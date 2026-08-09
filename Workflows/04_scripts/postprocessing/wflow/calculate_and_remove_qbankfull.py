#%%
# In this script, we calculate the 2-year return period of a 30yr wflow discharge simulation for different gauges, which are coupled later to SFINCS
# This bankfull discharge estimate is removed from the TC event discharge as an approximation for streamflow
# This script is based on: https://scaling-robot-wgkjqqr.pages.github.io/notebooks/Fit_univariate.html

# First, load the packages
import os
from datetime import datetime as datetime
from os.path import join
from hydromt_wflow import WflowSbmModel
import matplotlib.pyplot as plt
import pandas as pd
# NOTE: pyextremes is imported lazily further down, only on the branch that actually fits the
# extreme-value distribution. Importing it at module level made this script fail even when
# use_bankfull_corr is False, where no fitting happens at all.

#%%
# Set up wflow run variables (v1: hydromt.log.setuplog removed; no logger needed)
if "snakemake" in locals():
    wflow_root_30yr   = snakemake.params.wflow_root_forcing_30yr
    wflow_root_event  = snakemake.params.wflow_root_forcing
    data_cats         = snakemake.params.data_cat
    use_bankfull_corr = snakemake.params.use_bankfull_corr
    landuse_30yr      = snakemake.params.landuse_30yr
    precip_forcing    = snakemake.wildcards.precip_forcing
else:
    region            = "sofala"
    TC_name           = "Idai"
    precip_forcing    = "era5_hourly_zarr"
    CF_rain           = 0
    CF_rain_txt       = "0"
    CF_landuse        = "vito"
    landuse_30yr      = "vito"
    use_bankfull_corr = True
    wflow_root_30yr   = f"p:/11210471-001-compass/03_Runs/{region}/{TC_name}"
    wflow_root_event  = f"p:/11210471-001-compass/03_Runs/{region}/{TC_name}/wflow_{CF_landuse}/event_precip_{precip_forcing}_CF{CF_rain_txt}"
    curdir            = '../../../'
    data_cats         = [
        join(curdir, "03_data_catalogs", "datacatalog_general.yml"),
        join(curdir, "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling.yml"),
        join(curdir, "03_data_catalogs", "datacatalog_SFINCS_obspoints.yml"),
        join(curdir, "03_data_catalogs", "datacatalog_CF_forcing.yml")
        ]

#%%
# The 30-yr run is decoupled: it lives under its own land-use directory and is not
# necessarily the same land use as the event run.
wflow_path_30yr = join(wflow_root_30yr, f"wflow_{landuse_30yr}", f"event_precip_{precip_forcing}_CF0_30yr")
dis_out = os.path.join(wflow_root_event, "events", "run_default", "wflow_dis_no_bankfull.csv")

if not use_bankfull_corr:
    # Still emit the file the workflow expects, so downstream rules have a stable input.
    print("Not using bankfull correction, writing empty placeholder discharge file...")
    os.makedirs(os.path.dirname(dis_out), exist_ok=True)
    pd.DataFrame(list()).to_csv(dis_out)
    raise SystemExit(0)

# check whether the bankfull calculations have already been done
wflow_bankfull = f"{wflow_path_30yr}/warmup/qbankfull_wflow_gauges.csv"

if not os.path.exists(wflow_bankfull):
    from pyextremes import EVA   # only needed when the distribution is actually fitted

    # Read ('r') the Wflow 30yr warm-up results
    mod = WflowSbmModel(
        root=join(wflow_path_30yr, "warmup"),
        data_libs=data_cats,
        mode="r",
    )
    mod.read()

    # Read in the wflow discharge (v1: results['netcdf'] -> output_scalar component)
    df = mod.output_scalar.data['Q'].to_pandas()

    # Now we calculate the bankfull discharge, based on a 2 yr return period by using block maxima and fitting the distribution using Akaike Information Criterion (AIC)
    qbankfull = []

    # Loop over the different wflow gauges (output points)
    for gauge in df.columns:
        data = df[gauge]
        # Initialize a model for block maxima where we will store the results
        model_bm = EVA(data=data)

        # Sampling Annual Maxima, therefore using a block size of 365D
        peaks = model_bm.get_extremes(
            method="BM",
            extremes_type="high",
            block_size="365.2425D",
            errors="raise",
        )

        # Plot the selected extremes
        model_bm.plot_extremes()

        # We fit the extreme value models
        model_bm.fit_model(model="Emcee")

        # Estimate of return periods
        rp = [2] # 2-year return period for bankfull discharge
        summary = model_bm.get_summary(return_period=rp, alpha=0.95)
        print(summary)

        # Remove index and add gauge number
        summary = summary.reset_index(drop=True)
        summary["gauge"] = gauge

        # Add summary information to one df
        qbankfull.append(summary)

        # Plotting the fitting annual maxima
        fig, ax = model_bm.plot_return_values(alpha=0.95)
        ax.set_title(f"Gauge: {gauge}")
        ax.set_ylabel("Discharge [m³/s]")

    # Combine all into one DataFrame
    qbankfull_df = pd.concat(qbankfull, ignore_index=True)

    # Save to CSV
    qbankfull_df.to_csv(wflow_bankfull, index=False)

else:
    qbankfull_df = pd.read_csv(wflow_bankfull)


# %% ---------------------------------------------------------------
# Check removing bankfull discharge from factual event simulations
# ------------------------------------------------------------------
# Read ('r') the Wflow 30yr warm-up results
mod_F = WflowSbmModel(
    root=join(wflow_root_event, "events"),
    data_libs=data_cats,
    mode="r"
)
mod_F.read()
# %%
# Check results (v1: results['netcdf'] -> output_scalar component)
df_F = mod_F.output_scalar.data['Q'].to_pandas()
df_F

# %%
# Remove the qbankfull from all discharge values and set to zero if discharge is below 0
qbankfull_df = qbankfull_df.set_index('gauge')
qbankfull_df.index = qbankfull_df.index.astype(str)

# Ensure timeseries_df columns are strings for matching
df_F.columns = df_F.columns.astype(str)

df_F_no_bankfull = df_F.copy(deep=True)
for gauge in df_F_no_bankfull.columns:
    if gauge in qbankfull_df.index:
        qbankfull_gauge = qbankfull_df.loc[gauge, "return value"]
        df_F_no_bankfull[gauge] = df_F_no_bankfull[gauge] - qbankfull_gauge
        df_F_no_bankfull[gauge] = df_F_no_bankfull[gauge].clip(lower=0)     # ensures all values below 0 are set to 0
        df_F_no_bankfull.to_csv(dis_out, index=True)

# %% BANKFULL FIGURES
# Off by default: the plotting stack does not work in the compass-wflow pixi env, and this
# script runs as a snakemake rule where the figures are not needed.
make_bankfull_figures = False  # Set to True to produce the bankfull calculation figures
if make_bankfull_figures:
    # We select the first discharge location and have a look at the data
    data_F = df_F["1"]
    plt.figure()
    ax = data_F.plot()
    plt.ylabel("Discharge (m³/s)")

    # %% Plot the masked discharge compared to the full discharge
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot both time series on same axis
    df_F['1'].plot(ax=ax, label='Original', color='blue')
    df_F_no_bankfull['1'].plot(ax=ax, label='Masked', color='orange')

    # Add horizontal bankfull line
    ax.axhline(qbankfull_df.loc['1', 'return value'], color='red', linestyle=':', linewidth=2, label='Bankfull Q')

    ax.set_ylabel("Discharge (m³/s)")
    ax.set_title("Discharge with Bankfull Threshold for Gauge 1")
    ax.legend()

    plt.tight_layout()
    plt.show()

# %%

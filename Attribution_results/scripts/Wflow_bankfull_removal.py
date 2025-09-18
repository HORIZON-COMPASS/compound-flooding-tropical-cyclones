#%%
# In this script, we calculate the 2-year return period of a 30yr wflow discharge simulation for different gauges, which are coupled later to SFINCS
# This bankfull discharge estimate is removed from the TC event discharge as an approximation for streamflow
# This script is based on: https://scaling-robot-wgkjqqr.pages.github.io/notebooks/Fit_univariate.html

# First, load the packages using the pixi environement compass-wflow
import os
from datetime import datetime as datetime
from os.path import join
from hydromt.log import setuplog
from hydromt_wflow import WflowModel
from pyextremes import EVA
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt

#%%
# Set up wflow run variables
logger = setuplog("update", "./hydromt.log", log_level=10)

wflow_root_30yr  = f"../data/wflow/event_precip_era5_hourly_zarr_CF0_30yr"
wflow_root_event = f"../data/wflow/event_precip_era5_hourly_zarr_CF0"
maindir           = '../../'
data_cats        = [
    join(maindir, "Workflows", "03_data_catalogs", "datacatalog_general.yml"), 
    join(maindir, "Workflows", "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling.yml"), 
    join(maindir, "Workflows", "03_data_catalogs", "datacatalog_SFINCS_obspoints.yml"),
    join(maindir, "Workflows", "03_data_catalogs", "datacatalog_CF_forcing.yml")
    ]

#%%
# check whether the bankfull calculations have already been done
wflow_bankfull = f"{wflow_root_30yr}/warmup/qbankfull_wflow_gauges.csv"

if not os.path.exists(wflow_bankfull):
   # Read ('r') the Wflow 30yr warm-up results 
    mod = WflowModel(
        root=join(wflow_root_30yr, "warmup"),
        data_libs=data_cats,
        mode="r",
        logger=logger,
    )
    mod.read()

    # Read in the wflow discharge 
    df = mod.results['netcdf']['Q'].to_pandas()

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

    # Save to CSV - input Table S1
    qbankfull_df.to_csv(wflow_bankfull, index=False)

else:
    qbankfull_df = pd.read_csv(wflow_bankfull)


# %% ---------------------------------------------------------------
# Check removing bankfull discharge from factual event simulations
# ------------------------------------------------------------------
# Read ('r') the Wflow 30yr warm-up results 
mod_F = WflowModel(
    root=join(wflow_root_event, "events"),
    data_libs=data_cats,
    mode="r"
)
mod_F.read()
# %%
# Check results
df_F = mod_F.results['netcdf']['Q'].to_pandas()
df_F

# %%
# We select the first discharge location
data_F = df_F["1"]
# And have a look at the data
plt.figure()
ax = data_F.plot()
plt.ylabel("Discharge (m³/s)")

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


# %%
# Read the model
mod = WflowModel(
    root=join(wflow_root_30yr, "warmup"),
    data_libs=data_cats,
    mode="r",
    logger=logger,
)
mod.read()

# Read in the wflow discharge 
df = mod.results['netcdf']['Q'].to_pandas()

# Select first two gauges
gauges = df[['1', '2']]
rivers = ["Buzi", "Pungwe"]

# Create subplots for two gauges
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

qbankfull = []

for i, gauge in enumerate(gauges):
    data = df[gauge]
    model_bm = EVA(data=data)

    peaks = model_bm.get_extremes(
        method="BM",
        extremes_type="high",
        block_size="365.2425D",
        errors="raise",
    )

    model_bm.fit_model(model="Emcee")

    rp = [2]
    summary = model_bm.get_summary(return_period=rp, alpha=0.95)
    summary = summary.reset_index(drop=True)
    summary["gauge"] = gauge
    qbankfull.append(summary)

    # Plot on respective subplot
    model_bm.plot_return_values(alpha=0.95, ax=axes[i])
    # axes[i].plot(ax_sub.lines[0].get_xdata(), ax_sub.lines[0].get_ydata())
    axes[i].set_title(f"{rivers[i]} River (Gauge {i+1})")
    axes[i].set_ylabel("Discharge [m³/s]")
    axes[i].set_xlim(1,30)
    axes[i].text(0.02, 1.06, f"({chr(97+i)})", transform=axes[i].transAxes,
                 fontsize=12, fontweight='bold', va='top')
    
fig.savefig(f"../figures/fS2.png", dpi=300, bbox_inches='tight')
fig.savefig(f"../figures/fS2.pdf", dpi=300, bbox_inches='tight')

plt.tight_layout()
plt.show()

# %%
# Plot the masked discharge compared to the full discharge
fig, ax = plt.subplots(figsize=(12, 6))

# Plot both time series
df_F['1'].plot(ax=ax, label='Raw wflow output', color='steelblue', linewidth=1.8)
df_F_no_bankfull['1'].plot(ax=ax, label='Effective discharge', color='darkorange', linewidth=1.8)

# Add horizontal bankfull line
ax.axhline(
    qbankfull_df.loc['1', 'return value'],
    color='crimson', linestyle='--', linewidth=2,
    label='Estimated bankfull Q'
)

# Labels and title
ax.set_xlabel("Time", fontsize=13)       # X-axis label
ax.set_ylabel("Discharge [m³/s]", fontsize=13)
ax.set_title("Buzi River (Gauge 1): Discharge with bankfull threshold removed", fontsize=14)

# Tick label sizes
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)

# Grid for clarity
ax.grid(True, linestyle="--", alpha=0.6)

# Legend outside plot
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0, fontsize=11)

# Save
fig.savefig("../figures/fS3.png", dpi=300, bbox_inches='tight')
fig.savefig("../figures/fS3.pdf", dpi=300, bbox_inches='tight')

# %%

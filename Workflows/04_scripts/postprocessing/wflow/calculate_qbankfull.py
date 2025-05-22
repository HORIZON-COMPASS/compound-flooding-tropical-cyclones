#%%
# In this script, we calculate the 2-year return period of a 30yr wflow discharge simulation for different gauges, which are coupled later to SFINCS
# This bankfull discharge estimate is removed from the TC event discharge as an approximation for streamflow
# This script is based on: https://scaling-robot-wgkjqqr.pages.github.io/notebooks/Fit_univariate.html

# First, load the packages
from datetime import datetime as datetime
from os.path import join
from hydromt.log import setuplog
from hydromt_wflow import WflowModel
from pyextremes import EVA
import matplotlib.pyplot as plt
import pandas as pd

# Set up wflow run variables
logger = setuplog("update", "./hydromt.log", log_level=10)
wflow_root = f"p:/11210471-001-compass/03_Runs/sofala/Idai/wflow/event_precip_era5_hourly_zarr_CF0_30yr"
curdir              = '../../../'
data_cats           = [
        join(curdir, "03_data_catalogs", "datacatalog_general.yml"), 
        join(curdir, "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling.yml"), 
        join(curdir, "03_data_catalogs", "datacatalog_SFINCS_obspoints.yml"),
        join(curdir, "03_data_catalogs", "datacatalog_CF_forcing.yml")
        ]

# Read ('r') the Wflow 30yr warm-up results 
mod = WflowModel(
    root=join(wflow_root, "warmup"),
    data_libs=data_cats,
    mode="r",
    logger=logger,
)
mod.read()
#%%
# Read in the wflow discharge 
df = mod.results['netcdf']['Q'].to_pandas()
df

# %%
# We select the first discharge location
data = df["1"]
# And have a look at the data
plt.figure()
ax = data.plot()
plt.ylabel("Discharge (m³/s)")

# %%
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
    
#%%
# Combine all into one DataFrame
qbankfull_df = pd.concat(qbankfull, ignore_index=True)
qbankfull_df
#%%
# Save to CSV
file_path = f"{wflow_root}/warmup"
qbankfull_df.to_csv(os.path.join(file_path,"qbankfull_wflow_gauges.csv"), index=False)



# %% ---------------------------------------------------------------
# Check removing bankfull discharge from factual event simulations
# ------------------------------------------------------------------
wflow_event_root = f"p:/11210471-001-compass/03_Runs/sofala/Idai/wflow/event_precip_era5_hourly_zarr_CF0"

# Read ('r') the Wflow 30yr warm-up results 
mod_F = WflowModel(
    root=join(wflow_event_root, "events"),
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
# Mask values below 2-year return value to only keep flood flows 
qbankfull_df = qbankfull_df.set_index('gauge')
qbankfull_df.index = qbankfull_df.index.astype(str)

# Ensure timeseries_df columns are strings for matching
df_F.columns = df_F.columns.astype(str)

masked_df = df_F.copy(deep=True)
for gauge in masked_df.columns:
    if gauge in qbankfull_df.index:
        threshold = qbankfull_df.loc[gauge, "return value"]
        masked_df[gauge] = masked_df[gauge].where(masked_df[gauge] > threshold)

# %%
# Plot the masked discharge compared to the full discharge
fig, ax = plt.subplots(figsize=(12, 6))

# Plot both time series on same axis
df_F['1'].plot(ax=ax, label='Original', color='blue')
masked_df['1'].plot(ax=ax, label='Masked', color='orange')

# Add horizontal bankfull line
ax.axhline(qbankfull_df.loc['1', 'return value'], color='red', linestyle=':', linewidth=2, label='Bankfull Q')

ax.set_ylabel("Discharge (m³/s)")
ax.set_title("Discharge with Bankfull Threshold for Gauge 1")
ax.legend()

plt.tight_layout()
plt.show()

# %%

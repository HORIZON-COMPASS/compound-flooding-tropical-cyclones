# in this test script, the factuals and counterfactual runs are compared
#%% Import packages
import os
from pathlib import Path
from os.path import join
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from hydromt_sfincs import SfincsModel, utils
import contextily as ctx  # For adding basemap tiles
import seaborn as sns
import yaml
import itertools
import cartopy.crs as ccrs
import math

# %% Load the SFINCS model runs
config_general  = '../../Workflows/01_config_snakemake/config_general_MZB.yml'
# Load the YAML file
with open(config_general, 'r') as file:
    config = yaml.safe_load(file)

#%%
# Access values
runname_id          = config['runname_ids']['Idai']
tc_name             = runname_id['tc_name']
sid                 = runname_id['sid']
region              = runname_id['region']
bbox_sfincs         = runname_id['bbox_sfincs']
bbox_dfm            = runname_id['bbox_dfm']
start_time          = runname_id['start_time']
end_time            = runname_id['end_time']
ext_days            = runname_id['ext_days']
precip_forcing      = runname_id['precip_forcing']
dfm_res             = runname_id['dfm_res']
bathy               = runname_id['bathy']
wind_forcing        = runname_id['wind_forcing']
tidemodel           = runname_id['tidemodel']
dfm_obs_file        = runname_id['dfm_obs_file']
sfincs_obs_file     = runname_id['sfincs_obs_file']
dfm_verification_points = runname_id['dfm_verification_points']
config_sfincs_base  = runname_id['config_sfincs_base']
utmzone             = runname_id['utmzone']

# Access CF values
CF_value_rain = runname_id['CF_value_rain']
CF_value_wind = runname_id['CF_value_wind']
CF_value_SLR  = runname_id['CF_value_SLR']
#%%
# Define model path
# model_name     = f'event_tp_{precip_forcing}_CF{CF_value_rain}_{tidemodel}_CF{CF_value_SLR}_{wind_forcing}_CF{CF_value_wind}'
base_path     = f'p:/11210471-001-compass/03_Runs/{region}/{tc_name}/sfincs/'

#%%
# Load model output
# List to store all models
models = []
# Generate all combinations of CF values
for rain, wind, slr in itertools.product(CF_value_rain, CF_value_wind, CF_value_SLR):
    model_name = f'event_tp_{precip_forcing}_CF{rain}_{tidemodel}_CF{slr}_{wind_forcing}_CF{wind}'
    model_path = f'{base_path}/{model_name}'
    
    print(f'Model Name: {model_name}')
    print(f'Model Path: {model_path}')

    # Categorization logic
    CF_values = [rain, wind, slr]
    num_CF_diff = sum(1 for v in CF_values if v != 0)

    if num_CF_diff == 0:
        category = "Factual"
        cat_short = "F"
    elif num_CF_diff == 1:
        category = "Single Driver Counterfactual"
        cat_short = "CF_DR_single"
    elif num_CF_diff == 2:
        category = "Counterfactual Driver Pair"
        cat_short = "CF_DR_pair"
    else:
        category = "Counterfactual Compound Driver"
        cat_short = "CF_DR_compound"

    CF_info = { "rain": rain, "wind": wind, "SLR": slr }
    CF_info = { param: value for param, value in CF_info.items() if value != '0' } 

    # Read Sfincs model and store in list
    model_obj = SfincsModel(model_path, mode="r")
    models.append({
        "model_name": model_name,
        "model_path": model_path,
        "sfincs_model": model_obj,
        "category": category,
        "CF_info": CF_info
    })

# Print results for verification
for model in models:
    print(f"{model['CF_info']}")
    print(f"{model['category']} -> {model['model_name']} -> {model['model_path']} -> CF_info: {model['CF_info']}")

# %%
# read global surface water occurance (GSWO) data to mask permanent water
models[0]["sfincs_model"].data_catalog.from_yml(os.path.join('../../Workflows/03_data_catalogs/datacatalog_general.yml'))
gswo = models[0]["sfincs_model"].data_catalog.get_rasterdataset("gswo", geom=models[0]["sfincs_model"].region, buffer=1000)

# %%
### loop over the different sfincs_model to compute hmax ###
# we set a threshold to mask minimum flood depth
hmin = 0.05

# Can be turned into a function!
for model in models:
    # if "hmax_masked" in models["sfincs_model"].results:
    #     print("hmax_masked exists!")
    # else:
    data        = model["sfincs_model"]  # Load the SFINCS model object
    sfincs_root = model["model_path"]
    title       = model["model_name"]
    print(f"Processing model: {title}")
    # first we are going to select our highest-resolution elevation dataset
    depfile = join(sfincs_root, "subgrid", "dep_subgrid.tif")
    da_dep = data.data_catalog.get_rasterdataset(depfile)

    # as we have a subgrid model, we don't have hmax available, so we are using zsmax (maximum water levels)
    # compute the maximum over all time steps
    da_zsmax = data.results["zsmax"].max(dim="timemax")
    # Fourthly, we downscale the floodmap
    da_hmax = utils.downscale_floodmap(
        zsmax=da_zsmax,
        dep=da_dep,
        hmin=hmin,
        # floodmap_fn=join(sfincs_root, "gis/floodmap.tif") # uncomment to save floodmap to <mod.root>/floodmap.tif
        )

    # we use the GSWO dataset to mask permanent water by first geprojecting it to the subgrid of hmax
    gswo_mask = gswo.raster.reproject_like(da_hmax, method="max")
    # permanent water where water occurence > 5%
    da_hmax_masked = da_hmax.where(gswo_mask <= 5)
    # Add the name attribute for identification
    data.results['hmax'] = da_hmax
    data.results['hmax_masked'] = da_hmax_masked

        # Open the existing NetCDF dataset in append mode
        # nc_file = join(sfincs_root, "sfincs_his.nc")
        # nc_file_out = join(sfincs_root, "sfincs_his_test.nc")
        
        # try:
        #     with xr.open_dataset(nc_file, mode="a") as ds:  
        #         # Convert DataArrays to Dataset with proper variable names
        #         ds_new = xr.Dataset({
        #             "hmax": da_hmax,
        #             "hmax_masked": da_hmax_masked
        #         })
                
        #         # Merge with existing dataset and save
        #         ds_updated = xr.merge([ds, ds_new])
        #         ds_updated.to_netcdf(nc_file, mode="w")  # Overwrite with new data

        #     print(f"Saved hmax and hmax_masked to {nc_file}")

        # except Exception as e:
        #     print(f"Error saving to {nc_file}: {e}")
        
    del da_dep, da_zsmax, da_hmax, da_hmax_masked

# %% Create subplots and set up figure title
projection = models[0]['sfincs_model'].crs.to_epsg()
# Determine subplot grid size
num_models = len(models)
num_cols = 2  # Adjust as needed
num_rows = math.ceil(num_models / num_cols)

fig, axes = plt.subplots(nrows=num_rows, 
                         ncols=num_cols, 
                         figsize=(8, 10), 
                         constrained_layout=True, 
                         subplot_kw={"projection": ccrs.epsg(projection)})
fig.suptitle("Masked hmax", fontsize=12, y=1.02)

# Loop through datasets and plot
for model, ax in zip(models, axes.flatten()):
    # Plot the masked water depth and add basemap
    im = model['sfincs_model'].results['hmax_masked'].plot.pcolormesh(
        ax=ax, cmap="Blues", vmin=0, vmax=3.0, add_colorbar=False
    )
    ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=12, crs=model['sfincs_model'].results['hmax_masked'].rio.crs, attribution=False)
    
    # Set axis labels and title
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(model['category'], fontsize=11)

# Add colorbar and label
fig.colorbar(im, ax=axes, orientation="vertical", fraction=0.02, pad=0.04).set_label('Masked hmax (m)', rotation=270, labelpad=15)

# Show the plot
plt.show()

#%% Create plots of the difference between F and CF scenarios
# Calculate the differences
factual_hmax = None 

for model in models:
        # Store factual hmax for comparison
    if model["category"] == "Factual":
        factual_hmax = model['sfincs_model'].results["hmax_masked"]

    # Compute difference for counterfactual models
    hmax_diff = None
    if factual_hmax is not None and model["category"] != "Factual":
        hmax_diff = factual_hmax - model['sfincs_model'].results["hmax_masked"]
        model['sfincs_model'].results["hmax_diff"] = hmax_diff
        print(f"hmax_diff calculated for {model['model_name']}")

#%%
# Create subplots
# Filter out only counterfactual models with hmax_diff
counterfactual_models = [m for m in models if "hmax_diff" in m["sfincs_model"].results]

# Determine subplot grid size
num_models = len(counterfactual_models)
num_cols = 2  # Adjust as needed
num_rows = math.ceil(num_models / num_cols)

projection = models[0]['sfincs_model'].crs.to_epsg()
# Create figure and subplots
fig, axes = plt.subplots(nrows=num_rows, 
                         ncols=num_cols, 
                         figsize=(14, 8), 
                         constrained_layout=True, 
                         subplot_kw={"projection": ccrs.epsg(projection)}
                         )
fig.suptitle("Difference Plots: Factual - Counterfactual Water Level", fontsize=16, y=1.06)

# Flatten axes array for easy indexing (works for any grid size)
axes = axes.flatten()

# Define colormap settings
cmap = "RdBu"
vmin, vmax = -0.3, 0.3

# Loop through counterfactual models and plot
for idx, model in enumerate(counterfactual_models):
    ax = axes[idx]
    hmax_diff = model['sfincs_model'].results["hmax_diff"]

    # Plot difference
    im = hmax_diff.plot.pcolormesh(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False)

    # Add basemap
    ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=12, crs=model['sfincs_model'].crs, attribution=False)

    # Set title and labels
    ax.set_title(f"{model['model_name']}", fontsize=12)
    ax.set_xlabel("x coordinate UTM zone 36S [m]", fontsize=10)
    ax.set_ylabel("y coordinate UTM zone 36S [m]", fontsize=10)
    ax.tick_params(axis="both", labelsize=9)

# Hide any unused subplots (if the grid is larger than the number of models)
for idx in range(len(counterfactual_models), len(axes)):
    fig.delaxes(axes[idx])

# Add shared colorbar
cbar = fig.colorbar(im, ax=axes, orientation="vertical", fraction=0.02, pad=0.04)
cbar.set_label('Difference in Maximum Water Level (m)', rotation=270, labelpad=20, fontsize=12)
cbar.ax.tick_params(labelsize=10)

# Show the plot
plt.show()


# %% Calculate the surface area of one grid cell
dx = abs(models[0]['sfincs_model'].results['hmax_masked'].x[1] - models[0]['sfincs_model'].results['hmax_masked'].x[0])  # Grid resolution in x-direction (meters)
dy = abs(models[0]['sfincs_model'].results['hmax_masked'].y[1] - models[0]['sfincs_model'].results['hmax_masked'].y[0])  # Grid resolution in y-direction (meters)
cell_area = dx * dy 

# %% Calculate the area of the flooded cells as flood extent with hmax_masked as variable
flood_characs = []

for model in models:
    flooded_cells = model['sfincs_model'].results['hmax_masked'] > 0  # Create a boolean mask
    flood_extent = (flooded_cells * cell_area).sum().compute()  # Compute the total flooded area

    # Convert the result to square kilometers
    flood_extent_km2 = flood_extent / 1e6  # Convert square meters to square kilometers
    
    # Store the result in the dictionary
    model['flood_extent_km2'] = flood_extent_km2
    print(f"for model {model['model_name']}, the flooded area: {flood_extent_km2}")


# %% Calculate the volume for each flooded cell (depth * area)
for model in models:
    flooded_cells = model['sfincs_model'].results['hmax_masked'] > 0  # Create a boolean mask
    flood_volume = (model['sfincs_model'].results['hmax_masked'] * flooded_cells * cell_area).sum().compute()  # Sum the volume of all flooded cells

    # Convert flood volume to cubic kilometers (optional, for readability)
    flood_volume_km3 = flood_volume / 1e9  # Convert cubic meters to cubic kilometers

    model['flood_volume_km3'] = flood_volume_km3
    print(f"for model {model['model_name']}, the flood volume is {flood_volume_km3}")

#%% Calculate the flood volume difference between the factual and counterfactual datasets:

factual_flood_volume = None 
factual_flood_extent = None 

for model in models:
        # Store factual hmax for comparison
    if model["category"] == "Factual":
        factual_flood_volume = model['flood_volume_km3']
        factual_flood_extent = model['flood_extent_km2']

    # Compute difference for counterfactual models
    flood_volume_diff = None

    if factual_flood_volume is not None and model["category"] != "Factual":
        flood_volume_diff = (factual_flood_volume - model['flood_volume_km3']) / factual_flood_volume * 100
        model['Volume_diff_from_F(%)'] = flood_volume_diff
        print(f"flood_volume_diff calculated for {model['model_name']}")

    flood_extent_diff = None
    if factual_flood_extent is not None and model["category"] != "Factual":
        flood_extent_diff = (factual_flood_extent - model['flood_extent_km2']) / factual_flood_extent * 100
        model['Extent_diff_from_F(%)'] = flood_extent_diff
        print(f"flood_extent_diff calculated for {model['model_name']}")

#%%
# Now create a DataFrame for plotting, excluding models with category "Factual"
model_names = []
volume_diffs = []
drivers = []

for model in models:
    if model['category'] != "Factual":
        model_names.append(model["model_name"])
        volume_diff_value = model["Volume_diff_from_F(%)"].values.flatten()[0]
        volume_diffs.append(volume_diff_value)

        CF_info = model['CF_info']
        # Check if each driver is non-zero and construct the description
        if CF_info.get("rain") != 0 and CF_info.get("wind") != 0 and CF_info.get("SLR") != 0:
            drivers.append("Compound")
        elif CF_info.get("rain") != 0 and CF_info.get("wind") != 0:
            drivers.append("Rain & wind")
        elif CF_info.get("rain") != 0 and CF_info.get("SLR") != 0:
            drivers.append("Rain & SLR")
        elif CF_info.get("wind") != 0 and CF_info.get("SLR") != 0:
            drivers.append("Wind & SLR")
        elif CF_info.get("rain") != 0:
            drivers.append("Rain")
        elif CF_info.get("wind") != 0:
            drivers.append("Wind")
        elif CF_info.get("SLR") != 0:
            drivers.append("SLR")

# Create a DataFrame for seaborn plotting
df = pd.DataFrame({
    "Model name": model_names,
    "Flood volume difference (F-CF in %)": volume_diffs,
    "CF driver": drivers
})

sns.catplot(data=df, 
            x="CF driver", 
            y="Flood volume difference (F-CF in %)", 
            kind="swarm", 
            height=6, 
            aspect=2)


# %%

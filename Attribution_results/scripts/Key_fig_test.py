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
        "cat_short": cat_short,
        "CF_info": CF_info
    })

extra_model_obj = SfincsModel("p:/11210471-001-compass/03_Runs/sofala/Idai/sfincs/event_tp_era5_hourly_CF0_GTSMv41opendap_CF-0.14_toSFINCSwaterlevel_spw_IBTrACS_CF0", mode="r")
extra_CF_info = {"rain": 0, "wind": 0, "SLR": -0.14 }
models.append({
    "model_name": "event_tp_era5_hourly_CF0_GTSMv41opendap_CF-0.14_toSFINCSwaterlevel_spw_IBTrACS_CF0",
    "model_path": "p:/11210471-001-compass/03_Runs/sofala/Idai/sfincs/event_tp_era5_hourly_CF0_GTSMv41opendap_CF-0.14_toSFINCSwaterlevel_spw_IBTrACS_CF0",
    "sfincs_model": extra_model_obj,
    "category": "Single Driver Counterfactual",
    "cat_short": "CF_DR_single",
    "CF_info": { param: value for param, value in extra_CF_info.items() if value != '0' } 
})

# Print results for verification
for model in models:
    print(f"{model['CF_info']}")
    print(f"{model['category']} -> {model['model_name']} -> {model['model_path']} -> CF_info: {model['CF_info']}")

# %%
# read global surface water occurance (GSWO) data to mask permanent water
models[0]["sfincs_model"].data_catalog.from_yml(os.path.join('../../Workflows/03_data_catalogs/datacatalog_general.yml'))
gswo = models[0]["sfincs_model"].data_catalog.get_rasterdataset("gswo", geom=models[0]["sfincs_model"].region, buffer=1000)

# %% Load the SFINCS models
root         = '../../sfincs_models'
TC           = 'Idai'
scenario_F   = 'factuals'
scenario_CF  = 'counterfactuals'
comp_drivers = 'compound_drivers'
sing_drivers = 'single_drivers'


### Factuals ###
# compound drivers
sfincs_root_F_PluvCoast = os.path.join(root,TC,scenario_F,comp_drivers,'sfincs_MZ_ERA5Land_nodischarge')
mod_F_PluvCoast = SfincsModel(sfincs_root_F_PluvCoast, mode="r")

# single drivers 
sfincs_root_F_Coast = os.path.join(root,TC,scenario_F,sing_drivers,'sfincs_MZ_coastal')
mod_F_Coast = SfincsModel(sfincs_root_F_Coast, mode="r")

sfincs_root_F_Coast_MDT = os.path.join(root,TC,scenario_F,sing_drivers,'sfincs_MZ_coastal_MDT')
mod_F_Coast_MDT = SfincsModel(sfincs_root_F_Coast_MDT, mode="r")


### Counterfactuals ###
# compound drivers
sfincs_root_CF_7precip_PluvCoast = os.path.join(root,TC,scenario_CF,comp_drivers,'sfincs_MZ_ERA5Land_CF7%_compd')
mod_CF_7precip_PluvCoast = SfincsModel(sfincs_root_CF_7precip_PluvCoast, mode="r")

# single drivers
# sfincs_root_CF_noSLR_Coast = os.path.join(root,TC,scenario_CF,sing_drivers,'sfincs_MZ_SLR_isimip2015')
# mod_CF_noSLR_Coast = SfincsModel(sfincs_root_CF_noSLR_Coast, mode="r")

sfincs_root_CF_noSLR_Coast_MDT = os.path.join(root,TC,scenario_CF,sing_drivers,'sfincs_MZ_coastal_MDT_noSLR')
mod_CF_noSLR_Coast_MDT = SfincsModel(sfincs_root_CF_noSLR_Coast_MDT, mode="r")

# %%
datasets = [
    (mod_F_PluvCoast, sfincs_root_F_PluvCoast, "F_PluvCoast", "Factual", "Pluvial", "precipitation"),
    (mod_CF_7precip_PluvCoast, sfincs_root_CF_7precip_PluvCoast, "CF_7precip_PluvCoast", "Counterfactual", "Pluvial", "precipitation"),
    (mod_F_Coast, sfincs_root_F_Coast, "F_SLR_Coast", "Factual", "Coastal", "SLR-noMDT"),
    (mod_F_Coast_MDT,sfincs_root_F_Coast_MDT, "F_SLR_Coast_MDT", "Factual", "Coastal", "SLR"),
    (mod_CF_noSLR_Coast_MDT, sfincs_root_CF_noSLR_Coast_MDT, "CF_noSLR_Coast_MDT", "Counterfactual", "Coastal", "noSLR")
]

# %%
# read global surface water occurance (GSWO) data to mask permanent water
mod_F_PluvCoast.data_catalog.from_yml(os.path.join('../../COMPASS/Workflows/03_data_catalogs/datacatalog_general.yml'))
gswo = mod_F_PluvCoast.data_catalog.get_rasterdataset("gswo", geom=mod_F_PluvCoast.region, buffer=1000)

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

for (data, sfincs_root, title, scenario, flood_type, driver) in (datasets):
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
                #   ds.close()

        #     ds_updated.to_netcdf(nc_file, mode="w")  # Overwrite with new data

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
    im = model['sfincs_model'].results()['hmax_masked'].plot.pcolormesh(
        ax=ax, cmap="Blues", vmin=0, vmax=3.0, add_colorbar=False
    )
    ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=12, crs=model['sfincs_model'].results['hmax_masked'].rio.crs, attribution=False)

    # Add the name attribute for identification
    data.results['hmax_masked'] = da_hmax_masked
    
del da_dep, da_zsmax, da_hmax, da_hmax_masked

# %% Create subplots and set up figure title
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(8, 15), constrained_layout=True)
fig.suptitle("Masked hmax: pluvial & coastal flooding", fontsize=16, y=1.02)

# Define high-resolution basemap source
basemap_source = ctx.providers.Esri.WorldImagery

# Loop through datasets and plot
for (data, sfincs_root, name, scenario, flood_type, driver), ax in zip(datasets, axes.flatten()):
    # Plot the masked water depth and add basemap
    im = data.results['hmax_masked'].plot.pcolormesh(
        ax=ax, cmap="Blues", vmin=0, vmax=3.0, add_colorbar=False
    )
    ctx.add_basemap(ax, source=basemap_source, zoom=12, crs=data.results['hmax_masked'].rio.crs, attribution=False)
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
    non_zero_CF_info = {key: value for key, value in model["CF_info"].items() if value != 0}
    cf_info_str = ", ".join(f"{key}: {value}" for key, value in non_zero_CF_info.items())

    ax.set_title(f"{model['cat_short']} ({cf_info_str})", fontsize=12)
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

# Extract individual datasets from the collection
mod_F_PluvCoast = datasets[0][0]
mod_CF_7precip_PluvCoast = datasets[1][0]
mod_F_SLR_Coast_MDT = datasets[3][0]
mod_CF_noSLR_Coast_MDT = datasets[4][0]

# Calculate differences
difference_pluv_coast = mod_F_PluvCoast.results['hmax_masked'] - mod_CF_7precip_PluvCoast.results['hmax_masked']
difference_slr_coast_mdt = mod_F_SLR_Coast_MDT.results['hmax_masked'] - mod_CF_noSLR_Coast_MDT.results['hmax_masked']

# Create subplots
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 8), constrained_layout=True)
fig.suptitle("Difference Plots: Factual - Counterfactual waterlevel", fontsize=16, y=1.06)

# Define high-resolution basemap source
basemap_source = ctx.providers.Esri.WorldImagery

# Plot Difference 1: Pluvial & Coastal Flooding
im1 = difference_pluv_coast.plot.pcolormesh(
    ax=axes[0], cmap="RdBu", vmin=-0.3, vmax=0.3, add_colorbar=False
)
ctx.add_basemap(
    axes[0],
    source=basemap_source,
    zoom=12,
    crs=difference_pluv_coast.rio.crs,
    attribution=False
)
axes[0].set_title(r"Pluvial & coastal flooding: CF -7% rainfall", fontsize=14)
axes[0].set_xlabel("x coordinate UTM zone 36S [m]", fontsize=12)
axes[0].set_ylabel("y coordinate UTM zone 36S [m]", fontsize=12)
axes[0].tick_params(axis="both", labelsize=11)

# Plot Difference 2: Coastal MDT
im2 = difference_slr_coast_mdt.plot.pcolormesh(
    ax=axes[1], cmap="RdBu", vmin=-0.3, vmax=0.3, add_colorbar=False
)
ctx.add_basemap(
    axes[1],
    source=basemap_source,
    zoom=12,
    crs=difference_slr_coast_mdt.rio.crs,
    attribution=False
)
axes[1].set_title("Coastal flooding: CF -0.14 m SLR", fontsize=14)
axes[1].set_xlabel("x coordinate UTM zone 36S [m]", fontsize=12)
axes[1].set_ylabel("y coordinate UTM zone 36S [m]", fontsize=12)
axes[1].tick_params(axis="both", labelsize=11)

# Add shared colorbar
cbar = fig.colorbar(im1, ax=axes, orientation="vertical", fraction=0.02, pad=0.04)
cbar.set_label('Difference in maximum waterlevel (m)', rotation=270, labelpad=20, fontsize=14)
cbar.ax.tick_params(labelsize=12)

# Show the plot
plt.show()


# %% Calculate the surface area of one grid cell
dx = abs(models[0]['sfincs_model'].results['hmax_masked'].x[1] - models[0]['sfincs_model'].results['hmax_masked'].x[0])  # Grid resolution in x-direction (meters)
dy = abs(models[0]['sfincs_model'].results['hmax_masked'].y[1] - models[0]['sfincs_model'].results['hmax_masked'].y[0])  # Grid resolution in y-direction (meters)

# %% Calculate the surface area of one grid cell
dx = abs(mod_F_PluvCoast.results['hmax_masked'].x[1] - mod_F_PluvCoast.results['hmax_masked'].x[0])  # Grid resolution in x-direction (meters)
dy = abs(mod_F_PluvCoast.results['hmax_masked'].y[1] - mod_F_PluvCoast.results['hmax_masked'].y[0])  # Grid resolution in y-direction (meters)
cell_area = dx * dy 

# %% Calculate the area of the flooded cells as flood extent with hmax_masked as variable
flood_characs = []

for model in models:
    flooded_cells = model['sfincs_model'].results['hmax_masked'] > 0  # Create a boolean mask

for (data, sfincs_root, name, scenario, flood_type, driver) in (datasets):
    flooded_cells = data.results['hmax_masked'] > 0  # Create a boolean mask
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
    flood_characs.append({
        'name': name,
        'scenario': scenario,
        'flood_type': flood_type,
        'driver': driver,
        'flood_extent_km2': flood_extent_km2
    })

#%% Print results
for i, data in enumerate(flood_characs):
    print(f"Dataset: {data['name']}, Flooded Area: {data['flood_extent_km2'].item()} km²")

# %% Calculate the volume for each flooded cell (depth * area)
for (data, sfincs_root, name, scenario, flood_type, driver), characs in zip(datasets, flood_characs):
    flooded_cells = data.results['hmax_masked'] > 0  # Create a boolean mask
    flood_volume = (data.results['hmax_masked'] * flooded_cells * cell_area).sum().compute()  # Sum the volume of all flooded cells

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


#%%
# List to store all data
data_list = []

for model in valid_models:
    model_path = model["model_path"]
    hisfile = os.path.join(model_path, "sfincs_his.nc")

    # Open dataset
    ds_his = xr.open_dataset(hisfile)
    ds_his["station_id"] = ds_his["station_id"].astype(int)

    # Convert to Pandas DataFrame
    df = ds_his[["point_zs"]].to_dataframe().reset_index()
    
    # Add model name for identification
    non_zero_CF_info = {key: value for key, value in model["CF_info"].items() if value != 0}
    cf_info_str = ", ".join(f"{key}: {value}" for key, value in non_zero_CF_info.items())
    
    df["model_name"]  = model["model_name"]
    df["category"]    = model["category"]
    df["cat_short"]   = model["cat_short"]
    df["CF_info"]     = model["CF_info"]
    df["CF_info_str"] = cf_info_str

    # Append to list
    data_list.append(df)

# Combine all models into one DataFrame
df_all = pd.concat(data_list, ignore_index=True)


#%%
# plot timeseries at output points 
# Get unique station IDs
station_ids = df_all["station_id"].unique()
nrows = len(station_ids)  # One row per station

# Create figure and axes with higher resolution and shared x-axis
fig, axes = plt.subplots(nrows=nrows, ncols=1, figsize=(14, 40), sharex=True, constrained_layout=True)

# Ensure axes is iterable even when there's only one subplot
if nrows == 1:
    axes = [axes]

# Plot data for each station
for ax, station_id in zip(axes, station_ids):
    df_station = df_all[df_all["station_id"] == station_id]
    
    for model_name in df_station["model_name"].unique():
        df_model = df_station[df_station["model_name"] == model_name]
        ax.plot(df_model["time"], df_model["point_zs"], label=model_name, linewidth=2)

    ax.set_title(f"Station {station_id}", fontsize=18)
    ax.set_ylabel("Water Level [m]", fontsize=14)
    ax.grid()
    
    # Increase y-axis tick size
    ax.tick_params(axis="y", labelsize=14)

# Shared x-axis label
axes[-1].set_xlabel("Time", fontsize=16)
axes[-1].tick_params(axis="x", rotation=60, labelsize=14)  # Rotate x-ticks for clarity

# Add a single legend at the top
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=16)

    if name == characs['name']:
        characs['flood_volume_km3'] = flood_volume_km3

#%% Print results
for i, data in enumerate(flood_characs):
    print(f"Dataset: {data['name']}, Flood Volume: {data['flood_volume_km3'].item()} km³")
#%% Calculate the flood volume difference between the factual and counterfactual datasets:
for data in flood_characs:
    if data['flood_type'] == 'Pluvial' and data['scenario'] == 'Counterfactual' and data['driver'] == 'precipitation':
        counterfactual_volume = data['flood_volume_km3']
        factual_data = next(d for d in flood_characs if d['flood_type'] == 'Pluvial' and d['scenario'] == 'Factual' and d['driver'] == 'precipitation')
        factual_volume = factual_data['flood_volume_km3']

        volume_dif_prc = (factual_volume - counterfactual_volume)/factual_volume * 100
        data['Volume_diff_from_F(%)'] = volume_dif_prc


    if data['flood_type'] == 'Coastal' and data['scenario'] == 'Counterfactual' and data['driver'] == 'noSLR':
        counterfactual_volume = data['flood_volume_km3']
        factual_data = next(d for d in flood_characs if d['flood_type'] == 'Coastal' and d['scenario'] == 'Factual' and d['driver'] == 'SLR')
        factual_volume = factual_data['flood_volume_km3']

        volume_dif_prc = (factual_volume - counterfactual_volume)/factual_volume * 100
        data['Volume_diff_from_F(%)'] = volume_dif_prc

#%%
# Sample Data
categories = ['Flood drivers']
data = {
    'Precipitation': [np.random.normal(loc=5, scale=2, size=100)],
    'Wind': [np.random.normal(loc=2, scale=1, size=100)],
    'SLR': [np.random.normal(loc=3, scale=1, size=100)],
    'Compound': [np.random.normal(loc=6, scale=1, size=100)],
}

# Define spacing parameters
n_sets = len(data)
n_categories = len(categories)
group_width = 0.8  # Reduce group width to tighten spacing
set_spacing = 0.2  # Less spacing between boxes
box_width = 0.12  # Keep boxes narrow
positions = []

# Compute positions for each box
for i, category in enumerate(categories):
    base = i * (group_width + 1.5)  # Adjust category spacing
    positions.extend([base + j * set_spacing for j in range(n_sets)])

# Flatten data for plotting
flattened_data = [item for sublist in zip(*data.values()) for item in sublist]

# Plot
fig, ax = plt.subplots(figsize=(8, 5))  # Smaller figure size
box = ax.boxplot(flattened_data, positions=positions, widths=box_width, patch_artist=True)

# Add colors for each set
colors = ['lightblue', 'lightgreen', 'coral', 'pink']
for patch, color in zip(box['boxes'], colors * n_categories):
    patch.set_facecolor(color)

# Custom x-ticks for categories
category_centers = [(i * (group_width + 1.5) + (group_width - set_spacing) / 2) for i in range(n_categories)]
ax.set_xticks(category_centers)
ax.set_xticklabels(categories, fontsize=14)

# Labels and grid
ax.set_ylabel('Diff in flood volume [%]', fontsize=14)
ax.set_title('Climate Attribution of Driver Contributions to Compound Flood from TC Idai', fontsize=16)
ax.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Zero line for positive/negative values
ax.grid(axis='y', linestyle='--', alpha=0.7)

# Ensure x-axis label is properly centered
ax.set_xlabel('')

# Legend
handles = [plt.Line2D([0], [0], color=color, lw=4, label=key) for key, color in zip(data.keys(), colors)]
ax.legend(handles, data.keys(), title='Driver', loc='upper left', fontsize=14, title_fontsize=14)

plt.tight_layout()
plt.show()


#%%
# Define flood types and their suffixes
flood_types = {
    'Pluvial': 'pluv',
    'Fluvial': 'fluv',
    'Coastal': 'coast',
    'Compound': 'cmpd'
}

# Define drivers
drivers = ['precipitation', 'temperature', 'wind', 'noSLR']

# Create dictionaries to store values for each flood type and driver
driver_contributions = {
    'precipitation': {},
    'temperature': {},
    'wind': {},
    'noSLR': {}
}

# Iterate through flood types and drivers to generate np.array datasets
for flood_type, suffix in flood_types.items():
    for driver in drivers:
        # Filter and create an array for the current flood type and driver
        driver_values = np.array([
            d['Volume_diff_from_F(%)'] for d in flood_characs
            if d['flood_type'] == flood_type and d['scenario'] == 'Counterfactual' and d['driver'] == driver
        ])
        
        # Ensure an empty array is created if no matching values are found
        if driver_values.size == 0:
            driver_values = np.array([])

        # Use the flood type suffix as the key
        dataset_name = f"{suffix}"
        driver_contributions[driver][dataset_name] = driver_values

#%%
### Make the plot ###
categories = ['Pluvial', 'Fluvial', 'Coastal', 'Compound']
data = {
    'Precipitation': [
        driver_contributions['precipitation']['pluv'],
        driver_contributions['precipitation']['fluv'],
        driver_contributions['precipitation']['coast'],
        driver_contributions['precipitation']['cmpd']
    ],
    'Temperature': [
        driver_contributions['temperature']['pluv'],
        driver_contributions['temperature']['fluv'],
        driver_contributions['temperature']['coast'],
        driver_contributions['temperature']['cmpd']
    ],
    'Wind': [
        driver_contributions['wind']['pluv'],
        driver_contributions['wind']['fluv'],
        driver_contributions['wind']['coast'],
        driver_contributions['wind']['cmpd']
    ],
    'SLR': [
        driver_contributions['noSLR']['pluv'],
        driver_contributions['noSLR']['fluv'],
        driver_contributions['noSLR']['coast'],
        driver_contributions['noSLR']['cmpd']
    ],
}

# Prepare the data for seaborn
flat_data = []
for category, values in zip(categories, zip(*data.values())):
    for driver, driver_values in zip(data.keys(), values):
        flat_data.extend([(category, driver, value) for value in driver_values])

df = pd.DataFrame(flat_data, columns=['Category', 'Driver', 'Value'])

# Plot
plt.figure(figsize=(12, 6))
sns.swarmplot(data=df, x='Category', y='Value', hue='Driver', dodge=True, palette='Set2')

# Customization
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Zero line for positive/negative values
plt.title('Climate Attribution of Driver Contributions to Compound Flood from TC Idai')
plt.ylabel('Diff in Flood Volume [%]')
plt.xlabel('Flood Type')
plt.legend(title='Driver Contributions', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.show()





# %% 
                      #####################
                    ###### EXAMPLES ########
                      #####################
# Sample data
categories = ['Pluvial', 'Fluvial', 'Coastal', 'Compound']
data = {
    'Pluvial': [np.random.normal(loc=5, scale=2, size=5),
                np.random.normal(loc=4, scale=2, size=5),
                np.random.normal(loc=6, scale=2, size=5),
                np.random.normal(loc=6, scale=2, size=5)],
    'Fluvial': [np.random.normal(loc=-3, scale=1.5, size=5),
                np.random.normal(loc=-4, scale=1.5, size=5),
                np.random.normal(loc=-2, scale=1.5, size=5),
                np.random.normal(loc=-2, scale=1.5, size=5)],
    'Coastal': [np.random.normal(loc=2, scale=1, size=5),
                np.random.normal(loc=3, scale=1, size=5),
                np.random.normal(loc=1, scale=1, size=5),
                np.random.normal(loc=1, scale=1, size=5)],
    'Compound': [np.random.normal(loc=2, scale=1, size=5),
                 np.random.normal(loc=3, scale=1, size=5),
                 np.random.normal(loc=1, scale=1, size=5),
                 np.random.normal(loc=1, scale=1, size=5)],
}

# Create positions for the boxplots
n_sets = len(data)
n_categories = len(categories)
group_width = 1  # Width of each group of boxes
set_spacing = group_width / n_sets
positions = []

for i, category in enumerate(categories):
    base = i * (group_width + 0.5)  # 0.5 adds spacing between categories
    positions.extend([base + j * set_spacing for j in range(n_sets)])

# Flatten the data for plotting
flattened_data = [item for sublist in zip(*data.values()) for item in sublist]

# Plot
fig, ax = plt.subplots(figsize=(10, 6))
box = ax.boxplot(flattened_data, positions=positions, widths=0.2, patch_artist=True)

# Add colors for each set
colors = ['lightblue', 'lightgreen', 'coral', 'gold']  # Extended to include 4 colors
for patch, color in zip(box['boxes'], colors * n_categories):
    patch.set_facecolor(color)

# Custom x-ticks for categories
category_centers = [(i * (group_width + 0.5) + (group_width - set_spacing) / 2) for i in range(n_categories)]
ax.set_xticks(category_centers)
ax.set_xticklabels(categories)

# Labels and grid
ax.set_ylabel('Diff in Flood Volume [%]')
ax.set_title('Boxplot with Grouped Boxes and Spaced Categories')
ax.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Zero line for positive/negative values
ax.grid(axis='y', linestyle='--', alpha=0.7)

# Legend
handles = [plt.Line2D([0], [0], color=color, lw=4, label=key) for key, color in zip(data.keys(), colors)]
ax.legend(handles, data.keys(), title='Data Sets', loc='upper left')

plt.tight_layout()
plt.show()


# %%
# MAKE A VIOLIN PLOT
# Prepare the data for seaborn
flat_data = []
for category, values in zip(categories, zip(*data.values())):
    for dataset_index, dataset_values in enumerate(values):
        flat_data.extend([(category, f"Dataset {dataset_index+1}", value) for value in dataset_values])

df = pd.DataFrame(flat_data, columns=['Category', 'Dataset', 'Value'])

# Plot
plt.figure(figsize=(12, 6))
sns.violinplot(data=df, x='Category', y='Value', hue='Dataset', dodge=True, palette='Set2', inner='point')

# Customization
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Zero line for positive/negative values
plt.title('Climate Attribution Data by Category (Violin Plot)')
plt.ylabel('Diff in Flood Volume [%]')
plt.xlabel('Flood Type')
plt.legend(title='Driver', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
>>>>>>> 09401fe (Adding attribution plotting script)

plt.show()

#%%
# Get unique station IDs
station_ids = df_all["station_id"].unique()
nrows = len(station_ids)  # One row per station

# Create figure and axes with higher resolution and shared x-axis
fig, axes = plt.subplots(nrows=nrows, ncols=1, figsize=(14, 40), sharex=True, constrained_layout=True)

# Ensure axes is iterable even when there's only one subplot
if nrows == 1:
    axes = [axes]

# Define the two models to plot
selected_models = [
    "event_tp_era5_hourly_CF0_GTSMv41opendap_CF0_spw_IBTrACS_CF0",
    "event_tp_era5_hourly_CF0_GTSMv41opendap_CF-0.14_spw_IBTrACS_CF0",
    # "event_tp_era5_hourly_CF0_GTSMv41opendap_CF-0.14_toSFINCSwaterlevel_spw_IBTrACS_CF0"
]

# Plot data for each station
for ax, station_id in zip(axes, station_ids):
    df_station = df_all[df_all["station_id"] == station_id]

    # Filter the DataFrame to include only the selected models
    df_selected_models = df_station[df_station["model_name"].isin(selected_models)]

    # Plot each model separately
    for model_name in selected_models:
        df_model = df_selected_models[df_selected_models["model_name"] == model_name]
        ax.plot(df_model["time"], df_model["point_zs"], label=model_name, linewidth=2)

    ax.set_title(f"Station {station_id}", fontsize=18)
    ax.set_ylabel("Water Level [m]", fontsize=14)
    ax.grid()
    
    # Increase y-axis tick size
    ax.tick_params(axis="y", labelsize=14)

# Shared x-axis label
axes[-1].set_xlabel("Time", fontsize=16)
axes[-1].tick_params(axis="x", rotation=60, labelsize=14)  # Rotate x-ticks for clarity

# Add a single legend at the top
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=16)

plt.show()

# %%
# Plot basemap with observation points 
fig, ax = models[0]['sfincs_model'].plot_basemap(fn_out=None, bmap="sat", figsize=(11, 7))
=======

# Sample data for drivers (precipitation, temperature, etc.)
categories = ['Pluvial', 'Fluvial', 'Coastal', 'Compound']
data = {
    'Pluvial': [np.random.normal(loc=5, scale=2, size=10),
                np.random.normal(loc=4, scale=2, size=10),
                np.random.normal(loc=6, scale=2, size=10),
                np.random.normal(loc=6, scale=2, size=10)],
    'Fluvial': [np.random.normal(loc=-3, scale=1.5, size=10),
                np.random.normal(loc=-4, scale=1.5, size=10),
                np.random.normal(loc=-2, scale=1.5, size=10),
                np.random.normal(loc=-2, scale=1.5, size=10)],
    'Coastal': [np.random.normal(loc=2, scale=1, size=10),
                np.random.normal(loc=3, scale=1, size=10),
                np.random.normal(loc=1, scale=1, size=10),
                np.random.normal(loc=1, scale=1, size=10)],
    'Compound': [np.random.normal(loc=2, scale=1, size=10),
                 np.random.normal(loc=3, scale=1, size=10),
                 np.random.normal(loc=1, scale=1, size=10),
                 np.random.normal(loc=1, scale=1, size=10)],
}

# Extra base scenario data for the Pluvial Precipitation
extra_precipitation_data = np.random.normal(loc=7, scale=2, size=100)  # Extra data for Pluvial + Precipitation

# Prepare the data
flat_data = []
for category, values in zip(categories, zip(*data.values())):
    for dataset_index, dataset_values in enumerate(values):
        if category == 'Pluvial' and dataset_index == 0:  # Targeting Pluvial + Precipitation
            # Add the original and the extra data for Pluvial and Precipitation
            flat_data.extend([(category, f"Precipitation {dataset_index + 1}", value) for value in dataset_values])
            flat_data.extend([(category, f"Precipitation Extra", value) for value in extra_precipitation_data])
        else:
            flat_data.extend([(category, f"Dataset {dataset_index + 1}", value) for value in dataset_values])

df = pd.DataFrame(flat_data, columns=['Category', 'Dataset', 'Value'])

# Plot
plt.figure(figsize=(12, 6))

# Plot the regular violin plots for all categories and drivers
sns.violinplot(data=df, x='Category', y='Value', hue='Dataset', dodge=True, palette='Set2', inner='point', scale='count')

# Overlay the extra violin plot only for Pluvial Precipitation
extra_data = df[df['Dataset'] == 'Precipitation Extra']
sns.violinplot(data=extra_data, x='Category', y='Value', dodge=True, color='lightgrey', inner='point', scale='count', alpha=0.5)

# Customization
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Zero line for positive/negative values
plt.title('Climate Attribution Data by Category (Violin Plot) with Extra Precipitation Data')
plt.ylabel('Diff in Flood Volume [%]')
plt.xlabel('Flood Type')

# Adjust plot appearance
plt.legend(title='Driver', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.show()


#%%
# Example data
data = [2.5, 3.0, 2.8, 3.2, 2.9]

sns.swarmplot(data=[data])
plt.ylabel("Diff in Flood Volume [%]")
plt.title("Bee Swarm Plot for Small Sample Size")
plt.show()

# %%
sns.violinplot(data=[data], inner=None, color="lightblue")
sns.swarmplot(data=[data], color="black", size=7)
plt.ylabel("Values")
plt.title("Violin Plot with Data Points")
plt.show()

# %%

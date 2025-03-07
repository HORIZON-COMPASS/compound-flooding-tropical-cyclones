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

# %% Load the SFINCS model runs
config_general  = '../../Workflows/01_config_snakemake/config_general.yml'
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
verification_points = runname_id['verification_points']
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
    num_CF_diff = sum(1 for v in CF_values if v != '0')

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
    print(f"{model['category']} -> {model['model_name']} -> {model['model_path']} -> CF_info: {model['CF_info']}")


# %%
# datasets = [
#     (mod_F_PluvCoast, sfincs_root_F_PluvCoast, "F_PluvCoast", "Factual", "Pluvial", "precipitation"),
#     (mod_CF_7precip_PluvCoast, sfincs_root_CF_7precip_PluvCoast, "CF_7precip_PluvCoast", "Counterfactual", "Pluvial", "precipitation"),
#     (mod_F_Coast, sfincs_root_F_Coast, "F_SLR_Coast", "Factual", "Coastal", "SLR-noMDT"),
#     (mod_F_Coast_MDT,sfincs_root_F_Coast_MDT, "F_SLR_Coast_MDT", "Factual", "Coastal", "SLR"),
#     (mod_CF_noSLR_Coast_MDT, sfincs_root_CF_noSLR_Coast_MDT, "CF_noSLR_Coast_MDT", "Counterfactual", "Coastal", "noSLR")
# ]

# %%
# read global surface water occurance (GSWO) data to mask permanent water
models[0]["sfincs_model"].data_catalog.from_yml(os.path.join('../../Workflows/03_data_catalogs/datacatalog_general.yml'))
gswo = models[0]["sfincs_model"].data_catalog.get_rasterdataset("gswo", geom=models[0]["sfincs_model"].region, buffer=1000)
# %%
### loop over the different sfincs_model to compute hmax ###
# we set a threshold to mask minimum flood depth
hmin = 0.05

for model in models:
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
    ax.set_title(name, fontsize=12)

# Add colorbar and label
fig.colorbar(im, ax=axes, orientation="vertical", fraction=0.02, pad=0.04).set_label('Masked hmax (m)', rotation=270, labelpad=15)

# Show the plot
plt.show()

#%% Create plots of the difference between F and CF scenarios
# Calculate the differences
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
dx = abs(mod_F_PluvCoast.results['hmax_masked'].x[1] - mod_F_PluvCoast.results['hmax_masked'].x[0])  # Grid resolution in x-direction (meters)
dy = abs(mod_F_PluvCoast.results['hmax_masked'].y[1] - mod_F_PluvCoast.results['hmax_masked'].y[0])  # Grid resolution in y-direction (meters)
cell_area = dx * dy 

# %% Calculate the area of the flooded cells as flood extent with hmax_masked as variable
flood_characs = []

for (data, sfincs_root, name, scenario, flood_type, driver) in (datasets):
    flooded_cells = data.results['hmax_masked'] > 0  # Create a boolean mask
    flood_extent = (flooded_cells * cell_area).sum().compute()  # Compute the total flooded area

    # Convert the result to square kilometers
    flood_extent_km2 = flood_extent / 1e6  # Convert square meters to square kilometers
    
    # Store the result in the dictionary
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

plt.show()

#%%

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

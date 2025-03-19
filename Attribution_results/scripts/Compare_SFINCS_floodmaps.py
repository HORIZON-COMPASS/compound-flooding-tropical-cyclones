#%%
# Load the necessary packages
import os
from os.path import join
import yaml
import matplotlib.pyplot as plt
import itertools
import math
from hydromt_sfincs import SfincsModel, utils
from hydromt import DataCatalog
import contextily as ctx
import cartopy.crs as ccrs
import pandas as pd
import seaborn as sns
from pyproj import Transformer

# Function to load YAML configuration file
def load_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)


# Load the SFINCS models and create a model dictonary
def load_sfincs_models(config):
    """Generates model paths and categories for SFINCS runs based on CF values."""
    run = config['runname_ids']['Idai']
    base_path = f"p:/11210471-001-compass/03_Runs/{run['region']}/{run['tc_name']}/sfincs/"
    
    models = []
    for rain, wind, slr in itertools.product(run['CF_value_rain'], run['CF_value_wind'], run['CF_value_SLR']):
        model_name = f"event_tp_{run['precip_forcing']}_CF{rain}_{run['tidemodel']}_CF{slr}_{run['wind_forcing']}_CF{wind}"
        model_path = os.path.join(base_path, model_name)
        utmzone    = run['utmzone']
        model_obj  = SfincsModel(model_path, mode="r")

        num_CF_diff      = sum(v != 0 for v in [rain, wind, slr])
        categories       = ["Factual", "Single Driver Counterfactual", "Counterfactual Driver Pair", "Counterfactual Compound Driver"]
        short_cats       = ["F", "CF_DR_single", "CF_DR_pair", "CF_DR_compound"]
        CF_info          = {k: v for k, v in zip(["rain", "wind", "SLR"], [rain, wind, slr]) if v != 0}
        non_zero_CF_info = {k: v for k, v in CF_info.items() if v != 0}
        CF_info_str      = ", ".join(f"{k}: {v}" for k, v in non_zero_CF_info.items())

        models.append({
            "model_name":     model_name,
            "model_path":     model_path,
            "utmzone":        utmzone,
            "sfincs_model":   model_obj,
            "sfincs_results": model_obj.results,
            "category":       categories[num_CF_diff],
            "cat_short":      short_cats[num_CF_diff],
            "CF_info":        CF_info,
            "CF_info_str":    CF_info_str
        })

    return models


# read global surface water occurance (GSWO) data to mask permanent water for the model region
def gwso_sfincs_region(model):
    datacat_path = os.path.abspath("../../Workflows/03_data_catalogs/datacatalog_general.yml")
    data_catalog = DataCatalog(data_libs = [datacat_path])
    sfincs_region = model["sfincs_model"].region
    gwso_region = data_catalog.get_rasterdataset("gswo", geom=sfincs_region, buffer=1000)
    return gwso_region


# Compute the maximum water level (hmax) and mask out permanent water
def compute_hmax_masked(models, gwso_region):
    # we set a threshold to mask minimum flood depth
    hmin = 0.05

    for model in models:
        # select the highest-resolution elevation dataset
        print(f"Processing model: {model['model_name']}")
        depfile = join(model["model_path"], "subgrid", "dep_subgrid.tif")
        da_dep = model["sfincs_model"].data_catalog.get_rasterdataset(depfile)

        # compute the maximum over all time steps
        da_zsmax = model["sfincs_results"]["zsmax"].max(dim="timemax")
     
        # downscale the floodmap
        da_hmax = utils.downscale_floodmap(
            zsmax=da_zsmax,
            dep=da_dep,
            hmin=hmin,
            # floodmap_fn=join(sfincs_root, "gis/floodmap.tif") # uncomment to save floodmap to <mod.root>/floodmap.tif
            )
    
        # GSWO dataset to mask permanent water by first geprojecting it to the subgrid of hmax
        gswo_mask = gwso_region.raster.reproject_like(da_hmax, method="max")
        # permanent water where water occurence > 5%
        da_hmax_masked = da_hmax.where(gswo_mask <= 5)

        # Add the name attribute for identification
        model["sfincs_results"]['hmax'] = da_hmax
        model["sfincs_results"]['hmax_masked'] = da_hmax_masked

        # Open the existing NetCDF dataset in append mode
        # nc_file = join(model["model_path"], "sfincs_his.nc")
        
        # try:
        #     with xr.open_dataset(nc_file, mode="a") as ds:  
        #         # Convert DataArrays to Dataset with proper variable names
        #         ds_new = xr.Dataset({
        #             "hmax": da_hmax,
        #             "hmax_masked": da_hmax_masked
        #         })
                
        #         # Merge with existing dataset and save
        #         ds_updated = xr.merge([ds, ds_new])
        #         ds.close()

        #         ds_updated.to_netcdf(nc_file, mode="w")  # Overwrite with new data

        #         print(f"Saved hmax and hmax_masked to {nc_file}")

        # except Exception as e:
        #     print(f"Error saving to {nc_file}: {e}")

    return models


# Plot the maximum water spatially
def plot_masked_hmax(models, num_cols=2, figsize=(8, 10)):
    """
    Creates a series of subplots based on the given models, plotting the masked hmax values.

    Parameters:
    - models: List of model dictionaries, each containing 'sfincs_model' and 'category' data.
    - num_cols: Number of columns for the subplots. Default is 2.
    - figsize: Tuple for the figure size. Default is (8, 10).
    """
    # Determine the projection from the first model
    projection = models[0]['sfincs_model'].crs.to_epsg()

    # Determine subplot grid size
    num_models = len(models)
    num_rows = math.ceil(num_models / num_cols)

    # Create the subplots
    fig, axes = plt.subplots(nrows=num_rows, 
                             ncols=num_cols, 
                             figsize=figsize, 
                             constrained_layout=True, 
                             subplot_kw={"projection": ccrs.epsg(projection)})

    # Add a super title
    fig.suptitle("Masked hmax", fontsize=12, y=1.02)

    # Loop through datasets and plot
    for model, ax in zip(models, axes.flatten()):
        # Plot the masked water depth and add basemap
        im = model['sfincs_results']['hmax_masked'].plot.pcolormesh(
            ax=ax, cmap="Blues", vmin=0, vmax=3.0, add_colorbar=False
        )
        
        # Add basemap
        ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=12, crs=model['sfincs_results']['hmax_masked'].rio.crs, attribution=False)
        
        # Set axis labels and title
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

        ax.set_title(f"{model['cat_short']} ({model['CF_info_str']})", fontsize=11)

    # Add colorbar and label
    fig.colorbar(im, ax=axes, orientation="vertical", fraction=0.02, pad=0.04).set_label('Masked hmax (m)', rotation=270, labelpad=15)

    # Show the plot
    plt.show()


# Compute the differences between the Factual and Counterfactual masked hmax variables
def compute_hmax_diff(models):
    factual_hmax = None
    
    # Check if "Factual" category exists
    for model in models:
        if model["category"] == "Factual":
            factual_hmax = model['sfincs_results'].get("hmax_masked", None)
            if factual_hmax is None:
                print(f"Error: 'hmax_masked' not found for factual model: {model['model_name']}")
                return models  # Exit early if 'hmax_masked' is missing

    # Compute difference for counterfactual models
    for model in models:
        if model["category"] != "Factual":
            hmax_masked = model['sfincs_results'].get("hmax_masked", None)
            if hmax_masked is not None:
                hmax_diff = hmax_masked - factual_hmax
                model['sfincs_results']["hmax_diff"] = hmax_diff
                print(f"hmax_diff calculated for {model['model_name']}")
            else:
                print(f"Warning: 'hmax_masked' not found for counterfactual model: {model['model_name']}")
    
    return models

# Function to plot hmax_diff for counterfactual models
def plot_hmax_diff(models, num_cols=2, zoom_region_latlon=None):
    """
    Plots the hmax_diff for counterfactual models with an optional zoom feature.
    Each subplot has its own x and y labels, and all subplots retain ticks.

    Parameters:
        models (list): List of model dictionaries.
        num_cols (int): Number of columns in the subplot grid.
        figsize (tuple): Size of the figure.
        zoom_region_latlon (tuple, optional): (lat_min, lat_max, lon_min, lon_max) in latitude/longitude.
    """
    # Filter out only counterfactual models with hmax_diff
    counterfactual_models = [m for m in models if "hmax_diff" in m["sfincs_results"]]
    num_models = len(counterfactual_models)

    # Get projection from the first model (assuming all models have the same CRS)
    model_crs = models[0]['sfincs_model'].crs
    model_epsg = model_crs.to_epsg()

    # Convert lat/lon zoom to the model’s coordinate system
    zoom_region = None
    if zoom_region_latlon:
        lat_min, lat_max, lon_min, lon_max = zoom_region_latlon
        transformer = Transformer.from_crs("EPSG:4326", model_epsg, always_xy=True)
        xmin, ymin = transformer.transform(lon_min, lat_min)
        xmax, ymax = transformer.transform(lon_max, lat_max)
        zoom_region = (xmin, xmax, ymin, ymax)

    # Define grid dimensions
    num_rows = math.ceil(num_models / num_cols)

    # Create figure and subplots
    fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, 
                             figsize=(num_cols * 5, num_rows * 6), 
                             constrained_layout=True,
                             subplot_kw={"projection": ccrs.epsg(model_epsg)})

    # Flatten axes array for easy indexing (handles both single and multiple rows)
    if num_models == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    # Define colormap settings
    cmap = "RdBu"
    vmin, vmax = -0.3, 0.3

    # Loop through counterfactual models and plot
    for idx, model in enumerate(counterfactual_models):
        ax = axes[idx]
        hmax_diff = model['sfincs_results'].get("hmax_diff", None)

        if hmax_diff is not None:
            # Plot difference
            im = hmax_diff.plot.pcolormesh(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False)

            # Add basemap
            ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=9, 
                            crs=model_crs, attribution=False)

            # Apply zoom if region is specified
            if zoom_region:
                xmin, xmax, ymin, ymax = zoom_region
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)

            # Construct title with CF info
            non_zero_CF_info = {key: value for key, value in model["CF_info"].items() if value != 0}
            cf_info_str = ", ".join(f"{key}: {value}" for key, value in non_zero_CF_info.items())
            ax.set_title(f"{model['cat_short']} ({cf_info_str})", fontsize=12)

            # Set individual x and y labels for each subplot
            ax.set_xlabel(f"x coordinate UTM zone {model['utmzone']} [m]", fontsize=10)
            ax.set_ylabel(f"y coordinate UTM zone {model['utmzone']} [m]", fontsize=10)
            ax.xaxis.set_visible(True)
            ax.yaxis.set_visible(True)

            # Ensure tick labels are visible
            ax.tick_params(axis="both", labelsize=9)
        else:
            print(f"Error: 'hmax_diff' not found for counterfactual model: {model['model_name']}")

    # Hide any unused subplots (if the grid is larger than the number of models)
    for idx in range(num_models, len(axes)):
        fig.delaxes(axes[idx])

    # Add shared colorbar
    cbar = fig.colorbar(im, ax=axes, orientation="vertical", fraction=0.02, pad=0.02)
    cbar.set_label('Difference in Maximum Water Level (m)', rotation=270, labelpad=15, fontsize=12)
    cbar.ax.tick_params(labelsize=10)

    # Reduce whitespace between subplots
    fig.subplots_adjust(wspace=0.02, hspace=0.15)

    # Show the plot
    plt.show()


# Function to calculate the surface area of one grid cell
def calculate_cell_area(model):
    # Error handling for missing 'hmax_masked'
    hmax_masked = model['sfincs_results'].get('hmax_masked', None)
    if hmax_masked is None:
        print(f"Error: 'hmax_masked' not found for model: {model['model_name']}")
        return None  # Return None if no data is available

    # Calculate the cell area
    dx = abs(hmax_masked.x[1] - hmax_masked.x[0])  # Grid resolution in x-direction (meters)
    dy = abs(hmax_masked.y[1] - hmax_masked.y[0])  # Grid resolution in y-direction (meters)
    return dx * dy  # Area of one grid cell (m²)


# Function to calculate the flooded area (extent) for each model
def calculate_flood_extent(models):
    for model in models:
        # Error handling for missing 'hmax_masked'
        hmax_masked = model['sfincs_results'].get('hmax_masked', None)
        if hmax_masked is None:
            print(f"Error: 'hmax_masked' not found for model: {model['model_name']}")
            continue  # Skip this model if 'hmax_masked' is missing

        # Create a boolean mask for flooded cells (hmax_masked > 0)
        flooded_cells = hmax_masked > 0
        # Compute the total flooded area (in square meters)
        flood_extent = (flooded_cells * calculate_cell_area(model)).sum().compute()
        # Convert to square kilometers
        flood_extent_km2 = flood_extent / 1e6
        # Store in the model dictionary
        model['sfincs_results']['flood_extent_km2'] = flood_extent_km2
        print(f"for model {model['model_name']}, the flooded area: {flood_extent_km2}")
    return models


# Function to calculate the flooded volume (depth * area) for each model
def calculate_flood_volume(models):
    for model in models:
        # Error handling for missing 'hmax_masked'
        hmax_masked = model['sfincs_results'].get('hmax_masked', None)
        if hmax_masked is None:
            print(f"Error: 'hmax_masked' not found for model: {model['model_name']}")
            continue  # Skip this model if 'hmax_masked' is missing

        # Create a boolean mask for flooded cells (hmax_masked > 0)
        flooded_cells = hmax_masked > 0
        # Compute the flooded volume (sum of depth * area for each flooded cell)
        flood_volume = (hmax_masked * flooded_cells * calculate_cell_area(model)).sum().compute()
        # Convert to cubic kilometers
        flood_volume_km3 = flood_volume / 1e9
        # Store in the model dictionary
        model['sfincs_results']['flood_volume_km3'] = flood_volume_km3
        print(f"for model {model['model_name']}, the flood volume is {flood_volume_km3}")
    return models


# Function to calculate the flood volume and extent differences between factual and counterfactual datasets
def calculate_flood_differences(models):
    factual_flood_volume = None
    factual_flood_extent = None
    
    for model in models:
        # Store factual flood volume and extent for comparison
        if model["category"] == "Factual":
            factual_flood_volume = model['sfincs_results']['flood_volume_km3']
            factual_flood_extent = model['sfincs_results']['flood_extent_km2']

            if factual_flood_volume is None:
                print(f"Error: 'flood_volume_km3' not found for factual model: {model['model_name']}")
                continue  # Skip this model if factual flood volume is missing
            
            if factual_flood_extent is None:
                print(f"Error: 'flood_extent_km2' not found for factual model: {model['model_name']}")
                continue  # Skip this model if factual flood extent is missing

        # Compute flood volume difference for counterfactual models
        if factual_flood_volume is not None and model["category"] != "Factual":
            # Error handling for missing flood extent in counterfactual models
            flood_volume = model['sfincs_results'].get('flood_volume_km3', None)
            if flood_volume is None:
                print(f"Error: 'flood_volume_km3' not found for counterfactual model: {model['model_name']}")
                continue  # Skip this model if flood volume is missing
            
            flood_volume_diff = (factual_flood_volume - flood_volume) / factual_flood_volume * 100
            model['sfincs_results']['Volume_diff_from_F(%)'] = flood_volume_diff
            print(f"flood_volume_diff calculated for {model['model_name']}")

        # Compute flood extent difference for counterfactual models
        if factual_flood_extent is not None and model["category"] != "Factual":
            # Error handling for missing flood extent in counterfactual models
            flood_extent = model['sfincs_results'].get('flood_extent_km2', None)
            if flood_extent is None:
                print(f"Error: 'flood_extent_km2' not found for counterfactual model: {model['model_name']}")
                continue  # Skip this model if flood extent is missing

            flood_extent_diff = (factual_flood_extent - flood_extent) / factual_flood_extent * 100
            model['sfincs_results']['Extent_diff_from_F(%)'] = flood_extent_diff
            print(f"flood_extent_diff calculated for {model['model_name']}")

    return models


# Function to compare flood characteristics between factual and counterfactual models
def compare_flood_characteristics(models):
    results = {}
    factual_area, factual_volume = None, None

    # Find the factual model and calculate its flood characteristics
    for model in models:
        if model["category"] == "Factual":
            factual_area, factual_volume = calculate_flood_characteristics(model)
            break  # No need to continue once factual model is found

    if factual_area is None or factual_volume is None:
        print("Error: No factual model found or missing flood characteristics")
        return None  # Exit if no factual model found

    # Loop through counterfactual models and compute differences
    for model in models:
        if model["category"] != "Factual":
            cf_area, cf_volume = calculate_flood_characteristics(model)
            area_diff = cf_area - factual_area
            volume_diff = cf_volume - factual_volume

            model['sfincs_results']['Extent_diff_from_F(%)'] = area_diff
            model['sfincs_results']['Volume_diff_from_F(%)'] = volume_diff
                
    return models


# Function to create a categorical plot comparing flood volume differences by CF driver
def plot_flood_volume_difference_by_driver(models):
    # Create lists for model names, volume differences, and drivers
    model_names = []
    volume_diffs = []
    drivers = []

    # Collect data from models excluding "Factual"
    for model in models:
        if model['category'] != "Factual":
            model_names.append(model["model_name"])
            volume_diff_value = model['sfincs_results'].get("Volume_diff_from_F(%)", None)
            if volume_diff_value is not None:
                volume_diffs.append(volume_diff_value.values.flatten()[0])

            CF_info = model["CF_info"]
            if all(k in CF_info and CF_info[k] != 0 for k in ["rain", "wind", "SLR"]):
                drivers.append("Compound")
            elif "rain" in CF_info and "wind" in CF_info and CF_info["rain"] != 0 and CF_info["wind"] != 0:
                drivers.append("Rain & wind")
            elif "rain" in CF_info and "SLR" in CF_info and CF_info["rain"] != 0 and CF_info["SLR"] != 0:
                drivers.append("Rain & SLR")
            elif "wind" in CF_info and "SLR" in CF_info and CF_info["wind"] != 0 and CF_info["SLR"] != 0:
                drivers.append("Wind & SLR")
            elif "rain" in CF_info and CF_info["rain"] != 0:
                drivers.append("Rain")
            elif "wind" in CF_info and CF_info["wind"] != 0:
                drivers.append("Wind")
            elif "SLR" in CF_info and CF_info["SLR"] != 0:
                drivers.append("SLR")
            print(drivers)

    # Create a DataFrame for plotting
    df = pd.DataFrame({
        "Model name": model_names,
        "Flood volume difference (F-CF in %)": volume_diffs,
        "CF driver": drivers
    })

    # Create the seaborn categorical plot (swarm plot)
    sns.set(style="whitegrid")
    catplot = sns.catplot(data=df, 
                          x="CF driver", 
                          y="Flood volume difference (F-CF in %)", 
                          kind="swarm", 
                          height=6, 
                          aspect=2)
    catplot.set_xticklabels(rotation=45)  # Rotate x-axis labels for better readability
    plt.title("Flood Volume Difference by CF Driver")
    plt.tight_layout()
    plt.show()


##############################################################
# Use of functions
#%%
# Load snakemake config file to construct the model paths
config_path  = '../../Workflows/01_config_snakemake/config_general_MZB.yml'
cfg = load_config(config_path)

#%%
# Load the SFINCS models in one dictonary
models = load_sfincs_models(cfg)
print(models)

#%%
# Calculate hmax and mask out permanent water
gwso = gwso_sfincs_region(models[0])
models = compute_hmax_masked(models, gwso)
models = compute_hmax_diff(models)

#%%
# Calculate flood characteristics 
models = calculate_flood_extent(models)
models = calculate_flood_volume(models)
models = calculate_flood_differences(models)

#%%
# SOME PLOTTING
plot_hmax_diff(models)
plot_masked_hmax(models)
plot_flood_volume_difference_by_driver(models)

# %%
zoom_region = (-19.95, -19.6, 34.55, 35.1)
plot_hmax_diff(models, zoom_region_latlon=zoom_region)

# %%

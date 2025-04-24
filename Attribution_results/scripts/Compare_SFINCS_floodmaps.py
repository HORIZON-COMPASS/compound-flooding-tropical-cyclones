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
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.patches import Wedge
import xarray as xr
import geopandas as gpd
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
    factual_model = None  # Initialize factual_model before the loop
    
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

        # Create model dictionary
        model_dict = {
            "model_name": model_name,
            "model_path": model_path,
            "utmzone": utmzone,
            "sfincs_model": model_obj,
            "sfincs_results": model_obj.results,
            "category": categories[num_CF_diff],
            "cat_short": short_cats[num_CF_diff],
            "CF_info": CF_info,
            "CF_info_str": CF_info_str
        }

        # If the model is factual (CF values are all 0), assign it to factual_model
        if num_CF_diff == 0:
            factual_model = model_dict
        else:
            models.append(model_dict)

    # Ensure the factual model is the first one in the list
    if factual_model is not None:
        models.insert(0, factual_model)

    return models


def load_fiat_models(config):
    run = config['runname_ids']['Idai']
    base_path = f"p:/11210471-001-compass/03_Runs/{run['region']}/{run['tc_name']}/fiat/"
    
    models = []
    factual_model = None  # Initialize factual_model before the loop
    
    for rain, wind, slr in itertools.product(run['CF_value_rain'], run['CF_value_wind'], run['CF_value_SLR']):
        model_name = f"event_tp_{run['precip_forcing']}_CF{rain}_{run['tidemodel']}_CF{slr}_{run['wind_forcing']}_CF{wind}"
        model_path = os.path.join(base_path, model_name)
        fiat_results = pd.read_csv(join(f"{model_path}", "output/output.csv"))
        fiat_results_spatial =  gpd.read_file(join(f"{model_path}", "output/spatial.fgb"))

        num_CF_diff      = sum(v != 0 for v in [rain, wind, slr])
        categories       = ["Factual", "Single Driver Counterfactual", "Counterfactual Driver Pair", "Counterfactual Compound Driver"]
        short_cats       = ["F", "CF_DR_single", "CF_DR_pair", "CF_DR_compound"]
        CF_info          = {k: v for k, v in zip(["rain", "wind", "SLR"], [rain, wind, slr]) if v != 0}
        non_zero_CF_info = {k: v for k, v in CF_info.items() if v != 0}
        CF_info_str      = ", ".join(f"{k}: {v}" for k, v in non_zero_CF_info.items())

        model_dict = {
            "model_name":             model_name,
            "model_path":             model_path,
            "fiat_results":           fiat_results,
            "fiat_results_spatial":   fiat_results_spatial,
            "category":               categories[num_CF_diff],
            "cat_short":              short_cats[num_CF_diff],
            "CF_info":                CF_info,
            "CF_info_str":            CF_info_str
        }

        # If the model is factual (CF values are all 0), assign it to factual_model
        if num_CF_diff == 0:
            factual_model = model_dict
        else:
            models.append(model_dict)

    # Ensure the factual model is the first one in the list
    if factual_model is not None:
        models.insert(0, factual_model)

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
def plot_hmax_diff(models, num_cols=1, zoom_region_latlon=None):
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
    if num_models > 1:
        num_cols = 2
    
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


def zoom_in_sfincs_data(models, zoom_region_latlon=None):
    """
    Filters each model's sfincs_results data to retain only the region within the zoom extent.

    Parameters:
        models (list): List of model dictionaries.
        zoom_region_latlon (tuple, optional): (lat_min, lat_max, lon_min, lon_max)
    """
    if zoom_region_latlon:
        lat_min, lat_max, lon_min, lon_max = zoom_region_latlon

        # Get CRS and transformer
        model_crs = models[0]['sfincs_model'].crs
        model_epsg = model_crs.to_epsg()
        transformer = Transformer.from_crs("EPSG:4326", model_epsg, always_xy=True)

        # Transform lat/lon to model coordinates
        xmin, ymin = transformer.transform(lon_min, lat_min)
        xmax, ymax = transformer.transform(lon_max, lat_max)

        for model in models:
            sfincs_results = model['sfincs_results']

            for key, val in sfincs_results.items():
                if isinstance(val, xr.DataArray) and {'x', 'y'}.issubset(val.dims):
                    sfincs_results[key] = val.sel(x=slice(xmin, xmax), y=slice(ymin, ymax))

            print(f"Model {model['model_name']} zoomed to region {zoom_region_latlon}")
    
    else:
        print("No zoom region specified. Returning unzoomed models.")
    
    return models


def load_and_zoom_fiat_fgb(fiat_models, zoom_region_latlon=None):
    """
    Filter spatial FIAT results within a zoom region for all fiat_models.
    
    Parameters:
        fiat_models (list): List of fiat model dictionaries with a 'fiat_results_spatial' GeoDataFrame.
        zoom_region_latlon (tuple): (lat_min, lat_max, lon_min, lon_max)
    
    Returns:
        list: Updated fiat_models list with 'fiat_results_spatial' cropped to the zoom region.
    """
    if zoom_region_latlon:
        lat_min, lat_max, lon_min, lon_max = zoom_region_latlon

        for fiat in fiat_models:
            gdf = fiat['fiat_results_spatial']

            # Ensure CRS is WGS84 (lat/lon)
            if gdf.crs != "EPSG:4326":
                gdf = gdf.to_crs("EPSG:4326")

            # Spatial filter using bounds
            gdf_filtered = gdf.cx[lon_min:lon_max, lat_min:lat_max]
            print(f"[{fiat['model_name']}] Filtered to {len(gdf_filtered)} rows from original {len(gdf)}")

            # Save the filtered GeoDataFrame back
            fiat['fiat_results_spatial'] = gdf_filtered

    else:
        print("No zoom region specified. Returning unfiltered fiat models.")

    return fiat_models


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


def calculate_damage_differences(fiat_models):
    factual_total_damage = None  # Variable to store the factual total damage

    for model in fiat_models:
        # Store factual total damage for comparison
        if model["category"] == "Factual":
            factual_total_damage = model['fiat_results'].get("total_damage", None)
            factual_total_damage_spatial = model['fiat_results_spatial'].get("total_damage", None)

            if factual_total_damage is None:
                print(f"Error: 'total_damage' not found for factual model: {model['model_name']}")
                continue  # Skip this model if factual total damage is missing

        # Compute damage difference for counterfactual models
        if factual_total_damage is not None and model["category"] != "Factual":
            total_damage = model['fiat_results'].get('total_damage', None)
            if total_damage is None:
                print(f"Error: 'total_damage' not found for counterfactual model: {model['model_name']}")
                continue  # Skip this model if total damage is missing
            
            # Calculate the damage difference from the factual model (in percentage)
            damage_diff = (factual_total_damage - total_damage) / factual_total_damage * 100
            model['fiat_results']['Damage_diff_from_F(%)'] = damage_diff
            print(f"damage_diff calculated for {model['model_name']}")
        
        if factual_total_damage_spatial is not None and model["category"] != "Factual":
            total_spatial_damage = model['fiat_results_spatial'].get('total_damage', None)
            if total_spatial_damage is None:
                print(f"Error: 'total_damage' not found for counterfactual model: {model['model_name']}")
                continue  # Skip this model if total damage is missing
            
            # Calculate the damage difference from the factual model (in percentage)
            damage_diff = (factual_total_damage_spatial - total_spatial_damage) / factual_total_damage_spatial * 100
            model['fiat_results_spatial']['Damage_diff_from_F(%)'] = damage_diff
            print(f"damage_diff calculated for {model['model_name']}")


    return fiat_models


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
    print(len(model_names))
    print(len(volume_diffs))
    print(len(drivers))
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


def plot_flood_and_damage_difference_by_driver(sfincs_models, fiat_models):
    # Create lists for model names, volume differences, damage differences, and drivers
    model_names = []
    volume_diffs = []
    damage_diffs = []
    drivers = []

    # Collect data from models excluding "Factual"
    for sfincs, fiat in zip(sfincs_models, fiat_models):
        if sfincs['category'] != "Factual":
            model_names.append(sfincs["model_name"])
            
            # Flood volume difference
            volume_diff_value = sfincs['sfincs_results'].get("Volume_diff_from_F(%)", None)
            if volume_diff_value is not None:
                volume_diffs.append(volume_diff_value.values.flatten()[0])
            
            # Damage difference (from fiat results)
            mean_damage_diff = fiat['fiat_results']['Damage_diff_from_F(%)'][np.isfinite(fiat['fiat_results']['Damage_diff_from_F(%)'])].mean()
            if mean_damage_diff is not None:
                damage_diffs.append(mean_damage_diff)
            
            # Determine CF driver
            CF_info = sfincs["CF_info"]
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

    # Create a DataFrame for plotting
    df = pd.DataFrame({
        "Model name": model_names,
        "Flood volume difference (F-CF in %)": volume_diffs,
        "Mean damage difference (F-CF in %)": damage_diffs,
        "CF driver": drivers
    })

    # Create the plot
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot flood volume difference
    scatter1 = ax1.scatter(df['CF driver'], df['Flood volume difference (F-CF in %)'], 
                           color='b', label="Flood Volume Difference", zorder=2)
    ax1.set_xlabel("Counterfactual Flood Driver", fontsize=14, labelpad=15)  # Increase font size and space
    ax1.set_ylabel("Flood volume difference (F-CF in %)", fontsize=14, color='b')  # Font size and color

    # Set the ax1 color to match the scatter points
    ax1.tick_params(axis='y', labelcolor='b', labelsize=12)  # Increase tick label size
    ax1.set_ylabel("Flood volume difference (F-CF in %)", color='b')

    # Create a second y-axis for damage difference
    ax2 = ax1.twinx()

    # Plot damage difference
    scatter2 = ax2.scatter(df['CF driver'], df['Mean damage difference (F-CF in %)'], 
                           color='r', label="Damage Difference", zorder=1, alpha=0.7, marker='o')
    ax2.set_ylabel("Damage difference (%)", fontsize=14, color='r')  # Font size and color

    # Set the ax2 color to match the scatter points
    ax2.tick_params(axis='y', labelcolor='r', labelsize=12)  # Increase tick label size
    ax2.set_ylabel("Damage difference (%)", color='r')

    # Rotate x-axis labels for better readability
    ax1.set_xticklabels(df["CF driver"].unique(), fontsize=12)

    # Calculate the y-limits for both axes based on data and add some padding
    padding = 0.25  # You can adjust this padding value as needed

    # For ax1 (Flood Volume Difference)
    ax1_min = df['Flood volume difference (F-CF in %)'].min()
    ax1_max = df['Flood volume difference (F-CF in %)'].max()
    
    # For ax2 (Damage Difference)
    ax2_min = df['Mean damage difference (F-CF in %)'].min()
    ax2_max = df['Mean damage difference (F-CF in %)'].max()

    # Set the same y-limits for both axes and ensure 0 is included
    y_min = min(ax1_min, ax2_min, 0) - padding
    y_max = max(ax1_max, ax2_max, 0) + padding

    # Set y-limits to include the same range for both axes
    ax1.set_ylim([y_min, y_max])
    ax2.set_ylim([y_min, y_max])

    # Add black horizontal line at zero
    ax2.axhline(0, color='k', linewidth=0.2, zorder=0)

    # Title and layout
    plt.title("Flood Volume and Damage Difference by CF Driver", fontsize=16)  # Title font size
    plt.tight_layout()

    # Synchronize the gridlines between the two axes
    ax1.grid(True)  # Enable grid on ax1
    ax2.grid(False)  # Disable grid on ax2
    ax1.set_axisbelow(True)  # Ensure gridlines are below the scatter points

    # Show the plot
    plt.show()


def plot_extent_and_damage_difference_by_driver(sfincs_models, fiat_models):
    # Create lists for model names, volume differences, damage differences, and drivers
    model_names = []
    extent_diffs = []
    damage_diffs = []
    drivers = []

    # Collect data from models excluding "Factual"
    for sfincs, fiat in zip(sfincs_models, fiat_models):
        if sfincs['category'] != "Factual":
            model_names.append(sfincs["model_name"])
            
            # Flood volume difference
            extent_diff_value = sfincs['sfincs_results'].get("Extent_diff_from_F(%)", None)
            if extent_diff_value is not None:
                extent_diffs.append(extent_diff_value.values.flatten()[0])
            
            # Damage difference (from fiat results)
            mean_damage_diff = fiat['fiat_results']['Damage_diff_from_F(%)'][np.isfinite(fiat['fiat_results']['Damage_diff_from_F(%)'])].mean()
            if mean_damage_diff is not None:
                damage_diffs.append(mean_damage_diff)
            
            # Determine CF driver
            CF_info = sfincs["CF_info"]
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

    # Create a DataFrame for plotting
    df = pd.DataFrame({
        "Model name": model_names,
        "Flood extent difference (F-CF in %)": extent_diffs,
        "Mean damage difference (F-CF in %)": damage_diffs,
        "CF driver": drivers
    })

    # Create the plot
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot flood volume difference
    scatter1 = ax1.scatter(df['CF driver'], df['Flood extent difference (F-CF in %)'], 
                           color='b', label="Flood Extent Difference", zorder=2)
    ax1.set_xlabel("Counterfactual Flood Driver", fontsize=14, labelpad=15)  # Increase font size and space
    ax1.set_ylabel("Flood extent difference (F-CF in %)", fontsize=14, color='b')  # Font size and color

    # Set the ax1 color to match the scatter points
    ax1.tick_params(axis='y', labelcolor='b', labelsize=12)  # Increase tick label size
    # ax1.set_ylabel("Flood extent difference (F-CF in %)", color='b')

    # Create a second y-axis for damage difference
    ax2 = ax1.twinx()

    # Plot damage difference
    scatter2 = ax2.scatter(df['CF driver'], df['Mean damage difference (F-CF in %)'], 
                           color='r', label="Damage Difference", zorder=1, alpha=0.7, marker='o')
    ax2.set_ylabel("Damage difference (%)", fontsize=14, color='r')  # Font size and color

    # Set the ax2 color to match the scatter points
    ax2.tick_params(axis='y', labelcolor='r', labelsize=12)  # Increase tick label size
    # ax2.set_ylabel("Damage difference (%)", color='r')

    # Rotate x-axis labels for better readability
    ax1.set_xticklabels(df["CF driver"].unique(), fontsize=12)

    # Calculate the y-limits for both axes based on data and add some padding
    padding = 0.25  # You can adjust this padding value as needed

    # For ax1 (Flood Volume Difference)
    ax1_min = df['Flood extent difference (F-CF in %)'].min()
    ax1_max = df['Flood extent difference (F-CF in %)'].max()
    
    # For ax2 (Damage Difference)
    ax2_min = df['Mean damage difference (F-CF in %)'].min()
    ax2_max = df['Mean damage difference (F-CF in %)'].max()

    # Set the same y-limits for both axes and ensure 0 is included
    y_min = min(ax1_min, ax2_min, 0) - padding
    y_max = max(ax1_max, ax2_max, 0) + padding

    # Set y-limits to include the same range for both axes
    ax1.set_ylim([y_min, y_max])
    ax2.set_ylim([y_min, y_max])

    # Add black horizontal line at zero
    ax2.axhline(0, color='k', linewidth=0.2, zorder=0)

    # Title and layout
    plt.title("Flood Extent and Damage Difference by CF Driver", fontsize=16)  # Title font size
    plt.tight_layout()

    # Synchronize the gridlines between the two axes
    ax1.grid(True)  # Enable grid on ax1
    ax2.grid(False)  # Disable grid on ax2
    ax1.set_axisbelow(True)  # Ensure gridlines are below the scatter points

    # Show the plot
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
#%%
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
fiat_models = load_fiat_models(cfg)
#%%
fiat_models = calculate_damage_differences(fiat_models)

#%%
plot_flood_and_damage_difference_by_driver(models, fiat_models)
plot_extent_and_damage_difference_by_driver(models, fiat_models)
#%%
urban_models = zoom_in_sfincs_data(models, zoom_region_latlon=zoom_region)

#%%
urban_fiat = load_and_zoom_fiat_fgb(fiat_models, zoom_region_latlon=zoom_region)
#%%
plot_flood_and_damage_difference_by_driver(urban_models, urban_fiat)
plot_extent_and_damage_difference_by_driver(urban_models, urban_fiat)

#%%





def plot_flood_and_damage_difference_by_driver_split(sfincs_models, fiat_models):
    
    model_names = []
    volume_diffs = []
    damage_diffs = []
    drivers = []
    driver_sets = []

    for sfincs, fiat in zip(sfincs_models, fiat_models):
        if sfincs['category'] != "Factual":
            model_names.append(sfincs["model_name"])

            volume_diff_value = sfincs['sfincs_results'].get("Volume_diff_from_F(%)", None)
            if volume_diff_value is not None:
                volume_diffs.append(volume_diff_value.values.flatten()[0])

            mean_damage_diff = fiat['fiat_results']['Damage_diff_from_F(%)'][np.isfinite(fiat['fiat_results']['Damage_diff_from_F(%)'])].mean()
            if mean_damage_diff is not None:
                damage_diffs.append(mean_damage_diff)

            CF_info = sfincs["CF_info"]
            active_drivers = [k.title() if k != "SLR" else "SLR" for k in ["rain", "wind", "SLR"] if CF_info.get(k, 0) != 0]
            drivers.append(" & ".join(active_drivers) if len(active_drivers) > 1 else active_drivers[0])
            driver_sets.append(active_drivers)

    df = pd.DataFrame({
        "Model name": model_names,
        "Flood volume difference (F-CF in %)": volume_diffs,
        "Mean damage difference (%)": damage_diffs,
        "CF driver": drivers,
        "Driver set": driver_sets
    })

    color_map = {
        'SLR': '#0072B2',
        'Wind': '#009E73',
        'Rain': '#E69F00',
    }

    sns.set(style="whitegrid")
    fig, ax1 = plt.subplots(figsize=(12, 6))

    ax1.set_xlabel("CF driver", fontsize=14)
    ax1.set_ylabel("Flood volume difference (F-CF in %)", fontsize=14)
    ax1.tick_params(axis='x', labelsize=12)
    ax1.tick_params(axis='y', labelsize=12)
    ax1.set_xticks(range(len(df["CF driver"].unique())))
    ax1.set_xticklabels(df["CF driver"].unique(), rotation=20)
    

    def plot_split_scatter(ax, x, y, driver_set, damage, color_map, radius_base=0.2, radius_scale=0.1):
        radius = radius_base + radius_scale
        angle_step = 360 / len(driver_set)
        for i, driver in enumerate(driver_set):
            start_angle = i * angle_step
            end_angle = (i + 1) * angle_step
            wedge = Wedge((x, y), radius, start_angle, end_angle,
                          facecolor=color_map[driver], edgecolor='black', linewidth=0.3)
            ax.add_patch(wedge)

    unique_x = list(df["CF driver"].unique())
    x_offsets = {driver: i for i, driver in enumerate(unique_x)}
    swarm_offsets = {key: 0.0 for key in unique_x}
    
    y_spread = 0.5  # tighter vertical spread

    for i, row in df.iterrows():
        x_pos = x_offsets[row["CF driver"]]
        y_pos = row["Flood volume difference (F-CF in %)"] + swarm_offsets[row["CF driver"]]
        plot_split_scatter(ax1, x_pos, y_pos, row["Driver set"], row["Mean damage difference (%)"], color_map)
        swarm_offsets[row["CF driver"]] += y_spread * (-1)**i * 1  # alternate up and down, tighter spread

    ax2 = ax1.twinx()
    ax2.set_ylabel("Damage difference (%)", fontsize=14)
    ax2.tick_params(axis='y', labelsize=12)
    ax2.set_yticks([])

    ax1.set_xlim(-0.5, len(unique_x) - 0.5)
    y_min = min(df["Flood volume difference (F-CF in %)"]) - 1.5
    y_max = max(df["Flood volume difference (F-CF in %)"]) + 1.5
    ax1.set_ylim(y_min, y_max)
    

    plt.title("Flood Volume and Damage Difference by CF Driver (Split Points)", fontsize=16)
    plt.tight_layout()
    plt.show()

plot_flood_and_damage_difference_by_driver_split(models, fiat_models)
# %%
# %%

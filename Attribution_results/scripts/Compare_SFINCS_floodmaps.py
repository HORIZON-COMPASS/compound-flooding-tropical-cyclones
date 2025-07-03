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
from matplotlib.patches import Patch
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.colors import LinearSegmentedColormap
import gc
import xarray as xr
import platform, os

prefix = "p:/" if platform.system() == "Windows" else "/p/"

# Function to load YAML configuration file
def load_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)


# Load the SFINCS models and create a model dictonary
def load_sfincs_models(config):
    """Generates model paths and categories for SFINCS runs based on CF values."""
    run = config['runname_ids']['Idai']
    base_path = join(prefix, "11210471-001-compass", "03_Runs", run["region"], run["tc_name"], "sfincs")
    

    models = []
    factual_model = None  # Initialize factual_model before the loop
    
    for rain, wind, slr in itertools.product(run['CF_value_rain'], run['CF_value_wind'], run['CF_value_SLR']):
        model_name = f"event_tp_{run['precip_forcing']}_CF{rain}_{run['tidemodel']}_CF{slr}_{run['wind_forcing']}_CF{wind}"
        model_path = join(base_path, model_name)
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
    base_path = os.path.join(prefix, "11210471-001-compass", "03_Runs", run['region'], run['tc_name'], "fiat")
    
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

        # Free memory where possible
        del fiat_results, fiat_results_spatial, CF_info_str, CF_info
        
    # Ensure the factual model is the first one in the list
    if factual_model is not None:
        models.insert(0, factual_model)
    
    gc.collect()
    return models


# read global surface water occurance (GSWO) data to mask permanent water for the model region
def gwso_sfincs_region(model):
    if platform.system() == "Windows":
        datacat_path = os.path.abspath("../../Workflows/03_data_catalogs/datacatalog_general.yml")
    else:
        datacat_path = os.path.abspath("../../Workflows/03_data_catalogs/datacatalog_general___linux.yml")
    data_catalog = DataCatalog(data_libs = [datacat_path])
    sfincs_region = model["sfincs_model"].region
    gwso_region = data_catalog.get_rasterdataset("gswo", geom=sfincs_region, buffer=1000)
    return gwso_region


# Compute the maximum water level (hmax) and mask out permanent water
def compute_hmax_masked(models, gwso_region):
    import gc
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

        del da_hmax, da_zsmax, da_dep, gswo_mask  # Clean up to free memory
        gc.collect()
        
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
    Efficiently plots masked hmax for each model, loading data from disk if needed.
    """
    projection = models[0]['sfincs_model'].crs.to_epsg()
    num_models = len(models)
    num_rows = math.ceil(num_models / num_cols)

    fig, axes = plt.subplots(nrows=num_rows, 
                             ncols=num_cols, 
                             figsize=figsize, 
                             constrained_layout=True, 
                             subplot_kw={"projection": ccrs.epsg(projection)})

    fig.suptitle("Masked hmax", fontsize=12, y=1.02)

    for model, ax in zip(models, axes.flatten()):
        # Load hmax_masked from file if only path is stored
        hmax_masked = model['sfincs_results'].get('hmax_masked', None)
        if hmax_masked is None:
            hmax_masked = xr.open_dataarray(model['sfincs_results']['hmax_masked_path'])
        # Plot and immediately close if loaded from file
        im = hmax_masked.plot.pcolormesh(
            ax=ax, cmap="Blues", vmin=0, vmax=5.0, add_colorbar=False
        )
        ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=12, crs=hmax_masked.rio.crs, attribution=False)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(f"{model['cat_short']} ({model['CF_info_str']})", fontsize=11)
        
        hmax_masked.close()
        del hmax_masked
       

    fig.colorbar(im, ax=axes, orientation="vertical", fraction=0.02, pad=0.04).set_label('Masked hmax (m)', rotation=270, labelpad=15)
    fig.savefig(join("../figures/all_masked_hmax.png"), dpi=300, bbox_inches='tight')
    plt.show()
    gc.collect()


# Compute the differences between the Factual and Counterfactual masked hmax variables
def compute_hmax_diff(models):
    factual_hmax = None
    hmax_masked = None
    hmax_diff   = None

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
            
    del factual_hmax, hmax_masked, hmax_diff  # Clean up to free memory
    gc.collect()       

    return models


def plot_hmax_diff(models, num_cols=1, zoom_region_latlon=None, add_basemap = True, save_path="../figures/all_hmax_diff.png"):
    """
    Plots hmax_diff for counterfactual models with optional zoom, using memory-efficient handling.
    """
    counterfactual_models = [m for m in models if "hmax_diff" in m["sfincs_results"]]
    if not counterfactual_models:
        print("No counterfactual models with 'hmax_diff' found.")
        return

    model_crs = counterfactual_models[0]['sfincs_model'].crs
    model_epsg = model_crs.to_epsg()

    # Convert zoom region to projected coords
    zoom_region = None
    if zoom_region_latlon:
        lat_min, lat_max, lon_min, lon_max = zoom_region_latlon
        transformer = Transformer.from_crs("EPSG:4326", model_epsg, always_xy=True)
        xmin, ymin = transformer.transform(lon_min, lat_min)
        xmax, ymax = transformer.transform(lon_max, lat_max)
        zoom_region = (xmin, xmax, ymin, ymax)

    num_models = len(counterfactual_models)
    if num_models > 1:
        num_cols = 2
    num_rows = math.ceil(num_models / num_cols)

    fig, axes = plt.subplots(
        nrows=num_rows, ncols=num_cols,
        figsize=(num_cols * 5, num_rows * 5),
        constrained_layout=True,
        subplot_kw={"projection": ccrs.epsg(model_epsg)}
    )

    # Flatten axes
    axes = [axes] if num_models == 1 else axes.flatten()
    cmap = "RdBu"
    vmin, vmax = -0.3, 0.3
    colorbar_added = False

    for idx, (model, ax) in enumerate(zip(counterfactual_models, axes)):
        hmax_diff = model["sfincs_results"].get("hmax_diff")
        if hmax_diff is None:
            print(f"Missing hmax_diff for model: {model['model_name']}")
            continue

        # Plot diff
        im = hmax_diff.plot.pcolormesh(
            ax=ax, cmap=cmap, vmin=vmin, vmax=vmax,
            add_colorbar=False
        )
        # or: ctx.providers.Esri.WorldImagery
        if add_basemap:
            try:
                ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik, zoom=8, crs=model_crs, attribution=False)
            except Exception as e:
                print(f"Warning: Basemap not added for model {model['model_name']}: {e}")

        # Apply zoom
        if zoom_region:
            ax.set_xlim(zoom_region[0], zoom_region[1])
            ax.set_ylim(zoom_region[2], zoom_region[3])

        # Title & axis labels
        non_zero_CF_info = {k: v for k, v in model["CF_info"].items() if v != 0}
        cf_info_str = ", ".join(f"{k}: {v}" for k, v in non_zero_CF_info.items())
        ax.set_title(f"{model['cat_short']} ({cf_info_str})", fontsize=12)

        utmzone = model.get("utmzone", "UTM")
        ax.set_xlabel(f"x [{utmzone}] (m)", fontsize=10)
        ax.set_ylabel(f"y [{utmzone}] (m)", fontsize=10)
        ax.tick_params(axis="both", labelsize=9)

        

        # Free memory (optional if images are large)
        del hmax_diff
    
    # Add colorbar only once
    if not colorbar_added:
        cbar = fig.colorbar(im, ax=axes, orientation="vertical", fraction=0.02, pad=0.02)
        cbar.set_label("Δ Max Water Level (m)", rotation=270, labelpad=15, fontsize=12)
        cbar.ax.tick_params(labelsize=10)
        colorbar_added = True

    # Hide unused axes
    for idx in range(len(counterfactual_models), len(axes)):
        fig.delaxes(axes[idx])

    # Save and close
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)  # Release memory


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
    
    del hmax_masked
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

    del hmax_masked, flooded_cells, flood_extent, flood_extent_km2  # Clean up to free memory
    gc.collect()
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
        
    del hmax_masked, flooded_cells, flood_volume, flood_volume_km3  # Clean up to free memory
    gc.collect()
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
        
    del factual_flood_volume, factual_flood_extent  # Clean up to free memory
    gc.collect()

    return models


def calculate_damage_differences(fiat_models):
    factual_total_damage = None  # Variable to store the factual total damage

    for model in fiat_models:
        # Store factual total damage for comparison
        if model["category"] == "Factual":
            factual_total_damage = model['fiat_results'].get("total_damage", None).sum()
            factual_total_damage_spatial = model['fiat_results_spatial'].get("total_damage", None).sum()

            if factual_total_damage is None:
                print(f"Error: 'total_damage' not found for factual model: {model['model_name']}")
                continue  # Skip this model if factual total damage is missing

        # Compute damage difference for counterfactual models
        if factual_total_damage is not None and model["category"] != "Factual":
            total_damage = model['fiat_results'].get('total_damage', None).sum()
            if total_damage is None:
                print(f"Error: 'total_damage' not found for counterfactual model: {model['model_name']}")
                continue  # Skip this model if total damage is missing
            
            # Calculate the damage difference from the factual model (in percentage)
            damage_diff = (factual_total_damage - total_damage) / factual_total_damage * 100
            model['Damage_diff_from_F(%)'] = damage_diff
            print(f"damage_diff calculated for {model['model_name']}")
        
        if factual_total_damage_spatial is not None and model["category"] != "Factual":
            total_spatial_damage = model['fiat_results_spatial'].get('total_damage', None).sum()
            if total_spatial_damage is None:
                print(f"Error: 'total_damage' not found for counterfactual model: {model['model_name']}")
                continue  # Skip this model if total damage is missing
            
            # Calculate the damage difference from the factual model (in percentage)
            damage_diff = (factual_total_damage_spatial - total_spatial_damage) / factual_total_damage_spatial * 100
            model['Damage_spatial_diff_from_F(%)'] = damage_diff
            print(f"damage_diff calculated for {model['model_name']}")

    del factual_total_damage, factual_total_damage_spatial  # Clean up to free memory
    gc.collect()

    return fiat_models

################# PLOTTING ###################
# Function to create a categorical plot comparing flood volume differences by CF driver

def plot_abs_flood_difference_by_driver(sfincs_models):
    model_dict = {}

    for sf in sfincs_models:
        vol = sf['sfincs_results'].get("flood_volume_km3", None)
        if vol is None:
            continue
        try:
            vol = float(vol)
        except Exception:
            continue

        CF_info = sf["CF_info"]
        drivers = tuple(sorted(k.upper() for k in ["rain", "wind", "SLR"] if CF_info.get(k, 0) != 0))
        if not drivers:
            drivers = ("FACTUAL",)

        model_dict[drivers] = vol  # Only one model per driver combo assumed

    # Sort keys by combination length then alphabetically
    sorted_keys = sorted(model_dict.keys(), key=lambda k: (len(k), k))

    # Setup plot
    fig, ax = plt.subplots(figsize=(10, 5))  # Slightly smaller figsize
    x = np.arange(len(sorted_keys))
    width = 0.4

    for i, key in enumerate(sorted_keys):
        vol = model_dict[key]
        label = " & ".join(key)
        if label == "RAIN & SLR & WIND":
            label = "ALL"

        color = "#b0b0b0" if label == "FACTUAL" else "#97a6c4"
        ax.bar(x[i], vol, color=color, width=width, edgecolor='black')
        ax.text(x[i], vol + 0.1, f"{vol:.2f}", ha='center', va='bottom', fontsize=11)

        # Remove label text immediately if not needed further
        del label

    ax.set_xticks(x)
    ax.set_xticklabels([" & ".join(k) if k != ("RAIN", "SLR", "WIND") else "ALL" for k in sorted_keys], fontsize=12)
    ax.set_ylabel("Flood volume (km³)", fontsize=14)
    ax.set_title("Factual and Counterfactual Flood Volume", fontsize=14)
    ax.tick_params(axis='y', labelsize=12)
    ax.set_xlim(left=-0.5)
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)

    ax.legend(handles=[
        Patch(facecolor='#b0b0b0', edgecolor='black', label='Factual'),
        Patch(facecolor='#97a6c4', edgecolor='black', label='Counterfactual')
    ], loc='upper left', bbox_to_anchor=(1, 1), fontsize=13)

    plt.tight_layout()
    plt.savefig("../figures/abs_flood_difference_by_driver.png", dpi=300, bbox_inches="tight")
    plt.close(fig)  # More memory-efficient than plt.show() for batch runs

    # Clean-up
    del model_dict, sorted_keys, x, fig, ax
    gc.collect()


def plot_abs_flood_ext_difference_by_driver(sfincs_models):
    # Collect data by driver combination
    model_dict = {}
    for sf in (sfincs_models):
        name = sf['model_name']
        
        # Flood extent
        ext = sf['sfincs_results'].get("flood_extent_km2", None)
        if ext is None:
            continue
        ext = ext.values.flatten()[0]

        # Drivers active
        CF_info = sf["CF_info"]
        drivers = [k.upper() for k in ["rain", "wind", "SLR"] if CF_info.get(k, 0) != 0]
        key = tuple(sorted(drivers)) if drivers else ("FACTUAL",)
        model_dict[key] = {'extent': ext}

    # Build data for plotting
    all_keys = sorted(model_dict.keys(), key=lambda k: (len(k), k))
    data_plot = []

    for key in all_keys:
        ext_total = model_dict[key]['extent']

        label = " & ".join(key)
        if label == "RAIN & SLR & WIND":
            label = "ALL"

        data_plot.append({
            'label': label,
            'extent': ext_total,
        })

    # Sort data by the total magnitude of impact (volume + damage)
    data_plot.sort(key=lambda d: abs(d['extent']))

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(data_plot))

    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)  # This makes grid lines render below plot elements
    width = 0.4

    # Iterate over each model and plot the bars
    for i, d in enumerate(data_plot):
        if d['label'] == "FACTUAL":
            color = "#b0b0b0"  # grey for factual
        else:
            color = "#97a6c4"  # consistent pastel blue for all others
        ax.bar(x[i], d['extent'], color=color, width=width, edgecolor='black')
        ax.text(x[i], d['extent'] + 0.1, f"{d['extent']:.2f}", ha='center', va='bottom', fontsize=13)
        
    # Axis setup
    ax.set_xticks(x)
    ax.set_xticklabels([d['label'] for d in data_plot], fontsize=14)
    ax.set_ylabel("Flood extent (km²)", fontsize=16)
    ax.set_title("Factual and Counterfactual Flood Extent", fontsize=16)
    ax.tick_params(axis='y', labelsize=13)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

    # Set x-axis to start from zero
    ax.set_xlim(left=-0.5)

    legend_elements = [
        Patch(facecolor='#b0b0b0', edgecolor='black', label='Factual'),
        Patch(facecolor='#97a6c4', edgecolor='black', label='Counterfactual'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1),  fontsize=16)

    # Layout
    plt.savefig("../figures/abs_flood_ext_difference_by_driver.png", dpi=300, bbox_inches="tight")
    
    plt.tight_layout()
    plt.show()


def plot_abs_damage_difference_by_driver(fiat_models):
    # Collect data by driver combination
    model_dict = {}
    for fiat in (fiat_models):
       
        # Flood damage
        dam = fiat['fiat_results']['total_damage'].sum()

        # Drivers active
        CF_info = fiat["CF_info"]
        drivers = [k.upper() for k in ["rain", "wind", "SLR"] if CF_info.get(k, 0) != 0]
        key = tuple(sorted(drivers)) if drivers else ("FACTUAL",)
        model_dict[key] = {'damage': dam}

    # Build data for plotting
    all_keys = sorted(model_dict.keys(), key=lambda k: (len(k), k))
    data_plot = []

    for key in all_keys:
        vol_total = model_dict[key]['damage']

        label = " & ".join(key)
        if label == "RAIN & SLR & WIND":
            label = "ALL"

        data_plot.append({
            'label': label,
            'damage': vol_total,
        })

    # Sort data by the total magnitude of impact (damage)
    data_plot.sort(key=lambda d: abs(d['damage']))

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(data_plot))
    
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)  # This makes grid lines render below plot elements
    width = 0.4

    # Iterate over each model and plot the bars
    for i, d in enumerate(data_plot):
        if d['label'] == "FACTUAL":
            color = "#b0b0b0"  # grey for factual
        else:
            color = "#384860"  # consistent pastel blue for all others
        ax.bar(x[i], d['damage'], color=color, width=width, edgecolor='black')

        ax.text(x[i], d['damage'] + 0.1, f"{d['damage']:.1f}", ha='center', va='bottom', fontsize=13)


    # Axis setup
    ax.set_xticks(x)
    ax.set_xticklabels([d['label'] for d in data_plot], fontsize=14)
    ax.set_ylabel("Flood damage ($)", fontsize=16)
    ax.set_title("Factual and Counterfactual Flood Damage", fontsize=16)
    ax.tick_params(axis='y', labelsize=13)
    ax.xaxis.grid(False)
    
    ax.yaxis.get_offset_text().set_fontsize(12)  # Set your desired font size

    # Set x-axis to start from zero
    ax.set_xlim(left=-0.5)

    legend_elements = [
        Patch(facecolor='#b0b0b0', edgecolor='black', label='Factual'),
        Patch(facecolor='#384860', edgecolor='black', label='Counterfactual'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1),  fontsize=16)
    
    plt.savefig("../figures/abs_damage_diff_by_driver.png", dpi=300, bbox_inches="tight")
    
    plt.tight_layout()
    plt.show()


def plot_driver_combination_volume_damage(sfincs_models, fiat_models, filter_keys=None, tolerance=1e-6):    
    model_dict = {}
    for sf, fiat in zip(sfincs_models, fiat_models):      
        # Volume
        vol = sf['sfincs_results'].get("Volume_diff_from_F(%)", None)
        if vol is None:
            continue
        vol = vol.values.flatten()[0]

        # Damage
        dam = fiat['Damage_diff_from_F(%)']

        CF_info = sf["CF_info"]
        drivers = [k.upper() for k in ["rain", "wind", "SLR"] if CF_info.get(k, 0) != 0]
        key = tuple(sorted(drivers)) if drivers else ("FACTUAL",)
        model_dict[key] = {'volume': vol, 'damage': dam}

    # Build and sort data
    all_keys = sorted(model_dict.keys(), key=lambda k: (len(k), k))
    data_plot = []
    for key in all_keys:
        vol_total = model_dict[key]['volume']
        dam_total = model_dict[key]['damage']
        if abs(vol_total) < tolerance and abs(dam_total) < tolerance:
            continue

        label = " & ".join(key)
        if label == "RAIN & SLR & WIND":
            label = "ALL"
        data_plot.append({'label': label, 'volume': vol_total, 'damage': dam_total, 'key': key})
    
    data_plot.sort(key=lambda d: d['volume'])

    # Apply filtering
    if filter_keys is not None:
        filter_keys_normalized = [tuple(sorted(fk.split(" & "))) if isinstance(fk, str) else tuple(sorted(fk)) for fk in filter_keys]
        data_plot = [d for d in data_plot if tuple(sorted(d['key'])) in filter_keys_normalized]

    if not data_plot:
        print("No data available for the selected driver combinations.")
        return

    # Plot
    fig, ax = plt.subplots(figsize=(14, 8), dpi=300)
    x = np.arange(len(data_plot))

    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)  # This makes grid lines render below plot elements
    width = 0.3

    max_y = 0
    for i, d in enumerate(data_plot):
        # Volume bar
        ax.bar(x[i] - width/2, d['volume'], width=width, color="#97a6c4", edgecolor='black')
        # Damage bar
        ax.bar(x[i] + width/2, d['damage'], width=width, color="#384860", edgecolor='black',)
        
        # Annotate the value on top of the bar
        ax.text(x[i] - width/2, d['volume'] + 0.1, f"{d['volume']:.2f}", ha='center', va='bottom', fontsize=13)
        ax.text(x[i] + width/2, d['damage'] + 0.1, f"{d['damage']:.2f}", ha='center', va='bottom', fontsize=13)

        max_y = max(max_y, d['volume'], d['damage'])

    ax.set_xticks(x)
    ax.set_xticklabels([d['label'] for d in data_plot], fontsize=14)
    ax.set_ylabel("Attribution (F-CF) (%)", fontsize=16)
    ax.set_title("Flood Volume and Damage Attribution by CF Driver (Set)", fontsize=16)
    ax.set_ylim(0, max_y + 5)
    ax.set_xlim(-0.5, len(data_plot) - 0.5)
    ax.tick_params(axis='y', labelsize=12)
    ax.xaxis.grid(False)

    legend_elements = [
        Patch(facecolor='#97a6c4', edgecolor='black', label='Flood Volume'),
        Patch(facecolor='#384860', edgecolor='black', label='Damage'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1), fontsize=16)

    plt.savefig("../figures/volume_damage_diff_combined.png", dpi=300, bbox_inches="tight")
    
    plt.tight_layout(pad=4.0)
    plt.show()


def plot_driver_combination_extent_damage(sfincs_models, fiat_models, filter_keys=None, tolerance=1e-6):    
    model_dict = {}
    for sf, fiat in zip(sfincs_models, fiat_models):      
        # Volume
        ext = sf['sfincs_results'].get("Extent_diff_from_F(%)", None)
        if ext is None:
            continue
        ext = ext.values.flatten()[0]

        # Damage
        dam = fiat['Damage_diff_from_F(%)']

        CF_info = sf["CF_info"]
        drivers = [k.upper() for k in ["rain", "wind", "SLR"] if CF_info.get(k, 0) != 0]
        key = tuple(sorted(drivers)) if drivers else ("FACTUAL",)
        model_dict[key] = {'extent': ext, 'damage': dam}

    # Build and sort data
    all_keys = sorted(model_dict.keys(), key=lambda k: (len(k), k))
    data_plot = []
    for key in all_keys:
        ext_total = model_dict[key]['extent']
        dam_total = model_dict[key]['damage']
        if abs(ext_total) < tolerance and abs(dam_total) < tolerance:
            continue

        label = " & ".join(key)
        if label == "RAIN & SLR & WIND":
            label = "ALL"
        data_plot.append({'label': label, 'extent': ext_total, 'damage': dam_total, 'key': key})
    
    data_plot.sort(key=lambda d: d['damage'])

    # Apply filtering
    if filter_keys is not None:
        filter_keys_normalized = [tuple(sorted(fk.split(" & "))) if isinstance(fk, str) else tuple(sorted(fk)) for fk in filter_keys]
        data_plot = [d for d in data_plot if tuple(sorted(d['key'])) in filter_keys_normalized]

    if not data_plot:
        print("No data available for the selected driver combinations.")
        return

    # Plot
    fig, ax = plt.subplots(figsize=(14, 8), dpi=300)
    x = np.arange(len(data_plot))

    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)  # This makes grid lines render below plot elements
    width = 0.3

    max_y = 0
    for i, d in enumerate(data_plot):
        # Extent bar
        ax.bar(x[i] - width/2, d['extent'], width=width, color="#97a6c4", edgecolor='black')
        # Damage bar
        ax.bar(x[i] + width/2, d['damage'], width=width, color="#384860", edgecolor='black',)
        
        # Annotate the value on top of the bar
        ax.text(x[i] - width/2, d['extent'] + 0.1, f"{d['extent']:.2f}", ha='center', va='bottom', fontsize=13)
        ax.text(x[i] + width/2, d['damage'] + 0.1, f"{d['damage']:.2f}", ha='center', va='bottom', fontsize=13)

        max_y = max(max_y, d['extent'], d['damage'])

    ax.set_xticks(x)
    ax.set_xticklabels([d['label'] for d in data_plot], fontsize=14)
    ax.set_ylabel("Attribution (F-CF) (%)", fontsize=16)
    ax.set_title("Flood Extent and Damage Attribution by CF Driver (Set)", fontsize=16)
    ax.set_ylim(0, max_y + 5)
    ax.set_xlim(-0.5, len(data_plot) - 0.5)
    ax.tick_params(axis='y', labelsize=12)
    ax.xaxis.grid(False)

    legend_elements = [
        Patch(facecolor='#97a6c4', edgecolor='black', label='Flood Extent'),
        Patch(facecolor='#384860', edgecolor='black', label='Damage'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1), fontsize=16)

    plt.savefig("../figures/extent_damage_diff_combined.png", dpi=300, bbox_inches="tight")
    
    plt.tight_layout(pad=4.0)
    plt.show()


# Plot functions with driver decomposition for EGU ppt
def plot_driver_decomposition_volume_only(sfincs_models, fiat_models):
    color_map = {
        'RAIN': '#56B4E9',       # turquoise
        'WIND': '#3AA17E',       # ocean green
        'SLR':  '#8E7CC3',        # purple
    }

    # Collect data by driver combination
    model_dict = {}
    for sf, fiat in zip(sfincs_models, fiat_models):
        name = sf['model_name']
        
        # Flood volume
        vol = sf['sfincs_results'].get("Volume_diff_from_F(%)", None)
        if vol is None:
            continue
        vol = vol.values.flatten()[0]

        # Damage (though not used in this plot, we still need to track the max_y)
        dam = fiat['fiat_results']['Damage_diff_from_F(%)']
        dam = dam[np.isfinite(dam)].mean()
        if dam is None:
            continue

        # Drivers active
        CF_info = sf["CF_info"]
        drivers = [k.upper() for k in ["rain", "wind", "SLR"] if CF_info.get(k, 0) != 0]
        key = tuple(sorted(drivers))
        model_dict[key] = {'volume': vol, 'damage': dam}

    # Build data for plotting
    all_keys = sorted(model_dict.keys(), key=lambda k: (len(k), k))
    data_plot = []

    for key in all_keys:
        vol_total = model_dict[key]['volume']
        dam_total = model_dict[key]['damage']

        # Use individual driver contributions to calculate ratios
        vol_parts = {d: model_dict.get((d,), {}).get('volume', 0) for d in key}
        dam_parts = {d: model_dict.get((d,), {}).get('damage', 0) for d in key}

        # Calculate the total of individual contributions
        vol_sum = sum(vol_parts.values())
        dam_sum = sum(dam_parts.values())

        # Calculate the ratios for coloring
        vol_ratios = {d: vol_parts[d] / vol_sum if vol_sum > 0 else 0 for d in key}
        dam_ratios = {d: dam_parts[d] / dam_sum if dam_sum > 0 else 0 for d in key}

        label = " & ".join(key)
        if label == "RAIN & SLR & WIND":
            label = "COMPOUND"

        data_plot.append({
            'label': label,
            'volume': vol_total,
            'damage': dam_total,
            'volume_ratios': vol_ratios,
            'damage_ratios': dam_ratios
        })

    # Sort data by the total magnitude of impact (volume + damage)
    data_plot.sort(key=lambda d: abs(d['volume']) + abs(d['damage']))

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(data_plot))
    width = 0.35

    max_y = 0  # To dynamically track the maximum y-value for proper axis limits

    # Iterate over each model and plot the bars
    for i, d in enumerate(data_plot):
        vol_bottom = 0
        dam_bottom = 0

        # Plot flood volume for the combination (single bar for the total value)
        for part, ratio in d['volume_ratios'].items():
            # Color the bars based on the ratio for each driver in the combination
            ax.bar(x[i] - width/2, d['volume'] * ratio, bottom=vol_bottom, width=width,
                   color=color_map.get(part, '#cccccc'), edgecolor='black')
            vol_bottom += d['volume'] * ratio

        # Track the maximum y-value considering both flood volume and damage
        max_y = max(max_y, vol_bottom, dam_bottom, d['damage'])

        # Add symbols (~) for flood volume
        ax.text(x[i] - width/2, vol_bottom + 0.5, "~", ha='center', va='bottom', fontsize=18)

    # Axis setup
    ax.set_xticks(x)
    ax.set_xticklabels([d['label'] for d in data_plot])
    ax.set_ylabel("Difference (%)")
    ax.set_title("Flood Volume Differences by CF Driver Set")

    # Set x-axis to start from zero
    ax.set_xlim(left=0)

    # Set y-axis to start at zero (removes extra space below bars) and dynamically adjust based on data
    ax.set_ylim(bottom=0, top=max_y + 5)

    # Legend
    legend_elements = [
        Patch(facecolor=color_map['RAIN'], edgecolor='black', label='Rain'),
        Patch(facecolor=color_map['WIND'], edgecolor='black', label='Wind'),
        Patch(facecolor=color_map['SLR'], edgecolor='black', label='SLR'),
        Patch(facecolor='white', edgecolor='black', label='Flood Volume'),
    ]
    ax.legend(handles=legend_elements, title="Legend", loc='upper left', bbox_to_anchor=(1, 1))

    # Layout
    plt.tight_layout()
    plt.show()


def plot_driver_decomposition_extent_only(sfincs_models, fiat_models):
    color_map = {
        'RAIN': '#56B4E9',       # turquoise
        'WIND': '#3AA17E',       # ocean green
        'SLR':  '#8E7CC3',        # purple
    }

    # Collect data by driver combination
    model_dict = {}
    for sf, fiat in zip(sfincs_models, fiat_models):
        name = sf['model_name']
        
        # Flood extent
        ext = sf['sfincs_results'].get("Extent_diff_from_F(%)", None)
        if ext is None:
            continue
        ext = ext.values.flatten()[0]

        # Damage (though not used in this plot, we still need to track the max_y)
        dam = fiat['fiat_results']['Damage_diff_from_F(%)']
        dam = dam[np.isfinite(dam)].mean()
        if dam is None:
            continue

        # Drivers active
        CF_info = sf["CF_info"]
        drivers = [k.upper() for k in ["rain", "wind", "SLR"] if CF_info.get(k, 0) != 0]
        key = tuple(sorted(drivers))
        model_dict[key] = {'extent': ext, 'damage': dam}

    # Build data for plotting
    all_keys = sorted(model_dict.keys(), key=lambda k: (len(k), k))
    data_plot = []

    for key in all_keys:
        ext_total = model_dict[key]['extent']
        dam_total = model_dict[key]['damage']

        # Use individual driver contributions to calculate ratios
        ext_parts = {d: model_dict.get((d,), {}).get('extent', 0) for d in key}
        dam_parts = {d: model_dict.get((d,), {}).get('damage', 0) for d in key}

        # Calculate the total of individual contributions
        ext_sum = sum(ext_parts.values())
        dam_sum = sum(dam_parts.values())

        # Calculate the ratios for coloring
        ext_ratios = {d: ext_parts[d] / ext_sum if ext_sum > 0 else 0 for d in key}
        dam_ratios = {d: dam_parts[d] / dam_sum if dam_sum > 0 else 0 for d in key}

        label = " & ".join(key)
        if label == "RAIN & SLR & WIND":
            label = "COMPOUND"

        data_plot.append({
            'label': label,
            'extent': ext_total,
            'damage': dam_total,
            'extent_ratios': ext_ratios,
            'damage_ratios': dam_ratios
        })

    # Sort data by the total magnitude of impact (extent + damage)
    data_plot.sort(key=lambda d: abs(d['extent']) + abs(d['damage']))

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(data_plot))
    width = 0.35

    max_y = 0  # To dynamically track the maximum y-value for proper axis limits

    # Iterate over each model and plot the bars
    for i, d in enumerate(data_plot):
        ext_bottom = 0
        dam_bottom = 0

        # Plot flood extent for the combination (single bar for the total value)
        for part, ratio in d['extent_ratios'].items():
            # Color the bars based on the ratio for each driver in the combination
            ax.bar(x[i] - width/2, d['extent'] * ratio, bottom=ext_bottom, width=width,
                   color=color_map.get(part, '#cccccc'), edgecolor='black')
            ext_bottom += d['extent'] * ratio

        # Track the maximum y-value considering both flood extent and damage
        max_y = max(max_y, ext_bottom, dam_bottom, d['damage'])

        # Add symbols (~) for flood extent
        ax.text(x[i] - width/2, ext_bottom + 0.5, "~", ha='center', va='bottom', fontsize=18)

    # Axis setup
    ax.set_xticks(x)
    ax.set_xticklabels([d['label'] for d in data_plot])
    ax.set_ylabel("Difference (%)")
    ax.set_title("Flood Extent Differences by CF Driver Set")

    # Set x-axis to start from zero
    ax.set_xlim(left=0)

    # Set y-axis to start at zero (removes extra space below bars) and dynamically adjust based on data
    ax.set_ylim(bottom=0, top=max_y + 5)

    # Legend
    legend_elements = [
        Patch(facecolor=color_map['RAIN'], edgecolor='black', label='Rain'),
        Patch(facecolor=color_map['WIND'], edgecolor='black', label='Wind'),
        Patch(facecolor=color_map['SLR'], edgecolor='black', label='SLR'),
        Patch(facecolor='white', edgecolor='black', label='Flood Extent'),
    ]
    ax.legend(handles=legend_elements, title="Legend", loc='upper left', bbox_to_anchor=(1, 1))

    # Layout
    plt.tight_layout()
    plt.show()


def plot_driver_decomposition_volume_damage(sfincs_models, fiat_models, filter_keys=None, tolerance=1e-6):
    color_map = {
        'RAIN': '#56B4E9',       # turquoise
        'WIND': '#3AA17E',       # ocean green
        'SLR':  '#8E7CC3',        # purple
    }

    # Collect data by driver combination
    model_dict = {}
    for sf, fiat in zip(sfincs_models, fiat_models):
        name = sf['model_name']
        
        # Flood volume
        vol = sf['sfincs_results'].get("Volume_diff_from_F(%)", None)
        if vol is None:
            continue
        vol = vol.values.flatten()[0]  # Flatten the array to get a scalar value

        # Damage
        dam = fiat['Damage_diff_from_F(%)']

        # Print values for debugging
        print(f"Model: {name}, Volume: {vol}, Damage: {dam}")

        # Drivers active
        CF_info = sf["CF_info"]
        drivers = [k.upper() for k in ["rain", "wind", "SLR"] if CF_info.get(k, 0) != 0]
        key = tuple(sorted(drivers))  # Sort the driver combination to ensure consistency
        model_dict[key] = {'volume': vol, 'damage': dam}

    # Build data for plotting
    all_keys = sorted(model_dict.keys(), key=lambda k: (len(k), k))
    data_plot = []

    for key in all_keys:
        vol_total = model_dict[key]['volume']
        dam_total = model_dict[key]['damage']

        # Use individual driver contributions to calculate ratios
        vol_parts = {d: model_dict.get((d,), {}).get('volume', 0) for d in key}
        dam_parts = {d: model_dict.get((d,), {}).get('damage', 0) for d in key}

        # Calculate the total of individual contributions
        vol_sum = sum(vol_parts.values())
        dam_sum = sum(dam_parts.values())

        # Apply tolerance: skip combinations with very small totals
        if abs(vol_sum) < tolerance and abs(dam_sum) < tolerance:
            continue

        # Calculate the ratios for coloring
        vol_ratios = {d: vol_parts[d] / vol_sum if vol_sum > 0 else 0 for d in key}
        dam_ratios = {d: dam_parts[d] / dam_sum if dam_sum > 0 else 0 for d in key}

        # Generate a label
        label = " & ".join(key)
        if label == "RAIN & SLR & WIND":
            label = "COMPOUND"

        data_plot.append({
            'label': label,
            'volume': vol_total,
            'damage': dam_total,
            'volume_ratios': vol_ratios,
            'damage_ratios': dam_ratios,
            'key': key  # Store the key for filtering purposes
        })

    # Sort data by the total magnitude of impact (volume + damage)
    data_plot.sort(key=lambda d: abs(d['volume']) + abs(d['damage']))

    # Optional filter: only plot selected combinations
    if filter_keys is not None:
        # Normalize filter_keys (handle string and tuple types)
        filter_keys_normalized = [tuple(sorted(fk.split(" & "))) if isinstance(fk, str) else tuple(sorted(fk)) for fk in filter_keys]
        
        # Filter the data_plot based on the normalized keys
        filtered_data = []
        for d in data_plot:
            if tuple(sorted(d['key'])) in filter_keys_normalized:
                filtered_data.append(d)
            else:
                print(f"Filtered out {d['label']} with volume: {d['volume']} and damage: {d['damage']}")  # Debugging info

        data_plot = filtered_data

    # Check if any valid data exists after filtering
    if not data_plot:
        print("No data available for the selected driver combinations.")
        return

    fig, ax = plt.subplots(figsize=(14, 8), dpi=300)
    x = np.arange(len(data_plot))
    width = 0.3

    max_y = 0  # To dynamically track the maximum y-value for proper axis limits

    for i, d in enumerate(data_plot):
        vol_bottom = 0
        dam_bottom = 0

        # Plot flood volume for the combination (single bar for the total value)
        for part, ratio in d['volume_ratios'].items():
            # Color the bars based on the ratio for each driver in the combination
            ax.bar(x[i] - width/2, d['volume'] * ratio, bottom=vol_bottom, width=width,
                   color=color_map.get(part, '#cccccc'), edgecolor='black')
            vol_bottom += d['volume'] * ratio

        # Plot damage for the combination (single bar for the total value)
        for part, ratio in d['damage_ratios'].items():
            ax.bar(x[i] + width/2, d['damage'] * ratio, bottom=dam_bottom, width=width,
                   color=color_map.get(part, '#cccccc'), edgecolor='black', hatch='//')
            dam_bottom += d['damage'] * ratio

        # Update the max_y to ensure the y-axis is properly adjusted
        max_y = max(max_y, vol_bottom, dam_bottom)

        # Add symbols (~) for flood volume and ($) for damage
        ax.text(x[i] - width/2, vol_bottom + 0.5, "~", ha='center', va='bottom', fontsize=22)
        ax.text(x[i] + width/2, dam_bottom + 0.5, "$", ha='center', va='bottom', fontsize=22)

    # Axis setup with larger font sizes
    ax.set_xticks(x)
    ax.set_xticklabels([d['label'] for d in data_plot], fontsize=16)
    ax.set_ylabel("Difference (%)", fontsize=22)
    ax.set_title("Flood Volume and Damage Differences by CF Driver Set", fontsize=22)

    # Set y-axis to start at zero (removes extra space below bars) and dynamically adjust based on data
    ax.set_ylim(bottom=0, top=max_y + 5)
    ax.set_xlim(-0.5, len(data_plot) - 0.5)

    # Increase the font size of the y-axis ticks
    ax.tick_params(axis='y', labelsize=16)

    # Legend with larger font size and title
    legend_elements = [
        Patch(facecolor=color_map['RAIN'], edgecolor='black', label='Rain'),
        Patch(facecolor=color_map['WIND'], edgecolor='black', label='Wind'),
        Patch(facecolor=color_map['SLR'], edgecolor='black', label='SLR'),
        Patch(facecolor='white', edgecolor='black', label='Flood Volume'),
        Patch(facecolor='white', edgecolor='black', hatch='//', label='Damage'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1), fontsize=16)


    # Layout with more space
    plt.tight_layout(pad=4.0)
    plt.show()


def plot_driver_decomposition_extent_damage(sfincs_models, fiat_models, filter_keys=None):
    color_map = {
        'RAIN': '#56B4E9',       # turquoise
        'WIND': '#3AA17E',       # ocean green
        'SLR':  '#8E7CC3',        # purple
    }

    # Collect data by driver combination
    model_dict = {}
    for sf, fiat in zip(sfincs_models, fiat_models):
        name = sf['model_name']
        
        # Flood extent
        ext = sf['sfincs_results'].get("Extent_diff_from_F(%)", None)
        if ext is None:
            continue
        ext = ext.values.flatten()[0]

        # Damage
        dam = fiat['Damage_diff_from_F(%)']

        # Drivers active
        CF_info = sf["CF_info"]
        drivers = [k.upper() for k in ["rain", "wind", "SLR"] if CF_info.get(k, 0) != 0]
        key = tuple(drivers)
        model_dict[key] = {'extent': ext, 'damage': dam}

    # Build data for plotting
    all_keys = sorted(model_dict.keys(), key=lambda k: (len(k), k))
    data_plot = []

    for key in all_keys:
        ext_total = model_dict[key]['extent']
        dam_total = model_dict[key]['damage']

        ext_parts = {d: 0 for d in key}
        dam_parts = {d: 0 for d in key}

        for d in key:
            ext_parts[d] = model_dict.get((d,), {}).get('extent', 0)
            dam_parts[d] = model_dict.get((d,), {}).get('damage', 0)

        ext_sum = sum(ext_parts.values())
        dam_sum = sum(dam_parts.values())

        ext_ratios = {d: ext_parts[d] / ext_sum if ext_sum > 0 else 0 for d in key}
        dam_ratios = {d: dam_parts[d] / dam_sum if dam_sum > 0 else 0 for d in key}

        sorted_key = tuple(sorted(key))
        label = " & ".join(sorted_key)
        if set(sorted_key) == {"RAIN", "SLR", "WIND"}:
            label = "COMPOUND"

        data_plot.append({
            'label': label,
            'extent': ext_total,
            'damage': dam_total,
            'extent_ratios': ext_ratios,
            'damage_ratios': dam_ratios,
            'key': sorted_key
        })

    # Filter data_plot based on filter_keys (only for plotting)
    if filter_keys:
        normalized_filter = [tuple(sorted(fk)) for fk in filter_keys]
        data_plot = [d for d in data_plot if d['key'] in normalized_filter]

    print("\nData to plot (after filtering and processing):")
    for entry in data_plot:
        print(entry)

    if not data_plot:
        print("No data to plot after filtering. Please check the filter_keys.")
        return

    data_plot.sort(key=lambda d: abs(d['extent']) + abs(d['damage']))

    fig, ax = plt.subplots(figsize=(14, 8), dpi=300)
    x = np.arange(len(data_plot))
    width = 0.3
    max_y = 0

    for i, d in enumerate(data_plot):
        ext_bottom = 0
        dam_bottom = 0

        for part, ratio in d['extent_ratios'].items():
            ax.bar(x[i] - width/2, d['extent'] * ratio, bottom=ext_bottom, width=width,
                   color=color_map.get(part, '#cccccc'), edgecolor='black')
            ext_bottom += d['extent'] * ratio

        for part, ratio in d['damage_ratios'].items():
            ax.bar(x[i] + width/2, d['damage'] * ratio, bottom=dam_bottom, width=width,
                   color=color_map.get(part, '#cccccc'), edgecolor='black', hatch='//')
            dam_bottom += d['damage'] * ratio

        max_y = max(max_y, ext_bottom, dam_bottom)

        ax.text(x[i] - width/2, ext_bottom + 0.5, "~", ha='center', va='bottom', fontsize=22)
        ax.text(x[i] + width/2, dam_bottom + 0.5, "$", ha='center', va='bottom', fontsize=22)

    ax.set_xticks(x)
    ax.set_xticklabels([d['label'] for d in data_plot], fontsize=16)
    ax.set_ylabel("Difference (%)", fontsize=20)
    ax.set_title("Flood Extent and Damage Differences by CF Driver Set", fontsize=22)

    ax.set_xlim(-0.5, len(data_plot) - 0.5)
    ax.set_ylim(bottom=0, top=max_y + 5)
    ax.tick_params(axis='y', labelsize=16)
    
    legend_elements = [
        Patch(facecolor=color_map['RAIN'], edgecolor='black', label='Rain'),
        Patch(facecolor=color_map['WIND'], edgecolor='black', label='Wind'),
        Patch(facecolor=color_map['SLR'], edgecolor='black', label='SLR'),
        Patch(facecolor='white', edgecolor='black', label='Flood Extent'),
        Patch(facecolor='white', edgecolor='black', hatch='//', label='Damage'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1), fontsize=16)

    plt.tight_layout(pad=4.0)
    plt.show()


def plot_driver_decomposition_volume_damage_interaction(sfincs_models, fiat_models):
    color_map = {
        'RAIN': '#56B4E9',       # turquoise
        'WIND': '#3AA17E',       # ocean green
        'SLR':  '#8E7CC3',       # purple
        'INTERACTION': '#999999' # gray
    }

    # Collect data by driver combination
    model_dict = {}
    for sf, fiat in zip(sfincs_models, fiat_models):
        name = sf['model_name']
        
        # Flood volume
        vol = sf['sfincs_results'].get("Volume_diff_from_F(%)", None)
        if vol is None:
            continue
        vol = vol.values.flatten()[0]

        # Damage
        dam = fiat['Damage_diff_from_F(%)']

        # Drivers active
        CF_info = sf["CF_info"]
        drivers = [k.upper() for k in ["rain", "wind", "SLR"] if CF_info.get(k, 0) != 0]
        key = tuple(sorted(drivers))
        model_dict[key] = {'volume': vol, 'damage': dam, 'CF_info': CF_info}

    # Build data for plotting
    all_keys = sorted(model_dict.keys(), key=lambda k: (len(k), k))
    data_plot = []

    for key in all_keys:
        vol_total = model_dict[key]['volume']
        dam_total = model_dict[key]['damage']
        CF_info = model_dict[key]['CF_info']

        # Individual driver
        if len(key) == 1:
            data_plot.append({
                'label': key[0],
                'volume_parts': {key[0]: vol_total},
                'damage_parts': {key[0]: dam_total},
                'CF_info': CF_info
            })
        else:
            # Sum of individual contributions
            vol_parts = {d: model_dict.get((d,), {}).get('volume', 0) for d in key}
            dam_parts = {d: model_dict.get((d,), {}).get('damage', 0) for d in key}

            vol_sum = sum(vol_parts.values())
            dam_sum = sum(dam_parts.values())

            # Add interaction effect
            vol_parts['INTERACTION'] = vol_total - vol_sum
            dam_parts['INTERACTION'] = dam_total - dam_sum

            label = " & ".join(key)
            data_plot.append({
                'label': label,
                'volume_parts': vol_parts,
                'damage_parts': dam_parts,
                'CF_info': CF_info
            })

    # Sort data by the combined magnitude of volume and damage (low to high)
    data_plot.sort(key=lambda d: abs(sum(d['volume_parts'].values())) + abs(sum(d['damage_parts'].values())))

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(data_plot))
    width = 0.35

    max_y = 0  # To dynamically track the maximum y-value for proper axis limits

    for i, d in enumerate(data_plot):
        vol_bottom = 0
        dam_bottom = 0

        # Flood volume bar (left)
        for part, val in d['volume_parts'].items():
            # If it's the interaction part, set transparency
            alpha = 0.2 if part == 'INTERACTION' else 1.0
            ax.bar(x[i] - width/2, val, bottom=vol_bottom, width=width,
                   color=color_map.get(part, '#cccccc'), edgecolor='black', alpha=alpha)
            vol_bottom += val

        # Damage bar (right)
        for part, val in d['damage_parts'].items():
            # If it's the interaction part, adjust hatch pattern and transparency
            alpha = 0.7 if part == 'INTERACTION' else 1.0  # Less transparency for negative interactions
            hatch_pattern = '\\\\' if part == 'INTERACTION' and val < 0 else '//'
            
            bar = ax.bar(x[i] + width/2, val, bottom=dam_bottom, width=width,
                         color=color_map.get(part, '#cccccc'), edgecolor='black', hatch=hatch_pattern, alpha=alpha)
            dam_bottom += val

            # For interaction, if negative, make only the bottom border bold
            if part == 'INTERACTION' and val < 0:
                # Get the current bar's rectangle and modify the bottom edge
                rect = bar[0]
                # The bottom edge is at y = dam_bottom (start position) and extends to the bottom
                rect.set_edgecolor('black')  # Set the edge color to black
                rect.set_linewidth(1)  # Set normal line width for the other edges
                
                # Draw the bottom line manually with thicker line
                ax.plot([rect.get_x(), rect.get_x() + rect.get_width()],
                        [dam_bottom, dam_bottom], color='black', linewidth=3)

        # Update the max_y to ensure the y-axis is properly adjusted
        max_y = max(max_y, vol_bottom, dam_bottom)

        # Add actual value text for flood volume above bars
        ax.text(x[i] - width/2, vol_bottom + 0.5, f"{vol_total:.2f}%", ha='center', va='bottom', fontsize=12)

        # Add actual value text for damage above bars
        ax.text(x[i] + width/2, dam_bottom + 0.5, f"${dam_total:.2f}", ha='center', va='bottom', fontsize=12)

        # Add wave symbol (~) above the flood volume bars
        ax.text(x[i] - width/2, vol_bottom + 1.5, "~", ha='center', va='bottom', fontsize=18)

        # Add dollar sign ($) above the damage bars
        ax.text(x[i] + width/2, dam_bottom + 1.5, "$", ha='center', va='bottom', fontsize=18)

    # Axis setup
    ax.set_xticks(x)
    ax.set_xticklabels([d['label'] for d in data_plot], rotation=20)
    ax.set_ylabel("Difference (%)")
    ax.set_title("Flood Volume and Damage Differences by CF Driver Set")

    # Set x-axis to start from zero
    ax.set_xlim(left=0)

    # Set y-axis to start at zero (removes extra space below bars) and dynamically adjust based on data
    ax.set_ylim(bottom=0, top=max_y + 5)

    # Legend
    legend_elements = [
        Patch(facecolor=color_map['RAIN'], edgecolor='black', label='Rain'),
        Patch(facecolor=color_map['WIND'], edgecolor='black', label='Wind'),
        Patch(facecolor=color_map['SLR'], edgecolor='black', label='SLR'),
        Patch(facecolor=color_map['INTERACTION'], edgecolor='black', label='Interaction'),
        Patch(facecolor='white', edgecolor='black', label='Flood Volume'),
        Patch(facecolor='white', edgecolor='black', hatch='//', label='Damage'),
    ]
    ax.legend(handles=legend_elements, title="Legend", loc='upper left', bbox_to_anchor=(1, 1))

    # Layout
    plt.tight_layout()
    plt.show()


def plot_masked_hmax_factual_only(models, num_cols=1, figsize=(6, 6), dpi=300):
    """
    Plots only the factual model's masked hmax flood map with coordinate labels.
    """
    import matplotlib.ticker as mticker
    # Filter for factual model only
    factual_models = [m for m in models if all(v == 0 for v in m['CF_info'].values())]

    if not factual_models:
        print("No factual model found.")
        return

    projection = factual_models[0]['sfincs_model'].crs.to_epsg()
    num_models = len(factual_models)
    num_rows = math.ceil(num_models / num_cols)

    fig, axes = plt.subplots(nrows=num_rows,
                             ncols=num_cols,
                             figsize=figsize,
                             dpi=dpi,
                             constrained_layout=True,
                             subplot_kw={"projection": ccrs.epsg(projection)})

    if num_models == 1:
        axes = [axes]  # Ensure iterable

    fig.suptitle("Factual Max Flood Depth", fontsize=11)

    for model, ax in zip(factual_models, axes):
        # Plot hmax masked
        im = model['sfincs_results']['hmax_masked'].plot.pcolormesh(
            ax=ax, cmap="Blues", vmin=0, vmax=5.0, add_colorbar=False
        )

        # Add basemap
        ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery,
                        zoom=10,
                        crs=model['sfincs_results']['hmax_masked'].rio.crs,
                        attribution=False)

        # Add gridlines and format tick labels
        gl = ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {'size': 8}
        gl.ylabel_style = {'size': 8}
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, .2))
        gl.ylocator = mticker.FixedLocator(np.arange(-90, 90, .2))


        # Titles
        ax.set_title(f"", fontsize=16)

        # Add shared colorbar with larger size and labels
        cbar = fig.colorbar(im, ax=axes, orientation="vertical", 
                            fraction=0.035, pad=0.05)
        cbar.set_label('Flood depth (m)', rotation=270, labelpad=10, fontsize=9)
        cbar.ax.tick_params(labelsize=8)
    
    fig.savefig("../figures/factual_hmax_masked.png", bbox_inches='tight', dpi=dpi)
    plt.show()


def plot_hmax_diff_rain_only(models, num_cols=1, zoom_region_latlon=None, savefig_path=None):
    """
    Plots the hmax_diff for RAIN-only counterfactual models with enhanced visuals,
    including latitude/longitude gridlines for presentations.
    """
    import cartopy.crs as ccrs
    import matplotlib.ticker as mticker
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

    rain_only_models = [
        m for m in models 
        if "hmax_diff" in m["sfincs_results"] and m["CF_info"].get("rain", 0) != 0 
        and all(m["CF_info"].get(k, 0) == 0 for k in ["wind", "SLR"])
    ]

    if not rain_only_models:
        print("No RAIN-only models with 'hmax_diff' found.")
        return

    model_crs = rain_only_models[0]['sfincs_model'].crs
    model_epsg = model_crs.to_epsg()

    zoom_region = None
    if zoom_region_latlon:
        transformer = Transformer.from_crs("EPSG:4326", model_epsg, always_xy=True)
        lat_min, lat_max, lon_min, lon_max = zoom_region_latlon
        xmin, ymin = transformer.transform(lon_min, lat_min)
        xmax, ymax = transformer.transform(lon_max, lat_max)
        zoom_region = (xmin, xmax, ymin, ymax)

    num_models = len(rain_only_models)
    num_cols = 2 if num_models > 1 else 1
    num_rows = math.ceil(num_models / num_cols)

    fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols,
                             figsize=(num_cols * 7, num_rows * 6),
                             constrained_layout=True,
                             subplot_kw={"projection": ccrs.epsg(model_epsg)},
                             dpi=300)

    axes = np.array(axes).flatten() if num_models > 1 else [axes]
    
    cmap = LinearSegmentedColormap.from_list("white_red", ["red", "white"])
    vmin, vmax = 0, -0.4  # Adjust range if needed

    for i, model in enumerate(rain_only_models):
        ax = axes[i]
        hmax_diff = model["sfincs_results"]["hmax_diff"]

        # Plot
        im = hmax_diff.plot.pcolormesh(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False)
        ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=10,
                        crs=model_crs, attribution=False)

        # Gridlines with lat/lon labels
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.6, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"size": 10}
        gl.ylabel_style = {"size": 10}
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlocator = mticker.MaxNLocator(nbins=5)
        gl.ylocator = mticker.MaxNLocator(nbins=5)

        # Apply zoom
        if zoom_region:
            ax.set_xlim(zoom_region[0], zoom_region[1])
            ax.set_ylim(zoom_region[2], zoom_region[3])

        # cf_str = ", ".join(f"{k}: {v}" for k, v in model["CF_info"].items() if v != 0)
        ax.set_title(f"Counterfactual Rain", fontsize=14)

    for j in range(num_models, len(axes)):
        fig.delaxes(axes[j])

    # Colorbar
    cbar = fig.colorbar(im, ax=axes, orientation="vertical", fraction=0.035, pad=0.04)
    cbar.set_label('Difference in Max Water Level (m)', fontsize=14, rotation=270, labelpad=20)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_ticks(np.arange(-0.4, 0.1, 0.1))  # From -0.4 to 0 with steps of 0.1

    # if savefig_path:
    #     plt.savefig(savefig_path, bbox_inches='tight', dpi=300)

    plt.show()


def plot_hmax_diff_slr_rain(models, zoom_region_latlon=None):
    """
    Plots the hmax_diff for RAIN-only and SLR-only models side by side with shared colorbar.
    """


    # Filter first RAIN-only and SLR-only models
    rain_model = next((m for m in models if m["CF_info"].get("rain", 0) != 0 and m["CF_info"].get("SLR", 0) == 0 and "hmax_diff" in m["sfincs_results"]), None)
    slr_model = next((m for m in models if m["CF_info"].get("SLR", 0) != 0 and m["CF_info"].get("rain", 0) == 0 and "hmax_diff" in m["sfincs_results"]), None)

    if not rain_model or not slr_model:
        print("Either RAIN-only or SLR-only model with 'hmax_diff' not found.")
        return

    model_epsg = rain_model['sfincs_model'].crs.to_epsg()


    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300, constrained_layout=True,
                             subplot_kw={"projection": ccrs.epsg(model_epsg)})

    cmap = LinearSegmentedColormap.from_list("white_red", ["red", "white"])
    vmin, vmax = -0.4, 0
    plots = [(rain_model, axes[0], "Rain"), (slr_model, axes[1], "SLR")]

    for model, ax, title in plots:
        hmax_diff = model["sfincs_results"]["hmax_diff"]

        im = hmax_diff.plot.pcolormesh(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False)
        ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=9,
                        crs=model['sfincs_model'].crs, attribution=False)

        # Add gridlines with lat/lon
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = gl.right_labels = False
        gl.xlabel_style = gl.ylabel_style = {'size': 11}

        ax.set_title(title, fontsize=16)  # You can increase to 16 or more

    # Shared colorbar
    cbar = fig.colorbar(im, ax=axes, orientation="vertical", fraction=0.02, pad=0.02)
    cbar.set_label('Difference in Maximum Water Level (m)', rotation=270, labelpad=20, fontsize=13)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_ticks(np.arange(-0.4, 0.1, 0.1))

    
    plt.show()


def plot_hmax_diff_slr_wind_rain(models, zoom_region_latlon=None):
    # Filter first sginle driver only models
    rain_model = next((m for m in models if m["CF_info"].get("rain", 0) != 0 and m["CF_info"].get("SLR", 0) == 0 and m["CF_info"].get("wind", 0) == 0 and "hmax_diff" in m["sfincs_results"]), None)
    slr_model  = next((m for m in models if m["CF_info"].get("SLR", 0) != 0 and m["CF_info"].get("rain", 0) == 0 and m["CF_info"].get("wind", 0) == 0 and "hmax_diff" in m["sfincs_results"]), None)
    wind_model = next((m for m in models if m["CF_info"].get("wind", 0) != 0 and m["CF_info"].get("rain", 0) == 0 and m["CF_info"].get("SLR", 0) == 0 and "hmax_diff" in m["sfincs_results"]), None)

    if not rain_model or not slr_model or not wind_model:
        print("Either RAIN-only or SLR-only or WIND-only model with 'hmax_diff' not found.")
        return

    model_epsg = rain_model['sfincs_model'].crs.to_epsg()


    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(14, 6), dpi=300, constrained_layout=True,
                             subplot_kw={"projection": ccrs.epsg(model_epsg)})

    cmap = LinearSegmentedColormap.from_list("white_red", ["red", "white"])
    vmin, vmax = -0.4, 0
    plots = [(rain_model, axes[0], "Rain"), (slr_model, axes[1], "SLR"), (wind_model, axes[1], "Wind")]

    for model, ax, title in plots:
        hmax_diff = model["sfincs_results"]["hmax_diff"]

        im = hmax_diff.plot.pcolormesh(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False)
        ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=9,
                        crs=model['sfincs_model'].crs, attribution=False)

        # Add gridlines with lat/lon
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = gl.right_labels = False
        gl.xlabel_style = gl.ylabel_style = {'size': 11}

        ax.set_title(title, fontsize=16)  # You can increase to 16 or more

    # Shared colorbar
    cbar = fig.colorbar(im, ax=axes, orientation="vertical", fraction=0.02, pad=0.02)
    cbar.set_label('Difference in Maximum Water Level (m)', rotation=270, labelpad=20, fontsize=13)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_ticks(np.arange(-0.4, 0.1, 0.1))

    fig.savefig("../figures/hmax_diff_slr_widn_rain.png", bbox_inches='tight', dpi=300)
    plt.show()


def plot_driver_decomposition_extent(sfincs_models, filter_keys=None):
    color_map = {
        'RAIN': '#56B4E9',       # turquoise
        'WIND': '#3AA17E',       # ocean green
        'SLR':  '#8E7CC3',        # purple
    }

    # Collect data by driver combination
    model_dict = {}
    for sf in (sfincs_models):
        name = sf['model_name']
        
        # Flood extent
        ext = sf['sfincs_results'].get("Extent_diff_from_F(%)", None)
        if ext is None:
            continue
        ext = ext.values.flatten()[0]

        # Drivers active
        CF_info = sf["CF_info"]
        drivers = [k.upper() for k in ["rain", "wind", "SLR"] if CF_info.get(k, 0) != 0]
        key = tuple(drivers)
        model_dict[key] = {'extent': ext}

    # Build data for plotting
    all_keys = sorted(model_dict.keys(), key=lambda k: (len(k), k))
    data_plot = []

    for key in all_keys:
        ext_total = model_dict[key]['extent']

        ext_parts = {d: 0 for d in key}

        for d in key:
            ext_parts[d] = model_dict.get((d,), {}).get('extent', 0)

        ext_sum = sum(ext_parts.values())

        ext_ratios = {d: ext_parts[d] / ext_sum if ext_sum > 0 else 0 for d in key}

        sorted_key = tuple(sorted(key))
        label = " & ".join(sorted_key)
        if set(sorted_key) == {"RAIN", "SLR", "WIND"}:
            label = "COMPOUND"

        data_plot.append({
            'label': label,
            'extent': ext_total,
            'extent_ratios': ext_ratios,
            'key': sorted_key
        })

    # Filter data_plot based on filter_keys (only for plotting)
    if filter_keys:
        normalized_filter = [tuple(sorted(fk)) for fk in filter_keys]
        data_plot = [d for d in data_plot if d['key'] in normalized_filter]

    print("\nData to plot (after filtering and processing):")
    for entry in data_plot:
        print(entry)

    if not data_plot:
        print("No data to plot after filtering. Please check the filter_keys.")
        return

    data_plot.sort(key=lambda d: abs(d['extent']))

    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)
    x = np.arange(len(data_plot))
    width = 0.5
    max_y = 0

    for i, d in enumerate(data_plot):
        ext_bottom = 0

        for part, ratio in d['extent_ratios'].items():
            ax.bar(x[i], d['extent'] * ratio, bottom=ext_bottom, width=width,
                color=color_map.get(part, '#cccccc'), edgecolor='black')
            ext_bottom += d['extent'] * ratio

        max_y = max(max_y, ext_bottom)

    ax.set_xticks(x)
    ax.set_xticklabels([d['label'] for d in data_plot], fontsize=18)
    ax.set_ylabel("Flood difference (%)", fontsize=20)
    ax.set_title("Flood Extent Differences by Counterfactual Driver", fontsize=22)

    ax.set_xlim(-0.5, len(data_plot) - 0.5)
    ax.set_ylim(bottom=0, top=max_y + 5)
    ax.tick_params(axis='y', labelsize=16)
    
    legend_elements = [
        Patch(facecolor=color_map['RAIN'], edgecolor='black', label='Rain'),
        Patch(facecolor=color_map['WIND'], edgecolor='black', label='Wind'),
        Patch(facecolor=color_map['SLR'], edgecolor='black', label='SLR'),
        ]
    # ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1), fontsize=18)

    plt.tight_layout(pad=4.0)
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

# %%
fiat_models = load_fiat_models(cfg)
#%%
fiat_models = calculate_damage_differences(fiat_models)

#%%
# SOME PLOTTING of spatial flood maps
# plot_hmax_diff(models)

# plot_masked_hmax_factual_only(models)
# plot_masked_hmax(models)

# # PLOTTING paper outline
# plot_abs_flood_difference_by_driver(models)
# plot_abs_flood_ext_difference_by_driver(models)
# plot_abs_damage_difference_by_driver(fiat_models)

# plot_driver_combination_volume_damage(models, fiat_models)
# plot_driver_combination_extent_damage(models, fiat_models)

#%%
plot_hmax_diff_slr_wind_rain(models)

# %%
################################################
########## PPT SLIDES FIGS EGU 2025 ############
################################################
# plotting all driver(s) combinations for hazard and impact
# plot_driver_decomposition_volume_damage(models, fiat_models)
# plot_driver_decomposition_volume_damage(models, fiat_models, filter_keys=[("WIND","SLR"), ("RAIN",), ("RAIN", "WIND", "SLR")])

# plot_driver_decomposition_extent_damage(models, fiat_models)

# # Plotting potential interaction of drivers, compared to their individual sum
# plot_driver_decomposition_volume_damage_interaction(models, fiat_models)

# # plotting the change in flood extent and damage for specific driver(s) combinations
# plot_driver_decomposition_extent_damage(models, fiat_models, filter_keys=[("WIND","SLR"), ("RAIN", "WIND", "SLR")])
# plot_driver_decomposition_extent_damage(models, fiat_models, filter_keys=[("WIND","SLR"), ("RAIN",), ("RAIN", "WIND", "SLR")])

# #%%
# # Spatial plot of factual simulated compound flooding
# plot_masked_hmax_factual_only(models)

# # Plot the difference in flooding due to changes in rain only - spatial map
# plot_hmax_diff_rain_only(models)

# # Plot the difference in flooding for rain and SLR in one plot spatially
# plot_hmax_diff_slr_rain(models)

# # Plot the change in flood extent for specific driver (combinations)
# plot_driver_decomposition_extent(models,  filter_keys=[("WIND","SLR"), ("RAIN",), ("RAIN", "WIND", "SLR")])


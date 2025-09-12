#%% use pixi environment compass-wflow
# Load the necessary packages
import os
from os.path import join
import yaml
import itertools
import gc
import xarray as xr
import pandas as pd
import numpy as np
import platform
import gc

from hydromt_sfincs import SfincsModel, utils
from hydromt import DataCatalog

import matplotlib.pyplot as plt
from pyproj import Transformer
import geopandas as gpd
import cartopy.crs as ccrs
from matplotlib.patches import Patch
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.colors import LinearSegmentedColormap
from shapely.geometry import box
from rasterio.features import shapes
from shapely.geometry import shape
from matplotlib.cm import ScalarMappable
import matplotlib.patheffects as path_effects
from matplotlib.colors import Normalize
from matplotlib.colors import PowerNorm
from matplotlib.colors import TwoSlopeNorm
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

prefix = "p:/" if platform.system() == "Windows" else "/p/"

# General functions
def get_driver_group_model(group, models):
    return next(
        (m for m in models
         if all(m["CF_info"].get(d, 0) != 0 for d in group)
         and all(m["CF_info"].get(d, 0) == 0 for d in {"rain", "SLR", "wind"} - set(group))
         and "hmax_diff" in m["sfincs_results"]), None)

def get_single_driver_model(driver, models):
    return next(
        (m for m in models
         if m["CF_info"].get(driver, 0) != 0
         and all(m["CF_info"].get(d, 0) == 0 for d in {"rain", "SLR", "wind"} - {driver})
         and "hmax_diff" in m["sfincs_results"]), None)

def lat_formatter(x, pos):
    direction = 'N' if x >= 0 else 'S'
    return f"{abs(x):.1f}°{direction}"

def lon_formatter(x, pos):
    direction = 'E' if x >= 0 else 'W'
    return f"{abs(x):.1f}°{direction}"

def format_driver_label(drivers) -> str:
    """Format driver labels consistently."""
    if isinstance(drivers, (list, tuple)):
        if set(drivers) == {"RAIN", "SLR", "WIND"}:
            return "All"
        return " & ".join(d.upper() if d.upper() == "SLR" else d.capitalize() for d in drivers)
    return drivers.upper() if drivers.upper() == "SLR" else drivers.capitalize()

def format_pct(val, show_zero=False, signed=False):
    if val == 0 or val is None:
        return "0 %" if show_zero else ""
    abs_val = abs(val)
    if abs_val < 1:
        return "<1 %"
    if signed:
        return f"{val:+.0f} %"
    else:
        return f"{int(round(val))} %"


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
        his_path   = os.path.join(model_path,"sfincs_his.nc")
        ds_his     = xr.open_dataset(his_path, engine="netcdf4")

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
            "sfincs_his": ds_his,
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
        fiat_results =  gpd.read_file(join(f"{model_path}", "output/spatial.fgb"))

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
        del fiat_results, CF_info_str, CF_info
        
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
def compute_hmax_masked(models, gwso_region, model_region_gdf):
    # we set a threshold to mask minimum flood depth
    hmin = 0.05

    for model in models:
        # select the highest-resolution elevation dataset
        print(f"Processing model: {model['model_name']}")
        depfile = join(model["model_path"], "subgrid", "dep_subgrid.tif")
        da_dep = model["sfincs_model"].data_catalog.get_rasterdataset(depfile)

        # compute the maximum over all time steps
        # First timestep leads to incorrect diiference values for permanent water cells that are incorrectly unmasked; requires filtering.
        da_zsmax = model["sfincs_results"]["zsmax"].isel(timemax=slice(1, None)).max(dim="timemax")
     
        # downscale the floodmap
        da_hmax = utils.downscale_floodmap(
            zsmax=da_zsmax,
            dep=da_dep,
            hmin=hmin,
            reproj_method = "bilinear"
            )
    
        # GSWO dataset to mask permanent water by first geprojecting it to the subgrid of hmax
        gswo_mask = gwso_region.raster.reproject_like(da_hmax, method="max")
        # permanent water where water occurence > 5%
        da_hmax_masked = da_hmax.where(gswo_mask <= 5)

        # Add the name attribute for identification
        model["sfincs_results"]['hmax'] = da_hmax
        model["sfincs_results"]['hmax_masked'] = da_hmax_masked

        valid_mask = (gswo_mask <= 5).astype("uint8").squeeze()

        # Extract shapes
        shapes_gen = shapes(valid_mask.values, transform=valid_mask.rio.transform())
        valid_polygons = [shape(geom) for geom, val in shapes_gen if val == 1]
        gdf_valid = gpd.GeoDataFrame(geometry=valid_polygons, crs=gswo_mask.rio.crs)
        gdf_valid = gdf_valid.to_crs(model_region_gdf.crs)


        del da_hmax, da_zsmax, da_dep, gswo_mask  # Clean up to free memory
        gc.collect()

    return models, gdf_valid


# Compute the differences between the Factual and Counterfactual masked hmax variables
def compute_hmax_diff(models):
    factual_hmax = None
    hmax_masked = None
    hmax_diff   = None

    # Check if "Factual" category exists
    for model in models:
        if model["category"] == "Factual":
            factual_hmax = model['sfincs_results'].get("hmax_masked", None)
            mask_F_valid = ~np.isnan(factual_hmax)
            if factual_hmax is None:
                print(f"Error: 'hmax_masked' not found for factual model: {model['model_name']}")
                return models  # Exit early if 'hmax_masked' is missing
            
    # Compute difference for counterfactual models
    for model in models:
        if model["category"] != "Factual":
            hmax_masked = model['sfincs_results'].get("hmax_masked", None)
            mask_CF_valid = ~np.isnan(hmax_masked)

            # For F: where F is NaN but CF has a value, set F to 0
            hmax_F = factual_hmax.where(mask_F_valid | ~mask_CF_valid, 0)
            # For CF: where CF is NaN but F has a value, set CF to 0
            hmax_CF = hmax_masked.where(mask_CF_valid | ~mask_F_valid, 0)
            if hmax_masked is not None:
                hmax_diff = hmax_F - hmax_CF 
                model['sfincs_results']["hmax_diff"] = hmax_diff
                print(f"hmax_diff calculated for {model['model_name']}")
            else:
                print(f"Warning: 'hmax_masked' not found for counterfactual model: {model['model_name']}")
            
    del factual_hmax, hmax_masked, hmax_diff  # Clean up to free memory
    gc.collect()       

    return models


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
    
    area = dx*dy
    print(f"Cell area: {area}")

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
        flooded_cells = hmax_masked > 0 # Necessary for boolean mask!!
        # Compute the total flooded area (in square meters)
        flood_extent = (flooded_cells * calculate_cell_area(model)).sum(dim=['x', 'y']).compute()
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

        # Compute the flooded volume (sum of depth * area for each flooded cell)
        flood_volume = (hmax_masked.where(hmax_masked > 0, 0) * calculate_cell_area(model)).sum().compute()
        # Convert to cubic kilometers
        flood_volume_km3 = flood_volume / 1e9
        # Store in the model dictionary
        model['sfincs_results']['flood_volume_km3'] = flood_volume_km3
        model['sfincs_results']['flood_volume_m3'] = flood_volume
        print(f"for model {model['model_name']}, the flood volume is {flood_volume_km3}")
        
    del hmax_masked, flood_volume, flood_volume_km3  # Clean up to free memory
    gc.collect()
    return models


def calculate_flood_differences(models):
    factual_flood_volume = None
    factual_flood_extent = None
    
    for model in models:
        if model["category"] == "Factual":
            factual_flood_volume = model['sfincs_results']['flood_volume_km3']
            factual_flood_extent = model['sfincs_results']['flood_extent_km2']

        # Compute flood volume difference for counterfactual models
        if model["category"] != "Factual":
            # Counterfactual
            hmax_diff = model['sfincs_results'].get('hmax_diff', None)
            if hmax_diff is None:
                print(f"hmax_diff missing for {model['model_name']}")
                continue

            # Flood volume 
            flood_volume = model['sfincs_results'].get('flood_volume_km3', None)
            flood_volume_diff_perct = (factual_flood_volume - flood_volume) / factual_flood_volume * 100
            model['sfincs_results']['Volume_diff_from_F(%)'] = flood_volume_diff_perct
            print(f"flood_volume_diff calculated for {model['model_name']}")

            # Flooded area 
            flood_extent = model['sfincs_results'].get('flood_extent_km2', None)
            flood_extent_diff = (factual_flood_extent - flood_extent) / factual_flood_extent * 100
            model['sfincs_results']['Extent_diff_from_F(%)'] = flood_extent_diff
            print(f"flood_extent_diff calculated for {model['model_name']}")

            positive_mask = (hmax_diff > 0.01)
            affected_area = (positive_mask * calculate_cell_area(model)).sum().compute() / 1e6  # km²
            print(f"The MOST affected area by climate-change-induced flooding: {round(float(affected_area), -1)} km²")

    del factual_flood_volume, factual_flood_extent  # Clean up to free memory
    gc.collect()

    return models


# # Function to calculate the flood volume and extent differences between factual and counterfactual datasets
# def calculate_flood_differences(models):
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
        
        # Calculate the area affected by the counterfactual flooding, 
        if model["category"] != "Factual":
            hmax_diff = model['sfincs_results']["hmax_diff"]
            dx = abs(hmax_diff.x[1] - hmax_diff.x[0])
            dy = abs(hmax_diff.y[1] - hmax_diff.y[0])
            cell_area = dx * dy  # in map units (e.g., m² if projected)
            total_flooded_area = (cell_area * (hmax_diff).sum().values) / 1e6
            total_flooded_area = round(float(total_flooded_area), -1)

            
            print(f"The affect area by climate-change-induced flooding: {total_flooded_area} km2")

    del factual_flood_volume, factual_flood_extent  # Clean up to free memory
    gc.collect()

    return models


def calculate_damage_differences(fiat_models):
    factual_total_damage = None  # Variable to store the factual total damage

    for model in fiat_models:
        # Store factual total damage for comparison
        if model["category"] == "Factual":
            factual_total_damage = model['fiat_results'].get("total_damage", None).sum()

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

    del factual_total_damage  # Clean up to free memory
    gc.collect()

    return fiat_models

################# PLOTTING ###################
# For paper
def plot_hmax_diff_rain_slrwind_all(models, model_region_gdf, background):
    # Get models
    driver_groups = {
    "Rain": ["rain"],
    "SLR & Wind": ["SLR", "wind"],
    "All": ["rain", "SLR", "wind"]
    }
    subplot_labels = ['(a)', '(b)', '(c)']

    # === Create figure ===
    fig, axes = plt.subplots(1, 3, figsize=(10, 5), dpi=300, constrained_layout=True,
                            subplot_kw={"projection": ccrs.PlateCarree()}, sharey=True)

    # Plot settings
    cmap = LinearSegmentedColormap.from_list("white_red", ["white", "red"])
    vmin, vmax = 0, 0.6
    utm_crs = ccrs.UTM(zone=36, southern_hemisphere=True)
    model_region_gdf = model_region_gdf.to_crs("EPSG:4326")
    background = background.to_crs("EPSG:4326")

    # === Loop through plots ===
    for i, (ax, (title, group)) in enumerate(zip(axes, driver_groups.items())):
        model = get_driver_group_model(group, models)
        if model is None:
            print(f"Missing model for group {group}")
            continue

        # Plotting hmax_diff
        hmax_diff = model["sfincs_results"]["hmax_diff"]
        im = hmax_diff.plot.pcolormesh(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False, transform=utm_crs, rasterized=True)

        # Plot background
        mask_box = box(34.8, -20.3, 35.3, -19.9)  # minx, miny, maxx, maxy
        background_outside_box = background[~background.intersects(mask_box)] # removing errorneous lines outside model region
        background_outside_box.plot(ax=ax, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)
        
        # Plot model region
        model_region_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())

        # Set extent (based on actual lat/lon coordinates)
        minx, miny, maxx, maxy = model_region_gdf.bounds.minx.item(), model_region_gdf.bounds.miny.item(), model_region_gdf.bounds.maxx.item(), model_region_gdf.bounds.maxy.item()
        ax.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

        # Set title and figure annotations
        ax.set_title(title, fontsize=10)
        ax.text(0, 1.02, subplot_labels[i], transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='bottom', ha='left')

        # Add gridlines and format tick labels
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.xlocator = mticker.FixedLocator(np.arange(minx, maxx + 0.1, 0.2))
        gl.ylocator = mticker.FixedLocator(np.arange(miny, maxy + 0.1, 0.2))
        gl.xformatter = mticker.FuncFormatter(lon_formatter)
        gl.yformatter = mticker.FuncFormatter(lat_formatter)
        gl.right_labels = False
        gl.top_labels = False
        gl.xlabel_style = {'size': 9}
        gl.ylabel_style = {'size': 9}
        if i != 0: 
            gl.left_labels = False

    # === Final touches ===
    cbar = fig.colorbar(im, ax=axes, orientation="vertical", shrink=0.5, pad=0.01)
    cbar.set_label('Attributable flood depth (m)', labelpad=10, fontsize=9)
    cbar.ax.tick_params(labelsize=9)
    cbar.set_ticks(np.arange(0, 0.7, 0.2))

    fig.savefig("../figures/f04.png", bbox_inches='tight', dpi=300)
    fig.savefig("../figures/f04.pdf", bbox_inches='tight', dpi=300)
    

def plot_driver_combination_volume_extent_damage(sfincs_models, fiat_models, filter_keys=None, tolerance=1e-6):
    model_dict = {}
    for sf, fiat in zip(sfincs_models, fiat_models):
        # Percent changes
        vol = sf['sfincs_results'].get("Volume_diff_from_F(%)", None)
        ext = sf['sfincs_results'].get("Extent_diff_from_F(%)", None)
        dam = fiat.get('Damage_diff_from_F(%)', None)
        
        if vol is None or ext is None or dam is None:
            continue
        vol = float(vol)
        ext = float(ext)
        dam = float(dam)

        CF_info = sf["CF_info"]
        drivers = [k.upper() for k in ["rain", "wind", "SLR"] if CF_info.get(k, 0) != 0]
        key = tuple(sorted(drivers)) if drivers else ("FACTUAL",)
        
        model_dict[key] = {
            'extent_pct': ext,
            'volume_pct': vol,
            'damage_pct': dam
        }

    # Build and sort data
    all_keys = sorted(model_dict.keys(), key=lambda k: (len(k), k))
    data_plot = []
    for key in all_keys:
        vals = model_dict[key]
        # Skip if all changes are below tolerance
        if all(abs(vals[v]) < tolerance for v in ['volume_pct', 'extent_pct', 'damage_pct']):
            continue
        label = format_driver_label(key)
        data_plot.append({'label': label, 'key': key, **vals})

    data_plot.sort(key=lambda d: d['damage_pct'])

    # Apply filtering
    if filter_keys is not None:
        filter_keys_normalized = [tuple(sorted(fk.split(" & "))) if isinstance(fk, str) else tuple(sorted(fk)) for fk in filter_keys]
        data_plot = [d for d in data_plot if tuple(sorted(d['key'])) in filter_keys_normalized]

    if not data_plot:
        print("No data available for the selected driver combinations.")
        return

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    x = np.arange(len(data_plot))
    width = 0.2

    max_pct = 0
    for i, d in enumerate(data_plot):
        # Bars for % changes
        ax.bar(x[i] - width, d['extent_pct'], width=width, color="#384860", edgecolor='black')
        ax.bar(x[i], d['volume_pct'], width=width, color="#5a7d9a", edgecolor='black')
        ax.bar(x[i] + width, d['damage_pct'], width=width, color="#c34a36", edgecolor='black')

        # Annotate % change bars
        if d['extent_pct'] != 0:
            ax.text(x[i] - width, d['extent_pct'] + 0.35, format_pct(d['extent_pct']), ha='center', va='bottom', fontsize=13)
        if d['volume_pct'] != 0:
            ax.text(x[i], d['volume_pct'] + 0.35, format_pct(d['volume_pct']), ha='center', va='bottom', fontsize=13)
        if d['damage_pct'] != 0:
            ax.text(x[i] + width, d['damage_pct'] + 0.35, format_pct(d['damage_pct']), ha='center', va='bottom', fontsize=13)

        max_pct = max(max_pct, abs(d['volume_pct']), abs(d['extent_pct']), abs(d['damage_pct']))

    # Layout settings
    ax.set_xticks(x)
    ax.set_xticklabels([d['label'] for d in data_plot], fontsize=14)
    ax.set_ylabel("Attributable relative change (%)", fontsize=16)
    ax.set_ylim(0, max_pct*1.15)
    ax.set_xlim(-0.5, len(data_plot) - 0.5)
    ax.tick_params(axis='y', labelsize=14)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)

    legend_elements = [
        Patch(facecolor='#384860', edgecolor='black', label='Flood extent'),
        Patch(facecolor='#5a7d9a', edgecolor='black', label='Flood volume'),
        Patch(facecolor='#c34a36', edgecolor='black', label='Damage')
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1), fontsize=14)

    plt.tight_layout()
    plt.savefig("../figures/f05.png", dpi=300, bbox_inches="tight")
    plt.savefig("../figures/f05.pdf", dpi=300, bbox_inches="tight")
    

def table_abs_and_rel_vol_ext_dam(sfincs_models, fiat_models):
    usd_2010_to_2019 = 1.172 # Convert US-Dollars (2010) to US-Dollars (2019) - annual averages: 255.657 / 218.056
    
    data_dict = {}
    factual_data = None

    for i, (sf, fiat) in enumerate(zip(sfincs_models, fiat_models)):

        category = sf.get('category', '').lower()
        if category == 'factual':
            print("passing factual data")
            factual_data = {
                "vol_abs": float(sf['sfincs_results'].get("flood_volume_m3", None)),
                "ext_abs": float(sf['sfincs_results'].get("flood_extent_km2", None)),
                "dam_abs": float(fiat['fiat_results'].get("total_damage", None).sum()) * usd_2010_to_2019
            }
            continue  # skip adding factual to data_dict

        # Extract metrics
        vol_abs = sf['sfincs_results'].get("flood_volume_m3", None)
        ext_abs = sf['sfincs_results'].get("flood_extent_km2", None)
        dam_abs = fiat['fiat_results'].get("total_damage", None).sum() * usd_2010_to_2019

        vol_pct = sf['sfincs_results'].get("Volume_diff_from_F(%)", None)
        ext_pct = sf['sfincs_results'].get("Extent_diff_from_F(%)", None)
        dam_pct = fiat.get("Damage_diff_from_F(%)", None)

        try:
            vol_abs = float(vol_abs)
            ext_abs = float(ext_abs)
            dam_abs = float(dam_abs)
            vol_pct = float(vol_pct)
            ext_pct = float(ext_pct)
            dam_pct = float(dam_pct)
        except Exception:
            continue

        CF_info = sf.get("CF_info", {})
        drivers = tuple(sorted(k.upper() for k in ["rain", "wind", "SLR"] if CF_info.get(k, 0) != 0)) 
        # Store counterfactual data keyed by drivers
        data_dict[drivers] = {
            "vol_abs": vol_abs,
            "ext_abs": ext_abs,
            "dam_abs": dam_abs,
            "vol_pct": vol_pct,
            "ext_pct": ext_pct,
            "dam_pct": dam_pct
        }

    # Prepare final keys, start with factual
    sorted_keys = [("FACTUAL",)] + sorted(data_dict.keys(), key=lambda k: (len(k), k))
    labels = [format_driver_label(k) for k in sorted_keys]

    # Build rows: factual first, then counterfactuals
    table_data = []

    # Factual row
    table_data.append([
        f"{factual_data['vol_abs']/1e6:.3f}", "-",  # no percent change for factual
        f"{factual_data['ext_abs']:.3f}", "-",
        f"{factual_data['dam_abs']/1e6:.3f}", "-"
    ])

    # Counterfactual rows
    for key in sorted_keys[1:]:
        d = data_dict[key]
        row = [
            f"{d['vol_abs']/1e6:.3f}", f"{d['vol_pct']:.3f}",
            f"{d['ext_abs']:.3f}", f"{d['ext_pct']:.3f}",
            f"{d['dam_abs']/1e6:.3f}", f"{d['dam_pct']:.3f}"
        ]
        table_data.append(row)

    # Now create the table with 6 columns
    fig, ax = plt.subplots(figsize=(13, 0.7 + 0.45 * len(table_data)))
    ax.axis("off")

    table = ax.table(
        cellText=table_data,
        rowLabels=labels,
        colLabels=[
            "Flood Volume \n[10⁶ m³]", "Volume Δ \n[%]",
            "Flood Extent \n[km²]", "Extent Δ \n[%]",
            "Damage \n[10⁶ USD]", "Damage Δ \n[%]"
        ],
        cellLoc='center',
        rowLoc='center',
        loc='center'
    )

    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.15, 1.4)

    plt.tight_layout()

    columns = [
        "Flood Volume [10^6 m³]", "Volume Δ [%]",
        "Flood Extent [km²]", "Extent Δ [%]",
        "Damage [10^6 USD]", "Damage Δ [%]"
    ]

    df = pd.DataFrame(table_data, columns=columns)
    df.insert(0, "Counterfactual", labels)

    df.to_csv("../figures/TS2.csv", index=False)
    fig.savefig("../figures/TS2.png", dpi=300, bbox_inches="tight")
    print("✅ Saved table as 'TS2.csv/.png'")


def plot_cf_timeseries_from_models(models, stations_list=[5, 40], gauges_list=[1,2]):
    colors = plt.get_cmap('tab10').colors  

    # Extract single-driver models
    rain_model = get_single_driver_model("rain", models)
    slr_model = get_single_driver_model("SLR", models)
    wind_model = get_single_driver_model("wind", models)

    if not all([rain_model, slr_model, wind_model]):
        print("Missing single-driver model (rain, SLR, or wind) with 'hmax_diff'.")
        return

    counterfactuals = {
        "SLR": slr_model,
        "Wind": wind_model,
        "Rain": rain_model
    }

    # Create figure
    fig, axs = plt.subplots(4, 1, figsize=(6, 8), dpi=300, constrained_layout=True)
    fig.suptitle("Factual vs. Counterfactual Forcing")

    # Define plotting ranges
    start, end = "2019-03-14", "2019-03-16"
    start_dis, end_dis = "2019-03-17", "2019-03-22"

    # -------------------
    # 1. SLR subplot
    axs[0].plot(models[0]['sfincs_model'].forcing['bzs'].time.sel(time=slice(start, end)), 
                models[0]['sfincs_model'].forcing['bzs'].sel(index=stations_list[1]).sel(time=slice(start, end)),
                color=colors[1], label=f'S{stations_list[1]} Factual')
    axs[0].plot(slr_model['sfincs_model'].forcing['bzs'].time.sel(time=slice(start, end)), 
                slr_model['sfincs_model'].forcing['bzs'].sel(index=stations_list[1]).sel(time=slice(start, end)),
                color="#206AAF", linestyle='-', label=f'S{stations_list[1]} -0.14 cm SLR')
    axs[0].set_ylabel("Water level\n[m]")
    axs[0].set_title("CF SLR")

    # -------------------
    # 2. Wind subplot
    axs[1].plot(models[0]['sfincs_model'].forcing['bzs'].time.sel(time=slice(start, end)), 
                models[0]['sfincs_model'].forcing['bzs'].sel(index=stations_list[1]).sel(time=slice(start, end)),
                color=colors[1], label=f'S{stations_list[1]} Factual')
    axs[1].plot(wind_model['sfincs_model'].forcing['bzs'].time.sel(time=slice(start, end)), 
                wind_model['sfincs_model'].forcing['bzs'].sel(index=stations_list[1]).sel(time=slice(start, end)),
                color="#36CEC6", linestyle='-', label=f'S{stations_list[1]} -10% Wind')
    axs[1].set_ylabel("Water level\n[m]")
    axs[1].set_title("CF Wind")

    # -------------------
    # 3. Rain - discharge subplot
    axs[2].plot(models[0]['sfincs_model'].forcing['dis'].time.sel(time=slice(start_dis, end_dis)), 
                models[0]['sfincs_model'].forcing['dis'].sel(index=gauges_list[0]).sel(time=slice(start_dis, end_dis)),
                color=colors[2], label=f'G{gauges_list[0]} Factual')
    axs[2].plot(models[0]['sfincs_model'].forcing['dis'].time.sel(time=slice(start_dis, end_dis)), 
                models[0]['sfincs_model'].forcing['dis'].sel(index=gauges_list[1]).sel(time=slice(start_dis, end_dis)),
                color=colors[3], label=f'G{gauges_list[1]} Factual')

    axs[2].plot(rain_model['sfincs_model'].forcing['dis'].time.sel(time=slice(start_dis, end_dis)), 
                rain_model['sfincs_model'].forcing['dis'].sel(index=gauges_list[0]).sel(time=slice(start_dis, end_dis)),
                color=colors[2], linestyle='--', label=f'G{gauges_list[0]} -8% Rain')
    axs[2].plot(rain_model['sfincs_model'].forcing['dis'].time.sel(time=slice(start_dis, end_dis)), 
                rain_model['sfincs_model'].forcing['dis'].sel(index=gauges_list[1]).sel(time=slice(start_dis, end_dis)),
                color=colors[3], linestyle='--', label=f'G{gauges_list[1]} -8% Rain')
    axs[2].set_ylabel("Discharge\n[m³/s]")
    axs[2].set_title("CF Rain")

    # -------------------
    # 4. Rain - precipitation subplot
    axs[3].step(models[0]['sfincs_model'].forcing['precip_2d'].time.sel(time=slice(start, end)),
                models[0]['sfincs_model'].forcing['precip_2d'].sum(dim=["x", "y"]).sel(time=slice(start, end)),
                where='post', color=colors[0], label='Factual')

    axs[3].step(rain_model['sfincs_model'].forcing['precip_2d'].time,
                rain_model['sfincs_model'].forcing['precip_2d'].sum(dim=["x", "y"]),
                where='post', color=colors[4], linestyle='-', label='-8% Rain')
    axs[3].set_ylabel("Accum. precipitation\n[mm/h]")
    axs[3].set_xlabel("Day in March 2019")
    axs[3].set_title("CF Rain")

    # -------------------
    # Style all axes
    for i, ax in enumerate(axs):
        ax.grid(True, linestyle="--", alpha=0.6)
        ax.legend(fontsize=8, loc="upper right")
        ax.tick_params(labelsize=9)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
        # Add subplot label (a), (b), ...
        ax.text(0.02, 0.9, f"({chr(97+i)})", transform=ax.transAxes, fontsize=10, fontweight="bold")

    # set xlim for different plots
    axs[0].set_xlim([np.datetime64(start), np.datetime64(end)])
    axs[1].set_xlim([np.datetime64(start), np.datetime64(end)])
    axs[2].set_xlim([np.datetime64(start_dis), np.datetime64(end_dis)])
    axs[3].set_xlim([np.datetime64(start), np.datetime64(end)])

    # Save
    fig.savefig("../figures/fS12.png", bbox_inches="tight", dpi=300)
    fig.savefig("../figures/fS12.pdf", bbox_inches="tight", dpi=300)



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
# Read model region
model_region = gpd.read_file(join(prefix, "11210471-001-compass", "03_Runs", "sofala", "Idai", "sfincs", 
                                  "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0",
                                  "gis", "region.geojson"))

#%%
# Calculate hmax and mask out permanent water
gwso = gwso_sfincs_region(models[0])
models, gdf_valid = compute_hmax_masked(models, gwso, model_region)
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
# PLOTTING for paper
# Figure 4
plot_hmax_diff_rain_slrwind_all(models, model_region, gdf_valid)

#%%
# Table 2 & S2
table_abs_and_rel_vol_ext_dam(models, fiat_models)

#%%
# Figure 5
plot_driver_combination_volume_extent_damage(models, fiat_models, filter_keys=["RAIN", "SLR & WIND", "RAIN & SLR & WIND"])

#%%
# Figure S11
plot_cf_timeseries_from_models(models)

# %%

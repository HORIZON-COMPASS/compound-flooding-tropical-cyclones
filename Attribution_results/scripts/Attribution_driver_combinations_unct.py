#%% use pixi environment compass-wflow
# Load the necessary packages
import os
from os.path import join
import yaml
import gc
import pandas as pd
import numpy as np
import platform
import gc
from hydromt_sfincs import SfincsModel
from hydromt import DataCatalog
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs
from matplotlib.patches import Patch
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap
from shapely.geometry import box
from rasterio.features import shapes
from shapely.geometry import shape
from matplotlib.colors import LinearSegmentedColormap
import rioxarray as rxr  # Required for reading TIFF files

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


def generate_sfincs_run_combinations(config, runname='Idai'):
    """Generate the specific CF run combinations defined in your Snakefile."""
    
    run = config['runname_ids'][runname]

    rain_low, rain_high = run['CF_rain_uncert']
    rain_med = run['CF_value_rain'][1]

    wind_low, wind_high = run['CF_wind_uncert']
    wind_med = run['CF_value_wind'][1]

    slr_low, slr_high = run['CF_SLR_uncert']
    slr_med = run['CF_value_SLR'][1]

    # Start with factual
    run_combinations = [(0, 0, 0)]

    # All drivers (low / med / high)
    run_combinations.extend([
        (rain_low, slr_low, wind_low),
        (rain_med, slr_med, wind_med),
        (rain_high, slr_high, wind_high),
    ])

    # Rain only
    run_combinations.extend([
        (rain_low, 0, 0),
        (rain_med, 0, 0),
        (rain_high, 0, 0),
    ])

    # Wind only
    run_combinations.extend([
        (0, 0, wind_low),
        (0, 0, wind_med),
        (0, 0, wind_high),
    ])

    # SLR only
    run_combinations.extend([
        (0, slr_low, 0),
        (0, slr_med, 0),
        (0, slr_high, 0),
    ])

    # Wind + SLR
    run_combinations.extend([
        (0, slr_low, wind_low),
        (0, slr_med, wind_med),
        (0, slr_high, wind_high),
    ])

    # Medium mixed rain combinations
    run_combinations.extend([
        (rain_med, 0, wind_med),
        (rain_med, slr_med, 0),
    ])

    return run_combinations


# Load the SFINCS models and create a model dictonary
def load_sfincs_models(config):
    """Generates model paths and categories for SFINCS runs based on CF values."""
    run = config['runname_ids']['Idai']
    # base_path = join("..", "data", "sfincs")
    base_path = join(prefix, '11210471-001-compass/03_Runs/sofala/Idai/sfincs')
    run_combinations = generate_sfincs_run_combinations(config, runname='Idai')
    
    models = []
    factual_model = None  # Initialize factual_model before the loop

    for rain, slr, wind in run_combinations:
        model_name = f"event_tp_{run['precip_forcing']}_CF{rain}_{run['tidemodel']}_CF{slr}_{run['wind_forcing']}_CF{wind}"
        model_path = join(base_path, model_name)
        utmzone    = run['utmzone']
        model_obj  = SfincsModel(model_path, mode="r")
        # flood_file = join(model_path,"floodmap.tif")
        flood_file = join(model_path,"plot_output","floodmap.tif")
        model_obj.results['flood_file'] = flood_file

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
    # base_path = os.path.join('..', 'data', "fiat")
    base_path = join(prefix, '11210471-001-compass/03_Runs/sofala/Idai/fiat')
    run_combinations = generate_sfincs_run_combinations(config, runname='Idai')

    models = []
    factual_model = None  # Initialize factual_model before the loop

    for rain, slr, wind in run_combinations:
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
def gwso_sfincs_region(model, region):
    if platform.system() == "Windows":
        datacat_path = os.path.abspath("../../Workflows/03_data_catalogs/datacatalog_general.yml")
    else:
        datacat_path = os.path.abspath("../../Workflows/03_data_catalogs/datacatalog_general___linux.yml")
    data_catalog = DataCatalog(data_libs = [datacat_path])
    sfincs_region = region
    gwso_region = data_catalog.get_rasterdataset("gswo", geom=sfincs_region, buffer=1000)
    return gwso_region


# Compute the maximum water level (hmax) and mask out permanent water
def compute_hmax_masked(models, gwso_region, model_region_gdf):
    for model in models:
        # select the highest-resolution elevation dataset
        print(f"Processing model: {model['model_name']}")
       
        # Load flood map which is already downscaled, masked for permanent water and represents a cells as flooded from 0.05 m or more
        print("Loading preprocessed flood map...")
        flood_file = model['sfincs_results']['flood_file']
        hmax = rxr.open_rasterio(flood_file)

        # Remove band dimension if present (common with TIFF files)
        if "band" in hmax.dims:
            hmax = hmax.squeeze("band", drop=True)
        
        # Add the name attribute for identification
        model["sfincs_results"]['hmax_masked'] = hmax    

        # Create background polygon that aligns with the permanent water mask for the model region
        # GSWO dataset to mask permanent water by first geprojecting it to the subgrid of hmax
        gswo_mask = gwso_region.raster.reproject_like(hmax, method="max")    

        # permanent water where water occurence > 5%
        valid_mask = (gswo_mask <= 5).astype("uint8").squeeze()

        # Extract shapes
        shapes_gen = shapes(valid_mask.values, transform=valid_mask.rio.transform())
        valid_polygons = [shape(geom) for geom, val in shapes_gen if val == 1]
        gdf_valid = gpd.GeoDataFrame(geometry=valid_polygons, crs=gswo_mask.rio.crs)
        gdf_valid = gdf_valid.to_crs(model_region_gdf.crs)

        del hmax, gswo_mask  # Clean up to free memory
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
    # Get models for specific scenarios
    scenario_filters = [
    {"name": "Rain",        "CF_rain": -8,    "CF_SLR": 0,      "CF_wind": 0},
    {"name": "SLR & Wind",  "CF_rain": 0,     "CF_SLR": -0.1,   "CF_wind": -5},
    {"name": "All",         "CF_rain": -8,    "CF_SLR": -0.1,   "CF_wind": -5},
    ]

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

    for i, (ax, scenario) in enumerate(zip(axes, scenario_filters)):
        # Select the model that matches the scenario
        model = next((m for m in models
                    if m["CF_info"].get("rain",0) == scenario["CF_rain"] and
                        m["CF_info"].get("SLR",0)  == scenario["CF_SLR"] and
                        m["CF_info"].get("wind",0) == scenario["CF_wind"]), None)
        if model is None:
            print(f"Missing model for scenario {scenario['name']}")
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
        ax.set_title(scenario["name"], fontsize=10)
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
    

def plot_driver_combination_volume_extent_damage(sfincs_models, fiat_models):
    # Define the 3 medium-value scenarios
    medium_runs = [
        {'rain': -8, 'wind': 0, 'SLR': 0},                # Rain only
        {'rain': 0, 'wind': -5, 'SLR': -0.1},            # SLR & Wind
        {'rain': -8, 'wind': -5, 'SLR': -0.1},           # All combined
    ]
    
    scenario_labels = ["Rain", "SLR & Wind", "All"]
    
    # Collect data for plotting
    plot_data = []
    
    for medium_cf in medium_runs:
        # Find medium run
        medium_sf = next(sf for sf in sfincs_models if all(
            sf['CF_info'].get(k, 0) == v for k, v in medium_cf.items()))
        medium_fiat = next(fiat for fiat in fiat_models if fiat['model_name'] == medium_sf['model_name'])
        
        # Find low and high uncertainty runs for this driver combination
        driver_keys = [k for k, v in medium_cf.items() if v != 0]
        low_cf = {k: 0 for k in ['rain', 'wind', 'SLR']}
        high_cf = {k: 0 for k in ['rain', 'wind', 'SLR']}
        # Assign low/high only for drivers involved
        for k in driver_keys:
            if k == 'rain':
                low_cf[k] = -4
                high_cf[k] = -16
            elif k == 'wind':
                low_cf[k] = -1
                high_cf[k] = -10
            elif k == 'SLR':
                low_cf[k] = -0.05
                high_cf[k] = -0.15
        
        # Get low/high runs
        low_sf = next(sf for sf in sfincs_models if all(
            sf['CF_info'].get(k, 0) == v for k, v in low_cf.items()))
        high_sf = next(sf for sf in sfincs_models if all(
            sf['CF_info'].get(k, 0) == v for k, v in high_cf.items()))
        low_fiat  = next(fiat for fiat in fiat_models if fiat['model_name'] == low_sf['model_name'])
        high_fiat = next(fiat for fiat in fiat_models if fiat['model_name'] == high_sf['model_name'])
        
        # Extract values
        def get_vals(sf, fiat):
            return {
                'extent': float(sf['sfincs_results'].get("Extent_diff_from_F(%)", 0)),
                'volume': float(sf['sfincs_results'].get("Volume_diff_from_F(%)", 0)),
                'damage': float(fiat.get("Damage_diff_from_F(%)", 0))
            }
        
        medium_vals = get_vals(medium_sf, medium_fiat)
        low_vals    = get_vals(low_sf, low_fiat)
        high_vals   = get_vals(high_sf, high_fiat)
        
        # Compute whiskers
        yerr = {
            'extent': [[medium_vals['extent'] - low_vals['extent']], [high_vals['extent'] - medium_vals['extent']]],
            'volume': [[medium_vals['volume'] - low_vals['volume']], [high_vals['volume'] - medium_vals['volume']]],
            'damage': [[medium_vals['damage'] - low_vals['damage']], [high_vals['damage'] - medium_vals['damage']]]
        }
        
        plot_data.append({'label': None, 'medium': medium_vals, 'yerr': yerr})
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    x = np.arange(len(plot_data))
    width = 0.2
    ax.set_axisbelow(True)    
    ax.grid(True, linestyle='--', alpha=0.6)    
    
    max_val = 0
    for i, d in enumerate(plot_data):
        ax.bar(x[i] - width, d['medium']['extent'], width=width, color="#384860", edgecolor='#555555', 
               yerr=d['yerr']['extent'], capsize=5, ecolor='darkgrey', error_kw=dict(elinewidth=1.2,
                                                                                     capthick=1.2, ecolor='darkgrey'))
        ax.bar(x[i], d['medium']['volume'], width=width, color="#5a7d9a", edgecolor='#555555', 
               yerr=d['yerr']['volume'], capsize=5, ecolor='darkgrey', error_kw=dict(elinewidth=1.2,
                                                                                     capthick=1.2, ecolor='darkgrey'))
        ax.bar(x[i] + width, d['medium']['damage'], width=width, color="#c34a36", edgecolor='#555555', 
               yerr=d['yerr']['damage'], capsize=5, ecolor='darkgrey', error_kw=dict(elinewidth=1.2,
                                                                                     capthick=1.2, ecolor='darkgrey'))
        
        max_val = max(max_val, 
                      d['medium']['extent'] + d['yerr']['extent'][1][0],
                      d['medium']['volume'] + d['yerr']['volume'][1][0],
                      d['medium']['damage'] + d['yerr']['damage'][1][0])

        # Annotate % change bars
        if d['medium']['extent'] != 0:
            ax.text(x[i] - width, d['medium']['extent'] + 0.35, format_pct(d['medium']['extent']), ha='center', va='bottom', fontsize=13)
        if d['medium']['volume'] != 0:
            ax.text(x[i], d['medium']['volume'] + 0.35, format_pct(d['medium']['volume']), ha='center', va='bottom', fontsize=13)
        if d['medium']['damage'] != 0:
            ax.text(x[i] + width, d['medium']['damage'] + 0.35, format_pct(d['medium']['damage']), ha='center', va='bottom', fontsize=13)

    # Labels and ticks
    ax.set_xticks(x)
    ax.set_xticklabels(scenario_labels, fontsize=14)
    ax.set_ylabel("Attributable relative change (%)", fontsize=16)
    ax.set_ylim(0, max_val * 1.15)
    ax.tick_params(axis='y', labelsize=14)
    
    # Legend
    legend_elements = [
        Patch(facecolor='#384860', edgecolor='black', label='Flood extent'),
        Patch(facecolor='#5a7d9a', edgecolor='black', label='Flood volume'),
        Patch(facecolor='#c34a36', edgecolor='black', label='Damage')
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1), fontsize=14)
    
    
    plt.tight_layout()
    plt.savefig("../figures/f05.png", dpi=300, bbox_inches="tight")
    plt.savefig("../figures/f05.pdf", dpi=300, bbox_inches="tight")


def plot_driver_combination_absolute(sfincs_models, fiat_models):
    # Define conversion factor from 2010 euros to 2019 USD
    usd_2010_to_2019 = 1.172 # Convert US-Dollars (2010) to US-Dollars (2019) - annual averages: 255.657 / 218.056 (https://www.bls.gov/cpi/tables/supplemental-files/)

    medium_runs = [
        {'rain': 0, 'wind': 0,  'SLR': 0},       # Factual
        {'rain': -8, 'wind': 0,  'SLR': 0},      # Rain only
        {'rain': 0,  'wind': -5, 'SLR': -0.1},   # SLR & Wind
        {'rain': -8, 'wind': -5, 'SLR': -0.1},   # All combined
    ]

    scenario_labels = ["Factual", "Rain", "SLR & Wind", "All"]

    plot_data = []

    for medium_cf in medium_runs:
        medium_sf = next(sf for sf in sfincs_models if all(
            sf['CF_info'].get(k, 0) == v for k, v in medium_cf.items()))
        medium_fiat = next(f for f in fiat_models if f['model_name'] == medium_sf['model_name'])

        is_factual = all(v == 0 for v in medium_cf.values())

        def get_abs_vals(sf, fiat):
            return {
                'extent': float(sf['sfincs_results'].get("flood_extent_km2", 0)),
                'volume': float(sf['sfincs_results'].get("flood_volume_m3", 0) / 1e6),  # M m³
                'damage': float(fiat['fiat_results'].get("total_damage", 0).sum()/1e6 * usd_2010_to_2019),  # M USD
            }

        # Get medium values
        medium_vals = get_abs_vals(medium_sf, medium_fiat)

        if not is_factual:
            # Determine drivers
            driver_keys = [k for k, v in medium_cf.items() if v != 0]

            low_cf  = {k: 0 for k in ['rain', 'wind', 'SLR']}
            high_cf = {k: 0 for k in ['rain', 'wind', 'SLR']}

            for k in driver_keys:
                if k == 'rain':
                    low_cf[k], high_cf[k] = -16, -4
                elif k == 'wind':
                    low_cf[k], high_cf[k] = -10, -1
                elif k == 'SLR':
                    low_cf[k], high_cf[k] = -0.15, -0.05

            low_sf  = next(sf for sf in sfincs_models if all(sf['CF_info'].get(k, 0) == v for k, v in low_cf.items()))
            high_sf = next(sf for sf in sfincs_models if all(sf['CF_info'].get(k, 0) == v for k, v in high_cf.items()))
            low_fiat  = next(f for f in fiat_models if f['model_name'] == low_sf['model_name'])
            high_fiat = next(f for f in fiat_models if f['model_name'] == high_sf['model_name'])

            low_vals  = get_abs_vals(low_sf, low_fiat)
            high_vals = get_abs_vals(high_sf, high_fiat)

            # Absolute asymmetric error
            yerr = {}
            for k in ['extent', 'volume', 'damage']:
                v_med = medium_vals[k]
                v_min = min(low_vals[k], high_vals[k])
                v_max = max(low_vals[k], high_vals[k])
                yerr[k] = [[v_med - v_min], [v_max - v_med]]
        else:
            # Factual has no uncertainty
            yerr = {k: None for k in ['extent', 'volume', 'damage']}

        plot_data.append({'medium': medium_vals, 'yerr': yerr, 'is_factual': is_factual})


    # Plotting
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), dpi=300, sharex=True)

    metrics = [
        ('extent', 'Flood extent (km²)', '#384860', '(a)'),
        ('volume', 'Flood volume (M m³)', '#5a7d9a', '(b)'),
        ('damage', 'Damage (M USD)', '#c34a36', '(c)')
    ]

    x = np.arange(len(plot_data))
    err_kw = dict(elinewidth=1.2, capthick=1.2, ecolor='darkgrey')

    for ax, (key, ylabel, color, panel) in zip(axes, metrics):
        ax.set_axisbelow(True)
        ax.grid(True, linestyle='--', alpha=0.6)

        max_val = 0
        for i, d in enumerate(plot_data):
            ax.bar(
                x[i],
                d['medium'][key],
                width=0.55,
                color='lightgrey' if d['is_factual'] else color,
                edgecolor='#555555',
                yerr=None if d['is_factual'] else d['yerr'][key],
                capsize=5,
                error_kw=err_kw,
                zorder=3
            )
            # update max_val for y-axis
            if d['is_factual']:
                max_val = max(max_val, d['medium'][key])
            else:
                max_val = max(max_val, d['medium'][key] + d['yerr'][key][1][0])

        ax.set_xticks(x)
        ax.set_xticklabels(scenario_labels, fontsize=13)
        ax.set_ylabel(ylabel, fontsize=14)
        ax.set_ylim(0, max_val * 1.15)
        ax.tick_params(axis='y', labelsize=12)

        # Panel label
        ax.text(0.02, 1.07, panel, transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')

    plt.tight_layout()
    plt.show()
    plt.savefig("../figures/f05_abs.png", dpi=300, bbox_inches="tight")
    plt.savefig("../figures/f05_abs.pdf", dpi=300, bbox_inches="tight")
    plt.close()


def table_abs_and_rel_vol_ext_dam(sfincs_models, fiat_models):
    usd_2010_to_2019 = 1.172
    plot_data = []

    # Define the scenarios you want in order
    medium_runs = [
        {'rain': 0,  'wind': 0,  'SLR': 0},      # Factual
        {'rain': -8, 'wind': 0,  'SLR': 0},      # Rain only
        {'rain': 0,  'wind': -5, 'SLR': 0},      # Wind only
        {'rain': 0,  'wind': 0,  'SLR': -0.1},   # SLR only
        {'rain': 0,  'wind': -5, 'SLR': -0.1},   # SLR & Wind
        {'rain': -8, 'wind': -5, 'SLR': -0.1},   # All combined
    ]
    scenario_labels = ["Factual", "Rain", "Wind", "SLR", "SLR & Wind", "All"]

    # Helper to extract absolute and % difference values
    def get_vals(sf, fiat):
        return {
            'extent_abs': float(sf['sfincs_results'].get("flood_extent_km2", 0)),
            'volume_abs': float(sf['sfincs_results'].get("flood_volume_m3", 0)/1e6),  # million m³
            'damage_abs': float(fiat['fiat_results'].get("total_damage", 0).sum()/1e6 * usd_2010_to_2019),
            'volume_pct': float(sf['sfincs_results'].get("Volume_diff_from_F(%)", 0)),
            'extent_pct': float(sf['sfincs_results'].get("Extent_diff_from_F(%)", 0)),
            'damage_pct': float(fiat.get("Damage_diff_from_F(%)", 0))
        }

    for cf in medium_runs:
        is_factual = all(v == 0 for v in cf.values())
        medium_sf = next(sf for sf in sfincs_models if all(sf['CF_info'].get(k, 0) == v for k, v in cf.items()))
        medium_fiat = next(f for f in fiat_models if f['model_name'] == medium_sf['model_name'])
        medium_vals = get_vals(medium_sf, medium_fiat)

        # Default low/high
        low_vals = high_vals = None

        if not is_factual:
            driver_keys = [k for k, v in cf.items() if v != 0]
            low_cf = {k: 0 for k in ['rain', 'wind', 'SLR']}
            high_cf = {k: 0 for k in ['rain', 'wind', 'SLR']}
            for k in driver_keys:
                if k == 'rain':
                    low_cf[k], high_cf[k] = -16, -4
                elif k == 'wind':
                    low_cf[k], high_cf[k] = -10, -1
                elif k == 'SLR':
                    low_cf[k], high_cf[k] = -0.15, -0.05

            # Find low/high runs
            low_sf = next(sf for sf in sfincs_models if all(sf['CF_info'].get(k, 0) == v for k, v in low_cf.items()))
            high_sf = next(sf for sf in sfincs_models if all(sf['CF_info'].get(k, 0) == v for k, v in high_cf.items()))
            low_fiat = next(f for f in fiat_models if f['model_name'] == low_sf['model_name'])
            high_fiat = next(f for f in fiat_models if f['model_name'] == high_sf['model_name'])

            low_vals = get_vals(low_sf, low_fiat)
            high_vals = get_vals(high_sf, high_fiat)

        plot_data.append({
            'medium': medium_vals,
            'low': low_vals,
            'high': high_vals,
            'is_factual': is_factual
        })

    # --- Build table rows ---
    table_rows = []
    for d in plot_data:
        row = []

        # Absolute values
        for metric_abs in ['volume_abs', 'extent_abs', 'damage_abs']:
            med = d['medium'][metric_abs]
            if d['low'] and d['high']:
                low = min(d['low'][metric_abs], d['high'][metric_abs])
                high = max(d['low'][metric_abs], d['high'][metric_abs])
                row.append(f"{med:.0f} ({low:.0f} - {high:.0f})")
            else:
                row.append(f"{med:.0f}")

        # Relative % values
        for metric_pct, metric_abs in [('volume_pct', 'volume_abs'),
                                       ('extent_pct', 'extent_abs'),
                                       ('damage_pct', 'damage_abs')]:
            med_pct = d['medium'][metric_pct]
            if d['is_factual']:
                row.append("-")
            elif d['low'] and d['high']:
                low_pct = d['low'][metric_pct]
                high_pct = d['high'][metric_pct]
                row.append(f"{med_pct:.2f} ({min(low_pct, high_pct):.2f} - {max(low_pct, high_pct):.2f})")
            else:
                row.append(f"{med_pct:.2f}")
        table_rows.append(row)

    # --- Plot table ---
    fig, ax = plt.subplots(figsize=(14, 0.7 + 0.45 * len(table_rows)))
    ax.axis("off")
    table = ax.table(
        cellText=table_rows,
        rowLabels=scenario_labels,
        colLabels=["Flood Volume [10⁶ m³]", "Flood Extent [km²]", "Damage [10⁶ USD]",
                   "Volume Δ [%]", "Extent Δ [%]", "Damage Δ [%]"],
        cellLoc='center',
        rowLoc='center',
        loc='center'
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.15, 1.4)
    plt.tight_layout()
    fig.savefig("../figures/TS2.png", dpi=300, bbox_inches="tight")
    print("Table saved as 'TS2.png'")

    # --- Save CSV ---
    columns = ["Flood Volume [10^6 m³]", "Flood Extent [km²]", "Damage [10^6 USD]",
               "Volume Δ [%]", "Extent Δ [%]", "Damage Δ [%]"]
    df = pd.DataFrame(table_rows, columns=columns)
    df.insert(0, "Scenario", scenario_labels)
    df.to_csv("../figures/TS2.csv", index=False)
    print("✅ Table saved as 'TS2.csv'")


def plot_cf_timeseries_all(models, stations_list=[5, 40], gauges_list=[1],
                           start="2019-03-14", end="2019-03-16",
                           start_coast="2019-03-14 16:00:00", end_coast="2019-03-15 02:00:00",
                           start_dis="2019-03-17", end_dis="2019-03-22"):
    # Define drivers
    drivers = [
        {"title": "SLR", "key": "SLR", "values": [-0.05, -0.1, -0.15]},
        {"title": "Wind", "key": "wind", "values": [-1, -5, -10]},
        {"title": "Rain", "key": "rain", "values": [-4, -8, -16]}
    ]

    # Colors for CF low, medium, high (same across all drivers)
    cf_colors = ["#1B9E77", "#D95F02", "#7570B3"]  # low, medium, high

    # -------------------
    fig, axs = plt.subplots(4, 1, figsize=(8, 9), dpi=300, constrained_layout=True)
    fig.suptitle("Factual vs Counterfactual Forcing", fontsize=14)

    # Helper: select model for a given driver and value, all others 0
    def get_model(driver_key, value):
        return next(
            (m for m in models
             if m['CF_info'].get(driver_key, 0) == value and
             all(m['CF_info'].get(k, 0) == 0 for k in ['rain', 'wind', 'SLR'] if k != driver_key)),
            None
        )

    # Factual model
    model_f = next(m for m in models if all(v == 0 for v in m['CF_info'].values()))
    t_bzs = model_f['sfincs_model'].forcing['bzs'].time.sel(time=slice(start_coast, end_coast))
    t_dis = model_f['sfincs_model'].forcing['dis'].time.sel(time=slice(start_dis, end_dis))
    t_prec = model_f['sfincs_model'].forcing['precip_2d'].time.sel(time=slice(start, end))

    # -------------------
    # Plot SLR and Wind water level subplots
    for i, drv in enumerate(drivers[:2]):  # SLR, Wind
        ax = axs[i]
        ax.set_title(f"CF {drv['title']} at Station {stations_list[1]}")
        ax.set_ylabel("Water level [m]")
        ax.grid(True, linestyle="--", alpha=0.6)

        # Factual
        bzs_factual = model_f['sfincs_model'].forcing['bzs'].sel(index=stations_list[1]).sel(time=slice(start_coast, end_coast))
        ax.plot(t_bzs, bzs_factual, color='gray', linewidth=1.5, label='Factual')

        # CF lines (low-medium-high)
        for idx, v in enumerate(drv['values']):
            model_cf = get_model(drv['key'], v)
            if model_cf is None:
                continue
            bzs = model_cf['sfincs_model'].forcing['bzs'].sel(index=stations_list[1]).sel(time=slice(start_coast, end_coast))

            # Label with units
            if drv['title'] == "SLR":
                label = f"{drv['title']} {v} m"
            elif drv['title'] == "Wind":
                label = f"{drv['title']} {v}%"

            ax.plot(t_bzs, bzs, color=cf_colors[idx], linewidth=1, label=label)

        ax.legend(fontsize=8, loc='upper right')
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %Hh'))
        ax.set_ylim(2,3.5)
        ax.tick_params(labelsize=9)
        ax.text(0.02, 0.9, f"({chr(97+i)})", transform=ax.transAxes, fontsize=10, fontweight="bold")
            
    # -------------------
    # Plot Rain: discharge and precipitation
    rain_drv = drivers[2]
    # Discharge
    ax_dis = axs[2]
    ax_dis.set_title(f"CF {rain_drv['title']} at Gauge {gauges_list[0]}")
    ax_dis.set_ylabel("Discharge [m³/s]")
    ax_dis.grid(True, linestyle="--", alpha=0.6)

    dis_factual = model_f['sfincs_model'].forcing['dis'].sel(index=gauges_list[0]).sel(time=slice(start_dis, end_dis))
    ax_dis.plot(t_dis, dis_factual, color='gray', linewidth=1.5, label='Factual')

    for idx, v in enumerate(rain_drv['values']):
        model_cf = get_model(rain_drv['key'], v)
        if model_cf is None:
            continue
        dis = model_cf['sfincs_model'].forcing['dis'].sel(index=gauges_list[0]).sel(time=slice(start_dis, end_dis))
        ax_dis.plot(t_dis, dis, color=cf_colors[idx], linewidth=1.5, label=f"{rain_drv['title']} {v} %")

    ax_dis.legend(fontsize=8, loc='upper right')
    ax_dis.xaxis.set_major_formatter(mdates.DateFormatter('%d %Hh'))
    ax_dis.tick_params(labelsize=9)
    ax_dis.text(0.02, 0.9, "(c)", transform=ax_dis.transAxes, fontsize=10, fontweight="bold")

    # Precipitation
    ax_prec = axs[3]
    ax_prec.set_title(f"CF {rain_drv['title']}")
    ax_prec.set_ylabel("Accum. precipitation [mm/h]")
    ax_prec.set_xlabel("Day and hour in March 2019")
    ax_prec.grid(True, linestyle="--", alpha=0.6)

    prec_factual = model_f['sfincs_model'].forcing['precip_2d'].sum(dim=['x','y']).sel(time=slice(start, end))
    ax_prec.step(t_prec, prec_factual, where='post', color='gray', linewidth=1.5, label='Factual')

    for idx, v in enumerate(rain_drv['values']):
        model_cf = get_model(rain_drv['key'], v)
        if model_cf is None:
            continue
        prec = model_cf['sfincs_model'].forcing['precip_2d'].sum(dim=['x','y']).sel(time=slice(start, end))
        ax_prec.step(t_prec, prec, where='post', color=cf_colors[idx], linewidth=1.5, label=f"{rain_drv['title']} {v} %")

    ax_prec.legend(fontsize=8, loc='upper right')
    ax_prec.tick_params(labelsize=9)
    ax_prec.xaxis.set_major_formatter(mdates.DateFormatter('%d %Hh'))
    ax_prec.text(0.02, 0.9, "(d)", transform=ax_prec.transAxes, fontsize=10, fontweight="bold")

    # -------------------
    # Set xlim
    axs[0].set_xlim([np.datetime64(start_coast), np.datetime64(end_coast)])
    axs[1].set_xlim([np.datetime64(start_coast), np.datetime64(end_coast)])
    axs[2].set_xlim([np.datetime64(start_dis), np.datetime64(end_dis)])
    axs[3].set_xlim([np.datetime64(start), np.datetime64(end)])

    # -------------------
    # Save
    fig.savefig("../figures/fS12.png", bbox_inches="tight", dpi=300)
    fig.savefig("../figures/fS12.pdf", bbox_inches="tight", dpi=300)
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
# Read model region
model_region = gpd.read_file(join("..", "data", "sfincs", "gis", "region.geojson"))

#%%
# Calculate hmax and mask out permanent water
gwso = gwso_sfincs_region(models[0], model_region)
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
# plot_hmax_diff_rain_slrwind_all(models, model_region, gdf_valid)

# Figure 5
# plot_driver_combination_volume_extent_damage(models, fiat_models)

# Figure 5 - adapted for absolute values
# plot_driver_combination_absolute(models, fiat_models)

# Table 2 & S2
# table_abs_and_rel_vol_ext_dam(models, fiat_models)

# # Figure S12
plot_cf_timeseries_all(models)



# %%
# Maximum flood difference for different scenarios
rain_low = models[4]['sfincs_results']['hmax_diff'].quantile(0.99).round(1).item()
rain_med = models[5]['sfincs_results']['hmax_diff'].quantile(0.99).round(1).item()
rain_high = models[6]['sfincs_results']['hmax_diff'].quantile(0.99).round(1).item()

slrwind_low = models[13]['sfincs_results']['hmax_diff'].quantile(0.99).round(1).item()
slrwind_med = models[14]['sfincs_results']['hmax_diff'].quantile(0.99).round(1).item()
slrwind_high = models[15]['sfincs_results']['hmax_diff'].quantile(0.99).round(1).item()

print("99th percentile of maximum flood depth difference (hmax_diff) for different scenarios:")
print(f"Rain low/med/high: {rain_low} / {rain_med} / {rain_high} m")
print(f"SLR & Wind low/med/high: {slrwind_low} / {slrwind_med} / {slrwind_high} m")
# %%

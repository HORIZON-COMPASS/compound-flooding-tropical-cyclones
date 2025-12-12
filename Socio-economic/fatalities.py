#%%
import os
import yaml
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from os.path import join
import rasterio
import geopandas as gpd
import itertools
import warnings
warnings.filterwarnings('ignore')
import platform
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from hydromt_sfincs import SfincsModel, utils
from hydromt import DataCatalog
import gc
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import matplotlib.colorbar as mcolorbar

prefix = "p:/" if platform.system() == "Windows" else "/p/"

# Reading shapefiles
# shapefile_fp = "p:/11210471-001-compass/03_Runs/sofala/Idai/sfincs/event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0/gis/region.geojson"   # replace with your region shapefile
shapefile_fp = "c:/Code/Paper_2/Data/region.geojson"
# background = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_region_background.geojson")
background = gpd.read_file("c:/Code/Paper_2/Data/sofala_region_background.geojson")
background_utm = background.to_crs("EPSG:32736")
region_utm = gpd.read_file(shapefile_fp)
region = region_utm.to_crs("EPSG:4326")
# Flood model subgrid
BASE_RUN_PATH = Path("C:/Code/Paper_1/Data_submission")
sfincs_subgrid = BASE_RUN_PATH / "sfincs" / "subgrid" / "dep_subgrid.tif"

#%%
# get extent from raster transform
def get_extent(transform, width, height):
    left = transform[2]
    right = left + width * transform[0]
    top = transform[5]
    bottom = top + height * transform[4]
    return [left, right, top, bottom]

with rasterio.open(sfincs_subgrid) as src:
    flood_grid_crs, flood_grid_transform, flood_grid_shape = src.crs, src.transform, (src.height, src.width)
    
flood_extent = get_extent(flood_grid_transform, flood_grid_shape[1], flood_grid_shape[0])

# %%
###################################################################################
############################### Flood fatalities ##################################
###################################################################################
# Eq. 3 from Jonkman et al (2008), based on Boyd (2005)

def load_raster(path):
    with rasterio.open(path) as src:
        return src.read(1), src.crs, src.transform
    
#### Load SFINCS model ####
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

# read global surface water occurance (GSWO) data to mask permanent water for the model region
def gwso_sfincs_region(model):
    if platform.system() == "Windows":
        datacat_path = os.path.abspath("../Workflows/03_data_catalogs/datacatalog_general.yml")
    else:
        datacat_path = os.path.abspath("../Workflows/03_data_catalogs/datacatalog_general___linux.yml")
    data_catalog = DataCatalog(data_libs = [datacat_path])
    sfincs_region = model["sfincs_model"].region
    gwso_region = data_catalog.get_rasterdataset("gswo", geom=sfincs_region, buffer=1000)
    return gwso_region

# compute maximum flood depth, maximum rise rate of flood depth, maximum velocity at maximum flood depth, and their masked versions excluding permanent water
def compute_hmax_masked_riserate(models, gwso_region):
    # we set a threshold to mask minimum flood depth
    hmin = 0.05

    for model in models:
        if os.path.exists(model['model_path'] + "/sfincs_derived.nc"):
            print(f"Skipping {model['model_name']} — results already exist.")
            dataset = xr.open_dataset(join(model['model_path'], "sfincs_derived.nc"))

            model["sfincs_results"]['hmax'] = dataset['hmax']
            model["sfincs_results"]['hmax_masked'] = dataset['hmax_masked']
            model["sfincs_results"]['max_rise_rate_h'] = dataset['max_rise_rate_h']
            model["sfincs_results"]['vmax_masked'] = dataset['vmax_masked']
            model["sfincs_results"]['hvmax_masked'] = dataset['hvmax_masked']
        
        else:
            print(f"Processing model: {model['model_name']}")

            # select the highest-resolution elevation dataset
            depfile = join(model["model_path"], "subgrid", "dep_subgrid.tif")
            da_dep = model["sfincs_model"].data_catalog.get_rasterdataset(depfile)

            # compute the maximum over all time steps
            # First timestep leads to incorrect diiference values for permanent water cells that are incorrectly unmasked; requires filtering.
            da_zsmax = model["sfincs_results"]["zsmax"].isel(timemax=slice(1, None)).max(dim="timemax")
            da_zs = model["sfincs_results"]["zs"]
        
            # downscale the floodmap for max water depth
            da_hmax = utils.downscale_floodmap(
                zsmax=da_zsmax,
                dep=da_dep,
                hmin=hmin,
                reproj_method = "bilinear"
                )
            
            # Calculate rise rate
            dzs_dt_h = da_zs.diff(dim="time") / (da_zs["time"].diff("time") / np.timedelta64(1, "h"))
            # Take the per-pixel maximum rise rate (over time)
            dzs_dt_max_h = dzs_dt_h.max(dim="time", skipna=True)
            # Downscale or reproject to match dep_subgrid
            dzs_dt_max_h.rio.write_crs("EPSG:32736")
            dzs_dt_max_h_sub = dzs_dt_max_h.rio.reproject_match(da_dep)
            # Mask cells where max water depth ≤ hmin
            mask_valid = da_hmax > 0.05
            dzs_dt_max_h_sub = dzs_dt_max_h_sub.where(mask_valid)
            
            # velocity at maximum water depth - also downscale and mask
            da_vmax = model['sfincs_results']['vmax'].max(dim='timemax')
            da_vmax = da_vmax.rio.reproject_match(da_dep)
            da_vmax_sub = da_vmax.where(mask_valid)

            # GSWO dataset to mask permanent water by first geprojecting it to the subgrid of hmax
            gswo_mask = gwso_region.raster.reproject_like(da_hmax, method="max")
            # permanent water where water occurence > 5%
            da_hmax_masked = da_hmax.where(gswo_mask <= 5)

            # Do the same for rise rate
            # gswo_mask = gwso_region.raster.reproject_like(dzs_dt_max_h, method="max")
            dzs_dt_max_h_masked = dzs_dt_max_h_sub.where(gswo_mask <= 5)

            # same for max velocity
            vmax_masked = da_vmax_sub.where(gswo_mask <= 5)  

            # calculate hv product
            hvmax_masked = (da_hmax_masked * vmax_masked).where((da_hmax_masked > 0) & (vmax_masked > 0))      

            # --- Save results ---
            model["sfincs_results"]['hmax'] = da_hmax
            model["sfincs_results"]['hmax_masked'] = da_hmax_masked
            model["sfincs_results"]['max_rise_rate_h'] = dzs_dt_max_h_masked
            model["sfincs_results"]['vmax_masked'] = vmax_masked
            model["sfincs_results"]['hvmax_masked'] = hvmax_masked

            dataset = xr.Dataset({
                "hmax": da_hmax,
                "hmax_masked": da_hmax_masked,
                "max_rise_rate_h": dzs_dt_max_h_masked,
                "vmax_masked": vmax_masked,
                "hvmax_masked": hvmax_masked
            })
            dataset.to_netcdf(join(model["model_path"], "sfincs_derived.nc"))

            del da_hmax, da_zsmax, da_zs, da_dep, gswo_mask  # Clean up to free memory
            gc.collect()

    return models


#%% FATALITY FUNCTIONS
# Calculate fatalities according to Jonkman et al (2008)
def fatality_fraction_Jonkman_etal2008(h, mu_N, sigma_N):
    """Compute fatality fraction FD(h) from water depth (m) using Jonkman (2008) lognormal fit."""
    from scipy.stats import norm

    h = np.maximum(h, 1e-3)  # avoid log(0)
    return norm.cdf((np.log(h) - mu_N) / sigma_N)

# Fatality function from Boyd (2005) as reported in Jonkman et al. (2008)
def fatality_curve_Boyd2005(h):
    """Fatality curve based on Boyd et al. (2005) as reported in Jonkman et al. (2008).
    
    Parameters
    ----------
    h : float or array-like
        Flood depth in meters.

    Returns
    -------
    float or array-like
        Fatality rate (between 0 and 1).
    """
    h = np.asarray(h)
    return 0.34 / (1 + np.exp(20.37 - 6.18 * h))

#%%
import numpy as np
import geopandas as gpd
from shapely.geometry import box

def aggregate_point_gdf(
    gdf,
    factor=100,
    cell_size=None,
    sum_cols=None,
    mean_cols=None,
    count=True
):
    """
    Aggregates a point-based GeoDataFrame into coarse grid polygons.

    Parameters
    ----------
    gdf : GeoDataFrame
        Must contain geometry points or x/y columns.
    factor : int
        Coarsening factor relative to fine cell size.
    cell_size : float or None
        Size of native fine grid cell (meters). If None, is estimated.
    sum_cols : list of str
        Columns to sum in each coarse cell (e.g. population, fatalities).
    mean_cols : list of str
        Columns to average (e.g. flood depth).
    count : bool
        Whether to include count of points per coarse cell.

    Returns
    -------
    GeoDataFrame
        Coarse grid polygons with aggregated values.
    """

    gdf = gdf.copy()

    # -------------------------------
    # 1. Ensure geometry exists
    # -------------------------------
    if "geometry" not in gdf.columns:
        if ("x" not in gdf.columns) or ("y" not in gdf.columns):
            raise ValueError("GDF must have geometry or x/y columns.")
        gdf["geometry"] = gpd.points_from_xy(gdf.x, gdf.y)

    # -------------------------------
    # 2. Estimate fine cell size if needed
    # -------------------------------
    if cell_size is None:
        # compute min spacing in x direction
        xs = np.sort(gdf.geometry.x.unique())
        diffs = np.diff(xs)
        diffs = diffs[diffs > 0]
        cell_size = np.median(diffs)
        print(f"Estimated fine cell size = {cell_size:.2f} m")

    coarse_size = factor * cell_size
    print(f"Coarse grid size = {coarse_size:.2f} m")

    # -------------------------------
    # 3. Build coarse grid
    # -------------------------------
    xmin, ymin, xmax, ymax = gdf.total_bounds

    xs = np.arange(xmin, xmax + coarse_size, coarse_size)
    ys = np.arange(ymin, ymax + coarse_size, coarse_size)

    grid_cells = [
        box(x, y, x + coarse_size, y + coarse_size)
        for x in xs[:-1]
        for y in ys[:-1]
    ]

    grid = gpd.GeoDataFrame({"geometry": grid_cells}, crs=gdf.crs)
    grid["grid_id"] = np.arange(len(grid))

    # -------------------------------
    # 4. Spatial join points -> grid
    # -------------------------------
    joined = gpd.sjoin(gdf, grid, predicate="within")

    # -------------------------------
    # 5. Aggregation
    # -------------------------------
    agg_dict = {}

    if sum_cols is not None:
        for col in sum_cols:
            agg_dict[col] = (col, "sum")

    if mean_cols is not None:
        for col in mean_cols:
            agg_dict[col] = (col, "mean")

    if count:
        agg_dict["n_points"] = ("geometry", "count")

    aggregated = joined.groupby("grid_id").agg(**agg_dict)

    # -------------------------------
    # 6. Reattach geometry
    # -------------------------------
    out = grid.join(aggregated, on="grid_id", how="left")

    return out



#%%
# Load exposed population dataframes
data_path = Path("data/preprocessed/population")

# Load exposed population dataframes
df_pop_2019_exposed_F  = pd.read_csv(join(data_path, "df_pop_2019_F.csv"))
df_pop_2019_exposed_CF = pd.read_csv(join(data_path, "df_pop_2019_CF.csv"))
df_pop_1990_exposed_F  = pd.read_csv(join(data_path, "df_pop_1990_F.csv"))
df_pop_1990_exposed_CF = pd.read_csv(join(data_path, "df_pop_1990_CF.csv"))
# unform pop growth scenario
df_pop_2019_exposed_F_uniform  = pd.read_csv(join(data_path, "df_pop_2019_uniform_F.csv"))
df_pop_2019_exposed_CF_uniform = pd.read_csv(join(data_path, "df_pop_2019_uniform_CF.csv"))

# into gdf
gdf_pop_2019_exposed_F = gpd.GeoDataFrame(df_pop_2019_exposed_F, geometry=gpd.points_from_xy(df_pop_2019_exposed_F["x"], df_pop_2019_exposed_F["y"]), crs="EPSG:32736")

# Load exposed population rasters
ra_exposed_pop_2019_F, crs, transform = load_raster(join(data_path, "exposed_pop_2019_F.tif"))
ra_exposed_pop_2019_F_uniform, crs, transform = load_raster(join(data_path, "exposed_pop_2019_F_uniform.tif"))
ra_exposed_pop_1990_F, crs, transform = load_raster(join(data_path, "exposed_pop_1990_F.tif"))
ra_exposed_pop_2019_CF, crs, transform = load_raster(join(data_path, "exposed_pop_2019_CF.tif"))
ra_exposed_pop_2019_CF_uniform, crs, transform = load_raster(join(data_path, "exposed_pop_2019_CF_uniform.tif"))
ra_exposed_pop_1990_CF, crs, transform = load_raster(join(data_path, "exposed_pop_1990_CF.tif"))

df_pop_2019_exposed_F['fatalities'] = df_pop_2019_exposed_F['population'] * fatality_curve_Boyd2005(df_pop_2019_exposed_F['flood_depth'])
df_pop_2019_exposed_CF['fatalities'] = df_pop_2019_exposed_CF['population'] * fatality_curve_Boyd2005(df_pop_2019_exposed_CF['flood_depth'])
df_pop_1990_exposed_F['fatalities'] = df_pop_1990_exposed_F['population'] * fatality_curve_Boyd2005(df_pop_1990_exposed_F['flood_depth'])
df_pop_1990_exposed_CF['fatalities'] = df_pop_1990_exposed_CF['population'] * fatality_curve_Boyd2005(df_pop_1990_exposed_CF['flood_depth'])


df_pop_2019_exposed_F_uniform['fatalities'] = df_pop_2019_exposed_F_uniform['population'] * fatality_curve_Boyd2005(df_pop_2019_exposed_F_uniform['flood_depth'])
df_pop_2019_exposed_CF_uniform['fatalities'] = df_pop_2019_exposed_CF_uniform['population'] * fatality_curve_Boyd2005(df_pop_2019_exposed_CF_uniform['flood_depth'])

print(f"Total fatalities acc. to Boyd et al. (2005):")
print(f"Factual climate 2019 population scenario:                {df_pop_2019_exposed_F['fatalities'].sum():.0f} people")
print(f"Factual climate 2019 UNIFORM population scenario:        {df_pop_2019_exposed_F_uniform['fatalities'].sum():.0f} people")
print(f"Counterfactual climate 2019 UNIFORM population scenario: {df_pop_2019_exposed_CF_uniform['fatalities'].sum():.0f} people")
print(f"Counterfactual climate 2019 population scenario:         {df_pop_2019_exposed_CF['fatalities'].sum():.0f} people")
print(f"Factual climate 1990 population scenario:                {df_pop_1990_exposed_F['fatalities'].sum():.0f} people")
print(f"Counterfactual climate 1990 population scenario:         {df_pop_1990_exposed_CF['fatalities'].sum():.0f} people")


#%%
# Aggregate to larger cells
gdf_pop_2019_exposed_F['fatalities'] = gdf_pop_2019_exposed_F['population'] * fatality_curve_Boyd2005(gdf_pop_2019_exposed_F['flood_depth'])
agg_pop_2019_exposed_F = aggregate_point_gdf(gdf_pop_2019_exposed_F, factor=100, sum_cols=["fatalities"])


#%%
# Load snakemake config file to construct the model paths
config_path  = '../Workflows/01_config_snakemake/config_general_MZB.yml'
cfg = load_config(config_path)

# Load the SFINCS models in one dictonary
models = load_sfincs_models(cfg)
print(models)

# Calculate hmax and mask out permanent water
gwso = gwso_sfincs_region(models[0])


#%%
# Create model object for SFINCS model incl velocity output
base_path = join(prefix, "11210471-001-compass", "03_Runs", "sofala", "Idai", "sfincs")
model_name = f"event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_maxvel"
model_path = join(base_path, model_name)
model_obj  = SfincsModel(model_path, mode="r")
his_path   = os.path.join(model_path,"sfincs_his.nc")
ds_his     = xr.open_dataset(his_path, engine="netcdf4")

# Create model dictionary
model_velocity_F = {
    "model_name": model_name,
    "model_path": model_path,
    "sfincs_model": model_obj,
    "sfincs_results": model_obj.results,
    "sfincs_his": ds_his
    }

model_velocity_F = compute_hmax_masked_riserate([model_velocity_F], gwso)

hvmax_F = model_velocity_F[0]['sfincs_results']['hvmax_masked']
vmax_F = model_velocity_F[0]['sfincs_results']['vmax_masked']
hmax_F = model_velocity_F[0]['sfincs_results']['hmax_masked']
rise_F = model_velocity_F[0]['sfincs_results']['max_rise_rate_h']

# Apply Jonkman et al. (2008) fatality function based on zones defined by hvmax, vmax, hmax, and rise rate
# Initialize with zeros
fatality_frac_F = xr.zeros_like(hmax_F)
# Define masks
F_mask_breach = (hvmax_F >= 7) & (vmax_F > 2)
F_mask_rapid = (hmax_F > 2.1) & (rise_F > 0.5)
F_mask_remaining = (hmax_F > 0) & ~F_mask_breach & ~F_mask_rapid

# Apply per zone
fatality_frac_F = xr.where(F_mask_breach, 1, fatality_frac_F)
fatality_frac_F = xr.where(F_mask_rapid, fatality_fraction_Jonkman_etal2008(hmax_F, mu_N=1.46, sigma_N=0.28), fatality_frac_F)
fatality_frac_F = xr.where(F_mask_remaining, fatality_fraction_Jonkman_etal2008(hmax_F, mu_N=7.60, sigma_N=2.75), fatality_frac_F)
# Store the result
model_velocity_F[0]['sfincs_results']['fatality_frac'] = fatality_frac_F

#%%
# Calculate exposed fatalities for exposed population to factual flooding
flood_fatalities_F_2019 = model_velocity_F[0]['sfincs_results']['fatality_frac'] * ra_exposed_pop_2019_F
flood_fatalities_F_1990 = model_velocity_F[0]['sfincs_results']['fatality_frac'] * ra_exposed_pop_1990_F
flood_fatalities_F_2019_uniform = model_velocity_F[0]['sfincs_results']['fatality_frac'] * ra_exposed_pop_2019_F_uniform

#%%
# Create model object for SFINCS model incl velocity output
base_path = join(prefix, "11210471-001-compass", "03_Runs", "sofala", "Idai", "sfincs")
model_name = f"event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10_maxvel"
model_path = join(base_path, model_name)
model_obj  = SfincsModel(model_path, mode="r")
his_path   = os.path.join(model_path,"sfincs_his.nc")
ds_his     = xr.open_dataset(his_path, engine="netcdf4")

# Create model dictionary
model_velocity_CF = {
    "model_name": model_name,
    "model_path": model_path,
    "sfincs_model": model_obj,
    "sfincs_results": model_obj.results,
    "sfincs_his": ds_his
    }

model_velocity_CF = compute_hmax_masked_riserate([model_velocity_CF], gwso)

hvmax_CF = model_velocity_CF[0]['sfincs_results']['hvmax_masked']
vmax_CF = model_velocity_CF[0]['sfincs_results']['vmax_masked']
hmax_CF = model_velocity_CF[0]['sfincs_results']['hmax_masked']
rise_CF = model_velocity_CF[0]['sfincs_results']['max_rise_rate_h']

# Apply Jonkman et al. (2008) fatality function based on zones defined by hvmax, vmax, hmax, and rise rate
# Initialize with zeros
fatality_frac_CF = xr.zeros_like(hmax_CF)

# Define masks
CF_mask_breach = (hvmax_CF >= 7) & (vmax_CF > 2)
CF_mask_rapid = (hmax_CF > 2.1) & (rise_CF > 0.5)
CF_mask_remaining = (hmax_CF > 0) & ~CF_mask_breach & ~CF_mask_rapid

# Apply per zone
fatality_frac_CF = xr.where(CF_mask_breach, 1, fatality_frac_CF)
fatality_frac_CF = xr.where(CF_mask_rapid, fatality_fraction_Jonkman_etal2008(hmax_CF, mu_N=1.46, sigma_N=0.28), fatality_frac_CF)
fatality_frac_CF = xr.where(CF_mask_remaining, fatality_fraction_Jonkman_etal2008(hmax_CF, mu_N=7.60, sigma_N=2.75), fatality_frac_CF)

# Store the result
model_velocity_CF[0]['sfincs_results']['fatality_frac'] = fatality_frac_CF

#%%
# Calculate exposed fatalities for exposed population to factual flooding
flood_fatalities_CF_2019 = model_velocity_CF[0]['sfincs_results']['fatality_frac'] * ra_exposed_pop_2019_CF
flood_fatalities_CF_1990 = model_velocity_CF[0]['sfincs_results']['fatality_frac'] * ra_exposed_pop_1990_CF
flood_fatalities_CF_2019_uniform = model_velocity_CF[0]['sfincs_results']['fatality_frac'] * ra_exposed_pop_2019_CF_uniform

#%%
print(f"Total fatalities acc. to Jonkman et al. (2008):")
print(f"Factual climate 2019 population scenario:                {flood_fatalities_F_2019.sum().values:.0f} people")
print(f"Factual climate 2019 UNIFORM population scenario:        {flood_fatalities_F_2019_uniform.sum().values:.0f} people")
print(f"Counterfactual climate 2019 population scenario:         {flood_fatalities_CF_2019.sum().values:.0f} people")
print(f"Counterfactual climate 2019 UNIFORM population scenario: {flood_fatalities_CF_2019_uniform.sum().values:.0f} people")
print(f"Factual climate 1990 population scenario:                {flood_fatalities_F_1990.sum().values:.0f} people")
print(f"Counterfactual climate 1990 population scenario:         {flood_fatalities_CF_1990.sum().values:.0f} people")

#%%
def aggregate_fatalities(total_fatalities_array, flood_raster, transform, crs, region=None, background=None, factor=100):
    print("Aggregating population raster to coarser grid polygons...")

    # Pixel size from transform
    pixel_width = transform.a
    pixel_height = -transform.e

    # Block (coarse) size in meters
    cell_width = factor * pixel_width
    cell_height = factor * pixel_height

    # --- Helper functions ---
    def block_sum(arr, factor):
        nrows, ncols = arr.shape
        nrows_crop = nrows - nrows % factor
        ncols_crop = ncols - ncols % factor
        arr_cropped = arr[:nrows_crop, :ncols_crop]
        return arr_cropped.reshape(nrows_crop//factor, factor, ncols_crop//factor, factor).sum(axis=(1,3))

    def block_mean(arr, factor):
        nrows, ncols = arr.shape
        nrows_crop = nrows - nrows % factor
        ncols_crop = ncols - ncols % factor
        arr_cropped = arr[:nrows_crop, :ncols_crop]
        return np.nanmean(arr_cropped.reshape(nrows_crop//factor, factor, ncols_crop//factor, factor), axis=(1,3))

    # def block_count_threshold(arr, factor, threshold=1):
    #     nrows, ncols = arr.shape
    #     nrows_crop = nrows - nrows % factor
    #     ncols_crop = ncols - ncols % factor
    #     arr_cropped = arr[:nrows_crop, :ncols_crop]
    #     return (arr_cropped.reshape(nrows_crop//factor, factor, ncols_crop//factor, factor) > threshold).sum(axis=(1,3))

    # --- Aggregate population and flood ---
    total_agg = block_sum(total_fatalities_array.values, factor)
    # exposed_agg = block_sum(np.where(flood_raster > 0, total_pop_array, 0), factor)
    avg_flood_depth = block_mean(flood_raster.values, factor)
    # pixels_high = block_count_threshold(flood_raster, factor=factor, threshold=1)

    # --- Build coarse grid ---
    nrows_coarse, ncols_coarse = total_agg.shape
    x0, y0 = transform * (0, 0)
    x_coords = x0 + np.arange(ncols_coarse) * cell_width
    y_coords = y0 + np.arange(nrows_coarse) * -cell_height

    grid_cells = [box(x, y - cell_height, x + cell_width, y) for y in y_coords for x in x_coords]

    gdf_grid = gpd.GeoDataFrame(
        {
            "total_fatalities": total_agg.flatten(),
            # "exposed_population": exposed_agg.flatten(),
            "avg_flood_depth": avg_flood_depth.flatten(),
            # "pct_cells_higher_1m": (pixels_high.flatten() / (factor**2)) * 100
        },
        geometry=grid_cells,
        crs=crs
    )

    # --- Clip to region/background ---
    region_bg = gpd.overlay(region, background, how="intersection") if background is not None else region.copy()

    # --- Redistribute population proportionally to area inside region ---
    gdf_grid = gdf_grid.reset_index(names="cell_id")
    intersections = gpd.overlay(gdf_grid, region_bg, how="intersection")
    intersections["intersect_area"] = intersections.geometry.area
    area_sum = intersections.groupby("cell_id")["intersect_area"].transform("sum")
    intersections["norm_fraction"] = intersections["intersect_area"] / area_sum

    # Redistribute full population across pieces (sum per cell preserved)
    for col in ["total_fatalities"]:
        intersections[col] = intersections.groupby("cell_id")[col].transform("first") * intersections["norm_fraction"]

    # Aggregate pieces back to one polygon per cell
    gdf_grid_masked = intersections.dissolve(by="cell_id", aggfunc="sum")
    gdf_grid_masked["geometry"] = intersections.dissolve(by="cell_id").geometry

    # Relative exposure
    # gdf_grid_masked["relative_fatalities"] = (
    #     gdf_grid_masked["exposed_population"] / gdf_grid_masked["total_population"] * 100
    # )

    print(f"Total fatalities (original): {np.nansum(total_fatalities_array):,.2f}")
    print(f"Total fatalities (aggregated): {gdf_grid_masked['total_fatalities'].sum():,.2f}")
    print(f"Diff: {np.nansum(total_fatalities_array) - gdf_grid_masked['total_fatalities'].sum():,.2f}")
    print(f"Diff %: {((np.nansum(total_fatalities_array) - gdf_grid_masked['total_fatalities'].sum()) / np.nansum(total_fatalities_array)) * 100:,.2f}")
    
    return gdf_grid_masked

#%%
gdf_fatalities_2019_exposed_F_coarse  = aggregate_fatalities(flood_fatalities_F_2019, hmax_F, flood_grid_transform, flood_grid_crs, region_utm, background_utm)


# %%
#######################################################################
################################## PLOTTING ###########################
#######################################################################
# General settings for plotting
minx, miny, maxx, maxy = region_utm.bounds.minx.item(), region_utm.bounds.miny.item(), region_utm.bounds.maxx.item(), region_utm.bounds.maxy.item()
utm_crs = ccrs.UTM(zone=36, southern_hemisphere=True)
subplot_labels_2 = ['(a)', '(b)']
subplot_labels_3 = ['(a)', '(b)', '(c)']
subplot_labels_4 = ['(a)', '(b)', '(c)', '(d)']

#%%
print("Plotting spatially aggregated exposed population change")
from matplotlib.colors import PowerNorm
from matplotlib.cm import ScalarMappable

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5), dpi=300, sharey=True, 
                         constrained_layout=True, subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
im = axes[0].imshow(hmax_F, cmap='viridis', extent=flood_extent, origin='lower', 
                    vmin=0, vmax=3.5, zorder=2)

agg_pop_2019_exposed_F[agg_pop_2019_exposed_F['fatalities'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
norm = PowerNorm(gamma=0.5, vmin=agg_pop_2019_exposed_F['fatalities'].min(), vmax=agg_pop_2019_exposed_F['fatalities'].max())

plot = agg_pop_2019_exposed_F[agg_pop_2019_exposed_F['fatalities'] > 0].plot(column='fatalities', cmap='Greens', edgecolor='grey',
                                    linewidth=0.2, ax=axes[1], legend=False, zorder=2, norm=norm, rasterized=True)
subplot_labels = ['(a)', '(b)']

for i, ax in enumerate(axes):
    # Add model region
    region_utm.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)

    # # Add background and set extent (based on actual lat/lon coordinates)
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    minx, miny, maxx, maxy = region.bounds.minx.item(), region.bounds.miny.item(), region.bounds.maxx.item(), region.bounds.maxy.item()
    ax.set_extent(flood_extent, crs=ccrs.UTM(36, southern_hemisphere=True))

    # Add gridlines and format tick labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.right_labels = False
    gl.top_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}
    if i != 0: 
        gl.left_labels = False
    
    ax.text(0, 1.02, subplot_labels[i], transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='bottom', ha='left')

 # Colorbars
cbar = plt.colorbar(im, ax=axes[0], shrink=0.8)
cbar.set_label("Flood depth (m)")
axes[0].set_title("Factual Flood Depth")

# Continuous scale using actual min/max of your population column
vmin, vmax = agg_pop_2019_exposed_F["fatalities"].min(), agg_pop_2019_exposed_F["fatalities"].max()
sm1 = ScalarMappable(cmap="Greens", norm=norm)
sm1._A = []  # required for colorbar without passing data
cbar = fig.colorbar(sm1, ax=axes[1], orientation="vertical", shrink=0.8)
cbar.set_label("Aggregated fatalities (# people)", fontsize=10)
cbar.ax.tick_params(labelsize=9)

axes[0].set_title("Factual flooding", fontsize=11)
axes[1].set_title("Factual fatalities acc. to Boyd", fontsize=11)


#%%
# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5), dpi=300, sharey=True, 
                         constrained_layout=True, subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
im = axes[0].imshow(hmax_F, cmap='viridis', extent=flood_extent, origin='lower', 
                    vmin=0, vmax=3.5, zorder=2)

gdf_fatalities_2019_exposed_F_coarse[gdf_fatalities_2019_exposed_F_coarse['total_fatalities'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
norm = PowerNorm(gamma=0.5, vmin=gdf_fatalities_2019_exposed_F_coarse['total_fatalities'].min(), vmax=gdf_fatalities_2019_exposed_F_coarse['total_fatalities'].max())

plot = gdf_fatalities_2019_exposed_F_coarse[gdf_fatalities_2019_exposed_F_coarse['total_fatalities'] > 0].plot(column='total_fatalities', cmap='Greens', edgecolor='grey',
                                    linewidth=0.2, ax=axes[1], legend=False, zorder=2, norm=norm, rasterized=True)
subplot_labels = ['(a)', '(b)']

for i, ax in enumerate(axes):
    # Add model region
    region_utm.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)

    # # Add background and set extent (based on actual lat/lon coordinates)
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    minx, miny, maxx, maxy = region.bounds.minx.item(), region.bounds.miny.item(), region.bounds.maxx.item(), region.bounds.maxy.item()
    ax.set_extent(flood_extent, crs=ccrs.UTM(36, southern_hemisphere=True))

    # Add gridlines and format tick labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.right_labels = False
    gl.top_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}
    if i != 0: 
        gl.left_labels = False
    
    ax.text(0, 1.02, subplot_labels[i], transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='bottom', ha='left')

 # Colorbars
cbar = plt.colorbar(im, ax=axes[0], shrink=0.8)
cbar.set_label("Flood depth (m)")
axes[0].set_title("Factual Flood Depth")

# Continuous scale using actual min/max of your population column
vmin, vmax = agg_pop_2019_exposed_F["fatalities"].min(), agg_pop_2019_exposed_F["fatalities"].max()
sm1 = ScalarMappable(cmap="Greens", norm=norm)
sm1._A = []  # required for colorbar without passing data
cbar = fig.colorbar(sm1, ax=axes[1], orientation="vertical", shrink=0.8)
cbar.set_label("Aggregated fatalities (# people)", fontsize=10)
cbar.ax.tick_params(labelsize=9)

axes[0].set_title("Factual flooding", fontsize=11)
axes[1].set_title("Factual fatalities acc. to Jonkman", fontsize=11)


#%%
#############################################################################
############# plot flood variables contributing to fatalities ###############
#############################################################################
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 6), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Plot the flood depth
hmax_F = model_velocity_F[0]['sfincs_results']['hmax_masked'].load()
hmax_CF = model_velocity_CF[0]['sfincs_results']['hmax_masked'].load()
im1 = hmax_CF.plot.pcolormesh(ax=axes[0], cmap="viridis", vmin=0, vmax=3.5, add_colorbar=False, transform=utm_crs, rasterized=True)

# Plot the maximum water velocity
maxvel_F = model_velocity_F[0]['sfincs_results']['vmax_masked'].load()
maxvel_CF = model_velocity_CF[0]['sfincs_results']['vmax_masked'].load()
im2 = maxvel_CF.plot.pcolormesh(ax=axes[1], cmap="viridis", vmin=0, vmax=maxvel_F.quantile(0.99).item(), add_colorbar=False, transform=utm_crs, rasterized=True)

# Plot the maximum rise rate of water depth
maxriserate_F = model_velocity_F[0]['sfincs_results']['max_rise_rate_h'].load()
maxriserate_CF = model_velocity_CF[0]['sfincs_results']['max_rise_rate_h'].load()
im3 = maxriserate_CF.plot.pcolormesh(ax=axes[2], cmap="viridis", vmin=0, vmax=maxriserate_F.quantile(0.99).item(), add_colorbar=False, transform=utm_crs, rasterized=True)

for i, ax in enumerate(axes):
    # Add model region
    region_utm.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)

    # Add background and set extent (based on actual lat/lon coordinates)
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)

    # Add gridlines and format tick labels
    ax.set_xticks([630000,650000, 670000, 690000, 710000])
    ax.set_yticks(np.arange(miny, maxy + 20000, 20000))  # every 1 km
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda y, pos: f"{y/1e6:.2f}"))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, pos: f"{y/1e6:.2f}"))
    ax.tick_params(labelsize=8)
    ax.grid(True, which='major', linestyle='--', color='lightgray', linewidth=0.5)
    ax.set_xlabel("x coordinate UTM zone 36S [×10⁶ m]", size=8)
    ax.set_ylabel("y coordinate UTM zone 36S [×10⁶ m]", size=8)
    ax.set_title("")
    if i != 0:  
        ax.left_labels = False  # disable y-axis labels
        ax.set_ylabel("", size=8)
    ax.text(0, 1.02, subplot_labels_3[i], transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='bottom', ha='left')
    
for ax in axes:
    ax.set_extent([minx, maxx, miny, maxy], crs= ccrs.UTM(36, southern_hemisphere=True))

# ==== Colorbar for flood depth ====
cbar1 = fig.colorbar(im1, ax=axes[0], orientation="vertical", 
                     fraction=0.035, aspect=20, pad=0.02)
cbar1.set_label("Maximum flood depth (m)", labelpad=5, fontsize=9)
cbar1.ax.tick_params(labelsize=8)

# ==== Colorbar for flood velocity ====
cbar2 = fig.colorbar(im2, ax=axes[1], orientation="vertical", 
                     fraction=0.035, aspect=20, pad=0.02)
cbar2.set_label("Maximum velocity (m/s)", labelpad=5, fontsize=9)
cbar2.ax.tick_params(labelsize=8)

# ==== Colorbar for rate of rising ====
cbar2 = fig.colorbar(im3, ax=axes[2], orientation="vertical", 
                     fraction=0.035, aspect=20, pad=0.02)
cbar2.set_label("Maximum rising rate (m/h)", labelpad=5, fontsize=9)
cbar2.ax.tick_params(labelsize=8)

# fig.savefig("../figures/f03.png", bbox_inches='tight', dpi=300)
# fig.savefig("../figures/f03.pdf", bbox_inches='tight', dpi=300)


# %%
#############################################################################
############ Plot Factual and Counterfactual fatalities risk zones ##########
#############################################################################
F_fatal_mask_plot = xr.full_like(hmax_F, np.nan)
CF_fatal_mask_plot = xr.full_like(hmax_CF, np.nan)

# --- FACTUAL ---
# remaining → 1
F_fatal_mask_plot = xr.where(F_mask_remaining, 1, F_fatal_mask_plot)
# rapid → 2
F_fatal_mask_plot = xr.where(F_mask_rapid, 2, F_fatal_mask_plot)
# breach → 3
F_fatal_mask_plot = xr.where(F_mask_breach, 3, F_fatal_mask_plot)
# --- COUNTERFACTUAL ---
# remaining → 1
CF_fatal_mask_plot = xr.where(CF_mask_remaining, 1, CF_fatal_mask_plot)
# rapid → 2
CF_fatal_mask_plot = xr.where(CF_mask_rapid, 2, CF_fatal_mask_plot)
# breach → 3
CF_fatal_mask_plot = xr.where(CF_mask_breach, 3, CF_fatal_mask_plot)

# lightest → darkest
cmap_fatal_risk = mcolors.ListedColormap(['#fee0d2', "#fc9272", "#de2d26"])
bounds = [0.5, 1.5, 2.5, 3.5]  # boundaries between mask categories
norm_fatal_risk = mcolors.BoundaryNorm(bounds, cmap_fatal_risk.N)

#%%
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 6), dpi=300, constrained_layout=True, 
                         subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

im1 = F_fatal_mask_plot.plot.pcolormesh(ax=axes[0], cmap=cmap_fatal_risk, norm=norm_fatal_risk, add_colorbar=False,
                                 rasterized=True)

im2 = CF_fatal_mask_plot.plot.pcolormesh(ax=axes[1], cmap=cmap_fatal_risk, norm=norm_fatal_risk, add_colorbar=False,
                                 rasterized=True)

for i, ax in enumerate(axes):
    # Add model region
    region_utm.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)

    # # Add background and set extent (based on actual lat/lon coordinates)
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)

    # Add gridlines and format tick labels
    ax.set_xticks([630000, 650000, 670000, 690000, 710000])
    ax.set_yticks(np.arange(miny, maxy + 20000, 20000))  # every 1 km
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda y, pos: f"{y/1e6:.2f}"))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, pos: f"{y/1e6:.2f}"))
    ax.tick_params(labelsize=8)
    ax.set_axisbelow(True)
    ax.grid(True, linewidth=0.3)
    ax.set_xlabel("x coordinate UTM zone 36S [×10⁶ m]", size=8)
    ax.set_ylabel("y coordinate UTM zone 36S [×10⁶ m]", size=8)
    ax.set_title("")
    if i != 0:  
        ax.left_labels = False  # disable y-axis labels
        ax.set_ylabel("", size=8)
    ax.text(0, 1.02, subplot_labels_2[i], transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='bottom', ha='left')
    
for ax in axes:
    ax.set_extent([minx, maxx, miny, maxy], crs=ccrs.UTM(36, southern_hemisphere=True))

# Titles
axes[0].set_title("Factual climate", fontsize=10)
axes[1].set_title("Counterfactual climate", fontsize=10)

# Add manual axis for the vertical legend (colorbar)
cax = fig.add_axes([1.01, 0.25, 0.02, 0.5])  
cb = mcolorbar.ColorbarBase(cax, cmap=cmap_fatal_risk, norm=norm_fatal_risk, boundaries=bounds, ticks=[1, 2, 3],
                             orientation="vertical")
cb.ax.set_yticklabels(['Low', 'Middle', 'High'])
cb.ax.tick_params(labelsize=8)
cb.set_label("Fatality risk", fontsize=9, labelpad=10)

for t in cb.ax.get_yticklabels():
    t.set_rotation(90)      # vertical text
    t.set_va("center")      # vertically centered
    t.set_ha("center")      # horizontally centered


# %%
########################################
############ plot difference ###########
########################################
# Plot the difference in fatality zone between Factual and Counterfactual climate scenarios
diff_sign_fatal = F_fatal_mask_plot - CF_fatal_mask_plot

# Only two colors: decreased risk (-1) → green, increased risk (+1) → purple
cmap_diff_fatal_zone = mcolors.ListedColormap(["forestgreen", "purple"])
bounds_diff = [-1.5, 0, 1.5]
ticks = [(-1.5 + 0)/2, (0 + 1.5)/2]  # midpoints: [-0.75, 0.75]
norm_diff_fatal_zone = mcolors.BoundaryNorm(bounds_diff, cmap_diff_fatal_zone.N)

# Mask no-change cells (0) so they are transparent
diff_plot_fatal = diff_sign_fatal.where(diff_sign_fatal != 0)

fig, ax = plt.subplots(figsize=(4, 6), dpi=300,
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Optionally overlay flooded domain lightly
CF_fatal_mask_plot.plot.pcolormesh(ax=ax, cmap=cmap_fatal_risk, norm=norm_fatal_risk, add_colorbar=False, 
                             rasterized=True, alpha=0.6)

# Plot differences
diff_plot_fatal.plot.pcolormesh(ax=ax, cmap=cmap_diff_fatal_zone,norm=norm_diff_fatal_zone, add_colorbar=False, 
                          rasterized=True)

# Extent, grid, labels
ax.set_extent([minx, maxx, miny, maxy], crs=ccrs.UTM(36, southern_hemisphere=True))
ax.set_xticks([630000, 650000, 670000, 690000, 710000])
ax.set_yticks(np.arange(miny, maxy, 20000))
ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda y, pos: f"{y/1e6:.2f}"))
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, pos: f"{y/1e6:.2f}"))
ax.tick_params(labelsize=8)
ax.grid(True, linewidth=0.3)
ax.set_xlabel("x coordinate UTM zone 36S [×10⁶ m]", size=8)
ax.set_ylabel("y coordinate UTM zone 36S [×10⁶ m]", size=8)
ax.set_title("Change in fatality risk (F – CF)", fontsize=10)
ax.text(0, 1.02, '(c)', transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='bottom', ha='left')

# Add background and region boundary
background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
region_utm.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)

# Add vertical colorbar
cax = fig.add_axes([0.95, 0.25, 0.02, 0.5])
cb = mcolorbar.ColorbarBase(cax, cmap=cmap_diff_fatal_zone, norm=norm_diff_fatal_zone, boundaries=bounds_diff, 
                            ticks=ticks, orientation="vertical") 

cb.ax.set_yticklabels(["Decrease", "Increase"])
cb.ax.set_ylabel("Fatality risk change", fontsize=9, labelpad=10)
cb.ax.tick_params(labelsize=8)
for t in cb.ax.get_yticklabels():
    t.set_rotation(90)
    t.set_va("center")
    t.set_ha("center")

# %%
########################################################################################
########## Plot difference in exceeding flood fatality conditions thresholds ###########
########################################################################################
# Plot the difference in exceeding fatality-relevant thresholds between Factual and Counterfactual climate scenarios
# Function to compute difference in exceeding threshold
def compute_diff(F, CF, threshold, geq=False):
    Fv = F.values
    CFv = CF.values
    valid = (~np.isnan(Fv)) & (~np.isnan(CFv))
    out = np.full_like(Fv, np.nan, dtype=float)  # default NaN
    
    if geq:
        mask_F = Fv >= threshold
        mask_CF = CFv >= threshold
    else:
        mask_F = Fv > threshold
        mask_CF = CFv > threshold
    
    # Only consider cells where at least one exceeds threshold
    compute_idx = valid & (mask_F | mask_CF)
    
    out[compute_idx] = mask_F[compute_idx].astype(int) - mask_CF[compute_idx].astype(int)
    
    return xr.DataArray(out, coords=F.coords, dims=F.dims)

hmax_threshold_diff  = compute_diff(hmax_F, hmax_CF, 2.1)
vmax_threshold_diff  = compute_diff(vmax_F, vmax_CF, 2)
hvmax_threshold_diff = compute_diff(hvmax_F, hvmax_CF, 7, geq=True)
rise_threshold_diff  = compute_diff(rise_F, rise_CF, 0.5)

# Define colormap for difference plots
cmap_fl_diff = mcolors.ListedColormap([
    "#66c2a5",            # soft green
    mcolors.to_rgba("#fdae61", alpha=0.5),  # transparent soft orange
    "#ca3329"             # soft red
])
bounds = [-1.5, -0.5, 0.5, 1.5]
norm_fl_diff = mcolors.BoundaryNorm(bounds, cmap_fl_diff.N)


# Setup figure and axes
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(9, 10), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Plot maps
im1 = axes[0, 0].imshow(hmax_threshold_diff, cmap=cmap_fl_diff, norm=norm_fl_diff, extent=[minx, maxx, miny, maxy],
                     transform=ccrs.UTM(36, southern_hemisphere=True), origin='lower', zorder=2)
axes[0, 0].set_title('Flood depth > 2.1 m', size=11)

im2 = axes[0, 1].imshow(vmax_threshold_diff, cmap=cmap_fl_diff, norm=norm_fl_diff, extent=[minx, maxx, miny, maxy],
                     transform=ccrs.UTM(36, southern_hemisphere=True), origin='lower', zorder=2)
axes[0, 1].set_title('Velocity > 2 m/s', size=11)

im3 = axes[1, 0].imshow(hvmax_threshold_diff, cmap=cmap_fl_diff, norm=norm_fl_diff, extent=[minx, maxx, miny, maxy],
                     transform=ccrs.UTM(36, southern_hemisphere=True), origin='lower', zorder=2)
axes[1, 0].set_title('Flood depth * velocity ≥ 7 m2/s', size=11)

im4 = axes[1, 1].imshow(rise_threshold_diff, cmap=cmap_fl_diff, norm=norm_fl_diff, extent=[minx, maxx, miny, maxy],
                     transform=ccrs.UTM(36, southern_hemisphere=True), origin='lower', zorder=2)
axes[1, 1].set_title('Rise rate > 0.5 m/h', size=11)

label_idx = 0
for row in range(axes.shape[0]):
    for col in range(axes.shape[1]):
        ax = axes[row, col]
        # Add background
        for geom in background_utm.geometry:
            ax.add_geometries([geom], crs=ccrs.UTM(36, southern_hemisphere=True),
                              facecolor='#E0E0E0', edgecolor='none', zorder=0)
        # Add region boundary
        for geom in region_utm.geometry:
            ax.add_geometries([geom], crs=ccrs.UTM(36, southern_hemisphere=True),
                              facecolor='none', edgecolor='black', linewidth=0.3, zorder=1)
        # Set extent
        ax.set_extent([minx, maxx, miny, maxy], crs=ccrs.UTM(36, southern_hemisphere=True))

        # Gridlines & ticks
        ax.set_xticks([630000,650000,670000,690000,710000])
        ax.set_yticks(np.arange(miny, maxy, 20000))
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda y, pos: f"{y/1e6:.2f}"))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, pos: f"{y/1e6:.2f}"))
        ax.grid(True, linestyle='--', color='lightgray', linewidth=0.5)

        # Only set x-label for bottom row
        if row == axes.shape[0] - 1:
            ax.set_xlabel("x coordinate UTM zone 36S [×10⁶ m]", size=10)
        # Only set y-label for first column
        if col == 0:
            ax.set_ylabel("y coordinate UTM zone 36S [×10⁶ m]", size=10)
        else:
            ax.set_ylabel("")  # remove y-label
            ax.tick_params(labelleft=False)  # remove y-tick labels
        # Add subplot label
        ax.text(0.02, 0.95, subplot_labels_4[label_idx], transform=ax.transAxes,
                fontsize=12, fontweight='bold', va='top', ha='left')
        label_idx += 1

# Add manual axis for the vertical legend (colorbar)
cax = fig.add_axes([1.01, 0.25, 0.02, 0.5])
cb = mcolorbar.ColorbarBase(cax, cmap=cmap_fl_diff, norm=norm_fl_diff,
                            boundaries=[-1.5, -0.5, 0.5, 1.5], ticks=[-1, 0, 1], orientation='vertical')
cb.ax.set_yticklabels(['Decrease', 'No change', 'Increase'])
cb.ax.tick_params(labelsize=10)
cb.set_label("Changed fatality risk (F-CF)", fontsize=11, labelpad=10)

for t in cb.ax.get_yticklabels():
    t.set_rotation(90)      # vertical text
    t.set_va("center")      # vertically centered
    t.set_ha("center")      # horizontally centered

plt.show()


# %%
####################################################################################
###### Spatially plot factual and counterfactual fatalities for both methods #######
####################################################################################
diff_flood_fatalities_clim    = flood_fatalities_F_2019 - flood_fatalities_CF_2019
diff_flood_fatalities_uniform = flood_fatalities_F_2019 - flood_fatalities_F_2019_uniform
diff_flood_fatalities_pop     = flood_fatalities_F_2019 - flood_fatalities_F_1990 

# Setup figure and axes
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(9, 10), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Plot maps
im1 = axes[0, 0].imshow(flood_fatalities_F_2019, cmap=cmap_fl_diff, norm=norm_fl_diff, extent=[minx, maxx, miny, maxy],
                     transform=ccrs.UTM(36, southern_hemisphere=True), origin='lower', zorder=2)
axes[0, 0].set_title('Factual fatalities Jonkman et al.', size=11)

im2 = axes[0, 1].imshow(flood_fatalities_CF_2019, cmap=cmap_fl_diff, norm=norm_fl_diff, extent=[minx, maxx, miny, maxy],
                     transform=ccrs.UTM(36, southern_hemisphere=True), origin='lower', zorder=2)
axes[0, 1].set_title('Factual fatalities ', size=11)

im3 = axes[1, 0].imshow(hmax_threshold_diff, cmap=cmap_fl_diff, norm=norm_fl_diff, extent=[minx, maxx, miny, maxy],
                     transform=ccrs.UTM(36, southern_hemisphere=True), origin='lower', zorder=2)
axes[1, 0].set_title('Flood depth * velocity ≥ 7 m2/s', size=11)

im4 = axes[1, 1].imshow(rise_threshold_diff, cmap=cmap_fl_diff, norm=norm_fl_diff, extent=[minx, maxx, miny, maxy],
                     transform=ccrs.UTM(36, southern_hemisphere=True), origin='lower', zorder=2)
axes[1, 1].set_title('Rise rate > 0.5 m/h', size=11)

for row in range(axes.shape[0]):
    for col in range(axes.shape[1]):
        ax = axes[row, col]
        # Add background
        for geom in background_utm.geometry:
            ax.add_geometries([geom], crs=ccrs.UTM(36, southern_hemisphere=True),
                              facecolor='#E0E0E0', edgecolor='none', zorder=0)
        # Add region boundary
        for geom in region_utm.geometry:
            ax.add_geometries([geom], crs=ccrs.UTM(36, southern_hemisphere=True),
                              facecolor='none', edgecolor='black', linewidth=0.3, zorder=1)
        # Set extent
        ax.set_extent([minx, maxx, miny, maxy], crs=ccrs.UTM(36, southern_hemisphere=True))

        # Gridlines & ticks
        ax.set_xticks([630000,650000,670000,690000,710000])
        ax.set_yticks(np.arange(miny, maxy, 20000))
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda y, pos: f"{y/1e6:.2f}"))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, pos: f"{y/1e6:.2f}"))
        ax.grid(True, linestyle='--', color='lightgray', linewidth=0.5)

        # Only set x-label for bottom row
        if row == axes.shape[0] - 1:
            ax.set_xlabel("x coordinate UTM zone 36S [×10⁶ m]", size=10)
        # Only set y-label for first column
        if col == 0:
            ax.set_ylabel("y coordinate UTM zone 36S [×10⁶ m]", size=10)
        else:
            ax.set_ylabel("")  # remove y-label
            ax.tick_params(labelleft=False)  # remove y-tick labels

# Add manual axis for the vertical legend (colorbar)
cax = fig.add_axes([1.01, 0.25, 0.02, 0.5])
cb = mcolorbar.ColorbarBase(cax, cmap=cmap_fl_diff, norm=norm_fl_diff,
                            boundaries=[-1.5, -0.5, 0.5, 1.5], ticks=[-1, 0, 1], orientation='vertical')
cb.ax.set_yticklabels(['Decrease', 'No change', 'Increase'])
cb.ax.tick_params(labelsize=10)
cb.set_label("Changed fatality risk (F-CF)", fontsize=11, labelpad=10)

for t in cb.ax.get_yticklabels():
    t.set_rotation(90)      # vertical text
    t.set_va("center")      # vertically centered
    t.set_ha("center")      # horizontally centered

plt.show()
# %%

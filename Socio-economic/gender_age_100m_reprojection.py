#%%
print("Loading packages...")
import os
import yaml
import json
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from os.path import join
import rasterio
from rasterio import features
import geopandas as gpd
import itertools
import warnings
warnings.filterwarnings('ignore')
import platform
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import rioxarray as rxr 
from hydromt import DataCatalog
from tqdm import tqdm
import xarray as xr
import re
import glob


prefix = "p:/" if platform.system() == "Windows" else "/p/"

#%%
# ===== CONFIGURATION =====
print("Setting up paths and parameters...")
# EVENT_NAME = "Idai"
BASE_RUN_PATH = Path(os.path.join(prefix,"11210471-001-compass","03_Runs","sofala","Idai"))
# SCENARIO_PATH_F = "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" # factual
# SCENARIO_PATH_CF = "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10" # counterfactual

# ===== FILE PATHS =====
# Base directory for the specific event and scenario
# fiat_dir_F    = BASE_RUN_PATH / "fiat" / SCENARIO_PATH_F
# fiat_dir_CF   = BASE_RUN_PATH / "fiat" / SCENARIO_PATH_CF
# sfincs_dir_F  = BASE_RUN_PATH / "sfincs" / SCENARIO_PATH_F
# sfincs_dir_CF = BASE_RUN_PATH / "sfincs" / SCENARIO_PATH_CF

# ===== DATA CATALOG =====
# if platform.system() == "Windows":
#     datacat_path = os.path.abspath("../Workflows/03_data_catalogs/datacatalog_general.yml")
# else:
#     datacat_path = os.path.abspath("../Workflows/03_data_catalogs/datacatalog_general___linux.yml")
# data_catalog = DataCatalog(data_libs = [datacat_path])

#%%
# Input files
shapefile_fp = BASE_RUN_PATH / "sfincs" / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" / "gis" / "region.geojson"   # replace with your region shapefile
background = gpd.read_file(os.path.join(prefix, "11210471-001-compass","01_Data","sofala_geoms","sofala_region_background.geojson"))
# shapefile_sofala = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_province.shp")

# Load the admin3 district in the case study region to validate exposed people
# districts = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_districts_study_region.shp")

# # population in provided inthousand persons per grid cell
# population_raster_path_2020 = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_2020_30.tif")  
# population_raster_path_1990 = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_1990_30.tif")  

# # flood raster
# F_flooding = sfincs_dir_F / "floodmap.tif"
# CF_flooding = sfincs_dir_CF / "floodmap.tif"

# Flood model subgrid
sfincs_subgrid = BASE_RUN_PATH / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" / "sfincs" / "subgrid" / "dep_subgrid.tif"


#%% Read flood data and background polygons
# --- Flood grid properties ---
with rasterio.open(sfincs_subgrid) as src:
    flood_grid_crs, flood_grid_transform, flood_grid_shape = src.crs, src.transform, (src.height, src.width)

# --- Setup region ---
region = gpd.read_file(shapefile_fp).to_crs(flood_grid_crs)
region_geom = [json.loads(region.to_json())["features"][0]["geometry"]]

# --- Read flood rasters ---
# hmax_F_da = rxr.open_rasterio(F_flooding).squeeze("band", drop=True)  # if single-band
# hmax_CF_da = rxr.open_rasterio(CF_flooding).squeeze("band", drop=True)  # if single-band
# hmax_F = hmax_F_da.values
# hmax_CF = hmax_CF_da.values

# # get extent from raster transform
# def get_extent(transform, width, height):
#     left = transform[2]
#     right = left + width * transform[0]
#     top = transform[5]
#     bottom = top + height * transform[4]
#     return [left, right, top, bottom]

# flood_extent = get_extent(flood_grid_transform, flood_grid_shape[1], flood_grid_shape[0])

# # --- Background layers for plotting ---
# # Define a polygon to remove/mask out land layer (incorrect boundary over the ocean)
# mask_poly = Polygon([(34.9,-20.3), (36,-20.3), (36,-19.9), (34.9,-19.9)])
# bg_filtered = background.copy()
# bg_filtered['geometry'] = bg_filtered.geometry.apply(lambda g: g.difference(mask_poly))

# Reproject background and region to flood grid CRS for consistent plotting
background_utm = background.to_crs(flood_grid_crs)
# bg_filtered_utm = bg_filtered.to_crs(flood_grid_crs)
region_utm = region.to_crs(flood_grid_crs)
# districts_utm = districts.to_crs(flood_grid_crs)

# # Remove districts that are not connecting to the region
# drop_districts = ["Muanza", "Gororngosa-Sede", "Galinha"]
# districts_filtered = districts_utm[~districts_utm['NAME_3'].isin(drop_districts)]

#%%
####################################################################
############### Load World Pop Age & Sex Structures ################
####################################################################
# WorldPop data source: https://hub.worldpop.org/project/categories?id=8
worldpop_folder_2020_100m = os.path.join(prefix,"11210471-001-compass","01_Data","population_data", "Worldpop","moz_agesex_structures_2020_CN_100m_R2025A_v1")
path_moz_agesex_combined_2020_100m = os.path.join(worldpop_folder_2020_100m, "moz_agesex_structures_2020_CN_100m_R2025A_v1.nc")

#%% ###################################################################################
######################### Overlay population with flood map ###########################
#######################################################################################
print("Reprojecting WorldPop 100m age-sex rasters to flood grid...")
def prepare_worldpop_to_flood_grid(
    folder_raster,
    path_output_nc,
    land_gdf,
    flood_crs,
    flood_transform,
    flood_shape,
    region_gdf=None,
    exclude_genders=None,
    year=None,
    chunk_size={'x': 1000, 'y': 1000},
    province_geom=None,
    region_name=None,
    districts=None
):
    """
    Combine, clip (by shapefile), reproject, and redistribute WorldPop 100m age-sex rasters to flood grid.

    Args:
        folder_raster (str): Folder with WorldPop rasters (e.g., moz_f_00.tif).
        path_output_nc (str): Output NetCDF for combined flood-grid dataset.
        land_gdf (GeoDataFrame): Land geometry for redistributing only over land.
        flood_crs, flood_transform, flood_shape: Flood grid definition.
        region_gdf (GeoDataFrame, optional): Shapefile of the region to clip to.
        exclude_genders (list, optional): e.g., ['t'] to skip total gender.
        year (int, optional): Dataset year.
        chunk_size (dict): Chunk size for lazy raster reading.
        province_geom, region_name, districts: Optional redistribution metadata.
    """
    from socio_eco_funcs import reproject_and_redistribute_population_over_land  # <-- your existing helper

    files = sorted(glob.glob(os.path.join(folder_raster, "moz_*.tif")))
    if exclude_genders is None:
        exclude_genders = []

    if region_gdf is not None:
        # Ensure CRS matches raster CRS
        with rxr.open_rasterio(files[0]) as sample_da:
            raster_crs = sample_da.rio.crs
        if region_gdf.crs != raster_crs:
            region_gdf = region_gdf.to_crs(raster_crs)
        print(f"✅ Clipping rasters using shapefile ({len(region_gdf)} feature(s))")

    print(f"Processing {len(files)} WorldPop rasters...")

    reprojected_dict = {}
    genders, ages = [], []

    for f in tqdm(files, desc="Reprojecting to flood grid"):
        match = re.search(r"_([fmt])_(\d+)", f)
        if not match:
            print(f"Skipping unrecognized filename: {os.path.basename(f)}")
            continue

        gender, age = match.group(1), int(match.group(2))
        if gender in exclude_genders:
            continue

        da = rxr.open_rasterio(f, chunks=chunk_size).squeeze()
        da = da.where(da != -99999, 0)

        # --- 🗺️ Clip raster using region shapefile if provided ---
        if region_gdf is not None:
            da = da.rio.clip(region_gdf.geometry, region_gdf.crs, drop=True)
        
        # Write temporary clipped raster for redistribution
        tmp_clip_path = Path(folder_raster) / f"_tmp_clipped_{gender}_{age:02d}.tif"
        da.rio.to_raster(tmp_clip_path)

        # Redistribute population to flood grid
        pop_fine, _, _, _, _ = reproject_and_redistribute_population_over_land(
            pop_path=tmp_clip_path,
            land_gdf=land_gdf,
            flood_crs=flood_crs,
            flood_transform=flood_transform,
            flood_shape=flood_shape,
            # province_geom=province_geom,
            region=region_name,
            # districts=districts,
            year=year
        )

        reprojected_dict[(gender, age)] = pop_fine
        genders.append(gender)
        ages.append(age)
        tmp_clip_path.unlink(missing_ok=True)

    genders = sorted(set(genders))
    ages = sorted(set(ages))

    # Combine all into 4D array
    arr = np.zeros((len(genders), len(ages), flood_shape[0], flood_shape[1]), dtype=np.float32)
    for (g, a), data in reprojected_dict.items():
        gi, ai = genders.index(g), ages.index(a)
        arr[gi, ai, :, :] = data

    # Create spatial coordinates
    x_coords = np.arange(flood_shape[1]) * flood_transform[0] + flood_transform[2] + flood_transform[0] / 2
    y_coords = np.arange(flood_shape[0]) * flood_transform[4] + flood_transform[5] + flood_transform[4] / 2

    ds = xr.Dataset(
        {'population': (['genders', 'ages', 'y', 'x'], arr)},
        coords={'genders': genders, 'ages': ages, 'y': y_coords, 'x': x_coords}
    )
    ds = ds.rio.write_crs(flood_crs)
    ds = ds.chunk(chunk_size)

    ds.attrs.update({
        'title': f'WorldPop Age-Sex Data reprojected to flood grid ({year})',
        'source': 'WorldPop',
        'resolution': '100m',
        'country': 'Mozambique',
        'year': year,
        'chunked': True,
        'chunk_size': chunk_size,
        'clipped_to_region': region_gdf is not None,
        'reprojected_to_flood_grid': True,
        'created_date': pd.Timestamp.now().isoformat(),
    })

    ds.to_netcdf(path_output_nc, encoding={
        'population': {
            'zlib': True,
            'complevel': 6,
            'chunksizes': (len(genders), len(ages),
                           min(1000, flood_shape[0]), min(1000, flood_shape[1]))
        }
    })
    print(f"\n✅ Saved reprojected dataset: {path_output_nc}")

    return ds

# %%

ds_reprojected = prepare_worldpop_to_flood_grid(
    folder_raster=worldpop_folder_2020_100m,
    path_output_nc=path_moz_agesex_combined_2020_100m,
    land_gdf=background,
    flood_crs=flood_grid_crs,
    flood_transform=flood_grid_transform,
    flood_shape=flood_grid_shape,
    region_gdf=region,
    exclude_genders=['t'],   # skip total gender if not needed
    year=2020,
    chunk_size={'x': 1000, 'y': 1000}
)

# %%

#%%
import os
from matplotlib import cm
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
import matplotlib.colors as mcolors
from matplotlib.colors import PowerNorm
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import Normalize
from rasterio.vrt import WarpedVRT
from rasterio.warp import reproject, Resampling
from rasterio.mask import mask
from shapely.geometry import box
from shapely.geometry import Polygon
import cartopy.crs as ccrs
import rioxarray as rxr 
from hydromt_sfincs import SfincsModel, utils
from hydromt import DataCatalog
from tqdm import tqdm
import gc
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks


prefix = "p:/" if platform.system() == "Windows" else "/p/"

# ===== CONFIGURATION =====
EVENT_NAME = "Idai"
BASE_RUN_PATH = Path("C:/Code/Paper_1/Data_submission")
SCENARIO_PATH_F = "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" # factual
SCENARIO_PATH_CF = "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10" # counterfactual

# ===== FILE PATHS =====
# Base directory for the specific event and scenario
fiat_dir_F    = BASE_RUN_PATH / "fiat" / SCENARIO_PATH_F
fiat_dir_CF   = BASE_RUN_PATH / "fiat" / SCENARIO_PATH_CF
sfincs_dir_F  = BASE_RUN_PATH / "sfincs" / SCENARIO_PATH_F
sfincs_dir_CF = BASE_RUN_PATH / "sfincs" / SCENARIO_PATH_CF

# ===== DATA CATALOG =====
if platform.system() == "Windows":
    datacat_path = os.path.abspath("../Workflows/03_data_catalogs/datacatalog_general.yml")
else:
    datacat_path = os.path.abspath("../Workflows/03_data_catalogs/datacatalog_general___linux.yml")
data_catalog = DataCatalog(data_libs = [datacat_path])

#%%
# ===== Input files ==== #
shapefile_fp = "p:/11210471-001-compass/03_Runs/sofala/Idai/sfincs/event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0/gis/region.geojson"   # replace with your region shapefile
background = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_region_background.geojson")
shapefile_sofala = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_province.shp")

# Load the admin3 district in the case study region to validate exposed people
districts = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_districts_study_region.shp")

# population in provided inthousand persons per grid cell
population_raster_path_2020 = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_2020_30.tif")  
population_raster_path_1990 = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_1990_30.tif")  

# flood raster
F_flooding = sfincs_dir_F / "floodmap.tif"
CF_flooding = sfincs_dir_CF / "floodmap.tif"

# Flood model subgrid
sfincs_subgrid = BASE_RUN_PATH / "sfincs" / "subgrid" / "dep_subgrid.tif"


#%% Read flood data and background polygons
# --- Flood grid properties ---
with rasterio.open(sfincs_subgrid) as src:
    flood_grid_crs, flood_grid_transform, flood_grid_shape = src.crs, src.transform, (src.height, src.width)

# --- Setup region ---
region = gpd.read_file(shapefile_fp).to_crs(flood_grid_crs)
region_geom = [json.loads(region.to_json())["features"][0]["geometry"]]

# --- Read flood rasters ---
hmax_F_da = rxr.open_rasterio(F_flooding).squeeze("band", drop=True)  # if single-band
hmax_CF_da = rxr.open_rasterio(CF_flooding).squeeze("band", drop=True)  # if single-band
hmax_F = hmax_F_da.values
hmax_CF = hmax_CF_da.values

# get extent from raster transform
def get_extent(transform, width, height):
    left = transform[2]
    right = left + width * transform[0]
    top = transform[5]
    bottom = top + height * transform[4]
    return [left, right, top, bottom]

flood_extent = get_extent(flood_grid_transform, flood_grid_shape[1], flood_grid_shape[0])

# --- Background layers for plotting ---
# Define a polygon to remove/mask out land layer (incorrect boundary over the ocean)
mask_poly = Polygon([(34.9,-20.3), (36,-20.3), (36,-19.9), (34.9,-19.9)])
bg_filtered = background.copy()
bg_filtered['geometry'] = bg_filtered.geometry.apply(lambda g: g.difference(mask_poly))

# Reproject background and region to flood grid CRS for consistent plotting
background_utm = background.to_crs(flood_grid_crs)
bg_filtered_utm = bg_filtered.to_crs(flood_grid_crs)
region_utm = region.to_crs(flood_grid_crs)
districts_utm = districts.to_crs(flood_grid_crs)

# Remove districts that are not connecting to the region
drop_districts = ["Muanza", "Gororngosa-Sede", "Galinha"]
districts_filtered = districts_utm[~districts_utm['NAME_3'].isin(drop_districts)]


#%% Read and regrid population data
# --- Function to redistribute population over land pixels on flood grid ---
def reproject_and_redistribute_population_over_land(pop_path, land_gdf, flood_crs, flood_transform, flood_shape, province_geom=None, region=None, districts=None, year=None, out_raster_path=None):    
    print(f"▶ Loading {year} population raster...")
    with rasterio.open(pop_path) as src:
        pop = src.read(1, masked=True)
        pop_affine = src.transform
        pop_crs = src.crs

        # Clip to province
        province_geom = province_geom.to_crs(src.crs)
        province_geom = [json.loads(province_geom.to_json())["features"][0]["geometry"]]

        pop_sofala, transform_sofala = mask(src, province_geom, crop=True, nodata=src.nodata)

        # Dissolve districts into one geometry
        districts_single = districts.dissolve().reset_index(drop=True)
        districts_single = districts_single.to_crs(src.crs)
        districts_geom = [districts_single.geometry.iloc[0].__geo_interface__]

        pop_districts, pop_affine_districts = mask(src, districts_geom, crop=True, nodata=src.nodata)

        if out_raster_path is not None and os.path.exists(out_raster_path):
            print(f"▶ Loading existing raster from {out_raster_path}")
            with rasterio.open(out_raster_path) as src:
                pop_fine = src.read(1)
                return pop_fine, pop_sofala, transform_sofala, pop_districts, pop_affine_districts

        # Clip to region if provided
        if region is not None:
            region_wsg = region.to_crs(src.crs)
            region_geom = [json.loads(region_wsg.to_json())["features"][0]["geometry"]] 
            pop, pop_affine = rasterio.mask.mask(src, region_geom, crop=True, nodata=src.nodata)

        pop = pop.squeeze()

    # Prepare empty high-resolution array
    pop_fine = np.zeros(flood_shape, dtype=np.float32)

    # Rasterize land mask to flood grid (True for land)
    land_mask = features.rasterize(
        [(geom, 1) for geom in land_gdf.geometry],
        out_shape=flood_shape,
        transform=flood_transform,
        fill=0,
        dtype=np.uint8
    ).astype(bool)

    print("  ✔ Land mask created on flood grid.")

    # Loop through each coarse pixel
    print("▶ Redistributing population to fine grid...")
    for row in tqdm(range(pop.shape[0]), desc="  Processing coarse cells"):
        for col in range(pop.shape[1]):
            pop_value = pop[row, col]
            if np.isnan(pop_value) or pop_value <= 0:
                continue

            # Get coarse pixel bounds (in coarse CRS)
            x_min, y_max = pop_affine * (col, row)
            x_max, y_min = pop_affine * (col + 1, row + 1)
            coarse_bounds = box(x_min, y_min, x_max, y_max)

            # Transform to flood CRS
            coarse_bounds_flood = gpd.GeoSeries([coarse_bounds], crs=pop_crs).to_crs(flood_crs).iloc[0]

            # Rasterize this coarse cell footprint to flood grid
            coarse_mask = features.rasterize(
                [(coarse_bounds_flood, 1)],
                out_shape=flood_shape,
                transform=flood_transform,
                fill=0,
                # all_touched=True,
                dtype=np.uint8
            ).astype(bool)

            # Identify land pixels within this coarse cell
            valid_mask = coarse_mask & land_mask
            n_valid = valid_mask.sum()

            # Distribute population evenly over valid pixels
            if n_valid > 0:
                pop_fine[valid_mask] += pop_value / n_valid

    # else: all water → skip or optionally add to nearest land (not done here)
    total_input_pop = float(np.nansum(pop))
    total_output_pop = float(pop_fine.sum())
    diff = abs(total_output_pop - total_input_pop)
    rel_diff = diff / total_input_pop * 100

    # --- 7️⃣ Validation printout ---
    print("  ✔ Redistribution done.")
    print(f"  🔹 Input population:  {total_input_pop:,.0f}")
    print(f"  🔹 Output population: {total_output_pop:,.0f}")
    print(f"  🔹 Difference:        {diff:,.2f} ({rel_diff:.4f} %)")

    if rel_diff > 0.01:
        print("  ⚠ WARNING: Population not perfectly preserved — check CRS or mask alignment!")

    print(f"  ✔ Redistribution done. Total population preserved: {pop_fine.sum():,.0f}")
    
    # Optional: save the result as a GeoTIFF
    if out_raster_path is not None:
        print(f"▶ Saving redistributed population raster to {out_raster_path}")
        new_profile = {
            "driver": "GTiff",
            "dtype": rasterio.float32,
            "count": 1,
            "height": pop_fine.shape[0],
            "width": pop_fine.shape[1],
            "crs": flood_crs,
            "transform": flood_transform,
            "compress": "deflate"
        }
        with rasterio.open(out_raster_path, "w", **new_profile) as dst:
            dst.write(pop_fine, 1)

    return pop_fine, pop_sofala, transform_sofala, pop_districts, pop_affine_districts

# --- Function to link population raster to flood depth as DataFrame ---
def pop_raster_to_df(pop_array, flood_array, transform):
    print("Linking population raster to flood depth as DataFrame...")
    # Flatten arrays
    pop_flat = pop_array.ravel()
    flood_flat = flood_array.ravel()

    # Mask zero-pop cells
    mask = pop_flat > 0
    pop_vals = pop_flat[mask]
    flood_vals = flood_flat[mask]

    # Pixel coordinates (centers)
    rows, cols = np.indices(pop_array.shape)
    xs, ys = transform * (cols.ravel()[mask] + 0.5, rows.ravel()[mask] + 0.5)

    df = pd.DataFrame({
        "population": pop_vals,
        "flood_depth": flood_vals,
        "x": xs,
        "y": ys
    })

    return df

# --- Function to aggregate population and flood depth to coarser grid ---
# Problem is that population is either lost when area weighted averaging or overcounted when summing and redistributing.
def aggregate_pop(total_pop_array, flood_raster, transform, crs, region=None, background=None, factor=100):
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

    def block_count_threshold(arr, factor, threshold=1):
        nrows, ncols = arr.shape
        nrows_crop = nrows - nrows % factor
        ncols_crop = ncols - ncols % factor
        arr_cropped = arr[:nrows_crop, :ncols_crop]
        return (arr_cropped.reshape(nrows_crop//factor, factor, ncols_crop//factor, factor) > threshold).sum(axis=(1,3))

    # --- Aggregate population and flood ---
    total_agg = block_sum(total_pop_array, factor)
    exposed_agg = block_sum(np.where(flood_raster > 0, total_pop_array, 0), factor)
    avg_flood_depth = block_mean(flood_raster, factor)
    pixels_high = block_count_threshold(flood_raster, factor=factor, threshold=1)

    # --- Build coarse grid ---
    nrows_coarse, ncols_coarse = total_agg.shape
    x0, y0 = transform * (0, 0)
    x_coords = x0 + np.arange(ncols_coarse) * cell_width
    y_coords = y0 + np.arange(nrows_coarse) * -cell_height

    grid_cells = [box(x, y - cell_height, x + cell_width, y) for y in y_coords for x in x_coords]

    gdf_grid = gpd.GeoDataFrame(
        {
            "total_population": total_agg.flatten(),
            "exposed_population": exposed_agg.flatten(),
            "avg_flood_depth": avg_flood_depth.flatten(),
            "pct_cells_higher_1m": (pixels_high.flatten() / (factor**2)) * 100
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
    for col in ["total_population", "exposed_population"]:
        intersections[col] = intersections.groupby("cell_id")[col].transform("first") * intersections["norm_fraction"]

    # Aggregate pieces back to one polygon per cell
    gdf_grid_masked = intersections.dissolve(by="cell_id", aggfunc="sum")
    gdf_grid_masked["geometry"] = intersections.dissolve(by="cell_id").geometry

    # Relative exposure
    gdf_grid_masked["relative_population"] = (
        gdf_grid_masked["exposed_population"] / gdf_grid_masked["total_population"] * 100
    )

    print(f"Total pop (original): {np.nansum(total_pop_array):,.2f}")
    print(f"Total pop (aggregated): {gdf_grid_masked['total_population'].sum():,.2f}")
    print(f"Diff: {np.nansum(total_pop_array) - gdf_grid_masked['total_population'].sum():,.2f}")
    print(f"Diff %: {((np.nansum(total_pop_array) - gdf_grid_masked['total_population'].sum()) / np.nansum(total_pop_array)) * 100:,.2f}")

    return gdf_grid_masked



# %%
# ============================================================================================= #
# # More accurate aggregation function but very slow
# def aggregate_pop_sjoin_slow(total_pop_array, flood_raster, transform, crs, region_utm, background_utm, cell_size=2500):
#     # 1️⃣ Clip your region to ensure bounds match data extent
#     clipped_region = gpd.overlay(region_utm, background_utm, how="intersection")

#     # 2️⃣ Define cell size in meters (since UTM)
#     cell_size = 2500  # 2 km × 2 km grid cells

#     minx, miny, maxx, maxy = clipped_region.total_bounds
#     cols = np.arange(minx, maxx, cell_size)
#     rows = np.arange(miny, maxy, cell_size)

#     # 3️⃣ Build grid cells
#     grid_cells = [box(x, y, x + cell_size, y + cell_size) for x in cols for y in rows]
#     gdf_grid = gpd.GeoDataFrame(geometry=grid_cells, crs=clipped_region.crs)

#     # 4️⃣ Clip grid to region (only keep intersecting cells)
#     gdf_grid_masked = gpd.overlay(gdf_grid, clipped_region, how="intersection")

#     # 5️⃣ Convert your fine-resolution raster data to points
#     # Example: convert flood/pop arrays to GeoDataFrame of pixel centers
#     rows, cols = np.where(~np.isnan(total_pop_array))
#     xs, ys = rasterio.transform.xy(flood_grid_transform, rows, cols)
#     gdf_points = gpd.GeoDataFrame(
#         {"pop": total_pop_array[rows, cols],
#         "flood": hmax_F[rows, cols]},
#         geometry=gpd.points_from_xy(xs, ys),
#         crs=flood_grid_crs
#     )

#     # 6️⃣ Spatial join: assign each point to its coarse grid cell
#     joined = gpd.sjoin(gdf_points, gdf_grid_masked, how="left", predicate="within")

#     # 7️⃣ Aggregate statistics per grid cell
#     agg_pop_total = joined.groupby("index_right")["pop"].sum()
#     agg_pop_exposed = joined.loc[joined["flood"] > 0].groupby("index_right")["pop"].sum()
#     agg_flood_mean = joined.groupby("index_right")["flood"].mean()
#     agg_pct_above_1m = (
#         joined.loc[joined["flood"] > 1].groupby("index_right").size() /
#         joined.groupby("index_right").size()
#     ) * 100

#     # 8️⃣ Assign to grid GeoDataFrame
#     gdf_grid_masked["total_population"] = agg_pop_total
#     gdf_grid_masked["exposed_population"] = agg_pop_exposed
#     gdf_grid_masked["avg_flood_depth"] = agg_flood_mean
#     gdf_grid_masked["pct_cells_higher_1m"] = agg_pct_above_1m
#     gdf_grid_masked = gdf_grid_masked.fillna(0)

#     # 9️⃣ Relative exposure
#     gdf_grid_masked["relative_population"] = (
#         gdf_grid_masked["exposed_population"] / gdf_grid_masked["total_population"] * 100
#     ).replace(np.inf, 0)

#     print("✔ Aggregation finished")
#     print(f"Total pop: {np.nansum(total_pop_array):,.0f}")
#     print(f"Total pop: {gdf_grid_masked['total_population'].sum():,.0f}")

#     return gdf_grid_masked



#%%
# ============================================================================================ #
# ====================== Process population directly into flood grid ========================= #
# ============================================================================================ #

pop_arrays = {}
pop_sofala_arrays = {}
pop_sofala_districts = {}
pop_affine_sofala_districts = {}
for year, path in [(1990, population_raster_path_1990),
                   (2020, population_raster_path_2020)]:
    pop_arrays[year], pop_sofala_arrays[year], transform_sofala_land, pop_sofala_districts[year], pop_affine_sofala_districts[year] = reproject_and_redistribute_population_over_land(
        pop_path=path, land_gdf=background_utm, flood_crs=flood_grid_crs, flood_transform=flood_grid_transform,
        flood_shape=flood_grid_shape, province_geom=shapefile_sofala, region=region, districts=districts_filtered, year=year,
        out_raster_path=f"c:/Code/COMPASS_exposure/Data/Modified/population_{year}_region_regrid.tif"
    )

# --- Compute exposed population GeoDataFrames ---
df_pop_2020_flood_depth_F  = pop_raster_to_df(pop_arrays[2020], hmax_F, flood_grid_transform)
df_pop_2020_flood_depth_CF = pop_raster_to_df(pop_arrays[2020], hmax_CF, flood_grid_transform)
df_pop_1990_flood_depth_F  = pop_raster_to_df(pop_arrays[1990], hmax_F, flood_grid_transform)
df_pop_1990_flood_depth_CF = pop_raster_to_df(pop_arrays[1990], hmax_CF, flood_grid_transform)

# simple rasters for fast plotting
ra_exposed_pop_2020_F  = np.where(hmax_F > 0, pop_arrays[2020], 0)
ra_exposed_pop_2020_CF = np.where(hmax_CF > 0, pop_arrays[2020], 0)
ra_exposed_pop_1990_F  = np.where(hmax_F > 0, pop_arrays[1990], 0)
ra_exposed_pop_1990_CF = np.where(hmax_CF > 0, pop_arrays[1990], 0)

# Compute exposed pop on fine grid
df_pop_2020_exposed_F  = df_pop_2020_flood_depth_F[df_pop_2020_flood_depth_F['flood_depth'] > 0]
df_pop_2020_exposed_CF = df_pop_2020_flood_depth_CF[df_pop_2020_flood_depth_CF['flood_depth'] > 0]
df_pop_1990_exposed_F  = df_pop_1990_flood_depth_F[df_pop_1990_flood_depth_F['flood_depth'] > 0]
df_pop_1990_exposed_CF = df_pop_1990_flood_depth_CF[df_pop_1990_flood_depth_CF['flood_depth'] > 0]

# --- Aggregate population to coarser grid ---
gdf_pop_2020_exposed_F_coarse  = aggregate_pop(pop_arrays[2020], hmax_F, flood_grid_transform, flood_grid_crs, region_utm, background_utm)
gdf_pop_2020_exposed_CF_coarse = aggregate_pop(pop_arrays[2020], hmax_CF, flood_grid_transform, flood_grid_crs, region_utm, background_utm)
gdf_pop_1990_exposed_F_coarse  = aggregate_pop(pop_arrays[1990], hmax_F, flood_grid_transform, flood_grid_crs, region_utm, background_utm)
gdf_pop_1990_exposed_CF_coarse = aggregate_pop(pop_arrays[1990], hmax_CF, flood_grid_transform, flood_grid_crs, region_utm, background_utm)

#%%
# Uniform population growth
pop_growth = np.nansum(pop_arrays[2020]) / np.nansum(pop_arrays[1990])
print(f"Uniform population growth from 1990 to 2020 in Sofala region: {(pop_growth*100):.2f}%")

pop_array_uniform_2020 = pop_arrays[1990] * (pop_growth)

# Sanity check
print(f"{np.nansum(pop_array_uniform_2020):,.0f} people in 2020 with uniform growth")
print(f"{np.nansum(pop_arrays[2020]):,.0f} people in 2020 actual")

# Calculate exposure
ra_exposed_pop_2020_F_uniform = np.where(hmax_F > 0, pop_array_uniform_2020, 0)

# Aggregate to coarser cells
gdf_pop_2020_exposed_F_uniform_coarse = aggregate_pop(pop_array_uniform_2020, hmax_F, flood_grid_transform, flood_grid_crs, region_utm, background_utm)

#%%# ============================================================================================ #
# ===================== Print summary statistics of exposed population ========================== #
# =============================================================================================== #
print("2020 Factual exposed population stats:")
print("Total population in region:", np.nansum(pop_arrays[2020]))
print("Total population in Sofala:", np.nansum(pop_sofala_arrays[2020]))
print("Total exposed people:", np.nansum(ra_exposed_pop_2020_F).astype(int))
print("Exposed people percentage of total population:", 100 * np.nansum(ra_exposed_pop_2020_F).astype(int) / np.nansum(pop_arrays[2020]))

print("\n2020 Counterfactual exposed population stats:")
print("Total population in region:", np.nansum(pop_arrays[2020]))
print("Total population in Sofala:", np.nansum(pop_sofala_arrays[2020]))
print("Total exposed people:", np.nansum(ra_exposed_pop_2020_CF).astype(int))
print("Exposed people percentage of total population:", 100 * np.nansum(ra_exposed_pop_2020_CF).astype(int) / np.nansum(pop_arrays[2020]))

print("\n1990 Factual exposed population stats:")
print("Total population in region:", np.nansum(pop_arrays[1990]))
print("Total population in Sofala:", np.nansum(pop_sofala_arrays[1990]))
print("Total exposed people:", np.nansum(ra_exposed_pop_1990_F).astype(int))
print("Exposed people percentage of total population:", 100 * np.nansum(ra_exposed_pop_1990_F).astype(int) / np.nansum(pop_arrays[1990]))

print("\n2020 UNIFORM exposed population stats:")
print("Total population in region:", np.nansum(pop_array_uniform_2020))
print("Total exposed people:", np.nansum(ra_exposed_pop_2020_F_uniform).astype(int))
print("Exposed people percentage of total population:", 100 * np.nansum(ra_exposed_pop_2020_F_uniform).astype(int) / np.nansum(pop_arrays[2020]))

print("\nOne-line attribution numbers:")
print(f"Exposed population in 2020 Factual: {int(np.nansum(ra_exposed_pop_2020_F).astype(int)):,}")
print(f"Exposed population attributable to climate change: {int(np.nansum(ra_exposed_pop_2020_F) - np.nansum(ra_exposed_pop_2020_CF)):,} {100 * (np.nansum(ra_exposed_pop_2020_F) - np.nansum(ra_exposed_pop_2020_CF)) / np.nansum(ra_exposed_pop_2020_F):.2f}%")
print(f"Exposed population attributable to population change (2020-1990): {int(np.nansum(ra_exposed_pop_2020_F) - np.nansum(ra_exposed_pop_1990_F)):,} {100 * (np.nansum(ra_exposed_pop_2020_F) - np.nansum(ra_exposed_pop_1990_F)) / np.nansum(ra_exposed_pop_2020_F):.2f}%")
print(f"Exposed population attributable to population change (uniform growth): {int(np.nansum(ra_exposed_pop_2020_F) - np.nansum(ra_exposed_pop_2020_F_uniform)):,} {100 * (np.nansum(ra_exposed_pop_2020_F) - np.nansum(ra_exposed_pop_2020_F_uniform)) / np.nansum(ra_exposed_pop_2020_F):.2f}%")
print(f"Population growth from 1990 to 2020 in the region: {int(np.nansum(pop_arrays[2020]) - np.nansum(pop_arrays[1990])):,} {100 * (np.nansum(pop_arrays[2020]) - np.nansum(pop_arrays[1990])) / np.nansum(pop_arrays[2020]):.2f}%")


# %%
# # ============================================================================================ #
# # ================== Plot comparison rasterized vs vectorized population ===================== #
# # ============================================================================================ #

# gdf_pop_2020_exposed_F_coarse_vector = aggregate_pop_sjoin_slow(pop_arrays[2020], hmax_F, flood_grid_transform, flood_grid_crs, region_utm, background_utm, factor=100)

# gdf_2020_F_vector = gdf_pop_2020_exposed_F_coarse_vector.copy()
# gdf_2020_F_vector.loc[gdf_2020_F_vector["total_population"] == 0] = np.nan

# gdf_2020_F = gdf_pop_2020_exposed_F_coarse.copy()
# gdf_2020_F.loc[gdf_2020_F["total_population"] == 0] = np.nan

# # Ensure vector and raster GDFs have only the necessary columns
# gdf_vector = gdf_2020_F_vector[['geometry', 'total_population']].copy()
# gdf_raster = gdf_2020_F[['geometry', 'total_population']].copy()

# # Add temporary IDs
# gdf_vector['tmp_id'] = np.arange(len(gdf_vector))
# gdf_raster['tmp_id'] = np.arange(len(gdf_raster))

# # Spatial join using centroid matching for simplicity
# gdf_vector['centroid'] = gdf_vector.geometry.centroid
# gdf_raster['centroid'] = gdf_raster.geometry.centroid

# # Merge on centroid coordinates (rounded to avoid floating point issues)
# gdf_vector['cx'] = gdf_vector.centroid.x.round(2)
# gdf_vector['cy'] = gdf_vector.centroid.y.round(2)
# gdf_raster['cx'] = gdf_raster.centroid.x.round(2)
# gdf_raster['cy'] = gdf_raster.centroid.y.round(2)

# gdf_compare = gdf_vector.merge(
#     gdf_raster[['cx','cy','total_population']],
#     on=['cx','cy'],
#     suffixes=('_vector','_raster')
# )

# # Compute difference
# gdf_compare['change_in_population'] = (
#     gdf_compare['total_population_vector'] - gdf_compare['total_population_raster']
# )

# # Plotting
# fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True) # 10,6

# gdf_2020_F_vector = gdf_pop_2020_exposed_F_coarse_vector.copy()
# gdf_2020_F_vector.loc[gdf_2020_F_vector["total_population"] == 0] = np.nan

# gdf_2020_F = gdf_pop_2020_exposed_F_coarse.copy()
# gdf_2020_F.loc[gdf_2020_F["total_population"] == 0] = np.nan

# gdf_2020_F['change_in_population_vect'] = gdf_2020_F_vector['total_population'] - gdf_2020_F['total_population']

# # Define colormap: from white to #67CBE4
# cmap = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#FC6F37"])
# norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_pop_2020_exposed_F_coarse_vector['total_population']))  
# cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#BD2A2A"])
# norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_compare['change_in_population']))

# plot = gdf_2020_F_vector.plot(column="total_population", cmap=cmap, norm=norm ,
#                        linewidth=0.1, edgecolor="grey", ax=axes[0], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

# plot = gdf_2020_F.plot(column="total_population", cmap=cmap, norm=norm,
#                         linewidth=0.1, edgecolor="grey", ax=axes[1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

# plot = gdf_compare.plot(column="change_in_population", cmap=cmap_change, norm=norm_change,
#                        linewidth=0.1, edgecolor="grey", ax=axes[2], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

# for ax in axes:
#     region.boundary.plot(ax=ax, color='black', linewidth=1)

# xmin, xmax, ymin, ymax = flood_extent
# for i, ax in enumerate(axes):
#     background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
#     bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
#     region_utm.boundary.plot(ax=ax, color='black', linewidth=0.5, zorder=2)
#     ax.set_xlim(xmin, xmax)
#     ax.set_ylim(ymin, ymax)

# for i, ax in enumerate(axes[:2]):
#     sm = ScalarMappable(cmap=cmap, norm=norm)
#     sm._A = []  
#     cbar = plt.colorbar(sm, ax=ax, shrink=0.8)
#     cbar.set_label("Population (people per cell)")

# sm = ScalarMappable(cmap=cmap_change, norm=norm_change)
# sm._A = []  
# cbar = plt.colorbar(sm, ax=axes[2], shrink=0.8)
# cbar.set_label("Difference in population (people per cell)")

# axes[0].set_title("Factual 2020 vectorized", fontsize=10)
# axes[1].set_title("Factual 2020 rasterized", fontsize=10)
# axes[2].set_title("Vectorized - Rasterized", fontsize=10)

# axes[0].set_ylabel("y coordinate UTM zone 36S [m]")
# axes[0].set_xlabel("x coordinate UTM zone 36S [m]")
# axes[1].set_xlabel("x coordinate UTM zone 36S [m]")
# axes[2].set_xlabel("x coordinate UTM zone 36S [m]")

# fig.suptitle("Total aggregated population in study region", fontsize=12)

# plt.tight_layout()
# plt.show()



#%% ============================================================================================ # 
# ================== Plot distribution of flood depth per exposed population =================== #
# ============================================================================================== #
def compute_cdf_and_bins(gdf, bins, depth_col="flood_depth", pop_col="population"):
    """Return sorted (depth, cdf) arrays and population counts per depth bin."""
    flood = gdf[depth_col].values
    pop   = gdf[pop_col].values

    # Mask invalid
    mask = ~np.isnan(flood) & (pop > 0)
    flood, pop = flood[mask], pop[mask]

    # Continuous CDF
    idx = np.argsort(flood)
    flood_sorted = flood[idx]
    pop_sorted   = pop[idx]
    cdf = np.cumsum(pop_sorted) / np.sum(pop_sorted)

    # Binned population
    pop_by_depth = pd.Series(pop).groupby(pd.cut(flood, bins)).sum()

    return flood_sorted, cdf, pop_by_depth

# --- Inputs ---
bins = np.arange(0, 3.5 + 0.25, 0.25)

# --- Compute for both datasets ---
flood_F_2020_sorted, cdf_F, pop_2020_by_depth_F   = compute_cdf_and_bins(df_pop_2020_exposed_F, bins)
flood_CF_2020_sorted, cdf_CF_clim, pop_2020_by_depth_CF = compute_cdf_and_bins(df_pop_2020_exposed_CF, bins)
flood_F_1990_sorted, cdf_CF_pop, pop_1990_by_depth_F = compute_cdf_and_bins(df_pop_1990_exposed_F, bins)
flood_CF_1990_sorted, cdf_CF_pop_clim, pop_1990_by_depth_CF = compute_cdf_and_bins(df_pop_1990_exposed_CF, bins)

bins_fine = np.arange(0, 3.5 + 0.1, 0.1)
flood_F_2020_sorted_fine, cdf_F_fine, pop_2020_by_depth_F_fine   = compute_cdf_and_bins(df_pop_2020_exposed_F, bins_fine)
flood_CF_2020_sorted_fine, cdf_CF_clim_fine, pop_2020_by_depth_CF_fine = compute_cdf_and_bins(df_pop_2020_exposed_CF, bins_fine)
flood_F_1990_sorted_fine, cdf_CF_pop_fine, pop_1990_by_depth_F_fine = compute_cdf_and_bins(df_pop_1990_exposed_F, bins_fine)
flood_CF_1990_sorted_fine, cdf_CF_pop_clim_fine, pop_1990_by_depth_CF_fine = compute_cdf_and_bins(df_pop_1990_exposed_CF, bins_fine)

# --- Plot continuous CDF ---
plt.figure(figsize=(8,5))
plt.plot(flood_F_2020_sorted, cdf_F, label="Factual")
plt.plot(flood_CF_2020_sorted, cdf_CF_clim, label="Counterfactual Climate")
plt.plot(flood_F_1990_sorted, cdf_CF_pop, label="Counterfactual Population")
plt.plot(flood_CF_1990_sorted, cdf_CF_pop_clim, label="Counterfactual Climate & Population")
plt.xlim(0, 3.5)
plt.xlabel("Flood depth (m)")
plt.ylabel("Cumulative fraction of exposed population")
plt.title("Continuous CDF of exposed population by flood depth")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.show()

# --- Plot binned bar chart ---
bin_centers = bins[:-1] + np.diff(bins) / 2
plt.figure(figsize=(8,5))
plt.bar(bin_centers - 0.05, pop_2020_by_depth_F.values, width=0.05, label="Factual", alpha=0.7)
plt.bar(bin_centers, pop_2020_by_depth_CF.values, width=0.05, label="Counterfactual Climate", alpha=0.7)
plt.bar(bin_centers + 0.05, pop_1990_by_depth_F.values, width=0.05, label="Counterfactual Population", alpha=0.7)
plt.bar(bin_centers + 0.1, pop_1990_by_depth_CF.values, width=0.05, label="Counterfactual Climate & Population", alpha=0.7)

plt.xlabel("Flood depth (m)")
plt.ylabel("Exposed population")
plt.title("Population exposed by flood depth bins")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.xlim(0, 3.5)
plt.show()


# --- PLot absolute line plot ---
bin_centers = bins_fine[:-1] + np.diff(bins_fine) / 2
plt.figure(figsize=(8,5))
plt.plot(bin_centers, pop_2020_by_depth_F_fine.values, label="Factual")
plt.plot(bin_centers, pop_2020_by_depth_CF_fine.values, label="Counterfactual Climate")
plt.plot(bin_centers, pop_1990_by_depth_F_fine.values, label="Counterfactual Population")
plt.plot(bin_centers, pop_1990_by_depth_CF_fine.values, label="Counterfactual Climate & Population")
plt.xlabel("Flood depth (m)")  
plt.ylabel("Exposed population")
plt.title("Population exposed by flood depth bins") 
plt.legend()
plt.xlim(0, 3.5)
plt.grid(True, linestyle="--", alpha=0.5)

#%%
# ================================================================================================== #
# =========================== Plot density-weighted population exposure ============================ #
# ================================================================================================== #


x = np.linspace(0, 3.5, 200)

kde = gaussian_kde(df_pop_2020_exposed_F['flood_depth'].values, weights=df_pop_2020_exposed_F['population'].values)
y_F = kde(x) # probability density
# Find peaks and sort by height
peaks, _ = find_peaks(y_F)
sorted_peaks_F = peaks[np.argsort(y_F[peaks])[::-1]]

kde = gaussian_kde(df_pop_2020_exposed_CF['flood_depth'].values, weights=df_pop_2020_exposed_CF['population'].values)
y_CF_clim = kde(x) # probability density
# Find peaks and sort by height
peaks, _ = find_peaks(y_CF_clim)
sorted_peaks_CF_clim = peaks[np.argsort(y_CF_clim[peaks])[::-1]]

kde = gaussian_kde(df_pop_1990_exposed_F['flood_depth'].values, weights=df_pop_1990_exposed_F['population'].values)
y_CF_pop = kde(x) # probability density
# Find peaks and sort by height
peaks, _ = find_peaks(y_CF_pop)
sorted_peaks_CF_pop = peaks[np.argsort(y_CF_pop[peaks])[::-1]]

kde = gaussian_kde(df_pop_1990_exposed_CF['flood_depth'].values, weights=df_pop_1990_exposed_CF['population'].values)
y_CF_clim_pop = kde(x) 

#%%
# Plot density-weighted population exposure per flood depth
plt.figure(figsize=(8,5))
plt.plot(x, y_F, label="Factual")
plt.plot(x, y_CF_clim, label="Counterfactual Climate")
plt.plot(x, y_CF_pop, label="Counterfactual Population")
plt.plot(x, y_CF_clim_pop, label="Counterfactual Climate & Population")

if len(sorted_peaks_F) > 1:
    second_peak_x = x[sorted_peaks_F[1]]
    plt.axvline(second_peak_x, linestyle="--", color = '#1f77b4', linewidth=0.8)
for i, idx in enumerate(sorted_peaks_F[:2]):  # only top 2 peaks
    plt.annotate(f"{x[idx]:.2f} m",
                 xy=(x[idx], y_F[idx]),
                 xytext=(x[idx]+0.1, y_F[idx]+0.001), color='#1f77b4', fontsize=9,
                 arrowprops=dict(arrowstyle="->", lw=1, color='#1f77b4'))

if len(sorted_peaks_CF_clim) > 1:
    second_peak_x = x[sorted_peaks_CF_clim[1]]
    plt.axvline(second_peak_x, linestyle="--", color='#ff7f0e', linewidth=0.8)
for i, idx in enumerate(sorted_peaks_CF_clim[:2]):  # only top 2 peaks
    plt.annotate(f"{x[idx]:.2f} m",
                 xy=(x[idx], y_CF_clim[idx]),
                 xytext=(x[idx]+0.1, y_CF_clim[idx]+0.001), color='#ff7f0e', fontsize=9,    
                 arrowprops=dict(arrowstyle="->", lw=1, color='#ff7f0e'))

for i, idx in enumerate(sorted_peaks_CF_pop[:2]):  # only top 2 peaks
    plt.annotate(f"{x[idx]:.2f} m",
                 xy=(x[idx], y_CF_pop[idx]),
                 xytext=(x[idx]+0.1, y_CF_pop[idx]+0.001), color='#2ca02c', fontsize=9,
                 arrowprops=dict(arrowstyle="->", lw=1, color='#2ca02c'))

plt.xlabel("Flood depth (m)")
plt.ylabel("Density-weighted population exposure")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.xlim(0, 3.5)


#%%
#%% ============================================================================================ # 
# =================== Plot the flood depth per exposed population spatially ==================== #
# ============================================================================================== #
# Plot average flood depth per exposed population cell
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True) # 10,6

# Plot average flood depth among exposed population
gdf_2020_F = gdf_pop_2020_exposed_F_coarse.copy()
gdf_2020_F.loc[gdf_2020_F["exposed_population"] <= 0, "avg_flood_depth"] = np.nan

gdf_2020_CF = gdf_pop_2020_exposed_CF_coarse.copy()
gdf_2020_CF.loc[gdf_2020_CF["exposed_population"] <= 0, "avg_flood_depth"] = np.nan

gdf_2020_CF['change_in_flood_depth'] = gdf_2020_F['avg_flood_depth'] - gdf_2020_CF['avg_flood_depth']

# Define colormap: from white to #67CBE4
cmap = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#67CBE4"])
cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#651F94"])

plot = gdf_2020_F.plot(column="avg_flood_depth", cmap=cmap, vmin=0, vmax=3.5, linewidth=0.1, 
                edgecolor="grey", ax=axes[0], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2020_CF.plot(column="avg_flood_depth", cmap=cmap, vmin=0, vmax=3.5, linewidth=0.1, 
                edgecolor="grey", ax=axes[1], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2020_CF.plot(column="change_in_flood_depth", cmap=cmap_change, vmin=0, vmax=0.5, linewidth=0.1, 
                edgecolor="grey", ax=axes[2], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

xmin, xmax, ymin, ymax = flood_extent
for i, ax in enumerate(axes):
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
    region_utm.boundary.plot(ax=ax, color='black', linewidth=0.5, zorder=2)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

for i, ax in enumerate(axes[:2]):
    sm = ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=3.5))
    sm._A = []  
    cbar = plt.colorbar(sm, ax=ax, shrink=0.8)
    cbar.set_label("Average flood depth among exposed population (m)")

sm = ScalarMappable(cmap=cmap_change, norm=plt.Normalize(vmin=0, vmax=0.5))
sm._A = []  
cbar = plt.colorbar(sm, ax=axes[2], shrink=0.8)
cbar.set_label("Difference in Average flood depth \namong exposed population (m)")

axes[0].set_title("Factual", fontsize=10)
axes[1].set_title("No Climate Change", fontsize=10)
axes[2].set_title("Factual - Counterfactual", fontsize=10)

axes[0].set_ylabel("y coordinate UTM zone 36S [m]")
axes[0].set_xlabel("x coordinate UTM zone 36S [m]")
axes[1].set_xlabel("x coordinate UTM zone 36S [m]")
axes[2].set_xlabel("x coordinate UTM zone 36S [m]")

fig.suptitle("Average flood depth among exposed population", fontsize=12)

plt.tight_layout()
plt.show()

#%% ============================================================================================ # 
# ============================= Plot the change in population spatially ======================== #
# ============================================================================================== #
# Plot average flood depth per exposed population cell
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True) # 10,6

gdf_2020_F = gdf_pop_2020_exposed_CF_coarse.copy()
gdf_2020_F.loc[gdf_2020_F["total_population"] == 0] = np.nan

gdf_1990_F = gdf_pop_1990_exposed_CF_coarse.copy()
gdf_1990_F.loc[gdf_1990_F["total_population"] == 0] = np.nan

gdf_1990_F['change_in_population'] = gdf_2020_F['total_population'] - gdf_1990_F['total_population']

# Define colormap: from white to #67CBE4
cmap = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#FC6F37"])
norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_pop_2020_exposed_F_coarse['total_population']))  
cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#BD2A2A"])
norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_1990_F['change_in_population']))

plot = gdf_2020_F.plot(column="total_population", cmap=cmap, norm=norm ,
                       linewidth=0.1, edgecolor="grey", ax=axes[0], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_1990_F.plot(column="total_population", cmap=cmap, norm=norm,
                        linewidth=0.1, edgecolor="grey", ax=axes[1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_1990_F.plot(column="change_in_population", cmap=cmap_change, norm=norm_change,
                       linewidth=0.1, edgecolor="grey", ax=axes[2], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

for ax in axes:
    region.boundary.plot(ax=ax, color='black', linewidth=1)

xmin, xmax, ymin, ymax = flood_extent
for i, ax in enumerate(axes):
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
    region_utm.boundary.plot(ax=ax, color='black', linewidth=0.5, zorder=2)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

for i, ax in enumerate(axes[:2]):
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []  
    cbar = plt.colorbar(sm, ax=ax, shrink=0.8)
    cbar.set_label("Population (people per cell)")

sm = ScalarMappable(cmap=cmap_change, norm=norm_change)
sm._A = []  
cbar = plt.colorbar(sm, ax=axes[2], shrink=0.8)
cbar.set_label("Difference in population (people per cell)")

axes[0].set_title("Factual", fontsize=10)
axes[1].set_title("No population growth", fontsize=10)
axes[2].set_title("Factual - Counterfactual", fontsize=10)

axes[0].set_ylabel("y coordinate UTM zone 36S [m]")
axes[0].set_xlabel("x coordinate UTM zone 36S [m]")
axes[1].set_xlabel("x coordinate UTM zone 36S [m]")
axes[2].set_xlabel("x coordinate UTM zone 36S [m]")

fig.suptitle("Total population in study region", fontsize=12)

plt.tight_layout()
plt.show()


#%%
#%% ============================================================================================ # 
# ================ Plot the diff in uniform and spatial change in population =================== #
# ============================================================================================== #
# Plot average flood depth per exposed population cell
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True) # 10,6

gdf_2020_F = gdf_pop_2020_exposed_CF_coarse.copy()
gdf_2020_F.loc[gdf_2020_F["total_population"] == 0] = np.nan

gdf_2020_F_uni = gdf_pop_2020_exposed_F_uniform_coarse.copy()
gdf_2020_F_uni.loc[gdf_2020_F_uni["total_population"] == 0] = np.nan

gdf_2020_F_uni['change_in_population'] = gdf_2020_F['total_population'] - gdf_2020_F_uni['total_population']

# Define colormap: from white to #67CBE4
cmap = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#FC6F37"])
norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_pop_2020_exposed_F_coarse['total_population']))  
cmap_change = plt.get_cmap('RdBu_r')
norm_change = mcolors.TwoSlopeNorm(vmin=np.nanmin(gdf_2020_F_uni['change_in_population']),
                                   vcenter=0,
                                   vmax=np.nanmax(gdf_2020_F_uni['change_in_population']))

plot = gdf_2020_F.plot(column="total_population", cmap=cmap, norm=norm ,
                       linewidth=0.1, edgecolor="grey", ax=axes[0], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2020_F_uni.plot(column="total_population", cmap=cmap, norm=norm,
                        linewidth=0.1, edgecolor="grey", ax=axes[1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2020_F_uni.plot(column="change_in_population", cmap=cmap_change, norm=norm_change,
                       linewidth=0.1, edgecolor="grey", ax=axes[2], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

for ax in axes:
    region.boundary.plot(ax=ax, color='black', linewidth=1)

xmin, xmax, ymin, ymax = flood_extent
for i, ax in enumerate(axes):
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
    region_utm.boundary.plot(ax=ax, color='black', linewidth=0.5, zorder=2)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

for i, ax in enumerate(axes[:2]):
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []  
    cbar = plt.colorbar(sm, ax=ax, shrink=0.8)
    cbar.set_label("Population (people per cell)")

sm = ScalarMappable(cmap=cmap_change, norm=norm_change)
sm._A = []  
cbar = plt.colorbar(sm, ax=axes[2], shrink=0.8)
cbar.set_label("Difference in population (people per cell)")

axes[0].set_title("Factual 2020", fontsize=10)
axes[1].set_title("Uniform population growth 2020", fontsize=10)
axes[2].set_title("Factual - Uniform Factual", fontsize=10)

axes[0].set_ylabel("y coordinate UTM zone 36S [m]")
axes[0].set_xlabel("x coordinate UTM zone 36S [m]")
axes[1].set_xlabel("x coordinate UTM zone 36S [m]")
axes[2].set_xlabel("x coordinate UTM zone 36S [m]")

fig.suptitle("Total population in study region", fontsize=12)

plt.tight_layout()
plt.show()



#%% ============================================================================================ # 
# =============== Plot the change in exposed population > 1 m flood depth spatially ============ #
# ============================================================================================== #
# Plot average flood depth per exposed population cell
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True) # 10,6

# Plot average flood depth among exposed population
gdf_2020_F = gdf_pop_2020_exposed_F_coarse.copy()
gdf_2020_F.loc[gdf_2020_F["exposed_population"] <= 1, "avg_flood_depth"] = np.nan

gdf_2020_CF = gdf_pop_2020_exposed_CF_coarse.copy()
gdf_2020_CF.loc[gdf_2020_CF["exposed_population"] <= 1, "avg_flood_depth"] = np.nan

gdf_2020_CF['change_in_exposed_population_>1m'] = gdf_2020_F['exposed_population'] - gdf_2020_CF['exposed_population']

# Define colormap: from white to #67CBE4
cmap = 'Blues'
norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2020_F['exposed_population']))
cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#651F94"])
norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2020_CF['change_in_exposed_population_>1m']))

plot = gdf_2020_F.plot(column="exposed_population", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[0], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2020_CF.plot(column="exposed_population", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[1], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2020_CF.plot(column="change_in_exposed_population_>1m", cmap=cmap_change, norm=norm_change, linewidth=0.1, 
                edgecolor="grey", ax=axes[2], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

xmin, xmax, ymin, ymax = flood_extent
for i, ax in enumerate(axes):
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
    region_utm.boundary.plot(ax=ax, color='black', linewidth=0.5, zorder=2)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

for i, ax in enumerate(axes[:2]):
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []  
    cbar = plt.colorbar(sm, ax=ax, shrink=0.8)
    cbar.set_label("Exposed population > 1 m flood depth")

sm = ScalarMappable(cmap=cmap_change, norm=norm_change)
sm._A = []  
cbar = plt.colorbar(sm, ax=axes[2], shrink=0.8)
cbar.set_label("Difference in Exposed population > 1 m flood depth")

axes[0].set_title("Factual", fontsize=10)
axes[1].set_title("No Climate Change", fontsize=10)
axes[2].set_title("Factual - Counterfactual", fontsize=10)

axes[0].set_ylabel("y coordinate UTM zone 36S [m]")
axes[0].set_xlabel("x coordinate UTM zone 36S [m]")
axes[1].set_xlabel("x coordinate UTM zone 36S [m]")
axes[2].set_xlabel("x coordinate UTM zone 36S [m]")

fig.suptitle("Population exposed to > 1 m flood depth", fontsize=12)

plt.tight_layout()
plt.show()

#%% ============================================================================================ # 
# ============================= Plot % cells with flood depth > 1 m ============================ #
# ============================================================================================== #
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot average flood depth among exposed population
gdf_2020_F = gdf_pop_2020_exposed_F_coarse.copy()
gdf_2020_CF = gdf_pop_2020_exposed_CF_coarse.copy()

gdf_2020_CF['change_in_%more_1m'] = gdf_2020_F['pct_cells_higher_1m'] - gdf_2020_CF['pct_cells_higher_1m']

# Define colormap: from white to #67CBE4
cmap = 'Blues'
norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2020_F['pct_cells_higher_1m']))
cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#651F94"])
norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2020_CF['change_in_%more_1m']))

plot = gdf_2020_F.plot(column="pct_cells_higher_1m", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[0], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2020_CF.plot(column="pct_cells_higher_1m", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[1], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2020_CF.plot(column="change_in_%more_1m", cmap=cmap_change, norm=norm_change, linewidth=0.1, 
                edgecolor="grey", ax=axes[2], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

xmin, xmax, ymin, ymax = flood_extent
for i, ax in enumerate(axes):
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
    region_utm.boundary.plot(ax=ax, color='black', linewidth=0.5, zorder=2)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

for i, ax in enumerate(axes[:2]):
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []  
    cbar = plt.colorbar(sm, ax=ax, shrink=0.8)
    cbar.set_label("% cells > 1 m flood depth")

sm = ScalarMappable(cmap=cmap_change, norm=norm_change)
sm._A = []  
cbar = plt.colorbar(sm, ax=axes[2], shrink=0.8)
cbar.set_label("Difference in % cells > 1 m flood depth")

axes[0].set_title("Factual", fontsize=10)
axes[1].set_title("No Climate Change", fontsize=10)
axes[2].set_title("Factual - Counterfactual", fontsize=10)

axes[0].set_ylabel("y coordinate UTM zone 36S [m]")
axes[0].set_xlabel("x coordinate UTM zone 36S [m]")
axes[1].set_xlabel("x coordinate UTM zone 36S [m]")
axes[2].set_xlabel("x coordinate UTM zone 36S [m]")

fig.suptitle("% cells > 1 m flood depth", fontsize=12)

plt.tight_layout()
plt.show()
#%% ============================================================================================ # 
# ======================== Plot the factual flood and exposed population ======================= #
# ============================================================================================== #
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

for i, ax in enumerate(axes):
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
    region_utm.boundary.plot(ax=ax, color='black', linewidth=1, zorder=2)

# Factual flood
im = axes[0].imshow(hmax_F, cmap='viridis', extent=flood_extent, origin='lower', vmin=0, vmax=3.5)
cbar = plt.colorbar(im, ax=axes[0], shrink=0.8)
cbar.set_label("Flood depth (m)")
axes[0].set_title("Factual Flood Depth")

# Exposed population
im = axes[1].imshow(ra_exposed_pop_2020_F, cmap='viridis', extent=flood_extent, origin='lower')
cbar = plt.colorbar(im, ax=axes[1], shrink=0.8)
cbar.set_label("Exposed population")
axes[1].set_title("Factual Exposed Population")

plt.tight_layout()
plt.show()


#%%
# ============================================================================================ #
# ======================= Compute total exposed population per district ====================== #
# ============================================================================================ #
pop_totals = []
pop_exposed = []
pop_per_district = []

for _, row in districts_filtered.iterrows():
    # Mask raster to district polygon
    from rasterio import features
    district_mask = features.geometry_mask([row.geometry],
                                           out_shape=ra_exposed_pop_2020_F.shape,
                                           transform=flood_grid_transform,
                                           invert=True)
    pop_exposed.append(ra_exposed_pop_2020_F[district_mask].sum())
    pop_totals.append(pop_arrays[2020][district_mask].sum())

districts_filtered['pop_exposed'] = pop_exposed
districts_filtered['pop_total'] = pop_totals

# Make sure districts are in the same CRS as the original pop raster
with rasterio.open(population_raster_path_2020) as src:
    pop_crs = src.crs

districts_native = districts_filtered.to_crs(pop_crs)

for _, row in districts_native.iterrows():
    district_mask = features.geometry_mask(
        [row.geometry],
        out_shape=pop_sofala_districts[2020][0].shape,
        transform=pop_affine_sofala_districts[2020],
        invert=True
    )
    pop_per_district.append(pop_sofala_districts[2020][0][district_mask].sum())

districts_filtered['pop_per_district'] = pop_per_district

#%%
# --- Plot ---
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

for ax in axes:
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
    districts_filtered.boundary.plot(ax=ax, color='orange', linewidth=2, zorder=2)
    region_utm.boundary.plot(ax=ax, color='lightblue', linewidth=0.5, zorder=2)

# Rasterize the case-study region polygon
region_mask = rasterio.features.rasterize(
    [(geom, 1) for geom in region.geometry],
    out_shape=ra_exposed_pop_2020_F.shape,
    transform=flood_grid_transform,
    fill=0,
    all_touched=True,
    dtype=np.uint8
).astype(bool)

# Mask raster outside region
ra_exposed_pop_masked = np.where(region_mask, ra_exposed_pop_2020_F, np.nan)
    
# Exposed population
im = axes[0].imshow(ra_exposed_pop_masked, cmap='viridis', extent=flood_extent, origin='lower')
cbar = plt.colorbar(im, ax=axes[0], shrink=0.8)
cbar.set_label("Exposed population")
axes[0].set_title("Factual Exposed Population")

# Exposed population
im = axes[1].imshow(pop_arrays[2020], cmap='Reds', extent=flood_extent, origin='lower')
cbar = plt.colorbar(im, ax=axes[1], shrink=0.8)
cbar.set_label("Total population")
axes[1].set_title("Total 2020 Population")

# Dictionary with (dx, dy) offsets in map units for specific districts
label_offsets = {
    "Sofala": (0, 10000),  # move 0 m east, 10000 m north
    "Nhamatanda": (18000, -22000),
    "Estaquinha": (20000, -5000),
    # add more as needed
}

outside_districts = ["Nhamatanda", "Estaquinha"]

for idx, row in districts_filtered.iterrows():
    x, y = row.geometry.centroid.x, row.geometry.centroid.y
    # Apply offset if district in dictionary
    if row['NAME_3'] in label_offsets:
        dx, dy = label_offsets[row['NAME_3']]
        x += dx
        y += dy
    axes[0].text(x, y, f"{row['pop_exposed']:,.0f}", fontsize=8, ha='center', va='center',
            color='white', fontweight='bold', zorder=5)
    
    axes[1].text(x, y, f"{row['pop_total']:,.0f}", fontsize=8, ha='center', va='center',
            color='black', fontweight='bold', zorder=5)
    
    # District name label
    name_x = x - 8000 if row['NAME_3'] in outside_districts else x
    name_y = y - 3000  # all 3000 south of pop label

    axes[0].text(
        name_x, name_y, row['NAME_3'], fontsize=8, ha='center', va='center',
        color='coral', fontweight='bold', zorder=5
    )
    axes[1].text(
        name_x, name_y, row['NAME_3'], fontsize=8, ha='center', va='center',
        color='coral', fontweight='bold', zorder=5
    )
    
plt.tight_layout()

#%%
# --- Plot ---
fig, ax = plt.subplots(1, 1, figsize=(12, 6))

# background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
# bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
districts_filtered.boundary.plot(ax=ax, color='orange', linewidth=2, zorder=2)
region_utm.boundary.plot(ax=ax, color='lightblue', linewidth=1, zorder=2)

# Dictionary with (dx, dy) offsets in map units for specific districts
label_offsets = {
    "Sofala": (0, 10000),  # move 0 m east, 10000 m north
    "Nhamatanda": (18000, -22000),
    "Estaquinha": (20000, -5000),
    # add more as needed
}
outside_districts = ["Nhamatanda", "Estaquinha"]

for idx, row in districts_filtered.iterrows():
    x, y = row.geometry.centroid.x, row.geometry.centroid.y
    ax.text(x, y, f"{row['pop_per_district']:,.0f}", fontsize=8, ha='center', va='center',
            color='black', fontweight='bold', zorder=5)
    
    # District name label
    name_x = x
    name_y = y - 4000  # all 3000 south of pop label

    ax.text(
        name_x, name_y, row['NAME_3'], fontsize=8, ha='center', va='center',
        color='grey', fontweight='bold', zorder=5
    )

ax.set_title("2020 District Population")

plt.tight_layout()

#%%
# table with numbers per district
df_district_summary = districts_filtered[['NAME_3', 'pop_per_district', 'pop_total', 'pop_exposed']].copy()
df_district_summary.columns = ['District', 'District Population (2020)', 'Total Population in Region', 'Exposed Population (2020 Factual)']
df_district_summary["% of district pop"] = (100 * df_district_summary['Total Population in Region'] / df_district_summary['District Population (2020)']).round(0)
df_district_summary["% exposed"] = (100 * df_district_summary['Exposed Population (2020 Factual)'] / df_district_summary['Total Population in Region']).round(0)
df_district_summary[['District', 'District Population (2020)', 'Total Population in Region', 'Exposed Population (2020 Factual)']] = df_district_summary[['District', 'District Population (2020)', 'Total Population in Region', 'Exposed Population (2020 Factual)']].round(0)

df_district_summary.to_csv("c:/Code/COMPASS_exposure/Data/Modified/sofala_district_exposed_population_summary.csv", index=False)

#%% ============================================================================================ # 
# --- Prepare colormap for absolute pop exposed ---
colors = plt.cm.Blues(np.linspace(0, 1, 256))
colors[0] = [1, 1, 1, 0]  # make first color transparent
pop_cmap = mcolors.ListedColormap(colors)

# --- Prepare colormap for change in pop exposed ---
colors = plt.cm.RdBu_r(np.linspace(0, 1, 256))
mid = 128  # midpoint index in 256
colors[mid] = [1, 1, 1, 0]  # RGBA, fully transparent
diff_cmap = mcolors.ListedColormap(colors)

pop_diff_PG = ra_exposed_pop_2020_F - ra_exposed_pop_1990_F
pop_diff_CC = ra_exposed_pop_2020_F - ra_exposed_pop_2020_CF

vmax_pg = np.nanmax(pop_diff_PG)
vmax_cc = np.nanmax(pop_diff_CC)

#%%
# --- Setup figure ---
fig, axes = plt.subplots(2, 3, figsize=(18, 12), sharex=True, sharey=True, dpi=300, constrained_layout=True)

for i, ax in enumerate(axes.flatten()):
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
    region_utm.boundary.plot(ax=ax, color='black', linewidth=1, zorder=2)

# --- Population colormap ---
pop_bins = np.arange(0, np.nanmax(ra_exposed_pop_2020_F)+1, 1)  # 0,1,2,... max
pop_colors = plt.cm.Blues(np.linspace(0, 1, len(pop_bins)-1))
pop_colors[0] = [1, 1, 1, 0]  # RGBA
pop_cmap_discrete = mcolors.ListedColormap(pop_colors)
pop_norm = BoundaryNorm(pop_bins, pop_cmap_discrete.N, extend='neither')

# --- Difference colormap ---
diff_bins = np.arange(-int(vmax_cc), int(vmax_cc)+1, 1)  # integer steps
diff_colors = plt.cm.RdBu_r(np.linspace(0, 1, len(diff_bins)-1))
mid_idx = np.where(diff_bins[:-1] == 0)[0]
if len(mid_idx) > 0:
    diff_colors[mid_idx[0]] = [1, 1, 1, 0]
diff_cmap_discrete = mcolors.ListedColormap(diff_colors)
diff_norm = BoundaryNorm(diff_bins, diff_cmap_discrete.N, extend='neither')

# --- TOP ROW ---
# Population change effect
im = axes[0,0].imshow(np.round(ra_exposed_pop_2020_F).astype(int), cmap=pop_cmap_discrete, norm=pop_norm,
                      extent=flood_extent, origin='lower', alpha=0.8, zorder=3)
axes[0,0].set_title("Exposed 2020 Population")

im = axes[0,1].imshow(np.round(ra_exposed_pop_1990_F).astype(int), cmap=pop_cmap_discrete, norm=pop_norm,
                      extent=flood_extent, origin='lower', alpha=0.8, zorder=3)
cbar = plt.colorbar(im, ax=axes[0,1], shrink=0.8)
cbar.set_label("Exposed population")
axes[0,1].set_title("Exposed 1990 Population")

im = axes[0,2].imshow(np.round(pop_diff_PG).astype(int), cmap=diff_cmap_discrete, norm=diff_norm,
                      extent=flood_extent, origin='lower', alpha=0.8, zorder=3)
cbar = plt.colorbar(im, ax=axes[0,2], shrink=0.8)
cbar.set_label("Attributable exposed population")
axes[0,2].set_title("Population Change (2020 - 1990)")

# # --- BOTTOM ROW --
# Climate change effect
im = axes[1,0].imshow(np.round(ra_exposed_pop_2020_F).astype(int), cmap=pop_cmap_discrete, norm=pop_norm, extent=flood_extent,
                      origin='lower', alpha=0.8, zorder=3)
axes[1,0].set_title("Factual Climate")

im = axes[1,1].imshow(np.round(ra_exposed_pop_2020_CF).astype(int), cmap=pop_cmap_discrete, norm=pop_norm, extent=flood_extent,
                      origin='lower', alpha=0.8, zorder=3)
cbar = plt.colorbar(im, ax=axes[1,1], shrink=0.8)
cbar.set_label("Exposed population")
axes[1,1].set_title("Counterfactual Climate")

im = axes[1,2].imshow(np.round(pop_diff_CC).astype(int), cmap=diff_cmap_discrete, norm=diff_norm, 
                      extent=flood_extent, origin='lower', alpha=0.8, zorder=3)
cbar = plt.colorbar(im, ax=axes[1,2], shrink=0.8)
cbar.set_label("Attributable exposed population")
axes[1,2].set_title("Climate Change")

plt.show()



#%%
# ============================================================================================ #
# ===================== Differences in AGGREGATED exposed population ========================= #
# ============================================================================================ #
print("Plotting spatially aggregated exposed population change")
# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5), dpi=300, sharey=True, 
                         constrained_layout=True, subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
im = axes[0].imshow(hmax_F, cmap='viridis', extent=flood_extent, origin='lower', 
                    vmin=0, vmax=3.5, zorder=2)

gdf_pop_2020_exposed_F_coarse[gdf_pop_2020_exposed_F_coarse['exposed_population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
norm = PowerNorm(gamma=0.5, vmin=gdf_pop_2020_exposed_F_coarse['exposed_population'].min(), vmax=gdf_pop_2020_exposed_F_coarse['exposed_population'].max())

plot = gdf_pop_2020_exposed_F_coarse[gdf_pop_2020_exposed_F_coarse['exposed_population'] > 0].plot(column='exposed_population', cmap='Blues', edgecolor='grey',
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
vmin, vmax = gdf_pop_2020_exposed_F_coarse["exposed_population"].min(), gdf_pop_2020_exposed_F_coarse["exposed_population"].max()
sm1 = ScalarMappable(cmap="Blues", norm=norm)
sm1._A = []  # required for colorbar without passing data
cbar = fig.colorbar(sm1, ax=axes[1], orientation="vertical", shrink=0.8)
cbar.set_label("Aggregated exposed population (# people)", fontsize=10)
cbar.ax.tick_params(labelsize=9)

axes[0].set_title("Factual flooding", fontsize=11)
axes[1].set_title("Factual population exposure", fontsize=11)



#%%
print("Plotting spatially aggregated RELATIVE exposed population change")

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5), dpi=300, sharey=True, 
                         constrained_layout=True, subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
im = axes[0].imshow(hmax_F, cmap='viridis', extent=flood_extent, origin='lower', 
                    vmin=0, vmax=3.5, zorder=2)

gdf_pop_2020_exposed_F_coarse[gdf_pop_2020_exposed_F_coarse['relative_population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2020_exposed_F_coarse[gdf_pop_2020_exposed_F_coarse['relative_population'] > 0].plot(column='relative_population', cmap='Blues', edgecolor='grey',
                                                                 linewidth=0.2, ax=axes[1], legend=False, zorder=2, rasterized=True)
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
vmin, vmax = gdf_pop_2020_exposed_F_coarse["relative_population"].min(), gdf_pop_2020_exposed_F_coarse["relative_population"].max()
sm1 = ScalarMappable(cmap="Blues", norm=Normalize(vmin=vmin, vmax=vmax))
sm1._A = []  # required for colorbar without passing data
cbar = fig.colorbar(sm1, ax=axes[1], orientation="vertical", shrink=0.8)
cbar.set_label("Relative exposed population (%)", fontsize=10)
cbar.ax.tick_params(labelsize=9)

axes[0].set_title("Factual flooding", fontsize=11)
axes[1].set_title("Factual population exposure", fontsize=11)



# %%
print("Plotting spatially aggregated exposed population damage for F, CF population and diff")

gdf_pop_1990_exposed_F_coarse['population_diff'] = (gdf_pop_2020_exposed_F_coarse['exposed_population'] - gdf_pop_1990_exposed_F_coarse['exposed_population'])

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 6), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Create colormap normalization
norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2020_exposed_F_coarse["exposed_population"].max())
norm_diff = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_1990_exposed_F_coarse["population_diff"].max())
red_half = LinearSegmentedColormap.from_list("bwr_red", plt.cm.bwr(np.linspace(0.5, 1, 256)))

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
gdf_pop_2020_exposed_F_coarse[gdf_pop_2020_exposed_F_coarse['exposed_population'] == 0].plot(ax=axes[0], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2020_exposed_F_coarse[gdf_pop_2020_exposed_F_coarse['exposed_population'] > 0].plot(column='exposed_population', cmap='Blues',  edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[0], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_1990_exposed_F_coarse[gdf_pop_1990_exposed_F_coarse['exposed_population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_1990_exposed_F_coarse[gdf_pop_1990_exposed_F_coarse['exposed_population'] > 0].plot(column='exposed_population', cmap='Blues', edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[1], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_1990_exposed_F_coarse[gdf_pop_1990_exposed_F_coarse['population_diff'] <= 0].plot(ax=axes[2], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_1990_exposed_F_coarse[gdf_pop_1990_exposed_F_coarse['population_diff'] > 0].plot(column='population_diff', cmap=red_half, norm=norm_diff, edgecolor='grey', 
                                                                 linewidth=0.2, ax=axes[2], legend=False, zorder=2, missing_kwds={"color": "white", "edgecolor": "none"}, rasterized=True)
subplot_labels = ['(a)', '(b)', '(c)']

for i, ax in enumerate(axes):
    # Add model region
    region_utm.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)

    # Add background and set extent (based on actual lat/lon coordinates)
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    minx, miny, maxx, maxy = region_utm.bounds.minx.item(), region_utm.bounds.miny.item(), region_utm.bounds.maxx.item(), region_utm.bounds.maxy.item()
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
sm1 = ScalarMappable(cmap="Blues", norm=norm_pop)
sm1._A = []  # required for colorbar without passing data
cbar = fig.colorbar(sm1, ax=axes[1], orientation="vertical", shrink=0.5)
cbar.set_label("Aggregated exposed population (# people)", fontsize=10)
cbar.ax.tick_params(labelsize=9)

sm2 = ScalarMappable(cmap=red_half, norm=norm_diff)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, ax=axes[2], orientation="vertical", shrink=0.5)
cbar2.set_label("Attributable affected population [# people]", fontsize = 9)
cbar2.ax.tick_params(labelsize=8)

axes[0].set_title("Factual population \n(2020)", fontsize=9)
axes[1].set_title("Counterfactual population \n(1990)", fontsize=9)
axes[2].set_title("Factual - Counterfactual", fontsize=9)

# %%
print("Plotting spatially aggregated exposed population for F, CF climate and diff")

gdf_pop_2020_exposed_CF_coarse['population_diff'] = (gdf_pop_2020_exposed_F_coarse['exposed_population'] - gdf_pop_2020_exposed_CF_coarse['exposed_population'])

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 5), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Create colormap normalization
norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2020_exposed_F_coarse["exposed_population"].max())
norm_diff = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2020_exposed_CF_coarse["population_diff"].max())
red_half = LinearSegmentedColormap.from_list("bwr_red", plt.cm.bwr(np.linspace(0.5, 1, 256)))

cmap = LinearSegmentedColormap.from_list("custom_cmap", ["white", "#67CBE4"])

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
gdf_pop_2020_exposed_F_coarse[gdf_pop_2020_exposed_F_coarse['exposed_population'] == 0].plot(ax=axes[0], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2020_exposed_F_coarse[gdf_pop_2020_exposed_F_coarse['exposed_population'] > 0].plot(column='exposed_population', cmap='Blues',  edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[0], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_2020_exposed_CF_coarse[gdf_pop_2020_exposed_CF_coarse['exposed_population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2020_exposed_CF_coarse[gdf_pop_2020_exposed_CF_coarse['exposed_population'] > 0].plot(column='exposed_population', cmap='Blues', edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[1], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_2020_exposed_CF_coarse[gdf_pop_2020_exposed_CF_coarse['population_diff'] <= 0].plot(ax=axes[2], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2020_exposed_CF_coarse[gdf_pop_2020_exposed_CF_coarse['population_diff'] > 0].plot(column='population_diff', cmap=red_half, norm=norm_diff, edgecolor='grey', 
                                                                 linewidth=0.2, ax=axes[2], legend=False, zorder=2, missing_kwds={"color": "white", "edgecolor": "none"}, rasterized=True)
subplot_labels = ['(a)', '(b)', '(c)']

for i, ax in enumerate(axes):
    # Add model region
    region_utm.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)

    # Add background and set extent (based on actual lat/lon coordinates)
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    minx, miny, maxx, maxy = region_utm.bounds.minx.item(), region_utm.bounds.miny.item(), region_utm.bounds.maxx.item(), region_utm.bounds.maxy.item()
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
sm1 = ScalarMappable(cmap="Blues", norm=norm_pop)
sm1._A = []  # required for colorbar without passing data
cbar = fig.colorbar(sm1, ax=axes[1], orientation="vertical", shrink=0.5)
cbar.set_label("Aggregated exposed population (# people)", fontsize=10)
cbar.ax.tick_params(labelsize=9)

sm2 = ScalarMappable(cmap=red_half, norm=norm_diff)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, ax=axes[2], orientation="vertical", shrink=0.5)
cbar2.set_label("Attributable affected population [# people]", fontsize = 9)
cbar2.ax.tick_params(labelsize=8)

axes[0].set_title("Factual climate", fontsize=9)
axes[1].set_title("Counterfactual climate", fontsize=9)
axes[2].set_title("Factual - Counterfactual", fontsize=9)



# %%
# Plotting uniform population 1990 vs spatially differing 1990
print("Plotting spatially aggregated exposed population damage for F, CF population and diff")

gdf_pop_2020_exposed_F_uniform_coarse['population_diff'] = (gdf_pop_2020_exposed_F_coarse['exposed_population'] - gdf_pop_2020_exposed_F_uniform_coarse['exposed_population'])

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 6), dpi=300, sharey=True, constrained_layout=True,
                         subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Create colormap normalization
norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2020_exposed_F_coarse["exposed_population"].max())
norm_2020_diff = mcolors.TwoSlopeNorm(vmin=(gdf_pop_2020_exposed_F_uniform_coarse["population_diff"].min()), vcenter=0, vmax=(gdf_pop_2020_exposed_F_uniform_coarse["population_diff"].max()))
cmap_2020_diff = plt.get_cmap('RdBu_r')


# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
gdf_pop_2020_exposed_F_coarse[gdf_pop_2020_exposed_F_coarse['exposed_population'] == 0].plot(ax=axes[0], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2020_exposed_F_coarse[gdf_pop_2020_exposed_F_coarse['exposed_population'] > 0].plot(column='exposed_population', cmap='Blues',  edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[0], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_2020_exposed_F_uniform_coarse[gdf_pop_2020_exposed_F_uniform_coarse['exposed_population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2020_exposed_F_uniform_coarse[gdf_pop_2020_exposed_F_uniform_coarse['exposed_population'] > 0].plot(column='exposed_population', cmap='Blues', edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[1], legend=False, 
                                                                               zorder=2, rasterized=True)
plot = gdf_pop_2020_exposed_F_uniform_coarse.plot(column='population_diff', cmap=cmap_2020_diff, norm=norm_2020_diff, edgecolor='grey',
                                                  linewidth=0.2, ax=axes[2], missing_kwds={"color": "white", "edgecolor": "none"}, 
                                                  zorder=2, rasterized=True)

subplot_labels = ['(a)', '(b)', '(c)']

for i, ax in enumerate(axes):
    # Add model region
    region_utm.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)

    # Add background and set extent (based on actual lat/lon coordinates)
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    minx, miny, maxx, maxy = region_utm.bounds.minx.item(), region_utm.bounds.miny.item(), region_utm.bounds.maxx.item(), region_utm.bounds.maxy.item()
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
sm1 = ScalarMappable(cmap="Blues", norm=norm_pop)
sm1._A = []  # required for colorbar without passing data
cbar = fig.colorbar(sm1, ax=axes[1], orientation="vertical", shrink=0.5)
cbar.set_label("Aggregated exposed population (# people)", fontsize=10)
cbar.ax.tick_params(labelsize=9)

sm2 = ScalarMappable(cmap=cmap_2020_diff, norm=norm_2020_diff)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, ax=axes[2], orientation="vertical", shrink=0.5)
cbar2.set_label("Attributable affected population [# people]", fontsize = 9)
cbar2.ax.tick_params(labelsize=8)

axes[0].set_title("Factual exposed population \n(2020)", fontsize=9)
axes[1].set_title("Factual exposed population \n(2020 - uniform)", fontsize=9)
axes[2].set_title("Spatial - Uniform", fontsize=9)

# %% #################################################################
##################### diff in relative population ####################
######################################################################
print("Plotting spatially aggregated RELATIVE exposed population damage for F, CF population and diff")

gdf_pop_1990_exposed_F_coarse['population_diff'] = (gdf_pop_2020_exposed_F_coarse['relative_population'] - gdf_pop_1990_exposed_F_coarse['relative_population'])

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 6), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Create colormap normalization
norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2020_exposed_F_coarse["relative_population"].max())
norm_diff = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_1990_exposed_F_coarse["population_diff"].max())
red_half = LinearSegmentedColormap.from_list("bwr_red", plt.cm.bwr(np.linspace(0.5, 1, 256)))

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
gdf_pop_2020_exposed_F_coarse[gdf_pop_2020_exposed_F_coarse['relative_population'] == 0].plot(ax=axes[0], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2020_exposed_F_coarse[gdf_pop_2020_exposed_F_coarse['relative_population'] > 0].plot(column='relative_population', cmap='Blues',  edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[0], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_1990_exposed_F_coarse[gdf_pop_1990_exposed_F_coarse['relative_population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_1990_exposed_F_coarse[gdf_pop_1990_exposed_F_coarse['relative_population'] > 0].plot(column='relative_population', cmap='Blues', edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[1], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_1990_exposed_F_coarse[gdf_pop_1990_exposed_F_coarse['population_diff'] <= 0].plot(ax=axes[2], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_1990_exposed_F_coarse[gdf_pop_1990_exposed_F_coarse['population_diff'] > 0].plot(column='population_diff', cmap=red_half, norm=norm_diff, edgecolor='grey', 
                                                                 linewidth=0.2, ax=axes[2], legend=False, zorder=2, missing_kwds={"color": "white", "edgecolor": "none"}, rasterized=True)
subplot_labels = ['(a)', '(b)', '(c)']

for i, ax in enumerate(axes):
    # Add model region
    region_utm.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)

    # Add background and set extent (based on actual lat/lon coordinates)
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    minx, miny, maxx, maxy = region_utm.bounds.minx.item(), region_utm.bounds.miny.item(), region_utm.bounds.maxx.item(), region_utm.bounds.maxy.item()
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
sm1 = ScalarMappable(cmap="Blues", norm=norm_pop)
sm1._A = []  # required for colorbar without passing data
cbar = fig.colorbar(sm1, ax=axes[1], orientation="vertical", shrink=0.5)
cbar.set_label("Aggregated relative exposed population (# people)", fontsize=10)
cbar.ax.tick_params(labelsize=9)

sm2 = ScalarMappable(cmap=red_half, norm=norm_diff)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, ax=axes[2], orientation="vertical", shrink=0.5)
cbar2.set_label("Attributable relative affected population [# people]", fontsize = 9)
cbar2.ax.tick_params(labelsize=8)

axes[0].set_title("Factual population \n(2020)", fontsize=9)
axes[1].set_title("Counterfactual population \n(1990)", fontsize=9)
axes[2].set_title("Factual - Counterfactual", fontsize=9)

#%%
print("Plotting spatially RELATIVE aggregated exposed population damage for F, CF climate and diff")

gdf_pop_2020_exposed_CF_coarse['population_diff'] = (gdf_pop_2020_exposed_F_coarse['relative_population'] - gdf_pop_2020_exposed_CF_coarse['relative_population'])

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 5), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Create colormap normalization
norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2020_exposed_F_coarse["relative_population"].max())
norm_diff = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_1990_exposed_F_coarse["population_diff"].max())
red_half = LinearSegmentedColormap.from_list("bwr_red", plt.cm.bwr(np.linspace(0.5, 1, 256)))

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
gdf_pop_2020_exposed_F_coarse[gdf_pop_2020_exposed_F_coarse['relative_population'] == 0].plot(ax=axes[0], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2020_exposed_F_coarse[gdf_pop_2020_exposed_F_coarse['relative_population'] > 0].plot(column='relative_population', cmap='Blues',  edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[0], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_2020_exposed_CF_coarse[gdf_pop_2020_exposed_CF_coarse['relative_population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2020_exposed_CF_coarse[gdf_pop_2020_exposed_CF_coarse['relative_population'] > 0].plot(column='relative_population', cmap='Blues', edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[1], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_2020_exposed_CF_coarse[gdf_pop_2020_exposed_CF_coarse['population_diff'] <= 0].plot(ax=axes[2], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2020_exposed_CF_coarse[gdf_pop_2020_exposed_CF_coarse['population_diff'] > 0].plot(column='population_diff', cmap=red_half, norm=norm_diff, edgecolor='grey', 
                                                                 linewidth=0.2, ax=axes[2], legend=False, zorder=2, missing_kwds={"color": "white", "edgecolor": "none"}, rasterized=True)
subplot_labels = ['(a)', '(b)', '(c)']

for i, ax in enumerate(axes):
    # Add model region
    region_utm.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)

    # Add background and set extent (based on actual lat/lon coordinates)
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    minx, miny, maxx, maxy = region_utm.bounds.minx.item(), region_utm.bounds.miny.item(), region_utm.bounds.maxx.item(), region_utm.bounds.maxy.item()
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
sm1 = ScalarMappable(cmap="Blues", norm=norm_pop)
sm1._A = []  # required for colorbar without passing data
cbar = fig.colorbar(sm1, ax=axes[1], orientation="vertical", shrink=0.5)
cbar.set_label("Aggregated relativeexposed population (# people)", fontsize=10)
cbar.ax.tick_params(labelsize=9)

sm2 = ScalarMappable(cmap=red_half, norm=norm_diff)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, ax=axes[2], orientation="vertical", shrink=0.5)
cbar2.set_label("Attributable relative affected population [# people]", fontsize = 9)
cbar2.ax.tick_params(labelsize=8)

axes[0].set_title("Factual climate", fontsize=9)
axes[1].set_title("Counterfactual climate", fontsize=9)
axes[2].set_title("Factual - Counterfactual", fontsize=9)



# %%
###################################################################################
############################### Flood fatalities ##################################
###################################################################################
# Eq. 3 from Jonkman et al (2008), based on Boyd (2005)

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

# Calculate fatalities according to Jonkman et al (2008)
def fatality_fraction_Jonkman_etal2008(h, mu_N, sigma_N):
    """Compute fatality fraction FD(h) from water depth (m) using Jonkman (2008) lognormal fit."""
    from scipy.stats import norm

    h = np.maximum(h, 1e-3)  # avoid log(0)
    return norm.cdf((np.log(h) - mu_N) / sigma_N)


#%%
# Fatality function from Boyd (2005) as reported in Jonkman et al. (2008)
def fatality_curve_Boyd2005(h):
    """Fatality curve based on Boyd (2005) as reported in Jonkman et al. (2008).
    
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
# Apply fatality function of Boyd (2005)
df_pop_2020_exposed_F['fatalities'] = df_pop_2020_exposed_F['population'] * fatality_curve_Boyd2005(df_pop_2020_exposed_F['flood_depth'])
df_pop_2020_exposed_CF['fatalities'] = df_pop_2020_exposed_CF['population'] * fatality_curve_Boyd2005(df_pop_2020_exposed_CF['flood_depth'])
df_pop_1990_exposed_F['fatalities'] = df_pop_1990_exposed_F['population'] * fatality_curve_Boyd2005(df_pop_1990_exposed_F['flood_depth'])

print(f"Total fatalities acc. to Boyd et al. (2005):")
print(f"Factual climate 2020 population scenario:        {df_pop_2020_exposed_F['fatalities'].sum():.0f} people")
print(f"Counterfactual climate 2020 population scenario: {df_pop_2020_exposed_CF['fatalities'].sum():.0f} people")
print(f"Factual climate 1990 population scenario:        {df_pop_1990_exposed_F['fatalities'].sum():.0f} people")

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
model_velocity = {
    "model_name": model_name,
    "model_path": model_path,
    "sfincs_model": model_obj,
    "sfincs_results": model_obj.results,
    "sfincs_his": ds_his
    }

model_velocity = compute_hmax_masked_riserate([model_velocity], gwso)

hvmax = model_velocity[0]['sfincs_results']['hvmax_masked']
vmax = model_velocity[0]['sfincs_results']['vmax_masked']
hmax = model_velocity[0]['sfincs_results']['hmax_masked']
rise = model_velocity[0]['sfincs_results']['max_rise_rate_h']

# Apply Jonkman et al. (2008) fatality function based on zones defined by hvmax, vmax, hmax, and rise rate
# Initialize with zeros
fatality_frac = xr.zeros_like(hmax)

# Define masks
mask_breach = (hvmax >= 7) & (vmax > 2)
mask_rapid = (hmax > 2.1) & (rise > 0.5)
mask_remaining = (hmax > 0) & ~mask_breach & ~mask_rapid

# Apply per zone
fatality_frac = xr.where(mask_breach, 1, fatality_frac)
fatality_frac = xr.where(mask_rapid, fatality_fraction_Jonkman_etal2008(hmax, mu_N=1.46, sigma_N=0.28), fatality_frac)
fatality_frac = xr.where(mask_remaining, fatality_fraction_Jonkman_etal2008(hmax, mu_N=7.60, sigma_N=2.75), fatality_frac)

# Store the result
model_velocity[0]['sfincs_results']['fatality_frac'] = fatality_frac


# Calculate exposed fatalities for exposed population to factual flooding
flood_fatalities_F = model_velocity[0]['sfincs_results']['fatality_frac'] * ra_exposed_pop_2020_F

print(f"Total fatalities acc. to Jonkman et al. (2008): {flood_fatalities_F.sum().values:.0f} people")
# %%

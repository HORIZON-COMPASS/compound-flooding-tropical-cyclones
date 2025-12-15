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
from affine import Affine
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
from matplotlib.colors import PowerNorm, TwoSlopeNorm
from rasterio.transform import rowcol



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
background = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_region_background.geojson", driver="GeoJSON")
region = gpd.read_file(shapefile_fp)
shapefile_sofala = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_province.shp")

# Load the admin3 district in the case study region to validate exposed people
districts_adm3 = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_districts_study_region.shp")
districts_adm2 = data_catalog.get_geodataframe("gadm_level2", geom=region, buffer=1000)

# population in provided inthousand persons per grid cell
population_raster_path_2019 = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_2019_30.tif")  
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
region = region.to_crs(flood_grid_crs)
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
districts_adm3_utm = districts_adm3.to_crs(flood_grid_crs)
districts_adm2_utm = districts_adm2.to_crs(flood_grid_crs)

# Remove districts that are not connecting to the region
drop_districts = ["Muanza", "Gororngosa-Sede", "Galinha"]
districts_adm3_filtered = districts_adm3_utm[~districts_adm3_utm['NAME_3'].isin(drop_districts)]


#%% Read and regrid population data
# --- Function to redistribute population over land pixels on flood grid ---
def reproject_and_redistribute_population_over_land(pop_path, land_gdf, flood_crs, flood_transform, flood_shape, province_geom=None, region=None, districts_adm3=None, districts_adm2=None, year=None, out_raster_path=None):    
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
        districts_adm3_single = districts_adm3.dissolve().reset_index(drop=True)
        districts_adm3_single = districts_adm3_single.to_crs(src.crs)
        districts_adm3_geom   = [districts_adm3_single.geometry.iloc[0].__geo_interface__]

        districts_adm2_single = districts_adm2.dissolve().reset_index(drop=True)
        districts_adm2_single = districts_adm2_single.to_crs(src.crs)
        districts_adm2_geom   = [districts_adm2_single.geometry.iloc[0].__geo_interface__]

        pop_districts_adm3, pop_affine_districts_adm3 = mask(src, districts_adm3_geom, crop=True, nodata=src.nodata)
        pop_districts_adm2, pop_affine_districts_adm2 = mask(src, districts_adm2_geom, crop=True, nodata=src.nodata)

        if out_raster_path is not None and os.path.exists(out_raster_path):
            print(f"▶ Loading existing raster from {out_raster_path}")
            with rasterio.open(out_raster_path) as src:
                pop_fine = src.read(1)
                return pop_fine, pop_sofala, transform_sofala, pop_districts_adm3, pop_affine_districts_adm3, pop_districts_adm2, pop_affine_districts_adm2

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
                # pop_fine[valid_mask] += pop_value / n_valid # old resulting in "people behind the comma"

                # --- Integer redistribution (whole people only) ---
                P = int(round(pop_value))
                N = n_valid

                base = P // N
                remainder = P % N

                # convert mask indices to linear index list
                valid_indices = np.where(valid_mask)

                # assign the base number to all pixels
                pop_fine[valid_indices] += base

                # assign remainder one-by-one to the first 'remainder' pixels
                # (or shuffle if you prefer randomness)
                if remainder > 0:
                    # deterministic: first R indices
                    # if random desired: uncomment the shuffle line below
                    perm = np.random.permutation(len(valid_indices[0]))
                    chosen = perm[:remainder]
                    # chosen = np.arange(remainder)
                    pop_fine[valid_indices[0][chosen], valid_indices[1][chosen]] += 1


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
        H, W = pop_fine.shape
        a, b, c, d, e, f = flood_transform  # unpack affine

        if e > 0:
            print("  ⚠ Detected positive y-resolution in transform → fixing for QGIS")

            # 1) flip array vertically
            pop_fine_to_write = np.flipud(pop_fine)

            # 2) fix y-scale sign and y-origin
            new_e = -abs(e)
            new_f = f + e * (H - 1)

            fixed_transform = Affine(a, b, c, d, new_e, new_f)

        else:
            pop_fine_to_write = pop_fine
            fixed_transform = flood_transform

        # ------------------------------------------------------------------

        new_profile = {
            "driver": "GTiff",
            "dtype": rasterio.float32,
            "count": 1,
            "height": H,
            "width": W,
            "crs": flood_crs,
            "transform": fixed_transform,
            "compress": "deflate"
        }

        with rasterio.open(out_raster_path, "w", **new_profile) as dst:
            dst.write(pop_fine_to_write, 1)

    return pop_fine, pop_sofala, transform_sofala, pop_districts_adm3, pop_affine_districts_adm3, pop_districts_adm2, pop_affine_districts_adm2

# --- Function to link population raster to flood depth as DataFrame ---
def pop_raster_to_gdf(pop_array, flood_array, transform, year, climate, export_df=True, export_path=None):
    print("Linking population raster to flood depth as DataFrame...")

    # Check if shapes match
    if pop_array.shape == flood_array.shape:
        print("✔ Shapes match")
    else:
        print("✖ Shapes do NOT match!", pop_array.shape, flood_array.shape)

    # Flatten arrays
    pop_flat = pop_array.ravel()
    flood_flat = flood_array.ravel()

    # Mask zero-pop cells
    mask = (pop_flat > 0) & (flood_flat > 0)
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

    if export_df:
        file_name = f"df_pop_{year}_{climate}.csv"
        df.to_csv(join(export_path, file_name), index=False)
        print(f"▶ Exported DataFrame to {join(export_path, file_name)}")

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

# --- Function to save raster only if it does not yet exist ---
def save_raster(array, out_path, transform, crs):
    # """Save raster only if it does NOT yet exist."""
    # if os.path.exists(out_path):
    #     print(f"✔ File already exists, skipping: {out_path}")
    #     return
    
    print(f"▶ Writing raster: {out_path}")
    profile = {
        "driver": "GTiff",
        "dtype": rasterio.float32,
        "count": 1,
        "height": array.shape[0],
        "width": array.shape[1],
        "crs": crs,
        "transform": transform,
        "compress": "deflate"
    }
    
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(array.astype("float32"), 1)



#%%
# ============================================================================================ #
# ====================== Process population directly into flood grid ========================= #
# ============================================================================================ #
export_path = "p:/11210471-001-compass/04_Results/Idai_socioeconomic/preprocessed/population/"

pop_arrays = {}
pop_sofala_arrays = {}
pop_sofala_districts_adm3 = {}
pop_affine_sofala_districts_adm3 = {}
pop_sofala_districts_adm2 = {}
pop_affine_sofala_districts_adm2 = {}
# --- Reproject to flood grid and redistribute population rasters over land ---
for year, path in [(1990, population_raster_path_1990),
                   (2019, population_raster_path_2019)]:
    pop_arrays[year], pop_sofala_arrays[year], transform_sofala_land, pop_sofala_districts_adm3[year], pop_affine_sofala_districts_adm3[year], pop_sofala_districts_adm2[year], pop_affine_sofala_districts_adm2[year] = reproject_and_redistribute_population_over_land(
        pop_path=path, land_gdf=background_utm, flood_crs=flood_grid_crs, flood_transform=flood_grid_transform,
        flood_shape=flood_grid_shape, province_geom=shapefile_sofala, region=region, districts_adm3=districts_adm3_filtered,
        districts_adm2=districts_adm2, year=year,
        out_raster_path=f"c:/Code/COMPASS_exposure/Data/Modified/population_{year}_region_regrid_wholepeople.tif"
    ) 

# --- Compute exposed population GeoDataFrames ---
gdf_pop_2019_flood_depth_F  = pop_raster_to_gdf(pop_arrays[2019], hmax_F, flood_grid_transform, year=2019, climate="F", export_df=True, export_path=export_path)
gdf_pop_2019_flood_depth_CF = pop_raster_to_gdf(pop_arrays[2019], hmax_CF, flood_grid_transform, year=2019, climate="CF", export_df=True, export_path=export_path)
gdf_pop_1990_flood_depth_F  = pop_raster_to_gdf(pop_arrays[1990], hmax_F, flood_grid_transform, year=1990, climate="F", export_df=True, export_path=export_path)
gdf_pop_1990_flood_depth_CF = pop_raster_to_gdf(pop_arrays[1990], hmax_CF, flood_grid_transform, year=1990, climate="CF", export_df=True, export_path=export_path)

# simple rasters for fast plotting
ra_exposed_pop_2019_F  = np.where(hmax_F > 0, pop_arrays[2019], 0)
ra_exposed_pop_2019_CF = np.where(hmax_CF > 0, pop_arrays[2019], 0)
ra_exposed_pop_1990_F  = np.where(hmax_F > 0, pop_arrays[1990], 0)
ra_exposed_pop_1990_CF = np.where(hmax_CF > 0, pop_arrays[1990], 0)

# Compute exposed pop on fine grid
gdf_pop_2019_exposed_F  = gdf_pop_2019_flood_depth_F[gdf_pop_2019_flood_depth_F['flood_depth'] > 0]
gdf_pop_2019_exposed_CF = gdf_pop_2019_flood_depth_CF[gdf_pop_2019_flood_depth_CF['flood_depth'] > 0]
gdf_pop_1990_exposed_F  = gdf_pop_1990_flood_depth_F[gdf_pop_1990_flood_depth_F['flood_depth'] > 0]
gdf_pop_1990_exposed_CF = gdf_pop_1990_flood_depth_CF[gdf_pop_1990_flood_depth_CF['flood_depth'] > 0]

# --- Aggregate population to coarser grid ---
gdf_pop_2019_exposed_F_coarse  = aggregate_pop(pop_arrays[2019], hmax_F, flood_grid_transform, flood_grid_crs, region_utm, background_utm)
gdf_pop_2019_exposed_CF_coarse = aggregate_pop(pop_arrays[2019], hmax_CF, flood_grid_transform, flood_grid_crs, region_utm, background_utm)
gdf_pop_1990_exposed_F_coarse  = aggregate_pop(pop_arrays[1990], hmax_F, flood_grid_transform, flood_grid_crs, region_utm, background_utm)
gdf_pop_1990_exposed_CF_coarse = aggregate_pop(pop_arrays[1990], hmax_CF, flood_grid_transform, flood_grid_crs, region_utm, background_utm)

#%%
# save raster for fatality analysis in other script
# save_raster(ra_exposed_pop_2019_F,  "p:/11210471-001-compass/04_Results/Idai_socioeconomic/preprocessed/population/exposed_pop_2019_F.tif",  flood_grid_transform, flood_grid_crs)
# save_raster(ra_exposed_pop_2019_CF, "p:/11210471-001-compass/04_Results/Idai_socioeconomic/preprocessed/population/exposed_pop_2019_CF.tif", flood_grid_transform, flood_grid_crs)
# save_raster(ra_exposed_pop_1990_F,  "p:/11210471-001-compass/04_Results/Idai_socioeconomic/preprocessed/population/exposed_pop_1990_F.tif",  flood_grid_transform, flood_grid_crs)
# save_raster(ra_exposed_pop_1990_CF, "p:/11210471-001-compass/04_Results/Idai_socioeconomic/preprocessed/population/exposed_pop_1990_CF.tif", flood_grid_transform, flood_grid_crs)

#%%
# Uniform population growth
pop_growth = np.nansum(pop_arrays[2019]) / np.nansum(pop_arrays[1990])
print(f"Uniform population growth from 1990 to 2019 in Sofala region: {(pop_growth*100):.2f}%")

pop_array_uniform_2019 = pop_arrays[1990] * (pop_growth)

# Sanity check
print(f"{np.nansum(pop_array_uniform_2019):,.0f} people in 2019 with uniform growth")
print(f"{np.nansum(pop_arrays[2019]):,.0f} people in 2019 actual")

# Calculate exposure
ra_exposed_pop_2019_F_uniform = np.where(hmax_F > 0, pop_array_uniform_2019, 0)
ra_exposed_pop_2019_CF_uniform = np.where(hmax_CF > 0, pop_array_uniform_2019, 0)

# get flood depth per affected population
gdf_pop_2019_flood_depth_F_uniform  = pop_raster_to_gdf(pop_array_uniform_2019, hmax_F, flood_grid_transform, year='2019_uniform', climate="F", export_df=True, export_path=export_path)
gdf_pop_2019_flood_depth_CF_uniform = pop_raster_to_gdf(pop_array_uniform_2019, hmax_CF, flood_grid_transform, year='2019_uniform', climate="CF", export_df=True, export_path=export_path)

gdf_pop_2019_exposed_F_uniform  = gdf_pop_2019_flood_depth_F_uniform[gdf_pop_2019_flood_depth_F_uniform['flood_depth'] > 0]
gdf_pop_2019_exposed_CF_uniform = gdf_pop_2019_flood_depth_CF_uniform[gdf_pop_2019_flood_depth_CF_uniform['flood_depth'] > 0]

# Aggregate to coarser cells
gdf_pop_2019_exposed_F_uniform_coarse = aggregate_pop(pop_array_uniform_2019, hmax_F, flood_grid_transform, flood_grid_crs, region_utm, background_utm)

# save_raster(ra_exposed_pop_2019_F_uniform,  "p:/11210471-001-compass/04_Results/Idai_socioeconomic/preprocessed/population/exposed_pop_2019_F_uniform.tif",  flood_grid_transform, flood_grid_crs)
# save_raster(ra_exposed_pop_2019_CF_uniform, "p:/11210471-001-compass/04_Results/Idai_socioeconomic/preprocessed/population/exposed_pop_2019_CF_uniform.tif", flood_grid_transform, flood_grid_crs)
# save_raster(pop_array_uniform_2019,         "p:/11210471-001-compass/04_Results/Idai_socioeconomic/preprocessed/population/population_2019_uniform.tif",      flood_grid_transform, flood_grid_crs)

#%% ============================================================================================= #
# ===================== Print summary statistics of exposed population ========================== #
# =============================================================================================== #
print("2019 Factual exposed population stats:")
print("Total population in region:", np.nansum(pop_arrays[2019]))
print("Total population in Sofala:", np.nansum(pop_sofala_arrays[2019]))
print("Total exposed people:", np.nansum(ra_exposed_pop_2019_F).astype(int))
print("Exposed people percentage of total population:", 100 * np.nansum(ra_exposed_pop_2019_F).astype(int) / np.nansum(pop_arrays[2019]))

print("\n2019 Counterfactual exposed population stats:")
print("Total population in region:", np.nansum(pop_arrays[2019]))
print("Total population in Sofala:", np.nansum(pop_sofala_arrays[2019]))
print("Total exposed people:", np.nansum(ra_exposed_pop_2019_CF).astype(int))
print("Exposed people percentage of total population:", 100 * np.nansum(ra_exposed_pop_2019_CF).astype(int) / np.nansum(pop_arrays[2019]))

print("\n1990 Factual exposed population stats:")
print("Total population in region:", np.nansum(pop_arrays[1990]))
print("Total population in Sofala:", np.nansum(pop_sofala_arrays[1990]))
print("Total exposed people:", np.nansum(ra_exposed_pop_1990_F).astype(int))
print("Exposed people percentage of total population:", 100 * np.nansum(ra_exposed_pop_1990_F).astype(int) / np.nansum(pop_arrays[1990]))

print("\n1990 Counterfactual exposed population stats:")
print("Total population in region:", np.nansum(pop_arrays[1990]))
print("Total population in Sofala:", np.nansum(pop_sofala_arrays[1990]))
print("Total exposed people:", np.nansum(ra_exposed_pop_1990_CF).astype(int))
print("Exposed people percentage of total population:", 100 * np.nansum(ra_exposed_pop_1990_CF).astype(int) / np.nansum(pop_arrays[1990]))

print("\n2019 UNIFORM exposed population stats:")
print("Total population in region:", np.nansum(pop_array_uniform_2019))
print("Total exposed people:", np.nansum(ra_exposed_pop_2019_F_uniform).astype(int))
print("Exposed people percentage of total population:", 100 * np.nansum(ra_exposed_pop_2019_F_uniform).astype(int) / np.nansum(pop_arrays[2019]))

print("\nOne-line attribution numbers:")
print(f"Exposed population in 2019 Factual: {int(np.nansum(ra_exposed_pop_2019_F).astype(int)):,}")
print(f"Exposed population attributable to climate change: {int(np.nansum(ra_exposed_pop_2019_F) - np.nansum(ra_exposed_pop_2019_CF)):,} {100 * (np.nansum(ra_exposed_pop_2019_F) - np.nansum(ra_exposed_pop_2019_CF)) / np.nansum(ra_exposed_pop_2019_F):.2f}%")
print(f"Exposed population attributable to population change (2019-1990): {int(np.nansum(ra_exposed_pop_2019_F) - np.nansum(ra_exposed_pop_1990_F)):,} {100 * (np.nansum(ra_exposed_pop_2019_F) - np.nansum(ra_exposed_pop_1990_F)) / np.nansum(ra_exposed_pop_2019_F):.2f}%")
print(f"Exposed population attributable to population change and climate change: {int(np.nansum(ra_exposed_pop_2019_F) - np.nansum(ra_exposed_pop_1990_CF)):,} {100 * (np.nansum(ra_exposed_pop_2019_F) - np.nansum(ra_exposed_pop_1990_CF)) / np.nansum(ra_exposed_pop_2019_F):.2f}%")

print(f"Exposed population attributable to population change (uniform growth): {int(np.nansum(ra_exposed_pop_2019_F) - np.nansum(ra_exposed_pop_2019_F_uniform)):,} {100 * (np.nansum(ra_exposed_pop_2019_F) - np.nansum(ra_exposed_pop_2019_F_uniform)) / np.nansum(ra_exposed_pop_2019_F):.2f}%")
print(f"Population growth from 1990 to 2019 in the region: {int(np.nansum(pop_arrays[2019]) - np.nansum(pop_arrays[1990])):,} {100 * (np.nansum(pop_arrays[2019]) - np.nansum(pop_arrays[1990])) / np.nansum(pop_arrays[1990]):.2f}%")


# %%
# # ============================================================================================ #
# # ================== Plot comparison rasterized vs vectorized population ===================== #
# # ============================================================================================ #

# gdf_pop_2019_exposed_F_coarse_vector = aggregate_pop_sjoin_slow(pop_arrays[2019], hmax_F, flood_grid_transform, flood_grid_crs, region_utm, background_utm, factor=100)

# gdf_2019_F_vector = gdf_pop_2019_exposed_F_coarse_vector.copy()
# gdf_2019_F_vector.loc[gdf_2019_F_vector["total_population"] == 0] = np.nan

# gdf_2019_F = gdf_pop_2019_exposed_F_coarse.copy()
# gdf_2019_F.loc[gdf_2019_F["total_population"] == 0] = np.nan

# # Ensure vector and raster GDFs have only the necessary columns
# gdf_vector = gdf_2019_F_vector[['geometry', 'total_population']].copy()
# gdf_raster = gdf_2019_F[['geometry', 'total_population']].copy()

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

# gdf_2019_F_vector = gdf_pop_2019_exposed_F_coarse_vector.copy()
# gdf_2019_F_vector.loc[gdf_2019_F_vector["total_population"] == 0] = np.nan

# gdf_2019_F = gdf_pop_2019_exposed_F_coarse.copy()
# gdf_2019_F.loc[gdf_2019_F["total_population"] == 0] = np.nan

# gdf_2019_F['change_in_population_vect'] = gdf_2019_F_vector['total_population'] - gdf_2019_F['total_population']

# # Define colormap: from white to #67CBE4
# cmap = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#FC6F37"])
# norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_pop_2019_exposed_F_coarse_vector['total_population']))  
# cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#BD2A2A"])
# norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_compare['change_in_population']))

# plot = gdf_2019_F_vector.plot(column="total_population", cmap=cmap, norm=norm ,
#                        linewidth=0.1, edgecolor="grey", ax=axes[0], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

# plot = gdf_2019_F.plot(column="total_population", cmap=cmap, norm=norm,
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

# axes[0].set_title("Factual 2019 vectorized", fontsize=10)
# axes[1].set_title("Factual 2019 rasterized", fontsize=10)
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
flood_F_2019_sorted, cdf_F, pop_2019_by_depth_F   = compute_cdf_and_bins(gdf_pop_2019_exposed_F, bins)
flood_CF_2019_sorted, cdf_CF_clim, pop_2019_by_depth_CF = compute_cdf_and_bins(gdf_pop_2019_exposed_CF, bins)
flood_F_1990_sorted, cdf_CF_pop, pop_1990_by_depth_F = compute_cdf_and_bins(gdf_pop_1990_exposed_F, bins)
flood_CF_1990_sorted, cdf_CF_pop_clim, pop_1990_by_depth_CF = compute_cdf_and_bins(gdf_pop_1990_exposed_CF, bins)

bins_fine = np.arange(0, 3.5 + 0.1, 0.1)
flood_F_2019_sorted_fine, cdf_F_fine, pop_2019_by_depth_F_fine   = compute_cdf_and_bins(gdf_pop_2019_exposed_F, bins_fine)
flood_CF_2019_sorted_fine, cdf_CF_clim_fine, pop_2019_by_depth_CF_fine = compute_cdf_and_bins(gdf_pop_2019_exposed_CF, bins_fine)
flood_F_1990_sorted_fine, cdf_CF_pop_fine, pop_1990_by_depth_F_fine = compute_cdf_and_bins(gdf_pop_1990_exposed_F, bins_fine)
flood_CF_1990_sorted_fine, cdf_CF_pop_clim_fine, pop_1990_by_depth_CF_fine = compute_cdf_and_bins(gdf_pop_1990_exposed_CF, bins_fine)

# --- Plot continuous CDF ---
plt.figure(figsize=(8,5))
plt.plot(flood_F_2019_sorted, cdf_F, label="Factual")
plt.plot(flood_CF_2019_sorted, cdf_CF_clim, label="Counterfactual Climate")
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
plt.bar(bin_centers - 0.05, pop_2019_by_depth_F.values, width=0.05, label="Factual", alpha=0.7)
plt.bar(bin_centers, pop_2019_by_depth_CF.values, width=0.05, label="Counterfactual Climate", alpha=0.7)
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
plt.plot(bin_centers, pop_2019_by_depth_F_fine.values, label="Factual")
plt.plot(bin_centers, pop_2019_by_depth_CF_fine.values, label="Counterfactual Climate")
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

kde = gaussian_kde(gdf_pop_2019_exposed_F['flood_depth'].values, weights=gdf_pop_2019_exposed_F['population'].values)
y_F = kde(x) # probability density
# Find peaks and sort by height
peaks, _ = find_peaks(y_F)
sorted_peaks_F = peaks[np.argsort(y_F[peaks])[::-1]]

kde = gaussian_kde(gdf_pop_2019_exposed_CF['flood_depth'].values, weights=gdf_pop_2019_exposed_CF['population'].values)
y_CF_clim = kde(x) # probability density
# Find peaks and sort by height
peaks, _ = find_peaks(y_CF_clim)
sorted_peaks_CF_clim = peaks[np.argsort(y_CF_clim[peaks])[::-1]]

kde = gaussian_kde(gdf_pop_1990_exposed_F['flood_depth'].values, weights=gdf_pop_1990_exposed_F['population'].values)
y_CF_pop = kde(x) # probability density
# Find peaks and sort by height
peaks, _ = find_peaks(y_CF_pop)
sorted_peaks_CF_pop = peaks[np.argsort(y_CF_pop[peaks])[::-1]]

kde = gaussian_kde(gdf_pop_1990_exposed_CF['flood_depth'].values, weights=gdf_pop_1990_exposed_CF['population'].values)
y_CF_clim_pop = kde(x) 

kde = gaussian_kde(gdf_pop_2019_exposed_F_uniform['flood_depth'].values, weights=gdf_pop_2019_exposed_F_uniform['population'].values)
y_F_uni = kde(x) # probability density
# Find peaks and sort by height
peaks, _ = find_peaks(y_F_uni)
sorted_peaks_F_uni = peaks[np.argsort(y_F_uni[peaks])[::-1]]

#%%
# Plot density-weighted population exposure per flood depth
# with arrows indicating peak differences
def match_peaks_by_x(x, peaksA, peaksB):
    """
    Match peaks between two scenarios by closest x position.
    Returns list of tuples: (peakA_idx, peakB_idx)
    where either can be None if unmatched.
    """
    xsA = x[peaksA]
    xsB = x[peaksB]

    matched = []
    used_B = set()

    # Match each A peak to nearest B peak
    for iA, xA in zip(peaksA, xsA):
        if len(xsB) > 0:
            # find closest B peak not yet used
            diffs = [(abs(xA - xB), iB) for xB, iB in zip(xsB, peaksB) if iB not in used_B]
            if diffs:
                _, bestB = min(diffs, key=lambda t: t[0])
                used_B.add(bestB)
                matched.append((iA, bestB))
            else:
                matched.append((iA, None))
        else:
            matched.append((iA, None))

    # Add B peaks with no match in A
    for iB in peaksB:
        if iB not in used_B:
            matched.append((None, iB))

    return matched

# Draw arrows between matched peaks
def draw_peak_arrows(x, yA, yB, peaksA, peaksB, color, label_prefix, label_offsets, arrow_offsets,
                     min_diff=0.01):

    pairs = match_peaks_by_x(x, peaksA, peaksB)

    for i, (iA, iB) in enumerate(pairs):

        # CASE 1: Peak exists in both lines --------------------------------------
        if iA is not None and iB is not None:
            xA, yA_val = x[iA], yA[iA]
            xB, yB_val = x[iB], yB[iB]

        # CASE 2: Peak in A but not in B -----------------------------------------
        elif iA is not None and iB is None:
            xA = x[iA]
            yA_val = yA[iA]
            xB = xA  # vertical arrow
            yB_val = yB[iA]  # value on other line at same x

        # CASE 3: Peak in B but not in A -----------------------------------------
        elif iB is not None and iA is None:
            xB = x[iB]
            yB_val = yB[iB]
            xA = xB  # vertical arrow
            yA_val = yA[iB]  # value on other line at same x

        # compute difference
        diff = abs(yA_val - yB_val)

        # FILTER 1 — small differences
        if diff < min_diff:
            continue

        # vertical arrow if xA == xB (peak ↔ no peak)
        offset_x_arrow, offset_y_arrow = arrow_offsets.get(i, (0.0, 0.0))
        plt.annotate(
            "",
            xy=(xB+offset_x_arrow, yB_val+offset_y_arrow),
            xytext=(xA, yA_val),
            arrowprops=dict(arrowstyle="->", lw=1.4, color=color))

        # add label
        offset_x, offset_y = label_offsets.get(i, (0.03, 0.0)) 
        xm = xA + offset_x  # midpoint in x (same if vertical)
        ym = (yA_val + yB_val)/2 + offset_y
        plt.text(xm, ym, f"{label_prefix}",
                 ha="left", va="center", fontsize=9, color=color, fontweight="bold",
                 bbox=dict(boxstyle="round",pad=0.15, facecolor="lightgrey", 
                           edgecolor="none", alpha=0.5))

# Colours based on conceptual figure
colours = ['#00B050', '#1E2E57', "#28C2E9", '#9B59B6']

fig, (ax_main, ax_zoom) = plt.subplots(2, 1, figsize=(9, 8), gridspec_kw={'height_ratios': [3, 1]}, dpi=300)

plt.sca(ax_main)

# Define the boundaries
x_bg = np.linspace(0, 3.5, 500)  # example x array

low_mask = x_bg < 0.5                
mid_mask = (x_bg >= 0.5) & (x_bg < 1.5)  
high_mask = x_bg >= 1.5                

# # Fill each region from 0 to the curve
ax_main.fill_between(x_bg[low_mask], 0, max(y_CF_clim_pop)*1.05, color="#d9d9d9", alpha=0.3, step='post')
ax_main.fill_between(x_bg[mid_mask], 0, max(y_CF_clim_pop)*1.05, color="#b3b3b3", alpha=0.3, step='post')
ax_main.fill_between(x_bg[high_mask], 0, max(y_CF_clim_pop)*1.05, color="#808080", alpha=0.3, step='post')

ax_main.plot(x, y_F, label=f"Factual ({np.nansum(pop_arrays[2019]):,.0f} people)", color=colours[0], linewidth=2)
ax_main.plot(x, y_CF_clim, label=f"Counterfactual Climate ({np.nansum(pop_arrays[2019]):,.0f} people)", color=colours[1], linewidth=1)
ax_main.plot(x, y_CF_pop, label=f"Counterfactual Population ({np.nansum(pop_arrays[1990]):,.0f} people)", color=colours[2], linewidth=1)
ax_main.plot(x, y_CF_clim_pop, label=f"Counterfactual Climate & Population ({np.nansum(pop_arrays[1990]):,.0f} people)", color=colours[3], linewidth=1)

ax_main.set_ylabel("Density-weighted affected population")
ax_main.grid(True, linestyle="--", alpha=0.5)
ax_main.set_xlim(0, 3.5)

# peak annotations 
if len(sorted_peaks_F) > 1:
    second_peak_x = x[sorted_peaks_F[1]]
    plt.annotate(f"{x[sorted_peaks_F[0]]:.2f} m",
                 xy=(x[sorted_peaks_F[0]]+0.02, y_F[sorted_peaks_F[0]]),
                 xytext=(x[sorted_peaks_F[0]]+0.15, y_F[sorted_peaks_F[0]]-0.05), color=colours[0], fontsize=8,
                 arrowprops=dict(arrowstyle="->", lw=0.8, color=colours[0], linestyle='--'))
    plt.annotate(f"{x[sorted_peaks_F[1]]:.2f} m", 
                 xy=(x[sorted_peaks_F[1]]+0.02, y_F[sorted_peaks_F[1]]),
                 xytext=(x[sorted_peaks_F[1]]+0.15, y_F[sorted_peaks_F[1]]), color=colours[0], fontsize=8,
                 arrowprops=dict(arrowstyle="->", lw=0.8, color=colours[0], linestyle='--'))

if len(sorted_peaks_CF_clim) > 1:
    second_peak_x = x[sorted_peaks_CF_clim[1]]
for i, idx in enumerate(sorted_peaks_CF_clim[:2]): 
    plt.annotate(f"{x[idx]:.2f} m",
                 xy=(x[idx]+0.02, y_CF_clim[idx]),
                 xytext=(x[idx]+0.15, y_CF_clim[idx]+0.05), color=colours[1], fontsize=8,
                 arrowprops=dict(arrowstyle="->", lw=0.8, color=colours[1], linestyle='--'))

for i, idx in enumerate(sorted_peaks_CF_pop[:2]): 
    plt.annotate(f"{x[idx]:.2f} m",
                 xy=(x[idx]+0.01, y_CF_pop[idx]),
                 xytext=(x[idx]+0.15, y_CF_pop[idx]+0.03), color=colours[2], fontsize=8,
                 arrowprops=dict(arrowstyle="->", lw=0.8, color=colours[2], linestyle='--'))

# draw arrows
draw_peak_arrows(x, y_CF_clim, y_F, sorted_peaks_CF_clim, sorted_peaks_F, 
                 label_offsets = {0: (0.1, 0.05), 1: (0.17, 0.03)},
                 arrow_offsets = {0: (0.0, 0.0), 1: (0.0, 0.00)},
                 color=colours[0], label_prefix="Climate change")
draw_peak_arrows(x, y_CF_pop, y_F, sorted_peaks_CF_pop, sorted_peaks_F,
                 label_offsets = {0: (0.1, 0.002), 1: (0.13, -0.01)},
                 arrow_offsets = {0: (0.0, -0.03), 1: (0.0, 0.00)},
                 color=colours[0], label_prefix="Population change")

leg = ax_main.legend(loc='upper right', fontsize=9)
leg.get_frame().set_facecolor('white')  # background color
leg.get_frame().set_alpha(0.4)            # transparency

ax_main.fill_between(x, y_F, y_CF_clim_pop, where=((y_CF_clim_pop > y_CF_clim) & (x < 0.5)),
                     alpha=0.15, color=colours[3])
ax_main.fill_between(x, y_CF_pop, y_F, where=((y_CF_pop > y_F)),
                     alpha=0.15, color=colours[2])
ax_main.fill_between(x, y_F, y_CF_clim, where=(y_CF_clim > y_F),
                     alpha=0.15, color=colours[1])
ax_main.fill_between(x, y_CF_pop, y_F, where=(y_F > y_CF_pop),
                     alpha=0.15, color=colours[0])

ax_main.text(0.25, max(y_F)*0.07, "Low", ha="center", va="center", fontsize=9, color="#5C5C5C", fontweight="bold")
ax_main.text(1.0, max(y_F)*0.07, "Medium", ha="center", va="center", fontsize=9, color="#5C5C5C", fontweight="bold")
ax_main.text(2.5, max(y_F)*0.07, "High", ha="center", va="center", fontsize=9, color="#5C5C5C", fontweight="bold")

ax_main.set_ylim(0, y_CF_clim_pop.max() * 1.05)
# -----------------
# ZOOMED PLOT
# -----------------
plt.sca(ax_zoom)

ax_zoom.fill_between(x_bg[high_mask], 0, max(y_CF_clim_pop)*1.05, color="#808080", alpha=0.3, step='post')
ax_zoom.text(2.5, max(y_F)*0.05, "High", ha="center", va="center", fontsize=9, color="#5C5C5C", fontweight="bold")

ax_zoom.plot(x, y_F, color=colours[0], linewidth=2)
ax_zoom.plot(x, y_CF_clim, color=colours[1], linewidth=1)
ax_zoom.plot(x, y_CF_pop, color=colours[2], linewidth=1)
ax_zoom.plot(x, y_CF_clim_pop, color=colours[3], linewidth=1)

ax_zoom.set_xlim(1.5, x.max())
ax_zoom.set_ylim(0, 0.15)
ax_zoom.set_xlabel("Flood depth (m)")
ax_zoom.set_ylabel("")
ax_zoom.set_title("Zoomed-in to tail")
ax_zoom.grid(True, linestyle="--", alpha=0.5)


ax_zoom.fill_between(x, y_CF_pop, y_F, where=(y_CF_pop > y_F),
                     alpha=0.15, color=colours[2])
ax_zoom.fill_between(x, y_F, y_CF_clim_pop, where=(y_CF_clim_pop > y_F),
                     alpha=0.15, color=colours[3])
ax_zoom.fill_between(x, y_CF_clim, y_F, where=(y_F > y_CF_clim),
                     alpha=0.15, color=colours[0])
ax_zoom.fill_between(x, y_CF_clim_pop, y_F, where=(y_F > y_CF_clim_pop),
                     alpha=0.15, color=colours[0])


ax_zoom.annotate("", xy=(1.62, 0.10), xytext=(1.62, 0.057), ha="center", va="center", 
                 arrowprops=dict(arrowstyle="->", lw=1.2, color=colours[0], shrinkA=0, shrinkB=0))
ax_zoom.annotate("", xy=(2.67, 0.032), xytext=(2.65, 0.022), ha="center", va="center", 
                 arrowprops=dict(arrowstyle="->", lw=1.2, color=colours[0], shrinkA=0, shrinkB=0))
ax_zoom.annotate("", xy=(2.15, 0.073), xytext=(2.15, 0.05), ha="center", va="center", 
                 arrowprops=dict(arrowstyle="->", lw=1.2, color=colours[3], shrinkA=0, shrinkB=0))
ax_zoom.annotate("", xy=(1.66, 0.115), xytext=(1.65, 0.1), ha="center", va="center", 
                 arrowprops=dict(arrowstyle="->", lw=1.2, color=colours[2], shrinkA=0, shrinkB=0))

ax_zoom.text(1.67, 0.04, f"Climate change", ha="center", va="center", fontsize=8, color=colours[0], fontweight="bold",
             bbox=dict(boxstyle="round",pad=0.15, facecolor="lightgrey", edgecolor="none", alpha=0.5))
ax_zoom.text(2.67, 0.04, f"Climate change", ha="left", va="center", fontsize=8, color=colours[0], fontweight="bold",
             bbox=dict(boxstyle="round",pad=0.15, facecolor="lightgrey", edgecolor="none", alpha=0.5))
ax_zoom.text(2.15, 0.087, f"Relative population change", ha="center", va="center", fontsize=8, color=colours[3], fontweight="bold",
             bbox=dict(boxstyle="round",pad=0.15, facecolor="lightgrey", edgecolor="none", alpha=0.5))
ax_zoom.text(1.68, 0.125, f"Relative population change \n         + climate change", ha="left", va="center", fontsize=8, color=colours[2], 
             fontweight="bold", bbox=dict(boxstyle="round",pad=0.15, facecolor="lightgrey", edgecolor="none", alpha=0.5))

#%%
# Check whether population count varies per cell
mask = (y_F > y_CF_pop) & (x < 1.5)

# X-range of the filled area
x_min = x[mask].min()
x_max = x[mask].max()

mask = (gdf_pop_2019_exposed_F['flood_depth'] > x_min) & \
       (gdf_pop_2019_exposed_F['flood_depth'] < x_max)

pop_in_range = gdf_pop_2019_exposed_F.loc[mask, 'population']

print(pop_in_range.describe())

# Histogram of values
plt.hist(pop_in_range, bins=30, edgecolor='k', alpha=0.7)
plt.xlabel("Population in cell")
plt.ylabel("Number of cells")
plt.title(f"Population distribution for flood depth {x_min:.2f}-{x_max:.2f} m")
plt.show()

# %%
# -----------------
# SUMMARY BOX
# -----------------
mask = (x > 0) & (x < 0.5)
x_min_min = x[mask].min()
x_max_min = x[mask].max()

mask = (x > 0.5) & (x < 1.5)
x_min_med = x[mask].min()
x_max_med = x[mask].max()

mask = (x > 1.5)
x_min_max = x[mask].min()
x_max_max = x[mask].max()

pop_affected_F_min_depth = np.nansum(gdf_pop_2019_exposed_F['population'][(gdf_pop_2019_exposed_F['flood_depth'] > x_min_min) & (gdf_pop_2019_exposed_F['flood_depth'] < x_max_min)])
pop_affected_CF_pop_min_depth = np.nansum(gdf_pop_1990_exposed_F['population'][(gdf_pop_1990_exposed_F['flood_depth'] > x_min_min) & (gdf_pop_1990_exposed_F['flood_depth'] < x_max_min)])
pop_affected_CF_clim_min_depth = np.nansum(gdf_pop_2019_exposed_CF['population'][(gdf_pop_2019_exposed_CF['flood_depth'] > x_min_min) & (gdf_pop_2019_exposed_CF['flood_depth'] < x_max_min)])
pop_affected_CF_clim_pop_min_depth = np.nansum(gdf_pop_1990_exposed_CF['population'][(gdf_pop_1990_exposed_CF['flood_depth'] > x_min_min) & (gdf_pop_1990_exposed_CF['flood_depth'] < x_max_min)])
pop_affected_F_min_depth_rel = pop_affected_F_min_depth / np.nansum(pop_arrays[2019]) * 100
pop_affected_CF_pop_min_depth_rel = pop_affected_CF_pop_min_depth / np.nansum(pop_arrays[1990]) * 100
pop_affected_CF_clim_min_depth_rel = pop_affected_CF_clim_min_depth / np.nansum(pop_arrays[2019]) * 100
pop_affected_CF_clim_pop_min_depth_rel = pop_affected_CF_clim_pop_min_depth / np.nansum(pop_arrays[1990]) * 100

pop_affected_F_med_depth = np.nansum(gdf_pop_2019_exposed_F['population'][(gdf_pop_2019_exposed_F['flood_depth'] > x_min_med) & (gdf_pop_2019_exposed_F['flood_depth'] < x_max_med)])
pop_affected_CF_pop_med_depth = np.nansum(gdf_pop_1990_exposed_F['population'][(gdf_pop_1990_exposed_F['flood_depth'] > x_min_med) & (gdf_pop_1990_exposed_F['flood_depth'] < x_max_med)])
pop_affected_CF_clim_med_depth = np.nansum(gdf_pop_2019_exposed_CF['population'][(gdf_pop_2019_exposed_CF['flood_depth'] > x_min_med) & (gdf_pop_2019_exposed_CF['flood_depth'] < x_max_med)])
pop_affected_CF_clim_pop_med_depth = np.nansum(gdf_pop_1990_exposed_CF['population'][(gdf_pop_1990_exposed_CF['flood_depth'] > x_min_med) & (gdf_pop_1990_exposed_CF['flood_depth'] < x_max_med)])
pop_affected_F_med_depth_rel = pop_affected_F_med_depth / np.nansum(pop_arrays[2019]) * 100
pop_affected_CF_pop_med_depth_rel = pop_affected_CF_pop_med_depth / np.nansum(pop_arrays[1990]) * 100
pop_affected_CF_clim_med_depth_rel = pop_affected_CF_clim_med_depth / np.nansum(pop_arrays[2019]) * 100
pop_affected_CF_clim_pop_med_depth_rel = pop_affected_CF_clim_pop_med_depth / np.nansum(pop_arrays[1990]) * 100

pop_affected_F_max_depth = np.nansum(gdf_pop_2019_exposed_F['population'][(gdf_pop_2019_exposed_F['flood_depth'] > x_min_max) & (gdf_pop_2019_exposed_F['flood_depth'] < x_max_max)])
pop_affected_CF_pop_max_depth = np.nansum(gdf_pop_1990_exposed_F['population'][(gdf_pop_1990_exposed_F['flood_depth'] > x_min_max) & (gdf_pop_1990_exposed_F['flood_depth'] < x_max_max)])
pop_affected_CF_clim_max_depth = np.nansum(gdf_pop_2019_exposed_CF['population'][(gdf_pop_2019_exposed_CF['flood_depth'] > x_min_max) & (gdf_pop_2019_exposed_CF['flood_depth'] < x_max_max)])
pop_affected_CF_clim_pop_max_depth = np.nansum(gdf_pop_1990_exposed_CF['population'][(gdf_pop_1990_exposed_CF['flood_depth'] > x_min_max) & (gdf_pop_1990_exposed_CF['flood_depth'] < x_max_max)])
pop_affected_F_max_depth_rel = pop_affected_F_max_depth / np.nansum(pop_arrays[2019]) * 100
pop_affected_CF_pop_max_depth_rel = pop_affected_CF_pop_max_depth / np.nansum(pop_arrays[1990]) * 100
pop_affected_CF_clim_max_depth_rel = pop_affected_CF_clim_max_depth / np.nansum(pop_arrays[2019]) * 100
pop_affected_CF_clim_pop_max_depth_rel = pop_affected_CF_clim_pop_max_depth / np.nansum(pop_arrays[1990]) * 100

total_affected_F = np.nansum(gdf_pop_2019_exposed_F['population'][gdf_pop_2019_exposed_F['flood_depth'] > 0])
total_affected_CF_clim = np.nansum(gdf_pop_2019_exposed_CF['population'][gdf_pop_2019_exposed_CF['flood_depth'] > 0])
total_affected_CF_pop = np.nansum(gdf_pop_1990_exposed_F['population'][gdf_pop_1990_exposed_F['flood_depth'] > 0])
total_affected_CF_clim_pop = np.nansum(gdf_pop_1990_exposed_CF['population'][gdf_pop_1990_exposed_CF['flood_depth'] > 0])

# Build summary table
# Depth ranges
depth_ranges = ["0-0.5 m", "0.5-1.5 m", ">1.5 m"]

# Labels for sub-rows
labels = ["# affected", "% of total", "# attributable", "% attributed"]

# Prepare values for each scenario
F_abs = [pop_affected_F_min_depth, pop_affected_F_med_depth, pop_affected_F_max_depth]
CF_pop_abs = [pop_affected_CF_pop_min_depth, pop_affected_CF_pop_med_depth, pop_affected_CF_pop_max_depth]
CF_clim_abs = [pop_affected_CF_clim_min_depth, pop_affected_CF_clim_med_depth, pop_affected_CF_clim_max_depth]
CF_clim_pop_abs = [pop_affected_CF_clim_pop_min_depth, pop_affected_CF_clim_pop_med_depth, pop_affected_CF_clim_pop_max_depth]

F_rel = [pop_affected_F_min_depth_rel, pop_affected_F_med_depth_rel, pop_affected_F_max_depth_rel]
CF_pop_rel = [pop_affected_CF_pop_min_depth_rel, pop_affected_CF_pop_med_depth_rel, pop_affected_CF_pop_max_depth_rel]
CF_clim_rel = [pop_affected_CF_clim_min_depth_rel, pop_affected_CF_clim_med_depth_rel, pop_affected_CF_clim_max_depth_rel]
CF_clim_pop_rel = [pop_affected_CF_clim_pop_min_depth_rel, pop_affected_CF_clim_pop_med_depth_rel, pop_affected_CF_clim_pop_max_depth_rel]

# Attributable values
attr_CF_pop_abs = [F_abs[i] - CF_pop_abs[i] for i in range(3)]
attr_CF_clim_abs = [F_abs[i] - CF_clim_abs[i] for i in range(3)]
attr_CF_clim_pop_abs = [F_abs[i] - CF_clim_pop_abs[i] for i in range(3)]

attr_CF_pop_pct = [attr_CF_pop_abs[i]/F_abs[i]*100 for i in range(3)]
attr_CF_clim_pct = [attr_CF_clim_abs[i]/F_abs[i]*100 for i in range(3)]
attr_CF_clim_pop_pct = [attr_CF_clim_pop_abs[i]/F_abs[i]*100 for i in range(3)]

# Build table rows
rows = []
for i, dr in enumerate(depth_ranges):
    rows.append({"Depth Range": dr, "Label": labels[0], "F": F_abs[i], "CF pop": CF_pop_abs[i], "CF clim": CF_clim_abs[i], "CF clim & pop": CF_clim_pop_abs[i]})
    rows.append({"Depth Range": dr, "Label": labels[1], "F": F_rel[i], "CF pop": CF_pop_rel[i], "CF clim": CF_clim_rel[i], "CF clim & pop": CF_clim_pop_rel[i]})
    rows.append({"Depth Range": dr, "Label": labels[2], "F": "-", "CF pop": attr_CF_pop_abs[i], "CF clim": attr_CF_clim_abs[i], "CF clim & pop": attr_CF_clim_pop_abs[i]})
    rows.append({"Depth Range": dr, "Label": labels[3], "F": "-", "CF pop": attr_CF_pop_pct[i], "CF clim": attr_CF_clim_pct[i], "CF clim & pop": attr_CF_clim_pop_pct[i]})

# Add extra row for total population
rows.append({"Depth Range": "Total population", "Label": "-", 
             "F": np.nansum(pop_arrays[2019]), 
             "CF pop": np.nansum(pop_arrays[1990]), 
             "CF clim": np.nansum(pop_arrays[2019]), 
             "CF clim & pop": np.nansum(pop_arrays[1990])})

# Create DataFrame
df_table = pd.DataFrame(rows)
df_table.to_csv("c:/Code/Paper_2/Output/flood_depth_impact_summary_table.csv", index=False)
# Display summary table
pd.set_option('display.max_colwidth', None)
print(df_table)


#%% ============================================================================================ # 
# =================== Plot the flood depth per exposed population spatially ==================== #
# ============================================================================================== #
# Plot average flood depth per exposed population cell
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True) # 10,6

# Plot average flood depth among exposed population
gdf_2019_F = gdf_pop_2019_exposed_F_coarse.copy()
gdf_2019_F.loc[gdf_2019_F["exposed_population"] <= 0, "avg_flood_depth"] = np.nan

gdf_2019_CF = gdf_pop_2019_exposed_CF_coarse.copy()
gdf_2019_CF.loc[gdf_2019_CF["exposed_population"] <= 0, "avg_flood_depth"] = np.nan

gdf_2019_CF['change_in_flood_depth'] = gdf_2019_F['avg_flood_depth'] - gdf_2019_CF['avg_flood_depth']

# Define colormap: from white to #67CBE4
cmap = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#67CBE4"])
cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#651F94"])

plot = gdf_2019_F.plot(column="avg_flood_depth", cmap=cmap, vmin=0, vmax=3.5, linewidth=0.1, 
                edgecolor="grey", ax=axes[0], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_CF.plot(column="avg_flood_depth", cmap=cmap, vmin=0, vmax=3.5, linewidth=0.1, 
                edgecolor="grey", ax=axes[1], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_CF.plot(column="change_in_flood_depth", cmap=cmap_change, vmin=0, vmax=0.5, linewidth=0.1, 
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

gdf_2019_F = gdf_pop_2019_exposed_F_coarse.copy()
gdf_2019_F.loc[gdf_2019_F["total_population"] == 0] = np.nan

gdf_1990_F = gdf_pop_1990_exposed_F_coarse.copy()
gdf_1990_F.loc[gdf_1990_F["total_population"] == 0] = np.nan

gdf_1990_F['change_in_population'] = gdf_2019_F['total_population'] - gdf_1990_F['total_population']

# Define colormap: from white to #67CBE4
cmap = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#FC6F37"])
norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_pop_2019_exposed_F_coarse['total_population']))  
cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#BD2A2A"])
norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_1990_F['change_in_population']))

plot = gdf_2019_F.plot(column="total_population", cmap=cmap, norm=norm ,
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
fig, axes = plt.subplots(1, 3, figsize=(16, 6), sharey=True, constrained_layout=True)

gdf_2019_F = gdf_pop_2019_exposed_F_coarse.copy()
gdf_2019_F.loc[gdf_2019_F["total_population"] == 0] = np.nan

gdf_2019_F_uni = gdf_pop_2019_exposed_F_uniform_coarse.copy()
gdf_2019_F_uni.loc[gdf_2019_F_uni["total_population"] == 0] = np.nan

gdf_2019_F_uni['change_in_population'] = gdf_2019_F_uni['total_population'] - gdf_2019_F['total_population']

# Define colormap: from white to #67CBE4
cmap = mcolors.LinearSegmentedColormap.from_list("white_to_orange", ["#ffffff", "#651F94"])
norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_pop_2019_exposed_F_coarse['total_population']))  
cmap_change = plt.get_cmap('RdBu_r')
norm_change = mcolors.TwoSlopeNorm(vmin=np.nanmin(gdf_2019_F_uni['change_in_population']),
                                   vcenter=0,
                                   vmax=np.nanmax(gdf_2019_F_uni['change_in_population']))

plot = gdf_2019_F.plot(column="total_population", cmap=cmap, norm=norm ,
                       linewidth=0.1, edgecolor="grey", ax=axes[0], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_F_uni.plot(column="total_population", cmap=cmap, norm=norm,
                        linewidth=0.1, edgecolor="grey", ax=axes[1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_F_uni.plot(column="change_in_population", cmap=cmap_change, norm=norm_change,
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

# --- Population colorbar spanning axes[0] and axes[1]
sm_pop = ScalarMappable(cmap=cmap, norm=norm)
sm_pop._A = []

cbar_pop = fig.colorbar(sm_pop, ax=[axes[0], axes[1]], shrink=0.8, pad=0.02)
cbar_pop.set_label("Population (people per cell)")

# --- Difference colorbar for axis[2] only
sm_change = ScalarMappable(cmap=cmap_change, norm=norm_change)
sm_change._A = []

cbar_change = fig.colorbar(sm_change, ax=axes[2], shrink=0.8, pad=0.02)
cbar_change.set_label("Difference in population (people per cell)")

axes[0].set_title("Factual 2019 population", fontsize=10)
axes[1].set_title("Uniform growth 2019 population", fontsize=10)
axes[2].set_title("Uniform growth - Factual", fontsize=10)

axes[0].set_ylabel("y coordinate UTM zone 36S [m]")
axes[0].set_xlabel("x coordinate UTM zone 36S [m]")
axes[1].set_xlabel("x coordinate UTM zone 36S [m]")
axes[2].set_xlabel("x coordinate UTM zone 36S [m]")

fig.suptitle("Total population in study region", fontsize=12)

plt.show()


#%%
# Do the same but zoom into Beira for the 25 m resolution rasters

# bbox in raster CRS
xmin, xmax = 690000, 702000
ymin, ymax = 7803000, 7818000

# Use rowcol safely
row_ul, col_ul = rowcol(flood_grid_transform, xmin, ymax, op=int)  # upper-left
row_lr, col_lr = rowcol(flood_grid_transform, xmax, ymin, op=int)  # lower-right

# Clip indices to raster shape
row_min = max(0, min(row_ul, row_lr))
row_max = min(ra_exposed_pop_2019_F.shape[0], max(row_ul, row_lr))
col_min = max(0, min(col_ul, col_lr))
col_max = min(ra_exposed_pop_2019_F.shape[1], max(col_ul, col_lr))

print("Row indices:", row_min, row_max)
print("Col indices:", col_min, col_max)

# Compute coordinates of the clipped raster edges
x_min_clip, y_max_clip = flood_grid_transform * (col_min, row_min)  # top-left
x_max_clip, y_min_clip = flood_grid_transform * (col_max, row_max)  # bottom-right

extent_beira = [x_min_clip, x_max_clip, y_max_clip, y_min_clip]
print("Beira extent:", extent_beira)

# Slice raster
raster_F_beira_pop    = pop_arrays[2019][row_min:row_max, col_min:col_max]
raster_UF_beira_pop   = pop_array_uniform_2019[row_min:row_max, col_min:col_max]
pop_diff_UG_beira_pop = raster_UF_beira_pop - raster_F_beira_pop

raster_F_beira_affected    = ra_exposed_pop_2019_F[row_min:row_max, col_min:col_max]
raster_UF_beira_affected   = ra_exposed_pop_2019_F_uniform[row_min:row_max, col_min:col_max]
pop_diff_UG_beira_affected = raster_UF_beira_affected - raster_F_beira_affected

print("Clipped raster shape:", raster_F_beira_affected.shape)


# Plotting
fig, axes = plt.subplots(2, 3, figsize=(16, 12), sharex=True, sharey=True, dpi=300, constrained_layout=True)

for i, ax in enumerate(axes.flatten()):
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
    region_utm.boundary.plot(ax=ax, color='black', linewidth=1, zorder=2)

# --- Population colormap ---
cmap_white_orange = mcolors.LinearSegmentedColormap.from_list("white_to_orange", ["#ffffff", "#651F94"])
pop_bins = np.arange(0, np.nanmax(raster_F_beira_pop)+1, 1) 
pop_colors = cmap_white_orange(np.linspace(0, 1, len(pop_bins)-1))
pop_colors[0] = [1, 1, 1, 0]  # RGBA
pop_cmap_discrete = mcolors.ListedColormap(pop_colors)
pop_norm = BoundaryNorm(pop_bins, pop_cmap_discrete.N, extend='neither')

# --- Affected population colormap ---
pop_bins = np.arange(0, np.nanmax(raster_F_beira_affected)+1, 1)
pop_colors = plt.cm.Blues(np.linspace(0, 1, len(pop_bins)-1))
pop_colors[0] = [1, 1, 1, 0]  # RGBA
pop_affc_cmap_discrete = mcolors.ListedColormap(pop_colors)
pop_affc_norm = BoundaryNorm(pop_bins, pop_affc_cmap_discrete.N, extend='neither')

# --- Difference colormap population ---
vmax_ug = np.nanmax(pop_diff_UG_beira_pop)
diff_bins = np.arange(-int(vmax_ug), int(vmax_ug)+1, 1)  # integer steps
diff_colors = plt.cm.RdBu_r(np.linspace(0, 1, len(diff_bins)-1))
mid_idx = np.where(diff_bins[:-1] == 0)[0]
if len(mid_idx) > 0:
    diff_colors[mid_idx[0]] = [1, 1, 1, 0]
diff_cmap_discrete = mcolors.ListedColormap(diff_colors)
diff_norm = BoundaryNorm(diff_bins, diff_cmap_discrete.N, extend='neither')

# --- Difference colormap affected ---
vmax_ug = np.nanmax(pop_diff_UG_beira_affected)
diff_bins = np.arange(-int(vmax_ug), int(vmax_ug)+1, 1)  # integer steps
diff_colors = plt.cm.RdBu_r(np.linspace(0, 1, len(diff_bins)-1))
mid_idx = np.where(diff_bins[:-1] == 0)[0]
if len(mid_idx) > 0:
    diff_colors[mid_idx[0]] = [1, 1, 1, 0]
diff_affc_cmap_discrete = mcolors.ListedColormap(diff_colors)
diff_affc_norm = BoundaryNorm(diff_bins, diff_affc_cmap_discrete.N, extend='neither')

# --- TOP ROW ---
# Population change effect
im = axes[0,0].imshow(raster_F_beira_pop, cmap=pop_cmap_discrete, norm=pop_norm,
                      extent=extent_beira, origin='lower', alpha=0.8, zorder=3)
axes[0,0].set_title("Factual 2019 Population")

im = axes[0,1].imshow(raster_UF_beira_pop, cmap=pop_cmap_discrete, norm=pop_norm,
                      extent=extent_beira, origin='lower', alpha=0.8, zorder=3)
cbar = plt.colorbar(im, ax=axes[0,1], shrink=0.8)
cbar.set_label("Exposed population")
axes[0,1].set_title("Uniform 2019 Population")

im = axes[0,2].imshow(pop_diff_UG_beira_pop, cmap=diff_cmap_discrete, norm=diff_norm,
                      extent=extent_beira, origin='lower', alpha=0.8, zorder=3)
cbar = plt.colorbar(im, ax=axes[0,2], shrink=0.8)
cbar.set_label("Difference in population")
axes[0,2].set_title("Uniform - Factual population")

# # --- BOTTOM ROW --
# Population change effect
im = axes[1,0].imshow(raster_F_beira_affected, cmap=pop_affc_cmap_discrete, norm=pop_affc_norm,
                      extent=extent_beira, origin='lower', alpha=0.8, zorder=3)
axes[1,0].set_title("Factual affected population")

im = axes[1,1].imshow(raster_UF_beira_affected, cmap=pop_affc_cmap_discrete, norm=pop_affc_norm,
                      extent=extent_beira, origin='lower', alpha=0.8, zorder=3)
cbar = plt.colorbar(im, ax=axes[1,1], shrink=0.8)
cbar.set_label("Affected population per 25 m grid cell")
axes[1,1].set_title("Factual affected population Uniform")
im = axes[1,2].imshow(pop_diff_UG_beira_affected, cmap=diff_affc_cmap_discrete, norm=diff_affc_norm,
                      extent=extent_beira, origin='lower', alpha=0.8, zorder=3)
cbar = plt.colorbar(im, ax=axes[1,2], shrink=0.8)
cbar.set_label("Difference in affected population")
axes[1,2].set_title("Uniform - Factual affected population")

axes[0,0].set_ylabel("y coordinate UTM zone 36S [m]")
axes[1,0].set_ylabel("y coordinate UTM zone 36S [m]")
axes[1,0].set_xlabel("x coordinate UTM zone 36S [m]")
axes[1,1].set_xlabel("x coordinate UTM zone 36S [m]")
axes[1,2].set_xlabel("x coordinate UTM zone 36S [m]")

fig.suptitle("Uniform vs. spatially differing population growth", fontsize=14)

plt.show()


#%%

# Colours based on conceptual figure
colours = ['#00B050', "#086632"]

fig, ax_main = plt.subplots(1, 1, figsize=(8, 6),dpi=300)


ax_main.plot(x, y_F, label=f"Factual population ({np.nansum(pop_arrays[2019]):,.0f} people)", color=colours[0], linewidth=2)
ax_main.plot(x, y_F_uni, label=f"Uniform population ({np.nansum(pop_array_uniform_2019):,.0f} people)", color=colours[1], linewidth=1)

ax_main.set_ylabel("Density-weighted population exposure")
ax_main.grid(True, linestyle="--", alpha=0.5)
ax_main.set_xlim(0, 3.5)
ax_main.set_xlabel("Flood depth")
ax_main.set_ylim(0, np.nanmax(y_CF_clim_pop)*1.05)

# peak annotations 
if len(sorted_peaks_F) > 1:
    second_peak_x = x[sorted_peaks_F[1]]
    plt.annotate(f"{x[sorted_peaks_F[0]]:.2f} m",
                 xy=(x[sorted_peaks_F[0]]+0.02, y_F[sorted_peaks_F[0]]),
                 xytext=(x[sorted_peaks_F[0]]+0.15, y_F[sorted_peaks_F[0]]-0.05), color=colours[0], fontsize=8,
                 arrowprops=dict(arrowstyle="->", lw=0.8, color=colours[0], linestyle='--'))
    plt.annotate(f"{x[sorted_peaks_F[1]]:.2f} m", 
                 xy=(x[sorted_peaks_F[1]]+0.02, y_F[sorted_peaks_F[1]]),
                 xytext=(x[sorted_peaks_F[1]]+0.15, y_F[sorted_peaks_F[1]]), color=colours[0], fontsize=8,
                 arrowprops=dict(arrowstyle="->", lw=0.8, color=colours[0], linestyle='--'))

if len(sorted_peaks_F_uni) > 1:
    second_peak_x = x[sorted_peaks_F_uni[1]]
for i, idx in enumerate(sorted_peaks_F_uni[:2]): 
    plt.annotate(f"{x[idx]:.2f} m",
                 xy=(x[idx]+0.02, y_F_uni[idx]),
                 xytext=(x[idx]+0.15, y_F_uni[idx]+0.05), color=colours[1], fontsize=8,
                 arrowprops=dict(arrowstyle="->", lw=0.8, color=colours[1], linestyle='--'))

ax_main.legend(loc='upper right', fontsize=9)


#%% ============================================================================================ # 
# =============== Plot the change in exposed population > 1 m flood depth spatially ============ #
# ============================================================================================== #
# Plot average flood depth per exposed population cell
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True) # 10,6

# Plot average flood depth among exposed population
gdf_2019_F = gdf_pop_2019_exposed_F_coarse.copy()
gdf_2019_F.loc[gdf_2019_F["exposed_population"] <= 1, "avg_flood_depth"] = np.nan

gdf_2019_CF = gdf_pop_2019_exposed_CF_coarse.copy()
gdf_2019_CF.loc[gdf_2019_CF["exposed_population"] <= 1, "avg_flood_depth"] = np.nan

gdf_2019_CF['change_in_exposed_population_>1m'] = gdf_2019_F['exposed_population'] - gdf_2019_CF['exposed_population']

# Define colormap: from white to #67CBE4
cmap = 'Blues'
norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2019_F['exposed_population']))
cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#651F94"])
norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2019_CF['change_in_exposed_population_>1m']))

plot = gdf_2019_F.plot(column="exposed_population", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[0], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_CF.plot(column="exposed_population", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[1], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_CF.plot(column="change_in_exposed_population_>1m", cmap=cmap_change, norm=norm_change, linewidth=0.1, 
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
gdf_2019_F = gdf_pop_2019_exposed_F_coarse.copy()
gdf_2019_CF = gdf_pop_2019_exposed_CF_coarse.copy()

gdf_2019_CF['change_in_%more_1m'] = gdf_2019_F['pct_cells_higher_1m'] - gdf_2019_CF['pct_cells_higher_1m']

# Define colormap: from white to #67CBE4
cmap = 'Blues'
norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2019_F['pct_cells_higher_1m']))
cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#651F94"])
norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2019_CF['change_in_%more_1m']))

plot = gdf_2019_F.plot(column="pct_cells_higher_1m", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[0], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_CF.plot(column="pct_cells_higher_1m", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[1], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_CF.plot(column="change_in_%more_1m", cmap=cmap_change, norm=norm_change, linewidth=0.1, 
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
im = axes[1].imshow(ra_exposed_pop_2019_F, cmap='viridis', extent=flood_extent, origin='lower')
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
pop_per_district_adm3 = []

for _, row in districts_adm3_filtered.iterrows():
    # Mask raster to district polygon
    from rasterio import features
    district_mask = features.geometry_mask([row.geometry],
                                           out_shape=ra_exposed_pop_2019_F.shape,
                                           transform=flood_grid_transform,
                                           invert=True)
    pop_exposed.append(ra_exposed_pop_2019_F[district_mask].sum())
    pop_totals.append(pop_arrays[2019][district_mask].sum())

districts_adm3_filtered['pop_exposed'] = pop_exposed
districts_adm3_filtered['pop_total'] = pop_totals

# Make sure districts are in the same CRS as the original pop raster
with rasterio.open(population_raster_path_2019) as src:
    pop_crs = src.crs

districts_native = districts_adm3_filtered.to_crs(pop_crs)

for _, row in districts_native.iterrows():
    district_mask = features.geometry_mask(
        [row.geometry],
        out_shape=pop_sofala_districts_adm3[2019][0].shape,
        transform=pop_affine_sofala_districts_adm3[2019],
        invert=True
    )
    pop_per_district_adm3.append(pop_sofala_districts_adm3[2019][0][district_mask].sum())

districts_adm3_filtered['pop_per_district'] = pop_per_district_adm3

#%%
# Do the same for admin 2 level
pop_totals = []
pop_exposed = []
pop_per_district_adm2 = []

for _, row in districts_adm2.iterrows():
    # Mask raster to district polygon
    from rasterio import features
    district_mask = features.geometry_mask([row.geometry],
                                           out_shape=ra_exposed_pop_2019_F.shape,
                                           transform=flood_grid_transform,
                                           invert=True)
    pop_exposed.append(ra_exposed_pop_2019_F[district_mask].sum())
    pop_totals.append(pop_arrays[2019][district_mask].sum())

districts_adm2['pop_exposed'] = pop_exposed
districts_adm2['pop_total'] = pop_totals

# Make sure districts are in the same CRS as the original pop raster
with rasterio.open(population_raster_path_2019) as src:
    pop_crs = src.crs

districts_native = districts_adm2.to_crs(pop_crs)

for _, row in districts_native.iterrows():
    district_mask = features.geometry_mask(
        [row.geometry],
        out_shape=pop_sofala_districts_adm2[2019][0].shape,
        transform=pop_affine_sofala_districts_adm2[2019],
        invert=True
    )
    pop_per_district_adm2.append(pop_sofala_districts_adm2[2019][0][district_mask].sum())

districts_adm2['pop_per_district'] = pop_per_district_adm2

#%%
# --- Plot ---
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

for ax in axes:
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
    districts_adm3_filtered.boundary.plot(ax=ax, color='orange', linewidth=2, zorder=2)
    region_utm.boundary.plot(ax=ax, color='lightblue', linewidth=0.5, zorder=2)

# Rasterize the case-study region polygon
region_mask = rasterio.features.rasterize(
    [(geom, 1) for geom in region.geometry],
    out_shape=ra_exposed_pop_2019_F.shape,
    transform=flood_grid_transform,
    fill=0,
    all_touched=True,
    dtype=np.uint8
).astype(bool)

# Mask raster outside region
ra_exposed_pop_masked = np.where(region_mask, ra_exposed_pop_2019_F, np.nan)
    
# Exposed population
im = axes[0].imshow(ra_exposed_pop_masked, cmap='viridis', extent=flood_extent, origin='lower')
cbar = plt.colorbar(im, ax=axes[0], shrink=0.8)
cbar.set_label("Exposed population")
axes[0].set_title("Factual Exposed Population")

# Exposed population
im = axes[1].imshow(pop_arrays[2019], cmap='Reds', extent=flood_extent, origin='lower')
cbar = plt.colorbar(im, ax=axes[1], shrink=0.8)
cbar.set_label("Total population")
axes[1].set_title("Total 2019 Population")

# Dictionary with (dx, dy) offsets in map units for specific districts
label_offsets = {
    "Sofala": (0, 10000),  # move 0 m east, 10000 m north
    "Nhamatanda": (18000, -22000),
    "Estaquinha": (20000, -5000),
    # add more as needed
}

outside_districts = ["Nhamatanda", "Estaquinha"]

for idx, row in districts_adm3_filtered.iterrows():
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
districts_adm3_filtered.boundary.plot(ax=ax, color='orange', linewidth=2, zorder=2)
region_utm.boundary.plot(ax=ax, color='lightblue', linewidth=1, zorder=2)

# Dictionary with (dx, dy) offsets in map units for specific districts
label_offsets = {
    "Sofala": (0, 10000),  # move 0 m east, 10000 m north
    "Nhamatanda": (18000, -22000),
    "Estaquinha": (20000, -5000),
    # add more as needed
}
outside_districts = ["Nhamatanda", "Estaquinha"]

for idx, row in districts_adm3_filtered.iterrows():
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

ax.set_title("2019 District Population")

plt.tight_layout()


#%%
# --- Plot ---
fig, ax = plt.subplots(1, 1, figsize=(12, 6))

# background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
# bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
districts_adm2_utm = districts_adm2.to_crs(region_utm.crs)
districts_adm2_utm.boundary.plot(ax=ax, color='orange', linewidth=2, zorder=2)
region.boundary.plot(ax=ax, color='lightblue', linewidth=1, zorder=2)

# Dictionary with (dx, dy) offsets in map units for specific districts
label_offsets = {
    "Sofala": (0, 10000),  # move 0 m east, 10000 m north
    "Nhamatanda": (18000, -22000),
    "Estaquinha": (20000, -5000),
    # add more as needed
}

for idx, row in districts_adm2_utm.iterrows():
    x, y = row.geometry.centroid.x, row.geometry.centroid.y
    ax.text(x, y, f"{row['pop_per_district']:,.0f}", fontsize=8, ha='center', va='center',
            color='black', fontweight='bold', zorder=5)
    
    # District name label
    name_x = x
    name_y = y - 7000  # all 3000 south of pop label

    ax.text(
        name_x, name_y, row['NAME_2'], fontsize=8, ha='center', va='center',
        color='grey', fontweight='bold', zorder=5
    )

ax.set_title("2019 District Population")

plt.tight_layout()


#%%
# table with numbers per district
df_district_summary = districts_adm3_filtered[['NAME_3', 'pop_per_district', 'pop_total', 'pop_exposed']].copy()
df_district_summary.columns = ['District', 'District Population (2019)', 'Total Population in Region', 'Exposed Population (2019 Factual)']
df_district_summary["% of district pop"] = (100 * df_district_summary['Total Population in Region'] / df_district_summary['District Population (2019)']).round(0)
df_district_summary["% exposed"] = (100 * df_district_summary['Exposed Population (2019 Factual)'] / df_district_summary['Total Population in Region']).round(0)
df_district_summary[['District', 'District Population (2019)', 'Total Population in Region', 'Exposed Population (2019 Factual)']] = df_district_summary[['District', 'District Population (2019)', 'Total Population in Region', 'Exposed Population (2019 Factual)']].round(0)

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

pop_diff_PG = ra_exposed_pop_2019_F - ra_exposed_pop_1990_F
pop_diff_CC = ra_exposed_pop_2019_F - ra_exposed_pop_2019_CF

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
pop_bins = np.arange(0, np.nanmax(ra_exposed_pop_2019_F)+1, 1)  # 0,1,2,... max
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
im = axes[0,0].imshow(np.round(ra_exposed_pop_2019_F).astype(int), cmap=pop_cmap_discrete, norm=pop_norm,
                      extent=flood_extent, origin='lower', alpha=0.8, zorder=3)
axes[0,0].set_title("Exposed 2019 Population")

im = axes[0,1].imshow(np.round(ra_exposed_pop_1990_F).astype(int), cmap=pop_cmap_discrete, norm=pop_norm,
                      extent=flood_extent, origin='lower', alpha=0.8, zorder=3)
cbar = plt.colorbar(im, ax=axes[0,1], shrink=0.8)
cbar.set_label("Exposed population")
axes[0,1].set_title("Exposed 1990 Population")

im = axes[0,2].imshow(np.round(pop_diff_PG).astype(int), cmap=diff_cmap_discrete, norm=diff_norm,
                      extent=flood_extent, origin='lower', alpha=0.8, zorder=3)
cbar = plt.colorbar(im, ax=axes[0,2], shrink=0.8)
cbar.set_label("Attributable exposed population")
axes[0,2].set_title("Population Change (2019 - 1990)")

# # --- BOTTOM ROW --
# Climate change effect
im = axes[1,0].imshow(np.round(ra_exposed_pop_2019_F).astype(int), cmap=pop_cmap_discrete, norm=pop_norm, extent=flood_extent,
                      origin='lower', alpha=0.8, zorder=3)
axes[1,0].set_title("Factual Climate")

im = axes[1,1].imshow(np.round(ra_exposed_pop_2019_CF).astype(int), cmap=pop_cmap_discrete, norm=pop_norm, extent=flood_extent,
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

gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['exposed_population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
norm = PowerNorm(gamma=0.5, vmin=gdf_pop_2019_exposed_F_coarse['exposed_population'].min(), vmax=gdf_pop_2019_exposed_F_coarse['exposed_population'].max())

plot = gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['exposed_population'] > 0].plot(column='exposed_population', cmap='Blues', edgecolor='grey',
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
vmin, vmax = gdf_pop_2019_exposed_F_coarse["exposed_population"].min(), gdf_pop_2019_exposed_F_coarse["exposed_population"].max()
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

gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['relative_population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['relative_population'] > 0].plot(column='relative_population', cmap='Blues', edgecolor='grey',
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
vmin, vmax = gdf_pop_2019_exposed_F_coarse["relative_population"].min(), gdf_pop_2019_exposed_F_coarse["relative_population"].max()
sm1 = ScalarMappable(cmap="Blues", norm=Normalize(vmin=vmin, vmax=vmax))
sm1._A = []  # required for colorbar without passing data
cbar = fig.colorbar(sm1, ax=axes[1], orientation="vertical", shrink=0.8)
cbar.set_label("Relative exposed population (%)", fontsize=10)
cbar.ax.tick_params(labelsize=9)

axes[0].set_title("Factual flooding", fontsize=11)
axes[1].set_title("Factual population exposure", fontsize=11)



# %%
print("Plotting spatially aggregated exposed population for F, CF population and diff")

gdf_pop_1990_exposed_F_coarse['population_diff'] = (gdf_pop_2019_exposed_F_coarse['exposed_population'] - gdf_pop_1990_exposed_F_coarse['exposed_population'])

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 6), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Create colormap normalization
norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2019_exposed_F_coarse["exposed_population"].max())
norm_diff = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_1990_exposed_F_coarse["population_diff"].max())
red_half = LinearSegmentedColormap.from_list("bwr_red", plt.cm.bwr(np.linspace(0.5, 1, 256)))

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['exposed_population'] == 0].plot(ax=axes[0], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['exposed_population'] > 0].plot(column='exposed_population', cmap='Blues',  edgecolor='grey',
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
cbar.set_label("Aggregated exposed population [# people]", fontsize=10)
cbar.ax.tick_params(labelsize=9)

sm2 = ScalarMappable(cmap=red_half, norm=norm_diff)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, ax=axes[2], orientation="vertical", shrink=0.5)
cbar2.set_label("Attributable exposed population [# people]", fontsize = 9)
cbar2.ax.tick_params(labelsize=8)

axes[0].set_title("Factual population \n(2019)", fontsize=9)
axes[1].set_title("Counterfactual population \n(1990)", fontsize=9)
axes[2].set_title("Factual - Counterfactual", fontsize=9)

# %%
print("Plotting spatially aggregated exposed population for F, CF climate and diff")

gdf_pop_2019_exposed_CF_coarse['population_diff'] = (gdf_pop_2019_exposed_F_coarse['exposed_population'] - gdf_pop_2019_exposed_CF_coarse['exposed_population'])

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 5), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Create colormap normalization
norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2019_exposed_F_coarse["exposed_population"].max())
norm_diff = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2019_exposed_CF_coarse["population_diff"].max())
red_half = LinearSegmentedColormap.from_list("bwr_red", plt.cm.bwr(np.linspace(0.5, 1, 256)))

cmap = LinearSegmentedColormap.from_list("custom_cmap", ["white", "#67CBE4"])

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['exposed_population'] == 0].plot(ax=axes[0], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['exposed_population'] > 0].plot(column='exposed_population', cmap='Blues',  edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[0], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_2019_exposed_CF_coarse[gdf_pop_2019_exposed_CF_coarse['exposed_population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2019_exposed_CF_coarse[gdf_pop_2019_exposed_CF_coarse['exposed_population'] > 0].plot(column='exposed_population', cmap='Blues', edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[1], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_2019_exposed_CF_coarse[gdf_pop_2019_exposed_CF_coarse['population_diff'] <= 0].plot(ax=axes[2], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2019_exposed_CF_coarse[gdf_pop_2019_exposed_CF_coarse['population_diff'] > 0].plot(column='population_diff', cmap=red_half, norm=norm_diff, edgecolor='grey', 
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
cbar.set_label("Aggregated exposed population [# people]", fontsize=10)
cbar.ax.tick_params(labelsize=9)

sm2 = ScalarMappable(cmap=red_half, norm=norm_diff)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, ax=axes[2], orientation="vertical", shrink=0.5)
cbar2.set_label("Attributable exposed population [# people]", fontsize = 9)
cbar2.ax.tick_params(labelsize=8)

axes[0].set_title("Factual climate", fontsize=9)
axes[1].set_title("Counterfactual climate", fontsize=9)
axes[2].set_title("Factual - Counterfactual", fontsize=9)



# %%
# Plotting uniform population growth vs spatially differing growth
print("Plotting spatially aggregated exposed population for F, CF population and diff")

gdf_pop_2019_exposed_F_uniform_coarse['population_diff'] = (gdf_pop_2019_exposed_F_coarse['exposed_population'] - gdf_pop_2019_exposed_F_uniform_coarse['exposed_population'])

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 6), dpi=300, sharey=True, constrained_layout=True,
                         subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Create colormap normalization
norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2019_exposed_F_coarse["exposed_population"].max())
norm_2019_diff = mcolors.TwoSlopeNorm(vmin=(gdf_pop_2019_exposed_F_uniform_coarse["population_diff"].min()), vcenter=0, vmax=(gdf_pop_2019_exposed_F_uniform_coarse["population_diff"].max()))
cmap_2019_diff = plt.get_cmap('RdBu_r')


# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['exposed_population'] == 0].plot(ax=axes[0], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['exposed_population'] > 0].plot(column='exposed_population', cmap='Blues',  edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[0], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_2019_exposed_F_uniform_coarse[gdf_pop_2019_exposed_F_uniform_coarse['exposed_population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2019_exposed_F_uniform_coarse[gdf_pop_2019_exposed_F_uniform_coarse['exposed_population'] > 0].plot(column='exposed_population', cmap='Blues', edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[1], legend=False, 
                                                                               zorder=2, rasterized=True)
plot = gdf_pop_2019_exposed_F_uniform_coarse.plot(column='population_diff', cmap=cmap_2019_diff, norm=norm_2019_diff, edgecolor='grey',
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

sm2 = ScalarMappable(cmap=cmap_2019_diff, norm=norm_2019_diff)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, ax=axes[2], orientation="vertical", shrink=0.5)
cbar2.set_label("Attributable exposed population [# people]", fontsize = 9)
cbar2.ax.tick_params(labelsize=8)

axes[0].set_title("Factual exposed population \n(2019)", fontsize=9)
axes[1].set_title("Factual exposed population \n(2019 - uniform)", fontsize=9)
axes[2].set_title("Spatial - Uniform", fontsize=9)

# %% #################################################################
##################### diff in relative population ####################
######################################################################
print("Plotting spatially aggregated RELATIVE exposed population damage for F, CF population and diff")
from matplotlib.colors import TwoSlopeNorm

gdf_pop_1990_exposed_F_coarse['population_diff'] = (gdf_pop_2019_exposed_F_coarse['relative_population'] - gdf_pop_1990_exposed_F_coarse['relative_population'])

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 6), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Create colormap normalization
vmin_diff = np.nanpercentile(gdf_pop_1990_exposed_F_coarse["population_diff"], 1)
vmax_diff = np.nanpercentile(gdf_pop_1990_exposed_F_coarse["population_diff"], 99)
norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2019_exposed_F_coarse["relative_population"].max())
norm_diff = TwoSlopeNorm(vmin=vmin_diff, vcenter=0, vmax=vmax_diff)
cmap_diff_pop = plt.get_cmap('RdBu_r')

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['relative_population'] == 0].plot(ax=axes[0], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['relative_population'] > 0].plot(column='relative_population', cmap='Blues',  edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[0], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_1990_exposed_F_coarse[gdf_pop_1990_exposed_F_coarse['relative_population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_1990_exposed_F_coarse[gdf_pop_1990_exposed_F_coarse['relative_population'] > 0].plot(column='relative_population', cmap='Blues', edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[1], legend=False, 
                                                                               zorder=2, rasterized=True)

plot = gdf_pop_1990_exposed_F_coarse.plot(
    column='population_diff',
    cmap=cmap_diff_pop,
    norm=norm_diff,
    edgecolor='grey',
    linewidth=0.2,
    ax=axes[2],
    legend=False,
    zorder=2,
    rasterized=True,
    missing_kwds={"color": "white", "edgecolor": "none"}
)

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
cbar.set_label("Aggregated relative exposed population [%]", fontsize=10)
cbar.ax.tick_params(labelsize=9)

sm2 = ScalarMappable(cmap=cmap_diff_pop, norm=norm_diff)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, ax=axes[2], orientation="vertical", shrink=0.5)
cbar2.set_label("Attributable relative exposed population [%]", fontsize = 9)
cbar2.ax.tick_params(labelsize=8)

axes[0].set_title("Factual population \n(2019)", fontsize=9)
axes[1].set_title("Counterfactual population \n(1990)", fontsize=9)
axes[2].set_title("Factual - Counterfactual", fontsize=9)

#%%
print("Plotting spatially RELATIVE aggregated exposed population for F, CF climate and diff")

gdf_pop_2019_exposed_CF_coarse['population_diff'] = (gdf_pop_2019_exposed_F_coarse['relative_population'] - gdf_pop_2019_exposed_CF_coarse['relative_population'])

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 5), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Create colormap normalization
norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2019_exposed_F_coarse["relative_population"].max())
norm_diff = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanpercentile(gdf_pop_2019_exposed_CF_coarse["population_diff"], 99))
red_half = LinearSegmentedColormap.from_list("bwr_red", plt.cm.bwr(np.linspace(0.5, 1, 256)))

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['relative_population'] == 0].plot(ax=axes[0], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['relative_population'] > 0].plot(column='relative_population', cmap='Blues',  edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[0], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_2019_exposed_CF_coarse[gdf_pop_2019_exposed_CF_coarse['relative_population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2019_exposed_CF_coarse[gdf_pop_2019_exposed_CF_coarse['relative_population'] > 0].plot(column='relative_population', cmap='Blues', edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[1], legend=False, 
                                                                               zorder=2, rasterized=True)

# gdf_pop_2019_exposed_CF_coarse[gdf_pop_2019_exposed_CF_coarse['population_diff'] <= 0].plot(ax=axes[2], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2019_exposed_CF_coarse[gdf_pop_2019_exposed_CF_coarse['population_diff'] > 0].plot(column='population_diff', cmap=red_half, norm=norm_diff, edgecolor='grey', 
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
cbar.set_label("Aggregated relative exposed population [%]", fontsize=10)
cbar.ax.tick_params(labelsize=9)

sm2 = ScalarMappable(cmap=red_half, norm=norm_diff)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, ax=axes[2], orientation="vertical", shrink=0.5)
cbar2.set_label("Attributable relative exposed population [%]", fontsize = 9)
cbar2.ax.tick_params(labelsize=8)

axes[0].set_title("Factual climate", fontsize=9)
axes[1].set_title("Counterfactual climate", fontsize=9)
axes[2].set_title("Factual - Counterfactual", fontsize=9)


# %%

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
from scipy.ndimage import gaussian_filter1d
from matplotlib.colors import TwoSlopeNorm


prefix = "p:/" if platform.system() == "Windows" else "/p/"

# ===== CONFIGURATION =====
EVENT_NAME = "Idai"
BASE_RUN_PATH = Path("p:/11210471-001-compass/03_Runs/sofala/Idai")
SCENARIO_PATH_F = "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" # factual
SCENARIO_PATH_CF = "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.1_era5_hourly_spw_IBTrACS_CF-5" # counterfactual

# ===== FILE PATHS =====
# Base directory for the specific event and scenario
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
region = gpd.read_file("p:/11210471-001-compass/03_Runs/sofala/Idai/sfincs/event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0/gis/region.geojson")
background = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_region_background.geojson", driver="GeoJSON")
shapefile_sofala = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_province.shp")

# Load the admin3 district in the case study region to validate exposed people
beira_district = gpd.read_file("../Attribution_results/data/gis/Beira_region.shp")
districts_adm3 = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_districts_study_region.shp")
districts_adm2 = data_catalog.get_geodataframe("gadm_level2", geom=region, buffer=1000)

# population in provided inthousand persons per grid cell
population_raster_path_2019 = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_2019_30.tif")  
population_raster_path_1990 = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_1990_30.tif")  

settlement_type_path = Path("C:/Code/COMPASS/Socio-economic/results/gis/avg_rural_per_grid.tif")

with rasterio.open(settlement_type_path) as src:
        settlement_type_grid = src.read(1, masked=True)
        settlement_type_grid_affine = src.transform
        settlement_type_grid_crs = src.crs

# Open original 1 km population rasters
with rasterio.open(population_raster_path_2019) as src_2019:
    region_proj = region.to_crs(src_2019.crs)
    pop_2019_1km, transform_pop_2019_1km = mask(src_2019, region_proj.geometry, crop=True)
    print("No-data value Population 2019:", src_2019.nodata)

with rasterio.open(population_raster_path_1990) as src_1990:
    region_proj = region.to_crs(src_1990.crs)
    pop_1990_1km, transform_pop_1990_1km = mask(src_1990, region_proj.geometry, crop=True)
    print("No-data value Population 1990:", src_1990.nodata)


# flood raster
F_flooding = sfincs_dir_F / "plot_output" / "floodmap.tif"
CF_flooding = sfincs_dir_CF / "plot_output" / "floodmap.tif"

# Flood model subgrid
sfincs_subgrid = join(sfincs_dir_F, "subgrid", "dep_subgrid.tif")


#%% Read flood data and background polygons
# --- Flood grid properties ---
with rasterio.open(sfincs_subgrid) as src:
    flood_grid_crs, flood_grid_transform, flood_grid_shape = src.crs, src.transform, (src.height, src.width)

# --- Setup region ---
region = region.to_crs(flood_grid_crs)
region_wsg84 = region.to_crs("EPSG:4326")
region_geom = [json.loads(region.to_json())["features"][0]["geometry"]]

# --- Read flood rasters ---
hmax_F_da = rxr.open_rasterio(F_flooding).squeeze("band", drop=True)  # if single-band
hmax_CF_da = rxr.open_rasterio(CF_flooding).squeeze("band", drop=True)  # if single-band
hmax_F = hmax_F_da.values
hmax_CF = hmax_CF_da.values
hmax_diff = hmax_F - hmax_CF

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
beira_utm = beira_district.to_crs(flood_grid_crs)

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

    # --- Validation printout ---
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
def pop_raster_to_gdf(pop_array, flood_array, settlement_type_array, transform, year, climate, export_df=True, export_path=None):
    print("Linking population raster to flood depth as DataFrame...")

    # Check if shapes match
    if pop_array.shape == flood_array.shape == settlement_type_array.shape:
        print("✔ Shapes match")
    else:
        print("✖ Shapes do NOT match!", pop_array.shape, flood_array.shape, settlement_type_array.shape)

    # Flatten arrays
    pop_flat = pop_array.ravel()
    flood_flat = flood_array.ravel()
    settlement_type_flat = settlement_type_array.ravel()

    # Mask zero-pop cells
    mask = (pop_flat > 0) & (flood_flat > 0)
    pop_vals = pop_flat[mask]
    flood_vals = flood_flat[mask]
    settlement_type_vals = settlement_type_flat[mask]

    # Pixel coordinates (centers)
    rows, cols = np.indices(pop_array.shape)
    xs, ys = transform * (cols.ravel()[mask] + 0.5, rows.ravel()[mask] + 0.5)

    df = pd.DataFrame({
        "population": pop_vals,
        "flood_depth": flood_vals,
        "settlement_type": settlement_type_vals,
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
    cells_flooded = block_count_threshold(flood_raster, factor=factor, threshold=0)
    pixels_high = block_count_threshold(flood_raster, factor=factor, threshold=1)
    pixels_higher = block_count_threshold(flood_raster, factor=factor, threshold=1.5)

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
            "pct_cells_flooded": (cells_flooded.flatten() / (factor**2)) * 100,
            "pct_cells_higher_1m": (pixels_high.flatten() / (factor**2)) * 100,
            "pct_cells_higher_1.5m": (pixels_higher.flatten() / (factor**2)) * 100
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
gdf_pop_2019_flood_depth_F  = pop_raster_to_gdf(pop_arrays[2019], hmax_F, settlement_type_grid, flood_grid_transform, year=2019, climate="F", export_df=True, export_path=export_path)
gdf_pop_2019_flood_depth_CF = pop_raster_to_gdf(pop_arrays[2019], hmax_CF, settlement_type_grid, flood_grid_transform, year=2019, climate="CF", export_df=True, export_path=export_path)
gdf_pop_1990_flood_depth_F  = pop_raster_to_gdf(pop_arrays[1990], hmax_F, settlement_type_grid, flood_grid_transform, year=1990, climate="F", export_df=True, export_path=export_path)
gdf_pop_1990_flood_depth_CF = pop_raster_to_gdf(pop_arrays[1990], hmax_CF, settlement_type_grid, flood_grid_transform, year=1990, climate="CF", export_df=True, export_path=export_path)

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
gdf_pop_2019_exposed_CF_uniform_coarse = aggregate_pop(pop_array_uniform_2019, hmax_CF, flood_grid_transform, flood_grid_crs, region_utm, background_utm)
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



#%% ======================================================================================= #
# ============= Compare estimated and reported totala nd exposed popuation ================ #
# ========================================================================================= #
# def get_extent(transform, width, height):
#     left = transform[2]
#     right = left + width * transform[0]
#     top = transform[5]
#     bottom = top + height * transform[4]
#     return [left, right, bottom, top]

# # Mask raster outside region (already cropped with rasterio.mask.mask)
# # pop_WP_masked = np.where(pop_WP_2020_coarse[0] == 0, np.nan, pop_WP_2020_coarse[0])
# pop_2019_1km_masked = np.where(pop_2019_1km[0] == 0, np.nan, pop_2019_1km[0])
# extent_2019_1km = get_extent(transform_pop_2019_1km, pop_2019_1km.shape[2], pop_2019_1km.shape[1])
# vmax_pop_2019_1km = np.percentile(pop_2019_1km_masked[~np.isnan(pop_2019_1km_masked)], 99.9)


# fig, axes = plt.subplots(1, 1, figsize=(6, 6))
# ax = axes

# # Increase extent
# xmin, ymin, xmax, ymax = region.total_bounds

# background.plot(ax=ax, color="#E0E0E0", zorder=0)
# bg_filtered.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
# region.boundary.plot(ax=ax, color='black', linewidth=1, zorder=2)
# districts_adm2_utm.boundary.plot(ax=ax, color='red', linewidth=0.5, zorder=5)
# districts_adm3_filtered.boundary.plot(ax=ax, color='blue', linewidth=0.5, zorder=4)
# ax.set_xlim(xmin - 0.2, xmax + 0.2)
# ax.set_ylim(ymin - 0.2, ymax + 0.2)

# for idx, row in districts_adm3_filtered.iterrows():
#     x = row.geometry.centroid.x
#     y = row.geometry.centroid.y
        
# ax.text(
#     x, y,
#     row["NAME_3"],
#     fontsize=7,
#     ha="center",
#     color="darkblue",
#     zorder=10,
#     fontweight="bold"
# )

# for idx, row in districts_adm2_utm.iterrows():
#     x = row.geometry.centroid.x
#     y = row.geometry.centroid.y
        
# ax.text(
#     x-0.15, y,
#     row["NAME_2"],
#     fontsize=7,
#     ha="center",
#     color="red",
#     zorder=10,
#     fontweight="bold"
# )

# # Plot Historical Exposure data 2020
# im = ax.imshow(pop_2019_1km_masked, cmap="Blues", extent=extent_2019_1km,  
#                     origin='upper', norm=PowerNorm(gamma=0.5, vmin=0, vmax=vmax_pop_2019_1km)) 
# cbar = plt.colorbar(im, ax=ax, shrink=0.8)
# cbar.set_label("Population (people per cell)")
# ax.set_title("Hist. Exposure 2020 \n[~1 km grid]")


#%% ============================================================================================ # 
# ================== Plot distribution of flood depth per exposed population =================== #
# ============================================================================================== #
def compute_cdf_and_bins(gdf, bins, depth_col="flood_depth", pop_col="population"):
    """Return population counts per flood depth bin."""
    flood = gdf[depth_col].values
    pop   = gdf[pop_col].values

    # Mask invalid
    mask = ~np.isnan(flood) & (pop > 0) # select only cells with flooding AND population
    flood, pop = flood[mask], pop[mask]

    # Binned population
    pop_by_depth = pd.Series(pop).groupby(pd.cut(flood, bins)).sum()

    return pop_by_depth
    
bins_fine = np.arange(0, 3.5 + 0.02, 0.01)
low_mask = bins_fine[:-1] < 0.5
mid_mask = (bins_fine[:-1] >= 0.5) & (bins_fine[:-1] < 1.5)
high_mask = bins_fine[:-1] >= 1.5
pop_2019_by_depth_F_fine  = compute_cdf_and_bins(gdf_pop_2019_exposed_F, bins_fine)
pop_2019_by_depth_CF_fine  = compute_cdf_and_bins(gdf_pop_2019_exposed_CF, bins_fine)
pop_1990_by_depth_F_fine   = compute_cdf_and_bins(gdf_pop_1990_exposed_F, bins_fine)
pop_1990_by_depth_CF_fine  = compute_cdf_and_bins(gdf_pop_1990_exposed_CF, bins_fine)

bins_coarse = np.arange(0, 3.5 + 0.2, 0.1)
pop_2019_by_depth_F_coarse = compute_cdf_and_bins(gdf_pop_2019_exposed_F, bins_coarse)
pop_2019_by_depth_CF_coarse = compute_cdf_and_bins(gdf_pop_2019_exposed_CF, bins_coarse)
pop_1990_by_depth_F_coarse  = compute_cdf_and_bins(gdf_pop_1990_exposed_F, bins_coarse)
pop_1990_by_depth_CF_coarse = compute_cdf_and_bins(gdf_pop_1990_exposed_CF, bins_coarse)

#%% --- settings for plotting ---
# Colours based on conceptual figure
colours = ['#00B050', '#1E2E57', "#28C2E9", '#9B59B6']

# Bin centers for plotting
bin_centers = bins_fine[:-1] + np.diff(bins_fine) / 2
bin_centers_coarse = bins_coarse[:-1] + np.diff(bins_coarse) / 2

# Plotting masks for different flood depth ranges (low, medium, high)
x_bg = np.linspace(0, 3.5, 500)  # example x array
low_mask_bg = x_bg < 0.5
mid_mask_bg = (x_bg >= 0.5) & (x_bg < 1.5)
high_mask_bg = x_bg >= 1.5



#%%
# ================== Plot attributable % of exposed population per flood depth =================== #
Change_per_flood_depth = pd.DataFrame({
    "Factual": pop_2019_by_depth_F_coarse,
    "CF_climate": pop_2019_by_depth_CF_coarse,
    "CF_population": pop_1990_by_depth_F_coarse,
    "CF_climate_population": pop_1990_by_depth_CF_coarse
})

Change_per_flood_depth["Rel_change_CF_climate"] = (Change_per_flood_depth["Factual"] - Change_per_flood_depth["CF_climate"]) / Change_per_flood_depth["Factual"] * 100
Change_per_flood_depth["Rel_change_CF_population"] = (Change_per_flood_depth["Factual"] - Change_per_flood_depth["CF_population"]) / Change_per_flood_depth["Factual"] * 100
Change_per_flood_depth["Rel_change_CF_climate_population"] = (Change_per_flood_depth["Factual"] - Change_per_flood_depth["CF_climate_population"]) / Change_per_flood_depth["Factual"] * 100
Change_per_flood_depth["Dominant_driver_cc"] = Change_per_flood_depth["Rel_change_CF_climate"] / Change_per_flood_depth["Rel_change_CF_population"]
Change_per_flood_depth["Dominant_driver_pc"] = Change_per_flood_depth["Rel_change_CF_population"] / Change_per_flood_depth["Rel_change_CF_climate"]

fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=300)
ax.plot(bin_centers_coarse, Change_per_flood_depth["Rel_change_CF_climate"], label="Climate change", color=colours[1])
ax.plot(bin_centers_coarse, Change_per_flood_depth["Rel_change_CF_population"], label="Population change", color=colours[2])
ax.plot(bin_centers_coarse, Change_per_flood_depth["Rel_change_CF_climate_population"], label="Climate & Population change", color=colours[3])

ymax = max(Change_per_flood_depth["Rel_change_CF_climate"].max(), Change_per_flood_depth["Rel_change_CF_population"].max(), Change_per_flood_depth["Rel_change_CF_climate_population"].max()) * 1.05
ymin = min(Change_per_flood_depth["Rel_change_CF_climate"].min(), Change_per_flood_depth["Rel_change_CF_population"].min(), Change_per_flood_depth["Rel_change_CF_climate_population"].min()) * 1.05
ax.fill_between(x_bg[low_mask_bg], ymin, ymax, color="#d9d9d9", alpha=0.3)
ax.fill_between(x_bg[mid_mask_bg], ymin, ymax, color="#b3b3b3", alpha=0.3)
ax.fill_between(x_bg[high_mask_bg], ymin, ymax, color="#808080", alpha=0.3)

ax.text(0.25, ymax*0.12, "Low", ha="center", fontweight="bold", color="#5C5C5C", fontsize=10)
ax.text(1.0, ymax*0.12, "Medium", ha="center", fontweight="bold", color="#5C5C5C", fontsize=10)
ax.text(2.5, ymax*0.12, "High", ha="center", fontweight="bold", color="#5C5C5C", fontsize=10)

ax.axhline(y=0, linestyle="--", color="black", alpha=0.7)
ax.set_xlabel("Flood depth (m)")  
ax.set_ylabel("Attributable exposed population (%)")
ax.set_title("Attributable population exposed by flood depth bins") 
ax.legend(loc="upper right")
ax.set_xlim(0, 3.5)
ax.set_ylim(ymin, ymax)
ax.grid(True, linestyle="--", alpha=0.5)


# ================== Plot dominant driver of attributable exposed population per flood depth =================== #
fig, ax = plt.subplots(1, 1, figsize=(8,5), dpi=300)
ax.plot(bin_centers_coarse, Change_per_flood_depth["Dominant_driver_cc"], color=colours[1])
ax.axhline(y=0, linestyle="--", color="black", alpha=0.7)
ax.set_xlabel("Flood depth (m)")  
ax.set_ylabel("Dominant driver of attributable exposed population \n(Climate change / Population change)")
ax.set_title("Dominant driver of attributable % of population exposed \n(Climate change / Population change)")
ax.set_xlim(0, 3.5)
ax.grid(True, linestyle="--", alpha=0.5)


fig, ax = plt.subplots(1, 1, figsize=(8,5), dpi=300)
ax.plot(bin_centers_coarse, Change_per_flood_depth["Dominant_driver_pc"], color=colours[2])
ax.axhline(y=0, linestyle="--", color="black", alpha=0.7)
ax.set_xlabel("Flood depth (m)")  
ax.set_ylabel("Dominant driver of attributable exposed population \n(Population change / Climate change)")
ax.set_title("Dominant driver of attributable % of population exposed \n(Population change / Climate change)")
ax.set_xlim(0, 3.5)
ax.grid(True, linestyle="--", alpha=0.5)

#%%
# =========================== Plotting exposed population per flood depth ============================ #
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



#%%
# Plotting absolute change in exposed population per water depth
fig, ax = plt.subplots(figsize=(8,5))

x = bin_centers
y_F = pop_2019_by_depth_F_fine.values
y_CF_clim = pop_2019_by_depth_CF_fine.values
y_CF_pop = pop_1990_by_depth_F_fine.values
y_CF_clim_pop = pop_1990_by_depth_CF_fine.values

# Possible smoothing of curves
# sigma = 1.2  
# y_F = gaussian_filter1d(y_F, sigma)
# y_CF_clim = gaussian_filter1d(y_CF_clim, sigma)
# y_CF_pop = gaussian_filter1d(y_CF_pop, sigma)
# y_CF_clim_pop = gaussian_filter1d(y_CF_clim_pop, sigma)

ymax = max(y_F.max(), y_CF_clim.max(), y_CF_pop.max(), y_CF_clim_pop.max()) * 1.05

ax.fill_between(x_bg[low_mask_bg], 0, ymax, color="#d9d9d9", alpha=0.3)
ax.fill_between(x_bg[mid_mask_bg], 0, ymax, color="#b3b3b3", alpha=0.3)
ax.fill_between(x_bg[high_mask_bg], 0, ymax, color="#808080", alpha=0.3)

ax.plot(x, y_F, label=f"Factual ({np.nansum(ra_exposed_pop_2019_F).astype(int):,.0f} people)", color=colours[0], linewidth=2)
ax.plot(x, y_CF_clim, label=f"Counterfactual Climate ({np.nansum(ra_exposed_pop_2019_CF).astype(int):,.0f} people)", color=colours[1], linewidth=1)
ax.plot(x, y_CF_pop, label=f"Counterfactual Population ({np.nansum(ra_exposed_pop_1990_F).astype(int):,.0f} people)", color=colours[2], linewidth=1)
ax.plot(x, y_CF_clim_pop, label=f"Counterfactual Climate & Population ({np.nansum(ra_exposed_pop_1990_CF).astype(int):,.0f} people)", color=colours[3], linewidth=1)

# For plotting background fills, we need to interpolate the y-values to match the finer x-grid
# x_fine = np.linspace(x.min(), x.max(), 500)
# y_F_fine = np.interp(x_fine, x, y_F)
# y_CF_clim_fine = np.interp(x_fine, x, y_CF_clim)
# y_CF_pop_fine = np.interp(x_fine, x, y_CF_pop)
# y_CF_clim_pop_fine = np.interp(x_fine, x, y_CF_clim_pop)
# ax.fill_between(x_fine, y_F_fine, y_CF_clim_fine, where=(y_CF_clim_fine > y_F_fine),
#                 alpha=0.5, color=colours[1])
# ax.fill_between(x_fine, y_CF_pop_fine, y_F_fine, where=(y_F_fine > y_CF_pop_fine),
#                 alpha=0.5, color=colours[2])
# ax.fill_between(x_fine, y_CF_clim_pop_fine, y_F_fine, where=(y_F_fine > y_CF_clim_pop_fine),
#                 alpha=0.5, color=colours[3])

# peaks_F, _ = find_peaks(y_F)

# if len(peaks_F) > 0:
#     idx = peaks_F[np.argmax(y_F[peaks_F])]  # highest peak
#     ax.annotate(f"{x[idx]:.2f} m",
#                 xy=(x[idx], y_F[idx]),
#                 xytext=(x[idx] + 0.2, y_F[idx] * 0.9),
#                 arrowprops=dict(arrowstyle="->"))


# Find peaks and sort by height
peaks_F, props_F = find_peaks(y_F, prominence=0.02, distance=5) 
sorted_peaks_F = peaks_F[np.argsort(props_F["prominences"])[::-1]]
peaks_CF_clim, props_CF_clim = find_peaks(y_CF_clim, prominence=0.02, distance=5)
sorted_peaks_CF_clim = peaks_CF_clim[np.argsort(props_CF_clim["prominences"])[::-1]]
peaks_CF_pop, props_CF_pop = find_peaks(y_CF_pop, prominence=0.02, distance=10)
sorted_peaks_CF_pop = peaks_CF_pop[np.argsort(props_CF_pop["prominences"])[::-1]]
# find second peak for CF_pop manually due to slowrise
mask = (x >= 1.0) & (x <= 1.5) # peak is located between 1 and 1.5 m flood depth
if np.any(mask):
    idx_peak_2nd_CF_pop = np.argmax(y_CF_pop[mask])
    idx_peak_2nd_CF_pop = np.where(mask)[0][idx_peak_2nd_CF_pop]

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

for i, idx in enumerate(sorted_peaks_CF_pop[:1]): 
    plt.annotate(f"{x[idx]:.2f} m",
                 xy=(x[idx]+0.01, y_CF_pop[idx]),
                 xytext=(x[idx]+0.15, y_CF_pop[idx]+0.03), color=colours[2], fontsize=8,
                 arrowprops=dict(arrowstyle="->", lw=0.8, color=colours[2], linestyle='--'))

plt.annotate(f"{x[idx_peak_2nd_CF_pop]:.2f} m",
                 xy=(x[idx_peak_2nd_CF_pop]+0.01, y_CF_pop[idx_peak_2nd_CF_pop]),
                 xytext=(x[idx_peak_2nd_CF_pop]+0.15, y_CF_pop[idx_peak_2nd_CF_pop]+0.03),
                 color=colours[2], fontsize=8,
                 arrowprops=dict(arrowstyle="->", lw=0.8, color=colours[2], linestyle='--'))
    
ax.text(0.25, ymax*0.03, "Low", ha="center", fontweight="bold", color="#5C5C5C", fontsize=10)
ax.text(1.0, ymax*0.03, "Medium", ha="center", fontweight="bold", color="#5C5C5C", fontsize=10)
ax.text(2.5, ymax*0.03, "High", ha="center", fontweight="bold", color="#5C5C5C", fontsize=10)

ax.set_xlabel("Flood depth (m)")
ax.set_ylabel("Exposed population")
ax.set_xlim(0.05, 3.5)
ax.set_ylim(0, ymax)
ax.grid(True, linestyle="--", alpha=0.5)
ax.legend()
# ax.set_title("Smoothed sigma=1.2")

#%%
# --- Plot absolute line plot COMPARED TO FACTUAL ---
fig, ax = plt.subplots(figsize=(8,5))

perct_attr_clim = (np.nansum(ra_exposed_pop_2019_F) - np.nansum(ra_exposed_pop_2019_CF)) / np.nansum(ra_exposed_pop_2019_F) * 100
perct_attr_pop = (np.nansum(ra_exposed_pop_2019_F) - np.nansum(ra_exposed_pop_1990_F)) / np.nansum(ra_exposed_pop_2019_F) * 100
perct_attr_clim_pop = (np.nansum(ra_exposed_pop_2019_F) - np.nansum(ra_exposed_pop_1990_CF)) / np.nansum(ra_exposed_pop_2019_F) * 100

diff_clim = (pop_2019_by_depth_F_fine.values - pop_2019_by_depth_CF_fine.values)
diff_pop = (pop_2019_by_depth_F_fine.values - pop_1990_by_depth_F_fine.values)
diff_clim_pop = (pop_2019_by_depth_F_fine.values - pop_1990_by_depth_CF_fine.values)

ax.plot(bin_centers, diff_clim, label=f"Climate change ({int(np.round(perct_attr_clim))} %)", color=colours[1])
ax.plot(bin_centers, diff_pop, label=f"Population change ({int(np.round(perct_attr_pop))} %)", color=colours[2])
ax.plot(bin_centers, diff_clim_pop, label=f"Climate & Population change ({int(np.round(perct_attr_clim_pop))} %)", color=colours[3])
ax.set_xlabel("Flood depth (m)")  
ax.set_ylabel("Absolute change in exposed population")
ax.axhline(y=0, linestyle="--", color="black", alpha=0.7)
ax.legend(loc='upper right', fontsize=9)
ax.set_xlim(0.05, 3.5)
ax.set_ylim(diff_clim.min()*1.1, diff_clim_pop.max()*1.1)
ax.grid(True, linestyle="--", alpha=0.5)

ax.fill_between(x_bg[low_mask_bg], diff_clim.min()*1.1, diff_clim_pop.max()*1.1, color="#d9d9d9", alpha=0.3)
ax.fill_between(x_bg[mid_mask_bg], diff_clim.min()*1.1, diff_clim_pop.max()*1.1, color="#b3b3b3", alpha=0.3)
ax.fill_between(x_bg[high_mask_bg], diff_clim.min()*1.1, diff_clim_pop.max()*1.1, color="#808080", alpha=0.3)
ax.text(0.25, diff_clim_pop.max()*0.1, "Low", ha="center", fontweight="bold", color="#5C5C5C", fontsize=10,
        bbox=dict(boxstyle="round",pad=0.15, facecolor="white", edgecolor="none", alpha=0.7))
ax.text(1.0, diff_clim_pop.max()*0.1, "Medium", ha="center", fontweight="bold", color="#5C5C5C", fontsize=10,
        bbox=dict(boxstyle="round",pad=0.15, facecolor="white", edgecolor="none", alpha=0.7))
ax.text(2.5, diff_clim_pop.max()*0.1, "High", ha="center", fontweight="bold", color="#5C5C5C", fontsize=10,
        bbox=dict(boxstyle="round",pad=0.15, facecolor="white", edgecolor="none", alpha=0.7))


#%%
# Plotting bar plot of absolute change in exposed population per water depth category
def compute_attr_per_flood_depth_mask(diff, baseline, mask):
    return np.nansum(diff[mask]) / np.nansum(baseline[mask]) * 100

# Absolute change in exposed population per flood depth category
low_vals_abs_diff = [diff_clim[low_mask].sum(), diff_pop[low_mask].sum(), diff_clim_pop[low_mask].sum()]
mid_vals_abs_diff = [diff_clim[mid_mask].sum(), diff_pop[mid_mask].sum(), diff_clim_pop[mid_mask].sum()]
high_vals_abs_diff = [diff_clim[high_mask].sum(), diff_pop[high_mask].sum(), diff_clim_pop[high_mask].sum()]

# Relative attributable change in exposed population per flood depth category
low_vals_attr  = [compute_attr_per_flood_depth_mask(diff_clim, pop_2019_by_depth_F_fine, low_mask),
                  compute_attr_per_flood_depth_mask(diff_pop, pop_2019_by_depth_F_fine, low_mask),
                  compute_attr_per_flood_depth_mask(diff_clim_pop, pop_2019_by_depth_F_fine, low_mask)]
mid_vals_attr  = [compute_attr_per_flood_depth_mask(diff_clim, pop_2019_by_depth_F_fine, mid_mask),
                  compute_attr_per_flood_depth_mask(diff_pop, pop_2019_by_depth_F_fine, mid_mask),
                  compute_attr_per_flood_depth_mask(diff_clim_pop, pop_2019_by_depth_F_fine, mid_mask)]
high_vals_attr = [compute_attr_per_flood_depth_mask(diff_clim, pop_2019_by_depth_F_fine, high_mask), 
                  compute_attr_per_flood_depth_mask(diff_pop, pop_2019_by_depth_F_fine, high_mask), 
                  compute_attr_per_flood_depth_mask(diff_clim_pop, pop_2019_by_depth_F_fine, high_mask)]

data_abs_diff = np.array([low_vals_abs_diff, mid_vals_abs_diff, high_vals_abs_diff])
data_attr = np.array([low_vals_attr, mid_vals_attr, high_vals_attr])


fig, axes = plt.subplots(1, 2, figsize=(10,5), sharex=True)

bar_width = 0.25
x_pos = np.arange(3)  # Low, Medium, High

labels = ["Low", "Medium", "High"]
subplot_labels_2 = ["(a)", "(b)"]

axes[0].bar(x_pos - bar_width, data_abs_diff[:,0], width=bar_width, 
       label=f"Climate change ({int(np.round(perct_attr_clim))} %)", 
       color=colours[1])
axes[0].bar(x_pos, data_abs_diff[:,1], width=bar_width, 
       label=f"Population change ({int(np.round(perct_attr_pop))} %)", 
       color=colours[2])
axes[0].bar(x_pos + bar_width, data_abs_diff[:,2], width=bar_width, 
       label=f"Climate & Population change ({int(np.round(perct_attr_clim_pop))} %)", 
       color=colours[3])

axes[1].bar(x_pos - bar_width, data_attr[:,0], width=bar_width, 
       label=f"Climate change ({int(np.round(perct_attr_clim))} %)", 
       color=colours[1])
axes[1].bar(x_pos, data_attr[:,1], width=bar_width, 
       label=f"Population change ({int(np.round(perct_attr_pop))} %)", 
       color=colours[2])
axes[1].bar(x_pos + bar_width, data_attr[:,2], width=bar_width, 
       label=f"Climate & population change ({int(np.round(perct_attr_clim_pop))} %)", 
       color=colours[3])

for i, ax in enumerate(axes):
    ax.axhline(0, linestyle="--", color="black", alpha=0.7)
    ax.set_axisbelow(True)
    ax.grid(True, axis='y', linestyle="--", alpha=0.5)
    ax.set_xlabel("Flood depth", fontsize=9)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontdict={'fontweight': 'bold', 'color': '#5C5C5C'})
    ax.text(0, 1.02, subplot_labels_2[i],
            transform=ax.transAxes,
            fontsize=10, fontweight="bold",
            va="bottom", ha="left")
    

axes[0].set_ylabel("Absolute change in exposed population (people)", fontsize=9)
axes[1].set_ylabel("Attributable exposed population (%)", fontsize=9)
axes[1].legend(fontsize=9, loc='upper right', bbox_to_anchor=(1, 1.2))

plt.tight_layout()
plt.show()

#%%
fig, ax = plt.subplots(1, 1, figsize=(11,5), sharex=True, dpi=300)

bar_width = 0.25
x_pos = np.arange(3)  # Low, Medium, High

labels = ["Low \n(<0.5 m)", "Medium \n(0.5–1.5 m)", "High \n(>1.5 m)"]

ax.bar(x_pos - bar_width, data_abs_diff[:,0], width=bar_width, 
       label=f"Climate change ({int(np.round(perct_attr_clim))} %)", 
       color=colours[1])
ax.bar(x_pos, data_abs_diff[:,1], width=bar_width, 
       label=f"Population change ({int(np.round(perct_attr_pop))} %)", 
       color=colours[2])
ax.bar(x_pos + bar_width, data_abs_diff[:,2], width=bar_width, 
       label=f"Climate & Population change ({int(np.round(perct_attr_clim_pop))} %)", 
       color=colours[3])

ax.axhline(0, linestyle="--", color="black", alpha=0.7)
ax.set_axisbelow(True)
ax.grid(True, axis='y', linestyle="--", alpha=0.5)
ax.set_xlabel("Flood depth", fontsize=9)
ax.set_xticks(x_pos)
ax.set_xticklabels(labels, fontdict={'fontweight': 'bold', 'color': '#5C5C5C'})
    
ax.set_ylabel("Absolute change in exposed population (people)", fontsize=9)
ax.legend(fontsize=9, loc='upper right', bbox_to_anchor=(1.5, 1))

plt.tight_layout()
plt.show()


#%%
fig, ax = plt.subplots(1, 1, figsize=(12,5), sharex=True, dpi=300)

bar_width = 0.25
x_pos = np.arange(3)  # Low, Medium, High

labels = ["Low \n(<0.5 m)", "Medium \n(0.5–1.5 m)", "High \n(>1.5 m)"]

ax.bar(x_pos - bar_width, data_attr[:,0], width=bar_width, 
       label=f"Climate change ({int(np.round(perct_attr_clim))} %)", 
       color=colours[1], edgecolor='grey', linewidth=0.3)
ax.bar(x_pos, data_attr[:,1], width=bar_width, 
       label=f"Population change ({int(np.round(perct_attr_pop))} %)", 
       color=colours[2], edgecolor='grey', linewidth=0.3)
ax.bar(x_pos + bar_width, data_attr[:,2], width=bar_width, 
       label=f"Climate & population change ({int(np.round(perct_attr_clim_pop))} %)", 
       color=colours[3], edgecolor='grey', linewidth=0.3)


ax.axhline(0, linestyle="-", color="black", alpha=0.7)
ax.set_axisbelow(True)
ax.grid(True, axis='y', linestyle="--", alpha=0.5)
ax.set_xlabel("Flood depth", fontsize=12)
ax.set_xticks(x_pos)
ax.set_xticklabels(labels, fontdict={'fontsize': 12, 'fontweight': 'bold', 'color': '#5C5C5C'})
    
ax.set_ylabel("Attributable exposed population (%)", fontsize=12)
ax.legend(fontsize=11, loc='upper right', bbox_to_anchor=(1.6, 1), frameon=False)

plt.tight_layout()
plt.show()


#%% ---------------------------------------------------------- #
# SUMMARY TABLE OF POPULATION EXPOSED PER FLOOD DEPTH CATEGORY #
# ------------------------------------------------------------ #
depth_bins = {"0-0.5 m": (0, 0.5), "0.5-1.5 m": (0.5, 1.5), ">1.5 m": (1.5, np.inf)}

scenarios = {
    "F": (gdf_pop_2019_exposed_F, 2019),
    "CF_clim": (gdf_pop_2019_exposed_CF, 2019),
    "CF_pop": (gdf_pop_1990_exposed_F, 1990),
    "CF_clim_pop": (gdf_pop_1990_exposed_CF, 1990)}

total_population = {
    2019: np.nansum(pop_arrays[2019]),
    1990: np.nansum(pop_arrays[1990])}

# -----------------
# COMPUTE ABSOLUTE & RELATIVE VALUES PER SCENARIO & DEPTH BIN
# -----------------
results_abs = {scenario: [] for scenario in scenarios}
results_rel = {scenario: [] for scenario in scenarios}

for depth_label, (dmin, dmax) in depth_bins.items():
    for scenario_name, (gdf, year) in scenarios.items():        
        mask = (gdf["flood_depth"] > dmin) & (gdf["flood_depth"] <= dmax)
        value = np.nansum(gdf.loc[mask, "population"])        
        results_abs[scenario_name].append(value)
        results_rel[scenario_name].append(value / total_population[year] * 100)

# -----------------
# ATTRIBUTABLE VALUES (relative to F)
# -----------------
attr_abs = {}
attr_pct = {}

for scenario in ["CF_pop", "CF_clim", "CF_clim_pop"]:
    attr_abs[scenario] = [results_abs["F"][i] - results_abs[scenario][i] for i in range(len(depth_bins))]
    attr_pct[scenario] = [
        (attr_abs[scenario][i] / results_abs["F"][i] * 100)
        if results_abs["F"][i] != 0 else np.nan
        for i in range(len(depth_bins))]

# -----------------
# BUILD TABLE
# -----------------
rows = []
labels = ["# affected", "% of total", "# attributable", "% attributed"]

for i, depth_label in enumerate(depth_bins.keys()):
    # Absolute
    rows.append({
        "Depth Range": depth_label,
        "Label": labels[0],
        **{s: results_abs[s][i] for s in scenarios}})

    # Relative
    rows.append({
        "Depth Range": depth_label,
        "Label": labels[1],
        **{s: results_rel[s][i] for s in scenarios}})

    # Attributable absolute
    rows.append({
        "Depth Range": depth_label,
        "Label": labels[2],
        "F": "-",
        **{s: attr_abs[s][i] for s in attr_abs}})

    # Attributable percent
    rows.append({
        "Depth Range": depth_label,
        "Label": labels[3],
        "F": "-",
        **{s: attr_pct[s][i] for s in attr_pct}})

# Add total population row
rows.append({
    "Depth Range": "Total population",
    "Label": "-",
    **{"F": total_population[2019],
       "CF_pop": total_population[1990],
       "CF_clim": total_population[2019],
       "CF_clim_pop": total_population[1990]}})

table_exposed_pop_per_flood_depth = pd.DataFrame(rows)

table_exposed_pop_per_flood_depth.to_csv("c:/Code/Paper_2/Output/flood_depth_impact_summary_table.csv",
                                         index=False)

print(table_exposed_pop_per_flood_depth)




#%%
# Fig 2
import matplotlib.patheffects as path_effects
print("Plotting attributable exposed population (three drivers)")

# --- Compute differences ---
gdf_CF_pop = gdf_pop_1990_exposed_F_coarse.copy()
gdf_CF_clim = gdf_pop_2019_exposed_CF_coarse.copy()
gdf_CF_clim_pop = gdf_pop_1990_exposed_CF_coarse.copy()

total_factual = gdf_pop_2019_exposed_F_coarse["exposed_population"].sum()
gdf_CF_clim["diff"] = (gdf_pop_2019_exposed_F_coarse["exposed_population"] - gdf_CF_clim["exposed_population"])
gdf_CF_pop["diff"] = (gdf_pop_2019_exposed_F_coarse["exposed_population"] - gdf_CF_pop["exposed_population"])
gdf_CF_clim_pop["diff"] = (gdf_pop_2019_exposed_F_coarse["exposed_population"] - gdf_CF_clim_pop["exposed_population"])

datasets = [
    ("Climate change", gdf_CF_clim),
    ("Population change", gdf_CF_pop),
    ("Climate & population change", gdf_CF_clim_pop)]

# --- Figure ---
fig, axes = plt.subplots(1, 3, figsize=(11, 5), dpi=300, constrained_layout=True,
                         subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Shared normalization
vmax = max(ds["diff"].max() for _, ds in datasets)
norm_diff = PowerNorm(gamma=0.5, vmin=0, vmax=vmax)

subplot_labels = ["(a)", "(b)", "(c)"]

# --- Plot loop ---
for i, (title, gdf) in enumerate(datasets):
    ax = axes[i]

    # Plot zero/negative as white
    gdf[gdf["diff"] <= 0].plot(ax=ax, color="white", edgecolor="grey", linewidth=0.2, zorder=1)

    # Plot positive difference
    gdf[gdf["diff"] > 0].plot(column="diff", cmap=plt.cm.Reds, norm=norm_diff, edgecolor="grey",
                              linewidth=0.2, ax=ax, legend=False, zorder=2, rasterized=True)

    # Region + background
    background_utm.plot(ax=ax, color="#E0E0E0", zorder=0)
    region_utm.boundary.plot(ax=ax, edgecolor="black", linewidth=0.3)
    ax.set_extent(flood_extent, crs=ccrs.UTM(36, southern_hemisphere=True))
    beira_utm.boundary.plot(ax=ax, edgecolor="black", linewidth=0.5, zorder=3, alpha=0.7)

    # Plot city and river locations and names
    ax.plot(34.862, -19.833, marker='o', color='black', markersize=3, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
    text = ax.text(34.852, -19.89, "Beira", transform=ccrs.PlateCarree(), fontsize=8, color='black', zorder=5)
    text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
    ax.text(34.975, -19.66, "Beira District", fontsize=8, ha="center", va="center", style="italic", 
            transform = ccrs.PlateCarree(), zorder=4, color="#5C5C5C")
    
    # Buzi River marker and label
    ax.plot(34.43, -19.89, marker='o', color='black', markersize=3, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
    text2 = ax.text(34.44, -19.87, "Buzi River", transform=ccrs.PlateCarree(),
                    fontsize=8, color='black', zorder=5)
    text2.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
    # Pungwe River marker and label
    ax.plot(34.543, -19.545, marker='o', color='black', markersize=3, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
    text3 = ax.text(34.554, -19.52, "Pungwe River", transform=ccrs.PlateCarree(),
                    fontsize=8, color='black', zorder=5)
    text3.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])

    # Gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color="gray", alpha=0.5, linestyle="--")
    gl.right_labels = False
    gl.top_labels = False
    if i != 0:
        gl.left_labels = False
    ax.set_title(title, fontsize=10)
    ax.text(0, 1.02, subplot_labels[i],
            transform=ax.transAxes,
            fontsize=10, fontweight="bold",
            va="bottom", ha="left")
    
    rel_change = (gdf["diff"].sum() / total_factual * 100)
    ax.text(
        0.98, 0.98,
        f"~{round(gdf['diff'].sum(), -3):,.0f} people\n({rel_change:.0f} %)",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.8)
    )

# --- Shared colorbar ---
sm = ScalarMappable(cmap=plt.cm.Reds, norm=norm_diff)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes, orientation="vertical", shrink=0.6, pad=0.02)
cbar.set_label("Attributable exposed population [# people]", fontsize=10)
cbar.ax.tick_params(labelsize=9)


# --- Stats of exposed population change in Beira district ---
# Select only cells that intersect Beira District
beira_cells_clim = gpd.overlay(gdf_CF_clim, beira_utm, how="intersection")
beira_cells_pop = gpd.overlay(gdf_CF_pop, beira_utm, how="intersection")
beira_cells_clim_pop = gpd.overlay(gdf_CF_clim_pop, beira_utm, how="intersection")

# Sum exposed population
print(f"Total exposed population change in Beira District (CF Clim): {round(beira_cells_clim['diff'].sum(), -3):,.0f}")
print(f"Total exposed population change in Beira District (CF Pop): {round(beira_cells_pop['diff'].sum(), -3):,.0f}")
print(f"Total exposed population change in Beira District (CF Clim & Pop): {round(beira_cells_clim_pop['diff'].sum(), -3):,.0f}")



# %%
print("Plotting attributable exposed population (three drivers)")

# --- Compute differences ---
gdf_CF_pop = gdf_pop_1990_exposed_F_coarse.copy()
gdf_CF_clim = gdf_pop_2019_exposed_CF_coarse.copy()
gdf_CF_clim_pop = gdf_pop_1990_exposed_CF_coarse.copy()

total_factual = gdf_pop_2019_exposed_F_coarse["exposed_population"].sum()
gdf_CF_clim["pct_attr"] = (gdf_pop_2019_exposed_F_coarse["exposed_population"] - gdf_CF_clim["exposed_population"]) / gdf_pop_2019_exposed_F_coarse["exposed_population"] * 100
gdf_CF_pop["pct_attr"] = (gdf_pop_2019_exposed_F_coarse["exposed_population"] - gdf_CF_pop["exposed_population"]) / gdf_pop_2019_exposed_F_coarse["exposed_population"] * 100
gdf_CF_clim_pop["pct_attr"] = (gdf_pop_2019_exposed_F_coarse["exposed_population"] - gdf_CF_clim_pop["exposed_population"]) / gdf_pop_2019_exposed_F_coarse["exposed_population"] * 100

datasets = [
    ("Climate change", gdf_CF_clim),
    ("Population change", gdf_CF_pop),
    ("Climate & population change", gdf_CF_clim_pop)]

# --- Figure ---
fig, axes = plt.subplots(1, 3, figsize=(11, 5), dpi=300, constrained_layout=True,
                         subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Shared normalization
vmax = max(ds["pct_attr"].max() for _, ds in datasets)
vmin = min(ds["pct_attr"].min() for _, ds in datasets)
norm_diff = TwoSlopeNorm(vmin=-100, vcenter=0, vmax=100)
cmap_diff = plt.cm.RdBu_r
subplot_labels = ["(a)", "(b)", "(c)"]

# --- Plot loop ---
for i, (title, gdf) in enumerate(datasets):
    ax = axes[i]

    # Plot zero/negative as white
    gdf[gdf["pct_attr"] <= 0].plot(ax=ax, color="white", edgecolor="grey", linewidth=0.2, zorder=1)

    # Plot positive difference
    gdf[gdf["pct_attr"] > 0].plot(column="pct_attr", cmap=cmap_diff, norm=norm_diff, edgecolor="grey",
                              linewidth=0.2, ax=ax, legend=False, zorder=2, rasterized=True)

    # Region + background
    background_utm.plot(ax=ax, color="#E0E0E0", zorder=0)
    region_utm.boundary.plot(ax=ax, edgecolor="black", linewidth=0.3)
    ax.set_extent(flood_extent, crs=ccrs.UTM(36, southern_hemisphere=True))
    beira_utm.boundary.plot(ax=ax, edgecolor="black", linewidth=0.5, zorder=3, alpha=0.7)

    # Plot city and river locations and names
    ax.plot(34.862, -19.833, marker='o', color='black', markersize=3, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
    text = ax.text(34.852, -19.89, "Beira", transform=ccrs.PlateCarree(), fontsize=8, color='black', zorder=5)
    text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
    ax.text(34.975, -19.66, "Beira District", fontsize=8, ha="center", va="center", style="italic", 
            transform = ccrs.PlateCarree(), zorder=4, color="#5C5C5C")
    
    # Buzi River marker and label
    ax.plot(34.43, -19.89, marker='o', color='black', markersize=3, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
    text2 = ax.text(34.44, -19.87, "Buzi River", transform=ccrs.PlateCarree(),
                    fontsize=8, color='black', zorder=5)
    text2.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
    # Pungwe River marker and label
    ax.plot(34.543, -19.545, marker='o', color='black', markersize=3, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
    text3 = ax.text(34.554, -19.52, "Pungwe River", transform=ccrs.PlateCarree(),
                    fontsize=8, color='black', zorder=5)
    text3.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])

    # Gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color="gray", alpha=0.5, linestyle="--")
    gl.right_labels = False
    gl.top_labels = False
    if i != 0:
        gl.left_labels = False
    ax.set_title(title, fontsize=10)
    ax.text(0, 1.02, subplot_labels[i],
            transform=ax.transAxes,
            fontsize=10, fontweight="bold",
            va="bottom", ha="left")
    
    # rel_change = (gdf["pct_attr"].sum() / total_factual * 100)
    # ax.text(
    #     0.98, 0.98,
    #     f"~{round(gdf['pct_attr'].sum(), -3):,.0f} people\n({rel_change:.0f} %)",
    #     transform=ax.transAxes,
    #     ha="right",
    #     va="top",
    #     fontsize=9,
    #     bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.8)
    # )

# --- Shared colorbar ---
sm = ScalarMappable(cmap=cmap_diff, norm=norm_diff)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes, orientation="vertical", shrink=0.6, pad=0.02)
cbar.set_label("Attributable exposed population [%]", fontsize=10)
cbar.ax.tick_params(labelsize=9)


# --- Stats of exposed population change in Beira district ---
# Select only cells that intersect Beira District
beira_cells_clim = gpd.overlay(gdf_CF_clim, beira_utm, how="intersection")
beira_cells_pop = gpd.overlay(gdf_CF_pop, beira_utm, how="intersection")
beira_cells_clim_pop = gpd.overlay(gdf_CF_clim_pop, beira_utm, how="intersection")

# Sum exposed population
# print(f"Total exposed population change in Beira District (CF Clim): {round(beira_cells_clim['diff'].sum(), -3):,.0f}")
# print(f"Total exposed population change in Beira District (CF Pop): {round(beira_cells_pop['diff'].sum(), -3):,.0f}")
# print(f"Total exposed population change in Beira District (CF Clim & Pop): {round(beira_cells_clim_pop['diff'].sum(), -3):,.0f}")


#%%
print("Plotting spatially aggregated RELATIVE exposed population damage for F, CF population and diff")

# Compute difference
gdf_pop_1990_exposed_F_coarse['population_diff'] = (gdf_pop_2019_exposed_F_coarse['relative_population'] -
                                                    gdf_pop_1990_exposed_F_coarse['relative_population'])
gdf_pop_2019_exposed_CF_coarse['population_diff'] = (gdf_pop_2019_exposed_F_coarse['relative_population'] -
                                                    gdf_pop_2019_exposed_CF_coarse['relative_population'])
gdf_pop_1990_exposed_CF_coarse['population_diff'] = (gdf_pop_2019_exposed_F_coarse['relative_population'] -
                                                    gdf_pop_1990_exposed_CF_coarse['relative_population'])

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(11, 6), dpi=300, sharey=True, constrained_layout=True,
                         subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Symmetric scale around zero for difference
diff1 = gdf_pop_1990_exposed_F_coarse["population_diff"].values
diff2 = gdf_pop_2019_exposed_CF_coarse["population_diff"].values
diff3 = gdf_pop_1990_exposed_CF_coarse["population_diff"].values

# Combine
all_diffs = np.concatenate([diff1, diff2, diff3])

# Remove NaNs
all_diffs = all_diffs[~np.isnan(all_diffs)]

# Symmetric absolute max
vabs = np.max(np.abs(all_diffs))

norm_diff = TwoSlopeNorm(vmin=-50, vcenter=0, vmax=50)
cmap_diff_pop = plt.get_cmap('RdBu_r')

# Plotting differences
gdf_pop_2019_exposed_CF_coarse.plot(column='population_diff', cmap=cmap_diff_pop, norm=norm_diff,
                                    edgecolor='grey', linewidth=0.2, ax=axes[0], legend=False, zorder=2,
                                    rasterized=True, missing_kwds={"color": "white"})
gdf_pop_1990_exposed_F_coarse.plot(column='population_diff', cmap=cmap_diff_pop, norm=norm_diff,
                                   edgecolor='grey', linewidth=0.2, ax=axes[1], legend=False, zorder=2,
                                   rasterized=True, missing_kwds={"color": "white"})
gdf_pop_1990_exposed_CF_coarse.plot(column='population_diff', cmap=cmap_diff_pop, norm=norm_diff,
                                    edgecolor='grey', linewidth=0.2, ax=axes[2], legend=False, zorder=2,
                                    rasterized=True, missing_kwds={"color": "white"})

# Formatting
subplot_labels = ['(a)', '(b)', '(c)']
scenarios = ["Climate change", "Population change", "Climate & Population change"]

for i, ax in enumerate(axes):
    region_utm.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    ax.set_extent(flood_extent, crs=ccrs.UTM(36, southern_hemisphere=True))

    gl = ax.gridlines(draw_labels=True, linewidth=0.4,
                      color='gray', alpha=0.5, linestyle='--')
    gl.right_labels = False
    gl.top_labels = False
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    if i != 0:
        gl.left_labels = False

    ax.text(0, 1.02, subplot_labels[i],
            transform=ax.transAxes,
            fontsize=10, fontweight='bold',
            va='bottom', ha='left')
    ax.set_title(f"{scenarios[i]}", fontsize=9)

# Difference colorbar
sm_diff = ScalarMappable(cmap=cmap_diff_pop, norm=norm_diff)
sm_diff.set_array([])
cbar2 = fig.colorbar(sm_diff, ax=axes[:3], shrink=0.6)
cbar2.set_label("Relative exposed population change [%]", fontsize=9)
cbar2.ax.tick_params(labelsize=8)


#%%
# ============================================================================= #
######## Calculate and plot attributable exposed population per district ########
# ============================================================================= #
# Fig 5
dataset_dict = {"F": {"exposed": ra_exposed_pop_2019_F,
                      "total": pop_arrays[2019]},
                "CF_clim": {"exposed": ra_exposed_pop_2019_CF,
                            "total": pop_arrays[2019]},
                "CF_pop": {"exposed": ra_exposed_pop_1990_F,
                           "total": pop_arrays[1990]},
                "CF_clim_pop": {"exposed": ra_exposed_pop_1990_CF,
                                "total": pop_arrays[1990]}}

results_exposed = {key: [] for key in dataset_dict.keys()}
results_total = {key: [] for key in dataset_dict.keys()}

for _, row in districts_adm3_filtered.iterrows():
    district_mask = features.geometry_mask([row.geometry], out_shape=ra_exposed_pop_2019_F.shape,
                                           transform=flood_grid_transform, invert=True)

    for key, data in dataset_dict.items():
        exposed_raster = data["exposed"]
        total_raster = data["total"]

        results_exposed[key].append(exposed_raster[district_mask].sum())
        results_total[key].append(total_raster[district_mask].sum())

# Assign to GeoDataFrame
for key in dataset_dict.keys():
    districts_adm3_filtered[f"pop_exposed_{key}"] = results_exposed[key]
    districts_adm3_filtered[f"pop_total_{key}"] = results_total[key]

    districts_adm3_filtered[f"relative_exposed_{key}"] = (districts_adm3_filtered[f"pop_exposed_{key}"] /
                                                          districts_adm3_filtered[f"pop_total_{key}"] * 100)

# Calculate attributable fractions
districts_adm3_filtered["attr_clim"] = ((districts_adm3_filtered["pop_exposed_F"] - 
                                        districts_adm3_filtered["pop_exposed_CF_clim"]) / 
                                        districts_adm3_filtered["pop_exposed_F"] * 100)

districts_adm3_filtered["attr_pop"] = ((districts_adm3_filtered["pop_exposed_F"] - 
                                       districts_adm3_filtered["pop_exposed_CF_pop"]) /
                                       districts_adm3_filtered["pop_exposed_F"] * 100)

districts_adm3_filtered["attr_clim_pop"] = ((districts_adm3_filtered["pop_exposed_F"] - 
                                             districts_adm3_filtered["pop_exposed_CF_clim_pop"]) / 
                                            districts_adm3_filtered["pop_exposed_F"] * 100)

# Combine all attributable values
all_attr = np.concatenate([districts_adm3_filtered["attr_clim"].values, 
                           districts_adm3_filtered["attr_pop"].values, 
                           districts_adm3_filtered["attr_clim_pop"].values])

all_attr = all_attr[~np.isnan(all_attr)]
vabs = np.nanpercentile(np.abs(all_attr), 99)
# norm_attr = Normalize(vmin=0, vmax=vabs)
norm_attr = PowerNorm(gamma=0.5, vmin=0, vmax=vabs)
cmap_attr = LinearSegmentedColormap.from_list("reds_from_rdbu", plt.cm.RdBu_r(np.linspace(0.5, 1, 256)))
columns = ["attr_clim", "attr_pop", "attr_clim_pop"]
titles = ["Climate change", "Population change", "Climate and population change"]

# Clip to region and remove permanent water bodies
districts_adm3_clipped = gpd.clip(districts_adm3_filtered, background_utm)
districts_adm3_clipped = gpd.clip(districts_adm3_clipped, region_utm)
districts_adm3_clipped = districts_adm3_clipped[~districts_adm3_clipped.geometry.is_empty]

fig, axes = plt.subplots(1, 3, figsize=(11, 6), dpi=300, sharey=True, constrained_layout=True, 
                         subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

for i, (ax, col, title) in enumerate(zip(axes, columns, titles)):
    districts_adm3_clipped.plot(column=col, cmap=cmap_attr, norm=norm_attr, edgecolor="grey", 
                                 linewidth=0.2, ax=ax, legend=False, rasterized=True, 
                                 missing_kwds={"color": "white"})

    region_utm.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)
    background_utm.plot(ax=ax, color='#E0E0E0', zorder=0)
    ax.set_extent(flood_extent, crs=ccrs.UTM(36, southern_hemisphere=True))
    ax.set_title(title, fontsize=9)

    gl = ax.gridlines(draw_labels=True, linewidth=0.4,
                      color='gray', alpha=0.5, linestyle='--')
    gl.right_labels = False
    gl.top_labels = False
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    if i != 0:
        gl.left_labels = False

    ax.text(0, 1.02, subplot_labels[i],
            transform=ax.transAxes,
            fontsize=10, fontweight='bold',
            va='bottom', ha='left')
    
    # Add values at the center of each district
    for _, row in districts_adm3_clipped.iterrows():
        # Get centroid coordinates
        x, y = row.geometry.centroid.x, row.geometry.centroid.y
        value = row[col]
        name = row["NAME_3"]

        # Plot text
        ax.text(x, y, f"{name}\n{int(round(value))} %", ha="center", va="center", fontsize=7, color="black",
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.3, boxstyle='round'))

# Shared colorbar
sm = ScalarMappable(cmap=cmap_attr, norm=norm_attr)
sm.set_array([])

cbar = fig.colorbar(sm, ax=axes, shrink=0.6)
cbar.set_label("Attributable relative exposed population [%]", fontsize=9)
cbar.ax.tick_params(labelsize=8)


#%%
# Plot effect of climate change alone
vmin = 0
vmax = np.nanpercentile(districts_adm3_clipped["attr_clim"].values, 99)
norm = Normalize(vmin=vmin, vmax=vmax)
reds = LinearSegmentedColormap.from_list("reds_only", plt.cm.RdBu_r(np.linspace(0.5, 1, 256)))

fig, ax = plt.subplots(1, 1, figsize=(8, 4), dpi=300, constrained_layout=True, 
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

districts_adm3_clipped.plot(column="attr_clim", cmap=reds, norm=norm, edgecolor="grey", 
                            linewidth=0.2, ax=ax, legend=False, rasterized=True, 
                            missing_kwds={"color": "white"})

region_utm.boundary.plot(ax=ax, edgecolor="black", linewidth=0.3, zorder=2)
background_utm.plot(ax=ax, color="#E0E0E0", zorder=0)

ax.set_extent(flood_extent, crs=ccrs.UTM(36, southern_hemisphere=True))
ax.set_title("Climate change", fontsize=11)
gl = ax.gridlines(draw_labels=True, linewidth=0.4,
                      color='gray', alpha=0.5, linestyle='--')
gl.right_labels = False
gl.top_labels = False
gl.xlabel_style = {'size': 8}
gl.ylabel_style = {'size': 8}

sm = ScalarMappable(cmap=reds, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, shrink=0.7)
cbar.set_label("Attributable relative exposed population [%]", fontsize=10)
cbar.ax.tick_params(labelsize=9)

# Add values at the center of each district
for _, row in districts_adm3_clipped.iterrows():
    # Get centroid coordinates
    x, y = row.geometry.centroid.x, row.geometry.centroid.y
    value = row["attr_clim"]
    name = row["NAME_3"]

    # Plot text
    ax.text(x, y, f"{name}\n{int(round(value))} %", ha="center", va="center", fontsize=7, color="black",
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.3, boxstyle='round'))







#%% #############################################################################################
#################################################################################################
#################################### SUPPLEMENTARY FIGURES ###################################### 
#################################################################################################
#################################################################################################

# ============================================================================================== # 
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

subplot_labels = ['(a)', '(b)', '(c)']
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
    ax.text(0, 1.02, subplot_labels[i], transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='bottom', ha='left')

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



#%% ============================================================================================ # 
# =============== Plot the change in exposed population > 1 m flood depth spatially ============ #
# ============================================================================================== #
# Plot average flood depth per exposed population cell
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True) # 10,6

# Plot average flood depth among exposed population
gdf_2019_F = gdf_pop_2019_exposed_F_coarse.copy()
gdf_2019_F.loc[gdf_2019_F["relative_population"] <= 1.5, "avg_flood_depth"] = np.nan

gdf_2019_CF = gdf_pop_2019_exposed_CF_coarse.copy()
gdf_2019_CF.loc[gdf_2019_CF["relative_population"] <= 1.5, "avg_flood_depth"] = np.nan

gdf_2019_CF['change_in_relative_population_>1m'] = gdf_2019_F['relative_population'] - gdf_2019_CF['relative_population']

# Define colormap: from white to #67CBE4
cmap = 'Blues'
norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2019_F['relative_population']))
cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#651F94"])
norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2019_CF['change_in_relative_population_>1m']))

plot = gdf_2019_F.plot(column="relative_population", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[0], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_CF.plot(column="relative_population", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[1], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_CF.plot(column="change_in_relative_population_>1m", cmap=cmap_change, norm=norm_change, linewidth=0.1, 
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
    cbar.set_label("Exposed population > 1.5 m flood depth")

sm = ScalarMappable(cmap=cmap_change, norm=norm_change)
sm._A = []  
cbar = plt.colorbar(sm, ax=axes[2], shrink=0.8)
cbar.set_label("Difference in relative xposed population > 1.5 m flood depth")

axes[0].set_title("Factual", fontsize=10)
axes[1].set_title("No Climate Change", fontsize=10)
axes[2].set_title("Factual - Counterfactual", fontsize=10)

axes[0].set_ylabel("y coordinate UTM zone 36S [m]")
axes[0].set_xlabel("x coordinate UTM zone 36S [m]")
axes[1].set_xlabel("x coordinate UTM zone 36S [m]")
axes[2].set_xlabel("x coordinate UTM zone 36S [m]")

fig.suptitle("Population exposed to > 1.5 m flood depth", fontsize=12)

plt.tight_layout()
plt.show()

#%% ============================================================================================ # 
# ===================== Plotting attributable % exposed to flood depth > 1.5 m ================= #
# ============================================================================================== #
# Plot average flood depth per exposed population cell
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True) # 10,6

# Plot average flood depth among exposed population
gdf_2019_F = gdf_pop_2019_exposed_F_coarse.copy()
gdf_2019_F.loc[gdf_2019_F["avg_flood_depth"] == 0, "exposed_population"] = 0

gdf_2019_CF = gdf_pop_2019_exposed_CF_coarse.copy()
gdf_2019_CF.loc[gdf_2019_CF["avg_flood_depth"] == 0, "exposed_population"] = 0

gdf_2019_CF['pect_change_in_exposed_population_>1.5m'] = (gdf_2019_F['exposed_population'] - gdf_2019_CF['exposed_population']) / gdf_2019_F['exposed_population'] * 100
# Define colormap: from white to #67CBE4
cmap = 'Blues'
norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2019_F['exposed_population']))
cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#651F94"])
norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2019_CF['pect_change_in_exposed_population_>1.5m']))

plot = gdf_2019_F.plot(column="exposed_population", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[0], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_CF.plot(column="exposed_population", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[1], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_CF.plot(column="pect_change_in_exposed_population_>1.5m", cmap=cmap_change, norm=norm_change, linewidth=0.1, 
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
    cbar.set_label("Exposed population")

sm = ScalarMappable(cmap=cmap_change, norm=norm_change)
sm._A = []  
cbar = plt.colorbar(sm, ax=axes[2], shrink=0.8)
cbar.set_label("Attributable exposed population")

axes[0].set_title("Factual", fontsize=10)
axes[1].set_title("No Climate Change", fontsize=10)
axes[2].set_title("(F - CF) / F * 100%", fontsize=10)

axes[0].set_ylabel("y coordinate UTM zone 36S [m]")
axes[0].set_xlabel("x coordinate UTM zone 36S [m]")
axes[1].set_xlabel("x coordinate UTM zone 36S [m]")
axes[2].set_xlabel("x coordinate UTM zone 36S [m]")

fig.suptitle("Total exposed population", fontsize=12)

plt.tight_layout()
plt.show()


#%% ============================================================================================ # 
# ============================= Plot % cells with flood depth > 1.5 m ============================ #
# ============================================================================================== #
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot average flood depth among exposed population
gdf_2019_F = gdf_pop_2019_exposed_F_coarse.copy()
gdf_2019_CF = gdf_pop_2019_exposed_CF_coarse.copy()

gdf_2019_CF['change_in_%more_1.5m'] = gdf_2019_F['pct_cells_higher_1.5m'] - gdf_2019_CF['pct_cells_higher_1.5m']

# Define colormap: from white to #67CBE4
cmap = 'Blues'
norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2019_F['pct_cells_higher_1.5m']))
cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#651F94"])
norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=(gdf_2019_CF['change_in_%more_1.5m']).quantile(0.99))

plot = gdf_2019_F.plot(column="pct_cells_higher_1.5m", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[0], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_CF.plot(column="pct_cells_higher_1.5m", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[1], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_CF.plot(column="change_in_%more_1.5m", cmap=cmap_change, norm=norm_change, linewidth=0.1, 
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
    cbar.set_label("% cells > 1.5 m flood depth")

sm = ScalarMappable(cmap=cmap_change, norm=norm_change)
sm._A = []  
cbar = plt.colorbar(sm, ax=axes[2], shrink=0.8)
cbar.set_label("Difference in % cells > 1.5 m flood depth")

axes[0].set_title("Factual", fontsize=10)
axes[1].set_title("No Climate Change", fontsize=10)
axes[2].set_title("Factual - Counterfactual", fontsize=10)

axes[0].set_ylabel("y coordinate UTM zone 36S [m]")
axes[0].set_xlabel("x coordinate UTM zone 36S [m]")
axes[1].set_xlabel("x coordinate UTM zone 36S [m]")
axes[2].set_xlabel("x coordinate UTM zone 36S [m]")

fig.suptitle("% cells > 1.5 m flood depth", fontsize=12)

plt.tight_layout()
plt.show()



#%% ============================================================================================ # 
# ================================ Plot % cells that are flooded  ============================== #
# ============================================================================================== #
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot average flood depth among exposed population
gdf_2019_F = gdf_pop_2019_exposed_F_coarse.copy()
gdf_2019_CF = gdf_pop_2019_exposed_CF_coarse.copy()

gdf_2019_CF['change_in_%_cells_flooded'] = gdf_2019_F['pct_cells_flooded'] - gdf_2019_CF['pct_cells_flooded']

# Define colormap: from white to #67CBE4
cmap = 'Blues'
norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2019_F['pct_cells_flooded']))
cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#651F94"])
norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=(gdf_2019_CF['change_in_%_cells_flooded']).quantile(0.99))

plot = gdf_2019_F.plot(column="pct_cells_flooded", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[0], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_CF.plot(column="pct_cells_flooded", cmap=cmap, norm=norm, linewidth=0.1, 
                edgecolor="grey", ax=axes[1], zorder=2,
                missing_kwds={"color": "none", "edgecolor": "none"})

plot = gdf_2019_CF.plot(column="change_in_%_cells_flooded", cmap=cmap_change, norm=norm_change, linewidth=0.1, 
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
    cbar.set_label("% cells flooded")

sm = ScalarMappable(cmap=cmap_change, norm=norm_change)
sm._A = []  
cbar = plt.colorbar(sm, ax=axes[2], shrink=0.8)
cbar.set_label("Difference in % cells flooded")

axes[0].set_title("Factual", fontsize=10)
axes[1].set_title("No Climate Change", fontsize=10)
axes[2].set_title("Factual - Counterfactual", fontsize=10)

axes[0].set_ylabel("y coordinate UTM zone 36S [m]")
axes[0].set_xlabel("x coordinate UTM zone 36S [m]")
axes[1].set_xlabel("x coordinate UTM zone 36S [m]")
axes[2].set_xlabel("x coordinate UTM zone 36S [m]")

fig.suptitle("% cells flooded", fontsize=12)

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
# ================================================================================== #
# SUPPLEMENTARY: Factual flood, population & exposed population + changes (6-panel)  #
# ================================================================================== #
def setup_map_axes(
    axes,
    region_utm,
    background_utm,
    flood_extent,
    subplot_labels=None,
    titles=None,
    show_left_labels_only=True,
    label_offset=(0, 1.02),
):
    """
    Apply standard map formatting to one or more Cartopy axes:
    region boundary, background, extent, gridlines, labels.
    """
    axes_arr = np.atleast_1d(axes)
    ncols = axes_arr.shape[-1] if axes_arr.ndim >= 2 else axes_arr.size
    nrows = axes_arr.shape[0] if axes_arr.ndim >= 2 else 1
    axes_flat = axes_arr.ravel()
    for i, ax in enumerate(axes_flat):
        background_utm.plot(ax=ax, color="#E0E0E0", zorder=0)
        region_utm.boundary.plot(ax=ax, edgecolor="black", linewidth=0.3)
        ax.set_extent(flood_extent, crs=ccrs.UTM(36, southern_hemisphere=True))

        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color="gray", alpha=0.5, linestyle="--")
        gl.right_labels = False
        gl.top_labels = False
        gl.xlabel_style = {"size": 9}
        gl.ylabel_style = {"size": 9}
        if show_left_labels_only and i % ncols != 0:
            gl.left_labels = False
        if i // ncols < nrows - 1:
            gl.bottom_labels = False
        

        if subplot_labels and i < len(subplot_labels):
            ax.text(
                label_offset[0],
                label_offset[1],
                subplot_labels[i],
                transform=ax.transAxes,
                fontsize=10,
                fontweight="bold",      
                va="bottom",
                ha="left",
            )
        if titles and i < len(titles):
            ax.set_title(titles[i], fontsize=10)

def plot_supp_factual_changes_overview(hmax_F, gdf_pop_2019_exposed_F_coarse, hmax_diff, 
                                       gdf_pop_1990_exposed_F_coarse, gdf_pop_1990_exposed_CF_coarse,
                                       region_utm, background_utm, flood_extent):
    # Data preparation for plotting
    gdf_2019 = gdf_pop_2019_exposed_F_coarse.copy()
    gdf_2019.loc[gdf_2019["total_population"] == 0] = np.nan
    gdf_1990 = gdf_pop_1990_exposed_F_coarse.copy()
    gdf_1990.loc[gdf_1990["total_population"] == 0] = np.nan
    gdf_2019["change_in_population"] = gdf_2019["total_population"] - gdf_1990["total_population"]
    gdf_attr = gdf_pop_2019_exposed_F_coarse.copy()
    gdf_attr["attr_exposed_pop"] = gdf_pop_2019_exposed_F_coarse["exposed_population"] - gdf_pop_1990_exposed_CF_coarse["exposed_population"]

    # colour maps and norms
    norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_2019["total_population"].max())
    cmap_pop = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#FC6F37"])
    norm_pop_exposed = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2019_exposed_F_coarse["exposed_population"].max())
    cmap_pop_exposed = mcolors.LinearSegmentedColormap.from_list("white_to_darkblue", ["#ffffff", "#67CBE4"])
    cmap_change = plt.cm.Reds
    norm_pop_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2019["change_in_population"]))
    norm_pop_exposed_change = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_1990_exposed_CF_coarse["attr_exposed_pop"].max())

    fig, axes = plt.subplots(2, 3, figsize=(16, 8), dpi=300, sharex=True, sharey=True, constrained_layout=True,
                             subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

    # Plot 1 - Factual fooding
    im_1 = axes[0,0].imshow(hmax_F, cmap="viridis", extent=flood_extent, origin="lower", vmin=0, vmax=3.5, zorder=2)

    # Plot 2 - Factual population
    gdf_2019.plot(column="total_population", cmap=cmap_pop, norm=norm_pop, linewidth=0.1,
                  edgecolor="grey", ax=axes[0,1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

    # Plot 3 - Factual exposed population
    gdf_exposed_F = gdf_pop_2019_exposed_F_coarse.copy()
    gdf_exposed_F.loc[gdf_exposed_F["exposed_population"] == 0, "exposed_population"] = np.nan
    gdf_exposed_F.plot(column="exposed_population", cmap=cmap_pop_exposed, edgecolor="grey", 
                       norm=norm_pop_exposed, linewidth=0.2, ax=axes[0,2], legend=False, 
                       zorder=2, rasterized=True,
                       missing_kwds={"color": "none", "edgecolor": "none"})
    
    # Plot 4 - Climate change
    im_2 = axes[1,0].imshow(hmax_diff, cmap=cmap_change, extent=flood_extent, origin="lower", vmin=0, vmax=0.5, zorder=2)

    # Plot 5 - Population change
    gdf_2019.plot(column="change_in_population", cmap=cmap_change, norm=norm_pop_change, linewidth=0.1,
                  edgecolor="grey", ax=axes[1,1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    
    # Plot 6 - Attributable exposed population
    gdf_attr.loc[gdf_attr["attr_exposed_pop"] <= 0, "attr_exposed_pop"] = np.nan
    gdf_attr.plot(column="attr_exposed_pop", cmap=cmap_change, edgecolor="grey", 
                  norm=norm_pop_exposed_change, linewidth=0.2, ax=axes[1,2], legend=False, 
                  zorder=2, rasterized=True,
                  missing_kwds={"color": "none", "edgecolor": "none"})

    setup_map_axes(axes, region_utm, background_utm, flood_extent,
                   subplot_labels=["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"],
                   titles=["Factual flooding", "Factual population", "Factual exposed population", 
                           "Climate change", "Population change", "Attributable exposed population"])

    # Colour bars top row
    plt.colorbar(im_1, ax=axes[0,0], shrink=0.8).set_label("Flood depth (m)")
    sm = ScalarMappable(cmap=cmap_pop, norm=norm_pop)
    sm._A = []
    fig.colorbar(sm, ax=axes[0,1], shrink=0.8).set_label("Aggregated population (# people)")
    sm = ScalarMappable(cmap=cmap_pop_exposed, norm=norm_pop_exposed)
    sm._A = []
    fig.colorbar(sm, ax=axes[0,2], shrink=0.8).set_label("Aggregated exposed population (# people)")

    # Colour bars for the bottom row
    plt.colorbar(im_2, ax=axes[1,0], shrink=0.8).set_label("Attributable flood depth (m)")
    sm = ScalarMappable(cmap=cmap_change, norm=norm_pop_change)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes[1,1], orientation="vertical", shrink=0.8)
    cbar.set_label("Change in population [# people]", fontsize=10)
    sm = ScalarMappable(cmap=cmap_change, norm=norm_pop_exposed_change)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes[1,2], orientation="vertical", shrink=0.8)
    cbar.set_label("Attributable exposed population [# people]", fontsize=10)

    return fig


fig = plot_supp_factual_changes_overview(hmax_F, gdf_pop_2019_exposed_F_coarse, hmax_diff, gdf_pop_1990_exposed_F_coarse, gdf_pop_1990_exposed_CF_coarse, region_utm, background_utm, flood_extent)
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
gdf_2019_F = gdf_pop_2019_exposed_F_coarse.copy()
gdf_2019_F.loc[gdf_2019_F["total_population"] == 0] = np.nan

gdf_1990_F = gdf_pop_1990_exposed_F_coarse.copy()
gdf_1990_F.loc[gdf_1990_F["total_population"] == 0] = np.nan

gdf_1990_F['change_in_population'] = gdf_2019_F['total_population'] - gdf_1990_F['total_population']

subplot_labels = ['(a)', '(b)', '(c)']
# Define colormap: from white to #67CBE4
cmap = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#FC6F37"])
norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_pop_2019_exposed_F_coarse['total_population']))  
cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#BD2A2A"])
norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_1990_F['change_in_population']))

plot = gdf_2019_F.plot(column="total_population", cmap=cmap, norm=norm ,
                       linewidth=0.1, edgecolor="grey", ax=axes[0], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

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
# Fig 3
print("Plotting spatially aggregated RELATIVE exposed population for F, and diff for CF clim and pop")

gdf_pop_2019_exposed_CF_coarse["rel_diff"] = (gdf_pop_2019_exposed_F_coarse["relative_population"] - 
                                              gdf_pop_2019_exposed_CF_coarse["relative_population"])

gdf_pop_1990_exposed_F_coarse["rel_diff"] = (gdf_pop_2019_exposed_F_coarse["relative_population"] - 
                                             gdf_pop_1990_exposed_F_coarse["relative_population"])


fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 5), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Create colormap normalization
vmax_diff_pop = (gdf_pop_1990_exposed_F_coarse["rel_diff"]).quantile(0.99)
vmin_diff_pop = (gdf_pop_1990_exposed_F_coarse["rel_diff"]).quantile(0.01)
vmax_diff_clim = (gdf_pop_2019_exposed_CF_coarse["rel_diff"]).quantile(0.99)
vmin_diff_clim = (gdf_pop_2019_exposed_CF_coarse["rel_diff"]).quantile(0.01)
norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2019_exposed_F_coarse["relative_population"].max())
norm_diff_pop = TwoSlopeNorm(vmin=vmin_diff_pop, vcenter=0, vmax=vmax_diff_pop)
norm_diff_clim = PowerNorm(gamma=0.5, vmin=vmin_diff_clim, vmax=vmax_diff_clim)
red_half = LinearSegmentedColormap.from_list("bwr_red", plt.cm.bwr(np.linspace(0.5, 1, 256)))

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['relative_population'] == 0].plot(ax=axes[0], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse['relative_population'] > 0].plot(column='relative_population', cmap='Blues',  edgecolor='grey',
                                                                               norm=norm_pop, linewidth=0.2, ax=axes[0], legend=False, 
                                                                               zorder=2, rasterized=True)

gdf_pop_2019_exposed_CF_coarse[gdf_pop_2019_exposed_CF_coarse['rel_diff'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_2019_exposed_CF_coarse.plot(column='rel_diff', cmap=red_half, norm=norm_diff_clim, edgecolor='grey', 
                                           linewidth=0.2, ax=axes[1], legend=False, zorder=2, 
                                           missing_kwds={"color": "white", "edgecolor": "none"}, rasterized=True)

gdf_pop_1990_exposed_F_coarse[gdf_pop_1990_exposed_F_coarse['rel_diff'] == 0].plot(ax=axes[2], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_pop_1990_exposed_F_coarse.plot(column='rel_diff', cmap='bwr', norm=norm_diff_pop, edgecolor='grey', 
                                          linewidth=0.2, ax=axes[2], legend=False, zorder=2, 
                                          missing_kwds={"color": "white", "edgecolor": "none"}, rasterized=True)

subplot_labels = ['(a)', '(b)', '(c)']

for i, ax in enumerate(axes):
    # Add model region
    region_utm.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)

    # Set extent 
    minx, miny, maxx, maxy = region_utm.bounds.minx.item(), region_utm.bounds.miny.item(), region_utm.bounds.maxx.item(), region_utm.bounds.maxy.item()
    ax.set_extent(flood_extent, crs=ccrs.UTM(36, southern_hemisphere=True))
    
    # Plot background
    mask_box = box(34.8, -20.3, 35.3, -19.9)  # minx, miny, maxx, maxy
    background_outside_box = background[~background.intersects(mask_box)] # removing errorneous lines outside model region
    background_outside_box.plot(ax=ax, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)
    background_outside_box.boundary.plot(ax=ax, color="#818181", linewidth=0.2, 
                                             transform=ccrs.PlateCarree(), zorder=1)
    
    # Add gridlines and format tick labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.right_labels = False
    gl.top_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}
    if i != 0: 
        gl.left_labels = False
    
    ax.text(-0.05, 1.02, subplot_labels[i], transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='bottom', ha='left')

# Colorbars
sm1 = ScalarMappable(cmap="Blues", norm=norm_pop)
sm1._A = []  # required for colorbar without passing data
cbar = fig.colorbar(sm1, ax=axes[0], orientation="vertical", shrink=0.5)
cbar.set_label("Aggregated relative exposed population [%]", fontsize=10)
cbar.ax.tick_params(labelsize=9)

sm2 = ScalarMappable(cmap=red_half, norm=norm_diff_clim)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, ax=axes[1], orientation="vertical", shrink=0.5)
cbar2.set_label("Attributable relative exposed population [%]", fontsize = 9)
cbar2.ax.tick_params(labelsize=8)

sm3 = ScalarMappable(cmap='bwr', norm=norm_diff_pop)
sm3.set_array([])
cbar3 = fig.colorbar(sm3, ax=axes[2], orientation="vertical", shrink=0.5)
cbar3.set_label("Attributable relative exposed population [%]", fontsize = 9)
cbar3.ax.tick_params(labelsize=8)

axes[0].set_title("Factual", fontsize=8)
axes[1].set_title("Factual - Counterfactual climate", fontsize=8)
axes[2].set_title("Factual - Counterfactual population", fontsize=8)


# %%

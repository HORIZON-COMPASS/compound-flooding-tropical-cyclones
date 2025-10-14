#%%
print("Loading packages...")
import os
import json
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from os.path import join
import rasterio
from rasterio import features
import geopandas as gpd
import warnings
warnings.filterwarnings('ignore')
import platform
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import rioxarray as rxr 
from hydromt import DataCatalog
import xarray as xr
import re
import glob
from tqdm import tqdm
import shutil
import re
from shapely.geometry import box
from numba import njit, prange
from concurrent.futures import ProcessPoolExecutor, as_completed


prefix = "p:/" if platform.system() == "Windows" else "/p/"

#%%
# ===== CONFIGURATION =====
print("Setting up paths and parameters...")
BASE_RUN_PATH = Path(os.path.join(prefix,"11210471-001-compass","03_Runs","sofala","Idai"))

#%%
# Input files
shapefile_fp = BASE_RUN_PATH / "sfincs" / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" / "gis" / "region.geojson"   # replace with your region shapefile
background = gpd.read_file(os.path.join(prefix, "11210471-001-compass","01_Data","sofala_geoms","sofala_region_background.geojson"))
region = gpd.read_file(shapefile_fp)

# Flood model subgrid
sfincs_subgrid = BASE_RUN_PATH / "sfincs" / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" / "subgrid" / "dep_subgrid.tif"

# --- Flood grid properties ---
with rasterio.open(sfincs_subgrid) as src:
    flood_grid_crs, flood_grid_transform, flood_grid_shape = src.crs, src.transform, (src.height, src.width)


#%%
####################################################################
############### Load World Pop Age & Sex Structures ################
####################################################################
# WorldPop data source: https://hub.worldpop.org/project/categories?id=8
worldpop_folder_2020_100m = os.path.join(prefix,"11210471-001-compass","01_Data","population_data", "Worldpop","moz_agesex_structures_2020_CN_100m_R2025A_v1")
path_moz_agesex_combined_2020_100m = os.path.join(worldpop_folder_2020_100m, "moz_agesex_structures_2020_CN_100m_R2025A_v1.nc")

#%%
print("Reprojecting WorldPop 100m age-sex rasters to flood grid...")

# %%
def reproject_and_redistribute_population_over_land(
    pop_data,  # can be a raster path or a rioxarray DataArray
    land_gdf,
    flood_crs,
    flood_transform,
    flood_shape,
    province_geom=None,
    region=None,
    districts=None,
    year=None,
):
    """Reproject and redistribute population DataArray or raster over land pixels."""

    print(f"▶ Processing {year} population data...")

    # --- 1️⃣ Load raster or use DataArray directly ---
    with rasterio.open(pop_data) as src:
        pop = src.read(1, masked=True)
        pop_affine = src.transform
        pop_crs = src.crs

        # --- 2️⃣ Clip or dissolve if needed ---
        if province_geom is not None:
            print("Clipping population data to province geometry...")
            province_geom = province_geom.to_crs(pop_crs)
            province_geom = [json.loads(province_geom.to_json())["features"][0]["geometry"]]
            pop, pop_affine = rasterio.mask.mask(pop, province_geom, crop=True, nodata=src.nodata)

        if districts is not None:
            print("Clipping population data to district geometry...")
            districts_single = districts.dissolve().reset_index(drop=True).to_crs(pop_crs)
            districts_geom = [districts_single.geometry.iloc[0].__geo_interface__]
            pop, pop_affine = rasterio.mask.mask(pop, districts_geom, crop=True, nodata=src.nodata)

        if region is not None:
            print("Clipping population data to region geometry...")
            region_wsg = region.to_crs(src.crs)
            region_geom = [json.loads(region_wsg.to_json())["features"][0]["geometry"]] 
            pop, pop_affine = rasterio.mask.mask(src, region_geom, crop=True, nodata=src.nodata)

        pop = pop.squeeze()

    # --- 3️⃣ Prepare land mask ---
    print("Creating land mask on flood grid...")
    land_mask = features.rasterize(
        [(geom, 1) for geom in land_gdf.geometry],
        out_shape=flood_shape,
        transform=flood_transform,
        fill=0,
        dtype=np.uint8,
    ).astype(bool)

    print("  ✔ Land mask created on flood grid.")

    # --- 4️⃣ Redistribute population ---
    pop_fine = np.zeros(flood_shape, dtype=np.float32)
    print("▶ Redistributing population to fine grid...")

    for row in tqdm(range(pop.shape[0]), desc="  Processing coarse cells"):
        for col in range(pop.shape[1]):
            pop_value = pop[row, col]
            if np.isnan(pop_value) or pop_value <= 0:
                continue

            x_min, y_max = pop_affine * (col, row)
            x_max, y_min = pop_affine * (col + 1, row + 1)
            coarse_bounds = box(x_min, y_min, x_max, y_max)

            # Transform to flood CRS
            coarse_bounds_flood = gpd.GeoSeries([coarse_bounds], crs=pop_crs).to_crs(flood_crs).iloc[0]

            coarse_mask = features.rasterize(
                [(coarse_bounds_flood, 1)],
                out_shape=flood_shape,
                transform=flood_transform,
                fill=0,
                dtype=np.uint8,
            ).astype(bool)

            valid_mask = coarse_mask & land_mask
            n_valid = valid_mask.sum()

            if n_valid > 0:
                pop_fine[valid_mask] += pop_value / n_valid

    # --- 5️⃣ Validation ---
    total_input = float(np.nansum(pop))
    total_output = float(pop_fine.sum())
    diff = abs(total_output - total_input)
    rel_diff = diff / total_input * 100
    print(f"  🔹 Input pop:  {total_input:,.0f}")
    print(f"  🔹 Output pop: {total_output:,.0f}")
    print(f"  🔹 Δ = {diff:,.2f} ({rel_diff:.4f} %)")

    if rel_diff > 0.01:
        print("  ⚠ WARNING: Population mismatch — check CRS or mask alignment!")

    return pop_fine


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
    chunk_size={'x': 1000, 'y': 1000}
):
    print(f"▶ Preparing WorldPop rasters from {folder_raster}")

    files = sorted(Path(folder_raster).glob("*.tif"))
    ds_vars = {}

    for file in tqdm(files, desc="Processing rasters"):
        name = file.stem
        if exclude_genders and any(g in name for g in exclude_genders):
            continue

        da = Path(file)

        pop_fine = reproject_and_redistribute_population_over_land(
            pop_data=da,
            land_gdf=land_gdf,
            flood_crs=flood_crs,
            flood_transform=flood_transform,
            flood_shape=flood_shape,
            region=region_gdf,
            year=year,
        )

        ds_vars[name] = (("y", "x"), pop_fine)

    ds = xr.Dataset(ds_vars)
    ds = ds.chunk(chunk_size)

    ds.attrs.update({
        'title': f'WorldPop Data reprojected to flood grid ({year})',
        'source': 'WorldPop',
        'year': year,
        'chunked': 'True',
        'chunk_size': str(chunk_size),
        'clipped_to_region': str(region_gdf is not None),
        'created_date': pd.Timestamp.now().isoformat(),
    })

    print(f"💾 Saving to {path_output_nc}")
    ds.to_netcdf(path_output_nc)

    return ds


#%%
# ===== RUN THE PROCESSING =====

# ds = prepare_worldpop_to_flood_grid(
#     folder_raster=worldpop_folder_2020_100m,
#     path_output_nc=path_moz_agesex_combined_2020_100m,
#     land_gdf=background,
#     flood_crs=flood_grid_crs,
#     flood_transform=flood_grid_transform,
#     flood_shape=flood_grid_shape,
#     region_gdf=region,
#     exclude_genders=['t'],
#     year=2020
# )

# %%
# ===== OPTIMIZED VERSION WITH NUMBA & VECTORIZATION =====

@njit(parallel=True)
def redistribute_population(pop_coarse, mapping, land_mask, pop_fine):
    """
    Fast redistribution using precomputed mask indices.
    
    pop_coarse: 2D array of coarse population
    mapping: 3D array of shape (coarse_rows, coarse_cols, max_pixels_per_cell) 
             containing flat indices in pop_fine for each coarse cell
    land_mask: flattened land mask
    pop_fine: flattened fine grid to accumulate population
    """
    for i in prange(pop_coarse.shape[0]):
        for j in range(pop_coarse.shape[1]):
            val = pop_coarse[i, j]
            if val <= 0 or np.isnan(val):
                continue
            cell_indices = mapping[i, j]
            n = len(cell_indices)
            if n > 0:
                portion = val / n
                for idx in cell_indices:
                    pop_fine[idx] += portion
    return pop_fine

def reproject_and_redistribute_population_over_land_vectorized(
    pop_data,  # raster path
    land_gdf,
    flood_crs,
    flood_transform,
    flood_shape,
    region=None,
    year=None
):
    """Vectorized and accelerated redistribution of population."""
    
    print(f"▶ Processing {year} population data (vectorized)...")

    # --- Load raster ---
    with rasterio.open(pop_data) as src:
        pop_coarse = src.read(1, masked=True)
        coarse_affine = src.transform
        coarse_crs = src.crs

    # --- Clip to region if provided ---
    if region is not None:
        region_proj = region.to_crs(coarse_crs)
        region_geom = [feature["geometry"] for feature in json.loads(region_proj.to_json())["features"]]
        with rasterio.open(pop_data) as src:
            pop_coarse, coarse_affine = rasterio.mask.mask(src, region_geom, crop=True, nodata=0)
        pop_coarse = pop_coarse.squeeze()

    # --- Precompute land mask on fine grid ---
    land_mask = features.rasterize(
        [(geom, 1) for geom in land_gdf.geometry],
        out_shape=flood_shape,
        transform=flood_transform,
        fill=0,
        dtype=np.uint8
    ).astype(bool)
    land_mask_flat = land_mask.ravel()
    pop_fine_flat = np.zeros(flood_shape[0] * flood_shape[1], dtype=np.float32)
    print("  ✔ Land mask created on flood grid.")

    # --- Create GeoDataFrame of coarse pixels ---
    coarse_rows, coarse_cols = pop_coarse.shape
    polys = []
    for i in range(coarse_rows):
        for j in range(coarse_cols):
            x0, y0 = coarse_affine * (j, i)
            x1, y1 = coarse_affine * (j+1, i+1)
            polys.append(box(x0, y0, x1, y1))
    gdf_coarse = gpd.GeoDataFrame({'row': np.repeat(np.arange(coarse_rows), coarse_cols),
                                   'col': np.tile(np.arange(coarse_cols), coarse_rows)},
                                  geometry=polys, crs=coarse_crs)
    
    # --- Transform all coarse polygons at once ---
    gdf_coarse = gdf_coarse.to_crs(flood_crs)
    
    # --- Rasterize coarse polygons to fine grid and store indices ---
    max_pixels_per_cell = 5000  # adjust depending on resolution
    mapping = np.empty((coarse_rows, coarse_cols), dtype=object)
    
    for idx, row in gdf_coarse.iterrows():
        mask = features.rasterize(
            [(row.geometry, 1)],
            out_shape=flood_shape,
            transform=flood_transform,
            fill=0,
            dtype=np.uint8
        ).astype(bool)
        valid_mask = np.flatnonzero(mask & land_mask_flat)
        mapping[row.row, row.col] = valid_mask

    # --- Flatten coarse pop and redistribute ---
    pop_fine_flat = redistribute_population(pop_coarse, mapping, land_mask_flat, pop_fine_flat)
    pop_fine = pop_fine_flat.reshape(flood_shape)

    # --- Validation ---
    total_input = float(np.nansum(pop_coarse))
    total_output = float(pop_fine.sum())
    diff = abs(total_output - total_input)
    rel_diff = diff / total_input * 100
    print(f"  🔹 Input pop:  {total_input:,.0f}")
    print(f"  🔹 Output pop: {total_output:,.0f}")
    print(f"  🔹 Δ = {diff:,.2f} ({rel_diff:.4f} %)")
    if rel_diff > 0.01:
        print("  ⚠ WARNING: Population mismatch — check CRS or mask alignment!")

    return pop_fine


def prepare_worldpop_to_flood_grid_slurm(
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
    n_workers=6  # adjust for number of CPUs on node or SLURM job
):
    """
    Process all WorldPop rasters and reproject to flood grid using vectorized redistribution.
    Designed for parallel execution (local or SLURM cluster).
    """
    print(f"▶ Preparing WorldPop rasters from {folder_raster} with {n_workers} workers")

    files = sorted(Path(folder_raster).glob("*.tif"))
    if exclude_genders:
        files = [f for f in files if not any(g in f.stem for g in exclude_genders)]

    if len(files) == 0:
        raise ValueError("No raster files found in folder or all excluded by filter.")

    ds_vars = {}

    # --- Helper to process one raster file ---
    def process_file(file):
        name = file.stem
        pop_fine = reproject_and_redistribute_population_over_land_vectorized(
            pop_data=file,
            land_gdf=land_gdf,
            flood_crs=flood_crs,
            flood_transform=flood_transform,
            flood_shape=flood_shape,
            region=region_gdf,
            year=year
        )
        return name, pop_fine

    # --- Parallel execution ---
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_file, f): f for f in files}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing rasters"):
            name, pop_fine = future.result()
            ds_vars[name] = (("y", "x"), pop_fine)

    # --- Combine all rasters into one xarray.Dataset ---
    ds = xr.Dataset(ds_vars)
    ds = ds.chunk(chunk_size)

    ds.attrs.update({
        'title': f'WorldPop Data reprojected to flood grid ({year})',
        'source': 'WorldPop',
        'year': year,
        'chunked': True,
        'chunk_size': str(chunk_size),
        'clipped_to_region': str(region_gdf is not None),
        'created_date': pd.Timestamp.now().isoformat(),
        'parallelized': True,
        'n_workers': n_workers
    })

    print(f"💾 Saving to {path_output_nc}")
    ds.to_netcdf(path_output_nc)

    return ds

#%%
test_folder = os.path.join(prefix,"11210471-001-compass","01_Data","population_data","Worldpop","moz_agesex_structures_2020_CN_100m_R2025A_v1","test_folder")
ds = prepare_worldpop_to_flood_grid_slurm(
    folder_raster=test_folder,
    path_output_nc=test_folder / "test_output.nc",
    land_gdf=background,
    flood_crs=flood_grid_crs,
    flood_transform=flood_grid_transform,
    flood_shape=flood_grid_shape,
    region_gdf=region,
    exclude_genders=['t'],
    year=2020
)
# %%

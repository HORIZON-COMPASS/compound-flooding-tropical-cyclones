#%%
import os
import numpy as np
import geopandas as gpd
import rasterio
from rasterio import features
from shapely.geometry import box
from tqdm import tqdm
import xarray as xr
import pandas as pd
from numba import njit
import platform
from pathlib import Path
import json
import numpy as np
from rasterio import mask



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


# -----------------------------
# Numba-accelerated redistribution
# -----------------------------
@njit
def redistribute_population_numba(pop_coarse, mapping, pop_fine_flat):
    for i in range(pop_coarse.shape[0]):
        for j in range(pop_coarse.shape[1]):
            val = pop_coarse[i, j]
            if val <= 0:
                continue
            idxs = mapping[i, j]
            n = len(idxs)
            if n > 0:
                for k in idxs:
                    pop_fine_flat[k] += val / n
    return pop_fine_flat

# -----------------------------
# Vectorized reproject + redistribute
# -----------------------------
def reproject_and_redistribute_population_over_land(
    pop_data,  # raster path
    land_gdf,
    flood_crs,
    flood_transform,
    flood_shape,
    region=None,
    year=None
):
    """Vectorized & accelerated redistribution of population onto flood grid."""

    print(f"▶ Processing {year} population data...")

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
            pop_coarse, coarse_affine = mask.mask(src, region_geom, crop=True, nodata=0)
        pop_coarse = pop_coarse.squeeze()

    # --- Prepare land mask on fine grid ---
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

    # --- Create all coarse cell polygons at once ---
    rows, cols = pop_coarse.shape
    polys = [box(*coarse_affine * (j, i), *coarse_affine * (j+1, i+1))
             for i in range(rows) for j in range(cols)]
    gdf_coarse = gpd.GeoDataFrame(
        {'row': np.repeat(np.arange(rows), cols),
         'col': np.tile(np.arange(cols), rows)},
        geometry=polys, crs=coarse_crs
    )

    # --- Transform all coarse polygons at once ---
    gdf_coarse = gdf_coarse.to_crs(flood_crs)

    # --- Rasterize each coarse cell and store indices of fine pixels ---
    mapping = np.empty((rows, cols), dtype=object)
    for idx, row in tqdm(gdf_coarse.iterrows(), total=len(gdf_coarse), desc="Rasterizing coarse cells"):
        mask = features.rasterize(
            [(row.geometry, 1)],
            out_shape=flood_shape,
            transform=flood_transform,
            fill=0,
            dtype=np.uint8
        ).astype(bool)
        mapping[row.row, row.col] = np.flatnonzero(mask & land_mask_flat)

    # --- Redistribute population with Numba ---
    pop_fine_flat = redistribute_population_numba(pop_coarse, mapping, pop_fine_flat)
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

# -----------------------------
# Wrapper for all rasters
# -----------------------------
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

        pop_fine = reproject_and_redistribute_population_over_land(
            pop_data=file,
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

test_folder = os.path.join(prefix,"11210471-001-compass","01_Data","population_data","Worldpop","moz_agesex_structures_2020_CN_100m_R2025A_v1","test_folder")
test_output = os.path.join(prefix,"11210471-001-compass","01_Data","population_data","Worldpop","moz_agesex_structures_2020_CN_100m_R2025A_v1","test_folder", "test_output.nc")
ds = prepare_worldpop_to_flood_grid(
    folder_raster=test_folder,
    path_output_nc=test_output,
    land_gdf=background,
    flood_crs=flood_grid_crs,
    flood_transform=flood_grid_transform,
    flood_shape=flood_grid_shape,
    region_gdf=region,
    exclude_genders=['t'],
    year=2020
)

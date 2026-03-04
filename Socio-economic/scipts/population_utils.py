"""
Core population processing functions:
- Reprojection and redistribution of population to flood grid
- Population-to-GeoDataFrame linking
- Spatial aggregation
- Raster I/O utilities
"""
import os
import json
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
import rasterio.mask
from rasterio import features
from rasterio.mask import mask
from affine import Affine
from shapely.geometry import box
from tqdm import tqdm

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Helper: get raster extent from transform
# ---------------------------------------------------------------------------
def get_extent(transform, width, height):
    """Return [left, right, top, bottom] from an affine transform."""
    left = transform[2]
    right = left + width * transform[0]
    top = transform[5]
    bottom = top + height * transform[4]
    return [left, right, top, bottom]


# ---------------------------------------------------------------------------
# Reproject and redistribute population onto flood grid
# ---------------------------------------------------------------------------
def reproject_and_redistribute_population_over_land(
    pop_path,
    land_gdf,
    flood_crs,
    flood_transform,
    flood_shape,
    province_geom=None,
    region=None,
    districts_adm3=None,
    districts_adm2=None,
    year=None,
    out_raster_path=None,
):
    """
    Read a coarse population raster, clip it, and redistribute population
    evenly across land pixels on the fine-resolution flood grid.

    Returns
    -------
    pop_fine : np.ndarray
        Redistributed population on the flood grid.
    pop_sofala, transform_sofala
        Population clipped to Sofala province.
    pop_districts_adm3, pop_affine_districts_adm3
        Population clipped to dissolved ADM3 districts.
    pop_districts_adm2, pop_affine_districts_adm2
        Population clipped to dissolved ADM2 districts.
    """
    print(f"▶ Loading {year} population raster...")
    with rasterio.open(pop_path) as src:
        pop = src.read(1, masked=True)
        pop_affine = src.transform
        pop_crs = src.crs

        # Clip to province
        province_proj = province_geom.to_crs(src.crs)
        province_json = [json.loads(province_proj.to_json())["features"][0]["geometry"]]
        pop_sofala, transform_sofala = mask(src, province_json, crop=True, nodata=src.nodata)

        # Dissolve districts and clip
        districts_adm3_single = districts_adm3.dissolve().reset_index(drop=True).to_crs(src.crs)
        districts_adm3_geom = [districts_adm3_single.geometry.iloc[0].__geo_interface__]

        districts_adm2_single = districts_adm2.dissolve().reset_index(drop=True).to_crs(src.crs)
        districts_adm2_geom = [districts_adm2_single.geometry.iloc[0].__geo_interface__]

        pop_districts_adm3, pop_affine_districts_adm3 = mask(
            src, districts_adm3_geom, crop=True, nodata=src.nodata
        )
        pop_districts_adm2, pop_affine_districts_adm2 = mask(
            src, districts_adm2_geom, crop=True, nodata=src.nodata
        )

        # Early return if pre-computed raster exists
        if out_raster_path is not None and os.path.exists(out_raster_path):
            print(f"▶ Loading existing raster from {out_raster_path}")
            with rasterio.open(out_raster_path) as cached:
                pop_fine = cached.read(1)
                return (
                    pop_fine,
                    pop_sofala,
                    transform_sofala,
                    pop_districts_adm3,
                    pop_affine_districts_adm3,
                    pop_districts_adm2,
                    pop_affine_districts_adm2,
                )

        # Clip to region
        if region is not None:
            region_wsg = region.to_crs(src.crs)
            region_geom = [json.loads(region_wsg.to_json())["features"][0]["geometry"]]
            pop, pop_affine = rasterio.mask.mask(src, region_geom, crop=True, nodata=src.nodata)

        pop = pop.squeeze()

    # --- Prepare fine grid ------------------------------------------------
    pop_fine = np.zeros(flood_shape, dtype=np.float32)

    land_mask = features.rasterize(
        [(geom, 1) for geom in land_gdf.geometry],
        out_shape=flood_shape,
        transform=flood_transform,
        fill=0,
        dtype=np.uint8,
    ).astype(bool)
    print("  ✔ Land mask created on flood grid.")

    # --- Redistribute each coarse pixel -----------------------------------
    print("▶ Redistributing population to fine grid...")
    for row in tqdm(range(pop.shape[0]), desc="  Processing coarse cells"):
        for col in range(pop.shape[1]):
            pop_value = pop[row, col]
            if np.isnan(pop_value) or pop_value <= 0:
                continue

            # Coarse pixel bounds (in coarse CRS)
            x_min, y_max = pop_affine * (col, row)
            x_max, y_min = pop_affine * (col + 1, row + 1)
            coarse_bounds = box(x_min, y_min, x_max, y_max)

            # Transform to flood CRS
            coarse_bounds_flood = gpd.GeoSeries([coarse_bounds], crs=pop_crs).to_crs(flood_crs).iloc[0]

            # Rasterize coarse cell onto flood grid
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
                P = int(round(pop_value))
                base = P // n_valid
                remainder = P % n_valid

                valid_indices = np.where(valid_mask)
                pop_fine[valid_indices] += base

                if remainder > 0:
                    perm = np.random.permutation(len(valid_indices[0]))
                    chosen = perm[:remainder]
                    pop_fine[valid_indices[0][chosen], valid_indices[1][chosen]] += 1

    # --- Validation -------------------------------------------------------
    total_input_pop = float(np.nansum(pop))
    total_output_pop = float(pop_fine.sum())
    diff = abs(total_output_pop - total_input_pop)
    rel_diff = diff / total_input_pop * 100

    print("  ✔ Redistribution done.")
    print(f"  🔹 Input population:  {total_input_pop:,.0f}")
    print(f"  🔹 Output population: {total_output_pop:,.0f}")
    print(f"  🔹 Difference:        {diff:,.2f} ({rel_diff:.4f} %)")

    if rel_diff > 0.01:
        print("  ⚠ WARNING: Population not perfectly preserved — check CRS or mask alignment!")

    # --- Optionally save --------------------------------------------------
    if out_raster_path is not None:
        _save_redistributed_raster(pop_fine, out_raster_path, flood_transform, flood_crs)

    return (
        pop_fine,
        pop_sofala,
        transform_sofala,
        pop_districts_adm3,
        pop_affine_districts_adm3,
        pop_districts_adm2,
        pop_affine_districts_adm2,
    )


def _save_redistributed_raster(pop_fine, out_path, flood_transform, flood_crs):
    """Write redistributed population to GeoTIFF, fixing y-axis if needed."""
    print(f"▶ Saving redistributed population raster to {out_path}")
    H, W = pop_fine.shape
    a, b, c, d, e, f = flood_transform

    if e > 0:
        print("  ⚠ Detected positive y-resolution in transform → fixing for QGIS")
        pop_fine_to_write = np.flipud(pop_fine)
        new_e = -abs(e)
        new_f = f + e * (H - 1)
        fixed_transform = Affine(a, b, c, d, new_e, new_f)
    else:
        pop_fine_to_write = pop_fine
        fixed_transform = flood_transform

    profile = {
        "driver": "GTiff",
        "dtype": rasterio.float32,
        "count": 1,
        "height": H,
        "width": W,
        "crs": flood_crs,
        "transform": fixed_transform,
        "compress": "deflate",
    }
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(pop_fine_to_write, 1)


# ---------------------------------------------------------------------------
# Link population raster to flood depth as DataFrame
# ---------------------------------------------------------------------------
def pop_raster_to_gdf(
    pop_array,
    flood_array,
    transform,
    year,
    climate,
    export_df=True,
    export_path=None,
):
    """
    Flatten population, flood-depth and settlement-type arrays into a
    DataFrame of pixel-level records (only cells with pop > 0 and flooding).
    """
    print("Linking population raster to flood depth as DataFrame...")

    if pop_array.shape == flood_array.shape:
        print("✔ Shapes match")
    else:
        print("✖ Shapes do NOT match!", pop_array.shape, flood_array.shape)

    pop_flat = pop_array.ravel()
    flood_flat = flood_array.ravel()

    valid = (pop_flat > 0) & (flood_flat > 0)
    pop_vals = pop_flat[valid]
    flood_vals = flood_flat[valid]

    rows, cols = np.indices(pop_array.shape)
    xs, ys = transform * (cols.ravel()[valid] + 0.5, rows.ravel()[valid] + 0.5)

    df = pd.DataFrame(
        {
            "population": pop_vals,
            "flood_depth": flood_vals,
            "x": xs,
            "y": ys,
        }
    )

    if export_df and export_path is not None:
        file_name = f"df_pop_{year}_{climate}.csv"
        out = os.path.join(export_path, file_name)
        df.to_csv(out, index=False)
        print(f"▶ Exported DataFrame to {out}")

    return df


# ---------------------------------------------------------------------------
# Aggregate population / flood to coarser grid
# ---------------------------------------------------------------------------
def aggregate_pop(
    total_pop_array,
    flood_raster,
    transform,
    crs,
    region=None,
    background=None,
    factor=100,
):
    """
    Block-aggregate population and flood statistics to a coarser grid,
    then clip to the region and redistribute proportionally by area.
    """
    print("Aggregating population raster to coarser grid polygons...")

    pixel_width = transform.a
    pixel_height = -transform.e
    cell_width = factor * pixel_width
    cell_height = factor * pixel_height

    def _block_sum(arr, f):
        nr, nc = arr.shape
        nr_c, nc_c = nr - nr % f, nc - nc % f
        return arr[:nr_c, :nc_c].reshape(nr_c // f, f, nc_c // f, f).sum(axis=(1, 3))

    def _block_mean(arr, f):
        nr, nc = arr.shape
        nr_c, nc_c = nr - nr % f, nc - nc % f
        return np.nanmean(arr[:nr_c, :nc_c].reshape(nr_c // f, f, nc_c // f, f), axis=(1, 3))

    def _block_count(arr, f, threshold=1):
        nr, nc = arr.shape
        nr_c, nc_c = nr - nr % f, nc - nc % f
        return (arr[:nr_c, :nc_c].reshape(nr_c // f, f, nc_c // f, f) > threshold).sum(axis=(1, 3))

    total_agg = _block_sum(total_pop_array, factor)
    exposed_agg = _block_sum(np.where(flood_raster > 0, total_pop_array, 0), factor)
    avg_flood = _block_mean(flood_raster, factor)
    cells_flooded = _block_count(flood_raster, factor, threshold=0)
    cells_high = _block_count(flood_raster, factor, threshold=1)
    cells_higher = _block_count(flood_raster, factor, threshold=1.5)

    # Build coarse grid polygons
    nrows_c, ncols_c = total_agg.shape
    x0, y0 = transform * (0, 0)
    x_coords = x0 + np.arange(ncols_c) * cell_width
    y_coords = y0 + np.arange(nrows_c) * -cell_height

    grid_cells = [
        box(x, y - cell_height, x + cell_width, y) for y in y_coords for x in x_coords
    ]

    gdf_grid = gpd.GeoDataFrame(
        {
            "total_population": total_agg.flatten(),
            "exposed_population": exposed_agg.flatten(),
            "avg_flood_depth": avg_flood.flatten(),
            "pct_cells_flooded": (cells_flooded.flatten() / (factor**2)) * 100,
            "pct_cells_higher_1m": (cells_high.flatten() / (factor**2)) * 100,
            "pct_cells_higher_1.5m": (cells_higher.flatten() / (factor**2)) * 100,
        },
        geometry=grid_cells,
        crs=crs,
    )

    # Clip to region/background
    region_bg = (
        gpd.overlay(region, background, how="intersection")
        if background is not None
        else region.copy()
    )

    gdf_grid = gdf_grid.reset_index(names="cell_id")
    intersections = gpd.overlay(gdf_grid, region_bg, how="intersection")
    intersections["intersect_area"] = intersections.geometry.area
    area_sum = intersections.groupby("cell_id")["intersect_area"].transform("sum")
    intersections["norm_fraction"] = intersections["intersect_area"] / area_sum

    for col in ["total_population", "exposed_population"]:
        intersections[col] = (
            intersections.groupby("cell_id")[col].transform("first")
            * intersections["norm_fraction"]
        )

    gdf_grid_masked = intersections.dissolve(by="cell_id", aggfunc="sum")
    gdf_grid_masked["geometry"] = intersections.dissolve(by="cell_id").geometry
    gdf_grid_masked["relative_population"] = (
        gdf_grid_masked["exposed_population"] / gdf_grid_masked["total_population"] * 100
    )

    print(f"Total pop (original): {np.nansum(total_pop_array):,.2f}")
    print(f"Total pop (aggregated): {gdf_grid_masked['total_population'].sum():,.2f}")
    diff_pct = (
        (np.nansum(total_pop_array) - gdf_grid_masked["total_population"].sum())
        / np.nansum(total_pop_array)
        * 100
    )
    print(f"Diff %: {diff_pct:,.2f}")

    return gdf_grid_masked


# ---------------------------------------------------------------------------
# Compute population counts per flood-depth bin
# ---------------------------------------------------------------------------
def compute_cdf_and_bins(gdf, bins, depth_col="flood_depth", pop_col="population"):
    """Return population counts per flood-depth bin."""
    flood = gdf[depth_col].values
    pop = gdf[pop_col].values
    valid = ~np.isnan(flood) & (pop > 0)
    flood, pop = flood[valid], pop[valid]
    return pd.Series(pop).groupby(pd.cut(flood, bins)).sum()


# ---------------------------------------------------------------------------
# Compute attributable fraction per flood-depth mask
# ---------------------------------------------------------------------------
def compute_attr_per_flood_depth_mask(diff, baseline, mask):
    """Return the attributable percentage for a given depth mask."""
    return np.nansum(diff[mask]) / np.nansum(baseline[mask]) * 100


# ---------------------------------------------------------------------------
# Save a raster to GeoTIFF
# ---------------------------------------------------------------------------
def save_raster(array, out_path, transform, crs):
    """Write a 2-D array as a single-band GeoTIFF."""
    print(f"▶ Writing raster: {out_path}")
    profile = {
        "driver": "GTiff",
        "dtype": rasterio.float32,
        "count": 1,
        "height": array.shape[0],
        "width": array.shape[1],
        "crs": crs,
        "transform": transform,
        "compress": "deflate",
    }
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(array.astype("float32"), 1)

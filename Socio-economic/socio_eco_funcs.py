import os
import json
import numpy as np
import pandas as pd
from os.path import join
import rasterio
from rasterio import features
import geopandas as gpd
import warnings
warnings.filterwarnings('ignore')
from rasterio.mask import mask
from shapely.geometry import box
from tqdm import tqdm


# --- Function to redistribute population over land pixels on flood grid ---
def reproject_and_redistribute_population_over_land(pop_path, land_gdf, flood_crs, flood_transform, flood_shape, province_geom=None, region=None, districts=None, year=None, out_raster_path=None):    
    print(f"▶ Loading {year} population raster...")
    with rasterio.open(pop_path) as src:
        pop = src.read(1, masked=True)
        pop_affine = src.transform
        pop_crs = src.crs

        if province_geom is not None:
            # Clip to province
            province_geom = province_geom.to_crs(src.crs)
            province_geom = [json.loads(province_geom.to_json())["features"][0]["geometry"]]

            pop_sofala, transform_sofala = mask(src, province_geom, crop=True, nodata=src.nodata)

        if districts is not None:
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
        "population": pop_vals.astype(int),
        "flood_depth": flood_vals,
        "x": xs,
        "y": ys
    })

    return df

# --- Function to aggregate population and flood depth to coarser grid ---
def aggregate_pop(total_pop_array, flood_raster, transform, crs, region=None, background=None, factor=100):
    print("Aggregating population raster to coarser grid polygons...")
    # Pixel size from transform
    pixel_width = transform.a
    pixel_height = -transform.e

    # Block (coarse) size in meters
    cell_width = factor * pixel_width
    cell_height = factor * pixel_height

    # Aggregate raster by factor
    def block_sum(arr, factor):
        nrows, ncols = arr.shape
        nrows_crop = nrows - nrows % factor
        ncols_crop = ncols - ncols % factor
        arr_cropped = arr[:nrows_crop, :ncols_crop]
        agg = arr_cropped.reshape(nrows_crop//factor, factor,
                                ncols_crop//factor, factor).sum(axis=(1,3))
        return agg

    def block_mean(arr, factor):
        nrows, ncols = arr.shape
        nrows_crop = nrows - nrows % factor
        ncols_crop = ncols - ncols % factor
        arr_cropped = arr[:nrows_crop, :ncols_crop]
        mean = arr_cropped.reshape(nrows_crop//factor, factor,
                                   ncols_crop//factor, factor).mean(axis=(1,3))
        return mean

    total_agg = block_sum(total_pop_array, factor)
    exposed_agg = block_sum(np.where(flood_raster > 0, total_pop_array, 0), factor)

    # Aggregate flood depth
    avg_flood_depth = block_mean(flood_raster, factor)

    # Coarse grid coordinates based on raster bounds
    nrows_coarse, ncols_coarse = total_agg.shape
    x0, y0 = transform * (0, 0)  # top-left pixel center
    x_coords = x0 + np.arange(ncols_coarse) * cell_width
    y_coords = y0 + np.arange(nrows_coarse) * -cell_height  # negative because raster Y decreases downward

    grid_cells = []
    for j, y in enumerate(y_coords):
        for i, x in enumerate(x_coords):
            grid_cells.append(box(x, y - cell_height, x + cell_width, y))

    # Create GeoDataFrame
    gdf_grid = gpd.GeoDataFrame(
        {
            "total_population": total_agg.flatten(),
            "exposed_population": exposed_agg.flatten(),
            "avg_flood_depth": avg_flood_depth.flatten()
        },
        geometry=grid_cells,
        crs=crs
    )

    # Clip to region/background
    region_bg = gpd.overlay(region, background, how="intersection") if background is not None else region.copy()
    gdf_grid_masked = gpd.overlay(gdf_grid, region_bg, how="intersection")

    # Relative population
    gdf_grid_masked["relative_population"] = (
        gdf_grid_masked["exposed_population"] / gdf_grid_masked["total_population"] * 100
    )
    return gdf_grid_masked


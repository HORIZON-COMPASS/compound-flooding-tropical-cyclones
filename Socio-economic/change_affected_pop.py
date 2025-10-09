#%%
import rasterio
import geopandas as gpd
import numpy as np
from shapely.geometry import box
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
import platform
import os
import matplotlib.pyplot as plt
from rasterio.warp import reproject, Resampling
from rasterio.mask import mask
import json
from matplotlib.colors import PowerNorm
import matplotlib.colors as colors
# import matplotlib.cm as cm
from shapely.geometry import Polygon
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap
from rasterio.vrt import WarpedVRT
import matplotlib.colors as mcolors
import rioxarray as rxr  # Required for reading TIFF files
import pandas as pd
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import geopandas as gpd
from shapely.geometry import box
from rasterio.transform import Affine
from rasterio.features import rasterize
from shapely.geometry import Point
import yaml
import xarray as xr
import os
from os.path import join
from hydromt_sfincs import SfincsModel, utils
from hydromt import DataCatalog
import itertools

prefix = "p:/" if platform.system() == "Windows" else "/p/"

def lat_formatter(x, pos):
    direction = 'N' if x >= 0 else 'S'
    return f"{abs(x):.1f}°{direction}"

def lon_formatter(x, pos):
    direction = 'E' if x >= 0 else 'W'
    return f"{abs(x):.1f}°{direction}"


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


#%%
# Input files
shapefile_fp = "p:/11210471-001-compass/03_Runs/sofala/Idai/sfincs/event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0/gis/region.geojson"   # replace with your region shapefile
background = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_region_background.geojson")

# population in provided inthousand persons per grid cell
population_raster_path_2020 = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_2020_30.tif")  
population_raster_path_1990 = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_1990_30.tif")  

# flood raster
F_flooding = sfincs_dir_F / "floodmap.tif"
CF_flooding = sfincs_dir_CF / "floodmap.tif"

# Flood model subgrid
sfincs_subgrid = BASE_RUN_PATH / "sfincs" / "subgrid" / "dep_subgrid.tif"



#%% Read flood data and backgorund polygons
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



#%% Read and regrid population data
def reproject_and_clip_population(pop_path, region_geom, flood_crs, flood_transform, flood_shape, year):
    print("Reproject population raster directly into flood grid CRS & resolution, clipped to the region of interest")
    with rasterio.open(pop_path) as src:
        # Build a WarpedVRT to act like the raster is already in flood CRS
        with WarpedVRT(src, crs=flood_crs, resampling=Resampling.sum) as vrt:
            # Clip to region (already in flood CRS)
            pop_clipped, transform_clipped = mask(vrt, region_geom, crop=True, nodata=src.nodata)

        # Prepare output array at flood grid resolution
        pop_regridded = np.zeros(flood_shape, dtype=pop_clipped.dtype)

        # Reproject clipped population → flood grid directly
        reproject(
            source=pop_clipped,
            destination=pop_regridded,
            src_transform=transform_clipped,
            src_crs=vrt.crs,
            dst_transform=flood_transform,
            dst_crs=flood_crs,
            resampling=Resampling.sum
        )

        print(f"[{year}] No-data value:", src.nodata)
        return pop_regridded
    

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


# --- Process population directly into flood grid ---
pop_arrays = {}
for year, path in [(1990, population_raster_path_1990),
                   (2020, population_raster_path_2020)]:
    pop_arrays[year] = reproject_and_clip_population(
        path, region_geom, flood_grid_crs, flood_grid_transform, flood_grid_shape, year
    )

# --- Compute exposed population GeoDataFrames ---
df_pop_2020_flood_depth_F  = pop_raster_to_df(pop_arrays[2020], hmax_F, flood_grid_transform)
df_pop_2020_flood_depth_CF = pop_raster_to_df(pop_arrays[2020], hmax_CF, flood_grid_transform)
df_pop_1990_flood_depth_F = pop_raster_to_df(pop_arrays[1990], hmax_F, flood_grid_transform)

# simple rasters for fast plotting
ra_exposed_pop_2020_F = np.where(hmax_F > 0, pop_arrays[2020], 0)
ra_exposed_pop_2020_CF = np.where(hmax_CF > 0, pop_arrays[2020], 0)
ra_exposed_pop_1990_F = np.where(hmax_F > 0, pop_arrays[1990], 0)

# Compute exposed pop on fine grid
df_pop_2020_exposed_F  = df_pop_2020_flood_depth_F[df_pop_2020_flood_depth_F['flood_depth'] > 0]
df_pop_2020_exposed_CF = df_pop_2020_flood_depth_CF[df_pop_2020_flood_depth_CF['flood_depth'] > 0]
df_pop_1990_exposed_F  = df_pop_1990_flood_depth_F[df_pop_1990_flood_depth_F['flood_depth'] > 0]

# --- Aggregate population to coarser grid ---
gdf_pop_2020_exposed_F_coarse = aggregate_pop(pop_arrays[2020], hmax_F, flood_grid_transform, flood_grid_crs, region_utm, background_utm)
gdf_pop_2020_exposed_CF_coarse = aggregate_pop(pop_arrays[2020], hmax_CF, flood_grid_transform, flood_grid_crs, region_utm, background_utm)
gdf_pop_1990_exposed_F_coarse = aggregate_pop(pop_arrays[1990], hmax_F, flood_grid_transform, flood_grid_crs, region_utm, background_utm)

# %%
print("2020 Factual exposed population stats:")
print("Max people exposed per cell:", np.nanmax(df_pop_2020_exposed_F["population"].astype(int)))
print("Total exposed people:", np.nansum(df_pop_2020_exposed_F["population"].astype(int)))
print("Count of exposed cells:", np.sum(~np.isnan(df_pop_2020_exposed_F['population'])))

print("2020 Counterfactual exposed population stats:")
print("Max people exposed per cell:", np.nanmax(df_pop_2020_exposed_CF["population"].astype(int)))
print("Total exposed people:", np.nansum(df_pop_2020_exposed_CF["population"].astype(int)))
print("Count of exposed cells:", np.sum(~np.isnan(df_pop_2020_exposed_CF['population'])))

print("1990 Factual exposed population stats:")
print("Max people exposed per cell:", np.nanmax(df_pop_1990_exposed_F["population"].astype(int)))
print("Total exposed people:", np.nansum(df_pop_1990_exposed_F["population"].astype(int)))
print("Count of exposed cells:", np.sum(~np.isnan(df_pop_1990_exposed_F['population'])))

print("One-line attribution numbers:")
print(f"Exposed population in 2020 Factual: {int(np.nansum(df_pop_2020_exposed_F['population'].astype(int))):,}")
print(f"Exposed population attributable to climate change: {int(np.nansum(df_pop_2020_exposed_F['population'].astype(int)) - np.nansum(df_pop_2020_exposed_CF['population'].astype(int))):,} {100 * (np.nansum(df_pop_2020_exposed_F['population'].astype(int)) - np.nansum(df_pop_2020_exposed_CF['population'].astype(int))) / np.nansum(df_pop_2020_exposed_F['population'].astype(int)):.2f}%")
print(f"Exposed population attributable to population change (2020-1990): {int(np.nansum(df_pop_2020_exposed_F['population'].astype(int)) - np.nansum(df_pop_1990_exposed_F['population'].astype(int))):,} {100 * (np.nansum(df_pop_2020_exposed_F['population'].astype(int)) - np.nansum(df_pop_1990_exposed_F['population'].astype(int))) / np.nansum(df_pop_2020_exposed_F['population'].astype(int)):.2f}%")
print(f"Population growth from 1990 to 2020 in Sofala region: {int(np.nansum(pop_arrays[2020]) - np.nansum(pop_arrays[1990])):,} {100 * (np.nansum(pop_arrays[2020]) - np.nansum(pop_arrays[1990])) / np.nansum(pop_arrays[2020]):.2f}%")


#%% ============================================================================================ # 
# Plot distribution of flood depth per exposed population
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

# --- Plot continuous CDF ---
plt.figure(figsize=(8,5))
plt.plot(flood_F_2020_sorted, cdf_F, label="Factual")
plt.plot(flood_CF_2020_sorted, cdf_CF_clim, label="Counterfactual Climate")
plt.plot(flood_F_1990_sorted, cdf_CF_pop, label="Counterfactual Population")
plt.xlim(0, 3.5)
plt.xlabel("Flood depth (m)")
plt.ylabel("Cumulative fraction of exposed population")
plt.title("Continuous CDF of exposed population by flood depth")
plt.legend()
plt.grid(True)
plt.show()

# --- Plot binned bar chart ---
bin_centers = bins[:-1] + np.diff(bins) / 2
plt.figure(figsize=(8,5))
plt.bar(bin_centers - 0.05, pop_2020_by_depth_F.values, width=0.05, label="Factual", alpha=0.7)
plt.bar(bin_centers, pop_2020_by_depth_CF.values, width=0.05, label="Counterfactual Climate", alpha=0.7)
plt.bar(bin_centers + 0.05, pop_1990_by_depth_F.values, width=0.05, label="Counterfactual Population", alpha=0.7)

plt.xlabel("Flood depth (m)")
plt.ylabel("Exposed population")
plt.title("Population exposed by flood depth bins")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.show()



#%% ============================================================================================ # 
# Plot the factual flood and exposed population
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


#%% ============================================================================================ # 
# Plot the difference between 1990 and 2020 of the HE dataset, and between WP 2020 and HE 2020
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
# ===================== Differences in AGGREGATED exposed population ======================== #
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
print("Plotting spatially aggregated exposed population damage for F, CF climate and diff")

gdf_pop_2020_exposed_CF_coarse['population_diff'] = (gdf_pop_2020_exposed_F_coarse['exposed_population'] - gdf_pop_2020_exposed_CF_coarse['exposed_population'])

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 5), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

# Create colormap normalization
norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2020_exposed_F_coarse["exposed_population"].max())
norm_diff = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2020_exposed_CF_coarse["population_diff"].max())
red_half = LinearSegmentedColormap.from_list("bwr_red", plt.cm.bwr(np.linspace(0.5, 1, 256)))

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


df_pop_2020_exposed_F['fatalities'] = df_pop_2020_exposed_F['population'] * fatality_curve_Boyd2005(df_pop_2020_exposed_F['flood_depth'])

# %%
CF_rain_flooding = BASE_RUN_PATH / "sfincs" / "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" / "floodmap.tif"
CF_SLR_wind_flooding = BASE_RUN_PATH / "sfincs" / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10" / "floodmap.tif"
hmax_CF_rain_da = rxr.open_rasterio(CF_rain_flooding).squeeze("band", drop=True)  # if single-band
hmax_CF_SLR_wind_da = rxr.open_rasterio(CF_SLR_wind_flooding).squeeze("band", drop=True)  # if single-band

# river flooding vs. coastal flooding
def calculate_flood_type_area(df, hmax_F_da, hmax_CF_da, hmax_CF_rain_da, hmax_CF_SLR_wind_da):
    # Masks for valid values
    mask_F_valid = ~np.isnan(hmax_F_da.values)
    mask_CF_total_valid = ~np.isnan(hmax_CF_da.values)
    mask_CF_rain_valid = ~np.isnan(hmax_CF_rain_da.values)
    mask_CF_SLR_wind_valid = ~np.isnan(hmax_CF_SLR_wind_da.values)

    # Replace .where() with np.where()
    hmax_F_filled = np.where(mask_F_valid | ~mask_CF_total_valid, hmax_F_da.values, 0)
    hmax_CF_total_filled = np.where(mask_CF_total_valid | ~mask_F_valid, hmax_CF_da.values, 0)
    hmax_diff_total = hmax_F_filled - hmax_CF_total_filled

    # Fluvial (rain) difference
    hmax_CF_fluvial_filled = np.where(mask_CF_rain_valid | ~mask_F_valid, hmax_CF_rain_da.values, 0)
    hmax_diff_fluvial = hmax_F_filled - hmax_CF_fluvial_filled

    # Coastal (SLR + wind) difference
    hmax_CF_coastal_filled = np.where(mask_CF_SLR_wind_valid | ~mask_F_valid, hmax_CF_SLR_wind_da.values, 0)
    hmax_diff_coastal = hmax_F_filled - hmax_CF_coastal_filled

    # Masks for each type
    fluvial_mask = (hmax_diff_fluvial > 0)
    coastal_mask = (hmax_diff_coastal > 0)

    # Interpolate flood type for df based on (x,y)
    from scipy.interpolate import RegularGridInterpolator

    # Create interpolators
    interp_fluvial = RegularGridInterpolator((hmax_F_da.y, hmax_F_da.x), fluvial_mask, bounds_error=False, fill_value=False)
    interp_coastal = RegularGridInterpolator((hmax_F_da.y, hmax_F_da.x), coastal_mask, bounds_error=False, fill_value=False)

    # Apply to df
    df = df.copy()
    df['is_fluvial'] = interp_fluvial(df[['y', 'x']])
    df['is_coastal'] = interp_coastal(df[['y', 'x']])

    # Assign flood type
    df['flood_type'] = 'none'
    df.loc[df['is_fluvial'] == 1, 'flood_type'] = 'fluvial'
    df.loc[df['is_coastal'] == 1, 'flood_type'] = 'coastal'

    return df

df_pop_2020_exposed_F_type = calculate_flood_type_area(df_pop_2020_exposed_F, hmax_F_da, hmax_CF_da, hmax_CF_rain_da, hmax_CF_SLR_wind_da)
# %%
def mortality_rate(depth, a, b):
    return 1 / (1 + np.exp(-(a + b * depth)))

# Jonkman parameters
params = {
    "fluvial": {"a": -6.4, "b": 2.0},
    "coastal": {"a": -5.4, "b": 2.0},
}

def compute_fatalities(row):
    ft = row['flood_type']
    if ft == 'none' or row['flood_depth'] <= 0:
        return 0
    a, b = params[ft]['a'], params[ft]['b']
    pf = mortality_rate(row['flood_depth'], a, b)
    return pf * row['population']

df_pop_2020_exposed_F_type['expected_fatalities'] = df_pop_2020_exposed_F_type.apply(compute_fatalities, axis=1)

# %%
import gc
from rasterio.features import shapes
from shapely.geometry import shape

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

#%%
# Compute the maximum water level (hmax) and mask out permanent water
def compute_hmax_masked(models, gwso_region):
    # we set a threshold to mask minimum flood depth
    hmin = 0.05

    for model in models:
        # select the highest-resolution elevation dataset
        print(f"Processing model: {model['model_name']}")
        depfile = join(model["model_path"], "subgrid", "dep_subgrid.tif")
        da_dep = model["sfincs_model"].data_catalog.get_rasterdataset(depfile)

        # compute the maximum over all time steps
        # First timestep leads to incorrect diiference values for permanent water cells that are incorrectly unmasked; requires filtering.
        da_zsmax = model["sfincs_results"]["zsmax"].isel(timemax=slice(1, None)).max(dim="timemax")
        da_zs = model["sfincs_results"]["zs"]
        da_zb = model["sfincs_results"]["zb"]
     
        # downscale the floodmap for max water depth
        da_hmax = utils.downscale_floodmap(
            zsmax=da_zsmax,
            dep=da_dep,
            hmin=hmin,
            reproj_method = "bilinear"
            )
        
        # downscale the floodmap for hourly water depth
        # da_h = utils.downscale_floodmap(
        #     zsmax=da_zs,
        #     dep=da_dep,
        #     hmin=hmin,
        #     reproj_method = "bilinear"
        #     )
        da_h = da_zs - da_zb
        da_h = da_h.where(da_h > hmin, np.nan)
    
        # GSWO dataset to mask permanent water by first geprojecting it to the subgrid of hmax
        gswo_mask = gwso_region.raster.reproject_like(da_hmax, method="max")
        # permanent water where water occurence > 5%
        da_hmax_masked = da_hmax.where(gswo_mask <= 5)
        da_h_masked = da_h.where(gswo_mask <= 5)

        # --- Compute velocity using Manning's equation ---
        manning_file = join(model["model_path"], "subgrid", "manning_subgrid.tif")
        manning_ra = rxr.open_rasterio(manning_file).squeeze("band", drop=True) 
        manning_subgrid = manning_ra.rio.reproject_match(da_h_masked)

        # Extract grid spacing from coordinates (assumes regular grid)
        dx = float(da_h_masked['x'][0,1] - da_h_masked['x'][0,0])
        dy = float(da_h_masked['y'][1,0] - da_h_masked['y'][0,0])

        # Compute slopes
        dH_dx = da_h_masked.differentiate('m') / dx  # slope in x
        dH_dy = da_h_masked.differentiate('n') / dy  # slope in y
        S = np.sqrt(dH_dx**2 + dH_dy**2)  # slope magnitude

        # Compute velocity
        da_v = (1 / manning_subgrid) * da_h_masked**(2/3) * S**0.5
        da_v = da_v.fillna(0).rename('velocity')
        da_v.attrs['units'] = 'm/s'

        # --- Save results ---
        model["sfincs_results"]['hmax'] = da_hmax
        model["sfincs_results"]['hmax_masked'] = da_hmax_masked
        model["sfincs_results"]['h_masked'] = da_h_masked
        model["sfincs_results"]['velocity'] = da_v

        del da_hmax, da_zsmax, da_zs, da_dep, gswo_mask  # Clean up to free memory
        gc.collect()

    return models

# %%
print("=== Step 1: Extract SFINCS results ===")
da_zs = model["sfincs_results"]["zs"]      # (time, y, x)
da_zb = model["sfincs_results"]["zb"]      # (y, x)
print("  - Extracted water surface and bed elevation.")

print("=== Step 2: Compute hourly depth ===")
da_h = da_zs - da_zb                        # (time, y, x)
da_h_valid = da_h.where(da_h > hmin)
valid_mask = ~da_h_valid.isnull().all(dim="time")
print(f"  - Applied hmin={hmin} mask. Valid pixels: {valid_mask.sum().values}/{valid_mask.size}")

print("=== Step 3: Compute per-pixel max depth and time index ===")
imax = da_h_valid.argmax(dim="time", skipna=True)
imax = imax.where(valid_mask)
da_hmax = da_h_valid.max(dim="time", skipna=True)
print("  - Computed max depth per pixel and index of max time.")

print("=== Step 4: Compute slope from water surface ===")
dx = float(da_zs['x'][1] - da_zs['x'][0])
dy = float(da_zs['y'][1] - da_zs['y'][0])
dZs_dx = da_zs.differentiate("x") / dx
dZs_dy = da_zs.differentiate("y") / dy
S = np.sqrt(dZs_dx**2 + dZs_dy**2)
print("  - Computed slope magnitude over time.")

print("=== Step 5: Extract slope at time of max depth ===")
# Bring slope to memory (or compute chunkwise if too big)
# Define a small function that selects per-pixel slope at the given time index
def select_slope_at_hmax(S_block, imax_block):
    out = np.take_along_axis(S_block, imax_block[None, ...], axis=0)
    return out.squeeze(0)

# Apply blockwise (Dask parallelized)
S_at_hmax = xr.map_blocks(
    select_slope_at_hmax,
    S,
    args=(imax),
    template=xr.zeros_like(imax, dtype=float)
)
S_at_hmax.name = "slope_at_hmax"
print("  - Extracted slope at max depth per pixel (chunked + parallelized).")


print("=== Step 6: Mask permanent water ===")
gswo_mask = gwso.raster.reproject_like(da_hmax, method="max")
mask_perm_water = gswo_mask > 5
S_masked = S_at_hmax.where(~mask_perm_water)
hmax_masked = da_hmax.where(~mask_perm_water)
print(f"  - Applied permanent water mask and clipped hmax to hmin={hmin}.")

print("=== Step 7: Reproject to Manning grid ===")
S_sub = S_masked.rio.reproject_match(manning_subgrid)
hmax_sub = hmax_masked.rio.reproject_match(manning_subgrid)
print("  - Reprojected slope and max depth to Manning grid.")

print("=== Step 8: Compute velocity using Manning's equation ===")
eps = 1e-9
da_v_hmax = (1.0 / (manning_subgrid + eps)) * (hmax_sub ** (2.0 / 3.0)) * np.sqrt(S_sub + eps)
da_v_hmax = da_v_hmax.fillna(0).rename("velocity")
da_v_hmax.attrs["units"] = "m/s"
print("  - Computed velocity at max water depth.")

print("=== Done: Velocity at max water depth computed ===")
da_v_hmax


#%%

# Function to compute hourly velocity using Manning's equation
def compute_hourly_velocity(ds, n_manning):
    # 1. Compute water depth
    h = ds['zs'] - ds['zb']  # (n, m, time)
    h = h.where(h > 0, 0.0)  # ensure non-negative depths

    # 2. Compute slopes in x and y directions
    # Get grid spacing from coordinates (assuming uniform grid)
    dx = float(ds['x'][0,1] - ds['x'][0,0])  # meters per cell
    dy = float(ds['y'][1,0] - ds['y'][0,0])

    # compute gradients
    dH_dx = h.differentiate('m') / dx  # approximate slope in x-direction
    dH_dy = h.differentiate('n') / dy  # approximate slope in y-direction

    # magnitude of slope
    S = np.sqrt(dH_dx**2 + dH_dy**2)
    
    # 3. Compute velocity using Manning's equation
    # v = (1/n) * h^(2/3) * S^(1/2)
    v = (1 / n_manning) * h**(2/3) * S**0.5
    v = v.fillna(0)  # fill NaNs with 0 for dry cells

    # 4. Attach coordinates and return as DataArray
    v = v.rename('velocity')
    v.attrs['units'] = 'm/s'
    
    return v


#%%
# Load snakemake config file to construct the model paths
config_path  = '../Workflows/01_config_snakemake/config_general_MZB.yml'
cfg = load_config(config_path)

#%%
# Load the SFINCS models in one dictonary
models = load_sfincs_models(cfg)
print(models)
# %%
# Calculate hmax and mask out permanent water
gwso = gwso_sfincs_region(models[0])

#%%
models, gdf_valid = compute_hmax_masked(models, gwso, model_region)
models = compute_hmax_diff(models)

# %%

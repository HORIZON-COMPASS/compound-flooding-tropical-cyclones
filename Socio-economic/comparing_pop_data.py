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
population_raster_path_WP = os.path.join(prefix,"11210471-001-compass","01_Data","population_data","moz_ppp_2020_UNadj_constrained.tif")

# flood raster
F_flooding = sfincs_dir_F / "floodmap.tif"
CF_flooding = sfincs_dir_CF / "floodmap.tif"

#%%
# --- Load region ---
region = gpd.read_file(shapefile_fp)
with rasterio.open(population_raster_path_1990) as src1990:
    region = region.to_crs(src1990.crs)  # match CRS

region_geom = [json.loads(region.to_json())["features"][0]["geometry"]]

background = background.to_crs("EPSG:4326")  # Do once

# --- Clip rasters ---
with rasterio.open(population_raster_path_1990) as src_HE_1990:
    pop_HE_1990, transform_HE_1990 = mask(src_HE_1990, region_geom, crop=True)
    print("No-data value Historical Exposure 1990:", src_HE_1990.nodata)

with rasterio.open(population_raster_path_2020) as src_HE_2020:
    pop_HE_2020, transform_HE_2020 = mask(src_HE_2020, region_geom, crop=True)
    print("No-data value Historical Exposure 2020:", src_HE_2020.nodata)

with rasterio.open(population_raster_path_WP) as src_WP_2020:
    pop_WP_2020, transform_WP_2020 = mask(src_WP_2020, region_geom, crop=True)
    print("No-data value World Pop 2020:", src_WP_2020.nodata)
    pop_WP_2020[pop_WP_2020 == src_WP_2020.nodata] = 0


# Example: get extent from raster transform
def get_extent(transform, width, height):
    left = transform[2]
    right = left + width * transform[0]
    top = transform[5]
    bottom = top + height * transform[4]
    return [left, right, bottom, top]

extent_WP = get_extent(transform_WP_2020, pop_WP_2020.shape[2], pop_WP_2020.shape[1])
extent_HE2020 = get_extent(transform_HE_2020, pop_HE_2020.shape[2], pop_HE_2020.shape[1])
extent_HE1990 = get_extent(transform_HE_1990, pop_HE_1990.shape[2], pop_HE_1990.shape[1])

#%%
# Plot the raw pop data for the region
vmax_WP = np.percentile(pop_WP_2020[0], 99.9)  # 99th percentile
vmax_HE2020 = np.percentile(pop_HE_2020[0], 99.9)
vmax_HE1990 = np.percentile(pop_HE_1990[0], 99.9)

# Mask raster outside region (already cropped with rasterio.mask.mask)
pop_WP_masked = np.where(pop_WP_2020[0] == 0, np.nan, pop_WP_2020[0])
pop_HE_2020_masked = np.where(pop_HE_2020[0] == 0, np.nan, pop_HE_2020[0])
pop_HE_1990_masked = np.where(pop_HE_1990[0] == 0, np.nan, pop_HE_1990[0])

# Define a polygon to remove/mask out
mask_poly = Polygon([(34.9,-20.3), (36,-20.3), (36,-19.8), (34.9,-19.8)])

# Subtract the polygon from all geometries
bg_filtered = background.copy()
bg_filtered['geometry'] = bg_filtered.geometry.apply(lambda g: g.difference(mask_poly))


fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for ax in axes:
    background.plot(ax=ax, color="#E0E0E0", zorder=0)
    bg_filtered.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
    region.boundary.plot(ax=ax, color='black', linewidth=1, zorder=2)

# Plot World pop data
im = axes[0].imshow(pop_WP_masked, cmap="Blues", extent=extent_WP,  origin='upper', vmin=0, vmax=vmax_WP) 
cbar = plt.colorbar(im, ax=axes[0], shrink=0.8)
cbar.set_label("Population (people per cell)")
axes[0].set_title("World Pop 2020 \n[100 m grid]")

# Plot Historical Exposure data 2020
im = axes[1].imshow(pop_HE_2020_masked, cmap="Blues", extent=extent_HE2020,  origin='upper', vmin=0, vmax=vmax_HE2020) 
cbar = plt.colorbar(im, ax=axes[1], shrink=0.8)
cbar.set_label("Population (people per cell)")
axes[1].set_title("Hist. Exposure 2020 \n[~1 km grid]")

# Plot Historical Exposure data 1990
im = axes[2].imshow(pop_HE_1990_masked, cmap="Blues", extent=extent_HE1990,  origin='upper', vmin=0, vmax=vmax_HE1990)
cbar = plt.colorbar(im, ax=axes[2], shrink=0.8)
cbar.set_label("Population (people per cell)")
axes[2].set_title("Hist. Exposure 1990 \n[~1 km grid]")

for ax in axes:
    region.boundary.plot(ax=ax, color='black', linewidth=1)


# %% Regrid WP data to HE for comparison
pop_WP_2020_coarse = np.zeros(pop_HE_2020.shape, dtype=np.float32)

# Reproject WP → HE grid
reproject(
    source=pop_WP_2020,
    destination=pop_WP_2020_coarse,
    src_transform=transform_WP_2020,
    src_crs=src_WP_2020.crs,
    dst_transform=transform_HE_2020,
    dst_crs=src_HE_2020.crs,
    resampling=Resampling.sum   # sum keeps population counts correct
)


#%% Plot all pop data on the same grid
# cmap = cm.get_cmap("Blues").copy()
# cmap.set_bad(color="#E0E0E0")  # grey

# Mask raster outside region (already cropped with rasterio.mask.mask)
pop_WP_masked = np.where(pop_WP_2020_coarse[0] == 0, np.nan, pop_WP_2020_coarse[0])
pop_HE_2020_masked = np.where(pop_HE_2020[0] == 0, np.nan, pop_HE_2020[0])
pop_HE_1990_masked = np.where(pop_HE_1990[0] == 0, np.nan, pop_HE_1990[0])

vmax_WP = np.percentile(pop_WP_masked[~np.isnan(pop_WP_masked)], 99.9)

# Define a polygon to remove/mask out
mask_poly = Polygon([(34.9,-20.3), (36,-20.3), (36,-19.8), (34.9,-19.8)])

# Subtract the polygon from all geometries
bg_filtered = background.copy()
bg_filtered['geometry'] = bg_filtered.geometry.apply(lambda g: g.difference(mask_poly))


fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for ax in axes:
    background.plot(ax=ax, color="#E0E0E0", zorder=0)
    bg_filtered.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
    region.boundary.plot(ax=ax, color='black', linewidth=1, zorder=2)

# Plot WorldPop 2020
im = axes[0].imshow(pop_WP_masked, cmap="Blues", extent=extent_WP, origin='upper', norm=PowerNorm(gamma=0.5, vmin=0, vmax=vmax_WP))
cbar = plt.colorbar(im, ax=axes[0], shrink=0.8)
cbar.set_label("Population (people per cell)")
axes[0].set_title("World Pop 2020 \n[~1 km grid]")

# Plot Historical Exposure data 2020
im = axes[1].imshow(pop_HE_2020_masked, cmap="Blues", extent=extent_HE2020,  origin='upper', norm=PowerNorm(gamma=0.5, vmin=0, vmax=vmax_HE2020)) 
cbar = plt.colorbar(im, ax=axes[1], shrink=0.8)
cbar.set_label("Population (people per cell)")
axes[1].set_title("Hist. Exposure 2020 \n[~1 km grid]")

# Plot Historical Exposure data 1990
im = axes[2].imshow(pop_HE_1990_masked, cmap="Blues", extent=extent_HE1990,  origin='upper', norm=PowerNorm(gamma=0.5, vmin=0, vmax=vmax_HE1990))
cbar = plt.colorbar(im, ax=axes[2], shrink=0.8)
cbar.set_label("Population (people per cell)")
axes[2].set_title("Hist. Exposure 1990 \n[~1 km grid]")



#%% Plot the difference between 1990 and 2020 of the HE dataset, and between WP 2020 and HE 2020
# --- Make sure dimensions match ---
fig, axes = plt.subplots(1, 2, figsize=(12, 6))  # 1 row, 3 columns

# --- Make sure dimensions match ---
if pop_WP_2020_coarse.shape != pop_HE_2020.shape:
    raise ValueError("WP 2020 and HE 2020 clipped rasters have different shapes!")
if pop_HE_1990.shape != pop_HE_2020.shape:
    raise ValueError("1990 and 2020 clipped rasters have different shapes!")

# --- Compute difference ---
pop_diff_ds = pop_HE_2020[0] - pop_WP_2020_coarse[0]   # increase (+) or decrease (-)
pop_diff_yrs = pop_HE_2020[0] - pop_HE_1990[0]   # increase (+) or decrease (-)

# --- Plot ---
for ax in axes:
    region.boundary.plot(ax=ax, color='black', linewidth=1)
    background.plot(ax=ax, color='#E0E0E0', zorder=0)
    bg_filtered.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)

# Create colormap that starts with white for zero
cmap = plt.cm.Reds
cmap = colors.ListedColormap(['white'] + [cmap(i) for i in range(1, cmap.N)])

pop_diff_ds_masked = np.where(pop_diff_ds == 0, np.nan, pop_diff_ds)
pop_diff_yrs_masked = np.where(pop_diff_yrs == 0, np.nan, pop_diff_yrs)

im = axes[0].imshow(pop_diff_ds_masked, cmap=cmap, extent=extent_WP, norm=PowerNorm(gamma=0.5, vmin=0, vmax=np.percentile(pop_diff_ds, 99.9)))
cbar = plt.colorbar(im, ax=axes[0], shrink=0.8)
cbar.set_label("Population difference (HE - WP)")
axes[0].set_title("HE pop diff between World Pop and Hist. Exposure")

im = axes[1].imshow(pop_diff_yrs_masked, cmap=cmap, extent=extent_HE2020, norm=PowerNorm(gamma=0.5, vmin=0, vmax=np.percentile(pop_diff_yrs, 99.9)))
cbar = plt.colorbar(im, ax=axes[1], shrink=0.8)
cbar.set_label("Population change (2020 - 1990)")
axes[1].set_title("HE pop diff between 2020 and 1990")

#%%
# Some numbers
print(f"Total population HE 2020: {pop_HE_2020[0].sum():.0f} people")
print(f"Total population HE 1990: {pop_HE_1990[0].sum():.0f} people")
print(f"Total population WP 2020: {pop_WP_2020[0].sum():.0f} people")
print(f"Total population WP 2020 (coarsened): {pop_WP_2020_coarse[0].sum():.0f} people") # sanity check

print(f"Difference in total population between HE 2020 and 1990: {pop_diff_yrs.sum():.0f} people")
print(f"Difference in total population between HE 2020 and WP 2020: {pop_diff_ds.sum():.0f} people")


#%%
# ============================================================================================ #
# ============================= Differences in exposed population ============================ #
# ============================================================================================ #
WP_2020 = BASE_RUN_PATH / "fiat" / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" / "output" / "spatial_with_pop_and_flood_worldpop.fgb"
HE_2020 = BASE_RUN_PATH / "fiat" / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" / "output" / "spatial_with_pop_and_flood_hist_exposure.fgb"
HE_1990 = BASE_RUN_PATH / "fiat" / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" / "output" / "spatial_with_pop_and_flood_hist_exposure_1990.fgb"

# Read the geodataframes
gdf_WP_2020 = gpd.read_file(WP_2020)
gdf_HE_2020 = gpd.read_file(HE_2020)
gdf_HE_1990 = gpd.read_file(HE_1990)


#%%
# Function to process a population GeoDataFrame
print("Extracting coordinates from population centroids...")
def process_gdf(gdf_path, value_column='population', crs=ccrs.PlateCarree()):
    gdf = gpd.read_file(gdf_path)
    crs = gdf.crs
    
    # Compute centroids
    gdf['centroid'] = gdf.geometry.centroid
    gdf['x'] = gdf['centroid'].x
    gdf['y'] = gdf['centroid'].y
    
    # Keep only relevant columns
    gdf[value_column] = gdf[value_column]
    gdf = gdf[['object_id', 'x', 'y', value_column]]
    
    # Remove rows with missing coordinates
    gdf = gdf.dropna(subset=['x', 'y'])
    
    # Convert back to GeoDataFrame with points
    gdf_points = gpd.GeoDataFrame(
        gdf,
        geometry=gpd.points_from_xy(gdf['x'], gdf['y']),
        crs=crs
    )
    
    # Reproject to target CRS
    gdf_points = gdf_points.to_crs(crs)

    print(gdf_points.head())
    
    return gdf_points

# Apply to all three datasets
filt_WP_2020 = process_gdf(WP_2020)
filt_HE_2020 = process_gdf(HE_2020)
filt_HE_1990 = process_gdf(HE_1990)


#%%
def aggregate_damage_to_grid(gdf_pop, region, background, cell_size=0.025):
    # Clip region to background
    clipped_region = gpd.overlay(region, background, how="intersection")
    
    # Bounds
    minx, miny, maxx, maxy = clipped_region.total_bounds
    cols = np.arange(minx, maxx, cell_size)
    rows = np.arange(miny, maxy, cell_size)
    
    # Build grid cells
    grid_cells = [box(x, y, x + cell_size, y + cell_size) for x in cols for y in rows]
    gdf_grid = gpd.GeoDataFrame(geometry=grid_cells, crs=clipped_region.crs)
    
    # Clip grid to region
    gdf_grid_masked = gpd.overlay(gdf_grid, clipped_region, how='intersection')
    
    # Spatial join: assign grid cell index to each point
    joined = gpd.sjoin(gdf_pop, gdf_grid_masked, how="left", predicate="within")
    
    # Aggregate damage/population per grid cell
    if "population" in gdf_pop.columns:
        agg_tot = joined.groupby(joined.index_right)["population"].sum()
        gdf_grid_masked["population"] = agg_tot
        gdf_grid_masked["population"] = gdf_grid_masked["population"].fillna(0)
    else:
        print("No population column found!")
    
    return gdf_grid_masked

gdf_grid_WP_2020 = aggregate_damage_to_grid(filt_WP_2020, region, background)
gdf_grid_HE_2020 = aggregate_damage_to_grid(filt_HE_2020, region, background)
gdf_grid_HE_1990 = aggregate_damage_to_grid(filt_HE_1990, region, background)

# %%
#%%
print("Plotting spatially aggregated total damage for WP and HE 2020")
# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5), dpi=300, sharey=True, constrained_layout=True,
                         subplot_kw={"projection": ccrs.PlateCarree()})

# Create colormap normalization
norm_WP = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_grid_WP_2020['population'].max())
norm_HE = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_grid_HE_2020["population"].max())
# red_half = LinearSegmentedColormap.from_list("bwr_red", plt.cm.bwr(np.linspace(0.5, 1, 256)))

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
gdf_grid_WP_2020[gdf_grid_WP_2020['population'] == 0].plot(ax=axes[0], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_grid_WP_2020[gdf_grid_WP_2020['population'] > 0].plot(column='population', cmap='Blues', norm=norm_WP, edgecolor='grey', 
                                                                linewidth=0.2, ax=axes[0], legend=False, zorder=2, rasterized=True)

gdf_grid_HE_2020[gdf_grid_HE_2020['population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_grid_HE_2020[gdf_grid_HE_2020['population'] > 0].plot(column='population', cmap='Blues', norm=norm_HE, edgecolor='grey', 
                                                                 linewidth=0.2, ax=axes[1], legend=False, zorder=2, rasterized=True)
subplot_labels = ['(a)', '(b)']

for i, ax in enumerate(axes):
    # Add model region
    region.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())

    # # Add background and set extent (based on actual lat/lon coordinates)
    background.plot(ax=ax, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)
    minx, miny, maxx, maxy = region.bounds.minx.item(), region.bounds.miny.item(), region.bounds.maxx.item(), region.bounds.maxy.item()
    ax.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

    # Add gridlines and format tick labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.FixedLocator(np.arange(minx, maxx + 0.1, 0.2))
    gl.ylocator = mticker.FixedLocator(np.arange(miny, maxy + 0.1, 0.2))
    gl.xformatter = mticker.FuncFormatter(lon_formatter)
    gl.yformatter = mticker.FuncFormatter(lat_formatter)
    gl.right_labels = False
    gl.top_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}
    if i != 0: 
        gl.left_labels = False
    
    ax.text(0, 1.02, subplot_labels[i], transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='bottom', ha='left')

 # Colorbars
sm1 = ScalarMappable(norm=norm_WP, cmap="Blues")
sm1.set_array([])
cbar = fig.colorbar(sm1, ax=axes[0], orientation="vertical", shrink=0.6, pad=0.01)
cbar.set_label("Aggregated exposed population", fontsize = 10)
cbar.ax.tick_params(labelsize=9)

sm1 = ScalarMappable(norm=norm_HE, cmap="Blues")
sm1.set_array([])
cbar = fig.colorbar(sm1, ax=axes[1], orientation="vertical", shrink=0.6, pad=0.01)
cbar.set_label("Aggregated exposed population", fontsize = 10)
cbar.ax.tick_params(labelsize=9)

axes[0].set_title("World Pop 2020", fontsize=11)
axes[1].set_title("Hist. Exposure", fontsize=11)

# fig.savefig("../figures/fS12.png", bbox_inches='tight', dpi=300)
# fig.savefig("../figures/fS12.pdf", bbox_inches='tight', dpi=300)
# %%
print("Plotting spatially aggregated total damage for F, CF and diff")

gdf_grid_HE_1990['population_diff'] = (gdf_grid_HE_2020['population'] - gdf_grid_HE_1990['population'])

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 5), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.PlateCarree()})

# Create colormap normalization
norm_HE_2020 = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_grid_HE_2020["population"].max())
norm_HE_1990 = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_grid_HE_1990["population"].max())
norm_diff = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_grid_HE_1990["population_diff"].max())
red_half = LinearSegmentedColormap.from_list("bwr_red", plt.cm.bwr(np.linspace(0.5, 1, 256)))

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
gdf_grid_HE_2020[gdf_grid_HE_2020['population'] == 0].plot(ax=axes[0], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_grid_HE_2020[gdf_grid_HE_2020['population'] > 0].plot(column='population', cmap='Blues', norm=norm_HE_2020, edgecolor='grey', 
                                                                linewidth=0.2, ax=axes[0], legend=False, zorder=2, rasterized=True)

gdf_grid_HE_1990[gdf_grid_HE_1990['population'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_grid_HE_1990[gdf_grid_HE_1990['population'] > 0].plot(column='population', cmap='Blues', norm=norm_HE_1990, edgecolor='grey', 
                                                                 linewidth=0.2, ax=axes[1], legend=False, zorder=2, rasterized=True)

gdf_grid_HE_1990[gdf_grid_HE_1990['population_diff'] <= 0].plot(ax=axes[2], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_grid_HE_1990[gdf_grid_HE_1990['population_diff'] > 0].plot(column='population_diff', cmap=red_half, norm=norm_diff, edgecolor='grey', 
                                                                 linewidth=0.2, ax=axes[2], legend=False, zorder=2, missing_kwds={"color": "white", "edgecolor": "none"}, rasterized=True)
subplot_labels = ['(a)', '(b)', '(c)']

for i, ax in enumerate(axes):
    # Add model region
    region.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())

    # # Add background and set extent (based on actual lat/lon coordinates)
    background.plot(ax=ax, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)
    minx, miny, maxx, maxy = region.bounds.minx.item(), region.bounds.miny.item(), region.bounds.maxx.item(), region.bounds.maxy.item()
    ax.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

    # Add gridlines and format tick labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.FixedLocator(np.arange(minx, maxx + 0.1, 0.2))
    gl.ylocator = mticker.FixedLocator(np.arange(miny, maxy + 0.1, 0.2))
    gl.xformatter = mticker.FuncFormatter(lon_formatter)
    gl.yformatter = mticker.FuncFormatter(lat_formatter)
    gl.right_labels = False
    gl.top_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}
    if i != 0: 
        gl.left_labels = False
    
    ax.text(0, 1.02, subplot_labels[i], transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='bottom', ha='left')

# Colorbars
sm1 = ScalarMappable(norm=norm_HE_2020, cmap="Blues")
sm1.set_array([])
cbar = fig.colorbar(sm1, ax=axes[0], orientation="vertical", shrink=0.4, pad=0.01)
cbar.set_label("Aggregated exposed population [# people]", fontsize = 9)
cbar.ax.tick_params(labelsize=9)

sm1 = ScalarMappable(norm=norm_HE_1990, cmap="Blues")
sm1.set_array([])
cbar = fig.colorbar(sm1, ax=axes[1], orientation="vertical", shrink=0.4, pad=0.01)
cbar.set_label("Aggregated exposed population [# people]", fontsize = 9)
cbar.ax.tick_params(labelsize=9)

sm2 = ScalarMappable(cmap=red_half, norm=norm_diff)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, ax=axes[2], orientation="vertical", shrink=0.4, pad=0.01)
cbar2.set_label("Attributable affected population [# people]", fontsize = 9)
cbar2.ax.tick_params(labelsize=8)

axes[0].set_title("Factual (2020)", fontsize=10)
axes[1].set_title("Counterfactual (1990)", fontsize=10)
axes[2].set_title("Factual - Counterfactual", fontsize=10)

# fig.savefig("../figures/fS12.png", bbox_inches='tight', dpi=300)
# fig.savefig("../figures/fS12.pdf", bbox_inches='tight', dpi=300)
# %%
# Statistics
total_2020_HE = gdf_grid_HE_2020["population"].sum()
total_1990_HE = gdf_grid_HE_1990["population"].sum()
total_2020_WP = gdf_grid_WP_2020["population"].sum()
total_diff = gdf_grid_HE_1990['population_diff'].sum()
perc_diff = (gdf_grid_HE_2020['population'].sum() - gdf_grid_HE_1990['population'].sum()) / gdf_grid_HE_2020['population'].sum() * 100

print(f"Total exposed population 2020 (WP): {total_2020_WP:.0f} people")
print(f"Total exposed population 2020 (HE): {total_2020_HE:.0f} people")
print(f"Total exposed population 1990 (HE): {total_1990_HE:.0f} people")
print(f"Attributable exposed population to population change: {total_diff:.0f} people")
print(f"Attributable exposed population to population change: {perc_diff:.0f} %")
# %%
# For "mock up" fig
gdf_grid_HE_2020.to_file("data/gdf_grid_HE_2020.gpkg", driver="GPKG")
gdf_grid_HE_1990.to_file("data/gdf_grid_HE_1990.gpkg", driver="GPKG")

# %%

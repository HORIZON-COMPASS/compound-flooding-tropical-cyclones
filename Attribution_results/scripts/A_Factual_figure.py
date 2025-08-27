#%%
print("Loading packages")
import os
from os.path import join
from pathlib import Path
import xarray as xr
import numpy as np
import geopandas as gpd

from hydromt_sfincs import SfincsModel, utils
from hydromt import DataCatalog

import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
from matplotlib.colors import PowerNorm
from matplotlib.colors import BoundaryNorm
import matplotlib.gridspec as gridspec
from matplotlib.cm import ScalarMappable
import matplotlib.patheffects as path_effects
from matplotlib.colors import TwoSlopeNorm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize

from pyproj import Transformer
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import contextily as ctx
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from shapely.geometry import box
from rasterio.features import shapes
from shapely.geometry import shape

import warnings
warnings.filterwarnings('ignore')

import platform
prefix = "p:/" if platform.system() == "Windows" else "/p/"

def lat_formatter(x, pos):
    direction = 'N' if x >= 0 else 'S'
    return f"{abs(x):.1f}°{direction}"

def lon_formatter(x, pos):
    direction = 'E' if x >= 0 else 'W'
    return f"{abs(x):.1f}°{direction}"

def custom_formatter(value, pos=None):
    return f"{value:.1f}°"

#%%
print("Loading Factual model")
# define model and data catalog file paths
factual_model_dir = os.path.join(prefix,"11210471-001-compass","03_Runs","sofala","Idai","sfincs","event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0")
cf_model_dir = os.path.join(prefix,"11210471-001-compass","03_Runs","sofala","Idai","sfincs","event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10")

if platform.system() == "Windows":
    datacat = [
        '../../Workflows/03_data_catalogs/datacatalog_general.yml'
        ]
else:
    datacat = [
        '../../Workflows/03_data_catalogs/datacatalog_general___linux.yml'
        ]

data_catalog = DataCatalog(data_libs = datacat)

#%%
# Load in model, model region, and buffer model region
mod = SfincsModel(factual_model_dir, data_libs=datacat, mode="r")

hisfile = os.path.join(factual_model_dir,"sfincs_his.nc")
ds_his = xr.open_dataset(hisfile, engine='netcdf4')

model_region_gdf = gpd.read_file(join(
    prefix, "11210471-001-compass", "03_Runs", "sofala", "Idai", "sfincs", 
    "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0", "gis", "region.geojson"
)).to_crs("EPSG:4326") 

#%%
print("Masking permanent water")
# we set a threshold to mask minimum flood depth
hmin = 0.05

# compute the maximum over all time steps
da_zsmax = mod.results["zsmax"].max(dim="timemax")
     
# downscale the floodmap
depfile  = join(factual_model_dir, "subgrid", "dep_subgrid.tif")
da_dep   = mod.data_catalog.get_rasterdataset(depfile)

da_hmax = utils.downscale_floodmap(
      zsmax=da_zsmax,
      dep=da_dep,
      hmin=hmin,
      )
    
# GSWO dataset to mask permanent water by first geprojecting it to the subgrid of hmax
sfincs_region = mod.region
projection    = mod.crs.to_epsg()
gwso_region   = data_catalog.get_rasterdataset("gswo", geom=sfincs_region, buffer=1000)
gswo_mask     = gwso_region.raster.reproject_like(da_hmax, method="max")
# permanent water where water occurence > 5%
da_hmax_masked = da_hmax.where(gswo_mask <= 5)

# Add the name attribute for identification
mod.results['hmax'] = da_hmax
mod.results['hmax_masked'] = da_hmax_masked

# Make own background shape
valid_mask = (gswo_mask <= 5).astype("uint8").squeeze()

# Extract shapes
shapes_gen = shapes(valid_mask.values, transform=valid_mask.rio.transform())
valid_polygons = [shape(geom) for geom, val in shapes_gen if val == 1]
gdf_valid = gpd.GeoDataFrame(geometry=valid_polygons, crs=gswo_mask.rio.crs)
gdf_valid = gdf_valid.to_crs(model_region_gdf.crs)

del da_hmax, da_zsmax, da_dep  # Clean up to free memory


# #%%
# print("Plotting factual hmax masked")
# fig, ax = plt.subplots(nrows=1,
#                        ncols=1,
#                        figsize=(4,6),
#                        dpi=300,
#                        constrained_layout=True,
#                        subplot_kw={"projection": ccrs.epsg(projection)})

# # fig.suptitle("Factual Max Flood Depth", fontsize=11)

# # Plot hmax masked
# hmax = mod.results['hmax_masked']

# vmin = float(hmax.min())
# vmax = float(hmax.max())

# im = hmax.plot.pcolormesh(
#     ax=ax, cmap="viridis", vmin=vmin, vmax=3.5, add_colorbar=False
# )

# gdf_valid.plot(ax=ax, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)

# # Add model region
# model_region_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())

# # Set extent (based on actual lat/lon coordinates)
# minx, miny, maxx, maxy = model_region_gdf.bounds.minx.item(), model_region_gdf.bounds.miny.item(), model_region_gdf.bounds.maxx.item(), model_region_gdf.bounds.maxy.item()
# ax.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

# # Add gridlines and format tick labels
# gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
# gl.xlocator = mticker.FixedLocator(np.arange(minx, maxx + 0.1, 0.2))
# gl.ylocator = mticker.FixedLocator(np.arange(miny, maxy + 0.1, 0.2))
# gl.xformatter = mticker.FuncFormatter(lon_formatter)
# gl.yformatter = mticker.FuncFormatter(lat_formatter)
# gl.right_labels = False
# gl.top_labels = False
# gl.xlabel_style = {'size': 9}
# gl.ylabel_style = {'size': 9}

# # Titles
# ax.set_title("Factual Max Flood Depth", fontsize=11)

# # Add shared colorbar with larger size and labels
# cbar = fig.colorbar(im, ax=ax, orientation="vertical", 
#                     fraction=0.035, pad=0.05)
# cbar.set_label('Flood depth (m)', labelpad=5, fontsize=9)
# cbar.ax.tick_params(labelsize=8)
    
# fig.savefig("../figures/factual_hmax_masked.png", bbox_inches='tight', dpi=300)
# plt.show()

# #%%
# ### PLOT TIMESERIES OUTPUT POINTS
# print("Plotting timeseries")
# fig, ax = plt.subplots(3,1,figsize=(4,6), sharex=True, dpi=300)
# ax[0].plot(ds_his.time,ds_his['point_zs'].isel(stations=5),color='r', label=f'Station 5')
# ax[0].plot(mod.forcing['bzs'].time, mod.forcing['bzs'].sel(index=40), label=f'Station 40')

# ax[1].plot(mod.forcing['dis'].time, mod.forcing['dis'].sel(index=1), label=f'Gauge 1')
# ax[1].plot(mod.forcing['dis'].time, mod.forcing['dis'].sel(index=2), label=f'Gauge 2')

# ax[2].step(mod.forcing['precip_2d'].time, mod.forcing['precip_2d'].sum(dim=["x", "y"]), where='post')

# for a in ax:
#     a.set_xlim(ds_his.time.min(), ds_his.time.max())
#     a.legend(fontsize=9, loc="upper right") 
#     a.tick_params(axis='both', labelsize=10)
#     a.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.7)

# ax[0].set_title(f"Timeseries of water levels _base run")
# # ax.set_title(f"Timeseries of water levels \n Location: {ds_his['station_id'].isel(stations=5).values}")
# ax[0].set_ylabel("Water level height \n [m]")
# ax[1].set_ylabel("Discharge \n [m3/s]")
# ax[2].set_ylabel("Accumulated precipitation \n [mm/h] [m]")
# ax[2].xaxis.set_major_formatter(mdates.DateFormatter('%d'))
# ax[2].set_xlabel("Day in March 2019")
# fig.savefig("../figures/sfincs_output_TS.png", bbox_inches='tight', dpi=300)
# fig.show()


#%%
# print("Plotting factual hmax masked and timeseries for specific gauges/stations")
# gauges_list   = [1,2]
# stations_list = [5, 40]

# fig = plt.figure(figsize=(12, 5), dpi=300, constrained_layout=True)

# projection = mod.crs.to_epsg()
# colors = plt.get_cmap('tab10').colors  # Tuple of 10 colors

# gs = gridspec.GridSpec(nrows=3, ncols=2, figure=fig,
#                        width_ratios=[2, 1],  # Left col 3x wider than right col
#                        height_ratios=[1, 1, 1],  # Three rows on right all equal height
#                        wspace=0.05, hspace=0.03)

# ax_map = fig.add_subplot(gs[:, 0], projection=ccrs.PlateCarree())

# # Plot hmax masked
# utm_crs = ccrs.UTM(zone=36, southern_hemisphere=True)
# im = mod.results['hmax_masked'].plot.pcolormesh(
#     ax=ax_map, cmap="viridis", vmin=0, vmax=3.5, add_colorbar=False, transform=utm_crs, rasterized=True)

# ax_map.set_title("Factual Max Flood Depth", fontsize=10)

# gdf_valid.plot(ax=ax_map, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)

# # Add model region
# model_region_gdf.boundary.plot(ax=ax_map, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())

# # Set extent (based on actual lat/lon coordinates)
# minx, miny, maxx, maxy = model_region_gdf.bounds.minx.item(), model_region_gdf.bounds.miny.item(), model_region_gdf.bounds.maxx.item(), model_region_gdf.bounds.maxy.item()
# ax_map.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

# # Add gridlines and format tick labels
# gl = ax_map.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
# gl.xlocator = mticker.FixedLocator(np.arange(minx, maxx + 0.1, 0.2))
# gl.ylocator = mticker.FixedLocator(np.arange(miny, maxy + 0.1, 0.2))
# gl.xformatter = mticker.FuncFormatter(lon_formatter)
# gl.yformatter = mticker.FuncFormatter(lat_formatter)
# gl.right_labels = False
# gl.top_labels = False
# gl.xlabel_style = {'size': 9}
# gl.ylabel_style = {'size': 9}

# # ==== Plot city and river names ====
# # Plot Beira location
# ax_map.plot(34.848, -19.832, marker='o', color='black', markersize=4, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
# text = ax_map.text(34.852, -19.81, "Beira", transform=ccrs.PlateCarree(),
#                     fontsize=8, color='black', zorder=5)
# text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])

# # Buzi River marker and label
# ax_map.plot(34.43, -19.89, marker='o', color='black', markersize=4, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
# text2 = ax_map.text(34.44, -19.87, "Buzi River", transform=ccrs.PlateCarree(),
#                 fontsize=8, color='black', zorder=5)
# text2.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
# # Pungwe River marker and label
# ax_map.plot(34.543, -19.545, marker='o', color='black', markersize=4, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
# text3 = ax_map.text(34.554, -19.52, "Pungwe River", transform=ccrs.PlateCarree(),
#                 fontsize=8, color='black', zorder=5)
# text3.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])

# # Colorbar for map
# cbar = fig.colorbar(im, ax=ax_map, orientation="vertical", shrink=0.6, pad=0.01)
# cbar.set_label('Flood Depth (m)', labelpad=10, fontsize=9)
# cbar.ax.tick_params(labelsize=9)

# # RIGHT: three stacked time series subplots
# ax0 = fig.add_subplot(gs[0, 1])
# ax1 = fig.add_subplot(gs[1, 1], sharex=ax0)
# ax2 = fig.add_subplot(gs[2, 1], sharex=ax0)

# # Plotting timeseries
# ax0.plot(ds_his.time, ds_his['point_zs'].isel(stations=stations_list[0]), color=colors[4], label=f'S{stations_list[0]}')
# ax0.plot(mod.forcing['bzs'].time, mod.forcing['bzs'].sel(index=stations_list[1]), color=colors[1], label=f'S{stations_list[1]}')
# ax0.set_ylabel("Water level height \n [m]", fontsize=9)
# ax0.legend(fontsize=8, loc="upper right", markerscale=0.8)
# ax0.grid(True, linestyle='--', alpha=0.6)
# ax0.tick_params(labelsize=9)
# ax0.set_xlim(ds_his.time.min(), ds_his.time.max())
# ax0.set_title("Factual Time Series", fontsize=10)
# ax0.tick_params(bottom=False, labelbottom=False)

# ax1.plot(mod.forcing['dis'].time, mod.forcing['dis'].sel(index=gauges_list[0]), color=colors[2], label=f'G{gauges_list[0]}')
# ax1.plot(mod.forcing['dis'].time, mod.forcing['dis'].sel(index=gauges_list[1]), color=colors[3], label=f'G{gauges_list[1]}')
# ax1.set_ylabel("Discharge \n [m³/s]", fontsize=9)
# ax1.legend(fontsize=8, loc="upper right", markerscale=0.8)
# ax1.grid(True, linestyle='--', alpha=0.6)
# ax1.tick_params(labelsize=9)
# ax1.set_xlim(ds_his.time.min(), ds_his.time.max())
# ax1.tick_params(bottom=False, labelbottom=False)

# ax2.step(mod.forcing['precip_2d'].time, mod.forcing['precip_2d'].sum(dim=["x", "y"]), where='post', color=colors[0])
# ax2.set_ylabel("Accum. precipitation \n [mm/h]", fontsize=10)
# ax2.set_xlabel("Day in March 2019")
# ax2.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
# ax2.grid(True, linestyle='--', alpha=0.6)
# ax2.tick_params(labelsize=9)
# ax2.set_xlim(ds_his.time.min(), ds_his.time.max())

# # Plot station and Gauge points on the flood map and annotate
# ax_map.scatter(ds_his['point_zs'].isel(stations=stations_list[0]).point_x.values, ds_his['point_zs'].isel(stations=stations_list[0]).point_y.values, color=colors[4], s=30, edgecolor='k', zorder=5, transform=ccrs.epsg(projection), label=f'Station {stations_list[0]}')
# ax_map.scatter(mod.forcing['bzs'].sel(index=stations_list[1]).geometry.item().x, mod.forcing['bzs'].sel(index=stations_list[1]).geometry.item().y, color=colors[1], s=30, edgecolor='k', zorder=5, transform=ccrs.epsg(projection), label=f'Station {stations_list[1]}')

# ax_map.scatter(mod.forcing['dis'].sel(index=gauges_list[0]).geometry.item().x, mod.forcing['dis'].sel(index=gauges_list[0]).geometry.item().y, color=colors[2], s=30, edgecolor='k', zorder=5, transform=ccrs.epsg(projection), label=f'Gauge {gauges_list[0]}')
# ax_map.scatter(mod.forcing['dis'].sel(index=gauges_list[1]).geometry.item().x, mod.forcing['dis'].sel(index=gauges_list[1]).geometry.item().y, color=colors[3], s=30, edgecolor='k', zorder=5, transform=ccrs.epsg(projection), label=f'Gauge {gauges_list[1]}')

# offset = 2000  # Adjust offset depending on map scale (in map projection units)
# ax_map.annotate(f'S{stations_list[0]}', xy=(ds_his['point_zs'].isel(stations=stations_list[0]).point_x.values, ds_his['point_zs'].isel(stations=stations_list[0]).point_y.values),
#                 xytext=(ds_his['point_zs'].isel(stations=stations_list[0]).point_x.values + 2000, ds_his['point_zs'].isel(stations=stations_list[0]).point_y.values - 4000),
#                 transform=ccrs.epsg(projection),
#                 fontsize=8, color=colors[4], fontweight='bold',
#                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=0.5, alpha=0.8)
# )
# ax_map.annotate(f'S{stations_list[1]}', xy=(mod.forcing['bzs'].sel(index=stations_list[1]).geometry.item().x, mod.forcing['bzs'].sel(index=stations_list[1]).geometry.item().y),
#                 xytext=(mod.forcing['bzs'].sel(index=stations_list[1]).geometry.item().x + 2000, mod.forcing['bzs'].sel(index=stations_list[1]).geometry.item().y - 4000),
#                 transform=ccrs.epsg(projection),
#                 fontsize=8, color=colors[1], fontweight='bold',
#                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=0.5, alpha=0.8)
# )
# ax_map.annotate(f'G{gauges_list[0]}', xy=(mod.forcing['dis'].sel(index=gauges_list[0]).geometry.item().x, mod.forcing['dis'].sel(index=gauges_list[0]).geometry.item().y),
#                 xytext=(mod.forcing['dis'].sel(index=gauges_list[0]).geometry.item().x - 6500, mod.forcing['dis'].sel(index=gauges_list[0]).geometry.item().y - 4000),
#                 transform=ccrs.epsg(projection),
#                 fontsize=8, color=colors[2], fontweight='bold',
#                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=0.5, alpha=0.8)
# )
# ax_map.annotate(f'G{gauges_list[1]}', xy=(mod.forcing['dis'].sel(index=gauges_list[1]).geometry.item().x, mod.forcing['dis'].sel(index=gauges_list[1]).geometry.item().y),
#                 xytext=(mod.forcing['dis'].sel(index=gauges_list[1]).geometry.item().x + 2500, mod.forcing['dis'].sel(index=gauges_list[1]).geometry.item().y),
#                 transform=ccrs.epsg(projection),
#                 fontsize=8, color=colors[3], fontweight='bold',
#                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=0.5, alpha=0.8)
# )
# props = dict( facecolor='w', lw=0, alpha=0.8)
# ax_map.text(0.03, 0.95, '(a)', transform=ax_map.transAxes, bbox=props, fontsize=10, fontweight='bold')
# ax0.text(0.02, 0.85, '(b)', transform=ax0.transAxes, bbox=props, fontsize=10, fontweight='bold')
# ax1.text(0.02, 0.85, '(c)', transform=ax1.transAxes, bbox=props, fontsize=10, fontweight='bold')
# ax2.text(0.02, 0.85, '(d)', transform=ax2.transAxes, bbox=props, fontsize=10, fontweight='bold')


# fig.savefig("../figures/factual_hmax_timeseries.png", bbox_inches='tight', dpi=300)
# fig.savefig("../figures/factual_hmax_timeseries.pdf", bbox_inches='tight', dpi=300)
# plt.show()





#%% ############################################
# ========== Factual damage plotting ===========
################################################
# Define conversion factor from 2010 euros to 2019 USD
eur_to_usd = 1.326 #
usd_2010_to_2019 = 1.152 

# Base paths - update these as needed
BASE_RUN_PATH = Path(os.path.join(prefix,"11210471-001-compass","03_Runs","sofala"))
OUTPUT_DIR = Path("../figures")

# ===== FILE PATHS =====
# buildings with no damage have already been masked out
file_cf0 = BASE_RUN_PATH / "Idai" / "fiat" / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" / "output" / "output_relative_damage.fgb"
file_cfall = BASE_RUN_PATH / "Idai" / "fiat" / "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10" / "output" / "output_relative_damage.fgb"

model_region_gdf = gpd.read_file(join(prefix, "11210471-001-compass", "03_Runs", "sofala", "Idai", "sfincs", 
                                      "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0", 
                                      "gis", "region.geojson")).to_crs("EPSG:4326") 

beira_region = gpd.read_file(join(prefix, "11210471-001-compass", "01_Data", "sofala_geoms", "Beira_region.shp")
                             ).to_crs("EPSG:4326")

# Create output directory if it doesn't exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ===== LOAD DAMAGE DATA =====
print("Loading damage data files...")
print(f"Loading CF0: {file_cf0}")

# Read the geodataframes
gdf_cf0 = gpd.read_file(file_cf0)
gdf_cfall = gpd.read_file(file_cfall)

#%%
# ===== EXTRACT COORDINATES =====
print("Extracting coordinates from building centroids...")
# Coordinate system: PlateCarree = lat/lon
crs = ccrs.PlateCarree()

# Extract x, y coordinates from geometry centroids (for polygon buildings)
gdf_cf0['centroid'] = gdf_cf0.geometry.centroid
gdf_cf0['x'] = gdf_cf0['centroid'].x
gdf_cf0['y'] = gdf_cf0['centroid'].y
gdf_cf0['total_damage_USD'] = gdf_cf0['total_damage'] * eur_to_usd * usd_2010_to_2019
gdf_cf0['max_total_damage_USD'] = gdf_cf0['max_damage_total'] * eur_to_usd * usd_2010_to_2019

gdf_cf0 = gdf_cf0[['object_id', 'x', 'y', 'total_damage_USD', 'max_total_damage_USD']]

# Remove buildings with no coordinates
gdf_cf0 = gdf_cf0.dropna(subset=['x', 'y'])

print(gdf_cf0[['object_id', 'x', 'y', 'total_damage_USD', 'max_total_damage_USD']].head())

# Convert to GeoDataFrame if not already
gdf_cf0_crs = gpd.GeoDataFrame(gdf_cf0,
                               geometry=gpd.points_from_xy(gdf_cf0["x"], gdf_cf0["y"]),
                               crs="EPSG:32736"  # or match whatever CRS your data is in
                               )

# Ensure same CRS
gdf_damage = gdf_cf0_crs.to_crs(crs)


# same for the CF
# Extract x, y coordinates from geometry centroids (for polygon buildings)
gdf_cfall['centroid'] = gdf_cfall.geometry.centroid
gdf_cfall['x'] = gdf_cfall['centroid'].x
gdf_cfall['y'] = gdf_cfall['centroid'].y
gdf_cfall['total_damage_USD'] = gdf_cfall['total_damage'] * eur_to_usd * usd_2010_to_2019
gdf_cfall['max_total_damage_USD'] = gdf_cfall['max_damage_total'] * eur_to_usd * usd_2010_to_2019

gdf_cfall = gdf_cfall[['object_id', 'x', 'y', 'total_damage_USD', 'max_total_damage_USD']]

# Remove buildings with no coordinates
gdf_cfall = gdf_cfall.dropna(subset=['x', 'y'])

print(gdf_cfall[['object_id', 'x', 'y', 'total_damage_USD', 'max_total_damage_USD']].head())

# Convert to GeoDataFrame if not already
gdf_cfall_crs = gpd.GeoDataFrame(gdf_cfall,
                               geometry=gpd.points_from_xy(gdf_cfall["x"], gdf_cfall["y"]),
                               crs="EPSG:32736"  # or match whatever CRS your data is in
                               )

# Ensure same CRS
gdf_cf_damage = gdf_cfall_crs.to_crs(crs)


#%%
# ===== STATISTICS =====
print(f"\nDamage Statistics for Idai (total_damage_USD):")
cf0_damage = gdf_cf0['total_damage_USD']
cfall_damage = gdf_cfall['total_damage_USD']

# Select rows where total_damage == max_damage_total
mask_equal = np.isclose(gdf_cf0["total_damage_USD"], gdf_cf0["max_total_damage_USD"], rtol=1e-5)
gdf_equal = gdf_cf0[mask_equal]

# Select rows where total_damage is at least 10% of max_damage_total
df_10pct_or_more = gdf_cf0[gdf_cf0["total_damage_USD"] >= 0.1 * gdf_cf0["max_total_damage_USD"]]

# Calculate total damage within Beira region
damage_in_beira = gpd.sjoin(gdf_damage, beira_region, how="inner", predicate="intersects")
total_damage_beira = damage_in_beira['total_damage_USD'].sum()
total_damage = cf0_damage.sum()
perct_dam_beira = (total_damage_beira / total_damage) * 100
print(f"Total damage in Beira region: ${total_damage_beira:.0f}, {perct_dam_beira:.0f}% of total")

print(f"CF0 max damage: ${cf0_damage.max():.0f}")
print(f"CF0 mean damage: ${cf0_damage.mean():.0f}")
print(f"CF0 total damage: ${cf0_damage.sum():.0f}")

print(f"# buildings damaged: {len(cf0_damage[cf0_damage>0])}")
print(f"# buildings totally destroyed: {len(gdf_equal)}")
print(f"# buildings >10% destroyed: {len(df_10pct_or_more)}")

print(f"CFall max damage: ${cfall_damage.max():.0f}")
print(f"CFall mean damage: ${cfall_damage.mean():.0f}")
print(f"CFall total damage: ${cfall_damage.sum():.0f}")



#%%
# === PLOTTING ===
# Create damage plotting settings
max_total = cf0_damage.quantile(0.9)
boundaries = np.linspace(0, max_total, 11)
damage_norm = BoundaryNorm(boundaries, ncolors=256, clip=True)
damage_cmap = plt.get_cmap('Reds')
damage_label = 'Total Damage [USD]'
# Plot settings
point_size = 10
alpha = 0.8

#%%
# print("Creating scattered factual damage plot...")
# fig, axes = plt.subplots(1, 1, figsize=(4, 6), subplot_kw={'projection': crs})


# # Plot damage points
# scatter1 = axes.scatter(gdf_damage.geometry.x, gdf_damage.geometry.y, c=gdf_damage['total_damage_USD'], 
#                         cmap=damage_cmap, norm=damage_norm, s=point_size, alpha=alpha, transform=crs)

# # Add model region
# model_region_gdf.boundary.plot(ax=axes, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())

# # # Add background and set extent (based on actual lat/lon coordinates)
# gdf_valid.plot(ax=axes, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)
# minx, miny, maxx, maxy = model_region_gdf.bounds.minx.item(), model_region_gdf.bounds.miny.item(), model_region_gdf.bounds.maxx.item(), model_region_gdf.bounds.maxy.item()
# axes.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

# # Add gridlines and format tick labels
# gl = axes.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
# gl.xlocator = mticker.FixedLocator(np.arange(minx, maxx + 0.1, 0.2))
# gl.ylocator = mticker.FixedLocator(np.arange(miny, maxy + 0.1, 0.2))
# gl.xformatter = mticker.FuncFormatter(lon_formatter)
# gl.yformatter = mticker.FuncFormatter(lat_formatter)
# gl.right_labels = False
# gl.top_labels = False
# gl.xlabel_style = {'size': 9}
# gl.ylabel_style = {'size': 9}


# # Add colorbar
# cbar1 = plt.colorbar(scatter1, ax=axes, shrink=0.8, pad=0.06)
# cbar1.set_label(damage_label, fontsize=10)

# # Titles and save
# plt.tight_layout()
# plt.suptitle(f'Factual Total Damage', fontsize=16, y=0.96)

# output_file_main = OUTPUT_DIR / f'damage_idai_total_damage_scatter.png'
# plt.savefig(output_file_main, dpi=300, bbox_inches='tight')
# plt.show()



# %% ##############################################################
# ============= Spatially aggregate damage damage  ================
###################################################################
print("Spatially aggregating damage data")
# clip permanent water shapes to model region
clipped_region = gpd.overlay(model_region_gdf, gdf_valid, how="intersection")

# Set resolution in degrees or meters (depending on CRS)
cell_size = 0.025  # e.g., in degrees for EPSG:4326

# Bounds
minx, miny, maxx, maxy = clipped_region.total_bounds
cols = np.arange(minx, maxx, cell_size)
rows = np.arange(miny, maxy, cell_size)

# Build grid cells
grid_cells = [box(x, y, x + cell_size, y + cell_size) for x in cols for y in rows]
gdf_grid = gpd.GeoDataFrame(geometry=grid_cells, crs=clipped_region.crs)

# Clip again to final shape
gdf_grid_masked = gpd.overlay(gdf_grid, clipped_region, how='intersection')
gdf_cf_grid_masked = gpd.overlay(gdf_grid, clipped_region, how='intersection')

# %%  
# Spatial join: assign grid cell index to each point
joined = gpd.sjoin(gdf_damage, gdf_grid_masked, how="left", predicate="within")

# Aggregate damage per grid cell
agg_tot_damage = joined.groupby(joined.index_right)["total_damage_USD"].sum()
agg_max_damage = joined.groupby(joined.index_right)["max_total_damage_USD"].sum()

# Assign damage to grid GeoDataFrame
gdf_grid_masked["total_damage_USD"] = agg_tot_damage
gdf_grid_masked["max_total_damage_USD"] = agg_max_damage
gdf_grid_masked["total_damage_USD"] = gdf_grid_masked["total_damage_USD"].fillna(0)
gdf_grid_masked['total_damage_M']   = gdf_grid_masked["total_damage_USD"]/1e6
gdf_grid_masked["max_total_damage_USD"] = gdf_grid_masked["max_total_damage_USD"].fillna(0)

# Relative percentual damage per grid
gdf_grid_masked['relative_aggr_damage'] = (agg_tot_damage / agg_max_damage) * 100 # %
gdf_grid_masked["relative_aggr_damage"] = gdf_grid_masked["relative_aggr_damage"].fillna(0)

# Same for CF data
# Spatial join: assign grid cell index to each point
joined = gpd.sjoin(gdf_cf_damage, gdf_cf_grid_masked, how="left", predicate="within")

# Aggregate damage per grid cell
agg_tot_cf_damage = joined.groupby(joined.index_right)["total_damage_USD"].sum()
agg_max_cf_damage = joined.groupby(joined.index_right)["max_total_damage_USD"].sum()

# Assign damage to grid GeoDataFrame
gdf_cf_grid_masked["total_damage_USD"]     = agg_tot_cf_damage
gdf_cf_grid_masked["max_total_damage_USD"] = agg_max_cf_damage
gdf_cf_grid_masked["total_damage_USD"]     = gdf_cf_grid_masked["total_damage_USD"].fillna(0)
gdf_cf_grid_masked['total_damage_M']       = gdf_cf_grid_masked["total_damage_USD"]/1e6
gdf_cf_grid_masked["max_total_damage_USD"] = gdf_cf_grid_masked["max_total_damage_USD"].fillna(0)

# Relative percentual damage per grid
gdf_cf_grid_masked['relative_aggr_damage']      = (agg_tot_cf_damage / agg_max_cf_damage) * 100 # %
gdf_cf_grid_masked["relative_aggr_damage"]      = gdf_cf_grid_masked["relative_aggr_damage"].fillna(0)
gdf_cf_grid_masked["relative_aggr_damage_diff"] = ((agg_tot_damage - agg_tot_cf_damage) / agg_max_damage) * 100 # %
gdf_cf_grid_masked["relative_aggr_damage_diff"] = gdf_cf_grid_masked["relative_aggr_damage_diff"].fillna(0)

# %% ##############################################################
# ======= Plot the factual aggregated damage and flooding =========
###################################################################
# print("Plotting spatially aggregated total damage")
# # plot the total_damage_SD, emphazizing lower values
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,6), dpi=300,
#                        subplot_kw={"projection": ccrs.PlateCarree()})

# # Create colormap normalization
# norm = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_grid_masked['total_damage_USD'].max())

# # Plot using GeoPandas, but draw to custom ax and return colorbar mappable
# plot = gdf_grid_masked.plot(column="total_damage_USD", cmap="Reds", norm=norm, edgecolor="grey",
#                             linewidth=0.2, ax=ax, legend=False)

# # Add colorbar
# sm = ScalarMappable(norm=norm, cmap="Reds")
# cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.06)
# cbar.set_label(damage_label, fontsize=10) 

# # Add model region
# model_region_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())

# # # Add background and set extent (based on actual lat/lon coordinates)
# gdf_valid.plot(ax=ax, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)
# minx, miny, maxx, maxy = model_region_gdf.bounds.minx.item(), model_region_gdf.bounds.miny.item(), model_region_gdf.bounds.maxx.item(), model_region_gdf.bounds.maxy.item()
# ax.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

# # Add gridlines and format tick labels
# gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
# gl.xlocator = mticker.FixedLocator(np.arange(minx, maxx + 0.1, 0.2))
# gl.ylocator = mticker.FixedLocator(np.arange(miny, maxy + 0.1, 0.2))
# gl.xformatter = mticker.FuncFormatter(lon_formatter)
# gl.yformatter = mticker.FuncFormatter(lat_formatter)
# gl.right_labels = False
# gl.top_labels = False
# gl.xlabel_style = {'size': 9}
# gl.ylabel_style = {'size': 9}

# ax.set_title("Total Aggregated Flood Damage", fontsize=12)
# fig.savefig("../figures/total_aggregated_damage.png", bbox_inches='tight', dpi=300)

# plt.show()



# %%
# print("Plotting spatially relative total damage")
# # plot the total_damage, emphazizing lower values
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,6), dpi=300,
#                        subplot_kw={"projection": ccrs.PlateCarree()})

# # Plot using GeoPandas, but draw to custom ax and return colorbar mappable
# plot = gdf_grid_masked.plot(column="relative_aggr_damage", cmap="Reds", edgecolor="grey",
#                             linewidth=0.2, ax=ax, legend=False, vmin = 0, vmax = 100)

# # Add colorbar
# sm = ScalarMappable(cmap="Reds", norm=plt.Normalize(vmin=0, vmax=100))
# sm.set_array([])  # doesn't need real data; just placeholder for colorbar
# cbar = fig.colorbar(sm, ax=ax, shrink=0.8, pad=0.06)
# cbar.set_label("Relative Damage [%]", fontsize=10) 

# # Add model region
# model_region_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())

# # # Add background and set extent (based on actual lat/lon coordinates)
# gdf_valid.plot(ax=ax, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)
# minx, miny, maxx, maxy = model_region_gdf.bounds.minx.item(), model_region_gdf.bounds.miny.item(), model_region_gdf.bounds.maxx.item(), model_region_gdf.bounds.maxy.item()
# ax.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

# # Add gridlines and format tick labels
# gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
# gl.xlocator = mticker.FixedLocator(np.arange(minx, maxx + 0.1, 0.2))
# gl.ylocator = mticker.FixedLocator(np.arange(miny, maxy + 0.1, 0.2))
# gl.xformatter = mticker.FuncFormatter(lon_formatter)
# gl.yformatter = mticker.FuncFormatter(lat_formatter)
# gl.right_labels = False
# gl.top_labels = False
# gl.xlabel_style = {'size': 9}
# gl.ylabel_style = {'size': 9}

# ax.set_title("Relative damage Map", fontsize=12)
# fig.savefig("../figures/relative_aggregated_damage.png", bbox_inches='tight', dpi=300)

# plt.show()


# %%
# Plot the damage and flooding as sub panels
# print("plotting factual flooding and aggregated total damage")
# fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 6), dpi=300, constrained_layout=True, 
#                          subplot_kw={"projection": ccrs.PlateCarree()})

# # Plot the flooding
# utm_crs = ccrs.UTM(zone=36, southern_hemisphere=True)
# hmax = mod.results['hmax_masked'].load()
# im = hmax.plot.pcolormesh(ax=axes[0], cmap="viridis", vmin=0, vmax=3.5, add_colorbar=False, transform=utm_crs, rasterized=True)

# # Plot the total damage
# norm = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_grid_masked['total_damage_M'].max())
# gdf_grid_masked[gdf_grid_masked['total_damage_M'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)

# plot = gdf_grid_masked[gdf_grid_masked['total_damage_M'] > 0].plot(column='total_damage_M', cmap='Reds', norm=norm, edgecolor='grey', 
#                                                                  linewidth=0.2, ax=axes[1], legend=False, zorder=2)

# background = gdf_valid.to_crs("EPSG:4326")  # Do once
# region_boundary = model_region_gdf.to_crs("EPSG:4326")
# subplot_labels = ['(a)', '(b)']

# for i, ax in enumerate(axes):
#     # Add model region
#     region_boundary.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)

#     # # Add background and set extent (based on actual lat/lon coordinates)
#     background.plot(ax=ax, color='#E0E0E0', zorder=0)

#     # Add gridlines and format tick labels
#     gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
#     gl.xlocator = mticker.FixedLocator(np.arange(minx, maxx + 0.1, 0.2))
#     gl.ylocator = mticker.FixedLocator(np.arange(miny, maxy + 0.1, 0.2))
#     gl.xformatter = mticker.FuncFormatter(lon_formatter)
#     gl.yformatter = mticker.FuncFormatter(lat_formatter)
#     gl.right_labels = False
#     gl.top_labels = False
#     gl.xlabel_style = {'size': 9}
#     gl.ylabel_style = {'size': 9}
#     if i == 1:  
#         gl.left_labels = False  # disable y-axis labels
            
#     # ==== Plot city and river names ====
#     # Plot Beira location
#     ax.plot(34.848, -19.832, marker='o', color='black', markersize=4, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
#     text = ax.text(34.84, -19.89, "Beira", transform=ccrs.PlateCarree(),
#                         fontsize=8, color='black', zorder=5)
#     text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
    
#     # Buzi River marker and label
#     ax.plot(34.43, -19.89, marker='o', color='black', markersize=4, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
#     text2 = ax.text(34.44, -19.87, "Buzi River", transform=ccrs.PlateCarree(),
#                     fontsize=8, color='black', zorder=5)
#     text2.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])

#     # Pungwe River marker and label
#     ax.plot(34.543, -19.545, marker='o', color='black', markersize=4, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
#     text3 = ax.text(34.554, -19.52, "Pungwe River", transform=ccrs.PlateCarree(),
#                     fontsize=8, color='black', zorder=5)
#     text3.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])

#     ax.text(0, 1.02, subplot_labels[i], transform=ax.transAxes,
#             fontsize=10, fontweight='bold', va='bottom', ha='left')

# minx, miny, maxx, maxy = region_boundary.bounds.minx.item(), region_boundary.bounds.miny.item(), region_boundary.bounds.maxx.item(), region_boundary.bounds.maxy.item()
# for ax in axes:
#     ax.set_extent([minx, maxx, miny, maxy], crs=ccrs.PlateCarree())

# # Titles
# axes[0].set_title("", fontsize=10)
# axes[1].set_title("", fontsize=10)

# # ==== Colorbar for Flood Depth ====
# cbar1 = fig.colorbar(im, ax=axes[0], orientation="vertical", 
#                      fraction=0.035, aspect=20, pad=0.01)
# cbar1.set_label("Maximum flood depth (m)", labelpad=6, fontsize=9)
# cbar1.ax.tick_params(labelsize=8)

# # ==== Colorbar for Damage ====
# sm = ScalarMappable(norm=norm, cmap="Reds")
# sm.set_array([])  # Required to avoid warning, even if dummy
# cbar2 = fig.colorbar(sm, ax=axes[1], orientation="vertical", 
#                      fraction=0.035, aspect=20, pad=0.01)
# cbar2.set_label('Aggregated total damage [M USD]', labelpad=6, fontsize=9)
# cbar2.ax.tick_params(labelsize=8)
# # Make the 1e7 offset text smaller
# cbar2.ax.yaxis.offsetText.set_fontsize(7)


# fig.savefig("../figures/f03.png", bbox_inches='tight', dpi=300)
# fig.savefig("../figures/f03.pdf", bbox_inches='tight', dpi=300)

# %%
# print("plotting factual flooding and aggregated relative damage")
# fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 6), dpi=300, constrained_layout=True, 
#                          subplot_kw={"projection": ccrs.PlateCarree()})

# # Plot the flooding
# utm_crs = ccrs.UTM(zone=36, southern_hemisphere=True)
# hmax = mod.results['hmax_masked'].load()
# im = hmax.plot.pcolormesh(ax=axes[0], cmap="viridis", vmin=0, vmax=3.5, add_colorbar=False, transform=utm_crs, rasterized=True)

# gdf_grid_masked['plot_rel_agg_dam'] = gdf_grid_masked['relative_aggr_damage'].replace(0, np.nan)

# plot = gdf_grid_masked.plot(column='plot_rel_agg_dam', cmap='Reds', edgecolor="grey", linewidth=0.2, 
#                             ax=axes[1], legend=False, vmin=0, vmax=100, missing_kwds={"color": "white"}, rasterized=True)

# background = gdf_valid.to_crs("EPSG:4326") 
# region_boundary = model_region_gdf.to_crs("EPSG:4326")
# subplot_labels = ['(a)', '(b)']

# for i, ax in enumerate(axes):
#     # Add model region
#     region_boundary.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)

#     # # Add background and set extent (based on actual lat/lon coordinates)
#     background.plot(ax=ax, color='#E0E0E0', zorder=0)

#     # Add gridlines and format tick labels
#     gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
#     gl.xlocator = mticker.FixedLocator(np.arange(minx, maxx + 0.1, 0.2))
#     gl.ylocator = mticker.FixedLocator(np.arange(miny, maxy + 0.1, 0.2))
#     gl.xformatter = mticker.FuncFormatter(lon_formatter)
#     gl.yformatter = mticker.FuncFormatter(lat_formatter)
#     gl.right_labels = False
#     gl.top_labels = False
#     gl.xlabel_style = {'size': 9}
#     gl.ylabel_style = {'size': 9}
#     if i == 1:  
#         gl.left_labels = False  # disable y-axis labels
    
#     ax.text(0, 1.02, subplot_labels[i], transform=ax.transAxes,
#             fontsize=10, fontweight='bold', va='bottom', ha='left')

#     # ==== Plot city and river names ====
#     # Plot Beira location
#     ax.plot(34.848, -19.832, marker='o', color='black', markersize=4, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
#     text = ax.text(34.84, -19.89, "Beira", transform=ccrs.PlateCarree(),
#                         fontsize=8, color='black', zorder=5)
#     text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
    
#     # Buzi River marker and label
#     ax.plot(34.43, -19.89, marker='o', color='black', markersize=4, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
#     text2 = ax.text(34.44, -19.87, "Buzi River", transform=ccrs.PlateCarree(),
#                     fontsize=8, color='black', zorder=5)
#     text2.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])

#     # Pungwe River marker and label
#     ax.plot(34.543, -19.545, marker='o', color='black', markersize=4, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
#     text3 = ax.text(34.554, -19.52, "Pungwe River", transform=ccrs.PlateCarree(),
#                     fontsize=8, color='black', zorder=5)
#     text3.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])

# minx, miny, maxx, maxy = region_boundary.bounds.minx.item(), region_boundary.bounds.miny.item(), region_boundary.bounds.maxx.item(), region_boundary.bounds.maxy.item()
# for ax in axes:
#     ax.set_extent([minx, maxx, miny, maxy], crs=ccrs.PlateCarree())

# # Titles
# axes[0].set_title("")
# # axes[1].set_title("(b) Aggregated Relative Damage", fontsize=10)

# # ==== Colorbar for Flood Depth ====
# cbar1 = fig.colorbar(im, ax=axes[0], orientation="vertical", 
#                      fraction=0.035, aspect=20, pad=0.01)
# cbar1.set_label("Maximum Flood Depth (m)", labelpad=6, fontsize=9)
# cbar1.ax.tick_params(labelsize=8)

# # ==== Colorbar for Damage ====
# sm = ScalarMappable(cmap="Reds", norm=plt.Normalize(vmin=0, vmax=100))
# sm.set_array([])
# cbar2 = fig.colorbar(sm, ax=axes[1], orientation="vertical", 
#                      fraction=0.035, aspect=20, pad=0.01)
# cbar2.set_label("Relative Damage [%]", labelpad=6, fontsize=9)
# cbar2.ax.tick_params(labelsize=8)
# # Make the 1e7 offset text smaller
# cbar2.ax.yaxis.offsetText.set_fontsize(7)


# fig.savefig("../figures/factual_flooding_and_aggregated_rel_damage.png", bbox_inches='tight', dpi=300)
# fig.savefig("../figures/factual_flooding_and_aggregated_rel_damage.pdf", bbox_inches='tight', dpi=300)


#%%
print("Plotting spatially aggregated total damage for F, CF and diff")

gdf_cf_grid_masked['total_damage_diff'] = (gdf_grid_masked['total_damage_M'] - gdf_cf_grid_masked['total_damage_M'])

# plot the total_damage, emphazizing lower values
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 5), dpi=300, sharey=True, constrained_layout=True,
                       subplot_kw={"projection": ccrs.PlateCarree()})

# Create colormap normalization
norm = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_grid_masked['total_damage_M'].max())
norm_diff = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_cf_grid_masked["total_damage_diff"].max())
red_half = LinearSegmentedColormap.from_list("bwr_red", plt.cm.bwr(np.linspace(0.5, 1, 256)))

# Plot using GeoPandas, but draw to custom ax and return colorbar mappable
gdf_grid_masked[gdf_grid_masked['total_damage_M'] == 0].plot(ax=axes[0], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_grid_masked[gdf_grid_masked['total_damage_M'] > 0].plot(column='total_damage_M', cmap='Reds', norm=norm, edgecolor='grey', 
                                                                 linewidth=0.2, ax=axes[0], legend=False, zorder=2, rasterized=True)

gdf_cf_grid_masked[gdf_cf_grid_masked['total_damage_M'] == 0].plot(ax=axes[1], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_cf_grid_masked[gdf_cf_grid_masked['total_damage_M'] > 0].plot(column='total_damage_M', cmap='Reds', norm=norm, edgecolor='grey', 
                                                                 linewidth=0.2, ax=axes[1], legend=False, zorder=2, rasterized=True)

gdf_cf_grid_masked[gdf_cf_grid_masked['total_damage_diff'] <= 0].plot(ax=axes[2], color='white', edgecolor='grey', linewidth=0.2, zorder=1)
plot = gdf_cf_grid_masked[gdf_cf_grid_masked['total_damage_diff'] > 0].plot(column='total_damage_diff', cmap=red_half, norm=norm_diff, edgecolor='grey', 
                                                                 linewidth=0.2, ax=axes[2], legend=False, zorder=2, missing_kwds={"color": "white", "edgecolor": "none"}, rasterized=True)
subplot_labels = ['(a)', '(b)', '(c)']

for i, ax in enumerate(axes):
    # Add model region
    model_region_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())

    # # Add background and set extent (based on actual lat/lon coordinates)
    gdf_valid.plot(ax=ax, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)
    minx, miny, maxx, maxy = model_region_gdf.bounds.minx.item(), model_region_gdf.bounds.miny.item(), model_region_gdf.bounds.maxx.item(), model_region_gdf.bounds.maxy.item()
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
sm1 = ScalarMappable(norm=norm, cmap="Reds")
sm1.set_array([])
cbar = fig.colorbar(sm1, ax=axes[0:2], orientation="vertical", shrink=0.4, pad=0.01)
cbar.set_label("Aggregated Total Damage [M USD]", fontsize = 9)
cbar.ax.tick_params(labelsize=8)

norm = Normalize(vmin=0, vmax=gdf_cf_grid_masked["total_damage_diff"].max())
sm2 = ScalarMappable(cmap=red_half, norm=norm_diff)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, ax=axes[2], orientation="vertical", shrink=0.4, pad=0.01)
cbar2.set_label("Attributable Total Damage [M USD]", fontsize = 9)
cbar2.ax.tick_params(labelsize=8)

axes[0].set_title("Factual", fontsize=10)
axes[1].set_title("Counterfactual", fontsize=10)
axes[2].set_title("Factual - Counterfactual", fontsize=10)

# fig.suptitle("Total Aggregated Flood Damage", fontsize=12)
fig.savefig("../figures/f_and_cf_total_aggregated_damage.png", bbox_inches='tight', dpi=300)
fig.savefig("../figures/f_and_cf_total_aggregated_damage.pdf", bbox_inches='tight', dpi=300)

# %%
gdf_cf_grid_masked['rel_dam_diff'] = (gdf_grid_masked['relative_aggr_damage'] - gdf_cf_grid_masked['relative_aggr_damage']) / gdf_grid_masked['relative_aggr_damage'] * 100

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10,5), dpi=300, constrained_layout=True, 
                         subplot_kw={"projection": ccrs.PlateCarree()})

# Create white-to-red colormap
white_red = LinearSegmentedColormap.from_list("white_red", ["white", "red"])

# Set normalization
norm_diff = Normalize(vmin=0, vmax=15)

gdf_grid_masked['plot_rel_agg_dam'] = gdf_grid_masked['relative_aggr_damage'].replace(0, np.nan)
plot = gdf_grid_masked.plot(column='plot_rel_agg_dam', cmap='Reds', edgecolor="grey", linewidth=0.2, 
                            ax=axes[0], legend=False, vmin=0, vmax=100, missing_kwds={"color": "white"}, rasterized=True)

gdf_cf_grid_masked['plot_rel_agg_dam'] = gdf_cf_grid_masked['relative_aggr_damage'].replace(0, np.nan)
plot = gdf_cf_grid_masked.plot(column='plot_rel_agg_dam', cmap='Reds', edgecolor="grey", linewidth=0.2, 
                            ax=axes[1], legend=False, vmin=0, vmax=100, missing_kwds={"color": "white"}, rasterized=True)

gdf_cf_grid_masked['plot_rel_dam_diff'] = gdf_cf_grid_masked['relative_aggr_damage_diff'].replace(0, np.nan)
plot = gdf_cf_grid_masked.plot(column='plot_rel_dam_diff', cmap=white_red, edgecolor="grey", linewidth=0.2, 
                            ax=axes[2], legend=False, norm=norm_diff, missing_kwds={"color": "white"}, rasterized=True)

background = gdf_valid.to_crs("EPSG:4326") 
region_boundary = model_region_gdf.to_crs("EPSG:4326")
minx, miny, maxx, maxy = region_boundary.bounds.minx.item(), region_boundary.bounds.miny.item(), region_boundary.bounds.maxx.item(), region_boundary.bounds.maxy.item()

for i, ax in enumerate(axes):
    # Add model region
    region_boundary.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3)

    # # Add background and set extent (based on actual lat/lon coordinates)
    background.plot(ax=ax, color='#E0E0E0', zorder=0)

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
        gl.left_labels = False  # disable y-axis labels
    ax.text(0, 1.02, subplot_labels[i], transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='bottom', ha='left')
    ax.set_extent([minx, maxx, miny, maxy], crs=ccrs.PlateCarree())

# Titles
axes[0].set_title("Factual", fontsize=10)
axes[1].set_title("Counterfactual", fontsize=10)
axes[2].set_title("Factual - Counterfactual", fontsize=10)

# ==== Colorbar for Damage ====
sm1 = ScalarMappable(cmap="Reds", norm=plt.Normalize(vmin=0, vmax=100))
sm1.set_array([])
cbar = fig.colorbar(sm1, ax=axes[0:2], orientation="vertical", shrink=0.4, pad=0.01)
cbar.set_label("Aggregated Relative Damage [%]", fontsize=9)
cbar.ax.tick_params(labelsize=8)

sm2 = ScalarMappable(cmap=white_red, norm=norm_diff)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, ax=axes[2], orientation="vertical", shrink=0.4, pad=0.01)
cbar2.set_label("Attributable Relative Damage [%]", fontsize=9)
cbar2.ax.tick_params(labelsize=8)

fig.savefig("../figures/f_and_cf_aggregated_rel_damage.png", bbox_inches='tight', dpi=300)
fig.savefig("../figures/f_and_cf_aggregated_rel_damage.pdf", bbox_inches='tight', dpi=300)

# %%

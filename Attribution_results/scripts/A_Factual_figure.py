#%% Use compass-wflow pixi environment
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
from matplotlib.colors import PowerNorm
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import ScalarMappable
import matplotlib.patheffects as path_effects
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
from matplotlib import gridspec
import cartopy.crs as ccrs
from shapely.geometry import box
from rasterio.features import shapes
from shapely.geometry import shape
import matplotlib.dates as mdates

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


# %%
print("Plotting factual hmax masked and timeseries for specific gauges/stations")
gauges_list   = [1,2]
stations_list = [5, 40]

fig = plt.figure(figsize=(12, 5), dpi=300, constrained_layout=True)

projection = mod.crs.to_epsg()
colors = plt.get_cmap('tab10').colors  # Tuple of 10 colors

gs = gridspec.GridSpec(nrows=3, ncols=2, figure=fig,
                       width_ratios=[2, 1],  # Left col 3x wider than right col
                       height_ratios=[1, 1, 1],  # Three rows on right all equal height
                       wspace=0.05, hspace=0.03)

ax_map = fig.add_subplot(gs[:, 0], projection=ccrs.PlateCarree())

# Plot hmax masked
utm_crs = ccrs.UTM(zone=36, southern_hemisphere=True)
im = mod.results['hmax_masked'].plot.pcolormesh(
    ax=ax_map, cmap="viridis", vmin=0, vmax=3.5, add_colorbar=False, transform=utm_crs, rasterized=True)

ax_map.set_title("Factual Max Flood Depth", fontsize=10)

gdf_valid.plot(ax=ax_map, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)

# Add model region
model_region_gdf.boundary.plot(ax=ax_map, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())

# Set extent (based on actual lat/lon coordinates)
minx, miny, maxx, maxy = model_region_gdf.bounds.minx.item(), model_region_gdf.bounds.miny.item(), model_region_gdf.bounds.maxx.item(), model_region_gdf.bounds.maxy.item()
ax_map.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

# Add gridlines and format tick labels
gl = ax_map.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.xlocator = mticker.FixedLocator(np.arange(minx, maxx + 0.1, 0.2))
gl.ylocator = mticker.FixedLocator(np.arange(miny, maxy + 0.1, 0.2))
gl.xformatter = mticker.FuncFormatter(lon_formatter)
gl.yformatter = mticker.FuncFormatter(lat_formatter)
gl.right_labels = False
gl.top_labels = False
gl.xlabel_style = {'size': 9}
gl.ylabel_style = {'size': 9}

# ==== Plot city and river names ====
# Plot Beira location
ax_map.plot(34.848, -19.832, marker='o', color='black', markersize=4, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
text = ax_map.text(34.852, -19.81, "Beira", transform=ccrs.PlateCarree(),
                    fontsize=8, color='black', zorder=5)
text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])

# Buzi River marker and label
ax_map.plot(34.43, -19.89, marker='o', color='black', markersize=4, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
text2 = ax_map.text(34.44, -19.87, "Buzi River", transform=ccrs.PlateCarree(),
                fontsize=8, color='black', zorder=5)
text2.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
# Pungwe River marker and label
ax_map.plot(34.543, -19.545, marker='o', color='black', markersize=4, markeredgecolor='white', transform=ccrs.PlateCarree(), zorder=5)
text3 = ax_map.text(34.554, -19.52, "Pungwe River", transform=ccrs.PlateCarree(),
                fontsize=8, color='black', zorder=5)
text3.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])

# Colorbar for map
cbar = fig.colorbar(im, ax=ax_map, orientation="vertical", shrink=0.6, pad=0.01)
cbar.set_label('Maximum flood depth (m)', labelpad=10, fontsize=9)
cbar.ax.tick_params(labelsize=9)

# RIGHT: three stacked time series subplots
ax0 = fig.add_subplot(gs[0, 1])
ax1 = fig.add_subplot(gs[1, 1], sharex=ax0)
ax2 = fig.add_subplot(gs[2, 1], sharex=ax0)

# Plotting timeseries
ax0.plot(ds_his.time, ds_his['point_zs'].isel(stations=stations_list[0]), color=colors[4], label=f'S{stations_list[0]}')
ax0.plot(mod.forcing['bzs'].time, mod.forcing['bzs'].sel(index=stations_list[1]), color=colors[1], label=f'S{stations_list[1]}')
ax0.set_ylabel("Water level height \n [m]", fontsize=9)
ax0.legend(fontsize=8, loc="upper right", markerscale=0.8)
ax0.grid(True, linestyle='--', alpha=0.6)
ax0.tick_params(labelsize=9)
ax0.set_xlim(ds_his.time.min(), ds_his.time.max())
ax0.set_title("Factual time series", fontsize=10)
ax0.tick_params(bottom=False, labelbottom=False)

ax1.plot(mod.forcing['dis'].time, mod.forcing['dis'].sel(index=gauges_list[0]), color=colors[2], label=f'G{gauges_list[0]}')
ax1.plot(mod.forcing['dis'].time, mod.forcing['dis'].sel(index=gauges_list[1]), color=colors[3], label=f'G{gauges_list[1]}')
ax1.set_ylabel("Discharge \n [m³/s]", fontsize=9)
ax1.legend(fontsize=8, loc="upper right", markerscale=0.8)
ax1.grid(True, linestyle='--', alpha=0.6)
ax1.tick_params(labelsize=9)
ax1.set_xlim(ds_his.time.min(), ds_his.time.max())
ax1.tick_params(bottom=False, labelbottom=False)

ax2.step(mod.forcing['precip_2d'].time, mod.forcing['precip_2d'].sum(dim=["x", "y"]), where='post', color=colors[0])
ax2.set_ylabel("Accum. precipitation \n [mm/h]", fontsize=10)
ax2.set_xlabel("Day in March 2019")
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
ax2.grid(True, linestyle='--', alpha=0.6)
ax2.tick_params(labelsize=9)
ax2.set_xlim(ds_his.time.min(), ds_his.time.max())

# Plot station and Gauge points on the flood map and annotate
ax_map.scatter(ds_his['point_zs'].isel(stations=stations_list[0]).point_x.values, ds_his['point_zs'].isel(stations=stations_list[0]).point_y.values, color=colors[4], s=30, edgecolor='k', zorder=5, transform=ccrs.epsg(projection), label=f'Station {stations_list[0]}')
ax_map.scatter(mod.forcing['bzs'].sel(index=stations_list[1]).geometry.item().x, mod.forcing['bzs'].sel(index=stations_list[1]).geometry.item().y, color=colors[1], s=30, edgecolor='k', zorder=5, transform=ccrs.epsg(projection), label=f'Station {stations_list[1]}')

ax_map.scatter(mod.forcing['dis'].sel(index=gauges_list[0]).geometry.item().x, mod.forcing['dis'].sel(index=gauges_list[0]).geometry.item().y, color=colors[2], s=30, edgecolor='k', zorder=5, transform=ccrs.epsg(projection), label=f'Gauge {gauges_list[0]}')
ax_map.scatter(mod.forcing['dis'].sel(index=gauges_list[1]).geometry.item().x, mod.forcing['dis'].sel(index=gauges_list[1]).geometry.item().y, color=colors[3], s=30, edgecolor='k', zorder=5, transform=ccrs.epsg(projection), label=f'Gauge {gauges_list[1]}')

offset = 2000  # Adjust offset depending on map scale (in map projection units)
ax_map.annotate(f'S{stations_list[0]}', xy=(ds_his['point_zs'].isel(stations=stations_list[0]).point_x.values, ds_his['point_zs'].isel(stations=stations_list[0]).point_y.values),
                xytext=(ds_his['point_zs'].isel(stations=stations_list[0]).point_x.values + 2000, ds_his['point_zs'].isel(stations=stations_list[0]).point_y.values - 4000),
                transform=ccrs.epsg(projection),
                fontsize=8, color=colors[4], fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=0.5, alpha=0.8)
)
ax_map.annotate(f'S{stations_list[1]}', xy=(mod.forcing['bzs'].sel(index=stations_list[1]).geometry.item().x, mod.forcing['bzs'].sel(index=stations_list[1]).geometry.item().y),
                xytext=(mod.forcing['bzs'].sel(index=stations_list[1]).geometry.item().x + 2000, mod.forcing['bzs'].sel(index=stations_list[1]).geometry.item().y - 4000),
                transform=ccrs.epsg(projection),
                fontsize=8, color=colors[1], fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=0.5, alpha=0.8)
)
ax_map.annotate(f'G{gauges_list[0]}', xy=(mod.forcing['dis'].sel(index=gauges_list[0]).geometry.item().x, mod.forcing['dis'].sel(index=gauges_list[0]).geometry.item().y),
                xytext=(mod.forcing['dis'].sel(index=gauges_list[0]).geometry.item().x - 6500, mod.forcing['dis'].sel(index=gauges_list[0]).geometry.item().y - 4000),
                transform=ccrs.epsg(projection),
                fontsize=8, color=colors[2], fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=0.5, alpha=0.8)
)
ax_map.annotate(f'G{gauges_list[1]}', xy=(mod.forcing['dis'].sel(index=gauges_list[1]).geometry.item().x, mod.forcing['dis'].sel(index=gauges_list[1]).geometry.item().y),
                xytext=(mod.forcing['dis'].sel(index=gauges_list[1]).geometry.item().x + 2500, mod.forcing['dis'].sel(index=gauges_list[1]).geometry.item().y),
                transform=ccrs.epsg(projection),
                fontsize=8, color=colors[3], fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=0.5, alpha=0.8)
)
props = dict( facecolor='w', lw=0, alpha=0.8)
ax_map.text(0.03, 0.95, '(a)', transform=ax_map.transAxes, bbox=props, fontsize=10, fontweight='bold')
ax0.text(0.02, 0.85, '(b)', transform=ax0.transAxes, bbox=props, fontsize=10, fontweight='bold')
ax1.text(0.02, 0.85, '(c)', transform=ax1.transAxes, bbox=props, fontsize=10, fontweight='bold')
ax2.text(0.02, 0.85, '(d)', transform=ax2.transAxes, bbox=props, fontsize=10, fontweight='bold')


fig.savefig("../figures/fS13.png", bbox_inches='tight', dpi=300)
fig.savefig("../figures/fS13.pdf", bbox_inches='tight', dpi=300)
plt.show()


# %%
# Some numbers
max_buzi = mod.forcing['dis'].sel(index=gauges_list[0]).max()
max_pungwe = mod.forcing['dis'].sel(index=gauges_list[1]).max()

max_buzi_1000 = np.round(max_buzi.values, -3)
max_pungwe_1000 = np.round(max_pungwe.values, -3)

print(f"Max discharge Buzi gauge: {max_buzi_1000} m3/s (rounded to 1000)")
print(f"Max discharge Pungwe gauge: {max_pungwe_1000} m3/s (rounded to 1000)")


#%%
# Coastal water level
max_S5 = ds_his['point_zs'].isel(stations=stations_list[0]).max().values
max_tide_S5 = ds_his['point_zs'].sel(time=slice("2019-03-19", "2019-03-25"), stations=stations_list[0]).max().values
max_storm_surge_S5 = max_S5 - max_tide_S5

print(f"Max coastal water level at S5: {max_S5:.1f} m")
print(f"Max storm surge at S5: {max_storm_surge_S5:.1f} m")
# %%

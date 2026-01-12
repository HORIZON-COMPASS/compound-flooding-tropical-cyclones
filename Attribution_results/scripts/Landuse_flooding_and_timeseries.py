#%% Use compass-wflow pixi environment
print("Loading packages")
import os
from os.path import join
import xarray as xr
import numpy as np
import geopandas as gpd
from hydromt_sfincs import SfincsModel
from hydromt import DataCatalog
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.dates as mdates
import rioxarray as rxr  # Required for reading TIFF files
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
from matplotlib.colors import TwoSlopeNorm


#%%
print("Loading model paths and data catalog")
# define file paths to directory when different model outputs are stored
models_dir = join("../data/sfincs/") # change this to your own path

datacat = ['../../Workflows/03_data_catalogs/datacatalog_general.yml']
data_catalog = DataCatalog(data_libs = datacat)

# Define different model directories
# Factual model from Doris' submitted paper with ESA worldcover in wflow domain and vito in sfincs domain
doris_factual_model_dir        = join(models_dir, "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_doris")
# New factual model with updated landuse for 2020
factual_model_dir            = join(models_dir, "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_lisboa_2020")
# Counterfactual climate model with 8% of precipitation removed acc. to Clausius Clapeyron at 1.1 degree of warming, 
# 14 cm of sea level rise (SLR) removed and with 10% reduced maximum tropical cyclone wind speeds
CF_climate_model_dir         = join(models_dir, "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10_lisboa_2020")
# Counterfactual landuse model
CF_landuse_model_dir         = join(models_dir, "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_lisboa_2000")
# Counterfactual climate and landuse model
CF_climate_landuse_model_dir = join(models_dir, "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10_lisboa_2000")

#%%
# Load in models and model region
print("Loading model data and model region")
# Factual models submitted for Doris paper and run for Poppy's paper
mod_F_doris = SfincsModel(doris_factual_model_dir, data_libs=datacat, mode="r")
ds_his_F_doris= xr.open_dataset(join(doris_factual_model_dir,"sfincs_his.nc"), engine='netcdf4')

mod_F  = SfincsModel(factual_model_dir, data_libs=datacat, mode="r")
ds_his_F = xr.open_dataset(join(factual_model_dir,"sfincs_his.nc"), engine='netcdf4')

# Counterfactual climate model
mod_CF_clim  = SfincsModel(CF_climate_model_dir, data_libs=datacat, mode="r")
ds_his_CF_clim = xr.open_dataset(join(CF_climate_model_dir,"sfincs_his.nc"), engine='netcdf4')

# Counterfactual landuse model
mod_CF_landuse  = SfincsModel(CF_landuse_model_dir, data_libs=datacat, mode="r")
ds_his_CF_landuse = xr.open_dataset(join(CF_landuse_model_dir,"sfincs_his.nc"), engine='netcdf4')

# Counterfactual climate and landuse model
mod_CF_clim_landuse  = SfincsModel(CF_climate_landuse_model_dir, data_libs=datacat, mode="r")
ds_his_CF_clim_landuse = xr.open_dataset(join(CF_climate_landuse_model_dir,"sfincs_his.nc"), engine='netcdf4')

# Model region
model_region_gdf = gpd.read_file(join("..", "data", "sfincs", "gis", "region.geojson")).to_crs("EPSG:4326") 

#%%
print("Loading flood maps")
# Load flood map which is already downscaled, masked for permanent water and represents a cells as flooded from 0.05 m or more
model_dirs = {
    "F_doris": doris_factual_model_dir,
    "F": factual_model_dir,
    "CF_clim": CF_climate_model_dir,
    "CF_landuse": CF_landuse_model_dir,
    "CF_clim_landuse": CF_climate_landuse_model_dir,
}

hmax = {}

for key, d in model_dirs.items():
    da = rxr.open_rasterio(
        Path(d) / "floodmap.tif",
        masked=True,
        chunks="auto"  # enables lazy evaluation
    )

    if "band" in da.dims:
        da = da.squeeze("band", drop=True)

    hmax[key] = da

#%%
# Load own background shape for plotting that aligns with permanent water mask: permanent water where water occurence > 5%
gdf_valid = gpd.read_file(join("..", "data", "valid_landmask.gpkg"), layer="valid_mask").to_crs(model_region_gdf.crs)

# settings for plotting
projection = 32736
utm_crs = ccrs.UTM(zone=36, southern_hemisphere=True)


#%%
# Plot max flood depth maps for comparison
fig, axes = plt.subplots(1,3, figsize=(15, 8), dpi=300, subplot_kw={"projection": utm_crs},
                       constrained_layout=True)

# Plot hmax masked
im1 = hmax['F_doris'].plot.pcolormesh(ax=axes[0], cmap="viridis", vmin=0, vmax=3.5, 
                               add_colorbar=False, transform=utm_crs, rasterized=True)
im2 = hmax['F'].plot.pcolormesh(ax=axes[1], cmap="viridis", vmin=0, vmax=3.5, 
                               add_colorbar=True, transform=utm_crs, rasterized=True)
hmax_diff = hmax['F'] - hmax['F_doris']
im3 = hmax_diff.plot.pcolormesh(ax=axes[2], cmap="bwr", vmin=-3.5, vmax=3.5, 
                               add_colorbar=True, transform=utm_crs, rasterized=True)

axes[0].set_title("Factual (Doris) Max Flood Depth", fontsize=10)
axes[1].set_title("Factual (Poppy) Max Flood Depth", fontsize=10)
axes[2].set_title("Difference in Max Flood Depth", fontsize=10)

minx, miny, maxx, maxy = model_region_gdf.bounds.minx.item(), model_region_gdf.bounds.miny.item(), model_region_gdf.bounds.maxx.item(), model_region_gdf.bounds.maxy.item()
for ax in axes:
    # Add gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.2, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    # Plot valid landmask
    gdf_valid.plot(ax=ax, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)
    # Add model region
    model_region_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())
    # Set extent (based on actual lat/lon coordinates)
    ax.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

#%%
print("Plotting factual hmax masked and timeseries for specific gauges/stations")
gauges_list   = [1,2]
stations_list = [5, 40]

fig, axes = plt.subplots(3, 1, figsize=(10, 10), dpi=300, constrained_layout=True)

projection = 32736
colors = plt.get_cmap('tab10').colors  # Tuple of 10 colors

# Plotting timeseries
axes[0].plot(ds_his_F_doris.time, ds_his_F_doris['point_zs'].isel(stations=stations_list[0]), color=colors[4], label=f'S{stations_list[0]}_F_doris')
axes[0].plot(mod_F_doris.forcing['bzs'].time, mod_F_doris.forcing['bzs'].sel(index=stations_list[1]), color=colors[1], label=f'S{stations_list[1]}_F_doris')
axes[0].plot(ds_his_F.time, ds_his_F['point_zs'].isel(stations=stations_list[0]), color=colors[4], label=f'S{stations_list[0]}_F', linestyle='--')
axes[0].plot(mod_F.forcing['bzs'].time, mod_F.forcing['bzs'].sel(index=stations_list[1]), color=colors[1], label=f'S{stations_list[1]}_F', linestyle='--')
axes[0].set_ylabel("Water level height \n [m]", fontsize=9)
axes[0].legend(fontsize=8, loc="upper right", markerscale=0.8)
axes[0].grid(True, linestyle='--', alpha=0.6)
axes[0].tick_params(labelsize=9)
axes[0].set_xlim(ds_his_F_doris.time.min(), ds_his_F_doris.time.max())
axes[0].set_title("Factual time series", fontsize=10)
axes[0].tick_params(bottom=False, labelbottom=False)

axes[1].plot(mod_F_doris.forcing['dis'].time, mod_F_doris.forcing['dis'].sel(index=gauges_list[0]), color=colors[2], label=f'G{gauges_list[0]}_F_doris')
axes[1].plot(mod_F_doris.forcing['dis'].time, mod_F_doris.forcing['dis'].sel(index=gauges_list[1]), color=colors[3], label=f'G{gauges_list[1]}_F_doris')
axes[1].plot(mod_F.forcing['dis'].time, mod_F.forcing['dis'].sel(index=gauges_list[0]), color=colors[2], label=f'G{gauges_list[0]}_F', linestyle='--')
axes[1].plot(mod_F.forcing['dis'].time, mod_F.forcing['dis'].sel(index=gauges_list[1]), color=colors[3], label=f'G{gauges_list[1]}_F', linestyle='--')
axes[1].set_ylabel("Discharge \n [m³/s]", fontsize=9)
axes[1].legend(fontsize=8, loc="upper right", markerscale=0.8)
axes[1].grid(True, linestyle='--', alpha=0.6)
axes[1].tick_params(labelsize=9)
axes[1].set_xlim(ds_his_F_doris.time.min(), ds_his_F_doris.time.max())
axes[1].tick_params(bottom=False, labelbottom=False)

axes[2].step(mod_F_doris.forcing['precip_2d'].time, mod_F_doris.forcing['precip_2d'].sum(dim=["x", "y"]), where='post', color=colors[0], label='Precipitation_F_doris')
axes[2].step(mod_F.forcing['precip_2d'].time, mod_F.forcing['precip_2d'].sum(dim=["x", "y"]), where='post', color=colors[1], label='Precipitation_F', linestyle='--')
axes[2].set_ylabel("Accum. precipitation \n [mm/h]", fontsize=10)
axes[2].set_xlabel("Day in March 2019")
axes[2].xaxis.set_major_formatter(mdates.DateFormatter('%d'))
axes[2].grid(True, linestyle='--', alpha=0.6)
axes[2].legend(fontsize=8, loc="upper right", markerscale=0.8)
axes[2].tick_params(labelsize=9)
axes[2].set_xlim(ds_his_F_doris.time.min(), ds_his_F_doris.time.max())

# %%
print("Plotting factual and counterfactual timeseries for specific gauges/stations")
gauges_list   = [1,2]
stations_list = [5, 40]

fig, axes = plt.subplots(3, 1, figsize=(15, 10), dpi=300, constrained_layout=True)

projection = 32736
colors = plt.get_cmap('tab10').colors  # Tuple of 10 colors

# Plotting timeseries
axes[0].plot(ds_his_F.time, ds_his_F['point_zs'].isel(stations=stations_list[0]), color=colors[4], label=f'S{stations_list[0]}_F')
axes[0].plot(mod_F.forcing['bzs'].time, mod_F.forcing['bzs'].sel(index=stations_list[1]), color=colors[1], label=f'S{stations_list[1]}_F')
axes[0].plot(ds_his_CF_clim.time, ds_his_CF_clim['point_zs'].isel(stations=stations_list[0]), color=colors[4], label=f'S{stations_list[0]}_CF_clim', linestyle='--')
axes[0].plot(mod_CF_clim.forcing['bzs'].time, mod_CF_clim.forcing['bzs'].sel(index=stations_list[1]), color=colors[1], label=f'S{stations_list[1]}_CF_clim', linestyle='--')
axes[0].set_ylabel("Water level height \n [m]", fontsize=9)
axes[0].legend(fontsize=8, loc="upper right", markerscale=0.8)
axes[0].grid(True, linestyle='--', alpha=0.6)
axes[0].tick_params(labelsize=9)
axes[0].set_xlim(ds_his_F.time.min(), ds_his_F.time.max())
axes[0].set_title("Factual time series", fontsize=10)
axes[0].tick_params(bottom=False, labelbottom=False)

axes[1].plot(mod_F.forcing['dis'].time, mod_F.forcing['dis'].sel(index=gauges_list[0]), color=colors[2], label=f'G{gauges_list[0]}_F')
axes[1].plot(mod_F.forcing['dis'].time, mod_F.forcing['dis'].sel(index=gauges_list[1]), color=colors[3], label=f'G{gauges_list[1]}_F')
axes[1].plot(mod_CF_clim.forcing['dis'].time, mod_CF_clim.forcing['dis'].sel(index=gauges_list[0]), color=colors[2], label=f'G{gauges_list[0]}_CF_clim', linestyle='--')
axes[1].plot(mod_CF_clim.forcing['dis'].time, mod_CF_clim.forcing['dis'].sel(index=gauges_list[1]), color=colors[3], label=f'G{gauges_list[1]}_CF_clim', linestyle='--')
axes[1].plot(mod_CF_landuse.forcing['dis'].time, mod_CF_landuse.forcing['dis'].sel(index=gauges_list[0]), color=colors[2], label=f'G{gauges_list[0]}_CF_lu', linestyle=':')
axes[1].plot(mod_CF_landuse.forcing['dis'].time, mod_CF_landuse.forcing['dis'].sel(index=gauges_list[1]), color=colors[3], label=f'G{gauges_list[1]}_CF_lu', linestyle=':')
axes[1].plot(mod_CF_clim_landuse.forcing['dis'].time, mod_CF_clim_landuse.forcing['dis'].sel(index=gauges_list[0]), color=colors[2], label=f'G{gauges_list[0]}_CF_clim_lu', linestyle='-.')
axes[1].plot(mod_CF_clim_landuse.forcing['dis'].time, mod_CF_clim_landuse.forcing['dis'].sel(index=gauges_list[1]), color=colors[3], label=f'G{gauges_list[1]}_CF_clim_lu', linestyle='-.')
axes[1].set_ylabel("Discharge \n [m³/s]", fontsize=9)
axes[1].legend(fontsize=8, loc="upper right", markerscale=0.8)
axes[1].grid(True, linestyle='--', alpha=0.6)
axes[1].tick_params(labelsize=9)
axes[1].set_xlim(ds_his_F.time.min(), ds_his_F.time.max())
axes[1].tick_params(bottom=False, labelbottom=False)

axes[2].step(mod_F.forcing['precip_2d'].time, mod_F.forcing['precip_2d'].sum(dim=["x", "y"]), where='post', color=colors[0], label='Precipitation_F')
axes[2].step(mod_CF_clim.forcing['precip_2d'].time, mod_CF_clim.forcing['precip_2d'].sum(dim=["x", "y"]), where='post', color=colors[1], label='Precipitation_CF', linestyle='--')
axes[2].set_ylabel("Accum. precipitation \n [mm/h]", fontsize=10)
axes[2].set_xlabel("Day in March 2019")
axes[2].xaxis.set_major_formatter(mdates.DateFormatter('%d'))
axes[2].grid(True, linestyle='--', alpha=0.6)
axes[2].legend(fontsize=8, loc="upper right", markerscale=0.8)
axes[2].tick_params(labelsize=9)
axes[2].set_xlim(ds_his_F.time.min(), ds_his_F.time.max())

# %%
print("Plotting counterfactual hmax masked maps")
# Plot max flood depth maps
fig, axes = plt.subplots(1,3, figsize=(15, 8), dpi=300, subplot_kw={"projection": utm_crs},
                       constrained_layout=True)

# Plot hmax masked
hmax_diff_CF_clim = hmax['F'] - hmax['CF_clim']
hmax_diff_CF_clim.load()
# norm_CF_clim = TwoSlopeNorm(vmin=0, vmax=hmax_diff_CF_clim.quantile(0.99).values)
im1 = hmax_diff_CF_clim.plot.pcolormesh(ax=axes[0], cmap="Reds", vmin=0, vmax=hmax_diff_CF_clim.quantile(0.99).values,
                               add_colorbar=True, transform=utm_crs, rasterized=True)

hmax_diff_CF_landuse = hmax['F'] - hmax['CF_landuse']
hmax_diff_CF_landuse.load()
norm_CF_landuse = TwoSlopeNorm(vmin=hmax_diff_CF_landuse.quantile(0.01).values, vcenter=0.0, vmax=hmax_diff_CF_landuse.quantile(0.99).values)
im2 = hmax_diff_CF_landuse.plot.pcolormesh(ax=axes[1], cmap="bwr", norm=norm_CF_landuse, 
                               add_colorbar=True, transform=utm_crs, rasterized=True)

hmax_diff_CF_clim_landuse = hmax['F'] - hmax['CF_clim_landuse']
hmax_diff_CF_clim_landuse.load()
norm_CF_clim_landuse = TwoSlopeNorm(vmin=hmax_diff_CF_clim_landuse.quantile(0.01).values, vcenter=0.0, vmax=hmax_diff_CF_clim_landuse.quantile(0.99).values)
im3 = hmax_diff_CF_clim_landuse.plot.pcolormesh(ax=axes[2], cmap="bwr", norm=norm_CF_clim_landuse, 
                               add_colorbar=True, transform=utm_crs, rasterized=True)

axes[0].set_title("F - CF climate", fontsize=10)
axes[1].set_title("F - CF landuse", fontsize=10)
axes[2].set_title("F - CF climate & landuse", fontsize=10)
fig.suptitle("Difference in max flood depth", fontsize=11)

minx, miny, maxx, maxy = model_region_gdf.bounds.minx.item(), model_region_gdf.bounds.miny.item(), model_region_gdf.bounds.maxx.item(), model_region_gdf.bounds.maxy.item()
for ax in axes:
    # Add gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.2, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    # Plot valid landmask
    gdf_valid.plot(ax=ax, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)
    # Add model region
    model_region_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())
    # Set extent (based on actual lat/lon coordinates)
    ax.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

# %%

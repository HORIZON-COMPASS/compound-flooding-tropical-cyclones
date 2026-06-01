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


#%% This script - shows maps of the scenarios for the whole model region. 


print("Loading model paths and data catalog")
# define file paths to directory when different model outputs are stored
models_dir = join("/Data/Scenario_floodmaps") # change this to your own path

#datacat = ['../../Workflows/03_data_catalogs/datacatalog_general.yml']
#data_catalog = DataCatalog(data_libs = datacat)

# Define different model directories
# Factual model from Doris' submitted paper with ESA worldcover in wflow domain and vito in sfincs domain
# doris_factual_model_dir        = join(models_dir, "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_doris")
# New factual model with updated landuse for 2020
factual_model_dir            = join(models_dir, "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_lisboa_2020")
# Counterfactual climate model with 8% of precipitation removed acc. to Clausius Clapeyron at 1.1 degree of warming, 
# 14 cm of sea level rise (SLR) removed and with 10% reduced maximum tropical cyclone wind speeds
CF_climate_model_dir         = join(models_dir, "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.1_era5_hourly_spw_IBTrACS_CF-5_lisboa_2020")
# Counterfactual landuse model
CF_landuse_model_dir         = join(models_dir, "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_lisboa_2000")
# Counterfactual climate and landuse model
CF_climate_landuse_model_dir = join(models_dir, "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.1_era5_hourly_spw_IBTrACS_CF-5_lisboa_2000")

#%%
# Load in models and model region
print("Loading model data and model region")
# Factual models submitted for Doris paper and run for Poppy's paper
#mod_F_doris = SfincsModel(doris_factual_model_dir, data_libs=datacat, mode="r")
#ds_his_F_doris= xr.open_dataset(join(doris_factual_model_dir,"sfincs_his.nc"), engine='netcdf4')

#mod_F  = SfincsModel(factual_model_dir, data_libs=datacat, mode="r")
ds_his_F = xr.open_dataset(join(factual_model_dir,"sfincs_his.nc"), engine='netcdf4')

# Counterfactual climate model
#mod_CF_clim  = SfincsModel(CF_climate_model_dir, data_libs=datacat, mode="r")
ds_his_CF_clim = xr.open_dataset(join(CF_climate_model_dir,"sfincs_his.nc"), engine='netcdf4')

# Counterfactual landuse model
#mod_CF_landuse  = SfincsModel(CF_landuse_model_dir, data_libs=datacat, mode="r")
ds_his_CF_landuse = xr.open_dataset(join(CF_landuse_model_dir,"sfincs_his.nc"), engine='netcdf4')

# Counterfactual climate and landuse model
#mod_CF_clim_landuse  = SfincsModel(CF_climate_landuse_model_dir, data_libs=datacat, mode="r")
ds_his_CF_clim_landuse = xr.open_dataset(join(CF_climate_landuse_model_dir,"sfincs_his.nc"), engine='netcdf4')

# Model region
model_region_gdf = gpd.read_file("Data/region.geojson").to_crs("EPSG:4326") 

#%%
print("Loading flood maps")
# Load flood map which is already downscaled, masked for permanent water and represents a cells as flooded from 0.05 m or more
model_dirs = {
    #"F_doris": doris_factual_model_dir,
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
gdf_valid = gpd.read_file(
    "/Data/gis_files/valid_landmask.gpkg",
    layer="valid_mask",
).to_crs(model_region_gdf.crs)

# settings for plotting
projection = 32736
utm_crs = ccrs.UTM(zone=36, southern_hemisphere=True)



#%%
# Plot max flood depth maps for comparison
fig, axes = plt.subplots(1,3, figsize=(15, 8), dpi=300, subplot_kw={"projection": utm_crs},
                       constrained_layout=True)

# Plot hmax masked
im1 = hmax['F'].plot.pcolormesh(ax=axes[0], cmap="viridis", vmin=0, vmax=3.5, 
                               add_colorbar=False, transform=utm_crs, rasterized=True)
im2 = hmax['CF_landuse'].plot.pcolormesh(ax=axes[1], cmap="viridis", vmin=0, vmax=3.5, 
                               add_colorbar=True, transform=utm_crs, rasterized=True)
hmax_diff = hmax['F'] - hmax['CF_landuse']
im3 = hmax_diff.plot.pcolormesh(ax=axes[2], cmap="bwr", vmin=-2, vmax=2, 
                               add_colorbar=True, transform=utm_crs, rasterized=True)

axes[0].set_title("Factual Max Flood Depth (2020)", fontsize=10)
axes[1].set_title("Counterfactual Max Flood Depth (2000)", fontsize=10)
axes[2].set_title("Difference in Max Flood Depth (2020-2000)", fontsize=10)

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
# %% Comparing ( Factual - Counterfactual) 
print("Plotting counterfactual hmax masked maps")
# Plot max flood depth maps
fig, axes = plt.subplots(1,3, figsize=(15, 8), dpi=300, subplot_kw={"projection": utm_crs},
                       constrained_layout=True)

# Plot hmax masked - Factual 2020 (observed climate change) - COunterfactual 2000 (removed climate change)
hmax_diff_CF_clim = hmax['F'] - hmax['CF_clim']
hmax_diff_CF_clim.load()
#norm_CF_clim = TwoSlopeNorm(vmin=0, vmax=hmax_diff_CF_clim.quantile(0.99).values, vcenter=0.0,)
im1 = hmax_diff_CF_clim.plot.pcolormesh(ax=axes[0], cmap="Reds",vmin=0, vmax=hmax_diff_CF_clim.quantile(0.99).values,
                               add_colorbar=True, transform=utm_crs, rasterized=True)

# Factual 2020 land-use (observed climate)  - Counterfactual 2000 land use (+climate change removed)
hmax_diff_CF_landuse = hmax['F'] - hmax['CF_landuse']
hmax_diff_CF_landuse.load()
norm_CF_landuse = TwoSlopeNorm(vmin=hmax_diff_CF_landuse.quantile(0.01).values, vcenter=0.0, vmax=hmax_diff_CF_landuse.quantile(0.99).values)
im2 = hmax_diff_CF_landuse.plot.pcolormesh(ax=axes[1], cmap="bwr", norm=norm_CF_landuse, 
                               add_colorbar=True, transform=utm_crs, rasterized=True)

# Factual 
hmax_diff_CF_clim_landuse = hmax['F'] - hmax['CF_clim_landuse']
hmax_diff_CF_clim_landuse.load()
norm_CF_clim_landuse = TwoSlopeNorm(vmin=hmax_diff_CF_clim_landuse.quantile(0.01).values, vcenter=0.0, vmax=hmax_diff_CF_clim_landuse.quantile(0.99).values)
im3 = hmax_diff_CF_clim_landuse.plot.pcolormesh(ax=axes[2], cmap="bwr", norm=norm_CF_clim_landuse, 
                               add_colorbar=True, transform=utm_crs, rasterized=True)

axes[0].set_title("CLimate Change", fontsize=10)
axes[1].set_title("LULC Change", fontsize=10)
axes[2].set_title("Climate & LULC Change", fontsize=10)
fig.suptitle("Difference in max flood depth between scenarios (m)", fontsize=14)

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
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs

# --- compute diffs ---
hmax_diff_CF_clim = (hmax["F"] - hmax["CF_clim"]).load()
hmax_diff_CF_landuse = (hmax["F"] - hmax["CF_landuse"]).load()
hmax_diff_CF_both = (hmax["F"] - hmax["CF_clim_landuse"]).load()

data_list = [hmax_diff_CF_clim, hmax_diff_CF_landuse, hmax_diff_CF_both]
titles = ["Climate change", "LULC Change", "Climate + LULC"]
panel_labels = ["(a)", "(b)", "(c)"]

# --- shared scale (pick one) ---
# If you want the "paper" style like your screenshot (positive-only Reds):
vmin = -0.1
vmax = float(max([da.quantile(0.99).values for da in data_list]))  # shared vmax

cmap = plt.get_cmap("Reds").copy()
cmap.set_bad(alpha=0)  # make NaNs transparent (optional)
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

# --- bigger plots, smaller margins ---
fig, axes = plt.subplots(
    1, 3,
    figsize=(15.5, 6.2),   # wider figure = bigger panels
    dpi=300,
    subplot_kw={"projection": utm_crs}
)

# Reserve only a small strip on the right for the colorbar
fig.subplots_adjust(left=0.04, right=0.90, bottom=0.06, top=0.92, wspace=0.04)

minx = model_region_gdf.bounds.minx.item()
miny = model_region_gdf.bounds.miny.item()
maxx = model_region_gdf.bounds.maxx.item()
maxy = model_region_gdf.bounds.maxy.item()

panel_labels = ["(a)", "(b)", "(c)"]

for ax, lab in zip(axes, panel_labels):
    bbox = ax.get_position()  # axes position in figure coords

    fig.text(
        bbox.x0,          # left edge of panel
        bbox.y1 + 0.01,   # slightly above panel
        lab,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="bold"
    )

for i, (ax, da, title, lab) in enumerate(zip(axes, data_list, titles)):
    ax.set_facecolor("#fffcfc")

    # context/landmask
    gdf_valid.plot(ax=ax, color="#e0e0e0", transform=ccrs.PlateCarree(), zorder=0)

    # OPTIONAL: positive-only like your screenshot
    da_plot = da.where(da > 0)

    # main layer (no per-axis colorbar)
    da_plot.plot.pcolormesh(
        ax=ax, cmap=cmap, norm=norm,
        add_colorbar=False,
        transform=utm_crs,
        rasterized=True,
        zorder=2
    )

    # domain outline
    model_region_gdf.boundary.plot(
        ax=ax, edgecolor="black", linewidth=0.6,
        transform=ccrs.PlateCarree(), zorder=3
    )

    # subtle gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.35, color="gray", alpha=0.35, linestyle="--")
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 8}
    gl.ylabel_style = {"size": 8}

    # Reduce label clutter: only left panel shows y labels; only middle shows x labels if you want
    if i != 0:
        gl.left_labels = False

    ax.set_title(title, fontsize=13)
    ax.text(0.02, 0.98, lab, transform=ax.transAxes, ha="left", va="top",
            fontsize=11, fontweight="bold")

    ax.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

# ---- ONE slim shared colorbar ----
# Make a dedicated tiny axis for the colorbar: [left, bottom, width, height]
cax = fig.add_axes([0.92, 0.16, 0.012, 0.68])  # thinner + tall, tweak to taste
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])

cbar = fig.colorbar(sm, cax=cax)
cbar.set_label("Difference in max flood depth between scenarios (m)", fontsize=10)
cbar.ax.tick_params(labelsize=8)

plt.show()

# %%
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import TwoSlopeNorm
import cartopy.crs as ccrs

data_list = [hmax_diff_CF_clim, hmax_diff_CF_landuse, hmax_diff_CF_both]
titles = ["Climate change", "LULC Change", "Climate + LULC"]
panel_labels = ["(a)", "(b)", "(c)"]

# --- shared diverging scale (robust + symmetric around 0) ---
q = 0.99
vabs = max(float(abs(da).quantile(q).values) for da in data_list)
vmin, vmax = -vabs, vabs

cmap = plt.get_cmap("RdBu_r").copy()   # diverging (blue=negative, red=positive)
cmap.set_bad(alpha=0)                  # NaNs transparent (optional)
norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

# --- bigger plots, smaller margins ---
fig, axes = plt.subplots(
    1, 3,
    figsize=(15.8, 6.2),
    dpi=300,
    subplot_kw={"projection": utm_crs}
)

# Reserve a small strip on the right for ONE shared colorbar
fig.subplots_adjust(left=0.04, right=0.90, bottom=0.06, top=0.92, wspace=0.04)

# extent (you already have these)
minx = model_region_gdf.bounds.minx.item()
miny = model_region_gdf.bounds.miny.item()
maxx = model_region_gdf.bounds.maxx.item()
maxy = model_region_gdf.bounds.maxy.item()

# --- plot panels ---
for i, (ax, da, title, lab) in enumerate(zip(axes, data_list, titles, panel_labels)):
    ax.set_facecolor("#eef3f6")  

    # context/landmask
    gdf_valid.plot(
        ax=ax,
        color="#d9d9d9",
        transform=ccrs.PlateCarree(),
        zorder=0
    )

    # main layer (NO per-axis colorbar)
    da.plot.pcolormesh(
        ax=ax,
        cmap=cmap,
        norm=norm,
        add_colorbar=False,
        transform=utm_crs,
        rasterized=True,
        zorder=2
    )

    # domain outline
    model_region_gdf.boundary.plot(
        ax=ax,
        edgecolor="black",
        linewidth=0.6,
        transform=ccrs.PlateCarree(),
        zorder=3
    )

    # subtle gridlines
    gl = ax.gridlines(
        draw_labels=True,
        linewidth=0.35,
        color="gray",
        alpha=0.35,
        linestyle="--"
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 8}
    gl.ylabel_style = {"size": 8}
    if i != 0:
        gl.left_labels = False  # reduce clutter

    # title (bigger)
    ax.set_title(title, fontsize=13, pad=8)

    ax.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

# --- panel labels ABOVE each plot (Doris style) ---
for ax, lab in zip(axes, panel_labels):
    bbox = ax.get_position()
    fig.text(
        bbox.x0,
        bbox.y1 + 0.01,
        lab,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="bold"
    )

# ---- ONE slim shared colorbar ----
cax = fig.add_axes([0.92, 0.16, 0.012, 0.68])  # [left, bottom, width, height]
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])

cbar = fig.colorbar(sm, cax=cax)
cbar.set_label("Difference in max flood depth between scenarios (m)", fontsize=10)
cbar.ax.tick_params(labelsize=8)

plt.show()

# %%

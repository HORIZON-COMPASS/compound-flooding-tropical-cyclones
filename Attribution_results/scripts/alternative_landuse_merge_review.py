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
from hydromt_wflow import WflowModel
import pandas as pd
import rioxarray as rxr
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colors import BoundaryNorm
from rasterio.enums import Resampling
from matplotlib.lines import Line2D 
from matplotlib.colors import TwoSlopeNorm


#%% This script - shows maps of the scenarios for the whole model region. 
print("Loading model paths and data catalog")
# define file paths to directory when different model outputs are stored
# models_dir = join("p:/11210471-001-compass/03_Runs/sofala/Idai/sfincs") # change this to your own path
runs_dir_sfincs = join("p:/11210471-001-compass/03_Runs/sofala/Idai/sfincs")
runs_dir_wflow = join("p:/11210471-001-compass/03_Runs/sofala/Idai/wflow")
models_dir_base = join("p:/11210471-001-compass/02_Models/sofala/Idai/")

curdir           = '../../Workflows/'
data_cats        = [join(curdir, "03_data_catalogs", "datacatalog_general.yml"), 
                    join(curdir, "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling.yml"), 
                    join(curdir, "03_data_catalogs", "datacatalog_SFINCS_obspoints.yml"),
                    join(curdir, "03_data_catalogs", "datacatalog_CF_forcing.yml")]

#datacat = ['../../Workflows/03_data_catalogs/datacatalog_general.yml']
data_catalog = DataCatalog(data_libs = data_cats)

wflow_region = gpd.read_file(os.path.join(models_dir_base, "wflow", "staticgeoms", "region.geojson")).to_crs("EPSG:4326")
wflow_basins = gpd.read_file(os.path.join(models_dir_base, "wflow", "staticgeoms", "basins.geojson")).to_crs("EPSG:4326")
wflow_basins_dissolved = wflow_basins.dissolve().reset_index(drop=True)
bbox_wflow_region = wflow_region.total_bounds  # [minx, miny, maxx, maxy]

# SFINCS model region
model_region_gdf = gpd.read_file(os.path.join(runs_dir_sfincs, "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_lisboa_2020",
                                              "gis/region.geojson")).to_crs("EPSG:4326") 

# Load own background shape for plotting that aligns with permanent water mask: permanent water where water occurence > 5%
gdf_background = gpd.read_file(
    "C:/Code/clim_pop_change_impact_TC_Idai/Socio-economic/data/gis/case_study_region_background.geojson",
    driver="GeoJSON").to_crs(model_region_gdf.crs)

# settings for plotting
projection = 32736
utm_crs = ccrs.UTM(zone=36, southern_hemisphere=True)


#%%
# Define different model directories# Factual model with updated landuse for 2020
factual_model_dir    = join(runs_dir_sfincs, "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_lisboa_2020")

# Counterfactual landuse model
CF_landuse_model_dir = join(runs_dir_sfincs, "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_lisboa_2000")

# Counterfactual climate model with 8% of precipitation removed acc. to Clausius Clapeyron at 1.1 degree of warming, 
# 14 cm of sea level rise (SLR) removed and with 10% reduced maximum tropical cyclone wind speeds
CF_climate_model_dir = join(runs_dir_sfincs, "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.1_era5_hourly_spw_IBTrACS_CF-5_lisboa_2020")

# Factual model with ESA Worldcover merged wiith Lisboa 2020
lisboa_esa_F_dir     = join(runs_dir_sfincs, "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_lisboa_2020_esa_merged")

#%%
# Load in models and model region
print("Loading model data and model region")
# Factual model
mod_F  = SfincsModel(factual_model_dir, data_libs=data_cats, mode="r")
ds_his_F = xr.open_dataset(join(factual_model_dir,"sfincs_his.nc"), engine='netcdf4')

# CF landuse models 
mod_CF  = SfincsModel(CF_landuse_model_dir, data_libs=data_cats, mode="r")
ds_his_CF = xr.open_dataset(join(CF_landuse_model_dir,"sfincs_his.nc"), engine='netcdf4')

# CF climate models 
mod_CF_clim = SfincsModel(CF_climate_model_dir, data_libs=data_cats, mode="r")
ds_his_CF_clim = xr.open_dataset(join(CF_climate_model_dir,"sfincs_his.nc"), engine='netcdf4')

# Factual model with ESA Worldcover merged with Lisboa 2020
mod_F_esa  = SfincsModel(lisboa_esa_F_dir, data_libs=data_cats, mode="r")
ds_his_F_esa = xr.open_dataset(join(lisboa_esa_F_dir,"sfincs_his.nc"), engine='netcdf4')


#%%
print("Loading flood maps")
# Load flood map which is already downscaled, masked for permanent water and represents a cells as flooded from 0.05 m or more
model_dirs = {
    "F": factual_model_dir,
    "CF_landuse": CF_landuse_model_dir,
    "CF_climate": CF_climate_model_dir,
    "F_esa": lisboa_esa_F_dir,
}

hmax = {}
for key, d in model_dirs.items():
    da = rxr.open_rasterio(Path(d) / "plot_output" /"floodmap.tif", masked=True, chunks="auto")

    if "band" in da.dims:
        da = da.squeeze("band", drop=True)

    hmax[key] = da


#%%
hmax_diff_lisboa = hmax['F'] - hmax['CF_landuse']
hmax_diff_esa = hmax['F'] - hmax['F_esa']

vmax = max(hmax_diff_lisboa.quantile(0.99).values, hmax_diff_esa.quantile(0.99).values)
vmin = min(hmax_diff_lisboa.quantile(0.01).values, hmax_diff_esa.quantile(0.01).values)
norm = TwoSlopeNorm(vcenter=0, vmin=vmin, vmax=vmax)

cmap = "RdBu_r"  # diverging colormap for differences

# Plot max flood depth maps for comparison
fig, axes = plt.subplots(1,2, figsize=(10, 8), dpi=300, subplot_kw={"projection": utm_crs},
                         constrained_layout=True)

# Plot hmax masked
im1 = hmax_diff_lisboa.plot.pcolormesh(ax=axes[0], cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False,
                                       norm=norm, transform=utm_crs, rasterized=True)
im2 = hmax_diff_esa.plot.pcolormesh(ax=axes[1], cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False,
                                     norm=norm, transform=utm_crs, rasterized=True)

axes[0].set_title("Difference in max flood depth Lisboa F - CF landuse", fontsize=10)
axes[1].set_title("Difference in max flood depth Lisboa F - ESA F", fontsize=10)

minx, miny, maxx, maxy = model_region_gdf.bounds.minx.item(), model_region_gdf.bounds.miny.item(), model_region_gdf.bounds.maxx.item(), model_region_gdf.bounds.maxy.item()
for ax in axes:
    # Add gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.2, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    # Plot valid landmask
    gdf_background.plot(ax=ax, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)
    # Add model region
    model_region_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())
    # Set extent (based on actual lat/lon coordinates)
    ax.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

cbar = fig.colorbar(im1, ax=axes, orientation="vertical", pad=0.04, shrink=0.6)
cbar.set_label("Difference in maximum flood depth [m]", fontsize=10)

plt.show()


#%%
hmax_diff_clim = hmax['F'] - hmax['CF_climate']

vmax = max(hmax_diff_lisboa.quantile(0.99).values, hmax_diff_esa.quantile(0.99).values, hmax_diff_clim.quantile(0.99).values)
vmin = min(hmax_diff_lisboa.quantile(0.01).values, hmax_diff_esa.quantile(0.01).values, hmax_diff_clim.quantile(0.01).values)
norm = TwoSlopeNorm(vcenter=0, vmin=vmin, vmax=vmax)

cmap = "RdBu_r"  # diverging colormap for differences

# Plot max flood depth maps for comparison
fig, axes = plt.subplots(1,3, figsize=(12, 8), dpi=300, subplot_kw={"projection": utm_crs},
                         constrained_layout=True, sharey=True)

# Plot hmax masked
im1 = hmax_diff_lisboa.plot.pcolormesh(ax=axes[0], cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False,
                                       norm=norm, transform=utm_crs, rasterized=True)
im2 = hmax_diff_clim.plot.pcolormesh(ax=axes[1], cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False,
                                     norm=norm, transform=utm_crs, rasterized=True)
im3 = hmax_diff_esa.plot.pcolormesh(ax=axes[2], cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False,
                                     norm=norm, transform=utm_crs, rasterized=True)

axes[0].set_title("Difference in max flood depth Lisboa F - CF landuse", fontsize=10)
axes[1].set_title("Difference in max flood depth Lisboa F - CF climate", fontsize=10)
axes[2].set_title("Difference in max flood depth Lisboa F - ESA CF", fontsize=10)

minx, miny, maxx, maxy = model_region_gdf.bounds.minx.item(), model_region_gdf.bounds.miny.item(), model_region_gdf.bounds.maxx.item(), model_region_gdf.bounds.maxy.item()
for i, ax in enumerate(axes):
    gl = ax.gridlines(draw_labels=True, linewidth=0.2, color='gray', alpha=0.5,
                      linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    # No longitude labels on top
    gl.bottom_labels = True
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}

    # Plot valid landmask
    gdf_background.plot(ax=ax, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)
    # Add model region
    model_region_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())
    # Set extent (based on actual lat/lon coordinates)
    ax.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

cbar = fig.colorbar(im1, ax=axes, orientation="vertical", pad=0.04, shrink=0.6)
cbar.set_label("Difference in maximum flood depth [m]", fontsize=10)

plt.show()



#%%
##################################################################
######################### WFLOW MODELS ###########################
##################################################################
# Loading wflow models for comparison
wflow_lisboa_vito_F_event_path            = os.path.join(runs_dir_wflow, "event_precip_era5_hourly_zarr_CF0_lisboa_2020")
wflow_lisboa_F_notmerged_small_event_path = os.path.join(runs_dir_wflow, "event_precip_era5_hourly_zarr_CF0_lisboa_2020_notmerged_small")
wflow_lisboa_esa_F_path                   = os.path.join(runs_dir_wflow, "event_precip_era5_hourly_zarr_CF0_lisboa_2020_esa_merged")

# Comparing wflow models
mod_lisboa_vito_F      = WflowModel(root=join(wflow_lisboa_vito_F_event_path, "events"), data_libs=data_cats, mode="r")
mod_lisboa_esa_F       = WflowModel(root=join(wflow_lisboa_esa_F_path, "events"), data_libs=data_cats, mode="r")
mod_lisboa_F_notmerged = WflowModel(root=join(wflow_lisboa_F_notmerged_small_event_path, "events"), data_libs=data_cats, mode="r")

dict_wflow_models = {
    "lisboa_vito_F": mod_lisboa_vito_F,
    "lisboa_esa_F": mod_lisboa_esa_F,
}

# %%
# Check results
df_dis_wflow_mod = {}

for name, mod in dict_wflow_models.items():
    mod.read()
    df_dis_wflow = mod.results['netcdf']['Q'].to_pandas()
    df_dis_wflow_mod[name] = df_dis_wflow


# # %%
# # Plot discharge at specific location
# plt.figure()

# for name, mod in dict_wflow_models.items():
#     df_wflow = df_dis_wflow_mod[name]
#     ax = df_wflow["3"].plot(label=name)

# plt.ylabel("Discharge (m³/s)")
# plt.legend()


#%%
# Get gauge IDs from first model
gauges = sorted(next(iter(df_dis_wflow_mod.values())).columns, key=int)

ncols = 2
nrows = int(np.ceil(len(gauges) / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(12, 3 * nrows), sharex=True,
                         constrained_layout=True)
fig.subplots_adjust(wspace=0.08, hspace=0.05)

axes = axes.flatten()

for ax, gauge in zip(axes, gauges):
    for model_name, df in df_dis_wflow_mod.items():
        ax.plot(df.index, df[gauge], label=model_name)

    ax.set_title(f"Gauge {gauge}")
    ax.grid(alpha=0.3)
    ax.set_xlim(df.index.min(), df.index.max())
    ax.xaxis.set_major_locator(mdates.DayLocator()) # interval=5
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d"))

for ax in axes[::2]:  # column 0
    ax.set_ylabel("Discharge [m³/s]")

# Remove unused axes if odd number of gauges
for ax in axes[len(gauges):]:
    fig.delaxes(ax)

# X-label only on bottom row
for ax in axes[-ncols:]:
    ax.set_xlabel("Day in March 2019")

axes[0].legend()

plt.show()


#%%
mod_ini_lisboa_2020 = WflowModel(root=os.path.join(models_dir_base, "wflow_lisboa_2020"), mode="r+", config_fn=os.path.join(models_dir_base, "wflow_lisboa_2020", "wflow_sbm.toml"))

gdf_rivers = mod_ini_lisboa_2020.geoms["rivers"]
q_locs = mod_ini_lisboa_2020.geoms["gauges_locs"]

# Plotting
fig, ax = plt.subplots()
mod_ini_lisboa_2020.geoms["rivers"].plot(ax=ax)
q_locs.plot(ax=ax, color="red")
for idx, row in q_locs.iterrows():
    plt.text(row.geometry.x, row.geometry.y, q_locs.loc[idx,'index'], fontsize=9, ha='right')
plt.show()


#%%
# Plotting the river network and discharge gauges
fig, ax = plt.subplots(figsize=(8, 8))

# River network
gdf_rivers.plot(ax=ax, color="steelblue", linewidth=1, alpha=0.8, zorder=1)
model_region_gdf.boundary.plot(ax=ax, color="grey", linewidth=1.5, zorder=0)
wflow_basins_dissolved.boundary.plot(ax=ax, color="blue", linewidth=1.2, linestyle="--", zorder=0)

# Gauges
q_locs_plot = q_locs.copy()
offset = 0.01  # adjust to your CRS

q_locs_plot.loc[q_locs_plot["index"] == 1, "geometry"] = (q_locs_plot.loc[q_locs_plot["index"] == 1, "geometry"
                                                                          ].translate(xoff=-offset, yoff=offset))
q_locs_plot.loc[q_locs_plot["index"] == 9, "geometry"] = (q_locs_plot.loc[q_locs_plot["index"] == 9, "geometry"
                                                                          ].translate(xoff=offset, yoff=0))
q_locs_plot.plot(ax=ax, color="crimson", markersize=20, edgecolor="black")

# Labels with white background
label_offsets = {1: (-10, 5), 9: (5, -5)}
for idx, row in q_locs.iterrows():
    gauge_id = row["index"]

    dx, dy = label_offsets.get(gauge_id, (5, 5))

    ax.annotate(str(gauge_id), (row.geometry.x, row.geometry.y), xytext=(dx, dy), 
                textcoords="offset points", fontsize=9, fontweight="bold",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=0.2),
                zorder=4)

ax.set_title("Wflow river network and discharge gauges", fontsize=12, fontweight="bold")
# ax.set_axis_off()

# Custom legend
legend_elements = [Line2D([0], [0], color="steelblue", lw=2, label="River network"),
                   Line2D([0], [0], color="blue", lw=2, linestyle="--", label="Wflow basin boundary"),
                   Line2D([0], [0], color="grey", lw=2, label="SFINCS region boundary"),
                   Line2D([0], [0], marker="o", color="w", markerfacecolor="crimson",
                          markeredgecolor="black", markersize=8, label="Gauges")]
ax.legend(handles=legend_elements, loc="lower right", frameon=True)

plt.tight_layout()
plt.show()


#%%
##################################################################
######################## LAND USE MAPS ###########################
##################################################################
# Loading wflow models for comparison
lisboa_vito_merged_F   = os.path.join(runs_dir_wflow, "event_precip_era5_hourly_zarr_CF0_lisboa_2020")
lisboa_notmerged_small = os.path.join(runs_dir_wflow, "event_precip_era5_hourly_zarr_CF0_lisboa_2020_notmerged_small")
lisboa_esa_merged_F    = os.path.join(runs_dir_wflow, "event_precip_era5_hourly_zarr_CF0_lisboa_2020_esa_merged")


# Comparing land use maps from and for the different wflow models
ds_lisboa_vito_merged_F = xr.open_dataset(os.path.join(lisboa_vito_merged_F, "staticmaps.nc"))
ds_lisboa_notmerged_small = xr.open_dataset(os.path.join(lisboa_notmerged_small, "staticmaps.nc"))
ds_lisboa_esa_merged_F = xr.open_dataset(os.path.join(lisboa_esa_merged_F, "staticmaps.nc"))

lulc_lisboa_vito_F = ds_lisboa_vito_merged_F["wflow_landuse"]
lulc_notmerged_small = ds_lisboa_notmerged_small["wflow_landuse"]
lulc_esa_merged_F = ds_lisboa_esa_merged_F["wflow_landuse"]

dict_landuse_wflow = {
    "lisboa_vito_F": lulc_lisboa_vito_F,
    "lisboa_notmerged_small": lulc_notmerged_small,
    "lisboa_esa_merged_F": lulc_esa_merged_F
}

# %%
# Read reclassified but unmerged Lisboa 2020 land use map
lisboa_2020_notmerged = data_catalog.get_rasterdataset("lisboa_2020_notmerged")

# Match Wflow grid
lisboa_match = lisboa_2020_notmerged.rio.reproject_match(
    dict_landuse_wflow["lisboa_vito_F"], resampling=Resampling.nearest)
lisboa_match = lisboa_match.where(lisboa_match != 255, np.nan)

# Make polygons of valid Lisboa coverage for plotting
valid_lisboa = lisboa_match.notnull().astype("uint8")

# Mark 0 as nodata so only valid-area polygons are returned
valid_lisboa = valid_lisboa.rio.write_nodata(0)
gdf_valid = valid_lisboa.raster.vectorize()

# Keep only polygons representing valid Lisboa coverage
gdf_valid = gdf_valid[gdf_valid["value"] == 1]


#%%
# Project to an equal-area or UTM CRS
gdf_valid_proj = gdf_valid.to_crs(32736)  # example UTM zone, adapt if needed
wflow_proj = wflow_basins_dissolved.to_crs(32736)

intersection = gdf_valid_proj.geometry.union_all().intersection(
    wflow_proj.geometry.iloc[0])

area_lisboa_km2 = intersection.area / 1e6
area_basin_km2 = wflow_proj.geometry.iloc[0].area / 1e6
pct_coverage = 100 * area_lisboa_km2 / area_basin_km2
pct_merged = 100 - pct_coverage

print(f"Lisboa coverage: {area_lisboa_km2:.1f} km²")
print(f"Basin area: {area_basin_km2:.1f} km²")
print(f"Coverage: {pct_coverage:.1f}%")
print(f"Percentage of basin area that is merged: {pct_merged:.1f}%")

#%%
# -----------------------
# Classes & colours
# -----------------------
reclass = {
    10: 0,   # Tree cover -> Forest
    114: 0,  # Forest -> Forest

    90: 6,   # Wetland -> Wetland & Mangroves
    95: 6,   # Mangroves -> Wetland & Mangroves
    122: 6,  # Mangrove -> Wetland & Mangroves

    20: 1,
    30: 2,
    40: 3,
    50: 4,
    60: 5,
    126: 0,
}

class_labels = {
    0:  "Forest",
    1:  "Shrubland",
    2:  "Grassland",
    3:  "Cropland",
    4:  "Settlement",
    5:  "Other Land",
    6:  "Wetland & Mangroves",
    # 7:  "Forestry Plantation",
}

class_colors = {
    0: "#1B7837",
    1:  "#A6B96F",
    2:  "#F3E2C7",
    3:  "#F1C40F",
    4:  "#D7191C",
    5:  "#BDBDBD",
    # 90:  "#2C7FB8",
    6: "#1FB5AA",
    # 126: "#8E63CE",
}

vals = np.array(sorted(class_labels.keys()))
bounds = np.concatenate([[vals.min()-0.5], vals[:-1]+0.5, [vals.max()+0.5]])
norm = BoundaryNorm(bounds, cmap.N)

colors = [class_colors[v] for v in vals]
cmap = ListedColormap(colors)

lisboa_match_plot = xr.apply_ufunc(lambda x: np.vectorize(lambda v: reclass.get(v, v))(x),
                                   lisboa_match, vectorize=True)

dict_landuse_wflow_plot = {}
for key, da in dict_landuse_wflow.items():
    dict_landuse_wflow_plot[key] = xr.apply_ufunc(lambda x: np.vectorize(lambda v: reclass.get(v, v))(x),
                                                  da, vectorize=True)

xmin, ymin, xmax, ymax = wflow_basins_dissolved.total_bounds  # [minx, miny, maxx, maxy]


#%%
fig, axes = plt.subplots(1, 3, figsize=(12, 6), sharey=True, constrained_layout=True)

im = lisboa_match_plot.plot(ax=axes[0], cmap=cmap, norm=norm, add_colorbar=False)
# dict_landuse_wflow_plot["lisboa_notmerged_small"].plot(ax=axes[1], cmap=cmap, norm=norm, add_colorbar=False)
dict_landuse_wflow_plot["lisboa_vito_F"].plot(ax=axes[1], cmap=cmap, norm=norm, add_colorbar=False)
dict_landuse_wflow_plot["lisboa_esa_merged_F"].plot(ax=axes[2], cmap=cmap, norm=norm, add_colorbar=False)

# cbar = fig.colorbar(im, ax=axes, orientation="vertical", pad=0.02, shrink=0.6)
# cbar.set_ticks(vals)
# cbar.set_ticklabels([class_labels[v] for v in vals])

for i, ax in enumerate(axes.flat):
    if i != 0:
        ax.set_ylabel("")

for ax in axes.flat:
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
for ax in [axes[1], axes[2]]:
    gdf_valid.boundary.plot(ax=ax, color="black", linewidth=1.2, linestyle="--")
    wflow_basins_dissolved.boundary.plot(ax=ax, color="blue", linewidth=1.2, linestyle="--")

wflow_basins.boundary.plot(ax=axes[0], color="blue", linewidth=1.2, linestyle="--")

axes[0].set_title("Lisboa Only")
# axes[1].set_title("Lisboa Only wflow")
axes[1].set_title("Lisboa + Vito Factual")
axes[2].set_title("Lisboa + ESA Factual")

for label, ax in zip(["(a)", "(b)", "(c)"], axes):
    ax.text(0.02, 0.98, label, transform=ax.transAxes, fontsize=12, fontweight="bold",
            va="top", ha="left", bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))

# Add legend for the lisboa and wflow boundaries
legend_elements = [Line2D([0], [0], color="black", linestyle="--", linewidth=1.5,
                          label="Original Lisboa coverage"), 
                   Line2D([0], [0], color="blue", linestyle="-", linewidth=1.5,
                          label="Wflow basin boundary"),]
axes[2].legend(handles=legend_elements, loc="upper right", frameon=True, fontsize=9)

cbar = fig.colorbar(im, ax=axes,
                    orientation="vertical",
                    pad=0.02,
                    shrink=0.6)

cbar.set_ticks(vals)
cbar.set_ticklabels([class_labels[v] for v in vals])

for tick in cbar.ax.get_yticklabels():
    tick.set_verticalalignment("center")

plt.show()




# %%
# Merging ESA Worldcover and Lisboa 2020 land use maps for the wflow model region
esa = data_catalog.get_rasterdataset("esa_worldcover", bbox=bbox_wflow_region, buffer=2)
lisboa_2020_notmerged = data_catalog.get_rasterdataset("lisboa_2020_notmerged", bbox=bbox_wflow_region, buffer=5)

esa_wflow_region_path = r"p:/11210471-001-compass/MSc_internship/Poppy/wflow_staticgeoms/esa_wflow_region.tif"
if not os.path.exists(esa_wflow_region_path):
    esa.raster.to_raster(esa_wflow_region_path)

# Save the merged land use map to a new file
lisboa_2020_esa_path = r"p:/11210471-001-compass/MSc_internship/Poppy/LULC_data/lisboa_2020_esa_merged.tif"

if not os.path.exists(lisboa_2020_esa_path):
    # Make an empty template raster with the Lisboa resolution and extent as the Wflow model region
    xmin, ymin, xmax, ymax = bbox_wflow_region

    # Optional buffer
    buffer = 0.05
    xmin -= buffer
    ymin -= buffer
    xmax += buffer
    ymax += buffer

    dx, dy = lisboa_2020_notmerged.rio.resolution()

    nx = int(np.ceil((xmax - xmin) / abs(dx)))
    ny = int(np.ceil((ymax - ymin) / abs(dy)))

    x = xmin + abs(dx) * (0.5 + np.arange(nx))
    y = ymax - abs(dy) * (0.5 + np.arange(ny))

    template = xr.DataArray(np.full((ny, nx), np.nan, dtype=np.float32), coords={"y": y, "x": x}, 
                            dims=("y", "x")).rio.write_crs(lisboa_2020_notmerged.rio.crs)

    lisboa_match = lisboa_2020_notmerged.rio.reproject_match(template,
                                                            resampling=Resampling.nearest)

    esa_match = esa.rio.reproject_match(template, resampling=Resampling.nearest)

    merged = xr.where(lisboa_match != lisboa_match.rio.nodata, lisboa_match, esa_match)
    merged = merged.rio.write_nodata(0)
    merged = merged.rio.write_crs(lisboa_2020_notmerged.rio.crs)
    merged.raster.to_raster(lisboa_2020_esa_path)

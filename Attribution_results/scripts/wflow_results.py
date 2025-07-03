# %%
import os
from os.path import join
import pandas as pd
import numpy as np
from datetime import datetime as datetime
import geopandas as gpd
import xarray as xr

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker

import contextily as ctx
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from hydromt import DataCatalog
from hydromt_sfincs import SfincsModel, utils
from hydromt_wflow import WflowModel

import platform
from pathlib import Path

#%%
if platform.system() == "Windows":
    P = Path("p:/")
else:
    # adjust this if your cluster mounts “P:” somewhere else, e.g. “/mnt/p”
    P = Path("/p/")

# Now build your base directories
BASE        = P / "11210471-001-compass"
BASE_RUNS   = P / "11210471-001-compass" / "03_Runs"  / "sofala" / "Idai" / "wflow"
BASE_SFINCS = P / "11210471-001-compass" / "03_Runs"  / "sofala" / "Idai" / "sfincs" 
BASE_MODELS = P / "11210471-001-compass" / "02_Models" / "sofala" / "Idai" / "wflow"
BASE_DATA   = P / "11210471-001-compass" / "01_Data"

def make_model(subfolder: str):
    root = BASE_RUNS / subfolder / "events"
    cfg  = root / "wflow_sbm.toml"
    return WflowModel(root=str(root), mode="r", config_fn=str(cfg))

#%%
#Load mdoels
mod_ini = WflowModel(root=str(BASE_MODELS), mode="r+", config_fn=str(BASE_MODELS / "wflow_sbm.toml"))

mod_base = make_model("event_precip_era5_hourly_zarr_CF0_oldsettings")
mod_fplns = make_model("event_precip_era5_hourly_zarr_CF0_floodplains")
mod_fplns_f_ = make_model("event_precip_era5_hourly_zarr_CF0_floodplains_f_")
mod_f_soilthick = make_model("event_precip_era5_hourly_zarr_CF0_f_soilthick2")
mod_f_soilthick_cnst = make_model("event_precip_era5_hourly_zarr_CF0_f_soilthick")
mod_fplns_f_soilthick = make_model("event_precip_era5_hourly_zarr_CF0_floodplain_f_soilthick")
mod_fplns_f_soilthick_chirps = make_model("event_precip_chirps_CF0")
mod_fplns_f_soilthick_era5_daily = make_model("event_precip_era5_daily_CF0")

# %%
# # Get station IDs
station_ids = list(mod_base.results['netcdf']['Q'].Q_gauges_locs.values)

# Dictionary of models to loop over
models = {
    "mod_base": mod_base,
    "mod_fplns": mod_fplns,
    "mod_fplns_f_": mod_fplns_f_,
    "mod_fplns_f_soilthick": mod_fplns_f_soilthick,
    "mod_f_soilthick": mod_f_soilthick,
    "mod_fplns_f_soilthick_chirps": mod_fplns_f_soilthick_chirps
}

#%%
# # Plot all models for all stations
# for model_name, model in models.items():
#     fig, ax = plt.subplots()
#     for st in station_ids:
#         model.results['netcdf']['Q'].sel(Q_gauges_locs=st).plot(ax=ax, label=f'Gauge {st}')
#     ax.legend(title='Discharge Gauges')
#     ax.set_title(model_name)

#%%
# Plot comparison for each gauge individually
compare_models = {
    "Base": mod_base,
    "Incl floodplains": mod_fplns,
    "Floodplains with f_": mod_fplns_f_,
    "f_ and soilthickness": mod_f_soilthick,
    "Floodplains, f_ and soilthickness": mod_fplns_f_soilthick,
    "Floodplains, f_ and soilthickness and CHIRPS": mod_fplns_f_soilthick_chirps,
    "Floodplains, f_ and soilthickness and ERA5 daily": mod_fplns_f_soilthick_era5_daily
    # Uncomment to include constant variant:
    # "f_ and soilthickness const": mod_f_soilthick_cnst
}

for gauge_id in ['1', '2']:
    fig, ax = plt.subplots()
    for label, model in compare_models.items():
        model.results['netcdf']['Q'].sel(Q_gauges_locs=gauge_id).plot(ax=ax, label=label)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.35), ncol=2)
    ax.set_title(f"Discharge Gauge {gauge_id}")


# %% ##############################################################
############################## GRDC DATA ##########################
## ################################################################
# Read in GRDC data
path = str(BASE_DATA / "GRDC" / "GRDC-Daily.nc")
GRDC_data = xr.open_dataset(path)

# Mod wflow 30 yr
wflow_30yr_new =  make_model("event_precip_era5_hourly_zarr_CF0_30yr")
# %%
# Get the indices where river_name matches your selection
names = GRDC_data['river_name'].values.astype(str)
wanted = ['BUDZI', 'RIOPUNGOE']
idx = [i for i, n in enumerate(names) if n in wanted]

# Select by index along the 'id' dimension
selected = GRDC_data.isel(id=idx)

# Now plot
lats = selected['geo_y'].values
lons = selected['geo_x'].values

plt.figure(figsize=(6, 6))
ax = plt.gca()
# Plot points (note: contextily expects lon=x, lat=y in Web Mercator)
plt.scatter(lons, lats, marker='o', color='red', s=10)
for i, name in enumerate(selected['river_name'].values):
    plt.text(lons[i] - 0.07, lats[i], name, fontsize=8, fontweight='bold', ha='right')

# Set axis labels and title
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Selected GRDC Rivers in Blue')

# Convert to Web Mercator for contextily
gdf_stations = gpd.GeoDataFrame(
    {'name': selected['river_name'].values},
    geometry=gpd.points_from_xy(lons, lats),
    crs="EPSG:4326"
)

# Getting the model regions
gdf_wflow = gpd.read_file(str(BASE_MODELS / "staticgeoms" / "basins.geojson"))
gdf_sfincs = gpd.read_file(str(BASE / "02_Models" / "sofala" / "Idai" / "sfincs" / "gis" / "region.geojson"))
gdf_wflow = gdf_wflow.to_crs("EPSG:4326")
gdf_sfincs = gdf_sfincs.to_crs("EPSG:4326")

# Get the SFINCS region geometry (assuming it's a single polygon)
sfincs_geom = gdf_sfincs.unary_union

# Calculate distance from each station to the SFINCS region (in degrees, since EPSG:4326)
gdf_stations['distance'] = gdf_stations.geometry.apply(lambda x: x.distance(sfincs_geom))

# Select the two closest stations
closest_stations = gdf_stations.nsmallest(2, 'distance')

# Plot SFINCS region and set up legend entry
gdf_sfincs.plot(ax=ax, edgecolor='pink', facecolor='pink', linewidth=1, alpha=0.5, zorder=3)
sfincs_patch = mpatches.Patch(facecolor='pink', edgecolor='pink', alpha=0.5, label="SFINCS Region")

# Plot wflow basins region and set up legend entry
gdf_wflow.plot(ax=ax, edgecolor='lightskyblue', facecolor='lightskyblue', linewidth=1, alpha=0.5, zorder=2)
wflow_patch = mpatches.Patch(facecolor='lightskyblue', edgecolor='lightskyblue', alpha=0.5, label="Wflow Basins")

# Plot all and the two closest stations
closest_stations.plot(ax=ax, marker='o', color='blue', zorder=5, label='Closest to SFINCS')
gdf_stations.plot(ax=ax, marker='o', color='red', zorder=2)
# wflow_30yr_new.results['netcdf']['Q'].sel(Q_gauges_locs='1').plot

#The folder is the location of the wflow model in 02_Models
gdf_gauges = mod_ini.geoms["gauges_locs"]

gauge1 = gdf_gauges[gdf_gauges['index'] == 1]
gauge2 = gdf_gauges[gdf_gauges['index'] == 2]

# Plot gauge 1 and 2
gauge1.plot(ax=ax, color='#1f77b4', markersize=30, label='wflow Gauge 1', zorder=10)
gauge2.plot(ax=ax, color='#2ca02c', markersize=30, label='wflow Gauge 2', zorder=10)
mod_ini.geoms["rivers"].plot(ax=ax, color='white')

# Add text labels at the gauge locations
ax.text(
    gauge1.geometry.x.iloc[0]+0.1, gauge1.geometry.y.iloc[0],
    str(gauge1['index'].iloc[0]), fontsize=9, ha='right', zorder=10,
    fontweight='bold')
ax.text(
    gauge2.geometry.x.iloc[0]+0.1, gauge2.geometry.y.iloc[0],
    str(gauge2['index'].iloc[0]), fontsize=9, ha='right', zorder=10,
    fontweight='bold')

# Set the extent and add OSM basemap
ax.set_xlim(gdf_stations.geometry.x.min()-0.2, gdf_stations.geometry.x.max()+1)
ax.set_ylim(gdf_stations.geometry.y.min()-0.5, gdf_stations.geometry.y.max()+0.5)
ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik, attribution=False, crs=gdf_stations.crs)

plt.show()
# %%
# Get the indices of the two closest stations in the original dataset
# Get coordinates from closest_stations
closest_coords = set(zip(closest_stations.geometry.x, closest_stations.geometry.y))

# Get coordinates from selected and their indices
selected_coords = list(zip(GRDC_data['geo_x'].values, GRDC_data['geo_y'].values))

# Find indices in 'selected' that match closest_stations
matching_indices = [i for i, coord in enumerate(selected_coords) if coord in closest_coords]

plt.figure(figsize=(10, 5))

# Plot GRDC stations
for idx in matching_indices:
    station_name = GRDC_data['river_name'].values[idx]
    plt.plot(GRDC_data['time'].values, GRDC_data['runoff_mean'].isel(id=idx), label=f"GRDC {station_name}")

# Plot wflow results (add labels and plot only once)
plt.plot(
    wflow_30yr_new.results['netcdf']['Q'].time,
    wflow_30yr_new.results['netcdf']['Q'].sel(Q_gauges_locs='1'),
    label='Q1', color = '#1f77b4'
)
plt.plot(
    wflow_30yr_new.results['netcdf']['Q'].time,
    wflow_30yr_new.results['netcdf']['Q'].sel(Q_gauges_locs='2'),
    label='Q2', color = '#2ca02c'
)

plt.xlabel('Time')
plt.ylabel('Discharge')
plt.title('Timeseries for Two Closest GRDC Stations and wflow')
plt.legend()
plt.tight_layout()
plt.show()

#%%
fig, axs = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

# Find indices for BUDZI and RIOPUNGOE
names = GRDC_data['river_name'].values.astype(str)
idx_budzi = [i for i, n in enumerate(names) if n == 'BUDZI']
idx_riopungoe = [i for i, n in enumerate(names) if n == 'RIOPUNGOE']

# Plot BUDZI and Q1
axs[0].plot(
    wflow_30yr_new.results['netcdf']['Q'].time,
    wflow_30yr_new.results['netcdf']['Q'].sel(Q_gauges_locs='1'),
    label='wflow Q1', color='#1f77b4'
)
station_name = GRDC_data['river_name'].values[45]
axs[0].plot(GRDC_data['time'].values, GRDC_data['runoff_mean'].isel(id=45),
            color='#ff7f0e', label=f"GRDC {station_name}")

# Add horizontal lines for max values
axs[0].axhline(GRDC_data['runoff_mean'].isel(id=45).max(), color='#ff7f0e', linestyle='--', linewidth=1, label='GRDC max')
axs[0].axhline(wflow_30yr_new.results['netcdf']['Q'].sel(Q_gauges_locs='1').max(), color='#1f77b4', linestyle='--', linewidth=1, label='wflow Q1 max')

axs[0].set_ylabel('Discharge')
axs[0].set_title('BUDZI and wflow Q1')
axs[0].legend()

# Plot RIOPUNGOE and Q2
axs[1].plot(
    wflow_30yr_new.results['netcdf']['Q'].time,
    wflow_30yr_new.results['netcdf']['Q'].sel(Q_gauges_locs='2'),
    label='wflow Q2', color='#2ca02c'
)

station_name = GRDC_data['river_name'].values[38]
axs[1].plot(GRDC_data['time'].values, GRDC_data['runoff_mean'].isel(id=38), 
            color='#9467bd',label=f"GRDC {station_name}")

# Add horizontal lines for max values
axs[1].axhline(GRDC_data['runoff_mean'].isel(id=38).max(), color='#9467bd', linestyle='--', linewidth=1, label='GRDC max')
axs[1].axhline(wflow_30yr_new.results['netcdf']['Q'].sel(Q_gauges_locs='2').max(), color='#2ca02c', linestyle='--', linewidth=1, label='wflow Q2 max')


axs[1].set_xlabel('Time')
axs[1].set_ylabel('Discharge')
axs[1].set_title('RIOPUNGOE and wflow Q2')
axs[1].legend()

plt.tight_layout()
plt.show()

# %%
# From Albrecht
# Read and drop the first row (header example)
# Q_Day_1494100 = pd.read_csv(
#     r"c:\Code\2025-07-01_12-51\1494100_Q_Day.Cmd.txt",
#     encoding="latin1",
#     delimiter=";",
#     skiprows=5,
#     names=["date", "time", "Q"],
#     na_values=["-999.000", "-999"],
#     comment="#"
# )
# Q_Day_1494100 = Q_Day_1494100.iloc[1:].reset_index(drop=True)
# %%




##########################################################
############### PLOTTING SFINCS RESULTS ##################
##########################################################
# %%
def make_sfincs_model(subfolder: str):
    root = BASE_SFINCS / subfolder
    return SfincsModel(root=str(root), mode="r", data_libs=datacat)

datacat = [
        '../../Workflows/03_data_catalogs/datacatalog_general.yml',
        '../../Workflows/03_data_catalogs/datacatalog_SFINCS_obspoints.yml',
        '../../Workflows/03_data_catalogs/datacatalog_SFINCS_coastal_coupling.yml',
        '../../Workflows/03_data_catalogs/datacatalog_CF_forcing.yml'
        ]

data_catalog = DataCatalog(data_libs = datacat)

#%%
# Load in the SFINCS models
mod_old        = make_sfincs_model("event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_base")
mod_wflow_all  = make_sfincs_model("event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_updatedwflow_all")
mod_noFpln     = make_sfincs_model("event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_updatedwflow_nofloodplains")
mod_ERA5_daily = make_sfincs_model("event_tp_era5_daily_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0")
mod_chirps     = make_sfincs_model("event_tp_chirps_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0")


#%%
# we set a threshold to mask minimum flood depth
hmin = 0.05

for mod, model_name in [
    (mod_old, "old"),
    (mod_wflow_all, "wflow all"),
    (mod_noFpln, "no floodplains"),
    (mod_ERA5_daily, "ERA5_daily"),
    (mod_chirps, "CHIRPS daily"),
]:
    # compute the maximum over all time steps
    da_zsmax = mod.results["zsmax"].max(dim="timemax")
        
    # downscale the floodmap
    depfile  = str(BASE_SFINCS / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_base" / "subgrid" / "dep_subgrid.tif")
    da_dep   = mod.data_catalog.get_rasterdataset(depfile)

    da_hmax = utils.downscale_floodmap(
        zsmax=da_zsmax,
        dep=da_dep,
        hmin=hmin,
        # floodmap_fn=join(sfincs_root, "gis/floodmap.tif") # uncomment to save floodmap to <mod.root>/floodmap.tif
        )
        
    # GSWO dataset to mask permanent water by first geprojecting it to the subgrid of hmax
    sfincs_region = mod.region
    gwso_region   = data_catalog.get_rasterdataset("gswo", geom=sfincs_region, buffer=1000)
    gswo_mask     = gwso_region.raster.reproject_like(da_hmax, method="max")
    # permanent water where water occurence > 5%
    da_hmax_masked = da_hmax.where(gswo_mask <= 5)

    # Add the name attribute for identification
    mod.results['hmax'] = da_hmax
    mod.results['hmax_masked'] = da_hmax_masked

    del da_hmax, da_zsmax, da_dep, gswo_mask  # Clean up to free memory


#%%
# Plot the actual flood maps
projection = mod.crs.to_epsg()

fig, ax = plt.subplots(nrows=3,
                       ncols=2,
                       figsize=(8,8),
                       dpi=300,
                       constrained_layout=True,
                       subplot_kw={"projection": ccrs.epsg(projection)})

fig.suptitle("Factual Max Flood Depth", fontsize=11)

# Plot hmax masked
im = mod_old.results['hmax_masked'].plot.pcolormesh(
     ax=ax[0,0], cmap="Blues", vmin=0, vmax=5.0, 
     add_colorbar=False)

im = mod_wflow_all.results['hmax_masked'].plot.pcolormesh(
     ax=ax[0,1], cmap="Blues", vmin=0, vmax=5.0, 
     add_colorbar=False)

im = mod_noFpln.results['hmax_masked'].plot.pcolormesh(
     ax=ax[1,0], cmap="Blues", vmin=0, vmax=5.0, 
     add_colorbar=False)

im = mod_ERA5_daily.results['hmax_masked'].plot.pcolormesh(
     ax=ax[1,1], cmap="Blues", vmin=0, vmax=5.0, 
     add_colorbar=False)

im = mod_chirps.results['hmax_masked'].plot.pcolormesh(
     ax=ax[2,0], cmap="Blues", vmin=0, vmax=5.0, 
     add_colorbar=False)

# Add basemap and gridlines for each subplot
for i, axis in enumerate(ax.flat):
    ctx.add_basemap(
        axis,
        source=ctx.providers.Esri.WorldImagery,
        zoom=10,
        crs=mod_old.results['hmax_masked'].rio.crs,
        attribution=False
    )
    gl = axis.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, .2))
    gl.ylocator = mticker.FixedLocator(np.arange(-90, 90, .2))

# Titles
ax[0,0].set_title(f"Old wflow settings", fontsize=12)
ax[0,1].set_title(f"New wflow settings", fontsize=12)
ax[1,0].set_title(f"No 1D floodplain", fontsize=12)
ax[1,1].set_title(f"ERA5 daily", fontsize=12)
ax[2,0].set_title(f"CHIRPS daily", fontsize=12)

# Add shared colorbar with larger size and labels
cbar = fig.colorbar(im, ax=ax, orientation="vertical", 
                    fraction=0.035, pad=0.05)
cbar.set_label('Flood depth (m)', rotation=270, labelpad=10, fontsize=9)
cbar.ax.tick_params(labelsize=8)
    
fig.savefig("../figures/wflow_test_hmax_masked.png", bbox_inches='tight', dpi=300)
plt.show()



# %%
# Calculate the difference between the two hmax_masked maps (mod_new - mod_old)
hmax_diff_all = mod_wflow_all.results['hmax_masked'] - mod_old.results['hmax_masked']
hmax_diff_noFpln = mod_noFpln.results['hmax_masked'] - mod_old.results['hmax_masked']
hmax_diff_ERA5_daily = mod_ERA5_daily.results['hmax_masked'] - mod_old.results['hmax_masked']
hmax_diff_chirps = mod_chirps.results['hmax_masked'] - mod_old.results['hmax_masked']

# Plot the difference map as a new subplot
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8, 8), dpi=300,
                       constrained_layout=True,
                       subplot_kw={"projection": ccrs.epsg(projection)})

fig.suptitle("Difference in Flood Depth", fontsize=11)

# Plot the difference
im0 = hmax_diff_all.plot.pcolormesh(
    ax=ax[0,0], cmap="RdBu", vmin=-2, vmax=2, add_colorbar=False,
    transform=ccrs.PlateCarree())

im1 = hmax_diff_noFpln.plot.pcolormesh(
    ax=ax[0,1], cmap="RdBu", vmin=-2, vmax=2, add_colorbar=False,
    transform=ccrs.PlateCarree())

im2 = hmax_diff_ERA5_daily.plot.pcolormesh(
    ax=ax[1,0], cmap="RdBu", vmin=-2, vmax=2, add_colorbar=False,
    transform=ccrs.PlateCarree())

im3 = hmax_diff_chirps.plot.pcolormesh(
    ax=ax[1,1], cmap="RdBu", vmin=-2, vmax=2, add_colorbar=False,
    transform=ccrs.PlateCarree())

ax[0,0].set_title("Difference (All settings - Old)", fontsize=12)
ax[0,1].set_title("Difference (No floodplains - Old)", fontsize=12)
ax[1,0].set_title("Difference (ERA5 daily - Old)", fontsize=12)
ax[1,1].set_title("Difference (CHIRPS - Old)", fontsize=12)


# Add basemap and gridlines for each subplot
for i, axis in enumerate(ax.flat):
    ctx.add_basemap(
        axis,
        source=ctx.providers.Esri.WorldImagery,
        zoom=10,
        crs=mod_old.results['hmax_masked'].rio.crs,
        attribution=False
    )
    gl = axis.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, .2))
    gl.ylocator = mticker.FixedLocator(np.arange(-90, 90, .2))

# Add shared colorbar
cbar = fig.colorbar(im3, ax=ax, orientation="vertical", fraction=0.04, pad=0.05)
cbar.set_label('Difference in Flood depth (m)', rotation=270, labelpad=10, fontsize=9)
cbar.ax.tick_params(labelsize=8)

fig.savefig("../figures/wflow_test_hmax_diff.png", bbox_inches='tight', dpi=300)

plt.show()

# %%

# %%
print("Loading packages...")
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
from shapely.geometry import Point

import contextily as ctx
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from hydromt import DataCatalog
from hydromt_sfincs import SfincsModel, utils
from hydromt_wflow import WflowModel

import platform
from pathlib import Path
from geopy.distance import geodesic

from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
from hydromt.stats import skills as skillstats  # NSE, KGE, etc.

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
    print(f"Loading model: {subfolder} from {BASE_RUNS}")
    return WflowModel(root=str(root), mode="r", config_fn=str(cfg))

#%%
#Load mdoels
mod_ini = WflowModel(root=str(BASE_MODELS), mode="r+", config_fn=str(BASE_MODELS / "wflow_sbm.toml"))

#%%
mod_base = make_model("event_precip_era5_hourly_zarr_CF0_oldsettings")
mod_fplns = make_model("event_precip_era5_hourly_zarr_CF0_floodplains")
mod_fplns_f_ = make_model("event_precip_era5_hourly_zarr_CF0_floodplains_f_")
mod_f_soilthick = make_model("event_precip_era5_hourly_zarr_CF0_f_soilthick2")
mod_f_soilthick_cnst = make_model("event_precip_era5_hourly_zarr_CF0_f_soilthick")
mod_fplns_f_soilthick = make_model("event_precip_era5_hourly_zarr_CF0_floodplain_f_soilthick")
mod_fplns_f_soilthick_maxlk = make_model("event_precip_era5_hourly_zarr_CF0")
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
    "mod_f_soilthick": mod_f_soilthick,
    "mod_fplns_f_soilthick": mod_fplns_f_soilthick,
    "mod_fplns_f_soilthick_maxlk": mod_fplns_f_soilthick_maxlk,
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
    # "Floodplains, f_, soilthickness and maxleakage": mod_fplns_f_soilthick_maxlk,
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
print("Loading GRDC dataset...")
grdc_path = BASE_DATA / "GRDC" / "GRDC-Daily.nc"
GRDC_data = xr.open_dataset(str(grdc_path))

print("Initializing 30-year WFLOW model...")
def make_model_warmup(subfolder: str):
    root = BASE_RUNS / subfolder / "warmup"
    cfg = root / "wflow_sbm.toml"
    return WflowModel(root=str(root), mode="r", config_fn=str(cfg))

wflow_30yr_new = make_model_warmup("event_precip_era5_daily_CF0_30yr_1954")
wflow_30yr_1989 = make_model_warmup("event_precip_era5_hourly_zarr_CF0_30yr")
wflow_30yr_maxL06 = make_model_warmup("event_precip_era5_daily_CF0_30yr")

# %% ##############################################################
######################### SELECT GRDC STATIONS ####################
###################################################################

print("Selecting GRDC stations...")
wanted_rivers = ['BUDZI', 'RIOPUNGOE']
names = GRDC_data['river_name'].values.astype(str)
idx = [i for i, n in enumerate(names) if n in wanted_rivers]
selected = GRDC_data.isel(id=idx)

print(f"Found {len(selected.id)} matching GRDC stations.")

lats, lons = selected['geo_y'].values, selected['geo_x'].values
gdf_stations = gpd.GeoDataFrame({'name': selected['river_name'].values},
                                 geometry=gpd.points_from_xy(lons, lats),
                                 crs="EPSG:4326")


# %% ##############################################################
######################## LOAD REGIONS & MATCH #####################
###################################################################

print("Loading WFLOW and SFINCS regions...")
gdf_wflow = gpd.read_file(BASE_MODELS / "staticgeoms" / "basins.geojson").to_crs("EPSG:4326")
gdf_sfincs = gpd.read_file(BASE / "02_Models" / "sofala" / "Idai" / "sfincs" / "gis" / "region.geojson").to_crs("EPSG:4326")

print("Finding closest GRDC stations to SFINCS domain...")
sfincs_geom = gdf_sfincs.unary_union
gdf_stations['distance'] = gdf_stations.geometry.distance(sfincs_geom)
closest_stations = gdf_stations.nsmallest(2, 'distance')
print(f"Closest stations: {list(closest_stations['name'])}")


#%%
# %% ##############################################################
############# MATCH WFLOW GRID POINTS TO CLOSEST STATIONS #########
###################################################################

print("Matching GRDC stations to WFLOW river grid points...")
q_river = wflow_30yr_new.results['output']['q_river']
lat_vals, lon_vals = q_river['lat'].values, q_river['lon'].values

station_coords = [(pt.y, pt.x) for pt in closest_stations.geometry]
station_series, matched_indices, matched_gridpoints = [], [], []

for lat, lon in station_coords:
    lat_idx, lon_idx = (abs(lat_vals - lat)).argmin(), (abs(lon_vals - lon)).argmin()
    matched_lat, matched_lon = float(lat_vals[lat_idx]), float(lon_vals[lon_idx])
    matched_gridpoints.append((matched_lat, matched_lon))
    matched_indices.append((lat_idx, lon_idx))
    station_series.append(q_river[:, lat_idx, lon_idx])

print("Matched WFLOW grid points to GRDC stations.")

matches_df = pd.DataFrame({
    'station_index': closest_stations.index,
    'GRDC_station_lat': [lat for lat, _ in station_coords],
    'GRDC_station_lon': [lon for _, lon in station_coords],
    'wflow_matched_grid_lat': [lat for lat, _ in matched_gridpoints],
    'wflow_matched_grid_lon': [lon for _, lon in matched_gridpoints],
    'lat_idx': [i for i, _ in matched_indices],
    'lon_idx': [j for _, j in matched_indices],
    'distance_km': [
        geodesic((a, b), (c, d)).km for (a, b), (c, d) in zip(station_coords, matched_gridpoints)
    ]
})
print(matches_df)


print("Matching GRDC stations to WFLOW river grid points (wflow_30yr_1989)...")


#%%
# Extract river discharge and grid coordinates
q_river_1989 = wflow_30yr_1989.results['output']['q_river']
lat_vals_1989, lon_vals_1989 = q_river_1989['lat'].values, q_river_1989['lon'].values

# Match GRDC stations to closest WFLOW grid cells
station_coords_1989 = [(pt.y, pt.x) for pt in closest_stations.geometry]
station_series_1989, matched_indices_1989, matched_gridpoints_1989 = [], [], []

for lat, lon in station_coords_1989:
    lat_idx = abs(lat_vals_1989 - lat).argmin()
    lon_idx = abs(lon_vals_1989 - lon).argmin()
    matched_lat, matched_lon = float(lat_vals_1989[lat_idx]), float(lon_vals_1989[lon_idx])
    matched_gridpoints_1989.append((matched_lat, matched_lon))
    matched_indices_1989.append((lat_idx, lon_idx))
    station_series_1989.append(q_river_1989[:, lat_idx, lon_idx])

print("Matched WFLOW grid points to GRDC stations for wflow_30yr_1989.")

matches_df_1989 = pd.DataFrame({
    'station_index': closest_stations.index,
    'GRDC_station_lat': [lat for lat, _ in station_coords_1989],
    'GRDC_station_lon': [lon for _, lon in station_coords_1989],
    'wflow_matched_grid_lat': [lat for lat, _ in matched_gridpoints_1989],
    'wflow_matched_grid_lon': [lon for _, lon in matched_gridpoints_1989],
    'lat_idx': [i for i, _ in matched_indices_1989],
    'lon_idx': [j for _, j in matched_indices_1989],
    'distance_km': [
        geodesic((a, b), (c, d)).km for (a, b), (c, d) in zip(station_coords_1989, matched_gridpoints_1989)
    ]
})

print(matches_df_1989)


#%%
# Extract river discharge and grid coordinates
q_river_max = wflow_30yr_maxL06.results['output']['q_river']
lat_vals_max, lon_vals_max = q_river_max['lat'].values, q_river_max['lon'].values

# Match GRDC stations to closest WFLOW grid cells
station_coords_max = [(pt.y, pt.x) for pt in closest_stations.geometry]
station_series_max, matched_indices_max, matched_gridpoints_max = [], [], []

for lat, lon in station_coords_max:
    lat_idx = abs(lat_vals_max - lat).argmin()
    lon_idx = abs(lon_vals_max - lon).argmin()
    matched_lat, matched_lon = float(lat_vals_max[lat_idx]), float(lon_vals_max[lon_idx])
    matched_gridpoints_max.append((matched_lat, matched_lon))
    matched_indices_max.append((lat_idx, lon_idx))
    station_series_max.append(q_river_max[:, lat_idx, lon_idx])

print("Matched WFLOW grid points to GRDC stations for wflow_30yr_maxL06.")

# matches_df_max = pd.DataFrame({
#     'station_index': closest_stations.index,
#     'GRDC_station_lat': [lat for lat, _ in station_coords_max],
#     'GRDC_station_lon': [lon for _, lon in station_coords_max],
#     'wflow_matched_grid_lat': [lat for lat, _ in matched_gridpoints_max],
#     'wflow_matched_grid_lon': [lon for _, lon in matched_gridpoints_max],
#     'lat_idx': [i for i, _ in matched_indices_max],
#     'lon_idx': [j for _, j in matched_indices_max],
#     'distance_km': [
#         geodesic((a, b), (c, d)).km for (a, b), (c, d) in zip(station_coords_max, matched_gridpoints_max)
#     ]
# })


# %% ##############################################################
################### PLOT SPATIAL MATCHING OVERVIEW ################
###################################################################

# print("Plotting spatial match of GRDC and WFLOW points...")
# matched_wflow_gdf = gpd.GeoDataFrame(
#     geometry=[Point(lon, lat) for lat, lon in matches_df[['wflow_matched_grid_lat', 'wflow_matched_grid_lon']].values],
#     crs='EPSG:4326'
# )

# fig, ax = plt.subplots(figsize=(10, 10))
# gdf_wflow.plot(ax=ax, edgecolor='skyblue', facecolor='skyblue', alpha=0.5, label="Wflow Basins")
# gdf_sfincs.plot(ax=ax, edgecolor='pink', facecolor='pink', alpha=0.5, label="SFINCS Region")
# mod_ini.geoms["rivers"].plot(ax=ax, color='white', zorder=1)
# gdf_stations.plot(ax=ax, color='red', label='All GRDC Stations', zorder=10)
# closest_stations.plot(ax=ax, color='blue', label='Closest GRDC Stations', zorder=10)
# matched_wflow_gdf.plot(ax=ax, marker='x', color='green', markersize=50, label='Matched WFLOW Grid', zorder=10)

# # Add text labels to matched WFLOW grid points
# for i, point in enumerate(matched_wflow_gdf.geometry):
#     ax.text(point.x + 0.02, point.y + 0.02, f"Q{i}", fontsize=9, fontweight='bold', color='green', ha='left', zorder=10)

# for i, point in enumerate(closest_stations.geometry):
#     ax.text(point.x + 0.02, point.y - 0.02, f"G{i}", fontsize=9, fontweight='bold', color='blue', ha='left', zorder=10)

# gdf_gauges = mod_ini.geoms["gauges_locs"]
# # for i in [1, 2]:
# #     gauge = gdf_gauges[gdf_gauges['index'] == i]
# #     gauge.plot(ax=ax, markersize=30, label=f'wflow-sfincs Gauge {i}', zorder=10)
# #     ax.text(gauge.geometry.x.iloc[0] + 0.1, gauge.geometry.y.iloc[0],
# #             str(i), fontsize=9, fontweight='bold', ha='right')


# ax.set_xlim(gdf_stations.geometry.x.min() - 0.2, gdf_stations.geometry.x.max() + 1)
# ax.set_ylim(gdf_stations.geometry.y.min() - 0.5, gdf_stations.geometry.y.max() + 0.5)
# ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik, attribution=False, crs=gdf_stations.crs)
# ax.legend(loc='lower left')
# plt.tight_layout()
# plt.savefig('../figures/GRDC_stations_wflow_gridpoints.png', dpi=300)
# plt.show()


# %% ##############################################################
############## COMPARE WFLOW VS. GRDC TIME SERIES DIRECTLY ########
###################################################################

# print("Comparing WFLOW time series to GRDC records...")
# fig, axs = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

# for i, river in enumerate(wanted_rivers):
#     idx_grdc = [i for i, n in enumerate(names) if n == river][0]
#     color = ["#0C75C0", "#4DA54E"]
#     ts = station_series[i]

#     axs[i].plot(ts.time, ts.values, label=f'wflow Q{i}', color=color[i])
#     axs[i].plot(GRDC_data['time'], GRDC_data['runoff_mean'].isel(id=idx[matches_df['station_index'].iloc[i]]), label=f'GRDC {river} (G{i})', color='orange')
    
#     axs[i].axhline(GRDC_data['runoff_mean'].isel(id=idx[matches_df['station_index'].iloc[i]]).max(), linestyle='--', color='orange', label='GRDC max')
#     axs[i].axhline(ts.values.max(), linestyle='--', color=color[i], label='wflow max')

#     axs[i].set_ylabel('Discharge [m3/s]')
#     axs[i].set_title(f"{river} and wflow Q{i+1}")
#     axs[i].legend()

# axs[1].set_xlabel('Time')
# plt.tight_layout()
# plt.savefig('../figures/GRDC_stations_wflow_timeseries_comparison.png', dpi=300)
# plt.show()


# %% ##############################################################
############## COMPARE WFLOW VS. GRDC STATISTICALLY ###############
###################################################################
# print("Comparing WFLOW time series to GRDC records with skill statistics...")

# fig, axs = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
# color = ["#0C75C0", "#4DA54E"]  # Blue for BUDZI, Green for RIOPUNGOE

# for i, river in enumerate(wanted_rivers):
#     print(f"Processing {river}...")

#     idx_grdc = [i for i, n in enumerate(names) if n == river][0]
#     grdc_ts = GRDC_data['runoff_mean'].isel(id=idx[matches_df['station_index'].iloc[i]])
#     wflow_ts = station_series[i]

#     # Stats
#     print(grdc_ts)
#     print(wflow_ts)
#     common_time = np.intersect1d(wflow_ts.time.values, grdc_ts.time.values)
#     wflow_aligned = wflow_ts.sel(time=common_time)
#     grdc_aligned = grdc_ts.sel(time=common_time)

#     wflow_stat = grdc_aligned.chunk(dict(time=-1))
#     grdc_stat = wflow_aligned.chunk(dict(time=-1))

#     kge = skillstats.kge(wflow_stat, grdc_stat, dim='time')["kge"].values.round(2)
#     # nse = skillstats.nashsutcliffe(wflow_stat, grdc_stat, dim='month').values.round(2)

#     # Plot time series
#     axs[i].plot(wflow_ts.time, wflow_ts.values, label=f'WFLOW Q{i},  KGE: {kge}', color=color[i])
#     axs[i].plot(grdc_ts.time, grdc_ts.values, label=f'GRDC {river} (G{i})', color='orange')

#     # Add max lines
#     axs[i].axhline(grdc_ts.max().values, linestyle='--', color='orange', label='GRDC max')
#     axs[i].axhline(wflow_ts.values.max(), linestyle='--', color=color[i], label='WFLOW max')

#     axs[i].set_ylabel('Discharge [m³/s]')
#     axs[i].set_title(f"{river} and WFLOW Q{i}")
#     axs[i].legend()
#     axs[i].grid(True)

# axs[1].set_xlabel('Time')
# plt.tight_layout()
# plt.savefig('../figures/GRDC_stations_wflow_timeseries_comparison_with_stats.png', dpi=300)
# plt.show()

# print("Daily time series comparison with statistics saved.")


#%%
# plot for 1974 & 1975
# print("Comparing WFLOW time series to GRDC records with skill statistics...")

# fig, axs = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
# color = ["#0C75C0", "#4DA54E"]  # Blue for BUDZI, Green for RIOPUNGOE

# # Define time range
# time_range = slice("1974-01-01", "1975-12-31")

# for i, river in enumerate(wanted_rivers):
#     print(f"Processing {river}...")

#     idx_grdc = [i for i, n in enumerate(names) if n == river][0]
#     grdc_ts = GRDC_data['runoff_mean'].isel(id=idx[matches_df['station_index'].iloc[i]])
#     wflow_ts = station_series[i]

#     # Filter both time series to 1974–1975
#     grdc_ts = grdc_ts.sel(time=time_range)
#     wflow_ts = wflow_ts.sel(time=time_range)

#     # Align times for fair comparison
#     common_time = np.intersect1d(wflow_ts.time.values, grdc_ts.time.values)
#     wflow_aligned = wflow_ts.sel(time=common_time)
#     grdc_aligned = grdc_ts.sel(time=common_time)

#     # Stats
#     wflow_stat = grdc_aligned.chunk(dict(time=-1))
#     grdc_stat = wflow_aligned.chunk(dict(time=-1))

#     kge = skillstats.kge(wflow_stat, grdc_stat, dim='time')["kge"].values.round(2)
#     # nse = skillstats.nashsutcliffe(wflow_stat, grdc_stat, dim='time').values.round(2)

#     # Plot time series
#     axs[i].plot(wflow_ts.time, wflow_ts.values, label=f'WFLOW Q{i},  KGE: {kge}', color=color[i])
#     axs[i].plot(grdc_ts.time, grdc_ts.values, label=f'GRDC {river} (G{i})', color='orange')

#     # Add max lines
#     axs[i].axhline(grdc_ts.max().values, linestyle='--', color='orange', label='GRDC max')
#     axs[i].axhline(wflow_ts.values.max(), linestyle='--', color=color[i], label='WFLOW max')

#     axs[i].set_ylabel('Discharge [m³/s]')
#     axs[i].set_title(f"{river} and WFLOW Q{i} (1974–1975)")
#     axs[i].legend()
#     axs[i].grid(True)

# axs[1].set_xlabel('Time')
# plt.tight_layout()
# plt.savefig('../figures/GRDC_stations_wflow_timeseries_comparison_1974_1975.png', dpi=300)
# plt.show()

# print("Daily time series comparison with statistics saved.")



#%%
print("Creating monthly climatology comparison between WFLOW and GRDC with uncertainty bands and statistics...")

fig, axs = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
color = ["#0C75C0", "#4DA54E"]  # Blue for BUDZI, Green for RIOPUNGOE

for i, river in enumerate(wanted_rivers):
    print(f"Processing {river}...")

    # Extract time series
    grdc_ts = GRDC_data['runoff_mean'].isel(id=idx[matches_df['station_index'].iloc[i]])
    wflow_ts = station_series[i]
    wflow_ts_1989 = station_series_1989[i]
    wfow_maxL06 = station_series_max[i]

    # Group by month
    grdc_monthly_mean = grdc_ts.groupby(grdc_ts['time.month']).mean()
    grdc_monthly_std = grdc_ts.groupby(grdc_ts['time.month']).std()

    wflow_monthly_mean = wflow_ts.groupby(wflow_ts['time.month']).mean()
    wflow_monthly_std = wflow_ts.groupby(wflow_ts['time.month']).std()

    wflow_1989_monthly_mean = wflow_ts_1989.groupby(wflow_ts_1989['time.month']).mean()
    wflow_1989_monthly_std = wflow_ts_1989.groupby(wflow_ts_1989['time.month']).std()

    wfow_maxL06_monthly_mean = wfow_maxL06.groupby(wfow_maxL06['time.month']).mean()
    wfow_maxL06_monthly_std = wfow_maxL06.groupby(wfow_maxL06['time.month']).std()

    months = np.arange(1, 13)

    # Stats
    wflow_kge = wflow_monthly_mean.drop_vars(['lat', 'lon'], errors='ignore').reset_coords(drop=True).chunk(dict(month=-1))
    wflow_1989_kge = wflow_1989_monthly_mean.drop_vars(['lat', 'lon'], errors='ignore').reset_coords(drop=True).chunk(dict(month=-1))
    wflow_max_kge = wfow_maxL06_monthly_mean.drop_vars(['lat', 'lon'], errors='ignore').reset_coords(drop=True).chunk(dict(month=-1))
    grdc_kge = grdc_monthly_mean.drop_vars(['id'], errors='ignore').reset_coords(drop=True).chunk(dict(month=-1))

    print(wflow_kge)
    print(grdc_kge)
    kge = skillstats.kge(wflow_kge, grdc_kge, dim='month')["kge"].values.round(2)
    kge_1989 = skillstats.kge(wflow_1989_kge, grdc_kge, dim='month')["kge"].values.round(2)
    kge_max = skillstats.kge(wflow_max_kge, grdc_kge, dim='month')["kge"].values.round(2)
    # nse = skillstats.nashsutcliffe(wflow_kge, grdc_kge, dim='month').values.round(2)

    # Plot monthly means
    axs[i].plot(months, grdc_monthly_mean.values, label=f"GRDC {river} 1954-1984 (G{i})", color='orange')
    axs[i].plot(months, wflow_monthly_mean.values, label=f"WFLOW Q{i} 1954-1984, KGE: {kge}", color=color[i])
    axs[i].plot(months, wflow_1989_monthly_mean.values, label=f"WFLOW Q{i} 1989-2019, KGE: {kge_1989}", color='#9E4DA5')
    axs[i].plot(months, wfow_maxL06_monthly_mean.values, label=f"WFLOW Q{i} 1989-2019 ML 0.6, KGE: {kge_max}", color='#9E4DA5')

    # Uncertainty bands (mean ± std)
    axs[i].fill_between(months,
                        (grdc_monthly_mean - grdc_monthly_std).clip(min=0),
                        grdc_monthly_mean + grdc_monthly_std,
                        color='orange', alpha=0.3)

    axs[i].fill_between(months,
                        (wflow_monthly_mean - wflow_monthly_std).clip(min=0),
                        wflow_monthly_mean + wflow_monthly_std,
                        color=color[i], alpha=0.3)
    
    axs[i].fill_between(months,
                        (wflow_1989_monthly_mean - wflow_1989_monthly_std).clip(min=0),
                        wflow_1989_monthly_mean + wflow_1989_monthly_std,
                        color="#9E4DA5", alpha=0.3)

    axs[i].set_ylabel('Monthly Mean Discharge [m³/s]')
    axs[i].set_title(f"Monthly Climatology: {river}")
    axs[i].set_xticks(months)
    axs[i].set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                            'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    axs[i].legend()
    axs[i].grid(True)

axs[1].set_xlabel('Month')

plt.tight_layout()
plt.savefig('../figures/monthly_climatology_GRDC_vs_WFLOW_with_uncertainty_stats_upd.png', dpi=300)
plt.show()

print("Monthly climatology with uncertainty bands and Hydromt stats saved.")



#%%


##########################################################
############### PLOTTING SFINCS RESULTS ##################
##########################################################
# %%
# def make_sfincs_model(subfolder: str):
#     root = BASE_SFINCS / subfolder
#     return SfincsModel(root=str(root), mode="r", data_libs=datacat)

# datacat = [
#         '../../Workflows/03_data_catalogs/datacatalog_general.yml',
#         '../../Workflows/03_data_catalogs/datacatalog_SFINCS_obspoints.yml',
#         '../../Workflows/03_data_catalogs/datacatalog_SFINCS_coastal_coupling.yml',
#         '../../Workflows/03_data_catalogs/datacatalog_CF_forcing.yml'
#         ]

# data_catalog = DataCatalog(data_libs = datacat)

# #%%
# # Load in the SFINCS models
# mod_old        = make_sfincs_model("event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_base")
# mod_wflow_all  = make_sfincs_model("event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_updatedwflow_all")
# mod_noFpln     = make_sfincs_model("event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_updatedwflow_nofloodplains")
# mod_ERA5_daily = make_sfincs_model("event_tp_era5_daily_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0")
# mod_chirps     = make_sfincs_model("event_tp_chirps_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0")


# #%%
# # we set a threshold to mask minimum flood depth
# hmin = 0.05

# for mod, model_name in [
#     (mod_old, "old"),
#     (mod_wflow_all, "wflow all"),
#     (mod_noFpln, "no floodplains"),
#     (mod_ERA5_daily, "ERA5_daily"),
#     (mod_chirps, "CHIRPS daily"),
# ]:
#     # compute the maximum over all time steps
#     da_zsmax = mod.results["zsmax"].max(dim="timemax")
        
#     # downscale the floodmap
#     depfile  = str(BASE_SFINCS / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_base" / "subgrid" / "dep_subgrid.tif")
#     da_dep   = mod.data_catalog.get_rasterdataset(depfile)

#     da_hmax = utils.downscale_floodmap(
#         zsmax=da_zsmax,
#         dep=da_dep,
#         hmin=hmin,
#         # floodmap_fn=join(sfincs_root, "gis/floodmap.tif") # uncomment to save floodmap to <mod.root>/floodmap.tif
#         )
        
#     # GSWO dataset to mask permanent water by first geprojecting it to the subgrid of hmax
#     sfincs_region = mod.region
#     gwso_region   = data_catalog.get_rasterdataset("gswo", geom=sfincs_region, buffer=1000)
#     gswo_mask     = gwso_region.raster.reproject_like(da_hmax, method="max")
#     # permanent water where water occurence > 5%
#     da_hmax_masked = da_hmax.where(gswo_mask <= 5)

#     # Add the name attribute for identification
#     mod.results['hmax'] = da_hmax
#     mod.results['hmax_masked'] = da_hmax_masked

#     del da_hmax, da_zsmax, da_dep, gswo_mask  # Clean up to free memory


# #%%
# # Plot the actual flood maps
# projection = mod.crs.to_epsg()

# fig, ax = plt.subplots(nrows=3,
#                        ncols=2,
#                        figsize=(8,8),
#                        dpi=300,
#                        constrained_layout=True,
#                        subplot_kw={"projection": ccrs.epsg(projection)})

# fig.suptitle("Factual Max Flood Depth", fontsize=11)

# # Plot hmax masked
# im = mod_old.results['hmax_masked'].plot.pcolormesh(
#      ax=ax[0,0], cmap="Blues", vmin=0, vmax=5.0, 
#      add_colorbar=False)

# im = mod_wflow_all.results['hmax_masked'].plot.pcolormesh(
#      ax=ax[0,1], cmap="Blues", vmin=0, vmax=5.0, 
#      add_colorbar=False)

# im = mod_noFpln.results['hmax_masked'].plot.pcolormesh(
#      ax=ax[1,0], cmap="Blues", vmin=0, vmax=5.0, 
#      add_colorbar=False)

# im = mod_ERA5_daily.results['hmax_masked'].plot.pcolormesh(
#      ax=ax[1,1], cmap="Blues", vmin=0, vmax=5.0, 
#      add_colorbar=False)

# im = mod_chirps.results['hmax_masked'].plot.pcolormesh(
#      ax=ax[2,0], cmap="Blues", vmin=0, vmax=5.0, 
#      add_colorbar=False)

# # Add basemap and gridlines for each subplot
# for i, axis in enumerate(ax.flat):
#     ctx.add_basemap(
#         axis,
#         source=ctx.providers.Esri.WorldImagery,
#         zoom=10,
#         crs=mod_old.results['hmax_masked'].rio.crs,
#         attribution=False
#     )
#     gl = axis.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
#     gl.top_labels = False
#     gl.right_labels = False
#     gl.xlabel_style = {'size': 8}
#     gl.ylabel_style = {'size': 8}
#     gl.xformatter = LONGITUDE_FORMATTER
#     gl.yformatter = LATITUDE_FORMATTER
#     gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, .2))
#     gl.ylocator = mticker.FixedLocator(np.arange(-90, 90, .2))

# # Titles
# ax[0,0].set_title(f"Old wflow settings", fontsize=12)
# ax[0,1].set_title(f"New wflow settings", fontsize=12)
# ax[1,0].set_title(f"No 1D floodplain", fontsize=12)
# ax[1,1].set_title(f"ERA5 daily", fontsize=12)
# ax[2,0].set_title(f"CHIRPS daily", fontsize=12)

# # Add shared colorbar with larger size and labels
# cbar = fig.colorbar(im, ax=ax, orientation="vertical", 
#                     fraction=0.035, pad=0.05)
# cbar.set_label('Flood depth (m)', rotation=270, labelpad=10, fontsize=9)
# cbar.ax.tick_params(labelsize=8)
    
# fig.savefig("../figures/wflow_test_hmax_masked.png", bbox_inches='tight', dpi=300)
# plt.show()



# # %%
# # Calculate the difference between the two hmax_masked maps (mod_new - mod_old)
# hmax_diff_all = mod_wflow_all.results['hmax_masked'] - mod_old.results['hmax_masked']
# hmax_diff_noFpln = mod_noFpln.results['hmax_masked'] - mod_old.results['hmax_masked']
# hmax_diff_ERA5_daily = mod_ERA5_daily.results['hmax_masked'] - mod_old.results['hmax_masked']
# hmax_diff_chirps = mod_chirps.results['hmax_masked'] - mod_old.results['hmax_masked']

# # Plot the difference map as a new subplot
# fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8, 8), dpi=300,
#                        constrained_layout=True,
#                        subplot_kw={"projection": ccrs.epsg(projection)})

# fig.suptitle("Difference in Flood Depth", fontsize=11)

# # Plot the difference
# im0 = hmax_diff_all.plot.pcolormesh(
#     ax=ax[0,0], cmap="RdBu", vmin=-2, vmax=2, add_colorbar=False,
#     transform=ccrs.PlateCarree())

# im1 = hmax_diff_noFpln.plot.pcolormesh(
#     ax=ax[0,1], cmap="RdBu", vmin=-2, vmax=2, add_colorbar=False,
#     transform=ccrs.PlateCarree())

# im2 = hmax_diff_ERA5_daily.plot.pcolormesh(
#     ax=ax[1,0], cmap="RdBu", vmin=-2, vmax=2, add_colorbar=False,
#     transform=ccrs.PlateCarree())

# im3 = hmax_diff_chirps.plot.pcolormesh(
#     ax=ax[1,1], cmap="RdBu", vmin=-2, vmax=2, add_colorbar=False,
#     transform=ccrs.PlateCarree())

# ax[0,0].set_title("Difference (All settings - Old)", fontsize=12)
# ax[0,1].set_title("Difference (No floodplains - Old)", fontsize=12)
# ax[1,0].set_title("Difference (ERA5 daily - Old)", fontsize=12)
# ax[1,1].set_title("Difference (CHIRPS - Old)", fontsize=12)


# # Add basemap and gridlines for each subplot
# for i, axis in enumerate(ax.flat):
#     ctx.add_basemap(
#         axis,
#         source=ctx.providers.Esri.WorldImagery,
#         zoom=10,
#         crs=mod_old.results['hmax_masked'].rio.crs,
#         attribution=False
#     )
#     gl = axis.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
#     gl.top_labels = False
#     gl.right_labels = False
#     gl.xlabel_style = {'size': 8}
#     gl.ylabel_style = {'size': 8}
#     gl.xformatter = LONGITUDE_FORMATTER
#     gl.yformatter = LATITUDE_FORMATTER
#     gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, .2))
#     gl.ylocator = mticker.FixedLocator(np.arange(-90, 90, .2))

# # Add shared colorbar
# cbar = fig.colorbar(im3, ax=ax, orientation="vertical", fraction=0.04, pad=0.05)
# cbar.set_label('Difference in Flood depth (m)', rotation=270, labelpad=10, fontsize=9)
# cbar.ax.tick_params(labelsize=8)

# fig.savefig("../figures/wflow_test_hmax_diff.png", bbox_inches='tight', dpi=300)

# plt.show()

# # %%

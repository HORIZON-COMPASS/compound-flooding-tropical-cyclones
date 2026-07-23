
#%% script from Natalia Aleksandrove and adapted by Doris Vertegaal to track TCs from ClimateDT data
# run with pixi env compass-snake-dfm
from datetime import datetime
import platform
from matplotlib import cm
from matplotlib.lines import Line2D
import contextily as ctx
import matplotlib.lines as mlines
import numpy as np
import xarray as xr
import os
import geopandas as gpd
import re
from rasterio.features import rasterize
import rasterio as rio
import hydromt
import pandas as pd
from scipy.ndimage import minimum_filter
import matplotlib.pyplot as plt
import dfm_tools as dfmt
import cartopy.crs as ccrs


#%%

prefix = "p:/" if platform.system() == "Windows" else "/p/"

def preprocess_add_r(ds):
    source = ds.encoding.get('source', '')
    match = re.search(r'_r(\d+)_', source)
    # realization = int(match.group(1)) if match else -1
    # ds = ds.assign_coords(realization=realization)
    return ds

def get_mask_from_gdf(gdf, ds):
    new_gdf = gpd.GeoDataFrame(geometry=[gdf.unary_union], crs=gdf.crs)
    tr = rio.transform.from_bounds(west=ds.longitude[0], south=ds.latitude[0], east=ds.longitude[-1], north=ds.latitude[-1], 
                                   width=len(ds.longitude), height=len(ds.latitude))
    mask = rasterize(new_gdf.geometry, out_shape=(len(ds.latitude), len(ds.longitude)), transform=tr)
    mask_da = xr.DataArray(mask[::-1, :], coords=[ds.latitude, ds.longitude], dims=['latitude', 'longitude'])

    return mask_da

# Open the DFM his files in the output folder
def open_ds_his(dir, model):
    for fname in os.listdir(os.path.join(dir,model,'output')):
        if fname.endswith('_his.nc'):
            print(fname)
            file_nc_his = os.path.join(dir,model,'output',fname)

    #open hisfile with xarray and print netcdf structure
    if file_nc_his is not None:
        ds = xr.open_mfdataset(file_nc_his, preprocess=dfmt.preprocess_hisnc)

    ds['windmag'] = np.sqrt(ds['windx']**2 + ds['windy']**2)
    ds['windmag'].attrs['long_name'] = 'wind speed'
    ds['windmag'].attrs['units'] = 'm/s'

    return ds

# Open the DFM map files in the output folder
def open_ds_map(dir, model): 
    file_nc_map = []
    for fname in os.listdir(os.path.join(dir,model,'output')):
        if fname.endswith("map.nc"):
            print(fname)
            file_nc_map.append(os.path.join(dir,model,'output',fname))

    ds = dfmt.open_partitioned_dataset(file_nc_map)

    # compute magnitude of wind
    ds['mesh2d_windmag'] = np.sqrt(ds['mesh2d_windx']**2 + ds['mesh2d_windy']**2)

    return ds

# Plot wind from DFM spatially for a time slice
def plot_wind_slice(ds, time_str, var='mesh2d_windmag', title=None):
    """Plot wind magnitude at a specific time from a UGRID dataset."""
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10, 6))

    da = ds[var].sel(time=time_str)
    da.ugrid.plot(ax=ax, transform=ccrs.PlateCarree(), cmap='viridis')

    ax.coastlines()
    ax.set_title(title or f"Wind magnitude at {time_str}")
    return ax

# Find local minima in the mean sea level pressure (MSLP) data
def find_mslp_minima(mslp):
    return mslp == minimum_filter(mslp, size=3, mode="nearest")

# Calculate a fast distance in kilometers between two lat/lon points
def fast_distance_km(lat1, lon1, lat2, lon2):
    return 111 * np.sqrt((lat1-lat2)**2 + (lon1-lon2)**2)

#%%
# Data directory
datadir = os.path.join(prefix, '11210471-001-compass/01_Data/ECMWF_ClimateDT/data/preprocessed/Gen2/Idai_full')

# Data catalog
data_cat = ['../../../03_data_catalogs/datacatalog_general.yml', '../../../03_data_catalogs/datacatalog_CF_forcing.yml'] 
data_catalog = hydromt.data_catalog.DataCatalog(data_libs = data_cat)

# Spatial and temporal domain
bbox_wflow_sfincs = [32.375,-20.941667,35.35,-17.525]
# bbox_idai = [30, -28, 45, -9] 
bbox_idai = [30, -22, 45, -12] 
areaname = 'wflow+sfincs'

time_min = "2019-03-1"
time_max = "2019-03-25"

# Load model shapefiles for wflow and sfincs
gdf_sfincs = gpd.read_file(os.path.join(prefix, "11210471-001-compass/02_Models/sofala/Idai/sfincs/gis/region.geojson"))
gdf_sfincs = gdf_sfincs.to_crs(epsg=4326) 

gdf_wflow = gpd.read_file(os.path.join(prefix, "11210471-001-compass/02_Models/sofala/Idai/wflow/staticgeoms/basins.geojson"))
gdf_wflow = gdf_wflow.to_crs(epsg=4326) 

# DFM data
dir_runs = os.path.join(prefix, '11210471-001-compass/03_Runs/sofala/Idai/dfm')

#%%
# Load TC Idai track from shapefile and convert to list of (lat, lon) tuples
data_base = os.path.join(prefix, "11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS")

shapefile_path = os.path.join(data_base, "IBTrACS_IDAI.shp")
gdf = gpd.read_file(shapefile_path)
tc_idai = gdf[gdf['SID'] == '2019063S18038']

tc_idai_points = list(zip(tc_idai.geometry.y, tc_idai.geometry.x))

# Point from where Idai was a TC1
first_tc1_idx = tc_idai[tc_idai['USA_SSHS'] == 1].index[0]
tc_idai_tc1 = tc_idai.loc[first_tc1_idx:]
tc_idai_points_TC = list(zip(tc_idai_tc1.geometry.y, tc_idai_tc1.geometry.x))


#%%
##### DestinE ClimateDT storyline data #####
# Reading in the data for the control and historical runs, and selecting the spatial domain of interest
ds_cont = xr.open_mfdataset(os.path.join(datadir, 'climateDT_*_cont_*_Idai_full_20190301_to_20190329.nc'),
                            preprocess=preprocess_add_r)

ds_hist = xr.open_mfdataset(os.path.join(datadir, 'climateDT_*_hist_*_Idai_full_20190301_to_20190329.nc'),
                            preprocess=preprocess_add_r)

ds_tp2k = xr.open_mfdataset(os.path.join(datadir, 'climateDT_*_Tplus2.0K_*_Idai_full_20190301_to_20190329.nc'),
                            preprocess=preprocess_add_r)

# Reading best TC like Idai per realisationt racked using TC_tracking_Idai.py
TC_like_Idai_climDT = gpd.read_file('data/best_track_per_realization.gpkg')

#%%
###### WIND ######
##### Observational data #######
# Reading in ERA5 data for the same time period and spatial domain
ds_era5_large = data_catalog.get_rasterdataset(
    data_like = "era5_hourly_zarr",
    time_tuple = ("2019-03-01","2019-03-29"),
    variables = ['precip', 'wind10_u', 'wind10_v', 'press_msl'],
    bbox=bbox_idai, 
    buffer = 1
)

ds_era5 = ds_era5_large.sel(latitude=slice(bbox_wflow_sfincs[3]+0.25, bbox_wflow_sfincs[1]-0.25),
                            longitude=slice(bbox_wflow_sfincs[0]-0.25, bbox_wflow_sfincs[2]+0.25))
ds_tp_era5 = ds_era5['precip'].sortby('latitude', ascending=True)
ds_tp_era5 = ds_tp_era5 / 1000 # as the data catalog converts it from m to mm. making it explcitit here for clarity.
ds_tp_era5.attrs['units'] = 'm'

ds_era5_wind = ds_era5[['wind10_u','wind10_v','press_msl']].sortby('latitude', ascending=True)
ds_era5_wind_large = ds_era5_large[['wind10_u','wind10_v','press_msl']].sortby('latitude', ascending=True)

#%%
# track Idai in ERA5 and extend until end of intense precipitation
RADIUS_KM = 250
era5_track = []
# -------------------------
# Follow observed Idai track
# -------------------------
for _, row in tc_idai.iterrows():
    t = np.datetime64(row["ISO_TIME"])
    ds_t = ds_era5_wind_large.sel(time=t, method="nearest")
    mslp = ds_t["press_msl"].values

    print(f"Tracking TC Idai at time {t} and finding mslp minima")
    minima = find_mslp_minima(mslp)
    yy, xx = np.where(minima)

    best = None
    best_score = None
    for i, j in zip(yy, xx):
        lat_c = float(ds_t.latitude.values[i])
        lon_c = float(ds_t.longitude.values[j])
        dist = fast_distance_km(row["LAT"], row["LON"], lat_c, lon_c)

        if dist > RADIUS_KM:
            continue

        score = (dist, mslp[i, j])
        if best_score is None or score < best_score:
            best_score = score
            best = (i, j)

    if best is None:
        continue

    i, j = best
    era5_track.append({"time": t, "lat": float(ds_t.latitude.values[i]), 
                       "lon": float(ds_t.longitude.values[j]), "mslp": float(mslp[i, j]),})
#%%
# -----------------------
# DAILY SUBPLOTS OF PRECIPITATION AND MSLP AROUND LANDFALL AND AFTER
# -----------------------
days = pd.date_range("2019-03-13", "2019-03-21", freq="D")

fig, axes = plt.subplots(3, 3, figsize=(15, 10), constrained_layout=True)

axes = axes.ravel()
for ax, day in zip(axes, days):
    # Daily precipitation (mm) over wflow domain
    ds_day = ds_era5_large.sel(time=slice(day, day + pd.Timedelta(days=1)))
    precip = ds_day["precip"].sum("time")
    precip_max = float(max(ds_day["precip"].sum("time").max() for d in days))

    # Mean sea-level pressure (hPa) over Idai domain
    mslp = ds_day["press_msl"].mean("time")

    # Precipitation shading
    precip.plot(ax=ax, cmap="Blues", vmin=0, vmax=precip_max, alpha=0.7, add_colorbar=False, zorder=1)

    # Pressure contours
    levels = np.arange(np.floor(float(mslp.min())), np.ceil(float(mslp.max())) + 2, 2)
    cs = ax.contour(ds_era5_large.longitude, ds_era5_large.latitude, mslp, levels=levels, colors="yellow", linewidths=1.5, zorder=10)
    ax.clabel(cs, fmt="%d", fontsize=7)

    # Optional: ERA5 track
    if "era5_track" in locals():
        ax.plot([p["lon"] for p in era5_track], [p["lat"] for p in era5_track], "ORANGE", lw=2, zorder=5)

    # Wflow boundary
    gdf_wflow.boundary.plot(ax=ax, color="lightblue", linewidth=2, zorder=4)

    # Set to larger extent    
    ax.set_xlim(ds_era5_large.longitude.min(), 38)
    ax.set_ylim(ds_era5_large.latitude.min(), -16)

    # Basemap
    ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=7, crs=gdf_wflow.crs,       attribution=False, zorder=0, rasterized=True)

    ax.set_title(day.strftime("%d %b %Y"))
    ax.set_xlabel("")
    ax.set_ylabel("")

# Shared colorbar
fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, precip_max), cmap="Blues"),
             ax=axes, shrink=0.8, label="Daily precipitation (mm)")
plt.suptitle("ERA5 Daily Precipitation (shading) and daily mean MSLP (contours)\n13-21 March 2019",
    fontsize=18)

plt.show()




#%%
# -----------------------
# DAILY MSLP ONLY
# -----------------------
days = pd.date_range("2019-03-04", "2019-03-21", freq="D")

fig, axes = plt.subplots(3, 6, figsize=(24, 12), constrained_layout=True)

axes = axes.ravel()
for ax, day in zip(axes, days):
    # 12 UTC snapshot
    mslp = (ds_era5_large.sel(time=f"{day:%Y-%m-%d}T12:00", method="nearest")["press_msl"]) 

    # Filled pressure field
    pcm = ax.contourf(ds_era5_large.longitude, ds_era5_large.latitude, mslp.values,
                      levels=np.arange(950, 1026, 2), cmap="viridis", alpha=0.6, zorder=1)

    # Pressure contours
    cs = ax.contour(ds_era5_large.longitude, ds_era5_large.latitude, mslp.values, 
                    levels=np.arange(950, 1026, 4), colors="yellow", linewidths=1.5, zorder=10)
    ax.clabel(cs, fmt="%d", fontsize=7)

    # ERA5-tracked Idai
    if len(era5_track) > 0:
        track_day = [p for p in era5_track if pd.Timestamp(p["time"]) <= day + pd.Timedelta(days=1)]

        if len(track_day) > 0:
            ax.plot([p["lon"] for p in track_day], [p["lat"] for p in track_day], color="red", 
                    lw=2, zorder=20)
            ax.scatter(track_day[-1]["lon"], track_day[-1]["lat"], c="red", s=40, 
                       edgecolors="white", zorder=21)
    # Wflow boundary
    gdf_wflow.boundary.plot(ax=ax, color="cyan", linewidth=2, zorder=30)

    ax.set_xlim(bbox_idai[0], bbox_idai[2])
    ax.set_ylim(bbox_idai[1], bbox_idai[3])

    ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=7, crs=gdf_wflow.crs, 
                    attribution=False, zorder=0, rasterized=True)

    ax.set_title(day.strftime("%d %b %Y"))
    ax.set_xlabel("")
    ax.set_ylabel("")

# remove unused panels
for ax in axes[len(days):]:
    ax.remove()

fig.colorbar(pcm, ax=axes.tolist(), shrink=0.8, label="MSLP (hPa)")
plt.suptitle("ERA5 Mean Sea-Level Pressure and Tracked Idai Position\n4-21 March 2019", fontsize=18)

plt.show()

#%%
# # -------------------------
# # Extend to 21 March
# # -------------------------
# end_time = np.datetime64("2019-03-21T23:00")
# current_lat = era5_track[-1]["lat"]
# current_lon = era5_track[-1]["lon"]
# all_times = ds_era5_wind.time.values
# future_times = all_times[all_times > era5_track[-1]["time"]]

# for t in future_times:
#     if t > end_time:
#         break

#     ds_t = ds_era5_wind.sel(time=t)
#     mslp = ds_t["press_msl"].values
#     minima = find_mslp_minima(mslp)
#     yy, xx = np.where(minima)

#     best = None
#     best_score = None
#     for i, j in zip(yy, xx):
#         lat_c = float(ds_t.latitude.values[i])
#         lon_c = float(ds_t.longitude.values[j])
#         dist = fast_distance_km(current_lat, current_lon, lat_c, lon_c)

#         # maximum translation speed constraint
#         if dist > 50:
#             continue

#         score = (mslp[i, j], dist)
#         if best_score is None or score < best_score:
#             best_score = score
#             best = (i, j)

#     if best is None:
#         break

#     i, j = best
#     current_lat = float(ds_t.latitude.values[i])
#     current_lon = float(ds_t.longitude.values[j])

#     era5_track.append({"time": t, "lat": current_lat, "lon": current_lon, "mslp": float(mslp[i, j]),})

# # %%
# # Get masks for wflow region, and Pungwe and Buzi basins 
# mask_da = get_mask_from_gdf(gdf_wflow, ds_cont)
# mask_da_era5 = get_mask_from_gdf(gdf_wflow, ds_tp_era5)

# gdf_Pungwe = gdf_wflow[gdf_wflow['geometry'].area>2][gdf_wflow['value']==42]
# gdf_Buzi = gdf_wflow[gdf_wflow['geometry'].area>2][gdf_wflow['value']==67]

# mask_da_pungwe = get_mask_from_gdf(gdf_Pungwe, ds_cont)
# mask_da_buzi = get_mask_from_gdf(gdf_Buzi, ds_cont)

# %%
# IBTrACS wind speed data merged with ERA5 in the background, from DFM
ds_dfm_SLR_0cm_wind_0_original  = open_ds_his(dir_runs, 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0')
wm_dfm_SLR_0cm_wind_0 = ds_dfm_SLR_0cm_wind_0_original["windmag"].sel(station='BEIRA IHO').compute()
pm_dfm_SLR_0cm_wind_0 = ds_dfm_SLR_0cm_wind_0_original["patm"].sel(station='BEIRA IHO').compute() / 100  # Convert from Pa to hPa


# Calculate wind speed magnitude from u and v components
ds_hist['umag'] = np.sqrt(ds_hist['10u']**2 + ds_hist['10v']**2)
ds_cont['umag'] = np.sqrt(ds_cont['10u']**2 + ds_cont['10v']**2)
ds_tp2k['umag'] = np.sqrt(ds_tp2k['10u']**2 + ds_tp2k['10v']**2)
ds_hist['umag'].attrs['units'] = 'm/s'
ds_cont['umag'].attrs['units'] = 'm/s'  
ds_tp2k['umag'].attrs['units'] = 'm/s'

# Same for ERA5 wind data
ds_era5_wind['umag'] = np.sqrt(ds_era5_wind['wind10_u']**2 + ds_era5_wind['wind10_v']**2)
ds_era5_wind['umag'].attrs['units'] = 'm/s'

# Correct pressure from Pa to hPa
ds_hist['msl'] = ds_hist['msl'] / 100
ds_cont['msl'] = ds_cont['msl'] / 100
ds_tp2k['msl'] = ds_tp2k['msl'] / 100

# open map files of the factual DFM run
# model_name = "event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0"
# ds_map_F = open_ds_map(dir_runs, model_name)

# model_name = "event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF-10"
# ds_map_CF10 = open_ds_map(dir_runs, model_name)

# #%%
# # Plot wind speed spatially
# plot_wind_slice(ds_map_F, time_str='2019-03-14 12:00:00')
# plot_wind_slice(ds_map_CF10, time_str='2019-03-14 12:00:00')

# # Actual plotting
# fig, ax  = plt.subplots(1, 1, figsize=(12,6))
# ax.plot(wm_dfm_SLR_0cm_wind_0["time"], wm_dfm_SLR_0cm_wind_0, label="wind SLR 0cm Wind 0")

# plt.xlabel("Time")
# plt.ylabel("Wind speed (m/s)")
# plt.legend()
# plt.tight_layout()
# plt.show()

#%%
coords = [wm_dfm_SLR_0cm_wind_0['station_x_coordinate'].values, 
          wm_dfm_SLR_0cm_wind_0['station_y_coordinate'].values]

fig, (ax0, ax1) = plt.subplots(figsize=(10, 7), nrows=2)
lims=[0, 1]
times=['2019-03-12','2019-03-17']

times = [datetime.strptime('2019-03-12', '%Y-%m-%d'),  datetime.strptime('2019-03-17', '%Y-%m-%d')]


for rr in ds_hist.realization.values:
    ds_hist.sel(latitude=coords[1], longitude=coords[0], method='nearest').sel(time=slice(times[0],times[1]), realization=rr).umag.plot(ax=ax0, linewidth=2, alpha=0.5, color='tab:blue')
    ds_cont.sel(latitude=coords[1], longitude=coords[0], method='nearest').sel(time=slice(times[0],times[1]), realization=rr).umag.plot(ax=ax0, linewidth=2, alpha=0.5, color='tab:orange')
    ds_tp2k.sel(latitude=coords[1], longitude=coords[0], method='nearest').sel(time=slice(times[0],times[1]), realization=rr).umag.plot(ax=ax0, linewidth=2, alpha=0.5, color='tab:green')
    
ds_hist.sel(latitude=coords[1], longitude=coords[0], method='nearest').sel(time=slice(times[0],times[1])).mean(dim='realization').umag.plot(ax=ax0, label='ClimateDT factual', linewidth=2, color='tab:blue')
ds_cont.sel(latitude=coords[1], longitude=coords[0], method='nearest').sel(time=slice(times[0],times[1])).mean(dim='realization').umag.plot(ax=ax0, label='ClimateDT counterfactual', linewidth=2, color='tab:orange')
ds_tp2k.sel(latitude=coords[1], longitude=coords[0], method='nearest').sel(time=slice(times[0],times[1])).mean(dim='realization').umag.plot(ax=ax0, label='ClimateDT TP2K', linewidth=2, color='tab:green')
ds_era5_wind.sel(latitude=coords[1], longitude=coords[0], method='nearest').sel(time=slice(times[0],times[1])).umag.plot(ax=ax0, label='ERA5 reanalysis', color='gray', linewidth=2)
wm_dfm_SLR_0cm_wind_0.sel(time=slice(times[0],times[1])).plot(ax=ax0, label='IBTrACS+ERA5', color='black', linewidth=1.5, linestyle='--')

for rr in ds_hist.realization.values:
    ds_hist.sel(latitude=coords[1], longitude=coords[0], method='nearest').sel(time=slice(times[0],times[1]), realization=rr).msl.plot(ax=ax1, linewidth=2, alpha=0.5, color='tab:blue')
    ds_cont.sel(latitude=coords[1], longitude=coords[0], method='nearest').sel(time=slice(times[0],times[1]), realization=rr).msl.plot(ax=ax1, linewidth=2, alpha=0.5, color='tab:orange')
    ds_tp2k.sel(latitude=coords[1], longitude=coords[0], method='nearest').sel(time=slice(times[0],times[1]), realization=rr).msl.plot(ax=ax1, linewidth=2, alpha=0.5, color='tab:green')

ds_hist.sel(latitude=coords[1], longitude=coords[0], method='nearest').sel(time=slice(times[0],times[1])).mean(dim='realization').msl.plot(ax=ax1, label='ClimateDT factual', linewidth=2, color='tab:blue')
ds_cont.sel(latitude=coords[1], longitude=coords[0], method='nearest').sel(time=slice(times[0],times[1])).mean(dim='realization').msl.plot(ax=ax1, label='ClimateDT counterfactual', linewidth=2, color='tab:orange')
ds_tp2k.sel(latitude=coords[1], longitude=coords[0], method='nearest').sel(time=slice(times[0],times[1])).mean(dim='realization').msl.plot(ax=ax1, label='ClimateDT TP2K', linewidth=2, color='tab:green')
ds_era5_wind.sel(latitude=coords[1], longitude=coords[0], method='nearest').sel(time=slice(times[0],times[1])).press_msl.plot(ax=ax1, label='ERA5 reanalysis', color='gray', linewidth=2)
pm_dfm_SLR_0cm_wind_0.sel(time=slice(times[0],times[1])).plot(ax=ax1, label='IBTrACS+ERA5', color='black', linewidth=1.5, linestyle='--')


for ax in [ax0, ax1]:
    ax.set_xlim(times[0], times[1])
    ax.grid()
    ax.set_title('')
    if ax == ax0:
        ax.set_ylabel('10-m wind magnitude [m/s]')
        ax.set_xlabel('')
        ax.legend()
    else:
        ax.set_ylabel('Mean sea level pressure [hPa]')
        ax.set_xlabel('Time (day)')
        # ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        

fig.suptitle(f'Comparing wind speed and pressure for BEIRA_IHO: lat {coords[1]}, lon {coords[0]}', fontweight='bold')

# %%
# Calculate the maximum wind speed and minimum pressure for each realization in the specified time range
hist_max_wind = [ds_hist.sel(latitude=coords[1], longitude=coords[0], method='nearest',
                             realization=r).sel(time=slice(*times)).umag.max().values
                             for r in ds_hist.realization.values]
cont_max_wind = [ds_cont.sel(latitude=coords[1], longitude=coords[0], method='nearest',
                             realization=r).sel(time=slice(*times)).umag.max().values
                             for r in ds_cont.realization.values]
tp2k_max_wind = [ds_tp2k.sel(latitude=coords[1], longitude=coords[0], method='nearest',
                             realization=r).sel(time=slice(*times)).umag.max().values
                             for r in ds_tp2k.realization.values]

hist_min_msl = [ds_hist.sel(latitude=coords[1], longitude=coords[0], method='nearest',
                            realization=r).sel(time=slice(*times)).msl.min().values
                            for r in ds_hist.realization.values]
cont_min_msl = [ds_cont.sel(latitude=coords[1], longitude=coords[0], method='nearest',
                            realization=r).sel(time=slice(*times)).msl.min().values
                            for r in ds_cont.realization.values]
tp2k_min_msl = [ds_tp2k.sel(latitude=coords[1], longitude=coords[0], method='nearest',
                            realization=r).sel(time=slice(*times)).msl.min().values
                            for r in ds_tp2k.realization.values]

#%%
# Create boxplots for maximum wind speed and minimum pressure
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(8, 4))

colors = ['tab:orange', 'tab:blue', 'tab:green']

# Peak wind
bp0 = ax0.boxplot([cont_max_wind, hist_max_wind, tp2k_max_wind],
                  labels=['CF_past', 'Factual', 'CF_future'], patch_artist=True)

# Minimum pressure
bp1 = ax1.boxplot([cont_min_msl, hist_min_msl, tp2k_min_msl],
                  labels=['CF_past', 'Factual', 'CF_future'], patch_artist=True)

for bp in [bp0, bp1]:
    for i, c in enumerate(colors):
        bp['boxes'][i].set_facecolor(c)
        bp['boxes'][i].set_alpha(0.4)
        bp['medians'][i].set_color(c)
        bp['medians'][i].set_linewidth(2)

ax0.set_ylabel('Maximum 10-m wind speed [m/s]')
ax1.set_ylabel('Minimum MSLP [Pa]')

ax0.set_title('Peak wind')
ax1.set_title('Minimum pressure')

ax0.grid(axis='y')
ax1.grid(axis='y')

plt.tight_layout()
# %%




# %%
data_base = "p:/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS"

shapefile_path = os.path.join(data_base, "IBTrACS_IDAI.shp")
gdf = gpd.read_file(shapefile_path)
tc_idai = gdf[gdf['SID'] == '2019063S18038']

# Normalize windspeed for point sizing (adjust scaling as needed)
tc_idai["size"] = (tc_idai["USA_WIND"] - tc_idai["USA_WIND"].min()) / (
    tc_idai["USA_WIND"].max() - tc_idai["USA_WIND"].min()) * 138 + 19  # Scale between 20 and 145

wind_min = tc_idai["USA_WIND"].min()
wind_max = tc_idai["USA_WIND"].max()

# create colormap and normalize function for colors (wind speed to 0-1)
cmap_idai = cm.get_cmap("Reds")
norm = plt.Normalize(wind_min, wind_max)



# %%
# Set up figure
fig, ax = plt.subplots(figsize=(9, 5))

# Plot SFINCS and wflow regions and set up legend entry
gdf_wflow.plot(ax=ax, edgecolor="gray", facecolor="none")
gdf_sfincs.plot(ax=ax, edgecolor="black", facecolor="none") 
    
# Plot TC Tracks
# tc_idai_filtered = tc_idai[tc_idai.geometry.y < -18]
tc_scatter = ax.scatter(tc_idai.geometry.x, tc_idai.geometry.y, s=tc_idai["size"]*0.5,  # Size based on wind speed
                        c=tc_idai["USA_WIND"], cmap=cmap_idai, alpha=0.7, label="Track TC Idai", zorder=6, rasterized=True)
ax.plot(tc_idai.geometry.x, tc_idai.geometry.y, color="grey", linewidth=1, alpha=0.7, linestyle="-", zorder=5, rasterized=True)
# Create a custom legend entry for the TC scatter plot (you can modify the color or markersize)
tc_marker = Line2D([0], [0], marker='o', color='darkgrey', markerfacecolor='red', markersize=5, label='Track TC Idai')

# Add basemap (LOWER zoom = faster)
ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=7, crs=tc_idai.crs, attribution=False, zorder=0, rasterized=True)

txt = ax.text(32.01, -21.99, "Sources: Esri, i-cubed, USDA, USGS, AEX, \nGeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, \nand the GIS User Community | Powered by Esri",
              fontsize=5.5, color='white', alpha=0.7, ha='left', va='bottom', zorder=20)

# Define ticks for x and y axes
xticks = np.arange(32, 42, 1)    
yticks = np.arange(-21, -14, 1) 
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_xticklabels([f"{x}°E" for x in xticks])
ax.set_yticklabels([f"{abs(y)}°S" for y in yticks])

# Set limits
extent = [30, -22, 45, -14.5]  # [lon_min, lat_min, lon_max, lat_max]
ax.set_xlim(extent[0], extent[2])
ax.set_ylim(extent[1], extent[3])

ax.annotate("Track TC Idai", xy=(36.02,-19.8), xytext=(37.2,-19.15), textcoords='data',
            arrowprops=dict(arrowstyle="->", color='darkred', lw=1.5),
            bbox=dict(boxstyle="round,pad=0.3", fc="#d3d3d3", ec="darkred", alpha=0.7),
            fontsize=10, color='darkred', fontweight='bold', zorder=10)

# Add Legend
size_labels = ['0-50 km/h', '50-100 km/h', '>100 km/h']
legend_colors = ['#ff9999', '#ff4d4d', '#b30000']  # light to dark red
size_handles = [
    mlines.Line2D([], [], marker='o', color='k', markerfacecolor=legend_colors[0], alpha=0.7, markersize=2.5, label=size_labels[0], linestyle=''),
    mlines.Line2D([], [], marker='o', color='k', markerfacecolor=legend_colors[1], alpha=0.7, markersize=3.75, label=size_labels[1], linestyle=''),
    mlines.Line2D([], [], marker='o', color='k', markerfacecolor=legend_colors[2], alpha=0.7, markersize=5, label=size_labels[2], linestyle='')
]


# %%

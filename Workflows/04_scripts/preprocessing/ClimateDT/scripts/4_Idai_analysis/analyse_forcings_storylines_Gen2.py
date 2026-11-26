#%% script from Natalia Aleksandrove and adapted by Doris Vertegaal to track TCs from ClimateDT data
# run with pixi env compass-snake-dfm
import platform
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import xarray as xr
import os
import geopandas as gpd
import re
import glob
from rasterio.features import rasterize
import rasterio as rio
import hydromt
import numpy as np
import pandas as pd
import geopandas as gpd
import os
import pickle
from scipy.ndimage import minimum_filter
from shapely.geometry import Point
from shapely.geometry import LineString
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import dfm_tools as dfmt
import os
import cartopy.crs as ccrs
import ast
import hydromt

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
gdf_sfincs = gpd.read_file(r"p:\11210471-001-compass\02_Models\sofala\Idai\sfincs\gis\region.geojson")
gdf_sfincs = gdf_sfincs.to_crs(epsg=4326) 

gdf_wflow = gpd.read_file(r"p:\11210471-001-compass\02_Models\sofala\Idai\wflow\staticgeoms\basins.geojson")
gdf_wflow = gdf_wflow.to_crs(epsg=4326) 

# DFM data
dir_runs = os.path.join(prefix, '11210471-001-compass/03_Runs/sofala/Idai/dfm')

#%%
##### DestinE ClimateDT storyline data #####
# Reading in the data for the control and historical runs, and selecting the spatial domain of interest
ds_cont = xr.open_mfdataset(os.path.join(datadir, 'climateDT_*_cont_*_Idai_full_20190301_to_20190329.nc'),
                            preprocess=preprocess_add_r)
# ds_cont_large = ds_cont_raw.sel(latitude=slice(bbox_idai[1], bbox_idai[3]),
#                       longitude=slice(bbox_idai[0], bbox_idai[2]))

ds_hist = xr.open_mfdataset(os.path.join(datadir, 'climateDT_*_hist_*_Idai_full_20190301_to_20190329.nc'),
                            preprocess=preprocess_add_r)
# ds_hist_large = ds_hist.sel(latitude=slice(bbox_idai[1], bbox_idai[3]),
#                       longitude=slice(bbox_idai[0], bbox_idai[2]))

ds_tp2k = xr.open_mfdataset(os.path.join(datadir, 'climateDT_*_Tplus2.0K_*_Idai_full_20190301_to_20190329.nc'),
                            preprocess=preprocess_add_r)

# Converting precipitation from mm to m
ds_cont['tp'] = ds_cont['tp']/1000
ds_cont['tp'].attrs['units'] = 'm'

ds_hist['tp'] = ds_hist['tp']/1000
ds_hist['tp'].attrs['units'] = 'm'

ds_tp2k['tp'] = ds_tp2k['tp']/1000
ds_tp2k['tp'].attrs['units'] = 'm'

#%%
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


# %%
# Get masks for wflow region, and Pungwe and Buzi basins 
mask_da = get_mask_from_gdf(gdf_wflow, ds_cont)
mask_da_era5 = get_mask_from_gdf(gdf_wflow, ds_tp_era5)

gdf_Pungwe = gdf_wflow[gdf_wflow['geometry'].area>2][gdf_wflow['value']==42]
gdf_Buzi = gdf_wflow[gdf_wflow['geometry'].area>2][gdf_wflow['value']==67]

mask_da_pungwe = get_mask_from_gdf(gdf_Pungwe, ds_cont)
mask_da_buzi = get_mask_from_gdf(gdf_Buzi, ds_cont)


#%%
#### Comparing meteo forcings #####
# PRECIPITATION
# preparing data
times=['2019-03-13','2019-03-21']

# Plot area-aggragated precipitation (conversion to m^3)
cdt_x = np.mean(ds_cont.longitude.values[1:]-ds_cont.longitude.values[:-1])*111.32*1000*np.cos(np.deg2rad(ds_cont.latitude.mean().values))  # approx. conversion to m^2
cdt_y = np.mean(ds_cont.latitude.values[1:]-ds_cont.latitude.values[:-1])*110.574*1000  # approx. conversion to m^2

era5_x = np.mean(ds_tp_era5.longitude.values[1:]-ds_tp_era5.longitude.values[:-1])*111.32*1000*np.cos(np.deg2rad(ds_tp_era5.latitude.mean().values))
era5_y = np.mean(ds_tp_era5.latitude.values[1:]-ds_tp_era5.latitude.values[:-1])*110.574*1000 

# --- TOTAL VOLUMES over selected period ---
hist_totals = [(ds_hist.sel(realization=rr) * cdt_x * cdt_y).where(mask_da).sel(time=slice(times[0], times[1])).tp.sum(dim=['time', 'latitude', 'longitude']).values
               for rr in ds_hist.realization.values]
cont_totals = [(ds_cont.sel(realization=rr) * cdt_x * cdt_y).where(mask_da).sel(time=slice(times[0], times[1])).tp.sum(dim=['time', 'latitude', 'longitude']).values
               for rr in ds_cont.realization.values]
tp2k_totals = [(ds_tp2k.sel(realization=rr) * cdt_x * cdt_y).where(mask_da).sel(time=slice(times[0], times[1])).tp.sum(dim=['time', 'latitude', 'longitude']).values
               for rr in ds_tp2k.realization.values]
era5_total = (ds_tp_era5 * era5_x * era5_y).where(mask_da_era5).sel(time=slice(times[0], times[1])).sum(dim=['time', 'latitude', 'longitude']).values

#%%
peak1 = ("2019-03-14", "2019-03-15 18:00:00")
peak2 = ("2019-03-17 13:00:00", "2019-03-18 19:00:00")
colors = ['tab:orange', 'tab:blue', 'tab:green']

hist_peak1 = [(ds_hist.sel(realization=rr) * cdt_x * cdt_y).where(mask_da).sel(time=slice(*peak1))
              .tp.sum(dim=['time', 'latitude', 'longitude']).values for rr in ds_hist.realization.values]
cont_peak1 = [(ds_cont.sel(realization=rr) * cdt_x * cdt_y).where(mask_da).sel(time=slice(*peak1))
              .tp.sum(dim=['time', 'latitude', 'longitude']).values for rr in ds_cont.realization.values]
tp2k_peak1 = [(ds_tp2k.sel(realization=rr)*cdt_x*cdt_y).where(mask_da).sel(time=slice(*peak1))
              .tp.sum(dim=['time', 'latitude', 'longitude']).values for rr in ds_tp2k.realization.values]
era5_peak1 = (ds_tp_era5*era5_x*era5_y).where(mask_da_era5).sel(time=slice(*peak1)).sum(dim=['time', 'latitude', 'longitude']).values
hist_peak2 = [(ds_hist.sel(realization=rr)*cdt_x*cdt_y).where(mask_da).sel(time=slice(*peak2))
              .tp.sum(dim=['time', 'latitude', 'longitude']).values for rr in ds_hist.realization.values]
tp2k_peak1 = [(ds_tp2k.sel(realization=rr)*cdt_x*cdt_y).where(mask_da).sel(time=slice(*peak1))
              .tp.sum(dim=['time', 'latitude', 'longitude']).values for rr in ds_tp2k.realization.values]
hist_peak2 = [(ds_hist.sel(realization=rr)*cdt_x*cdt_y).where(mask_da).sel(time=slice(*peak2))
              .tp.sum(dim=['time', 'latitude', 'longitude']).values for rr in ds_hist.realization.values]
cont_peak2 = [(ds_cont.sel(realization=rr)*cdt_x*cdt_y).where(mask_da).sel(time=slice(*peak2))
              .tp.sum(dim=['time', 'latitude', 'longitude']).values for rr in ds_cont.realization.values]
tp2k_peak2 = [(ds_tp2k.sel(realization=rr)*cdt_x*cdt_y).where(mask_da).sel(time=slice(*peak2))
              .tp.sum(dim=['time', 'latitude', 'longitude']).values for rr in ds_tp2k.realization.values]
era5_peak2 = (ds_tp_era5*era5_x*era5_y).where(mask_da_era5).sel(time=slice(*peak2)).sum(dim=['time', 'latitude', 'longitude']).values

#%%
#### plotting #####
fig = plt.figure(figsize=(14, 8))
gs = fig.add_gridspec(2, 3, height_ratios=[2, 1.3], hspace=0.25)

ax0 = fig.add_subplot(gs[0, :])  # full top row
ax1 = fig.add_subplot(gs[1, 0])
ax2 = fig.add_subplot(gs[1, 1])
ax3 = fig.add_subplot(gs[1, 2])

p1 = (ds_cont.median(dim='realization')*cdt_x*cdt_y).where(mask_da).sel(time=slice(times[0],times[1])).tp.sum(dim=['latitude','longitude']).plot(ax=ax0, label=f'CF_past', color='tab:orange', linewidth=2)
p2 = (ds_hist.median(dim='realization')*cdt_x*cdt_y).where(mask_da).sel(time=slice(times[0],times[1])).tp.sum(dim=['latitude','longitude']).plot(ax=ax0, label=f'Factual', color='tab:blue', linewidth=2)
p3 = (ds_tp2k.median(dim='realization')*cdt_x*cdt_y).where(mask_da).sel(time=slice(times[0],times[1])).tp.sum(dim=['latitude','longitude']).plot(ax=ax0, label=f'CF_future', color='tab:green', linewidth=2)
p4 = (ds_tp_era5*era5_x*era5_y).where(mask_da_era5).sel(time=slice(times[0],times[1])).sum(dim=['latitude','longitude']).plot(ax=ax0, label='ERA5 reanalysis', color='k')

# Calculate the total precipitation volume over the masked area as a function of time
for rr in ds_hist.realization.values:
    (ds_cont.sel(realization=rr)*cdt_x*cdt_y).where(mask_da).sel(time=slice(times[0],times[1])).tp.sum(dim=['latitude','longitude']).plot(ax=ax0, color='tab:orange', alpha=0.3)
    (ds_hist.sel(realization=rr)*cdt_x*cdt_y).where(mask_da).sel(time=slice(times[0],times[1])).tp.sum(dim=['latitude','longitude']).plot(ax=ax0, color='tab:blue', alpha=0.3)
    (ds_tp2k.sel(realization=rr)*cdt_x*cdt_y).where(mask_da).sel(time=slice(times[0],times[1])).tp.sum(dim=['latitude','longitude']).plot(ax=ax0, color='tab:green', alpha=0.3)

# Highlighting the two peaks in precipitation
ax0.axvspan(np.datetime64(peak1[0]), np.datetime64(peak1[1]), color="grey", alpha=0.2, label="Peak 1")
ax0.axvspan(np.datetime64(peak2[0]), np.datetime64(peak2[1]), color="grey", alpha=0.4, label="Peak 2")

ax0.set_title("")
ax0.legend()
ax0.grid()
ax0.set_ylabel('Volume (m$^3$)')
ax0.set_xlim([np.datetime64(times[0]), np.datetime64(times[1])])


axes = [ax1, ax2, ax3]
titles = ["Total event", "First rainfall peak", "Second rainfall peak"]

datasets = [[cont_totals, hist_totals, tp2k_totals], [cont_peak1, hist_peak1, tp2k_peak1],
            [cont_peak2, hist_peak2, tp2k_peak2]]

for ax, data, title in zip(axes, datasets, titles):
    bp = ax.boxplot(data, labels=['CF_past', 'Factual', 'CF_future'], patch_artist=True)
 
    for i, color in enumerate(colors):
        bp['boxes'][i].set_facecolor(color)
        bp['boxes'][i].set_alpha(0.4)

        bp['medians'][i].set_color(color)
        bp['medians'][i].set_linewidth(2)

        bp['whiskers'][2*i].set_color(color)
        bp['whiskers'][2*i+1].set_color(color)

        bp['caps'][2*i].set_color(color)
        bp['caps'][2*i+1].set_color(color)

    ax.set_title(title)
    ax.grid(axis='y')

    # Median percentage changes relative to Factual
    med_cf_past = np.median(data[0])
    med_hist = np.median(data[1])
    med_cf_future = np.median(data[2])

    change_cf_past = 100 * (med_cf_past - med_hist) / med_hist
    change_cf_future = 100 * (med_cf_future - med_hist) / med_hist

    y_cf_past = np.max(data[0]) * 1.05
    y_cf_future = np.max(data[2]) * 1.05

    ax.text(
        1, y_cf_past,
        f"{change_cf_past:+.1f}%",
        ha="center",
        color="tab:orange",
        fontweight="bold"
    )

    ax.text(
        3, y_cf_future,
        f"{change_cf_future:+.1f}%",
        ha="center",
        color="tab:green",
        fontweight="bold"
    )

    ymax = max(np.max(d) for d in data)
    ax.set_ylim(top=ymax * 1.2)


ax1.axhline(era5_total, color='k', linestyle='--')
ax2.axhline(era5_peak1, color='k', linestyle='--')
ax3.axhline(era5_peak2, color='k', linestyle='--')
ax3.set_ylim(top=era5_peak2 * 1.05)

ax1.set_ylabel('Volume (m$^3$)')

fig.suptitle(f'Area-aggregated precipitation over the wflow+SFINCS model regions', 
             fontweight='bold', y=0.91)


#%%
# %%
# Spatial comaprison of precipitation for the selected period (Factual)
lims=[0, 1.1]

fig, axs = plt.subplots(figsize=(10, 8), ncols=3, nrows=2, sharex=True, sharey=True, constrained_layout=True)

# Plot ERA5 
p = ds_tp_era5.sel(time=slice(times[0], times[1])).sum(dim='time').plot(
    ax=axs[0, 0], vmin=lims[0], vmax=lims[1], cmap='turbo', add_colorbar=False)
axs[0, 0].set_xlim([bbox_wflow_sfincs[0], bbox_wflow_sfincs[2]])
axs[0, 0].set_ylim([bbox_wflow_sfincs[1], bbox_wflow_sfincs[3]])
axs[0, 0].set_title('ERA5 reanalysis')

# Plot ClimateDT historical scenario
for rr,ax in enumerate(axs.ravel()[1:]):
    ds_hist.sel(time=slice(times[0], times[1]), realization=rr+1).tp.sum(dim='time').plot(
        ax=ax, vmin=lims[0], vmax=lims[1], cmap='turbo', add_colorbar=False)
    ax.set_title(f'ClimateDT factual scenario, r{rr+1}')
    if rr != 2:
        ax.set_ylabel('')

for ax in axs.ravel():
    gdf_wflow.plot(ax=ax, edgecolor="gray", facecolor="none")
    gdf_sfincs.plot(ax=ax, edgecolor="black", facecolor="none")
    ax.set_xlim(bbox_wflow_sfincs[0], bbox_wflow_sfincs[2])
    ax.set_ylim(bbox_wflow_sfincs[1], bbox_wflow_sfincs[3])

cbar = fig.colorbar(p, ax=axs.ravel().tolist(), orientation="horizontal", pad=0.01, 
                    aspect=40, fraction=0.1, shrink=0.8)
cbar.set_label("Cumulative precipitation [m]", labelpad=6)
fig.set_constrained_layout_pads(wspace=0.1)
fig.suptitle(f'Cumulative precipitation from {times[0]} to {times[1]}')


#%%
# Total TC period - ERA5 and Factual realizations
lims=[0, 1.1]

fig, axs = plt.subplots(figsize=(10, 8), ncols=3, nrows=2, sharex=True, sharey=True, constrained_layout=True)

# Plot ERA5 
p = ds_tp_era5.sel(time=slice(times[0], times[1])).sum(dim='time').plot(
    ax=axs[0, 0], vmin=lims[0], vmax=lims[1], cmap='turbo', add_colorbar=False)
axs[0, 0].set_xlim([bbox_wflow_sfincs[0], bbox_wflow_sfincs[2]])
axs[0, 0].set_ylim([bbox_wflow_sfincs[1], bbox_wflow_sfincs[3]])
axs[0, 0].set_title('ERA5 reanalysis')

# Plot ClimateDT historical scenario
for rr,ax in enumerate(axs.ravel()[1:]):
    ds_hist.sel(time=slice(times[0], times[1]), realization=rr+1).tp.sum(dim='time').plot(
        ax=ax, vmin=lims[0], vmax=lims[1], cmap='turbo', add_colorbar=False)
    ax.set_title(f'ClimateDT factual scenario, r{rr+1}')
    if rr != 2:
        ax.set_ylabel('')

for ax in axs.ravel():
    gdf_wflow.plot(ax=ax, edgecolor="gray", facecolor="none")
    gdf_sfincs.plot(ax=ax, edgecolor="black", facecolor="none")
    ax.set_xlim(bbox_wflow_sfincs[0], bbox_wflow_sfincs[2])
    ax.set_ylim(bbox_wflow_sfincs[1], bbox_wflow_sfincs[3])

cbar = fig.colorbar(p, ax=axs.ravel().tolist(), orientation="horizontal", pad=0.01, 
                    aspect=40, fraction=0.1, shrink=0.8)
cbar.set_label("Cumulative precipitation [m]", labelpad=6)
fig.set_constrained_layout_pads(wspace=0.1)
fig.suptitle(f'Cumulative precipitation from {times[0]} to {times[1]}')


#%%
periods = {"Total_event": (times[0], times[1]), "Peak1": peak1, "Peak2": peak2}

stats_tp = {}
for name, period in periods.items():
    hist = ds_hist.sel(time=slice(*period)).tp.sum("time")
    cont = ds_cont.sel(time=slice(*period)).tp.sum("time")
    tp2k = ds_tp2k.sel(time=slice(*period)).tp.sum("time")

    stats_tp[name] = {
        "hist_mean": hist.mean("realization"),
        "cont_mean": cont.mean("realization"),
        "tp2k_mean": tp2k.mean("realization"),

        "hist_median": hist.median("realization"),
        "cont_median": cont.median("realization"),
        "tp2k_median": tp2k.median("realization"),

        "diff_mean_CF_past": cont.mean("realization") - hist.mean("realization"),
        "diff_mean_CF_future": tp2k.mean("realization") - hist.mean("realization"),
        "diff_median_CF_past": cont.median("realization") - hist.median("realization"),
        "diff_median_CF_future": tp2k.median("realization") - hist.median("realization"),
        }
    
#%%
summary_tp_change = pd.DataFrame({
        "CF_past_mean [%]": [((stats_tp[p]["cont_mean"].sum() - stats_tp[p]["hist_mean"].sum())
                              / stats_tp[p]["hist_mean"].sum() * 100).values for p in periods],
        "CF_future_mean [%]": [((stats_tp[p]["tp2k_mean"].sum() - stats_tp[p]["hist_mean"].sum())
                                / stats_tp[p]["hist_mean"].sum() * 100).values for p in periods],
        "CF_past_median [%]": [((stats_tp[p]["cont_median"].sum() - stats_tp[p]["hist_median"].sum())
                                / stats_tp[p]["hist_median"].sum() * 100).values for p in periods        ],
        "CF_future_median [%]": [((stats_tp[p]["tp2k_median"].sum() - stats_tp[p]["hist_median"].sum())
                / stats_tp[p]["hist_median"].sum() * 100).values for p in periods],},
                index=periods.keys())

print(summary_tp_change.round(1))

#%%
# Total TC period - Factual mean, CF past and CF future mean
lims=[0, 1]

fig, axs = plt.subplots(figsize=(10, 5), ncols=3, nrows=1, sharex=True, sharey=True, constrained_layout=True)

# Plot ERA5 
p1 = stats_tp['Total_event']['hist_mean'].plot(ax=axs[0], vmin=lims[0], vmax=lims[1], cmap='turbo', add_colorbar=False)
p1 = stats_tp['Total_event']['cont_mean'].plot(ax=axs[1], vmin=lims[0], vmax=lims[1], cmap='turbo', add_colorbar=False)
p1 = stats_tp['Total_event']['tp2k_mean'].plot(ax=axs[2], vmin=lims[0], vmax=lims[1], cmap='turbo', add_colorbar=False)

for i, ax in enumerate(axs):
    gdf_wflow.plot(ax=ax, edgecolor="gray", facecolor="none")
    gdf_sfincs.plot(ax=ax, edgecolor="black", facecolor="none")
    ax.set_xlim(bbox_wflow_sfincs[0], bbox_wflow_sfincs[2])
    ax.set_ylim(bbox_wflow_sfincs[1], bbox_wflow_sfincs[3])
    if i != 0:
        ax.set_ylabel('')

axs[0].set_title('Factual mean')
axs[1].set_title('CF past mean')
axs[2].set_title('CF future mean')

cbar = fig.colorbar(p1, ax=axs[:], orientation="vertical", pad=0.01, 
                    aspect=40, fraction=0.1, shrink=0.8)
cbar.set_label("Cumulative precipitation [m]", labelpad=6)
fig.set_constrained_layout_pads(wspace=0.1)
fig.suptitle(f'Mean cumulative precipitation from {times[0]} to {times[1]}', fontweight='bold')


#%%
# Total TC period - Factual mean and absolute difference from CF scenario means
fig, axs = plt.subplots(figsize=(10, 4), ncols=3, nrows=1, sharex=True, sharey=True, 
                        constrained_layout=True)

# Plot mean precipitation from the historical scenario
p1 = stats_tp['Total_event']['hist_mean'].plot(ax=axs[0], vmin=0, vmax=stats_tp['Total_event']['hist_mean'].max().values, cmap='turbo', add_colorbar=False)
axs[0].set_title('Factual mean')

# Plot absolute difference
p2 = stats_tp['Total_event']['diff_mean_CF_past'].plot(ax=axs[1], cmap='RdBu', vmin=-0.2, vmax=0.2, add_colorbar=False)
p3 = stats_tp['Total_event']['diff_mean_CF_future'].plot(ax=axs[2], cmap='RdBu', vmin=-0.2, vmax=0.2, add_colorbar=False)
axs[1].set_title('CF past mean')
axs[2].set_title('CF future mean')

for i, ax in enumerate(axs):
    gdf_wflow.plot(ax=ax, edgecolor="gray", facecolor="none")
    gdf_sfincs.plot(ax=ax, edgecolor="black", facecolor="none")
    ax.set_xlim(bbox_wflow_sfincs[0], bbox_wflow_sfincs[2])
    ax.set_ylim(bbox_wflow_sfincs[1], bbox_wflow_sfincs[3])
    if i != 0:
        ax.set_ylabel('')

cbar1 = fig.colorbar(p1, ax=axs[0], orientation="vertical", pad=0.05, aspect=20, 
                    fraction=0.1, shrink=0.6)
cbar1.set_label("Cumulative mean precipitation [m]", labelpad=6)

cbar2 = fig.colorbar(p2, ax=axs[1:], orientation="vertical", pad=0.05, aspect=20, 
                    fraction=0.1, shrink=0.6)
cbar2.set_label("Absolute difference [m]", labelpad=6)

fig.set_constrained_layout_pads(wspace=0.1)
fig.suptitle(f'Mean cumulative precipitation from {times[0]} to {times[1]}', fontweight='bold')

# Total TC period - Factual mean and absolute difference from CF scenario means
fig, axs = plt.subplots(figsize=(10, 4), ncols=3, nrows=1, sharex=True, sharey=True, 
                        constrained_layout=True)

#%%
# Total TC period - Factual median and absolute difference from CF scenario medians
fig, axs = plt.subplots(figsize=(10, 4), ncols=3, nrows=1, sharex=True, sharey=True, 
                        constrained_layout=True)

# Plot median precipitation from the historical scenario
p1 = stats_tp['Total_event']['hist_median'].plot(ax=axs[0], vmin=0, vmax=stats_tp['Total_event']['hist_median'].max().values, cmap='turbo', add_colorbar=False)
axs[0].set_title('Factual median')

# Plot absolute difference
p2 = stats_tp['Total_event']['diff_median_CF_past'].plot(ax=axs[1], cmap='RdBu', vmin=-0.2, vmax=0.2, add_colorbar=False)
p3 = stats_tp['Total_event']['diff_median_CF_future'].plot(ax=axs[2], cmap='RdBu', vmin=-0.2, vmax=0.2, add_colorbar=False)
axs[1].set_title('CF past median')
axs[2].set_title('CF future median')

for i, ax in enumerate(axs):
    gdf_wflow.plot(ax=ax, edgecolor="gray", facecolor="none")
    gdf_sfincs.plot(ax=ax, edgecolor="black", facecolor="none")
    ax.set_xlim(bbox_wflow_sfincs[0], bbox_wflow_sfincs[2])
    ax.set_ylim(bbox_wflow_sfincs[1], bbox_wflow_sfincs[3])
    if i != 0:
        ax.set_ylabel('')

cbar1 = fig.colorbar(p1, ax=axs[0], orientation="vertical", pad=0.05, aspect=20, 
                    fraction=0.1, shrink=0.6)
cbar1.set_label("Cumulative median precipitation [m]", labelpad=6)

cbar2 = fig.colorbar(p2, ax=axs[1:], orientation="vertical", pad=0.05, aspect=20, 
                    fraction=0.1, shrink=0.6)
cbar2.set_label("Absolute difference [m]", labelpad=6)

fig.set_constrained_layout_pads(wspace=0.1)
fig.suptitle(f'Median cumulative precipitation from {times[0]} to {times[1]}', fontweight='bold')


#%%
# Peak 1 of TC period - Factual mean and absolute difference from CF scenario means
fig, axs = plt.subplots(figsize=(10, 4), ncols=3, nrows=1, sharex=True, sharey=True, 
                        constrained_layout=True)

# Plot mean precipitation from the historical scenario
p1 = stats_tp['Peak1']['hist_mean'].plot(ax=axs[0], vmin=0, vmax=stats_tp['Peak1']['hist_mean'].max().values, cmap='turbo', add_colorbar=False)
axs[0].set_title('Factual mean')

# Plot absolute difference
p2 = stats_tp['Peak1']['diff_mean_CF_past'].plot(ax=axs[1], cmap='RdBu', vmin=-0.1, vmax=0.1, add_colorbar=False)
p3 = stats_tp['Peak1']['diff_mean_CF_future'].plot(ax=axs[2], cmap='RdBu', vmin=-0.1, vmax=0.1, add_colorbar=False)
axs[1].set_title('CF past mean')
axs[2].set_title('CF future mean')

for i, ax in enumerate(axs):
    gdf_wflow.plot(ax=ax, edgecolor="gray", facecolor="none")
    gdf_sfincs.plot(ax=ax, edgecolor="black", facecolor="none")
    ax.set_xlim(bbox_wflow_sfincs[0], bbox_wflow_sfincs[2])
    ax.set_ylim(bbox_wflow_sfincs[1], bbox_wflow_sfincs[3])
    if i != 0:
        ax.set_ylabel('')

cbar1 = fig.colorbar(p1, ax=axs[0], orientation="vertical", pad=0.05, aspect=20, 
                    fraction=0.1, shrink=0.6)
cbar1.set_label("Cumulative mean precipitation [m]", labelpad=6)

cbar2 = fig.colorbar(p2, ax=axs[1:], orientation="vertical", pad=0.05, aspect=20, 
                    fraction=0.1, shrink=0.6)
cbar2.set_label("Absolute difference [m]", labelpad=6)

fig.set_constrained_layout_pads(wspace=0.1)
fig.suptitle(f'Mean cumulative precipitation from {peak1[0]} to {peak1[1]}', fontweight='bold')


#%%
# Peak 1 of TC period - Factual median and absolute difference from CF scenario medians
fig, axs = plt.subplots(figsize=(8, 8), ncols=2, nrows=2, sharex=True, sharey=True, 
                        constrained_layout=True)

# Plot median precipitation from the historical scenario
p0 = ds_tp_era5.sel(time=slice(*peak1)).sum(dim='time').plot(ax=axs[0, 0], vmin=0, vmax=stats_tp['Peak1']['hist_median'].max().values, cmap='turbo', add_colorbar=False)
axs[0, 0].set_title('ERA5 reanalysis')
p1 = stats_tp['Peak1']['hist_median'].plot(ax=axs[0, 1], vmin=0, vmax=stats_tp['Peak1']['hist_median'].max().values, cmap='turbo', add_colorbar=False)
axs[0, 1].set_title('Factual median')

# Plot absolute difference
p2 = stats_tp['Peak1']['diff_median_CF_past'].plot(ax=axs[1, 0], cmap='RdBu', vmin=-0.1, vmax=0.1, add_colorbar=False)
p3 = stats_tp['Peak1']['diff_median_CF_future'].plot(ax=axs[1, 1], cmap='RdBu', vmin=-0.1, vmax=0.1, add_colorbar=False)
axs[1, 0].set_title('CF past median')
axs[1, 1].set_title('CF future median')

for i, ax in enumerate(axs.flatten()):
    gdf_wflow.plot(ax=ax, edgecolor="gray", facecolor="none")
    gdf_sfincs.plot(ax=ax, edgecolor="black", facecolor="none")
    ax.set_xlim(bbox_wflow_sfincs[0], bbox_wflow_sfincs[2])
    ax.set_ylim(bbox_wflow_sfincs[1], bbox_wflow_sfincs[3])
    if i != 0:
        ax.set_ylabel('')

cbar1 = fig.colorbar(p1, ax=axs[:1], orientation="vertical", pad=0.05, aspect=20, 
                    fraction=0.1, shrink=0.6)
cbar1.set_label("Cumulative median precipitation [m]", labelpad=6)

cbar2 = fig.colorbar(p2, ax=axs[1:], orientation="vertical", pad=0.05, aspect=20, 
                    fraction=0.1, shrink=0.6)
cbar2.set_label("Absolute difference [m]", labelpad=6)

fig.set_constrained_layout_pads(wspace=0.1)
fig.suptitle(f'Median cumulative precipitation from {peak1[0]} to {peak1[1]}', fontweight='bold')


#%%
# Peak 2 of TC period - Factual mean and absolute difference from CF scenario means
fig, axs = plt.subplots(figsize=(10, 4), ncols=3, nrows=1, sharex=True, sharey=True, 
                        constrained_layout=True)

# Plot mean precipitation from the historical scenario
p1 = stats_tp['Peak2']['hist_mean'].plot(ax=axs[0], vmin=0, vmax=stats_tp['Peak2']['hist_mean'].max().values, cmap='turbo', add_colorbar=False)
axs[0].set_title('Factual mean')

# Plot absolute difference
p2 = stats_tp['Peak2']['diff_mean_CF_past'].plot(ax=axs[1], cmap='RdBu', vmin=-0.1, vmax=0.1, add_colorbar=False)
p3 = stats_tp['Peak2']['diff_mean_CF_future'].plot(ax=axs[2], cmap='RdBu', vmin=-0.1, vmax=0.1, add_colorbar=False)
axs[1].set_title('CF past mean')
axs[2].set_title('CF future mean')

for i, ax in enumerate(axs):
    gdf_wflow.plot(ax=ax, edgecolor="gray", facecolor="none")
    gdf_sfincs.plot(ax=ax, edgecolor="black", facecolor="none")
    ax.set_xlim(bbox_wflow_sfincs[0], bbox_wflow_sfincs[2])
    ax.set_ylim(bbox_wflow_sfincs[1], bbox_wflow_sfincs[3])
    if i != 0:
        ax.set_ylabel('')

cbar1 = fig.colorbar(p1, ax=axs[0], orientation="vertical", pad=0.05, aspect=20, 
                    fraction=0.1, shrink=0.6)
cbar1.set_label("Cumulative mean precipitation [m]", labelpad=6)

cbar2 = fig.colorbar(p2, ax=axs[1:], orientation="vertical", pad=0.05, aspect=20, 
                    fraction=0.1, shrink=0.6)
cbar2.set_label("Absolute difference [m]", labelpad=6)

fig.set_constrained_layout_pads(wspace=0.1)
fig.suptitle(f'Mean cumulative precipitation from {peak2[0]} to {peak2[1]}', fontweight='bold')

#%%
# Peak 2 of TC period - Factual median and absolute difference from CF scenario means
fig, axs = plt.subplots(figsize=(8, 8), ncols=2, nrows=2, sharex=True, sharey=True, 
                        constrained_layout=True)

# Plot mean precipitation from the historical scenario
p0 = ds_tp_era5.sel(time=slice(*peak2)).sum(dim='time').plot(ax=axs[0, 0], vmin=0, vmax=stats_tp['Peak2']['hist_median'].max().values, cmap='turbo', add_colorbar=False)
axs[0, 0].set_title('ERA5 reanalysis')
p1 = stats_tp['Peak2']['hist_median'].plot(ax=axs[0, 1], vmin=0, vmax=stats_tp['Peak2']['hist_median'].max().values, cmap='turbo', add_colorbar=False)
axs[0, 1].set_title('Factual median')

# Plot absolute difference
p2 = stats_tp['Peak2']['diff_median_CF_past'].plot(ax=axs[1, 0], cmap='RdBu', vmin=-0.1, vmax=0.1, add_colorbar=False)
p3 = stats_tp['Peak2']['diff_median_CF_future'].plot(ax=axs[1, 1], cmap='RdBu', vmin=-0.1, vmax=0.1, add_colorbar=False)
axs[1, 0].set_title('CF past median')
axs[1, 1].set_title('CF future median')

for i, ax in enumerate(axs.flatten()):
    gdf_wflow.plot(ax=ax, edgecolor="gray", facecolor="none")
    gdf_sfincs.plot(ax=ax, edgecolor="black", facecolor="none")
    ax.set_xlim(bbox_wflow_sfincs[0], bbox_wflow_sfincs[2])
    ax.set_ylim(bbox_wflow_sfincs[1], bbox_wflow_sfincs[3])
    if i != 0:
        ax.set_ylabel('')

cbar1 = fig.colorbar(p1, ax=axs[:1], orientation="vertical", pad=0.05, aspect=20, 
                    fraction=0.1, shrink=0.6)
cbar1.set_label("Cumulative median precipitation [m]", labelpad=6)

cbar2 = fig.colorbar(p2, ax=axs[1:], orientation="vertical", pad=0.05, aspect=20, 
                    fraction=0.1, shrink=0.6)
cbar2.set_label("Absolute difference [m]", labelpad=6)

fig.set_constrained_layout_pads(wspace=0.1)
fig.suptitle(f'Median cumulative precipitation from {peak2[0]} to {peak2[1]}', fontweight='bold')

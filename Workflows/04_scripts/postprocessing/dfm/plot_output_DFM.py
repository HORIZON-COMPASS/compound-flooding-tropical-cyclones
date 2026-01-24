#%% Import the necessary packages using pixi environment compass-snake-dfm
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import dfm_tools as dfmt
import os
import cartopy.crs as ccrs
import ast
import hydromt

if "snakemake" in locals():
    dir_runs = os.path.abspath(snakemake.output.dir_event_model) 
    model = snakemake.params.model_name
    sfincs_bbox = ast.literal_eval(snakemake.params.sfincs_bbox)
    dfm_bbox = ast.literal_eval(snakemake.params.dfm_bbox)
else:
    # dir_runs = r'p:\11210471-001-compass\02_Models\Delft3DFM\mozambique_model\computations'
    # model = 'mozambique_spw_Idai_area_MZB_500m_gebco2024' 
    region = "sofala"
    tc_name = "Idai"
    dfm_res = "450"
    bathy = "gebco2024_MZB"
    tidemodel = 'GTSMv41' # tidemodel: FES2014, FES2012, EOT20, GTSMv41, GTSMv41opendap
    wind_forcing = "era5_hourly_spw_IBTrACS"
    CF_SLR = 0
    CF_SLR_txt = "0"
    CF_wind = 0
    CF_wind_txt = "0"
    dir_runs = f'p:/11210471-001-compass/03_Runs/{region}/{tc_name}/dfm'
    # model_name = f'event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}_nochnk'
    model_name = f'event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}'
    dfm_bbox_model = ast.literal_eval("[32.3,42.5,-27.4,-9.5]")   
    crop_bbox = ast.literal_eval("[34, -20.2, 35.6, -19.2]")
    sfincs_bbox = ast.literal_eval("[34.33,-20.12,34.95,-19.30]")
    dfm_obs_file = "p:/11210471-001-compass/01_Data/Coastal_boundary/points/coastal_bnd_MZB_5mMSL_points_1km.shp"
    

# Open the his files in the output folder
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


# Open the map files in the output folder
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

# Some functions to plot the model results
def plot_timeserie(data, variable, station_name='BEIRA IHO'):
    # Plot the stations
    fig, ax = plt.subplots(1,1,figsize=(8,5))
    # Find the index of the station
    station_idx  = data.station.values.tolist().index(station_name)

    # Plot the water level for the selected station over time for the different simulations
    ax.plot(data.time.values, data[variable][:, station_idx], label=f'{station_name}')
    print(data[variable].attrs['long_name'])
    # Set labels and title
    unit = data[variable].attrs['units']
    long_name = data[variable].attrs['long_name']
    ax.set_ylabel(f'{long_name} ({unit})')
    ax.set_xlabel('Time')
    ax.set_title(f'{variable} for Station {station_name}')
    # Add legend
    ax.legend(loc=1, fontsize=8)


def plot_timeseries_multi(datasets, variable, station_name='BEIRA IHO', labels=None, start_date=None, end_date=None, y_max=None, y_min=None):
    """
    Plot a time series for the same variable and station across multiple datasets,
    with varying line thicknesses for better visibility.
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    if labels is not None and len(labels) != len(datasets):
        raise ValueError("Length of labels must match length of datasets.")

    for i, data in enumerate(datasets):
        if station_name not in data.station.values:
            print(f"Station '{station_name}' not found in dataset {i}")
            continue
        station_idx = data.station.values.tolist().index(station_name)

        label = labels[i] if labels is not None else f'Dataset {i+1}'
        line_width = 6 - i * 1  # Increase line width for each dataset
        ax.plot(data.time.values, data[variable][:, station_idx], label=label)

    # Use metadata from the first dataset for labels
    unit = datasets[0][variable].attrs.get('units', '')
    long_name = datasets[0][variable].attrs.get('long_name', variable)

    ax.set_ylabel(f'{long_name} ({unit})')
    ax.set_xlabel('Time')
    ax.set_title(f'{long_name} at Station {station_name}')
    ax.set_xlim(start_date, end_date)
    ax.set_ylim(y_min, y_max)
    ax.legend(loc='upper right', fontsize=8)
    plt.tight_layout()
    plt.show()


def plot_wind_slice(ds, time_str, var='mesh2d_windmag', title=None):
    """Plot wind magnitude at a specific time from a UGRID dataset."""
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10, 6))

    da = ds[var].sel(time=time_str)
    da.ugrid.plot(ax=ax, transform=ccrs.PlateCarree(), cmap='viridis')

    ax.coastlines()
    ax.set_title(title or f"Wind magnitude at {time_str}")
    return ax


def plot_wind_difference(ds1, ds2, time_str, var='mesh2d_windmag', title=None):
    """Plot wind magnitude difference between two datasets at a specific time."""
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10, 6))

    da1 = ds1[var].sel(time=time_str)
    da2 = ds2[var].sel(time=time_str)

    da_diff = da1 - da2
    da_diff.ugrid.plot(ax=ax, transform=ccrs.PlateCarree(), cmap='RdBu_r')

    ax.coastlines()
    ax.set_title(title or f"Wind magnitude difference at {time_str}")
    return ax


def compute_and_plot_max_waterlevel(ds, var='mesh2d_s1', title="Maximum water level", add_coastlines=True):
    """
    Efficiently compute and plot the maximum water level from a UGRID dataset.

    Parameters:
        ds (xarray.Dataset): Input dataset with UGRID mesh.
        var (str): Variable name for water level.
        title (str): Title of the plot.
        add_coastlines (bool): If True, adds coastlines to the plot.
    
    Returns:
        matplotlib.axes.Axes: The plot axis.
    """
    max_var = f'max_{var}'

    # Cache max value if not already in dataset
    if max_var not in ds:
        ds[max_var] = ds[var].max(dim='time', keep_attrs=True)

    da_max = ds[max_var]

    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    da_max.ugrid.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='Blues',
        rasterized=True
    )

    if add_coastlines:
        ax.coastlines(resolution='10m')

    ax.set_title(title)
    return ax


def plot_waterlevel_difference(ds1, ds2, var='mesh2d_s1', title='Water level difference', time_dim='time'):
    # Compute max water level over time
    da1_max = ds1[var].max(dim=time_dim)
    da2_max = ds2[var].max(dim=time_dim)

    # Calculate the difference
    diff = da1_max - da2_max

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    diff.ugrid.plot(ax=ax, transform=ccrs.PlateCarree(), cmap='RdBu_r', center=0, rasterized=True)

    ax.coastlines(resolution='10m')
    ax.set_title(title)
    return ax


#%%
# Load the data catalog with all the different DFM outputs with different SLR and wind CFs
path_data_cat       = os.path.abspath("../../../03_data_catalogs/datacatalog_CF_forcing.yml")
datacatalog = hydromt.DataCatalog(data_libs=[path_data_cat])


#%%
# Load all the different DFM outputs with different SLR and wind CFs
dfm_run = 'dfm_output_event_450_gebco2024_MZB_GTSMv41_CF-0.05_era5_hourly_spw_IBTrACS_CF-1_waves'
ds_dfm_SLR_0cm_wind_0 = datacatalog.get_geodataframe('dfm_output_event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_waves')
ds_dfm_SLR_5cm_wind_0 = datacatalog.get_geodataframe('dfm_output_event_450_gebco2024_MZB_GTSMv41_CF-0.05_era5_hourly_spw_IBTrACS_CF0_waves')
ds_dfm_SLR_10cm_wind_0 = datacatalog.get_geodataframe('dfm_output_event_450_gebco2024_MZB_GTSMv41_CF-0.1_era5_hourly_spw_IBTrACS_CF0_waves')
ds_dfm_SLR_15cm_wind_0 = datacatalog.get_geodataframe('dfm_output_event_450_gebco2024_MZB_GTSMv41_CF-0.15_era5_hourly_spw_IBTrACS_CF0_waves')
ds_dfm_SLR_0cm_wind_1 = datacatalog.get_geodataframe('dfm_output_event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF-1_waves')
ds_dfm_SLR_0cm_wind_5 = datacatalog.get_geodataframe('dfm_output_event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF-5_waves')
ds_dfm_SLR_0cm_wind_10 = datacatalog.get_geodataframe('dfm_output_event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF-10_waves')
ds_dfm_SLR_5cm_wind_1 = datacatalog.get_geodataframe('dfm_output_event_450_gebco2024_MZB_GTSMv41_CF-0.05_era5_hourly_spw_IBTrACS_CF-1_waves')
ds_dfm_SLR_10cm_wind_5 = datacatalog.get_geodataframe('dfm_output_event_450_gebco2024_MZB_GTSMv41_CF-0.1_era5_hourly_spw_IBTrACS_CF-5_waves')
ds_dfm_SLR_15cm_wind_10 = datacatalog.get_geodataframe('dfm_output_event_450_gebco2024_MZB_GTSMv41_CF-0.15_era5_hourly_spw_IBTrACS_CF-10_waves')

#%%
########################################################################################################
########################################## PLOT WATERLEVEL #############################################
########################################################################################################
# Plot timeseries for station 40 under different scenarios for the peak of the storm surge
wl_dfm_SLR_0cm_wind_0 = ds_dfm_SLR_0cm_wind_0["waterlevel"].sel(station=40).compute()

wl_dfm_SLR_5cm_wind_0 = ds_dfm_SLR_5cm_wind_0["waterlevel"].sel(station=40).compute()
wl_dfm_SLR_10cm_wind_0 = ds_dfm_SLR_10cm_wind_0["waterlevel"].sel(station=40).compute()
wl_dfm_SLR_15cm_wind_0 = ds_dfm_SLR_15cm_wind_0["waterlevel"].sel(station=40).compute()

wl_dfm_SLR_0cm_wind_1 = ds_dfm_SLR_0cm_wind_1["waterlevel"].sel(station=40).compute()
wl_dfm_SLR_0cm_wind_5 = ds_dfm_SLR_0cm_wind_5["waterlevel"].sel(station=40).compute()
wl_dfm_SLR_0cm_wind_10 = ds_dfm_SLR_0cm_wind_10["waterlevel"].sel(station=40).compute()

wl_dfm_SLR_5cm_wind_1 = ds_dfm_SLR_5cm_wind_1["waterlevel"].sel(station=40).compute()
wl_dfm_SLR_10cm_wind_5 = ds_dfm_SLR_10cm_wind_5["waterlevel"].sel(station=40).compute()
wl_dfm_SLR_15cm_wind_10 = ds_dfm_SLR_15cm_wind_10["waterlevel"].sel(station=40).compute()


# Actual plotting
fig, ax  = plt.subplots(1, 1, figsize=(12,6))
ax.plot(wl_dfm_SLR_0cm_wind_0["time"], wl_dfm_SLR_0cm_wind_0, label="Waterlevel SLR 0cm Wind 0")

ax.plot(wl_dfm_SLR_5cm_wind_0["time"], wl_dfm_SLR_5cm_wind_0, label="Waterlevel SLR 5cm Wind 0")
ax.plot(wl_dfm_SLR_10cm_wind_0["time"], wl_dfm_SLR_10cm_wind_0, label="Waterlevel SLR 10cm Wind 0")
ax.plot(wl_dfm_SLR_15cm_wind_0["time"], wl_dfm_SLR_15cm_wind_0, label="Waterlevel SLR 15cm Wind 0")

ax.plot(wl_dfm_SLR_0cm_wind_1["time"], wl_dfm_SLR_0cm_wind_1, label="Waterlevel SLR 0cm Wind 1")
ax.plot(wl_dfm_SLR_0cm_wind_5["time"], wl_dfm_SLR_0cm_wind_5, label="Waterlevel SLR 0cm Wind 5")
ax.plot(wl_dfm_SLR_0cm_wind_10["time"], wl_dfm_SLR_0cm_wind_10, label="Waterlevel SLR 0cm Wind 10")

ax.plot(wl_dfm_SLR_5cm_wind_1["time"], wl_dfm_SLR_5cm_wind_1, label="Waterlevel SLR 5cm Wind 1")
ax.plot(wl_dfm_SLR_10cm_wind_5["time"], wl_dfm_SLR_10cm_wind_5, label="Waterlevel SLR 10cm Wind 5")
ax.plot(wl_dfm_SLR_15cm_wind_10["time"], wl_dfm_SLR_15cm_wind_10, label="Waterlevel SLR 15cm Wind 10")


ax.set_xlim([np.datetime64('2019-03-14 18:00'), np.datetime64('2019-03-15 01:00')])
ax.set_ylim([2,3.5])
plt.xlabel("Time")
plt.ylabel("Water Level (m)")
plt.legend()
plt.tight_layout()
plt.show()


#%%
########################################################################################################
############################################# PLOT WIND ################################################
########################################################################################################
# Plot timeseries for wind speed at station 40 under different scenarios for the peak of the storm surge
ds_dfm_SLR_0cm_wind_0_original  = open_ds_his(dir_runs, 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0')
ds_dfm_SLR_0cm_wind_0_05mergefrac  = open_ds_his(dir_runs, 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_0.5mergefrac')
ds_dfm_SLR_0cm_wind_1_original  = open_ds_his(dir_runs, 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF-1')
ds_dfm_SLR_0cm_wind_5_original  = open_ds_his(dir_runs, 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF-5')
ds_dfm_SLR_0cm_wind_10_original = open_ds_his(dir_runs, 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF-10')

wm_dfm_SLR_0cm_wind_0 = ds_dfm_SLR_0cm_wind_0_original["windmag"].sel(station='BEIRA IHO').compute()
wm_dfm_SLR_0cm_wind_0_05mergefrac = ds_dfm_SLR_0cm_wind_0_05mergefrac["windmag"].sel(station='BEIRA IHO').compute()
wm_dfm_SLR_0cm_wind_1 = ds_dfm_SLR_0cm_wind_1_original["windmag"].sel(station='BEIRA IHO').compute()
wm_dfm_SLR_0cm_wind_5 = ds_dfm_SLR_0cm_wind_5_original["windmag"].sel(station='BEIRA IHO').compute()
wm_dfm_SLR_0cm_wind_10 = ds_dfm_SLR_0cm_wind_10_original["windmag"].sel(station='BEIRA IHO').compute()

# Actual plotting
fig, ax  = plt.subplots(1, 1, figsize=(12,6))
ax.plot(wm_dfm_SLR_0cm_wind_0["time"], wm_dfm_SLR_0cm_wind_0, label="wind SLR 0cm Wind 0")

ax.plot(wm_dfm_SLR_0cm_wind_0_05mergefrac["time"], wm_dfm_SLR_0cm_wind_0_05mergefrac, label="wind SLR 0cm Wind 0 0.5 mergefrac")

# ax.plot(wm_dfm_SLR_0cm_wind_1["time"], wm_dfm_SLR_0cm_wind_1, label="Wind SLR 0cm Wind 1")
# ax.plot(wm_dfm_SLR_0cm_wind_5["time"], wm_dfm_SLR_0cm_wind_5, label="Wind SLR 0cm Wind 5")
# ax.plot(wm_dfm_SLR_0cm_wind_10["time"], wm_dfm_SLR_0cm_wind_10, label="Wind SLR 0cm Wind 10")

plt.xlabel("Time")
plt.ylabel("Wind speed (m/s)")
plt.legend()
plt.tight_layout()
plt.show()


#%%
###############################################################################################
######################## SPATIAL PLOTS FOR WIND AND WATER LEVEL ###############################
###############################################################################################
# open map files of Factual and CF-10 wind
model_name = "event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0"
ds_map_F = open_ds_map(dir_runs, model_name)

model_name = "event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF-10"
ds_map_CF10 = open_ds_map(dir_runs, model_name)

#%%
# Plot wind speed spatially
plot_wind_slice(ds_map_F, time_str='2019-03-14 12:00:00')
plot_wind_slice(ds_map_CF10, time_str='2019-03-14 12:00:00')

#%%
# Plot the wind speed difference between two model runs
plot_wind_difference(ds_map_F, ds_map_CF10, time_str='2019-03-14 12:00:00')

#%%
# Calculate the maximum water levels spatially and plot it
compute_and_plot_max_waterlevel(ds_map_F)
compute_and_plot_max_waterlevel(ds_map_CF10)

#%%
# Plot the water level difference between two model runs
plot_waterlevel_difference(ds_map_F, ds_map_CF10)


#%%
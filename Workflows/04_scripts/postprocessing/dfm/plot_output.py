#%% Import the necessary packages using pixi environment compass-snake-dfm
import xarray as xr
import matplotlib.pyplot as plt
import contextily as ctx
import numpy as np
import dfm_tools as dfmt
import os
import glob
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import contextily as ctx
import geopandas as gpd
from shapely.geometry import LineString
import ast
import cartopy.crs as ccrs
import pyproj
from pyproj import Transformer

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
    wind_forcing = "era5_hourly"
    CF_SLR = 0
    CF_SLR_txt = "0"
    CF_wind = 0
    CF_wind_txt = "0"
    dir_runs = f'p:/11210471-001-compass/03_Runs/{region}/{tc_name}/dfm'
    # model_name = f'event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}_nochnk'
    model_name = 'event_450_gebco2024_MZB_GTSMv41_CF0_spw_IBTrACS_CF0'
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


def plot_timeseries_multi(datasets, variable, station_name='BEIRA IHO', labels=None):
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
    ax.legend(loc='upper right', fontsize=8)
    plt.tight_layout()
    plt.show()


def plot_stations_and_grid(ds_his, ds_map, dfm_bbox, crop_bbox_str="[34.8, -20.5, 35.4, -19.3]", crs='EPSG:4326', input_crs="EPSG:32736"):
    crop_bbox = ast.literal_eval(crop_bbox_str)
    transformer = Transformer.from_crs(input_crs, crs, always_xy=True)
    lon, lat = transformer.transform(ds_his.station_x_coordinate.values, ds_his.station_y_coordinate.values)
    ds_his = ds_his.assign_coords({"station_lon": ("station", lon), "station_lat": ("station", lat)})
    
    fig, ax = plt.subplots(1, 2, figsize=(15, 7), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Full domain plot
    ds_map.grid.plot(ax=ax[0], edgecolor='grey', linewidth=0.5, alpha=0.5)
    ctx.add_basemap(ax=ax[0], source=ctx.providers.Esri.WorldImagery, zoom=7, crs=crs, attribution=False)
    ax[0].set_title("Full DFM model domain")
    
    letter_stations = [s for s in ds_his.station.values if not s.isdigit() and s.lower() not in ["airport", "beira off"]]
    ds_his.sel(station=letter_stations).plot.scatter(ax=ax[0], x='station_x_coordinate', y='station_y_coordinate', marker="o")
    for txt in letter_stations:
        ax[0].text(ds_his.station_x_coordinate.sel(station=txt), ds_his.station_y_coordinate.sel(station=txt), txt, size=7, color='white', fontweight='bold')
    ax[0].set_xticks(range(int(dfm_bbox[0]), int(dfm_bbox[2]) + 1, 2))
    ax[0].set_yticks(range(int(dfm_bbox[1]), int(dfm_bbox[3]) + 1, 2))
    ax[0].set_xticklabels([str(i) for i in range(int(dfm_bbox[0]), int(dfm_bbox[2]) + 1, 5)])
    ax[0].set_yticklabels([str(i) for i in range(int(dfm_bbox[1]), int(dfm_bbox[3]) + 1, 5)])
    
    # Zoomed domain
    stations_in_bbox = ds_his.where(
        (ds_his.station_lon >= crop_bbox[0]) & (ds_his.station_lon <= crop_bbox[2]) &
        (ds_his.station_lat >= crop_bbox[1]) & (ds_his.station_lat <= crop_bbox[3]),
        drop=True
    )
    stations_in_bbox_letter = ds_his.where(
        (ds_his.station_x_coordinate >= crop_bbox[0]) & (ds_his.station_x_coordinate <= crop_bbox[2]) &
        (ds_his.station_y_coordinate >= crop_bbox[1]) & (ds_his.station_y_coordinate <= crop_bbox[3]),
        drop=True
    )
    no_airport = [s for s in stations_in_bbox_letter.station.values if s.lower() != "airport"]

    ds_map.grid.plot(ax=ax[1], edgecolor='grey', linewidth=0.5, alpha=0.5)
    ctx.add_basemap(ax=ax[1], source=ctx.providers.Esri.WorldImagery, zoom=9, crs=crs, attribution=False)
    ax[1].plot(ds_his['station_lon'], ds_his['station_lat'], 'xc')
    ax[1].set_xlim(crop_bbox[0], crop_bbox[2])
    ax[1].set_ylim(crop_bbox[1], crop_bbox[3])
    ax[1].set_title("Model domain zoomed into SFINCS bbox")

    stations_in_bbox.plot.scatter(ax=ax[1], x='station_lon', y='station_lat')
    for txt in stations_in_bbox.station.values[::5]:
        ax[1].text(stations_in_bbox.station_lon.sel(station=txt), stations_in_bbox.station_lat.sel(station=txt), txt, size=7, color='white', fontweight='bold')

    stations_in_bbox_letter.plot.scatter(ax=ax[1], x='station_x_coordinate', y='station_y_coordinate', marker="o")
    for txt in no_airport:
        ax[1].text(stations_in_bbox_letter.station_x_coordinate.sel(station=txt), stations_in_bbox_letter.station_y_coordinate.sel(station=txt), txt, size=7, color='white', fontweight='bold')

    ax[1].set_xticks(range(int(crop_bbox[0]), int(crop_bbox[2]) + 1, 2))
    ax[1].set_yticks(range(int(crop_bbox[1]), int(crop_bbox[3]) + 1, 2))
    ax[1].set_xticklabels([str(i) for i in range(int(crop_bbox[0]), int(crop_bbox[2]) + 1, 5)])
    ax[1].set_yticklabels([str(i) for i in range(int(crop_bbox[1]), int(crop_bbox[3]) + 1, 5)])

    for a in ax:
        dfmt.plot_coastlines(ax=a, min_area=1000, linewidth=0.5, zorder=0)
        dfmt.plot_borders(ax=a, zorder=0)
        a.set_xlabel("Longitude")
        a.set_ylabel("Latitude")

    plt.tight_layout()
    return fig, ax


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
    """
    Plot the spatial difference in maximum water level between two datasets over time.

    Automatically caches max values as 'max_<var>' in each dataset.
    """
    max_var = f'max_{var}'

    # Compute max only if not already cached
    if max_var not in ds1:
        ds1[max_var] = ds1[var].max(dim=time_dim, keep_attrs=True)
    if max_var not in ds2:
        ds2[max_var] = ds2[var].max(dim=time_dim, keep_attrs=True)

    # Compute spatial difference
    diff = ds1[max_var] - ds2[max_var]

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    diff.ugrid.plot(ax=ax, transform=ccrs.PlateCarree(), cmap='RdBu_r', center=0, rasterized=True)

    ax.coastlines(resolution='10m')
    ax.set_title(title)
    return ax


def plot_waterlevel_difference(ds1, ds2, var='mesh2d_s1', title='Water level difference', time_dim='time'):
    """
    Plot the spatial difference in maximum water level between two datasets over time.
    
    Parameters:
        ds1, ds2 (xarray.Dataset): Datasets with UGRID structure and water level variable.
        var (str): Variable name for water level.
        title (str): Title for the plot.
        time_dim (str): Name of the time dimension to reduce over.
    """
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

#%% Plot the timeserie for the specified variable and staion (default is BEIRA IHO)
model_name = 'event_450_gebco2024_MZB_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF0'
ds_his = open_ds_his(dir_runs, model_name)
plot_timeserie(ds_his, variable='waterlevel', station_name='1821')
plot_timeserie(ds_his, variable='windmag')


#%%
import hydromt
path_data_cat       = os.path.abspath("../../../03_data_catalogs/datacatalog_SFINCS_coastal_coupling.yml")
datacatalog = hydromt.DataCatalog(data_libs=[path_data_cat])
dfm_run = 'dfm_output_event_450_gebco2024_MZB_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF0_waves'
ds_dfm = datacatalog.get_geodataframe(dfm_run)

#%%
    # tide_surge = ds_combined["tide_surge"].sel(match=station_to_plot).compute()
    # wave_wl = ds_combined["wave_setup"].sel(match=station_to_plot).compute()
waterlevel = ds_dfm["waterlevel"].sel(station_name='1821').compute()

plt.figure(figsize=(12,6))
plt.plot(waterlevel["time"], waterlevel, label="Waterlevel")
    # plt.plot(wave_wl["time"], wave_wl, label="Wave Induced WL")
    # plt.plot(tide_surge["time"], tide_surge, label="Tide + Surge Induced WL")

    # wave_max = ds_combined["wave_setup"].max(dim='time').values[station_idx]

plt.xlabel("Time")
plt.ylabel("Water Level (m)")
# plt.title(f"Water levels for station {station_to_plot} with max wave setup of {wave_max:.2f} m")
plt.legend()
plt.tight_layout()
plt.show()

#%%

test = xr.open_dataset(datacatalog[dfm_run].path)

 #%%
# Compare the timeseries of different models
# model_name = 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_trad'
# ds_his_ERA5_spw_chnk_noOperand = open_ds_his(dir_runs, model_name)

model_name = 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0'
ds_his_ERA5_spw = open_ds_his(dir_runs, model_name)

model_name = 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_CF0'
ds_his_ERA5 = open_ds_his(dir_runs, model_name)

model_name = 'event_450_gebco2024_MZB_GTSMv41_CF0_spw_IBTrACS_CF0'
ds_his_spw = open_ds_his(dir_runs, model_name)

model_name = 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_plusoperand_500km'
ds_his_ERA5_spw_plusoperand_500km = open_ds_his(dir_runs, model_name)

model_name = 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_plusoperand'
ds_his_ERA5_spw_plusoperand = open_ds_his(dir_runs, model_name)

#%%
ds_his_list = [ds_his_ERA5_spw_plusoperand, ds_his_ERA5_spw_plusoperand_500km, ds_his_spw, ds_his_ERA5]
annot = ['ERA5_spw_merged0.75_rad900km', 'ERA5_spw_merged0.75_rad500km', 'spw_only', 'ERA5_only_dynChnk']

plot_timeseries_multi(ds_his_list, variable='waterlevel', labels=annot)
plot_timeseries_multi(ds_his_list, variable='windmag', labels=annot)

 #%%
# Compare the timeseries of different CFs
model_name = 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0'
Fact = open_ds_his(dir_runs, model_name)

model_name = 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF-10'
CF_wind10 = open_ds_his(dir_runs, model_name)

model_name = 'event_450_gebco2024_MZB_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF0'
CF_SLR14cm = open_ds_his(dir_runs, model_name)

model_name = 'event_450_gebco2024_MZB_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10'
CF_wind10_SLR14cm = open_ds_his(dir_runs, model_name)

#%%
ds_his_list = [Fact, CF_wind10, CF_SLR14cm, CF_wind10_SLR14cm]
annot = ['Fact', 'CF_wind10', 'CF_SLR14cm', 'CF_wind10_SLR14cm']

ds_his_list = [CF_SLR14cm, CF_wind10_SLR14cm]
annot = ['CF_SLR14cm', 'CF_wind10_SLR14cm']

plot_timeseries_multi(ds_his_list, variable='waterlevel', labels=annot)
plot_timeseries_multi(ds_his_list, variable='windmag', labels=annot)


#%%
# open map files 
model_name = "event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_plusoperand"
ds_map_ERA5_spw_plusoperand = open_ds_map(dir_runs, model_name)

#%%
model_name = "event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0_plusoperand_500km"
ds_map_ERA5_spw_plusoperand_500km = open_ds_map(dir_runs, model_name)

model_name = 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_CF0' 
ds_map_ERA5 = open_ds_map(dir_runs, model_name)

model_name = 'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0' 
ds_map_ERA5_spw = open_ds_map(dir_runs, model_name)

#%%
model_name = 'event_450_gebco2024_MZB_GTSMv41_CF0_spw_IBTrACS_CF0' 
ds_map_spw = open_ds_map(dir_runs, model_name)

#%%
# Plot wind speed spatially
plot_wind_slice(ds_map_ERA5_spw, time_str='2019-03-14 12:00:00')

#%%
# Plot the wind speed difference between two model runs
plot_wind_difference(ds_map_spw, ds_map_ERA5_spw_plusoperand, time_str='2019-03-14 12:00:00')

#%%
# Calculate the maximum water levels spatially and plot it
compute_and_plot_max_waterlevel(ds_map_ERA5_spw_plusoperand)

compute_and_plot_max_waterlevel(ds_map_spw)


#%%
plot_waterlevel_difference(ds_map_ERA5_spw_plusoperand, ds_map_spw)

#%% Load the DFM model grid for visualisation
grid_ds = xr.open_dataset(os.path.join(dir_runs,model_name,"grid_network.nc"))
dfm_obs = gpd.read_file(dfm_obs_file)

#%% Load the TC track as shapefile
shapefile_path = "p:/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS/IBTrACS_IDAI.shp"
tc_track = gpd.read_file(shapefile_path)

line = LineString(tc_track.geometry)
line_gdf = gpd.GeoDataFrame(geometry=[line], crs=tc_track.crs)

#%%
# Plot spatially the different station and the model grid
plot_stations_and_grid(ds_his_ERA5_spw, ds_map_ERA5_spw, dfm_bbox_model)



#%%
# Load topography/bathymetry data
import xarray as xr
topo_bathy_ds = xr.open_dataset("path_to_topography_or_bathymetry_file.nc")  # Replace with your dataset path
elevation = topo_bathy_ds["elevation"]  # Replace 'elevation' with the variable name in the dataset

# Define a function to plot the background
def add_topo_bathy_background(ax, extent=None):
    if extent:
        elevation_subset = elevation.sel(
            lon=slice(extent[0], extent[2]), lat=slice(extent[1], extent[3])
        )
    else:
        elevation_subset = elevation
    
    # Plot bathymetry and topography as an image
    im = ax.pcolormesh(
        elevation_subset["lon"], elevation_subset["lat"], elevation_subset,
        cmap="terrain", transform=ccrs.PlateCarree(), shading='auto'
    )
    return im

# Plotting
fig, axes = plt.subplots(1, 2, figsize=(15, 7), subplot_kw={'projection': ccrs.PlateCarree()})  # Two subplots side by side

# Define extent for zoomed-in plot (optional)
sfincs_extent = [sfincs_bbox[0], sfincs_bbox[2], sfincs_bbox[1], sfincs_bbox[3]]

# Add topography and bathymetry to both plots
add_topo_bathy_background(axes[0])
im = add_topo_bathy_background(axes[1], extent=sfincs_extent)

# Plot other elements like grid, bounding boxes, and tc_track on top of the background
for ax in axes:
    ax.scatter(filtered_tc_x_coords, filtered_tc_y_coords,
               s=point_sizes, c='orange', alpha=transparencies, edgecolor='black', label="TC Track", transform=ccrs.PlateCarree())
    ax.plot(filtered_tc_x_coords, filtered_tc_y_coords, color="black", linewidth=0.5, label="TC Track Line", transform=ccrs.PlateCarree())

# Grid network points
axes[0].scatter(grid_ds['mesh2d_node_x'].values, grid_ds['mesh2d_node_y'].values,
                s=2, color=(0, 0, 1, 0.5), label="Grid Network", transform=ccrs.PlateCarree())
axes[1].scatter(grid_ds['mesh2d_node_x'].values, grid_ds['mesh2d_node_y'].values,
                s=2, color=(0, 0, 1, 0.5), label="Grid Network", transform=ccrs.PlateCarree())

# Plot bounding boxes
axes[0].plot([lon_min_dfm, lon_max_dfm, lon_max_dfm, lon_min_dfm, lon_min_dfm],
             [lat_min_dfm, lat_min_dfm, lat_max_dfm, lat_max_dfm, lat_min_dfm],
             color="orange", label="BBox 1", transform=ccrs.PlateCarree())
axes[0].plot([sfincs_bbox[0], sfincs_bbox[2], sfincs_bbox[2], sfincs_bbox[0], sfincs_bbox[0]],
             [sfincs_bbox[1], sfincs_bbox[1], sfincs_bbox[3], sfincs_bbox[3], sfincs_bbox[1]],
             color="red", label="BBox 2", transform=ccrs.PlateCarree())

axes[1].plot(x_coords, y_coords, color="yellow", linewidth=2, label="Observation Line", transform=ccrs.PlateCarree())

# Add colorbar for topography/bathymetry
cbar = fig.colorbar(im, ax=axes, orientation="horizontal", fraction=0.046, pad=0.1)
cbar.set_label("Elevation (m)")

# Labels, titles, and layout
axes[0].set_title("Stations with Names (Letters)")
axes[1].set_title("Stations with Numeric Names")
axes[0].set_xlabel("Longitude")
axes[0].set_ylabel("Latitude")
axes[1].set_xlabel("Longitude")
axes[1].set_ylabel("Latitude")

plt.tight_layout()
plt.show()


#%%
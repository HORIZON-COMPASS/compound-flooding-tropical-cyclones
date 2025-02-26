#%% Import the necessary packages
import xarray as xr
import matplotlib.pyplot as plt
import contextily as ctx
import numpy as np
import dfm_tools as dfmt
import os
import glob
from matplotlib.gridspec import GridSpec
#import cartopy.crs as ccrs
import contextily as ctx
import geopandas as gpd
from shapely.geometry import LineString
import ast
import cartopy.crs as ccrs

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
    tidemodel = 'GTSMv41opendap' # tidemodel: FES2014, FES2012, EOT20, GTSMv41, GTSMv41opendap
    wind_forcing = "spw_IBTrACS"
    CF_SLR = 0
    CF_SLR_txt = "0"
    CF_wind = 0
    CF_wind_txt = "0"
    dir_runs = f'p:/11210471-001-compass/03_Runs/{region}/{tc_name}/dfm'
    model = f'event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}'
    dfm_bbox = ast.literal_eval("[32.3,42.5,-27.4,-9.5]")   
    crop_bbox = ast.literal_eval("[34, -20.2, 35.6, -19.2]")
    sfincs_bbox = ast.literal_eval("[34.33,-20.12,34.95,-19.30]")
    dfm_obs_file = "p:/11210471-001-compass/01_Data/Coastal_boundary/points/coastal_bnd_MZB_5mMSL_points_1km.shp"
    
#%% Open the his files in the output folder
for fname in os.listdir(os.path.join(dir_runs,model,'output')):
    if fname.endswith('_his.nc'):
        file_nc_his = os.path.join(dir_runs,model,'output',fname)

#open hisfile with xarray and print netcdf structure
if file_nc_his is not None:
    ds_his = xr.open_mfdataset(file_nc_his, preprocess=dfmt.preprocess_hisnc)

    # locate map files 
file_nc_map = []
for fname in os.listdir(os.path.join(dir_runs,model,'output')):
    if fname.endswith("map.nc"):
        print(fname)
        file_nc_map.append(os.path.join(dir_runs,model,'output',fname))

ds_map = dfmt.open_partitioned_dataset(file_nc_map)

# compute magnitude of wind
ds_map['mesh2d_windmag'] = np.sqrt(ds_map['mesh2d_windx']**2 + ds_map['mesh2d_windy']**2)

#%% Load the DFM model grid for visualisation
grid_ds = xr.open_dataset(os.path.join(dir_runs,model,"grid_network.nc"))
dfm_obs = gpd.read_file(dfm_obs_file)

#%% Load the TC track as shapefile
shapefile_path = "p:/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS/IBTrACS_IDAI.shp"
tc_track = gpd.read_file(shapefile_path)

line = LineString(tc_track.geometry)
line_gdf = gpd.GeoDataFrame(geometry=[line], crs=tc_track.crs)

#%%
dfm_bbox_strp = [float(x) for x in dfm_bbox.strip("[]").split(",")]
lon_min_dfm, lon_max_dfm, lat_min_dfm, lat_max_dfm = dfm_bbox_strp
#%% Plot the stations
fig, ax = plt.subplots(1,1,figsize=(10,5))
# Find the index of the station 'BEIRA IHO'
station_name = 'BEIRA IHO'
station_idx  = ds_his.station.values.tolist().index(station_name)

# Plot the water level for the selected station over time for the different simulations
ax.plot(ds_his.time.values, ds_his.waterlevel[:, station_idx], label='BEIRA IHO - F')

# Set labels and title
ax.set_xlabel('Time')
ax.set_ylabel('Water Level (m)')
ax.set_title(f'Water Level for Station {station_name}')
# Add legend
ax.legend(loc=1, fontsize=8)

# Show the plot
plt.tight_layout()
plt.show()


#%%
# plot net/grid for the whole DFM domain and zoomed into the SFINCS domain
fig, ax = plt.subplots(1, 2, figsize=(15, 7), subplot_kw={'projection': ccrs.PlateCarree()})

pc = ds_map.grid.plot(ax=ax[0], edgecolor='white', linewidth=0.5, alpha=0.5)
ctx.add_basemap(ax=ax[0], source=ctx.providers.Esri.WorldImagery, zoom=7, crs=ccrs.PlateCarree(), attribution=False)
ax[0].plot(ds_his['station_x_coordinate'], ds_his['station_y_coordinate'], 'xc')
ax[0].set_title("Full DFM model domain")

# Set x and y ticks
ax[0].set_xticks(range(int(dfm_bbox[0]), int(dfm_bbox[2])+1, 2))  # Adjust tick interval as needed
ax[0].set_yticks(range(int(dfm_bbox[1]), int(dfm_bbox[3])+1, 2))  # Adjust tick interval as needed
ax[0].set_xticklabels([str(i) for i in range(int(dfm_bbox[0]), int(dfm_bbox[2])+1, 5)])  # Longitude labels
ax[0].set_yticklabels([str(i) for i in range(int(dfm_bbox[1]), int(dfm_bbox[3])+1, 5)])  # Latitude labels

# Plotting the second map with the cropped bbox to SFINCS region
pc2 = ds_map.grid.plot(ax=ax[1], edgecolor='white', linewidth=0.5, alpha=0.5)
ctx.add_basemap(ax=ax[1], source=ctx.providers.Esri.WorldImagery, zoom=9, crs=ccrs.PlateCarree(), attribution=False)
ax[1].plot(ds_his['station_x_coordinate'], ds_his['station_y_coordinate'], 'xc')
ax[1].set_xlim([crop_bbox[0], crop_bbox[2]])
ax[1].set_ylim([crop_bbox[1], crop_bbox[3]])  # Set the extent for the bounding box
ax[1].set_title("Model zomain zoomed into SFINCS bbox")
ax[1].plot([sfincs_bbox[0], sfincs_bbox[2], sfincs_bbox[2], sfincs_bbox[0], sfincs_bbox[0]],
           [sfincs_bbox[1], sfincs_bbox[1], sfincs_bbox[3], sfincs_bbox[3], sfincs_bbox[1]],color="red", label="BBox 1", transform=ccrs.PlateCarree())

# Set x and y ticks for the cropped region
ax[1].set_xticks(range(int(crop_bbox[0]), int(crop_bbox[2])+1, 2))  # Adjust tick interval as needed
ax[1].set_yticks(range(int(crop_bbox[1]), int(crop_bbox[3])+1, 2))  # Adjust tick interval as needed
ax[1].set_xticklabels([str(i) for i in range(int(crop_bbox[0]), int(crop_bbox[2])+1, 2)])  # Longitude labels
ax[1].set_yticklabels([str(i) for i in range(int(crop_bbox[1]), int(crop_bbox[3])+1, 2)])  # Latitude labels

# Show coastlines and borders for both plots
for a in ax:
    dfmt.plot_coastlines(ax=a, min_area=1000, linewidth=0.5, zorder=0)
    dfmt.plot_borders(ax=a, zorder=0)
    # ax.legend(loc='lower right', fontsize=8)
    a.set_xlabel("Longitude")
    a.set_ylabel("Latitude")
    a.set_xticks
    a.set_yticks



#%% 
# Find time of max WL at the output location
id_ts_max = ds_his['waterlevel'].sel(station=['BEIRA IHO']).argmax().values.tolist()
time_ts_max = ds_his['time'].isel(time=id_ts_max).values

#%%
# Plot the TC track and max total water level
fig = plt.figure(layout="constrained",figsize=(12,4))
gs = GridSpec(1, 3, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])#projection = ccrs.epsg(32736))

ax1.set_title('Total water level [mMSL]')
#sc = ax1.scatter(ds_dfm_map['mesh2d_face_x'].values, ds_dfm_map['mesh2d_face_y'].values, c=ds_dfm_map['mesh2d_s1'].sel(time=time_ts_max,method='nearest').values,vmin=-3,vmax=3,  cmap='viridis',s=5)
sc = ds_map['mesh2d_s1'].sel(time=time_ts_max,method='nearest').ugrid.plot(ax=ax1,vmin=-3,vmax=3,cmap='viridis')
plot_loc = ax1.scatter(ds_his.sel(station=['BEIRA IHO']).station_x_coordinate, ds_his.sel(station=['BEIRA IHO']).station_y_coordinate,c='k',label='Model output point')
ax1.set_aspect('equal')
ax1.set_ylim([-22,-19]); ax1.set_xlim([34.2,37])
ctx.add_basemap(ax=ax1, source=ctx.providers.Esri.WorldTopoMap, crs='EPSG:4326', attribution=False)
tc_track.plot(ax=ax1,color='mediumblue', linewidth=1)
line_gdf.plot(ax=ax1, color='mediumblue', linewidth=1, label='TC track Idai (IBTrACS)')
ax1.legend(loc='lower right')

# plot the waterlevel over time at the output point
ax2 = fig.add_subplot(gs[0, 1:])
ax2.plot(ds_his.sel(time=slice('2019-03-11','2019-03-20')).sel(station='BEIRA IHO').time,ds_his.sel(time=slice('2019-03-11','2019-03-20')).sel(station='BEIRA IHO').waterlevel,color='k')
lims = ax2.get_ylim()
ax2.plot([time_ts_max,time_ts_max],lims,'k--')
ax2.set_ylim(lims)
ax2.grid()
ax2.set_title('Timeseries of water levels at model output point')

fig.suptitle('Delft3D-FM model output for tropical cyclone Idai based on IBTrACS')
# %%

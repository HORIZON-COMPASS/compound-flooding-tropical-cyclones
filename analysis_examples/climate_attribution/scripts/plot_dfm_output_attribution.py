# %% imports
import xarray as xr
import matplotlib.pyplot as plt
import contextily as ctx
import numpy as np
import dfm_tools as dfmt
import ast
import cartopy.crs as ccrs

#%%

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
    CF_SLR = -0.14
    CF_SLR_txt = "-0.14"
    CF_wind = -10
    CF_wind_txt = "-10"
    dir_runs = f'p:/11210471-001-compass/03_Runs/{region}/{tc_name}/dfm'
    CF_SLR_model = f'event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF0'
    CF_SLR_wind_model = f'event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}'
    CF_wind_model = f'event_{dfm_res}_{bathy}_{tidemodel}_CF0_{wind_forcing}_CF{CF_wind_txt}'
    F__model = f'event_{dfm_res}_{bathy}_{tidemodel}_CF0_{wind_forcing}_CF0'
    dfm_bbox = ast.literal_eval("[32.3,42.5,-27.4,-9.5]") 
    crop_bbox = ast.literal_eval("[34, -20.5, 35.6, -19.3]")
    sfincs_bbox = ast.literal_eval("[34.33,-20.12,34.95,-19.30]")
    dfm_obs_file = "p:/11210471-001-compass/01_Data/Coastal_boundary/points/coastal_bnd_MZB_5mMSL_points_1km.shp"
    crs = 'EPSG:4326' # coordinate reference system

#%%
# Open the his files in the output folder for the counterfactual (CF) and factual (F) model runs
for fname in os.listdir(os.path.join(dir_runs,CF_SLR_model,'output')):
    if fname.endswith('_his.nc'):
        CF_SLR_file_nc_his = os.path.join(dir_runs,CF_SLR_model,'output',fname)

for fname in os.listdir(os.path.join(dir_runs,CF_SLR_wind_model,'output')):
    if fname.endswith('_his.nc'):
        CF_SLR_wind_file_nc_his = os.path.join(dir_runs,CF_SLR_wind_model,'output',fname)

for fname in os.listdir(os.path.join(dir_runs,CF_wind_model,'output')):
    if fname.endswith('_his.nc'):
        CF_wind_file_nc_his = os.path.join(dir_runs,CF_wind_model,'output',fname)

for fname in os.listdir(os.path.join(dir_runs,F__model,'output')):
    if fname.endswith('_his.nc'):
        F__file_nc_his = os.path.join(dir_runs,F__model,'output',fname)
        
# open hisfile with xarray and print netcdf structure
if CF_SLR_file_nc_his is not None:
    ds_his_CF_SLR = xr.open_mfdataset(CF_SLR_file_nc_his, preprocess=dfmt.preprocess_hisnc)

if CF_SLR_wind_file_nc_his is not None:
    ds_his_CF_SLR_wind = xr.open_mfdataset(CF_SLR_wind_file_nc_his, preprocess=dfmt.preprocess_hisnc)

if CF_wind_file_nc_his is not None:
    ds_his_CF_wind = xr.open_mfdataset(CF_wind_file_nc_his, preprocess=dfmt.preprocess_hisnc)

if F__file_nc_his is not None:
    ds_his_F = xr.open_mfdataset(F__file_nc_his, preprocess=dfmt.preprocess_hisnc)

#%%
# locate map files 
file_nc_map_F = []
for fname in os.listdir(os.path.join(dir_runs,F__model,'output')):
    if fname.endswith("map.nc"):
        print(fname)
        file_nc_map_F.append(os.path.join(dir_runs,F__model,'output',fname))

file_nc_map_CF_SLR_wind = []
for fname in os.listdir(os.path.join(dir_runs,CF_SLR_wind_model,'output')):
    if fname.endswith("map.nc"):
        print(fname)
        file_nc_map_CF_SLR_wind.append(os.path.join(dir_runs,CF_SLR_wind_model,'output',fname))

file_nc_map_CF_wind = []
for fname in os.listdir(os.path.join(dir_runs,CF_wind_model,'output')):
    if fname.endswith("map.nc"):
        print(fname)
        file_nc_map_CF_wind.append(os.path.join(dir_runs,CF_wind_model,'output',fname))

file_nc_map_CF_SLR = []
for fname in os.listdir(os.path.join(dir_runs,CF_SLR_model,'output')):
    if fname.endswith("map.nc"):
        print(fname)
        file_nc_map_CF_SLR.append(os.path.join(dir_runs,CF_SLR_model,'output',fname))

ds_map_F           = dfmt.open_partitioned_dataset(file_nc_map_F)
ds_map_CF_SLR_wind = dfmt.open_partitioned_dataset(file_nc_map_CF_SLR_wind)
ds_map_CF_wind     = dfmt.open_partitioned_dataset(file_nc_map_CF_wind)
ds_map_CF_SLR      = dfmt.open_partitioned_dataset(file_nc_map_CF_SLR)

#%%
# get and print variable properties
# vars_pd_F  = dfmt.get_ncvarproperties(ds_his_F)
# vars_pd_CF = dfmt.get_ncvarproperties(ds_his_CF_SLR)
# print(vars_pd_F)
# print(vars_pd_CF)

#%%
# plot his data: waterlevel at the BEIRA IHO station
fig, ax = plt.subplots(1,1,figsize=(10,5))
# Find the index of the station 'BEIRA IHO'
station_name = 'BEIRA IHO'
station_idx_F           = ds_his_F.station.values.tolist().index(station_name)
station_idx_CF_SLR      = ds_his_CF_SLR.station.values.tolist().index(station_name)
station_idx_CF_wind     = ds_his_CF_wind.station.values.tolist().index(station_name)
station_idx_CF_SLR_wind = ds_his_CF_SLR_wind.station.values.tolist().index(station_name)

# Plot the water level for the selected station over time for the different simulations
ax.plot(ds_his_F.time.values, ds_his_F.waterlevel[:, station_idx_F], label='BEIRA IHO - F')
ax.plot(ds_his_CF_SLR.time.values, ds_his_CF_SLR.waterlevel[:, station_idx_CF_SLR], label='BEIRA IHO - CF_SLR -0.14m')
ax.plot(ds_his_CF_wind.time.values, ds_his_CF_wind.waterlevel[:, station_idx_CF_wind], label='BEIRA IHO - CF_wind -10%')
ax.plot(ds_his_CF_SLR_wind.time.values, ds_his_CF_SLR_wind.waterlevel[:, station_idx_CF_SLR_wind], label='BEIRA IHO - CF_SLR -0.14m & CF_wind -10%')

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
crop_bbox = ast.literal_eval("[34, -20.5, 35.6, -19]")
# plot net/grid for the whole DFM domain and zoomed into the SFINCS domain
fig, ax = plt.subplots(1, 2, figsize=(15, 7), subplot_kw={'projection': ccrs.PlateCarree()})

pc = ds_map_F.grid.plot(ax=ax[0], edgecolor='white', linewidth=0.5, alpha=0.5)
ctx.add_basemap(ax=ax[0], source=ctx.providers.Esri.WorldImagery, zoom=7, crs=crs, attribution=False)
ax[0].plot(ds_his_F['station_x_coordinate'], ds_his_F['station_y_coordinate'], 'xc')
ax[0].set_title("Full DFM model domain")

# Set x and y ticks
ax[0].set_xticks(range(int(dfm_bbox[0]), int(dfm_bbox[2])+1, 2))  # Adjust tick interval as needed
ax[0].set_yticks(range(int(dfm_bbox[1]), int(dfm_bbox[3])+1, 2))  # Adjust tick interval as needed
ax[0].set_xticklabels([str(i) for i in range(int(dfm_bbox[0]), int(dfm_bbox[2])+1, 5)])  # Longitude labels
ax[0].set_yticklabels([str(i) for i in range(int(dfm_bbox[1]), int(dfm_bbox[3])+1, 5)])  # Latitude labels

# Plotting the second map with the cropped bbox to SFINCS region
pc2 = ds_map_F.grid.plot(ax=ax[1], edgecolor='white', linewidth=0.5, alpha=0.5)
ctx.add_basemap(ax=ax[1], source=ctx.providers.Esri.WorldImagery, zoom=9, crs=crs, attribution=False)
ax[1].plot(ds_his_F['station_x_coordinate'], ds_his_F['station_y_coordinate'], 'xc')
ax[1].set_xlim([crop_bbox[0], crop_bbox[2]])
ax[1].set_ylim([crop_bbox[1], crop_bbox[3]])  # Set the extent for the bounding box
ax[1].set_title("Model zomain zoomed into SFINCS bbox")
ax[1].plot([sfincs_bbox[0], sfincs_bbox[2], sfincs_bbox[2], sfincs_bbox[0], sfincs_bbox[0]],
           [sfincs_bbox[1], sfincs_bbox[1], sfincs_bbox[3], sfincs_bbox[3], sfincs_bbox[1]],color="red", label="BBox 1", transform=ccrs.PlateCarree())

# Set x and y ticks for the cropped region
ax[1].set_xticks(range(int(crop_bbox[0]), int(crop_bbox[2])+1, 2))  # Adjust tick interval as needed
ax[1].set_yticks(range(int(crop_bbox[1]), int(crop_bbox[3])+1, 2))  # Adjust tick interval as needed
ax[1].set_xticklabels([str(i) for i in range(int(crop_bbox[0]), int(crop_bbox[2])+1, 5)])  # Longitude labels
ax[1].set_yticklabels([str(i) for i in range(int(crop_bbox[1]), int(crop_bbox[3])+1, 5)])  # Latitude labels

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
# Compute magnitude of wind by combining X and Y wind speeds
ds_his_F['windmag'] = np.sqrt(ds_his_F['windx']**2 + ds_his_F['windy']**2)
ds_map_F['mesh2d_windmag'] = np.sqrt(ds_map_F['mesh2d_windx']**2 + ds_map_F['mesh2d_windy']**2)

ds_his_CF_SLR_wind['windmag'] = np.sqrt(ds_his_CF_SLR_wind['windx']**2 + ds_his_CF_SLR_wind['windy']**2)
ds_map_CF_SLR_wind['mesh2d_windmag'] = np.sqrt(ds_map_CF_SLR_wind['mesh2d_windx']**2 + ds_map_CF_SLR_wind['mesh2d_windy']**2)

ds_his_CF_wind['windmag'] = np.sqrt(ds_his_CF_wind['windx']**2 + ds_his_CF_wind['windy']**2)
ds_map_CF_wind['mesh2d_windmag'] = np.sqrt(ds_map_CF_wind['mesh2d_windx']**2 + ds_map_CF_wind['mesh2d_windy']**2)

ds_his_CF_SLR['windmag'] = np.sqrt(ds_his_CF_SLR['windx']**2 + ds_his_CF_SLR['windy']**2)
ds_map_CF_SLR['mesh2d_windmag'] = np.sqrt(ds_map_CF_SLR['mesh2d_windx']**2 + ds_map_F['mesh2d_windy']**2)

# Find time of max wind at the BEIRA IHO station
id_ts_max = ds_his_F['windmag'].sel(station=['BEIRA IHO']).argmax().values.tolist()
time_ts_max = ds_his_F['time'].isel(time=id_ts_max).values

# id_ts_max_CF = ds_his_CF_SLR_wind['windmag'].sel(station=['BEIRA IHO']).argmax().values.tolist()
# time_ts_max_CF = ds_his_CF_SLR_wind['time'].isel(time=id_ts_max_CF).values
#%%
# Plot the max wind speed at BEIRA IHO for the Factual simulation
fig, ax = plt.subplots(1, 1, figsize=(15, 7), subplot_kw={'projection': ccrs.PlateCarree()})

sc = ds_map_F['mesh2d_windmag'].sel(time=time_ts_max,method='nearest').ugrid.plot(ax=ax,cmap='viridis')
plot_loc = ax.scatter(ds_his_F.sel(station=['BEIRA IHO']).station_x_coordinate, ds_his_F.sel(station=['BEIRA IHO']).station_y_coordinate,c='red',s=5,label='BEIRA IHO')
ctx.add_basemap(ax=ax, source=ctx.providers.Esri.WorldTopoMap, crs=ccrs.PlateCarree(), attribution=False)

# Add lat/lon gridlines and labels
gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5, alpha=0.7)
gl.right_labels = False  # Hide labels on the right
gl.top_labels = False  # Hide labels on top
# add legend
ax.legend(loc='lower right', fontsize=8)
# Access the colorbar and set label
cbar = sc.colorbar
cbar.set_label("Wind speed (m/s)", fontsize=9)
# Set title for the specific axis (recommended for single plots)
ax.set_title(f"Maximum wind speed at BEIRA IHO at {time_ts_max}", fontsize=10)

#%%
# Plot the max wind speed at BEIRA IHO for the F & CF simulation, and the difference
fig, axes = plt.subplots(1, 3, figsize=(18, 7), subplot_kw={'projection': ccrs.PlateCarree()})

# Original Wind Speed Plot
sc1 = ds_map_F['mesh2d_windmag'].sel(time=time_ts_max, method='nearest').ugrid.plot(ax=axes[0], cmap='viridis')
axes[0].scatter(ds_his_F.sel(station=['BEIRA IHO']).station_x_coordinate, 
                ds_his_F.sel(station=['BEIRA IHO']).station_y_coordinate, 
                c='red', s=5, label='BEIRA IHO')
ctx.add_basemap(axes[0], source=ctx.providers.Esri.WorldTopoMap, crs=ccrs.PlateCarree(), attribution=False)
axes[0].set_title(f"Maximum F Wind Speed at BEIRA IHO at {time_ts_max}", fontsize=10)
axes[0].legend(loc='lower right', fontsize=8)

sc2 = ds_map_CF_wind['mesh2d_windmag'].sel(time=time_ts_max, method='nearest').ugrid.plot(ax=axes[1], cmap='viridis')
axes[1].scatter(ds_his_F.sel(station=['BEIRA IHO']).station_x_coordinate, 
                ds_his_F.sel(station=['BEIRA IHO']).station_y_coordinate, 
                c='red', s=5, label='BEIRA IHO')
ctx.add_basemap(axes[1], source=ctx.providers.Esri.WorldTopoMap, crs=ccrs.PlateCarree(), attribution=False)
axes[1].set_title(f"Maximum CF Wind Speed at BEIRA IHO at {time_ts_max}", fontsize=10)
axes[1].legend(loc='lower right', fontsize=8)

# Difference Plot (ds_map_F - ds_map_CF_SLR_wind)
wind_diff = ds_map_F['mesh2d_windmag'].sel(time=time_ts_max, method='nearest') - \
            ds_map_CF_wind['mesh2d_windmag'].sel(time=time_ts_max, method='nearest')
vmax = wind_diff.max().compute()  # Computes max() first
sc3 = wind_diff.ugrid.plot(ax=axes[2], cmap='RdBu_r', vmin=-vmax, vmax=vmax)
axes[2].scatter(ds_his_F.sel(station=['BEIRA IHO']).station_x_coordinate, 
                ds_his_F.sel(station=['BEIRA IHO']).station_y_coordinate, 
                c='red', s=5, label='BEIRA IHO')
ctx.add_basemap(axes[2], source=ctx.providers.Esri.WorldTopoMap, crs=ccrs.PlateCarree(), attribution=False)
axes[2].set_title(f"Max wind Speed Difference at {time_ts_max}", fontsize=10)
axes[2].legend(loc='lower right', fontsize=8)

# Add Gridlines to both subplots
for ax in axes:
    gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5, alpha=0.7)
    gl.right_labels = False  
    gl.top_labels = False  

# Set Colorbar Labels
sc1.colorbar.set_label("Wind Speed (m/s)", fontsize=9)
sc2.colorbar.set_label("Wind Speed (m/s)", fontsize=9)
sc3.colorbar.set_label("Wind Speed Difference (m/s)", fontsize=9)

# Set a figure-wide title
fig.suptitle("Max wind Speed at BEIRA IHO", fontsize=14)

plt.tight_layout()
plt.show()

# %%
# Plot the wind magnitude over time
fig, ax = plt.subplots(1,1,figsize=(10,5))

# Plot the water level for the selected station over time
ax.plot(ds_his_F.time.values, ds_his_F.windmag[:, station_idx_F], label='BEIRA IHO - F')
ax.plot(ds_his_CF_SLR.time.values, ds_his_CF_SLR.windmag[:, station_idx_CF_SLR], label='BEIRA IHO - CF_SLR -0.14m')
ax.plot(ds_his_CF_wind.time.values, ds_his_CF_wind.windmag[:, station_idx_CF_SLR_wind], label='BEIRA IHO - CF_wind -10%')
ax.plot(ds_his_CF_SLR_wind.time.values, ds_his_CF_SLR_wind.windmag[:, station_idx_CF_SLR_wind], label='BEIRA IHO - CF_SLR -0.14m & CF_wind -10%')

# Set labels and title
ax.set_xlabel('Time')
ax.set_ylabel('Wind Speed (m/s)')
ax.set_title(f'Wind speed at {station_name}')

# Add legend
ax.legend(loc=1, fontsize=8)

# Show the plot
plt.tight_layout()
plt.show()


# %%
# Compare coordinate values
print(ds_map_F['mesh2d_windmag'].shape)
print(ds_map_CF_SLR_wind['mesh2d_windmag'].shape)

# %%
print(ds_map_F['mesh2d_windmag'].sel(time=time_ts_max, method='nearest').time.values)
print(ds_map_CF_SLR_wind['mesh2d_windmag'].sel(time=time_ts_max, method='nearest').time.values)

# %%

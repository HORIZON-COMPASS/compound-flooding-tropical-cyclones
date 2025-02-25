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
    tidemodel = 'FES2014' # tidemodel: FES2014, FES2012, EOT20, GTSMv41, GTSMv41opendap
    wind_forcing = "spw_IBTrACS"
    CF_SLR = -0.14
    CF_SLR_txt = "-0.14"
    CF_wind = 0
    CF_wind_txt = "0"
    dir_runs = f'p:/11210471-001-compass/03_Runs/{region}/{tc_name}/dfm'
    CF_model = f'event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}'
    F__model = f'event_{dfm_res}_{bathy}_{tidemodel}_CF0_{wind_forcing}_CF0'
    dfm_bbox = ast.literal_eval("[32.3,42.5,-27.4,-9.5]") 
    crop_bbox = ast.literal_eval("[34, -20.5, 35.6, -19.5]")
    sfincs_bbox = ast.literal_eval("[34.33,-20.12,34.95,-19.30]")
    dfm_obs_file = "p:/11210471-001-compass/01_Data/Coastal_boundary/points/coastal_bnd_MZB_5mMSL_points_1km.shp"
    crs = 'EPSG:4326' # coordinate reference system

#%%
# Open the his files in the output folder for the counterfactual (CF) and factual (F) model runs
for fname in os.listdir(os.path.join(dir_runs,CF_model,'output')):
    if fname.endswith('_his.nc'):
        CF_file_nc_his = os.path.join(dir_runs,CF_model,'output',fname)

for fname in os.listdir(os.path.join(dir_runs,F__model,'output')):
    if fname.endswith('_his.nc'):
        F__file_nc_his = os.path.join(dir_runs,F__model,'output',fname)

# open hisfile with xarray and print netcdf structure
if CF_file_nc_his is not None:
    ds_his_CF = xr.open_mfdataset(CF_file_nc_his, preprocess=dfmt.preprocess_hisnc)

# open hisfile with xarray and print netcdf structure
if F__file_nc_his is not None:
    ds_his_F = xr.open_mfdataset(F__file_nc_his, preprocess=dfmt.preprocess_hisnc)

print(ds_his_F)
print(ds_his_CF)

#%%
# get and print variable properties
vars_pd_F  = dfmt.get_ncvarproperties(ds_his_F)
vars_pd_CF = dfmt.get_ncvarproperties(ds_his_CF)
print(vars_pd_F)
print(vars_pd_CF)

#%%
# plot his data: waterlevel at stations
fig, ax = plt.subplots(1,1,figsize=(10,5))
# Find the index of the station 'BEIRA IHO'
station_name = 'BEIRA IHO'
station_idx_F  = ds_his_F.station.values.tolist().index(station_name)
station_idx_CF = ds_his_CF.station.values.tolist().index(station_name)

# Plot the water level for the selected station over time
ax.plot(ds_his_F.time.values, ds_his_F.waterlevel[:, station_idx_F], label='BEIRA IHO - F')
ax.plot(ds_his_CF.time.values, ds_his_CF.waterlevel[:, station_idx_CF], label='BEIRA IHO - CF_SLR -0.14m')

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
# plot his data: waterlevel at stations
fig, ax = plt.subplots(1,1,figsize=(10,5))
ds_his['BEIRA IHO'].waterlevel.plot.line(ax=ax, x='time')
# ax.legend(ds_his.station.to_series(), loc=1, fontsize=8)

#%%
    # locate map files 
file_nc_map = []
for fname in os.listdir(os.path.join(dir_runs,model,'output')):
    if fname.endswith("map.nc"):
        print(fname)
        file_nc_map.append(os.path.join(dir_runs,model,'output',fname))

ds_map = dfmt.open_partitioned_dataset(file_nc_map)
ds_map
#%%
# get and print variable properties
vars_pd = dfmt.get_ncvarproperties(ds_map)
vars_pd.head(10)

#%%
# plot net/grid. use random variable and plot line to get grid
fig, ax = plt.subplots(1, 2, figsize=(15, 7), subplot_kw={'projection': ccrs.PlateCarree()})

pc = ds_map.grid.plot(ax=ax[0], edgecolor='white', linewidth=0.5, alpha=0.5)
ctx.add_basemap(ax=ax[0], source=ctx.providers.Esri.WorldImagery, zoom=7, crs=crs, attribution=False)
ax[0].plot(ds_his['station_x_coordinate'], ds_his['station_y_coordinate'], 'xc')
ax[0].set_title("Full DFM model domain")

# Set x and y ticks for the cropped region
ax[0].set_xticks(range(int(dfm_bbox[0]), int(dfm_bbox[2])+1, 2))  # Adjust tick interval as needed
ax[0].set_yticks(range(int(dfm_bbox[1]), int(dfm_bbox[3])+1, 2))  # Adjust tick interval as needed
ax[0].set_xticklabels([str(i) for i in range(int(dfm_bbox[0]), int(dfm_bbox[2])+1, 5)])  # Longitude labels
ax[0].set_yticklabels([str(i) for i in range(int(dfm_bbox[1]), int(dfm_bbox[3])+1, 5)])  # Latitude labels

# Plotting the second map with the defined bbox
pc2 = ds_map.grid.plot(ax=ax[1], edgecolor='white', linewidth=0.5, alpha=0.5)
ctx.add_basemap(ax=ax[1], source=ctx.providers.Esri.WorldImagery, zoom=9, crs=crs, attribution=False)
ax[1].plot(ds_his['station_x_coordinate'], ds_his['station_y_coordinate'], 'xc')
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
# compute magnitude of wind
ds_map['mesh2d_windmag'] = np.sqrt(ds_map['mesh2d_windx']**2 + ds_map['mesh2d_windy']**2)
# %%

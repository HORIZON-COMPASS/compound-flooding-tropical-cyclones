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
    bathy = "gebco2024"
    tidemodel = 'GTSMv41opendap' # tidemodel: FES2014, FES2012, EOT20, GTSMv41, GTSMv41opendap
    wind_forcing = "spw_IBTrACS_ext"
    dir_runs = f'p:/11210471-001-compass/03_Runs/{region}/{tc_name}/dfm'
    model = f'event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}'
    dfm_bbox = "[32.3,42.5,-27.4,-9.5]"   
    crop_bbox = ast.literal_eval("[34, -20.5, 35.6, -19.5]")
    sfincs_bbox = ast.literal_eval("[34.33,-20.12,34.95,-19.30]")
    dfm_obs_file = "p:/11210471-001-compass/01_Data/Coastal_boundary/points/coastal_bnd_MZB_5mMSL_points_1km.shp"
    
#%% Open the his files in the output folder
for fname in os.listdir(os.path.join(dir_runs,model,'output')):
    if fname.endswith('_his.nc'):
        file_nc_his = os.path.join(dir_runs,model,'output',fname)

#open hisfile with xarray and print netcdf structure
if file_nc_his is not None:
    ds_his = xr.open_mfdataset(file_nc_his, preprocess=dfmt.preprocess_hisnc)

#%% Load the DFM model grid for visualisation
grid_ds = xr.open_dataset(os.path.join(dir_runs,model,"grid_network.nc"))
dfm_obs = gpd.read_file(dfm_obs_file)

#%% Load the TC track as shapefile
shapefile_path = "p:/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS/IBTrACS_IDAI.shp"
tc_track = gpd.read_file(shapefile_path)

#%%
dfm_bbox_strp = [float(x) for x in dfm_bbox.strip("[]").split(",")]
lon_min_dfm, lon_max_dfm, lat_min_dfm, lat_max_dfm = dfm_bbox_strp
#%% Plot the stations
# Extract the coordinates from the 'geometry' column
x_coords = [point.x for point in dfm_obs.geometry if point is not None]
y_coords = [point.y for point in dfm_obs.geometry if point is not None]

# Plotting
fig, axes = plt.subplots(1, 2, figsize=(15, 7), subplot_kw={'projection': ccrs.PlateCarree()})  # Two subplots side by side

# Plot 1: Stations with letter names
letter_stations = [s for s in ds_his.station.values if not s.isdigit()]
ds_his.sel(station=letter_stations).plot.scatter(ax=axes[0], x='station_x_coordinate', y='station_y_coordinate', marker="o")
for tt, txt in enumerate(letter_stations):
    axes[0].text(ds_his.station_x_coordinate.sel(station=txt), ds_his.station_y_coordinate.sel(station=txt), txt, size=7)
axes[0].set_title("Stations with Names (Letters)")
axes[0].scatter(grid_ds['mesh2d_node_x'].values, grid_ds['mesh2d_node_y'].values, s=2, color="blue", label="Grid Network", transform=ccrs.PlateCarree())
axes[0].plot([lon_min_dfm, lon_max_dfm, lon_max_dfm, lon_min_dfm, lon_min_dfm],
             [lat_min_dfm, lat_min_dfm, lat_max_dfm, lat_max_dfm, lat_min_dfm],
             color="orange", label="BBox 1", transform=ccrs.PlateCarree())
axes[0].plot([sfincs_bbox[0], sfincs_bbox[2], sfincs_bbox[2], sfincs_bbox[0], sfincs_bbox[0]],
             [sfincs_bbox[1], sfincs_bbox[1], sfincs_bbox[3], sfincs_bbox[3], sfincs_bbox[1]],color="red", label="BBox 1", transform=ccrs.PlateCarree())

# Plot 2: Stations with numeric names at the SFINCS bbox
numeric_stations = [s for s in ds_his.station.values if s.isdigit()]
ds_his.sel(station=numeric_stations).plot.scatter(ax=axes[1], x='station_x_coordinate', y='station_y_coordinate', marker="o")
for tt, txt in enumerate(numeric_stations):
    axes[1].text(ds_his.station_x_coordinate.sel(station=txt), ds_his.station_y_coordinate.sel(station=txt), txt, size=7)
axes[1].set_title("Stations with Numeric Names")
axes[1].scatter(grid_ds['mesh2d_node_x'].values, grid_ds['mesh2d_node_y'].values, s=2, color="blue", label="Grid Network", transform=ccrs.PlateCarree())
axes[1].plot([sfincs_bbox[0], sfincs_bbox[2], sfincs_bbox[2], sfincs_bbox[0], sfincs_bbox[0]],
             [sfincs_bbox[1], sfincs_bbox[1], sfincs_bbox[3], sfincs_bbox[3], sfincs_bbox[1]],color="red", label="BBox 1", transform=ccrs.PlateCarree())

# Plot the observation line (connecting the points from dfm_obs)
axes[1].plot([point.x for point in dfm_obs.geometry if point is not None], [point.y for point in dfm_obs.geometry if point is not None], color="yellow", linewidth=2, label="Observation Line", transform=ccrs.PlateCarree())

# Set limits for Plot 2 (zooming into the SFINCS bbox)
axes[1].set_xlim([crop_bbox[0], crop_bbox[2]])
axes[1].set_ylim([crop_bbox[1], crop_bbox[3]])

# Show coastlines and borders for both plots
for ax in axes:
    dfmt.plot_coastlines(ax=ax, min_area=1000, linewidth=0.5, zorder=0)
    dfmt.plot_borders(ax=ax, zorder=0)

# Display the plot
plt.tight_layout()
plt.show()

#%%
# Extract x and y coordinates from the 'geometry' column of dfm_obs
x_coords = [point.x for point in dfm_obs.geometry if point is not None]
y_coords = [point.y for point in dfm_obs.geometry if point is not None]

# Plotting
fig, axes = plt.subplots(1, 2, figsize=(15, 7), subplot_kw={'projection': ccrs.PlateCarree()})  # Two subplots side by side

# Plot 1: Stations with letter names
letter_stations = [s for s in ds_his.station.values if not s.isdigit()]
ds_his.sel(station=letter_stations).plot.scatter(ax=axes[0], x='station_x_coordinate', y='station_y_coordinate', marker="o")
for tt, txt in enumerate(letter_stations):
    axes[0].text(ds_his.station_x_coordinate.sel(station=txt), ds_his.station_y_coordinate.sel(station=txt), txt, size=7)
axes[0].set_title("Stations with Names (Letters)")

# Grid network points in a nicer blue with transparency
axes[0].scatter(grid_ds['mesh2d_node_x'].values, grid_ds['mesh2d_node_y'].values, 
                s=2, color=(0, 0, 1, 0.5), label="Grid Network", transform=ccrs.PlateCarree())

# Add background features
axes[0].set_facecolor('lightblue')
axes[0].coastlines(resolution='50m', color='black', linewidth=1)
axes[0].add_feature(cfeature.LAND, facecolor='lightgray', zorder=1)
axes[0].add_feature(cfeature.OCEAN, facecolor='lightblue', zorder=0)

# Plot bounding boxes
axes[0].plot([lon_min_dfm, lon_max_dfm, lon_max_dfm, lon_min_dfm, lon_min_dfm],
             [lat_min_dfm, lat_min_dfm, lat_max_dfm, lat_max_dfm, lat_min_dfm],
             color="orange", label="BBox 1", transform=ccrs.PlateCarree())
axes[0].plot([sfincs_bbox[0], sfincs_bbox[2], sfincs_bbox[2], sfincs_bbox[0], sfincs_bbox[0]],
             [sfincs_bbox[1], sfincs_bbox[1], sfincs_bbox[3], sfincs_bbox[3], sfincs_bbox[1]],
             color="red", label="BBox 2", transform=ccrs.PlateCarree())

axes[0].set_xlabel("Longitude")  # Add x-axis label
axes[0].set_ylabel("Latitude")  # Add y-axis label

# Plot 2: Stations with numeric names at the SFINCS bbox
numeric_stations = [s for s in ds_his.station.values if s.isdigit()]
ds_his.sel(station=numeric_stations).plot.scatter(ax=axes[1], x='station_x_coordinate', y='station_y_coordinate', marker="o")
for tt, txt in enumerate(numeric_stations):
    axes[1].text(ds_his.station_x_coordinate.sel(station=txt), ds_his.station_y_coordinate.sel(station=txt), txt, size=7)
axes[1].set_title("Stations with Numeric Names")

# Grid network points in a nicer blue with transparency
axes[1].scatter(grid_ds['mesh2d_node_x'].values, grid_ds['mesh2d_node_y'].values, 
                s=2, color=(0, 0, 1, 0.5), label="Grid Network", transform=ccrs.PlateCarree())

# Plot the observation line (connecting the points from dfm_obs)
axes[1].plot(x_coords, y_coords, color="yellow", linewidth=2, label="Observation Line", transform=ccrs.PlateCarree())

# Set limits for Plot 2 (zoom into the SFINCS bbox)
axes[1].set_xlim([sfincs_bbox[0], sfincs_bbox[2]])
axes[1].set_ylim([sfincs_bbox[1], sfincs_bbox[3]])

# Add background features
axes[1].set_facecolor('lightblue')
axes[1].coastlines(resolution='50m', color='black', linewidth=1)
axes[1].add_feature(cfeature.LAND, facecolor='lightgray', zorder=1)
axes[1].add_feature(cfeature.OCEAN, facecolor='lightblue', zorder=0)

axes[1].set_xlabel("Longitude")  # Add x-axis label
axes[1].set_ylabel("Latitude")  # Add y-axis label

# Show coastlines and borders for both plots
for ax in axes:
    dfmt.plot_coastlines(ax=ax, min_area=1000, linewidth=0.5, zorder=0)
    dfmt.plot_borders(ax=ax, zorder=0)

# Display the plot
plt.tight_layout()
plt.show()

#%%
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from cartopy.io.shapereader import natural_earth

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










#%% Plot the water level at the Beira IHO and offshore station, and the obs_points
if file_nc_his is not None:
    fig, ax = plt.subplots(1,1,figsize=(10,5))
    ds_his_sel.sel(stations=['BEIRA IHO']).waterlevel.plot.line(ax=ax, x='time')
    ax.legend(ds_his_sel.stations.to_series(),bbox_to_anchor=(1.04, 1),loc="upper left",fontsize=8) 
    plt.grid()

#%% Open the produced waterlevel maps
files_nc_map = glob.glob(os.path.join(dir_runs,model,'output','*map.nc'))
ds_dfm_map = dfmt.open_partitioned_dataset(files_nc_map)

# Find time of max WL at the output location
id_ts_max = ds_his_sel['waterlevel'].sel(stations=['BEIRA IHO']).argmax().values.tolist()
time_ts_max = ds_his_sel['time'].isel(time=id_ts_max).values

#%%

# Check the contents of the GeoDataFrame
print(gdf.head())
line = LineString(gdf.geometry)
line_gdf = gpd.GeoDataFrame(geometry=[line], crs=gdf.crs)

# Plot the TC track and max total water level
fig = plt.figure(layout="constrained",figsize=(12,4))
gs = GridSpec(1, 3, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])#projection = ccrs.epsg(32736))

ax1.set_title('Total water level [mMSL]')
#sc = ax1.scatter(ds_dfm_map['mesh2d_face_x'].values, ds_dfm_map['mesh2d_face_y'].values, c=ds_dfm_map['mesh2d_s1'].sel(time=time_ts_max,method='nearest').values,vmin=-3,vmax=3,  cmap='viridis',s=5)
sc = ds_dfm_map['mesh2d_s1'].sel(time=time_ts_max,method='nearest').ugrid.plot(ax=ax1,vmin=-3,vmax=3,cmap='viridis')
plot_loc = ax1.scatter(ds_his_sel.sel(stations=['BEIRA IHO']).station_x_coordinate, ds_his_sel.sel(stations=['BEIRA IHO']).station_y_coordinate,c='k',label='Model output point')
ax1.set_aspect('equal')
ax1.set_ylim([-22,-19]); ax1.set_xlim([34.2,37])
ctx.add_basemap(ax=ax1, source=ctx.providers.Esri.WorldTopoMap, crs='EPSG:4326', attribution=False)
gdf.plot(ax=ax1,color='mediumblue', linewidth=1)
line_gdf.plot(ax=ax1, color='mediumblue', linewidth=1, label='TC track Idai (IBTrACS)')
ax1.legend(loc='lower right')

# plot the waterlevel over time at the output point
ax2 = fig.add_subplot(gs[0, 1:])
ax2.plot(ds_his_sel.sel(time=slice('2019-03-11','2019-03-20')).sel(stations='BEIRA IHO').time,ds_his_sel.sel(time=slice('2019-03-11','2019-03-20')).sel(stations='BEIRA IHO').waterlevel,color='k')
lims = ax2.get_ylim()
ax2.plot([time_ts_max,time_ts_max],lims,'k--')
ax2.set_ylim(lims)
ax2.grid()
ax2.set_title('Timeseries of water levels at model output point')

fig.suptitle('Delft3D-FM model output for tropical cyclone Idai based on IBTrACS')
# %%

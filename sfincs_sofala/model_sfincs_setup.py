#%% Import the correct packages
import os
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import geopandas as gpd
import matplotlib.pyplot as plt

import hydromt
from hydromt.log import setuplog
from hydromt_sfincs import SfincsModel


#%% Sources to help create this script
# https://deltares.github.io/hydromt_sfincs/latest/_examples/build_from_script.html#
# Pizza course: p:\11208235-013-egyptian-tsunami-mode\SFINCS_course_notebooks\

#%% We select the area for the location 

#For now we select the same initial bbox as Dirk's paper. Later this should be more flexible
bbox = [34.33,-20.12,34.95,-19.30] 
model_res = 100 #By defaulft
data_catalog  = hydromt.DataCatalog(data_libs = ['datacatalog.yml']) #To correct for the location of the GTSM data

#%% Specify root_folder and logger_name
root_folder  = Path('sfincs_sofala')
logger_name = 'SFINCS_log_sofala'
logger = setuplog(logger_name, log_level=10)

#Define library path - for example if we merge deltares library with us. 

# initialize model
sf = SfincsModel(
    root=root_folder,
    mode="w+",
    data_libs = ['datacatalog.yml'],
    logger=logger,
)

# %% We use the same setup as in the .ini file from Dirk's paper
sf.setup_grid_from_region(
    region = {'bbox': bbox},  # make sure the region is in lon/lat
    res = model_res, #By default 100 m
    rotated=False, # non-rotated grid. 
    crs='utm' # automatically the closest UTM zone is selected (unit is in meters), 
)
#%% Plot model region
fig, ax = sf.plot_basemap(plot_region=True,bmap='sat')

# %% We follow the steps from the ini file from Dirk's paper (now called a yml file)
sf.setup_dep(datasets_dep= [{'elevtn': 'merit_hydro'}]) #Setup topobathy --- here bathymetry is not used!
# datasets_dep = [{"elevtn": "merit_hydro", "zmin": 0.001}, {"elevtn": "gebco"}]
_ = sf.plot_basemap(variable='dep',bmap='sat', plot_region=True) #Plotting the outcome
#%% We call osm - to be used later to define the waterlevel boundary conditions
gdf_include = sf.data_catalog.get_geodataframe('osm_coastlines', bbox=bbox)

#Plotting osm there
fig, ax = sf.plot_basemap(plot_region=True,bmap='sat')
gdf_include.to_crs(sf.crs).boundary.plot(ax=ax, color="b", lw=1, ls="--")
#fig, ax  = sf.plot_basemap(variable='msk', bmap='sat', zoomlevel=10)
#%% Set up the mask
sf.setup_mask_active(zmin=0, reset_mask=True)
sf.setup_mask_active(
    include_mask = None, # change to None if you don't 
    exclude_mask = gdf_include, 
    drop_area = 1000, 
    fill_area = 0, # not filling anything? 
    reset_mask = False
)
# Plot the mask. Using variable='msk' will display the mask values for the active cells.
_ = sf.plot_basemap(variable='msk', bmap='sat')
# %% Assigning the waterlevel downstream boundary condition
sf.setup_mask_bounds(btype = 'waterlevel',
                     zmin = 0,
                     #zmax = 1,
                     include_mask = 'osm_coastlines',
                     include_mask_buffer = 200, 
                     reset_bounds = True)
# Inspect the updated mask. Mask value 2 means waterlevel boundary, mask value 3 means outflow boundary
_ = sf.plot_basemap(variable='msk', bmap='sat')

#%% Setup river inflow and outflow
river_inflow_kwargs = dict(
    hydrography='merit_hydro',
    river_upa=500,
    river_len=1000,
    keep_rivers_geom=True,
)

sf.setup_river_inflow(**river_inflow_kwargs)


#SETUP THE RIVER OUTFLOW
sf.setup_river_outflow(
                      hydrography="merit_hydro", #rivers = 'rivers_lin2019_v1',
                      river_upa = 10, # Is don't get this
                      river_width = 4e3,
                      keep_rivers_geom= True)

sf.plot_basemap('basemap.png', bmap='sat')
#Q: what means src in the plotted map? --> discharge points

#%% We try to get the river bathymetry
hydro = sf.data_catalog.get_rasterdataset('merit_hydro', bbox=bbox)
rivers = sf.data_catalog.get_geodataframe('rivers_lin2019_v1', bbox=bbox)
#%%
#We export it to check it
source_names=["merit_hydro", "rivers_lin2019_v1"]
folder_name = "tmp_data_export"
data_catalog.export_data(
    data_root=folder_name,
    bbox=bbox,
    source_names=source_names,
    meta={"version": "1"},
)
#%% For now making dummy rivers -- why dummy as it copies rivers_inflow based on merit_hydro?
gdf_riv = sf.geoms["rivers_inflow"].copy()
gdf_riv["rivwth"] = 100 # width [m]
gdf_riv["rivdph"] = 1.5  # depth [m]
gdf_riv["manning"] = 0.03  # manning coefficient [s.m-1/3]
# gdf_riv[["geometry", "rivwth", "rivdph", "manning"]]

datasets_riv = [{'centerlines': gdf_riv}]

#%% If we use setup_subgrid, then surface roughness, more detailed elevation and rivers to be assigned with 
# setup_subgrid
# datasets_rgh = [
#     {'lulc': 'esa_worldcover', 'reclass_table': 'esa_worldcover_mapping'}
# ] #To check later - in Dirk's paper uses vito - we can change that later - sf.setup_manning_roughness(datasets_rgh = [{'manning':'vito'}],  manning_sea = 0.02)

datasets_rgh = [{"lulc": "vito"}]

#sf.setup_manning_roughness(datasets_rgh = datasets_rgh,  manning_sea = 0.02)
#%%
# Does the order determine which dataset to prioritize? 
datasets_dep = [
    {'elevtn': 'merit_hydro'}, 
    {'elevtn': 'copdem30', 'zmin' : 0.001}, # what means zmin argument here?
    {'elevtn': 'gebco'}
]

#Create the subgrid where we burn the river bathymetry as well
sf.setup_subgrid(
    datasets_dep= datasets_dep,
    datasets_rgh = datasets_rgh,
    datasets_riv = datasets_riv,
    nr_subgrid_pixels = 3,  #Fill in number of subgrid files -- can this be uneven number?
    write_dep_tif=True,  # save a cloud-optimized geotiff of the subgrid topography
    write_man_tif=True,
    nrmax=5000,  # set tile size a bit larger speed up processing (default 2000)
)

# we can plot the 2D subgrid variables
_ = sf.plot_basemap(variable="subgrid.z_zmin", plot_bounds=False, bmap="sat", zoomlevel=12)

# Use predefined plotting function 'plot_basemap' to show your full SFINCS model setup
_ = sf.plot_basemap(fn_out="basemap.png", bmap="sat", zoomlevel=12)

#%% Infiltration data
sf.setup_cn_infiltration(
    'gcn250', antecedent_moisture='dry'
)

#%% Add forcing

# Specify the simulation time in the model config - not done now
model_time_config = {
    "tref": "20190314 000000", #FILL IN THE REFERENCE TIME (can be any date)
    "tstart": "20190314 000000", #FILL IN THE START TIME OF THE SIMULATION
    "tstop": "20190323 000000", #FILL IN THE END TIME OF THE SIMULATION
    "dtout": 3600, #FILL IN THE TIMESTEP OF THE MAP OUTPUT
    "dthisout" : 3600, #FILL IN THE TIMESTEP OF THE SCALAR OUTPUT
}
sf.setup_config(**model_time_config)

#%% Set up rainfall forcing from ERA5
sf.setup_precip_forcing_from_grid(
    precip='era5_hourly',
    aggregate=False
)


#%% Set up wind forcing from ERA5
sf.setup_wind_forcing_from_grid(
    wind = 'era5_hourly'
)

#%% Set up pressure forcing from ERA5
sf.setup_pressure_forcing_from_grid(
    press='era5_hourly'
)

#%% Set up coastal water level forcing
# change to locations and timeseries
sf.setup_waterlevel_forcing(
    geodataset='gtsm_codec_reanalysis_hourly_v1',
    buffer=200
)
#%% plot all forcings
_ = sf.plot_forcing()

#%% SETUP THE DISCHARGE
#Create synthetic discharge timeseries
q1 = 100
q2 = 1000
q3 = 500
q4 = 2500
discharge = sf.forcing['dis'].copy()
discharge[:,0] = q1
discharge[:,1] = q2
discharge[:,2] = q3
discharge[:,3] = q4

sf.setup_discharge_forcing(discharge) #set discharge in model using hydroMT function
sf.forcing['dis']  # print discharge

# only write forcing
#sf.write_forcing()

#%% Plot model summary so far
_ = sf.plot_basemap(variable='dep', bmap='sat')

#%% Saving the model for now
sf.write()
#%% To reload a model already existing
# logger = setuplog('SFINCS_log_sofala', log_level=20)
# sf = SfincsModel(data_libs=['datacatalog.yml'], root=root_folder, mode='r', logger=logger)
# sf.read()

# ... do stuff

# # change the model root (to not overwrite existing model)
# model_root = 'sfincs_sofala_subgrid'
# sf.set_root(model_root, mode= 'r+')
# %%

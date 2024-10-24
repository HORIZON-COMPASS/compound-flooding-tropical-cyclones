#%% Import the correct packages
import os
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xarray as xr

import hydromt
from hydromt.log import setuplog
from hydromt_sfincs import SfincsModel

#%% Sources to help create this script
# https://deltares.github.io/hydromt_sfincs/latest/_examples/build_from_script.html#

#%% User settings

bbox = [34.33,-20.12,34.95,-19.30] # Sofala region 
model_res = 100 # resolution
datacat = os.path.join('..','..','datacatalog.yml')
modelname = 'sfincs_sofala_test'
coupling_mask = 'coastal_coupling_msk'

#For now we select the same initial bbox as Dirk's paper. Later this should be more flexible
bbox = [34.33,-20.12,34.95,-19.30] 
model_res = 100 #By defaulft
datacat = os.path.join('..','boundary_conditions','datacatalog.yml')
data_catalog  = hydromt.DataCatalog(data_libs = [datacat]) #To correct for the location of the GTSM data

#%% Specify root_folder and logger_name
root_folder  = os.path.join('..','computations','sfincs_sofala_test')
logger_name = 'SFINCS_log_sofala'
logger = setuplog(logger_name, log_level=10)

#Define library path - for example if we merge deltares library with us. 

# initialize model
sf = SfincsModel(
    root=root_folder,
    mode="w+",
    data_libs = [datacat],
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
sf.setup_dep(datasets_dep= [{'elevtn': 'merit_hydro', 'zmin': 0.001}, 
                            {'elevtn': 'gebco_v2024', 'reproj_method': 'bilinear'}]) #'offset': 'mdt_cnes_cls18',

_ = sf.plot_basemap(variable='dep',bmap='sat', plot_region=True) #Plotting the outcome
#%% We call osm - to be used later to define the waterlevel boundary conditions
gdf_include = sf.data_catalog.get_geodataframe(coupling_mask, bbox=bbox) # 'osm_coastlines' can also be used

#Plotting osm there
fig, ax = sf.plot_basemap(plot_region=True,bmap='sat')
gdf_include.to_crs(sf.crs).boundary.plot(ax=ax, color="b", lw=1, ls="--")
#fig, ax  = sf.plot_basemap(variable='msk', bmap='sat', zoomlevel=10)
#%% Set up the mask
sf.setup_mask_active(zmin=-10, reset_mask=True)
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
                     zmin = -10,
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

#%%
# To check the river network interactively: 
# Here we can see which segments are present in the river network
sf.geoms['rivers_inflow'].explore()


#%% We try to get the river bathymetry
hydro = sf.data_catalog.get_rasterdataset('merit_hydro', bbox=bbox)
rivers = sf.data_catalog.get_geodataframe('rivers_lin2019_v1', bbox=bbox)
#%%
#We export it to check it
source_names=["merit_hydro", "rivers_lin2019_v1"]
folder_name = os.path.join('..','boundary_conditions',"tmp_data_export")
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

datasets_rgh = [{"lulc": "vito", 'reclass_table': 'vito_mapping'}]

#sf.setup_manning_roughness(datasets_rgh = datasets_rgh,  manning_sea = 0.02)
#%%
# Does the order determine which dataset to prioritize? 
datasets_dep = [{'elevtn': 'merit_hydro', 'zmin': 0.001}, 
                 {'elevtn': 'gebco_v2024', 'offset': 'mdt_cnes_cls18', 'reproj_method': 'bilinear'}
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
    'gcn250', antecedent_moisture='dry' # can be changed to 'avg' and 'wet'
)

#%% Add forcing

# Specify the simulation time in the model config - 14 - 23 march - but 9-16 for now from IbTracks)
model_time_config = {
    "tref": "20190309 000000", #FILL IN THE REFERENCE TIME (can be any date)
    "tstart": "20190309 000000", #FILL IN THE START TIME OF THE SIMULATION
    "tstop": "20190316 000000", #FILL IN THE END TIME OF THE SIMULATION
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
    geodataset='dfm_output_MZ_doris', # 'gtsm_codec_reanalysis_hourly_v1'
    buffer=5000
)

# #%% SETUP THE DISCHARGE - Synthetic discharge for now
# Generate new time range for discharge that matches the model settings
tstart = pd.to_datetime(model_time_config["tstart"], format="%Y%m%d %H%M%S")
tstop = pd.to_datetime(model_time_config["tstop"], format="%Y%m%d %H%M%S")
new_time = pd.date_range(start=tstart, end=tstop, freq=f'{model_time_config["dtout"]}S')

# Extend discharge data to the time length of the model
dis = sf.forcing['dis']
extended_data = np.zeros((len(new_time), dis.sizes['index'])) # extend data length
extended_data[:dis.sizes['time'], :] = dis.values # add correct time values

# Create new DataArray with extended time
sf.forcing['dis'] = xr.DataArray(
    extended_data,
    coords={
        'time': new_time,
        'index': dis.coords['index'],
        'geometry': dis.coords['geometry'],
        'uparea': dis.coords['uparea']
    },
    dims=['time', 'index']
)

# Set custom discharge values
sf.forcing['dis'][:, 0] = 100    # For index 1
sf.forcing['dis'][:, 1] = 10000  # For index 2
sf.forcing['dis'][:, 2] = 5000   # For index 3
sf.forcing['dis'][:, 3] = 2500   # For index 4

# Print to verify
print(sf.forcing['dis'])

#%% plot all forcings
sf.write_forcing()
_ = sf.plot_forcing()

#%%
# Add observation points from Eilander et al. (2022)
# https://zenodo.org/records/7274465 - 2_code/2_experiment/obs_locs.geojson
obs_points = Path('p:/11210471-001-compass/01_Data/obs_locs.geojson')
sf.setup_observation_points(locations=obs_points, merge=False)

#%% Plot model summary so far
_ = sf.plot_basemap(variable='dep', bmap='sat')

#%% Saving the model for now
sf.write()

#%%
def print_directory_tree(directory):
    for root, dirs, files in os.walk(directory):
        level = root.replace(directory, '').count(os.sep)
        indent = ' ' * 2 * (level)
        print('{}{}/'.format(indent, os.path.basename(root)))
        subindent = ' ' * 2 * (level + 1)
        for f in files:
            print(f'{subindent}+ {f}')

print_directory_tree(sf.root)

#%%
# make a sfincs_log.txt file to write output when running SFINCS
# Create an empty text file
file_name = "sfincs_log.txt"
with open(os.path.join(root_folder, file_name), "w") as file:
    pass

#%%
#  create a bat file to run the model on windows
batch_content = 'call "p:/11210471-001-compass/02_Models/00_executables/SFINCS_v2.1.1_Dollerup_release_exe/sfincs.exe" > sfincs_log.txt'

file_name = 'run_sfincs.bat'
with open(os.path.join(root_folder, file_name), "w") as file:
    file.write(batch_content)
# Now you can run the SFINCS model by opening the bat file on a windows computer

#%% To reload a model already existing
# logger = setuplog('SFINCS_log_sofala', log_level=20)
# sf = SfincsModel(data_libs=['datacatalog.yml'], root=root_folder, mode='r', logger=logger)
# sf.read()

# ... do stuff

# # change the model root (to not overwrite existing model)
# model_root = 'sfincs_sofala_subgrid'
# sf.set_root(model_root, mode= 'r+')
# %%

#%% Checking outcome of wflow

#Import packages
#%%
import pandas as pd
import matplotlib.pyplot as plt
#from matplotlib import colors
import os
from hydromt_wflow import WflowModel
import xarray as xr
#import geopandas as gpd


#%%
#We check where the stations are:
#The folder is the location of the wflow model in 02_Models
model_path = r"p:\11210471-001-compass\02_Models\sofala\Idai\wflow"
mod_ini = WflowModel(root=os.path.join(model_path), mode="r+", config_fn=os.path.join(model_path, "wflow_sbm.toml"))
q_locs = mod_ini.geoms["gauges_locs"]

#Plotting
fig, ax = plt.subplots()
mod_ini.geoms["rivers"].plot(ax=ax)
q_locs.plot(ax=ax, color="red")
for idx, row in q_locs.iterrows():
    plt.text(row.geometry.x, row.geometry.y, q_locs.loc[idx,'index'], fontsize=9, ha='right')
plt.show()

#%%
#We check where the stations are:
from os.path import join
from hydromt.log import setuplog
curdir           = '../../../'
data_cats        = [
        join(curdir, "03_data_catalogs", "datacatalog_general.yml"), 
        join(curdir, "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling.yml"), 
        join(curdir, "03_data_catalogs", "datacatalog_SFINCS_obspoints.yml"),
        join(curdir, "03_data_catalogs", "datacatalog_CF_forcing.yml")
        ]
logger = setuplog("update", "./hydromt.log", log_level=10)

#The folder is the location of the wflow model in 02_Models
wflow_root_30yr = r"p:\11210471-001-compass\03_Runs\sofala\Idai\wflow\event_precip_era5_hourly_zarr_CF0_30yr"
mod = WflowModel(
        root=join(wflow_root_30yr, "warmup", "run_default"),
        data_libs=data_cats,
        mode="r+",
        logger=logger,
        config_fn=os.path.join(wflow_root_30yr, "warmup", "run_default", "wflow_sbm.toml")
    )
q_locs = mod.geoms["gauges_locs"]

#Plotting
fig, ax = plt.subplots()
mod.geoms["rivers"].plot(ax=ax)
q_locs.plot(ax=ax, color="red")
for idx, row in q_locs.iterrows():
    plt.text(row.geometry.x, row.geometry.y, q_locs.loc[idx,'index'], fontsize=9, ha='right')
plt.show()
#%%We import the modelled data - from 03_Runs
Folder_start = "event_precip_era5_hourly_zarr_CF0" 
# Folder_events = os.path.join(Folder_start, "events")

fn_config = r"c:\CODE\COMPASS\wflow_checks\event_precip_era5_hourly_zarr_CF0\events\wflow_sbm.toml"
#fn_config = os.path.join(Folder_start, "events", "wflow_sbm.toml")
mod_event = WflowModel(root=os.path.join(Folder_start, "events"), mode="r+", config_fn=fn_config)
mod_event.results.keys()

Q = mod_event.results['netcdf']

# # Alternative is that we load the results directly from the netcdf file
# ds = xr.open_dataset(os.path.join( Folder_p, "run_default", "output_scalar.nc"))
# Q = ds["Q"]


#%%
#We plot the results
plt.figure()
ds["Q"].sel(Q_gauges_locs = "1").plot(label="Gauge 1")
ds["Q"].sel(Q_gauges_locs = "2").plot(label="Gauge 2")
plt.legend()

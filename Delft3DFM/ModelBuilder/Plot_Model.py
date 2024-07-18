# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:12:56 2023

@author: Tammo Zijlker
"""
#%% 
import os
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
from dfm_tools import modelbuilder as mb
import hydrolib.core.dflowfm as hcdfm
import xarray as xr
import pandas as pd
import contextily as ctx
import getpass
from pathlib import Path
import xugrid as xu

#%% Settings

#Bounding box for model cutout
# lon_min, lon_max, lat_min, lat_max = -15.0, 13.0, 43.116760, 64
# lon_min, lon_max, lat_min, lat_max = 52.70, 58.38,-24.43, -19.20 # Reunion
lon_min, lon_max, lat_min, lat_max = 52.70, 58.38,-24.43, -19.20 # Tacloban

## input
# dir_output = Path(r'p:\11208614-de-370a\01_models\Basque\DFLOWFM\Models\DCSM_GTSM_cutout_2013-2017')
# dir_output = Path(r'p:\11208614-de-370a\01_models\Reunion\DFLOWFM\Models\reunion_custom_grid')
# dir_output = Path(r'p:\11208614-de-370a\01_models\Philippines\DFLOWFM\Tacloban_Bay')
# dir_output = Path(r'p:\11208614-de-370a\01_models\Philippines\DFLOWFM\Model_Version0')
dir_output = Path(r'p:\11208614-de-370a\01_models\BES\DFLOWFM\Models\st_maarten')

path_style = 'unix' # windows / unix
overwrite = False # used for downloading of forcing data. Always set to True when changing the domain
paths_relative = False #TODO: currently only works with path_style='windows' (same OS as IDE)

#Model files (existing)
# poly_file = os.path.join(dir_output, 'humber_seaboundary.pli')
# poly_file = os.path.join(dir_output, 'tacloban_bay.pli')
# network_file = 'humber_net.nc'
# network_file = 'DFM_interpreted_idomain_humber_net.nc'
# network_file = 'reunion_net.nc'
# network_file = 'tacloban_bay_net.nc'
network_file = 'st_maarten_net.nc'
# network_file = 'phil_gtsm_cutout_v1_net.nc'

# local_grid = dfmt.open_partitioned_dataset(str(dir_output / network_file))
local_grid = xu.open_dataset(str(dir_output / network_file))
# seaboundary = hcdfm.PolyFile(poly_file)
# bnd  = dfmt.pointlike_to_DataFrame(seaboundary.objects[0])  # convert to dataframe

#observation points sfincs_boundary
# obs = hcdfm.XYNModel(dir_output / 'sfincs_obs.xyn') #load obsfile
# obs = hcdfm.XYNModel(dir_output / 'sfincs_basque_moved_obs.xyn') #load obsfile
# obs = hcdfm.XYNModel(dir_output / 'reunion_SFINCS_obs.xyn') #load obsfile
# obs = dfmt.pointlike_to_DataFrame(obs)  # convert to dataframe
lon_min, lat_min, lon_max, lat_max = local_grid.grid.bounds
#%%
fig, ax = plt.subplots()
local_grid.grid.plot(ax=ax,linewidth=1)
# ax.scatter(bnd['x'], bnd['y'], color = 'black', marker = 'o', label=None)
ax.set_xlim([lon_min-0.1, lon_max+0.1])
ax.set_ylim([lat_min-0.1, lat_max+0.1])

ax.set_xlabel(r'Longitude [$^\circ$]')
ax.set_ylabel(r'Latitude [$^\circ$]')

# Hide X and Y axes label marks
# ax.xaxis.set_tick_params(labelbottom=False)
# ax.yaxis.set_tick_params(labelleft=False)

# Hide X and Y axes tick marks
# ax.set_xticks([])
# ax.set_yticks([])
f = 1.0/np.cos(np.mean([lat_min, lat_max])*np.pi/180)

plt.gca().set_aspect(f)

ctx.add_basemap(ax=ax, crs='EPSG:4326', attribution=False,source=ctx.providers.Esri.WorldImagery)
# plt.savefig(dir_output / 'figures' / 'network' / 'basque_cutout_grid.png', dpi=300, bbox_inches='tight')	
# plt.savefig(dir_output / 'figures' / 'network' / 'reunion_cutout_grid.png', dpi=300, bbox_inches='tight')	
plt.savefig(dir_output / 'figures' / 'network' / 'grid.png', dpi=300, bbox_inches='tight')	
#%% Make zoom over SFINCS domain
fig, ax = plt.subplots()
# lon_min, lon_max, lat_min, lat_max = -0.8, 0.5, 53.2, 54 # Humber
# lon_min, lat_min, lon_max, lat_max = -4.328613,431,-0.986023,44.549378 # Basque
# lon_min, lon_max, lat_min, lat_max =  54.7,56.4,-21.9,-20.5 # Reunion
# lon_min, lat_min, lon_max, lat_max = 124.832153,10.606620,125.565491,11.383109 # Tacloban
lon_min, lat_min, lon_max, lat_max = -63.197479,17.98,-62.903595,18.145199 # ST. Maarten
local_grid.ugrid.sel(x=slice(lon_min, lon_max), y=slice(lat_min, lat_max))\
    .grid.plot(ax=ax,linewidth=1, zorder = 1)

# ax.scatter(obs['x'], obs['y'], color = 'red', marker = '*', label='SFINCS boundary', zorder = 2)

ax.set_xlabel(r'Longitude [$^\circ$]')
ax.set_ylabel(r'Latitude [$^\circ$]')

ax.set_xlim([lon_min, lon_max])
ax.set_ylim([lat_min, lat_max])

ax.set_xlabel(r'Longitude [$^\circ$]')
ax.set_ylabel(r'Latitude [$^\circ$]')
# ax.set_xlabel('')
# ax.set_ylabel('')

# Hide X and Y axes label marks
# ax.xaxis.set_tick_params(labelbottom=False)
# ax.yaxis.set_tick_params(labelleft=False)

# Hide X and Y axes tick marks
# ax.set_xticks([])
# ax.set_yticks([])
f = 1.0/np.cos(np.mean([lat_min, lat_max])*np.pi/180)

plt.gca().set_aspect(f)
# ax.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0), ncol=1)

ctx.add_basemap(ax=ax, crs='EPSG:4326', attribution=False,source=ctx.providers.Esri.WorldImagery, zorder = 0)
# plt.savefig(dir_output / 'figures' / 'network' / 'humber_cutout_grid_zoom_humber.png', dpi=300, bbox_inches='tight')	
plt.savefig(dir_output / 'figures' / 'network' / 'grid_zoom.png', dpi=300, bbox_inches='tight')

#%% Plot bathymetry map
fig, ax = plt.subplots()

local_grid['mesh2d_node_z']\
    .ugrid.sel(x=slice(lon_min, lon_max), y=slice(lat_min, lat_max))\
    .ugrid.plot(ax=ax, cmap = 'jet', vmax = 0 )
local_grid.ugrid.sel(x=slice(lon_min, lon_max), y=slice(lat_min, lat_max))\
    .grid.plot(ax=ax,linewidth=0.2, zorder = 1)
    
ax.set_xlabel(r'Longitude [$^\circ$]')
ax.set_ylabel(r'Latitude [$^\circ$]')

f = 1.0/np.cos(np.mean([lat_min, lat_max])*np.pi/180)

plt.gca().set_aspect(f)
ctx.add_basemap(ax=ax, crs='EPSG:4326', attribution=False,source=ctx.providers.Esri.WorldImagery, zorder = 0)
plt.savefig(dir_output / 'figures' / 'network' / 'zoom_grid_bathymetry.png', dpi=300, bbox_inches='tight')

# %%

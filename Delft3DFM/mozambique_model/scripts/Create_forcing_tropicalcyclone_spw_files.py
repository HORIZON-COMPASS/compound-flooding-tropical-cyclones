# This script makes use of the CoastalHazardsToolkit (CHT): https://github.com/Deltares-research/CoastalHazardsToolkit/ 
# To use the underlying scripts, the CHT needs to be installed locally (pulled from GitHub and pip install -e . in the folder of the repository)

# Script by: n-aleksandrova


# Bugs and to be fixed:
# script says: "self.background_pressure = 1012 Pa" I think the units should not be Pa, but mbar, atm pressure is in the order of 10^5 Pa.

#%%
# Load packages
from cht_cyclones.tropical_cyclone import TropicalCyclone
import matplotlib.pyplot as plt
import os
import xarray as xr
import pandas as pd
import numpy as np

#%%
# specify cyclone name
tc_name = 'Freddy'#'Kenneth' #'Idai'
tc_year = 2023

dir_base = r'p:\11210471-001-compass\02_Models\Delft3DFM\mozambique_model'

#%%
# Define directory
#dir = os.path.dirname(r'c:\Users\aleksand\git_projects\CoastalHazardsToolkit\src\cht\tropical_cyclone')

# %%
# directory where the IBTRACS database is stored
print('Open IBTrACS dataset...')
ds_ibtracs = xr.open_dataset(r'p:\11210471-001-compass\01_Data\IBTrACS\IBTrACS.ALL.v04r01.nc')

# %%
def create_track(ds_tc):
    # Initialize the tropical cyclone object
    print('- Initializing tc object...')
    tc = TropicalCyclone(name=tc_name)
    tc.nr_radial_bins = 600
    tc.phi_spiral = 22.6
    tc.spiderweb_radius = 900
    tc.extend_track = 3

    # Only keep the data that is not NaN (filtered based on rmw availability)
    tmp = ds_tc.time.where(~ds_tc.usa_rmw.isnull(),drop=True).values
    data_time = [pd.to_datetime(i) for i in tmp]
    data_lon = ds_tc.usa_lon.where(~ds_tc.usa_rmw.isnull(),drop=True).fillna(-999)
    data_lat = ds_tc.usa_lat.where(~ds_tc.usa_rmw.isnull(),drop=True).fillna(-999)
    data_wind = ds_tc.usa_wind.where(~ds_tc.usa_rmw.isnull(),drop=True).fillna(-999)
    data_pres = ds_tc.usa_pres.where(~ds_tc.usa_rmw.isnull(),drop=True).fillna(-999)
    data_rmw = ds_tc.usa_rmw.where(~ds_tc.usa_rmw.isnull(),drop=True).fillna(-999)
    data_r34 = ds_tc.usa_r34.where(~ds_tc.usa_rmw.isnull(),drop=True).fillna(-999)

    # fill in track data
    print('- Filling in track data...')
    tc.provide_track(datetimes = data_time, lons = data_lon.values, lats = data_lat.values,
                 winds = data_wind.values, pressures = data_pres.values,
                 rmw = data_rmw.values, r35 = data_r34.values)
    
    # functions to improve the data
    print('- Improve schematization...')
    tc.account_for_forward_speed()
    tc.estimate_missing_values()

    return tc

#%%
# Find the cyclone in the IBTRACS database
# Note: the record for one cyclone can contain several parts, this is the case for Freddy (2023)
print(f'Find TC {tc_name} ({tc_year}) in the database...')

ds_ibtracs = ds_ibtracs.where(ds_ibtracs.season == tc_year,drop=True)
id = []
for ii,ids in enumerate(ds_ibtracs.storm):
    yr = pd.to_datetime(ds_ibtracs.isel(date_time=0).time.values[ii]).year
    if (ds_ibtracs.name.values[ii].decode(encoding="utf-8") == tc_name.upper()) & (yr == tc_year):
        id.append(ids.item())

# %%

for id_track in id:
    # select only the data for this cyclone
    ds_tc = ds_ibtracs.isel(storm=id_track,drop=True)

    print(f'Processing TC {tc_name} ({tc_year}),  track {ds_tc.sid} ...') 
    
    # crop the length of the dataset for specific cyclones
    # e.g. Freddy (2023) record is very long, we do not need the part that is far beyond the model domain
    if (tc_name == 'Freddy') & (tc_year == 2023):
        ds_tc = ds_tc.where(ds_tc.time > np.datetime64('2023-02-21'),drop=True)

    # create spw file for this specific track
    tc = create_track(ds_tc)

    # export to spiderweb
    print('- Saving track...')
    tc.to_spiderweb(os.path.join(dir_base,'boundary_conditions','meteo','TC',f'tc_{tc_name.upper()}_{ds_tc.isel(date_time=0).sid.item().decode(encoding="utf-8")}.spw'))

    del tc

#ds_ibtracs.close(); del ds_ibtracs
# %%

# This script makes use of the CoastalHazardsToolkit (CHT): https://github.com/Deltares-research/CoastalHazardsToolkit/ (as of 2024 only available for the Deltares research organization on GitHub)
# To use the underlying scripts, the CHT needs to be installed locally (pulled from GitHub and pip install -e . in the folder of the repository)

# Script by: n-aleksandrova


# Bugs and to be fixed:
# - cht_cyclones expects a package called 'fiona', is not installed by default
# np.NaN needs to be replaced by np.nan for compatibility with numpy 2.0 (in tropical_cyclone.py) - for now downgraded to numpy 1.26.4
# script says: "self.background_pressure = 1012 Pa" I think the units should not be Pa, but mbar, atm pressure is in the order of 10^5 Pa.

#%%
# Load packages
from cht_cyclones.tropical_cyclone import TropicalCyclone
import matplotlib.pyplot as plt
import os
import xarray as xr
import pandas as pd

#%%
# specify cyclone name
tc_name = 'Idai'

#%%
# Define directory
#dir = os.path.dirname(r'c:\Users\aleksand\git_projects\CoastalHazardsToolkit\src\cht\tropical_cyclone')

# %%
# directory where the IBTRACS database is stored
ds_ibtracs = xr.open_dataset(r'p:\archivedprojects\11205281-tcwise\00_data\IBTrACS.ALL.v04r00_download240223.nc')

# %%
# Initialize the tropical cyclone object
tc = TropicalCyclone(name=tc_name)

# %%
#Change default settings
#tc.background_pressure = 1020
tc.nr_radial_bins = 600
tc.phi_spiral = 22.6
tc.spiderweb_radius = 900

# merge factor

#%%
# Find the cyclone in the IBTRACS database
for ii,ids in enumerate(ds_ibtracs.storm):
    if ds_ibtracs.name.values[ii].decode(encoding="utf-8") == tc_name.upper():
        id = ids.item()
        break

# %%
# select only the data for this cyclone
ds_tc = ds_ibtracs.isel(storm=id,drop=True)
ds_ibtracs.close(); del ds_ibtracs

# %%
# Only keep the data that is not NaN (filtered based on rmw availability)
tmp = ds_tc.time.where(~ds_tc.usa_rmw.isnull(),drop=True).values
data_time = [pd.to_datetime(i) for i in tmp]
data_lon = ds_tc.usa_lon.where(~ds_tc.usa_rmw.isnull(),drop=True).fillna(-999)
data_lat = ds_tc.usa_lat.where(~ds_tc.usa_rmw.isnull(),drop=True).fillna(-999)
data_wind = ds_tc.usa_wind.where(~ds_tc.usa_rmw.isnull(),drop=True).fillna(-999)
data_pres = ds_tc.usa_pres.where(~ds_tc.usa_rmw.isnull(),drop=True).fillna(-999)
data_rmw = ds_tc.usa_rmw.where(~ds_tc.usa_rmw.isnull(),drop=True).fillna(-999)
data_r34 = ds_tc.usa_r34.where(~ds_tc.usa_rmw.isnull(),drop=True).fillna(-999)
# %%
# fill in track data
tc.provide_track(datetimes = data_time, lons = data_lon.values, lats = data_lat.values,
                 winds = data_wind.values, pressures = data_pres.values,
                 rmw = data_rmw.values, r35 = data_r34.values)
# %%
# functions to improve the data
tc.account_for_forward_speed()
tc.estimate_missing_values()
# %%
# export to spiderweb
tc.to_spiderweb(os.path.join('..','boundary_conditions','meteo','TC',f'tc_{tc_name.upper()}.spw'))
# %%
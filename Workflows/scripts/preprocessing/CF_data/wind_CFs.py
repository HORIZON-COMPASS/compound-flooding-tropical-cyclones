# This script is adapted from the Create_forcing_tropicalcyclone_spw_files.py
#%%
# Importing the necessary packages
import os
import numpy as np
import xarray as xr
import pandas as pd
# import hydromt
import geopandas as gpd
import copy
import yaml
from datetime import datetime
from cht_cyclones.tropical_cyclone import TropicalCyclone

# Read snakemake values if run in workflow, 
# otherwise use absolute values to test script
if "snakemake" in locals():
    start_date = np.datetime64(snakemake.params.start_date) 
    end_date = np.datetime64(snakemake.params.end_date) 
    # bbox = ast.literal_eval(snakemake.params.bbox)
    # path_data_cat = snakemake.params.data_cat
    # CF_catalog_path = snakemake.params.CF_data_cat
    tc_name = snakemake.wildcards.tc_name
    CF_value = float(snakemake.wildcards.CF_value_wind)
    CF_value_txt = snakemake.wildcards.CF_value_wind
    output_CF_wind = os.path.abspath(snakemake.output.CF_wind)
    ibtracs_path = os.path.abspath(snakemake.input.ibtracs_path)
else:
    start_date = np.datetime64("2019-03-09") 
    end_date = np.datetime64("2019-03-24") 
    # bbox = [34.33,-20.12,34.95,-19.30]
    # path_data_cat = os.path.abspath("../../../data_catalogs/datacatalog_general.yml")
    tc_name = "Idai"
    CF_value = -10
    CF_value_txt = "-10"
    output_CF_wind = os.path.abspath(f"p:/11210471-001-compass/01_Data/counterfactuals/wind/{tc_name}_{CF_value}.nc")
    # CF_catalog_path = os.path.abspath("../../../data_catalogs/datacatalog_CF_forcing.yml")
    ibtracs_path = 'p:/11210471-001-compass/01_Data/IBTrACS/IBTrACS.ALL.v04r01.nc'
#%%
# specify cyclone name
tc_name = 'Idai' #'Kenneth', 'Idai', 'Freddy'
tc_year = start_date.astype('datetime64[ms]').astype(datetime).year
extend_days = 9 # change to check with how many days it need to be extended

#%%
# directory where the IBTRACS database is stored
print('Open IBTrACS dataset...')
ds_ibtracs = xr.open_dataset(ibtracs_path)

#%%
def create_track(ds_tc):
    # Initialize the tropical cyclone object
    print('- Initializing tc object...')
    tc = TropicalCyclone(name=tc_name)
    tc.nr_radial_bins = 600
    tc.phi_spiral = 22.6
    tc.spiderweb_radius = 900
    tc.extend_track = extend_days+1

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

    # drop unnecessary variables
    sid = ds_tc.sid.item().decode(encoding="utf-8")
    keys = list(ds_tc.keys()) # remove all variables except water level
    keyslist = ['usa_lon','usa_lat','usa_wind','usa_pres','usa_rmw','usa_r34'] # these are the variables to keep
    for kk in keyslist:
        keys.remove(kk)
    ds_tc = ds_tc.drop_vars(keys)
    print(ds_tc)
    
    # crop the length of the dataset for specific cyclones
    # e.g. Freddy (2023) record is very long, we do not need the part that is far beyond the model domain
    if (tc_name == 'Freddy') & (tc_year == 2023):
        ds_tc = ds_tc.where(ds_tc.time > np.datetime64('2023-02-21'),drop=True)

    # TO DO: define extend_days based on datetime


    # create spw file for this specific track
    # tc = create_track(ds_tc)

    # export to spiderweb
    # print('- Saving track...')
    # tc.to_spiderweb(os.path.join(dir_base,'boundary_conditions','meteo','TC',f'tc_{tc_name.upper()}_{sid}.spw'))

    # del tc

#ds_ibtracs.close(); del ds_ibtracs
#%%
ds_tc["usa_r34"]


#%%
import matplotlib.pyplot as plt
ds = ds_tc
# Assuming your dataset is loaded as 'ds'
lat = ds['lat'].values
lon = ds['lon'].values
wind = ds['usa_pres'].values

# Remove NaN values to avoid plotting issues
valid_data = ~np.isnan(lat) & ~np.isnan(lon) & ~np.isnan(wind)
lat = lat[valid_data]
lon = lon[valid_data]
wind = wind[valid_data]

# Create the scatter plot
plt.figure(figsize=(10, 6))
sc = plt.scatter(lon, lat, c=wind, cmap="viridis", s=20, edgecolor="k")
plt.colorbar(sc, label="Wind Speed (m/s)")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.title("Spatial Distribution of Tropical Cyclone Wind")
plt.grid()
plt.show()



# %%
# create spw file for this specific track
tc = create_track(ds_tc)

# export to spiderweb
print('- Saving track...')
tc.to_spiderweb(os.path.join(dir_base,'boundary_conditions','meteo','TC',f'tc_{tc_name.upper()}_{sid}.spw'))

del tc


# %%

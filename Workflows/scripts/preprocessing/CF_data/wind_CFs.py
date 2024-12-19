# This script is adapted from the Create_forcing_tropicalcyclone_spw_files.py
#%%
# Importing the necessary packages
import os
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from datetime import datetime
from cht_cyclones.tropical_cyclone import TropicalCyclone

# Read snakemake values if run in workflow, 
# otherwise use absolute values to test script
if "snakemake" in locals():
    start_date = np.datetime64(snakemake.params.start_date) 
    end_date = np.datetime64(snakemake.params.end_date) 
    tc_name = snakemake.wildcards.tc_name
    CF_value = float(snakemake.wildcards.CF_value_wind)
    CF_value_txt = snakemake.wildcards.CF_value_wind
    output_CF_wind = os.path.abspath(snakemake.output.CF_wind)
else:
    start_date = np.datetime64("2019-03-09") 
    end_date = np.datetime64("2019-03-24") 
    tc_name = "Idai"
    CF_value = -10
    CF_value_txt = "-10"
    output_CF_wind = f"p:/11210471-001-compass/01_Data/counterfactuals/wind/tc_{tc_name}_{CF_value_txt}.spw"
    
#%%
#extract TC year
tc_year = start_date.astype('datetime64[ms]').astype(datetime).year

#%%
# directory where the IBTRACS database is stored
ibtracs_path = 'p:/11210471-001-compass/01_Data/IBTrACS/IBTrACS.ALL.v04r01.nc'
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
    # if (tc_name == 'Freddy') & (tc_year == 2023):
    #     ds_tc = ds_tc.where(ds_tc.time > np.datetime64('2023-02-21'),drop=True)

ds_ibtracs.close(); del ds_ibtracs

#%%
# Calculated extended days necessary of the track
# Filter out NaT values
valid_times = ds_tc.time.values[~np.isnat(ds_tc.time.values)]

# Get the last valid time
if len(valid_times) > 0:
    last_valid_time = valid_times[-1]

    # Calculate the difference with end_date defined by snakemake
    difference = (end_date - last_valid_time).astype('timedelta64[D]').item()

    # Define extend_days if the difference is greater than zero
    if difference.days > 0:
        extend_days = difference.days
        print(f"extend_days: {extend_days}")
    else:
        print("No extension needed.")
else:
    print("No valid time values found in the dataset.")

#%%
# Load the era5 dataset for environmental pressure estimate 
# era5_ds = xr.open_dataset("p:/wflow_global/hydromt/meteo/era5_daily/nc_merged/era5_2019_daily.nc")
# era5_ds = era5_ds.rename({"latitude": "lat", "longitude": "lon"})
# #%%
# # Extract data from ERA5 datasets to obtain an estimate for environmental pressure
# tc_lat, tc_lon, usa_pres, tc_time = ds_tc['lat'].values, ds_tc['lon'].values, ds_tc['usa_pres'].values, ds_tc['time'].values
# era5_lat, era5_lon, era5_time, era5_pressure = era5_ds['lat'].values, era5_ds['lon'].values, era5_ds['time'].values, era5_ds['msl'].values # convert Pa to hPa

# # Remove NaN values
# valid_data = ~np.isnan(tc_lat) & ~np.isnan(tc_lon) & ~np.isnan(usa_pres)
# tc_lat, tc_lon, usa_pres, tc_time = tc_lat[valid_data], tc_lon[valid_data], usa_pres[valid_data], tc_time[valid_data]

# # Match TC time with ERA5 time (find closest time step in ERA5)
# era5_time_index = np.abs(era5_time - tc_time[:, None]).argmin(axis=1)

# # Find the closest grid points in ERA5 and extract msl values
# matched_lat, matched_lon, matched_time, matched_era5_pressure = [], [], [], []
# for i, time_index in enumerate(era5_time_index):
#     lat_idx, lon_idx = np.abs(era5_lat - tc_lat[i]).argmin(), np.abs(era5_lon - tc_lon[i]).argmin()
#     matched_lat.append(era5_lat[lat_idx])
#     matched_lon.append(era5_lon[lon_idx])
#     matched_time.append(era5_time[time_index])
#     matched_era5_pressure.append(era5_pressure[time_index, lat_idx, lon_idx]/100) # divide by 100 to convert from Pa to hPa = mb

#     print(matched_era5_pressure)
    

# matched_era5_ds = xr.Dataset(
#     {
#         "era5_msl": ("date_time", matched_era5_pressure),
#     },
#     coords={
#         "time": ("date_time", tc_time),
#         "lat": ("date_time", matched_lat),
#         "lon": ("date_time", matched_lon),
#     }
# )

# Plot the matched ERA5 pressure data as a scatter plot over the TC points
# fig, ax = plt.subplots(figsize=(10, 6))
# scatter = ax.scatter(matched_ds['lon'], matched_ds['lat'], c=matched_ds['matched_pressure'], cmap='coolwarm', edgecolors='k', s=50)
# fig.colorbar(scatter, label="Matched ERA5 Pressure (hPa)")
# ax.set_title("Matched ERA5 Pressure at TC Points")
# ax.set_xlabel("Longitude")
# ax.set_ylabel("Latitude")
# ax.grid(True)
# plt.show()

#%%
#### CF calculations ####
# Create counterfactual wind based on CF_value
ds_tc["usa_wind"] = ds_tc["usa_wind"] * ((100 + CF_value)/100)


# Correct for the cooresponsing (small) change in pressure:
# The central pressure at each track position is increased by CF_value times 
# the difference between central pressure and environmental/background pressure,
# defined in cht-cyclones as self.background_pressure = 1012 Pa
ds_tc["usa_pres"] = ds_tc["usa_pres"] + ((100 - CF_value)/100) * (1012 - ds_tc["usa_pres"])

# # Mask for non-NaN values in the 'usa_pres' column
# valid_usa_pres_mask = ~np.isnan(ds_tc["usa_pres"])

# # Apply the operation only for valid 'usa_pres' values
# ds_tc["usa_pres"][valid_usa_pres_mask] = ds_tc["usa_pres"][valid_usa_pres_mask] + ((100 - CF_value)/100) * (matched_era5_ds["era5_msl"] - ds_tc["usa_pres"][valid_usa_pres_mask])

#%%
# create spw file for this specific track
tc = create_track(ds_tc)

#%%
# export to spiderweb
print('- Saving track...')
tc.to_spiderweb(output_CF_wind)

del tc


# %%

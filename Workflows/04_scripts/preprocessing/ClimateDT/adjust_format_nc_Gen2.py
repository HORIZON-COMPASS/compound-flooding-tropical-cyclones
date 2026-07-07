
#%%
import matplotlib.pyplot as plt
import pandas as pd
import os
from datetime import datetime
import numpy as np
import xarray as xr
from glob import glob
import rioxarray as rxr

#%%
main_dir = "p:/11210471-001-compass/01_Data/ECMWF_ClimateDT"
# data_dir = os.path.join(main_dir, 'data', 'nc_Sofala')
data_dir = os.path.join(main_dir, 'data', 'nc_Idai_full')

out_dir = os.path.join(main_dir, 'data','preprocessed','Gen2','Idai_full')
os.makedirs(out_dir, exist_ok=True)

#%%
varnames = ['avg_tprate', 'msl','10u','10v', '2t', 'avg_sdswrf', 'tisr']
experiments = ["cont","hist","Tplus2.0K"]
areaname = 'Idai_full'
bbox = [32, -23, 36, -17]
bbox = [30, -22, 45, -12]     # Idai
# incl madagascar piece [30, -22, 45, -14.5] 
renamevars = {'avg_tprate':'tp', 'msl':'msl','10u':'10u','10v':'10v', '2t':'2t', 'avg_sdswrf':'ssrd', 'avg_tdswrf':'tisr', 'avg_snswrf': 'ssr'}


standard_names = {'10u':'eastward_wind', '10v':'northward_wind', 'msl':'air_pressure', 
                  '2t':'air_temperature', 'ssrd':'surface_solar_radiation_downwards',
                  'tisr':'top_incoming_solar_radiation',
                  'tp':'precipitation', 'ssr':'surface_net_shortwave_radiation'}
units = {'10u':'m/s', '10v':'m/s', 'msl':'Pa', '2t':'K', 'ssrd':'W/m2', 'tisr':'W/m2', 'tp':'mm/hr'}

realizations = ['1', '2', '3', '4', '5']

#%%
for experiment in experiments[1:2]:
    for varname in varnames[1:4]:
        for realization in realizations:
    
            files = glob(os.path.join(data_dir, f'climateDT_{varname}_{experiment}_{areaname}_*_r{realization}.nc'))
            print(files)
            ds = xr.open_mfdataset(files)
            ds.load()

            ds = ds.sortby('latitude', ascending=True) # sort latitude in ascending order
            ds = ds.sel(latitude=slice(bbox[1], bbox[3]), longitude=slice(bbox[0], bbox[2]))

            if varname == 'avg_tprate':
                ds[varname] = ds[varname] * 3600 # convert from kg/m2/s to mm/hr

            ds = ds.rename({varname: renamevars[varname]})

            ds['time'].attrs['standard_name'] = 'time'    
            ds['time'].attrs['long_name'] = 'time'    
            ds[renamevars[varname]].attrs['standard_name'] = standard_names[renamevars[varname]]
            ds[renamevars[varname]].attrs['long_name'] = standard_names[renamevars[varname]]
            ds[renamevars[varname]].attrs['units'] = units[renamevars[varname]]
            ds[renamevars[varname]].encoding["dtype"] = "float32"

            ds['longitude'].encoding["dtype"] = "float32"
            ds['latitude'].encoding["dtype"] = "float32"

            starttime = np.datetime_as_string(ds.time.values[0], unit='D').replace('-','')
            endtime = np.datetime_as_string(ds.time.values[-1], unit='D').replace('-','')
            ds.rio.write_crs("epsg:4326", inplace="True")

            ds = ds.expand_dims(realization=[int(realization)])
            print(int(realization), ds.dims)

            outfile = os.path.join(
                out_dir, f'climateDT_{renamevars[varname]}_{experiment}_r{realization}_{areaname}_{starttime}_to_{endtime}.nc')
            
            if os.path.exists(outfile):    
                print(f"Skipping existing file: {outfile}")    
                continue
            
            ds.to_netcdf(outfile)


    # # Concat wind
    # ds1 = xr.open_dataset(os.path.join(out_dir, f'climateDT_10v_{experiment}_{areaname}_{starttime}_to_{endtime}.nc'))
    # ds2 = xr.open_dataset(os.path.join(out_dir, f'climateDT_10u_{experiment}_{areaname}_{starttime}_to_{endtime}.nc'))
    # ds_ws = xr.merge([ds1, ds2])

    # ds_ws = ds_ws.rename({'10v':'northward_wind', '10u':'eastward_wind'})
    # ds_ws['northward_wind'].attrs['standard_name'] = 'northward_wind'
    # ds_ws['eastward_wind'].attrs['standard_name'] = 'eastward_wind'
    
    # ds_ws.to_netcdf(os.path.join(out_dir, f'climateDT_windxy_{experiment}_{areaname}_{starttime}_to_{endtime}.nc'))

# %%

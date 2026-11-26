
#%%
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import os
from datetime import datetime
import numpy as np
import xarray as xr
from glob import glob
import rioxarray as rxr

#%%
SCRATCH = Path("/scratch-shared/dvertegaal1")
main_dir = SCRATCH / "data" / "ECMWF_ClimateDT"
# data_dir = os.path.join(main_dir, 'data', 'nc_Sofala')
data_dir = main_dir 

out_dir = main_dir / 'preprocessed' / 'Idai_full'
os.makedirs(out_dir, exist_ok=True)

#%%
# varnames = ['avg_tprate', 'msl','10u','10v', '2t', 'avg_sdswrf', 'tisr']
params_pl = {60: "pv", 129: "z", 131: "u", 132: "v"}
experiments = ["cont","hist","Tplus2.0K"]
date_range = "20190303/to/20190318"
date_string = date_range.replace("/to/", "_to_")
areaname = 'Idai_full'
# bbox = [32, -23, 36, -17]
bbox = [30, -22, 45, -12]     # Idai
# incl madagascar piece [30, -22, 45, -14.5] 
renamevars = {'avg_tprate':'tp', 'msl':'msl','10u':'10u','10v':'10v',
            #   '2t':'2t', 'avg_sdswrf':'ssrd', 'avg_tdswrf':'tisr', 'avg_snswrf': 'ssr'
              }


standard_names = {'pv': 'potential_vorticity', 'z': 'geopotential', 'u':'eastward_wind', 'v':'northward_wind'}
units = {'pv': 'K m**2 kg**-1 s**-1', 'z': 'm**2 s**-2', 'u':'m/s', 'v':'m/s'}

realizations = ['1', '2', '3', '4', '5']

#%%
for experiment in experiments[1:2]:
    for realization in realizations[0:1]:
    
            files = glob(str(data_dir / (f"climateDT_{'_'.join(params_pl.values())}_{experiment}_r{realization}_"f"{date_string}.nc")))
            print(files)
            ds = xr.open_mfdataset(files)
            ds.load()

            ds = ds.sortby('latitude', ascending=True) # sort latitude in ascending order
            # ds = ds.sel(latitude=slice(bbox[1], bbox[3]), longitude=slice(bbox[0], bbox[2]))

            # if varname == 'avg_tprate':
            #     ds[varname] = ds[varname] * 3600 # convert from kg/m2/s to mm/hr

            # ds = ds.rename({varname: renamevars[varname]})

            ds['time'].attrs['standard_name'] = 'time'    
            ds['time'].attrs['long_name'] = 'time'   
            ds['longitude'].encoding["dtype"] = "float32"
            ds['latitude'].encoding["dtype"] = "float32"

            for var in ["pv", "z", "u", "v"]:
                            ds[var].attrs["standard_name"] = standard_names[var]
                            ds[var].attrs["long_name"] = standard_names[var]
                            ds[var].attrs["units"] = units[var]
                            ds[var].encoding["dtype"] = "float32"

            starttime = np.datetime_as_string(ds.time.values[0], unit='D').replace('-','')
            endtime = np.datetime_as_string(ds.time.values[-1], unit='D').replace('-','')
            ds.rio.write_crs("EPSG:4326", inplace=True)

            # ds = ds.expand_dims(realization=[int(realization)])
            # print(int(realization), ds.dims)

            outfile = os.path.join(
                out_dir, f'climateDT_{'_'.join(params_pl.values())}_{experiment}_r{realization}_{areaname}_{starttime}_to_{endtime}.nc')
            
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

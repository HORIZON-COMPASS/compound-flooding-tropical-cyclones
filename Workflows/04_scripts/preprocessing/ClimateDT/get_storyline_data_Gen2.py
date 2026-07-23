'''
Script to access and view the storyline data from Climate DT.

Author: aleksand

Based on examples from: polytope-examples GitHub page

To run this script you need to install a Python environment with (pip):
polytope-client
earthkit==0.11.2
earthkit-data==0.17.0
covjsonkit==0.2.1
lxml 
conflator

You also need an access token in the home directory, file ~/.polytopeapirc
This is created by running the desp-authentication.py script from https://github.com/destination-earth-digital-twins/polytope-examples/

Available period: 2017 to today 

Guide to scenarios:
- "hist" - observed event
- "cont" - control scenario (mid-20th century)
- "Tplus2.0K" - 2 degree warming, SSP3.7.0

More info: https://destine.ecmwf.int/news/the-fast-development-of-destines-climate-change-adaptation-digital-twin/
'''

#%%
import earthkit.data
from earthkit.plots.interactive import Chart
from polytope.api import Client
from earthkit.plots.geo import domains
# from earthkit.geo import nearest_point_haversine
import matplotlib.pyplot as plt
import pandas as pd
import os
from datetime import datetime
import numpy as np
import time

import os
import tempfile
import shutil

print("TMPDIR:", os.environ.get("TMPDIR"))
print("tempdir:", tempfile.gettempdir())

for p in ["/tmp", tempfile.gettempdir(), "/var/lib/containers"]:
    try:
        usage = shutil.disk_usage(p)
        print(
            p,
            f"total={usage.total/1e9:.1f} GB",
            f"free={usage.free/1e9:.1f} GB"
        )
    except Exception as e:
        print(p, e)
        

def is_valid_earthkit_data(data):
    try:
        return len(data) > 0
    except Exception:
        return False

earthkit.data.config.set("maximum-cache-size", "1G") 

#%%

bbox_wflowsfincs = [32, -23, 36, -17] # from https://boundingbox.klokantech.com/, CSV format
bbox_dfm = [32, -28, 43, -9]
bbox_idai = [30, -22, 45, -12] 
areaname = 'Idai_full'
bbox_zoomin = [34.725665,-19.869092,35.060748,-19.630819] # Beira
areaname_zoomin = 'Beira'
p_ref = (-19.8, 34.84) # Beira
locname = 'Beira'
date = 20190315
hour = 500

params = [235055, 151, 165, 166, 167, 235035, 235053] # total precip, msl, u10, v10, 2t, ssrd, tisr
varnames = {235055:'avg_tprate', 151:'msl', 165:'10u', 166:'10v', 
            167:'2t', 235035:'avg_sdswrf', 235053:'avg_tdswrf'}
units = {235055:'kg m**-2 s**-1', 151:'Pa', 165:'m/s', 166:'m/s',
         167:'K', 235035:'W m**-2', 235053:'W m**-2'}
dates = ["20190301/to/20190309", "20190310/to/20190319", "20190320/to/20190329"]
experiments = ["cont","hist","Tplus2.0K"]

plots = False

# Parameters available at: https://confluence.ecmwf.int/display/DDCZ/Climate+DT+Phase+1+data+catalogue
request = {
     "activity": "story-nudging",
     "class": "d1",
     "dataset": "climate-dt",
     "generation": "2",
     "levtype": "sfc",
     "model": "ifs-fesom",
     "expver": "0001",
     "resolution": "high",
     "stream": "clte",
     "type": "fc",
     "time": "/".join([f"{h:02d}00" for h in range(0, 24)]),
     "grid": "0.1/0.1",
     "area": "-12/30/-22/45",
    #  'feature' : {                 
    #     "type" : "boundingbox",
    #     "points" : [[bbox[3], bbox[0]], [bbox[1], bbox[2]]],
    #     "axes" : ["latitude", "longitude"]
    #     }
 }

#%% Trying to access via GRIB file
for param in params:
    varname = varnames[param]
    request['param'] = param
    bbox = bbox_idai

    for experiment in experiments:
        request['experiment'] = experiment
        for realization in range(1,6):
            request['realization'] = str(realization)
            for dd, daterange in enumerate(dates):
                request['date'] = daterange

                data_file_nc = os.path.join('/p/11210471-001-compass/01_Data/ECMWF_ClimateDT/data', f'nc_{areaname}', f'climateDT_{varname}_{experiment}_{areaname}_{daterange.replace("/to/", "_to_")}_r{realization}.nc')
                
                if os.path.exists(data_file_nc):
                    print(f'{data_file_nc} already exists.')
                    continue
                else:
                    # Access data
                    max_retries = 3
                    for attempt in range(max_retries):
                        try:
                            data = earthkit.data.from_source("polytope", "destination-earth", request,
                                                             address='polytope.mn5.apps.dte.destination-earth.eu', stream=False)

                            # check data integrity
                            if is_valid_earthkit_data(data):
                                break
                            else:
                                print(f"⚠️ Empty dataset (attempt {attempt+1})")

                        except Exception as e:
                            print(f"⚠️ Download error (attempt {attempt+1}): {e}")

                        time.sleep(5)  # wait before retry

                    else:
                        # log failure and skip gracefully
                        with open("failed_requests.txt", "a") as f:
                            f.write(f"{request}\n")
                        print("Skipping request after retries:", request)
                        continue
                    
                    if len(data) == 0:
                        print("❌ No valid GRIB messages, skipping:", request)
                        continue

                    # Subset to bbox and save to disk
                    ds = data.to_xarray(add_earthkit_attrs=False)
                    ds_sel = ds.where((ds.longitude < bbox[2]) & (ds.longitude > bbox[0]) &
                                    (ds.latitude < bbox[3]) & (ds.latitude > bbox[1]), 
                                    drop=True)
                    ds_sel = ds_sel.rename({'forecast_reference_time':'time'})
                    ds_sel.attrs['experiment'] = experiment
                    ds_sel.attrs['bbox'] = bbox
                    
                    ds_sel.to_netcdf(data_file_nc)

                    time.sleep(2)
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

earthkit.data.config.set("maximum-cache-size", "5G") 

#%%

bbox_wflowsfincs = [32, -23, 36, -17] # from https://boundingbox.klokantech.com/, CSV format
bbox_dfm = [32, -28, 43, -9]
bbox_idai = [30, -22, 45, -14.5] 
areaname = 'Idai_full'
bbox_zoomin = [34.725665,-19.869092,35.060748,-19.630819] # Beira
areaname_zoomin = 'Beira'
p_ref = (-19.8, 34.84) # Beira
locname = 'Beira'
date = 20190315
hour = 500

params = [235055, 151, 165, 166, 167, 235035]#, 212] # total precip, msl, u10, v10, 2t, ssrd, tisr
varnames = {235055:'avg_tprate', 151:'msl', 165:'10u', 166:'10v', 
            167:'2t', 235035:'avg_sdswrf'}#, 212:'tisr'}
units = {235055:'kg m**-2 s**-1', 151:'Pa', 165:'m/s', 166:'m/s',
         167:'K', 235035:'W m**-2'}#, 212:'J/m2'}
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
    #  'feature' : {                 
    #     "type" : "boundingbox",
    #     "points" : [[bbox[3], bbox[0]], [bbox[1], bbox[2]]],
    #     "axes" : ["latitude", "longitude"]
    #     }
 }
#%%

# for param in params[:2]:
#     for experiment in experiments:
#         request['date'] = dates
#         request['param'] = param
#         request["experiment"] = experiment
#         data_file = f"data/climate-dt-story-nudging_{experiment}_{dates.replace('/','_')}_{param}.covjson"

#         if LIVE_REQUEST:
#             data = earthkit.data.from_source("polytope", "destination-earth", request, address="polytope.lumi.apps.dte.destination-earth.eu", stream=False)
#             data.to_target("file", data_file)
#         else:
#             data = earthkit.data.from_source("file", data_file) 

# #%%
# data_file = f"data/climate-dt-story-nudging_{experiments[1]}_{dates.replace('/','_')}_{params[0]}.covjson"

# data = earthkit.data.from_source("file", data_file) 

# da = data.to_xarray(add_earthkit_attrs=False)

# data_file_nc = data_file.replace('.covjson','.nc')
# da.to_netcdf(data_file_nc)

#%% Trying to access via GRIB file

#request.pop('feature')

for param in params[1:4]:

    varname = varnames[param]
    request['param'] = param

    if param in [151, 165, 166]:
        bbox = bbox_idai
    else:
        bbox = bbox_wflowsfincs

    if varname=='tp':
        factor = 1000 # to convert from m to mm
    else:
        factor = 1

    for experiment in experiments[1:2]:
        request['experiment'] = experiment

        for realization in range(1,6): # range(1,11):
            request['realization'] = str(realization)

            for dd, daterange in enumerate(dates):
                request['date'] = daterange

                data_file_nc = os.path.join('/p/11210471-001-compass/01_Data/ECMWF_ClimateDT/data', f'nc_{areaname}', f'climateDT_{varname}_{experiment}_{areaname}_{daterange.replace("/to/", "_to_")}_r{realization}.nc')
                
                if os.path.exists(data_file_nc):
                    print(f'{data_file_nc} already exists.')
                    continue
                else:
                    # Access data
                    data = earthkit.data.from_source("polytope", "destination-earth", 
                                                    request, 
                                                    address='polytope.mn5.apps.dte.destination-earth.eu',
                                                    stream=False)
                    
                    # check temp path where the data is downloaded to:
                    #data.path 

                    # Save to disk - global file
                    #data_file_nc = os.path.join('data', f'climateDT_{experiment}_{varnames[param]}_{areaname}.nc')
                    #data_latlon.to_netcdf(data_file_nc)

                    # Subset to bbox and save to disk
                    ds = data.to_xarray(add_earthkit_attrs=False)
                    ds_sel = ds.where((ds.longitude < bbox[2]) & (ds.longitude > bbox[0]) &
                                    (ds.latitude < bbox[3]) & (ds.latitude > bbox[1]), 
                                    drop=True)
                    ds_sel = ds_sel.rename({'forecast_reference_time':'time'})
                    ds_sel.attrs['experiment'] = experiment
                    ds_sel.attrs['bbox'] = bbox
                    
                    ds_sel.to_netcdf(data_file_nc)

#                 # Plotting single timestep on a map
#                 if plots & dd==2:
#                     domain_selection = domains.Domain.from_bbox(
#                         bbox=[bbox[0], bbox[2], bbox[1], bbox[3]],
#                         name=areaname)
#                     chart = earthkit.plots.Map(domain=domain_selection)
#                     if param==228:
#                         chart.contourf(data.sel(dataDate=date,dataTime=hour), 
#                                     units=units[param], auto_style=True)
#                         chart.grid_points(data.sel(dataDate=20190317,dataTime=0))
#                     elif param==151:
#                         chart.quickplot(data.sel(dataDate=date,dataTime=hour), 
#                                         style = earthkit.plots.styles.Contour(linecolors='blue',levels={"step": 4},
#                                                                             units="hPa", legend_style=None,labels=True))
#                     elif param==165 or param==166:
#                         chart.quickplot(data.sel(dataDate=date,dataTime=hour))
#                     chart.coastlines()
#                     chart.borders()
#                     chart.gridlines()
#                     chart.legend()
#                     chart.title(f"Variable: {varname}, experiment {experiment}\n time: {date} {hour:04.0f}")
#                     chart.fig.savefig(os.path.join('data', 'figures', f'climateDT_{varname}_{experiment}_{date}_{hour:04.0f}_{areaname}_r{realization}.png'), format='png')


#                 # Plot regridded data sample for a single timestep
#                 if plots & dd==2:
#                     chart = earthkit.plots.Map(domain=domain_selection)
#                     if param==228:
#                         chart.contourf(data_latlon.sel(dataDate=date,dataTime=hour), 
#                                     units=units[param], auto_style=True)
#                         chart.grid_points(data_latlon.sel(dataDate=20190317,dataTime=0))
#                     elif param==151:
#                         chart.quickplot(data_latlon.sel(dataDate=date,dataTime=hour), 
#                                         style = earthkit.plots.styles.Contour(linecolors='blue',levels={"step": 4},
#                                                                             units="hPa", legend_style=None,labels=True))
#                     elif param==165 or param==166:
#                         chart.quickplot(data_latlon.sel(dataDate=date,dataTime=hour))
#                     chart.coastlines()
#                     chart.borders()
#                     chart.gridlines()
#                     chart.legend()
#                     chart.title(f"Variable: {varname}, experiment {experiment}\n time: {date} {hour:04.0f} - Interpolated")
#                     chart.fig.savefig(os.path.join('data', 'figures', f'climateDT_{varname}_{experiment}_{date}_{hour:04.0f}_{areaname} - interpolated_r{realization}.png'), format='png')

#                 # extracting a single location - plotting timeseries
#                 # example of extracting directly from lazy-loaded .grib file
#                 # latlon = data.to_latlon()
#                 # lat = latlon["lat"]
#                 # lon = latlon["lon"]
#                 # idx, dist = nearest_point_haversine(p_ref, (lat, lon))
                
#                 # v = data_sel.values[:,idx].squeeze()
#                 # v=v * factor
#                 # t = data_sel.metadata("valid_datetime")
#                 # time = pd.to_datetime(t)
#                 # # df = pd.DataFrame({
#                 # #     "time": time,
#                 # #     "value": v
#                 # # })
#                 # fig, ax = plt.subplots()
#                 # ax.plot(time,v)
#                 # ax.set_ylim([0,40])
#                 # ax.set_title(f'Variable: {varname}, location {p_ref}, experiment {experiment}')
#                 # fig.savefig(os.path.join('data', 'figures', f'climateDT_{experiment}_{varname}_timeseries_{locname}.png'))

#                 # extracting mean over a small area (bbox_small)
#                 # extracting from the netCDF file
#                 if plots & dd==2:
#                     ds_sel_small = ds_sel.where((ds_sel.longitude < bbox_zoomin[2]) &
#                                     (ds_sel.longitude > bbox_zoomin[0]) &
#                                     (ds_sel.latitude < bbox_zoomin[3]) &
#                                     (ds_sel.latitude > bbox_zoomin[1]), drop=True)
                    
#                     fig, ax = plt.subplots()
#                     (ds_sel_small.mean(dim=['latitude','longitude'])[varname]*factor).plot(ax=ax)
#                     if param==228:
#                         ax.set_ylim([0,40])
#                     elif param==151:
#                         ax.set_ylim([97000, 102000])
#                     elif param==165 or param==166:
#                         ax.set_ylim([-20, 20])
#                     if param==228:
#                         ax.set_title(f'Variable: {varname} \n area {areaname_zoomin} \n experiment {experiment} \n accumulated sum: {ds_sel_small['tp'].sum().item()*1000:.1f} mm')
#                     else:
#                         ax.set_title(f'Variable: {varname} \n area {areaname_zoomin} \n experiment {experiment}')
#                     fig.savefig(os.path.join('data', 'figures', f'climateDT_{varname}_{experiment}_timeseries_{areaname_zoomin}_r{realization}.png'))

#                 # Plot accumulated rainfall - average value per hour
#                 if plots and param == 228 and dd==2:
#                     datestr = datetime.strptime(str(date), "%Y%m%d").strftime("%Y-%m-%d")
#                     ds_sel_1day = ds_sel.sel(time=datestr)

#                     ds_sel_av = ds_sel_1day.sum(dim='time')/24
#                     ds_sel_av['tp'].attrs = ds_sel_1day['tp'].attrs

#                     # Plot
#                     chart = earthkit.plots.Map(domain=domain_selection)
#                     cc = chart.contourf(ds_sel_av['tp'], 
#                                 units="mm",
#                                 auto_style=True)
#                     chart.grid_points(data_latlon.sel(dataDate=20190317,dataTime=0))
#                     chart.coastlines(linewidth=2)
#                     chart.borders()
#                     chart.gridlines()
#                     chart.legend()
#                     chart.title(f"Variable: {varname}, date {date} - hourly average \n experiment {experiment}")
#                     chart.fig.savefig(os.path.join('data', 'figures', f'climateDT_{varname}_{experiment}_{date}_{areaname}_hourly_average_r{realization}.png'), format='png')

#                 ds_sel.close(); del ds_sel
#                 ds.close(); del ds
#                 del data_latlon
#                 del data


# # %%

# %%

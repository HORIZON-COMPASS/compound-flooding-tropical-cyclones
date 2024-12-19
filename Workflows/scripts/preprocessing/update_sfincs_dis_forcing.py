# %%
from datetime import datetime as datetime
from datetime import timedelta
from os.path import basename, join
from shutil import copy

import pandas as pd
import xarray as xr
from hydromt.config import configread
from hydromt.log import setuplog
from hydromt_wflow import WflowModel
from hydromt_sfincs import SfincsModel
from hydromt_sfincs.sfincs_input import SfincsInput

# %%
# model and data paths/
logger = setuplog("update", "./hydromt.log", log_level=10)
if "snakemake" in locals():
    sfincs_model_folder = snakemake.params.dir_run_with_forcing
    wflow_root = snakemake.params.wflow_root_forcing
    data_cats = snakemake.params.data_cats
else:
    curdir = r'c:\Git_repos\COMPASS\Workflows'
    wflow_root = r"p:\11210471-001-compass\03_Runs\quelimane\Freddy2\era5_hourly\wflow"
    event =  "freddy2"
    sfincs_model_folder = r"p:\11210471-001-compass\03_Runs\quelimane\Freddy2\era5_hourly\sfincs"
    data_cats = [
            join(curdir, "data_catalogs", "datacatalog_general___linux.yml"), 
            join(curdir, "data_catalogs", "datacatalog_SFINCS_coastal_coupling___linux.yml"), 
            join(curdir, "data_catalogs", "datacatalog_SFINCS_obspoints___linux.yml")
        ]



#%%
mod = WflowModel(
    root=join(wflow_root, 'events'),
    data_libs=data_cats,
    mode="r",
    logger=logger,
)
mod.read()
#%%
# reftime_object = datetime.strptime(ref_time, "%Y-%m-%dT%H:%M:%S")

inp = SfincsInput.from_file(join(sfincs_model_folder,"sfincs.inp"))
config = inp.to_dict()


# Write dis file to sfincs event model folder
reftime_object = config["tref"]
df = mod.results['netcdf']['Q'].to_pandas()
df.index = (df.index - reftime_object).total_seconds()
df.to_csv(
    join(sfincs_model_folder, "sfincs.dis"),
    sep=" ",
    header=False,
)


    # precip = mod.data_catalog.get_rasterdataset(
    #     'era5_hourly',
    #     geom = mod.region,
    #     buffer = 5,
    #     variables = 'precip',
    #     time_tuple = (config['tstart'], config['tstop']),                     
    # ).to_dataset().rename({
    #     'precip' : 'Precipitation',
    #     'longitude' : 'x',
    #     'latitude' : 'y'
    #     })

    # #Convert time units to minutes
    # encoding = dict(
    #             time={"units": f"minutes since 1900-01-01", "dtype": "float64", '_FillValue': None}
    #         )

    # precip.to_netcdf(join(sfincs_model_folder,"precip.nc"), encoding = encoding)
    # Update inp file
    # config.pop('amprfile', None)
config.update({'netamprfile': "precip.nc"})
config.update({"disfile": "sfincs.dis"})
config.update({"srcfile": "sfincs.src"})
inp = SfincsInput.from_dict(config)
inp.write(inp_fn=join(sfincs_model_folder, "sfincs.inp"))

# %%

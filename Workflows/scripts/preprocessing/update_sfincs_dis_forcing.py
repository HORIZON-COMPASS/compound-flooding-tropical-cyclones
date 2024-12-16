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
# if "snakemake" in locals():
#     wflow_root = snakemake.params.wflow_root
#     event = snakemake.params.event
#     cosmos_event_root = snakemake.params.cosmos_event_root
#     use_case = snakemake.params.use_case
#     sfincs_root = snakemake.params.sfincs_root
#     sf_model_name = snakemake.params.sfincs_event_name
#     meteo_option = snakemake.params.meteo_option
#     meteo_fn = snakemake.params.meteo_fn
# else
wflow_root = "c:/Git_repos/COMPASS/Wflow/wflow_sofala"
event =  "freddy2"
sfincs_model_folder = "c:/Git_repos/COMPASS/SFINCS/sfincs_sofala/computations/sfincs_CLI_Freddy2"




#%%
mod = WflowModel(
    root=join(wflow_root, "event", event),
    data_libs=["deltares_data", r"c:/Git_repos/COMPASS/SFINCS/datacatalog_general.yml"],
    mode="r",
    logger=logger,
)

# reftime_object = datetime.strptime(ref_time, "%Y-%m-%dT%H:%M:%S")
for var in mod.results['netcdf'].data_vars:
    inp = SfincsInput.from_file(join(sfincs_model_folder,"sfincs.inp"))
    config = inp.to_dict()


    # Write dis file to sfincs event model folder
    reftime_object = config["tref"]
    df = mod.results['netcdf'][var].to_pandas()
    df.index = (df.index - reftime_object).total_seconds()
    df.to_csv(
        join(sfincs_model_folder, "sfincs.dis"),
        sep=" ",
        header=False,
    )


    precip = mod.data_catalog.get_rasterdataset(
        'era5_hourly',
        geom = mod.region,
        buffer = 5,
        variables = 'precip',
        time_tuple = (config['tstart'], config['tstop']),                     
    ).to_dataset().rename({
        'precip' : 'Precipitation',
        'longitude' : 'x',
        'latitude' : 'y'
        })


    #Convert time units to minutes
    encoding = dict(
                time={"units": f"minutes since 1900-01-01", "dtype": "float64", '_FillValue': None}
            )

    precip.to_netcdf(join(sfincs_model_folder,"precip.nc"), encoding = encoding)
    # Update inp file
    config.pop('amprfile', None)
    config.update({'netamprfile': "precip.nc"})
    config.update({"disfile": "sfincs.dis"})
    config.update({"srcfile": "sfincs.src"})
    inp = SfincsInput.from_dict(config)
    inp.write(inp_fn=join(sfincs_model_folder, "sfincs.inp"))

# %%

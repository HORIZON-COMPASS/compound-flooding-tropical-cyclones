# %%
from os.path import basename, join

import pandas as pd
from hydromt.config import configread
from hydromt.log import setuplog
from hydromt_wflow import WflowModel
from datetime import datetime as datetime
from datetime import timedelta

# %%
# model and data paths/
logger = setuplog("update", "./hydromt.log", log_level=10)
# wflow_root = r"p:\11208614-de-370a\01_models\Humber\wflow"
if "snakemake" in locals():
    wflow_root_noforcing = snakemake.params.wflow_root_noforcing
    wflow_root_forcing = snakemake.params.wflow_root_forcing
    start_time = snakemake.params.start_time
    end_time = snakemake.params.end_time
    data_cat = snakemake.params.data_cat
    meteo_fn = snakemake.params.forcing

# else:
#     script_root = r"p:\11208614-de-370a\02_scripts\wflow_sfincs_snake"
#     wflow_root = r"p:\11208614-de-370a\01_models\Reunion\wflow"
#     use_case = "Reunion"
#     start_time = "2024-01-12T00:00:00"
#     end_time = "2024-01-19T00:00:00"
#     event = "2024-01-12"
#     meteo_fn = "era5_hourly_local"
#     meteo_option = "ERA5"

# %% Setup forcing
start_time_object = datetime.strptime(start_time, "%Y%m%d %H%M%S") - timedelta(days=2)
end_time_object = datetime.strptime(end_time, "%Y%m%d %H%M%S")
start_time = datetime.strftime(start_time_object, "%Y-%m-%dT%H:%M:%S")
end_time = datetime.strftime(end_time_object, "%Y-%m-%dT%H:%M:%S")

mod = WflowModel(
    root=wflow_root_noforcing,
    data_libs=[data_cat],
    mode="r",
    logger=logger,
)
mod.read()

opt = {
    "setup_config": {
        "starttime": start_time,
        "endtime": end_time,
        "timestepsecs": 3600,
        "model.reinit": False,
        "input.path_static": join("..","staticmaps.nc"),
        "input.path_forcing":"inmaps.nc",
    },
    "setup_precip_forcing": {
        "precip_fn": meteo_fn,
        "precip_clim_fn": None,
    },
    "setup_temp_pet_forcing": {
        #"temp_pet_fn": meteo_fn,
        "temp_pet_fn": "era5_hourly_zarr", #"era5_hourly_zarr"
        "press_correction": True,
        "temp_correction": True,
        "pet_method": "debruin",
        "skip_pet": False,
    },
}

mod.set_root(join(wflow_root_forcing, "events"), mode="w+")
mod.update(opt=opt, write=False)
mod.write_forcing()
mod.write_config()
# %%

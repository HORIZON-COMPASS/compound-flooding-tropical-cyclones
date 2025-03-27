# %%
from datetime import datetime as datetime
from datetime import timedelta
from os.path import basename, join

import pandas as pd
from hydromt.config import configread
from hydromt.log import setuplog
from hydromt_wflow import WflowModel

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
else:
    precip_forcing = "era5_hourly"
    CF_rain = 0
    CF_rain_txt = "0"
    wflow_root_noforcing = "p:/11210471-001-compass/02_Models/sofala/Idai/wflow_test"
    wflow_root_forcing = f"p:/11210471-001-compass/03_Runs/sofala/Idai/wflow_test/event_precip_{precip_forcing}_CF{CF_rain_txt}"
    start_time = "20190309 000000"
    end_time = "20190325 060000"
    data_cat = [
        '../../../03_data_catalogs/datacatalog_general.yml',
        '../../../03_data_catalogs/datacatalog_CF_forcing.yml',
    ] 

# %% Setup forcing Warmup
mod = WflowModel(
    root=wflow_root_noforcing,
    data_libs=data_cat,
    mode="r",
    logger=logger,
)
mod.read()
start_time_object = datetime.strptime(start_time, "%Y%m%d %H%M%S") - timedelta(days=2) #Start wflow 2 days before sfincs
start_time_warmup = datetime.strftime(
    start_time_object - timedelta(days=365), "%Y-%m-%dT%H:%M:%S"
)
end_time_warmup = datetime.strftime(start_time_object, "%Y-%m-%dT%H:%M:%S")
opt = {
    "setup_config": {
        "starttime": start_time_warmup,
        "endtime": end_time_warmup,
        "timestepsecs": 86400,
        "model.reinit": True,
        "state.path_output": join(
           "..", "..", "events", "instate", "instates.nc"
        ),
        "input.path_static": join("..","staticmaps.nc"),
        "input.path_forcing":"inmaps.nc",
    },
    "setup_precip_forcing": {
        "precip_fn": "era5_daily",
        "precip_clim_fn": None,
        "chunksize": 10,
    },
    "setup_temp_pet_forcing": {
        "temp_pet_fn": "era5_daily",
        "press_correction": True,
        "temp_correction": True,
        "dem_forcing_fn": "era5_orography",
        "pet_method": "debruin",
        "skip_pet": False,
        "chunksize": 10,
    },
}

mod.set_root(join(wflow_root_forcing, "warmup"), mode="w+")
mod.setup_config(**opt["setup_config"])
mod.setup_precip_forcing(**opt["setup_precip_forcing"])
mod.setup_temp_pet_forcing(**opt["setup_temp_pet_forcing"])
mod.write_forcing(chunksize=10)
mod.write_grid()
mod.write_config()

# %%

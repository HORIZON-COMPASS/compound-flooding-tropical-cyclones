# %%
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
    wflow_root = snakemake.params.wflow_root
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

mod = WflowModel(
    root=wflow_root,
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
        "input.path_static": "..\staticmaps.nc",
    },
    "setup_precip_forcing": {
        "precip_fn": meteo_fn,
        "precip_clim_fn": None,
    },
    "setup_temp_pet_forcing": {
        "temp_pet_fn": meteo_fn,
        "press_correction": True,
        "temp_correction": True,
        "pet_method": "debruin",
        "skip_pet": False,
    },
}

mod.set_root(join(wflow_root, "events"), mode="w+")
mod.update(opt=opt, write=False)
mod.write_forcing(fn_out=join(mod.root, "inmaps.nc"))
mod.write_config()
# %%

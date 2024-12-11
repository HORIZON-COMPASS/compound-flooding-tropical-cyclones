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
    wflow_root = snakemake.params.wflow_root
    start_time = snakemake.params.start_time
    end_time = snakemake.params.end_time
    data_cat = snakemake.params.data_cat
# else:
#     script_root = "p:/11208614-de-370a/02_scripts/wflow_sfincs_snake"
#     wflow_root = (
#         "p:/11208614-de-370a/01_models/Philippines/wflow/wflow_philippines_all_islands"
#     )
#     use_case = "Philippines"
#     start_time = "2013-11-06T00:00:00"
#     end_time = "2013-11-12T00:00:00"
#     event = "2013-11-06"

# %% Setup forcing Warmup
mod = WflowModel(
    root=wflow_root,
    data_libs=[data_cat],
    mode="r",
    logger=logger,
)
mod.read()
start_time_object = datetime.strptime(start_time, "%Y-%m-%d %H:%M:%S")
start_time_warmup = datetime.strftime(
    start_time_object - timedelta(days=365), "%Y-%m-%d %H:%M:%S"
)
end_time_warmup = start_time
opt = {
    "setup_config": {
        "starttime": start_time_warmup,
        "endtime": end_time_warmup,
        "timestepsecs": 86400,
        "model.reinit": True,
        "state.path_output": join(
            wflow_root, "events", "instate", "instates.nc"
        ),
        "input.path_static": "..\staticmaps.nc",
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

mod.set_root(join(wflow_root, "warmup"), mode="w+")
mod.setup_config(**opt["setup_config"])
mod.setup_precip_forcing(**opt["setup_precip_forcing"])
mod.setup_temp_pet_forcing(**opt["setup_temp_pet_forcing"])
mod.write_forcing(fn_out=join(mod.root, "inmaps.nc"), chunksize=10)
mod.write_config()

# %%

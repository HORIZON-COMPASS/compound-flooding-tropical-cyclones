# %%
# Add warmup forcing to a Wflow base model — HydroMT v1 / hydromt_wflow 1.0.x (WflowSbmModel).
# Migrated from the v0 set_root/setup_*/write_* API: in v1 we relocate + apply the setup steps
# via mod.update(model_out=..., steps=..., write=True). The setup_config keys use the Wflow.jl
# v1 TOML layout (time.* section); see migration_ref_files/wflow_build.yml.
from datetime import datetime as datetime
from datetime import timedelta
from os.path import basename, join
import pandas as pd
from hydromt_wflow import WflowSbmModel

# %%
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
    wflow_root_noforcing = "/p/11210471-001-compass/02_Models/sofala/Idai/wflow"
    wflow_root_forcing = f"/p/11210471-001-compass/03_Runs/sofala/Idai/wflow/event_precip_{precip_forcing}_CF{CF_rain_txt}"
    start_time = "20190309 000000"
    end_time = "20190325 060000"
    data_cat = [
        '../../../03_data_catalogs/datacatalog_general_v1___linux.yml',
        '../../../03_data_catalogs/datacatalog_CF_forcing_v1___linux.yml',
    ]

# %% Read the base model
mod = WflowSbmModel(root=wflow_root_noforcing, data_libs=data_cat, mode="r")
mod.read()

start_time_object = datetime.strptime(start_time, "%Y%m%d %H%M%S") - timedelta(days=2)  # start 2 days before sfincs
start_time_warmup = datetime.strftime(start_time_object - timedelta(days=20), "%Y-%m-%dT%H:%M:%S")
end_time_warmup = datetime.strftime(start_time_object, "%Y-%m-%dT%H:%M:%S")

# %% v1 steps (setup_config takes a `data` dict; time settings live under [time] in Wflow.jl v1)
steps = [
    {"setup_config": {"data": {
        "time.starttime": start_time_warmup,
        "time.endtime": end_time_warmup,
        "time.timestepsecs": 86400,
        "model.cold_start__flag": True,   # v0 model.reinit=True (cold start); renamed in Wflow.jl v1
        "state.path_output": join("..", "..", "events", "instate", "instates.nc"),
        "input.path_static": join("..", "staticmaps.nc"),
        "input.path_forcing": "inmaps.nc",
    }}},
    {"setup_precip_forcing": {"precip_fn": "era5_daily", "precip_clim_fn": None, "chunksize": 10}},
    {"setup_temp_pet_forcing": {
        "temp_pet_fn": "era5_daily",
        "press_correction": True,
        "temp_correction": True,
        "dem_forcing_fn": "era5_orography",
        "pet_method": "debruin",
        "skip_pet": False,
        "chunksize": 10,
    }},
]

# Relocate to the warmup run dir and write the full model (forcing + staticmaps + config)
mod.update(model_out=join(wflow_root_forcing, "warmup"), steps=steps, write=True, forceful_overwrite=True)

# %%

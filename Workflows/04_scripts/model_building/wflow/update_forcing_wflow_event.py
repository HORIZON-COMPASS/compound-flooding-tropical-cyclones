# %%
# Add event forcing to a Wflow base model — HydroMT v1 / hydromt_wflow 1.0.x (WflowSbmModel).
# Migrated from the v0 set_root/opt/write_* API: in v1 we relocate + apply the steps via
# mod.update(model_out=..., steps=..., write=False) then write forcing + config (the event reuses
# the warmup staticmaps via input.path_static, so we do NOT rewrite the grid). setup_config uses
# the Wflow.jl v1 TOML keys (time.* section; model.reinit -> model.cold_start__flag).
from os.path import join
from hydromt_wflow import WflowSbmModel
from datetime import datetime as datetime
from datetime import timedelta

# %%
if "snakemake" in locals():
    tc_name              = snakemake.params.tc_name
    wflow_root_noforcing = snakemake.params.wflow_root_noforcing
    wflow_root_forcing   = snakemake.params.wflow_root_forcing
    start_time           = snakemake.params.start_time
    end_time             = snakemake.params.end_time
    data_cat             = snakemake.params.data_cat
    precip_forcing       = snakemake.wildcards.precip_forcing
    CF_rain              = float(snakemake.wildcards.CF_rain)
    CF_rain_txt          = snakemake.wildcards.CF_rain
    meteo_fn             = snakemake.params.meteo_forcing
else:
    tc_name              = "Idai"
    precip_forcing       = "era5_hourly"
    CF_rain              = -7
    CF_rain_txt          = f"{CF_rain}"
    wflow_root_noforcing = "/p/11210471-001-compass/02_Models/sofala/Idai/wflow"
    wflow_root_forcing   = f"/p/11210471-001-compass/03_Runs/sofala/Idai/wflow/event_precip_{precip_forcing}_CF{CF_rain_txt}"
    start_time           = "20190309 000000"
    end_time             = "20190325 060000"
    data_cat             = ['../../../03_data_catalogs/datacatalog_general_v1___linux.yml',
                            '../../../03_data_catalogs/datacatalog_CF_forcing_v1___linux.yml',
                            ]
    meteo_fn             = "era5_hourly"

# %%
# Setup forcing time window
start_time_object = datetime.strptime(start_time, "%Y%m%d %H%M%S") - timedelta(days=2)
end_time_object = datetime.strptime(end_time, "%Y%m%d %H%M%S")
start_time = datetime.strftime(start_time_object, "%Y-%m-%dT%H:%M:%S")
end_time = datetime.strftime(end_time_object, "%Y-%m-%dT%H:%M:%S")

# %% Read the base model
mod = WflowSbmModel(root=wflow_root_noforcing, data_libs=data_cat, mode="r")
mod.read()

# %% Build the v1 steps list
steps = [
    {"setup_config": {"data": {
        "time.starttime": start_time,
        "time.endtime": end_time,
        "time.timestepsecs": 3600,
        "model.cold_start__flag": False,   # v0 model.reinit=False (warm start from warmup states)
        "input.path_static": join("..", "staticmaps.nc"),
        "input.path_forcing": "inmaps.nc",
    }}},
    {"setup_temp_pet_forcing": {
        "temp_pet_fn": meteo_fn,
        "press_correction": True,
        "temp_correction": True,
        "pet_method": "debruin",
        "skip_pet": False,
    }},
]

# Rainfall forcing (factual uses the raw product; counterfactual uses the CF-shifted product)
if CF_rain is None:
    print("Error: CF_rain value not found")
elif CF_rain == 0:
    steps.append({"setup_precip_forcing": {"precip_fn": precip_forcing, "precip_clim_fn": None}})
else:
    steps.append({"setup_precip_forcing": {"precip_fn": f'{precip_forcing}_CF{CF_rain_txt}_{tc_name}', "precip_clim_fn": None}})

# %%
# Relocate to the events run dir, apply steps without writing, then write forcing + config only
# (the event reuses the warmup staticmaps, so we deliberately do not rewrite the grid).
mod.update(model_out=join(wflow_root_forcing, "events"), steps=steps, write=False, forceful_overwrite=True)
mod.forcing.write()
mod.config.write()
# %%

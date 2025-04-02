# %%
from os.path import join
from hydromt.log import setuplog
from hydromt_wflow import WflowModel
from datetime import datetime as datetime
from datetime import timedelta

# %%
# model and data paths/
logger = setuplog("update", "./hydromt.log", log_level=10)

if "snakemake" in locals():
    tc_name              = snakemake.wildcards.runname
    wflow_root_noforcing = snakemake.params.wflow_root_noforcing
    wflow_root_forcing   = snakemake.params.wflow_root_forcing
    start_time           = snakemake.params.start_time
    end_time             = snakemake.params.end_time
    data_cat             = snakemake.params.data_cat
    precip_forcing       = snakemake.wildcards.precip_forcing
    CF_rain              = float(snakemake.wildcards.CF_rain)
    CF_rain_txt          = snakemake.wildcards.CF_rain
else:
    tc_name              = "Idai"
    precip_forcing       = "era5_hourly"
    CF_rain              = -7
    CF_rain_txt          = f"{CF_rain}"
    wflow_root_noforcing = "p:/11210471-001-compass/02_Models/sofala/Idai/wflow"
    wflow_root_forcing   = f"p:/11210471-001-compass/03_Runs/sofala/Idai/wflow/event_precip_{precip_forcing}_CF{CF_rain_txt}"
    start_time           = "20190309 000000"
    end_time             = "20190325 060000"
    data_cat             = ['../../../03_data_catalogs/datacatalog_general.yml',
                            '../../../03_data_catalogs/datacatalog_CF_forcing.yml',
                            ] 

# %% 
# Setup forcing
start_time_object = datetime.strptime(start_time, "%Y%m%d %H%M%S") - timedelta(days=2)
end_time_object = datetime.strptime(end_time, "%Y%m%d %H%M%S")
start_time = datetime.strftime(start_time_object, "%Y-%m-%dT%H:%M:%S")
end_time = datetime.strftime(end_time_object, "%Y-%m-%dT%H:%M:%S")

mod = WflowModel(
    root=wflow_root_noforcing,
    data_libs=data_cat,
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
    "setup_temp_pet_forcing": {
        "temp_pet_fn": "era5_hourly",
        "press_correction": True,
        "temp_correction": True,
        "pet_method": "debruin",
        "skip_pet": False,
    },
}

# Add rainfall forcing
if CF_rain is None:
    print(f"Error: CF_rain value not found")
elif CF_rain == 0:
    opt["setup_precip_forcing"] = {
        "precip_fn": precip_forcing,
        "precip_clim_fn": None,  # Use a different forcing file
    }
else:
    opt["setup_precip_forcing"] = {
        "precip_fn": f'{precip_forcing}_CF{CF_rain_txt}_{tc_name}',
        "precip_clim_fn": None,  # Use a different forcing file
    }

#%%
mod.set_root(join(wflow_root_forcing, "events"), mode="w+")
mod.update(opt=opt, write=False)
mod.write_forcing()
mod.write_config()
# %%

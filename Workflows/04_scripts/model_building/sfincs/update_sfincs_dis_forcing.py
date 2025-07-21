# %%
from datetime import datetime as datetime
from os.path import join
import os
from hydromt.log import setuplog
from hydromt_wflow import WflowModel
from hydromt_sfincs.sfincs_input import SfincsInput
import pandas as pd

# %%
# model and data paths/
logger = setuplog("update", "./hydromt.log", log_level=10)
if "snakemake" in locals():
    sfincs_model_folder   = snakemake.params.dir_run_with_forcing
    wflow_root            = snakemake.params.wflow_root_forcing
    data_cats             = snakemake.params.data_cats
    wflow_dis_no_bankfull = snakemake.input.wflow_dis_no_bankfull
    wflow_base            = snakemake.params.wflow_base
else:
    curdir              = '../../../'
    region              = "sofala"
    tc_name             = "Idai"
    precip_forcing      = "era5_hourly_zarr"
    wind_forcing        = 'era5_hourly_spw_IBTrACS'
    tidemodel           = 'GTSMv41' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap, tpxo80_opendap
    CF_rain_txt         = "0"
    CF_SLR_txt          = "0"
    CF_wind_txt         = "0"
    wflow_root          = f"p:/11210471-001-compass/03_Runs/{region}/{tc_name}/wflow/event_precip_{precip_forcing}_CF{CF_rain_txt}"
    wflow_base          = f"p:/11210471-001-compass/02_Models/{region}/{tc_name}/wflow"
    sfincs_model_folder = f"p:/11210471-001-compass/03_Runs/{region}/{tc_name}/sfincs/event_tp_{precip_forcing}_CF{CF_rain_txt}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}"
    data_cats           = [
        join(curdir, "03_data_catalogs", "datacatalog_general.yml"), 
        join(curdir, "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling.yml"), 
        join(curdir, "03_data_catalogs", "datacatalog_SFINCS_obspoints.yml"),
        join(curdir, "03_data_catalogs", "datacatalog_CF_forcing.yml")
        ]
    wflow_dis_no_bankfull = f"{wflow_root}/events/run_default/wflow_dis_no_qbankfull.csv"

#%%
# Read config rile from sfincs model with coastal and meteo forcing
inp = SfincsInput.from_file(join(sfincs_model_folder,"sfincs.inp"))
config = inp.to_dict()

#%%
# Write dis file to sfincs event model folder
reftime_object = config["tref"]

#%%
# mod = WflowModel(
#     root=join(wflow_root, 'events'),
#     data_libs=data_cats,
#     mode="r",
#     logger=logger,
# )
# mod.read()
# #%%
# df_check = mod.results['netcdf']['Q'].to_pandas()
# df_check.index = (df_check.index - reftime_object).total_seconds()
# df_check.to_csv(
#     join(sfincs_model_folder, "sfincs_compare.dis"),
#     sep=" ",
#     header=False,
# )

# Read the original wflow gauge order to match that of sfincs
mod_ini = WflowModel(root=os.path.join(wflow_base), mode="r+", config_fn=os.path.join(wflow_base, "wflow_sbm.toml"))
q_locs = mod_ini.geoms["gauges_locs"]

# Read wflow qbankfull corrected discharge
df = pd.read_csv(wflow_dis_no_bankfull)
df.set_index('time', inplace=True)
df.columns.name = "Q_gauges_locs"
df.index = pd.to_datetime(df.index)
df.index = (df.index - reftime_object).total_seconds()
df = df[q_locs['index'].astype(str).values] # adjust gauge order to match the gauge location in the sfincs.src file
df.to_csv(
    join(sfincs_model_folder, "sfincs.dis"),
    sep=" ",
    header=False,
)

#%%
config.update({"disfile": "sfincs.dis"})
config.update({"srcfile": "sfincs.src"})
inp = SfincsInput.from_dict(config)
inp.write(inp_fn=join(sfincs_model_folder, "sfincs.inp"))

# %%

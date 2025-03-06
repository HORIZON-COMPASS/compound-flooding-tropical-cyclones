# %%
from datetime import datetime as datetime
from os.path import join
from hydromt.log import setuplog
from hydromt_wflow import WflowModel
from hydromt_sfincs.sfincs_input import SfincsInput

# %%
# model and data paths/
logger = setuplog("update", "./hydromt.log", log_level=10)
if "snakemake" in locals():
    sfincs_model_folder = snakemake.params.dir_run_with_forcing
    wflow_root = snakemake.params.wflow_root_forcing
    data_cats = snakemake.params.data_cats
else:
    curdir = '../../../'
    region = "sofala"
    tc_name = "Idai"
    precip_forcing = "era5_hourly"
    wind_forcing = 'spw_IBTrACS'
    tidemodel = 'GTSMv41opendap' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap, tpxo80_opendap
    CF_SLR_txt = "0"
    CF_wind_txt = "0"
    CF_rain_txt = "0"
    wflow_root = f"p:/11210471-001-compass/03_Runs/{region}/{tc_name}/wflow/event_precip_{precip_forcing}"
    sfincs_model_folder = f"p:/11210471-001-compass/03_Runs/{region}/{tc_name}/sfincs/event_tp_{precip_forcing}_CF{CF_rain_txt}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}"
    data_cats = [
            join(curdir, "03_data_catalogs", "datacatalog_general.yml"), 
            join(curdir, "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling.yml"), 
            join(curdir, "03_data_catalogs", "datacatalog_SFINCS_obspoints.yml"),
            join(curdir, "03_data_catalogs", "datacatalog_CF_forcing.yml")
        ]


#%%
mod = WflowModel(
    root=join(wflow_root, "events"),
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

config.update({"disfile": "sfincs.dis"})
config.update({"srcfile": "sfincs.src"})
inp = SfincsInput.from_dict(config)
inp.write(inp_fn=join(sfincs_model_folder, "sfincs.inp"))

# %%

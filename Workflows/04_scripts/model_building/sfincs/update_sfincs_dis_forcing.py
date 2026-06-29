# %%
# Couple Wflow event discharge into the SFINCS model as sfincs.dis — HydroMT v1.
# Migrated from the v0 API: hydromt_sfincs.sfincs_input.SfincsInput is gone in v2 (the sfincs.inp
# is the `config` component now); WflowModel -> WflowSbmModel; results['netcdf'] -> output_scalar.
from datetime import datetime as datetime
from os.path import join
import pandas as pd
from hydromt_sfincs import SfincsModel
from hydromt_wflow import WflowSbmModel

# %%
if "snakemake" in locals():
    sfincs_model_folder   = snakemake.params.dir_run_with_forcing
    wflow_root            = snakemake.params.wflow_root_forcing
    data_cats             = snakemake.params.data_cats
    wflow_base            = snakemake.params.wflow_base
else:
    curdir              = '../../../'
    region              = "sofala"
    tc_name             = "Idai"
    precip_forcing      = "era5_hourly_zarr"
    wind_forcing        = 'spw_IBTrACS'
    tidemodel           = 'GTSMv41opendap'
    CF_rain_txt         = "0"
    CF_SLR_txt          = "0"
    CF_wind_txt         = "0"
    wflow_root          = f"/p/11210471-001-compass/03_Runs/{region}/{tc_name}/wflow/event_precip_{precip_forcing}_CF{CF_rain_txt}"
    wflow_base          = f"/p/11210471-001-compass/02_Models/{region}/{tc_name}/wflow"
    sfincs_model_folder = f"/p/11210471-001-compass/03_Runs/{region}/{tc_name}/sfincs/event_tp_{precip_forcing}_CF{CF_rain_txt}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}_nobankfull"
    data_cats           = [
        join(curdir, "03_data_catalogs", "datacatalog_general_v1___linux.yml"),
        join(curdir, "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling_v1___linux.yml"),
        join(curdir, "03_data_catalogs", "datacatalog_SFINCS_obspoints_v1___linux.yml"),
        join(curdir, "03_data_catalogs", "datacatalog_CF_forcing_v1___linux.yml"),
    ]

# %%
# Read the SFINCS config (sfincs.inp) via the v2 config component
sf = SfincsModel(root=sfincs_model_folder, mode="r+")
sf.config.read()
reftime_object = sf.config.get("tref")   # v2: config.data is not subscriptable; use .get()
if not isinstance(reftime_object, datetime):
    reftime_object = datetime.strptime(str(reftime_object), "%Y%m%d %H%M%S")

# %%
# Read the original wflow gauge order so the .dis columns match the sfincs source order
mod_ini = WflowSbmModel(root=wflow_base, mode="r", config_filename="wflow_sbm.toml")
mod_ini.geoms.read()
q_locs = mod_ini.geoms.data["gauges_locs"]   # gauge geom from setup_gauges(basename="locs")

# Read the wflow event discharge output (v1: results['netcdf'] -> output_scalar)
mod = WflowSbmModel(root=join(wflow_root, 'events'), data_libs=data_cats, mode="r")
mod.read()
df = mod.output_scalar.data['Q'].to_pandas()
df.index = (df.index - reftime_object).total_seconds()

# Order columns to match the sfincs source points, write the .dis file
df = df[q_locs['fid'].astype(str).values]   # v1 gauge id column is 'fid' (was 'index'); matches output_scalar Q columns
df.to_csv(join(sfincs_model_folder, "sfincs.dis"), sep=" ", header=False)

# %%
# Point the SFINCS config at the discharge + source files and write it back
sf.config.set("disfile", "sfincs.dis")
sf.config.set("srcfile", "sfincs.src")
sf.config.write()
# %%

# %%
from datetime import datetime as datetime

from hydromt.config import configread
from hydromt.log import setuplog
from hydromt_sfincs import SfincsModel

# %%
# model and data paths/
logger = setuplog("update", "./hydromt.log", log_level=10)
if "snakemake" in locals():
    sfincs_mod_no_forcing = snakemake.params.dir_run_no_forcing
    sfincs_mod_with_forcing = snakemake.params.dir_run_with_forcing
    forcing_yml = snakemake.params.forcing_yml
    data_cats = snakemake.params.data_cats
    spw_file = snakemake.input.spw_file_in



opt = configread(forcing_yml, abs_path=True)  # read settings from ini file
kwargs = opt.pop("global", {})

#%%
mod = SfincsModel(
    root=sfincs_mod_no_forcing,
    data_libs=data_cats,
    mode="r",
    logger=logger,
)

opt["setup_config"]["spwfile"] = spw_file

mod.update(
    model_out = sfincs_mod_with_forcing,
    write=True,
    forceful_overwrite=True,
    opt=opt
)

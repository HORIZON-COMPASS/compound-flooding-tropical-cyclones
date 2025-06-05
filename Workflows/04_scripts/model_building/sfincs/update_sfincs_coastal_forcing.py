# %%
from datetime import datetime as datetime

from hydromt.config import configread
from hydromt.log import setuplog
from hydromt import data_catalog
import shutil
import os
from os.path import join
from hydromt_sfincs import SfincsModel

# %%
# model and data paths/
logger = setuplog("update", "./hydromt.log", log_level=10)
if "snakemake" in locals():
    sfincs_mod_no_forcing = snakemake.params.dir_run_no_forcing
    sfincs_mod_with_forcing = snakemake.params.dir_run_with_forcing
    data_cats = snakemake.params.data_cats
    wind_forcing = snakemake.params.wind_forcing
    start_time = snakemake.params.start_time
    end_time = snakemake.params.end_time
    precip_forcing = snakemake.wildcards.forcing
    use_dfm = snakemake.params.use_dfm
    coastal_ts = snakemake.params.coastal_ts
    dfm_output = snakemake.params.dfm_output
    utmzone = snakemake.params.utmzone

data_catalog = data_catalog.DataCatalog(data_cats)

#%% Create a configuration
opt = {
    'setup_config': {
        'dtout': 3600,
        'dthisout': 3600,
        'storemeteo': 1,
        'utmzone':utmzone
        }
}

#  Set starttime and endtime
opt["setup_config"]["tref"] = start_time
opt["setup_config"]["tstart"] = start_time
opt["setup_config"]["tstop"] = end_time

# Add rainfall forcing
opt['setup_precip_forcing_from_grid'] = dict(precip=precip_forcing,aggregate=False)

#Add coastal water level forcing - either from DFM or from an existing time series such as GTSM
if use_dfm:
    # Add coastal water level forcing from Delft3D-FM model
    opt['setup_waterlevel_forcing'] = dict(geodataset=dfm_output,buffer=1000,merge=False)
else:
    # Add coastal water level forcing from an existing time series
    opt['setup_waterlevel_forcing'] = dict(geodataset=coastal_ts,buffer=1000,merge=False)

#%%
mod = SfincsModel(
    root=sfincs_mod_no_forcing,
    data_libs=data_cats,
    mode="r",
    logger=logger,
)

wind_forcing_str = str(wind_forcing).lower() if wind_forcing is not None else "none"

SKIP_WIND_KEYWORDS = ["no_wind", "none", "false", ""]

if wind_forcing_str not in SKIP_WIND_KEYWORDS:
    logger.info(f"Adding wind forcing using: {wind_forcing}")
    if 'spw' in wind_forcing_str: # Check if it's a spiderweb file type
        logger.info(f"Setting up SPIDERWEB wind forcing for: {wind_forcing}")
        spw_input = data_catalog[wind_forcing].path
        spw_file = os.path.basename(spw_input)
        spw_copy = os.path.join(sfincs_mod_with_forcing,spw_file)
        shutil.copyfile(spw_input, spw_copy)
        opt["setup_config"]["spwfile"] =  os.path.basename(spw_file)

        logger.info(f"Set SFINCS config: spwfile='{os.path.basename(spw_file)}', meteotype='spiderweb'")

    else: # Assuming gridded data like ERA5
        logger.info(f"Setting up gridded wind forcing using data catalog entry: {wind_forcing}")
        opt["setup_wind_forcing_from_grid"] = dict(wind=wind_forcing)
else:
    logger.info(f"Skipping wind forcing based on configuration value: '{wind_forcing}'")

mod.update(
    model_out = sfincs_mod_with_forcing,
    write=True,
    forceful_overwrite=True,
    opt=opt
)

# %%
from datetime import datetime as datetime
from os.path import basename, join, exists
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
    tc_name                 = snakemake.params.tc_name
    sfincs_mod_no_forcing   = snakemake.params.dir_run_no_forcing
    sfincs_mod_with_forcing = snakemake.params.dir_run_with_forcing
    data_cats               = snakemake.params.data_cats
    wind_forcing            = snakemake.wildcards.wind_forcing
    start_time              = snakemake.params.start_time
    end_time                = snakemake.params.end_time
    precip_forcing          = snakemake.wildcards.precip_forcing
    use_dfm                 = snakemake.params.use_dfm
    coastal_ts              = snakemake.params.coastal_ts
    dfm_output              = snakemake.params.dfm_output
    utmzone                 = snakemake.params.utmzone
    obs_points              = snakemake.params.sfincs_obs_points
    CF_wind_txt             = snakemake.wildcards.CF_wind
    CF_rain                 = float(snakemake.wildcards.CF_rain)
    CF_rain_txt             = snakemake.wildcards.CF_rain
else:
    region                  = "test"
    utmzone                 = '36s'
    tc_name                 = "Kenneth"
    wind_forcing            = 'no_wind'
    precip_forcing          = 'era5_hourly_zarr'
    dfm_res                 = "450"
    bathy                   = "gebco2024_MZB"
    tidemodel               = 'GTSMv41opendap' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap, tpxo80_opendap
    data_cats               = [
        '../../../03_data_catalogs/datacatalog_general___linux.yml',
        '../../../03_data_catalogs/datacatalog_SFINCS_obspoints___linux.yml',
        '../../../03_data_catalogs/datacatalog_SFINCS_coastal_coupling___linux.yml',
        '../../../03_data_catalogs/datacatalog_CF_forcing___linux.yml',
    ]    
    CF_rain                 = -8
    CF_rain_txt             = f"{CF_rain}"
    CF_SLR_txt              = "0"
    CF_wind_txt             = "0"
    start_time = "20190425 000000"                         # Start time of the SFINCS model run in format: YYYYMMDD HHMMSS           
    end_time = "20190430 000000"     
    use_dfm                 = False
    use_waves               = False
    # dfm_model               = f"event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}_test"
    # dfm_output              = f"dfm_output_{dfm_model}"
    sfincs_mod_no_forcing   = os.path.join(f"/p/11210471-001-compass/02_Models/{region}/{tc_name}/sfincs")
    sfincs_mod_with_forcing = os.path.join(f"/p/11210471-001-compass/03_Runs/{region}/{tc_name}/sfincs/event_tp_{precip_forcing}_CF{CF_rain_txt}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}")
    obs_points              = os.path.join("/p/11210471-001-compass/01_Data/sfincs_obs_points/obs_locs_Kenneth.geojson")
    coastal_ts = "gtsm_codec_reanalysis_hourly_v3"                   # DataCatalog handle for the coastal boundary conditions if D-FM is NOT used


#%%
data_cat = data_catalog.DataCatalog(data_cats)

if not exists(sfincs_mod_with_forcing):
    os.mkdir(sfincs_mod_with_forcing)
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
if CF_rain is None:
    print(f"Error: CF_rain value not found")
elif CF_rain == 0:
    opt['setup_precip_forcing_from_grid'] = dict(precip=f'{precip_forcing}', aggregate=False)
else:
    opt['setup_precip_forcing_from_grid'] = dict(precip=f'{precip_forcing}_CF{CF_rain_txt}_{tc_name}', aggregate=False)


#Add coastal water level forcing - either from DFM or from an existing time series such as GTSM
if use_dfm and use_waves:
    # Add coastal water level forcing from Delft3D-FM model
    opt['setup_waterlevel_forcing'] = dict(geodataset=f'{dfm_output}_waves',buffer=1000,merge=False)
elif use_dfm:
    # Add coastal water level forcing from Delft3D-FM model
    opt['setup_waterlevel_forcing'] = dict(geodataset=dfm_output,buffer=1000,merge=False)
else:
    # Add coastal water level forcing from an existing time series
    opt['setup_waterlevel_forcing'] = dict(geodataset=coastal_ts,buffer=1000,merge=False)

# Add observation points for timeserie output
opt['setup_observation_points'] = dict(locations=obs_points, merge=False)

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
mod.plot_forcing()

# %%
# Copy the subgrid folder from the base model to the event model for postprocessing of the results
if not os.path.exists(os.path.join(sfincs_mod_with_forcing, 'subgrid')):
    shutil.copytree(os.path.join(sfincs_mod_no_forcing, 'subgrid'), os.path.join(sfincs_mod_with_forcing, 'subgrid'))
else:
    print(f"Folder already exists: {os.path.join(sfincs_mod_with_forcing, 'subgrid')}")

# %%

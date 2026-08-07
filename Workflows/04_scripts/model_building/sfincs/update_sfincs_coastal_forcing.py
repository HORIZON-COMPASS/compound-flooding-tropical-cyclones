# %%
# Add event forcing (precip / coastal water level / wind / discharge / obs points) to a SFINCS
# base model — HydroMT v1 / hydromt_sfincs v2. Migrated from the v0 setup_*/opt API: forcing is
# now applied as a `steps` list of component methods via mod.update(steps=...).
from datetime import datetime as datetime
from os.path import basename, join, exists
import logging
import shutil
import os
import hydromt
from hydromt.data_catalog import DataCatalog
from hydromt_sfincs import SfincsModel

# v1: hydromt.log.setuplog is removed; use stdlib logging so the existing logger.info calls work
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("compass")

# %%
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
    use_waves               = snakemake.params.get('use_waves', False)
    coastal_ts              = snakemake.params.coastal_ts
    dfm_output              = snakemake.params.dfm_output
    utmzone                 = snakemake.params.utmzone
    obs_points              = snakemake.params.sfincs_obs_points
    skip_coastal_forcing    = snakemake.params.get('skip_coastal_forcing', False)
    skip_discharge_forcing  = snakemake.params.get('skip_discharge_forcing', False)
    # Optional wildcards for counterfactual scenarios (may not exist in simplified workflows)
    CF_wind_txt             = getattr(snakemake.wildcards, 'CF_wind', '0')
    CF_rain                 = float(snakemake.wildcards.CF_rain)
    CF_rain_txt             = snakemake.wildcards.CF_rain
    CF_SLR_txt              = getattr(snakemake.wildcards, 'CF_SLR', '0')
else:
    region                  = "durban"
    utmzone                 = '36s'
    tc_name                 = "Durban_April2022"
    wind_forcing            = 'no_wind'
    precip_forcing          = 'era5_hourly'
    dfm_res                 = "450"
    bathy                   = "gebco2024_MZB"
    tidemodel               = 'GTSMv41opendap'
    data_cats               = [
        '../../../03_data_catalogs/datacatalog_general_v1___linux.yml',
        '../../../03_data_catalogs/datacatalog_SFINCS_obspoints_v1___linux.yml',
        '../../../03_data_catalogs/datacatalog_SFINCS_coastal_coupling_v1___linux.yml',
        '../../../03_data_catalogs/datacatalog_CF_forcing_v1___linux.yml',
    ]
    CF_rain                 = 0
    CF_rain_txt             = f"{CF_rain}"
    CF_SLR_txt              = "0"
    CF_wind_txt             = "0"
    start_time = "20220410 000000"
    end_time = "20220413 000000"
    use_dfm                 = False
    use_waves               = False
    skip_coastal_forcing    = True
    skip_discharge_forcing  = True
    sfincs_mod_no_forcing   = f"/p/11210471-001-compass/02_Models/{region}/Durban2022/sfincs"
    sfincs_mod_with_forcing = f"/p/11210471-001-compass/03_Runs/{region}/Durban2022/sfincs/event_precip_{precip_forcing}_CF{CF_rain_txt}_{wind_forcing}"
    obs_points              = "/p/11210471-001-compass/01_Data/sfincs_obs_points/obs_locs_durban.geojson"
    coastal_ts = "gtsm_codec_reanalysis_hourly_v3"


# %%
data_cat = DataCatalog(data_cats)

if not exists(sfincs_mod_with_forcing):
    os.makedirs(sfincs_mod_with_forcing)

# %%
# Build the v1 steps list. Order matters: config first, then the forcing components.
steps = []

# Base config (replaces setup_config): time + output settings
config_data = {
    'dtout': 3600,
    'dthisout': 3600,
    'storemeteo': 1,
    'utmzone': utmzone,
    'tref': start_time,
    'tstart': start_time,
    'tstop': end_time,
}

# Rainfall forcing (setup_precip_forcing_from_grid -> precipitation.create)
if CF_rain is None:
    print("Error: CF_rain value not found")
elif CF_rain == 0:
    precip_name = precip_forcing
else:
    precip_name = f'{precip_forcing}_CF{CF_rain_txt}_{tc_name}'
steps.append({"precipitation.create": dict(precip=precip_name, aggregate=False)})

# Coastal water level forcing (setup_waterlevel_forcing -> water_level.create)
if not skip_coastal_forcing:
    if use_dfm and use_waves:
        logger.info(f"Adding coastal water level forcing from D-FM with waves: {dfm_output}_waves")
        steps.append({"water_level.create": dict(geodataset=f'{dfm_output}_waves', buffer=1000, merge=False)})
    elif use_dfm:
        logger.info(f"Adding coastal water level forcing from D-FM: {dfm_output}")
        steps.append({"water_level.create": dict(geodataset=dfm_output, buffer=1000, merge=False)})
    else:
        logger.info(f"Adding coastal water level forcing from time series: {coastal_ts}")
        steps.append({"water_level.create": dict(geodataset=coastal_ts, buffer=1000, merge=False)})
else:
    logger.info("Skipping coastal water level forcing (skip_coastal_forcing=True)")

# Observation points (setup_observation_points -> observation_points.create)
if exists(obs_points):
    try:
        import json
        with open(obs_points, 'r') as f:
            obs_data = json.load(f)
        if obs_data.get('features') and len(obs_data['features']) > 0:
            steps.append({"observation_points.create": dict(locations=obs_points, merge=False)})
            logger.info(f"Adding {len(obs_data['features'])} observation points from {obs_points}")
        else:
            logger.info(f"Skipping observation points (empty GeoJSON file: {obs_points})")
    except Exception as e:
        logger.warning(f"Could not read observation points file {obs_points}: {e}")
else:
    logger.info(f"Skipping observation points (file not found: {obs_points})")

# Wind forcing (spiderweb via config.spwfile, or gridded via wind.create)
wind_forcing_str = str(wind_forcing).lower() if wind_forcing is not None else "none"
SKIP_WIND_KEYWORDS = ["no_wind", "none", "false", ""]
if wind_forcing_str not in SKIP_WIND_KEYWORDS:
    logger.info(f"Adding wind forcing using: {wind_forcing}")
    if 'spw' in wind_forcing_str:  # spiderweb file
        logger.info(f"Setting up SPIDERWEB wind forcing for: {wind_forcing}")
        # The spiderweb handle is CONSTRUCTED from the CF wind value and the TC name;
        # `wind_forcing` only selects this branch and is not itself a catalog key
        # (no catalog defines e.g. 'era5_hourly_spw_IBTrACS'). This mirrors v0's
        # data_cat[f"spw_IBTrACS_CF{CF_wind_txt}_{tc_name}"].
        spw_key = f"spw_IBTrACS_CF{CF_wind_txt}_{tc_name}"
        # .full_uri, not .uri: .uri is the raw catalog string ("11210471-.../x.spw"),
        # whereas .full_uri applies the catalog root ("/p/") as v0's .path did.
        spw_input = data_cat.get_source(spw_key).full_uri
        spw_file = os.path.basename(spw_input)
        spw_copy = os.path.join(sfincs_mod_with_forcing, spw_file)
        shutil.copyfile(spw_input, spw_copy)
        config_data["spwfile"] = os.path.basename(spw_file)
        logger.info(f"Set SFINCS config: spwfile='{os.path.basename(spw_file)}', meteotype='spiderweb'")
    else:  # gridded wind (e.g. ERA5)
        logger.info(f"Setting up gridded wind forcing using data catalog entry: {wind_forcing}")
        steps.append({"wind.create": dict(wind=wind_forcing)})
else:
    logger.info(f"Skipping wind forcing based on configuration value: '{wind_forcing}'")

# Discharge forcing from gridded data (setup_discharge_forcing_from_grid ->
# discharge_points.create_from_grid)
if not skip_discharge_forcing:
    if "snakemake" in locals():
        discharge_forcing = snakemake.params.get('discharge_forcing', None)
        discharge_uparea = snakemake.params.get('discharge_uparea', 'glofas_uparea')
    else:
        discharge_forcing = 'glofas_v4_durban_apr2022'
        discharge_uparea = 'glofas_uparea'
    if discharge_forcing:
        logger.info(f"Adding discharge forcing from gridded data: {discharge_forcing}")
        steps.append({"discharge_points.create_from_grid": dict(
            discharge=discharge_forcing,
            uparea=discharge_uparea,
            wdw=1,
            rel_error=0.1,
            abs_error=100,
        )})
    else:
        logger.info("No discharge_forcing specified, skipping discharge forcing")
else:
    logger.info("Skipping discharge forcing (skip_discharge_forcing=True)")

# config.update goes first in the steps list
steps.insert(0, {"config.update": config_data})

# %%
# Read the base model and apply the forcing steps, writing to the event run directory
mod = SfincsModel(root=sfincs_mod_no_forcing, data_libs=data_cats, mode="r")
mod.update(
    model_out=sfincs_mod_with_forcing,
    write=True,
    forceful_overwrite=True,
    steps=steps,
)
mod.plot_forcing()

# %%
# Copy the subgrid folder from the base model to the event model for result postprocessing
if not os.path.exists(os.path.join(sfincs_mod_with_forcing, 'subgrid')):
    shutil.copytree(os.path.join(sfincs_mod_no_forcing, 'subgrid'),
                    os.path.join(sfincs_mod_with_forcing, 'subgrid'))
else:
    print(f"Folder already exists: {os.path.join(sfincs_mod_with_forcing, 'subgrid')}")

# %%

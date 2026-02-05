### Simplified Snakemake workflow for SFINCS with precipitation-only forcing
### No Wflow coupling, no coastal forcing, no discharge forcing
### Designed for rapid precipitation-driven flood modeling

import os
import re
from snakemake.io import Wildcards
from snakemake import shell
from os.path import join
from itertools import product

curdir = os.getcwd()

# Include the SFINCS build rules (reuse existing base model setup)
include: 'snakefile_sfincs_build.smk'

# Platform-specific path handling
if os.name == 'nt':  # Running on Windows
    root_dir = join("p:/", config['root_dir'])
elif os.name == "posix":  # Running on Linux
    root_dir = join("/p", config['root_dir'])

# Define directories
dir_runs = config['dir_runs']
dir_data = config["dir_data"]
dir_models = config["dir_models"]

# Extract configuration for wildcard combinations
runname_ids = list(config['runname_ids'].keys())
region = [value['region'] for key, value in config['runname_ids'].items()]
precip_forcing = [value['precip_forcing'] for key, value in config['runname_ids'].items()]
wind_forcing = [value['wind_forcing'] for key, value in config['runname_ids'].items()]
CF_rain = [value['CF_value_rain'] for key, value in config['runname_ids'].items()]

# Wildcard constraints
wildcard_constraints:
    precip_forcing='|'.join([re.escape(x) for x in precip_forcing]),
    wind_forcing='|'.join([re.escape(x) for x in wind_forcing]),
    CF_rain=r"-?\d*\.?\d+",

# Create run combinations (simplified - only precipitation CF, no tidemodel/SLR/wind CF)
run_combinations = []
for key, value in config['runname_ids'].items():
    for tp in value['CF_value_rain']:
        run_combinations.append((
            value['region'],
            key,
            value['precip_forcing'],
            value['wind_forcing'],
            tp
        ))

# Unpack into separate wildcard lists
region, runname_ids, precip_forcing, wind_forcing, CF_rain = zip(*run_combinations)

# ==================== Helper Functions ====================

def get_bbox(wildcards):
    """Get bounding box for SFINCS model"""
    return config["runname_ids"][wildcards.runname]["bbox_sfincs"]

def get_starttime(wildcards):
    """Get simulation start time"""
    return config["runname_ids"][wildcards.runname]["start_time"]

def get_endtime(wildcards):
    """Get simulation end time"""
    return config["runname_ids"][wildcards.runname]["end_time"]

def get_tcname(wildcards):
    """Get event name"""
    return config["runname_ids"][wildcards.runname]['tc_name']

def get_utmzone(wildcards):
    """Get UTM zone for projection"""
    return config['runname_ids'][wildcards.runname]['utmzone']

def get_datacatalog(wildcards):
    """Get data catalog files based on platform"""
    if os.name == 'nt':  # Windows
        return [
            join(curdir, '..', "03_data_catalogs", "datacatalog_general.yml"),
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling.yml"),
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_obspoints.yml"),
            join(curdir, '..', "03_data_catalogs", "datacatalog_CF_forcing.yml")
        ]
    elif os.name == "posix":  # Linux
        return [
            join(curdir, '..', "03_data_catalogs", "datacatalog_general___linux.yml"),
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling___linux.yml"),
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_obspoints___linux.yml"),
            join(curdir, '..', "03_data_catalogs", "datacatalog_CF_forcing___linux.yml")
        ]

def get_skip_coastal_forcing(wildcards):
    """Check if coastal forcing should be skipped"""
    return config['runname_ids'][wildcards.runname].get('skip_coastal_forcing', False)

def get_skip_discharge_forcing(wildcards):
    """Check if discharge forcing should be skipped"""
    return config['runname_ids'][wildcards.runname].get('skip_discharge_forcing', False)

def get_wind_forcing(wildcards):
    """Get wind forcing configuration"""
    return config['runname_ids'][wildcards.runname]['wind_forcing']

def get_use_dfm(wildcards):
    """Check if D-FM is used"""
    return config['runname_ids'][wildcards.runname].get('use_dfm', False)

def get_coastal_ts(wildcards):
    """Get coastal time series dataset"""
    return config['runname_ids'][wildcards.runname].get('coastal_ts', 'gtsm_codec_reanalysis_hourly_v3')

# ==================== Main Workflow Rules ====================

# Main rule - final output target
rule all_sfincs_precip_only:
    input:
        expand(
            join(root_dir, dir_runs, "{region}", "{runname}", "sfincs",
                 "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}",
                 "plot_output", "sfincs_basemap.png"),
            zip,
            region=region,
            runname=runname_ids,
            precip_forcing=precip_forcing,
            CF_rain=CF_rain,
            wind_forcing=wind_forcing
        )

# Rule 1: Add precipitation forcing (skip coastal and discharge)
rule add_precip_forcing_sfincs:
    input:
        msk_file = join(root_dir, dir_models, "{region}", "{runname}", "sfincs", "sfincs.msk"),
    params:
        tc_name = get_tcname,
        dir_run_no_forcing = directory(join(root_dir, dir_models, "{region}", "{runname}", "sfincs")),
        dir_run_with_forcing = lambda wildcards: directory(
            join(root_dir, dir_runs, wildcards.region, wildcards.runname, "sfincs",
                 f"event_precip_{wildcards.precip_forcing}_CF{wildcards.CF_rain}_{wildcards.wind_forcing}")
        ),
        data_cats = get_datacatalog,
        wind_forcing = get_wind_forcing,
        start_time = get_starttime,
        end_time = get_endtime,
        use_dfm = get_use_dfm,
        coastal_ts = get_coastal_ts,
        dfm_output = lambda wildcards: f"dfm_output_dummy",  # Placeholder, not used
        utmzone = get_utmzone,
        sfincs_obs_points = lambda wildcards: join(
            root_dir, dir_data, "sfincs_obs_points",
            config["runname_ids"][wildcards.runname]["sfincs_obs_file"]
        ),
        skip_coastal_forcing = get_skip_coastal_forcing,
        skip_discharge_forcing = get_skip_discharge_forcing,
    output:
        # Output is sfincs.inp (not bzs since no coastal forcing)
        inp_file = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs",
                       "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}",
                       "sfincs.inp"),
    script:
        join('..', "04_scripts", "model_building", "sfincs", "update_sfincs_coastal_forcing.py")

# Rule 2: Run SFINCS model
rule run_sfincs_model:
    threads: 16
    input:
        inp_file = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs",
                       "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}",
                       "sfincs.inp"),
    params:
        dir_run_with_forcing = lambda wildcards: directory(
            join(root_dir, dir_runs, wildcards.region, wildcards.runname, "sfincs",
                 f"event_precip_{wildcards.precip_forcing}_CF{wildcards.CF_rain}_{wildcards.wind_forcing}")
        ),
        exe = join(root_dir, dir_models, "00_executables", "SFINCS_v2.1.1_Dollerup_release_exe", 'sfincs.exe'),
    output:
        mapout = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs",
                     "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}",
                     "sfincs_map.nc"),
    run:
        if os.name == 'nt':
            import subprocess
            print(f"Running SFINCS model at {params.dir_run_with_forcing}")
            with open(join(params.dir_run_with_forcing, "sfincs.log"), "w") as f:
                subprocess.run([str(params.exe)], stdout=f, cwd=params.dir_run_with_forcing)
        if os.name == 'posix':
            shell("docker image ls")
            shell("docker run --mount src={params.dir_run_with_forcing},target=/data,type=bind deltares/sfincs-cpu:latest sfincs")

# Rule 3: Post-process and visualize flood maps
rule sfincs_plot_floodmap:
    input:
        mapout = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs",
                     "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}",
                     "sfincs_map.nc"),
    params:
        dir_run = directory(
            join(root_dir, dir_runs, "{region}", "{runname}", "sfincs",
                 "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}")
        ),
        dir_model_no_forcing = directory(join(root_dir, dir_models, "{region}", "{runname}", "sfincs")),
        datacat = get_datacatalog
    output:
        figure = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs",
                     "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}",
                     "plot_output", "sfincs_basemap.png"),
        floodmap = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs",
                       "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}",
                       "plot_output", "floodmap.tif")
    script:
        join(curdir, '..', "04_scripts", "postprocessing", "sfincs", "sfincs_postprocess.py")

### Snakemake workflow for FIAT flood impact assessment
### Compatible with the precip-only SFINCS workflow (snakefile_sfincs_precip_only.smk)
### Expects SFINCS outputs at: event_precip_{precip}_CF{rain}_{wind}/

import os
import re
from snakemake.io import Wildcards
from snakemake import shell
from os.path import join
from itertools import product

curdir = os.getcwd()
if os.name == 'nt':
    root_dir = join("p:/", config['root_dir'])
elif os.name == "posix":
    root_dir = join("/p", config['root_dir'])

dir_data   = config["dir_data"]
dir_models = config["dir_models"]
dir_runs   = config["dir_runs"]


def get_country(wildcards):
    return config['runname_ids'][wildcards.runname]['country']

def get_continent(wildcards):
    return config['runname_ids'][wildcards.runname]['continent']

def get_config(wildcards):
    config_fiat_base = config["runname_ids"][wildcards.runname]['config_fiat_base']
    return join(curdir, '..', "05_config_models", "03_fiat", config_fiat_base)

def get_datacatalog(wildcards):
    if os.name == 'nt':
        return join('..', "03_data_catalogs", "datacatalog_fiat.yml")
    elif os.name == "posix":
        return join('..', "03_data_catalogs", "datacatalog_fiat___linux.yml")


# Build run combinations from config: one entry per (runname, CF_rain) pair
run_combinations = []
for key, value in config['runname_ids'].items():
    for tp in value['CF_value_rain']:
        run_combinations.append((
            value['region'],
            key,
            value['precip_forcing'],
            value['wind_forcing'],
            tp,
        ))

region, runname_ids, precip_forcing, wind_forcing, CF_rain = zip(*run_combinations)

wildcard_constraints:
    precip_forcing='|'.join([re.escape(x) for x in set(precip_forcing)]),
    wind_forcing='|'.join([re.escape(x) for x in set(wind_forcing)]),
    CF_rain=r"-?\d*\.?\d+",


rule all_fiat_precip_only:
    input:
        expand(
            join(root_dir, dir_runs, "{region}", "{runname}", "fiat",
                 "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}",
                 "output", "spatial.fgb"),
            zip,
            region=region,
            runname=runname_ids,
            precip_forcing=precip_forcing,
            CF_rain=CF_rain,
            wind_forcing=wind_forcing,
        )


rule build_fiat_model_precip_only:
    input:
        floodmap = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs",
                        "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}",
                        "plot_output", "floodmap.tif"),
    params:
        dir_run_with_forcing = directory(
            join(root_dir, dir_runs, "{region}", "{runname}", "sfincs",
                 "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}")
        ),
        datacat_fiat  = get_datacatalog,
        model_folder  = join(root_dir, dir_runs, "{region}", "{runname}", "fiat",
                             "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}"),
        continent     = get_continent,
        country       = get_country,
        config        = get_config,
    output:
        fiat_settings = join(root_dir, dir_runs, "{region}", "{runname}", "fiat",
                             "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}",
                             "settings.toml"),
    script:
        join('..', "04_scripts", "model_building", "fiat", "setup_fiat.py")


rule run_fiat_model_precip_only:
    input:
        fiat_settings = join(root_dir, dir_runs, "{region}", "{runname}", "fiat",
                             "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}",
                             "settings.toml"),
    params:
        dir_run_with_forcing = lambda wildcards: directory(
            join(root_dir, dir_runs, wildcards.region, wildcards.runname, "fiat",
                 f"event_precip_{wildcards.precip_forcing}_CF{wildcards.CF_rain}_{wildcards.wind_forcing}")
        ),
    output:
        out = join(root_dir, dir_runs, "{region}", "{runname}", "fiat",
                   "event_precip_{precip_forcing}_CF{CF_rain}_{wind_forcing}",
                   "output", "spatial.fgb"),
    run:
        if os.name == 'nt':
            import subprocess
            print(f"Running FIAT model with settings: {input.fiat_settings}")
            with open(join(params.dir_run_with_forcing, "fiat_log.txt"), "w") as log_file:
                subprocess.run(["fiat", "run", input.fiat_settings], stdout=log_file, stderr=subprocess.PIPE, cwd=params.dir_run_with_forcing)
            print("Finished running FIAT on Windows.")
        if os.name == 'posix':
            print(f"Running FIAT model with settings: {input.fiat_settings}")
            shell(f"fiat run {input.fiat_settings}")
            print("Finished running FIAT on Linux.")

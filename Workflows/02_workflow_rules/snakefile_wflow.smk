workflow.global_resources["io_heavy"] = 1  # Only allow 1 job at a time

#%%### Import some useful python libraries
import os
from snakemake.io import Wildcards
from os.path import join
from itertools import product

curdir = os.getcwd()
if os.name == 'nt': #Running on windows
    root_dir = join("p:/",config['root_dir'])
elif os.name == "posix": #Running on linux
    root_dir = join("/p", config['root_dir'])
dir_runs = config['dir_runs']
dir_models = config['dir_models']

def get_region(wildcards):
    print(test)
    return config["runname_ids"][wildcards.runname]['region']

def get_tcname(wildcards):
    return config["runname_ids"][wildcards.runname]['tc_name']

def get_start_time(wildcards):
    return config["runname_ids"][wildcards.runname]['start_time']

def get_end_time(wildcards):
    return config["runname_ids"][wildcards.runname]['end_time']

def get_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return [
            join(curdir, '..', "03_data_catalogs", "datacatalog_general.yml"), 
            join(curdir, '..', "03_data_catalogs", "datacatalog_CF_forcing.yml")
        ]
    elif os.name == "posix": #Running on linux
        return [
            join(curdir, '..', "03_data_catalogs", "datacatalog_general___linux.yml"),
            join(curdir, '..', "03_data_catalogs", "datacatalog_CF_forcing___linux.yml")
        ]

runname_ids = list(config['runname_ids'].keys())
region = [value['region'] for key, value in config['runname_ids'].items()]
precip_forcing = [value['precip_forcing'] for key, value in config['runname_ids'].items()]
CF_rain = [value['CF_value_rain'] for key, value in config['runname_ids'].items()]
CF_landuse = [value['CF_landuse'] for key, value in config['runname_ids'].items()]

# To prevent unwanted wildcard underscore splitting
wildcard_constraints:
    precip_forcing='|'.join([re.escape(x) for x in precip_forcing]),
    CF_rain=r"-?\d*\.?\d+", # Matches integer and floating-point numbers (positive and negative)

run_combinations = []
for key, value in config['runname_ids'].items():
    for tp, lulc in product(value['CF_value_rain'], value['CF_landuse']):
        run_combinations.append((value['region'], key, value['precip_forcing'], tp, lulc))

# Unpack into separate wildcard lists
region, runname_ids, precip_forcing, CF_rain, CF_landuse = zip(*run_combinations)

# Uncommend the first "expand" line, and commend the other two, when running the first rule 'make_base_model_wflow' only, which is necessary before runnign the snakefile_wflow_30yr.smk once.
rule all_wflow:
    input:
        # expand(join(root_dir, dir_models, "{region}", "{runname}", "wflow", 'staticmaps.nc'), zip, region=region, runname=runname_ids, precip_forcing=precip_forcing, CF_rain=CF_rain),
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "wflow", "event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "events", "run_default", "output_scalar.nc"), zip, region=region, runname=runname_ids, precip_forcing=precip_forcing, CF_rain=CF_rain, landuse=CF_landuse),
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "events", "run_default", "wflow_dis_no_qbankfull.csv"), zip, region=region, runname=runname_ids, precip_forcing=precip_forcing, CF_rain=CF_rain, landuse=CF_landuse),

rule make_base_model_wflow:
    input:
        config_file = join(curdir,'..', "05_config_models", "01_wflow", "config_wflow_{landuse}.yml"),
        region_geom = join(root_dir, dir_models, "{region}", "{runname}", "sfincs_{landuse}", "gis", "region.geojson"),
        dir_sfincs_model = join(root_dir, dir_models, "{region}", "{runname}", "sfincs_{landuse}"),
        src_file = join(root_dir, dir_models, "{region}", "{runname}", "sfincs_{landuse}", "gis", "src.geojson")
    params:
        dir_model = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{landuse}"),
        data_cat = get_datacatalog,
    output: 
        toml_file = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{landuse}", 'wflow_sbm.toml'),
        staticmaps = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{landuse}", 'staticmaps.nc'), 
    script:
        join(curdir, '..', "04_scripts", "model_building", "wflow", "setup_wflow_base.py")

# update wflow forcing for warmup
rule update_forcing_wflow_warmup:
    input: 
        toml_file = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{landuse}", 'wflow_sbm.toml'),
        staticmaps = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{landuse}", 'staticmaps.nc'), 
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow", "event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "warmup", "inmaps.nc"),
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow", "event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "warmup", "wflow_sbm.toml"),
    params:
        wflow_root_noforcing = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{landuse}"),
        wflow_root_forcing= join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}"),
        start_time = get_start_time,
        end_time = get_end_time,
        data_cat = get_datacatalog,
    script:
        join(curdir, '..', "04_scripts", "model_building", "wflow", "update_forcing_wflow_warmup.py")

# update wflow forcing for event
rule update_forcing_wflow_event:
    input: 
        toml_file = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{landuse}", 'wflow_sbm.toml'),
        staticmaps = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{landuse}", 'staticmaps.nc'), 
        previous_rule = join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "warmup", "inmaps.nc")
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "events", "inmaps.nc"),
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "events", "wflow_sbm.toml"),
    params:
        wflow_root_noforcing = directory(join(root_dir, dir_models, "{region}", "{runname}", "wflow_{landuse}")),
        wflow_root_forcing= directory(join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}")),
        start_time = get_start_time,
        end_time = get_end_time,
        forcing = "{precip_forcing}",
        data_cat = get_datacatalog,
        tc_name = get_tcname
    script:
        join(curdir, '..',  "04_scripts", "model_building", "wflow", "update_forcing_wflow_event.py")

rule run_wflow_warmup:
    threads: 16
    input:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "warmup", "inmaps.nc"),
        toml = join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "warmup", "wflow_sbm.toml"),
        previous_rule = join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "events", "inmaps.nc"),  
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "events", "instate", "instates.nc"),
    params:
        exe = join(root_dir, dir_models, "00_executables", "wflow0.8.1", "wflow_cli", "bin", "wflow_cli.exe"),
        julia_env_fn = "~/.julia/environments/v1.9"
    shell:
        """
        {params.exe} {input.toml} || julia +1.9 --threads 4 --project={params.julia_env_fn} -e "using Wflow; Wflow.run()" "{input.toml}"
        """

rule run_wflow_event:
    threads: 16
    input:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "events", "instate", "instates.nc"),
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "events", "inmaps.nc"),
        toml = join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "events", "wflow_sbm.toml"),
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "events", "run_default", "output_scalar.nc"),
    params:
        exe = join(root_dir, dir_models, "00_executables", "wflow0.8.1", "wflow_cli", "bin", "wflow_cli.exe"),
        julia_env_fn = "~/.julia/environments/v1.9",
    shell:
        """
        {params.exe} {input.toml} || julia +1.9 --threads 4 --project={params.julia_env_fn} -e "using Wflow; Wflow.run()" "{input.toml}"
        """

# remove bankfull discharge 
rule postprocess_discharge:
    input:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "events", "run_default", "output_scalar.nc"),
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF0_{landuse}_30yr", "warmup", "run_default", "output_scalar.nc"),
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}", "events", "run_default", "wflow_dis_no_qbankfull.csv")
    params:
        wflow_root_forcing_30yr = directory(join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF0_{landuse}_30yr")),
        wflow_root_forcing = directory(join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}_{landuse}")),
        data_cat = get_datacatalog,
    script: join(curdir, '..',  "04_scripts", "postprocessing", "wflow", "calculate_and_remove_qbankfull.py")

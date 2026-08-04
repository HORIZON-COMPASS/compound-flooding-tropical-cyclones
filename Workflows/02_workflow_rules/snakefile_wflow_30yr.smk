workflow.global_resources["io_heavy"] = 1  # Only allow 1 job at a time

# %%
### Import some useful python libraries
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
        ]
    elif os.name == "posix": #Running on linux
        return [
            join(curdir, '..', "03_data_catalogs", "datacatalog_general___linux.yml"),
        ]

runname_ids = list(config['runname_ids'].keys())
region = [value['region'] for key, value in config['runname_ids'].items()]
precip_forcing = [value['precip_forcing'] for key, value in config['runname_ids'].items()]
landuse = [value['bankfull_corr_lulc'] for key, value in config['runname_ids'].items()]

# To prevent unwanted wildcard underscore splitting
wildcard_constraints:
    precip_forcing='|'.join([re.escape(x) for x in precip_forcing]),

# If one would want to analyse the effect of landuse on the bankfull discharge over 30 years, one can run a warmup run of wflow for 30 years for every landuse configuration (not implemented here)
run_combinations = []
for key, value in config['runname_ids'].items():
    # bankfull_corr_lulc is a single land-use handle (not a list), so there is exactly one
    # 30-yr warmup combination per run.
    lulc = value['bankfull_corr_lulc']
    run_combinations.append((value['region'], key, value['precip_forcing'], lulc))

# Unpack into separate wildcard lists
region, runname_ids, precip_forcing, landuse = zip(*run_combinations)

rule all_wflow:
    input:
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{landuse}","event_precip_{precip_forcing}_CF0_30yr", "warmup", "run_default", "output_scalar.nc"), zip, region=region, runname=runname_ids, precip_forcing=precip_forcing, landuse=landuse),

# update wflow forcing for warmup
rule update_forcing_wflow_warmup:
    input: 
        toml_file = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{landuse}", 'wflow_sbm.toml'),
        staticmaps = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{landuse}", 'staticmaps.nc'), 
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{landuse}","event_precip_{precip_forcing}_CF0_30yr", "warmup", "inmaps.nc"),
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{landuse}","event_precip_{precip_forcing}_CF0_30yr", "warmup", "wflow_sbm.toml"),
    params:
        wflow_root_noforcing = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{landuse}"),
        wflow_root_forcing= join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{landuse}","event_precip_{precip_forcing}_CF0_30yr"),
        start_time = get_start_time,
        end_time = get_end_time,
        data_cat = get_datacatalog,
    script:
        join(curdir, '..', "04_scripts", "model_building", "wflow", "update_forcing_wflow_warmup_30yr.py")

rule run_wflow_warmup:
    input:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{landuse}","event_precip_{precip_forcing}_CF0_30yr", "warmup", "inmaps.nc"),
        toml = join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{landuse}","event_precip_{precip_forcing}_CF0_30yr", "warmup", "wflow_sbm.toml"),
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{landuse}","event_precip_{precip_forcing}_CF0_30yr", "warmup", "run_default", "output_scalar.nc"),
    params:
        exe = join(root_dir, dir_models, "00_executables", "wflow0.8.1", "wflow_cli", "bin", "wflow_cli.exe"),
        julia_env_fn = "~/.julia/environments/v1.9"
    shell:
        """
        julia +1.9 --threads 4 --project={params.julia_env_fn} -e "using Wflow; Wflow.run()" "{input.toml}"
        """
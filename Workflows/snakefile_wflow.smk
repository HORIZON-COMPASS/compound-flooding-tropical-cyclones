#%%### Import some useful python libraries
import os
from snakemake.io import Wildcards
from os.path import join

curdir = os.getcwd()
if os.name == 'nt': #Running on windows
    root_dir = join("p:/",config['root_dir'])
elif os.name == "posix": #Running on linux
    root_dir = join("/p", config['root_dir'])
dir_runs = config['dir_runs']

def get_region(test):
    print(test)
    return config["runname_ids"][test]['region']

def get_experiment_id(wildcards):
    return config["runname_ids"][wildcards.runname]

def get_forcing(wildcards):
    return config["runname_ids"][wildcards.runname]['forcing']

def get_start_time(wildcards):
    return config["runname_ids"][wildcards.runname]['start_time']

def get_end_time(wildcards):
    return config["runname_ids"][wildcards.runname]['end_time']

def get_bbox(wildcards):
    prebbox = config["runname_ids"][wildcards.runname]["bbox"]
    arg_bbox = "{" + "'bbox': "+ prebbox + "}"
    return arg_bbox

# def get_dir_model_base(wildcards):
#     print(wildcards)
#     return join(root_dir, "02_Models", config["runname_ids"][wildcards.runname]['region'], config["runname_ids"][wildcards.runname], "wflow")

def get_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return "data_catalogs/datacatalog_general.yml"
    elif os.name == "posix": #Running on linux
        return "data_catalogs/datacatalog_general___linux.yml"

regions = [value['region'] for key, value in config['runname_ids'].items()]
runname_ids = list(config['runname_ids'].keys())  #
forcing = [value['forcing'] for key, value in config['runname_ids'].items()]


#TODO Update rules with expand functions
rule all_wflow:
    # input:
        # expand(join(root_dir, "02_Models", get_region("{runname}"), "{runname}", "wflow", 'wflow_sbm.toml'), runname=runname_ids)
    input:
        expand(join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "wflow", "events", "run_default", "output_scalar.nc"), zip, region=regions, runname=runname_ids, forcing=forcing),
        

rule make_base_model_wflow:
    input:
        config_file = join(curdir, "config_wflow", "wflow_build_{region}.yml"),
        region_geom = join(curdir, "..", "03_Runs", "sfincs_{runname}", "gis", "region.geojson"),
        dir_sfincs_model = "c:/Git_repos/COMPASS/03_Runs/sfincs_{runname}",
    params:
        dir_model = join(root_dir, "02_Models", "{region}", "{runname}", "wflow"),
        data_cat = get_datacatalog,
        arg_bbox = get_bbox,
    output: 
        toml_file = join(root_dir, "02_Models", "{region}", "{runname}", "wflow", 'wflow_sbm.toml'),
        staticmaps = join(root_dir, "02_Models", "{region}", "{runname}", "wflow", 'staticmaps.nc'), 
    script:
        join("scripts", "model_building", "wflow", "setup_wflow_base.py")



# update wflow forcing for warmup
rule update_forcing_wflow_warmup:
    input: 
        toml_file = join(root_dir, "02_Models", "{region}", "{runname}", "wflow", 'wflow_sbm.toml'),
        staticmaps = join(root_dir, "02_Models", "{region}", "{runname}", "wflow", 'staticmaps.nc'), 
    output:
        join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "wflow", "warmup", "inmaps.nc"),
        join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "wflow", "warmup", "wflow_sbm.toml"),
    params:
        wflow_root_noforcing = join(root_dir, "02_Models", "{region}", "{runname}", "wflow"),
        wflow_root_forcing= join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "wflow"),
        start_time = get_start_time,
        end_time = get_end_time,
        data_cat = get_datacatalog,
    script:
        join(curdir, "scripts", "preprocessing", "update_forcing_wflow_warmup.py")

# update wflow forcing for event
rule update_forcing_wflow_event:
    input: 
        toml_file = join(root_dir, "02_Models", "{region}", "{runname}", "wflow", 'wflow_sbm.toml'),
        staticmaps = join(root_dir, "02_Models", "{region}", "{runname}", "wflow", 'staticmaps.nc'), 
    output:
        join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "wflow", "events", "inmaps.nc"),
        join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "wflow", "events", "wflow_sbm.toml"),
    params:
        wflow_root_noforcing = join(root_dir, "02_Models", "{region}", "{runname}", "wflow"),
        wflow_root_forcing= join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "wflow"),
        start_time = get_start_time,
        end_time = get_end_time,
        forcing = "{forcing}",
        data_cat = get_datacatalog,
    script:
        join(curdir, "scripts","preprocessing", "update_forcing_wflow_event.py")

rule run_wflow_warmup:
    input:
        join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "wflow", "warmup", "inmaps.nc"),
        toml = join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "wflow", "warmup", "wflow_sbm.toml"),
    output:
        join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "wflow", "events", "instate", "instates.nc"),
    params:
        exe = join(root_dir, "02_Models", "00_executables", "wflow0.8.1", "wflow_cli", "bin", "wflow_cli.exe"),
        julia_env_fn = "~/.julia/environments/v1.9"
    shell:
        """
        julia --threads 4 --project={params.julia_env_fn} -e "using Wflow; Wflow.run()" "{input.toml}"
        """

rule run_wflow_event:
    input:
        join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "wflow", "events", "instate", "instates.nc"),
        join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "wflow", "events", "inmaps.nc"),
        toml = join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "wflow", "events", "wflow_sbm.toml"),
    output:
        join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "wflow", "events", "run_default", "output_scalar.nc"),
    params:
        exe = join(root_dir, "02_Models", "00_executables", "wflow0.8.1", "wflow_cli", "bin", "wflow_cli.exe"),
        julia_env_fn = "~/.julia/environments/v1.9",
    shell:
        """
        julia --threads 4 --project={params.julia_env_fn} -e "using Wflow; Wflow.run()" "{input.toml}"
        """
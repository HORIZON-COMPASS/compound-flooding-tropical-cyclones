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
dir_models = config['dir_models']

def get_region(wildcards):
    print(test)
    return config["runname_ids"][wildcards.runname]['region']

def get_forcing(wildcards):
    return config["runname_ids"][wildcards.runname]['forcing']

def get_start_time(wildcards):
    return config["runname_ids"][wildcards.runname]['start_time']

def get_end_time(wildcards):
    return config["runname_ids"][wildcards.runname]['end_time']

def get_bbox(wildcards):
    prebbox = config["runname_ids"][wildcards.runname]["bbox_sfincs"]
    arg_bbox = "{" + "'bbox': "+ prebbox + "}"
    return arg_bbox

def get_config(wildcards):
    config_wflow_base = config["runname_ids"][wildcards.runname]['config_wflow_base']
    return join(curdir, '..', "05_config_models", "01_wflow", config_wflow_base)

# def get_dir_model_base(wildcards):
#     print(wildcards)
#     return join(root_dir, dir_models, config["runname_ids"][wildcards.runname]['region'], config["runname_ids"][wildcards.runname], "wflow")

def get_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return "../03_data_catalogs/datacatalog_general.yml"
    elif os.name == "posix": #Running on linux
        return "../03_data_catalogs/datacatalog_general___linux.yml"

runname_ids = list(config['runname_ids'].keys())
regions = [value['region'] for key, value in config['runname_ids'].items()]
forcing = [value['forcing'] for key, value in config['runname_ids'].items()]

# To prevent unwanted wildcard underscore splitting
wildcard_constraints:
    forcing='|'.join([re.escape(x) for x in forcing]),

rule all_wflow:
    input:
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "wflow", "event_precip_{forcing}", "events", "run_default", "output_scalar.nc"), zip, region=regions, runname=runname_ids, forcing=forcing),
        

rule make_base_model_wflow:
    input:
        #config_file = join(curdir,'..', "05_config_models", "01_wflow", "config_wflow.yml"),
        region_geom = join(root_dir, dir_models, "{region}", "{runname}", "sfincs", "gis", "region.geojson"),
        dir_sfincs_model = join(root_dir, dir_models, "{region}", "{runname}", "sfincs"),
        src_file = join(root_dir, dir_models, "{region}", "{runname}", "sfincs", "gis", "src.geojson"),
        config_file = get_config  
    params:
        dir_model = join(root_dir, dir_models, "{region}", "{runname}", "wflow"),
        data_cat = get_datacatalog,
        arg_bbox = get_bbox,
    output: 
        toml_file = join(root_dir, dir_models, "{region}", "{runname}", "wflow", 'wflow_sbm.toml'),
        staticmaps = join(root_dir, dir_models, "{region}", "{runname}", "wflow", 'staticmaps.nc'), 
    script:
        join(curdir, '..', "04_scripts", "model_building", "wflow", "setup_wflow_base.py")

# update wflow forcing for warmup
rule update_forcing_wflow_warmup:
    input: 
        toml_file = join(root_dir, dir_models, "{region}", "{runname}", "wflow", 'wflow_sbm.toml'),
        staticmaps = join(root_dir, dir_models, "{region}", "{runname}", "wflow", 'staticmaps.nc'), 
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}", "warmup", "inmaps.nc"),
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}", "warmup", "wflow_sbm.toml"),
    params:
        wflow_root_noforcing = join(root_dir, dir_models, "{region}", "{runname}", "wflow"),
        wflow_root_forcing= join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}"),
        start_time = get_start_time,
        end_time = get_end_time,
        data_cat = get_datacatalog,
    script:
        join(curdir, '..', "04_scripts", "model_building", "wflow", "update_forcing_wflow_warmup.py")

# update wflow forcing for event
rule update_forcing_wflow_event:
    input: 
        toml_file = join(root_dir, dir_models, "{region}", "{runname}", "wflow", 'wflow_sbm.toml'),
        staticmaps = join(root_dir, dir_models, "{region}", "{runname}", "wflow", 'staticmaps.nc'), 
        #Creates a serial dependency with the previous rule to avoid error when running the workflow
        previous_rule = join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}", "warmup", "inmaps.nc") 
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}", "events", "inmaps.nc"),
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}", "events", "wflow_sbm.toml"),
    params:
        wflow_root_noforcing = directory(join(root_dir, dir_models, "{region}", "{runname}", "wflow")),
        wflow_root_forcing= directory(join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}")),
        start_time = get_start_time,
        end_time = get_end_time,
        forcing = "{forcing}",
        data_cat = get_datacatalog,
    script:
        join(curdir, '..',  "04_scripts", "model_building", "wflow", "update_forcing_wflow_event.py")

rule run_wflow_warmup:
    input:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}", "warmup", "inmaps.nc"),
        toml = join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}", "warmup", "wflow_sbm.toml"),
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}", "events", "instate", "instates.nc"),
    params:
        exe = join(root_dir, dir_models, "00_executables", "wflow0.8.1", "wflow_cli", "bin", "wflow_cli.exe"),
        julia_env_fn = "~/.julia/environments/v1.9"
    shell:
        """
        {params.exe} {input.toml} || julia +1.9 --threads 4 --project={params.julia_env_fn} -e "using Wflow; Wflow.run()" "{input.toml}"
        """

rule run_wflow_event:
    input:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}", "events", "instate", "instates.nc"),
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}", "events", "inmaps.nc"),
        toml = join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}", "events", "wflow_sbm.toml"),
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}", "events", "run_default", "output_scalar.nc"),
    params:
        exe = join(root_dir, dir_models, "00_executables", "wflow0.8.1", "wflow_cli", "bin", "wflow_cli.exe"),
        julia_env_fn = "~/.julia/environments/v1.9",
    shell:
        """
        {params.exe} {input.toml} || julia +1.9 --threads 4 --project={params.julia_env_fn} -e "using Wflow; Wflow.run()" "{input.toml}"
        """
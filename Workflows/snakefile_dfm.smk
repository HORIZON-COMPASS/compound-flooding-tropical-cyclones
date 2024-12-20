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

# define other directories:
dir_data   = config["dir_data"]
dir_models = config["dir_models"]
dir_runs   = config["dir_runs"]

def get_forcing(wildcards):
    return config["tc_name"][wildcards.tc_name]['forcing']

def get_start_time(wildcards):
    return config["tc_name"][wildcards.tc_name]['start_time']

def get_end_time(wildcards):
    return config["tc_name"][wildcards.tc_name]['end_time']

def get_bbox(wildcards):
    prebbox = config["tc_name"][wildcards.tc_name]["bbox"]
    arg_bbox = "{" + "'bbox': "+ prebbox + "}"
    return arg_bbox

def get_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return "data_catalogs/datacatalog_general.yml"
    elif os.name == "posix": #Running on linux
        return "data_catalogs/datacatalog_general___linux.yml"

# Define wildcards for path names
region = [value['region'] for key, value in config['tc_name'].items()]
tc_name = list(config['tc_name'].keys())  #
wind_forcing = [value['wind_forcing'] for key, value in config['tc_name'].items()]
dfm_res = [value['dfm_res'] for key, value in config['tc_name'].items()]
bathy = [value['bathy'] for key, value in config['tc_name'].items()]
tidemodel = [value['tidemodel'] for key, value in config['tc_name'].items()]

rule all_dfm:
    # input:
    #     expand(join(root_dir, "dir_runs", "{region}", "{tc_name}", "{forcing}", "dfm", "events", "run_default", "output_scalar.nc"), zip, region=region, tc_name=tc_name, forcing=wind_forcing),
    input:
        expand(join(root_dir, dir_models, "mozambique", "dfm", 'base_{dfm_res}_{bathy}_{tidemodel}'), dfm_res=dfm_res, bathy=bathy, tidemodel=tidemodel)

rule make_base_model_dfm:
    params:
        data_cat = get_datacatalog,
        dfm_bbox = get_dfm_bbox,
        output_bbox = get_bbox,
    output: 
        dir_model = join(root_dir, dir_models, "mozambique", "dfm", 'base_{dfm_res}_{bathy}_{tidemodel}'),
    script:
        join("scripts", "model_building", "dfm", "setup_dfm_base.py")

# rule make_dfm_model:
#      input:
#         config_file = join(curdir, "config_wflow", "wflow_build_{region}.yml"),
#         dir_dfm_model = join(root_dir, "02_Models","{region}", "{tc_name}", "wflow",)
#     params:
#         dir_model = join(root_dir, "02_Models", "{region}", "{tc_name}", "wflow"),
#         data_cat = get_datacatalog,
#         arg_bbox = get_bbox,
#     output: 
#         toml_file = join(root_dir, "02_Models", "{region}", "{tc_name}", "wflow", 'wflow_sbm.toml'),
#         staticmaps = join(root_dir, "02_Models", "{region}", "{tc_name}", "wflow", 'staticmaps.nc'), 
#     script:
#         join("scripts", "model_building", "wflow", "setup_wflow_base.py")

# rule run_dfm:
#     input:
#         join(root_dir, "03_Runs", "{region}", "{tc_name}", "{forcing}", "wflow", "events", "instate", "instates.nc"),
#         join(root_dir, "03_Runs", "{region}", "{tc_name}", "{forcing}", "wflow", "events", "inmaps.nc"),
#         toml = join(root_dir, "03_Runs", "{region}", "{tc_name}", "{forcing}", "wflow", "events", "wflow_sbm.toml"),
#     output:
#         join(root_dir, "03_Runs", "{region}", "{tc_name}", "{forcing}", "wflow", "events", "run_default", "output_scalar.nc"),
#     params:
#         exe = join(root_dir, "02_Models", "00_executables", "wflow0.8.1", "wflow_cli", "bin", "wflow_cli.exe"),
#         julia_env_fn = "~/.julia/environments/v1.9",
#     shell:
#         """
#         {params.exe} {input.toml} || julia --threads 4 --project={params.julia_env_fn} -e "using Wflow; Wflow.run()" "{input.toml}"
#         """
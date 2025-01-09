#%%### Import some useful python libraries
import os
from snakemake.io import Wildcards
from os.path import join

curdir = os.getcwd()
if os.name == 'nt': #Running on windows
    p_dir = join("p:/")
elif os.name == "posix": #Running on linux
    p_dir = join("/p")

root_dir = join(p_dir,config['root_dir'])
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
    bbox = config["tc_name"][wildcards.tc_name]["bbox"]
    return bbox

def get_dfm_bbox(wildcards):
    bbox = config["tc_name"][wildcards.tc_name]["bbox_dfm"]
    return bbox

def get_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return "data_catalogs/datacatalog_general.yml"
    elif os.name == "posix": #Running on linux
        return "data_catalogs/datacatalog_general___linux.yml"

# need to be adjusted for Linux
def get_obs_file(wildcards):
    obs_file = config["tc_name"][wildcards.tc_name]["dfm_obs_file"]
    return obs_file

# need to be adjusted for Linux
def get_verification_points(wildcards):
    verification_points = config["tc_name"][wildcards.tc_name]["verification_points"]
    return verification_points

# Define wildcards for path names
region = [value['region'] for key, value in config['tc_name'].items()]
dfm_res = [value['dfm_res'] for key, value in config['tc_name'].items()]
bathy = [value['bathy'] for key, value in config['tc_name'].items()]
tidemodel = [value['tidemodel'] for key, value in config['tc_name'].items()]
tc_name = list(config['tc_name'].keys())
wind_forcing = [value['wind_forcing'] for key, value in config['tc_name'].items()]

# To prevent unwanted wildcard underscore splitting
wildcard_constraints:
    wind_forcing='|'.join([re.escape(x) for x in wind_forcing]),

rule all_dfm:
    input:
        # expand(join(root_dir, dir_models, "{region}", "{tc_name}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}", "ext_file_new.ext"), region=region, tc_name=tc_name, dfm_res=dfm_res, bathy=bathy, tidemodel=tidemodel)
        expand(join(root_dir, dir_models, "{region}", "{tc_name}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "settings.mdu"), region=region, tc_name=tc_name, dfm_res=dfm_res, bathy=bathy, tidemodel=tidemodel, wind_forcing=wind_forcing)

rule make_model_dfm_base:
    params:
        data_cat = get_datacatalog,
        dfm_bbox = get_dfm_bbox,
        output_bbox = get_bbox,
    output: 
        dir_model = directory(join(root_dir, dir_models, "{region}", "{tc_name}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}")),
        ext_file_new = join(root_dir, dir_models, "{region}", "{tc_name}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}", "ext_file_new.ext"),
    script:
        join("scripts", "model_building", "dfm", "setup_dfm_base.py")

rule make_dfm_model_event:
    input:
        ext_file_new = join(root_dir, dir_models, "{region}", "{tc_name}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}", 'ext_file_new.ext'),
        dimrset      = join(root_dir, "/d-hydro/dimrset/weekly/2.25.17.78708"),
        base_mdu     = join("scripts", "model_building", "dfm", "base_model_settings.mdu"),
        batchfile_h7 = join("scripts", "model_building", "dfm", "submit_singularity_h7.sh"),
    params:
        dir_base_model = directory(join(root_dir, dir_models, "{region}", "{tc_name}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}")),
        start_time   = get_start_time,
        end_time     = get_end_time,
        dfm_bbox     = get_dfm_bbox,
        output_bbox  = get_bbox,
        dfm_obs_file = get_obs_file,
        verification_points = get_verification_points,
        data_cat     = get_datacatalog,
    output: 
        dir_event_model = directory(join(root_dir, dir_models, "{region}", "{tc_name}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}")),
        mdu_file = join(root_dir, dir_models, "{region}", "{tc_name}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "settings.mdu"),
    script:
        join("scripts", "model_building", "dfm", "setup_dfm_event.py")

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
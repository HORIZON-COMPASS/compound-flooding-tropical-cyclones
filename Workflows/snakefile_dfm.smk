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

def get_sfincs_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return "data_catalogs/datacatalog_sfincs.yml"
    elif os.name == "posix": #Running on linux
        return "data_catalogs/datacatalog_sfincs___linux.yml"

# # need to be adjusted for Linux
# def get_obs_file(wildcards):
#     return config["tc_name"][wildcards.tc_name]["dfm_obs_file"]

# # need to be adjusted for Linux
# def get_verification_points(wildcards):
#     verification_points = config["tc_name"][wildcards.tc_name]["verification_points"]
#     return verification_points

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
        expand(join(root_dir, dir_runs, "{region}", "{tc_name}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}.mdu"), region=region, tc_name=tc_name, dfm_res=dfm_res, bathy=bathy, tidemodel=tidemodel, wind_forcing=wind_forcing)
        # expand(join(root_dir, dir_runs, "{region}", "{tc_name}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "output", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}_his.nc"), region=region, tc_name=tc_name, dfm_res=dfm_res, bathy=bathy, tidemodel=tidemodel, wind_forcing=wind_forcing)
        # "data_catalogs/datacatalog_sfincs.yml"

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
    params:
        dir_base_model = directory(join(root_dir, dir_models, "{region}", "{tc_name}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}")),
        start_time   = get_start_time,
        end_time     = get_end_time,
        dfm_bbox     = get_dfm_bbox,
        output_bbox  = get_bbox,
        # dfm_obs      = get_obs_file,
        dfm_obs_file = lambda wildcards: join(root_dir, dir_data, "Coastal_boundary", "points", config["tc_name"][wildcards.tc_name]["dfm_obs_file"]),
        # verif_points = get_verification_points,
        verif_points_file = lambda wildcards: join(root_dir, dir_data, "Coastal_boundary", "points", config["tc_name"][wildcards.tc_name]["verification_points"]),
        data_cat     = get_datacatalog,
        dimrset      = join(p_dir, "d-hydro", "dimrset", "weekly", "2.28.06"),
        uniformwind  = join(root_dir, dir_data, "uniformwind0.wnd"),
        model_name   = "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}",
    output: 
        dir_event_model = directory(join(root_dir, dir_runs, "{region}", "{tc_name}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}")),
        mdu_file = join(root_dir, dir_runs, "{region}", "{tc_name}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}.mdu"),
        submit_script_linux = join(root_dir, dir_runs, "{region}", "{tc_name}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "submit_singularity_h7.sh"),
        submit_script_windows = join(root_dir, dir_runs, "{region}", "{tc_name}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "run_parallel.bat"),
    script:
        join("scripts", "model_building", "dfm", "setup_dfm_event.py")

rule run_dfm:
    input:
        submit_script_linux = join(root_dir, dir_runs, "{region}", "{tc_name}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "submit_singularity_h7.sh"),
        submit_script_windows = join(root_dir, dir_runs, "{region}", "{tc_name}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "run_parallel.bat"),
    params:
        windows_output = directory(join(root_dir, dir_runs, "{region}", "{tc_name}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "output")),
    output:
        his_file = join(root_dir, dir_runs, "{region}", "{tc_name}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "output", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}_his.nc"),        
    run:
        if os.name == 'nt':
            print("Executing DFM...")
            shell("cmd /c {input.submit_script_windows}")
            print("Finished running")
            # Move the output folder one level up on Windows
            # os.rename("{params.windows_output}", "{params.model_output}")  # Move folder one level up
        if os.name == 'posix':
            shell("sbatch {input.submit_script_linux}")

rule add_dfm_output_to_catalog:
    input:
        his_file = join(root_dir, dir_runs, "{region}", "{tc_name}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "his.nc"),
    params:
        model_name      = "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}",
        dir_event_model = directory(join(root_dir, dir_runs, "{region}", "{tc_name}", "dfm", "{model_name}")),
        sfincs_data_cat = get_sfincs_datacatalog,
    output:
        sfincs_data_cat_update = "{params.sfincs_data_cat}"
    script:
        join("scripts", "postprocessing", "dfm", "output_to_catalog.py")

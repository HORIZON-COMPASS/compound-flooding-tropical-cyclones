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
    return config["runname_ids"][wildcards.runname]['forcing']

def get_start_time(wildcards):
    return config["runname_ids"][wildcards.runname]['start_time']

def get_end_time(wildcards):
    return config["runname_ids"][wildcards.runname]['end_time']

def get_bbox(wildcards):
    bbox = config["runname_ids"][wildcards.runname]["bbox"]
    return bbox

def get_dfm_bbox(wildcards):
    bbox = config["runname_ids"][wildcards.runname]["bbox_dfm"]
    return bbox

def get_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return "data_catalogs/datacatalog_general.yml"
    elif os.name == "posix": #Running on linux
        return "data_catalogs/datacatalog_general___linux.yml"

def get_sfincs_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return "data_catalogs/datacatalog_SFINCS_coastal_coupling.yml"
    elif os.name == "posix": #Running on linux
        return "data_catalogs/datacatalog_SFINCS_coastal_coupling___linux.yml"

# Define wildcards for path names
runname_ids = list(config['runname_ids'].keys())
tc_name = [value['tc_name'] for key, value in config['runname_ids'].items()]
region = [value['region'] for key, value in config['runname_ids'].items()]
dfm_res = [value['dfm_res'] for key, value in config['runname_ids'].items()]
bathy = [value['bathy'] for key, value in config['runname_ids'].items()]
tidemodel = [value['tidemodel'] for key, value in config['runname_ids'].items()]
wind_forcing = [value['wind_forcing'] for key, value in config['runname_ids'].items()]

# To prevent unwanted wildcard underscore splitting
wildcard_constraints:
    wind_forcing='|'.join([re.escape(x) for x in wind_forcing]),

# Define the script dynamically based on OS before the rule
submit_script_system = "run_parallel.bat" if os.name == 'nt' else "submit_singularity_h7.sh"

rule all_dfm:
    input:
        # expand(join(root_dir, dir_models, "{region}", "{runname}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}", "ext_file_new.ext"), zip, region=region, runname=runname_ids, dfm_res=dfm_res, bathy=bathy, tidemodel=tidemodel)
        # expand(join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}.mdu"), zip, region=region, runname=runname_ids, tc_name=tc_name, dfm_res=dfm_res, bathy=bathy, tidemodel=tidemodel, wind_forcing=wind_forcing)
        # expand(join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "output", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}_his.nc"), zip, region=region, runname=runname_ids, tc_name=tc_name, dfm_res=dfm_res, bathy=bathy, tidemodel=tidemodel, wind_forcing=wind_forcing),
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "postprocessing_done.txt"), zip, region=region, runname=runname_ids, tc_name=tc_name, dfm_res=dfm_res, bathy=bathy, tidemodel=tidemodel, wind_forcing=wind_forcing),   

rule make_model_dfm_base:
    params:
        data_cat = get_datacatalog,
        dfm_bbox = get_dfm_bbox,
        output_bbox = get_bbox,
    output: 
        dir_model = directory(join(root_dir, dir_models, "{region}", "{runname}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}")),
        ext_file_new = join(root_dir, dir_models, "{region}", "{runname}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}", "ext_file_new.ext"),
    script:
        join("scripts", "model_building", "dfm", "setup_dfm_base.py")

rule make_dfm_model_event:
    input:
        ext_file_new = join(root_dir, dir_models, "{region}", "{runname}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}", "ext_file_new.ext"),
    params:
        dir_base_model = directory(join(root_dir, dir_models, "{region}", "{runname}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}")),
        start_time   = get_start_time,
        end_time     = get_end_time,
        dfm_bbox     = get_dfm_bbox,
        output_bbox  = get_bbox,
        dfm_obs_file = lambda wildcards: join(root_dir, dir_data, "Coastal_boundary", "points", config["runname_ids"][wildcards.runname]["dfm_obs_file"]),
        verif_points_file = lambda wildcards: join(root_dir, dir_data, "Coastal_boundary", "points", config["runname_ids"][wildcards.runname]["verification_points"]),
        data_cat     = get_datacatalog,
        dimrset      = join(p_dir, "d-hydro", "dimrset", "weekly", "2.28.06"),
        uniformwind  = join(root_dir, dir_data, "uniformwind0.wnd"),
        model_name   = "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}",
    output: 
        dir_event_model = directory(join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}")),
        mdu_file = join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}.mdu"),
        submit_script = join(root_dir,dir_runs,"{region}", "{runname}","dfm", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}",submit_script_system),
    script:
        join("scripts", "model_building", "dfm", "setup_dfm_event.py")

rule run_dfm:
    input:
        submit_script = join(root_dir,dir_runs,"{region}", "{runname}","dfm", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}",submit_script_system),
    output:
        his_file = join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "output", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}_his.nc"),        
    run:
        if os.name == 'nt':
            print("Executing DFM...")
            shell("cmd /c {input.submit_script}")
            print("Finished running")
        if os.name == 'posix':
            shell("sbatch {input.submit_script}")

rule add_dfm_output_to_catalog:
    input:
        his_file = join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "output", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}_his.nc"),
    params:
        model_name       = "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}",
        sfincs_data_cat  = get_sfincs_datacatalog,
        root_dir         = p_dir,
    output:
        done_file = join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "postprocessing_done.txt"),
    script:
        join("scripts", "postprocessing", "dfm", "output_to_catalog.py")
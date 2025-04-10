#%%### Import some useful python libraries
import os
from snakemake.io import Wildcards
from os.path import join

curdir = os.getcwd()
if os.name == 'nt': #Running on windows
    disk_dir = join("p:/")
elif os.name == "posix": #Running on linux
    disk_dir = join("/p")

root_dir = join(disk_dir,config['root_dir'])

# define other directories:
dir_models = config["dir_models"]
dir_runs   = config["dir_runs"]

def get_forcing(wildcards):
    return config["runname_ids"][wildcards.runname]['forcing']

def get_tcname(wildcards):
    return config["runname_ids"][wildcards.runname]['tc_name']

def get_start_time(wildcards):
    return config["runname_ids"][wildcards.runname]['start_time']

def get_end_time(wildcards):
    return config["runname_ids"][wildcards.runname]['end_time']

def get_sfincs_bbox(wildcards):
    bbox = config["runname_ids"][wildcards.runname]["bbox_sfincs"]
    return bbox

def get_dfm_bbox(wildcards):
    bbox = config["runname_ids"][wildcards.runname]["bbox_dfm"]
    return bbox

def get_dfm_dxy_base(wildcards):
    dfm_dxy_base = config["runname_ids"][wildcards.runname]["dfm_dxy_base"]
    return dfm_dxy_base

def get_dfm_obs_points(wildcards):
    dfm_obs_file = config["runname_ids"][wildcards.runname]["dfm_obs_file"]
    return dfm_obs_file
    
def get_dfm_verification_points(wildcards):
    verification_points = config["runname_ids"][wildcards.runname]["dfm_verification_points"]
    return verification_points

def get_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return join(curdir, '..', "03_data_catalogs", "datacatalog_general.yml")
    elif os.name == "posix": #Running on linux
        return join(curdir, '..', "03_data_catalogs", "datacatalog_general___linux.yml")

def get_sfincs_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return  join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling.yml")
    elif os.name == "posix": #Running on linux
        return join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling___linux.yml")

# Define wildcards for path names
runname_ids = list(config['runname_ids'].keys())
region = [value['region'] for key, value in config['runname_ids'].items()]
dfm_res = [value['dfm_res'] for key, value in config['runname_ids'].items()]
bathy = [value['bathy'] for key, value in config['runname_ids'].items()]
tidemodel = [value['tidemodel'] for key, value in config['runname_ids'].items()]
wind_forcing = [value['wind_forcing'] for key, value in config['runname_ids'].items()]

# To prevent unwanted wildcard underscore splitting
wildcard_constraints:
    wind_forcing='|'.join([re.escape(x) for x in wind_forcing]),
    bathy='|'.join([re.escape(x) for x in bathy]),

# Define the script dynamically based on OS before the rule
submit_script_system = "run_parallel.bat" if os.name == 'nt' else "run_singularity_h7.sh"
submit_script_system_copy = "run_parallel_copy.bat" if os.name == 'nt' else "run_singularity_h7_copy.sh"

rule all_dfm:
    input:
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "postprocessing_done.txt"), zip, region=region, runname=runname_ids, dfm_res=dfm_res, bathy=bathy, tidemodel=tidemodel, wind_forcing=wind_forcing),   
        

rule make_model_dfm_base:
    params:
        data_cat = get_datacatalog,
        dfm_bbox = get_dfm_bbox,
        output_bbox = get_sfincs_bbox,
        dfm_dxy_base = get_dfm_dxy_base,
    output: 
        dir_model = directory(join(root_dir, dir_models, "{region}", "{runname}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}")),
        ext_file_new = join(root_dir, dir_models, "{region}", "{runname}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}", "ext_file_new.ext"),
    resources:
        partition = '4vcpu',
        time = '0-0:30:00',
        jobname = 'dfm_base',
        taskspernode = 4,
    script:
        join("..", "04_scripts", "model_building", "dfm", "setup_dfm_base.py")
        
rule make_dfm_model_event:
    input:
        ext_file_new = join(root_dir, dir_models, "{region}", "{runname}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}", "ext_file_new.ext"),
    params:
        tc_name = get_tcname,
        dir_base_model = directory(join(root_dir, dir_models, "{region}", "{runname}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}")),
        start_time   = get_start_time,
        end_time     = get_end_time,
        dfm_bbox     = get_dfm_bbox,
        output_bbox  = get_sfincs_bbox,
        verif_points = get_dfm_verification_points,
        data_cat     = get_datacatalog,
        sfincs_data_cat = get_sfincs_datacatalog,
        dimrset      = join(disk_dir, "d-hydro", "dimrset", "weekly", "2.28.06"),
        model_name   = "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}",
        dfm_obs_file = get_dfm_obs_points,
        submit_script_file = join(root_dir,dir_runs,"{region}", "{runname}","dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}",submit_script_system),
    output:  
        dir_event_model = directory(join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}")),
        mdu_file = join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}.mdu"),
        submit_script = join(root_dir,dir_runs,"{region}", "{runname}","dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}",submit_script_system),
    resources:
        partition = '4vcpu',
        time = '0-0:30:00',
        jobname = 'dfm_event',
        taskspernode = 4,
    script:
        join("..", "04_scripts", "model_building", "dfm", "setup_dfm_event.py")


rule run_dfm:
    input:
        submit_script = join(root_dir,dir_runs,"{region}", "{runname}","dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}",submit_script_system),
    params:
        dir_event_model = directory(join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}")),
        submit_script_copy = join(root_dir,dir_runs,"{region}", "{runname}","dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}",submit_script_system_copy),
    output:
        his_file = join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "output", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}_0000_his.nc"),        
    resources:
        partition = '16vcpu',
        time =  '0-3:00:00',
        jobname = 'dfm_run',    
        taskspernode = 16,
    run:
        if os.name == 'nt':
            print("Executing DFM...")
            shell("cmd /c {input.submit_script}")
            print("Finished running")
        if os.name == 'posix':
            shell("cp {input.submit_script} {params.submit_script_copy}")
            shell("chmod +x {params.submit_script_copy}")
            shell("{params.submit_script_copy}")

rule add_dfm_output_to_catalog:
    input:
        his_file = join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "output", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}_0000_his.nc"),
    params:
        model_name       = "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}",
        sfincs_data_cat  = get_sfincs_datacatalog,
        root_dir         = disk_dir,
    output:
        done_file = join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}", "postprocessing_done.txt"),
    resources:
        partition = '4vcpu',
        time = '0-0:30:00',
        jobname = 'dfm_out',
        taskspernode = 4,
    script:
        join("..", "04_scripts", "postprocessing", "dfm", "output_to_catalog.py")
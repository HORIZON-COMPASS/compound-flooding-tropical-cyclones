#%%### Import some useful python libraries
import os
from snakemake.io import Wildcards
from os.path import join
from itertools import product

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

def get_tc_ext_days(wildcards):
    return config["runname_ids"][wildcards.runname]['ext_days']

def get_tc_sid(wildcards):
    return config["runname_ids"][wildcards.runname]['sid']

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

def get_sfincs_datacatalog2(wildcards):
    if os.name == 'nt': #Running on windows
        return  join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_obspoints.yml")
    elif os.name == "posix": #Running on linux
        return join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_obspoints___linux.yml")

# Define wildcards for path names
runname_ids = list(config['runname_ids'].keys())
region = [value['region'] for key, value in config['runname_ids'].items()]
dfm_res = [value['dfm_res'] for key, value in config['runname_ids'].items()]
bathy = [value['bathy'] for key, value in config['runname_ids'].items()]
tidemodel = [value['tidemodel'] for key, value in config['runname_ids'].items()]
# CF_SLR = [value['CF_value_SLR'] for key, value in config['runname_ids'].items()]
wind_forcing = [value['wind_forcing'] for key, value in config['runname_ids'].items()]
# CF_wind = [value['CF_value_wind'] for key, value in config['runname_ids'].items()]

# To prevent unwanted wildcard underscore splitting
wildcard_constraints:
    wind_forcing='|'.join([re.escape(x) for x in wind_forcing]),
    bathy='|'.join([re.escape(x) for x in bathy]),
    CF_wind=r"-?\d+"  # Ensures only numbers are captured (prevents '10_his.nc')


# activate when having multiple CF values!!
#  Generate all combinations of CF_SLR and CF_wind for each runname
run_combinations = []
for key, value in config['runname_ids'].items():
    for slr, wind in product(value['CF_value_SLR'], value['CF_value_wind']):
        run_combinations.append((value['region'], key, value['dfm_res'], value['bathy'], value['tidemodel'], slr, value['wind_forcing'], wind))

# Unpack into separate wildcard lists
region, runname_ids, dfm_res, bathy, tidemodel, CF_SLR, wind_forcing, CF_wind = zip(*run_combinations)

# Define the script dynamically based on OS before the rule
submit_script_system = "run_parallel.bat" if os.name == 'nt' else "run_singularity_h7.sh"
submit_script_system_copy = "run_parallel_copy.bat" if os.name == 'nt' else "run_singularity_h7_copy.sh"

rule all_dfm:
    input:
        # expand(join(root_dir, dir_models, "{region}", "{runname}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}", "ext_file_new.ext"), zip, region=region, runname=runname_ids, dfm_res=dfm_res, bathy=bathy, tidemodel=tidemodel, CF_SLR=CF_SLR,) 
        # expand(join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "settings.mdu"), zip, region=region, runname=runname_ids, dfm_res=dfm_res, bathy=bathy, tidemodel=tidemodel, CF_SLR=CF_SLR, wind_forcing=wind_forcing, CF_wind=CF_wind)
        # expand(join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "output", "settings_0000_his.nc"), zip, region=region, runname=runname_ids, dfm_res=dfm_res, bathy=bathy, tidemodel=tidemodel, CF_SLR=CF_SLR, wind_forcing=wind_forcing, CF_wind=CF_wind),
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "postprocessing_done.txt"), zip, region=region, runname=runname_ids, dfm_res=dfm_res, bathy=bathy, tidemodel=tidemodel, CF_SLR=CF_SLR, wind_forcing=wind_forcing, CF_wind=CF_wind),   

rule make_model_dfm_base:
    params:
        data_cat = get_datacatalog,
        dfm_bbox = get_dfm_bbox,
        output_bbox = get_sfincs_bbox,
        dfm_dxy_base = get_dfm_dxy_base,
    output: 
        dir_model = directory(join(root_dir, dir_models, "{region}", "{runname}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}")),
        ext_file_new = join(root_dir, dir_models, "{region}", "{runname}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}", "ext_file_new.ext"),
    script:
        join("..", "04_scripts", "model_building", "dfm", "setup_dfm_base.py")

# Modify rule to also create CF for other datasets than spw files 
# rule create_CF_wind_forcing:
#     params:
#         start_date = get_start_time,
#         end_date   = get_end_time,
#         tc_name    = get_tcname,
#         data_cat   = get_datacatalog,
#         root_dir   = p_dir,
#         sid        = get_tc_sid,
#         ext_days   = get_tc_ext_days,
#     output:
#         path_CF_wind = join(root_dir, dir_data, "SPW_forcing_files", "tc_{tc_name}_CF{CF_wind}_{sid}_ext{ext_days}d.spw")
#     script:
#         join("..", "04_scripts", "preprocessing", "create_tc_IBTrACS_wind.py")

rule make_dfm_model_event:
    input:
        ext_file_new = join(root_dir, dir_models, "{region}", "{runname}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}", "ext_file_new.ext"),
    params:
        tc_name = get_tcname,
        dir_base_model = directory(join(root_dir, dir_models, "{region}", "{runname}", "dfm", "base_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}")),
        dir_event_model = directory(join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}")),
        start_time   = get_start_time,
        end_time     = get_end_time,
        dfm_bbox     = get_dfm_bbox,
        output_bbox  = get_sfincs_bbox,
        verif_points = get_dfm_verification_points,
        data_cat     = get_datacatalog,
        sfincs_data_cat = get_sfincs_datacatalog,
        dimrset      = join(p_dir, "d-hydro", "dimrset", "weekly", "2.28.06"),
        uniformwind  = join(root_dir, dir_data, "uniformwind0.wnd"),
        model_name   = "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}",
        dfm_obs_file = get_dfm_obs_points,
    output: 
        mdu_file = join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "settings.mdu"),
        submit_script = join(root_dir,dir_runs,"{region}", "{runname}","dfm", "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}",submit_script_system),
    script:
        join("..", "04_scripts", "model_building", "dfm", "setup_dfm_event.py")
        join("..", "04_scripts", "model_building", "dfm", "setup_dfm_event.py")

rule run_dfm:
    input:
        mdu_file = join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "settings.mdu"),
        submit_script = join(root_dir,dir_runs,"{region}", "{runname}","dfm", "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}",submit_script_system),
    params:
        dir_event_model = directory(join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}")),
        submit_script_copy = join(root_dir,dir_runs,"{region}", "{runname}","dfm", "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}",submit_script_system_copy),    
    output:
        his_file = join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "output", "settings_0000_his.nc"),        
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
        his_file = join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "output", "settings_0000_his.nc"),
    params:
        model_name       = "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}",
        sfincs_data_cat  = get_sfincs_datacatalog,
        root_dir         = disk_dir,
    output:
        done_file = join(root_dir, dir_runs, "{region}", "{runname}", "dfm", "event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "postprocessing_done.txt"),
    script:
        join("..", "04_scripts", "postprocessing", "dfm", "output_to_catalog.py")
        join("..", "04_scripts", "postprocessing", "dfm", "output_to_catalog.py")
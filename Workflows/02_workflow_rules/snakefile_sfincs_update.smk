### Import some useful python libraries
import os
from snakemake.io import Wildcards
from snakemake import shell
from os.path import join

curdir = os.getcwd()
if os.name == 'nt': #Running on windows
    root_dir = join("p:/",config['root_dir'])
elif os.name == "posix": #Running on linux
    root_dir = join("/p", config['root_dir'])

# define other directories:
dir_data   = config["dir_data"]
dir_models = config["dir_models"]
dir_runs   = config["dir_runs"]

def get_bbox(wildcards):
    bbox = config["runname_ids"][wildcards.runname]["bbox_sfincs"]
    return bbox

def get_starttime(wildcards):
    start_time = config["runname_ids"][wildcards.runname]["start_time"]
    return start_time

def get_endtime(wildcards):
    end_time = config["runname_ids"][wildcards.runname]["end_time"]
    return end_time

def get_wind_forcing(wildcards):
    # Returns the data catalog handle for the .spw file 
    return config['runname_ids'][wildcards.runname]['wind_forcing']

def get_utmzone(wildcards):
    return config['runname_ids'][wildcards.runname]['utmzone']

def get_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return [
            join(curdir, '..', "03_data_catalogs", "data_catalog_MO.yml"), 
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling.yml"), 
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_obspoints.yml")
        ]
    elif os.name == "posix": #Running on linux
        return [
            join(curdir, '..', "03_data_catalogs", "data_catalog_MO.yml"), 
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling___linux.yml"), 
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_obspoints___linux.yml")
        ]


def get_use_dfm(wildcards):
    return config['runname_ids'][wildcards.runname]['use_dfm']

def get_coastal_ts(wildcards):
    return config['runname_ids'][wildcards.runname]['coastal_ts']

runname_ids = list(config['runname_ids'].keys())
regions = [value['region'] for key, value in config['runname_ids'].items()]
wind_forcing = [value['wind_forcing'] for key, value in config['runname_ids'].items()]
forcing = [value['forcing'] for key, value in config['runname_ids'].items()]

# To prevent unwanted wildcard underscore splitting
wildcard_constraints:
    forcing='|'.join([re.escape(x) for x in forcing]),
    wind_forcing='|'.join([re.escape(x) for x in wind_forcing]),

rule all_sfincs_update:
    input:
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "plot_output", "sfincs_basemap.png"), zip, region=regions, runname=runname_ids, forcing=forcing)


rule add_forcing_coastal_meteo_sfincs:
    input:
        msk_file = join(root_dir, dir_models, "{region}", "{runname}", "sfincs", "sfincs.msk"),
    params:
        dir_run_no_forcing = directory(join(root_dir, dir_models, "{region}", "{runname}", "sfincs")),
        dir_run_with_forcing = directory(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}")),
        data_cats = get_datacatalog,
        wind_forcing = get_wind_forcing,
        start_time = get_starttime, 
        end_time = get_endtime,
        use_dfm = get_use_dfm,
        coastal_ts = get_coastal_ts,
        dfm_output = lambda wildcards: "dfm_output_event_"+ config['runname_ids'][wildcards.runname]["dfm_res"] + "_" + config['runname_ids'][wildcards.runname]["bathy"] + "_" + config['runname_ids'][wildcards.runname]["tidemodel"] + "_" + config['runname_ids'][wildcards.runname]["wind_forcing"],
        utmzone = get_utmzone,
    output:
        bzs_file = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "sfincs.bzs"),
    script:
        join( '..', "04_scripts", "model_building", "sfincs", "update_sfincs_coastal_forcing.py")

rule update_dis_forcing_sfincs:
    input:
        bzs_file = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "sfincs.bzs"),
        wflow_output = join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}", "events", "run_default", "output_scalar.nc"),
    params:
        dir_run_with_forcing = directory(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}")),
        wflow_root_forcing = directory(join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}")),
        data_cats = get_datacatalog,
    output:
        dis_file = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "sfincs.dis"),
    script:
        join( '..', "04_scripts", "model_building", "sfincs", "update_sfincs_dis_forcing.py")


rule run_sfincs_model:
    input:
#        batchfile = "{dir_run}"+"/sfincs_"+"{runname}"+"/run_sfincs.bat"
        dis_file = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "sfincs.dis"),
    params:
        dir_run_with_forcing = directory(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}")),
        exe = join(root_dir, dir_models, "00_executables", "SFINCS_v2.1.1_Dollerup_release_exe", 'sfincs.exe'),
        currentdir = curdir
    output:
        mapout = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "sfincs_map.nc"),
    run:
        if os.name == 'nt':
            import subprocess
            print(f"Running Sfincs model at {params.dir_run_with_forcing} with bin {params.exe}")
            print("Executing SFINCS...")
            with open (join(params.dir_run_with_forcing,"sfincs.log"), "w") as f:
                subprocess.run([str(params.exe)], stdout=f, cwd=params.dir_run_with_forcing)
                print("Finished running")
        if os.name == 'posix':
            shell("docker image ls")
            shell("docker run --mount src={params.dir_run_with_forcing},target=/data,type=bind deltares/sfincs-cpu:latest sfincs")



rule sfincs_plot_floodmap:
    input:
        mapout = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "sfincs_map.nc"),
    params:
        dir_run = directory(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}")),
        dir_model_no_forcing = directory(join(root_dir, dir_models, "{region}", "{runname}", "sfincs")),
        datacat = get_datacatalog
    output:
        figure = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "plot_output", "sfincs_basemap.png")  
    script:
        join(curdir,  '..', "04_scripts", "postprocessing", "sfincs", "sfincs_postprocess.py")

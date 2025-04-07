### Import some useful python libraries
import os
from snakemake.io import Wildcards
from snakemake import shell
from os.path import join
from itertools import product

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

def get_tcname(wildcards):
    return config["runname_ids"][wildcards.runname]['tc_name']

def get_wind_forcing(wildcards):
    # Returns the data catalog handle for the .spw file 
    return config['runname_ids'][wildcards.runname]['wind_forcing']

def get_utmzone(wildcards):
    return config['runname_ids'][wildcards.runname]['utmzone']

def get_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return [
            join(curdir, '..', "03_data_catalogs", "datacatalog_general.yml"), 
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling.yml"), 
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_obspoints.yml"),
            join(curdir, '..', "03_data_catalogs", "datacatalog_CF_forcing.yml")
        ]
    elif os.name == "posix": #Running on linux
        return [
            join(curdir, '..', "03_data_catalogs", "datacatalog_general___linux.yml"), 
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling___linux.yml"), 
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_obspoints___linux.yml"),
            join(curdir, '..', "03_data_catalogs", "datacatalog_CF_forcing___linux.yml")
        ]

def get_use_dfm(wildcards):
    return config['runname_ids'][wildcards.runname]['use_dfm']

def get_coastal_ts(wildcards):
    return config['runname_ids'][wildcards.runname]['coastal_ts']

def get_cf_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return join(curdir, '..', '03_data_catalogs', 'datacatalog_CF_forcing.yml')
    elif os.name == "posix": #Running on linux
        return join(curdir, '..', '03_data_catalogs', 'datacatalog_CF_forcing___linux.yml')

runname_ids = list(config['runname_ids'].keys())
region = [value['region'] for key, value in config['runname_ids'].items()]
precip_forcing = [value['precip_forcing'] for key, value in config['runname_ids'].items()]
tidemodel = [value['tidemodel'] for key, value in config['runname_ids'].items()]
wind_forcing = [value['wind_forcing'] for key, value in config['runname_ids'].items()]
CF_rain = [value['CF_value_rain'] for key, value in config['runname_ids'].items()]
CF_SLR = [value['CF_value_SLR'] for key, value in config['runname_ids'].items()]
CF_wind = [value['CF_value_wind'] for key, value in config['runname_ids'].items()]

# To prevent unwanted wildcard underscore splitting
wildcard_constraints:
    precip_forcing='|'.join([re.escape(x) for x in precip_forcing]),
    wind_forcing='|'.join([re.escape(x) for x in wind_forcing]),
    CF_SLR=r"-?\d*\.?\d+",  # Matches integer and floating-point numbers (positive and negative)
    CF_wind=r"-?\d*\.?\d+",
    CF_rain=r"-?\d*\.?\d+",

run_combinations = []
for key, value in config['runname_ids'].items():
    for tp, slr, wind in product(value['CF_value_rain'], value['CF_value_SLR'], value['CF_value_wind']):
        run_combinations.append((value['region'], key, value['dfm_res'], value['bathy'], value['precip_forcing'], tp, value['tidemodel'], slr, value['wind_forcing'], wind))

# Unpack into separate wildcard lists
region, runname_ids, dfm_res, bathy, precip_forcing, CF_rain, tidemodel, CF_SLR, wind_forcing, CF_wind = zip(*run_combinations)

rule all_sfincs_update:
    input:
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "plot_output", "sfincs_basemap.png"), zip, region=region, runname=runname_ids, precip_forcing=precip_forcing, CF_rain=CF_rain, tidemodel=tidemodel, CF_SLR=CF_SLR, wind_forcing=wind_forcing, CF_wind=CF_wind),
        # expand(join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "sfincs.dis"), zip, region=region, runname=runname_ids, precip_forcing=precip_forcing, CF_rain=CF_rain, tidemodel=tidemodel, CF_SLR=CF_SLR, wind_forcing=wind_forcing, CF_wind=CF_wind),

rule add_forcing_coastal_meteo_sfincs:
    input:
        msk_file = join(root_dir, dir_models, "{region}", "{runname}", "sfincs", "sfincs.msk"),
    params:
        tc_name = get_tcname,
        dir_run_no_forcing = directory(join(root_dir, dir_models, "{region}", "{runname}", "sfincs")),
        dir_run_with_forcing = lambda wildcards: directory(join(root_dir, dir_runs, wildcards.region, wildcards.runname, "sfincs", 
                                                                  f"event_tp_{wildcards.precip_forcing}_CF{wildcards.CF_rain}_{wildcards.tidemodel}_CF{wildcards.CF_SLR}_{wildcards.wind_forcing}_CF{wildcards.CF_wind}")),
        data_cats = get_datacatalog,
        wind_forcing = get_wind_forcing,
        start_time = get_starttime, 
        end_time = get_endtime,
        use_dfm = get_use_dfm,
        coastal_ts = get_coastal_ts,
        dfm_output = lambda wildcards: f"dfm_output_event_{config['runname_ids'][wildcards.runname]['dfm_res']}_{config['runname_ids'][wildcards.runname]['bathy']}_{wildcards.tidemodel}_CF{wildcards.CF_SLR}_{wildcards.wind_forcing}_CF{wildcards.CF_wind}",
        utmzone = get_utmzone,
        sfincs_obs_points = lambda wildcards: join(root_dir, dir_data, "sfincs_obs_points", config["runname_ids"][wildcards.runname]["sfincs_obs_file"]),
    output:
        bzs_file = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "sfincs.bzs"),
    script:
        join( '..', "04_scripts", "model_building", "sfincs", "update_sfincs_coastal_forcing.py")

rule update_dis_forcing_sfincs:
    input:
        bzs_file = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "sfincs.bzs"),
        wflow_output = join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}", "events", "run_default", "output_scalar.nc"), # do not change!
    params:
        dir_run_with_forcing = lambda wildcards: directory(join(root_dir, dir_runs, wildcards.region, wildcards.runname, "sfincs", 
                                                                  f"event_tp_{wildcards.precip_forcing}_CF{wildcards.CF_rain}_{wildcards.tidemodel}_CF{wildcards.CF_SLR}_{wildcards.wind_forcing}_CF{wildcards.CF_wind}")),
        wflow_root_forcing = directory(join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}")),
        data_cats = get_datacatalog,
    output:
        dis_file = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "sfincs.dis"),
    script:
        join( '..', "04_scripts", "model_building", "sfincs", "update_sfincs_dis_forcing.py")


rule run_sfincs_model:
    input:
        dis_file = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "sfincs.dis"),
    params:
        dir_run_with_forcing = lambda wildcards: directory(join(root_dir, dir_runs, wildcards.region, wildcards.runname, "sfincs", 
                                                                  f"event_tp_{wildcards.precip_forcing}_CF{wildcards.CF_rain}_{wildcards.tidemodel}_CF{wildcards.CF_SLR}_{wildcards.wind_forcing}_CF{wildcards.CF_wind}")),
        exe = join(root_dir, dir_models, "00_executables", "SFINCS_v2.1.1_Dollerup_release_exe", 'sfincs.exe'),
    output:
        mapout = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "sfincs_map.nc"),
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
        mapout = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "sfincs_map.nc"),
    params:
        dir_run = directory(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}")),
        dir_model_no_forcing = directory(join(root_dir, dir_models, "{region}", "{runname}", "sfincs")),
        datacat = get_datacatalog
    output:
        figure   = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "plot_output", "sfincs_basemap.png"),  
        floodmap = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "plot_output", "floodmap.tif") 
    script:
        join(curdir,  '..', "04_scripts", "postprocessing", "sfincs", "sfincs_postprocess.py")

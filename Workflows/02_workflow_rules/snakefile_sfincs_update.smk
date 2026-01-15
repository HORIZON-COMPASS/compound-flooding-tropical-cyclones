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

def get_use_waves(wildcards):
    return config['runname_ids'][wildcards.runname]['use_waves']

def get_coastal_ts(wildcards):
    return config['runname_ids'][wildcards.runname]['coastal_ts']

def get_cf_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return join(curdir, '..', '03_data_catalogs', 'datacatalog_CF_forcing.yml')
    elif os.name == "posix": #Running on linux
        return join(curdir, '..', '03_data_catalogs', 'datacatalog_CF_forcing___linux.yml')


# Generate all combinations of runs based on config parameters
run_combinations = []

for runname, value in config['runname_ids'].items():
    region         = value['region']
    precip_forcing = value['precip_forcing']
    tidemodel      = value['tidemodel']
    wind_forcing   = value['wind_forcing']

    # Extract CF values from config
    rain_low,  rain_high  = value['CF_rain_uncert']
    rain_med              = value['CF_value_rain'][1]

    wind_low,  wind_high  = value['CF_wind_uncert']
    wind_med              = value['CF_value_wind'][1]

    slr_low,   slr_high   = value['CF_SLR_uncert']
    slr_med               = value['CF_value_SLR'][1]

    # Factual
    run_combinations.append(
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, 0, 0)
    )

    # All drivers (low / medium / high)
    run_combinations.extend([
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_low,  slr_low,  wind_low),
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_med,  slr_med,  wind_med),
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_high, slr_high, wind_high),
    ])

    # Rain only (low / medium / high)
    run_combinations.extend([
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_low,  0, 0),
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_med,  0, 0),
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_high, 0, 0),
    ])

    # Wind only (low / medium / high)
    run_combinations.extend([
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, 0, wind_low),
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, 0, wind_med),
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, 0, wind_high),
    ])

    # SLR only (low / medium / high)
    run_combinations.extend([
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, slr_low,  0),
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, slr_med,  0),
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, slr_high, 0),
    ])

    # Wind + SLR (low / medium / high)
    run_combinations.extend([
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, slr_low,  wind_low),
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, slr_med,  wind_med),
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, slr_high, wind_high),
    ])

    # Medium mixed rain combinations only
    run_combinations.extend([
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_med, 0, wind_med),
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_med, slr_med,  0),
    ])


# Unpack into wildcard lists
region, runname_ids, precip_forcing, tidemodel, wind_forcing, CF_rain, CF_SLR, CF_wind = map(list, zip(*run_combinations))

wildcard_constraints:
    precip_forcing='|'.join([re.escape(x) for x in set(precip_forcing)]),
    wind_forcing='|'.join([re.escape(x) for x in set(wind_forcing)]),
    CF_SLR=r"-?\d*\.?\d+",
    CF_wind=r"-?\d*\.?\d+",
    CF_rain=r"-?\d*\.?\d+",

rule all_sfincs_update:
    input:
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "plot_output", "sfincs_basemap.png"), zip, region=region, runname=runname_ids, precip_forcing=precip_forcing, CF_rain=CF_rain, tidemodel=tidemodel, CF_SLR=CF_SLR, wind_forcing=wind_forcing, CF_wind=CF_wind),
        # expand(join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "sfincs.dis"), zip, region=region, runname=runname_ids, precip_forcing=precip_forcing, CF_rain=CF_rain, tidemodel=tidemodel, CF_SLR=CF_SLR, wind_forcing=wind_forcing, CF_wind=CF_wind),
        
rule add_forcing_coastal_meteo_sfincs:
    input:
        msk_file             = join(root_dir, dir_models, "{region}", "{runname}", "sfincs", "sfincs.msk"),
    params:
        tc_name              = get_tcname,
        dir_run_no_forcing   = directory(join(root_dir, dir_models, "{region}", "{runname}", "sfincs")),
        dir_run_with_forcing = lambda wildcards: directory(join(root_dir, dir_runs, wildcards.region, wildcards.runname, "sfincs", 
                                                                  f"event_tp_{wildcards.precip_forcing}_CF{wildcards.CF_rain}_{wildcards.tidemodel}_CF{wildcards.CF_SLR}_{wildcards.wind_forcing}_CF{wildcards.CF_wind}")),
        data_cats            = get_datacatalog,
        wind_forcing         = get_wind_forcing,
        start_time           = get_starttime, 
        end_time             = get_endtime,
        use_dfm              = get_use_dfm,
        use_waves            = get_use_waves,
        coastal_ts           = get_coastal_ts,
        dfm_output           = lambda wildcards: f"dfm_output_event_{config['runname_ids'][wildcards.runname]['dfm_res']}_{config['runname_ids'][wildcards.runname]['bathy']}_{wildcards.tidemodel}_CF{wildcards.CF_SLR}_{wildcards.wind_forcing}_CF{wildcards.CF_wind}",
        utmzone              = get_utmzone,
        sfincs_obs_points    = lambda wildcards: join(root_dir, dir_data, "sfincs_obs_points", config["runname_ids"][wildcards.runname]["sfincs_obs_file"]),
    output:
        bzs_file             = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "sfincs.bzs"),
    script:
        join( '..', "04_scripts", "model_building", "sfincs", "update_sfincs_coastal_forcing.py")

rule update_dis_forcing_sfincs:
    input:
        bzs_file              = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "sfincs.bzs"),
        wflow_output          = join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}", "events", "run_default", "output_scalar.nc"), # do not change!
        wflow_dis_no_bankfull = join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}", "events", "run_default", "wflow_dis_no_qbankfull.csv"),
    params:
        dir_run_with_forcing  = lambda wildcards: directory(join(root_dir, dir_runs, wildcards.region, wildcards.runname, "sfincs", 
                                                                  f"event_tp_{wildcards.precip_forcing}_CF{wildcards.CF_rain}_{wildcards.tidemodel}_CF{wildcards.CF_SLR}_{wildcards.wind_forcing}_CF{wildcards.CF_wind}")),
        wflow_root_forcing    = directory(join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{precip_forcing}_CF{CF_rain}")),
        wflow_base            = directory(join(root_dir, dir_models, "{region}", "{runname}", "wflow")),
        data_cats             = get_datacatalog,
    output:
        dis_file              = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "sfincs.dis"),
    script:
        join( '..', "04_scripts", "model_building", "sfincs", "update_sfincs_dis_forcing.py")


rule run_sfincs_model:
    threads: 16 # increase when using more vcpu's
    input:
        dis_file             = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "sfincs.dis"),
    params:
        dir_run_with_forcing = lambda wildcards: directory(join(root_dir, dir_runs, wildcards.region, wildcards.runname, "sfincs", 
                                                                  f"event_tp_{wildcards.precip_forcing}_CF{wildcards.CF_rain}_{wildcards.tidemodel}_CF{wildcards.CF_SLR}_{wildcards.wind_forcing}_CF{wildcards.CF_wind}")),
        exe                  = join(root_dir, dir_models, "00_executables", "SFINCS_v2.1.1_Dollerup_release_exe", 'sfincs.exe'),
        sif                  = join(root_dir, dir_models, "00_executables", "sfincs-cpu_latest.sif"),
        sif_dir              = join(root_dir, dir_models, "00_executables"),
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
            shell("docker run --mount src={params.dir_run_with_forcing},target=/data,type=bind deltares/sfincs-cpu:sfincs-v2.2.0-col-dEze-Release")

rule sfincs_plot_floodmap:
    input:
        mapout = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "sfincs_map.nc"),
    params:
        dir_run              = directory(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}")),
        dir_model_no_forcing = directory(join(root_dir, dir_models, "{region}", "{runname}", "sfincs")),
        datacat              = get_datacatalog
    output:
        figure   = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "plot_output", "sfincs_basemap.png"),  
        floodmap = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "plot_output", "floodmap.tif") 
    script:
        join(curdir,  '..', "04_scripts", "postprocessing", "sfincs", "sfincs_postprocess.py")

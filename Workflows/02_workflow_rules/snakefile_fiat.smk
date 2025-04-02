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


def get_tcname(wildcards):
    return config["runname_ids"][wildcards.runname]['tc_name']

def get_wind_forcing(wildcards):
    # Returns the data catalog handle for the .spw file 
    return config['runname_ids'][wildcards.runname]['wind_forcing']

def get_country(wildcards):
    return config['runname_ids'][wildcards.runname]['country']

def get_continent(wildcards):
    return config['runname_ids'][wildcards.runname]['continent']
def get_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return [
            join(curdir, '..', "03_data_catalogs", "datacatalog_fiat.yml")
        ]
    elif os.name == "posix": #Running on linux
        return [
            join(curdir, '..', "03_data_catalogs", "datacatalog_fiat___linux.yml"), 
        ]


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


rule all_fiat_model:
    input:
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "fiat","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "settings.toml"), zip, region=region, runname=runname_ids, precip_forcing=precip_forcing, CF_rain=CF_rain, tidemodel=tidemodel, CF_SLR=CF_SLR, wind_forcing=wind_forcing, CF_wind=CF_wind),


rule build_fiat_model:
    input:
        floodmap             = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "plot_output", "floodmap.tif")
    params:
        dir_run_with_forcing = lambda wildcards: directory(join(root_dir, dir_runs, wildcards.region, wildcards.runname, "sfincs", 
                                                                  f"event_tp_{wildcards.precip_forcing}_CF{wildcards.CF_rain}_{wildcards.tidemodel}_CF{wildcards.CF_SLR}_{wildcards.wind_forcing}_CF{wildcards.CF_wind}")),
        datacat_fiat         = get_datacatalog
        model_folder         = directory(join(root_dir, dir_runs, "{region}", "{runname}", "fiat")),
        continent            = get_continent
        country              = get_country
        config               = config_file = join(curdir,'..', "05_config_models", "03_fiat", "config_fiat.yml"),
    output:
        fiat_settings        = join(root_dir,  dir_runs, "{region}", "{runname}", "fiat","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "settings.toml"),
    script:
        join( '..', "04_scripts", "model_building", "fiat", "setup_fiat.py")


rule run_fiat_model:
    input:
        fiat_settings = join(root_dir, dir_runs, "{region}", "{runname}", "fiat",
                             "event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "settings.toml"),
    params:
        dir_run_with_forcing = lambda wildcards: directory(join(root_dir, dir_runs, wildcards.region, wildcards.runname, "sfincs", 
                                                                  f"event_tp_{wildcards.precip_forcing}_CF{wildcards.CF_rain}_{wildcards.tidemodel}_CF{wildcards.CF_SLR}_{wildcards.wind_forcing}_CF{wildcards.CF_wind}")),
    output:
        mapout = join(root_dir, dir_runs, "{region}", "{runname}", "fiat", 
                      f"event_tp_{wildcards.precip_forcing}_CF{wildcards.CF_rain}_{wildcards.tidemodel}_CF{wildcards.CF_SLR}_{wildcards.wind_forcing}_CF{wildcards.CF_wind}", "fiat_output.nc"),
    run:
        if os.name == 'nt':  # For Windows
            import subprocess
            print(f"Running FIAT model with settings: {input.fiat_settings}")
            with open(join(params.dir_run_with_forcing, "fiat_log.txt"), "w") as log_file:
                subprocess.run(["fiat", "run", input.fiat_settings], stdout=log_file, stderr=subprocess.PIPE, cwd=params.dir_run_with_forcing)
            print("Finished running FIAT on Windows.")
        
        elif os.name == 'posix':  # For Linux
            print(f"Running FIAT model with settings: {input.fiat_settings}")
            shell(f"fiat run {input.fiat_settings}")
            print("Finished running FIAT on Linux.")
#%%### Import some useful python libraries
import os
from snakemake.io import Wildcards
from os.path import join
from itertools import product

curdir = os.getcwd()
if os.name == 'nt': #Running on windows
    root_dir = join("p:/",config['root_dir'])
elif os.name == "posix": #Running on linux
    root_dir = join("/p", config['root_dir'])
dir_runs = config['dir_runs']
dir_models = config['dir_models']

def get_region(wildcards):
    print(test)
    return config["runname_ids"][wildcards.runname]['region']

def get_tcname(wildcards):
    return config["runname_ids"][wildcards.runname]['tc_name']

def get_start_time(wildcards):
    return config["runname_ids"][wildcards.runname]['start_time']

def get_end_time(wildcards):
    return config["runname_ids"][wildcards.runname]['end_time']

def get_meteo_forcing(wildcards):
    # dataset used for wflow temperature/PET in the event. Defaults to precip_forcing (works when
    # that is a full met dataset like era5), but can be set separately when the precip product is
    # precip-only (e.g. ceh_gear radar) via a `meteo_forcing` config field.
    rid = config["runname_ids"][wildcards.runname]
    return rid.get('meteo_forcing', rid['precip_forcing'])

def get_bbox(wildcards):
    prebbox = config["runname_ids"][wildcards.runname]["bbox_sfincs"]
    arg_bbox = "{" + "'bbox': "+ prebbox + "}"
    return arg_bbox

def get_config_wflow(wildcards):
    config_wflow_base = config["runname_ids"][wildcards.runname]['config_wflow_base']
    return join(curdir, '..', "05_config_models", "01_wflow", config_wflow_base)

def get_river_upa(wildcards):
    return config["runname_ids"][wildcards.runname]["river_upa"]

def get_lulc_mapping(wildcards):
    # Reclassification table mapping land-use classes to wflow parameters. Distinct from the
    # SFINCS table: wflow needs the *_hydromtwflow variant (many parameters, not just Manning).
    # The land-use dataset itself is the {CF_landuse} path wildcard.
    return config["runname_ids"][wildcards.runname].get("lulc_mapping_wflow")

# def get_dir_model_base(wildcards):
#     print(wildcards)
#     return join(root_dir, dir_models, config["runname_ids"][wildcards.runname]['region'], config["runname_ids"][wildcards.runname], "wflow")

def get_use_bankfull_corr(wildcards):
    return config["runname_ids"][wildcards.runname]["use_bankfull_corr"]

def get_landuse_30yr_wflow(wildcards):
    # Land use of the decoupled 30-yr run used for the bankfull estimate. Falls back to the
    # event's own land use when the case does not define a separate one.
    rv = config["runname_ids"][wildcards.runname]
    return rv.get("bankfull_corr_lulc", rv["CF_landuse"][0])

def get_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return [
            join(curdir, '..', "03_data_catalogs", "datacatalog_general.yml"), 
            join(curdir, '..', "03_data_catalogs", "datacatalog_CF_forcing.yml")
        ]
    elif os.name == "posix": #Running on linux
        return [
            join(curdir, '..', "03_data_catalogs", "datacatalog_general_v1___linux.yml"),
            join(curdir, '..', "03_data_catalogs", "datacatalog_CF_forcing_v1___linux.yml")
        ]

runname_ids = list(config['runname_ids'].keys())
region = [value['region'] for key, value in config['runname_ids'].items()]
precip_forcing = [value['precip_forcing'] for key, value in config['runname_ids'].items()]
CF_rain = [value['CF_value_rain'] for key, value in config['runname_ids'].items()]

# To prevent unwanted wildcard underscore splitting
wildcard_constraints:
    precip_forcing='|'.join([re.escape(x) for x in precip_forcing]),
    CF_rain=r"-?\d*\.?\d+", # Matches integer and floating-point numbers (positive and negative)

run_combinations = []
for key, value in config['runname_ids'].items():
    # Model/run directories are suffixed with the land-use scenario (wflow_{CF_landuse}),
    # so several land-use variants of the same run can coexist.
    for tp, lulc in product(value['CF_value_rain'], value['CF_landuse']):
        run_combinations.append((value['region'], key, value['precip_forcing'], tp, lulc,
                                 value['use_bankfull_corr']))

# Unpack into separate wildcard lists
region, runname_ids, precip_forcing, CF_rain, CF_landuse, bankfull_corr = zip(*run_combinations)

rule all_wflow:
    input:
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}", "event_precip_{precip_forcing}_CF{CF_rain}", "events", "run_default", "output_scalar.nc"), zip, region=region, runname=runname_ids, precip_forcing=precip_forcing, CF_rain=CF_rain, CF_landuse=CF_landuse),
        # Only runs with use_bankfull_corr enabled request the corrected discharge series.
        # The file is wflow_dis_no_bankfull.csv and is what
        # rule postprocess_discharge produces.
        [fn for fn, use_bf in zip(
            expand(join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}", "event_precip_{precip_forcing}_CF{CF_rain}", "events", "run_default", "wflow_dis_no_bankfull.csv"), zip, region=region, runname=runname_ids, precip_forcing=precip_forcing, CF_rain=CF_rain, CF_landuse=CF_landuse),
            bankfull_corr) if use_bf],

rule make_base_model_wflow:
    input:
        #config_file = join(curdir,'..', "05_config_models", "01_wflow", "config_wflow.yml"),
        region_geom = join(root_dir, dir_models, "{region}", "{runname}", "sfincs_{CF_landuse}", "gis", "region.geojson"),
        dir_sfincs_model = join(root_dir, dir_models, "{region}", "{runname}", "sfincs_{CF_landuse}"),
        src_file = join(root_dir, dir_models, "{region}", "{runname}", "sfincs_{CF_landuse}", "gis", "dis.geojson"),  # v2 renamed src.geojson -> dis.geojson
        config_file = get_config_wflow  
    params:
        dir_model = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{CF_landuse}"),
        data_cat = get_datacatalog,
        arg_bbox = get_bbox,
        river_upa = get_river_upa,
        lulc_mapping_wflow = get_lulc_mapping
    output: 
        toml_file = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{CF_landuse}", 'wflow_sbm.toml'),
        staticmaps = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{CF_landuse}", 'staticmaps.nc'), 
    script:
        join(curdir, '..', "04_scripts", "model_building", "wflow", "setup_wflow_base.py")

# update wflow forcing for warmup
rule update_forcing_wflow_warmup:
    input: 
        toml_file = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{CF_landuse}", 'wflow_sbm.toml'),
        staticmaps = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{CF_landuse}", 'staticmaps.nc'), 
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "warmup", "inmaps.nc"),
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "warmup", "wflow_sbm.toml"),
    params:
        wflow_root_noforcing = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{CF_landuse}"),
        wflow_root_forcing= join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}"),
        start_time = get_start_time,
        end_time = get_end_time,
        data_cat = get_datacatalog,
    script:
        join(curdir, '..', "04_scripts", "model_building", "wflow", "update_forcing_wflow_warmup.py")

# update wflow forcing for event
rule update_forcing_wflow_event:
    input: 
        toml_file = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{CF_landuse}", 'wflow_sbm.toml'),
        staticmaps = join(root_dir, dir_models, "{region}", "{runname}", "wflow_{CF_landuse}", 'staticmaps.nc'), 
        previous_rule = join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "warmup", "inmaps.nc")
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "events", "inmaps.nc"),
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "events", "wflow_sbm.toml"),
    params:
        wflow_root_noforcing = directory(join(root_dir, dir_models, "{region}", "{runname}", "wflow_{CF_landuse}")),
        wflow_root_forcing= directory(join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}")),
        start_time = get_start_time,
        end_time = get_end_time,
        forcing = get_meteo_forcing,
        data_cat = get_datacatalog,
        tc_name = get_tcname
    script:
        join(curdir, '..',  "04_scripts", "model_building", "wflow", "update_forcing_wflow_event.py")

rule run_wflow_warmup:
    threads: 16
    input:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "warmup", "inmaps.nc"),
        toml = join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "warmup", "wflow_sbm.toml"),
        previous_rule = join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "events", "inmaps.nc"),  
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "events", "instate", "instates.nc"),
    params:
        exe = join(root_dir, dir_models, "00_executables", "wflow0.8.1", "wflow_cli", "bin", "wflow_cli.exe"),
        julia_env_fn = "~/.julia/environments/wflow1"  # Wflow.jl 1.0.2 (v1); v1.9 env has v0.8
    shell:
        # The wflow0.8.1 exe cannot run a Wflow.jl v1 TOML, so on Linux this falls through to
        # Julia with the wflow1 env (Wflow 1.0.2). Uses the juliaup default Julia (1.11.3, which
        # matches the wflow1 manifest); the old `+1.9` channel is not installed here.
        """
        export PATH="$HOME/.juliaup/bin:$PATH"
        {params.exe} {input.toml} || julia --threads 4 --project={params.julia_env_fn} -e "using Wflow; Wflow.run()" "{input.toml}"
        """

rule run_wflow_event:
    threads: 16
    input:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "events", "instate", "instates.nc"),
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "events", "inmaps.nc"),
        toml = join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "events", "wflow_sbm.toml"),
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "events", "run_default", "output_scalar.nc"),
    params:
        exe = join(root_dir, dir_models, "00_executables", "wflow0.8.1", "wflow_cli", "bin", "wflow_cli.exe"),
        julia_env_fn = "~/.julia/environments/wflow1",  # Wflow.jl 1.0.2 (v1); v1.9 env has v0.8
    shell:
        """
        export PATH="$HOME/.juliaup/bin:$PATH"
        {params.exe} {input.toml} || julia --threads 4 --project={params.julia_env_fn} -e "using Wflow; Wflow.run()" "{input.toml}"
        """

# remove bankfull discharge.
# The 30-yr run is decoupled: it is located under its own land-use directory
# (bankfull_corr_lulc), which need not match the event run's CF_landuse, so it is passed
# as a param rather than declared as an input.
rule postprocess_discharge:
    input:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "events", "run_default", "output_scalar.nc"),
    output:
        join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}", "events", "run_default", "wflow_dis_no_bankfull.csv")
    params:
        wflow_root_forcing_30yr = directory(join(root_dir, dir_runs, "{region}", "{runname}")),
        wflow_root_forcing = directory(join(root_dir, dir_runs, "{region}", "{runname}", "wflow_{CF_landuse}","event_precip_{precip_forcing}_CF{CF_rain}")),
        landuse_30yr = get_landuse_30yr_wflow,
        use_bankfull_corr = get_use_bankfull_corr,
        data_cat = get_datacatalog,
    script: join(curdir, '..',  "04_scripts", "postprocessing", "wflow", "calculate_and_remove_qbankfull.py")

### Import some useful python libraries
import os
from snakemake.io import Wildcards
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

def get_bathy(wildcards):
    bathy = config["runname_ids"][wildcards.runname]["bathy"]
    return bathy

def get_dfm_coastal_mask(wildcards):
    dfm_coastal_mask = config["runname_ids"][wildcards.runname]["dfm_coastal_mask"]
    return dfm_coastal_mask 

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

run_combinations = []
for key, value in config['runname_ids'].items():
    for lulc in value['CF_landuse']:
        run_combinations.append({
            "region": value['region'],
            "runname": key,
            "landuse": lulc
        })

regions     = [c["region"] for c in run_combinations]
runname_ids = [c["runname"] for c in run_combinations]
landuses    = [c["landuse"] for c in run_combinations]

rule all_sfincs_build:
    input:
        expand(join(root_dir, dir_models, "{region}", "{runname}", "sfincs_{landuse}", "sfincs.msk"), zip, region=regions, runname=runname_ids, landuse=landuses),
        expand(join(root_dir, dir_models, "{region}", "{runname}", "sfincs_{landuse}", "gis", "src.geojson"), zip, region=regions, runname=runname_ids, landuse=landuses)

rule make_base_model_sfincs:
    params:
        arg_bbox = get_bbox,
        dir_model_sfincs = join(root_dir, dir_models, "{region}", "{runname}", "sfincs_{landuse}"),
        data_cats = get_datacatalog,
        bathy = get_bathy,
        dfm_coastal_mask = get_dfm_coastal_mask,
    input:
        config_file = join(curdir,'..', "05_config_models", "02_sfincs", "sfincs_base_build_{landuse}.yml"),
    output: 
        dir_sfincs_model = directory(join(root_dir, dir_models, "{region}", "{runname}", "sfincs_{landuse}")),
        msk_file = join(root_dir, dir_models, "{region}", "{runname}", "sfincs_{landuse}" , "sfincs.msk"),
        src_file = join(root_dir, dir_models, "{region}", "{runname}", "sfincs_{landuse}", "gis", "src.geojson"),
        region_geom = join(root_dir, dir_models, "{region}", "{runname}", "sfincs_{landuse}", "gis", "region.geojson"),
    script:
        join("..","04_scripts", "model_building", "sfincs", "setup_sfincs_base.py")

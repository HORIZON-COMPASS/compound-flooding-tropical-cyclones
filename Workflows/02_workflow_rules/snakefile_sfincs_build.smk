### Import some useful python libraries
import os
from snakemake.io import Wildcards
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

def get_config(wildcards):
    config_sfincs_base = config["runname_ids"][wildcards.runname]['config_sfincs_base']
    return join(curdir, '..', "05_config_models", "02_sfincs", config_sfincs_base)

def get_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return [
            join(curdir, '..', "03_data_catalogs", "datacatalog_general.yml"), 
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling.yml"), 
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_obspoints.yml")
        ]
    elif os.name == "posix": #Running on linux
        return [
            join(curdir, '..', "03_data_catalogs", "datacatalog_general___linux.yml"), 
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_coastal_coupling___linux.yml"), 
            join(curdir, '..', "03_data_catalogs", "datacatalog_SFINCS_obspoints___linux.yml")
        ]

runname_ids = list(config['runname_ids'].keys())
regions = [value['region'] for key, value in config['runname_ids'].items()]

rule all_sfincs_build:
    input:
        expand(join(root_dir, dir_models, "{region}", "{runname}", "sfincs" , "sfincs.msk"), zip, region=regions, runname=runname_ids),
        expand(join(root_dir, dir_models, "{region}", "{runname}", "sfincs", "gis", "src.geojson"), zip, region=regions, runname=runname_ids)

rule make_base_model_sfincs:
    params:
        arg_bbox = get_bbox,
        dir_model_sfincs = join(root_dir, dir_models, "{region}", "{runname}", "sfincs"),
        data_cats = get_datacatalog
    input:
        config_file = get_config,
    output: 
        dir_sfincs_model = directory(join(root_dir, dir_models, "{region}", "{runname}", "sfincs")),
        msk_file = join(root_dir, dir_models, "{region}", "{runname}", "sfincs" , "sfincs.msk"),
        src_file = join(root_dir, dir_models, "{region}", "{runname}", "sfincs", "gis", "src.geojson"),
        region_geom = join(root_dir, dir_models, "{region}", "{runname}", "sfincs", "gis", "region.geojson"),
    script:
        join("..","04_scripts", "model_building", "sfincs", "setup_sfincs_base.py")

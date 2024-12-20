### Import some useful python libraries
import os
from snakemake.io import Wildcards
from os.path import join

curdir = os.getcwd()
if os.name == 'nt': #Running on windows
    root_dir = join("p:/",config['root_dir'])
elif os.name == "posix": #Running on linux
    root_dir = join("/p", config['root_dir'])

def get_bbox(wildcards):
    bbox = config["runname_ids"][wildcards.runname]["bbox"]
    return bbox

def get_path_spw_ori(wildcards):
    file_spw = config["runname_ids"][wildcards.runname]["file_spw"]
    path_spw_ori = join(root_dir, config['dir_spw'], file_spw)
    return path_spw_ori


def get_file_spw(wildcards):
    # Returns the file name for the .spw file (e.g., tc_FREDDY_2023061S22036_ext9d.spw)
    return config['runname_ids'][wildcards.runname]['file_spw']

def get_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return [
            join(curdir, "data_catalogs", "datacatalog_general.yml"), 
            join(curdir, "data_catalogs", "datacatalog_SFINCS_coastal_coupling.yml"), 
            join(curdir, "data_catalogs", "datacatalog_SFINCS_obspoints.yml")
        ]
    elif os.name == "posix": #Running on linux
        return [
            join(curdir, "data_catalogs", "datacatalog_general___linux.yml"), 
            join(curdir, "data_catalogs", "datacatalog_SFINCS_coastal_coupling___linux.yml"), 
            join(curdir, "data_catalogs", "datacatalog_SFINCS_obspoints___linux.yml")
        ]

regions = [value['region'] for key, value in config['runname_ids'].items()]
spw_files = [value['file_spw'] for key, value in config['runname_ids'].items()]
runname_ids = list(config['runname_ids'].keys())  #
forcing = [value['forcing'] for key, value in config['runname_ids'].items()]

rule all_sfincs_build:
    input:
        expand(join(root_dir, "02_Models", "{region}", "{runname}", "sfincs" , "sfincs.msk"), zip, region=regions, runname=runname_ids),
        expand(join(root_dir, "02_Models", "{region}", "{runname}", "sfincs", "gis", "src.geojson"), zip, region=regions, runname=runname_ids)

rule make_base_model_sfincs:
    params:
        arg_bbox = get_bbox,
        dir_model_sfincs = join(root_dir, "02_Models", "{region}", "{runname}", "sfincs"),
        data_cats = get_datacatalog
    input:
        config_file = join(curdir, "config_sfincs", "sfincs_base_build.yml"),
    output: 
        msk_file = join(root_dir, "02_Models", "{region}", "{runname}", "sfincs" , "sfincs.msk"),
        src_file = join(root_dir, "02_Models", "{region}", "{runname}", "sfincs", "gis", "src.geojson")
    script:
        join("scripts", "model_building", "sfincs", "setup_sfincs_base.py")

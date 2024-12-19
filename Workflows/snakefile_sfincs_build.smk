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
        expand(join(root_dir, "02_Models", "{region}", "{runname}", "sfincs" , "sfincs.msk"), zip, region=regions, runname=runname_ids)

rule make_base_model_sfincs:
    params:
        arg_bbox = get_bbox,
        dir_model_sfincs = join(root_dir, "02_Models", "{region}", "{runname}", "sfincs"),
        data_cats = get_datacatalog
    input:
        config_file = join(curdir, "config_sfincs", "sfincs_base_build.yml"),
    output: 
        msk_file = join(root_dir, "02_Models", "{region}", "{runname}", "sfincs" , "sfincs.msk")
    script:
        join("scripts", "model_building", "sfincs", "setup_sfincs_base.py")


# rule add_forcing:
#     input:
#         msk_file = join(root_dir, "02_Models", "{region}", "{runname}", "sfincs" , "sfincs.msk"),
#         spw_file_in = get_path_spw_ori,
#     params:
#         dir_run_no_forcing = join(root_dir, "02_Models", "{region}", "{runname}", "sfincs"),
#         dir_run_with_forcing = join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "sfincs"),
#         forcing_yml = "config_sfincs/sfincs_"+"{runname}"+"_forcing.yml",
#         data_cats = get_datacatalog,
#     output:
#         bzs_file = join(root_dir,  "03_Runs", "{region}", "{runname}", "{forcing}", "sfincs", "sfincs.bzs"),
#     script:
#         join("scripts", "preprocessing", "update_sfincs_coastal_forcing.py")


# # rule - check the inp file? for e.g: formatting in linux


# # rule add_batchfile:
# #     input:
# #         bzs_file = "{dir_run}"+"/sfincs_"+"{runname}"+"/sfincs.bzs" 
# #     params:
# #         dir_run = "{dir_run}"+"/sfincs_"+"{runname}"
# #     output:
# #         ("{dir_run}"+"/sfincs_"+"{runname}"+"/run_sfincs.bat")
# #     script:
# #         'add_batchfile.py'

# rule run_sfincs_model:
#     input:
# #        batchfile = "{dir_run}"+"/sfincs_"+"{runname}"+"/run_sfincs.bat"
#         bzs_file = join(root_dir,  "03_Runs", "{region}", "{runname}", "{forcing}", "sfincs", "sfincs.bzs") 
#     params:
#         dir_run_with_forcing = join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "sfincs"),
#         currentdir = curdir
#     output:
#         mapout = join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "sfincs", "sfincs_map.nc"),
#         hisout =join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "sfincs", "sfincs_his.nc"),     
#     shell:
#         '''
#         docker image ls 
#         docker run --mount src={params.dir_run_with_forcing},target=/data,type=bind deltares/sfincs-cpu:latest sfincs
#         '''
# #        '(cd {params.dir_run} && run_sfincs.bat)

# rule sfincs_plot_floodmap:
#     input:
#         mapout = join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "sfincs", "sfincs_map.nc"),
#         hisout =join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "sfincs", "sfincs_his.nc"), 
#     params:
#         dir_run = join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "sfincs"),
#         datacat = get_datacatalog
#     output:
#         figure = join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "sfincs", "plot_output", "sfincs_basemap.png")  
#     script:
#         join(curdir, "scripts", "postprocessing", "sfincs_postprocess.py")

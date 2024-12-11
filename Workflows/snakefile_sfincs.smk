### Import some useful python libraries
import os
from snakemake.io import Wildcards

curdir = os.getcwd()

def get_bbox(wildcards):
    prebbox = config["runname_ids"][wildcards.runname]["bbox"]
    arg_bbox = "{" + "'bbox': "+ prebbox + "}"
    return arg_bbox

def get_file_spw(wildcards):
    return config["runname_ids"][wildcards.runname]["file_spw"]

def get_path_spw_ori(wildcards):
    file_spw = config["runname_ids"][wildcards.runname]["file_spw"]
    path_spw_ori = config['dir_spw']+"/"+ file_spw
    return path_spw_ori

def get_path_spw_run(wildcards):
    file_spw = config["runname_ids"][wildcards.runname]["file_spw"]
    path_spw_run = config['dir_runs'] + "/sfincs_" + config["runname_ids"][wildcards.runname]["tc_name"] + "/"+ file_spw
    return path_spw_run

def get_datacatalog(wildcards):
    return config['datacatalog']

rule all:
    input:
        expand(("{dir_run}"+"/sfincs_"+"{runname}"+"/plot_output/sfincs_basemap.png"), dir_run=config['dir_runs'], runname= config["runname_ids"])

rule make_base_model:
    params:
        dir_run = "{dir_run}"+"/sfincs_"+"{runname}",
        rmfile = "{dir_run}"+"/sfincs_"+"{runname}/"+"hydromt_data.yml",
        arg_bbox = get_bbox
    output: 
        msk_file = "{dir_run}"+"/sfincs_{runname}"+"/sfincs.msk"
    shell:
        '''
        hydromt build sfincs {params.dir_run} --region "{params.arg_bbox}" -i config_sfincs/sfincs_base_build.yml --force-overwrite -v 
        del {params.rmfile} || rm {params.rmfile}
        '''

rule add_forcing:
    input:
        msk_file = "{dir_run}"+"/sfincs_"+"{runname}"+"/sfincs.msk"
    params:
        dir_run = "{dir_run}"+"/sfincs_"+"{runname}",
        forcing_yml = "config_sfincs/sfincs_"+"{runname}"+"_forcing.yml",
        path_spw_ori = get_path_spw_ori,
        path_spw_run = get_path_spw_run
    output:
        bzs_file = "{dir_run}"+"/sfincs_"+"{runname}"+"/sfincs.bzs" 
    shell:
        '''
        cp {params.path_spw_ori} {params.path_spw_run}
        hydromt update sfincs {params.dir_run} -i {params.forcing_yml} -v
        '''

# rule - check the inp file? for e.g: formatting in linux


# rule add_batchfile:
#     input:
#         bzs_file = "{dir_run}"+"/sfincs_"+"{runname}"+"/sfincs.bzs" 
#     params:
#         dir_run = "{dir_run}"+"/sfincs_"+"{runname}"
#     output:
#         ("{dir_run}"+"/sfincs_"+"{runname}"+"/run_sfincs.bat")
#     script:
#         'add_batchfile.py'

rule run_sfincs_model:
    input:
#        batchfile = "{dir_run}"+"/sfincs_"+"{runname}"+"/run_sfincs.bat"
        bzs_file = "{dir_run}"+"/sfincs_"+"{runname}"+"/sfincs.bzs"
    params:
        dir_run = "{dir_run}"+"/sfincs_"+"{runname}",
        currentdir = curdir
    output:
        mapout = "{dir_run}"+"/sfincs_"+"{runname}"+"/sfincs_map.nc",
        hisout = "{dir_run}"+"/sfincs_"+"{runname}"+"/sfincs_his.nc"     
    shell:
        '''
        docker image ls 
        docker run --mount src={params.dir_run},target=/data,type=bind deltares/sfincs-cpu:latest sfincs
        '''
#        '(cd {params.dir_run} && run_sfincs.bat)

rule sfincs_plot_floodmap:
    input:
        mapfile = "{dir_run}"+"/sfincs_"+"{runname}"+"/sfincs_map.nc",
        hisfile = "{dir_run}"+"/sfincs_"+"{runname}"+"/sfincs_his.nc"
    params:
        dir_run = "{dir_run}"+"/sfincs_"+"{runname}",
        datacat = get_datacatalog
    output:
        figure = "{dir_run}"+"/sfincs_"+"{runname}"+"/plot_output/sfincs_basemap.png"  
    script:
        "postprocess/sfincs_postprocess.py" 

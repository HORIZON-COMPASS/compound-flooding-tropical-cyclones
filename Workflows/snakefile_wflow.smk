### Import some useful python libraries
import os
from snakemake.io import Wildcards
from os.path import join

curdir = os.getcwd()
root_dir = config['root_dir']
event = config['event_name']


def get_bbox(wildcards):
    prebbox = config["runname_ids"][wildcards.runname]["bbox"]
    arg_bbox = "{" + "'bbox': "+ prebbox + "}"
    return arg_bbox

def get_datacatalog(wildcards):
    return config['datacatalog']

#TODO Update rules with expand functions
rule all:
    input:
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "event", "run_default", "output_scalar.nc"),

rule make_base_model:
    input:
        config_file = join(curdir, "config_wflow", "wflow_build_sofala.yml"),
        region_geom = join(root_dir, "02_Models", region, experiment_name, "SFINCS_noforcing", "gis", "basin.geojson")
    params:
        dir_model = join(root_dir, "02_Models", region, experiment_name, "wflow"),
        data_cat = join(cur_dir, "data_catalogs", "datacatalog_general.yml")
        arg_bbox = get_bbox,
    output: 
        toml_file = join(dir_model, 'wflow_sbm.toml'),
        staticmaps = join(dir_model, 'staticmaps.nc'),
    script:
        join("scripts", "model_building", "wflow")



# update wflow forcing for warmup
rule update_forcing_wflow_warmup:
    input: 
        toml_file = join(dir_model, 'wflow_sbm.toml'),
        staticmaps = join(dir_model, 'staticmaps.nc'),
    output:
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "warmup", "inmaps_warmup.nc"),
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "warmup", "wflow_sbm.toml"),
    params:
	    wflow_root = wflow_root,
    script:
        join(script_root, "wflow", "update_forcing_wflow_warmup.py")

# update wflow forcing for event
rule update_forcing_wflow_event:
    input: 
        toml_file = join(dir_model, 'wflow_sbm.toml'),
        staticmaps = join(dir_model, 'staticmaps.nc'),
    output:
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "event", "inmaps.nc"),
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "event", "wflow_sbm.toml"),
    params:
	    wflow_root = wflow_root,
    script:
        join(script_root, "wflow", "update_forcing_wflow_event.py")

rule run_wflow_warmup:
    input:
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "warmup", "inmaps_warmup.nc"),
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "warmup", "wflow_sbm.toml"),
    output:
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "warmup", "wflow_sbm.toml"),
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "warmup", "outstate", "outstates.nc"),
    params:
        exe = join(root_dir, "02_Models", "wflow0.8.1", "wflow_cli", "bin", "wflow_cli.exe"),
    shell:
        "{params.exe} {input.toml}"

rule copy_state_wflow:
    input:
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "warmup", "outstate", "outstates.nc"),
    output:
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "event", "instate", "instates.nc"),
    shell:
       "copy {input} {output}"

rule run_wflow_event:
    input:
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "event", "instate", "instates.nc"),
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "event", "inmaps.nc"),
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "event", "wflow_sbm.toml"),
    output:
        join(root_dir, "03_Runs", region, experiment_name, forcing, "wflow", "event", "run_default", "output_scalar.nc"),
    params:
        exe = join(root_dir, "02_Models", "wflow0.8.1", "wflow_cli", "bin", "wflow_cli.exe"),
    shell:
        "{params.exe} {input.toml}"
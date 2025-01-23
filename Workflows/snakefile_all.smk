### Import some useful python libraries
import os
from snakemake.io import Wildcards
from os.path import join

include: 'snakefile_sfincs_build.smk'
include: 'snakefile_wflow.smk'
# include: 'snakefile_dfm.smk' # to be included and tested
include: 'snakefile_sfincs_update.smk'
curdir = os.getcwd()
if os.name == 'nt': #Running on windows
    root_dir = join("p:/",config['root_dir'])
elif os.name == "posix": #Running on linux
    root_dir = join("/p", config['root_dir'])
dir_runs = config['dir_runs']



regions = [value['region'] for key, value in config['runname_ids'].items()]
runname_ids = list(config['runname_ids'].keys())  #
forcing = [value['forcing'] for key, value in config['runname_ids'].items()]

rule all_general:
    input:
        expand(join(root_dir, "03_Runs", "{region}", "{runname}", "{forcing}", "sfincs", "plot_output", "sfincs_basemap.png"), zip, region=regions, runname=runname_ids, forcing=forcing, spw_file = spw_files),


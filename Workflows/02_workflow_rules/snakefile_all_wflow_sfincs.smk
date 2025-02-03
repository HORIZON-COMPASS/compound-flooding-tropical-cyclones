### Import some useful python libraries
import os
from snakemake.io import Wildcards
from os.path import join

curdir = os.getcwd()

include: 'snakefile_sfincs_build.smk'
include: 'snakefile_wflow.smk'
include: 'snakefile_sfincs_update.smk'

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
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "plot_output", "sfincs_basemap.png"), zip, region=regions, runname=runname_ids, forcing=forcing),


### Import some useful python libraries
import os
from snakemake.io import Wildcards
from os.path import join
from itertools import product

curdir = os.getcwd()

include: 'snakefile_sfincs_build.smk'
include: 'snakefile_wflow.smk'
include: 'snakefile_sfincs_update.smk'
include: 'snakefile_fiat.smk'

if os.name == 'nt': #Running on windows
    root_dir = join("p:/",config['root_dir'])
elif os.name == "posix": #Running on linux
    root_dir = join("/p", config['root_dir'])
dir_runs = config['dir_runs']

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

rule all_general:
    input:
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "fiat","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "output", "spatial.fgb"), zip, region=region, runname=runname_ids, precip_forcing=precip_forcing, CF_rain=CF_rain, tidemodel=tidemodel, CF_SLR=CF_SLR, wind_forcing=wind_forcing, CF_wind=CF_wind),

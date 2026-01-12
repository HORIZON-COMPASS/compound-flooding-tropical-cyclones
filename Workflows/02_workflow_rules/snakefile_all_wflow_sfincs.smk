### Import some useful python libraries
import os
from snakemake.io import Wildcards
from os.path import join
from itertools import product

curdir = os.getcwd()

include: 'snakefile_sfincs_build.smk'
include: 'snakefile_wflow.smk'
include: 'snakefile_sfincs_update.smk'

if os.name == 'nt': #Running on windows
    root_dir = join("p:/",config['root_dir'])
elif os.name == "posix": #Running on linux
    root_dir = join("/p", config['root_dir'])
dir_runs = config['dir_runs']

run_combinations = []

for key, value in config['runname_ids'].items():
    for landuse in value['CF_landuse']:
        for rain, slr, wind in product(
            value['CF_value_rain'],
            value['CF_value_SLR'],
            value['CF_value_wind']
        ):
            run_combinations.append({
                "region": value["region"],
                "runname": key,
                "precip_forcing": value["precip_forcing"],
                "tidemodel": value["tidemodel"],
                "wind_forcing": value["wind_forcing"],
                "CF_rain": rain,
                "CF_SLR": slr,
                "CF_wind": wind,
                "landuse": landuse,
            })

region         = [c["region"] for c in run_combinations]
runname_ids    = [c["runname"] for c in run_combinations]
precip_forcing = [c["precip_forcing"] for c in run_combinations]
tidemodel      = [c["tidemodel"] for c in run_combinations]
wind_forcing   = [c["wind_forcing"] for c in run_combinations]
CF_rain        = [c["CF_rain"] for c in run_combinations]
CF_SLR         = [c["CF_SLR"] for c in run_combinations]
CF_wind        = [c["CF_wind"] for c in run_combinations]
landuse        = [c["landuse"] for c in run_combinations]

wildcard_constraints:
    precip_forcing='|'.join(map(re.escape, set(precip_forcing))),
    wind_forcing='|'.join(map(re.escape, set(wind_forcing))),
    CF_SLR=r"-?\d*\.?\d+",
    CF_wind=r"-?\d*\.?\d+",
    CF_rain=r"-?\d*\.?\d+",

rule all_general:
    input:
        expand(
            join(root_dir, dir_runs, "{region}", "{runname}", "sfincs", "event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}_{landuse}", "plot_output", "sfincs_basemap.png"),
            zip, region=region, runname=runname_ids, precip_forcing=precip_forcing, CF_rain=CF_rain, tidemodel=tidemodel, CF_SLR=CF_SLR, wind_forcing=wind_forcing, CF_wind=CF_wind, landuse=landuse)

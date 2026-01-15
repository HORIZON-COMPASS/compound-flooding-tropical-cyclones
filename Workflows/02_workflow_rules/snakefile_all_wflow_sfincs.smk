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

# Generate all combinations of runs based on config parameters
run_combinations = []

for runname, value in config['runname_ids'].items():
    region         = value['region']
    precip_forcing = value['precip_forcing']
    tidemodel      = value['tidemodel']
    wind_forcing   = value['wind_forcing']

    # Extract CF values from config
    rain_low,  rain_high  = value['CF_rain_uncert']
    rain_med              = value['CF_value_rain'][1]

    wind_low,  wind_high  = value['CF_wind_uncert']
    wind_med              = value['CF_value_wind'][1]

    slr_low,   slr_high   = value['CF_SLR_uncert']
    slr_med               = value['CF_value_SLR'][1]

    # Factual
    run_combinations.append(
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, 0, 0)
    )

    # All drivers (low / medium / high)
    run_combinations.extend([
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_low,  slr_low,  wind_low),
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_med,  slr_med,  wind_med),
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_high, slr_high, wind_high),
    ])

    # Rain only (low / medium / high)
    run_combinations.extend([
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_low,  0, 0),
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_med,  0, 0),
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_high, 0, 0),
    ])

    # Wind only (low / medium / high)
    run_combinations.extend([
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, 0, wind_low),
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, 0, wind_med),
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, 0, wind_high),
    ])

    # SLR only (low / medium / high)
    run_combinations.extend([
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, slr_low,  0),
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, slr_med,  0),
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, slr_high, 0),
    ])

    # Wind + SLR (low / medium / high)
    run_combinations.extend([
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, slr_low,  wind_low),
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, slr_med,  wind_med),
        (region, runname, precip_forcing, tidemodel, wind_forcing, 0, slr_high, wind_high),
    ])

    # Medium mixed rain combinations only
    run_combinations.extend([
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_med, 0, wind_med),
        (region, runname, precip_forcing, tidemodel, wind_forcing, rain_med, slr_med,  0),
    ])

# Unpack into wildcard lists
region, runname_ids, precip_forcing, tidemodel, wind_forcing, CF_rain, CF_SLR, CF_wind = map(
    list, zip(*run_combinations)
)

wildcard_constraints:
    precip_forcing='|'.join([re.escape(x) for x in set(precip_forcing)]),
    wind_forcing='|'.join([re.escape(x) for x in set(wind_forcing)]),
    CF_SLR=r"-?\d*\.?\d+",
    CF_wind=r"-?\d*\.?\d+",
    CF_rain=r"-?\d*\.?\d+",

rule all_general:
    input:
        expand(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_tp_{precip_forcing}_CF{CF_rain}_{tidemodel}_CF{CF_SLR}_{wind_forcing}_CF{CF_wind}", "plot_output", "sfincs_basemap.png"), zip, region=region, runname=runname_ids, precip_forcing=precip_forcing, CF_rain=CF_rain, tidemodel=tidemodel, CF_SLR=CF_SLR, wind_forcing=wind_forcing, CF_wind=CF_wind),

#!/bin/bash
#SBATCH --job-name=clim_download
#SBATCH --output=/gpfs/home4/dvertegaal1/logs/clim_%j.out
#SBATCH --time=24:00:00
#SBATCH --partition=staging
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

experiment="cont"

cd ~/git_repos/compound-flooding-tropical-cyclones/Workflows/04_scripts/preprocessing/ClimateDT/get_and_process_data || exit 1

pixi run -e compass-climatedt python -u get_storyline_data_Gen2_simple_levels.py "$experiment"
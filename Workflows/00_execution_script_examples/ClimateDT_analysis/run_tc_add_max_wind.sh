#!/bin/bash
#SBATCH --job-name=tc_wind_strategy
#SBATCH --output=/gpfs/home4/dvertegaal1/logs/tc_wind_strategy_%j.out
#SBATCH --time=06:00:00
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

set -eo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd /gpfs/home4/dvertegaal1/git_repos/compound-flooding-tropical-cyclones

pixi run -e compass-climatedt python -u \
    Workflows/04_scripts/preprocessing/ClimateDT/copilots/add_wind_to_strategy_tracks.py \
    --strategy all_ibtracs \
    --workers 5

echo "Finished at $(date)"
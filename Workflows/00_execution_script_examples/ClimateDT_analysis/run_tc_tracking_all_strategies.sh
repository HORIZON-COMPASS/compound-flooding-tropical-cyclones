#!/bin/bash
#SBATCH --job-name=tc_tracking_strategies
#SBATCH --output=/gpfs/home4/dvertegaal1/logs/tc_tracking_strategies_%j.out
#SBATCH --time=01:00:00
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
    Workflows/04_scripts/preprocessing/ClimateDT/copilots/run_tracker_strategies.py

echo "Finished at $(date)"
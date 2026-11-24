#!/bin/bash
#SBATCH --job-name=tc_tracking
#SBATCH --output=/gpfs/home4/dvertegaal1/logs/tc_tracking_%j.out
#SBATCH --time=02:00:00
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

set -eo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

SCEN=$1
REALIZATION=$2

cd /gpfs/home4/dvertegaal1/git_repos/compound-flooding-tropical-cyclones

echo "Scenario: $SCEN"
echo "Realization: $REALIZATION"

pixi run -e compass-climatedt python -u \
    Workflows/04_scripts/preprocessing/ClimateDT/tracking_TCs_SI.py \
    --scen "$SCEN" \
    --realization "$REALIZATION"

echo "Finished at $(date)"
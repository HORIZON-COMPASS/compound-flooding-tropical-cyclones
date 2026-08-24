#!/bin/bash
#SBATCH --job-name=pyjob            # give the job a name
#SBATCH --output=00_execution_script_examples/logs/slurm/pyjob_%j.out       # STDOUT → this file (%j is job ID)
#SBATCH --time=00:30:00             # hh:mm:ss wall‑time
#SBATCH --partition=test          # or your queue name
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

echo "Starting Python script on $(hostname) at $(date)"

module load pixi
# module load python/3.10

# Installing pixi environment
pixi install --environment compass-climatedt
pixi shell-hook --environment compass-climatedt > hook.sh
source hook.sh

# Run your script
echo "Running Python script..."
echo "Current directory: $(pwd)"

cd 04_scripts/preprocessing/ClimateDT

python -u get_storyline_data_Gen2.py

echo "Finished at $(date)"

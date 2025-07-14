#!/bin/bash
#SBATCH --job-name=pyjob            # give the job a name
#SBATCH --output=00_execution_script_examples/logs/slurm/pyjob_%j.out       # STDOUT → this file (%j is job ID)
#SBATCH --time=00:30:00             # hh:mm:ss wall‑time
#SBATCH --partition=test          # or your queue name
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --exclusive

set -euo pipefail

echo "Starting Python script on $(hostname) at $(date)"

module load pixi
# module load python/3.10

# Installing pixi environment
pixi install --environment compass-wflow
pixi shell-hook --environment compass-wflow > hook.sh
source hook.sh

# Run your script
cd ../Attribution_results/scripts

python plot_cf_pop_changes.py

echo "Finished at $(date)"

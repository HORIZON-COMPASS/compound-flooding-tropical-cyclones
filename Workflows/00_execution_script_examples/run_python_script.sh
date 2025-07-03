#!/bin/bash
#
#SBATCH --job-name=pyjob            # give the job a name
#SBATCH --output=pyjob_%j.out       # STDOUT → this file (%j is job ID)
#SBATCH --error=pyjob_%j.err        # STDERR → this file
#SBATCH --time=00:30:00             # hh:mm:ss wall‑time
#SBATCH --partition=test          # or your queue name
#SBATCH --nodes=1
#SBATCH --exclusive

set -euo pipefail

echo "Starting Python script on $(hostname) at $(date)"

module load pixi
module load python/3.10

# Installing pixi environment
pixi install --environment compass-wflow
pixi shell-hook --environment compass-wflow > hook.sh
source hook.sh

# Run your script
python ../Attribution_results/scripts/wflow_results.py

echo "Finished at $(date)"

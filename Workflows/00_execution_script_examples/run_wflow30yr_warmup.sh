#!/bin/bash
#SBATCH --job-name=compass-wflow-30yr         # Job name
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_wflow_%j_30yr.log     # Standard output and error log
#SBATCH --time=0-8:00:00           # Job duration (hh:mm:ss)
#SBATCH --partition 16vcpu
#SBATCH --exclusive 
#SBATCH --ntasks=1                  # Number of tasks (analyses) to run

set -e

module load pixi
module load julia

echo "Starting job on $(hostname) at $(date)"

#Going to the folder where git checkout is
ROOT="/u/morenodu/git_repos/compound-flooding-tropical-cyclones/cd "${ROOT}"

# Installing pixi environment
echo "Setting up pixi environment"
pixi install --environment compass-wflow
pixi shell-hook --environment compass-wflow > hook.sh
source hook.sh

# Install Julia environment
echo "Installing Julia packages (may be skipped if already installed)"
julia +1.9 -e 'using Pkg; Pkg.instantiate(); Pkg.add("Wflow")'

# Build the path to the input toml file
INPUT_TOML="/p/11210471-001-compass/03_Runs/sofala/Idai/wflow/event_precip_era5_hourly_zarr_CF0_30yr_copy/warmup/wflow_sbm.toml"
LOG_FILE="/p/11210471-001-compass/03_Runs/sofala/Idai/wflow/event_precip_era5_hourly_zarr_CF0_30yr_copy/warmup/run_default/log.txt"
# Run the Wflow CLI (this corresponds to the Snakemake shell command)
echo "Running Wflow"
julia +1.9 --threads 4 --project=~/.julia/environments/v1.9 -e "using Wflow; Wflow.run()" "${INPUT_TOML}" 2>&1 | tee "${LOG_FILE}"

echo "Job finished at $(date)"

exit 0

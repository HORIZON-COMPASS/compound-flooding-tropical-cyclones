#!/bin/bash
#SBATCH --job-name=compass-sfincs                                                      # Job name
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_wflow_sfincs_%j.log     # Standard output and error log
#SBATCH --time=0-3:00:00                                                               # Job duration (hh:mm:ss)
#SBATCH --partition 16vcpu
#SBATCH --exclusive 
#SBATCH --ntasks=1                                                                     # Number of tasks (analyses) to run

echo "=== Job started on $(date) ==="

# ────────────────────────────────────────────────
# Set temp directories for Pixi to avoid slow I/O
# ────────────────────────────────────────────────
export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
mkdir -p "$PIXI_CACHE_DIR"

echo "Loading modules..."
module load pixi
module load julia
module load apptainer

# Navigate to the repo directory
ROOT="/u/vertegaa/git_repos/COMPASS"
echo "Changing to ROOT directory: $ROOT"
cd "${ROOT}" || { echo "Failed to cd to ROOT directory!"; exit 1; }

# Install pixi environment if not already installed
echo "Installing pixi environment..."
pixi install --environment compass-wflow

# Activate pixi environment
echo "Activating pixi shell environment..."
eval "$(pixi shell-hook --environment compass-wflow)"

# Install Julia packages
echo "Setting up Julia environment..."
julia +1.9 -e 'using Pkg; Pkg.instantiate(); Pkg.add("Wflow")'

# Navigate to Snakemake workflow directory
echo "Changing to workflows directory..."
cd Workflows/02_workflow_rules || { echo "Failed to cd to workflow directory!"; exit 1; }

# Unlock Snakemake directory
echo "Unlocking Snakemake directory..."
snakemake --unlock -s snakefile_all_wflow_sfincs.smk --configfile ../01_config_snakemake/config_general_MZB_CF.yml 

# Dry-run of the workflow
echo "Running Snakemake dry-run..."
snakemake -s snakefile_all_wflow_sfincs.smk \
  --configfile ../01_config_snakemake/config_general_MZB_CF.yml \
  --cores 'all' --latency-wait 60 --wait-for-files --forceall

echo "=== Job finished on $(date) ==="
exit

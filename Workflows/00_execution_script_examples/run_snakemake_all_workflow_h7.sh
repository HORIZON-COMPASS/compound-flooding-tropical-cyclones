#!/bin/bash
#SBATCH --job-name=compass-sfincs                                                      # Job name
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_wflow_sfincs_%j.log     # Standard output and error log
#SBATCH --time=0-12:00:00                                                              # Job duration (hh:mm:ss)
#SBATCH --partition 16vcpu                                                             # Partition to use (e.g., 1vcpu, 4vcpu, test, etc.)
#SBATCH --exclusive                                                                    # Request exclusive access to the node
#SBATCH --ntasks=1                                                                     # Number of tasks (analyses) to run

# Make sure to adapt configuration for the specific use case (e.g., Idai, Freddy, Kenneth, etc.) by modifying the config_general_*.yml file accordingly.
echo "=== Job started on $(date) ==="

# Set temp directories for Pixi to avoid slow I/O
export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
mkdir -p "$PIXI_CACHE_DIR"

# /p is NFS, where the default HDF5 and SQLite file locking is unreliable. Without these,
# NetCDF reads fail with "HDF error" and GeoPackage reads with "database is locked".
export HDF5_USE_FILE_LOCKING=FALSE
export SQLITE_USE_OGR_VFS=YES

echo "Loading modules..."
module load pixi
module load julia
module load apptainer

# Navigate to the repo directory.
# Resolved automatically so this script works for any checkout; override by exporting
# COMPASS_ROOT before submitting.
ROOT="${COMPASS_ROOT:-$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel 2>/dev/null)}"
if [ -z "$ROOT" ]; then
    echo "Could not determine repo root; export COMPASS_ROOT=/path/to/compound-flooding-tropical-cyclones"; exit 1
fi
echo "Changing to ROOT directory: $ROOT"
cd "${ROOT}" || { echo "Failed to cd to ROOT directory!"; exit 1; }

# Install pixi environment if not already installed.
# compass-v1 is the HydroMT v1 stack (hydromt 1.4 / hydromt_sfincs 2.0 / hydromt_wflow 1.0).
echo "Installing pixi environment..."
pixi install --environment compass-v1

# Activate pixi environment
echo "Activating pixi shell environment..."
eval "$(pixi shell-hook --environment compass-v1)"

# Wflow.jl 1.0.2 lives in ~/.julia/environments/wflow1 and is invoked by the run rules with
# the juliaup default Julia. Nothing to install here; just make juliaup visible.
export PATH="$HOME/.juliaup/bin:$PATH"

# Navigate to Snakemake workflow directory
echo "Changing to workflows directory..."
cd Workflows/02_workflow_rules || { echo "Failed to cd to workflow directory!"; exit 1; }

# Set CONFIG to the use case you want to run.
CONFIG=../01_config_snakemake/config_general_mzb_Idai.yml

# Unlock Snakemake directory
echo "Unlocking Snakemake directory..."
snakemake --unlock -s snakefile_all_wflow_sfincs.smk --configfile "$CONFIG"

# Running of the workflow
echo "Running Snakemake..."
snakemake -s snakefile_all_wflow_sfincs.smk --configfile "$CONFIG" --cores 'all' --latency-wait 60 --wait-for-files --rerun-incomplete

echo "=== Job finished on $(date) ==="
exit
#!/bin/bash
#SBATCH --job-name=compass-sfincs          # Job name
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_sfincs_%j.log     # Standard output and error log
#SBATCH --time=0-0:30:00           # Job duration (hh:mm:ss)
#SBATCH --partition test
#SBATCH --exclusive 
#SBATCH --ntasks=1                  # Number of tasks (analyses) to run

# ────────────────────────────────────────────────
# Set temp directories for Pixi to avoid slow I/O
# ────────────────────────────────────────────────
export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
mkdir -p "$PIXI_CACHE_DIR"

echo "=== Loading modules ==="
module load pixi
module load apptainer 

echo "=== Changing directory to COMPASS root ==="
# Going to the folder where git checkout is
ROOT="/u/vertegaa/git_repos/COMPASS"
cd "${ROOT}"

# Installing pixi environment
echo "=== Setting up Pixi environment: compass-snake-sfincs ==="
pixi install --environment compass-snake-sfincs
eval "$(pixi shell-hook --environment compass-snake-sfincs)"

# Navigate to directory where the scripts are
echo "=== Changing to workflow rules directory ==="
cd Workflows/02_workflow_rules

 #Unlocking the directory for snakemake
echo "=== Unlocking Snakemake working directory ==="
snakemake --unlock -s snakefile_sfincs_update.smk --configfile ../01_config_snakemake/config_general_MZB.yml 

# running workflow with snakemake
echo "=== Performing dry run of Snakemake ==="
snakemake -s snakefile_sfincs_update.smk --configfile ../01_config_snakemake/config_general_MZB.yml --cores 'all' --latency-wait 60 --wait-for-files --rerun-incomplete

exit

#!/bin/bash
#SBATCH --job-name=sfincs_durban_precip_dis
#SBATCH --output=output_log_durban_%j.log      # Standard output and error log
#SBATCH --time=0-2:00:00                       # Job duration (2 hours should be sufficient for precip-only)
#SBATCH --partition=4vcpu                      # Partition
#SBATCH --exclusive
#SBATCH --ntasks=1                             # Number of tasks

# ==============================================================================
# SFINCS Workflow: Precipitation + Discharge Forcing (Durban April 2022)
# ==============================================================================
# This script runs SFINCS with:
# - ERA5 hourly precipitation
# - GloFAS v4.0 daily discharge forcing
# - No coastal forcing, no wind forcing
# ==============================================================================

echo "=========================================="
echo "SLURM_JOB_ID = $SLURM_JOB_ID"
echo "SLURM_NODELIST = $SLURM_NODELIST"
echo "=========================================="

# Activate the pixi environment
module load pixi
eval "$(pixi shell-hook -e compass-snake-sfincs)"

# Navigate to the workflow rules directory
cd /u/morenodu/git_repos/compound-flooding-tropical-cyclones/Workflows/02_workflow_rules

# Unlock the working directory if locked from a previous run
snakemake -s snakefile_sfincs_precip_discharge.smk \
    --configfile ../01_config_snakemake/config_durban_floods_2022_discharge.yml \
    --unlock

# Run the workflow
# --forceall: Force rerun of all rules
# --cores 4: Use 4 CPU cores
# --keep-going: Continue with independent jobs if one fails
snakemake -s snakefile_sfincs_precip_discharge.smk \
    --configfile ../01_config_snakemake/config_durban_floods_2022_discharge.yml \
    --cores 4 \
    --forceall \
    --keep-going

echo "=========================================="
echo "Workflow completed"
echo "=========================================="
#!/bin/bash
# Somerset quick SFINCS+Wflow test (CF0, 1-day event) — validates wflow/Julia + discharge coupling.
# Submit with: sbatch run_snakemake_all_workflow_somerset_wflow_quicktest.sh
#SBATCH --job-name=somerset-wqt
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_somerset_wqt_%j.log
#SBATCH --time=0-1:00:00
#SBATCH --partition=4vcpu
#SBATCH --exclusive
#SBATCH --ntasks=1

export HDF5_USE_FILE_LOCKING=FALSE
module load pixi
# Wflow runs use juliaup 1.11.3 (set on PATH inside the snakefile run rules) to match the wflow1 env.

ROOT="/u/morenodu/git_repos/compound-flooding-tropical-cyclones"
cd "${ROOT}"
eval "$(pixi shell-hook -e compass-v1)"

cd Workflows/02_workflow_rules
CONFIG=../01_config_snakemake/config_general_somerset_wflow_quicktest.yml
SMK=snakefile_all_wflow_sfincs.smk

snakemake --unlock -s $SMK --configfile $CONFIG
snakemake -s $SMK --configfile $CONFIG --cores 'all' --latency-wait 120 --wait-for-files --forceall

exit

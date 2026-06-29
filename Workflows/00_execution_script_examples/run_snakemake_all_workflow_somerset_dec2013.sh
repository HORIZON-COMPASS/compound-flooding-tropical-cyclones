#!/bin/bash
# Somerset Dec-2013 compound test (CF0): GTSM archive coastal + CEH precip + Wflow + SFINCS (v1).
# Submit with: sbatch run_snakemake_all_workflow_somerset_dec2013.sh
#SBATCH --job-name=somerset-dec
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_somerset_dec_%j.log
#SBATCH --time=0-1:00:00
#SBATCH --partition=4vcpu
#SBATCH --exclusive
#SBATCH --ntasks=1

export HDF5_USE_FILE_LOCKING=FALSE
module load pixi
ROOT="/u/morenodu/git_repos/compound-flooding-tropical-cyclones"
cd "${ROOT}"
eval "$(pixi shell-hook -e compass-v1)"
cd Workflows/02_workflow_rules
CONFIG=../01_config_snakemake/config_general_somerset_dec2013.yml
SMK=snakefile_all_wflow_sfincs.smk
snakemake --unlock -s $SMK --configfile $CONFIG
snakemake -s $SMK --configfile $CONFIG --cores 'all' --latency-wait 120 --wait-for-files --forceall
exit

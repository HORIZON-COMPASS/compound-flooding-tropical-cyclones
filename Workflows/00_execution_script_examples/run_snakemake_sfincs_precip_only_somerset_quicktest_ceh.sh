#!/bin/bash
# Somerset quick end-to-end test with local CEH-GEAR forcing (CF0, 3-hour Dec-2013 window).
# Submit with: sbatch run_snakemake_sfincs_precip_only_somerset_quicktest_ceh.sh
#SBATCH --job-name=somerset-qt-ceh
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_somerset_qt_ceh_%j.log
#SBATCH --time=0-0:40:00
#SBATCH --partition=4vcpu
#SBATCH --exclusive
#SBATCH --ntasks=1

export HDF5_USE_FILE_LOCKING=FALSE
module load pixi

ROOT="/u/morenodu/git_repos/compound-flooding-tropical-cyclones"
cd "${ROOT}"
eval "$(pixi shell-hook -e compass-v1)"

cd Workflows/02_workflow_rules
CONFIG=../01_config_snakemake/config_general_somerset_quicktest_ceh.yml
SMK=snakefile_sfincs_precip_only.smk

snakemake --unlock -s $SMK --configfile $CONFIG
snakemake -s $SMK --configfile $CONFIG --cores 'all' --latency-wait 120 --wait-for-files --forceall

exit

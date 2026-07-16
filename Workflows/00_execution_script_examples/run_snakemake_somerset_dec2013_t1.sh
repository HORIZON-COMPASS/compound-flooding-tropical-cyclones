#!/bin/bash
# Somerset Dec-2013 — T1 counterfactual (past/pre-industrial rainfall). Full Wflow + SFINCS pipeline.
# Submit with: sbatch run_snakemake_somerset_dec2013_t1.sh
#SBATCH --job-name=somerset-t1
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_somerset_t1_%j.log
#SBATCH --time=0-20:00:00
#SBATCH --partition=4vcpu
#SBATCH --exclusive
#SBATCH --ntasks=1

export HDF5_USE_FILE_LOCKING=FALSE
module load pixi          # Deltares HPC — replace with your env manager if needed

# Resolve repo root from script location (Workflows/00_execution_script_examples/)
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"
eval "$(pixi shell-hook -e compass-v1)"
cd Workflows/02_workflow_rules
CONFIG=../01_config_snakemake/config_general_somerset_dec2013_t1.yml
SMK=snakefile_all_wflow_sfincs.smk
snakemake --unlock -s $SMK --configfile $CONFIG
snakemake -s $SMK --configfile $CONFIG --cores 'all' --latency-wait 120 --wait-for-files --forceall
exit

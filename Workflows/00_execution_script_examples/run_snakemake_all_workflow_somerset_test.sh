#!/bin/bash
# Somerset Levels — full HydroMT v1 pipeline (SFINCS + Wflow), CF0 baseline, _test runname.
# Submit with: sbatch run_snakemake_all_workflow_somerset_test.sh
#SBATCH --job-name=somerset-test-v1
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_somerset_test_%j.log
#SBATCH --time=0-4:00:00
#SBATCH --partition 4vcpu
#SBATCH --exclusive
#SBATCH --ntasks=1

export HDF5_USE_FILE_LOCKING=FALSE

module load pixi
module load julia

ROOT="/u/morenodu/git_repos/compound-flooding-tropical-cyclones"
cd "${ROOT}"

# Activate the HydroMT v1 environment (hydromt 1.4 / hydromt_sfincs 2.0.0rc3 / hydromt_wflow 1.0.2)
eval "$(pixi shell-hook -e compass-v1)"

# Wflow.jl 1.0.2 is already installed in ~/.julia/environments/wflow1 (used by the run rules).

cd Workflows/02_workflow_rules

CONFIG=../01_config_snakemake/config_general_somerset_test.yml
SMK=snakefile_all_wflow_sfincs.smk

# Unlock + DAG + run
snakemake --unlock -s $SMK --configfile $CONFIG
snakemake -s $SMK --configfile $CONFIG --forceall --rulegraph | dot -Tpng > dag_somerset_test.png
snakemake -s $SMK --configfile $CONFIG --cores 'all' --latency-wait 180 --wait-for-files --forceall

exit

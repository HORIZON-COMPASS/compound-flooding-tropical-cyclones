#!/bin/bash
#SBATCH --job-name=compass-fiat-durban                                         # Job name
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_fiat_%j.log    # Standard output and error log
#SBATCH --time=0-2:00:00                                                        # Job duration (hh:mm:ss)
#SBATCH --partition 4vcpu
#SBATCH --exclusive
#SBATCH --ntasks=1                                                              # Number of tasks to run

module load pixi

ROOT="/u/morenodu/git_repos/compound-flooding-tropical-cyclones/"
cd "${ROOT}"

# Install pixi environment and required packages
pixi install --environment compass-fiat
pixi run --environment compass-fiat pip install "hydromt_fiat==0.5.10" "dask==2024.12.1"
pixi shell-hook --environment compass-fiat > hook.sh
source hook.sh

# Make sure the correct version of libstdc++ is being loaded
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH set to: $LD_LIBRARY_PATH"

cd Workflows/02_workflow_rules

# Unlock snakemake working directory
snakemake --unlock -s snakefile_fiat_precip_only.smk \
    --configfile ../01_config_snakemake/config_fiat_durban_floods_2022.yml

# Generate workflow DAG visualization
snakemake -s snakefile_fiat_precip_only.smk \
    --configfile ../01_config_snakemake/config_fiat_durban_floods_2022.yml \
    --forceall --rulegraph | dot -Tpng > dag_smk_fiat_durban.png

# Dry run to validate before submission
snakemake -n -s snakefile_fiat_precip_only.smk \
    --configfile ../01_config_snakemake/config_fiat_durban_floods_2022.yml

# Run the workflow
snakemake -s snakefile_fiat_precip_only.smk \
    --configfile ../01_config_snakemake/config_fiat_durban_floods_2022.yml \
    --cores all --latency-wait 60 --wait-for-files

exit

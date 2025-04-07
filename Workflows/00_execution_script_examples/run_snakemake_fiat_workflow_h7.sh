#!/bin/bash
#SBATCH --job-name=compass-fiat                                              # Job name
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_fiat_%j.log     # Standard output and error log
#SBATCH --time=0-2:00:00                                                       # Job duration (hh:mm:ss)
#SBATCH --partition 4vcpu
#SBATCH --exclusive 
#SBATCH --ntasks=1                                                             # Number of tasks (analyses) to run

module load pixi

#Going to the folder where git checkout is
#ROOT="/u/couasnon/git_repos/COMPASS/COMPASS"
#ROOT="/u/bovensch/git_repos/COMPASS"
ROOT="/u/vertegaa/git_repos/COMPASS"
# ROOT="/u/aleksand/compound-flooding-tropical-cyclones/"
cd "${ROOT}"

# Installing pixi environment
pixi install --environment compass-fiat
pixi run --environment compass-fiat pip install "hydromt_fiat @ git+https://github.com/Deltares/hydromt_fiat.git"
pixi run --environment compass-fiat conda install ibstdcxx-ng=12 gcc
pixi run --environment compass-fiat conda install --force-reinstall pandas xarray
pixi shell-hook --environment compass-fiat > hook.sh
source hook.sh

# Make sure the correct version of libstdc++ is being loaded
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH set to: $LD_LIBRARY_PATH"

# Navigate to directory where the scripts are
cd Workflows/02_workflow_rules

#Unlocking the directory for snakemake
snakemake --unlock -s snakefile_fiat.smk --configfile ../01_config_snakemake/config_general_MZB.yml 

# # running workflow with snakemake
snakemake -s snakefile_fiat.smk --configfile ../01_config_snakemake/config_general_MZB.yml --forceall --rulegraph | dot -Tpng > dag_smk_fiat.png
snakemake -n -s snakefile_fiat.smk --configfile ../01_config_snakemake/config_general_MZB.yml --cores 'all' --latency-wait 60 --wait-for-files

exit

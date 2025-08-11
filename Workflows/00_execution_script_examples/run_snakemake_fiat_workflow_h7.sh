#!/bin/bash
#SBATCH --job-name=compass-fiat                                              # Job name
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_fiat_%j.log     # Standard output and error log
#SBATCH --time=0-0:30:00                                                       # Job duration (hh:mm:ss)
#SBATCH --partition test
#SBATCH --exclusive 
#SBATCH --ntasks=1                                                             # Number of tasks (analyses) to run

export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
mkdir -p "$PIXI_CACHE_DIR"

# point to _all_ .so in the env's lib directory
export LD_LIBRARY_PATH="${ROOT}/.pixi/envs/compass-fiat/lib:${LD_LIBRARY_PATH:-}"

module load pixi

#Going to the folder where git checkout is
#ROOT="/u/couasnon/git_repos/COMPASS/COMPASS"
#ROOT="/u/bovensch/git_repos/COMPASS"
ROOT="/u/vertegaa/git_repos/COMPASS"
# ROOT="/u/aleksand/compound-flooding-tropical-cyclones/"
cd "${ROOT}"

# Installing pixi environment
# pixi install --environment compass-fiat
eval "$(pixi shell-hook --environment compass-fiat)"
pixi run --environment compass-fiat pip install "hydromt_fiat @ git+https://github.com/Deltares/hydromt_fiat.git"
# # Add required packages dynamically
# pixi add libstdcxx-ng=12
# pixi add gcc
# pixi reinstall pandas xarray

# Explicitly set LD_LIBRARY_PATH so it prefers env libs over system libs
export LD_LIBRARY_PATH="$ROOT/.pixi/envs/compass-fiat/lib:$LD_LIBRARY_PATH"

# Navigate to directory where the scripts are
cd Workflows/02_workflow_rules

#Unlocking the directory for snakemake
snakemake --unlock -s snakefile_fiat.smk --configfile ../01_config_snakemake/config_general_MZB.yml 

# # running workflow with snakemake
# snakemake -s snakefile_fiat.smk --configfile ../01_config_snakemake/config_general_MZB.yml --forceall --rulegraph | dot -Tpng > dag_smk_fiat.png
snakemake -s snakefile_fiat.smk --configfile ../01_config_snakemake/config_general_MZB.yml --cores 'all' --latency-wait 60 --wait-for-files --rerun-incomplete

exit

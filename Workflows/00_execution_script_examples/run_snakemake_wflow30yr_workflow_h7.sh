#!/bin/bash
#SBATCH --job-name=compass-wflow         # Job name
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_wflow_%j.log     # Standard output and error log
#SBATCH --time=0-4:00:00           # Job duration (hh:mm:ss)
#SBATCH --partition 16vcpu
#SBATCH --exclusive 
#SBATCH --ntasks=1                  # Number of tasks (analyses) to run

export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
mkdir -p "$PIXI_CACHE_DIR"

module load pixi
module load julia

#Going to the folder where git checkout is
#ROOT="/u/couasnon/git_repos/COMPASS/COMPASS"
#ROOT="/u/bovensch/git_repos/COMPASS"
ROOT="/u/vertegaa/git_repos/COMPASS"
# ROOT="/u/aleksand/compound-flooding-tropical-cyclones/"
cd "${ROOT}"

# Install pixi environment if not already installed
echo "Installing pixi environment..."
pixi install --environment compass-wflow

# Activate pixi environment
echo "Activating pixi shell environment..."
eval "$(pixi shell-hook --environment compass-wflow)"

# Install Julia environment
echo "Setting up Julia environment..."
julia +1.9 -e 'using Pkg; Pkg.instantiate(); Pkg.add("Wflow")'

# Navigate to directory where the scripts are
echo "Changing to workflows directory..."
cd Workflows/02_workflow_rules

#Unlocking the directory for snakemake
snakemake --unlock -s snakefile_wflow_30yr.smk --configfile ../01_config_snakemake/config_general_MZB.yml 

# running workflow with snakemake
# snakemake -s snakefile_wflow_30yr.smk --configfile ../01_config_snakemake/config_general_MZB.yml --forceall --rulegraph | dot -Tpdf > wflow_30yr.pdf
snakemake -s snakefile_wflow_30yr.smk --configfile ../01_config_snakemake/config_general_MZB.yml --cores 'all' --latency-wait 60 --wait-for-files --rerun-incomplete # --cores 4

exit

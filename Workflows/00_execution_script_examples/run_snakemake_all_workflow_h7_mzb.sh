#!/bin/bash
#SBATCH --job-name=compass-sfincs                                                      # Job name
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_wflow_sfincs_%j.log     # Standard output and error log
#SBATCH --time=0-2:00:00                                                               # Job duration (hh:mm:ss)
#SBATCH --partition 4vcpu
#SBATCH --exclusive 
#SBATCH --ntasks=1                                                                     # Number of tasks (analyses) to run

module load pixi
module load julia

#Going to the folder where git checkout is
#ROOT="/u/couasnon/git_repos/COMPASS/COMPASS"
#ROOT="/u/bovensch/git_repos/COMPASS"
ROOT="/u/vertegaa/git_repos/COMPASS"
# ROOT="/u/aleksand/compound-flooding-tropical-cyclones/"
cd "${ROOT}"

# Installing pixi environment
pixi install --environment compass-wflow
pixi shell-hook --environment compass-wflow > hook.sh
source hook.sh

# Install Julia environment
julia +1.9 -e 'using Pkg; Pkg.instantiate(); Pkg.add("Wflow")'


# Navigate to directory where the scripts are
cd Workflows/02_workflow_rules

#Unlocking the directory for snakemake
snakemake --unlock -s snakefile_all_wflow_sfincs.smk --configfile ../01_config_snakemake/config_general_MZB.yml 

# # running workflow with snakemake
snakemake -s snakefile_all_wflow_sfincs.smk --configfile ../01_config_snakemake/config_general_MZB.yml --forceall --rulegraph | dot -Tpng > dag_smk_all.png
snakemake -s snakefile_all_wflow_sfincs.smk --configfile ../01_config_snakemake/config_general_MZB.yml --cores 'all' --latency-wait 60 --wait-for-files

exit

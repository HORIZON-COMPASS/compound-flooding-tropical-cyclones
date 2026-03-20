#!/bin/bash
#SBATCH --job-name=compass-dfm         # Job name
#SBATCH --output=output_log_%j.log     # Standard output and error log
#SBATCH --time=0-03:00:00           # Job duration (hh:mm:ss)
#SBATCH --partition 16vcpu #test #16vcpu
#SBATCH --exclusive 
#SBATCH --ntasks=1                  # Number of tasks (analyses) to run

module load pixi

#Going to the folder where git checkout is
ROOT="/u/couasnon/git_repos/compound-flooding-tropical-cyclones/"
# ROOT="/u/bovensch/git_repos/COMPASS"
#ROOT="/u/vertegaa/git_repos/COMPASS"df 
#ROOT="/u/aleksand/compound-flooding-tropical-cyclones/"
cd "${ROOT}"

# Installing pixi environment
pixi install --environment compass-snake-dfm
pixi shell-hook --environment compass-snake-dfm > hook.sh
source hook.sh

# Navigate to directory where the scripts are
cd Workflows/02_workflow_rules

#Unlocking the directory for snakemake
snakemake --unlock -s snakefile_dfm.smk --configfile ../01_config_snakemake/config_general_somerset.yml 

# running workflow with snakemake
snakemake -s snakefile_dfm.smk --configfile ../01_config_snakemake/config_general_somerset.yml --forceall --rulegraph | dot -Tpdf > dag_dfm.pdf
snakemake -s snakefile_dfm.smk --configfile ../01_config_snakemake/config_general_somerset.yml --cores 'all' --latency-wait 60 --wait-for-files --keep-incomplete --rerun-triggers input # --forceall # --cores 4

exit

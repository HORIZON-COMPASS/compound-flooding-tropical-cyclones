#!/bin/bash
#SBATCH --job-name=compass-dfm         # Job name
#SBATCH --output=output_log_%j.log     # Standard output and error log
#SBATCH --time=0-0:30:00           # Job duration (hh:mm:ss)
#SBATCH --partition test
#SBATCH --exclusive 
#SBATCH --ntasks=1                  # Number of tasks (analyses) to run

module load pixi

#Going to the folder where git checkout is
#ROOT="/u/couasnon/git_repos/COMPASS/COMPASS"
# ROOT="/u/bovensch/git_repos/COMPASS"
ROOT="/u/vertegaa/git_repos/COMPASS"
cd "${ROOT}"

# Installing pixi environment
pixi install --environment compass-dfm
pixi shell-hook --environment compass-dfm > hook.sh
source hook.sh

# Navigate to directory where the scripts are
cd Workflows

#Unlocking the directory for snakemake
snakemake --unlock -s snakefile_dfm.smk --configfile config_snakemake/config_general_dfm.yml 

# running workflow with snakemake
snakemake -s snakefile_dfm.smk --configfile config_snakemake/config_general_dfm.yml --forceall --rulegraph | dot -Tpdf > dag.pdf
snakemake -s snakefile_dfm.smk --configfile config_snakemake/config_general_dfm.yml --cores 'all' --latency-wait 60 --wait-for-files  --forceall # --cores 4

exit

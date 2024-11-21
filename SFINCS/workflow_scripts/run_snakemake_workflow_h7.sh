#!/bin/bash
#SBATCH --job-name=compass-sfincs          # Job name
#SBATCH --output=output_log_%j.log     # Standard output and error log
#SBATCH --time=0-0:30:00 #0-2:00:00           # Job duration (hh:mm:ss)
#SBATCH --partition test #4vcpu
#SBATCH --exclusive 
#SBATCH --ntasks=1                  # Number of tasks (analyses) to run

module load pixi

#Going to the folder where git checkout is
ROOT="/u/couasnon/git_repos/COMPASS/COMPASS"
cd "${ROOT}"

# Installing pixi environment
pixi install --environment compass-snake-sfincs
pixi shell-hook --environment compass-snake-sfincs > hook.sh
source hook.sh

# Navigate to directory where the scripts are
cd SFINCS/workflow_scripts

#Unlocking the directory for snakemake
snakemake --unlock -s snakefile_dict --configfile config/config_general.yml 

# running workflow with snakemake
snakemake -s snakefile_dict --configfile config/config_general.yml --forceall --rulegraph | dot -Tpdf > dag.pdf
snakemake -n -s snakefile_dict --configfile config/config_general.yml --cores 4 --forcerun

exit

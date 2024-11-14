#!/bin/bash
#SBATCH --job-name=compass-sfincs          # Job name
#SBATCH --output=output_log.log     # Standard output and error log
#SBATCH --time=0-0:30:00 #0-2:00:00           # Job duration (hh:mm:ss)
#SBATCH --partition test #4vcpu
#SBATCH --ntasks=1                  # Number of tasks (analyses) to run
##SBATCH --mail-user=natalia.aleksandrova@deltares.nl
##SBATCH --mail-type=ALL

module load pixi

#Going to the folder where git checkout is
ROOT="/p/11210471-001-compass/02_Models/workflow_test/COMPASS"
cd "${ROOT}"

# Installing pixi environment
pixi install --environment compass-snake-sfincs
pixi shell-hook --environment compass-snake-sfincs > hook.sh
source hook.sh

# Navigate to directory where the scripts are
cd SFINCS/workflow_scripts

#Unlocking the directory for snakemake
snakemake --unlock -s snakefile --configfile config/config_Idai.yml 

# running workflow with snakemake
snakemake --configfile config/config_Idai.yml --forceall --rulegraph | dot -Tpdf > dag.pdf
snakemake -s snakefile --configfile config/config_Idai.yml --cores 4 --rerun-triggers mtime 

exit

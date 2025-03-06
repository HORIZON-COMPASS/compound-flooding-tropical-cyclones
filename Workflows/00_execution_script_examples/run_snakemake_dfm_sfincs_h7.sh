#!/bin/bash
#SBATCH --job-name=compass-sfincs         # Job name
#SBATCH --output=00_execution_script_examples/logs/slurm/output_log_%j.log     # Standard output and error log
#SBATCH --time=0-0:30:00           # Job duration (hh:mm:ss)
#SBATCH --partition test    # 16vcpu
#SBATCH --exclusive 
#SBATCH --ntasks=1                  # Number of tasks (analyses) to run

module load pixi
module load apptainer/1.2.5     # Load the Apptainer container system software.
module load intelmpi/2021.11.0   # Load the  message-passing library for parallel simulations.

#Going to the folder where git checkout is
#ROOT="/u/couasnon/git_repos/COMPASS/COMPASS"
#ROOT="/u/bovensch/git_repos/COMPASS"
ROOT="/u/vertegaa/git_repos/COMPASS"
# ROOT="/u/aleksand/compound-flooding-tropical-cyclones/"
cd "${ROOT}"

# Installing pixi environment
pixi install --environment compass-snake-dfm-sfincs
pixi shell-hook --environment compass-snake-dfm-sfincs > hook.sh
source hook.sh

# Navigate to directory where the scripts are
cd Workflows/02_workflow_rules

#Unlocking the directory for snakemake
snakemake --unlock -s 02_workflow_rules/snakefile_all_dfm_sfincs.smk --configfile ../01_config_snakemake/config_general.yml 

# running workflow with snakemake
snakemake -s snakefile_all_dfm_sfincs.smk.smk --configfile ../01_config_snakemake/config_general.yml --forceall --rulegraph | dot -Tpdf > dag.pdf
# snakemake -s snakefile_all_dfm_sfincs.smk.smk --configfile ../01_config_snakemake/config_general.yml --cores 'all' --latency-wait 60 --wait-for-files  --forceall # --cores 4
snakemake -s snakefile_all_dfm_sfincs.smk.smk --configfile ../01_config_snakemake/config_general.yml--cores 'all' --latency-wait 60 --wait-for-files--jobs 10 --executor cluster-generic --cluster-generic-submit-cmd "sbatch --job-name {resources.jobname} --time {resources.time} --partition {resources.partition} --ntasks-per-node={resources.taskspernode} --nodes=1 --parsable"

exit

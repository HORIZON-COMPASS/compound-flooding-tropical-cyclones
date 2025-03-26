#!/bin/bash
#SBATCH --job-name=compass-dfm         # Job name
#SBATCH --output=output_log_%j.log     # Standard output and error log
#SBATCH --time=0-03:00:00           # Job duration (hh:mm:ss)
#SBATCH --partition 1vcpu
#SBATCH --exclusive 
#SBATCH --ntasks=1                  # Number of tasks (analyses) to run

module load pixi
module load apptainer/1.2.5     # Load the Apptainer container system software.
module load intelmpi/2021.11.0   # Load the  message-passing library for parallel simulations.
 

#Going to the folder where git checkout is
ROOT="/u/couasnon/git_repos/compound-flooding-tropical-cyclones/"
# ROOT="/u/bovensch/git_repos/COMPASS"
#ROOT="/u/vertegaa/git_repos/COMPASS"
#ROOT="/u/aleksand/compound-flooding-tropical-cyclones/"
cd "${ROOT}"

# Installing pixi environment
pixi install --environment compass-snake-dfm
pixi shell-hook --environment compass-snake-dfm > hook.sh
source hook.sh

# Navigate to directory where the scripts are
cd Workflows/02_workflow_rules

#Unlocking the directory for snakemake
snakemake --unlock -s snakefile_dfm_cluster.smk --configfile ../01_config_snakemake/config_general_somerset.yml 

# Making the DAG
snakemake -s snakefile_dfm_cluster.smk --configfile ../01_config_snakemake/config_general_somerset.yml --rulegraph | dot -Tpng > dag_dfm_somerset2.png

# Snakemake using dedicated slurm functionality - does not work fully on h7 yet
#snakemake -s snakefile_dfm_cluster.smk --configfile config_snakemake/config_general.yml --jobs 10 --executor slurm --default-resources "slurm_account='hot'"  "slurm_partition='4vcpu'" "runtime=30" "tasks=1" --set-resources "run_dfm:slurm_partition='16vcpu'" "run_dfm:time=180"

# Snakemake using generic cluster functionality
snakemake -s snakefile_dfm_cluster.smk  --configfile ../01_config_snakemake/config_general_somerset.yml --jobs 10 --executor cluster-generic --cluster-generic-submit-cmd "sbatch --job-name {resources.jobname} --time {resources.time} --partition {resources.partition} --ntasks-per-node={resources.taskspernode} --nodes=1 --parsable"
exit

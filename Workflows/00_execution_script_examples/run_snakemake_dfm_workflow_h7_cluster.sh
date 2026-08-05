#!/bin/bash
#SBATCH --job-name=compass-dfm                                                # Job name
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_dfm_%j.log     # Standard output and error log
#SBATCH --time=0-3:00:00                                                      # Job duration (hh:mm:ss)
#SBATCH --partition 1vcpu                                                     # Partition to use (e.g., 1vcpu, 4vcpu, test, etc.)
#SBATCH --exclusive                                                           # Request exclusive access to the node
#SBATCH --ntasks=1                                                            # Number of tasks (analyses) to run

# Make sure to adapt configuration for the specific use case (e.g., Idai, Freddy, Kenneth, etc.) by modifying the config_general_*.yml file accordingly.
module load pixi
module load apptainer/1.2.5                                                   # Load the Apptainer container system software.
module load intelmpi/2021.11.0                                                # Load the  message-passing library for parallel simulations.
 

#Going to the folder where git checkout is
#ROOT="/u/couasnon/git_repos/COMPASS/COMPASS"
# ROOT="/u/bovensch/git_repos/COMPASS"
# ROOT="/u/morenodu/git_repos/compound-flooding-tropical-cyclones/
# ROOT="/u/aleksand/compound-flooding-tropical-cyclones/"
ROOT="/u/vertegaa/git_repos/COMPASS"
cd "${ROOT}"

# Installing pixi environment
pixi install --environment compass-snake-dfm
pixi shell-hook --environment compass-snake-dfm > hook.sh
source hook.sh

# Navigate to directory where the scripts are
cd Workflows/02_workflow_rules

#Unlocking the directory for snakemake
snakemake --unlock -s snakefile_dfm_cluster.smk --configfile ../01_config_snakemake/config_general_MZB.yml 

# running workflow with snakemake
snakemake -s snakefile_dfm_cluster.smk --configfile ../01_config_snakemake/config_general_MZB.yml --forceall --rulegraph | dot -Tpng > dag_dfm_cluster.png

# Snakemake using generic cluster functionality
snakemake -s snakefile_dfm_cluster.smk  --configfile ../01_config_snakemake/config_general_MZB.yml --jobs 10 --executor cluster-generic --scheduler greedy \
  --cluster-generic-submit-cmd "sbatch --job-name {resources.jobname} --time {resources.time} --partition {resources.partition} --ntasks-per-node={resources.taskspernode} --nodes=1 --parsable"

exit

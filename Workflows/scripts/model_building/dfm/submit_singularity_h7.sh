#! /bin/bash

# Usage: 
#   - Place this script in the same folder as the model
#   - Modify this script where needed (e.g. number of nodes, number of tasks per node, singularity version, model folder)
#   - Execute this script using:
#     sbatch ./submit_singularity_h7_mv.sh
#
# This is an h7 specific script for single- or multi-node simulations

#SBATCH --nodes=1               #-N, -n is total numer of nodes. $SLURM_NTASKS = "--nodes" times "--ntasks-per-node"
#SBATCH --ntasks-per-node=4     #you pay for a minimum of 1/4 of the cores on a node
#SBATCH --job-name=JOBNAME     #-J
#SBATCH --time 1-00:00:00        #-t, reduce the expected time if possible to increase your priority
#SBATCH --partition=4vcpu      #Partition name

#---You will need to modify the input below this line---

#Set MDU file
mduFile=MDUFILE

#---You do not need to modify anything below this line---

# Set the location of the Singularity container.
export singularityFolder=/p/11210471-001-compass/02_Models/00_executables/Apptainer_DFLOWFM_2023_02_copy/

#
#
# --- You shouldn't need to change the lines below ------------------------

# stop after an error occurred:
set -e
nPart=$SLURM_NTASKS*$SLURM_NNODES
echo ""
echo "nPart" $nPart

echo ""
echo "Partitioning..."
echo "SLURM_NTASKS: $SLURM_NTASKS"
srun -n 1 $singularityFolder/execute_singularity_h7.sh -p 0 run_dflowfm.sh --partition:ndomains=$SLURM_NTASKS:icgsolver=6 $mduFile

echo ""
echo "Simulation..."
srun $singularityFolder/execute_singularity_h7.sh -p 0 run_dflowfm.sh $mduFile
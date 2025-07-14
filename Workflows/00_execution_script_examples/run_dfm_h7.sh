#!/bin/bash
#SBATCH --job-name=compass-dfm          # Job name
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_dfm_%j.log     # Standard output and error log
#SBATCH --time=0-4:00:00           # Job duration (hh:mm:ss)
#SBATCH --partition 16vcpu
#SBATCH --nodes=1                    # Adjust as needed
#SBATCH --ntasks-per-node=16    
#SBATCH --exclusive 

module load apptainer/1.2.5     # Load the Apptainer container system software.
module load intelmpi/2021.11.0   # Load the  message-passing library for parallel simulations.
 
set -euo pipefail

echo "Running DFM model manually..."

EVENT_DIR="/p/11210471-001-compass/03_Runs/sofala/Idai/dfm/event_450_gebco2024_MZB_GTSMv41_CF0_spw_IBTrACS_CF0"
SCRIPT="${EVENT_DIR}/run_singularity_h7.sh"

# Optional: copy script and make it executable if needed
# cp "$SCRIPT" "${SCRIPT}_copy"
# chmod +x "${SCRIPT}_copy"
# "${SCRIPT}_copy"

chmod +x "$SCRIPT"
"$SCRIPT"

echo "Finished running DFM"

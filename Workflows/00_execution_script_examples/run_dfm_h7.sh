#!/bin/bash
#SBATCH --job-name=compass-dfm
#SBATCH --output=output_log_%j.log
#SBATCH --time=0-4:00:00
#SBATCH --partition=16vcpu
#SBATCH --nodes=1                    # Adjust as needed
#SBATCH --ntasks-per-node=16        # Match Snakemake's `taskspernode`
#SBATCH --exclusive

set -euo pipefail

echo "Running DFM model manually..."

EVENT_DIR="/p/11210471-001-compass/03_Runs/sofala/Idai/dfm/event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_zarr_CF0"
SCRIPT="${EVENT_DIR}/run_singularity_h7.sh"

# Optional: copy script and make it executable if needed
# cp "$SCRIPT" "${SCRIPT}_copy"
# chmod +x "${SCRIPT}_copy"
# "${SCRIPT}_copy"

chmod +x "$SCRIPT"
"$SCRIPT"

echo "Finished running DFM"

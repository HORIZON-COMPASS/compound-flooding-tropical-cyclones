#!/bin/bash
#SBATCH --job-name=compass-sfincs
#SBATCH --output=00_execution_script_examples/logs/slurm/output_sfincs_%j.log
#SBATCH --time=0-0:30:00
#SBATCH --partition=test
#SBATCH --nodes=1                    # Adjust as needed
#SBATCH --exclusive

# Run SFINCS using Docker
echo "Running SFINCS model using Docker..."

RUN_DIR="/p/11210471-001-compass/03_Runs/sofala/Idai/sfincs/event_tp_era5_hourly_CF0_GTSMv41_CF0_spw_IBTrACS_CF0_waves"

echo "Working directory: ${RUN_DIR}"

docker run --rm \
    --mount type=bind,src="${RUN_DIR}",target=/data \
    deltares/sfincs-cpu:latest sfincs

echo "Finished running SFINCS"

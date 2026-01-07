#!/bin/bash
#SBATCH --job-name=compass-sfincs
#SBATCH --output=00_execution_script_examples/logs/slurm/output_sfincs_%j.log
#SBATCH --time=0-2:00:00
#SBATCH --partition=16vcpu
#SBATCH --nodes=1                    # Adjust as needed
#SBATCH --exclusive

# Run SFINCS using Docker
echo "Running SFINCS model using Docker..."

RUN_DIR="/p/11210471-001-compass/03_Runs/sofala/Idai/sfincs/event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0"

echo "Working directory: ${RUN_DIR}"

docker run --rm \
    --mount type=bind,src="${RUN_DIR}",target=/data \
    deltares/sfincs-cpu:sfincs-v2.2.0-col-dEze-Release

echo "Finished running SFINCS"

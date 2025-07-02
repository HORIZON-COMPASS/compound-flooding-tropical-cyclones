#!/bin/bash
#SBATCH --job-name=compass-wflow
#SBATCH --output=00_execution_script_examples/logs/slurm/wflow_output_log_%j.log
#SBATCH --time=0-0:30:00
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --exclusive

set -euo pipefail

module load julia

echo "Running Wflow model manually..."

# Set run configuration
RUN_DIR="/p/11210471-001-compass/03_Runs/sofala/Idai/wflow/event_precip_era5_hourly_zarr_CF0_f_soilthick"

EXE="/p/11210471-001-compass/01_Models/00_executables/wflow0.8.1/wflow_cli/bin/wflow_cli.exe"
JULIA_ENV_FN="$HOME/.julia/environments/v1.9"

# 1. Warmup step
echo "Running warmup..."

TOML_WARMUP="${RUN_DIR}/warmup/wflow_sbm.toml"
if [[ -f "$EXE" ]]; then
    "$EXE" "$TOML_WARMUP" || julia +1.9 --threads 4 --project="$JULIA_ENV_FN" -e "using Wflow; Wflow.run()" "$TOML_WARMUP"
else
    julia +1.9 --threads 4 --project="$JULIA_ENV_FN" -e "using Wflow; Wflow.run()" "$TOML_WARMUP"
fi

# 2. Event step
echo "Running event..."

TOML_EVENT="${RUN_DIR}/events/wflow_sbm.toml"
if [[ -f "$EXE" ]]; then
    "$EXE" "$TOML_EVENT" || julia +1.9 --threads 4 --project="$JULIA_ENV_FN" -e "using Wflow; Wflow.run()" "$TOML_EVENT"
else
    julia +1.9 --threads 4 --project="$JULIA_ENV_FN" -e "using Wflow; Wflow.run()" "$TOML_EVENT"
fi

echo "Finished Wflow run."

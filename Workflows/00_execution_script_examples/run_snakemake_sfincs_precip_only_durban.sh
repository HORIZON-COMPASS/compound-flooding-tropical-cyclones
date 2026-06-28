#!/bin/bash
# Durban 2022 precipitation-only SFINCS workflow
# Forcing: ERA5 hourly precipitation (factual) + -10% counterfactual (climate attribution)
# Produces: event_precip_era5_hourly_CF0_no_wind  (factual)
#           event_precip_era5_hourly_CF-10_no_wind (counterfactual)
#SBATCH --job-name=durban-sfincs-precip-cf10  # Job name
#SBATCH --output=output_log_durban_%j.log      # Standard output and error log
#SBATCH --time=0-4:00:00                       # Job duration (4 hours for both factual + CF-10 runs)
#SBATCH --partition=4vcpu                      # Partition
#SBATCH --exclusive
#SBATCH --ntasks=1                             # Number of tasks

export HDF5_USE_FILE_LOCKING=FALSE

module load pixi

# Project root directory
ROOT="/u/morenodu/git_repos/compound-flooding-tropical-cyclones/"
cd "${ROOT}"

# Activate pixi environment (SFINCS only, no Wflow/Julia needed)
eval "$(pixi shell-hook -e compass-snake-sfincs)"

# Navigate to workflow directory
cd Workflows/02_workflow_rules

# Unlock directory for snakemake
snakemake --unlock \
    -s snakefile_sfincs_precip_only.smk \
    --configfile ../01_config_snakemake/config_durban_floods_2022.yml

# Generate workflow DAG (optional, for visualization)
snakemake \
    -s snakefile_sfincs_precip_only.smk \
    --configfile ../01_config_snakemake/config_durban_floods_2022.yml \
    --forceall \
    --rulegraph | dot -Tpng > dag_durban_precip_only.png

# Run the workflow
snakemake \
    -s snakefile_sfincs_precip_only.smk \
    --configfile ../01_config_snakemake/config_durban_floods_2022.yml \
    --cores 'all' \
    --latency-wait 180 \
    --wait-for-files \
    --forceall

exit

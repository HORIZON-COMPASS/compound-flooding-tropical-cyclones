#!/bin/bash
# v1 integration crash test: Somerset, the precip-only chain on the new land-use layout.
#SBATCH --job-name=v1int-durban
#SBATCH --output=/u/morenodu/git_repos/compass-v1-integrate/slurm_logs/v1int_durban_%j.log
#SBATCH --time=0-04:00:00
#SBATCH --partition=4vcpu
#SBATCH --ntasks=1
set -o pipefail
echo "=== Job started on $(date) on $(hostname) ==="
export HDF5_USE_FILE_LOCKING=FALSE
export SQLITE_USE_OGR_VFS=YES
export PATH="$HOME/.juliaup/bin:$PATH"      # v1 run rules use Wflow.jl 1.0.2 from ~/.julia/environments/wflow1
ROOT="/u/morenodu/git_repos/compass-v1-integrate"
ENVDIR="/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-v1"
export PATH="$ENVDIR/bin:$PATH"; export PROJ_DATA="$ENVDIR/share/proj"; export GDAL_DATA="$ENVDIR/share/gdal"
cd "$ROOT/Workflows/02_workflow_rules" || exit 1
python -c "import hydromt, hydromt_sfincs, hydromt_wflow; print('hydromt',hydromt.__version__,'| sfincs',hydromt_sfincs.__version__,'| wflow',hydromt_wflow.__version__)" 2>&1 | grep -v PROJ
snakemake --unlock -s snakefile_sfincs_precip_only.smk --configfile ../../smoketest/config_v1int_durban.yml >/dev/null 2>&1
echo "=== Running the precip-only chain (--keep-going) ==="
snakemake -s snakefile_sfincs_precip_only.smk --configfile ../../smoketest/config_v1int_durban.yml \
    --cores 4 --latency-wait 60 --wait-for-files --rerun-incomplete --keep-going --printshellcmds
echo "=== Snakemake exit code: $? ==="
echo "=== Finished on $(date) ==="

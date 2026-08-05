#!/bin/bash
#SBATCH --job-name=diag-smk
#SBATCH --output=/u/morenodu/git_repos/compass-v1-integrate/slurm_logs/diag_smk_%j.log
#SBATCH --time=0-00:25:00
#SBATCH --partition=4vcpu
#SBATCH --ntasks=1
export HDF5_USE_FILE_LOCKING=FALSE
export SQLITE_USE_OGR_VFS=YES
ENV=/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-v1
export PATH="$ENV/bin:$PATH" PROJ_DATA="$ENV/share/proj" GDAL_DATA="$ENV/share/gdal"
echo "=== node: $(hostname)  cores: $(nproc) ==="
R=/p/11210471-001-compass/03_Runs/somerset/SomersetLevels_v1int/wflow_vito/event_precip_ceh_gear_CF0/warmup
rm -rf "$R"
cd /u/morenodu/git_repos/compass-v1-integrate/Workflows/02_workflow_rules
rm -rf .snakemake/locks
echo "########## snakemake, single rule, on a COMPUTE node ##########"
timeout 600 snakemake -s snakefile_all_wflow_sfincs.smk \
  --configfile ../../smoketest/config_v1int_somerset.yml \
  --until update_forcing_wflow_warmup --cores 4 --latency-wait 60 --rerun-incomplete
echo "  exit=$? (124 = HUNG -> reproduced)"

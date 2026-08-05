#!/bin/bash
#SBATCH --job-name=diag-warmup
#SBATCH --output=/u/morenodu/git_repos/compass-v1-integrate/slurm_logs/diag_warmup_%j.log
#SBATCH --time=0-00:20:00
#SBATCH --partition=4vcpu
#SBATCH --ntasks=1
export HDF5_USE_FILE_LOCKING=FALSE
export SQLITE_USE_OGR_VFS=YES
ENV=/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-v1
export PATH="$ENV/bin:$PATH" PROJ_DATA="$ENV/share/proj" GDAL_DATA="$ENV/share/gdal"
cd /u/morenodu/git_repos/compass-v1-integrate
echo "=== node: $(hostname)  cores visible: $(nproc) ==="
rm -rf /p/11210471-001-compass/03_Runs/somerset/DIAG_warmup
echo "########## threaded (as the rule runs) ##########"
timeout 300 python smoketest/diag_warmup_hang.py threaded; echo "  exit=$? (124=HUNG)"
rm -rf /p/11210471-001-compass/03_Runs/somerset/DIAG_warmup
echo "########## synchronous ##########"
timeout 300 python smoketest/diag_warmup_hang.py synchronous; echo "  exit=$? (124=HUNG)"

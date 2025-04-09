#!/bin/bash
#SBATCH --job-name=compass-wlf-sfnc-ft                                       # Job name
#SBATCH --output=00_execution_script_examples/logs/slurm/slurm_fiat_%j.log   # Standard output and error log
#SBATCH --time=0-3:00:00                                                     # Job duration (hh:mm:ss)
#SBATCH --partition 4vcpu                                                    # test # 4vcpu
#SBATCH --exclusive 
#SBATCH --ntasks=1                                                           # Number of tasks (analyses) to run

module load pixi
module load julia

#Going to the folder where git checkout is
# ROOT="/u/couasnon/git_repos/compound-flooding-tropical-cyclones/"
# ROOT="/u/bovensch/git_repos/COMPASS"
ROOT="/u/vertegaa/git_repos/COMPASS"
# ROOT="/u/aleksand/compound-flooding-tropical-cyclones/"
cd "${ROOT}"

# Installing pixi environment
pixi install --environment compass-wflow-sfincs-fiat
pixi run --environment compass-wflow-sfincs-fiat pip install "hydromt_fiat @ git+https://github.com/Deltares/hydromt_fiat.git"
pixi run --environment compass-wflow-sfincs-fiat conda install ibstdcxx-ng=12 
pixi run --environment compass-wflow-sfincs-fiat conda install gcc
pixi run --environment compass-wflow-sfincs-fiat conda install --force-reinstall pandas xarray
pixi shell-hook --environment compass-wflow-sfincs-fiat > hook.sh
source hook.sh

# Make sure the correct version of libstdc++ is being loaded
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH set to: $LD_LIBRARY_PATH"

# Install Julia environment
julia +1.9 -e 'using Pkg; Pkg.instantiate(); Pkg.add("Wflow")'


# Navigate to directory where the scripts are
cd Workflows/02_workflow_rules

#Unlocking the directory for snakemake
snakemake --unlock -s snakefile_all_wflow_sfincs_fiat.smk --configfile ../01_config_snakemake/config_general_MZB.yml 

#running workflow with snakemake
snakemake -s snakefile_all_wflow_sfincs_fiat.smk --configfile ../01_config_snakemake/config_general_MZB.yml --forceall --rulegraph | dot -Tpng > dag_smk_all_mzb.png
snakemake -s snakefile_all_wflow_sfincs_fiat.smk --configfile ../01_config_snakemake/config_general_MZB.yml --cores 'all' --latency-wait 180 --wait-for-files

exit

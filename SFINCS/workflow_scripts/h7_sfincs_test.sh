#!/bin/bash
#SBATCH --job-name=test_sfincs         # Job name
#SBATCH --output=test_sfincs_%j.log     # Standard output and error log
#SBATCH --time=0-00:10:00           # Job duration (hh:mm:ss)
#SBATCH --partition test
#SBATCH --mail-user=anais.couasnon@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env
#SBATCH --exclusive


# SFINCS - example taken from: https://publicwiki.deltares.nl/display/Deltareken/SFINCS

#Going to the folder where SFINCS .inp script is
ROOT="/u/couasnon/git_repos/COMPASS/COMPASS/SFINCS/sfincs_sofala/computations/sfincs_Idai"
#cd "${ROOT}"

echo "STARTING $(date)"
docker image ls
#docker run --mount src=${PWD},target=/data,type=bind deltares/sfincs-cpu:sfincs-v2.0.3-Cauberg sfincs
docker run --mount src=${ROOT},target=/data,type=bind deltares/sfincs-cpu:sfincs-v2.0.3-Cauberg sfincs
wait
echo "ENDED $(date)"
  
#wait
# All is good
#exit 0





#Running the first sscript
#python run_test.py 
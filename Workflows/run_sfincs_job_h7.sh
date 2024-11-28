#!/bin/bash

MODELDIR=$1

module purge
module load apptainer
module load intelmpi
module load singularity

singularity exec --bind ${MODELDIR} /p/11202255-sfincs/executables/singularity/SFINCS_2023_release/SFINCS-v2.0.2-Blockhaus-Release-Q2-2023.sif sfincs

exit




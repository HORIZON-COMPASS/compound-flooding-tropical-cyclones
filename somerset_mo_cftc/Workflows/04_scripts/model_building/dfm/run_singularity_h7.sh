#! /bin/bash
 
#--- Setup the container ------------------------------------------------------------------------------------
# For use within Deltares, Delft3D FM Apptainer containers are available here: P:\d-hydro\delft3dfm_containers\
# Specify the folder that contains the required version of the Apptainer container
containerFolder=/p/d-hydro/delft3dfm_containers/delft3dfm_2025.01/
 
 
#--- Setup the model ----------------------------------------------------------------------------------------
# Specify the ROOT folder of your model, i.e. the folder that contains ALL of the input files and sub-folders, e.g:
modelFolder=MDUFOLDER

# Specify the folder containing your model's MDU file.
mdufileFolder=$modelFolder/
 
# Specify the folder containing your DIMR configuration file.
dimrconfigFolder=$modelFolder
 
# The name of the DIMR configuration file. The default name is dimr_config.xml. This file must already exist!
dimrFile=dimr_config.xml
 
# This setting might help to prevent errors due to temporary locking of NetCDF files. 
export HDF5_USE_FILE_LOCKING=FALSE


# Stop the computation after an error occurs.
set -e
 
# For parallel processes, the lines below update the <process> element in the DIMR configuration file.
# The updated list of numbered partitions is calculated from the user specified number of nodes and cores.
# You DO NOT need to modify the lines below.
PROCESSSTR="$(seq -s " " 0 $((SLURM_NTASKS-1)))"
sed -i "s/\(<process.*>\)[^<>]*\(<\/process.*\)/\1$PROCESSSTR\2/" $dimrconfigFolder/$dimrFile
 
# The name of the MDU file is read from the DIMR configuration file.
# You DO NOT need to modify the line below.
mduFile="$(sed -n 's/\r//; s/<inputFile>\(.*\).mdu<\/inputFile>/\1/p' $dimrconfigFolder/$dimrFile)".mdu
 
 
#--- Partition by calling the dflowfm executable -------------------------------------------------------------
if [ "$SLURM_NTASKS" -gt 1 ]; then 
    echo ""
    echo "Partitioning parallel model..."
    cd "$mdufileFolder"
    echo "Partitioning in folder ${PWD}"
    srun -n 1 -N 1 $containerFolder/execute_singularity_h7.sh -c $containerFolder -m $modelFolder dflowfm --nodisplay --autostartstop --partition:ndomains="$SLURM_NTASKS":icgsolver=6 "$mduFile"
else 
    #--- No partitioning ---
    echo ""
    echo "Sequential model..."
fi 
 
#--- Simulation by calling the dimr executable ----------------------------------------------------------------
echo ""
echo "Simulation..."
cd $dimrconfigFolder
echo "Computing in folder ${PWD}"
srun $containerFolder/execute_singularity_h7.sh -c $containerFolder -m $modelFolder dimr "$dimrFile"

@ echo off
rem When using mpich2 for the first time on a machine:
rem Execute "smpd -install" as administrator:
rem     Preparation: Check that your Delft3D installation contains "...\x64\share\bin\smpd.exe". Optionally copy it to a local directory (it will run as a service).
rem     "Start" -> "All programs" -> "Accessories", right-click "Command Prompt", "Run as Administrator"
rem     In this command box:
rem         cd ...\x64\share\bin
rem (i.e.   cd "c:\Program Files (x86)\Deltares\Delft3D FM Suite 2020.02 HMWQ (1.6.1.47098)\plugins\DeltaShell.Dimr\kernels\x64\share\bin\"   )
rem         smpd -install
rem     When there is an smpd already running on the machine, it must be ended first, using the Microsoft Task Manager, 
rem     or in the command  box: smpd -uninstall

rem Usage instructions for partitioned model (parallel computation)
rem  - Update the path to your Delft3D FM installation
rem  - Update the name of your mdu
rem  - Change the number of partitions
rem  - Update the settings in dimr_config.xml
rem     - Name of mdu file (e.g. <inputFile>Vietnam.mdu</inputFile>)
rem     - One process per partitions counting upward from 0 (e.g. <process>0 1</process> for 2 partitions)

rem User input
set D3D_folder="c:\Program Files\Deltares\Delft3D FM Suite 2023.02 HMWQ"
set MDU_file=Vietnam.mdu
set partitions=2

rem Partition the network and mdu
call %D3D_folder%\plugins\DeltaShell.Dimr\kernels\x64\dflowfm\scripts\run_dflowfm.bat "--partition:ndomains=%partitions%:icgsolver=6" %MDU_file%

rem Execute the simulation
call %D3D_folder%\plugins\DeltaShell.Dimr\kernels\x64\dimr\scripts\run_dimr_parallel.bat %partitions% dimr_config.xml

rem To prevent the DOS box from disappearing immediately: enable pause on the following line
pause

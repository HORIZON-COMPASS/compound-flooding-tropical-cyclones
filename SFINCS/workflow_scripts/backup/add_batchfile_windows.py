import os

dir_run = snakemake.params.dir_run

batch_content = 'call "p:/11210471-001-compass/02_Models/00_executables/SFINCS_v2.1.1_Dollerup_release_exe/sfincs.exe" > sfincs_log.txt'

file_name = 'run_sfincs.bat'
with open(os.path.join(dir_run, file_name), "w") as file:
    file.write(batch_content)
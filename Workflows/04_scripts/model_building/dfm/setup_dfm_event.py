# This script is based on the ModelBuilder template from dfm_tools v0.23.0 (accessed on 18/07/2024). 
# Adapted by Natalia Aleksandrova & Doris Vertegaal
# Use pixi environment compass-snake-dfm
#%%
# import packages
import os
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm
import xarray as xr
import pandas as pd
import geopandas as gpd
import shutil
from datetime import datetime, timedelta
import hydromt
import ast
import pyproj

#%%
if "snakemake" in locals():
    region               = snakemake.wildcards.region
    tc_name              = snakemake.params.tc_name
    dfm_res              = snakemake.wildcards.dfm_res
    bathy                = snakemake.wildcards.bathy
    tidemodel            = snakemake.wildcards.tidemodel
    bathy                = snakemake.wildcards.wind_forcing
    wind_forcing         = snakemake.wildcards.wind_forcing
    start_time           = snakemake.params.start_time
    end_time             = snakemake.params.end_time
    bbox_dfm             = ast.literal_eval(snakemake.params.dfm_bbox)
    sfincs_region        = os.path.abspath(snakemake.params.sfincs_region)
    dfm_obs_file         = snakemake.params.dfm_obs_file
    verification_points  = snakemake.params.verif_points
    path_data_cat        = snakemake.params.data_cat
    model_name           = snakemake.params.model_name
    dir_base_model       = os.path.abspath(snakemake.params.dir_base_model)
    dir_output_main      = os.path.abspath(snakemake.params.dir_event_model)
    dimrset_folder       = os.path.abspath(snakemake.params.dimrset)
    uniformwind_filename = os.path.abspath(snakemake.params.uniformwind)
    era5wind_path        = os.path.abspath(snakemake.params.era5wind)
    submit_script_file   = os.path.abspath(snakemake.output.submit_script)
    CF_SLR               = float(snakemake.wildcards.CF_SLR)
    CF_SLR_txt           = snakemake.wildcards.CF_SLR
    CF_wind              = float(snakemake.wildcards.CF_wind)
    CF_wind_txt          = snakemake.wildcards.CF_wind
else:
    region               = "sofala"
    tc_name              = "Idai"
    dfm_res              = "450"
    bathy                = "gebco2024_MZB"
    tidemodel            = 'GTSMv41' # tidemodel: FES2014, FES2012, EOT20, GTSMv41, GTSMv41opendap
    wind_forcing         = "spw_IBTrACS"
    CF_SLR               = 0
    CF_SLR_txt           = f"{CF_SLR}"
    CF_wind              = 0
    CF_wind_txt          = f"{CF_wind}"
    bbox_dfm             = ast.literal_eval("[32.3,42.5,-27.4,-9.5]")   
    sfincs_region        = f"p:/11210471-001-compass/02_Models/{region}/{tc_name}/sfincs/gis/region.geojson"
    start_time           = "20190309 000000"
    end_time             = "20190325 060000"
    dfm_obs_file         = "coastal_coupling_DFM_obs_points_MZB"
    verification_points  = "dfm_verif_points_MZB"
    path_data_cat        = [
        '../../../03_data_catalogs/datacatalog_general.yml',
        '../../../03_data_catalogs/datacatalog_SFINCS_obspoints.yml',
        '../../../03_data_catalogs/datacatalog_SFINCS_coastal_coupling.yml',
        '../../../03_data_catalogs/datacatalog_CF_forcing.yml',
        ]
    model_name           = f'event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}'
    base_model           = f'base_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR_txt}'
    dir_base_model       = f'p:/11210471-001-compass/02_Models/{region}/{tc_name}/dfm/{base_model}'
    dir_output_main      = f'p:/11210471-001-compass/03_Runs/{region}/{tc_name}/dfm/{model_name}'
    dimrset_folder       = "p:/d-hydro/dimrset/weekly/2.28.06/" # alternatively r"c:\Program Files\Deltares\Delft3D FM Suite 2023.03 HMWQ\plugins\DeltaShell.Dimr\kernels" #alternatively r"p:\d-hydro\dimrset\weekly\2.25.17.78708"
    submit_script_file   = 'run_parallel.bat'
    

#%%
# Define hydromt datacatalog
#data_catalog = hydromt.data_catalog.DataCatalog(path_data_cat)
data_catalog = hydromt.data_catalog.DataCatalog(data_libs=path_data_cat)

# Get base mdu and batchfile
script_path = os.path.dirname(os.path.realpath(__file__))
base_mdu = os.path.abspath(os.path.join(script_path, 'base_model_settings.mdu'))

#%%
# define model name, general settings and output location
dir_output_geom = os.path.join(dir_output_main,'geometry')
dir_output_bc = os.path.join(dir_output_main,'boundary_conditions')
# dir_windows_simulation = os.path.join(dir_output_main,'windows_simulation')

# make directories, if not yet present
os.makedirs(dir_output_main, exist_ok=True)

#%% Copy all files from the base model directory (excl. hidden files like .git files)
for item in os.listdir(dir_base_model):
    # Skip hidden files and directories (names starting with '.')
    if item.startswith('.'):
        continue

    src_path = os.path.join(dir_base_model, item)
    dest_path = os.path.join(dir_output_main, item)

    try:
        # Check if the item is a file or directory
        if os.path.isfile(src_path):
            shutil.copy2(src_path, dest_path)  # Copy file with metadata
        elif os.path.isdir(src_path):
            shutil.copytree(src_path, dest_path, dirs_exist_ok=True)  # Recursively copy directories
    except PermissionError:
        print(f"Permission denied: {src_path}")
    except Exception as e:
        print(f"Error copying {src_path}: {e}")
#%%
os.makedirs(dir_output_geom, exist_ok=True)
# os.makedirs(dir_windows_simulation, exist_ok=True)

# generate_grid = False # option to skip grid generation if this was already done.
overwrite = False # used for downloading of forcing data. Always set to True when changing the domain
crs = 'EPSG:4326' # coordinate reference system

#%%
# Define model domain
lon_min, lon_max, lat_min, lat_max = bbox_dfm
 
#%%
# Define spin up and the model dates in the correct format
spin_up =  3 # days

start_datetime = datetime.strptime(start_time, '%Y%m%d %H%M%S')
end_datetime = datetime.strptime(end_time, '%Y%m%d %H%M%S')
date_min = (start_datetime - timedelta(days=spin_up)).strftime('%Y-%m-%d %H:%M:%S')
date_max = (end_datetime).strftime('%Y-%m-%d %H:%M:%S')
ref_date = datetime(start_datetime.year, 1, 1).strftime('%Y-%m-%d %H:%M:%S')

#%% #########################################
###### Grid generation and refinement #######
#############################################
netfile = os.path.join(dir_output_main, 'grid_network.nc')
poly_file = os.path.join(dir_output_main, 'pli_file.pli')
pathfile_illegalcells = os.path.join(dir_output_main, "illegalcells.pol")

# Load the base model grid
xu_grid_uds = dfmt.open_partitioned_dataset(netfile)

#%%##################################################
#### Define boundary conditions from tidal model ####
#####################################################

# Modify the ext_file_new to contain only file names and not full paths for Linux
ext_file_new = os.path.join(dir_output_main,'ext_file_new.ext')

#%%#####################################
######### Define meteo forcing #########
########################################
# Function to more easily create meteo forcing with full paths and only file names
def create_forcing(quantity, filename, filetype, method, operand, use_basename=False):
    # Optionally use only the file name
    if use_basename:
        filename = os.path.basename(filename)
    return hcdfm.ExtOldForcing(quantity=quantity,
                               filename=filename,
                               filetype=filetype,
                               method=method,
                               operand=operand)

# generate old format external forcings file (.ext): spatial data
ext_file_old = os.path.join(dir_output_main, f'ext_file_old.ext')

# Initialise focring files for Linux and Windows separately
ext_old = hcdfm.ExtOldModel()


# Define model forcing
if 'spw' in wind_forcing and 'era5' in wind_forcing:
    meteo_type = 'spiderweb_era5_merged'
    spw = 1 
    spw_input = data_catalog[f"spw_IBTrACS_CF{CF_wind_txt}_{tc_name}"].path
    spw_file_origin = spw_input # change to path from datacatalog
elif 'spw' in wind_forcing:
    meteo_type = 'spiderweb'
    spw = 1 
    spw_input = data_catalog[f"spw_IBTrACS_CF{CF_wind_txt}_{tc_name}"].path
    spw_file_origin = spw_input # change to path from datacatalog
else:
    meteo_type = wind_forcing
    spw = 0
# Can we add option for spiderweb+ERA5? Could not make it work yet. -> Natalia: as far as I know that is not possible

# Create the forcing file (ext_old) depending on the wind data
if meteo_type == 'spiderweb_era5_merged':
    spw_file = os.path.basename(spw_file_origin)
    spw_copy = os.path.join(dir_output_main,spw_file)
    shutil.copyfile(spw_file_origin, spw_copy)

    # Add uniform wind file to set background wind speed to 0 and enable blending of the spw file with the background wind
    preprocess_era5_Idai_path = era5wind_path
    preprocess_era5_Idai_file = os.path.basename(preprocess_era5_Idai_path)
    preprocess_era5_Idai_dest = os.path.join(dir_output_main, preprocess_era5_Idai_file)
    shutil.copyfile(preprocess_era5_Idai_path, preprocess_era5_Idai_dest)

    # Create forcing with only file names for Linux
    ext_old.forcing.append(hcdfm.ExtOldForcing(quantity='airpressure_windx_windy',
                                                filename=preprocess_era5_Idai_file,
                                                varname='msl u10n v10n',
                                                filetype=hcdfm.ExtOldFileType.NetCDFGridData, #11
                                                method=hcdfm.ExtOldMethod.InterpolateTimeAndSpaceSaveWeights, #3
                                                operand=hcdfm.Operand.override, #O
                                                ))
    # Operand add (+) has been tested for desired spw+ERA5 behaviour. Override was not successful for ERA5 merging
    ext_old.forcing.append(create_forcing('airpressure_windx_windy', spw_copy, hcdfm.ExtOldFileType.SpiderWebData, 
                                          hcdfm.ExtOldMethod.PassThrough, hcdfm.Operand.add, use_basename=True))
    ext_old.save(filepath=ext_file_old) # save the file

    # Add setting to spw to be merge with era5 at the radius of 0.75
    with open(spw_copy, "r") as f:
        lines = f.readlines()

    has_merge_frac = any(line.strip().startswith("spw_merge_frac") for line in lines)
    new_lines = []

    for i, line in enumerate(lines):
        new_lines.append(line)
        if line.strip().startswith("spw_radius") and not has_merge_frac:
            # Add the line right after spw_radius
            new_lines.append("spw_merge_frac = 0.75\n")
            has_merge_frac = True  # To avoid adding multiple times

    # Write to new file or overwrite original
    with open(spw_copy, "w") as f:
        f.writelines(new_lines)

    print(f"Added spw_merge_frac to: {spw_copy}")
    

elif 'era5' in meteo_type: # ERA5 - download spatial fields of air pressure, wind speeds and Charnock coefficient
    dir_output_data_era5 = os.path.join(dir_output_bc,'meteo', 'ERA5')
    os.makedirs(dir_output_data_era5, exist_ok=True)
        
    varlist_list = [['msl','u10n','v10n']]
    for varlist in varlist_list:
        for varkey in varlist:
            #TODO change to use data catalog?
            dfmt.download_ERA5(varkey, 
                            longitude_min=lon_min, longitude_max=lon_max, latitude_min=lat_min, latitude_max=lat_max,
                            date_min=date_min, date_max=date_max,
                            dir_output=dir_output_data_era5, overwrite=overwrite)

    # ERA5 meteo - convert to netCDF for usage in Delft3D FM
    ext_old = dfmt.preprocess_merge_meteofiles_era5(ext_old=ext_old,
                                                    varkey_list=varlist_list,
                                                    dir_data=dir_output_data_era5,
                                                    dir_output=dir_output_main,
                                                    time_slice=slice(date_min, date_max))
    ext_old.save(filepath=ext_file_old) # save the file

elif meteo_type == 'spiderweb':
    spw_file = os.path.basename(spw_file_origin)
    spw_copy = os.path.join(dir_output_main,spw_file)
    shutil.copyfile(spw_file_origin, spw_copy)

    # Add uniform wind file to set background wind speed to 0 and enable blending of the spw file with the background wind
    uniwind_path = os.path.join(dir_output_main,"uniformwind0.wnd")
    with open(uniwind_path, 'w', encoding='utf-8') as file:
        file.write('0 0 0\n 1728000 0 0')    

    # Create forcing with only file names for Linux
    ext_old.forcing.append(create_forcing('windxy', uniwind_path, hcdfm.ExtOldFileType.TimeSeries, 
                                      hcdfm.ExtOldMethod.PassThrough, hcdfm.Operand.override, use_basename=True))
    ext_old.forcing.append(create_forcing('airpressure_windx_windy', spw_copy, hcdfm.ExtOldFileType.SpiderWebData, 
                                      hcdfm.ExtOldMethod.PassThrough, hcdfm.Operand.add, use_basename=True))
    ext_old.save(filepath=ext_file_old) # save the file

#%%##############################################
############## Generate obs file ################
#################################################
# The D-Flow FM model wil have mapoutput and hisoutput. 
# A file with coordinates of obs stations will be generated.

# Read shp file of points along the coastline and the sfincs region
gdfp = gpd.read_file(data_catalog[dfm_obs_file].path)
region_poly = gpd.read_file(sfincs_region)
region_poly = region_poly.to_crs(crs)

# Ensure both have the same CRS
gdfp = gdfp.to_crs(crs)

# Buffer the region_poly and clip the points along the coastline by the buffered region
buffer_size = 0.01  # Set buffer distance (adjust as needed)
region_poly_buffered = region_poly.buffer(buffer_size)
gdfp_clipped = gpd.clip(gdfp, region_poly_buffered)

# Extract x and y coordinates
xcor = gdfp_clipped.geometry.x; xcor.name = 'x'
ycor = gdfp_clipped.geometry.y; ycor.name = 'y'

# Combine into a DataFrame
tmp = pd.concat([xcor, ycor], axis=1).dropna()
tmp['names'] = tmp.index

try:
    if verification_points:  # Ensures 'verification_points' is not empty or None
        tmp2 = pd.read_table(data_catalog[verification_points].path, sep=" ", names=['x', 'y', 'names'])
        pd_obs = tmp2._append(tmp, ignore_index=True)  
        del tmp, tmp2
    else:
        raise ValueError("verification_points is empty.")  # Handle empty case explicitly
except NameError:
    # If 'verification_points' is not defined, set pd_obs to tmp
    pd_obs = tmp
    print("verification_points is not defined. Using 'tmp' as pd_obs.")
except ValueError as e:
    # If 'verification_points' is defined but empty, set pd_obs to tmp
    pd_obs = tmp
    print(f"{e} Using 'tmp' as pd_obs.")

# save obs points
file_obs = os.path.join(dir_output_main, f'obs_points.xyn')
pd_obs.to_csv(file_obs, sep=' ', header=False, index=False, float_format='%.6f')

# plot obs points
fig, ax = plt.subplots(figsize=(8,4))
xu_grid_uds.grid.plot(ax=ax,linewidth=0.5,color='k',alpha=0.2)
ax.plot(pd_obs['x'],pd_obs['y'],'rx')
region_poly.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=1)
dfmt.plot_coastlines(ax=ax, crs=region_poly.crs)

# Save the figure
output_path = os.path.join(dir_output_geom, 'grid_obs_points.png')
fig.savefig(output_path, dpi=300, bbox_inches='tight')

#%%#########################################
############ Generate mdu file #############
############################################
# In order for the model to run, we need a model definition file, i.e., a *.mdu file

# initialize mdu file and update settings
mdu_file = os.path.join(dir_output_main, f'settings.mdu')


# use mdu file from GTSM
base_mdu = base_mdu
mdu = hcdfm.FMModel(base_mdu)

if os.path.exists(pathfile_illegalcells):
    mdu.geometry.drypointsfile = pathfile_illegalcells


# add the grid (grid_network.nc, network file)
mdu.geometry.netfile = netfile
# Adjust initial water in case of SLR for Idai
if CF_SLR == 0: # factual
    mdu.geometry.waterlevini = 0.04
elif CF_SLR == -0.05:
    mdu.geometry.waterlevini = -0.01
elif CF_SLR == -0.10:
    mdu.geometry.waterlevini = -0.06
elif CF_SLR == -0.15:
    mdu.geometry.waterlevini = -0.11
else:
    pass

# add the external forcing files (.ext)
mdu.external_forcing.extforcefile = ext_file_old
mdu.external_forcing.extforcefilenew = ext_file_new

# Define drag coefficient 
mdu.wind.icdtyp = 3 
mdu.wind.cdbreakpoints = [0.001, 0.003, 0.0015]
mdu.wind.windspeedbreakpoints = [0, 25, 50]

# update time settings
mdu.time.refdate = pd.Timestamp(ref_date).strftime('%Y%m%d')
mdu.time.tunit = 'S'
mdu.time.dtmax = 30
mdu.time.startdatetime = pd.Timestamp(date_min).strftime('%Y%m%d%H%M%S')
mdu.time.stopdatetime = pd.Timestamp(date_max).strftime('%Y%m%d%H%M%S')
mdu.time.autotimestep = 3

# update output settings
mdu.output.obsfile = file_obs
mdu.output.hisinterval = [600] #s
mdu.output.mapinterval = [3600]
mdu.output.rstinterval = [0]
mdu.output.statsinterval = [3600]

# output wind
mdu.output.wrimap_wind = 1

# save .mdu file
mdu.save(mdu_file) 

#%% 
# Modify the ext_new file for Linux simulation (only containing file names and not full paths)

# make all paths relative (might be properly implemented in https://github.com/Deltares/HYDROLIB-core/issues/532)
dfmt.make_paths_relative(mdu_file)

#%%####################################################
############# Generate DIMR and bat file ##############
#######################################################
# In order to run the model via DIMR we need a `dimr_config.xml` file. 
# If you are running this notebook on a Windows platform, a *.bat file will also be created

nproc = 4 # number of processes

# Making bat file and dimr_config.xml file
dfmt.create_model_exec_files(file_mdu=mdu_file, nproc=nproc, dimrset_folder=dimrset_folder.replace('/', '\\'))
# default_bat_path = os.path.join(dir_output_main, "run_parallel.bat")
bat_file_path = os.path.join(dir_output_main, "run_parallel.bat")

# Remove "pause" from bat file and update MDU_file, dimr_config
if os.path.exists(bat_file_path):
    print(f"Found .bat file: {bat_file_path}. Modifying...")

    # Open the .bat file and read its lines
    with open(bat_file_path, "r") as infile:
        lines = infile.readlines()

    # Create a new file (or overwrite the existing file) with the changes
    with open(bat_file_path, "w") as outfile:
        # Add the working directory as the first line
        outfile.write(f'Rem Set working directory\n')
        outfile.write(f'cd /d "{os.path.normpath(dir_output_main)}"\n')
        outfile.write(f'\n')

        dimr_config_line_added = False  # Flag to check if the dimr_config line is added

        for line in lines:
            # Remove the pause command
            if line.strip().lower() == "pause":
                continue  # Skip writing this line

            # Replace every "/p/" with "p:\" in the line (happens when run on linux)
            line = line.replace("/p/", "p:\\")

            # Write any other lines as they are
            outfile.write(line)

    print(f"Updated the .bat file: {bat_file_path}")
else:
    print(f".bat file not found: {bat_file_path}. No changes made.")


# making singularity .sh file
if submit_script_file[-3:]=='.sh':

    replacements = {'JOBNAME': region, 'MDUFOLDER':os.path.dirname(mdu_file).replace('\\', '/').replace('p:/', '/p/')}

    with open(os.path.join(script_path,os.path.basename(submit_script_file))) as infile, open(submit_script_file, 'w',newline='\n') as outfile:
        for line in infile:
            for src, target in replacements.items():
                line = line.replace(src, target)
            outfile.write(line)

#%%############################################
############ Visualize model tree #############
###############################################

# visualize the model tree, show_tree is available for all HYDROLIB-core model components
mdu_obj = hcdfm.FMModel(mdu_file)
mdu_obj.show_tree()
# %%

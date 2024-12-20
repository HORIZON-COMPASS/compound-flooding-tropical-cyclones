# This script is based on the ModelBuilder template from dfm_tools v0.23.0 (accessed on 18/07/2024). 
# Adapted by Natalia Aleksandrova
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

#%%
if "snakemake" in locals():
    model_dir = snakemake.params.dir_model
    # config_file = snakemake.input.config_file
    # data_cat = snakemake.params.data_cat
    # bbox = snakemake.params.arg_bbox
    # region_geom = snakemake.input.region_geom
    # dir_sfincs_model = snakemake.input.dir_sfincs_model
else:
    bbox_dfm = "[32.3,42.5,-27.4,-9.5]"
    region = "sofala"
    tc_name = "Idai"
    dfm_res = "450"
    bathy = "gebco2024"
    tidemodel = 'GTSMv4.1_opendap' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap
    wind_forcing = "spw_IBTrACS_ext"
    start_time = "20190309 000000"
    end_time = "20190325 060000"
    model_name = f'mozambique_{region}_{tc_name}_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}'
    base_model = model_name = f'base_{dfm_res}_{bathy}_{tidemodel}'
    dir_base_model = f'p:/11210471-001-compass/02_Models/mozambique/dfm/{model_name}'
    dir_output_main = f'p:/11210471-001-compass/02_Models/{region}/{tc_name}/dfm/{model_name}'
    spw_input = 'p:/11210471-001-compass/02_Models/Delft3DFM/mozambique_model/boundary_conditions/meteo/TC/tc_IDAI_2019063S18038_ext9d.spw'
    output_bbox = [34, -20.5, 35.6, -19.5] 
    path_data_cat = os.path.abspath("../../../data_catalogs/datacatalog_general.yml")
    dfm_obs_file = 'p:/11210471-001-compass/01_Data/Coastal_boundary/points/coastal_bnd_MZB_5mMSL_points_1km.shp'
    verification_points = 'p:/11210471-001-compass/02_Models\Delft3DFM/mozambique_model/geometry/output_locations/MZB_Sofala_IHO_obs.xyn'

#%%
# Define hydromt datacatalog
data_catalog = hydromt.DataCatalog(data_libs = [path_data_cat])

# needed to define wind_forcing and bathy paths
#%%
# define model name, general settings and output location
dir_output_run = os.path.join(dir_output_main,'computations')
dir_output_geom = os.path.join(dir_output_main,'geometry')
dir_output_bc = os.path.join(dir_output_main,'boundary_conditions')

# make directories, if not yet present
os.makedirs(dir_output_run, exist_ok=True)

# generate_grid = False # option to skip grid generation if this was already done.
overwrite = False # used for downloading of forcing data. Always set to True when changing the domain
crs = 'EPSG:4326' # coordinate reference system

#%%
# domain and resolution
bbox_list = [float(x) for x in bbox_dfm.strip("[]").split(",")]
lon_min, lon_max, lat_min, lat_max = bbox_list

# dxy is the base grid resolution - i.e. the coarsest grid size in your grid. It is in degrees.
# dxy = 0.2 # degrees
# # Defines the minimum length of the edges of the computational mesh cells.
# min_edge_size = dfm_res # m 
#%%
# Define spin up and the model dates in the correct format
spin_up =  4 # days

start_datetime = datetime.strptime(start_time, '%Y%m%d %H%M%S')
end_datetime = datetime.strptime(end_time, '%Y%m%d %H%M%S')
date_min = (start_datetime - timedelta(days=spin_up)).strftime('%Y-%m-%d %H:%M:%S')
date_max = (end_datetime).strftime('%Y-%m-%d %H:%M:%S')
ref_date = datetime(start_datetime.year, 1, 1).strftime('%Y-%m-%d %H:%M:%S')


#%% #########################################
###### Grid generation and refinement #######
#############################################

# Copy the correct base model files and define them
shutil.copyfile(os.path.join(dir_base_model, base_model, 'grid_network.nc'), os.path.join(dir_output_run,'grid_network.nc'))
shutil.copyfile(os.path.join(dir_base_model, base_model, 'pli_file.pli'), os.path.join(dir_output_run,'pli_file.pli'))
shutil.copyfile(os.path.join(dir_base_model, base_model, 'illegalcells.pol'), os.path.join(dir_output_run,'illegalcells.pol'))

netfile = os.path.join(dir_output_run, 'grid_network.nc')
poly_file = os.path.join(dir_base_model, base_model, 'pli_file.pli')
pathfile_illegalcells = os.path.join(dir_base_model, base_model,"illegalcells.pol")

# Load the base model grid
xu_grid_uds = dfmt.open_partitioned_dataset(netfile)

#%%##################################################
#### Define boundary conditions from tidal model ####
#####################################################

# Copy the correct base model external forcings file (.ext): initial and open boundary condition 
shutil.copyfile(os.path.join(dir_base_model, base_model, 'ext_file_new.ext'), os.path.join(dir_output_run,'ext_file_new.ext'))

# and define it
ext_file_new = os.path.join(dir_output_run, f'ext_file_new.ext')

#%%#####################################
######### Define meteo forcing #########
########################################

# generate old format external forcings file (.ext): spatial data
ext_file_old = os.path.join(dir_output_run, f'ext_file_old.ext')

ext_old = hcdfm.ExtOldModel()

# Define model forcing
if 'spw' in wind_forcing:
    meteo_type = 'spiderweb'
    spw = 1
    spw_input = data_catalog[f'{wind_forcing}_{tc_name}'].path
    spw_file_origin = spw_input # change to path from datacatalog
else:
    meteo_type = wind_forcing
    spw = 0
# Can we add option for spiderweb+ERA5? Could not make it work yet.

if meteo_type=='ERA5': # ERA5 - download spatial fields of air pressure, wind speeds and Charnock coefficient
    dir_output_data_era5 = os.path.join(dir_output_bc,'meteo', 'ERA5')
    os.makedirs(dir_output_data_era5, exist_ok=True)
        
    varlist_list = [['msl','u10n','v10n','chnk']]
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
                                                    dir_output=dir_output_run,
                                                    time_slice=slice(date_min, date_max))
elif meteo_type == 'spiderweb':
    spw_file = os.path.basename(spw_file_origin)
    shutil.copyfile(spw_file_origin, os.path.join(dir_output_run,spw_file))

    uniformwind_filename = "p:/11210471-001-compass/01_Data/uniformwind0.wnd"
    shutil.copyfile(uniformwind_filename,os.path.join(dir_output_run,"uniformwind0.wnd"))

    forcing_uniformwind = hcdfm.ExtOldForcing(quantity='windxy',
                                        filename=uniformwind_filename,
                                        filetype=hcdfm.ExtOldFileType.TimeSeries, 
                                        method=hcdfm.ExtOldMethod.PassThrough, 
                                        operand=hcdfm.Operand.override)
    ext_old.forcing.append(forcing_uniformwind)
    
    forcing_spw = hcdfm.ExtOldForcing(quantity='airpressure_windx_windy',
                                        filename=spw_file,
                                        filetype=hcdfm.ExtOldFileType.SpiderWebData, 
                                        method=hcdfm.ExtOldMethod.PassThrough, 
                                        operand=hcdfm.Operand.add)
    ext_old.forcing.append(forcing_spw)

ext_old.save(filepath=ext_file_old) # , path_style=path_style)

#%%##############################################
############## Generate obs file ################
#################################################
# The D-Flow FM model wil have mapoutput and hisoutput. 
# A file with coordinates of obs stations will be generated.

# Read shp file of points along the MZB coastline
gdfp = gpd.read_file(dfm_obs_file)

# crop output points to the area where output is needed for the flood model
gdfp = gdfp.cx[output_bbox[0]:output_bbox[2],output_bbox[1]:output_bbox[3]]

# Convert points to the xyn file format
xcor = gdfp['geometry'][:].x; xcor.name = 'x'
ycor = gdfp['geometry'][:].y; ycor.name = 'y'
tmp = pd.concat([xcor,ycor],axis=1)
tmp = tmp.dropna()
tmp['names'] = tmp.index

# Attempt to check if 'verification_points' is defined
try:
    # Try to use 'verification_points'
    tmp2 = pd.read_table(verification_points, sep=" ", names=['x', 'y', 'names'])
    pd_obs = tmp2._append(tmp, ignore_index=True)  
    del tmp, tmp2
except NameError:
    # If 'verification_points' is not defined, set pd_obs to tmp
    pd_obs = tmp
    print("verification_points is not defined. Using 'tmp' as pd_obs.")

# save obs points
file_obs = os.path.join(dir_output_run, f'obs_points.xyn')
pd_obs.to_csv(file_obs, sep=' ', header=False, index=False, float_format='%.6f')

# plot obs points
fig, ax = plt.subplots(figsize=(8,4))
xu_grid_uds.grid.plot(ax=ax,linewidth=0.5,color='k',alpha=0.2)
ax.plot(pd_obs['x'],pd_obs['y'],'rx')
dfmt.plot_coastlines(ax=ax, crs=crs)


#%%#########################################
############ Generate mdu file #############
############################################
# In order for the model to run, we need a model definition file, i.e., a *.mdu file

# initialize mdu file and update settings
mdu_file = os.path.join(dir_output_run, f'settings.mdu')

#mdu = hcdfm.FMModel()

# use mdu file from GTSM
base_mdu = os.path.abspath('base_model_settings.mdu')
mdu = hcdfm.FMModel(base_mdu)

if os.path.exists(pathfile_illegalcells):
    mdu.geometry.drypointsfile = pathfile_illegalcells

# add the grid (grid_network.nc, network file)
mdu.geometry.netfile = netfile

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
mdu.save(mdu_file) # ,path_style=path_style)

# make all paths relative (might be properly implemented in https://github.com/Deltares/HYDROLIB-core/issues/532)
dfmt.make_paths_relative(mdu_file)


#%%####################################################
############# Generate DIMR and bat file ##############
#######################################################
# In order to run the model via DIMR we need a `dimr_config.xml` file. 
# If you are running this notebook on a Windows platform, a *.bat file will also be created

nproc = 4 # number of processes
dimrset_folder = r"p:\d-hydro\dimrset\weekly\2.25.17.78708" # alternatively r"c:\Program Files\Deltares\Delft3D FM Suite 2023.03 HMWQ\plugins\DeltaShell.Dimr\kernels" #alternatively r"p:\d-hydro\dimrset\weekly\2.25.17.78708"
dfmt.create_model_exec_files(file_mdu=mdu_file, nproc=nproc, dimrset_folder=dimrset_folder)

# add a batch file to submit job to h7 cluster (as al alternative to running it on your node or computer)
batchfile_h7 = os.path.abspath("submit_singularity_h7.sh")

# maybe not necessary?
pathfile_h7 = os.path.join(dir_output_run,'submit_singularity_h7.sh')

replacements = {'JOBNAME':'MZB', 'MDUFILE':mdu_file}

with open(batchfile_h7) as infile, open(pathfile_h7, 'w',newline='\n') as outfile:
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
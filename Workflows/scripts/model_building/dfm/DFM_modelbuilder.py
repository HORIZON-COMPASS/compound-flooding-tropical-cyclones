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

#%%
if "snakemake" in locals():
    # model_dir = snakemake.params.dir_model
    # config_file = snakemake.input.config_file
    # data_cat = snakemake.params.data_cat
    # bbox = snakemake.params.arg_bbox
    # region_geom = snakemake.input.region_geom
    # dir_sfincs_model = snakemake.input.dir_sfincs_model
else:
    dir_output_main = r'p:\11210471-001-compass\02_Models\Delft3DFM\mozambique_model'
    model_name = 'mozambique_spw_Freddy2_area_MZB_500m_gebco2024'

#%%
# user input - to be converted to snakemake inputs and params


generate_grid = False # option to skip grid generation if this was already done.
overwrite = False # used for downloading of forcing data. Always set to True when changing the domain
crs = 'EPSG:4326' # coordinate reference system

#dir_output_main = os.path.abspath('..')

dir_output_run = os.path.join(dir_output_main,'computations',model_name)
dir_output_geom = os.path.join(dir_output_main,'geometry')
dir_output_bc = os.path.join(dir_output_main,'boundary_conditions')

# path_style = 'windows' # windows / unix

# domain and resolution
# dxy is the base grid resolution - i.e. the coarsest grid size in your grid. It is in degrees.
# min_edge_size is the 
if (model_name=='mozambique_large') | (model_name=='mozambique_spw_Idai') | (model_name=='mozambique_spw_Idai_differentspwfile'):
    lon_min, lon_max, lat_min, lat_max =  32.55, 47.22, -27.84, -9.15
    dxy = 0.24; min_edge_size = 1500 # in meters
elif model_name=='mozambique_spw_Idai_highresgrid':
    lon_min, lon_max, lat_min, lat_max =  32.55, 47.22, -27.84, -9.6
    dxy = 0.1; min_edge_size = 400 # in meters

elif 'areaBeira_500m' in model_name:
    lon_min, lon_max, lat_min, lat_max = 34, 40.3, -22.79,-15.98
    dxy = 0.2 ; min_edge_size = 450 # m 
elif 'area_MZBlarge_500m' in model_name:
    lon_min, lon_max, lat_min, lat_max = 32.3, 45.5, -27.4,-9.5
    dxy = 0.2 # degrees
    min_edge_size = 450 # m 
elif 'area_MZBlarge_200m' in model_name:
    lon_min, lon_max, lat_min, lat_max = 32.3, 45.5, -27.4,-9.5
    dxy = 0.2 # degrees
    min_edge_size = 200 # m    
elif 'area_MZB_500m' in model_name:
    lon_min, lon_max, lat_min, lat_max = 32.3, 42.5, -27.4,-9.5
    dxy = 0.2 # degrees
    min_edge_size = 450 # m  
elif 'area_MZB_200m' in model_name:
    lon_min, lon_max, lat_min, lat_max = 32.3, 42.5, -27.4,-9.5
    dxy = 0.2 # degrees
    min_edge_size = 200 # m      

#dates as understood by pandas.period_range(). 
if 'Idai' in model_name:
    date_min = '2019-03-05 00:00:00' 
    date_max = '2019-03-25 06:00:00'
    ref_date = '2019-01-01'
elif 'Kenneth' in model_name:
    date_min = '2019-04-21 21:00:00' 
    date_max = '2019-05-05 00:00:00'
    ref_date = '2019-01-01' 
elif 'Freddy1' in model_name:
    date_min = '2023-02-21 00:00:00'
    date_max = '2023-03-05 06:00:00'
    ref_date = '2023-01-01'    
elif 'Freddy2' in model_name:
    date_min = '2023-03-04 06:00:00'
    date_max = '2023-03-23 00:00:00'
    ref_date = '2023-01-01'    

# forcing
if model_name=='mozambique_large':
    meteo_type = 'ERA5'
    spw = 0
elif 'spw' in model_name:
    meteo_type = 'spiderweb'
    spw = 1
# Can we add option for spiderweb+ERA5? Could not make it work yet.

if 'spw_Idai' in model_name:    
    spw_name = 'tc_Idai'
    spw_file_origin = r'p:\11210471-001-compass\02_Models\Delft3DFM\mozambique_model\boundary_conditions\meteo\TC\tc_IDAI_2019063S18038_ext9d.spw'
    output_bbox = [34, -20.5, 35.6, -19.5] 
elif 'spw_Kenneth' in model_name:    
    spw_name = 'tc_Kenneth'
    spw_file_origin = r'p:\11210471-001-compass\02_Models\Delft3DFM\mozambique_model\boundary_conditions\meteo\TC\tc_KENNETH_2019112S10053_ext9d.spw'
    output_bbox = [40, -13.5, 42, -11.2] 
elif 'spw_Freddy1' in model_name:    
    spw_name = 'tc_Freddy1'
    spw_file_origin = r'p:\11210471-001-compass\02_Models\Delft3DFM\mozambique_model\boundary_conditions\meteo\TC\tc_FREDDY_2023036S12118_ext9d.spw'
    output_bbox = [34, -23.3, 36, -21] 
elif 'spw_Freddy2' in model_name:    
    spw_name = 'tc_Freddy2'
    spw_file_origin = r'p:\11210471-001-compass\02_Models\Delft3DFM\mozambique_model\boundary_conditions\meteo\TC\tc_FREDDY_2023061S22036_ext9d.spw'
    output_bbox = [36.3, -19.7, 38, -16.7] 

# Choose tidal model
tidemodel = 'GTSMv4.1_opendap' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap

#%%    
# make directories, if not yet present
os.makedirs(dir_output_run, exist_ok=True)


#%%
#############################################
###### Grid generation and refinement #######
#############################################

# File name for the grid network
netfile = os.path.join(dir_output_run, f'{model_name}_net.nc')

# Load grid if it was already generated
if not generate_grid:
    xu_grid_uds = dfmt.open_partitioned_dataset(netfile)

#%%
if generate_grid:
    # generate spherical regular grid
    mk_object = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy, crs=crs)

    # generate plifile from grid extent and coastlines
    bnd_gdf = dfmt.generate_bndpli_cutland(mk=mk_object, res='h', buffer=0.2)
    bnd_gdf_interp = dfmt.interpolate_bndpli(bnd_gdf, res=0.03)

#%%
poly_file = os.path.join(dir_output_run, f'{model_name}.pli')
if generate_grid:
    # filter out boundary sections that are very short (<5km in this case)
    ids = []
    for ii in range(len(bnd_gdf_interp.index)):
        if bnd_gdf_interp.iloc[ii].values.to_crs('EPSG:3857').length < 5000: # approx. 111 km in 1 degree
            ids.append(ii)

    bnd_gdf_interp = bnd_gdf_interp.drop(ids)
    bnd_gdf_interp

    pli_polyfile = dfmt.geodataframe_to_PolyFile(bnd_gdf_interp)
    pli_polyfile.save(poly_file)

#%%
if generate_grid:
    # plot basegrid and polyline
    fig, ax = plt.subplots()
    mk_object.mesh2d_get().plot_edges(ax,zorder=1)
    bnd_gdf_interp.plot(ax=ax, edgecolor='r')
    # ctx.add_basemap(ax=ax, crs=crs, attribution=False)
    dfmt.plot_coastlines(ax=ax, crs=crs)

#%%
# connect to a coarse version of the GEBCO_2022 dataset on OPeNDAP
# alternatively download your own full resolution cutout from https://download.gebco.net (use a buffer of e.g. 1 degree)
#file_nc_bathy = "https://opendap.deltares.nl/thredds/dodsC/opendap/deltares/Delft3D/netcdf_example_files/GEBCO_2022/GEBCO_2022_coarsefac08.nc"
if generate_grid:
    if 'gebco2020' in model_name:
        file_bathy = r'p:\metocean-data\open\GEBCO\2020\GEBCO_2020.nc'
        data_bathy = xr.open_dataset(file_bathy).elevation
        data_bathy_sel = data_bathy.sel(lon=slice(lon_min-1, lon_max+1), lat=slice(lat_min-1, lat_max+1),drop=True)
    elif 'gebco2023' in model_name:
        file_nc_bathy_sel = r'p:\11210471-001-compass\01_Data\Gebco\gebco_2023_n-8.0_s-29.0_w31.5_e48.5.nc'
        data_bathy_sel = xr.open_dataset(file_nc_bathy_sel).elevation
    elif 'gebco2024' in model_name:
        file_nc_bathy_sel = r'p:\11210471-001-compass\01_Data\Gebco\gebco_2024_n-8.0_s-30.0_w30.0_e50.0.nc'
        data_bathy_sel = xr.open_dataset(file_nc_bathy_sel).elevation
        
    data_bathy_sel.load()

# alternatively you can connect to ETOPO 30s, for which there is also a 15s (15 arcseconds) resolution dataset available
# file_nc_bathy = "https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/30s/30s_surface_elev_netcdf/ETOPO_2022_v1_30s_N90W180_surface.nc"
# data_bathy = xr.open_dataset(file_nc_bathy).z

# subset to area of interest
#data_bathy_sel = data_bathy.sel(lon=slice(lon_min-1, lon_max+1), lat=slice(lat_min-1, lat_max+1))

#%%
if generate_grid:
    # refine grid
    dfmt.refine_basegrid(mk=mk_object, data_bathy_sel=data_bathy_sel, min_edge_size=min_edge_size)

    # plot
    fig, ax = plt.subplots()
    mk_object.mesh2d_get().plot_edges(ax,zorder=1)
    # ctx.add_basemap(ax=ax, crs=crs, attribution=False)
    dfmt.plot_coastlines(ax=ax, crs=crs)

#%%
if generate_grid:
    # remove land with GSHHS coastlines
    dfmt.meshkernel_delete_withcoastlines(mk=mk_object,res='i',min_area=50,crs=crs)

#%%
if generate_grid:
    # plot
    fig, ax = plt.subplots(figsize=(12,7))
    mk_object.mesh2d_get().plot_edges(ax,zorder=1)
    # ctx.add_basemap(ax=ax, crs=crs, attribution=False)
    dfmt.plot_coastlines(ax=ax,res='i',min_area=50,crs=crs)
    #ax.set_xlim(46, 46.75)
    #ax.set_ylim(-16.25,-15.50)

#%%
pathfile_illegalcells = os.path.join(dir_output_run,"illegalcells.pol")
if generate_grid:
    # derive illegalcells geodataframe
    illegalcells_gdf = dfmt.meshkernel_get_illegalcells(mk=mk_object)
    # create and add drypointsfile if there are any cells generated that will result in high orthogonality
    if len(illegalcells_gdf) > 0:
        illegalcells_polyfile = dfmt.geodataframe_to_PolyFile(illegalcells_gdf)
        illegalcells_polyfile.save(pathfile_illegalcells)

#%%
if generate_grid:
    # convert to xugrid
    xu_grid_uds = dfmt.meshkernel_to_UgridDataset(mk=mk_object, crs=crs)

    # interpolate bathymetry onto the grid
    data_bathy_interp = data_bathy_sel.interp(lon=xu_grid_uds.obj.mesh2d_node_x, lat=xu_grid_uds.obj.mesh2d_node_y)
    xu_grid_uds['mesh2d_node_z'] = data_bathy_interp.clip(max=10)

    # write xugrid grid to netcdf
    xu_grid_uds.ugrid.to_netcdf(netfile)

#%%
# plot bathymetry and grid
fig, ax = plt.subplots(figsize=(8,4))
xu_grid_uds.mesh2d_node_z.ugrid.plot(ax=ax,center=False)
xu_grid_uds.grid.plot(ax=ax,linewidth=0.5,color='white',alpha=0.2)
# ctx.add_basemap(ax=ax, crs=crs, attribution=False)
dfmt.plot_coastlines(ax=ax, crs=crs)


#%%
#####################################################
### Generate boundary conditions from tidal model ###
#####################################################

# generate new format external forcings file (.ext): initial and open boundary condition
ext_file_new = os.path.join(dir_output_run, f'{model_name}_new.ext')
ext_new = hcdfm.ExtModel()

#%%
# interpolate tidal components to boundary conditions file (.bc)
dfmt.interpolate_tide_to_bc(ext_new=ext_new, tidemodel=tidemodel, file_pli=poly_file, component_list=None)

#%%
#save new ext file
ext_new.save(filepath=ext_file_new) # ,path_style=path_style)


#%%
########################################
######### Define meteo forcing #########
########################################

# generate old format external forcings file (.ext): spatial data
ext_file_old = os.path.join(dir_output_run, f'{model_name}_old.ext')

#%%
ext_old = hcdfm.ExtOldModel()

if meteo_type=='ERA5': # ERA5 - download spatial fields of air pressure, wind speeds and Charnock coefficient
    dir_output_data_era5 = os.path.join(dir_output_bc,'meteo', 'ERA5')
    os.makedirs(dir_output_data_era5, exist_ok=True)
        
    varlist_list = [['msl','u10n','v10n','chnk']]
    for varlist in varlist_list:
        for varkey in varlist:
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
    spw_file = spw_name + '.spw'
    shutil.copyfile(spw_file_origin, os.path.join(dir_output_run,spw_file))

    uniformwind_filename = 'uniformwind0.wnd'
    shutil.copyfile(os.path.join(dir_output_bc,'meteo',uniformwind_filename),os.path.join(dir_output_run,uniformwind_filename))

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

elif meteo_type == 'spiderweb+ERA5': # does not work so far
    # ERA5
    dir_output_data_era5 = os.path.join(dir_output_bc,'meteo', 'ERA5')
    os.makedirs(dir_output_data_era5, exist_ok=True)
    
    varlist_list = [['msl','u10n','v10n']]
    for varlist in varlist_list:
        for varkey in varlist:
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
    
    # spiderweb
    spw_file = spw_name + '.spw'
    shutil.copyfile(spw_file_origin, os.path.join(dir_output_run,spw_file))
    forcing_spw = hcdfm.ExtOldForcing(quantity='airpressure_windx_windy',
                                        filename=spw_file,
                                        filetype=hcdfm.ExtOldFileType.SpiderWebData, 
                                        method=hcdfm.ExtOldMethod.PassThrough, 
                                        operand=hcdfm.Operand.override)
    ext_old.forcing.append(forcing_spw)

ext_old.save(filepath=ext_file_old) # , path_style=path_style)

#%%
if 'ERA5' in meteo_type:
    # plot converted ERA5 data
    file_era5 = os.path.join(dir_output_bc,'meteo','ERA5',f'*{date_min[:4]}*.nc')
    ds_era5 = xr.open_mfdataset(file_era5)
    ds_era5

    # plot
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
    ds_era5.u10n.isel(time=0).plot(ax=ax1)
    dfmt.plot_coastlines(ax=ax1, crs=crs)
    ds_era5.v10n.isel(time=0).plot(ax=ax2)
    dfmt.plot_coastlines(ax=ax2, crs=crs)
    fig.tight_layout()


#%%
#################################################
############## Generate obs file ################
#################################################
# The D-Flow FM model wil have mapoutput and hisoutput. 
# A file with coordinates of obs stations will be generated.

# Read shp file of points along the MZB coastline
gdfp = gpd.read_file(r'p:\11210471-001-compass\01_Data\Coastal_boundary\points\coastal_bnd_MZB_5mMSL_points_1km.shp')

#%%
# crop output points to the area where output is needed for the flood model
gdfp = gdfp.cx[output_bbox[0]:output_bbox[2],output_bbox[1]:output_bbox[3]]

#%%
# Convert points to the xyn file format
xcor = gdfp['geometry'][:].x; xcor.name = 'x'
ycor = gdfp['geometry'][:].y; ycor.name = 'y'
tmp = pd.concat([xcor,ycor],axis=1)
tmp = tmp.dropna()
tmp['names'] = tmp.index

#%%
# Load extra points that are used for verification
tmp2 = pd.read_table(r'p:\11210471-001-compass\02_Models\Delft3DFM\mozambique_model\geometry\output_locations\MZB_Sofala_IHO_obs.xyn',sep=" ",names=['x', 'y', 'names'])
pd_obs = tmp2._append(tmp,ignore_index=True); del tmp, tmp2

# save obs points
file_obs = os.path.join(dir_output_run, f'{model_name}_obs.xyn')
pd_obs.to_csv(file_obs, sep=' ', header=False, index=False, float_format='%.6f')

# plot obs points
fig, ax = plt.subplots(figsize=(8,4))
xu_grid_uds.grid.plot(ax=ax,linewidth=0.5,color='k',alpha=0.2)
ax.plot(pd_obs['x'],pd_obs['y'],'rx')
dfmt.plot_coastlines(ax=ax, crs=crs)


#%%
############################################
############ Generate mdu file #############
############################################
# In order for the model to run, we need a model definition file, i.e., a *.mdu file

# initialize mdu file and update settings
mdu_file = os.path.join(dir_output_run, f'{model_name}.mdu')

#mdu = hcdfm.FMModel()

# use mdu file from GTSM
base_mdu = os.path.join(dir_output_main,'general','base_model_settings.mdu')
mdu = hcdfm.FMModel(base_mdu)

if os.path.exists(pathfile_illegalcells):
    mdu.geometry.drypointsfile = pathfile_illegalcells

# add the grid (_net.nc, network file)
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


#%%
#######################################################
############# Generate DIMR and bat file ##############
#######################################################
# In order to run the model via DIMR we need a `dimr_config.xml` file. 
# If you are running this notebook on a Windows platform, a *.bat file will also be created

nproc = 4 # number of processes
dimrset_folder = r"p:\d-hydro\dimrset\weekly\2.25.17.78708" # alternatively r"c:\Program Files\Deltares\Delft3D FM Suite 2023.03 HMWQ\plugins\DeltaShell.Dimr\kernels" #alternatively r"p:\d-hydro\dimrset\weekly\2.25.17.78708"
dfmt.create_model_exec_files(file_mdu=mdu_file, nproc=nproc, dimrset_folder=dimrset_folder)

#%%
# add a batch file to submit job to h7 cluster (as al alternative to running it on your node or computer)
batchfile_h7 = r'p:\11210471-001-compass\02_Models\Delft3DFM\mozambique_model\scripts\submit_singularity_h7.sh'

pathfile_h7 = os.path.join(dir_output_run,'submit_singularity_h7.sh')

replacements = {'JOBNAME':'MZB', 'MDUFILE':f'{model_name}.mdu'}

with open(batchfile_h7) as infile, open(pathfile_h7, 'w',newline='\n') as outfile:
    for line in infile:
        for src, target in replacements.items():
            line = line.replace(src, target)
        outfile.write(line)

#%%
###############################################
############ Visualize model tree #############
###############################################

# visualize the model tree, show_tree is available for all HYDROLIB-core model components
mdu_obj = hcdfm.FMModel(mdu_file)
mdu_obj.show_tree()
# %%

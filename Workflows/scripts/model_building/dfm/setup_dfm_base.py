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
    bbox_dfm = snakemake.params.dfm_bbox
    dfm_res = snakemake.wildcards.dfm_res
    bathy = snakemake.wildcards.bathy
    tidemodel = snakemake.wildcards.tidemodel
    model_name = snakemake.params.model_name
    dir_output_main = snakemake.params.dir_model
    path_data_cat = snakemake.params.data_cat
    # grid_network = snakemake.output.grid_network
    # illigalcells = snakemake.output.illigalcells
    # pli_file = snakemake.output.pli_file
    # ext_file_new = snakemake.output.ext_file_new
else:
    bbox_dfm = "[32.3,42.5,-27.4,-9.5]"
    dfm_res = "450"
    bathy = "gebco2024"
    tidemodel = 'GTSMv4.1_opendap' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap
    model_name = f'base_{dfm_res}_{bathy}_{tidemodel}'
    dir_output_main = f'p:/11210471-001-compass/02_Models/mozambique/dfm/'
    path_data_cat = os.path.abspath("../../../data_catalogs/datacatalog_general.yml")
    # grid_network = os.path.join(dir_output_run, 'grid_network.nc')
    # illigalcells = os.path.join(dir_output_run, 'pli_file.pli')
    # pli_file = os.path.join(dir_output_run,"illegalcells.pol")
    # ext_file_new =

#%%
# Define hydromt datacatalog
data_catalog = hydromt.DataCatalog(data_libs = [path_data_cat])

# needed to define wind_forcing and bathy paths
#%%
# define model name, general settings and output location
dir_output_run = os.path.join(dir_output_main, model_name)

# make directories, if not yet present
os.makedirs(dir_output_run, exist_ok=True)

generate_grid = False # option to skip grid generation if this was already done.
overwrite = False # used for downloading of forcing data. Always set to True when changing the domain
crs = 'EPSG:4326' # coordinate reference system

#%%
# domain and resolution
bbox_list = [float(x) for x in bbox_dfm.strip("[]").split(",")]
lon_min, lon_max, lat_min, lat_max = bbox_list

# dxy is the base grid resolution - i.e. the coarsest grid size in your grid. It is in degrees.
dxy = 0.2 # degrees
# Defines the minimum length of the edges of the computational mesh cells.
min_edge_size = dfm_res # m 

#%%##########################################
###### Grid generation and refinement #######
#############################################

# File name for the grid network
netfile = os.path.join(dir_output_run, 'grid_network.nc')
poly_file = os.path.join(dir_output_run, 'pli_file.pli')
pathfile_illegalcells = os.path.join(dir_output_run,"illegalcells.pol")

# Load grid if it was already generated
if os.path.exists(netfile):
    # Load the grid if the file does not exist
    xu_grid_uds = dfmt.open_partitioned_dataset(netfile)

else:
    # generate spherical regular grid
    mk_object = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy, crs=crs)

    # generate plifile from grid extent and coastlines
    bnd_gdf = dfmt.generate_bndpli_cutland(mk=mk_object, res='h', buffer=0.2)
    bnd_gdf_interp = dfmt.interpolate_bndpli(bnd_gdf, res=0.03)

    # filter out boundary sections that are very short (<5km in this case)
    ids = []
    for ii in range(len(bnd_gdf_interp.index)):
        if bnd_gdf_interp.iloc[ii].values.to_crs('EPSG:3857').length < 5000: # approx. 111 km in 1 degree
            ids.append(ii)

    bnd_gdf_interp = bnd_gdf_interp.drop(ids)

    pli_polyfile = dfmt.geodataframe_to_PolyFile(bnd_gdf_interp)
    pli_polyfile.save(poly_file)

    # plot basegrid and polyline
    fig, ax = plt.subplots()
    mk_object.mesh2d_get().plot_edges(ax,zorder=1)
    bnd_gdf_interp.plot(ax=ax, edgecolor='r')
    dfmt.plot_coastlines(ax=ax, crs=crs)

    # Define bathymetry
    file_nc_bathy_sel = data_catalog[bathy].path
    data_bathy_sel = xr.open_dataset(file_nc_bathy_sel).elevation # what does elevation do?
        
    data_bathy_sel.load()

    # subset to area of interest
    # data_bathy_sel = data_bathy.sel(lon=slice(lon_min-1, lon_max+1), lat=slice(lat_min-1, lat_max+1))

    # refine grid
    dfmt.refine_basegrid(mk=mk_object, data_bathy_sel=data_bathy_sel, min_edge_size=min_edge_size)

    # plot
    fig, ax = plt.subplots()
    mk_object.mesh2d_get().plot_edges(ax,zorder=1)
    dfmt.plot_coastlines(ax=ax, crs=crs)

    # remove land with GSHHS coastlines
    dfmt.meshkernel_delete_withcoastlines(mk=mk_object,res='i',min_area=50,crs=crs)

    # plot
    fig, ax = plt.subplots(figsize=(12,7))
    mk_object.mesh2d_get().plot_edges(ax,zorder=1)
    dfmt.plot_coastlines(ax=ax,res='i',min_area=50,crs=crs)

    # derive illegalcells geodataframe
    illegalcells_gdf = dfmt.meshkernel_get_illegalcells(mk=mk_object)
    # create and add drypointsfile if there are any cells generated that will result in high orthogonality
    if len(illegalcells_gdf) > 0:
        illegalcells_polyfile = dfmt.geodataframe_to_PolyFile(illegalcells_gdf)
        illegalcells_polyfile.save(pathfile_illegalcells)

    # convert to xugrid
    xu_grid_uds = dfmt.meshkernel_to_UgridDataset(mk=mk_object, crs=crs)

    # interpolate bathymetry onto the grid
    data_bathy_interp = data_bathy_sel.interp(lon=xu_grid_uds.obj.mesh2d_node_x, lat=xu_grid_uds.obj.mesh2d_node_y)
    xu_grid_uds['mesh2d_node_z'] = data_bathy_interp.clip(max=10)

    # write xugrid grid to netcdf
    xu_grid_uds.ugrid.to_netcdf(netfile)

# plot bathymetry and grid
fig, ax = plt.subplots(figsize=(8,4))
xu_grid_uds.mesh2d_node_z.ugrid.plot(ax=ax,center=False)
xu_grid_uds.grid.plot(ax=ax,linewidth=0.5,color='white',alpha=0.2)
# ctx.add_basemap(ax=ax, crs=crs, attribution=False)
dfmt.plot_coastlines(ax=ax, crs=crs)


#%%##################################################
### Generate boundary conditions from tidal model ###
#####################################################

# generate new format external forcings file (.ext): initial and open boundary condition
ext_file_new = os.path.join(dir_output_run, 'ext_file_new.ext')
ext_new = hcdfm.ExtModel()

# interpolate tidal components to boundary conditions file (.bc)
dfmt.interpolate_tide_to_bc(ext_new=ext_new, tidemodel=tidemodel, file_pli=poly_file, component_list=None)

#save new ext file
ext_new.save(filepath=ext_file_new) # ,path_style=path_style)

# %%
# This script is based on the ModelBuilder template from dfm_tools v0.31.0 (accessed on 20/01/2025). 
# Adapted by Natalia Aleksandrova & Doris Vertegaal
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
from hydromt import DataCatalog
import numpy as np

#%%
if "snakemake" in locals():
    dfm_res = float(snakemake.wildcards.dfm_res)
    dxy_base = float(snakemake.params.dfm_dxy_base)
    bathy = snakemake.wildcards.bathy
    tidemodel = snakemake.wildcards.tidemodel
    dir_output_main = os.path.abspath(snakemake.output.dir_model)
    bbox_dfm = snakemake.params.dfm_bbox
    path_data_cat = os.path.abspath(snakemake.params.data_cat)
    CF_value = float(snakemake.wildcards.CF_SLR)
    CF_value_txt = snakemake.wildcards.CF_SLR
else:
    dfm_res_txt = "450"
    dfm_res = 450 # m
    dxy_base = 0.02 # degrees
    bathy = "gebco2024_MZB"
    tidemodel = 'FES2014' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap, tpxo80_opendap
    CF_value = -0.14
    CF_value_txt = "0.14"
    dir_output_main = f'p:/11210471-001-compass/02_Models/sofala/Idai/dfm/base_{dfm_res_txt}_{bathy}_{tidemodel}_CF{CF_value_txt}'
    bbox_dfm = "[32.3,42.5,-27.4,-9.5]"
    path_data_cat = os.path.abspath("../../../03_data_catalogs/datacatalog_general.yml")

# Correct for the missing . in the snake that snakemake cannot read
if tidemodel == "GTSMv41opendap":
    tidemodel = "GTSMv4.1_opendap"
else:
    pass

if tidemodel == "GTSMv41":
    tidemodel = "GTSMv4.1"
else:
    pass
#%%
# Define hydromt datacatalog
data_catalog = DataCatalog(data_libs = [path_data_cat])

# needed to define wind_forcing and bathy paths
#%%
# define model name, general settings and output location
dir_output_run = os.path.join(dir_output_main)
dir_output_geometry = os.path.join(dir_output_main, "geometry")

# make directories, if not yet present
os.makedirs(dir_output_run, exist_ok=True)
os.makedirs(dir_output_geometry, exist_ok=True)

crs = 'EPSG:4326' # coordinate reference system

#%%
# domain and resolution
bbox_list = [float(x) for x in bbox_dfm.strip("[]").split(",")]
lon_min, lon_max, lat_min, lat_max = bbox_list
print(bbox_dfm)

# dxy is the base grid resolution - i.e. the coarsest grid size in your grid. It is in degrees.
dxy = dxy_base # degrees
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
    bnd_gdf = dfmt.generate_bndpli_cutland(mk=mk_object, res='h', buffer=0.05)
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
    # Save the figure
    output_path = os.path.join(dir_output_geometry, 'basegrid_polyline.png')
    fig.savefig(output_path, dpi=300, bbox_inches='tight')

    # Define bathymetry
    data_bathy_sel = data_catalog.get_rasterdataset(bathy)
    data_bathy_sel = data_bathy_sel.rename({'x':'lon','y':'lat'})
    data_bathy_sel = data_bathy_sel.sortby("lat")
    data_bathy_sel = data_bathy_sel.sortby("lon")
    data_bathy_sel = data_bathy_sel.sel(lon=slice(lon_min-1,lon_max+1),lat=slice(lat_min-1,lat_max+1), drop=True)
    data_bathy_sel.load()
    if '_FillValue' in data_bathy_sel.attrs:
        data_bathy_sel_masked = xr.where(data_bathy_sel == data_bathy_sel.attrs['_FillValue'], 0, data_bathy_sel)
        data_bathy_sel = data_bathy_sel_masked
    data_bathy_sel = data_bathy_sel.fillna(0)

    # subset to area of interest
    # data_bathy_sel = data_bathy.sel(lon=slice(lon_min-1, lon_max+1), lat=slice(lat_min-1, lat_max+1))

    # refine grid
    dfmt.refine_basegrid(mk=mk_object, data_bathy_sel=data_bathy_sel, min_edge_size=min_edge_size)

    # plot
    fig, ax = plt.subplots()
    mk_object.mesh2d_get().plot_edges(ax,zorder=1)
    dfmt.plot_coastlines(ax=ax, crs=crs)
    # Save the figure
    output_path = os.path.join(dir_output_geometry, 'basegrid_refined.png')
    fig.savefig(output_path, dpi=300, bbox_inches='tight')

    # remove land with GSHHS coastlines
    dfmt.meshkernel_delete_withcoastlines(mk=mk_object,res='i',min_area=50,crs=crs)

    # plot
    fig, ax = plt.subplots(figsize=(12,7))
    mk_object.mesh2d_get().plot_edges(ax,zorder=1)
    dfmt.plot_coastlines(ax=ax,res='i',min_area=50,crs=crs)
    # Save the figure
    output_path = os.path.join(dir_output_geometry, 'basegrid_refined_noland.png')
    fig.savefig(output_path, dpi=300, bbox_inches='tight')

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

# Save the figure
output_path = os.path.join(dir_output_geometry, 'grid_refined_bathy.png')
fig.savefig(output_path, dpi=300, bbox_inches='tight')


#%%##################################################
### Generate boundary conditions from tidal model ###
#####################################################

# Generate new format external forcings file (.ext): initial and open boundary condition
ext_new_path = os.path.join(dir_output_run, 'ext_file_new.ext')
ext_new = hcdfm.ExtModel()

# Interpolate tidal components to boundary conditions file (.bc)
dfmt.interpolate_tide_to_bc(ext_new=ext_new, tidemodel=tidemodel, file_pli=poly_file, component_list=None)

# Save new ext file
ext_new.save(filepath=ext_new_path)

# Make the paths relative
with open(ext_new_path, 'r') as file:
    lines = file.readlines()
# Modify the lines that contain file paths
modified_lines = []
for line in lines:
    if 'locationFile' in line or 'forcingFile' in line:
        # Split the line by the first '=' and get the key and path
        key, path = line.split('=', 1)

        # Remove leading/trailing spaces from the path and get just the file name
        path = path.strip()
        file_name = os.path.basename(path)
        # Replace the path with just the file name and retain the comment
        modified_line = f'{key.strip()} = {file_name}\n'
        modified_lines.append(modified_line)
    else:
        modified_lines.append(line)

# Write the modified lines to the output file
with open(ext_new_path, 'w') as file:
    file.writelines(modified_lines)
print(f'Modified file saved to: {ext_new_path}')

# %%
# Add the A0 constituent to include SLR in the tidal boundary
if CF_value != 0:
    file_bc = os.path.join(dir_output_main, f"tide_{tidemodel}_pli_file.bc")
    forcingmodel_object = hcdfm.ForcingModel(file_bc)

    for forcing in forcingmodel_object.forcing:
        if forcing.datablock:
            last_entry = forcing.datablock[-1]  # Get the last row
            
            # Check if "A0" already exists and update it if so
            if last_entry[0] == "A0":
                last_entry[1] = CF_value  # Update CF_value
            else:
                forcing.datablock.append(["A0", CF_value, 0])  # Append if A0 doesn't exist

    # Save back to the same file
    forcingmodel_object.save(file_bc)

else:
    pass

# %%

# -*- coding: utf-8 -*-
"""
Created on Tue May  2 16:02:14 2023

@author: zijlker
"""
from pathlib import Path
import contextily as ctx
import geopandas
import pandas as pd
import geopandas as gpd
import numpy as np
import pyproj
import xarray as xr
import dfm_tools as dfmt
# from global_land_mask import globe
# import cartopy.io.shapereader as shpreader
import shapely
from shapely.geometry import LineString, MultiPoint, Point
# import cartopy.crs as ccrs
from scipy.spatial import KDTree
import hydrolib.core.dflowfm as hcdfm

def crs_conversion(df, inProj, outProj):
    # Create the transformer object from the source and target CRS
     transformer = pyproj.Transformer.from_crs(inProj, outProj, always_xy=True)
    
     # Transform the coordinates using the transformer object
     transformed_coords = [transformer.transform(x, y) for x, y in zip(df['x'], df['y'])]
    
     # Create a new DataFrame with the transformed coordinates
     transformed_df = pd.DataFrame(transformed_coords, columns=['x', 'y'])
    
     return transformed_df

def generate_sfincs_bnd_from_delft3d_xyn(obs_file_input,
                                         obs_file_input_crs,
                                         bnd_file_output,
                                         bnd_file_output_crs):
    obs = hcdfm.XYNModel(obs_file_input)
    gdf = dfmt.pointlike_to_geodataframe_points(obs)
    gdf['x'] = gdf.geometry.x
    gdf['y'] = gdf.geometry.y 
    gdf = gdf[['x','y']]
    
    #Write .bnd file
    bnd_file_transf = crs_conversion(gdf,obs_file_input_crs, bnd_file_output_crs)
    bnd_file_transf.to_csv(bnd_file_output, sep = ' ', index=False, header =False)
    
    
def generate_delft3d_obsfile_from_sfincs_bnd(bnd_file_input,
                                             bnd_file_input_crs,
                                             obsfile_output,
                                             obsfile_output_crs):
    #Load sfincs bnd file 
    bnd_file = pd.read_csv(bnd_file_input, sep = ' ', names = ['x', 'y'])
    
    #Transform coordinates
    bnd_file_transf = crs_conversion(bnd_file,bnd_file_input_crs, obsfile_output_crs)
    
    #Add column with names to for DFLOWFM model
    names = [f'sfincs_{i:03d}' for i in range(len(bnd_file_transf))]
    bnd_file_transf['name'] = names
    
    # Create obs file for dflowfm
    obs_sfincs = hcdfm.XYNModel()
    points = [hcdfm.XYNPoint(x=x, y=y, n=f'sfincs_bnd_{i:03d}') \
        for i, (x, y) in enumerate(zip(bnd_file_transf['x'].values, 
                                       bnd_file_transf['y'].values))]
    obs_sfincs.points = points
    obs_sfincs.save(obsfile_output)

    #Write observations file
    # bnd_file_transf.to_csv(obsfile_output, sep = ' ', index=False, header =False)
    
def generate_sfincs_bzs_from_dflowfm_his(reftime, 
                                         obsfile_dflowfm_sfincs_stations,
                                         hisfile_dflowfm,
                                         bzs_outputfile):
    
    his = xr.open_mfdataset(hisfile_dflowfm, preprocess=dfmt.preprocess_hisnc)
    
    obsfile = pd.read_csv(obsfile_dflowfm_sfincs_stations, sep = ' ', names = ['x', 'y', 'name'])
    
    bzs = his['waterlevel'].drop(['station_x_coordinate', 'station_y_coordinate'])\
                           .to_dataframe()\
                           .unstack('stations').loc[:, ('waterlevel',obsfile['name'].values)]
    
    bzs.index = (bzs.index - reftime).total_seconds()
    bzs.to_csv(bzs_outputfile, sep = ' ', header = False, index=True)

def generate_sfincs_bzs_from_dflowfm_his_humber(reftime, 
                                         obsfile_dflowfm_sfincs_stations,
                                         hisfile_dflowfm,
                                         bzs_outputfile):
    
    his = xr.open_mfdataset(hisfile_dflowfm, preprocess=dfmt.preprocess_hisnc)
    
    obsfile = pd.read_csv(obsfile_dflowfm_sfincs_stations, sep = ' ', names = ['x', 'y', 'name'])
    
    bzs = his['waterlevel'].drop(['station_x_coordinate', 'station_y_coordinate'])\
                           .to_dataframe()\
                           .unstack('stations')
    # print(bzs)
    bzs.loc[:, ('waterlevel','sfincs_000')] = bzs.loc[:, ('waterlevel','sfincs_001')]
    bzs.loc[:, ('waterlevel','sfincs_009')] = bzs.loc[:, ('waterlevel','sfincs_010')]
    
    bzs = bzs.loc[:, ('waterlevel',obsfile['name'].values)]
    print(bzs)
    bzs.index = (bzs.index - reftime).total_seconds()
    bzs.to_csv(bzs_outputfile, sep = ' ', header = False, index=True)
    
def generate_coastline_obspoints(lon_min : float, # degrees
                                     lon_max : float, # degrees
                                     lat_min : float, # degrees
                                     lat_max : float, # degrees
                                     interval : float,# degrees
                                     resolution : str ,# degrees
                                     threshold_mindepth : float, # meters
                                     file_nc : Path, # path to netcdf file
                                     output_file : Path): # path to output file,):

    #%% create obs file and snap to grid 

    # Points along coastline with dfm_tools
    coastline = dfmt.get_coastlines_gdb(res = resolution,
                             bbox = [lon_min, lat_min, lon_max, lat_max],
                             crs = "EPSG:4326")
    linestrings = [polygon.boundary for polygon in coastline['geometry']]
    linstrings_cropped = [shapely.ops.clip_by_rect(linestring, lon_min, lat_min, lon_max, lat_max) for linestring in linestrings]
    shp_clipped = gpd.GeoDataFrame(geometry = gpd.GeoSeries(linstrings_cropped))
    shp_clipped = shp_clipped[~shp_clipped.is_empty]
    
    #Build GeoSeries with points along coastline 
    cpoints = gpd.GeoSeries()
    for index, line in shp_clipped.iterrows():
        if 'MULTILINESTRING' in str(line.values):       
            shapes = str(line.values[0]).strip().split("\n")       
            gdf = gpd.GeoDataFrame({'geometry': shapes})
            gdf['geometry'] = gpd.GeoSeries.from_wkt(gdf.geometry)
            gdf = gdf.set_geometry('geometry').explode()
            for iindex, iline in  gdf.iterrows():
                distances = np.arange(0, iline.geometry.length, interval)
                points = MultiPoint([LineString(iline.geometry).interpolate(distance) for distance in distances])
                gs =gpd.GeoSeries(Point(pnt.x,pnt.y) for pnt in points.geoms)
                cpoints=pd.concat([cpoints, gs])
        else:
            distances = np.arange(0, line.geometry.length, interval)
            points = MultiPoint([LineString(line.geometry).interpolate(distance) for distance in distances])
            gs =gpd.GeoSeries(Point(pnt.x,pnt.y) for pnt in points.geoms)
            cpoints=pd.concat([cpoints, gs])
    locs_coast = gpd.GeoDataFrame(cpoints, geometry=cpoints.geometry, crs="EPSG:4326")

    # Snap to model points and -5 m depth ## TO DO
    obs=pd.DataFrame()
    obs['x'] = locs_coast.geometry.x.values
    obs['y'] = locs_coast.geometry.y.values

    # #Create obs file for dflowfm
    # obs_sfincs = hcdfm.XYNModel()
    # points = [hcdfm.XYNPoint(x=x, y=y, n=f'sfincs_bnd_{i:03d}') \
    #     for i, (x, y) in enumerate(zip(obs['x'].values, 
    #                                    obs['y'].values))]
    # obs_sfincs.points = points
    # obs_sfincs.save('sealevel_not_snapped_obs.xyn')

    # retrieve modeldata

    ds_net = xr.open_dataset(file_nc)

    #%% calculate face center x, y and z
    # net_face_x = []
    # net_face_y = []
    # net_face_z = []
    # for face in ds_net.mesh2d_nFaces:
    #     idx_node = ds_net.mesh2d_face_nodes\
    #         .sel(mesh2d_nFaces = face)\
    #         .dropna(dim='mesh2d_nMax_face_nodes').data - 1
    #     net_face_x_sel = ds_net.mesh2d_node_x.sel(mesh2d_nNodes = idx_node).data.mean()
    #     net_face_y_sel = ds_net.mesh2d_node_y.sel(mesh2d_nNodes = idx_node).data.mean()
    #     net_face_z_sel = ds_net.mesh2d_node_z.sel(mesh2d_nNodes = idx_node).data.mean()
    #     net_face_x.append(net_face_x_sel)
    #     net_face_y.append(net_face_y_sel)
    #     net_face_z.append(net_face_z_sel)


    #Make dictionary with cell centers
    net_faces = pd.DataFrame({'x':ds_net['mesh2d_face_x'].data,
                              'y':ds_net['mesh2d_face_y'].data,
                              'z':ds_net['mesh2d_flowelem_bl'].data})

    #Select only cells with z < threshold_mindepth
    bool_valid_cells = (net_faces['z']<-threshold_mindepth)

    #creating kdtree with valid cell centers (cartesian coordinates)
    def xlonylat2xyzcartesian(data):
        """
        necessary to calculate cartesian distances, otherwise nearest neigbour can fail.
        https://stackoverflow.com/questions/45127141/find-the-nearest-point-in-distance-for-all-the-points-in-the-dataset-python
        """
        R = 6367
        phi = np.deg2rad(data['y'])
        theta = np.deg2rad(data['x'])
        data = pd.DataFrame()
        data['x_cart'] = R * np.cos(phi) * np.cos(theta)
        data['y_cart'] = R * np.cos(phi) * np.sin(theta)
        data['z_cart'] = R * np.sin(phi)
        return data

    data_celcenxy_valid = net_faces.loc[bool_valid_cells, ['x', 'y']].reset_index() #pd.DataFrame({'x':net_node_x[bool_valid_cells],'y':net_node_y[bool_valid_cells]})#,'area':data_cellarea[bool_valid_cells]})
    data_celcenxy_valid_cart = xlonylat2xyzcartesian(data_celcenxy_valid)
    tree = KDTree(data_celcenxy_valid_cart[['x_cart','y_cart','z_cart']])

    def dist_to_arclength(chord_length):
        """
        https://stackoverflow.com/questions/45127141/find-the-nearest-point-in-distance-for-all-the-points-in-the-dataset-python
        """
        R = 6367 # earth radius
        central_angle = 2*np.arcsin(chord_length/(2.0*R)) 
        arclength = R*central_angle*1000
        return arclength

    #finding nearest cellcenter-neighbors of each obspoint in file
    data_obsorg_cart = xlonylat2xyzcartesian(obs)
    distance_cart, index = tree.query(data_obsorg_cart, k=1)
    data_celcenxy_validsel = data_celcenxy_valid.loc[index,:]

    # write 
    obs_snapped=pd.DataFrame()
    obs_snapped['x'] = data_celcenxy_validsel['x'].values
    obs_snapped['y'] = data_celcenxy_validsel['y'].values
    obs_snapped = obs_snapped.drop_duplicates()

    #Create obs file for dflowfm
    obs_sfincs = hcdfm.XYNModel()
    points = [hcdfm.XYNPoint(x=x, y=y, n=f'sfincs_bnd_{i:03d}') for i, (x, y) in enumerate(zip(obs_snapped['x'].values, obs_snapped['y'].values))]
    obs_sfincs.points = points
    obs_sfincs.save(output_file)
    

if __name__=='__main__':
    # bnd_file_input = Path(r'p:\11208614-de-370a\01_models\Humber\sfincs\sfincs_humber\sfincs.bnd')
    # bnd_file_input_crs = 'epsg:32630' 
    # obsfile_output = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\Humber_DFLOWFM_cutout_DCSM\sfincs_obs.xyn')
    # obsfile_output_crs = 'epsg:4326'
    # generate_delft3d_obsfile_from_sinfcs_bnd(bnd_file_input,
    #                                          bnd_file_input_crs,
    #                                          obsfile_output,
    #                                          obsfile_output_crs)
    
    # bnd_file_input = Path(r'p:\11208614-de-370a\01_models\Basque\sfincs\sfincs_basque\sfincs.bnd')
    # bnd_file_input_crs = 'epsg:32630'
    # obsfile_output_crs = 'epsg:4326'
    # generate_delft3d_obsfile_from_sinfcs_bnd(bnd_file_input,
    #                                          bnd_file_input_crs,
    #                                          obsfile_output,
    #                                          obsfile_output_crs)
    
    
    # obsfile_output = Path(r'p:\11208614-de-370a\01_models\Basque\DFLOWFM\Models\Humber_Model_Copy\sfincs_basque_moved_obs.xyn')
    # bzs_outputfile = Path(r'p:\11208614-de-370a\01_models\Basque\sfincs\DFLOWFM_seaboundary\basque.bzs')
    # tref = pd.Timestamp('2021-11-25')
    # hisfile_dflowfm = Path(r'p:\11208614-de-370a\01_models\Basque\DFLOWFM\Models\Humber_Model_Copy\output\humber_0000_his.nc')                                         
    # generate_sfincs_bzs_from_dflowfm_his(tref,
    #                                      obsfile_output,
    #                                      hisfile_dflowfm,
    #                                      bzs_outputfile)
    
    # obsfile_output = Path(r'p:\11208614-de-370a\01_models\Reunion\DFLOWFM\Models\reunion_custom_grid2\reunion_SFINCS_obs.xyn')
    # bzs_outputfile = Path(r'p:\11208614-de-370a\01_models\Reunion\DFLOWFM\Models\reunion_custom_grid2\sfincs\reunion.bzs')
    # tref = pd.Timestamp('2018-04-01')
    # hisfile_dflowfm = Path(r'p:\11208614-de-370a\01_models\Reunion\DFLOWFM\Models\reunion_custom_grid2\output\reunion_his.nc')                                         
    # generate_sfincs_bzs_from_dflowfm_his(tref,
    #                                      obsfile_output,
    #                                      hisfile_dflowfm,
    #                                      bzs_outputfile)
    
    # obsfile_output = Path(r'p:\11208614-de-370a\01_models\BES\DFLOWFM\Models\bonaire\bonaire_SFINCS_obs.xyn')
    # bzs_outputfile = Path(r'p:\11208614-de-370a\01_models\BES\DFLOWFM\Models\bonaire\sfincs\bonaire.bzs')
    # tref = pd.Timestamp('2016-09-25')
    # hisfile_dflowfm = Path(r'p:\11208614-de-370a\01_models\BES\DFLOWFM\Models\bonaire\output\bonaire_his.nc')                                         
    # generate_sfincs_bzs_from_dflowfm_his(tref,
    #                                      obsfile_output,
    #                                      hisfile_dflowfm,
    #                                      bzs_outputfile)
    
    # obsfile_output = Path(r'p:\11208614-de-370a\01_models\BES\DFLOWFM\Models\saba_st_eustatius\saba_SFINCS_obs.xyn')
    # bzs_outputfile = Path(r'p:\11208614-de-370a\01_models\BES\DFLOWFM\Models\saba_st_eustatius\sfincs\saba.bzs')
    # tref = pd.Timestamp('2017-09-01')
    # hisfile_dflowfm = Path(r'p:\11208614-de-370a\01_models\BES\DFLOWFM\Models\saba_st_eustatius\output\saba_st_eustatius_his.nc')                                         
    # generate_sfincs_bzs_from_dflowfm_his(tref,
    #                                      obsfile_output,
    #                                      hisfile_dflowfm,
    #                                      bzs_outputfile)
    
    # obsfile_output = Path(r'p:\11208614-de-370a\01_models\BES\DFLOWFM\Models\saba_st_eustatius\st_eustatius_SFINCS_obs.xyn')
    # bzs_outputfile = Path(r'p:\11208614-de-370a\01_models\BES\DFLOWFM\Models\saba_st_eustatius\sfincs\st_eustatius.bzs')
    # tref = pd.Timestamp('2017-09-01')
    # hisfile_dflowfm = Path(r'p:\11208614-de-370a\01_models\BES\DFLOWFM\Models\saba_st_eustatius\output\saba_st_eustatius_his.nc')                                         
    # generate_sfincs_bzs_from_dflowfm_his(tref,
    #                                      obsfile_output,
    #                                      hisfile_dflowfm,
    #                                      bzs_outputfile)

    # obsfile_output = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm_2013\geometry\output_locations\sfincs_obs.xyn')
    # bzs_outputfile = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm_2013\output\sfincs\humber.bzs')
    # tref = pd.Timestamp('2013-11-13')
    # hisfile_dflowfm = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm_2013\output\DCSM-FM_0_5nm_0000_his.nc')                                         
    # generate_sfincs_bzs_from_dflowfm_his(tref,
    #                                      obsfile_output,
    #                                      hisfile_dflowfm,
    #                                      bzs_outputfile)
    
    # obsfile_output = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm_2007\geometry\output_locations\sfincs_obs.xyn')
    # bzs_outputfile = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm_2007\output\sfincs\humber.bzs')
    # tref = pd.Timestamp('2006-01-01')
    # hisfile_dflowfm = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm_2007\output\DCSM-FM_0_5nm_0000_his.nc')                                         
    # generate_sfincs_bzs_from_dflowfm_his(tref,
    #                                      obsfile_output,
    #                                      hisfile_dflowfm,
    #                                      bzs_outputfile)
    
    obsfile_output = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\geometry\output_locations\sfincs_obs.xyn')
    bzs_outputfile = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\computations\Xavier_ERA5\output\sfincs\humber.bzs')
    tref = pd.Timestamp('2013-11-13')
    hisfile_dflowfm = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\computations\Xavier_ERA5\output\DCSM-FM_0_5nm_0000_his.nc')                                         
    generate_sfincs_bzs_from_dflowfm_his_humber(tref,
                                         obsfile_output,
                                         hisfile_dflowfm,
                                         bzs_outputfile)
    
    obsfile_output = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\geometry\output_locations\sfincs_obs.xyn')
    bzs_outputfile = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\computations\Xavier_IFS\output\sfincs\humber.bzs')
    tref = pd.Timestamp('2013-11-13')
    hisfile_dflowfm = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\computations\Xavier_IFS\output\DCSM-FM_0_5nm_0000_his.nc')                                         
    generate_sfincs_bzs_from_dflowfm_his_humber(tref,
                                         obsfile_output,
                                         hisfile_dflowfm,
                                         bzs_outputfile)

    obsfile_output = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\geometry\output_locations\sfincs_obs.xyn')
    bzs_outputfile = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\computations\Chiara_IFS\output\sfincs\humber.bzs')
    tref = pd.Timestamp('2020-02-05')
    hisfile_dflowfm = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\computations\Chiara_IFS\output\DCSM-FM_0_5nm_0000_his.nc')                                         
    generate_sfincs_bzs_from_dflowfm_his_humber(tref,
                                         obsfile_output,
                                         hisfile_dflowfm,
                                         bzs_outputfile)

    obsfile_output = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\geometry\output_locations\sfincs_obs.xyn')
    bzs_outputfile = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\computations\Chiara_DTExtremes\output\sfincs\humber.bzs')
    tref = pd.Timestamp('2020-02-05')
    hisfile_dflowfm = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\computations\Chiara_DTExtremes\output\DCSM-FM_0_5nm_0000_his.nc')                                         
    generate_sfincs_bzs_from_dflowfm_his_humber(tref,
                                         obsfile_output,
                                         hisfile_dflowfm,
                                         bzs_outputfile)

    obsfile_output = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\geometry\output_locations\sfincs_obs.xyn')
    bzs_outputfile = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\computations\Chiara_ERA5\output\sfincs\humber.bzs')
    tref = pd.Timestamp('2020-02-05')
    hisfile_dflowfm = Path(r'p:\11208614-de-370a\01_models\Humber\DFLOWFM\Models\DCSM-FM_0_5nm\computations\Chiara_ERA5\output\DCSM-FM_0_5nm_0000_his.nc')                                         
    generate_sfincs_bzs_from_dflowfm_his_humber(tref,
                                         obsfile_output,
                                         hisfile_dflowfm,
                                         bzs_outputfile)
    
    
    # bnd_file_input_crs = 'epsg:32630'
    # obsfile_output_crs = 'epsg:4326'
    # generate_sfincs_bnd_from_delft3d_xyn(r'p:\11208614-de-370a\03_modelruns\cosmos\run_folder\scenarios\basque_dec2021_ifs\20211209_00z\models\northeast_atlantic\delft3dfm\delft3dfm_basque\input\sfincs_basque_moved_obs.xyn',
    #                                      obsfile_output_crs,
    #                                         r'p:\11208614-de-370a\03_modelruns\cosmos\run_folder\scenarios\basque_dec2021_ifs\20211209_00z\models\northeast_atlantic\delft3dfm\delft3dfm_basque\input\sfincs.bnd',
    #                                         bnd_file_input_crs)
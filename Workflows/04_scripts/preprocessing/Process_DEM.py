#%%
from hydromt_wflow import WflowModel
from hydromt.config import configread
from hydromt.log import setuplog
from os.path import join, basename
import pandas as pd
from hydromt import DataCatalog, flw
import matplotlib.pyplot as plt
import geopandas as gpd
import pyflwdir
import numpy as np
import os
import xarray as xr

def add_hydrography_layers(DEM, output_path, write_flowdir_vector=True):
    """Function that derives Flowdir, uparea, streamorder and basin layers from a DEM
    All layers are written to a tiff file and a basins geojson is created"""

    DEM["flwdir"] = flw.d8_from_dem(
        da_elv=DEM['elevtn'],
        gdf_stream=None,
        max_depth=-1,  # no local pits
        outlets="edges",
        idxs_pit=None)
    dims = DEM.raster.dims
    flwdir = flw.flwdir_from_da(
        DEM["flwdir"], ftype="infer", check_ftype=True)
    if write_flowdir_vector:
        feate = flwdir.streams(min_sto=6)  # minimum stream order for plotting
        gdfe = gpd.GeoDataFrame.from_features(feate, crs=DEM.raster.crs)
        gdfe.to_file(os.path.join(output_path, 'flowdir.geojson'))
    #Derive uparea
    uparea = flwdir.upstream_area(unit="km2")
    DEM["uparea"] = xr.Variable(dims, uparea, attrs=dict(_FillValue=-9999))

    #Streamorder
    strord = flwdir.stream_order()
    DEM["strord"] = xr.Variable(dims, strord)
    DEM["strord"].raster.set_nodata(255)

    # slope
    slope = pyflwdir.dem.slope(
        elevtn=DEM["elevtn"].values,
        nodata=DEM["elevtn"].raster.nodata,
        latlon=False,  # True if geographic crs, False if projected crs
        transform=DEM["elevtn"].raster.transform,
    )
    DEM["slope"] = xr.Variable(dims, slope)
    DEM["slope"].raster.set_nodata(DEM["elevtn"].raster.nodata)

    # basin at the pits locations
    basins = flwdir.basins(idxs=flwdir.idxs_pit).astype(np.int32)
    DEM["basins"] = xr.Variable(dims, basins, attrs=dict(_FillValue=0))

    # basin index file
    gdf_basins = DEM["basins"].raster.vectorize()

    #Export the gridded data as tif files in a new folder
    # export the hydrography data as tif files (one per variable)
    DEM.raster.to_mapstack(
        root= output_path,
        driver="GTiff")
    # export the basin index as geosjon
    gdf_basins.to_file(
        os.path.join(output_path, "DEM_Bonaire_10m_processed_basins.geojson"), driver="GeoJSON"
    )
#%% Create Flowdir map Bonaire

#Get DEM data
data_catalog = DataCatalog(r"../SFINCS/datacatalog_general.yml")
dem = data_catalog.get_rasterdataset(
    "fabdem")
dem = dem.to_dataset()
dem.load()
output_path_dem = r"DEM_processed"



#%%
add_hydrography_layers(dem, output_path_dem)


# %%

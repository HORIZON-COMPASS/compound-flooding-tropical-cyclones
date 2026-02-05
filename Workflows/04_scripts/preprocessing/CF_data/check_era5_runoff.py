# Script to check and process era5 and era5-land run off
#%%
# Importing the necessary packages
import os
import numpy as np
import ast
import xarray as xr
import hydromt
import geopandas as gpd
import copy
import yaml
from datetime import datetime


#%%
tc_name = "Idai"
start_date = "20190309 000000"
end_date = "20190325 060000"
wflow_region = "p:/11210471-001-compass/02_Models/sofala/Idai/wflow/staticgeoms/region.geojson"
data_cat = [
    '../../../03_data_catalogs/datacatalog_general.yml',
    '../../../03_data_catalogs/datacatalog_CF_forcing.yml',
] 
precip_name = "era5_hourly"
CF_value = -7
CF_value_txt = f"{CF_value}"
output_CF_rainfall = f"p:/11210471-001-compass/01_Data/counterfactuals/precipitation/{precip_name}_CF{CF_value}_{tc_name}.nc"
CF_catalog_path = "../../../03_data_catalogs/datacatalog_CF_forcing.yml"

#%%
# Read data catalog
data_catalog = hydromt.data_catalog.DataCatalog(data_libs = data_cat)

# Read the region that needs precipitation input
region = gpd.read_file(wflow_region)

# Convert to datetime objects
start_dt = datetime.strptime(start_date, "%Y%m%d %H%M%S")
end_dt = datetime.strptime(end_date, "%Y%m%d %H%M%S")

# Pass as a time tuple for HydroMT
time_range = (start_dt, end_dt)

#%%
# Load raster data for specified region and time range
precip_data = data_catalog.get_rasterdataset(
    data_like = "era5_hourly",
    time_tuple = time_range,
    variables = 'precip',
    geom=region, 
    buffer = 2 # cells
)

precip_data_zarr = data_catalog.get_rasterdataset(
    data_like = "era5_hourly_zarr",
    time_tuple = time_range,
    variables = 'precip',
    geom=region, 
    buffer = 2 # cells
)

#%%
# Load raster data for specified region and time range
total_runoff_data = data_catalog.get_rasterdataset(
    data_like = "era5_hourly",
    time_tuple = time_range,
    variables = 'ro',
    geom=region, 
    buffer = 2 # cells
)

runoff_coef = data_catalog.get_rasterdataset(
    data_like = "era5_hourly",
    time_tuple = time_range,
    variables = 'src',
    geom=region, 
    buffer = 2 # cells
)
#%%
# Adjust precipitation values acc to the CF_value expressed as a factor
precip_data.load()
# %%
total_runoff_data.load()
# %%
total_runoff_data.plot()

# %%
lat = -17.75
lon = 32
precip_data.sel(latitude=lat, longitude=lon).plot()

total_runoff_data.sel(latitude=lat, longitude=lon).plot()
# %%

precip_data.to_netcdf("c:/Code/Data_Poppy/era5_hourly_MZB.nc")
# %%

#%%
# Importing the necessary packages
import hydromt
import geopandas as gpd
from datetime import datetime
import matplotlib.pyplot as plt

tc_name = "Idai"
start_date = "20190309 000000"
end_date = "20190325 060000"
wflow_region = "p:/11210471-001-compass/02_Models/sofala/Idai/wflow/staticgeoms/region.geojson"

data_cat = ['../../../03_data_catalogs/datacatalog_general.yml'] 

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
# Original data catalog entry for era5_hourly using lock:false
era5_tp = data_catalog.get_rasterdataset(
    data_like = 'era5_hourly',
    time_tuple = time_range,
    variables = 'precip',
    geom=region, 
)

# era5_hourly_zarr has auto chunking, therefore also lock: true
era5_tp_zarr = data_catalog.get_rasterdataset(
    data_like = 'era5_hourly_zarr',
    time_tuple = time_range,
    variables = 'precip',
    geom=region, 
)

# Using Deltares datacatalog (where lock is true as default)
data_catalog_Delt = hydromt.DataCatalog(data_libs=["deltares_data"])
era5_tp_Delt = data_catalog_Delt.get_rasterdataset(
    data_like = 'era5_hourly',
    time_tuple = time_range,
    variables = 'precip',
    geom=region, 
)

#%%
# Loading era5 precip data using our own datacatalog without lock: false
# reload data catalog after removing lock:false
data_cat = ['../../../03_data_catalogs/datacatalog_general.yml'] 
data_catalog = hydromt.data_catalog.DataCatalog(data_libs = data_cat)

# Load raster data for specified region and time range
era5_tp_lock_TRUE = data_catalog.get_rasterdataset(
    data_like = 'era5_hourly',
    time_tuple = time_range,
    variables = 'precip',
    geom=region, 
)
# %%
lat = -17.75
lon = 32.5

# era5_tp.sel(longitude=lon, latitude=lat).plot()
era5_tp_zarr.sel(longitude=lon, latitude=lat).plot()
era5_tp_Delt.sel(longitude=lon, latitude=lat).plot()
era5_tp_lock_TRUE.sel(longitude=lon, latitude=lat).plot()
# %%
lon, lat = 32.5, -17.75  # Example coordinates

# Create a figure and axis
plt.figure(figsize=(10, 5))

# Plot each dataset with labels
# plt.plot(era5_tp.sel(longitude=lon, latitude=lat), label="ERA5 hourly Lock False")
plt.plot(era5_tp_zarr.sel(longitude=lon, latitude=lat), label="ERA5 hourly Zarr")
plt.plot(era5_tp_Delt.sel(longitude=lon, latitude=lat), label="ERA5 hourly from Deltares datacatalog")
plt.plot(era5_tp_lock_TRUE.sel(longitude=lon, latitude=lat), label="ERA5 hourly Lock True")

# Add title and labels
plt.title("Precipitation Time Series at ({}, {})".format(lon, lat))
plt.xlabel("Time")
plt.ylabel("Precipitation (mm)")
plt.legend()  # Add a legend

# Show the plot
plt.show()

# %%

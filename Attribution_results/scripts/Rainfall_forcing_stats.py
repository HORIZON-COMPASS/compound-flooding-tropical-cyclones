
# ## Compare rainfall forcing
#%%
# Import the correct packages, use compass-wflow pixi environment
import geopandas as gpd
import numpy as np
import hydromt

# In[2]:
# Read data catalog
path_data_cat = "../../Workflows/03_data_catalogs/datacatalog_general.yml"
data_catalog = hydromt.data_catalog.DataCatalog(data_libs = path_data_cat)

sfincs        = "p:/11210471-001-compass/02_Models/sofala/Idai/sfincs/gis/region.geojson"
wflow_region  = "p:/11210471-001-compass/02_Models/sofala/Idai/wflow/staticgeoms/region.geojson"
wflow_basins  = "p:/11210471-001-compass/02_Models/sofala/Idai/wflow/staticgeoms/basins.geojson"

# region of eastern Africa
bbox_big = (29,-27,46,-9)
region = gpd.read_file(wflow_region)
sfincs_region = gpd.read_file(sfincs).to_crs(epsg=4326).total_bounds
basins = gpd.read_file(wflow_basins)

# time range of the event
start_date = np.datetime64('2019-03-09T00:00') 
end_date = np.datetime64('2019-03-25T06:00') 
time_range = (start_date, end_date)

# In[3]:
# Load raster data for specified region and time range for ERA5
era5_tp = data_catalog.get_rasterdataset(
    data_like = 'era5_hourly_zarr',
    time_tuple = time_range,
    variables = 'precip',
    geom=region,
    buffer=2
)

# and for whole of Mozambique (MZ)
era5_tp_MZ = data_catalog.get_rasterdataset(
    data_like = 'era5_hourly_zarr',
    time_tuple = time_range,
    variables = 'precip',
    bbox=bbox_big,
    buffer=2
)

# and for SFINCS domain
era5_tp_sfincs = data_catalog.get_rasterdataset(
    data_like = 'era5_hourly_zarr',
    time_tuple = time_range,
    variables = 'precip',
    bbox=sfincs_region,
    buffer=2
)

# ERA5 rainfall statistics
min = era5_tp.sum(dim='time').min().values
max = era5_tp.sum(dim='time').max().values
print(f"The accumulated rainfall over the wflow region ranged from {min:.1f} and {max:.1f} mm (min and max).")

mean = era5_tp.sum(dim='time').mean().values
print(f"The mean accumulated rainfall over the wflow region is {mean:.1f} mm")

lower = era5_tp.sum(dim='time').quantile(0.05).values
upper = era5_tp.sum(dim='time').quantile(0.95).values
print(f"Rainfall over the wflow region ranged between {lower:.1f} and {upper:.1f} mm (5th–95th percentile).")

mean_sfincs = era5_tp_sfincs.sum(dim='time').mean().values
mean_sfincs_10 = np.round(mean_sfincs, -1)
print(f"The mean accumulated rainfall over the sfincs domain is {mean_sfincs_10:.0f} mm")


# %%

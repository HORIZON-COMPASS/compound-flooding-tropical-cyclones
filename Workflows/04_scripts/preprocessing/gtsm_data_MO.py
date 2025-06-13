#%%
# Importing the necessary packages
import hydromt
import geopandas as gpd
from datetime import datetime
import matplotlib.pyplot as plt
import yaml
import ast
import xarray as xr

#%%
# Load config file
config_file = "../../01_config_snakemake/config_general_somerset.yml"
with open(config_file, "r") as f:
    config = yaml.safe_load(f)
    config = config['runname_ids']['SomersetLevels']

# DFM model region 
bbox_dfm = ast.literal_eval(config['bbox_dfm'])
start_date = config['start_time']
end_date =config['end_time']


# Read data catalog
data_cat = ['../../03_data_catalogs/datacatalog_general.yml'] 
data_catalog = hydromt.data_catalog.DataCatalog(data_libs = data_cat)


# Convert to datetime objects
start_dt = datetime.strptime(start_date, "%Y%m%d %H%M%S")
end_dt = datetime.strptime(end_date, "%Y%m%d %H%M%S")

# Pass as a time tuple for HydroMT
time_range = (start_dt, end_dt)
#%%
# # Load raster data for specified region and time range
# gtsm_somerset = data_catalog.get_geodataset(
#     data_like = 'gtsm_reanalysis_waterlevel_hourly',
#     time_tuple = time_range,
#     # variables = 'precip',
#     bbox=bbox_dfm, 
#     buffer=2000, # meters
# )

# # %%
# # export to netcdf 
# gtsm_somerset.to_netcdf(path="C:/Code/gtsm_somerset.nc")


# %%
# Do the same for the modis_lai data
# bbox = [-4, -2, 50, 52]

# # Load raster data for specified region and time range
# modisLAI_somerset = data_catalog.get_rasterdataset(
#     data_like = 'modis_lai',
#     time_tuple = time_range,
#     bbox=bbox, 
#     buffer=2000, # meters
# )


# %%
# Export to netcdf
modisLAI_somerset.to_netcdf(path="C:/Code/modis_lai_somerset.nc")
# %%
# Do the same for the soilgrid data
bbox = [-4, -2, 50, 52]

# Load raster data for specified region and time range
soilgrids_somerset = data_catalog.get_rasterdataset(
    data_like = 'soilgrids',
    time_tuple = time_range,
    bbox=bbox, 
    buffer=1000, # meters
)


# %%
# Export to netcdf
# soilgrids_somerset.to_netcdf(path="C:/Code/soilgrids_somerset.nc")
# %%
soilgrids_somerset.to_netcdf(
    "p:/11210471-001-compass/01_Data/soilgrids_somerset_MO.nc",
    engine="netcdf4",
    encoding={v: {"zlib": True, "complevel": 4} for v in soilgrids_somerset.data_vars}
)

with zipfile.ZipFile("p:/11210471-001-compass/01_Data/soilgrids_somerset_MO.zip", 'w', compression=zipfile.ZIP_DEFLATED) as zf:
    zf.write("p:/11210471-001-compass/01_Data/soilgrids_somerset_MO.nc", arcname="p:/11210471-001-compass/01_Data/soilgrids_somerset_MO.nc")

# %%

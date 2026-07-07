#%% Works with hydromt-sfincs environment but not compass-snake-sfincs
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

# Read snakemake values if run in workflow, 
# otherwise use absolute values to test script
if "snakemake" in locals():
    tc_name = snakemake.wildcards.runname
    bbox = os.path.abspath(snakemake.input.wflow_bbox)
    start_date = snakemake.params.start_date
    end_date = snakemake.params.end_date
    data_cat = snakemake.params.data_cat
    CF_catalog_path = snakemake.params.CF_data_cat
    precip_name = snakemake.wildcards.precip_forcing
    CF_value = float(snakemake.wildcards.CF_rain)
    CF_value_txt = snakemake.wildcards.CF_rain
    output_CF_rainfall = os.path.abspath(snakemake.output.CF_rainfall)
    CF_data_cat = os.path.abspath(snakemake.params.CF_data_cat) 
else:
    tc_name = "Idai"
    start_date = "20190306 000000" # Ensure at least 2 days prior for wflow warm up
    end_date = "20190325 060000"
    wflow_region = f"p:/11210471-001-compass/02_Models/sofala/{tc_name}/wflow/staticgeoms/region.geojson"
    data_cat = [
        '../../../03_data_catalogs/datacatalog_ClimateDT.yml',
    ] 
    ds_name = "ClimateDT_hist_r1"
    # CF_value = -16
    # CF_value_txt = f"{CF_value}"
    # output_CF_rainfall = f"p:/11210471-001-compass/01_Data/counterfactuals/precipitation/{precip_name}_CF{CF_value}_{tc_name}.nc"
    CF_catalog_path = "../../../03_data_catalogs/datacatalog_CF_forcing.yml"


#%%
test = xr.open_dataset("p:/11210471-001-compass/01_Data/ECMWF_ClimateDT/data/preprocessed/Gen2/climateDT_msl_hist_r1_Sofala_20190301_to_20190329.nc")


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
test_data = data_catalog.get_rasterdataset(
    data_like = ds_name,
    time_tuple = time_range,
    # variables = 'precip',
    geom=region, 
    buffer = 2 # cells
)

# %%

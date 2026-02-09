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
        '../../../03_data_catalogs/datacatalog_general.yml',
        '../../../03_data_catalogs/datacatalog_CF_forcing.yml',
    ] 
    precip_name = "era5_hourly_zarr"
    CF_value = -16
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
    data_like = precip_name,
    time_tuple = time_range,
    variables = 'precip',
    geom=region, 
    buffer = 2 # cells
)

#%%
# Adjust precipitation values acc to the CF_value expressed as a factor
precip_data_CF = precip_data * ((100 + CF_value)/ 100)


#%%
# Save the data in as NetCDF
ds_out = precip_data_CF.to_dataset(name='precip')  # gives a Dataset with one variable
ds_out = ds_out.drop_vars("spatial_ref")

ds_out.to_netcdf(output_CF_rainfall)


#%%
# Adding modified rainfall dataset to the CF data catalog
CF_datacatalog = hydromt.data_catalog.DataCatalog(data_libs=[CF_catalog_path])
CF_dataset_name = f"{precip_name}_CF{CF_value_txt}_{tc_name}"

 # Add the CF data to the CF data catalog
CF_datacatalog[CF_dataset_name] = copy.deepcopy(data_catalog[f"{precip_name}"])
CF_datacatalog[CF_dataset_name].path = os.path.abspath(output_CF_rainfall)
CF_datacatalog[CF_dataset_name].driver = "netcdf"

# Check if the 'rename' attribute exists and delete it
if all(hasattr(CF_datacatalog[CF_dataset_name], attr) for attr in ['rename', 'unit_mult', 'unit_add']):
    del CF_datacatalog[CF_dataset_name].rename
    del CF_datacatalog[CF_dataset_name].unit_mult
    del CF_datacatalog[CF_dataset_name].unit_add

# Check if 'meta' exists and is a dictionary
if hasattr(CF_datacatalog[CF_dataset_name], 'meta') and isinstance(CF_datacatalog[CF_dataset_name].meta, dict):
    # Remove all metadata except for 'notes' and add CF info
    notes = CF_datacatalog[CF_dataset_name].meta.pop('notes', None)  # Save the current 'notes' if exists
    CF_datacatalog[CF_dataset_name].meta.clear()

    if notes is not None:
        CF_datacatalog[CF_dataset_name].meta['notes'] = f"Copied {precip_name} dataset and adjusted with CF value {CF_value_txt}%"
    
# Save the updated catalog
CF_datacatalog.to_yml(CF_catalog_path)


#%%

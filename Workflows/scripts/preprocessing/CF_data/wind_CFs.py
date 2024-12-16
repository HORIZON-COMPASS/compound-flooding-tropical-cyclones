# Importing the necessary packages
import os
import numpy as np
import xarray as xr
# import hydromt
import geopandas as gpd
import copy
import yaml

# Read snakemake values if run in workflow, 
# otherwise use absolute values to test script
if "snakemake" in locals():
    start_date = np.datetime64(snakemake.params.start_date) 
    end_date = np.datetime64(snakemake.params.end_date) 
    bbox = ast.literal_eval(snakemake.params.bbox)
    path_data_cat = snakemake.params.data_cat
    CF_catalog_path = snakemake.params.CF_data_cat
    wind_name = snakemake.wildcards.wind_name
    CF_value = float(snakemake.wildcards.CF_value_wind)
    CF_value_txt = snakemake.wildcards.CF_value_wind
    output_CF_wind = os.path.abspath(snakemake.output.CF_wind)
else:
    start_date = np.datetime64("2019-03-09") 
    end_date = np.datetime64("2019-03-24") 
    bbox = [34.33,-20.12,34.95,-19.30]
    path_data_cat = os.path.abspath("../../../data_catalogs/datacatalog_general.yml")
    wind_name = "ERA5land_Idai"
    CF_value = -10
    CF_value_txt = "-10"
    output_CF_rainfall = os.path.abspath(f"p:/11210471-001-compass/01_Data/counterfactuals/wind/{wind_name}_{CF_value}.nc")
    CF_catalog_path = os.path.abspath("../../../data_catalogs/datacatalog_CF_forcing.yml")
    
# Read data catalog
data_catalog = hydromt.DataCatalog(data_libs = [path_data_cat])

# Load raster data
precip_data = data_catalog.get_rasterdataset(
    data_like = precip_name,
    variables = 'precip',
    geom=region, 
)

# Adjust precipitation values acc to the CF_value expressed as a factor
precip_data_CF = precip_data.copy()
precip_data_CF = precip_data_CF * ((100 + CF_value)/100)

# Saving CF precip dataset
precip_data_CF = precip_data_CF.to_dataset(name='precip')
precip_data_CF.to_netcdf(output_CF_rainfall)
precip_data_CF.close()


# Adding modified rainfall dataset to the CF data catalog
CF_datacatalog = hydromt.DataCatalog(data_libs=[CF_catalog_path])
CF_dataset_name = f"{precip_name}_{CF_value_txt}"

# # Add the CF data to the CF data catalog
if CF_dataset_name not in CF_datacatalog:
    CF_datacatalog[CF_dataset_name] = copy.deepcopy(data_catalog[f"{precip_name}"])
    CF_datacatalog[CF_dataset_name].path = os.path.abspath(output_CF_rainfall)

    # Check if the 'rename' attribute exists and delete it
    if hasattr(CF_datacatalog[CF_dataset_name], 'rename'):
        del CF_datacatalog[CF_dataset_name].rename

    # Check if 'meta' exists and is a dictionary
    if hasattr(CF_datacatalog[CF_dataset_name], 'meta') and isinstance(CF_datacatalog[CF_dataset_name].meta, dict):
        # Remove all metadata except for 'notes' and add CF info
        notes = CF_datacatalog[CF_dataset_name].meta.pop('notes', None)  # Save the current 'notes' if exists
        CF_datacatalog[CF_dataset_name].meta.clear()

        if notes is not None:
            CF_datacatalog[CF_dataset_name].meta['notes'] = f"Copied {precip_name} dataset and adjusted with CF value {CF_value_txt}%"
    
    # Save the updated catalog
    CF_datacatalog.to_yml(CF_catalog_path)

else:
    print(f"Dataset '{CF_dataset_name}' is already in the catalog.")

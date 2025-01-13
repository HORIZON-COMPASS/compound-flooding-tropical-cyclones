# The dfm output is added to a dfm specific data catalog
# Importing the necessary packages
import os
import numpy as np
import ast
import xarray as xr
import hydromt
import geopandas as gpd
import copy
import yaml

#%%
if "snakemake" in locals():
    path_data_cat = os.path.abspath(snakemake.params.dfm_data_cat)
    his_path = 
    run_dir = 
else:
    path_data_cat = os.path.abspath("../../../data_catalogs/datacatalog_dfm_output.yml")

#%%
# Adding modified rainfall dataset to the CF data catalog
dfm_datacatalog = hydromt.DataCatalog(data_libs=[path_data_cat])
dfm_run = f"{precip_name}_{CF_value_txt}"

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

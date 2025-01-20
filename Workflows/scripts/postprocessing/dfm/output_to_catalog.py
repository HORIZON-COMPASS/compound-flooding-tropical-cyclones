#%% The dfm output is added to a dfm specific data catalog
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
    his_path = os.path.abspath(snakemake.input.his_file)
    path_data_cat = os.path.abspath(snakemake.output.sfincs_data_cat)
    run_dir = os.path.abspath(snakemake.params.dir_event_model)
    model_name = snakemake.params.model_name
    tc_name = snakemake.wildcards.tc_name
else:
    region = "sofala"
    tc_name = "Idai"
    dfm_res = "450"
    bathy = "gebco2024"
    tidemodel = 'GTSMv41opendap' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap
    wind_forcing = "spw_IBTrACS_ext"
    model_name = f'event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}_TEST2'
    path_data_cat = os.path.abspath("../../../data_catalogs/datacatalog_SFINCS_coastal_coupling.yml")
    run_dir = f'p:/11210471-001-compass/03_Runs/{region}/{tc_name}/dfm/{model_name}'
    his_path = f"{model_name}_his.nc"

#%%
# Adding modified rainfall dataset to the CF data catalog
datacatalog = hydromt.DataCatalog(data_libs=[path_data_cat])
dfm_run = f"dfm_output_{model_name}"

# # Add the CF data to the CF data catalog
if dfm_run not in datacatalog:
    datacatalog[dfm_run] = copy.deepcopy(datacatalog[f"dfm_output_MZB_{tc_name}"])
    datacatalog[dfm_run].path = os.path.abspath(his_path)
    
    # Save the updated catalog
    datacatalog.to_yml(path_data_cat)

else:
    print(f"Dataset '{dfm_run}' is already in the catalog.")

# %%

# %%

#%%
# Importing the necessary packages
import os
import numpy as np
import xarray as xr
import hydromt
import geopandas as gpd
import copy
import yaml
import dfm_tools as dfmt

#%%
# if "snakemake" in locals():
#     start_date = np.datetime64(snakemake.params.start_date) 
#     end_date = np.datetime64(snakemake.params.end_date) 
#     domain = os.path.abspath(snakemake.input.model_domain)
#     path_data_cat = os.path.abspath(snakemake.input.data_cat)
#     CF_catalog_path = os.path.abspath(snakemake.input.CF_data_cat)
#     precip_name = snakemake.wildcards.precip_name
#     CF_value = float(snakemake.wildcards.CF_value_rain)
#     CF_value_txt = snakemake.wildcards.CF_value_rain
#     output_CF_rainfall = os.path.abspath(snakemake.output.CF_rainfall)
# else:
#     start_date = np.datetime64("2019-03-09") 
#     end_date = np.datetime64("2019-03-24") 
#     domain = os.path.abspath("../../runs/SFINCS/base_model/gis/region.geojson")
#     path_data_cat = os.path.abspath("../../input_data/factual/datacatalog_general.yml")
#     precip_name = "ERA5land_Idai"
#     CF_value = 10
#     output_CF_rainfall = os.path.abspath(f"../../input_data/counterfactual/{precip_name}_{CF_value}.nc")
#     CF_catalog_path = os.path.abspath("../../input_data/counterfactual/datacatalog_CF.yml")

# interpolate tidal components to boundary conditions file (.bc)

file_path_bc = (r"p:\11210471-001-compass\02_Models\Delft3DFM\mozambique_model\computations\mozambique_spw_Idai_areaBeira_500m_gebco2023_doris\tide_FES2014.bc")
# ds = xr.open_dataset(file_path)
# %%
from hydrolib.core.dflowfm.bc.models import ForcingModel

# Path to your .bc file
# file_bc_path = r"p:\...\tide_FES2014.bc"

# Load the .bc file
forcing_model = ForcingModel(filepath=file_path)

# Inspect the contents
print(forcing_model)
# %%
import hydrolib.core.dflowfm as hcdfm
file_pli = r"p:\11210471-001-compass\02_Models\Delft3DFM\mozambique_model\computations\mozambique_spw_Idai_areaBeira_500m_gebco2023_doris\mozambique_spw_Idai_areaBeira_500m_gebco2023_doris.pli"
boundary_object = hcdfm.Boundary(quantity='waterlevelbnd', 
                                 locationfile=file_pli,
                                 forcingfile=forcing_model)
# %%
# from dfm_tools import Boundary

# Path to your boundary condition file (.bc)
bc_file_path = file_path_bc
# Load the boundary condition file
bc = boundary_object

# Define the sea level rise (e.g., 0.5 meters)
sea_level_rise = 0.5  # meters

# Modify the water level boundary condition (assuming it's named 'waterlevelbnd')
# if 'waterlevelbnd' in bc.boundary_names:
    # Access the boundary and add the sea level rise
bc.quantity['waterlevelbnd'].data = bc.quantity['waterlevelbnd'].data + sea_level_rise


# Save the updated boundary condition file
# bc.save(r"path_to_updated_boundary_condition_file.bc")


# %%




# %%

#%%
# Activate the compass-fiat-hydromt environment & import required packages 
import os
from hydromt_fiat.fiat import FiatModel
from hydromt.log import setuplog
from pathlib import Path
import geopandas as gpd
import matplotlib.pyplot as plt
import shutil
import pandas as pd

#%%
# Define output folder and model name
# model_root = "computations"  # where to save the FIAT model
model_name = "sfincs_idai_test" # name of the case study

#%%
# Set up data catalog and hydromt logger
model_folder = (Path(os.path.join("p:/11210471-001-compass/03_Runs/FIAT/", model_name)))  # path to model folder
data_catalog = (Path(os.path.abspath("")).resolve().parent / "hydromt_fiat_catalog_global.yml") # path to data catalog
logger_name = "hydromt_fiat"  # name of the logger
logger = setuplog(logger_name, log_level=10) # setup logger

#%%
# Read in the region of interest (e.g. from SFINCS model output)
region = gpd.read_file(Path("p:/11210471-001-compass/03_Runs/sfincs_Idai/gis/region.geojson"))
continent = "Africa"
country = "Mozambique"

# Uncomment to create a region from lat lon coordinates
# area_of_interest = {
#     "type": "FeatureCollection",
#     "features": [
#         {
#             "type": "Feature",
#             "properties": {},
#             "geometry": {
#                 "coordinates": [
#                      [
#             [34.37939695948343, -20.07031561711855],
#             [34.387653838143834, -19.377711500714003],
#             [34.92179500387799, -19.368799942263994],
#             [34.92293033010172, -20.07136792215392],
#             [34.37939695948343, -20.07031561711855]
#           ]
#         ],
#                 "type": "Polygon",
#             },
#         }
#     ],
# }
# region = gpd.GeoDataFrame.from_features(area_of_interest, crs="EPSG:4326")

#%%
### Setup vulnerability parameters ###
vulnerability_fn = "jrc_vulnerability_curves"
vulnerability_identifiers_and_linking_fn = "jrc_vulnerability_curves_linking"
unit = "meters"

### Setup exposure parameters ###
asset_locations = "OSM"
occupancy_type = "OSM"
max_potential_damage = "jrc_damage_values"
ground_floor_height = 0
damage_types = ["total"]
unit = "meters"
crs = "EPSG:4326"

# Set up hazard parameters
path_map = (r"p:\11210471-001-compass\03_Runs\sfincs_Idai\gis\floodmap.tif") # change this to your own floodmap
crs_flood = "EPSG:32736"    # original crs of floodmap
map_type = "max_depth"
# not sure about vertical ref but the default "datum" does not give different result than "DEM"

### Setup output parameters ###
output_dir = "output"
output_csv_name = "output.csv"
output_vector_name = "spatial.gpkg"

#%%
# Define model configurations
configuration = {
    "setup_output": {
        "output_dir": output_dir,
        "output_csv_name": output_csv_name,
        "output_vector_name": output_vector_name,
    },
    "setup_vulnerability": {
        "vulnerability_fn": vulnerability_fn,
        "vulnerability_identifiers_and_linking_fn": vulnerability_identifiers_and_linking_fn,
        "continent": continent,
        "unit": unit,
    },
    "setup_exposure_buildings": {
        "asset_locations": asset_locations,
        "occupancy_type": occupancy_type,
        "max_potential_damage": max_potential_damage,
        "ground_floor_height": ground_floor_height,
        "unit": unit,
        "dst_crs": crs,
        "damage_types": damage_types,
        "country": country,
    },
    "setup_hazard": {
        "map_fn": path_map,
        "map_type": map_type,
        "crs": crs_flood, 
    },
    "setup_global_settings": {
        "crs": 4326,
    },
}

#%%
# Set up model
if model_folder.exists():
    shutil.rmtree(model_folder)
fiat_model = FiatModel(root=model_folder, mode="w", data_libs=[data_catalog], logger=logger)

#%%
# Build the model
fiat_model.build(region={"geom": region}, opt=configuration, write=False)

#%%
# Get the geodataframe with exposure data to check
# gdf = fiat_model.exposure.get_full_gdf(fiat_model.exposure.exposure_db)
# gdf.to_file('exposure_test.geojson', driver='GeoJSON')

#%%
# Get the range of (possible) water depths
water_depths = fiat_model.vulnerability.hazard_values
# Plot damage curves for some occupancy types
line_styles = ['--', '-', ':']
for function_name, ls in zip(fiat_model.vulnerability.functions.keys(), line_styles):
    dmg = [float(i) for i in fiat_model.vulnerability.functions[function_name]]
    plt.plot(water_depths, dmg, label=function_name, ls=ls)
plt.xlabel('depth (m)')
plt.ylabel('damage fraction (-)')
plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
plt.show()

#%%
# Write the model
fiat_model.write()

#%%
fiat_model_new = FiatModel(root=model_folder, mode="r", data_libs=[data_catalog], logger=logger)
fiat_model_new.read()

#%% Some bug fixes to make the model run
# First for the buildings file
# Load the GeoPackage file
gdf = gpd.read_file(f'{model_folder}/exposure/buildings.gpkg')

# Rename the 'Object ID' column to 'object_id'
gdf.rename(columns={'Object ID': 'object_id'}, inplace=True)

# Save the changes back to the file
gdf.to_file(f'{model_folder}/exposure/buildings.gpkg', driver="GPKG")


# And secondly for the exposure csv file
# Load the CSV file
df = pd.read_csv(f"{model_folder}/exposure/exposure.csv")

# Rename columns
df = df.rename(columns={
    "Object ID": "object_id",
    "Primary Object Type": "primary_object_type",
    "Secondary Object Type": "secondary_object_type",
    "Max Potential Damage: Total": "max_damage_total",
    "Ground Floor Height": "ground_flht",
    "Extraction Method": "extract_method",
    "Ground Elevation": "ground_elevtn",
    "Damage Function: Total": "fn_damage_total"
})

# Save the updated DataFrame to a new CSV
df.to_csv(f"{model_folder}/exposure/exposure.csv", index=False)
#%%
# To run the model, use the "execute_fiat_example.ipynb" script

#%%
# Activate the compass-fiat-hydromt environment & import required packages 
import os
from hydromt_fiat.fiat import FiatModel
from hydromt.log import setuplog
from pathlib import Path
import geopandas as gpd
import shutil
from hydromt.config import configread
import rasterio
import toml

#%%
if "snakemake" in locals():
    sfincs_mod_with_forcing = Path(os.path.abspath(snakemake.params.dir_run_with_forcing))
    data_catalog            = snakemake.params.datacat_fiat
    model_folder            = Path(os.path.abspath(snakemake.params.model_folder))
    continent               = snakemake.params.continent
    country                 = snakemake.params.country
    config_file             = snakemake.params.config
    floodmap                = snakemake.input.floodmap
else:
    continent               = "Africa"
    country                 = "Mozambique"
    region                  = "sofala"
    tc_name                 = "Idai"
    wind_forcing            = 'spw_IBTrACS'
    precip_forcing          = 'era5_hourly_zarr'
    bathy                   = "gebco2024_MZB"
    tidemodel               = 'GTSMv41opendap' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap, tpxo80_opendap
    data_catalog            = ['../../../03_data_catalogs/datacatalog_fiat.yml']    
    CF_rain_txt             = "0"
    CF_SLR_txt              = "0"
    CF_wind_txt             = "0"
    model_name              = f"event_tp_{precip_forcing}_CF{CF_rain_txt}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}"
    sfincs_mod_with_forcing = os.path.join(f"p:/11210471-001-compass/03_Runs/{region}/{tc_name}/sfincs/{model_name}")
    model_folder            = (Path(os.path.join("p:/11210471-001-compass/03_Runs/sofala/Idai/fiat", model_name)))  # path to model folder
    # model_folder            = Path(os.path.join("c:/Code/Delft-FIAT/tests", model_name))
    config_file             = '../../../05_config_models/03_fiat/config_fiat.yml'
    floodmap                = f"p:/11210471-001-compass/03_Runs/{region}/{tc_name}/sfincs/{model_name}/plot_output/floodmap.tif"

# %%
# Set up data catalog and hydromt logger
logger_name = "hydromt_fiat"  # name of the logger
logger = setuplog(logger_name, log_level=10) # setup logger
config = configread(config_file, abs_path=True)  # read settings from ini file

#%%
# Read in the region of interest (e.g. from SFINCS model output)
region = gpd.read_file(os.path.join(sfincs_mod_with_forcing, "gis/region.geojson"))

#%%
# Get the crs of the flood map
with rasterio.open(floodmap) as src:
    crs_flood = src.crs  # Get the CRS
    crs_flood = crs_flood.to_string()

#%%
# Set up region specific parameters
config["setup_vulnerability"]["continent"]    = continent
config["setup_exposure_buildings"]["country"] = country
config["setup_hazard"]["map_fn"]              = floodmap
config["setup_hazard"]["crs"]                 = crs_flood

#%%
# Set up model
# if model_folder.exists():
#     shutil.rmtree(model_folder)
fiat_model = FiatModel(root=model_folder, mode="w", data_libs=data_catalog, logger=logger)

#%%
# Build and write the model
fiat_model.build(region={"geom": region}, opt=config, write=True)


#%%
# Debugging to allow running the model from a different location than python environment is stored
# Load the buildings.gpkg file
gdf = gpd.read_file(f"{model_folder}/exposure/buildings.gpkg")

# Save as .fgb
gdf.to_file(f"{model_folder}/exposure/buildings.fgb", driver="FlatGeobuf")

#%%
# Refer to the new file in the settings.toml
with open(f"{model_folder}/settings.toml", "r") as f:
    settings = toml.load(f)

# Update the file path
settings["exposure"]["geom"]["file1"] = "exposure/buildings.fgb"
settings["output"]["geom"]["name1"] = "spatial.fgb"

# Save the updated TOML file
with open(f"{model_folder}/settings.toml", "w") as f:
    toml.dump(settings, f)
# %%
# To run the model, use the "execute_fiat_example.ipynb" script

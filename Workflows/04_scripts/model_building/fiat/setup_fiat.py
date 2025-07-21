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
import tempfile
import pandas as pd
import logging
from pathlib import Path
import platform
prefix = "p:/" if platform.system() == "Windows" else "/p/"

#%%
if "snakemake" in locals():
    sfincs_mod_with_forcing = Path(os.path.abspath(snakemake.params.dir_run_with_forcing))
    data_catalog            = Path(os.path.abspath(snakemake.params.datacat_fiat))
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
    wind_forcing            = 'era5_hourly_spw_IBTrACS'
    precip_forcing          = 'era5_hourly_zarr'
    bathy                   = "gebco2024_MZB"
    tidemodel               = 'GTSMv41' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap, tpxo80_opendap
    data_catalog            = '../../../03_data_catalogs/datacatalog_fiat.yml'  
    CF_rain_txt             = "-8"
    CF_SLR_txt              = "0"
    CF_wind_txt             = "0"
    model_name              = f"event_tp_{precip_forcing}_CF{CF_rain_txt}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}"
    sfincs_mod_with_forcing = os.path.join(f"p:/11210471-001-compass/03_Runs/{region}/{tc_name}/sfincs/{model_name}")
    model_folder            = (Path(os.path.join("p:/11210471-001-compass/03_Runs/sofala/Idai/fiat", model_name)))  # path to model folder
    # model_folder            = Path(os.path.join("c:/Code/Delft-FIAT/tests", model_name))
    config_file             = '../../../05_config_models/03_fiat/fiat_base_build.yml'
    floodmap                = f"p:/11210471-001-compass/03_Runs/{region}/{tc_name}/sfincs/{model_name}/plot_output/floodmap.tif"
    sofala_exposure         = 'p:/11210471-001-compass/01_Data/fiat/sofala/exposure'

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
# Define temporary directory for building the model
with tempfile.TemporaryDirectory() as tmpdir:
    print("temporary directory for building the model")
    tmp_model_base = Path(tmpdir)
    tmp_model_folder = tmp_model_base / "fiat_model_temp"

    # Set up the model in temp directory
    print("Build the model")
    fiat_model = FiatModel(root=tmp_model_folder, mode="w+", data_libs=[data_catalog], logger=logger)
    fiat_model.build(region={"geom": region}, opt=config, write=True)

    # Try .fgb first, fallback to .gpkg
    fgb_path = tmp_model_folder / "exposure" / "buildings.fgb"
    gpkg_path = tmp_model_folder / "exposure" / "buildings.gpkg"

    if fgb_path.exists():
        pass
    elif gpkg_path.exists():
        gdf = gpd.read_file(gpkg_path)
        gdf.to_file(f"{tmp_model_folder}/exposure/buildings.fgb", driver="FlatGeobuf")
    else:
        raise FileNotFoundError("Neither buildings.fgb nor buildings.gpkg found in the exposure folder.")

    # Update settings.toml
    toml_file = tmp_model_folder / "settings.toml"
    with open(toml_file, "r") as f:
        settings = toml.load(f)
    settings["exposure"]["geom"]["file1"] = "exposure/buildings.fgb"
    settings["output"]["geom"]["name1"] = "spatial.fgb"
    with open(toml_file, "w") as f:
        toml.dump(settings, f)

    # Ensure Linux readability for CSV
    vuln_csv = tmp_model_folder / "vulnerability/vulnerability_curves.csv"
    with open(vuln_csv, 'r') as f:
        lines = f.read().splitlines()
    with open(vuln_csv, 'w') as f:
        for line in lines:
            f.write(line + '\n')
    
    # Remove unrealistic building for sofala
    if "sofala" in str(model_folder).lower():
        # Load buildings and exposure
        gdf_buildings = gpd.read_file(tmp_model_folder / "exposure" / "buildings.fgb")
        df_exposure = pd.read_csv(tmp_model_folder / "exposure" / "exposure.csv")

        # ---- Identify the "unrealistic" building ----
        gdf_utm = gdf.to_crs(epsg=32736)

        # Add area column in square meters
        gdf_utm["area_m2"] = gdf_utm.geometry.area
        big_buildings = gdf_utm[gdf_utm['area_m2'] > 1000000]

        if not big_buildings.empty:
            print("Deleting these buildings:")
            print(big_buildings)

            exclude_ids = big_buildings["object_id"].tolist()
            gdf_buildings = gdf_buildings[~gdf_buildings["object_id"].isin(exclude_ids)].copy()

            # Drop corresponding entries from df_exposure
            df_exposure = df_exposure[~df_exposure['object_id'].isin(big_buildings['object_id'])]

            # Save the cleaned files
            gdf_buildings.to_file(os.path.join(tmp_model_folder,"exposure","buildings.fgb"), driver="FlatGeobuf")
            df_exposure.to_csv(os.path.join(tmp_model_folder,"exposure","exposure.csv"), index=False)
        else:
            print("No buildings above the area threshold were found.")

    # Close all handlers attached to the logger
    for handler in logging.root.handlers[:]:
        handler.close()
        logging.root.removeHandler(handler)
        
    # Move CONTENTS of tmp_model_folder into model_folder
    print("Move CONTENTS of tmp_model_folder into model_folder")
    if model_folder.exists():
        shutil.rmtree(model_folder)
    # Create target directory if not exists
    model_folder.parent.mkdir(parents=True, exist_ok=True)

    for item in tmp_model_folder.iterdir():
        shutil.move(str(item), str(model_folder / item.name))

#%%
# To run the model, use the "execute_fiat_example.ipynb" script

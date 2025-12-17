# %% Use pixi environment compass-wflow
from os.path import join, exists
import os
from hydromt.config import configread
from hydromt.log import setuplog
from hydromt_wflow import WflowModel
import geopandas as gpd

# %%
if "snakemake" in locals():
    model_dir        = snakemake.params.dir_model
    config_file      = snakemake.input.config_file
    data_cat         = snakemake.params.data_cat
    region_geom      = snakemake.input.region_geom
    dir_sfincs_model = snakemake.input.dir_sfincs_model
else:
    landuse           = "lisboa_2000"
    model_dir        = f"/p/11210471-001-compass/02_Models/sofala/Idai/wflow_{landuse}"
    config_file      = f"../../../05_config_models/01_wflow/config_wflow_{landuse}.yml"
    data_cat         = [
        '../../../03_data_catalogs/datacatalog_general___linux.yml',
        '../../../03_data_catalogs/datacatalog_CF_forcing___linux.yml'
        ] 
    region_geom      = f'/p/11210471-001-compass/02_Models/sofala/Idai/sfincs_{landuse}/gis/region.geojson'
    dir_sfincs_model = f'/p/11210471-001-compass/02_Models/sofala/Idai/sfincs_{landuse}'

#%%
# Check whether model folder exists
if not exists(model_dir):
    os.mkdir(model_dir)

# model and data paths/
logger = setuplog("update", join(model_dir, "hydromt.log"), log_level=10)

# read settings from ini file
opt = configread(config_file, abs_path=True)  
kwargs = opt.pop("global", {})

# Set up output points based on SFINCS inflow river points
opt['setup_gauges'] = {
    'gauges_fn': join(dir_sfincs_model, "gis", "src.geojson"),
    'snap_to_river': True,
    'snap_uparea': True,
    'rel_error': 0.2,
    'derive_subcatch': False,
    'index_col': 'index',
    'basename': 'locs'
}

#%%
# Read SFINCS region
region = gpd.read_file(region_geom).to_crs(epsg = '4326')

#%%
# Set up wflow 
mod = WflowModel(
    root=model_dir, data_libs=data_cat, mode="w+", logger=logger, **kwargs
)
mod.config['csv'] = None
# %% BUILD MODEL
mod.build(region={"basin": region}, opt=opt)
mod.config['csv'] = None
mod.write_config()

# %%

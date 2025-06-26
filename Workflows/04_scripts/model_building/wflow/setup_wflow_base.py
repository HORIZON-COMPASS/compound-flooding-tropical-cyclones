# %%
from os.path import basename, join, exists

import os
from hydromt.config import configread
from hydromt.log import setuplog
from hydromt_wflow import WflowModel
import geopandas as gpd
# %%
if "snakemake" in locals():
    model_dir = snakemake.params.dir_model
    config_file = snakemake.input.config_file
    data_cat = snakemake.params.data_cat
    bbox = snakemake.params.arg_bbox
    region_geom = snakemake.input.region_geom
    dir_sfincs_model = snakemake.input.dir_sfincs_model
    river_upa = snakemake.params.river_upa
else:
    model_dir = r"p:/11210471-001-compass\02_Models\quelimane\Freddy2\wflow"
    config_file = r"c:\Git_repos\COMPASS\Workflows\config_wflow\wflow_build_quelimane.yml"
    data_cat = r"c:\Git_repos\COMPASS\Workflows\data_catalogs/datacatalog_general.yml"
    bbox = [34.33,-20.12,34.95,-19.30]
    region_geom = r'p:\11210471-001-compass\02_Models\quelimane\Freddy2\sfincs\gis\region.geojson'
    dir_sfincs_model = r'p:\11210471-001-compass\02_Models\quelimane\Freddy2\sfincs'

if not exists(model_dir):
    os.mkdir(model_dir)
# model and data paths/
logger = setuplog("update", join(model_dir, "hydromt.log"), log_level=10)

print("Reading config file:", model_dir)
print("dir_sfincs_model:", dir_sfincs_model)
print("region_geom:", region_geom)
print("data_cat:", data_cat)


opt = configread(config_file, abs_path=True)  # read settings from yml file
kwargs = opt.pop("global", {})

opt['setup_gauges'] = {
    'gauges_fn': join(dir_sfincs_model, "gis", "src.geojson"),
    'snap_to_river': True,
    'snap_uparea': True,
    'rel_error': 0.2,
    'derive_subcatch': False,
    'index_col': 'index',
    'basename': 'locs'
}

opt['setup_rivers']['river_upa'] = river_upa

#%%

region = gpd.read_file(region_geom).to_crs(epsg = '4326')


#%%
mod = WflowModel(
    root=model_dir, data_libs=[data_cat], mode="w+", logger=logger, **kwargs
)
mod.config['csv'] = None
# %% BUILD MODEL
mod.build(region={"basin": region}, opt=opt)
mod.config['csv'] = None
mod.write_config()

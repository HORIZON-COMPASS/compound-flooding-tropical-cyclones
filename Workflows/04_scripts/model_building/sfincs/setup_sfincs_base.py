# %% Use pixi environment compass-snake-sfincs
from os.path import join, exists
import geopandas as gpd
from hydromt.config import configread
import ast
from hydromt.log import setuplog
import hydromt
from hydromt_sfincs import SfincsModel
#%%
def get_local_vector_data(file, bbox, data_cat):
    dataCat = hydromt.data_catalog.DataCatalog(data_cat)
    vector = dataCat.get_geodataframe(
        data_like = file,
        bbox = bbox)
    return vector
# %%
if "snakemake" in locals():
    model_dir        = snakemake.params.dir_model_sfincs
    config_file      = snakemake.input.config_file
    data_cats        = snakemake.params.data_cats
    bbox             = ast.literal_eval(snakemake.params.arg_bbox)
    bathy            = snakemake.params.bathy
    dfm_coastal_mask = snakemake.params.dfm_coastal_mask
    region_name      = snakemake.wildcards.region
else:
    model_dir        = 'p:/11210471-001-compass/02_Models/sofala/Idai/sfincs'
    config_file      = '../../../05_config_models/02_sfincs/sfincs_base_build.yml'
    data_cats        = ['../../../03_data_catalogs/datacatalog_general.yml',
                        '../../../03_data_catalogs/datacatalog_SFINCS_obspoints.yml',
                        '../../../03_data_catalogs/datacatalog_SFINCS_coastal_coupling.yml']
    bbox             = ast.literal_eval("[34.33,-20.12,34.95,-19.30]")
    bathy            = 'gebco2024_MZB'
    dfm_coastal_mask = 'coastal_coupling_msk_MZB'
    region_name      = 'sofala'
    

# Check whether model folder exists. If not, make one
if not exists(model_dir):
    os.mkdir(model_dir)
#%%
# model and data paths/
logger = setuplog("update", join(model_dir, "hydromt.log"), log_level=10)
opt = configread(config_file, abs_path=True)  # read settings from ini file
kwargs = opt.pop("global", {})

# fill in the configuration for SFINCS with arguments from the snakemake config file
opt['setup_dep']['datasets_dep'] = opt['setup_dep']['datasets_dep'] + [{'elevtn': bathy, 'reproj_method': 'bilinear'}]   
opt['setup_subgrid']['datasets_dep'] = opt['setup_subgrid']['datasets_dep'] + [{'elevtn': bathy, 'reproj_method': 'bilinear'}]   

#%%
region = get_local_vector_data(
    file = 'basin_atlas_level12_v10',
    bbox = bbox,
    data_cat = data_cats[0],
)

#%%
# Set up model region
opt['setup_mask_active']['mask'] = region
# opt['setup_mask_active']['mask_buffer'] = 1000
opt['setup_mask_active']['exclude_mask'] = dfm_coastal_mask

if region_name == 'sofala':
    opt['setup_mask_active2']['include_mask'] = 'sofala_incl_polygon'
    opt['setup_mask_active2']['exclude_mask'] = 'sofala_excl_northern_basin'
    opt['setup_mask_active3']['include_mask'] = 'sofala_incl_polygon_2'
else:
    pass

opt['setup_mask_bounds']['include_mask'] = dfm_coastal_mask

#%%
# Initialise model object
mod = SfincsModel(
    root=model_dir, data_libs=data_cats, mode="w+", logger=logger, **kwargs
)

# %% BUILD MODEL
mod.build(region={"geom": region}, opt=opt)


#%% Plot the region and boundaries
fig, ax = mod.plot_basemap(
    fn_out=model_dir, 
    variable="dep",
    plot_bounds=True,
    plot_geoms=True, 
    plot_region=True,
    bmap="sat",
    zoomlevel=12,
    figsize=(8, 6))

# %%

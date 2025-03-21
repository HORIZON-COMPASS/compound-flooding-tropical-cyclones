# %%
from os.path import join, exists

from hydromt.config import configread
import ast
from hydromt.log import setuplog
import hydromt
from hydromt_sfincs import SfincsModel

def get_local_vector_data(file, bbox, data_cat):
    dataCat = hydromt.data_catalog.DataCatalog(data_cat)
    vector = dataCat.get_geodataframe(
        data_like = file,
        bbox=  bbox)
    return vector
# %%
if "snakemake" in locals():
    model_dir = snakemake.params.dir_model_sfincs
    config_file = snakemake.input.config_file
    data_cats = snakemake.params.data_cats
    bbox = ast.literal_eval(snakemake.params.arg_bbox)
    bathy = snakemake.params.bathy
    dfm_coastal_mask = snakemake.params.dfm_coastal_mask
else:
    model_dir = 'p:/11210471-001-compass/02_Models/sofala/Idai/sfincs_test'
    config_file = '../../../05_config_models/02_sfincs/sfincs_base_build.yml'
    data_cats = [
        '../../../03_data_catalogs/datacatalog_general.yml',
        '../../../03_data_catalogs/datacatalog_SFINCS_obspoints.yml',
        '../../../03_data_catalogs/datacatalog_SFINCS_coastal_coupling.yml',
    ]
    bbox = ast.literal_eval("[34.33,-20.12,34.95,-19.30]")
    bathy = 'gebco2024_MZB'
    dfm_coastal_mask = 'coastal_coupling_msk_MZB'

if not exists(model_dir):
    os.mkdir(model_dir)

# model and data paths/
logger = setuplog("update", join(model_dir, "hydromt.log"), log_level=10)
opt = configread(config_file, abs_path=True)  # read settings from ini file
kwargs = opt.pop("global", {})

# fill in the configuration for SFINCS with arguments from the snakemake config file
opt['setup_dep']['datasets_dep'] = opt['setup_dep']['datasets_dep'] + [{'elevtn': bathy, 'reproj_method': 'bilinear'}]   
opt['setup_mask_active']['exclude_mask'] = dfm_coastal_mask
opt['setup_subgrid']['datasets_dep'] = opt['setup_subgrid']['datasets_dep'] + [{'elevtn': bathy, 'reproj_method': 'bilinear'}]   

#%%
region = get_local_vector_data(
    file = 'basin_atlas_level12_v10', #TODO change to 'merit_hydro_basins....' 
    bbox = bbox,
    data_cat = data_cats[0],
)
#%%
# Set up model region
opt['setup_mask_active']['mask'] = region

#%%
# Initialise model object
mod = SfincsModel(
    root=model_dir, data_libs=data_cats, mode="w+", logger=logger, **kwargs
)

# %% BUILD MODEL
mod.build(region={"geom": region}, opt=opt)


# %%
# Include extra polygon
opt2 = {
    'setup_mask_active': {
        'include_mask': 'sofala_incl_polygon',
        'reset_mask': False
        }
}

mod.update(
    write=True,
    # forceful_overwrite=True,
    opt=opt2
)
# %%


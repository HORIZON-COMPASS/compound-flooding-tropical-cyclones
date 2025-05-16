# %%
from os.path import basename, join, exists

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
    model_dir = r'p:\11210471-001-compass\02_Models\somerset\SomersetLevels\sfincs'
    config_file = r'c:\CODE\COMPASS\compound-flooding-tropical-cyclones\Workflows\05_config_models\02_sfincs\sfincs_base_build.yml'
    data_cats = [
        r'c:\CODE\COMPASS\compound-flooding-tropical-cyclones\Workflows\03_data_catalogs\datacatalog_general.yml',
        r'c:\CODE\COMPASS\compound-flooding-tropical-cyclones\Workflows\03_data_catalogs\datacatalog_SFINCS_obspoints.yml',
        r'c:\CODE\COMPASS\compound-flooding-tropical-cyclones\Workflows\03_data_catalogs\datacatalog_SFINCS_coastal_coupling.yml',
        ]
    #bbox =[-3.2913,50.9637,-2.5063,51.3508]
    bbox =[-3.16169,51.06687,-2.867119,51.258058]
    #bbox =[36.7,-18.35,37.41,-17.64]
    bathy = 'gebco'
    #bathy = 'emodnet_bathy_E4_2018_msl'
    dfm_coastal_mask = 'coastal_coupling_msk_SMST'

    # data_cats = [
    #     r'c:\Git_repos\COMPASS\Workflows\data_catalogs\datacatalog_general.yml',
    #     r'c:\Git_repos\COMPASS\Workflows\data_catalogs\datacatalog_SFINCS_obspoints.yml',
    #     r'c:\Git_repos\COMPASS\Workflows\data_catalogs\datacatalog_SFINCS_coastal_coupling.yml',
    # ]

#%%
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

opt['setup_mask_active']['mask'] = region
opt['setup_mask_active']['mask_buffer'] = 2000
opt['setup_mask_active']['exclude_mask'] = dfm_coastal_mask

opt['setup_mask_bounds']['include_mask'] = dfm_coastal_mask

#%%
mod = SfincsModel(
    root=model_dir, data_libs=data_cats, mode="w+", logger=logger, **kwargs
)

# %% BUILD MODEL
mod.build(region={"geom": region}, opt=opt)

# %%

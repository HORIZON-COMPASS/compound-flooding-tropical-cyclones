# %%
from os.path import basename, join, exists

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
else:
    model_dir = r'p:\11210471-001-compass\02_Models\quelimane\Freddy2\sfincs'
    config_file = r'c:\Git_repos\COMPASS\Workflows\config_sfincs\sfincs_base_build.yml'
    data_cats = [
        r'c:\Git_repos\COMPASS\Workflows\data_catalogs\datacatalog_general.yml',
        r'c:\Git_repos\COMPASS\Workflows\data_catalogs\datacatalog_SFINCS_obspoints.yml',
        r'c:\Git_repos\COMPASS\Workflows\data_catalogs\datacatalog_SFINCS_coastal_coupling.yml',
    ]
    bbox ="[36.7,-18.35,37.41,-17.64]"


if not exists(model_dir):
    os.mkdir(model_dir)

# model and data paths/
logger = setuplog("update", join(model_dir, "hydromt.log"), log_level=10)
opt = configread(config_file, abs_path=True)  # read settings from ini file
kwargs = opt.pop("global", {})

#%%
region = get_local_vector_data(
    file = 'basin_atlas_level12_v10', #TODO change to 'merit_hydro_basins....' 
    bbox = bbox,
    data_cat = data_cats[0],
)
#%%
opt['setup_mask_active']['mask'] = region
#%%
mod = SfincsModel(
    root=model_dir, data_libs=data_cats, mode="w+", logger=logger, **kwargs
)

# %% BUILD MODEL
mod.build(region={"geom": region}, opt=opt)

# %%

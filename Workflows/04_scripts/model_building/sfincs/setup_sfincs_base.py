# %%
from os.path import join, exists
import geopandas as gpd
from hydromt.config import configread
import ast
from hydromt.log import setuplog
import hydromt
from hydromt_sfincs import SfincsModel

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
    file = 'basin_atlas_level12_v10',
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
# Include extra polygon for Sofala region
if region_name == 'sofala':
    opt2 = {
        'setup_mask_active': {
        'include_mask': 'sofala_incl_polygon',
        'reset_mask': False
        }}
    
    mod.update(
        write=True,
        opt=opt2
        )
else:
    pass

# %%
# Set up river depths using r+ mode 
# mod = SfincsModel(
#     root=model_dir, data_libs=data_cats, mode="r+", logger=logger, **kwargs
# ) 
# mod.read()
# gdf_powlaw = gpd.read_file(join(model_dir, "gis", "rivers_inflow.geojson"))
# depth_powlaw = hydromt.workflows.rivers.river_depth(data= gdf_powlaw,
#                                                     method = "powlaw",
#                                                     qbankfull_name= "qbankfull",
#                                                     hc = 0.27,
#                                                     hp = 0.30,)
# gdf_powlaw["rivdph"] = depth_powlaw
# gdf_powlaw.to_file(join(model_dir, "gis", "rivers_inflow_with_depth.geojson"))

# #%%
# # Can this be based more on the config?
# mod.setup_subgrid(datasets_dep=opt['setup_subgrid']['datasets_dep'] + [{'elevtn': bathy, 'reproj_method': 'bilinear'}],
#                   datasets_rgh=opt['setup_subgrid']['datasets_rgh'],
#                   datasets_riv=[{"centerlines": gdf_powlaw, # does this increase both width and depth?
#                                  "manning": 0.035, 
#                                  # "rivdph": 3 # to fill missing values
#                                  }],
#                   write_dep_tif = True,
#                   write_man_tif = True,
#                   nr_subgrid_pixels = 10, # increase when burning in rivers
#                   nrmax = 5000, 
#                   )

# #%%
# mod.write()

# %%

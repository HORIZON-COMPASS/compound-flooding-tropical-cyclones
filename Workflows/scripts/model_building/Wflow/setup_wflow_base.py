# %%
from os.path import basename, join

import pandas as pd
from hydromt.config import configread
from hydromt.log import setuplog
from hydromt_wflow import WflowModel

# %%
# model and data paths/
logger = setuplog("update", "./hydromt.log", log_level=10)
yml_file = "wflow_build_sofala.yml"

opt = configread(yml_file, abs_path=True)  # read settings from ini file
kwargs = opt.pop("global", {})
model_root = "wflow_sofala"
bbox = [36.68,-18.35,37.41,-17.64]
data_cat = r'c:\Git_repos\COMPASS\SFINCS\datacatalog_general.yml'

mod = WflowModel(
    root=model_root, data_libs=[data_cat], mode="w+", logger=logger, **kwargs
)

# %% BUILD MODEL
mod.build(region={"basin": r"c:\Git_repos\COMPASS\SFINCS\sfincs_sofala\computations\sfincs_CLI_Freddy2\gis\region.geojson"}, opt=opt)

# %%

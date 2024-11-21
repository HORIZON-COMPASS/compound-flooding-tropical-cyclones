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
bbox = [36.7,-18.35,37.41,-17.64]

mod = WflowModel(
    root=model_root, data_libs=["deltares_data"], mode="w+", logger=logger, **kwargs
)

# %% BUILD MODEL
mod.build(region={"basin": bbox, "outlets": True}, opt=opt)

# %%

# %%
from os.path import basename, join, exists

import pandas as pd
from hydromt.config import configread
from hydromt.log import setuplog
from hydromt_wflow import WflowModel

# %%
if "snakemake" in locals():
    model_dir = snakemake.params.model_dir
    config_file = snakemake.inputs.config_file
    data_cat = snakemake.params.data_cat
    bbox = snakemake.params.arg_bbox
    region_geom = snakemake.inputs.region_geom

if not exists(model_dir):
    os.mkdir(model_dir)
# model and data paths/
logger = setuplog("update", join(model_dir, "hydromt.log"), log_level=10)


opt = configread(config_file, abs_path=True)  # read settings from ini file
kwargs = opt.pop("global", {})


mod = WflowModel(
    root=model_dir, data_libs=[data_cat], mode="w+", logger=logger, **kwargs
)

# %% BUILD MODEL
mod.build(region={"basin": region_geom}, opt=opt)

# %%

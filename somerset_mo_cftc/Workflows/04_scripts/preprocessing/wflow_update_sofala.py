# %%
from os.path import basename, join

import pandas as pd
from hydromt.config import configread
from hydromt.log import setuplog
from hydromt_wflow import WflowModel

# %%
# model and data paths/
logger = setuplog("update", "./hydromt.log", log_level=10)
wflow_root = r"c:\Git_repos\COMPASS\Wflow\wflow_sofala"
event = "freddy2"
# %% Setup forcing
mod = WflowModel(
    root=wflow_root,
    data_libs=[r"c:\Git_repos\COMPASS\SFINCS\datacatalog_general.yml"],
    mode="r",
    logger=logger,
)
yml_file_forcing = r"c:\Git_repos\COMPASS\Wflow\wflow_update_sofala.yml"
opt = configread(yml_file_forcing, abs_path=True)  # read settings from ini file
mod.update(join(wflow_root, "event", event), forceful_overwrite=True, opt=opt)

# %% change output settings
# setting_toml = {
#     "netcdf.path": f"output_scalar.nc",
#     "netcdf.variable": [
#         {"name": "Q_src", "map": "gauges_src", "parameter": "lateral.river.q_av"}
# #     ],
# # }
# for option in setting_toml:
#     mod.set_config(option, setting_toml[option])
# mod.write_config()

#

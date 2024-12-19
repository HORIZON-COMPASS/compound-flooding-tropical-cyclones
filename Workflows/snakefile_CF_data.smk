#Creating a start for workflow
import os
from os.path import join
from snakemake.io import Wildcards
import pandas as pd
# from snakemake.utils import Paramspace


### Create the Paramspace object ###
# paramspace = Paramspace(pd.read_csv("config/run_names.csv", sep=","), filename_params="*", param_sep="~", filename_sep="-")

curdir = os.getcwd()
if os.name == 'nt': #Running on windows
    root_dir = join("p:/",config['root_dir'])
elif os.name == "posix": #Running on linux
    root_dir = join("/p", config['root_dir'])

def get_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return "data_catalogs/datacatalog_general.yml"
    elif os.name == "posix": #Running on linux
        return "data_catalogs/datacatalog_general___linux.yml"

def get_CF_datacatalog(wildcards):
    if os.name == 'nt': #Running on windows
        return "data_catalogs/datacatalog_CF_forcing.yml"
    elif os.name == "posix": #Running on linux
        return "data_catalogs/datacatalog_CF_forcing___linux.yml"

# Define wildcards from config file
precip_name=config["precip_name"]
CF_value_rain=config["CF_value_rain"]
tc_name=config["tc_name"]
CF_value_wind=config["CF_value_wind"]
# tidemodel = config["tidemodel"]
# SLR_value = config["SLR_value"]

rule all: 
    input: 
        expand(join(root_dir, "01_Data", "counterfactuals", "wind", "tc_{tc_name}_{CF_value_wind}.spw"), 
               tc_name=tc_name, CF_value_wind=CF_value_wind)
        # expand(join(root_dir, "01_Data", "counterfactuals", "precipitation", "{precip_name}_{CF_value_rain}.nc"), 
        #        precip_name=precip_name, CF_value_rain=CF_value_rain)

# rule create_SFINCS_base_model:
    # input:
    # output:
        # fn_base_model = "../Idai/runs/SFINCS/" + "base_model"

# Creat CF rainfall dataset and add to CF data catalog
# rule create_CF_rainfall:
#     params:
#         start_date=config["start_date"],
#         end_date=config["end_date"],
#         bbox=config["bbox"],
#         data_cat=get_datacatalog,
#         CF_data_cat=get_CF_datacatalog,
#         precip_name=lambda wildcards: wildcards.precip_name,
#         CF_value_rain=lambda wildcards: wildcards.CF_value_rain
#     output:
#         CF_rainfall=join(root_dir, "01_Data", "counterfactuals", "precipitation", "{precip_name}_{CF_value_rain}.nc")
#     script:
#         join(curdir,"scripts","preprocessing","CF_data","precip_CFs.py")  # Path to the CF script

### Other CF rules for wflow & DFM ###
# rule create_CF_noSLR:
#     input:
#         data_cat="../Idai/input_data/factual/datacatalog_general.yml",
#         CF_data_cat="../Idai/input_data/counterfactual/datacatalog_CF.yml"
#     params:
#         tidemodel = lambda wildcards: wildcards.tidemodel
#         SLR_value = lambda wildcards: wildcards.SLR_value # Remove SLR from Z0 (= Mean Sea Level offset) tide component
#     output:
#         CF_SLR="../Idai/input_data/counterfactual/SLR/{tidemodel}_{SLR_value}.nc"
#     script:
#         "../Idai/preprocessing/scripts/sealevel_CFs.py" 

# rule create_CF_wind:
#     params:
#         start_date=config["start_date"],
#         end_date=config["end_date"],
#         # tc_name=lambda wildcards: wildcards.tc_name,
#         tc_name=lambda wildcards: wildcards.tc_name,
#         # CF_value_wind=lambda wildcards: wildcards.CF_value_wind,
#         CF_value_wind=lambda wildcards: wildcards.CF_value_wind,
#     output:
#         CF_wind=join(root_dir, "01_Data", "counterfactuals", "wind", "tc_{tc_name}_{CF_value_wind}.spw")
#     script:
#         join(curdir,"scripts","preprocessing","CF_data","wind_CFs.py") 

# rule update_SFINCS_precip: #Updating the yml file #Creating the new rainfall
#     input:
#         fn_base_model = "../Idai/runs/SFINCS/" + "base_model" + "sfincs.inp"
#         updated_yml
#     output:
#         fn_updated_model_precip
#     script:
#         #

# rule update_SFINCS_waterlevel: 
#     input:
#         fn_updated_model_precip 
#         updated_yml
#     output:
#         fn_updated_model_precip_wl
#     script:
#         #

# rule update_SFINCS_discharge: #Updating the yml file #Creating the new rainfall
#     input:
#         fn_updated_model_precip_wl 
#         updated_yml
#     output:
#         fn_updated_model_precip_wl_dis
#     script:
#         #

# rule run_updated_SFINCS:
#     input:
#         fn_updated_model_precip_wl_dis
#     params:
#         sfincs_params = paramspace.instance
#     output:
#         "runs/SFINCS/" f"{paramspace.wildcard_pattern}"

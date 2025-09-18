#%% use pixi environment compass-wflow
# Load modules
import os
from hydromt_sfincs import SfincsModel

wind_forcing         = 'era5_hourly_spw_IBTrACS'
precip_forcing       = 'era5_hourly_zarr'
tidemodel            = 'GTSMv41' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap, tpxo80_opendap
datacat              = [
    '../../Workflows/03_data_catalogs/datacatalog_general.yml',
    '../../Workflows/03_data_catalogs/datacatalog_SFINCS_obspoints.yml',
    '../../Workflows/03_data_catalogs/datacatalog_SFINCS_coastal_coupling.yml',
    '../../Workflows/03_data_catalogs/datacatalog_CF_forcing.yml'
    ]
CF_SLR_txt           = "0"
CF_wind_txt          = "0"
CF_rain_txt          = "0"

# Factual model
model_name           = f"event_tp_{precip_forcing}_CF{CF_rain_txt}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}"
dir_run              = f"../data/sfincs/{model_name}"
outfile_png          = f"../figures/fS1.png"
outfile_pdf          = f"../figures/fS1.pdf"

# select the model and datacatalog
sfincs_root = dir_run
mod = SfincsModel(sfincs_root, data_libs=datacat, mode="r")

# reading in the model results
mod.read_results()

# Save as png
fig, ax = mod.plot_basemap(
    fn_out=os.path.join(os.path.abspath(os.path.dirname(outfile_png)),os.path.basename(outfile_png)), 
    plot_geoms=True, 
    figsize=(8, 6))

# And as pdf
fig, ax = mod.plot_basemap(
    fn_out=os.path.join(os.path.abspath(os.path.dirname(outfile_pdf)),os.path.basename(outfile_pdf)), 
    plot_geoms=True, 
    figsize=(8, 6))
# %%

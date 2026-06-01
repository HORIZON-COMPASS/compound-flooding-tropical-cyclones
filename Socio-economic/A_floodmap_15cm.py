#%% Use pixi environment compass-sfincs 
# Load modules
import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from os.path import join
from hydromt_sfincs import SfincsModel, utils

#%% Set parameters
region               = "sofala"
tc_name              = "Idai"
wind_forcing         = 'era5_hourly_spw_IBTrACS'
precip_forcing       = 'era5_hourly_zarr'
tidemodel            = 'GTSMv41' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap, tpxo80_opendap
datacat              = [
    '../Workflows/03_data_catalogs/datacatalog_general.yml',
    '../Workflows/03_data_catalogs/datacatalog_SFINCS_obspoints.yml',
    '../Workflows/03_data_catalogs/datacatalog_SFINCS_coastal_coupling.yml',
    '../Workflows/03_data_catalogs/datacatalog_CF_forcing.yml'
    ]
CF_SLR_txt           = "-0.1"
CF_wind_txt          = "-5"
CF_rain_txt          = "-8"
model_name           = f"event_tp_{precip_forcing}_CF{CF_rain_txt}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}"
dir_run              = f"p:/11210471-001-compass/03_runs/{region}/{tc_name}/sfincs/{model_name}"
mapfile              = f"{dir_run}/sfincs_map.nc"
outfile              = f"{dir_run}/plot_output/sfincs_basemap.png"
floodmap             = f"{dir_run}/plot_output/floodmap_15cm.tif"

#%%
print("------- Checking what we got ------")
print("Model run directory: ", dir_run)
print("mapfile: ", mapfile)
print("Output figure basemap: ", outfile)
print("Floodmap output: ", floodmap)

#%%
# select the model and datacatalog
sfincs_root = dir_run
mod = SfincsModel(sfincs_root, data_libs=datacat, mode="r")

# reading in the model results
mod.read_results()

#%%
### PLOT BASEMAP
fig, ax = mod.plot_basemap(
    fn_out=os.path.join(os.path.abspath(os.path.dirname(outfile)),os.path.basename(outfile)), 
    plot_geoms=True, 
    figsize=(8, 6))

#%%
### PLOT FORCING
_ = mod.plot_forcing(
    fn_out = os.path.join(os.path.abspath(os.path.dirname(outfile)),'sfincs_forcing.png'))

#%%
# compute the maximum water level over all time steps
da_zsmax = mod.results["zsmax"].max(dim="timemax")

# select our highest-resolution elevation dataset
depfile = join(dir_run, "subgrid", "dep_subgrid.tif")
da_dep = mod.data_catalog.get_rasterdataset(depfile)

# we set a threshold to mask minimum flood depth
hmin = 0.15

# Downscale the floodmap
da_hmax = utils.downscale_floodmap(
    zsmax=da_zsmax,
    dep=da_dep,
    hmin=hmin,
    reproj_method = "bilinear",
    # floodmap_fn=floodmap # uncomment to save to <mod.root>/floodmap.tif)
)

# we use the GSWO dataset to mask permanent water by first reprojecting it to the subgrid of hmax
gswo = mod.data_catalog.get_rasterdataset("gswo", geom=mod.region, buffer=1000)
gswo_mask = gswo.raster.reproject_like(da_hmax, method="max")

# permanent water where water occurence > 5%
da_hmax_masked = da_hmax.where(gswo_mask <= 5)

# save the masked floodmap
# da_hmax_masked.raster.to_raster(os.path.join(os.path.abspath(os.path.dirname(outfile)),'sfincs_output_hmax_AllTime.tif'))
da_hmax_masked.raster.to_raster(os.path.join(os.path.abspath(os.path.dirname(outfile)), floodmap))

# %%

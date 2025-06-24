#%% 
# Load modules
import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from os.path import join
from hydromt_sfincs import SfincsModel, utils

#%%
# We read the snakemake parameters
if "snakemake" in locals():
    mapfile              = snakemake.input.mapout
    outfile              = snakemake.output.figure
    floodmap             = snakemake.output.floodmap
    dir_model_no_forcing = snakemake.params.dir_model_no_forcing
    dir_run              = snakemake.params.dir_run
    datacat              = snakemake.params.datacat
else:
    region               = "sofala"
    tc_name              = "Idai"
    wind_forcing         = 'spw_IBTrACS'
    precip_forcing       = 'era5_hourly'
    tidemodel            = 'GTSMv41opendap' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap, tpxo80_opendap
    datacat              = [
        '../../../03_data_catalogs/datacatalog_general.yml',
        '../../../03_data_catalogs/datacatalog_SFINCS_obspoints.yml',
        '../../../03_data_catalogs/datacatalog_SFINCS_coastal_coupling.yml',
        '../../../03_data_catalogs/datacatalog_CF_forcing.yml'
        ]
    CF_SLR_txt           = "0"
    CF_wind_txt          = "-10"
    CF_rain_txt          = "0"
    model_name           = f"event_tp_{precip_forcing}_CF{CF_rain_txt}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}"
    dir_run            = f"p:/11210471-001-compass/03_runs/{region}/{tc_name}/sfincs/{model_name}"
    mapfile              = f"{dir_run}/sfincs_map.nc"
    outfile              = f"{dir_run}/plot_output/sfincs_basemap.png"
    floodmap             = f"{dir_run}/plot_output/floodmap.tif"

#%%
print("------- Checking what we got ------")
print("Model run directory: ", dir_run)
print("mapfile: ", mapfile)
print("Output figure basemap: ", outfile)
#%%
# select the model and datacatalog
sfincs_root = dir_run
mod = SfincsModel(sfincs_root, data_libs=datacat, mode="r")

#%%
# mod.data_catalog.from_yml(datacat)
#mod.data_catalog.from_yml(datacat,root='p:/') # ---> FOR WINDOWS
#%%
# select our highest-resolution elevation dataset
depfile = join(dir_run, "subgrid", "dep_subgrid.tif")
da_dep = mod.data_catalog.get_rasterdataset(depfile)

# read global surface water occurance (GSWO) data to mask permanent water
gswo = mod.data_catalog.get_rasterdataset("gswo", geom=mod.region, buffer=1000)

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
### PLOT MAX INUNDATION
# compute the maximum water level over all time steps
da_zsmax = mod.results["zsmax"].max(dim="timemax")

# we set a threshold to mask minimum flood depth
hmin = 0.05

# Downscale the floodmap
da_hmax = utils.downscale_floodmap(
    zsmax=da_zsmax,
    dep=da_dep,
<<<<<<< HEAD
    hmin=hmin,
    floodmap_fn=floodmap # uncomment to save to <mod.root>/floodmap.tif)
)
=======
    hmin=hmin, 
    reproj_method ="bilinear")
>>>>>>> origin/main

# we use the GSWO dataset to mask permanent water by first reprojecting it to the subgrid of hmax
gswo_mask = gswo.raster.reproject_like(da_hmax, method="max")

# permanent water where water occurence > 5%
da_hmax_masked = da_hmax.where(gswo_mask <= 5)

# basemap plot the masked hmax on top
fig, ax = mod.plot_basemap(
    fn_out=None,
    figsize=(8, 6),
    variable=da_hmax_masked,
    plot_bounds=False,
    plot_geoms=False,
    bmap="sat",
    zoomlevel=11,
    vmin=0,
    vmax=2.5,
    cbar_kwargs={"shrink": 0.6, "anchor": (0, 0)},
)

tstart = np.datetime_as_string(mod.results['zs'].time.values[0],'m').replace('T',' ')
tend = np.datetime_as_string(mod.results['zs'].time.values[-1],'m').replace('T',' ')
ax.set_title(f"SFINCS masked maximum water depth \n Period: {tstart} to {tend}")
fig.savefig(os.path.join(os.path.abspath(os.path.dirname(outfile)),f'sfincs_output_hmax_AllTime.png'))
del da_zsmax

<<<<<<< HEAD
#%%
=======
# save as raster
da_hmax_masked.raster.to_raster(os.path.join(os.path.abspath(os.path.dirname(outfile)),'sfincs_output_hmax_AllTime.tif'))


>>>>>>> origin/main
### PLOT MAX INUNDATION PER TIMEMAX TIMESTAMP
# repeat the same steps as above, but for individual timesteps of timemax variable
for ii,timestamp in enumerate(mod.results['zsmax'].timemax.values):
    da_zsmax = mod.results["zsmax"].isel(timemax=ii)
    da_hmax = utils.downscale_floodmap(
        zsmax=da_zsmax,
        dep=da_dep,
        hmin=hmin)
    
    # permanent water where water occurence > 5%
    da_hmax_masked = da_hmax.where(gswo_mask <= 5)

    # basemap plot the masked hmax on top
    fig, ax = mod.plot_basemap(
        fn_out=None,
        figsize=(8, 6),
        variable=da_hmax_masked,
        plot_bounds=False,
        plot_geoms=False,
        bmap="sat",
        zoomlevel=11,
        vmin=0,
        vmax=2.5,
        cbar_kwargs={"shrink": 0.6, "anchor": (0, 0)},
    )
    if ii==0:
        tstart = np.datetime_as_string(mod.results['zs'].time.values[0],'m').replace('T',' ')
    else:
        tstart = np.datetime_as_string(mod.results['zsmax'].timemax.values[ii-1],'m').replace('T',' ')
    tend = np.datetime_as_string(mod.results['zsmax'].timemax.values[ii],'m').replace('T',' ')
    ax.set_title(f"SFINCS masked maximum water depth \n Period: {tstart} to {tend}")

    figname_ext = f'period_{tstart}_to_{tend}'
    for r in ((":", ""), (" ", "_"), ("-", "")):
        figname_ext = figname_ext.replace(*r)
    fig.savefig(os.path.join(os.path.abspath(os.path.dirname(outfile)),f'sfincs_output_hmax_{figname_ext}.png'))

<<<<<<< HEAD
#%%
=======
    # save as raster
    da_hmax_masked.raster.to_raster(os.path.join(os.path.abspath(os.path.dirname(outfile)),f'sfincs_output_hmax_{figname_ext}.tif'))


>>>>>>> origin/main
### PLOT TIMESERIES OUTPUT POINTS
if os.path.exists(join(dir_run,"sfincs_his.nc")):
    hisfile = os.path.join(dir_run,"sfincs_his.nc")
    ds_his = xr.open_dataset(hisfile)
    ds_his["station_id"] = ds_his["station_id"].astype(int)

    for ii,loc_id in enumerate(ds_his['station_id'].values):
        fig, ax = plt.subplots(1,1,figsize=(10,5))
        ax.plot(ds_his.time,ds_his['point_zs'].isel(stations=ii),color='r',label=f'')
        ax.grid()
        ax.set_title(f"Timeseries of water levels \n Location: {loc_id}")
        ax.set_ylabel("Inundation height [m]")
        fig.savefig(os.path.join(os.path.abspath(os.path.dirname(outfile)),f'sfincs_output_TS_loc_{loc_id}.png'))
    else:
        print("No sfincs_his.nc file found in model run directory. Skipping timeseries plots.")
#%%
# This code is based on that of DirkEilander. (2022). DirkEilander/compound_flood_modelling: revised paper (Version v2). Zenodo. https://doi.org/10.5281/zenodo.7274465

# Load modules
import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from hydromt_sfincs import SfincsModel, utils
from os.path import join, isfile, basename, isdir, dirname
import glob
import hydromt

#%%
# Set variable of the runs
region               = "sofala"
tc_name              = "Idai"
wind_forcing         = 'spw_IBTrACS'
precip_forcing       = 'era5_hourly_zarr'
tidemodel            = 'GTSMv41opendap' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap, tpxo80_opendap
datacat              = [
    '../../../03_data_catalogs/datacatalog_general.yml',
    '../../../03_data_catalogs/datacatalog_SFINCS_obspoints.yml',
    '../../../03_data_catalogs/datacatalog_SFINCS_coastal_coupling.yml',
    '../../../03_data_catalogs/datacatalog_CF_forcing.yml'
    ]
CF_SLR_txt           = "0"
CF_wind_txt          = "0"
CF_rain_txt          = "0"
model_name           = f"event_tp_{precip_forcing}_CF{CF_rain_txt}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}"
root                 = f"p:/11210471-001-compass/03_runs/{region}/{tc_name}/sfincs/{model_name}"
mapfile              = f"{root}/sfincs_map.nc"
outfile              = f"{root}/plot_output/sfincs_basemap.png"
floodmap             = f"{root}/plot_output/floodmap.tif"
subgrid              = f"{root}/subgrid/dep_subgrid.tif"


#%%
mod0 = SfincsModel(root, mode='r')
mod0.read_results()

#%%
# read and resample images to model domain
ddir = r'p:/11210471-001-compass/01_Data/Satellite/Idai/3_eo_rapid'
sat_dir = r'p:/11210471-001-compass/01_Data/Satellite/Idai/3_eo_rapid/reprojected'
dates = ['20190319', '20190320']
types = ['flooding', 'Dry']

hmin = 0.05

da_dep = mod0.data_catalog.get_rasterdataset(subgrid)
# compute the maximum over all time steps
da_zsmax = mod0.results["zsmax"].max(dim="timemax")
     
# downscale the floodmap
da_hmax = utils.downscale_floodmap(
    zsmax=da_zsmax,
    dep=da_dep,
    hmin=hmin,
    # floodmap_fn=join(sfincs_root, "gis/floodmap.tif") # uncomment to save floodmap to <mod.root>/floodmap.tif
    )

rm_fns = []
for date in dates:
    for name in types:
        fns = glob.glob(join(ddir, date, f'{name}*.tif'))
        da_lst = []
        for fn in fns:
            print(fn)
            da_obs0 = hydromt.open_raster(fn).load().astype(np.int8)
            da_obs0.raster.set_nodata(-1)
            try:
                da_obs0 = da_obs0.raster.reproject_like(da_hmax, method='max')
                da_lst.append(da_obs0)
            except IndexError:
                da_obs0.close()
                rm_fns.append(fn)
                print('out of bounds')
                pass
            
        if len(da_lst) > 0:
            print(f'concatenate {len(da_lst)} files')
            da = xr.concat(da_lst, dim='img').max('img').load().astype(np.uint8)
            da = da.where(da==1,0)
            da.raster.set_nodata(0)
            da.raster.to_raster(join(sat_dir, f'{name.lower()}_{date}.tif'), compress='deflate')

if len(rm_fns) > 0:
    print(f'deleting {len(rm_fns)} files')
    for fn in rm_fns:
        os.unlink(fn)


# %% TODO ask Dirk about the data
# read and resample processed sentinel 1 images to model domain
# ddir = r'p:/11210471-001-compass/01_Data/Satellite/Idai/3_eo_sentinel1'
# events = {
#     'idai': ['20190319', '20190320']
# }
# fns = glob.glob(join(ddir, f'sofala*.tif'))
# for event, dates in events.items():
#     da_obs_lst = []
#     for date in dates:
#         name = f'sofala_floodmap_RefinedLee_{date}_10m.tif'
#         fn = join(ddir, name)
#         da_obs0 = hydromt.open_raster(fn).load().astype(np.int8)
#         da_obs0.raster.set_nodata(-1)
#         da_obs0 = da_obs0.raster.reproject_like(mod0, method='max')
#         # da_obs0 = da_obs0.where(da_obs0==1,0)
#         # da_obs0.raster.set_nodata(-1)
#         da_obs_lst.append(da_obs0)
#         da_obs0.raster.to_raster(join(root, 'gis', f'sofala_floodmap_{date}.tif'), compress='deflate')
#     da_obs_max = xr.concat(da_obs_lst, dim='time').max('time')
#     da_obs0.raster.to_raster(join(root, 'gis', f'sofala_floodmap_{event}_max.tif'), compress='deflate')

#%%
event = {
    'idai': ['20190319', '20190320']
}
sfx_skill = {}
# postfix = '_100m'

# select our highest-resolution elevation dataset
depfile = join(root, "subgrid", "dep_subgrid.tif")
da_dep = mod0.data_catalog.get_rasterdataset(depfile)

# read global surface water occurance (GSWO) data to mask permanent water
gswo = mod0.data_catalog.get_rasterdataset("gswo", geom=mod0.region, buffer=1000)

gswo_mask = gswo.raster.reproject_like(da_hmax, method="max")

# permanent water where water occurence > 5%
msk = da_hmax.where(gswo_mask <= 5)


for event, dates in event.items():
    da_sfx = xr.open_dataarray(join(root, f'sfincs_map.nc'))
    # da_cmf = xr.open_dataarray(join(mdir1, f'flddph_{event}_v2.nc'))

    da_obs_lst = []
    for date in dates:
        # read observations
        da_obs = hydromt.open_raster(join(sat_dir, f'flooding_{date}.tif'))
        # da_obs = hydromt.open_raster(join(roots[postfix], 'gis', f'sofala_floodmap_{date}.tif'))
        da_obs_lst.append(da_obs)

        # validate SFINCS
        da_sim = da_sfx.sel(time=date).squeeze().raster.flipud()
        da_skill, da_cm = skill(da_sim, da_obs, msk, hmin=hmin)
        df_skill = da_skill.reset_coords(drop=True).to_dataframe()
        sfx_skill[date] = df_skill.loc[sfx_rm.keys()].rename(sfx_rm).drop(columns='E')

    # max exent
    da_obs = hydromt.open_raster(join(sat_dir, f'sofala_floodmap_{event}_max.tif'))
    # validate SFINCS
    da_sim = da_sfx.sel(time=slice(*tslice_max[event])).max('time').raster.flipud()
    da_skill, da_cm = skill(da_sim, da_obs, msk[postfix], hmin=hmin)
    df_skill = da_skill.reset_coords(drop=True).to_dataframe()
    sfx_skill[event] = df_skill.loc[sfx_rm.keys()].rename(sfx_rm).drop(columns='E')

    # validate CMF
    da_sim = da_cmf.sel(time=slice(*tslice_max[event])).max('time')
    da_skill, da_cm = skill(da_sim, da_obs, msk[postfix], hmin=hmin)
    df_skill = da_skill.reset_coords(drop=True).to_dataframe()
    cmf_skill[event] = df_skill.loc[cmf_rm.keys()].rename(cmf_rm).drop(columns='E')
   

dfs = []
for date in sfx_skill.keys():
    df1 = pd.concat([cmf_skill[date],sfx_skill[date]],axis=1,keys=['CMF', 'SFINCS']).swaplevel(0,1,axis=1).sort_index(axis=0).sort_index(axis=1)
    dfs.append(df1)
df1 = pd.concat(dfs, axis=1, keys=sfx_skill.keys())
# df1
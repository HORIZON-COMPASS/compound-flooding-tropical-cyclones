# %% Script of jupyter notebook made by Natalia Aleksandrova - edits by Doris Vertegaal
import sys
from pathlib import Path

# Add current script directory to sys.path
sys.path.append(str(Path(__file__).resolve().parent))

from skill import skill #source: @DirkEilander, https://github.com/DirkEilander/compound_flood_modelling , Nov 2022
import numpy as np
from os.path import join
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
from string import ascii_lowercase as abcd
import hydromt
import rioxarray as rxr
import rasterio as rio
from rasterio import features
from rasterio.features import sieve
from rasterio.crs import CRS
import geopandas as gpd
import platform
from shapely.geometry import box
from rasterio.features import shapes
from shapely.geometry import shape
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import os

import platform
prefix = "p:/" if platform.system() == "Windows" else "/p/"

def lat_formatter(x, pos):
    direction = 'N' if x >= 0 else 'S'
    return f"{abs(x):.1f}°{direction}"
def lon_formatter(x, pos):
    direction = 'E' if x >= 0 else 'W'
    return f"{abs(x):.1f}°{direction}"

def custom_formatter(value, pos=None):
    return f"{value:.1f}°"

# %% [markdown]
event = 'idai'

# %%
dir_obs_unosat = join(prefix,'11210471-001-compass','01_Data','Validation_UNOSAT')
dir_obs_cems = join(prefix,'11210471-001-compass', '01_Data','Validation_GFM')

rundirs = join(prefix,'11210471-001-compass','03_Runs')

# %%
event == 'idai'
mdir = join(rundirs, 'sofala','Idai','sfincs','event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_spw_IBTrACS_CF0')
file_floodmap_unosat = join(dir_obs_unosat, 'TC20190312MOZ_SHP', 'ST1_20190319_WaterExtent_ManicaSofalaProvinces.shp')
file_floodmap_cems = join(dir_obs_cems, 'Idai_2019','maximum_2019_03','maximum_flood_extent_2019-03-01_2019-03-30_beira_2025_06_26T11_41_09_531327.geojson')
sfincs_crs = '32736'
sfincs_utm = '36S'

# %%
hmin=0.05

# %%
print(os.getcwd())
if platform.system() == "Windows":
    dataCat = hydromt.data_catalog.DataCatalog(join('..','..','Workflows', "03_data_catalogs", "datacatalog_general.yml"))
else:
    dataCat = hydromt.data_catalog.DataCatalog(join('..','..','Workflows', "03_data_catalogs", "datacatalog_general___linux.yml"))

# %%
# Read model output data 
da_model = rxr.open_rasterio(join(prefix,'11210471-001-compass', '03_Runs', 'sofala', 'Idai', 'sfincs', 
                                  'event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0',
                                  'plot_output','floodmap.tif'))

# %%
# add crs information and convert to rioxarray dataset
da_model = da_model.sortby(["x", "y"])
da_model.rio.write_crs(f"epsg:{sfincs_crs}", inplace=True)
ds = da_model.to_dataset(name='model')

# %%
# read region mask
region = gpd.read_file(join(prefix,'11210471-001-compass', '03_Runs', 'sofala', 'Idai', 'sfincs',
                            'event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0', 'gis', 'region.geojson'))
region_wsg = region.to_crs("EPSG:4326")
# %%
# read permanent water mask
gswo = dataCat.get_rasterdataset("gswo", geom=region, buffer=1000)
gswo_mask = gswo.raster.reproject_like(da_model, method="max")
valid_mask = (gswo_mask <= 5).astype("uint8").squeeze()

# Extract shapes
shapes_gen = shapes(valid_mask.values, transform=valid_mask.rio.transform())
valid_polygons = [shape(geom) for geom, val in shapes_gen if val == 1]
gdf_valid = gpd.GeoDataFrame(geometry=valid_polygons, crs=gswo_mask.rio.crs)
gdf_valid = gdf_valid.to_crs(region_wsg.crs)

# %%
# Load raster flood map with rasterio as a basis for rasterization of the region and EO flood map vectors
raster = rio.open(join(prefix, '11210471-001-compass', '03_Runs', 'sofala', 'Idai', 'sfincs', 
                       'event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0', 'plot_output', 'floodmap.tif'), 'r+')
raster.crs = CRS.from_epsg(sfincs_crs)

# %%
# clean flood map - remove pixelation
msk = raster.read_masks()
msk2 = sieve(msk, size=150, connectivity=4)
ds['model_clean'] = ds['model'].where(msk2) 

# %%
# rasterize SFINCS region mask
rasterized = features.rasterize(region.geometry,
                                out_shape = raster.shape,
                                fill = 0,
                                out = None,
                                transform = raster.transform,
                                all_touched = False,
                                default_value = 1,
                                dtype = None)

ds['sfincs_mask'] = (('y', 'x'), rasterized)

# %%
if file_floodmap_unosat is not None:
    
    gpd_obs = gpd.read_file(file_floodmap_unosat)
    gpd_obs = gpd_obs.to_crs(sfincs_crs)

    # Get list of geometries for all features in vector file
    geom = [shapes for shapes in gpd_obs.geometry]; del gpd_obs


    # rasterize vector data based on the flood map from SFINCS
    rasterized = features.rasterize(geom,
                                    out_shape = raster.shape,
                                    fill = 0,
                                    out = None,
                                    transform = raster.transform,
                                    all_touched = False,
                                    default_value = 1,
                                    dtype = None)

    ds['obs_unosat'] = (('y', 'x'), rasterized)
    del rasterized, geom

    # Apply water mask
    ds['obs_unosat'] = ds['obs_unosat'].where(gswo_mask <= 5)

    # Apply region mask
    ds['obs_unosat'] = ds['obs_unosat'].where(ds.sfincs_mask) 

    print('UNOSAT flood map loaded.')
else:
    print('UNOSAT flood maps are not included in the analysis')

# %%
if file_floodmap_cems is not None:
    gpd_obs = gpd.read_file(file_floodmap_cems)
    gpd_obs = gpd_obs.to_crs(sfincs_crs)

    # Get list of geometries for all features in vector file
    geom = [shapes for shapes in gpd_obs.geometry]; del gpd_obs

    # rasterize vector data based on the flood map from SFINCS
    rasterized = features.rasterize(geom,
                                    out_shape = raster.shape,
                                    fill = 0,
                                    out = None,
                                    transform = raster.transform,
                                    all_touched = False,
                                    default_value = 1,
                                    dtype = None)

    ds['obs_cems'] = (('y', 'x'), rasterized)
    del rasterized

    # Apply water mask
    ds['obs_cems'] = ds['obs_cems'].where(gswo_mask <= 5)

    # Apply region mask
    ds['obs_cems'] = ds['obs_cems'].where(ds.sfincs_mask) 

    print('UNOSAT flood map loaded.')
else:
    print('CEMS flood maps are not included in the analysis')

# %%
# Plot settings
cm_dict = {
    1: ('false neg.', '#dd8452'),
    2: ('false pos.', '#c44e52'),
    3: ('true pos.', '#4c72b0'),
}
levels = [k for k,v in cm_dict.items()] + [4]
colors = [v[1] for k,v in cm_dict.items()]
cmap, norm = clrs.from_levels_and_colors(levels, colors)
ticklabs = [v[0] for k,v in cm_dict.items()]
ticks = np.array(levels[:-1])+np.diff(levels)/2.

props = dict( facecolor='w', lw=0, alpha=0.8)

# %%
results = []  # store tuples (source, da_skill, da_cm)

source_list = [x.replace('obs_','') for x in list(ds.keys()) if 'obs' in x]

for source in source_list:
    da_skill, da_cm = skill(ds['model_clean'], ds[f'obs_{source}'], ds['sfincs_mask'], hmin=0.3)
    da_cm = da_cm.load()
    da_skill = da_skill.load()
    results.append((source, da_skill, da_cm))


# %%
# Plot with validation metrics
source_list = [x.replace('obs_','') for x in list(ds.keys()) if 'obs' in x]

utm_crs = ccrs.UTM(zone=36, southern_hemisphere=True)  # or False if in north
ll_crs = ccrs.PlateCarree()

fig, axs = plt.subplots(
    figsize=(len(source_list)*4,6),
    nrows=1, ncols=len(source_list),
    subplot_kw={'projection': ll_crs},
    sharex = True, sharey=True
)
if len(source_list) > 1: axs = axs.flatten()

for ii, source in enumerate(source_list):
    if len(source_list) > 1:
        ax=axs[ii]
    else: 
        ax=axs

    # Calculate skill
    da_skill, da_cm = skill(ds['model_clean'], ds[f'obs_{source}'], ds['sfincs_mask'], hmin=0.3)
    df_skill = da_skill.reset_coords(drop=True).to_dataframe()

    da_cm = da_cm.load()
    da_skill = da_skill.load()
    hr, csi, fr = np.round(da_skill['H'].item(),2), np.round(da_skill['C'].item(),2), np.round(da_skill['F'].item(),2)
    cs = da_cm.where(da_cm>0).plot(ax=ax, cmap=cmap, norm=norm, transform=utm_crs, add_colorbar=False)
    ax.text(0.8, 0.85, f'C: {csi:.2f}\nH: {hr:.2f}\nF: {fr:.2f}', transform=ax.transAxes, bbox=props, fontsize=10)

    # Add model region
    gdf_valid.plot(ax=ax, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)
    region_wsg.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())
    
    # Set extent (based on actual lat/lon coordinates)
    minx, miny, maxx, maxy = region_wsg.bounds.minx.item(), region_wsg.bounds.miny.item(), region_wsg.bounds.maxx.item(), region_wsg.bounds.maxy.item()
    ax.set_extent([minx, maxx, miny, maxy], ccrs.PlateCarree())

    # Add gridlines and format tick labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.FixedLocator(np.arange(minx, maxx + 0.1, 0.2))
    gl.ylocator = mticker.FixedLocator(np.arange(miny, maxy + 0.1, 0.2))
    gl.xformatter = mticker.FuncFormatter(lon_formatter)
    gl.yformatter = mticker.FuncFormatter(lat_formatter)
    gl.right_labels = False
    gl.top_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}
    if ii == 1:  
        gl.left_labels = False  # disable y-axis labels

fig.subplots_adjust(wspace=0.04, hspace=0.06)
axs[0].set_title(f'{source_list[0].upper()}')
axs[1].set_title(f'{source_list[1].upper()}')

# Add a colorbar axis 
cbar_ax = fig.add_axes([0.93, 0.2, 0.015, 0.5])

# Draw the colorbar
cbar=fig.colorbar(cs, cax=cbar_ax, orientation='vertical', ticks=ticks)
cbar_ax.set_yticklabels(ticklabs, va='center', rotation=90)

# Subplot (a) - first plot
axs[0].text(0.0, 1.06, "(a)", transform=axs[0].transAxes,
             fontsize=11, fontweight='bold', va='top', ha='left')

# Subplot (b) - second plot
axs[1].text(0.0, 1.06, "(b)", transform=axs[1].transAxes,
             fontsize=11, fontweight='bold', va='top', ha='left')

fig.savefig("../figures/fS7.png",dpi=300, bbox_inches='tight')
fig.savefig("../figures/fS7.pdf",dpi=300, bbox_inches='tight')
# %%

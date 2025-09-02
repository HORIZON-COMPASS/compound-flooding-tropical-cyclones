#%%
print("Loading packages")
import os
from os.path import join
import xarray as xr
import numpy as np
import geopandas as gpd
from hydromt_sfincs import SfincsModel, utils
from hydromt import DataCatalog
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from pyproj import Transformer
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import contextily as ctx
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from shapely.geometry import box
from rasterio.features import shapes
from shapely.geometry import shape

import warnings
warnings.filterwarnings('ignore')

import platform
prefix = "p:/" if platform.system() == "Windows" else "/p/"

def lat_formatter(x, pos):
    return f"{abs(x):.1f}°"
def lon_formatter(x, pos):
    return f"{abs(x):.1f}°"

def custom_formatter(value, pos=None):
    return f"{value:.1f}°"

#%%
print("Loading Factual model")
# define model and data catalog file paths
factual_model_dir = os.path.join(prefix,"11210471-001-compass","03_Runs","sofala","Idai","sfincs","event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0")
cf_model_dir = os.path.join(prefix,"11210471-001-compass","03_Runs","sofala","Idai","sfincs","event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10")
# factual_model_dir = r"c:\Code\Paper_1\Tests\event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0"
if platform.system() == "Windows":
    datacat = [
        '../../Workflows/03_data_catalogs/datacatalog_general.yml'
        ]
else:
    datacat = [
        '../../Workflows/03_data_catalogs/datacatalog_general___linux.yml'
        ]

data_catalog = DataCatalog(data_libs = datacat)

#%%
# Load in model, model region, and buffer model region
mod = SfincsModel(factual_model_dir, data_libs=datacat, mode="r")

hisfile = os.path.join(factual_model_dir,"sfincs_his.nc")
ds_his = xr.open_dataset(hisfile, engine='netcdf4')

model_region_gdf = gpd.read_file(join(
    prefix, "11210471-001-compass", "03_Runs", "sofala", "Idai", "sfincs", 
    "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0", "gis", "region.geojson"
)).to_crs("EPSG:4326") 

#%%
print("Masking permanent water")
# we set a threshold to mask minimum flood depth
hmin = 0.05

# compute the maximum over all time steps
da_zsmax = mod.results["zsmax"].max(dim="timemax")
     
# downscale the floodmap
depfile  = join(factual_model_dir, "subgrid", "dep_subgrid.tif")
da_dep   = mod.data_catalog.get_rasterdataset(depfile)

da_hmax = utils.downscale_floodmap(
      zsmax=da_zsmax,
      dep=da_dep,
      hmin=hmin,
      )
    
# GSWO dataset to mask permanent water by first geprojecting it to the subgrid of hmax
sfincs_region = mod.region
projection    = mod.crs.to_epsg()
gwso_region   = data_catalog.get_rasterdataset("gswo", geom=sfincs_region, buffer=1000)
gswo_mask     = gwso_region.raster.reproject_like(da_hmax, method="max")
# permanent water where water occurence > 5%
da_hmax_masked = da_hmax.where(gswo_mask <= 5)

# Add the name attribute for identification
mod.results['hmax'] = da_hmax
mod.results['hmax_masked'] = da_hmax_masked

# Make own background shape
valid_mask = (gswo_mask <= 5).astype("uint8").squeeze()

# Extract shapes
shapes_gen = shapes(valid_mask.values, transform=valid_mask.rio.transform())
valid_polygons = [shape(geom) for geom, val in shapes_gen if val == 1]
gdf_valid = gpd.GeoDataFrame(geometry=valid_polygons, crs=gswo_mask.rio.crs)
gdf_valid = gdf_valid.to_crs(model_region_gdf.crs)

del da_hmax, da_zsmax, da_dep  # Clean up to free memory


#%%
print("Plotting factual hmax masked")
fig, ax = plt.subplots(nrows=1,
                       ncols=1,
                       figsize=(4,6),
                       dpi=300,
                       constrained_layout=True,
                       subplot_kw={"projection": ccrs.epsg(projection)})

# fig.suptitle("Factual Max Flood Depth", fontsize=11)

# Plot hmax masked
hmax = mod.results['hmax_masked']

vmin = float(hmax.min())
vmax = float(hmax.max())

im = hmax.plot.pcolormesh(
    ax=ax, cmap="Blues", vmin=vmin, vmax=3.5, add_colorbar=False
)

# gdf_valid.plot(ax=ax, color='#E0E0E0', transform=ccrs.PlateCarree(), zorder=0)
ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=10, crs=hmax.rio.crs, attribution=False)

# Add model region
# model_region_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=0.3, transform=ccrs.PlateCarree())

# Set extent (based on actual lat/lon coordinates)
minx, miny, maxx, maxy = model_region_gdf.bounds.minx.item(), model_region_gdf.bounds.miny.item(), model_region_gdf.bounds.maxx.item(), model_region_gdf.bounds.maxy.item()
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

# Titles
ax.set_title("", fontsize=8)

# Add shared colorbar with larger size and labels
cbar = fig.colorbar(im, ax=ax, orientation="vertical", 
                    fraction=0.035, pad=0.05)
cbar.set_label('Flood depth (m)', labelpad=5, fontsize=9)
cbar.ax.tick_params(labelsize=8)
    
fig.savefig("../figures/factual_hmax_masked.png", bbox_inches='tight', dpi=300)
plt.show()
# %%

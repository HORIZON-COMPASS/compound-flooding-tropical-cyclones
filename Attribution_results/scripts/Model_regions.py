#%% Import the necessary packages
# Use pixi environment compass-snake-dfm
import os
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx
import numpy as np
from matplotlib import cm
import matplotlib.lines as mlines
import cartopy.crs as ccrs
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from matplotlib.patches import ConnectionPatch

#%% Load TC track shapefiles as geopandas geodataframe - obtain from https://www.ncei.noaa.gov/products/international-best-track-archive (v4r01 SI points)
data_base = "../data/"

shapefile_path = os.path.join(data_base, "ibtracs/IBTrACS.SI.list.v04r01.points.shp")
gdf = gpd.read_file(shapefile_path)
tc_idai = gdf[gdf['SID'] == '2019063S18038']

# Normalize windspeed for point sizing (adjust scaling as needed)
tc_idai["size"] = (tc_idai["USA_WIND"] - tc_idai["USA_WIND"].min()) / (
    tc_idai["USA_WIND"].max() - tc_idai["USA_WIND"].min()
) * 138 + 19  # Scale between 20 and 145

cmap_idai = cm.get_cmap("Reds")

wind_min = tc_idai["USA_WIND"].min()
wind_max = tc_idai["USA_WIND"].max()

# The same colormap
cmap_idai = cm.get_cmap("Reds")

# Normalize function for colors (wind speed to 0-1)
norm = plt.Normalize(wind_min, wind_max)

# %% PLOTTING MODEL DOMAIN FIGURES FOR PAPER
# Getting the model regions
gdf_wflow = gpd.read_file(os.path.join(data_base, "wflow/gis/basins.geojson"))
gdf_sfincs = gpd.read_file(os.path.join(data_base, "sfincs/gis/region.geojson"))
gdf_sfincs = gdf_sfincs.to_crs("EPSG:4326")
gdf_snapwave = gpd.read_file(os.path.join(data_base, "snapwave/gis/SnapWave_region_sofala_only.shp"))
gdf_snapwave = gdf_snapwave.to_crs("EPSG:4326")
dfm_grid = gpd.read_file(os.path.join(data_base, "dfm/dfm_grid.gpkg"))


#%%
# Set up figure
fig, ax = plt.subplots(figsize=(9, 5))

# Plot SFINCS region and set up legend entry
gdf_sfincs.plot(ax=ax, edgecolor='pink', facecolor='pink', linewidth=1, alpha=0.5, zorder=4, rasterized=True)
sfincs_patch = mpatches.Patch(facecolor='pink', edgecolor='pink', alpha=0.5, label="SFINCS Region")

# Plot wflow basins region and set up legend entry
gdf_wflow.plot(ax=ax, edgecolor='lightskyblue', facecolor='lightskyblue', linewidth=1, alpha=0.5, zorder=3, rasterized=True)
wflow_patch = mpatches.Patch(facecolor='lightskyblue', edgecolor='lightskyblue', alpha=0.5, label="Wflow Basins")

# Plot snapwave region and set up legend entry
gdf_snapwave.plot(ax=ax, facecolor='#FFFF99', edgecolor='#FFFF99', linewidth=1, alpha=0.5, zorder=2, rasterized=True)
snap_patch = mpatches.Patch(facecolor='#FFFF99', edgecolor='#FFFF99', alpha=0.5, label="SnapWave Domain")

# Plotting of DFM Grid
dfm_grid.plot(ax=ax, edgecolor='white', linewidth=0.5, alpha=0.5, zorder=1, rasterized=True)

# Plot TC Tracks
tc_idai_filtered = tc_idai[tc_idai.geometry.y < -18]
tc_scatter = ax.scatter(tc_idai_filtered.geometry.x, tc_idai_filtered.geometry.y, s=tc_idai_filtered["size"]*0.5,  # Size based on wind speed
                        c=tc_idai_filtered["USA_WIND"], cmap=cmap_idai, alpha=0.7, label="Track TC Idai", zorder=6, rasterized=True)
ax.plot(tc_idai_filtered.geometry.x, tc_idai_filtered.geometry.y, color="grey", linewidth=1, alpha=0.7, linestyle="-", zorder=5, rasterized=True)
# Create a custom legend entry for the TC scatter plot (you can modify the color or markersize)
tc_marker = Line2D([0], [0], marker='o', color='darkgrey', markerfacecolor='red', markersize=5, label='Track TC Idai')

# Add basemap (LOWER zoom = faster)
ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=7, crs=tc_idai_filtered.crs, attribution=False, zorder=0, rasterized=True)

txt = ax.text(
    32.01, -21.99,  # x, y in figure coordinates (0=left/bottom, 1=right/top)
    "Tiles © Esri -- Source: Esri, i-cubed, USDA, USGS, AEX, \nGeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and \nthe GIS User Community",
    fontsize=5.5,
    color='white',
    alpha=0.7,
    ha='left',
    va='bottom',
    zorder=20,
)

# Define ticks for x and y axes
xticks = np.arange(32, 42, 1)    
yticks = np.arange(-21, -14, 1) 

ax.set_xticks(xticks)
ax.set_yticks(yticks)

ax.set_xticklabels([f"{x}°E" for x in xticks])
ax.set_yticklabels([f"{abs(y)}°S" for y in yticks])

# Set limits
extent = [32, -22, 39, -16.5]  # [lon_min, lat_min, lon_max, lat_max]
ax.set_xlim(extent[0], extent[2])
ax.set_ylim(extent[1], extent[3])

ax.text(33.85, -18.92, "SFINCS Domain", color='pink', fontsize=10, fontweight='bold',
    bbox=dict(boxstyle="round,pad=0.3", facecolor='grey', alpha=0.7, edgecolor='pink'),
    zorder=10)

ax.text(32.5, -17.3, "Wflow Domain", color='lightskyblue', fontsize=10, fontweight='bold',
        bbox=dict(boxstyle="round,pad=0.3", facecolor='grey', alpha=0.7, edgecolor='lightskyblue'),
        zorder=10)

ax.text(36.6, -21.7, "D-Flow FM Domain", color='white', fontsize=10, fontweight='bold',
        bbox=dict(boxstyle="round,pad=0.3", facecolor='grey', alpha=0.7, edgecolor='white'),
        zorder=10)

ax.text(35, -20.7, "SnapWave Domain", color='#FFFF99', fontsize=10, fontweight='bold',
        bbox=dict(boxstyle="round,pad=0.3", facecolor='grey', alpha=0.7, edgecolor='#FFFF99'),
        zorder=10)

ax.annotate("Track TC Idai", xy=(36.02,-19.8), xytext=(37.2,-19.15), textcoords='data',
            arrowprops=dict(arrowstyle="->", color='darkred', lw=1.5),
            bbox=dict(boxstyle="round,pad=0.3", fc="#d3d3d3", ec="darkred", alpha=0.7),
            fontsize=10, color='darkred', fontweight='bold', zorder=10)

# Add Legend
size_labels = ['0-50 km/h', '50-100 km/h', '>100 km/h']
legend_colors = ['#ff9999', '#ff4d4d', '#b30000']  # light to dark red
size_handles = [
    mlines.Line2D([], [], marker='o', color='k', markerfacecolor=legend_colors[0], alpha=0.7, markersize=2.5, label=size_labels[0], linestyle=''),
    mlines.Line2D([], [], marker='o', color='k', markerfacecolor=legend_colors[1], alpha=0.7, markersize=3.75, label=size_labels[1], linestyle=''),
    mlines.Line2D([], [], marker='o', color='k', markerfacecolor=legend_colors[2], alpha=0.7, markersize=5, label=size_labels[2], linestyle='')
]

# --- Create legend ---
leg = ax.legend(handles=size_handles, loc="upper right", fontsize=8, 
                edgecolor='grey', title="Wind speed", title_fontsize=8,
                bbox_to_anchor=(1, 1))

# --- Inset: East Africa context map placed outside
inset_position = [0.65, 0.11, 0.6, 0.6]  # [left, bottom, width, height]
inset_ax = fig.add_axes(inset_position, projection=ccrs.PlateCarree())
inset_ax.set_extent([20, 50.2, -35, 0], crs=ccrs.PlateCarree())
inset_ax.add_feature(cfeature.LAND, facecolor='lightgrey')
inset_ax.add_feature(cfeature.BORDERS, linewidth=0.4, color='grey')
inset_ax.add_feature(cfeature.COASTLINE, linewidth=0.3)

# --- Red bounding box on inset map for main map extent ---
rect = mpatches.Rectangle((extent[0], extent[1]),  # (x, y)
                          extent[2] - extent[0],   # width
                          extent[3] - extent[1],   # height
                          linewidth=1.5, edgecolor='red', facecolor='none', zorder=5)
inset_ax.add_patch(rect)

# --- Add country borders and names ---
# Read country shapefile from Natural Earth
shapefile = shpreader.natural_earth(resolution='110m',
                                     category='cultural',
                                     name='admin_0_countries')
gdf_countries = gpd.read_file(shapefile)

countries_to_label = [
    "Malawi", "Tanzania", "Botswana", "Eswatini", "Lesotho", "Madagascar", "Rwanda", "Burundi"]
gdf_filtered = gdf_countries[gdf_countries['NAME'].isin(countries_to_label)]

for _, row in gdf_filtered.iterrows():
    # Use representative point inside the country
    lon, lat = row.geometry.representative_point().coords[0]
    name = row['NAME']
    inset_ax.text(lon, lat, name, fontsize=7, transform=ccrs.PlateCarree(), 
                  ha='center', va='center', color='black', zorder=12)

# Add label for Mozambique with offset
inset_ax.text(32, -15.2, 'Mozambique', transform=ccrs.PlateCarree(),
              fontsize=8, weight='bold', ha='left', va='top', zorder=12)

# Correct other country label locations
inset_ax.text(25, -13.5, 'Zambia', transform=ccrs.PlateCarree(),
              fontsize=7, ha='left', va='top', zorder=12)
inset_ax.text(26, -18.2, 'Zimbabwe', transform=ccrs.PlateCarree(),
              fontsize=7, ha='left', va='top', zorder=12)
inset_ax.text(22.4, -27.2, 'South Africa', transform=ccrs.PlateCarree(),
              fontsize=7, ha='left', va='top', zorder=12)

inset_ax.set_title("East Africa", fontsize=8)

# Coordinates of the red box in inset (data coords)
x0, x1, y0, y1 = extent[0], extent[2], extent[1], extent[3]

# Top connection: from top-right of main → top-left of inset bbox
con_top = ConnectionPatch(
    xyA=(x0, y1), coordsA=inset_ax.transData,  # top-left in inset
    xyB=(x1, y1), coordsB=ax.transData,        # top-right in main
    color="red", lw=1, zorder=2)

# Bottom connection: from bottom-right of main → bottom-left of inset bbox
con_bottom = ConnectionPatch(
    xyA=(x0, y0), coordsA=inset_ax.transData,  # bottom-left in inset
    xyB=(x1, y0), coordsB=ax.transData,        # bottom-right in main
    color="red", lw=1, zorder=2)

# Add to the figure so it can span both axes
fig.add_artist(con_top)
fig.add_artist(con_bottom)

fig.savefig("../figures/f02.png", dpi=300, bbox_inches="tight", transparent=False)
fig.savefig("../figures/f02.pdf", dpi=300, bbox_inches="tight", transparent=False)
plt.show()
# %%

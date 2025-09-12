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
import dfm_tools as dfmt
import xarray as xr
from shapely.geometry import MultiPoint

#%% Load TC track shapefiles as geopandas geodataframe
shapefile_path = "p:/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS/IBTrACS_IDAI.shp"
tc_idai = gpd.read_file(shapefile_path)

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
gdf_wflow = gpd.read_file(r"p:\11210471-001-compass\02_Models\sofala\Idai\wflow\staticgeoms\basins.geojson")
gdf_sfincs = gpd.read_file(r"p:\11210471-001-compass\02_Models\sofala\Idai\sfincs\gis\region.geojson")
gdf_sfincs = gdf_sfincs.to_crs("EPSG:4326")
gdf_snapwave = gpd.read_file(r"P:\11210471-001-compass\01_Data\sofala_geoms\SnapWave_region_sofala_only.shp")
gdf_snapwave = gdf_snapwave.to_crs("EPSG:4326")

# Getting DFM grid
dir_runs = f'p:/11210471-001-compass/03_Runs/sofala/Idai/dfm'
F__model = f'event_450_gebco2024_MZB_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0'

file_nc_map_F = []
for fname in os.listdir(os.path.join(dir_runs,F__model,'output')):
    if fname.endswith("map.nc"):
        print(fname)
        file_nc_map_F.append(os.path.join(dir_runs,F__model,'output',fname))
ds_map_F = dfmt.open_partitioned_dataset(file_nc_map_F)


#%%
# Calculate surface area of DFM grid
# Extract node coordinates
x = ds_map_F['mesh2d_node_x'].values
y = ds_map_F['mesh2d_node_y'].values

# Create a convex hull polygon of all nodes (outer edge)
points = MultiPoint(list(zip(x, y)))
boundary_polygon = points.convex_hull  # outer edge

# Save to GeoDataFrame
gdf = gpd.GeoDataFrame(index=[0], crs="EPSG:4326", geometry=[boundary_polygon])

# Remove land from DFM grid boundary
land_gdf = gpd.read_file(r"p:\11210471-001-compass\01_Data\land_polygon\ne_10m_land\ne_10m_land.shp").to_crs(gdf.crs)
ocean_only = gpd.overlay(gdf, land_gdf, how='difference')

# Your existing DFM grid plot
fig, ax = plt.subplots(figsize=(8, 8))
ds_map_F.grid.plot(ax=ax, edgecolor='white', linewidth=0.5, alpha=0.5, zorder=1)
ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=7, crs=gdf.crs, attribution=False, zorder=0)

# Plot the boundary polygon
gdf.boundary.plot(ax=ax, edgecolor='red', linewidth=2, zorder=2)  # adjust color/width
ocean_only.plot(ax=ax, edgecolor='orange', alpha=0.5, linewidth=2, zorder=2)  # adjust color/width
ax.set_title("DFM Grid with Outer Boundary")
plt.show()

# Reproject to UTM 36S
ocean_utm = ocean_only.to_crs("EPSG:32736")
# Calculate area in m²
ocean_utm["area_m2"] = ocean_utm.geometry.area
# Sum total ocean area
total_ocean_area = ocean_utm["area_m2"].sum()
total_ocean_area_rounded = round((total_ocean_area/1e6), -2)
print(f"Total ocean area: {total_ocean_area_rounded} km²")


#%%
# Set up figure
fig, ax = plt.subplots(figsize=(9, 5))

# Plot SFINCS region and set up legend entry
gdf_sfincs.plot(ax=ax, edgecolor='pink', facecolor='pink', linewidth=1, alpha=0.5, zorder=4, rasterized=True)
sfincs_patch = mpatches.Patch(facecolor='pink', edgecolor='pink', alpha=0.5, label="SFINCS Region")

# Plot wflow basins region and set up legend entry
gdf_wflow.plot(ax=ax, edgecolor='lightskyblue', facecolor='lightskyblue', linewidth=1, alpha=0.5, zorder=3, rasterized=True)
wflow_patch = mpatches.Patch(facecolor='lightskyblue', edgecolor='lightskyblue', alpha=0.5, label="Wflow Basins")

# Plot wflow basins region and set up legend entry
gdf_snapwave.plot(ax=ax, facecolor='#FFFF99', edgecolor='#FFFF99', linewidth=1, alpha=0.5, zorder=2, rasterized=True)
snap_patch = mpatches.Patch(facecolor='#FFFF99', edgecolor='#FFFF99', alpha=0.5, label="SnapWave Domain")

# Plotting of DFM Grid
ds_map_F.grid.plot(ax=ax, edgecolor='white', linewidth=0.5, alpha=0.5, zorder=1, rasterized=True)

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
inset_position = [0.63, 0.11, 0.6, 0.6]  # [left, bottom, width, height]
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

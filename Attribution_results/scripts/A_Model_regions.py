#%% Import the necessary packages
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx
import numpy as np
from matplotlib import cm
import matplotlib.lines as mlines
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import cartopy.feature as cfeature
from shapely.geometry import box
from cartopy.feature import ShapelyFeature
import cartopy.io.shapereader as shpreader

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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Getting the model regions
gdf_wflow = gpd.read_file(r"p:\11210471-001-compass\02_Models\sofala\Idai\wflow\staticgeoms\basins.geojson")
gdf_sfincs = gpd.read_file(r"p:\11210471-001-compass\02_Models\sofala\Idai\sfincs\gis\region.geojson")
gdf_sfincs = gdf_sfincs.to_crs("EPSG:4326")
grid_dfm = xr.open_dataset(r"p:\11210471-001-compass\02_Models\sofala\Idai\dfm\base_450_gebco2024_MZB_GTSMv41opendap\grid_network.nc")

# Set up figure
fig, ax = plt.subplots(figsize=(9, 5))

# Plot SFINCS region and set up legend entry
gdf_sfincs.plot(ax=ax, edgecolor='pink', facecolor='pink', linewidth=1, alpha=0.5, zorder=3)
sfincs_patch = mpatches.Patch(facecolor='pink', edgecolor='pink', alpha=0.5, label="SFINCS Region")

# Plot wflow basins region and set up legend entry
gdf_wflow.plot(ax=ax, edgecolor='lightskyblue', facecolor='lightskyblue', linewidth=1, alpha=0.5, zorder=2)
wflow_patch = mpatches.Patch(facecolor='lightskyblue', edgecolor='lightskyblue', alpha=0.5, label="Wflow Basins")

# Efficient Plotting of DFM Grid (reduce marker size or plot subset of points)
dfm_scatter = ax.scatter(grid_dfm['mesh2d_node_x'], grid_dfm['mesh2d_node_y'], 
                         s=1, color="white", alpha=0.7, label="DFM Grid", zorder=1)

# Plot TC Tracks
tc_idai_filtered = tc_idai[tc_idai.geometry.y < -18]
tc_scatter = ax.scatter(tc_idai_filtered.geometry.x, tc_idai_filtered.geometry.y, s=tc_idai_filtered["size"]*0.5,  # Size based on wind speed
                        c=tc_idai_filtered["USA_WIND"], cmap=cmap_idai, alpha=0.7, label="Track TC Idai", zorder=6)
ax.plot(tc_idai_filtered.geometry.x, tc_idai_filtered.geometry.y, color="grey", linewidth=1, alpha=0.7, linestyle="-", zorder=5)
# Create a custom legend entry for the TC scatter plot (you can modify the color or markersize)
tc_marker = Line2D([0], [0], marker='o', color='darkgrey', markerfacecolor='red', markersize=5, label='Track TC Idai')

# Add basemap (LOWER zoom = faster)
ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=7, crs=tc_idai_filtered.crs, attribution=False, zorder=0)

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

ax.text(36.3, -21, "DFM Domain", color='white', fontsize=10, fontweight='bold',
        bbox=dict(boxstyle="round,pad=0.3", facecolor='grey', alpha=0.7, edgecolor='white'),
        zorder=10)

ax.annotate("Track TC Idai", xy=(36.5,-19.8), xytext=(37.2,-19.15), textcoords='data',
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
ax.legend(handles=size_handles, loc="upper right", fontsize=8, 
          edgecolor='grey', title="Wind Speed", title_fontsize=8, bbox_to_anchor=(1.27, 1))

ax.set_title("Model Domains and Track TC Idai", fontsize=11)

# --- Inset: East Africa context map placed outside
inset_position = [0.65, 0.1, 0.6, 0.6]  # [left, bottom, width, height]
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
    "Zimbabwe", "Malawi", "Tanzania", "Botswana", "Eswatini", "Lesotho", "Madagascar", "Rwanda", "Burundi"]
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
inset_ax.text(22, -27.2, 'South Africa', transform=ccrs.PlateCarree(),
              fontsize=7, ha='left', va='top', zorder=12)

# Optional: remove axis ticks and labels
# inset_ax.set_xticks([])
# inset_ax.set_yticks([])
inset_ax.set_title("East Africa", fontsize=8)

# Define corners of the extent in (lon, lat)
# inset_pts = [
#     (extent[0], extent[1]),  # bottom-left
#     (extent[0], extent[3]),  # top-left
# ]

# # Connect to the corresponding right corners on main map
# main_pts = [
#     (extent[2], extent[1]),  # bottom-right
#     (extent[2], extent[3]),  # top-right
# ]

# # Coordinate transforms
# # 1. Transform from data to display
# display_coords_inset = [inset_ax.projection.transform_point(x, y, src_crs=ccrs.PlateCarree()) for x, y in inset_pts]
# display_coords_main = main_pts  # plain lon/lat, already PlateCarree

# # 2. Convert display coords to pixel coordinates
# trans_inset_data = inset_ax.transData.transform
# trans_main_data = ax.transData.transform
# inv_fig = fig.transFigure.inverted()

# for (ix, iy), (mx, my) in zip(display_coords_inset, display_coords_main):
#     # Convert data to pixel space, then to figure space
#     p1_fig = inv_fig.transform(trans_inset_data((ix, iy)))
#     p2_fig = inv_fig.transform(trans_main_data((mx, my)))

#     line = mlines.Line2D([p1_fig[0], p2_fig[0]], [p1_fig[1], p2_fig[1]],
#                          transform=fig.transFigure,
#                          color='red', linewidth=1)
#     fig.lines.append(line)


plt.savefig("../figures/model_domains_sofala.png", dpi=300, bbox_inches="tight", transparent=False)
plt.show()
# %%

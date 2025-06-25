#%% Import the necessary packages
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx
import numpy as np
from matplotlib import cm
import matplotlib.lines as mlines
import xarray as xr
import cartopy.crs as ccrs
# import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

#%% Load TC track shapefiles as geopandas geodataframe
shapefile_path = "p:/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS/IBTrACS_IDAI.shp"
tc_idai = gpd.read_file(shapefile_path)

shapefile_path = "p:/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS/IBTrACS_KENNETH.shp"
tc_kenneth = gpd.read_file(shapefile_path)

shapefile_path = "p:/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS/IBTrACS_FREDDY_part1.shp"
tc_freddy1 = gpd.read_file(shapefile_path)

shapefile_path = "p:/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS/IBTrACS_FREDDY_part2.shp"
tc_freddy2 = gpd.read_file(shapefile_path)

# %%
# Normalize windspeed for point sizing (adjust scaling as needed)
tc_idai["size"] = (tc_idai["USA_WIND"] - tc_idai["USA_WIND"].min()) / (
    tc_idai["USA_WIND"].max() - tc_idai["USA_WIND"].min()
) * 138 + 19  # Scale between 20 and 145

# Normalize windspeed for tc_kenneth wind sizing
tc_kenneth["size"] = (tc_kenneth["USA_WIND"] - tc_kenneth["USA_WIND"].min()) / (
    tc_kenneth["USA_WIND"].max() - tc_kenneth["USA_WIND"].min()
) * 138 + 19  # Scale between 20 and 145

tc_freddy1["size"] = (tc_freddy1["USA_WIND"] - tc_freddy1["USA_WIND"].min()) / (
    tc_freddy1["USA_WIND"].max() - tc_freddy1["USA_WIND"].min()
) * 138 + 19  # Scale between 20 and 145

# Normalize windspeed for tc_kenneth wind sizing
tc_freddy2["size"] = (tc_freddy2["USA_WIND"] - tc_freddy2["USA_WIND"].min()) / (
    tc_freddy2["USA_WIND"].max() - tc_freddy2["USA_WIND"].min()
) * 138 + 19  # Scale between 20 and 145

#%%
# Define bounding box for Mozambique region
bbox = [31, -27, 50, -10]  # [min_lon, min_lat, max_lon, max_lat]

colors = {"TC Idai": "darkred", "TC Kenneth": "darkorange", "TC Freddy - 1": "dodgerblue", "TC Freddy - 2": "dodgerblue"}

# Maintain correct aspect ratio
fig_width = 8  # Fixed width
aspect_ratio = abs((bbox[3] - bbox[1]) / (bbox[2] - bbox[0]))  # lat range / lon range
fig_height = fig_width * aspect_ratio  # Adjust height

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Define colors and colormap
cmap_idai = cm.get_cmap("Reds")
cmap_kenneth = cm.get_cmap("Oranges")
cmap_freddy1 = cm.get_cmap("Blues")
cmap_freddy2 = cm.get_cmap("Blues")

# Normalize windspeed for color mapping
norm_idai = plt.Normalize(vmin=tc_idai["USA_WIND"].min(), vmax=tc_idai["USA_WIND"].max())
norm_kenneth = plt.Normalize(vmin=tc_kenneth["USA_WIND"].min(), vmax=tc_kenneth["USA_WIND"].max())
norm_freddy1 = plt.Normalize(vmin=tc_freddy1["USA_WIND"].min(), vmax=tc_freddy1["USA_WIND"].max())
norm_freddy2 = plt.Normalize(vmin=tc_freddy2["USA_WIND"].min(), vmax=tc_freddy2["USA_WIND"].max())

# Plot grey track lines for each cyclone
for tc in [tc_idai, tc_kenneth, tc_freddy1, tc_freddy2]:
    ax.plot(
        tc.geometry.x, 
        tc.geometry.y, 
        color="grey",  
        linewidth=1,  # Adjust line thickness
        alpha=0.7,  # Make it slightly transparent
        linestyle="-",
        zorder=1
    )

# Plot TC tracks with different colors and sizes
for tc, name, cmap, norm in zip([tc_idai, tc_kenneth, tc_freddy1, tc_freddy2], colors.keys(), [cmap_idai, cmap_kenneth, cmap_freddy1, cmap_freddy2], [norm_idai, norm_kenneth, norm_freddy1, norm_freddy2]):
    sc = ax.scatter(
        tc.geometry.x, 
        tc.geometry.y, 
        s=tc["size"],  # Size based on wind speed
        c=tc["USA_WIND"],  # Color based on wind speed
        cmap=cmap, 
        norm=norm,
        label=name,  # Add to legend
        edgecolor="black",
        alpha=0.7,
        zorder=2
    )

# Add Google Maps-style basemap (add zoom=10, )
ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=7, crs=tc_idai.crs, attribution=False)

# Load the world shapefile (update the path to the downloaded shapefile)
world_shapefile = "C:/Code/processing/naturalearthdata/ne_50m_admin_0_countries.shp"
world = gpd.read_file(world_shapefile)

# Extract Mozambique border
world.boundary.plot(ax=ax, linewidth=1, color="white", alpha=0.5)

# Set lat/lon limits
ax.set_xlim([bbox[0], bbox[2]])
ax.set_ylim([bbox[1], bbox[3]])

# Axis labels and ticks
ax.set_xlabel("Longitude", fontsize=10)
ax.set_ylabel("Latitude", fontsize=10)
ax.tick_params(axis="both", which="major", labelsize=10)

# Grid
ax.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)

# Add legend
# Manually create legend for TC colors
handles = [
    mlines.Line2D([], [], color="darkred", label="TC Idai"),
    mlines.Line2D([], [], color="darkorange", label="TC Kenneth"),
    mlines.Line2D([], [], color="dodgerblue", label="TC Freddy"),
]

# Manually create legend for windspeed categories (size categories)
size_labels = ['0-50 km/h', '50-100 km/h', '>100 km/h']
size_handles = [
    mlines.Line2D([], [], marker='o', color='w', markerfacecolor='white', markeredgecolor='black', markersize=5, label=size_labels[0]),
    mlines.Line2D([], [], marker='o', color='w', markerfacecolor='white', markeredgecolor='black',markersize=7.5, label=size_labels[1]),
    mlines.Line2D([], [], marker='o', color='w', markerfacecolor='white', markeredgecolor='black',markersize=10, label=size_labels[2])
]

# Add both legend entries to the plot
ax.legend(handles=handles + size_handles, fontsize=10, loc="upper left")

# Save figure with transparent background (optional)
# plt.savefig("mozambique_tc_tracks.png", dpi=300, bbox_inches="tight", transparent=True)

plt.show()
# %% PLOTTING MODEL DOMAIN FIGURES FOR PAPER
# Getting the model regions
gdf_wflow = gpd.read_file(r"p:\11210471-001-compass\02_Models\sofala\Idai\wflow\staticgeoms\basins.geojson")
gdf_sfincs = gpd.read_file(r"p:\11210471-001-compass\02_Models\sofala\Idai\sfincs\gis\region.geojson")
gdf_sfincs = gdf_sfincs.to_crs("EPSG:4326")
grid_dfm = xr.open_dataset(r"p:\11210471-001-compass\02_Models\sofala\Idai\dfm\base_450_gebco2024_MZB_GTSMv41opendap\grid_network.nc")

# Set up figure
fig, ax = plt.subplots(figsize=(8, 6))
# fig.set_size_inches(10, 8, forward=True)

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
tc_scatter = ax.scatter(tc_idai.geometry.x, tc_idai.geometry.y, s=tc_idai["size"]*0.5,  # Size based on wind speed
                        c=tc_idai["USA_WIND"], cmap=cmap_idai, alpha=0.7, label="Track TC Idai", zorder=6)
ax.plot(tc_idai.geometry.x, tc_idai.geometry.y, color="grey", linewidth=1, alpha=0.7, linestyle="-", zorder=5)
# Create a custom legend entry for the TC scatter plot (you can modify the color or markersize)
tc_marker = Line2D([0], [0], marker='o', color='darkgrey', markerfacecolor='red', markersize=5, label='Track TC Idai')

# Add basemap (LOWER zoom = faster)
ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=7, crs=tc_idai.crs, attribution=False, zorder=0)

# Set limits
ax.set_xlim([32, 41])
ax.set_ylim([-21.2, -15])

# Add Legend
size_handles = [
    mlines.Line2D([], [], marker='o', color='w', markerfacecolor='white', alpha=0.5, markersize=2.5, label=size_labels[0], linestyle=''),
    mlines.Line2D([], [], marker='o', color='w', markerfacecolor='white', alpha=0.5, markersize=3.75, label=size_labels[1], linestyle=''),
    mlines.Line2D([], [], marker='o', color='w', markerfacecolor='white', alpha=0.5, markersize=5, label=size_labels[2], linestyle='')
]

# Combine all legend handles into one list
combined_handles = size_handles + [sfincs_patch, wflow_patch, dfm_scatter, tc_marker]

# Add the combined legend with the title for wind speed size
ax.legend(handles=combined_handles, loc="upper left", fontsize=8, facecolor='grey', edgecolor='white', title="Wind Speed", title_fontsize=8)

ax.set_xlabel("Longitude", fontsize=10)
ax.set_ylabel("Latitude", fontsize=10)


# plt.savefig(r"p:\11210471-001-compass\04_Results\model_regions\sofala_Idai_dfm_legends.png", dpi=300, bbox_inches="tight", transparent=False)

# Show plot
plt.show()
# %%

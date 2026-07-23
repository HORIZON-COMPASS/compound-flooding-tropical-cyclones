# %% Import the necessary packages
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx
import numpy as np
from matplotlib import cm
import matplotlib.lines as mlines
import xarray as xr
import cartopy.crs as ccrs
import seaborn as sns

# %% Load TC track shapefiles as geopandas geodataframe
shapefile_path = (
    "p:/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS/IBTrACS_IDAI.shp"
)
tc_idai = gpd.read_file(shapefile_path)

shapefile_path = (
    "p:/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS/IBTrACS_KENNETH.shp"
)
tc_kenneth = gpd.read_file(shapefile_path)

shapefile_path = (
    "p:/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS/IBTrACS_FREDDY_part1.shp"
)
tc_freddy1 = gpd.read_file(shapefile_path)

shapefile_path = (
    "p:/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS/IBTrACS_FREDDY_part2.shp"
)
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

# %%
# Define bounding box for Mozambique region
bbox = [31, -27, 50, -10]  # [min_lon, min_lat, max_lon, max_lat]

colors = {
    "TC Idai": "darkred",
    "TC Kenneth": "darkorange",
    "TC Freddy - 1": "dodgerblue",
    "TC Freddy - 2": "dodgerblue",
}

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
norm_idai = plt.Normalize(
    vmin=tc_idai["USA_WIND"].min(), vmax=tc_idai["USA_WIND"].max()
)
norm_kenneth = plt.Normalize(
    vmin=tc_kenneth["USA_WIND"].min(), vmax=tc_kenneth["USA_WIND"].max()
)
norm_freddy1 = plt.Normalize(
    vmin=tc_freddy1["USA_WIND"].min(), vmax=tc_freddy1["USA_WIND"].max()
)
norm_freddy2 = plt.Normalize(
    vmin=tc_freddy2["USA_WIND"].min(), vmax=tc_freddy2["USA_WIND"].max()
)

# Plot grey track lines for each cyclone
for tc in [tc_idai, tc_kenneth, tc_freddy1, tc_freddy2]:
    ax.plot(
        tc.geometry.x,
        tc.geometry.y,
        color="grey",
        linewidth=1,  # Adjust line thickness
        alpha=0.7,  # Make it slightly transparent
        linestyle="-",
        zorder=1,
    )

# Plot TC tracks with different colors and sizes
for tc, name, cmap, norm in zip(
    [tc_idai, tc_kenneth, tc_freddy1, tc_freddy2],
    colors.keys(),
    [cmap_idai, cmap_kenneth, cmap_freddy1, cmap_freddy2],
    [norm_idai, norm_kenneth, norm_freddy1, norm_freddy2],
):
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
        zorder=2,
    )

# Add Google Maps-style basemap (add zoom=10, )
ctx.add_basemap(
    ax,
    source=ctx.providers.Esri.WorldImagery,
    zoom=7,
    crs=tc_idai.crs,
    attribution=False,
)

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
size_labels = ["0-50 km/h", "50-100 km/h", ">100 km/h"]
size_handles = [
    mlines.Line2D(
        [],
        [],
        marker="o",
        color="w",
        markerfacecolor="white",
        markeredgecolor="black",
        markersize=5,
        label=size_labels[0],
    ),
    mlines.Line2D(
        [],
        [],
        marker="o",
        color="w",
        markerfacecolor="white",
        markeredgecolor="black",
        markersize=7.5,
        label=size_labels[1],
    ),
    mlines.Line2D(
        [],
        [],
        marker="o",
        color="w",
        markerfacecolor="white",
        markeredgecolor="black",
        markersize=10,
        label=size_labels[2],
    ),
]

# Add both legend entries to the plot
ax.legend(handles=handles + size_handles, fontsize=10, loc="upper left")

# Save figure with transparent background (optional)
# plt.savefig("mozambique_tc_tracks.png", dpi=300, bbox_inches="tight", transparent=True)

plt.show()

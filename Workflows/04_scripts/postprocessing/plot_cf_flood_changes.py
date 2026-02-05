import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import BoundaryNorm, ListedColormap
import rioxarray as rxr  # Required for reading TIFF files
import contextily as ctx  # For adding basemap tiles
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

# City markers for different regions
try:
    from mozambique_cities import mozambique_cities
except ImportError:
    mozambique_cities = []

# South Africa cities (Durban area)
south_africa_cities = [
    {"name": "Durban", "lat": -29.8587, "lon": 31.0218},
    # {"name": "Pinetown", "lat": -29.8167, "lon": 30.8667},
    # {"name": "Umlazi", "lat": -29.9667, "lon": 30.8833},
    # {"name": "Amanzimtoti", "lat": -30.0500, "lon": 30.8833},
]

# Map events to their city lists
EVENT_CITIES = {
    "Freddy": mozambique_cities,
    "Kenneth": mozambique_cities,
    "Idai": mozambique_cities,
    "Durban_April2022": south_africa_cities,
}


# ===== HELPER FUNCTIONS =====
def add_city_markers(ax, extent, cities, fontsize=8, markersize=50):
    """
    Add city markers to a map if they fall within the plot extent.

    Parameters
    ----------
    ax : matplotlib axis
        The axis to plot on
    extent : tuple
        (minx, maxx, miny, maxy) extent of the plot in EPSG:4326
    cities : list
        List of city dictionaries with 'name', 'lat', 'lon' keys
    fontsize : int
        Font size for city labels
    markersize : int
        Size of city markers
    """
    minx, maxx, miny, maxy = extent

    cities_in_extent = []
    for city in cities:
        if minx <= city["lon"] <= maxx and miny <= city["lat"] <= maxy:
            cities_in_extent.append(city)

    if not cities_in_extent:
        return 0

    # Plot city markers
    for city in cities_in_extent:
        ax.scatter(
            city["lon"],
            city["lat"],
            s=markersize,
            c="red",
            marker="o",
            edgecolors="white",
            linewidths=1.5,
            zorder=10,
            alpha=0.8,
        )

        # Add city label with background for better visibility
        # Offset label slightly up and to the right of the marker
        label_offset_x = 0.02  # Degrees longitude
        label_offset_y = 0.02  # Degrees latitude
        ax.text(
            city["lon"] + label_offset_x,
            city["lat"] + label_offset_y,
            city["name"],
            fontsize=fontsize,
            ha="left",
            va="bottom",
            zorder=11,
            bbox=dict(
                boxstyle="round,pad=0.3", facecolor="white", alpha=0.7, edgecolor="none"
            ),
        )

    return len(cities_in_extent)


# File paths
# TODO: this is a preliminary script that is still not fully embedded in the snakemake workflow. The user is required to change the events & input files for now.
# For idai runs, I tested: CF: event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10; F: event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0

# ===== CONFIGURATION =====
# Set your event name here
EVENT_NAME = "Durban_April2022"  # Change this to: "Kenneth", "Freddy", "Idai", "Durban_April2022"

# Base paths - update these as needed
OUTPUT_DIR = Path("/p/11210471-001-compass/04_Results/CF_figs")

# ===== EVENT-SPECIFIC CONFIGURATION =====
# Maps event names to their specific folder paths and base directories
EVENT_CONFIG = {
    "Freddy": {
        "base_path": Path("/p/11210471-001-compass/03_Runs/test"),
        "factual": "event_tp_era5_hourly_CF0_GTSMv41opendap_CF0_no_wind_CF0",
        "counterfactual": "event_tp_era5_hourly_CF-8_GTSMv41opendap_CF0_no_wind_CF0",
        "folder_name": "Freddy",  # Folder name in base_path
    },
    "Kenneth": {
        "base_path": Path("/p/11210471-001-compass/03_Runs/test"),
        "factual": "event_tp_era5_hourly_zarr_CF0_GTSMv41opendap_CF0_no_wind_CF0",
        "counterfactual": "event_tp_era5_hourly_zarr_CF-8_GTSMv41opendap_CF0_no_wind_CF0",
        "folder_name": "Kenneth",  # Folder name in base_path
    },
    "Idai": {
        "base_path": Path("/p/11210471-001-compass/03_Runs/sofala"),
        "factual": "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0",
        "counterfactual": "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10",
        "folder_name": "Idai",  # Folder name in base_path
    },
    "Durban_April2022": {
        "base_path": Path("/p/11210471-001-compass/03_Runs/durban"),
        "factual": "event_precip_era5_hourly_CF0_no_wind",
        "counterfactual": "event_precip_era5_hourly_CF-8_no_wind",
        "folder_name": "Durban2022",  # Folder name in base_path
    },
}

# Validate event name
if EVENT_NAME not in EVENT_CONFIG:
    raise ValueError(
        f"Unknown event: {EVENT_NAME}. Valid options: {list(EVENT_CONFIG.keys())}"
    )

# Get event-specific configuration
event_cfg = EVENT_CONFIG[EVENT_NAME]
BASE_RUN_PATH = event_cfg["base_path"]
FOLDER_NAME = event_cfg.get(
    "folder_name", EVENT_NAME
)  # Use folder_name if provided, else EVENT_NAME
EVENT_CITY_LIST = EVENT_CITIES.get(EVENT_NAME, [])  # Get cities for this event

# ===== DYNAMIC FILE PATHS =====
# Construct file paths based on event name
file_cf0 = (
    BASE_RUN_PATH
    / FOLDER_NAME
    / "sfincs"
    / event_cfg["factual"]
    / "plot_output"
    / "sfincs_output_hmax_AllTime.tif"
)
file_cf8 = (
    BASE_RUN_PATH
    / FOLDER_NAME
    / "sfincs"
    / event_cfg["counterfactual"]
    / "plot_output"
    / "sfincs_output_hmax_AllTime.tif"
)

# Create output directory if it doesn't exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Processing event: {EVENT_NAME}")

# ===== LOAD PREPROCESSED TIFF FILES =====
print("Loading preprocessed flood depth files...")
zsmax_cf0 = rxr.open_rasterio(file_cf0)
zsmax_cf8 = rxr.open_rasterio(file_cf8)

# Remove band dimension if present (common with TIFF files)
if "band" in zsmax_cf0.dims:
    zsmax_cf0 = zsmax_cf0.squeeze("band", drop=True)
if "band" in zsmax_cf8.dims:
    zsmax_cf8 = zsmax_cf8.squeeze("band", drop=True)

# Store original CRS for reference
original_crs = zsmax_cf0.rio.crs
print(f"Original CRS: {original_crs}")

# ===== CALCULATE METRICS IN ORIGINAL CRS (meters) =====
# This ensures accurate area/volume calculations before any reprojection
print("Calculating flood metrics in original projected CRS...")

# Handle NaN values for difference calculation (in original CRS)
mask_cf0_valid = ~np.isnan(zsmax_cf0)
mask_cf8_valid = ~np.isnan(zsmax_cf8)

# For CF0: where CF0 is NaN but CF8 has a value, set CF0 to 0
zsmax_cf0_calc = zsmax_cf0.where(mask_cf0_valid | ~mask_cf8_valid, 0)
# For CF8: where CF8 is NaN but CF0 has a value, set CF8 to 0
zsmax_cf8_calc = zsmax_cf8.where(mask_cf8_valid | ~mask_cf0_valid, 0)

# Calculate difference in original CRS
diff_calc = zsmax_cf0_calc - zsmax_cf8_calc

# Flood extent (area) change - calculated in original projected CRS (meters)
flood_threshold = 0.05  # meters
flood_mask_cf0 = zsmax_cf0_calc > flood_threshold
flood_mask_cf8 = zsmax_cf8_calc > flood_threshold

# Get cell size in meters from the original projected CRS
try:
    dx = float(np.abs(zsmax_cf0.x.diff("x").median()))
    dy = float(np.abs(zsmax_cf0.y.diff("y").median()))
    cell_area_m2 = dx * dy
    print(f"Cell size: {dx:.2f} x {dy:.2f} m = {cell_area_m2:.2f} m²")
except Exception:
    cell_area_m2 = np.nan
    print("Warning: Could not determine cell size")

area_cf0_m2 = flood_mask_cf0.sum().values * cell_area_m2
area_cf8_m2 = flood_mask_cf8.sum().values * cell_area_m2
extent_area_diff_m2 = area_cf0_m2 - area_cf8_m2
extent_area_pct = (
    (extent_area_diff_m2 / area_cf8_m2 * 100.0) if area_cf8_m2 > 0 else np.nan
)

# Flood volume (sum depth * cell area over flooded cells) - in original CRS
if not np.isnan(cell_area_m2):
    volume_cf0_m3 = float(
        zsmax_cf0_calc.where(flood_mask_cf0).sum(skipna=True) * cell_area_m2
    )
    volume_cf8_m3 = float(
        zsmax_cf8_calc.where(flood_mask_cf8).sum(skipna=True) * cell_area_m2
    )
    volume_diff_m3 = volume_cf0_m3 - volume_cf8_m3
    volume_pct = (
        (volume_diff_m3 / volume_cf8_m3 * 100.0) if volume_cf8_m3 > 0 else np.nan
    )
else:
    volume_cf0_m3 = volume_cf8_m3 = volume_diff_m3 = volume_pct = np.nan

# Mean differences restricted to flooded cells (using same threshold)
flood_mask_union = (zsmax_cf0_calc > flood_threshold) | (
    zsmax_cf8_calc > flood_threshold
)
flood_mask_intersection = (zsmax_cf0_calc > flood_threshold) & (
    zsmax_cf8_calc > flood_threshold
)
mean_diff_flooded_union = diff_calc.where(flood_mask_union).mean(skipna=True)

# Percent difference (relative change) - calculated before reprojection
percent_diff = xr.where(
    zsmax_cf8_calc > 0, (diff_calc / zsmax_cf8_calc) * 100.0, np.nan
)

# ===== REPROJECT TO LAT/LON FOR PLOTTING =====
# Only reproject AFTER metrics are calculated
if original_crs != "EPSG:4326":
    print("Reprojecting to EPSG:4326 (lat/lon) for plotting...")
    zsmax_cf0 = zsmax_cf0_calc.rio.reproject("EPSG:4326")
    zsmax_cf8 = zsmax_cf8_calc.rio.reproject("EPSG:4326")
    diff = diff_calc.rio.reproject("EPSG:4326")
    print("Reprojection complete")
else:
    print("Already in EPSG:4326")
    zsmax_cf0 = zsmax_cf0_calc
    zsmax_cf8 = zsmax_cf8_calc
    diff = diff_calc

# Only show areas with positive water depth (above ground) - for plotting
zsmax_cf0_plot = zsmax_cf0.where(zsmax_cf0 > 0.05)
zsmax_cf8_plot = zsmax_cf8.where(zsmax_cf8 > 0.05)

# Summary metrics
max_depth_cf0 = float(zsmax_cf0.max())
max_depth_cf8 = float(zsmax_cf8.max())
mean_diff_all = float(diff.mean())
mean_diff_union = float(mean_diff_flooded_union)

# Final concise summary
print("\n================ SUMMARY ================")
print(f"Event: {EVENT_NAME}")

print(f"\n[1] Flood Extent (threshold > {flood_threshold} m)")
print(f"  Factual area:        {area_cf0_m2/1e6:.3f} km^2")
print(f"  Counterfactual area: {area_cf8_m2/1e6:.3f} km^2")
print(f"  Absolute change:     {extent_area_diff_m2/1e6:.3f} km^2")
print(f"  Percent change:      {extent_area_pct:.2f} %")

print("\n[2] Flood Depth")
print(f"  Max depth factual:        {max_depth_cf0:.3f} m")
print(f"  Max depth counterfactual: {max_depth_cf8:.3f} m")
print(f"  Mean depth diff (all cells):        {mean_diff_all:.3f} m")
print(f"  Mean depth diff (flooded union):    {mean_diff_union:.3f} m")
print(
    f"  Mean percent depth change (where CF-8 > 0): {percent_diff.mean(skipna=True).values:.2f} %"
)

print(f"\n[3] Flood Volume (threshold > {flood_threshold} m)")
print(f"  Factual volume:        {volume_cf0_m3/1e6:.3f} Mm^3")
print(f"  Counterfactual volume: {volume_cf8_m3/1e6:.3f} Mm^3")
print(f"  Absolute change:       {volume_diff_m3/1e6:.3f} Mm^3")
print(f"  Percent change:        {volume_pct:.2f} %")
print("=========================================\n")

# Remove previous detailed exceedance and cell count prints
# ===== DETERMINE COORDINATE SYSTEM =====
# Try to determine UTM zone from the data coordinates
try:
    # Get approximate center coordinates
    center_x = float(zsmax_cf0.x.mean())
    center_y = float(zsmax_cf0.y.mean())

    # Guess UTM zone based on coordinates (this is approximate)
    if center_x > 300000 and center_x < 800000:  # Typical UTM coordinate range
        if center_y < 0:  # Southern hemisphere
            utm_zone = 37  # Adjust based on your region
            southern = True
        else:
            utm_zone = 37  # Adjust based on your region
            southern = False
    else:
        utm_zone = 37  # Default for your region
        southern = True

    print(f"Using UTM zone {utm_zone}, Southern: {southern}")
except Exception as e:
    print(f"Could not determine UTM zone: {e}")
    utm_zone = 37
    southern = True

# ===== CREATE MAIN COMPARISON PLOT =====
print("Creating plots...")

# Data is now in EPSG:4326
data_crs = "EPSG:4326"
print(f"Data CRS for plotting: {data_crs}")

# Use regular matplotlib axes (not Cartopy) for contextily compatibility
fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# Define colormaps and levels
flood_levels = np.arange(0, 4, 0.2)
flood_cmap = plt.cm.viridis

# Difference levels (starting from 0 for Reds colormap)
dif_min = 0
dif_max = 0.5
diff_levels = np.linspace(dif_min, dif_max, 11)
diff_cmap = plt.cm.Reds  # Reds colormap with 0 as minimum

# Plot CF0
ax1 = axes[0]
im1 = zsmax_cf0_plot.plot(
    ax=ax1,
    levels=flood_levels,
    cmap=flood_cmap,
    add_colorbar=False,
    x="x",
    y="y",
    alpha=0.7,
    zorder=2,
)
ax1.set_aspect("equal")
ax1.set_title("Factual", fontsize=12, fontweight="bold")

# Plot CF-8s
ax2 = axes[1]
im2 = zsmax_cf8_plot.plot(
    ax=ax2,
    levels=flood_levels,
    cmap=flood_cmap,
    add_colorbar=False,
    x="x",
    y="y",
    alpha=0.7,
    zorder=2,
)
ax2.set_aspect("equal")
ax2.set_title("Counterfactual", fontsize=12, fontweight="bold")

# Plot difference
ax3 = axes[2]
im3 = diff.plot(
    ax=ax3,
    levels=diff_levels,
    cmap=diff_cmap,
    add_colorbar=False,
    x="x",
    y="y",
    alpha=0.7,
    zorder=2,
)
ax3.set_aspect("equal")
ax3.set_title(
    "Flood Depth Changes: Factual vs Counterfactual", fontsize=12, fontweight="bold"
)

# Add satellite basemap to all axes
for ax in axes:
    try:
        ctx.add_basemap(
            ax=ax,
            source=ctx.providers.OpenStreetMap.Mapnik,
            crs=data_crs,
            attribution=False,
            zorder=1,
            zoom=12,  # Higher zoom = higher resolution
        )
    except Exception as e:
        print(f"Could not add basemap: {e}")

# Add city markers to all axes
extent = (
    zsmax_cf0_plot.x.min().item(),
    zsmax_cf0_plot.x.max().item(),
    zsmax_cf0_plot.y.min().item(),
    zsmax_cf0_plot.y.max().item(),
)
print(
    f"Plot extent: lon=[{extent[0]:.2f}, {extent[1]:.2f}], lat=[{extent[2]:.2f}, {extent[3]:.2f}]"
)

for ax in axes:
    num_cities = add_city_markers(
        ax, extent, EVENT_CITY_LIST, fontsize=7, markersize=40
    )

if num_cities > 0:
    print(f"Added {num_cities} city markers to plots")
else:
    print("No cities fall within the plot extent")

# Format axis labels as lat/lon
for ax in axes:
    ax.set_xlabel("Longitude [°]", fontsize=10)
    ax.set_ylabel("Latitude [°]", fontsize=10)

# Add colorbars
cbar1 = plt.colorbar(im1, ax=ax1, shrink=0.8, pad=0.1)
cbar1.set_label("Max Flood Depth [m]", fontsize=10)

cbar2 = plt.colorbar(im2, ax=ax2, shrink=0.8, pad=0.1)
cbar2.set_label("Max Flood Depth [m]", fontsize=10)

cbar3 = plt.colorbar(im3, ax=ax3, shrink=0.8, pad=0.1)
cbar3.set_label("Depth Difference [m]", fontsize=10)

# Adjust layout and add main title
plt.tight_layout()
plt.suptitle(
    f"Flood Depth Analysis - {EVENT_NAME}: Factual vs Counterfactual",
    fontsize=16,
    fontweight="bold",
    y=1.02,
)

# Save the figure
output_file_main = (
    OUTPUT_DIR / f"sfincs_{EVENT_NAME.lower()}_flood_depth_comparison.png"
)
plt.savefig(output_file_main, dpi=300, bbox_inches="tight")
plt.show()

# ===== CREATE BAR CHARTS FOR AGGREGATED FLOOD METRICS =====
print("Creating aggregated flood metrics bar charts...")


def create_bar_chart_with_annotation(
    ax, scenarios, totals, ylabel, title, unit_prefix="", unit_suffix=""
):
    """Helper function to create a bar chart with climate change attribution annotation."""
    bars = ax.bar(
        scenarios, totals, color=["steelblue", "steelblue"], alpha=1, width=0.4
    )

    ax.set_ylabel(ylabel, fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_axisbelow(True)

    # Calculate difference for annotation
    difference = abs(totals[1] - totals[0])
    percentage_diff = (
        (difference / min(totals[0], totals[1])) * 100
        if min(totals[0], totals[1]) > 0
        else 0
    )

    # Get bar positions
    bar_positions = [bar.get_x() + bar.get_width() / 2 for bar in bars]
    bar_heights = [bar.get_height() for bar in bars]

    # Determine which bar is higher
    higher_bar_idx = 0 if totals[0] > totals[1] else 1
    lower_bar_idx = 1 - higher_bar_idx

    # Annotation parameters for vertical difference line
    line_x_position = bar_positions[1] + (bar_positions[1] - bar_positions[0]) * 0.3
    line_extension = (bar_positions[1] - bar_positions[0]) * 0.04

    # Draw the vertical difference line
    ax.plot(
        [line_x_position, line_x_position],
        [bar_heights[lower_bar_idx], bar_heights[higher_bar_idx]],
        "k-",
        linewidth=2,
    )

    # Add extended dashed horizontal lines
    left_edge = ax.get_xlim()[0]
    for total in totals:
        ax.plot(
            [left_edge, line_x_position],
            [total, total],
            linestyle="--",
            color="gray",
            alpha=0.7,
            linewidth=1,
        )

    # Small horizontal ticks at both ends
    ax.plot(
        [line_x_position - line_extension, line_x_position + line_extension],
        [bar_heights[lower_bar_idx], bar_heights[lower_bar_idx]],
        "k-",
        linewidth=1.5,
    )
    ax.plot(
        [line_x_position - line_extension, line_x_position + line_extension],
        [bar_heights[higher_bar_idx], bar_heights[higher_bar_idx]],
        "k-",
        linewidth=1.5,
    )

    # Add text annotation
    text_x = line_x_position + (bar_positions[1] - bar_positions[0]) * 0.1

    # Calculate midpoint between the two bar heights
    midpoint_y = (bar_heights[lower_bar_idx] + bar_heights[higher_bar_idx]) / 2

    # Calculate spacing based on the difference between bars (with a minimum)
    bar_diff = bar_heights[higher_bar_idx] - bar_heights[lower_bar_idx]
    text_spacing = max(bar_diff * 0.15, max(totals) * 0.03)  # Ensure minimum spacing

    # Format the difference text
    diff_text = f"{unit_prefix}{difference:.2f}{unit_suffix} ({percentage_diff:.0f}%)"

    # Place the value text slightly below midpoint
    ax.text(text_x, midpoint_y - text_spacing, diff_text, ha="left", va="center", fontsize=10, rotation=0)

    # Place the label text slightly above midpoint
    ax.text(
        text_x,
        midpoint_y + text_spacing,
        "Climate change\nattribution:",
        ha="left",
        va="center",
        fontsize=9,
        style="italic",
        rotation=0,
    )

    # Set y-axis to start from 0 with appropriate margin
    ax.set_ylim(0, max(totals) * 1.15)

    # Remove top and right spines for cleaner look
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ===== FLOOD EXTENT BAR CHART =====
print("Creating flood extent bar chart...")
fig_extent, ax_extent = plt.subplots(1, 1, figsize=(6, 6))

scenarios = ["Counterfactual", "Factual"]
extent_totals = [area_cf8_m2 / 1e6, area_cf0_m2 / 1e6]  # Convert to km²

create_bar_chart_with_annotation(
    ax_extent,
    scenarios,
    extent_totals,
    ylabel="Flood Extent [km²]",
    title=f"Flood Extent - {EVENT_NAME}",
    unit_prefix="",
    unit_suffix=" km²",
)

plt.tight_layout()
output_file_extent_bar = (
    OUTPUT_DIR / f"sfincs_{EVENT_NAME.lower()}_flood_extent_barchart.png"
)
plt.savefig(output_file_extent_bar, dpi=300, bbox_inches="tight")
plt.show()

# ===== FLOOD VOLUME BAR CHART =====
print("Creating flood volume bar chart...")
fig_volume, ax_volume = plt.subplots(1, 1, figsize=(6, 6))

volume_totals = [volume_cf8_m3 / 1e6, volume_cf0_m3 / 1e6]  # Convert to Mm³

create_bar_chart_with_annotation(
    ax_volume,
    scenarios,
    volume_totals,
    ylabel="Flood Volume [Mm³]",
    title=f"Flood Volume - {EVENT_NAME}",
    unit_prefix="",
    unit_suffix=" Mm³",
)

plt.tight_layout()
output_file_volume_bar = (
    OUTPUT_DIR / f"sfincs_{EVENT_NAME.lower()}_flood_volume_barchart.png"
)
plt.savefig(output_file_volume_bar, dpi=300, bbox_inches="tight")
plt.show()

# ===== EXPORT AGGREGATED FLOOD DATA TO CSV =====
print("Exporting aggregated flood data to CSV...")

# Create DataFrame with flood extent and volume data
flood_data = pd.DataFrame(
    {
        "Scenario": ["Counterfactual", "Factual", "Difference", "Percent_Change"],
        "Flood_Extent_km2": [
            area_cf8_m2 / 1e6,
            area_cf0_m2 / 1e6,
            extent_area_diff_m2 / 1e6,
            extent_area_pct,
        ],
        "Flood_Volume_Mm3": [
            volume_cf8_m3 / 1e6,
            volume_cf0_m3 / 1e6,
            volume_diff_m3 / 1e6,
            volume_pct,
        ],
    }
)

# Add metadata columns
flood_data["Event"] = EVENT_NAME
flood_data["Flood_Threshold_m"] = flood_threshold

# Reorder columns
flood_data = flood_data[
    ["Event", "Scenario", "Flood_Extent_km2", "Flood_Volume_Mm3", "Flood_Threshold_m"]
]

# Save to CSV
output_file_csv = OUTPUT_DIR / f"sfincs_{EVENT_NAME.lower()}_flood_aggregated_data.csv"
flood_data.to_csv(output_file_csv, index=False)
print(f"Aggregated flood data saved to: {output_file_csv}")
print(flood_data.to_string(index=False))

# ===== DETAILED DIFFERENCE PLOT =====
print("Creating detailed difference plot...")
fig2, ax = plt.subplots(1, 1, figsize=(12, 8))

# Plot only the difference with better resolution
im_diff = diff.plot(
    ax=ax,
    levels=diff_levels,
    cmap="Reds",
    add_colorbar=False,
    x="x",
    y="y",
    alpha=0.7,
    zorder=2,
)
ax.set_aspect("equal")

# Add satellite basemap
try:
    ctx.add_basemap(
        ax=ax,
        source=ctx.providers.OpenStreetMap.Mapnik,
        crs=data_crs,
        attribution=False,
        zorder=1,
        zoom=12,  # Higher zoom = higher resolution
    )
except Exception as e:
    print(f"Could not add basemap: {e}")

# Add city markers
add_city_markers(ax, extent, EVENT_CITY_LIST, fontsize=9, markersize=60)

# Format axis labels as lat/lon
ax.set_xlabel("Longitude [°]", fontsize=12)
ax.set_ylabel("Latitude [°]", fontsize=12)

# Colorbar
cbar = plt.colorbar(im_diff, ax=ax, shrink=0.8, pad=0.1)
cbar.set_label("Flood Depth Difference [m]", fontsize=12)

# Title
ax.set_title(
    f"Flood Depth Changes - {EVENT_NAME}: Factual vs Counterfactual",
    fontsize=14,
    fontweight="bold",
    pad=20,
)

# Save
output_file_diff = OUTPUT_DIR / f"sfincs_{EVENT_NAME.lower()}_difference_only.png"
plt.savefig(output_file_diff, dpi=300, bbox_inches="tight")
plt.show()

# ===== INDIVIDUAL FACTUAL PLOT =====
print("Creating individual factual plot...")
fig3, ax_f = plt.subplots(1, 1, figsize=(12, 8))

im_factual = zsmax_cf0_plot.plot(
    ax=ax_f,
    levels=flood_levels,
    cmap=flood_cmap,
    add_colorbar=False,
    x="x",
    y="y",
    alpha=0.7,
    zorder=2,
)
ax_f.set_aspect("equal")

# Add satellite basemap
try:
    ctx.add_basemap(
        ax=ax_f,
        source=ctx.providers.OpenStreetMap.Mapnik,
        crs=data_crs,
        attribution=False,
        zorder=1,
        zoom=12,  # Higher zoom = higher resolution
    )
except Exception as e:
    print(f"Could not add basemap: {e}")

# Add city markers
add_city_markers(ax_f, extent, EVENT_CITY_LIST, fontsize=9, markersize=60)

# Format axis labels as lat/lon
ax_f.set_xlabel("Longitude [°]", fontsize=12)
ax_f.set_ylabel("Latitude [°]", fontsize=12)

# Colorbar
cbar_f = plt.colorbar(im_factual, ax=ax_f, shrink=0.8, pad=0.1)
cbar_f.set_label("Max Flood Depth [m]", fontsize=12)

# Title
ax_f.set_title(
    f"Factual Flood Depth - {EVENT_NAME}",
    fontsize=14,
    fontweight="bold",
    pad=20,
)

# Save
output_file_factual = OUTPUT_DIR / f"sfincs_{EVENT_NAME.lower()}_factual.png"
plt.savefig(output_file_factual, dpi=300, bbox_inches="tight")
plt.show()

# ===== INDIVIDUAL COUNTERFACTUAL PLOT =====
print("Creating individual counterfactual plot...")
fig4, ax_cf = plt.subplots(1, 1, figsize=(12, 8))

im_counterfactual = zsmax_cf8_plot.plot(
    ax=ax_cf,
    levels=flood_levels,
    cmap=flood_cmap,
    add_colorbar=False,
    x="x",
    y="y",
    alpha=0.7,
    zorder=2,
)
ax_cf.set_aspect("equal")

# Add satellite basemap
try:
    ctx.add_basemap(
        ax=ax_cf,
        source=ctx.providers.OpenStreetMap.Mapnik,
        crs=data_crs,
        attribution=False,
        zorder=1,
        zoom=12,  # Higher zoom = higher resolution
    )
except Exception as e:
    print(f"Could not add basemap: {e}")

# Add city markers
add_city_markers(ax_cf, extent, EVENT_CITY_LIST, fontsize=9, markersize=60)

# Format axis labels as lat/lon
ax_cf.set_xlabel("Longitude [°]", fontsize=12)
ax_cf.set_ylabel("Latitude [°]", fontsize=12)

# Colorbar
cbar_cf = plt.colorbar(im_counterfactual, ax=ax_cf, shrink=0.8, pad=0.1)
cbar_cf.set_label("Max Flood Depth [m]", fontsize=12)

# Title
ax_cf.set_title(
    f"Counterfactual Flood Depth - {EVENT_NAME}",
    fontsize=14,
    fontweight="bold",
    pad=20,
)

# Save
output_file_counterfactual = (
    OUTPUT_DIR / f"sfincs_{EVENT_NAME.lower()}_counterfactual.png"
)
plt.savefig(output_file_counterfactual, dpi=300, bbox_inches="tight")
plt.show()

# Close datasets
zsmax_cf0.close()
zsmax_cf8.close()

print(f"Analysis complete. Output figures:")
print(f"  - {output_file_main}")
print(f"  - {output_file_extent_bar}")
print(f"  - {output_file_volume_bar}")
print(f"  - {output_file_diff}")
print(f"  - {output_file_factual}")
print(f"  - {output_file_counterfactual}")
plt.show()

# ===== SUMMARY STATISTICS =====
print(f"\nAreas with significant changes:")
print(f"Grid cells with >0.1m increase: {(diff > 0.1).sum().values}")
print(f"Grid cells with >0.1m decrease: {(diff < -0.1).sum().values}")
print(f"Grid cells with >0.5m increase: {(diff > 0.5).sum().values}")
print(f"Grid cells with >0.5m decrease: {(diff < -0.5).sum().values}")

print("\nPercent change exceedances (same mask as above):")
print(f"Cells with >10% increase: {(percent_diff > 10).sum().values}")
print(f"Cells with >50% increase: {(percent_diff > 50).sum().values}")
print(f"Cells with >100% increase: {(percent_diff > 100).sum().values}")
print(f"Cells with >10% decrease: {(percent_diff < -10).sum().values}")
print(f"Cells with >50% decrease: {(percent_diff < -50).sum().values}")

# Close datasets
zsmax_cf0.close()
zsmax_cf8.close()

print(
    f"\nAnalysis complete for {EVENT_NAME}! Check the saved PNG files in: {OUTPUT_DIR}"
)
print("Files created:")
print(f"  - {output_file_main}")
print(f"  - {output_file_extent_bar}")
print(f"  - {output_file_volume_bar}")
print(f"  - {output_file_csv}")
print(f"  - {output_file_diff}")
print(f"  - {output_file_factual}")
print(f"  - {output_file_counterfactual}")

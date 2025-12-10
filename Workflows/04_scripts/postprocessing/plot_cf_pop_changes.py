import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import BoundaryNorm, ListedColormap
import contextily as ctx
import warnings

warnings.filterwarnings("ignore")

# File paths
from pathlib import Path

# by @dumontgoulart - adapted for population exposure analysis
# TODO: this is a preliminary script that is still not fully embedded in the snakemake workflow. The user is required to change the events & input files for now.

# ===== CONFIGURATION =====
# Set your event name here
EVENT_NAME = "Kenneth"  # Change this to: "Kenneth", "Freddy", etc.

# Choose population column to analyze
POPULATION_COLUMN = "population"  # Main population column to analyze

# Inundation depth threshold (in meters)
INUNDATION_THRESHOLD = 0.2  # Only consider areas with >0.2m flooding

# Base paths - update these as needed
BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs/test")
OUTPUT_DIR = Path("/p/11210471-001-compass/04_Results/CF_figs")

# ===== DYNAMIC FILE PATHS =====
# Construct file paths based on event name
file_cf0 = (
    BASE_RUN_PATH
    / EVENT_NAME
    / "fiat"
    / "event_tp_era5_hourly_zarr_CF0_GTSMv41opendap_CF0_no_wind_CF0"
    / "output"
    / "spatial_with_pop_and_flood.fgb"
)
file_cf8 = (
    BASE_RUN_PATH
    / EVENT_NAME
    / "fiat"
    / "event_tp_era5_hourly_zarr_CF-8_GTSMv41opendap_CF0_no_wind_CF0"
    / "output"
    / "spatial_with_pop_and_flood.fgb"
)

# Create output directory if it doesn't exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Processing population exposure for event: {EVENT_NAME}")
print(f"Using population column: {POPULATION_COLUMN}")
print(f"Inundation depth threshold: {INUNDATION_THRESHOLD}m")

# ===== LOAD POPULATION EXPOSURE DATA =====
print("Loading population exposure data files...")
print(f"Loading CF0: {file_cf0}")
print(f"Loading CF-8: {file_cf8}")

# Read the geodataframes
gdf_cf0 = gpd.read_file(file_cf0)
gdf_cf8 = gpd.read_file(file_cf8)

print(f"CF0 data shape: {gdf_cf0.shape}")
print(f"CF-8 data shape: {gdf_cf8.shape}")
print(f"CF0 columns: {list(gdf_cf0.columns)}")
print(f"CF-8 columns: {list(gdf_cf8.columns)}")

# Check if the population column exists
if POPULATION_COLUMN not in gdf_cf0.columns:
    print(
        f"Warning: {POPULATION_COLUMN} not found in CF0 data. Available columns: {list(gdf_cf0.columns)}"
    )
    # Try alternative column names
    if "pop" in gdf_cf0.columns:
        POPULATION_COLUMN = "pop"
    elif "total_population" in gdf_cf0.columns:
        POPULATION_COLUMN = "total_population"
    else:
        raise ValueError("No suitable population column found!")

print(f"Final population column used: {POPULATION_COLUMN}")

# Check if inundation depth column exists
if "inun_depth" not in gdf_cf0.columns:
    print(
        "Warning: 'inun_depth' column not found. Available columns with 'depth' or 'inun':"
    )
    depth_cols = [
        col
        for col in gdf_cf0.columns
        if "depth" in col.lower() or "inun" in col.lower()
    ]
    print(depth_cols)
    if depth_cols:
        inun_col = depth_cols[0]
        print(f"Using column: {inun_col}")
    else:
        raise ValueError("No suitable inundation depth column found!")
else:
    inun_col = "inun_depth"

# ===== FILTER BY INUNDATION DEPTH =====
print(f"Filtering data by inundation depth > {INUNDATION_THRESHOLD}m...")
print(f"CF0 records before filtering: {len(gdf_cf0)}")
print(f"CF-8 records before filtering: {len(gdf_cf8)}")

# Apply inundation depth filter
gdf_cf0_filtered = gdf_cf0[gdf_cf0[inun_col] > INUNDATION_THRESHOLD].copy()
gdf_cf8_filtered = gdf_cf8[gdf_cf8[inun_col] > INUNDATION_THRESHOLD].copy()

print(f"CF0 records after filtering: {len(gdf_cf0_filtered)}")
print(f"CF-8 records after filtering: {len(gdf_cf8_filtered)}")

# ===== CONVERT TO LAT/LON =====
print("Converting data to lat/lon coordinates (EPSG:4326)...")
print(f"Original CRS: {gdf_cf0.crs}")

# Convert to EPSG:4326 (WGS84 lat/lon) if not already
if gdf_cf0.crs != "EPSG:4326":
    gdf_cf0_filtered = gdf_cf0_filtered.to_crs("EPSG:4326")
    gdf_cf8_filtered = gdf_cf8_filtered.to_crs("EPSG:4326")
    print("Converted to EPSG:4326")
else:
    print("Already in EPSG:4326")

# ===== EXTRACT COORDINATES =====
print("Extracting coordinates from geometry centroids...")
# Extract x, y coordinates from geometry centroids
gdf_cf0_filtered["centroid"] = gdf_cf0_filtered.geometry.centroid
gdf_cf0_filtered["x"] = gdf_cf0_filtered["centroid"].x
gdf_cf0_filtered["y"] = gdf_cf0_filtered["centroid"].y

gdf_cf8_filtered["centroid"] = gdf_cf8_filtered.geometry.centroid
gdf_cf8_filtered["x"] = gdf_cf8_filtered["centroid"].x
gdf_cf8_filtered["y"] = gdf_cf8_filtered["centroid"].y

print(
    f"Sample coordinates CF0: lon={gdf_cf0_filtered['x'].iloc[0]:.4f}°, lat={gdf_cf0_filtered['y'].iloc[0]:.4f}°"
)
print(
    f"Sample coordinates CF-8: lon={gdf_cf8_filtered['x'].iloc[0]:.4f}°, lat={gdf_cf8_filtered['y'].iloc[0]:.4f}°"
)

# ===== MERGE DATA FOR DIFFERENCE CALCULATION =====
print("Merging data for difference calculation...")
# Merge on object_id to compare same locations
# Use coordinates from CF0 data, but handle cases where locations might only exist in one scenario
merged = gdf_cf0_filtered[["object_id", "x", "y", POPULATION_COLUMN, inun_col]].merge(
    gdf_cf8_filtered[["object_id", "x", "y", POPULATION_COLUMN, inun_col]],
    on="object_id",
    how="outer",
    suffixes=("_cf0", "_cf8"),
)

# Handle coordinates: use CF0 coordinates where available, otherwise CF8
merged["x"] = merged["x_cf0"].fillna(merged["x_cf8"])
merged["y"] = merged["y_cf0"].fillna(merged["y_cf8"])

# Handle inundation depth: use CF0 values where available, otherwise CF8
merged["inun_depth"] = merged[f"{inun_col}_cf0"].fillna(merged[f"{inun_col}_cf8"])

# Drop the temporary coordinate columns
merged = merged.drop(
    [col for col in merged.columns if col.endswith("_cf0") or col.endswith("_cf8")],
    axis=1,
)

# Re-add the population columns with proper names
merged[f"{POPULATION_COLUMN}_cf0"] = (
    gdf_cf0_filtered.set_index("object_id")[POPULATION_COLUMN]
    .reindex(merged["object_id"])
    .fillna(0)
    .values
)
merged[f"{POPULATION_COLUMN}_cf8"] = (
    gdf_cf8_filtered.set_index("object_id")[POPULATION_COLUMN]
    .reindex(merged["object_id"])
    .fillna(0)
    .values
)

# Calculate difference in exposed population (CF0 - CF8)
merged["pop_exposure_diff"] = (
    merged[f"{POPULATION_COLUMN}_cf0"] - merged[f"{POPULATION_COLUMN}_cf8"]
)

# Remove locations with no coordinates
merged = merged.dropna(subset=["x", "y"])

print(f"Merged data shape: {merged.shape}")
print(f"Sample merged data:")
print(
    merged[
        [
            "object_id",
            "x",
            "y",
            f"{POPULATION_COLUMN}_cf0",
            f"{POPULATION_COLUMN}_cf8",
            "pop_exposure_diff",
        ]
    ].head()
)

# ===== STATISTICS =====
print(f"\nPopulation Exposure Statistics for {EVENT_NAME}:")
cf0_pop = merged[f"{POPULATION_COLUMN}_cf0"]
cf8_pop = merged[f"{POPULATION_COLUMN}_cf8"]
pop_diff = merged["pop_exposure_diff"]

print(f"CF0 max exposed population: {cf0_pop.max():.0f} people")
print(f"CF0 mean exposed population: {cf0_pop.mean():.1f} people")
print(f"CF0 total exposed population: {cf0_pop.sum():.0f} people")
print(f"CF0 locations with >0 exposed population: {(cf0_pop > 0).sum()}")
print(f"CF-8 max exposed population: {cf8_pop.max():.0f} people")
print(f"CF-8 mean exposed population: {cf8_pop.mean():.1f} people")
print(f"CF-8 total exposed population: {cf8_pop.sum():.0f} people")
print(f"CF-8 locations with >0 exposed population: {(cf8_pop > 0).sum()}")
print(
    f"Difference in total exposed population: {cf0_pop.sum() - cf8_pop.sum():.0f} people"
)
if cf8_pop.sum() > 0:
    print(
        f"Percentage change: {((cf0_pop.sum() - cf8_pop.sum()) / cf8_pop.sum() * 100):.1f}%"
    )
else:
    print("Percentage change: Cannot calculate (CF-8 has no exposed population)")
print(f"Max increase in exposed population (CF0 vs CF-8): {pop_diff.max():.0f} people")
print(f"Max decrease in exposed population (CF0 vs CF-8): {pop_diff.min():.0f} people")
print(f"Mean difference in exposed population: {pop_diff.mean():.1f} people")

# ===== DETERMINE COORDINATE SYSTEM =====
try:
    # Get approximate center coordinates
    center_x = merged["x"].mean()
    center_y = merged["y"].mean()

    print(f"Center coordinates: {center_x:.0f}, {center_y:.0f}")

    # Determine coordinate system based on coordinate values
    if abs(center_x) > 180 or abs(center_y) > 90:
        # Likely projected coordinates (UTM, etc.)
        if center_y < 0:
            # Southern hemisphere
            utm_zone = int((center_x + 180) / 6) + 1
            southern = True
        else:
            # Northern hemisphere
            utm_zone = int((center_x + 180) / 6) + 1
            southern = False
        print(
            f"Detected projected coordinates - UTM zone {utm_zone}, Southern: {southern}"
        )
    else:
        # Likely geographic coordinates (lat/lon)
        print("Detected geographic coordinates (lat/lon)")
        utm_zone = None
        southern = None

except Exception as e:
    print(f"Could not determine coordinates: {e}")
    utm_zone = 37
    southern = True

# ===== PLOTTING SETUP =====
print("Setting up plots...")

# Data is now in EPSG:4326
data_crs = "EPSG:4326"
print(f"Data CRS for plotting: {data_crs}")

# Define colormaps and levels with discrete intervals - using blues for population
# Population exposure colormap
max_pop = max(cf0_pop.quantile(0.95), cf8_pop.quantile(0.95))
if max_pop > 1000:
    # Create 10 intervals from 0 to max
    boundaries = np.linspace(0, max_pop, 11)
elif max_pop > 100:
    # For smaller populations, use smaller intervals
    boundaries = np.linspace(0, max_pop, 11)
else:
    # For very small populations
    boundaries = np.linspace(0, max_pop, 11)

pop_norm = BoundaryNorm(boundaries, ncolors=256, clip=True)
pop_cmap = plt.get_cmap("Blues")  # Blue colormap for population
pop_label = "Exposed Population [people]"
vmin, vmax = 0, max_pop

# Difference colormap with discrete intervals
diff_max = max(abs(pop_diff.quantile(0.05)), abs(pop_diff.quantile(0.95)))
diff_boundaries = np.linspace(-diff_max, diff_max, 11)
diff_norm = BoundaryNorm(diff_boundaries, ncolors=256, clip=True)
diff_cmap = plt.get_cmap("RdBu_r")  # Red for increases, Blue for decreases

print(f"Population boundaries: {boundaries}")
print(f"Difference boundaries: {diff_boundaries}")

# ===== CREATE MAIN COMPARISON PLOT (3-PANEL) =====
print("Creating population exposure comparison plots...")

# Use regular matplotlib axes (not Cartopy) for contextily compatibility
fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# Plot settings
point_size = 20  # Point size for visibility
alpha = 0.8  # Transparency for better colors

# Plot CF0 (Factual)
scatter1 = axes[0].scatter(
    merged["x"],
    merged["y"],
    c=merged[f"{POPULATION_COLUMN}_cf0"],
    cmap=pop_cmap,
    norm=pop_norm,
    s=point_size,
    alpha=alpha,
    zorder=2,
)
axes[0].set_aspect("equal")
axes[0].set_title("Factual", fontsize=12, fontweight="bold")

# Plot CF-8 (Counterfactual)
scatter2 = axes[1].scatter(
    merged["x"],
    merged["y"],
    c=merged[f"{POPULATION_COLUMN}_cf8"],
    cmap=pop_cmap,
    norm=pop_norm,
    s=point_size,
    alpha=alpha,
    zorder=2,
)
axes[1].set_aspect("equal")
axes[1].set_title("Counterfactual", fontsize=12, fontweight="bold")

# Plot Difference
scatter3 = axes[2].scatter(
    merged["x"],
    merged["y"],
    c=merged["pop_exposure_diff"],
    cmap=diff_cmap,
    norm=diff_norm,
    s=point_size,
    alpha=alpha,
    zorder=2,
)
axes[2].set_aspect("equal")
axes[2].set_title(
    "Population Exposure Changes: Factual vs Counterfactual",
    fontsize=12,
    fontweight="bold",
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

# Format axis labels as lat/lon
for ax in axes:
    ax.set_xlabel("Longitude [°]", fontsize=10)
    ax.set_ylabel("Latitude [°]", fontsize=10)

# Add colorbars
cbar1 = plt.colorbar(scatter1, ax=axes[0], shrink=0.8, pad=0.1)
cbar1.set_label(pop_label, fontsize=10)

cbar2 = plt.colorbar(scatter2, ax=axes[1], shrink=0.8, pad=0.1)
cbar2.set_label(pop_label, fontsize=10)

cbar3 = plt.colorbar(scatter3, ax=axes[2], shrink=0.8, pad=0.1)
cbar3.set_label(f"{pop_label} Difference", fontsize=10)

# Adjust layout and add main title
plt.tight_layout()
plt.suptitle(
    f"Population Exposure Analysis - {EVENT_NAME}: Factual vs Counterfactual (depth > {INUNDATION_THRESHOLD}m)",
    fontsize=16,
    fontweight="bold",
    y=1.02,
)

# Save the figure
output_file_main = (
    OUTPUT_DIR / f"population_exposure_{EVENT_NAME.lower()}_comparison.png"
)
plt.savefig(output_file_main, dpi=300, bbox_inches="tight")
plt.show()

# ===== CREATE SEPARATE BAR CHART FOR TOTAL EXPOSED POPULATION =====
print("Creating separate total exposed population bar chart...")
# Calculate total exposed population for bar chart
total_cf0 = cf0_pop.sum()
total_cf8 = cf8_pop.sum()
bar_label = "Total Exposed Population [people]"

# Create separate bar chart figure
fig_bar, ax_bar = plt.subplots(1, 1, figsize=(6, 6))

# Create bar chart
scenarios = ["Counterfactual", "Factual"]
totals = [total_cf8, total_cf0]
bars = ax_bar.bar(
    scenarios, totals, color=["lightblue", "lightblue"], alpha=1, width=0.4
)

# Style the bar chart
ax_bar.set_ylabel(bar_label, fontsize=12, fontweight="bold")
ax_bar.set_title(
    f"Total Exposed Population - {EVENT_NAME}", fontsize=14, fontweight="bold"
)
ax_bar.set_axisbelow(True)

# Calculate difference for annotation
difference = abs(total_cf0 - total_cf8)
percentage_diff = (
    (difference / min(total_cf0, total_cf8)) * 100
    if min(total_cf0, total_cf8) > 0
    else 0
)

# Add annotation showing climate change attribution
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
ax_bar.plot(
    [line_x_position, line_x_position],
    [bar_heights[lower_bar_idx], bar_heights[higher_bar_idx]],
    "k-",
    linewidth=2,
)

# Add extended dashed horizontal lines from y-axis to the vertical difference line
left_edge = ax_bar.get_xlim()[0]
for total in totals:
    ax_bar.plot(
        [left_edge, line_x_position],
        [total, total],
        linestyle="--",
        color="gray",
        alpha=0.7,
        linewidth=1,
    )

# Small horizontal ticks at both ends
ax_bar.plot(
    [line_x_position - line_extension, line_x_position + line_extension],
    [bar_heights[lower_bar_idx], bar_heights[lower_bar_idx]],
    "k-",
    linewidth=1.5,
)
ax_bar.plot(
    [line_x_position - line_extension, line_x_position + line_extension],
    [bar_heights[higher_bar_idx], bar_heights[higher_bar_idx]],
    "k-",
    linewidth=1.5,
)

# Add text annotation positioned at the center between bars
text_x = line_x_position + (bar_positions[1] - bar_positions[0]) * 0.1
text_y = ((bar_heights[lower_bar_idx] + bar_heights[higher_bar_idx]) / 2) * 0.97

# Format the difference text for population
if difference >= 1e6:
    diff_text = f"{difference/1e6:.1f}M people ({percentage_diff:.0f}%)"
elif difference >= 1e3:
    diff_text = f"{difference/1e3:.1f}K people ({percentage_diff:.0f}%)"
else:
    diff_text = f"{difference:.0f} people ({percentage_diff:.0f}%)"

# Add the annotation text
ax_bar.text(text_x, text_y, diff_text, ha="left", va="center", fontsize=10, rotation=0)
text_offset = (bar_heights[higher_bar_idx] - bar_heights[lower_bar_idx]) * 0.35
ax_bar.text(
    text_x,
    text_y + text_offset,
    "Climate change \n attribution:",
    ha="left",
    va="center",
    fontsize=9,
    style="italic",
    rotation=0,
)

# Set y-axis to start from 0 with appropriate margin
ax_bar.set_ylim(0, max(totals) * 1.1)

# Remove top and right spines for cleaner look
ax_bar.spines["top"].set_visible(False)
ax_bar.spines["right"].set_visible(False)

plt.tight_layout()

# Save the bar chart
output_file_bar = (
    OUTPUT_DIR / f"population_exposure_{EVENT_NAME.lower()}_totals_barchart.png"
)
plt.savefig(output_file_bar, dpi=300, bbox_inches="tight")
plt.show()

# ===== DETAILED DIFFERENCE PLOT =====
print("Creating detailed population exposure difference plot...")

fig2, ax = plt.subplots(1, 1, figsize=(12, 8))

# Plot only the difference with larger points
scatter_diff = ax.scatter(
    merged["x"],
    merged["y"],
    c=merged["pop_exposure_diff"],
    cmap=diff_cmap,
    norm=diff_norm,
    s=point_size * 2,
    alpha=alpha,
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
        zoom=14,
    )
except Exception as e:
    print(f"Could not add basemap: {e}")

# Format axis labels as lat/lon
ax.set_xlabel("Longitude [°]", fontsize=12)
ax.set_ylabel("Latitude [°]", fontsize=12)

# Colorbar
cbar = plt.colorbar(scatter_diff, ax=ax, shrink=0.8, pad=0.1)
cbar.set_label(f"{pop_label} Difference", fontsize=12)

# Title
ax.set_title(
    f"Population Exposure Changes - {EVENT_NAME}: Factual vs Counterfactual\n(Red = More exposed people, Blue = Fewer exposed people)",
    fontsize=14,
    fontweight="bold",
    pad=20,
)

# Save
output_file_diff = (
    OUTPUT_DIR / f"population_exposure_{EVENT_NAME.lower()}_difference_only.png"
)
plt.savefig(output_file_diff, dpi=300, bbox_inches="tight")
plt.show()

# ===== CREATE FACTUAL-ONLY PLOT =====
print("Creating factual-only population exposure plot...")

fig_factual, ax_factual = plt.subplots(1, 1, figsize=(10, 8))

# Plot only the factual (CF0) data
scatter_factual = ax_factual.scatter(
    merged["x"],
    merged["y"],
    c=merged[f"{POPULATION_COLUMN}_cf0"],
    cmap=pop_cmap,
    norm=pop_norm,
    s=point_size * 1.5,
    alpha=alpha,
    zorder=2,
)
ax_factual.set_aspect("equal")

# Add satellite basemap
try:
    ctx.add_basemap(
        ax=ax_factual,
        source=ctx.providers.OpenStreetMap.Mapnik,
        crs=data_crs,
        attribution=False,
        zorder=1,
        zoom=14,
    )
except Exception as e:
    print(f"Could not add basemap: {e}")

# Format axis labels as lat/lon
ax_factual.set_xlabel("Longitude [°]", fontsize=12)
ax_factual.set_ylabel("Latitude [°]", fontsize=12)

# Colorbar
cbar_factual = plt.colorbar(scatter_factual, ax=ax_factual, shrink=0.8, pad=0.1)
cbar_factual.set_label(pop_label, fontsize=12)

# Title
ax_factual.set_title(
    f"Population Exposure - {EVENT_NAME}: Factual Scenario\n(depth > {INUNDATION_THRESHOLD}m)",
    fontsize=14,
    fontweight="bold",
    pad=20,
)

# Save the factual-only plot
output_file_factual = (
    OUTPUT_DIR / f"population_exposure_{EVENT_NAME.lower()}_factual_only.png"
)
plt.savefig(output_file_factual, dpi=300, bbox_inches="tight")
plt.show()

# ===== CREATE COUNTERFACTUAL-ONLY PLOT =====
print("Creating counterfactual-only population exposure plot...")

fig_cf, ax_cf = plt.subplots(1, 1, figsize=(10, 8))

# Plot only the counterfactual (CF-8) data
scatter_cf = ax_cf.scatter(
    merged["x"],
    merged["y"],
    c=merged[f"{POPULATION_COLUMN}_cf8"],
    cmap=pop_cmap,
    norm=pop_norm,
    s=point_size * 1.5,
    alpha=alpha,
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
        zoom=14,
    )
except Exception as e:
    print(f"Could not add basemap: {e}")

# Format axis labels as lat/lon
ax_cf.set_xlabel("Longitude [°]", fontsize=12)
ax_cf.set_ylabel("Latitude [°]", fontsize=12)

# Colorbar
cbar_cf = plt.colorbar(scatter_cf, ax=ax_cf, shrink=0.8, pad=0.1)
cbar_cf.set_label(pop_label, fontsize=12)

# Title
ax_cf.set_title(
    f"Population Exposure - {EVENT_NAME}: Counterfactual Scenario\n(depth > {INUNDATION_THRESHOLD}m)",
    fontsize=14,
    fontweight="bold",
    pad=20,
)

# Save the counterfactual-only plot
output_file_counterfactual = (
    OUTPUT_DIR / f"population_exposure_{EVENT_NAME.lower()}_counterfactual_only.png"
)
plt.savefig(output_file_counterfactual, dpi=300, bbox_inches="tight")
plt.show()

# ===== SUMMARY STATISTICS =====
print(f"\nLocations with significant population exposure changes:")
# Define thresholds based on population scale
low_threshold = max(10, diff_max * 0.1)
high_threshold = max(50, diff_max * 0.25)

print(
    f"Locations with >{low_threshold:.0f} people more exposed: {(pop_diff > low_threshold).sum()}"
)
print(
    f"Locations with >{low_threshold:.0f} people less exposed: {(pop_diff < -low_threshold).sum()}"
)
print(
    f"Locations with >{high_threshold:.0f} people more exposed: {(pop_diff > high_threshold).sum()}"
)
print(
    f"Locations with >{high_threshold:.0f} people less exposed: {(pop_diff < -high_threshold).sum()}"
)

# Inundation depth statistics
print(f"\nInundation depth statistics (for flooded areas > {INUNDATION_THRESHOLD}m):")
print(f"Mean inundation depth: {merged['inun_depth'].mean():.2f}m")
print(f"Max inundation depth: {merged['inun_depth'].max():.2f}m")
print(f"Median inundation depth: {merged['inun_depth'].median():.2f}m")

# ===== EXPORT SUMMARY DATA TO CSV =====
print("\nExporting summary data to CSV...")

# Calculate total exposed population for each scenario
total_cf0 = cf0_pop.sum()
total_cf8 = cf8_pop.sum()
difference = total_cf0 - total_cf8
percentage_change = ((total_cf0 - total_cf8) / total_cf8 * 100) if total_cf8 > 0 else 0

# Create summary DataFrame
import pandas as pd

summary_data = pd.DataFrame(
    {
        "Scenario": [
            "Factual (CF0)",
            "Counterfactual (CF-8)",
            "Difference",
            "Percentage Change",
        ],
        "Population_Column": [
            POPULATION_COLUMN,
            POPULATION_COLUMN,
            POPULATION_COLUMN,
            POPULATION_COLUMN,
        ],
        "Value": [total_cf0, total_cf8, difference, percentage_change],
        "Unit": ["people", "people", "people", "%"],
        "Inundation_Threshold_m": [
            INUNDATION_THRESHOLD,
            INUNDATION_THRESHOLD,
            INUNDATION_THRESHOLD,
            INUNDATION_THRESHOLD,
        ],
    }
)

# Save to CSV
csv_output_file = (
    OUTPUT_DIR / f"population_{EVENT_NAME.lower()}_{POPULATION_COLUMN}_summary.csv"
)
summary_data.to_csv(csv_output_file, index=False)
print(f"  Saved summary CSV: {csv_output_file}")

print(f"\nAnalysis complete for {EVENT_NAME}! Check the saved files in: {OUTPUT_DIR}")
print("Files created:")
print(f"  - {output_file_main}")
print(f"  - {output_file_bar}")
print(f"  - {output_file_diff}")
print(f"  - {output_file_factual}")
print(f"  - {output_file_counterfactual}")
print(f"  - {csv_output_file}")
print(f"\nTotal locations analyzed: {len(merged)}")
print(
    f"Locations with exposed population in CF0: {(merged[f'{POPULATION_COLUMN}_cf0'] > 0).sum()}"
)
print(
    f"Locations with exposed population in CF-8: {(merged[f'{POPULATION_COLUMN}_cf8'] > 0).sum()}"
)

"""
Plot counterfactual vs factual damage changes - RASTERIZED VERSION

This script rasterizes building-level damage data to a regular grid (default 0.025°)
by summing all damage values from buildings within each grid cell.

This provides a cleaner visualization for large datasets and allows for
easier comparison with other raster datasets.

by @dumontgoulart
"""

import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.colors import BoundaryNorm, ListedColormap
import contextily as ctx
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

# ===== CONFIGURATION =====
# Set your event name here
EVENT_NAME = "Freddy"  # Change this to: "Kenneth", "Freddy", etc.

# Choose damage column to plot: "total_damage" or "relative_damage"
DAMAGE_COLUMN = "total_damage"  # Change this to "total_damage" if preferred

# Grid resolution in degrees
GRID_RESOLUTION = 0.025  # degrees (~2.5 km at equator)

# Base paths - update these as needed
BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs/test")
OUTPUT_DIR = Path("/p/11210471-001-compass/04_Results/CF_figs")

# ===== DYNAMIC FILE PATHS =====
# Construct file paths based on event name
file_cf0 = (
    BASE_RUN_PATH
    / EVENT_NAME
    / "fiat"
    / "event_tp_era5_hourly_CF0_GTSMv41opendap_CF0_no_wind_CF0"
    / "output"
    / "output_relative_damage.fgb"
)
file_cf8 = (
    BASE_RUN_PATH
    / EVENT_NAME
    / "fiat"
    / "event_tp_era5_hourly_CF-8_GTSMv41opendap_CF0_no_wind_CF0"
    / "output"
    / "output_relative_damage.fgb"
)

# Create output directory if it doesn't exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Processing economic damage (RASTERIZED) for event: {EVENT_NAME}")
print(f"Using damage column: {DAMAGE_COLUMN}")
print(f"Grid resolution: {GRID_RESOLUTION}°")


# ===== RASTERIZATION FUNCTION =====
def rasterize_damage_to_grid(
    gdf: gpd.GeoDataFrame,
    value_column: str,
    resolution: float = 0.025,
    bounds: tuple = None,
    agg_func: str = "sum",
) -> xr.DataArray:
    """
    Rasterize point/polygon damage data to a regular grid.

    Parameters
    ----------
    gdf : GeoDataFrame
        GeoDataFrame with damage data
    value_column : str
        Column name containing the values to aggregate
    resolution : float
        Grid resolution in degrees
    bounds : tuple
        (minx, miny, maxx, maxy) bounds for the grid. If None, uses data bounds.
    agg_func : str
        Aggregation function: 'sum', 'mean', 'max', 'count'

    Returns
    -------
    xr.DataArray
        Rasterized damage data
    """
    # Get centroids if geometry is polygon
    if gdf.geometry.iloc[0].geom_type != "Point":
        centroids = gdf.geometry.centroid
        x_coords = centroids.x.values
        y_coords = centroids.y.values
    else:
        x_coords = gdf.geometry.x.values
        y_coords = gdf.geometry.y.values

    values = gdf[value_column].values

    # Determine bounds
    if bounds is None:
        minx, miny, maxx, maxy = (
            x_coords.min(),
            y_coords.min(),
            x_coords.max(),
            y_coords.max(),
        )
    else:
        minx, miny, maxx, maxy = bounds

    # Add small buffer to include edge points
    minx -= resolution / 2
    maxx += resolution / 2
    miny -= resolution / 2
    maxy += resolution / 2

    # Create grid edges
    x_edges = np.arange(minx, maxx + resolution, resolution)
    y_edges = np.arange(miny, maxy + resolution, resolution)

    # Create grid centers
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2

    # Bin the data
    # Find which bin each point belongs to
    x_bins = np.digitize(x_coords, x_edges) - 1
    y_bins = np.digitize(y_coords, y_edges) - 1

    # Clip to valid range
    x_bins = np.clip(x_bins, 0, len(x_centers) - 1)
    y_bins = np.clip(y_bins, 0, len(y_centers) - 1)

    # Create output array
    grid = np.zeros((len(y_centers), len(x_centers)))
    count_grid = np.zeros((len(y_centers), len(x_centers)))

    # Aggregate values into grid cells
    for i in range(len(values)):
        if not np.isnan(values[i]):
            if agg_func == "sum":
                grid[y_bins[i], x_bins[i]] += values[i]
            elif agg_func == "max":
                grid[y_bins[i], x_bins[i]] = max(grid[y_bins[i], x_bins[i]], values[i])
            elif agg_func in ["mean", "count"]:
                grid[y_bins[i], x_bins[i]] += values[i]
            count_grid[y_bins[i], x_bins[i]] += 1

    # Calculate mean if requested
    if agg_func == "mean":
        with np.errstate(divide="ignore", invalid="ignore"):
            grid = np.where(count_grid > 0, grid / count_grid, np.nan)
    elif agg_func == "count":
        grid = count_grid

    # Set cells with no data to NaN
    if agg_func != "count":
        grid = np.where(count_grid > 0, grid, np.nan)

    # Create xarray DataArray
    da = xr.DataArray(
        data=grid,
        dims=["y", "x"],
        coords={"y": y_centers, "x": x_centers},
        attrs={
            "resolution": resolution,
            "units": "USD" if "damage" in value_column.lower() else "",
            "aggregation": agg_func,
        },
    )

    return da


# ===== LOAD DAMAGE DATA =====
print("Loading damage data files...")
print(f"Loading CF0: {file_cf0}")
print(f"Loading CF-8: {file_cf8}")

# Read the geodataframes
gdf_cf0 = gpd.read_file(file_cf0)
gdf_cf8 = gpd.read_file(file_cf8)

print(f"CF0 data shape: {gdf_cf0.shape}")
print(f"CF-8 data shape: {gdf_cf8.shape}")
print(f"CF0 columns: {list(gdf_cf0.columns)}")

# Check if the damage column exists
if DAMAGE_COLUMN not in gdf_cf0.columns:
    print(
        f"Warning: {DAMAGE_COLUMN} not found in CF0 data. Available columns: {list(gdf_cf0.columns)}"
    )
    if "relative_damage" in gdf_cf0.columns:
        DAMAGE_COLUMN = "relative_damage"
    elif "total_damage" in gdf_cf0.columns:
        DAMAGE_COLUMN = "total_damage"
    else:
        raise ValueError("No suitable damage column found!")

print(f"Final damage column used: {DAMAGE_COLUMN}")

# ===== CONVERT TO LAT/LON =====
print("Converting data to lat/lon coordinates (EPSG:4326)...")
print(f"Original CRS: {gdf_cf0.crs}")

if gdf_cf0.crs != "EPSG:4326":
    gdf_cf0 = gdf_cf0.to_crs("EPSG:4326")
    gdf_cf8 = gdf_cf8.to_crs("EPSG:4326")
    print("Converted to EPSG:4326")
else:
    print("Already in EPSG:4326")

# ===== DETERMINE COMMON BOUNDS =====
print("Determining common bounds for rasterization...")
bounds_cf0 = gdf_cf0.total_bounds
bounds_cf8 = gdf_cf8.total_bounds

# Use the union of both bounds
common_bounds = (
    min(bounds_cf0[0], bounds_cf8[0]),  # minx
    min(bounds_cf0[1], bounds_cf8[1]),  # miny
    max(bounds_cf0[2], bounds_cf8[2]),  # maxx
    max(bounds_cf0[3], bounds_cf8[3]),  # maxy
)
print(f"Common bounds: {common_bounds}")

# ===== RASTERIZE DAMAGE DATA =====
print(f"Rasterizing damage data to {GRID_RESOLUTION}° grid...")

# Rasterize CF0 (Factual)
print("  Rasterizing CF0 (Factual)...")
da_cf0 = rasterize_damage_to_grid(
    gdf_cf0,
    value_column=DAMAGE_COLUMN,
    resolution=GRID_RESOLUTION,
    bounds=common_bounds,
    agg_func="sum",
)
print(f"  CF0 raster shape: {da_cf0.shape}")

# Rasterize CF8 (Counterfactual)
print("  Rasterizing CF-8 (Counterfactual)...")
da_cf8 = rasterize_damage_to_grid(
    gdf_cf8,
    value_column=DAMAGE_COLUMN,
    resolution=GRID_RESOLUTION,
    bounds=common_bounds,
    agg_func="sum",
)
print(f"  CF-8 raster shape: {da_cf8.shape}")

# ===== CALCULATE DIFFERENCE =====
print("Calculating difference raster...")

# Handle NaN values: treat as 0 for difference calculation
da_cf0_filled = da_cf0.fillna(0)
da_cf8_filled = da_cf8.fillna(0)

# Calculate difference (Factual - Counterfactual)
da_diff = da_cf0_filled - da_cf8_filled

# Mask difference where both are 0/NaN
mask_any_damage = (da_cf0_filled > 0) | (da_cf8_filled > 0)
da_diff = da_diff.where(mask_any_damage)

print(f"Difference raster shape: {da_diff.shape}")

# ===== STATISTICS =====
print(f"\nRasterized Damage Statistics for {EVENT_NAME} ({DAMAGE_COLUMN}):")

total_cf0 = float(da_cf0_filled.sum())
total_cf8 = float(da_cf8_filled.sum())
total_diff = total_cf0 - total_cf8

if DAMAGE_COLUMN == "total_damage":
    print(f"CF0 (Factual) total damage: ${total_cf0:,.0f}")
    print(f"CF-8 (Counterfactual) total damage: ${total_cf8:,.0f}")
    print(f"Difference (climate attribution): ${total_diff:,.0f}")
    if total_cf8 > 0:
        print(f"Percentage change: {(total_diff / total_cf8 * 100):.1f}%")
else:
    print(f"CF0 (Factual) total: {total_cf0:.2f}")
    print(f"CF-8 (Counterfactual) total: {total_cf8:.2f}")
    print(f"Difference: {total_diff:.2f}")

print(f"Grid cells with CF0 damage: {int((da_cf0_filled > 0).sum())}")
print(f"Grid cells with CF-8 damage: {int((da_cf8_filled > 0).sum())}")
print(f"Grid cells with difference: {int(mask_any_damage.sum())}")

# ===== PLOTTING SETUP =====
print("\nSetting up plots...")

data_crs = "EPSG:4326"

# Define colormaps and levels
if DAMAGE_COLUMN == "total_damage":
    # Use quantile-based boundaries for better visualization
    max_val = max(
        float(da_cf0.quantile(0.95, skipna=True)),
        float(da_cf8.quantile(0.95, skipna=True)),
    )
    if max_val > 0:
        boundaries = np.linspace(0, max_val, 11)
    else:
        boundaries = np.linspace(0, 1, 11)
    damage_norm = BoundaryNorm(boundaries, ncolors=256, clip=True)
    damage_cmap = plt.get_cmap("Reds")
    damage_label = f"Total Damage per {GRID_RESOLUTION}° cell [USD]"
else:
    boundaries = np.linspace(0, 100, 11)
    damage_norm = BoundaryNorm(boundaries, ncolors=256, clip=True)
    damage_cmap = plt.get_cmap("Reds")
    damage_label = f"Relative Damage per {GRID_RESOLUTION}° cell [%]"

# Difference colormap
diff_max = float(da_diff.quantile(0.9, skipna=True))
diff_min = float(da_diff.quantile(0.1, skipna=True))
diff_abs_max = max(abs(diff_max), abs(diff_min))
if diff_abs_max > 0:
    diff_boundaries = np.linspace(-diff_abs_max, diff_abs_max, 11)
else:
    diff_boundaries = np.linspace(-1, 1, 11)
diff_norm = BoundaryNorm(diff_boundaries, ncolors=256, clip=True)
diff_cmap = plt.get_cmap("RdBu_r")

print(f"Damage boundaries: {boundaries}")
print(f"Difference boundaries: {diff_boundaries}")

# ===== CREATE MAIN COMPARISON PLOT (3-PANEL) =====
print("Creating rasterized damage comparison plots...")

fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# Plot CF0 (Factual)
im1 = da_cf0.plot(
    ax=axes[0],
    cmap=damage_cmap,
    norm=damage_norm,
    add_colorbar=False,
    x="x",
    y="y",
    alpha=0.8,
    zorder=2,
)
axes[0].set_aspect("equal")
axes[0].set_title("Factual", fontsize=12, fontweight="bold")

# Plot CF-8 (Counterfactual)
im2 = da_cf8.plot(
    ax=axes[1],
    cmap=damage_cmap,
    norm=damage_norm,
    add_colorbar=False,
    x="x",
    y="y",
    alpha=0.8,
    zorder=2,
)
axes[1].set_aspect("equal")
axes[1].set_title("Counterfactual", fontsize=12, fontweight="bold")

# Plot Difference
im3 = da_diff.plot(
    ax=axes[2],
    cmap=diff_cmap,
    norm=diff_norm,
    add_colorbar=False,
    x="x",
    y="y",
    alpha=0.8,
    zorder=2,
)
axes[2].set_aspect("equal")
axes[2].set_title(
    "Damage Changes: Factual vs Counterfactual", fontsize=12, fontweight="bold"
)

# Add basemap to all axes
for ax in axes:
    try:
        ctx.add_basemap(
            ax=ax,
            source=ctx.providers.OpenStreetMap.Mapnik,
            crs=data_crs,
            attribution=False,
            zorder=1,
            zoom=12,
        )
    except Exception as e:
        print(f"Could not add basemap: {e}")

# Format axis labels
for ax in axes:
    ax.set_xlabel("Longitude [°]", fontsize=10)
    ax.set_ylabel("Latitude [°]", fontsize=10)

# Add colorbars
cbar1 = plt.colorbar(im1, ax=axes[0], shrink=0.8, pad=0.1)
cbar1.set_label(damage_label, fontsize=10)

cbar2 = plt.colorbar(im2, ax=axes[1], shrink=0.8, pad=0.1)
cbar2.set_label(damage_label, fontsize=10)

cbar3 = plt.colorbar(im3, ax=axes[2], shrink=0.8, pad=0.1)
cbar3.set_label(f"{damage_label} Difference", fontsize=10)

plt.tight_layout()
plt.suptitle(
    f"Economic Damage Analysis (Rasterized {GRID_RESOLUTION}°) - {EVENT_NAME}: Factual vs Counterfactual",
    fontsize=16,
    fontweight="bold",
    y=1.02,
)

# Save the figure
output_file_main = (
    OUTPUT_DIR
    / f"damage_{EVENT_NAME.lower()}_{DAMAGE_COLUMN}_raster_{GRID_RESOLUTION}deg_comparison.png"
)
plt.savefig(output_file_main, dpi=300, bbox_inches="tight")
plt.show()

# ===== CREATE BAR CHART (TOTAL DAMAGE) =====
if DAMAGE_COLUMN == "total_damage":
    print("Creating total damage bar chart...")

    fig_bar, ax_bar = plt.subplots(1, 1, figsize=(6, 6))

    scenarios = ["Counterfactual", "Factual"]
    totals = [total_cf8, total_cf0]
    bars = ax_bar.bar(
        scenarios, totals, color=["lightcoral", "lightcoral"], alpha=1, width=0.4
    )

    ax_bar.set_ylabel("Total Damage [USD]", fontsize=12, fontweight="bold")
    ax_bar.set_title(f"Total Damage - {EVENT_NAME}", fontsize=14, fontweight="bold")
    ax_bar.set_axisbelow(True)

    # Calculate difference for annotation
    difference = abs(total_cf0 - total_cf8)
    percentage_diff = (
        (difference / min(total_cf0, total_cf8)) * 100
        if min(total_cf0, total_cf8) > 0
        else 0
    )

    # Add annotation
    bar_positions = [bar.get_x() + bar.get_width() / 2 for bar in bars]
    bar_heights = [bar.get_height() for bar in bars]
    higher_bar_idx = 0 if totals[0] > totals[1] else 1
    lower_bar_idx = 1 - higher_bar_idx

    line_x_position = bar_positions[1] + (bar_positions[1] - bar_positions[0]) * 0.3
    line_extension = (bar_positions[1] - bar_positions[0]) * 0.04

    # Draw vertical difference line
    ax_bar.plot(
        [line_x_position, line_x_position],
        [bar_heights[lower_bar_idx], bar_heights[higher_bar_idx]],
        "k-",
        linewidth=2,
    )

    # Dashed horizontal lines
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

    # Ticks
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

    # Text annotation
    text_x = line_x_position + (bar_positions[1] - bar_positions[0]) * 0.1
    text_y = ((bar_heights[lower_bar_idx] + bar_heights[higher_bar_idx]) / 2) * 0.97

    if difference >= 1e9:
        diff_text = f"${difference/1e9:.1f}B ({percentage_diff:.0f}%)"
    elif difference >= 1e6:
        diff_text = f"${difference/1e6:.1f}M ({percentage_diff:.0f}%)"
    elif difference >= 1e3:
        diff_text = f"${difference/1e3:.1f}K ({percentage_diff:.0f}%)"
    else:
        diff_text = f"${difference:.0f} ({percentage_diff:.0f}%)"

    ax_bar.text(text_x, text_y, diff_text, ha="left", va="center", fontsize=10)
    text_offset = (bar_heights[higher_bar_idx] - bar_heights[lower_bar_idx]) * 0.35
    ax_bar.text(
        text_x,
        text_y + text_offset,
        "Climate change \n attribution:",
        ha="left",
        va="center",
        fontsize=9,
        style="italic",
    )

    ax_bar.set_ylim(0, max(totals) * 1.1)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    plt.tight_layout()

    output_file_bar = (
        OUTPUT_DIR
        / f"damage_{EVENT_NAME.lower()}_{DAMAGE_COLUMN}_raster_totals_barchart.png"
    )
    plt.savefig(output_file_bar, dpi=300, bbox_inches="tight")
    plt.show()
else:
    output_file_bar = None

# ===== DETAILED DIFFERENCE PLOT =====
print("Creating detailed difference plot...")

fig2, ax = plt.subplots(1, 1, figsize=(12, 8))

im_diff = da_diff.plot(
    ax=ax,
    cmap=diff_cmap,
    norm=diff_norm,
    add_colorbar=False,
    x="x",
    y="y",
    alpha=0.8,
    zorder=2,
)
ax.set_aspect("equal")

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

ax.set_xlabel("Longitude [°]", fontsize=12)
ax.set_ylabel("Latitude [°]", fontsize=12)

cbar = plt.colorbar(im_diff, ax=ax, shrink=0.8, pad=0.1)
cbar.set_label(f"{damage_label} Difference", fontsize=12)

ax.set_title(
    f"Economic Damage Changes (Rasterized {GRID_RESOLUTION}°) - {EVENT_NAME}\n(Red = Higher damage, Blue = Lower damage)",
    fontsize=14,
    fontweight="bold",
    pad=20,
)

output_file_diff = (
    OUTPUT_DIR
    / f"damage_{EVENT_NAME.lower()}_{DAMAGE_COLUMN}_raster_{GRID_RESOLUTION}deg_difference_only.png"
)
plt.savefig(output_file_diff, dpi=300, bbox_inches="tight")
plt.show()

# ===== INDIVIDUAL FACTUAL PLOT =====
print("Creating factual-only plot...")

fig3, ax_f = plt.subplots(1, 1, figsize=(12, 8))

im_factual = da_cf0.plot(
    ax=ax_f,
    cmap=damage_cmap,
    norm=damage_norm,
    add_colorbar=False,
    x="x",
    y="y",
    alpha=0.8,
    zorder=2,
)
ax_f.set_aspect("equal")

try:
    ctx.add_basemap(
        ax=ax_f,
        source=ctx.providers.OpenStreetMap.Mapnik,
        crs=data_crs,
        attribution=False,
        zorder=1,
        zoom=14,
    )
except Exception as e:
    print(f"Could not add basemap: {e}")

ax_f.set_xlabel("Longitude [°]", fontsize=12)
ax_f.set_ylabel("Latitude [°]", fontsize=12)

cbar_f = plt.colorbar(im_factual, ax=ax_f, shrink=0.8, pad=0.1)
cbar_f.set_label(damage_label, fontsize=12)

ax_f.set_title(
    f"Economic Damage (Rasterized {GRID_RESOLUTION}°) - {EVENT_NAME}: Factual",
    fontsize=14,
    fontweight="bold",
    pad=20,
)

output_file_factual = (
    OUTPUT_DIR
    / f"damage_{EVENT_NAME.lower()}_{DAMAGE_COLUMN}_raster_{GRID_RESOLUTION}deg_factual_only.png"
)
plt.savefig(output_file_factual, dpi=300, bbox_inches="tight")
plt.show()

# ===== INDIVIDUAL COUNTERFACTUAL PLOT =====
print("Creating counterfactual-only plot...")

fig4, ax_cf = plt.subplots(1, 1, figsize=(12, 8))

im_cf = da_cf8.plot(
    ax=ax_cf,
    cmap=damage_cmap,
    norm=damage_norm,
    add_colorbar=False,
    x="x",
    y="y",
    alpha=0.8,
    zorder=2,
)
ax_cf.set_aspect("equal")

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

ax_cf.set_xlabel("Longitude [°]", fontsize=12)
ax_cf.set_ylabel("Latitude [°]", fontsize=12)

cbar_cf = plt.colorbar(im_cf, ax=ax_cf, shrink=0.8, pad=0.1)
cbar_cf.set_label(damage_label, fontsize=12)

ax_cf.set_title(
    f"Economic Damage (Rasterized {GRID_RESOLUTION}°) - {EVENT_NAME}: Counterfactual",
    fontsize=14,
    fontweight="bold",
    pad=20,
)

output_file_counterfactual = (
    OUTPUT_DIR
    / f"damage_{EVENT_NAME.lower()}_{DAMAGE_COLUMN}_raster_{GRID_RESOLUTION}deg_counterfactual_only.png"
)
plt.savefig(output_file_counterfactual, dpi=300, bbox_inches="tight")
plt.show()

# ===== EXPORT RASTER DATA TO NETCDF =====
print("\nExporting rasterized data to NetCDF...")

# Create dataset with all rasters
ds_output = xr.Dataset(
    {
        "damage_factual": da_cf0,
        "damage_counterfactual": da_cf8,
        "damage_difference": da_diff,
    },
    attrs={
        "event": EVENT_NAME,
        "damage_column": DAMAGE_COLUMN,
        "grid_resolution_deg": GRID_RESOLUTION,
        "aggregation_method": "sum",
        "crs": "EPSG:4326",
        "units": "USD" if DAMAGE_COLUMN == "total_damage" else "percent",
    },
)

output_file_nc = (
    OUTPUT_DIR
    / f"damage_{EVENT_NAME.lower()}_{DAMAGE_COLUMN}_raster_{GRID_RESOLUTION}deg.nc"
)
ds_output.to_netcdf(output_file_nc)
print(f"  Saved NetCDF: {output_file_nc}")

# ===== EXPORT SUMMARY DATA TO CSV =====
print("Exporting summary data to CSV...")

summary_data = pd.DataFrame(
    {
        "Scenario": [
            "Factual (CF0)",
            "Counterfactual (CF-8)",
            "Difference",
            "Percentage Change",
        ],
        "Damage_Type": [DAMAGE_COLUMN] * 4,
        "Value": [
            total_cf0,
            total_cf8,
            total_diff,
            (total_diff / total_cf8 * 100) if total_cf8 > 0 else 0,
        ],
        "Unit": [
            "USD" if DAMAGE_COLUMN == "total_damage" else "%",
            "USD" if DAMAGE_COLUMN == "total_damage" else "%",
            "USD" if DAMAGE_COLUMN == "total_damage" else "%",
            "%",
        ],
        "Grid_Resolution_deg": [GRID_RESOLUTION] * 4,
    }
)

output_file_csv = (
    OUTPUT_DIR
    / f"damage_{EVENT_NAME.lower()}_{DAMAGE_COLUMN}_raster_{GRID_RESOLUTION}deg_summary.csv"
)
summary_data.to_csv(output_file_csv, index=False)
print(f"  Saved CSV: {output_file_csv}")

# ===== FINAL SUMMARY =====
print(f"\nAnalysis complete for {EVENT_NAME}! Check the saved files in: {OUTPUT_DIR}")
print("Files created:")
print(f"  - {output_file_main}")
if output_file_bar:
    print(f"  - {output_file_bar}")
print(f"  - {output_file_diff}")
print(f"  - {output_file_factual}")
print(f"  - {output_file_counterfactual}")
print(f"  - {output_file_nc}")
print(f"  - {output_file_csv}")
print(f"\nGrid statistics:")
print(f"  - Resolution: {GRID_RESOLUTION}° (~{GRID_RESOLUTION * 111:.1f} km)")
print(f"  - Grid size: {da_cf0.shape[1]} x {da_cf0.shape[0]} cells")
print(f"  - Total cells with data: {int(mask_any_damage.sum())}")

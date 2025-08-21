import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import BoundaryNorm, ListedColormap
import rioxarray as rxr  # Required for reading TIFF files
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

# File paths
# TODO: this is a preliminary script that is still not fully embedded in the snakemake workflow. The user is required to change the events & input files for now.
# For idai runs, I tested: CF: event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10; F: event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0

# ===== CONFIGURATION =====
# Set your event name here
EVENT_NAME = "Kenneth"  # Change this to: "Kenneth", "Freddy", etc.

# Base paths - update these as needed
BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs/test")
OUTPUT_DIR = Path("/p/11210471-001-compass/04_Results/CF_figs")

# ===== DYNAMIC FILE PATHS =====
# Construct file paths based on event name
file_cf0 = (
    BASE_RUN_PATH
    / EVENT_NAME
    / "sfincs"
    / "event_tp_era5_hourly_zarr_CF0_GTSMv41opendap_CF0_no_wind_CF0"
    / "plot_output"
    / "sfincs_output_hmax_AllTime.tif"
)
file_cf8 = (
    BASE_RUN_PATH
    / EVENT_NAME
    / "sfincs"
    / "event_tp_era5_hourly_zarr_CF-8_GTSMv41opendap_CF0_no_wind_CF0"
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

# ===== HANDLE NAN VALUES FOR DIFFERENCE CALCULATION =====
print("Handling NaN values for difference calculation...")
# Create masks for where each dataset has valid values
mask_cf0_valid = ~np.isnan(zsmax_cf0)
mask_cf8_valid = ~np.isnan(zsmax_cf8)

# For CF0: where CF0 is NaN but CF8 has a value, set CF0 to 0
zsmax_cf0 = zsmax_cf0.where(mask_cf0_valid | ~mask_cf8_valid, 0)
# For CF8: where CF8 is NaN but CF0 has a value, set CF8 to 0
zsmax_cf8 = zsmax_cf8.where(mask_cf8_valid | ~mask_cf0_valid, 0)

# ===== CALCULATE DIFFERENCE =====
print("Calculating differences...")
diff = zsmax_cf0 - zsmax_cf8  # Difference in maximum flood depth (CF0 - CF-8)

# ---- Percent difference (relative change) ----
# (Factual - Counterfactual) / Counterfactual * 100
# Mask where counterfactual depth is zero or NaN to avoid misleading huge %.
percent_diff = xr.where(zsmax_cf8 > 0, (diff / zsmax_cf8) * 100.0, np.nan)

# Flood extent (area) change (using same threshold as plots)
flood_threshold = 0.05  # meters
flood_mask_cf0 = zsmax_cf0 > flood_threshold
flood_mask_cf8 = zsmax_cf8 > flood_threshold
try:
    dx = float(np.abs(zsmax_cf0.x.diff("x").median()))
    dy = float(np.abs(zsmax_cf0.y.diff("y").median()))
    cell_area_m2 = dx * dy
except Exception:
    cell_area_m2 = np.nan
area_cf0_m2 = flood_mask_cf0.sum().values * cell_area_m2
area_cf8_m2 = flood_mask_cf8.sum().values * cell_area_m2
extent_area_diff_m2 = area_cf0_m2 - area_cf8_m2
extent_area_pct = (
    (extent_area_diff_m2 / area_cf8_m2 * 100.0) if area_cf8_m2 > 0 else np.nan
)

# Flood volume (sum depth * cell area over flooded cells)
if not np.isnan(cell_area_m2):
    volume_cf0_m3 = float(
        zsmax_cf0.where(flood_mask_cf0).sum(skipna=True) * cell_area_m2
    )
    volume_cf8_m3 = float(
        zsmax_cf8.where(flood_mask_cf8).sum(skipna=True) * cell_area_m2
    )
    volume_diff_m3 = volume_cf0_m3 - volume_cf8_m3
    volume_pct = (
        (volume_diff_m3 / volume_cf8_m3 * 100.0) if volume_cf8_m3 > 0 else np.nan
    )
else:
    volume_cf0_m3 = volume_cf8_m3 = volume_diff_m3 = volume_pct = np.nan

# Only show areas with positive water depth (above ground)
zsmax_cf0_plot = zsmax_cf0.where(zsmax_cf0 > 0.05)
zsmax_cf8_plot = zsmax_cf8.where(zsmax_cf8 > 0.05)

# Optional: only show differences where there's significant flooding
# diff = diff.where((zsmax_cf0 > 0.01) | (zsmax_cf8 > 0.01))

# ===== PRINT STATISTICS (REDUCED) =====
# Remove previous verbose block
# print(f"\nStatistics for {EVENT_NAME}:")
# ...previous detailed prints removed...

# Mean differences restricted to flooded cells (using same threshold)
flood_mask_union = (zsmax_cf0 > flood_threshold) | (zsmax_cf8 > flood_threshold)
flood_mask_intersection = (zsmax_cf0 > flood_threshold) & (zsmax_cf8 > flood_threshold)
mean_diff_flooded_union = diff.where(flood_mask_union).mean(skipna=True)
# (intersection retained if needed later, but not printed)
# mean_diff_flooded_intersection = diff.where(flood_mask_intersection).mean(skipna=True)

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
try:
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(20, 6),
        subplot_kw={
            "projection": ccrs.UTM(zone=utm_zone, southern_hemisphere=southern)
        },
    )
    use_cartopy = True
except Exception as e:
    print(f"Cartopy projection failed: {e}")
    print("Falling back to simple plots without projection...")
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    use_cartopy = False

# Define colormaps and levels
flood_levels = np.arange(0, 4, 0.2)
flood_cmap = plt.cm.Blues

# Difference levels (symmetric around 0)
dif_min = -0.5
dif_max = 0.5
diff_levels = np.linspace(dif_min, dif_max, 21)
diff_cmap = plt.cm.RdBu  # Red for increases, Blue for decreases

# Plot CF0
ax1 = axes[0]
if use_cartopy:
    im1 = zsmax_cf0_plot.plot(
        ax=ax1,
        levels=flood_levels,
        cmap=flood_cmap,
        add_colorbar=False,
        transform=ccrs.UTM(zone=utm_zone, southern_hemisphere=southern),
        x="x",
        y="y",
    )
    ax1.coastlines(resolution="10m")
else:
    im1 = zsmax_cf0_plot.plot(
        ax=ax1, levels=flood_levels, cmap=flood_cmap, add_colorbar=False, x="x", y="y"
    )
ax1.set_title("Factual", fontsize=12, fontweight="bold")

# Plot CF-8
ax2 = axes[1]
if use_cartopy:
    im2 = zsmax_cf8_plot.plot(
        ax=ax2,
        levels=flood_levels,
        cmap=flood_cmap,
        add_colorbar=False,
        transform=ccrs.UTM(zone=utm_zone, southern_hemisphere=southern),
        x="x",
        y="y",
    )
    ax2.coastlines(resolution="10m")
else:
    im2 = zsmax_cf8_plot.plot(
        ax=ax2, levels=flood_levels, cmap=flood_cmap, add_colorbar=False, x="x", y="y"
    )
ax2.set_title("Counterfactual", fontsize=12, fontweight="bold")

# Plot difference
ax3 = axes[2]
if use_cartopy:
    im3 = diff.plot(
        ax=ax3,
        levels=diff_levels,
        cmap=diff_cmap,
        add_colorbar=False,
        transform=ccrs.UTM(zone=utm_zone, southern_hemisphere=southern),
        x="x",
        y="y",
    )
    ax3.coastlines(resolution="10m")
else:
    im3 = diff.plot(
        ax=ax3, levels=diff_levels, cmap=diff_cmap, add_colorbar=False, x="x", y="y"
    )
ax3.set_title(
    "Flood Depth Changes: Factual vs Counterfactual", fontsize=12, fontweight="bold"
)

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

# ===== DETAILED DIFFERENCE PLOT =====
print("Creating detailed difference plot...")
try:
    fig2, ax = plt.subplots(
        1,
        1,
        figsize=(12, 8),
        subplot_kw={
            "projection": ccrs.UTM(zone=utm_zone, southern_hemisphere=southern)
        },
    )
    use_cartopy2 = True
except:
    fig2, ax = plt.subplots(1, 1, figsize=(12, 8))
    use_cartopy2 = False

# Plot only the difference with better resolution
if use_cartopy2:
    im_diff = diff.plot(
        ax=ax,
        levels=diff_levels,
        cmap="RdBu",
        add_colorbar=False,
        transform=ccrs.UTM(zone=utm_zone, southern_hemisphere=southern),
        x="x",
        y="y",
    )
    ax.coastlines(resolution="10m")
else:
    im_diff = diff.plot(
        ax=ax, levels=diff_levels, cmap="RdBu", add_colorbar=False, x="x", y="y"
    )

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

# Close datasets
zsmax_cf0.close()
zsmax_cf8.close()

print(f"Analysis complete. Output figures:")
print(f"  - {output_file_main}")
print(f"  - {output_file_diff}")
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
print(f"  - {output_file_diff}")

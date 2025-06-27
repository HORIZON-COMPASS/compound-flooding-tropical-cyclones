import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import BoundaryNorm, ListedColormap
import warnings
warnings.filterwarnings('ignore')

# File paths
from pathlib import Path

# ===== CONFIGURATION =====
# Set your event name here
EVENT_NAME = "Kenneth"  # Change this to: "Kenneth", "Freddy", etc.

# Choose damage column to plot: "total_damage" or "relative_damage"
DAMAGE_COLUMN = "total_damage"  # Change this to "total_damage" if preferred

# Base paths - update these as needed
BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs/test")
OUTPUT_DIR = Path("/p/11210471-001-compass/04_Results/CF_figs")

# ===== DYNAMIC FILE PATHS =====
# Construct file paths based on event name
file_cf0 = BASE_RUN_PATH / EVENT_NAME / "fiat" / "event_tp_era5_hourly_zarr_CF0_GTSMv41opendap_CF0_no_wind_CF0" / "output" / "output_relative_damage.fgb"
file_cf8 = BASE_RUN_PATH / EVENT_NAME / "fiat" / "event_tp_era5_hourly_zarr_CF-8_GTSMv41opendap_CF0_no_wind_CF0" / "output" / "output_relative_damage.fgb"

# Create output directory if it doesn't exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Processing economic damage for event: {EVENT_NAME}")
print(f"Using damage column: {DAMAGE_COLUMN}")

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
print(f"CF-8 columns: {list(gdf_cf8.columns)}")

# Check if the damage column exists
if DAMAGE_COLUMN not in gdf_cf0.columns:
    print(f"Warning: {DAMAGE_COLUMN} not found in CF0 data. Available columns: {list(gdf_cf0.columns)}")
    # Try alternative column names
    if 'relative_damage' in gdf_cf0.columns:
        DAMAGE_COLUMN = 'relative_damage'
    elif 'total_damage' in gdf_cf0.columns:
        DAMAGE_COLUMN = 'total_damage'
    else:
        raise ValueError("No suitable damage column found!")

print(f"Final damage column used: {DAMAGE_COLUMN}")

# ===== EXTRACT COORDINATES =====
print("Extracting coordinates from building centroids...")
# Extract x, y coordinates from geometry centroids (for polygon buildings)
gdf_cf0['centroid'] = gdf_cf0.geometry.centroid
gdf_cf0['x'] = gdf_cf0['centroid'].x
gdf_cf0['y'] = gdf_cf0['centroid'].y

gdf_cf8['centroid'] = gdf_cf8.geometry.centroid
gdf_cf8['x'] = gdf_cf8['centroid'].x
gdf_cf8['y'] = gdf_cf8['centroid'].y

print(f"Sample coordinates CF0: x={gdf_cf0['x'].iloc[0]:.2f}, y={gdf_cf0['y'].iloc[0]:.2f}")
print(f"Sample coordinates CF-8: x={gdf_cf8['x'].iloc[0]:.2f}, y={gdf_cf8['y'].iloc[0]:.2f}")

# ===== MERGE DATA FOR DIFFERENCE CALCULATION =====
print("Merging data for difference calculation...")
# Merge on object_id to compare same buildings
# Use coordinates from CF0 data, but handle cases where buildings might only exist in one scenario
merged = gdf_cf0[['object_id', 'x', 'y', DAMAGE_COLUMN]].merge(
    gdf_cf8[['object_id', 'x', 'y', DAMAGE_COLUMN]], 
    on='object_id', 
    how='outer', 
    suffixes=('_cf0', '_cf8')
)

# Handle coordinates: use CF0 coordinates where available, otherwise CF8
merged['x'] = merged['x_cf0'].fillna(merged['x_cf8'])
merged['y'] = merged['y_cf0'].fillna(merged['y_cf8'])

# Drop the temporary coordinate columns
merged = merged.drop(['x_cf0', 'y_cf0', 'x_cf8', 'y_cf8'], axis=1)

# Fill NaN values with 0 for buildings that don't exist in one scenario
merged[f'{DAMAGE_COLUMN}_cf0'] = merged[f'{DAMAGE_COLUMN}_cf0'].fillna(0)
merged[f'{DAMAGE_COLUMN}_cf8'] = merged[f'{DAMAGE_COLUMN}_cf8'].fillna(0)

# Convert relative damage from 0-1 scale to 0-100% scale if needed
if DAMAGE_COLUMN == 'relative_damage':
    print("Converting relative damage from 0-1 scale to 0-100% scale...")
    merged[f'{DAMAGE_COLUMN}_cf0'] = merged[f'{DAMAGE_COLUMN}_cf0'] * 100
    merged[f'{DAMAGE_COLUMN}_cf8'] = merged[f'{DAMAGE_COLUMN}_cf8'] * 100

# Calculate difference (CF0 - CF8)
merged['damage_diff'] = merged[f'{DAMAGE_COLUMN}_cf0'] - merged[f'{DAMAGE_COLUMN}_cf8']

# Remove buildings with no coordinates
merged = merged.dropna(subset=['x', 'y'])

print(f"Merged data shape: {merged.shape}")
print(f"Sample merged data (after scaling):")
print(merged[['object_id', 'x', 'y', f'{DAMAGE_COLUMN}_cf0', f'{DAMAGE_COLUMN}_cf8', 'damage_diff']].head())

# ===== STATISTICS =====
print(f"\nDamage Statistics for {EVENT_NAME} ({DAMAGE_COLUMN}):")
cf0_damage = merged[f'{DAMAGE_COLUMN}_cf0']
cf8_damage = merged[f'{DAMAGE_COLUMN}_cf8']
damage_diff = merged['damage_diff']

if DAMAGE_COLUMN == 'relative_damage':
    print(f"CF0 max damage: {cf0_damage.max():.1f}%")
    print(f"CF0 mean damage: {cf0_damage.mean():.1f}%")
    print(f"CF0 buildings with >0% damage: {(cf0_damage > 0).sum()}")
    print(f"CF-8 max damage: {cf8_damage.max():.1f}%")
    print(f"CF-8 mean damage: {cf8_damage.mean():.1f}%")
    print(f"CF-8 buildings with >0% damage: {(cf8_damage > 0).sum()}")
    print(f"Max increase (CF0 vs CF-8): {damage_diff.max():.1f}%")
    print(f"Max decrease (CF0 vs CF-8): {damage_diff.min():.1f}%")
    print(f"Mean difference: {damage_diff.mean():.1f}%")
else:
    print(f"CF0 max damage: {cf0_damage.max():.3f}")
    print(f"CF0 mean damage: {cf0_damage.mean():.3f}")
    print(f"CF0 total damage: {cf0_damage.sum():.0f}")
    print(f"CF-8 max damage: {cf8_damage.max():.3f}")
    print(f"CF-8 mean damage: {cf8_damage.mean():.3f}")
    print(f"CF-8 total damage: {cf8_damage.sum():.0f}")
    print(f"Max increase (CF0 vs CF-8): {damage_diff.max():.3f}")
    print(f"Max decrease (CF0 vs CF-8): {damage_diff.min():.3f}")
    print(f"Mean difference: {damage_diff.mean():.3f}")

# ===== DETERMINE COORDINATE SYSTEM =====
try:
    # Get approximate center coordinates
    center_x = merged['x'].mean()
    center_y = merged['y'].mean()
    
    print(f"Center coordinates: {center_x:.0f}, {center_y:.0f}")
    
    # Determine UTM zone based on region (adjust for your area)
    utm_zone = 37  # Adjust based on your region
    southern = True  # Set to False if in northern hemisphere
    
    print(f"Using UTM zone {utm_zone}, Southern: {southern}")
except Exception as e:
    print(f"Could not determine coordinates: {e}")
    utm_zone = 37
    southern = True

# ===== PLOTTING SETUP =====
print("Setting up plots...")

# Create projection for cartopy
try:
    crs = ccrs.UTM(zone=utm_zone, southern_hemisphere=southern)
    use_cartopy = True
except Exception as e:
    print(f"Cartopy projection failed: {e}")
    use_cartopy = False

# Define colormaps and levels with discrete intervals
if DAMAGE_COLUMN == 'relative_damage':
    # For relative damage (0-100%) with 10% intervals
    boundaries = np.linspace(0, 100, 11)  # Creates 10 intervals from 0 to 100
    damage_norm = BoundaryNorm(boundaries, ncolors=256, clip=True)
    damage_cmap = plt.get_cmap('Reds')  # Darker red colormap
    damage_label = 'Relative Damage [%]'
    vmin, vmax = 0, 100
elif DAMAGE_COLUMN == 'total_damage':
    # For total damage (currency) with appropriate intervals
    max_total = max(cf0_damage.quantile(0.9), cf8_damage.quantile(0.9))
    # Create 10 intervals from 0 to max
    boundaries = np.linspace(0, max_total, 11)
    damage_norm = BoundaryNorm(boundaries, ncolors=256, clip=True)
    damage_cmap = plt.get_cmap('Reds')
    damage_label = 'Total Damage [USD]'
    vmin, vmax = 0, max_total
else:
    # Generic fallback
    max_val = max(cf0_damage.max(), cf8_damage.max())
    boundaries = np.linspace(0, max_val, 11)
    damage_norm = BoundaryNorm(boundaries, ncolors=256, clip=True)
    damage_cmap = plt.get_cmap('Reds')
    damage_label = f'{DAMAGE_COLUMN}'
    vmin, vmax = 0, max_val

# Difference colormap with discrete intervals
diff_max = damage_diff.quantile(0.6)
if DAMAGE_COLUMN == 'relative_damage':
    # For relative damage differences, use smaller intervals
    diff_boundaries = np.linspace(-20, 20, 11)  # 5% intervals from -50% to +50%
else:
    # For absolute damage differences
    diff_boundaries = np.linspace(-diff_max, diff_max, 11)

diff_norm = BoundaryNorm(diff_boundaries, ncolors=256, clip=True)
diff_cmap = plt.get_cmap('RdBu_r')  # Red for increases, Blue for decreases

print(f"Damage boundaries: {boundaries}")
print(f"Difference boundaries: {diff_boundaries}")

# ===== CREATE MAIN COMPARISON PLOT =====
print("Creating damage comparison plots...")

if use_cartopy:
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), 
                            subplot_kw={'projection': crs})
else:
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# Plot settings
point_size = 17  # Increased point size for better visibility
alpha = 0.9      # Increased transparency for better colors

# Plot CF0 (Factual)
ax1 = axes[0]
if use_cartopy:
    scatter1 = ax1.scatter(merged['x'], merged['y'], c=merged[f'{DAMAGE_COLUMN}_cf0'], 
                          cmap=damage_cmap, norm=damage_norm, s=point_size, alpha=alpha, 
                          transform=crs)
    # ax1.coastlines(resolution='10m')
    # ax1.add_feature(cfeature.BORDERS)
    ax1.add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
else:
    scatter1 = ax1.scatter(merged['x'], merged['y'], c=merged[f'{DAMAGE_COLUMN}_cf0'], 
                          cmap=damage_cmap, norm=damage_norm, s=point_size, alpha=alpha)
ax1.set_title('Factual', fontsize=12, fontweight='bold')

# Plot CF-8 (Counterfactual)
ax2 = axes[1]
if use_cartopy:
    scatter2 = ax2.scatter(merged['x'], merged['y'], c=merged[f'{DAMAGE_COLUMN}_cf8'], 
                          cmap=damage_cmap, norm=damage_norm, s=point_size, alpha=alpha, 
                          transform=crs)
    # ax2.coastlines(resolution='10m')
    # ax2.add_feature(cfeature.BORDERS)
    ax2.add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
else:
    scatter2 = ax2.scatter(merged['x'], merged['y'], c=merged[f'{DAMAGE_COLUMN}_cf8'], 
                          cmap=damage_cmap, norm=damage_norm, s=point_size, alpha=alpha)
ax2.set_title('Counterfactual', fontsize=12, fontweight='bold')

# Plot Difference
ax3 = axes[2]
if use_cartopy:
    scatter3 = ax3.scatter(merged['x'], merged['y'], c=merged['damage_diff'], 
                          cmap=diff_cmap, norm=diff_norm, s=point_size, alpha=alpha, 
                          transform=crs)
    # ax3.coastlines(resolution='10m')
    # ax3.add_feature(cfeature.BORDERS)
    ax3.add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
else:
    scatter3 = ax3.scatter(merged['x'], merged['y'], c=merged['damage_diff'], 
                          cmap=diff_cmap, norm=diff_norm, s=point_size, alpha=alpha)
ax3.set_title('Damage Changes: Factual vs Counterfactual', fontsize=12, fontweight='bold')

# Add colorbars
cbar1 = plt.colorbar(scatter1, ax=ax1, shrink=0.8, pad=0.1)
cbar1.set_label(damage_label, fontsize=10)

cbar2 = plt.colorbar(scatter2, ax=ax2, shrink=0.8, pad=0.1)
cbar2.set_label(damage_label, fontsize=10)

cbar3 = plt.colorbar(scatter3, ax=ax3, shrink=0.8, pad=0.1)
cbar3.set_label(f'{damage_label} Difference', fontsize=10)

# Adjust layout and add main title
plt.tight_layout()
plt.suptitle(f'Economic Damage Analysis - {EVENT_NAME}: Factual vs Counterfactual ({DAMAGE_COLUMN})', 
             fontsize=16, fontweight='bold', y=1.02)

# Save the figure
output_file_main = OUTPUT_DIR / f'damage_{EVENT_NAME.lower()}_{DAMAGE_COLUMN}_comparison.png'
plt.savefig(output_file_main, dpi=300, bbox_inches='tight')
plt.show()

# ===== DETAILED DIFFERENCE PLOT =====
print("Creating detailed damage difference plot...")

if use_cartopy:
    fig2, ax = plt.subplots(1, 1, figsize=(12, 8), subplot_kw={'projection': crs})
else:
    fig2, ax = plt.subplots(1, 1, figsize=(12, 8))

# Plot only the difference with larger points
if use_cartopy:
    scatter_diff = ax.scatter(merged['x'], merged['y'], c=merged['damage_diff'], 
                             cmap=diff_cmap, norm=diff_norm, s=point_size*2, alpha=alpha, 
                             transform=crs)
    # ax.coastlines(resolution='10m')
    # ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
else:
    scatter_diff = ax.scatter(merged['x'], merged['y'], c=merged['damage_diff'], 
                             cmap=diff_cmap, norm=diff_norm, s=point_size*2, alpha=alpha)

# Colorbar
cbar = plt.colorbar(scatter_diff, ax=ax, shrink=0.8, pad=0.1)
cbar.set_label(f'{damage_label} Difference', fontsize=12)

# Title
ax.set_title(f'Economic Damage Changes - {EVENT_NAME}: Factual vs Counterfactual\n(Red = Higher damage, Blue = Lower damage)', 
             fontsize=14, fontweight='bold', pad=20)

# Save
output_file_diff = OUTPUT_DIR / f'damage_{EVENT_NAME.lower()}_{DAMAGE_COLUMN}_difference_only.png'
plt.savefig(output_file_diff, dpi=300, bbox_inches='tight')
plt.show()

# ===== SUMMARY STATISTICS =====
print(f"\nBuildings with significant damage changes:")
if DAMAGE_COLUMN == 'relative_damage':
    print(f"Buildings with >10% damage increase: {(damage_diff > 10).sum()}")
    print(f"Buildings with >10% damage decrease: {(damage_diff < -10).sum()}")
    print(f"Buildings with >25% damage increase: {(damage_diff > 25).sum()}")
    print(f"Buildings with >25% damage decrease: {(damage_diff < -25).sum()}")
    print(f"Buildings with >5% damage increase: {(damage_diff > 5).sum()}")
    print(f"Buildings with >5% damage decrease: {(damage_diff < -5).sum()}")
else:
    threshold_low = diff_max * 0.1
    threshold_high = diff_max * 0.25
    print(f"Buildings with >{threshold_low:.0f} damage increase: {(damage_diff > threshold_low).sum()}")
    print(f"Buildings with >{threshold_low:.0f} damage decrease: {(damage_diff < -threshold_low).sum()}")
    print(f"Buildings with >{threshold_high:.0f} damage increase: {(damage_diff > threshold_high).sum()}")
    print(f"Buildings with >{threshold_high:.0f} damage decrease: {(damage_diff < -threshold_high).sum()}")

print(f"\nAnalysis complete for {EVENT_NAME}! Check the saved PNG files in: {OUTPUT_DIR}")
print("Files created:")
print(f"  - {output_file_main}")
print(f"  - {output_file_diff}")
print(f"\nTotal buildings analyzed: {len(merged)}")
if DAMAGE_COLUMN == 'relative_damage':
    print(f"Buildings with damage in CF0: {(merged[f'{DAMAGE_COLUMN}_cf0'] > 0).sum()}")
    print(f"Buildings with damage in CF-8: {(merged[f'{DAMAGE_COLUMN}_cf8'] > 0).sum()}")
else:
    print(f"Buildings with damage in CF0: {(merged[f'{DAMAGE_COLUMN}_cf0'] > 0).sum()}")
    print(f"Buildings with damage in CF-8: {(merged[f'{DAMAGE_COLUMN}_cf8'] > 0).sum()}")
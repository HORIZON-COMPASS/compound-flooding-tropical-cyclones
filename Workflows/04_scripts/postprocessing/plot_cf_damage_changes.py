#%%
import geopandas as gpd
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import BoundaryNorm, ListedColormap
import warnings
warnings.filterwarnings('ignore')
import platform
# File paths
from pathlib import Path
# by @dumontgoulart
#TODO: this is a preliminary script that is still not fully embedded in the snakemake workflow. The user is required to change the events & input files for now.

# ===== CONFIGURATION =====
# Set your event name here
EVENT_NAME = "Idai"  # Change this to: "Kenneth", "Freddy", etc.

# Choose damage column to plot: "total_damage" or "relative_damage"
DAMAGE_COLUMN = "total_damage"  # Change this to "total_damage" if preferred

# Base paths - update these as needed
prefix = "p:/" if platform.system() == "Windows" else "/p/"

# BASE_RUN_PATH = Path(os.path.join(prefix,"11210471-001-compass","03_Runs","test"))
BASE_RUN_PATH = Path(os.path.join(prefix,"11210471-001-compass","03_Runs","sofala"))
OUTPUT_DIR = Path(os.path.join(prefix,"11210471-001-compass","04_Results","CF_figs", "redone"))

# ===== DYNAMIC FILE PATHS =====
# Construct file paths based on event name
# file_cf0 = BASE_RUN_PATH / EVENT_NAME / "fiat" / "event_tp_era5_hourly_zarr_CF0_GTSMv41opendap_CF0_no_wind_CF0" / "output" / "output_relative_damage.fgb"
# file_cf8 = BASE_RUN_PATH / EVENT_NAME / "fiat" / "event_tp_era5_hourly_zarr_CF-8_GTSMv41opendap_CF0_no_wind_CF0" / "output" / "output_relative_damage.fgb"

file_cf0 = BASE_RUN_PATH / EVENT_NAME / "fiat" / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" / "output" / "output_relative_damage.fgb"
file_cf8 = BASE_RUN_PATH / EVENT_NAME / "fiat" / "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10" / "output" / "output_relative_damage.fgb"

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

#%%
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


# ===== CORRECT CURRENCY =====
eur_to_usd = 1.326 # Convert JRC Damage Values (Euro 2010) into US-Dollars (2010)
usd_2010_to_2019 = 1.172 # Convert US-Dollars (2010) to US-Dollars (2019) - annual averages: 255.657 / 218.056
usd_2010_to_2023 = 1.397 # Convert US-Dollars (2010) to US-Dollars (2019) - annual averages: 304.702 / 218.056

if EVENT_NAME == "Idai" or EVENT_NAME == "Kenneth":
    gdf_cf0[DAMAGE_COLUMN] = gdf_cf0[DAMAGE_COLUMN] * eur_to_usd * usd_2010_to_2019
    gdf_cf8[DAMAGE_COLUMN] = gdf_cf8[DAMAGE_COLUMN] * eur_to_usd * usd_2010_to_2019
elif EVENT_NAME == "Freddy":
    gdf_cf0[DAMAGE_COLUMN] = gdf_cf0[DAMAGE_COLUMN] * eur_to_usd * usd_2010_to_2023
    gdf_cf8[DAMAGE_COLUMN] = gdf_cf8[DAMAGE_COLUMN] * eur_to_usd * usd_2010_to_2023
else:
    print("No currency conversion applied.")
    

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
    print(f"CF0 total damage: {cf0_damage.sum():.0f}%")
    print(f"CF0 buildings with >0% damage: {(cf0_damage > 0).sum()}")
    print(f"CF-8 max damage: {cf8_damage.max():.1f}%")
    print(f"CF-8 mean damage: {cf8_damage.mean():.1f}%")
    print(f"CF-8 total damage: {cf8_damage.sum():.0f}%")
    print(f"CF-8 buildings with >0% damage: {(cf8_damage > 0).sum()}")
    print(f"Difference in total damage: {cf0_damage.sum() - cf8_damage.sum():.0f}%")
    print(f"Percentage change: {((cf0_damage.sum() - cf8_damage.sum()) / cf8_damage.sum() * 100):.1f}%")
    print(f"Max increase (CF0 vs CF-8): {damage_diff.max():.1f}%")
    print(f"Max decrease (CF0 vs CF-8): {damage_diff.min():.1f}%")
    print(f"Mean difference: {damage_diff.mean():.1f}%")
else:
    print(f"CF0 max damage: ${cf0_damage.max():.0f}")
    print(f"CF0 mean damage: ${cf0_damage.mean():.0f}")
    print(f"CF0 total damage: ${cf0_damage.sum():.0f}")
    print(f"CF-8 max damage: ${cf8_damage.max():.0f}")
    print(f"CF-8 mean damage: ${cf8_damage.mean():.0f}")
    print(f"CF-8 total damage: ${cf8_damage.sum():.0f}")
    print(f"Difference in total damage: ${cf0_damage.sum() - cf8_damage.sum():.0f}")
    print(f"Percentage change: {((cf0_damage.sum() - cf8_damage.sum()) / cf8_damage.sum() * 100):.1f}%")
    print(f"Max increase (CF0 vs CF-8): ${damage_diff.max():.0f}")
    print(f"Max decrease (CF0 vs CF-8): ${damage_diff.min():.0f}")
    print(f"Mean difference: ${damage_diff.mean():.0f}")

#%%
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
use_cartopy=True
# Create projection for cartopy
try:
    crs = ccrs.UTM(zone=utm_zone, southern_hemisphere=southern)
    use_cartopy = True
except Exception as e:
    print(f"Cartopy projection failed: {e}")
    use_cartopy = False

use_cartopy=True

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
    if EVENT_NAME == "Idai" or EVENT_NAME == "Kenneth":
        damage_label = 'Total Damage [USD 2019]'
    elif EVENT_NAME == "Freddy":
        damage_label = 'Total Damage [USD 2023]'
    else:
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
diff_max = damage_diff.quantile(0.9)
diff_min = damage_diff.quantile(0.1)
if DAMAGE_COLUMN == 'relative_damage':
    # For relative damage differences, use smaller intervals
    diff_boundaries = np.linspace(-20, 20, 11)  # intervals from -20% to +20%
else:
    # For absolute damage differences
    diff_boundaries = np.linspace(diff_min, diff_max, 11)

diff_norm = BoundaryNorm(diff_boundaries, ncolors=256, clip=True)
diff_cmap = plt.get_cmap('RdBu_r')  # Red for increases, Blue for decreases

print(f"Damage boundaries: {boundaries}")
print(f"Difference boundaries: {diff_boundaries}")

# ===== CREATE MAIN COMPARISON PLOT (3-PANEL) =====
# print("Creating damage comparison plots...")
# Convert merged DataFrame to GeoDataFrame
gdf = gpd.GeoDataFrame(
    merged,
    geometry=gpd.points_from_xy(merged['x'], merged['y']),
    crs=f"EPSG:32736"  # UTM Zone 36S (adjust if needed)
)

# Convert to lat/lon
gdf = gdf.to_crs("EPSG:4326")

# Use these for plotting
merged['lon'] = gdf.geometry.x
merged['lat'] = gdf.geometry.y


crs = ccrs.PlateCarree()

if use_cartopy:
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), 
                            subplot_kw={'projection': crs})
else:
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# Plot settings
point_size = 20  # Point size for visibility
alpha = 0.8      # Transparency for better colors

# Plot CF0 (Factual)
if use_cartopy:
    scatter1 = axes[0].scatter(merged['lon'], merged['lat'], c=merged[f'{DAMAGE_COLUMN}_cf0'], 
                          cmap=damage_cmap, norm=damage_norm, s=point_size, alpha=alpha, 
                          transform=crs)
    axes[0].add_feature(cfeature.LAND, color='lightgray', alpha=0.5)

else:
    # axes[0].set_facecolor('lightgray')
    scatter1 = axes[0].scatter(merged['lon'], merged['lat'], c=merged[f'{DAMAGE_COLUMN}_cf0'], 
                          cmap=damage_cmap, norm=damage_norm, s=point_size, alpha=alpha)
axes[0].set_title('Factual', fontsize=12, fontweight='bold')

# Plot CF-8 (Counterfactual)
if use_cartopy:
    scatter2 = axes[1].scatter(merged['lon'], merged['lat'], c=merged[f'{DAMAGE_COLUMN}_cf8'], 
                          cmap=damage_cmap, norm=damage_norm, s=point_size, alpha=alpha, 
                          transform=crs)
    axes[1].add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
else:
    # axes[1].set_facecolor('lightgray')
    scatter2 = axes[1].scatter(merged['lon'], merged['lat'], c=merged[f'{DAMAGE_COLUMN}_cf8'], 
                          cmap=damage_cmap, norm=damage_norm, s=point_size, alpha=alpha)
axes[1].set_title('Counterfactual', fontsize=12, fontweight='bold')

# Plot Difference
# Remove invalid values
plot_data_diff = merged.copy()
plot_data_diff = plot_data_diff[
    np.isfinite(plot_data_diff['lon']) &
    np.isfinite(plot_data_diff['lat']) &
    np.isfinite(plot_data_diff['damage_diff']) &
    (plot_data_diff['lat'] >= -90) & (plot_data_diff['lat'] <= 90) &
    (plot_data_diff['lon'] >= -180) & (plot_data_diff['lon'] <= 180)
]

if use_cartopy:
    scatter3 = axes[2].scatter(plot_data_diff['lon'], plot_data_diff['lat'], c=plot_data_diff['damage_diff'], 
                          cmap=diff_cmap, norm=diff_norm, s=point_size, alpha=alpha, 
                          transform=crs)
    axes[2].add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
else:
    # axes[2].set_facecolor('lightgray')
    scatter3 = axes[2].scatter(plot_data_diff['lon'], plot_data_diff['lat'], c=plot_data_diff['damage_diff'], 
                          cmap=diff_cmap, norm=diff_norm, s=point_size, alpha=alpha)
axes[2].set_title('Damage Changes: Factual vs Counterfactual', fontsize=12, fontweight='bold')

# Add colorbars
cbar1 = plt.colorbar(scatter1, ax=axes[0], shrink=0.8, pad=0.1)
cbar1.set_label(damage_label, fontsize=10)

cbar2 = plt.colorbar(scatter2, ax=axes[1], shrink=0.8, pad=0.1)
cbar2.set_label(damage_label, fontsize=10)

cbar3 = plt.colorbar(scatter3, ax=axes[2], shrink=0.8, pad=0.1)
cbar3.set_label(f'{damage_label} Difference', fontsize=10)

# for ax in axes:
#     ax.add_feature(cfeature.OCEAN, color='white')  # or lightblue

# Adjust layout and add main title
plt.tight_layout()
# plt.suptitle(f'Economic Damage Analysis - {EVENT_NAME}: Factual vs Counterfactual ({DAMAGE_COLUMN})', 
#              fontsize=16, fontweight='bold', y=1.02)

# Save the figure
output_file_main = OUTPUT_DIR / f'damage_{EVENT_NAME.lower()}_{DAMAGE_COLUMN}_comparison.png'
plt.savefig(output_file_main, dpi=300, bbox_inches='tight')

#%%
# ===== CREATE SEPARATE BAR CHART (ONLY FOR TOTAL DAMAGE) =====
if DAMAGE_COLUMN == 'total_damage':
    print("Creating separate total damage bar chart...")
    # Calculate total damage for bar chart
    total_cf0 = cf0_damage.sum()
    total_cf8 = cf8_damage.sum()
    if EVENT_NAME == "Idai" or EVENT_NAME == "Kenneth":
        bar_label = 'Total Damage [USD 2019]'
    elif EVENT_NAME == "Freddy":
        bar_label = 'Total Damage [USD 2023]'
    else:
        bar_label = 'Total Damage [USD]'

    # Create separate bar chart figure
    fig_bar, ax_bar = plt.subplots(1, 1, figsize=(6, 6))

    # Create bar chart
    scenarios = ['Counterfactual', 'Factual']
    totals = [total_cf8, total_cf0]
    bars = ax_bar.bar(scenarios, totals, color=['lightcoral', 'lightcoral'],
                      alpha=1, width=0.4)

    # Style the bar chart
    ax_bar.set_ylabel(bar_label, fontsize=12, fontweight='bold')
    ax_bar.set_title(f'Total Damage - {EVENT_NAME}', fontsize=14, fontweight='bold')
    ax_bar.set_axisbelow(True)

    # Calculate difference for annotation
    difference = abs(total_cf0 - total_cf8)
    percentage_diff = (difference / total_cf0) * 100

    # Add annotation showing climate change attribution
    # Get bar positions
    bar_positions = [bar.get_x() + bar.get_width()/2 for bar in bars]
    bar_heights = [bar.get_height() for bar in bars]

    # Determine which bar is higher
    higher_bar_idx = 0 if totals[0] > totals[1] else 1
    lower_bar_idx = 1 - higher_bar_idx

    # Annotation parameters for vertical difference line
    line_x_position = bar_positions[1] + (bar_positions[1] - bar_positions[0]) * 0.3  # Position to the right of the rightmost bar
    line_extension = (bar_positions[1] - bar_positions[0]) * 0.04  # Small horizontal ticks

    # Draw the vertical difference line
    # Main vertical line showing the height difference
    ax_bar.plot([line_x_position, line_x_position], 
               [bar_heights[lower_bar_idx], bar_heights[higher_bar_idx]], 
               'k-', linewidth=2)

    # Add extended dashed horizontal lines from y-axis to the vertical difference line
    left_edge = ax_bar.get_xlim()[0]
    for total in totals:
        ax_bar.plot([left_edge, line_x_position], [total, total], 
                   linestyle='--', color='gray', alpha=0.7, linewidth=1)

    # Small horizontal ticks at both ends
    ax_bar.plot([line_x_position - line_extension, line_x_position + line_extension], 
               [bar_heights[lower_bar_idx], bar_heights[lower_bar_idx]], 
               'k-', linewidth=1.5)
    ax_bar.plot([line_x_position - line_extension, line_x_position + line_extension], 
               [bar_heights[higher_bar_idx], bar_heights[higher_bar_idx]], 
               'k-', linewidth=1.5)

    # Add text annotation positioned at the center between bars
    text_x = line_x_position + (bar_positions[1] - bar_positions[0]) * 0.1  # Position to the right of the line
    text_y = ((bar_heights[lower_bar_idx] + bar_heights[higher_bar_idx]) / 2) *0.97  # Center vertically on the line

    # Format the difference text based on damage type
    # Format large numbers nicely
    if difference >= 1e9:
        diff_text = f'${difference/1e9:.1f}B ({percentage_diff:.0f}%)'
    elif difference >= 1e6:
        diff_text = f'${difference/1e6:.1f}M ({percentage_diff:.0f}%)'
    elif difference >= 1e3:
        diff_text = f'${difference/1e3:.1f}K ({percentage_diff:.0f}%)'
    else:
        diff_text = f'${difference:.0f} ({percentage_diff:.0f}%)'

    # Add the annotation text
    ax_bar.text(text_x, text_y, diff_text, 
               ha='left', va='center', fontsize=10, rotation=0)
    # Use a small vertical offset for the second text line
    text_offset = (bar_heights[higher_bar_idx] - bar_heights[lower_bar_idx]) * 0.35  # Smaller vertical offset
    ax_bar.text(text_x, text_y + text_offset, 'Climate change \n attribution:', 
               ha='left', va='center', fontsize=9, style='italic', rotation=0)

    # Set y-axis to start from 0 with appropriate margin
    ax_bar.set_ylim(0, max(totals) * 1.1)

    # Remove top and right spines for cleaner look
    ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)

    plt.tight_layout()

    # Save the bar chart
    output_file_bar = OUTPUT_DIR / f'damage_{EVENT_NAME.lower()}_{DAMAGE_COLUMN}_totals_barchart.png'
    plt.savefig(output_file_bar, dpi=300, bbox_inches='tight')
else:
    print("Skipping bar chart generation for relative damage visualization...")
    output_file_bar = None

#%%
# ===== DETAILED DIFFERENCE PLOT =====
print("Creating detailed damage difference plot...")

if use_cartopy:
    fig2, ax = plt.subplots(1, 1, figsize=(12, 8), subplot_kw={'projection': crs})
else:
    fig2, ax = plt.subplots(1, 1, figsize=(12, 8))

# Plot only the difference with larger points
if use_cartopy:
    scatter_diff = ax.scatter(merged['lon'], merged['lat'], c=merged['damage_diff'], 
                             cmap=diff_cmap, norm=diff_norm, s=point_size*2, alpha=alpha, 
                             transform=crs)
    ax.add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
else:
    ax.set_facecolor('lightgray')
    scatter_diff = ax.scatter(merged['lon'], merged['lat'], c=merged['damage_diff'], 
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

print("Creating factual-only damage plot...")

if use_cartopy:
    fig_factual, ax_factual = plt.subplots(1, 1, figsize=(10, 8), subplot_kw={'projection': crs})
else:
    fig_factual, ax_factual = plt.subplots(1, 1, figsize=(10, 8))

# Plot only the factual (CF0) data
if use_cartopy:
    scatter_factual = ax_factual.scatter(merged['lon'], merged['lat'], c=merged[f'{DAMAGE_COLUMN}_cf0'], 
                                       cmap=damage_cmap, norm=damage_norm, s=point_size*1.5, alpha=alpha)
    ax_factual.add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
else:
    ax_factual.set_facecolor('lightgray')
    scatter_factual = ax_factual.scatter(merged['lon'], merged['lat'], c=merged[f'{DAMAGE_COLUMN}_cf0'], 
                                       cmap=damage_cmap, norm=damage_norm, s=point_size*1.5, alpha=alpha)

# Colorbar
cbar_factual = plt.colorbar(scatter_factual, ax=ax_factual, shrink=0.8, pad=0.1)
cbar_factual.set_label(damage_label, fontsize=12)

# Title
ax_factual.set_title(f'Economic Damage - {EVENT_NAME}: Factual Scenario\n({DAMAGE_COLUMN})', 
                    fontsize=14, fontweight='bold', pad=20)

# Save the factual-only plot
output_file_factual = OUTPUT_DIR / f'damage_{EVENT_NAME.lower()}_{DAMAGE_COLUMN}_factual_only.png'
plt.savefig(output_file_factual, dpi=200, bbox_inches='tight')

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
    print(f"Buildings with >${threshold_low:.0f} damage increase: {(damage_diff > threshold_low).sum()}")
    print(f"Buildings with >${threshold_low:.0f} damage decrease: {(damage_diff < -threshold_low).sum()}")
    print(f"Buildings with >${threshold_high:.0f} damage increase: {(damage_diff > threshold_high).sum()}")
    print(f"Buildings with >${threshold_high:.0f} damage decrease: {(damage_diff < -threshold_high).sum()}")

print(f"\nAnalysis complete for {EVENT_NAME}! Check the saved PNG files in: {OUTPUT_DIR}")
# print("Files created:")
# print(f"  - {output_file_main}")
# if DAMAGE_COLUMN == 'total_damage' and output_file_bar:
#     print(f"  - {output_file_bar}")
# print(f"  - {output_file_diff}")
print(f"\nTotal buildings analyzed: {len(merged)}")
if DAMAGE_COLUMN == 'relative_damage':
    print(f"Buildings with damage in CF0: {(merged[f'{DAMAGE_COLUMN}_cf0'] > 0).sum()}")
    print(f"Buildings with damage in CF-8: {(merged[f'{DAMAGE_COLUMN}_cf8'] > 0).sum()}")
else:
    print(f"Buildings with damage in CF0: {(merged[f'{DAMAGE_COLUMN}_cf0'] > 0).sum()}")
    print(f"Buildings with damage in CF-8: {(merged[f'{DAMAGE_COLUMN}_cf8'] > 0).sum()}")
# %%

# Save statistics to a text file
with open(f"summary_damage_({DAMAGE_COLUMN})_{EVENT_NAME}.txt", "w") as f:
    if DAMAGE_COLUMN == 'relative_damage':
        print(f"Damage Statistics for {EVENT_NAME} ({DAMAGE_COLUMN}):", file=f)
        print(f"CF0 max damage: {cf0_damage.max():.1f}%", file=f)
        print(f"CF0 mean damage: {cf0_damage.mean():.1f}%", file=f)
        print(f"CF0 total damage: {cf0_damage.sum():.0f}%", file=f)
        print(f"CF0 buildings with >0% damage: {(cf0_damage > 0).sum()}", file=f)
        print(f"CF-8 max damage: {cf8_damage.max():.1f}%", file=f)
        print(f"CF-8 mean damage: {cf8_damage.mean():.1f}%", file=f)
        print(f"CF-8 total damage: {cf8_damage.sum():.0f}%", file=f)
        print(f"CF-8 buildings with >0% damage: {(cf8_damage > 0).sum()}", file=f)
        print(f"Difference in total damage: {cf0_damage.sum() - cf8_damage.sum():.0f}%", file=f)
        print(f"Percentage change: {((cf0_damage.sum() - cf8_damage.sum()) / cf8_damage.sum() * 100):.1f}%", file=f)
        print(f"Max increase (CF0 vs CF-8): {damage_diff.max():.1f}%", file=f)
        print(f"Max decrease (CF0 vs CF-8): {damage_diff.min():.1f}%", file=f)
        print(f"Mean difference: {damage_diff.mean():.1f}%", file=f)

        print(f"Buildings with >10% damage increase: {(damage_diff > 10).sum()}", file=f)
        print(f"Buildings with >10% damage decrease: {(damage_diff < -10).sum()}", file=f)
        print(f"Buildings with >25% damage increase: {(damage_diff > 25).sum()}", file=f)
        print(f"Buildings with >25% damage decrease: {(damage_diff < -25).sum()}", file=f)
        print(f"Buildings with >5% damage increase: {(damage_diff > 5).sum()}", file=f)
        print(f"Buildings with >5% damage decrease: {(damage_diff < -5).sum()}", file=f)

    else:
        print(f"Damage Statistics for {EVENT_NAME} ({DAMAGE_COLUMN}):", file=f)
        print(f"CF0 max damage: ${cf0_damage.max():.0f}", file=f)
        print(f"CF0 mean damage: ${cf0_damage.mean():.0f}", file=f)
        print(f"CF0 total damage: ${cf0_damage.sum():.0f}", file=f)
        print(f"CF-8 max damage: ${cf8_damage.max():.0f}", file=f)
        print(f"CF-8 mean damage: ${cf8_damage.mean():.0f}", file=f)
        print(f"CF-8 total damage: ${cf8_damage.sum():.0f}", file=f)
        print(f"Difference in total damage: ${cf0_damage.sum() - cf8_damage.sum():.0f}", file=f)
        print(f"Percentage change: {((cf0_damage.sum() - cf8_damage.sum()) / cf8_damage.sum() * 100):.1f}%", file=f)
        print(f"Max increase (CF0 vs CF-8): ${damage_diff.max():.0f}", file=f)
        print(f"Max decrease (CF0 vs CF-8): ${damage_diff.min():.0f}", file=f)
        print(f"Mean difference: ${damage_diff.mean():.0f}", file=f)

        threshold_low = diff_max * 0.1
        threshold_high = diff_max * 0.25
        print(f"Buildings with >${threshold_low:.0f} damage increase: {(damage_diff > threshold_low).sum()}", file=f)
        print(f"Buildings with >${threshold_low:.0f} damage decrease: {(damage_diff < -threshold_low).sum()}", file=f)
        print(f"Buildings with >${threshold_high:.0f} damage increase: {(damage_diff > threshold_high).sum()}", file=f)
        print(f"Buildings with >${threshold_high:.0f} damage decrease: {(damage_diff < -threshold_high).sum()}", file=f)


    print(f"\nTotal buildings analyzed: {len(merged)}", file=f)

    if DAMAGE_COLUMN == 'relative_damage':
        print(f"Damage Statistics for {EVENT_NAME} ({DAMAGE_COLUMN}):", file=f)
        print(f"Buildings with damage in CF0: {(merged[f'{DAMAGE_COLUMN}_cf0'] > 0).sum()}", file=f)
        print(f"Buildings with damage in CF-8: {(merged[f'{DAMAGE_COLUMN}_cf8'] > 0).sum()}", file=f)
    else:
        print(f"Damage Statistics for {EVENT_NAME} ({DAMAGE_COLUMN}):", file=f)
        print(f"Buildings with damage in CF0: {(merged[f'{DAMAGE_COLUMN}_cf0'] > 0).sum()}", file=f)
        print(f"Buildings with damage in CF-8: {(merged[f'{DAMAGE_COLUMN}_cf8'] > 0).sum()}", file=f)

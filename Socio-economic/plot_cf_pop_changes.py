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
from pathlib import Path


# ===== CONFIGURATION =====
# Set your event name here
EVENT_NAME = "Idai"  # Change this to: "Kenneth", "Freddy", etc.

# Choose population column to analyze
POPULATION_COLUMN = "population"  # Main population column to analyze

# Inundation depth threshold (in meters)
INUNDATION_THRESHOLD = 0  # Only consider areas with >0.05m flooding

# Base paths - update these as needed
prefix = "p:/" if platform.system() == "Windows" else "/p/"

BASE_RUN_PATH = Path("C:/Code/Paper_1/Data_submission") # For Idai
OUTPUT_DIR = Path("p:/11210471-001-compass/04_Results/socioeconomic/pop")

# ===== DYNAMIC FILE PATHS =====
file_cf0 = BASE_RUN_PATH / "fiat" / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" / "output" / "spatial_with_pop_and_flood_worldpop.fgb"
file_cf8 = BASE_RUN_PATH / "fiat" / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" / "output" / "spatial_with_pop_and_flood_hist_exposure.fgb"

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
    print(f"Warning: {POPULATION_COLUMN} not found in CF0 data. Available columns: {list(gdf_cf0.columns)}")
    # Try alternative column names
    if 'pop' in gdf_cf0.columns:
        POPULATION_COLUMN = 'pop'
    elif 'total_population' in gdf_cf0.columns:
        POPULATION_COLUMN = 'total_population'
    else:
        raise ValueError("No suitable population column found!")

print(f"Final population column used: {POPULATION_COLUMN}")

# Check if inundation depth column exists
if 'inun_depth' not in gdf_cf0.columns:
    print("Warning: 'inun_depth' column not found. Available columns with 'depth' or 'inun':")
    depth_cols = [col for col in gdf_cf0.columns if 'depth' in col.lower() or 'inun' in col.lower()]
    print(depth_cols)
    if depth_cols:
        inun_col = depth_cols[0]
        print(f"Using column: {inun_col}")
    else:
        raise ValueError("No suitable inundation depth column found!")
else:
    inun_col = 'inun_depth'

#%%
# ===== FILTER BY INUNDATION DEPTH =====
print(f"Filtering data by inundation depth > {INUNDATION_THRESHOLD}m...")
print(f"CF0 records before filtering: {len(gdf_cf0)}")
print(f"CF-8 records before filtering: {len(gdf_cf8)}")

# Apply inundation depth filter
gdf_cf0_filtered = gdf_cf0[gdf_cf0[inun_col] > INUNDATION_THRESHOLD].copy()
gdf_cf8_filtered = gdf_cf8[gdf_cf8[inun_col] > INUNDATION_THRESHOLD].copy()

print(f"CF0 records after filtering: {len(gdf_cf0_filtered)}")
print(f"CF-8 records after filtering: {len(gdf_cf8_filtered)}")

# ===== EXTRACT COORDINATES =====
print("Extracting coordinates from geometry centroids...")
# Extract x, y coordinates from geometry centroids
gdf_cf0_filtered['centroid'] = gdf_cf0_filtered.geometry.centroid
gdf_cf0_filtered['x'] = gdf_cf0_filtered['centroid'].x
gdf_cf0_filtered['y'] = gdf_cf0_filtered['centroid'].y

gdf_cf8_filtered['centroid'] = gdf_cf8_filtered.geometry.centroid
gdf_cf8_filtered['x'] = gdf_cf8_filtered['centroid'].x
gdf_cf8_filtered['y'] = gdf_cf8_filtered['centroid'].y

print(f"Sample coordinates CF0: x={gdf_cf0_filtered['x'].iloc[0]:.2f}, y={gdf_cf0_filtered['y'].iloc[0]:.2f}")
print(f"Sample coordinates CF-8: x={gdf_cf8_filtered['x'].iloc[0]:.2f}, y={gdf_cf8_filtered['y'].iloc[0]:.2f}")

# ===== MERGE DATA FOR DIFFERENCE CALCULATION =====
print("Merging data for difference calculation...")
# Merge on object_id to compare same locations
# Use coordinates from CF0 data, but handle cases where locations might only exist in one scenario
merged = gdf_cf0_filtered[['object_id', 'x', 'y', POPULATION_COLUMN, inun_col]].merge(
    gdf_cf8_filtered[['object_id', 'x', 'y', POPULATION_COLUMN, inun_col]], 
    on='object_id', 
    how='outer', 
    suffixes=('_cf0', '_cf8')
)

# Handle coordinates: use CF0 coordinates where available, otherwise CF8
merged['x'] = merged['x_cf0'].fillna(merged['x_cf8'])
merged['y'] = merged['y_cf0'].fillna(merged['y_cf8'])

# Handle inundation depth: use CF0 values where available, otherwise CF8
merged['inun_depth'] = merged[f'{inun_col}_cf0'].fillna(merged[f'{inun_col}_cf8'])

# Drop the temporary coordinate columns
merged = merged.drop([col for col in merged.columns if col.endswith('_cf0') or col.endswith('_cf8')], axis=1)

# Re-add the population columns with proper names
merged[f'{POPULATION_COLUMN}_cf0'] = gdf_cf0_filtered.set_index('object_id')[POPULATION_COLUMN].reindex(merged['object_id']).fillna(0).values
merged[f'{POPULATION_COLUMN}_cf8'] = gdf_cf8_filtered.set_index('object_id')[POPULATION_COLUMN].reindex(merged['object_id']).fillna(0).values

# Calculate difference in exposed population (CF0 - CF8)
merged['pop_exposure_diff'] = merged[f'{POPULATION_COLUMN}_cf0'] - merged[f'{POPULATION_COLUMN}_cf8']

# Remove locations with no coordinates
merged = merged.dropna(subset=['x', 'y'])

print(f"Merged data shape: {merged.shape}")
print(f"Sample merged data:")
print(merged[['object_id', 'x', 'y', f'{POPULATION_COLUMN}_cf0', f'{POPULATION_COLUMN}_cf8', 'pop_exposure_diff']].head())

# ===== STATISTICS =====
print(f"\nPopulation Exposure Statistics for {EVENT_NAME}:")
cf0_pop = merged[f'{POPULATION_COLUMN}_cf0']
cf8_pop = merged[f'{POPULATION_COLUMN}_cf8']
pop_diff = merged['pop_exposure_diff']

print(f"CF0 max exposed population: {cf0_pop.max():.0f} people")
print(f"CF0 mean exposed population: {cf0_pop.mean():.1f} people")
print(f"CF0 total exposed population: {cf0_pop.sum():.0f} people")
print(f"CF0 locations with >0 exposed population: {(cf0_pop > 0).sum()}")
print(f"CF-8 max exposed population: {cf8_pop.max():.0f} people")
print(f"CF-8 mean exposed population: {cf8_pop.mean():.1f} people")
print(f"CF-8 total exposed population: {cf8_pop.sum():.0f} people")
print(f"CF-8 locations with >0 exposed population: {(cf8_pop > 0).sum()}")
print(f"Difference in total exposed population: {cf0_pop.sum() - cf8_pop.sum():.0f} people")
if cf8_pop.sum() > 0:
    print(f"Percentage change: {((cf0_pop.sum() - cf8_pop.sum()) / cf8_pop.sum() * 100):.1f}%")
else:
    print("Percentage change: Cannot calculate (CF-8 has no exposed population)")
print(f"Max increase in exposed population (CF0 vs CF-8): {pop_diff.max():.0f} people")
print(f"Max decrease in exposed population (CF0 vs CF-8): {pop_diff.min():.0f} people")
print(f"Mean difference in exposed population: {pop_diff.mean():.1f} people")

# ===== DETERMINE COORDINATE SYSTEM =====
try:
    # Get approximate center coordinates
    center_x = merged['x'].mean()
    center_y = merged['y'].mean()
    
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
        print(f"Detected projected coordinates - UTM zone {utm_zone}, Southern: {southern}")
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

# Create projection for cartopy with better fallback options
use_cartopy = True
crs = None

try:
    if utm_zone is not None:
        # Try UTM projection first
        crs = ccrs.UTM(zone=utm_zone, southern_hemisphere=southern)
        print(f"Using UTM projection: zone {utm_zone}, southern={southern}")
    else:
        # Try PlateCarree for lat/lon data
        crs = ccrs.PlateCarree()
        print("Using PlateCarree projection for lat/lon data")
except Exception as e:
    print(f"UTM/PlateCarree projection failed: {e}")
    try:
        # Fallback to PlateCarree
        crs = ccrs.PlateCarree()
        print("Using PlateCarree projection as fallback")
    except Exception as e2:
        print(f"All cartopy projections failed: {e2}")
        use_cartopy = False

if not use_cartopy:
    print("Cartopy disabled - using matplotlib only")

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
pop_cmap = plt.get_cmap('Blues')  # Blue colormap for population
pop_label = 'Exposed Population [people]'
vmin, vmax = 0, max_pop

# Difference colormap with discrete intervals
diff_max = max(abs(pop_diff.quantile(0.05)), abs(pop_diff.quantile(0.95)))
diff_min = min(abs(pop_diff.quantile(0.05)), abs(pop_diff.quantile(0.95)))
diff_boundaries = np.linspace(diff_min, diff_max, 11)
diff_norm = BoundaryNorm(diff_boundaries, ncolors=256, clip=True)
diff_cmap = plt.get_cmap('Reds')  # Red for increases, Blue for decreases

print(f"Population boundaries: {boundaries}")
print(f"Difference boundaries: {diff_boundaries}")

# ===== CREATE MAIN COMPARISON PLOT (3-PANEL) =====
print("Creating population exposure comparison plots...")
print(f"Using cartopy: {use_cartopy}")

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
    scatter1 = axes[0].scatter(merged['x'], merged['y'], c=merged[f'{POPULATION_COLUMN}_cf0'], 
                          cmap=pop_cmap, norm=pop_norm, s=point_size, alpha=alpha, 
                          transform=crs)
    # Add multiple cartopy features for better visualization
    try:
        axes[0].add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
        axes[0].add_feature(cfeature.OCEAN, color='white', alpha=0.5)
        axes[0].add_feature(cfeature.COASTLINE, color='black', alpha=0.3, linewidth=0.5)
    except:
        print("Some cartopy features failed to load")
else:
    scatter1 = axes[0].scatter(merged['x'], merged['y'], c=merged[f'{POPULATION_COLUMN}_cf0'], 
                          cmap=pop_cmap, norm=pop_norm, s=point_size, alpha=alpha)
axes[0].set_title('Factual', fontsize=12, fontweight='bold')

# Plot CF-8 (Counterfactual)
if use_cartopy:
    scatter2 = axes[1].scatter(merged['x'], merged['y'], c=merged[f'{POPULATION_COLUMN}_cf8'], 
                          cmap=pop_cmap, norm=pop_norm, s=point_size, alpha=alpha, 
                          transform=crs)
    try:
        axes[1].add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
        axes[1].add_feature(cfeature.OCEAN, color='white', alpha=0.5)
        axes[1].add_feature(cfeature.COASTLINE, color='black', alpha=0.3, linewidth=0.5)
    except:
        print("Some cartopy features failed to load")
else:
    scatter2 = axes[1].scatter(merged['x'], merged['y'], c=merged[f'{POPULATION_COLUMN}_cf8'], 
                          cmap=pop_cmap, norm=pop_norm, s=point_size, alpha=alpha)
axes[1].set_title('Counterfactual', fontsize=12, fontweight='bold')

# Plot Difference
if use_cartopy:
    scatter3 = axes[2].scatter(merged['x'], merged['y'], c=merged['pop_exposure_diff'], 
                          cmap=diff_cmap, norm=diff_norm, s=point_size, alpha=alpha, 
                          transform=crs)
    try:
        axes[2].add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
        axes[2].add_feature(cfeature.OCEAN, color='white', alpha=0.5)
        axes[2].add_feature(cfeature.COASTLINE, color='black', alpha=0.3, linewidth=0.5)
    except:
        print("Some cartopy features failed to load")
else:
    scatter3 = axes[2].scatter(merged['x'], merged['y'], c=merged['pop_exposure_diff'], 
                          cmap=diff_cmap, norm=diff_norm, s=point_size, alpha=alpha)
axes[2].set_title('Population Exposure Changes: Factual vs Counterfactual', fontsize=12, fontweight='bold')

# Add colorbars
from mpl_toolkits.axes_grid1 import make_axes_locatable

# === Colorbar for scatter1
cbar1 = fig.colorbar(scatter1, ax=axes[0], orientation='vertical', shrink=0.7, pad=0.02)
cbar1.set_label(pop_label, fontsize=10)

# === Colorbar for scatter2
cbar2 = fig.colorbar(scatter2, ax=axes[1], orientation='vertical', shrink=0.7, pad=0.02)
cbar2.set_label(pop_label, fontsize=10)

# === Colorbar for scatter3 (difference)
cbar3 = fig.colorbar(scatter3, ax=axes[2], orientation='vertical', shrink=0.7, pad=0.02)
cbar3.set_label(f'{pop_label} Difference', fontsize=10)

# Adjust layout and add main title
plt.tight_layout()
plt.suptitle(f'Population Exposure Analysis - {EVENT_NAME}: Factual vs Counterfactual (depth > {INUNDATION_THRESHOLD}m)', 
             fontsize=16, fontweight='bold', y=1.02)

# Save the figure
output_file_main = OUTPUT_DIR / f'population_exposure_{EVENT_NAME.lower()}_comparison.png'
plt.savefig(output_file_main, dpi=300, bbox_inches='tight')

#%%
# ===== CREATE SEPARATE BAR CHART FOR TOTAL EXPOSED POPULATION =====
print("Creating separate total exposed population bar chart...")
# Calculate total exposed population for bar chart
total_cf0 = cf0_pop.sum()
total_cf8 = cf8_pop.sum()
bar_label = 'Total Exposed Population [people]'

# Create separate bar chart figure
fig_bar, ax_bar = plt.subplots(1, 1, figsize=(6, 6))

# Create bar chart
scenarios = ['Counterfactual', 'Factual']
totals = [total_cf8, total_cf0]
bars = ax_bar.bar(scenarios, totals, color=['lightblue', 'lightblue'],
                  alpha=1, width=0.4)

# Style the bar chart
ax_bar.set_ylabel(bar_label, fontsize=12, fontweight='bold')
ax_bar.set_title(f'Total Exposed Population - {EVENT_NAME}', fontsize=14, fontweight='bold')
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
line_x_position = bar_positions[1] + (bar_positions[1] - bar_positions[0]) * 0.3
line_extension = (bar_positions[1] - bar_positions[0]) * 0.04

# Draw the vertical difference line
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
text_x = line_x_position + (bar_positions[1] - bar_positions[0]) * 0.1
text_y = ((bar_heights[lower_bar_idx] + bar_heights[higher_bar_idx]) / 2) * 0.97

# Format the difference text for population
if difference >= 1e6:
    diff_text = f'{difference/1e6:.1f}M people ({percentage_diff:.0f}%)'
elif difference >= 1e3:
    diff_text = f'{difference/1e3:.1f}K people ({percentage_diff:.0f}%)'
else:
    diff_text = f'{difference:.0f} people ({percentage_diff:.0f}%)'

# Add the annotation text
ax_bar.text(text_x, text_y, diff_text, 
           ha='left', va='center', fontsize=10, rotation=0)
text_offset = (bar_heights[higher_bar_idx] - bar_heights[lower_bar_idx]) * 0.8
ax_bar.text(text_x, text_y + text_offset, 'Climate change \n attribution:', 
           ha='left', va='center', fontsize=9, style='italic', rotation=0)

# Set y-axis to start from 0 with appropriate margin
ax_bar.set_ylim(0, max(totals) * 1.1)

# Remove top and right spines for cleaner look
ax_bar.spines['top'].set_visible(False)
ax_bar.spines['right'].set_visible(False)

plt.tight_layout()

# Save the bar chart
output_file_bar = OUTPUT_DIR / f'population_exposure_{EVENT_NAME.lower()}_totals_barchart.png'
plt.savefig(output_file_bar, dpi=300, bbox_inches='tight')

#%%
# ===== DETAILED DIFFERENCE PLOT =====
print("Creating detailed population exposure difference plot...")

if use_cartopy:
    fig2, ax = plt.subplots(1, 1, figsize=(12, 8), subplot_kw={'projection': crs})
else:
    fig2, ax = plt.subplots(1, 1, figsize=(12, 8))

# Plot only the difference with larger points
if use_cartopy:
    scatter_diff = ax.scatter(merged['x'], merged['y'], c=merged['pop_exposure_diff'], 
                             cmap=diff_cmap, norm=diff_norm, s=point_size*2, alpha=alpha, 
                             transform=crs)
    ax.add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
else:
    scatter_diff = ax.scatter(merged['x'], merged['y'], c=merged['pop_exposure_diff'], 
                             cmap=diff_cmap, norm=diff_norm, s=point_size*2, alpha=alpha)

# Colorbar
cbar = plt.colorbar(scatter_diff, ax=ax, shrink=0.8, pad=0.1)
cbar.set_label(f'{pop_label} Difference', fontsize=12)

# Title
ax.set_title(f'Population Exposure Changes - {EVENT_NAME}: Factual vs Counterfactual\n(Red = More exposed people, Blue = Fewer exposed people)', 
             fontsize=14, fontweight='bold', pad=20)

# Save
output_file_diff = OUTPUT_DIR / f'population_exposure_{EVENT_NAME.lower()}_difference_only.png'
plt.savefig(output_file_diff, dpi=300, bbox_inches='tight')

# ===== CREATE FACTUAL-ONLY PLOT =====
print("Creating factual-only population exposure plot...")
print(merged[['x', 'y']].head())
use_cartopy = True

# print(f"merged.crs: {merged.crs}")

gdf = gpd.GeoDataFrame(
    merged,
    geometry=gpd.points_from_xy(merged['x'], merged['y']),
    crs=f"EPSG:4326"  # UTM Zone 36S (adjust if needed)
)

# Convert to lat/lon
gdf = gdf.to_crs("EPSG:4326")

print(f"gdf.crs: {gdf.crs}")

# Use these for plotting
merged['lon'] = gdf.geometry.x
merged['lat'] = gdf.geometry.y

crs = ccrs.PlateCarree()

if use_cartopy:
    fig_factual, ax_factual = plt.subplots(1, 1, figsize=(10, 8), subplot_kw={'projection': crs})
else:
    fig_factual, ax_factual = plt.subplots(1, 1, figsize=(10, 8))

# Plot only the factual (CF0) data
if use_cartopy:
    scatter_factual = ax_factual.scatter(merged['lon'], merged['lat'], c=merged[f'{POPULATION_COLUMN}_cf0'], 
                                        cmap=pop_cmap, norm=pop_norm, s=point_size*1.5, alpha=alpha, 
                                        transform=ccrs.PlateCarree())
    try:
        ax_factual.add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
        ax_factual.add_feature(cfeature.OCEAN, color='white', alpha=0.5)
        ax_factual.add_feature(cfeature.COASTLINE, color='black', alpha=0.3, linewidth=0.5)
    except:
        print("Some cartopy features failed to load")
else:
    scatter_factual = ax_factual.scatter(merged['lon'], merged['lat'], c=merged[f'{POPULATION_COLUMN}_cf0'], 
                                       cmap=pop_cmap, norm=pop_norm, s=point_size*1.5, alpha=alpha)

# Colorbar
cbar_factual = plt.colorbar(scatter_factual, ax=ax_factual, shrink=0.8, pad=0.1)
cbar_factual.set_label(pop_label, fontsize=12)

# Title
ax_factual.set_title(f'Population Exposure - {EVENT_NAME}: Factual Scenario\n(depth > {INUNDATION_THRESHOLD}m)', 
                    fontsize=14, fontweight='bold', pad=20)

# Save the factual-only plot
output_file_factual = OUTPUT_DIR / f'population_exposure_{EVENT_NAME.lower()}_factual_only.png'
plt.savefig(output_file_factual, dpi=300, bbox_inches='tight')

# ===== SUMMARY STATISTICS =====
print(f"\nLocations with significant population exposure changes:")
# Define thresholds based on population scale
low_threshold = max(10, diff_max * 0.1)
high_threshold = max(50, diff_max * 0.25)

print(f"Locations with >{low_threshold:.0f} people more exposed: {(pop_diff > low_threshold).sum()}")
print(f"Locations with >{low_threshold:.0f} people less exposed: {(pop_diff < -low_threshold).sum()}")
print(f"Locations with >{high_threshold:.0f} people more exposed: {(pop_diff > high_threshold).sum()}")
print(f"Locations with >{high_threshold:.0f} people less exposed: {(pop_diff < -high_threshold).sum()}")

# Inundation depth statistics
print(f"\nInundation depth statistics (for flooded areas > {INUNDATION_THRESHOLD}m):")
print(f"Mean inundation depth: {merged['inun_depth'].mean():.2f}m")
print(f"Max inundation depth: {merged['inun_depth'].max():.2f}m")
print(f"Median inundation depth: {merged['inun_depth'].median():.2f}m")

print(f"\nAnalysis complete for {EVENT_NAME}! Check the saved PNG files in: {OUTPUT_DIR}")
print("Files created:")
# print(f"  - {output_file_main}")
# print(f"  - {output_file_bar}")
# print(f"  - {output_file_diff}")
print(f"\nTotal locations analyzed: {len(merged)}")
print(f"Locations with exposed population in CF0: {(merged[f'{POPULATION_COLUMN}_cf0'] > 0).sum()}")
print(f"Locations with exposed population in CF-8: {(merged[f'{POPULATION_COLUMN}_cf8'] > 0).sum()}")
# %%

with open(f"{OUTPUT_DIR}/summary_population_{EVENT_NAME}.txt", "w") as f:
    print(f"\nPopulation Exposure Statistics for {EVENT_NAME}:", file=f)
    print(f"\nCF0 is Factual, CF8 is Counterfactual", file=f)

    print(f"CF0 max exposed population: {cf0_pop.max():.0f} people", file=f)
    print(f"CF0 mean exposed population: {cf0_pop.mean():.1f} people", file=f)
    print(f"CF0 total exposed population: {cf0_pop.sum():.0f} people", file=f)
    print(f"CF0 locations with >0 exposed population: {(cf0_pop > 0).sum()}", file=f)
    print(f"CF-8 max exposed population: {cf8_pop.max():.0f} people", file=f)
    print(f"CF-8 mean exposed population: {cf8_pop.mean():.1f} people", file=f)
    print(f"CF-8 total exposed population: {cf8_pop.sum():.0f} people", file=f)
    print(f"CF-8 locations with >0 exposed population: {(cf8_pop > 0).sum()}", file=f)
    print(f"Difference in total exposed population: {cf0_pop.sum() - cf8_pop.sum():.0f} people", file=f)
    if cf8_pop.sum() > 0:
        print(f"Percentage change: {((cf0_pop.sum() - cf8_pop.sum()) / cf8_pop.sum() * 100):.1f}%", file=f)
    else:
        print("Percentage change: Cannot calculate (CF-8 has no exposed population)", file=f)

    print(f"Max increase in exposed population (CF0 vs CF-8): {pop_diff.max():.0f} people", file=f)
    print(f"Max decrease in exposed population (CF0 vs CF-8): {pop_diff.min():.0f} people", file=f)
    print(f"Mean difference in exposed population: {pop_diff.mean():.1f} people", file=f)

    print(f"Locations with >{low_threshold:.0f} people more exposed: {(pop_diff > low_threshold).sum()}", file=f)
    print(f"Locations with >{low_threshold:.0f} people less exposed: {(pop_diff < -low_threshold).sum()}", file=f)
    print(f"Locations with >{high_threshold:.0f} people more exposed: {(pop_diff > high_threshold).sum()}", file=f)
    print(f"Locations with >{high_threshold:.0f} people less exposed: {(pop_diff < -high_threshold).sum()}", file=f)

    # Inundation depth statistics
    print(f"\nInundation depth statistics (for flooded areas > {INUNDATION_THRESHOLD}m):", file=f)
    print(f"Mean inundation depth: {merged['inun_depth'].mean():.2f}m", file=f)
    print(f"Max inundation depth: {merged['inun_depth'].max():.2f}m", file=f)
    print(f"Median inundation depth: {merged['inun_depth'].median():.2f}m", file=f)

    print(f"\nAnalysis complete for {EVENT_NAME}! Check the saved PNG files in: {OUTPUT_DIR}", file=f)
    print("Files created:", file=f)
    print(f"\nTotal locations analyzed: {len(merged)}", file=f)
    print(f"Locations with exposed population in CF0: {(merged[f'{POPULATION_COLUMN}_cf0'] > 0).sum()}", file=f)
    print(f"Locations with exposed population in CF-8: {(merged[f'{POPULATION_COLUMN}_cf8'] > 0).sum()}", file=f)


# %%

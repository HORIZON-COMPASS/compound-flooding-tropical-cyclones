import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import BoundaryNorm, ListedColormap
import warnings
warnings.filterwarnings('ignore')

# File paths
from pathlib import Path

file_cf0 = "/p/11210471-001-compass/03_Runs/test/Kenneth/sfincs/event_tp_era5_hourly_zarr_CF0_GTSMv41opendap_CF0_no_wind_CF0/sfincs_map.nc"
file_cf8 = "/p/11210471-001-compass/03_Runs/test/Kenneth/sfincs/event_tp_era5_hourly_zarr_CF-8_GTSMv41opendap_CF0_no_wind_CF0/sfincs_map.nc"
output_dir = Path("/p/11210471-001-compass/04_Results/CF_figs")

# Create output directory if it doesn't exist
output_dir.mkdir(parents=True, exist_ok=True)

# Create output directory if it doesn't exist
output_dir.mkdir(parents=True, exist_ok=True)

# Load datasets (with decode_times=False to avoid time parsing issues)
print("Loading CF0 dataset...")
ds_cf0 = xr.open_dataset(file_cf0, decode_times=False)

print("Loading CF-8 dataset...")
ds_cf8 = xr.open_dataset(file_cf8, decode_times=False)

print("Dataset variables available:")
print(f"CF0: {list(ds_cf0.data_vars.keys())}")
print(f"CF-8: {list(ds_cf8.data_vars.keys())}")

# Extract water depth (zs - zb) and find maximum over time
print("Extracting water depths...")
print(f"CF0 zs shape: {ds_cf0.zs.shape}")
print(f"CF0 zb shape: {ds_cf0.zb.shape}")
print(f"CF-8 zs shape: {ds_cf8.zs.shape}")
print(f"CF-8 zb shape: {ds_cf8.zb.shape}")

# Calculate water depth (water level minus bed level) for each time step
print("Calculating water depths (zs - zb)...")
water_depth_cf0 = ds_cf0.zs - ds_cf0.zb  # Water depth at each time step for CF0
water_depth_cf8 = ds_cf8.zs - ds_cf8.zb  # Water depth at each time step for CF-8

# Take maximum across time dimension to get overall maximum flood depth
print("Finding maximum water depths over time...")
zsmax_cf0 = water_depth_cf0.max(dim='time')  # Maximum water depth for CF0
zsmax_cf8 = water_depth_cf8.max(dim='time')  # Maximum water depth for CF-8

print(f"After time max - CF0 shape: {zsmax_cf0.shape}")
print(f"After time max - CF-8 shape: {zsmax_cf8.shape}")

# Calculate difference (CF-8 minus CF0)
print("Calculating differences...")
diff = zsmax_cf8 - zsmax_cf0

# Mask areas with no flooding in either scenario and apply _FillValue
print("Applying masks and cleaning data...")
# Replace fill values with NaN
zsmax_cf0 = zsmax_cf0.where(zsmax_cf0 != -99999.0)
zsmax_cf8 = zsmax_cf8.where(zsmax_cf8 != -99999.0)

# Only show areas with positive water depth (above ground)
zsmax_cf0 = zsmax_cf0.where(zsmax_cf0 > 0)
zsmax_cf8 = zsmax_cf8.where(zsmax_cf8 > 0)

# Optional: only show differences where there's significant flooding
# diff = diff.where((zsmax_cf0 > 0.01) | (zsmax_cf8 > 0.01))

# Print some statistics
print(f"\nStatistics:")
print(f"CF0 max flood depth: {zsmax_cf0.max().values:.3f} m")
print(f"CF-8 max flood depth: {zsmax_cf8.max().values:.3f} m")
print(f"Maximum increase (CF-8 vs CF0): {diff.max().values:.3f} m")
print(f"Maximum decrease (CF-8 vs CF0): {diff.min().values:.3f} m")
print(f"Mean difference: {diff.mean().values:.3f} m")

# Create the plot
print("Creating plots...")
try:
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), 
                            subplot_kw={'projection': ccrs.UTM(zone=37, southern_hemisphere=True)})
    use_cartopy = True
except Exception as e:
    print(f"Cartopy projection failed: {e}")
    print("Falling back to simple plots without projection...")
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    use_cartopy = False

# Define colormaps and levels
flood_levels = np.arange(0, 3.1, 0.1)
flood_cmap = plt.cm.Blues

# Difference levels (symmetric around 0)
dif_min = -0.5
dif_max = 0.5
diff_levels = np.linspace(dif_min, dif_max, 21)
diff_cmap = plt.cm.RdBu  # Red for increases, Blue for decreases

# Plot CF0
ax1 = axes[0]
if use_cartopy:
    im1 = zsmax_cf0.plot(ax=ax1, levels=flood_levels, cmap=flood_cmap, 
                         add_colorbar=False, transform=ccrs.UTM(zone=37, southern_hemisphere=True),
                         x='x', y='y')
    ax1.coastlines(resolution='10m')
else:
    im1 = zsmax_cf0.plot(ax=ax1, levels=flood_levels, cmap=flood_cmap, 
                         add_colorbar=False, x='x', y='y')
ax1.set_title('CF0: Maximum Flood Depth Above Ground [m]', fontsize=12, fontweight='bold')

# Plot CF-8
ax2 = axes[1]
if use_cartopy:
    im2 = zsmax_cf8.plot(ax=ax2, levels=flood_levels, cmap=flood_cmap, 
                         add_colorbar=False, transform=ccrs.UTM(zone=37, southern_hemisphere=True),
                         x='x', y='y')
    ax2.coastlines(resolution='10m')
else:
    im2 = zsmax_cf8.plot(ax=ax2, levels=flood_levels, cmap=flood_cmap, 
                         add_colorbar=False, x='x', y='y')
ax2.set_title('CF-8: Maximum Flood Depth Above Ground [m]', fontsize=12, fontweight='bold')

# Plot difference
ax3 = axes[2]
if use_cartopy:
    im3 = diff.plot(ax=ax3, levels=diff_levels, cmap=diff_cmap, 
                    add_colorbar=False, transform=ccrs.UTM(zone=37, southern_hemisphere=True),
                    x='x', y='y')
    ax3.coastlines(resolution='10m')
else:
    im3 = diff.plot(ax=ax3, levels=diff_levels, cmap=diff_cmap, 
                    add_colorbar=False, x='x', y='y')
ax3.set_title('Difference (CF-8 minus CF0) [m]', fontsize=12, fontweight='bold')

# Add colorbars
cbar1 = plt.colorbar(im1, ax=ax1, shrink=0.8, pad=0.1)
cbar1.set_label('Max Flood Depth Above Ground [m]', fontsize=10)

cbar2 = plt.colorbar(im2, ax=ax2, shrink=0.8, pad=0.1)
cbar2.set_label('Max Flood Depth Above Ground [m]', fontsize=10)

cbar3 = plt.colorbar(im3, ax=ax3, shrink=0.8, pad=0.1)
cbar3.set_label('Depth Difference [m]', fontsize=10)

# Adjust layout and add main title
plt.tight_layout()
plt.suptitle('SFINCS Flood Depth Above Ground Comparison: CF0 vs CF-8', 
             fontsize=16, fontweight='bold', y=1.02)

# Save the figure
plt.savefig(output_dir / 'sfincs_flood_depth_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

# Additional analysis: create a simple difference-only plot
print("Creating detailed difference plot...")
try:
    fig2, ax = plt.subplots(1, 1, figsize=(12, 8),
                           subplot_kw={'projection': ccrs.UTM(zone=37, southern_hemisphere=True)})
    use_cartopy2 = True
except:
    fig2, ax = plt.subplots(1, 1, figsize=(12, 8))
    use_cartopy2 = False

# Plot only the difference with better resolution
if use_cartopy2:
    im_diff = diff.plot(ax=ax, levels=diff_levels, cmap='RdBu', 
                       add_colorbar=False, transform=ccrs.UTM(zone=37, southern_hemisphere=True),
                       x='x', y='y')
else:
    im_diff = diff.plot(ax=ax, levels=diff_levels, cmap='RdBu', 
                       add_colorbar=False, x='x', y='y')

# Colorbar
cbar = plt.colorbar(im_diff, ax=ax, shrink=0.8, pad=0.1)
cbar.set_label('Flood Depth Difference (CF-8 minus CF0) [m]', fontsize=12)

# Title
ax.set_title('Flood Depth Changes: CF-8 vs CF0\n(Red = Increased flooding, Blue = Decreased flooding)', 
             fontsize=14, fontweight='bold', pad=20)

# Save
plt.savefig(output_dir / 'sfincs_difference_only.png', dpi=300, bbox_inches='tight')
plt.show()

# Print summary of areas with significant changes
print(f"\nAreas with significant changes:")
print(f"Grid cells with >0.1m increase: {(diff > 0.1).sum().values}")
print(f"Grid cells with >0.1m decrease: {(diff < -0.1).sum().values}")
print(f"Grid cells with >0.5m increase: {(diff > 0.5).sum().values}")
print(f"Grid cells with >0.5m decrease: {(diff < -0.5).sum().values}")

# Close datasets
ds_cf0.close()
ds_cf8.close()

print(f"\nAnalysis complete! Check the saved PNG files in: {output_dir}")
print("Files created:")
print(f"  - {output_dir / 'sfincs_flood_depth_comparison.png'}")
print(f"  - {output_dir / 'sfincs_difference_only.png'}")
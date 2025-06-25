# %% In this script the regional SLR is calculate by using ISIMIP data from Treu et al. (2023): https://data.isimip.org/search/query/10.48364/ISIMIP.749905.1/
# Use the 'compass-snake-dfm' python environment
# Import the necessary packages
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import yaml
import ast
from shapely.geometry import box
import matplotlib.patches as mpatches
from scipy.stats import linregress
import pandas as pd

# Load config file
config_file = "../../../01_config_snakemake/config_general_MZB.yml"
with open(config_file, "r") as f:
    config = yaml.safe_load(f)
    config = config['runname_ids']['Idai']

# %%
# Region of case study area in Mozambique
lat_MZB = [-20.12, -19.30]
lon_MZB = [34.33, 34.95]

# DFM model region 
lon_min, lon_max, lat_min, lat_max = ast.literal_eval(config['bbox_dfm'])
lat_DFM = [lat_min, lat_max]
lon_DFM = [lon_min, lon_max]

# case study settings
start_date = np.datetime64('2019-03-09T00:00') 
end_date = np.datetime64('2019-03-24T00:00') 
     

# %%
# Reading in Factual (F) and Counterfactual (CF) data
SLR_hist_F  = xr.open_dataset(r"c:\Code\Paper_1\ISIMIP_SLR\hcc_obsclim_geocentricwaterlevel_global_hourly_2015.nc", engine='netcdf4')
SLR_hist_CF = xr.open_dataset(r"c:\Code\Paper_1\ISIMIP_SLR\hcc_counterclim_geocentricwaterlevel_global_hourly_2015.nc", engine='netcdf4')
                                                         
# %%
# Load lat and lon variables into memory
lat = SLR_hist_F['lat'].values
lon = SLR_hist_F['lon'].values

# Create a mask for stations within the region
mask_MZB = (lat >= lat_MZB[0]) & (lat <= lat_MZB[1]) & (lon >= lon_MZB[0]) & (lon <= lon_MZB[1])
mask_DFM = (lat >= lat_DFM[0]) & (lat <= lat_DFM[1]) & (lon >= lon_DFM[0]) & (lon <= lon_DFM[1])

# Apply the mask to filter stations within the region
stations_MZB_F  = SLR_hist_F.sel(stations=mask_MZB)
stations_DFM_F  = SLR_hist_F.sel(stations=mask_DFM)

## Do the same for the counterclim ##
# Load lat and lon variables into memory
lat = SLR_hist_CF['lat'].values
lon = SLR_hist_CF['lon'].values

# Create a mask for stations within the region
mask_MZB = (lat >= lat_MZB[0]) & (lat <= lat_MZB[1]) & (lon >= lon_MZB[0]) & (lon <= lon_MZB[1])
mask_DFM = (lat >= lat_DFM[0]) & (lat <= lat_DFM[1]) & (lon >= lon_DFM[0]) & (lon <= lon_DFM[1])

# Apply the mask to filter stations within the region
stations_MZB_CF  = SLR_hist_CF.sel(stations=mask_MZB)
stations_DFM_CF  = SLR_hist_CF.sel(stations=mask_DFM)

# %%
# Define the zoomed-in extent for the stations
# Adjust these limits as needed based on the locations of your stations
zoom_extent_MZB = [stations_MZB_F['lon'].min() - 1, stations_MZB_F['lon'].max() + 1,
                   stations_MZB_F['lat'].min() - 1, stations_MZB_F['lat'].max() + 1]
zoom_extent_DFM = [stations_DFM_F['lon'].min() - 1, stations_DFM_F['lon'].max() + 1,
                  stations_DFM_F['lat'].min() - 1, stations_DFM_F['lat'].max() + 1]

# Create a figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})

# List of station datasets, colors, titles, and zoom extents for each plot
stations_data = [stations_MZB_F, stations_DFM_F]
colors = ['red', 'orange']
titles = ['Zoomed Map with stations_MZB', 'Zoomed Map with stations_DFM']
zoom_extents = [zoom_extent_MZB, zoom_extent_DFM]

# Loop over both subplots to apply the same settings
for ax, stations, color, title, zoom_extent in zip(axes, stations_data, colors, titles, zoom_extents):
    ax.set_extent(zoom_extent, crs=ccrs.PlateCarree())  # Set the zoomed extent for each subplot
    ax.coastlines()                                    # Add coastlines
    ax.add_feature(cfeature.BORDERS)                   # Add country borders
    ax.add_feature(cfeature.LAND)                      # Add land feature
    ax.add_feature(cfeature.OCEAN)                     # Add ocean feature
    ax.gridlines(draw_labels=True)                     # Add gridlines with labels

    # Scatter plot for the stations
    ax.scatter(stations['lon'], stations['lat'], color=color, edgecolors='black', label=title, transform=ccrs.PlateCarree())

    # Set the title for each plot
    ax.set_title(title)

# %%
# Load lat and lon variables into memory
lon = stations_DFM_F['lon'].values

# Filter stations within the region and select 5 evenly spaced points
filtered_stations = stations_DFM_F.sel(stations=(lon >= lon.min()) & (lon <= lon.max()))
selected_stations_DFM_F = filtered_stations.isel(stations=np.round(np.linspace(0, len(filtered_stations['lat']) - 1, 5)).astype(int))

# Get coordinates for selected stations in the CF dataset
selected_lats, selected_lons = selected_stations_DFM_F['lat'].values, selected_stations_DFM_F['lon'].values
combined_mask = np.isin(stations_DFM_CF['lat'].values, selected_lats) & np.isin(stations_DFM_CF['lon'].values, selected_lons)
selected_stations_DFM_CF = stations_DFM_CF.sel(stations=combined_mask)

sel_stations_DFM = [selected_stations_DFM_F, selected_stations_DFM_CF]
# Define zoom extent for the map
zoom_extent = [32, 42, filtered_stations['lat'].min() - 1, filtered_stations['lat'].max() + 1]


# Create a figure with higher resolution and shared y-axis
fig, axs = plt.subplots(
    1, 2, 
    figsize=(8, 8), 
    dpi=300,  # ⬅️ increase resolution
    sharey=True, 
    subplot_kw={'projection': ccrs.PlateCarree()}
)

# Prepare the datasets and titles for looping
titles = ['Selected Stations F','Selected Stations CF']
colors = ['orange','blue']
bbox_polygon = box(lon_min, lat_min, lon_max, lat_max)
bbox_patch = mpatches.Patch(edgecolor='red', facecolor='none', linewidth=1, linestyle='--', label='DFM domain')

# Loop through each axis and corresponding dataset
for i, (ax, data, title, color) in enumerate(zip(axs, sel_stations_DFM, titles, colors)):
    ax.set_extent(zoom_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)

    ax.add_geometries(
        [bbox_polygon],
        crs=ccrs.PlateCarree(),
        edgecolor='red',
        facecolor='None',
        linewidth=1.5,
        linestyle='--'
    )
    
    # Setup gridlines with control over labels
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = i == 0        # only show y-axis on the first subplot
    gl.bottom_labels = True        # show x-axis on both; set to i == 1 if only on right

    ax.scatter(data['lon'], data['lat'], color=color, edgecolors='black', label=title, marker='o')

# Create combined legend on figure level
handles = [mpatches.Patch(color=c, label=t) for c, t in zip(colors, titles)] + [bbox_patch]
fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(1, 1.03), borderaxespad=0.)

fig.subplots_adjust(wspace=0.05)  # Smaller wspace = less space between columns

plt.tight_layout()
plt.show()

# %%
# To be sure, reindex one dataset to the other for the stations aound Beira (MZB), and the whole of Mozambique (MZ)
stations_MZB_CF_aligned = stations_MZB_CF['geocentricwaterlevel'].sel(time=stations_MZB_F['time'])
stations_DFM_CF_aligned = sel_stations_DFM[1]['geocentricwaterlevel'].sel(time=sel_stations_DFM[0]['time'])

# Subtract the CounterFactual (CF) from the Factual (CF) dataset for the stations near Beira
water_level_difference_MZB = stations_MZB_F['geocentricwaterlevel'] - stations_MZB_CF_aligned

# Subtract the CounterFactual (CF) from the Factual (CF) dataset for the five selected stations in the whole of Mozambique
water_level_difference_DFM  = sel_stations_DFM[0]['geocentricwaterlevel'] - stations_DFM_CF_aligned

# %%
# Calculate the mean and standard deviation across stations
mean_water_level_difference_MZB = water_level_difference_MZB.mean(dim='stations')
std_water_level_difference_MZB = water_level_difference_MZB.std(dim='stations')

mean_water_level_difference_DFM = water_level_difference_DFM.mean(dim='stations')
std_water_level_difference_DFM = water_level_difference_DFM.std(dim='stations')

# Create a figure with two subplots side-by-side
fig, axes = plt.subplots(2, 1, figsize=(8, 8))

# --- Plot 1: Mean Water Level Difference for MZB and MZ with Uncertainty Bounds ---
axes[0].plot(mean_water_level_difference_MZB['time'], mean_water_level_difference_MZB, label="Beira stations, MZ", color='orange')
axes[0].fill_between(mean_water_level_difference_MZB['time'],
                     mean_water_level_difference_MZB - std_water_level_difference_MZB,
                     mean_water_level_difference_MZB + std_water_level_difference_MZB,
                     color='orange', alpha=0.2)

axes[0].plot(mean_water_level_difference_DFM['time'], mean_water_level_difference_DFM, label="Five stations along MZ coast", color='blue')
axes[0].fill_between(mean_water_level_difference_DFM['time'],
                     mean_water_level_difference_DFM - std_water_level_difference_DFM,
                     mean_water_level_difference_DFM + std_water_level_difference_DFM,
                     color='blue', alpha=0.2)

# Customize the first plot
axes[0].set_title('Mean Difference of 2015 Geocentric Water Levels (Factual - Counterfactual)')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Mean Water Level Difference (mm)')
axes[0].grid()
axes[0].legend()

# --- Plot 2: Only Mean and Standard Deviation for MZB ---

axes[1].plot(mean_water_level_difference_MZB['time'], mean_water_level_difference_MZB, label="Beira stations, MZ", color='orange')
axes[1].fill_between(mean_water_level_difference_MZB['time'],
                     mean_water_level_difference_MZB - std_water_level_difference_MZB,
                     mean_water_level_difference_MZB + std_water_level_difference_MZB,
                     color='orange', alpha=0.2)

# Customize the second plot
axes[1].set_title('Mean Difference of 2015 Geocentric Water Levels for Beira, MZ (Factual - Counterfactual)')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('Mean Water Level Difference (mm)')
axes[1].grid()
axes[1].legend()

# Display the plots
plt.tight_layout()
plt.show()

# %%
# Extend and plot the timeseries until the time of the TC
# Extract 2015 time and values (convert to numeric for regression)
time_2015 = pd.to_datetime(mean_water_level_difference_DFM['time'].values)
x_2015 = (time_2015 - time_2015[0]).total_seconds() / (3600 * 24)  # days since start
y_2015 = mean_water_level_difference_MZB.values

# Fit linear trend
slope, intercept, r_value, p_value, std_err = linregress(x_2015, y_2015)

# Create new time range up to 2019
start_date = time_2015[0]
end_date = pd.Timestamp(end_date) # defined in the beginning
extended_time = pd.date_range(start=start_date, end=end_date, freq='D')
x_extended = (extended_time - start_date).days
y_extended = intercept + slope * x_extended

# Plot original and extrapolated
fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
ax.plot(time_2015, y_2015, label="2015 data", color='blue',linewidth=2)
ax.plot(extended_time, y_extended, label="Extrapolated trend to 2019", color='red', linestyle='--')

ax.set_title("Extrapolated Trend of 2015 Sea Level Rise for Five DFM stations along MZ coast")
ax.set_xlabel("Time")
ax.set_ylabel("Sea level rise (mm)")
ax.legend()
ax.grid()
plt.tight_layout()
plt.show()




# %%
#################################################################################
### Do the same for the long dataset to get the correct SLR for the Bathy ref ###
#################################################################################
# Load the 1990 & 2000 F and CF dataset to calculte the average SLR in this period
SLR_hist_F_1990 = xr.open_dataset(r"C:\Code\Paper_1\ISIMIP_SLR\hcc_obsclim_geocentricwaterlevel_global_hourly_1990.nc", engine='netcdf4')
SLR_hist_CF_1990 = xr.open_dataset(r"C:\Code\Paper_1\ISIMIP_SLR\hcc_counterclim_geocentricwaterlevel_global_hourly_1990.nc", engine='netcdf4')

SLR_hist_F_2000 = xr.open_dataset(r"C:\Code\Paper_1\ISIMIP_SLR\hcc_obsclim_geocentricwaterlevel_global_hourly_2000.nc", engine='netcdf4')
SLR_hist_CF_2000 = xr.open_dataset(r"C:\Code\Paper_1\ISIMIP_SLR\hcc_counterclim_geocentricwaterlevel_global_hourly_2000.nc", engine='netcdf4')

# Apply spatial mask to the long term monthly data
SLR_hist_F_1990 = SLR_hist_F_1990.sel(stations=mask_DFM)
SLR_hist_CF_1990 = SLR_hist_CF_1990.sel(stations=mask_DFM)
SLR_hist_F_2000 = SLR_hist_F_2000.sel(stations=mask_DFM)
SLR_hist_CF_2000 = SLR_hist_CF_2000.sel(stations=mask_DFM)

# Select only five stations in the region
sel_stations_F_1990  = SLR_hist_F_1990.sel(stations=combined_mask)
sel_stations_CF_1990 = SLR_hist_CF_1990.sel(stations=combined_mask)
sel_stations_F_2000  = SLR_hist_F_2000.sel(stations=combined_mask)
sel_stations_CF_2000 = SLR_hist_CF_2000.sel(stations=combined_mask)

sel_stations_DFM_1990 = [sel_stations_F_1990, sel_stations_CF_1990]
sel_stations_DFM_2000 = [sel_stations_F_2000, sel_stations_CF_2000]

# Create a figure with higher resolution and shared y-axis
fig, axs = plt.subplots(
    1, 2, 
    figsize=(8, 8), 
    dpi=300,  # ⬅️ increase resolution
    sharey=True, 
    subplot_kw={'projection': ccrs.PlateCarree()}
)

# Prepare the datasets and titles for looping
titles = ['Selected Stations F long','Selected Stations CF long']

# Loop through each axis and corresponding dataset
for i, (ax, data, title, color) in enumerate(zip(axs, sel_stations_DFM_2000, titles, colors)):
    ax.set_extent(zoom_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)

    ax.add_geometries(
        [bbox_polygon],
        crs=ccrs.PlateCarree(),
        edgecolor='red',
        facecolor='None',
        linewidth=1.5,
        linestyle='--'
    )
    
    # Setup gridlines with control over labels
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = i == 0        # only show y-axis on the first subplot
    gl.bottom_labels = True        # show x-axis on both; set to i == 1 if only on right

    ax.scatter(data['lon'], data['lat'], color=color, edgecolors='black', label=title, marker='o')

# Create combined legend on figure level
handles = [mpatches.Patch(color=c, label=t) for c, t in zip(colors, titles)] + [bbox_patch]
fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(1, 1.03), borderaxespad=0.)

fig.subplots_adjust(wspace=0.05)  # Smaller wspace = less space between columns

plt.tight_layout()
plt.show()


#%%
# To be sure, reindex one dataset to the other for the stations of the whole of Mozambique (MZ)
stations_DFM_CF_1990_aligned = sel_stations_DFM_1990[1]['geocentricwaterlevel'].sel(time=sel_stations_DFM_1990[0]['time'])
stations_DFM_CF_2000_aligned = sel_stations_DFM_2000[1]['geocentricwaterlevel'].sel(time=sel_stations_DFM_2000[0]['time'])

# Subtract the CounterFactual (CF) from the Factual (F) dataset for the five selected stations in the whole of Mozambique
water_level_difference_1990  = sel_stations_DFM_1990[0]['geocentricwaterlevel'] - stations_DFM_CF_1990_aligned
water_level_difference_2000  = sel_stations_DFM_2000[0]['geocentricwaterlevel'] - stations_DFM_CF_2000_aligned


# %%
# Calculate the mean and standard deviation across stations
mean_water_level_difference_1990 = water_level_difference_1990.mean(dim='stations')
std_water_level_difference_1990 = water_level_difference_1990.std(dim='stations')
mean_water_level_difference_2000 = water_level_difference_2000.mean(dim='stations')
std_water_level_difference_2000 = water_level_difference_2000.std(dim='stations')

mean_water_level_difference_1990 = water_level_difference_1990.mean(dim='stations')
std_water_level_difference_1990 = water_level_difference_1990.std(dim='stations')
mean_water_level_difference_2000 = water_level_difference_2000.mean(dim='stations')
std_water_level_difference_2000 = water_level_difference_2000.std(dim='stations')

#%%
bathy_mean_SLR_ref = mean_water_level_difference_1990.mean() + (mean_water_level_difference_2000.mean() - mean_water_level_difference_1990.mean())/2

#%%
# Create a figure with two subplots side-by-side
fig, axes = plt.subplots(1, 1, figsize=(8, 5))

# --- Plot 1: Mean Water Level Difference for MZB and MZ with Uncertainty Bounds ---
axes.plot(mean_water_level_difference_1990['time'], mean_water_level_difference_1990, label="SLR 1990 along five stations of MZ coast", color='orange')
axes.fill_between(mean_water_level_difference_1990['time'],
                     mean_water_level_difference_1990 - std_water_level_difference_1990,
                     mean_water_level_difference_1990 + std_water_level_difference_1990,
                     color='orange', alpha=0.2)

axes.plot(mean_water_level_difference_2000['time'], mean_water_level_difference_2000, label="SLR 2000 along five stations of MZ coast", color='blue')
axes.fill_between(mean_water_level_difference_2000['time'],
                     mean_water_level_difference_2000 - std_water_level_difference_2000,
                     mean_water_level_difference_2000 + std_water_level_difference_2000,
                     color='blue', alpha=0.2)

axes.axhline(y=bathy_mean_SLR_ref, color='red', linestyle='--', linewidth=1, label = "Mean SLR between 1990 & 2000")
# Annotate the line
axes.annotate(f'{bathy_mean_SLR_ref.values:.2f}', xy=((np.datetime64('2000-01-01T00:00')), bathy_mean_SLR_ref), xytext=((np.datetime64('2000-01-01T00:00')), bathy_mean_SLR_ref - 2.5),
            arrowprops=dict(arrowstyle='->'), fontsize=8)

# Customize the first plot
axes.set_title('Mean Difference of 2015 Geocentric Water Levels (Factual - Counterfactual)')
axes.set_xlabel('Time')
axes.set_ylabel('Mean Water Level Difference (mm)')
axes.grid()
axes.legend()


# %%

#%%
# Import the correct packages
import os
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature
import geopandas as gpd
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import xarray as xr
import netCDF4 as nc


#%%
# Region of case study area in Mozambique
lat_MZB = [-20.12, -19.30]
lon_MZB = [34.33, 34.95]

# region of eastern Africa
lat_MZ = [-27, -9]
lon_MZ = [29, 46]

# case study settings
start_date = np.datetime64('2019-03-09T00:00') 
end_date = np.datetime64('2019-03-24T00:00') 


# %%
# Reading in ERA5 the data from wflow folder
SLR_hist_F  = xr.open_dataset(r"z:\Code\Paper_1\ISIMIP_SLR\hcc_obsclim_geocentricwaterlevel_global_hourly_2015.nc", engine='netcdf4')
SLR_hist_CF = xr.open_dataset(r"z:\Code\Paper_1\ISIMIP_SLR\hcc_counterclim_geocentricwaterlevel_global_hourly_2015.nc", engine='netcdf4')


# In[4]:


# Load lat and lon variables into memory
lat = SLR_hist_F['lat'].values
lon = SLR_hist_F['lon'].values

# Create a mask for stations within the region
mask_MZB = (lat >= lat_MZB[0]) & (lat <= lat_MZB[1]) & (lon >= lon_MZB[0]) & (lon <= lon_MZB[1])
mask_MZ = (lat >= lat_MZ[0]) & (lat <= lat_MZ[1]) & (lon >= lon_MZ[0]) & (lon <= lon_MZ[1])

# Apply the mask to filter stations within the region
stations_MZB_F = SLR_hist_F.sel(stations=mask_MZB)
stations_MZ_F  = SLR_hist_F.sel(stations=mask_MZ)

## Do the same for the counterclim ##
# Load lat and lon variables into memory
lat = SLR_hist_CF['lat'].values
lon = SLR_hist_CF['lon'].values

# Create a mask for stations within the region
mask_MZB = (lat >= lat_MZB[0]) & (lat <= lat_MZB[1]) & (lon >= lon_MZB[0]) & (lon <= lon_MZB[1])
mask_MZ = (lat >= lat_MZ[0]) & (lat <= lat_MZ[1]) & (lon >= lon_MZ[0]) & (lon <= lon_MZ[1])

# Apply the mask to filter stations within the region
stations_MZB_CF = SLR_hist_CF.sel(stations=mask_MZB)
stations_MZ_CF  = SLR_hist_CF.sel(stations=mask_MZ)


# In[5]:


# Create a figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(15, 8), subplot_kw={'projection': ccrs.PlateCarree()})

# List of station datasets and colors for each plot
stations_data = [stations_MZB_F, stations_MZ_F]
colors = ['red', 'blue']
titles = ['Map with stations_MZB', 'Map with stations_MZ']

# Loop over both subplots to apply the same settings
for ax, stations, color, title in zip(axes, stations_data, colors, titles):
    ax.set_global()                        # Set global extent for the plot
    ax.coastlines()                        # Add coastlines
    ax.add_feature(cfeature.BORDERS)       # Add country borders
    ax.add_feature(cfeature.LAND)          # Add land feature
    ax.add_feature(cfeature.OCEAN)         # Add ocean feature
    ax.gridlines(draw_labels=True)         # Add gridlines with labels

    # Scatter plot for the stations
    ax.scatter(stations['lon'], stations['lat'], color=color, label=title, transform=ccrs.PlateCarree(), s=8)

    # Set the title for each plot
    ax.set_title(title)

# Display the figure
plt.show()


# In[6]:


# Define the zoomed-in extent for the stations
# Adjust these limits as needed based on the locations of your stations
zoom_extent_MZB = [stations_MZB_F['lon'].min() - 1, stations_MZB_F['lon'].max() + 1,
                   stations_MZB_F['lat'].min() - 1, stations_MZB_F['lat'].max() + 1]
zoom_extent_MZ = [stations_MZ_F['lon'].min() - 1, stations_MZ_F['lon'].max() + 1,
                  stations_MZ_F['lat'].min() - 1, stations_MZ_F['lat'].max() + 1]

# Create a figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(15, 8), subplot_kw={'projection': ccrs.PlateCarree()})

# List of station datasets, colors, titles, and zoom extents for each plot
stations_data = [stations_MZB_F, stations_MZ_F]
colors = ['red', 'blue']
titles = ['Zoomed Map with stations_MZB', 'Zoomed Map with stations_MZ']
zoom_extents = [zoom_extent_MZB, zoom_extent_MZ]

# Loop over both subplots to apply the same settings
for ax, stations, color, title, zoom_extent in zip(axes, stations_data, colors, titles, zoom_extents):
    ax.set_extent(zoom_extent, crs=ccrs.PlateCarree())  # Set the zoomed extent for each subplot
    ax.coastlines()                                    # Add coastlines
    ax.add_feature(cfeature.BORDERS)                   # Add country borders
    ax.add_feature(cfeature.LAND)                      # Add land feature
    ax.add_feature(cfeature.OCEAN)                     # Add ocean feature
    ax.gridlines(draw_labels=True)                     # Add gridlines with labels

    # Scatter plot for the stations
    ax.scatter(stations['lon'], stations['lat'], color=color, label=title, transform=ccrs.PlateCarree())

    # Set the title for each plot
    ax.set_title(title)

# Display the figure
plt.show()


# In[17]:


# Load lat and lon variables into memory
lon = stations_MZ_F['lon'].values

# Filter stations within the region and select 5 evenly spaced points
filtered_stations = stations_MZ_F.sel(stations=(lon >= 32) & (lon <= 42))
selected_stations_MZ_F = filtered_stations.isel(stations=np.round(np.linspace(0, len(filtered_stations['lat']) - 1, 5)).astype(int))

# Get coordinates for selected stations in the CF dataset
selected_lats, selected_lons = selected_stations_MZ_F['lat'].values, selected_stations_MZ_F['lon'].values
combined_mask = np.isin(stations_MZ_CF['lat'].values, selected_lats) & np.isin(stations_MZ_CF['lon'].values, selected_lons)
selected_stations_MZ_CF = stations_MZ_CF.sel(stations=combined_mask)

# Define zoom extent for the map
zoom_extent = [32, 42, filtered_stations['lat'].min() - 1, filtered_stations['lat'].max() + 1]

# Create a figure with subplots for both datasets
fig, axs = plt.subplots(1, 2, figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()}, tight_layout=True)

# Prepare the datasets and titles for looping
MZ_sel_datasets = [selected_stations_MZ_F, selected_stations_MZ_CF]
titles = ['Selected Stations F', 'Selected Stations CF']
colors = ['red', 'blue']

# Loop through each axis and corresponding dataset
for ax, data, title, color in zip(axs, MZ_sel_datasets, titles, colors):
    ax.set_extent(zoom_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.gridlines(draw_labels=True)
    ax.scatter(data['lon'], data['lat'], color=color, label=f'{title}', marker='o')
    ax.set_title(title, fontsize=14)
    ax.legend()

# Adjust the layout for a tighter fit
plt.tight_layout(pad=3.0)  # Adjust pad as needed for spacing
plt.show()


# In[8]:


# To be sure, reindex one dataset to the other for the stations aound Beira (MZB), and the whole of Mozambique (MZ)
stations_MZB_CF_aligned = stations_MZB_CF['geocentricwaterlevel'].sel(time=stations_MZB_F['time'])
stations_MZ_CF_aligned = MZ_sel_datasets[1]['geocentricwaterlevel'].sel(time=MZ_sel_datasets[0]['time'])

# Subtract the CounterFactual (CF) from the Factual (CF) dataset for the stations near Beira
water_level_difference_MZB = stations_MZB_F['geocentricwaterlevel'] - stations_MZB_CF_aligned

# Subtract the CounterFactual (CF) from the Factual (CF) dataset for the five selected stations in the whole of Mozambique
water_level_difference_MZ  = MZ_sel_datasets[0]['geocentricwaterlevel'] - stations_MZ_CF_aligned


# In[9]:


# Calculate the mean and standard deviation across stations
mean_water_level_difference_MZB = water_level_difference_MZB.mean(dim='stations')
std_water_level_difference_MZB = water_level_difference_MZB.std(dim='stations')

mean_water_level_difference_MZ = water_level_difference_MZ.mean(dim='stations')
std_water_level_difference_MZ = water_level_difference_MZ.std(dim='stations')

# Create a figure with two subplots side-by-side
fig, axes = plt.subplots(1, 2, figsize=(18, 6))

# --- Plot 1: Mean Water Level Difference for MZB and MZ with Uncertainty Bounds ---
axes[0].plot(mean_water_level_difference_MZB['time'], mean_water_level_difference_MZB, label="Beira stations, MZ", color='blue')
axes[0].fill_between(mean_water_level_difference_MZB['time'],
                     mean_water_level_difference_MZB - std_water_level_difference_MZB,
                     mean_water_level_difference_MZB + std_water_level_difference_MZB,
                     color='blue', alpha=0.2)

axes[0].plot(mean_water_level_difference_MZ['time'], mean_water_level_difference_MZ, label="Five stations along MZ coast", color='orange')
axes[0].fill_between(mean_water_level_difference_MZ['time'],
                     mean_water_level_difference_MZ - std_water_level_difference_MZ,
                     mean_water_level_difference_MZ + std_water_level_difference_MZ,
                     color='orange', alpha=0.2)

# Customize the first plot
axes[0].set_title('Mean Difference of 2015 Geocentric Water Levels (Factual - Counterfactual)')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Mean Water Level Difference (mm)')
axes[0].grid()
axes[0].legend()

# --- Plot 2: Only Mean and Standard Deviation for MZB ---

axes[1].plot(mean_water_level_difference_MZB['time'], mean_water_level_difference_MZB, label="Beira stations, MZ", color='blue')
axes[1].fill_between(mean_water_level_difference_MZB['time'],
                     mean_water_level_difference_MZB - std_water_level_difference_MZB,
                     mean_water_level_difference_MZB + std_water_level_difference_MZB,
                     color='blue', alpha=0.2)

# Customize the second plot
axes[1].set_title('Mean Difference of 2015 Geocentric Water Levels for Beira, MZ (Factual - Counterfactual)')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('Mean Water Level Difference (mm)')
axes[1].grid()
axes[1].legend()

# Display the plots
plt.tight_layout()
plt.show()


# In[10]:


# Substract the two datasets for the five selected stations inthe whole of  Mozambique
# water_level_difference_MZ  = MZ_sel_datasets[0]['geocentricwaterlevel'] - MZ_sel_datasets[1]['geocentricwaterlevel']

# Calculate the mean difference over the specified time range
mean_water_level_difference_MZB = water_level_difference_MZB.mean(dim=['stations'])
mean_water_level_difference_MZ  = water_level_difference_MZ.mean(dim=['stations'])

# Plot the mean result
plt.figure(figsize=(12, 6))
mean_water_level_difference_MZB.plot() 
mean_water_level_difference_MZ.plot() 

plt.title('Mean Difference of 2015 Geocentric Water Levels for Beira, MZ (Factual - Counterfactual)')
plt.xlabel('Time')
plt.ylabel('Mean Water Level Difference (mm)')
plt.grid()
plt.show()


# In[11]:


# Plot the mean result
plt.figure(figsize=(12, 6))

# Plot each dataset with labels
mean_water_level_difference_MZB.plot(label="Beira stations, MZ")
mean_water_level_difference_MZ.plot(label="five stations along MZ coast")

# Set the title and axis labels
plt.title('Mean Difference of 2015 Geocentric Water Levels for Beira, MZ (Factual - Counterfactual)')
plt.xlabel('Time')
plt.ylabel('Mean Water Level Difference (mm)')

# Add grid and legend
plt.grid()
plt.legend()  # Display the legend with specified labels
plt.show()


# In[ ]:


# Ensure you have your aligned datasets for geocentric water levels
stations_MZB_CF_aligned = stations_MZB_CF['geocentricwaterlevel'].sel(time=stations_MZB_F['time'])

# Subtract the two datasets to get the difference
water_level_difference = stations_MZB_F['geocentricwaterlevel'] - stations_MZB_CF_aligned

# Calculate the mean difference over the specified time range
mean_water_level_difference = water_level_difference.mean(dim='stations')

# Set up the plot
plt.figure(figsize=(12, 6))

# Plot the mean water level difference
mean_water_level_difference.plot(label='Mean Water Level Difference', color='blue')  # Using Xarray's built-in plotting

# Add titles and labels
plt.title('Mean Difference of Geocentric Water Levels for Beira, MZ (Factual - Counterfactual)')
plt.xlabel('Time')
plt.ylabel('Mean Water Level Difference (mm)')
plt.grid()
plt.legend()

# Plot differences per station with correct time axis
for station in water_level_difference.stations.values:
    # Extract the time and the difference data for the specific station
    station_data = water_level_difference.sel(stations=station)
    
    # Plot with the correct time axis
    plt.plot(station_data['time'], station_data, label=f'Station {station}', alpha=0.5)  # Alpha for transparency

# Optionally, adjust the legend to avoid overlap
plt.legend(loc='upper right', bbox_to_anchor=(1.2, 1), fontsize='small')

plt.tight_layout()  # Adjust layout for better spacing
plt.show()


# In[ ]:





# In[ ]:





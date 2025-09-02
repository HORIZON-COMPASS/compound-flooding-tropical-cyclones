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
from datetime import datetime
from matplotlib.patches import Patch
import matplotlib.dates as mdates

# Load config file
config_file = "../../Workflows/01_config_snakemake/config_general_MZB.yml"
with open(config_file, "r") as f:
    config = yaml.safe_load(f)
    config = config['runname_ids']['Idai']

# %%
# DFM model region 
lon_min, lon_max, lat_min, lat_max = ast.literal_eval(config['bbox_dfm'])
lat_DFM = [lat_min, lat_max]
lon_DFM = [lon_min, lon_max]

# case study settings
start_date = np.datetime64(datetime.strptime(config['start_time'], "%Y%m%d %H%M%S"))
end_date = np.datetime64(datetime.strptime(config['end_time'], "%Y%m%d %H%M%S"))

# %%
# Reading in Factual (F) and Counterfactual (CF) data
SLR_hist_F  = xr.open_dataset(r"c:\Code\Paper_1\ISIMIP_SLR\hcc_obsclim_geocentricwaterlevel_global_hourly_2015.nc", engine='netcdf4')
SLR_hist_CF = xr.open_dataset(r"c:\Code\Paper_1\ISIMIP_SLR\hcc_counterclim_geocentricwaterlevel_global_hourly_2015.nc", engine='netcdf4')
                                                         
# %%
# Load lat and lon variables into memory
lat = SLR_hist_F['lat'].values
lon = SLR_hist_F['lon'].values

# Create a mask for stations within the DFM model domain
mask_DFM = (lat >= lat_DFM[0]) & (lat <= lat_DFM[1]) & (lon >= lon_DFM[0]) & (lon <= lon_DFM[1])
# Apply the mask to filter stations within the region
stations_DFM_F  = SLR_hist_F.sel(stations=mask_DFM)

## Do the same for the counterclim ##
lat = SLR_hist_CF['lat'].values
lon = SLR_hist_CF['lon'].values

# Create a mask for stations within the region
mask_DFM = (lat >= lat_DFM[0]) & (lat <= lat_DFM[1]) & (lon >= lon_DFM[0]) & (lon <= lon_DFM[1])
# Apply the mask to filter stations within the region
stations_DFM_CF  = SLR_hist_CF.sel(stations=mask_DFM)


# %%
# Filter stations within the region and select 5 evenly spaced points
selected_stations_DFM_F = stations_DFM_F.isel(stations=np.round(np.linspace(0, len(stations_DFM_F['lat']) - 1, 5)).astype(int))

# Get coordinates for selected stations
selected_lats, selected_lons = selected_stations_DFM_F['lat'].values, selected_stations_DFM_F['lon'].values

# To be sure, reindex one dataset to the other
combined_mask = np.isin(stations_DFM_CF['lat'].values, selected_lats) & np.isin(stations_DFM_CF['lon'].values, selected_lons)
selected_stations_DFM_CF = stations_DFM_CF.sel(stations=combined_mask)
stations_DFM_CF_aligned = selected_stations_DFM_CF['geocentricwaterlevel'].sel(time=selected_stations_DFM_F['time'])

# Calculate the water level difference between factual and counterfactual datasets
water_level_difference_DFM  = selected_stations_DFM_F['geocentricwaterlevel'] - stations_DFM_CF_aligned

# Average over the 5 selected stations
mean_water_level_difference_DFM = water_level_difference_DFM.mean(dim='stations')
std_water_level_difference_DFM = water_level_difference_DFM.std(dim='stations')


#%%
# Plot all ISIMIP stations in the DFM region and select five evenly spaced stations along the coast
sel_stations_DFM = [stations_DFM_F, selected_stations_DFM_F]

# Define zoom extent for the map
zoom_extent = [32, 42, stations_DFM_F['lat'].min() - 1, stations_DFM_F['lat'].max() + 1]

# Create a figure with higher resolution and shared y-axis
fig, axs = plt.subplots(
    1, 2, 
    figsize=(8, 8), 
    dpi=300,  # ⬅️ increase resolution
    sharey=True, 
    subplot_kw={'projection': ccrs.PlateCarree()},
    constrained_layout=True
)

# Prepare the datasets and titles for looping
titles = ['All ISIMIP stations within DFM domain','Selected Stations']
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
fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(1, 1.04), borderaxespad=0.)

# Subplot (a) - first plot
axs[0].text(0.0, 1.03, "(a)", transform=axs[0].transAxes,
             fontsize=14, fontweight='bold', va='top', ha='left')

# Subplot (b) - second plot
axs[1].text(0.0, 1.03, "(b)", transform=axs[1].transAxes,
             fontsize=14, fontweight='bold', va='top', ha='left')

fig.savefig("../figures/fS8.png",dpi=300, bbox_inches='tight')
fig.savefig("../figures/fS8.pdf",dpi=300, bbox_inches='tight')
plt.show()


# %%
# Extend and plot the timeseries until the time of the TC, and the std deviation between the 5 selected stations
# --- Prepare data ---
time_2015 = pd.to_datetime(mean_water_level_difference_DFM['time'].values)
x_2015 = (time_2015 - time_2015[0]).total_seconds() / (3600 * 24)  # days since start
y_2015 = mean_water_level_difference_DFM.values

# Fit linear trend
slope, intercept, r_value, p_value, std_err = linregress(x_2015, y_2015)

# Extend time to 2019
start_date = time_2015[0]
extended_time = pd.date_range(start=start_date, end=pd.Timestamp(end_date), freq='D')
x_extended = (extended_time - start_date).days
y_extended = intercept + slope * x_extended

# --- Create figure ---
fig, axes = plt.subplots(
    1, 2, figsize=(12, 5), dpi=300, sharey=True,
    gridspec_kw={'width_ratios': [1, 1.3], 'wspace': 0},  # no horizontal space
    constrained_layout=True
)

# --- Plot 1: Mean water level difference with uncertainty ---
axes[0].plot(time_2015, y_2015, label="Mean of 5 stations along MZ coast", color='blue', linewidth=2)
axes[0].fill_between(time_2015,
                     y_2015 - std_water_level_difference_DFM,
                     y_2015 + std_water_level_difference_DFM,
                     color='blue', alpha=0.2)

axes[0].set_title('Mean Difference of 2015 Geocentric Water Levels \n(Factual - Counterfactual)', fontsize=12)
axes[0].set_xlabel('Month (2015)', fontsize=11)
axes[0].set_ylabel("Sea level rise (mm)", fontsize=11)
axes[0].grid(True)
axes[0].set_xlim([time_2015[0], time_2015[-1]])
axes[0].xaxis.set_major_locator(mdates.MonthLocator())
axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
std_patch = Patch(facecolor='blue', alpha=0.2, label='±1 Std Dev')
axes[0].legend(handles=[axes[0].lines[0], std_patch], fontsize=10)

# --- Plot 2: Original 2015 data and extrapolated trend ---
axes[1].plot(time_2015, y_2015, label="Mean of 5 stations along MZ coast", color='blue', linewidth=2)
print_time = end_date.strftime("%Y-%m-%d")
axes[1].plot(extended_time, y_extended, label=f"Extrapolated trend to {print_time}", color='red', linestyle='--', linewidth=2)

axes[1].set_title("Extrapolated Trend of 2015 Sea Level Rise", fontsize=12)
axes[1].set_xlabel("Time", fontsize=11)
axes[1].set_xlim([time_2015[0], extended_time[-1]])
axes[1].grid(True)
axes[1].legend(fontsize=10)

# Subplot (a) - first plot
axes[0].text(0.0, 1.06, "(a)", transform=axes[0].transAxes,
             fontsize=14, fontweight='bold', va='top', ha='left')

# Subplot (b) - second plot
axes[1].text(0.0, 1.06, "(b)", transform=axes[1].transAxes,
             fontsize=14, fontweight='bold', va='top', ha='left')

fig.savefig("../figures/fS9.png",dpi=300, bbox_inches='tight')
fig.savefig("../figures/fS9.pdf",dpi=300, bbox_inches='tight')
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
combined_mask = np.isin(stations_DFM_F['lat'].values, selected_lats) & np.isin(stations_DFM_F['lon'].values, selected_lons)
selected_stations_DFM_F = stations_DFM_F.sel(stations=combined_mask)
sel_stations_F_1990  = SLR_hist_F_1990.sel(stations=combined_mask)
sel_stations_CF_1990 = SLR_hist_CF_1990.sel(stations=combined_mask)
sel_stations_F_2000  = SLR_hist_F_2000.sel(stations=combined_mask)
sel_stations_CF_2000 = SLR_hist_CF_2000.sel(stations=combined_mask)

sel_stations_DFM_1990 = [sel_stations_F_1990, sel_stations_CF_1990]
sel_stations_DFM_2000 = [sel_stations_F_2000, sel_stations_CF_2000]

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
            arrowprops=dict(arrowstyle='->'), fontsize=10)

# Customize the first plot
axes.set_title('Mean Difference of 2015 Geocentric Water Levels (Factual - Counterfactual)', fontsize=12)
axes.set_xlabel('Time', fontsize=11)
axes.set_ylabel('Mean Water Level Difference (mm)', fontsize=11)
axes.grid()
axes.legend()
axes.set_xlim([mean_water_level_difference_1990['time'][0], mean_water_level_difference_2000['time'][-1]])
axes.xaxis.set_major_locator(mdates.YearLocator(2))

fig.savefig("../figures/fS10.png",dpi=300, bbox_inches='tight')
fig.savefig("../figures/fS10.pdf",dpi=300, bbox_inches='tight')

# %%

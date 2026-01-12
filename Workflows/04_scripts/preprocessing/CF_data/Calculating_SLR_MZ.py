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
import platform
import os
import statsmodels.api as sm
from scipy.interpolate import interp1d
import requests
from io import StringIO

prefix = "p:/" if platform.system() == "Windows" else "/p/"

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
start_time = config["start_time"]  # "20190309 000000"
end_time = config["end_time"]  # "20190325 000000"
start_time_iso = f"{start_time[:4]}-{start_time[4:6]}-{start_time[6:8]}T{start_time[9:11]}:{start_time[11:13]}:{start_time[13:15]}"
end_time_iso = f"{end_time[:4]}-{end_time[4:6]}-{end_time[6:8]}T{end_time[9:11]}:{end_time[11:13]}:{end_time[13:15]}"
start_date_event = np.datetime64(start_time_iso) 
end_date_event = np.datetime64(end_time_iso) 

# %%
# Reading in Factual (F) and Counterfactual (CF) data
all_year_dir = os.path.join(prefix, "11210471-001-compass/01_Data/ISIMIP/SLR/DFM_MZ_subset/MZB_combined")

SL_hist_F_gc = xr.open_dataset(os.path.join(all_year_dir, "hcc_obsclim_geocentricwaterlevel_global_hourly_all_years.nc"), engine='netcdf4')
SL_hist_CF_gc = xr.open_dataset(os.path.join(all_year_dir, "hcc_counterclim_geocentricwaterlevel_global_hourly_all_years.nc"), engine='netcdf4')
SLR_hist_gc = SL_hist_F_gc['geocentricwaterlevel'] - SL_hist_CF_gc['geocentricwaterlevel']

SL_hist_F_wl  = xr.open_dataset(os.path.join(all_year_dir, "hcc_obsclim_waterlevel_global_hourly_all_years.nc"), engine='netcdf4')
SL_hist_CF_wl = xr.open_dataset(os.path.join(all_year_dir, "hcc_counterclim_waterlevel_global_hourly_all_years.nc"), engine='netcdf4')
SLR_hist_wl = SL_hist_F_wl['waterlevel'] - SL_hist_CF_wl['waterlevel']

# %%
# Some helper functions
def select_stations(ds, lat_bounds, lon_bounds):
    lat = ds["lat"].values
    lon = ds["lon"].values
    mask = (
        (lat >= lat_bounds[0]) & (lat <= lat_bounds[1]) &
        (lon >= lon_bounds[0]) & (lon <= lon_bounds[1])
    )
    return ds.sel(stations=mask)

def select_evenly_spaced(ds, n=5):
    """Select n evenly spaced stations based on index order."""
    idx = np.round(np.linspace(0, ds.sizes["stations"] - 1, n)).astype(int)
    return ds.isel(stations=idx)

def match_stations_by_coords(ds_ref, ds_other):
    """Select stations in ds_other that match lat/lon of ds_ref."""
    lat, lon = ds_ref["lat"].values, ds_ref["lon"].values
    mask = np.isin(ds_other["lat"].values, lat) & np.isin(ds_other["lon"].values, lon)
    return ds_other.sel(stations=mask)

def align_cf_to_f(f_ds, cf_ds, var):
    """Align CF dataset to F time axis and subtract."""
    cf_aligned = cf_ds[var].sel(time=f_ds["time"])
    return f_ds[var] - cf_aligned

# ------------------------------------------------------------------
### Select stations in the MZB and DFM regions ###
# Geocentric water levels
SLR_MZB_gc = select_stations(SLR_hist_gc,  lat_MZB, lon_MZB)
SLR_DFM_gc = select_stations(SLR_hist_gc,  lat_DFM, lon_DFM)

# Non-geocentric water levels
SLR_MZB_wl = select_stations(SLR_hist_wl,  lat_MZB, lon_MZB)
SLR_DFM_wl = select_stations(SLR_hist_wl,  lat_DFM, lon_DFM)

# ------------------------------------------------------------------
# Stats across stations
mean_SLR_MZB_gc, std_SLR_MZB_gc = SLR_MZB_gc.mean("stations"), SLR_MZB_gc.std("stations")
mean_SLR_DFM_gc, std_SLR_DFM_gc = SLR_DFM_gc.mean("stations"), SLR_DFM_gc.std("stations")

mean_SLR_MZB_wl, std_SLR_MZB_wl = SLR_MZB_wl.mean("stations"), SLR_MZB_wl.std("stations")
mean_SLR_DFM_wl, std_SLR_DFM_wl = SLR_DFM_wl.mean("stations"), SLR_DFM_wl.std("stations")


# %%
# Extend and plot the timeseries until the time of the TC
time_DFM = pd.to_datetime(mean_SLR_DFM_gc['time'].values)
x_DFM = (time_DFM - time_DFM[0]).total_seconds() / (3600 * 24)  # days since start
y_DFM_gc = mean_SLR_DFM_gc.values
y_DFM_wl = mean_SLR_DFM_wl.values

# --- Linear Trendline ---
slope_gc, intercept_gc, r_value_gc, p_value_gc, std_err_gc = linregress(x_DFM, y_DFM_gc)
slope_wl, intercept_wl, r_value_wl, p_value_wl, std_err_wl = linregress(x_DFM, y_DFM_wl)

# --- LOWESS Trendline ---
y_series_gc = pd.Series(y_DFM_gc, index=time_DFM)
y_daily_gc = y_series_gc.resample("1D").mean()
x_daily_gc = (y_daily_gc.index - time_DFM[0]).days

lowess_daily_gc = sm.nonparametric.lowess(y_daily_gc.values, x_daily_gc, frac=0.3, return_sorted=False)
interp_lowess_gc = interp1d(x_daily_gc, lowess_daily_gc, kind="linear", fill_value="extrapolate")

# Same for water level only
y_series_wl = pd.Series(y_DFM_wl, index=time_DFM)
y_daily_wl = y_series_wl.resample("1D").mean()
x_daily_wl = (y_daily_wl.index - time_DFM[0]).days

lowess_daily_wl = sm.nonparametric.lowess(y_daily_wl.values, x_daily_wl, frac=0.3, return_sorted=False)
interp_lowess_wl = interp1d(x_daily_wl, lowess_daily_wl, kind="linear", fill_value="extrapolate")



# Create new time range up to 2019
start_date = time_DFM[0]
end_date = pd.Timestamp(end_date_event) # defined in the beginning
extended_time = pd.date_range(start=start_date, end=end_date, freq='D')

# Create DataFrame for extended time and calculated trend values
df_extended = pd.DataFrame({'time': extended_time})
x_extended = (extended_time - start_date).days  # days since start
df_extended['linear_gc'] = intercept_gc + slope_gc * x_extended
df_extended['lowess_gc'] = interp_lowess_gc(x_extended)
df_extended['linear_wl'] = intercept_wl + slope_wl * x_extended
df_extended['lowess_wl'] = interp_lowess_wl(x_extended)


#%%###################################################################################
#################################### PLOTS ###########################################
######################################################################################

############################### PLOT LONG-TERM TRENDS ################################
fig, axes = plt.subplots(2, 1, figsize=(8, 8))

# --- Plot 1: Mean Water Level Difference for MZB and MZ with Uncertainty Bounds ---
axes[0].plot(mean_SLR_MZB_gc['time'], mean_SLR_MZB_gc, label="Beira stations (gc)", color='orange')
axes[0].fill_between(mean_SLR_MZB_gc['time'],
                     mean_SLR_MZB_gc - std_SLR_MZB_gc,
                     mean_SLR_MZB_gc + std_SLR_MZB_gc,
                     color='orange', alpha=0.2)

axes[0].plot(mean_SLR_DFM_gc['time'], mean_SLR_DFM_gc, label="DFM stations (gc)", color='red')
axes[0].fill_between(mean_SLR_DFM_gc['time'],
                     mean_SLR_DFM_gc - std_SLR_DFM_gc,
                     mean_SLR_DFM_gc + std_SLR_DFM_gc,
                     color='red', alpha=0.2)

# same for waterlevel only
axes[0].plot(mean_SLR_MZB_wl['time'], mean_SLR_MZB_wl, label="Beira stations (wl)", color='lightblue')
axes[0].fill_between(mean_SLR_MZB_wl['time'],
                     mean_SLR_MZB_wl - std_SLR_MZB_wl,
                     mean_SLR_MZB_wl + std_SLR_MZB_wl,
                     color='lightblue', alpha=0.2)

axes[0].plot(mean_SLR_DFM_wl['time'], mean_SLR_DFM_wl, label="DFM stations (wl)", color='blue')
axes[0].fill_between(mean_SLR_DFM_wl['time'],
                     mean_SLR_DFM_wl - std_SLR_DFM_wl,
                     mean_SLR_DFM_wl + std_SLR_DFM_wl,
                     color='blue', alpha=0.2)

# Customize the first plot
axes[0].set_title('Mean SLR for Beira and DFM stations with and without VLM (Factual - Counterfactual)')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Mean SLR (mm)')
axes[0].grid()
axes[0].legend()

# --- Plot 2: Only Mean and Standard Deviation for MZB ---
axes[1].plot(mean_SLR_MZB_gc['time'], mean_SLR_MZB_gc, label="Beira stations (gc)", color='orange')
axes[1].fill_between(mean_SLR_MZB_gc['time'],
                     mean_SLR_MZB_gc - std_SLR_MZB_gc,
                     mean_SLR_MZB_gc + std_SLR_MZB_gc,
                     color='orange', alpha=0.2)

# same for waterlevel only
axes[1].plot(mean_SLR_MZB_wl['time'], mean_SLR_MZB_wl, label="Beira stations (wl)", color='lightblue')
axes[1].fill_between(mean_SLR_MZB_wl['time'],
                     mean_SLR_MZB_wl - std_SLR_MZB_wl,
                     mean_SLR_MZB_wl + std_SLR_MZB_wl,
                     color='lightblue', alpha=0.2)

# Customize the second plot
axes[1].set_title('Mean Difference 1985-2015 Water Levels for Beira, MZ (Factual - Counterfactual)')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('Mean Water Level Difference (mm)')
axes[1].grid()
axes[1].legend()

# Display the plots
plt.tight_layout()
plt.show()

# %%
############################### PLOT EXTRAPOLATED TRENDS ################################
# Plot original and extrapolated
fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
ax.plot(time_DFM, y_DFM_gc, label="ISIMIP SLR data (gc)", color='orange',linewidth=2)
ax.fill_between(mean_SLR_DFM_gc['time'],
                     mean_SLR_DFM_gc - std_SLR_DFM_gc,
                     mean_SLR_DFM_gc + std_SLR_DFM_gc,
                     color='orange', alpha=0.2)
ax.plot(extended_time, df_extended['linear_gc'], label=f"Extrapolated linear trend (gc)", color='grey', linestyle='--')
ax.plot(extended_time, df_extended['lowess_gc'], label=f"Extrapolated LOWESS trend (gc)", color='black', linestyle='--')

ax.plot(time_DFM, y_DFM_wl, label="ISIMIP SLR data (wl)", color='blue',linewidth=2)
ax.fill_between(mean_SLR_DFM_wl['time'],
                     mean_SLR_DFM_wl - std_SLR_DFM_wl,
                     mean_SLR_DFM_wl + std_SLR_DFM_wl,
                     color='blue', alpha=0.2)
ax.plot(extended_time, df_extended['linear_wl'], label=f"Extrapolated linear trend (wl)", color='grey', linestyle='--')
ax.plot(extended_time, df_extended['lowess_wl'], label=f"Extrapolated LOWESS trend (wl)", color='black', linestyle='--')

ax.set_title(f"Extrapolated ISIMIP SLR trend to {end_date.date()} for stations within DFM domain")
ax.set_xlabel("Time")
ax.set_ylabel("Sea level rise (mm)")
ax.legend(fontsize=8)
ax.grid()
plt.tight_layout()
plt.show()


# %%
############################### PLOT STATIONS SPATIALLY ################################
# Define the zoomed-in extent for the stations
zoom_extent = [32, 42, lat_DFM[0] - 1, lat_DFM[1] + 1]

# Create a figure with two subplots
fig, ax = plt.subplots(1, 1, figsize=(8, 8), subplot_kw={'projection': ccrs.PlateCarree()})

# List of station datasets, colors, titles, and zoom extents for each plot
stations_data = [SLR_MZB_wl, SLR_DFM_wl]
colors = ['red', 'orange']
bbox_polygon = box(lon_min, lat_min, lon_max, lat_max)
bbox_patch = mpatches.Patch(edgecolor='red', facecolor='none', linewidth=1, linestyle='--', label='DFM domain')

# Loop over both subplots to apply the same settings
ax.set_extent(zoom_extent, crs=ccrs.PlateCarree())  # Set the zoomed extent for each subplot
ax.coastlines()                                    # Add coastlines
ax.add_feature(cfeature.BORDERS)                   # Add country borders
ax.add_feature(cfeature.LAND)                      # Add land feature
ax.add_feature(cfeature.OCEAN)                     # Add ocean feature
ax.gridlines(draw_labels=True)                     # Add gridlines with labels

ax.add_geometries(
        [bbox_polygon],
        crs=ccrs.PlateCarree(),
        edgecolor='red',
        facecolor='None',
        linewidth=1.5,
        linestyle='--'
    )
# Scatter plot for the stations
ax.scatter(SLR_DFM_wl['lon'], SLR_DFM_wl['lat'], color='orange', edgecolors='black', label="DFM stations", transform=ccrs.PlateCarree())
ax.scatter(SLR_MZB_wl['lon'], SLR_MZB_wl['lat'], color='red', edgecolors='black', label="Beira stations", transform=ccrs.PlateCarree())

# Create combined legend on figure level
titles = ['Beira stations','DFM stations']
handles = [mpatches.Patch(color=c, label=t) for c, t in zip(colors, titles)] + [bbox_patch]
fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(1, 1.03), borderaxespad=0.)

# Set the title for each plot
ax.set_title("Zoomed map with ISIMIP stations")


# %%
#################################################################################
###################### Calculate the SLR for the Bathy ref ######################
#################################################################################
# use lowess fit and calculate mean SLR between 1990 and 2000
bathy_mean_SLR_ref_wl = df_extended[df_extended["time"].dt.year.between(1990, 2000)]['lowess_wl'].mean()
bathy_mean_SLR_ref_gc = df_extended[df_extended["time"].dt.year.between(1990, 2000)]['lowess_gc'].mean()
mean_time = df_extended[df_extended["time"].dt.year.between(1990, 2000)]['time'].mean()


#%%
# Plot the mean SLR between 1990 and 2000 and its annotated average
fig, axes = plt.subplots(1, 1, figsize=(8, 5))

axes.plot(df_extended["time"][df_extended["time"].dt.year.between(1990, 2000)], df_extended[df_extended["time"].dt.year.between(1990, 2000)]['lowess_gc'], color='orange', label='LOWESS SLR (gc)')
axes.plot(df_extended["time"][df_extended["time"].dt.year.between(1990, 2000)], df_extended[df_extended["time"].dt.year.between(1990, 2000)]['lowess_wl'], color='blue', label='LOWESS SLR (wl)')

# Annotate the line
axes.annotate(f'{bathy_mean_SLR_ref_wl:.2f}', xy=(mean_time, bathy_mean_SLR_ref_wl), xytext=(mean_time, bathy_mean_SLR_ref_wl - 5),
            arrowprops=dict(arrowstyle='->'), fontsize=8)
axes.annotate(f'{bathy_mean_SLR_ref_gc:.2f}', xy=(mean_time, bathy_mean_SLR_ref_gc), xytext=(mean_time, bathy_mean_SLR_ref_gc - 5),
            arrowprops=dict(arrowstyle='->'), fontsize=8)

# Customize the first plot
axes.set_title('Mean SLR 1990 - 2000 within DFM domain and its annotate average')
axes.set_xlabel('Time')
axes.set_ylabel('Mean Water Level Difference (mm)')
axes.grid()
axes.legend()


#%%
SLR_start_event_wl = df_extended['lowess_wl'][df_extended['time'] == pd.to_datetime(start_date_event)].values[0]
SLR_start_event_wl_minstd = SLR_start_event_wl - std_SLR_DFM_wl.mean()
SLR_start_event_wl_maxstd = SLR_start_event_wl + std_SLR_DFM_wl.mean()

SLR_start_event_gc = df_extended['lowess_gc'][df_extended['time'] == pd.to_datetime(start_date_event)].values[0]
SLR_start_event_gc_minstd = SLR_start_event_gc - std_SLR_DFM_gc.mean()
SLR_start_event_gc_maxstd = SLR_start_event_gc + std_SLR_DFM_gc.mean()

print("ISIMIP geocentric water level ds, excl VLM")
print(f"SLR at start of event (2019-03-09): {SLR_start_event_gc:.2f} mm")
print(f"SLR at start of event -std: {SLR_start_event_gc_minstd:.2f} mm")
print(f"SLR at start of event +std: {SLR_start_event_gc_maxstd:.2f} mm")

print("\nISIMIP water level ds, incl VLM")
print(f"SLR at start of event (2019-03-09): {SLR_start_event_wl:.2f} mm")
print(f"SLR at start of event -std: {SLR_start_event_wl_minstd:.2f} mm")
print(f"SLR at start of event +std: {SLR_start_event_wl_maxstd:.2f} mm")





# %%
#################################################################################
###################### Validate with PSMSL RLR data ##############################  
#################################################################################
# Function to calculate linear and LOWESS trends and plot
def lin_lowess_trend_plot(data, date_col=None, waterlevel_col='sla',
                          remove_mean=True, figure_plotting=True, lowess_frac=1):
    # -----------------------------
    # Prepare series
    # -----------------------------
    if isinstance(data, pd.DataFrame):
        if date_col is not None:
            data = data.set_index(date_col)
        y_series = data[waterlevel_col]
    else:
        y_series = data

    # x in days since start
    x_days = (y_series.index - y_series.index[0]).total_seconds() / (3600*24)
    y = y_series.values
    not_nan = ~np.isnan(y)

    # -----------------------------
    # Linear trend
    # -----------------------------
    slope_lin, intercept_lin, _, _, _ = linregress(x_days[not_nan], y[not_nan])
    linear_fit = intercept_lin + slope_lin * x_days
    slope_lin_per_year = slope_lin * 365.25

    # -----------------------------
    # LOWESS trend
    # -----------------------------
    # Resample daily to reduce noise if needed
    y_daily = y_series.resample('D').mean().interpolate()
    x_daily = (y_daily.index - y_series.index[0]).days

    lowess_fit = sm.nonparametric.lowess(y_daily.values, x_daily, frac=lowess_frac, return_sorted=False)
    # Interpolate back to original dates
    interp_lowess = interp1d(x_daily, lowess_fit, kind='linear', fill_value='extrapolate')
    lowess_original = interp_lowess(x_days)

    # Slope of LOWESS trend (per year)
    slope_lowess, _, _, _, _ = linregress(x_days, lowess_original)
    slope_lowess_per_year = slope_lowess * 365.25

    # -----------------------------
    # Plotting
    # -----------------------------
    if figure_plotting:
        plt.figure(figsize=(12,5))
        plt.plot(x_days, y, label='Original data', alpha=0.5)
        plt.plot(y_series.index, linear_fit, '--', color='green', label=f'Linear Trend ({slope_lin_per_year:.2f}/yr)')
        plt.plot(y_series.index, lowess_original, '-', color='red', label=f'LOWESS Trend ({slope_lowess_per_year:.2f}/yr)')
        plt.xlabel('Time (monthly)')
        plt.ylabel('Sea level (mm)')
        plt.title('Monthly PSMSL tide gauge data for Durban (SA)')
        plt.legend()
        plt.show()

    # -----------------------------
    # Return results
    # -----------------------------
    return {
        'slope_linear_per_year': slope_lin_per_year,
        'slope_lowess_per_year': slope_lowess_per_year,
        'linear_fit': linear_fit,
        'lowess_fit': lowess_original,
        'y': y,
        'x_days': x_days,
        'time_index': y_series.index
    }

# Convert decimal year to datetime
def decimal_year_to_datetime(dec_year):
    year = np.floor(dec_year).astype(int)
    rem = dec_year - year
    dt = pd.to_datetime(year, format="%Y") + pd.to_timedelta(rem*365.25, unit='D')
    return dt

#%%
# Load PSMSL RLR data for Durban (station ID 284)
url = "https://psmsl.org/data/obtaining/rlr.monthly.data/284.rlrdata"
r = requests.get(url)
text = r.text

# Read the data into a DataFrame
df_psmsl_durban = pd.read_csv(StringIO(text), sep=";", header=None, 
                              names=["dec_year", "sla", "flag", "other"], usecols=[0,1,2])

# Clean the data
df_psmsl_durban = df_psmsl_durban[df_psmsl_durban['flag'] == 0]
df_psmsl_durban['sla'] = pd.to_numeric(df_psmsl_durban['sla'], errors='coerce')
df_psmsl_durban['sla'].replace(-99999, np.nan, inplace=True)
df_psmsl_durban['time'] = decimal_year_to_datetime(df_psmsl_durban['dec_year'])

# Run the trend plotting function
results = lin_lowess_trend_plot(df_psmsl_durban, date_col='time', waterlevel_col='sla')
print("Linear slope (mm/yr):", results['slope_linear_per_year'])
print("LOWESS slope (mm/yr):", results['slope_lowess_per_year'])


# %%

slope_gc_per_year = slope_gc * 365.25
slope_wl_per_year = slope_wl * 365.25

print(f"Geocentric trend: {slope_gc_per_year:.3f} units/year")
print(f"Water level trend: {slope_wl_per_year:.3f} units/year")



slope_lowess_gc, intercept_lowess_gc, r_val, p_val, std_err = linregress(x_extended, df_extended['lowess_gc'])
slope_lowess_gc_per_year = slope_lowess_gc * 365.25  # convert from per day to per year

slope_lowess_wl, intercept_lowess_wl, r_val, p_val, std_err = linregress(x_extended, df_extended['lowess_wl'])
slope_lowess_wl_per_year = slope_lowess_wl * 365.25

print(f"LOWESS trend geocentric: {slope_lowess_gc_per_year:.3f} units/year")
print(f"LOWESS trend waterlevel: {slope_lowess_wl_per_year:.3f} units/year")



# %%

# %% Use pixi evironment: compass-wflow
print("Loading packages...")
import pandas as pd
import numpy as np
from datetime import datetime as datetime
import geopandas as gpd
import xarray as xr
import matplotlib.pyplot as plt
from shapely.geometry import Point
import contextily as ctx
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from hydromt_wflow import WflowModel
import platform
from pathlib import Path
from geopy.distance import geodesic
from hydromt.stats import skills as skillstats  # NSE, KGE, etc.
import calendar

#%%
if platform.system() == "Windows":
    P = Path("p:/")
else:
    # adjust this if your cluster mounts “P:” somewhere else, e.g. “/mnt/p”
    P = Path("/p/")

# Now build your base directories
BASE        = P / "11210471-001-compass"
BASE_RUNS   = P / "11210471-001-compass" / "03_Runs"  / "sofala" / "Idai" / "wflow"
BASE_SFINCS = P / "11210471-001-compass" / "03_Runs"  / "sofala" / "Idai" / "sfincs" 
BASE_MODELS = P / "11210471-001-compass" / "02_Models" / "sofala" / "Idai" / "wflow"
BASE_DATA   = P / "11210471-001-compass" / "01_Data"

# Load base model for shapefiles
mod_ini = WflowModel(root=str(BASE_MODELS), mode="r+", config_fn=str(BASE_MODELS / "wflow_sbm.toml"))


# %% ##############################################################
############################## GRDC DATA ##########################
## ################################################################
print("Loading GRDC dataset...")
# Retrieved from https://grdc.bafg.de/data/data_portal/

grdc_path = BASE_DATA / "GRDC" / "GRDC-Daily.nc"
GRDC_data = xr.open_dataset(str(grdc_path))

print("Initializing 30-year WFLOW model...")
def make_model_warmup(subfolder: str):
    root = BASE_RUNS / subfolder / "warmup"
    cfg = root / "wflow_sbm.toml"
    return WflowModel(root=str(root), mode="r", config_fn=str(cfg))

wflow_30yr_hist = make_model_warmup("event_precip_era5_daily_CF0_30yr_1954")

# %% ##############################################################
######################### SELECT GRDC STATIONS ####################
###################################################################
print("Selecting GRDC stations...")
wanted_rivers = ['BUDZI', 'RIOPUNGOE']
names = GRDC_data['river_name'].values.astype(str)
idx = [i for i, n in enumerate(names) if n in wanted_rivers]
selected = GRDC_data.isel(id=idx)

print(f"Found {len(selected.id)} matching GRDC stations.")

lats, lons = selected['geo_y'].values, selected['geo_x'].values
gdf_stations = gpd.GeoDataFrame({'name': selected['river_name'].values},
                                 geometry=gpd.points_from_xy(lons, lats),
                                 crs="EPSG:4326")

# %% ##############################################################
######################## LOAD REGIONS & MATCH #####################
###################################################################
print("Loading WFLOW and SFINCS regions...")
gdf_wflow = gpd.read_file(BASE_MODELS / "staticgeoms" / "basins.geojson").to_crs("EPSG:4326")
gdf_sfincs = gpd.read_file(BASE / "02_Models" / "sofala" / "Idai" / "sfincs" / "gis" / "region.geojson").to_crs("EPSG:4326")

print("Finding closest GRDC stations to SFINCS domain...")
sfincs_geom = gdf_sfincs.unary_union
gdf_stations['distance'] = gdf_stations.geometry.distance(sfincs_geom)
closest_stations = gdf_stations.nsmallest(2, 'distance')
print(f"Closest stations: {list(closest_stations['name'])}")


# %% ##############################################################
############# MATCH WFLOW GRID POINTS TO CLOSEST STATIONS #########
###################################################################
print("Matching GRDC stations to WFLOW river grid points...")
q_river = wflow_30yr_hist.results['output']['q_river']
lat_vals, lon_vals = q_river['lat'].values, q_river['lon'].values

station_coords = [(pt.y, pt.x) for pt in closest_stations.geometry]
station_series, matched_indices, matched_gridpoints = [], [], []

for lat, lon in station_coords:
    lat_idx, lon_idx = (abs(lat_vals - lat)).argmin(), (abs(lon_vals - lon)).argmin()
    matched_lat, matched_lon = float(lat_vals[lat_idx]), float(lon_vals[lon_idx])
    matched_gridpoints.append((matched_lat, matched_lon))
    matched_indices.append((lat_idx, lon_idx))
    station_series.append(q_river[:, lat_idx, lon_idx])

print("Matched WFLOW grid points to GRDC stations.")

matches_df = pd.DataFrame({
    'station_index': closest_stations.index,
    'GRDC_station_lat': [lat for lat, _ in station_coords],
    'GRDC_station_lon': [lon for _, lon in station_coords],
    'wflow_matched_grid_lat': [lat for lat, _ in matched_gridpoints],
    'wflow_matched_grid_lon': [lon for _, lon in matched_gridpoints],
    'lat_idx': [i for i, _ in matched_indices],
    'lon_idx': [j for _, j in matched_indices],
    'distance_km': [
        geodesic((a, b), (c, d)).km for (a, b), (c, d) in zip(station_coords, matched_gridpoints)
    ]
})

print(matches_df)


# %% ##############################################################
################### PLOT SPATIAL MATCHING OVERVIEW ################
###################################################################

print("Plotting spatial match of GRDC and WFLOW points...")
matched_wflow_gdf = gpd.GeoDataFrame(
    geometry=[Point(lon, lat) for lat, lon in matches_df[['wflow_matched_grid_lat', 'wflow_matched_grid_lon']].values],
    crs='EPSG:4326'
)

fig, ax = plt.subplots(figsize=(12,7))
gdf_wflow.plot(ax=ax, edgecolor='skyblue', facecolor='skyblue', alpha=0.8, label="Wflow basins")
gdf_sfincs.plot(ax=ax, edgecolor='pink', facecolor='pink', alpha=0.5, label="SFINCS region")
mod_ini.geoms["rivers"].plot(ax=ax, color='white', zorder=1)
gdf_stations.plot(ax=ax, color='red', label='All GRDC stations', zorder=10)
closest_stations.plot(ax=ax, color='orange', label='Closest GRDC stations (G)', zorder=10)
matched_wflow_gdf.plot(ax=ax, marker='x', color="blue", linewidth=1.5, markersize=50, label='Matched wflow grid goint (Q)', zorder=10)

# Add text labels to matched WFLOW grid points
for i, point in enumerate(matched_wflow_gdf.geometry):
    ax.text(point.x + 0.02, point.y + 0.05, f"Q{i+1}", fontsize=12, fontweight='bold', color='blue', ha='left', zorder=10, 
            bbox=dict(boxstyle="round,pad=0.1", facecolor='grey', alpha=0.5, edgecolor='grey'))

for i, point in enumerate(closest_stations.geometry):
    ax.text(point.x + 0.02, point.y - 0.1, f"G{i+1}", fontsize=12, fontweight='bold', color='orange', ha='left', zorder=10, 
            bbox=dict(boxstyle="round,pad=0.1", facecolor='grey', alpha=0.5, edgecolor='grey'))

gdf_gauges = mod_ini.geoms["gauges_locs"]

ax.set_xlim(gdf_stations.geometry.x.min() - 0.2, gdf_stations.geometry.x.max() + 1)
ax.set_ylim(gdf_stations.geometry.y.min() - 0.5, gdf_stations.geometry.y.max() + 0.5)

# Add basemap (LOWER zoom = faster)
ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=9, crs=gdf_stations.crs, attribution=False, zorder=0)

txt = ax.text(
    33.9, -20.75,  # x, y in figure coordinates (0=left/bottom, 1=right/top)
    "Tiles © Esri -- Source: Esri, i-cubed, USDA, USGS, AEX, \nGeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and \nthe GIS User Community",
    fontsize=7,
    color='white',
    alpha=0.7,
    ha='left',
    va='bottom',
    zorder=20,
)

x = gdf_stations.geometry.x
y = gdf_stations.geometry.y
xticks = np.arange(32.5, 36, 0.5)    
yticks = np.arange(-21, -17.5, 0.5) 
ax.set_xticklabels([f"{x}°E" for x in xticks], fontsize=12)
ax.set_yticklabels([f"{abs(y)}°S" for y in yticks], fontsize=12)

ax.legend(loc='upper right')
plt.tight_layout()
plt.savefig('../figures/fS4.png', dpi=300, bbox_inches="tight")
plt.savefig('../figures/fS4.pdf', dpi=300, bbox_inches="tight")
plt.show()


# %% ##############################################################
############## COMPARE WFLOW VS. GRDC STATISTICALLY ###############
###################################################################
print("Creating monthly climatology comparison between WFLOW and GRDC with uncertainty bands and statistics...")
river_names = ['Buzi', 'Pungwe'] # English names for titles and annotations

print("Precomputing monthly climatologies...")

months = np.arange(1, 13)
month_labels = [calendar.month_abbr[m] for m in months]

# Storage
climatology = {}

# --- Precompute climatologies for all rivers ---
for i, river in enumerate(wanted_rivers):
    grdc_ts = GRDC_data['runoff_mean'].isel(id=idx[matches_df['station_index'].iloc[i]])
    wflow_ts = station_series[i]

    climatology[river] = {
        "grdc_mean": grdc_ts.groupby("time.month").mean().values,
        "grdc_std": grdc_ts.groupby("time.month").std().values,
        "wflow_mean": wflow_ts.groupby("time.month").mean().values,
        "wflow_std": wflow_ts.groupby("time.month").std().values,
    }

#%%
print("Plotting...")

# --- Plotting ---
fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True, sharey=True)
colors = ["#0C75C0", "#4DA54E"]

for i, (river, river_name) in enumerate(zip(wanted_rivers, river_names)):
    clim = climatology[river]

    # Compute KGE once per river
    kge = skillstats.kge(
        xr.DataArray(clim["wflow_mean"], dims=["month"], coords={"month": months}),
        xr.DataArray(clim["grdc_mean"], dims=["month"], coords={"month": months}),
        dim="month"
    )["kge"].values.round(2)

    # Plot GRDC
    axs[i].plot(months, clim["grdc_mean"], label=f"GRDC (G{i+1})", color="orange")
    axs[i].fill_between(months,
                        np.clip(clim["grdc_mean"] - clim["grdc_std"], 0, None),
                        clim["grdc_mean"] + clim["grdc_std"],
                        color="orange", alpha=0.3)

    # Plot WFLOW
    axs[i].plot(months, clim["wflow_mean"], label=f"wflow (Q{i+1}), KGE={kge}", color=colors[i])
    axs[i].fill_between(months,
                        np.clip(clim["wflow_mean"] - clim["wflow_std"], 0, None),
                        clim["wflow_mean"] + clim["wflow_std"],
                        color=colors[i], alpha=0.3)

    # Style
    axs[i].set_ylabel("Discharge [m³/s]", fontsize=12)
    axs[i].set_title(f"{river_name}: Monthly Climatology", fontsize=13)
    axs[i].set_xticks(months)
    axs[i].set_xticklabels(month_labels)
    axs[i].set_xlim(1, 12)
    axs[i].legend(loc="upper right", fontsize=10)
    axs[i].grid(True, linestyle="--", alpha=0.6)
    axs[i].tick_params(axis="both", labelsize=11)

axs[0].text(0.0, 1.03, "(a)", transform=axs[0].transAxes,
             fontsize=14, fontweight='bold', va='top', ha='left')

# Subplot (b) - second plot
axs[1].text(0.0, 1.03, "(b)", transform=axs[1].transAxes,
             fontsize=14, fontweight='bold', va='top', ha='left')

axs[-1].set_xlabel("Month", fontsize=12)
fig.suptitle("Observed (GRDC) vs. simulated (wflow) monthly discharge (1954–1984)", fontsize=14, y=0.95)

plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig("../figures/fS5.png", dpi=300, bbox_inches="tight")
plt.savefig("../figures/fS5.pdf", dpi=300, bbox_inches="tight")
plt.show()


# %%
print("Comparing wflow time series to GRDC records with skill statistics...")

fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
color = ["#0C75C0", "#4DA54E"]  # Blue for BUDZI, Green for RIOPUNGOE

for i, (river, river_name) in enumerate(zip(wanted_rivers, river_names)):
    print(f"Processing {river}...")

    idx_grdc = [i for i, n in enumerate(names) if n == river][0]
    grdc_ts = GRDC_data['runoff_mean'].isel(id=idx[matches_df['station_index'].iloc[i]])
    wflow_ts = station_series[i]

    # Stats
    print(grdc_ts)
    print(wflow_ts)
    common_time = np.intersect1d(wflow_ts.time.values, grdc_ts.time.values)
    wflow_aligned = wflow_ts.sel(time=common_time)
    grdc_aligned = grdc_ts.sel(time=common_time)

    wflow_stat = grdc_aligned.chunk(dict(time=-1))
    grdc_stat = wflow_aligned.chunk(dict(time=-1))

    kge = skillstats.kge(wflow_stat, grdc_stat, dim='time')["kge"].values.round(2)

    # Plot time series
    axs[i].plot(wflow_ts.time, wflow_ts.values, label=f'wflow (Q{i+1}),  KGE: {kge}', color=color[i])
    axs[i].plot(grdc_ts.time, grdc_ts.values, label=f'GRDC (G{i+1})', color='orange')

    # Add max lines
    axs[i].axhline(grdc_ts.max().values, linestyle='--', color='orange', label='GRDC max')
    axs[i].axhline(wflow_ts.values.max(), linestyle='--', color=color[i], label='wflow max')

    axs[i].set_ylabel('Discharge [m³/s]', fontsize=12)
    axs[i].set_title(f"The {river_name} River; Gauge Q{i+1}", fontsize=13)
    axs[i].legend(loc='upper right', fontsize=10)
    axs[i].grid(True)
    axs[i].tick_params(axis='both', which='major', labelsize=11)
    axs[i].set_xlim([np.datetime64('1954-01-01'), np.datetime64('1984-12-31')])

axs[0].text(0, 1.03, "(a)", transform=axs[0].transAxes, 
            fontsize=14, fontweight='bold', va='top', ha='left')
axs[1].text(0, 1.03, "(b)", transform=axs[1].transAxes, 
            fontsize=14, fontweight='bold', va='top', ha='left')

fig.suptitle(f"Observed (GRDC) vs. simulated (wflow) discharge (1954–1984)", fontsize=14) 
axs[1].set_xlabel('Time', fontsize=12)
fig.savefig('../figures/fS6.png', dpi=300)
fig.savefig('../figures/fS6.pdf', dpi=300)
plt.show()

print("Daily time series comparison with statistics saved.")


# %%

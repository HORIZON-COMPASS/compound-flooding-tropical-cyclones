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
from hydromt_wflow import WflowModel
from pathlib import Path
from geopy.distance import geodesic
from hydromt.stats import skills as skillstats  # NSE, KGE, etc.
import calendar
from hydromt import data_catalog
import string  
import matplotlib.dates as mdates
import matplotlib.patheffects as pe
from matplotlib.ticker import FuncFormatter
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

#%%
# Set directories
BASE        = Path("p:/11210471-001-compass/")
BASE_RUNS   = BASE / "03_Runs/sofala/Idai/wflow"
BASE_SFINCS = BASE / "02_Models/sofala/Idai/sfincs" 
BASE_MODELS = Path("p:/11210471-001-compass/02_Models/sofala/Idai/wflow")
BASE_DATA   = Path("p:/11210471-001-compass/01_Data/GRDC")

# Load base model for shapefiles
mod_ini = WflowModel(root=str(BASE_MODELS), mode="r+", config_fn=str(BASE_MODELS / "wflow_sbm.toml"))

# Load wflow and SFINCS regions and basins
gdf_wflow = gpd.read_file(BASE_MODELS / "staticgeoms" / "basins.geojson").to_crs("EPSG:4326")
gdf_sfincs = gpd.read_file(BASE_SFINCS / "gis" / "region.geojson", driver="GeoJSON").to_crs("EPSG:4326")

# Make sure to change "C:/Code/COMPASS/" to directory where to github repo is stored
print("Initializing 30-year WFLOW model...")
def make_model_warmup(subfolder: str):
    root = (Path(BASE_RUNS) / subfolder / "warmup").resolve()
    cfg = (Path(BASE_RUNS) / subfolder / "warmup" / "wflow_sbm.toml").resolve()
    return WflowModel(root=str(root), mode="r", config_fn=str(cfg))

wflow_30yr_hist = make_model_warmup("event_precip_era5_daily_CF0_30yr_1954")
wflow_30yr = make_model_warmup("event_precip_era5_hourly_zarr_CF0_30yr")

#%%
# Function to plot gauges, basins, rivers, and SFINCS region
def plot_gauges_with_basins_and_sfincs(
    rivers, gauges, basins_gdf, sfincs_region,
    offsets=(None), figsize=(8, 8)
):
    if offsets is None:
        offsets = {}

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Plot wflow basins
    gdf_wflow.plot(ax=ax, edgecolor='skyblue', facecolor='skyblue', alpha=0.8, label="Wflow basins")
    basins_gdf.plot(ax=ax, facecolor="none", edgecolor="#555555", linewidth=1, zorder=1, label="Wflow basins")

    # Plot rivers
    rivers.plot(ax=ax, linewidth=0.8, color="white", zorder=1, label="Rivers")

    # Plot SFINCS region outline
    gdf_sfincs.plot(ax=ax, edgecolor='pink', facecolor='pink', alpha=0.5, label="SFINCS region")

    # Plot gauges
    gauges.plot(ax=ax, color="red", markersize=15, zorder=3, label="Gauges")

    # Add gauge numbers with offsets
    for _, r in gauges.iterrows():
        g = int(r["index"])
        dx, dy = offsets.get(g, (0.0, -0.07))
        ax.text(
            r.geometry.x + dx,
            r.geometry.y + dy,
            str(g),
            fontsize=8,
            weight="bold",
            ha="center",
            va="center",
            zorder=4,
            path_effects=[
                pe.Stroke(linewidth=2.5, foreground="white"),
                pe.Normal(),
            ],
        )

    # Add basemap (LOWER zoom = faster)
    ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=9, crs=gdf_wflow.crs, attribution=False, zorder=0)

    txt = ax.text(
        33.9, -21.07,  # x, y in figure coordinates (0=left/bottom, 1=right/top)
        "Tiles © Esri -- Source: Esri, i-cubed, USDA, USGS, AEX, \nGeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and \nthe GIS User Community",
        fontsize=6,
        color='white',
        alpha=0.7,
        ha='left',
        va='bottom',
        zorder=20,
    )
    # Axis labels & formatting
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.2f}°E"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.2f}°N"))

    ax.set_title("Gauges within wflow and SFINCS model domains")

    legend_handles = [
        mlines.Line2D([], [], color="white", markeredgecolor='black', linewidth=1.5, label="Rivers"),
        mpatches.Patch(facecolor="none", edgecolor="#555555", linewidth=1, label="wflow basins"),
        mpatches.Patch(facecolor="skyblue", label="wflow domain"),
        mpatches.Patch(facecolor="pink", label="SFINCS region"),
        mlines.Line2D([], [], color="red", marker="o", linestyle="None", markersize=6, label="Gauges"),
    ]
    ax.legend(handles=legend_handles, loc="upper right")

    fig.savefig("../figures/fS4.png", dpi=300, bbox_inches="tight")

    return fig, ax


fig, ax = plot_gauges_with_basins_and_sfincs(
    rivers=mod_ini.geoms["rivers"],
    gauges=mod_ini.geoms["gauges_locs"],
    basins_gdf=gdf_wflow,
    sfincs_region=gdf_sfincs,
    offsets={1: (-0.05, 0), 10: (0.09, -0.01)}
)


# %% ##############################################################
############################## GRDC DATA ##########################
## ################################################################
print("Loading GRDC dataset...")
# Retrieved from https://grdc.bafg.de/data/data_portal/

# Make sure to retrieve the data and store it in the correct folder
grdc_path = BASE_DATA / "GRDC-Daily.nc"
GRDC_data = xr.open_dataset(str(grdc_path))

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
########################## MATCH REGIONS   ########################
###################################################################

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

axs[0].text(0.0, 1.07, "(a)", transform=axs[0].transAxes,
             fontsize=14, fontweight='bold', va='top', ha='left')

# Subplot (b) - second plot
axs[1].text(0.0, 1.07, "(b)", transform=axs[1].transAxes,
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

axs[0].text(0, 1.07, "(a)", transform=axs[0].transAxes, 
            fontsize=14, fontweight='bold', va='top', ha='left')
axs[1].text(0, 1.07, "(b)", transform=axs[1].transAxes, 
            fontsize=14, fontweight='bold', va='top', ha='left')

fig.suptitle(f"Observed (GRDC) vs. simulated (wflow) discharge (1954–1984)", fontsize=14) 
axs[1].set_xlabel('Time', fontsize=12)
fig.savefig('../figures/fS6.png', dpi=300)
fig.savefig('../figures/fS6.pdf', dpi=300)
plt.show()

print("Daily time series comparison with statistics saved.")








# %% ##############################################################
################# COMPARE WFLOW VS. GLOFAS DATA ###################
###################################################################
import rioxarray 

# Load data catalog
data_cats = ['../../Workflows/03_data_catalogs/datacatalog_general.yml']    
data_cat = data_catalog.DataCatalog(data_cats)

# Time range
start_date = "20190307 000000"
end_date = "20190325 000000"
start_dt = datetime.strptime(start_date, "%Y%m%d %H%M%S")
end_dt = datetime.strptime(end_date, "%Y%m%d %H%M%S")

# Load GLOFAS dataset (version 3.1 and 4.0 with 0.1 and 0.05 degree spatial resolution, respectively)
# Can be downloaded here: https://ewds.climate.copernicus.eu/datasets/cems-glofas-historical
glofas_ds_v31 = data_cat.get_rasterdataset('glofas_era5_v31', geom=gdf_wflow)
# glofas_ds_v40 = data_cat.get_rasterdataset('glofas_era5_v40', geom=gdf_wflow)
# glofas_ds_v40_event = data_cat.get_rasterdataset('glofas_era5_v40_wflowMZ_Idai', geom=gdf_wflow)
# Rename valid_time to time
# glofas_ds_v40 = glofas_ds_v40.rename({"valid_time": "time"})

glofas_uparea_v40 = xr.open_dataset("p:/11210471-001-compass/01_Data/glofas_v4/uparea_glofas_v4_0.nc")
glofas_uparea_v31 = rioxarray.open_rasterio("p:/wflow_global/hydromt/hydro/glofas_era5/glofas_v31_uparea.tif")
glofas_uparea_v40 = glofas_uparea_v40.rio.write_crs("EPSG:4326", inplace=True)
glofas_uparea_v40_clipped = glofas_uparea_v40.rio.clip(gdf_wflow.geometry, gdf_wflow.crs)

glofas_uparea_v31 = glofas_uparea_v31.rio.write_crs("EPSG:4326")  # if not set
glofas_uparea_v31_clipped = glofas_uparea_v31.rio.clip(gdf_wflow.geometry, gdf_wflow.crs)

# Interpolate upstream area to GloFAS v4.0 grid
glofas_uparea_interp = glofas_uparea_v40_clipped['uparea'].interp(
    latitude=glofas_ds_v31['latitude'],
    longitude=glofas_ds_v31['longitude'],
    method="nearest" 
)
glofas_ds_v31 = glofas_ds_v31.assign_coords(uparea=glofas_uparea_interp)

# %% ##############################################################
##################### SELECT GLOFAS GRID CELLS ####################
###################################################################
# Convert GloFAS uparea to GeoDataFrame
# df = (glofas_ds_v31.to_dataframe().dropna().reset_index())

# # Explicitly store original WGS84 coordinates
# df["lat_wgs84"] = df["latitude"]
# df["lon_wgs84"] = df["longitude"]

# gdf_glofas_uparea_v31 = gpd.GeoDataFrame(df, geometry=[Point(lon, lat) for lat, lon in zip(df['latitude'], df['longitude'])],
#                        crs="EPSG:4326")

# # Reproject to metric CRS for distance calculation
# gdf_glofas_uparea_v31 = gdf_glofas_uparea_v31.to_crs("EPSG:32736")
# gdf_sfincs_wsg84 = gdf_sfincs.to_crs("EPSG:32736")

# # Exclude cells inside SFINCS
# gdf_outside = gdf_glofas_uparea_v31[~gdf_glofas_uparea_v31.within(gdf_sfincs_wsg84.unary_union)]
# # Compute distance to SFINCS boundary
# sfincs_boundary = gdf_sfincs_wsg84.unary_union.boundary
# gdf_outside["dist_to_sfincs_m"] = gdf_outside.distance(sfincs_boundary)

# # Keep only cells within 2 km
# gdf_near = gdf_outside[gdf_outside["dist_to_sfincs_m"] <= 2000]

# # Select top 3 by upstream area
# gdf_top3 = gdf_near.sort_values("uparea", ascending=False).head(3)  
# gdf_top3 = gdf_near.sort_values("uparea", ascending=False).head(3)

# # Export for QGIS (back to lat/lon)
# gdf_top3.to_crs("EPSG:4326").to_file("C:/Code/Paper_1/QGIS/top3_uparea_near_sfincs_v31.geojson", driver="GeoJSON")

# print("Top 3 upstream GloFAS cells within 2 km of SFINCS exported.")

#%%
# ------------------------------
# 1. Sum discharge over the event
# ------------------------------
da_event = glofas_ds_v40.sel(time=slice(start_dt, end_dt))  # select event period
event_total = da_event.sum(dim="time")  # shape: (lat, lon)

# ------------------------------
# 2. Flatten to 2D arrays for GeoDataFrame
# ------------------------------
lats = glofas_ds_v40.latitude.values
lons = glofas_ds_v40.longitude.values

lat_grid, lon_grid = np.meshgrid(lats, lons, indexing="ij")  # shape: (lat, lon)
uparea = glofas_ds_v40.uparea.values  # shape: (lat, lon)
event_discharge = event_total.values  # shape: (lat, lon)

# ------------------------------
# 3. Create GeoDataFrame
# ------------------------------
df = pd.DataFrame({
    "latitude": lat_grid.ravel(),
    "longitude": lon_grid.ravel(),
    "uparea": uparea.ravel(),
    "event_discharge": event_discharge.ravel()
})

# Drop NaNs if necessary
df = df.dropna(subset=["uparea", "event_discharge"])

# Make GeoDataFrame
gdf = gpd.GeoDataFrame(df, geometry=[Point(xy) for xy in zip(df["longitude"], df["latitude"])], crs="EPSG:4326")

# ------------------------------
# 4. Filter by SFINCS buffer
# ------------------------------
gdf_utm = gdf.to_crs("EPSG:32736")
gdf_sfincs_utm = gdf_sfincs.to_crs("EPSG:32736")

sfincs_boundary = gdf_sfincs_utm.unary_union.boundary
gdf_utm["dist_to_sfincs_m"] = gdf_utm.geometry.distance(sfincs_boundary)

# Keep only cells outside SFINCS and within 2 km
gdf_near = gdf_utm[(gdf_utm["dist_to_sfincs_m"] <= 2000) & (~gdf_utm.within(gdf_sfincs_utm.unary_union))]

# ------------------------------
# 5. Select top 3 by event discharge, break ties with upstream area
# ------------------------------
gdf_top3 = gdf_near.sort_values(["event_discharge", "uparea"], ascending=[False, False]).head(3)

# Export for QGIS if needed
gdf_top3.to_crs("EPSG:4326").to_file("c:/Code/Paper_1/QGIS/top3_glofas_high_discharge_v31.geojson", driver="GeoJSON")

print(gdf_top3[["latitude", "longitude", "uparea", "event_discharge", "dist_to_sfincs_m"]])




# %%
print("Selecting discharge from wflow and add coordinates...")
gauges = mod_ini.geoms["gauges_locs"]
gauges_wgs84 = gauges.to_crs(epsg=4326)
gauges_wgs84["lon"] = gauges_wgs84.geometry.x
gauges_wgs84["lat"] = gauges_wgs84.geometry.y

Q_wflow_30yr = wflow_30yr.results['netcdf']['Q']
df_wflow_event = pd.read_csv(BASE_RUNS / "event_precip_era5_hourly_zarr_CF0" / "events" / "run_default" / "wflow_dis_no_qbankfull.csv", parse_dates=['time'], index_col='time')
Q_wflow_event = xr.DataArray(df_wflow_event.values, dims=["time", "Q_gauges_locs"], 
                             coords={"time": df_wflow_event.index, "Q_gauges_locs": df_wflow_event.columns.astype(str)}, 
                             name="Q").to_dataset()
Q_wflow_event = Q_wflow_event['Q']


# Based on the resolution of GloFAS data, we only select our gauges with a upstream area > 500 km2 - similar to the approach for in Harrigan et al. (2020)
print("Selecting only gauges with upstream area >= 500 km2...")
gauges_filtered = mod_ini.geoms['gauges_locs'][mod_ini.geoms['gauges_locs'].uparea >= 500]

Q_wflow_30yr_filtered = Q_wflow_30yr.sel(Q_gauges_locs=gauges_filtered['index'].astype(str).values)
Q_wflow_event_filtered = Q_wflow_event.sel(Q_gauges_locs=gauges_filtered['index'].astype(str).values)

# Add correct lat and lon to gauges in wflow output
gauges_lookup = (gauges_wgs84.set_index("index"))
gauges_lookup.index = gauges_lookup.index.astype(str)
gauge_ids = Q_wflow_30yr_filtered.coords["Q_gauges_locs"].values.astype(str)
lon = xr.DataArray(gauges_lookup.loc[gauge_ids, "lon"].values, dims="Q_gauges_locs", 
                   coords={"Q_gauges_locs": gauge_ids}, name="lon")
lat = xr.DataArray(gauges_lookup.loc[gauge_ids, "lat"].values, dims="Q_gauges_locs",
                   coords={"Q_gauges_locs": gauge_ids}, name="lat")

Q_wflow_30yr_filtered = Q_wflow_30yr_filtered.assign_coords(lon=lon, lat=lat)
Q_wflow_event_filtered = Q_wflow_event_filtered.assign_coords(lon=lon, lat=lat)


print("Merge gauges when located within 10 km of each others...")
# Merge gauges when located within 10 km of each other
gdf_wflow_gauges = gpd.GeoDataFrame({"Q_gauges_locs": Q_wflow_30yr_filtered.Q_gauges_locs.values},
                       geometry=[Point(xy) for xy in zip(Q_wflow_30yr_filtered.lon.values,
                                                         Q_wflow_30yr_filtered.lat.values)], crs="EPSG:4326")
gdf_wflow_gauges_utm = gdf_wflow_gauges.to_crs("EPSG:32736")

cluster_id = np.full(len(gdf_wflow_gauges_utm), -1, dtype=int)
current_cluster = 0

# Simple clustering based on distance threshold
for i in range(len(gdf_wflow_gauges_utm)):
    if cluster_id[i] != -1:
        continue

    distances = gdf_wflow_gauges_utm.geometry.distance(gdf_wflow_gauges_utm.geometry.iloc[i])
    members = np.where(distances <= 10000)[0]

    cluster_id[members] = current_cluster
    current_cluster += 1

gdf_wflow_gauges["cluster"] = cluster_id

rep_station = (
    gdf_wflow_gauges.groupby("cluster")["Q_gauges_locs"]
    .min()
    .rename("rep_station")
)

gdf_wflow_gauges = gdf_wflow_gauges.merge(rep_station, left_on="cluster", right_index=True)

# Combine time series by cluster
def combine_by_cluster(ds):
    df = ds.to_dataframe().reset_index()

    df = df.merge(gdf_wflow_gauges[["Q_gauges_locs", "cluster", "rep_station"]], on="Q_gauges_locs")

    df_comb = (df.groupby(["time", "cluster"]).sum(numeric_only=True).reset_index())

    df_comb = df_comb.merge(rep_station, on="cluster")

    return (df_comb.set_index(["time", "rep_station"]).to_xarray().rename({"rep_station": "Q_gauges_locs"}))

Q_wflow_30yr_filt_combined = combine_by_cluster(Q_wflow_30yr_filtered)
Q_wflow_event_filt_combined = combine_by_cluster(Q_wflow_event_filtered)

gdf_rep = (gdf_wflow_gauges[gdf_wflow_gauges.Q_gauges_locs == gdf_wflow_gauges.rep_station].drop_duplicates("cluster"))

# gdf_rep.to_file("C:/Code/Paper_1/QGIS/wflow_filt_combined_gauges.geojson", driver="GeoJSON")

# Only select the top three combined gauges
Q_wflow_30yr_filt_combined = Q_wflow_30yr_filt_combined.sel(Q_gauges_locs=gdf_rep.head(3).rep_station.values)
Q_wflow_event_filt_combined = Q_wflow_event_filt_combined.sel(Q_gauges_locs=gdf_rep.head(3).rep_station.values)

#%%
# Select the top three GloFAS gauges/grid cells near the sfincs domain with the largest upstream area
glofas_v40_sel_list = []

for lat, lon in zip(gdf_top3["latitude"], gdf_top3["longitude"]):
    sel_point = glofas_ds_v40.sel(latitude=lat, longitude=lon, method="nearest")
    glofas_v40_sel_list.append(sel_point)

# Concatenate along a new 'points' dimension
glofas_v40_sel = xr.concat(glofas_v40_sel_list, dim="points")

# Export points for QGIS for sanity check
geometry = [Point(lon, lat) for lat, lon in zip(gdf_top3["latitude"], gdf_top3["longitude"])]
gdf_top3_export = gpd.GeoDataFrame(
    gdf_top3[['uparea', 'dist_to_sfincs_m']],
    geometry=geometry,
    crs="EPSG:4326"
)

glofas_v40_sel = glofas_v40_sel.assign_coords(
    Q_gauges_locs=("points", [str(i) for i in range(1, len(gdf_top3)+1)])  # '1', '2', '3'
)
glofas_v40_sel = glofas_v40_sel.swap_dims({"points": "Q_gauges_locs"})

# Export as GeoJSON
# gdf_top3_export.to_file("C:/Code/Paper_1/QGIS/top3_glofas_nearest_sfincs.geojson", driver="GeoJSON")

# Q_wflow_event_daily = (
#     Q_wflow_event_filt_combined
#     .resample(time="1D")
#     .mean()   # or .sum() depending on your definition
# )

# Align time series
Q_wflow_30yr_aligned, glofas_30yr_aligned = xr.align(Q_wflow_30yr_filt_combined, glofas_v40_sel, join="inner")
Q_wflow_event_aligned, glofas_event_aligned = xr.align(Q_wflow_event_filt_combined, glofas_v40_sel, join="inner")


#%%
# Plot all gauges time series for the event
print("Plotting time series comparison for all gauges...")

# Get all gauge IDs
gauge_ids = Q_wflow_event_aligned['Q_gauges_locs'].values

# Prepare DataFrame to store stats
df_stats = pd.DataFrame(columns=['Gauge #', 'wflow [m³]', 'GloFAS [m³]', 'wflow / GloFAS [%]'])

# Plotting setup
n = len(gauge_ids)
ncols = 2
nrows = int(np.ceil(n / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(12, nrows*3), sharex=True, sharey=False)
axes = axes.flatten()

letters = list(string.ascii_lowercase)

for i, g in enumerate(gauge_ids):
    ax = axes[i]
    print(g)
    # Select timeseries
    wflow_ts = Q_wflow_event_aligned.sel(Q_gauges_locs=g, time=slice(start_dt, end_dt))
    glofas_ts = glofas_event_aligned.sel(Q_gauges_locs=g, time=slice(start_dt, end_dt))
    
    # Plot
    ax.plot(wflow_ts.time, wflow_ts['Q'], label="WFLOW", lw=1.5)
    ax.plot(glofas_ts.time, glofas_ts, label="GloFAS", lw=1.5, alpha=0.7)
    
    # Time differences in seconds
    dt_wflow = np.diff(wflow_ts.time.values).astype("timedelta64[s]").astype(float)
    dt_wflow = np.append(dt_wflow, dt_wflow[-1])
    
    dt_glofas = np.diff(glofas_ts.time.values).astype("timedelta64[s]").astype(float)
    dt_glofas = np.append(dt_glofas, dt_glofas[-1])
    
    # Total volume
    total_Q = (wflow_ts['Q'].values * dt_wflow).sum()
    total_G = (glofas_ts.values * dt_glofas).sum()
    
    # Save stats
    df_stats.loc[i] = [g, total_Q, total_G, total_Q / total_G * 100]
    
    # Plot styling
    ax.set_title(f"Gauge {g}")
    ax.set_ylabel("Discharge [m³/s]")
    ax.grid(alpha=0.3)
    if i == 0:
        ax.legend()
    # Letter annotation
    ax.text(0.02, 1.08, f'({letters[i]})', transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left')
    
    # Formatting x-axis as days in 2019
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))  # every 2 days
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d"))
    ax.tick_params(axis='x')
    ax.set_xlabel(" Day in March 2019")

fig.suptitle("Discharge Comparison at Gauges for GloFAS v4.0", fontsize=15, fontweight='bold')

# Hide empty axes
for j in range(i+1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.show()

print(df_stats)


#%%
# Compare historical 30 years of wflow run and GloFAS v3.1 data and calculate KGE
print("Calculating KGE for 30-year comparison...")
kge_all = skillstats.kge(
    sim=glofas_30yr_aligned.chunk(dict(time=-1)),
    obs=Q_wflow_30yr_aligned['Q'].chunk(dict(time=-1)),
    dim="time"
)

# Convert to DataFrame
df_kge = pd.DataFrame({
    'gauge': kge_all["Q_gauges_locs"].values,
    'kge': kge_all["kge"].values
})


#%%
# Plot function for gauges comparison of 30yr wflow and GloFAS v3.1 data 
def plot_gauges_comparison(gauges=["1", "2"], kge_all=None,
                           wflow_ds=Q_wflow_30yr_aligned['Q'],
                           glofas_ds=glofas_30yr_aligned):
    # Prepare figure
    fig, axes = plt.subplots(len(gauges), 1, figsize=(10, 6), sharex=True)
    
    if len(gauges) == 1:
        axes = [axes]

    for i, g in enumerate(gauges):
        ax = axes[i]

        # Compute gauge timeseries if needed
        wflow_g = wflow_ds.sel(Q_gauges_locs=g).compute()
        glofas_g = glofas_ds.sel(Q_gauges_locs=g).compute()

        # Plot
        ax.plot(wflow_g.time.values, wflow_g.values, label="WFLOW", lw=1.5)
        ax.plot(glofas_g.time.values, glofas_g.values, label="GloFAS", lw=1.5)

        # KGE value
        if kge_all is not None:
            kge_val = float(kge_all["kge"].sel(Q_gauges_locs=g).compute().values)
            title = f"Gauge {g} — KGE = {kge_val:.2f}"
        else:
            title = f"Gauge {g}"
        ax.set_title(title, fontsize=12)

        # Label a/b
        ax.text(-0.05, 1.05, f"({chr(97+i)})", transform=ax.transAxes,
                fontsize=14, fontweight="bold")

        ax.set_ylabel("Discharge [m³/s]")
        ax.grid(alpha=0.3)
        ax.legend()

    # X-axis as years
    axes[-1].xaxis.set_major_locator(mdates.YearLocator(5))
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    axes[-1].set_xlabel("Year")

    fig.tight_layout()
    plt.show()

# Usage
plot_gauges_comparison(gauges=["1", "2"], kge_all=kge_all)


# %%
# Overview DataFrame with KGE and event volume comparison
print("Creating overview DataFrame with KGE and event volume comparison...")
df_overview_glofas = df_kge.copy(deep=True)
df_overview_glofas['wflow / glofas event discharge [%]'] = df_stats['wflow / GloFAS [%]']
df_overview_glofas['absolute difference [m³ * 1e8]'] = (df_stats['wflow [m³]'] - df_stats['GloFAS [m³]']) / 1e8
df_overview_glofas = df_overview_glofas.set_index('gauge')

df_overview_glofas
# %%

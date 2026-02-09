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
from hydromt.stats import skills as skillstats  # NSE, KGE, etc.
from hydromt import data_catalog
import string  
import matplotlib.dates as mdates
import matplotlib.patheffects as pe
from matplotlib.ticker import FuncFormatter
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import rioxarray

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
################# COMPARE WFLOW VS. GLOFAS DATA ###################
###################################################################
# Load data catalog
data_cats = ['../../Workflows/03_data_catalogs/datacatalog_general.yml']
data_cat = data_catalog.DataCatalog(data_cats)

# Event time range
start_date = "20190307 000000"
end_date = "20190325 000000"
start_dt = datetime.strptime(start_date, "%Y%m%d %H%M%S")
end_dt = datetime.strptime(end_date, "%Y%m%d %H%M%S")

# Load GloFAS datasets
glofas_ds_v31 = data_cat.get_rasterdataset('glofas_era5_v31', geom=gdf_wflow)
glofas_ds_v40 = data_cat.get_rasterdataset('glofas_era5_v40', geom=gdf_wflow)
glofas_ds_v40 = glofas_ds_v40.rename({"valid_time": "time"})

# --------------------------------------------------------------------
# Load upstream area datasets
# v4.0
glofas_uparea_v40 = xr.open_dataset("p:/11210471-001-compass/01_Data/glofas_v4/uparea_glofas_v4_0.nc")
glofas_uparea_v40 = glofas_uparea_v40.rio.write_crs("EPSG:4326", inplace=True)
glofas_uparea_v40_clipped = glofas_uparea_v40.rio.clip(gdf_wflow.geometry, gdf_wflow.crs)

# v3.1
glofas_uparea_v31 = rioxarray.open_rasterio("p:/wflow_global/hydromt/hydro/glofas_era5/glofas_v31_uparea.tif")
glofas_uparea_v31 = glofas_uparea_v31.rio.write_crs("EPSG:4326")  # ensure CRS is set
glofas_uparea_v31_clipped = glofas_uparea_v31.rio.clip(gdf_wflow.geometry, gdf_wflow.crs)

# Make sure both uparea datasets are xarray Datasets with proper names
glofas_uparea_v31_clipped.name = "uparea"

# --------------------------------------------------------------------
# Convert DataArray to Dataset if needed
if isinstance(glofas_uparea_v31_clipped, xr.DataArray):
    glofas_uparea_v31_ds = glofas_uparea_v31_clipped.to_dataset(name="uparea")
else:
    glofas_uparea_v31_ds = glofas_uparea_v31_clipped

# Rename dims to match v4.0
glofas_uparea_v31_ds = glofas_uparea_v31_ds.rename({"x": "longitude", "y": "latitude"})

# If lat/lon are 2D, reduce to 1D
if glofas_uparea_v31_ds.latitude.ndim == 2:
    glofas_uparea_v31_ds = glofas_uparea_v31_ds.assign_coords(
        latitude=glofas_uparea_v31_ds.latitude[:, 0],
        longitude=glofas_uparea_v31_ds.longitude[0, :]
    )

# Add spatial_ref coordinate
glofas_uparea_v31_ds = glofas_uparea_v31_ds.assign_coords(spatial_ref=0)

# Interpolate upstream area to GloFAS grid
glofas_uparea_interp_v40 = glofas_uparea_v40_clipped['uparea'].interp(latitude=glofas_ds_v40['latitude'], longitude=glofas_ds_v40['longitude'], method="nearest")

glofas_uparea_interp_v31 = glofas_uparea_v31_ds['uparea'].interp(latitude=glofas_ds_v31['latitude'], longitude=glofas_ds_v31['longitude'], method="nearest")

# --------------------------------------------------------------------
# Convert discharge DataArrays to Datasets and add uparea as a data variable
glofas_ds_v40 = glofas_ds_v40.to_dataset(name="discharge")
glofas_ds_v40["uparea"] = glofas_uparea_interp_v40

glofas_ds_v31 = glofas_ds_v31.to_dataset(name="discharge")
glofas_ds_v31["uparea"] = glofas_uparea_interp_v31



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
def combine_by_cluster(ds, gdf_wflow_gauges, rep_station):
    # Convert to DataFrame
    df = ds.to_dataframe().reset_index()
    
    # Merge cluster & rep_station info
    df = df.merge(
        gdf_wflow_gauges[["Q_gauges_locs", "cluster", "rep_station"]],
        on="Q_gauges_locs",
        how="left"
    )
    
    # Group by cluster and time, sum discharge
    df_comb = df.groupby(["time", "cluster"]).sum(numeric_only=True).reset_index()
    df_comb = df_comb.merge(rep_station, on="cluster")
    ds_comb = df_comb.set_index(["time", "rep_station"]).to_xarray().rename({"rep_station": "Q_gauges_locs"})
    
    # 6. Assign lat/lon coordinates from representative gauges
    reps = gdf_wflow_gauges.set_index("Q_gauges_locs").loc[ds_comb.Q_gauges_locs.values]
    ds_comb = ds_comb.assign_coords(
        lon=("Q_gauges_locs", reps.geometry.x.values),
        lat=("Q_gauges_locs", reps.geometry.y.values)
    )
    
    return ds_comb

# Combine by cluster and keep coordinates
Q_wflow_30yr_filt_combined = combine_by_cluster(Q_wflow_30yr_filtered, gdf_wflow_gauges=gdf_wflow_gauges, rep_station=rep_station)
Q_wflow_event_filt_combined = combine_by_cluster(Q_wflow_event_filtered, gdf_wflow_gauges=gdf_wflow_gauges, rep_station=rep_station)

gdf_rep = (gdf_wflow_gauges[gdf_wflow_gauges.Q_gauges_locs == gdf_wflow_gauges.rep_station].drop_duplicates("cluster"))

gdf_rep = gdf_rep.head(2)
# gdf_rep.to_file("C:/Code/Paper_1/QGIS/wflow_filt_combined_gauges_v3.geojson", driver="GeoJSON")

# Only select the top three combined gauges
Q_wflow_30yr_filt_combined = Q_wflow_30yr_filt_combined.sel(Q_gauges_locs=gdf_rep.head(3).rep_station.values)
Q_wflow_event_filt_combined = Q_wflow_event_filt_combined.sel(Q_gauges_locs=gdf_rep.head(3).rep_station.values)

# Resample to daily means to match GloFAS temporal resolution
Q_wflow_event_filt_resampled = Q_wflow_event_filt_combined.resample(time="1D").mean()

#%%
print("Selecting GloFAS stations closest to specified gauges...")
gauges_v4 = pd.DataFrame({
    "Q_gauges_locs": ["1", "2"],
    "latitude": [-19.9318, -19.2229],
    "longitude": [34.2329, 34.5354],
    "uparea": [2.674231e10, 2.5425527e10],
}).set_index("Q_gauges_locs")

discharge_sel_v4 = []
uparea_sel_v4 = []

for gauge, row in gauges_v4.iterrows():
    ds_point = glofas_ds_v40.sel(latitude=row.latitude, longitude=row.longitude, method="nearest")
    discharge_sel_v4.append(ds_point["discharge"])
    uparea_sel_v4.append(ds_point["uparea"])

glofas_ds_v40_sel = xr.Dataset(data_vars={"discharge": (("Q_gauges_locs", "time"), xr.concat(discharge_sel_v4, dim="Q_gauges_locs").data),
                                          "uparea": (("Q_gauges_locs",), xr.concat(uparea_sel_v4, dim="Q_gauges_locs").data),},
                                          coords={"Q_gauges_locs": gauges_v4.index.values, "time": glofas_ds_v40.time,
                                                  "latitude": ("Q_gauges_locs", gauges_v4.latitude.values), 
                                                  "longitude": ("Q_gauges_locs", gauges_v4.longitude.values)})


gauges_v31 = pd.DataFrame({
    "Q_gauges_locs": ["1", "2", "3"],
    "latitude": [-19.985, -19.143, -20.278],
    "longitude": [34.24, 34.558, 34.544],
    "uparea": [2.71818e10, 2.5136839e10, 2.8109243e9],
}).set_index("Q_gauges_locs")

discharge_sel_v31 = []
uparea_sel_v31 = []

for gauge, row in gauges_v31.iterrows():
    ds_point = glofas_ds_v31.sel(latitude=row.latitude, longitude=row.longitude, method="nearest")
    discharge_sel_v31.append(ds_point["discharge"])
    uparea_sel_v31.append(ds_point["uparea"])

glofas_ds_v31_sel = xr.Dataset(data_vars={"discharge": (("Q_gauges_locs", "time"), xr.concat(discharge_sel_v31, dim="Q_gauges_locs").data),
                                          "uparea": (("Q_gauges_locs",), xr.concat([ua.squeeze() for ua in uparea_sel_v31], dim="Q_gauges_locs").data),},
                                          coords={"Q_gauges_locs": gauges_v31.index.values, "time": glofas_ds_v31.time,
                                                  "latitude": ("Q_gauges_locs", gauges_v31.latitude.values), 
                                                  "longitude": ("Q_gauges_locs", gauges_v31.longitude.values)})


# Align time series
Q_wflow_30yr_aligned, glofas_30yr_aligned_v40 = xr.align(Q_wflow_30yr_filt_combined, glofas_ds_v40_sel, join="inner")
Q_wflow_event_aligned, glofas_event_aligned_v40 = xr.align(Q_wflow_event_filt_resampled, glofas_ds_v40_sel, join="inner")

Q_wflow_30yr_aligned, glofas_30yr_aligned_v31 = xr.align(Q_wflow_30yr_filt_combined, glofas_ds_v31_sel, join="inner")
Q_wflow_event_aligned, glofas_event_aligned_v31 = xr.align(Q_wflow_event_filt_resampled, glofas_ds_v31_sel, join="inner")

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
fig, axes = plt.subplots(nrows, ncols, figsize=(10, nrows*3), sharex=True, sharey=False)
axes = axes.flatten()

letters = list(string.ascii_lowercase)

for i, g in enumerate(gauge_ids):
    ax = axes[i]
    print(g)
    # Select timeseries
    wflow_ts = Q_wflow_event_aligned.sel(Q_gauges_locs=g, time=slice(start_dt, end_dt))
    glofas_ts = glofas_event_aligned_v40.sel(Q_gauges_locs=g, time=slice(start_dt, end_dt))
    
    # Plot
    ax.plot(wflow_ts.time, wflow_ts['Q'], label="WFLOW", lw=1.5)
    ax.plot(glofas_ts.time, glofas_ts['discharge'], label="GloFAS", lw=1.5, alpha=0.7)
    
    # Time differences in seconds
    dt_wflow = np.diff(wflow_ts.time.values).astype("timedelta64[s]").astype(float)
    dt_wflow = np.append(dt_wflow, dt_wflow[-1])
    
    dt_glofas = np.diff(glofas_ts.time.values).astype("timedelta64[s]").astype(float)
    dt_glofas = np.append(dt_glofas, dt_glofas[-1])
    
    # Total volume
    total_Q = (wflow_ts['Q'].values * dt_wflow).sum()
    total_G = (glofas_ts['discharge'].values * dt_glofas).sum()
    
    # Save stats
    df_stats.loc[i] = [g, total_Q, total_G, total_Q / total_G * 100]
    
    # Plot styling
    axes[0].set_title(f"Buzi river (gauge 1)")
    axes[1].set_title(f"Pungwe river (gauge 2)")
    ax.set_ylabel("Discharge [m³/s]")
    ax.grid(alpha=0.3)
    if i == 0:
        ax.legend()
    # Letter annotation
    ax.text(0.0, 1.1, f'({letters[i]})', transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left')
    
    # Formatting x-axis as days in 2019
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))  # every 2 days
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d"))
    ax.tick_params(axis='x')
    ax.set_xlabel(" Day in March 2019")

fig.suptitle("Discharge Comparison at Gauges for GloFAS v4.0", fontsize=13, fontweight='bold')

# Hide empty axes
for j in range(i+1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.show()

print(df_stats)



#%%
# Compare historical 30 years of wflow run and GloFAS v4.0 data and calculate KGE
print("Calculating KGE for 30-year comparison...")
kge_all_v40 = skillstats.kge(
    sim=glofas_30yr_aligned_v40['discharge'].chunk(dict(time=-1)),
    obs=Q_wflow_30yr_aligned['Q'].chunk(dict(time=-1)),
    dim="time"
)

# Convert to DataFrame
df_kge_v40 = pd.DataFrame({
    'gauge': kge_all_v40["Q_gauges_locs"].values,
    'kge': kge_all_v40["kge"].values
})


#%%
# Plot function for gauges comparison of 30yr wflow and GloFAS v4.0 data 
def plot_gauges_comparison(gauges=["1", "2"], kge_all=None, wflow_ds=Q_wflow_30yr_aligned['Q'], glofas_ds=None):
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
plot_gauges_comparison(gauges=["1", "2"], kge_all=kge_all_v40, glofas_ds=glofas_30yr_aligned_v40['discharge'])


# %%
# Overview DataFrame with KGE and event volume comparison
print("Creating overview DataFrame with KGE and event volume comparison...")
df_overview_glofas_v40 = df_kge_v40.copy(deep=True)
df_overview_glofas_v40['wflow / glofas event discharge [%]'] = df_stats['wflow / GloFAS [%]']
df_overview_glofas_v40['absolute difference [m³ * 1e8]'] = (df_stats['wflow [m³]'] - df_stats['GloFAS [m³]']) / 1e8
df_overview_glofas_v40 = df_overview_glofas_v40.set_index('gauge')

df_overview_glofas_v40


# %%
#### same but for GloFAS v3.1 compared to GloFAS v4.0 ####

# Plot all gauges time series for the event
print("Plotting time series comparison for all gauges...")

# Get all gauge IDs
gauge_ids = Q_wflow_event_aligned['Q_gauges_locs'].values

# Prepare DataFrame to store stats
df_stats_31 = pd.DataFrame(columns=['Gauge #', 'wflow [m³]', 'GloFAS [m³]', 'wflow / GloFAS [%]', 'GloFAS / wflow [%]', '(GloFAS - wflow) / wflow [%]'])
df_stats_40 = pd.DataFrame(columns=['Gauge #', 'wflow [m³]', 'GloFAS [m³]', 'wflow / GloFAS [%]', 'GloFAS / wflow [%]', '(GloFAS - wflow) / wflow [%]'])

# Plotting setup
n = len(gauge_ids)
ncols = 2
nrows = int(np.ceil(n / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(10, nrows*3), sharex=True, sharey=False)
axes = axes.flatten()

letters = list(string.ascii_lowercase)

for i, g in enumerate(gauge_ids):
    ax = axes[i]
    print(g)
    # Select timeseries
    wflow_ts = Q_wflow_event_aligned.sel(Q_gauges_locs=g, time=slice(start_dt, end_dt))
    glofas_ts_v40 = glofas_event_aligned_v40.sel(Q_gauges_locs=g, time=slice(start_dt, end_dt))
    glofas_ts_v31 = glofas_event_aligned_v31.sel(Q_gauges_locs=g, time=slice(start_dt, end_dt))
    
    # Plot
    ax.plot(wflow_ts.time, wflow_ts['Q'], label="wflow", lw=1.5)
    ax.plot(glofas_ts_v40.time, glofas_ts_v40['discharge'], label="GloFAS v4.0", lw=1.5, alpha=0.7)
    ax.plot(glofas_ts_v31.time, glofas_ts_v31['discharge'], label="GloFAS v3.1", lw=1.5, alpha=0.7)
    
    # Time differences in seconds
    dt_wflow = np.diff(wflow_ts.time.values).astype("timedelta64[s]").astype(float)
    dt_wflow = np.append(dt_wflow, dt_wflow[-1])
    
    dt_glofas_v40 = np.diff(glofas_ts_v40.time.values).astype("timedelta64[s]").astype(float)
    dt_glofas_v40 = np.append(dt_glofas_v40, dt_glofas_v40[-1])
    
    dt_glofas_v31 = np.diff(glofas_ts_v31.time.values).astype("timedelta64[s]").astype(float)
    dt_glofas_v31 = np.append(dt_glofas_v31, dt_glofas_v31[-1])
    
    # Total volume
    total_Q = (wflow_ts['Q'].values * dt_wflow).sum()
    total_G_v40 = (glofas_ts_v40['discharge'].values * dt_glofas_v40).sum()
    total_G_v31 = (glofas_ts_v31['discharge'].values * dt_glofas_v31).sum()
    
    # Save stats
    df_stats_40.loc[i] = [g, total_Q, total_G_v40, total_Q / total_G_v40 * 100, total_G_v40 / total_Q * 100, (total_G_v40 - total_Q) / total_Q * 100]
    df_stats_31.loc[i] = [g, total_Q, total_G_v31, total_Q / total_G_v31 * 100, total_G_v31 / total_Q * 100, (total_G_v31 - total_Q) / total_Q * 100]
    
    # Plot styling
    axes[0].set_title(f"Buzi river (gauge 1)")
    axes[1].set_title(f"Pungwe river (gauge 2)")
    ax.set_ylabel("Discharge [m³/s]")
    ax.grid(alpha=0.3)
    if i == 0:
        ax.legend()
    # Letter annotation
    ax.text(0.0, 1.1, f'({letters[i]})', transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left')
    
    # Formatting x-axis as days in 2019
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))  # every 2 days
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d"))
    ax.tick_params(axis='x')
    ax.set_xlabel(" Day in March 2019")
    ax.set_xlim(start_dt, end_dt)

axes[0].legend()

fig.suptitle("Discharge comparison of wflow with GloFAS during TC Idai", fontsize=13, fontweight='bold')
fig.savefig("../figures/fS5.png", dpi=300, bbox_inches="tight")

# Hide empty axes
for j in range(i+1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.show()

print(df_stats_40)
print(df_stats_31)



# %%
# Compare historical 30 years of wflow run and GloFAS v3.1 data and calculate KGE
print("Calculating KGE for 30-year comparison...")
kge_all_v31 = skillstats.kge(
    sim=glofas_30yr_aligned_v31['discharge'].chunk({'time': -1}),
    obs=Q_wflow_30yr_aligned['Q'].chunk({'time': -1}),
    dim="time"
)

# Convert to DataFrame
df_kge_v31 = pd.DataFrame({
    'gauge': kge_all_v31["Q_gauges_locs"].values,
    'kge': kge_all_v31["kge"].values
})


#%%
def plot_gauges_comparison_versions(gauges=["1", "2"], kge_all_v40=None, kge_all_v31=None, wflow_ds=Q_wflow_30yr_aligned['Q'], glofas_ds_v40=None, glofas_ds_v31=None):
    # Prepare figure
    fig, axes = plt.subplots(len(gauges), 1, figsize=(10, 6), sharex=True)
    
    if len(gauges) == 1:
        axes = [axes]

    for i, g in enumerate(gauges):
        ax = axes[i]

        # Compute gauge timeseries if needed
        wflow_g = wflow_ds.sel(Q_gauges_locs=g).compute()
        glofas_g_v40 = glofas_ds_v40.sel(Q_gauges_locs=g).compute()
        glofas_g_v31 = glofas_ds_v31.sel(Q_gauges_locs=g).compute()

        # Plot
        ax.plot(wflow_g.time.values, wflow_g.values, label="wflow", lw=1.5, zorder=1)
        ax.plot(glofas_g_v40.time.values, glofas_g_v40.values, label="GloFAS v4.0", lw=1.5, zorder=3)
        ax.plot(glofas_g_v31.time.values, glofas_g_v31.values, label="GloFAS v3.1", lw=1.5, zorder=2)

        # KGE value
        if kge_all_v40 is not None and kge_all_v31 is not None:
            kge_val_v40 = float(kge_all_v40["kge"].sel(Q_gauges_locs=g).compute().values)
            kge_val_v31 = float(kge_all_v31["kge"].sel(Q_gauges_locs=g).compute().values)
            title = f"Gauge {g} — KGE v4.0 = {kge_val_v40:.2f}, KGE v3.1 = {kge_val_v31:.2f}"
        else:
            title = f"Gauge {g}"
        ax.set_title(title, fontsize=12)

        # Label a/b
        ax.text(-0.05, 1.05, f"({chr(97+i)})", transform=ax.transAxes,
                fontsize=14, fontweight="bold")

        ax.set_ylabel("Discharge [m³/s]")
        ax.grid(alpha=0.3)
    
    axes[0].legend()

    # X-axis as years
    axes[-1].xaxis.set_major_locator(mdates.YearLocator(5))
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    axes[-1].set_xlabel("Year")

    fig.tight_layout()
    plt.show()

# Usage
plot_gauges_comparison_versions(gauges=["1", "2"], kge_all_v40=kge_all_v40, kge_all_v31=kge_all_v31, glofas_ds_v40=glofas_30yr_aligned_v40['discharge'], glofas_ds_v31=glofas_30yr_aligned_v31['discharge'])


# #%%
# # Plot function for gauges comparison of 30yr wflow and GloFAS v3.1 data 
# plot_gauges_comparison(gauges=["1", "2"], kge_all=kge_all_v31, glofas_ds=glofas_30yr_aligned_v31['discharge'])


# # %%
# # Overview DataFrame with KGE and event volume comparison
# print("Creating overview DataFrame with KGE and event volume comparison...")
# df_overview_glofas_v31 = df_kge_v31.copy(deep=True)
# df_overview_glofas_v31['wflow / glofas event discharge [%]'] = df_stats['wflow / GloFAS [%]']
# df_overview_glofas_v31['absolute difference [m³ * 1e8]'] = (df_stats['wflow [m³]'] - df_stats['GloFAS [m³]']) / 1e8
# df_overview_glofas_v31 = df_overview_glofas_v31.set_index('gauge')

# df_overview_glofas_v31
# # %%

# %%

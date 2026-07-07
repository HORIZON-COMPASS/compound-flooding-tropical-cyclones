#%% script from Natalia Aleksandrove and adapted by Doris Vertegaal to track TCs from ClimateDT data
import platform
import numpy as np
import xarray as xr
import os
import geopandas as gpd
import re
import glob
from rasterio.features import rasterize
import rasterio as rio
import hydromt
import numpy as np
import pandas as pd
import geopandas as gpd
import os
import pickle
from scipy.ndimage import minimum_filter
from shapely.geometry import Point

#%%
def preprocess_add_r(ds):
    source = ds.encoding.get('source', '')
    match = re.search(r'_r(\d+)_', source)
    # realization = int(match.group(1)) if match else -1
    # ds = ds.assign_coords(realization=realization)
    return ds

prefix = "p:/" if platform.system() == "Windows" else "/p/"

# Data directory
datadir = os.path.join(prefix, '11210471-001-compass/01_Data/ECMWF_ClimateDT/data/preprocessed/Gen2/Idai_full')

# # Data catalog
# data_cat = ['../data_catalogs/datacatalog_general.yml', '../data_catalogs/datacatalog_CF_forcing.yml'] 
# data_catalog = hydromt.data_catalog.DataCatalog(data_libs = data_cat)

# Spatial and temporal domain
# bbox_wflow_sfincs = [32,-20.6,35.6,-17.5]
# bbox_idai = [30, -28, 45, -9] 
bbox_idai = [30, -22, 45, -12] 
# areaname = 'wflow+sfincs'

time_min = "2019-03-1"
time_max = "2019-03-25"

# gdf_sfincs = gpd.read_file(r"p:\11210471-001-compass\02_Models\sofala\Idai\sfincs\gis\region.geojson")
# gdf_sfincs = gdf_sfincs.to_crs(epsg=4326) 

# gdf_wflow = gpd.read_file(r"p:\11210471-001-compass\02_Models\sofala\Idai\wflow\staticgeoms\basins.geojson")
# gdf_wflow = gdf_wflow.to_crs(epsg=4326) 



#%%
##### DestinE ClimateDT storyline data #####
# Reading in the data for the control and historical runs, and selecting the spatial domain of interest
# ds_cont_raw = xr.open_mfdataset(os.path.join(datadir, 'climateDT_*_cont_*_Idai_full_20190301_to_20190329.nc'),
#                             preprocess=preprocess_add_r)
# ds_cont_large = ds_cont_raw.sel(latitude=slice(bbox_idai[1], bbox_idai[3]),
#                       longitude=slice(bbox_idai[0], bbox_idai[2]))

ds_hist_raw = xr.open_mfdataset(os.path.join(datadir, 'climateDT_*_hist_*_Idai_full_20190301_to_20190329.nc'),
                            preprocess=preprocess_add_r)
ds_hist_large = ds_hist_raw.sel(latitude=slice(bbox_idai[1], bbox_idai[3]),
                      longitude=slice(bbox_idai[0], bbox_idai[2]))


#%%
##### Observational data #######
# Reading in ERA5 data for the same time period and spatial domain
# ds_era5_large = data_catalog.get_rasterdataset(
#     data_like = "era5_hourly_zarr",
#     time_tuple = ("2019-03-01","2019-03-29"),
#     variables = ['precip', 'wind10_u', 'wind10_v', 'press_msl'],
#     bbox=bbox_idai, 
#     buffer = 1
# )

#%%
# Creating TC track from raster data
# ---------------------------
# SETTINGS
# ---------------------------
MSLP_THRESHOLD = 105000  # Pa (same threshold as IBTrACS)
raw_track_file = "tracks_all_realizations.pkl"

# ---------------------------
# FUNCTIONS
# ---------------------------
def find_mslp_minima(mslp_2d):
    return (mslp_2d == minimum_filter(mslp_2d, size=3, mode='nearest'))

def fast_distance_km(lat1, lon1, lat2, lon2):
    return 111 * np.sqrt((lat1 - lat2)**2 + (lon1 - lon2)**2)

# ---------------------------
# PRECOMPUTE WIND
# ---------------------------
print("Precomputing wind...")
ds_hist_large['wind'] = np.sqrt(ds_hist_large['10u']**2 + ds_hist_large['10v']**2)

tracks_all_realizations = {}

lat = ds_hist_large.latitude.values
lon = ds_hist_large.longitude.values

# ---------------------------
# LOOP PER REALIZATION
# ---------------------------
for r in ds_hist_large.realization.values:
    print(f"Processing realization {r}...")

    ds_r = ds_hist_large.sel(realization=r)

    mslp_array = ds_r['msl']
    wind_array = ds_r['wind']

    # Cache wind per timestep
    wind_cache = {t: wind_array.sel(time=t).values for t in ds_r.time.values}

    # ---------------------------
    # FIND MINIMA
    # ---------------------------
    candidates = []

    for t in ds_r.time.values:

        mslp = mslp_array.sel(time=t).values

        minima_mask = find_mslp_minima(mslp)
        yy, xx = np.where(minima_mask)

        for i, j in zip(yy, xx):
            # Optional filtering
            if mslp[i, j] > MSLP_THRESHOLD:
                continue

            candidates.append({
                "time": t,
                "lat": lat[i],
                "lon": lon[j],
                "i": i,
                "j": j,
                "mslp": float(mslp[i, j]),
                "realization": r
            })

    print(f"  Candidates found: {len(candidates)}")

    # ---------------------------
    # TRACK BUILDING
    # ---------------------------
    tracks = []

    for point in candidates:
        assigned = False

        for track in tracks:
            last = track[-1]

            dt = (point["time"] - last["time"]) / np.timedelta64(1, 'h')

            if dt == 1:
                if fast_distance_km(point["lat"], point["lon"], last["lat"], last["lon"]) < 300:
                    track.append(point)
                    assigned = True
                    break

        if not assigned:
            tracks.append([point])

    print(f"  Tracks before filtering: {len(tracks)}")

    # ---------------------------
    # ADD WIND
    # ---------------------------
    for track in tracks:
        for pt in track:
            wind_t = wind_cache[pt["time"]]
            pt["wind"] = float(wind_t[pt["i"], pt["j"]])

    # ---------------------------
    # FILTER TRACKS
    # ---------------------------
    valid_tracks = []
    invalid_tracks = []

    for track in tracks:

        winds = [pt["wind"] for pt in track]
        mslps = [pt["mslp"] for pt in track]

        duration_hours = len(track)
        max_wind = np.max(winds)
        pressure_drop = np.max(mslps) - np.min(mslps)

        if (
            max_wind >= 25 and
            # pressure_drop >= 27 and #TODO: should be >27 below the yearly maximum over its lifetime
            duration_hours >= 72
        ):
            valid_tracks.append(track)
        else:
            invalid_tracks.append(track)

    print(f"  Valid tracks: {len(valid_tracks)}")
    print(f"  Invalid tracks: {len(invalid_tracks)}")

    # Store per realization
    tracks_all_realizations[r] = {
        "valid": valid_tracks,
        "invalid": invalid_tracks
    }

print("Tracking complete!")

#%%
# Load TC Idai track from shapefile and convert to list of (lat, lon) tuples
data_base = os.path.join(prefix, "11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS")

shapefile_path = os.path.join(data_base, "IBTrACS_IDAI.shp")
gdf = gpd.read_file(shapefile_path)
tc_idai = gdf[gdf['SID'] == '2019063S18038']

tc_idai_points = list(zip(tc_idai.geometry.y, tc_idai.geometry.x))


def track_distance(track, ref_points):
    track_points = [(pt["lat"], pt["lon"]) for pt in track]

    # match shortest length
    n = min(len(track_points), len(ref_points))
    if n < 10:
        return np.inf  # too short, discard

    dists = []
    for i in range(n):
        lat1, lon1 = track_points[i]
        lat2, lon2 = ref_points[i]

        d = fast_distance_km(lat1, lon1, lat2, lon2)
        dists.append(d)

    return np.mean(dists)

best_matches_per_realization = {}

for r, data in tracks_all_realizations.items():
    all_tracks = data["valid"] + data["invalid"]
    matches = []

    for tid, track in enumerate(all_tracks):
        dist = track_distance(track, tc_idai_points)

        if np.isfinite(dist):
            matches.append({
                "realization": int(r),
                "track_id": tid,
                "distance": dist,
                "track": track
            })

    # sort within this realization
    matches = sorted(matches, key=lambda x: x["distance"])

    # keep top 5 for *this* realization
    best_matches_per_realization[r] = matches[:5]

rows = []

for r, matches in best_matches_per_realization.items():
    for rank, item in enumerate(matches):
        track = item["track"]

        for pt in track:
            rows.append({
                "realization": r,
                "track_rank": rank,          # rank within realization
                "distance": item["distance"],
                "time": pt["time"],
                "lat": pt["lat"],
                "lon": pt["lon"],
                "geometry": Point(pt["lon"], pt["lat"]),
                "wind": pt["wind"],
                "mslp": pt["mslp"]
            })

gdf_best = gpd.GeoDataFrame(rows, crs="EPSG:4326")
gdf_best.to_file("tracks_best_per_realization.gpkg", driver="GPKG")


from shapely.geometry import LineString

rows_lines = []

for r, matches in best_matches_per_realization.items():
    for rank, item in enumerate(matches):

        coords = [(pt["lon"], pt["lat"]) for pt in item["track"]]

        rows_lines.append({
            "realization": r,
            "track_rank": rank,
            "distance": item["distance"],
            "geometry": LineString(coords)
        })

gdf_lines = gpd.GeoDataFrame(rows_lines, crs="EPSG:4326")
gdf_lines.to_file("tracks_best_lines.gpkg", driver="GPKG")



#%%
# ---------------------------
# SAVE RAW TRACKS
# ---------------------------
with open(raw_track_file, "wb") as f:
    pickle.dump(tracks_all_realizations, f)

# ---------------------------
# FLATTEN TO TABLE
# ---------------------------
rows = []

for r, data in tracks_all_realizations.items():

    for is_valid, tracks in zip(
        [True, False],
        [data["valid"], data["invalid"]]
    ):
        for tid, track in enumerate(tracks):
            for pt in track:
                rows.append({
                    "realization": r,
                    "track_id": tid,
                    "time": pt["time"],
                    "lat": pt["lat"],
                    "lon": pt["lon"],
                    "i": pt["i"],
                    "j": pt["j"],
                    "wind": pt["wind"],
                    "mslp": pt["mslp"],
                    "is_valid": is_valid
                })

# ✅ FAST save (recommended)
df_tracks = pd.DataFrame(rows)
df_tracks.to_parquet("tracks.parquet")

# ---------------------------
# (OPTIONAL) SAVE GEO VERSION
# ---------------------------
gdf_tracks = gpd.GeoDataFrame(
    df_tracks,
    geometry=gpd.points_from_xy(df_tracks.lon, df_tracks.lat),
    crs="EPSG:4326"
)

gdf_tracks.to_parquet("tracks_geo.parquet")

print("✅ Tracks saved (raw + parquet + geo)!")
# %%

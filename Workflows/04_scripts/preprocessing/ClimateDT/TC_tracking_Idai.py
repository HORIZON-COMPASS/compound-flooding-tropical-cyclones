'''
Script to track TCs in the storyline data from ClimateDT data.

Author: dvertegaal
'''

#%% script from Natalia Aleksandrove and adapted by Doris Vertegaal to track TCs from ClimateDT data
import platform
import numpy as np
import xarray as xr
import os
import geopandas as gpd
import re
import numpy as np
import pandas as pd
import geopandas as gpd
import os
import pickle
from scipy.ndimage import minimum_filter
from shapely.geometry import Point
from shapely.geometry import LineString

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

# Spatial and temporal domain 
bbox_idai = [30, -22, 45, -12] 

time_min = "2019-03-1"
time_max = "2019-03-25"


# Load TC Idai track from shapefile and convert to list of (lat, lon) tuples
data_base = os.path.join(prefix, "11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS")

shapefile_path = os.path.join(data_base, "IBTrACS_IDAI.shp")
gdf = gpd.read_file(shapefile_path)
tc_idai = gdf[gdf['SID'] == '2019063S18038']

tc_idai_points = list(zip(tc_idai.geometry.y, tc_idai.geometry.x))

# Point from where Idai was a TC1
first_tc1_idx = tc_idai[tc_idai['USA_SSHS'] == 1].index[0]
tc_idai_tc1 = tc_idai.loc[first_tc1_idx:]
tc_idai_points_TC = list(zip(tc_idai_tc1.geometry.y, tc_idai_tc1.geometry.x))


#%%
##### DestinE ClimateDT storyline data #####
# Present-day
ds_hist = xr.open_mfdataset(os.path.join(datadir, 'climateDT_*_hist_*_Idai_full_20190301_to_20190329.nc'),
                            preprocess=preprocess_add_r)


#%%
# Creating TC track from raster data
# ---------------------------
# SETTINGS
# ---------------------------
# MSLP_THRESHOLD = 105000  # Pa (same threshold as IBTrACS)
# raw_track_file = "tracks_all_realizations.pkl"

# # ---------------------------
# # FUNCTIONS
# # ---------------------------
# def find_mslp_minima(mslp_2d):
#     return (mslp_2d == minimum_filter(mslp_2d, size=3, mode='nearest'))

# def fast_distance_km(lat1, lon1, lat2, lon2):
#     return 111 * np.sqrt((lat1 - lat2)**2 + (lon1 - lon2)**2)

# # ---------------------------
# # PRECOMPUTE WIND
# # ---------------------------
# print("Precomputing wind...")
# ds_hist['wind'] = np.sqrt(ds_hist['10u']**2 + ds_hist['10v']**2)

# tracks_all_realizations = {}

# lat = ds_hist.latitude.values
# lon = ds_hist.longitude.values

# # ---------------------------
# # LOOP PER REALIZATION
# # ---------------------------
# for r in ds_hist.realization.values:
#     print(f"Processing realization {r}...")

#     ds_r = ds_hist.sel(realization=r)

#     mslp_array = ds_r['msl']
#     wind_array = ds_r['wind']

#     # Cache wind per timestep
#     wind_cache = {t: wind_array.sel(time=t).values for t in ds_r.time.values}

#     # ---------------------------
#     # FIND MINIMA
#     # ---------------------------    
#     times = ds_r.time.values
#     mslp_values = mslp_array.values

#     candidates = []
#     for it, t in enumerate(times):

#         mslp = mslp_values[it]

#         minima_mask = find_mslp_minima(mslp)

#         yy, xx = np.where(minima_mask & (mslp <= MSLP_THRESHOLD))

#         for i, j in zip(yy, xx):
#             candidates.append({
#                 "time": t,
#                 "lat": lat[i],
#                 "lon": lon[j],
#                 "i": i,
#                 "j": j,
#                 "mslp": float(mslp[i, j]),
#                 "realization": r
#             })

#     print(f"  Candidates found: {len(candidates)}")

#     # ---------------------------
#     # TRACK BUILDING
#     # ---------------------------
#     tracks = []
#     active_tracks = []

#     # group candidates by time
#     candidates_by_time = {}
#     for c in candidates:
#         candidates_by_time.setdefault(c["time"], []).append(c)

#     times = sorted(candidates_by_time.keys())

#     for t in times:
#         current_candidates = candidates_by_time[t]
#         used = set()
#         new_active_tracks = []

#         # try to extend existing tracks
#         for track in active_tracks:      
#             last = track[-1]
#             dt = (t - last["time"]) / np.timedelta64(1, "h")
    
#             best_idx = None
#             best_cand = None
#             best_score = None

#             # skip if candidates are further than 1 h apart
#             if dt != 1:
#                 continue

#             # best = None
#             for idx, cand in enumerate(current_candidates): 
#                 # Quick checks to skip candidates that are already used or too far away
#                 if idx in used:
#                     continue
#                 if abs(cand["lat"] - last["lat"]) > 0.5:
#                     continue
#                 if abs(cand["lon"] - last["lon"]) > 0.5:
#                     continue
                
#                 # Compute distance between candidates in km
#                 dist = fast_distance_km(cand["lat"], cand["lon"], last["lat"], last["lon"])

#                 # Max translation speed of Idai was 11 knots (20 km/h), so 50 km is reasonable
#                 if dist < 50: # km                    
#                     score = (cand["mslp"], dist)
#                     if best_score is None or score < best_score:
#                         best_score = score
#                         best_idx = idx
#                         best_cand = cand

#             if best_cand is not None:    
#                 track.append(best_cand)
#                 used.add(best_idx)
#                 new_active_tracks.append(track)
#                 print(f"  Extended track: " f"mslp={best_cand['mslp']:.0f}, "f"dist={best_score[1]:.1f} km"            )

#         # start new tracks from unused candidates
#         for idx, cand in enumerate(current_candidates):
#             if idx in used:
#                 continue
#             track = [cand]
#             tracks.append(track)
#             new_active_tracks.append(track)

#         active_tracks = new_active_tracks

#     print(f"  Tracks before filtering: {len(tracks)}")
    
#     # ---------------------------
#     # ADD WIND
#     # ---------------------------    
#     # RMW from IBTrACS (miles) with a 10% buffer
#     rmw_miles = tc_idai['USA_RMW'].max()
#     search_radius_km = rmw_miles * 1.60934 * 2

#     print(f"Searching max wind within {search_radius_km:.1f} km")
    
#     # Calculate the max wind speed within the 1.1* radius of max wind speed of Idai for each point in the track
#     # approximate grid spacing in km
#     dlat_km = 111 * float(np.mean(np.diff(lat)))
#     dlon_km = 111 * np.cos(np.deg2rad(np.mean(lat))) * float(np.mean(np.diff(lon)))

#     window = 30

#     ii, jj = np.meshgrid(np.arange(-window, window + 1), np.arange(-window, window + 1),
#                          indexing="ij")
#     dist_km = np.sqrt((ii * dlat_km)**2 + (jj * dlon_km)**2)
#     dist_mask = dist_km <= search_radius_km

#     for track in tracks:
#         for pt in track:
#             wind_t = wind_cache[pt["time"]]

#             i0 = pt["i"]
#             j0 = pt["j"]

#             i1 = max(0, i0 - window)
#             i2 = min(len(lat), i0 + window + 1)

#             j1 = max(0, j0 - window)
#             j2 = min(len(lon), j0 + window + 1)

#             subset = wind_t[i1:i2, j1:j2]

#             mask = dist_mask[(i1 - i0 + window):(i2 - i0 + window),
#                              (j1 - j0 + window):(j2 - j0 + window)]

#             pt["wind"] = float(np.nanmax(subset[mask]))

#     # ---------------------------
#     # FILTER TRACKS
#     # ---------------------------
#     valid_tracks = []
#     invalid_tracks = []

#     for track in tracks:
#         winds = [pt["wind"] for pt in track]
#         mslps = [pt["mslp"] for pt in track]

#         duration_hours = len(track)
#         max_wind = np.max(winds)
#         pressure_drop = np.max(mslps) - np.min(mslps)

#         if (
#             max_wind >= 25 and
#             # pressure_drop >= 27 and #TODO: should be >27 below the yearly maximum over its lifetime
#             duration_hours >= 72
#         ):
#             valid_tracks.append(track)
#         else:
#             invalid_tracks.append(track)

#     print(f"  Valid tracks: {len(valid_tracks)}")
#     print(f"  Invalid tracks: {len(invalid_tracks)}")

#     # Store per realization
#     tracks_all_realizations[r] = {
#         "valid": valid_tracks,
#         "invalid": invalid_tracks
#     }

# print("Tracking complete!")

# #%%
# # ----------------------------------
# # SAVE BEST MATCHES PER REALIZATION
# # ----------------------------------
# def track_distance(track, ref_points):
#     track_points = [(pt["lat"], pt["lon"]) for pt in track]

#     # match shortest length
#     n = min(len(track_points), len(ref_points))
#     if n < 10:
#         return np.inf  # too short, discard

#     dists = []
#     for i in range(n):
#         lat1, lon1 = track_points[i]
#         lat2, lon2 = ref_points[i]

#         d = fast_distance_km(lat1, lon1, lat2, lon2)
#         dists.append(d)

#     return np.mean(dists)

# best_match_per_realization = {}

# for r, data in tracks_all_realizations.items():
#     all_tracks = data["valid"] + data["invalid"]
    
#     best_match = None
#     best_distance = np.inf

#     for tid, track in enumerate(all_tracks):

#         dist = track_distance(track, tc_idai_points_TC)

#         if np.isfinite(dist) and dist < best_distance:
#             best_distance = dist
#             best_match = {
#                 "realization": int(r),
#                 "track_id": tid,
#                 "distance": dist,
#                 "track": track
#             }

#     if best_match is not None:
#         best_match_per_realization[r] = best_match


# rows = []
# for r, item in best_match_per_realization.items():
#     track = item["track"]

#     for pt in track:
#         rows.append({
#             "realization": r,
#             "distance": item["distance"],
#             "time": pt["time"],
#             "lat": pt["lat"],
#             "lon": pt["lon"],
#             "wind": pt["wind"],
#             "mslp": pt["mslp"],
#             "geometry": Point(pt["lon"], pt["lat"])
#         })


# gdf_best = gpd.GeoDataFrame(rows, crs="EPSG:4326")
# gdf_best.to_file("best_track_per_realization.gpkg", driver="GPKG")


# # rows_lines = []
# # for r, match in best_match_per_realization.items():
# #     for rank, item in enumerate([match]):

# #         coords = [(pt["lon"], pt["lat"]) for pt in item["track"]]

# #         rows_lines.append({
# #             "realization": r,
# #             "track_rank": rank,
# #             "distance": item["distance"],
# #             "geometry": LineString(coords)
# #         })

# # gdf_lines = gpd.GeoDataFrame(rows_lines, crs="EPSG:4326")
# # gdf_lines.to_file("tracks_best_lines.gpkg", driver="GPKG")



# #%%
# # ---------------------------
# # SAVE RAW TRACKS
# # ---------------------------
# with open(raw_track_file, "wb") as f:
#     pickle.dump(tracks_all_realizations, f)

# # ---------------------------
# # FLATTEN TO TABLE
# # ---------------------------
# rows = []

# for r, data in tracks_all_realizations.items():

#     for is_valid, tracks in zip(
#         [True, False],
#         [data["valid"], data["invalid"]]
#     ):
#         for tid, track in enumerate(tracks):
#             for pt in track:
#                 rows.append({
#                     "realization": r,
#                     "track_id": tid,
#                     "time": pt["time"],
#                     "lat": pt["lat"],
#                     "lon": pt["lon"],
#                     "i": pt["i"],
#                     "j": pt["j"],
#                     "wind": pt["wind"],
#                     "mslp": pt["mslp"],
#                     "is_valid": is_valid
#                 })

# # Save
# df_tracks = pd.DataFrame(rows)
# df_tracks.to_parquet("tracks.parquet")

# # ---------------------------
# # SAVE GEO VERSION
# # ---------------------------
# gdf_tracks = gpd.GeoDataFrame(
#     df_tracks,
#     geometry=gpd.points_from_xy(df_tracks.lon, df_tracks.lat),
#     crs="EPSG:4326"
# )

# gdf_tracks.to_parquet("tracks_geo.parquet")

# print("Tracks saved (raw + parquet + geo)!")




# %%
# starting from Idai
# ---------------------------
# IDAI-GUIDED TRACKING
# One Idai analogue track per realization
# ---------------------------
MSLP_THRESHOLD = 105000  # Pa
SEARCH_RADIUS_KM = 200
raw_track_file = "tracks_all_realizations.pkl"

def find_mslp_minima(mslp_2d):
    return (mslp_2d == minimum_filter(mslp_2d, size=3, mode="nearest"))

def fast_distance_km(lat1, lon1, lat2, lon2):
    return 111 * np.sqrt((lat1 - lat2) ** 2 + (lon1 - lon2) ** 2)

print("Precomputing wind...")
ds_hist["wind"] = np.sqrt(ds_hist["10u"]**2 + ds_hist["10v"]**2)

lat = ds_hist.latitude.values
lon = ds_hist.longitude.values

# ---------------------------
# PRECOMPUTE WIND SEARCH MASK
# ---------------------------
rmw_miles = tc_idai["USA_RMW"].max()
search_radius_wind_km = rmw_miles * 1.60934 * 2

dlat_km = 111 * float(np.mean(np.diff(lat)))
dlon_km = 111 * np.cos(np.deg2rad(np.mean(lat))) * float(np.mean(np.diff(lon)))

window = int(np.ceil(search_radius_wind_km / min(dlat_km, dlon_km)))

ii, jj = np.meshgrid(
    np.arange(-window, window + 1),
    np.arange(-window, window + 1),
    indexing="ij"
)

dist_km = np.sqrt((ii * dlat_km) ** 2 + (jj * dlon_km) ** 2)
dist_mask = dist_km <= search_radius_wind_km

tracks_all_realizations = {}

# ---------------------------
# LOOP OVER REALIZATIONS
# ---------------------------
for r in ds_hist.realization.values:

    print(f"Processing realization {r}...")

    ds_r = ds_hist.sel(realization=r)

    mslp_array = ds_r["msl"]
    wind_array = ds_r["wind"]

    wind_cache = {
        t: wind_array.sel(time=t).values
        for t in ds_r.time.values
    }

    # ---------------------------
    # FIND MINIMA
    # ---------------------------
    candidates_by_time = {}

    times = ds_r.time.values
    mslp_values = mslp_array.values

    for it, t in enumerate(times):

        mslp = mslp_values[it]

        minima_mask = find_mslp_minima(mslp)

        yy, xx = np.where(
            minima_mask &
            (mslp <= MSLP_THRESHOLD)
        )

        candidates = []

        for i, j in zip(yy, xx):
            candidates.append({
                "time": t,
                "lat": lat[i],
                "lon": lon[j],
                "i": i,
                "j": j,
                "mslp": float(mslp[i, j]),
                "realization": int(r)
            })

        candidates_by_time[t] = candidates

    # ---------------------------
    # FOLLOW OBSERVED IDAI TRACK
    # ---------------------------
    idai_track = []

    for obs_pt in tc_idai_points_TC:

        t = np.datetime64(obs_pt["time"])

        if t not in candidates_by_time:
            continue

        best_cand = None
        best_score = None

        for cand in candidates_by_time[t]:
            dist = fast_distance_km(
                obs_pt["lat"],
                obs_pt["lon"],
                cand["lat"],
                cand["lon"]
            )

            if dist > SEARCH_RADIUS_KM:
                continue

            # prioritize proximity, then lower pressure
            score = (dist, cand["mslp"])

            if best_score is None or score < best_score:
                best_score = score
                best_cand = cand

        if best_cand is None:
            continue

        # ---------------------------
        # ADD MAX WIND
        # ---------------------------
        wind_t = wind_cache[best_cand["time"]]

        i0 = best_cand["i"]
        j0 = best_cand["j"]

        i1 = max(0, i0 - window)
        i2 = min(len(lat), i0 + window + 1)

        j1 = max(0, j0 - window)
        j2 = min(len(lon), j0 + window + 1)

        subset = wind_t[i1:i2, j1:j2]

        mask = dist_mask[
            (i1 - i0 + window):(i2 - i0 + window),
            (j1 - j0 + window):(j2 - j0 + window)
        ]

        best_cand["wind"] = float(np.nanmax(subset[mask]))

        idai_track.append(best_cand)

    tracks_all_realizations[int(r)] = idai_track

    print(f"  Analogue track length: {len(idai_track)}")

print("Tracking complete!")

# ---------------------------
# SAVE RAW TRACKS
# ---------------------------
# with open(raw_track_file, "wb") as f:
#     pickle.dump(tracks_all_realizations, f)

# ---------------------------
# SAVE POINTS
# ---------------------------
rows = []

for r, track in tracks_all_realizations.items():

    for pt in track:

        rows.append({
            "realization": r,
            "time": pt["time"],
            "lat": pt["lat"],
            "lon": pt["lon"],
            "i": pt["i"],
            "j": pt["j"],
            "wind": pt["wind"],
            "mslp": pt["mslp"],
            "geometry": Point(pt["lon"], pt["lat"])
        })

gdf_points = gpd.GeoDataFrame(
    rows,
    crs="EPSG:4326"
)

gdf_points.to_file(
    "idai_analogues_points.gpkg",
    driver="GPKG"
)

# ---------------------------
# SAVE LINES
# ---------------------------
rows_lines = []

for r, track in tracks_all_realizations.items():

    if len(track) < 2:
        continue

    rows_lines.append({
        "realization": r,
        "geometry": LineString(
            [(pt["lon"], pt["lat"]) for pt in track]
        )
    })

gdf_lines = gpd.GeoDataFrame(
    rows_lines,
    crs="EPSG:4326"
)

gdf_lines.to_file("idai_analogues_lines.gpkg", driver="GPKG")

# ---------------------------
# PARQUET OUTPUTS
# ---------------------------
# df_tracks = pd.DataFrame(rows)

# df_tracks.to_parquet(
#     "idai_analogues.parquet"
# )

# gdf_points.to_parquet(
#     "idai_analogues_geo.parquet"
# )

print("Saved:")
print("  idai_analogues_points.gpkg")
print("  idai_analogues_lines.gpkg")
# print("  idai_analogues.parquet")
# print("  idai_analogues_geo.parquet")
# print("  tracks_all_realizations.pkl")
# %%

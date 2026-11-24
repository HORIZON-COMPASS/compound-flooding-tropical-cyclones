"""
Add maximum wind speed (within a radius of the track centre) to every point of the
tracks of a chosen tracking strategy, using the pickle master file created by
copilots/run_tracker_strategies.py (keys: (strategy, scenario, realization, storm_name)).

Reuses process_member_wind / add_max_wind_within_radius from tracking_TCs_SI.py.
"""
import os
import pickle
import platform
import argparse
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd
import xarray as xr
import glob
import traceback

from general_TC_funcs import import_ibtracs_period


def save_pickle_atomic(data, output_file):
    output_file = Path(output_file)
    temp_path = output_file.with_name(output_file.name + ".tmp")
    with open(temp_path, "wb") as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(temp_path, output_file)


def add_max_wind_within_radius(track, wind_speed, radius_km=200, time_col="time",
                               lat_col="lat", lon_col="lon"):
    wind_col = "max_wind_model"
    rmw_col = "rmw_model_km"

    track[wind_col] = np.nan
    track[rmw_col] = np.nan
    track["rmw_lat"] = np.nan
    track["rmw_lon"] = np.nan

    r_deg = radius_km / 111.0

    for idx, point in track.iterrows():
        try:
            time = pd.Timestamp(point[time_col])
            lat0 = point[lat_col]
            lon0 = point[lon_col]

            if pd.isna(lat0) or pd.isna(lon0):
                continue

            wind_ts = wind_speed.sel(time=time, method="nearest")

            subset = wind_ts.sel(latitude=slice(lat0 + r_deg, lat0 - r_deg),
                                 longitude=slice(lon0 - r_deg, lon0 + r_deg))
            if subset.size == 0:
                print(f"No grid cells found "
                      f"time={time} "
                      f"lat={lat0:.2f} "
                      f"lon={lon0:.2f}", flush=True)
                continue

            lat2d, lon2d = xr.broadcast(subset.latitude, subset.longitude)

            dx = (lon2d - lon0) * np.cos(np.deg2rad(lat0)) * 111.0  # approximate conversion to km
            dy = (lat2d - lat0) * 111.0  # approximate conversion to km

            distance_km = np.sqrt(dx**2 + dy**2)

            masked_values = subset.where(distance_km <= radius_km).values

            if np.all(np.isnan(masked_values)):
                print(f"All winds NaN "
                      f"time={time} "
                      f"lat={lat0:.2f} "
                      f"lon={lon0:.2f}", flush=True)
                continue

            iy, ix = np.unravel_index(np.nanargmax(masked_values), masked_values.shape)

            track.loc[idx, wind_col] = float(masked_values[iy, ix])
            track.loc[idx, rmw_col] = float(distance_km.values[iy, ix])
            track.loc[idx, "rmw_lat"] = float(subset.latitude.values[iy])
            track.loc[idx, "rmw_lon"] = float(subset.longitude.values[ix])

        except Exception as error:
            print(f"FAILED timestep {idx} "
                  f"time={time} "
                  f"lat={lat0} "
                  f"lon={lon0}: {error}", flush=True)
            raise

    return track


def process_member_wind(args):
    scen, realization, input_file, data_path, output_path, radius_km, use_ibtracs_rmw, storms = args

    output_file = os.path.join(output_path,
                               f"storms_{scen}_r{realization}_wind_{radius_km}km.pkl")

    if os.path.exists(output_file):
        print(f"SKIP {scen} r{realization}: wind file exists", flush=True)
        return output_file

    print(f"START {scen} r{realization}", flush=True)

    with open(input_file, "rb") as f:
        member_tracks = pickle.load(f)

    u_files = sorted(glob.glob(
        os.path.join(data_path, "10u", scen, f"r{realization}",
                     f"climateDT_10u_{scen}_r{realization}_*.nc")))

    v_files = sorted(glob.glob(
        os.path.join(data_path, "10v", scen, f"r{realization}",
                     f"climateDT_10v_{scen}_r{realization}_*.nc")))

    if not u_files or not v_files:
        print(f"Missing wind files: {scen} r{realization}", flush=True)
        return None

    data = xr.open_mfdataset(u_files + v_files, combine="by_coords",
                             data_vars="minimal", coords="minimal",
                             compat="override", parallel=False)

    try:
        wind_speed = np.hypot(data["u10"], data["v10"])

        results = {}
        for key, track in member_tracks.items():
            storm_name = key[2] if isinstance(key, tuple) else str(key)

            if use_ibtracs_rmw and storm_name in storms:
                print(f"Using IBTrACS RMW for wind radius for {storm_name}", flush=True)
                rmw = storms[storm_name]["rmw"].max()

                if not np.isnan(rmw):
                    search_radius = rmw * 1.852 + 20 # from nautical miles to km + buffer of 20 km
                else:
                    print(f"IBTrACS RMW is NaN for {storm_name}, using default radius of {radius_km} km", flush=True)
                    search_radius = radius_km
            else:
                print(f"Using default radius of {radius_km} km for {storm_name}", flush=True)
                search_radius = radius_km

            try:
                results[key] = add_max_wind_within_radius(track=track,
                                                          wind_speed=wind_speed,
                                                          radius_km=search_radius)
                print(f"DONE {scen} r{realization}: {storm_name}", flush=True)

            except Exception as error:
                print(f"FAILED {scen} r{realization} {storm_name}: {error}", flush=True)
                results[key] = track.copy()

    except Exception:
        print(f"REALIZATION FAILED {scen} r{realization}", flush=True)
        traceback.print_exc()
        raise

    finally:
        data.close()

    save_pickle_atomic(results, output_file)

    size_gb = os.path.getsize(output_file) / 1024**3
    print(f"SAVED {scen} r{realization}: {len(results)} storms, {size_gb:.3f} GB", flush=True)

    return output_file




def add_wind_to_strategy_tracks(strategy="all_ibtracs", master_file=None, scenarios=None,
                                realizations=None, wind_radius_km=200, use_ibtracs_rmw=True,
                                start_time="2017-01-01", end_time="2026-08-01",
                                num_workers=1):
    """
    Add max wind speed to all tracked points of one tracking strategy.

    Parameters
    ----------
    strategy : str
        Tracking strategy in the master file (e.g. "all_ibtracs", "2point", "5point").
    master_file : path-like, optional
        Master pickle file. Defaults to storms_all_strategies_{start}_{end}.pkl.
    scenarios, realizations : list of str, optional
        Subset of members to process. Default: all present for the strategy.
    wind_radius_km : float
        Default search radius; used when IBTrACS RMW is unavailable or not requested.
    use_ibtracs_rmw : bool
        Use max IBTrACS RMW (nm -> km, + 20 km) as search radius per storm.

    Returns
    -------
    dict
        {(strategy, scenario, realization, storm_name): track_df with wind columns}
    """
    if platform.system() == "Windows":
        data_path_base = Path("c:/Code/test/ClimateDT/")
        processed_path = Path("P:/11210471-001-compass/01_Data/ECMWF_ClimateDT/analysis_output/")
    else:
        data_path_base = Path("/projects/prjs2226/data/ClimateDT/sfc/raw/")
        processed_path = Path("/projects/prjs2226/data/ClimateDT/sfc/processed/")

    tracks_path = processed_path / "msl" / "tracked"
    master_file = Path(master_file) if master_file else (
        tracks_path / f"storms_all_strategies_{start_time}_{end_time}.pkl")

    rmw_tag = "rmw" if use_ibtracs_rmw else f"{wind_radius_km}km"
    out_path = processed_path / "wind" / f"tracked_{strategy}_{rmw_tag}"
    member_path = out_path / "members"
    os.makedirs(member_path, exist_ok=True)

    print(f"Loading master file: {master_file}", flush=True)
    with open(master_file, "rb") as f:
        master = pickle.load(f)

    # Group the selected strategy's tracks per (scenario, realization)
    members = {}
    for (strat, scen, real, storm), track in master.items():
        if strat != strategy or track is None:
            continue
        if scenarios and scen not in scenarios:
            continue
        if realizations and str(real) not in [str(r) for r in realizations]:
            continue
        members.setdefault((scen, str(real)), {})[storm] = track
    del master

    if not members:
        raise ValueError(f"No tracks found for strategy '{strategy}' in {master_file}")

    storms = import_ibtracs_period(start=start_time, end=end_time) if use_ibtracs_rmw else {}

    # process_member_wind reads per-member pickles keyed by storm name
    jobs = []
    for (scen, real), tracks in members.items():
        input_file = member_path / f"storms_{scen}_r{real}_{strategy}.pkl"
        save_pickle_atomic(tracks, input_file)
        jobs.append((scen, real, str(input_file), data_path_base, str(out_path),
                     wind_radius_km, use_ibtracs_rmw, storms))

    with Pool(processes=num_workers) as pool:
        wind_files = [r for r in pool.imap_unordered(process_member_wind, jobs) if r is not None]

    results = {}
    for wind_file in sorted(wind_files):
        scen, real = os.path.basename(wind_file).split("_wind_")[0][len("storms_"):].rsplit("_r", 1)
        with open(wind_file, "rb") as f:
            for storm, track in pickle.load(f).items():
                results[(strategy, scen, real, storm)] = track

    out_file = out_path / f"storms_{strategy}_max_wind_{rmw_tag}_{start_time}_{end_time}.pkl"
    save_pickle_atomic(results, out_file)
    print(f"Saved {len(results)} tracks to {out_file}", flush=True)

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--strategy", default="all_ibtracs")
    parser.add_argument("--master-file", default=None)
    parser.add_argument("--scen", nargs="*", default=None)
    parser.add_argument("--realization", nargs="*", default=None)
    parser.add_argument("--radius-km", type=float, default=200)
    parser.add_argument("--no-ibtracs-rmw", action="store_true")
    parser.add_argument("--workers", type=int, default=1)
    a = parser.parse_args()

    add_wind_to_strategy_tracks(strategy=a.strategy, master_file=a.master_file,
                                scenarios=a.scen, realizations=a.realization,
                                wind_radius_km=a.radius_km,
                                use_ibtracs_rmw=a.no_ibtracs_rmw,
                                num_workers=a.workers)

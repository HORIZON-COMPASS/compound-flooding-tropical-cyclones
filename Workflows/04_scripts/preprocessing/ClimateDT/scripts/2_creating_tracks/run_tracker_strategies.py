"""
Parallel execution wrapper for different storm tracking strategies.

Allows running multiple tracker strategies on the same data for comparison.
"""

import os
import pickle
import pandas as pd
import xarray as xr
from pathlib import Path
from multiprocessing import Pool
import platform
import glob

from tracker_strategies import (
    strategy_2point_init,
    strategy_5point_init,
    strategy_all_ibtracs_init,
)


STRATEGY_FUNCTIONS = {
    "2point": strategy_2point_init,
    "5point": strategy_5point_init,
    "all_ibtracs": strategy_all_ibtracs_init,
}


def strategy_file(processed_path_tracks, scen, realization, strategy_name):
    """Path of the result file for one scenario/realization/strategy."""
    return os.path.join(processed_path_tracks,
                        f"storms_{scen}_r{realization}_strategy_{strategy_name}.pkl")


def parse_strategy_file(path):
    """Return (scen, realization, strategy) from a per-strategy result filename."""
    name = os.path.basename(path)[len("storms_"):-len(".pkl")]
    head, strategy = name.split("_strategy_", 1)
    scen, realization = head.rsplit("_r", 1)
    return scen, realization, strategy


def process_member_tracking_strategies(args):
    """
    Track a single ensemble member using multiple strategies.
    
    Parameters
    ----------
    args : tuple
        (scen, realization, storms, data_path_msl, processed_path_tracks)
    
    Returns
    -------
    list
        Paths of result files written by this call
    """
    scen, realization, storms, data_path_msl, processed_path_tracks = args
    
    print(f"Processing {scen} r{realization}", flush=True)

    pending = [name for name in STRATEGY_FUNCTIONS
               if not os.path.exists(strategy_file(processed_path_tracks, scen, realization, name))]
    if not pending:
        print(f"SKIP {scen} r{realization}: all strategies done", flush=True)
        return []

    data_path = data_path_msl / scen / f"r{realization}"

    if not data_path.exists():
        print(f"Path does not exist: {data_path}", flush=True)
        return None

    files = sorted(data_path.glob("*.nc"))
    if not files:
        print(f"No NetCDF files found in {data_path}", flush=True)
        return None

    data = xr.open_mfdataset(files, combine="by_coords", parallel=False)

    # Rename coordinates if necessary (latitude -> lat, longitude -> lon)
    if "latitude" in data.coords:
        data = data.rename({"latitude": "lat", "longitude": "lon"})

    saved = []
    try:
        for strategy_name in pending:
            strategy_func = STRATEGY_FUNCTIONS[strategy_name]
            tracked_storms = {}
            for storm_name, best_track in storms.items():
                try:
                    tracked_storms[storm_name] = strategy_func(data, best_track, s_size=1.5)
                    print(f"{storm_name} tracked with {strategy_name}", flush=True)
                except Exception as e:
                    print(f"Failed to track {storm_name} with {strategy_name}: {e}", flush=True)
                    tracked_storms[storm_name] = None

            n_ok = len([s for s in tracked_storms.values() if s is not None])
            print(f"  {strategy_name}: completed {n_ok}/{len(storms)} storms", flush=True)

            # Save each strategy immediately (write-then-rename) so work survives a job timeout
            output_file = strategy_file(processed_path_tracks, scen, realization, strategy_name)
            tmp_file = output_file + ".tmp"
            with open(tmp_file, "wb") as f:
                pickle.dump(tracked_storms, f, protocol=pickle.HIGHEST_PROTOCOL)
            os.replace(tmp_file, output_file)
            saved.append(output_file)
            print(f"Saved {os.path.basename(output_file)}", flush=True)
    finally:
        data.close()

    return saved


def run_strategy_comparison(experiments, realizations, storms, data_path_msl, 
                            processed_path_tracks, num_workers=5):
    """
    Run all tracking strategies on specified experiments and realizations.
    
    Parameters
    ----------
    experiments : list
        List of scenario names (e.g., ['hist', 'cont', 'Tplus2.0K'])
    realizations : list
        List of realization numbers (e.g., ['1', '2', '3', '4', '5'])
    storms : dict
        Dictionary of storms from IBTrACS
    data_path_msl : Path
        Base path to MSL data
    processed_path_tracks : Path
        Output directory for results
    num_workers : int
        Number of parallel workers
    """
    jobs = []
    for scen in experiments:
        for real_idx, realization in enumerate(realizations):
            jobs.append((scen, realization, storms, data_path_msl, processed_path_tracks))
    
    print(f"Submitting {len(jobs)} jobs for strategy comparison", flush=True)
    
    with Pool(processes=num_workers) as pool:
        for result in pool.imap_unordered(process_member_tracking_strategies, jobs):
            if result:
                print(f"Completed: {[os.path.basename(r) for r in result]}", flush=True)



def combine_strategy_results(processed_path_tracks, start_time, end_time):
    """
    Combine all individual strategy pickle files into a single master file.
    
    Creates a flattened dictionary with tuple keys:
    {(strategy, scenario, realization, storm_name): track_df}
    
    This format is directly compatible with compare_climatedt_tracks()
    and plot_tracks_all() functions.
    
    Parameters
    ----------
    processed_path_tracks : Path or str
        Path to directory containing individual result files
    start_time : str
        Start time for naming the combined file (e.g., "2017-01-01")
    end_time : str
        End time for naming the combined file (e.g., "2026-08-01")
    
    Returns
    -------
    str
        Path to the combined master pickle file
    """
    processed_path_tracks = Path(processed_path_tracks)
    
    # Master dictionary with tuple keys for efficient access
    master_results = {}
    
    # Load all individual pickle files
    pattern = str(processed_path_tracks / "storms_*_r*_strategy_*.pkl")
    pkl_files = sorted(glob.glob(pattern))

    if not pkl_files:
        print("No strategy result files found to combine", flush=True)
        return None

    print(f"Combining {len(pkl_files)} strategy result files...", flush=True)
    print(f"Creating flattened tuple-key structure...", flush=True)
    
    total_storms = 0

    for pkl_file in pkl_files:
        scenario, realization, strategy = parse_strategy_file(pkl_file)
        print(f"  Loading {scenario} r{realization} {strategy}...", flush=True)
        
        with open(pkl_file, "rb") as f:
            strategy_data = pickle.load(f)
        
        # Flatten the strategy data into tuple keys
        # Input: {storm_name: track_df}
        # Output: {(strategy, scenario, realization, storm_name): track_df}
        for storm_name, track_df in strategy_data.items():
            key = (strategy, scenario, realization, storm_name)
            master_results[key] = track_df
            total_storms += 1

    # Save combined result to master file
    master_file = processed_path_tracks / f"storms_all_strategies_{start_time}_{end_time}.pkl"
    
    with open(master_file, "wb") as f:
        pickle.dump(master_results, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    print(f"\n{'='*80}")
    print(f"Successfully combined {len(pkl_files)} strategy result files")
    print(f"Total storms: {total_storms}")
    print(f"Unique strategies: {len(set(k[0] for k in master_results.keys()))}")
    print(f"Unique scenarios: {len(set(k[1] for k in master_results.keys()))}")
    print(f"Unique realizations: {len(set(k[2] for k in master_results.keys()))}")
    print(f"\nSaved to: {master_file}")
    print(f"{'='*80}\n")
    
    return str(master_file)


def load_strategy_results(processed_path_tracks, scen=None, realization=None):
    """
    Load tracking results from pickle files.
    
    Parameters
    ----------
    processed_path_tracks : Path
        Path to directory containing result files
    scen : str, optional
        Filter by scenario. If None, loads all.
    realization : str, optional
        Filter by realization. If None, loads all.
    
    Returns
    -------
    dict
        Nested dictionary: {scen: {realization: {strategy: tracked_storms}}}
    """
    results = {}
    
    import glob

    pattern = f"storms_{scen or '*'}_r{realization or '*'}_strategy_*.pkl"
    files = sorted(glob.glob(os.path.join(processed_path_tracks, pattern)))

    for pkl_file in files:
        file_scen, file_real, strategy = parse_strategy_file(pkl_file)
        with open(pkl_file, "rb") as f:
            results.setdefault(file_scen, {}).setdefault(file_real, {})[strategy] = pickle.load(f)

    return results


def summary_stats_by_strategy(strategy_results):
    """
    Generate summary statistics for each strategy.
    
    Parameters
    ----------
    strategy_results : dict
        Nested dictionary of tracking results by strategy
    
    Returns
    -------
    pd.DataFrame
        Summary statistics with rows for each strategy
    """
    summary = {}
    
    for strategy_name, tracked_storms in strategy_results.items():
        if tracked_storms is None:
            summary[strategy_name] = {'status': 'FAILED',
                                      'n_storms_attempted': 0,
                                      'n_storms_successful': 0,
                                      'success_rate': 0.0}
            continue
        
        n_attempted = len(tracked_storms)
        n_successful = sum(1 for v in tracked_storms.values() if v is not None)
        
        summary[strategy_name] = {'status': 'SUCCESS',
                                  'n_storms_attempted': n_attempted,
                                  'n_storms_successful': n_successful,
                                  'success_rate': n_successful / n_attempted if n_attempted > 0 else 0.0}
    
    return pd.DataFrame(summary).T


if __name__ == "__main__":
    import sys
    from general_TC_funcs import import_ibtracs_period
    
    if platform.system() == "Windows":
        data_path_base = Path("c:/Code/test/ClimateDT/")
        processed_path = Path("P:/11210471-001-compass/01_Data/ECMWF_ClimateDT/analysis_output/")
    else:
        data_path_base = Path("/projects/prjs2226/data/ClimateDT/sfc/raw/")
        processed_path = Path("/projects/prjs2226/data/ClimateDT/sfc/processed/")
    
    processed_path_tracks = Path(os.path.join(processed_path, "msl", "tracked"))
    
    os.makedirs(processed_path_tracks, exist_ok=True)
    
    start_time = "2017-01-01"
    end_time = "2026-08-01"
    
    # Load IBTrACS reference
    print("Loading IBTrACS data...", flush=True)
    storms = import_ibtracs_period(start=start_time, end=end_time)
    
    # Run strategy comparison
    experiments = [
        "hist", 
        "cont", 
        "Tplus2.0K"
        ]  # Adjust as needed
    realizations = ["1", "2", "3", "4", "5"]  # Adjust as needed
    
    run_strategy_comparison(experiments, 
                            realizations, 
                            storms, 
                            data_path_base / "msl",  
                            processed_path_tracks,
                            num_workers=5)
    
    print("Strategy comparison complete!", flush=True)
    
    # Combine all individual results into a master file
    print("\nCombining strategy results into master file...", flush=True)
    master_file = combine_strategy_results(processed_path_tracks, start_time, end_time)
    
    if master_file:
        print(f"Master file created: {master_file}", flush=True)


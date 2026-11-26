
import os
import pickle
from pathlib import Path
import platform
import numpy as np
import pandas as pd
import xarray as xr
import warnings
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from general_TC_funcs import import_ibtracs_period


def haversine_distance(lat1, lon1, lat2, lon2):
    """
    Calculate haversine distance between two points in km.
    
    Parameters
    ----------
    lat1, lon1 : float or array
        First point(s)
    lat2, lon2 : float or array
        Second point(s)
    
    Returns
    -------
    float or array
        Distance in km
    """
    R = 6371  # Earth radius in km
    
    lat1_rad = np.radians(lat1)
    lat2_rad = np.radians(lat2)
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    
    a = np.sin(dlat/2)**2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    
    return R * c


def normalize_pressure_column(df, pres_col_name):
    """
    Normalize pressure column to hPa units and rename to 'msl'.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with pressure column
    pres_col_name : str
        Name of the pressure column ('msl', 'min_pres', or 'pressure')
    
    Returns
    -------
    pd.DataFrame
        DataFrame with pressure column renamed to 'msl' and converted to hPa if needed
    """
    df = df.copy()

    # Detect units and convert to hPa if necessary
    # If max value > 10000, likely in Pa (100000 Pa = 1000 hPa typical for storms)
    pressure_values = df[pres_col_name].dropna()
    
    if len(pressure_values) > 0:
        max_pres = pressure_values.max()
        if max_pres > 10000:  # Likely in Pa
            df[pres_col_name] = df[pres_col_name] / 100.0
    
    # Rename to standard 'msl' name
    if pres_col_name != 'msl':
        df = df.rename(columns={pres_col_name: 'msl'})
    
    return df



# load and align tracks
def load_existing_tracks(track_file):
    print(f"Loading existing tracked results: {track_file}", flush=True)
    with open(track_file, "rb") as f:
        tracks = pickle.load(f)
    
    return tracks


def align_tracks(tracked_track, reference_track, time_col='time'):
    """
    Align two tracks by time and return matched point pairs.
    
    Normalizes pressure units to hPa and standardizes column names.
    Filters reference data to only include complete observations (lat, lon, msl all present).
    Tracked data is assumed to already have complete observations.
    
    Parameters
    ----------
    tracked_track : pd.DataFrame
        Tracked storm with columns: time, lat, lon, msl/min_pres (all complete)
    reference_track : pd.DataFrame
        Reference track (IBTrACS) with same columns
    time_col : str
        Name of time column
    
    Returns
    -------
    pd.DataFrame
        Aligned tracks with tracked and reference position and pressure columns
        next to each other. Model step and meanSea coordinates are omitted.
        Only rows where reference has all three fields (lat, lon, msl) present.
        All pressures in hPa.
    """
    if tracked_track is None or len(tracked_track) == 0:
        return None
    
    # Normalize pressure columns in both datasets
    tracked_track = tracked_track.copy()
    reference_track = reference_track.copy()
    
    # Find and normalize reference pressure column
    tracked_pres_col = 'msl'
    tracked_track = normalize_pressure_column(tracked_track, tracked_pres_col)
  
    ref_pres_col = 'min_pres'
    reference_track = normalize_pressure_column(reference_track, ref_pres_col)

    # Filter reference to only include complete observations (lat, lon, msl all present)
    reference_complete = reference_track[reference_track[['lat', 'lon', 'msl']].notna().all(axis=1)].copy()
    
    # Ensure time columns are datetime (reference already floored to hourly in import_ibtracs_period)
    tracked_time = pd.to_datetime(tracked_track[time_col])
    reference_time = pd.to_datetime(reference_complete[time_col])
    
    # Find overlap
    time_min = max(tracked_time.min(), reference_time.min())
    time_max = min(tracked_time.max(), reference_time.max())
    
    if time_min > time_max or pd.isna(time_min) or pd.isna(time_max):
        return None
    
    # Merge on closest times
    tracked_subset = tracked_track[(tracked_time >= time_min) & (tracked_time <= time_max)].copy()
    reference_subset = reference_complete.copy()
    
    if len(tracked_subset) == 0 or len(reference_subset) == 0:
        return None
    
    # Rename columns before merge to ensure suffixes are applied properly
    tracked_cols = {c: f'{c}_tracked' if c != time_col else c 
                    for c in tracked_subset.columns}
    ref_cols = {c: f'{c}_ref' if c != time_col else c 
                for c in reference_subset.columns}
    
    tracked_subset = tracked_subset.rename(columns=tracked_cols)
    reference_subset = reference_subset.rename(columns=ref_cols)
    
    aligned = pd.merge_asof(tracked_subset.sort_values(time_col),
                            reference_subset.sort_values(time_col),
                            on=time_col, direction='nearest')

    aligned = aligned.drop(columns=['step_tracked', 'meanSea_tracked'], errors='ignore')
    paired_cols = [time_col, 'valid_time_tracked', 'lat_tracked', 'lat_ref', 'lon_tracked', 'lon_ref',
                   'msl_tracked', 'msl_ref']
    other_cols = [column for column in aligned.columns if column not in paired_cols]
    aligned = aligned[paired_cols + other_cols]

    return aligned if not aligned.empty else None


def create_aligned_tracks_dataframe(storms, tracks, tracking_strategy=None):
    """Build an unfiltered table of matched points for inspection.

    ``should_exclude`` marks matches whose IBTrACS latitude is south of 40S.
    The tracked and reference timestamps are stored separately.
    """
    aligned_tracks = []
    for key, track in tracks.items():
        if track is None:
            continue
        if len(key) == 4:
            strategy, scenario, realization, name = key
        elif len(key) == 3:
            scenario, realization, name = key
            strategy = tracking_strategy if tracking_strategy is not None else "2-point-init"
        else:
            print(f"Unexpected key format (len={len(key)}): {key}, skipping")
            continue

        reference = storms[name].copy()
        reference["reference_time"] = reference["time"]
        aligned = align_tracks(track, reference)
        if aligned is None:
            print(f"No aligned points for {key}, skipping")
            continue

        aligned["strategy"] = strategy
        aligned["scenario"] = scenario
        aligned["realization"] = realization
        aligned["storm"] = name
        aligned["should_exclude"] = aligned["lat_ref"] < -40
        aligned_tracks.append(aligned)

    if not aligned_tracks:
        print("No aligned tracks found")
        return pd.DataFrame()
    return pd.concat(aligned_tracks, ignore_index=True)


def add_aligned_wind(aligned_tracks, wind_tracks, label):
    """Add per-point model wind (m/s) using the exact track key and time.

    The wind tracks must contain unique times and matching track coordinates.
    Rows without a matching wind point retain NaN and are reported.
    """
    keys = ["strategy", "scenario", "realization", "storm", "time"]
    wind_column = f"max_wind_model_{label}"
    if wind_column in aligned_tracks.columns:
        raise ValueError(f"Column already exists: {wind_column}")

    wind_frames = []
    for key, track in wind_tracks.items():
        if len(key) != 4:
            raise ValueError(f"Expected a four-part wind track key, got {key!r}")
        if track is None:
            print(f"No wind track for {key}")
            continue
        strategy, scenario, realization, storm = key
        wind = track[["time", "lat", "lon", "max_wind_model"]].copy()
        wind["strategy"] = strategy
        wind["scenario"] = scenario
        wind["realization"] = realization
        wind["storm"] = storm
        wind_frames.append(wind)

    if not wind_frames:
        raise ValueError(f"No wind tracks found for {label}")

    wind_data = pd.concat(wind_frames, ignore_index=True)
    if wind_data[keys].isna().any().any():
        raise ValueError(f"Missing track keys or times in {label} wind data")
    if wind_data.duplicated(keys).any():
        raise ValueError(f"Duplicate track keys and times in {label} wind data")

    wind_data = wind_data.rename(columns={
        "lat": "_wind_lat",
        "lon": "_wind_lon",
        "max_wind_model": wind_column,
    })
    merged = aligned_tracks.merge(wind_data, on=keys, how="left",
                                  validate="many_to_one", indicator="_wind_match")
    matched = merged["_wind_match"].eq("both")
    same_position = (
        np.isclose(merged.loc[matched, "lat_tracked"], merged.loc[matched, "_wind_lat"],
                   rtol=0, atol=1e-5)
        & np.isclose(merged.loc[matched, "lon_tracked"], merged.loc[matched, "_wind_lon"],
                     rtol=0, atol=1e-5)
    )
    if not same_position.all():
        raise ValueError(f"{label} wind points do not match tracked coordinates")

    missing_count = (~matched).sum()
    nan_count = merged.loc[matched, wind_column].isna().sum()
    print(f"{label} wind: {matched.sum()} matched points, {missing_count} unmatched, "
          f"{nan_count} matched with missing wind")
    return merged.drop(columns=["_wind_lat", "_wind_lon", "_wind_match"])



# Calculate and plot track comparison (metrics)
def compute_position_metrics(tracked_track, reference_track, aligned=None):
    """
    Compute position error metrics.
    
    Parameters
    ----------
    tracked_track : pd.DataFrame
        Tracked positions
    reference_track : pd.DataFrame
        Reference positions
    aligned : pd.DataFrame, optional
        Already aligned point pairs to use instead of aligning again.
    
    Returns
    -------
    dict
        Dictionary with metrics: MDE, CEP, MXE, etc.
    """
    if aligned is None:
        aligned = align_tracks(tracked_track, reference_track)
    if aligned is None:
        raise ValueError("No aligned track points available for position metrics")
    
    # Calculate distances
    distances = haversine_distance(aligned['lat_tracked'].values, aligned['lon_tracked'].values,
                                   aligned['lat_ref'].values, aligned['lon_ref'].values)
    
    # Compute metrics
    metrics = {
        'n_points': len(aligned),
        'MDE': np.mean(distances),  # Mean Distance Error
        'CEP': np.percentile(distances, 50),  # Circular Error Probable (50th percentile)
        'MXE': np.max(distances),  # Maximum Error
        'P95E': np.percentile(distances, 95),  # 95th percentile
        'P90E': np.percentile(distances, 90),  # 90th percentile
        'std_error': np.std(distances),
        'RMSE': np.sqrt(np.mean(distances**2)),
    }
    
    return metrics


def compute_intensity_metrics(tracked_track, reference_track, storm_name=None, debug=False, aligned=None):
    """
    Compute intensity (MSLP) error metrics.
    
    Parameters
    ----------
    tracked_track : pd.DataFrame
        Tracked pressures
    reference_track : pd.DataFrame
        Reference pressures
    storm_name : str, optional
        Storm name for debugging output
    debug : bool, optional
        Print debug information
    aligned : pd.DataFrame, optional
        Already aligned point pairs to use instead of aligning again.
    
    Returns
    -------
    dict
        Dictionary with metrics: RMSE, bias, correlation, etc.
    """
    if aligned is None:
        aligned = align_tracks(tracked_track, reference_track)
    if aligned is None:
        raise ValueError("No aligned track points available for intensity metrics")
    
    # Get pressure columns (could be named msl or different)
    pres_cols = [c for c in aligned.columns if 'msl' in c.lower() or 'min_pres' in c.lower()]
    
    if debug:
        print(f"[{storm_name}] Pressure columns found: {pres_cols}")
    
    # Find tracked and reference pressure columns
    pres_tracked_cols = [c for c in pres_cols if 'tracked' in c.lower()]
    pres_ref_cols = [c for c in pres_cols if 'ref' in c.lower()]

    pres_tracked = aligned[pres_tracked_cols[0]].values
    pres_ref = aligned[pres_ref_cols[0]].values
    
    # Compute metrics
    errors = pres_tracked - pres_ref
    metrics = {
        'n_points': len(aligned),
        'RMSE': np.sqrt(np.mean(errors**2)),
        'bias': np.mean(errors),
        'std_error': np.std(errors),
        'correlation': np.corrcoef(pres_tracked, pres_ref)[0, 1] if len(pres_tracked) > 1 else np.nan,
    }
    
    return metrics


def compute_track_progression_metrics(tracked_track, reference_track):
    """
    Compute metrics related to track progression and longevity.
    
    Normalizes pressure units to hPa and standardizes column names.
    Filters reference data to only include complete observations (lat, lon, msl).
    Tracked data is assumed to already have complete observations.
    
    Parameters
    ----------
    tracked_track : pd.DataFrame
        Tracked positions (all complete)
    reference_track : pd.DataFrame
        Reference positions
    
    Returns
    -------
    dict
        Dictionary with metrics: duration, completion rate, etc.
        
    Notes
    -----
    completion_rate : ratio of tracked duration to reference duration
        - 1.0 means the tracked storm spans the entire reference period
        - < 1.0 means the tracked storm was shorter (earlier start or earlier end)
    point_coverage : ratio of tracked points to reference points
        - Indicates data frequency/resolution match between tracked and reference
        - Only counts complete observations in reference
    """
    if tracked_track is None or len(tracked_track) == 0:
        reference_track = reference_track.copy()
        
        # Find and normalize reference pressure column
        ref_pres_col = None
        if 'msl' in reference_track.columns:
            ref_pres_col = 'msl'
        elif 'min_pres' in reference_track.columns:
            ref_pres_col = 'min_pres'
        elif 'pressure' in reference_track.columns:
            ref_pres_col = 'pressure'
        
        if ref_pres_col:
            reference_track = normalize_pressure_column(reference_track, ref_pres_col)
        
        # Filter reference to complete observations only
        reference_complete = reference_track[
            reference_track[['lat', 'lon', 'msl']].notna().all(axis=1)
        ]
        
        reference_duration = (reference_complete['time'].max() - reference_complete['time'].min()).total_seconds() / 3600 if len(reference_complete) > 1 else 0
        
        return {
            'track_duration_hours': 0,
            'reference_duration_hours': reference_duration,
            'n_tracked_points': 0,
            'n_reference_points': len(reference_complete),
            'completion_rate': 0.0,
            'point_coverage': 0.0,
        }
    
    # Normalize pressure in both datasets
    tracked_track = tracked_track.copy()
    reference_track = reference_track.copy()
    
    # Find and normalize tracked pressure column
    tracked_pres_col = None
    if 'msl' in tracked_track.columns:
        tracked_pres_col = 'msl'
    elif 'min_pres' in tracked_track.columns:
        tracked_pres_col = 'min_pres'
    elif 'pressure' in tracked_track.columns:
        tracked_pres_col = 'pressure'
    
    if tracked_pres_col:
        tracked_track = normalize_pressure_column(tracked_track, tracked_pres_col)
    
    # Find and normalize reference pressure column
    ref_pres_col = None
    if 'msl' in reference_track.columns:
        ref_pres_col = 'msl'
    elif 'min_pres' in reference_track.columns:
        ref_pres_col = 'min_pres'
    elif 'pressure' in reference_track.columns:
        ref_pres_col = 'pressure'
    
    if ref_pres_col:
        reference_track = normalize_pressure_column(reference_track, ref_pres_col)
    
    # Filter reference to complete observations only (lat, lon, msl all present)
    reference_complete = reference_track[
        reference_track[['lat', 'lon', 'msl']].notna().all(axis=1)
    ]
    
    if len(reference_complete) == 0:
        return {
            'track_duration_hours': 0,
            'reference_duration_hours': 0,
            'n_tracked_points': len(tracked_track),
            'n_reference_points': 0,
            'completion_rate': 0.0,
            'point_coverage': 0.0,
        }
    
    track_duration = (tracked_track['time'].max() - tracked_track['time'].min()).total_seconds() / 3600
    reference_duration = (reference_complete['time'].max() - reference_complete['time'].min()).total_seconds() / 3600
    
    # Completion rate: ratio of duration tracked vs reference
    # 1.0 = tracked the entire storm, < 1.0 = tracked only part of the storm
    completion_rate = track_duration / reference_duration if reference_duration > 0 else 0.0
    
    # Point coverage: ratio of point density
    # Indicates if temporal resolution is similar
    point_coverage = len(tracked_track) / len(reference_complete) if len(reference_complete) > 0 else 0.0
    
    return {
        'track_duration_hours': track_duration,
        'reference_duration_hours': reference_duration,
        'n_tracked_points': len(tracked_track),
        'n_reference_points': len(reference_complete),
        'completion_rate': completion_rate,
        'point_coverage': point_coverage,
    }


def compare_all_strategies(storm_name, strategies_tracks, reference_track, debug=False,
                           min_ref_lat=None, min_ref_sshs=None, filter_metrics=False,
                           land_points=None):
    """
    Compare all strategies for a single storm.

    Duration warnings are printed if debug=True. A warning is emitted when a
    strategy has no aligned points after filtering.

    Parameters
    ----------
    storm_name : str
        Name of the storm
    strategies_tracks : dict
        Dictionary mapping strategy name to tracked DataFrame
    reference_track : pd.DataFrame
        Reference track from IBTrACS
    debug : bool, optional
        Print debug information and warnings
    min_ref_lat : float, optional
        Minimum latitude of matched IBTrACS points, inclusive.
    min_ref_sshs : int, optional
        Minimum usa_sshs category of matched IBTrACS points, inclusive.
        When filters are enabled, track duration spans the first and last
        retained matched points, including any gaps between them.
    filter_metrics : bool, optional
        If True, apply latitude and SSHS filters to ALL metrics (position, intensity, duration).
        If False, only apply to duration metrics (default behavior).
    land_points : str, optional
        How to handle landfall points: 'excl' excludes reference points with landfall_ref=0,
        None includes all points.

    Returns
    -------
    pd.DataFrame
        Comparison table with all metrics for all strategies
    """
    results = []
    warning_messages = []
    
    # Filter reference track for landfall if specified
    reference_to_use = reference_track.copy()
        
    for strategy_name, tracked_track in strategies_tracks.items():
        if tracked_track is None:
            continue
        
        aligned = align_tracks(tracked_track, reference_to_use)
        if aligned is None:
            warnings.warn(f"{storm_name} ({strategy_name}): no aligned points after filtering; skipping",
                          stacklevel=2)
            continue

        # Apply filters to ALL metrics if filter_metrics=True
        if filter_metrics and (min_ref_lat is not None or min_ref_sshs is not None or land_points == 'excl'):
            # Filter aligned points by criteria
            filtered_aligned = aligned.copy()
            if min_ref_lat is not None:
                filtered_aligned = filtered_aligned[filtered_aligned['lat_ref'] >= min_ref_lat]
            if min_ref_sshs is not None:
                filtered_aligned = filtered_aligned[filtered_aligned['usa_sshs_ref'] >= min_ref_sshs]
            if land_points == 'excl':
                if 'landfall_ref' in filtered_aligned.columns:
                    filtered_aligned = filtered_aligned[filtered_aligned['landfall_ref'] == 0].copy()
                else:
                    print(f"Warning: 'landfall_ref' column not found, skipping landfall filter")
            
            # Skip if no points remain after filtering
            if len(filtered_aligned) == 0:
                warning_messages.append(f"{storm_name} ({strategy_name}): no points after filtering; skipping")
                continue
            
            # Compute metrics on filtered aligned points
            pos_metrics = compute_position_metrics(tracked_track, reference_to_use, aligned=filtered_aligned)
            int_metrics = compute_intensity_metrics(tracked_track, reference_to_use, storm_name=storm_name,
                                                    debug=debug, aligned=filtered_aligned)
            progression_track = tracked_track[tracked_track['time'].isin(filtered_aligned['time'])]
            progression_reference = reference_to_use
            if min_ref_lat is not None:
                progression_reference = progression_reference[progression_reference['lat'] >= min_ref_lat]
            if min_ref_sshs is not None:
                progression_reference = progression_reference[progression_reference['usa_sshs'] >= min_ref_sshs]
        else:
            # Original behavior: only filter duration metrics
            pos_metrics = compute_position_metrics(tracked_track, reference_to_use, aligned=aligned)
            int_metrics = compute_intensity_metrics(tracked_track, reference_to_use, storm_name=storm_name,
                                                    debug=debug, aligned=aligned)
            if min_ref_lat is not None or min_ref_sshs is not None:
                progression_track = tracked_track[tracked_track['time'].isin(aligned['time'])]
                progression_reference = reference_to_use
                if min_ref_lat is not None:
                    progression_reference = progression_reference[progression_reference['lat'] >= min_ref_lat]
                if min_ref_sshs is not None:
                    progression_reference = progression_reference[progression_reference['usa_sshs'] >= min_ref_sshs]
            else:
                progression_track = tracked_track
                progression_reference = reference_to_use
        
        prog_metrics = compute_track_progression_metrics(progression_track, progression_reference)
        
        # Check for warnings       
        if prog_metrics['completion_rate'] < 1.0:
            warning_messages.append(f"{storm_name}: Incomplete tracking ({prog_metrics['completion_rate']:.1%} of storm duration)")
        
        if prog_metrics['track_duration_hours'] != prog_metrics['reference_duration_hours']:
            warning_messages.append(f"{storm_name}: Duration mismatch (tracked={prog_metrics['track_duration_hours']:.1f}h, ref={prog_metrics['reference_duration_hours']:.1f}h)")
        
        result = {
            'storm': storm_name,
            'strategy': strategy_name,
            'n_points_compared': pos_metrics['n_points'],
            'MDE_km': pos_metrics['MDE'],
            'CEP_km': pos_metrics['CEP'],
            'MXE_km': pos_metrics['MXE'],
            'P90E_km': pos_metrics['P90E'],
            'P95E_km': pos_metrics['P95E'],
            'position_RMSE_km': pos_metrics['RMSE'],
            'RMSE_pres_hPa': int_metrics['RMSE'],
            'bias_pres_hPa': int_metrics['bias'],
            'corr_pres': int_metrics['correlation'],
            'track_duration_hrs': prog_metrics['track_duration_hours'],
            'n_tracked_points': prog_metrics['n_tracked_points'],
            'n_reference_points': prog_metrics['n_reference_points'],
            'data_coverage': prog_metrics['point_coverage'],
        }
        
        results.append(result)
    
    df = pd.DataFrame(results)
    
    # Print warnings if debug enabled
    if debug and warning_messages:
        print("\nWarnings:")
        for w in warning_messages:
            print(f"{w}")
    
    return df


def compare_climatedt_tracks(storms, tracks, tracking_strategy=None, results_dir=None,
                             plot_mode='none', min_ref_lat=-40, min_ref_sshs=None,
                             scen=None, land_points=None):
    """
    Compare tracked storms against IBTrACS reference.

    Parameters
    ----------
    storms : dict
        Dictionary of reference storm tracks (IBTrACS)
    tracks : dict
        Dictionary of tracked storm data with keys:
        - (scenario, realization, storm_name) for single strategy
        - (strategy, scenario, realization, storm_name) for multiple strategies
    tracking_strategy : str, optional
        Name of the tracking strategy. Only used if tracks keys are 3-tuples.
        If None and tracks keys are 4-tuples, strategies are auto-detected.
        (default: None)
    results_dir : Path or str, optional
        Directory to save results. If None, uses ClimateDT/results/tracking relative to script.
        (default: None)
    plot_mode : str, optional
        Legacy plotting option. Prefer calling plot_comparison_tracks separately.
        Accepts 'none', 'all', or 'low_correlation'. (default: 'none')
    min_ref_lat : float or None, optional
        Minimum matched IBTrACS latitude, inclusive. Defaults to -40 (40°S);
        use None to disable the latitude filter.
    min_ref_sshs : int or None, optional
        Minimum matched IBTrACS usa_sshs category, inclusive. None disables
        the filter; missing categories are excluded when it is enabled.
    scen : str, optional
        Scenario to filter ('hist', 'cont', 'Tplus2.0K', etc.). If None, processes all scenarios.
    land_points : str, optional
        How to handle landfall points: 'excl' excludes points over land, None includes all.
        Requires 'should_exclude' column in tracks data.
    
    Returns
    -------
    pd.DataFrame
        Combined comparison results for all storms with scenario, realization, 
        and strategy columns. One row per (scenario, realization, strategy, storm) combination.
    Path
        Directory where results were saved
    """
    if results_dir is None:
        # Try to use __file__ if available (script mode)
        try:
            results_dir = Path(__file__).parent.parent.parent / "results" / "tracking"
        except NameError:
            # Fallback for Jupyter/IPython (no __file__)
            results_dir = Path.cwd() / "results" / "tracking"
            print(f"  (Running in notebook mode, using: {results_dir})")
    else:
        results_dir = Path(results_dir)
    
    results_dir.mkdir(parents=True, exist_ok=True)
    
    # Auto-detect if multiple strategies in tracks by checking sample key length
    if tracks:
        sample_key = next(iter(tracks.keys()))
        is_multi_strategy = len(sample_key) == 4  # 4-tuple = multi-strategy
    else:
        is_multi_strategy = False
    
    if is_multi_strategy and tracking_strategy is not None:
        print(f"Warning: tracking_strategy parameter ignored (auto-detected from tracks)")
    
    print("=" * 80)
    if is_multi_strategy:
        print("TRACKER COMPARISON - MULTIPLE STRATEGIES")
    else:
        if tracking_strategy is None:
            tracking_strategy = "2-point-init"
        print(f"TRACKER COMPARISON - {tracking_strategy.upper()}")
    print("=" * 80)

    all_results = []
    n_storms_compared = 0

    for key, tracked_track in tracks.items():
        if tracked_track is None:
            continue
        
        # Parse key based on actual format (handle both 3-tuple and 4-tuple)
        if len(key) == 4:
            # Multi-strategy: (strategy, scenario, realization, storm_name)
            strategy, scen_key, realization, storm_name = key
        elif len(key) == 3:
            # Single strategy: (scenario, realization, storm_name)
            scen_key, realization, storm_name = key
            strategy = tracking_strategy
        else:
            print(f"Unexpected key format (len={len(key)}): {key}, skipping")
            continue

        # Filter by scenario if specified
        if scen is not None and scen_key != scen:
            continue

        n_storms_compared += 1
        reference = storms[storm_name]

        comparison_df = compare_all_strategies(
            storm_name, {strategy: tracked_track}, reference,
            min_ref_lat=min_ref_lat, min_ref_sshs=min_ref_sshs,
            filter_metrics=True,  # Pass flag to filter metrics
            land_points=land_points,  # Pass land filtering option
        )

        if len(comparison_df) > 0:
            comparison_df['scenario'] = scen_key
            comparison_df['realization'] = realization
            all_results.append(comparison_df)
            
    print(f"\nAttempted {n_storms_compared} track comparisons")
        
    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)
    
    # Combine all results
    if all_results:
        df_combined = pd.concat(all_results, ignore_index=True)
        
        # Reorder columns: scenario, realization, strategy first
        cols = df_combined.columns.tolist()
        for col in ['scenario', 'realization', 'strategy']:
            if col in cols:
                cols.remove(col)
        df_combined = df_combined[['scenario', 'realization', 'strategy'] + cols]
        
        # Remove redundant columns
        cols_to_keep = [c for c in df_combined.columns 
                       if c not in ['n_tracked_points', 'n_reference_points', 'data_coverage']]
        df_combined = df_combined[cols_to_keep]
        
        print("\n" + "-" * 80)
        print("METRICS SUMMARY")
        print("-" * 80)
        
        # Summary statistics
        summary_cols = ['MDE_km', 'CEP_km', 'RMSE_pres_hPa', 'bias_pres_hPa']
        for col in summary_cols:
            if col in df_combined.columns:
                data = df_combined[col].dropna()
                # if len(data) > 0:
                #     print(f"{col:20s}: mean={data.mean():8.2f}, std={data.std():8.2f}, "
                #           f"min={data.min():8.2f}, max={data.max():8.2f}")
        
        print(f"\nTotal comparisons: {len(df_combined)}")
        
        if plot_mode != 'none':
            plot_comparison_tracks(
                df_combined, storms, tracks, results_dir, plot_mode,
                min_ref_lat=min_ref_lat, min_ref_sshs=min_ref_sshs,
            )
        
        return df_combined
    
    return df_combined


def plot_tracks_simple(track, tc_idai_climatedt, scenarios):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.pyplot as plt
    bbox_idai = [30, 45, -22, -12] 

    fig, axes = plt.subplots(2, 2, figsize=(8, 6), dpi=300, constrained_layout=True,
                             subplot_kw={"projection": ccrs.PlateCarree()}, 
                             sharex=True, sharey=True)

    axes = axes.flatten()
    for i, (ax, scen) in enumerate(zip(axes, scenarios)):
        # background map
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax.add_feature(cfeature.LAND, facecolor="lightgray")
        ax.add_feature(cfeature.OCEAN, facecolor="lightblue")

        # IBTrACS track
        ax.plot(track["lon"], track["lat"], "k-", lw=2, transform=ccrs.PlateCarree(),
                label="IBTrACS")

        # realizations
        for (scenario, realization), min_p in tc_idai_climatedt.items():
            if scenario == scen:
                ax.plot(min_p["lon"], min_p["lat"], alpha=0.5, lw=1, transform=ccrs.PlateCarree(),
                        label=f"r{realization}")

        ax.set_extent(bbox_idai, crs=ccrs.PlateCarree())
        ax.set_title(scen)

        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        # Left labels only on first column (plots 0 and 2)
        gl.left_labels = (i % 2 == 0)
        gl.bottom_labels = (i != 0)

    # Use unused axes for legend
    if len(scenarios) < len(axes):
        legend_ax = axes[len(scenarios)]
        legend_ax.axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    legend_ax.legend(handles, labels, loc="upper left", frameon=True)

    plt.tight_layout()
    plt.show()


def plot_tracks_all(storms, tc_climatedt, figure_path="figures"):
    """
    Plot spatial tracks for storms across scenarios and strategies.
    
    Supports both single and multiple strategies:
    - Single strategy: keys = (scenario, realization, storm_name) → 1 row × 3 scenarios
    - Multiple strategies: keys = (strategy, scenario, realization, storm_name) → n_strategies rows × 3 scenarios
    
    Parameters
    ----------
    storms : dict
        Dictionary of reference storm tracks (IBTrACS)
    tc_climatedt : dict
        Dictionary of tracked storm data with keys:
        - (scenario, realization, storm_name) for single strategy
        - (strategy, scenario, realization, storm_name) for multiple strategies
    figure_path : str
        Directory to save figures
    """
    os.makedirs(figure_path, exist_ok=True)
    
    scenarios = ["cont", "hist", "Tplus2.0K"]
    colors = plt.cm.tab10.colors
    
    # Auto-detect if multiple strategies
    if tc_climatedt:
        sample_key = next(iter(tc_climatedt.keys()))
        is_multi_strategy = len(sample_key) == 4
    else:
        is_multi_strategy = False
    
    # Extract unique storms and strategies
    unique_storms = set()
    strategies = set()
    
    for key in tc_climatedt.keys():
        if is_multi_strategy:
            strategy, scen, realization, name = key
            strategies.add(strategy)
        else:
            scen, realization, name = key
        unique_storms.add(name)
    
    strategies = sorted(list(strategies)) if strategies else [None]
    
    for storm_name in sorted(unique_storms):
        print(f"Plotting tracks for {storm_name}")
        
        # Skip if not in reference
        if storm_name not in storms:
            # print(f"  {storm_name}: not in IBTrACS reference, skipping.")
            continue
        
        track = storms[storm_name]
        
        # Organize tracks by strategy and scenario/realization
        tracks_by_strategy = {}
        for key, model_track in tc_climatedt.items():
            if is_multi_strategy:
                strategy, scen, realization, name = key
            else:
                scen, realization, name = key
                strategy = strategies[0]
            
            if name != storm_name:
                continue
            
            if strategy not in tracks_by_strategy:
                tracks_by_strategy[strategy] = {}
            
            tracks_by_strategy[strategy][(scen, realization)] = model_track
        
        if not tracks_by_strategy:
            # print(f"  No model tracks found for {storm_name}, skipping.")
            continue
        
        # Shared extent based on IBTrACS and all ClimateDT tracks
        all_lon = [track["lon"].to_numpy()]
        all_lat = [track["lat"].to_numpy()]
        
        for strategy_tracks in tracks_by_strategy.values():
            for model_track in strategy_tracks.values():
                if model_track is not None:
                    all_lon.append(model_track["lon"].to_numpy())
                    all_lat.append(model_track["lat"].to_numpy())
        
        all_lon = np.concatenate(all_lon)
        all_lat = np.concatenate(all_lat)
        
        # Check if there are any valid coordinates
        valid_lon = all_lon[np.isfinite(all_lon)]
        valid_lat = all_lat[np.isfinite(all_lat)]
        
        if len(valid_lon) == 0 or len(valid_lat) == 0:
            print(f"  Skipping {storm_name}: no valid coordinates")
            continue
        
        margin = 3
        extent = [np.nanmin(all_lon) - margin,
                  np.nanmax(all_lon) + margin,
                  np.nanmin(all_lat) - margin,
                  np.nanmax(all_lat) + margin]
        
        # Create subplot grid: one row per strategy, 3 columns for scenarios
        n_strategies = len(tracks_by_strategy)
        figsize = (15, 4 + (n_strategies - 1) * 4)  # Adjust height based on number of strategies
        
        fig, axes = plt.subplots(n_strategies, 3, figsize=figsize,
                                 subplot_kw={"projection": ccrs.PlateCarree()},
                                 constrained_layout=True)
        
        # Handle case where there's only one strategy (axes is 1D)
        if n_strategies == 1:
            axes = axes.reshape(1, -1)
        
        for strategy_idx, (strategy, model_tracks) in enumerate(sorted(tracks_by_strategy.items())):
            for scen_idx, scen in enumerate(scenarios):
                ax = axes[strategy_idx, scen_idx]
                
                ax.coastlines()
                ax.add_feature(cfeature.BORDERS, linewidth=0.5)
                ax.add_feature(cfeature.LAND, facecolor="lightgray")
                ax.add_feature(cfeature.OCEAN, facecolor="lightblue")
                
                # Plot IBTrACS reference
                ax.plot(track["lon"], track["lat"], "k-", lw=2.5, transform=ccrs.PlateCarree(),
                       label="IBTrACS")
                
                # Start of track
                ax.plot(track.lon.iloc[0], track.lat.iloc[0], marker="o", color="red", markersize=3,
                       transform=ccrs.PlateCarree(), zorder=10, label="Track start")
                
                # Plot ensemble members
                for i, realization in enumerate(["1", "2", "3", "4", "5"]):
                    key = (scen, realization)
                    
                    if key not in model_tracks:
                        continue
                    
                    model_track = model_tracks[key]
                    
                    if model_track is not None:
                        ax.plot(model_track["lon"], model_track["lat"], color=colors[i], lw=1.3,
                               alpha=0.8, transform=ccrs.PlateCarree(), label=f"r{realization}")
                
                ax.set_extent(extent, crs=ccrs.PlateCarree())
                
                # Set title: scenario name on top row, strategy name on left column
                if strategy_idx == 0:
                    ax.set_title(scen, fontweight='bold')
                if scen_idx == 0:
                    ax.text(-0.45, 0.5, strategy, transform=ax.transAxes,
                           fontsize=11, fontweight='bold', va='center', rotation=90)
                
                gl = ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5)
                gl.top_labels = False
                gl.right_labels = False
                if strategy_idx < n_strategies - 1:
                    gl.bottom_labels = False
        
        # Add metadata
        max_wind = track.wind.max() * 0.514444  # convert from knots to m/s
        
        # TC category at maximum wind
        idx_max = track.wind.idxmax()
        tc_strength = track.usa_sshs.loc[idx_max]
        start_date = pd.to_datetime(track.time.min()).strftime("%Y-%m-%d")
        
        # Max intensity point (mark on middle subplot of first strategy row)
        axes[0, 1].scatter(track.lon.loc[idx_max], track.lat.loc[idx_max], color="white", s=30,
                          marker="x", edgecolor="black", linewidth=0.5, transform=ccrs.PlateCarree(),
                          zorder=5, label=f"Max intensity (SSHS={tc_strength})")
        
        fig.suptitle(f"{storm_name} | {start_date} | Max wind: {max_wind:.1f} m/s | {tc_strength}",
                    fontsize=14, fontweight='bold')
        
        # Unified legend at bottom
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=7, frameon=True, fontsize=9)
        
        safe_name = storm_name.replace(" ", "_").replace("/", "_")
        outfile = os.path.join(figure_path, f"{safe_name}_tracked.png")
        
        fig.savefig(outfile, bbox_inches="tight", dpi=150)
        print(f"  Saved to: {outfile}")
        
        if platform.system() == "Windows":
            plt.show()
        
        plt.close(fig)


def plot_comparison_tracks(results_df, storms, tracks, results_dir, plot_mode='none',
                           min_ref_lat=-40, min_ref_sshs=None):
    """Plot all or low-correlation storms, filtering tracked points by IBTrACS.

    The full reference track is shown for context. Figures are saved under
    ``results_dir/figures/tracking``. Use ``plot_mode='none'`` to skip plotting.
    """
    if plot_mode not in ('none', 'all', 'low_correlation'):
        raise ValueError(f"Unknown plot_mode: {plot_mode!r}")
    if plot_mode == 'none' or results_df is None or results_df.empty:
        return

    print("\n" + "-" * 80)
    print("GENERATING SPATIAL TRACK PLOTS")
    print("-" * 80)

    if plot_mode == 'low_correlation' and 'corr_pres' in results_df.columns:
        threshold = results_df['corr_pres'].quantile(0.1)
        plotted_storms = results_df[results_df['corr_pres'] <= threshold]['storm'].unique()
        print(f"Found {len(plotted_storms)} storms with low correlation (r <= {threshold:.3f})")
    else:
        plotted_storms = results_df['storm'].unique()
        if plot_mode == 'low_correlation':
            print(f"'corr_pres' column not found, plotting all {len(plotted_storms)} storms")
        else:
            print(f"Plotting all {len(plotted_storms)} storms")

    tracks_to_plot = {}
    for key, model_track in tracks.items():
        if model_track is None or len(key) not in (3, 4):
            continue
        storm_name = key[-1]
        if storm_name not in plotted_storms:
            continue
        if min_ref_lat is not None or min_ref_sshs is not None:
            aligned = align_tracks(model_track, storms[storm_name])
            if aligned is None:
                continue
            if min_ref_lat is not None:
                aligned = aligned[aligned['lat_ref'] >= min_ref_lat]
            if min_ref_sshs is not None:
                aligned = aligned[aligned['usa_sshs_ref'] >= min_ref_sshs]
            if aligned.empty:
                continue
            model_track = model_track[model_track["time"].isin(aligned["time"])].copy()
        tracks_to_plot[key] = model_track

    if tracks_to_plot:
        figure_dir = Path(results_dir) / "figures" / "tracking"
        figure_dir.mkdir(parents=True, exist_ok=True)
        print(f"Plotting {len(plotted_storms)} storms...")
        plot_tracks_all(storms, tracks_to_plot, str(figure_dir))
        print(f"Saved track plots to: {figure_dir}")


def save_comparison_figures(figures, results_dir, filename):
    """Save a single figure or the per-strategy figures returned by plot helpers."""
    if figures is None:
        return
    if not isinstance(figures, dict):
        figures = {None: figures}
    for strategy, figure in figures.items():
        safe_strategy = str(strategy).replace('/', '_').replace('\\', '_')
        suffix = f"_{safe_strategy}" if strategy is not None else ""
        plot_file = results_dir / f"{filename}{suffix}.png"
        figure.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close(figure)
        print(f"Saved to: {plot_file}")

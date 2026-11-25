"""
Comparison utilities for different storm tracking strategies.

Provides functions to compare tracking results and compute performance metrics.
"""

import numpy as np
import pandas as pd
import xarray as xr
import warnings
from scipy.spatial.distance import cdist


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


def aggregate_comparison_by_storm(all_results):
    """
    Aggregate comparison results across multiple storms.
    
    Parameters
    ----------
    all_results : dict
        Dictionary mapping storm name to comparison DataFrame
    
    Returns
    -------
    dict
        Dictionary with aggregated statistics by strategy
    """
    # Concatenate all results
    df_all = pd.concat(all_results.values(), ignore_index=True)
    
    # Group by strategy and compute statistics
    by_strategy = df_all.groupby('strategy').agg({
        'n_points': ['count', 'mean'],
        'MDE_km': ['mean', 'std', 'min', 'max'],
        'CEP_km': ['mean', 'std'],
        'MXE_km': ['mean', 'std', 'max'],
        'P90E_km': ['mean', 'std'],
        'RMSE_pres_hPa': ['mean', 'std'],
        'bias_pres_hPa': ['mean', 'std'],
        'corr_pres': ['mean', 'std'],
        'completion_rate': ['mean', 'std', 'min', 'max'],
    })
    
    return by_strategy


def plot_metrics_distribution(results_df, figsize=None):
    """
    Plot distribution of key tracking metrics for easy analysis.
    Creates separate plots for each strategy if multiple strategies present.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Results from compare_all_strategies with columns: strategy, MDE_km, MXE_km,
        position_RMSE_km, RMSE_pres_hPa, bias_pres_hPa, corr_pres
    figsize : tuple, optional
        Figure size per strategy (width, height). Default: (15, 8)
    
    Returns
    -------
    dict or matplotlib.figure.Figure
        Single figure if one strategy, dict of figures if multiple strategies
    """
    import matplotlib.pyplot as plt
    
    if figsize is None:
        figsize = (15, 8)
    
    # Key metrics to plot (6 total)
    metrics = [
        ('MDE_km', 'Mean Position Error (km)', 15),
        ('MXE_km', 'Maximum Position Error (km)', 15),
        ('position_RMSE_km', 'Position RMSE (km)', 10),
        ('RMSE_pres_hPa', 'Pressure RMSE (hPa)', 15),
        ('bias_pres_hPa', 'Pressure Bias (hPa)', 15),
        ('corr_pres', 'Pressure Correlation', 20),
    ]
    
    # Detect strategies
    strategies = results_df['strategy'].unique() if 'strategy' in results_df.columns else [None]
    
    figs = {}
    
    for strategy in strategies:
        # Filter data by strategy
        if strategy is not None:
            data_to_plot = results_df[results_df['strategy'] == strategy]
            title_prefix = f'{strategy}'
        else:
            data_to_plot = results_df
            title_prefix = 'All Data'
        
        if len(data_to_plot) == 0:
            print(f"  No data for strategy: {strategy}")
            continue
        
        # Create figure with 2 rows x 3 columns
        fig, axes = plt.subplots(2, 3, figsize=figsize, constrained_layout=True)
        fig.suptitle(f'Tracking Metrics Distribution - {title_prefix}', fontsize=14, fontweight='bold')
        
        axes = axes.flatten()
        
        for idx, (metric, label, bins) in enumerate(metrics):
            ax = axes[idx]
            
            if metric not in data_to_plot.columns:
                ax.text(0.5, 0.5, f'{metric}\nNot available', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=10)
                ax.set_title(label)
                continue
            
            data = data_to_plot[metric].dropna()
            
            if len(data) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(label)
                continue
            
            # Histogram with statistics
            ax.hist(data, bins=bins, alpha=0.7, edgecolor='black', color='steelblue')
            ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2, 
                      label=f'Mean: {data.mean():.2f}')
            ax.axvline(data.median(), color='green', linestyle='--', linewidth=2, 
                      label=f'Median: {data.median():.2f}')
            
            ax.set_xlabel(label)
            ax.set_ylabel('Frequency')
            ax.set_title(label, fontsize=11, fontweight='bold')
            ax.legend(fontsize=8, loc='upper right')
            ax.grid(axis='y', alpha=0.3)
            
            # Statistics text box
            stats_text = (f'N={len(data)}\nMean={data.mean():.2f}\n'
                         f'Std={data.std():.2f}\nMin={data.min():.2f}\nMax={data.max():.2f}')
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                   fontsize=8, verticalalignment='top', horizontalalignment='left',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
        
        figs[strategy if strategy is not None else 'all'] = fig
    
    # Return single figure or dict of figures
    if len(figs) == 1:
        return list(figs.values())[0]
    else:
        return figs


def plot_metrics_distribution_filtered(results_df, storms, min_sshs=1, figsize=None):
    """
    Plot metrics distribution for storms with usa_sshs >= threshold in reference.
    Creates separate plots for each strategy if multiple strategies present.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Results from compare_all_strategies with strategy column
    storms : dict
        Dictionary of reference storm tracks (IBTrACS) with 'usa_sshs' data
    min_sshs : int, optional
        Minimum SSHS category to include (default: 1, hurricane strength)
    figsize : tuple, optional
        Figure size per strategy (width, height). Default: (15, 8)
    
    Returns
    -------
    dict or matplotlib.figure.Figure
        Single figure if one strategy, dict of figures if multiple strategies
    """
    import matplotlib.pyplot as plt
    
    if figsize is None:
        figsize = (15, 8)
    
    # Filter to storms with usa_sshs >= min_sshs
    filtered_results = []
    for _, row in results_df.iterrows():
        storm_name = row['storm']
        if storm_name not in storms:
            continue
        
        track = storms[storm_name]
        if 'usa_sshs' in track.columns:
            max_sshs = track['usa_sshs'].max()
            if pd.notna(max_sshs) and max_sshs >= min_sshs:
                filtered_results.append(row)
    
    if not filtered_results:
        print(f"No storms found with usa_sshs >= {min_sshs}")
        return None
    
    filtered_df = pd.DataFrame(filtered_results)
    
    # Key metrics to plot (6 total)
    metrics = [
        ('MDE_km', 'Mean Position Error (km)', 15),
        ('MXE_km', 'Maximum Position Error (km)', 15),
        ('position_RMSE_km', 'Position RMSE (km)', 10),
        ('RMSE_pres_hPa', 'Pressure RMSE (hPa)', 15),
        ('bias_pres_hPa', 'Pressure Bias (hPa)', 15),
        ('corr_pres', 'Pressure Correlation', 20),
    ]
    
    # Detect strategies
    strategies = filtered_df['strategy'].unique() if 'strategy' in filtered_df.columns else [None]
    
    figs = {}
    
    for strategy in strategies:
        # Filter data by strategy
        if strategy is not None:
            data_to_plot = filtered_df[filtered_df['strategy'] == strategy]
            title_prefix = f'{strategy}'
        else:
            data_to_plot = filtered_df
            title_prefix = 'All Data'
        
        if len(data_to_plot) == 0:
            print(f"  No data for strategy: {strategy}")
            continue
        
        # Create figure with 2 rows x 3 columns
        fig, axes = plt.subplots(2, 3, figsize=figsize, constrained_layout=True)
        fig.suptitle(f'Tracking Metrics (SSHS ≥ {min_sshs}) - {title_prefix} (N={len(data_to_plot)} storms)', 
                    fontsize=14, fontweight='bold')
        
        axes = axes.flatten()
        
        for idx, (metric, label, bins) in enumerate(metrics):
            ax = axes[idx]
            
            if metric not in data_to_plot.columns:
                ax.text(0.5, 0.5, f'{metric}\nNot available', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=10)
                ax.set_title(label)
                continue
            
            data = data_to_plot[metric].dropna()
            
            if len(data) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(label)
                continue
            
            # Histogram with statistics
            ax.hist(data, bins=bins, alpha=0.7, edgecolor='black', color='coral')
            ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2, 
                      label=f'Mean: {data.mean():.2f}')
            ax.axvline(data.median(), color='green', linestyle='--', linewidth=2, 
                      label=f'Median: {data.median():.2f}')
            
            ax.set_xlabel(label)
            ax.set_ylabel('Frequency')
            ax.set_title(label, fontsize=11, fontweight='bold')
            ax.legend(fontsize=8, loc='upper right')
            ax.grid(axis='y', alpha=0.3)
            
            # Statistics text box
            stats_text = (f'N={len(data)}\nMean={data.mean():.2f}\n'
                         f'Std={data.std():.2f}\nMin={data.min():.2f}\nMax={data.max():.2f}')
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                   fontsize=8, verticalalignment='top', horizontalalignment='left',
                   bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))
        
        figs[strategy if strategy is not None else 'all'] = fig
    
    # Return single figure or dict of figures
    if len(figs) == 1:
        return list(figs.values())[0]
    else:
        return figs


def plot_strategy_comparison(results_df, figsize=(14, 6)):
    """
    Compare metrics across different tracking strategies.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Results from compare_all_strategies
    figsize : tuple, optional
        Figure size (width, height)
    
    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    import matplotlib.pyplot as plt
    
    if 'strategy' not in results_df.columns:
        print("No strategy column found in results")
        return None
    
    strategies = results_df['strategy'].unique()
    if len(strategies) == 1:
        print(f"Only one strategy found: {strategies[0]}")
        return None
    
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    fig.suptitle(f'Strategy Comparison: {", ".join(strategies)}', fontsize=14, fontweight='bold')
    
    # Position error comparison
    ax = axes[0, 0]
    for strategy in strategies:
        data = results_df[results_df['strategy'] == strategy]['MDE_km']
        ax.hist(data.dropna(), bins=15, alpha=0.6, label=strategy, edgecolor='black')
    ax.set_xlabel('Position Error MDE (km)')
    ax.set_ylabel('Frequency')
    ax.set_title('Position Error Distribution')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # Pressure error comparison
    ax = axes[0, 1]
    for strategy in strategies:
        data = results_df[results_df['strategy'] == strategy]['RMSE_pres_hPa']
        ax.hist(data.dropna(), bins=15, alpha=0.6, label=strategy, edgecolor='black')
    ax.set_xlabel('Pressure RMSE (hPa)')
    ax.set_ylabel('Frequency')
    ax.set_title('Pressure Error Distribution')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # Box plot comparison - position metrics
    ax = axes[1, 0]
    position_metrics = ['MDE_km', 'CEP_km', 'position_RMSE_km']
    for i, metric in enumerate(position_metrics):
        if metric in results_df.columns:
            data = results_df[metric].dropna()
            ax.boxplot(data, positions=[i], widths=0.6, patch_artist=True,
                      boxprops=dict(facecolor='lightblue', alpha=0.7))
    ax.set_xticklabels(position_metrics)
    ax.set_ylabel('Error (km)')
    ax.set_title('Position Metrics Comparison')
    ax.grid(axis='y', alpha=0.3)
    
    # Box plot comparison - intensity metrics
    ax = axes[1, 1]
    intensity_metrics = ['RMSE_pres_hPa', 'bias_pres_hPa', 'corr_pres']
    positions = []
    for i, metric in enumerate(intensity_metrics):
        if metric in results_df.columns:
            data = results_df[metric].dropna()
            ax.boxplot(data, positions=[i], widths=0.6, patch_artist=True,
                      boxprops=dict(facecolor='lightcoral', alpha=0.7))
            positions.append(i)
    ax.set_xticklabels([intensity_metrics[i] for i in positions])
    ax.set_ylabel('Magnitude')
    ax.set_title('Intensity Metrics Comparison')
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    return fig


if __name__ == "__main__":
    # Example usage
    print("Tracker comparison utilities module")
    print("Import and use functions for analyzing tracking results")

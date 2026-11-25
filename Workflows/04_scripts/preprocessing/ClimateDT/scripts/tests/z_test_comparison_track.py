#%%
"""Compare existing ClimateDT tracks against IBTrACS reference data."""

import os
import pickle
from pathlib import Path
import platform
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from z_tracker_comparison import align_tracks, compare_all_strategies
from z_general_TC_funcs import import_ibtracs_period


def load_existing_tracks(track_file):
    print(f"Loading existing tracked results: {track_file}", flush=True)
    with open(track_file, "rb") as f:
        tracks = pickle.load(f)
    
    return tracks


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
                all_lon.append(model_track["lon"].to_numpy())
                all_lat.append(model_track["lat"].to_numpy())
        
        all_lon = np.concatenate(all_lon)
        all_lat = np.concatenate(all_lat)
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


def plot_comparison_metrics(results_df, storms, results_dir):
    """Generate and save metric distributions and strategy comparison."""
    from tracker_comparison import (
        plot_metrics_distribution,
        plot_metrics_distribution_filtered,
        plot_strategy_comparison,
    )

    print("\nGenerating metrics distribution plots...")
    save_comparison_figures(plot_metrics_distribution(results_df), results_dir,
                            "metrics_distribution")
    save_comparison_figures(
        plot_metrics_distribution_filtered(results_df, storms, min_sshs=1),
        results_dir, "metrics_distribution_sshs1plus",
    )
    if results_df['strategy'].nunique() > 1:
        save_comparison_figures(plot_strategy_comparison(results_df), results_dir,
                                "strategy_comparison")


if __name__ == "__main__":
    if platform.system() == "Windows":
        # processed_path = Path("P:/11210471-001-compass/01_Data/ECMWF_ClimateDT/analysis_output/SI_tracked_storms/")
        processed_path = Path("C:/Code/COMPASS/Workflows/04_scripts/preprocessing/ClimateDT/data/processed/")
        track_file = (processed_path / "storms_all_strategies_2017-01-01_2026-08-01.pkl")
        wind_files = {
            "200km": processed_path / "refreshed_winds" / "storms_all_ibtracs_max_wind_200km_2017-01-01_2026-08-01.pkl",
            "rmw": processed_path / "refreshed_winds" / "storms_all_ibtracs_max_wind_rmw_2017-01-01_2026-08-01.pkl",
        }
    else:
        track_file = Path("/projects/prjs2226/data/ClimateDT/sfc/processed/msl/tracked/"
                          "storms_all_strategies_2017-01-01_2026-08-01.pkl")
        wind_files = {}  # Set paths here to include wind data.
    
    start_time = "2017-01-01"
    end_time = "2026-08-01"
    plot_mode = 'none'  # 'none', 'all', or 'low_correlation' (spatial tracks)
    plot_metrics = True  # Generate metrics distributions and strategy comparison
    min_ref_lat = -40  # Inclusive latitude cutoff; None disables it
    min_ref_sshs = 1  # Set to 1 to compare only category 1+ points

    # Load IBTrACS reference
    print("\nLoading IBTrACS reference data...", flush=True)
    storms = import_ibtracs_period(start=start_time, end=end_time)
    print(f"Loaded {len(storms)} reference storms")

    # Load existing tracked results
    print("\nLoading existing tracked results...", flush=True)
    tracks = load_existing_tracks(track_file)

    # Keep all aligned points available for inspection before comparison.
    aligned_tracks = create_aligned_tracks_dataframe(storms, tracks)
    
    print(f"Aligned track points available for inspection: {len(aligned_tracks)}")
    if not aligned_tracks.empty:
        if not wind_files:
            print("No wind files configured; aligned tracks will not include model wind")
        for label, wind_file in wind_files.items():
            aligned_tracks = add_aligned_wind(aligned_tracks, load_existing_tracks(wind_file), label)
            aligned_tracks.to_csv(Path(track_file).parent / "aligned_tracks_inspection.csv", index=False)

    # Get results as dataframe
    results_df, results_dir = compare_climatedt_tracks(
        storms, tracks, 
        # tracking_strategy="2-point-init", # disabled if multiple strategies are present
        min_ref_lat=min_ref_lat,
        min_ref_sshs=min_ref_sshs,
    )
    
    if results_df is not None and len(results_df) > 0:
        plot_comparison_tracks(
            results_df, storms, tracks, results_dir, plot_mode,
            min_ref_lat=min_ref_lat, min_ref_sshs=min_ref_sshs,
        )
        print(f"\nResults dataframe shape: {results_df.shape}")
        print(f"Columns: {results_df.columns.tolist()}")
        
        # Save results to results directory
        output_csv = results_dir / "tracking_comparison_results.csv"
        results_df.to_csv(output_csv, index=False)
        print(f"\nDetailed results saved to: {output_csv}")
        
        if plot_metrics:
            plot_comparison_metrics(results_df, storms, results_dir)
        
        print(f"\nAll results saved to: {results_dir}")
    else:
        print("No comparison results were produced")

# %%

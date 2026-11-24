"""
Quick start examples for tracker strategies.

Shows how to use the different tracking strategies and comparison utilities.
"""

import pandas as pd
import xarray as xr
from pathlib import Path
import pickle

# Import tracking strategies
from tracker_strategies import (
    strategy_2point_init,
    strategy_5point_init,
    strategy_all_ibtracs_init
)

# Import comparison utilities
from tracker_comparison import (
    compare_all_strategies,
    compute_position_metrics,
    compute_intensity_metrics,
    aggregate_comparison_by_storm,
)


# ============================================================================
# EXAMPLE 1: Track a single storm with all strategies
# ============================================================================

def example_single_storm_comparison():
    """
    Track one storm using all four strategies and compare results.
    """
    print("\n" + "="*70)
    print("EXAMPLE 1: Single Storm Comparison")
    print("="*70)
    
    # Load data (adjust paths as needed)
    data = xr.open_mfdataset(
        "/projects/prjs2226/data/ClimateDT/sfc/raw/hist/r1/msl/*.nc",
        combine="by_coords"
    )
    
    # Load IBTrACS reference
    from tracking_TCs_SI import import_ibtracs_period
    storms = import_ibtracs_period(start="2020-01-01", end="2020-12-31")
    
    # Select one storm
    storm_name = list(storms.keys())[0]
    best_track = storms[storm_name]
    
    print(f"\nTracking storm: {storm_name}")
    print(f"IBTrACS track points: {len(best_track)}")
    
    # Track with each strategy
    results = {}
    strategies = {
        '2point': strategy_2point_init,
        '5point': strategy_5point_init,
        'all_ibtracs': strategy_all_ibtracs_init,
        'intensity': strategy_intensity_init,
    }
    
    for strategy_name, strategy_func in strategies.items():
        try:
            print(f"\n  Running {strategy_name}...", end=" ", flush=True)
            track = strategy_func(data, best_track, s_size=5)
            results[strategy_name] = track
            print(f"SUCCESS ({len(track)} points)")
        except Exception as e:
            print(f"FAILED: {e}")
            results[strategy_name] = None
    
    # Compare results
    print("\n" + "-"*70)
    print("COMPARISON RESULTS")
    print("-"*70)
    
    comparison_df = compare_all_strategies(storm_name, results, best_track)
    print("\n" + comparison_df.to_string())
    
    return comparison_df


# ============================================================================
# EXAMPLE 2: Batch processing with results summary
# ============================================================================

def example_batch_processing():
    """
    Process multiple storms and create summary comparison.
    """
    print("\n" + "="*70)
    print("EXAMPLE 2: Batch Processing Summary")
    print("="*70)
    
    # Load data
    data = xr.open_mfdataset(
        "/projects/prjs2226/data/ClimateDT/sfc/raw/hist/r1/msl/*.nc",
        combine="by_coords"
    )
    
    # Load IBTrACS
    from tracking_TCs_SI import import_ibtracs_period
    storms = import_ibtracs_period(start="2020-01-01", end="2020-12-31")
    
    print(f"\nProcessing {len(storms)} storms")
    
    # Track with each strategy
    strategies = {
        '2point': strategy_2point_init,
        '5point': strategy_5point_init,
        'all_ibtracs': strategy_all_ibtracs_init,
        'intensity': strategy_intensity_init,
    }
    
    all_comparisons = {}
    
    for i, (storm_name, best_track) in enumerate(storms.items(), 1):
        if i > 5:  # Process only first 5 for demo
            break
        
        print(f"\n  [{i}/{min(5, len(storms))}] {storm_name}...", end=" ", flush=True)
        
        results = {}
        for strategy_name, strategy_func in strategies.items():
            try:
                track = strategy_func(data, best_track, s_size=5)
                results[strategy_name] = track
            except Exception as e:
                results[strategy_name] = None
        
        comp_df = compare_all_strategies(storm_name, results, best_track)
        all_comparisons[storm_name] = comp_df
        print(f"DONE")
    
    # Aggregate results
    print("\n" + "-"*70)
    print("AGGREGATED RESULTS")
    print("-"*70)
    
    summary = aggregate_comparison_by_storm(all_comparisons)
    print("\n" + summary.to_string())
    
    return all_comparisons, summary


# ============================================================================
# EXAMPLE 3: Load and analyze pre-computed results
# ============================================================================

def example_load_results():
    """
    Load previously computed strategy comparison results.
    """
    print("\n" + "="*70)
    print("EXAMPLE 3: Load Pre-computed Results")
    print("="*70)
    
    from run_tracker_strategies import load_strategy_results
    
    processed_path_tracks = Path("/projects/prjs2226/data/ClimateDT/sfc/processed/msl/tracked")
    
    # Load results for hist r1
    results = load_strategy_results(processed_path_tracks, scen="hist", realization="1")
    
    print(f"\nLoaded results for:")
    for scen, scen_data in results.items():
        for real, real_data in scen_data.items():
            print(f"  {scen} r{real}")
            for strategy_name, tracks in real_data.items():
                n_success = sum(1 for v in tracks.values() if v is not None)
                print(f"    {strategy_name}: {n_success}/{len(tracks)} storms")
    
    return results


# ============================================================================
# EXAMPLE 4: Detailed metric analysis
# ============================================================================

def example_detailed_metrics():
    """
    Compute detailed metrics for individual tracks.
    """
    print("\n" + "="*70)
    print("EXAMPLE 4: Detailed Metric Analysis")
    print("="*70)
    
    # Load data
    data = xr.open_mfdataset(
        "/projects/prjs2226/data/ClimateDT/sfc/raw/hist/r1/msl/*.nc",
        combine="by_coords"
    )
    
    # Load IBTrACS
    from tracking_TCs_SI import import_ibtracs_period
    storms = import_ibtracs_period(start="2020-01-01", end="2020-12-31")
    
    # Get first storm
    storm_name = list(storms.keys())[0]
    best_track = storms[storm_name]
    
    print(f"\nAnalyzing storm: {storm_name}")
    
    # Track with one strategy
    track = strategy_5point_init(data, best_track, s_size=5)
    
    # Compute detailed metrics
    print("\nPosition Metrics:")
    pos_metrics = compute_position_metrics(track, best_track)
    for key, value in pos_metrics.items():
        if isinstance(value, float):
            print(f"  {key:20s}: {value:8.2f}")
        else:
            print(f"  {key:20s}: {value}")
    
    print("\nIntensity Metrics:")
    int_metrics = compute_intensity_metrics(track, best_track)
    for key, value in int_metrics.items():
        if isinstance(value, float):
            print(f"  {key:20s}: {value:8.2f}")
        else:
            print(f"  {key:20s}: {value}")


# ============================================================================
# EXAMPLE 5: Save comparison results
# ============================================================================

def example_save_results():
    """
    Save comparison results to CSV for further analysis.
    """
    print("\n" + "="*70)
    print("EXAMPLE 5: Save Results to CSV")
    print("="*70)
    
    # Get example data (using example 2)
    all_comparisons, summary = example_batch_processing()
    
    # Concatenate all comparisons
    all_df = pd.concat(all_comparisons.values(), ignore_index=True)
    
    # Save to CSV
    output_file = "tracker_comparison_all_storms.csv"
    all_df.to_csv(output_file, index=False)
    print(f"\nSaved to {output_file}")
    
    # Save summary
    summary_file = "tracker_comparison_summary.csv"
    summary.to_csv(summary_file)
    print(f"Saved to {summary_file}")
    
    return all_df


# ============================================================================
# Strategy Descriptions
# ============================================================================

def print_strategy_descriptions():
    """
    Print descriptions of each tracking strategy.
    """
    descriptions = {
        'strategy_2point_init': """
        STRATEGY 1: 2-Point Initialization (CURRENT DEFAULT)
        ────────────────────────────────────────────────────
        • Uses first 2 IBTrACS observation points
        • Initializes search boxes around point 1
        • After point 2, extrapolates forward using previous position
        • Good for: Balanced initialization without over-constraining
        • Pros: Simple, established method, less dependent on early obs errors
        • Cons: May miss initial position if model has bias
        """,
        
        'strategy_5point_init': """
        STRATEGY 2: 5-Point Initialization
        ──────────────────────────────────
        • Uses first 5 IBTrACS observation points
        • Establishes more robust initial track direction
        • Provides better estimate of forward velocity
        • Good for: Storms with noisy early observations
        • Pros: Better initialization, captures track curvature early
        • Cons: More data-dependent, may fail if early points are bad
        """,
        
        'strategy_all_ibtracs_init': """
        STRATEGY 3: All IBTrACS Points Guidance
        ──────────────────────────────────────
        • Uses all available IBTrACS points
        • At each model timestep, searches around nearest IBTrACS point
        • Effectively allows full track to guide the search
        • Good for: Best-case scenario with perfect knowledge
        • Pros: Should give best possible tracking results
        • Cons: Unrealistic - assumes knowledge of future; use as benchmark
        """,
        
        'strategy_intensity_init': """
        STRATEGY 4: Intensity-Based Initialization
        ──────────────────────────────────────────
        • Finds point of maximum intensity (min MSLP) in IBTrACS
        • Tracks backward from peak intensity to genesis
        • Tracks forward from peak intensity to dissipation
        • Good for: Capturing well-developed storm core
        • Pros: Strong initialization at known intensity peak
        • Cons: May miss early/late track portions; intensity shifts
        """,
    }
    
    print("\n" + "="*70)
    print("STRATEGY DESCRIPTIONS")
    print("="*70)
    
    for name, desc in descriptions.items():
        print(desc)


if __name__ == "__main__":
    import sys
    
    print("\n" + "="*70)
    print("TRACKER STRATEGIES - QUICK START EXAMPLES")
    print("="*70)
    
    print_strategy_descriptions()
    
    print("\n" + "="*70)
    print("TO RUN EXAMPLES, CALL:")
    print("="*70)
    print("""
    from tracker_examples import *
    
    # Single storm comparison
    comp_df = example_single_storm_comparison()
    
    # Batch processing
    all_comp, summary = example_batch_processing()
    
    # Load pre-computed results
    results = example_load_results()
    
    # Detailed metrics
    example_detailed_metrics()
    
    # Save results
    all_df = example_save_results()
    """)

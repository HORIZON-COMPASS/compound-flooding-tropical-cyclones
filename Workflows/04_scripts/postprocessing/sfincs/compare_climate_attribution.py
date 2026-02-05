"""
Climate Attribution Analysis: Factual vs Counterfactual Flooding.

This script compares flood outputs between:
- Factual (CF0): Current climate precipitation (ERA5)
- Counterfactual (CF-8): Pre-industrial proxy (-8% precipitation)

The difference quantifies the contribution of climate change to flooding.

Produces:
1. Bilateral hmax difference map showing climate-attributable flooding
2. Summary statistics with percentage attribution
3. Bar chart comparing total flood metrics

Author: Generated for Durban 2022 climate attribution analysis
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import rioxarray as rxr
import contextily as ctx
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")


# ===== CONFIGURATION =====
EVENT_NAME = "Durban_April2022"
REGION = "durban"
RUNNAME = "Durban2022"

# Base paths
BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs")
OUTPUT_DIR = Path("/p/11210471-001-compass/04_Results/climate_attribution")

# Run configurations
RUNS = {
    "Factual (CF0)": {
        "folder": "event_precip_era5_hourly_CF0_no_wind",
        "cf_value": 0,
    },
    "Counterfactual (CF-8)": {
        "folder": "event_precip_era5_hourly_CF-8_no_wind",
        "cf_value": -8,
    },
}

# Plotting thresholds
FLOOD_THRESHOLD = 0.05  # meters
HMAX_VMAX = 2.0
DIFF_VMAX = 0.5


# ===== HELPER FUNCTIONS =====
def get_run_path(run_name: str) -> Path:
    """Get the run directory path."""
    config = RUNS[run_name]
    return BASE_RUN_PATH / REGION / RUNNAME / "sfincs" / config["folder"]


def load_hmax_tif(run_name: str) -> xr.DataArray:
    """Load the hmax TIFF file."""
    run_path = get_run_path(run_name)
    tif_path = run_path / "plot_output" / "sfincs_output_hmax_AllTime.tif"

    if not tif_path.exists():
        raise FileNotFoundError(f"TIFF file not found: {tif_path}")

    da = rxr.open_rasterio(tif_path)

    if "band" in da.dims:
        da = da.squeeze("band", drop=True)

    if da.rio.crs != "EPSG:4326":
        da = da.rio.reproject("EPSG:4326")

    return da


def calculate_flood_volume(hmax: xr.DataArray, threshold: float = FLOOD_THRESHOLD) -> float:
    """Calculate total flood volume in m³."""
    try:
        dx = float(np.abs(hmax.x.diff("x").median()))
        dy = float(np.abs(hmax.y.diff("y").median()))
        dx_m = dx * 96000  # Approximate for Durban latitude (~30°S)
        dy_m = dy * 111000
        cell_area_m2 = dx_m * dy_m
    except Exception:
        cell_area_m2 = 100 * 100

    flood_mask = hmax > threshold
    volume_m3 = float(hmax.where(flood_mask).sum(skipna=True) * cell_area_m2)

    return volume_m3


def calculate_flood_extent(hmax: xr.DataArray, threshold: float = FLOOD_THRESHOLD) -> float:
    """Calculate total flood extent in km²."""
    try:
        dx = float(np.abs(hmax.x.diff("x").median()))
        dy = float(np.abs(hmax.y.diff("y").median()))
        dx_m = dx * 96000
        dy_m = dy * 111000
        cell_area_m2 = dx_m * dy_m
    except Exception:
        cell_area_m2 = 100 * 100

    flood_mask = hmax > threshold
    extent_m2 = float(flood_mask.sum() * cell_area_m2)
    extent_km2 = extent_m2 / 1e6

    return extent_km2


def calculate_mean_depth(hmax: xr.DataArray, threshold: float = FLOOD_THRESHOLD) -> float:
    """Calculate mean flood depth in flooded areas."""
    flood_mask = hmax > threshold
    return float(hmax.where(flood_mask).mean(skipna=True))


# ===== PLOTTING FUNCTIONS =====
def plot_hmax_comparison(hmax_factual: xr.DataArray, hmax_counterfactual: xr.DataArray, output_path: Path):
    """Create comparison plot showing factual, counterfactual, and climate-attributable flooding."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    flood_levels = np.linspace(0, HMAX_VMAX, 21)
    diff_levels = np.linspace(-DIFF_VMAX, DIFF_VMAX, 21)

    # Calculate difference: factual - counterfactual (positive = additional flooding from climate change)
    hmax_factual_aligned = hmax_factual.interp_like(hmax_counterfactual, method='nearest')
    diff = hmax_factual_aligned - hmax_counterfactual

    # Panel 1: Factual (current climate)
    hmax1_plot = hmax_factual.where(hmax_factual > FLOOD_THRESHOLD)
    im1 = hmax1_plot.plot(
        ax=axes[0], levels=flood_levels, cmap="viridis",
        add_colorbar=False, x="x", y="y", alpha=0.7, zorder=2
    )
    axes[0].set_title("Factual (Current Climate)", fontsize=12, fontweight="bold")

    # Panel 2: Counterfactual (-8% precipitation)
    hmax2_plot = hmax_counterfactual.where(hmax_counterfactual > FLOOD_THRESHOLD)
    im2 = hmax2_plot.plot(
        ax=axes[1], levels=flood_levels, cmap="viridis",
        add_colorbar=False, x="x", y="y", alpha=0.7, zorder=2
    )
    axes[1].set_title("Counterfactual (-8% Precip)", fontsize=12, fontweight="bold")

    # Panel 3: Difference (climate-attributable)
    im3 = diff.plot(
        ax=axes[2], levels=diff_levels, cmap="RdBu_r",
        add_colorbar=False, x="x", y="y", alpha=0.7, zorder=2
    )
    axes[2].set_title("Climate-Attributable Flooding", fontsize=12, fontweight="bold")

    # Format axes and add basemaps
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel("Longitude [°]")
        ax.set_ylabel("Latitude [°]")
        try:
            ctx.add_basemap(
                ax=ax, source=ctx.providers.OpenStreetMap.Mapnik,
                crs="EPSG:4326", attribution=False, zorder=1, zoom=11
            )
        except Exception:
            pass

    # Add colorbars
    cbar1 = fig.colorbar(im1, ax=axes[0], shrink=0.8, pad=0.02)
    cbar1.set_label("Max Depth [m]")
    cbar2 = fig.colorbar(im2, ax=axes[1], shrink=0.8, pad=0.02)
    cbar2.set_label("Max Depth [m]")
    cbar3 = fig.colorbar(im3, ax=axes[2], shrink=0.8, pad=0.02)
    cbar3.set_label("Difference [m]")

    # Statistics annotation
    mean_diff = float(diff.mean(skipna=True))
    max_diff = float(diff.max(skipna=True))

    plt.suptitle(
        f"Climate Attribution Analysis - {EVENT_NAME}\n"
        f"Mean additional depth: {mean_diff:.3f}m | Max additional depth: {max_diff:.3f}m",
        fontsize=13, fontweight="bold", y=1.02
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def plot_attribution_metrics(metrics: dict, output_path: Path):
    """Create bar chart comparing flood metrics with climate attribution percentage."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Calculate attribution percentages
    factual = metrics["Factual (CF0)"]
    counterfactual = metrics["Counterfactual (CF-8)"]

    vol_attributable = factual['volume'] - counterfactual['volume']
    ext_attributable = factual['extent'] - counterfactual['extent']
    depth_attributable = factual['mean_depth'] - counterfactual['mean_depth']

    vol_attr_pct = (vol_attributable / factual['volume']) * 100
    ext_attr_pct = (ext_attributable / factual['extent']) * 100
    depth_attr_pct = (depth_attributable / factual['mean_depth']) * 100

    # Colors: factual (current), counterfactual (pre-industrial), attributable (difference)
    colors = ['#FF6B6B', '#4ECDC4', '#FFE66D']  # Coral, Teal, Yellow
    labels = ['Factual\n(Current)', 'Counterfactual\n(-8% Precip)', 'Climate-\nAttributable']

    # Panel 1: Volume
    volumes = [factual['volume'] / 1e6, counterfactual['volume'] / 1e6, vol_attributable / 1e6]
    bars1 = axes[0].bar(range(3), volumes, color=colors, edgecolor='black', linewidth=1)
    axes[0].set_ylabel("Flood Volume [Mm³]", fontsize=11, fontweight="bold")
    axes[0].set_title(f"Flood Volume\n(Attribution: {vol_attr_pct:.1f}%)", fontsize=12, fontweight="bold")
    axes[0].set_xticks(range(3))
    axes[0].set_xticklabels(labels, fontsize=9)
    axes[0].set_ylim(0, max(volumes) * 1.2)
    axes[0].spines["top"].set_visible(False)
    axes[0].spines["right"].set_visible(False)
    axes[0].grid(axis='y', alpha=0.3)
    for bar, val in zip(bars1, volumes):
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(volumes)*0.02,
                     f'{val:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Panel 2: Extent
    extents = [factual['extent'], counterfactual['extent'], ext_attributable]
    bars2 = axes[1].bar(range(3), extents, color=colors, edgecolor='black', linewidth=1)
    axes[1].set_ylabel("Flood Extent [km²]", fontsize=11, fontweight="bold")
    axes[1].set_title(f"Flood Extent\n(Attribution: {ext_attr_pct:.1f}%)", fontsize=12, fontweight="bold")
    axes[1].set_xticks(range(3))
    axes[1].set_xticklabels(labels, fontsize=9)
    axes[1].set_ylim(0, max(extents) * 1.2)
    axes[1].spines["top"].set_visible(False)
    axes[1].spines["right"].set_visible(False)
    axes[1].grid(axis='y', alpha=0.3)
    for bar, val in zip(bars2, extents):
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(extents)*0.02,
                     f'{val:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Panel 3: Mean Depth
    depths = [factual['mean_depth'], counterfactual['mean_depth'], depth_attributable]
    bars3 = axes[2].bar(range(3), depths, color=colors, edgecolor='black', linewidth=1)
    axes[2].set_ylabel("Mean Depth [m]", fontsize=11, fontweight="bold")
    axes[2].set_title(f"Mean Flood Depth\n(Attribution: {depth_attr_pct:.1f}%)", fontsize=12, fontweight="bold")
    axes[2].set_xticks(range(3))
    axes[2].set_xticklabels(labels, fontsize=9)
    axes[2].set_ylim(0, max(depths) * 1.2)
    axes[2].spines["top"].set_visible(False)
    axes[2].spines["right"].set_visible(False)
    axes[2].grid(axis='y', alpha=0.3)
    for bar, val in zip(bars3, depths):
        axes[2].text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(depths)*0.02,
                     f'{val:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

    plt.suptitle(
        f"Climate Change Attribution - {EVENT_NAME}\n"
        f"Comparing current climate vs pre-industrial proxy (-8% precipitation)",
        fontsize=14, fontweight="bold", y=1.02
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def export_summary(metrics: dict, output_path: Path):
    """Export summary statistics to CSV with attribution analysis."""
    factual = metrics["Factual (CF0)"]
    counterfactual = metrics["Counterfactual (CF-8)"]

    # Calculate attributable amounts and percentages
    vol_attributable = factual['volume'] - counterfactual['volume']
    ext_attributable = factual['extent'] - counterfactual['extent']
    depth_attributable = factual['mean_depth'] - counterfactual['mean_depth']

    vol_attr_pct = (vol_attributable / factual['volume']) * 100
    ext_attr_pct = (ext_attributable / factual['extent']) * 100
    depth_attr_pct = (depth_attributable / factual['mean_depth']) * 100

    data = [
        {
            "Scenario": "Factual (CF0 - Current Climate)",
            "Flood_Volume_Mm3": factual['volume'] / 1e6,
            "Flood_Extent_km2": factual['extent'],
            "Mean_Depth_m": factual['mean_depth'],
        },
        {
            "Scenario": "Counterfactual (CF-8 - Pre-industrial proxy)",
            "Flood_Volume_Mm3": counterfactual['volume'] / 1e6,
            "Flood_Extent_km2": counterfactual['extent'],
            "Mean_Depth_m": counterfactual['mean_depth'],
        },
        {
            "Scenario": "Climate-Attributable (Factual - Counterfactual)",
            "Flood_Volume_Mm3": vol_attributable / 1e6,
            "Flood_Extent_km2": ext_attributable,
            "Mean_Depth_m": depth_attributable,
        },
        {
            "Scenario": "Attribution Percentage (%)",
            "Flood_Volume_Mm3": vol_attr_pct,
            "Flood_Extent_km2": ext_attr_pct,
            "Mean_Depth_m": depth_attr_pct,
        },
    ]

    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")

    # Print summary
    print("\n" + "="*80)
    print("CLIMATE ATTRIBUTION ANALYSIS SUMMARY")
    print(f"Event: {EVENT_NAME}")
    print("="*80)
    print(f"\n{'Metric':<25} {'Factual':>15} {'Counterfactual':>18} {'Attributable':>15} {'%':>8}")
    print("-"*80)
    print(f"{'Flood Volume [Mm³]':<25} {factual['volume']/1e6:>15.2f} "
          f"{counterfactual['volume']/1e6:>18.2f} {vol_attributable/1e6:>15.2f} {vol_attr_pct:>7.1f}%")
    print(f"{'Flood Extent [km²]':<25} {factual['extent']:>15.2f} "
          f"{counterfactual['extent']:>18.2f} {ext_attributable:>15.2f} {ext_attr_pct:>7.1f}%")
    print(f"{'Mean Depth [m]':<25} {factual['mean_depth']:>15.2f} "
          f"{counterfactual['mean_depth']:>18.2f} {depth_attributable:>15.2f} {depth_attr_pct:>7.1f}%")
    print("="*80)
    print(f"\nInterpretation: Approximately {vol_attr_pct:.1f}% of the flood volume is")
    print(f"attributable to climate change (the +8% precipitation from warming).")
    print("="*80 + "\n")


# ===== MAIN EXECUTION =====
def main():
    """Main execution function."""
    print(f"\n{'='*60}")
    print(f"Climate Attribution Analysis - {EVENT_NAME}")
    print(f"{'='*60}\n")

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    print("Loading data...")
    hmax_data = {}
    metrics = {}

    for name in RUNS.keys():
        try:
            print(f"  Loading {name}...")
            hmax = load_hmax_tif(name)
            hmax_data[name] = hmax

            # Calculate metrics
            metrics[name] = {
                'volume': calculate_flood_volume(hmax),
                'extent': calculate_flood_extent(hmax),
                'mean_depth': calculate_mean_depth(hmax),
            }

            print(f"    Volume: {metrics[name]['volume']/1e6:.2f} Mm³, "
                  f"Extent: {metrics[name]['extent']:.2f} km², "
                  f"Mean Depth: {metrics[name]['mean_depth']:.2f} m")

        except Exception as e:
            print(f"  ERROR loading {name}: {e}")
            print(f"  Make sure the SFINCS run has completed and post-processing is done.")
            return

    # Generate plots
    print("\nGenerating plots...")

    # 1. Hmax comparison
    print("  Creating hmax comparison plot...")
    plot_hmax_comparison(
        hmax_data["Factual (CF0)"],
        hmax_data["Counterfactual (CF-8)"],
        OUTPUT_DIR / f"hmax_climate_attribution_{EVENT_NAME.lower()}.png"
    )

    # 2. Attribution metrics bar chart
    print("  Creating attribution metrics chart...")
    plot_attribution_metrics(
        metrics,
        OUTPUT_DIR / f"metrics_climate_attribution_{EVENT_NAME.lower()}.png"
    )

    # 3. Export summary
    print("  Exporting summary CSV...")
    export_summary(
        metrics,
        OUTPUT_DIR / f"summary_climate_attribution_{EVENT_NAME.lower()}.csv"
    )

    print(f"\n{'='*60}")
    print(f"Analysis complete! Output files saved to: {OUTPUT_DIR}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()

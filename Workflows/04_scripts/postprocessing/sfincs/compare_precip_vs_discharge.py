"""
Compare SFINCS flood model outputs: Precipitation-only vs Precipitation + Discharge.

This script compares flood outputs between:
- Precipitation-only forcing (ERA5)
- Precipitation + River Discharge forcing (ERA5 + GloFAS)

Produces:
1. Bilateral hmax difference map showing where discharge adds flooding
2. Summary statistics with percentage change in flood volume and extent
3. Bar chart comparing total flood metrics

Author: Generated for Durban 2022 compound flooding analysis
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

# Base paths
BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs")
OUTPUT_DIR = Path("/p/11210471-001-compass/04_Results/precip_vs_discharge_comparison")

# Run configurations
RUNS = {
    "Precip Only (ERA5)": {
        "runname": "Durban2022",
        "folder": "event_precip_era5_hourly_CF0_no_wind",
    },
    "Precip + Discharge (ERA5 + GloFAS)": {
        "runname": "Durban2022_dis",
        "folder": "event_precip_era5_hourly_CF0_no_wind_dis",
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
    return BASE_RUN_PATH / REGION / config["runname"] / "sfincs" / config["folder"]


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


def load_discharge_forcing(run_name: str) -> pd.DataFrame:
    """Load the discharge forcing file if it exists."""
    run_path = get_run_path(run_name)
    dis_path = run_path / "sfincs.dis"

    if not dis_path.exists():
        return None

    # Read discharge file (space-separated, no header)
    df = pd.read_csv(dis_path, sep=r'\s+', header=None)
    df.columns = ['time_sec'] + [f'src_{i+1}' for i in range(len(df.columns) - 1)]
    df['time_hours'] = df['time_sec'] / 3600

    return df


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


def calculate_max_depth(hmax: xr.DataArray) -> float:
    """Calculate maximum flood depth."""
    return float(hmax.max(skipna=True))


def calculate_mean_depth(hmax: xr.DataArray, threshold: float = FLOOD_THRESHOLD) -> float:
    """Calculate mean flood depth in flooded areas."""
    flood_mask = hmax > threshold
    return float(hmax.where(flood_mask).mean(skipna=True))


# ===== PLOTTING FUNCTIONS =====
def plot_hmax_comparison(hmax_precip: xr.DataArray, hmax_discharge: xr.DataArray, output_path: Path):
    """Create comparison plot showing both runs and their difference."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    flood_levels = np.linspace(0, HMAX_VMAX, 21)
    diff_levels = np.linspace(-DIFF_VMAX, DIFF_VMAX, 21)

    # Calculate difference: discharge - precip (positive = more flooding with discharge)
    # Align grids first
    hmax_precip_aligned = hmax_precip.interp_like(hmax_discharge, method='nearest')
    diff = hmax_discharge - hmax_precip_aligned

    # Panel 1: Precip only
    hmax1_plot = hmax_precip.where(hmax_precip > FLOOD_THRESHOLD)
    im1 = hmax1_plot.plot(
        ax=axes[0], levels=flood_levels, cmap="viridis",
        add_colorbar=False, x="x", y="y", alpha=0.7, zorder=2
    )
    axes[0].set_title("Precipitation Only (ERA5)", fontsize=12, fontweight="bold")

    # Panel 2: Precip + Discharge
    hmax2_plot = hmax_discharge.where(hmax_discharge > FLOOD_THRESHOLD)
    im2 = hmax2_plot.plot(
        ax=axes[1], levels=flood_levels, cmap="viridis",
        add_colorbar=False, x="x", y="y", alpha=0.7, zorder=2
    )
    axes[1].set_title("Precip + Discharge (ERA5 + GloFAS)", fontsize=12, fontweight="bold")

    # Panel 3: Difference
    im3 = diff.plot(
        ax=axes[2], levels=diff_levels, cmap="RdBu_r",
        add_colorbar=False, x="x", y="y", alpha=0.7, zorder=2
    )
    axes[2].set_title("Difference (Discharge - Precip Only)", fontsize=12, fontweight="bold")

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
    max_increase = float(diff.max(skipna=True))

    plt.suptitle(
        f"Impact of River Discharge on Flood Depth - {EVENT_NAME}\n"
        f"Mean diff: {mean_diff:.3f}m | Max increase: {max_increase:.3f}m",
        fontsize=13, fontweight="bold", y=1.02
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def plot_discharge_timeseries(discharge_df: pd.DataFrame, output_path: Path):
    """Plot discharge time series from GloFAS."""
    if discharge_df is None:
        print("No discharge data available")
        return

    fig, ax = plt.subplots(figsize=(10, 5))

    colors = plt.cm.Set1(np.linspace(0, 1, len(discharge_df.columns) - 2))

    for i, col in enumerate([c for c in discharge_df.columns if c.startswith('src_')]):
        values = discharge_df[col].values
        if values.max() > 0:  # Only plot non-zero sources
            ax.plot(discharge_df['time_hours'], values, label=col, linewidth=2, color=colors[i])

    ax.set_xlabel("Time [hours from start]", fontsize=11)
    ax.set_ylabel("Discharge [m³/s]", fontsize=11)
    ax.set_title(f"GloFAS River Discharge Forcing - {EVENT_NAME}", fontsize=12, fontweight="bold")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, discharge_df['time_hours'].max())

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def plot_metrics_comparison(metrics: dict, output_path: Path):
    """Create bar chart comparing flood metrics with percentage change."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    run_names = list(metrics.keys())
    x = np.arange(len(run_names))
    width = 0.5

    # Colors
    colors = ['#4ECDC4', '#FF6B6B']  # Teal for precip-only, coral for with discharge

    # Calculate percentage changes
    precip_key = "Precip Only (ERA5)"
    discharge_key = "Precip + Discharge (ERA5 + GloFAS)"

    vol_change = ((metrics[discharge_key]['volume'] - metrics[precip_key]['volume']) /
                  metrics[precip_key]['volume'] * 100)
    ext_change = ((metrics[discharge_key]['extent'] - metrics[precip_key]['extent']) /
                  metrics[precip_key]['extent'] * 100)
    mean_change = ((metrics[discharge_key]['mean_depth'] - metrics[precip_key]['mean_depth']) /
                   metrics[precip_key]['mean_depth'] * 100)

    # Panel 1: Volume
    volumes = [metrics[n]['volume'] / 1e6 for n in run_names]
    bars1 = axes[0].bar(x, volumes, width, color=colors, edgecolor='black', linewidth=1)
    axes[0].set_ylabel("Flood Volume [Mm³]", fontsize=11, fontweight="bold")
    axes[0].set_title(f"Total Flood Volume\n({vol_change:+.1f}% change)", fontsize=12, fontweight="bold")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(['Precip\nOnly', 'Precip +\nDischarge'], fontsize=10)
    axes[0].set_ylim(0, max(volumes) * 1.2)
    axes[0].spines["top"].set_visible(False)
    axes[0].spines["right"].set_visible(False)
    axes[0].grid(axis='y', alpha=0.3)
    for bar, val in zip(bars1, volumes):
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(volumes)*0.02,
                     f'{val:.2f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Panel 2: Extent
    extents = [metrics[n]['extent'] for n in run_names]
    bars2 = axes[1].bar(x, extents, width, color=colors, edgecolor='black', linewidth=1)
    axes[1].set_ylabel("Flood Extent [km²]", fontsize=11, fontweight="bold")
    axes[1].set_title(f"Total Flood Extent\n({ext_change:+.1f}% change)", fontsize=12, fontweight="bold")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(['Precip\nOnly', 'Precip +\nDischarge'], fontsize=10)
    axes[1].set_ylim(0, max(extents) * 1.2)
    axes[1].spines["top"].set_visible(False)
    axes[1].spines["right"].set_visible(False)
    axes[1].grid(axis='y', alpha=0.3)
    for bar, val in zip(bars2, extents):
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(extents)*0.02,
                     f'{val:.2f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Panel 3: Mean depth
    mean_depths = [metrics[n]['mean_depth'] for n in run_names]
    bars3 = axes[2].bar(x, mean_depths, width, color=colors, edgecolor='black', linewidth=1)
    axes[2].set_ylabel("Mean Depth [m]", fontsize=11, fontweight="bold")
    axes[2].set_title(f"Mean Flood Depth\n({mean_change:+.1f}% change)", fontsize=12, fontweight="bold")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(['Precip\nOnly', 'Precip +\nDischarge'], fontsize=10)
    axes[2].set_ylim(0, max(mean_depths) * 1.2)
    axes[2].spines["top"].set_visible(False)
    axes[2].spines["right"].set_visible(False)
    axes[2].grid(axis='y', alpha=0.3)
    for bar, val in zip(bars3, mean_depths):
        axes[2].text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(mean_depths)*0.02,
                     f'{val:.2f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

    plt.suptitle(
        f"Impact of River Discharge on Flood Metrics - {EVENT_NAME}",
        fontsize=14, fontweight="bold", y=1.02
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def export_summary(metrics: dict, output_path: Path):
    """Export summary statistics to CSV."""
    precip_key = "Precip Only (ERA5)"
    discharge_key = "Precip + Discharge (ERA5 + GloFAS)"

    # Calculate percentage changes
    vol_change = ((metrics[discharge_key]['volume'] - metrics[precip_key]['volume']) /
                  metrics[precip_key]['volume'] * 100)
    ext_change = ((metrics[discharge_key]['extent'] - metrics[precip_key]['extent']) /
                  metrics[precip_key]['extent'] * 100)
    max_change = ((metrics[discharge_key]['max_depth'] - metrics[precip_key]['max_depth']) /
                  metrics[precip_key]['max_depth'] * 100)
    mean_change = ((metrics[discharge_key]['mean_depth'] - metrics[precip_key]['mean_depth']) /
                   metrics[precip_key]['mean_depth'] * 100)

    data = []
    for name, m in metrics.items():
        data.append({
            "Run": name,
            "Flood_Volume_Mm3": m['volume'] / 1e6,
            "Flood_Extent_km2": m['extent'],
            "Max_Depth_m": m['max_depth'],
            "Mean_Depth_m": m['mean_depth'],
        })

    # Add change row
    data.append({
        "Run": "% Change (Discharge vs Precip)",
        "Flood_Volume_Mm3": vol_change,
        "Flood_Extent_km2": ext_change,
        "Max_Depth_m": max_change,
        "Mean_Depth_m": mean_change,
    })

    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")

    # Print summary
    print("\n" + "="*70)
    print("SUMMARY: Impact of River Discharge on Flooding")
    print("="*70)
    print(f"\n{'Metric':<25} {'Precip Only':>15} {'Precip+Discharge':>18} {'Change':>12}")
    print("-"*70)
    print(f"{'Flood Volume [Mm³]':<25} {metrics[precip_key]['volume']/1e6:>15.2f} "
          f"{metrics[discharge_key]['volume']/1e6:>18.2f} {vol_change:>+11.1f}%")
    print(f"{'Flood Extent [km²]':<25} {metrics[precip_key]['extent']:>15.2f} "
          f"{metrics[discharge_key]['extent']:>18.2f} {ext_change:>+11.1f}%")
    print(f"{'Max Depth [m]':<25} {metrics[precip_key]['max_depth']:>15.2f} "
          f"{metrics[discharge_key]['max_depth']:>18.2f} {max_change:>+11.1f}%")
    print(f"{'Mean Depth [m]':<25} {metrics[precip_key]['mean_depth']:>15.2f} "
          f"{metrics[discharge_key]['mean_depth']:>18.2f} {mean_change:>+11.1f}%")
    print("="*70 + "\n")


# ===== MAIN EXECUTION =====
def main():
    """Main execution function."""
    print(f"\n{'='*60}")
    print(f"Precip vs Precip+Discharge Comparison - {EVENT_NAME}")
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
                'max_depth': calculate_max_depth(hmax),
                'mean_depth': calculate_mean_depth(hmax),
            }

            print(f"    Volume: {metrics[name]['volume']/1e6:.2f} Mm³, "
                  f"Extent: {metrics[name]['extent']:.2f} km², "
                  f"Max: {metrics[name]['max_depth']:.2f} m")

        except Exception as e:
            print(f"  ERROR loading {name}: {e}")
            return

    # Load discharge data
    discharge_df = load_discharge_forcing("Precip + Discharge (ERA5 + GloFAS)")

    # Generate plots
    print("\nGenerating plots...")

    # 1. Hmax comparison
    print("  Creating hmax comparison plot...")
    plot_hmax_comparison(
        hmax_data["Precip Only (ERA5)"],
        hmax_data["Precip + Discharge (ERA5 + GloFAS)"],
        OUTPUT_DIR / f"hmax_precip_vs_discharge_{EVENT_NAME.lower()}.png"
    )

    # 2. Discharge time series
    if discharge_df is not None:
        print("  Creating discharge time series plot...")
        plot_discharge_timeseries(
            discharge_df,
            OUTPUT_DIR / f"discharge_timeseries_{EVENT_NAME.lower()}.png"
        )

    # 3. Metrics comparison bar chart
    print("  Creating metrics comparison chart...")
    plot_metrics_comparison(
        metrics,
        OUTPUT_DIR / f"metrics_comparison_{EVENT_NAME.lower()}.png"
    )

    # 4. Export summary
    print("  Exporting summary CSV...")
    export_summary(
        metrics,
        OUTPUT_DIR / f"summary_precip_vs_discharge_{EVENT_NAME.lower()}.csv"
    )

    print(f"\n{'='*60}")
    print(f"Analysis complete! Output files saved to: {OUTPUT_DIR}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
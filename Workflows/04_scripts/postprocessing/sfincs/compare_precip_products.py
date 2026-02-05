"""
Compare SFINCS flood model outputs from different precipitation products.

This script compares flood outputs (hmax, precipitation forcing, total water volume)
across multiple precipitation forcing datasets (e.g., ERA5, GPM IMERG, MSWEP).

Produces:
1. Bilateral hmax difference maps between all pairs of precipitation products
2. Precipitation forcing time series comparison
3. Total water volume bar chart comparison

Author: Generated for Durban 2022 precipitation-only SFINCS workflow
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import rioxarray as rxr
import contextily as ctx
from pathlib import Path
from itertools import combinations
import warnings

warnings.filterwarnings("ignore")


# ===== CONFIGURATION =====
# Event settings
EVENT_NAME = "Durban_April2022"
REGION = "durban"
RUNNAME = "Durban2022"

# Base paths
BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs")
OUTPUT_DIR = Path("/p/11210471-001-compass/04_Results/precip_comparison")

# Precipitation products to compare
# Format: {display_name: folder_suffix}
PRECIP_PRODUCTS = {
    "ERA5": "era5_hourly",
    "GPM IMERG": "gpm_imerg_durban",
    "MSWEP": "mswep_v280_3h_durban_apr2022",
}

# Plotting thresholds
FLOOD_THRESHOLD = 0.05  # meters - minimum depth to consider as flooded
HMAX_VMIN = 0.0
HMAX_VMAX = 2.0
DIFF_VMIN = -0.5
DIFF_VMAX = 0.5


# ===== HELPER FUNCTIONS =====
def get_run_path(precip_name: str) -> Path:
    """Get the run directory path for a given precipitation product."""
    precip_suffix = PRECIP_PRODUCTS[precip_name]
    return BASE_RUN_PATH / REGION / RUNNAME / "sfincs" / f"event_precip_{precip_suffix}_CF0_no_wind"


def load_hmax_tif(precip_name: str) -> xr.DataArray:
    """Load the hmax TIFF file for a precipitation product."""
    run_path = get_run_path(precip_name)
    tif_path = run_path / "plot_output" / "sfincs_output_hmax_AllTime.tif"

    if not tif_path.exists():
        raise FileNotFoundError(f"TIFF file not found: {tif_path}")

    da = rxr.open_rasterio(tif_path)

    # Remove band dimension if present
    if "band" in da.dims:
        da = da.squeeze("band", drop=True)

    # Reproject to EPSG:4326 if needed
    if da.rio.crs != "EPSG:4326":
        da = da.rio.reproject("EPSG:4326")

    return da


def load_precip_forcing(precip_name: str) -> xr.DataArray:
    """Load the precipitation forcing data for a precipitation product."""
    run_path = get_run_path(precip_name)
    precip_path = run_path / "precip_2d.nc"

    if not precip_path.exists():
        raise FileNotFoundError(f"Precipitation file not found: {precip_path}")

    ds = xr.open_dataset(precip_path)
    # Get precipitation variable (usually 'Precipitation' or 'precip')
    precip_var = [v for v in ds.data_vars if 'precip' in v.lower()][0]
    return ds[precip_var]


def get_time_step_hours(precip: xr.DataArray) -> float:
    """Calculate time step in hours from precipitation data."""
    if len(precip.time) < 2:
        return 1.0

    # Use timedelta64 with minutes for better precision with sub-hourly data
    time_diff = np.diff(precip.time.values).astype('timedelta64[m]').astype(float) / 60.0
    dt_hours = np.median(time_diff)

    # Sanity check
    if dt_hours <= 0 or dt_hours > 24:
        dt_hours = 1.0

    return dt_hours


def load_sfincs_output(precip_name: str) -> xr.Dataset:
    """Load the SFINCS map output for volume calculations."""
    run_path = get_run_path(precip_name)
    map_path = run_path / "sfincs_map.nc"

    if not map_path.exists():
        raise FileNotFoundError(f"SFINCS map file not found: {map_path}")

    return xr.open_dataset(map_path)


def calculate_flood_volume(hmax: xr.DataArray, threshold: float = FLOOD_THRESHOLD) -> float:
    """Calculate total flood volume in m³ from hmax raster."""
    # Get cell area from coordinates
    try:
        dx = float(np.abs(hmax.x.diff("x").median()))
        dy = float(np.abs(hmax.y.diff("y").median()))
        # Convert degrees to approximate meters (at ~30°S latitude)
        # 1 degree longitude ≈ 96 km, 1 degree latitude ≈ 111 km
        dx_m = dx * 96000  # Approximate for Durban latitude
        dy_m = dy * 111000
        cell_area_m2 = dx_m * dy_m
    except Exception:
        # Fallback to rough estimate if coordinate calculation fails
        cell_area_m2 = 100 * 100  # 100m x 100m default

    # Calculate volume: sum of (depth * cell_area) for flooded cells
    flood_mask = hmax > threshold
    volume_m3 = float(hmax.where(flood_mask).sum(skipna=True) * cell_area_m2)

    return volume_m3


def calculate_flood_extent(hmax: xr.DataArray, threshold: float = FLOOD_THRESHOLD) -> float:
    """Calculate total flood extent in km² from hmax raster."""
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


# ===== PLOTTING FUNCTIONS =====
def plot_hmax_comparison(hmax_data: dict, output_path: Path):
    """
    Create a comparison plot showing hmax for all precipitation products.

    Parameters
    ----------
    hmax_data : dict
        Dictionary mapping product names to hmax DataArrays
    output_path : Path
        Output path for the figure
    """
    n_products = len(hmax_data)
    fig, axes = plt.subplots(1, n_products, figsize=(6 * n_products, 5))

    if n_products == 1:
        axes = [axes]

    flood_levels = np.linspace(HMAX_VMIN, HMAX_VMAX, 21)

    for ax, (name, hmax) in zip(axes, hmax_data.items()):
        # Mask low values for plotting
        hmax_plot = hmax.where(hmax > FLOOD_THRESHOLD)

        im = hmax_plot.plot(
            ax=ax,
            levels=flood_levels,
            cmap="viridis",
            add_colorbar=False,
            x="x", y="y",
            alpha=0.7,
            zorder=2,
        )
        ax.set_aspect("equal")
        ax.set_title(f"{name}", fontsize=12, fontweight="bold")
        ax.set_xlabel("Longitude [°]")
        ax.set_ylabel("Latitude [°]")

        # Add basemap
        try:
            ctx.add_basemap(
                ax=ax,
                source=ctx.providers.OpenStreetMap.Mapnik,
                crs="EPSG:4326",
                attribution=False,
                zorder=1,
                zoom=11,
            )
        except Exception as e:
            print(f"Could not add basemap: {e}")

    # Add shared colorbar
    cbar = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
    cbar.set_label("Max Flood Depth [m]", fontsize=11)

    plt.suptitle(
        f"Maximum Flood Depth Comparison - {EVENT_NAME}",
        fontsize=14, fontweight="bold", y=1.02
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def plot_bilateral_differences(hmax_data: dict, output_dir: Path):
    """
    Create bilateral difference plots between all pairs of precipitation products.

    Parameters
    ----------
    hmax_data : dict
        Dictionary mapping product names to hmax DataArrays
    output_dir : Path
        Output directory for figures
    """
    product_names = list(hmax_data.keys())
    pairs = list(combinations(product_names, 2))

    for name1, name2 in pairs:
        hmax1 = hmax_data[name1]
        hmax2 = hmax_data[name2]

        # Handle NaN values - set to 0 where other has data
        mask1_valid = ~np.isnan(hmax1)
        mask2_valid = ~np.isnan(hmax2)
        hmax1_filled = hmax1.where(mask1_valid | ~mask2_valid, 0)
        hmax2_filled = hmax2.where(mask2_valid | ~mask1_valid, 0)

        # Calculate difference: name1 - name2
        diff = hmax1_filled - hmax2_filled

        # Create figure with 3 panels
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        flood_levels = np.linspace(HMAX_VMIN, HMAX_VMAX, 21)
        diff_levels = np.linspace(DIFF_VMIN, DIFF_VMAX, 21)

        # Plot product 1
        hmax1_plot = hmax1.where(hmax1 > FLOOD_THRESHOLD)
        im1 = hmax1_plot.plot(
            ax=axes[0], levels=flood_levels, cmap="viridis",
            add_colorbar=False, x="x", y="y", alpha=0.7, zorder=2
        )
        axes[0].set_title(f"{name1}", fontsize=12, fontweight="bold")

        # Plot product 2
        hmax2_plot = hmax2.where(hmax2 > FLOOD_THRESHOLD)
        im2 = hmax2_plot.plot(
            ax=axes[1], levels=flood_levels, cmap="viridis",
            add_colorbar=False, x="x", y="y", alpha=0.7, zorder=2
        )
        axes[1].set_title(f"{name2}", fontsize=12, fontweight="bold")

        # Plot difference
        im3 = diff.plot(
            ax=axes[2], levels=diff_levels, cmap="RdBu_r",
            add_colorbar=False, x="x", y="y", alpha=0.7, zorder=2
        )
        axes[2].set_title(f"Difference ({name1} - {name2})", fontsize=12, fontweight="bold")

        # Add basemaps and format axes
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

        # Calculate statistics
        mean_diff = float(diff.mean(skipna=True))
        max_diff = float(diff.max(skipna=True))
        min_diff = float(diff.min(skipna=True))

        plt.suptitle(
            f"Flood Depth Comparison: {name1} vs {name2}\n"
            f"Mean diff: {mean_diff:.3f}m | Max diff: {max_diff:.3f}m | Min diff: {min_diff:.3f}m",
            fontsize=13, fontweight="bold", y=1.05
        )

        plt.tight_layout()
        filename = f"hmax_diff_{name1.lower().replace(' ', '_')}_vs_{name2.lower().replace(' ', '_')}.png"
        output_path = output_dir / filename
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Saved: {output_path}")


def plot_precip_forcing_comparison(precip_data: dict, output_path: Path):
    """
    Create a comparison plot of precipitation forcing time series.

    Parameters
    ----------
    precip_data : dict
        Dictionary mapping product names to precipitation DataArrays
    output_path : Path
        Output path for the figure
    """
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    colors = plt.cm.Set1(np.linspace(0, 1, len(precip_data)))

    # Top panel: Time series of spatial mean precipitation
    ax1 = axes[0]
    for (name, precip), color in zip(precip_data.items(), colors):
        # Calculate spatial mean precipitation rate over time
        precip_mean = precip.mean(dim=['x', 'y']) if 'x' in precip.dims else precip.mean(dim=['m', 'n'])
        ax1.plot(precip.time, precip_mean, label=name, color=color, linewidth=1.5)

    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mean Precipitation Rate [mm/hr]")
    ax1.set_title("Spatial Mean Precipitation Rate Over Time", fontsize=12, fontweight="bold")
    ax1.legend(loc="upper right")
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(axis='x', rotation=45)

    # Bottom panel: Cumulative precipitation
    ax2 = axes[1]
    for (name, precip), color in zip(precip_data.items(), colors):
        precip_mean = precip.mean(dim=['x', 'y']) if 'x' in precip.dims else precip.mean(dim=['m', 'n'])

        # Get time step in hours using helper function
        dt_hours = get_time_step_hours(precip)

        # Calculate cumulative (rate * time)
        cumsum = np.cumsum(precip_mean.values * dt_hours)
        ax2.plot(precip.time, cumsum, label=name, color=color, linewidth=2)

    ax2.set_xlabel("Time")
    ax2.set_ylabel("Cumulative Precipitation [mm]")
    ax2.set_title("Cumulative Precipitation Over Time", fontsize=12, fontweight="bold")
    ax2.legend(loc="upper left")
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='x', rotation=45)

    plt.suptitle(
        f"Precipitation Forcing Comparison - {EVENT_NAME}",
        fontsize=14, fontweight="bold", y=1.02
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def plot_volume_comparison(volumes: dict, extents: dict, output_path: Path):
    """
    Create bar charts comparing total flood volume and extent.

    Parameters
    ----------
    volumes : dict
        Dictionary mapping product names to flood volumes (m³)
    extents : dict
        Dictionary mapping product names to flood extents (km²)
    output_path : Path
        Output path for the figure
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    product_names = list(volumes.keys())
    x = np.arange(len(product_names))
    width = 0.6
    colors = plt.cm.Set2(np.linspace(0, 1, len(product_names)))

    # Volume bar chart
    ax1 = axes[0]
    volume_values = [volumes[name] / 1e6 for name in product_names]  # Convert to Mm³
    bars1 = ax1.bar(x, volume_values, width, color=colors, edgecolor='black', linewidth=1)

    ax1.set_ylabel("Flood Volume [Mm³]", fontsize=11, fontweight="bold")
    ax1.set_title("Total Flood Volume", fontsize=12, fontweight="bold")
    ax1.set_xticks(x)
    ax1.set_xticklabels(product_names, fontsize=10)
    ax1.set_ylim(0, max(volume_values) * 1.15)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.grid(axis='y', alpha=0.3)

    # Add value labels on bars
    for bar, val in zip(bars1, volume_values):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(volume_values)*0.02,
                 f'{val:.2f}', ha='center', va='bottom', fontsize=10)

    # Extent bar chart
    ax2 = axes[1]
    extent_values = [extents[name] for name in product_names]
    bars2 = ax2.bar(x, extent_values, width, color=colors, edgecolor='black', linewidth=1)

    ax2.set_ylabel("Flood Extent [km²]", fontsize=11, fontweight="bold")
    ax2.set_title("Total Flood Extent", fontsize=12, fontweight="bold")
    ax2.set_xticks(x)
    ax2.set_xticklabels(product_names, fontsize=10)
    ax2.set_ylim(0, max(extent_values) * 1.15)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.grid(axis='y', alpha=0.3)

    # Add value labels on bars
    for bar, val in zip(bars2, extent_values):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(extent_values)*0.02,
                 f'{val:.2f}', ha='center', va='bottom', fontsize=10)

    plt.suptitle(
        f"Flood Impact Comparison by Precipitation Product - {EVENT_NAME}",
        fontsize=14, fontweight="bold", y=1.02
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def export_summary_csv(volumes: dict, extents: dict, precip_totals: dict, output_path: Path):
    """Export summary statistics to CSV."""
    data = []
    for name in volumes.keys():
        data.append({
            "Precipitation_Product": name,
            "Flood_Volume_Mm3": volumes[name] / 1e6,
            "Flood_Extent_km2": extents[name],
            "Total_Precipitation_mm": precip_totals.get(name, np.nan),
        })

    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")
    print(df.to_string(index=False))


# ===== MAIN EXECUTION =====
def main():
    """Main execution function."""
    print(f"\n{'='*60}")
    print(f"Precipitation Product Comparison - {EVENT_NAME}")
    print(f"{'='*60}\n")

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load all data
    print("Loading data...")
    hmax_data = {}
    precip_data = {}
    volumes = {}
    extents = {}
    precip_totals = {}

    for name in PRECIP_PRODUCTS.keys():
        try:
            print(f"  Loading {name}...")

            # Load hmax
            hmax_data[name] = load_hmax_tif(name)

            # Load precipitation forcing
            precip_data[name] = load_precip_forcing(name)

            # Calculate volume and extent
            volumes[name] = calculate_flood_volume(hmax_data[name])
            extents[name] = calculate_flood_extent(hmax_data[name])

            # Calculate total precipitation
            precip = precip_data[name]
            precip_mean = precip.mean(dim=['x', 'y']) if 'x' in precip.dims else precip.mean(dim=['m', 'n'])
            dt_hours = get_time_step_hours(precip)
            precip_totals[name] = float(np.sum(precip_mean.values * dt_hours))

            print(f"    Volume: {volumes[name]/1e6:.2f} Mm³, Extent: {extents[name]:.2f} km², "
                  f"Total Precip: {precip_totals[name]:.1f} mm")

        except Exception as e:
            print(f"  ERROR loading {name}: {e}")
            continue

    if len(hmax_data) < 2:
        print("ERROR: Need at least 2 products to compare. Exiting.")
        return

    # Generate plots
    print("\nGenerating plots...")

    # 1. Overall hmax comparison
    print("  Creating hmax comparison plot...")
    plot_hmax_comparison(hmax_data, OUTPUT_DIR / f"hmax_comparison_{EVENT_NAME.lower()}.png")

    # 2. Bilateral differences
    print("  Creating bilateral difference plots...")
    plot_bilateral_differences(hmax_data, OUTPUT_DIR)

    # 3. Precipitation forcing comparison
    print("  Creating precipitation forcing comparison...")
    plot_precip_forcing_comparison(precip_data, OUTPUT_DIR / f"precip_forcing_comparison_{EVENT_NAME.lower()}.png")

    # 4. Volume and extent bar charts
    print("  Creating volume/extent bar charts...")
    plot_volume_comparison(volumes, extents, OUTPUT_DIR / f"volume_extent_comparison_{EVENT_NAME.lower()}.png")

    # 5. Export summary CSV
    print("  Exporting summary CSV...")
    export_summary_csv(volumes, extents, precip_totals, OUTPUT_DIR / f"summary_{EVENT_NAME.lower()}.csv")

    print(f"\n{'='*60}")
    print(f"Analysis complete! Output files saved to: {OUTPUT_DIR}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()

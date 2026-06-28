"""
Compare Climate DT Factual (Historical) vs Counterfactual (Control) scenarios.

This script compares flood outputs and precipitation forcing between:
- Factual (hist): Climate DT historical scenario (current climate)
- Counterfactual (cont): Climate DT control scenario (pre-industrial proxy)

Produces:
1. Flood depth (hmax) maps: Factual | Counterfactual | Difference
2. Precipitation time series comparison (rate and cumulative)
3. Precipitation spatial snapshots at multiple timesteps
4. Precipitation difference maps
5. Flood volume/extent bar charts with attribution
6. Summary statistics CSV

Author: Generated for Durban 2022 Climate DT climate attribution analysis
"""

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
EVENT_NAME = "Durban_April2022_ClimateDT"
REGION = "durban"

# Base paths
BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs")
OUTPUT_DIR = Path("/p/11210471-001-compass/04_Results/climateDT_attribution")

# Scenario configuration
# Factual = Historical (current climate), Counterfactual = Control (pre-industrial proxy)
SCENARIOS = {
    "Factual (hist)": {
        "runname": "Durban2022_ClimateDT",
        "suffix": "climateDT_tp_hist_Durban",
    },
    "Counterfactual (cont)": {
        "runname": "Durban2022_ClimateDT_cont",
        "suffix": "climateDT_tp_cont_Durban",
    },
}

# Plotting thresholds
FLOOD_THRESHOLD = 0.05  # meters - minimum depth to consider as flooded
HMAX_VMIN = 0.0
HMAX_VMAX = 2.0
DIFF_VMIN = -0.5
DIFF_VMAX = 0.5
PRECIP_VMAX = 50  # mm/hr for precipitation plots


# ===== HELPER FUNCTIONS =====
def get_run_path(scenario_name: str) -> Path:
    """Get the run directory path for a given scenario."""
    scenario_config = SCENARIOS[scenario_name]
    return (
        BASE_RUN_PATH
        / REGION
        / scenario_config["runname"]
        / "sfincs"
        / f"event_precip_{scenario_config['suffix']}_CF0_no_wind"
    )


def load_hmax_tif(scenario_name: str) -> xr.DataArray:
    """Load the hmax TIFF file for a scenario."""
    run_path = get_run_path(scenario_name)
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


def load_precip_forcing(scenario_name: str) -> xr.DataArray:
    """Load the precipitation forcing data for a scenario."""
    run_path = get_run_path(scenario_name)
    precip_path = run_path / "precip_2d.nc"

    if not precip_path.exists():
        raise FileNotFoundError(f"Precipitation file not found: {precip_path}")

    ds = xr.open_dataset(precip_path)
    # Get precipitation variable (usually 'Precipitation' or 'precip')
    precip_var = [v for v in ds.data_vars if "precip" in v.lower()][0]
    return ds[precip_var]


def get_time_step_hours(precip: xr.DataArray) -> float:
    """Calculate time step in hours from precipitation data."""
    if len(precip.time) < 2:
        return 1.0

    time_diff = (
        np.diff(precip.time.values).astype("timedelta64[m]").astype(float) / 60.0
    )
    dt_hours = np.median(time_diff)

    if dt_hours <= 0 or dt_hours > 24:
        dt_hours = 1.0

    return dt_hours


def calculate_flood_volume(
    hmax: xr.DataArray, threshold: float = FLOOD_THRESHOLD
) -> float:
    """Calculate total flood volume in m³ from hmax raster."""
    try:
        dx = float(np.abs(hmax.x.diff("x").median()))
        dy = float(np.abs(hmax.y.diff("y").median()))
        # Convert degrees to approximate meters (at ~30°S latitude)
        dx_m = dx * 96000
        dy_m = dy * 111000
        cell_area_m2 = dx_m * dy_m
    except Exception:
        cell_area_m2 = 100 * 100

    flood_mask = hmax > threshold
    volume_m3 = float(hmax.where(flood_mask).sum(skipna=True) * cell_area_m2)

    return volume_m3


def calculate_flood_extent(
    hmax: xr.DataArray, threshold: float = FLOOD_THRESHOLD
) -> float:
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
    Create a 3-panel comparison: Factual | Counterfactual | Difference.
    """
    names = list(hmax_data.keys())
    hmax_fact = hmax_data[names[0]]
    hmax_cf = hmax_data[names[1]]

    # Handle NaN values for difference
    mask_fact = ~np.isnan(hmax_fact)
    mask_cf = ~np.isnan(hmax_cf)
    hmax_fact_filled = hmax_fact.where(mask_fact | ~mask_cf, 0)
    hmax_cf_filled = hmax_cf.where(mask_cf | ~mask_fact, 0)
    diff = hmax_fact_filled - hmax_cf_filled

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    flood_levels = np.linspace(HMAX_VMIN, HMAX_VMAX, 21)
    diff_levels = np.linspace(DIFF_VMIN, DIFF_VMAX, 21)

    # Plot Factual
    hmax_fact_plot = hmax_fact.where(hmax_fact > FLOOD_THRESHOLD)
    im1 = hmax_fact_plot.plot(
        ax=axes[0],
        levels=flood_levels,
        cmap="viridis",
        add_colorbar=False,
        x="x",
        y="y",
        alpha=0.8,
        zorder=2,
    )
    axes[0].set_title(f"{names[0]}", fontsize=12, fontweight="bold")

    # Plot Counterfactual
    hmax_cf_plot = hmax_cf.where(hmax_cf > FLOOD_THRESHOLD)
    im2 = hmax_cf_plot.plot(
        ax=axes[1],
        levels=flood_levels,
        cmap="viridis",
        add_colorbar=False,
        x="x",
        y="y",
        alpha=0.8,
        zorder=2,
    )
    axes[1].set_title(f"{names[1]}", fontsize=12, fontweight="bold")

    # Plot Difference (Factual - Counterfactual = Climate Change Attribution)
    im3 = diff.plot(
        ax=axes[2],
        levels=diff_levels,
        cmap="RdBu_r",
        add_colorbar=False,
        x="x",
        y="y",
        alpha=0.8,
        zorder=2,
    )
    axes[2].set_title("Difference (Factual - Counterfactual)\n= Climate Attribution", fontsize=11, fontweight="bold")

    # Add basemaps and format axes
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel("Longitude [°]")
        ax.set_ylabel("Latitude [°]")
        try:
            ctx.add_basemap(
                ax=ax,
                source=ctx.providers.OpenStreetMap.Mapnik,
                crs="EPSG:4326",
                attribution=False,
                zorder=1,
                zoom=11,
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

    plt.suptitle(
        f"Climate DT Flood Depth Comparison - {EVENT_NAME}\n"
        f"Mean attribution: {mean_diff:.3f}m | Max attribution: {max_diff:.3f}m",
        fontsize=13,
        fontweight="bold",
        y=1.05,
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def plot_precip_timeseries(precip_data: dict, output_path: Path):
    """
    Create precipitation time series comparison with rate and cumulative panels.
    """
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))

    colors = {"Factual (hist)": "#1f77b4", "Counterfactual (cont)": "#ff7f0e"}

    # Top panel: Time series of spatial mean precipitation rate
    ax1 = axes[0]
    for name, precip in precip_data.items():
        precip_mean = (
            precip.mean(dim=["x", "y"])
            if "x" in precip.dims
            else precip.mean(dim=["m", "n"])
        )
        ax1.plot(
            precip.time,
            precip_mean,
            label=name,
            color=colors.get(name, "gray"),
            linewidth=2,
        )

    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mean Precipitation Rate [mm/hr]")
    ax1.set_title(
        "Spatial Mean Precipitation Rate Over Time", fontsize=12, fontweight="bold"
    )
    ax1.legend(loc="upper right", fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(axis="x", rotation=45)

    # Fill between to show difference
    names = list(precip_data.keys())
    precip1 = precip_data[names[0]]
    precip2 = precip_data[names[1]]
    mean1 = (
        precip1.mean(dim=["x", "y"])
        if "x" in precip1.dims
        else precip1.mean(dim=["m", "n"])
    )
    mean2 = (
        precip2.mean(dim=["x", "y"])
        if "x" in precip2.dims
        else precip2.mean(dim=["m", "n"])
    )
    ax1.fill_between(
        precip1.time.values,
        mean1.values,
        mean2.values,
        alpha=0.3,
        color="green",
        label="Difference (attribution)",
    )
    ax1.legend(loc="upper right", fontsize=10)

    # Bottom panel: Cumulative precipitation
    ax2 = axes[1]
    cumsum_data = {}
    for name, precip in precip_data.items():
        precip_mean = (
            precip.mean(dim=["x", "y"])
            if "x" in precip.dims
            else precip.mean(dim=["m", "n"])
        )
        dt_hours = get_time_step_hours(precip)
        cumsum = np.cumsum(precip_mean.values * dt_hours)
        cumsum_data[name] = cumsum
        ax2.plot(
            precip.time,
            cumsum,
            label=name,
            color=colors.get(name, "gray"),
            linewidth=2.5,
        )

    ax2.set_xlabel("Time")
    ax2.set_ylabel("Cumulative Precipitation [mm]")
    ax2.set_title(
        "Cumulative Precipitation Over Time", fontsize=12, fontweight="bold"
    )
    ax2.legend(loc="upper left", fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis="x", rotation=45)

    # Add annotation for total difference
    total_diff = cumsum_data[names[0]][-1] - cumsum_data[names[1]][-1]
    pct_diff = (total_diff / cumsum_data[names[1]][-1]) * 100 if cumsum_data[names[1]][-1] > 0 else 0
    ax2.annotate(
        f"Total difference:\n{total_diff:.1f} mm ({pct_diff:.1f}%)",
        xy=(0.98, 0.5),
        xycoords="axes fraction",
        ha="right",
        va="center",
        fontsize=11,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", edgecolor="orange"),
    )

    plt.suptitle(
        f"Climate DT Precipitation Forcing Comparison - {EVENT_NAME}",
        fontsize=14,
        fontweight="bold",
        y=1.02,
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def plot_precip_snapshots(precip_data: dict, output_path: Path, n_snapshots: int = 4):
    """
    Create spatial snapshots of precipitation at multiple timesteps.
    Shows both scenarios side-by-side with difference.
    """
    names = list(precip_data.keys())
    precip_fact = precip_data[names[0]]
    precip_cf = precip_data[names[1]]

    # Find peak timestep and create snapshots around it
    precip_mean = (
        precip_fact.mean(dim=["x", "y"])
        if "x" in precip_fact.dims
        else precip_fact.mean(dim=["m", "n"])
    )
    peak_idx = int(precip_mean.argmax().values)

    # Select snapshot indices (before peak, at peak, after peak)
    n_times = len(precip_fact.time)
    snapshot_indices = np.linspace(0, n_times - 1, n_snapshots, dtype=int)

    # Make sure peak is included
    if peak_idx not in snapshot_indices:
        snapshot_indices = sorted(list(set(list(snapshot_indices) + [peak_idx])))[:n_snapshots]

    fig, axes = plt.subplots(n_snapshots, 3, figsize=(15, 4 * n_snapshots))

    for row, idx in enumerate(snapshot_indices):
        # Get data at this timestep
        fact_snap = precip_fact.isel(time=idx)
        cf_snap = precip_cf.isel(time=idx)
        diff_snap = fact_snap - cf_snap

        time_str = np.datetime_as_string(precip_fact.time.values[idx], unit="h")

        # Determine coordinate names
        x_coord = "x" if "x" in fact_snap.dims else "n"
        y_coord = "y" if "y" in fact_snap.dims else "m"

        # Plot Factual
        im1 = fact_snap.plot(
            ax=axes[row, 0],
            cmap="Blues",
            vmin=0,
            vmax=PRECIP_VMAX,
            add_colorbar=False,
            x=x_coord,
            y=y_coord,
        )
        axes[row, 0].set_title(f"{names[0]}\n{time_str}", fontsize=10)

        # Plot Counterfactual
        im2 = cf_snap.plot(
            ax=axes[row, 1],
            cmap="Blues",
            vmin=0,
            vmax=PRECIP_VMAX,
            add_colorbar=False,
            x=x_coord,
            y=y_coord,
        )
        axes[row, 1].set_title(f"{names[1]}\n{time_str}", fontsize=10)

        # Plot Difference
        diff_max = max(abs(float(diff_snap.min())), abs(float(diff_snap.max())), 5)
        im3 = diff_snap.plot(
            ax=axes[row, 2],
            cmap="RdBu_r",
            vmin=-diff_max,
            vmax=diff_max,
            add_colorbar=False,
            x=x_coord,
            y=y_coord,
        )
        mean_diff = float(diff_snap.mean())
        axes[row, 2].set_title(f"Difference (mean: {mean_diff:.2f} mm/hr)\n{time_str}", fontsize=10)

        # Format axes
        for ax in axes[row]:
            ax.set_aspect("equal")
            if row == n_snapshots - 1:
                ax.set_xlabel("Longitude" if x_coord == "x" else "X")
            else:
                ax.set_xlabel("")
            ax.set_ylabel("Latitude" if y_coord == "y" else "Y")

    # Add colorbars
    cbar1 = fig.colorbar(im1, ax=axes[:, 0], shrink=0.6, pad=0.02)
    cbar1.set_label("Precip Rate [mm/hr]")
    cbar2 = fig.colorbar(im2, ax=axes[:, 1], shrink=0.6, pad=0.02)
    cbar2.set_label("Precip Rate [mm/hr]")
    cbar3 = fig.colorbar(im3, ax=axes[:, 2], shrink=0.6, pad=0.02)
    cbar3.set_label("Difference [mm/hr]")

    plt.suptitle(
        f"Precipitation Spatial Snapshots - {EVENT_NAME}",
        fontsize=14,
        fontweight="bold",
        y=1.01,
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def plot_precip_difference_map(precip_data: dict, output_path: Path):
    """
    Create a map showing the cumulative precipitation difference.
    """
    names = list(precip_data.keys())
    precip_fact = precip_data[names[0]]
    precip_cf = precip_data[names[1]]

    dt_hours = get_time_step_hours(precip_fact)

    # Calculate cumulative precipitation (sum over time * dt)
    x_coord = "x" if "x" in precip_fact.dims else "n"
    y_coord = "y" if "y" in precip_fact.dims else "m"

    cumsum_fact = precip_fact.sum(dim="time") * dt_hours
    cumsum_cf = precip_cf.sum(dim="time") * dt_hours
    cumsum_diff = cumsum_fact - cumsum_cf

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Determine common vmax
    vmax_cum = max(float(cumsum_fact.max()), float(cumsum_cf.max()))

    # Plot Factual cumulative
    im1 = cumsum_fact.plot(
        ax=axes[0],
        cmap="Blues",
        vmin=0,
        vmax=vmax_cum,
        add_colorbar=False,
        x=x_coord,
        y=y_coord,
    )
    total_fact = float(cumsum_fact.mean())
    axes[0].set_title(f"{names[0]}\nMean: {total_fact:.1f} mm", fontsize=12, fontweight="bold")

    # Plot Counterfactual cumulative
    im2 = cumsum_cf.plot(
        ax=axes[1],
        cmap="Blues",
        vmin=0,
        vmax=vmax_cum,
        add_colorbar=False,
        x=x_coord,
        y=y_coord,
    )
    total_cf = float(cumsum_cf.mean())
    axes[1].set_title(f"{names[1]}\nMean: {total_cf:.1f} mm", fontsize=12, fontweight="bold")

    # Plot Difference
    diff_max = max(abs(float(cumsum_diff.min())), abs(float(cumsum_diff.max())))
    im3 = cumsum_diff.plot(
        ax=axes[2],
        cmap="RdBu_r",
        vmin=-diff_max,
        vmax=diff_max,
        add_colorbar=False,
        x=x_coord,
        y=y_coord,
    )
    mean_diff = float(cumsum_diff.mean())
    pct_diff = (mean_diff / total_cf) * 100 if total_cf > 0 else 0
    axes[2].set_title(
        f"Difference (Attribution)\nMean: {mean_diff:.1f} mm ({pct_diff:.1f}%)",
        fontsize=12,
        fontweight="bold",
    )

    # Format axes
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel("Longitude" if x_coord == "x" else "X")
        ax.set_ylabel("Latitude" if y_coord == "y" else "Y")

    # Add colorbars
    cbar1 = fig.colorbar(im1, ax=axes[0], shrink=0.8, pad=0.02)
    cbar1.set_label("Cumulative Precip [mm]")
    cbar2 = fig.colorbar(im2, ax=axes[1], shrink=0.8, pad=0.02)
    cbar2.set_label("Cumulative Precip [mm]")
    cbar3 = fig.colorbar(im3, ax=axes[2], shrink=0.8, pad=0.02)
    cbar3.set_label("Difference [mm]")

    plt.suptitle(
        f"Cumulative Precipitation Comparison - {EVENT_NAME}",
        fontsize=14,
        fontweight="bold",
        y=1.02,
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def plot_attribution_bars(
    volumes: dict,
    extents: dict,
    precip_totals: dict,
    output_path: Path,
):
    """
    Create bar charts showing flood volume, extent, and precipitation with attribution.
    """
    names = list(volumes.keys())
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    colors = ["#1f77b4", "#ff7f0e"]
    x = np.arange(len(names))
    width = 0.5

    # Volume bar chart
    ax1 = axes[0]
    volume_values = [volumes[name] / 1e6 for name in names]  # Convert to Mm³
    bars1 = ax1.bar(x, volume_values, width, color=colors, edgecolor="black", linewidth=1)

    vol_diff = volume_values[0] - volume_values[1]
    vol_pct = (vol_diff / volume_values[1]) * 100 if volume_values[1] > 0 else 0

    ax1.set_ylabel("Flood Volume [Mm³]", fontsize=11, fontweight="bold")
    ax1.set_title(f"Total Flood Volume\nAttribution: {vol_diff:.2f} Mm³ ({vol_pct:.1f}%)", fontsize=11, fontweight="bold")
    ax1.set_xticks(x)
    ax1.set_xticklabels(names, fontsize=9, rotation=15, ha="right")
    ax1.set_ylim(0, max(volume_values) * 1.2)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.grid(axis="y", alpha=0.3)

    for bar, val in zip(bars1, volume_values):
        ax1.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(volume_values) * 0.02,
            f"{val:.2f}",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    # Extent bar chart
    ax2 = axes[1]
    extent_values = [extents[name] for name in names]
    bars2 = ax2.bar(x, extent_values, width, color=colors, edgecolor="black", linewidth=1)

    ext_diff = extent_values[0] - extent_values[1]
    ext_pct = (ext_diff / extent_values[1]) * 100 if extent_values[1] > 0 else 0

    ax2.set_ylabel("Flood Extent [km²]", fontsize=11, fontweight="bold")
    ax2.set_title(f"Total Flood Extent\nAttribution: {ext_diff:.2f} km² ({ext_pct:.1f}%)", fontsize=11, fontweight="bold")
    ax2.set_xticks(x)
    ax2.set_xticklabels(names, fontsize=9, rotation=15, ha="right")
    ax2.set_ylim(0, max(extent_values) * 1.2)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.grid(axis="y", alpha=0.3)

    for bar, val in zip(bars2, extent_values):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(extent_values) * 0.02,
            f"{val:.2f}",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    # Precipitation bar chart
    ax3 = axes[2]
    precip_values = [precip_totals[name] for name in names]
    bars3 = ax3.bar(x, precip_values, width, color=colors, edgecolor="black", linewidth=1)

    precip_diff = precip_values[0] - precip_values[1]
    precip_pct = (precip_diff / precip_values[1]) * 100 if precip_values[1] > 0 else 0

    ax3.set_ylabel("Total Precipitation [mm]", fontsize=11, fontweight="bold")
    ax3.set_title(f"Total Precipitation\nAttribution: {precip_diff:.1f} mm ({precip_pct:.1f}%)", fontsize=11, fontweight="bold")
    ax3.set_xticks(x)
    ax3.set_xticklabels(names, fontsize=9, rotation=15, ha="right")
    ax3.set_ylim(0, max(precip_values) * 1.2)
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    ax3.grid(axis="y", alpha=0.3)

    for bar, val in zip(bars3, precip_values):
        ax3.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(precip_values) * 0.02,
            f"{val:.1f}",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    plt.suptitle(
        f"Climate Attribution Summary - {EVENT_NAME}",
        fontsize=14,
        fontweight="bold",
        y=1.02,
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def export_summary_csv(
    volumes: dict,
    extents: dict,
    precip_totals: dict,
    output_path: Path,
):
    """Export summary statistics to CSV with attribution calculations."""
    names = list(volumes.keys())

    data = []
    for name in names:
        data.append(
            {
                "Scenario": name,
                "Flood_Volume_Mm3": volumes[name] / 1e6,
                "Flood_Extent_km2": extents[name],
                "Total_Precipitation_mm": precip_totals[name],
            }
        )

    # Add attribution row
    vol_attr = (volumes[names[0]] - volumes[names[1]]) / 1e6
    ext_attr = extents[names[0]] - extents[names[1]]
    precip_attr = precip_totals[names[0]] - precip_totals[names[1]]

    vol_pct = (vol_attr / (volumes[names[1]] / 1e6)) * 100 if volumes[names[1]] > 0 else 0
    ext_pct = (ext_attr / extents[names[1]]) * 100 if extents[names[1]] > 0 else 0
    precip_pct = (precip_attr / precip_totals[names[1]]) * 100 if precip_totals[names[1]] > 0 else 0

    data.append(
        {
            "Scenario": "Attribution (Factual - Counterfactual)",
            "Flood_Volume_Mm3": vol_attr,
            "Flood_Extent_km2": ext_attr,
            "Total_Precipitation_mm": precip_attr,
        }
    )
    data.append(
        {
            "Scenario": "Attribution (%)",
            "Flood_Volume_Mm3": vol_pct,
            "Flood_Extent_km2": ext_pct,
            "Total_Precipitation_mm": precip_pct,
        }
    )

    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")
    print("\n" + df.to_string(index=False))


# ===== MAIN EXECUTION =====
def main():
    """Main execution function."""
    print(f"\n{'='*70}")
    print(f"Climate DT Factual vs Counterfactual Comparison - {EVENT_NAME}")
    print(f"{'='*70}\n")

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load all data
    print("Loading data...")
    hmax_data = {}
    precip_data = {}
    volumes = {}
    extents = {}
    precip_totals = {}

    for name in SCENARIOS.keys():
        try:
            print(f"  Loading {name}...")
            run_path = get_run_path(name)
            print(f"    Path: {run_path}")

            # Load hmax
            hmax_data[name] = load_hmax_tif(name)

            # Load precipitation forcing
            precip_data[name] = load_precip_forcing(name)

            # Calculate volume and extent
            volumes[name] = calculate_flood_volume(hmax_data[name])
            extents[name] = calculate_flood_extent(hmax_data[name])

            # Calculate total precipitation
            precip = precip_data[name]
            precip_mean = (
                precip.mean(dim=["x", "y"])
                if "x" in precip.dims
                else precip.mean(dim=["m", "n"])
            )
            dt_hours = get_time_step_hours(precip)
            precip_totals[name] = float(np.sum(precip_mean.values * dt_hours))

            print(
                f"    Volume: {volumes[name]/1e6:.2f} Mm³, Extent: {extents[name]:.2f} km², "
                f"Total Precip: {precip_totals[name]:.1f} mm"
            )

        except Exception as e:
            print(f"  ERROR loading {name}: {e}")
            import traceback
            traceback.print_exc()
            return

    # Generate plots
    print("\nGenerating plots...")

    # 1. Hmax comparison (3-panel)
    print("  Creating hmax comparison plot...")
    plot_hmax_comparison(
        hmax_data, OUTPUT_DIR / f"hmax_comparison_{EVENT_NAME.lower()}.png"
    )

    # 2. Precipitation time series
    print("  Creating precipitation time series plot...")
    plot_precip_timeseries(
        precip_data, OUTPUT_DIR / f"precip_timeseries_{EVENT_NAME.lower()}.png"
    )

    # 3. Precipitation spatial snapshots
    print("  Creating precipitation spatial snapshots...")
    plot_precip_snapshots(
        precip_data, OUTPUT_DIR / f"precip_snapshots_{EVENT_NAME.lower()}.png"
    )

    # 4. Cumulative precipitation difference map
    print("  Creating cumulative precipitation difference map...")
    plot_precip_difference_map(
        precip_data, OUTPUT_DIR / f"precip_cumulative_diff_{EVENT_NAME.lower()}.png"
    )

    # 5. Attribution bar charts
    print("  Creating attribution bar charts...")
    plot_attribution_bars(
        volumes,
        extents,
        precip_totals,
        OUTPUT_DIR / f"attribution_bars_{EVENT_NAME.lower()}.png",
    )

    # 6. Export summary CSV
    print("  Exporting summary CSV...")
    export_summary_csv(
        volumes,
        extents,
        precip_totals,
        OUTPUT_DIR / f"summary_{EVENT_NAME.lower()}.csv",
    )

    print(f"\n{'='*70}")
    print(f"Analysis complete! Output files saved to: {OUTPUT_DIR}")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
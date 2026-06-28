"""
ClimateDT Climate Attribution Analysis Script

This script compares factual (historical) and counterfactual (control) climate scenarios
from ECMWF ClimateDT data to assess climate change attribution for:
- Total Precipitation (tp)
- 2m Temperature (2t)
- Mean Sea Level Pressure (msl)

The analysis includes:
- Accumulated precipitation totals
- Mean, min, max values over the entire period
- Spatial statistics (domain-averaged, point extremes)
- Temporal evolution and peak analysis
- Percentage/absolute differences between scenarios
- Comprehensive visualizations

Usage:
    python climateDT_attribution_analysis.py --data_dir /path/to/data --output_dir /path/to/output

Author: Generated for COMPASS project
"""

# %%
import os
import argparse
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple
import warnings

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.gridspec import GridSpec
from scipy.ndimage import label as ndimage_label, gaussian_filter, shift as ndimage_shift

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("Warning: cartopy not available. Spatial plots will use plain matplotlib.")

warnings.filterwarnings("ignore")


# =============================================================================
# Configuration Classes
# =============================================================================


@dataclass
class AnalysisConfig:
    """Configuration for the climate attribution analysis."""

    data_dir: Path  # Directory containing the ClimateDT NetCDF files
    output_dir: Path  # Directory for saving outputs
    region_name: str = "Durban"  # Region name for labels
    date_range: str = "20220409_to_20220414"  # Date range in filename format
    scenarios: List[str] = field(default_factory=lambda: ["hist", "cont"])
    scenario_labels: Dict[str, str] = field(
        default_factory=lambda: {
            "hist": "Factual (Historical)",
            "cont": "Counterfactual (Pre-industrial)",
        }
    )
    variables: List[str] = field(default_factory=lambda: ["tp", "msl", "2t"])

    # Storm core detection parameters (used by compute_storm_core_analysis)
    core_gaussian_sigma: float = 1.5        # Gaussian smoothing in grid cells
    core_percentile_threshold: float = 75.0  # Percentile of active pixels for threshold
    core_abs_threshold_mm: float = 1.0      # Hard minimum mm below which pixels excluded
    core_quiet_threshold_mm: float = 0.5    # If domain P75 < this, timestep is "quiet"
    core_peak_window_hours: float = 24.0    # Rolling window for peak accumulation (Method 1)

    variable_info: Dict[str, Dict] = field(
        default_factory=lambda: {
            "tp": {
                "name": "Total Precipitation",
                "unit": "mm",
                "unit_conversion": 1000,  # m to mm
                "colormap": "Blues",
                "diff_colormap": "BrBG",
                "aggregation": "sum",  # How to aggregate over time
            },
            "msl": {
                "name": "Mean Sea Level Pressure",
                "unit": "hPa",
                "unit_conversion": 0.01,  # Pa to hPa
                "colormap": "RdYlBu_r",
                "diff_colormap": "RdBu_r",
                "aggregation": "mean",
            },
            "2t": {
                "name": "2m Temperature",
                "unit": "°C",
                "unit_conversion": 1,  # Keep in K, convert later
                "kelvin_offset": -273.15,  # K to °C
                "colormap": "RdYlBu_r",
                "diff_colormap": "RdBu_r",
                "aggregation": "mean",
            },
        }
    )


# =============================================================================
# Data Loading Functions
# =============================================================================


def get_file_path(config: AnalysisConfig, variable: str, scenario: str) -> Path:
    """Construct the file path for a given variable and scenario."""
    filename = (
        f"climateDT_{variable}_{scenario}_{config.region_name}_{config.date_range}.nc"
    )
    return config.data_dir / filename


def load_dataset(config: AnalysisConfig, variable: str, scenario: str) -> xr.Dataset:
    """Load a single NetCDF dataset."""
    filepath = get_file_path(config, variable, scenario)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    ds = xr.open_dataset(filepath)
    return ds


def load_all_scenarios(
    config: AnalysisConfig, variable: str
) -> Dict[str, xr.DataArray]:
    """Load all scenarios for a given variable and apply unit conversions."""
    data = {}
    var_info = config.variable_info[variable]

    for scenario in config.scenarios:
        ds = load_dataset(config, variable, scenario)
        da = ds[variable]

        # Apply unit conversions
        da = da * var_info["unit_conversion"]
        if "kelvin_offset" in var_info:
            da = da + var_info["kelvin_offset"]

        # Update units attribute
        da.attrs["units"] = var_info["unit"]
        da.attrs["scenario"] = scenario

        data[scenario] = da

    return data


# =============================================================================
# Statistical Analysis Functions
# =============================================================================


def compute_temporal_statistics(da: xr.DataArray, var_info: Dict) -> Dict:
    """Compute temporal statistics for a data array."""
    stats = {}

    # Domain-averaged time series
    domain_mean = da.mean(dim=["latitude", "longitude"])

    # Overall statistics
    if var_info["aggregation"] == "sum":
        # For precipitation: accumulated total
        stats["accumulated_total"] = float(da.sum(dim="time").mean())
        stats["max_accumulated_point"] = float(da.sum(dim="time").max())
        stats["accumulated_spatial"] = da.sum(dim="time")
    else:
        # For temperature/pressure: period mean
        stats["period_mean"] = float(da.mean())
        stats["mean_spatial"] = da.mean(dim="time")

    # Extreme values
    stats["max_value"] = float(da.max())
    stats["max_time"] = str(da.where(da == da.max(), drop=True).time.values[0])
    stats["min_value"] = float(da.min())
    stats["min_time"] = str(da.where(da == da.min(), drop=True).time.values[0])

    # Domain-averaged statistics
    stats["domain_mean_max"] = float(domain_mean.max())
    stats["domain_mean_min"] = float(domain_mean.min())
    stats["domain_mean_mean"] = float(domain_mean.mean())
    stats["domain_mean_std"] = float(domain_mean.std())

    # Time series for plotting
    stats["domain_mean_timeseries"] = domain_mean

    return stats


def compute_percentile_statistics(
    da: xr.DataArray, percentiles: List[int] = [10, 25, 50, 75, 90, 95, 99]
) -> Dict:
    """Compute percentile statistics for extreme value analysis."""
    flat_data = da.values.flatten()
    flat_data = flat_data[~np.isnan(flat_data)]

    return {f"p{p}": float(np.percentile(flat_data, p)) for p in percentiles}


def compute_comparison_statistics(
    data_factual: xr.DataArray, data_counterfactual: xr.DataArray, var_info: Dict
) -> Dict:
    """Compute comparison statistics between factual and counterfactual scenarios."""
    comparison = {}

    # Absolute difference
    diff = data_factual - data_counterfactual

    # Relative difference (percentage)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_diff = (diff / np.abs(data_counterfactual)) * 100
        rel_diff = xr.where(np.isinf(rel_diff), np.nan, rel_diff)

    if var_info["aggregation"] == "sum":
        # For precipitation
        acc_factual = data_factual.sum(dim="time")
        acc_counter = data_counterfactual.sum(dim="time")
        acc_diff = acc_factual - acc_counter
        acc_rel_diff = (acc_diff / acc_counter) * 100

        comparison["accumulated_diff_spatial"] = acc_diff
        comparison["accumulated_rel_diff_spatial"] = acc_rel_diff
        comparison["accumulated_diff_mean"] = float(acc_diff.mean())
        comparison["accumulated_rel_diff_mean"] = float(acc_rel_diff.mean())
        comparison["accumulated_diff_max"] = float(acc_diff.max())
        comparison["accumulated_rel_diff_max"] = float(acc_rel_diff.max())
    else:
        # For temperature/pressure - mean over time
        mean_factual = data_factual.mean(dim="time")
        mean_counter = data_counterfactual.mean(dim="time")
        mean_diff = mean_factual - mean_counter
        mean_rel_diff = (mean_diff / np.abs(mean_counter)) * 100

        comparison["mean_diff_spatial"] = mean_diff
        comparison["mean_rel_diff_spatial"] = mean_rel_diff
        comparison["mean_diff_mean"] = float(mean_diff.mean())
        comparison["mean_rel_diff_mean"] = float(mean_rel_diff.mean())
        comparison["mean_diff_max"] = float(mean_diff.max())
        comparison["mean_diff_min"] = float(mean_diff.min())

    # Extreme values comparison
    comparison["max_diff"] = float(data_factual.max() - data_counterfactual.max())
    comparison["min_diff"] = float(data_factual.min() - data_counterfactual.min())

    # Time-varying domain-mean difference
    domain_mean_diff = data_factual.mean(
        dim=["latitude", "longitude"]
    ) - data_counterfactual.mean(dim=["latitude", "longitude"])
    comparison["domain_mean_diff_timeseries"] = domain_mean_diff
    comparison["domain_mean_diff_max"] = float(domain_mean_diff.max())
    comparison["domain_mean_diff_min"] = float(domain_mean_diff.min())

    return comparison


# =============================================================================
# Storm Attribution Analysis Functions (Methods 1, 2, 3)
# =============================================================================


# ---------------------------------------------------------------------------
# Method 1 — Temporal Relaxation: Peak Accumulation Window
# ---------------------------------------------------------------------------

def find_peak_accumulation_window(
    da: xr.DataArray,
    window_hours: float = 24.0,
) -> Dict:
    """
    Find the contiguous time window of a given length that maximises
    domain-averaged accumulated precipitation.

    Parameters
    ----------
    da : xr.DataArray
        Precipitation DataArray (time, latitude, longitude) in mm.
    window_hours : float
        Duration of the rolling window in hours (default 24).

    Returns
    -------
    dict with keys:
        start_idx, end_idx   : int — slice indices into the time axis
        start_time, end_time : str — ISO timestamps of the window
        window_accum         : xr.DataArray (lat, lon) — accumulated precip
                               over the window
        domain_mean_mm       : float — domain-averaged total over the window
    """
    # Infer timestep in hours
    if len(da.time) < 2:
        dt_hours = 1.0
    else:
        dt_ns = np.diff(da.time.values).astype("timedelta64[s]").astype(float)
        dt_hours = float(np.median(dt_ns)) / 3600.0
    dt_hours = max(dt_hours, 1e-3)

    N = max(1, int(round(window_hours / dt_hours)))  # window in timesteps

    # Domain-mean time series → rolling sum to find peak window
    domain_mean_ts = da.mean(dim=["latitude", "longitude"])
    rolling_sum = domain_mean_ts.rolling(time=N, min_periods=N).sum()

    # argmax gives the LAST index of the best window (rolling is right-aligned)
    end_idx = int(rolling_sum.argmax(dim="time").values)
    start_idx = max(0, end_idx - N + 1)

    window_slice = da.isel(time=slice(start_idx, end_idx + 1))
    window_accum = window_slice.sum(dim="time")

    return {
        "start_idx": start_idx,
        "end_idx": end_idx,
        "start_time": str(da.time.values[start_idx])[:16],
        "end_time": str(da.time.values[end_idx])[:16],
        "window_accum": window_accum,
        "domain_mean_mm": float(window_accum.mean()),
    }


def compute_peak_window_comparison(
    data_factual: xr.DataArray,
    data_counterfactual: xr.DataArray,
    window_hours: float = 24.0,
) -> Dict:
    """
    Compare accumulated precipitation over each scenario's own peak window.

    Each scenario independently selects its best N-hour window, so a storm
    that peaks 6 hours later in the counterfactual is still measured at its
    own best moment (not penalised for being late).

    Parameters
    ----------
    data_factual, data_counterfactual : xr.DataArray
        Precipitation DataArrays (time, latitude, longitude) in mm.
    window_hours : float
        Window duration (hours). Configurable via AnalysisConfig.core_peak_window_hours.

    Returns
    -------
    dict with keys:
        hist_window      : dict from find_peak_accumulation_window()
        cont_window      : dict from find_peak_accumulation_window()
        diff_accum       : xr.DataArray (lat, lon) = factual - counterfactual window accum
        window_offset_hours : float — timing offset between the two peak windows
        scalar_stats     : dict
    """
    hist_win = find_peak_accumulation_window(data_factual, window_hours)
    cont_win = find_peak_accumulation_window(data_counterfactual, window_hours)

    diff_accum = hist_win["window_accum"] - cont_win["window_accum"]

    # Timing offset between the two window centres
    hist_centre_idx = (hist_win["start_idx"] + hist_win["end_idx"]) / 2.0
    cont_centre_idx = (cont_win["start_idx"] + cont_win["end_idx"]) / 2.0

    dt_ns = np.diff(data_factual.time.values).astype("timedelta64[s]").astype(float)
    dt_hours = float(np.median(dt_ns)) / 3600.0 if len(dt_ns) > 0 else 1.0
    window_offset_hours = (hist_centre_idx - cont_centre_idx) * dt_hours

    abs_diff = hist_win["domain_mean_mm"] - cont_win["domain_mean_mm"]
    rel_diff_pct = (
        abs_diff / cont_win["domain_mean_mm"] * 100
        if cont_win["domain_mean_mm"] > 0 else np.nan
    )

    return {
        "hist_window": hist_win,
        "cont_window": cont_win,
        "diff_accum": diff_accum,
        "window_offset_hours": window_offset_hours,
        "window_hours": window_hours,
        "scalar_stats": {
            "factual_window_accum_mm": hist_win["domain_mean_mm"],
            "counterfactual_window_accum_mm": cont_win["domain_mean_mm"],
            "absolute_diff_mm": abs_diff,
            "relative_diff_pct": rel_diff_pct,
            "window_offset_hours": window_offset_hours,
            "factual_window_start": hist_win["start_time"],
            "factual_window_end": hist_win["end_time"],
            "counterfactual_window_start": cont_win["start_time"],
            "counterfactual_window_end": cont_win["end_time"],
        },
    }


# ---------------------------------------------------------------------------
# Method 2 — Object-Based Decomposition: Intensity × Area × Volume
# ---------------------------------------------------------------------------

def compute_object_decomposition(
    data_factual: xr.DataArray,
    data_counterfactual: xr.DataArray,
    storm_core_results: Dict,
) -> Dict:
    """
    Decompose the volume attribution into Intensity and Area contributions.

    At each active (non-quiet) timestep:
    - Intensity = P95 of precipitation pixels inside the core mask
    - Area      = core footprint in km²  (from core_pixel_count × cell_area)

    Both are compared between factual and counterfactual to show *why* the
    volume changed: did the storm get wetter, grow larger, or both?

    Parameters
    ----------
    data_factual, data_counterfactual : xr.DataArray
        Precipitation DataArrays (time, latitude, longitude) in mm.
    storm_core_results : dict
        Returned by compute_storm_core_analysis().

    Returns
    -------
    dict with keys:
        time                 : np.ndarray of timestamps
        intensity_hist_ts    : np.ndarray — P95-in-core per timestep (factual)
        intensity_cont_ts    : np.ndarray — P95-in-core per timestep (counterfactual)
        area_hist_ts_km2     : np.ndarray — core area per timestep (factual)
        area_cont_ts_km2     : np.ndarray — core area per timestep (counterfactual)
        scalar_stats         : dict (event-aggregated metrics)
    """
    hist_mask = storm_core_results["factual_core_mask"]
    cont_mask = storm_core_results["counterfactual_core_mask"]
    hist_track = storm_core_results["factual_track_df"]
    cont_track = storm_core_results["counterfactual_track_df"]

    # Cell area for the domain (use mean lat for approximation)
    dlat = float(np.abs(np.diff(data_factual.latitude.values).mean()))
    dlon = float(np.abs(np.diff(data_factual.longitude.values).mean()))
    km_per_deg = 111.32
    mean_lat_rad = float(np.deg2rad(data_factual.latitude.values.mean()))
    cell_area_km2 = dlat * dlon * (km_per_deg ** 2) * np.cos(mean_lat_rad)

    n_time = len(data_factual.time)
    intensity_hist = np.full(n_time, np.nan)
    intensity_cont = np.full(n_time, np.nan)
    area_hist_km2 = np.full(n_time, np.nan)
    area_cont_km2 = np.full(n_time, np.nan)

    for i in range(n_time):
        h_mask = hist_mask.isel(time=i).values
        c_mask = cont_mask.isel(time=i).values
        h_field = data_factual.isel(time=i).values
        c_field = data_counterfactual.isel(time=i).values

        if h_mask.any():
            vals = h_field[h_mask]
            intensity_hist[i] = float(np.percentile(vals[~np.isnan(vals)], 95)) if len(vals) > 0 else np.nan
            area_hist_km2[i] = float(h_mask.sum()) * cell_area_km2

        if c_mask.any():
            vals = c_field[c_mask]
            intensity_cont[i] = float(np.percentile(vals[~np.isnan(vals)], 95)) if len(vals) > 0 else np.nan
            area_cont_km2[i] = float(c_mask.sum()) * cell_area_km2

    # Event-aggregated stats (mean over active timesteps only)
    mean_int_hist = float(np.nanmean(intensity_hist))
    mean_int_cont = float(np.nanmean(intensity_cont))
    int_diff = mean_int_hist - mean_int_cont
    int_rel_pct = (int_diff / mean_int_cont * 100) if mean_int_cont > 0 else np.nan

    mean_area_hist = float(np.nanmean(area_hist_km2))
    mean_area_cont = float(np.nanmean(area_cont_km2))
    area_diff = mean_area_hist - mean_area_cont
    area_rel_pct = (area_diff / mean_area_cont * 100) if mean_area_cont > 0 else np.nan

    sc = storm_core_results["scalar_stats"]

    return {
        "time": data_factual.time.values,
        "intensity_hist_ts": intensity_hist,
        "intensity_cont_ts": intensity_cont,
        "area_hist_ts_km2": area_hist_km2,
        "area_cont_ts_km2": area_cont_km2,
        "scalar_stats": {
            "mean_intensity_hist_mm": mean_int_hist,
            "mean_intensity_cont_mm": mean_int_cont,
            "intensity_diff_mm": int_diff,
            "intensity_rel_diff_pct": int_rel_pct,
            "mean_area_hist_km2": mean_area_hist,
            "mean_area_cont_km2": mean_area_cont,
            "area_diff_km2": area_diff,
            "area_rel_diff_pct": area_rel_pct,
            "volume_factual_mm_km2": sc["factual_total_core_volume_mm_km2"],
            "volume_cont_mm_km2": sc["counterfactual_total_core_volume_mm_km2"],
            "volume_rel_diff_pct": sc["relative_diff_pct"],
        },
    }


# ---------------------------------------------------------------------------
# Method 3 — Storm-Relative Compositing: Centroid-Aligned Difference Map
# ---------------------------------------------------------------------------

def _shift_field_by_latlon(
    field_2d: np.ndarray,
    lat: np.ndarray,
    lon: np.ndarray,
    delta_lat: float,
    delta_lon: float,
) -> np.ndarray:
    """
    Translate a 2D field by (delta_lat, delta_lon) degrees using a discrete
    pixel shift (scipy.ndimage.shift with spline interpolation, order=1).

    The shift in pixels is computed from the grid spacing, so fractional
    pixel offsets are handled by linear interpolation between grid points.
    Pixels that fall outside the domain after shifting are filled with 0.

    Parameters
    ----------
    field_2d : np.ndarray  (lat x lon)
    lat, lon : np.ndarray  — 1D coordinate arrays
    delta_lat, delta_lon : float  — shift to apply in degrees
        Positive delta_lat moves features northward (towards higher lat index
        if lat is ascending, towards lower index if descending).

    Returns
    -------
    shifted : np.ndarray (lat x lon)
    """
    dlat = float(np.median(np.diff(lat)))   # degrees per pixel (signed)
    dlon = float(np.median(np.diff(lon)))

    # pixels to shift: delta / spacing  (positive = shift array in that axis)
    shift_lat = delta_lat / dlat
    shift_lon = delta_lon / dlon

    # Fill NaN with 0 before shifting to avoid artefacts at boundary
    field_clean = np.where(np.isnan(field_2d), 0.0, field_2d)
    shifted = ndimage_shift(
        field_clean,
        shift=(shift_lat, shift_lon),
        order=1,          # bilinear — fast and clean
        mode="constant",
        cval=0.0,
    )
    return shifted


def _compute_radial_profile(
    field_2d: np.ndarray,
    lat: np.ndarray,
    lon: np.ndarray,
    centroid_lat: float,
    centroid_lon: float,
    n_bins: int = 30,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the azimuthally-averaged precipitation profile as a function
    of distance from the centroid (in km).

    Uses Haversine distance for accurate km conversion at any latitude.

    Returns
    -------
    r_km : np.ndarray (n_bins,) — bin centre distances in km
    mean_precip : np.ndarray (n_bins,) — mean precipitation per bin (mm)
    """
    R_EARTH = 6371.0  # km
    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)
    clat_rad = np.deg2rad(centroid_lat)
    clon_rad = np.deg2rad(centroid_lon)

    LAT_r, LON_r = np.meshgrid(lat_rad, lon_rad, indexing="ij")

    dlat = LAT_r - clat_rad
    dlon = LON_r - clon_rad
    a = np.sin(dlat / 2) ** 2 + np.cos(clat_rad) * np.cos(LAT_r) * np.sin(dlon / 2) ** 2
    dist_km = 2 * R_EARTH * np.arcsin(np.sqrt(np.clip(a, 0, 1)))

    # Bin precipitation by distance
    max_dist = dist_km.max()
    bin_edges = np.linspace(0, max_dist, n_bins + 1)
    r_km = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    mean_precip = np.full(n_bins, np.nan)
    flat_dist = dist_km.ravel()
    flat_field = field_2d.ravel()
    valid = ~np.isnan(flat_field)

    for k in range(n_bins):
        in_bin = valid & (flat_dist >= bin_edges[k]) & (flat_dist < bin_edges[k + 1])
        if in_bin.any():
            mean_precip[k] = float(flat_field[in_bin].mean())

    return r_km, mean_precip


def compute_storm_relative_composite(
    data_factual: xr.DataArray,
    data_counterfactual: xr.DataArray,
    storm_core_results: Dict,
    composite_halfwidth_hours: float = 3.0,
) -> Dict:
    """
    Compute a storm-relative difference map by shifting the counterfactual
    precipitation field so its centroid aligns with the factual centroid.

    This removes the spatial displacement signal and reveals the pure
    intensity/structure change as a (hopefully) concentric pattern around
    a shared centre, rather than a noisy red/blue dipole.

    Compositing is done over ±composite_halfwidth_hours around the factual
    volume peak for stability (avoids noise from a single timestep).

    Parameters
    ----------
    data_factual, data_counterfactual : xr.DataArray
        Precipitation DataArrays (time, latitude, longitude) in mm.
    storm_core_results : dict
        Returned by compute_storm_core_analysis().
    composite_halfwidth_hours : float
        Half-width of the compositing window in hours (default 3).

    Returns
    -------
    dict with keys:
        factual_composite      : np.ndarray (lat, lon) — factual avg over window
        cf_composite_raw       : np.ndarray (lat, lon) — CF avg before shifting
        cf_composite_shifted   : np.ndarray (lat, lon) — CF shifted to factual centroid
        storm_relative_diff    : np.ndarray (lat, lon) = factual - shifted CF
        lat, lon               : np.ndarray — coordinate arrays
        radial_profile_hist    : dict  r_km, precip_mean
        radial_profile_cont    : dict  r_km, precip_mean (from shifted CF)
        factual_centroid_lat, factual_centroid_lon   : float
        cf_centroid_lat, cf_centroid_lon             : float (before shift)
        displacement_lat_deg, displacement_lon_deg   : float
        displacement_km        : float
        peak_window_str        : str (e.g. "2022-04-11T06 ± 3h")
        scalar_stats           : dict
    """
    hist_vol_ts = storm_core_results["factual_core_volume_ts"]
    hist_track = storm_core_results["factual_track_df"]
    cont_track = storm_core_results["counterfactual_track_df"]

    lat = data_factual.latitude.values
    lon = data_factual.longitude.values

    # Infer dt in hours
    dt_ns = np.diff(data_factual.time.values).astype("timedelta64[s]").astype(float)
    dt_hours = float(np.median(dt_ns)) / 3600.0 if len(dt_ns) > 0 else 1.0

    half_n = max(0, int(round(composite_halfwidth_hours / dt_hours)))
    peak_idx = int(hist_vol_ts.argmax().values)
    win_start = max(0, peak_idx - half_n)
    win_end = min(len(data_factual.time) - 1, peak_idx + half_n)

    # Average factual and counterfactual over the window
    factual_composite = data_factual.isel(time=slice(win_start, win_end + 1)).mean(dim="time").values
    cf_composite_raw = data_counterfactual.isel(time=slice(win_start, win_end + 1)).mean(dim="time").values

    # Centroid: mean over window rows that are not quiet
    win_h_track = hist_track.iloc[win_start:win_end + 1]
    win_c_track = cont_track.iloc[win_start:win_end + 1]

    h_active = win_h_track[~win_h_track["is_quiet"]]
    c_active = win_c_track[~win_c_track["is_quiet"]]

    clat_f = float(h_active["centroid_lat"].mean()) if len(h_active) > 0 else float(hist_track["centroid_lat"].mean(skipna=True))
    clon_f = float(h_active["centroid_lon"].mean()) if len(h_active) > 0 else float(hist_track["centroid_lon"].mean(skipna=True))
    clat_c = float(c_active["centroid_lat"].mean()) if len(c_active) > 0 else float(cont_track["centroid_lat"].mean(skipna=True))
    clon_c = float(c_active["centroid_lon"].mean()) if len(c_active) > 0 else float(cont_track["centroid_lon"].mean(skipna=True))

    delta_lat = clat_f - clat_c
    delta_lon = clon_f - clon_c

    # Shift counterfactual composite to factual centroid
    cf_composite_shifted = _shift_field_by_latlon(cf_composite_raw, lat, lon, delta_lat, delta_lon)
    storm_relative_diff = factual_composite - cf_composite_shifted

    # Radial profiles (both from factual centroid)
    r_h, mp_h = _compute_radial_profile(factual_composite, lat, lon, clat_f, clon_f)
    r_c, mp_c = _compute_radial_profile(cf_composite_shifted, lat, lon, clat_f, clon_f)

    # Displacement magnitude (Haversine)
    R_EARTH = 6371.0
    dlat_r = np.deg2rad(delta_lat)
    dlon_r = np.deg2rad(delta_lon)
    clat_r = np.deg2rad(clat_c)
    a = np.sin(dlat_r / 2) ** 2 + np.cos(clat_r) * np.cos(clat_r + dlat_r) * np.sin(dlon_r / 2) ** 2
    displacement_km = float(2 * R_EARTH * np.arcsin(np.sqrt(np.clip(a, 0, 1))))

    peak_time_str = str(data_factual.time.values[peak_idx])[:16]

    return {
        "factual_composite": factual_composite,
        "cf_composite_raw": cf_composite_raw,
        "cf_composite_shifted": cf_composite_shifted,
        "storm_relative_diff": storm_relative_diff,
        "lat": lat,
        "lon": lon,
        "radial_profile_hist": {"r_km": r_h, "precip_mean": mp_h},
        "radial_profile_cont": {"r_km": r_c, "precip_mean": mp_c},
        "factual_centroid_lat": clat_f,
        "factual_centroid_lon": clon_f,
        "cf_centroid_lat": clat_c,
        "cf_centroid_lon": clon_c,
        "displacement_lat_deg": delta_lat,
        "displacement_lon_deg": delta_lon,
        "displacement_km": displacement_km,
        "peak_window_str": f"{peak_time_str} ± {composite_halfwidth_hours:.0f}h",
        "scalar_stats": {
            "displacement_km": displacement_km,
            "displacement_lat_deg": delta_lat,
            "displacement_lon_deg": delta_lon,
            "factual_composite_mean_mm": float(np.nanmean(factual_composite)),
            "cf_composite_mean_mm": float(np.nanmean(cf_composite_raw)),
            "storm_relative_diff_mean_mm": float(np.nanmean(storm_relative_diff)),
            "storm_relative_diff_max_mm": float(np.nanmax(storm_relative_diff)),
            "storm_relative_diff_min_mm": float(np.nanmin(storm_relative_diff)),
            "peak_window_str": f"{peak_time_str} ± {composite_halfwidth_hours:.0f}h",
        },
    }


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def compute_lagrangian_attribution(
    data_factual: xr.DataArray,
    data_counterfactual: xr.DataArray,
    storm_core_results: Dict,
    config: "AnalysisConfig",
) -> Dict:
    """
    Run all three storm attribution methods and bundle into lagrangian_results.

    Parts:
        1. peak_window    — temporal relaxation: each scenario's own best window
        2. object_decomp  — intensity (P95) × area × volume decomposition
        3. storm_relative — centroid-aligned difference map + radial profiles

    Parameters
    ----------
    data_factual, data_counterfactual : xr.DataArray
    storm_core_results : dict  — from compute_storm_core_analysis()
    config : AnalysisConfig

    Returns
    -------
    lagrangian_results : dict with keys "peak_window", "object_decomp",
                         "storm_relative"
    """
    window_hours = getattr(config, "core_peak_window_hours", 24.0)

    print("  Part 1: peak accumulation window comparison...")
    peak_window = compute_peak_window_comparison(
        data_factual, data_counterfactual, window_hours=window_hours
    )

    print("  Part 2: object-based intensity × area decomposition...")
    object_decomp = compute_object_decomposition(
        data_factual, data_counterfactual, storm_core_results
    )

    print("  Part 3: storm-relative centroid-aligned composite...")
    storm_relative = compute_storm_relative_composite(
        data_factual, data_counterfactual, storm_core_results
    )

    pw = peak_window["scalar_stats"]
    od = object_decomp["scalar_stats"]
    sr = storm_relative["scalar_stats"]
    print(
        f"  Peak window: F={pw['factual_window_accum_mm']:.1f} mm  "
        f"CF={pw['counterfactual_window_accum_mm']:.1f} mm  "
        f"Δ={pw['relative_diff_pct']:.1f}%  "
        f"offset={pw['window_offset_hours']:.1f}h\n"
        f"  Intensity: {od['intensity_rel_diff_pct']:.1f}%  "
        f"Area: {od['area_rel_diff_pct']:.1f}%  "
        f"Volume: {od['volume_rel_diff_pct']:.1f}%\n"
        f"  Centroid displacement: {sr['displacement_km']:.1f} km"
    )

    lagrangian_results = {
        "peak_window": peak_window,
        "object_decomp": object_decomp,
        "storm_relative": storm_relative,
    }
    return lagrangian_results


# ---------------------------------------------------------------------------
# Attribution Visualization Functions
# ---------------------------------------------------------------------------

def plot_peak_window_comparison(
    peak_window_results: Dict,
    config: "AnalysisConfig",
    save_path: Optional[Path] = None,
):
    """
    Three-panel map comparing accumulated precipitation over each scenario's
    own peak window — removing any timing offset bias.

    Panels: factual window accum | counterfactual window accum | difference.
    Panel titles include the exact time window for each scenario and their
    domain-mean accumulations for quick reading.

    Parameters
    ----------
    peak_window_results : dict  — from compute_peak_window_comparison()
    config : AnalysisConfig
    save_path : Path, optional

    Returns
    -------
    fig, axes
    """
    hw = peak_window_results["hist_window"]
    cw = peak_window_results["cont_window"]
    diff = peak_window_results["diff_accum"]
    sc = peak_window_results["scalar_stats"]

    lon = hw["window_accum"].longitude.values
    lat = hw["window_accum"].latitude.values

    if HAS_CARTOPY:
        fig, axes = plt.subplots(
            1, 3, figsize=(16, 5), subplot_kw={"projection": ccrs.PlateCarree()}
        )
    else:
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    plot_kwargs = {"transform": ccrs.PlateCarree()} if HAS_CARTOPY else {}

    vmax = float(max(hw["window_accum"].max(), cw["window_accum"].max()))

    im1 = axes[0].pcolormesh(lon, lat, hw["window_accum"], cmap="Blues",
                              vmin=0, vmax=vmax, **plot_kwargs)
    axes[0].set_title(
        f"{config.scenario_labels['hist']}\n"
        f"{hw['start_time']} → {hw['end_time']}\n"
        f"Domain mean: {sc['factual_window_accum_mm']:.1f} mm",
        fontsize=9,
    )
    plt.colorbar(im1, ax=axes[0], label="Accum. precip [mm]", shrink=0.8)

    im2 = axes[1].pcolormesh(lon, lat, cw["window_accum"], cmap="Blues",
                              vmin=0, vmax=vmax, **plot_kwargs)
    axes[1].set_title(
        f"{config.scenario_labels['cont']}\n"
        f"{cw['start_time']} → {cw['end_time']}\n"
        f"Domain mean: {sc['counterfactual_window_accum_mm']:.1f} mm",
        fontsize=9,
    )
    plt.colorbar(im2, ax=axes[1], label="Accum. precip [mm]", shrink=0.8)

    abs_max = float(np.abs(diff).max())
    abs_max = abs_max if abs_max > 0 else 1.0
    im3 = axes[2].pcolormesh(lon, lat, diff, cmap="RdBu_r",
                              vmin=-abs_max, vmax=abs_max, **plot_kwargs)
    offset_str = (
        f"CF peaks {abs(sc['window_offset_hours']):.1f}h "
        f"{'earlier' if sc['window_offset_hours'] > 0 else 'later'}"
        if sc["window_offset_hours"] != 0 else "same timing"
    )
    axes[2].set_title(
        f"Difference (Factual − Counterfactual)\n"
        f"Δ = {sc['absolute_diff_mm']:.1f} mm  ({sc['relative_diff_pct']:.1f}%)\n"
        f"{offset_str}",
        fontsize=9,
    )
    plt.colorbar(im3, ax=axes[2], label="Δ Accum. precip [mm]", shrink=0.8)

    for ax in axes:
        if HAS_CARTOPY:
            ax.coastlines(resolution="10m", linewidth=0.8)
            ax.add_feature(cfeature.BORDERS, linestyle=":", linewidth=0.6)
        else:
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")

    fig.suptitle(
        f"Peak {sc['window_offset_hours'] and peak_window_results['window_hours']:.0f}h "
        f"Accumulation Window Comparison — {config.region_name}\n"
        f"Each scenario uses its own best {peak_window_results['window_hours']:.0f}h window",
        fontsize=12, fontweight="bold",
    )
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig, axes


def plot_object_decomposition(
    decomp_results: Dict,
    storm_core_results: Dict,
    config: "AnalysisConfig",
    save_path: Optional[Path] = None,
):
    """
    Three-panel time series plot decomposing the storm attribution into
    Intensity (P95 in core), Area (core footprint km²), and Volume (mm·km²).

    Each panel shows the factual (blue) and counterfactual (red dashed) time
    series with fill_between, plus a text annotation of the event-mean
    relative difference.

    Parameters
    ----------
    decomp_results : dict  — from compute_object_decomposition()
    storm_core_results : dict  — from compute_storm_core_analysis()
    config : AnalysisConfig
    save_path : Path, optional

    Returns
    -------
    fig, axes
    """
    time = decomp_results["time"]
    od = decomp_results["scalar_stats"]
    sc = storm_core_results["scalar_stats"]
    vol_hist = storm_core_results["factual_core_volume_ts"]
    vol_cont = storm_core_results["counterfactual_core_volume_ts"]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=False)

    def _ts_panel(ax, ts_hist, ts_cont, ylabel, mean_hist, mean_cont, rel_pct, title):
        ax.plot(time, ts_hist, color="#1f77b4", lw=2,
                label=config.scenario_labels["hist"])
        ax.plot(time, ts_cont, color="#d62728", lw=2, ls="--",
                label=config.scenario_labels["cont"])
        ax.fill_between(time, ts_hist, ts_cont,
                        where=~(np.isnan(ts_hist) | np.isnan(ts_cont)),
                        alpha=0.2, color="purple")
        ax.axhline(mean_hist, color="#1f77b4", ls=":", lw=1.2, alpha=0.7)
        ax.axhline(mean_cont, color="#d62728", ls=":", lw=1.2, alpha=0.7)
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Time")
        ax.tick_params(axis="x", rotation=30)
        ax.grid(True, alpha=0.3)
        sign = "+" if rel_pct >= 0 else ""
        ax.text(
            0.03, 0.97,
            f"F mean: {mean_hist:.1f}\nCF mean: {mean_cont:.1f}\nΔ: {sign}{rel_pct:.1f}%",
            transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8),
            fontfamily="monospace",
        )
        ax.legend(fontsize=7, loc="upper right")

    _ts_panel(
        axes[0],
        decomp_results["intensity_hist_ts"],
        decomp_results["intensity_cont_ts"],
        "P95 precip in core [mm]",
        od["mean_intensity_hist_mm"], od["mean_intensity_cont_mm"],
        od["intensity_rel_diff_pct"],
        "Intensity (P95 in core)",
    )
    _ts_panel(
        axes[1],
        decomp_results["area_hist_ts_km2"],
        decomp_results["area_cont_ts_km2"],
        "Core area [km²]",
        od["mean_area_hist_km2"], od["mean_area_cont_km2"],
        od["area_rel_diff_pct"],
        "Area (core footprint)",
    )
    _ts_panel(
        axes[2],
        vol_hist.values,
        vol_cont.values,
        "Core volume [mm·km²]",
        od["volume_factual_mm_km2"] / max(sc["n_active_timesteps_factual"], 1),
        od["volume_cont_mm_km2"] / max(sc["n_active_timesteps_counterfactual"], 1),
        od["volume_rel_diff_pct"],
        "Volume (Intensity × Area)",
    )

    fig.suptitle(
        f"Storm Attribution Decomposition — {config.region_name}\n"
        f"Intensity: {'+' if od['intensity_rel_diff_pct'] >= 0 else ''}{od['intensity_rel_diff_pct']:.1f}%  |  "
        f"Area: {'+' if od['area_rel_diff_pct'] >= 0 else ''}{od['area_rel_diff_pct']:.1f}%  |  "
        f"Volume: {'+' if od['volume_rel_diff_pct'] >= 0 else ''}{od['volume_rel_diff_pct']:.1f}%",
        fontsize=12, fontweight="bold",
    )
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig, axes


def plot_storm_relative_composite(
    composite_results: Dict,
    config: "AnalysisConfig",
    save_path: Optional[Path] = None,
):
    """
    Four-panel figure for storm-relative compositing.

    [0,0] Factual composite — precipitation field at peak window
    [0,1] Counterfactual shifted to factual centroid (same scale)
    [1,0] Storm-relative difference = factual - shifted CF
          Expected: concentric pattern, not a dipole
    [1,1] Radial profiles — azimuthal mean precipitation vs distance (km)
          from the shared factual centroid; shows rain-shield expansion/contraction

    Parameters
    ----------
    composite_results : dict  — from compute_storm_relative_composite()
    config : AnalysisConfig
    save_path : Path, optional

    Returns
    -------
    fig, axes (2 x 2)
    """
    cr = composite_results
    lat, lon = cr["lat"], cr["lon"]
    sc = cr["scalar_stats"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    vmax = float(max(
        np.nanmax(cr["factual_composite"]),
        np.nanmax(cr["cf_composite_shifted"]),
    ))
    vmax = max(vmax, 1.0)

    # [0,0] Factual composite
    ax = axes[0, 0]
    im1 = ax.pcolormesh(lon, lat, cr["factual_composite"],
                        cmap="Blues", vmin=0, vmax=vmax, shading="auto")
    ax.scatter(cr["factual_centroid_lon"], cr["factual_centroid_lat"],
               s=250, marker="*", color="yellow", edgecolors="black",
               linewidths=1.0, zorder=10, label="Factual centroid")
    ax.set_title(
        f"{config.scenario_labels['hist']}\nComposite ({cr['peak_window_str']})",
        fontsize=9, fontweight="bold",
    )
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    ax.legend(fontsize=7, loc="upper right")
    plt.colorbar(im1, ax=ax, label="Precip [mm]", shrink=0.85)

    # [0,1] Counterfactual shifted
    ax = axes[0, 1]
    im2 = ax.pcolormesh(lon, lat, cr["cf_composite_shifted"],
                        cmap="Blues", vmin=0, vmax=vmax, shading="auto")
    ax.scatter(cr["factual_centroid_lon"], cr["factual_centroid_lat"],
               s=250, marker="*", color="yellow", edgecolors="black",
               linewidths=1.0, zorder=10, label="Shared centroid (factual)")
    ax.scatter(cr["cf_centroid_lon"], cr["cf_centroid_lat"],
               s=120, marker="o", color="red", edgecolors="white",
               linewidths=0.8, zorder=9, label=f"CF original centroid")
    disp_dir = []
    if abs(cr["displacement_lat_deg"]) > 0.01:
        disp_dir.append("N" if cr["displacement_lat_deg"] > 0 else "S")
    if abs(cr["displacement_lon_deg"]) > 0.01:
        disp_dir.append("E" if cr["displacement_lon_deg"] > 0 else "W")
    disp_str = "".join(disp_dir) if disp_dir else "~0"
    ax.set_title(
        f"{config.scenario_labels['cont']} — shifted to factual centroid\n"
        f"Displacement: {sc['displacement_km']:.1f} km {disp_str}",
        fontsize=9, fontweight="bold",
    )
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    ax.legend(fontsize=7, loc="upper right")
    plt.colorbar(im2, ax=ax, label="Precip [mm]", shrink=0.85)

    # [1,0] Storm-relative difference
    ax = axes[1, 0]
    abs_max = max(abs(float(np.nanmax(cr["storm_relative_diff"]))),
                  abs(float(np.nanmin(cr["storm_relative_diff"]))), 1.0)
    im3 = ax.pcolormesh(lon, lat, cr["storm_relative_diff"],
                        cmap="RdBu_r", vmin=-abs_max, vmax=abs_max, shading="auto")
    ax.scatter(cr["factual_centroid_lon"], cr["factual_centroid_lat"],
               s=250, marker="*", color="yellow", edgecolors="black",
               linewidths=1.0, zorder=10)
    ax.set_title(
        f"Storm-Relative Difference (Factual − Shifted CF)\n"
        f"Mean Δ = {sc['storm_relative_diff_mean_mm']:.2f} mm  "
        f"Max = {sc['storm_relative_diff_max_mm']:.1f}  "
        f"Min = {sc['storm_relative_diff_min_mm']:.1f}",
        fontsize=9, fontweight="bold",
    )
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    plt.colorbar(im3, ax=ax, label="Δ Precip [mm]", shrink=0.85)

    # [1,1] Radial profiles
    ax = axes[1, 1]
    rp_h = cr["radial_profile_hist"]
    rp_c = cr["radial_profile_cont"]
    ax.plot(rp_h["r_km"], rp_h["precip_mean"], color="#1f77b4", lw=2.5,
            label=config.scenario_labels["hist"])
    ax.plot(rp_c["r_km"], rp_c["precip_mean"], color="#d62728", lw=2.5, ls="--",
            label=f"{config.scenario_labels['cont']} (shifted)")
    ax.fill_between(
        rp_h["r_km"],
        np.where(np.isnan(rp_h["precip_mean"]), 0, rp_h["precip_mean"]),
        np.where(np.isnan(rp_c["precip_mean"]), 0, rp_c["precip_mean"]),
        alpha=0.2, color="purple",
    )
    ax.axvline(50, color="gray", ls=":", lw=1, alpha=0.7,
               label="~50 km (convective-core scale)")
    ax.set_xlabel("Distance from centroid [km]")
    ax.set_ylabel("Azimuthal mean precip [mm]")
    ax.set_title(
        "Radial Precipitation Profile\n(both referenced to factual centroid)",
        fontsize=9, fontweight="bold",
    )
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)

    fig.suptitle(
        f"Storm-Relative Composite — {config.region_name}\n"
        f"Peak: {cr['peak_window_str']}  |  "
        f"Centroid displacement removed: {sc['displacement_km']:.1f} km {disp_str}",
        fontsize=12, fontweight="bold",
    )
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig, axes


def plot_lagrangian_dashboard(
    lagrangian_results: Dict,
    config: "AnalysisConfig",
    save_path: Optional[Path] = None,
):
    """
    Single consolidated dashboard combining all three attribution methods.

    Layout (4 rows via GridSpec):

    Row 1 — Temporal Relaxation (Part 1)
        [Factual 24h peak accum] [CF 24h peak accum] [Difference map]
        Panel titles carry exact time windows; difference panel annotates
        the timing offset between scenarios.

    Row 2 — Object Decomposition (Part 2)
        [Intensity time series] [Area time series] [Attribution text summary]
        Each time series shows factual (blue) vs counterfactual (red dashed)
        with fill_between.  The text box reports Intensity %, Area %, Volume %.

    Row 3 — Storm-Relative Compositing (Part 3, maps)
        [Factual ±3h composite] [CF shifted to factual centroid] [Diff map RdBu]
        Centroids marked with stars; CF pre-shift centroid shown as red circle.

    Row 4 — Radial Profile (Part 3, profile)
        Wide single panel: azimuthal mean precipitation (mm) vs distance from
        factual centroid (km, 0–150 km) for factual (blue) and shifted CF (red).
        fill_between shows structural expansion / contraction of the rain shield.

    Parameters
    ----------
    lagrangian_results : dict  — from compute_lagrangian_attribution()
    config : AnalysisConfig
    save_path : Path, optional

    Returns
    -------
    fig
    """
    pw_res = lagrangian_results["peak_window"]
    od_res = lagrangian_results["object_decomp"]
    sr_res = lagrangian_results["storm_relative"]

    pw = pw_res["scalar_stats"]
    od = od_res["scalar_stats"]
    sr = sr_res["scalar_stats"]


    clr_hist = "#1f77b4"
    clr_cont = "#d62728"
    lbl_hist = config.scenario_labels["hist"]
    lbl_cont = config.scenario_labels["cont"]

    fig = plt.figure(figsize=(18, 22))
    gs = GridSpec(
        4, 3,
        figure=fig,
        height_ratios=[1, 1, 1, 0.75],
        hspace=0.42,
        wspace=0.32,
    )

    # ------------------------------------------------------------------
    # ROW 1 — Peak window accumulation maps
    # ------------------------------------------------------------------
    hw = pw_res["hist_window"]
    cw = pw_res["cont_window"]
    diff_accum = pw_res["diff_accum"]

    lon_pw = hw["window_accum"].longitude.values
    lat_pw = hw["window_accum"].latitude.values
    vmax_pw = float(max(hw["window_accum"].max(), cw["window_accum"].max()))

    ax_r1 = [fig.add_subplot(gs[0, c]) for c in range(3)]

    im = ax_r1[0].pcolormesh(lon_pw, lat_pw, hw["window_accum"],
                              cmap="Blues", vmin=0, vmax=vmax_pw, shading="auto")
    plt.colorbar(im, ax=ax_r1[0], label="mm", shrink=0.85)
    ax_r1[0].set_title(
        f"{lbl_hist}\n{hw['start_time']} → {hw['end_time']}\n"
        f"Domain mean: {pw['factual_window_accum_mm']:.1f} mm",
        fontsize=8, fontweight="bold",
    )

    im2 = ax_r1[1].pcolormesh(lon_pw, lat_pw, cw["window_accum"],
                               cmap="Blues", vmin=0, vmax=vmax_pw, shading="auto")
    plt.colorbar(im2, ax=ax_r1[1], label="mm", shrink=0.85)
    ax_r1[1].set_title(
        f"{lbl_cont}\n{cw['start_time']} → {cw['end_time']}\n"
        f"Domain mean: {pw['counterfactual_window_accum_mm']:.1f} mm",
        fontsize=8, fontweight="bold",
    )

    abs_max_pw = max(float(np.abs(diff_accum).max()), 1.0)
    im3 = ax_r1[2].pcolormesh(lon_pw, lat_pw, diff_accum,
                               cmap="RdBu_r", vmin=-abs_max_pw, vmax=abs_max_pw, shading="auto")
    plt.colorbar(im3, ax=ax_r1[2], label="Δ mm", shrink=0.85)
    sign_str = (f"CF peaks {abs(pw['window_offset_hours']):.1f}h "
                f"{'earlier' if pw['window_offset_hours'] > 0 else 'later'}"
                if pw["window_offset_hours"] != 0 else "same timing")
    ax_r1[2].set_title(
        f"Difference (F − CF)\nΔ = {pw['absolute_diff_mm']:.1f} mm "
        f"({pw['relative_diff_pct']:+.1f}%)\n{sign_str}",
        fontsize=8, fontweight="bold",
    )
    for ax in ax_r1:
        ax.set_xlabel("Lon"); ax.set_ylabel("Lat")

    # ------------------------------------------------------------------
    # ROW 2 — Object decomposition time series + summary text
    # ------------------------------------------------------------------
    ax_int = fig.add_subplot(gs[1, 0])
    ax_area = fig.add_subplot(gs[1, 1])
    ax_text = fig.add_subplot(gs[1, 2])

    time = od_res["time"]

    def _ts(ax, ts_h, ts_c, ylabel, mean_h, mean_c, rel_pct, title):
        mask = ~(np.isnan(ts_h) | np.isnan(ts_c))
        ax.plot(time, ts_h, color=clr_hist, lw=2, label=lbl_hist)
        ax.plot(time, ts_c, color=clr_cont, lw=2, ls="--", label=lbl_cont)
        if mask.any():
            ax.fill_between(time, ts_h, ts_c, where=mask, alpha=0.2, color="purple")
        ax.axhline(mean_h, color=clr_hist, ls=":", lw=1.2, alpha=0.7)
        ax.axhline(mean_c, color=clr_cont, ls=":", lw=1.2, alpha=0.7)
        ax.set_title(title, fontsize=9, fontweight="bold")
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(axis="x", rotation=30, labelsize=7)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, loc="upper right")
        sign = "+" if rel_pct >= 0 else ""
        ax.text(0.03, 0.97,
                f"F: {mean_h:.1f}  CF: {mean_c:.1f}\nΔ: {sign}{rel_pct:.1f}%",
                transform=ax.transAxes, fontsize=8, va="top", fontfamily="monospace",
                bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    _ts(ax_int,
        od_res["intensity_hist_ts"], od_res["intensity_cont_ts"],
        "P95 in core [mm]",
        od["mean_intensity_hist_mm"], od["mean_intensity_cont_mm"],
        od["intensity_rel_diff_pct"],
        "Intensity (P95 within core)")

    _ts(ax_area,
        od_res["area_hist_ts_km2"], od_res["area_cont_ts_km2"],
        "Core area [km²]",
        od["mean_area_hist_km2"], od["mean_area_cont_km2"],
        od["area_rel_diff_pct"],
        "Area (core footprint)")

    # Text summary panel
    ax_text.axis("off")
    cc_expected = 7.0  # % per degree warming (Clausius-Clapeyron)
    summary = (
        "ATTRIBUTION BREAKDOWN\n"
        + "─" * 28 + "\n"
        f"  Intensity (thermo.):\n"
        f"    {od['intensity_rel_diff_pct']:+.1f}%\n"
        f"    (C-C expectation: +{cc_expected:.0f}%/K)\n\n"
        f"  Area (dynamic):\n"
        f"    {od['area_rel_diff_pct']:+.1f}%\n\n"
        f"  Volume (total):\n"
        f"    {od['volume_rel_diff_pct']:+.1f}%\n"
        + "─" * 28 + "\n"
        f"  Peak window Δ:\n"
        f"    {pw['relative_diff_pct']:+.1f}% ({pw['absolute_diff_mm']:+.1f} mm)\n"
        f"  Timing offset:\n"
        f"    {pw['window_offset_hours']:+.1f} h\n"
        + "─" * 28 + "\n"
        f"  Centroid shift:\n"
        f"    {sr['displacement_km']:.1f} km\n"
        f"    (Δlat={sr['displacement_lat_deg']:+.2f}°,\n"
        f"     Δlon={sr['displacement_lon_deg']:+.2f}°)"
    )
    ax_text.text(
        0.08, 0.95, summary,
        transform=ax_text.transAxes, fontsize=9, va="top",
        fontfamily="monospace",
        bbox=dict(boxstyle="round", facecolor="#f0f4ff", edgecolor="#3a5a9e",
                  linewidth=1.5, alpha=0.95),
    )

    # ------------------------------------------------------------------
    # ROW 3 — Storm-relative compositing maps
    # ------------------------------------------------------------------
    lat_sr = sr_res["lat"]
    lon_sr = sr_res["lon"]
    clat_f = sr_res["factual_centroid_lat"]
    clon_f = sr_res["factual_centroid_lon"]
    clat_c = sr_res["cf_centroid_lat"]
    clon_c = sr_res["cf_centroid_lon"]

    vmax_sr = max(
        float(np.nanmax(sr_res["factual_composite"])),
        float(np.nanmax(sr_res["cf_composite_shifted"])),
        1.0,
    )

    ax_r3 = [fig.add_subplot(gs[2, c]) for c in range(3)]

    im_f = ax_r3[0].pcolormesh(lon_sr, lat_sr, sr_res["factual_composite"],
                                cmap="Blues", vmin=0, vmax=vmax_sr, shading="auto")
    ax_r3[0].scatter(clon_f, clat_f, s=220, marker="*", color="yellow",
                     edgecolors="black", linewidths=0.9, zorder=10)
    ax_r3[0].set_title(f"{lbl_hist}\nPeak composite ({sr_res['peak_window_str']})",
                       fontsize=8, fontweight="bold")
    plt.colorbar(im_f, ax=ax_r3[0], label="mm", shrink=0.85)

    im_c = ax_r3[1].pcolormesh(lon_sr, lat_sr, sr_res["cf_composite_shifted"],
                                cmap="Blues", vmin=0, vmax=vmax_sr, shading="auto")
    ax_r3[1].scatter(clon_f, clat_f, s=220, marker="*", color="yellow",
                     edgecolors="black", linewidths=0.9, zorder=10,
                     label="Shared centroid")
    ax_r3[1].scatter(clon_c, clat_c, s=100, marker="o", color="red",
                     edgecolors="white", linewidths=0.8, zorder=9,
                     label="CF original centroid")
    # Arrow showing the shift applied
    ax_r3[1].annotate(
        "", xy=(clon_f, clat_f), xytext=(clon_c, clat_c),
        arrowprops=dict(arrowstyle="->", color="red", lw=1.5),
        zorder=11,
    )
    ax_r3[1].set_title(
        f"{lbl_cont} — shifted to factual centroid\n"
        f"Displacement: {sr['displacement_km']:.1f} km removed",
        fontsize=8, fontweight="bold",
    )
    ax_r3[1].legend(fontsize=6, loc="upper right")
    plt.colorbar(im_c, ax=ax_r3[1], label="mm", shrink=0.85)

    diff_sr = sr_res["storm_relative_diff"]
    abs_max_sr = max(abs(float(np.nanmax(diff_sr))), abs(float(np.nanmin(diff_sr))), 1.0)
    im_d = ax_r3[2].pcolormesh(lon_sr, lat_sr, diff_sr,
                                cmap="RdBu_r", vmin=-abs_max_sr, vmax=abs_max_sr, shading="auto")
    ax_r3[2].scatter(clon_f, clat_f, s=220, marker="*", color="yellow",
                     edgecolors="black", linewidths=0.9, zorder=10)
    ax_r3[2].set_title(
        f"Storm-Relative Difference (F − shifted CF)\n"
        f"Mean Δ = {sr['storm_relative_diff_mean_mm']:.2f} mm  "
        f"Max = {sr['storm_relative_diff_max_mm']:.1f}",
        fontsize=8, fontweight="bold",
    )
    plt.colorbar(im_d, ax=ax_r3[2], label="Δ mm", shrink=0.85)

    for ax in ax_r3:
        ax.set_xlabel("Lon"); ax.set_ylabel("Lat")

    # ------------------------------------------------------------------
    # ROW 4 — Radial profiles (wide, spanning all 3 columns)
    # ------------------------------------------------------------------
    ax_rad = fig.add_subplot(gs[3, :])

    rp_h = sr_res["radial_profile_hist"]
    rp_c = sr_res["radial_profile_cont"]
    r_h, mp_h = rp_h["r_km"], rp_h["precip_mean"]
    r_c, mp_c = rp_c["r_km"], rp_c["precip_mean"]

    # Clip to 150 km as requested
    mask_h = r_h <= 150
    mask_c = r_c <= 150

    ax_rad.plot(r_h[mask_h], mp_h[mask_h], color=clr_hist, lw=2.5, label=lbl_hist)
    ax_rad.plot(r_c[mask_c], mp_c[mask_c], color=clr_cont, lw=2.5, ls="--",
                label=f"{lbl_cont} (shifted)")

    # fill_between on the common r grid
    r_common = r_h[mask_h]
    mp_h_cl = mp_h[mask_h]
    mp_c_cl = np.interp(r_common, r_c[mask_c], mp_c[mask_c],
                        left=np.nan, right=np.nan)
    valid = ~(np.isnan(mp_h_cl) | np.isnan(mp_c_cl))
    if valid.any():
        ax_rad.fill_between(
            r_common[valid], mp_h_cl[valid], mp_c_cl[valid],
            where=mp_h_cl[valid] >= mp_c_cl[valid],
            alpha=0.25, color=clr_hist, label="Factual > CF (intensification)",
        )
        ax_rad.fill_between(
            r_common[valid], mp_h_cl[valid], mp_c_cl[valid],
            where=mp_h_cl[valid] < mp_c_cl[valid],
            alpha=0.25, color=clr_cont, label="CF > Factual",
        )

    ax_rad.axvline(50, color="gray", ls=":", lw=1.2, alpha=0.7,
                   label="50 km (~convective core scale)")
    ax_rad.set_xlabel("Distance from factual centroid [km]", fontsize=10)
    ax_rad.set_ylabel("Azimuthal mean precipitation [mm]", fontsize=10)
    ax_rad.set_title(
        "Radial Precipitation Profile  —  Both referenced to factual centroid\n"
        "Blue fill = factual enhancement; red fill = CF stronger at that radius",
        fontsize=9, fontweight="bold",
    )
    ax_rad.set_xlim(0, 150)
    ax_rad.set_ylim(bottom=0)
    ax_rad.legend(fontsize=8, loc="upper right")
    ax_rad.grid(True, alpha=0.3)

    # ------------------------------------------------------------------
    # Overall title
    # ------------------------------------------------------------------
    fig.suptitle(
        f"Lagrangian Storm Attribution Dashboard — {config.region_name}\n"
        f"Factual vs Counterfactual  |  "
        f"Intensity: {od['intensity_rel_diff_pct']:+.1f}%  "
        f"Area: {od['area_rel_diff_pct']:+.1f}%  "
        f"Volume: {od['volume_rel_diff_pct']:+.1f}%  "
        f"|  Peak-window Δ: {pw['relative_diff_pct']:+.1f}%",
        fontsize=13, fontweight="bold", y=1.005,
    )

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig


def run_full_analysis(config: AnalysisConfig, variable: str) -> Dict:
    """Run the complete analysis for a given variable."""
    print(f"\n{'='*60}")
    print(f"Analyzing: {config.variable_info[variable]['name']}")
    print(f"{'='*60}")

    # Load data
    data = load_all_scenarios(config, variable)
    var_info = config.variable_info[variable]

    results = {
        "variable": variable,
        "variable_info": var_info,
        "data": data,
        "scenarios": {},
    }

    # Compute statistics for each scenario
    for scenario, da in data.items():
        print(f"\nComputing statistics for {config.scenario_labels[scenario]}...")
        scenario_stats = {
            "temporal": compute_temporal_statistics(da, var_info),
            "percentiles": compute_percentile_statistics(da),
        }
        results["scenarios"][scenario] = scenario_stats

    # Compute comparison statistics
    print("\nComputing comparison statistics...")
    results["comparison"] = compute_comparison_statistics(
        data["hist"], data["cont"], var_info
    )

    # Dynamic storm core analysis (precipitation only)
    if variable == "tp":
        print("\nComputing dynamic storm core analysis...")
        try:
            results["storm_core"] = compute_storm_core_analysis(
                data["hist"], data["cont"], config
            )
            sc = results["storm_core"]["scalar_stats"]
            print(
                f"  Factual core volume:        {sc['factual_total_core_volume_mm_km2']:.0f} mm·km²\n"
                f"  Counterfactual core volume: {sc['counterfactual_total_core_volume_mm_km2']:.0f} mm·km²\n"
                f"  Relative difference:        {sc['relative_diff_pct']:.1f}%"
            )
        except Exception as exc:
            print(f"  WARNING: Storm core analysis failed: {exc}")
            results["storm_core"] = None

    # Three-method storm attribution (precipitation only, requires storm core)
    if variable == "tp" and results.get("storm_core") is not None:
        print("\nComputing Lagrangian storm attribution (Methods 1–3)...")
        try:
            results["storm_attribution"] = compute_lagrangian_attribution(
                data["hist"], data["cont"], results["storm_core"], config
            )
            # Print a compact summary immediately so it's visible in the terminal
            attr = results["storm_attribution"]
            pw = attr["peak_window"]["scalar_stats"]
            od = attr["object_decomp"]["scalar_stats"]
            sr = attr["storm_relative"]["scalar_stats"]
            print("\n" + "=" * 60)
            print("  LAGRANGIAN ATTRIBUTION SUMMARY")
            print("=" * 60)
            print(f"  Method 1 — Peak accumulation window ({attr['peak_window']['window_hours']:.0f} h)")
            print(f"    Factual:        {pw['factual_window_accum_mm']:.1f} mm  "
                  f"({pw['factual_window_start']} → {pw['factual_window_end']})")
            print(f"    Counterfactual: {pw['counterfactual_window_accum_mm']:.1f} mm  "
                  f"({pw['counterfactual_window_start']} → {pw['counterfactual_window_end']})")
            print(f"    Change: {pw['absolute_diff_mm']:+.1f} mm  ({pw['relative_diff_pct']:+.1f}%)"
                  f"   |  Timing offset: {pw['window_offset_hours']:.1f} h")
            print(f"  Method 2 — Object decomposition")
            print(f"    Intensity (P95 in core): {od['intensity_rel_diff_pct']:+.1f}%"
                  f"  ({od['mean_intensity_hist_mm']:.1f} → {od['mean_intensity_cont_mm']:.1f} mm)")
            print(f"    Area (core footprint):   {od['area_rel_diff_pct']:+.1f}%"
                  f"  ({od['mean_area_hist_km2']:.0f} → {od['mean_area_cont_km2']:.0f} km²)")
            print(f"    Volume (mm·km²):         {od['volume_rel_diff_pct']:+.1f}%")
            print(f"  Method 3 — Storm-relative composite")
            print(f"    Centroid displacement removed: {sr['displacement_km']:.1f} km"
                  f"  (Δlat={sr['displacement_lat_deg']:+.3f}°, Δlon={sr['displacement_lon_deg']:+.3f}°)")
            print(f"    Storm-relative Δ precip:  mean={sr['storm_relative_diff_mean_mm']:+.2f} mm"
                  f"  max={sr['storm_relative_diff_max_mm']:+.2f} mm"
                  f"  min={sr['storm_relative_diff_min_mm']:+.2f} mm")
            print("=" * 60)
        except Exception as exc:
            print(f"  WARNING: Storm attribution analysis failed: {exc}")
            import traceback
            traceback.print_exc()
            results["storm_attribution"] = None

    return results


# =============================================================================
# Storm Core Detection Functions
# =============================================================================


def detect_storm_core_single_timestep(
    field_2d: np.ndarray,
    lat_1d: np.ndarray,
    lon_1d: np.ndarray,
    gaussian_sigma: float = 1.5,
    percentile_threshold: float = 75.0,
    abs_threshold_mm: float = 1.0,
    quiet_threshold_mm: float = 0.5,
) -> Tuple[np.ndarray, Optional[Tuple[float, float]], Dict]:
    """
    Detect the storm precipitation core in a single 2D field.

    Uses Gaussian smoothing to suppress sub-grid noise, a percentile threshold
    to identify heavy precipitation pixels, and connected-component labeling
    to isolate the dominant storm blob. The largest connected component
    (by pixel count) is selected as the storm core.

    Parameters
    ----------
    field_2d : np.ndarray
        2D precipitation field (lat x lon) in mm. May contain NaN.
    lat_1d : np.ndarray
        1D latitude coordinate array.
    lon_1d : np.ndarray
        1D longitude coordinate array.
    gaussian_sigma : float
        Smoothing standard deviation in grid cells (default 1.5).
    percentile_threshold : float
        Percentile of non-zero smoothed pixels used as the lower threshold
        (default 75.0 = P75).
    abs_threshold_mm : float
        Hard minimum precipitation value in mm (default 1.0).
    quiet_threshold_mm : float
        If the P{percentile_threshold} of the smoothed field is below this,
        the timestep is flagged as quiet and an empty mask is returned.

    Returns
    -------
    core_mask_2d : np.ndarray (bool, shape lat x lon)
    centroid : tuple (lat, lon) or None
        Precipitation-weighted centroid in geographic coordinates.
    diagnostics : dict
        Keys: is_quiet, n_components, core_pixel_count, core_mean_precip,
              threshold_used.
    """
    nan_mask = np.isnan(field_2d)
    field_clean = np.where(nan_mask, 0.0, field_2d)

    # Gaussian smoothing to suppress isolated convective noise
    smoothed = gaussian_filter(field_clean.astype(float), sigma=gaussian_sigma)
    smoothed[nan_mask] = 0.0  # keep NaN regions clean

    # Determine threshold from active (non-zero) smoothed pixels
    active_vals = smoothed[smoothed > 0]
    if len(active_vals) == 0:
        return (
            np.zeros_like(field_2d, dtype=bool),
            None,
            {"is_quiet": True, "n_components": 0, "core_pixel_count": 0,
             "core_mean_precip": 0.0, "threshold_used": 0.0},
        )

    pct_val = float(np.percentile(active_vals, percentile_threshold))
    if pct_val < quiet_threshold_mm:
        return (
            np.zeros_like(field_2d, dtype=bool),
            None,
            {"is_quiet": True, "n_components": 0, "core_pixel_count": 0,
             "core_mean_precip": 0.0, "threshold_used": pct_val},
        )

    threshold = max(abs_threshold_mm, pct_val)

    # Binary mask and connected-component labeling
    binary = (smoothed >= threshold) & (~nan_mask)
    labeled_array, n_components = ndimage_label(binary)

    if n_components == 0:
        return (
            np.zeros_like(field_2d, dtype=bool),
            None,
            {"is_quiet": False, "n_components": 0, "core_pixel_count": 0,
             "core_mean_precip": 0.0, "threshold_used": threshold},
        )

    # Select the largest connected component
    component_sizes = np.bincount(labeled_array.ravel())[1:]  # exclude background (label 0)
    largest_label = int(np.argmax(component_sizes)) + 1
    core_mask_2d = labeled_array == largest_label

    # Precipitation-weighted centroid (use original field, not smoothed)
    row_idx, col_idx = np.where(core_mask_2d)
    weights = field_clean[row_idx, col_idx]
    weight_sum = weights.sum()
    if weight_sum > 0:
        centroid_lat = float(np.average(lat_1d[row_idx], weights=weights))
        centroid_lon = float(np.average(lon_1d[col_idx], weights=weights))
        centroid = (centroid_lat, centroid_lon)
    else:
        centroid = (float(lat_1d[row_idx].mean()), float(lon_1d[col_idx].mean()))

    core_precip_vals = field_clean[core_mask_2d]
    diagnostics = {
        "is_quiet": False,
        "n_components": n_components,
        "core_pixel_count": int(core_mask_2d.sum()),
        "core_mean_precip": float(core_precip_vals.mean()) if len(core_precip_vals) else 0.0,
        "threshold_used": threshold,
    }
    return core_mask_2d, centroid, diagnostics


def compute_core_mask_timeseries(
    da: xr.DataArray,
    gaussian_sigma: float = 1.5,
    percentile_threshold: float = 75.0,
    abs_threshold_mm: float = 1.0,
    quiet_threshold_mm: float = 0.5,
) -> Tuple[xr.DataArray, pd.DataFrame]:
    """
    Compute the dynamic storm core mask for every timestep in a DataArray.

    At each timestep the precipitation field is analysed independently via
    detect_storm_core_single_timestep(). The resulting boolean masks are
    stacked back into a DataArray with the same coordinates as the input.

    Parameters
    ----------
    da : xr.DataArray
        Precipitation DataArray with dims (time, latitude, longitude) in mm.
    gaussian_sigma, percentile_threshold, abs_threshold_mm, quiet_threshold_mm
        Passed directly to detect_storm_core_single_timestep().

    Returns
    -------
    core_mask_da : xr.DataArray (bool, time x latitude x longitude)
    track_df : pd.DataFrame
        Columns: time, centroid_lat, centroid_lon, core_pixel_count,
                 core_mean_precip_mm, is_quiet, threshold_used, n_components.
    """
    lat_1d = da.latitude.values
    lon_1d = da.longitude.values
    n_time, n_lat, n_lon = da.shape

    core_mask_arr = np.zeros((n_time, n_lat, n_lon), dtype=bool)
    records = []

    for i in range(n_time):
        field_2d = da.isel(time=i).values
        mask, centroid, diags = detect_storm_core_single_timestep(
            field_2d, lat_1d, lon_1d,
            gaussian_sigma=gaussian_sigma,
            percentile_threshold=percentile_threshold,
            abs_threshold_mm=abs_threshold_mm,
            quiet_threshold_mm=quiet_threshold_mm,
        )
        core_mask_arr[i] = mask
        records.append({
            "time": da.time.values[i],
            "centroid_lat": centroid[0] if centroid is not None else np.nan,
            "centroid_lon": centroid[1] if centroid is not None else np.nan,
            "core_pixel_count": diags["core_pixel_count"],
            "core_mean_precip_mm": diags["core_mean_precip"],
            "is_quiet": diags["is_quiet"],
            "threshold_used": diags["threshold_used"],
            "n_components": diags["n_components"],
        })

    core_mask_da = xr.DataArray(
        core_mask_arr,
        coords=da.coords,
        dims=da.dims,
        attrs={"description": "Dynamic storm core boolean mask"},
    )
    track_df = pd.DataFrame(records)
    return core_mask_da, track_df


def compute_core_volume_timeseries(
    da: xr.DataArray,
    core_mask_da: xr.DataArray,
) -> xr.DataArray:
    """
    Compute the spatial integral of precipitation within the storm core
    at each timestep, with cosine-latitude area weighting.

    Parameters
    ----------
    da : xr.DataArray
        Precipitation DataArray (time, latitude, longitude) in mm.
    core_mask_da : xr.DataArray (bool)
        Dynamic core mask from compute_core_mask_timeseries(), same shape.

    Returns
    -------
    core_volume_ts : xr.DataArray (time,)
        Core precipitation volume in mm·km². Zero for quiet timesteps.
        1 mm·km² = 1000 m³ of water.
    """
    dlat = float(np.abs(np.diff(da.latitude.values).mean()))
    dlon = float(np.abs(np.diff(da.longitude.values).mean()))
    km_per_deg = 111.32

    # Per-latitude cell area (varies with cos(lat))
    lat_rad = np.deg2rad(da.latitude.values)
    cell_area_km2_1d = dlat * dlon * (km_per_deg ** 2) * np.cos(lat_rad)

    # Broadcast to 2D (lat, lon) then align with DataArray
    cell_area_2d = np.outer(cell_area_km2_1d, np.ones(len(da.longitude)))
    cell_area_da = xr.DataArray(
        cell_area_2d,
        coords={"latitude": da.latitude, "longitude": da.longitude},
        dims=["latitude", "longitude"],
    )

    masked_precip = da.where(core_mask_da, other=0.0)
    core_volume_ts = (masked_precip * cell_area_da).sum(dim=["latitude", "longitude"])
    core_volume_ts.attrs["units"] = "mm*km2"
    core_volume_ts.attrs["description"] = "Storm core precipitation volume (1 mm*km2 = 1000 m3)"
    return core_volume_ts


def compute_storm_core_analysis(
    data_factual: xr.DataArray,
    data_counterfactual: xr.DataArray,
    config: "AnalysisConfig",
) -> Dict:
    """
    Compute the difference in storm core precipitation volume between
    factual and counterfactual scenarios.

    The storm core is detected independently for each scenario at every
    timestep (Option B: independent masks). This is the most physically
    honest approach — each storm is measured on its own terms, avoiding
    the inflation that would result from imposing the (larger) factual
    core mask onto the (potentially smaller) counterfactual field.

    Parameters
    ----------
    data_factual : xr.DataArray
        Historical precipitation DataArray (time, latitude, longitude) in mm.
    data_counterfactual : xr.DataArray
        Control/pre-industrial precipitation DataArray, same shape.
    config : AnalysisConfig
        Reads core_gaussian_sigma, core_percentile_threshold,
        core_abs_threshold_mm, core_quiet_threshold_mm.

    Returns
    -------
    dict with keys:
        factual_core_mask, counterfactual_core_mask   : xr.DataArray (bool)
        factual_track_df, counterfactual_track_df     : pd.DataFrame
        factual_core_volume_ts, counterfactual_core_volume_ts : xr.DataArray (mm*km2)
        core_volume_diff_ts                           : xr.DataArray (mm*km2)
        core_fraction_hist, core_fraction_cont        : xr.DataArray (0-1, lat x lon)
        scalar_stats                                  : dict
    """
    sigma = getattr(config, "core_gaussian_sigma", 1.5)
    pct = getattr(config, "core_percentile_threshold", 75.0)
    abs_thr = getattr(config, "core_abs_threshold_mm", 1.0)
    quiet_thr = getattr(config, "core_quiet_threshold_mm", 0.5)

    print("  Computing factual storm core masks...")
    hist_mask, hist_track = compute_core_mask_timeseries(
        data_factual, gaussian_sigma=sigma, percentile_threshold=pct,
        abs_threshold_mm=abs_thr, quiet_threshold_mm=quiet_thr,
    )
    print("  Computing counterfactual storm core masks...")
    cont_mask, cont_track = compute_core_mask_timeseries(
        data_counterfactual, gaussian_sigma=sigma, percentile_threshold=pct,
        abs_threshold_mm=abs_thr, quiet_threshold_mm=quiet_thr,
    )

    hist_vol_ts = compute_core_volume_timeseries(data_factual, hist_mask)
    cont_vol_ts = compute_core_volume_timeseries(data_counterfactual, cont_mask)

    # Align on time (inner join) in case of any 1-step offset
    hist_vol_ts, cont_vol_ts = xr.align(hist_vol_ts, cont_vol_ts, join="inner")
    core_volume_diff_ts = hist_vol_ts - cont_vol_ts

    # Time-integrated core frequency (fraction of timesteps in core)
    core_fraction_hist = hist_mask.mean(dim="time")
    core_fraction_cont = cont_mask.mean(dim="time")

    # Grid cell area for area statistics
    dlat = float(np.abs(np.diff(data_factual.latitude.values).mean()))
    dlon = float(np.abs(np.diff(data_factual.longitude.values).mean()))
    km_per_deg = 111.32
    mean_lat_rad = float(np.deg2rad(data_factual.latitude.values.mean()))
    cell_area_km2 = dlat * dlon * (km_per_deg ** 2) * np.cos(mean_lat_rad)

    # Scalar statistics
    n_active_hist = int((~hist_track["is_quiet"]).sum())
    n_active_cont = int((~cont_track["is_quiet"]).sum())

    hist_total_vol = float(hist_vol_ts.sum())
    cont_total_vol = float(cont_vol_ts.sum())
    abs_diff = hist_total_vol - cont_total_vol
    rel_diff_pct = (abs_diff / cont_total_vol * 100) if cont_total_vol > 0 else np.nan

    hist_mean_area = (
        float(hist_track.loc[~hist_track["is_quiet"], "core_pixel_count"].mean()) * cell_area_km2
        if n_active_hist > 0 else 0.0
    )
    cont_mean_area = (
        float(cont_track.loc[~cont_track["is_quiet"], "core_pixel_count"].mean()) * cell_area_km2
        if n_active_cont > 0 else 0.0
    )

    scalar_stats = {
        "factual_total_core_volume_mm_km2": hist_total_vol,
        "counterfactual_total_core_volume_mm_km2": cont_total_vol,
        "absolute_diff_mm_km2": abs_diff,
        "relative_diff_pct": rel_diff_pct,
        "factual_mean_core_area_km2": hist_mean_area,
        "counterfactual_mean_core_area_km2": cont_mean_area,
        "factual_mean_centroid_lat": float(hist_track["centroid_lat"].mean(skipna=True)),
        "factual_mean_centroid_lon": float(hist_track["centroid_lon"].mean(skipna=True)),
        "counterfactual_mean_centroid_lat": float(cont_track["centroid_lat"].mean(skipna=True)),
        "counterfactual_mean_centroid_lon": float(cont_track["centroid_lon"].mean(skipna=True)),
        "n_active_timesteps_factual": n_active_hist,
        "n_active_timesteps_counterfactual": n_active_cont,
        "gaussian_sigma_used": sigma,
        "percentile_threshold_used": pct,
    }

    return {
        "factual_core_mask": hist_mask,
        "counterfactual_core_mask": cont_mask,
        "factual_track_df": hist_track,
        "counterfactual_track_df": cont_track,
        "factual_core_volume_ts": hist_vol_ts,
        "counterfactual_core_volume_ts": cont_vol_ts,
        "core_volume_diff_ts": core_volume_diff_ts,
        "core_fraction_hist": core_fraction_hist,
        "core_fraction_cont": core_fraction_cont,
        "scalar_stats": scalar_stats,
    }


# =============================================================================
# Visualization Functions
# =============================================================================


def plot_spatial_comparison(
    results: Dict,
    config: AnalysisConfig,
    ax_list: List = None,
    save_path: Optional[Path] = None,
):
    """Create spatial comparison plots for accumulated/mean values."""
    var_info = results["variable_info"]
    variable = results["variable"]

    if ax_list is None:
        if HAS_CARTOPY:
            fig, axes = plt.subplots(
                1, 3, figsize=(15, 5), subplot_kw={"projection": ccrs.PlateCarree()}
            )
        else:
            fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    else:
        axes = ax_list
        fig = axes[0].get_figure()

    # Get spatial data based on aggregation type
    if var_info["aggregation"] == "sum":
        factual_spatial = results["scenarios"]["hist"]["temporal"][
            "accumulated_spatial"
        ]
        counter_spatial = results["scenarios"]["cont"]["temporal"][
            "accumulated_spatial"
        ]
        diff_spatial = results["comparison"]["accumulated_diff_spatial"]
        title_suffix = "Accumulated"
    else:
        factual_spatial = results["scenarios"]["hist"]["temporal"]["mean_spatial"]
        counter_spatial = results["scenarios"]["cont"]["temporal"]["mean_spatial"]
        diff_spatial = results["comparison"]["mean_diff_spatial"]
        title_suffix = "Mean"

    # Common extent
    lon = factual_spatial.longitude.values
    lat = factual_spatial.latitude.values
    extent = [lon.min(), lon.max(), lat.min(), lat.max()]

    # Plot factual
    vmin = min(factual_spatial.min(), counter_spatial.min())
    vmax = max(factual_spatial.max(), counter_spatial.max())

    plot_kwargs = {"transform": ccrs.PlateCarree()} if HAS_CARTOPY else {}

    im1 = axes[0].pcolormesh(
        lon,
        lat,
        factual_spatial,
        cmap=var_info["colormap"],
        vmin=vmin,
        vmax=vmax,
        **plot_kwargs,
    )
    axes[0].set_title(
        f'{config.scenario_labels["hist"]}\n{title_suffix} {var_info["name"]}'
    )
    if HAS_CARTOPY:
        axes[0].coastlines()
        axes[0].add_feature(cfeature.BORDERS, linestyle=":")
    axes[0].set_xlabel("Longitude")
    axes[0].set_ylabel("Latitude")
    plt.colorbar(im1, ax=axes[0], label=var_info["unit"], shrink=0.8)

    # Plot counterfactual
    im2 = axes[1].pcolormesh(
        lon,
        lat,
        counter_spatial,
        cmap=var_info["colormap"],
        vmin=vmin,
        vmax=vmax,
        **plot_kwargs,
    )
    axes[1].set_title(
        f'{config.scenario_labels["cont"]}\n{title_suffix} {var_info["name"]}'
    )
    if HAS_CARTOPY:
        axes[1].coastlines()
        axes[1].add_feature(cfeature.BORDERS, linestyle=":")
    axes[1].set_xlabel("Longitude")
    plt.colorbar(im2, ax=axes[1], label=var_info["unit"], shrink=0.8)

    # Plot difference
    abs_max = np.abs(diff_spatial).max()
    im3 = axes[2].pcolormesh(
        lon,
        lat,
        diff_spatial,
        cmap=var_info["diff_colormap"],
        vmin=-abs_max,
        vmax=abs_max,
        **plot_kwargs,
    )
    axes[2].set_title(
        f'Difference (Factual - Counterfactual)\n{title_suffix} {var_info["name"]}'
    )
    if HAS_CARTOPY:
        axes[2].coastlines()
        axes[2].add_feature(cfeature.BORDERS, linestyle=":")
    axes[2].set_xlabel("Longitude")
    plt.colorbar(im3, ax=axes[2], label=f"Δ {var_info['unit']}", shrink=0.8)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig, axes


def plot_timeseries_comparison(
    results: Dict, config: AnalysisConfig, save_path: Optional[Path] = None
):
    """Create time series comparison plots."""
    var_info = results["variable_info"]

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    # Get domain-mean time series
    ts_factual = results["scenarios"]["hist"]["temporal"]["domain_mean_timeseries"]
    ts_counter = results["scenarios"]["cont"]["temporal"]["domain_mean_timeseries"]
    ts_diff = results["comparison"]["domain_mean_diff_timeseries"]

    time = ts_factual.time.values

    # Plot absolute values
    axes[0].plot(
        time, ts_factual.values, "b-", linewidth=2, label=config.scenario_labels["hist"]
    )
    axes[0].plot(
        time,
        ts_counter.values,
        "r--",
        linewidth=2,
        label=config.scenario_labels["cont"],
    )
    axes[0].set_ylabel(f'{var_info["name"]} [{var_info["unit"]}]')
    axes[0].set_title(f'Domain-Averaged {var_info["name"]} Time Series')
    axes[0].legend(loc="upper right")
    axes[0].grid(True, alpha=0.3)

    # Fill between to show difference
    axes[0].fill_between(
        time, ts_factual.values, ts_counter.values, alpha=0.3, color="purple"
    )

    # Plot difference
    axes[1].plot(time, ts_diff.values, "k-", linewidth=2)
    axes[1].axhline(y=0, color="gray", linestyle="--", alpha=0.7)
    axes[1].fill_between(
        time,
        0,
        ts_diff.values,
        where=ts_diff.values > 0,
        color="blue",
        alpha=0.3,
        label="Factual higher",
    )
    axes[1].fill_between(
        time,
        0,
        ts_diff.values,
        where=ts_diff.values < 0,
        color="red",
        alpha=0.3,
        label="Counterfactual higher",
    )
    axes[1].set_ylabel(f'Difference [{var_info["unit"]}]')
    axes[1].set_xlabel("Time")
    axes[1].set_title("Difference (Factual - Counterfactual)")
    axes[1].legend(loc="upper right")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig, axes


def plot_extreme_value_comparison(
    results: Dict, config: AnalysisConfig, save_path: Optional[Path] = None
):
    """Create extreme value comparison plots with percentiles."""
    var_info = results["variable_info"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Get percentile data
    percentiles_factual = results["scenarios"]["hist"]["percentiles"]
    percentiles_counter = results["scenarios"]["cont"]["percentiles"]

    percentile_keys = sorted(
        [k for k in percentiles_factual.keys()], key=lambda x: int(x[1:])
    )
    percentile_vals = [int(k[1:]) for k in percentile_keys]

    factual_vals = [percentiles_factual[k] for k in percentile_keys]
    counter_vals = [percentiles_counter[k] for k in percentile_keys]

    # Bar comparison
    x = np.arange(len(percentile_keys))
    width = 0.35

    axes[0].bar(
        x - width / 2,
        factual_vals,
        width,
        label=config.scenario_labels["hist"],
        color="blue",
        alpha=0.7,
    )
    axes[0].bar(
        x + width / 2,
        counter_vals,
        width,
        label=config.scenario_labels["cont"],
        color="red",
        alpha=0.7,
    )
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([f"P{p}" for p in percentile_vals])
    axes[0].set_ylabel(f'{var_info["name"]} [{var_info["unit"]}]')
    axes[0].set_title(f'Percentile Comparison: {var_info["name"]}')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3, axis="y")

    # Difference plot
    diff_vals = np.array(factual_vals) - np.array(counter_vals)
    colors_bar = ["blue" if d > 0 else "red" for d in diff_vals]
    axes[1].bar(x, diff_vals, color=colors_bar, alpha=0.7)
    axes[1].axhline(y=0, color="black", linestyle="-", linewidth=0.5)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([f"P{p}" for p in percentile_vals])
    axes[1].set_ylabel(f'Difference [{var_info["unit"]}]')
    axes[1].set_title("Percentile Differences (Factual - Counterfactual)")
    axes[1].grid(True, alpha=0.3, axis="y")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig, axes


def plot_storm_core_snapshots(
    data_factual: xr.DataArray,
    data_counterfactual: xr.DataArray,
    storm_core_results: Dict,
    config: "AnalysisConfig",
    n_snapshots: int = 4,
    save_path: Optional[Path] = None,
):
    """
    Multi-panel snapshot plot for visual inspection of the detected core cluster.

    For each selected timestep (peak + evenly distributed) shows a row of 4 panels:
      [Factual precip + core contour] [Counterfactual precip + core contour]
      [Core mask factual]             [Core mask counterfactual]

    The core boundary is drawn as a thick red contour and the precipitation-
    weighted centroid is marked with a yellow star, so you can judge at a
    glance whether the cluster detection looks physically reasonable.

    Parameters
    ----------
    data_factual : xr.DataArray
        Precipitation DataArray (time, latitude, longitude) in mm.
    data_counterfactual : xr.DataArray
        Same for counterfactual scenario.
    storm_core_results : dict
        Returned by compute_storm_core_analysis().
    config : AnalysisConfig
    n_snapshots : int
        Number of timesteps to show (default 4). The timestep with the
        highest factual domain-mean precipitation is always included.
    save_path : Path, optional

    Returns
    -------
    fig, axes (n_snapshots x 4 array)
    """
    hist_mask = storm_core_results["factual_core_mask"]
    cont_mask = storm_core_results["counterfactual_core_mask"]
    hist_track = storm_core_results["factual_track_df"]
    cont_track = storm_core_results["counterfactual_track_df"]

    lon = data_factual.longitude.values
    lat = data_factual.latitude.values
    n_times = len(data_factual.time)

    # Choose snapshot timesteps: evenly spaced, but always include the peak
    domain_mean = data_factual.mean(dim=["latitude", "longitude"])
    peak_idx = int(domain_mean.argmax().values)
    candidate_indices = np.linspace(0, n_times - 1, n_snapshots, dtype=int).tolist()
    if peak_idx not in candidate_indices:
        # Replace the candidate closest to peak with the peak itself
        distances = [abs(i - peak_idx) for i in candidate_indices]
        candidate_indices[int(np.argmin(distances))] = peak_idx
    snapshot_indices = sorted(set(candidate_indices))[:n_snapshots]

    # Shared colour scale across all timesteps for precipitation
    vmax_precip = float(max(data_factual.max(), data_counterfactual.max()))
    vmax_precip = max(vmax_precip, 1.0)

    fig, axes = plt.subplots(
        len(snapshot_indices), 4,
        figsize=(18, 4.5 * len(snapshot_indices)),
    )
    # Ensure axes is always 2D even for a single row
    if len(snapshot_indices) == 1:
        axes = axes[np.newaxis, :]

    col_titles = [
        f"{config.scenario_labels['hist']}\nPrecip + Core boundary",
        f"{config.scenario_labels['cont']}\nPrecip + Core boundary",
        f"{config.scenario_labels['hist']}\nCore mask (1 = in core)",
        f"{config.scenario_labels['cont']}\nCore mask (1 = in core)",
    ]
    for col, title in enumerate(col_titles):
        axes[0, col].set_title(title, fontsize=9, fontweight="bold")

    for row, tidx in enumerate(snapshot_indices):
        h_field = data_factual.isel(time=tidx).values
        c_field = data_counterfactual.isel(time=tidx).values
        h_mask_2d = hist_mask.isel(time=tidx).values.astype(float)
        c_mask_2d = cont_mask.isel(time=tidx).values.astype(float)

        time_str = str(data_factual.time.values[tidx])[:16]
        is_peak = (tidx == peak_idx)
        row_label = f"{time_str}{'  ← PEAK' if is_peak else ''}"

        h_row = hist_track.iloc[tidx]
        c_row = cont_track.iloc[tidx]

        # --- Panel 0: factual precip + core contour ---
        ax = axes[row, 0]
        pcm = ax.pcolormesh(lon, lat, h_field, cmap="Blues",
                            vmin=0, vmax=vmax_precip, shading="auto")
        if h_mask_2d.any():
            ax.contour(lon, lat, h_mask_2d, levels=[0.5],
                       colors="red", linewidths=2.0)
        if not h_row["is_quiet"] and not np.isnan(h_row["centroid_lat"]):
            ax.scatter(h_row["centroid_lon"], h_row["centroid_lat"],
                       s=180, marker="*", color="yellow",
                       edgecolors="black", linewidths=0.8, zorder=10)
        ax.set_ylabel(row_label, fontsize=8)
        ax.set_xlabel("Lon" if row == len(snapshot_indices) - 1 else "")
        plt.colorbar(pcm, ax=ax, label="mm", shrink=0.85, pad=0.02)

        # --- Panel 1: counterfactual precip + core contour ---
        ax = axes[row, 1]
        pcm2 = ax.pcolormesh(lon, lat, c_field, cmap="Blues",
                             vmin=0, vmax=vmax_precip, shading="auto")
        if c_mask_2d.any():
            ax.contour(lon, lat, c_mask_2d, levels=[0.5],
                       colors="red", linewidths=2.0)
        if not c_row["is_quiet"] and not np.isnan(c_row["centroid_lat"]):
            ax.scatter(c_row["centroid_lon"], c_row["centroid_lat"],
                       s=180, marker="*", color="yellow",
                       edgecolors="black", linewidths=0.8, zorder=10)
        ax.set_xlabel("Lon" if row == len(snapshot_indices) - 1 else "")
        plt.colorbar(pcm2, ax=ax, label="mm", shrink=0.85, pad=0.02)

        # --- Panel 2: factual core mask (binary) ---
        ax = axes[row, 2]
        ax.pcolormesh(lon, lat, h_mask_2d, cmap="Reds",
                      vmin=0, vmax=1, shading="auto", alpha=0.7)
        # Overlay light precipitation field as contourf for context
        ax.contourf(lon, lat, h_field,
                    levels=np.linspace(0, vmax_precip, 8),
                    cmap="Blues", alpha=0.35)
        if not h_row["is_quiet"] and not np.isnan(h_row["centroid_lat"]):
            ax.scatter(h_row["centroid_lon"], h_row["centroid_lat"],
                       s=180, marker="*", color="gold",
                       edgecolors="black", linewidths=0.8, zorder=10)
        n_px = int(h_row["core_pixel_count"])
        thr = h_row["threshold_used"]
        ax.set_xlabel("Lon" if row == len(snapshot_indices) - 1 else "")
        ax.set_title(f"n_px={n_px}  thr={thr:.1f}mm", fontsize=7, color="darkred",
                     pad=2)

        # --- Panel 3: counterfactual core mask (binary) ---
        ax = axes[row, 3]
        ax.pcolormesh(lon, lat, c_mask_2d, cmap="Reds",
                      vmin=0, vmax=1, shading="auto", alpha=0.7)
        ax.contourf(lon, lat, c_field,
                    levels=np.linspace(0, vmax_precip, 8),
                    cmap="Blues", alpha=0.35)
        if not c_row["is_quiet"] and not np.isnan(c_row["centroid_lat"]):
            ax.scatter(c_row["centroid_lon"], c_row["centroid_lat"],
                       s=180, marker="*", color="gold",
                       edgecolors="black", linewidths=0.8, zorder=10)
        n_px_c = int(c_row["core_pixel_count"])
        thr_c = c_row["threshold_used"]
        ax.set_xlabel("Lon" if row == len(snapshot_indices) - 1 else "")
        ax.set_title(f"n_px={n_px_c}  thr={thr_c:.1f}mm", fontsize=7, color="darkred",
                     pad=2)

    # Common y-axis labels for lat
    for row_axes in axes:
        for ax in row_axes:
            ax.set_ylabel(ax.get_ylabel() or "Lat")

    sc = storm_core_results["scalar_stats"]
    fig.suptitle(
        f"Storm Core Cluster Inspection — {config.region_name}  "
        f"(red contour = core boundary, ★ = precip-weighted centroid)\n"
        f"Core vol attribution: {sc['absolute_diff_mm_km2']:.0f} mm·km²  "
        f"({sc['relative_diff_pct']:.1f}%)",
        fontsize=11, fontweight="bold", y=1.01,
    )
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig, axes


def plot_centroid_tracks(
    storm_core_results: Dict,
    config: "AnalysisConfig",
    save_path: Optional[Path] = None,
):
    """
    Map showing the precipitation-weighted centroid tracks for both scenarios.

    The track is colored by normalized time (blue=early, red=late) and
    connected with a thin line to reveal the storm trajectory. Start and
    end positions are annotated with distinct markers.

    Parameters
    ----------
    storm_core_results : dict
        Returned by compute_storm_core_analysis().
    config : AnalysisConfig
    save_path : Path, optional

    Returns
    -------
    fig, ax
    """
    hist_track = storm_core_results["factual_track_df"]
    cont_track = storm_core_results["counterfactual_track_df"]
    sc = storm_core_results["scalar_stats"]

    # Filter to active (non-quiet) timesteps
    h_active = hist_track[~hist_track["is_quiet"]].reset_index(drop=True)
    c_active = cont_track[~cont_track["is_quiet"]].reset_index(drop=True)

    if HAS_CARTOPY:
        fig, ax = plt.subplots(
            figsize=(10, 8), subplot_kw={"projection": ccrs.PlateCarree()}
        )
    else:
        fig, ax = plt.subplots(figsize=(10, 8))

    plot_kwargs = {"transform": ccrs.PlateCarree()} if HAS_CARTOPY else {}
    cmap_time = plt.cm.coolwarm

    def _plot_track(track_df, label, linestyle, marker_color):
        if len(track_df) == 0:
            return
        norm_time = np.linspace(0, 1, len(track_df))
        sc_plot = ax.scatter(
            track_df["centroid_lon"],
            track_df["centroid_lat"],
            c=norm_time,
            cmap=cmap_time,
            s=60,
            zorder=5,
            label=label,
            edgecolors=marker_color,
            linewidths=0.8,
            **plot_kwargs,
        )
        # Connecting line
        ax.plot(
            track_df["centroid_lon"],
            track_df["centroid_lat"],
            color=marker_color,
            lw=1.0,
            ls=linestyle,
            alpha=0.6,
            zorder=4,
            **plot_kwargs,
        )
        # Start marker
        ax.scatter(
            track_df["centroid_lon"].iloc[0],
            track_df["centroid_lat"].iloc[0],
            s=150, marker="o", color=marker_color, zorder=6, edgecolors="white",
            linewidths=1.5, **plot_kwargs,
        )
        # End marker
        ax.scatter(
            track_df["centroid_lon"].iloc[-1],
            track_df["centroid_lat"].iloc[-1],
            s=200, marker="*", color=marker_color, zorder=6, edgecolors="white",
            linewidths=1.0, **plot_kwargs,
        )
        return sc_plot

    _plot_track(h_active, config.scenario_labels["hist"], "-", "#1f77b4")
    _plot_track(c_active, config.scenario_labels["cont"], "--", "#d62728")

    if HAS_CARTOPY:
        ax.coastlines(resolution="10m", linewidth=0.8)
        ax.add_feature(cfeature.BORDERS, linestyle=":", linewidth=0.6)
        ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False,
                     alpha=0.4, linestyle="--")
    else:
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.grid(True, alpha=0.3)

    ax.legend(loc="upper left", fontsize=9)
    ax.set_title(
        f"Storm Core Centroid Track — {config.region_name}\n"
        f"Circle = start, Star = end; color = time (blue→red)",
        fontsize=11,
    )

    # Stats text box
    stats_str = (
        f"Factual core vol:      {sc['factual_total_core_volume_mm_km2']:.0f} mm·km²\n"
        f"Counterfact core vol:  {sc['counterfactual_total_core_volume_mm_km2']:.0f} mm·km²\n"
        f"Difference:            {sc['absolute_diff_mm_km2']:.0f} mm·km² "
        f"({sc['relative_diff_pct']:.1f}%)\n"
        f"Active timesteps — F:{sc['n_active_timesteps_factual']}  "
        f"CF:{sc['n_active_timesteps_counterfactual']}"
    )
    ax.text(
        0.02, 0.02, stats_str,
        transform=ax.transAxes, fontsize=8, verticalalignment="bottom",
        fontfamily="monospace",
        bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8),
    )

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig, ax


def plot_core_volume_timeseries(
    storm_core_results: Dict,
    config: "AnalysisConfig",
    save_path: Optional[Path] = None,
):
    """
    Two-panel time series of storm core precipitation volume.

    Panel 1: Core volume time series for both scenarios with fill_between.
    Panel 2: Difference (factual - counterfactual) with signed colour fill.

    Parameters
    ----------
    storm_core_results : dict
        Returned by compute_storm_core_analysis().
    config : AnalysisConfig
    save_path : Path, optional

    Returns
    -------
    fig, axes (array of 2)
    """
    hist_vol = storm_core_results["factual_core_volume_ts"]
    cont_vol = storm_core_results["counterfactual_core_volume_ts"]
    diff_vol = storm_core_results["core_volume_diff_ts"]
    sc = storm_core_results["scalar_stats"]

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    time = hist_vol.time.values

    # Panel 1 — absolute volumes
    axes[0].plot(time, hist_vol.values, "b-", linewidth=2,
                 label=config.scenario_labels["hist"])
    axes[0].plot(time, cont_vol.values, "r--", linewidth=2,
                 label=config.scenario_labels["cont"])
    axes[0].fill_between(time, hist_vol.values, cont_vol.values,
                          alpha=0.25, color="purple")
    axes[0].set_ylabel("Core precipitation volume [mm·km²]")
    axes[0].set_title(
        f"Storm Core Precipitation Volume — {config.region_name}\n"
        f"Total: F={sc['factual_total_core_volume_mm_km2']:.0f}  "
        f"CF={sc['counterfactual_total_core_volume_mm_km2']:.0f}  "
        f"Δ={sc['absolute_diff_mm_km2']:.0f} mm·km² ({sc['relative_diff_pct']:.1f}%)"
    )
    axes[0].legend(loc="upper right")
    axes[0].grid(True, alpha=0.3)

    # Panel 2 — difference
    axes[1].plot(time, diff_vol.values, "k-", linewidth=2)
    axes[1].axhline(y=0, color="gray", linestyle="--", alpha=0.7)
    axes[1].fill_between(
        time, 0, diff_vol.values,
        where=diff_vol.values > 0, color="blue", alpha=0.3,
        label="Factual higher",
    )
    axes[1].fill_between(
        time, 0, diff_vol.values,
        where=diff_vol.values < 0, color="red", alpha=0.3,
        label="Counterfactual higher",
    )
    axes[1].set_ylabel("Δ Core volume [mm·km²]")
    axes[1].set_xlabel("Time")
    axes[1].set_title("Difference (Factual − Counterfactual)")
    axes[1].legend(loc="upper right")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig, axes


def plot_time_integrated_core_mask(
    storm_core_results: Dict,
    config: "AnalysisConfig",
    save_path: Optional[Path] = None,
):
    """
    Three-panel spatial map of time-integrated storm core frequency.

    Panels: factual core frequency | counterfactual core frequency | difference.
    Each pixel shows the fraction of active timesteps it was inside the core
    (0 = never, 1 = always).

    Parameters
    ----------
    storm_core_results : dict
        Returned by compute_storm_core_analysis().
    config : AnalysisConfig
    save_path : Path, optional

    Returns
    -------
    fig, axes (array of 3)
    """
    hist_frac = storm_core_results["core_fraction_hist"]
    cont_frac = storm_core_results["core_fraction_cont"]
    diff_frac = hist_frac - cont_frac
    sc = storm_core_results["scalar_stats"]

    if HAS_CARTOPY:
        fig, axes = plt.subplots(
            1, 3, figsize=(16, 5),
            subplot_kw={"projection": ccrs.PlateCarree()},
        )
    else:
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    lon = hist_frac.longitude.values
    lat = hist_frac.latitude.values
    plot_kwargs = {"transform": ccrs.PlateCarree()} if HAS_CARTOPY else {}

    # Panel 1 — factual
    im1 = axes[0].pcolormesh(
        lon, lat, hist_frac, cmap="Blues", vmin=0, vmax=1, **plot_kwargs
    )
    axes[0].set_title(
        f'{config.scenario_labels["hist"]}\nCore Frequency\n'
        f'(n={sc["n_active_timesteps_factual"]} active steps)'
    )
    plt.colorbar(im1, ax=axes[0], label="Fraction of timesteps in core", shrink=0.8)

    # Panel 2 — counterfactual
    im2 = axes[1].pcolormesh(
        lon, lat, cont_frac, cmap="Blues", vmin=0, vmax=1, **plot_kwargs
    )
    axes[1].set_title(
        f'{config.scenario_labels["cont"]}\nCore Frequency\n'
        f'(n={sc["n_active_timesteps_counterfactual"]} active steps)'
    )
    plt.colorbar(im2, ax=axes[1], label="Fraction of timesteps in core", shrink=0.8)

    # Panel 3 — difference
    abs_max = float(np.abs(diff_frac).max())
    abs_max = abs_max if abs_max > 0 else 0.01
    im3 = axes[2].pcolormesh(
        lon, lat, diff_frac, cmap="RdBu_r", vmin=-abs_max, vmax=abs_max,
        **plot_kwargs,
    )
    axes[2].set_title(
        "Difference (Factual − Counterfactual)\nCore Frequency"
    )
    plt.colorbar(im3, ax=axes[2], label="Δ Fraction in core", shrink=0.8)

    for ax in axes:
        if HAS_CARTOPY:
            ax.coastlines(resolution="10m", linewidth=0.8)
            ax.add_feature(cfeature.BORDERS, linestyle=":", linewidth=0.6)
        else:
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")

    fig.suptitle(
        f"Time-Integrated Storm Core Mask — {config.region_name}\n"
        f"Core vol attribution: {sc['absolute_diff_mm_km2']:.0f} mm·km² "
        f"({sc['relative_diff_pct']:.1f}%)",
        fontsize=13, fontweight="bold",
    )
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig, axes


def plot_storm_core_animation(
    data_factual: xr.DataArray,
    data_counterfactual: xr.DataArray,
    storm_core_results: Dict,
    config: "AnalysisConfig",
    save_path: Optional[Path] = None,
    n_frames_max: int = 72,
):
    """
    Animated GIF of the evolving storm core mask at each timestep.

    Each frame shows a 2-panel map (factual | counterfactual) with:
    - Precipitation field as pcolormesh background (Blues)
    - Core mask boundary as a thick red contour
    - Centroid marked with a yellow star

    This function is opt-in and NOT called by the default pipeline.
    Pass a .gif save_path to produce output.

    Parameters
    ----------
    data_factual, data_counterfactual : xr.DataArray
        Full precipitation DataArrays (time, latitude, longitude).
    storm_core_results : dict
        Returned by compute_storm_core_analysis().
    config : AnalysisConfig
    save_path : Path, optional
        Must end in .gif or .mp4.
    n_frames_max : int
        Maximum number of animation frames (subsampled if needed).
    """
    from matplotlib.animation import FuncAnimation

    hist_mask = storm_core_results["factual_core_mask"]
    cont_mask = storm_core_results["counterfactual_core_mask"]
    hist_track = storm_core_results["factual_track_df"]
    cont_track = storm_core_results["counterfactual_track_df"]

    n_times = len(data_factual.time)
    if n_times > n_frames_max:
        frame_indices = np.linspace(0, n_times - 1, n_frames_max, dtype=int)
    else:
        frame_indices = np.arange(n_times)

    lon = data_factual.longitude.values
    lat = data_factual.latitude.values

    vmax_precip = float(
        max(data_factual.max(), data_counterfactual.max()) * 0.8
    )
    vmax_precip = max(vmax_precip, 1.0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    def _init_panel(ax, title):
        pcm = ax.pcolormesh(lon, lat, np.zeros((len(lat), len(lon))),
                            cmap="Blues", vmin=0, vmax=vmax_precip)
        ax.set_title(title)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        star = ax.scatter([], [], s=200, marker="*", color="yellow",
                          edgecolors="black", linewidths=0.8, zorder=10)
        return pcm, star

    pcm1, star1 = _init_panel(ax1, config.scenario_labels["hist"])
    pcm2, star2 = _init_panel(ax2, config.scenario_labels["cont"])
    contour_artists = []
    title_text = fig.suptitle("", fontsize=11)

    def update(frame_idx):
        # Remove old contours
        for coll in contour_artists:
            for c in coll.collections if hasattr(coll, "collections") else [coll]:
                c.remove()
        contour_artists.clear()

        i = frame_idx
        h_field = data_factual.isel(time=i).values
        c_field = data_counterfactual.isel(time=i).values
        h_mask_2d = hist_mask.isel(time=i).values.astype(float)
        c_mask_2d = cont_mask.isel(time=i).values.astype(float)

        pcm1.set_array(h_field.ravel())
        pcm2.set_array(c_field.ravel())

        # Core mask contour
        if h_mask_2d.any():
            cs1 = ax1.contour(lon, lat, h_mask_2d, levels=[0.5], colors="red", linewidths=2)
            contour_artists.append(cs1)
        if c_mask_2d.any():
            cs2 = ax2.contour(lon, lat, c_mask_2d, levels=[0.5], colors="red", linewidths=2)
            contour_artists.append(cs2)

        # Centroid
        clat = hist_track["centroid_lat"].iloc[i]
        clon = hist_track["centroid_lon"].iloc[i]
        if not np.isnan(clat):
            star1.set_offsets([[clon, clat]])
        else:
            star1.set_offsets(np.empty((0, 2)))

        clat2 = cont_track["centroid_lat"].iloc[i]
        clon2 = cont_track["centroid_lon"].iloc[i]
        if not np.isnan(clat2):
            star2.set_offsets([[clon2, clat2]])
        else:
            star2.set_offsets(np.empty((0, 2)))

        time_str = str(data_factual.time.values[i])[:16]
        title_text.set_text(f"Storm Core — {config.region_name}  |  {time_str}")

    anim = FuncAnimation(fig, update, frames=frame_indices, interval=200, blit=False)

    if save_path:
        writer = "pillow" if str(save_path).endswith(".gif") else "ffmpeg"
        anim.save(str(save_path), writer=writer, fps=6, dpi=100)
        print(f"Saved animation: {save_path}")
        plt.close(fig)
    else:
        plt.show()


def plot_summary_dashboard(
    all_results: Dict, config: AnalysisConfig, save_path: Optional[Path] = None
):
    """Create a comprehensive summary dashboard for all variables."""
    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(3, 4, figure=fig, hspace=0.35, wspace=0.3)

    variables = list(all_results.keys())

    for i, variable in enumerate(variables):
        results = all_results[variable]
        var_info = results["variable_info"]

        # Spatial difference map
        if HAS_CARTOPY:
            ax_map = fig.add_subplot(gs[i, 0], projection=ccrs.PlateCarree())
        else:
            ax_map = fig.add_subplot(gs[i, 0])

        if var_info["aggregation"] == "sum":
            diff_spatial = results["comparison"]["accumulated_diff_spatial"]
            title = f"{var_info['name']}\nAccum. Difference"
        else:
            diff_spatial = results["comparison"]["mean_diff_spatial"]
            title = f"{var_info['name']}\nMean Difference"

        lon = diff_spatial.longitude.values
        lat = diff_spatial.latitude.values
        abs_max = np.abs(diff_spatial).max()

        dash_plot_kwargs = {"transform": ccrs.PlateCarree()} if HAS_CARTOPY else {}
        im = ax_map.pcolormesh(
            lon,
            lat,
            diff_spatial,
            cmap=var_info["diff_colormap"],
            vmin=-abs_max,
            vmax=abs_max,
            **dash_plot_kwargs,
        )
        if HAS_CARTOPY:
            ax_map.coastlines()
        ax_map.set_title(title, fontsize=10)
        plt.colorbar(im, ax=ax_map, label=f"Δ{var_info['unit']}", shrink=0.7)

        # Time series comparison
        ax_ts = fig.add_subplot(gs[i, 1:3])
        ts_factual = results["scenarios"]["hist"]["temporal"]["domain_mean_timeseries"]
        ts_counter = results["scenarios"]["cont"]["temporal"]["domain_mean_timeseries"]

        ax_ts.plot(
            ts_factual.time.values,
            ts_factual.values,
            "b-",
            linewidth=1.5,
            label="Factual",
        )
        ax_ts.plot(
            ts_counter.time.values,
            ts_counter.values,
            "r--",
            linewidth=1.5,
            label="Counterfactual",
        )
        ax_ts.set_ylabel(var_info["unit"])
        ax_ts.set_title(f'{var_info["name"]} Time Series', fontsize=10)
        ax_ts.legend(loc="upper right", fontsize=8)
        ax_ts.grid(True, alpha=0.3)
        ax_ts.tick_params(axis="x", rotation=45)

        # Statistics summary
        ax_stats = fig.add_subplot(gs[i, 3])
        ax_stats.axis("off")

        stats_text = f"**{var_info['name']}**\n"
        stats_text += "-" * 25 + "\n"

        fact_stats = results["scenarios"]["hist"]["temporal"]
        cont_stats = results["scenarios"]["cont"]["temporal"]
        comp = results["comparison"]

        if var_info["aggregation"] == "sum":
            stats_text += f"Factual Accum: {fact_stats['accumulated_total']:.1f} {var_info['unit']}\n"
            stats_text += f"Counterfact Accum: {cont_stats['accumulated_total']:.1f} {var_info['unit']}\n"
            stats_text += (
                f"Difference: {comp['accumulated_diff_mean']:.1f} {var_info['unit']}\n"
            )
            stats_text += f"Rel. Change: {comp['accumulated_rel_diff_mean']:.1f}%\n"
        else:
            stats_text += (
                f"Factual Mean: {fact_stats['period_mean']:.2f} {var_info['unit']}\n"
            )
            stats_text += f"Counterfact Mean: {cont_stats['period_mean']:.2f} {var_info['unit']}\n"
            stats_text += (
                f"Difference: {comp['mean_diff_mean']:.2f} {var_info['unit']}\n"
            )

        stats_text += "-" * 25 + "\n"
        stats_text += f"Factual Max: {fact_stats['max_value']:.2f}\n"
        stats_text += f"Counterfact Max: {cont_stats['max_value']:.2f}\n"
        stats_text += f"Max Difference: {comp['max_diff']:.2f}\n"
        stats_text += "-" * 25 + "\n"
        stats_text += f"Factual Min: {fact_stats['min_value']:.2f}\n"
        stats_text += f"Counterfact Min: {cont_stats['min_value']:.2f}\n"
        stats_text += f"Min Difference: {comp['min_diff']:.2f}\n"

        ax_stats.text(
            0.05,
            0.95,
            stats_text,
            transform=ax_stats.transAxes,
            fontsize=9,
            verticalalignment="top",
            fontfamily="monospace",
            bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.5),
        )

    fig.suptitle(
        f"ClimateDT Climate Attribution Analysis: {config.region_name}\n"
        f"Factual vs Counterfactual Comparison",
        fontsize=14,
        fontweight="bold",
    )

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig


# =============================================================================
# Report Generation
# =============================================================================


def generate_text_report(
    all_results: Dict, config: AnalysisConfig, save_path: Optional[Path] = None
) -> str:
    """Generate a comprehensive text report of the analysis."""
    lines = []
    lines.append("=" * 80)
    lines.append("CLIMATEDT CLIMATE ATTRIBUTION ANALYSIS REPORT")
    lines.append("=" * 80)
    lines.append(f"\nRegion: {config.region_name}")
    lines.append(f"Date Range: {config.date_range}")
    lines.append(f"Data Directory: {config.data_dir}")
    lines.append(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("\n")

    for variable, results in all_results.items():
        var_info = results["variable_info"]
        fact_stats = results["scenarios"]["hist"]["temporal"]
        cont_stats = results["scenarios"]["cont"]["temporal"]
        fact_pct = results["scenarios"]["hist"]["percentiles"]
        cont_pct = results["scenarios"]["cont"]["percentiles"]
        comp = results["comparison"]

        lines.append("-" * 80)
        lines.append(f"\n{var_info['name'].upper()} ({variable})")
        lines.append("-" * 40)

        lines.append("\n### SUMMARY STATISTICS ###")
        lines.append(
            f"\n{'Metric':<35} {'Factual':>15} {'Counterfact':>15} {'Difference':>15}"
        )
        lines.append("-" * 80)

        if var_info["aggregation"] == "sum":
            lines.append(
                f"{'Accumulated Total (domain mean)':<35} "
                f"{fact_stats['accumulated_total']:>15.2f} "
                f"{cont_stats['accumulated_total']:>15.2f} "
                f"{comp['accumulated_diff_mean']:>15.2f}"
            )
            lines.append(
                f"{'Max Accumulated (point)':<35} "
                f"{fact_stats['max_accumulated_point']:>15.2f} "
                f"{cont_stats['max_accumulated_point']:>15.2f} "
                f"{comp['accumulated_diff_max']:>15.2f}"
            )
            lines.append(
                f"{'Relative Change (%)':<35} "
                f"{'':>15} "
                f"{'':>15} "
                f"{comp['accumulated_rel_diff_mean']:>15.1f}%"
            )
        else:
            lines.append(
                f"{'Period Mean':<35} "
                f"{fact_stats['period_mean']:>15.2f} "
                f"{cont_stats['period_mean']:>15.2f} "
                f"{comp['mean_diff_mean']:>15.2f}"
            )

        lines.append(
            f"{'Maximum Value':<35} "
            f"{fact_stats['max_value']:>15.2f} "
            f"{cont_stats['max_value']:>15.2f} "
            f"{comp['max_diff']:>15.2f}"
        )
        lines.append(
            f"{'Minimum Value':<35} "
            f"{fact_stats['min_value']:>15.2f} "
            f"{cont_stats['min_value']:>15.2f} "
            f"{comp['min_diff']:>15.2f}"
        )

        lines.append(
            f"{'Domain Mean - Max':<35} "
            f"{fact_stats['domain_mean_max']:>15.2f} "
            f"{cont_stats['domain_mean_max']:>15.2f} "
            f"{comp['domain_mean_diff_max']:>15.2f}"
        )
        lines.append(
            f"{'Domain Mean - Min':<35} "
            f"{fact_stats['domain_mean_min']:>15.2f} "
            f"{cont_stats['domain_mean_min']:>15.2f} "
            f"{comp['domain_mean_diff_min']:>15.2f}"
        )
        lines.append(
            f"{'Domain Mean - Std Dev':<35} "
            f"{fact_stats['domain_mean_std']:>15.2f} "
            f"{cont_stats['domain_mean_std']:>15.2f} "
            f"{'':>15}"
        )

        lines.append("\n### EXTREME VALUE ANALYSIS (Percentiles) ###")
        lines.append(
            f"\n{'Percentile':<15} {'Factual':>15} {'Counterfact':>15} {'Difference':>15}"
        )
        lines.append("-" * 60)

        for key in sorted(fact_pct.keys(), key=lambda x: int(x[1:])):
            pct_label = key.upper()
            diff = fact_pct[key] - cont_pct[key]
            lines.append(
                f"{pct_label:<15} "
                f"{fact_pct[key]:>15.2f} "
                f"{cont_pct[key]:>15.2f} "
                f"{diff:>15.2f}"
            )

        lines.append("\n### TIME OF EXTREMES ###")
        lines.append(f"Factual Maximum at: {fact_stats['max_time']}")
        lines.append(f"Counterfactual Maximum at: {cont_stats['max_time']}")
        lines.append(f"Factual Minimum at: {fact_stats['min_time']}")
        lines.append(f"Counterfactual Minimum at: {cont_stats['min_time']}")

        # Storm core section (tp only)
        if variable == "tp" and results.get("storm_core") is not None:
            sc = results["storm_core"]["scalar_stats"]
            lines.append("\n### DYNAMIC STORM CORE PRECIPITATION ANALYSIS ###")
            lines.append(
                f"  Method: Gaussian-smoothed (σ={sc['gaussian_sigma_used']}) "
                f"connected-component tracking, P{sc['percentile_threshold_used']:.0f} threshold"
            )
            lines.append(
                f"  Core mask: independent per scenario (Option B — each storm on its own terms)"
            )
            lines.append(
                f"\n{'Metric':<45} {'Factual':>15} {'Counterfact':>15} {'Difference':>15}"
            )
            lines.append("-" * 90)
            lines.append(
                f"{'Total core volume (mm·km²)':<45} "
                f"{sc['factual_total_core_volume_mm_km2']:>15.0f} "
                f"{sc['counterfactual_total_core_volume_mm_km2']:>15.0f} "
                f"{sc['absolute_diff_mm_km2']:>15.0f}"
            )
            lines.append(
                f"{'Relative change (%)':<45} "
                f"{'':>15} {'':>15} "
                f"{sc['relative_diff_pct']:>14.1f}%"
            )
            lines.append(
                f"{'Mean core footprint area (km²)':<45} "
                f"{sc['factual_mean_core_area_km2']:>15.0f} "
                f"{sc['counterfactual_mean_core_area_km2']:>15.0f}"
            )
            lines.append(
                f"{'Active (non-quiet) timesteps':<45} "
                f"{sc['n_active_timesteps_factual']:>15} "
                f"{sc['n_active_timesteps_counterfactual']:>15}"
            )
            lines.append(
                f"{'Mean centroid latitude (°)':<45} "
                f"{sc['factual_mean_centroid_lat']:>15.3f} "
                f"{sc['counterfactual_mean_centroid_lat']:>15.3f} "
                f"{sc['factual_mean_centroid_lat'] - sc['counterfactual_mean_centroid_lat']:>15.3f}"
            )
            lines.append(
                f"{'Mean centroid longitude (°)':<45} "
                f"{sc['factual_mean_centroid_lon']:>15.3f} "
                f"{sc['counterfactual_mean_centroid_lon']:>15.3f} "
                f"{sc['factual_mean_centroid_lon'] - sc['counterfactual_mean_centroid_lon']:>15.3f}"
            )
            lines.append(
                f"\n  Note: 1 mm·km² = 1000 m³ of water. Core volume represents only the "
                f"dominant precipitation cluster (largest connected component) at each timestep."
            )

        # Three-method attribution sub-report (tp only)
        if variable == "tp" and results.get("storm_attribution") is not None:
            attr = results["storm_attribution"]
            pw = attr["peak_window"]["scalar_stats"]
            od = attr["object_decomp"]["scalar_stats"]
            sr = attr["storm_relative"]["scalar_stats"]

            lines.append("\n### METHOD 1 — PEAK ACCUMULATION WINDOW ###")
            lines.append(f"  Window duration: {attr['peak_window']['window_hours']:.0f} h (independent per scenario)")
            lines.append(
                f"\n{'Metric':<45} {'Factual':>15} {'Counterfact':>15} {'Difference':>15}"
            )
            lines.append("-" * 90)
            lines.append(
                f"{'Domain-mean accumulation (mm)':<45} "
                f"{pw['factual_window_accum_mm']:>15.1f} "
                f"{pw['counterfactual_window_accum_mm']:>15.1f} "
                f"{pw['absolute_diff_mm']:>15.1f}"
            )
            lines.append(
                f"{'Relative change (%)':<45} {'':>15} {'':>15} "
                f"{pw['relative_diff_pct']:>14.1f}%"
            )
            lines.append(f"  Factual peak window:        {pw['factual_window_start']} → {pw['factual_window_end']}")
            lines.append(f"  Counterfactual peak window: {pw['counterfactual_window_start']} → {pw['counterfactual_window_end']}")
            sign = "earlier" if pw["window_offset_hours"] > 0 else "later"
            lines.append(f"  Timing offset: {abs(pw['window_offset_hours']):.1f} h  (CF peaks {sign} than factual)")

            lines.append("\n### METHOD 2 — OBJECT DECOMPOSITION (Intensity × Area × Volume) ###")
            lines.append(
                f"\n{'Property':<35} {'Factual':>12} {'Counterfact':>12} {'Δ':>10} {'Δ%':>8}"
            )
            lines.append("-" * 77)
            lines.append(
                f"{'Intensity — P95 in core (mm)':<35} "
                f"{od['mean_intensity_hist_mm']:>12.2f} "
                f"{od['mean_intensity_cont_mm']:>12.2f} "
                f"{od['intensity_diff_mm']:>10.2f} "
                f"{od['intensity_rel_diff_pct']:>7.1f}%"
            )
            lines.append(
                f"{'Area — mean core footprint (km²)':<35} "
                f"{od['mean_area_hist_km2']:>12.0f} "
                f"{od['mean_area_cont_km2']:>12.0f} "
                f"{od['area_diff_km2']:>10.0f} "
                f"{od['area_rel_diff_pct']:>7.1f}%"
            )
            lines.append(
                f"{'Volume — total core (mm·km²)':<35} "
                f"{od['volume_factual_mm_km2']:>12.0f} "
                f"{od['volume_cont_mm_km2']:>12.0f} "
                f"{od['volume_factual_mm_km2'] - od['volume_cont_mm_km2']:>10.0f} "
                f"{od['volume_rel_diff_pct']:>7.1f}%"
            )

            lines.append("\n### METHOD 3 — STORM-RELATIVE COMPOSITE ###")
            lines.append(f"  Composite window: {sr['peak_window_str']}")
            lines.append(f"  Centroid displacement removed: {sr['displacement_km']:.1f} km")
            lines.append(f"    (Δlat={sr['displacement_lat_deg']:+.3f}°, Δlon={sr['displacement_lon_deg']:+.3f}°)")
            lines.append(f"  Storm-relative mean Δ: {sr['storm_relative_diff_mean_mm']:.2f} mm")
            lines.append(f"  Storm-relative max Δ:  {sr['storm_relative_diff_max_mm']:.2f} mm")
            lines.append(f"  Storm-relative min Δ:  {sr['storm_relative_diff_min_mm']:.2f} mm")

        lines.append("\n")

    lines.append("=" * 80)
    lines.append("END OF REPORT")
    lines.append("=" * 80)

    report = "\n".join(lines)

    if save_path:
        with open(save_path, "w") as f:
            f.write(report)
        print(f"Report saved: {save_path}")

    return report


def export_statistics_csv(
    all_results: Dict, config: AnalysisConfig, save_path: Optional[Path] = None
) -> pd.DataFrame:
    """Export key statistics to a CSV file."""
    records = []

    for variable, results in all_results.items():
        var_info = results["variable_info"]
        fact_stats = results["scenarios"]["hist"]["temporal"]
        cont_stats = results["scenarios"]["cont"]["temporal"]
        comp = results["comparison"]

        record = {
            "variable": variable,
            "variable_name": var_info["name"],
            "unit": var_info["unit"],
            "factual_max": fact_stats["max_value"],
            "counterfactual_max": cont_stats["max_value"],
            "max_diff": comp["max_diff"],
            "factual_min": fact_stats["min_value"],
            "counterfactual_min": cont_stats["min_value"],
            "min_diff": comp["min_diff"],
            "factual_domain_mean": fact_stats["domain_mean_mean"],
            "counterfactual_domain_mean": cont_stats["domain_mean_mean"],
        }

        if var_info["aggregation"] == "sum":
            record["factual_accumulated"] = fact_stats["accumulated_total"]
            record["counterfactual_accumulated"] = cont_stats["accumulated_total"]
            record["accumulated_diff"] = comp["accumulated_diff_mean"]
            record["accumulated_rel_diff_pct"] = comp["accumulated_rel_diff_mean"]
        else:
            record["factual_period_mean"] = fact_stats["period_mean"]
            record["counterfactual_period_mean"] = cont_stats["period_mean"]
            record["period_mean_diff"] = comp["mean_diff_mean"]

        # Storm core columns (tp only)
        if variable == "tp" and results.get("storm_core") is not None:
            sc = results["storm_core"]["scalar_stats"]
            record["factual_core_volume_mm_km2"] = sc["factual_total_core_volume_mm_km2"]
            record["counterfactual_core_volume_mm_km2"] = sc["counterfactual_total_core_volume_mm_km2"]
            record["core_volume_diff_mm_km2"] = sc["absolute_diff_mm_km2"]
            record["core_volume_rel_diff_pct"] = sc["relative_diff_pct"]
            record["factual_mean_core_area_km2"] = sc["factual_mean_core_area_km2"]
            record["counterfactual_mean_core_area_km2"] = sc["counterfactual_mean_core_area_km2"]
            record["factual_centroid_lat"] = sc["factual_mean_centroid_lat"]
            record["factual_centroid_lon"] = sc["factual_mean_centroid_lon"]
            record["counterfactual_centroid_lat"] = sc["counterfactual_mean_centroid_lat"]
            record["counterfactual_centroid_lon"] = sc["counterfactual_mean_centroid_lon"]

        # Three-method attribution columns (tp only)
        if variable == "tp" and results.get("storm_attribution") is not None:
            attr = results["storm_attribution"]
            pw = attr["peak_window"]["scalar_stats"]
            od = attr["object_decomp"]["scalar_stats"]
            sr = attr["storm_relative"]["scalar_stats"]
            record["factual_peak_window_start"] = pw["factual_window_start"]
            record["factual_peak_window_end"] = pw["factual_window_end"]
            record["cont_peak_window_start"] = pw["counterfactual_window_start"]
            record["cont_peak_window_end"] = pw["counterfactual_window_end"]
            record["window_offset_hours"] = pw["window_offset_hours"]
            record["peak_window_accum_diff_mm"] = pw["absolute_diff_mm"]
            record["peak_window_accum_rel_diff_pct"] = pw["relative_diff_pct"]
            record["intensity_diff_pct"] = od["intensity_rel_diff_pct"]
            record["area_diff_pct"] = od["area_rel_diff_pct"]
            record["centroid_displacement_km"] = sr["displacement_km"]
            record["storm_relative_diff_mean_mm"] = sr["storm_relative_diff_mean_mm"]

        records.append(record)

    df = pd.DataFrame(records)

    if save_path:
        df.to_csv(save_path, index=False)
        print(f"Statistics exported: {save_path}")

    return df


# =============================================================================
# Main Analysis Function
# =============================================================================


def run_attribution_analysis(
    data_dir: str,
    output_dir: str,
    region_name: str = "Durban",
    date_range: str = "20220409_to_20220414",
    save_figures: bool = True,
    show_figures: bool = False,
) -> Dict:
    """
    Run the complete climate attribution analysis.

    Parameters
    ----------
    data_dir : str
        Path to directory containing ClimateDT NetCDF files
    output_dir : str
        Path to directory for saving outputs
    region_name : str
        Region name (used in filenames and labels)
    date_range : str
        Date range string as in filenames
    save_figures : bool
        Whether to save figures to disk
    show_figures : bool
        Whether to display figures interactively

    Returns
    -------
    Dict
        Dictionary containing all analysis results
    """
    # Initialize configuration
    config = AnalysisConfig(
        data_dir=Path(data_dir),
        output_dir=Path(output_dir),
        region_name=region_name,
        date_range=date_range,
    )

    # Create output directory if needed
    config.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'#'*60}")
    print("CLIMATEDT CLIMATE ATTRIBUTION ANALYSIS")
    print(f"{'#'*60}")
    print(f"\nData directory: {config.data_dir}")
    print(f"Output directory: {config.output_dir}")
    print(f"Region: {config.region_name}")
    print(f"Date range: {config.date_range}")

    # Check which variables are available
    available_vars = []
    for var in config.variables:
        filepath = get_file_path(config, var, "hist")
        if filepath.exists():
            available_vars.append(var)
            print(f"Found: {var}")
        else:
            print(f"Not found: {var} ({filepath})")

    if not available_vars:
        raise FileNotFoundError("No valid variable files found in the data directory.")

    # Run analysis for each variable
    all_results = {}
    for variable in available_vars:
        results = run_full_analysis(config, variable)
        all_results[variable] = results

        # Generate individual variable plots
        if save_figures or show_figures:
            # Spatial comparison
            fig, _ = plot_spatial_comparison(
                results,
                config,
                save_path=(
                    config.output_dir / f"spatial_{variable}_{region_name}.png"
                    if save_figures
                    else None
                ),
            )
            if show_figures:
                plt.show()
            else:
                plt.close(fig)

            # Time series comparison
            fig, _ = plot_timeseries_comparison(
                results,
                config,
                save_path=(
                    config.output_dir / f"timeseries_{variable}_{region_name}.png"
                    if save_figures
                    else None
                ),
            )
            if show_figures:
                plt.show()
            else:
                plt.close(fig)

            # Extreme value comparison
            fig, _ = plot_extreme_value_comparison(
                results,
                config,
                save_path=(
                    config.output_dir / f"extremes_{variable}_{region_name}.png"
                    if save_figures
                    else None
                ),
            )
            if show_figures:
                plt.show()
            else:
                plt.close(fig)

    # Storm core figures (precipitation only)
    sc_results = all_results.get("tp", {}).get("storm_core")
    if sc_results is not None and (save_figures or show_figures):
        print("\nGenerating storm core figures...")
        sc_prefix = config.output_dir / f"storm_core_{config.region_name}"
        tp_data = all_results["tp"]["data"]

        fig, _ = plot_storm_core_snapshots(
            tp_data["hist"], tp_data["cont"], sc_results, config,
            n_snapshots=4,
            save_path=Path(f"{sc_prefix}_cluster_snapshots.png") if save_figures else None,
        )
        if show_figures:
            plt.show()
        else:
            plt.close(fig)

        fig, _ = plot_centroid_tracks(
            sc_results, config,
            save_path=Path(f"{sc_prefix}_centroid_tracks.png") if save_figures else None,
        )
        if show_figures:
            plt.show()
        else:
            plt.close(fig)

        fig, _ = plot_core_volume_timeseries(
            sc_results, config,
            save_path=Path(f"{sc_prefix}_volume_timeseries.png") if save_figures else None,
        )
        if show_figures:
            plt.show()
        else:
            plt.close(fig)

        fig, _ = plot_time_integrated_core_mask(
            sc_results, config,
            save_path=Path(f"{sc_prefix}_integrated_mask.png") if save_figures else None,
        )
        if show_figures:
            plt.show()
        else:
            plt.close(fig)

    # Lagrangian attribution dashboard (all three methods in one figure)
    lagrangian = all_results.get("tp", {}).get("storm_attribution")
    if lagrangian is not None and (save_figures or show_figures):
        print("\nGenerating Lagrangian attribution dashboard...")
        dashboard_path = (
            config.output_dir / f"storm_attr_{config.region_name}_lagrangian_dashboard.png"
            if save_figures
            else None
        )
        fig = plot_lagrangian_dashboard(lagrangian, config, save_path=dashboard_path)
        if show_figures:
            plt.show()
        else:
            plt.close(fig)

    # Generate summary dashboard
    if save_figures or show_figures:
        fig = plot_summary_dashboard(
            all_results,
            config,
            save_path=(
                config.output_dir / f"summary_dashboard_{region_name}.png"
                if save_figures
                else None
            ),
        )
        if show_figures:
            plt.show()
        else:
            plt.close(fig)

    # Generate text report
    report = generate_text_report(
        all_results,
        config,
        save_path=config.output_dir / f"attribution_report_{region_name}.txt",
    )
    print("\n" + report)

    # Export statistics to CSV
    export_statistics_csv(
        all_results,
        config,
        save_path=config.output_dir / f"attribution_statistics_{region_name}.csv",
    )

    print(f"\n{'#'*60}")
    print("ANALYSIS COMPLETE")
    print(f"{'#'*60}")
    print(f"\nOutputs saved to: {config.output_dir}")

    return all_results


# =============================================================================
# Command Line Interface
# =============================================================================


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="ClimateDT Climate Attribution Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze larger area data
  python climateDT_attribution_analysis.py \\
      --data_dir /p/11210471-001-compass/01_Data/ECMWF_ClimateDT/data/nc_Durban \\
      --output_dir ./output/larger_area

  # Analyze smaller area data
  python climateDT_attribution_analysis.py \\
      --data_dir /p/11210471-001-compass/01_Data/ECMWF_ClimateDT/data/nc_Durban/smaller_area \\
      --output_dir ./output/smaller_area

  # Custom region name and date range
  python climateDT_attribution_analysis.py \\
      --data_dir /path/to/data \\
      --output_dir /path/to/output \\
      --region_name "MyRegion" \\
      --date_range "20220409_to_20220414"
        """,
    )

    parser.add_argument(
        "--data_dir",
        "-d",
        type=str,
        required=True,
        help="Directory containing ClimateDT NetCDF files",
    )

    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        required=True,
        help="Directory for saving output files",
    )

    parser.add_argument(
        "--region_name",
        "-r",
        type=str,
        default="Durban",
        help="Region name for labeling (default: Durban)",
    )

    parser.add_argument(
        "--date_range",
        type=str,
        default="20220409_to_20220414",
        help="Date range string as in filenames (default: 20220409_to_20220414)",
    )

    parser.add_argument(
        "--no_save_figures", action="store_true", help="Do not save figures to disk"
    )

    parser.add_argument(
        "--show_figures", action="store_true", help="Display figures interactively"
    )

    return parser.parse_args()


# =============================================================================
# Main Entry Point
# =============================================================================

if __name__ == "__main__":
    # Check if running interactively (e.g., in Jupyter or IDE)
    import sys

    if len(sys.argv) > 1 and sys.argv[1] != "-f":  # -f is used by Jupyter
        # Command line execution
        args = parse_arguments()

        results = run_attribution_analysis(
            data_dir=args.data_dir,
            output_dir=args.output_dir,
            region_name=args.region_name,
            date_range=args.date_range,
            save_figures=not args.no_save_figures,
            show_figures=args.show_figures,
        )
    else:
        # Interactive/testing execution with default values
        print("Running in interactive mode with default settings...")
        print(
            "For command-line usage, run: python climateDT_attribution_analysis.py --help"
        )

        # Default paths for testing - adjust these as needed
        DEFAULT_DATA_DIR_LARGE = (
            "/p/11210471-001-compass/01_Data/ECMWF_ClimateDT/data/nc_Durban"
        )
        DEFAULT_DATA_DIR_SMALL = "/p/11210471-001-compass/01_Data/ECMWF_ClimateDT/data/nc_Durban/smaller_area"
        DEFAULT_OUTPUT_DIR = (
            "/p/11210471-001-compass/01_Data/ECMWF_ClimateDT/analysis_output"
        )

        # Example: Run analysis for larger area
        # results = run_attribution_analysis(
        #     data_dir=DEFAULT_DATA_DIR_LARGE,
        #     output_dir=f"{DEFAULT_OUTPUT_DIR}/larger_area",
        #     region_name="Durban",
        #     save_figures=True,
        #     show_figures=True
        # )

        pass

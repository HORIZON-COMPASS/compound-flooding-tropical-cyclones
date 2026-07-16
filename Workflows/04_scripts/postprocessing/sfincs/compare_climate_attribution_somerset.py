"""
Climate Attribution Analysis: Somerset Dec 2013 — Factual vs T1 (past) vs T3 (future).

Compares max flood depth outputs across three CEH-GEAR rainfall scenarios:
- Factual:     ceh_gear_compass  — observed hourly rainfall
- T1 (past):   ceh_t1_compass    — pre-industrial / past climate (less intense)
- T3 (future): ceh_t3_compass    — future climate signal (more intense)

Produces:
1. 3-panel hmax map: Factual | T1 | T3
2. Side-by-side attribution difference maps (Factual−T1, T3−Factual)
3. Bar chart for flood metrics (volume, extent, mean depth) with attribution brackets
4. Summary CSV

Adapted from compare_climate_attribution.py (Durban 2022).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import rioxarray as rxr
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")


# ===== CONFIGURATION =====
EVENT_NAME    = "SomersetLevels_Dec2013"
EVENT_DISPLAY = "Somerset Levels Dec 2013"
REGION        = "somerset"

BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs")
OUTPUT_DIR    = Path("/p/11210471-001-compass/04_Results/climate_attribution_somerset")

# Each entry: run name → scenario folder (under BASE_RUN_PATH/REGION/<scenario>/sfincs/)
RUNS = {
    "Factual": {
        "scenario": "SomersetLevels_dec_factual",
        "folder":   "event_tp_ceh_gear_compass_CF0_GTSMv41opendap_CF0_no_wind_CF0",
        "color":    "#4472C4",
        "label":    "Factual\n(observed)",
    },
    "T1 (past)": {
        "scenario": "SomersetLevels_dec_t1",
        "folder":   "event_tp_ceh_t1_compass_CF0_GTSMv41opendap_CF0_no_wind_CF0",
        "color":    "#ED7D31",
        "label":    "T1\n(past climate)",
    },
    "T3 (future)": {
        "scenario": "SomersetLevels_dec_t3",
        "folder":   "event_tp_ceh_t3_compass_CF0_GTSMv41opendap_CF0_no_wind_CF0",
        "color":    "#A9D18E",
        "label":    "T3\n(future climate)",
    },
}

FACTUAL_KEY     = "Factual"
COUNTERFACTUALS = ["T1 (past)", "T3 (future)"]

FLOOD_THRESHOLD = 0.05   # m
HMAX_VMAX       = 2.5    # m  (Somerset Levels — relatively shallow but wide)
DIFF_VMAX       = 0.5    # m

# Somerset ~51°N: 1° longitude ≈ 111000 × cos(51°) ≈ 69800 m
LAT_REF = 51.0
DX_PER_DEG = 111_000 * np.cos(np.radians(LAT_REF))   # ≈ 69 830 m/deg
DY_PER_DEG = 111_000                                   # m/deg


# ===== HELPER FUNCTIONS =====
def get_hmax_path(run_name: str) -> Path:
    r = RUNS[run_name]
    return (BASE_RUN_PATH / REGION / r["scenario"] / "sfincs"
            / r["folder"] / "plot_output" / "sfincs_output_hmax_AllTime.tif")


def load_hmax_tif(run_name: str) -> xr.DataArray:
    path = get_hmax_path(run_name)
    if not path.exists():
        raise FileNotFoundError(f"TIFF not found: {path}")
    da = rxr.open_rasterio(path)
    if "band" in da.dims:
        da = da.squeeze("band", drop=True)
    if da.rio.crs and str(da.rio.crs) != "EPSG:4326":
        da = da.rio.reproject("EPSG:4326")
    da = da.where(da > 0)   # mask no-data (0 or negative = not flooded)
    return da


def _cell_area(da: xr.DataArray):
    dx_deg = float(np.abs(da.x.diff("x").median()))
    dy_deg = float(np.abs(da.y.diff("y").median()))
    return dx_deg * DX_PER_DEG * dy_deg * DY_PER_DEG   # m²


def calculate_flood_volume(hmax: xr.DataArray) -> float:
    cell_m2 = _cell_area(hmax)
    return float(hmax.where(hmax > FLOOD_THRESHOLD).sum(skipna=True) * cell_m2)


def calculate_flood_extent(hmax: xr.DataArray) -> float:
    cell_m2 = _cell_area(hmax)
    return float((hmax > FLOOD_THRESHOLD).sum() * cell_m2) / 1e6   # km²


def calculate_mean_depth(hmax: xr.DataArray) -> float:
    return float(hmax.where(hmax > FLOOD_THRESHOLD).mean(skipna=True))


# ===== PLOTTING FUNCTIONS =====
def plot_hmax_threepanel(hmax_data: dict, output_path: Path):
    """3-panel map: Factual | T1 | T3 (all same colour scale)."""
    flood_levels = np.linspace(0, HMAX_VMAX, 21)

    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    for ax, run_name in zip(axes, [FACTUAL_KEY] + COUNTERFACTUALS):
        da = hmax_data[run_name].where(hmax_data[run_name] > FLOOD_THRESHOLD)
        im = da.plot(ax=ax, levels=flood_levels, cmap="Blues",
                     add_colorbar=False, x="x", y="y", alpha=0.85, zorder=2)
        ax.set_aspect("equal")
        ax.set_title(RUNS[run_name]["label"].replace("\n", " — "), fontsize=13, fontweight="bold")
        ax.set_xlabel("Longitude [°]", fontsize=11)
        ax.set_ylabel("Latitude [°]", fontsize=11)
        fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02).set_label("Max depth [m]", fontsize=10)

    plt.suptitle(
        f"Maximum Flood Depth — {EVENT_DISPLAY}\nFactual | T1 past climate | T3 future climate",
        fontsize=14, fontweight="bold", y=1.01,
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def plot_hmax_diff_panels(hmax_data: dict, output_path: Path):
    """2-panel attribution difference map: (Factual−T1) left, (T3−Factual) right."""
    hmax_fct = hmax_data[FACTUAL_KEY]
    hmax_t1  = hmax_data["T1 (past)"].interp_like(hmax_fct, method="nearest")
    hmax_t3  = hmax_data["T3 (future)"].interp_like(hmax_fct, method="nearest")

    diff_past   = hmax_fct - hmax_t1    # positive = more flooding in factual than past
    diff_future = hmax_t3 - hmax_fct    # positive = more flooding in future than factual

    all_vals = np.concatenate([
        diff_past.values[np.isfinite(diff_past.values)],
        diff_future.values[np.isfinite(diff_future.values)],
    ])
    vlim = float(np.quantile(np.abs(all_vals), 0.95)) if len(all_vals) else 1.0
    vlim = max(vlim, 0.05)
    levels = np.linspace(-vlim, vlim, 21)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    titles = [
        f"Factual − T1 (past climate)\n+ve = climate change increased depth",
        f"T3 (future) − Factual\n+ve = future climate increases depth",
    ]
    for ax, da, title in zip(axes, [diff_past, diff_future], titles):
        im = da.plot(ax=ax, levels=levels, cmap="RdBu_r",
                     add_colorbar=False, x="x", y="y", alpha=0.85, zorder=2)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xlabel("Longitude [°]", fontsize=11)
        ax.set_ylabel("Latitude [°]", fontsize=11)

    plt.suptitle(
        f"Climate-Attributable Flood Depth Change — {EVENT_DISPLAY}",
        fontsize=14, fontweight="bold", y=1.01,
    )
    plt.tight_layout(rect=[0, 0, 0.88, 1])
    cax = fig.add_axes([0.91, 0.15, 0.02, 0.68])
    fig.colorbar(im, cax=cax).set_label("Depth difference [m]", fontsize=10)
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def _bracket(ax, x_line, v_a, v_b, ymax, label):
    x_ext = 0.05
    lo, hi = min(v_a, v_b), max(v_a, v_b)
    if hi - lo < 1e-9:
        return
    ax.plot([x_line, x_line], [lo, hi], "k-", linewidth=1.6)
    for y in [lo, hi]:
        ax.plot([x_line - x_ext, x_line + x_ext], [y, y], "k-", linewidth=1.2)
    pct = (v_a - v_b) / v_b * 100 if v_b else 0.0
    sign = "+" if pct >= 0 else ""
    ax.text(x_line + 0.07, (lo + hi) / 2,
            f"{sign}{pct:.1f}%", va="center", ha="left", fontsize=8, style="italic")
    ax.text(x_line + 0.07, hi + ymax * 0.03, label,
            va="bottom", ha="left", fontsize=7, color="gray", style="italic")


def plot_attribution_metrics(metrics: dict, output_path: Path):
    """1×3 bar chart with all 3 scenarios and 2 attribution brackets."""
    specs = [
        ("volume",     "Flood Volume [Mm³]", 1e6),
        ("extent",     "Flood Extent [km²]", 1),
        ("mean_depth", "Mean Depth [m]",      1),
    ]
    colors = [RUNS[k]["color"] for k in [FACTUAL_KEY] + COUNTERFACTUALS]
    bar_labels = [RUNS[k]["label"] for k in [FACTUAL_KEY] + COUNTERFACTUALS]

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))

    for ax, (col, ylabel, scale) in zip(axes, specs):
        v_fct = metrics[FACTUAL_KEY][col] / scale
        v_t1  = metrics["T1 (past)"][col]  / scale
        v_t3  = metrics["T3 (future)"][col] / scale

        bars = ax.bar([0, 1, 2], [v_fct, v_t1, v_t3],
                      color=colors, edgecolor="black", linewidth=0.8, width=0.5)
        ax.set_ylabel(ylabel, fontsize=11, fontweight="bold")
        ax.set_title(ylabel.split(" [")[0], fontsize=12, fontweight="bold")
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(bar_labels, fontsize=9)
        ymax = max(v_fct, v_t1, v_t3) * 1.6 if max(v_fct, v_t1, v_t3) > 0 else 1.0
        ax.set_ylim(0, ymax)
        ax.set_xlim(-0.5, 3.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)

        for bar, val in zip(bars, [v_fct, v_t1, v_t3]):
            fmt = f"{val:.3f}" if col == "mean_depth" else f"{val:.2f}"
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + ymax * 0.01,
                    fmt, ha="center", va="bottom", fontsize=9, fontweight="bold")

        # Brackets: Factual vs T1 (climate change signal — past to present)
        _bracket(ax, 0.65, v_fct, v_t1, ymax, "Factual vs T1")
        # Bracket: T3 vs Factual (future change signal)
        _bracket(ax, 1.65, v_t3, v_fct, ymax, "T3 vs Factual")

    plt.suptitle(
        f"Flood Attribution — {EVENT_DISPLAY}\nFactual vs T1 (past) vs T3 (future) CEH-GEAR rainfall",
        fontsize=14, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def export_summary(metrics: dict, output_path: Path):
    w = 100
    print("\n" + "=" * w)
    print(f"  SFINCS FLOOD ATTRIBUTION SUMMARY — {EVENT_DISPLAY}")
    print("=" * w)
    print(f"{'Metric':<22} {'Factual':>14} {'T1 (past)':>14} {'T3 (future)':>14}"
          f" {'CC Fct/T1%':>11} {'CC T3/Fct%':>11}")
    print("-" * w)

    rows = []
    for col, label, scale, fmt in [
        ("volume",     "Flood Volume [Mm³]", 1e6, ".2f"),
        ("extent",     "Flood Extent [km²]", 1,   ".2f"),
        ("mean_depth", "Mean Depth [m]",      1,   ".4f"),
    ]:
        v_fct = metrics[FACTUAL_KEY][col] / scale
        v_t1  = metrics["T1 (past)"][col]  / scale
        v_t3  = metrics["T3 (future)"][col] / scale
        pct_t1  = (v_fct - v_t1) / v_t1  * 100 if v_t1  else 0.0   # +ve = fct wetter than past
        pct_t3  = (v_t3  - v_fct) / v_fct * 100 if v_fct else 0.0  # +ve = future wetter than fct
        f = f"{{:{fmt}}}"
        print(f"{label:<22} {f.format(v_fct):>14} {f.format(v_t1):>14} {f.format(v_t3):>14}"
              f" {pct_t1:>+10.1f}% {pct_t3:>+10.1f}%")
        rows.append({
            "Metric": label,
            "Factual": v_fct,
            "T1_past": v_t1,
            "T3_future": v_t3,
            "pct_change_Factual_vs_T1": round(pct_t1, 2),
            "pct_change_T3_vs_Factual": round(pct_t3, 2),
        })
    print("=" * w + "\n")
    pd.DataFrame(rows).to_csv(output_path, index=False)
    print(f"Saved: {output_path}")


# ===== MAIN =====
def main():
    print(f"\n{'='*65}")
    print(f"SFINCS Climate Attribution — {EVENT_DISPLAY}")
    print(f"{'='*65}\n")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading hmax rasters...")
    hmax_data = {}
    metrics   = {}
    for name in RUNS:
        print(f"  {name} ...")
        hmax = load_hmax_tif(name)
        hmax_data[name] = hmax
        metrics[name] = {
            "volume":     calculate_flood_volume(hmax),
            "extent":     calculate_flood_extent(hmax),
            "mean_depth": calculate_mean_depth(hmax),
        }
        m = metrics[name]
        print(f"    Volume: {m['volume']/1e6:.3f} Mm³  |  "
              f"Extent: {m['extent']:.2f} km²  |  Mean depth: {m['mean_depth']:.4f} m")

    print("\nGenerating figures...")

    plot_hmax_threepanel(
        hmax_data,
        OUTPUT_DIR / f"hmax_threepanel_{EVENT_NAME}.png",
    )
    plot_hmax_diff_panels(
        hmax_data,
        OUTPUT_DIR / f"hmax_attribution_diff_{EVENT_NAME}.png",
    )
    plot_attribution_metrics(
        metrics,
        OUTPUT_DIR / f"metrics_attribution_{EVENT_NAME}.png",
    )
    export_summary(
        metrics,
        OUTPUT_DIR / f"summary_attribution_{EVENT_NAME}.csv",
    )

    print(f"\n{'='*65}")
    print(f"Done. Outputs: {OUTPUT_DIR}")
    print(f"{'='*65}\n")


if __name__ == "__main__":
    main()

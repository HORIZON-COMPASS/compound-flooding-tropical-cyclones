"""
Climate Attribution Analysis: Factual vs Counterfactual Flooding (SFINCS).

Compares max flood depth outputs between:
- Factual (CF0):               Current climate precipitation (ERA5)
- Clausius-Clapeyron (CF-8):   Pre-industrial proxy via CC scaling (-8%)
- ClimateDT (CF-10):           Pre-industrial proxy via ClimateDT model (-10%)

Produces per counterfactual:
1. 3-panel hmax map (Factual | Counterfactual | Difference)

Produces combined:
2. Bar chart for flood metrics (volume, extent, mean depth) — 3 bars, 2 brackets
3. Summary CSV

Author: Generated for Durban 2022 climate attribution analysis
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
EVENT_NAME    = "Durban_April2022"
EVENT_DISPLAY = "Durban 2022"
REGION        = "durban"
RUNNAME       = "Durban2022"

BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs")
OUTPUT_DIR    = Path("/p/11210471-001-compass/04_Results/climate_attribution")

RUNS = {
    "Factual": {
        "folder":   "event_precip_era5_hourly_CF0_no_wind",
        "cf_value": 0,
        "color":    "#FF6B6B",
    },
    "Clausius-Clapeyron": {
        "folder":   "event_precip_era5_hourly_CF-8_no_wind",
        "cf_value": -8,
        "color":    "#4ECDC4",
    },
    "ClimateDT": {
        "folder":   "event_precip_era5_hourly_CF-10_no_wind",
        "cf_value": -10,
        "color":    "#FFE66D",
    },
}

FACTUAL_KEY      = "Factual"
COUNTERFACTUALS  = ["Clausius-Clapeyron", "ClimateDT"]

FLOOD_THRESHOLD = 0.05   # m
HMAX_VMAX       = 2.0    # m
DIFF_VMAX       = 0.5    # m


# ===== HELPER FUNCTIONS =====
def get_run_path(run_name: str) -> Path:
    return BASE_RUN_PATH / REGION / RUNNAME / "sfincs" / RUNS[run_name]["folder"]


def load_hmax_tif(run_name: str) -> xr.DataArray:
    path = get_run_path(run_name) / "plot_output" / "sfincs_output_hmax_AllTime.tif"
    if not path.exists():
        raise FileNotFoundError(f"TIFF not found: {path}")
    da = rxr.open_rasterio(path)
    if "band" in da.dims:
        da = da.squeeze("band", drop=True)
    if da.rio.crs != "EPSG:4326":
        da = da.rio.reproject("EPSG:4326")
    return da


def calculate_flood_volume(hmax: xr.DataArray) -> float:
    dx_m = float(np.abs(hmax.x.diff("x").median())) * 96_000
    dy_m = float(np.abs(hmax.y.diff("y").median())) * 111_000
    return float(hmax.where(hmax > FLOOD_THRESHOLD).sum(skipna=True) * dx_m * dy_m)


def calculate_flood_extent(hmax: xr.DataArray) -> float:
    dx_m = float(np.abs(hmax.x.diff("x").median())) * 96_000
    dy_m = float(np.abs(hmax.y.diff("y").median())) * 111_000
    return float((hmax > FLOOD_THRESHOLD).sum() * dx_m * dy_m) / 1e6


def calculate_mean_depth(hmax: xr.DataArray) -> float:
    return float(hmax.where(hmax > FLOOD_THRESHOLD).mean(skipna=True))


# ===== PLOTTING FUNCTIONS =====
def plot_hmax_comparison(
    hmax_fct: xr.DataArray,
    hmax_cft: xr.DataArray,
    cft_label: str,
    output_path: Path,
):
    """3-panel map: Factual | Counterfactual | Difference."""
    flood_levels = np.linspace(0, HMAX_VMAX, 21)
    diff_levels  = np.linspace(-DIFF_VMAX, DIFF_VMAX, 21)

    hmax_fct_aligned = hmax_fct.interp_like(hmax_cft, method="nearest")
    diff = hmax_fct_aligned - hmax_cft

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    panels = [
        (hmax_fct.where(hmax_fct > FLOOD_THRESHOLD), flood_levels, "viridis",
         "Factual", "Max Depth [m]"),
        (hmax_cft.where(hmax_cft > FLOOD_THRESHOLD), flood_levels, "viridis",
         cft_label, "Max Depth [m]"),
        (diff, diff_levels, "RdBu_r",
         f"Climate-Attributable\n(Factual − {cft_label})", "Difference [m]"),
    ]
    for ax, (da, levels, cmap, title, cbar_lbl) in zip(axes, panels):
        im = da.plot(ax=ax, levels=levels, cmap=cmap,
                     add_colorbar=False, x="x", y="y", alpha=0.75, zorder=2)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.set_xlabel("Longitude [°]", fontsize=11)
        ax.set_ylabel("Latitude [°]", fontsize=11)
        try:
            ctx.add_basemap(ax=ax, source=ctx.providers.OpenStreetMap.Mapnik,
                            crs="EPSG:4326", attribution=False, zorder=1, zoom=11)
        except Exception:
            pass
        fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02).set_label(cbar_lbl, fontsize=10)

    mean_diff = float(diff.mean(skipna=True))
    max_diff  = float(diff.max(skipna=True))
    plt.suptitle(
        f"Flood Attribution — {EVENT_DISPLAY}  |  Factual vs {cft_label}\n"
        f"Mean additional depth: {mean_diff:.3f} m  |  Max additional depth: {max_diff:.3f} m",
        fontsize=13, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def _bracket(ax, x_line, v_fct, v_cft, attr_pct, ymax, unit, label):
    """Draw an attribution bracket between two bar heights."""
    x_ext = 0.05
    lo, hi = min(v_fct, v_cft), max(v_fct, v_cft)
    ax.plot([x_line, x_line], [lo, hi], "k-", linewidth=1.6)
    for y in [lo, hi]:
        ax.plot([x_line - x_ext, x_line + x_ext], [y, y], "k-", linewidth=1.2)
    for v in [v_fct, v_cft]:
        ax.plot([ax.get_xlim()[0], x_line], [v, v],
                linestyle="--", color="gray", alpha=0.4, linewidth=0.7)
    sign = "+" if (v_fct - v_cft) >= 0 else ""
    diff = v_fct - v_cft
    txt = f"{sign}{diff:.2f}{unit}\n({attr_pct:.0f}%)" if unit else \
          f"{sign}{diff:.2f}\n({attr_pct:.0f}%)"
    ax.text(x_line + 0.07, (lo + hi) / 2, txt,
            va="center", ha="left", fontsize=8, style="italic")
    ax.text(x_line + 0.07, hi + ymax * 0.04, label,
            va="bottom", ha="left", fontsize=7, style="italic", color="gray")


def plot_attribution_metrics(metrics: dict, output_path: Path):
    """1×3 bar chart — 3 bars per metric, 2 attribution brackets."""
    fct = metrics[FACTUAL_KEY]
    cc  = metrics[COUNTERFACTUALS[0]]
    cdt = metrics[COUNTERFACTUALS[1]]

    specs = [
        ("volume", "Flood Volume [Mm³]",  1e6,  ""),
        ("extent", "Flood Extent [km²]",  1,    ""),
        ("mean_depth", "Mean Depth [m]",  1,    ""),
    ]
    colors = [RUNS[FACTUAL_KEY]["color"],
              RUNS[COUNTERFACTUALS[0]]["color"],
              RUNS[COUNTERFACTUALS[1]]["color"]]
    bar_labels = ["Factual", "Clausius-\nClapeyron", "ClimateDT"]

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))

    for ax, (col, ylabel, scale, unit) in zip(axes, specs):
        v_fct = fct[col] / scale
        v_cc  = cc[col]  / scale
        v_cdt = cdt[col] / scale
        pct_cc  = (v_fct - v_cc)  / v_fct * 100 if v_fct else 0.0
        pct_cdt = (v_fct - v_cdt) / v_fct * 100 if v_fct else 0.0

        bars = ax.bar([0, 1, 2], [v_fct, v_cc, v_cdt],
                      color=colors, edgecolor="black", linewidth=1, width=0.5)
        ax.set_ylabel(ylabel, fontsize=11, fontweight="bold")
        ax.set_title(ylabel.split(" [")[0], fontsize=12, fontweight="bold")
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(bar_labels, fontsize=10)
        ymax = max(v_fct, v_cc, v_cdt) * 1.55 if max(v_fct, v_cc, v_cdt) > 0 else 1.0
        ax.set_ylim(0, ymax)
        ax.set_xlim(-0.5, 3.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)

        for bar, val in zip(bars, [v_fct, v_cc, v_cdt]):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + ymax * 0.012,
                        f"{val:.2f}", ha="center", va="bottom",
                        fontsize=9, fontweight="bold")

        if v_fct > 0 and v_cc > 0:
            _bracket(ax, 1.45, v_fct, v_cc, pct_cc, ymax, unit,
                     "Clausius-Clapeyron\nattribution:")
        if v_fct > 0 and v_cdt > 0:
            _bracket(ax, 2.55, v_fct, v_cdt, pct_cdt, ymax, unit,
                     "ClimateDT\nattribution:")

    plt.suptitle(
        f"Flood Attribution — {EVENT_DISPLAY}\n"
        f"Factual vs Clausius-Clapeyron & ClimateDT counterfactuals",
        fontsize=14, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def export_summary(metrics: dict, output_path: Path):
    fct = metrics[FACTUAL_KEY]
    cc  = metrics[COUNTERFACTUALS[0]]
    cdt = metrics[COUNTERFACTUALS[1]]

    w = 95
    print("\n" + "=" * w)
    print(f"  SFINCS FLOOD ATTRIBUTION SUMMARY — {EVENT_DISPLAY}")
    print("=" * w)
    print(f"{'Metric':<22} {'Factual':>13} {'CC (CF-8)':>13} {'CDT (CF-10)':>13}"
          f" {'Attr CC%':>9} {'Attr CDT%':>10}")
    print("-" * w)

    rows = []
    for col, label, scale, fmt in [
        ("volume",     "Flood Volume [Mm³]", 1e6, ".2f"),
        ("extent",     "Flood Extent [km²]", 1,   ".2f"),
        ("mean_depth", "Mean Depth [m]",      1,   ".3f"),
    ]:
        v_fct = fct[col] / scale
        v_cc  = cc[col]  / scale
        v_cdt = cdt[col] / scale
        pct_cc  = (v_fct - v_cc)  / v_fct * 100 if v_fct else 0.0
        pct_cdt = (v_fct - v_cdt) / v_fct * 100 if v_fct else 0.0
        f = f"{{:{fmt}}}"
        print(f"{label:<22} {f.format(v_fct):>13} {f.format(v_cc):>13} {f.format(v_cdt):>13}"
              f" {pct_cc:>8.1f}% {pct_cdt:>9.1f}%")
        rows.append({
            "Metric": label,
            "Factual": fct[col] / scale,
            "Clausius_Clapeyron_CF-8": v_cc,
            "ClimateDT_CF-10": v_cdt,
            "Attribution_CC_pct": round(pct_cc, 2),
            "Attribution_CDT_pct": round(pct_cdt, 2),
        })

    print("=" * w + "\n")
    pd.DataFrame(rows).to_csv(output_path, index=False)
    print(f"Saved: {output_path}")


def plot_hmax_attribution_comparison(hmax_data: dict, output_path: Path):
    """
    2-panel map: climate-attributable hmax difference for CC (left) and CDT (right).

    Both panels share the same diverging colour scale so the two attribution
    methods are directly comparable at a glance.
    """
    hmax_fct = hmax_data[FACTUAL_KEY]
    cc_label  = COUNTERFACTUALS[0].split(" (")[0] if " (" in COUNTERFACTUALS[0] else COUNTERFACTUALS[0]
    cdt_label = COUNTERFACTUALS[1].split(" (")[0] if " (" in COUNTERFACTUALS[1] else COUNTERFACTUALS[1]

    diffs = {}
    for cft_key in COUNTERFACTUALS:
        hmax_cft = hmax_data[cft_key]
        hmax_fct_aligned = hmax_fct.interp_like(hmax_cft, method="nearest")
        diffs[cft_key] = hmax_fct_aligned - hmax_cft

    # Shared symmetric colour limits across both panels
    all_vals = np.concatenate([
        diffs[k].values[~np.isnan(diffs[k].values)] for k in COUNTERFACTUALS
    ])
    diff_abs = float(np.quantile(np.abs(all_vals), 0.95)) if len(all_vals) else 1.0
    if diff_abs == 0:
        diff_abs = 1.0
    diff_levels = np.linspace(-diff_abs, diff_abs, 21)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    panels = [
        (diffs[COUNTERFACTUALS[0]], f"Factual − {cc_label}"),
        (diffs[COUNTERFACTUALS[1]], f"Factual − {cdt_label}"),
    ]
    ims = []
    for ax, (da, title) in zip(axes, panels):
        im = da.plot(ax=ax, levels=diff_levels, cmap="RdBu_r",
                     add_colorbar=False, x="x", y="y", alpha=0.75, zorder=2)
        ims.append(im)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.set_xlabel("Longitude [°]", fontsize=11)
        ax.set_ylabel("Latitude [°]", fontsize=11)
        try:
            ctx.add_basemap(ax=ax, source=ctx.providers.OpenStreetMap.Mapnik,
                            crs="EPSG:4326", attribution=False, zorder=1, zoom=11)
        except Exception:
            pass

    plt.suptitle(
        f"Climate-Attributable Flood Depth — {EVENT_DISPLAY}\n"
        f"Clausius-Clapeyron vs ClimateDT attribution",
        fontsize=14, fontweight="bold", y=1.02,
    )
    plt.tight_layout(rect=[0, 0, 0.88, 1])
    cax = fig.add_axes([0.91, 0.15, 0.02, 0.68])
    fig.colorbar(ims[0], cax=cax).set_label("Depth Difference [m]", fontsize=10)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# ===== MAIN =====
def main():
    print(f"\n{'='*60}")
    print(f"SFINCS Climate Attribution — {EVENT_DISPLAY}")
    print(f"{'='*60}\n")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading hmax rasters...")
    hmax_data = {}
    metrics   = {}
    for name in RUNS:
        try:
            print(f"  Loading {name}...")
            hmax = load_hmax_tif(name)
            hmax_data[name] = hmax
            metrics[name] = {
                "volume":     calculate_flood_volume(hmax),
                "extent":     calculate_flood_extent(hmax),
                "mean_depth": calculate_mean_depth(hmax),
            }
            m = metrics[name]
            print(f"    Volume: {m['volume']/1e6:.2f} Mm³  |  "
                  f"Extent: {m['extent']:.2f} km²  |  "
                  f"Mean depth: {m['mean_depth']:.3f} m")
        except Exception as e:
            print(f"  ERROR loading {name}: {e}")
            return

    print("\nGenerating plots...")

    for cft in COUNTERFACTUALS:
        slug = cft.lower().replace(" ", "_").replace("-", "")
        print(f"  hmax map: Factual vs {cft}...")
        plot_hmax_comparison(
            hmax_data[FACTUAL_KEY], hmax_data[cft], cft,
            OUTPUT_DIR / f"hmax_attribution_{slug}_{EVENT_NAME.lower()}.png",
        )

    print("  hmax attribution comparison (CC vs CDT side-by-side)...")
    plot_hmax_attribution_comparison(
        hmax_data,
        OUTPUT_DIR / f"hmax_attribution_comparison_{EVENT_NAME.lower()}.png",
    )

    print("  Attribution metrics bar chart...")
    plot_attribution_metrics(
        metrics,
        OUTPUT_DIR / f"metrics_attribution_{EVENT_NAME.lower()}.png",
    )

    print("  Exporting summary CSV...")
    export_summary(
        metrics,
        OUTPUT_DIR / f"summary_attribution_{EVENT_NAME.lower()}.csv",
    )

    print(f"\n{'='*60}")
    print(f"Done. Outputs saved to: {OUTPUT_DIR}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()

"""
Climate Attribution Analysis: Factual vs Counterfactual Flood Impact (Delft-FIAT).

Compares building-level flood impact outputs between:
- Factual (CF0):                       Current climate precipitation (ERA5)
- Clausius-Clapeyron (CF-8):           Pre-industrial proxy via CC scaling (-8%)
- ClimateDT (CF-10):                   Pre-industrial proxy via ClimateDT model (-10%)

Produces:
1. Rasterized damage maps per counterfactual (Factual | CFX | Difference)
2. Rasterized population-exposed maps per counterfactual
3. 2-panel attribution comparison map — damage difference CC vs CDT side-by-side
4. 2-panel attribution comparison map — population difference CC vs CDT side-by-side
5. Attribution bar chart — 3 bars per metric, 2 attribution brackets
6. Printed summary table and CSV export

Population metric requires running fiat_add_pop_exposed_metric.py first to generate
spatial_with_pop_and_flood.fgb.  Falls back to 0 if that file is not yet available.

by @dumontgoulart
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.colors import BoundaryNorm
import contextily as ctx
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")


# ===== CONFIGURATION =====
EVENT_NAME = "Durban_April2022"   # used in output file names — keep as-is
EVENT_DISPLAY = "Durban 2022"     # used in plot titles
REGION = "durban"
RUNNAME = "Durban2022"

BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs")
OUTPUT_DIR = Path("/p/11210471-001-compass/04_Results/fiat_climate_attribution")

RUNS = {
    "Factual (CF0)": {
        "folder": "event_precip_era5_hourly_CF0_no_wind",
        "cf_value": 0,
        "color": "#FF6B6B",  # coral
    },
    "Clausius-Clapeyron (CF-8)": {
        "folder": "event_precip_era5_hourly_CF-8_no_wind",
        "cf_value": -8,
        "color": "#4ECDC4",  # teal
    },
    "ClimateDT (CF-10)": {
        "folder": "event_precip_era5_hourly_CF-10_no_wind",
        "cf_value": -10,
        "color": "#FFE66D",  # gold
    },
}

FACTUAL_KEY = "Factual (CF0)"
COUNTERFACTUALS = ["Clausius-Clapeyron (CF-8)", "ClimateDT (CF-10)"]

GRID_RESOLUTION = 0.01  # degrees (~1 km at Durban latitude)


# ===== HELPER FUNCTIONS =====
def get_fiat_path(run_name: str) -> Path:
    """Return the FIAT output directory for the given scenario name."""
    return BASE_RUN_PATH / REGION / RUNNAME / "fiat" / RUNS[run_name]["folder"]


def load_spatial_fgb(run_name: str) -> gpd.GeoDataFrame:
    """Load spatial output and reproject to EPSG:4326.

    Prefers spatial_with_pop_and_flood.fgb (has 'population' column) when
    available; falls back to spatial.fgb otherwise.
    """
    base = get_fiat_path(run_name) / "output"
    path = base / "spatial_with_pop_and_flood.fgb"
    if not path.exists():
        path = base / "spatial.fgb"
    if not path.exists():
        raise FileNotFoundError(f"FIAT output not found in {base}")
    gdf = gpd.read_file(path)
    if str(gdf.crs) != "EPSG:4326":
        gdf = gdf.to_crs("EPSG:4326")
    return gdf


def compute_impact_metrics(gdf: gpd.GeoDataFrame) -> dict:
    """Compute aggregate impact metrics from a FIAT spatial output GeoDataFrame."""
    affected = gdf[gdf["total_damage"] > 0]
    pop_exposed = float(affected["population"].sum()) if "population" in affected.columns else 0.0
    return {
        "total_damage": float(affected["total_damage"].sum()),
        "n_affected": int(len(affected)),
        "mean_damage_per_building": float(affected["total_damage"].mean()) if len(affected) > 0 else 0.0,
        "pop_exposed": pop_exposed,
    }


def rasterize_damage(
    gdf: gpd.GeoDataFrame,
    value_col: str = "total_damage",
    resolution: float = GRID_RESOLUTION,
    bounds: tuple = None,
) -> xr.DataArray:
    """
    Aggregate building-level values onto a regular lat/lon grid (sum per cell).

    Parameters
    ----------
    gdf : GeoDataFrame in EPSG:4326
    value_col : column to rasterize
    resolution : grid cell size in degrees
    bounds : (minx, miny, maxx, maxy); uses data extent if None

    Returns
    -------
    xr.DataArray with NaN where no buildings are present
    """
    gdf = gdf[gdf[value_col] > 0].copy()
    if len(gdf) == 0:
        return None

    centroids = gdf.geometry.centroid
    xs = centroids.x.values
    ys = centroids.y.values
    vals = gdf[value_col].values

    if bounds is None:
        bounds = (xs.min(), ys.min(), xs.max(), ys.max())
    minx, miny, maxx, maxy = bounds
    minx -= resolution / 2
    maxx += resolution / 2
    miny -= resolution / 2
    maxy += resolution / 2

    x_edges = np.arange(minx, maxx + resolution, resolution)
    y_edges = np.arange(miny, maxy + resolution, resolution)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2

    x_bins = np.clip(np.digitize(xs, x_edges) - 1, 0, len(x_centers) - 1)
    y_bins = np.clip(np.digitize(ys, y_edges) - 1, 0, len(y_centers) - 1)

    grid = np.zeros((len(y_centers), len(x_centers)))
    count = np.zeros_like(grid)
    for i in range(len(vals)):
        if not np.isnan(vals[i]):
            grid[y_bins[i], x_bins[i]] += vals[i]
            count[y_bins[i], x_bins[i]] += 1
    grid = np.where(count > 0, grid, np.nan)

    return xr.DataArray(
        data=grid,
        dims=["y", "x"],
        coords={"y": y_centers, "x": x_centers},
    )


def _common_bounds(gdfs: dict) -> tuple:
    """Compute lon/lat bounding box covering all scenarios' affected buildings."""
    all_minx, all_miny, all_maxx, all_maxy = [], [], [], []
    for gdf in gdfs.values():
        affected = gdf[gdf["total_damage"] > 0]
        if len(affected) == 0:
            continue
        c = affected.geometry.centroid
        all_minx.append(c.x.min())
        all_miny.append(c.y.min())
        all_maxx.append(c.x.max())
        all_maxy.append(c.y.max())
    return (min(all_minx), min(all_miny), max(all_maxx), max(all_maxy))


# ===== PLOTTING FUNCTIONS =====
def _three_panel_spatial(
    da_fct: xr.DataArray,
    da_cft: xr.DataArray,
    fct_label: str,
    cft_label: str,
    value_cbar: str,
    diff_cbar: str,
    scenario_cmap,
    suptitle: str,
    output_path: Path,
):
    """
    Generic 3-panel spatial comparison: Factual | Counterfactual | Difference.

    Used by both the damage and population map functions.
    """
    da_fct_f = da_fct.fillna(0)
    da_cft_f = da_cft.fillna(0)
    da_diff = (da_fct_f - da_cft_f).where((da_fct_f > 0) | (da_cft_f > 0))

    all_vals = np.concatenate([
        da_fct.values[~np.isnan(da_fct.values)],
        da_cft.values[~np.isnan(da_cft.values)],
    ])
    vmax = float(np.quantile(all_vals, 0.95)) if len(all_vals) > 0 else 1.0
    val_norm = BoundaryNorm(np.linspace(0, vmax, 11), ncolors=256, clip=True)

    diff_vals = da_diff.values[~np.isnan(da_diff.values)]
    diff_abs = float(np.quantile(np.abs(diff_vals), 0.95)) if len(diff_vals) > 0 else 1.0
    if diff_abs == 0:
        diff_abs = 1.0
    diff_norm = BoundaryNorm(np.linspace(-diff_abs, diff_abs, 11), ncolors=256, clip=True)
    diff_cmap = plt.get_cmap("RdBu_r")

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    panels = [
        (da_fct, val_norm,  scenario_cmap, fct_label,                                  value_cbar),
        (da_cft, val_norm,  scenario_cmap, cft_label,                                  value_cbar),
        (da_diff, diff_norm, diff_cmap,    f"Climate-Attributable\n(Factual − {cft_label})", diff_cbar),
    ]
    for ax, (da, norm, cmap, title, cbar_label) in zip(axes, panels):
        im = da.plot(ax=ax, cmap=cmap, norm=norm, add_colorbar=False,
                     x="x", y="y", alpha=0.85, zorder=2)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xlabel("Longitude [°]", fontsize=10)
        ax.set_ylabel("Latitude [°]", fontsize=10)
        try:
            ctx.add_basemap(ax=ax, source=ctx.providers.OpenStreetMap.Mapnik,
                            crs="EPSG:4326", attribution=False, zorder=1, zoom=12)
        except Exception:
            pass
        plt.colorbar(im, ax=ax, shrink=0.8, pad=0.1).set_label(cbar_label, fontsize=9)

    plt.suptitle(suptitle, fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def plot_damage_spatial_comparison(
    fct_gdf: gpd.GeoDataFrame,
    cft_gdf: gpd.GeoDataFrame,
    cft_label: str,
    bounds: tuple,
    output_path: Path,
):
    """3-panel rasterized damage map for one factual–counterfactual pair."""
    da_fct = rasterize_damage(fct_gdf, bounds=bounds)
    da_cft = rasterize_damage(cft_gdf, bounds=bounds)
    if da_fct is None or da_cft is None:
        print(f"  Insufficient damage data for {cft_label}, skipping.")
        return
    cft_display = cft_label.split(" (")[0]  # strip "(CF-X)" for display
    _three_panel_spatial(
        da_fct, da_cft,
        fct_label="Factual",
        cft_label=cft_display,
        value_cbar=f"Total Damage per {GRID_RESOLUTION}° cell [USD]",
        diff_cbar=f"Damage Difference per {GRID_RESOLUTION}° cell [USD]",
        scenario_cmap=plt.get_cmap("Reds"),
        suptitle=f"Economic Flood Damage (Rasterized {GRID_RESOLUTION}°) — {EVENT_DISPLAY}\nFactual vs {cft_display}",
        output_path=output_path,
    )


def plot_pop_spatial_comparison(
    fct_gdf: gpd.GeoDataFrame,
    cft_gdf: gpd.GeoDataFrame,
    cft_label: str,
    bounds: tuple,
    output_path: Path,
):
    """3-panel rasterized population-exposed map for one factual–counterfactual pair."""
    if "population" not in fct_gdf.columns or "population" not in cft_gdf.columns:
        print(f"  No population column for {cft_label} — run fiat_add_pop_exposed_metric.py first.")
        return
    da_fct = rasterize_damage(fct_gdf, value_col="population", bounds=bounds)
    da_cft = rasterize_damage(cft_gdf, value_col="population", bounds=bounds)
    if da_fct is None or da_cft is None:
        print(f"  Insufficient population data for {cft_label}, skipping.")
        return
    cft_display = cft_label.split(" (")[0]  # strip "(CF-X)" for display
    _three_panel_spatial(
        da_fct, da_cft,
        fct_label="Factual",
        cft_label=cft_display,
        value_cbar=f"Population Exposed per {GRID_RESOLUTION}° cell",
        diff_cbar=f"Population Difference per {GRID_RESOLUTION}° cell",
        scenario_cmap=plt.get_cmap("YlOrRd"),
        suptitle=f"Population Exposed to Flooding (Rasterized {GRID_RESOLUTION}°) — {EVENT_DISPLAY}\nFactual vs {cft_display}",
        output_path=output_path,
    )


def _two_panel_attribution_diff(
    da_cc_diff: xr.DataArray,
    da_cdt_diff: xr.DataArray,
    diff_cbar: str,
    suptitle: str,
    output_path: Path,
):
    """
    2-panel map showing the climate-attributable difference for CC (left) and CDT (right).

    Both panels share the same diverging colour scale so the two methods are
    directly comparable at a glance.
    """
    # Shared symmetric colour limits across both panels
    all_diffs = []
    for da in [da_cc_diff, da_cdt_diff]:
        v = da.values[~np.isnan(da.values)]
        if len(v):
            all_diffs.append(v)
    if all_diffs:
        combined = np.concatenate(all_diffs)
        diff_abs = float(np.quantile(np.abs(combined), 0.95))
        if diff_abs == 0:
            diff_abs = 1.0
    else:
        diff_abs = 1.0

    diff_norm = BoundaryNorm(np.linspace(-diff_abs, diff_abs, 11), ncolors=256, clip=True)
    diff_cmap = plt.get_cmap("RdBu_r")

    cc_label  = COUNTERFACTUALS[0].split(" (")[0]
    cdt_label = COUNTERFACTUALS[1].split(" (")[0]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    panels = [
        (da_cc_diff,  f"Factual − {cc_label}"),
        (da_cdt_diff, f"Factual − {cdt_label}"),
    ]
    ims = []
    for ax, (da, title) in zip(axes, panels):
        im = da.plot(ax=ax, cmap=diff_cmap, norm=diff_norm, add_colorbar=False,
                     x="x", y="y", alpha=0.85, zorder=2)
        ims.append(im)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.set_xlabel("Longitude [°]", fontsize=11)
        ax.set_ylabel("Latitude [°]", fontsize=11)
        try:
            ctx.add_basemap(ax=ax, source=ctx.providers.OpenStreetMap.Mapnik,
                            crs="EPSG:4326", attribution=False, zorder=1, zoom=12)
        except Exception:
            pass

    plt.suptitle(suptitle, fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout(rect=[0, 0, 0.88, 1])
    cax = fig.add_axes([0.91, 0.15, 0.02, 0.68])
    cbar = fig.colorbar(ims[0], cax=cax)
    cbar.set_label(diff_cbar, fontsize=10)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def plot_damage_attribution_comparison(
    gdfs: dict,
    bounds: tuple,
    output_path: Path,
):
    """2-panel figure: climate-attributable damage difference for CC vs CDT."""
    da_fct = rasterize_damage(gdfs[FACTUAL_KEY], bounds=bounds)
    if da_fct is None:
        print("  No factual damage data for attribution comparison, skipping.")
        return

    da_fct_f = da_fct.fillna(0)
    diffs = {}
    for cft_key in COUNTERFACTUALS:
        da_cft = rasterize_damage(gdfs[cft_key], bounds=bounds)
        if da_cft is None:
            print(f"  No data for {cft_key}, skipping attribution comparison.")
            return
        da_cft_f = da_cft.fillna(0)
        diffs[cft_key] = (da_fct_f - da_cft_f).where((da_fct_f > 0) | (da_cft_f > 0))

    _two_panel_attribution_diff(
        da_cc_diff=diffs[COUNTERFACTUALS[0]],
        da_cdt_diff=diffs[COUNTERFACTUALS[1]],
        diff_cbar=f"Damage Difference per {GRID_RESOLUTION}° cell [USD]",
        suptitle=(
            f"Climate-Attributable Economic Damage — {EVENT_DISPLAY}\n"
            f"Clausius-Clapeyron vs ClimateDT attribution"
        ),
        output_path=output_path,
    )


def plot_pop_attribution_comparison(
    gdfs: dict,
    bounds: tuple,
    output_path: Path,
):
    """2-panel figure: climate-attributable population-exposed difference for CC vs CDT."""
    for key in [FACTUAL_KEY] + COUNTERFACTUALS:
        if "population" not in gdfs[key].columns:
            print(f"  No population column for {key} — run fiat_add_pop_exposed_metric.py first.")
            return

    da_fct = rasterize_damage(gdfs[FACTUAL_KEY], value_col="population", bounds=bounds)
    if da_fct is None:
        print("  No factual population data for attribution comparison, skipping.")
        return

    da_fct_f = da_fct.fillna(0)
    diffs = {}
    for cft_key in COUNTERFACTUALS:
        da_cft = rasterize_damage(gdfs[cft_key], value_col="population", bounds=bounds)
        if da_cft is None:
            print(f"  No population data for {cft_key}, skipping attribution comparison.")
            return
        da_cft_f = da_cft.fillna(0)
        diffs[cft_key] = (da_fct_f - da_cft_f).where((da_fct_f > 0) | (da_cft_f > 0))

    _two_panel_attribution_diff(
        da_cc_diff=diffs[COUNTERFACTUALS[0]],
        da_cdt_diff=diffs[COUNTERFACTUALS[1]],
        diff_cbar=f"Population Difference per {GRID_RESOLUTION}° cell",
        suptitle=(
            f"Climate-Attributable Population Exposed — {EVENT_DISPLAY}\n"
            f"Clausius-Clapeyron vs ClimateDT attribution"
        ),
        output_path=output_path,
    )


def _attribution_bracket(ax, x_line, v_fct, v_cft, v_attr, attr_pct, ymax, unit, label):
    """Draw a single attribution bracket at x_line between factual and one counterfactual."""
    x_ext = 0.05
    lo, hi = min(v_fct, v_cft), max(v_fct, v_cft)
    ax.plot([x_line, x_line], [lo, hi], "k-", linewidth=1.6)
    for y_tick in [lo, hi]:
        ax.plot([x_line - x_ext, x_line + x_ext], [y_tick, y_tick], "k-", linewidth=1.2)
    # Dashed horizontal guides
    for val in [v_fct, v_cft]:
        ax.plot([ax.get_xlim()[0], x_line], [val, val],
                linestyle="--", color="gray", alpha=0.4, linewidth=0.7)
    mid_y = (lo + hi) / 2
    sign = "+" if v_attr >= 0 else ""
    if unit:
        attr_txt = f"{sign}{v_attr:.1f}{unit}\n({attr_pct:.0f}%)"
    else:
        attr_txt = f"{sign}{v_attr:,.0f}\n({attr_pct:.0f}%)"
    ax.text(x_line + 0.07, mid_y, attr_txt,
            va="center", ha="left", fontsize=8, style="italic", color="black")
    ax.text(x_line + 0.07, hi + ymax * 0.04, label,
            va="bottom", ha="left", fontsize=7, style="italic", color="gray")


def plot_attribution_metrics(metrics: dict, output_path: Path):
    """
    2×2 bar chart with 3 bars per subplot (Factual | CC | ClimateDT).

    Two attribution brackets are drawn per subplot — one for each counterfactual
    vs the factual — labelled with their approach name.
    """
    fct = metrics[FACTUAL_KEY]
    cc  = metrics[COUNTERFACTUALS[0]]
    cdt = metrics[COUNTERFACTUALS[1]]

    metric_specs = [
        ("total_damage",            "Total Damage [USD]",              1e6,  "M"),
        ("n_affected",              "Affected Buildings [#]",          1,    ""),
        ("mean_damage_per_building","Mean Damage per Building [USD]",  1e3,  "K"),
        ("pop_exposed",             "Population Exposed [#]",          1,    ""),
    ]

    bar_colors  = [RUNS[FACTUAL_KEY]["color"],
                   RUNS[COUNTERFACTUALS[0]]["color"],
                   RUNS[COUNTERFACTUALS[1]]["color"]]
    bar_labels  = [
        "Factual",
        "Clausius-\nClapeyron",
        "ClimateDT",
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for ax, (col, ylabel, scale, unit) in zip(axes, metric_specs):
        v_fct = fct[col] / scale
        v_cc  = cc[col]  / scale
        v_cdt = cdt[col] / scale

        attr_cc  = v_fct - v_cc
        attr_cdt = v_fct - v_cdt
        pct_cc   = (attr_cc  / v_fct * 100) if v_fct != 0 else 0.0
        pct_cdt  = (attr_cdt / v_fct * 100) if v_fct != 0 else 0.0

        bars = ax.bar([0, 1, 2], [v_fct, v_cc, v_cdt],
                      color=bar_colors, edgecolor="black", linewidth=1, width=0.5)

        ax.set_ylabel(ylabel, fontsize=10, fontweight="bold")
        ax.set_title(ylabel.split(" [")[0], fontsize=11, fontweight="bold")
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(bar_labels, fontsize=9)
        ymax = max(v_fct, v_cc, v_cdt) * 1.55 if max(v_fct, v_cc, v_cdt) > 0 else 1.0
        ax.set_ylim(0, ymax)
        ax.set_xlim(-0.5, 3.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)

        # Value labels on bars
        for bar, val in zip(bars, [v_fct, v_cc, v_cdt]):
            if val > 0:
                lbl = f"{val:.1f}{unit}" if unit else f"{val:,.0f}" if val >= 10 else f"{val:.2f}"
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + ymax * 0.012,
                        lbl, ha="center", va="bottom", fontsize=9, fontweight="bold")

        # Attribution bracket: Factual vs CC (between bars 1 and 2)
        if v_fct > 0 and v_cc > 0:
            _attribution_bracket(ax, x_line=1.45, v_fct=v_fct, v_cft=v_cc,
                                 v_attr=attr_cc, attr_pct=pct_cc, ymax=ymax,
                                 unit=unit, label="Clausius-Clapeyron\nattribution:")

        # Attribution bracket: Factual vs ClimateDT (to the right of bar 2)
        if v_fct > 0 and v_cdt > 0:
            _attribution_bracket(ax, x_line=2.55, v_fct=v_fct, v_cft=v_cdt,
                                 v_attr=attr_cdt, attr_pct=pct_cdt, ymax=ymax,
                                 unit=unit, label="ClimateDT\nattribution:")

    plt.suptitle(
        f"Climate Change Attribution — Flood Impact Metrics\n"
        f"{EVENT_DISPLAY}  |  Factual vs Clausius-Clapeyron & ClimateDT counterfactuals",
        fontsize=13, fontweight="bold", y=1.01,
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def export_summary(metrics: dict, output_path: Path):
    """Print attribution table to console and write it to CSV."""
    fct = metrics[FACTUAL_KEY]
    cc  = metrics[COUNTERFACTUALS[0]]
    cdt = metrics[COUNTERFACTUALS[1]]

    col_specs = [
        ("total_damage",            "Total Damage [USD]",          lambda v: f"${v:,.0f}"),
        ("n_affected",              "Affected Buildings",           lambda v: f"{v:,.0f}"),
        ("mean_damage_per_building","Mean Damage/Building [USD]",  lambda v: f"${v:,.0f}"),
        ("pop_exposed",             "Population Exposed",           lambda v: f"{v:,.0f}"),
    ]

    w = 100
    print("\n" + "=" * w)
    print(f"FIAT CLIMATE ATTRIBUTION SUMMARY  —  {EVENT_NAME}")
    print("=" * w)
    hdr = (f"{'Metric':<32} {'Factual':>13} {'CC (CF-8)':>13} {'CDT (CF-10)':>13}"
           f" {'Attr CC%':>9} {'Attr CDT%':>10}")
    print(hdr)
    print("-" * w)

    rows = []
    for col, label, fmt in col_specs:
        v_fct = fct[col]
        v_cc  = cc[col]
        v_cdt = cdt[col]
        pct_cc  = (v_fct - v_cc)  / v_fct * 100 if v_fct != 0 else 0.0
        pct_cdt = (v_fct - v_cdt) / v_fct * 100 if v_fct != 0 else 0.0
        print(f"{label:<32} {fmt(v_fct):>13} {fmt(v_cc):>13} {fmt(v_cdt):>13}"
              f" {pct_cc:>8.1f}% {pct_cdt:>9.1f}%")
        rows.append({
            "Metric": label,
            "Factual_CF0": v_fct,
            "Clausius_Clapeyron_CF-8": v_cc,
            "ClimateDT_CF-10": v_cdt,
            "Attribution_CC_pct": round(pct_cc, 2),
            "Attribution_CDT_pct": round(pct_cdt, 2),
        })

    print("=" * w + "\n")
    pd.DataFrame(rows).to_csv(output_path, index=False)
    print(f"Saved: {output_path}")


# ===== MAIN =====
def main():
    print(f"\n{'='*60}")
    print(f"FIAT Climate Attribution Analysis — {EVENT_NAME}")
    print(f"{'='*60}\n")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading FIAT spatial outputs...")
    gdfs = {}
    metrics = {}
    for name in RUNS:
        try:
            print(f"  Loading {name}...")
            gdf = load_spatial_fgb(name)
            gdfs[name] = gdf
            metrics[name] = compute_impact_metrics(gdf)
            m = metrics[name]
            print(
                f"    Total damage: ${m['total_damage']:,.0f} | "
                f"Affected: {m['n_affected']:,} | "
                f"Population exposed: {m['pop_exposed']:,.0f}"
            )
        except Exception as e:
            print(f"  ERROR loading {name}: {e}")
            return

    # Common spatial bounds across all three scenarios
    bounds = _common_bounds(gdfs)

    print("\nGenerating plots...")

    for cft_key in COUNTERFACTUALS:
        slug = cft_key.lower().replace(" ", "_").replace("(", "").replace(")", "").replace("-", "")
        print(f"  Damage map: Factual vs {cft_key}...")
        plot_damage_spatial_comparison(
            gdfs[FACTUAL_KEY], gdfs[cft_key], cft_key, bounds,
            OUTPUT_DIR / f"fiat_damage_spatial_{slug}_{EVENT_NAME.lower()}.png",
        )
        print(f"  Population map: Factual vs {cft_key}...")
        plot_pop_spatial_comparison(
            gdfs[FACTUAL_KEY], gdfs[cft_key], cft_key, bounds,
            OUTPUT_DIR / f"fiat_pop_spatial_{slug}_{EVENT_NAME.lower()}.png",
        )

    print("  Damage attribution comparison (CC vs CDT side-by-side)...")
    plot_damage_attribution_comparison(
        gdfs, bounds,
        OUTPUT_DIR / f"fiat_damage_attribution_comparison_{EVENT_NAME.lower()}.png",
    )
    print("  Population attribution comparison (CC vs CDT side-by-side)...")
    plot_pop_attribution_comparison(
        gdfs, bounds,
        OUTPUT_DIR / f"fiat_pop_attribution_comparison_{EVENT_NAME.lower()}.png",
    )

    print("  Creating attribution metrics bar charts...")
    plot_attribution_metrics(
        metrics,
        OUTPUT_DIR / f"fiat_metrics_attribution_{EVENT_NAME.lower()}.png",
    )

    print("  Exporting summary...")
    export_summary(
        metrics,
        OUTPUT_DIR / f"fiat_summary_attribution_{EVENT_NAME.lower()}.csv",
    )

    print(f"\n{'='*60}")
    print(f"Analysis complete!  Outputs saved to: {OUTPUT_DIR}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()

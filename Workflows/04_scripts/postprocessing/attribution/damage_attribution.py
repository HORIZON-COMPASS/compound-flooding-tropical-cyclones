"""
FIAT Economic Damage Attribution Analysis — generic multi-event script.

Compares building-level economic damage across N scenarios (factual + counterfactuals)
for any supported event. Produces:
  1. Spatial damage maps (rasterized or scatter) per factual–counterfactual pair
  2. Climate-attributable damage difference maps (shared colour scale across counterfactuals)
  3. Attribution metrics bar chart (total damage, affected buildings, mean damage, population)
  4. Summary CSV

Select the event by setting EVENT_NAME and DAMAGE_COLUMN below.

RASTERIZE = True  → aggregate buildings onto a regular lat/lon grid (better for dense cities)
RASTERIZE = False → scatter each building as a point (better for sparse rural areas)

Population columns require `fiat_add_pop_exposed_metric.py` to have been run first
to produce `spatial_with_pop_and_flood.fgb`. The script falls back to `spatial.fgb`
and skips population figures gracefully if the enriched file is absent.

Run:
    pixi run -e compass-v1 python attribution/damage_attribution.py

Supported events
----------------
  Durban_April2022   3 scenarios: factual vs Clausius-Clapeyron & ClimateDT
  Freddy             2 scenarios: factual vs CC −8%
  Kenneth            2 scenarios: factual vs CC −8%
  Idai               2 scenarios: factual vs counterfactual
"""

import warnings
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.colors import BoundaryNorm

try:
    import contextily as ctx
    _HAS_CTX = True
except ImportError:
    _HAS_CTX = False

warnings.filterwarnings("ignore")

# ── SELECT EVENT ──────────────────────────────────────────────────────────────
EVENT_NAME     = "Durban_April2022"
DAMAGE_COLUMN  = "total_damage"   # "total_damage" or "relative_damage"
RASTERIZE      = True             # True → grid, False → scatter points

BASE_RUN_PATH  = Path("/p/11210471-001-compass/03_Runs")
GRID_RESOLUTION = 0.01            # degrees (~1 km at Durban, ~0.7 km at Somerset)

# ── EVENT CONFIGURATION ───────────────────────────────────────────────────────
# path_pattern resolves to:
#   BASE_RUN_PATH / region / runname / fiat / folder / output / fgb_file
#
# Each scenario entry:
#   runname    : subdirectory under region/
#   folder     : FIAT event subfolder under runname/fiat/
#   color      : bar chart colour
#   label      : display label
#   is_factual : True for exactly one entry per event

EVENT_CONFIG = {
    "Durban_April2022": {
        "display_name": "Durban April 2022",
        "region":       "durban",
        "output_dir":   Path("/p/11210471-001-compass/04_Results/fiat_climate_attribution"),
        "fgb_file":     "spatial.fgb",
        "scenarios": {
            "Factual (CF0)": {
                "runname":    "Durban2022",
                "folder":     "event_precip_era5_hourly_CF0_no_wind",
                "color":      "#FF6B6B",
                "label":      "Factual\n(ERA5)",
                "is_factual": True,
            },
            "Clausius-Clapeyron": {
                "runname": "Durban2022",
                "folder":  "event_precip_era5_hourly_CF-8_no_wind",
                "color":   "#4ECDC4",
                "label":   "CC\n(−8%)",
            },
            "ClimateDT": {
                "runname": "Durban2022",
                "folder":  "event_precip_era5_hourly_CF-10_no_wind",
                "color":   "#FFE66D",
                "label":   "ClimateDT\n(−10%)",
            },
        },
    },

    "Freddy": {
        "display_name": "Tropical Cyclone Freddy",
        "region":       "test",
        "output_dir":   Path("/p/11210471-001-compass/04_Results/CF_figs"),
        "fgb_file":     "output_relative_damage.fgb",
        "scenarios": {
            "Factual": {
                "runname":    "Freddy",
                "folder":     "event_tp_era5_hourly_CF0_GTSMv41opendap_CF0_no_wind_CF0",
                "color":      "#4472C4",
                "label":      "Factual",
                "is_factual": True,
            },
            "Counterfactual": {
                "runname": "Freddy",
                "folder":  "event_tp_era5_hourly_CF-8_GTSMv41opendap_CF0_no_wind_CF0",
                "color":   "#ED7D31",
                "label":   "Counterfactual\n(CC −8%)",
            },
        },
    },

    "Kenneth": {
        "display_name": "Tropical Cyclone Kenneth",
        "region":       "test",
        "output_dir":   Path("/p/11210471-001-compass/04_Results/CF_figs"),
        "fgb_file":     "output_relative_damage.fgb",
        "scenarios": {
            "Factual": {
                "runname":    "Kenneth",
                "folder":     "event_tp_era5_hourly_zarr_CF0_GTSMv41opendap_CF0_no_wind_CF0",
                "color":      "#4472C4",
                "label":      "Factual",
                "is_factual": True,
            },
            "Counterfactual": {
                "runname": "Kenneth",
                "folder":  "event_tp_era5_hourly_zarr_CF-8_GTSMv41opendap_CF0_no_wind_CF0",
                "color":   "#ED7D31",
                "label":   "Counterfactual\n(CC −8%)",
            },
        },
    },

    "Idai": {
        "display_name": "Tropical Cyclone Idai",
        "region":       "sofala",
        "output_dir":   Path("/p/11210471-001-compass/04_Results/CF_figs"),
        "fgb_file":     "spatial.fgb",
        "scenarios": {
            "Factual": {
                "runname":    "Idai",
                "folder":     "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0",
                "color":      "#4472C4",
                "label":      "Factual",
                "is_factual": True,
            },
            "Counterfactual": {
                "runname": "Idai",
                "folder":  "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10",
                "color":   "#ED7D31",
                "label":   "Counterfactual\n(CC −8%)",
            },
        },
    },
}

# ── DERIVED CONFIG ────────────────────────────────────────────────────────────
cfg          = EVENT_CONFIG[EVENT_NAME]
DISPLAY_NAME = cfg["display_name"]
REGION       = cfg["region"]
SCENARIOS    = cfg["scenarios"]
OUTPUT_DIR   = cfg["output_dir"]
FGB_FILE     = cfg["fgb_file"]

FACTUAL_KEY     = next(k for k, v in SCENARIOS.items() if v.get("is_factual"))
COUNTERFACTUALS = [k for k in SCENARIOS if k != FACTUAL_KEY]

ZOOM_LEVEL = 12


# ── DATA LOADING ──────────────────────────────────────────────────────────────
def fgb_path(name: str) -> Path:
    sc = SCENARIOS[name]
    base = BASE_RUN_PATH / REGION / sc["runname"] / "fiat" / sc["folder"] / "output"
    enriched = base / "spatial_with_pop_and_flood.fgb"
    return enriched if enriched.exists() else base / FGB_FILE


def load_gdf(name: str) -> gpd.GeoDataFrame:
    path = fgb_path(name)
    if not path.exists():
        raise FileNotFoundError(f"FIAT output not found:\n  {path}")
    gdf = gpd.read_file(path)
    if str(gdf.crs) != "EPSG:4326":
        gdf = gdf.to_crs("EPSG:4326")
    centroids = gdf.geometry.centroid
    gdf = gdf.copy()
    gdf["x"] = centroids.x
    gdf["y"] = centroids.y
    return gdf


def compute_metrics(gdf: gpd.GeoDataFrame) -> dict:
    affected = gdf[gdf[DAMAGE_COLUMN] > 0] if DAMAGE_COLUMN in gdf.columns else gdf.iloc[:0]
    total_dam = float(affected[DAMAGE_COLUMN].sum()) if DAMAGE_COLUMN in affected.columns else 0.0
    pop = float(affected["population"].sum()) if "population" in affected.columns else 0.0
    return {
        "total_damage":             total_dam,
        "n_affected":               int(len(affected)),
        "mean_damage_per_building": float(affected[DAMAGE_COLUMN].mean()) if len(affected) > 0 else 0.0,
        "pop_exposed":              pop,
    }


# ── RASTERIZATION ─────────────────────────────────────────────────────────────
def _common_bounds(gdfs: dict) -> tuple:
    mins, maxs = {"x": [], "y": []}, {"x": [], "y": []}
    for gdf in gdfs.values():
        sub = gdf[gdf[DAMAGE_COLUMN] > 0] if DAMAGE_COLUMN in gdf.columns else gdf
        if len(sub) == 0:
            continue
        for ax in ("x", "y"):
            mins[ax].append(sub[ax].min())
            maxs[ax].append(sub[ax].max())
    return (min(mins["x"]), min(mins["y"]), max(maxs["x"]), max(maxs["y"]))


def rasterize(gdf: gpd.GeoDataFrame, col: str, bounds: tuple) -> xr.DataArray:
    sub = gdf[gdf[col] > 0].copy() if col in gdf.columns else gdf.iloc[:0].copy()
    if len(sub) == 0:
        return None
    xs, ys = sub["x"].values, sub["y"].values
    vals   = sub[col].values
    minx, miny, maxx, maxy = bounds
    r = GRID_RESOLUTION
    x_edges = np.arange(minx - r / 2, maxx + r, r)
    y_edges = np.arange(miny - r / 2, maxy + r, r)
    xc = (x_edges[:-1] + x_edges[1:]) / 2
    yc = (y_edges[:-1] + y_edges[1:]) / 2
    xi = np.clip(np.digitize(xs, x_edges) - 1, 0, len(xc) - 1)
    yi = np.clip(np.digitize(ys, y_edges) - 1, 0, len(yc) - 1)
    grid  = np.zeros((len(yc), len(xc)))
    count = np.zeros_like(grid)
    for i in range(len(vals)):
        if not np.isnan(vals[i]):
            grid[yi[i], xi[i]]  += vals[i]
            count[yi[i], xi[i]] += 1
    return xr.DataArray(
        np.where(count > 0, grid, np.nan), dims=["y", "x"],
        coords={"y": yc, "x": xc},
    )


# ── BASEMAP HELPER ────────────────────────────────────────────────────────────
def _basemap(ax):
    if _HAS_CTX:
        try:
            ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik,
                            crs="EPSG:4326", attribution=False, zorder=1, zoom=ZOOM_LEVEL)
        except Exception:
            pass


# ── PLOT: SPATIAL MAPS PER SCENARIO PAIR ──────────────────────────────────────
def plot_spatial_comparison(gdfs: dict, bounds: tuple):
    cmap_val  = plt.get_cmap("Reds")
    cmap_diff = plt.get_cmap("RdBu_r")
    fct_gdf   = gdfs[FACTUAL_KEY]

    for cft_key in COUNTERFACTUALS:
        cft_gdf    = gdfs[cft_key]
        cft_label  = SCENARIOS[cft_key]["label"].replace("\n", " ")
        slug = cft_key.lower().replace(" ", "_").replace("(", "").replace(")", "").replace("-", "")
        out  = OUTPUT_DIR / f"damage_spatial_{slug}_{EVENT_NAME.lower()}.png"

        if RASTERIZE:
            da_fct = rasterize(fct_gdf, DAMAGE_COLUMN, bounds)
            da_cft = rasterize(cft_gdf, DAMAGE_COLUMN, bounds)
            if da_fct is None or da_cft is None:
                print(f"  Insufficient data for {cft_key}, skipping spatial map.")
                continue
            da_fct_f = da_fct.fillna(0)
            da_cft_f = da_cft.fillna(0)
            da_diff  = (da_fct_f - da_cft_f).where((da_fct_f > 0) | (da_cft_f > 0))
            all_vals = np.concatenate([
                da_fct.values[~np.isnan(da_fct.values)],
                da_cft.values[~np.isnan(da_cft.values)],
            ])
            vmax  = float(np.quantile(all_vals, 0.95)) if len(all_vals) > 0 else 1.0
            vnorm = BoundaryNorm(np.linspace(0, vmax, 11), ncolors=256, clip=True)
            dvals = da_diff.values[~np.isnan(da_diff.values)]
            dabs  = max(float(np.quantile(np.abs(dvals), 0.95)), 1e-6) if len(dvals) > 0 else 1.0
            dnorm = BoundaryNorm(np.linspace(-dabs, dabs, 11), ncolors=256, clip=True)

            fig, axes = plt.subplots(1, 3, figsize=(20, 6))
            for ax, (da, norm, cmap, title, clabel) in zip(axes, [
                (da_fct,  vnorm, cmap_val,  "Factual",            f"Damage per {GRID_RESOLUTION}° cell [USD]"),
                (da_cft,  vnorm, cmap_val,  cft_label,            f"Damage per {GRID_RESOLUTION}° cell [USD]"),
                (da_diff, dnorm, cmap_diff, f"Factual − {cft_label}", f"Damage difference [USD]"),
            ]):
                im = da.plot(ax=ax, cmap=cmap, norm=norm, add_colorbar=False,
                             x="x", y="y", alpha=0.85, zorder=2)
                ax.set_aspect("equal")
                ax.set_title(title, fontsize=12, fontweight="bold")
                ax.set_xlabel("Longitude [°]", fontsize=10)
                ax.set_ylabel("Latitude [°]", fontsize=10)
                _basemap(ax)
                plt.colorbar(im, ax=ax, shrink=0.8, pad=0.1).set_label(clabel, fontsize=9)
        else:
            merged_fct  = fct_gdf[fct_gdf[DAMAGE_COLUMN] > 0]
            merged_cft  = cft_gdf[cft_gdf[DAMAGE_COLUMN] > 0]
            all_vals = np.concatenate([merged_fct[DAMAGE_COLUMN].values, merged_cft[DAMAGE_COLUMN].values])
            vmax  = float(np.quantile(all_vals, 0.95)) if len(all_vals) > 0 else 1.0
            vnorm = BoundaryNorm(np.linspace(0, vmax, 11), ncolors=256, clip=True)

            merged = merged_fct[["object_id", "x", "y", DAMAGE_COLUMN]].merge(
                merged_cft[["object_id", "x", "y", DAMAGE_COLUMN]],
                on="object_id", how="outer", suffixes=("_fct", "_cft"),
            )
            merged["x"] = merged["x_fct"].fillna(merged["x_cft"])
            merged["y"] = merged["y_fct"].fillna(merged["y_cft"])
            merged[f"{DAMAGE_COLUMN}_fct"] = merged[f"{DAMAGE_COLUMN}_fct"].fillna(0)
            merged[f"{DAMAGE_COLUMN}_cft"] = merged[f"{DAMAGE_COLUMN}_cft"].fillna(0)
            merged["diff"] = merged[f"{DAMAGE_COLUMN}_fct"] - merged[f"{DAMAGE_COLUMN}_cft"]
            dabs = max(float(merged["diff"].abs().quantile(0.95)), 1e-6)
            dnorm = BoundaryNorm(np.linspace(-dabs, dabs, 11), ncolors=256, clip=True)

            fig, axes = plt.subplots(1, 3, figsize=(20, 6))
            for ax, (col, norm, cmap, title) in zip(axes, [
                (f"{DAMAGE_COLUMN}_fct", vnorm,  cmap_val,  "Factual"),
                (f"{DAMAGE_COLUMN}_cft", vnorm,  cmap_val,  cft_label),
                ("diff",                 dnorm,  cmap_diff, f"Factual − {cft_label}"),
            ]):
                sc = ax.scatter(merged["x"], merged["y"], c=merged[col],
                                cmap=cmap, norm=norm, s=18, alpha=0.8, zorder=2)
                ax.set_aspect("equal")
                ax.set_title(title, fontsize=12, fontweight="bold")
                ax.set_xlabel("Longitude [°]", fontsize=10)
                ax.set_ylabel("Latitude [°]", fontsize=10)
                _basemap(ax)
                plt.colorbar(sc, ax=ax, shrink=0.8, pad=0.1).set_label(DAMAGE_COLUMN, fontsize=9)

        plt.suptitle(f"Economic Damage — {DISPLAY_NAME}\nFactual vs {cft_label}",
                     fontsize=14, fontweight="bold", y=1.02)
        plt.tight_layout()
        plt.savefig(out, dpi=200, bbox_inches="tight")
        plt.close()
        print(f"Saved: {out}")


# ── PLOT: ATTRIBUTION COMPARISON MAP (shared scale across counterfactuals) ────
def plot_attribution_map(gdfs: dict, bounds: tuple):
    if len(COUNTERFACTUALS) < 2:
        return  # nothing to compare when only one counterfactual
    if not RASTERIZE:
        return  # scatter mode doesn't easily support shared-scale attribution maps

    da_fct = rasterize(gdfs[FACTUAL_KEY], DAMAGE_COLUMN, bounds)
    if da_fct is None:
        return
    da_fct_f = da_fct.fillna(0)

    diffs = {}
    for cft in COUNTERFACTUALS:
        da_cft = rasterize(gdfs[cft], DAMAGE_COLUMN, bounds)
        if da_cft is None:
            continue
        da_cft_f   = da_cft.fillna(0)
        diffs[cft] = (da_fct_f - da_cft_f).where((da_fct_f > 0) | (da_cft_f > 0))

    if not diffs:
        return

    all_diffs = np.concatenate([d.values[~np.isnan(d.values)] for d in diffs.values()])
    dabs      = max(float(np.quantile(np.abs(all_diffs), 0.95)), 1e-6) if len(all_diffs) > 0 else 1.0
    dnorm     = BoundaryNorm(np.linspace(-dabs, dabs, 11), ncolors=256, clip=True)
    cmap_diff = plt.get_cmap("RdBu_r")

    n   = len(diffs)
    fig, axes = plt.subplots(1, n, figsize=(8 * n, 6))
    if n == 1:
        axes = [axes]
    ims = []
    for ax, cft in zip(axes, diffs):
        cft_label = SCENARIOS[cft]["label"].replace("\n", " ")
        im = diffs[cft].plot(ax=ax, cmap=cmap_diff, norm=dnorm, add_colorbar=False,
                             x="x", y="y", alpha=0.85, zorder=2)
        ims.append(im)
        ax.set_aspect("equal")
        ax.set_title(f"Factual − {cft_label}", fontsize=12, fontweight="bold")
        ax.set_xlabel("Longitude [°]", fontsize=10)
        ax.set_ylabel("Latitude [°]", fontsize=10)
        _basemap(ax)

    plt.suptitle(f"Climate-Attributable Damage — {DISPLAY_NAME}", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout(rect=[0, 0, 0.88, 1])
    cax = fig.add_axes([0.91, 0.15, 0.02, 0.68])
    fig.colorbar(ims[0], cax=cax).set_label(f"Damage per {GRID_RESOLUTION}° cell [USD]", fontsize=9)
    out = OUTPUT_DIR / f"damage_attribution_comparison_{EVENT_NAME.lower()}.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


# ── PLOT: METRICS BAR CHART ───────────────────────────────────────────────────
def _draw_bracket(ax, x_line, v_fct, v_cft, ymax, label):
    if v_fct == 0 or v_cft == 0 or abs(v_fct - v_cft) < 1e-9:
        return
    lo, hi = min(v_fct, v_cft), max(v_fct, v_cft)
    ext = 0.05
    ax.plot([x_line, x_line], [lo, hi], "k-", linewidth=1.6)
    for y in [lo, hi]:
        ax.plot([x_line - ext, x_line + ext], [y, y], "k-", linewidth=1.2)
    for val in [v_fct, v_cft]:
        ax.plot([ax.get_xlim()[0], x_line], [val, val],
                linestyle="--", color="gray", alpha=0.35, linewidth=0.7)
    attr = v_fct - v_cft
    pct  = attr / v_cft * 100 if v_cft else 0.0
    sign = "+" if attr >= 0 else ""
    ax.text(x_line + 0.07, (lo + hi) / 2,
            f"{sign}{pct:.1f}%", va="center", ha="left", fontsize=8, style="italic")
    ax.text(x_line + 0.07, hi + ymax * 0.04,
            label, va="bottom", ha="left", fontsize=7, style="italic", color="gray")


def plot_metrics_barchart(metrics: dict):
    n      = len(SCENARIOS)
    names  = list(SCENARIOS.keys())
    colors = [SCENARIOS[k]["color"] for k in names]
    labels = [SCENARIOS[k]["label"] for k in names]
    x      = list(range(n))

    specs = [
        ("total_damage",             "Total Damage [USD]",          1e6,  "M"),
        ("n_affected",               "Affected Buildings [#]",      1,    ""),
        ("mean_damage_per_building", "Mean Damage per Building [$]", 1e3, "K"),
        ("pop_exposed",              "Population Exposed [#]",       1,    ""),
    ]
    has_pop = any(metrics[k]["pop_exposed"] > 0 for k in metrics)
    if not has_pop:
        specs = specs[:3]
        n_plots = 3
        figsize = (16, 5)
    else:
        n_plots = 4
        figsize = (14, 10)

    ncols = 3 if n_plots == 3 else 2
    nrows = 1 if n_plots == 3 else 2
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = np.array(axes).flatten()

    for ax, (col, ylabel, scale, unit) in zip(axes, specs):
        vals = [metrics[k][col] / scale for k in names]
        bars = ax.bar(x, vals, color=colors, edgecolor="black", linewidth=0.8, width=0.5)
        ymax = max(vals) * 1.65 if max(vals) > 0 else 1.0
        ax.set_ylabel(ylabel, fontsize=10, fontweight="bold")
        ax.set_title(ylabel.split(" [")[0], fontsize=11, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylim(0, ymax)
        ax.set_xlim(-0.5, n - 0.5 + 1.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)
        for bar, val in zip(bars, vals):
            if val > 0:
                lbl = f"{val:.1f}{unit}" if unit else f"{val:,.0f}"
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + ymax * 0.01,
                        lbl, ha="center", va="bottom", fontsize=9, fontweight="bold")

        v_fct  = metrics[FACTUAL_KEY][col] / scale
        fct_xi = names.index(FACTUAL_KEY)
        for j, cft in enumerate(COUNTERFACTUALS):
            cft_xi  = names.index(cft)
            x_brace = (fct_xi + cft_xi) / 2 + 0.5 + j * 0.6
            _draw_bracket(ax, x_brace, v_fct, metrics[cft][col] / scale, ymax,
                          f"vs {SCENARIOS[cft]['label'].replace(chr(10), ' ')}")

    cft_str = " & ".join(SCENARIOS[k]["label"].replace("\n", " ") for k in COUNTERFACTUALS)
    plt.suptitle(f"Damage Attribution — {DISPLAY_NAME}\nFactual vs {cft_str}",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    out = OUTPUT_DIR / f"damage_metrics_attribution_{EVENT_NAME.lower()}.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


# ── SUMMARY CSV + PRINT ───────────────────────────────────────────────────────
def export_summary(metrics: dict):
    col_specs = [
        ("total_damage",             "Total Damage [USD]",          lambda v: f"${v:,.0f}"),
        ("n_affected",               "Affected Buildings",           lambda v: f"{v:,}"),
        ("mean_damage_per_building", "Mean Damage/Building [USD]",  lambda v: f"${v:,.0f}"),
        ("pop_exposed",              "Population Exposed",           lambda v: f"{v:,.0f}"),
    ]
    w = 26 + 16 * len(SCENARIOS)
    print("\n" + "=" * w)
    print(f"  DAMAGE ATTRIBUTION SUMMARY — {DISPLAY_NAME}")
    print("=" * w)
    hdr = f"{'Metric':<30}" + "".join(f" {k:>14}" for k in SCENARIOS)
    print(hdr)
    print("-" * w)

    rows = []
    for col, label, fmt in col_specs:
        vals = {k: metrics[k][col] for k in SCENARIOS}
        row  = {"Metric": label}
        line = f"{label:<30}"
        for k in SCENARIOS:
            row[k] = vals[k]
            line  += f" {fmt(vals[k]):>14}"
        print(line)
        for cft in COUNTERFACTUALS:
            v_fct = vals[FACTUAL_KEY]
            v_cft = vals[cft]
            pct   = (v_fct - v_cft) / v_cft * 100 if v_cft else 0.0
            row[f"pct_Factual_vs_{cft.replace(' ', '_')}"] = round(pct, 2)
        rows.append(row)

    print("=" * w + "\n")
    out = OUTPUT_DIR / f"damage_summary_attribution_{EVENT_NAME.lower()}.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"Saved: {out}")


# ── MAIN ──────────────────────────────────────────────────────────────────────
def main():
    print(f"\n{'='*65}")
    print(f"FIAT Damage Attribution — {DISPLAY_NAME}")
    print(f"EVENT_NAME={EVENT_NAME!r}  DAMAGE_COLUMN={DAMAGE_COLUMN!r}  RASTERIZE={RASTERIZE}")
    print(f"{'='*65}\n")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading FIAT outputs...")
    gdfs, metrics = {}, {}
    for name in SCENARIOS:
        print(f"  {name} ...")
        gdf = load_gdf(name)
        gdfs[name]     = gdf
        metrics[name]  = compute_metrics(gdf)
        m = metrics[name]
        print(f"    damage=${m['total_damage']:,.0f}  affected={m['n_affected']:,}  "
              f"pop={m['pop_exposed']:,.0f}")

    bounds = _common_bounds(gdfs)

    print("\nGenerating figures...")
    plot_spatial_comparison(gdfs, bounds)
    plot_attribution_map(gdfs, bounds)
    plot_metrics_barchart(metrics)
    export_summary(metrics)

    print(f"\n{'='*65}")
    print(f"Done. Outputs: {OUTPUT_DIR}")
    print(f"{'='*65}\n")


if __name__ == "__main__":
    main()

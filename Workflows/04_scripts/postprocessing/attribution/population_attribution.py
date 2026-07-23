"""
FIAT Population Exposure Attribution Analysis — generic multi-event script.

Compares building-level population exposure to flooding across N scenarios
(factual + counterfactuals) for any supported event. Requires FIAT output
files enriched with population data — run `fiat_add_pop_exposed_metric.py`
first to produce `spatial_with_pop_and_flood.fgb`.

Produces:
  1. Spatial scatter map of exposed population per scenario pair (factual | CFX | diff)
  2. Bar chart of total exposed population with attribution brackets
  3. Summary CSV

Select the event by setting EVENT_NAME and INUNDATION_THRESHOLD below.

Run:
    pixi run -e compass-v1 python attribution/population_attribution.py

Supported events
----------------
  Freddy    2 scenarios: factual vs CC −8%
  Kenneth   2 scenarios: factual vs CC −8%
  Idai      2 scenarios: factual vs counterfactual
"""

import warnings
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm

try:
    import contextily as ctx
    _HAS_CTX = True
except ImportError:
    _HAS_CTX = False

warnings.filterwarnings("ignore")

# ── SELECT EVENT ──────────────────────────────────────────────────────────────
EVENT_NAME          = "Freddy"
INUNDATION_THRESHOLD = 0.2   # m — only count locations with flood depth > this
POPULATION_COLUMN   = "population"

BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs")
ZOOM_LEVEL    = 12

# ── EVENT CONFIGURATION ───────────────────────────────────────────────────────
EVENT_CONFIG = {
    "Freddy": {
        "display_name": "Tropical Cyclone Freddy",
        "region":       "test",
        "output_dir":   Path("/p/11210471-001-compass/04_Results/CF_figs"),
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

FACTUAL_KEY     = next(k for k, v in SCENARIOS.items() if v.get("is_factual"))
COUNTERFACTUALS = [k for k in SCENARIOS if k != FACTUAL_KEY]


# ── DATA LOADING ──────────────────────────────────────────────────────────────
def _fgb_path(name: str) -> Path:
    sc   = SCENARIOS[name]
    base = BASE_RUN_PATH / REGION / sc["runname"] / "fiat" / sc["folder"] / "output"
    enriched = base / "spatial_with_pop_and_flood.fgb"
    if not enriched.exists():
        raise FileNotFoundError(
            f"Population-enriched FIAT output not found:\n  {enriched}\n"
            "Run fiat_add_pop_exposed_metric.py first."
        )
    return enriched


def load_gdf(name: str) -> gpd.GeoDataFrame:
    path = _fgb_path(name)
    gdf  = gpd.read_file(path)
    if str(gdf.crs) != "EPSG:4326":
        gdf = gdf.to_crs("EPSG:4326")
    centroids    = gdf.geometry.centroid
    gdf          = gdf.copy()
    gdf["x"]     = centroids.x
    gdf["y"]     = centroids.y

    # Resolve population column
    if POPULATION_COLUMN in gdf.columns:
        gdf["_pop"] = gdf[POPULATION_COLUMN]
    elif "pop" in gdf.columns:
        gdf["_pop"] = gdf["pop"]
    elif "total_population" in gdf.columns:
        gdf["_pop"] = gdf["total_population"]
    else:
        raise ValueError(f"No population column found in {path}. Columns: {list(gdf.columns)}")

    # Resolve inundation depth column
    depth_cols = [c for c in gdf.columns if "inun" in c.lower() or "depth" in c.lower()]
    if depth_cols:
        gdf["_depth"] = gdf[depth_cols[0]]
    else:
        gdf["_depth"] = 0.0  # if missing, keep all records

    return gdf[gdf["_depth"] > INUNDATION_THRESHOLD].copy()


def compute_metrics(gdf: gpd.GeoDataFrame) -> dict:
    has_pop = (gdf["_pop"] > 0).any()
    return {
        "total_pop":     float(gdf["_pop"].sum()),
        "n_locations":   int(len(gdf)),
        "mean_pop":      float(gdf[gdf["_pop"] > 0]["_pop"].mean()) if has_pop else 0.0,
    }


# ── BASEMAP ───────────────────────────────────────────────────────────────────
def _basemap(ax):
    if _HAS_CTX:
        try:
            ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik,
                            crs="EPSG:4326", attribution=False, zorder=1, zoom=ZOOM_LEVEL)
        except Exception:
            pass


# ── PLOT: SPATIAL COMPARISON ──────────────────────────────────────────────────
def plot_spatial_comparison(gdfs: dict):
    cmap_pop  = plt.get_cmap("Blues")
    cmap_diff = plt.get_cmap("RdBu_r")
    fct_gdf   = gdfs[FACTUAL_KEY]

    for cft_key in COUNTERFACTUALS:
        cft_gdf   = gdfs[cft_key]
        cft_label = SCENARIOS[cft_key]["label"].replace("\n", " ")
        slug      = cft_key.lower().replace(" ", "_").replace("(", "").replace(")", "").replace("-", "")

        all_pop = np.concatenate([fct_gdf["_pop"].values, cft_gdf["_pop"].values])
        vmax    = float(np.quantile(all_pop[all_pop > 0], 0.95)) if (all_pop > 0).any() else 1.0
        vnorm   = BoundaryNorm(np.linspace(0, vmax, 11), ncolors=256, clip=True)

        merged = fct_gdf[["object_id", "x", "y", "_pop"]].merge(
            cft_gdf[["object_id", "x", "y", "_pop"]],
            on="object_id", how="outer", suffixes=("_fct", "_cft"),
        )
        merged["x"]      = merged["x_fct"].fillna(merged["x_cft"])
        merged["y"]      = merged["y_fct"].fillna(merged["y_cft"])
        merged["_pop_fct"] = merged["_pop_fct"].fillna(0)
        merged["_pop_cft"] = merged["_pop_cft"].fillna(0)
        merged["diff"]     = merged["_pop_fct"] - merged["_pop_cft"]
        dabs  = max(float(merged["diff"].abs().quantile(0.95)), 1e-6)
        dnorm = BoundaryNorm(np.linspace(-dabs, dabs, 11), ncolors=256, clip=True)

        fig, axes = plt.subplots(1, 3, figsize=(20, 6))
        for ax, (col, norm, cmap, title, clabel) in zip(axes, [
            ("_pop_fct", vnorm,  cmap_pop,  "Factual",   "Exposed population"),
            ("_pop_cft", vnorm,  cmap_pop,  cft_label,   "Exposed population"),
            ("diff",     dnorm,  cmap_diff, f"Factual − {cft_label}", "Population difference"),
        ]):
            sc = ax.scatter(merged["x"], merged["y"], c=merged[col],
                            cmap=cmap, norm=norm, s=18, alpha=0.8, zorder=2)
            ax.set_aspect("equal")
            ax.set_title(title, fontsize=12, fontweight="bold")
            ax.set_xlabel("Longitude [°]", fontsize=10)
            ax.set_ylabel("Latitude [°]", fontsize=10)
            _basemap(ax)
            plt.colorbar(sc, ax=ax, shrink=0.8, pad=0.1).set_label(clabel, fontsize=9)

        plt.suptitle(
            f"Population Exposed to Flooding — {DISPLAY_NAME}  (depth > {INUNDATION_THRESHOLD} m)\n"
            f"Factual vs {cft_label}",
            fontsize=14, fontweight="bold", y=1.02,
        )
        plt.tight_layout()
        out = OUTPUT_DIR / f"pop_exposure_spatial_{slug}_{EVENT_NAME.lower()}.png"
        plt.savefig(out, dpi=200, bbox_inches="tight")
        plt.close()
        print(f"Saved: {out}")


# ── PLOT: BAR CHART ───────────────────────────────────────────────────────────
def _bracket(ax, x_pos, v_fct, v_cft, ymax, label):
    if v_fct == 0 or v_cft == 0 or abs(v_fct - v_cft) < 1e-9:
        return
    lo, hi = min(v_fct, v_cft), max(v_fct, v_cft)
    ax.plot([x_pos, x_pos], [lo, hi], "k-", linewidth=1.6)
    for y in [lo, hi]:
        ax.plot([x_pos - 0.05, x_pos + 0.05], [y, y], "k-", linewidth=1.2)
    attr = v_fct - v_cft
    pct  = attr / v_cft * 100 if v_cft else 0.0
    sign = "+" if attr >= 0 else ""
    ax.text(x_pos + 0.08, (lo + hi) / 2,
            f"{sign}{pct:.1f}%", va="center", ha="left", fontsize=9, style="italic")
    ax.text(x_pos + 0.08, hi + ymax * 0.04,
            label, va="bottom", ha="left", fontsize=7, style="italic", color="gray")


def plot_barchart(metrics: dict):
    names  = list(SCENARIOS.keys())
    colors = [SCENARIOS[k]["color"] for k in names]
    labels = [SCENARIOS[k]["label"] for k in names]
    x      = list(range(len(names)))
    vals   = [metrics[k]["total_pop"] for k in names]

    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    bars = ax.bar(x, vals, color=colors, edgecolor="black", linewidth=0.8, width=0.4)
    ymax = max(vals) * 1.6 if max(vals) > 0 else 1.0
    ax.set_ylabel("Total Exposed Population [people]", fontsize=11, fontweight="bold")
    ax.set_title(f"Population Exposure — {DISPLAY_NAME}\n(depth > {INUNDATION_THRESHOLD} m)",
                 fontsize=12, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylim(0, ymax)
    ax.set_xlim(-0.5, len(names) - 0.5 + 1.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3)
    ax.set_axisbelow(True)

    for bar, val in zip(bars, vals):
        if val > 0:
            if val >= 1e6:
                lbl = f"{val/1e6:.2f}M"
            elif val >= 1e3:
                lbl = f"{val/1e3:.1f}K"
            else:
                lbl = f"{val:.0f}"
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + ymax * 0.01,
                    lbl, ha="center", va="bottom", fontsize=10, fontweight="bold")

    v_fct  = metrics[FACTUAL_KEY]["total_pop"]
    fct_xi = names.index(FACTUAL_KEY)
    for j, cft in enumerate(COUNTERFACTUALS):
        cft_xi  = names.index(cft)
        x_brace = (fct_xi + cft_xi) / 2 + 0.5 + j * 0.5
        _bracket(ax, x_brace, v_fct, metrics[cft]["total_pop"], ymax,
                 f"vs {SCENARIOS[cft]['label'].replace(chr(10), ' ')}")

    plt.tight_layout()
    out = OUTPUT_DIR / f"pop_exposure_barchart_{EVENT_NAME.lower()}.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


# ── SUMMARY CSV + PRINT ───────────────────────────────────────────────────────
def export_summary(metrics: dict):
    w = 26 + 16 * len(SCENARIOS)
    print("\n" + "=" * w)
    print(f"  POPULATION EXPOSURE ATTRIBUTION SUMMARY — {DISPLAY_NAME}")
    print(f"  Inundation threshold: {INUNDATION_THRESHOLD} m")
    print("=" * w)
    hdr = f"{'Metric':<28}" + "".join(f" {k:>14}" for k in SCENARIOS)
    print(hdr)
    print("-" * w)

    specs = [
        ("total_pop",   "Total Exposed Population"),
        ("n_locations", "Flooded Locations"),
        ("mean_pop",    "Mean Pop/Location"),
    ]
    rows = []
    for col, label in specs:
        vals = {k: metrics[k][col] for k in SCENARIOS}
        row  = {"Metric": label}
        line = f"{label:<28}"
        for k in SCENARIOS:
            row[k] = vals[k]
            line  += f" {vals[k]:>14,.1f}"
        print(line)
        for cft in COUNTERFACTUALS:
            v_fct = vals[FACTUAL_KEY]
            v_cft = vals[cft]
            pct   = (v_fct - v_cft) / v_cft * 100 if v_cft else 0.0
            row[f"pct_Factual_vs_{cft.replace(' ', '_')}"] = round(pct, 2)
        rows.append(row)

    print("=" * w + "\n")
    out = OUTPUT_DIR / f"pop_exposure_summary_{EVENT_NAME.lower()}.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"Saved: {out}")


# ── MAIN ──────────────────────────────────────────────────────────────────────
def main():
    print(f"\n{'='*65}")
    print(f"Population Exposure Attribution — {DISPLAY_NAME}")
    print(f"EVENT_NAME={EVENT_NAME!r}  threshold={INUNDATION_THRESHOLD} m")
    print(f"{'='*65}\n")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading FIAT population outputs...")
    gdfs, metrics = {}, {}
    for name in SCENARIOS:
        print(f"  {name} ...")
        gdf = load_gdf(name)
        gdfs[name]    = gdf
        metrics[name] = compute_metrics(gdf)
        m = metrics[name]
        print(f"    exposed={m['total_pop']:,.0f} people  locations={m['n_locations']:,}")

    print("\nGenerating figures...")
    plot_spatial_comparison(gdfs)
    plot_barchart(metrics)
    export_summary(metrics)

    print(f"\n{'='*65}")
    print(f"Done. Outputs: {OUTPUT_DIR}")
    print(f"{'='*65}\n")


if __name__ == "__main__":
    main()

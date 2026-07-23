"""
SFINCS Flood Attribution Analysis — generic multi-event script.

Compares maximum flood depth (hmax) across N scenarios (factual + counterfactuals)
for any supported event. Produces:
  1. N-panel hmax map (one panel per scenario, shared colour scale)
  2. Attribution difference maps (factual − each counterfactual), shared diverging scale
  3. Bar chart: flood volume / extent / mean depth across all scenarios
  4. Summary CSV

Select the event by setting EVENT_NAME below. Paths for each scenario's
sfincs_output_hmax_AllTime.tif are resolved as:
    BASE_RUN_PATH / region / runname / sfincs / folder / plot_output / sfincs_output_hmax_AllTime.tif

Run:
    pixi run -e compass-v1 python attribution/flood_attribution.py

Supported events
----------------
  Durban_April2022       Durban 2022 — ERA5 factual vs Clausius-Clapeyron & ClimateDT
  Somerset_Dec2013       Somerset Dec 2013 — factual (CEH-GEAR) vs T1 (past) & T3 (future)
  Freddy                 Tropical Cyclone Freddy — ERA5 factual vs CC counterfactual
  Kenneth                Tropical Cyclone Kenneth — ERA5 factual vs CC counterfactual
  Idai                   Tropical Cyclone Idai — ERA5 factual vs CC counterfactual
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rioxarray as rxr
import xarray as xr

try:
    import contextily as ctx
    _HAS_CTX = True
except ImportError:
    _HAS_CTX = False

warnings.filterwarnings("ignore")

# ── SELECT EVENT ──────────────────────────────────────────────────────────────
EVENT_NAME = "Somerset_Dec2013"

BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs")

# ── EVENT CONFIGURATION ───────────────────────────────────────────────────────
# Each entry defines:
#   display_name : human-readable title used in figure labels
#   region       : subdirectory under BASE_RUN_PATH
#   lat_ref      : reference latitude (°) for degree→metre area conversion
#   output_dir   : where figures and CSV are saved
#   scenarios    : ordered dict — first entry with is_factual=True is the reference
#                  runname : scenario folder under region/
#                  folder  : SFINCS event subfolder under runname/sfincs/
#                  color   : bar/line colour
#                  label   : display label (use \n for line break)

EVENT_CONFIG = {
    "Somerset_Dec2013": {
        "display_name": "Somerset Levels Dec 2013",
        "region":       "somerset",
        "lat_ref":      51.0,
        "output_dir":   Path("/p/11210471-001-compass/04_Results/climate_attribution_somerset"),
        "scenarios": {
            "Factual": {
                "runname":    "SomersetLevels_dec_factual",
                "folder":     "event_tp_ceh_gear_compass_CF0_GTSMv41opendap_CF0_no_wind_CF0",
                "color":      "#2166ac",
                "label":      "Factual\n(observed)",
                "is_factual": True,
            },
            "T1 (past)": {
                "runname": "SomersetLevels_dec_t1",
                "folder":  "event_tp_ceh_t1_compass_CF0_GTSMv41opendap_CF0_no_wind_CF0",
                "color":   "#d6604d",
                "label":   "T1\n(past climate)",
            },
            "T3 (future)": {
                "runname": "SomersetLevels_dec_t3",
                "folder":  "event_tp_ceh_t3_compass_CF0_GTSMv41opendap_CF0_no_wind_CF0",
                "color":   "#4dac26",
                "label":   "T3\n(future climate)",
            },
        },
    },

    "Durban_April2022": {
        "display_name": "Durban April 2022",
        "region":       "durban",
        "lat_ref":      -29.8,
        "output_dir":   Path("/p/11210471-001-compass/04_Results/climate_attribution"),
        "scenarios": {
            "Factual": {
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
        "lat_ref":      -20.0,
        "output_dir":   Path("/p/11210471-001-compass/04_Results/climate_attribution"),
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
        "lat_ref":      -12.0,
        "output_dir":   Path("/p/11210471-001-compass/04_Results/climate_attribution"),
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
        "lat_ref":      -20.0,
        "output_dir":   Path("/p/11210471-001-compass/04_Results/climate_attribution"),
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
FLOOD_THRESHOLD = 0.05   # m — minimum depth counted as flooded
HMAX_VMAX       = 2.5    # m — colour scale max for hmax maps
ZOOM_LEVEL      = 11     # contextily zoom (reduce for large domains)

cfg          = EVENT_CONFIG[EVENT_NAME]
DISPLAY_NAME = cfg["display_name"]
REGION       = cfg["region"]
SCENARIOS    = cfg["scenarios"]
OUTPUT_DIR   = cfg["output_dir"]
LAT_REF      = cfg["lat_ref"]
DX_PER_DEG   = 111_000 * np.cos(np.radians(LAT_REF))
DY_PER_DEG   = 111_000

FACTUAL_KEY     = next(k for k, v in SCENARIOS.items() if v.get("is_factual"))
COUNTERFACTUALS = [k for k in SCENARIOS if k != FACTUAL_KEY]


# ── HELPERS ───────────────────────────────────────────────────────────────────
def hmax_path(name: str) -> Path:
    sc = SCENARIOS[name]
    return (BASE_RUN_PATH / REGION / sc["runname"] / "sfincs"
            / sc["folder"] / "plot_output" / "sfincs_output_hmax_AllTime.tif")


def load_hmax(name: str) -> xr.DataArray:
    path = hmax_path(name)
    if not path.exists():
        raise FileNotFoundError(f"hmax TIF not found:\n  {path}")
    da = rxr.open_rasterio(path)
    if "band" in da.dims:
        da = da.squeeze("band", drop=True)
    if da.rio.crs and str(da.rio.crs) != "EPSG:4326":
        da = da.rio.reproject("EPSG:4326")
    return da.where(da > 0)


def _cell_area_m2(da: xr.DataArray) -> float:
    dx = float(np.abs(da.x.diff("x").median())) * DX_PER_DEG
    dy = float(np.abs(da.y.diff("y").median())) * DY_PER_DEG
    return dx * dy


def flood_volume_m3(hmax: xr.DataArray) -> float:
    return float(hmax.where(hmax > FLOOD_THRESHOLD).sum(skipna=True) * _cell_area_m2(hmax))


def flood_extent_km2(hmax: xr.DataArray) -> float:
    return float((hmax > FLOOD_THRESHOLD).sum() * _cell_area_m2(hmax)) / 1e6


def flood_mean_depth_m(hmax: xr.DataArray) -> float:
    return float(hmax.where(hmax > FLOOD_THRESHOLD).mean(skipna=True))


def _add_basemap(ax):
    if _HAS_CTX:
        try:
            ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik,
                            crs="EPSG:4326", attribution=False,
                            zorder=1, zoom=ZOOM_LEVEL)
        except Exception:
            pass


# ── PLOT: N-PANEL HMAX MAP ────────────────────────────────────────────────────
def plot_hmax_panels(hmax_data: dict, output_path: Path):
    n = len(SCENARIOS)
    flood_levels = np.linspace(0, HMAX_VMAX, 21)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 7))
    if n == 1:
        axes = [axes]

    for ax, name in zip(axes, SCENARIOS):
        da = hmax_data[name].where(hmax_data[name] > FLOOD_THRESHOLD)
        im = da.plot(ax=ax, levels=flood_levels, cmap="Blues",
                     add_colorbar=False, x="x", y="y", alpha=0.85, zorder=2)
        ax.set_aspect("equal")
        ax.set_title(SCENARIOS[name]["label"].replace("\n", " — "),
                     fontsize=13, fontweight="bold")
        ax.set_xlabel("Longitude [°]", fontsize=11)
        ax.set_ylabel("Latitude [°]", fontsize=11)
        _add_basemap(ax)
        fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02).set_label("Max depth [m]", fontsize=10)

    scenario_str = " | ".join(SCENARIOS[k]["label"].replace("\n", " ") for k in SCENARIOS)
    plt.suptitle(f"Maximum Flood Depth — {DISPLAY_NAME}\n{scenario_str}",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# ── PLOT: ATTRIBUTION DIFFERENCE MAPS ─────────────────────────────────────────
def plot_attribution_diffs(hmax_data: dict, output_path: Path):
    hmax_fct = hmax_data[FACTUAL_KEY]
    diffs    = {}
    for cft in COUNTERFACTUALS:
        aligned     = hmax_fct.interp_like(hmax_data[cft], method="nearest")
        diffs[cft]  = aligned - hmax_data[cft]

    all_vals = np.concatenate([
        diffs[k].values[np.isfinite(diffs[k].values)] for k in diffs
    ])
    vlim = float(np.quantile(np.abs(all_vals), 0.95)) if len(all_vals) else 1.0
    vlim = max(vlim, 0.05)
    levels = np.linspace(-vlim, vlim, 21)

    n = len(COUNTERFACTUALS)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 7))
    if n == 1:
        axes = [axes]

    ims = []
    for ax, cft in zip(axes, COUNTERFACTUALS):
        lbl = SCENARIOS[cft]["label"].replace("\n", " ")
        im  = diffs[cft].plot(ax=ax, levels=levels, cmap="RdBu_r",
                               add_colorbar=False, x="x", y="y", alpha=0.85, zorder=2)
        ims.append(im)
        ax.set_aspect("equal")
        fct_lbl = SCENARIOS[FACTUAL_KEY]["label"].replace("\n", " ")
        ax.set_title(f"{fct_lbl} − {lbl}\n+ve = more flooding in factual",
                     fontsize=12, fontweight="bold")
        ax.set_xlabel("Longitude [°]", fontsize=11)
        ax.set_ylabel("Latitude [°]", fontsize=11)
        _add_basemap(ax)

    plt.suptitle(f"Climate-Attributable Flood Depth Change — {DISPLAY_NAME}",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout(rect=[0, 0, 0.9, 1])
    cax = fig.add_axes([0.92, 0.15, 0.015, 0.68])
    fig.colorbar(ims[0], cax=cax).set_label("Depth difference [m]", fontsize=10)
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# ── PLOT: METRICS BAR CHART ───────────────────────────────────────────────────
def _bracket(ax, x_pos, v_a, v_b, ymax, label):
    if abs(v_a - v_b) < 1e-9 or v_b == 0:
        return
    lo, hi = min(v_a, v_b), max(v_a, v_b)
    ax.plot([x_pos, x_pos], [lo, hi], "k-", linewidth=1.6)
    for y in [lo, hi]:
        ax.plot([x_pos - 0.05, x_pos + 0.05], [y, y], "k-", linewidth=1.2)
    pct  = (v_a - v_b) / v_b * 100
    sign = "+" if pct >= 0 else ""
    ax.text(x_pos + 0.08, (lo + hi) / 2,
            f"{sign}{pct:.1f}%", va="center", ha="left", fontsize=8, style="italic")
    ax.text(x_pos + 0.08, hi + ymax * 0.03,
            label, va="bottom", ha="left", fontsize=7, color="gray", style="italic")


def plot_metrics_barchart(metrics: dict, output_path: Path):
    names  = list(SCENARIOS.keys())
    colors = [SCENARIOS[k]["color"] for k in names]
    labels = [SCENARIOS[k]["label"] for k in names]
    x      = list(range(len(names)))

    specs = [
        ("volume",     "Flood Volume [Mm³]",  1e6),
        ("extent",     "Flood Extent [km²]",  1),
        ("mean_depth", "Mean Depth [m]",       1),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))

    for ax, (col, ylabel, scale) in zip(axes, specs):
        vals = [metrics[k][col] / scale for k in names]
        bars = ax.bar(x, vals, color=colors, edgecolor="black",
                      linewidth=0.8, width=0.5)
        ax.set_ylabel(ylabel, fontsize=11, fontweight="bold")
        ax.set_title(ylabel.split(" [")[0], fontsize=12, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=9)
        ymax = max(vals) * 1.6 if max(vals) > 0 else 1.0
        ax.set_ylim(0, ymax)
        ax.set_xlim(-0.5, len(names) - 0.5 + 1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)

        for bar, val in zip(bars, vals):
            fmt = f"{val:.3f}" if col == "mean_depth" else f"{val:.2f}"
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + ymax * 0.01,
                    fmt, ha="center", va="bottom", fontsize=9, fontweight="bold")

        v_fct = metrics[FACTUAL_KEY][col] / scale
        fct_i = names.index(FACTUAL_KEY)
        for j, cft in enumerate(COUNTERFACTUALS):
            cft_i   = names.index(cft)
            x_brace = (fct_i + cft_i) / 2 + 0.4
            _bracket(ax, x_brace, v_fct,
                     metrics[cft][col] / scale, ymax,
                     f"Factual vs\n{SCENARIOS[cft]['label'].replace(chr(10),' ')}")

    cft_str = " & ".join(SCENARIOS[k]["label"].replace("\n", " ") for k in COUNTERFACTUALS)
    plt.suptitle(f"Flood Attribution — {DISPLAY_NAME}\n"
                 f"Factual vs {cft_str}",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# ── SUMMARY CSV ───────────────────────────────────────────────────────────────
def export_summary(metrics: dict, output_path: Path):
    rows = []
    v_fct_dict = {col: metrics[FACTUAL_KEY][col] for col in ("volume", "extent", "mean_depth")}

    specs = [
        ("volume",     "Flood Volume [Mm³]", 1e6, ".2f"),
        ("extent",     "Flood Extent [km²]", 1,   ".2f"),
        ("mean_depth", "Mean Depth [m]",      1,   ".4f"),
    ]

    header = f"{'Metric':<22}" + "".join(f" {k:>14}" for k in SCENARIOS)
    print("\n" + "=" * (22 + 15 * len(SCENARIOS)))
    print(f"  SFINCS FLOOD ATTRIBUTION — {DISPLAY_NAME}")
    print("=" * (22 + 15 * len(SCENARIOS)))
    print(header)
    print("-" * (22 + 15 * len(SCENARIOS)))

    for col, label, scale, fmt in specs:
        vals = {k: metrics[k][col] / scale for k in SCENARIOS}
        row  = {"Metric": label}
        line = f"{label:<22}"
        for k in SCENARIOS:
            v = vals[k]
            row[k] = v
            line   += f" {f'{v:{fmt}}':>14}"
        print(line)
        for cft in COUNTERFACTUALS:
            v_fct = vals[FACTUAL_KEY]
            v_cft = vals[cft]
            pct   = (v_fct - v_cft) / v_cft * 100 if v_cft else 0.0
            row[f"pct_Factual_vs_{cft.replace(' ', '_')}"] = round(pct, 2)
        rows.append(row)

    print("=" * (22 + 15 * len(SCENARIOS)) + "\n")
    pd.DataFrame(rows).to_csv(output_path, index=False)
    print(f"Saved: {output_path}")


# ── MAIN ──────────────────────────────────────────────────────────────────────
def main():
    print(f"\n{'='*65}")
    print(f"SFINCS Flood Attribution — {DISPLAY_NAME}")
    print(f"EVENT_NAME = {EVENT_NAME!r}")
    print(f"{'='*65}\n")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading hmax rasters...")
    hmax_data, metrics = {}, {}
    for name in SCENARIOS:
        print(f"  {name} ...")
        hmax = load_hmax(name)
        hmax_data[name] = hmax
        metrics[name] = {
            "volume":     flood_volume_m3(hmax),
            "extent":     flood_extent_km2(hmax),
            "mean_depth": flood_mean_depth_m(hmax),
        }
        m = metrics[name]
        print(f"    vol={m['volume']/1e6:.3f} Mm³  ext={m['extent']:.2f} km²  "
              f"depth={m['mean_depth']:.4f} m")

    slug = EVENT_NAME.lower()
    print("\nGenerating figures...")
    plot_hmax_panels(hmax_data,
                     OUTPUT_DIR / f"hmax_panels_{slug}.png")
    plot_attribution_diffs(hmax_data,
                           OUTPUT_DIR / f"hmax_attribution_diff_{slug}.png")
    plot_metrics_barchart(metrics,
                          OUTPUT_DIR / f"metrics_attribution_{slug}.png")
    export_summary(metrics,
                   OUTPUT_DIR / f"summary_attribution_{slug}.csv")

    print(f"\n{'='*65}")
    print(f"Done. Outputs: {OUTPUT_DIR}")
    print(f"{'='*65}\n")


if __name__ == "__main__":
    main()

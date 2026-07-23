"""
SFINCS Flood Time Series — 3-panel plot from sfincs_map.nc.

Computes at every output timestep:
  (a) Flood extent   [km²] — area where depth > FLOOD_THRESHOLD
  (b) Flood volume   [Mm³] — sum of depth × cell_area
  (c) Depth stats    [m]   — mean, 95th percentile, and max over flooded cells

Water depth:
    h = max(zs − zb, 0)
where zb (n, m) is the subgrid bed level and zs (time, n, m) is the water
surface elevation output by SFINCS.

Set USE_HOURLY = False (default) to use zsmax (88 daily-max steps per 3-month
run) and keep memory low.  Set True to use zs (2112 hourly steps).

Compares multiple scenarios on the same axes when SCENARIOS has >1 entry.

Run:
    pixi run -e compass-v1 python sfincs/plot_sfincs_timeseries.py
"""

import glob
import warnings
from pathlib import Path

import cftime
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

warnings.filterwarnings("ignore")

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
BASE  = Path("/p/11210471-001-compass/03_Runs/somerset")
OUT   = Path("/p/11210471-001-compass/04_Results/climate_attribution_somerset")

# Scenarios to plot.  Each entry: (display_label, hex_colour).
# Set to a single entry for a single-scenario plot.
SCENARIOS = {
    "Factual":     ("SomersetLevels_dec_factual", "#2166ac"),
    "T1 (past)":   ("SomersetLevels_dec_t1",      "#d6604d"),
    "T3 (future)": ("SomersetLevels_dec_t3",       "#4dac26"),
}

FLOOD_THRESHOLD = 0.05   # m — minimum depth to count as flooded
USE_HOURLY      = False  # True → hourly zs (more memory); False → daily zsmax


# ── HELPERS ───────────────────────────────────────────────────────────────────
def _find_nc(scenario_dir: Path) -> str:
    matches = glob.glob(str(scenario_dir / "sfincs/event_tp_*/sfincs_map.nc"))
    if not matches:
        raise FileNotFoundError(
            f"No sfincs_map.nc found under:\n  {scenario_dir / 'sfincs/event_tp_*'}"
        )
    return matches[0]


def _decode_times(ds: xr.Dataset) -> list:
    if USE_HOURLY:
        raw, units = ds["time"].values,    ds["time"].attrs["units"]
        calendar   = ds["time"].attrs.get("calendar", "standard")
    else:
        raw, units = ds["timemax"].values, ds["timemax"].attrs["units"]
        calendar   = ds["timemax"].attrs.get("calendar", "standard")
    return cftime.num2pydate(raw, units=units, calendar=calendar)


def _cell_areas(ds: xr.Dataset) -> np.ndarray:
    cx = ds["corner_x"].values   # (n+1, m+1)
    cy = ds["corner_y"].values
    v1x = cx[:-1, 1:]  - cx[:-1, :-1]
    v1y = cy[:-1, 1:]  - cy[:-1, :-1]
    v2x = cx[1:,  :-1] - cx[:-1, :-1]
    v2y = cy[1:,  :-1] - cy[:-1, :-1]
    return np.abs(v1x * v2y - v1y * v2x)   # (n, m)  m²


def load_scenario(scenario_dir: Path) -> tuple:
    """Return (times, extent_km2, volume_Mm3, mean_h, p95_h, max_h)."""
    nc = _find_nc(scenario_dir)
    print(f"    {nc}")
    ds = xr.open_dataset(nc, decode_times=False)

    msk    = ds["msk"].values          # (n, m)
    zb     = ds["zb"].values           # (n, m)
    area   = _cell_areas(ds)           # (n, m)  m²
    active = msk > 0
    times  = _decode_times(ds)

    zs_data = (ds["zs"] if USE_HOURLY else ds["zsmax"]).values   # (T, n, m)
    ds.close()

    zs_clean = np.where(np.isnan(zs_data), zb[np.newaxis], zs_data)
    depth    = np.maximum(zs_clean - zb[np.newaxis], 0.0)
    depth[:, ~active] = 0.0
    flooded  = depth > FLOOD_THRESHOLD

    n_steps    = depth.shape[0]
    extent_km2 = np.empty(n_steps)
    volume_Mm3 = np.empty(n_steps)
    mean_h     = np.full(n_steps, np.nan)
    p95_h      = np.full(n_steps, np.nan)
    max_h      = np.full(n_steps, np.nan)

    for i in range(n_steps):
        mask = flooded[i] & active
        extent_km2[i] = (flooded[i] * area).sum() / 1e6
        volume_Mm3[i] = (depth[i]   * area).sum() / 1e6
        if mask.any():
            vals      = depth[i][mask]
            mean_h[i] = vals.mean()
            p95_h[i]  = np.percentile(vals, 95)
            max_h[i]  = vals.max()

    return times, extent_km2, volume_Mm3, mean_h, p95_h, max_h


# ── MAIN ──────────────────────────────────────────────────────────────────────
OUT.mkdir(parents=True, exist_ok=True)

results = {}
for label, (scenario, _) in SCENARIOS.items():
    print(f"Loading {label} ...")
    results[label] = load_scenario(BASE / scenario)
    t, ext, vol, mh, p95, mxh = results[label]
    t_lbl = "hourly" if USE_HOURLY else "daily-max"
    print(f"  {len(t)} {t_lbl} steps  "
          f"peak ext={np.nanmax(ext):.1f} km²  "
          f"peak vol={np.nanmax(vol):.3f} Mm³  "
          f"peak depth={np.nanmax(mxh):.2f} m")

fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

panel_labels = ["a", "b", "c"]
for label, (scenario, color) in SCENARIOS.items():
    times, ext, vol, mh, p95, mxh = results[label]
    kw  = dict(linewidth=1.8, label=label, color=color)
    kw2 = dict(linewidth=1.4, color=color)

    # (a) flood extent
    axes[0].plot(times, ext, **kw)
    axes[0].fill_between(times, ext, alpha=0.10, color=color)

    # (b) flood volume
    axes[1].plot(times, vol, **kw)
    axes[1].fill_between(times, vol, alpha=0.10, color=color)

    # (c) depth stats — solid=max, dashed=95th, dotted=mean
    axes[2].plot(times, mxh, linestyle="-",  **kw2, label=f"{label} max")
    axes[2].plot(times, p95, linestyle="--", **kw2, label=f"{label} 95th pct")
    axes[2].plot(times, mh,  linestyle=":",  **kw2, label=f"{label} mean")

ylabels = [
    "Flood extent (km²)",
    "Flood volume (Mm³)",
    "Water depth h = zs − zb  (m)",
]
for ax, ylabel, plbl in zip(axes, ylabels, panel_labels):
    ax.set_ylabel(ylabel, fontsize=11, fontweight="bold")
    ax.set_title(plbl, loc="left", fontsize=11, fontweight="bold")
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax.xaxis.set_minor_locator(mdates.WeekdayLocator(byweekday=0))
    ax.grid(axis="both", alpha=0.3, linestyle="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

axes[0].legend(fontsize=10, framealpha=0.7)
axes[2].legend(fontsize=8, framealpha=0.7, ncol=len(SCENARIOS))

plt.setp(axes[2].xaxis.get_majorticklabels(), rotation=25, ha="right")
axes[2].set_xlabel("Time", fontsize=11)

scenario_str = " | ".join(SCENARIOS.keys())
t_lbl = "hourly" if USE_HOURLY else "daily-max"
fig.suptitle(
    f"SFINCS Flood Time Series — {scenario_str}  [{t_lbl}]\n"
    f"h = max(zs − zb, 0)  |  threshold = {FLOOD_THRESHOLD} m",
    fontsize=12, fontweight="bold",
)
plt.tight_layout()

slug    = "_".join(s.lower().replace(" ", "_").replace("(", "").replace(")", "")
                   for s in SCENARIOS)
out_path = OUT / f"timeseries_flood_{slug}_{t_lbl}.png"
plt.savefig(out_path, dpi=200, bbox_inches="tight")
print(f"\nSaved: {out_path}")
plt.show()

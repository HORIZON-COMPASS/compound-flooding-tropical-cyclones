"""
SFINCS Flood Time Series — 3-panel plot.

Computes at every output timestep:
  (a) Flood extent   [km²] — area where depth > FLOOD_THRESHOLD
  (b) Flood volume   [Mm³] — sum of depth × cell_area over flooded cells
  (c) Depth stats    [m]   — mean, 95th percentile, and max over flooded cells

Two sources are available, selected with SOURCE:

  "downscaled"  (default, recommended for reported numbers)
      Reads the per-timestep floodmaps that sfincs_postprocess.py already
      writes: plot_output/sfincs_output_hmax_period_*.tif
      These are downscaled onto the subgrid DEM (dep_subgrid.tif) and masked
      for permanent water (GSWO > 5%), so extents are directly consistent
      with sfincs_output_hmax_AllTime.tif.

  "coarse"
      Reads sfincs_map.nc and computes h = max(zs − zb, 0) on the SFINCS
      computational grid.  Cells are all-or-nothing at that resolution and
      zb is an effective bed level biased toward the low-lying part of each
      cell, so extents come out systematically HIGHER than the downscaled
      floodmaps.  Useful for temporal diagnostics, not for absolute extent.

  "both"
      Overlays the two on panels (a) and (b) — solid = downscaled,
      dashed = coarse — to make the resolution effect visible.
      Panel (c) shows the downscaled stats only, to stay readable.

Note that the peak of the extent time series is always ≤ the extent of
sfincs_output_hmax_AllTime.tif: the tif is an envelope (a cell counts if it
was ever wet), whereas the time series is the largest single-moment extent.
SHOW_ALLTIME_REF draws that envelope as a reference line on panel (a).

Volume and extent both use FLOOD_THRESHOLD, so the two panels describe the
same set of cells.

Run:
    pixi run -e compass-v1 python sfincs/plot_sfincs_timeseries.py
"""

import glob
import re
import warnings
from datetime import datetime
from pathlib import Path

import cftime
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import rasterio
import xarray as xr
import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs
from shapely import contains_xy

warnings.filterwarnings("ignore")

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
BASE  = Path("/home/azureuser/cloudfiles/code/Users/eloise.matthews/data/output_jul26")
OUT   = Path("/home/azureuser/cloudfiles/code/Users/eloise.matthews/hydromt_v1_somerset/analysis/plots")
SHP_PATH = Path("/home/azureuser/cloudfiles/code/Users/eloise.matthews/local_data/static/somerset_shape/new_shape2.shp")

# Scenarios to plot.  Each entry: (display_label, hex_colour).
# Set to a single entry for a single-scenario plot.
SCENARIOS = {
    "Counterfactual":   ("t1_downscaled","#1f77b4"),
    "Factual":     ("factual_downscaled","#ff7f0e"),
    "Future": ("t3_downscaled","#2ca02c"),
}

SOURCE            = "downscaled"  # "downscaled" | "coarse" | "both"
FLOOD_THRESHOLD   = 0.3         # m — minimum depth to count as flooded
SHOW_ALLTIME_REF  = True          # draw AllTime envelope extent on panel (a)
USE_HOURLY        = False         # coarse only: True → hourly zs, False → daily zsmax


# ── PATH HELPERS ──────────────────────────────────────────────────────────────
def _find_nc(scenario_dir: Path) -> str:
    matches = glob.glob(str(scenario_dir / "sfincs_map.nc"))
    if not matches:
        raise FileNotFoundError(
            f"No sfincs_map.nc found under:\n  {scenario_dir}"
        )
    return matches[0]


def _find_period_tifs(scenario_dir: Path) -> list:
    """Per-timestep downscaled floodmaps, sorted by the end of their period."""
    pattern = str(
        scenario_dir / "sfincs_output_hmax_period_*.tif"
    )
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(
            f"No per-period floodmaps found under:\n  {pattern}\n"
            "Run sfincs_postprocess.py for this scenario first, or set "
            'SOURCE = "coarse".'
        )
    stamped = [(_period_end(f), f) for f in matches]
    stamped.sort()
    return stamped


def _find_alltime_tif(scenario_dir: Path):
    matches = glob.glob(
        str(scenario_dir / "sfincs_output_hmax_AllTime.tif")
    )
    return matches[0] if matches else None


def _period_end(fn: str) -> datetime:
    """Parse the period END timestamp out of a floodmap filename.

    sfincs_postprocess.py names these `period_{tstart}_to_{tend}`, and the
    end is what lines up with the `timemax` coordinate in sfincs_map.nc.
    """
    m = re.search(r"period_\d{8}_\d{4}_to_(\d{8})_(\d{4})\.tif$", fn)
    if not m:
        raise ValueError(f"Cannot parse a period end timestamp from: {fn}")
    return datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H%M")


# ── GEOMETRY HELPERS ──────────────────────────────────────────────────────────
def _decode_times(ds: xr.Dataset) -> list:
    if USE_HOURLY:
        raw, units = ds["time"].values,    ds["time"].attrs["units"]
        calendar   = ds["time"].attrs.get("calendar", "standard")
    else:
        raw, units = ds["timemax"].values, ds["timemax"].attrs["units"]
        calendar   = ds["timemax"].attrs.get("calendar", "standard")
    return cftime.num2pydate(raw, units=units, calendar=calendar)


def _cell_areas(ds: xr.Dataset) -> np.ndarray:
    """Coarse-grid cell areas [m²] from the corner coordinates."""
    cx = ds["corner_x"].values   # (n+1, m+1)
    cy = ds["corner_y"].values
    v1x = cx[:-1, 1:]  - cx[:-1, :-1]
    v1y = cy[:-1, 1:]  - cy[:-1, :-1]
    v2x = cx[1:,  :-1] - cx[:-1, :-1]
    v2y = cy[1:,  :-1] - cy[:-1, :-1]
    return np.abs(v1x * v2y - v1y * v2x)   # (n, m)  m²


def _pixel_areas(src) -> np.ndarray:
    """Raster pixel areas [m²].

    Scalar for a projected CRS; a per-row column vector (broadcasts over the
    raster) for a geographic one, where pixel area shrinks with latitude.
    """
    res_x, res_y = abs(src.res[0]), abs(src.res[1])
    if src.crs is not None and src.crs.is_projected:
        return np.float64(res_x * res_y)

    rows = np.arange(src.height) + 0.5
    _, lats = src.transform * (np.zeros_like(rows), rows)
    metres_per_deg_lat = 110_540.0
    metres_per_deg_lon = 111_320.0 * np.cos(np.deg2rad(np.asarray(lats)))
    return (res_x * metres_per_deg_lon * res_y * metres_per_deg_lat)[:, np.newaxis]


def _stats_from_depth(h: np.ndarray, area) -> tuple:
    """Extent [km²], volume [Mm³] and depth stats over cells above threshold."""
    flooded = h > FLOOD_THRESHOLD
    extent  = float((flooded * area).sum() / 1e6)
    volume  = float((h * flooded * area).sum() / 1e6)
    if flooded.any():
        vals = h[flooded]
        return extent, volume, vals.mean(), np.percentile(vals, 95), vals.max()
    return extent, volume, np.nan, np.nan, np.nan


# ── LOADERS ───────────────────────────────────────────────────────────────────
def load_scenario_downscaled(scenario_dir: Path, flood_mask_shp: list) -> tuple:
    """Extent/volume/depth stats from the per-period downscaled floodmaps."""
    stamped = _find_period_tifs(scenario_dir)
    print(f"    {len(stamped)} downscaled floodmaps")

    times = [t for t, _ in stamped]
    n     = len(stamped)
    extent_km2 = np.empty(n)
    volume_Mm3 = np.empty(n)
    mean_h     = np.full(n, np.nan)
    p95_h      = np.full(n, np.nan)
    max_h      = np.full(n, np.nan)
    # crs = None
    # xs, ys = None
    area_mask = None

    area = None
    for i, (_, fn) in enumerate(stamped):
        with rasterio.open(fn) as src:
            if area is None:
                area = _pixel_areas(src)
            arr = src.read(1, masked=True)

            if area_mask is None:
                cols, rows = np.meshgrid(np.arange(src.width), np.arange(src.height))
                shape = cols.shape
                xs, ys = rasterio.transform.xy(src.transform, rows, cols)
        
        if area_mask is None:
            x, y = xs.reshape(shape), ys.reshape(shape)
            lon_lat_area = ccrs.TransverseMercator(central_longitude=-2.0, central_latitude=49.0, false_easting=400000.0, false_northing=-100000.0, scale_factor=0.9996012717).transform_points("EPSG:32630", x, y)
            masks_area = [contains_xy(flood_mask_shp[i], lon_lat_area[:, :, 0], lon_lat_area[:, :, 1]) for i in range(len(flood_mask_shp))]
            area_mask = np.logical_or.reduce(masks_area)

        arr = np.where(~area_mask, 0.0, arr)
        # Dry / permanently-masked pixels carry no water.
        h = np.ma.filled(arr.astype("float64"), 0.0)
        h[~np.isfinite(h)] = 0.0
        
        extent_km2[i], volume_Mm3[i], mean_h[i], p95_h[i], max_h[i] = _stats_from_depth(
            h, area
        )

    return times, extent_km2, volume_Mm3, mean_h, p95_h, max_h


def load_scenario_coarse(scenario_dir: Path, flood_mask_shp: list) -> tuple:
    """Extent/volume/depth stats from sfincs_map.nc on the computational grid."""
    nc = _find_nc(scenario_dir)
    print(f"    {nc}")
    ds = xr.open_dataset(nc, decode_times=False)

    # reproject shape file
    x, y = np.meshgrid(ds.x[0], ds.y[:,0])
    lon_lat = ccrs.TransverseMercator(central_longitude=-2.0, central_latitude=49.0, false_easting=400000.0, false_northing=-100000.0, scale_factor=0.9996012717).transform_points(ccrs.TransverseMercator(central_longitude=-3.0, false_easting=500000.0, scale_factor=0.9996), x, y)
    print(f"lon_lat shape: {lon_lat.shape}")

    masks_area = [contains_xy(flood_mask_shp[i], lon_lat[:, :, 0], lon_lat[:, :, 1]) for i in range(len(flood_mask_shp))]
    area_mask = np.logical_or.reduce(masks_area)
    # Save the mask as a dataset
    mask_ds = xr.DataArray(area_mask, name='flood_area', coords={'lat': lon_lat[:, 0, 1], 'lon': lon_lat[0, :, 0]})

    msk    = ds["msk"].values          # (n, m)
    zb     = ds["zb"].where(mask_ds.values).values           # (n, m)
    area   = _cell_areas(ds)           # (n, m)  m²
    active = msk > 0
    times  = _decode_times(ds)

    zs_data = (ds["zs"] if USE_HOURLY else ds["zsmax"]).where(mask_ds.values).values   # (T, n, m)
    ds.close()

    zs_clean = np.where(np.isnan(zs_data), zb[np.newaxis], zs_data)
    depth    = np.maximum(zs_clean - zb[np.newaxis], 0.0)
    depth[:, ~active] = 0.0

    n_steps    = depth.shape[0]
    extent_km2 = np.empty(n_steps)
    volume_Mm3 = np.empty(n_steps)
    mean_h     = np.full(n_steps, np.nan)
    p95_h      = np.full(n_steps, np.nan)
    max_h      = np.full(n_steps, np.nan)

    for i in range(n_steps):
        extent_km2[i], volume_Mm3[i], mean_h[i], p95_h[i], max_h[i] = _stats_from_depth(
            depth[i], area
        )

    return times, extent_km2, volume_Mm3, mean_h, p95_h, max_h


def load_alltime_extent(scenario_dir: Path, flood_mask_shp: list):
    """Extent [km²] of the AllTime envelope floodmap, or None if absent."""
    fn = _find_alltime_tif(scenario_dir)
    if fn is None:
        return None
    with rasterio.open(fn) as src:
        area = _pixel_areas(src)
        arr  = src.read(1, masked=True)
        cols, rows = np.meshgrid(np.arange(src.width), np.arange(src.height))
        shape = cols.shape
        xs, ys = rasterio.transform.xy(src.transform, rows, cols)

    x, y = xs.reshape(shape), ys.reshape(shape)
    lon_lat_area = ccrs.TransverseMercator(central_longitude=-2.0, central_latitude=49.0, false_easting=400000.0, false_northing=-100000.0, scale_factor=0.9996012717).transform_points("EPSG:32630", x, y)
    masks_area = [contains_xy(flood_mask_shp[i], lon_lat_area[:, :, 0], lon_lat_area[:, :, 1]) for i in range(len(flood_mask_shp))]
    area_mask = np.logical_or.reduce(masks_area)
    arr = np.where(~area_mask, 0.0, arr)
    h = np.ma.filled(arr.astype("float64"), 0.0)
    h[~np.isfinite(h)] = 0.0
    return float(((h > FLOOD_THRESHOLD) * area).sum() / 1e6)


LOADERS = {"downscaled": load_scenario_downscaled, "coarse": load_scenario_coarse}


# ── MAIN ──────────────────────────────────────────────────────────────────────
if SOURCE not in ("downscaled", "coarse", "both"):
    raise ValueError(f'SOURCE must be "downscaled", "coarse" or "both" — got {SOURCE!r}')

sources = ["downscaled", "coarse"] if SOURCE == "both" else [SOURCE]
primary = sources[0]

OUT.mkdir(parents=True, exist_ok=True)

# load in shape file for cropping
area_reader = shpreader.Reader(SHP_PATH)
flood_mask_shp = [record.geometry for record in area_reader.records()]

results   = {}
alltime   = {}
for label, (scenario, _) in SCENARIOS.items():
    print(f"Loading {label} ...")
    for src_name in sources:
        results[(label, src_name)] = LOADERS[src_name](BASE / scenario, flood_mask_shp)
        t, ext, vol, mh, p95, mxh = results[(label, src_name)]
        if src_name == "coarse":
            step_lbl = "hourly" if USE_HOURLY else "daily-max"
        else:
            step_lbl = "daily-max, subgrid"
        print(f"  [{src_name:10s}] {len(t)} {step_lbl} steps  "
              f"peak ext={np.nanmax(ext):.1f} km²  "
              f"peak vol={np.nanmax(vol):.3f} Mm³  "
              f"peak depth={np.nanmax(mxh):.2f} m")

    if SHOW_ALLTIME_REF:
        alltime[label] = load_alltime_extent(BASE / scenario, flood_mask_shp)
        if alltime[label] is not None:
            print(f"  [AllTime   ] envelope extent={alltime[label]:.1f} km²")

fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

panel_labels = ["a", "b"]
linestyles   = {"downscaled": "-", "coarse": "--"}
# Shading only helps when there is a single curve to read against the axis;
# with several scenarios the overlapping fills obscure each other.
fill = len(SCENARIOS) == 1

for label, (scenario, color) in SCENARIOS.items():
    for src_name in sources:
        times, ext, vol, mh, p95, mxh = results[(label, src_name)]
        leg = label if len(sources) == 1 else f"{label} ({src_name})"
        ls  = linestyles[src_name]

        # (a) flood extent
        axes[0].plot(times, ext, linewidth=1.8, label=leg, color=color, linestyle=ls)
        if fill and src_name == primary:
            axes[0].fill_between(times, ext, alpha=0.10, color=color)

        # (b) flood volume
        axes[1].plot(times, vol, linewidth=1.8, label=leg, color=color, linestyle=ls)
        if fill and src_name == primary:
            axes[1].fill_between(times, vol, alpha=0.10, color=color)

    # (c) depth stats for the primary source only — solid=max, dashed=95th, dotted=mean
    # times, _, _, mh, p95, mxh = results[(label, primary)]
    # kw = dict(linewidth=1.4, color=color)
    # axes[2].plot(times, mxh, linestyle="-",  **kw, label=f"{label} max")
    # axes[2].plot(times, p95, linestyle="--", **kw, label=f"{label} 95th pct")
    # axes[2].plot(times, mh,  linestyle=":",  **kw, label=f"{label} mean")

    # AllTime envelope reference on panel (a)
    if SHOW_ALLTIME_REF and alltime.get(label) is not None:
        axes[0].axhline(
            alltime[label], color=color, linestyle=":", linewidth=1.2, alpha=0.7,
            label="AllTime envelope" if label == list(SCENARIOS)[0] else None,
        )

ylabels = [
    "Flood extent (km²)",
    "Flood volume (Mm³)",
    # "Water depth (m)",
]
for ax, ylabel, plbl in zip(axes, ylabels, panel_labels):
    ax.set_ylabel(ylabel, fontsize=11, fontweight="bold")
    ax.set_title(plbl, loc="left", fontsize=11, fontweight="bold")
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    ax.xaxis.set_minor_locator(mdates.WeekdayLocator(byweekday=0))
    ax.grid(axis="both", alpha=0.3, linestyle="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

axes[0].legend(fontsize=9, framealpha=0.7, title="Scenario", title_fontsize=9, loc="upper left", bbox_to_anchor=(0.01, 0.90))
axes[0].set_ylim(top=140)

plt.setp(axes[1].xaxis.get_majorticklabels(), rotation=25, ha="right")
axes[1].set_xlabel("Time", fontsize=11)

# if primary == "downscaled":
#     src_note = "subgrid-downscaled floodmaps, GSWO-masked (h from sfincs_output_hmax_period_*.tif)"
# else:
#     src_note = "computational grid (h = max(zs − zb, 0) from sfincs_map.nc)"
# if SOURCE == "both":
#     src_note += "   |   dashed = computational grid"

scenario_str = " | ".join(SCENARIOS.keys())
fig.suptitle(
    f"SFINCS Flood Time Series - Somerset Winter 2013/14", # {scenario_str}\n"
    # f"{src_note}  |  threshold = {FLOOD_THRESHOLD} m",
    fontsize=12, fontweight="bold",
)
plt.tight_layout()

slug    = "_".join(s.lower().replace(" ", "_").replace("(", "").replace(")", "")
                   for s in SCENARIOS)
suffix  = SOURCE if SOURCE != "coarse" else ("hourly" if USE_HOURLY else "daily-max")
out_path = OUT / f"timeseries_flood_{slug}_{suffix}.png"
plt.savefig(out_path, dpi=400, bbox_inches="tight")
print(f"\nSaved: {out_path}")
plt.show()

"""
Report figures — Durban 2022 Factual scenario.

Produces three publication-ready single-panel maps:
  1. Flood hazard map   (SFINCS max water depth)
  2. Economic damage map (Delft-FIAT, rasterized total damage per building)
  3. Population exposure map (Delft-FIAT + WorldPop, rasterized population)

Also prints a concise table of key result numbers.

Output directory: /p/11210471-001-compass/04_Results/report_figures/

by @dumontgoulart
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import xarray as xr
import rioxarray as rxr
from matplotlib.colors import BoundaryNorm
import contextily as ctx
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")


# ===== CONFIGURATION =====
EVENT_DISPLAY = "Durban 2022"
REGION        = "durban"
RUNNAME       = "Durban2022"
SCENARIO      = "event_precip_era5_hourly_CF0_no_wind"

BASE_RUNS = Path("/p/11210471-001-compass/03_Runs")
OUTPUT_DIR = Path("/p/11210471-001-compass/04_Results/report_figures")

SFINCS_DIR = BASE_RUNS / REGION / RUNNAME / "sfincs" / SCENARIO
FIAT_DIR   = BASE_RUNS / REGION / RUNNAME / "fiat"  / SCENARIO

FLOOD_MAP_PATH  = SFINCS_DIR / "plot_output" / "floodmap.tif"
FIAT_OUTPUT     = FIAT_DIR   / "output" / "spatial_with_pop_and_flood.fgb"
FIAT_FALLBACK   = FIAT_DIR   / "output" / "spatial.fgb"

FLOOD_THRESHOLD = 0.05   # m — minimum depth counted as flooded
HMAX_VMAX       = 2.0    # m — cap for colorbar (deeper values still shown at max colour)
GRID_RESOLUTION = 0.005  # degrees (~500 m) — raster grid for damage/pop maps
BASEMAP_ZOOM    = 12


# ===== HELPERS =====
def _load_fiat() -> gpd.GeoDataFrame:
    path = FIAT_OUTPUT if FIAT_OUTPUT.exists() else FIAT_FALLBACK
    gdf = gpd.read_file(path)
    if str(gdf.crs) != "EPSG:4326":
        gdf = gdf.to_crs("EPSG:4326")
    return gdf


def _rasterize(gdf: gpd.GeoDataFrame, value_col: str, bounds: tuple) -> xr.DataArray:
    """Sum value_col per GRID_RESOLUTION cell; NaN where no buildings."""
    sub = gdf[gdf[value_col] > 0].copy()
    if len(sub) == 0:
        return None
    c  = sub.geometry.centroid
    xs, ys, vals = c.x.values, c.y.values, sub[value_col].values

    minx, miny, maxx, maxy = bounds
    r = GRID_RESOLUTION
    x_edges = np.arange(minx - r/2, maxx + r, r)
    y_edges = np.arange(miny - r/2, maxy + r, r)
    xc = (x_edges[:-1] + x_edges[1:]) / 2
    yc = (y_edges[:-1] + y_edges[1:]) / 2

    xb = np.clip(np.digitize(xs, x_edges) - 1, 0, len(xc) - 1)
    yb = np.clip(np.digitize(ys, y_edges) - 1, 0, len(yc) - 1)

    grid  = np.zeros((len(yc), len(xc)))
    count = np.zeros_like(grid)
    for i in range(len(vals)):
        if not np.isnan(vals[i]):
            grid[yb[i], xb[i]]  += vals[i]
            count[yb[i], xb[i]] += 1
    grid = np.where(count > 0, grid, np.nan)
    return xr.DataArray(data=grid, dims=["y", "x"],
                        coords={"y": yc, "x": xc})


def _basemap(ax):
    try:
        ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik,
                        crs="EPSG:4326", attribution=False, zorder=1,
                        zoom=BASEMAP_ZOOM)
    except Exception:
        pass


def _style(ax, title, fontsize=14):
    ax.set_title(title, fontsize=fontsize, fontweight="bold", pad=10)
    ax.set_xlabel("Longitude [°]", fontsize=fontsize - 3)
    ax.set_ylabel("Latitude [°]", fontsize=fontsize - 3)
    ax.set_aspect("equal")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ===== KEY STATISTICS =====
def print_key_numbers(da_flood: xr.DataArray, gdf: gpd.GeoDataFrame):
    affected = gdf[gdf["total_damage"] > 0]
    has_pop  = "population" in gdf.columns

    # Flood extent from raster
    flooded = da_flood.where(da_flood > FLOOD_THRESHOLD)
    # Approximate cell area at ~30°S
    dx_m = float(np.abs(da_flood.x.diff("x").median())) * 96_000
    dy_m = float(np.abs(da_flood.y.diff("y").median())) * 111_000
    cell_m2 = dx_m * dy_m
    flood_extent_km2 = float((da_flood > FLOOD_THRESHOLD).sum() * cell_m2) / 1e6
    flood_vol_Mm3    = float(flooded.sum(skipna=True) * cell_m2) / 1e6

    total_dmg = affected["total_damage"].sum()
    mean_dmg  = affected["total_damage"].mean()
    max_dmg   = affected["total_damage"].max()
    pop_exp   = affected["population"].sum() if has_pop else float("nan")

    w = 60
    print("\n" + "=" * w)
    print(f"  KEY RESULTS — {EVENT_DISPLAY} | Factual scenario")
    print("=" * w)
    print(f"  {'Flood hazard (SFINCS)'}")
    print(f"    Flood extent (>{FLOOD_THRESHOLD} m):  {flood_extent_km2:.1f} km²")
    print(f"    Flood volume:                {flood_vol_Mm3:.2f} Mm³")
    print(f"    Max water depth:             {float(da_flood.max(skipna=True)):.2f} m")
    print(f"    Mean depth (flooded cells):  {float(flooded.mean(skipna=True)):.2f} m")
    print()
    print(f"  {'Flood impact (Delft-FIAT)'}")
    print(f"    Affected buildings:          {len(affected):,}")
    print(f"    Total economic damage:       USD {total_dmg:,.0f}")
    print(f"                                 (~USD {total_dmg/1e6:.0f}M)")
    print(f"    Mean damage / building:      USD {mean_dmg:,.0f}")
    print(f"    Max damage (single bldg):    USD {max_dmg:,.0f}")
    if has_pop:
        print(f"    Population exposed:          {pop_exp:,.0f}")
    print("=" * w + "\n")

    return {
        "flood_extent_km2": round(flood_extent_km2, 2),
        "flood_volume_Mm3": round(flood_vol_Mm3, 3),
        "max_depth_m": round(float(da_flood.max(skipna=True)), 2),
        "mean_depth_flooded_m": round(float(flooded.mean(skipna=True)), 2),
        "n_affected_buildings": int(len(affected)),
        "total_damage_USD": round(total_dmg, 0),
        "mean_damage_per_building_USD": round(mean_dmg, 0),
        "max_damage_single_building_USD": round(max_dmg, 0),
        "population_exposed": round(pop_exp, 0) if has_pop else None,
    }


# ===== PANEL RENDERERS (usable standalone or inside a combined figure) =====
def _render_flood_panel(da: xr.DataArray, ax, fontsize: int = 11):
    """Draw the flood hazard panel onto ax; returns the mappable."""
    norm = BoundaryNorm(np.linspace(0, HMAX_VMAX, 21), ncolors=256, clip=True)
    cmap = plt.get_cmap("Blues")
    im = da.where(da > FLOOD_THRESHOLD).plot(
        ax=ax, cmap=cmap, norm=norm, add_colorbar=False,
        x="x", y="y", alpha=0.85, zorder=2,
    )
    _basemap(ax)
    _style(ax, "Maximum Flood Depth", fontsize=fontsize)
    return im, f"Max Water Depth [m]  (capped at {HMAX_VMAX} m)", "max"


def _render_damage_panel(da: xr.DataArray, ax, fontsize: int = 11):
    """Draw the economic damage panel onto ax; returns the mappable."""
    vmax = float(np.nanquantile(da.values, 0.95))
    norm = BoundaryNorm(np.linspace(0, vmax, 11), ncolors=256, clip=True)
    cmap = plt.get_cmap("Reds")
    im = da.plot(ax=ax, cmap=cmap, norm=norm, add_colorbar=False,
                 x="x", y="y", alpha=0.85, zorder=2)
    _basemap(ax)
    _style(ax, "Economic Flood Damage", fontsize=fontsize)
    return im, f"Total Damage per {GRID_RESOLUTION}° cell [USD]", "max"


def _render_pop_panel(da: xr.DataArray, ax, fontsize: int = 11):
    """Draw the population exposure panel onto ax; returns the mappable."""
    vmax = float(np.nanquantile(da.values, 0.95))
    norm = BoundaryNorm(np.linspace(0, vmax, 11), ncolors=256, clip=True)
    cmap = plt.get_cmap("YlOrRd")
    im = da.plot(ax=ax, cmap=cmap, norm=norm, add_colorbar=False,
                 x="x", y="y", alpha=0.85, zorder=2)
    _basemap(ax)
    _style(ax, "Population Exposed to Flooding", fontsize=fontsize)
    return im, f"Population per {GRID_RESOLUTION}° cell", "max"


def _add_cbar(fig, im, ax, label, extend, fmt=None):
    cbar = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02, extend=extend)
    cbar.set_label(label, fontsize=10)
    if fmt:
        cbar.ax.yaxis.set_major_formatter(fmt)
    cbar.ax.tick_params(labelsize=9)


# ===== FIGURE 1 — FLOOD HAZARD MAP =====
def plot_flood_map(da: xr.DataArray, output_path: Path):
    fig, ax = plt.subplots(figsize=(10, 9))
    im, cbar_label, extend = _render_flood_panel(da, ax)
    _add_cbar(fig, im, ax, cbar_label, extend)
    ax.set_title(f"Maximum Flood Depth — {EVENT_DISPLAY}", fontsize=14,
                 fontweight="bold", pad=10)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# ===== FIGURE 2 — ECONOMIC DAMAGE MAP =====
def plot_damage_map(gdf: gpd.GeoDataFrame, output_path: Path):
    da = _rasterize(gdf, "total_damage",
                    gdf[gdf["total_damage"] > 0].geometry.centroid.total_bounds)
    if da is None:
        print("  No damage data, skipping damage map.")
        return
    fig, ax = plt.subplots(figsize=(10, 9))
    im, cbar_label, extend = _render_damage_panel(da, ax)
    _add_cbar(fig, im, ax, cbar_label, extend,
              fmt=ticker.FuncFormatter(
                  lambda v, _: f"${v/1e3:.0f}K" if v < 1e6 else f"${v/1e6:.1f}M"))
    ax.set_title(f"Economic Flood Damage — {EVENT_DISPLAY}", fontsize=14,
                 fontweight="bold", pad=10)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# ===== FIGURE 3 — POPULATION EXPOSURE MAP =====
def plot_pop_map(gdf: gpd.GeoDataFrame, output_path: Path):
    if "population" not in gdf.columns:
        print("  No population column — run fiat_add_pop_exposed_metric.py first.")
        return
    da = _rasterize(gdf, "population",
                    gdf[gdf["population"] > 0].geometry.centroid.total_bounds)
    if da is None:
        print("  No population data, skipping pop map.")
        return
    fig, ax = plt.subplots(figsize=(10, 9))
    im, cbar_label, extend = _render_pop_panel(da, ax)
    _add_cbar(fig, im, ax, cbar_label, extend)
    ax.set_title(f"Population Exposed to Flooding — {EVENT_DISPLAY}", fontsize=14,
                 fontweight="bold", pad=10)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# ===== FIGURE 4 — COMBINED THREE-PANEL =====
def plot_combined(da_flood: xr.DataArray, gdf: gpd.GeoDataFrame, output_path: Path):
    """Single figure with flood | damage | population panels left to right.

    All panels share the same geographic bounding box (derived from FIAT affected
    buildings) so they render at the same size and aspect ratio.
    """
    has_pop = "population" in gdf.columns
    fs = 15  # font size for combined figure

    # Common bounds from FIAT affected buildings extent
    aff = gdf[gdf["total_damage"] > 0]
    c = aff.geometry.centroid
    pad = 0.02  # degrees of padding around the data
    common_bounds = (c.x.min() - pad, c.y.min() - pad,
                     c.x.max() + pad, c.y.max() + pad)
    xlim = (common_bounds[0], common_bounds[2])
    ylim = (common_bounds[1], common_bounds[3])

    # Rasterise damage and population to the same common bounds
    da_dmg = _rasterize(gdf, "total_damage", common_bounds)
    da_pop = _rasterize(gdf, "population",   common_bounds) if has_pop else None

    fig, axes = plt.subplots(1, 3, figsize=(24, 8))

    # Panel 1 — flood (render full raster; xlim/ylim below will zoom to common extent)
    im1, lbl1, ext1 = _render_flood_panel(da_flood, axes[0], fontsize=fs)
    _add_cbar(fig, im1, axes[0], lbl1, ext1)

    # Panel 2 — damage
    if da_dmg is not None:
        im2, lbl2, ext2 = _render_damage_panel(da_dmg, axes[1], fontsize=fs)
        _add_cbar(fig, im2, axes[1], lbl2, ext2,
                  fmt=ticker.FuncFormatter(
                      lambda v, _: f"${v/1e3:.0f}K" if v < 1e6 else f"${v/1e6:.1f}M"))
    else:
        axes[1].set_visible(False)

    # Panel 3 — population
    if da_pop is not None:
        im3, lbl3, ext3 = _render_pop_panel(da_pop, axes[2], fontsize=fs)
        _add_cbar(fig, im3, axes[2], lbl3, ext3)
    else:
        axes[2].set_visible(False)

    # Enforce identical axes limits on all visible panels
    for ax in axes:
        if ax.get_visible():
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

    plt.suptitle(f"Flood Hazard, Damage & Population Exposure — {EVENT_DISPLAY}",
                 fontsize=18, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# ===== MAIN =====
def main():
    print(f"\n{'='*60}")
    print(f"  Report figures — {EVENT_DISPLAY} | Factual")
    print(f"{'='*60}\n")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load flood map
    print("Loading SFINCS flood map...")
    if not FLOOD_MAP_PATH.exists():
        print(f"  ERROR: {FLOOD_MAP_PATH} not found.")
        return
    da_flood = rxr.open_rasterio(FLOOD_MAP_PATH).squeeze("band", drop=True)
    if da_flood.rio.crs and str(da_flood.rio.crs) != "EPSG:4326":
        da_flood = da_flood.rio.reproject("EPSG:4326")
    # Replace nodata with NaN
    nodata = da_flood.rio.nodata
    if nodata is not None:
        da_flood = da_flood.where(da_flood != nodata)

    # Load FIAT output
    print("Loading FIAT spatial output...")
    path = FIAT_OUTPUT if FIAT_OUTPUT.exists() else FIAT_FALLBACK
    if not path.exists():
        print(f"  ERROR: FIAT output not found at {path}.")
        return
    gdf = _load_fiat()
    print(f"  Loaded {len(gdf):,} buildings from {path.name}")

    # Key numbers
    stats = print_key_numbers(da_flood, gdf)

    # Export stats to CSV
    stats_path = OUTPUT_DIR / f"key_numbers_{EVENT_DISPLAY.lower().replace(' ', '_')}_factual.csv"
    pd.DataFrame([stats]).T.rename(columns={0: "Value"}).to_csv(stats_path)
    print(f"Saved: {stats_path}")

    # Figures
    print("\nGenerating report figures...")

    print("  Figure 1: Flood hazard map...")
    plot_flood_map(
        da_flood,
        OUTPUT_DIR / f"report_flood_map_{EVENT_DISPLAY.lower().replace(' ', '_')}_factual.png",
    )

    print("  Figure 2: Economic damage map...")
    plot_damage_map(
        gdf,
        OUTPUT_DIR / f"report_damage_map_{EVENT_DISPLAY.lower().replace(' ', '_')}_factual.png",
    )

    print("  Figure 3: Population exposure map...")
    plot_pop_map(
        gdf,
        OUTPUT_DIR / f"report_pop_map_{EVENT_DISPLAY.lower().replace(' ', '_')}_factual.png",
    )

    print("  Figure 4: Combined three-panel figure...")
    plot_combined(
        da_flood, gdf,
        OUTPUT_DIR / f"report_combined_{EVENT_DISPLAY.lower().replace(' ', '_')}_factual.png",
    )

    print(f"\n{'='*60}")
    print(f"  Done. Outputs saved to: {OUTPUT_DIR}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()

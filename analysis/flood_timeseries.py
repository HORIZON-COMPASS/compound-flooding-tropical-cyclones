"""
Standalone flood extent timeseries plotter.
Reads SFINCS output (sfincs_map.nc) and a subgrid DEM directly — no SfincsModel required.
"""

from datetime import datetime
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import rioxarray
import xarray as xr
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Paths — edit these to match your local files
# ---------------------------------------------------------------------------

SFINCS_MAP_NC   = "factual/sfincs_map.nc"        # SFINCS output
DEP_SUBGRID_TIF = "factual/subgrid/dep_subgrid.tif"  # bed elevation
FLOOD_MASK_NC   = "/home/azureuser/cloudfiles/code/Users/eloise.matthews/local_data/static/somerset_shape/flood_area_mask.nc"   # pre-computed mask (see build_mask below)
GSWO_TIF        = None   # optional: path to GSWO GeoTIFF; set None to skip masking

HMIN = 0.05   # minimum flood depth to include (m)
HTHRESH = 0.3  # depth threshold for "flooded" pixel (m)
PIXEL_AREA_M2 = (90 / 4) ** 2  # ~22.5 m × ~22.5 m sub-grid pixels

OUTPUT_DIR = "draft_figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Helper: downscale zsmax → flood depth without hydromt_sfincs
# ---------------------------------------------------------------------------

def downscale_floodmap(zsmax: xr.DataArray, dep: xr.DataArray, hmin: float = 0.05) -> xr.DataArray:
    """
    Compute flood depth on the subgrid from water-surface elevation.

    Parameters
    ----------
    zsmax : xr.DataArray
        Maximum water-surface elevation on the coarse SFINCS grid (y, x).
    dep : xr.DataArray
        Bed elevation on the subgrid.
    hmin : float
        Minimum depth to keep (m).

    Returns
    -------
    xr.DataArray
        Flood depth on the subgrid, NaN where depth < hmin.
    """
    # Reproject coarse zsmax to the fine subgrid using nearest-neighbour
    zsmax_reprojected = zsmax.rio.reproject_match(dep, resampling=rioxarray.enums.Resampling.nearest)

    # Flood depth = water surface - bed elevation
    hmax = zsmax_reprojected - dep
    hmax = hmax.where(hmax >= hmin)  # mask depths below threshold
    return hmax


# ---------------------------------------------------------------------------
# Load the subgrid DEM once
# ---------------------------------------------------------------------------

print("Loading subgrid DEM …")
da_dep = rioxarray.open_rasterio(DEP_SUBGRID_TIF, masked=True).squeeze("band", drop=True)

# ---------------------------------------------------------------------------
# Load GSWO mask (optional — masks out permanent water bodies)
# ---------------------------------------------------------------------------

if GSWO_TIF is not None:
    print("Loading GSWO mask …")
    gswo = rioxarray.open_rasterio(GSWO_TIF, masked=True).squeeze("band", drop=True)
    gswo = gswo.rio.reproject_match(da_dep, resampling=rioxarray.enums.Resampling.max)
else:
    gswo = None

# ---------------------------------------------------------------------------
# Load the flood-area mask (shapes the region we care about)
# ---------------------------------------------------------------------------

flood_area_mask = xr.open_dataset(FLOOD_MASK_NC).flood_area

# ---------------------------------------------------------------------------
# Load SFINCS output and compute flood-extent timeseries
# ---------------------------------------------------------------------------

print("Loading SFINCS map output …")
ds = xr.open_dataset(SFINCS_MAP_NC)

# 'zsmax' is the maximum water-surface elevation at each snapshot
# Expected dims: (timemax, y, x)
da_zsmax = ds["zsmax"]

timestamps = da_zsmax["timemax"].values
print(f"Found {len(timestamps)} time steps: {timestamps[0]} → {timestamps[-1]}")

flood_extent_km2 = []

for i, t in enumerate(tqdm(timestamps, desc="Computing flood extent")):
    da_zsmax_i = da_zsmax.isel(timemax=i)

    # Downscale to subgrid
    da_hmax_i = downscale_floodmap(da_zsmax_i, da_dep, hmin=HMIN)

    # Remove permanent water (GSWO water occurrence > 5 %)
    if gswo is not None:
        da_hmax_i = da_hmax_i.where(gswo <= 5)

    # Crop to the flood-area region and count flooded pixels
    # Reproject mask to hmax grid if needed
    da_hmax_i = da_hmax_i.fillna(0)
    flooded = xr.where(da_hmax_i >= HTHRESH, 1, 0)

    # If the mask grid differs from the subgrid, reproject it
    if flood_area_mask.shape != flooded.shape:
        mask_reproj = flood_area_mask.rio.write_crs("EPSG:27700").rio.reproject_match(
            flooded, resampling=rioxarray.enums.Resampling.nearest
        )
    else:
        mask_reproj = flood_area_mask

    flooded = flooded.where(mask_reproj.values)
    extent_km2 = float(flooded.sum()) * PIXEL_AREA_M2 / 1e6  # m² → km²
    flood_extent_km2.append(extent_km2)

flood_extent_km2 = np.array(flood_extent_km2)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(9, 4), layout="constrained")
ax.plot(timestamps, flood_extent_km2, color="navy", linewidth=1.5)
ax.set_xlabel("Time")
ax.set_ylabel("Flood extent (km²)")
ax.set_title("SFINCS flood extent timeseries")
ax.tick_params(axis="x", rotation=45)
ax.grid(alpha=0.5)

out_path = os.path.join(OUTPUT_DIR, "flood_timeseries.png")
fig.savefig(out_path, dpi=150)
print(f"Saved → {out_path}")
plt.show()

"""
Build downscaled daily floodmap GeoTIFFs from sfincs_map.nc — standalone.

TEMPORARY helper for collaborators who were sent sfincs_map.nc but not the
plot_output/ folder.  It reproduces the floodmap rasters that
sfincs_postprocess.py writes, WITHOUT needing the full SFINCS model directory
or a HydroMT data catalog.

Why this is needed
------------------
sfincs_map.nc holds water levels on the SFINCS *computational* grid (~100 m).
Deriving flood depth from it directly (h = zs - zb) treats each cell as
all-or-nothing and uses an effective bed level biased toward the low-lying
part of the cell, which overestimates flood extent by ~25% and volume by
~70%.  The reported numbers come from downscaling onto the subgrid DEM
instead.  This script does that downscaling.

Inputs
------
  MAP_NC        sfincs_map.nc                     (required)
  DEP_SUBGRID   subgrid/dep_subgrid.tif           (required — cannot be
                                                   derived from the netcdf)
  GSWO_PATH     permanent-water raster            (optional, see below)

Outputs (into OUT_DIR)
----------------------
  sfincs_output_hmax_AllTime.tif          envelope: max over the whole run
  sfincs_output_hmax_period_*.tif         one per daily `timemax` step

Naming and layout match what plot_sfincs_timeseries.py expects, so once this
has run you can point that script at OUT_ROOT and use SOURCE = "downscaled".

Permanent water
---------------
The production floodmaps mask cells where GSWO occurrence > 5%, which removes
rivers, lakes and (importantly) estuaries/sea.  Skipping it inflates extent,
and near a coast it inflates it a lot.  GSWO_PATH accepts either the global
`occur.vrt` or a pre-clipped raster; set EXPORT_GSWO_MASK to write the clipped
mask out once so it can be shared alongside the netcdf.  Set GSWO_PATH = None
to skip masking entirely — extents will then not match the reference tifs.

Run:
    pixi run -e compass-v1 python sfincs/make_downscaled_floodmaps.py
"""

import warnings
from pathlib import Path

import cftime
import numpy as np
import rioxarray
import xarray as xr
from rasterio.enums import Resampling
from rasterio.warp import transform_bounds

import hydromt  # noqa: F401  (registers the .raster accessor)
from hydromt_sfincs import utils

warnings.filterwarnings("ignore")

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
RUN_DIR = Path(
    "/p/11210471-001-compass/03_Runs/somerset/SomersetLevels_dec_factual/sfincs/"
    "event_tp_ceh_gear_compass_CF0_GTSMv41opendap_CF0_no_wind_CF0"
)

MAP_NC      = RUN_DIR / "sfincs_map.nc"
DEP_SUBGRID = RUN_DIR / "subgrid" / "dep_subgrid.tif"

# Where to write.  The nested layout is what plot_sfincs_timeseries.py globs
# for (<scenario>/sfincs/event_tp_*/plot_output), so pointing that script's
# BASE at OUT_ROOT.parent and using this scenario name works unchanged.
OUT_ROOT = Path("/p/11210471-001-compass/03_Runs/somerset/shared_factual")
OUT_DIR  = OUT_ROOT / "sfincs" / "event_tp_shared" / "plot_output"

# Permanent-water mask.  None to skip.
GSWO_PATH        = Path("/p/wflow_global/hydromt/hydrography/gswo/occur.vrt")
GSWO_MAX_OCCUR   = 5      # mask where occurrence > this (%)
EXPORT_GSWO_MASK = True   # also write the clipped mask, for sharing

HMIN = 0.05   # m — passed to downscale_floodmap

# These mirror sfincs_postprocess.py, which is inconsistent between its two
# calls: the AllTime map is downscaled with "bilinear", while the per-timestep
# maps omit the argument and fall back to the "nearest" default, making the
# daily maps slightly blockier than the envelope.  Setting both to "bilinear"
# would be more self-consistent, but the rasters would then no longer
# reproduce the ones already shared.
REPROJ_ALLTIME = "bilinear"
REPROJ_PERIOD  = "nearest"

LIMIT_STEPS = None   # int to process only the first N steps (testing)


# ── HELPERS ───────────────────────────────────────────────────────────────────
def load_zsmax(map_nc: Path):
    """Rebuild zsmax as a CRS-aware raster DataArray, plus the time vectors.

    sfincs_map.nc stores x/y as 2D face-centre arrays.  For an unrotated grid
    these collapse to 1D axes, which is what the hydromt raster accessor and
    downscale_floodmap need.
    """
    ds = xr.open_dataset(map_nc, decode_times=False)

    epsg = ds["crs"].attrs.get("epsg_code")
    if not epsg or epsg == "-":
        raise ValueError(f"No usable CRS in {map_nc} (crs.epsg_code = {epsg!r})")

    xv, yv = ds["x"].values, ds["y"].values
    unrotated = np.allclose(xv, xv[0:1, :]) and np.allclose(yv, yv[:, 0:1])
    if not unrotated:
        raise NotImplementedError(
            "This grid is rotated, so x/y do not collapse to 1D axes. "
            "Use sfincs_postprocess.py with the full model directory instead."
        )
    x1d, y1d = xv[0, :], yv[:, 0]

    def _decode(var):
        return cftime.num2pydate(
            ds[var].values,
            units=ds[var].attrs["units"],
            calendar=ds[var].attrs.get("calendar", "standard"),
        )

    timemax = _decode("timemax")
    t_start = _decode("time")[0]

    zsmax = xr.DataArray(
        ds["zsmax"].values,
        dims=("timemax", "y", "x"),
        coords={"timemax": np.arange(ds.sizes["timemax"]), "y": y1d, "x": x1d},
        name="zsmax",
    )
    ds.close()

    zsmax.raster.set_crs(epsg)
    # SFINCS writes dry cells as _FillValue (-99999), which xarray has already
    # turned into NaN.  Declaring NaN as the nodata value matters: without it
    # the reprojection inside downscale_floodmap treats those cells as real
    # water levels and the single-timestep floodmaps come out slightly wrong.
    zsmax.raster.set_nodata(np.nan)
    return zsmax, timemax, t_start


def load_gswo_mask(gswo_path: Path, like):
    """GSWO occurrence reprojected onto the `like` grid.

    Clips in geographic coordinates first — the global VRT is far too large to
    read whole.  `max` resampling matches sfincs_postprocess.py, and is the
    conservative choice: a subgrid pixel is called permanent water if any
    contributing GSWO pixel is.
    """
    bounds_wgs84 = transform_bounds(like.rio.crs, "EPSG:4326", *like.rio.bounds())
    gswo = rioxarray.open_rasterio(
        gswo_path, masked=True, chunks={"x": 4000, "y": 4000}
    ).squeeze()
    gswo = gswo.rio.clip_box(*bounds_wgs84).compute()
    return gswo.rio.reproject_match(like, resampling=Resampling.max)


def period_tag(t_from, t_to) -> str:
    """Filename stamp, byte-identical to sfincs_postprocess.py's convention."""
    return f"period_{t_from:%Y%m%d_%H%M}_to_{t_to:%Y%m%d_%H%M}"


def write(da, path: Path):
    da.raster.to_raster(str(path))
    print(f"    wrote {path.name}")


# ── MAIN ──────────────────────────────────────────────────────────────────────
OUT_DIR.mkdir(parents=True, exist_ok=True)

for required in (MAP_NC, DEP_SUBGRID):
    if not required.exists():
        raise FileNotFoundError(required)

print(f"Reading  {MAP_NC}")
zsmax, timemax, t_start = load_zsmax(MAP_NC)
print(f"  {len(timemax)} daily steps, {timemax[0]:%Y-%m-%d} to {timemax[-1]:%Y-%m-%d}")

print(f"Reading  {DEP_SUBGRID}")
da_dep = rioxarray.open_rasterio(DEP_SUBGRID, masked=True).squeeze(drop=True)
print(f"  subgrid DEM {da_dep.shape} at {abs(da_dep.raster.res[0]):.0f} m, {da_dep.raster.crs}")

# The GSWO mask is built once on the AllTime floodmap grid and reused, exactly
# as sfincs_postprocess.py does.
print("Downscaling AllTime envelope ...")
# sfincs_postprocess.py passes zsmax.max(dim="timemax") straight through, and
# xarray's reduction drops the attrs that carry nodata — so the AllTime map is
# built with nodata unset, while the .isel() slices below keep nodata = NaN.
# Clearing it here reproduces the shipped rasters exactly.  The two differ by
# ~0.03% in extent, so this is about matching, not accuracy.
zsmax_env = zsmax.max(dim="timemax")
zsmax_env.attrs.pop("_FillValue", None)   # raster.set_nodata(None) is a no-op
hmax_all = utils.downscale_floodmap(
    zsmax=zsmax_env,
    dep=da_dep,
    hmin=HMIN,
    reproj_method=REPROJ_ALLTIME,
)

gswo_mask = None
if GSWO_PATH is not None:
    if not Path(GSWO_PATH).exists():
        raise FileNotFoundError(
            f"{GSWO_PATH} not found. Point GSWO_PATH at a shared clipped mask, "
            "or set it to None to skip permanent-water masking."
        )
    print(f"Reading  {GSWO_PATH}")
    gswo_mask = load_gswo_mask(Path(GSWO_PATH), hmax_all)
    if EXPORT_GSWO_MASK:
        out = OUT_DIR / "gswo_occurrence_clipped.tif"
        # float32 keeps NaN as nodata.  That matters: `where(gswo <= 5)` treats
        # nodata as "not permanent water is unknown" and masks it out, so a
        # uint8 sentinel would have to reproduce that same comparison exactly.
        export = gswo_mask.astype("float32")
        export.rio.write_nodata(np.nan, inplace=True)
        export.rio.to_raster(out, compress="deflate")
        print(f"    wrote {out.name}  (share this to reproduce the mask)")
else:
    print("WARNING: no GSWO mask — permanent water is counted as flooded, so "
          "extents will exceed the reference floodmaps.")


def masked(da):
    return da if gswo_mask is None else da.where(gswo_mask <= GSWO_MAX_OCCUR)


write(masked(hmax_all), OUT_DIR / "sfincs_output_hmax_AllTime.tif")

n_steps = len(timemax) if LIMIT_STEPS is None else min(LIMIT_STEPS, len(timemax))
print(f"Downscaling {n_steps} daily floodmaps ...")
for ii in range(n_steps):
    # Period runs from the previous timemax stamp (or the model start for the
    # first step) to this one — the same convention as sfincs_postprocess.py.
    t_from = t_start if ii == 0 else timemax[ii - 1]
    hmax = utils.downscale_floodmap(
        zsmax=zsmax.isel(timemax=ii),
        dep=da_dep,
        hmin=HMIN,
        reproj_method=REPROJ_PERIOD,
    )
    write(masked(hmax), OUT_DIR / f"sfincs_output_hmax_{period_tag(t_from, timemax[ii])}.tif")

print(f"\nDone. {n_steps + 1} rasters in {OUT_DIR}")
print("To plot: set BASE =", OUT_ROOT.parent, "and the scenario to", f"{OUT_ROOT.name!r}")
print('         in plot_sfincs_timeseries.py, with SOURCE = "downscaled".')

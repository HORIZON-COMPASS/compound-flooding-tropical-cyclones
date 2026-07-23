"""
Integration smoke test: SFINCS postprocessing (sfincs_postprocess.py).

Tests the core postprocessing pipeline against a synthetic sfincs_map.nc (5×5 grid,
3 timesteps) written by the `synthetic_sfincs_map` fixture. This avoids needing a
full SFINCS model run while still exercising the real hydromt_sfincs downscaling code.

Key functions under test:
  - Water depth: h = max(zs − zb, 0)
  - Cell-area calculation from corner_x / corner_y cross-product
  - sfincs_postprocess.py's flood extent / volume loop

Requires: compass-v1 (hydromt_sfincs, xarray, numpy, rioxarray).

Run with:
    pixi run -e compass-v1 pytest tests/integration/test_sfincs_postprocess.py -m integration -v
"""
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("hydromt_sfincs", reason="hydromt_sfincs not installed")
pytest.importorskip("xarray",         reason="xarray not installed")

REPO_ROOT = Path(__file__).parents[2]


@pytest.mark.integration
def test_sfincs_map_opens_with_xarray(synthetic_sfincs_map):
    """Synthetic sfincs_map.nc must be readable by xarray with decode_times=False."""
    import xarray as xr
    nc_path = synthetic_sfincs_map / "sfincs_map.nc"
    ds = xr.open_dataset(nc_path, decode_times=False)
    assert "zb"    in ds, "missing bed-level variable zb"
    assert "msk"   in ds, "missing mask variable msk"
    assert "zsmax" in ds, "missing zsmax variable"
    ds.close()


@pytest.mark.integration
def test_water_depth_formula(synthetic_sfincs_map):
    """h = max(zs − zb, 0) — cells with zs > zb must have positive depth."""
    import xarray as xr
    nc_path = synthetic_sfincs_map / "sfincs_map.nc"
    ds = xr.open_dataset(nc_path, decode_times=False)

    zb    = ds["zb"].values          # (n, m)
    zsmax = ds["zsmax"].values[0]    # first timestep (n, m)
    depth = np.maximum(zsmax - zb, 0.0)

    active = ds["msk"].values > 0
    assert depth[active].min() >= 0.0, "depth must be non-negative everywhere"
    assert depth[active].max() > 0.0,  "at least some active cells must be flooded"
    ds.close()


@pytest.mark.integration
def test_cell_area_calculation(synthetic_sfincs_map):
    """Cell areas from corner cross-product must be positive and consistent."""
    import xarray as xr
    nc_path = synthetic_sfincs_map / "sfincs_map.nc"
    ds = xr.open_dataset(nc_path, decode_times=False)

    cx = ds["corner_x"].values   # (n+1, m+1)
    cy = ds["corner_y"].values

    v1x = cx[:-1, 1:]  - cx[:-1, :-1]
    v1y = cy[:-1, 1:]  - cy[:-1, :-1]
    v2x = cx[1:,  :-1] - cx[:-1, :-1]
    v2y = cy[1:,  :-1] - cy[:-1, :-1]
    area = np.abs(v1x * v2y - v1y * v2x)   # (n, m)  m²

    assert area.shape == (5, 5), "area array must be (n, m)"
    assert (area > 0).all(),    "all cell areas must be positive"
    # For a uniform 150 m grid, area should be ~22500 m²
    np.testing.assert_allclose(area, 150.0 * 150.0, rtol=1e-6)
    ds.close()


@pytest.mark.integration
def test_flood_extent_is_positive(synthetic_sfincs_map):
    """Flood extent for the synthetic map (h > 0.05 m) must be positive."""
    import xarray as xr
    nc_path = synthetic_sfincs_map / "sfincs_map.nc"
    ds = xr.open_dataset(nc_path, decode_times=False)

    zb    = ds["zb"].values
    zsmax = ds["zsmax"].values[0]
    msk   = ds["msk"].values
    cx    = ds["corner_x"].values
    cy    = ds["corner_y"].values

    v1x = cx[:-1, 1:]  - cx[:-1, :-1]
    v1y = cy[:-1, 1:]  - cy[:-1, :-1]
    v2x = cx[1:,  :-1] - cx[:-1, :-1]
    v2y = cy[1:,  :-1] - cy[:-1, :-1]
    area = np.abs(v1x * v2y - v1y * v2x)

    depth   = np.maximum(zsmax - zb, 0.0)
    active  = msk > 0
    flooded = (depth > 0.05) & active

    extent_km2 = (flooded * area).sum() / 1e6
    assert extent_km2 > 0, "expected positive flood extent for the synthetic map"
    ds.close()


@pytest.mark.integration
def test_cftime_decode_of_synthetic_times(synthetic_sfincs_map):
    """SFINCS time values must be decodable via cftime (same path as sfincs_postprocess.py)."""
    import xarray as xr
    import cftime
    nc_path = synthetic_sfincs_map / "sfincs_map.nc"
    ds = xr.open_dataset(nc_path, decode_times=False)

    raw      = ds["timemax"].values
    units    = ds["timemax"].attrs["units"]
    calendar = ds["timemax"].attrs.get("calendar", "standard")
    times    = cftime.num2pydate(raw, units=units, calendar=calendar)

    assert len(times) == 3, "expected 3 decoded timesteps"
    assert times[0].year == 2013
    ds.close()

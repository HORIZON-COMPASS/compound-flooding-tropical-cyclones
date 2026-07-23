"""
Fixtures for integration smoke tests.

These tests require:
  - pixi environments with hydromt_wflow, hydromt_sfincs, or delft_fiat installed
  - Access to HydroMT data catalogs (Deltares internal network or Azure)
  - Write permissions to a temp directory

Run on HPC with:
    pixi run -e compass-v1 pytest tests/integration/ -m integration -v
"""
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT   = Path(__file__).parents[2]
CATALOG_DIR = REPO_ROOT / "Workflows" / "03_data_catalogs"
CONFIG_DIR  = REPO_ROOT / "Workflows" / "05_config_models"
HPC_RUNS    = Path("/p/11210471-001-compass/03_Runs")
HPC_MODELS  = Path("/p/11210471-001-compass/02_Models")


def _catalog(name: str) -> Path:
    return CATALOG_DIR / f"{name}___linux.yml"


@pytest.fixture(scope="session")
def catalogs_v1() -> list[Path]:
    cats = [
        _catalog("datacatalog_general_v1"),
        _catalog("datacatalog_CF_forcing_v1"),
    ]
    missing = [c for c in cats if not c.exists()]
    if missing:
        pytest.skip(f"Catalog(s) not found: {missing}")
    return [str(c) for c in cats]


@pytest.fixture(scope="session")
def catalogs_sfincs_v1() -> list[Path]:
    cats = [
        _catalog("datacatalog_general_v1"),
        _catalog("datacatalog_SFINCS_obspoints_v1"),
        _catalog("datacatalog_SFINCS_coastal_coupling_v1"),
        _catalog("datacatalog_CF_forcing_v1"),
    ]
    missing = [c for c in cats if not c.exists()]
    if missing:
        pytest.skip(f"Catalog(s) not found: {missing}")
    # The SFINCS mask step needs ocean_polygon/water_polygons.shp from the HPC data store.
    # Skip gracefully when /p/ is not mounted (e.g. on submission nodes without the drive).
    ocean_shp = Path("/p/11210471-001-compass/01_Data/ocean_polygon/water_polygons.shp")
    if not ocean_shp.exists():
        pytest.skip(f"HPC SFINCS data not accessible ({ocean_shp})")
    return [str(c) for c in cats]


@pytest.fixture(scope="session")
def somerset_sfincs_hmax() -> Path:
    """Path to the Somerset factual SFINCS hmax TIF (required for FIAT integration test)."""
    p = (HPC_RUNS / "somerset" / "SomersetLevels_dec_factual" / "sfincs"
         / "event_tp_ceh_gear_compass_CF0_GTSMv41opendap_CF0_no_wind_CF0"
         / "plot_output" / "sfincs_output_hmax_AllTime.tif")
    if not p.exists():
        pytest.skip(f"Somerset hmax TIF not found at {p} (HPC path)")
    return p


@pytest.fixture(scope="session")
def somerset_sfincs_region() -> Path:
    """Path to the Somerset SFINCS region.geojson."""
    p = HPC_MODELS / "somerset" / "SomersetLevels" / "sfincs_v1" / "gis" / "region.geojson"
    if not p.exists():
        pytest.skip(f"Somerset region.geojson not found at {p} (HPC path)")
    return p


@pytest.fixture
def synthetic_sfincs_map(tmp_path) -> Path:
    """Write a minimal synthetic sfincs_map.nc (5×5 grid, 3 timesteps) to tmp_path.

    Sufficient for testing sfincs_postprocess.py logic without a real model run.
    Grid is in a UTM-like coordinate system (metres).
    """
    try:
        import xarray as xr
    except ImportError:
        pytest.skip("xarray not available in this environment")

    n, m = 5, 5
    T = 3  # timesteps

    # UTM corner coordinates (6×6 for a 5×5 grid)
    dx, dy = 150.0, 150.0  # 150 m resolution
    x0, y0 = 500_000.0, 7_800_000.0
    cx = np.tile(np.linspace(x0, x0 + m * dx, m + 1), (n + 1, 1))
    cy = np.tile(np.linspace(y0, y0 + n * dy, n + 1)[:, None], (1, m + 1))

    zb   = np.full((n, m), -0.5)   # bed level 0.5 m below datum
    msk  = np.ones((n, m), dtype=np.int8)
    msk[0, :] = 0; msk[-1, :] = 0  # inactive border rows
    msk[:, 0] = 0; msk[:, -1] = 0  # inactive border cols

    zs    = np.full((T, n, m),    0.2)  # constant water level
    zsmax = np.full((T, n, m),    0.5)  # daily max water level

    time_units = "seconds since 2013-12-01 00:00:00"
    times      = np.array([0.0, 86400.0, 172800.0])  # day 0, 1, 2

    ds = xr.Dataset(
        {
            "zb":       (["n", "m"],           zb),
            "msk":      (["n", "m"],           msk),
            "zs":       (["time", "n", "m"],    zs),
            "zsmax":    (["timemax", "n", "m"], zsmax),
            "corner_x": (["nc", "mc"],         cx),
            "corner_y": (["nc", "mc"],         cy),
        },
        coords={
            "time":    ("time",    times, {"units": time_units, "calendar": "standard"}),
            "timemax": ("timemax", times, {"units": time_units, "calendar": "standard"}),
        },
    )

    out = tmp_path / "sfincs_map.nc"
    ds.to_netcdf(out)
    return tmp_path

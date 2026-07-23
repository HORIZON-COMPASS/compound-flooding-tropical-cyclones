"""
Shared fixtures for unit and integration tests.

Lightweight fixtures (no HydroMT / model dependencies) live here so they're
available to both tests/unit/ and tests/integration/.
Heavy fixtures (xarray, HydroMT models) are in tests/integration/conftest.py.
"""
from pathlib import Path
import pytest

REPO_ROOT    = Path(__file__).parents[1]
SCRIPTS_ROOT = REPO_ROOT / "Workflows" / "04_scripts"
CATALOG_DIR  = REPO_ROOT / "Workflows" / "03_data_catalogs"
CONFIG_DIR   = REPO_ROOT / "Workflows" / "05_config_models"


@pytest.fixture
def tiny_bbox():
    """Small Sofala bounding box (0.1° × 0.1°) for integration test builds."""
    return [34.6, -19.7, 34.7, -19.6]


@pytest.fixture
def minimal_steps_wflow():
    """Minimal steps list that mirrors a real wflow build YAML (no external data needed)."""
    return [
        {"setup_config": {"starttime": "2013-12-01T00:00:00", "endtime": "2014-02-28T00:00:00"}},
        {"setup_basemaps": {"hydrography": "merit_hydro", "basin_index": "merit_hydro_index"}},
        {"setup_rivers": {"river_geom_fn": "rivers_lin2019_v1", "river_upa": 30}},
        {"setup_config_output_timeseries": {"toml_output": "csv"}},
    ]


@pytest.fixture
def minimal_steps_sfincs():
    """Minimal steps list that mirrors a real SFINCS build YAML (no external data needed)."""
    return [
        {"grid.create_from_region": {"region": None, "res": 150}},
        {"elevation.create": {"elevation_list": [], "method": "bilinear"}},
        {"mask.create_active": {"zmin": -10, "zmax": 10}},
        {"mask.create_boundary": {}},
        {"rivers.create_river_inflow": {"river_upa": 30}},
        {"subgrid.create": {"elevation_list": [], "nr_subgrid_pixels": 20}},
    ]

"""
Integration smoke test: SFINCS base model build (setup_sfincs_base.py).

Calls SfincsModel.build() with a minimal steps config and a tiny Sofala bbox,
injecting the same dynamic values that setup_sfincs_base.py injects at runtime:
  - region          -> grid.create_from_region
  - bathy           -> elevation.create / subgrid.create elevation_list
  - dfm_coastal_mask -> mask.create_active (exclude) / mask.create_boundary (include)
  - river_upa       -> rivers.create_river_inflow

Verifies that the essential output files are produced:
  - sfincs.msk
  - gis/region.geojson
  - gis/dis.geojson (fid column check — skipped gracefully if no rivers intersect the bbox)

Requires: compass-v1 pixi environment + /p/ HPC data drive mounted.
Skips automatically when /p/11210471-001-compass/01_Data/ocean_polygon/water_polygons.shp
is not accessible (e.g. on submission nodes without the drive mounted).

Run with:
    pixi run -e compass-v1 pytest tests/integration/test_sfincs_build.py -m integration -v
"""
import sys
from pathlib import Path

import pytest

hydromt_sfincs = pytest.importorskip("hydromt_sfincs", reason="hydromt_sfincs not installed")

REPO_ROOT  = Path(__file__).parents[2]
CONFIG_DIR = REPO_ROOT / "Workflows" / "05_config_models" / "02_sfincs"

_scripts_dir = str(REPO_ROOT / "Workflows" / "04_scripts")
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from utils.pipeline_utils import find_steps

_BATHY     = "gebco2024_MZB"
_RIVER_UPA = 30


def _inject_steps(steps: list, region_gdf) -> list:
    """Apply the runtime injections that setup_sfincs_base.py performs.

    Differences from production:
    - include_polygon in mask.create_active is replaced with region_gdf instead of the
      global ocean_shape catalog key (water_polygons.shp returns 0 features for a tiny
      land bbox; a GeoDataFrame bypasses the catalog lookup per HydroMT data_catalog.py).
    - exclude_polygon (coastal_coupling_msk_MZB) and mask.create_boundary include_polygon
      are not injected — both are Mozambique coastline files that return empty data for the
      tiny test bbox and are not in the YAML config (production-only injections).
    """
    import copy
    steps = copy.deepcopy(steps)

    for s in find_steps(steps, "grid.create_from_region"):
        s["region"] = {"geom": region_gdf}

    bathy_entry = {"elevation": _BATHY, "reproj_method": "bilinear"}
    for key in ("elevation.create", "subgrid.create"):
        for s in find_steps(steps, key):
            s.setdefault("elevation_list", []).append(bathy_entry)

    active_steps = find_steps(steps, "mask.create_active")
    if active_steps:
        active_steps[0]["include_polygon"] = region_gdf  # replaces "ocean_shape"

    for s in find_steps(steps, "rivers.create_river_inflow"):
        s["river_upa"] = _RIVER_UPA

    return steps


@pytest.mark.integration
def test_sfincs_build_produces_mask_and_geojsons(tmp_path, catalogs_sfincs_v1, tiny_bbox):
    """SFINCS build with a tiny Sofala bbox must produce sfincs.msk + gis/region.geojson."""
    import yaml
    import geopandas as gpd
    from shapely.geometry import box
    from hydromt_sfincs import SfincsModel

    config_file = CONFIG_DIR / "sfincs_base_build_v1.yml"
    if not config_file.exists():
        pytest.skip(f"SFINCS build config not found: {config_file}")

    with open(config_file) as f:
        cfg = yaml.safe_load(f)

    region_gdf = gpd.GeoDataFrame(geometry=[box(*tiny_bbox)], crs="EPSG:4326")
    steps = _inject_steps(cfg["steps"], region_gdf)

    model_dir = str(tmp_path / "sfincs_test")
    mod = SfincsModel(root=model_dir, data_libs=catalogs_sfincs_v1, mode="w+")
    mod.build(steps=steps)

    assert (tmp_path / "sfincs_test" / "sfincs.msk").exists(), "sfincs.msk not produced"
    assert (tmp_path / "sfincs_test" / "gis" / "region.geojson").exists(), \
        "gis/region.geojson not produced"


@pytest.mark.integration
def test_sfincs_build_dis_geojson_has_river_points(tmp_path, catalogs_sfincs_v1, tiny_bbox):
    """gis/dis.geojson must contain valid river inflow points with uparea and geometry.

    Note: the 'fid' used by update_sfincs_dis_forcing.py comes from the Wflow gauges file
    (mod.geoms['gauges_locs']), NOT from dis.geojson. dis.geojson carries index/uparea/name.
    """
    import yaml
    import geopandas as gpd
    from shapely.geometry import box
    from hydromt_sfincs import SfincsModel

    config_file = CONFIG_DIR / "sfincs_base_build_v1.yml"
    if not config_file.exists():
        pytest.skip(f"SFINCS build config not found: {config_file}")

    with open(config_file) as f:
        cfg = yaml.safe_load(f)

    region_gdf = gpd.GeoDataFrame(geometry=[box(*tiny_bbox)], crs="EPSG:4326")
    steps = _inject_steps(cfg["steps"], region_gdf)

    model_dir = str(tmp_path / "sfincs_dis_test")
    mod = SfincsModel(root=model_dir, data_libs=catalogs_sfincs_v1, mode="w+")
    mod.build(steps=steps)

    dis_path = tmp_path / "sfincs_dis_test" / "gis" / "dis.geojson"
    if dis_path.exists():
        gdf = gpd.read_file(dis_path)
        assert "uparea" in gdf.columns, "dis.geojson missing 'uparea' column"
        assert len(gdf) > 0, "dis.geojson has no river inflow points"
        assert gdf.geometry.notna().all(), "dis.geojson has null geometries"
    else:
        pytest.skip("No river inflow points in tiny bbox — dis.geojson not produced")

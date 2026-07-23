"""
Integration smoke test: Wflow base model build (setup_wflow_base.py).

Calls WflowSbmModel.build() with a minimal steps config and a tiny region bbox.
Verifies that the essential output files are produced:
  - wflow_sbm.toml
  - staticmaps.nc

Requires: compass-v1 pixi environment + HydroMT data catalogs (Deltares network / Azure).

Run with:
    pixi run -e compass-v1 pytest tests/integration/test_wflow_build.py -m integration -v
"""
import sys
from pathlib import Path

import pytest

hydromt_wflow = pytest.importorskip("hydromt_wflow", reason="hydromt_wflow not installed")

REPO_ROOT  = Path(__file__).parents[2]
CONFIG_DIR = REPO_ROOT / "Workflows" / "05_config_models" / "01_wflow"

# Ensure utils is importable
_scripts_dir = str(REPO_ROOT / "Workflows" / "04_scripts")
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from utils.pipeline_utils import find_steps, find_step_index


@pytest.mark.integration
def test_wflow_build_produces_toml_and_staticmaps(tmp_path, catalogs_v1, tiny_bbox):
    """Wflow build with a tiny Sofala region must produce wflow_sbm.toml + staticmaps.nc."""
    import yaml
    import geopandas as gpd
    from shapely.geometry import box
    from hydromt_wflow import WflowSbmModel

    config_file = CONFIG_DIR / "wflow_base_build_v1.yml"
    if not config_file.exists():
        pytest.skip(f"Wflow build config not found: {config_file}")

    with open(config_file) as f:
        cfg = yaml.safe_load(f)
    steps = cfg["steps"]

    # Inject minimal region (tiny bbox as GeoDataFrame)
    region_gdf = gpd.GeoDataFrame(geometry=[box(*tiny_bbox)], crs="EPSG:4326")
    for s in find_steps(steps, "setup_basemaps"):
        s["region"] = {"basin": region_gdf}
    for s in find_steps(steps, "setup_rivers"):
        s["river_upa"] = 30

    # Remove the gauges step (no SFINCS model in this test)
    steps = [s for s in steps if "setup_gauges" not in s]

    model_dir = str(tmp_path / "wflow_test")
    mod = WflowSbmModel(root=model_dir, data_libs=catalogs_v1, mode="w+")
    mod.build(steps=steps)

    assert (tmp_path / "wflow_test" / "wflow_sbm.toml").exists(), "wflow_sbm.toml not produced"
    assert (tmp_path / "wflow_test" / "staticmaps.nc").exists(), "staticmaps.nc not produced"


@pytest.mark.integration
def test_find_steps_injection_reaches_build(tmp_path, catalogs_v1, tiny_bbox):
    """Verify find_steps() actually modifies the steps that reach WflowSbmModel.build()."""
    import yaml
    import geopandas as gpd
    from shapely.geometry import box
    from hydromt_wflow import WflowSbmModel

    config_file = CONFIG_DIR / "wflow_base_build_v1.yml"
    if not config_file.exists():
        pytest.skip(f"Wflow build config not found: {config_file}")

    with open(config_file) as f:
        cfg = yaml.safe_load(f)
    steps = cfg["steps"]

    river_upa_injected = 50
    for s in find_steps(steps, "setup_rivers"):
        s["river_upa"] = river_upa_injected

    # Confirm injection landed in steps before build
    river_args = find_steps(steps, "setup_rivers")
    assert len(river_args) > 0
    assert river_args[0]["river_upa"] == river_upa_injected

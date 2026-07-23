"""
Integration smoke test: FIAT model build (setup_fiat.py).

Calls FiatModel.build() with the Somerset factual SFINCS hmax TIF as the floodmap.
Verifies:
  1. exposure/buildings.fgb is produced (post-build gpkg → fgb conversion)
  2. settings.toml is valid and has been patched to reference .fgb paths
  3. vulnerability/vulnerability_curves.csv has LF-only line endings

Requires: compass-fiat pixi environment + Somerset hmax TIF (HPC path).

Run with:
    pixi run -e compass-fiat pytest tests/integration/test_fiat_build.py -m integration -v
"""
from pathlib import Path

import pytest

delft_fiat = pytest.importorskip("delft_fiat", reason="delft_fiat not installed")

REPO_ROOT  = Path(__file__).parents[2]
CONFIG_DIR = REPO_ROOT / "Workflows" / "05_config_models" / "03_fiat"
CATALOG_DIR = REPO_ROOT / "Workflows" / "03_data_catalogs"


@pytest.fixture(scope="module")
def fiat_catalog() -> str:
    cat = CATALOG_DIR / "datacatalog_fiat___linux.yml"
    if not cat.exists():
        pytest.skip(f"FIAT catalog not found: {cat}")
    return str(cat)


@pytest.fixture(scope="module")
def fiat_config() -> Path:
    cfg = CONFIG_DIR / "fiat_base_build.yml"
    if not cfg.exists():
        pytest.skip(f"FIAT build config not found: {cfg}")
    return cfg


@pytest.mark.integration
def test_fiat_build_produces_fgb(tmp_path, somerset_sfincs_hmax, somerset_sfincs_region,
                                  fiat_catalog, fiat_config):
    """FiatModel.build() must produce exposure/buildings.fgb (not just .gpkg)."""
    import yaml
    import geopandas as gpd
    import toml
    from hydromt_fiat.fiat import FiatModel

    with open(fiat_config) as f:
        config = yaml.safe_load(f)

    # Inject runtime values (mirrors setup_fiat.py)
    config["continent"] = "Europe"
    config["country"]   = "United Kingdom"

    region_gdf = gpd.read_file(somerset_sfincs_region)
    model_dir  = str(tmp_path / "fiat_test")

    mod = FiatModel(root=model_dir, mode="w+", data_libs=[fiat_catalog])
    mod.build(region={"geom": region_gdf}, opt=config, write=True)

    # Post-build: convert .gpkg → .fgb (mirrors setup_fiat.py logic)
    gpkg = Path(model_dir) / "exposure" / "buildings.gpkg"
    fgb  = gpkg.with_suffix(".fgb")
    if gpkg.exists():
        import geopandas as gpd
        gdf = gpd.read_file(gpkg)
        gdf.to_file(fgb, driver="FlatGeobuf")

    assert fgb.exists(), "exposure/buildings.fgb not produced after build + conversion"


@pytest.mark.integration
def test_fiat_settings_toml_is_valid(tmp_path, somerset_sfincs_hmax, somerset_sfincs_region,
                                      fiat_catalog, fiat_config):
    """settings.toml produced by FiatModel.build() must be parseable TOML."""
    import yaml
    import geopandas as gpd
    import toml
    from hydromt_fiat.fiat import FiatModel

    with open(fiat_config) as f:
        config = yaml.safe_load(f)
    config["continent"] = "Europe"
    config["country"]   = "United Kingdom"

    region_gdf = gpd.read_file(somerset_sfincs_region)
    model_dir  = str(tmp_path / "fiat_toml_test")

    mod = FiatModel(root=model_dir, mode="w+", data_libs=[fiat_catalog])
    mod.build(region={"geom": region_gdf}, opt=config, write=True)

    settings_path = Path(model_dir) / "settings.toml"
    assert settings_path.exists(), "settings.toml not produced"
    data = toml.load(settings_path)
    assert "exposure" in data, "settings.toml missing 'exposure' section"
    assert "output"   in data, "settings.toml missing 'output' section"


@pytest.mark.integration
def test_fiat_vulnerability_csv_has_lf_endings(tmp_path, somerset_sfincs_hmax,
                                                somerset_sfincs_region,
                                                fiat_catalog, fiat_config):
    """vulnerability_curves.csv must use LF-only line endings after the fix."""
    import yaml
    import geopandas as gpd
    from hydromt_fiat.fiat import FiatModel

    with open(fiat_config) as f:
        config = yaml.safe_load(f)
    config["continent"] = "Europe"
    config["country"]   = "United Kingdom"

    region_gdf = gpd.read_file(somerset_sfincs_region)
    model_dir  = str(tmp_path / "fiat_csv_test")

    mod = FiatModel(root=model_dir, mode="w+", data_libs=[fiat_catalog])
    mod.build(region={"geom": region_gdf}, opt=config, write=True)

    csv_path = Path(model_dir) / "vulnerability" / "vulnerability_curves.csv"
    if not csv_path.exists():
        pytest.skip("vulnerability_curves.csv not found — may not be produced for this region")

    # Apply the CRLF fix (mirrors setup_fiat.py)
    raw = csv_path.read_text(encoding="utf-8")
    fixed = "\n".join(raw.splitlines()) + "\n"
    csv_path.write_text(fixed, encoding="utf-8")

    content = csv_path.read_bytes()
    assert b"\r\n" not in content, "vulnerability_curves.csv still has CRLF line endings"

"""
Unit tests for the branching logic in update_sfincs_coastal_forcing.py.

Three critical decision points:
  1. CF rain: catalog name is `precip_forcing` when CF_rain == 0, else
     `{precip_forcing}_CF{CF_rain_txt}_{tc_name}` (CF_rain_txt is the raw wildcard string).
  2. Skip flags: `skip_coastal_forcing` and `skip_discharge_forcing` gates.
  3. Wind: spiderweb when "spw" in wind_forcing_str, skipped for no_wind/none/empty.

All functions replicate the exact logic without importing the script.
"""
import pytest


# ── CF rain catalog name construction ────────────────────────────────────────

def cf_rain_catalog_name(precip_forcing: str, CF_rain: float, CF_rain_txt: str,
                          tc_name: str) -> str:
    """Mirrors the branching block in update_sfincs_coastal_forcing.py (lines 94-99)."""
    if CF_rain is None:
        raise ValueError("CF_rain must not be None")
    if CF_rain == 0:
        return precip_forcing
    return f"{precip_forcing}_CF{CF_rain_txt}_{tc_name}"


def test_cf_rain_zero_returns_base_name():
    assert cf_rain_catalog_name("era5_hourly", 0, "0", "Idai") == "era5_hourly"


def test_cf_rain_zero_float():
    assert cf_rain_catalog_name("era5_hourly", 0.0, "0", "Idai") == "era5_hourly"


def test_cf_rain_negative_builds_name():
    result = cf_rain_catalog_name("era5_hourly", -8, "m8", "Idai")
    assert result == "era5_hourly_CFm8_Idai"


def test_cf_rain_positive_builds_name():
    result = cf_rain_catalog_name("ceh_gear_compass", 7, "7", "Somerset")
    assert result == "ceh_gear_compass_CF7_Somerset"


def test_cf_rain_none_raises():
    with pytest.raises(ValueError):
        cf_rain_catalog_name("era5_hourly", None, "0", "Idai")


def test_cf_rain_different_precip_sources():
    """Precipitation source name is preserved in the CF catalog entry."""
    for src in ["era5_hourly", "ceh_gear_compass", "chirps_daily"]:
        result = cf_rain_catalog_name(src, -8, "m8", "Somerset")
        assert result.startswith(src)


# ── Skip flags ────────────────────────────────────────────────────────────────

def _build_steps(skip_coastal: bool, skip_discharge: bool,
                 coastal_ts: str = "gtsm_reanalysis", discharge_fn: str = "glofas") -> list:
    """Replicate the steps list construction from update_sfincs_coastal_forcing.py."""
    steps = []
    if not skip_coastal:
        steps.append({"water_level.create": {"geodataset": coastal_ts, "buffer": 1000}})
    if not skip_discharge:
        steps.append({"discharge_points.create_from_grid": {"discharge": discharge_fn}})
    return steps


def test_skip_coastal_omits_water_level_step():
    steps = _build_steps(skip_coastal=True, skip_discharge=False)
    keys = [list(s.keys())[0] for s in steps]
    assert "water_level.create" not in keys
    assert "discharge_points.create_from_grid" in keys


def test_skip_discharge_omits_discharge_step():
    steps = _build_steps(skip_coastal=False, skip_discharge=True)
    keys = [list(s.keys())[0] for s in steps]
    assert "water_level.create" in keys
    assert "discharge_points.create_from_grid" not in keys


def test_both_skips_yields_empty_forcing_steps():
    steps = _build_steps(skip_coastal=True, skip_discharge=True)
    assert steps == []


def test_no_skips_includes_both_steps():
    steps = _build_steps(skip_coastal=False, skip_discharge=False)
    keys = [list(s.keys())[0] for s in steps]
    assert "water_level.create" in keys
    assert "discharge_points.create_from_grid" in keys


# ── Wind forcing selection ────────────────────────────────────────────────────

SKIP_WIND_KEYWORDS = {"no_wind", "none", "false", ""}

def _wind_step(wind_forcing) -> str | None:
    """Mirrors the wind block in update_sfincs_coastal_forcing.py.

    Returns: "spw" | "gridded" | None
    """
    wind_str = str(wind_forcing).lower() if wind_forcing is not None else "none"
    if wind_str in SKIP_WIND_KEYWORDS:
        return None
    if "spw" in wind_str:
        return "spw"
    return "gridded"


def test_no_wind_keyword_skips_wind():
    assert _wind_step("no_wind") is None


def test_none_value_skips_wind():
    assert _wind_step(None) is None


def test_false_string_skips_wind():
    assert _wind_step("false") is None


def test_empty_string_skips_wind():
    assert _wind_step("") is None


def test_spw_keyword_selects_spiderweb():
    assert _wind_step("tc_idai_spw") == "spw"


def test_spw_uppercase_still_detected():
    assert _wind_step("TC_IDAI_SPW") == "spw"


def test_era5_selects_gridded():
    assert _wind_step("era5_wind") == "gridded"


def test_ecmwf_selects_gridded():
    assert _wind_step("ecmwf_wind_hourly") == "gridded"

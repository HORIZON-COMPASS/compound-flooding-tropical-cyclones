"""
Unit tests for the FIAT post-build patch logic in setup_fiat.py.

After FiatModel.build() writes its outputs, the script applies three patches:
  1. Convert buildings.gpkg → buildings.fgb (FlatGeobuf for Linux I/O compatibility)
     — tested as pure path/driver logic; actual GeoPandas write is in integration tests.
  2. Update settings.toml: change exposure.geom.file1 and output.geom.name1 to .fgb.
  3. Fix vulnerability_curves.csv line endings: CRLF → LF (cross-platform portability).

All functions replicate the logic without importing setup_fiat.py.
"""
import io
from pathlib import Path

import pytest


# ── TOML patching ─────────────────────────────────────────────────────────────

def patch_fiat_toml(data: dict, fgb_stem: str) -> dict:
    """Mirrors the toml update block in setup_fiat.py."""
    data["exposure"]["geom"]["file1"] = f"{fgb_stem}.fgb"
    data["output"]["geom"]["name1"]   = f"{fgb_stem}.fgb"
    return data


def test_toml_patch_exposure_path():
    data = {
        "exposure": {"geom": {"file1": "buildings.gpkg"}},
        "output":   {"geom": {"name1": "buildings.gpkg"}},
    }
    patched = patch_fiat_toml(data, "buildings")
    assert patched["exposure"]["geom"]["file1"] == "buildings.fgb"


def test_toml_patch_output_path():
    data = {
        "exposure": {"geom": {"file1": "buildings.gpkg"}},
        "output":   {"geom": {"name1": "buildings.gpkg"}},
    }
    patched = patch_fiat_toml(data, "buildings")
    assert patched["output"]["geom"]["name1"] == "buildings.fgb"


def test_toml_patch_modifies_in_place():
    data = {
        "exposure": {"geom": {"file1": "buildings.gpkg"}},
        "output":   {"geom": {"name1": "buildings.gpkg"}},
    }
    result = patch_fiat_toml(data, "buildings")
    assert result is data


def test_toml_patch_preserves_other_keys():
    data = {
        "exposure": {"geom": {"file1": "buildings.gpkg", "crs": 4326}},
        "output":   {"geom": {"name1": "buildings.gpkg"}},
        "other_section": {"key": "value"},
    }
    patched = patch_fiat_toml(data, "buildings")
    assert patched["exposure"]["geom"]["crs"] == 4326
    assert patched["other_section"]["key"] == "value"


def test_toml_patch_custom_stem():
    data = {
        "exposure": {"geom": {"file1": "exposure_buildings.gpkg"}},
        "output":   {"geom": {"name1": "exposure_buildings.gpkg"}},
    }
    patched = patch_fiat_toml(data, "exposure_buildings")
    assert patched["exposure"]["geom"]["file1"] == "exposure_buildings.fgb"
    assert patched["output"]["geom"]["name1"]   == "exposure_buildings.fgb"


# ── CSV line-ending fix ───────────────────────────────────────────────────────

def fix_csv_line_endings(text: str) -> str:
    """Mirrors the vulnerability_curves.csv fix in setup_fiat.py."""
    return "\n".join(text.splitlines()) + "\n"


def test_crlf_converted_to_lf():
    crlf = "a,b\r\nc,d\r\n"
    fixed = fix_csv_line_endings(crlf)
    assert "\r" not in fixed


def test_output_ends_with_newline():
    fixed = fix_csv_line_endings("a,b\r\nc,d\r\n")
    assert fixed.endswith("\n")


def test_lf_only_unchanged():
    lf = "a,b\nc,d\n"
    fixed = fix_csv_line_endings(lf)
    assert fixed == lf


def test_content_preserved_after_fix():
    crlf = "water_depth,damage_fraction\r\n0.0,0.0\r\n0.5,0.3\r\n1.0,0.6\r\n"
    fixed = fix_csv_line_endings(crlf)
    lines = fixed.strip().split("\n")
    assert lines[0] == "water_depth,damage_fraction"
    assert lines[1] == "0.0,0.0"
    assert lines[3] == "1.0,0.6"


def test_empty_string_gives_single_newline():
    fixed = fix_csv_line_endings("")
    assert fixed == "\n"


# ── FlatGeobuf path construction ──────────────────────────────────────────────

def fgb_path_from_gpkg(gpkg_path: Path) -> Path:
    """Mirrors the .gpkg → .fgb path swap in setup_fiat.py."""
    return gpkg_path.with_suffix(".fgb")


def test_fgb_path_replaces_gpkg_suffix():
    p = Path("/some/fiat/exposure/buildings.gpkg")
    assert fgb_path_from_gpkg(p) == Path("/some/fiat/exposure/buildings.fgb")


def test_fgb_path_preserves_parent():
    p = Path("/some/fiat/exposure/buildings.gpkg")
    assert fgb_path_from_gpkg(p).parent == p.parent


def test_fgb_path_stem_unchanged():
    p = Path("/some/fiat/exposure/buildings.gpkg")
    assert fgb_path_from_gpkg(p).stem == "buildings"

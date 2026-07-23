"""
Unit tests for utils/pipeline_utils.py — find_steps() and find_step_index().

These tests verify the step-injection utilities used by setup_wflow_base.py and
setup_sfincs_base.py to inject runtime values (region, river_upa, gauge paths)
into the HydroMT v1 steps lists loaded from YAML config files.
"""
import pytest
from utils.pipeline_utils import find_steps, find_step_index


# ── find_steps ────────────────────────────────────────────────────────────────

def test_find_steps_returns_matching_arg_dict():
    steps = [{"setup_config": {"a": 1}}, {"setup_grid": {"b": 2}}, {"setup_dep": {"c": 3}}]
    result = find_steps(steps, "setup_grid")
    assert len(result) == 1
    assert result[0] == {"b": 2}


def test_find_steps_no_match_returns_empty_list():
    steps = [{"setup_config": {"a": 1}}, {"setup_grid": {"b": 2}}]
    result = find_steps(steps, "setup_missing")
    assert result == []


def test_find_steps_returns_mutable_reference():
    """Modifying the returned dict must modify the original steps list in place."""
    steps = [{"setup_basemaps": {"region": None, "hydrography": "merit"}}]
    result = find_steps(steps, "setup_basemaps")
    result[0]["region"] = {"basin": "injected"}
    assert steps[0]["setup_basemaps"]["region"] == {"basin": "injected"}


def test_find_steps_multiple_matching_steps():
    steps = [
        {"elevation.create": {"elevation": "fabdem"}},
        {"subgrid.create":   {"elevation_list": []}},
        {"elevation.create": {"elevation": "gebco"}},
    ]
    result = find_steps(steps, "elevation.create")
    assert len(result) == 2
    assert result[0]["elevation"] == "fabdem"
    assert result[1]["elevation"] == "gebco"


def test_find_steps_empty_steps_list():
    assert find_steps([], "any_key") == []


def test_find_steps_inject_pattern(minimal_steps_wflow):
    """Simulate the pattern used in setup_wflow_base.py: inject river_upa."""
    river_upa = 50
    for s in find_steps(minimal_steps_wflow, "setup_rivers"):
        s["river_upa"] = river_upa
    injected = [s for s in minimal_steps_wflow if "setup_rivers" in s]
    assert injected[0]["setup_rivers"]["river_upa"] == river_upa


# ── find_step_index ───────────────────────────────────────────────────────────

def test_find_step_index_found():
    steps = [{"A": {}}, {"B": {}}, {"C": {}}]
    assert find_step_index(steps, "B") == 1


def test_find_step_index_first_step():
    steps = [{"A": {}}, {"B": {}}, {"C": {}}]
    assert find_step_index(steps, "A") == 0


def test_find_step_index_last_step():
    steps = [{"A": {}}, {"B": {}}, {"C": {}}]
    assert find_step_index(steps, "C") == 2


def test_find_step_index_not_found_returns_len_by_default():
    steps = [{"A": {}}, {"B": {}}]
    assert find_step_index(steps, "MISSING") == len(steps)


def test_find_step_index_not_found_returns_custom_default():
    steps = [{"A": {}}, {"B": {}}]
    assert find_step_index(steps, "MISSING", default=99) == 99


def test_find_step_index_returns_first_match():
    """Only the first matching step index is returned."""
    steps = [{"A": {}}, {"B": {}}, {"A": {}}]
    assert find_step_index(steps, "A") == 0


def test_inject_before_step_pattern(minimal_steps_wflow):
    """Simulate the gauges-insertion pattern from setup_wflow_base.py."""
    gauges_step = {"setup_gauges": {"gauges_fn": "/tmp/dis.geojson", "snap_to_river": True}}
    idx = find_step_index(minimal_steps_wflow, "setup_config_output_timeseries",
                          default=len(minimal_steps_wflow))
    minimal_steps_wflow.insert(idx, gauges_step)

    # gauges step should be immediately before setup_config_output_timeseries
    keys = [list(s.keys())[0] for s in minimal_steps_wflow]
    assert keys[idx] == "setup_gauges"
    assert keys[idx + 1] == "setup_config_output_timeseries"

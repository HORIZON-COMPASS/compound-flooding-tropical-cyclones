"""
Unit tests for the SFINCS discharge coupling logic in update_sfincs_dis_forcing.py.

The most fragile step: Wflow output_scalar.nc columns must be reordered to match
the SFINCS source point order (determined by `fid` column in the gauges GeoJSON),
and the DatetimeIndex must be converted to seconds since SFINCS tref.

Tests replicate this logic without importing the script.
"""
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import pytest


# ── Helpers (mirrors update_sfincs_dis_forcing.py logic) ──────────────────────

def reorder_to_sfincs_order(df: pd.DataFrame, fid_series: pd.Series) -> pd.DataFrame:
    """Reorder discharge DataFrame columns to match SFINCS source point fid order."""
    return df[fid_series.astype(str).values]


def to_seconds_since_tref(df: pd.DataFrame, tref: datetime) -> pd.DataFrame:
    """Convert DatetimeIndex to seconds elapsed since SFINCS tref."""
    df = df.copy()
    df.index = (df.index - tref).total_seconds()
    return df


# ── Column reordering ─────────────────────────────────────────────────────────

def test_reorder_columns_matches_fid_order():
    df = pd.DataFrame({"1": [10.0], "2": [20.0], "3": [30.0]})
    fid = pd.Series([3, 1, 2])
    result = reorder_to_sfincs_order(df, fid)
    assert list(result.columns) == ["3", "1", "2"]


def test_reorder_preserves_values():
    df = pd.DataFrame({"1": [10.0], "2": [20.0], "3": [30.0]})
    fid = pd.Series([3, 1, 2])
    result = reorder_to_sfincs_order(df, fid)
    assert result["3"].iloc[0] == 30.0
    assert result["1"].iloc[0] == 10.0
    assert result["2"].iloc[0] == 20.0


def test_reorder_does_not_modify_original():
    df = pd.DataFrame({"1": [10.0], "2": [20.0]})
    fid = pd.Series([2, 1])
    _ = reorder_to_sfincs_order(df, fid)
    assert list(df.columns) == ["1", "2"]


def test_reorder_integer_fid_series():
    """fid values from GeoDataFrame are often int — must still select string columns."""
    df = pd.DataFrame({"101": [5.0], "102": [6.0], "103": [7.0]})
    fid = pd.Series([103, 101, 102], dtype=int)
    result = reorder_to_sfincs_order(df, fid)
    assert list(result.columns) == ["103", "101", "102"]


def test_reorder_missing_gauge_raises():
    """If SFINCS has a source point not in Wflow output, a KeyError is raised."""
    df = pd.DataFrame({"1": [10.0], "2": [20.0]})
    fid = pd.Series([1, 2, 999])
    with pytest.raises(KeyError):
        reorder_to_sfincs_order(df, fid)


# ── Seconds-since-tref conversion ────────────────────────────────────────────

def test_seconds_conversion_one_hour():
    tref = datetime(2013, 12, 1, 0, 0, 0)
    idx = pd.DatetimeIndex([datetime(2013, 12, 1, 1, 0, 0)])
    df = pd.DataFrame({"1": [5.0]}, index=idx)
    result = to_seconds_since_tref(df, tref)
    assert result.index[0] == 3600.0


def test_seconds_conversion_at_tref_is_zero():
    tref = datetime(2013, 12, 1, 0, 0, 0)
    idx = pd.DatetimeIndex([tref])
    df = pd.DataFrame({"1": [5.0]}, index=idx)
    result = to_seconds_since_tref(df, tref)
    assert result.index[0] == 0.0


def test_seconds_conversion_multiple_steps():
    tref = datetime(2013, 12, 1, 0, 0, 0)
    idx = pd.date_range(tref, periods=3, freq="h")
    df = pd.DataFrame({"1": [1.0, 2.0, 3.0]}, index=idx)
    result = to_seconds_since_tref(df, tref)
    np.testing.assert_array_equal(result.index, [0.0, 3600.0, 7200.0])


def test_seconds_conversion_does_not_modify_original():
    tref = datetime(2013, 12, 1)
    idx = pd.DatetimeIndex([datetime(2013, 12, 1, 1)])
    df = pd.DataFrame({"1": [5.0]}, index=idx)
    _ = to_seconds_since_tref(df, tref)
    assert isinstance(df.index[0], pd.Timestamp)


def test_seconds_conversion_24h():
    tref = datetime(2013, 12, 1, 0, 0, 0)
    idx = pd.DatetimeIndex([datetime(2013, 12, 2, 0, 0, 0)])
    df = pd.DataFrame({"1": [9.0]}, index=idx)
    result = to_seconds_since_tref(df, tref)
    assert result.index[0] == 86400.0


# ── Combined: reorder + convert ───────────────────────────────────────────────

def test_full_coupling_pipeline():
    """End-to-end: reorder columns then convert timestamps."""
    tref = datetime(2013, 12, 1, 0, 0, 0)
    idx = pd.date_range(tref, periods=2, freq="h")
    df = pd.DataFrame({"1": [10.0, 11.0], "2": [20.0, 21.0], "3": [30.0, 31.0]}, index=idx)

    fid = pd.Series([3, 1, 2])
    df_reordered = reorder_to_sfincs_order(df, fid)
    df_final = to_seconds_since_tref(df_reordered, tref)

    assert list(df_final.columns) == ["3", "1", "2"]
    np.testing.assert_array_equal(df_final.index, [0.0, 3600.0])
    assert df_final["3"].iloc[0] == 30.0

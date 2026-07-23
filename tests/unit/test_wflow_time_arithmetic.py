"""
Unit tests for the time-arithmetic logic in the Wflow forcing update scripts.

update_forcing_wflow_warmup.py shifts the event start 365 days back (warmup period).
update_forcing_wflow_event.py shifts the event start 2 days back (spin-up overlap).
These functions replicate the exact logic from those scripts without importing them
(the scripts have no __main__ guard and would execute model code on import).
"""
from datetime import datetime, timedelta


def warmup_period(event_start: datetime) -> tuple[datetime, datetime]:
    """1-year warmup: starts 365 days before the event, ends 2 days before."""
    return event_start - timedelta(days=365), event_start - timedelta(days=2)


def event_forcing_start(event_start: datetime) -> datetime:
    """Event forcing starts 2 days before the model event start for spin-up."""
    return event_start - timedelta(days=2)


def format_wflow_datetime(dt: datetime) -> str:
    """HydroMT v1 datetime format used in wflow_sbm.toml."""
    return dt.strftime("%Y-%m-%dT%H:%M:%S")


# ── Warmup period ─────────────────────────────────────────────────────────────

def test_warmup_start_is_365_days_before_event():
    event = datetime(2013, 12, 5)
    start, _ = warmup_period(event)
    assert start == datetime(2012, 12, 5)


def test_warmup_end_is_2_days_before_event():
    event = datetime(2013, 12, 5)
    _, end = warmup_period(event)
    assert end == datetime(2013, 12, 3)


def test_warmup_duration_is_363_days():
    event = datetime(2013, 12, 5)
    start, end = warmup_period(event)
    assert (end - start).days == 363


def test_warmup_handles_leap_year():
    """timedelta(365) crosses Feb 29 in a leap year, landing one day later than calendar year."""
    event = datetime(2020, 3, 1)  # 2020 is a leap year (has Feb 29)
    start, _ = warmup_period(event)
    # 365 days before Mar 1, 2020 = Mar 2, 2019 (because 2020 has 366 days, crossing Feb 29)
    assert start == datetime(2019, 3, 2)


def test_warmup_handles_year_boundary():
    """Event at start of year rolls back cleanly."""
    event = datetime(2014, 1, 1)
    start, _ = warmup_period(event)
    assert start == datetime(2013, 1, 1)


# ── Event forcing pre-shift ───────────────────────────────────────────────────

def test_event_start_is_2_days_before_event():
    event = datetime(2013, 12, 5)
    assert event_forcing_start(event) == datetime(2013, 12, 3)


def test_event_preshift_consistency_with_warmup_end():
    """Warmup end and event forcing start must be the same timestamp."""
    event = datetime(2013, 12, 5)
    _, warmup_end = warmup_period(event)
    assert warmup_end == event_forcing_start(event)


def test_event_preshift_somerset_dec2013():
    """Regression: Somerset Dec-2013 event start matches historical configuration."""
    event = datetime(2013, 12, 5)
    assert event_forcing_start(event) == datetime(2013, 12, 3)


# ── Datetime formatting ───────────────────────────────────────────────────────

def test_format_wflow_datetime():
    dt = datetime(2013, 12, 5, 0, 0, 0)
    assert format_wflow_datetime(dt) == "2013-12-05T00:00:00"


def test_format_wflow_datetime_non_midnight():
    dt = datetime(2013, 12, 5, 6, 30, 0)
    assert format_wflow_datetime(dt) == "2013-12-05T06:30:00"

"""
Shared utilities for the COMPASS pipeline build scripts.

These functions are used across model_building/wflow/ and model_building/sfincs/
to inject runtime values into HydroMT v1 steps lists loaded from YAML config files.
"""


def find_steps(steps: list, key: str) -> list:
    """Return the argument dicts of every step whose single key == `key`.

    Each step in `steps` is a single-key dict like {"component.method": {args...}}.
    Returns a list of the *mutable* arg dicts so callers can modify them in place
    to inject runtime values (region bbox, river_upa, etc.) before calling build/update.

    Returns an empty list if no step matches `key`.
    """
    return [list(s.values())[0] for s in steps if list(s.keys())[0] == key]


def find_step_index(steps: list, key: str, default: int | None = None) -> int:
    """Return the index of the first step whose key == `key`.

    If not found, returns `default` when given, or len(steps) otherwise.
    Useful for inserting a new step before an existing one via steps.insert(idx, new_step).
    """
    for i, s in enumerate(steps):
        if list(s.keys())[0] == key:
            return i
    return default if default is not None else len(steps)

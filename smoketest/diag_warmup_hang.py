"""Pinpoint the update_forcing_wflow_warmup hang.

Replicates exactly what the snakemake rule does, but:
  - forces the dask *synchronous* scheduler, so there is no thread pool that can deadlock;
    if the hang disappears here, the threaded scheduler is the culprit
  - arms faulthandler so a hang prints the live Python stack instead of sitting silent
  - runs each step separately, so we learn which one blocks

Usage:  python diag_warmup_hang.py [threaded]
        (default is synchronous; pass "threaded" to reproduce the workflow's scheduler)
"""
import faulthandler
import sys
import time
from datetime import datetime, timedelta
from os.path import join

import dask

MODE = sys.argv[1] if len(sys.argv) > 1 else "synchronous"
SCHED = "threads" if MODE == "threaded" else "synchronous"

# Dump every live thread's stack after 120 s and keep repeating, so a block is visible.
faulthandler.dump_traceback_later(120, repeat=True, exit=False)

WFLOW_ROOT = "/p/11210471-001-compass/02_Models/somerset/SomersetLevels_v1int/wflow_vito"
CATS = [
    "/u/morenodu/git_repos/compass-v1-integrate/Workflows/03_data_catalogs/datacatalog_general_v1___linux.yml",
    "/u/morenodu/git_repos/compass-v1-integrate/Workflows/03_data_catalogs/datacatalog_CF_forcing_v1___linux.yml",
]
import os
# Target matters: the rule writes to /p (NFS). Pass "local" to write to /tmp instead.
OUT = ("/tmp/claude-1028338/-u-morenodu-git-repos-compound-flooding-tropical-cyclones/978daef8-ed38-49cd-ac5c-1dea811b52d0/scratchpad/diag_warmup_out"
       if os.environ.get("DIAG_TARGET") == "local"
       else "/p/11210471-001-compass/03_Runs/somerset/DIAG_warmup/wflow_vito/warmup")

# Same window the crash test used: 20-day warmup ending 2 days before the sfincs start.
start_time = "20140101 000000"
t0 = datetime.strptime(start_time, "%Y%m%d %H%M%S") - timedelta(days=2)
start_warmup = datetime.strftime(t0 - timedelta(days=20), "%Y-%m-%dT%H:%M:%S")
end_warmup = datetime.strftime(t0, "%Y-%m-%dT%H:%M:%S")


def step(label, fn):
    print(f"\n--- {label} ---", flush=True)
    t = time.time()
    try:
        fn()
        print(f"    OK in {time.time()-t:.1f}s", flush=True)
    except Exception as e:
        print(f"    FAILED after {time.time()-t:.1f}s: {type(e).__name__}: {str(e)[:200]}", flush=True)
        raise


def main():
    from hydromt_wflow import WflowSbmModel

    print(f"dask scheduler = {SCHED}", flush=True)
    print(f"warmup window  = {start_warmup} -> {end_warmup}", flush=True)

    with dask.config.set(scheduler=SCHED):
        mod = WflowSbmModel(root=WFLOW_ROOT, data_libs=CATS, mode="r")
        step("read base model", mod.read)

        # Relocate before the setup steps so writes land in the scratch dir, not on /p.
        # hydromt 1.x: model.root is a ModelRoot object with .set(); set_root() is gone.
        step("root.set -> scratch", lambda: mod.root.set(OUT, mode="w+"))

        step("setup_config", lambda: mod.setup_config(data={
            "time.starttime": start_warmup,
            "time.endtime": end_warmup,
            "time.timestepsecs": 86400,
            "model.cold_start__flag": True,
            "state.path_output": join("..", "..", "events", "instate", "instates.nc"),
            "input.path_static": join("..", "staticmaps.nc"),
            "input.path_forcing": "inmaps.nc",
        }))

        # These only build a lazy graph; the cost lands in the write below.
        step("setup_precip_forcing (lazy)", lambda: mod.setup_precip_forcing(
            precip_fn="era5_daily", precip_clim_fn=None, chunksize=10))

        step("setup_temp_pet_forcing (lazy)", lambda: mod.setup_temp_pet_forcing(
            temp_pet_fn="era5_daily",
            press_correction=True,
            temp_correction=True,
            dem_forcing_fn="era5_orography",
            pet_method="debruin",
            skip_pet=False,
            chunksize=10,
        ))

        # Faithful replication of the rule: one mod.update() that writes the FULL model
        # (forcing + staticmaps + config) to the target root.
        step("mod.update(write=True) -> full model write", lambda: mod.update(
            model_out=OUT, write=True, forceful_overwrite=True, steps=[]))


if __name__ == "__main__":
    main()

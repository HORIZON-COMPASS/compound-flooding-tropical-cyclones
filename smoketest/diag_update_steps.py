"""Replicate the rule EXACTLY: steps passed inside mod.update(), not pre-applied."""
import faulthandler, sys, time
from datetime import datetime, timedelta
from os.path import join
import dask
from hydromt_wflow import WflowSbmModel

SCHED = "threads" if (len(sys.argv) > 1 and sys.argv[1] == "threaded") else "synchronous"
faulthandler.dump_traceback_later(90, repeat=True, exit=False)

WFLOW_ROOT = "/p/11210471-001-compass/02_Models/somerset/SomersetLevels_v1int/wflow_vito"
OUT = "/p/11210471-001-compass/03_Runs/somerset/DIAG_steps/wflow_vito/warmup"
CATS = ["/u/morenodu/git_repos/compass-v1-integrate/Workflows/03_data_catalogs/datacatalog_general_v1___linux.yml",
        "/u/morenodu/git_repos/compass-v1-integrate/Workflows/03_data_catalogs/datacatalog_CF_forcing_v1___linux.yml"]

t0 = datetime.strptime("20140101 000000", "%Y%m%d %H%M%S") - timedelta(days=2)
sw = datetime.strftime(t0 - timedelta(days=20), "%Y-%m-%dT%H:%M:%S")
ew = datetime.strftime(t0, "%Y-%m-%dT%H:%M:%S")

steps = [
    {"setup_config": {"data": {
        "time.starttime": sw, "time.endtime": ew, "time.timestepsecs": 86400,
        "model.cold_start__flag": True,
        "state.path_output": join("..", "..", "events", "instate", "instates.nc"),
        "input.path_static": join("..", "staticmaps.nc"),
        "input.path_forcing": "inmaps.nc"}}},
    {"setup_precip_forcing": {"precip_fn": "era5_daily", "precip_clim_fn": None, "chunksize": 10}},
    {"setup_temp_pet_forcing": {"temp_pet_fn": "era5_daily", "press_correction": True,
        "temp_correction": True, "dem_forcing_fn": "era5_orography",
        "pet_method": "debruin", "skip_pet": False, "chunksize": 10}},
]

print(f"scheduler={SCHED}  window={sw} -> {ew}", flush=True)
with dask.config.set(scheduler=SCHED):
    mod = WflowSbmModel(root=WFLOW_ROOT, data_libs=CATS, mode="r")
    mod.read(); print("  read OK", flush=True)
    t = time.time()
    mod.update(model_out=OUT, steps=steps, write=True, forceful_overwrite=True)
    print(f"  update(steps=...) OK in {time.time()-t:.1f}s", flush=True)

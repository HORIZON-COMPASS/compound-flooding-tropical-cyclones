from pathlib import Path
from polytope.api import Client
import xarray as xr

SCRATCH = Path("/scratch-shared/dvertegaal1")

OUTDIR = SCRATCH / "data" / "ECMWF_ClimateDT"
OUTDIR.mkdir(parents=True, exist_ok=True)

TMPDIR = Path("/scratch-shared/dvertegaal1/tmp")
TMPDIR.mkdir(exist_ok=True)

client = Client(address="polytope.mn5.apps.dte.destination-earth.eu")

# Params at pressue levels (pl)
params_pl = [60, 129, 131, 132] # potential vorticity, geopotential, U component wind, V component wind
varnames_pl = {60:'pv', 129:'z', 131:'u', 132:'v'}
experiments = ['hist']
realizations = [1]
dates = ['20190303/to/20190318']

request = {
    "activity": "story-nudging",
    "class": "d1",
    "dataset": "climate-dt",
    "generation": "2",
    "levtype": "pl",
    "levelist": "1000/925/850/700/500/300/200/100",
    "model": "ifs-fesom",
    "expver": "0001",
    "resolution": "high",
    "stream": "clte",
    "type": "fc",
    "grid": "0.1/0.1",
    "area": "-12/30/-22/45",
    "param": [str(p) for p in params_pl],
    "experiment": experiments[0],
    "realization": str(realizations[0]),
    "date": dates[0],
    "time": "/".join([f"{h:02d}00" for h in range(24)]),
}

date_str = dates[0].replace("/to/", "_to_")

fname = (
    f"climateDT_"
    f"{'_'.join(varnames_pl[p] for p in params_pl)}_"
    f"{experiments[0]}_"
    f"r{realizations[0]}_"
    f"{date_str}"
)

tempfile = TMPDIR / f"{fname}.grib"
outfile = OUTDIR / f"{fname}.nc"

print("Downloading to:", tempfile)

client.retrieve("destination-earth", request, str(tempfile))

print("Download finished")

try:
    ds = xr.open_dataset(tempfile, engine="cfgrib", backend_kwargs={"indexpath": ""})
    print("Opened GRIB file:", tempfile)
    print("GRIB file:", ds)

    encoding = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
    ds.to_netcdf(outfile, encoding=encoding)
finally:
    tempfile.unlink(missing_ok=True)

print("Saved:", outfile)
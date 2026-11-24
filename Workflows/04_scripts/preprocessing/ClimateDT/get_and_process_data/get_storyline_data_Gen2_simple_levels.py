'''
Script to access and view the storyline data from Climate DT.

Author: Doris Vertegaal 

Based on examples from: polytope-examples GitHub page

To run this script you need to install a Python environment with (pip):
polytope-client

Each request downloads one variable and pressure level for the configured
date range. NetCDF files are stored under variable/scenario/realization;
existing files are skipped so an interrupted run can be restarted.
Pass a scenario (hist, cont, or Tplus2.0K) to download only that scenario;
without an argument, all scenarios are downloaded.

You also need an access token in the home directory, file ~/.polytopeapirc
This is created by running the desp-authentication.py script from https://github.com/destination-earth-digital-twins/polytope-examples/

Available period: 2017 to today 

Guide to scenarios:
- "hist" - observed event
- "cont" - control scenario (mid-20th century)
- "Tplus2.0K" - 2 degree warming, SSP3.7.0

More info: https://destine.ecmwf.int/news/the-fast-development-of-destines-climate-change-adaptation-digital-twin/
and https://platform.destine.eu/services/documents-and-api/doc/?service_name=climate-dt-user-guide
'''
#%%
import argparse
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory

import xarray as xr
from polytope.api import Client

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
PATH = Path("/projects/prjs2226")
OUTDIR = PATH / "data" / "ClimateDT" / "pl" / "raw"
TMPDIR = PATH / "data" / "ClimateDT" / "pl" / "tmp"

params_pl = {60: "pv", 129: "z", 131: "u", 132: "v", 130: "t", 133: "q", 135: "w", 157: "r", 246: "clwc"}
levels = [1000, 925, 850, 700, 600, 500, 400, 300, 200, 100, 50, 1]

experiments = ["hist", "cont", "Tplus2.0K"]
realizations = ["1", "2", "3", "4", "5"]
date_range = "20190301/to/20190322"
date_string = date_range.replace("/to/", "_to_")
area = "-12/30/-22/45"

# ------------------------------------------------------------------
# Base Polytope request
# ------------------------------------------------------------------
base_request = {
    "activity": "story-nudging",
    "class": "d1",
    "dataset": "climate-dt",
    "generation": "2",
    "levtype": "pl",
    "model": "ifs-fesom",
    "expver": "0001",
    "resolution": "high",
    "stream": "clte",
    "type": "fc",
    "grid": "0.1/0.1",
    "area": area,
    "date": date_range,
    "time": "/".join(f"{hour:02d}00" for hour in range(24)),
}


def main(selected_experiments=experiments):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    TMPDIR.mkdir(parents=True, exist_ok=True)
    client = Client(address="polytope.mn5.apps.dte.destination-earth.eu")

    for experiment in selected_experiments:
        for realization in realizations:
            for param, variable_name in params_pl.items():
                variable_outdir = OUTDIR / variable_name / experiment / f"r{realization}"
                variable_outdir.mkdir(parents=True, exist_ok=True)

                for level in levels:
                    output_file = variable_outdir / (
                        f"climateDT_{variable_name}_{level}_{experiment}_"
                        f"r{realization}_{date_string}.nc"
                    )
                    if output_file.exists():
                        print(f"Already exists, skipping: {output_file}", flush=True)
                        continue

                    request = base_request.copy()
                    request.update({
                        "param": str(param),
                        "levelist": str(level),
                        "experiment": experiment,
                        "realization": realization,
                    })

                    print(
                        f"Downloading {variable_name} at {level} hPa, "
                        f"{experiment}, r{realization}", flush=True,
                    )
                    with TemporaryDirectory(dir=TMPDIR) as tmpdir:
                        grib_file = Path(tmpdir) / "download.grib"
                        client.retrieve("destination-earth", request, str(grib_file))

                        with xr.open_dataset(
                            grib_file, engine="cfgrib", backend_kwargs={"indexpath": ""}
                        ) as ds:
                            ds = ds.load()

                        if "isobaricInhPa" not in ds.dims:
                            ds = ds.expand_dims(isobaricInhPa=[level])
                        ds.attrs.update({
                            "experiment": experiment,
                            "realization": realization,
                            "date_range": date_range,
                            "area": area,
                            "source": "Destination Earth Climate DT",
                        })
                        encoding = {
                            name: {"zlib": True, "complevel": 4}
                            for name in ds.data_vars
                        }

                        # Write next to the destination so replacement is atomic.
                        with NamedTemporaryFile(
                            dir=variable_outdir, suffix=".nc", delete=False
                        ) as tmp:
                            tmp_output = Path(tmp.name)
                        try:
                            ds.to_netcdf(tmp_output, encoding=encoding)
                            tmp_output.replace(output_file)
                        finally:
                            tmp_output.unlink(missing_ok=True)

                    print(f"Saved: {output_file}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download Climate DT pressure-level data.")
    parser.add_argument(
        "experiment", nargs="?", choices=experiments,
        help="Scenario to download (default: all scenarios).",
    )
    args = parser.parse_args()
    main([args.experiment] if args.experiment else experiments)

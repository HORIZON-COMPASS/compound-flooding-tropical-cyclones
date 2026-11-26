'''
Script to access and view the storyline data from Climate DT.

Author: Doris Vertegaal 

Based on examples from: polytope-examples GitHub page

To run this script you need to install a Python environment with (pip):
polytope-client

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
from pathlib import Path
import xarray as xr
from polytope.api import Client
from calendar import monthrange
import sys

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
PATH = Path("/projects/prjs2226")
OUTDIR = PATH / "data" / "ClimateDT" / "sfc" / "raw"
TMPDIR = PATH / "data" / "ClimateDT" / "sfc" / "tmp"

OUTDIR.mkdir(parents=True, exist_ok=True)
TMPDIR.mkdir(parents=True, exist_ok=True)

client = Client(address="polytope.mn5.apps.dte.destination-earth.eu")

params_sfc = {235055: "avg_tprate", 151: "msl", 165: "10u", 166: "10v"} # total precip, msl, u10, v10

# experiments = ["hist", "cont", "Tplus2.0K"]
# realizations = ["1", "2", "3", "4", "5"]
experiment = sys.argv[1]
realizations = ["1", "2", "3", "4", "5"]

start_year, start_month = 2017, 1
end_year, end_month = 2026, 7
times = "/".join(f"{hour:02d}00" for hour in range(0, 24, 3))

area = "0/10/-40/135"

def monthly_periods(start_year, start_month, end_year, end_month):
    year, month = start_year, start_month

    while (year, month) <= (end_year, end_month):
        last_day = monthrange(year, month)[1]

        start_date = f"{year:04d}{month:02d}01"
        end_date = f"{year:04d}{month:02d}{last_day:02d}"

        yield year, month, start_date, end_date

        if month == 12:
            year += 1
            month = 1
        else:
            month += 1

# date_string = date_range.replace("/to/", "_to_")


# ------------------------------------------------------------------
# Base Polytope request
# ------------------------------------------------------------------
base_request = {
    "activity": "story-nudging",
    "class": "d1",
    "dataset": "climate-dt",
    "generation": "2",
    "levtype": "sfc", # pl
    "model": "ifs-fesom",
    "expver": "0001",
    "resolution": "high",
    "stream": "clte",
    "type": "fc",
    "grid": "0.1/0.1",
    "area": area,
    "time": times,
}


# ------------------------------------------------------------------
# Download one variable/scenario/realization/month at a time
# ------------------------------------------------------------------
for realization in realizations:
    for param, variable_name in params_sfc.items():
        variable_outdir = (OUTDIR / variable_name / experiment / f"r{realization}")

        variable_outdir.mkdir(parents=True, exist_ok=True)

        for year, month, start_date, end_date in monthly_periods(start_year, start_month, end_year, end_month):
            date_range = f"{start_date}/to/{end_date}"

            output_file = variable_outdir / (f"climateDT_{variable_name}_{experiment}_" 
                                            f"r{realization}_{year:04d}{month:02d}.nc")

            grib_file = TMPDIR / (f"climateDT_{variable_name}_{experiment}_"
                                f"r{realization}_{year:04d}{month:02d}.grib")

            # Allows safely restarting the script.
            if output_file.exists():
                print(f"Already exists, skipping: {output_file}")
                continue

            request = base_request.copy()
            request.update({"param": str(param), "experiment": experiment, "realization": realization, "date": date_range})

            try:
                print(f"Downloading {variable_name}, {experiment}, "f"r{realization}, {year:04d}-{month:02d}", flush=True)

                client.retrieve("destination-earth", request, str(grib_file))

                # Open only the current monthly file.
                with xr.open_dataset(grib_file, engine="cfgrib", backend_kwargs={"indexpath": ""}) as ds:
                    ds = ds.load()

                ds.attrs.update({
                        "experiment": experiment,
                        "realization": realization,
                        "date_range": date_range,
                        "area": area,
                        "source": "Destination Earth Climate DT",
                    })

                encoding = {name: {"zlib": True, "complevel": 4, "shuffle": True} for name in ds.data_vars}
                ds.to_netcdf(output_file, encoding=encoding)

                print(f"Saved: {output_file}", flush=True)

            except Exception as e:
                print(f"FAILED: {output_file}")
                print(e)
                output_file.unlink(missing_ok=True)
                continue

            finally:
                grib_file.unlink(missing_ok=True)

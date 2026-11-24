from pathlib import Path
import xarray as xr
from dask.diagnostics import ProgressBar

PATH = Path("/projects/prjs2226")
DATA_DIR = PATH / "data" / "ClimateDT" / "sfc" / "raw" / "msl"
OUT_DIR = PATH / "data" / "ClimateDT" / "sfc" / "processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)

experiments = ["hist", 
            #    "cont", "Tplus2.0K"
               ]
realizations = ["1",
                #  "2", "3", "4", "5"
                 ]

variables = ["msl"]
# rename_vars = {"avg_tprate": "tp", "msl": "msl", "u10": "u10",  "v10": "v10"}
standard_names = {'msl': 'mean_sea_level_pressure'}
units = {'msl': 'Pa'}

expected_months = (2026 - 2017) * 12 + 7

for experiment in experiments:
    for realization in realizations:

        outfile = (OUT_DIR / f"climateDT_msl_{experiment}_r{realization}_20170101_to_20260731.nc")

        if outfile.exists():
            print(f"Skipping {outfile}")
            continue

        files = sorted((DATA_DIR / experiment / f"r{realization}").glob(f"climateDT_msl_{experiment}_r{realization}_*.nc")) 
        

        print(f"{experiment} r{realization}: "
              f"{len(files)} monthly MSL files")

        if len(files) == 0:
            raise FileNotFoundError(f"No files found for {experiment} r{realization}")

        if len(files) != expected_months:
            print(f"WARNING: found {len(files)} files, "
                  f"expected {expected_months}")

        ds = xr.open_mfdataset(files, combine="by_coords", parallel=False, data_vars="minimal", coords="minimal", compat="override").sortby("latitude")
        print("Dataset opened")
        print(ds)

        ds.attrs.update({"experiment": experiment, "realization": realization, "area": "0/10/-40/135", "source": "Destination Earth ClimateDT"})

        encoding = {"msl": {"dtype": "float32", "zlib": True, "complevel": 4, "shuffle": True}, "latitude": {"dtype": "float32"}, "longitude": {"dtype": "float32"}}

        print(f"Writing {outfile}")

        write_job = ds.to_netcdf(outfile, engine="netcdf4", format="NETCDF4", unlimited_dims=["time"], encoding=encoding, compute=False)

        with ProgressBar():
            write_job.compute()

        ds.close()

        print(f"Saved {outfile}")
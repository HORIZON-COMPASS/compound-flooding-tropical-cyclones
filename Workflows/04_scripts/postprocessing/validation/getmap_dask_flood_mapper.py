from importlib.resources import files

#import hvplot.xarray  # noqa
import dask
import xarray as xr
from dask.distributed import Client
from dask_flood_mapper import flood


client = Client(processes=False, threads_per_worker=2, n_workers=1, memory_limit="15GB")

time_range = "2019-03-09/2019-03-27"
bbox = [34.1, -20.3, 35.2, -19.05]

fd = flood.decision(bbox=bbox, datetime=time_range).compute()
fd

client.close()
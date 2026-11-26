# Based on script from Marjolein Ribberink: 
# https://github.com/MRibberink/OpheliaPaperCode/blob/main/hurricane_funcs.py
#%%

def import_ibtracs_period(start="2017-01-01", end="2026-07-01"):
    import xarray as xr
    import pandas as pd
    from os.path import join
    import numpy as np
    import platform
    from pathlib import Path

    if platform.system() == "Windows":
        path = Path("P:/11210471-001-compass/01_Data/IBTrACS/")
    else: # Linux (e.g. Snellius)
        path = Path("/projects/prjs2226/data/IBTrACS/")

    tracks = xr.open_dataset(join(path, "IBTrACS.SI.v04r01.nc"))

    storms = {}
    for i in range(tracks.sizes["storm"]):

        storm = tracks.isel(storm=i)
        time = pd.to_datetime(storm.time.values)
        valid = ((time >= pd.Timestamp(start)) & (time < pd.Timestamp(end)) & ~pd.isna(time))

        if valid.sum() == 0:
            continue

        name = np.asarray(storm.name.values).item()

        if isinstance(name, (bytes, np.bytes_)):
            name = name.decode("utf-8").strip()
        else:
            name = str(name).strip()

        tc = {
            "time": time[valid],
            "lat": storm.reunion_lat.values[valid],
            "lon": storm.reunion_lon.values[valid],
            "min_pres": storm.reunion_pres.values[valid],
            "wind": storm.reunion_wind.values[valid],
            "eye": storm.usa_eye.values[valid],
            "nature": storm.nature.values[valid],
            "usa_sshs": storm.usa_sshs.values[valid],
            "rmw": storm.reunion_rmw.values[valid],
        }

        storms[name] = pd.DataFrame(tc)

    print(f"Found {len(storms)} storms")

    return storms


# %%

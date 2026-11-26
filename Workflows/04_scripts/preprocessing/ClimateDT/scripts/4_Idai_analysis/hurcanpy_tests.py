# https://github.com/Huracan-project/huracanpy

#%%
import huracanpy
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from os.path import join
import platform

prefix = "p:/" if platform.system() == "Windows" else "/p/"
path = ("C:/Code/COMPASS/Workflows/04_scripts/preprocessing/ClimateDT/data/")
tracks=xr.open_dataset(join(path, 'IBTrACS.SI.v04r01.nc'))
track_year=tracks.where(tracks.season==2019,drop=True)
data=track_year.where(track_year.name==b'IDAI',drop=True)
tc={'time':data.time.values[data.time.values==data.time.values],
    'lat':data.reunion_lat.values[data.time.values==data.time.values],
    'lon':data.reunion_lon.values[data.time.values==data.time.values],
    'min_pres':data.reunion_pres.values[data.time.values==data.time.values],
    'wind':data.reunion_wind.values[data.time.values==data.time.values],
    'eye':data.usa_eye.values[data.time.values==data.time.values],
    'nature':data.nature.values[data.time.values==data.time.values],
    'rmw':data.reunion_rmw.values[data.time.values==data.time.values]
    }
tc_idai=pd.DataFrame(tc)
tc_idai['track_id'] = data.sid.values.item().decode()
tc_idai["time"] = pd.to_datetime(tc_idai["time"].dt.floor("h"))


#%%

ib = huracanpy.load(source="ibtracs", ibtracs_subset="wmo")
ib_idai = ib.where(ib.track_id==data.sid.values.item().decode(), drop=True)
ib_idai

# wmo wind are likely 10-min sustained winds
huracanpy.plot.tracks(ib_idai.lon, ib_idai.lat, intensity_var=ib_idai.wind)

#%%
# Needs track_id, lon, lat, time
# # tracks = huracanpy.load(huracanpy.example_csv_file)
# track = huracanpy.load(data)
# print(track)

# Iterating over all tracks
# Each track will be a subset of the xarray Dataset with a unique track_id
# The track_id is not necessarily an integer, it follows whatever you have loaded
# e.g. could be a string for IBTrACS
for track_id, ib_idai in ib_idai.groupby("track_id"):
    # Do something with the track
    print(track_id, len(ib_idai.time))

    

# Quickly view the tracks
ib_idai.hrcn.plot_tracks(intensity_var_name="wind")

# Add a new variable to the tracks and plot this instead
tracks = ib_idai.hrcn.add_is_land()
tracks.hrcn.plot_tracks(intensity_var_name="is_land")

plt.show()
# %%
import xarray as xr
import numpy as np


# inspecitive vertical profile along tracked Idai in hist scenario, realization 1
hist_pl_r1 = xr.open_dataset(r"p:\11210471-001-compass\01_Data\ECMWF_ClimateDT\data\preprocessed\Gen2\Idai_full\climateDT_pv_z_u_v_hist_r1_20190303_to_20190318.nc")
hist_pl_r1

#%%
ds=hist_pl_r1
# ds: ClimateDT pressure-level dataset
# ib_idai: HuracanPy/IBTrACS track

times = []
profiles = []

for i in range(len(ib_idai.record)):
    t = ib_idai.time.values[i]
    lat = ib_idai.lat.values[i]
    lon = ib_idai.lon.values[i]

    prof = (ds.sel(time=t, latitude=lat, longitude=lon, method="nearest"))
    profiles.append(prof)
    times.append(prof.time.values)

track_ds = xr.concat(profiles, dim="time")
track_ds = track_ds.assign_coords(time=times)

track_ds['wind'] = np.sqrt(track_ds['u']**2 + track_ds['v']**2)



#%%
levels = [1000, 850, 700, 500, 200]

fig, ax = plt.subplots()
for lev in levels:
    track_ds["wind"].sel(isobaricInhPa=lev).plot(ax=ax, label=f"{lev} hPa")

ax.legend()


# %%
track_ds["pv"].plot.contour(x="time", y="isobaricInhPa", levels=20, yincrease=False)

#%%
track_daily = track_ds.resample(time="1D").mean()

track_daily["wind"] = np.sqrt(track_daily["u"]**2 + track_daily["v"]**2)

fig, ax = plt.subplots(figsize=(8, 5))

track_daily["wind"].plot(x="time", y="isobaricInhPa", yincrease=False, cmap="viridis", ax=ax)

ax.set_ylabel("Pressure (hPa)")
ax.set_title("Daily mean wind speed along Idai track")
plt.show()


# %%
import matplotlib.pyplot as plt
import numpy as np

# Daily averages (optional)
wind_ds = track_ds.resample(time="1D").mean()

# Coordinates
times = wind_ds.time.values
levels = wind_ds.isobaricInhPa.values

# Make 2D coordinate arrays
X, Y = np.meshgrid(times, levels)

fig, ax = plt.subplots(figsize=(10, 6))

q = ax.quiver(X, Y, wind_ds["u"].values, wind_ds["v"].values, 
              np.sqrt(wind_ds["u"].values**2 + wind_ds["v"].values**2), cmap="viridis", scale=500)

ax.invert_yaxis()
ax.set_ylabel("Pressure (hPa)")
ax.set_xlabel("Time")
ax.set_title("Daily mean wind vectors along TC Idai track")

plt.colorbar(q, ax=ax, label="Wind speed (m/s)")
plt.show()

# %%

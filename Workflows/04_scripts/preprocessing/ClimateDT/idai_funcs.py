# idai_funcs.py
# Based on script from Marjolein Ribberink: 
# https://github.com/MRibberink/OpheliaPaperCode/blob/main/hurricane_funcs.py
from ntpath import join
import platform
import os

prefix = "p:/" if platform.system() == "Windows" else "/p/"


def import_ibtracs(storm='IDAI',year=2019):
    #Imports the IBTRACS data and selects out the most important
    import xarray as xr
    import pandas as pd
    from os.path import join
    
    path = join(prefix, "11210471-001-compass/01_Data/IBTrACS/")
    tracks=xr.open_dataset(join(path, 'IBTrACS.SI.v04r01.nc'))
    track_year=tracks.where(tracks.season==year,drop=True)
    data=track_year.where(track_year.name==storm.encode(),drop=True)
    tc={
        'time':data.time.values[data.time.values==data.time.values],
        'lat':data.reunion_lat.values[data.time.values==data.time.values],
        'lon':data.reunion_lon.values[data.time.values==data.time.values],
        'min_pres':data.reunion_pres.values[data.time.values==data.time.values],
        'wind':data.reunion_wind.values[data.time.values==data.time.values],
        'eye':data.usa_eye.values[data.time.values==data.time.values],
        'nature':data.nature.values[data.time.values==data.time.values],
        'rmw':data.reunion_rmw.values[data.time.values==data.time.values]
        }

    new_track=pd.DataFrame(data=tc)
    return(new_track)


def import_data(path=os.path.join(prefix,'11210471-001-compass/01_Data/ECMWF_ClimateDT/data/preprocessed/Gen2/Idai_full')):
    import xarray as xr

    datadir    = path
    ds_hist    = xr.open_mfdataset(os.path.join(datadir, 'climateDT_*_hist_*_Idai_full_20190301_to_20190329.nc'),
                                  combine='by_coords', data_vars='all')
    ds_cont    = xr.open_mfdataset(os.path.join(datadir, 'climateDT_*_cont_*_Idai_full_20190301_to_20190329.nc'),
                                  combine='by_coords', data_vars='all')
    ds_Tplus2k = xr.open_mfdataset(os.path.join(datadir, 'climateDT_*_Tplus2.0K_*_Idai_full_20190301_to_20190329.nc'),
                                  combine='by_coords', data_vars='all')

    ds_hist    = ds_hist.expand_dims(scenario=["hist"])
    ds_cont    = ds_cont.expand_dims(scenario=["cont"])
    ds_Tplus2k = ds_Tplus2k.expand_dims(scenario=["Tplus2.0K"])

    ds_all = xr.concat([ds_hist, ds_cont, ds_Tplus2k], dim="scenario")
    return ds_all


def dist(a, b):
    #Calculates the magnitude of two variables in an xarray-compliant way
    import xarray as xr
    import numpy as np
    func = lambda x, y: np.sqrt(x**2 + y**2)
    return xr.apply_ufunc(func, a, b)


def storm_tracker(data,best_track,model,s_size):
    import numpy as np
    import xarray as xr
    import pandas as pd
    import datetime as dt
    import bisect

    if model=='RACMO' or model=='RV':
        variables=data.mslp
    elif model=='OPER':
        variables=data.msl
    elif model=='ERA5' or model=='ClimateDT':
        data=data.rename({'latitude':'lat','longitude':'lon'})
        variables=data.msl
    elif model=="GFS":
        data=data.rename({'latitude':'lat','longitude':'lon'})
        variables=data.mslet
    else:
        raise Exception(model+" not recognized. Options are: 'fcst','OPER','HRES','RACMO','deg2','ERA5','ClimateDT'.")

    best_track["time"] = pd.to_datetime(best_track["time"].dt.floor("h"))
    enter=max([variables.time[0], best_track.loc[best_track[["lat", "lon", "min_pres"]].notna().all(axis=1), "time"].iloc[0]])
    exit=min([variables.time[-1], best_track.loc[best_track[["lat", "lon", "min_pres"]].notna().all(axis=1), "time"].iloc[-1]])
    timesteps=pd.to_datetime(data.time.values)
    idx_enter=bisect.bisect_left(timesteps, enter) #Index of entrance time
    idx_exit=bisect.bisect(timesteps, exit)        #Index of dissipation time

    mod_time=timesteps[idx_enter:idx_exit]         #Times actually plotted

    # find mslp minima in box around first point. In case of multiple minima, find the closest to the IBTrACS point

    idx_0=bisect.bisect_left(best_track['time'],enter)

    c_frame=variables.sel(time=enter).compute()
    search_area=c_frame.where((c_frame.lat>(float(best_track.lat[idx_0])-s_size)) &
                              (c_frame.lat<(float(best_track.lat[idx_0])+s_size)) &
                              (c_frame.lon>(float(best_track.lon[idx_0])-s_size)) &
                              (c_frame.lon<(float(best_track.lon[idx_0])+s_size)),drop=True).compute()
    #return c_frame
    #relv=c_frame.where(c_frame.vo==search_area.vo.max(), drop=True).to_dataframe().reset_index()

    min_list=c_frame.where(c_frame==search_area.min(), drop=True)
    dist_list=dist(min_list.lat-float(best_track.lat[idx_0]),min_list.lon-float(best_track.lon[idx_0]))
    min_p=c_frame.where(dist_list==dist_list.min(),drop=True).to_dataframe().reset_index()


    s_a=search_area
    #second point is also different since we can't extrapolate yet, also can't use the next track point
    #Instead, we just up the time 1, which works since it isn't moving fast enough yet, and compare with point 1

    c_frame2=variables.sel(time=enter+(data.time[1]-data.time[0])).compute()
    search_area2=c_frame2.where((c_frame2.lat>(float(best_track.lat[idx_0])-s_size)) &
                              (c_frame2.lat<(float(best_track.lat[idx_0])+s_size)) &
                              (c_frame2.lon>(float(best_track.lon[idx_0])-s_size)) &
                              (c_frame2.lon<(float(best_track.lon[idx_0])+s_size)),drop=True).compute()

    #relv=pd.concat([relv,c_frame2.where(c_frame2.vo==search_area2.vo.max(), drop=True).to_dataframe().reset_index()],ignore_index=True)

    min_list2=c_frame2.where(c_frame2==search_area2.min(), drop=True)
    dist_list2=dist(min_list2.lat-float(best_track.lat[idx_0]),min_list2.lon-float(best_track.lon[idx_0]))
    min_p=pd.concat([min_p,c_frame2.where(dist_list2==dist_list2.min(),drop=True).to_dataframe().reset_index()],ignore_index=True)

    #After this we extrapolate based on the previous 2 points
    s_a=xr.concat([s_a,search_area2], dim='time', join="outer")
    for i in np.arange(2,len(mod_time)):

        #print(i)
        c_frame=variables.sel(time=enter+(data.time[i]-data.time[0])).compute()
        guess_lat=min_p['lat'][i-1]
        guess_lon=min_p['lon'][i-1]
        if guess_lat.size>1:
            guess_lat=guess_lat[0]
        if guess_lon.size>1:
            guess_lon=guess_lon[0]

        guess_p=c_frame.sel(lat=guess_lat,method='nearest').sel(lon=guess_lon,method='nearest')
        search_area=c_frame.where((c_frame.lat>(float(guess_p.lat.values))-s_size) &
                                  (c_frame.lat<(float(guess_p.lat.values))+s_size) &
                                  (c_frame.lon>(float(guess_p.lon.values))-s_size) &
                                  (c_frame.lon<(float(guess_p.lon.values))+s_size),drop=True)

        #relv=pd.concat([relv,search_area.where(search_area.vo==search_area.vo.max(), drop=True).to_dataframe().reset_index()],ignore_index=True)

        min_list=c_frame.where(c_frame==search_area.min(), drop=True)
        s_a=xr.concat([s_a, search_area],dim="time", join="outer")
        dist_list=dist(min_list.lat-guess_lat,min_list.lon-guess_lon)
        min_p=pd.concat([min_p,search_area.where(dist_list==dist_list.min(),drop=True).to_dataframe().reset_index()],ignore_index=True)
    print("done")
    return min_p



def plot_tracks(track, tc_idai_climatedt, scenarios):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.pyplot as plt
    bbox_idai = [30, 45, -22, -12] 

    fig, axes = plt.subplots(2, 2, figsize=(8, 6), dpi=300, constrained_layout=True,
                             subplot_kw={"projection": ccrs.PlateCarree()}, 
                             sharex=True, sharey=True)

    axes = axes.flatten()
    for i, (ax, scen) in enumerate(zip(axes, scenarios)):
        # background map
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax.add_feature(cfeature.LAND, facecolor="lightgray")
        ax.add_feature(cfeature.OCEAN, facecolor="lightblue")

        # IBTrACS track
        ax.plot(track["lon"], track["lat"], "k-", lw=2, transform=ccrs.PlateCarree(),
                label="IBTrACS")

        # realizations
        for (scenario, realization), min_p in tc_idai_climatedt.items():
            if scenario == scen:
                ax.plot(min_p["lon"], min_p["lat"], alpha=0.5, lw=1, transform=ccrs.PlateCarree(),
                        label=f"r{realization}")

        ax.set_extent(bbox_idai, crs=ccrs.PlateCarree())
        ax.set_title(scen)

        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        # Left labels only on first column (plots 0 and 2)
        gl.left_labels = (i % 2 == 0)
        gl.bottom_labels = (i != 0)

    # Use unused axes for legend
    if len(scenarios) < len(axes):
        legend_ax = axes[len(scenarios)]
        legend_ax.axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    legend_ax.legend(handles, labels, loc="upper left", frameon=True)

    plt.tight_layout()
    plt.show()



def main(run_tracking=False):
    import os
    import numpy as np
    import pickle
    print('hello!')
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    processed_path = ('data/processed')
    data_path=os.path.join(prefix,'11210471-001-compass/01_Data/ECMWF_ClimateDT/data/preprocessed/Gen2/Idai_full')
    print("Importing IBTrACS data...")
    new_track=import_ibtracs()
    print("Importing ClimateDT data...")
    data_all=import_data(path=data_path)

    if run_tracking:
        print("Running storm tracking...")
        tc_idai_climatedt = {}
        for scen in data_all.scenario.values:
            for realization in data_all.realization.values:
                print(f"Scenario: {scen}, " f"Realization: {realization}")

                data_subset = data_all.sel(scenario=scen, realization=realization)
                min_p = storm_tracker(data_subset, new_track, model="ClimateDT", s_size=1.5)
                tc_idai_climatedt[(scen, realization)] = min_p

        with open(os.path.join(processed_path, "tc_idai_tracked_climatedt.pkl"), "wb") as f:
            pickle.dump(tc_idai_climatedt, f)

    else:
        print("Loading precomputed tracks...")
        with open(os.path.join(processed_path, "tc_idai_tracked_climatedt.pkl"), "rb") as f:
            tc_idai_climatedt = pickle.load(f)

    return new_track, data_all, tc_idai_climatedt
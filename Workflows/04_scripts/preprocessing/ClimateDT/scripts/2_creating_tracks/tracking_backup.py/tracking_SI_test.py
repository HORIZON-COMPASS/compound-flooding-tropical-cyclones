#%%
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

def import_ibtracs_period(start="2019-01-01", end="2019-02-01"):
    import xarray as xr
    import pandas as pd
    from os.path import join
    import numpy as np

    path = join(prefix, "11210471-001-compass/01_Data/IBTrACS/")
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
            "rmw": storm.reunion_rmw.values[valid],
        }

        storms[name] = pd.DataFrame(tc)

    print(f"Found {len(storms)} storms")

    return storms

# def import_data(path=('c:/Code/test/ClimateDT')):
#     import xarray as xr

#     datadir    = path
#     ds_hist    = xr.open_mfdataset(os.path.join(datadir, 'climateDT_msl_hist_r1_2019*.nc'),
#                                   combine='by_coords', data_vars='all')
#     ds_cont    = xr.open_mfdataset(os.path.join(datadir, 'climateDT_*_cont_*_Idai_full_20190301_to_20190329.nc'),
#                                   combine='by_coords', data_vars='all')
#     ds_Tplus2k = xr.open_mfdataset(os.path.join(datadir, 'climateDT_*_Tplus2.0K_*_Idai_full_20190301_to_20190329.nc'),
#                                   combine='by_coords', data_vars='all')

#     ds_hist    = ds_hist.expand_dims(scenario=["hist"])
#     ds_cont    = ds_cont.expand_dims(scenario=["cont"])
#     ds_Tplus2k = ds_Tplus2k.expand_dims(scenario=["Tplus2.0K"])

#     ds_all = xr.concat([ds_hist, ds_cont, ds_Tplus2k], dim="scenario")
#     ds = ds_hist
#     return ds


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

def plot_tracks_all(storms, tc_climatedt, figure_path="figures"):
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    os.makedirs(figure_path, exist_ok=True)

    scenarios = ["hist", "cont", "Tplus2.0K"]
    colors = plt.cm.tab10.colors

    for storm_name, track in storms.items():
        print(f"Plotting {storm_name}")
        model_tracks = {(scen, realization): model_track for (scen, realization, name), model_track 
                        in tc_climatedt.items() if name == storm_name}

        if not model_tracks:
            print(f"No model tracks found for {storm_name}, skipping.")
            continue

        # Shared extent based on IBTrACS and all ClimateDT tracks
        all_lon = [track["lon"].to_numpy()]
        all_lat = [track["lat"].to_numpy()]

        for model_track in model_tracks.values():
            all_lon.append(model_track["lon"].to_numpy())
            all_lat.append(model_track["lat"].to_numpy())

        all_lon = np.concatenate(all_lon)
        all_lat = np.concatenate(all_lat)
        margin = 3

        extent = [np.nanmin(all_lon) - margin,
                  np.nanmax(all_lon) + margin,
                  np.nanmin(all_lat) - margin,
                  np.nanmax(all_lat) + margin]

        fig, axes = plt.subplots(1, 3, figsize=(15, 5), subplot_kw={"projection": ccrs.PlateCarree()}, constrained_layout=True)

        for ax, scen in zip(axes, scenarios):
            ax.coastlines()
            ax.add_feature(cfeature.BORDERS, linewidth=0.5)
            ax.add_feature(cfeature.LAND, facecolor="lightgray")
            ax.add_feature(cfeature.OCEAN, facecolor="lightblue")

            ax.plot(track["lon"], track["lat"], "k-", lw=2.5, transform=ccrs.PlateCarree(), label="IBTrACS")

            for i, realization in enumerate(["1", "2", "3", "4", "5"]):
                key = (scen, realization)

                if key not in model_tracks:
                    continue

                model_track = model_tracks[key]

                ax.plot(model_track["lon"], model_track["lat"], color=colors[i], lw=1.3, alpha=0.8, 
                        transform=ccrs.PlateCarree(), label=f"r{realization}")

            ax.set_extent(extent, crs=ccrs.PlateCarree())
            ax.set_title(scen)

            gl = ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5)
            gl.top_labels = False
            gl.right_labels = False

        fig.suptitle(storm_name, fontsize=14)

        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=6, frameon=True)

        safe_name = storm_name.replace(" ", "_").replace("/", "_")

        outfile = os.path.join(figure_path, f"{safe_name}_tracked.png")

        fig.savefig(outfile, bbox_inches="tight")

        print(f"Saved {outfile}")

        plt.show()
        plt.close(fig)


def main_backup(run_tracking=False):
    import os
    import numpy as np
    import pickle
    import glob
    import xarray as xr
    print('hello!')
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    processed_path = ('processed/')
    data_path=('c:/Code/test/ClimateDT')

    print("Importing IBTrACS data...")
    start_time = "2019-01-01"
    end_time = "2019-03-01"
    storms = import_ibtracs_period(start=start_time, end=end_time)

    print("Importing ClimateDT data...")
    experiments = ["hist", "cont", "Tplus2.0K"]
    realizations = ["1", "2", "3", "4", "5"]

    for scen in experiments:
        for realization in realizations:
            print(f"Processing scenario: {scen}, realization: {realization}")

            files = sorted(glob.glob(os.path.join(data_path, f"climateDT_msl_{scen}_r{realization}_*.nc")))
            data = xr.open_mfdataset(files, combine="by_coords", data_vars="all")

            if run_tracking:
                print("Running storm tracking...")
                # scen = 'hist'
                # realization = 1
                tc_climatedt = {}

                for storm_name, track in storms.items():
                    print(f"Tracking {storm_name}")

                    if len(track) < 2:
                        print(f"Skipping {storm_name}: too few track points")
                        continue

                    try:
                        min_p = storm_tracker(data, track, model="ClimateDT", s_size=1.5)
                        tc_climatedt[scen, realization, storm_name] = min_p
                    except Exception as e:
                        print(f"Failed for {scen} {realization} {storm_name}: {e}")

                with open(os.path.join(data_path, processed_path, f"storms_tracked_climatedt_{scen}_{realization}_{start_time}_{end_time}.pkl"), "wb") as f:
                    pickle.dump(tc_climatedt, f)

            else:
                print("Loading precomputed tracks...")
                with open(os.path.join(data_path, processed_path, f"storms_tracked_climatedt_{scen}_{realization}_{start_time}_{end_time}.pkl"), "rb") as f:
                    tc_climatedt = pickle.load(f)

        return storms, data, tc_climatedt



def main(run_tracking=False):
    import os
    import glob
    import pickle
    import pandas as pd
    import xarray as xr

    print("hello!")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    data_path = "c:/Code/test/ClimateDT"
    processed_path = os.path.join(data_path, "processed")

    start_time = "2017-01-01"
    end_time   = "2026-08-01"

    print("Importing IBTrACS data...")
    storms = import_ibtracs_period(start=start_time, end=end_time)

    track_file = os.path.join(processed_path, f"storms_tracked_climatedt_{start_time}_{end_time}.pkl")

    table_file = os.path.join(processed_path, f"storms_tracked_climatedt_{start_time}_{end_time}.parquet")

    if run_tracking:
        experiments = ["hist", "cont", "Tplus2.0K"]
        realizations = ["1",
                        #  "2", "3", "4", "5"
                         ]

        tc_climatedt = {}
        all_tracks = []

        for scen in experiments:
            for realization in realizations:
                print(f"\nProcessing {scen} r{realization}")

                files = sorted(glob.glob(os.path.join(data_path, f"climateDT_msl_{scen}_r{realization}_*.nc")))

                if len(files) == 0:
                    print(f"No files found for {scen} r{realization}")
                    continue

                data = xr.open_mfdataset(files, combine="by_coords", data_vars="all")

                for storm_name, track in storms.items():
                    print(f"Tracking {storm_name}")
                    valid = (track[["lat", "lon", "min_pres"]].notna().all(axis=1))

                    if valid.sum() < 2:
                        continue

                    try:
                        min_p = storm_tracker(data, track.copy(), model="ClimateDT", s_size=1.5)
                        tc_climatedt[(scen, realization, storm_name)] = min_p

                        tmp = min_p.copy()
                        tmp["scenario"] = scen
                        tmp["realization"] = realization
                        tmp["storm_name"] = storm_name

                        all_tracks.append(tmp)

                    except Exception as e:
                        print(f"Failed: {scen} r{realization} {storm_name}")

                data.close()

        with open(track_file, "wb") as f:
            pickle.dump(tc_climatedt, f, protocol=pickle.HIGHEST_PROTOCOL)

        if all_tracks:
            tracks_df = pd.concat(all_tracks, ignore_index=True)
            tracks_df.to_parquet(table_file, index=False)
            print(f"Saved {len(tracks_df)} track points")

    else:
        with open(track_file, "rb") as f:
            tc_climatedt = pickle.load(f)

        tracks_df = pd.read_parquet(table_file)

    return storms, tc_climatedt


# track = import_ibtracs(storm="IDAI", year=2019)

if __name__=='__main__':
    storms, tc_climatedt = main(run_tracking=False)

plot_tracks_all(storms, tc_climatedt)

# %%

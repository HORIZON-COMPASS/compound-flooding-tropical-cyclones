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
    import pandas as pd
    import platform

    os.makedirs(figure_path, exist_ok=True)

    scenarios = ["cont", "hist", "Tplus2.0K"]
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

            # Start of track
            ax.plot(track.lon.iloc[0], track.lat.iloc[0], marker="o", color="red", markersize=3,
                    transform=ccrs.PlateCarree(), zorder=10, label="Track start")

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

        max_wind = track.wind.max() * 0.514444  # convert from knots to m/s

        # TC category at maximum wind
        idx_max = track.wind.idxmax()
        tc_strength = track.usa_sshs.loc[idx_max]
        start_date = pd.to_datetime(track.time.min()).strftime("%Y-%m-%d")

        # Max intensity point
        axes[1].scatter(track.lon.loc[idx_max], track.lat.loc[idx_max], color="blue", s=30,
                        marker="o", edgecolor="black", linewidth=0.5, transform=ccrs.PlateCarree(),
                        zorder=5, label=f"Max intensity (SSHS={tc_strength})")

        fig.suptitle(f"{storm_name} | {start_date} | Max wind: {max_wind:.1f} m/s | {tc_strength}")

        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=6, frameon=True)

        safe_name = storm_name.replace(" ", "_").replace("/", "_")

        outfile = os.path.join(figure_path, f"{safe_name}_tracked.png")

        fig.savefig(outfile, bbox_inches="tight")

        print(f"Saved {outfile}")
        if platform.system() == "Windows":
            plt.show()

        plt.close(fig)


def plot_ibtracs_context(storms, tc_climatedt):
    import platform
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import pandas as pd

    # Only storms that have ClimateDT tracks
    model_storms = {key[-1] for key in tc_climatedt.keys()}

    for storm_name in model_storms:
        fig = plt.figure(figsize=(10, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())

        ax.set_extent([10, 130, -40, 0])
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)

        # Plot all IBTrACS tracks in grey
        for name, track in storms.items():
            track = track.dropna(subset=["lon", "lat"])

            if len(track) == 0:
                continue

            ax.plot(track.lon, track.lat, color="#CCCCCC", linewidth=1.5,
                    transform=ccrs.PlateCarree(), zorder=1)

        # Highlight selected storm
        if storm_name in storms:
            track = storms[storm_name].dropna(subset=["lon", "lat"]) 

            ax.plot(track.lon, track.lat, color="red", linewidth=2.5,
                    transform=ccrs.PlateCarree(), zorder=3)

            # Annotation
            ax.annotate(storm_name, (track.lon.iloc[0], track.lat.iloc[0]),
                        xytext=(5, 10), textcoords="offset points", fontsize=12,
                        color="red", fontweight="bold")

        ax.set_title(f"TC {storm_name} within South Indian Ocean Basin (IBTrACS tracks)")

        fig.tight_layout()
        if platform.system() == "Windows":
            plt.show()
        fig.savefig(f"ibtracs_context_{storm_name}.png", dpi=300)
        plt.close(fig)

        
def process_member_tracking(args):
    import os
    import glob
    import pickle
    import time
    import pandas as pd
    import xarray as xr

    scen, realization, storms, data_path, processed_path = args

    pkl_file = os.path.join(processed_path, f"storms_{scen}_r{realization}.pkl")

    if os.path.exists(pkl_file):
        print(f"SKIP {scen} r{realization}: already processed", flush=True)
        return pkl_file
    
    worker_start = time.time()

    print(f"START {scen} r{realization}", flush=True)

    files = sorted(glob.glob(os.path.join(data_path, scen, f"r{realization}", 
                                          f"climateDT_msl_{scen}_r{realization}_*.nc")))

    if not files:
        print(f"No files found: {scen} r{realization}", flush=True)
        return None

    t0 = time.time()

    data = xr.open_mfdataset(files, combine="by_coords", data_vars="minimal", coords="minimal", 
                             compat="override", parallel=False)
    print(f"{scen} r{realization}: opened dataset in {time.time()-t0:.1f}s", flush=True)

    results = {}
    try:
        for storm_name, track in storms.items():
            valid = track[["lat", "lon", "min_pres"]].notna().all(axis=1)

            if valid.sum() < 2:
                continue

            storm_start = time.time()
            try:
                print(f"START {scen} r{realization}: {storm_name}", flush=True)

                results[(scen, realization, storm_name)
                ] = storm_tracker(data, track.copy(), model="ClimateDT", s_size=1.5)

                print(f"DONE {scen} r{realization}: {storm_name} ({time.time()-storm_start:.1f}s)", 
                      flush=True)

            except Exception as error:
                print(f"FAILED {scen} r{realization} {storm_name}: {error}", flush=True)

    finally:
        data.close()
        del data

    pkl_file = os.path.join(processed_path, f"storms_{scen}_r{realization}.pkl")

    with open(pkl_file, "wb") as f:
        pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"Saved {len(results)} storms to {os.path.basename(pkl_file)}", flush=True)

    # free memory before worker exits
    del results

    print(f"FINISHED {scen} r{realization} in {(time.time()-worker_start)/60:.1f} min",
          flush=True)

    return pkl_file


def process_member_wind(args):
    import os
    import glob
    import pickle
    import numpy as np
    import xarray as xr
    import traceback

    scen, realization, input_file, data_path, output_path, radius_km, use_ibtracs_rmw = args

    output_file = os.path.join(output_path,
        f"storms_{scen}_r{realization}_wind_{radius_km}km.pkl")

    if os.path.exists(output_file):
        print(f"SKIP {scen} r{realization}: wind file exists", flush=True)
        return output_file

    print(f"START {scen} r{realization}", flush=True)

    with open(input_file, "rb") as f:
        member_tracks = pickle.load(f)

    u_files = sorted(glob.glob(
                     os.path.join(data_path, "10u", scen, f"r{realization}",
                     f"climateDT_10u_{scen}_r{realization}_*.nc",)))

    v_files = sorted(glob.glob(
                     os.path.join(data_path, "10v", scen, f"r{realization}",
                     f"climateDT_10v_{scen}_r{realization}_*.nc")))

    if not u_files or not v_files:
        print(f"Missing wind files: {scen} r{realization}", flush=True)
        return None

    data = xr.open_mfdataset(u_files + v_files, combine="by_coords", 
                             data_vars="minimal", coords="minimal",
                             compat="override", parallel=False)

    try:
        wind_speed = np.hypot(data["u10"], data["v10"])

        results = {}
        for key, track in member_tracks.items():
            storm_name = key[2] if isinstance(key, tuple) else str(key)

            if use_ibtracs_rmw and storm_name in storms:
                print(f"Using IBTrACS RMW for wind radius for {storm_name}", flush=True)
                rmw = storms[storm_name]["rmw"].max()

                if not np.isnan(rmw):
                    search_radius = rmw * 1.852 + 20
                else:
                    print(f"IBTrACS RMW is NaN for {storm_name}, using default radius of {radius_km} km", flush=True)
                    search_radius = radius_km

            else:
                print(f"Using default radius of {radius_km} km for {storm_name}", flush=True)
                search_radius = radius_km

            try:
                results[key] = add_max_wind_within_radius(track=track,
                                                          wind_speed=wind_speed,
                                                          radius_km=search_radius)

                print(f"DONE {scen} r{realization}: {storm_name}", flush=True)

            except Exception as error:
                print(f"FAILED {scen} r{realization} {storm_name}: {error}",
                      flush=True)

                results[key] = track.copy()

    except Exception:
        print(f"REALIZATION FAILED {scen} r{realization}", flush=True)
        traceback.print_exc()
        raise

    finally:
        data.close()

    with open(output_file, "wb") as f:
        pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)

    size_gb = os.path.getsize(output_file) / 1024**3

    print(f"SAVED {scen} r{realization}: {len(results)} storms, {size_gb:.3f} GB",
          flush=True)

    del wind_speed
    del member_tracks
    del results
    del data

    return output_file


def add_max_wind_within_radius(track, wind_speed, radius_km=200, time_col="time", 
                               lat_col="lat", lon_col="lon"):
    import numpy as np
    import pandas as pd
    import xarray as xr

    wind_col = f"max_wind_model"
    rmw_col = "rmw_model_km"

    track[wind_col] = np.nan
    track[rmw_col] = np.nan
    track["rmw_lat"] = np.nan
    track["rmw_lon"] = np.nan

    r_deg = radius_km / 111.0

    for idx, point in track.iterrows():
        try:
            time = pd.Timestamp(point[time_col])
            lat0 = point[lat_col]
            lon0 = point[lon_col]

            if pd.isna(lat0) or pd.isna(lon0):
                continue
            
            wind_ts = wind_speed.sel(time=time, method="nearest")

            subset = wind_ts.sel(latitude=slice(lat0 + r_deg, lat0 - r_deg),
                                longitude=slice(lon0 - r_deg, lon0 + r_deg))
            if subset.size == 0:
                print(f"No grid cells found "
                      f"time={time} "
                      f"lat={lat0:.2f} "
                      f"lon={lon0:.2f}", flush=True)
                continue

            lat2d, lon2d = xr.broadcast(subset.latitude, subset.longitude)

            dx = (lon2d - lon0) * np.cos(np.deg2rad(lat0)) * 111.0 # approximate conversion to km
            dy = (lat2d - lat0) * 111.0 # approximate conversion to km

            distance_km = np.sqrt(dx**2 + dy**2)

            masked = subset.where(distance_km <= radius_km)

            masked_values = masked.values

            if np.all(np.isnan(masked_values)):
                print(f"All winds NaN "
                      f"time={time} "
                      f"lat={lat0:.2f} "
                      f"lon={lon0:.2f}", flush=True)
                continue

            if np.all(np.isnan(masked_values)):
                continue

            iy, ix = np.unravel_index(np.nanargmax(masked_values), masked_values.shape)

            # Write the maximum wind and RMW back to the original track DataFrame
            track.loc[idx, wind_col] = float(masked_values[iy, ix])
            track.loc[idx, rmw_col] = float(distance_km.values[iy, ix])

            track.loc[idx, "rmw_lat"] = float(subset.latitude.values[iy])
            track.loc[idx, "rmw_lon"] = float(subset.longitude.values[ix])

        except Exception as error:
            print(f"FAILED timestep {idx} "
                  f"time={time} "
                  f"lat={lat0} "
                  f"lon={lon0}: {error}", flush=True)
            raise

    return track


def main(run_tracking=False, add_wind=False, plot_tracks=False, plot_track_context=False, 
         wind_radius_km=200, use_ibtracs_rmw=False):
    import os
    import glob
    import pickle
    import pandas as pd
    from multiprocessing import Pool
    import platform
    from pathlib import Path

    if platform.system() == "Windows":
        data_path_base = Path("c:/Code/test/ClimateDT/")
        processed_path = Path("P:/11210471-001-compass/01_Data/ECMWF_ClimateDT/analysis_output/")
        processed_path_tracks = os.path.join(processed_path, "SI_tracked_storms")
        processed_path_tracked_wind = os.path.join(processed_path, "SI_tracked_storms")
        figure_path = processed_path / "figures"
    else:
        data_path_base = Path("/projects/prjs2226/data/ClimateDT/sfc/raw/")
        processed_path = Path("/projects/prjs2226/data/ClimateDT/sfc/processed/")
        processed_path_tracks = Path(os.path.join(processed_path, "msl", "tracked"))
        processed_path_tracked_wind = Path(os.path.join(processed_path, "wind", "tracked"))
        figure_path = processed_path / "figures"

    os.makedirs(processed_path, exist_ok=True)
    os.makedirs(figure_path, exist_ok=True)
    os.makedirs(os.path.join(processed_path_tracks), exist_ok=True)
    os.makedirs(os.path.join(processed_path_tracked_wind), exist_ok=True)


    start_time = "2017-01-01"
    end_time = "2026-08-01"

    # output files
    track_file = os.path.join(processed_path_tracks,
                              f"storms_tracked_{start_time}_{end_time}.pkl")

    wind_track_file = os.path.join(processed_path_tracked_wind,
        f"storms_tracked_max_wind_{wind_radius_km}km_{start_time}_{end_time}.pkl")

    # Loading IBTrACS data in SI basin for the specified period
    print("Importing IBTrACS data...")
    storms = import_ibtracs_period(start=start_time, end=end_time)

    # --------------------------------------------------
    # Track storms in ClimateDT data
    # --------------------------------------------------
    if run_tracking:
        experiments = ["hist", "cont", "Tplus2.0K"]
        realizations = ["1", "2", "3", "4", "5"]

        jobs = []
        for scen in experiments:
            for realization in realizations:
                pkl_file = os.path.join(processed_path_tracks, f"storms_{scen}_r{realization}.pkl")

                if os.path.exists(pkl_file):
                    print(f"SKIP {scen} r{realization}", flush=True)
                    continue

                data_path_msl = data_path_base / "msl"
                jobs.append((scen, realization, storms, data_path_msl, processed_path_tracks))

        print(f"Submitting {len(jobs)} jobs", flush=True)

        saved_files = []
        with Pool(processes=5) as pool:
            for result in pool.imap_unordered(process_member_tracking, jobs):
                print(f"worker returned: {os.path.basename(result)}", flush=True,)

                saved_files.append(result)

        # Load completed results
        tc_climatedt = {}
        for pkl_file in sorted(glob.glob(os.path.join(processed_path_tracks, "storms_*_r*.pkl",))):   
            print(f"Loading {os.path.basename(pkl_file)}", flush=True,)

            with open(pkl_file, "rb") as f:
                tc_climatedt.update(pickle.load(f))

        print(f"Loaded {len(tc_climatedt)} tracks", flush=True,)

        # PLOT FIRST
        print("Plotting all tracks...", flush=True)
        if plot_tracks:
            plot_tracks_all(storms, tc_climatedt, figure_path=figure_path)

        # Plot context of IBTrACS tracks and highlight the selected storm
        if plot_track_context:
            print("Plotting IBTrACS context...")
            plot_ibtracs_context(storms, tc_climatedt)

        # SAVE COMBINED PICKLE
        with open(track_file, "wb") as f:
            pickle.dump(tc_climatedt, f, protocol=pickle.HIGHEST_PROTOCOL)

        print(f"Saved {len(tc_climatedt)} tracks", flush=True)

    else:
        print(f"Loading existing track file: {track_file}", flush=True)
        with open(track_file, "rb") as f:
            tc_climatedt = pickle.load(f)

        if plot_tracks:
            print("Plotting all tracks...", flush=True)
            plot_tracks_all(storms, tc_climatedt, figure_path=figure_path)

        # Plot context of IBTrACS tracks and highlight the selected storm
        if plot_track_context:
            print("Plotting IBTrACS context...")
            plot_ibtracs_context(storms, tc_climatedt)

    # --------------------------------------------------
    # Add wind speed within tracked TC centre radius
    # --------------------------------------------------
    if add_wind:
        experiments = ["hist", "cont", "Tplus2.0K"]
        realizations = ["1", "2", "3", "4", "5"]

        jobs = []
        for scen in experiments:
            for realization in realizations:
                input_file = os.path.join(processed_path_tracks,
                             f"storms_{scen}_r{realization}.pkl")

                if not os.path.exists(input_file):
                    print(f"Missing tracking file: {input_file}", flush=True)
                    continue

                jobs.append((scen, realization, input_file, data_path_base,
                             processed_path_tracked_wind, wind_radius_km, use_ibtracs_rmw))

        wind_files = []
        with Pool(processes=5) as pool:
            for result in pool.imap_unordered(process_member_wind, jobs):
                if result is not None:
                    wind_files.append(result)

        tc_climatedt_wind = {}
        for wind_file in sorted(wind_files):
            with open(wind_file, "rb") as f:
                tc_climatedt_wind.update(pickle.load(f))

        wind_track_file = os.path.join(processed_path_tracked_wind,
                                       f"storms_tracked_max_wind_{radius_km}km_"
                                       f"{start_time}_{end_time}.pkl")

        with open(wind_track_file, "wb") as f:
            pickle.dump(tc_climatedt_wind, f, protocol=pickle.HIGHEST_PROTOCOL)

        size_gb = os.path.getsize(wind_track_file) / 1024**3

        print(f"Saved {len(tc_climatedt_wind)} tracks to "
              f"{wind_track_file} ({size_gb:.3f} GB)", flush=True)

        tc_climatedt = tc_climatedt_wind

    else:
        print(f"Loading existing wind track file: {wind_track_file}", flush=True)
        with open(wind_track_file, "rb") as f:
            tc_climatedt = pickle.load(f)


    return storms, tc_climatedt


if __name__ == "__main__":
    storms, tc_climatedt = main(run_tracking=False, add_wind=True, 
                                plot_tracks=False, plot_track_context=False, 
                                wind_radius_km=200)


# %%






#%%
# import os
# import glob
# import xarray as xr
# from pathlib import Path
# import numpy as np

# data_path = Path("c:/Code/test/ClimateDT/")
# processed_path = Path("P:/11210471-001-compass/01_Data/ECMWF_ClimateDT/analysis_output/SI_tracked_storms")
# figure_path = processed_path / "figures"

# ds_10u = os.path.join(data_path, f"climateDT_10u_hist_r1_201901.nc")
# ds_10v = os.path.join(data_path, f"climateDT_10v_hist_r1_201901.nc")
# ds_tp = os.path.join(data_path, f"climateDT_avg_tprate_hist_r1_201901.nc")
# files = sorted(glob.glob(ds_10u)) + sorted(glob.glob(ds_10v)) + sorted(glob.glob(ds_tp))
# data = xr.open_mfdataset(files, combine="by_coords", data_vars="minimal", coords="minimal", 
#                              compat="override", parallel=False)
# data['wind_speed'] = np.sqrt(data['u10']**2 + data['v10']**2)
# data['tp'] = data['avg_tprate'] * 3600 # convert from kg/m2/s to mm/hr




# %%
# storm_name = "DESMOND"
# scenario = "hist"
# realization = "1"

# tc_desmond = storms["DESMOND"]
# rmw = storms["DESMOND"].rmw.max() * 1.852 # convert from nautical miles to km
# search_radius = rmw + 20 # adding 20 km for safety
# model_tracks = tc_climatedt[(scenario, realization, "DESMOND")]

# track_with_wind_200km = add_max_wind_within_radius(model_tracks, data['wind_speed'], radius_km=200)
# track_with_wind_rmw = add_max_wind_within_radius(model_tracks, data['wind_speed'], radius_km=search_radius)

# %%
# import matplotlib.pyplot as plt
# fig = plt.figure(figsize=(10, 5))

# fig.suptitle("Comparison of Maximum Wind Speed")

# plt.plot(track_with_wind_200km['time'], track_with_wind_200km['max_wind_model'].values)  # Model wind speed
# # plt.plot(track_with_wind_rmw['time'], track_with_wind_rmw['max_wind_model'].values)  # Model wind speed
# plt.plot(storms['DESMOND']['time'], storms['DESMOND']['wind'].values * 0.514444)  # Convert knots to m/s

# plt.xlabel("Time")
# plt.ylabel("Wind Speed (m/s)")
# plt.legend(["Model (200km)", "Model (RMW)", "Observed"])




# %%

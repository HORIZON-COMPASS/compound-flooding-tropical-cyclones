"""
Storm tracker strategies module.

Provides different initialization and tracking strategies for tropical cyclone tracking.
Each strategy uses different IBTrACS information to initialize the search and guide the tracking.
"""

import numpy as np
import xarray as xr
import pandas as pd
import bisect


def dist(a, b):
    """Calculates the magnitude of two variables in an xarray-compliant way."""
    func = lambda x, y: np.sqrt(x**2 + y**2)
    return xr.apply_ufunc(func, a, b)


def _find_closest_min_pressure(c_frame, guess_lat, guess_lon, s_size):
    """
    Find the minimum MSLP in search box and return closest to guess point.
    
    Parameters
    ----------
    c_frame : xr.DataArray
        Current time frame of MSLP data
    guess_lat : float
        Reference latitude
    guess_lon : float
        Reference longitude
    s_size : float
        Search box half-width (degrees)
    
    Returns
    -------
    pd.DataFrame
        DataFrame with lat, lon, msl of tracked point
    """
    search_area = c_frame.where(
        (c_frame.lat > (float(guess_lat) - s_size)) &
        (c_frame.lat < (float(guess_lat) + s_size)) &
        (c_frame.lon > (float(guess_lon) - s_size)) &
        (c_frame.lon < (float(guess_lon) + s_size)),
        drop=True
    ).compute() # selects the search area around the guessed or observed track point
    
    if search_area.size == 0:
        return None
    
    min_list  = search_area.where(search_area == search_area.min(), drop=True) # finds the cells equal to minumum pressure in the search area
    dist_list = dist(min_list.lat - float(guess_lat), min_list.lon - float(guess_lon)) # calculates the Euclidian distance from the guessed point to each of the minimum pressure points
    result    = min_list.where(dist_list == dist_list.min(), drop=True).to_dataframe().reset_index().iloc[[0]] # selects the closest minimum pressure point to the guessed point and returns it as a DataFrame, if multiple points tie, the first one is selected
    
    return result


def _track_forward(data, variables, min_p, track_time, best_track, s_size, enter):
    """
    Track forward from initialized points using pressure minima.
    
    Parameters
    ----------
    data : xr.Dataset
        Full dataset
    variables : xr.DataArray
        MSLP data
    min_p : pd.DataFrame
        Initial track points
    track_time : pd.DatetimeIndex
        Times to track
    best_track : pd.DataFrame
        IBTrACS reference
    s_size : float
        Search box size
    enter : pd.Timestamp
        Start time
    
    Returns
    -------
    pd.DataFrame
        Complete tracked storm path
    """
    idx_0 = bisect.bisect_left(best_track['time'], enter)
    n_init_points = len(min_p) # number of points already tracked
    
    for i in np.arange(n_init_points, len(track_time)):
        c_frame = variables.sel(time=enter + (data.time[i] - data.time[0])).compute() # time step of the model data to be tracked

        # guess_lat and guess_lon are the last tracked points, which are used as the reference for the next search
        guess_lat = min_p['lat'].iloc[i - 1] 
        guess_lon = min_p['lon'].iloc[i - 1]
        
        result = _find_closest_min_pressure(c_frame, guess_lat, guess_lon, s_size)
        if result is not None:
            min_p = pd.concat([min_p, result], ignore_index=True)
    
    return min_p


def strategy_2point_init(data, best_track, s_size):
    """
    STRATEGY 1 (CURRENT): Initialize with first 2 IBTrACS points, then extrapolate forward.
    
    Uses the first two observation points from IBTrACS to set up initial conditions,
    then tracks forward using pressure minimum search in expanding search area.
    
    Parameters
    ----------
    data : xr.Dataset
        Dataset with MSLP data (variables.msl)
    best_track : pd.DataFrame
        IBTrACS track with columns: time, lat, lon, min_pres
    s_size : float
        Search box half-width in degrees
    
    Returns
    -------
    pd.DataFrame
        Tracked storm positions (lat, lon, time, msl)
    """
    variables = data.msl
    
    best_track["time"] = pd.to_datetime(best_track["time"].dt.floor("h")) # rounding to hours to match time steps in the model data
    enter = max([variables.time[0], best_track.loc[best_track[["lat", "lon", "min_pres"]].notna().all(axis=1), "time"].iloc[0]])   # start of track with valid data
    exit  = min([variables.time[-1], best_track.loc[best_track[["lat", "lon", "min_pres"]].notna().all(axis=1), "time"].iloc[-1]]) # end of track with valid data
    
    timesteps = pd.to_datetime(data.time.values)
    idx_enter = bisect.bisect_left(timesteps, enter) # index of entrance time
    idx_exit = bisect.bisect(timesteps, exit)        # index of dissipation time
    track_time = timesteps[idx_enter:idx_exit]       # times to be tracked
    
    idx_0 = bisect.bisect_left(best_track['time'], enter) # index of first valid IBTrACS point
    
    # First point
    c_frame = variables.sel(time=enter).compute()
    result1 = _find_closest_min_pressure(c_frame, best_track.lat[idx_0], best_track.lon[idx_0], s_size)
    min_p = result1
    
    # Second point
    c_frame2 = variables.sel(time=enter + (data.time[1] - data.time[0])).compute()
    result2 = _find_closest_min_pressure(c_frame2, best_track.lat[idx_0], best_track.lon[idx_0], s_size)
    min_p = pd.concat([min_p, result2], ignore_index=True)
    
    # Track forward from point 2
    min_p = _track_forward(data, variables, min_p, track_time, best_track, s_size, enter)
    
    print("Strategy 2-point initialization: done")
    return min_p


def strategy_5point_init(data, best_track, s_size):
    """
    STRATEGY 2: Initialize with first 5 IBTrACS points, then track forward.
    
    Uses the first five observation points from IBTrACS to establish initial track direction
    and position, providing more robust initialization for subsequent tracking.
    
    Parameters
    ----------
    data : xr.Dataset
        Dataset with MSLP data (variables.msl)
    best_track : pd.DataFrame
        IBTrACS track with columns: time, lat, lon, min_pres
    s_size : float
        Search box half-width in degrees
    
    Returns
    -------
    pd.DataFrame
        Tracked storm positions (lat, lon, time, msl)
    """
    variables = data.msl
    
    best_track["time"] = pd.to_datetime(best_track["time"].dt.floor("h"))
    enter = max([variables.time[0], best_track.loc[best_track[["lat", "lon", "min_pres"]].notna().all(axis=1), "time"].iloc[0]])
    exit = min([variables.time[-1], best_track.loc[best_track[["lat", "lon", "min_pres"]].notna().all(axis=1), "time"].iloc[-1]])
    
    timesteps = pd.to_datetime(data.time.values)
    idx_enter = bisect.bisect_left(timesteps, enter)
    idx_exit = bisect.bisect(timesteps, exit)
    track_time = timesteps[idx_enter:idx_exit]
    
    idx_0 = bisect.bisect_left(best_track['time'], enter)
    
    # Initialize with first 5 points
    min_p = None
    for point_idx in range(min(5, len(track_time))):
        c_frame = variables.sel(time=enter + (data.time[point_idx] - data.time[0])).compute()
        result = _find_closest_min_pressure(c_frame, best_track.lat[idx_0], best_track.lon[idx_0], s_size)
        
        if result is not None:
            if min_p is None:
                min_p = result
            else:
                min_p = pd.concat([min_p, result], ignore_index=True)
    
    if min_p is None:
        raise ValueError("Could not initialize tracker with first 5 points")
    
    # Track forward from point 5 onwards
    min_p = _track_forward(data, variables, min_p, track_time, best_track, s_size, enter)
    
    print("Strategy 5-point initialization: done")
    return min_p


def strategy_all_ibtracs_init(data, best_track, s_size):
    """
    STRATEGY 3: Use all IBTrACS points to guide initialization, find global MSLP minimum.
    
    For each time step, finds the minimum MSLP in model data closest to the IBTrACS
    observation. This effectively allows the full IBTrACS track to guide the initialization
    and all subsequent tracking steps.
    
    Parameters
    ----------
    data : xr.Dataset
        Dataset with MSLP data (variables.msl)
    best_track : pd.DataFrame
        IBTrACS track with columns: time, lat, lon, min_pres
    s_size : float
        Search box half-width in degrees
    
    Returns
    -------
    pd.DataFrame
        Tracked storm positions (lat, lon, time, msl)
    """
    variables = data.msl
    
    best_track["time"] = pd.to_datetime(best_track["time"].dt.floor("h"))
    enter = max([variables.time[0], best_track.loc[best_track[["lat", "lon", "min_pres"]].notna().all(axis=1), "time"].iloc[0]])
    exit  = min([variables.time[-1], best_track.loc[best_track[["lat", "lon", "min_pres"]].notna().all(axis=1), "time"].iloc[-1]])
    
    timesteps = pd.to_datetime(data.time.values)
    idx_enter = bisect.bisect_left(timesteps, enter)
    idx_exit  = bisect.bisect(timesteps, exit)
    track_time  = timesteps[idx_enter:idx_exit]
    
    min_p = None
    
    # Track all points using complete IBTrACS guidance
    for i in range(len(track_time)):
        c_frame = variables.sel(time=enter + (data.time[i] - data.time[0])).compute()
        
        # Find closest IBTrACS point in time
        time_diff = np.abs(best_track['time'] - track_time[i])
        closest_ibtracs_idx = time_diff.argmin()
        
        ibtracs_lat = best_track.lat.iloc[closest_ibtracs_idx]
        ibtracs_lon = best_track.lon.iloc[closest_ibtracs_idx]
        
        result = _find_closest_min_pressure(c_frame, ibtracs_lat, ibtracs_lon, s_size)
        
        if result is not None:
            if min_p is None:
                min_p = result
            else:
                min_p = pd.concat([min_p, result], ignore_index=True)
    
    if min_p is None:
        raise ValueError("Could not track storm using all IBTrACS points")
    
    print("Strategy all-IBTrACS initialization: done")
    return min_p

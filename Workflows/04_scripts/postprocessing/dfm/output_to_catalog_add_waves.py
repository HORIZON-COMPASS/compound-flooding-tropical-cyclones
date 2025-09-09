#%% Adding DFM output to the SFINCS coastal coupling data catalog
# use compass-wflow environment
# Importing the necessary packages
import os
import hydromt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
import contextily as ctx
from pyproj import Transformer
import geopandas as gpd
import pandas as pd
import numpy as np
import yaml
from datetime import datetime
import xarray as xr
import cartopy.crs as ccrs
from pathlib import Path
from scipy.interpolate import PchipInterpolator
import matplotlib.dates as mdates

#%%
if "snakemake" in locals():
    his_path            = os.path.abspath(snakemake.input.his_file)
    path_data_cat       = os.path.abspath(snakemake.params.cf_data_cat)
    model_name          = snakemake.params.model_name
    root_dir            = os.path.abspath(snakemake.params.root_dir)
    snake_done          = os.path.abspath(snakemake.output.done_file)
    use_wave            = snakemake.params.use_waves
    path_data_cat_coast = os.path.abspath(snakemake.params.coast_data_cat)
    wave_output         = snakemake.params.wave_output
    start_date          = snakemake.params.start_time
    end_date            = snakemake.params.start_time
else:
    region              = "sofala"
    tc_name             = "Idai"
    dfm_res             = "450"
    bathy               = "gebco2024_MZB"
    tidemodel           = 'GTSMv41' # tidemodel: FES2014, FES2012, EOT20, GTSMv41, GTSMv41opendap
    wind_forcing        = "era5_hourly_spw_IBTrACS"
    CF_SLR_txt          = "0"
    CF_wind_txt         = "0"
    start_date          = "20190309 000000"
    end_date            = "20190325 060000"
    model_name          = f'event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}'
    path_data_cat       = "C:/Code/COMPASS/Workflows/03_data_catalogs/datacatalog_CF_forcing.yml"
    run_dir             = f'p:/11210471-001-compass/03_Runs/{region}/{tc_name}/dfm/{model_name}'
    root_dir            = 'p:\\'
    snake_done          = os.path.join(run_dir, "postprocessing_done.txt")
    his_path1           = os.path.join(run_dir, "output", f"{model_name}_0000_his.nc")
    his_path2           = os.path.join(run_dir, "output", f"settings_0000_his.nc")  # Alternative file name
    # Use the first path that exists
    his_path            = his_path1 if os.path.exists(his_path1) else his_path2
    wave_output         = 'SNAPWAVE_Idai_setup'
    path_data_cat_coast = "C:/Code/COMPASS/Workflows/03_data_catalogs/datacatalog_SFINCS_coastal_coupling.yml"
    use_wave            = True

#%%
script_dir = Path(__file__).resolve().parent
config_file = script_dir / "../../../01_config_snakemake/config_general_MZB.yml"

with open(config_file, "r") as f:
    config = yaml.safe_load(f)
    config = config['runname_ids']['Idai']

# %% Loading the SFINCS coastal coupling data catalog & DFM output path
datacatalog = hydromt.DataCatalog(data_libs=[path_data_cat])
dfm_run = f"dfm_output_{model_name}"

#%% Specifying the DFM output information for the data catalog entry
adapter = hydromt.data_adapter.GeoDatasetAdapter(
    path=os.path.abspath(his_path),
    driver="netcdf",
    driver_kwargs={
        "chunks": {
            "station": 10,
            "time": -1
        }
    },
    rename={
        "station_x_coordinate": "lon",
        "station_y_coordinate": "lat"
    },
    meta={
        "category": "ocean"
    },
    crs=4326
    )

#%% Add the DFM output to the catalog and save
if dfm_run not in datacatalog:
    # Add new source if it doesn't exist
    datacatalog.add_source(dfm_run, adapter)
    print(f"Dataset '{dfm_run}' has been added to the catalog.")
else:
    # Update the existing source
    datacatalog[dfm_run] = adapter
    print(f"Dataset '{dfm_run}' has been updated in the catalog.")

# Save the updated catalog
datacatalog.to_yml(path_data_cat, root=root_dir)
   
# %% ------------------------------------------ #
# Add simulated SNAPWAVE waves to the DFM output
# Using model output and code developed by Fernaldi Gradiyanto & Tim Leijnse
# --------------------------------------------- #

# Add wave component to DFM derived waterlevel
if use_wave:
    # Convert to datetime objects
    start_dt = datetime.strptime(start_date, "%Y%m%d %H%M%S")
    end_dt = datetime.strptime(end_date, "%Y%m%d %H%M%S")
    time_range = (start_dt, end_dt)

    # Load the right datacatalog
    datacatalog_coast = hydromt.DataCatalog(data_libs=[path_data_cat_coast])
    print("Loading wave data")
    ds_wave = datacatalog_coast.get_geodataset(wave_output)
    # Transform to WGS84 (lat/lon)
    ds_wave.rio.set_crs("EPSG:32736", inplace=True) # same crs as in datacatalog
    transformer = Transformer.from_crs(32736, 4326, always_xy=True)
    x_coords = ds_wave['x'].values 
    y_coords = ds_wave['y'].values
    ds_wave = ds_wave.rename({'transects': 'station'})
    lon, lat = transformer.transform(x_coords, y_coords)
    ds_wave = ds_wave.assign_coords(lon=("station", lon), lat=("station", lat))
    ds_wave = ds_wave.drop_vars(["x", "y"])  # Drop x and y coordinates if they exist
    
    # Plot to check data for station 40 as example
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    ds_wave.sel(station=40).plot(ax=axes[0])
    ds_wave.sel(time=slice(pd.to_datetime("2019-03-24"), pd.to_datetime("2019-03-25 06:00:00")), station=40).plot(ax=axes[1])
    axes[0].set_title("Full time series (station 40)")
    axes[1].set_title("Subset time series (station 40)")
    plt.tight_layout()
    plt.show()

    # Set wave values to 0 after cutoff for all stations due to incorrect data
    cutoff = pd.to_datetime("2019-03-24 21:00:00")
    # Replace values after cutoff with these averages
    ds_wave = ds_wave.where(ds_wave.time <= cutoff, np.nan)

    # Check correction
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    ds_wave.sel(station=40).plot(ax=axes[0])
    ds_wave.sel(time=slice(pd.to_datetime("2019-03-24"), pd.to_datetime("2019-03-25 06:00:00")), station=40).plot(ax=axes[1])
    axes[0].set_title("Full time series (station 40)")
    axes[1].set_title("Subset time series (station 40)")
    plt.tight_layout()
    plt.show()

    # Get the dfm run output
    print("Loading DFM data")
    ds_dfm = datacatalog.get_geodataset(dfm_run)
    ds_dfm = ds_dfm.sel(time=slice(*time_range))

    # --- Filter DFM DataArray (remove IHO stations) ---
    ds_filtered = ds_dfm['waterlevel'].sel(
        station=~ds_dfm['station_name'].astype(str).str.match(r'^[A-Za-z]')
    )

    # --- Convert DFM to GeoDataFrame ---
    print("Convert DFM to GeoDataFrame")
    df_dfm = pd.DataFrame({
        "lon": ds_filtered.lon.values,
        "lat": ds_filtered.lat.values,
        "max_component": ds_filtered.max(dim="time").values,
        "station_name": ds_filtered["station_name"].values,
        "station": ds_filtered["station"].values
    })

    gdf_dfm = gpd.GeoDataFrame(
        df_dfm,
        geometry=gpd.points_from_xy(df_dfm.lon, df_dfm.lat),
        crs="EPSG:4326"
    ).to_crs(epsg=32736)

    # --- Convert Wave DataArray to GeoDataFrame ---
    print("Convert Wave DataArray to GeoDataFrame")
    df_wave = pd.DataFrame({
        "lon": ds_wave.lon.values,
        "lat": ds_wave.lat.values,
        "station": ds_wave.station.values,
        "max_component": ds_wave.max(dim="time").values
    })

    gdf_wave = gpd.GeoDataFrame(
        df_wave,
        geometry=gpd.points_from_xy(df_wave.lon, df_wave.lat),
        crs="EPSG:4326"
    ).to_crs(epsg=32736)


    # ---- MATCHING DFM TIDE & SURE TO WAVE OUTPUT ----
    # Find the nearest wave point to each DFM point
    print("MATCHING DFM TIDE & SURE TO WAVE OUTPUT")
    gdf_dfm_matched = gdf_dfm.sjoin_nearest(
        gdf_wave,  # Only join necessary columns
        how='left',
        lsuffix='dfm',
        rsuffix='wave',
        distance_col='distance_to_wave')
    
    # Plot which wave transects are matched to which DFM 5 m depth-contour points
    # Wave points in light blue
    gdf_wave_latlon = gdf_wave.to_crs("EPSG:4326")
    ax = gdf_wave_latlon.plot(color='lightblue', markersize=20, label='Wave points', figsize=(8, 8))

    # Plot DFM points in red
    gdf_dfm_latlon = gdf_dfm.to_crs("EPSG:4326")
    gdf_dfm_latlon.plot(ax=ax, color='red', markersize=20, label='DFM stations')

    # Optionally, plot lines connecting DFM to nearest wave point
    for _, row in gdf_dfm_matched.iterrows():
        x_vals = [row['lon_dfm'], row['lon_wave']]
        y_vals = [row['lat_dfm'], row['lat_wave']]
        plt.plot(x_vals, y_vals, color='gray', linestyle='--', linewidth=0.7)

    # Annotate DFM stations
    for _, row in gdf_dfm_latlon.iterrows():
        ax.text(row['lon']+0.01, row['lat'], row['station'], fontsize=6, ha='left', va='bottom', color='black')

    for _, row in gdf_wave_latlon.iterrows():
        ax.text(row['lon'], row['lat'], row['station'], fontsize=7, ha='right', va='bottom', color='black')

    # Labels and legend
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.title("DFM stations matched to nearest wave points")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Combine wave and DFM data
    # Select wave data for each matched DFM station
    wave_for_dfm = ds_wave.sel(station=xr.DataArray(gdf_dfm_matched["station_wave"].values, dims="station"))
    ds_dfm_matched = ds_filtered.sel(station=gdf_dfm_matched["station_dfm"].values)
    wave_for_dfm = wave_for_dfm.assign_coords(station=ds_dfm_matched.station)

    # Now combine
    ds_combined = xr.Dataset({
        "tide_surge": ds_dfm_matched,
        "wave_setup": wave_for_dfm
    })

    # Compute total water level
    ds_combined["waterlevel"] = ds_combined["tide_surge"] + ds_combined["wave_setup"].fillna(0)
    ds_combined["waterlevel"].attrs = ds_dfm["waterlevel"].attrs.copy()
   
    print("Export the combined dataset")
    output_path = Path(datacatalog[dfm_run].path).parent
    # ds_combined.to_netcdf(os.path.join(output_path, "settings_0000_his_WAVES.nc"))

    # Add the DFM output with the added waves to the catalog and save
    print("Add the combined dataset to the datacatalog ")
    dfm_run_waves = f"dfm_output_{model_name}_waves"

    adapter = hydromt.data_adapter.GeoDatasetAdapter(
        path=os.path.join(output_path, "settings_0000_his_WAVES.nc"),
        driver="netcdf",
        driver_kwargs={
            "chunks": {
                "station": 10,
                "time": -1
            }
        },
        rename={
            "station_x_coordinate": "lon",
            "station_y_coordinate": "lat"
        },
        meta={
            "category": "ocean",
            "notes": "DFM output combined with wave data from SNAPWAVE provides by Fernaldi Gradyanto",
        },
        crs=4326
        )

    if dfm_run_waves not in datacatalog:
        # Add new source if it doesn't exist
        datacatalog.add_source(dfm_run_waves, adapter)
        print(f"Dataset '{dfm_run_waves}' has been added to the catalog.")
    else:
        # Update the existing source
        datacatalog[dfm_run_waves] = adapter
        print(f"Dataset '{dfm_run_waves}' has been updated in the catalog.")

    # Save the updated catalog
    datacatalog.to_yml(path_data_cat, root=root_dir)

else:
    pass


#%% Make a file for snakemake to track the cata,log update
with open(snake_done, 'w') as file:
    # Write content to the file
    file.write("# Empty file used to make the snakemake dfm workflow add_data_to_catalog rule work\n")

#%% ###########################################
#############    SOME PLOTTING    #############
###############################################
# plot variable for one of the matched stations
def plot_station_wave_components(station_idx=30):
    station_to_plot = ds_combined["station"].values[station_idx]

    tide_surge = ds_combined["tide_surge"].sel(station=station_to_plot).compute()
    wave_wl = ds_combined["wave_setup"].sel(station=station_to_plot).compute()
    waterlevel = ds_combined["waterlevel"].sel(station=station_to_plot).compute()

    plt.figure(figsize=(9,6))
    plt.plot(waterlevel["time"], waterlevel, label="Total water level")
    plt.plot(wave_wl["time"], wave_wl, label="Wave Induced WL")
    # plt.plot(tide_surge["time"], tide_surge, label="Tide + Surge Induced WL")

    wave_max = ds_combined["wave_setup"].max(dim='time').values[station_idx]

    plt.xlabel("Time")
    plt.ylabel("Water level (m)")
    plt.title(f"Water levels for station {station_to_plot} with max wave setup of {wave_max:.2f} m")
    plt.legend()
    plt.tight_layout()
    plt.xlim(tide_surge["time"].min().values, tide_surge["time"].max().values)
    plt.grid(True)
    plt.show()

# Plot the selected staions and their tide+surge, waves only and combined
def plot_wave_impact(station_idx = 30):
    # Plot the selected staions and their tide+surge, waves only and combined
    # ------------------------------------------------------------------
    # MAP: max(wave_setup)
    # ------------------------------------------------------------------
    # Create figure and two axes
    fig = plt.figure(figsize=(10, 6), dpi=300)
    # First subplot with projection
    ax_map = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
    # Second subplot without projection (regular 2D plot)
    ax_abs = fig.add_subplot(1, 2, 2)
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:4326", always_xy=True)

    # Plot which wave transects are matched to which DFM 5 m depth-contour points
    # Wave points in light blue
    gdf_wave_subset = gdf_wave_latlon[
        (gdf_wave_latlon['lat'] >= -20.3) & (gdf_wave_latlon['lat'] <= -19.58) &
        (gdf_wave_latlon['lon'] >= 34.7) & (gdf_wave_latlon['lon'] <= 35.22)]    
    gdf_wave_subset.plot(ax=ax_map, color='lightblue', markersize=40, edgecolor="k",
                         linewidth=0.2, label='Wave output points')

    vmin = 0.0
    vmax = (ds_combined["waterlevel"].max(dim='time').values).max()
    vcenter = (vmin + vmax) / 2.0
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

    sc1 = ax_map.scatter(
            ds_combined.lon,
            ds_combined.lat,
            c=ds_combined["waterlevel"].max(dim='time').values,
            cmap="viridis",
            norm=norm,
            s=40,
            edgecolor="k",
            linewidth=0.2,
            label = "DFM output points"
        )
    
    # Plot lines connecting DFM to nearest wave point
    label_added = False
    for _, row in gdf_dfm_matched.iterrows():
        x_vals = [row['lon_dfm'], row['lon_wave']]
        y_vals = [row['lat_dfm'], row['lat_wave']]
        
        if not label_added:
            ax_map.plot(x_vals, y_vals, color='gray', linestyle='--', linewidth=0.7, label='Match line')
            label_added = True
        else:
            ax_map.plot(x_vals, y_vals, color='gray', linestyle='--', linewidth=0.7)

    # Extract coordinates and value
    station = ds_combined["station"].values[station_idx]
    lon = ds_combined.lon.sel(station=station).values
    lat = ds_combined.lat.sel(station=station).values
    ax_map.annotate(
        f"Station {station}",           # text
        xy=(lon, lat),               # point coordinates
        xytext=(lon + 0.03, lat),  # offset in lon/lat
        arrowprops=dict(arrowstyle='->', color='white', lw=1),
        fontsize=8,
        color='white',
        fontweight='bold'
    )

    # Labels and legend
    ax_map.set_title("DFM stations matched to nearest \nwave points", fontsize=11)

    ax_map.set_extent([ds_combined['lon'].min().values - 0.1, 
                       ds_combined['lon'].max().values + 0.025,
                       ds_combined['lat'].min().values - 0.03, 
                    ds_combined['lat'].max().values + 0.03], crs=ccrs.PlateCarree())
  
    # Add basemap (LOWER zoom = faster)
    ctx.add_basemap(ax_map, source=ctx.providers.Esri.WorldImagery, zoom=10, 
                    crs="EPSG:4326", attribution=False, zorder=0)

    txt = ax_map.text(
        34.94, -20.3,  # x, y in figure coordinates (0=left/bottom, 1=right/top)
        "Tiles © Esri -- Source: Esri, i-cubed, USDA, USGS, " \
        "\nAEX, GeoEye, Getmapping, Aerogrid, IGN, IGP," \
        "\nUPR-EGP, and the GIS User Community",
        fontsize=5.5,
        color='white',
        alpha=0.7,
        ha='left',
        va='bottom',
        zorder=20,
    )

    lon_ticks = np.linspace(*ax_map.get_xlim(), 5)
    lat_ticks = np.linspace(*ax_map.get_ylim(), 5)
    ax_map.set_xticks(lon_ticks)
    ax_map.set_yticks(lat_ticks)
    ax_map.set_xticklabels([f"{transformer.transform(x, df_dfm.lat.min())[0]:.1f}°E" for x in lon_ticks], fontsize=11)
    ax_map.set_yticklabels([f"{transformer.transform(df_dfm.lon.min(), y)[1]:.1f}°S" for y in lat_ticks], fontsize=11)
    ax_map.legend(fontsize=9)
    sm = plt.cm.ScalarMappable(norm=norm, cmap="viridis")
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax_map, orientation="vertical", fraction=0.04, pad=0.02)
    cbar.set_label("DFM max combined water level (m)")

    # ------------------------------------------------------------------
    # BAR CHART: Combined dataset
    # ------------------------------------------------------------------
    df_combined = pd.DataFrame({
            'lon': ds_combined.lon.values,
            'lat': ds_combined.lat.values,
            'wl_max_component': ds_combined["waterlevel"].max(dim='time').values,
            'wave_max_component': ds_combined["wave_setup"].max(dim='time').values,
            'tide_surge_max_component': ds_combined["tide_surge"].max(dim='time').values
        })
       
    ax_abs.barh(df_combined['lat'], df_combined["wl_max_component"],
                height=0.005, color="green", label="Max combined water level")
    
    ax_abs.barh(df_combined['lat'], df_combined["tide_surge_max_component"],
                height=0.005, color="steelblue", label="Max water level from tide & surge")

    ax_abs.barh(df_combined['lat'], df_combined["wave_max_component"],
                height=0.005, color="orange", label="Max wave setup")

    ax_abs.set_xlabel("Max water level component (m)")
    ax_abs.set_title("Max water level component at DFM points", fontsize=11)
    ax_abs.grid(axis="x", alpha=0.3)
    ax_abs.legend(fontsize=9, loc="lower left")
    ax_abs.set_yticklabels([])

    # Add subplot labels
    ax_map.text(0.07, 1.04, "(a)", transform=ax_map.transAxes, fontsize=12, 
                fontweight='bold', va='top', ha='right')
    ax_abs.text(0.02, 1.04, "(b)", transform=ax_abs.transAxes, fontsize=12, 
                fontweight='bold', va='top', ha='right')
    
    plt.show()


def max_wave_ratio(ds_combined):
    # 1. Find the time index where waterlevel is maximum per station
    idx_max_wl = ds_combined["waterlevel"].argmax(dim="time").compute()  # numpy array now

    # 2. Select wave setup at those times per station
    wave_setup_at_max_wl = ds_combined["wave_setup"].isel(
        time=xr.DataArray(idx_max_wl, dims="station")
    )

    # 3. Get summary stats (compute for Dask arrays)
    max_wl = ds_combined["waterlevel"].max(dim="time").max().compute().item()
    max_wave_setup_at_max_wl = wave_setup_at_max_wl.max().compute().item()
    min_wave_setup_at_max_wl = wave_setup_at_max_wl.min().compute().item()

    # 4. Ratio at max WL
    wl_at_max = ds_combined["waterlevel"].isel(
        time=xr.DataArray(idx_max_wl, dims="station")
    )
    ratio_at_max_wl = (wave_setup_at_max_wl / wl_at_max).max().compute().item()

    print(f"Maximum total water level: {max_wl:.2f} m")
    print(f"Wave setup during max WL: {min_wave_setup_at_max_wl:.2f}–{max_wave_setup_at_max_wl:.2f} m")
    print(f"Max ratio (wave_setup / waterlevel at max WL): {ratio_at_max_wl*100:.0f}%")


# Plot the selected staions and their tide+surge, waves only and combined
def plot_wave_overview(station_idx = 30):
    # Plot the selected staions and their tide+surge, waves only and combined
    # ------------------------------------------------------------------
    # MAP: max(wave_setup)
    # ------------------------------------------------------------------
    # Define figure with unequal subplot widths
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=(12, 6), dpi=300)

# Map subplot: full height, left side
    ax_map = fig.add_axes([0.05, 0.05, 0.5, 0.8], projection=ccrs.PlateCarree())

    # Time series subplot: 2/3 height, vertically centered, right side
    ax_abs = fig.add_axes([0.65, 0.2, 0.3, 0.45])
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:4326", always_xy=True)

    # Plot which wave transects are matched to which DFM 5 m depth-contour points
    # Wave points in light blue
    gdf_wave_subset = gdf_wave_latlon[
        (gdf_wave_latlon['lat'] >= -20.3) & (gdf_wave_latlon['lat'] <= -19.58) &
        (gdf_wave_latlon['lon'] >= 34.7) & (gdf_wave_latlon['lon'] <= 35.22)]    
    gdf_wave_subset.plot(ax=ax_map, color='lightblue', markersize=40, edgecolor="k",
                         linewidth=0.2, label='Wave output points')

    vmin = 0.0
    vmax = (ds_combined["waterlevel"].max(dim='time').values).max()
    vcenter = (vmin + vmax) / 2.0
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

    sc1 = ax_map.scatter(
            ds_combined.lon,
            ds_combined.lat,
            c=ds_combined["waterlevel"].max(dim='time').values,
            cmap="viridis",
            norm=norm,
            s=40,
            edgecolor="k",
            linewidth=0.2,
            label = "DFM output points"
        )
    
    # Plot lines connecting DFM to nearest wave point
    label_added = False
    for _, row in gdf_dfm_matched.iterrows():
        x_vals = [row['lon_dfm'], row['lon_wave']]
        y_vals = [row['lat_dfm'], row['lat_wave']]
        
        if not label_added:
            ax_map.plot(x_vals, y_vals, color='gray', linestyle='--', linewidth=0.7, label='Match line')
            label_added = True
        else:
            ax_map.plot(x_vals, y_vals, color='gray', linestyle='--', linewidth=0.7)

    # Extract coordinates and value
    station = ds_combined["station"].values[station_idx]
    lon = ds_combined.lon.sel(station=station).values
    lat = ds_combined.lat.sel(station=station).values
    ax_map.annotate(
        f"Station {station}",           # text
        xy=(lon, lat),               # point coordinates
        xytext=(lon + 0.03, lat),  # offset in lon/lat
        arrowprops=dict(arrowstyle='->', color='white', lw=1),
        fontsize=8,
        color='white',
        fontweight='bold'
    )

    # Labels and legend
    ax_map.set_title("DFM stations matched to nearest \nwave points", fontsize=11)

    ax_map.set_extent([ds_combined['lon'].min().values - 0.1, 
                       ds_combined['lon'].max().values + 0.025,
                       ds_combined['lat'].min().values - 0.03, 
                    ds_combined['lat'].max().values + 0.03], crs=ccrs.PlateCarree())
  
    # Add basemap (LOWER zoom = faster)
    ctx.add_basemap(ax_map, source=ctx.providers.Esri.WorldImagery, zoom=10, 
                    crs="EPSG:4326", attribution=False, zorder=0)

    txt = ax_map.text(
        34.94, -20.3,  # x, y in figure coordinates (0=left/bottom, 1=right/top)
        "Tiles © Esri -- Source: Esri, i-cubed, USDA, USGS, " \
        "\nAEX, GeoEye, Getmapping, Aerogrid, IGN, IGP," \
        "\nUPR-EGP, and the GIS User Community",
        fontsize=5.5,
        color='white',
        alpha=0.7,
        ha='left',
        va='bottom',
        zorder=20,
    )

    lon_ticks = np.linspace(*ax_map.get_xlim(), 5)
    lat_ticks = np.linspace(*ax_map.get_ylim(), 5)
    ax_map.set_xticks(lon_ticks)
    ax_map.set_yticks(lat_ticks)
    ax_map.set_xticklabels([f"{transformer.transform(x, df_dfm.lat.min())[0]:.1f}°E" for x in lon_ticks], fontsize=10)
    ax_map.set_yticklabels([f"{transformer.transform(df_dfm.lon.min(), y)[1]:.1f}°S" for y in lat_ticks], fontsize=10)
    ax_map.legend(fontsize=9)
    sm = plt.cm.ScalarMappable(norm=norm, cmap="viridis")
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax_map, orientation="vertical", fraction=0.04, pad=0.02, shrink=0.7)
    cbar.set_label("DFM max combined water level (m)", fontsize=11)

    # ---- Timeseries plotting (ax_abs) ----
    station_to_plot = ds_combined["station"].values[station_idx]
    tide_surge = ds_combined["tide_surge"].sel(station=station_to_plot, time=slice("2019-03-12","2019-03-17")).compute()
    wave_wl = ds_combined["wave_setup"].sel(station=station_to_plot, time=slice("2019-03-12","2019-03-17")).fillna(0).compute()
    waterlevel = ds_combined["waterlevel"].sel(station=station_to_plot, time=slice("2019-03-12","2019-03-17")).compute()

    ax_abs.plot(waterlevel["time"], waterlevel, label="Total water level")
    ax_abs.plot(wave_wl["time"], wave_wl, label="Wave setup")
    ax_abs.plot(tide_surge["time"], tide_surge, label="Tide + surge")

    wave_max = ds_combined["wave_setup"].max(dim='time').values[station_idx]

    ax_abs.set_ylabel("Water level (m)", fontsize=11)
    ax_abs.set_title(f"Water levels at station {station_to_plot}\nMax wave setup = {wave_max:.2f} m", fontsize=11)
    ax_abs.legend(fontsize=9)
    ax_abs.set_xlim(tide_surge["time"].min().values, tide_surge["time"].max().values)
    ax_abs.set_xlabel("Day in March 2019", fontsize=10)
    ax_abs.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    ax_abs.tick_params(axis='both', labelsize=9)
    ax_abs.grid(True)   

    ax_map.text(0.0, 1.05, "(a)", transform=ax_map.transAxes,
             fontsize=12, fontweight='bold', va='top', ha='left')

    # Subplot (b) - second plot
    ax_abs.text(0.0, 1.08, "(b)", transform=ax_abs.transAxes,
                fontsize=12, fontweight='bold', va='top', ha='left')  


    fig.tight_layout()

    fig.savefig('../../../../Attribution_results/figures/fS7_waves.png', dpi=300, bbox_inches = 'tight')
    fig.savefig('../../../../Attribution_results/figures/fS7_waves.pdf', dpi=300, bbox_inches = 'tight')
    return fig, ax_map, ax_abs




if use_wave:
    # plot_wave_impact()
    # plot_station_wave_components()
    max_wave_ratio(ds_combined)
    # plot_wave_overview()


# %%

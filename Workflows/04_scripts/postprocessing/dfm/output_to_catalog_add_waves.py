#%% Adding DFM output to the SFINCS coastal coupling data catalog
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

#%%
if "snakemake" in locals():
    his_path = os.path.abspath(snakemake.input.his_file)
    path_data_cat = os.path.abspath(snakemake.params.cf_data_cat)
    model_name = snakemake.params.model_name
    root_dir = os.path.abspath(snakemake.params.root_dir)
    snake_done = os.path.abspath(snakemake.output.done_file)
    use_wave = snakemake.params.use_waves
    path_data_cat_coast = os.path.abspath(snakemake.params.coast_data_cat)
    wave_output = snakemake.params.wave_output
else:
    region = "sofala"
    tc_name = "Idai"
    dfm_res = "450"
    bathy = "gebco2024_MZB"
    tidemodel = 'GTSMv41' # tidemodel: FES2014, FES2012, EOT20, GTSMv41, GTSMv41opendap
    wind_forcing = "spw_IBTrACS"
    CF_SLR_txt = "0"
    CF_wind_txt = "0"
    start_date = "20190309 000000"
    end_date = "20190325 060000"
    model_name = f'event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}'
    path_data_cat = os.path.abspath("../../../03_data_catalogs/datacatalog_CF_forcing.yml")
    run_dir = f'p:/11210471-001-compass/03_Runs/{region}/{tc_name}/dfm/{model_name}'
    root_dir = 'p:\\'
    snake_done = os.path.join(run_dir, "postprocessing_done.txt")
    his_path1 = os.path.join(run_dir, "output", f"{model_name}_0000_his.nc")
    his_path2 = os.path.join(run_dir, "output", f"settings_0000_his.nc")  # Alternative file name
    # Use the first path that exists
    his_path = his_path1 if os.path.exists(his_path1) else his_path2
    wave_output = 'SNAPWAVE_Idai_setup'
    path_data_cat_coast = os.path.abspath("../../../03_data_catalogs/datacatalog_SFINCS_coastal_coupling.yml")
    use_wave = True

#%%
config_file = "../../../01_config_snakemake/config_general_MZB.yml"
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
        "station_y_coordinate": "lat",
        "station": "index"
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
# Using model output and code developed by Fernaldi Gradiyanto
# --------------------------------------------- #
if use_wave:
    # Convert to datetime objects
    start_dt = datetime.strptime(start_date, "%Y%m%d %H%M%S")
    end_dt = datetime.strptime(end_date, "%Y%m%d %H%M%S")
    time_range = (start_dt, end_dt)

    # Load the right datacatalog
    datacatalog_coast = hydromt.DataCatalog(data_libs=[path_data_cat_coast])
    print("Loading wave data")
    ds_wave = datacatalog_coast.get_geodataset(wave_output)

    # Get the dfm run output
    print("Loading DFM data")
    ds_dfm = datacatalog.get_geodataframe(dfm_run, variables='waterlevel')
    ds_dfm = ds_dfm.sel(time=slice(*time_range))

    # --- Filter DFM DataArray (remove IHO stations) ---
    ds_filtered = ds_dfm.sel(
        station=~ds_dfm['station_name'].astype(str).str.match(r'^[A-Za-z]')
    )

    # --- Convert DFM to GeoDataFrame ---
    print("Convert DFM to GeoDataFrame")
    df_dfm = pd.DataFrame({
        "lon": ds_filtered.lon.values,
        "lat": ds_filtered.lat.values,
        "max_component": ds_filtered.max(dim="time").values,
        "station_name": ds_filtered["station_name"].values
    })
    gdf_dfm = gpd.GeoDataFrame(
        df_dfm,
        geometry=gpd.points_from_xy(df_dfm.lon, df_dfm.lat),
        crs="EPSG:4326"
    ).to_crs(epsg=3857)

    # --- Convert Wave DataArray to GeoDataFrame ---
    print("Convert Wave DataArray to GeoDataFrame")
    df_wave = pd.DataFrame({
        "lon": ds_wave.lon.values,
        "lat": ds_wave.lat.values,
        "max_component": ds_wave.max(dim="time").values
    })
    gdf_wave = gpd.GeoDataFrame(
        df_wave,
        geometry=gpd.points_from_xy(df_wave.lon, df_wave.lat),
        crs="EPSG:4326"
    ).to_crs(epsg=3857)


    # ---- MATCHING DFM TIDE & SURE TO WAVE OUTPUT ----
    # Find the nearest wave point to each DFM point
    print("MATCHING DFM TIDE & SURE TO WAVE OUTPUT")
    gdf_dfm_matched = gdf_dfm.sjoin_nearest(
        gdf_wave[['geometry', 'max_component']],  # Only join necessary columns
        how='left',
        distance_col='distance_to_wave')

    # Rename joined wave column for clarity
    gdf_dfm_matched = gdf_dfm_matched.rename(columns={"max_component_right": "wave_max_component", 
                                                      "max_component_left": "wl_max_component",
                                                      "index_right": "wave_index",
                                                      })

    # Filter to keep only those with distance <= 50 meters
    gdf_dfm_matched = gdf_dfm_matched[gdf_dfm_matched["distance_to_wave"] <= 50]
    gdf_dfm_matched.index.name = "station"

    # Subset the xarray DataArrays using the index pairs
    ds_dfm_matched = ds_filtered.isel(station=xr.DataArray(gdf_dfm_matched.index.values, dims="match"))
    ds_wave_matched = ds_wave.isel(station=xr.DataArray(gdf_dfm_matched["wave_index"].values, dims="match"))

    # Add wave_induced_wl from ds_wave_matched to ds_dfm_matched
    print("Combining the matched wave and DFM data in one dataset")
    ds_combined = xr.merge([ds_dfm_matched, ds_wave_matched], compat="override")

    # rename variables and sum them to get the waterlevel incl. waves
    ds_combined = ds_combined.rename({"waterlevel": "tide_surge", "wave_induced_wl": "wave_setup"})
    ds_combined["waterlevel"] = ds_combined["tide_surge"] + ds_combined["wave_setup"].where(~ds_combined["wave_setup"].isnull(), 0)

    print("Export the combined dataset")
    output_path = Path(datacatalog[dfm_run].path).parent
    ds_combined.to_netcdf(os.path.join(output_path, "settings_0000_his_WAVES.nc"))

    # Add the DFM output with the added waves to the catalog and save
    print("Add the combined dataset to the datacatalog ")
    dfm_run_waves = f"dfm_output_{model_name}_waves"

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
            "station_y_coordinate": "lat",
            "station": "index"
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

#%% ###########################################
#############    SOME PLOTTING    #############
###############################################
# plot variable for one of the matched stations
def plot_station_wave_components(station_idx=79):
    station_to_plot = ds_combined["match"].values[station_idx]

    tide_surge = ds_combined["tide_surge"].sel(match=station_to_plot).compute()
    wave_wl = ds_combined["wave_setup"].sel(match=station_to_plot).compute()
    waterlevel = ds_combined["waterlevel"].sel(match=station_to_plot).compute()

    plt.figure(figsize=(12,6))
    plt.plot(waterlevel["time"], waterlevel, label="Waterlevel")
    plt.plot(wave_wl["time"], wave_wl, label="Wave Induced WL")
    plt.plot(tide_surge["time"], tide_surge, label="Tide + Surge Induced WL")

    wave_max = ds_combined["wave_setup"].max(dim='time').values[station_idx]

    plt.xlabel("Time")
    plt.ylabel("Water Level (m)")
    plt.title(f"Water levels for station {station_to_plot} with max wave setup of {wave_max:.2f} m")
    plt.legend()
    plt.tight_layout()
    plt.show()

# Plot the selected staions and their tide+surge, waves only and combined
def plot_wave_impact():
    # Plot the selected staions and their tide+surge, waves only and combined
    # ------------------------------------------------------------------
    # MAP: max(wave_setup)
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(14, 6), dpi=150)
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.2, 1, 1], wspace=0.3)
    ax_map = fig.add_subplot(gs[0, 0], projection = ccrs.PlateCarree())

    transformer = Transformer.from_crs("EPSG:4326", "EPSG:4326", always_xy=True)

    vmin = 0.0
    vmax = df_dfm["max_component"].max()
    vcenter = (vmin + vmax) / 2.0
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

    sc1 = ax_map.scatter(
        df_dfm.lon,
        df_dfm.lat,
        c=df_dfm["max_component"],
        cmap="viridis",
        norm=norm,
        s=40,
        edgecolor="k",
        linewidth=0.2,
        label="DFM"
    )

    df_wave_cut = df_wave[
        (df_wave.lon >= df_dfm.lon.min()) & (df_wave.lon <= df_dfm.lon.max()) &
        (df_wave.lat >= df_dfm.lat.min()) & (df_wave.lat <= df_dfm.lat.min())
    ]

    ctx.add_basemap(ax_map, source=ctx.providers.OpenStreetMap.Mapnik, crs="EPSG:4326", attribution=False)

    lon_ticks = np.linspace(*ax_map.get_xlim(), 5)
    lat_ticks = np.linspace(*ax_map.get_ylim(), 5)
    ax_map.set_xticks(lon_ticks)
    ax_map.set_yticks(lat_ticks)
    ax_map.set_xticklabels([f"{transformer.transform(x, df_dfm.lat.min())[0]:.1f}°E" for x in lon_ticks], fontsize=11)
    ax_map.set_yticklabels([f"{transformer.transform(df_dfm.lon.min(), y)[1]:.1f}°S" for y in lat_ticks], fontsize=11)
    ax_map.set_xlabel("Longitude", fontsize=11)
    ax_map.set_ylabel("Latitude", fontsize=11)
    ax_map.set_title("DFM output points", fontsize=12)

    sm = plt.cm.ScalarMappable(norm=norm, cmap="viridis")
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax_map, orientation="vertical", fraction=0.04, pad=0.02)
    cbar.set_label("Max Water Level Component (m)")

    ax_map.legend(loc="upper left", fontsize=11)

    # ------------------------------------------------------------------
    # BAR CHART 1: Original datasets
    # ------------------------------------------------------------------
    ax_abs = fig.add_subplot(gs[0, 1])
    lat = ds_filtered.lat.values
    y_min, y_max = lat.min() - 0.01, lat.max() + 0.01
    bar_h = (y_max - y_min) / len(lat) * 0.6

    ax_abs.barh(gdf_dfm_matched['lat'], gdf_dfm_matched["wl_max_component"],
                height=0.005, color="steelblue", label="Tide + Surge")

    ax_abs.barh(gdf_dfm_matched['lat'], gdf_dfm_matched["wave_max_component"],
                left=gdf_dfm_matched["wl_max_component"],
                height=0.005, color="orange", label="Wave Setup")

    ax_abs.set_ylim(y_min, y_max)
    ax_abs.set_xlabel("Max Water Level Component (m)")
    ax_abs.set_title("Max of Tide & Surge + Wave")
    ax_abs.grid(axis="x", alpha=0.3)
    ax_abs.legend(fontsize=8, loc="lower right")
    ax_abs.set_yticklabels([])

    # ------------------------------------------------------------------
    # BAR CHART 2: Combined dataset
    # ------------------------------------------------------------------
    ax_abs = fig.add_subplot(gs[0, 2])
    df_combined = pd.DataFrame({
        'lon': ds_combined.lon.values,
        'lat': ds_combined.lat.values,
        'wl_max_component': ds_combined["waterlevel"].max(dim='time').values,
        'wave_max_component': ds_combined["wave_setup"].max(dim='time').values,
        'tide_surge_max_component': ds_combined["tide_surge"].max(dim='time').values
    })

    ax_abs.barh(df_combined['lat'], df_combined["wave_max_component"],
                left=df_combined["tide_surge_max_component"],
                height=0.005, color="orange", label="Wave Setup (on top of T+S)")

    ax_abs.barh(df_combined['lat'], df_combined["tide_surge_max_component"],
                height=0.005, color="steelblue", label="Tide + Surge")

    ax_abs.barh(df_combined['lat'], df_combined["wl_max_component"],
                height=0.005, color="green", label="Total waterlevel")

    ax_abs.set_ylim(y_min, y_max)
    ax_abs.set_xlabel("Max Water Level Component (m)")
    ax_abs.set_title("Max of Tide & Surge + Wave, and Combined")
    ax_abs.grid(axis="x", alpha=0.3)
    ax_abs.legend(fontsize=8, loc="lower right")
    ax_abs.set_yticklabels([])

    plt.show()


if use_wave:
    plot_wave_impact()
    plot_station_wave_components()

# %%

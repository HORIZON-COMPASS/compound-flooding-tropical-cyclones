#%% 
# import necessary packages
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import xarray as xr
import contextily as ctx
import cartopy.crs as ccrs
from scipy.spatial import cKDTree
from shapely.geometry import Polygon
import dfm_tools as dfmt
from adjustText import adjust_text

#%% 
# Load IHO station locations
IHO_stations = pd.read_csv("p:/11210471-001-compass/01_Data/Coastal_boundary/points/MZB_Sofala_IHO_obs.xyn", delim_whitespace=True, header=None, names=["X", "Y", "Station"])

# Load ERA5 data (not time corrected!!)
filename = r"p:\11210471-001-compass\01_Data\ERA5\Idai\wave\data_stream-wave_stepType-instant.nc"
era5 = xr.open_dataset(filename)
era5=era5.sel(valid_time=slice('2019-03-09','2019-03-25')) # this is specific to the current model setup
era5.load()

#%%
# Convert station data to a GeoDataFrame
gdf = gpd.GeoDataFrame(IHO_stations, 
                        geometry=gpd.points_from_xy(IHO_stations["X"], IHO_stations["Y"]), 
                        crs="EPSG:4326") 

#%%
# Extract ERA5 lat/lon grid
era5_lats = era5.latitude.values
era5_lons = era5.longitude.values
era5_grid_points = np.array([(lat, lon) for lat in era5_lats for lon in era5_lons])

# Create KDTree for nearest neighbor search
tree = cKDTree(era5_grid_points)

# Find nearest ERA5 grid point for each station
station_coords = IHO_stations[['Y', 'X']].values
distances, indices = tree.query(station_coords)

# Store nearest lat/lon
IHO_stations['nearest_lat'], IHO_stations['nearest_lon'] = zip(*[era5_grid_points[idx] for idx in indices])

# Check for stations without a valid nearest point
missing_stations = IHO_stations[distances > 1.0]  # Adjust threshold if needed
if not missing_stations.empty:
    print("Warning: No valid grid point found for these stations:\n", missing_stations)

# Plot all stations' SWH time series in a single plot
plt.figure(figsize=(12, 6))

for i, (index, station) in enumerate(IHO_stations.iterrows()):
    try:
        # Extract nearest grid SWH time series
        nearest_swh = era5.sel(latitude=station['nearest_lat'], longitude=station['nearest_lon'], method="nearest")['swh']

        # Plot with label
        plt.plot(nearest_swh.valid_time, nearest_swh, label=f"{station['Station']} ({station['nearest_lat']:.2f}, {station['nearest_lon']:.2f})")

    except KeyError:
        print(f"Skipping {station['Station']} – No valid nearest ERA5 grid point.")

plt.xlabel("Time")
plt.ylabel("SWH (m)")
plt.title("SWH Time Series for Stations (no data for BEIRA IHO & airport)")
plt.legend()
plt.grid(True)
plt.show()

#%%
# Create a figure with a standard Matplotlib axis (NO Cartopy)
fig, ax = plt.subplots(1, 1, figsize=(10, 7), subplot_kw={'projection': ccrs.PlateCarree()})

# Plot stations
ax.scatter(gdf.geometry.x, gdf.geometry.y, color="red", label="Stations")

# Create a list to store annotation text objects
texts = []

# Annotate stations
for i, txt in enumerate(IHO_stations["Station"]):
    text = ax.annotate(txt, (gdf.geometry.x[i], gdf.geometry.y[i]), fontsize=8, ha="left", color="lightgrey")
    
    texts.append(text)
    
    # Get the grid box for the station
    lat = IHO_stations["nearest_lat"].iloc[i]
    lon = IHO_stations["nearest_lon"].iloc[i]
    
    # Define the grid box boundaries (assuming 0.5 degree resolution, adjust accordingly)
    lat_min = lat - 0.25  # 0.5 degrees box around the station
    lat_max = lat + 0.25
    lon_min = lon - 0.25
    lon_max = lon + 0.25

    # Create a polygon for the ERA5 grid box
    grid_box = Polygon([(lon_min, lat_min), (lon_min, lat_max), (lon_max, lat_max), (lon_max, lat_min)])

    # Plot the grid box
    ax.plot(*grid_box.exterior.xy, color="blue", linestyle="--", alpha=0.7, label="ERA5 Grid Box" if i == 0 else "")

# Adjust the position of the annotations to avoid overlap
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="->", color='grey'))

# Add lat/lon ticks based on your data
lat_ticks = np.arange(np.floor(gdf.geometry.y.min()), np.ceil(gdf.geometry.y.max()) + 1, 1)
lon_ticks = np.arange(np.floor(gdf.geometry.x.min()), np.ceil(gdf.geometry.x.max()) + 1, 1)

ax.set_xticks(lon_ticks)
ax.set_yticks(lat_ticks)

# Label the axes
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")

# Add basemap
ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=9, crs=ccrs.PlateCarree(), attribution=False, zorder=0)

plt.legend(fontsize=8)
plt.show()

# %%

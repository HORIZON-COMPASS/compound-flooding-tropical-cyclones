"""
SFINCS Flood Animation Script
==============================
Creates an animation (MP4) of flood depth evolution over time from SFINCS model output.

Requires:
- ffmpeg installed: conda install ffmpeg -c conda-forge
- sfincs_map.nc file in the event folder

Author: @dumontgoulart
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import contextily as ctx
import geopandas as gpd
import pandas as pd
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

# ===== CONFIGURATION =====
# Set your event name here
EVENT_NAME = "Kenneth"  # Change this to: "Kenneth", "Freddy", etc.

# Scenario to animate (factual or counterfactual)
SCENARIO = "event_tp_era5_hourly_zarr_CF0_GTSMv41opendap_CF0_no_wind_CF0"  # Factual
# SCENARIO = "event_tp_era5_hourly_zarr_CF-8_GTSMv41opendap_CF0_no_wind_CF0"  # Counterfactual

# Animation settings
HMIN = 0.05  # Minimum water depth threshold (m)
VMIN = 0  # Colorbar minimum (m)
VMAX = 3  # Colorbar maximum (m)
STEP = 1  # One frame every <step> timesteps (increase to speed up)
INTERVAL = 250  # Milliseconds between frames
FPS = 4  # Frames per second for saved video
DPI = 200  # Resolution of saved video
ZOOM = 11  # Basemap zoom level (higher = more detail, lower for zoomed out view)
EXTENT_BUFFER = 0.30  # Buffer around data extent (0.30 = 30% zoom out)

# Base paths - update these as needed
BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs/test")
OUTPUT_DIR = Path("/p/11210471-001-compass/04_Results/CF_figs")

# TC Track paths
TC_TRACKS_BASE = Path("/p/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS")
EVENT_TO_TC_TRACKS = {
    "Freddy": [
        TC_TRACKS_BASE / "IBTrACS_FREDDY_part1.shp",
        TC_TRACKS_BASE / "IBTrACS_FREDDY_part2.shp",
    ],
    "Kenneth": [TC_TRACKS_BASE / "IBTrACS_KENNETH.shp"],
    "Idai": [TC_TRACKS_BASE / "IBTrACS_IDAI.shp"],
}

# ===== DYNAMIC FILE PATHS =====
# Construct file path to sfincs_map.nc
sfincs_map_file = BASE_RUN_PATH / EVENT_NAME / "sfincs" / SCENARIO / "sfincs_map.nc"

# Create output directory if it doesn't exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Processing SFINCS animation for event: {EVENT_NAME}")
print(f"Scenario: {SCENARIO}")
print(f"Input file: {sfincs_map_file}")

# ===== LOAD SFINCS OUTPUT =====
print("Loading SFINCS map output...")

if not sfincs_map_file.exists():
    raise FileNotFoundError(f"SFINCS map file not found: {sfincs_map_file}")

ds = xr.open_dataset(sfincs_map_file)
print(f"Dataset variables: {list(ds.data_vars)}")
print(f"Dataset dimensions: {dict(ds.dims)}")

# Get water depth variable (usually 'h' or calculate from zs - zb)
if "h" in ds.data_vars:
    da_h = ds["h"].copy()
    print("Using 'h' variable for water depth")
elif "zs" in ds.data_vars and "zb" in ds.data_vars:
    # Calculate water depth as zs (water surface level) - zb (bed level)
    da_h = (ds["zs"] - ds["zb"]).copy()
    print("Calculating water depth as zs - zb (water surface - bed level)")
elif "zs" in ds.data_vars:
    da_h = ds["zs"].copy()
    print("Warning: Using 'zs' directly (no 'zb' found for water depth calculation)")
else:
    # Try to find a suitable variable
    depth_vars = [
        v
        for v in ds.data_vars
        if "h" in v.lower() or "depth" in v.lower() or "zs" in v.lower()
    ]
    if depth_vars:
        da_h = ds[depth_vars[0]].copy()
        print(f"Using '{depth_vars[0]}' variable for water depth")
    else:
        raise ValueError(
            f"Could not find water depth variable. Available: {list(ds.data_vars)}"
        )

# Check time dimension
if "time" not in da_h.dims:
    raise ValueError("No time dimension found in data. Cannot create animation.")

print(f"Time steps available: {da_h.time.size}")
print(f"Time range: {da_h.time.values[0]} to {da_h.time.values[-1]}")

# ===== LOAD TC TRACK DATA =====
tc_track_gdf = None
tc_track_times = None

if EVENT_NAME in EVENT_TO_TC_TRACKS:
    tc_track_files = EVENT_TO_TC_TRACKS[EVENT_NAME]
    print(f"\nLoading TC track data for {EVENT_NAME}...")

    # Load all track parts
    track_parts = []
    for track_file in tc_track_files:
        if track_file.exists():
            gdf_track = gpd.read_file(track_file)
            # Convert to EPSG:4326 if needed
            if gdf_track.crs != "EPSG:4326":
                gdf_track = gdf_track.to_crs("EPSG:4326")
            track_parts.append(gdf_track)
            print(f"  Loaded: {track_file.name}")
        else:
            print(f"  Warning: Track file not found: {track_file}")

    if track_parts:
        # Combine all track parts
        tc_track_gdf = pd.concat(track_parts, ignore_index=True)
        print(f"  Total track points: {len(tc_track_gdf)}")

        # Extract timestamps from TC track
        # Look for time field in the shapefile
        time_fields = [
            col
            for col in tc_track_gdf.columns
            if "time" in col.lower() or "date" in col.lower() or "ISO_TIME" in col
        ]
        if time_fields:
            time_field = time_fields[0]
            print(f"  Using time field: {time_field}")
            tc_track_times = pd.to_datetime(tc_track_gdf[time_field])
            print(
                f"  TC track time range: {tc_track_times.min()} to {tc_track_times.max()}"
            )
        else:
            print(
                f"  Warning: No time field found in TC track. Available columns: {tc_track_gdf.columns.tolist()}"
            )
    else:
        print("  No valid track files found")
else:
    print(f"\nNo TC track configured for event: {EVENT_NAME}")

# ===== MASK WATER DEPTH =====
print(f"Masking water depth < {HMIN}m...")
da_h = da_h.where(da_h > HMIN)

# Drop spatial_ref if present (causes issues with plotting)
if "spatial_ref" in da_h.coords:
    da_h = da_h.drop("spatial_ref")

da_h.attrs.update(long_name="flood depth", unit="m")

# ===== DETERMINE COORDINATE NAMES =====
# SFINCS output may use different coordinate names
x_coord = None
y_coord = None

for possible_x in ["xc", "x", "lon", "longitude"]:
    if possible_x in da_h.coords or possible_x in da_h.dims:
        x_coord = possible_x
        break

for possible_y in ["yc", "y", "lat", "latitude"]:
    if possible_y in da_h.coords or possible_y in da_h.dims:
        y_coord = possible_y
        break

if x_coord is None or y_coord is None:
    print(f"Available coords: {list(da_h.coords)}")
    print(f"Available dims: {list(da_h.dims)}")
    raise ValueError("Could not determine x/y coordinate names")

print(f"Using coordinates: x={x_coord}, y={y_coord}")

# ===== GET CRS FOR BASEMAP =====
try:
    if "crs" in ds.attrs:
        data_crs = ds.attrs["crs"]
    elif "spatial_ref" in ds:
        data_crs = ds["spatial_ref"].attrs.get("crs_wkt", None)
    else:
        # Assume UTM based on coordinate values
        center_x = float(da_h[x_coord].mean())
        center_y = float(da_h[y_coord].mean())
        if center_x > 100000:  # Likely UTM
            # Try to infer EPSG from coordinate range
            data_crs = "EPSG:32737"  # Default to UTM 37S for Mozambique region
            print(f"Assuming CRS: {data_crs}")
        else:
            data_crs = "EPSG:4326"  # WGS84
    print(f"Data CRS: {data_crs}")
except Exception as e:
    print(f"Could not determine CRS: {e}")
    data_crs = "EPSG:32737"  # Default

# ===== CREATE ANIMATION =====
print("Creating animation...")

# Get extent for basemap (with buffer to show more context)
x_min_data, x_max_data = float(da_h[x_coord].min()), float(da_h[x_coord].max())
y_min_data, y_max_data = float(da_h[y_coord].min()), float(da_h[y_coord].max())

# Apply buffer to extent for zoom out (to show TC track movement better)
x_buffer = (x_max_data - x_min_data) * EXTENT_BUFFER
y_buffer = (y_max_data - y_min_data) * EXTENT_BUFFER
x_min = x_min_data - x_buffer
x_max = x_max_data + x_buffer
y_min = y_min_data - y_buffer
y_max = y_max_data + y_buffer

print(
    f"Data extent: X [{x_min_data:.2f}, {x_max_data:.2f}], Y [{y_min_data:.2f}, {y_max_data:.2f}]"
)
print(
    f"Buffered extent ({EXTENT_BUFFER*100:.0f}%): X [{x_min:.2f}, {x_max:.2f}], Y [{y_min:.2f}, {y_max:.2f}]"
)

# Calculate figure size based on data aspect ratio
data_width = x_max - x_min
data_height = y_max - y_min
aspect_ratio = data_width / data_height

# Set figure size to match data aspect ratio (with some padding for colorbar)
fig_height = 8
fig_width = fig_height * aspect_ratio + 2  # Add space for colorbar
fig_width = min(max(fig_width, 8), 16)  # Clamp between 8 and 16 inches

# Create figure with constrained layout to minimize white space
fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height), constrained_layout=True)

# Set axis limits FIRST (before adding basemap so tiles cover full extent)
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_aspect("equal")

# Add basemap (must be done after setting extent so tiles cover buffered area)
try:
    ctx.add_basemap(
        ax=ax,
        source=ctx.providers.OpenStreetMap.Mapnik,
        crs=data_crs,
        attribution=False,
        zorder=1,
        zoom=ZOOM,
    )
    print("Basemap added successfully")
except Exception as e:
    print(f"Could not add basemap: {e}")

# Plot initial frame
da_h0 = da_h.isel(time=0)
t0 = str(da_h0.time.values)[:19]  # Format timestamp

# Create the plot
cax_h = da_h0.plot(
    x=x_coord,
    y=y_coord,
    ax=ax,
    vmin=VMIN,
    vmax=VMAX,
    cmap=plt.cm.viridis,
    alpha=0.7,
    zorder=2,
    add_colorbar=True,
    cbar_kwargs={"shrink": 0.8, "pad": 0.02, "label": "Flood Depth [m]"},
)

ax.set_title(f"SFINCS water depth {t0}")

# Format axis labels with scientific notation
ax.ticklabel_format(style="scientific", axis="both", scilimits=(0, 0))
ax.xaxis.get_offset_text().set_fontsize(10)
ax.yaxis.get_offset_text().set_fontsize(10)

# Ensure axis limits are maintained after plotting
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

# ===== PLOT TC TRACK =====
tc_track_line = None
tc_track_dot = None

if tc_track_gdf is not None:
    print("Adding TC track to animation...")

    # Extract coordinates from track geometry
    # Convert to data CRS if needed (for plotting on the map)
    tc_track_plot = tc_track_gdf.to_crs(data_crs)

    # Get x, y coordinates from geometry
    track_x = tc_track_plot.geometry.x.values
    track_y = tc_track_plot.geometry.y.values

    # Plot full track as dashed blue line
    tc_track_line = ax.plot(
        track_x,
        track_y,
        color="blue",
        linestyle="--",
        linewidth=2,
        zorder=5,
        label="TC Track",
        alpha=0.8,
    )[0]

    # Create scatter plot for moving dot (initially at first position)
    tc_track_dot = ax.scatter(
        [track_x[0]],
        [track_y[0]],
        color="black",
        alpha=0.7,
        s=100,
        zorder=6,
        edgecolors="white",
        linewidths=2,
        label="TC Position",
    )

    # Add legend
    ax.legend(loc="upper right", fontsize=10, framealpha=0.9)

    print("  TC track added successfully")


def update_plot(
    i,
    da_h,
    cax_h,
    ax,
    tc_track_gdf=None,
    tc_track_times=None,
    tc_track_dot=None,
    data_crs=None,
):
    """Update function for animation."""
    da_hi = da_h.isel(time=i)
    t = str(da_hi.time.values)[:19]  # Format timestamp

    # Try to format nicely if it's a datetime
    try:
        t = da_hi.time.dt.strftime("%d-%B-%Y %H:%M:%S").item()
    except:
        pass

    ax.set_title(f"SFINCS water depth {t}")

    # Update the data
    # Get the QuadMesh object and update its array
    cax_h.set_array(da_hi.values.ravel())

    # Update TC track dot position if available
    if (
        tc_track_dot is not None
        and tc_track_times is not None
        and tc_track_gdf is not None
    ):
        # Get current flood timestep
        current_time = pd.to_datetime(da_hi.time.values)

        # Find closest TC track point in time
        time_diffs = (tc_track_times - current_time).abs()
        closest_idx = time_diffs.argmin()

        # Get the coordinates of the closest track point in the data CRS
        tc_track_plot = tc_track_gdf.to_crs(data_crs)
        closest_x = tc_track_plot.geometry.x.values[closest_idx]
        closest_y = tc_track_plot.geometry.y.values[closest_idx]

        # Update dot position
        tc_track_dot.set_offsets([[closest_x, closest_y]])

    return [cax_h]


# Create animation
print(f"Generating {len(range(0, da_h.time.size, STEP))} frames...")
ani = animation.FuncAnimation(
    fig,
    update_plot,
    frames=np.arange(0, da_h.time.size, STEP),
    interval=INTERVAL,
    fargs=(da_h, cax_h, ax, tc_track_gdf, tc_track_times, tc_track_dot, data_crs),
    blit=False,
)

# ===== SAVE ANIMATION =====
output_file = (
    OUTPUT_DIR / f"sfincs_{EVENT_NAME.lower()}_{SCENARIO}_flood_animation2.mp4"
)
print(f"Saving animation to: {output_file}")
print("This may take a few minutes...")

try:
    ani.save(str(output_file), fps=FPS, dpi=DPI, writer="ffmpeg")
    print(f"Animation saved successfully: {output_file}")
except Exception as e:
    print(f"Error saving animation: {e}")
    print("Make sure ffmpeg is installed: conda install ffmpeg -c conda-forge")

    # Try saving as GIF as fallback
    output_file_gif = (
        OUTPUT_DIR / f"sfincs_{EVENT_NAME.lower()}_{SCENARIO}_flood_animation.gif"
    )
    try:
        print(f"Trying to save as GIF instead: {output_file_gif}")
        ani.save(str(output_file_gif), fps=FPS, dpi=DPI // 2, writer="pillow")
        print(f"GIF saved successfully: {output_file_gif}")
    except Exception as e2:
        print(f"Could not save as GIF either: {e2}")

# Close the dataset
ds.close()

plt.close()

print(f"\nAnimation complete for {EVENT_NAME}!")
print(f"Output file: {output_file}")

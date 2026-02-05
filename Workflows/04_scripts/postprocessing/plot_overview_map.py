"""
Plot overview map showing the region of analysis in context.

This script creates a zoomed-out map with an OSM basemap showing:
1. The broader geographic context (country/regional level)
2. The region of analysis highlighted as a semi-transparent polygon

The region is read from the SFINCS model's gis/region.geojson file.

by @dumontgoulart
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import contextily as ctx
from pathlib import Path
import warnings
import xarray as xr
import numpy as np
from shapely.geometry import box
from shapely.ops import unary_union
from rasterio import features
from affine import Affine

warnings.filterwarnings("ignore")

# ===== CONFIGURATION =====
# Set your event name here
EVENT_NAME = "Freddy"  # Change this to: "Kenneth", "Freddy", "Idai", etc.

# Base paths - update these as needed
BASE_SFINCS_PATH = Path("/p/11210471-001-compass/03_Runs")
OUTPUT_DIR = Path("/p/11210471-001-compass/04_Results/CF_figs")

# Map styling options
ZOOM_BUFFER = 4.0  # Buffer around region in degrees (adjust for zoom level)
REGION_ALPHA = 0.4  # Transparency of region polygon (0-1)
REGION_COLOR = "red"  # Color of region polygon
REGION_EDGE_COLOR = "darkred"  # Edge color of region polygon
REGION_EDGE_WIDTH = 2  # Edge width of region polygon
BASEMAP_ZOOM = 10  # Basemap zoom level (lower = more zoomed out)

# ===== DYNAMIC FILE PATHS =====
# Map event names to their SFINCS run folders
# Update this dictionary as needed for your events
EVENT_TO_SFINCS_FOLDER = {
    "Freddy": "sfincs_Freddy2",
    "Kenneth": "sfincs_Kenneth",
    "Idai": "sfincs_Idai",
}

# Map event names to their TC track files
TC_TRACKS_BASE = Path("/p/11210471-001-compass/01_Data/IBTrACS/SELECTED_TRACKS")
EVENT_TO_TC_TRACKS = {
    "Freddy": [
        TC_TRACKS_BASE / "IBTrACS_FREDDY_part1.shp",
        TC_TRACKS_BASE / "IBTrACS_FREDDY_part2.shp",
    ],
    "Kenneth": [TC_TRACKS_BASE / "IBTrACS_KENNETH.shp"],
    "Idai": [TC_TRACKS_BASE / "IBTrACS_IDAI.shp"],
}

# Get SFINCS folder name (use mapping or default to sfincs_{EVENT_NAME})
sfincs_folder = EVENT_TO_SFINCS_FOLDER.get(EVENT_NAME, f"sfincs_{EVENT_NAME}")

# Construct path to staticmaps.nc from wflow
# Use the staticmaps.nc file from wflow to define the region extent
BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs/test")
staticmaps_file = (
    BASE_RUN_PATH
    / EVENT_NAME
    / "wflow"
    / "event_precip_era5_hourly_CF0"
    / "staticmaps.nc"
)

# Get TC track files for this event
tc_track_files = EVENT_TO_TC_TRACKS.get(EVENT_NAME, [])

# Set region_file to staticmaps_file (used by create_overview_map)
region_file = staticmaps_file

# Create output directory if it doesn't exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Creating overview map for event: {EVENT_NAME}")
print(f"Staticmaps file: {staticmaps_file}")


def create_overview_map(
    region_path: Path,
    event_name: str,
    output_dir: Path,
    zoom_buffer: float = 0.5,
    region_alpha: float = 0.4,
    region_color: str = "red",
    region_edge_color: str = "darkred",
    region_edge_width: float = 2,
    basemap_zoom: int = 10,
    tc_track_files: list = None,
):
    """
    Create a zoomed-out overview map showing the region of analysis.

    Parameters
    ----------
    region_path : Path
        Path to the region.geojson file
    event_name : str
        Name of the event (used for title and filename)
    output_dir : Path
        Directory to save the output figure
    zoom_buffer : float
        Buffer around region bounds in degrees
    region_alpha : float
        Transparency of region polygon (0-1)
    region_color : str
        Fill color of region polygon
    region_edge_color : str
        Edge color of region polygon
    region_edge_width : float
        Edge width of region polygon
    basemap_zoom : int
        Zoom level for basemap tiles
    tc_track_files : list
        List of paths to TC track shapefiles
    """
    # ===== LOAD REGION DATA FROM STATICMAPS.NC =====
    print(f"Loading region extent from: {region_path}")

    if not region_path.exists():
        raise FileNotFoundError(f"Staticmaps file not found: {region_path}")

    # Load the netCDF file and extract extent from dem_subgrid
    ds = xr.open_dataset(region_path)

    # Get the spatial coordinates - check multiple possible names
    if "x" in ds.coords and "y" in ds.coords:
        x_coords = ds["x"].values
        y_coords = ds["y"].values
    elif "lon" in ds.coords and "lat" in ds.coords:
        x_coords = ds["lon"].values
        y_coords = ds["lat"].values
    elif "longitude" in ds.coords and "latitude" in ds.coords:
        x_coords = ds["longitude"].values
        y_coords = ds["latitude"].values
    elif "longitude" in ds.dims and "latitude" in ds.dims:
        x_coords = ds["longitude"].values
        y_coords = ds["latitude"].values
    else:
        raise ValueError(
            f"Could not find x/y, lon/lat, or longitude/latitude coordinates in {region_path}. Available coords: {list(ds.coords)}, dims: {list(ds.dims)}"
        )

    # ===== CREATE APPROXIMATE DEM SILHOUETTE =====
    print("Creating approximate DEM silhouette...")

    # Try to load dem_subgrid to get the valid data mask
    if "dem_subgrid" in ds.data_vars:
        dem = ds["dem_subgrid"].values
        print(f"DEM shape: {dem.shape}")
    elif "wflow_dem" in ds.data_vars:
        dem = ds["wflow_dem"].values
        print(f"DEM shape: {dem.shape}")
    else:
        # Fallback to bounding box if no DEM variable found
        print("No DEM variable found, using bounding box")
        dem = None

    # Get the bounding box (needed for both approaches)
    minx, maxx = float(x_coords.min()), float(x_coords.max())
    miny, maxy = float(y_coords.min()), float(y_coords.max())

    # Determine CRS
    if abs(minx) <= 180 and abs(maxx) <= 180 and abs(miny) <= 90 and abs(maxy) <= 90:
        data_crs = "EPSG:4326"
        print("Detected CRS: EPSG:4326 (WGS84)")
    else:
        if "crs" in ds.attrs:
            data_crs = ds.attrs["crs"]
        elif hasattr(ds, "spatial_ref"):
            data_crs = ds.spatial_ref.attrs.get("crs_wkt", "EPSG:32737")
        else:
            data_crs = "EPSG:32737"
        print(f"Using CRS: {data_crs}")

    if dem is not None:
        # Create a valid data mask (where DEM is not NaN and not fill value)
        # Coarsen the mask for faster processing (use every Nth pixel)
        coarsen_factor = max(
            1, min(dem.shape) // 100
        )  # Aim for ~100 pixels on smallest dimension
        print(f"Coarsening by factor {coarsen_factor} for faster processing")

        # Create binary mask of valid data
        valid_mask = ~np.isnan(dem) & (dem != -9999) & (dem != -9999.0)

        # Coarsen the mask
        if coarsen_factor > 1:
            # Use simple slicing for coarsening
            valid_mask_coarse = valid_mask[::coarsen_factor, ::coarsen_factor]
            x_coarse = x_coords[::coarsen_factor]
            y_coarse = y_coords[::coarsen_factor]
        else:
            valid_mask_coarse = valid_mask
            x_coarse = x_coords
            y_coarse = y_coords

        print(f"Coarsened mask shape: {valid_mask_coarse.shape}")

        # Create transform for rasterio
        dx = abs(x_coarse[1] - x_coarse[0]) if len(x_coarse) > 1 else 1
        dy = abs(y_coarse[1] - y_coarse[0]) if len(y_coarse) > 1 else 1

        # Determine if y is ascending or descending
        y_ascending = y_coarse[0] < y_coarse[-1] if len(y_coarse) > 1 else True

        if y_ascending:
            transform = Affine.translation(
                x_coarse.min() - dx / 2, y_coarse.max() + dy / 2
            ) * Affine.scale(dx, -dy)
        else:
            transform = Affine.translation(
                x_coarse.min() - dx / 2, y_coarse.min() - dy / 2
            ) * Affine.scale(dx, dy)

        # Convert mask to polygons using rasterio.features
        try:
            shapes = list(
                features.shapes(
                    valid_mask_coarse.astype(np.uint8),
                    mask=valid_mask_coarse,
                    transform=transform,
                )
            )

            if shapes:
                # Create polygons from shapes (only where mask is True)
                from shapely.geometry import shape

                polygons = [shape(geom) for geom, value in shapes if value == 1]

                if polygons:
                    # Merge all polygons into one
                    region_polygon = unary_union(polygons)

                    # Simplify the polygon for faster rendering
                    # Use tolerance based on pixel size
                    simplify_tolerance = max(dx, dy) * 2
                    region_polygon = region_polygon.simplify(
                        simplify_tolerance, preserve_topology=True
                    )

                    print(
                        f"Created simplified DEM silhouette with {region_polygon.geom_type}"
                    )
                else:
                    print("No valid polygons found, using bounding box")
                    region_polygon = box(minx, miny, maxx, maxy)
            else:
                print("No shapes extracted, using bounding box")
                region_polygon = box(minx, miny, maxx, maxy)

        except Exception as e:
            print(f"Error creating silhouette: {e}, using bounding box")
            region_polygon = box(minx, miny, maxx, maxy)
    else:
        # Fallback to bounding box
        region_polygon = box(minx, miny, maxx, maxy)

    # Create a GeoDataFrame with the region polygon
    gdf_region = gpd.GeoDataFrame({"geometry": [region_polygon]}, crs=data_crs)

    print(f"Region CRS: {gdf_region.crs}")
    print(f"Region bounds: {gdf_region.total_bounds}")

    ds.close()

    # Convert to WGS84 (lat/lon) for display
    gdf_region_ll = gdf_region.to_crs("EPSG:4326")

    # ===== CALCULATE MAP EXTENT =====
    # Get bounds in lat/lon
    minx, miny, maxx, maxy = gdf_region_ll.total_bounds
    print(f"Original bounds: x=[{minx:.4f}, {maxx:.4f}], y=[{miny:.4f}, {maxy:.4f}]")

    # Apply buffer in degrees
    minx_buf = minx - zoom_buffer
    maxx_buf = maxx + zoom_buffer
    miny_buf = miny - zoom_buffer
    maxy_buf = maxy + zoom_buffer
    print(
        f"Buffered bounds: x=[{minx_buf:.4f}, {maxx_buf:.4f}], y=[{miny_buf:.4f}, {maxy_buf:.4f}]"
    )

    # ===== CREATE FIGURE =====
    # Calculate aspect ratio for figure sizing
    data_width = maxx_buf - minx_buf
    data_height = maxy_buf - miny_buf
    aspect_ratio = data_width / data_height
    print(
        f"Data extent: width={data_width:.4f}, height={data_height:.4f}, aspect={aspect_ratio:.4f}"
    )

    # Set figure size based on aspect ratio to match data aspect
    fig_height = 8
    fig_width = fig_height * aspect_ratio
    print(f"Figure size: {fig_width:.2f} x {fig_height}")

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)

    # ===== SET MAP EXTENT FIRST =====
    # Set extent BEFORE adding basemap so contextily knows what tiles to fetch
    ax.set_xlim(minx_buf, maxx_buf)
    ax.set_ylim(miny_buf, maxy_buf)

    # ===== ADD BASEMAP =====
    print("Adding basemap...")
    try:
        ctx.add_basemap(
            ax=ax,
            source=ctx.providers.OpenStreetMap.Mapnik,
            crs="EPSG:4326",
            attribution=False,
            zorder=1,
            zoom=basemap_zoom,
        )
        print("Basemap added successfully")
    except Exception as e:
        print(f"Warning: Could not add basemap: {e}")

    # ===== PLOT REGION ON TOP =====
    gdf_region_ll.plot(
        ax=ax,
        facecolor=region_color,
        edgecolor=region_edge_color,
        alpha=region_alpha,
        linewidth=region_edge_width,
        zorder=3,
    )

    # ===== PLOT TC TRACKS =====
    if tc_track_files:
        print(f"Loading TC track files for {event_name}...")
        for i, track_file in enumerate(tc_track_files):
            if track_file.exists():
                print(f"  Loading: {track_file.name}")
                try:
                    gdf_track = gpd.read_file(track_file)
                    # Convert to EPSG:4326 if needed
                    if gdf_track.crs != "EPSG:4326":
                        gdf_track = gdf_track.to_crs("EPSG:4326")

                    # Plot the track line connecting the points
                    # Extract coordinates for line plot
                    track_x = gdf_track.geometry.x.values
                    track_y = gdf_track.geometry.y.values

                    # Plot thin connecting line
                    ax.plot(
                        track_x,
                        track_y,
                        color="black",
                        linewidth=1,
                        zorder=4,
                        alpha=0.6,
                    )

                    # Plot the dots on top
                    gdf_track.plot(
                        ax=ax,
                        color="black",
                        alpha=0.4,
                        markersize=60,
                        zorder=5,
                        label=(
                            f"TC Track {i+1}" if len(tc_track_files) > 1 else "TC Track"
                        ),
                    )
                    print(f"  ✓ Plotted TC track: {track_file.name}")
                except Exception as e:
                    print(f"  ✗ Could not load TC track {track_file.name}: {e}")
            else:
                print(f"  ✗ TC track file not found: {track_file}")
    else:
        print("No TC track files specified")

    # ===== STYLING =====
    # Don't set aspect - figure size already matches data aspect ratio
    # Setting aspect='equal' with constrained_layout can cause extent issues
    ax.set_title(
        f"Region of Analysis: {event_name}",
        fontsize=14,
        fontweight="bold",
        pad=10,
    )

    # Add axis labels
    ax.set_xlabel("Longitude [°]", fontsize=10)
    ax.set_ylabel("Latitude [°]", fontsize=10)

    # ===== ADD LEGEND =====
    legend_elements = [
        mpatches.Patch(
            facecolor=region_color,
            edgecolor=region_edge_color,
            alpha=region_alpha,
            linewidth=region_edge_width,
            label="Region of Analysis",
        )
    ]

    # Add TC track to legend if plotted
    if tc_track_files and any(f.exists() for f in tc_track_files):
        from matplotlib.lines import Line2D

        legend_elements.append(
            Line2D([0], [0], color="black", linewidth=2, label="TC Track")
        )

    ax.legend(
        handles=legend_elements,
        loc="upper right",
        fontsize=10,
        frameon=True,
        facecolor="white",
        edgecolor="gray",
        framealpha=0.9,
    )

    # ===== SAVE FIGURE =====
    output_file = output_dir / f"{event_name}_overview_map.png"
    fig.savefig(output_file, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Overview map saved to: {output_file}")

    plt.close(fig)

    return output_file


def create_mozambique_africa_overview(
    event_name: str,
    output_dir: Path,
    region_alpha: float = 0.4,
    region_color: str = "red",
    region_edge_color: str = "darkred",
    region_edge_width: float = 2,
    basemap_zoom: int = 5,
    borders_file: Path = Path(
        "/p/11210471-001-compass/01_Data/admin_borders/CNTR_RG_60M_2024_4326.gpkg"
    ),
):
    """
    Create an overview map showing Mozambique highlighted within Africa.

    Parameters
    ----------
    event_name : str
        Name of the event (used for filename)
    output_dir : Path
        Directory to save the output figure
    region_alpha : float
        Transparency of Mozambique polygon (0-1)
    region_color : str
        Fill color of Mozambique polygon
    region_edge_color : str
        Edge color of Mozambique polygon
    region_edge_width : float
        Edge width of Mozambique polygon
    basemap_zoom : int
        Zoom level for basemap tiles
    borders_file : Path
        Path to the country borders geopackage file
    """
    print("\n" + "=" * 60)
    print("Creating Mozambique-Africa overview map")
    print("=" * 60)

    # Load country borders and extract Mozambique
    print(f"Loading country borders from: {borders_file}")
    gdf_borders = gpd.read_file(borders_file)
    print(f"Loaded {len(gdf_borders)} country polygons")

    # Filter for Mozambique
    gdf_mozambique = gdf_borders[gdf_borders["NAME_ENGL"] == "Mozambique"].copy()

    if len(gdf_mozambique) == 0:
        print("ERROR: Mozambique not found in borders file!")
        print(
            f"Available countries: {sorted(gdf_borders['NAME_ENGL'].unique())[:10]}..."
        )
        raise ValueError("Mozambique not found in borders file")

    print(f"Found Mozambique polygon")

    # Ensure it's in EPSG:4326
    if gdf_mozambique.crs != "EPSG:4326":
        gdf_mozambique = gdf_mozambique.to_crs("EPSG:4326")

    # Get Mozambique bounds
    moz_bounds = gdf_mozambique.total_bounds
    print(
        f"Mozambique bounds: x=[{moz_bounds[0]:.2f}, {moz_bounds[2]:.2f}], y=[{moz_bounds[1]:.2f}, {moz_bounds[3]:.2f}]"
    )

    # Africa extent with some buffer
    # Showing southern and eastern Africa prominently
    africa_minx, africa_maxx = 10, 52
    africa_miny, africa_maxy = -35, 5

    print(
        f"Africa extent: x=[{africa_minx}, {africa_maxx}], y=[{africa_miny}, {africa_maxy}]"
    )

    # ===== CREATE FIGURE =====
    data_width = africa_maxx - africa_minx
    data_height = africa_maxy - africa_miny
    aspect_ratio = data_width / data_height
    print(
        f"Map extent: width={data_width:.1f}°, height={data_height:.1f}°, aspect={aspect_ratio:.2f}"
    )

    fig_height = 10
    fig_width = fig_height * aspect_ratio
    print(f"Figure size: {fig_width:.2f} x {fig_height}")

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)

    # ===== SET MAP EXTENT FIRST =====
    ax.set_xlim(africa_minx, africa_maxx)
    ax.set_ylim(africa_miny, africa_maxy)

    # ===== ADD BASEMAP =====
    print("Adding basemap...")
    try:
        ctx.add_basemap(
            ax=ax,
            source=ctx.providers.OpenStreetMap.Mapnik,
            crs="EPSG:4326",
            attribution=False,
            zorder=1,
            zoom=basemap_zoom,
        )
        print("Basemap added successfully")
    except Exception as e:
        print(f"Warning: Could not add basemap: {e}")

    # ===== PLOT MOZAMBIQUE HIGHLIGHT =====
    # Plot the actual Mozambique polygon
    gdf_mozambique.plot(
        ax=ax,
        facecolor=region_color,
        edgecolor=region_edge_color,
        alpha=region_alpha,
        linewidth=region_edge_width,
        zorder=3,
    )

    # ===== STYLING =====
    ax.set_title(
        f"Regional Context: Mozambique - {event_name}",
        fontsize=16,
        fontweight="bold",
        pad=15,
    )

    ax.set_xlabel("Longitude [°]", fontsize=12)
    ax.set_ylabel("Latitude [°]", fontsize=12)

    # ===== ADD LEGEND =====
    legend_patch = mpatches.Patch(
        facecolor=region_color,
        edgecolor=region_edge_color,
        alpha=region_alpha,
        linewidth=region_edge_width,
        label="Mozambique",
    )
    ax.legend(
        handles=[legend_patch],
        loc="upper left",
        fontsize=11,
        frameon=True,
        facecolor="white",
        edgecolor="gray",
        framealpha=0.9,
    )

    # ===== ADD NORTH ARROW =====
    arrow_x, arrow_y = 0.95, 0.15
    arrow_length = 0.08

    ax.annotate(
        "",
        xy=(arrow_x, arrow_y + arrow_length),
        xytext=(arrow_x, arrow_y),
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="->", color="black", lw=2.5),
    )
    ax.text(
        arrow_x,
        arrow_y + arrow_length + 0.02,
        "N",
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        ha="center",
        va="bottom",
    )

    # ===== SAVE FIGURE =====
    output_file = output_dir / f"{event_name}_mozambique_africa_overview.png"
    fig.savefig(output_file, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Mozambique-Africa overview map saved to: {output_file}")

    plt.close(fig)

    return output_file


# ===== MAIN EXECUTION =====
if __name__ == "__main__":
    try:
        # Create regional overview map (zoomed into analysis region)
        output_file1 = create_overview_map(
            region_path=region_file,
            event_name=EVENT_NAME,
            output_dir=OUTPUT_DIR,
            zoom_buffer=ZOOM_BUFFER,
            region_alpha=REGION_ALPHA,
            region_color=REGION_COLOR,
            region_edge_color=REGION_EDGE_COLOR,
            region_edge_width=REGION_EDGE_WIDTH,
            basemap_zoom=BASEMAP_ZOOM,
            tc_track_files=tc_track_files,
        )
        print(f"\n✓ Regional overview map created: {output_file1}")

        # Create Mozambique-Africa overview map
        output_file2 = create_mozambique_africa_overview(
            event_name=EVENT_NAME,
            output_dir=OUTPUT_DIR,
            region_alpha=REGION_ALPHA,
            region_color=REGION_COLOR,
            region_edge_color=REGION_EDGE_COLOR,
            region_edge_width=REGION_EDGE_WIDTH,
            basemap_zoom=5,  # More zoomed out for Africa view
        )
        print(f"✓ Mozambique-Africa overview map created: {output_file2}")

        print("\n" + "=" * 60)
        print("All overview maps created successfully!")
        print("=" * 60)

    except FileNotFoundError as e:
        print(f"\nError: {e}")
        print("Please check the EVENT_NAME and BASE_SFINCS_PATH configuration.")
        print(
            "Make sure the region.geojson file exists in the SFINCS model's gis folder."
        )

    except Exception as e:
        print(f"\nUnexpected error: {e}")
        import traceback

        traceback.print_exc()

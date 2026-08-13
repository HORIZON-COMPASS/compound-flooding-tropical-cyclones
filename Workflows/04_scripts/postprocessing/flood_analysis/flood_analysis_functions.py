"""
Functions for flood_analysis_script.py
"""

import itertools

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import contextily as ctx
import matplotlib.pyplot as plt
import numpy as np
import rioxarray as rxr  # Required for reading TIFF files
import xarray as xr
from cartopy.io.img_tiles import OSM
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from pyproj import Transformer
from shapely import contains_xy


def add_city_markers(ax, extent, cities, fontsize=8, markersize=50):
    """
    Add city markers to a map if they fall within the plot extent.

    Parameters:
    ax (matplotlib axis): The axis to plot on
    extent (tuple): (minx, maxx, miny, maxy) extent of the plot in EPSG:4326
    cities (list): List of city dictionaries with 'name', 'lat', 'lon' keys
    fontsize (int): Font size for city labels
    markersize (int): Size of city markers

    Returns:
    int: Number of cities plotted
    """
    minx, maxx, miny, maxy = extent

    cities_in_extent = []
    for city in cities:
        if minx <= city["lon"] <= maxx and miny <= city["lat"] <= maxy:
            cities_in_extent.append(city)

    if not cities_in_extent:
        return 0

    # Plot city markers
    for city in cities_in_extent:
        ax.scatter(
            city["lon"],
            city["lat"],
            s=markersize,
            c="red",
            marker="o",
            edgecolors="white",
            linewidths=1.5,
            zorder=10,
            alpha=0.8,
        )

        # Add city label with background for better visibility
        # Offset label slightly up and to the right of the marker
        label_offset_x = 0.005  # Degrees longitude
        label_offset_y = 0.005  # Degrees latitude
        ax.text(
            city["lon"] + label_offset_x,
            city["lat"] + label_offset_y,
            city["name"],
            fontsize=fontsize,
            ha="left",
            va="bottom",
            zorder=11,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7, edgecolor="none"),
        )

    return len(cities_in_extent)


def _degree_ticks(vmin, vmax, n_ticks=5):
    """
    Create evenly spaced degree ticks for axis labelling.

    Parameters:
    vmin (float): Minimum value of the axis
    vmax (float): Maximum value of the axis
    n_ticks (int): Number of ticks to generate

    Returns:
    np.ndarray: Array of tick values
    """
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        return np.array([])
    if np.isclose(vmin, vmax):
        delta = 0.1
        return np.array([vmin - delta, vmin, vmax + delta])
    
    return np.round(np.linspace(vmin, vmax, n_ticks), 2)


def format_lonlat_axes(ax, extent, label_fontsize=10, tick_fontsize=9):
    """
    Apply longitude/latitude ticks and degree-format labels to a map axis.
    
    Parameters:
    ax (matplotlib axis): The axis to format
    extent (tuple): (minx, maxx, miny, maxy) extent of the plot in EPSG:4326
    label_fontsize (int): Font size for axis labels
    tick_fontsize (int): Font size for tick labels
    """
    minx, maxx, miny, maxy = extent
    ax.set_extent([minx, maxx, miny, maxy], crs=ccrs.PlateCarree())

    x_ticks = _degree_ticks(minx, maxx)
    y_ticks = _degree_ticks(miny, maxy)

    if x_ticks.size:
        ax.set_xticks(x_ticks, crs=ccrs.PlateCarree())
    if y_ticks.size:
        ax.set_yticks(y_ticks, crs=ccrs.PlateCarree())

    ax.xaxis.set_major_formatter(LongitudeFormatter(number_format=".2f", degree_symbol="°"))
    ax.yaxis.set_major_formatter(LatitudeFormatter(number_format=".2f", degree_symbol="°"))
    ax.tick_params(axis="both", labelsize=tick_fontsize)
    ax.set_xlabel("Longitude", fontsize=label_fontsize)
    ax.set_ylabel("Latitude", fontsize=label_fontsize)


def load_file(scenario, file):
    """
    Load a preprocessed flood file for a given scenario.

    Parameters:
    scenario (str): The scenario name e.g. "factual"
    file (pathlib.Path): Path to the .tif file

    Returns:
    zsmax (xarray.DataArray): The loaded flood data
    original_crs (pyproj.CRS): The original coordinate reference system of the data
    """
    if not file.exists():
        raise FileNotFoundError(f"Required file not found for scenario '{scenario}': {file}")

    print(f"Loading preprocessed flood file for scenario '{scenario}'...")
    zsmax = rxr.open_rasterio(file)

    if "band" in zsmax.dims:
        zsmax = zsmax.squeeze("band", drop=True)

    original_crs = zsmax.rio.crs
    print(f"Original CRS: {original_crs}")

    return zsmax, original_crs


def crop_to_lonlat_domain(da, domain):
    """
    Crop the data to the given domain.

    If the data is in a projected CRS (metres), the lat/lon domain bounds are
    reprojected to match before clipping.  If already in EPSG:4326, direct
    coordinate selection is used with automatic handling of axis ordering.

    Parameters:
    da (xarray.DataArray): The data to crop
    domain (dict): Dictionary with keys 'lat_min', 'lat_max', 'lon_min', 'lon_max'

    Returns:
    xarray.DataArray: The cropped data
    """
    crs = da.rio.crs
    if crs is not None and str(crs) != "EPSG:4326":
        # Reproject the four corners of the bounding box into the data CRS.
        transformer = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
        xmin, ymin = transformer.transform(domain["lon_min"], domain["lat_min"])
        xmax, ymax = transformer.transform(domain["lon_max"], domain["lat_max"])
        return da.rio.clip_box(minx=xmin, miny=ymin, maxx=xmax, maxy=ymax)

    # Already in lat/lon — use .sel() with ordering-aware slices.
    x_name = "x" if "x" in da.coords else "lon"
    y_name = "y" if "y" in da.coords else "lat"

    x_vals = da.coords[x_name].values
    y_vals = da.coords[y_name].values

    # Respect coordinate ordering (ascending vs descending).
    x_slice = (
        slice(domain["lon_min"], domain["lon_max"])
        if x_vals[0] <= x_vals[-1]
        else slice(domain["lon_max"], domain["lon_min"])
    )
    y_slice = (
        slice(domain["lat_min"], domain["lat_max"])
        if y_vals[0] <= y_vals[-1]
        else slice(domain["lat_max"], domain["lat_min"])
    )

    return da.sel({x_name: x_slice, y_name: y_slice})


def apply_shapefile(zsmax, shapefile):
    """
    Apply a shapefile mask to the given data array.

    Parameters:
    zsmax (xarray.DataArray): The data to mask
    shapefile (str): Path to the shapefile

    Returns:
    xarray.DataArray: The masked data
    """
    # load in shape file for cropping
    area_reader = shpreader.Reader(shapefile)
    flood_mask_shp = [record.geometry for record in area_reader.records()]

    # reproject shape file
    x, y = np.meshgrid(zsmax.x, zsmax.y)
    lon_lat_area = ccrs.TransverseMercator(
        central_longitude=-2.0,
        central_latitude=49.0,
        false_easting=400000.0,
        false_northing=-100000.0,
        scale_factor=0.9996012717,
    ).transform_points(zsmax.rio.crs, x, y)
    print(f"lon_lat shape: {lon_lat_area.shape}")

    masks_area = [
        contains_xy(flood_mask_shp[i], lon_lat_area[:, :, 0], lon_lat_area[:, :, 1])
        for i in range(len(flood_mask_shp))
    ]
    area_mask = np.logical_or.reduce(masks_area)

    # Save the mask as a dataset
    mask_ds = xr.DataArray(
        area_mask,
        name="flood_area",
        coords={"lat": lon_lat_area[:, 0, 1], "lon": lon_lat_area[0, :, 0]},
    )

    zsmax = zsmax.where(mask_ds.values)

    return zsmax


def handle_nans(zsmax_dict):
    """
    Set NaNs values to 0 for difference calculations in original CRS, but only where other scenario(s) has valid data.

    Parameters:
    zsmax_dict (dict): Dictionary of xarray.DataArray objects for each scenario

    Returns:
    zsmax_calc_dict (dict): Dictionary of xarray.DataArray objects with NaNs handled
    """
    mask_valids = {}
    for key in zsmax_dict:
        mask_valid = ~np.isnan(zsmax_dict[key])
        mask_valids[key] = mask_valid

    if "factual" in zsmax_dict:
        zsmax_factual = zsmax_dict["factual"]
        mask_factual_valid = mask_valids["factual"]
    if "counterfactual" in zsmax_dict:
        zsmax_counterfactual = zsmax_dict["counterfactual"]
        mask_counterfactual_valid = mask_valids["counterfactual"]
    if "future" in zsmax_dict:
        zsmax_future = zsmax_dict["future"]
        mask_future_valid = mask_valids["future"]

    zsmax_calc_dict = {}

    if len(zsmax_dict) == 3:  # all three scenarios
        zsmax_factual_calc = zsmax_factual.where(
            mask_factual_valid | ~mask_counterfactual_valid | ~mask_future_valid, 0
        )
        zsmax_counterfactual_calc = zsmax_counterfactual.where(
            mask_counterfactual_valid | ~mask_factual_valid | ~mask_future_valid, 0
        )
        zsmax_future_calc = zsmax_future.where(
            mask_future_valid | ~mask_factual_valid | ~mask_counterfactual_valid, 0
        )
        zsmax_calc_dict["factual"] = zsmax_factual_calc
        zsmax_calc_dict["counterfactual"] = zsmax_counterfactual_calc
        zsmax_calc_dict["future"] = zsmax_future_calc

    elif len(zsmax_dict) == 2:
        if "factual" in zsmax_dict and "counterfactual" in zsmax_dict:
            zsmax_factual_calc = zsmax_factual.where(
                mask_factual_valid | ~mask_counterfactual_valid, 0
            )
            zsmax_counterfactual_calc = zsmax_counterfactual.where(
                mask_counterfactual_valid | ~mask_factual_valid, 0
            )
            zsmax_calc_dict["factual"] = zsmax_factual_calc
            zsmax_calc_dict["counterfactual"] = zsmax_counterfactual_calc
        elif "factual" in zsmax_dict and "future" in zsmax_dict:
            zsmax_factual_calc = zsmax_factual.where(mask_factual_valid | ~mask_future_valid, 0)
            zsmax_future_calc = zsmax_future.where(mask_future_valid | ~mask_factual_valid, 0)
            zsmax_calc_dict["factual"] = zsmax_factual_calc
            zsmax_calc_dict["future"] = zsmax_future_calc
        elif "counterfactual" in zsmax_dict and "future" in zsmax_dict:
            zsmax_counterfactual_calc = zsmax_counterfactual.where(
                mask_counterfactual_valid | ~mask_future_valid, 0
            )
            zsmax_future_calc = zsmax_future.where(
                mask_future_valid | ~mask_counterfactual_valid, 0
            )
            zsmax_calc_dict["counterfactual"] = zsmax_counterfactual_calc
            zsmax_calc_dict["future"] = zsmax_future_calc

    elif len(zsmax_dict) == 1:
        print(
            "No other scenarios available to compare with for NaN handling, replacing NaNs with 0s."
        )
        key = next(iter(zsmax_dict))
        zsmax_calc = zsmax_dict[key].where(mask_valids[key], 0)
        zsmax_calc_dict[key] = zsmax_calc

    return zsmax_calc_dict


def get_cell_size(zsmax_dict):
    """
    Get cell size in meters from the original projected CRS.

    Parameters:
    zsmax_dict (dict): Dictionary of xarray.DataArray objects for each scenario

    Returns:
    cell_area_m2 (float): Cell area in square meters
    """
    try:
        any_zsmax = next(iter(zsmax_dict.values()))
        dx = float(np.abs(any_zsmax.x.diff("x").median()))
        dy = float(np.abs(any_zsmax.y.diff("y").median()))
        cell_area_m2 = dx * dy
        print(f"Cell size: {dx:.2f} x {dy:.2f} m = {cell_area_m2:.2f} m²")
    except Exception:
        cell_area_m2 = np.nan
        print("Warning: Could not determine cell size")

    return cell_area_m2


def pairwise_extent(pairs, zsmax_calc_dict, area_m2, FLOOD_THRESHOLD):
    """
    Pairwise flood extent metrics, keyed by "{b}_vs_{a}"

    Parameters:
    pairs (list of tuples): List of scenario pairs to compare
    zsmax_calc_dict (dict): Dictionary of xarray.DataArray objects for each scenario
    area_m2 (dict): Dictionary of flood extent areas for each scenario
    FLOOD_THRESHOLD (float): Threshold for considering a cell as flooded

    Returns:
    diff_calc_dict (dict): Difference in flood depth between scenarios
    mean_diff_flooded_union_dict (dict): Mean depth difference over the union of flooded areas between scenarios
    percent_diff_dict (dict): Percent difference in flood depth between scenarios
    extent_area_diff_m2_dict (dict): Difference in flood extent area between scenarios
    extent_area_pct_dict (dict): Percent difference in flood extent area between scenarios
    """
    diff_calc_dict = {}
    flood_mask_union_dict = {}
    mean_diff_flooded_union_dict = {}
    percent_diff_dict = {}
    extent_area_diff_m2_dict = {}
    extent_area_pct_dict = {}

    for a, b in pairs:
        key = f"{b}_vs_{a}"
        diff_calc_dict[key] = zsmax_calc_dict[b] - zsmax_calc_dict[a]
        flood_mask_union_dict[key] = (zsmax_calc_dict[a] > FLOOD_THRESHOLD) | (
            zsmax_calc_dict[b] > FLOOD_THRESHOLD
        )
        mean_diff_flooded_union_dict[key] = (
            diff_calc_dict[key].where(flood_mask_union_dict[key]).mean(skipna=True)
        )
        percent_diff_dict[key] = xr.where(
            zsmax_calc_dict[a] > 0,
            (diff_calc_dict[key] / zsmax_calc_dict[a]) * 100.0,
            np.nan,
        )
        extent_area_diff_m2_dict[key] = area_m2[b] - area_m2[a]
        extent_area_pct_dict[key] = (
            (extent_area_diff_m2_dict[key] / area_m2[a] * 100.0) if area_m2[a] > 0 else np.nan
        )

    return (
        diff_calc_dict,
        mean_diff_flooded_union_dict,
        percent_diff_dict,
        extent_area_diff_m2_dict,
        extent_area_pct_dict,
    )


def pairwise_volume(pairs, volume_m3):
    """
    Pairwise flood volume metrics, keyed by "{b}_vs_{a}"

    Parameters:
    pairs (list of tuples): List of scenario pairs to compare
    volume_m3 (dict): Dictionary of flood volumes for each scenario

    Returns:
    volume_diff_m3_dict (dict): Difference in flood volume between scenarios
    volume_pct_dict (dict): Percent difference in flood volume between scenarios
    """
    volume_diff_m3_dict = {}
    volume_pct_dict = {}

    for a, b in pairs:
        key = f"{b}_vs_{a}"
        volume_diff_m3_dict[key] = volume_m3[b] - volume_m3[a]
        volume_pct_dict[key] = (
            (volume_diff_m3_dict[key] / volume_m3[a] * 100.0) if volume_m3[a] > 0 else np.nan
        )

    return volume_diff_m3_dict, volume_pct_dict


def determine_coord_system(zsmax_calc_dict):
    """
    Determine the UTM zone and hemisphere (northern/southern) based on the data coordinates.

    Parameters:
    zsmax_calc_dict (dict): Dictionary of flood depth data arrays for each scenario

    Returns:
    utm_zone (int): UTM zone number
    southern (bool): True if southern hemisphere, False if northern hemisphere
    """
    try:
        # Get approximate center coordinates from an entry in zsmax_calc_dict
        center_x = float(zsmax_calc_dict[list(zsmax_calc_dict.keys())[0]].x.mean())
        center_y = float(zsmax_calc_dict[list(zsmax_calc_dict.keys())[0]].y.mean())

        # Guess UTM zone based on coordinates (this is approximate)
        if center_x > 300000 and center_x < 800000:  # Typical UTM coordinate range
            if center_y < 0:  # Southern hemisphere
                utm_zone = 37  # Adjust based on your region
                southern = True
            else:
                utm_zone = 37  # Adjust based on your region
                southern = False
        else:
            utm_zone = 37  # Default for your region
            southern = True

    except Exception as e:
        print(f"Could not determine UTM zone: {e}")
        utm_zone = 37
        southern = True

    return utm_zone, southern


def main_comparison_plot(zsmax_plot_dict, diff_calc_dict, scenarios, pairs, 
                         event_city_list, event_title, event_name, output_dir, 
                         flood_levels, flood_cmap, diff_levels, diff_cmap):
    """
    Create a single figure containing individual scenario maps and pairwise difference maps.

    Parameters:
    zsmax_plot_dict (dict): Dictionary of flood depth data arrays for each scenario
    diff_calc_dict (dict): Dictionary of difference data arrays for each scenario pair
    scenarios (list): List of scenario names
    pairs (list of tuples): List of scenario pairs to compare
    event_city_list (list): List of city dictionaries with 'name', 'lat', 'lon' keys
    event_title (str): Title of the event for the figure
    event_name (str): Name of the event for output file naming
    output_dir (pathlib.Path): Directory to save the output figure
    flood_levels (np.ndarray): Levels for flood depth plotting
    flood_cmap (matplotlib colormap): Colormap for flood depth plotting
    diff_levels (np.ndarray): Levels for difference plotting
    diff_cmap (matplotlib colormap): Colormap for difference plotting

    Returns:
    extent (tuple): (minx, maxx, miny, maxy) extent of the plot in EPSG:4326
    output_file_main (pathlib.Path): Path to the saved output figure
    """
    panel_specs = []
    for scenario in scenarios:
        panel_specs.append(
            {
                "data": zsmax_plot_dict[scenario],
                "title": scenario.capitalize(),
                "is_difference": False,
            }
        )

    for a, b in pairs:
        key = f"{b}_vs_{a}"
        if key in diff_calc_dict:
            panel_specs.append(
                {
                    "data": diff_calc_dict[key],
                    "title": f"Flood Depth Changes: {b.capitalize()} vs {a.capitalize()}",
                    "is_difference": True,
                }
            )

    n_panels = len(panel_specs)
    ncols = min(3, max(1, n_panels))
    nrows = int(np.ceil(n_panels / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(6.5 * ncols, 4.5 * nrows),  subplot_kw={"projection": ccrs.PlateCarree()})
    axes = np.atleast_1d(axes).ravel()

    first_scenario = scenarios[0]
    extent = (
        zsmax_plot_dict[first_scenario].x.min().item(),
        zsmax_plot_dict[first_scenario].x.max().item(),
        zsmax_plot_dict[first_scenario].y.min().item(),
        zsmax_plot_dict[first_scenario].y.max().item(),
    )
    print(
        f"Plot extent: lon=[{extent[0]:.2f}, {extent[1]:.2f}], lat=[{extent[2]:.2f}, {extent[3]:.2f}]"
    )

    for idx, spec in enumerate(panel_specs):
        ax = axes[idx]
        im = spec["data"].plot(
            ax=ax,
            levels=diff_levels if spec["is_difference"] else flood_levels,
            cmap=diff_cmap if spec["is_difference"] else flood_cmap,
            add_colorbar=False,
            x="x",
            y="y",
            alpha=0.7,
            zorder=2,
        )
        imagery = OSM()
        ax.add_image(imagery, 12)
        ax.set_aspect("equal")
        ax.set_title(spec["title"], fontsize=12, fontweight="bold", pad=10)

        add_city_markers(ax, extent, event_city_list, fontsize=7, markersize=40)
        format_lonlat_axes(ax, extent, label_fontsize=10, tick_fontsize=8)

        cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.1)
        cbar.set_label(
            "Depth Difference [m]" if spec["is_difference"] else "Max Flood Depth [m]",
            fontsize=10,
        )

    for ax in axes[n_panels:]:
        ax.axis("off")

    plt.tight_layout(h_pad=3.0, w_pad=2.0)
    plt.suptitle(
        f"Flood Depth Analysis - {event_title}",
        fontsize=16,
        fontweight="bold",
        y=1.02,
    )

    output_file_main = output_dir / f"sfincs_{event_name.lower()}_flood_depth_comparison.png"
    plt.savefig(output_file_main, dpi=400, bbox_inches="tight")
    plt.close()

    return extent, output_file_main


def individual_plot(scenario, zsmax_plot_dict, extent, EVENT_CITY_LIST, EVENT_TITLE, 
                    EVENT_NAME, OUTPUT_DIR, flood_levels, flood_cmap, output_file_individuals):
    """
    Create an individual flood depth map for a given scenario.

    Parameters:
    scenario (str): The scenario name e.g. "factual"
    zsmax_plot_dict (dict): Dictionary of flood depth data arrays for each scenario
    extent (tuple): (minx, maxx, miny, maxy) extent of the plot in EPSG:4326
    EVENT_CITY_LIST (list): List of city dictionaries with 'name', 'lat', 'lon' keys
    EVENT_TITLE (str): Title of the event for the figure
    EVENT_NAME (str): Name of the event for output file naming
    OUTPUT_DIR (pathlib.Path): Directory to save the output figure
    flood_levels (np.ndarray): Levels for flood depth plotting
    flood_cmap (matplotlib colormap): Colormap for flood depth plotting
    output_file_individuals (list): List to store output file paths of individual scenario plots

    Returns:
    output_file_individuals (list): Updated list of output file paths of individual scenario plots
    """
    print(f"Creating individual {scenario} plot...")

    zsmax_plot = zsmax_plot_dict.get(scenario, None)
    if zsmax_plot is None:
        print(f"No data available for scenario '{scenario}'. Skipping plot.")
        return

    fig, ax = plt.subplots(1, 1, figsize=(12, 8),  subplot_kw={"projection": ccrs.PlateCarree()})

    im_plot = zsmax_plot.plot(
        ax=ax,
        levels=flood_levels,
        cmap=flood_cmap,
        add_colorbar=False,
        x="x",
        y="y",
        alpha=0.7,
        zorder=2,
    )
    imagery = OSM()
    ax.add_image(imagery, 12)
    ax.set_aspect("equal")

    # Add city markers
    add_city_markers(ax, extent, EVENT_CITY_LIST, fontsize=9, markersize=60)

    # Format axis labels and ticks as lon/lat degrees
    format_lonlat_axes(ax, extent, label_fontsize=12, tick_fontsize=10)

    # Colorbar
    cbar = plt.colorbar(im_plot, ax=ax, shrink=0.8, pad=0.1)
    cbar.set_label("Max Flood Depth [m]", fontsize=12)

    # Title
    ax.set_title(
        f"{scenario.capitalize()} Flood Depth - {EVENT_TITLE}",
        fontsize=14,
        fontweight="bold",
        pad=20,
    )

    # Save
    output_file = OUTPUT_DIR / f"sfincs_{EVENT_NAME.lower()}_{scenario}.png"
    output_file_individuals.append(output_file)
    plt.savefig(output_file, dpi=400, bbox_inches="tight")
    plt.close()

    return output_file_individuals


def difference_plot(a, b, diff_calc_dict, extent, EVENT_CITY_LIST, EVENT_TITLE, 
                    EVENT_NAME, OUTPUT_DIR, diff_levels, diff_cmap, output_file_diffs):
    """
    Create a depth difference plot for a given pair of scenarios.

    Parameters:
    a (str): The first scenario name e.g. "factual"
    b (str): The second scenario name e.g. "counterfactual"
    diff_calc_dict (dict): Dictionary of difference data arrays for each scenario pair
    extent (tuple): (minx, maxx, miny, maxy) extent of the plot in EPSG:4326
    EVENT_CITY_LIST (list): List of city dictionaries with 'name', 'lat', 'lon' keys
    EVENT_TITLE (str): Title of the event for the figure
    EVENT_NAME (str): Name of the event for output file naming
    OUTPUT_DIR (pathlib.Path): Directory to save the output figure
    diff_levels (np.ndarray): Levels for difference plotting
    diff_cmap (matplotlib colormap): Colormap for difference plotting
    output_file_diffs (list): List to store output file paths of difference plots

    Returns:
    output_file_diffs (list): Updated list of output file paths of difference plots
    """
    key = f"{b}_vs_{a}"
    if key not in diff_calc_dict:
        return None

    print(f"Creating detailed {b.capitalize()}-{a.capitalize()} difference plot...")
    fig_diff, ax = plt.subplots(1, 1, figsize=(12, 8),  subplot_kw={"projection": ccrs.PlateCarree()})

    im_diff = diff_calc_dict[key].plot(
        ax=ax,
        levels=diff_levels,
        cmap=diff_cmap,
        add_colorbar=False,
        x="x",
        y="y",
        alpha=0.7,
        zorder=2,
    )
    imagery = OSM()
    ax.add_image(imagery, 12)
    ax.set_aspect("equal")

    add_city_markers(ax, extent, EVENT_CITY_LIST, fontsize=9, markersize=60)
    format_lonlat_axes(ax, extent, label_fontsize=12, tick_fontsize=10)

    cbar = plt.colorbar(im_diff, ax=ax, shrink=0.8, pad=0.1)
    cbar.set_label("Flood Depth Difference [m]", fontsize=12)

    ax.set_title(
        f"Flood Depth Changes - {EVENT_TITLE}: {b.capitalize()} vs {a.capitalize()}",
        fontsize=14,
        fontweight="bold",
        pad=20,
    )

    output_file_diff = OUTPUT_DIR / f"sfincs_{EVENT_NAME.lower()}_difference_{b}-{a}.png"
    output_file_diffs.append(output_file_diff)
    plt.savefig(output_file_diff, dpi=400, bbox_inches="tight")
    plt.close()

    return output_file_diffs


def create_bar_chart_with_annotation(ax, scenarios, totals, ylabel, title, unit_prefix="", unit_suffix=""):
    """
    Helper function to create a bar chart with climate change attribution annotation.

    Parameters:
    ax (matplotlib axis): The axis to plot on
    scenarios (list): List of scenario names
    totals (list): List of total values corresponding to each scenario
    ylabel (str): Label for the y-axis
    title (str): Title of the bar chart
    unit_prefix (str): Prefix for the units in the annotation (default: "")
    unit_suffix (str): Suffix for the units in the annotation (default: "")
    """
    if len(scenarios) != len(totals):
        raise ValueError("'scenarios' and 'totals' must have the same length.")

    bar_colors = ["steelblue"] * len(scenarios)
    bars = ax.bar(scenarios, totals, color=bar_colors, alpha=1, width=0.4)

    ax.set_ylabel(ylabel, fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_axisbelow(True)

    # Annotate all pairwise changes in the provided scenario order.
    annotation_pairs = list(itertools.combinations(range(len(totals)), 2))

    max_total = max(totals) if totals else 0
    y_pad = max(max_total * 0.06, 0.02)

    for i, (idx_a, idx_b) in enumerate(annotation_pairs):
        a = totals[idx_a]
        b = totals[idx_b]
        diff_val = b - a
        pct = (diff_val / a * 100.0) if a != 0 else np.nan
        label = f"{scenarios[idx_b]} - {scenarios[idx_a]}"

        x0 = bars[idx_a].get_x() + bars[idx_a].get_width() / 2
        x1 = bars[idx_b].get_x() + bars[idx_b].get_width() / 2
        y0 = bars[idx_a].get_height()
        y1 = bars[idx_b].get_height()
        y_top = max(y0, y1) + y_pad * (1.1 + i * 1.3)

        # Bracket line linking the compared bars.
        ax.plot(
            [x0, x0, x1, x1],
            [y_top - y_pad * 0.15, y_top, y_top, y_top - y_pad * 0.15],
            color="black",
            linewidth=1.2,
        )

        label_text = (
            f"{label}: {unit_prefix}{diff_val:+.2f}{unit_suffix} ({pct:+.1f}%)"
            if np.isfinite(pct)
            else f"{label}: {unit_prefix}{diff_val:+.2f}{unit_suffix}"
        )

        ax.text(
            (x0 + x1) / 2,
            y_top + y_pad * 0.1,
            label_text,
            ha="center",
            va="bottom",
            fontsize=9,
        )

    # Set y-axis to start from 0 with appropriate margin
    extra_top = y_pad * (1.8 + 1.3 * len(annotation_pairs))
    ax.set_ylim(0, max_total + extra_top)

    # Get bar positions
    bar_positions = [bar.get_x() + bar.get_width() / 2 for bar in bars]

    # Annotation parameters for vertical difference line
    if len(bar_positions) >= 2:
        x_span = bar_positions[-1] - bar_positions[0]
        line_x_position = bar_positions[-1] + max(x_span * 0.3, 0.2)
    elif len(bar_positions) == 1:
        line_x_position = bar_positions[0] + 0.2
    else:
        line_x_position = ax.get_xlim()[1]

    # Add extended dashed horizontal lines
    left_edge = ax.get_xlim()[0]
    for total in totals:
        ax.plot(
            [left_edge, line_x_position],
            [total, total],
            linestyle="--",
            color="gray",
            alpha=0.7,
            linewidth=1,
        )

    # Remove top and right spines for cleaner look
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def bar_charts(SCENARIOS, area_m2, volume_m3, EVENT_TITLE, EVENT_NAME, OUTPUT_DIR):
    """
    Create bar charts for flood extent and volume for each scenario. Annotates with 
    pairwise comparisons for attribution results.

    Parameters:
    SCENARIOS (list): List of scenario names
    area_m2 (dict): Dictionary of flood extent in m² for each scenario
    volume_m3 (dict): Dictionary of flood volume in m³ for each scenario
    EVENT_TITLE (str): Title of the event
    EVENT_NAME (str): Name of the event
    OUTPUT_DIR (Path): Directory to save the output bar chart

    Returns:
    output_file_agg_bar (Path): Path to the saved bar chart image
    """
    fig_bars, (ax_extent, ax_volume) = plt.subplots(1, 2, figsize=(14, 6))

    scenarios = [scenario.capitalize() for scenario in SCENARIOS]
    extent_totals = [
        area_m2.get(scenario, np.nan) / 1e6 for scenario in SCENARIOS
    ]  # Convert to km²
    volume_totals = [
        volume_m3.get(scenario, np.nan) / 1e6 for scenario in SCENARIOS
    ]  # Convert to Mm³

    create_bar_chart_with_annotation(
        ax_extent,
        scenarios,
        extent_totals,
        ylabel="Flood Extent [km²]",
        title=f"(a) Flood Extent - {EVENT_TITLE}",
        unit_prefix="",
        unit_suffix=" km²",
    )

    create_bar_chart_with_annotation(
        ax_volume,
        scenarios,
        volume_totals,
        ylabel="Flood Volume [Mm³]",
        title=f"(b) Flood Volume - {EVENT_TITLE}",
        unit_prefix="",
        unit_suffix=" Mm³",
    )

    plt.tight_layout()
    output_file_agg_bar = (
        OUTPUT_DIR / f"sfincs_{EVENT_NAME.lower()}_flood_extent_volume_barchart.png"
    )
    plt.savefig(output_file_agg_bar, dpi=400, bbox_inches="tight")
    plt.close()

    return output_file_agg_bar

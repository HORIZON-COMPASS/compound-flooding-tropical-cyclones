"""
This script processes SFINCS output files (sfincs_output_hmax_AllTime.tif) and compares model runs 
for attribution studies. It will calculate flood extent, volume, and depth metrics, as well as 
pairwise differences between scenarios provided, and output the results to a CSV.

It will also create plots of the flood depth for each scenario and the differences between 
scenarios, and create bar charts of flood extent and volume.

To use, update flood_analysis_config.py with parameters for your event. Also update the 
configuration section below with the output directory, scenarios to process, flood depth threshold, 
and event name.

Based on plot_counterfactual_flood_changes.py from compound_flooding_tropical_cyclones repository, 
but updated to be more flexible and handle any number of scenarios for a given event.
"""

import itertools
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from flood_analysis_config import (EVENT_CITIES, EVENT_CONFIG, diff_cmap,
                                   diff_levels, flood_cmap, flood_levels)
from flood_analysis_functions import (apply_shapefile, bar_charts,
                                      crop_to_lonlat_domain,
                                      determine_coord_system, difference_plot,
                                      get_cell_size, handle_nans,
                                      individual_plot, load_file,
                                      main_comparison_plot, pairwise_extent,
                                      pairwise_volume)

warnings.filterwarnings("ignore")


# ===== CONFIGURATION =====
# Base path - update as needed
OUTPUT_DIR = Path(
    "/data/users/ukcr.compass/COMPASS/somerset-compass/storyline_results_jul26/reformatted_cf_flood_changes_output"
)

# Depth threshold over which to consider flood extent
FLOOD_THRESHOLD = 0.3  # meters

# Scenarios running for - will analyse each individually and do pairwise comparisons
SCENARIOS = ["counterfactual", "factual", "future"]
# The code will handle any number of scenarios, but list order matters for pairwise comparisons
# Should be ordered as past -> future, earlier scenarios will be subtracted from later ones (b - a)
pairs = list(itertools.combinations(SCENARIOS, 2))
# This means you could do e.g. two counterfactual and one factual, or one factual and two future
# e.g. if SCENARIOS = ["counterfactual", "factual", "future"], then pairwise comparisons will be:
# factual - counterfactual, future - factual, future - counterfactual
# Make sure scenario names match keys for scenarios' data folders in EVENT_CONFIG dictionary

# Event to select from config
EVENT_NAME = "Somerset_Winter201314"


# ===== MAIN SCRIPT =====
print(f"Processing event: {EVENT_NAME}")

# Create output directory if it doesn't exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Validate event name
if EVENT_NAME not in EVENT_CONFIG:
    raise ValueError(
        f"Unknown event: {EVENT_NAME}. Valid options: {list(EVENT_CONFIG.keys())}"
    )

# Get event-specific configuration
event_cfg = EVENT_CONFIG[EVENT_NAME]
BASE_RUN_PATH = event_cfg["base_path"]  # Get base path where event data is stored
EVENT_TITLE = event_cfg["title"]  # Get event title for plots
EVENT_DOMAIN = event_cfg["domain"]  # Get domain for this event
EVENT_SHAPEFILE = event_cfg["shapefile"]  # Get shapefile path for this event, if being used
EVENT_CITY_LIST = EVENT_CITIES.get(EVENT_NAME, [])  # Get cities for this event

# Gather files for scenarios for event being processed
files_to_process = []

for scenario in SCENARIOS:
    if scenario not in event_cfg:
        print(
            f"Warning: Scenario '{scenario}' not defined for event '{EVENT_NAME}' in config. Skipping."
        )
        continue

    file = BASE_RUN_PATH / event_cfg[scenario] / "sfincs_output_hmax_AllTime.tif"
    files_to_process.append((scenario, file))

if len(files_to_process) == 0:
    raise ValueError(
        "No valid scenarios selected for processing. Please check SCENARIOS list."
    )

# Load and process each scenario file
zsmax_dict = {}
for scenario, file in files_to_process:
    zsmax, original_crs = load_file(scenario, file)
    zsmax_cropped = crop_to_lonlat_domain(zsmax, EVENT_DOMAIN)
    if EVENT_SHAPEFILE.exists():
        zsmax_cropped = apply_shapefile(zsmax_cropped, EVENT_SHAPEFILE)
    zsmax_dict[scenario] = zsmax_cropped

cell_area_m2 = get_cell_size(zsmax_dict)
zsmax_calc_dict = handle_nans(zsmax_dict)

# Calculate area/volume before any reprojection
print("Calculating flood metrics in original projected CRS...")

# ===== FLOOD EXTENT =====
# Build flood masks and areas
flood_mask_dict = {
    scenario: zsmax_calc_dict[scenario] > FLOOD_THRESHOLD for scenario in SCENARIOS
}
# Flood extent area per scenario (m²)
area_m2 = {
    scenario: float(flood_mask_dict[scenario].sum().values) * cell_area_m2
    for scenario in SCENARIOS
}

if len(SCENARIOS) > 1:
    (   diff_calc_dict,
        mean_diff_flooded_union_dict,
        percent_diff_dict,
        extent_area_diff_m2_dict,
        extent_area_pct_dict,
    ) = pairwise_extent(pairs, zsmax_calc_dict, area_m2, FLOOD_THRESHOLD)
else:
    print("Only one scenario present — skipping pairwise comparisons.")

# ===== FLOOD VOLUME =====
# Flood volume (sum depth * cell area over flooded cells)
if not np.isnan(cell_area_m2):
    volume_m3 = {
        scenario: float(
            zsmax_calc_dict[scenario].where(flood_mask_dict[scenario]).sum(skipna=True)
            * cell_area_m2
        )
        for scenario in SCENARIOS
    }
    if len(SCENARIOS) > 1:
        volume_diff_m3_dict, volume_pct_dict = pairwise_volume(pairs, volume_m3)
else:
    print("Warning: Could not determine cell size, skipping volume calculations.")
    volume_m3 = {scenario: np.nan for scenario in SCENARIOS}

# ===== REPROJECT TO LAT/LON FOR LATER PLOTTING =====
# Only reproject AFTER extent and volume metrics are calculated
if original_crs != "EPSG:4326":
    print("Reprojecting to EPSG:4326 (lat/lon) for plotting...")
    for scenario in SCENARIOS:
        zsmax_calc_dict[scenario] = zsmax_calc_dict[scenario].rio.reproject("EPSG:4326")
    for a, b in pairs:
        key = f"{b}_vs_{a}"
        diff_calc_dict[key] = diff_calc_dict[key].rio.reproject("EPSG:4326")
    print("Reprojection complete")
else:
    print("Already in EPSG:4326")

# ===== FLOOD DEPTH =====
max_depth = {scenario: float(zsmax_calc_dict[scenario].max()) for scenario in SCENARIOS}
mean_depth_diff_all = {
    f"{b}_vs_{a}": float(diff_calc_dict[f"{b}_vs_{a}"].mean()) for a, b in pairs
}
mean_depth_diff_union = {
    f"{b}_vs_{a}": float(mean_diff_flooded_union_dict[f"{b}_vs_{a}"]) for a, b in pairs
}

# ===== PRINT SUMMARY OF FLOOD METRICS =====
print("\n================ SUMMARY ================")
print(f"Event: {EVENT_NAME}")

print(f"\n[1] Flood Extent (threshold > {FLOOD_THRESHOLD} m)")
for scenario in SCENARIOS:
    print(f"  {scenario.capitalize()} area: {area_m2[scenario]/1e6:.3f} km^2")
for a, b in pairs:
    key = f"{b}_vs_{a}"
    print(f"{key}")
    print(f"  Absolute change:     {extent_area_diff_m2_dict[key]/1e6:.3f} km^2")
    print(f"  Percent change:      {extent_area_pct_dict[key]:.2f} %")

print("\n[2] Flood Depth")
for scenario in SCENARIOS:
    print(f"  {scenario.capitalize()} max depth: {max_depth[scenario]:.3f} m")
for a, b in pairs:
    key = f"{b}_vs_{a}"
    print(f"{key}")
    print(f"  Mean depth diff (all cells):        {mean_depth_diff_all[key]:.3f} m")
    print(f"  Mean depth diff (flooded union):    {mean_depth_diff_union[key]:.3f} m")
    print(
        f"  Mean percent depth change (where {a} > 0): {percent_diff_dict[key].mean(skipna=True).values:.2f} %"
    )

print(f"\n[3] Flood Volume (threshold > {FLOOD_THRESHOLD} m)")
for scenario in SCENARIOS:
    print(f"  {scenario.capitalize()} volume: {volume_m3[scenario]/1e6:.3f} Mm^3")
for a, b in pairs:
    key = f"{b}_vs_{a}"
    print(f"{key}")
    print(f"  Absolute change:     {volume_diff_m3_dict[key]/1e6:.3f} Mm^3")
    print(f"  Percent change:      {volume_pct_dict[key]:.2f} %")
print("=========================================\n")

# ===== EXPORT AGGREGATED FLOOD DATA TO CSV =====
print("Exporting aggregated flood data to CSV...")

# Create DataFrame rows from scenario and pairwise keys
csv_rows = []
volume_diff_m3_for_csv = locals().get("volume_diff_m3_dict", {})
volume_pct_for_csv = locals().get("volume_pct_dict", {})

# Scenario totals
for scenario in SCENARIOS:
    csv_rows.append(
        {
            "Scenario": scenario.capitalize(),
            "Flood_Extent_km2": area_m2.get(scenario, np.nan) / 1e6,
            "Flood_Volume_Mm3": volume_m3.get(scenario, np.nan) / 1e6,
        }
    )

# Pairwise rows: one for absolute difference, one for percent change
for a, b in pairs:
    key = f"{b}_vs_{a}"
    csv_rows.append(
        {
            "Scenario": f"{key}_difference",
            "Flood_Extent_km2": extent_area_diff_m2_dict.get(key, np.nan) / 1e6,
            "Flood_Volume_Mm3": volume_diff_m3_for_csv.get(key, np.nan) / 1e6,
        }
    )
    csv_rows.append(
        {
            "Scenario": f"{key}_percent_change",
            "Flood_Extent_km2": extent_area_pct_dict.get(key, np.nan),
            "Flood_Volume_Mm3": volume_pct_for_csv.get(key, np.nan),
        }
    )

flood_data = pd.DataFrame(csv_rows)

# Add metadata columns
flood_data["Event"] = EVENT_NAME
flood_data["Flood_Threshold_m"] = FLOOD_THRESHOLD

# Reorder columns
flood_data = flood_data[
    ["Event", "Scenario", "Flood_Extent_km2", "Flood_Volume_Mm3", "Flood_Threshold_m"]
]

# Save to CSV
output_file_csv = (
    OUTPUT_DIR / f"sfincs_{EVENT_NAME.lower()}_flood_aggregated_data_reformatted.csv"
)
flood_data.to_csv(output_file_csv, index=False)
print(f"Aggregated flood data saved to: {output_file_csv}")
print(flood_data.to_string(index=False))


# ===== PLOTTING =====
zsmax_plot_dict = {}
for scenario in SCENARIOS:
    zsmax = zsmax_calc_dict.get(scenario, None)
    # Only show areas with positive water depth (above ground)
    zsmax_plot = zsmax.where(zsmax > 0.05)
    zsmax_plot_dict[scenario] = zsmax_plot

# Data is now in EPSG:4326
data_crs = "EPSG:4326"
print(f"Creating plots in {data_crs}...")
utm_zone, southern = determine_coord_system(zsmax_calc_dict)
print(f"Using UTM zone {utm_zone}, Southern: {southern}")

extent, output_file_main = main_comparison_plot(
    zsmax_plot_dict,
    diff_calc_dict,
    SCENARIOS,
    pairs,
    EVENT_CITY_LIST,
    EVENT_TITLE,
    EVENT_NAME,
    OUTPUT_DIR,
    flood_levels,
    flood_cmap,
    diff_levels,
    diff_cmap,
)

output_file_individuals = []
for scenario in SCENARIOS:
    individual_plot(
        scenario,
        zsmax_plot_dict,
        extent,
        EVENT_CITY_LIST,
        EVENT_TITLE,
        EVENT_NAME,
        OUTPUT_DIR,
        flood_levels,
        flood_cmap,
        output_file_individuals
    )

output_file_diffs = []
for a, b in pairs:
    difference_plot(
        a,
        b,
        diff_calc_dict,
        extent,
        EVENT_CITY_LIST,
        EVENT_TITLE,
        EVENT_NAME,
        OUTPUT_DIR,
        diff_levels,
        diff_cmap,
        output_file_diffs,
    )

print("Creating aggregated flood metrics bar charts...")
output_file_agg_bar = bar_charts(
    SCENARIOS, area_m2, volume_m3, EVENT_TITLE, EVENT_NAME, OUTPUT_DIR
)


# ===== CLEANUP =====
print(
    f"\nAnalysis complete for {EVENT_NAME}! Check the saved PNG files in: {OUTPUT_DIR}"
)
print("Files created:")
print(f"  - {output_file_main}")
for output_file_individual in output_file_individuals:
    print(f"  - {output_file_individual}")
for output_file_diff in output_file_diffs:
    print(f"  - {output_file_diff}")
print(f"  - {output_file_agg_bar}")
print(f"  - {output_file_csv}")

for ds in zsmax_dict.values():
    try:
        ds.close()
    except Exception:
        pass

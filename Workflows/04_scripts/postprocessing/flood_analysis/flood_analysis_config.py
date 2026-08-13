"""
Configuration file for flood analysis and plotting.
EVENT_CONFIG is used for mapping event names to their parameters e.g. 
specific folder paths where sfincs_output_hmax_AllTime.tif is located.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Event cities and domain for plots
# Somerset - domain matched to that of precip inputs
somerset_domain = {
    "lat_min": 50.9,
    "lat_max": 51.2,
    "lon_min": -3.1,
    "lon_max": -2.6,
}
somerset_cities = [
    {"name": "Taunton", "lat": 51.0153, "lon": -3.1060},
    {"name": "Bridgwater", "lat": 51.1250, "lon": -3.0000},
    {"name": "Glastonbury", "lat": 51.1470, "lon": -2.7170},
]

# Map event name to specific folder paths and base directory
EVENT_CONFIG = {
    "Somerset_Winter201314": {
        "base_path": Path(
            "/data/users/ukcr.compass/COMPASS/somerset-compass/storyline_results_jul26"
        ),
        "factual": "factual",  # subfolder name for each scenario provided in SCENARIOS list - keys must match scenario names
        "counterfactual": "t1",
        "future": "t3",
        "title": "Somerset Winter 2013/14",
        "domain": somerset_domain,
        "shapefile": Path(
            "/data/users/ukcr.compass/COMPASS/somerset-compass/storyline_results_jul26/new_shape/new_shape2.shp"
        ),  # leave blank if not using shapefile
    },
}

# Map events to their city lists
EVENT_CITIES = {
    "Somerset_Winter201314": somerset_cities,
}

# Colormaps and levels for plots
flood_levels = np.arange(0, 5, 0.5)
flood_cmap = plt.cm.Blues

# Difference levels (starting from 0 for Blues colormap)
dif_min = 0
dif_max = 0.1
diff_levels = np.linspace(dif_min, dif_max, 11)
diff_cmap = plt.cm.Blues

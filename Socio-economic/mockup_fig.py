#%%
import geopandas as gpd
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import BoundaryNorm, ListedColormap
import warnings
warnings.filterwarnings('ignore')
import platform
from pathlib import Path
import rasterio
from rasterio.mask import mask
from shapely.geometry import box as shapely_box

prefix = "p:/" if platform.system() == "Windows" else "/p/"

# ===== CONFIGURATION =====
EVENT_NAME = "Idai"
BASE_RUN_PATH = Path("C:/Code/Paper_1/Data_submission")
SCENARIO_PATH = "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" # factual


#%% ######################################################################
############################### POPULATION ###############################
##########################################################################
HE_2020 = "data/gdf_grid_HE_2020.gpkg"
HE_1990 ="data/gdf_grid_HE_1990.gpkg"

# Read the geodataframes
gdf_HE_2020 = gpd.read_file(HE_2020)
gdf_HE_1990 = gpd.read_file(HE_1990)


# %%
# Statistics
total_2020_HE = gdf_HE_2020["population"].sum()
total_1990_HE = gdf_HE_1990["population"].sum()
total_diff = gdf_HE_1990['population_diff'].sum()
perc_diff = (gdf_HE_2020['population'].sum() - gdf_HE_1990['population_diff'].sum()) / gdf_HE_2020['population'].sum() * 100

print(f"Total exposed population 2020 (HE): {total_2020_HE:.0f} people")
print(f"Total exposed population 1990 (HE): {total_1990_HE:.0f} people")
print(f"Attributable exposed population to population change: {total_diff:.0f} people")
print(f"Attributable exposed population to population change: {perc_diff:.0f} %")




# %%
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

# Example categories
categories = ["Climate change", "GDP change", "FA change", "Population change"]

# Example values (two scenarios per category)
scenario1 = [30, 47, 25, 0]
scenario2 = [20, 0, 0, 15]
colors1 = ["#1E2E57", "#67CBE4", "#67CBE4",  "#67CBE4"]
colors2 = ["#495B88", "#67CBE4", "#67CBE4",  "#67CBE4"]

x = np.arange(len(categories))  # category positions
width = 0.35  # bar width

fig, ax1 = plt.subplots(figsize=(6,4))

# Bars
ax2 = ax1.twinx()
ax1.bar(x - width/2, scenario1, width, color=colors1, edgecolor="black")
ax2.bar(x + width/2, scenario2, width, color=colors2, edgecolor="black", hatch="//")

# Formatting
ax1.set_ylabel("Attributable damage [%]")
ax1.set_xlabel("Driver")
ax1.set_title("Mock-up bar plot")
ax1.set_xticks(x)
ax1.set_xticklabels(categories)

damage_patch = mpatches.Patch(facecolor="white", edgecolor="black", label="Damage")
population_patch = mpatches.Patch(facecolor="white", edgecolor="black", hatch="//", label="Affected population")
ax1.legend(handles=[damage_patch, population_patch], loc="upper right")

ax2.set_ylabel("Attributable affected population [%]")
plt.tight_layout()
plt.show()

"""
Configuration for the socio-economic attribution analysis.
All paths, constants, and event-specific settings are centralized here.
"""
import os
import platform
from pathlib import Path

import numpy as np

# ===== PLATFORM =====
IS_WINDOWS = platform.system() == "Windows"
PREFIX = "p:/" if IS_WINDOWS else "/p/"

# ===== EVENT CONFIGURATION =====
EVENT_NAME = "Idai"
BASE_RUN_PATH = Path("p:/11210471-001-compass/03_Runs/sofala/Idai")
SCENARIO_PATH_F = "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0"
SCENARIO_PATH_CF = "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.1_era5_hourly_spw_IBTrACS_CF-5"

# ===== SFINCS DIRECTORIES =====
SFINCS_DIR_F = BASE_RUN_PATH / "sfincs" / SCENARIO_PATH_F
SFINCS_DIR_CF = BASE_RUN_PATH / "sfincs" / SCENARIO_PATH_CF

# ===== DATA CATALOG =====
if IS_WINDOWS:
    DATACAT_PATH = os.path.abspath("../../Workflows/03_data_catalogs/datacatalog_general.yml")
else:
    DATACAT_PATH = os.path.abspath("../../Workflows/03_data_catalogs/datacatalog_general___linux.yml")

# ===== INPUT FILE PATHS =====
REGION_PATH = SFINCS_DIR_F / "gis" / "region.geojson"
BACKGROUND_PATH = Path("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_region_background.geojson")
SOFALA_PROVINCE_PATH = Path("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_province.shp")
BEIRA_DISTRICT_PATH = Path("../../Attribution_results/data/gis/Beira_region.shp")
DISTRICTS_ADM3_PATH = Path("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_districts_study_region.shp")

POPULATION_2019_PATH = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_2019_30.tif")
POPULATION_1990_PATH = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_1990_30.tif")
SETTLEMENT_TYPE_PATH = Path("C:/Code/COMPASS/Socio-economic/results/gis/avg_rural_per_grid.tif")

F_FLOODING_PATH = SFINCS_DIR_F / "plot_output" / "floodmap.tif"
CF_FLOODING_PATH = SFINCS_DIR_CF / "plot_output" / "floodmap.tif"
SFINCS_SUBGRID_PATH = SFINCS_DIR_F / "subgrid" / "dep_subgrid.tif"

# ===== OUTPUT PATHS =====
EXPORT_PATH = Path("p:/11210471-001-compass/04_Results/Idai_socioeconomic/preprocessed/population/")
SUMMARY_TABLE_PATH = Path("c:/Code/Paper_2/Output/flood_depth_impact_summary_table.csv")
DISTRICT_SUMMARY_PATH = Path("c:/Code/COMPASS_exposure/Data/Modified/sofala_district_exposed_population_summary.csv")

# ===== REGRIDDED POPULATION OUTPUT PATHS =====
def regridded_pop_path(year):
    return Path(f"c:/Code/COMPASS_exposure/Data/Modified/population_{year}_region_regrid_wholepeople.tif")

# ===== DISTRICTS TO DROP =====
DROP_DISTRICTS = ["Muanza", "Gororngosa-Sede", "Galinha"]

# ===== PLOTTING CONSTANTS =====
COLOURS = ['#00B050', '#1E2E57', "#28C2E9", '#9B59B6']

DEPTH_BINS = {"0-0.5 m": (0, 0.5), "0.5-1.5 m": (0.5, 1.5), ">1.5 m": (1.5, np.inf)}

# Mask polygon for removing incorrect coastline boundary
MASK_POLY_COORDS = [(34.9, -20.3), (36, -20.3), (36, -19.9), (34.9, -19.9)]

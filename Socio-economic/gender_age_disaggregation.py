#%%
import os
import json
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
import rasterio
import geopandas as gpd
import warnings
warnings.filterwarnings('ignore')
import platform
import matplotlib.pyplot as plt
from rasterio.warp import reproject, Resampling
from shapely.geometry import Polygon
import rioxarray as rxr 
from hydromt import DataCatalog
import re
import glob


prefix = "p:/" if platform.system() == "Windows" else "/p/"

#%%
# ===== CONFIGURATION =====
EVENT_NAME = "Idai"
BASE_RUN_PATH = Path("C:/Code/Paper_1/Data_submission")
SCENARIO_PATH_F = "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0" # factual
SCENARIO_PATH_CF = "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10" # counterfactual

# ===== FILE PATHS =====
# Base directory for the specific event and scenario
fiat_dir_F    = BASE_RUN_PATH / "fiat" / SCENARIO_PATH_F
fiat_dir_CF   = BASE_RUN_PATH / "fiat" / SCENARIO_PATH_CF
sfincs_dir_F  = BASE_RUN_PATH / "sfincs" / SCENARIO_PATH_F
sfincs_dir_CF = BASE_RUN_PATH / "sfincs" / SCENARIO_PATH_CF

# ===== DATA CATALOG =====
if platform.system() == "Windows":
    datacat_path = os.path.abspath("../Workflows/03_data_catalogs/datacatalog_general.yml")
else:
    datacat_path = os.path.abspath("../Workflows/03_data_catalogs/datacatalog_general___linux.yml")
data_catalog = DataCatalog(data_libs = [datacat_path])

#%%
# Input files
shapefile_fp = "p:/11210471-001-compass/03_Runs/sofala/Idai/sfincs/event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0/gis/region.geojson"   # replace with your region shapefile
background = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_region_background.geojson")
shapefile_sofala = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_province.shp")

# Load the admin3 district in the case study region to validate exposed people
districts = gpd.read_file("p:/11210471-001-compass/01_Data/sofala_geoms/sofala_districts_study_region.shp")

# population in provided inthousand persons per grid cell
population_raster_path_2020 = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_2020_30.tif")  
population_raster_path_1990 = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_1990_30.tif")  

# flood raster
F_flooding = sfincs_dir_F / "floodmap.tif"
CF_flooding = sfincs_dir_CF / "floodmap.tif"

# Flood model subgrid
sfincs_subgrid = BASE_RUN_PATH / "sfincs" / "subgrid" / "dep_subgrid.tif"


#%% Read flood data and background polygons
# --- Flood grid properties ---
with rasterio.open(sfincs_subgrid) as src:
    flood_grid_crs, flood_grid_transform, flood_grid_shape = src.crs, src.transform, (src.height, src.width)

# --- Setup region ---
region = gpd.read_file(shapefile_fp).to_crs(flood_grid_crs)
region_geom = [json.loads(region.to_json())["features"][0]["geometry"]]

# --- Read flood rasters ---
hmax_F_da = rxr.open_rasterio(F_flooding).squeeze("band", drop=True)  # if single-band
hmax_CF_da = rxr.open_rasterio(CF_flooding).squeeze("band", drop=True)  # if single-band
hmax_F = hmax_F_da.values
hmax_CF = hmax_CF_da.values

# get extent from raster transform
def get_extent(transform, width, height):
    left = transform[2]
    right = left + width * transform[0]
    top = transform[5]
    bottom = top + height * transform[4]
    return [left, right, top, bottom]

flood_extent = get_extent(flood_grid_transform, flood_grid_shape[1], flood_grid_shape[0])

# --- Background layers for plotting ---
# Define a polygon to remove/mask out land layer (incorrect boundary over the ocean)
mask_poly = Polygon([(34.9,-20.3), (36,-20.3), (36,-19.9), (34.9,-19.9)])
bg_filtered = background.copy()
bg_filtered['geometry'] = bg_filtered.geometry.apply(lambda g: g.difference(mask_poly))

# Reproject background and region to flood grid CRS for consistent plotting
background_utm = background.to_crs(flood_grid_crs)
bg_filtered_utm = bg_filtered.to_crs(flood_grid_crs)
region_utm = region.to_crs(flood_grid_crs)
districts_utm = districts.to_crs(flood_grid_crs)

# Remove districts that are not connecting to the region
drop_districts = ["Muanza", "Gororngosa-Sede", "Galinha"]
districts_filtered = districts_utm[~districts_utm['NAME_3'].isin(drop_districts)]

#%%
####################################################################
############### Load World Pop Age & Sex Structures ################
####################################################################
# WorldPop data source: https://hub.worldpop.org/project/categories?id=8

# Function to read, clip, and combine age-sex raster files
def read_and_combine_age_sex_rasters(path_combined_ds, folder_raster, region_gdf=None,
                                     year=None, exclude_genders=None, chunksize=1000):
    """
    Read, clip, and combine WorldPop age-sex raster files into a single 4D xarray Dataset:
        dimensions: genders, ages, y, x
    """

    if os.path.exists(path_combined_ds):
        print(f"📂 Loaded existing combined dataset: {path_combined_ds}")
        return xr.open_dataset(path_combined_ds)

    print(f"🔄 Creating combined dataset from {folder_raster}")

    files = sorted(glob.glob(os.path.join(folder_raster, "moz_*.tif")))
    if not files:
        raise FileNotFoundError(f"No raster files found in {folder_raster}")

    data_by_gender_age = {}

    for f in files:
        match = re.search(r"_([fmt])_(\d+)", f)
        if not match:
            print(f"⚠️ Skipping unrecognized filename: {os.path.basename(f)}")
            continue

        gender, age = match.group(1), int(match.group(2))
        if exclude_genders and gender in exclude_genders:
            print(f"🚫 Skipping excluded gender '{gender}': {os.path.basename(f)}")
            continue

        # Lazy read
        da = rxr.open_rasterio(f, chunks={'x': chunksize, 'y': chunksize}).squeeze()
        da = da.where(da != -99999, 0)

        # --- Clip to region if provided ---
        if region_gdf is not None:
            region_proj = region_gdf.to_crs(da.rio.crs)
            geoms = list(region_proj.geometry)
            da = da.rio.clip(geoms, da.rio.crs, drop=True)

        data_by_gender_age.setdefault(gender, {})[age] = da

    # --- Convert to DataArrays and combine ---
    all_genders = sorted(data_by_gender_age.keys())
    all_ages = sorted({age for gdict in data_by_gender_age.values() for age in gdict.keys()})

    da_list, gender_list, age_list = [], [], []

    for gender in all_genders:
        for age in all_ages:
            if age not in data_by_gender_age[gender]:
                continue
            da_list.append(data_by_gender_age[gender][age])
            gender_list.append(gender)
            age_list.append(age)

    # Build MultiIndex for gender-age
    multi_index = pd.MultiIndex.from_arrays([gender_list, age_list], names=('gender', 'age'))

    # Concatenate along this combined dimension
    ds_all = xr.concat(da_list, dim=pd.Index(range(len(multi_index)), name='z'))
    ds_all = ds_all.assign_coords(z=multi_index).unstack('z')

    ds_all.name = 'population'

    # Set nodata values to 0 instead of -99999 
    ds_all = ds_all.where(ds_all != -99999, 0)

    # Add metadata
    ds_all.attrs.update({
        'title': f'WorldPop Age-Sex Population Data {year or ""}',
        'source': 'WorldPop',
        'resolution': '100m',
        'country': 'Mozambique',
        'created_date': pd.Timestamp.now().isoformat()
    })

    if year is not None:
        ds_all.attrs['year'] = year

    # Optionally save to NetCDF (uncomment if desired)
    # ds_all.to_netcdf(path_combined_ds)
    # print(f"💾 Saved combined dataset to {path_combined_ds}")

    return ds_all


# Function to plot population pyramid from multi-dimensional dataset
def plot_population_pyramid_multidim(ds, combine_ages=None):
    year = ds.attrs['year']
    # Extract data by gender
    female_pop = ds.sel(gender='f').sum(dim=['x', 'y'])
    male_pop = ds.sel(gender='m').sum(dim=['x', 'y'])
    
    ages = female_pop.age.values
    
    if combine_ages:
        # Create new arrays with combined age groups
        new_ages = []
        new_female_pop = []
        new_male_pop = []
        
        used_ages = set()
        
        # Add combined age groups
        for new_age_label, age_list in combine_ages.items():
            available_ages = [age for age in age_list if age in ages]
            if available_ages:
                # Sum population for these ages
                female_combined = female_pop.sel(age=available_ages).sum().values
                male_combined = male_pop.sel(age=available_ages).sum().values
                
                new_ages.append(new_age_label)
                new_female_pop.append(female_combined)
                new_male_pop.append(male_combined)
                
                used_ages.update(available_ages)
        
        # Add remaining individual ages
        for age in ages:
            if age not in used_ages:
                new_ages.append(str(age))
                new_female_pop.append(female_pop.sel(age=age).values)
                new_male_pop.append(male_pop.sel(age=age).values)
        
        ages = new_ages
        female_values = np.array(new_female_pop)
        male_values = np.array(new_male_pop)
    else:
        female_values = female_pop.values
        male_values = male_pop.values
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    y_pos = np.arange(len(ages))
    bar_height = 0.6
    
    # Plot female (negative/left side) and male (positive/right side)
    ax.barh(y_pos, -female_values, bar_height, 
            label='Female', color='#ff6b9d', alpha=0.8)
    ax.barh(y_pos, male_values, bar_height, 
            label='Male', color='#4ecdc4', alpha=0.8)
    
    # Customize
    ax.set_yticks(y_pos)
    ax.set_yticklabels(ages)
    ax.set_xlabel('Population')
    ax.set_ylabel('Age')
    ax.set_title(f'Population Pyramid by Age and Gender - Study region ({year})')
    ax.axvline(x=0, color='black', linewidth=0.8)
    
    # Format x-axis
    ax_ticks = ax.get_xticks()
    ax.set_xticklabels([f'{abs(int(x/1000)):.0f}K' for x in ax_ticks])
    
    ax.legend()
    plt.tight_layout()
    plt.show()


#%%
# Path to all TIFFs
worldpop_folder_2020_1km  = "C:/Code/Paper_2/Data/World_pop/moz_agesex_structures_2020_CN_1km_R2025A_UA_v1"
worldpop_folder_2020_100m = "C:/Code/Paper_2/Data/World_pop/moz_agesex_structures_2020_CN_100m_R2025A_v1"
worldpop_folder_2015_1km  = "C:/Code/Paper_2/Data/World_pop/moz_agesex_structures_2015_CN_1km_R2025A_UA_v1"
worldpop_folder_2030_1km  = "C:/Code/Paper_2/Data/World_pop/moz_agesex_structures_2030_CN_1km_R2025A_UA_v1"

path_moz_agesex_combined_2020_1km  = os.path.join(worldpop_folder_2020_1km, "moz_agesex_structures_2020_CN_1km_R2025A_UA_v1.nc")
path_moz_agesex_combined_2020_100m = os.path.join(worldpop_folder_2020_100m, "moz_agesex_structures_2020_CN_100m_R2025A_v1.nc")
path_moz_agesex_combined_2015_1km  = os.path.join(worldpop_folder_2015_1km, "moz_agesex_structures_2015_CN_1km_R2025A_UA_v1.nc")
path_moz_agesex_combined_2030_1km  = os.path.join(worldpop_folder_2030_1km, "moz_agesex_structures_2030_CN_1km_R2025A_UA_v1.nc")

#%%
# Read and combine WorldPop rasters into a multi-dimensional dataset for 1 km resolution
ds_agesex_MZB_2020_1km = read_and_combine_age_sex_rasters(path_moz_agesex_combined_2020_1km, worldpop_folder_2020_1km, region, year=2020, exclude_genders=['t'])
ds_agesex_MZB_2015_1km = read_and_combine_age_sex_rasters(path_moz_agesex_combined_2015_1km, worldpop_folder_2015_1km, region, year=2015, exclude_genders=['t'])
ds_agesex_MZB_2030_1km = read_and_combine_age_sex_rasters(path_moz_agesex_combined_2030_1km, worldpop_folder_2030_1km, region, year=2030, exclude_genders=['t'])

# Do the same for 100 m resolution dataset
ds_agesex_MZB_2020_100m = read_and_combine_age_sex_rasters(path_moz_agesex_combined_2020_100m, worldpop_folder_2020_100m, region, year=2020, exclude_genders=['t'])

#%%
# Use with combined age groups
plot_population_pyramid_multidim(ds_agesex_MZB_2020_1km, combine_ages={'0-1': [0, 1]})
plot_population_pyramid_multidim(ds_agesex_MZB_2015_1km, combine_ages={'0-1': [0, 1]})
plot_population_pyramid_multidim(ds_agesex_MZB_2030_1km, combine_ages={'0-1': [0, 1]})


#%% ###################################################################################
######################### Overlay population with flood map ###########################
#######################################################################################
# Function to upscale flood raster to population raster grid
def upscale_floodmap_to_pop(pop_raster, flood_array, output_path, region=None):
    if output_path.exists():
        print(f"File {output_path} already exists. Skipping processing.")
        upscaled_flood_array = rxr.open_rasterio(output_path).squeeze("band", drop=True)
        return upscaled_flood_array
    else:   
        pop = pop_raster.squeeze()

        # inputs: hmax_F_da (source), pop_da (target/reference)
        src_da = flood_array
        dst_transform = pop.rio.transform()
        dst_crs = pop.rio.crs
        dst_shape = pop.shape  # (y, x)

        # prepare empty destination array
        dst = np.full(dst_shape, np.nan, dtype=src_da.dtype)

        reproject(
            source=src_da.values,
            destination=dst,
            src_transform=src_da.rio.transform(),
            src_crs=src_da.rio.crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            resampling=Resampling.bilinear,   # better than average!
            src_nodata=getattr(src_da.rio, "nodata", None),
            dst_nodata=np.nan
        )

        upscaled_flood_array = xr.DataArray(
            dst,
            dims=("y", "x"),
            coords={"y": pop_raster.y, "x": pop_raster.x},
            attrs={"crs": dst_crs.to_string()}
        )

        return upscaled_flood_array

# Function to compute population exposed to flooding
def population_exposure_to_flooding(ds_pop, flood_raster_F=None, flood_raster_CF=None, compute_percent=True):
    if flood_raster_F is None and flood_raster_CF is None:
        raise ValueError("At least one of flood_raster_F or flood_raster_CF must be provided")

    # Ensure Dataset
    if isinstance(ds_pop, xr.DataArray):
        ds_pop = ds_pop.to_dataset(name="population")

    new_vars = {}

    # Masks
    if flood_raster_F is not None:
        mask_F = flood_raster_F > 0
        new_vars[f'population_exposed_F'] = xr.where(mask_F, ds_pop['population'], 0)

    if flood_raster_CF is not None:
        mask_CF = flood_raster_CF > 0
        new_vars[f'population_exposed_CF'] = xr.where(mask_CF, ds_pop['population'], 0)

    # Assign exposed population
    ds_pop = ds_pop.assign(**new_vars)

    # --- Compute % exposed and differences ---
    total_pop = ds_pop["population"]
    pop_F = ds_pop[f"population_exposed_F"] if flood_raster_F is not None else None
    pop_CF = ds_pop[f"population_exposed_CF"] if flood_raster_CF is not None else None

    def aggregate_gender_age(data):
        female = data.sel(gender="f").sum(dim=["x","y"])
        male   = data.sel(gender="m").sum(dim=["x","y"])
        return female, male

    female_pop, male_pop = aggregate_gender_age(total_pop)

    if pop_F is not None:
        female_F, male_F = aggregate_gender_age(pop_F)
    if pop_CF is not None:
        female_CF, male_CF = aggregate_gender_age(pop_CF)

    # Optional: compute % exposed
    if compute_percent:
        if pop_F is not None:
            ds_pop[f'population_exposed_F_pct'] = xr.concat(
                [100 * female_F / female_pop, 100 * male_F / male_pop], 
                dim=pd.Index(["f","m"], name="gender")
            )
        if pop_CF is not None:
            ds_pop[f'population_exposed_CF_pct'] = xr.concat(
                [100 * female_CF / female_pop, 100 * male_CF / male_pop], 
                dim=pd.Index(["f","m"], name="gender")
            )
        if pop_F is not None and pop_CF is not None:
            # Absolute difference
            ds_pop[f'population_exposed_diff'] = xr.concat(
                [female_F - female_CF, male_F - male_CF],
                dim=pd.Index(["f","m"], name="gender")
            )
            # % difference
            ds_pop[f'population_exposed_diff_pct'] = xr.concat(
                [100 * (female_F - female_CF) / female_pop, 100 * (male_F - male_CF) / male_pop],
                dim=pd.Index(["f","m"], name="gender")
            )

    print(f"✅ Added exposed population variables (absolute, % and difference)")
    return ds_pop

# Function to plot exposed population pyramids for factual and counterfactual flooding
def plot_exposed_population_pyramids(ds, combine_ages=None, percent=True):
    # Extract data arrays
    total_pop = ds["population"]
    pop_exp_F = ds["population_exposed_F"]
    pop_exp_CF = ds["population_exposed_CF"]

    # Aggregate to gender-age groups
    def aggregate_pop(data, combine_ages=None):
        female = data.sel(gender="f").sum(dim=["x", "y"])
        male = data.sel(gender="m").sum(dim=["x", "y"])
        age = female.age.values

        if combine_ages:
            new_ages, f_vals, m_vals = [], [], []
            used_ages = set()
            for label, age_list in combine_ages.items():
                valid_ages = [a for a in age_list if a in age]
                if valid_ages:
                    new_ages.append(label)
                    f_vals.append(female.sel(age=valid_ages).sum().values)
                    m_vals.append(male.sel(age=valid_ages).sum().values)
                    used_ages.update(valid_ages)
            # Add remaining individual ages
            for a in age:
                if a not in used_ages:
                    new_ages.append(str(a))
                    f_vals.append(female.sel(age=a).values)
                    m_vals.append(male.sel(age=a).values)
            return new_ages, np.array(f_vals), np.array(m_vals)
        else:
            return age, female.values, male.values

    ages, female_pop, male_pop = aggregate_pop(total_pop, combine_ages)
    _, female_F, male_F = aggregate_pop(pop_exp_F, combine_ages)
    _, female_CF, male_CF = aggregate_pop(pop_exp_CF, combine_ages)

    # Compute percentages if requested
    if percent:
        female_F = 100 * female_F / female_pop
        male_F   = 100 * male_F / male_pop
        female_CF = 100 * female_CF / female_pop
        male_CF   = 100 * male_CF / male_pop

    # --- PLOTTING ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharey=True)
    bar_height = 0.6
    y_pos = np.arange(len(ages))

    # Left: factual
    axes[0].barh(y_pos, -female_F, bar_height, color="#ff6b9d", label="Female")
    axes[0].barh(y_pos,  male_F, bar_height, color="#4ecdc4", label="Male")
    axes[0].set_title("Factual exposure" + (" (%)" if percent else ""))

    # Right: counterfactual
    axes[1].barh(y_pos, -female_CF, bar_height, color="#ff6b9d")
    axes[1].barh(y_pos,  male_CF, bar_height, color="#4ecdc4")
    axes[1].set_title("Counterfactual exposure" + (" (%)" if percent else ""))

    # Shared formatting
    for ax in axes:
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(ages)
        ax.set_xlabel("% Exposed" if percent else "Exposed population")
        ax.invert_yaxis()
        xticks = ax.get_xticks()
        ax.set_xticklabels([abs(x) for x in xticks])

    axes[0].legend(loc="lower right")
    fig.suptitle("Exposure Population Pyramids", fontsize=14)
    plt.tight_layout()
    plt.show()

# Function to plot difference in exposed population between factual and counterfactual flooding
def plot_exposure_difference_pyramid(ds, combine_ages=None, percent=False):
    total_pop = ds["population"]
    pop_F = ds["population_exposed_F"]
    pop_CF = ds["population_exposed_CF"]

    # Aggregate to gender-age groups
    def aggregate_pop(data):
        female = data.sel(gender="f").sum(dim=["x", "y"])
        male = data.sel(gender="m").sum(dim=["x", "y"])
        return female, male

    female_pop, male_pop = aggregate_pop(total_pop)
    female_F, male_F = aggregate_pop(pop_F)
    female_CF, male_CF = aggregate_pop(pop_CF)

    ages = female_pop.age.values

    # Optional combine ages
    if combine_ages:
        new_ages, f_vals, m_vals, f_F, m_F, f_CF, m_CF = [], [], [], [], [], [], []
        used_ages = set()
        for label, age_list in combine_ages.items():
            valid_ages = [a for a in age_list if a in ages]
            if valid_ages:
                new_ages.append(label)
                f_vals.append(female_pop.sel(age=valid_ages).sum().values)
                m_vals.append(male_pop.sel(age=valid_ages).sum().values)
                f_F.append(female_F.sel(age=valid_ages).sum().values)
                m_F.append(male_F.sel(age=valid_ages).sum().values)
                f_CF.append(female_CF.sel(age=valid_ages).sum().values)
                m_CF.append(male_CF.sel(age=valid_ages).sum().values)
                used_ages.update(valid_ages)
        # Add remaining ages
        for a in ages:
            if a not in used_ages:
                new_ages.append(str(a))
                f_vals.append(female_pop.sel(age=a).values)
                m_vals.append(male_pop.sel(age=a).values)
                f_F.append(female_F.sel(age=a).values)
                m_F.append(male_F.sel(age=a).values)
                f_CF.append(female_CF.sel(age=a).values)
                m_CF.append(male_CF.sel(age=a).values)
        ages = new_ages
        female_pop, male_pop = np.array(f_vals), np.array(m_vals)
        female_F, male_F = np.array(f_F), np.array(m_F)
        female_CF, male_CF = np.array(f_CF), np.array(m_CF)

    # Compute differences
    diff_female = female_F - female_CF
    diff_male   = male_F - male_CF

    if percent:
        diff_female = 100 * diff_female / female_pop
        diff_male   = 100 * diff_male / male_pop

    # --- PLOTTING ---
    fig, ax = plt.subplots(figsize=(8, 7))
    bar_height = 0.6
    y_pos = np.arange(len(ages))

    ax.barh(y_pos, -diff_female, bar_height, color="#ff6b9d", label="Female")
    ax.barh(y_pos,  diff_male,   bar_height, color="#4ecdc4", label="Male")

    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(ages)
    ax.set_xlabel("% Difference" if percent else "Absolute Difference")
    ax.set_title("Difference in Exposed Population (F - CF)")
    ax.legend(loc="lower right")

    plt.tight_layout()
    plt.show()

#%%
# floodmap to pop 100 m grid
output_F_tif = Path(sfincs_dir_F / "floodmap_100m.tif")
output_CF_tif = Path(sfincs_dir_CF / "floodmap_100m.tif")

# Use the first population raster as reference
first_pop_raster = ds_agesex_MZB_2020_100m.isel(gender=0, age=0)

# Upscale flood rasters to population raster grid (100 m)
hmax_F_100m = upscale_floodmap_to_pop(first_pop_raster, hmax_F_da, output_F_tif, region=region)
hmax_CF_100m = upscale_floodmap_to_pop(first_pop_raster, hmax_CF_da, output_CF_tif, region=region)

#%%
# overlay every 100 m clipped population raster with flood map
ds_agesex_MZB_2020_100m = population_exposure_to_flooding(ds_agesex_MZB_2020_100m, flood_raster_F=hmax_F_100m, flood_raster_CF=hmax_CF_100m)

#%%
# Plot the exposed population pyramids for factual and counterfactual flooding
plot_exposed_population_pyramids(ds_agesex_MZB_2020_100m, combine_ages={'0-1': [0, 1]}, percent=True)


#%%
# Plot the difference in exposed population between factual and counterfactual flooding
plot_exposure_difference_pyramid(ds_agesex_MZB_2020_100m, combine_ages={'0-1': [0, 1]}, percent=False)


# %%
def mask_population_per_district(ds_pop, district_geoms, district_names=None):
    masked_districts = {}
    for i, geom in enumerate(district_geoms):
        masked_ds = mask_population_per_district_single(ds_pop, geom)
        name = district_names[i] if district_names is not None else f"district_{i}"
        masked_districts[name] = masked_ds
    return masked_districts

def mask_population_per_district_single(ds_pop, district_geom):
    transform = ds_pop.rio.transform()
    mask = rasterio.features.rasterize(
        [(district_geom, 1)],
        out_shape=(len(ds_pop.y), len(ds_pop.x)),
        transform=transform,
        fill=0,
        dtype=np.uint8
    )
    mask_da = xr.DataArray(mask, coords={'y': ds_pop.y, 'x': ds_pop.x}, dims=['y', 'x'])

    pop_vars = [v for v in ds_pop.data_vars if 'population' in v]
    masked_vars = {v: ds_pop[v] * mask_da for v in pop_vars}

    ds_masked = ds_pop.copy()
    ds_masked.update(masked_vars)
    return ds_masked


#%%# Mask population data per district
# Remove districts that are not connecting to the region
districts_clean = districts[['NAME_3', 'geometry']][~districts['NAME_3'].isin(["Nhamatanda", "Estaquinha"])]


masked_districts = mask_population_per_district(
    ds_agesex_MZB_2020_100m,
    district_geoms=districts_clean.geometry.values,
    district_names=districts_clean.NAME_3.values
)

# %%
plot_exposed_population_pyramids(masked_districts['Cidade Da Beira'], combine_ages={'0-1': [0, 1]}, percent=True)
# %%
plot_exposed_population_pyramids(masked_districts['Sofala'], combine_ages={'0-1': [0, 1]}, percent=True)

# %%

plot_exposed_population_pyramids(masked_districts['Buzi'], combine_ages={'0-1': [0, 1]}, percent=True)


plot_exposed_population_pyramids(masked_districts['Dondo'], combine_ages={'0-1': [0, 1]}, percent=True)


plot_exposed_population_pyramids(masked_districts['Mafambisse'], combine_ages={'0-1': [0, 1]}, percent=True)

#%%
# plot_exposed_population_pyramids(masked_districts['Nhamatanda'], combine_ages={'0-1': [0, 1]}, percent=True)

plot_exposed_population_pyramids(masked_districts['Tica'], combine_ages={'0-1': [0, 1]}, percent=True)

plot_exposed_population_pyramids(masked_districts['Gororngosa-Sede'], combine_ages={'0-1': [0, 1]}, percent=True)


# %%
print("Plot difference in exposed population pyramids per district")
print("for Beira")
plot_exposure_difference_pyramid(masked_districts['Cidade Da Beira'], combine_ages={'0-1': [0, 1]}, percent=True)
# %%
print("for Sofala")
plot_exposure_difference_pyramid(masked_districts['Sofala'], combine_ages={'0-1': [0, 1]}, percent=True)
# %%
print("for Buzi")
plot_exposure_difference_pyramid(masked_districts['Buzi'], combine_ages={'0-1': [0, 1]}, percent=True)
# %%
print("for Dondo")
plot_exposure_difference_pyramid(masked_districts['Dondo'], combine_ages={'0-1': [0, 1]}, percent=True)
# %%
print("for Mafambisse")
plot_exposure_difference_pyramid(masked_districts['Mafambisse'], combine_ages={'0-1': [0, 1]}, percent=True)
# %%    
print("for Tica")
plot_exposure_difference_pyramid(masked_districts['Tica'], combine_ages={'0-1': [0, 1]}, percent=True)
# %%

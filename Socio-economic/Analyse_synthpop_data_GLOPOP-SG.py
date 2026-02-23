#%% 
# First use read_synthpop_data_GLOPOP-SG.py to read the .dat file and combine it with the .tif file to add coordinates. 
# Then use this script to analyze the combined data using pixi env compass-socio.
import os
from pathlib import Path
import platform
import numpy as np
import pandas as pd
import gzip
import rasterio
import rioxarray as rxr
import geopandas as gpd
import json
from shapely import Point
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
from hydromt import DataCatalog
import rasterio.features as features
from rasterio.mask import mask
from affine import Affine
from shapely.geometry import box
from shapely.geometry import Polygon
from tqdm import tqdm
from rasterio.warp import reproject, Resampling

prefix = "p:/" if platform.system() == "Windows" else "/p/"

#%%
# ===== FILE PATHS =====
# Base directory for the specific event and scenario
BASE_RUN_PATH = Path("p:/11210471-001-compass/03_Runs/sofala/Idai")
BASE_DATA_PATH = Path("p:/11210471-001-compass/01_Data")
sfincs_dir_F  = BASE_RUN_PATH / "sfincs" / "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0"
sfincs_dir_CF = BASE_RUN_PATH / "sfincs" / "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.1_era5_hourly_spw_IBTrACS_CF-5"


# ===== DATA CATALOG =====
if platform.system() == "Windows":
    datacat_path = os.path.abspath("../Workflows/03_data_catalogs/datacatalog_general.yml")
else:
    datacat_path = os.path.abspath("../Workflows/03_data_catalogs/datacatalog_general___linux.yml")
data_catalog = DataCatalog(data_libs = [datacat_path])


# ===== INPUT FILES =====
# Flood model subgrid
sfincs_subgrid = os.path.join(sfincs_dir_F, "subgrid", "dep_subgrid.tif")
with rasterio.open(sfincs_subgrid) as src:
    flood_grid_crs, flood_grid_transform, flood_grid_shape = src.crs, src.transform, (src.height, src.width)

# flood rasters
hmax_F = rxr.open_rasterio(sfincs_dir_F / "plot_output" / "floodmap.tif").squeeze("band", drop=True).values
hmax_CF = rxr.open_rasterio(sfincs_dir_CF / "plot_output" / "floodmap.tif").squeeze("band", drop=True).values

# get extent from raster transform
def get_extent(transform, width, height):
    left = transform[2]
    right = left + width * transform[0]
    top = transform[5]
    bottom = top + height * transform[4]
    return [left, right, top, bottom]

flood_extent = get_extent(flood_grid_transform, flood_grid_shape[1], flood_grid_shape[0])


# Background layers for plotting
mask_poly = Polygon([(34.9,-20.3), (36,-20.3), (36,-19.9), (34.9,-19.9)])
region = gpd.read_file(os.path.join(sfincs_dir_F, "gis/region.geojson")).to_crs(flood_grid_crs)
region_geom = [json.loads(region.to_json())["features"][0]["geometry"]]

background = gpd.read_file(os.path.join(BASE_DATA_PATH, "sofala_geoms", "sofala_region_background.geojson"), driver="GeoJSON")
bg_filtered = background.copy()
bg_filtered['geometry'] = bg_filtered.geometry.apply(lambda g: g.difference(mask_poly))
background_utm = background.to_crs(flood_grid_crs)
bg_filtered_utm = bg_filtered.to_crs(flood_grid_crs)


# Load administrative boundaries
shapefile_sofala = gpd.read_file(os.path.join(BASE_DATA_PATH, "sofala_geoms", "sofala_province.shp"))
districts_adm3 = gpd.read_file(os.path.join(BASE_DATA_PATH, "sofala_geoms", "sofala_districts_study_region.shp"))
districts_adm2 = data_catalog.get_geodataframe("gadm_level2", geom=region, buffer=1000)
districts_adm3_utm = districts_adm3.to_crs(flood_grid_crs)
districts_adm2_utm = districts_adm2.to_crs(flood_grid_crs)

# Remove districts that are not connecting to the region
drop_districts = ["Muanza", "Gororngosa-Sede", "Galinha"]
districts_adm3_filtered = districts_adm3_utm[~districts_adm3_utm['NAME_3'].isin(drop_districts)]


# ===== HOUSEHOLD AND POPULATION DATA =====
# Dominik's data: population in provided in thousand persons per grid cell
population_raster_path_2019 = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_2019_30.tif")  
population_raster_path_1990 = Path("c:/Code/COMPASS_exposure/Data/Outputs/Population/Pop_1990_30.tif")  

# GLOPOP-SG: Read population and their socio-economic characteristics 
df_pop_charac = pd.read_csv(Path('data', 'GLOPOP-SG', 'synthpop_MOZr107_grid_combined.csv'))
glopop_path = Path('data', 'GLOPOP-SG', 'MOZr107_population.tif')
pop_grid = rxr.open_rasterio(glopop_path).squeeze("band", drop=True)  
year_glopop = 2015

# And its corresponding .tif file to add coordinates of the grid cells
grid_nr_filepath = Path('data', 'GLOPOP-SG', 'MOZr107_grid_nr.tif')
with rasterio.open(grid_nr_filepath) as src:
    grid_ids = src.read(1)                   # read first band
    grid_id_transform = src.transform      
    grid_id_crs = src.crs      
    nodata = src.nodata
    profile = src.profile





#%%  ===== FUNCTIONS FOR POPULATION EXPOSURE =====
# Redistribute population and Grid IDs to flood grid and land mask
def reproject_and_redistribute_population_over_land(pop_path, grid_id_path, land_gdf, flood_crs, flood_transform, flood_shape, province_geom=None, region=None, districts_adm3=None, districts_adm2=None, year=None, out_raster_path=None):    
    print(f"▶ Loading {year} population raster...")
    with rasterio.open(pop_path) as src:
        pop = src.read(1, masked=True)
        pop_affine = src.transform
        pop_crs = src.crs

        print(f"Population raster CRS: {pop_crs}")

        # Clip to province
        province_geom = province_geom.to_crs(src.crs)
        province_geom = [json.loads(province_geom.to_json())["features"][0]["geometry"]]

        pop_sofala, transform_sofala = mask(src, province_geom, crop=True, nodata=src.nodata)

        # Dissolve districts into one geometry
        districts_adm3_single = districts_adm3.dissolve().reset_index(drop=True)
        districts_adm3_single = districts_adm3_single.to_crs(src.crs)
        districts_adm3_geom   = [districts_adm3_single.geometry.iloc[0].__geo_interface__]

        districts_adm2_single = districts_adm2.dissolve().reset_index(drop=True)
        districts_adm2_single = districts_adm2_single.to_crs(src.crs)
        districts_adm2_geom   = [districts_adm2_single.geometry.iloc[0].__geo_interface__]

        pop_districts_adm3, pop_affine_districts_adm3 = mask(src, districts_adm3_geom, crop=True, nodata=src.nodata)
        pop_districts_adm2, pop_affine_districts_adm2 = mask(src, districts_adm2_geom, crop=True, nodata=src.nodata)

        if out_raster_path is not None and os.path.exists(out_raster_path) and os.path.exists(out_raster_path.replace(".tif", "_grid_ID.tif")):
            print(f"▶ Loading existing raster from {out_raster_path}")
            with rasterio.open(out_raster_path) as src:
                pop_fine = src.read(1)
            print(f"▶ Loading existing raster from {out_raster_path.replace('.tif', '_grid_ID.tif')}")
            with rasterio.open(out_raster_path.replace(".tif", "_grid_ID.tif"), "r") as src_id:
                grid_id_fine = src_id.read(1)
            return pop_fine, grid_id_fine, pop_sofala, transform_sofala, pop_districts_adm3, pop_affine_districts_adm3, pop_districts_adm2, pop_affine_districts_adm2

        # Clip to region if provided
        if region is not None:
            region_wsg = region.to_crs(src.crs)
            region_geom = [json.loads(region_wsg.to_json())["features"][0]["geometry"]] 
            pop, pop_affine = rasterio.mask.mask(src, region_geom, crop=True, nodata=src.nodata)

        pop = pop.squeeze()

    # Prepare empty high-resolution array
    pop_fine = np.zeros(flood_shape, dtype=np.float32)

    # Grid nr raster to link population counts to characteristics again
    with rasterio.open(grid_id_path) as src_id:
        grid_id = src_id.read(1)
        grid_id_transform = src_id.transform
        grid_id_crs = src_id.crs

    # ---- FINE grid ID (aligned with flood grid) ----
    grid_id_fine = np.full(flood_shape, 0, dtype=np.int32)

    reproject(
        source=grid_id,
        destination=grid_id_fine,
        src_transform=grid_id_transform,
        src_crs=grid_id_crs,
        dst_transform=flood_transform,
        dst_crs=flood_crs,
        resampling=Resampling.nearest
    )

    # Rasterize land mask to flood grid (True for land)
    land_mask = features.rasterize(
        [(geom, 1) for geom in land_gdf.geometry],
        out_shape=flood_shape,
        transform=flood_transform,
        fill=0,
        dtype=np.uint8
    ).astype(bool)

    print("  ✔ Land mask created on flood grid.")

    # Loop through each coarse pixel
    print("▶ Redistributing population to fine grid...")
    for row in tqdm(range(pop.shape[0]), desc="  Processing coarse cells"):
        for col in range(pop.shape[1]):
            pop_value = pop[row, col]
            if np.isnan(pop_value) or pop_value <= 0:
                continue

            # Get coarse pixel bounds (in coarse CRS)
            x_min, y_max = pop_affine * (col, row)
            x_max, y_min = pop_affine * (col + 1, row + 1)
            coarse_bounds = box(x_min, y_min, x_max, y_max)

            # Transform to flood CRS
            coarse_bounds_flood = gpd.GeoSeries([coarse_bounds], crs=pop_crs).to_crs(flood_crs).iloc[0]

            # Rasterize this coarse cell footprint to flood grid
            coarse_mask = features.rasterize(
                [(coarse_bounds_flood, 1)],
                out_shape=flood_shape,
                transform=flood_transform,
                fill=0,
                # all_touched=True,
                dtype=np.uint8
            ).astype(bool)

            # Identify land pixels within this coarse cell
            valid_mask = coarse_mask & land_mask
            n_valid = valid_mask.sum()

            # Distribute population evenly over valid pixels
            if n_valid > 0:
                # pop_fine[valid_mask] += pop_value / n_valid # old resulting in "people behind the comma"

                # --- Integer redistribution (whole people only) ---
                P = int(round(pop_value))
                N = n_valid

                base = P // N
                remainder = P % N

                # convert mask indices to linear index list
                valid_indices = np.where(valid_mask)

                # assign the base number to all pixels
                pop_fine[valid_indices] += base

                # assign remainder randomly
                if remainder > 0:
                    perm = np.random.permutation(len(valid_indices[0]))
                    chosen = perm[:remainder]
                    pop_fine[valid_indices[0][chosen], valid_indices[1][chosen]] += 1

    # else: all water → skip or optionally add to nearest land (not done here)
    total_input_pop = float(np.nansum(pop))
    total_output_pop = float(pop_fine.sum())
    diff = abs(total_output_pop - total_input_pop)
    rel_diff = diff / total_input_pop * 100

    # --- 7️⃣ Validation printout ---
    print("  ✔ Redistribution done.")
    print(f"  🔹 Input population:  {total_input_pop:,.0f}")
    print(f"  🔹 Output population: {total_output_pop:,.0f}")
    print(f"  🔹 Difference:        {diff:,.2f} ({rel_diff:.4f} %)")

    if rel_diff > 0.01:
        print("  ⚠ WARNING: Population not perfectly preserved — check CRS or mask alignment!")

    print(f"  ✔ Redistribution done. Total population preserved: {pop_fine.sum():,.0f}")
    
    # Optional: save the result as a GeoTIFF
    if out_raster_path is not None:
        H_pop, W_pop = pop_fine.shape
        H_id, W_id = grid_id_fine.shape
        a = flood_transform.a
        b = flood_transform.b
        c = flood_transform.c
        d = flood_transform.d
        e = flood_transform.e
        f = flood_transform.f

        if e > 0:
            print("  ⚠ Detected positive y-resolution in transform → fixing for QGIS")

            # 1) flip array vertically
            pop_fine_to_write = np.flipud(pop_fine)
            grid_id_fine_to_write = np.flipud(grid_id_fine)

            # 2a) fix y-scale sign and y-origin for pop grid
            new_e = -abs(e)
            new_f = f + e * (H_pop - 1)
            fixed_transform_pop = Affine(a, b, c, d, new_e, new_f)

            # 2b) fix y-scale sign and y-origin for grid ID grid
            new_e = -abs(e)
            new_f = f + e * (H_id - 1)
            fixed_transform_grid_id = Affine(a, b, c, d, new_e, new_f)

        else:
            pop_fine_to_write = pop_fine
            grid_id_fine_to_write = grid_id_fine
            fixed_transform_pop = flood_transform
            fixed_transform_grid_id = flood_transform

        # ------------------------------------------------------------------
        new_profile_pop = {
            "driver": "GTiff",
            "dtype": rasterio.float32,
            "count": 1,
            "height": H_pop,
            "width": W_pop,
            "crs": flood_crs,
            "transform": fixed_transform_pop,
            "compress": "deflate"
        }

        new_profile_grid_id = {
            "driver": "GTiff",
            "dtype": rasterio.int32,
            "count": 1,
            "height": H_id,
            "width": W_id,
            "crs": flood_crs,
            "transform": fixed_transform_grid_id,
            "compress": "deflate"
        }

        print(f"▶ Saving redistributed population raster to {out_raster_path}")
        with rasterio.open(out_raster_path, "w", **new_profile_pop) as dst:
            dst.write(pop_fine_to_write.astype(rasterio.float32), 1)

        print(f"▶ Saving redistributed grid ID raster to {out_raster_path.replace('.tif', '_grid_ID.tif')}")
        with rasterio.open(out_raster_path.replace(".tif", "_grid_ID.tif"), "w", **new_profile_grid_id) as dst:
            dst.write(grid_id_fine_to_write.astype(rasterio.int32), 1)

    return pop_fine, grid_id_fine, pop_sofala, transform_sofala, pop_districts_adm3, pop_affine_districts_adm3, pop_districts_adm2, pop_affine_districts_adm2

# aggregate exposed population and average flood depth to original grid based on grid ID raster
def aggregate_exposed_pop_to_grid_id(pop_array_fine, flood_array, grid_id_array_fine):
    # Mask population to flooded cells
    pop_exposed_array = np.where(flood_array > 0, pop_array_fine, 0)
    flood_depth_exposed_array = np.where(flood_array > 0, flood_array, 0)

    # Flatten arrays for aggregation
    df = pd.DataFrame({
        "grid_id": grid_id_array_fine.flatten(),
        "pop_total": pop_array_fine.flatten(),
        "pop_exposed": pop_exposed_array.flatten(),
        "flood_depth": flood_depth_exposed_array.flatten()
    })

    df = df[df["grid_id"] != 0]

    agg = df.groupby("grid_id").agg(
        total_population=("pop_total", "sum"),
        exposed_population=("pop_exposed", "sum"),
        avg_flood_depth=("flood_depth", lambda x: x[x > 0].mean() if (x > 0).any() else 0.0)
    ).reset_index()

    agg["exposure_ratio"] = agg["exposed_population"] / agg["total_population"]

    print("\n--- SANITY CHECKS ---")

    # Global population conservation
    print("Total fine population:", pop_array_fine.sum())
    print("Total aggregated population:", agg["total_population"].sum())

    # Global exposed population conservation
    print("Total fine exposed population:", pop_exposed_array.sum())
    print("Total aggregated exposed population:", agg["exposed_population"].sum())
    
    # Population preservation percentages
    print(f"Global population preservation: {agg['total_population'].sum() / pop_array_fine.sum() * 100:.4f} %")
    print(f"Global exposed population preservation: {agg['exposed_population'].sum() / pop_exposed_array.sum() * 100:.4f} %")

    # Exposed never exceeds total
    n_exceed = (agg["exposed_population"] > agg["total_population"]).sum()
    print("Cells where exposed > total:", n_exceed, "(should be 0)")

    # Cells with flood but zero exposed population
    n_flood_no_exposure = (
        (agg["avg_flood_depth"] > 0) &
        (agg["exposed_population"] == 0)
    ).sum()
    print("Cells with flood depth > 0 but exposed = 0:", n_flood_no_exposure, "(can be zero when pop lives in other 25 m cells than the flooding, but check)")

    print("----------------------\n")


    return agg


#%%
# Step 1: Downscale total population and grid IDs to the flood raster (25 m).
export_path = "p:/11210471-001-compass/04_Results/Idai_socioeconomic/preprocessed/population_characteristics/"

pop_arrays = {}
pop_sofala_arrays = {}
pop_sofala_districts_adm3 = {}
pop_affine_sofala_districts_adm3 = {}
pop_sofala_districts_adm2 = {}
pop_affine_sofala_districts_adm2 = {}
# --- Reproject to flood grid and redistribute population rasters over land ---
pop_fine_glopop, grid_id_fine_glopop, pop_sofala_arrays_glopop, transform_sofala_land, pop_districts_adm3_glopop, pop_affine_districts_adm3_glopop, pop_districts_adm2_glopop, pop_affine_districts_adm2_glopop = reproject_and_redistribute_population_over_land(
        pop_path=glopop_path, grid_id_path=grid_nr_filepath, land_gdf=background_utm, flood_crs=flood_grid_crs, flood_transform=flood_grid_transform,
        flood_shape=flood_grid_shape, province_geom=shapefile_sofala, region=region, districts_adm3=districts_adm3_filtered,
        districts_adm2=districts_adm2, year=year_glopop,
        out_raster_path=f"{export_path}population_GLOPOP_SG_MOZr107_regrid.tif") 


#%%
# Step 2 & 3) Compute exposed population at 25 m resolution and aggregate it and avg flood depth to grid ID raster (1 km)
gdf_pop_glopop_exposed_F_coarse  = aggregate_exposed_pop_to_grid_id(pop_fine_glopop, hmax_F, grid_id_fine_glopop)
gdf_pop_glopop_exposed_CF_coarse = aggregate_exposed_pop_to_grid_id(pop_fine_glopop, hmax_CF, grid_id_fine_glopop)



# %%
pop_charac_exposed_F = df_pop_charac.merge(gdf_pop_glopop_exposed_F_coarse[["grid_id", "avg_flood_depth", "exposure_ratio"]], on="grid_id", how="left")
pop_charac_exposed_CF = df_pop_charac.merge(gdf_pop_glopop_exposed_CF_coarse[["grid_id", "avg_flood_depth", "exposure_ratio"]], on="grid_id", how="left")


# aggregate individual-level df to grid-cell mean characteristics
grid_stats = df_pop_charac.groupby('grid_id').agg({
    'WEALTH': 'mean',
    'AGE': 'mean',
    'EDUC': 'mean',
    'GENDER': 'mean',
    'HHSIZE_CAT': 'mean',
    # add more columns if needed
}).reset_index()


#%%
# Get mean social characteristics per disctrict
# get raster dimensions
H, W = grid_id_fine_glopop.shape
transform = grid_id_transform

# create coordinates of cell centers
rows, cols = np.indices((H, W))
xs, ys = rasterio.transform.xy(transform, rows, cols, offset='center')
xs = np.array(xs).flatten()
ys = np.array(ys).flatten()
grid_ids = grid_id_fine_glopop.flatten()

# make a GeoDataFrame
gdf_grid = gpd.GeoDataFrame({
    'grid_id': grid_ids,
    'geometry': [Point(x, y) for x, y in zip(xs, ys)]
}, crs=districts_adm3.crs)

# join the grid stats to the GeoDataFrame
gdf_grid = gdf_grid.merge(grid_stats, on='grid_id', how='left')

# this assigns each grid cell to a district
gdf_grid_districts = gpd.sjoin(
    gdf_grid, 
    districts_adm3[['geometry', 'district_id']],  # keep only ID
    how='inner',
    predicate='intersects'
)

district_stats = gdf_grid_districts.groupby('district_id').agg({
    'WEALTH': 'mean',
    'AGE': 'mean',
    'EDUC': 'mean',
    'GENDER': 'mean',
    'HHSIZE_CAT': 'mean',
    # you can also compute sum if needed, e.g., population
}).reset_index()


# %% #########################################################################
######################### ======== PLOTTING ======== #########################
##############################################################################
# helper functions
def recode_educ(x):
    if 1 <= x <= 2:
        return "E1-2"
    elif 3 <= x <= 5:
        return "E3-5"

def recode_wealth(x):
    if x in [1, 2]:
        return "W1-2"
    elif x == 3:
        return "W3"
    elif x in [4, 5]:
        return "W4-5"

def recode_age(x):
    if 2 <= x <= 7:
        return "A2-7"
    elif x == 1 or x == 8:
        return "A1/8"

def recode_age_coarse(x):
    if 1 <= x <= 8:
        return "A1-8"

def recode_gender(x):
    if 0 <= x <= 1:
        return "G0-1"
        
def recode_household(x):
    if x == 1:
        return "H1"
    elif 2 <= x <= 6:
        return "H2-6"

def recode_household_coarse(x):
    if 1 <= x <= 6:
        return "H1-6"

def make_group_label(row):
    return f"W{int(row['WEALTH'])}_R{int(row['RURAL'])}_A{int(row['AGE'])}_E{int(row['EDUC'])}_G{int(row['GENDER'])}"

def make_group_label_med(row):
    return (
        f"{recode_wealth(int(row['WEALTH']))}_"
        f"R{int(row['RURAL'])}_"
        f"{recode_age(int(row['AGE']))}_"
        f"{recode_educ(int(row['EDUC']))}_"
        f"G{int(row['GENDER'])}_"
        f"{recode_household(int(row['HHSIZE_CAT']))}"
    )

def make_group_label_coarse(row):
    return (
        f"{recode_wealth(int(row['WEALTH']))}_"
        f"R{int(row['RURAL'])}_"
        f"{recode_age_coarse(int(row['AGE']))}_"
        f"{recode_educ(int(row['EDUC']))}_"
        f"{recode_gender(int(row['GENDER']))}_"
        f"{recode_household_coarse(int(row['HHSIZE_CAT']))}"
    )

def categorize_flood_depth(df):
    df['flood_category'] = pd.cut(df['avg_flood_depth'], bins=flood_bins_edges, labels=flood_bins_labels, right=False)
    return df

def aggregate_exposure(df, var):
    return (df.groupby([var, 'flood_category'])['exposure_ratio'].sum().reset_index()
            .rename(columns={'exposure_ratio': 'weighted_exposed'}))


# Available socio-economic characteristics with data
socio_vars = ['WEALTH', 'RURAL', 'AGE', 'EDUC', 'HHSIZE_CAT', 'GENDER']

pop_charac_exposed_F['group_label'] = pop_charac_exposed_F.apply(make_group_label, axis=1)
pop_charac_exposed_CF['group_label'] = pop_charac_exposed_CF.apply(make_group_label, axis=1)
pop_charac_exposed_F['group_label_med'] = pop_charac_exposed_F.apply(make_group_label_med, axis=1)
pop_charac_exposed_CF['group_label_med'] = pop_charac_exposed_CF.apply(make_group_label_med, axis=1)
pop_charac_exposed_F['group_label_coarse'] = pop_charac_exposed_F.apply(make_group_label_coarse, axis=1)
pop_charac_exposed_CF['group_label_coarse'] = pop_charac_exposed_CF.apply(make_group_label_coarse, axis=1)

# --- Define flood depth categories ---
flood_bins_labels = ['Low', 'Medium', 'High']
flood_bins_edges = [0, 0.5, 1.5, 3.5]  # Low <0.5m, Medium 0.5–1.5m, High >1.5m

pop_charac_exposed_F = categorize_flood_depth(pop_charac_exposed_F)
pop_charac_exposed_CF = categorize_flood_depth(pop_charac_exposed_CF)

# Define flood depth bins (0 to 3.5 m with 0.01 m intervals)
flood_bins = np.arange(0, 3.5 + 0.1, 0.02)
pop_charac_exposed_F['flood_bin'] = pd.cut(pop_charac_exposed_F['avg_flood_depth'], bins=flood_bins, right=False)
pop_charac_exposed_CF['flood_bin'] = pd.cut(pop_charac_exposed_CF['avg_flood_depth'], bins=flood_bins, right=False)




#%%
# Plotting flood depth bins separately for each socio-economic variable
for var in socio_vars:
    # Create labels for this variable
    pop_charac_exposed_F[f'{var}_label'] = pop_charac_exposed_F[var].astype(int).apply(lambda x: f"{var[0]}{x}")
    pop_charac_exposed_CF[f'{var}_label'] = pop_charac_exposed_CF[var].astype(int).apply(lambda x: f"{var[0]}{x}")
    
    # Aggregate factual and counterfactual
    agg_F = (pop_charac_exposed_F.groupby([f'{var}_label', 'flood_bin'])['exposure_ratio'].sum().reset_index().rename(columns={'exposure_ratio': 'weighted_pop_F'}))
    agg_CF = (pop_charac_exposed_CF.groupby([f'{var}_label', 'flood_bin'])['exposure_ratio'].sum().reset_index().rename(columns={'exposure_ratio': 'weighted_pop_CF'}))
    # Merge F and CF
    agg = pd.merge(agg_F, agg_CF, on=[f'{var}_label', 'flood_bin'], how='outer')
    agg[['weighted_pop_F', 'weighted_pop_CF']] = agg[['weighted_pop_F', 'weighted_pop_CF']].fillna(0)
    
    # Compute midpoint of flood bins
    agg['flood_mid'] = agg['flood_bin'].apply(lambda x: (x.left + x.right)/2)

    groups = agg[f'{var}_label'].unique()
    n_groups = len(groups)
    
    fig, axes = plt.subplots(n_groups, 1, figsize=(7, n_groups*1), sharex=True)

    tick_bins = np.arange(0.05, 3.5 + 0.5, 0.5)
    
    for i, grp in enumerate(groups):
        ax = axes[i] if n_groups > 1 else axes
        sub = agg[agg[f'{var}_label'] == grp]
        
        ax.plot(sub['flood_mid'], sub['weighted_pop_F'], color='#1f77b4', linewidth=1, label='Factual')
        ax.plot(sub['flood_mid'], sub['weighted_pop_CF'], color='#ff7f0e', linewidth=1, label='CF')
        ax.set_ylabel(grp)
        ax.set_ylim(0, max(agg['weighted_pop_F'].max(), agg['weighted_pop_CF'].max())*1.1)
        ax.set_xlim(0.05, 3.5)
        ax.set_xticks(tick_bins)
        ax.set_xticklabels([f"{x:.1f}" for x in tick_bins], rotation=45)
        
        if i == 0:
            ax.legend()
        if i == n_groups - 1:
            ax.set_xlabel('Flood depth [m]')
        else:
            ax.set_xlabel('')
    
    plt.suptitle(f"Population exposure per {var}")
    plt.tight_layout()
    plt.show()




# %%
# --- Choose socio-economic variable separately ---
var = 'AGE_label'

# --- Aggregate absolute exposed population per category and flood level ---
agg_F = aggregate_exposure(pop_charac_exposed_F, var)
agg_CF = aggregate_exposure(pop_charac_exposed_CF, var)

# --- Calculate % change ---
agg = pd.merge(agg_F, agg_CF, on=[var, 'flood_category'], how='outer', suffixes=('_F', '_CF'))
agg['pct_change'] = ((agg['weighted_exposed_F'] - agg['weighted_exposed_CF']) / agg['weighted_exposed_F'] * 100).round(1)

# --- Pivot to MultiIndex columns: top = category, sub = F / CF / % ---
table = agg.pivot(index='flood_category', columns=var, values=['weighted_exposed_F', 'weighted_exposed_CF', 'pct_change'])

# Reorder levels: top = category, second = F/CF/%
table = table.reorder_levels([1,0], axis=1)

# Rename metrics to F / CF / %
flat_cols = []
for cat, metric in table.columns:
    if metric == 'weighted_exposed_F':
        flat_cols.append((cat, 'F'))
    elif metric == 'weighted_exposed_CF':
        flat_cols.append((cat, 'CF'))
    else:
        flat_cols.append((cat, '%'))

table.columns = pd.MultiIndex.from_tuples(flat_cols)

# Sort categories alphabetically (W1, W2, ...)
table = table.sort_index(axis=1, level=0)

# Optional: reset index for printing
# table = table.reset_index()

table_swapped = table.T

table_swapped

# # Flatten MultiIndex for CSV export
# table_csv = table.copy()
# table_csv.columns = ['_'.join(col).strip() if isinstance(col, tuple) else col for col in table_csv.columns.values]

# # Export to CSV
# table_csv.to_csv("results/exposed_population_by_category.csv", index=False)



# %%
# All individual variables separately 
variables = {
    'AGE': 'A',
    'WEALTH': 'W',
    'EDUC': 'E',
    'RURAL': 'R',
    'HHSIZE_CAT': 'H',
    'GENDER': 'G'
}

for col, prefix in variables.items():
    pop_charac_exposed_F[f"{col}_label"] = (
        pop_charac_exposed_F[col].astype(int).apply(lambda x: f"{prefix}{x}"))
    pop_charac_exposed_CF[f"{col}_label"] = (
        pop_charac_exposed_CF[col].astype(int).apply(lambda x: f"{prefix}{x}"))

all_tables = []

for col, prefix in variables.items():
    var = f"{col}_label"

    agg_F = aggregate_exposure(pop_charac_exposed_F, var)
    agg_CF = aggregate_exposure(pop_charac_exposed_CF, var)
    agg = pd.merge(agg_F, agg_CF, on=[var, 'flood_category'], how='outer', suffixes=('_F', '_CF'))

    # % change
    agg['pct_change'] = ((agg['weighted_exposed_CF'] - agg['weighted_exposed_F'])
                        / agg['weighted_exposed_F'] * 100)

    # Pivot
    table = agg.pivot(index='flood_category', columns=var,
                      values=['weighted_exposed_F', 'weighted_exposed_CF', 'pct_change'])

    # Ensure correct flood category order
    table = table.reindex(flood_bins_labels)
    table = table.reorder_levels([1, 0], axis=1)

    # Rename metrics
    renamed_cols = []
    for cat, metric in table.columns:
        if metric == 'weighted_exposed_F':
            renamed_cols.append((cat, 'F'))
        elif metric == 'weighted_exposed_CF':
            renamed_cols.append((cat, 'CF'))
        else:
            renamed_cols.append((cat, '%'))

    table.columns = pd.MultiIndex.from_tuples(renamed_cols)
    table = table.sort_index(axis=1, level=0)

    # Transpose → categories become rows
    table_t = table.T

    # Force row order inside each category
    metric_order = ['%', 'CF', 'F']
    new_index = []
    for cat in sorted(set([i[0] for i in table_t.index])):
        for m in metric_order:
            if (cat, m) in table_t.index:
                new_index.append((cat, m))

    table_t = table_t.loc[new_index]

    # Make clean MultiIndex row structure
    table_t.index = pd.MultiIndex.from_tuples(table_t.index, names=[col, 'metric'])
    all_tables.append(table_t)

# Combine all variables vertically
final_table = pd.concat(all_tables)
final_table = final_table.reset_index()
final_table.columns.name = None
final_table = final_table.rename(columns={
    final_table.columns[0]: "Category",
    final_table.columns[1]: "Metric"
})

final_table


# %%
# --- Aggregate by group_label and flood category ---
agg_exposed_F = (pop_charac_exposed_F
    .groupby(['group_label_med', 'flood_category'])['exposure_ratio']
    .sum().reset_index().rename(columns={'exposure_ratio': 'F'}))

agg_exposed_CF = (pop_charac_exposed_CF
    .groupby(['group_label_med', 'flood_category'])['exposure_ratio']
    .sum().reset_index().rename(columns={'exposure_ratio': 'CF'}))

agg_pop = (pop_charac_exposed_F
    .groupby(['group_label_med'])
    .size().reset_index(name='total_population'))

# Merge F and CF exposed and total pop counts
agg = pd.merge(agg_exposed_F, agg_exposed_CF, on=['group_label_med', 'flood_category'], how='outer')
agg = agg.merge(agg_pop, on=['group_label_med'], how='left')

# --- Calculate % change only - per flood depth ---
agg['pct_change'] = np.where(agg['F'] > 0, (agg['F'] - agg['CF']) / agg['F'] * 100, 
                             np.nan).round(2)
agg['F_rate'] = agg['F'] / agg['total_population'] * 100
agg['CF_rate'] = agg['CF'] / agg['total_population'] * 100
agg['rate_diff'] = agg['F_rate'] - agg['CF_rate']

# --- Pivot: rows = group_label_med, columns = flood depth ---
pct_table = agg.pivot(
    index='group_label_med',
    columns='flood_category',
    values=['F', 'CF', 'pct_change', 'F_rate', 'CF_rate', 'rate_diff', 'total_population'])

table = pct_table.reorder_levels([1, 0], axis=1)

# Rename metrics for clarity
new_cols = []
for flood_cat, metric in table.columns:
    if metric == 'F':
        new_cols.append((flood_cat, 'F'))
    elif metric == 'CF':
        new_cols.append((flood_cat, 'CF'))
    elif metric == 'pct_change':
        new_cols.append((flood_cat, '%'))
    elif metric == 'total_population':
        new_cols.append((flood_cat, 'Total_pop'))
    elif metric == 'F_rate':
        new_cols.append((flood_cat, 'F_rate'))
    elif metric == 'CF_rate':   
        new_cols.append((flood_cat, 'CF_rate'))
    elif metric == 'rate_diff':
        new_cols.append((flood_cat, 'rate_diff'))

table.columns = pd.MultiIndex.from_tuples(new_cols)

# Sort flood categories
table = table.sort_index(axis=1, level=0)

# Optional: round numeric values safely
for col in table.columns.levels[0]:  # top-level: flood categories
    # Make a copy first to avoid pandas recursion issues
    table[(col, 'F')] = table[(col, 'F')].copy().round(0)       # population → 0 decimals
    table[(col, 'CF')] = table[(col, 'CF')].copy().round(0)
    table[(col, '%')] = pd.to_numeric(table[(col, '%')].copy(), errors='coerce').round(2)  # % → 2 decimals
    table[(col, 'F_rate')] = pd.to_numeric(table[(col, 'F_rate')].copy(), errors='coerce').round(4)  # rates → 4 decimals
    table[(col, 'CF_rate')] = pd.to_numeric(table[(col, 'CF_rate')].copy(), errors='coerce').round(4)
    table[(col, 'rate_diff')] = pd.to_numeric(table[(col, 'rate_diff')].copy(), errors='coerce').round(4)

export_path = "results/exposed_population_by_group_and_flood_category.csv"
table.to_csv(export_path)

table


#%%
agg_exposed_F_coarse = (pop_charac_exposed_F
    .groupby(['group_label_med', 'flood_category'])['exposure_ratio']
    .sum().reset_index().rename(columns={'exposure_ratio': 'F'}))

agg_exposed_CF_coarse = (pop_charac_exposed_CF
    .groupby(['group_label_med', 'flood_category'])['exposure_ratio']
    .sum().reset_index().rename(columns={'exposure_ratio': 'CF'}))

agg_pop_coarse = (pop_charac_exposed_F
    .groupby(['group_label_med'])
    .size().reset_index(name='total_population'))

# Merge F and CF exposed and total pop counts
agg_coarse = pd.merge(agg_exposed_F_coarse, agg_exposed_CF_coarse, on=['group_label_med', 'flood_category'], how='outer')
agg_coarse = agg_coarse.merge(agg_pop_coarse, on=['group_label_med'], how='left')

# Sum over all flood depths for each group_label_med
agg_total = agg_coarse.groupby('group_label_med')[['F', 'CF']].sum().reset_index()

# Calculate % difference across all depths
agg_total['pct_change_total'] = np.where(
    agg_total['F'] > 0,
    (agg_total['F'] - agg_total['CF']) / agg_total['F'] * 100,
    np.nan
).round(2)

agg_coarse['F_rate'] = agg_coarse['F'] / agg_coarse['total_population'] * 100
agg_coarse['CF_rate'] = agg_coarse['CF'] / agg_coarse['total_population'] * 100
agg_coarse['rate_diff'] = agg_coarse['F_rate'] - agg_coarse['CF_rate']

# Optional: also compute total rates
agg_total = agg_total.merge(agg_pop_coarse, on='group_label_med', how='left')
agg_total['F_rate_total'] = agg_total['F'] / agg_total['total_population'] * 100
agg_total['CF_rate_total'] = agg_total['CF'] / agg_total['total_population'] * 100
agg_total['rate_diff_total'] = agg_total['F_rate_total'] - agg_total['CF_rate_total']

MIN_POP = 100  # choose something meaningful
filtered_total = agg_total[agg_total['total_population'] >= MIN_POP]

top10_total = filtered_total.sort_values(['pct_change_total'], ascending=False).head(10)
bottom10_total = filtered_total.sort_values(['pct_change_total'], ascending=True).head(10)


# %%# Sort by Low flood category (descending = highest first)
MIN_POP = 100  # choose something meaningful
filtered = table[table[('High', 'Total_pop')] >= MIN_POP]

top10_low = filtered.sort_values([('Low', '%')], ascending=False).head(10)
bottom10_low = filtered.sort_values([('Low', '%')], ascending=True).head(10)

# Sort by Medium flood category
top10_medium = filtered.sort_values([('Medium', '%')], ascending=False).head(10)
bottom10_medium = filtered.sort_values([('Medium', '%')], ascending=True).head(10)

# Sort by High flood category
top10_high = filtered.sort_values([('High', '%')], ascending=False).head(10)
bottom10_high = filtered.sort_values([('High', '%')], ascending=True).head(10)

print("Top 10 groups with highest % change in exposed population for Low flood category:")
print(top10_low)
print("\nBottom 10 groups with lowest % change in exposed population for Low flood category:")
print(bottom10_low)

print("\nTop 10 groups with highest % change in exposed population for Medium flood category:")
print(top10_medium)
print("\nBottom 10 groups with lowest % change in exposed population for Medium flood category:")
print(bottom10_medium)

print("\nTop 10 groups with highest % change in exposed population for High flood category:")
print(top10_high)
print("\nBottom 10 groups with lowest % change in exposed population for High flood category:")
print(bottom10_high)

# %%

# Suppose your table has these columns:
# 'group_label_coarse', 'Low', 'Medium', 'High' (absolute exposed population)
# Example: table_coarse = df_grouped[['group_label_coarse', 'Low', 'Medium', 'High']]

# Data for plotting
labels = ['Low', 'Medium', 'High']
bar_width = 0.25
# x_pos = np.arange(len(table))  # one tick per coarse group

fig, ax = plt.subplots(figsize=(12,6))

# For each flood depth, plot bars side by side
for i, flood in enumerate(labels):
    ax.bar(x_pos + i*bar_width, table[flood], width=bar_width, label=flood)

# X-axis labels = coarse groups
ax.set_xticks(x_pos + bar_width)  # center the ticks
ax.set_xticklabels(table['group_label_coarse'], rotation=45, ha='right', fontsize=9)

ax.set_ylabel("Exposed population")
ax.set_xlabel("Flood depth")
ax.set_title("Exposed population per flood depth and coarse group")
ax.legend(title="Flood depth")
ax.grid(True, linestyle="--", alpha=0.5)

plt.tight_layout()
plt.show()

# %%
import matplotlib.pyplot as plt

# Make sure your table has '%', for example:
# table.columns = MultiIndex with (flood_depth, metric)
# row index = coarse group label

flood_bins_labels = ['Low', 'Medium', 'High']

fig, ax = plt.subplots(figsize=(10,6))

# Loop over coarse groups (rows)
for group_label in filtered.index:
    # Extract % values for this group across flood depths
    y_values = [filtered.loc[group_label, (flood, '%')] for flood in flood_bins_labels]
    ax.plot(flood_bins_labels, y_values, marker='o', label=group_label)

ax.set_xlabel("Flood depth")
ax.set_ylabel("Percentage change in exposed population")
ax.set_title("Percentage change per coarse social group and flood depth")
ax.axhline(0, color='black', linestyle='--', alpha=0.7)
ax.legend(title="Coarse groups", bbox_to_anchor=(1.05,1), loc='upper left')
ax.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.show()

# %%
flood_bins_labels = ['Low', 'Medium', 'High']
n_groups = len(filtered.index)
groups_per_plot = 100

for i in range(0, n_groups, groups_per_plot):
    fig, ax = plt.subplots(figsize=(10,6))
    
    subgroups = filtered.index[i:i+groups_per_plot]
    
    for group_label in subgroups:
        y_values = [filtered.loc[group_label, (flood, '%')] for flood in flood_bins_labels]
        ax.plot(flood_bins_labels, y_values, marker='o', label=group_label)
    
    ax.set_xlabel("Flood depth")
    ax.set_ylabel("Percentage change in exposed population")
    ax.set_title(f"Percentage change per coarse group (categories {i+1}-{i+len(subgroups)})")
    ax.axhline(0, color='black', linestyle='--', alpha=0.7)
    ax.legend(title="Coarse groups", bbox_to_anchor=(1.05,1), loc='upper left')
    ax.grid(True, linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.show()


# %%
flood_bins_labels = ['Low', 'Medium', 'High']

selected_groups = [
    "W4-5_R0_A2-7_E3-5_G1_H2-6",  
    "W1-2_R1_A1/8_E1-2_G0_H1",  
    "W4-5_R1_A1/8_E1-2_G0_H1",  
    "W1-2_R0_A1/8_E1-2_G0_H1",  
    "W1-2_R1_A1/8_E1-2_G0_H2-6",  
    "W1-2_R1_A1/8_E1-2_G1_H1",  
    # "W1-2_R1_A1/8_E4-5_G0_H1",    
    "W1-2_R1_A2-7_E1-2_G0_H1",     
]

x = np.arange(len(flood_bins_labels))  # positions for flood depths
width = 0.1  # width of each bar

fig, ax = plt.subplots(figsize=(10,6))

for i, group in enumerate(selected_groups):
    y_values = [filtered.loc[group, (flood, '%')] for flood in flood_bins_labels]
    total_attr_values = filtered_total[filtered_total['group_label_med'] == group]['pct_change_total'].values
    total_exposed_people_F = filtered_total[filtered_total['group_label_med'] == group]['F'].values
    ax.bar(x + i*width, y_values, width, label=f'{group} ({total_attr_values[0]:.1f}% - {total_exposed_people_F[0]:.0f} people)')

ax.set_xticks(x + width*1.5)
ax.set_xticklabels(flood_bins_labels)

ax.set_xlabel("Flood depth")
ax.set_ylabel("Percentage change in exposed population")
ax.set_title("Comparison of % change per flood depth\nSelected vulnerability groups")

ax.axhline(0, color='black', linestyle='--', alpha=0.7)
ax.legend(title="Coarse groups", bbox_to_anchor=(1.05,1), loc='upper left')
ax.grid(axis='y', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.show()


# change in fraction of population that is exposed due to climate change
fig, ax = plt.subplots(figsize=(10,6))

for i, group in enumerate(selected_groups):
    y_values = [filtered.loc[group, (flood, 'rate_diff')] for flood in flood_bins_labels]
    total_attr_values = filtered_total[filtered_total['group_label_med'] == group]['rate_diff_total'].values
    total_exposed_people_F = filtered_total[filtered_total['group_label_med'] == group]['F'].values
    total_exposed_people_CF = filtered_total[filtered_total['group_label_med'] == group]['CF'].values
    ax.bar(x + i*width, y_values, width, label=f'{group} ({total_attr_values[0]:.1f}% - {total_exposed_people_CF[0]:.0f}/{total_exposed_people_F[0]:.0f} CF/F people)')

ax.set_xticks(x + width*1.5)
ax.set_xticklabels(flood_bins_labels)

ax.set_xlabel("Flood depth")
ax.set_ylabel("Percentage change in exposed population")
ax.set_title("Comparison of % change per flood depth\nSelected vulnerability groups")

ax.axhline(0, color='black', linestyle='--', alpha=0.7)
ax.legend(title="Coarse groups", bbox_to_anchor=(1.05,1), loc='upper left')
ax.grid(axis='y', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.show()


# %%
age_groups = {
    "Young": "A1-2",
    "Medium": "A3-6",
    "Old": "A8"
}

fig, axes = plt.subplots(1, 3, figsize=(15,6), sharey=True)

for ax, (age_label, age_code) in zip(axes, age_groups.items()):
    
    gender_diff = []
    
    for flood in flood_bins_labels:
        
        # Filter men in this age group
        men = filtered[
            filtered.index.str.contains(age_code) &
            filtered.index.str.contains("G1")
        ][(flood, '%')].mean()
        
        # Filter women in this age group
        women = filtered[
            filtered.index.str.contains(age_code) &
            filtered.index.str.contains("G0")
        ][(flood, '%')].mean()
        
        gender_diff.append(men - women)
    
    ax.bar(flood_bins_labels, gender_diff)
    ax.set_title(age_label)
    ax.axhline(0, color='black', linestyle='--')
    ax.set_xlabel("Flood depth")

axes[0].set_ylabel("Men − Women (% change)")
plt.suptitle("Gender difference in % change per flood depth and age group")
plt.tight_layout()
plt.show()


# %%
wealth_groups = {
    "Poor": "W1-2",
    "Rich": "W4-5"
}


fig, axes = plt.subplots(1, 3, figsize=(16,6), sharey=True)

width = 0.35
x = np.arange(len(flood_bins_labels))

for ax, (age_label, age_code) in zip(axes, age_groups.items()):
    
    for i, (wealth_label, wealth_code) in enumerate(wealth_groups.items()):
        
        gender_diff = []
        
        for flood in flood_bins_labels:
            
            # MEN
            men = filtered[
                filtered.index.str.contains(age_code) &
                filtered.index.str.contains(wealth_code) &
                filtered.index.str.contains("G1")
            ][(flood, '%')].mean()
            
            # WOMEN
            women = filtered[
                filtered.index.str.contains(age_code) &
                filtered.index.str.contains(wealth_code) &
                filtered.index.str.contains("G0")
            ][(flood, '%')].mean()
            
            gender_diff.append(men - women)
        
        ax.bar(x + i*width, gender_diff, width, label=wealth_label)
    
    ax.set_title(age_label)
    ax.set_xticks(x + width/2)
    ax.set_xticklabels(flood_bins_labels)
    ax.axhline(0, color='black', linestyle='--')

axes[0].set_ylabel("Men − Women (% change)")
axes[0].legend(title="Wealth group")

plt.suptitle("Gender difference in % change by Age and Wealth")
plt.tight_layout()
plt.show()

# %%
hh_groups = {
    "Single": "H1",
    "Family": "H2-6"
}

fig, axes = plt.subplots(1, 3, figsize=(18,6), sharey=True)

width = 0.18
x = np.arange(len(flood_bins_labels))

for ax, (age_label, age_code) in zip(axes, age_groups.items()):
    
    bar_index = 0
    
    for wealth_label, wealth_code in wealth_groups.items():
        for hh_label, hh_code in hh_groups.items():
            
            gender_diff = []
            
            for flood in flood_bins_labels:
                
                men = filtered[
                    filtered.index.str.contains(age_code) &
                    filtered.index.str.contains(wealth_code) &
                    filtered.index.str.contains(hh_code) &
                    filtered.index.str.contains("G1")
                ][(flood, '%')].mean()
                
                women = filtered[
                    filtered.index.str.contains(age_code) &
                    filtered.index.str.contains(wealth_code) &
                    filtered.index.str.contains(hh_code) &
                    filtered.index.str.contains("G0")
                ][(flood, '%')].mean()
                
                gender_diff.append(men - women)
            
            ax.bar(
                x + bar_index*width,
                gender_diff,
                width,
                label=f"{wealth_label} – {hh_label}"
            )
            
            bar_index += 1
    
    ax.set_title(age_label)
    ax.set_xticks(x + width*1.5)
    ax.set_xticklabels(flood_bins_labels)
    ax.axhline(0, color='black', linestyle='--')

axes[0].set_ylabel("Men − Women (% change)")
axes[0].legend(title="Wealth – HH size", bbox_to_anchor=(1.05,1))

plt.suptitle("Gender difference in % change by Age, Wealth, and Household Size")
plt.tight_layout()
plt.show()

# %%
education_groups = {
    "Low education": "E1-3",
    "High education": "E4-5"
}

for edu_label, edu_code in education_groups.items():
    
    fig, axes = plt.subplots(1, 3, figsize=(18,6), sharey=True)
    
    width = 0.18
    x = np.arange(len(flood_bins_labels))
    
    for ax, (age_label, age_code) in zip(axes, age_groups.items()):
        
        bar_index = 0
        
        for wealth_label, wealth_code in wealth_groups.items():
            for hh_label, hh_code in hh_groups.items():
                
                gender_diff = []
                
                for flood in flood_bins_labels:
                    
                    men = filtered[
                        filtered.index.str.contains(age_code) &
                        filtered.index.str.contains(wealth_code) &
                        filtered.index.str.contains(hh_code) &
                        filtered.index.str.contains(edu_code) &
                        filtered.index.str.contains("G1")
                    ][(flood, '%')].mean()
                    
                    women = filtered[
                        filtered.index.str.contains(age_code) &
                        filtered.index.str.contains(wealth_code) &
                        filtered.index.str.contains(hh_code) &
                        filtered.index.str.contains(edu_code) &
                        filtered.index.str.contains("G0")
                    ][(flood, '%')].mean()
                    
                    gender_diff.append(men - women)
                
                ax.bar(
                    x + bar_index*width,
                    gender_diff,
                    width,
                    label=f"{wealth_label} – {hh_label}"
                )
                
                bar_index += 1
        
        ax.set_title(age_label)
        ax.set_xticks(x + width*1.5)
        ax.set_xticklabels(flood_bins_labels)
        ax.axhline(0, color='black', linestyle='--')
        ax.grid(axis='y', linestyle='--', alpha=0.4)
        ax.set_axisbelow(True)   # ensures grid is behind bars
    
    axes[0].set_ylabel("Men − Women (% change)")
    axes[0].legend(title="Wealth – HH size", bbox_to_anchor=(1.05,1))
    
    plt.suptitle(f"Gender difference by Age, Wealth & HH size\n({edu_label})")
    plt.tight_layout()
    plt.show()

# %%


rural_groups = {
    "Rural": "R1",
    "Urban": "R0"
}


for edu_label, edu_code in education_groups.items():
    for rural_label, rural_code in rural_groups.items():
        
        fig, axes = plt.subplots(1, 3, figsize=(18,6), sharey=True)
        width = 0.18
        x = np.arange(len(flood_bins_labels))
        
        for ax, (age_label, age_code) in zip(axes, age_groups.items()):
            
            bar_index = 0
            
            for wealth_label, wealth_code in wealth_groups.items():
                for hh_label, hh_code in hh_groups.items():
                    
                    gender_diff = []
                    
                    for flood in flood_bins_labels:
                        
                        men = filtered[
                            filtered.index.str.contains(age_code) &
                            filtered.index.str.contains(wealth_code) &
                            filtered.index.str.contains(hh_code) &
                            filtered.index.str.contains(edu_code) &
                            filtered.index.str.contains(rural_code) &
                            filtered.index.str.contains("G1")
                        ][(flood, '%')].mean()
                        
                        women = filtered[
                            filtered.index.str.contains(age_code) &
                            filtered.index.str.contains(wealth_code) &
                            filtered.index.str.contains(hh_code) &
                            filtered.index.str.contains(edu_code) &
                            filtered.index.str.contains(rural_code) &
                            filtered.index.str.contains("G0")
                        ][(flood, '%')].mean()
                        
                        gender_diff.append(men - women)
                    
                    ax.bar(
                        x + bar_index*width,
                        gender_diff,
                        width,
                        label=f"{wealth_label} – {hh_label}"
                    )
                    
                    bar_index += 1
            
            ax.set_title(age_label)
            ax.set_xticks(x + width*1.5)
            ax.set_xticklabels(flood_bins_labels)
            ax.axhline(0, color='black', linestyle='--')
            ax.grid(axis='y', linestyle='--', alpha=0.4)
            ax.set_axisbelow(True)
        
        axes[0].set_ylabel("Men − Women (% change)")
        axes[0].legend(title="Wealth – HH size", bbox_to_anchor=(1.05,1))
        
        plt.suptitle(
            f"{edu_label} – {rural_label}\nGender difference by Age, Wealth and HH size"
        )
        
        plt.tight_layout()
        plt.show()

# %%


flood_bins_labels = ['Low', 'Medium', 'High']

age_groups = {
    "Young": "A1-2",
    "Medium": "A3-6",
    "Old": "A8"
}

wealth_groups = {
    "Poor": "W1-2",
    "Rich": "W4-5"
}

hh_groups = {
    "Single": "H1",
    "Family": "H2-6"
}

education_groups = {
    "Low education": "E1-3",
    "High education": "E4-5"
}

rural_groups = {
    "Rural": "R1",
    "Urban": "R0"
}

gender_groups = {
    "Men": "G1",
    "Women": "G0"
}

for edu_label, edu_code in education_groups.items():
    for rural_label, rural_code in rural_groups.items():
        
        fig, axes = plt.subplots(1, 3, figsize=(20,6), sharey=True)
        
        x = np.arange(len(flood_bins_labels))
        width = 0.1  # smaller width since many bars
        
        for ax, (age_label, age_code) in zip(axes, age_groups.items()):
            
            bar_index = 0
            
            for wealth_label, wealth_code in wealth_groups.items():
                for hh_label, hh_code in hh_groups.items():
                    for gender_label, gender_code in gender_groups.items():
                        
                        values = []
                        
                        for flood in flood_bins_labels:
                            
                            val = filtered[
                                filtered.index.str.contains(age_code) &
                                filtered.index.str.contains(wealth_code) &
                                filtered.index.str.contains(hh_code) &
                                filtered.index.str.contains(edu_code) &
                                filtered.index.str.contains(rural_code) &
                                filtered.index.str.contains(gender_code)
                            ][(flood, '%')].mean()
                            
                            values.append(val)
                        
                        ax.bar(
                            x + bar_index*width,
                            values,
                            width,
                            label=f"{wealth_label}-{hh_label}-{gender_label}"
                        )
                        
                        bar_index += 1
            
            ax.set_title(age_label)
            ax.set_xticks(x + width*4)
            ax.set_xticklabels(flood_bins_labels)
            ax.axhline(0, color='black', linestyle='--')
            ax.grid(axis='y', linestyle='--', alpha=0.4)
            ax.set_axisbelow(True)
        
        axes[0].set_ylabel("% change in exposed population")
        
        axes[0].legend(
            title="Wealth – HH – Gender",
            bbox_to_anchor=(1.05,1),
            fontsize=8
        )
        
        plt.suptitle(
            f"{edu_label} – {rural_label}\n% Change by Age, Wealth, HH size and Gender"
        )
        
        plt.tight_layout()
        plt.show()

# %%
# Average change you are exposed to flooding due to climate change
exposed_frac_total_pop_F = pop_charac_exposed_F['exposure_ratio'].sum() / len(pop_charac_exposed_F) * 100
exposed_frac_total_pop_CF = pop_charac_exposed_CF['exposure_ratio'].sum() / len(pop_charac_exposed_CF) * 100
change_exposed_frac = exposed_frac_total_pop_F - exposed_frac_total_pop_CF

print(f"Average % of population exposed to flooding in Factual: {exposed_frac_total_pop_F:.2f}%")
print(f"Average % of population exposed to flooding in Climate Future: {exposed_frac_total_pop_CF:.2f}%")
print(f"Change in exposure: {change_exposed_frac:.2f}%")

# % of factual exposed population that can be attributed to climate change
change_attributed_to_climate = (pop_charac_exposed_F['exposure_ratio'].sum() - pop_charac_exposed_CF['exposure_ratio'].sum()) / pop_charac_exposed_F['exposure_ratio'].sum() * 100
print(f"Change in exposure attributed to climate change: {change_attributed_to_climate:.2f}% (based on Dominik's data is 7%)")
# %%

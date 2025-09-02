#%%
import rasterio
import geopandas as gpd
import numpy as np
from shapely.geometry import Point, box
from pyproj import Transformer
import pandas as pd
from pathlib import Path
import rasterio.windows
import xarray as xr
import warnings
warnings.filterwarnings('ignore')
import platform
import os

prefix = "p:/" if platform.system() == "Windows" else "/p/"

# ===== CONFIGURATION =====
EVENT_NAME = "Idai"
# BASE_RUN_PATH = Path("/p/11210471-001-compass/03_Runs/test")
BASE_RUN_PATH = Path(os.path.join(prefix, "11210471-001-compass","03_Runs","sofala"))
SCENARIO_PATH = "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0"

# ===== FILE PATHS =====
# Base directory for the specific event and scenario
base_dir = BASE_RUN_PATH / EVENT_NAME / "fiat" / SCENARIO_PATH

# Input files
# Use the spatial.fgb from output directory (pre-filtered buildings with damage)
buildings_fp_path = base_dir / "output" / "spatial.fgb"
exposure_csv_path = base_dir / "exposure" / "exposure.csv"

# Population raster path
population_raster_path = os.path.join(prefix,"11210471-001-compass","01_Data","population_data","moz_ppp_2020_UNadj_constrained.tif")  

# Flood depth NetCDF path
flood_depth_path = base_dir / "hazard" / "floodmap.nc"

# Output files
output_dir = base_dir / "output"
output_dir.mkdir(exist_ok=True)

output_buildings_path = output_dir / "spatial_with_pop_and_flood.fgb"
output_exposure_path = output_dir / "exposure_with_pop_and_flood.csv"

print(f"Processing population and flood depth allocation for event: {EVENT_NAME}")
print(f"Base directory: {base_dir}")

# ===== LOAD DATA =====
print("Loading building footprints from spatial.fgb...")
try:
    buildings_fp = gpd.read_file(buildings_fp_path)
    
    # Filter to only buildings with damage (if not already filtered)
    if 'total_damage' in buildings_fp.columns:
        original_count = len(buildings_fp)
        buildings_fp = buildings_fp[buildings_fp['total_damage'] > 0]
        filtered_count = len(buildings_fp)
        print(f"Filtered buildings with damage > 0: {filtered_count} out of {original_count}")
    
    print(f"Buildings footprints shape: {buildings_fp.shape}")
    print(f"Buildings footprints columns: {list(buildings_fp.columns)}")
    print(f"Buildings footprints CRS: {buildings_fp.crs}")
    
except FileNotFoundError:
    print(f"ERROR: spatial.fgb not found at {buildings_fp_path}")
    print("Make sure the FIAT model has been run and spatial.fgb has been generated.")
    exit(1)

print("\nLoading exposure data...")
exposure_df = pd.read_csv(exposure_csv_path)
print(f"Exposure data shape: {exposure_df.shape}")
print(f"Exposure data columns: {list(exposure_df.columns)}")

# ===== PREPARE BUILDING DATA =====
print("\nPreparing building data...")

# Check if Shape_Area exists, if not calculate it
if 'Shape_Area' not in buildings_fp.columns:
    # Calculate building areas
    if buildings_fp.crs == 'EPSG:4326':
        print("Converting to projected CRS for area calculation...")
        # Convert to a suitable UTM zone for accurate area calculation
        # For Mozambique, UTM zone 36S (EPSG:32736) is appropriate
        buildings_fp_projected = buildings_fp.to_crs('EPSG:32736')
        buildings_fp_projected['Shape_Area'] = buildings_fp_projected.geometry.area
        # Keep the area in the original dataframe
        buildings_fp['Shape_Area'] = buildings_fp_projected['Shape_Area']
    else:
        buildings_fp['Shape_Area'] = buildings_fp.geometry.area
else:
    print("Using existing Shape_Area column")

# Create building centroids for point-in-polygon operations
print("Creating building centroids...")
buildings_fp['centroid'] = buildings_fp.geometry.centroid
buildings_centroids = buildings_fp.copy()
buildings_centroids['geometry'] = buildings_centroids['centroid']

# Store original CRS for later
original_crs = buildings_centroids.crs

# Ensure we're working in the correct CRS for the population raster
# Most population rasters are in EPSG:4326
if buildings_centroids.crs != 'EPSG:4326':
    print("Converting centroids to EPSG:4326 for raster operations...")
    buildings_centroids = buildings_centroids.to_crs('EPSG:4326')

# ===== PREPARE RESULT DATAFRAME =====
print("Preparing result dataframe...")
result_df = buildings_centroids.copy()
result_df['population'] = 0.0

# ===== POPULATION ALLOCATION =====
print("\n===== STARTING POPULATION ALLOCATION =====")
print(f"Population raster path: {population_raster_path}")

# Check if population raster exists
if not Path(population_raster_path).exists():
    print(f"ERROR: Population raster not found at {population_raster_path}")
    print("Please update the population_raster_path variable with the correct path to your population raster.")
    exit(1)

# ===== CLIP POPULATION RASTER TO BUILDING EXTENT =====
print("Clipping population raster to building extent...")

# Get the bounds of all buildings (with some buffer)
buildings_bounds = buildings_centroids.total_bounds
buffer_deg = 0.01  # ~1km buffer in degrees
clip_bounds = [
    buildings_bounds[0] - buffer_deg,  # minx
    buildings_bounds[1] - buffer_deg,  # miny  
    buildings_bounds[2] + buffer_deg,  # maxx
    buildings_bounds[3] + buffer_deg   # maxy
]

print(f"Building bounds: {buildings_bounds}")
print(f"Clipping bounds (with buffer): {clip_bounds}")

# Load and clip the population raster
with rasterio.open(population_raster_path) as src:
    print(f"Population raster CRS: {src.crs}")
    print(f"Population raster shape: {src.shape}")
    print(f"Population raster bounds: {src.bounds}")
    
    # Create a window for clipping
    from rasterio.windows import from_bounds
    from rasterio.mask import mask
    from shapely.geometry import box as shapely_box
    
    # Create clipping geometry
    clip_geom = [shapely_box(*clip_bounds)]
    
    # Clip the raster
    try:
        clipped_data, clipped_transform = mask(src, clip_geom, crop=True, nodata=src.nodata)
        clipped_data = clipped_data[0]  # Get first band
        
        print(f"Original raster shape: {src.shape}")
        print(f"Clipped raster shape: {clipped_data.shape}")
        
        # Create a temporary profile for the clipped raster
        clipped_profile = src.profile.copy()
        clipped_profile.update({
            'height': clipped_data.shape[0],
            'width': clipped_data.shape[1],
            'transform': clipped_transform
        })
        
    except Exception as e:
        print(f"Warning: Could not clip raster ({e}). Using full raster...")
        clipped_data = src.read(1)
        clipped_transform = src.transform
        clipped_profile = src.profile.copy()
    
    # Get raster statistics
    valid_cells = clipped_data[clipped_data > 0]
    print(f"Total population in clipped raster: {valid_cells.sum():.0f}")
    print(f"Number of cells with population: {len(valid_cells)}")
    
    processed_cells = 0
    total_cells_with_pop = len(valid_cells)
    
    # ===== VECTORIZED PROCESSING (MUCH FASTER) =====
    print("Using vectorized processing for speed...")
    
    # Create spatial index for buildings (major speed improvement)
    print("Creating spatial index for buildings...")
    building_sindex = result_df.sindex
    
    # Get all cells with population > 0
    pop_cells = np.where(clipped_data > 0)
    cell_rows, cell_cols = pop_cells
    cell_populations = clipped_data[pop_cells]
    
    print(f"Processing {len(cell_rows)} cells with population in batches...")
    
    # Process in batches for memory efficiency
    batch_size = 1000
    total_batches = len(cell_rows) // batch_size + (1 if len(cell_rows) % batch_size else 0)
    
    for batch_idx in range(total_batches):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, len(cell_rows))
        
        print(f"Processing batch {batch_idx + 1}/{total_batches} (cells {start_idx}-{end_idx})...")
        
        # Get batch of cells
        batch_rows = cell_rows[start_idx:end_idx]
        batch_cols = cell_cols[start_idx:end_idx]
        batch_populations = cell_populations[start_idx:end_idx]
        
        # Create all cell polygons for this batch
        batch_polygons = []
        for i, (row, col) in enumerate(zip(batch_rows, batch_cols)):
            # Convert to geographic coordinates
            window = rasterio.windows.Window(col, row, 1, 1)
            cell_bounds = rasterio.windows.bounds(window, clipped_transform)
            cell_polygon = box(*cell_bounds)
            batch_polygons.append(cell_polygon)
        
        # Create GeoDataFrame of cell polygons
        cells_gdf = gpd.GeoDataFrame({
            'cell_id': range(start_idx, end_idx),
            'population': batch_populations,
            'geometry': batch_polygons
        }, crs=result_df.crs)
        
        # Spatial join to find buildings in each cell (vectorized!)
        print(f"  Performing spatial join for batch {batch_idx + 1}...")
        buildings_in_cells = gpd.sjoin(result_df, cells_gdf, how='inner', predicate='intersects')
        
        if not buildings_in_cells.empty:
            # Group by cell to calculate proportions
            cell_groups = buildings_in_cells.groupby('cell_id')
            
            for cell_id, cell_buildings in cell_groups:
                cell_population = batch_populations[cell_id - start_idx]
                
                # Calculate total area for buildings in this cell
                total_area_in_cell = cell_buildings['Shape_Area'].sum()
                
                if total_area_in_cell > 0:
                    # Calculate proportional population for each building
                    proportions = cell_buildings['Shape_Area'] / total_area_in_cell
                    pop_allocation = proportions * cell_population
                    
                    # Add to result (using building indices)
                    building_indices = cell_buildings.index
                    result_df.loc[building_indices, 'population'] += pop_allocation.values

print(f"Population allocation completed!")
print(f"Total population allocated: {result_df['population'].sum():.0f}")

# ===== MERGE WITH EXPOSURE DATA =====
print("\n===== MERGING WITH EXPOSURE DATA =====")

# Determine the common key for merging
# Common keys might be: 'object_id', 'Object ID', 'id', 'ID', etc.
building_keys = [col for col in result_df.columns if 'object' in col.lower() or col.lower() in ['id', 'fid']]
exposure_keys = [col for col in exposure_df.columns if 'object' in col.lower() or col.lower() in ['id', 'fid']]

print(f"Potential building keys: {building_keys}")
print(f"Potential exposure keys: {exposure_keys}")

# Try to find a common key
common_key = None
for bkey in building_keys:
    for ekey in exposure_keys:
        if bkey.lower() == ekey.lower():
            common_key = (bkey, ekey)
            break
    if common_key:
        break

if common_key:
    print(f"Using common key: {common_key[0]} (buildings) -> {common_key[1]} (exposure)")
    
    # Select columns to merge
    merge_columns = ['population', common_key[0]]
    if 'relative_damage' in result_df.columns:
        merge_columns.append('relative_damage')
    
    # Merge the data with exposure data
    exposure_with_updates = exposure_df.merge(
        result_df[merge_columns], 
        left_on=common_key[1], 
        right_on=common_key[0], 
        how='left'
    )
    
    # Fill NaN values with 0
    exposure_with_updates['population'] = exposure_with_updates['population'].fillna(0)
    if 'relative_damage' in exposure_with_updates.columns:
        exposure_with_updates['relative_damage'] = exposure_with_updates['relative_damage'].fillna(0)
    
else:
    print("WARNING: Could not find a common key for merging!")
    print("You may need to manually specify the key columns.")
    
    # Create a simple merge based on index if no common key is found
    exposure_with_updates = exposure_df.copy()
    if len(exposure_df) == len(result_df):
        exposure_with_updates['population'] = result_df['population'].values
        if 'relative_damage' in result_df.columns:
            exposure_with_updates['relative_damage'] = result_df['relative_damage'].values
    else:
        exposure_with_updates['population'] = 0

# ===== SAVE RESULTS =====
print("\n===== SAVING RESULTS =====")

# Clean up extra geometry columns before saving
print("Cleaning up geometry columns...")
columns_to_keep = []
geometry_columns = []

for col in result_df.columns:
    if hasattr(result_df[col], 'geom_type'):  # This is a geometry column
        geometry_columns.append(col)
    else:
        columns_to_keep.append(col)

print(f"Found geometry columns: {geometry_columns}")

# Keep only the main geometry column for saving
if 'geometry' in result_df.columns:
    # Use the original geometry column (building footprints)
    result_df_clean = result_df[columns_to_keep + ['geometry']].copy()
    result_df_clean = result_df_clean.set_geometry('geometry')
else:
    # Use the first geometry column found
    main_geom_col = geometry_columns[0]
    result_df_clean = result_df[columns_to_keep + [main_geom_col]].copy()
    result_df_clean = result_df_clean.set_geometry(main_geom_col)

print(f"Final columns for saving: {list(result_df_clean.columns)}")

# Save buildings with all updates
result_df_clean.to_file(output_buildings_path, driver='FlatGeobuf')
print(f"Spatial data with population and flood depth saved to: {output_buildings_path}")

# Save exposure with all updates
exposure_with_updates.to_csv(output_exposure_path, index=False)
print(f"Exposure with population and flood depth saved to: {output_exposure_path}")

# ===== SANITY CHECK =====
print("\n===== SANITY CHECK =====")
# Use clipped data for sanity check
total_population_raster = valid_cells.sum()  # This was calculated from clipped_data above
total_population_allocated = result_df['population'].sum()
ratio = total_population_allocated / total_population_raster if total_population_raster > 0 else 0

print(f"Total population from clipped raster: {total_population_raster:.0f}")
print(f"Total population allocated: {total_population_allocated:.0f}")
print(f"Allocation ratio: {ratio:.4f}")

if 0.9 <= ratio <= 1.1:
    print("Population sanity check PASSED: Population allocation is within acceptable range")
else:
    print("Population sanity check WARNING: Population allocation differs from raster total.")
    print("This is expected when using pre-filtered buildings (spatial.fgb) that only include flood-affected areas.")

# ===== SUMMARY =====
print("\n===== FINAL SUMMARY =====")
print(f"Buildings processed: {len(result_df)}")
print(f"Buildings with population > 0: {len(result_df[result_df['population'] > 0])}")
print(f"Buildings with both population and flood depth > 0: {len(result_df[(result_df['population'] > 0) & (result_df['inun_depth'] > 0)])}")

# Damage statistics (if available)
if 'total_damage' in result_df.columns:
    print(f"\nDamage statistics:")
    print(f"Total damage: ${result_df['total_damage'].sum():,.0f}")
    print(f"Average damage per building: ${result_df['total_damage'].mean():,.0f}")

# Population statistics
print(f"\nPopulation statistics:")
print(f"Average population per building: {result_df['population'].mean():.2f}")
print(f"Maximum population per building: {result_df['population'].max():.2f}")
print(f"Total affected population: {result_df['population'].sum():,.0f}")

# Flood statistics  
print(f"\nFlood depth statistics:")
print(f"Average flood depth per building (where > 0): {result_df[result_df['inun_depth'] > 0]['inun_depth'].mean():.2f} m")
print(f"Maximum flood depth: {result_df['inun_depth'].max():.2f} m")

print(f"\nOutput files saved to: {output_dir}")
# %%

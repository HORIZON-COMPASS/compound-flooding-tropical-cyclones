#%% -*- coding: utf-8 -*-
# Obtained from https://github.com/VU-IVM/GLOPOP-S/blob/main/READ_SYNTHPOP_DATA/read_synthpop_data.py
"""
Created on Mon Dec 18 17:32:30 2023

@author: Marijn Ton
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import gzip
import rasterio

#%%
# Data available on  https://doi.org/10.7910/DVN/KJC3RH

# Filename .dat file.
gdlcode = 'MOZr107' #example AFGr101. For all gdlcodes see Nr_individuals_data_availability.csv
filename = 'synthpop_' + gdlcode + '_grid.dat.gz' 
# filepath = Path('data', filename)
filepath = Path('data', 'GLOPOP-SG', filename)

with gzip.open(filepath, 'rb') as f:
    # Read the binary content of the file
    binary_content = f.read()

    data_np = np.frombuffer(binary_content, dtype=np.int32)

n_columns = 16

n_people = data_np.size// n_columns
# reshapa data
data_reshaped = np.reshape(data_np, (n_columns, n_people)).transpose()

attribute_names = ['HID', 'RELATE_HEAD', 'INCOME', 'WEALTH', 'RURAL', 'AGE', 'GENDER', 'EDUC', 
                   'HHTYPE', 'HHSIZE_CAT','AGRI_OWNERSHIP', 'FLOOR', 'WALL', 'ROOF', 'SOURCE', 'GRID_NR']

df = pd.DataFrame(data_reshaped, columns=attribute_names)

# column names explained:
# HID: household ID. Be aware that household ID's are unique within a region, not within a country. 
# RELATE_HEAD: relationship to household head. 1: head. 2: partner. 3: child. 4: relative. 5: non-relative. 
# INCOME: 1: poorest 20% individuals, 2: poorer 20%, 3: middle 20%, 4: richer 20%, 5: richest 20% individuals, -1: unavailable for country 
# WEALTH: 1: poorest 20% individuals, 2: poorer 20%, 3: middle 20%, 4: richer 20%, 5: richest 20% individuals, -1: unavailable for country 
# RURAL: settlement type. 1: rural. 0: urban. 
# AGE: age groups. 1:0-4, 2: 5-14, 3: 15-24, 4: 25-34, 5: 35-44, 6: 45-54, 7:55-64, 8:65+
# GENDER: gender. 1: male. 0: female. 
# EDUC: highest achieved education level. 1: less than primary. 2: primary completed. 
# 3: less than secondary. 4: secondary completed. 5: higher education
# HHTYPE: household types. 1: single. 2: couple. 3: couple with children. 4: one parent with children
# 5: couple with (non-)relatives. 6: couple with children and (non-)relatives. 
# 7: one parent with children and (non-)relatives. 8: other
# HHSIZE_CAT: household size categories. 1: 1 person. 2: 2 persons. 3: 3-4 persons. 4: 5-6. 5:7-10. 6: 10+. 
# AGRI_OWNERSHIP: Household owns land for agriculture. 1: yes. 0: no, -1: unavailable for country
# FLOOR: floor material. 1: natural, 2: rudimentary, 3: finished, -1: unavailable for country
# WALL: wall material. 1: natural, 2: rudimentary, 3: finished, -1: unavailable for country
# ROOF: roof material. 1: natural, 2: rudimentary, 3: finished, -1: unavailable for country
# SOURCE: data source of synthetic population. 1: LIS, 2: LIS survey, 3: LIS marginals, 
# 4: Modeled by LIS data, 5: DHS, 6: Modeled by DHS data
# GRID_NR: grid cell number


# %%
# Read corresponding .tif file to add coordinates of the grid cells
tif_file = gdlcode + '_grid_nr.tif'
tif_filepath = Path('data', 'GLOPOP-SG', tif_file)

# %%
# Open the TIFF
with rasterio.open(tif_filepath) as src:
    grid_ids = src.read(1)                   # read first band
    transform = src.transform            
    nodata = src.nodata
    profile = src.profile

# Generate row/col indices
rows, cols = np.indices(grid_ids.shape)
xs, ys = rasterio.transform.xy(transform, rows, cols)

# Flatten arrays
df_grid = pd.DataFrame({
    "grid_id": grid_ids.ravel(),
    "x": np.array(xs).ravel(),
    "y": np.array(ys).ravel()
})

# %%
# Now add the coordinate information to the synthetic pop data
df_merged = df.merge(
    df_grid[["grid_id", "x", "y"]],
    how="left",
    left_on="GRID_NR",
    right_on="grid_id"
)

# Save to csv
df_merged.to_csv(Path('data', 'GLOPOP-SG', f'synthpop_{gdlcode}_grid_combined.csv'), index=False)

# Aggregate population count per grid cell
df_pop = df_merged.groupby("GRID_NR").size().reset_index(name="pop_count").rename(columns={"GRID_NR": "grid_id"})


# %%
# Create population raster
pop_raster = np.zeros_like(grid_ids, dtype="float32")

# Create lookup dictionary
pop_dict = dict(zip(df_pop["grid_id"], df_pop["pop_count"]))

# Fill raster
# Vectorized assignment: replace grid_ids with population counts
vectorized_pop = np.vectorize(lambda gid: pop_dict.get(gid, 0))
pop_raster = vectorized_pop(grid_ids)

# Optionally set 0 for NoData pixels (grid_id == 0)
pop_raster[grid_ids == 0] = 0

# 5. Write population raster
profile.update(dtype="float32", nodata=0)

output_path = f"data/GLOPOP-SG/{gdlcode}_population.tif"
with rasterio.open(output_path, "w", **profile) as dst:
    dst.write(pop_raster, 1)

# %%

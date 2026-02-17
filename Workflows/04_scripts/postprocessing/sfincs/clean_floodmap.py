# Author: n-aleksandrova
# Created: 22-05-2025
# This script cleans the raster file of the flood map  to remove clusters of pixels
#  (10 pixels, but can be changed in the command below: -st 10)

# To run this script you need to install a python environment with rasterio:
#    conda create -n raster_env python=3.8 rasterio
#    conda activate raster_env

# %%
import os
import rasterio
from rasterio.features import sieve
import numpy as np

input_dir = r"P:\11210471-001-compass\03_Runs\test\Kenneth\sfincs\event_precip_era5_hourly\plot_output"
out_dir = os.path.join(input_dir, "floodmaps_cleaned")
os.makedirs(out_dir, exist_ok=True)

# %%
for file in os.listdir(input_dir):
    if file.endswith(".tif"):
        input = os.path.join(input_dir, file)
        filename = os.path.splitext(os.path.basename(file))[0]

        src = rasterio.open(input, mode="r")
        msk = src.read_masks()
        msk2 = sieve(msk, size=150, connectivity=4)

        raster = src.read()
        raster[0] = np.where(msk2 == 0, np.nan, raster[0])

        profile = src.meta
        file_out = os.path.join(
            input_dir, "floodmaps_cleaned", filename + "_filtered.tif"
        )
        with rasterio.open(file_out, "w", **profile) as dst:
            dst.write(raster)

        src.close()


# %%

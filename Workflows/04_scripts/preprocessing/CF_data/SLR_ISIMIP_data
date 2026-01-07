# %% In this script the regional SLR is calculate by using ISIMIP data from Treu et al. (2023): https://data.isimip.org/search/query/10.48364/ISIMIP.749905.1/
# Use the 'compass-snake-dfm' python environment
# Import the necessary packages
import numpy as np
import xarray as xr
import yaml
import os
import requests
import platform

prefix = "p:/" if platform.system() == "Windows" else "/p/"

# Load config file
config_file = "../../../01_config_snakemake/config_general_MZB.yml"
with open(config_file, "r") as f:
    config = yaml.safe_load(f)
    config = config['runname_ids']['Idai']

# %%
# Region of case study area in Mozambique
lat_MZB = [-20.12, -19.30]
lon_MZB = [34.33, 34.95]

# Time period for linear trend calculation
years = range(1985, 2016)

#%%
datasets = {
    "obsclim_geo": "hcc_obsclim_geocentricwaterlevel_global_hourly",
    "counterclim_geo": "hcc_counterclim_geocentricwaterlevel_global_hourly",
    "obsclim_wl": "hcc_obsclim_waterlevel_global_hourly",
    "counterclim_wl": "hcc_counterclim_waterlevel_global_hourly",
}


#%%
# Download and subset ISIMIP SLR data for MZB region
BASE_URL = "https://files.isimip.org/ISIMIP3a/InputData/climate/sealevel/"
OUTDIR = os.path.join(prefix, "11210471-001-compass/01_Data/ISIMIP/SLR/MZB_subset")
TMPDIR = os.path.join(prefix, "11210471-001-compass/01_Data/ISIMIP/SLR/tmp_global")

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(TMPDIR, exist_ok=True)

years = (2013, 2014)  # only download until 2015 for trend calculation

for key, prefix in datasets.items():
    # Determine clim for the URL path
    if "obsclim" in key:
        clim = "obsclim"
    elif "counterclim" in key:
        clim = "counterclim"
    else:
        raise ValueError(f"Cannot determine clim for {key}")

    for year in years:
        # Build URL and filenames
        fname = f"{clim}/global/hourly/historical/HCC/{prefix}_{year}.nc"
        url = BASE_URL + fname
        out_file = os.path.join(OUTDIR, f"{prefix}_MZB_{year}.nc")

        # --- skip if already exists ---
        if os.path.exists(out_file):
            print(f"Skipping {out_file}, already exists.")
            continue

        tmp_file = os.path.join(TMPDIR, f"{prefix}_{year}.nc")
        print(f"Processing {url}")

        # --- download using requests ---
        response = requests.get(url, stream=True)
        if response.status_code != 200:
            raise ValueError(f"Failed to download {url}, status code: {response.status_code}")

        os.makedirs(os.path.dirname(tmp_file), exist_ok=True)
        with open(tmp_file, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        # --- subset to your region ---
        ds = xr.open_dataset(tmp_file)
        mask = ((ds['lat'] >= lat_MZB[0]) & (ds['lat'] <= lat_MZB[1]) &
                (ds['lon'] >= lon_MZB[0]) & (ds['lon'] <= lon_MZB[1]))
        ds_sel = ds.sel(stations=mask)

        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        ds_sel.to_netcdf(out_file)

        ds.close()
        ds_sel.close()

        # --- delete the temporary global file ---
        os.remove(tmp_file)
        print(f"Saved subset to {out_file}")

# %%
# Combine yearly files into a single file for each dataset
def combine_yearly_files(prefix, years, indir, outdir, out_name=None, time_dim='time'):
    files = [os.path.join(indir, f"{prefix}_MZB_{year}.nc") for year in years]
    files = [f for f in files if os.path.exists(f)]

    if not files:
        print(f"No files found for {prefix} in {indir}")
        return None

    print(f"Combining {len(files)} files for {prefix}...")
    ds_list = [xr.open_dataset(f) for f in files]
    combined = xr.concat(ds_list, dim=time_dim)
    
    # Close the original datasets
    for ds in ds_list:
        ds.close()
    
    os.makedirs(outdir, exist_ok=True)
    outfile = out_name or os.path.join(outdir, f"{prefix}_all_years.nc")
    combined.to_netcdf(outfile)
    combined.close()
    print(f"Saved combined file to {outfile}")
    return outfile

years = range(1985, 2016)
indir = OUTDIR
outdir = os.path.join(OUTDIR, "MZB_combined")

for key, prefix in datasets.items():
    combine_yearly_files(prefix, years, indir, outdir)
#%% Use pixi environment compass-dataviz
import os
import cdsapi
import platform

#%% --------------------------------------------------
# Settings
bbox = [-17, 32, -21, 36]   # Mozambique region
years = range(1988, 2020)
months = range(1, 13)
days = [f"{d:02d}" for d in range(1, 32)]

prefix = "p:/" if platform.system() == "Windows" else "/p/"
outdir = os.path.join(prefix, "11210471-001-compass", "01_Data", "glofas_v4")
os.makedirs(outdir, exist_ok=True)

dataset = "cems-glofas-historical"

client = cdsapi.Client(url="https://ewds.climate.copernicus.eu/api")

#%% --------------------------------------------------
# functions
def yearly_request(year):
    return {
        "system_version": ["version_4_0"],
        "hydrological_model": ["lisflood"],
        "product_type": ["consolidated"],
        "variable": ["river_discharge_in_the_last_24_hours"],
        "hyear": [str(year)],
        "hmonth": [f"{m:02d}" for m in months],
        "hday": days,
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": bbox,
    }

def monthly_request(year, month):
    return {
        "system_version": ["version_4_0"],
        "hydrological_model": ["lisflood"],
        "product_type": ["consolidated"],
        "variable": ["river_discharge_in_the_last_24_hours"],
        "hyear": [str(year)],
        "hmonth": [f"{month:02d}"],
        "hday": days,
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": bbox,
    }

#%% --------------------------------------------------
# Download loop-
for year in years:
    print(f"\n📦 Processing {year}")

    year_file = os.path.join(outdir, f"glofas4_reanalysis_{year}.nc")

    # --- Try yearly download ---
    try:
        print("  → Trying yearly download")
        client.retrieve(dataset, yearly_request(year), year_file)
        print(f"Saved {year_file}")
        continue

    except Exception as e:
        print(f"Yearly download failed: {e}")
        print("Falling back to monthly downloads")

    # --- Monthly fallback ---
    for month in months:
        month_file = os.path.join(
            outdir, f"glofas4_reanalysis_{year}_{month:02d}.nc"
        )

        if os.path.exists(month_file):
            print(f"{month_file} already exists")
            continue

        try:
            print(f"Downloading {year}-{month:02d}")
            client.retrieve(dataset, monthly_request(year, month), month_file)
            print(f"Saved {month_file}")

        except Exception as e:
            print(f"Failed {year}-{month:02d}: {e}")

# %%

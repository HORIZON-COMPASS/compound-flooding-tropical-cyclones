# %%
# Clip the CEH-GEAR hourly rainfall data (factual + counterfactuals T1/T3) to the
# Somerset region and export as plain NetCDF files for the Met Office.
# This is the forcing-data companion to export_data_somerset_metoffice_v1.py.
# Source files are the local P: drive copies; no Azure resolver needed.
# Output lands in forMetOfficeUK_v1/data_forcing/ alongside data_static/.
#
# The catalog root is set to the ceh-gear source directory so that each per-file
# entry produces a clean output filename (e.g. CEH-GEAR-1hr-v2_201312.nc) rather
# than a glob-named merged file. After export, the 9 per-month catalog entries are
# consolidated back into 3 glob-based entries matching the main catalog structure.
import os
import tempfile
from os.path import join
from pathlib import Path

import yaml
import hydromt

path_repo = str(Path(__file__).resolve().parents[2])
path_dataexport = "/p/11210471-001-compass/01_Data/forMetOfficeUK_v1"
config_file = join(path_repo, "Workflows", "01_config_snakemake", "config_general_somerset_dec2013.yml")

# %%
with open(config_file) as f:
    cfg = yaml.safe_load(f)
bbox_sfincs = eval(cfg["runname_ids"]["SomersetLevels_dec"]["bbox_sfincs"])

# Same ±1° widening as the static export script
bbox_wide = [bbox_sfincs[0] - 1, bbox_sfincs[1] - 1, bbox_sfincs[2] + 1, bbox_sfincs[3] + 1]
print(f"Export bbox: {bbox_wide}")

# %%
# Build a temporary catalog with one entry per file (no globs) so that
# export_data() writes 9 cleanly named output files, not 3 glob-named ones.
src_dir = "/p/11210471-001-compass/01_Data/ceh-gear"
months = ["201312", "201401", "201402"]

variants = {
    "ceh_gear_compass": {
        "prefix": "CEH-GEAR-1hr-v2",
        "notes": "CEH-GEAR 1hr v2 gridded hourly rainfall over the UK (COMPASS local copy); rainfall_amount in kg m-2 (= mm) renamed to precip.",
    },
    "ceh_t1_compass": {
        "prefix": "ceh_t1",
        "notes": "CEH-GEAR 1hr v2 counterfactual variant T1 (COMPASS local copy); rainfall_amount in kg m-2 (= mm) renamed to precip.",
    },
    "ceh_t3_compass": {
        "prefix": "ceh_t3",
        "notes": "CEH-GEAR 1hr v2 counterfactual variant T3 (COMPASS local copy); rainfall_amount in kg m-2 (= mm) renamed to precip.",
    },
}

driver_options = {
    "chunks": {"latitude": 300, "longitude": 200, "time": 48},
    "combine": "by_coords",
    "decode_times": True,
    "parallel": True,
}

temp_catalog = {"meta": {"root": src_dir}}
entry_names = []

for key, info in variants.items():
    for month in months:
        entry_key = f"{key}_{month}"
        entry_names.append(entry_key)
        temp_catalog[entry_key] = {
            "data_type": "RasterDataset",
            "uri": f"{info['prefix']}_{month}.nc",
            "data_adapter": {"rename": {"rainfall_amount": "precip"}},
            "driver": {"name": "raster_xarray", "options": driver_options},
            "metadata": {
                "crs": 4326,
                "unit": "mm",
                "category": "meteo",
                "notes": info["notes"],
                "source_url": "https://catalogue.ceh.ac.uk/documents/fc9423d6-3d54-467f-bb2b-fc7357a3941f",
                "source_version": "CEH-GEAR-1hr-v2",
            },
        }

# %%
# Write the temporary catalog to a file, run the export, then clean up.
with tempfile.NamedTemporaryFile(mode="w", suffix=".yml", delete=False) as f:
    yaml.dump(temp_catalog, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
    temp_catalog_path = f.name

try:
    dc = hydromt.data_catalog.DataCatalog(data_libs=[temp_catalog_path])
    dc.export_data(
        new_root=join(path_dataexport, "data_forcing"),
        source_names=entry_names,
        bbox=bbox_wide,
    )
finally:
    os.unlink(temp_catalog_path)

# %%
# Post-process: consolidate the 9 per-month entries in the exported catalog into
# 3 glob-based entries matching the main catalog structure (same pattern as the
# basin_atlas layer-name fix in export_data_somerset_metoffice_v1.py).
catalog_path = join(path_dataexport, "data_forcing", "data_catalog.yml")
with open(catalog_path) as f:
    exported = yaml.safe_load(f)

consolidated = {"meta": exported["meta"]}
for key, info in variants.items():
    template = dict(exported[f"{key}_201312"])
    # Replace the month-specific filename with a glob pattern
    template["uri"] = template["uri"].replace("201312", "*")
    consolidated[key] = template

with open(catalog_path, "w") as f:
    yaml.dump(consolidated, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

print("Forcing export complete.")
print(f"Output: {join(path_dataexport, 'data_forcing')}")

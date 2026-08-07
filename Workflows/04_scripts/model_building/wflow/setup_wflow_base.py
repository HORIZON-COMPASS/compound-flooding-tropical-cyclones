# %%
# Build a Wflow base model with HydroMT v1 / hydromt_wflow 1.0.x (WflowSbmModel).
# Migrated from the v0 setup_*/opt API: the build config is now a `steps` list
# (wflow_base_build_v1.yml); region, river_upa and the SFINCS-derived gauges are injected
# into the steps at runtime before mod.build(steps=...).
from os.path import join, exists
import os
import sys
from pathlib import Path
import yaml
import geopandas as gpd
from hydromt_wflow import WflowSbmModel

# Ensure the shared scripts root is importable regardless of Snakemake CWD
_scripts_dir = str(Path(__file__).parents[2])
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from utils.pipeline_utils import find_steps, find_step_index


# %%
if "snakemake" in locals():
    model_dir        = snakemake.params.dir_model
    config_file      = snakemake.input.config_file
    data_cat         = snakemake.params.data_cat
    region_geom      = snakemake.input.region_geom
    dir_sfincs_model = snakemake.input.dir_sfincs_model
    river_upa        = snakemake.params.river_upa
    # CF_landuse is a path wildcard here (models live in wflow_{CF_landuse}/), so the land-use
    # handle comes from the wildcard; only the reclass table is passed as a param.
    landuse          = getattr(snakemake.wildcards, 'CF_landuse', None)
    lulc_mapping     = snakemake.params.get('lulc_mapping_wflow', None)
else:
    model_dir        = "/p/11210471-001-compass/02_Models/somerset/SomersetLevels/wflow_v1_test"
    config_file      = "../../../05_config_models/01_wflow/wflow_base_build_v1.yml"
    data_cat         = [
        "../../../03_data_catalogs/datacatalog_general_v1___linux.yml",
        "../../../03_data_catalogs/datacatalog_CF_forcing_v1___linux.yml",
    ]
    region_geom      = "/p/11210471-001-compass/02_Models/somerset/SomersetLevels/sfincs_v1/gis/region.geojson"
    dir_sfincs_model = "/p/11210471-001-compass/02_Models/somerset/SomersetLevels/sfincs_v1"
    river_upa        = 30
    landuse          = "vito"
    lulc_mapping     = "vito_mapping_wflow"

# Check whether model folder exists
if not exists(model_dir):
    os.makedirs(model_dir)

# %%
# Read the v1 build config (modeltype / global / steps)
with open(config_file) as f:
    cfg = yaml.safe_load(f)
steps = cfg["steps"]

# Read SFINCS region
region = gpd.read_file(region_geom).to_crs(epsg="4326")

# %%
# Inject dynamic values into the matching steps
for s in find_steps(steps, "setup_basemaps"):
    s["region"] = {"basin": region}
for s in find_steps(steps, "setup_rivers"):
    s["river_upa"] = river_upa

# land use -> setup_lulcmaps. Without lulc_mapping_fn hydromt_wflow falls back to its own built-in
# reclassification table rather than the project's calibrated one, so pass it whenever the config
# declares it. Entries the config leaves unset keep the build yml's defaults.
for s in find_steps(steps, "setup_lulcmaps"):
    if landuse:
        s["lulc_fn"] = landuse
    if lulc_mapping:
        s["lulc_mapping_fn"] = lulc_mapping

# Add a setup_gauges step based on the SFINCS inflow river points. Insert it before
# setup_config_output_timeseries (which references the gauge map by name).
# NOTE: in hydromt_sfincs v2 the river-inflow source points are written to gis/dis.geojson
# (the v0 gis/src.geojson was renamed).
gauges_step = {
    "setup_gauges": {
        "gauges_fn": join(dir_sfincs_model, "gis", "dis.geojson"),
        "snap_to_river": True,
        "snap_uparea": True,
        "rel_error": 0.2,
        "derive_subcatch": False,
        "index_col": "index",
        "basename": "locs",
    }
}
insert_at = find_step_index(steps, "setup_config_output_timeseries", default=len(steps))
steps.insert(insert_at, gauges_step)

# %%
# Build the model (v1: WflowSbmModel, no logger kwarg, build takes steps=)
mod = WflowSbmModel(root=model_dir, data_libs=data_cat, mode="w+")
mod.build(steps=steps)

# Disable CSV output (v0 did mod.config['csv'] = None) and persist the config
if "csv" in mod.config.data:
    mod.config.data.pop("csv", None)
mod.config.write()

# %%

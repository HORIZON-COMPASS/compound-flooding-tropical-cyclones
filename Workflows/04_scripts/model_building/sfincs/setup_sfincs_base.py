# %%
# Build a SFINCS base model with HydroMT v1 / hydromt_sfincs v2 (component-based steps).
# Migrated from the v0 setup_*/opt API: the build config is now a `steps` list
# (sfincs_base_build_v1.yml) and we inject the dynamic values into the matching steps
# before calling mod.build(steps=...).
from os.path import join, exists
import os
import ast
import sys
from pathlib import Path
import yaml
import geopandas as gpd
import hydromt
from hydromt_sfincs import SfincsModel

# Ensure the shared scripts root is importable regardless of Snakemake CWD
_scripts_dir = str(Path(__file__).parents[2])
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from utils.pipeline_utils import find_steps


def get_local_vector_data(file, bbox, data_cat):
    dataCat = hydromt.data_catalog.DataCatalog(data_cat)
    return dataCat.get_geodataframe(data_like=file, bbox=bbox)


# %%
if "snakemake" in locals():
    model_dir        = snakemake.params.dir_model_sfincs
    config_file      = snakemake.input.config_file
    data_cats        = snakemake.params.data_cats
    bbox             = ast.literal_eval(snakemake.params.arg_bbox)
    bathy            = snakemake.params.bathy
    dfm_coastal_mask = snakemake.params.dfm_coastal_mask
    river_upa        = snakemake.params.river_upa
    landuse          = snakemake.params.get('landuse', None)
    lulc_mapping     = snakemake.params.get('lulc_mapping_sfincs', None)
else:
    # Durban precip-only test case (config_durban_floods_2022.yml)
    model_dir = "/p/11210471-001-compass/02_Models/durban/Durban2022/sfincs"
    config_file = "../../../05_config_models/02_sfincs/sfincs_base_build_v1.yml"
    data_cats = [
        "../../../03_data_catalogs/datacatalog_general_v1___linux.yml",
        "../../../03_data_catalogs/datacatalog_SFINCS_obspoints_v1___linux.yml",
        "../../../03_data_catalogs/datacatalog_SFINCS_coastal_coupling_v1___linux.yml",
    ]
    bbox = [30.659688, -29.978273, 31.076825, -29.740075]
    bathy = "gebco2024_MZB"
    dfm_coastal_mask = "coastal_coupling_msk_MZB"
    river_upa = 30
    landuse = "vito"
    lulc_mapping = "vito_mapping_sfincs"

# Check whether model folder exists. If not, make one
if not exists(model_dir):
    os.makedirs(model_dir)

# %%
# Read the v1 build config (modeltype / global / steps)
with open(config_file) as f:
    cfg = yaml.safe_load(f)
steps = cfg["steps"]

# %%
# Build the model region from the basin atlas, clipped to the bbox
region = get_local_vector_data(
    file="basin_atlas_level12_v10",
    bbox=bbox,
    data_cat=data_cats[0],
)

# %%
# Inject the dynamic values into the matching steps
# region -> grid.create_from_region
for s in find_steps(steps, "grid.create_from_region"):
    s["region"] = {"geom": region}

# bathy -> elevation.create and subgrid.create elevation_list
bathy_entry = {"elevation": bathy, "reproj_method": "bilinear"}
for key in ("elevation.create", "subgrid.create"):
    for s in find_steps(steps, key):
        s.setdefault("elevation_list", []).append(bathy_entry)

# dfm coastal mask -> first mask.create_active (exclude) and mask.create_boundary (include)
active_steps = find_steps(steps, "mask.create_active")
if active_steps:
    active_steps[0]["exclude_polygon"] = dfm_coastal_mask
for s in find_steps(steps, "mask.create_boundary"):
    s["include_polygon"] = dfm_coastal_mask

# land use -> subgrid.create roughness_list (v2 renamed v0's datasets_rgh). Only the entries the
# config actually declares are overridden, so a config without CF_landuse / lulc_mapping_sfincs
# keeps the defaults written in the build yml.
if landuse or lulc_mapping:
    for s in find_steps(steps, "subgrid.create"):
        rgh = s.setdefault("roughness_list", [])
        if not rgh:
            rgh.append({})
        if landuse:
            rgh[0]["lulc"] = landuse
        if lulc_mapping:
            rgh[0]["reclass_table"] = lulc_mapping

# river upstream-area threshold -> river inflow (outflow handled after build, see below)
for s in find_steps(steps, "rivers.create_river_inflow"):
    s["river_upa"] = river_upa

# %%
# Initialise and build the model (v1: no logger kwarg, build takes steps=)
mod = SfincsModel(root=model_dir, data_libs=data_cats, mode="w+")
mod.build(steps=steps)

# River outflow is not a registered build step in hydromt_sfincs 2.0.0rc3, so call the
# component method directly on the built model, then re-write the affected outputs.
# rc3 bug: create_river_outflow references self.logger (missing) in its "no points" branch;
# give the component the attribute so it behaves, and guard the call defensively.
import logging
try:
    mod.rivers.logger = logging.getLogger("hydromt_sfincs")
except Exception:
    pass
try:
    mod.rivers.create_river_outflow(
        hydrography="merit_hydro",
        river_len=5000,
        river_upa=river_upa,
        keep_rivers_geom=True,
    )
    mod.write()
except Exception as e:
    print(f"WARNING: river outflow step skipped (hydromt_sfincs 2.0.0rc3 limitation): {e}")

# %%
# Plot the region and boundaries (model-level method, unchanged)
fig, ax = mod.plot_basemap(
    fn_out=model_dir,
    variable="dep",
    plot_bounds=True,
    plot_geoms=True,
    plot_region=True,
    bmap="sat",
    zoomlevel=12,
    figsize=(8, 6),
)

# %%

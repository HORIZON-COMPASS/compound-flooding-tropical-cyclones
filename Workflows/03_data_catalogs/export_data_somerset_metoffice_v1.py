# %%
# Clip the global static input datasets to the Somerset region and export a
# self-contained data package (+ catalog) for sharing with the Met Office.
# v1 equivalent of forMetOfficeUK_other/subset_data_hydroMT.py (HydroMT v0).
# Static data only for now; forcing data (CEH-GEAR, ERA5, GTSM) is a follow-up.
from os.path import join

import yaml
import hydromt

path_repo = "/u/morenodu/git_repos/compound-flooding-tropical-cyclones"
path_dataexport = "/p/11210471-001-compass/01_Data/forMetOfficeUK_v1"
path_datacat = join(path_repo, "Workflows", "03_data_catalogs")
config_file = join(path_repo, "Workflows", "01_config_snakemake", "config_general_somerset_dec2013.yml")

dc = hydromt.data_catalog.DataCatalog(
    data_libs=[join(path_datacat, "datacatalog_general_v1___linux.yml")]
)

# %%
with open(config_file) as f:
    cfg = yaml.safe_load(f)
bbox_sfincs = eval(cfg["runname_ids"]["SomersetLevels_dec"]["bbox_sfincs"])

# enlarge bbox (same margin as the v0 export script)
bbox_wide = [bbox_sfincs[0] - 1, bbox_sfincs[1] - 1, bbox_sfincs[2] + 1, bbox_sfincs[3] + 1]

# %%
static_data = [
    "fabdem",
    "merit_hydro",
    "merit_hydro_index",
    "rivers_lin2019_v1",
    "hydro_lakes",
    "hydro_reservoirs",
    "rgi",
    "modis_lai",
    "soilgrids",
    "ksathorfrac_global",
    "vito",
    "vito_mapping",
    "gcn250",
    "osm_coastlines",
    "ocean_shape",
    "basin_atlas_level12_v10",
    "gswo",
]

dc.export_data(
    new_root=join(path_dataexport, "data_static"),
    source_names=static_data,
    bbox=bbox_wide,
)

dc.export_data(
    new_root=join(path_dataexport, "data_static"),
    source_names=["gebco"],
    bbox=bbox_wide,
    append=True,
)

# %%
# export_data() writes the basin_atlas gpkg with a generic layer name ("basin_atlas_v10")
# but leaves the catalog's driver option pointing at the original source layer name
# ("BasinATLAS_v10_lev12"), which then fails to open. Patch it post-export.
catalog_path = join(path_dataexport, "data_static", "data_catalog.yml")
with open(catalog_path) as f:
    catalog_text = f.read()
catalog_text = catalog_text.replace(
    "layer: BasinATLAS_v10_lev12", "layer: basin_atlas_v10"
)
with open(catalog_path, "w") as f:
    f.write(catalog_text)

"""
Data loading, population redistribution, exposure calculation, and summary statistics.

Run this script to produce all intermediate arrays and GeoDataFrames
that the plotting script consumes.

Usage:
    python run_population_analysis.py
    (or run cell-by-cell in an interactive editor)
"""
#%%
import json
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
import rioxarray as rxr
from rasterio import features
from rasterio.mask import mask
from rasterio.transform import rowcol
from hydromt import DataCatalog
from shapely.geometry import Polygon

from config import (
    DATACAT_PATH,
    REGION_PATH, BACKGROUND_PATH, SOFALA_PROVINCE_PATH,
    BEIRA_DISTRICT_PATH, DISTRICTS_ADM3_PATH,
    POPULATION_2019_PATH, POPULATION_1990_PATH, SETTLEMENT_TYPE_PATH,
    F_FLOODING_PATH, CF_FLOODING_PATH, SFINCS_SUBGRID_PATH,
    EXPORT_PATH, SUMMARY_TABLE_PATH, DISTRICT_SUMMARY_PATH,
    SFINCS_DIR_F,
    DROP_DISTRICTS, MASK_POLY_COORDS, DEPTH_BINS,
    regridded_pop_path,
)
from population_utils import (
    get_extent,
    reproject_and_redistribute_population_over_land,
    pop_raster_to_gdf,
    aggregate_pop,
    compute_cdf_and_bins,
)

warnings.filterwarnings("ignore")

# ====================================================================== #
# 1. LOAD INPUT DATA                                                      #
# ====================================================================== #
#%%
data_catalog = DataCatalog(data_libs=[DATACAT_PATH])

# --- Vector layers ---
region = gpd.read_file(REGION_PATH)
background = gpd.read_file(BACKGROUND_PATH, driver="GeoJSON")
shapefile_sofala = gpd.read_file(SOFALA_PROVINCE_PATH)

beira_district = gpd.read_file(BEIRA_DISTRICT_PATH)
districts_adm3 = gpd.read_file(DISTRICTS_ADM3_PATH)
districts_adm2 = data_catalog.get_geodataframe("gadm_level2", geom=region, buffer=1000)

# --- Original 1 km population rasters (for comparison / validation) ---
with rasterio.open(POPULATION_2019_PATH) as src_2019:
    region_proj = region.to_crs(src_2019.crs)
    pop_2019_1km, transform_pop_2019_1km = mask(src_2019, region_proj.geometry, crop=True)
    print("No-data value Population 2019:", src_2019.nodata)

with rasterio.open(POPULATION_1990_PATH) as src_1990:
    region_proj = region.to_crs(src_1990.crs)
    pop_1990_1km, transform_pop_1990_1km = mask(src_1990, region_proj.geometry, crop=True)
    print("No-data value Population 1990:", src_1990.nodata)

# ====================================================================== #
# 2. FLOOD GRID SETUP                                                     #
# ====================================================================== #
#%%
with rasterio.open(SFINCS_SUBGRID_PATH) as src:
    flood_grid_crs = src.crs
    flood_grid_transform = src.transform
    flood_grid_shape = (src.height, src.width)

region = region.to_crs(flood_grid_crs)
region_wsg84 = region.to_crs("EPSG:4326")
region_geom = [json.loads(region.to_json())["features"][0]["geometry"]]

# --- Flood rasters ---
hmax_F_da = rxr.open_rasterio(F_FLOODING_PATH).squeeze("band", drop=True)
hmax_CF_da = rxr.open_rasterio(CF_FLOODING_PATH).squeeze("band", drop=True)
hmax_F = hmax_F_da.values
hmax_CF = hmax_CF_da.values
hmax_diff = hmax_F - hmax_CF

flood_extent = get_extent(flood_grid_transform, flood_grid_shape[1], flood_grid_shape[0])

# --- Background layers ---
mask_poly = Polygon(MASK_POLY_COORDS)
bg_filtered = background.copy()
bg_filtered["geometry"] = bg_filtered.geometry.apply(lambda g: g.difference(mask_poly))

background_utm = background.to_crs(flood_grid_crs)
bg_filtered_utm = bg_filtered.to_crs(flood_grid_crs)
region_utm = region.to_crs(flood_grid_crs)
districts_adm3_utm = districts_adm3.to_crs(flood_grid_crs)
districts_adm2_utm = districts_adm2.to_crs(flood_grid_crs)
beira_utm = beira_district.to_crs(flood_grid_crs)

districts_adm3_filtered = districts_adm3_utm[~districts_adm3_utm["NAME_3"].isin(DROP_DISTRICTS)]

# ====================================================================== #
# 3. REDISTRIBUTE POPULATION ONTO FLOOD GRID                              #
# ====================================================================== #
#%%
pop_arrays = {}
pop_sofala_arrays = {}
pop_sofala_districts_adm3 = {}
pop_affine_sofala_districts_adm3 = {}
pop_sofala_districts_adm2 = {}
pop_affine_sofala_districts_adm2 = {}

for year, path in [(1990, POPULATION_1990_PATH), (2019, POPULATION_2019_PATH)]:
    (
        pop_arrays[year],
        pop_sofala_arrays[year],
        transform_sofala_land,
        pop_sofala_districts_adm3[year],
        pop_affine_sofala_districts_adm3[year],
        pop_sofala_districts_adm2[year],
        pop_affine_sofala_districts_adm2[year],
    ) = reproject_and_redistribute_population_over_land(
        pop_path=path,
        land_gdf=background_utm,
        flood_crs=flood_grid_crs,
        flood_transform=flood_grid_transform,
        flood_shape=flood_grid_shape,
        province_geom=shapefile_sofala,
        region=region,
        districts_adm3=districts_adm3_filtered,
        districts_adm2=districts_adm2,
        year=year,
        out_raster_path=str(regridded_pop_path(year)),
    )

# ====================================================================== #
# 4. COMPUTE EXPOSURE                                                      #
# ====================================================================== #
#%%
# --- Pixel-level DataFrames ---
gdf_pop_2019_flood_depth_F = pop_raster_to_gdf(
    pop_arrays[2019], hmax_F, flood_grid_transform,
    year=2019, climate="F", export_df=True, export_path=str(EXPORT_PATH),
)
gdf_pop_2019_flood_depth_CF = pop_raster_to_gdf(
    pop_arrays[2019], hmax_CF, flood_grid_transform,
    year=2019, climate="CF", export_df=True, export_path=str(EXPORT_PATH),
)
gdf_pop_1990_flood_depth_F = pop_raster_to_gdf(
    pop_arrays[1990], hmax_F, flood_grid_transform,
    year=1990, climate="F", export_df=True, export_path=str(EXPORT_PATH),
)
gdf_pop_1990_flood_depth_CF = pop_raster_to_gdf(
    pop_arrays[1990], hmax_CF, flood_grid_transform,
    year=1990, climate="CF", export_df=True, export_path=str(EXPORT_PATH),
)

# --- Simple exposure rasters ---
ra_exposed_pop_2019_F = np.where(hmax_F > 0, pop_arrays[2019], 0)
ra_exposed_pop_2019_CF = np.where(hmax_CF > 0, pop_arrays[2019], 0)
ra_exposed_pop_1990_F = np.where(hmax_F > 0, pop_arrays[1990], 0)
ra_exposed_pop_1990_CF = np.where(hmax_CF > 0, pop_arrays[1990], 0)

# --- Exposed-only subsets ---
gdf_pop_2019_exposed_F = gdf_pop_2019_flood_depth_F[gdf_pop_2019_flood_depth_F["flood_depth"] > 0]
gdf_pop_2019_exposed_CF = gdf_pop_2019_flood_depth_CF[gdf_pop_2019_flood_depth_CF["flood_depth"] > 0]
gdf_pop_1990_exposed_F = gdf_pop_1990_flood_depth_F[gdf_pop_1990_flood_depth_F["flood_depth"] > 0]
gdf_pop_1990_exposed_CF = gdf_pop_1990_flood_depth_CF[gdf_pop_1990_flood_depth_CF["flood_depth"] > 0]

# --- Coarse aggregated grids ---
gdf_pop_2019_exposed_F_coarse = aggregate_pop(
    pop_arrays[2019], hmax_F, flood_grid_transform, flood_grid_crs, region_utm, background_utm
)
gdf_pop_2019_exposed_CF_coarse = aggregate_pop(
    pop_arrays[2019], hmax_CF, flood_grid_transform, flood_grid_crs, region_utm, background_utm
)
gdf_pop_1990_exposed_F_coarse = aggregate_pop(
    pop_arrays[1990], hmax_F, flood_grid_transform, flood_grid_crs, region_utm, background_utm
)
gdf_pop_1990_exposed_CF_coarse = aggregate_pop(
    pop_arrays[1990], hmax_CF, flood_grid_transform, flood_grid_crs, region_utm, background_utm
)

# ====================================================================== #
# 5. UNIFORM POPULATION GROWTH SCENARIO                                    #
# ====================================================================== #
#%%
pop_growth = np.nansum(pop_arrays[2019]) / np.nansum(pop_arrays[1990])
print(f"Uniform population growth from 1990 to 2019 in Sofala region: {pop_growth * 100:.2f}%")

pop_array_uniform_2019 = pop_arrays[1990] * pop_growth

print(f"{np.nansum(pop_array_uniform_2019):,.0f} people in 2019 with uniform growth")
print(f"{np.nansum(pop_arrays[2019]):,.0f} people in 2019 actual")

ra_exposed_pop_2019_F_uniform = np.where(hmax_F > 0, pop_array_uniform_2019, 0)
ra_exposed_pop_2019_CF_uniform = np.where(hmax_CF > 0, pop_array_uniform_2019, 0)

gdf_pop_2019_flood_depth_F_uniform = pop_raster_to_gdf(
    pop_array_uniform_2019, hmax_F, flood_grid_transform,
    year="2019_uniform", climate="F", export_df=True, export_path=str(EXPORT_PATH),
)
gdf_pop_2019_flood_depth_CF_uniform = pop_raster_to_gdf(
    pop_array_uniform_2019, hmax_CF, flood_grid_transform,
    year="2019_uniform", climate="CF", export_df=True, export_path=str(EXPORT_PATH),
)

gdf_pop_2019_exposed_F_uniform = gdf_pop_2019_flood_depth_F_uniform[
    gdf_pop_2019_flood_depth_F_uniform["flood_depth"] > 0
]
gdf_pop_2019_exposed_CF_uniform = gdf_pop_2019_flood_depth_CF_uniform[
    gdf_pop_2019_flood_depth_CF_uniform["flood_depth"] > 0
]

gdf_pop_2019_exposed_F_uniform_coarse = aggregate_pop(
    pop_array_uniform_2019, hmax_F, flood_grid_transform, flood_grid_crs, region_utm, background_utm
)
gdf_pop_2019_exposed_CF_uniform_coarse = aggregate_pop(
    pop_array_uniform_2019, hmax_CF, flood_grid_transform, flood_grid_crs, region_utm, background_utm
)

# ====================================================================== #
# 6. PRINT SUMMARY STATISTICS                                             #
# ====================================================================== #
#%%
def _print_stats(label, year, pop_arrays, pop_sofala, exposed):
    print(f"\n{label}:")
    print(f"  Total population in region: {np.nansum(pop_arrays[year]):,.0f}")
    if pop_sofala is not None:
        print(f"  Total population in Sofala: {np.nansum(pop_sofala[year]):,.0f}")
    n_exposed = int(np.nansum(exposed))
    pct = 100 * n_exposed / np.nansum(pop_arrays[year])
    print(f"  Total exposed people: {n_exposed:,}")
    print(f"  Exposed %: {pct:.2f}%")

_print_stats("2019 Factual", 2019, pop_arrays, pop_sofala_arrays, ra_exposed_pop_2019_F)
_print_stats("2019 Counterfactual", 2019, pop_arrays, pop_sofala_arrays, ra_exposed_pop_2019_CF)
_print_stats("1990 Factual", 1990, pop_arrays, pop_sofala_arrays, ra_exposed_pop_1990_F)
_print_stats("1990 Counterfactual", 1990, pop_arrays, pop_sofala_arrays, ra_exposed_pop_1990_CF)
_print_stats("2019 Uniform", 2019, pop_arrays, None, ra_exposed_pop_2019_F_uniform)

print("\n--- Attribution summary ---")
exp_F = np.nansum(ra_exposed_pop_2019_F)
exp_CF_clim = np.nansum(ra_exposed_pop_2019_CF)
exp_CF_pop = np.nansum(ra_exposed_pop_1990_F)
exp_CF_both = np.nansum(ra_exposed_pop_1990_CF)
exp_uni = np.nansum(ra_exposed_pop_2019_F_uniform)

def _attr(label, diff, base):
    print(f"  {label}: {int(diff):,}  ({100 * diff / base:.2f}%)")

_attr("Climate change", exp_F - exp_CF_clim, exp_F)
_attr("Population change (2019−1990)", exp_F - exp_CF_pop, exp_F)
_attr("Climate & population change", exp_F - exp_CF_both, exp_F)
_attr("Population redistribution (uniform)", exp_F - exp_uni, exp_F)
print(f"  Population growth 1990→2019: {int(np.nansum(pop_arrays[2019]) - np.nansum(pop_arrays[1990])):,}  "
      f"({100 * (np.nansum(pop_arrays[2019]) - np.nansum(pop_arrays[1990])) / np.nansum(pop_arrays[1990]):.2f}%)")

# ====================================================================== #
# 7. SUMMARY TABLE — EXPOSED POP PER FLOOD-DEPTH CATEGORY                 #
# ====================================================================== #
#%%
scenarios = {
    "F": (gdf_pop_2019_exposed_F, 2019),
    "CF_clim": (gdf_pop_2019_exposed_CF, 2019),
    "CF_pop": (gdf_pop_1990_exposed_F, 1990),
    "CF_clim_pop": (gdf_pop_1990_exposed_CF, 1990),
}

total_population = {
    2019: np.nansum(pop_arrays[2019]),
    1990: np.nansum(pop_arrays[1990]),
}

results_abs = {s: [] for s in scenarios}
results_rel = {s: [] for s in scenarios}

for depth_label, (dmin, dmax) in DEPTH_BINS.items():
    for sname, (gdf, year) in scenarios.items():
        m = (gdf["flood_depth"] > dmin) & (gdf["flood_depth"] <= dmax)
        value = np.nansum(gdf.loc[m, "population"])
        results_abs[sname].append(value)
        results_rel[sname].append(value / total_population[year] * 100)

attr_abs, attr_pct = {}, {}
for s in ["CF_pop", "CF_clim", "CF_clim_pop"]:
    attr_abs[s] = [results_abs["F"][i] - results_abs[s][i] for i in range(len(DEPTH_BINS))]
    attr_pct[s] = [
        (attr_abs[s][i] / results_abs["F"][i] * 100) if results_abs["F"][i] != 0 else np.nan
        for i in range(len(DEPTH_BINS))
    ]

rows = []
labels = ["# affected", "% of total", "# attributable", "% attributed"]
for i, depth_label in enumerate(DEPTH_BINS.keys()):
    rows.append({"Depth Range": depth_label, "Label": labels[0], **{s: results_abs[s][i] for s in scenarios}})
    rows.append({"Depth Range": depth_label, "Label": labels[1], **{s: results_rel[s][i] for s in scenarios}})
    rows.append({"Depth Range": depth_label, "Label": labels[2], "F": "-", **{s: attr_abs[s][i] for s in attr_abs}})
    rows.append({"Depth Range": depth_label, "Label": labels[3], "F": "-", **{s: attr_pct[s][i] for s in attr_pct}})

rows.append({
    "Depth Range": "Total population", "Label": "-",
    "F": total_population[2019], "CF_pop": total_population[1990],
    "CF_clim": total_population[2019], "CF_clim_pop": total_population[1990],
})

table_exposed_pop = pd.DataFrame(rows)
table_exposed_pop.to_csv(str(SUMMARY_TABLE_PATH), index=False)
print(table_exposed_pop)

# ====================================================================== #
# 8. DISTRICT-LEVEL EXPOSED POPULATION                                    #
# ====================================================================== #
#%%
dataset_dict = {
    "F": {"exposed": ra_exposed_pop_2019_F, "total": pop_arrays[2019]},
    "CF_clim": {"exposed": ra_exposed_pop_2019_CF, "total": pop_arrays[2019]},
    "CF_pop": {"exposed": ra_exposed_pop_1990_F, "total": pop_arrays[1990]},
    "CF_clim_pop": {"exposed": ra_exposed_pop_1990_CF, "total": pop_arrays[1990]},
}

results_exposed = {k: [] for k in dataset_dict}
results_total = {k: [] for k in dataset_dict}

for _, row in districts_adm3_filtered.iterrows():
    district_mask = features.geometry_mask(
        [row.geometry], out_shape=ra_exposed_pop_2019_F.shape,
        transform=flood_grid_transform, invert=True,
    )
    for key, data in dataset_dict.items():
        results_exposed[key].append(data["exposed"][district_mask].sum())
        results_total[key].append(data["total"][district_mask].sum())

for key in dataset_dict:
    districts_adm3_filtered[f"pop_exposed_{key}"] = results_exposed[key]
    districts_adm3_filtered[f"pop_total_{key}"] = results_total[key]
    districts_adm3_filtered[f"relative_exposed_{key}"] = (
        districts_adm3_filtered[f"pop_exposed_{key}"]
        / districts_adm3_filtered[f"pop_total_{key}"]
        * 100
    )

districts_adm3_filtered["attr_clim"] = (
    (districts_adm3_filtered["pop_exposed_F"] - districts_adm3_filtered["pop_exposed_CF_clim"])
    / districts_adm3_filtered["pop_exposed_F"] * 100
)
districts_adm3_filtered["attr_pop"] = (
    (districts_adm3_filtered["pop_exposed_F"] - districts_adm3_filtered["pop_exposed_CF_pop"])
    / districts_adm3_filtered["pop_exposed_F"] * 100
)
districts_adm3_filtered["attr_clim_pop"] = (
    (districts_adm3_filtered["pop_exposed_F"] - districts_adm3_filtered["pop_exposed_CF_clim_pop"])
    / districts_adm3_filtered["pop_exposed_F"] * 100
)

# --- Total & exposed population per district (for validation) ---
pop_per_district_adm3 = []
pop_exposed_list = []
pop_totals_list = []

for _, row in districts_adm3_filtered.iterrows():
    district_mask = features.geometry_mask(
        [row.geometry], out_shape=ra_exposed_pop_2019_F.shape,
        transform=flood_grid_transform, invert=True,
    )
    pop_exposed_list.append(ra_exposed_pop_2019_F[district_mask].sum())
    pop_totals_list.append(pop_arrays[2019][district_mask].sum())

districts_adm3_filtered["pop_exposed"] = pop_exposed_list
districts_adm3_filtered["pop_total"] = pop_totals_list

with rasterio.open(POPULATION_2019_PATH) as src:
    pop_crs = src.crs

districts_native = districts_adm3_filtered.to_crs(pop_crs)
for _, row in districts_native.iterrows():
    district_mask = features.geometry_mask(
        [row.geometry],
        out_shape=pop_sofala_districts_adm3[2019][0].shape,
        transform=pop_affine_sofala_districts_adm3[2019],
        invert=True,
    )
    pop_per_district_adm3.append(pop_sofala_districts_adm3[2019][0][district_mask].sum())

districts_adm3_filtered["pop_per_district"] = pop_per_district_adm3

# --- Same for ADM2 ---
pop_per_district_adm2 = []
pop_exposed_adm2 = []
pop_totals_adm2 = []

for _, row in districts_adm2.iterrows():
    district_mask = features.geometry_mask(
        [row.geometry], out_shape=ra_exposed_pop_2019_F.shape,
        transform=flood_grid_transform, invert=True,
    )
    pop_exposed_adm2.append(ra_exposed_pop_2019_F[district_mask].sum())
    pop_totals_adm2.append(pop_arrays[2019][district_mask].sum())

districts_adm2["pop_exposed"] = pop_exposed_adm2
districts_adm2["pop_total"] = pop_totals_adm2

districts_native_adm2 = districts_adm2.to_crs(pop_crs)
for _, row in districts_native_adm2.iterrows():
    district_mask = features.geometry_mask(
        [row.geometry],
        out_shape=pop_sofala_districts_adm2[2019][0].shape,
        transform=pop_affine_sofala_districts_adm2[2019],
        invert=True,
    )
    pop_per_district_adm2.append(pop_sofala_districts_adm2[2019][0][district_mask].sum())

districts_adm2["pop_per_district"] = pop_per_district_adm2

# --- District summary table ---
df_district_summary = districts_adm3_filtered[
    ["NAME_3", "pop_per_district", "pop_total", "pop_exposed"]
].copy()
df_district_summary.columns = [
    "District", "District Population (2019)", "Total Population in Region",
    "Exposed Population (2019 Factual)",
]
df_district_summary["% of district pop"] = (
    100 * df_district_summary["Total Population in Region"]
    / df_district_summary["District Population (2019)"]
).round(0)
df_district_summary["% exposed"] = (
    100 * df_district_summary["Exposed Population (2019 Factual)"]
    / df_district_summary["Total Population in Region"]
).round(0)
df_district_summary.to_csv(str(DISTRICT_SUMMARY_PATH), index=False)
print(df_district_summary)

# ====================================================================== #
# 9. BINNED FLOOD-DEPTH DISTRIBUTIONS                                     #
# ====================================================================== #
#%%
bins_fine = np.arange(0, 3.5 + 0.02, 0.01)
low_mask = bins_fine[:-1] < 0.5
mid_mask = (bins_fine[:-1] >= 0.5) & (bins_fine[:-1] < 1.5)
high_mask = bins_fine[:-1] >= 1.5

pop_2019_by_depth_F_fine = compute_cdf_and_bins(gdf_pop_2019_exposed_F, bins_fine)
pop_2019_by_depth_CF_fine = compute_cdf_and_bins(gdf_pop_2019_exposed_CF, bins_fine)
pop_1990_by_depth_F_fine = compute_cdf_and_bins(gdf_pop_1990_exposed_F, bins_fine)
pop_1990_by_depth_CF_fine = compute_cdf_and_bins(gdf_pop_1990_exposed_CF, bins_fine)

bins_coarse = np.arange(0, 3.5 + 0.2, 0.1)
pop_2019_by_depth_F_coarse = compute_cdf_and_bins(gdf_pop_2019_exposed_F, bins_coarse)
pop_2019_by_depth_CF_coarse = compute_cdf_and_bins(gdf_pop_2019_exposed_CF, bins_coarse)
pop_1990_by_depth_F_coarse = compute_cdf_and_bins(gdf_pop_1990_exposed_F, bins_coarse)
pop_1990_by_depth_CF_coarse = compute_cdf_and_bins(gdf_pop_1990_exposed_CF, bins_coarse)

bin_centers = bins_fine[:-1] + np.diff(bins_fine) / 2
bin_centers_coarse = bins_coarse[:-1] + np.diff(bins_coarse) / 2

print("✔ Analysis pipeline complete. Ready for plotting.")

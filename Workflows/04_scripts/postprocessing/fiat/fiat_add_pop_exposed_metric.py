"""
Add population-exposed metric to Delft-FIAT spatial output.

For each scenario listed in EVENT_CONFIG, allocates WorldPop grid-cell population
to flood-affected buildings proportionally by building footprint area.  Writes:
  - spatial_with_pop_and_flood.fgb  (spatial.fgb + 'population' column)
  - exposure_with_pop_and_flood.csv (exposure.csv + 'population' column)

by @dumontgoulart
"""

import rasterio
import geopandas as gpd
import numpy as np
from shapely.geometry import box
import pandas as pd
from pathlib import Path
import rasterio.windows
from rasterio.mask import mask as rasterio_mask
import warnings

warnings.filterwarnings("ignore")

# ===== CONFIGURATION =====
EVENT_NAME = "Durban2022"  # Change to: "Idai", "Durban2022"

EVENT_CONFIG = {
    "Idai": {
        "base_path": Path("/p/11210471-001-compass/03_Runs/sofala"),
        "scenarios": [
            "event_tp_era5_hourly_zarr_CF0_GTSMv41_CF0_era5_hourly_spw_IBTrACS_CF0",
            "event_tp_era5_hourly_zarr_CF-8_GTSMv41_CF-0.14_era5_hourly_spw_IBTrACS_CF-10",
        ],
        "population_raster": "/p/11210471-001-compass/01_Data/population_data/Worldpop/moz_ppp_2020_UNadj_constrained.tif",
        "utm_crs": "EPSG:32736",
    },
    "Durban2022": {
        "base_path": Path("/p/11210471-001-compass/03_Runs/durban"),
        "scenarios": [
            "event_precip_era5_hourly_CF0_no_wind",
            "event_precip_era5_hourly_CF-8_no_wind",
            "event_precip_era5_hourly_CF-10_no_wind",
        ],
        "population_raster": "/p/11210471-001-compass/01_Data/population_data/Worldpop/zaf_pop_2021_CN_100m_R2025A_v1.tif",
        "utm_crs": "EPSG:32736",
    },
}


# ===== CORE PROCESSING FUNCTION =====
def process_scenario(base_dir: Path, population_raster_path: str, utm_crs: str):
    """
    Allocate WorldPop population to flood-affected buildings for one scenario.

    Parameters
    ----------
    base_dir : Path
        Root of the FIAT scenario directory (contains output/, exposure/, hazard/).
    population_raster_path : str
        Path to the WorldPop GeoTIFF (EPSG:4326, population per pixel).
    utm_crs : str
        EPSG code for a projected CRS used to calculate building areas.
    """
    buildings_fp_path = base_dir / "output" / "spatial.fgb"
    exposure_csv_path = base_dir / "exposure" / "exposure.csv"
    output_dir = base_dir / "output"
    output_dir.mkdir(exist_ok=True)
    output_buildings_path = output_dir / "spatial_with_pop_and_flood.fgb"
    output_exposure_path = output_dir / "exposure_with_pop_and_flood.csv"

    print(f"\n  Base dir:   {base_dir}")
    print(f"  Pop raster: {population_raster_path}")

    # ----- Load buildings -----
    print("  Loading spatial.fgb...")
    if not buildings_fp_path.exists():
        print(f"  ERROR: {buildings_fp_path} not found — skipping.")
        return

    buildings_fp = gpd.read_file(buildings_fp_path)
    if "total_damage" in buildings_fp.columns:
        buildings_fp = buildings_fp[buildings_fp["total_damage"] > 0].copy()
    print(f"  Affected buildings: {len(buildings_fp):,}")

    # ----- Load exposure CSV -----
    print("  Loading exposure.csv...")
    exposure_df = pd.read_csv(exposure_csv_path) if exposure_csv_path.exists() else None

    # ----- Building areas -----
    if "Shape_Area" not in buildings_fp.columns:
        proj = buildings_fp.to_crs(utm_crs)
        buildings_fp["Shape_Area"] = proj.geometry.area
    else:
        print("  Using existing Shape_Area column.")

    # ----- Centroids in EPSG:4326 -----
    if str(buildings_fp.crs) != "EPSG:4326":
        buildings_fp = buildings_fp.to_crs("EPSG:4326")
    buildings_fp["centroid"] = buildings_fp.geometry.centroid
    result_df = buildings_fp.copy()
    result_df = result_df.set_geometry("centroid")
    result_df["population"] = 0.0

    # ----- Clip population raster to building extent -----
    if not Path(population_raster_path).exists():
        print(f"  ERROR: population raster not found at {population_raster_path}")
        return

    bounds = result_df.total_bounds
    buf = 0.01
    clip_geom = [
        box(bounds[0] - buf, bounds[1] - buf, bounds[2] + buf, bounds[3] + buf)
    ]

    print("  Clipping population raster...")
    with rasterio.open(population_raster_path) as src:
        print(f"    Raster CRS: {src.crs}  |  shape: {src.shape}")
        try:
            clipped_data, clipped_transform = rasterio_mask(
                src, clip_geom, crop=True, nodata=src.nodata
            )
            clipped_data = clipped_data[0]
        except Exception as e:
            print(f"    Warning: could not clip raster ({e}). Using full raster.")
            clipped_data = src.read(1)
            clipped_transform = src.transform

    nodata = -99999
    valid_mask = (clipped_data > 0) & np.isfinite(clipped_data)
    valid_cells = clipped_data[valid_mask]
    print(
        f"    Population in clipped area: {valid_cells.sum():,.0f}  |  cells: {len(valid_cells):,}"
    )

    # ----- Vectorised population allocation -----
    print("  Allocating population to buildings (vectorised)...")
    building_sindex = result_df.sindex

    pop_cells = np.where(valid_mask)
    cell_rows, cell_cols = pop_cells
    cell_populations = clipped_data[pop_cells]

    batch_size = 1000
    n_batches = len(cell_rows) // batch_size + (1 if len(cell_rows) % batch_size else 0)

    for b in range(n_batches):
        s, e = b * batch_size, min((b + 1) * batch_size, len(cell_rows))
        if b % 10 == 0:
            print(f"    Batch {b + 1}/{n_batches}...")

        polygons = []
        for row, col in zip(cell_rows[s:e], cell_cols[s:e]):
            w = rasterio.windows.Window(col, row, 1, 1)
            polygons.append(box(*rasterio.windows.bounds(w, clipped_transform)))

        cells_gdf = gpd.GeoDataFrame(
            {
                "cell_id": range(s, e),
                "population": cell_populations[s:e],
                "geometry": polygons,
            },
            crs="EPSG:4326",
        )
        joined = gpd.sjoin(result_df, cells_gdf, how="inner", predicate="intersects")
        if joined.empty:
            continue

        for cell_id, grp in joined.groupby("cell_id"):
            cell_pop = cell_populations[cell_id - s]
            total_area = grp["Shape_Area"].sum()
            if total_area > 0:
                result_df.loc[grp.index, "population"] += (
                    grp["Shape_Area"] / total_area * cell_pop
                ).values

    print(f"  Total population allocated: {result_df['population'].sum():,.0f}")

    # ----- Restore original geometry for saving -----
    result_df = result_df.set_geometry(buildings_fp.geometry.name)
    result_df = result_df.drop(columns=["centroid"], errors="ignore")

    # ----- Save buildings output -----
    result_df.to_file(output_buildings_path, driver="FlatGeobuf")
    print(f"  Saved: {output_buildings_path}")

    # ----- Merge with exposure CSV and save -----
    if exposure_df is not None:
        bkeys = [
            c
            for c in result_df.columns
            if "object" in c.lower() or c.lower() in ["id", "fid"]
        ]
        ekeys = [
            c
            for c in exposure_df.columns
            if "object" in c.lower() or c.lower() in ["id", "fid"]
        ]
        common_key = next(
            ((b, e) for b in bkeys for e in ekeys if b.lower() == e.lower()), None
        )
        if common_key:
            exp_out = exposure_df.merge(
                result_df[[common_key[0], "population"]],
                left_on=common_key[1],
                right_on=common_key[0],
                how="left",
            )
            exp_out["population"] = exp_out["population"].fillna(0)
        else:
            exp_out = exposure_df.copy()
            exp_out["population"] = 0
        exp_out.to_csv(output_exposure_path, index=False)
        print(f"  Saved: {output_exposure_path}")

    # ----- Summary -----
    print(f"  Buildings with population > 0: {(result_df['population'] > 0).sum():,}")
    if "total_damage" in result_df.columns:
        print(
            f"  Total damage:                  ${result_df['total_damage'].sum():,.0f}"
        )
    print(f"  Total population exposed:      {result_df['population'].sum():,.0f}")


# ===== MAIN =====
def main():
    cfg = EVENT_CONFIG[EVENT_NAME]
    print(f"\n{'='*60}")
    print(f"Population allocation — {EVENT_NAME}")
    print(f"{'='*60}")

    for scenario in cfg["scenarios"]:
        print(f"\nScenario: {scenario}")
        base_dir = cfg["base_path"] / EVENT_NAME / "fiat" / scenario
        process_scenario(base_dir, cfg["population_raster"], cfg["utm_crs"])

    print(f"\n{'='*60}")
    print(f"Done — {EVENT_NAME}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()

# Postprocessing examples

Tutorial notebooks for exploring SFINCS and FIAT outputs. No simulations are run — all notebooks read existing result files from completed pipeline runs.

**Environment**: `pixi run -e compass-v1 jupyter lab`  
**Reference event**: TC Idai / Sofala (`/p/11210471-001-compass/03_Runs/sofala/Idai/`)

Each notebook has a `# --- Paths ---` cell at the top — edit `REGION`, `EVENT`, and `SCENARIO` there to switch events.

---

## Notebooks

### `01_explore_flood_map.ipynb`
Load the maximum flood depth raster (`sfincs_output_hmax_AllTime.tif`) and compute basic flood metrics: extent (km²), volume (m³), and mean depth (m). Good starting point for any new run.

### `02_flood_timeseries.ipynb`
Open the raw gridded NetCDF (`sfincs_map.nc`) and reconstruct how the flood evolved day by day. Covers time decoding with `cftime`, exact cell-area calculation from corner coordinates, and plotting the flood hydrograph.

### `03_climate_attribution.ipynb`
Compare the factual scenario against one or more counterfactuals (e.g. −8% precipitation, −0.1 m sea level) to quantify the influence of climate change. Produces side-by-side flood maps, a diverging difference map, and an attribution bar chart.

### `04_damage_analysis.ipynb`
Load FIAT building-level damage output (`output/spatial.fgb`), explore the damage distribution, map losses spatially, and compare total damage across two scenarios.

---

## Scenario folder naming

```
event_tp_{precip}_{CF_rain}_{tide}_{CF_SLR}_{wind}_{CF_wind}
```

`CF0` = factual (no adjustment). `CF-8` on the precipitation term = −8% precipitation counterfactual.

## Related scripts

The production attribution scripts in `../attribution/` automate these analyses across all events and save summary CSVs to `/p/11210471-001-compass/04_Results/`.

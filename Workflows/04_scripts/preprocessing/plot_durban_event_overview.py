"""Create an overview figure for the Durban April 2022 precipitation event.

The figure contains:
1. A basemap with the analysed SFINCS domain highlighted.
2. A precipitation panel showing the domain-mean forcing rate and cumulative
   precipitation over the SFINCS domain.

By default, the script reads the Durban precipitation-only Snakemake config and
uses the first configured rainfall counterfactual value, which is the factual
run for the current Durban setup.
"""

from __future__ import annotations

import argparse
import os
import warnings
from pathlib import Path

import contextily as ctx
import geopandas as gpd
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import yaml
from shapely.geometry import box

warnings.filterwarnings("ignore")


DEFAULT_CONFIG = (
    Path(__file__).resolve().parents[2]
    / "01_config_snakemake"
    / "config_durban_floods_2022.yml"
)
DEFAULT_PRECIP_FORCING = "era5_hourly"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a Durban event overview figure from SFINCS inputs."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Path to the Snakemake config YAML.",
    )
    parser.add_argument(
        "--runname",
        default="Durban2022",
        help="Runname key inside the Snakemake config.",
    )
    parser.add_argument(
        "--precip-forcing",
        default=DEFAULT_PRECIP_FORCING,
        help="Precipitation forcing folder suffix to plot, for example era5_hourly or mswep_v316_3h_durban_apr2022.",
    )
    parser.add_argument(
        "--cf-rain",
        type=str,
        default=None,
        help="Optional rainfall counterfactual value overriding the first configured value.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output PNG path. Defaults to /p/<root>/04_Results/event_overview/.",
    )
    return parser.parse_args()


def load_config(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def get_root_path(config: dict) -> Path:
    if os.name == "nt":
        return Path("p:/") / config["root_dir"]
    return Path("/p") / config["root_dir"]


def format_cf_value(value: object) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        if float(value).is_integer():
            return str(int(value))
        return f"{value:g}"
    return str(value)


def get_run_settings(
    config: dict,
    runname: str,
    precip_forcing_override: str | None,
    cf_rain_override: str | None,
) -> dict:
    run_settings = config["runname_ids"][runname]
    cf_values = run_settings.get("CF_value_rain", [0])
    cf_rain = (
        cf_rain_override
        if cf_rain_override is not None
        else format_cf_value(cf_values[0])
    )
    return {
        "runname": runname,
        "event_name": run_settings["tc_name"],
        "region": run_settings["region"],
        "precip_forcing": precip_forcing_override or run_settings["precip_forcing"],
        "wind_forcing": run_settings["wind_forcing"],
        "cf_rain": cf_rain,
        "start_time": run_settings["start_time"],
        "end_time": run_settings["end_time"],
    }


def get_paths(config: dict, settings: dict, output_override: Path | None) -> dict:
    root_path = get_root_path(config)
    model_root = (
        root_path
        / config["dir_models"]
        / settings["region"]
        / settings["runname"]
        / "sfincs"
    )
    run_root = (
        root_path
        / config["dir_runs"]
        / settings["region"]
        / settings["runname"]
        / "sfincs"
        / f"event_precip_{settings['precip_forcing']}_CF{settings['cf_rain']}_{settings['wind_forcing']}"
    )
    output_path = output_override
    if output_path is None:
        output_path = (
            root_path
            / "04_Results"
            / "event_overview"
            / (
                f"durban_event_overview_{settings['runname'].lower()}_"
                f"{settings['precip_forcing'].lower()}.png"
            )
        )

    return {
        "region_geojson": model_root / "gis" / "region.geojson",
        "precip_netcdf": run_root / "precip_2d.nc",
        "sfincs_inp": run_root / "sfincs.inp",
        "output_png": output_path,
    }


def validate_inputs(paths: dict) -> None:
    missing = [
        str(path)
        for path in paths.values()
        if path.suffix != ".png" and not path.exists()
    ]
    if missing:
        raise FileNotFoundError(
            "Missing required input files:\n- " + "\n- ".join(missing)
        )


def load_region(region_path: Path) -> gpd.GeoDataFrame:
    region = gpd.read_file(region_path)
    if region.empty:
        raise ValueError(f"Region file is empty: {region_path}")
    if region.crs is None:
        region = region.set_crs("EPSG:4326")
    return region


def load_precipitation(precip_path: Path) -> tuple[xr.DataArray, str]:
    dataset = xr.open_dataset(precip_path)
    precip_vars = [name for name in dataset.data_vars if "precip" in name.lower()]
    if not precip_vars:
        raise ValueError(
            f"No precipitation variable found in {precip_path}. Available variables: {list(dataset.data_vars)}"
        )
    return dataset[precip_vars[0]], precip_vars[0]


def load_sfincs_input(sfincs_inp_path: Path) -> dict:
    sfincs_config: dict[str, str] = {}
    with sfincs_inp_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = [item.strip() for item in line.split("=", 1)]
            sfincs_config[key] = value
    return sfincs_config


def get_time_step_hours(precipitation: xr.DataArray) -> float:
    if precipitation.sizes.get("time", 0) < 2:
        return 1.0

    deltas = (
        np.diff(precipitation.time.values).astype("timedelta64[m]").astype(float) / 60.0
    )
    dt_hours = float(np.median(deltas))
    if dt_hours <= 0 or dt_hours > 24:
        return 1.0
    return dt_hours


def get_spatial_dims(precipitation: xr.DataArray) -> list[str]:
    return [dim for dim in precipitation.dims if dim != "time"]


def summarise_precipitation(
    precipitation: xr.DataArray,
) -> tuple[np.ndarray, np.ndarray, float]:
    spatial_dims = get_spatial_dims(precipitation)
    if not spatial_dims:
        raise ValueError(
            "Expected precipitation data with time and spatial dimensions."
        )

    mean_rate = precipitation.mean(dim=spatial_dims, skipna=True)
    dt_hours = get_time_step_hours(precipitation)
    cumulative_depth = np.cumsum(mean_rate.values * dt_hours)
    return mean_rate.values, cumulative_depth, dt_hours


def build_map_extent(region: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    region_geo = region.to_crs("EPSG:4326")
    minx, miny, maxx, maxy = region_geo.total_bounds
    dx = maxx - minx
    dy = maxy - miny
    buffer_x = max(dx * 1.5, 0.08)
    buffer_y = max(dy * 1.5, 0.08)
    extent_geom = box(
        minx - buffer_x, miny - buffer_y, maxx + buffer_x, maxy + buffer_y
    )
    return gpd.GeoDataFrame(geometry=[extent_geom], crs="EPSG:4326")


def get_projected_region(region: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    if region.crs is not None and not region.crs.is_geographic:
        return region
    projected_crs = region.estimate_utm_crs() or "EPSG:32736"
    return region.to_crs(projected_crs)


def get_resolution_from_coords(values: np.ndarray) -> float | None:
    if values.size < 2:
        return None
    deltas = np.diff(values.astype(float))
    if deltas.size == 0:
        return None
    return float(np.abs(np.median(deltas)))


def get_forcing_grid_resolution_m(
    precipitation: xr.DataArray,
) -> tuple[float | None, float | None]:
    spatial_dims = get_spatial_dims(precipitation)
    if len(spatial_dims) < 2:
        return None, None

    first_dim, second_dim = spatial_dims[:2]
    first_coords = precipitation.coords.get(first_dim)
    second_coords = precipitation.coords.get(second_dim)
    if first_coords is None or second_coords is None:
        return None, None

    dx = get_resolution_from_coords(np.asarray(first_coords.values))
    dy = get_resolution_from_coords(np.asarray(second_coords.values))
    if dx is None or dy is None:
        return None, None

    first_name = first_dim.lower()
    second_name = second_dim.lower()
    if first_name in {"lon", "longitude", "x"} and second_name in {
        "lat",
        "latitude",
        "y",
    }:
        if "lon" in first_name or "lat" in second_name:
            lat_mean = float(np.nanmean(np.asarray(second_coords.values)))
            dx = dx * 111320.0 * np.cos(np.deg2rad(lat_mean))
            dy = dy * 111320.0
    return dx, dy


def collect_metadata(
    region: gpd.GeoDataFrame,
    precipitation: xr.DataArray,
    precip_var_name: str,
    precip_path: Path,
    sfincs_input: dict,
) -> dict:
    mean_rate, cumulative_depth, dt_hours = summarise_precipitation(precipitation)
    timestamps = precipitation.time.values
    peak_index = int(np.nanargmax(mean_rate))
    projected_region = get_projected_region(region)
    area_km2 = float(projected_region.area.sum() / 1e6)
    bounds = projected_region.total_bounds
    forcing_dx_m, forcing_dy_m = get_forcing_grid_resolution_m(precipitation)

    dx_m = float(sfincs_input.get("dx", "nan"))
    dy_m = float(sfincs_input.get("dy", "nan"))
    grid_cols = int(float(sfincs_input.get("mmax", "0")))
    grid_rows = int(float(sfincs_input.get("nmax", "0")))
    active_cells = (
        int(round((area_km2 * 1e6) / (dx_m * dy_m))) if dx_m > 0 and dy_m > 0 else None
    )

    return {
        "forcing_file": precip_path.name,
        "forcing_variable": precip_var_name,
        "forcing_time_start": np.datetime_as_string(timestamps[0], unit="m"),
        "forcing_time_end": np.datetime_as_string(timestamps[-1], unit="m"),
        "forcing_time_step_hours": dt_hours,
        "forcing_time_steps": int(precipitation.sizes.get("time", 0)),
        "forcing_grid_shape": " x ".join(
            str(precipitation.sizes[dim]) for dim in get_spatial_dims(precipitation)
        ),
        "forcing_resolution_m": (forcing_dx_m, forcing_dy_m),
        "peak_mean_rate_mm_per_hr": float(np.nanmax(mean_rate)),
        "peak_time": np.datetime_as_string(timestamps[peak_index], unit="m"),
        "total_cumulative_mm": float(cumulative_depth[-1]),
        "domain_area_km2": area_km2,
        "domain_width_km": float((bounds[2] - bounds[0]) / 1e3),
        "domain_height_km": float((bounds[3] - bounds[1]) / 1e3),
        "model_resolution_m": (dx_m, dy_m),
        "model_grid_shape": f"{grid_cols} x {grid_rows}",
        "model_crs": sfincs_input.get("epsg", "unknown"),
        "model_utm_zone": sfincs_input.get("utmzone", "unknown"),
        "active_cells_approx": active_cells,
    }


def create_figure(
    region: gpd.GeoDataFrame,
    precipitation: xr.DataArray,
    metadata: dict,
    output_path: Path,
    event_name: str,
    runname: str,
    precip_label: str,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    mean_rate, cumulative_depth, dt_hours = summarise_precipitation(precipitation)
    timestamps = precipitation.time.values

    region_extent = build_map_extent(region)
    region_web = region.to_crs(epsg=3857)
    extent_web = region_extent.to_crs(epsg=3857)

    fig = plt.figure(figsize=(14, 7))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.1, 1.0], wspace=0.2)
    ax_map = fig.add_subplot(gs[0, 0])
    ax_ts = fig.add_subplot(gs[0, 1])

    extent_web.boundary.plot(ax=ax_map, alpha=0)
    try:
        ctx.add_basemap(
            ax_map,
            source=ctx.providers.Esri.WorldImagery,
            crs=extent_web.crs,
            attribution=False,
        )
    except Exception as exc:
        print(f"Could not add basemap: {exc}")

    region_web.plot(
        ax=ax_map,
        facecolor="#ef476f",
        edgecolor="#8d1428",
        linewidth=2.0,
        alpha=0.35,
        zorder=3,
    )
    region_web.boundary.plot(ax=ax_map, color="#8d1428", linewidth=2.0, zorder=4)

    minx, miny, maxx, maxy = extent_web.total_bounds
    ax_map.set_xlim(minx, maxx)
    ax_map.set_ylim(miny, maxy)
    ax_map.set_title("Analysed SFINCS domain", fontsize=13, fontweight="bold")
    ax_map.set_axis_off()

    bar_width_days = max(dt_hours / 24.0 * 0.85, 1 / 24.0)
    ax_ts.bar(
        timestamps,
        mean_rate,
        width=bar_width_days,
        color="#3a86ff",
        alpha=0.55,
        label="Domain-mean precipitation rate",
    )
    ax_ts.set_ylabel("Mean precipitation rate [mm/hr]", color="#1d4e89")
    ax_ts.tick_params(axis="y", labelcolor="#1d4e89")
    ax_ts.grid(True, alpha=0.25)

    ax_cumulative = ax_ts.twinx()
    ax_cumulative.plot(
        timestamps,
        cumulative_depth,
        color="#d1495b",
        linewidth=2.4,
        label="Cumulative precipitation",
    )
    ax_cumulative.set_ylabel("Cumulative precipitation [mm]", color="#8d1428")
    ax_cumulative.tick_params(axis="y", labelcolor="#8d1428")

    ax_ts.set_title(
        "SFINCS-domain precipitation forcing", fontsize=13, fontweight="bold"
    )
    ax_ts.set_xlabel("Time")
    ax_ts.xaxis.set_major_formatter(mdates.DateFormatter("%d %b\n%H:%M"))
    ax_ts.tick_params(axis="x", rotation=0)

    peak_index = int(np.nanargmax(mean_rate))
    ax_ts.axvline(
        timestamps[peak_index],
        color="#1d3557",
        linestyle="--",
        linewidth=1.2,
        alpha=0.7,
    )

    left_handles, left_labels = ax_ts.get_legend_handles_labels()
    right_handles, right_labels = ax_cumulative.get_legend_handles_labels()
    ax_ts.legend(
        left_handles + right_handles,
        left_labels + right_labels,
        loc="upper left",
        bbox_to_anchor=(0.02, 0.98),
        frameon=False,
        borderaxespad=0.0,
        ncol=1,
    )

    fig.suptitle(
        f"{event_name} overview ({runname})",
        fontsize=16,
        fontweight="bold",
        y=0.98,
    )
    fig.tight_layout(rect=(0.0, 0.03, 1.0, 0.95))
    fig.savefig(str(output_path), dpi=300, bbox_inches="tight")
    plt.close(fig)


def print_metadata_summary(metadata: dict, output_path: Path) -> None:
    forcing_dx_m, forcing_dy_m = metadata["forcing_resolution_m"]
    forcing_resolution_text = "n/a"
    if forcing_dx_m is not None and forcing_dy_m is not None:
        forcing_resolution_text = f"{forcing_dx_m:.1f} x {forcing_dy_m:.1f} m"

    print("Overview metadata")
    print(f"  Output figure: {output_path}")
    print(f"  Forcing file: {metadata['forcing_file']}")
    print(f"  Forcing variable: {metadata['forcing_variable']}")
    print(
        "  Forcing time window: "
        f"{metadata['forcing_time_start']} to {metadata['forcing_time_end']}"
    )
    print(
        "  Forcing timestep: "
        f"{metadata['forcing_time_step_hours']:.1f} h across {metadata['forcing_time_steps']} steps"
    )
    print(
        f"  Forcing grid: {metadata['forcing_grid_shape']} at {forcing_resolution_text}"
    )
    print(
        f"  Peak domain-mean rain rate: {metadata['peak_mean_rate_mm_per_hr']:.2f} mm/hr"
    )
    print(f"  Peak time: {metadata['peak_time']}")
    print(f"  Total cumulative precipitation: {metadata['total_cumulative_mm']:.2f} mm")
    print(f"  SFINCS domain area: {metadata['domain_area_km2']:.2f} km2")
    print(
        f"  SFINCS domain extent: {metadata['domain_width_km']:.2f} x {metadata['domain_height_km']:.2f} km"
    )
    print(
        "  SFINCS model resolution: "
        f"{metadata['model_resolution_m'][0]:.1f} x {metadata['model_resolution_m'][1]:.1f} m"
    )
    print(f"  SFINCS grid size: {metadata['model_grid_shape']}")
    print(f"  Approximate active cells: {metadata['active_cells_approx']:,}")
    print(f"  SFINCS CRS: EPSG:{metadata['model_crs']} ({metadata['model_utm_zone']})")


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    settings = get_run_settings(
        config,
        args.runname,
        args.precip_forcing,
        args.cf_rain,
    )
    paths = get_paths(config, settings, args.output)
    validate_inputs(paths)

    region = load_region(paths["region_geojson"])
    precipitation, precip_var_name = load_precipitation(paths["precip_netcdf"])
    sfincs_input = load_sfincs_input(paths["sfincs_inp"])
    metadata = collect_metadata(
        region=region,
        precipitation=precipitation,
        precip_var_name=precip_var_name,
        precip_path=paths["precip_netcdf"],
        sfincs_input=sfincs_input,
    )
    create_figure(
        region=region,
        precipitation=precipitation,
        metadata=metadata,
        output_path=paths["output_png"],
        event_name=settings["event_name"],
        runname=settings["runname"],
        precip_label=settings["precip_forcing"],
    )
    print_metadata_summary(metadata, paths["output_png"])
    print(f"Saved overview figure to: {paths['output_png']}")


if __name__ == "__main__":
    main()

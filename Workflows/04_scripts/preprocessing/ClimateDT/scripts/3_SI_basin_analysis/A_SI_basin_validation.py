#%%
import sys
from pathlib import Path

# Handle both script and interactive environments
if '__file__' in globals():
    sys.path.insert(0, str(Path(__file__).parent.parent))
else:
    # In interactive environment, add the scripts directory manually
    sys.path.insert(0, r"C:\Code\COMPASS\Workflows\04_scripts\preprocessing\ClimateDT\scripts")

import os
import pickle
import platform
from networkx import display
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from track_comparison_funcs import align_tracks, normalize_pressure_column
from general_TC_funcs import import_ibtracs_period
from track_comparison_funcs import add_aligned_wind, create_aligned_tracks_dataframe, load_existing_tracks, compare_climatedt_tracks, plot_tracks_all
import geopandas as gdf

processed_path = Path("C:/Code/COMPASS/Workflows/04_scripts/preprocessing/ClimateDT/data/processed/")
track_file = (processed_path / "storms_all_strategies_2017-01-01_2026-08-01.pkl")
wind_files = {"200km": processed_path / "storms_all_ibtracs_max_wind_200km_2017-01-01_2026-08-01.pkl",
              "rmw": processed_path / "storms_all_ibtracs_max_wind_rmw_2017-01-01_2026-08-01.pkl"}
precip_file = processed_path / "precip"/ "landfall_precip.nc"

start_time = "2017-01-01"
end_time = "2026-08-01"

storms = import_ibtracs_period(start=start_time, end=end_time)
tracks = load_existing_tracks(track_file)

# Keep all aligned points available for inspection before comparison.
aligned_tracks = create_aligned_tracks_dataframe(storms, tracks)

for wind_key in wind_files:
    aligned_tracks = add_aligned_wind(aligned_tracks, load_existing_tracks(wind_files[wind_key]), wind_key)


# %%
########################################################################
################## INSPECTING EXCLUDED IBTRACS STORMS ##################
########################################################################
def plot_storms_by_year(storms):
    fig, ax = plt.subplots(figsize=(12, 6))

    year_valid = {}
    year_invalid = {}

    for storm_name, df in storms.items():
        # Normalize pressure
        ref_pres_col = 'min_pres'
        df = normalize_pressure_column(df, ref_pres_col)
        
        # Extract year from storm's first timestamp
        year = pd.to_datetime(df['time']).dt.year.iloc[0]
        
        # Check if storm has ANY valid data points
        has_valid = (df['lat'].notna() & df['lon'].notna() & df['msl'].notna()).any()
        
        if year not in year_valid:
            year_valid[year] = 0
            year_invalid[year] = 0
        
        if has_valid:
            year_valid[year] += 1
        else:
            year_invalid[year] += 1

    # Sort by year and prepare data
    years = sorted(year_valid.keys())
    valid_counts = [year_valid[y] for y in years]
    invalid_counts = [year_invalid[y] for y in years]

    # Create stacked bar chart
    x_pos = range(len(years))
    ax.bar(x_pos, valid_counts, label=f'Storms with valid data (n={sum(valid_counts)})', alpha=0.8)
    ax.bar(x_pos, invalid_counts, bottom=valid_counts, label=f'Storms without valid data (n={sum(invalid_counts)})', alpha=0.5)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(years)
    ax.set_xlabel('Year')
    ax.set_ylabel('Number of Storms')
    ax.set_title('Storm Count by Year (colored by data validity)')
    ax.legend()
    plt.tight_layout()
    plt.show()

def filter_storms_without_valid_data(storms):
    storms_valid = {}

    for storm_name, df in storms.items():
        # Normalize pressure
        ref_pres_col = 'min_pres'
        df = normalize_pressure_column(df, ref_pres_col)
        
        # Check if storm has ANY valid data points
        has_valid = (df['lat'].notna() & df['lon'].notna() & df['msl'].notna()).any()
        
        if has_valid:
            storms_valid[storm_name] = df
        else:
            storms_without_valid[storm_name] = df

    # Display summary
    print(f"Total storms: {len(storms)}")
    print(f"Storms with valid data: {len(storms_valid)}")
    print("\nStorms with valid data:")
    for storm_name in sorted(storms_valid.keys()):
        df = storms_valid[storm_name]
        year = pd.to_datetime(df['time']).dt.year.iloc[0]
        print(f"  {storm_name} ({year}): {len(df)} rows")

    return storms_valid, storms_without_valid

def inspect_storms_without_valid_data(storms_without_valid):
    # Plot distribution of storms without valid data by year
    # Confirm: storms with PARTIAL valid data are still included in the valid group
    # (only storms with ZERO valid rows are in storms_without_valid)

    # Plot usa_sshs distribution for storms WITHOUT valid data
    fig, ax = plt.subplots(figsize=(10, 6))

    all_sshs = []
    for storm_name, df in storms_without_valid.items():
        # Collect all usa_sshs values (including NaN)
        all_sshs.extend(df['usa_sshs'].dropna().values)

    # Plot distribution
    ax.hist(all_sshs, bins=range(int(min(all_sshs))-1, int(max(all_sshs))+2), 
            edgecolor='black', alpha=0.7)

    ax.set_xlabel('USA SSHS Category')
    ax.set_ylabel('Count')
    ax.set_title(f'USA SSHS Distribution - Storms WITHOUT Valid Data (n={len(all_sshs)} points)')
    ax.set_xticks(range(int(min(all_sshs)), int(max(all_sshs))+1))
    plt.tight_layout()
    plt.show()

    # Also show breakdown by storm
    ssh_storms_not_valid = pd.DataFrame()
    
    for storm_name in sorted(storms_without_valid.keys()):
        df = storms_without_valid[storm_name]
        sshs_counts = df['usa_sshs'].value_counts().sort_index()
        ssh_storms_not_valid = pd.concat([ssh_storms_not_valid, pd.DataFrame({
            'storm': [storm_name] * len(sshs_counts),
            'sshs': sshs_counts.index,
            'count': sshs_counts.values
        })], ignore_index=True)
    return ssh_storms_not_valid

# applying the functions to inspect storms without valid data
plot_storms_by_year(storms)

storms_valid, storms_without_valid = filter_storms_without_valid_data(storms)

ssh_storms_not_valid = inspect_storms_without_valid_data(storms_without_valid)

#%%
# Display IBTrACS points south of 40S
south_of_40 = pd.concat([track.loc[track["lat"] < -40].assign(storm=name) 
                        for name, track in storms.items()], ignore_index=True)

display(south_of_40[["storm", "time", "lat", "lon", "usa_sshs"]])

# SSHS counts for each storm, including missing values
display(south_of_40.groupby("storm")["usa_sshs"].value_counts(dropna=False))

# %%
#######################################################################
################## COMPARING TRACKING STRATEGIES ######################
#######################################################################
# Calculate comparison metrics for track and msl, with and without filters
# Calculate with no filters at all
results_no_filt = compare_climatedt_tracks(storms, tracks, 
                                        #    tracking_strategy="all_ibtracs",
                                           scen='hist',
                                           min_ref_lat=None,
                                           min_ref_sshs=None,
                                           land_points=None,
                                           )

# Calculate metrics only for hist scenario with all filters applied to metrics
results_filt_all = compare_climatedt_tracks(storms, tracks, 
                                            # tracking_strategy="all_ibtracs",
                                            scen='hist',
                                            min_ref_lat=-40,
                                            min_ref_sshs=1,
                                            land_points='excl',
                                            )

# Calculate with no landfall exclusion but other filters
results_filt_40S = compare_climatedt_tracks(storms, tracks, 
                                            # tracking_strategy="all_ibtracs",
                                            scen='hist',
                                            min_ref_lat=-40,
                                            min_ref_sshs=None,
                                            land_points=None,
                                            )

results_filt_TC1 = compare_climatedt_tracks(storms, tracks, 
                                            # tracking_strategy="all_ibtracs",
                                            scen='hist',
                                            min_ref_lat=None,
                                            min_ref_sshs=1,
                                            land_points=None,
                                            )

results_filt_no_land = compare_climatedt_tracks(storms, tracks, 
                                                # tracking_strategy="all_ibtracs",
                                                scen='hist',
                                                min_ref_lat=None,
                                                min_ref_sshs=None,
                                                land_points='excl',
                                                )

#%%
# Comparing msl and track of the factual (hist) storms with all strategies
def plot_violin_metrics_strategies(results_no_filt):
    # Key metrics to plot (6 total)
    metrics = [
        ('MDE_km', 'Mean Position Error (km)', 15),
        ('MXE_km', 'Maximum Position Error (km)', 15),
        ('position_RMSE_km', 'Position RMSE (km)', 10),
        ('RMSE_pres_hPa', 'Pressure RMSE (hPa)', 15),
        ('bias_pres_hPa', 'Pressure Bias (hPa)', 15),
        ('corr_pres', 'Pressure Correlation', 20),
    ]
        
    # Color each violin by strategy and customize lines
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    axes = axes.flatten()
    strategies = sorted(results_no_filt['strategy'].unique())
    strategy_colors = {'2point': '#A1C9F4', '5point': '#FFB482', 'all_ibtracs': '#CC78BC'}

    # Selecting the factual scenario (called hist)
    data_F = results_no_filt[results_no_filt['scenario'] == 'hist']

    for idx, (metric, label, _) in enumerate(metrics):
        ax = axes[idx]
        data_to_plot = []
        
        for s in strategies:
            data = data_F[data_F['strategy'] == s][metric].dropna().values
            data_to_plot.append(data)
        
        # Create violins
        parts = ax.violinplot(data_to_plot, positions=range(len(strategies)), 
                            showmeans=True, showmedians=True)
        
        # Customize all parts
        for pc, strategy in zip(parts['bodies'], strategies):
            pc.set_facecolor(strategy_colors[strategy])
            pc.set_alpha(0.6)
            pc.set_edgecolor('black')
        
        parts['cmeans'].set_color('black')
        parts['cmeans'].set_linewidth(2)
        parts['cmeans'].set_linestyle('--')
        
        # Make whiskers/mins/maxes subtle
        parts['cbars'].set_color('gray')
        parts['cbars'].set_linewidth(0.8)
        parts['cmaxes'].set_color('gray')
        parts['cmaxes'].set_linewidth(0.8)
        parts['cmins'].set_color('gray')
        parts['cmins'].set_linewidth(0.8)
        
        # Now only the manual colored medians are drawn
        for pos, strategy in enumerate(strategies):
            data = data_to_plot[pos]
            median_val = np.median(data)
            mean_val = np.mean(data)
            
            ax.plot([pos - 0.1, pos + 0.1], [median_val, median_val], 
                color=strategy_colors[strategy], linewidth=3, linestyle='-')
            
            # Add mean label
            ax.text(pos, ax.get_ylim()[1] * 0.95, f'{mean_val:.2f}', 
                ha='center', fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor=strategy_colors[strategy], alpha=0.3))
            
        ax.set_xticks(range(len(strategies)))
        ax.set_xticklabels(strategies)
        ax.set_ylabel(label)
        ax.set_title(label, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        
        # Legend (only on first plot)
        if idx == 0:
            from matplotlib.lines import Line2D
            legend_elements = [Line2D([0], [0], color='black', linestyle='--', linewidth=2, label='Mean')]
            legend_elements += [Line2D([0], [0], color=strategy_colors[s], linewidth=2.5, label=f'{s} Median') 
                            for s in strategies]
            ax.legend(handles=legend_elements, fontsize=8, loc='upper right',
                    bbox_to_anchor=(0.98, 0.95))

    plt.suptitle('Tracking Metrics Comparison', fontsize=14, fontweight='bold')
    plt.show()

# Comparing msl and track mertrics for factual using all_ibtracs strategy with different filters applied
def plot_violin_metrics_filters(results_no_filt, results_filt_40S, results_filt_TC1, results_filt_no_land, results_filt_all):
    # Only use the all_ibtracs strategy for this plot
    # Use the pre-computed results dataframes instead of lambda filters
    filters = {'No filters': results_no_filt[results_no_filt['strategy'] == 'all_ibtracs'],
            'Lat > -40': results_filt_40S[results_filt_40S['strategy'] == 'all_ibtracs'],
            'SSHS >= 1': results_filt_TC1[results_filt_TC1['strategy'] == 'all_ibtracs'],
            'No landfall': results_filt_no_land[results_filt_no_land['strategy'] == 'all_ibtracs'],
            'All filters': results_filt_all[results_filt_all['strategy'] == 'all_ibtracs']}

    filter_colors = {'No filters': '#CC78BC',
                    'Lat > -40': '#fb8072',
                    'SSHS >= 1': '#80b1d3',
                    'No landfall': '#8dd3c7',
                    'All filters': '#ffffb3'}

    # Create violin plots for position RMSE and pressure correlation
    metrics_to_plot = [('position_RMSE_km', 'Position RMSE (km)'),
                    ('corr_pres', 'Pressure Correlation'),
                    ('RMSE_pres_hPa', 'Pressure RMSE (hPa)')]

    fig, axes = plt.subplots(3, 5, figsize=(18, 15), constrained_layout=True)

    # Calculate y-axis limits for each metric (row) across all filters
    ylim_per_metric = {}
    for metric_col, metric_label in metrics_to_plot:
        all_values = []
        for filtered_data in filters.values():
            data_to_plot = filtered_data[metric_col].dropna().values
            if len(data_to_plot) > 0:
                all_values.extend(data_to_plot)
        
        if all_values:
            ymin = np.percentile(all_values, 1)  # 1st percentile for better visualization
            ymax = np.percentile(all_values, 99)  # 99th percentile
            ylim_per_metric[metric_col] = (ymin, ymax)

    for metric_idx, (metric_col, metric_label) in enumerate(metrics_to_plot):
        for filter_idx, (filter_name, filtered_data) in enumerate(filters.items()):
            ax = axes[metric_idx, filter_idx]
            
            # Extract metric values from the already-filtered dataframe
            data_to_plot = filtered_data[metric_col].dropna().values
            
            # Create violin plot
            if len(data_to_plot) > 0:
                parts = ax.violinplot([data_to_plot], positions=[0], showmeans=True, showmedians=True)
                
                # Customize violin
                for pc in parts['bodies']:
                    pc.set_facecolor(filter_colors[filter_name])
                    pc.set_alpha(0.7)
                    pc.set_edgecolor('black')
                    pc.set_linewidth(1.5)
                
                # Customize mean and median
                parts['cmeans'].set_color('red')
                parts['cmeans'].set_linewidth(2)
                parts['cmeans'].set_linestyle('--')
                
                parts['cmedians'].set_color('darkblue')
                parts['cmedians'].set_linewidth(2)
                
                # Customize whiskers
                parts['cbars'].set_color('gray')
                parts['cbars'].set_linewidth(0.8)
                parts['cmaxes'].set_color('gray')
                parts['cmaxes'].set_linewidth(0.8)
                parts['cmins'].set_color('gray')
                parts['cmins'].set_linewidth(0.8)
                
                # Add statistics text
                mean_val = np.mean(data_to_plot)
                median_val = np.median(data_to_plot)
                ax.text(0.5, 0.95, f'mean={mean_val:.2f}\nmed={median_val:.2f}', 
                    transform=ax.transAxes, ha='center', va='top', fontsize=9,
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8))
                
                # Set y-axis limits per metric (same for all filters in the same row)
                if metric_col in ylim_per_metric:
                    ax.set_ylim(ylim_per_metric[metric_col])
                
                ax.set_xticks([])
                ax.grid(axis='y', alpha=0.3)
            else:
                ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center')
            
            # Labels
            if metric_idx == 0:
                ax.set_title(filter_name, fontweight='bold')
            if filter_idx == 0:
                ax.set_ylabel(metric_label, fontweight='bold')

    plt.suptitle('Tracking Metrics: Position RMSE and Pressure Correlation (all_ibtracs)', 
                fontsize=14, fontweight='bold')
    plt.show()

    fig.savefig("results/tracking/figures/tracking/metrics_violin_all_ibtracs.png", dpi=300)


plot_violin_metrics_strategies(results_no_filt)

plot_violin_metrics_filters(results_no_filt, results_filt_40S, results_filt_TC1, 
                            results_filt_no_land, results_filt_all)


#%%
# Plotting tracks for different strategies per storm
plot_tracks_all(storms, tracks, figure_path="../../figures/tracking/all_strategies")


# %%
#######################################################################
####################### VALIDATING WIND SPEED #########################
#######################################################################
wind_column = 'max_wind_model_200km'  # or 'max_wind_model_rmw'

df = aligned_tracks[aligned_tracks["strategy"] == "all_ibtracs"].copy()

# Compute per-storm metrics by scenario
metrics = []
for scenario in df['scenario'].unique():
    scenario_data = df[df['scenario'] == scenario]
    
    for storm_id in scenario_data['storm'].unique():
        storm_df = scenario_data[scenario_data['storm'] == storm_id]
        
        metrics.append({
            'scenario': scenario,
            'storm': storm_id,
            'wind_max': storm_df[wind_column].max(),
            'pres_min': storm_df['msl_tracked'].min(),
        })

metrics_df = pd.DataFrame(metrics)

# Plot PDF curves per scenario
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

scenarios = sorted(metrics_df['scenario'].unique())
colors = {'cont': 'blue', 'hist': 'orange', 'Tplus2.0K': 'red'}

# Max wind speed
for scenario in scenarios:
    data = metrics_df[metrics_df['scenario'] == scenario]['wind_max'].dropna()
    data.plot.kde(ax=axes[0], label=scenario, color=colors.get(scenario), linewidth=2)

axes[0].set_xlabel('Max Wind Speed (m/s)')
axes[0].set_ylabel('Probability Density')
axes[0].set_title('Max Wind Speed Distribution per Scenario')
axes[0].legend()
axes[0].grid(alpha=0.3)

# Min pressure
for scenario in scenarios:
    data = metrics_df[metrics_df['scenario'] == scenario]['pres_min'].dropna()
    data.plot.kde(ax=axes[1], label=scenario, color=colors.get(scenario), linewidth=2)

axes[1].set_xlabel('Min Pressure (hPa)')
axes[1].set_ylabel('Probability Density')
axes[1].set_title('Min Pressure Distribution per Scenario')
axes[1].legend()
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.show()

# Print confidence intervals (95%)
print("\nMax Wind Speed - 95% CI by Scenario:")
for scenario in scenarios:
    data = metrics_df[metrics_df['scenario'] == scenario]['wind_max'].dropna()
    ci = data.quantile([0.025, 0.975])
    print(f"{scenario}: {ci[0.025]:.1f} - {ci[0.975]:.1f} m/s (median: {data.median():.1f})")

print("\nMin Pressure - 95% CI by Scenario:")
for scenario in scenarios:
    data = metrics_df[metrics_df['scenario'] == scenario]['pres_min'].dropna()
    ci = data.quantile([0.025, 0.975])
    print(f"{scenario}: {ci[0.025]:.1f} - {ci[0.975]:.1f} hPa (median: {data.median():.1f})")




# %%
########################################################################
####################### VALIDATING PRECIPITATION #######################
########################################################################
def plot_landfall_precip(tracked_precip_file, scenario, realization, storm):
    import xarray as xr
    import pandas as pd
    import netCDF4
    import re

    """Plot all saved landfall candidates, with one panel per precipitation time."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.pyplot as plt

    def group_name(name):
        name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")
        if not name:
            raise ValueError("Storm name cannot be represented as a NetCDF group")
        return name

    parent = f"{scenario}/r{realization}/{group_name(storm)}"
    with netCDF4.Dataset(tracked_precip_file) as root:
        node = root
        for part in parent.split("/"):
            node = node.groups[part]
        events = sorted(node.groups)
    figures = []
    for event in events:
        with xr.open_dataset(tracked_precip_file, group=f"{parent}/{event}") as source:
            climdt_precip_data = source.load()
        count = climdt_precip_data.sizes["time"]
        columns = min(3, count)
        rows = (count + columns - 1) // columns
        projection = ccrs.PlateCarree()
        fig, axes = plt.subplots(rows, columns, figsize=(5 * columns, 4 * rows),
                                 subplot_kw={"projection": projection}, squeeze=False, layout="constrained")
        vmax = float(climdt_precip_data["tp_rate"].max(skipna=True))
        image = None
        for i, ax in enumerate(axes.flat):
            if i >= count:
                ax.set_visible(False)
                continue
            # background map
            ax.coastlines()
            ax.add_feature(cfeature.BORDERS, linewidth=0.5)
            ax.add_feature(cfeature.LAND, facecolor="lightgray")
            ax.add_feature(cfeature.OCEAN, facecolor="lightblue")

            image = ax.pcolormesh(climdt_precip_data.longitude, climdt_precip_data.latitude, climdt_precip_data.tp_rate.isel(time=i),
                                  transform=projection, cmap="Blues", vmin=0, vmax=vmax)
            ax.plot(climdt_precip_data.track_lon, climdt_precip_data.track_lat, color="red", marker=".",
                    label="24-hour track", transform=projection)
            ax.plot(climdt_precip_data.track_lon.isel(time=i), climdt_precip_data.track_lat.isel(time=i),
                    marker="o", color="orange", label="Current position", transform=projection)
            ax.plot(climdt_precip_data.attrs["model_landfall_longitude"], climdt_precip_data.attrs["model_landfall_latitude"],
                    marker="x", color="black", markersize=8, label="Model landfall position",
                    transform=projection)
            ax.set_title(pd.Timestamp(climdt_precip_data.time.values[i]).strftime("%Y-%m-%d %H:%M"))
        fig.colorbar(image, ax=list(axes.flat), label="Precipitation rate (mm/hour)", shrink=0.8)
        axes.flat[0].legend(loc="best")
        fig.suptitle(f"{storm} | {scenario} r{realization} | IBTrACS landfall {event}"
                     f"| USA_SSHS {climdt_precip_data.attrs['usa_sshs']:g}")

        output_path = Path(f"results/tracking/figures/precip/landfall_precip_{storm}_{scenario}_r{realization}.png")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path)
        plt.close(fig)  # Close figure to free memory
    return figures


def check_precip_data_availability(storms):
    # Filter each storm to complete observations only
    storms_valid = {}
    for storm_name, storm_df in storms.items():
        storms_valid[storm_name] = storm_df[storm_df[['lat', 'lon', 'min_pres']].notna().all(axis=1)]

    # Drop storms without data
    storms_valid = {name: df for name, df in storms_valid.items() if not df.empty}

    # Filter to landfall storms only (landfall == 0)
    storms_valid_landfall = {}
    for storm_name, storm_df in storms_valid.items():
        storms_valid_landfall[storm_name] = storm_df[storm_df['landfall'].eq(0)]

    # Then drop any that became empty after filtering
    storms_valid_landfall = {name: df for name, df in storms_valid_landfall.items() if not df.empty}
    
    # Print summary statistics
    print("\n" + "="*60)
    print("LANDFALL STORM DATA AVAILABILITY")
    print("="*60)
    print(f"Total storms: {len(storms.keys())}")
    print(f"Storms with valid tracking data: {len(storms_valid)}")
    print(f"Storms with landfall data: {len(storms_valid_landfall)}")
    print("="*60)
    return 

# Get statistics without plotting
stats = check_precip_data_availability(storms)

# Uncomment below to plot all storms with precipitation data
# for key, data in tracks.items():
#     strat, scen, real, storm = key

#     # Only plot storms with valid data
#     if data is not None and len(data) > 0:
#         print(storm)
#         try:
#             plot_landfall_precip(precip_file, scenario="hist", realization=1, storm=storm)
#         except (KeyError, FileNotFoundError) as e:
#             print(f"  → Skipped: {e}")


# %%

#%%
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx
import matplotlib.patches as mpatches
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
from hydromt_wflow import WflowModel


def plot_region_with_satellite(region, runname, discharge_locations_gdf, discharge_timeseries, save=True, custom_offsets = None):
    geojson_path = f'p:/11210471-001-compass/02_Models/{region}/{runname}/wflow/staticgeoms/basins.geojson'
    geojson_path_sfincs = f'p:/11210471-001-compass/02_Models/{region}/{runname}/sfincs/gis/region.geojson'
    geojson_path_rivers = f'p:/11210471-001-compass/02_Models/{region}/{runname}/wflow/staticgeoms/rivers.geojson'
    
    # Read the GeoJSON files
    gdf = gpd.read_file(geojson_path)
    gdf2 = gpd.read_file(geojson_path_sfincs)
    gdf2 = gdf2.to_crs(gdf.crs)  # Reproject to match gdf
    rivers = gpd.read_file(geojson_path_rivers).to_crs(gdf.crs) 

    # Plot the merged region
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal', 'box') 
    gdf.plot(ax=ax, edgecolor='black', facecolor='lightblue', alpha = 0.7)
    gdf2.plot(ax=ax, edgecolor='black', facecolor='red', alpha = 0.7)
    rivers.plot(ax=ax, edgecolor='blue', facecolor='none', alpha = 1, lw=0.4)

    # Adjust plot limits
    minx, miny, maxx, maxy = gdf.total_bounds
    buffer = 1  # Zoom-out buffer
    ax.set_xlim(minx - buffer, maxx + buffer)
    ax.set_ylim(miny - buffer, maxy + buffer)  
    ctx.add_basemap(ax, crs=gdf.crs.to_string(), source=ctx.providers.Esri.WorldImagery)

    # Add North Arrow
    north_arrow = mpatches.FancyArrowPatch(
        (0.9, 0.05), (0.9, 0.1), transform=ax.transAxes, color='white', arrowstyle='->', lw=2, mutation_scale=15
    )
    ax.add_patch(north_arrow)
    ax.text(0.91, 0.07, "N", transform=ax.transAxes, fontsize=14, color='white', va='bottom', ha='left', fontweight='bold')

    # Add the legend
    legend_elements = [
        mpatches.Patch(facecolor='lightblue', label='Wflow region', edgecolor='black'),
        mpatches.Patch(facecolor='red', label='SFINCS region', edgecolor='black')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=12, title="Regions", frameon=True, facecolor='white', edgecolor='black')

    # Add discharge time series graphs at given locations
    for idx, loc in discharge_locations_gdf.iterrows():
        # Extract the point geometry (longitude, latitude)
        point_geom = loc['geometry']
        point_x, point_y = point_geom.x, point_geom.y

        # Get the corresponding timeseries (from Q_gauges_locs)
        gauge_id = loc['index']  # assuming 'uparea' corresponds to the gauge ID
        timeseries = discharge_timeseries['Q'].sel(Q_gauges_locs=str(gauge_id))

        # Create inset axes for discharge time series plot
        # We'll use map coordinates (longitude, latitude) to position the inset
        inset_width = 1.9  # width of the inset in degrees
        inset_height = 1.9  # height of the inset in degrees

        # Position of the inset (in map coordinates)
        x_pos_inset = point_x + 0.02  # adjust position slightly to the right
        y_pos_inset = point_y + 0.02  # adjust position slightly up
        if custom_offsets and str(gauge_id) in custom_offsets:
            dx, dy = custom_offsets[str(gauge_id)]
       # Position of the inset (in map coordinates)
        x_pos_inset = point_x + dx  # adjust position based on custom offset
        y_pos_inset = point_y + dy  # adjust position based on custom offset

        # Ensure inset stays within map bounds
        # x_pos_inset = min(max(x_pos_inset, minx), maxx - inset_width)
        # y_pos_inset = min(max(y_pos_inset, miny), maxy - inset_height)

        # Create inset axes for timeseries plot
        ax_inset = inset_axes(ax, width=inset_width, height=inset_height, loc='lower left',
                              bbox_to_anchor=(x_pos_inset, y_pos_inset), bbox_transform=ax.transData)

        # Plot the discharge time series for the current location
        timeseries.plot(ax=ax_inset, color='green', lw=2)
        # Set all text in the inset to white
        ax_inset.set_title(f"Discharge", fontsize=10, color='white')
        ax_inset.set_xlabel('Time', fontsize=8, color='white')
        ax_inset.set_ylabel('Discharge (m³/s)', fontsize=8, color='white')
        ax_inset.tick_params(axis='both', which='major', labelsize=8, colors='white')


        # Add dashed lines from the point on the map to the inset plot
        ax.plot([point_x, x_pos_inset], [point_y, y_pos_inset], color='orange', linestyle='--', linewidth=2, alpha=0.7, transform=ax.transData)

        # Add the point on the map as a small marker
        ax.scatter(point_x, point_y, color='orange', s=100, edgecolor='black', zorder=5, transform=ax.transData)

    # Adjust plot appearance
    ax.set_title(f"Wflow results {region.capitalize()} for Hurricane {runname}", fontsize=20)
    ax.set_axis_off()  # Turn off the axis
    plt.tight_layout()
    
    if save:
        if not os.path.exists(f'p:/11210471-001-compass/04_Results/wflow'):
            os.mkdir(f'p:/11210471-001-compass/04_Results/wflow')
        fig.savefig(f'p:/11210471-001-compass/04_Results/wflow/{region}_{runname}_wflow_results.png')
    else:
        plt.show()


def plot_region_with_gauges_only(region, runname, discharge_locations_gdf, save=True, custom_offsets=None):
    geojson_path = f'p:/11210471-001-compass/02_Models/{region}/{runname}/wflow/staticgeoms/basins.geojson'
    geojson_path_sfincs = f'p:/11210471-001-compass/02_Models/{region}/{runname}/sfincs/gis/region.geojson'
    geojson_path_rivers = f'p:/11210471-001-compass/02_Models/{region}/{runname}/wflow/staticgeoms/rivers.geojson'
    
    # Read GeoJSONs
    gdf = gpd.read_file(geojson_path)
    gdf2 = gpd.read_file(geojson_path_sfincs).to_crs(gdf.crs)
    rivers = gpd.read_file(geojson_path_rivers).to_crs(gdf.crs)

    # Plot base map
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal', 'box')
    gdf.plot(ax=ax, edgecolor='black', facecolor='lightblue', alpha=0.7)
    gdf2.plot(ax=ax, edgecolor='black', facecolor='red', alpha=0.7)
    rivers.plot(ax=ax, edgecolor='blue', facecolor='none', alpha=1, lw=0.4)

    # Zoom to SFINCS region instead of Wflow
    minx, miny, maxx, maxy = gdf2.total_bounds
    buffer = 0.2  # Smaller buffer for a tighter zoom
    ax.set_xlim(minx - buffer, maxx + buffer)
    ax.set_ylim(miny - buffer, maxy + buffer)

    ctx.add_basemap(ax, crs=gdf.crs.to_string(), source=ctx.providers.Esri.WorldImagery, attribution=False)

    # Add north arrow
    north_arrow = mpatches.FancyArrowPatch((0.9, 0.05), (0.9, 0.1), transform=ax.transAxes, color='white', arrowstyle='->', lw=2, mutation_scale=15)
    ax.add_patch(north_arrow)
    ax.text(0.91, 0.07, "N", transform=ax.transAxes, fontsize=14, color='white', va='bottom', ha='left', fontweight='bold')

    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor='lightblue', label='Wflow region', edgecolor='black'),
        mpatches.Patch(facecolor='red', label='SFINCS region', edgecolor='black')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=12, title="Regions", frameon=True, facecolor='white', edgecolor='black')

    custom_offsets = {
    '1': (-0.02, 0.02),  
    '12': (-0.02, -0.02), 
    '14': (-0.02, -0.03), 
    '15': (-0.01, 0.035), 
    '16': (-0.01, 0.03), 
}
    # Annotate gauge locations with gauge IDs
    for idx, loc in discharge_locations_gdf.iterrows():
        point_geom = loc['geometry']
        point_x, point_y = point_geom.x, point_geom.y

        gauge_id = loc['index']

        # Plot gauge point
        ax.scatter(point_x, point_y, color='orange', s=30, edgecolor='black', zorder=5)

        # Annotate to the left with staggered vertical spacing
        dx = -0.02  # shift left
        dy = 0

        # Override if custom offset is provided
        if str(gauge_id) in custom_offsets:
            dx, dy = custom_offsets[str(gauge_id)]

        ax.text(point_x + dx, point_y + dy, str(gauge_id),
                fontsize=10, color='white', ha='right', va='center', fontweight='bold',
                bbox=dict(facecolor='black', edgecolor='none', pad=1.5, alpha=0.7))
    
    ax.set_title(f"Wflow gauge locations in {region.capitalize()} for Hurricane {runname}", fontsize=18)
    ax.set_axis_off()
    plt.tight_layout()

    if save:
        os.makedirs('p:/11210471-001-compass/04_Results/wflow', exist_ok=True)
        fig.savefig(f'p:/11210471-001-compass/04_Results/wflow/{region}_{runname}_wflow_gauge_map.png', dpi=300)
    else:
        plt.show()


#%%
# Define custom offsets for each location (you can specify the offset per gauge ID)
custom_offsets = {
    '1': (0.5, 0.5),  # Custom offset for gauge 1
    '2': (0.5, -1),  # Custom offset for gauge 2
    '3': (-2.5, 0.8),  # Custom offset for gauge 3
    '4': (-2.5, -1.3), # Custom offset for gauge 4
}
region = 'sofala'
event = 'Idai'

# Example discharge locations: [(x_pos, y_pos, (width, height), name)]
discharge_locations_gdf = gpd.read_file(f'p:/11210471-001-compass/02_Models/{region}/{event}/wflow/staticgeoms/gauges_locs.geojson')
mod = WflowModel(f'p:/11210471-001-compass/03_Runs/{region}/{event}/wflow/event_precip_era5_hourly_zarr_CF0/events', mode='r')
mod.read_results()


# Example timeseries (use actual discharge data for these)
discharge_timeseries = discharge_timeseries = mod.results['netcdf']

# Plot region with timeseries insets
plot_region_with_satellite('sofala', 'Idai', discharge_locations_gdf, discharge_timeseries, save=True, custom_offsets=custom_offsets)

#%%

plot_region_with_gauges_only(region="sofala", runname="Idai", discharge_locations_gdf=discharge_locations_gdf)


# %%

def plot_max_discharge_map(dataset, title="Maximum Discharge Map", save=False, filename=None):
    """
    Function to plot the maximum river discharge on a map using data from a xarray.Dataset.
    The basemap will be added using contextily for better background imagery.
    
    Parameters:
    - dataset: xarray Dataset with river discharge data
    - title: Title for the plot
    - save: Whether to save the plot
    - filename: If save is True, the filename to save the plot
    """
    # Extract max discharge over time
    max_discharge = dataset['q_river'].max(dim='time')  # Max across the time dimension
    
    # Plotting
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot the maximum discharge data as a heatmap
    im = ax.pcolormesh(dataset['lon'], dataset['lat'], max_discharge, cmap='YlGnBu', shading='auto')
    
    # Add a colorbar
    cbar = plt.colorbar(im, ax=ax, orientation='vertical', shrink=0.8)
    cbar.set_label('Max Discharge (m³/s)', fontsize=12)

    # Set the title
    ax.set_title(title, fontsize=16, fontweight='bold')

    # Add gridlines for latitude and longitude
    ax.grid(True, which='both', color='black', linewidth=0.5, linestyle='--')
    
    # Set the extent of the map (you can adjust the boundaries based on your data)
    ax.set_xlim([32.3, 35.4])  # Longitude range for Sofala region
    ax.set_ylim([-21.0, -17.5])  # Latitude range for Sofala region
    
    # Add basemap using contextily (OSM)
    ctx.add_basemap(ax, crs='EPSG:4326', source=ctx.providers.Esri.WorldImagery)

    # Optional: Add a marker or annotation for specific locations, like cities or river gauging points.
    # For example, you can add a text label on the plot:
    ax.text(34.8, -19.8, "Beira", fontsize=12, ha='center', color='white', fontweight='bold')

    # Save the plot if requested
    if save:
        if filename is None:
            filename = "max_discharge_map.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Saved plot as {filename}")

    # Show the plot
    plt.show()

# Example of how to call the function:
# Assuming 'dataset' is the xarray Dataset you provided
# plot_max_discharge_map(dataset, title="Maximum Discharge During Hurricane Idai", save=True, filename="max_discharge_idai.png")
# %%
dataset= mod.results['output']
plot_max_discharge_map(dataset)
# %%

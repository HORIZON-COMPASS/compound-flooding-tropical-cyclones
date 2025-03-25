#%%
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx
import matplotlib.patches as mpatches
import os
#%%



def plot_region_with_satellite(region, runname, save=True):
    geojson_path = f'p:/11210471-001-compass/02_Models/{region}/{runname}/wflow/staticgeoms/basins.geojson'
    geojson_path_sfincs = f'p:/11210471-001-compass/03_Runs/sfincs_{runname}/gis/region.geojson'
    geojson_path_rivers = f'p:/11210471-001-compass/02_Models/{region}/{runname}/wflow/staticgeoms/rivers.geojson'
    # Read the GeoJSON file
    gdf = gpd.read_file(geojson_path)
    
    
    gdf2 = gpd.read_file(geojson_path_sfincs)
    gdf2 = gdf2.to_crs(gdf.crs)  # Reprojecting to match gdf_1's CRS
    rivers = gpd.read_file(geojson_path_rivers).to_crs(gdf.crs) 
    print(gdf2)
    # Plot the merged region
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal', 'box') 
    gdf.plot(ax=ax, edgecolor='black', facecolor='lightblue', alpha = 0.7)
    gdf2.plot(ax=ax, edgecolor='black', facecolor='red', alpha = 0.7)
    rivers.plot(ax=ax, edgecolor='blue', facecolor='none', alpha = 1, lw=0.4)

    # Adjust the plot limits to zoom out a bit
    minx, miny, maxx, maxy = gdf.total_bounds
    buffer = 1  # Zoom-out buffer in degrees (adjust this as needed)
    ax.set_xlim(minx - buffer, maxx + buffer)
    ax.set_ylim(miny - buffer, maxy + buffer)  
    # Add the satellite background (OSM in this case)
    ctx.add_basemap(ax, crs=gdf.crs.to_string(), source=ctx.providers.Esri.WorldImagery)
    

    


 
    # Add a North Arrow
    north_arrow = mpatches.FancyArrowPatch(
        (0.9, 0.05), (0.9, 0.1), transform=ax.transAxes, color='white', arrowstyle='->', lw=2, mutation_scale=15
    )
    ax.add_patch(north_arrow)
    ax.text(0.91, 0.07, "N", transform=ax.transAxes, fontsize=14, color='white', va='bottom', ha='left', fontweight='bold')
    
    legend_elements = [
        mpatches.Patch(facecolor='lightblue', label='Wflow region', edgecolor='black'),
        mpatches.Patch(facecolor='red', label='SFINCS region', edgecolor='black')
    ]
    
     # Add the legend to the plot
    ax.legend(handles=legend_elements, loc='upper right', fontsize=12, title="Regions", frameon=True, facecolor='white', edgecolor='black')
    
    # Adjust plot appearance
    ax.set_title(f"Model domains {region.capitalize()}", fontsize=20)
    ax.set_axis_off()  # Turn off the axis
    plt.tight_layout()
    if save:
        if not os.path.exists(f'p:/11210471-001-compass/04_Results/model_regions'):
            os.mkdir(f'p:/11210471-001-compass/04_Results/model_regions')
        fig.savefig(f'p:/11210471-001-compass/04_Results/model_regions/{region}_{runname}_regions.png')
    # Show the plot
    else:
        plt.show()


# Replace this with the path to your GeoJSON file
#%%
plot_region_with_satellite('sofala', 'Idai')
#%%
plot_region_with_satellite('quelimane', 'Freddy2')

# %%

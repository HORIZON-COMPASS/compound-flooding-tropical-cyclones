import geopandas as gpd
import pandas as pd
from pathlib import Path
# by @dumontgoulart
#TODO: this is a preliminary script that is still not fully embedded in the snakemake workflow. The user is required to change the events & input files for now.

def main():
    """
    Main function to calculate and save relative damage.
    """
    event_name = "Freddy"  # Change this to your event name
    print(f"Processing event: {event_name}")
    directory = Path(f"/p/11210471-001-compass/03_Runs/test/{event_name}/fiat/event_tp_era5_hourly_CF0_GTSMv41opendap_CF0_no_wind_CF0")
    spatial_path = directory / "output" / "spatial.fgb"
    exposure_path = directory / "exposure" / "exposure.csv"
    output_dir = directory / "output"

    try:
        # Load the spatial data
        print(f"Reading spatial data from: {spatial_path}")
        spatial_gdf = gpd.read_file(spatial_path)
        spatial_gdf = spatial_gdf[spatial_gdf['total_damage'] > 0]

        # Load the exposure data
        print(f"Reading exposure data from: {exposure_path}")
        exposure_df = pd.read_csv(exposure_path)

        # Select only necessary columns from exposure data to avoid memory issues
        exposure_df_subset = exposure_df[['object_id', 'max_damage_total']]

        # Merge the GeoDataFrame with the exposure DataFrame
        print("Merging dataframes on 'object_id'...")
        merged_gdf = spatial_gdf.merge(exposure_df_subset, on="object_id", how="left")

        # Calculate 'relative_damage'
        print("Calculating 'relative_damage'...")
        # Handle cases where max_damage_total is 0 or NaN to avoid errors
        merged_gdf['max_damage_total'] = pd.to_numeric(merged_gdf['max_damage_total'], errors='coerce').fillna(0)
        merged_gdf['total_damage'] = pd.to_numeric(merged_gdf['total_damage'], errors='coerce').fillna(0)

        # Avoid division by zero
        merged_gdf['relative_damage'] = merged_gdf.apply(
            lambda row: row['total_damage'] / row['max_damage_total'] if row['max_damage_total'] != 0 else 0,
            axis=1
        )

        # Save the result
        output_path = Path(output_dir) / "output_relative_damage.fgb"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        print(f"Saving updated spatial data to: {output_path}")
        merged_gdf.to_file(output_path, driver='FlatGeobuf')

        print("Script finished successfully.")

    except FileNotFoundError as e:
        print(f"Error: {e}. Please check your file paths.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()

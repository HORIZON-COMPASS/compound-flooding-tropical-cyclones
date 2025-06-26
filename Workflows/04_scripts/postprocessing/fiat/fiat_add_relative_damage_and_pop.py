import geopandas as gpd
import pandas as pd
import argparse
from pathlib import Path

def main():
    """
    Main function to calculate and save relative damage.
    """
    parser = argparse.ArgumentParser(
        description="Calculates relative damage by merging FIAT's spatial output with exposure data."
    )
    parser.add_argument(
        "--spatial_path",
        type=str,
        required=True,
        help="Path to the spatial file (e.g., spatial.fgb) with 'object_id' and 'total_damage'.",
    )
    parser.add_argument(
        "--exposure_path",
        type=str,
        required=True,
        help="Path to the exposure CSV file with 'object_id' and 'max_damage_total'.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Path to the directory to save the output file.",
    )

    args = parser.parse_args()

    try:
        # Load the spatial data
        print(f"Reading spatial data from: {args.spatial_path}")
        spatial_gdf = gpd.read_file(args.spatial_path)

        # Load the exposure data
        print(f"Reading exposure data from: {args.exposure_path}")
        exposure_df = pd.read_csv(args.exposure_path)

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
        output_path = Path(args.output_dir) / "output_relative_damage.fgb"
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

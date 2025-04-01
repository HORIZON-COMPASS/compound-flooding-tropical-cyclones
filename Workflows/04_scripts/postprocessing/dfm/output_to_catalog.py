#%% Adding DFM output to the SFINCS coastal coupling data catalog
# Importing the necessary packages
import os
import hydromt

#%%
if "snakemake" in locals():
    his_path = os.path.abspath(snakemake.input.his_file)
    path_data_cat = os.path.abspath(snakemake.params.cf_data_cat)
    model_name = snakemake.params.model_name
    root_dir = os.path.abspath(snakemake.params.root_dir)
    snake_done = os.path.abspath(snakemake.output.done_file)
else:
    region = "sofala"
    tc_name = "Idai"
    dfm_res = "450"
    bathy = "gebco2024_MZB"
    tidemodel = 'GTSMv41opendap' # tidemodel: FES2014, FES2012, EOT20, GTSMv41, GTSMv41opendap
    wind_forcing = "spw_IBTrACS"
    CF_SLR_txt = "-0.14"
    CF_wind_txt = "0"
    model_name = f'event_{dfm_res}_{bathy}_{tidemodel}_CF{CF_SLR_txt}_{wind_forcing}_CF{CF_wind_txt}'
    path_data_cat = os.path.abspath("../../../03_data_catalogs/datacatalog_CF_forcing.yml")
    run_dir = f'p:/11210471-001-compass/03_Runs/{region}/{tc_name}/dfm/{model_name}'
    root_dir = 'p:\\'
    snake_done = os.path.join(run_dir, "postprocessing_done.txt")
    his_path1 = os.path.join(run_dir, "output", f"{model_name}_0000_his.nc")
    his_path2 = os.path.join(run_dir, "output", f"settings_0000_his.nc")  # Alternative file name
    # Use the first path that exists
    his_path = his_path1 if os.path.exists(his_path1) else his_path2
    

#%% Loading the SFINCS coastal coupling data catalog & DFM output path
datacatalog = hydromt.DataCatalog(data_libs=[path_data_cat])
dfm_run = f"dfm_output_{model_name}"

#%% Specifying the DFM output information for the data catalog entry
adapter = hydromt.data_adapter.GeoDatasetAdapter(
    path=os.path.abspath(his_path),
    driver="netcdf",
    driver_kwargs={
        "chunks": {
            "stations": 10,
            "time": -1
        }
    },
    rename={
        "station_x_coordinate": "lon",
        "station_y_coordinate": "lat",
        "stations": "index"
    },
    meta={
        "category": "ocean"
    },
    crs=4326
    )

#%% Add the DFM output to the catalog and save
if dfm_run not in datacatalog:
    # Add new source if it doesn't exist
    datacatalog.add_source(dfm_run, adapter)
    print(f"Dataset '{dfm_run}' has been added to the catalog.")
else:
    # Update the existing source
    datacatalog[dfm_run] = adapter
    print(f"Dataset '{dfm_run}' has been updated in the catalog.")

# Save the updated catalog
datacatalog.to_yml(path_data_cat, root=root_dir)

#%% Make a file for snakemake to track the cata,log update
with open(snake_done, 'w') as file:
    # Write content to the file
    file.write("# Empty file used to make the snakemake dfm workflow add_data_to_catalog rule work\n")
# %%

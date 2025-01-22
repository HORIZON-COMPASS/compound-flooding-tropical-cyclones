#%% Adding DFM output to the SFINCS coastal coupling data catalog
# Importing the necessary packages
import os
import hydromt

#%%
if "snakemake" in locals():
    his_path = os.path.abspath(snakemake.input.his_file)
    path_data_cat = os.path.abspath(snakemake.output.sfincs_data_cat)
    run_dir = os.path.abspath(snakemake.params.dir_event_model)
    model_name = snakemake.params.model_name
    tc_name = snakemake.wildcards.tc_name
    root_dir = os.path.abspath(snakemake.params.root_dir)
else:
    region = "sofala"
    tc_name = "Idai"
    dfm_res = "450"
    bathy = "gebco2024"
    tidemodel = 'GTSMv41opendap' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap
    wind_forcing = "spw_IBTrACS_ext"
    model_name = f'event_{dfm_res}_{bathy}_{tidemodel}_{wind_forcing}'
    path_data_cat = os.path.abspath("../../../data_catalogs/datacatalog_SFINCS_coastal_coupling.yml")
    run_dir = f'p:/11210471-001-compass/03_Runs/{region}/{tc_name}/dfm/{model_name}'
    his_path = os.path.join(run_dir, "output", f"{model_name}_his.nc")
    root_dir = 'p:\\'

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
# %%

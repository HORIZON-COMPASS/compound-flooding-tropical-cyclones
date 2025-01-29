
The structure of the folder is as follows:
- Workflows
    - Config - per model  : Model configuration files
    - Config_snakemake : Snakemake configuration files
    - data_catalogs : The data catalogs used for model building and used by HydroMT
    - scripts
        - model_building : Scripts used to build the models
        - preprocessing : Scripts used for preprocessing of data, updating models etc.
        - postprocessing : Scripts used for postprocessing of model results


## **Workflows**

Four snakemake workflow files are present. All workflows work both on Linux as well as on Windows. All snakemake workflows use the same config yml file: config_snakemake/config_general.yml. 

- [Workflows](#workflows)
  - [snakefile_all.smk](#snakefile_allsmk): This workflow combines all other snakemake workflows into one large workflow. The sequence of the workflows are: snakefile_sfincs_build.smk > snakefile_wflow.smk > snakefile_sfincs_update.smk
  - [snakefile_sfincs_build.smk](#snakefile_sfincs_buildsmk): This workflow builds a SFINCS model without adding any forcing yet. The workflow consists of just one rule.
  - [snakefile_wflow.smk](#snakefile_wflowsmk): This workflow builds a wflow model, using the SFINCS region as input. Gauges are added on the SFINCS inflow points, in order to generate output at the correct locations. Precipitation forcing is added based on the event start and end time. First, the wflow model is warmed up for a period of 1 year, with daily ERA5 data. The event itself is run with the forcing data as given in the snakemake configuration file.
  - [snakefile_sfincs_update](#snakefile_sfincs_updatesmk): This workflow updates the SFINCS model by adding forcing data and running the model simulations. It handles the addition of both meteorological and WFlow forcing data, executes the model, and generates the output.
   - [snakefile_dfm.smk](#snakefile_dfmsmk):This workflow creates a dfm base model and updates it by adding forcing data and running the model simulations. It also add the output to a data catalog which can be used as SFINCS waterlevel forcing.

In the section below, the different workflows are described:

### Workflow --- snakefile_all.smk ---
This workflow combines all other snakemake workflows into one large workflow. The sequence of the workflows are: 

snakefile_sfincs_build.smk > snakefile_wflow.smk > snakefile_sfincs_update.smk

### Workflow --- snakefile_sfincs_build.smk ---

This workflow builds a SFINCS model without adding any forcing yet. The workflow consists of just one rule.

#### Rule: make_base_model_sfincs
- **Script**: `scripts/model_building/sfincs/setup_sfincs_base.py`
- **Input**: `config_sfincs/sfincs_base_build.yml`
- **Output**: SFINCS model for the region (bounding box defined in `config_snakemake/config_general.yml`).
- **Description**: This rule creates the SFINCS model  using the bbox given in the snakemakeThe model domain is based on the bounding box as given in the config_snakemake/config_general.yml file.
Further options for building the model are given in the *config_sfincs\sfincs_base_build.yml* input file.
The bounding box is then  used to find all intersecting subbasins (hydroAtlas level 12). These intersecting subbasins form the model region. Also, the model region extents to the -5???m contour line. 


### Workflow --- snakefile_wflow.smk ---

This workflow builds a wflow model, using the SFINCS region as input. Gauges are added on the SFINCS inflow points, in order to generate output at the correct locations. Precipitation forcing is added based on the event start and end time. First, the model is warmed up for a period of 1 year, with daily ERA5 data. The event itself is run with the forcing data as given in the snakemake config file.

#### Rule: make_base_model_wflow
- **Script**: `scripts/model_building/sfincs/setup_wflow_base.py`
- **Input**: 
  - `config_wflow/wflow_build_{REGION}.yml`
  - `{SFINCS_DIR}/gis/region.geojson`
  - `{SFINCS_DIR}/gis/src.geojson`
- **Output**: WFlow model for the SFINCS region.
- **Description**: This rule creates the wflow model  using the SFINCS region 
Further options for building the model are given in the input yml file. The model is created for the entire upstream basin of the SFINCS region.
The model is created on the P-drive in the *p:\11210471-001-compass\02_Models* folder

#### Rule: update_forcing_wflow_warmup
- **Script**: `scripts/preprocessing/update_forcing_wflow_warmup.py`
- **Input**:
  - `toml_file`
  - `staticmaps`
- **Output**: Updated WFlow base model with ERA5 forcing data.
- **Description**: This rule uses the wflow base model and updates it with daily forcing data from ERA5. The start time is 1 year before the start time of the event (given in the snakemake config yml file) and the end time is 2 days before the event. The updated model is saved to the P-drive in the *p:\11210471-001-compass\03_Runs* folder

#### Rule: update_forcing_wflow_event
- **Script**: `scripts/preprocessing/update_forcing_wflow_event.py`
- **Input**:
  - `toml_file`
  - `staticmaps`
- **Output**: Updated WFlow base model for the event with the configured forcing data.
- **Description**: This rule uses the wflow base model and updates it with forcing data as configured in the snakemake config file. The start time is 2 days before the date in the config file and it uses the initial state from the warmup run. The updated model is saved to the P-drive in the *p:\11210471-001-compass\03_Runs* folder

#### Rule: run_wflow_warmup
- **Input**:
  - `inmaps.nc` (from warmup run)
  - `toml_file` (from warmup run)
- **Output**: WFlow warmup run output.
- **Description**: This rule runs the wflow warmup model. The final state is saved into the folder of the 'event' run. 

#### Rule: run_wflow_event
- **Input**:
  - `inmaps.nc` (from event run)
  - `toml_file` (from event run)
  - `instates.nc` (from warmup run)
- **Output**: WFlow event run output, saved as `output_scalar.nc`.
- **Description**: This rule runs the wflow event run. Output is saved as *output_scalar.nc*


### Workflow --- snakefile_dfm.smk ---

This workflow creates a dfm base model and updates it by adding forcing data and running the model simulations. It also add the output to a data catalog which can be used as SFINCS waterlevel forcing.

#### Rule: make_model_dfm_base
- **Script**: `scripts/model_building/dfm/setup_dfm_base.py`
- **Params**:
  - `bbox_dfm` (from config file)
  - `bbox_output` (from config file)
  - `data_cat` (from snake file)
- **Output**: DFM base model configuration, saved as multiple files in the base model folder:
  - `geometry` folder with base model grid figures
  - `pli_file.pli`
  - `illigalcells.pol`
  - `grid_network.nc`
  - `tide_{bathy}.bc`
  - `L*.pli` where * is a number, depending on how many boundary lines there are
  - `ext_file_new.ext`
- **Description**: This rule makes generates and refines base model grid and adds boundary conditions from a tidal model.

#### Rule: make_dfm_model_event
- **Script**: `scripts/model_building/dfm/setup_dfm_event.py`
- **Input**: 
  - `ext_file_new.ext` (from the base model)
  - (and all other files in the base model folder, which are copied inside the script to the event folder)
- **Params**:
  - `dir_base_model` (from make_model_dfm_base rule)
  - `start_time` (from config file)
  - `end_time` (from config file)
  - `bbox_dfm` (from config file)
  - `bbox_output` (from config run)
  - `dfm_obs_file` (from config file)
  - `verif_points_file` (from config file)
  - `data_cat` (from snake file)
  - `dimrset` (from Deltares dfm directory)
  - `uniformwind` (from data directory)
  - `model_name` (from snake rule)
- **Output**: 
  - `geometry` folder with base grid ánd event model figures
  - `{model_name}.mdu`
  - `ext_file_old.ext`
  - `dimr_config.xml`
  - `submit_singularity_h7.sh` (if run on linux)
  - `run_parallel.bat` (if run on windows)
- **Description**: This rule copies the base model configuration and adds the TC meteo forcing. It also set the model output points (through an obs file) and generates the files required to run the model (mdu, dimr, bat and sh)

#### Rule: run_dfm
- **Script**: `scripts/model_building/dfm/setup_dfm_event.py`
- **Input**: 
  - `submit_singularity_h7.sh` (if run on linux)
  - `run_parallel.bat` (if run on windows)
  - (and all other files in the event model folder, needed for model execution)
- **Run**:
  - Runs the batch file depending on the system that is used
- **Output**: 
  - `output/{model_name}_his.nc` (needed as sfincs water level forcing)
  - other output files, not needed for sfincs and therefore not included in the snake rule
- **Description**: This rule runs the dfm event model configuration, possible on windows and linux systems. It generates an output folder containing a his.nc, map.nc and .dia file. 

#### Rule: add_dfm_output_to_catalog
- **Script**: `scripts/postprocessing/dfm/output_to_catalog.py`
- **Input**: 
  - `output/{model_name}_his.nc` (from the event model directory)
- **Params**:
  - `model_name` (from snake rule)
  - `sfincs_data_cat` (from snake rule)
  - `root_dir` (from snake rule, differs for windows or linux)
- **Output**: 
  - `done_file` tracks whether the sfincs_data_cat is updated
- **Description**: This rule adds the event model his.nc output to the SFINCS data catalog, to be used as water level forcing for the SFINCS model. A 'done_file' is used for Snakemake to track whether the data catalog has been updated, but has no other purpose.


### Workflow --- snakefile_sfincs_update.smk ---

This workflow updates the SFINCS model by adding forcing data and running the model simulations. It handles the addition of both meteorological and WFlow forcing data, executes the model, and generates the output.

#### Rule: `add_forcing_coastal_meteo_sfincs`
- **Script**: `scripts/preprocessing/update_sfincs_coastal_forcing.py`
- **Input**: 
  - msk file
  - spw file
- **Output**: 
  - bzs file
- **Description**: This rule adds meteorological forcing from spiderweb data to the SFINCS model. The forcing data is incorporated into the model configuration (`sfincs.bzs`), enabling simulation based on coastal meteorological inputs.

#### Rule: `update_dis_forcing_sfincs`
- **Script**: `scripts/preprocessing/update_sfincs_dis_forcing.py`
- **Input**: 
  - bzs file
  - wflow output file
- **Output**: 
  - dis file
- **Description**: This rule updates the SFINCS model by adding forcing data derived from WFlow outputs. The `.dis` file is created, adding the discharge from wflow as forcing for SFINCS

#### Rule: `run_sfincs_model`
- **Script**: N/A (The model is executed via a subprocess or Docker)
- **Input**: 
  - dis file
- **Output**: 
  - sfincs_map.nc
  - sfincs_his.nc
- **Description**: This rule runs the SFINCS model with the updated `.dis` file. Depending on the operating system, the model is executed either via a subprocess (on Windows) or through Docker (on Linux).

#### Rule: `sfincs_plot_floodmap`
- **Script**: `scripts/postprocessing/sfincs_postprocess.py`
- **Input**: 
  - sfincs_map.nc
  - sfincs_his.nc
- **Output**: 
  - sfincs_basemap.png
- **Description**: This rule generates the final flood map visualization from the model outputs (`sfincs_map.nc` and `sfincs_his.nc`). The postprocessing script is used to produce the flood map as well as some other images.

---

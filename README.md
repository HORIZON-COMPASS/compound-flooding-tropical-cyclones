# Compound Flooding Modeling for Tropical Cyclones

This work is part of Work Package 1 of the COMPASS project whose overarching objective is to characterise compound extremes in current and future climates. COMPASS (COMPound extremes Attribution of climate change: towardS an operational Service) aims to develop a harmonized, yet flexible, methodological framework for **climate and impact attribution** of various complex **extremes** that include compound, sequential and cascading hazard events. For more information and useful links about the project, have a look at the introduction on the [COMPASS Github repository](https://github.com/HORIZON-COMPASS)

![logoCOMPASS](https://github.com/user-attachments/assets/4c3b95d4-bfc0-4727-a1e8-ee6653a03b5e)

## Description
This repository contains workflows to setup a global-to-local modelling chain for the hazard and impact modelling of tropical cyclones (TCs). More specifically, it contains specific workflows to setup a hydrological model, a coastal hydrodynamic model and a local compound flooding model for a user-defined region. Each model is setup through a separate workflow and provide boundary conditions for a specific TC driver to be considered to model local flooding. A final workflow combines all the separate workflows together for the local compound flood model setup. The workflows are created using the [snakemake](https://snakemake.readthedocs.io/en/stable/) package and [hydroMT](https://deltares.github.io/hydromt/stable/) model builders for wflow ([hydroMT-Wflow](https://deltares.github.io/hydromt_wflow/latest/)) and SFINCS ([hydromt-sfincs](https://deltares.github.io/hydromt_sfincs/latest/)) and the Python [dfm-tools](https://deltares.github.io/dfm_tools/) for the Delft3D-FM model. All the input datasets  are defined in a data catalog and specific model parameters and other pre-processing steps in a configuration file, both are .yml files. 

![image](https://github.com/user-attachments/assets/d50faecf-2b06-4193-8780-150476cb8315)

The code is under development and this file will be updated to reflect the current updates on the project. 

## Installation instructions
### Prerequisite
Make sure you have [Git](https://github.com/git-guides/install-git) installed. You will also need [pixi](https://pixi.sh/latest/#installation) to install the required dependencies (defined in specific environments). We suggest using [Visual Studio Code](https://code.visualstudio.com/) to browse through the code, develop and open the Jupyter notebooks (optional).

### Installation
In order to use the workflows, you will need to clone this repository and install the dependencies required to run the code. 
1. Open a terminal and clone this repository: `git clone https://github.com/HORIZON-COMPASS/compound-flooding-tropical-cyclones.git`
2. Navigate to the project directory: `cd compound-flooding-tropical-cyclones`
3. Install dependencies: `pixi install`. This will install all the environments required to run all the workflows. To install only specific environments, mention it, for example : `pixi install compass-sfincs`. All the enviroments available are listed [environments] in the pixi.toml file
   
## How to run
At the moment, four snakemake workflow files are present. All workflows work both on Linux as well as on Windows. A workflow consists of a set of rules to be executed in a specific order. All specific workflow configuration settings are prescribed in a configuration file (a .yml file).

- **snakefile_sfincs_build.smk**: This workflow builds a SFINCS model without adding any forcing yet. The workflow consists of just one rule.
- **snakefile_wflow.smk**: This workflow builds a wflow model, using the SFINCS region as input. Gauges are added on the SFINCS inflow points, in order to generate output at the correct locations. Precipitation forcing is added based on the event start and end time. First, the wflow model is warmed up for a period of 1 year, with daily ERA5 data. The event itself is run with the forcing data as given in the snakemake configuration file.
- **snakefile_dfm.smk**: This workflow creates a dfm base model and updates it by adding forcing data and running the model simulations. It also add the output to a data catalog which can be used as SFINCS waterlevel forcing.
- **snakefile_sfincs_update.smk**: This workflow updates the SFINCS model by adding forcing data and running the model simulations. It handles the addition of both meteorological and WFlow forcing data, executes the model, and generates the output.
- **snakefile_all.smk**: This workflow combines all other snakemake workflows into one large workflow. The sequence of the workflows are: 
snakefile_sfincs_build.smk > snakefile_wflow.smk > snakefile_sfincs_update.smk

All snakemake workflows use the same configuration file: config_snakemake/config_general.yml.

We provide examples on how to run each workflow in specific Jupyter notebook in the docs folder. Building the workflow requires the same general steps in snakemake, summarized below.
 1. Activating the environments to load all the required dependencies
 2. Navigate to the folder where the workflows are located (or make sure that you use the correct path afterwards)
 3. Creating a picture showing the workflow with the rules and their order to make sure the workflow is doing what is intended (optional)
 4. Unlocking the working directory in order to save the results of the workflow
 5. Running the workflow

This translate to the following command lines in the terminal, taking as an example the snakefile_wflow.smk workflow:
```
conda activate compass-wflow
cd Workflows
snakemake -s snakefile_wflow --configfile config_snakemake/config_general.yml  --dag | dot -Tpng > dag_all.png
snakemake --unlock -s snakefile_wflow --configfile config_snakemake/config_general.yml
snakemake all -c 1 -s snakefile_wflow --configfile config_snakemake/config_general.yml
```

There exists many snakemake commad line options that are worth exploring. The complete list is is on the [CLI documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) but below are a few useful ones:
- **-s**: selection of the snakefile workflow to run.
- **--configfile**: name of the config file with the model and climate options.
- **-c**: number of cores to use to run the workflows (if more than 1, the workflow will be parallelized).
- **--dry-run**: returns the list of steps (rules) in the workflow that will be run, **without actually running it**. This is useful to test whether the workflow will work as intended. 

## How to contribute
We welcome contributions to improve this project! Here are some ways you can help:
**Report Bugs**: If you find a bug, please open an issue with detailed information about the problem and how to reproduce it.
**Submit Pull Requests**: If you want to fix a bug or implement a feature, follow these steps:
1. Fork the repository.
2. Create a new branch (`git checkout -b feature/YourFeatureName`).
3. Make your changes.
4. Commit your changes (`git commit -m 'Add some feature'`).
5. Push to the branch (`git push origin feature/YourFeatureName`).
6. Open a pull request.
**Suggest Features**: Have an idea for a new feature? Open an issue to discuss it.

## Acknowledgements

![EU_logo](https://github.com/user-attachments/assets/e2fad699-697e-43fd-84be-032447d6dd21) The COMPASS project has received funding from the European Union’s HORIZON Research and Innovation Actions Programme under Grant Agreement No. 101135481

Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or of the European Health and Digital Executive Agency (HADEA). Neither the European Union nor the granting authority HADEA can be held responsible for them.

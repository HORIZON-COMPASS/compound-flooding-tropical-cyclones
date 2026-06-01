# Scripts and data to reproduce the figures from Vertegaal et al. (submitted)

## Data
For data from the models, the files are provided in the *data* folder that first has to be unzipped. 

For data from the data catalogs (e.g. gswo), links to download the data are provided in the data catalogs (see Workflows/03_data_catalogs).

For the data of 30 year wflow runs, the authors can be contacted (see email below), since the files are quite large. This data is necessary for the scripts *Wflow_bankfull_removal.py* & *Wflow_GRDC_comparison.py*.

For questions, you can contact the authors (main contact: doris.vertegaal@deltares.nl).


## Scripts
The script sin the folder *scripts* can be used to reproduce the figures and numbers presented in Vertegaal et al. (submitted). On top of every script is stated which Pixi environment from the pixi.toml file can be used to run the scripts. See the main *README.md* on how to use Pixi for Python environments. The figures and table are saved according to their figure or table number in the paper. See the Table below for which script is used for the production of which asset:

| Asset                 | Path                                                                   |
|-----------------------|------------------------------------------------------------------------|
| F01                   | *NA*                                                                   |
| F02                   | lulc_processing.ipynb                                                  |
| F03                   | Kumu                                                                   |
| F04                   | factual_counterfactual_maps.py                                         |                                                                               
|                                                                                                |
| *Supplement*                                                                                   |
| FS1                   |  *NA*                                                                  |
| FS2 (request data)    | River_discharge_analysis.ipynb                                         |

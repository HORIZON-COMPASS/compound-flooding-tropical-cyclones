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
| F02                   | Attribution_results\scripts\Model_regions.py                           |
| F03                   | Attribution_results\scripts\Flood_damage_and_attribution.py            |
| F04                   | Attribution_results\scripts\Attribution_driver_combinations.py         |
| F05                   | Attribution_results\scripts\Attribution_driver_combinations.py         |
| F06                   | Attribution_results\scripts\Flood_damage_and_attribution.py            |
| Table 1               | *NA*                                                                   |
| Table 2               | Attribution_results\scripts\Attribution_driver_combinations.py         |
|                                                                                                |
| *Supplement*                                                                                   |
| FS1                   | Attribution_results\scripts\SFINCS_basemap.py                          |
| FS2                   | Attribution_results\scripts\Wflow_bankfull_removal                     |
| FS3                   | Attribution_results\scripts\Wflow_bankfull_removal                     |
| FS4                   | Attribution_results\scripts\Wflow_GloFAS_comparison.py                 |
| FS5                   | Attribution_results\scripts\Wflow_GloFAS_comparison.py                 |
| FS6                   | Attribution_results\scripts\Wflow_GloFAS_comparison.py                 |
| FS7                   | Workflows\04_scripts\postprocessing\dfm\output_to_catalog_add_waves.py |
| FS8                   | Attribution_results\scripts\Validation_floodmaps_EO_data.py            |
| FS9                   | Attribution_results\scripts\Calculating_SLR.py                         |                       |
| FS10                  | Attribution_results\scripts\Attribution_driver_combinations.py         |
| FS11                  | Attribution_results\scripts\Factual_flooding_and_timeseries.py         |
| FS12                  | Attribution_results\scripts\Flood_damage_and_attribution.py            |
| Table S1              | Attribution_results\scripts\Wflow_bankfull_removal.py                  |
| Table S2              | Attribution_results\scripts\Attribution_driver_combinations.py         |
|                                                                                                |
| *In-text numbers*                                                                              |
| Fact flood extent     | Table 2                                                                |
| Fact discharge        | Attribution_results\scripts\Factual_flooding_and_timeseries.py         |
| Discharge normal      | Table S1                                                               |
| Mean accum rain       | Attribution_results\scripts\Rainfall_forcing_stats.py                  |
| Max coast water level | Attribution_results\scripts\Factual_flooding_and_timeseries.py         |
| Wave setup            | Workflows\04_scripts\postprocessing\dfm\output_to_catalog_add_waves.py |
| Flooded buildings     | Attribution_results\scripts\Flood_damage_and_attribution.py            |
| Total damage          | Attribution_results\scripts\Flood_damage_and_attribution.py            |
| Damage in Beira       | Attribution_results\scripts\Flood_damage_and_attribution.py            |
| Affected CF SLR & Wind area | Attribution_results\scripts\Attribution_driver_combinations.py   |
| Affected CF Rain area | Attribution_results\scripts\Attribution_driver_combinations.py         |

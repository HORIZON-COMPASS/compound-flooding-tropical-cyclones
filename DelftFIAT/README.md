# **A Delft-FIAT model for COMPASS**
Within COMPASS, the Delft-FIAT model can be used to calculate damages of compound flooding from TCs in Mozambique

*Made by Doris Vertegaal for the COMPASS project on 12-11-2024*

These scripts are based on the online examples of HydroMT-FIAT and Delft-FIAT 
HydroMT-FIAT: https://deltares.github.io/hydromt_fiat/latest/index.html
Delft-FIAT: https://deltares.github.io/Delft-FIAT/stable/

#### First, install the compass-fiat-hydromt and compass-fiat environments from Pixi:
* install pixi on your computer: open Windows PowerShell, and run (source: https://pixi.sh/latest/):
`iwr -useb https://pixi.sh/install.ps1 | iex`
* Navigate to your repository folder in PowerShell (e.g., C:/COMPASS/)
* Run the following in PowerShell:
`pixi install --environment compass-fiat-hydromt` & `pixi install --environment compass-fiat`


#### Second, build the fiat model using HydroMT-FIAT using: 
* scripts/fiat_hydromt.py 
* use the pixi environment compass-fiat-hydromt
* To run this model for your own event, give the path to your floodmap.tif and the correct region. Also make sure to change the country and continent.

This script is developed using the global example from the Hydromt_FIAT GitHub:
https://github.com/Deltares/hydromt_fiat/blob/main/examples/global_OSM_JRC.ipynb

The data described in the data catalog is copied from the data folder in the same GitHub


#### Third, run your fiat model using:
* scripts/execute_fiat_example.ipyn & activate your compass-fiat environment
Alternatively: You can also run fiat in Powershell by activating your environment and running:
! fiat run "{path to COMPASS repository folder}/computations/{model_name}/settings.toml"


#### Last, you can check the results your model run using:
* scripts/check_output.ipyn

Contributions:
Sarah Rautenbach for providing support in building and running the model for Mozambique 
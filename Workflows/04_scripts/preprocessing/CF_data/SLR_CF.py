# In this script, the DFM waterlevel forcing is adjusted to a 
# counterfactual scenario using ISIMIP SLR data  
#%%
# Import the correct packages
import hydromt 
import matplotlib.pyplot as plt
import numpy as np

#%% Region of case study area in Mozambique
lat_MZB = [-20.12, -19.30]
lon_MZB = [34.33, 34.95]

# region of eastern Africa
lat_MZ = [-27, -9]
lon_MZ = [29, 46]

# case study settings
start_date = np.datetime64('2019-03-09T00:00') 
end_date = np.datetime64('2019-03-24T00:00') 
#%% Reading in ISIMIP data
cat = hydromt.DataCatalog("../../../03_data_catalogs/datacatalog_SFINCS_coastal_coupling.yml")
# %%
dfm_idai = cat.get_dataset("dfm_output_MZB_Idai")
dfm_idai_CF = dfm_idai.copy()

# %%
dfm_idai_CF['waterlevel'] = dfm_idai['waterlevel'] - 0.141
# %% Export the counterfactual ds
outfile_path = "p:/11210471-001-compass/01_Data/counterfactuals/SLR/ISIMIP/"

dfm_idai_CF.to_netcdf(outfile_path + "dfm_output_MZB_Idai_SLR2015.his")
dfm_idai_CF.close()

# %%
# Plot the mean data
# Set up the subplots
fig, ax = plt.subplots(figsize=(15, 10))
# mean_data_F.plot(ax=ax)
dfm_idai_CF['waterlevel'].plot(ax=ax)

# %% Select water level data for one station
# Create the plot
fig, ax = plt.subplots(figsize=(15, 10))
ax.plot(dfm_idai_CF['time'], dfm_idai_CF['waterlevel'].sel(stations=[4]), label=f"CF")
ax.plot(dfm_idai['time'], dfm_idai['waterlevel'].sel(stations=[4]), label=f"F")

# Add labels and title
ax.set_xlabel("Time")
ax.set_ylabel("Water Level")
# ax.set_title(f"Water Level Forcing Over Time for {station_name}")
ax.legend()
plt.show()
# %%

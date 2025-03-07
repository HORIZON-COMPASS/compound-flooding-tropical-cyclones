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
cat = hydromt.DataCatalog("../../../03_data_catalogs/datacatalog_CF_forcing.yml")
# %%
# dfm_idai = cat.get_dataset("dfm_output_MZB_Idai")
dfm_Idai_CF0 = cat.get_dataset("dfm_output_event_450_gebco2024_MZB_GTSMv41opendap_CF0_spw_IBTrACS_CF0")
dfm_Idai_CF0_copy = dfm_Idai_CF0.copy()

# %%
dfm_Idai_CF0_copy['waterlevel'] = dfm_Idai_CF0_copy['waterlevel'] - 0.14
# %% Export the counterfactual ds
outfile_path = "p:/11210471-001-compass/01_Data/counterfactuals/SLR/ISIMIP/"

dfm_Idai_CF0_copy.to_netcdf(outfile_path + "dfm_output_event_450_gebco2024_MZB_GTSMv41opendap_CF0_spw_IBTrACS_CF0_SLR2015.his")
dfm_Idai_CF0_copy.close()

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
<<<<<<< HEAD
dfm_Idai_CF0_copy['waterlevel'].plot(ax=ax)
=======
dfm_idai_CF['waterlevel'].plot(ax=ax)
>>>>>>> 6064949 (adding SLR processing scripts)

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

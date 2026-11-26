#%%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt










#####################################################################
###############################      WIND      ######################
#####################################################################
fig, ax = plt.subplots(figsize=(8, 5))

for exp in ["cont", "hist", "Tplus2.0K"]:
    sns.kdeplot(data=df_climatedt_wind_max[df_climatedt_wind_max["experiment"] == exp], 
                x="max_wind", fill=False, common_norm=False, label=exp, ax=ax)

sns.kdeplot(data=ibracs_wind_max, x="max_wind_ms", fill=False, common_norm=False, 
            label="IBTrACS", ax=ax)

ax.set_xlabel("Maximum wind speed [m/s]")
ax.set_ylabel("Probability density")
ax.legend()
plt.show()
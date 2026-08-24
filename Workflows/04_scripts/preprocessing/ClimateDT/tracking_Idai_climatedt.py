
#%%
import xarray as xr
from idai_funcs import (import_ibtracs, plot_tracks, main)

track = import_ibtracs(storm="IDAI", year=2019)

if __name__=='__main__':
    track, data_all, tc_idai_climatedt = main(run_tracking=True)

plot_tracks(track, tc_idai_climatedt, data_all.scenario.values)


# %% Step 1: Run the full analysis pipeline
from run_population_analysis import *  # noqa: F403

# %% Step 2: Import all plot functions
from plot_attribution_figures import *  # noqa: F403

# %% Fig 2: Attribution maps
fig = plot_fig2_attribution_maps(
    gdf_pop_2019_exposed_F_coarse, gdf_pop_2019_exposed_CF_coarse,
    gdf_pop_1990_exposed_F_coarse, gdf_pop_1990_exposed_CF_coarse,
    region_utm, background_utm, beira_utm, flood_extent,
)
plt.show()

#%%
# Fig Sx: Factual flood, population, exposed population and changes (6-panel)
fig = plot_supp_factual_changes_overview(hmax_F, gdf_pop_2019_exposed_F_coarse, hmax_diff, 
                                         gdf_pop_1990_exposed_F_coarse, gdf_pop_1990_exposed_CF_coarse, 
                                         region_utm, background_utm, flood_extent)
plt.show()

# %% Attributable % per flood depth
fig = plot_attributable_pct_per_depth(
    bin_centers_coarse,
    pop_2019_by_depth_F_coarse, pop_2019_by_depth_CF_coarse,
    pop_1990_by_depth_F_coarse, pop_1990_by_depth_CF_coarse,
)
plt.show()

# %% Bar chart attribution
fig, diff_clim, diff_pop, diff_clim_pop, pct_cc, pct_pop, pct_both = plot_absolute_diff_per_depth(
    bin_centers,
    pop_2019_by_depth_F_fine, pop_2019_by_depth_CF_fine,
    pop_1990_by_depth_F_fine, pop_1990_by_depth_CF_fine,
    ra_exposed_pop_2019_F, ra_exposed_pop_2019_CF,
    ra_exposed_pop_1990_F, ra_exposed_pop_1990_CF,
)
plt.show()

# %% ... add more cells as needed
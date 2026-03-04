"""
All figure-generating code for the socio-economic attribution analysis.

Expects all variables from run_population_analysis.py to be available
(run that script first, or use `from run_population_analysis import *`
in an interactive session).

Usage:
    # In an interactive session after running run_population_analysis.py:
    exec(open("plot_attribution_figures.py").read())
    # Or run cell-by-cell in VS Code / Jupyter.
"""
#%%
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as path_effects
import cartopy.crs as ccrs
import rasterio
from rasterio import features
from rasterio.transform import rowcol
from matplotlib.cm import ScalarMappable
from matplotlib.colors import (
    PowerNorm, TwoSlopeNorm, BoundaryNorm, Normalize, LinearSegmentedColormap,
)
from scipy.signal import find_peaks
from shapely.geometry import box

from config import COLOURS, POPULATION_2019_PATH
from population_utils import compute_attr_per_flood_depth_mask
from plot_helpers import (
    setup_map_axes,
    setup_plain_axes,
    annotate_landmarks,
    add_flood_depth_bg,
    make_discrete_pop_cmap,
    make_discrete_diff_cmap,
)

# NOTE: This script assumes the following variables are in scope from
# run_population_analysis.py.  If running standalone, import them first.
# Key variables needed:
#   pop_arrays, pop_array_uniform_2019,
#   ra_exposed_pop_2019_F, ra_exposed_pop_2019_CF,
#   ra_exposed_pop_1990_F, ra_exposed_pop_1990_CF,
#   ra_exposed_pop_2019_F_uniform,
#   gdf_pop_2019_exposed_F, gdf_pop_2019_exposed_CF,
#   gdf_pop_1990_exposed_F, gdf_pop_1990_exposed_CF,
#   gdf_pop_2019_exposed_F_coarse, gdf_pop_2019_exposed_CF_coarse,
#   gdf_pop_1990_exposed_F_coarse, gdf_pop_1990_exposed_CF_coarse,
#   gdf_pop_2019_exposed_F_uniform_coarse,
#   districts_adm3_filtered, districts_adm2,
#   hmax_F, hmax_CF, hmax_diff
#   flood_grid_transform, flood_grid_crs, flood_grid_shape, flood_extent,
#   region, region_utm, background, background_utm, bg_filtered_utm, beira_utm,
#   bin_centers, bin_centers_coarse, bins_fine, bins_coarse,
#   low_mask, mid_mask, high_mask,
#   pop_2019_by_depth_F_fine, pop_2019_by_depth_CF_fine,
#   pop_1990_by_depth_F_fine, pop_1990_by_depth_CF_fine,
#   pop_2019_by_depth_F_coarse, pop_2019_by_depth_CF_coarse,
#   pop_1990_by_depth_F_coarse, pop_1990_by_depth_CF_coarse,


colours = COLOURS


# ===================================================================== #
# FIG: Attributable % exposed population per flood depth (coarse bins)   #
# ===================================================================== #
#%%
def plot_attributable_pct_per_depth(bin_centers_coarse, pop_2019_by_depth_F_coarse,
                                    pop_2019_by_depth_CF_coarse, pop_1990_by_depth_F_coarse,
                                    pop_1990_by_depth_CF_coarse):
    Change = pd.DataFrame({
        "Factual": pop_2019_by_depth_F_coarse,
        "CF_climate": pop_2019_by_depth_CF_coarse,
        "CF_population": pop_1990_by_depth_F_coarse,
        "CF_climate_population": pop_1990_by_depth_CF_coarse,
    })
    for col, cf in [("Rel_change_CF_climate", "CF_climate"),
                    ("Rel_change_CF_population", "CF_population"),
                    ("Rel_change_CF_climate_population", "CF_climate_population")]:
        Change[col] = (Change["Factual"] - Change[cf]) / Change["Factual"] * 100

    fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=300)
    ax.plot(bin_centers_coarse, Change["Rel_change_CF_climate"], label="Climate change", color=colours[1])
    ax.plot(bin_centers_coarse, Change["Rel_change_CF_population"], label="Population change", color=colours[2])
    ax.plot(bin_centers_coarse, Change["Rel_change_CF_climate_population"], label="Climate & Population change", color=colours[3])

    ymax = max(Change["Rel_change_CF_climate"].max(), Change["Rel_change_CF_population"].max(),
               Change["Rel_change_CF_climate_population"].max()) * 1.05
    ymin = min(Change["Rel_change_CF_climate"].min(), Change["Rel_change_CF_population"].min(),
               Change["Rel_change_CF_climate_population"].min()) * 1.05

    add_flood_depth_bg(ax, ymin, ymax, label_y_frac=0.1)
    ax.axhline(y=0, linestyle="--", color="black", alpha=0.7)
    ax.set_xlabel("Flood depth (m)")
    ax.set_ylabel("Attributable exposed population (%)")
    ax.set_title("Attributable population exposed by flood depth bins")
    ax.legend(loc="upper right")
    ax.set_xlim(0, 3.5)
    ax.set_ylim(ymin, ymax)
    ax.grid(True, linestyle="--", alpha=0.5)
    return fig


# ===================================================================== #
# FIG: Exposed population per flood depth (fine bins, absolute)          #
# ===================================================================== #
#%%
def plot_exposed_pop_per_depth(bin_centers, pop_2019_by_depth_F_fine, pop_2019_by_depth_CF_fine,
                               pop_1990_by_depth_F_fine, pop_1990_by_depth_CF_fine,
                               ra_exposed_pop_2019_F, ra_exposed_pop_2019_CF,
                               ra_exposed_pop_1990_F, ra_exposed_pop_1990_CF):
    fig, ax = plt.subplots(figsize=(8, 5))
    x = bin_centers
    y_F = pop_2019_by_depth_F_fine.values
    y_CF_clim = pop_2019_by_depth_CF_fine.values
    y_CF_pop = pop_1990_by_depth_F_fine.values
    y_CF_clim_pop = pop_1990_by_depth_CF_fine.values

    ymax = max(y_F.max(), y_CF_clim.max(), y_CF_pop.max(), y_CF_clim_pop.max()) * 1.05
    add_flood_depth_bg(ax, 0, ymax)

    ax.plot(x, y_F, label=f"Factual ({np.nansum(ra_exposed_pop_2019_F).astype(int):,.0f} people)",
            color=colours[0], linewidth=2)
    ax.plot(x, y_CF_clim, label=f"CF Climate ({np.nansum(ra_exposed_pop_2019_CF).astype(int):,.0f} people)",
            color=colours[1], linewidth=1)
    ax.plot(x, y_CF_pop, label=f"CF Population ({np.nansum(ra_exposed_pop_1990_F).astype(int):,.0f} people)",
            color=colours[2], linewidth=1)
    ax.plot(x, y_CF_clim_pop, label=f"CF Both ({np.nansum(ra_exposed_pop_1990_CF).astype(int):,.0f} people)",
            color=colours[3], linewidth=1)

    # Peak annotations
    for y_arr, c in [(y_F, colours[0]), (y_CF_clim, colours[1]), (y_CF_pop, colours[2])]:
        peaks, props = find_peaks(y_arr, prominence=0.02, distance=5)
        sorted_peaks = peaks[np.argsort(props["prominences"])[::-1]]
        for idx in sorted_peaks[:2]:
            plt.annotate(
                f"{x[idx]:.2f} m", xy=(x[idx] + 0.02, y_arr[idx]),
                xytext=(x[idx] + 0.15, y_arr[idx]), color=c, fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.8, color=c, linestyle="--"),
            )

    ax.set_xlabel("Flood depth (m)")
    ax.set_ylabel("Exposed population")
    ax.set_xlim(0.05, 3.5)
    ax.set_ylim(0, ymax)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend()
    return fig


# ===================================================================== #
# FIG: Absolute change compared to factual (diff curves)                 #
# ===================================================================== #
#%%
def plot_absolute_diff_per_depth(bin_centers, pop_2019_by_depth_F_fine, pop_2019_by_depth_CF_fine,
                                  pop_1990_by_depth_F_fine, pop_1990_by_depth_CF_fine,
                                  ra_exposed_pop_2019_F, ra_exposed_pop_2019_CF,
                                  ra_exposed_pop_1990_F, ra_exposed_pop_1990_CF):
    perct_attr_clim = (np.nansum(ra_exposed_pop_2019_F) - np.nansum(ra_exposed_pop_2019_CF)) / np.nansum(ra_exposed_pop_2019_F) * 100
    perct_attr_pop = (np.nansum(ra_exposed_pop_2019_F) - np.nansum(ra_exposed_pop_1990_F)) / np.nansum(ra_exposed_pop_2019_F) * 100
    perct_attr_clim_pop = (np.nansum(ra_exposed_pop_2019_F) - np.nansum(ra_exposed_pop_1990_CF)) / np.nansum(ra_exposed_pop_2019_F) * 100

    diff_clim = pop_2019_by_depth_F_fine.values - pop_2019_by_depth_CF_fine.values
    diff_pop = pop_2019_by_depth_F_fine.values - pop_1990_by_depth_F_fine.values
    diff_clim_pop = pop_2019_by_depth_F_fine.values - pop_1990_by_depth_CF_fine.values

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(bin_centers, diff_clim, label=f"Climate change ({int(np.round(perct_attr_clim))} %)", color=colours[1])
    ax.plot(bin_centers, diff_pop, label=f"Population change ({int(np.round(perct_attr_pop))} %)", color=colours[2])
    ax.plot(bin_centers, diff_clim_pop, label=f"Climate & Population change ({int(np.round(perct_attr_clim_pop))} %)", color=colours[3])

    ymin_val = diff_clim.min() * 1.1
    ymax_val = diff_clim_pop.max() * 1.1
    add_flood_depth_bg(ax, ymin_val, ymax_val)

    ax.set_xlabel("Flood depth (m)")
    ax.set_ylabel("Absolute change in exposed population")
    ax.axhline(y=0, linestyle="--", color="black", alpha=0.7)
    ax.legend(loc="upper right", fontsize=9)
    ax.set_xlim(0.05, 3.5)
    ax.set_ylim(ymin_val, ymax_val)
    ax.grid(True, linestyle="--", alpha=0.5)

    return fig, diff_clim, diff_pop, diff_clim_pop, perct_attr_clim, perct_attr_pop, perct_attr_clim_pop


# ===================================================================== #
# FIG: Bar plot — absolute & attributable change per depth category       #
# ===================================================================== #
#%%
def plot_bar_attribution(diff_clim, diff_pop, diff_clim_pop, low_mask, mid_mask, high_mask,
                          pop_2019_by_depth_F_fine, perct_attr_clim, perct_attr_pop, perct_attr_clim_pop):
    low_abs = [diff_clim[low_mask].sum(), diff_pop[low_mask].sum(), diff_clim_pop[low_mask].sum()]
    mid_abs = [diff_clim[mid_mask].sum(), diff_pop[mid_mask].sum(), diff_clim_pop[mid_mask].sum()]
    high_abs = [diff_clim[high_mask].sum(), diff_pop[high_mask].sum(), diff_clim_pop[high_mask].sum()]

    low_attr = [compute_attr_per_flood_depth_mask(d, pop_2019_by_depth_F_fine, low_mask)
                for d in [diff_clim, diff_pop, diff_clim_pop]]
    mid_attr = [compute_attr_per_flood_depth_mask(d, pop_2019_by_depth_F_fine, mid_mask)
                for d in [diff_clim, diff_pop, diff_clim_pop]]
    high_attr = [compute_attr_per_flood_depth_mask(d, pop_2019_by_depth_F_fine, high_mask)
                 for d in [diff_clim, diff_pop, diff_clim_pop]]

    data_abs = np.array([low_abs, mid_abs, high_abs])
    data_attr = np.array([low_attr, mid_attr, high_attr])

    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharex=True)
    bar_width = 0.25
    x_pos = np.arange(3)
    labels_x = ["Low\n(<0.5 m)", "Medium\n(0.5–1.5 m)", "High\n(>1.5 m)"]
    subplot_labels = ["(a)", "(b)"]

    for j, (ax, data, ylabel) in enumerate(zip(
            axes, [data_abs, data_attr],
            ["Absolute change in exposed population (people)", "Attributable exposed population (%)"])):
        ax.bar(x_pos - bar_width, data[:, 0], width=bar_width,
               label=f"Climate change ({int(np.round(perct_attr_clim))} %)", color=colours[1])
        ax.bar(x_pos, data[:, 1], width=bar_width,
               label=f"Population change ({int(np.round(perct_attr_pop))} %)", color=colours[2])
        ax.bar(x_pos + bar_width, data[:, 2], width=bar_width,
               label=f"Climate & population ({int(np.round(perct_attr_clim_pop))} %)", color=colours[3])
        ax.axhline(0, linestyle="--", color="black", alpha=0.7)
        ax.set_axisbelow(True)
        ax.grid(True, axis="y", linestyle="--", alpha=0.5)
        ax.set_xlabel("Flood depth", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels_x, fontdict={"fontweight": "bold", "color": "#5C5C5C"})
        ax.text(0, 1.02, subplot_labels[j], transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="bottom", ha="left")

    axes[1].legend(fontsize=9, loc="upper right", bbox_to_anchor=(1, 1.2))
    plt.tight_layout()
    return fig


# ===================================================================== #
# FIG 2: Attributable exposed population maps (3 drivers, absolute)      #
# ===================================================================== #
#%%
def plot_fig2_attribution_maps(gdf_pop_2019_exposed_F_coarse, gdf_pop_2019_exposed_CF_coarse,
                                gdf_pop_1990_exposed_F_coarse, gdf_pop_1990_exposed_CF_coarse,
                                region_utm, background_utm, beira_utm, flood_extent):
    print("Plotting Fig 2: attributable exposed population (three drivers)")
    gdf_CF_pop = gdf_pop_1990_exposed_F_coarse.copy()
    gdf_CF_clim = gdf_pop_2019_exposed_CF_coarse.copy()
    gdf_CF_clim_pop = gdf_pop_1990_exposed_CF_coarse.copy()

    total_factual = gdf_pop_2019_exposed_F_coarse["exposed_population"].sum()
    gdf_CF_clim["diff"] = gdf_pop_2019_exposed_F_coarse["exposed_population"] - gdf_CF_clim["exposed_population"]
    gdf_CF_pop["diff"] = gdf_pop_2019_exposed_F_coarse["exposed_population"] - gdf_CF_pop["exposed_population"]
    gdf_CF_clim_pop["diff"] = gdf_pop_2019_exposed_F_coarse["exposed_population"] - gdf_CF_clim_pop["exposed_population"]

    datasets = [
        ("Climate change", gdf_CF_clim),
        ("Population change", gdf_CF_pop),
        ("Climate & population change", gdf_CF_clim_pop),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(11, 5), dpi=300, constrained_layout=True,
                             subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

    vmax = max(ds["diff"].max() for _, ds in datasets)
    norm_diff = PowerNorm(gamma=0.5, vmin=0, vmax=vmax)
    subplot_labels = ["(a)", "(b)", "(c)"]

    for i, (title, gdf) in enumerate(datasets):
        ax = axes[i]
        gdf[gdf["diff"] <= 0].plot(ax=ax, color="white", edgecolor="grey", linewidth=0.2, zorder=1)
        gdf[gdf["diff"] > 0].plot(column="diff", cmap=plt.cm.Reds, norm=norm_diff, edgecolor="grey",
                                  linewidth=0.2, ax=ax, legend=False, zorder=2, rasterized=True)

        setup_map_axes(ax, region_utm, background_utm, flood_extent,
                       subplot_labels=[subplot_labels[i]], titles=[title])
        beira_utm.boundary.plot(ax=ax, edgecolor="black", linewidth=0.5, zorder=3, alpha=0.7)
        annotate_landmarks(ax)

        rel_change = gdf["diff"].sum() / total_factual * 100
        ax.text(0.98, 0.98, f"~{round(gdf['diff'].sum(), -3):,.0f} people\n({rel_change:.0f} %)",
                transform=ax.transAxes, ha="right", va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.8))

    sm = ScalarMappable(cmap=plt.cm.Reds, norm=norm_diff)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation="vertical", shrink=0.6, pad=0.02)
    cbar.set_label("Attributable exposed population [# people]", fontsize=10)
    return fig


# ===================================================================== #
# FIG 5: Attributable exposed pop per district                           #
# ===================================================================== #
#%%
def plot_fig5_district_attribution(districts_adm3_clipped, region_utm, background_utm, flood_extent):
    columns = ["attr_clim", "attr_pop", "attr_clim_pop"]
    titles = ["Climate change", "Population change", "Climate and population change"]
    subplot_labels = ["(a)", "(b)", "(c)"]

    all_attr = np.concatenate([districts_adm3_clipped[c].values for c in columns])
    all_attr = all_attr[~np.isnan(all_attr)]
    vabs = np.nanpercentile(np.abs(all_attr), 99)
    norm_attr = PowerNorm(gamma=0.5, vmin=0, vmax=vabs)
    cmap_attr = LinearSegmentedColormap.from_list("reds_from_rdbu", plt.cm.RdBu_r(np.linspace(0.5, 1, 256)))

    fig, axes = plt.subplots(1, 3, figsize=(11, 6), dpi=300, sharey=True, constrained_layout=True,
                             subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

    for i, (ax, col, title) in enumerate(zip(axes, columns, titles)):
        districts_adm3_clipped.plot(column=col, cmap=cmap_attr, norm=norm_attr, edgecolor="grey",
                                    linewidth=0.2, ax=ax, legend=False, rasterized=True,
                                    missing_kwds={"color": "white"})
        setup_map_axes(ax, region_utm, background_utm, flood_extent,
                       subplot_labels=[subplot_labels[i]], titles=[title])

        for _, row in districts_adm3_clipped.iterrows():
            x, y = row.geometry.centroid.x, row.geometry.centroid.y
            value = row[col]
            name = row["NAME_3"]
            ax.text(x, y, f"{name}\n{int(round(value))} %", ha="center", va="center", fontsize=7,
                    bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=0.3, boxstyle="round"))

    sm = ScalarMappable(cmap=cmap_attr, norm=norm_attr)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, shrink=0.6)
    cbar.set_label("Attributable relative exposed population [%]", fontsize=9)
    return fig


# ===================================================================== #
# FIG 3: Relative exposed pop — factual + diffs for clim & pop          #
# ===================================================================== #
#%%
def plot_fig3_relative_exposure(gdf_pop_2019_exposed_F_coarse, gdf_pop_2019_exposed_CF_coarse,
                                 gdf_pop_1990_exposed_F_coarse, region_utm, background,
                                 flood_extent):
    print("Plotting Fig 3: relative exposed population")

    gdf_pop_2019_exposed_CF_coarse["rel_diff"] = (
        gdf_pop_2019_exposed_F_coarse["relative_population"]
        - gdf_pop_2019_exposed_CF_coarse["relative_population"]
    )
    gdf_pop_1990_exposed_F_coarse["rel_diff"] = (
        gdf_pop_2019_exposed_F_coarse["relative_population"]
        - gdf_pop_1990_exposed_F_coarse["relative_population"]
    )

    fig, axes = plt.subplots(1, 3, figsize=(10, 5), dpi=300, sharey=True, constrained_layout=True,
                             subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

    norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2019_exposed_F_coarse["relative_population"].max())

    vmax_diff_clim = gdf_pop_2019_exposed_CF_coarse["rel_diff"].quantile(0.99)
    vmin_diff_clim = gdf_pop_2019_exposed_CF_coarse["rel_diff"].quantile(0.01)
    norm_diff_clim = PowerNorm(gamma=0.5, vmin=vmin_diff_clim, vmax=vmax_diff_clim)
    red_half = LinearSegmentedColormap.from_list("bwr_red", plt.cm.bwr(np.linspace(0.5, 1, 256)))

    vmax_diff_pop = gdf_pop_1990_exposed_F_coarse["rel_diff"].quantile(0.99)
    vmin_diff_pop = gdf_pop_1990_exposed_F_coarse["rel_diff"].quantile(0.01)
    norm_diff_pop = TwoSlopeNorm(vmin=vmin_diff_pop, vcenter=0, vmax=vmax_diff_pop)

    # Factual relative
    gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse["relative_population"] == 0].plot(
        ax=axes[0], color="white", edgecolor="grey", linewidth=0.2, zorder=1)
    gdf_pop_2019_exposed_F_coarse[gdf_pop_2019_exposed_F_coarse["relative_population"] > 0].plot(
        column="relative_population", cmap="Blues", edgecolor="grey", norm=norm_pop,
        linewidth=0.2, ax=axes[0], legend=False, zorder=2, rasterized=True)

    # Climate diff
    gdf_pop_2019_exposed_CF_coarse.plot(
        column="rel_diff", cmap=red_half, norm=norm_diff_clim, edgecolor="grey",
        linewidth=0.2, ax=axes[1], legend=False, zorder=2,
        missing_kwds={"color": "white", "edgecolor": "none"}, rasterized=True)

    # Population diff
    gdf_pop_1990_exposed_F_coarse.plot(
        column="rel_diff", cmap="bwr", norm=norm_diff_pop, edgecolor="grey",
        linewidth=0.2, ax=axes[2], legend=False, zorder=2,
        missing_kwds={"color": "white", "edgecolor": "none"}, rasterized=True)

    subplot_labels = ["(a)", "(b)", "(c)"]
    titles = ["Factual", "Factual - CF climate", "Factual - CF population"]

    # Background with coastline fix
    mask_box = box(34.8, -20.3, 35.3, -19.9)
    bg_outside = background[~background.intersects(mask_box)]

    for i, ax in enumerate(axes):
        region_utm.boundary.plot(ax=ax, edgecolor="black", linewidth=0.3)
        ax.set_extent(flood_extent, crs=ccrs.UTM(36, southern_hemisphere=True))
        bg_outside.plot(ax=ax, color="#E0E0E0", transform=ccrs.PlateCarree(), zorder=0)
        bg_outside.boundary.plot(ax=ax, color="#818181", linewidth=0.2,
                                 transform=ccrs.PlateCarree(), zorder=1)
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color="gray", alpha=0.5, linestyle="--")
        gl.right_labels = gl.top_labels = False
        gl.xlabel_style = gl.ylabel_style = {"size": 9}
        if i != 0:
            gl.left_labels = False
        ax.text(-0.05, 1.02, subplot_labels[i], transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="bottom", ha="left")
        ax.set_title(titles[i], fontsize=8)

    # Colorbars
    sm1 = ScalarMappable(cmap="Blues", norm=norm_pop)
    sm1._A = []
    fig.colorbar(sm1, ax=axes[0], orientation="vertical", shrink=0.5).set_label(
        "Aggregated relative exposed population [%]", fontsize=10)

    sm2 = ScalarMappable(cmap=red_half, norm=norm_diff_clim)
    sm2.set_array([])
    fig.colorbar(sm2, ax=axes[1], orientation="vertical", shrink=0.5).set_label(
        "Attributable relative exposed population [%]", fontsize=9)

    sm3 = ScalarMappable(cmap="bwr", norm=norm_diff_pop)
    sm3.set_array([])
    fig.colorbar(sm3, ax=axes[2], orientation="vertical", shrink=0.5).set_label(
        "Attributable relative exposed population [%]", fontsize=9)

    return fig


# ===================================================================== #
# SUPPLEMENTARY: Flood depth among exposed population                    #
# ===================================================================== #
#%%
def plot_supp_flood_depth_exposed(gdf_pop_2019_exposed_F_coarse, gdf_pop_2019_exposed_CF_coarse,
                                   region_utm, background_utm, bg_filtered_utm, flood_extent):
    gdf_F = gdf_pop_2019_exposed_F_coarse.copy()
    gdf_F.loc[gdf_F["exposed_population"] <= 0, "avg_flood_depth"] = np.nan
    gdf_CF = gdf_pop_2019_exposed_CF_coarse.copy()
    gdf_CF.loc[gdf_CF["exposed_population"] <= 0, "avg_flood_depth"] = np.nan
    gdf_CF["change_in_flood_depth"] = gdf_F["avg_flood_depth"] - gdf_CF["avg_flood_depth"]

    cmap_depth = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#67CBE4"])
    cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_purple", ["#ffffff", "#651F94"])

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    gdf_F.plot(column="avg_flood_depth", cmap=cmap_depth, vmin=0, vmax=3.5, linewidth=0.1,
               edgecolor="grey", ax=axes[0], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    gdf_CF.plot(column="avg_flood_depth", cmap=cmap_depth, vmin=0, vmax=3.5, linewidth=0.1,
                edgecolor="grey", ax=axes[1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    gdf_CF.plot(column="change_in_flood_depth", cmap=cmap_change, vmin=0, vmax=0.5, linewidth=0.1,
                edgecolor="grey", ax=axes[2], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

    setup_plain_axes(axes, region_utm, background_utm, bg_filtered_utm, flood_extent)

    for ax in axes[:2]:
        sm = ScalarMappable(cmap=cmap_depth, norm=plt.Normalize(vmin=0, vmax=3.5))
        sm._A = []
        plt.colorbar(sm, ax=ax, shrink=0.8).set_label("Average flood depth among exposed population (m)")

    sm = ScalarMappable(cmap=cmap_change, norm=plt.Normalize(vmin=0, vmax=0.5))
    sm._A = []
    plt.colorbar(sm, ax=axes[2], shrink=0.8).set_label("Difference in average flood depth (m)")

    axes[0].set_title("Factual")
    axes[1].set_title("No Climate Change")
    axes[2].set_title("Factual - Counterfactual")
    fig.suptitle("Average flood depth among exposed population", fontsize=12)
    plt.tight_layout()
    return fig


# ===================================================================== #
# SUPPLEMENTARY: Population change map (1990 vs 2019)                    #
# ===================================================================== #
#%%
def plot_supp_pop_change_map(gdf_pop_2019_exposed_F_coarse, gdf_pop_1990_exposed_F_coarse,
                              region_utm, background_utm, bg_filtered_utm, flood_extent):
    gdf_2019 = gdf_pop_2019_exposed_F_coarse.copy()
    gdf_2019.loc[gdf_2019["total_population"] == 0] = np.nan
    gdf_1990 = gdf_pop_1990_exposed_F_coarse.copy()
    gdf_1990.loc[gdf_1990["total_population"] == 0] = np.nan
    gdf_1990["change_in_population"] = gdf_2019["total_population"] - gdf_1990["total_population"]

    cmap = mcolors.LinearSegmentedColormap.from_list("white_to_orange", ["#ffffff", "#FC6F37"])
    norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2019["total_population"]))
    cmap_change = mcolors.LinearSegmentedColormap.from_list("white_to_red", ["#ffffff", "#BD2A2A"])
    norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_1990["change_in_population"]))

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    gdf_2019.plot(column="total_population", cmap=cmap, norm=norm, linewidth=0.1, edgecolor="grey",
                  ax=axes[0], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    gdf_1990.plot(column="total_population", cmap=cmap, norm=norm, linewidth=0.1, edgecolor="grey",
                  ax=axes[1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    gdf_1990.plot(column="change_in_population", cmap=cmap_change, norm=norm_change, linewidth=0.1,
                  edgecolor="grey", ax=axes[2], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

    setup_plain_axes(axes, region_utm, background_utm, bg_filtered_utm, flood_extent)

    for ax in axes[:2]:
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm._A = []
        plt.colorbar(sm, ax=ax, shrink=0.8).set_label("Population (people per cell)")

    sm = ScalarMappable(cmap=cmap_change, norm=norm_change)
    sm._A = []
    plt.colorbar(sm, ax=axes[2], shrink=0.8).set_label("Difference in population")

    axes[0].set_title("2019 Population")
    axes[1].set_title("1990 Population")
    axes[2].set_title("2019 − 1990")
    fig.suptitle("Total population in study region", fontsize=12)
    plt.tight_layout()
    return fig


# ===================================================================== #
# SUPPLEMENTARY: Uniform vs spatial pop growth (coarse + Beira zoom)     #
# ===================================================================== #
#%%
def plot_supp_uniform_vs_spatial(gdf_pop_2019_exposed_F_coarse, gdf_pop_2019_exposed_F_uniform_coarse,
                                  region_utm, background_utm, flood_extent):
    gdf_F = gdf_pop_2019_exposed_F_coarse.copy()
    gdf_U = gdf_pop_2019_exposed_F_uniform_coarse.copy()
    gdf_U["population_diff"] = gdf_F["exposed_population"] - gdf_U["exposed_population"]

    fig, axes = plt.subplots(1, 3, figsize=(10, 6), dpi=300, sharey=True, constrained_layout=True,
                             subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

    norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_F["exposed_population"].max())
    norm_diff = mcolors.TwoSlopeNorm(vmin=gdf_U["population_diff"].min(), vcenter=0,
                                     vmax=gdf_U["population_diff"].max())

    gdf_F[gdf_F["exposed_population"] > 0].plot(column="exposed_population", cmap="Blues", edgecolor="grey",
                                                  norm=norm_pop, linewidth=0.2, ax=axes[0], legend=False,
                                                  zorder=2, rasterized=True)
    gdf_U[gdf_U["exposed_population"] > 0].plot(column="exposed_population", cmap="Blues", edgecolor="grey",
                                                  norm=norm_pop, linewidth=0.2, ax=axes[1], legend=False,
                                                  zorder=2, rasterized=True)
    gdf_U.plot(column="population_diff", cmap="RdBu_r", norm=norm_diff, edgecolor="grey",
               linewidth=0.2, ax=axes[2], missing_kwds={"color": "white", "edgecolor": "none"},
               zorder=2, rasterized=True)

    setup_map_axes(axes, region_utm, background_utm, flood_extent,
                   subplot_labels=["(a)", "(b)", "(c)"],
                   titles=["Factual (2019)", "Uniform (2019)", "Spatial − Uniform"])

    sm1 = ScalarMappable(cmap="Blues", norm=norm_pop)
    sm1._A = []
    fig.colorbar(sm1, ax=axes[1], shrink=0.5).set_label("Exposed population (# people)")

    sm2 = ScalarMappable(cmap="RdBu_r", norm=norm_diff)
    sm2.set_array([])
    fig.colorbar(sm2, ax=axes[2], shrink=0.5).set_label("Attributable exposed population [# people]")
    return fig


# ===================================================================== #
# SUPPLEMENTARY: Beira zoom — 25 m resolution rasters                   #
# ===================================================================== #
#%%
def plot_supp_beira_zoom(pop_arrays, pop_array_uniform_2019, ra_exposed_pop_2019_F,
                          ra_exposed_pop_2019_F_uniform, flood_grid_transform,
                          background_utm, bg_filtered_utm, region_utm):
    xmin, xmax = 690000, 702000
    ymin, ymax = 7803000, 7818000

    row_ul, col_ul = rowcol(flood_grid_transform, xmin, ymax, op=int)
    row_lr, col_lr = rowcol(flood_grid_transform, xmax, ymin, op=int)

    row_min = max(0, min(row_ul, row_lr))
    row_max = min(ra_exposed_pop_2019_F.shape[0], max(row_ul, row_lr))
    col_min = max(0, min(col_ul, col_lr))
    col_max = min(ra_exposed_pop_2019_F.shape[1], max(col_ul, col_lr))

    x_min_clip, y_max_clip = flood_grid_transform * (col_min, row_min)
    x_max_clip, y_min_clip = flood_grid_transform * (col_max, row_max)
    extent_beira = [x_min_clip, x_max_clip, y_max_clip, y_min_clip]

    raster_F = pop_arrays[2019][row_min:row_max, col_min:col_max]
    raster_UF = pop_array_uniform_2019[row_min:row_max, col_min:col_max]
    diff_pop = raster_UF - raster_F

    raster_F_aff = ra_exposed_pop_2019_F[row_min:row_max, col_min:col_max]
    raster_UF_aff = ra_exposed_pop_2019_F_uniform[row_min:row_max, col_min:col_max]
    diff_aff = raster_UF_aff - raster_F_aff

    pop_cmap, pop_norm, _ = make_discrete_pop_cmap(
        mcolors.LinearSegmentedColormap.from_list("w2p", ["#ffffff", "#651F94"]),
        np.nanmax(raster_F))
    aff_cmap, aff_norm, _ = make_discrete_pop_cmap(plt.cm.Blues, np.nanmax(raster_F_aff))
    diff_cmap_pop, diff_norm_pop, _ = make_discrete_diff_cmap(np.nanmax(np.abs(diff_pop)))
    diff_cmap_aff, diff_norm_aff, _ = make_discrete_diff_cmap(np.nanmax(np.abs(diff_aff)))

    fig, axes = plt.subplots(2, 3, figsize=(16, 12), sharex=True, sharey=True, dpi=300, constrained_layout=True)

    for ax in axes.flatten():
        background_utm.plot(ax=ax, color="#E0E0E0", zorder=0)
        bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
        region_utm.boundary.plot(ax=ax, color="black", linewidth=1, zorder=2)

    # Top row: total population
    axes[0, 0].imshow(raster_F, cmap=pop_cmap, norm=pop_norm, extent=extent_beira, origin="lower", alpha=0.8, zorder=3)
    axes[0, 0].set_title("Factual 2019 Population")
    im = axes[0, 1].imshow(raster_UF, cmap=pop_cmap, norm=pop_norm, extent=extent_beira, origin="lower", alpha=0.8, zorder=3)
    plt.colorbar(im, ax=axes[0, 1], shrink=0.8).set_label("Population")
    axes[0, 1].set_title("Uniform 2019 Population")
    im = axes[0, 2].imshow(diff_pop, cmap=diff_cmap_pop, norm=diff_norm_pop, extent=extent_beira, origin="lower", alpha=0.8, zorder=3)
    plt.colorbar(im, ax=axes[0, 2], shrink=0.8).set_label("Difference")
    axes[0, 2].set_title("Uniform − Factual")

    # Bottom row: affected population
    axes[1, 0].imshow(raster_F_aff, cmap=aff_cmap, norm=aff_norm, extent=extent_beira, origin="lower", alpha=0.8, zorder=3)
    axes[1, 0].set_title("Factual affected")
    im = axes[1, 1].imshow(raster_UF_aff, cmap=aff_cmap, norm=aff_norm, extent=extent_beira, origin="lower", alpha=0.8, zorder=3)
    plt.colorbar(im, ax=axes[1, 1], shrink=0.8).set_label("Affected population")
    axes[1, 1].set_title("Uniform affected")
    im = axes[1, 2].imshow(diff_aff, cmap=diff_cmap_aff, norm=diff_norm_aff, extent=extent_beira, origin="lower", alpha=0.8, zorder=3)
    plt.colorbar(im, ax=axes[1, 2], shrink=0.8).set_label("Difference")
    axes[1, 2].set_title("Uniform − Factual affected")

    fig.suptitle("Uniform vs. spatially differing population growth (Beira zoom)", fontsize=14)
    return fig


# ===================================================================== #
# SUPPLEMENTARY: Relative exposure diffs (3-panel maps)                  #
# ===================================================================== #
#%%
def plot_supp_relative_diff_maps(gdf_pop_2019_exposed_F_coarse, gdf_pop_2019_exposed_CF_coarse,
                                  gdf_pop_1990_exposed_F_coarse, region_utm, background_utm,
                                  bg_filtered_utm, flood_extent):
    gdf_CF_clim = gdf_pop_2019_exposed_CF_coarse.copy()
    gdf_CF_pop = gdf_pop_1990_exposed_F_coarse.copy()

    gdf_CF_clim["population_diff"] = (gdf_pop_2019_exposed_F_coarse["relative_population"]
                                      - gdf_CF_clim["relative_population"])
    gdf_CF_pop["population_diff"] = (gdf_pop_2019_exposed_F_coarse["relative_population"]
                                     - gdf_CF_pop["relative_population"])

    norm_diff = TwoSlopeNorm(vmin=-50, vcenter=0, vmax=50)
    cmap_diff = plt.get_cmap("RdBu_r")

    fig, axes = plt.subplots(1, 2, figsize=(10, 6), dpi=300, sharey=True, constrained_layout=True,
                             subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

    gdf_CF_clim.plot(column="population_diff", cmap=cmap_diff, norm=norm_diff, edgecolor="grey",
                     linewidth=0.2, ax=axes[0], legend=False, zorder=2, rasterized=True,
                     missing_kwds={"color": "white"})
    gdf_CF_pop.plot(column="population_diff", cmap=cmap_diff, norm=norm_diff, edgecolor="grey",
                    linewidth=0.2, ax=axes[1], legend=False, zorder=2, rasterized=True,
                    missing_kwds={"color": "white"})

    setup_map_axes(axes, region_utm, background_utm, flood_extent,
                   subplot_labels=["(a)", "(b)"],
                   titles=["Climate change", "Population change"])

    sm = ScalarMappable(cmap=cmap_diff, norm=norm_diff)
    sm.set_array([])
    fig.colorbar(sm, ax=axes, shrink=0.6).set_label("Relative exposed population change [%]")
    return fig


# ================================================================================== #
# SUPPLEMENTARY: Factual flood, population, exposed population and changes (6-panel) #
# ================================================================================== #
#%%
def plot_supp_factual_changes_overview(hmax_F, gdf_pop_2019_exposed_F_coarse, hmax_diff, 
                                       gdf_pop_1990_exposed_F_coarse, gdf_pop_1990_exposed_CF_coarse,
                                       region_utm, background_utm, flood_extent):
    # Data preparation for plotting
    gdf_2019 = gdf_pop_2019_exposed_F_coarse.copy()
    gdf_2019.loc[gdf_2019["total_population"] == 0] = np.nan
    gdf_1990 = gdf_pop_1990_exposed_F_coarse.copy()
    gdf_1990.loc[gdf_1990["total_population"] == 0] = np.nan
    gdf_2019["change_in_population"] = gdf_2019["total_population"] - gdf_1990["total_population"]
    gdf_attr = gdf_pop_2019_exposed_F_coarse.copy()
    gdf_attr["attr_exposed_pop"] = gdf_pop_2019_exposed_F_coarse["exposed_population"] - gdf_pop_1990_exposed_CF_coarse["exposed_population"]

    # colour maps and norms
    norm_pop = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_2019["total_population"].max())
    cmap_pop = mcolors.LinearSegmentedColormap.from_list("white_to_blue", ["#ffffff", "#FC6F37"])
    norm_pop_exposed = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_2019_exposed_F_coarse["exposed_population"].max())
    cmap_pop_exposed = mcolors.LinearSegmentedColormap.from_list("white_to_darkblue", ["#ffffff", "#67CBE4"])
    cmap_change = plt.cm.Reds
    norm_pop_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_2019["change_in_population"]))
    norm_pop_exposed_change = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_pop_1990_exposed_CF_coarse["attr_exposed_pop"].max())

    fig, axes = plt.subplots(2, 3, figsize=(16, 8), dpi=300, sharex=True, sharey=True, constrained_layout=True,
                             subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

    # Plot 1 - Factual fooding
    im_1 = axes[0,0].imshow(hmax_F, cmap="viridis", extent=flood_extent, origin="lower", vmin=0, vmax=3.5, zorder=2)

    # Plot 2 - Factual population
    gdf_2019.plot(column="total_population", cmap=cmap_pop, norm=norm_pop, linewidth=0.1,
                  edgecolor="grey", ax=axes[0,1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

    # Plot 3 - Factual exposed population
    gdf_exposed_F = gdf_pop_2019_exposed_F_coarse.copy()
    gdf_exposed_F.loc[gdf_exposed_F["exposed_population"] == 0, "exposed_population"] = np.nan
    gdf_exposed_F.plot(column="exposed_population", cmap=cmap_pop_exposed, edgecolor="grey", 
                       norm=norm_pop_exposed, linewidth=0.2, ax=axes[0,2], legend=False, 
                       zorder=2, rasterized=True,
                       missing_kwds={"color": "none", "edgecolor": "none"})
    
    # Plot 4 - Climate change
    im_2 = axes[1,0].imshow(hmax_diff, cmap=cmap_change, extent=flood_extent, origin="lower", vmin=0, vmax=0.5, zorder=2)

    # Plot 5 - Population change
    gdf_2019.plot(column="change_in_population", cmap=cmap_change, norm=norm_pop_change, linewidth=0.1,
                  edgecolor="grey", ax=axes[1,1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    
    # Plot 6 - Attributable exposed population
    gdf_attr.loc[gdf_attr["attr_exposed_pop"] <= 0, "attr_exposed_pop"] = np.nan
    gdf_attr.plot(column="attr_exposed_pop", cmap=cmap_change, edgecolor="grey", 
                  norm=norm_pop_exposed_change, linewidth=0.2, ax=axes[1,2], legend=False, 
                  zorder=2, rasterized=True,
                  missing_kwds={"color": "none", "edgecolor": "none"})

    setup_map_axes(axes, region_utm, background_utm, flood_extent,
                   subplot_labels=["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"],
                   titles=["Factual flooding", "Factual population", "Factual exposed population", 
                           "Climate change", "Population change", "Attributable exposed population"])

    # Colour bars top row
    plt.colorbar(im_1, ax=axes[0,0], shrink=0.8).set_label("Flood depth (m)")
    sm = ScalarMappable(cmap=cmap_pop, norm=norm_pop)
    sm._A = []
    fig.colorbar(sm, ax=axes[0,1], shrink=0.8).set_label("Aggregated population (# people)")
    sm = ScalarMappable(cmap=cmap_pop_exposed, norm=norm_pop_exposed)
    sm._A = []
    fig.colorbar(sm, ax=axes[0,2], shrink=0.8).set_label("Aggregated exposed population (# people)")

    # Colour bars for the bottom row
    plt.colorbar(im_2, ax=axes[1,0], shrink=0.8).set_label("Attributable flood depth (m)")
    sm = ScalarMappable(cmap=cmap_change, norm=norm_pop_change)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes[1,1], orientation="vertical", shrink=0.8)
    cbar.set_label("Change in population [# people]", fontsize=10)
    sm = ScalarMappable(cmap=cmap_change, norm=norm_pop_exposed_change)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes[1,2], orientation="vertical", shrink=0.8)
    cbar.set_label("Attributable exposed population [# people]", fontsize=10)

    return fig



# ===================================================================== #
# FIG: Dominant driver of attributable exposed population                #
# ===================================================================== #
#%%
def plot_dominant_driver(bin_centers_coarse, pop_2019_by_depth_F_coarse,
                         pop_2019_by_depth_CF_coarse, pop_1990_by_depth_F_coarse,
                         pop_1990_by_depth_CF_coarse):
    Change = pd.DataFrame({
        "Factual": pop_2019_by_depth_F_coarse,
        "CF_climate": pop_2019_by_depth_CF_coarse,
        "CF_population": pop_1990_by_depth_F_coarse,
        "CF_climate_population": pop_1990_by_depth_CF_coarse,
    })
    Change["Rel_change_CF_climate"] = (Change["Factual"] - Change["CF_climate"]) / Change["Factual"] * 100
    Change["Rel_change_CF_population"] = (Change["Factual"] - Change["CF_population"]) / Change["Factual"] * 100
    Change["Dominant_driver_cc"] = Change["Rel_change_CF_climate"] / Change["Rel_change_CF_population"]
    Change["Dominant_driver_pc"] = Change["Rel_change_CF_population"] / Change["Rel_change_CF_climate"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), dpi=300)

    axes[0].plot(bin_centers_coarse, Change["Dominant_driver_cc"], color=colours[1])
    axes[0].axhline(y=0, linestyle="--", color="black", alpha=0.7)
    axes[0].set_xlabel("Flood depth (m)")
    axes[0].set_ylabel("Climate change / Population change")
    axes[0].set_title("Dominant driver (CC / PC)")
    axes[0].set_xlim(0, 3.5)
    axes[0].grid(True, linestyle="--", alpha=0.5)

    axes[1].plot(bin_centers_coarse, Change["Dominant_driver_pc"], color=colours[2])
    axes[1].axhline(y=0, linestyle="--", color="black", alpha=0.7)
    axes[1].set_xlabel("Flood depth (m)")
    axes[1].set_ylabel("Population change / Climate change")
    axes[1].set_title("Dominant driver (PC / CC)")
    axes[1].set_xlim(0, 3.5)
    axes[1].grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    return fig


# ===================================================================== #
# FIG 2b: Percentage attributable exposure maps (3-panel)                #
# ===================================================================== #
#%%
def plot_fig2b_pct_attribution_maps(gdf_pop_2019_exposed_F_coarse, gdf_pop_2019_exposed_CF_coarse,
                                     gdf_pop_1990_exposed_F_coarse, gdf_pop_1990_exposed_CF_coarse,
                                     region_utm, background_utm, beira_utm, flood_extent):
    print("Plotting Fig 2b: percentage attributable exposed population")
    gdf_CF_pop = gdf_pop_1990_exposed_F_coarse.copy()
    gdf_CF_clim = gdf_pop_2019_exposed_CF_coarse.copy()
    gdf_CF_clim_pop = gdf_pop_1990_exposed_CF_coarse.copy()

    gdf_CF_clim["pct_attr"] = ((gdf_pop_2019_exposed_F_coarse["exposed_population"]
                                 - gdf_CF_clim["exposed_population"])
                                / gdf_pop_2019_exposed_F_coarse["exposed_population"] * 100)
    gdf_CF_pop["pct_attr"] = ((gdf_pop_2019_exposed_F_coarse["exposed_population"]
                                - gdf_CF_pop["exposed_population"])
                               / gdf_pop_2019_exposed_F_coarse["exposed_population"] * 100)
    gdf_CF_clim_pop["pct_attr"] = ((gdf_pop_2019_exposed_F_coarse["exposed_population"]
                                     - gdf_CF_clim_pop["exposed_population"])
                                    / gdf_pop_2019_exposed_F_coarse["exposed_population"] * 100)

    datasets = [
        ("Climate change", gdf_CF_clim),
        ("Population change", gdf_CF_pop),
        ("Climate & population change", gdf_CF_clim_pop),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(11, 5), dpi=300, constrained_layout=True,
                             subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

    norm_diff = TwoSlopeNorm(vmin=-100, vcenter=0, vmax=100)
    cmap_diff = plt.cm.RdBu_r
    subplot_labels = ["(a)", "(b)", "(c)"]

    for i, (title, gdf) in enumerate(datasets):
        ax = axes[i]
        gdf[gdf["pct_attr"] <= 0].plot(ax=ax, color="white", edgecolor="grey", linewidth=0.2, zorder=1)
        gdf[gdf["pct_attr"] > 0].plot(column="pct_attr", cmap=cmap_diff, norm=norm_diff, edgecolor="grey",
                                       linewidth=0.2, ax=ax, legend=False, zorder=2, rasterized=True)

        setup_map_axes(ax, region_utm, background_utm, flood_extent,
                       subplot_labels=[subplot_labels[i]], titles=[title])
        beira_utm.boundary.plot(ax=ax, edgecolor="black", linewidth=0.5, zorder=3, alpha=0.7)
        annotate_landmarks(ax)

    sm = ScalarMappable(cmap=cmap_diff, norm=norm_diff)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation="vertical", shrink=0.6, pad=0.02)
    cbar.set_label("Attributable exposed population [%]", fontsize=10)
    return fig


# ===================================================================== #
# SUPPLEMENTARY: Climate-only district attribution map                   #
# ===================================================================== #
#%%
def plot_supp_climate_only_district(districts_adm3_clipped, region_utm, background_utm, flood_extent):
    from matplotlib.colors import Normalize as Norm
    vmin = 0
    vmax = np.nanpercentile(districts_adm3_clipped["attr_clim"].values, 99)
    norm = Norm(vmin=vmin, vmax=vmax)
    reds = LinearSegmentedColormap.from_list("reds_only", plt.cm.RdBu_r(np.linspace(0.5, 1, 256)))

    fig, ax = plt.subplots(1, 1, figsize=(8, 4), dpi=300, constrained_layout=True,
                           subplot_kw={"projection": ccrs.UTM(36, southern_hemisphere=True)})

    districts_adm3_clipped.plot(column="attr_clim", cmap=reds, norm=norm, edgecolor="grey",
                                linewidth=0.2, ax=ax, legend=False, rasterized=True,
                                missing_kwds={"color": "white"})

    setup_map_axes(ax, region_utm, background_utm, flood_extent, titles=["Climate change"])

    sm = ScalarMappable(cmap=reds, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, shrink=0.7).set_label("Attributable relative exposed population [%]", fontsize=10)

    for _, row in districts_adm3_clipped.iterrows():
        x, y = row.geometry.centroid.x, row.geometry.centroid.y
        ax.text(x, y, f"{row['NAME_3']}\n{int(round(row['attr_clim']))} %",
                ha="center", va="center", fontsize=7,
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=0.3, boxstyle="round"))
    return fig


# ===================================================================== #
# SUPPLEMENTARY: Exposed pop > 1.5 m flood depth (3-panel)              #
# ===================================================================== #
#%%
def plot_supp_exposed_gt15m(gdf_pop_2019_exposed_F_coarse, gdf_pop_2019_exposed_CF_coarse,
                             region_utm, background_utm, bg_filtered_utm, flood_extent):
    gdf_F = gdf_pop_2019_exposed_F_coarse.copy()
    gdf_F.loc[gdf_F["relative_population"] <= 1.5, "avg_flood_depth"] = np.nan
    gdf_CF = gdf_pop_2019_exposed_CF_coarse.copy()
    gdf_CF.loc[gdf_CF["relative_population"] <= 1.5, "avg_flood_depth"] = np.nan
    gdf_CF["change_gt15m"] = gdf_F["relative_population"] - gdf_CF["relative_population"]

    norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_F["relative_population"]))
    cmap_change = mcolors.LinearSegmentedColormap.from_list("w2p", ["#ffffff", "#651F94"])
    norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_CF["change_gt15m"]))

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    gdf_F.plot(column="relative_population", cmap="Blues", norm=norm, linewidth=0.1,
               edgecolor="grey", ax=axes[0], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    gdf_CF.plot(column="relative_population", cmap="Blues", norm=norm, linewidth=0.1,
                edgecolor="grey", ax=axes[1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    gdf_CF.plot(column="change_gt15m", cmap=cmap_change, norm=norm_change, linewidth=0.1,
                edgecolor="grey", ax=axes[2], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

    setup_plain_axes(axes, region_utm, background_utm, bg_filtered_utm, flood_extent)

    for ax in axes[:2]:
        sm = ScalarMappable(cmap="Blues", norm=norm)
        sm._A = []
        plt.colorbar(sm, ax=ax, shrink=0.8).set_label("Exposed population > 1.5 m flood depth")
    sm = ScalarMappable(cmap=cmap_change, norm=norm_change)
    sm._A = []
    plt.colorbar(sm, ax=axes[2], shrink=0.8).set_label("Difference")

    axes[0].set_title("Factual")
    axes[1].set_title("No Climate Change")
    axes[2].set_title("Factual - Counterfactual")
    fig.suptitle("Population exposed to > 1.5 m flood depth", fontsize=12)
    plt.tight_layout()
    return fig


# ===================================================================== #
# SUPPLEMENTARY: Attributable % of exposed pop (3-panel)                 #
# ===================================================================== #
#%%
def plot_supp_attr_pct_exposed(gdf_pop_2019_exposed_F_coarse, gdf_pop_2019_exposed_CF_coarse,
                                region_utm, background_utm, bg_filtered_utm, flood_extent):
    gdf_F = gdf_pop_2019_exposed_F_coarse.copy()
    gdf_F.loc[gdf_F["avg_flood_depth"] == 0, "exposed_population"] = 0
    gdf_CF = gdf_pop_2019_exposed_CF_coarse.copy()
    gdf_CF.loc[gdf_CF["avg_flood_depth"] == 0, "exposed_population"] = 0
    gdf_CF["pct_change"] = ((gdf_F["exposed_population"] - gdf_CF["exposed_population"])
                             / gdf_F["exposed_population"] * 100)

    norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_F["exposed_population"]))
    cmap_change = mcolors.LinearSegmentedColormap.from_list("w2p", ["#ffffff", "#651F94"])
    norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_CF["pct_change"]))

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    gdf_F.plot(column="exposed_population", cmap="Blues", norm=norm, linewidth=0.1,
               edgecolor="grey", ax=axes[0], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    gdf_CF.plot(column="exposed_population", cmap="Blues", norm=norm, linewidth=0.1,
                edgecolor="grey", ax=axes[1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    gdf_CF.plot(column="pct_change", cmap=cmap_change, norm=norm_change, linewidth=0.1,
                edgecolor="grey", ax=axes[2], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

    setup_plain_axes(axes, region_utm, background_utm, bg_filtered_utm, flood_extent)

    for ax in axes[:2]:
        sm = ScalarMappable(cmap="Blues", norm=norm)
        sm._A = []
        plt.colorbar(sm, ax=ax, shrink=0.8).set_label("Exposed population")
    sm = ScalarMappable(cmap=cmap_change, norm=norm_change)
    sm._A = []
    plt.colorbar(sm, ax=axes[2], shrink=0.8).set_label("Attributable exposed population")

    axes[0].set_title("Factual")
    axes[1].set_title("No Climate Change")
    axes[2].set_title("(F - CF) / F * 100%")
    fig.suptitle("Total exposed population", fontsize=12)
    plt.tight_layout()
    return fig


# ===================================================================== #
# SUPPLEMENTARY: % cells > 1.5 m flood depth (3-panel)                  #
# ===================================================================== #
#%%
def plot_supp_pct_cells_gt15m(gdf_pop_2019_exposed_F_coarse, gdf_pop_2019_exposed_CF_coarse,
                               region_utm, background_utm, bg_filtered_utm, flood_extent):
    gdf_F = gdf_pop_2019_exposed_F_coarse.copy()
    gdf_CF = gdf_pop_2019_exposed_CF_coarse.copy()
    gdf_CF["change_pct15m"] = gdf_F["pct_cells_higher_1.5m"] - gdf_CF["pct_cells_higher_1.5m"]

    norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_F["pct_cells_higher_1.5m"]))
    cmap_change = mcolors.LinearSegmentedColormap.from_list("w2p", ["#ffffff", "#651F94"])
    norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_CF["change_pct15m"].quantile(0.99))

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    gdf_F.plot(column="pct_cells_higher_1.5m", cmap="Blues", norm=norm, linewidth=0.1,
               edgecolor="grey", ax=axes[0], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    gdf_CF.plot(column="pct_cells_higher_1.5m", cmap="Blues", norm=norm, linewidth=0.1,
                edgecolor="grey", ax=axes[1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    gdf_CF.plot(column="change_pct15m", cmap=cmap_change, norm=norm_change, linewidth=0.1,
                edgecolor="grey", ax=axes[2], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

    setup_plain_axes(axes, region_utm, background_utm, bg_filtered_utm, flood_extent)

    for ax in axes[:2]:
        sm = ScalarMappable(cmap="Blues", norm=norm)
        sm._A = []
        plt.colorbar(sm, ax=ax, shrink=0.8).set_label("% cells > 1.5 m flood depth")
    sm = ScalarMappable(cmap=cmap_change, norm=norm_change)
    sm._A = []
    plt.colorbar(sm, ax=axes[2], shrink=0.8).set_label("Difference in % cells > 1.5 m")

    axes[0].set_title("Factual")
    axes[1].set_title("No Climate Change")
    axes[2].set_title("Factual - Counterfactual")
    fig.suptitle("% cells > 1.5 m flood depth", fontsize=12)
    plt.tight_layout()
    return fig


# ===================================================================== #
# SUPPLEMENTARY: % cells flooded (3-panel)                               #
# ===================================================================== #
#%%
def plot_supp_pct_cells_flooded(gdf_pop_2019_exposed_F_coarse, gdf_pop_2019_exposed_CF_coarse,
                                 region_utm, background_utm, bg_filtered_utm, flood_extent):
    gdf_F = gdf_pop_2019_exposed_F_coarse.copy()
    gdf_CF = gdf_pop_2019_exposed_CF_coarse.copy()
    gdf_CF["change_pct_flooded"] = gdf_F["pct_cells_flooded"] - gdf_CF["pct_cells_flooded"]

    norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_F["pct_cells_flooded"]))
    cmap_change = mcolors.LinearSegmentedColormap.from_list("w2p", ["#ffffff", "#651F94"])
    norm_change = PowerNorm(gamma=0.5, vmin=0, vmax=gdf_CF["change_pct_flooded"].quantile(0.99))

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    gdf_F.plot(column="pct_cells_flooded", cmap="Blues", norm=norm, linewidth=0.1,
               edgecolor="grey", ax=axes[0], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    gdf_CF.plot(column="pct_cells_flooded", cmap="Blues", norm=norm, linewidth=0.1,
                edgecolor="grey", ax=axes[1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    gdf_CF.plot(column="change_pct_flooded", cmap=cmap_change, norm=norm_change, linewidth=0.1,
                edgecolor="grey", ax=axes[2], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

    setup_plain_axes(axes, region_utm, background_utm, bg_filtered_utm, flood_extent)

    for ax in axes[:2]:
        sm = ScalarMappable(cmap="Blues", norm=norm)
        sm._A = []
        plt.colorbar(sm, ax=ax, shrink=0.8).set_label("% cells flooded")
    sm = ScalarMappable(cmap=cmap_change, norm=norm_change)
    sm._A = []
    plt.colorbar(sm, ax=axes[2], shrink=0.8).set_label("Difference in % cells flooded")

    axes[0].set_title("Factual")
    axes[1].set_title("No Climate Change")
    axes[2].set_title("Factual - Counterfactual")
    fig.suptitle("% cells flooded", fontsize=12)
    plt.tight_layout()
    return fig



# ===================================================================== #
# SUPPLEMENTARY: District validation plot (raster overlay)               #
# ===================================================================== #
#%%
def plot_supp_district_validation(districts_adm3_filtered, pop_arrays,
                                   ra_exposed_pop_2019_F, flood_extent,
                                   flood_grid_transform, region_utm,
                                   background_utm, bg_filtered_utm, region):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    for ax in axes:
        background_utm.plot(ax=ax, color="#E0E0E0", zorder=0)
        bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
        districts_adm3_filtered.boundary.plot(ax=ax, color="orange", linewidth=2, zorder=2)
        region_utm.boundary.plot(ax=ax, color="lightblue", linewidth=0.5, zorder=2)

    # Mask raster to region
    region_mask = features.rasterize(
        [(geom, 1) for geom in region.geometry],
        out_shape=ra_exposed_pop_2019_F.shape,
        transform=flood_grid_transform,
        fill=0, all_touched=True, dtype=np.uint8,
    ).astype(bool)
    ra_exposed_pop_masked = np.where(region_mask, ra_exposed_pop_2019_F, np.nan)

    im = axes[0].imshow(ra_exposed_pop_masked, cmap="viridis", extent=flood_extent, origin="lower")
    plt.colorbar(im, ax=axes[0], shrink=0.8).set_label("Exposed population")
    axes[0].set_title("Factual Exposed Population")

    im = axes[1].imshow(pop_arrays[2019], cmap="Reds", extent=flood_extent, origin="lower")
    plt.colorbar(im, ax=axes[1], shrink=0.8).set_label("Total population")
    axes[1].set_title("Total 2019 Population")

    label_offsets = {"Sofala": (0, 10000), "Nhamatanda": (18000, -22000), "Estaquinha": (20000, -5000)}
    outside_districts = ["Nhamatanda", "Estaquinha"]

    for _, row in districts_adm3_filtered.iterrows():
        x, y = row.geometry.centroid.x, row.geometry.centroid.y
        if row["NAME_3"] in label_offsets:
            dx, dy = label_offsets[row["NAME_3"]]
            x += dx
            y += dy
        axes[0].text(x, y, f"{row['pop_exposed']:,.0f}", fontsize=8, ha="center", va="center",
                     color="white", fontweight="bold", zorder=5)
        axes[1].text(x, y, f"{row['pop_total']:,.0f}", fontsize=8, ha="center", va="center",
                     color="black", fontweight="bold", zorder=5)
        name_x = x - 8000 if row["NAME_3"] in outside_districts else x
        name_y = y - 3000
        for ax in axes:
            ax.text(name_x, name_y, row["NAME_3"], fontsize=8, ha="center", va="center",
                    color="coral", fontweight="bold", zorder=5)

    plt.tight_layout()
    return fig



# ===================================================================== #
# SUPPLEMENTARY: Uniform vs spatially differing pop — plain axes         #
# ===================================================================== #
#%%
def plot_supp_uniform_vs_spatial_plain(gdf_pop_2019_exposed_F_coarse, gdf_pop_2019_exposed_F_uniform_coarse,
                                        region_utm, background_utm, bg_filtered_utm, region, flood_extent):
    gdf_F = gdf_pop_2019_exposed_F_coarse.copy()
    gdf_F.loc[gdf_F["total_population"] == 0] = np.nan
    gdf_U = gdf_pop_2019_exposed_F_uniform_coarse.copy()
    gdf_U.loc[gdf_U["total_population"] == 0] = np.nan
    gdf_U["change_in_population"] = gdf_U["total_population"] - gdf_F["total_population"]

    cmap = mcolors.LinearSegmentedColormap.from_list("w2p", ["#ffffff", "#651F94"])
    norm = PowerNorm(gamma=0.5, vmin=0, vmax=np.nanmax(gdf_F["total_population"]))
    cmap_change = plt.get_cmap("RdBu_r")
    norm_change = mcolors.TwoSlopeNorm(vmin=np.nanmin(gdf_U["change_in_population"]),
                                       vcenter=0, vmax=np.nanmax(gdf_U["change_in_population"]))

    fig, axes = plt.subplots(1, 3, figsize=(16, 6), sharey=True, constrained_layout=True)

    gdf_F.plot(column="total_population", cmap=cmap, norm=norm, linewidth=0.1, edgecolor="grey",
               ax=axes[0], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    gdf_U.plot(column="total_population", cmap=cmap, norm=norm, linewidth=0.1, edgecolor="grey",
               ax=axes[1], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})
    gdf_U.plot(column="change_in_population", cmap=cmap_change, norm=norm_change, linewidth=0.1,
               edgecolor="grey", ax=axes[2], zorder=2, missing_kwds={"color": "none", "edgecolor": "none"})

    for ax in axes:
        region.boundary.plot(ax=ax, color="black", linewidth=1)

    setup_plain_axes(axes, region_utm, background_utm, bg_filtered_utm, flood_extent)

    sm_pop = ScalarMappable(cmap=cmap, norm=norm)
    sm_pop._A = []
    fig.colorbar(sm_pop, ax=[axes[0], axes[1]], shrink=0.8, pad=0.02).set_label("Population (people per cell)")
    sm_change = ScalarMappable(cmap=cmap_change, norm=norm_change)
    sm_change._A = []
    fig.colorbar(sm_change, ax=axes[2], shrink=0.8, pad=0.02).set_label("Difference in population")

    axes[0].set_title("Factual 2019 population")
    axes[1].set_title("Uniform growth 2019 population")
    axes[2].set_title("Uniform growth - Factual")
    fig.suptitle("Total population in study region", fontsize=12)
    return fig

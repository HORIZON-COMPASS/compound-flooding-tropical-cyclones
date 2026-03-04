"""
Shared plotting helpers for the socio-economic attribution figures.
Reduces boilerplate for map setup, colorbars, landmarks, and flood-depth backgrounds.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as path_effects
import cartopy.crs as ccrs
from matplotlib.cm import ScalarMappable
from matplotlib.colors import PowerNorm, TwoSlopeNorm, BoundaryNorm, Normalize, LinearSegmentedColormap

from config import COLOURS


# ---------------------------------------------------------------------------
# Standard map formatting
# ---------------------------------------------------------------------------
def setup_map_axes(
    axes,
    region_utm,
    background_utm,
    flood_extent,
    subplot_labels=None,
    titles=None,
    show_left_labels_only=True,
    label_offset=(0, 1.02),
):
    """
    Apply standard map formatting to one or more Cartopy axes:
    region boundary, background, extent, gridlines, labels.
    """
    axes_arr = np.atleast_1d(axes)
    ncols = axes_arr.shape[-1] if axes_arr.ndim >= 2 else axes_arr.size
    nrows = axes_arr.shape[0] if axes_arr.ndim >= 2 else 1
    axes_flat = axes_arr.ravel()
    for i, ax in enumerate(axes_flat):
        background_utm.plot(ax=ax, color="#E0E0E0", zorder=0)
        region_utm.boundary.plot(ax=ax, edgecolor="black", linewidth=0.3)
        ax.set_extent(flood_extent, crs=ccrs.UTM(36, southern_hemisphere=True))

        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color="gray", alpha=0.5, linestyle="--")
        gl.right_labels = False
        gl.top_labels = False
        gl.xlabel_style = {"size": 9}
        gl.ylabel_style = {"size": 9}
        if show_left_labels_only and i % ncols != 0:
            gl.left_labels = False
        if i // ncols < nrows - 1:
            gl.bottom_labels = False

        if subplot_labels and i < len(subplot_labels):
            ax.text(
                label_offset[0],
                label_offset[1],
                subplot_labels[i],
                transform=ax.transAxes,
                fontsize=10,
                fontweight="bold",
                va="bottom",
                ha="left",
            )
        if titles and i < len(titles):
            ax.set_title(titles[i], fontsize=10)


# ---------------------------------------------------------------------------
# Standard non-projected map formatting (plain matplotlib, no Cartopy)
# ---------------------------------------------------------------------------
def setup_plain_axes(
    axes,
    region_utm,
    background_utm,
    bg_filtered_utm,
    flood_extent,
):
    """Apply standard formatting to non-Cartopy axes (simple imshow plots)."""
    xmin, xmax, ymin, ymax = flood_extent
    axes_flat = np.atleast_1d(axes).ravel()
    for ax in axes_flat:
        background_utm.plot(ax=ax, color="#E0E0E0", zorder=0)
        bg_filtered_utm.boundary.plot(ax=ax, color="#B0B0B0", linewidth=0.5, zorder=1)
        region_utm.boundary.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)


# ---------------------------------------------------------------------------
# Landmark annotations (Beira, Buzi River, Pungwe River)
# ---------------------------------------------------------------------------
def annotate_landmarks(ax, include_beira_district=True):
    """Add city and river markers to a Cartopy axis."""
    landmarks = [
        (34.862, -19.833, 34.852, -19.89, "Beira"),
        (34.43, -19.89, 34.44, -19.87, "Buzi River"),
        (34.543, -19.545, 34.554, -19.52, "Pungwe River"),
    ]
    for mx, my, tx, ty, label in landmarks:
        ax.plot(
            mx, my,
            marker="o", color="black", markersize=3, markeredgecolor="white",
            transform=ccrs.PlateCarree(), zorder=5,
        )
        text = ax.text(tx, ty, label, transform=ccrs.PlateCarree(), fontsize=8, color="black", zorder=5)
        text.set_path_effects([
            path_effects.Stroke(linewidth=3, foreground="white"),
            path_effects.Normal(),
        ])

    if include_beira_district:
        ax.text(
            34.975, -19.66, "Beira District",
            fontsize=8, ha="center", va="center", style="italic",
            transform=ccrs.PlateCarree(), zorder=4, color="#5C5C5C",
        )


# ---------------------------------------------------------------------------
# Flood-depth background shading (Low / Medium / High)
# ---------------------------------------------------------------------------
def add_flood_depth_bg(ax, ymin, ymax, label_y_frac=0.03):
    """Add Low / Medium / High shaded flood-depth backgrounds to a line-plot axis."""
    x_bg = np.linspace(0, 3.5, 500)
    low = x_bg < 0.5
    mid = (x_bg >= 0.5) & (x_bg < 1.5)
    high = x_bg >= 1.5

    ax.fill_between(x_bg[low], ymin, ymax, color="#d9d9d9", alpha=0.3)
    ax.fill_between(x_bg[mid], ymin, ymax, color="#b3b3b3", alpha=0.3)
    ax.fill_between(x_bg[high], ymin, ymax, color="#808080", alpha=0.3)

    label_y = ymin + (ymax - ymin) * label_y_frac
    bbox = dict(boxstyle="round", pad=0.15, facecolor="white", edgecolor="none", alpha=0.7)
    ax.text(0.25, label_y, "Low", ha="center", fontweight="bold", color="#5C5C5C", fontsize=10, bbox=bbox)
    ax.text(1.0, label_y, "Medium", ha="center", fontweight="bold", color="#5C5C5C", fontsize=10, bbox=bbox)
    ax.text(2.5, label_y, "High", ha="center", fontweight="bold", color="#5C5C5C", fontsize=10, bbox=bbox)


# ---------------------------------------------------------------------------
# Peak-matching utilities for flood-depth curves
# ---------------------------------------------------------------------------
def match_peaks_by_x(x, peaksA, peaksB):
    """Match peaks between two scenarios by closest x position."""
    xsA = x[peaksA]
    xsB = x[peaksB]
    matched = []
    used_B = set()

    for iA, xA in zip(peaksA, xsA):
        if len(xsB) > 0:
            diffs = [(abs(xA - xB), iB) for xB, iB in zip(xsB, peaksB) if iB not in used_B]
            if diffs:
                _, bestB = min(diffs, key=lambda t: t[0])
                used_B.add(bestB)
                matched.append((iA, bestB))
            else:
                matched.append((iA, None))
        else:
            matched.append((iA, None))

    for iB in peaksB:
        if iB not in used_B:
            matched.append((None, iB))

    return matched


def draw_peak_arrows(
    x, yA, yB, peaksA, peaksB, color, label_prefix,
    label_offsets, arrow_offsets, min_diff=0.01,
):
    """Draw arrows between matched peaks of two curves."""
    pairs = match_peaks_by_x(x, peaksA, peaksB)

    for i, (iA, iB) in enumerate(pairs):
        if iA is not None and iB is not None:
            xA, yA_val = x[iA], yA[iA]
            xB, yB_val = x[iB], yB[iB]
        elif iA is not None:
            xA, yA_val = x[iA], yA[iA]
            xB, yB_val = xA, yB[iA]
        elif iB is not None:
            xB, yB_val = x[iB], yB[iB]
            xA, yA_val = xB, yA[iB]
        else:
            continue

        diff = abs(yA_val - yB_val)
        if diff < min_diff:
            continue

        offset_x_arrow, offset_y_arrow = arrow_offsets.get(i, (0.0, 0.0))
        plt.annotate(
            "",
            xy=(xB + offset_x_arrow, yB_val + offset_y_arrow),
            xytext=(xA, yA_val),
            arrowprops=dict(arrowstyle="->", lw=1.4, color=color),
        )

        offset_x, offset_y = label_offsets.get(i, (0.03, 0.0))
        xm = xA + offset_x
        ym = (yA_val + yB_val) / 2 + offset_y
        plt.text(
            xm, ym, f"{label_prefix}",
            ha="left", va="center", fontsize=9, color=color, fontweight="bold",
            bbox=dict(boxstyle="round", pad=0.15, facecolor="lightgrey", edgecolor="none", alpha=0.5),
        )


# ---------------------------------------------------------------------------
# Discrete colormap builders
# ---------------------------------------------------------------------------
def make_discrete_pop_cmap(base_cmap, vmax, transparent_zero=True):
    """Build a discrete integer-step colormap for population counts."""
    bins = np.arange(0, int(vmax) + 1, 1)
    colors = base_cmap(np.linspace(0, 1, len(bins) - 1))
    if transparent_zero and len(colors) > 0:
        colors[0] = [1, 1, 1, 0]
    cmap = mcolors.ListedColormap(colors)
    norm = BoundaryNorm(bins, cmap.N, extend="neither")
    return cmap, norm, bins


def make_discrete_diff_cmap(vmax, transparent_zero=True):
    """Build a discrete integer-step diverging colormap for differences."""
    bins = np.arange(-int(vmax), int(vmax) + 1, 1)
    colors = plt.cm.RdBu_r(np.linspace(0, 1, len(bins) - 1))
    if transparent_zero:
        mid_idx = np.where(bins[:-1] == 0)[0]
        if len(mid_idx) > 0:
            colors[mid_idx[0]] = [1, 1, 1, 0]
    cmap = mcolors.ListedColormap(colors)
    norm = BoundaryNorm(bins, cmap.N, extend="neither")
    return cmap, norm, bins

"""
Inductively Coupled Plasma (ICP-OES / ICP-MS) Quantitative Analysis and Publication Plotting.

Implements:
- Multi-point calibration regression, blank subtraction, LOD / LOQ estimation
- Replicate statistical analysis (mean, standard deviation, RSD%)
- Dilution factor corrections and catalyst metal leaching quantification
  (dissolved mass, mass-normalized loss, percentage of initial inventory leached)
- Publication-quality calibration curves and leaching comparison figures.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats

from parsers.characterization.common import CharacterizationData
from modules.characterization.elemental_analysis import (
    calculate_catalyst_leaching_common,
    calculate_detection_limits,
    convert_concentration_units,
    linear_calibration_fit,
)
from scripts.utils import finalize_figure_layout, apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


def fit_icp_calibration(
    standards_concentration: Sequence[float],
    standards_intensity: Sequence[float],
    blank_intensities: Optional[Sequence[float]] = None,
) -> Dict[str, Any]:
    """
    Perform linear regression on ICP analytical calibration standards.

    Parameters
    ----------
    standards_concentration : sequence of float
        Known standard concentrations (e.g. in mg/L or ppm).
    standards_intensity : sequence of float
        Measured signal intensity/counts.
    blank_intensities : sequence of float, optional
        Replicate blank readings to evaluate background noise and detection limits.

    Returns
    -------
    Dict[str, Any]
        Regression parameters, R^2, slope, intercept, LOD, LOQ, residuals.
    """
    c_arr = np.asarray(standards_concentration, dtype=float)
    i_arr = np.asarray(standards_intensity, dtype=float)

    if len(c_arr) < 3:
        raise ValueError("At least 3 calibration standards required for regression.")

    fit_res = linear_calibration_fit(c_arr, i_arr, zero_intercept=False)

    lod = None
    loq = None
    if blank_intensities is not None and len(blank_intensities) >= 2:
        b_arr = np.asarray(blank_intensities, dtype=float)
        b_std = float(np.std(b_arr, ddof=1))
        det_limits = calculate_detection_limits(
            slope=fit_res["slope"],
            blank_stdev=b_std,
            method="blank_stdev",
            k_lod=3.3,
            k_loq=10.0,
        )
        lod = det_limits["lod_conc"]
        loq = det_limits["loq_conc"]

    return {
        "slope": fit_res["slope"],
        "intercept": fit_res["intercept"],
        "r_squared": fit_res["r_squared"],
        "std_err": fit_res["std_err_slope"],
        "lod_conc": lod,
        "loq_conc": loq,
        "standards_c": c_arr,
        "standards_i": i_arr,
        "fit_i": fit_res["y_pred"],
        "residuals": fit_res["residuals"],
    }


def calculate_catalyst_leaching(
    measured_concentration_mg_L: float,
    electrolyte_volume_mL: float,
    catalyst_mass_mg: float,
    initial_metal_wt_fraction: float,
    dilution_factor: float = 1.0,
    blank_concentration_mg_L: float = 0.0,
) -> Dict[str, Any]:
    """
    Calculate quantitative catalyst metal dissolution metrics.

    Parameters
    ----------
    measured_concentration_mg_L : float
        Measured concentration in solution (mg/L).
    electrolyte_volume_mL : float
        Total volume of testing electrolyte in milliliters.
    catalyst_mass_mg : float
        Total initial catalyst mass deposited/tested in milligrams.
    initial_metal_wt_fraction : float
        Initial weight fraction of the target metal in catalyst (e.g., 0.20 for 20 wt% Pt/C).
    dilution_factor : float, default 1.0
        Sample dilution factor applied prior to analysis.
    blank_concentration_mg_L : float, default 0.0
        Analytical blank concentration in mg/L.

    Returns
    -------
    Dict[str, Any]
        - 'corrected_concentration_mg_L': Corrected solution concentration.
        - 'dissolved_metal_mass_ug': Total dissolved mass of metal in micrograms.
        - 'mass_normalized_leaching_ug_per_mg': Leached metal per mg total catalyst.
        - 'percentage_metal_leached': Fraction of initial target metal leached (%).
    """
    return calculate_catalyst_leaching_common(
        measured_concentration_mg_L=measured_concentration_mg_L,
        electrolyte_volume_mL=electrolyte_volume_mL,
        catalyst_mass_mg=catalyst_mass_mg,
        initial_metal_wt_fraction=initial_metal_wt_fraction,
        dilution_factor=dilution_factor,
        blank_concentration_mg_L=blank_concentration_mg_L,
    )


def plot_icp_calibration(
    cal_result: Dict[str, Any],
    element: str = "Analyte",
    unit: str = "mg/L",
    title: Optional[str] = None,
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    """
    Generate publication 2-panel ICP calibration curve with fit line, R^2, LOD, and residuals.
    """
    apply_publication_style()
    fig, (ax_top, ax_bot) = create_publication_figure(
        columns="single",
        aspect_ratio=1.35,
        nrows=2,
        ncols=1,
        sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.0], "hspace": 0.08},
    )

    c = cal_result["standards_c"]
    i_obs = cal_result["standards_i"]
    i_fit = cal_result["fit_i"]
    res = cal_result["residuals"]

    # Top panel
    ax_top.scatter(c, i_obs, color="#0072B2", s=25, label="Standards", zorder=3)
    ax_top.plot(c, i_fit, color="#D55E00", linewidth=1.4, label="Linear Fit", zorder=2)

    r2 = cal_result["r_squared"]
    slope = cal_result["slope"]
    annot = f"$R^2 = {r2:.5f}$\nSlope: {slope:.1f}"
    if cal_result.get("lod_conc") is not None:
        annot += f"\nLOD: {cal_result['lod_conc']:.3f} {unit}"

    ax_top.text(
        0.05,
        0.80,
        annot,
        transform=ax_top.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.9),
    )

    ax_top.set_ylabel(r"Intensity $[\mathrm{counts}]$")
    ttl = title or f"ICP Calibration: {element}"
    ax_top.set_title(ttl, fontsize=9, pad=6)
    ax_top.legend(loc="lower right", frameon=False, fontsize=8)

    # Bottom panel: Residuals
    ax_bot.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax_bot.scatter(c, res, color="#009E73", s=18, zorder=3)
    ax_bot.set_xlabel(rf"Concentration $[\mathrm{{{unit}}}]$")
    ax_bot.set_ylabel(r"$\Delta I$")

    finalize_figure_layout(fig)
    return fig, (ax_top, ax_bot)


def plot_icp_leaching_comparison(
    samples: Sequence[str],
    leached_percentages: Sequence[float],
    error_bars: Optional[Sequence[float]] = None,
    element: str = "Metal",
    title: Optional[str] = None,
    color: str = "#D55E00",
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate publication bar chart comparing catalyst metal leaching percentages across materials.
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        created_fig = True
    else:
        fig = ax.get_figure()

    x_pos = np.arange(len(samples))
    y_vals = np.asarray(leached_percentages, dtype=float)
    y_err = np.asarray(error_bars, dtype=float) if error_bars is not None else None

    bars = ax.bar(x_pos, y_vals, yerr=y_err, capsize=3, color=color, edgecolor="black", linewidth=0.8, alpha=0.85, zorder=3)

    # Annotate values above bars
    for idx, rect in enumerate(bars):
        h = rect.get_height()
        err_h = y_err[idx] if y_err is not None else 0.0
        ax.annotate(
            f"{h:.1f}%",
            xy=(rect.get_x() + rect.get_width() / 2.0, h + err_h),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=7,
        )

    ax.set_xticks(x_pos)
    ax.set_xticklabels(samples, rotation=15, ha="right", fontsize=8)
    ax.set_ylabel(rf"{element} Leached $[\%]$")
    if title:
        ax.set_title(title, fontsize=9, pad=6)

    ax.set_ylim(0, max(y_vals) * 1.25 if len(y_vals) > 0 and max(y_vals) > 0 else 10.0)

    if created_fig:
        finalize_figure_layout(fig)
    return fig, ax

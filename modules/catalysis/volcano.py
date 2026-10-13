"""
Catalysis Volcano Plot Visualization and Microkinetic Descriptor Modeling.

Provides Sabatier volcano curves (activity or limiting potential versus descriptor
adsorption energies), two-leg theoretical volcano envelopes, optimal apex annotation,
and publication-quality catalyst point labeling.
"""

from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.config import get_figure_dimensions
from scripts.preprocess import KB_EV_PER_K
from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


def evaluate_volcano_model(
    x_range: np.ndarray,
    left_slope: float,
    left_intercept: float,
    right_slope: float,
    right_intercept: float,
) -> Tuple[np.ndarray, Tuple[float, float]]:
    """
    Evaluate theoretical two-leg Sabatier volcano relation:
    y = min(left_leg, right_leg) for limiting overpotential / max for activity.

    Parameters
    ----------
    x_range : np.ndarray
        Array of descriptor energy values.
    left_slope, left_intercept : float
        Parameters for the strong-binding leg.
    right_slope, right_intercept : float
        Parameters for the weak-binding leg.

    Returns
    -------
    Tuple[np.ndarray, Tuple[float, float]]
        (volcano_curve_y, (apex_x, apex_y))
    """
    left_leg = left_slope * x_range + left_intercept
    right_leg = right_slope * x_range + right_intercept

    # Intersection point (apex)
    if left_slope == right_slope:
        raise ValueError("Slopes must differ to form an apex.")

    apex_x = (right_intercept - left_intercept) / (left_slope - right_slope)
    apex_y = left_slope * apex_x + left_intercept

    # Typically in volcano curves, activity is bounded from above by the two legs (or overpotential from below)
    volcano_y = np.minimum(left_leg, right_leg)

    return volcano_y, (apex_x, apex_y)


def evaluate_orr_volcano(
    descriptor_delta_g_oh: Union[np.ndarray, Sequence[float]],
    scaling_ooh_intercept: float = 3.20,
    scaling_ooh_slope: float = 1.0,
    scaling_o_intercept: float = -0.20,
    scaling_o_slope: float = 2.0,
    two_leg_envelope: bool = True,
) -> Tuple[np.ndarray, np.ndarray, Tuple[float, float]]:
    """
    Evaluate associative 4-electron Oxygen Reduction Reaction (ORR) volcano model.

    Reaction steps (associative 4e- ORR mechanism):
    1. O2 + (H+ + e-) + * -> *OOH        (Delta G1 = Delta G_OOH - 4.92 eV)
    2. *OOH + (H+ + e-) -> *O + H2O       (Delta G2 = Delta G_O - Delta G_OOH)
    3. *O + (H+ + e-) -> *OH              (Delta G3 = Delta G_OH - Delta G_O)
    4. *OH + (H+ + e-) -> H2O + *         (Delta G4 = -Delta G_OH)

    Scaling relations (Nørskov et al., Viswanathan et al.):
    Delta G_OOH = scaling_ooh_slope * Delta G_OH + scaling_ooh_intercept  (typically 1.0 * dG_OH + 3.20 eV)
    Delta G_O   = scaling_o_slope * Delta G_OH + scaling_o_intercept      (typically 2.0 * dG_OH - 0.20 eV)

    Thermodynamic limiting potential:
    Two-leg Sabatier envelope: U_L = min(U_1, U_4) = min(4.92 - Delta G_OOH, Delta G_OH)
    All-step envelope: U_L = min(U_1, U_2, U_3, U_4)
    Theoretical overpotential:
    eta_ORR = 1.23 V - U_L

    Parameters
    ----------
    descriptor_delta_g_oh : array-like
        Adsorption free energy descriptor Delta G_OH* in eV.
    scaling_ooh_intercept : float, default 3.20
        Linear scaling intercept for Delta G_OOH* = slope * Delta G_OH* + intercept.
    scaling_ooh_slope : float, default 1.0
        Linear scaling slope for Delta G_OOH*.
    scaling_o_intercept : float, default -0.20
        Linear scaling intercept for Delta G_O*.
    scaling_o_slope : float, default 2.0
        Linear scaling slope for Delta G_O*.
    two_leg_envelope : bool, default True
        If True, uses the canonical two-leg Sabatier volcano envelope (U_1 and U_4).
        If False, includes all four elementary step potentials.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, Tuple[float, float]]
        (limiting_potentials_U_L, overpotentials_eta, (apex_dG_OH, apex_U_L))
    """
    dG_oh = np.asarray(descriptor_delta_g_oh, dtype=float)
    dG_ooh = scaling_ooh_slope * dG_oh + scaling_ooh_intercept
    dG_o = scaling_o_slope * dG_oh + scaling_o_intercept

    dG1 = dG_ooh - 4.92
    dG2 = dG_o - dG_ooh
    dG3 = dG_oh - dG_o
    dG4 = -dG_oh

    # Step potentials: U_i = -Delta G_i / e (in V)
    u1 = -dG1
    u2 = -dG2
    u3 = -dG3
    u4 = -dG4

    # Limiting potential is minimum of step potentials
    if two_leg_envelope:
        u_lim = np.minimum(u1, u4)
    else:
        u_lim = np.minimum(np.minimum(u1, u2), np.minimum(u3, u4))
    eta = 1.23 - u_lim

    # Apex is intersection of strong-binding leg (step 4: U_4 = dG_oh)
    # and weak-binding leg (step 1: U_1 = 4.92 - scaling_ooh_slope * dG_oh - scaling_ooh_intercept)
    apex_x = (4.92 - scaling_ooh_intercept) / (1.0 + scaling_ooh_slope)
    apex_y = apex_x  # since U_4 = dG_oh

    return u_lim, eta, (float(apex_x), float(apex_y))


def evaluate_her_volcano(
    descriptor_delta_g_h: Union[np.ndarray, Sequence[float]],
    temperature_k: float = 298.15,
    log_j0_max: float = 0.0,
    transfer_coefficient: float = 0.5,
) -> Tuple[np.ndarray, np.ndarray, Tuple[float, float]]:
    """
    Evaluate 2-electron Hydrogen Evolution Reaction (HER) volcano model.

    Reaction mechanism (Volmer-Heyrovsky / Volmer-Tafel):
    Volmer:    H+ + e- + * -> *H       (Delta G_H*)
    Heyrovsky: *H + H+ + e- -> H2 + *  (-Delta G_H*)

    Theoretical limiting overpotential:
    eta_HER = |Delta G_H*| / e

    Microkinetic exchange current density envelope:
    log10(j0) = log10(j0_max) - (transfer_coefficient * |Delta G_H*| / (k_B * T * ln(10)))

    Parameters
    ----------
    descriptor_delta_g_h : array-like
        Adsorption free energy descriptor Delta G_H* in eV.
    temperature_k : float, default 298.15
        Temperature in Kelvin.
    log_j0_max : float, default 0.0
        Log10 of maximum theoretical exchange current density at the apex.
    transfer_coefficient : float, default 0.5
        Electrochemical transfer coefficient alpha.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, Tuple[float, float]]
        (log10_j0, overpotentials_eta, (apex_dG_H, apex_log_j0))
    """
    dG_h = np.asarray(descriptor_delta_g_h, dtype=float)
    eta = np.abs(dG_h)

    # Slope in log10(j0) per eV: alpha / (k_B * T * ln(10))
    slope_ev = transfer_coefficient / (KB_EV_PER_K * float(temperature_k) * np.log(10.0))
    log_j0 = float(log_j0_max) - slope_ev * eta

    apex_x = 0.0
    apex_y = float(log_j0_max)

    return log_j0, eta, (apex_x, apex_y)



def plot_volcano(
    descriptors: Union[np.ndarray, Sequence[float], pd.Series],
    activities: Union[np.ndarray, Sequence[float], pd.Series],
    catalyst_names: Sequence[str],
    theoretical_model: Optional[Tuple[float, float, float, float]] = None,
    highlight_catalysts: Optional[Sequence[str]] = None,
    ax: Optional[plt.Axes] = None,
    x_label: str = r"Descriptor $\Delta G_{\mathrm{adsorbate}}$ (eV)",
    y_label: str = r"Theoretical Overpotential $-\eta$ (V) or Activity",
    title: Optional[str] = None,
    show_apex: bool = True,
    annotate_points: bool = True,
    marker_color: str = "#0072B2",
    highlight_color: str = "#D55E00",
    model_line_color: str = "#333333",
) -> Tuple[plt.Axes, Optional[Tuple[float, float]]]:
    """
    Plot catalytic volcano relationship with theoretical activity envelope and catalyst annotations.

    Parameters
    ----------
    descriptors : array-like
        Adsorption descriptor values (x-axis) in eV.
    activities : array-like
        Catalytic activity or limiting potential (y-axis).
    catalyst_names : sequence of str
        Labels for catalysts.
    theoretical_model : tuple of (left_slope, left_intercept, right_slope, right_intercept), optional
        Parameters defining the theoretical Sabatier volcano legs.
    highlight_catalysts : sequence of str, optional
        Names of top candidates to emphasize with distinct coloring.
    ax : plt.Axes, optional
        Matplotlib Axes to plot on.
    x_label, y_label : str
        Axis labels with units.
    title : str, optional
        Figure title.
    show_apex : bool, default True
        Whether to mark and label the volcano peak optimum.
    annotate_points : bool, default True
        Whether to print catalyst labels next to data points.
    marker_color, highlight_color, model_line_color : str
        Color styling.

    Returns
    -------
    Tuple[plt.Axes, Optional[Tuple[float, float]]]
        Matplotlib Axes and optimum apex coordinates (apex_x, apex_y).
    """
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    x_arr = np.asarray(descriptors, dtype=float)
    y_arr = np.asarray(activities, dtype=float)
    names = list(catalyst_names)
    highlights = set(highlight_catalysts) if highlight_catalysts else set()

    apex_coords: Optional[Tuple[float, float]] = None

    # Plot theoretical model if provided
    if theoretical_model is not None:
        l_s, l_i, r_s, r_i = theoretical_model
        x_min = min(x_arr.min() - 0.2, (r_i - l_i) / (l_s - r_s) - 0.5)
        x_max = max(x_arr.max() + 0.2, (r_i - l_i) / (l_s - r_s) + 0.5)
        x_dense = np.linspace(x_min, x_max, 250)

        y_model, apex_coords = evaluate_volcano_model(x_dense, l_s, l_i, r_s, r_i)
        ax.plot(x_dense, y_model, color=model_line_color, linestyle="-", linewidth=1.3, label="Microkinetic Model", zorder=2)

        if show_apex and apex_coords is not None:
            ax.scatter([apex_coords[0]], [apex_coords[1]], color="#E69F00", s=60, marker="*", edgecolor="#222222", zorder=5, label="Sabatier Optimum")
            ax.axvline(apex_coords[0], color="#999999", linestyle=":", linewidth=0.8, zorder=1)

    # Plot catalyst data points
    standard_x, standard_y, standard_names = [], [], []
    high_x, high_y, high_names = [], [], []

    for xi, yi, name in zip(x_arr, y_arr, names):
        if name in highlights:
            high_x.append(xi)
            high_y.append(yi)
            high_names.append(name)
        else:
            standard_x.append(xi)
            standard_y.append(yi)
            standard_names.append(name)

    if standard_x:
        ax.scatter(standard_x, standard_y, color=marker_color, s=26, edgecolor="#333333", linewidth=0.6, label="Catalysts", zorder=4)

    if high_x:
        ax.scatter(high_x, high_y, color=highlight_color, s=38, edgecolor="#222222", linewidth=0.8, label="Active Candidates", zorder=5)

    if annotate_points:
        for xi, yi, name in zip(x_arr, y_arr, names):
            is_high = name in highlights
            ax.annotate(
                name,
                xy=(xi, yi),
                xytext=(3.0, 3.0),
                textcoords="offset points",
                fontsize=6.8,
                fontweight="bold" if is_high else "normal",
                color="#000000" if is_high else "#333333",
                zorder=6,
            )

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if title:
        ax.set_title(title)

    ax.legend(loc="upper left", framealpha=0.85, fontsize=7.0)
    return ax, apex_coords


def plot_orr_volcano(
    descriptors_dG_OH: Union[np.ndarray, Sequence[float], pd.Series],
    catalyst_names: Sequence[str],
    highlight_catalysts: Optional[Sequence[str]] = None,
    ax: Optional[plt.Axes] = None,
    plot_quantity: str = "limiting_potential",
    title: Optional[str] = "Oxygen Reduction Reaction (ORR) Volcano",
    **kwargs: Any,
) -> Tuple[plt.Axes, Tuple[float, float]]:
    """
    Plot 4-electron ORR volcano curve and candidate catalysts.

    Parameters
    ----------
    descriptors_dG_OH : array-like
        Adsorption free energy Delta G_OH* in eV.
    catalyst_names : sequence of str
        Names of catalysts.
    highlight_catalysts : sequence of str, optional
        Names of top candidates to emphasize.
    ax : plt.Axes, optional
        Matplotlib Axes.
    plot_quantity : str, default 'limiting_potential'
        'limiting_potential' (U_L in V vs RHE) or 'overpotential' (eta in V).
    title : str, optional
        Figure title.

    Returns
    -------
    Tuple[plt.Axes, Tuple[float, float]]
    """
    dG_arr = np.asarray(descriptors_dG_OH, dtype=float)
    u_lim, eta, apex = evaluate_orr_volcano(dG_arr)

    # Dense model curve
    x_dense = np.linspace(min(dG_arr.min() - 0.2, 0.2), max(dG_arr.max() + 0.2, 1.6), 250)
    u_model, eta_model, apex_dense = evaluate_orr_volcano(x_dense)

    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    if plot_quantity == "limiting_potential":
        y_pts = u_lim
        y_model = u_model
        y_label = r"Limiting Potential $U_{\mathrm{L}}$ (V vs RHE)"
        apex_coords = apex_dense
    else:
        y_pts = eta
        y_model = eta_model
        y_label = r"Theoretical Overpotential $\eta^{\mathrm{ORR}}$ (V)"
        apex_coords = (apex_dense[0], 1.23 - apex_dense[1])

    ax.plot(x_dense, y_model, color="#333333", linewidth=1.3, label="ORR Sabatier Envelope", zorder=2)
    ax.scatter([apex_coords[0]], [apex_coords[1]], color="#E69F00", s=65, marker="*", edgecolor="#222222", zorder=5, label="Sabatier Optimum")
    ax.axvline(apex_coords[0], color="#999999", linestyle=":", linewidth=0.8, zorder=1)

    # Plot points
    highlights = set(highlight_catalysts) if highlight_catalysts else set()
    std_x, std_y, std_n = [], [], []
    hi_x, hi_y, hi_n = [], [], []

    for xi, yi, n in zip(dG_arr, y_pts, catalyst_names):
        if n in highlights:
            hi_x.append(xi)
            hi_y.append(yi)
            hi_n.append(n)
        else:
            std_x.append(xi)
            std_y.append(yi)
            std_n.append(n)

    if std_x:
        ax.scatter(std_x, std_y, color="#0072B2", s=28, edgecolor="#333333", linewidth=0.6, label="Catalysts", zorder=4)
    if hi_x:
        ax.scatter(hi_x, hi_y, color="#D55E00", s=40, edgecolor="#222222", linewidth=0.8, label="Active Candidates", zorder=5)

    for xi, yi, n in zip(dG_arr, y_pts, catalyst_names):
        is_hi = n in highlights
        ax.annotate(
            n,
            xy=(xi, yi),
            xytext=(3.0, 3.0),
            textcoords="offset points",
            fontsize=6.8,
            fontweight="bold" if is_hi else "normal",
            color="#000000" if is_hi else "#333333",
            zorder=6,
        )

    ax.set_xlabel(r"Descriptor $\Delta G_{\mathrm{OH}^*}$ (eV)")
    ax.set_ylabel(y_label)
    if title:
        ax.set_title(title)
    ax.legend(loc="lower left" if plot_quantity == "overpotential" else "lower right", framealpha=0.85, fontsize=7.0)
    return ax, apex_coords


def plot_her_volcano(
    descriptors_dG_H: Union[np.ndarray, Sequence[float], pd.Series],
    catalyst_names: Sequence[str],
    highlight_catalysts: Optional[Sequence[str]] = None,
    ax: Optional[plt.Axes] = None,
    temperature_k: float = 298.15,
    title: Optional[str] = "Hydrogen Evolution Reaction (HER) Volcano",
    **kwargs: Any,
) -> Tuple[plt.Axes, Tuple[float, float]]:
    """
    Plot 2-electron HER volcano curve and candidate catalysts.

    Parameters
    ----------
    descriptors_dG_H : array-like
        Adsorption free energy Delta G_H* in eV.
    catalyst_names : sequence of str
        Names of catalysts.
    highlight_catalysts : sequence of str, optional
        Names of top candidates to emphasize.
    ax : plt.Axes, optional
        Matplotlib Axes.
    temperature_k : float, default 298.15
        Temperature in Kelvin.
    title : str, optional
        Figure title.

    Returns
    -------
    Tuple[plt.Axes, Tuple[float, float]]
    """
    dG_arr = np.asarray(descriptors_dG_H, dtype=float)
    log_j0, eta, apex = evaluate_her_volcano(dG_arr, temperature_k=temperature_k)

    # Dense model curve
    x_dense = np.linspace(min(dG_arr.min() - 0.15, -0.8), max(dG_arr.max() + 0.15, 0.8), 250)
    log_j0_model, _, apex_dense = evaluate_her_volcano(x_dense, temperature_k=temperature_k)

    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    ax.plot(x_dense, log_j0_model, color="#333333", linewidth=1.3, label="Microkinetic HER Envelope", zorder=2)
    ax.scatter([apex_dense[0]], [apex_dense[1]], color="#E69F00", s=65, marker="*", edgecolor="#222222", zorder=5, label="Sabatier Optimum")
    ax.axvline(apex_dense[0], color="#999999", linestyle=":", linewidth=0.8, zorder=1)

    highlights = set(highlight_catalysts) if highlight_catalysts else set()
    std_x, std_y, std_n = [], [], []
    hi_x, hi_y, hi_n = [], [], []

    for xi, yi, n in zip(dG_arr, log_j0, catalyst_names):
        if n in highlights:
            hi_x.append(xi)
            hi_y.append(yi)
            hi_n.append(n)
        else:
            std_x.append(xi)
            std_y.append(yi)
            std_n.append(n)

    if std_x:
        ax.scatter(std_x, std_y, color="#0072B2", s=28, edgecolor="#333333", linewidth=0.6, label="Catalysts", zorder=4)
    if hi_x:
        ax.scatter(hi_x, hi_y, color="#D55E00", s=40, edgecolor="#222222", linewidth=0.8, label="Active Candidates", zorder=5)

    for xi, yi, n in zip(dG_arr, log_j0, catalyst_names):
        is_hi = n in highlights
        ax.annotate(
            n,
            xy=(xi, yi),
            xytext=(3.0, 3.0),
            textcoords="offset points",
            fontsize=6.8,
            fontweight="bold" if is_hi else "normal",
            color="#000000" if is_hi else "#333333",
            zorder=6,
        )

    ax.set_xlabel(r"Descriptor $\Delta G_{\mathrm{H}^*}$ (eV)")
    ax.set_ylabel(r"Exchange Current Density $\log_{10}(j_0 / j_0^{\mathrm{max}})$")
    if title:
        ax.set_title(title)
    ax.legend(loc="lower center", framealpha=0.85, fontsize=7.0)
    return ax, apex_dense


"""
Reaction Kinetics, Arrhenius Activation Energy Fitting, and Reactor Performance Visualization.

Implements rigorous Arrhenius parameter estimation (Ea, pre-exponential factor A)
with parameter uncertainties, residual diagnostics, conversion-selectivity curves,
and experimental-versus-model parity plots.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from scripts.config import get_figure_dimensions
from scripts.preprocess import GAS_CONSTANT_R_J_MOL_K
from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


def fit_arrhenius(
    temperatures_k: Union[np.ndarray, Sequence[float], pd.Series],
    rate_constants: Union[np.ndarray, Sequence[float], pd.Series],
    r_constant: float = GAS_CONSTANT_R_J_MOL_K,
) -> Dict[str, Any]:
    """
    Perform linear regression on Arrhenius equation: ln(k) = ln(A) - Ea / (R * T).

    Parameters
    ----------
    temperatures_k : array-like
        Absolute temperatures in Kelvin (T > 0).
    rate_constants : array-like
        Experimentally determined rate constants k (k > 0).
    r_constant : float, default 8.314462618 J/(mol*K)
        Universal gas constant.

    Returns
    -------
    Dict[str, Any]
        Dictionary with fitted parameters:
        - 'ea_j_mol': Activation energy in J/mol
        - 'ea_kj_mol': Activation energy in kJ/mol
        - 'ea_stderr_kj_mol': Standard error of Ea in kJ/mol
        - 'ln_a': Natural logarithm of pre-exponential factor ln(A)
        - 'ln_a_stderr': Standard error of ln(A)
        - 'pre_exponential_a': Pre-exponential factor A
        - 'r_squared': Coefficient of determination R^2
        - 'p_value': Two-sided p-value for the slope
        - 'residuals': Array of observed ln(k) - predicted ln(k)
        - 'inv_t': 1 / T values in K^-1
        - 'ln_k': ln(k) values
    """
    t_arr = np.asarray(temperatures_k, dtype=float)
    k_arr = np.asarray(rate_constants, dtype=float)

    if (t_arr <= 0).any():
        raise ValueError("Temperatures must be strictly positive in Kelvin.")
    if (k_arr <= 0).any():
        raise ValueError("Rate constants k must be strictly positive to compute ln(k).")
    if len(t_arr) != len(k_arr):
        raise ValueError(f"Lengths differ: {len(t_arr)} temperatures vs {len(k_arr)} rate constants.")
    if len(t_arr) < 3:
        raise ValueError("At least 3 temperature points are required for Arrhenius statistics.")

    inv_t = 1.0 / t_arr
    ln_k = np.log(k_arr)

    # Ordinary Least Squares regression
    reg = stats.linregress(inv_t, ln_k)

    slope = reg.slope          # -Ea / R
    intercept = reg.intercept  # ln(A)
    stderr_slope = reg.stderr
    stderr_intercept = reg.intercept_stderr

    ea_j = -slope * r_constant
    ea_kj = ea_j / 1000.0
    ea_stderr_kj = (stderr_slope * r_constant) / 1000.0

    predicted_ln_k = slope * inv_t + intercept
    residuals = ln_k - predicted_ln_k
    r2 = float(reg.rvalue ** 2)

    return {
        "ea_j_mol": float(ea_j),
        "ea_kj_mol": float(ea_kj),
        "ea_stderr_kj_mol": float(ea_stderr_kj),
        "ln_a": float(intercept),
        "ln_a_stderr": float(stderr_intercept),
        "pre_exponential_a": float(np.exp(intercept)),
        "r_squared": float(r2),
        "p_value": float(reg.pvalue),
        "residuals": residuals,
        "inv_t": inv_t,
        "ln_k": ln_k,
    }


def plot_arrhenius(
    temperatures_k: Union[np.ndarray, Sequence[float], pd.Series],
    rate_constants: Union[np.ndarray, Sequence[float], pd.Series],
    rate_constant_units: str = r"s$^{-1}$",
    energy_unit: str = "kJ/mol",
    show_residuals: bool = True,
    fig_axes: Optional[Tuple[plt.Figure, Union[plt.Axes, Sequence[plt.Axes]]]] = None,
    marker_color: str = "#0072B2",
    line_color: str = "#D55E00",
    title: Optional[str] = None,
) -> Tuple[plt.Figure, Dict[str, Any]]:
    """
    Generate publication-ready Arrhenius plot with regression line and residual diagnostics.

    Parameters
    ----------
    temperatures_k : array-like
        Absolute temperatures in Kelvin.
    rate_constants : array-like
        Rate constants k.
    rate_constant_units : str
        Units of rate constant for axis label.
    energy_unit : str, default 'kJ/mol'
        Unit for reporting Ea ('kJ/mol' or 'kcal/mol').
    show_residuals : bool, default True
        Whether to generate a companion residual plot underneath.
    fig_axes : tuple, optional
        Pre-existing (fig, axes) to plot on.
    marker_color, line_color : str
        Hex color strings.
    title : str, optional
        Figure title.

    Returns
    -------
    Tuple[plt.Figure, Dict[str, Any]]
        Matplotlib Figure and Arrhenius fitting results dictionary.
    """
    fit_res = fit_arrhenius(temperatures_k, rate_constants)

    inv_t_1000 = fit_res["inv_t"] * 1000.0  # Common scientific standard: 1000 / T (K^-1)
    ln_k = fit_res["ln_k"]

    if fig_axes is None:
        if show_residuals:
            fig, (ax_main, ax_res) = create_publication_figure(
                columns="single",
                aspect_ratio="1.1",
                nrows=2,
                ncols=1,
                sharex=True,
                gridspec_kw={"height_ratios": [3.0, 1.0], "hspace": 0.08},
            )
        else:
            fig, ax_main = create_publication_figure(columns="single", aspect_ratio="golden")
            ax_res = None
    else:
        fig, axes = fig_axes
        if isinstance(axes, (list, np.ndarray)):
            ax_main, ax_res = axes[0], axes[1]
        else:
            ax_main = axes
            ax_res = None

    # Main Arrhenius Plot
    ax_main.scatter(inv_t_1000, ln_k, color=marker_color, s=26, edgecolor="#222222", linewidth=0.6, zorder=4)

    # Fitted line
    x_dense_1000 = np.linspace(inv_t_1000.min() - 0.05, inv_t_1000.max() + 0.05, 100)
    x_dense_inv_t = x_dense_1000 / 1000.0
    slope = -fit_res["ea_j_mol"] / GAS_CONSTANT_R_J_MOL_K
    y_dense = slope * x_dense_inv_t + fit_res["ln_a"]

    ea_val = fit_res["ea_kj_mol"]
    ea_err = fit_res["ea_stderr_kj_mol"]
    r2_val = fit_res["r_squared"]

    label_text = (
        f"$E_{{\\mathrm{{a}}}} = {ea_val:.1f} \\pm {ea_err:.1f}\\text{{ {energy_unit}}}$\n"
        f"$\\ln(A) = {fit_res['ln_a']:.2f}$\n"
        f"$R^2 = {r2_val:.4f}$"
    )

    ax_main.plot(x_dense_1000, y_dense, color=line_color, linewidth=1.2, label=label_text, zorder=3)
    ax_main.set_ylabel(rf"$\ln(k / [{rate_constant_units}])$")
    ax_main.legend(loc="upper right", framealpha=0.9, fontsize=6.8)

    if title:
        ax_main.set_title(title)

    # Residuals Plot
    if ax_res is not None:
        residuals = fit_res["residuals"]
        ax_res.scatter(inv_t_1000, residuals, color=marker_color, s=20, edgecolor="#222222", linewidth=0.5)
        ax_res.axhline(0, color="#666666", linestyle="--", linewidth=0.7)
        ax_res.set_ylabel("Residuals")
        ax_res.set_xlabel(r"$1000 / T$ ($\mathrm{K}^{-1}$)")
        # Ticks for residuals
        res_limit = max(0.1, np.max(np.abs(residuals)) * 1.3)
        ax_res.set_ylim(-res_limit, res_limit)
    else:
        ax_main.set_xlabel(r"$1000 / T$ ($\mathrm{K}^{-1}$)")

    return fig, fit_res


def plot_conversion_selectivity(
    conversion: Union[np.ndarray, Sequence[float], pd.Series],
    selectivity: Union[np.ndarray, Sequence[float], pd.Series],
    yield_data: Optional[Union[np.ndarray, Sequence[float], pd.Series]] = None,
    catalyst_label: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    as_percentage: bool = True,
    marker_color: str = "#0072B2",
    line_color: str = "#0072B2",
) -> plt.Axes:
    """
    Plot catalyst Selectivity versus Conversion profile.

    Parameters
    ----------
    conversion : array-like
        Reactant conversion (0 to 1 fraction or 0 to 100%).
    selectivity : array-like
        Target product selectivity.
    yield_data : array-like, optional
        Target product yield (Y = X * S).
    catalyst_label : str, optional
        Label for legend.
    ax : plt.Axes, optional
        Matplotlib Axes.
    as_percentage : bool, default True
        Whether to format axes as 0-100% or 0-1.
    marker_color, line_color : str
        Hex color codes.

    Returns
    -------
    plt.Axes
    """
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    conv = np.asarray(conversion, dtype=float)
    sel = np.asarray(selectivity, dtype=float)

    # Normalize to 0-100% if requested
    if as_percentage:
        if conv.max() <= 1.0:
            conv = conv * 100.0
        if sel.max() <= 1.0:
            sel = sel * 100.0
        unit_str = " (%)"
    else:
        unit_str = ""

    # Sort by conversion for clean connecting curve
    order = np.argsort(conv)
    ax.plot(conv[order], sel[order], color=line_color, linestyle="-", linewidth=1.1, zorder=3)
    ax.scatter(conv, sel, color=marker_color, s=24, edgecolor="#222222", linewidth=0.5, label=catalyst_label, zorder=4)

    ax.set_xlabel(f"Conversion{unit_str}")
    ax.set_ylabel(f"Selectivity{unit_str}")

    if as_percentage:
        ax.set_xlim(0, 105)
        ax.set_ylim(0, 105)
    else:
        ax.set_xlim(0, 1.05)
        ax.set_ylim(0, 1.05)

    if catalyst_label:
        ax.legend(loc="best", framealpha=0.85)

    return ax


def plot_reactor_profiles(
    temperature_c: Union[np.ndarray, Sequence[float], pd.Series],
    conversion: Union[np.ndarray, Sequence[float], pd.Series],
    selectivity: Optional[Union[np.ndarray, Sequence[float], pd.Series]] = None,
    yield_val: Optional[Union[np.ndarray, Sequence[float], pd.Series]] = None,
    ax: Optional[plt.Axes] = None,
    title: Optional[str] = None,
) -> plt.Axes:
    """
    Plot reactor conversion, selectivity, and product yield versus reaction temperature.

    Parameters
    ----------
    temperature_c : array-like
        Reaction temperature in °C.
    conversion : array-like
        Reactant conversion (%).
    selectivity : array-like, optional
        Product selectivity (%).
    yield_val : array-like, optional
        Product yield (%). If None and selectivity provided, computed as X * S / 100.
    ax : plt.Axes, optional
        Matplotlib Axes.
    title : str, optional
        Title string.

    Returns
    -------
    plt.Axes
    """
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    t = np.asarray(temperature_c, dtype=float)
    x = np.asarray(conversion, dtype=float)

    # Conversion curve
    ax.plot(t, x, color="#0072B2", marker="o", markersize=3.8, linewidth=1.2, label="Conversion ($X$)")

    if selectivity is not None:
        s = np.asarray(selectivity, dtype=float)
        ax.plot(t, s, color="#009E73", marker="s", markersize=3.8, linewidth=1.2, label="Selectivity ($S$)")

        if yield_val is not None:
            y = np.asarray(yield_val, dtype=float)
        else:
            y = (x * s) / 100.0

        ax.plot(t, y, color="#D55E00", marker="^", markersize=3.8, linewidth=1.2, label="Yield ($Y$)")

    ax.set_xlabel(r"Temperature ($^\circ\mathrm{C}$)")
    ax.set_ylabel("Metric (%)")
    ax.set_ylim(0, 105)

    if title:
        ax.set_title(title)

    ax.legend(loc="center right", framealpha=0.85)
    return ax


def plot_parity_and_residuals(
    experimental: Union[np.ndarray, Sequence[float], pd.Series],
    predicted: Union[np.ndarray, Sequence[float], pd.Series],
    label: str = "Reaction Rate",
    units: str = r"mol\,g$^{-1}$\,s$^{-1}$",
    relative_error_band: float = 0.10,
    fig_axes: Optional[Tuple[plt.Figure, Union[plt.Axes, Sequence[plt.Axes]]]] = None,
) -> Tuple[plt.Figure, Dict[str, float]]:
    """
    Generate parity plot (predicted vs experimental) with +/- tolerance error cones and R^2.

    Parameters
    ----------
    experimental : array-like
        Observed experimental values.
    predicted : array-like
        Model predicted values.
    label : str
        Quantity name for axis labels.
    units : str
        Physical units string.
    relative_error_band : float, default 0.10
        Fractional error band (e.g. 0.10 = +/- 10% bounds).
    fig_axes : tuple, optional
        Existing Figure and Axes.

    Returns
    -------
    Tuple[plt.Figure, Dict[str, float]]
        Figure and validation statistics dictionary.
    """
    exp = np.asarray(experimental, dtype=float)
    pred = np.asarray(predicted, dtype=float)

    if len(exp) != len(pred):
        raise ValueError("Experimental and predicted arrays must have equal length.")

    if fig_axes is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="square")
    else:
        fig, ax = fig_axes

    val_min = min(float(exp.min()), float(pred.min())) * 0.9
    val_max = max(float(exp.max()), float(pred.max())) * 1.1

    # Parity line (y = x)
    line_x = np.linspace(val_min, val_max, 100)
    ax.plot(line_x, line_x, color="#333333", linestyle="-", linewidth=1.0, label="Parity ($y = x$)", zorder=2)

    # Relative error bounds (+/- band)
    upper_bound = line_x * (1.0 + relative_error_band)
    lower_bound = line_x * (1.0 - relative_error_band)
    ax.plot(line_x, upper_bound, color="#888888", linestyle="--", linewidth=0.7, label=rf"$\pm {int(relative_error_band*100)}\%$ Band", zorder=1)
    ax.plot(line_x, lower_bound, color="#888888", linestyle="--", linewidth=0.7, zorder=1)
    ax.fill_between(line_x, lower_bound, upper_bound, color="#EEEEEE", alpha=0.5, zorder=1)

    # Data points
    ax.scatter(exp, pred, color="#0072B2", s=28, edgecolor="#222222", linewidth=0.6, zorder=4)

    # Statistics
    ss_res = np.sum((exp - pred) ** 2)
    ss_tot = np.sum((exp - np.mean(exp)) ** 2)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
    mae = np.mean(np.abs(exp - pred))
    rmse = np.sqrt(np.mean((exp - pred) ** 2))

    stats_text = f"$R^2 = {r2:.4f}$\n$\\mathrm{{RMSE}} = {rmse:.3e}$\n$\\mathrm{{MAE}} = {mae:.3e}$"
    ax.annotate(
        stats_text,
        xy=(0.05, 0.75),
        xycoords="axes fraction",
        fontsize=6.8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#FAFAFA", edgecolor="#CCCCCC", alpha=0.9),
    )

    ax.set_xlim(val_min, val_max)
    ax.set_ylim(val_min, val_max)
    ax.set_xlabel(f"Experimental {label} ({units})")
    ax.set_ylabel(f"Predicted {label} ({units})")
    ax.legend(loc="lower right", framealpha=0.85, fontsize=6.8)

    return fig, {"r_squared": float(r2), "rmse": float(rmse), "mae": float(mae)}

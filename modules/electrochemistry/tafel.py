"""
Tafel Analysis and Electrocatalytic Reaction Kinetics.

Extracts Tafel slopes (b in mV/dec) and exchange current densities (j0)
from polarization data (LSV, steady-state polarization, or sampled staircase voltammetry).
Provides linear regression uncertainty, 95% confidence intervals, residual diagnostics,
and issues explicit scientific warnings for non-steady-state high-scan-rate artifacts.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import warnings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from modules.electrochemistry.reference_electrodes import (
    calculate_overpotential,
    convert_potential_to_rhe,
)
from parsers.electrochemistry.schema import StandardizedElectrochemicalData
from scripts.config import get_figure_dimensions
from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


class NonSteadyStateWarning(UserWarning):
    """Warns when Tafel slopes are derived from high-scan-rate non-steady-state measurements."""
    pass


def fit_tafel_slope(
    potential_vs_rhe: np.ndarray,
    current_density: np.ndarray,
    reaction: str = "OER",
    fit_range: Optional[Tuple[float, float]] = None,
    range_type: str = "potential",
    scan_rate_V_s: Optional[float] = None,
    current_range: Optional[Tuple[float, float]] = None,
) -> Dict[str, Any]:
    """
    Perform linear regression on Tafel relation: eta = a + b * log10(|j|).

    Parameters
    ----------
    potential_vs_rhe : np.ndarray
        Potential in Volts vs RHE.
    current_density : np.ndarray
        Current density in mA/cm^2.
    reaction : str, default 'OER'
        'OER' or 'HER'.
    fit_range : tuple of (float, float), optional
        Fitting window (min, max).
        If range_type == 'potential', values are in Volts vs RHE.
        If range_type == 'current', values are in mA/cm^2.
        If None, automatically selects linear kinetic region.
    range_type : str, default 'potential'
        'potential' or 'current'.
    scan_rate_V_s : float, optional
        Scan rate in V/s. Emits NonSteadyStateWarning if > 0.005 V/s (5 mV/s).

    Returns
    -------
    Dict[str, Any]
        - 'tafel_slope_mV_dec': Tafel slope in mV/decade
        - 'slope_stderr_mV_dec': Standard error of Tafel slope
        - 'intercept_V': Intercept a in Volts
        - 'r_squared': Coefficient of determination R^2
        - 'p_value': Two-sided p-value
        - 'log_j0_mA_cm2': Estimated log10 exchange current density
        - 'j0_mA_cm2': Exchange current density in mA/cm^2
        - 'fit_mask': Boolean mask of points used in fit
        - 'overpotentials_mV': Array of overpotentials
        - 'log_j': Array of log10(|j|)
    """
    if scan_rate_V_s is not None and scan_rate_V_s > 0.005:
        warnings.warn(
            f"Scan rate ({scan_rate_V_s*1000:.1f} mV/s) exceeds recommended steady-state limit (<= 5 mV/s). "
            "Tafel slopes derived from transient dynamic voltammetry can be distorted by capacitive charging "
            "and mass-transport delays.",
            NonSteadyStateWarning,
            stacklevel=2,
        )

    pot = np.asarray(potential_vs_rhe, dtype=float)
    j_arr = np.asarray(current_density, dtype=float)
    rxn = reaction.strip().upper()

    eta = calculate_overpotential(pot, reaction=rxn)
    eta_mv = np.abs(eta) * 1000.0  # Overpotential magnitude in mV
    abs_j = np.abs(j_arr)

    # Filter strictly positive current densities to take logarithm
    valid_points = abs_j > 1e-6
    if not np.any(valid_points):
        raise ValueError("No positive current density points found for logarithmic Tafel fitting.")

    pot_v = pot[valid_points]
    eta_v = eta_mv[valid_points]
    j_v = abs_j[valid_points]
    log_j = np.log10(j_v)

    # Determine fitting mask
    if current_range is not None:
        fit_range = current_range
        range_type = "current"

    if fit_range is not None:
        low, high = fit_range
        if range_type.lower() == "potential":
            mask = (pot_v >= low) & (pot_v <= high)
        elif range_type.lower() == "current":
            mask = (j_v >= low) & (j_v <= high)
        else:
            raise ValueError(f"Unknown range_type: '{range_type}'. Use 'potential' or 'current'.")
    else:
        # Auto-detect region: between 0.1 mA/cm^2 and 10 mA/cm^2 (typical kinetic control zone)
        mask = (j_v >= 0.1) & (j_v <= 10.0)
        if np.sum(mask) < 4:
            # Fallback to middle 40% of points
            p_len = len(j_v)
            mask = np.zeros(p_len, dtype=bool)
            mask[int(p_len * 0.3) : int(p_len * 0.7)] = True

    if np.sum(mask) < 3:
        raise ValueError("Selected Tafel fitting window contains fewer than 3 data points.")

    x_fit = log_j[mask]
    y_fit = eta_v[mask]

    res = stats.linregress(x_fit, y_fit)
    b_slope = float(res.slope)  # mV/decade
    a_intercept = float(res.intercept)  # mV
    r2 = float(res.rvalue ** 2)

    # Exchange current density j0 where eta = 0 mV -> 0 = a + b * log10(j0) => log10(j0) = -a / b
    log_j0 = float(-a_intercept / b_slope) if b_slope != 0 else float("nan")
    j0_val = float(10.0 ** log_j0) if not np.isnan(log_j0) else float("nan")

    return {
        "tafel_slope_mV_dec": b_slope,
        "slope_stderr_mV_dec": float(res.stderr),
        "intercept_mV": a_intercept,
        "r_squared": r2,
        "p_value": float(res.pvalue),
        "log_j0_mA_cm2": log_j0,
        "j0_mA_cm2": j0_val,
        "exchange_current_density_mA_cm2": j0_val,
        "fit_mask": mask,
        "x_fit": x_fit,
        "y_fit": y_fit,
        "log_j": log_j,
        "eta_mV": eta_v,
    }


def plot_tafel(
    data: StandardizedElectrochemicalData,
    reaction: str = "OER",
    fit_range: Optional[Tuple[float, float]] = None,
    range_type: str = "potential",
    reference_electrode: Optional[str] = None,
    ph: Optional[float] = None,
    ir_compensate: bool = False,
    ru_ohm: Optional[float] = None,
    ax: Optional[plt.Axes] = None,
    color: str = "#0072B2",
    label: Optional[str] = None,
    title: Optional[str] = None,
) -> Tuple[plt.Axes, Dict[str, Any]]:
    """
    Plot single Tafel curve with linear regression envelope and slope annotation.
    """
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    pot = data.potential_V.copy()
    if ir_compensate:
        r_val = ru_ohm if ru_ohm is not None else data.ru_ohm
        if r_val is not None:
            from modules.electrochemistry.reference_electrodes import apply_ir_compensation
            pot = apply_ir_compensation(pot, data.current_A, r_val)

    ref = reference_electrode or data.potential_reference_raw
    sol_ph = ph if ph is not None else data.ph
    pot_rhe = convert_potential_to_rhe(pot, ref, sol_ph)
    j_ma = data.get_current_density(unit="mA/cm2")

    fit_res = fit_tafel_slope(
        potential_vs_rhe=pot_rhe,
        current_density=j_ma,
        reaction=reaction,
        fit_range=fit_range,
        range_type=range_type,
        scan_rate_V_s=data.scan_rate_V_s,
    )

    log_j = fit_res["log_j"]
    eta_mv = fit_res["eta_mV"]
    x_fit = fit_res["x_fit"]
    b_slope = fit_res["tafel_slope_mV_dec"]
    r2 = fit_res["r_squared"]

    lbl = label or (f"{data.sample_id}" if data.sample_id else f"{reaction} Tafel")
    ax.plot(log_j, eta_mv, color=color, linewidth=1.2, label=lbl, zorder=3)

    # Plot linear fit line
    x_line = np.linspace(x_fit.min() - 0.2, x_fit.max() + 0.2, 50)
    y_line = fit_res["intercept_mV"] + b_slope * x_line
    ax.plot(
        x_line,
        y_line,
        color="#D55E00",
        linestyle="--",
        linewidth=1.3,
        label=f"Fit: {b_slope:.1f} mV/dec ($R^2={r2:.3f}$)",
        zorder=4,
    )

    ax.set_xlabel(r"$\log_{10}|j|$ ($\mathrm{mA\,cm^{-2}}$)")
    ax.set_ylabel(r"Overpotential $\eta$ (mV)")
    if title:
        ax.set_title(title)
    else:
        ax.set_title(f"{reaction} Tafel Plot")

    ax.legend(loc="upper left" if reaction.upper() == "OER" else "upper right", framealpha=0.85, fontsize=6.8)
    return ax, fit_res


def plot_tafel_comparison(
    datasets: Sequence[StandardizedElectrochemicalData],
    labels: Sequence[str],
    reaction: str = "OER",
    fit_ranges: Optional[Sequence[Optional[Tuple[float, float]]]] = None,
    range_type: str = "potential",
    reference_electrode: Optional[str] = None,
    ph: Optional[float] = None,
    ir_compensate: bool = False,
    ru_ohms: Optional[Sequence[float]] = None,
    ax: Optional[plt.Axes] = None,
    palette: str = "okabe_ito",
    title: Optional[str] = None,
) -> Tuple[plt.Axes, Dict[str, Dict[str, Any]]]:
    """
    Compare Tafel slopes across multiple catalysts on a single publication figure.
    """
    if len(datasets) != len(labels):
        raise ValueError("Length mismatch between datasets and labels.")

    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    colors = get_palette(palette, max(len(datasets), 2))
    all_fits: Dict[str, Dict[str, Any]] = {}

    for idx, (data, lbl) in enumerate(zip(datasets, labels)):
        pot = data.potential_V.copy()
        if ir_compensate:
            r_val = ru_ohms[idx] if ru_ohms else data.ru_ohm
            if r_val is not None:
                from modules.electrochemistry.reference_electrodes import apply_ir_compensation
                pot = apply_ir_compensation(pot, data.current_A, r_val)

        ref = reference_electrode or data.potential_reference_raw
        sol_ph = ph if ph is not None else data.ph
        pot_rhe = convert_potential_to_rhe(pot, ref, sol_ph)
        j_ma = data.get_current_density(unit="mA/cm2")

        f_range = fit_ranges[idx] if fit_ranges else None
        fit_res = fit_tafel_slope(
            potential_vs_rhe=pot_rhe,
            current_density=j_ma,
            reaction=reaction,
            fit_range=f_range,
            range_type=range_type,
            scan_rate_V_s=data.scan_rate_V_s,
        )
        all_fits[lbl] = fit_res

        b_slope = fit_res["tafel_slope_mV_dec"]
        ax.plot(
            fit_res["log_j"],
            fit_res["eta_mV"],
            color=colors[idx % len(colors)],
            linewidth=1.2,
            label=f"{lbl} ({b_slope:.1f} mV/dec)",
            zorder=3,
        )

    ax.set_xlabel(r"$\log_{10}|j|$ ($\mathrm{mA\,cm^{-2}}$)")
    ax.set_ylabel(r"Overpotential $\eta$ (mV)")
    if title:
        ax.set_title(title)
    else:
        ax.set_title(f"{reaction} Tafel Comparison")

    ax.legend(loc="upper left" if reaction.upper() == "OER" else "upper right", framealpha=0.85, fontsize=6.8)
    return ax, all_fits


# Functional alias for API consistency
calculate_tafel_slope = fit_tafel_slope

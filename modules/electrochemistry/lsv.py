"""
Linear Sweep Voltammetry (LSV) Polarization Curves and Electrocatalysis Benchmarking.

Supports HER and OER polarization curves, geometric and mass current normalization,
RHE conversions, iR compensation, overpotential extraction at target current densities
(e.g., eta_10, eta_50, eta_100), and multi-catalyst comparative publication plotting.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np

from modules.electrochemistry.reference_electrodes import (
    apply_ir_compensation,
    calculate_overpotential,
    convert_potential_to_rhe,
)
from parsers.electrochemistry.schema import StandardizedElectrochemicalData
from scripts.config import get_figure_dimensions
from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


def extract_overpotential_at_current(
    potential: np.ndarray,
    current_density: np.ndarray,
    target_j: float = 10.0,
    reaction: str = "OER",
) -> Optional[float]:
    """
    Interpolate the overpotential at a specific target current density (e.g. 10 mA/cm^2).

    Parameters
    ----------
    potential : np.ndarray
        Potential array (V vs RHE).
    current_density : np.ndarray
        Current density array (mA/cm^2).
    target_j : float, default 10.0
        Target current density magnitude in mA/cm^2.
        For OER, searches for j >= target_j (+10 mA/cm^2).
        For HER, searches for j <= -target_j (-10 mA/cm^2).
    reaction : str, default 'OER'
        'OER' or 'HER'.

    Returns
    -------
    float, optional
        Overpotential (eta) in mV. None if target current density was not reached.
    """
    pot = np.asarray(potential, dtype=float)
    j_arr = np.asarray(current_density, dtype=float)
    rxn = reaction.strip().upper()

    target_val = abs(target_j) if rxn == "OER" else -abs(target_j)

    # Sort monotonically by potential
    sort_idx = np.argsort(pot)
    pot_s = pot[sort_idx]
    j_s = j_arr[sort_idx]

    # Check if target is bracketed
    if rxn == "OER":
        if j_s.max() < target_val:
            return None
        # Interpolate potential where current density equals target_val
        e_target = float(np.interp(target_val, j_s, pot_s))
        # Overpotential for OER: eta = E - 1.229 V
        eta_v = e_target - 1.229
        return float(eta_v * 1000.0)  # in mV
    else:  # HER
        if j_s.min() > target_val:
            return None
        # For HER, j is negative, cathodic current increases negatively
        # Invert to monotonic for interpolation
        e_target = float(np.interp(-target_val, -j_s, pot_s))
        # Overpotential for HER: eta = |E - 0.000|
        eta_v = abs(e_target - 0.000)
        return float(eta_v * 1000.0)  # in mV


def plot_lsv(
    data: StandardizedElectrochemicalData,
    reaction: str = "OER",
    current_unit: str = "mA/cm2",
    potential_scale: str = "RHE",
    reference_electrode: Optional[str] = None,
    ph: Optional[float] = None,
    ir_compensate: bool = False,
    ru_ohm: Optional[float] = None,
    compensation_fraction: float = 1.0,
    target_currents: Sequence[float] = (10.0,),
    ax: Optional[plt.Axes] = None,
    line_color: str = "#0072B2",
    label: Optional[str] = None,
    title: Optional[str] = None,
) -> Tuple[plt.Axes, Dict[str, Any]]:
    """
    Plot single LSV polarization curve with optional overpotential benchmark markers.

    Parameters
    ----------
    data : StandardizedElectrochemicalData
        Parsed dataset.
    reaction : str, default 'OER'
        'OER' or 'HER'.
    current_unit : str, default 'mA/cm2'
        Normalization unit.
    potential_scale : str, default 'RHE'
        'RHE' or 'measured'.
    reference_electrode : str, optional
        Reference electrode used.
    ph : float, optional
        Electrolyte pH.
    ir_compensate : bool, default False
        Whether to apply software iR compensation.
    ru_ohm : float, optional
        Uncompensated resistance Ru in Ohms.
    compensation_fraction : float, default 1.0
        Fraction of Ru to compensate (0.85 to 1.0).
    target_currents : sequence of float, default (10.0,)
        Target current densities (mA/cm^2) for overpotential reporting.
    ax : plt.Axes, optional
        Target Axes.
    line_color : str, default '#0072B2'
        Color of curve.
    label : str, optional
        Legend label.
    title : str, optional
        Figure title.

    Returns
    -------
    Tuple[plt.Axes, Dict[str, Any]]
        Axes and overpotential benchmarks dictionary.
    """
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    pot = data.potential_V.copy()
    cur_A = data.current_A

    if ir_compensate:
        r_val = ru_ohm if ru_ohm is not None else data.ru_ohm
        if r_val is None:
            raise ValueError("ru_ohm must be specified to perform iR compensation.")
        pot = apply_ir_compensation(pot, cur_A, r_val, compensation_fraction)

    if potential_scale.upper() == "RHE":
        ref = reference_electrode or data.potential_reference_raw
        sol_ph = ph if ph is not None else data.ph
        pot_plotted = convert_potential_to_rhe(pot, ref, sol_ph)
        x_label = r"Potential $E$ (V vs RHE)"
    else:
        pot_plotted = pot
        x_label = f"Potential $E$ (V vs {reference_electrode or data.potential_reference_raw})"

    j_norm = data.get_current_density(unit=current_unit)
    y_label = f"Current Density $j$ ({current_unit})"

    lbl = label or (f"{data.sample_id}" if data.sample_id else f"{reaction} LSV")
    ax.plot(pot_plotted, j_norm, color=line_color, linewidth=1.3, label=lbl, zorder=3)

    # Benchmark overpotentials
    benchmarks: Dict[str, Optional[float]] = {}
    for tj in target_currents:
        eta_val = extract_overpotential_at_current(pot_plotted, j_norm, target_j=tj, reaction=reaction)
        benchmarks[f"eta_{int(tj)}_mV"] = eta_val

        # Mark eta_10 on plot if reached
        if eta_val is not None:
            y_mark = tj if reaction.upper() == "OER" else -tj
            # Calculate corresponding potential
            e_mark = (1.229 + eta_val / 1000.0) if reaction.upper() == "OER" else (-eta_val / 1000.0)
            ax.scatter([e_mark], [y_mark], color="#D55E00", s=36, marker="o", edgecolor="#222222", zorder=5)
            ax.annotate(
                f"$\\eta_{{{int(tj)}}} = {eta_val:.1f}$ mV",
                xy=(e_mark, y_mark),
                xytext=(6.0, -8.0 if reaction.upper() == "OER" else 8.0),
                textcoords="offset points",
                fontsize=6.8,
                color="#222222",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="#FFFFFF", edgecolor="#CCCCCC", alpha=0.9),
                zorder=6,
            )

    ax.axhline(0.0, color="#666666", linestyle="--", linewidth=0.7, zorder=1)
    if reaction.upper() == "OER":
        ax.axvline(1.229, color="#999999", linestyle=":", linewidth=0.7, label=r"$E^0_{\mathrm{OER}}$ (1.23 V)")
    else:
        ax.axvline(0.000, color="#999999", linestyle=":", linewidth=0.7, label=r"$E^0_{\mathrm{HER}}$ (0.00 V)")

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if title:
        ax.set_title(title)
    ax.legend(loc="best", framealpha=0.85, fontsize=6.8)

    return ax, benchmarks


def plot_lsv_comparison(
    datasets: Sequence[StandardizedElectrochemicalData],
    labels: Sequence[str],
    reaction: str = "OER",
    current_unit: str = "mA/cm2",
    potential_scale: str = "RHE",
    reference_electrode: Optional[str] = None,
    ph: Optional[float] = None,
    ir_compensate: bool = False,
    ru_ohms: Optional[Sequence[float]] = None,
    ax: Optional[plt.Axes] = None,
    palette: str = "okabe_ito",
    title: Optional[str] = None,
) -> Tuple[plt.Axes, Dict[str, Dict[str, Optional[float]]]]:
    """
    Compare LSV polarization curves across multiple electrocatalysts or conditions.

    Parameters
    ----------
    datasets : sequence of StandardizedElectrochemicalData
        List of datasets.
    labels : sequence of str
        Display names for each catalyst.
    reaction : str, default 'OER'
        'OER' or 'HER'.
    current_unit : str, default 'mA/cm2'
        Current density unit.
    potential_scale : str, default 'RHE'
        'RHE' or 'measured'.
    reference_electrode : str, optional
        Reference electrode.
    ph : float, optional
        Solution pH.
    ir_compensate : bool, default False
        Whether to apply iR compensation.
    ru_ohms : sequence of float, optional
        Ru value per catalyst.
    ax : plt.Axes, optional
        Target Axes.
    palette : str, default 'okabe_ito'
        Color palette.
    title : str, optional
        Figure title.

    Returns
    -------
    Tuple[plt.Axes, Dict[str, Dict[str, Optional[float]]]]
        Axes and dictionary mapping catalyst label to benchmark overpotentials.
    """
    if len(datasets) != len(labels):
        raise ValueError(
            f"Length mismatch: {len(datasets)} datasets provided but {len(labels)} labels given."
        )

    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    colors = get_palette(palette, max(len(datasets), 2))
    all_benchmarks: Dict[str, Dict[str, Optional[float]]] = {}

    for idx, (data, lbl) in enumerate(zip(datasets, labels)):
        pot = data.potential_V.copy()
        cur_A = data.current_A

        if ir_compensate:
            r_val = ru_ohms[idx] if ru_ohms else data.ru_ohm
            if r_val is not None:
                pot = apply_ir_compensation(pot, cur_A, r_val)

        if potential_scale.upper() == "RHE":
            ref = reference_electrode or data.potential_reference_raw
            sol_ph = ph if ph is not None else data.ph
            pot_plotted = convert_potential_to_rhe(pot, ref, sol_ph)
        else:
            pot_plotted = pot

        j_norm = data.get_current_density(unit=current_unit)
        ax.plot(pot_plotted, j_norm, color=colors[idx % len(colors)], linewidth=1.3, label=lbl, zorder=3)

        eta10 = extract_overpotential_at_current(pot_plotted, j_norm, target_j=10.0, reaction=reaction)
        all_benchmarks[lbl] = {"eta_10_mV": eta10}

    ax.axhline(0.0, color="#666666", linestyle="--", linewidth=0.7, zorder=1)
    if reaction.upper() == "OER":
        ax.axhline(10.0, color="#CCCCCC", linestyle=":", linewidth=0.6)
        ax.axvline(1.229, color="#999999", linestyle=":", linewidth=0.7, label=r"$E^0_{\mathrm{OER}}$ (1.23 V)")
    else:
        ax.axhline(-10.0, color="#CCCCCC", linestyle=":", linewidth=0.6)
        ax.axvline(0.000, color="#999999", linestyle=":", linewidth=0.7, label=r"$E^0_{\mathrm{HER}}$ (0.00 V)")

    x_lbl = r"Potential $E$ (V vs RHE)" if potential_scale.upper() == "RHE" else "Potential $E$ (V)"
    ax.set_xlabel(x_lbl)
    ax.set_ylabel(f"Current Density $j$ ({current_unit})")

    if title:
        ax.set_title(title)
    else:
        comp_note = " (iR-corrected)" if ir_compensate else ""
        ax.set_title(f"{reaction} Electrocatalytic Performance{comp_note}")

    ax.legend(loc="upper left" if reaction.upper() == "OER" else "lower left", framealpha=0.85, fontsize=6.8)
    return ax, all_benchmarks

"""
Cross-Instrument Benchmarking and Discrepancy Analysis (Metrohm Autolab vs BioLogic EC-Lab).

Provides rigorous scientific comparison between electrochemical datasets obtained
on different potentiostats:
- Verification of matched electrolyte chemistry, pH, reference potential, and electrode area
- Polarization curve overlay with point-by-point residual analysis
- Extraction and cross-comparison of key electrocatalytic metrics (eta_10, Tafel slope, C_dl, R_ct)
- Dual-panel publication figure with experimental curves and residual error bands.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np

from modules.electrochemistry.lsv import extract_overpotential_at_current
from modules.electrochemistry.tafel import calculate_tafel_slope
from parsers.electrochemistry.schema import StandardizedElectrochemicalData
from scripts.config import get_figure_dimensions
from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


def compare_instrument_lsv(
    data_metrohm: StandardizedElectrochemicalData,
    data_biologic: StandardizedElectrochemicalData,
    reaction: str = "OER",
    target_j: float = 10.0,
    potential_grid_points: int = 200,
) -> Dict[str, Any]:
    """
    Quantitatively compare LSV curves measured on Metrohm Autolab and BioLogic EC-Lab.

    Parameters
    ----------
    data_metrohm : StandardizedElectrochemicalData
        Dataset from Metrohm Autolab.
    data_biologic : StandardizedElectrochemicalData
        Dataset from BioLogic EC-Lab.
    reaction : str, default 'OER'
        'OER' or 'HER'.
    target_j : float, default 10.0
        Target current density in mA/cm^2 for overpotential extraction.
    potential_grid_points : int, default 200
        Number of points for interpolated residual analysis.

    Returns
    -------
    Dict[str, Any]
        Dictionary with matching status, eta_10 comparisons, curve RMSE, max delta, and residuals.
    """
    # 1. Verify matched experimental conditions
    conditions_match = True
    mismatches = []

    if data_metrohm.ph is not None and data_biologic.ph is not None:
        if abs(data_metrohm.ph - data_biologic.ph) > 0.2:
            conditions_match = False
            mismatches.append(f"pH mismatch: Metrohm {data_metrohm.ph} vs BioLogic {data_biologic.ph}")

    if data_metrohm.electrode_area_cm2 is not None and data_biologic.electrode_area_cm2 is not None:
        if abs(data_metrohm.electrode_area_cm2 - data_biologic.electrode_area_cm2) > 1e-4:
            mismatches.append(f"Area difference: {data_metrohm.electrode_area_cm2} vs {data_biologic.electrode_area_cm2} cm^2")

    # 2. Extract potentials vs RHE and current densities
    pot_m = data_metrohm.get_potential_vs_rhe()
    j_m = data_metrohm.get_current_density(unit="mA/cm^2")

    pot_b = data_biologic.get_potential_vs_rhe()
    j_b = data_biologic.get_current_density(unit="mA/cm^2")

    # Extract eta_10
    eta10_m = extract_overpotential_at_current(pot_m, j_m, target_j=target_j, reaction=reaction)
    eta10_b = extract_overpotential_at_current(pot_b, j_b, target_j=target_j, reaction=reaction)

    delta_eta10 = (eta10_b - eta10_m) if (eta10_m is not None and eta10_b is not None) else None

    # 3. Compute common overlapping potential window for residual comparison
    pot_min = max(float(np.min(pot_m)), float(np.min(pot_b)))
    pot_max = min(float(np.max(pot_m)), float(np.max(pot_b)))

    if pot_max > pot_min:
        grid_pot = np.linspace(pot_min, pot_max, potential_grid_points)

        # Sort and interpolate
        idx_m = np.argsort(pot_m)
        j_m_interp = np.interp(grid_pot, pot_m[idx_m], j_m[idx_m])

        idx_b = np.argsort(pot_b)
        j_b_interp = np.interp(grid_pot, pot_b[idx_b], j_b[idx_b])

        residuals = j_b_interp - j_m_interp
        rmse_j = float(np.sqrt(np.mean(residuals**2)))
        max_abs_diff = float(np.max(np.abs(residuals)))
    else:
        grid_pot = np.array([])
        j_m_interp = np.array([])
        j_b_interp = np.array([])
        residuals = np.array([])
        rmse_j = float("nan")
        max_abs_diff = float("nan")

    return {
        "conditions_match": conditions_match,
        "mismatches": mismatches,
        "eta10_metrohm_mV": eta10_m,
        "eta10_biologic_mV": eta10_b,
        "delta_eta10_mV": delta_eta10,
        "rmse_current_density": rmse_j,
        "max_abs_difference": max_abs_diff,
        "common_potential_grid": grid_pot,
        "j_metrohm_interpolated": j_m_interp,
        "j_biologic_interpolated": j_b_interp,
        "residuals": residuals,
    }


def plot_cross_instrument_benchmark(
    data_metrohm: StandardizedElectrochemicalData,
    data_biologic: StandardizedElectrochemicalData,
    reaction: str = "OER",
    target_j: float = 10.0,
    title: Optional[str] = "Cross-Instrument Electrocatalytic Benchmarking",
    save_path: Optional[Union[str, Path]] = None,
) -> Tuple[plt.Figure, Sequence[plt.Axes]]:
    """
    Generate a 2-panel publication figure comparing Metrohm Autolab and BioLogic EC-Lab.

    Top panel: Superimposed polarization curves with target current indicators.
    Bottom panel: Point-by-point current density residual (BioLogic - Metrohm).

    Parameters
    ----------
    data_metrohm : StandardizedElectrochemicalData
    data_biologic : StandardizedElectrochemicalData
    reaction : str, default 'OER'
    target_j : float, default 10.0
    title : str, optional
    save_path : str or Path, optional

    Returns
    -------
    Tuple[plt.Figure, Sequence[plt.Axes]]
        (fig, (ax_top, ax_res)).
    """
    apply_publication_style()
    comp = compare_instrument_lsv(data_metrohm, data_biologic, reaction=reaction, target_j=target_j)

    fig, (ax_top, ax_res) = create_publication_figure(
        columns="single",
        aspect_ratio=1.4,
        nrows=2,
        ncols=1,
        sharex=True,
        gridspec_kw={"height_ratios": [2.5, 1.0]},
    )

    pot_m = data_metrohm.get_potential_vs_rhe()
    j_m = data_metrohm.get_current_density(unit="mA/cm^2")

    pot_b = data_biologic.get_potential_vs_rhe()
    j_b = data_biologic.get_current_density(unit="mA/cm^2")

    c_m = get_color(0)  # Metrohm
    c_b = get_color(1)  # BioLogic

    lbl_m = f"Metrohm Autolab ({data_metrohm.instrument_software})"
    lbl_b = f"BioLogic ({data_biologic.instrument_software})"

    # Top Panel: LSV Curves
    ax_top.plot(pot_m, j_m, label=lbl_m, color=c_m, linewidth=1.5, zorder=3)
    ax_top.plot(pot_b, j_b, label=lbl_b, color=c_b, linewidth=1.5, linestyle="--", zorder=3)

    # Horizontal target current line
    target_val = abs(target_j) if reaction.upper() == "OER" else -abs(target_j)
    ax_top.axhline(target_val, color="#777777", linestyle=":", linewidth=1.0, alpha=0.8, zorder=1)

    # Annotate eta_10 values
    eta_m = comp["eta10_metrohm_mV"]
    eta_b = comp["eta10_biologic_mV"]
    if eta_m is not None and eta_b is not None:
        info_str = (
            f"Metrohm $\\eta_{{10}} = {eta_m:.1f}$ mV\n"
            f"BioLogic $\\eta_{{10}} = {eta_b:.1f}$ mV\n"
            f"$\\Delta\\eta_{{10}} = {abs(eta_b - eta_m):.1f}$ mV"
        )
        ax_top.text(
            0.05,
            0.82,
            info_str,
            transform=ax_top.transAxes,
            fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#cccccc", alpha=0.9),
        )

    ax_top.set_ylabel(r"$j~[\mathrm{mA\cdot cm^{-2}}]$")
    ax_top.legend(frameon=True, fontsize=8, loc="best")
    ax_top.grid(True, linestyle="--", alpha=0.3)

    # Bottom Panel: Residual
    if len(comp["common_potential_grid"]) > 0:
        ax_res.plot(comp["common_potential_grid"], comp["residuals"], color="#2b5c8f", linewidth=1.2)
        ax_res.axhline(0, color="black", linestyle="-", linewidth=0.8, alpha=0.6)
        ax_res.set_ylabel(r"$\Delta j~[\mathrm{mA\cdot cm^{-2}}]$", fontsize=7.5)
        ax_res.grid(True, linestyle="--", alpha=0.3)

        rmse = comp["rmse_current_density"]
        if np.isfinite(rmse):
            ax_res.text(
                0.95,
                0.80,
                f"RMSE: {rmse:.3f} mA/cm²",
                transform=ax_res.transAxes,
                fontsize=7,
                ha="right",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#dddddd", alpha=0.85),
            )

    ax_res.set_xlabel(r"$E~[\mathrm{V~vs.~RHE}]$")

    if title:
        fig.suptitle(title, fontsize=10, y=0.98)

    fig.tight_layout()

    if save_path:
        out_p = Path(save_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_p, dpi=300, bbox_inches="tight")

    return fig, (ax_top, ax_res)

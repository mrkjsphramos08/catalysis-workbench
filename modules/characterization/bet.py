"""
Brunauer-Emmett-Teller (BET) and Barrett-Joyner-Halenda (BJH) Gas Physisorption Analysis.

Implements:
- BET surface area linear transform with Rouquerol consistency criteria evaluation
- Monolayer capacity (q_m) and BET C constant calculation
- Specific surface area S_BET for N2 at 77.35 K (molecular area sigma = 0.162 nm^2)
- BJH mesopore size distribution (pore diameter D_p, incremental volume, dV/dD)
- Publication-quality isotherm, BET transform, and pore distribution plots.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

from parsers.characterization.common import CharacterizationData
from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette

# N2 physical constants at 77.35 K
N2_MOLAR_VOLUME_GAS_STP_CM3: float = 22414.0       # cm^3/mol
N2_LIQUID_MOLAR_VOLUME_CM3: float = 34.68           # cm^3/mol
N2_SURFACE_TENSION_MN_M: float = 8.72               # mN/m (dyn/cm)
N2_CROSS_SECTIONAL_AREA_NM2: float = 0.162          # nm^2 (16.2 A^2)
AVOGADRO_CONSTANT: float = 6.02214076e23
GAS_CONSTANT_J_MOL_K: float = 8.314462618


def calculate_bet_surface_area(
    p_p0: np.ndarray,
    quantity_adsorbed: np.ndarray,
    p_range: Tuple[float, float] = (0.05, 0.30),
    adsorbate: str = "N2",
    molecular_area_nm2: float = 0.162,
) -> Dict[str, Any]:
    """
    Calculate specific surface area (S_BET) using the linear BET equation:
    (p/p0) / [q * (1 - p/p0)] = (1 / (q_m * C)) + ((C - 1) / (q_m * C)) * (p/p0)

    Evaluates Rouquerol consistency criteria:
    1. q * (1 - p/p0) must continuously increase with p/p0.
    2. The BET constant C must be strictly positive (C > 0).
    3. The relative pressure corresponding to monolayer capacity, 1 / (sqrt(C) + 1),
       must fall within the selected fitting window.

    Parameters
    ----------
    p_p0 : np.ndarray
        Relative pressure array.
    quantity_adsorbed : np.ndarray
        Adsorbed volume in cm^3/g STP.
    p_range : tuple of (float, float), default (0.05, 0.30)
        Relative pressure fitting window [p_min, p_max].
    adsorbate : str, default 'N2'
        Adsorbate gas.
    molecular_area_nm2 : float, default 0.162
        Cross-sectional area of adsorbate molecule in nm^2 (16.2 Å^2 for N2).

    Returns
    -------
    Dict[str, Any]
        Dictionary of S_BET (m^2/g), C constant, q_m (cm^3/g STP), R^2, and Rouquerol checks.
    """
    p_arr = np.asarray(p_p0, dtype=float)
    q_arr = np.asarray(quantity_adsorbed, dtype=float)

    p1, p2 = sorted(p_range)
    mask = (p_arr >= p1) & (p_arr <= p2) & (p_arr < 1.0) & (q_arr > 0.0)

    if np.sum(mask) < 3:
        raise ValueError(
            f"Fewer than 3 valid adsorption data points in relative pressure window [{p1}, {p2}]."
        )

    p_fit = p_arr[mask]
    q_fit = q_arr[mask]

    # BET ordinate: y = (p/p0) / (q * (1 - p/p0))
    y_bet = p_fit / (q_fit * (1.0 - p_fit))

    slope, intercept, r_value, p_value, std_err = stats.linregress(p_fit, y_bet)

    # Monolayer capacity: q_m = 1 / (slope + intercept)
    q_m = 1.0 / (slope + intercept) if (slope + intercept) > 0 else np.nan
    # BET constant: C = (slope / intercept) + 1
    c_const = (slope / intercept) + 1.0 if intercept > 0 else np.nan

    # Specific surface area: S_BET = (q_m * N_A * sigma) / V_m
    # Convert sigma nm^2 -> m^2: 1 nm^2 = 1e-18 m^2
    area_factor = (AVOGADRO_CONSTANT * (molecular_area_nm2 * 1e-18)) / N2_MOLAR_VOLUME_GAS_STP_CM3
    s_bet = q_m * area_factor if np.isfinite(q_m) else np.nan

    # Rouquerol Consistency Checks
    # 1. q * (1 - p/p0) continuously increasing with p/p0
    q_rouq = q_fit * (1.0 - p_fit)
    rouq_diff = np.diff(q_rouq)
    rouq_monotonic = bool(np.all(rouq_diff > 0))

    # 2. C > 0
    c_positive = bool(c_const > 0) if np.isfinite(c_const) else False

    # 3. p_p0 at monolayer coverage: 1 / (sqrt(C) + 1) within fit range
    p_qm = 1.0 / (np.sqrt(c_const) + 1.0) if (c_positive and np.isfinite(c_const)) else np.nan
    p_qm_in_range = bool(p1 <= p_qm <= p2) if np.isfinite(p_qm) else False

    all_rouq_valid = rouq_monotonic and c_positive and p_qm_in_range

    return {
        "s_bet_m2_g": float(s_bet),
        "c_constant": float(c_const),
        "q_m_cm3_g_stp": float(q_m),
        "slope": float(slope),
        "intercept": float(intercept),
        "r_squared": float(r_value ** 2),
        "p_fit_range": (float(p1), float(p2)),
        "p_qm": float(p_qm) if np.isfinite(p_qm) else None,
        "rouquerol_criteria": {
            "monotonic_increasing": rouq_monotonic,
            "c_positive": c_positive,
            "qm_within_fit_range": p_qm_in_range,
            "overall_valid": all_rouq_valid,
        },
        "fit_p": p_fit,
        "fit_y": y_bet,
        "calc_y": slope * p_fit + intercept,
    }


def calculate_bjh_pore_distribution(
    p_p0: np.ndarray,
    quantity_adsorbed: np.ndarray,
    branch: str = "desorption",
) -> Dict[str, Any]:
    """
    Calculate mesopore size distribution using the Barrett-Joyner-Halenda (BJH) algorithm.

    Uses the cylindrical pore model with Halsey statistical thickness equation:
    t = 0.354 * [5.0 / -ln(p/p0)]^(1/3)  [nm]
    and the Kelvin capillary condensation equation:
    r_k = -0.953 / ln(p/p0)  [nm]
    Total pore radius: r_p = r_k + t, diameter D_p = 2 * r_p.

    Notes
    -----
    BJH assumes non-interacting cylindrical pores. Desorption branch data in materials with
    pore network connectivity may exhibit cavitation artifacts (~3.8 nm for N2 at 77 K).

    Parameters
    ----------
    p_p0 : np.ndarray
        Relative pressure array for the chosen branch.
    quantity_adsorbed : np.ndarray
        Adsorbed volume in cm^3/g STP.
    branch : str, default 'desorption'
        'desorption' or 'adsorption'.

    Returns
    -------
    Dict[str, Any]
        Pore diameters (nm), incremental pore volumes (cm^3/g), and dV/dD (cm^3/(g*nm)).
    """
    p_arr = np.asarray(p_p0, dtype=float)
    q_arr = np.asarray(quantity_adsorbed, dtype=float)

    # Mask valid mesoporous relative pressure range (0.35 <= p/p0 <= 0.99)
    mask = (p_arr >= 0.35) & (p_arr <= 0.99) & (q_arr > 0)
    if np.sum(mask) < 4:
        raise ValueError("Insufficient data points in mesoporous range (p/p0 between 0.35 and 0.99).")

    p_sub = p_arr[mask]
    q_sub = q_arr[mask]

    # Convert gas volume at STP (cm^3/g) to liquid volume of N2: V_liq = V_gas * 0.0015468 cm^3/cm^3
    v_liq = q_sub * (N2_LIQUID_MOLAR_VOLUME_CM3 / N2_MOLAR_VOLUME_GAS_STP_CM3)

    # Sort strictly descending for desorption, ascending for adsorption
    if branch.lower() == "desorption":
        sort_idx = np.argsort(p_sub)[::-1]
    else:
        sort_idx = np.argsort(p_sub)

    p_s = p_sub[sort_idx]
    v_s = v_liq[sort_idx]

    # Kelvin radius in nm: r_k = -2 * gamma * V_m / (R * T * ln(p/p0))
    # For N2 at 77.35 K: r_k = -0.953 / ln(p/p0) [nm]
    ln_p = np.log(p_s)
    r_k_nm = -0.953 / ln_p

    # Statistical thickness t using Halsey equation: t = 0.354 * (5 / -ln(p/p0))^(1/3) [nm]
    t_nm = 0.354 * ((5.0 / (-ln_p)) ** (1.0 / 3.0))

    # Total pore radius and diameter in nm
    r_p_nm = r_k_nm + t_nm
    d_p_nm = 2.0 * r_p_nm

    # Incremental pore volume between consecutive steps
    delta_v = np.abs(np.diff(v_s))
    delta_d = np.abs(np.diff(d_p_nm))

    mid_d_nm = (d_p_nm[:-1] + d_p_nm[1:]) / 2.0
    # Avoid zero division
    delta_d_safe = np.where(delta_d == 0, 1e-4, delta_d)
    dv_dd = delta_v / delta_d_safe

    return {
        "branch": branch,
        "pore_diameter_nm": mid_d_nm,
        "incremental_volume_cm3_g": delta_v,
        "cumulative_volume_cm3_g": np.cumsum(delta_v),
        "dv_dd_cm3_g_nm": dv_dd,
        "total_pore_volume_cm3_g": float(np.sum(delta_v)),
        "mean_pore_diameter_nm": float(np.average(mid_d_nm, weights=delta_v)) if np.sum(delta_v) > 0 else 0.0,
    }


def plot_isotherm(
    data: CharacterizationData,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot gas adsorption-desorption isotherm with distinct markers for branches.
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        created_fig = True
    else:
        fig = ax.get_figure()

    df = data.dataframe
    if df is not None and "branch" in df.columns:
        ads = df[df["branch"] == "adsorption"]
        des = df[df["branch"] == "desorption"]
        ax.plot(ads["p_p0"], ads["quantity_adsorbed"], marker="o", markersize=4, linestyle="-", color="#0072B2", label="Adsorption", zorder=3)
        if len(des) > 0:
            ax.plot(des["p_p0"], des["quantity_adsorbed"], marker="s", markersize=4, markerfacecolor="none", linestyle="--", color="#D55E00", label="Desorption", zorder=3)
    else:
        ax.plot(data.x, data.y, marker="o", markersize=4, color="#0072B2", label="Isotherm", zorder=3)

    ax.set_xlabel(r"Relative Pressure $p/p_0$")
    ax.set_ylabel(r"Quantity Adsorbed $[\mathrm{cm^3\cdot g^{-1}~STP}]$")
    if title:
        ax.set_title(title, fontsize=9, pad=6)

    ax.set_xlim(0, 1.0)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left", frameon=False, fontsize=8)

    if created_fig:
        fig.tight_layout()
    return fig, ax


def plot_bet_linear(
    bet_result: Dict[str, Any],
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot BET linear regression transform with annotations for S_BET, C constant, and R^2.
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        created_fig = True
    else:
        fig = ax.get_figure()

    p_fit = bet_result["fit_p"]
    y_fit = bet_result["fit_y"]
    calc_y = bet_result["calc_y"]

    ax.scatter(p_fit, y_fit, color="#0072B2", s=25, label="Data Points", zorder=3)
    ax.plot(p_fit, calc_y, color="#D55E00", linewidth=1.4, label="Linear Fit", zorder=2)

    s_bet = bet_result["s_bet_m2_g"]
    c_val = bet_result["c_constant"]
    r2 = bet_result["r_squared"]
    rouq_ok = bet_result["rouquerol_criteria"]["overall_valid"]

    status_tag = "Valid" if rouq_ok else "Caution (Rouquerol check failed)"
    annot = (
        rf"$S_{{\mathrm{{BET}}}} = {s_bet:.1f}~\mathrm{{m^2\cdot g^{{-1}}}}$"
        f"\n$C = {c_val:.1f}$"
        f"\n$R^2 = {r2:.5f}$"
        f"\nRouquerol: {status_tag}"
    )

    ax.text(
        0.05,
        0.75,
        annot,
        transform=ax.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.9),
    )

    ax.set_xlabel(r"Relative Pressure $p/p_0$")
    ax.set_ylabel(r"$\frac{p/p_0}{q(1 - p/p_0)}~[\mathrm{g\cdot cm^{-3}}]$")
    if title:
        ax.set_title(title, fontsize=9, pad=6)

    ax.legend(loc="lower right", frameon=False, fontsize=8)

    if created_fig:
        fig.tight_layout()
    return fig, ax


def plot_bjh_pore_distribution(
    bjh_result: Dict[str, Any],
    title: Optional[str] = None,
    color: str = "#009E73",
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot BJH pore size distribution (dV/dD vs D_p) on semi-logarithmic x-axis.
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        created_fig = True
    else:
        fig = ax.get_figure()

    d_p = bjh_result["pore_diameter_nm"]
    dv_dd = bjh_result["dv_dd_cm3_g_nm"]

    ax.plot(d_p, dv_dd, marker="o", markersize=3.5, color=color, linewidth=1.3, zorder=3)

    v_tot = bjh_result["total_pore_volume_cm3_g"]
    d_mean = bjh_result["mean_pore_diameter_nm"]
    branch_lbl = bjh_result.get("branch", "desorption").capitalize()

    annot = (
        f"Branch: {branch_lbl}"
        f"\nTotal Pore Vol: {v_tot:.3f} cm$^3$/g"
        f"\nMean Diameter: {d_mean:.1f} nm"
    )

    ax.text(
        0.95,
        0.90,
        annot,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.9),
    )

    ax.set_xscale("log")
    ax.set_xlabel(r"Pore Diameter $D_p~[\mathrm{nm}]$")
    ax.set_ylabel(r"$\mathrm{d}V / \mathrm{d}D_p~[\mathrm{cm^3\cdot g^{-1}\cdot nm^{-1}}]$")
    if title:
        ax.set_title(title, fontsize=9, pad=6)

    ax.set_ylim(bottom=0)

    if created_fig:
        fig.tight_layout()
    return fig, ax

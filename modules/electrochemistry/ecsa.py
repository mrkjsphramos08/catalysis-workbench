"""
Electrochemical Double-Layer Capacitance (Cdl) and ECSA Determination.

Evaluates double-layer capacitance (Cdl) from scan-rate series CVs in the non-Faradaic potential window:
    Delta j = j_a - j_c = 2 * C_dl * nu
Computes Electrochemical Surface Area (ECSA = C_dl * Area / C_s) and Roughness Factor (RF = C_dl / C_s).
Provides explicit scientific warnings regarding the invalidity of inferring active-site density directly from Cdl.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import warnings
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from parsers.electrochemistry.schema import StandardizedElectrochemicalData
from scripts.config import get_figure_dimensions
from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


class CdlInterpretationWarning(UserWarning):
    """Warns against conflating Cdl / ECSA with active-site density or intrinsic turnover frequency."""
    pass


def extract_cdl_from_scan_rates(
    cv_datasets: Sequence[StandardizedElectrochemicalData],
    potential_sample_V: float,
    specific_capacitance_mF_cm2: float = 0.040,
    area_cm2: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Extract double-layer capacitance (Cdl) from non-Faradaic CV scan-rate series.

    Parameters
    ----------
    cv_datasets : sequence of StandardizedElectrochemicalData
        CV datasets at distinct scan rates recorded within a purely non-Faradaic window.
    potential_sample_V : float
        Midpoint potential (V) within the non-Faradaic window where charging current is sampled.
    specific_capacitance_mF_cm2 : float, default 0.040 (40 uF/cm^2)
        Specific capacitance Cs of flat reference surface in mF/cm^2.
    area_cm2 : float, optional
        Electrode geometric area in cm^2.

    Returns
    -------
    Dict[str, Any]
        - 'cdl_mF_cm2': Double-layer capacitance in mF/cm^2
        - 'cdl_stderr': Fit standard error
        - 'r_squared': Fit R^2
        - 'ecsa_cm2': Estimated electrochemical surface area
        - 'roughness_factor': RF = ECSA / Area = Cdl / Cs
        - 'scan_rates_V_s': Array of scan rates (V/s)
        - 'delta_j_half_mA_cm2': Array of Delta j / 2 values (mA/cm^2)
    """
    warnings.warn(
        "Cdl and ECSA reflect electrostatically accessible electrochemical double-layer capacitance. "
        "Active-site density cannot be directly inferred from Cdl without independent selective probes "
        "(e.g., CO stripping, underpotential deposition, or redox titration).",
        CdlInterpretationWarning,
        stacklevel=2,
    )

    scan_rates = []
    delta_j_half = []

    for d in cv_datasets:
        if d.scan_rate_V_s is None or d.scan_rate_V_s <= 0:
            raise ValueError(f"Dataset from {d.filepath.name} lacks valid scan_rate_V_s.")

        pot = d.potential_V
        j_norm = d.get_current_density(area_cm2=area_cm2, unit="mA/cm2")

        # Use last cycle if multi-cycle
        if d.cycle_number is not None:
            last_c = d.unique_cycles[-1]
            mask = (d.cycle_number == last_c)
            pot = pot[mask]
            j_norm = j_norm[mask]

        # Find forward (anodic) and reverse (cathodic) current at sample potential
        # Separate half-cycles
        d_pot = np.diff(pot)
        forward_mask = np.append(d_pot >= 0, True)
        reverse_mask = np.append(d_pot < 0, True)

        pot_fwd = pot[forward_mask]
        j_fwd = j_norm[forward_mask]

        pot_rev = pot[reverse_mask]
        j_rev = j_norm[reverse_mask]

        if len(pot_fwd) < 2 or len(pot_rev) < 2:
            raise ValueError(f"Insufficient points in forward/reverse sweeps at {d.scan_rate_V_s} V/s.")

        # Interpolate currents at target potential
        # Ensure monotonic for interp
        idx_fwd = np.argsort(pot_fwd)
        j_anodic = float(np.interp(potential_sample_V, pot_fwd[idx_fwd], j_fwd[idx_fwd]))

        idx_rev = np.argsort(pot_rev)
        j_cathodic = float(np.interp(potential_sample_V, pot_rev[idx_rev], j_rev[idx_rev]))

        # Delta j / 2 = (j_a - j_c) / 2 = Cdl * nu
        dj_2 = (j_anodic - j_cathodic) / 2.0
        scan_rates.append(d.scan_rate_V_s)
        delta_j_half.append(dj_2)

    nu_arr = np.array(scan_rates, dtype=float)
    dj_arr = np.array(delta_j_half, dtype=float)

    # Linear fit: Delta j / 2 (mA/cm^2) = Cdl (mF/cm^2) * nu (V/s)
    # Note: 1 mA/cm^2 / (1 V/s) = 1 mF/cm^2 (exact SI unit alignment)
    res = stats.linregress(nu_arr, dj_arr)
    cdl_val = float(res.slope)
    r2 = float(res.rvalue ** 2)

    geom_area = area_cm2 or cv_datasets[0].electrode_area_cm2 or 1.0
    rf = float(cdl_val / specific_capacitance_mF_cm2)
    ecsa = float(rf * geom_area)

    return {
        "cdl_mF_cm2": cdl_val,
        "cdl_stderr": float(res.stderr),
        "r_squared": r2,
        "roughness_factor": rf,
        "ecsa_cm2": ecsa,
        "specific_capacitance_mF_cm2": specific_capacitance_mF_cm2,
        "scan_rates_V_s": nu_arr,
        "delta_j_half_mA_cm2": dj_arr,
        "potential_sample_V": potential_sample_V,
    }


def plot_cdl_analysis(
    cv_datasets: Sequence[StandardizedElectrochemicalData],
    potential_sample_V: float,
    specific_capacitance_mF_cm2: float = 0.040,
    area_cm2: Optional[float] = None,
    ax_cv: Optional[plt.Axes] = None,
    ax_fit: Optional[plt.Axes] = None,
    palette: str = "okabe_ito",
    title: Optional[str] = "ECSA / Cdl Determination",
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes], Dict[str, Any]]:
    """
    Generate dual-panel publication figure:
    (a) Non-Faradaic CV overlays at various scan rates.
    (b) Linear regression of capacitive current (Delta j / 2) vs scan rate to extract Cdl and ECSA.

    Parameters
    ----------
    cv_datasets : sequence of StandardizedElectrochemicalData
        CV series at different scan rates.
    potential_sample_V : float
        Midpoint potential (V) within non-Faradaic window.
    specific_capacitance_mF_cm2 : float, default 0.040 mF/cm^2
        Specific capacitance assumption.
    area_cm2 : float, optional
        Electrode geometric area in cm^2.
    ax_cv, ax_fit : plt.Axes, optional
        Subplot axes. If None, instantiates a double-panel figure.
    palette : str, default 'okabe_ito'
        Palette identifier.
    title : str, optional
        Overall title.

    Returns
    -------
    Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes], Dict[str, Any]]
        Figure, (ax_cv, ax_fit), and analysis result dictionary.
    """
    res = extract_cdl_from_scan_rates(
        cv_datasets=cv_datasets,
        potential_sample_V=potential_sample_V,
        specific_capacitance_mF_cm2=specific_capacitance_mF_cm2,
        area_cm2=area_cm2,
    )

    if ax_cv is None or ax_fit is None:
        fig, (ax_cv, ax_fit) = create_publication_figure(
            columns="double", aspect_ratio="2.2", nrows=1, ncols=2
        )
    else:
        fig = ax_cv.figure

    colors = get_palette(palette, max(len(cv_datasets), 2))

    # Panel (a): CV Overlays
    for idx, d in enumerate(cv_datasets):
        pot = d.potential_V
        j_norm = d.get_current_density(area_cm2=area_cm2, unit="mA/cm2")

        if d.cycle_number is not None:
            last_c = d.unique_cycles[-1]
            mask = (d.cycle_number == last_c)
            pot = pot[mask]
            j_norm = j_norm[mask]

        sr_mv = (d.scan_rate_V_s or 0.0) * 1000.0
        ax_cv.plot(pot, j_norm, color=colors[idx % len(colors)], linewidth=1.2, label=f"{sr_mv:.0f} mV/s")

    ax_cv.axvline(potential_sample_V, color="#D55E00", linestyle=":", linewidth=0.9, label=f"$E = {potential_sample_V:.2f}$ V")
    ax_cv.axhline(0.0, color="#666666", linestyle="--", linewidth=0.7)
    ax_cv.set_xlabel(r"Potential $E$ (V)")
    ax_cv.set_ylabel(r"Current Density $j$ ($\mathrm{mA\,cm^{-2}}$)")
    ax_cv.set_title("(a) Non-Faradaic CVs")
    ax_cv.legend(loc="upper left", framealpha=0.85, fontsize=6.8)

    # Panel (b): Cdl Linear Fit
    nu_mv_s = res["scan_rates_V_s"] * 1000.0  # mV/s for clean x-axis
    dj_vals = res["delta_j_half_mA_cm2"]
    cdl_val = res["cdl_mF_cm2"]
    r2 = res["r_squared"]
    rf = res["roughness_factor"]

    ax_fit.scatter(nu_mv_s, dj_vals, color="#0072B2", s=36, edgecolor="#222222", zorder=4, label="Sampled $\\Delta j / 2$")

    nu_fit_dense = np.linspace(0, max(nu_mv_s) * 1.05, 50)
    # y = Cdl * nu => (Cdl in mF/cm2) * (nu in V/s) = (Cdl) * (nu_mv / 1000)
    dj_fit_dense = cdl_val * (nu_fit_dense / 1000.0)

    ax_fit.plot(
        nu_fit_dense,
        dj_fit_dense,
        color="#D55E00",
        linestyle="-",
        linewidth=1.3,
        label=f"$C_{{\\mathrm{{dl}}}} = {cdl_val*1000:.1f}\\,\\mu\\mathrm{{F\\,cm^{{-2}}}}$\n$R^2 = {r2:.3f}$, $RF = {rf:.2f}$",
        zorder=3,
    )

    ax_fit.set_xlabel(r"Scan Rate $\nu$ ($\mathrm{mV\,s^{-1}}$)")
    ax_fit.set_ylabel(r"$\Delta j / 2$ ($\mathrm{mA\,cm^{-2}}$)")
    ax_fit.set_title("(b) Capacitive Current Scaling")
    ax_fit.legend(loc="upper left", framealpha=0.85, fontsize=6.8)

    return fig, (ax_cv, ax_fit), res


def calculate_double_layer_capacitance(
    scan_rates: Sequence[float],
    charging_currents: Sequence[float],
    scan_rate_unit: str = "mV/s",
    current_unit: str = "uA/cm2",
) -> Dict[str, Any]:
    """
    Fit charging current (Delta j / 2) vs scan rate (nu): Delta j / 2 = Cdl * nu.
    """
    sr = np.asarray(scan_rates, dtype=float)
    if scan_rate_unit.lower() in ("mv/s", "mv_s"):
        nu_v_s = sr * 1e-3
    else:
        nu_v_s = sr

    cur = np.asarray(charging_currents, dtype=float)
    if current_unit.lower() in ("ua/cm2", "ua/cm^2", "ua_cm2"):
        j_ma = cur * 1e-3
    elif current_unit.lower() in ("a/cm2", "a/cm^2", "a_cm2"):
        j_ma = cur * 1e3
    else:
        j_ma = cur

    res = stats.linregress(nu_v_s, j_ma)
    cdl_mF = float(res.slope)
    cdl_uF = cdl_mF * 1000.0

    return {
        "c_dl_mF_cm2": cdl_mF,
        "c_dl_uF_cm2": cdl_uF,
        "cdl_mF_cm2": cdl_mF,
        "cdl_uF_cm2": cdl_uF,
        "slope_stderr": float(res.stderr),
        "r_squared": float(res.rvalue ** 2),
        "intercept": float(res.intercept),
    }


def estimate_ecsa(
    cdl: float,
    cs_uF_cm2: float = 40.0,
    geometric_area_cm2: float = 1.0,
    warn_active_sites: bool = False,
) -> Dict[str, float]:
    """
    Estimate ECSA and roughness factor from double layer capacitance Cdl.
    """
    if warn_active_sites:
        warnings.warn(
            "Cdl and ECSA reflect electrostatically accessible electrochemical double-layer capacitance. "
            "Active-site density cannot be directly inferred from Cdl without independent selective probes "
            "(e.g., CO stripping, redox titrations, or in situ spectroscopic site counting).",
            CdlInterpretationWarning,
            stacklevel=2,
        )

    if cdl < 10.0 and cs_uF_cm2 > 10.0:
        cdl_uF = cdl * 1000.0
    else:
        cdl_uF = cdl

    rf = cdl_uF / float(cs_uF_cm2)
    ecsa = rf * float(geometric_area_cm2)

    return {
        "ecsa_cm2": ecsa,
        "roughness_factor": rf,
        "cdl_uF_cm2": cdl_uF,
        "cs_uF_cm2": cs_uF_cm2,
    }


plot_ecsa_analysis = plot_cdl_analysis

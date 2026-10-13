"""
Cyclic Voltammetry (CV) Analysis and Publication Visualization.

Supports multi-cycle plotting, scan-rate series overlays, redox peak detection
(anodic/cathodic peak potentials, peak currents, peak separations, formal potentials),
capacitive charging current integration, and publication figure generation.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

from parsers.electrochemistry.schema import StandardizedElectrochemicalData
from scripts.config import get_figure_dimensions
from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


def find_cv_peaks(
    potential: np.ndarray,
    current: np.ndarray,
    prominence: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Detect anodic and cathodic redox peaks from a single CV cycle.

    Parameters
    ----------
    potential : np.ndarray
        Potential array in Volts.
    current : np.ndarray
        Current or current density array.
    prominence : float, optional
        Required peak prominence. If None, auto-estimated from signal range.

    Returns
    -------
    Dict[str, Any]
        - 'e_pa': Anodic peak potential (V)
        - 'i_pa': Anodic peak current
        - 'e_pc': Cathodic peak potential (V)
        - 'i_pc': Cathodic peak current
        - 'delta_e_p': Peak separation |E_pa - E_pc| (V)
        - 'e_half': Formal / half-wave potential (E_pa + E_pc)/2 (V)
    """
    e_arr = np.asarray(potential, dtype=float)
    i_arr = np.asarray(current, dtype=float)

    i_range = np.ptp(i_arr)
    prom = prominence if prominence is not None else (i_range * 0.05 if i_range > 0 else 1e-6)

    # Anodic peaks (positive current maximum)
    peaks_anodic, _ = find_peaks(i_arr, prominence=prom)
    # Cathodic peaks (negative current minimum -> peak in -I)
    peaks_cathodic, _ = find_peaks(-i_arr, prominence=prom)

    e_pa = float(e_arr[peaks_anodic[np.argmax(i_arr[peaks_anodic])]]) if len(peaks_anodic) > 0 else None
    i_pa = float(np.max(i_arr[peaks_anodic])) if len(peaks_anodic) > 0 else None

    e_pc = float(e_arr[peaks_cathodic[np.argmin(i_arr[peaks_cathodic])]]) if len(peaks_cathodic) > 0 else None
    i_pc = float(np.min(i_arr[peaks_cathodic])) if len(peaks_cathodic) > 0 else None

    delta_ep = abs(e_pa - e_pc) if (e_pa is not None and e_pc is not None) else None
    e_half = (e_pa + e_pc) / 2.0 if (e_pa is not None and e_pc is not None) else None

    return {
        "e_pa": e_pa,
        "i_pa": i_pa,
        "e_pc": e_pc,
        "i_pc": i_pc,
        "delta_e_p": delta_ep,
        "e_half": e_half,
        "anodic_peak_indices": peaks_anodic.tolist(),
        "cathodic_peak_indices": peaks_cathodic.tolist(),
    }


def integrate_cv_charge(
    potential: np.ndarray,
    current_A: np.ndarray,
    scan_rate_V_s: float,
) -> Dict[str, float]:
    """
    Integrate total charge passed during a CV cycle: Q = (1/v) * integral(I dE).

    Parameters
    ----------
    potential : np.ndarray
        Potential array in Volts.
    current_A : np.ndarray
        Measured current in Amperes.
    scan_rate_V_s : float
        Scan rate in V/s (> 0).

    Returns
    -------
    Dict[str, float]
        - 'q_anodic_C': Total anodic charge in Coulombs (I > 0)
        - 'q_cathodic_C': Total cathodic charge in Coulombs (I < 0)
        - 'q_total_C': Net charge in Coulombs
    """
    if scan_rate_V_s <= 0:
        raise ValueError("Scan rate must be positive in V/s.")

    # Time step estimation from potential: dt = |dE| / v
    de = np.abs(np.diff(potential))
    dt = de / float(scan_rate_V_s)
    # Midpoint currents
    i_mid = (current_A[:-1] + current_A[1:]) / 2.0

    anodic_mask = i_mid > 0
    cathodic_mask = i_mid < 0

    q_anodic = float(np.sum(i_mid[anodic_mask] * dt[anodic_mask]))
    q_cathodic = float(np.sum(np.abs(i_mid[cathodic_mask]) * dt[cathodic_mask]))
    q_net = q_anodic - q_cathodic

    return {
        "q_anodic_C": q_anodic,
        "q_cathodic_C": q_cathodic,
        "q_total_C": q_net,
    }


def extract_redox_peaks(
    potential: np.ndarray,
    current: np.ndarray,
    prominence: Optional[float] = None,
) -> Dict[str, Any]:
    """Detect redox peaks with both abbreviated and descriptive metric keys."""
    res = find_cv_peaks(potential, current, prominence=prominence)
    res["anodic_peak_potential_V"] = res["e_pa"]
    res["cathodic_peak_potential_V"] = res["e_pc"]
    res["anodic_peak_current_A"] = res["i_pa"]
    res["cathodic_peak_current_A"] = res["i_pc"]
    res["delta_ep_mV"] = (res["delta_e_p"] * 1000.0) if res["delta_e_p"] is not None else None
    return res


def calculate_voltammetric_charge(
    time_or_potential: np.ndarray,
    current_A: np.ndarray,
    scan_rate_V_s: Optional[float] = None,
) -> Dict[str, float]:
    """
    Calculate voltammetric charge passed during a voltammetric cycle.
    Supports both (time_s, current_A) and (potential_V, current_A, scan_rate_V_s).
    """
    t_or_e = np.asarray(time_or_potential, dtype=float)
    i_arr = np.asarray(current_A, dtype=float)

    if scan_rate_V_s is not None and scan_rate_V_s > 0:
        res = integrate_cv_charge(t_or_e, i_arr, scan_rate_V_s)
        q_anodic = res["q_anodic_C"]
        q_cathodic = res["q_cathodic_C"]
        q_total = q_anodic + q_cathodic
        q_net = res["q_total_C"]
    else:
        dt = np.abs(np.diff(t_or_e))
        i_mid = (i_arr[:-1] + i_arr[1:]) / 2.0
        anodic_mask = i_mid > 0
        cathodic_mask = i_mid < 0
        q_anodic = float(np.sum(i_mid[anodic_mask] * dt[anodic_mask]))
        q_cathodic = float(np.sum(np.abs(i_mid[cathodic_mask]) * dt[cathodic_mask]))
        q_total = q_anodic + q_cathodic
        q_net = q_anodic - q_cathodic

    return {
        "q_anodic_C": q_anodic,
        "q_cathodic_C": q_cathodic,
        "q_total_C": q_net,
        "anodic_charge_C": q_anodic,
        "cathodic_charge_C": q_cathodic,
        "total_charge_C": q_total,
        "net_charge_C": q_net,
    }


def plot_cv(
    data: StandardizedElectrochemicalData,
    cycles: Optional[Sequence[int]] = None,
    current_unit: str = "mA/cm2",
    potential_scale: str = "RHE",
    reference_electrode: Optional[str] = None,
    ph: Optional[float] = None,
    ax: Optional[plt.Axes] = None,
    palette: str = "okabe_ito",
    mark_peaks: bool = False,
    title: Optional[str] = None,
) -> Tuple[plt.Axes, Dict[str, Any]]:
    """
    Plot Cyclic Voltammogram (CV) with multi-cycle support and publication styling.

    Parameters
    ----------
    data : StandardizedElectrochemicalData
        Parsed electrochemical dataset.
    cycles : sequence of int, optional
        Specific cycles to render (e.g. [1, 5, 10, 50]). If None, plots all available cycles.
    current_unit : str, default 'mA/cm2'
        'mA/cm2', 'A/cm2', 'mA', or 'A'.
    potential_scale : str, default 'RHE'
        'RHE' or 'measured'.
    reference_electrode : str, optional
        Reference electrode used if converting to RHE.
    ph : float, optional
        Electrolyte pH.
    ax : plt.Axes, optional
        Target Matplotlib Axes.
    palette : str, default 'okabe_ito'
        Palette identifier.
    mark_peaks : bool, default False
        Whether to highlight redox peaks on the primary cycle.
    title : str, optional
        Plot title.

    Returns
    -------
    Tuple[plt.Axes, Dict[str, Any]]
        Axes and analysis summary dictionary (peaks, scan rate, cycles).
    """
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    # Determine potential values
    if potential_scale.upper() == "RHE" and (data.ph is not None or ph is not None):
        pot = data.get_potential_vs_rhe(reference=reference_electrode, ph=ph)
        x_label = r"Potential $E$ (V vs RHE)"
    else:
        pot = data.potential_V
        ref_label = reference_electrode or data.potential_reference_raw or "Ref"
        x_label = f"Potential $E$ (V vs {ref_label})"

    # Determine current values
    u_clean = current_unit.strip().lower()
    if "cm" in u_clean or "area" in u_clean:
        cur = data.get_current_density(unit=current_unit)
        y_label = f"Current Density $j$ ({current_unit})"
    else:
        if u_clean in ("ma", "milliamp"):
            cur = data.current_A * 1000.0
            y_label = "Current $I$ (mA)"
        else:
            cur = data.current_A
            y_label = "Current $I$ (A)"

    all_cycles = data.unique_cycles
    target_cycles = [c for c in cycles if c in all_cycles] if cycles else all_cycles

    colors = get_palette(palette, max(len(target_cycles), 2))
    cycle_mask_base = data.cycle_number if data.cycle_number is not None else np.ones(len(pot), dtype=int)

    analysis_res: Dict[str, Any] = {
        "technique": "CV",
        "scan_rate_V_s": data.scan_rate_V_s,
        "cycles_plotted": target_cycles,
    }

    for idx, c_num in enumerate(target_cycles):
        c_mask = (cycle_mask_base == c_num)
        c_pot = pot[c_mask]
        c_cur = cur[c_mask]

        if len(c_pot) == 0:
            continue

        lbl = f"Cycle {c_num}" if len(target_cycles) > 1 else None
        ax.plot(c_pot, c_cur, color=colors[idx % len(colors)], linewidth=1.2, label=lbl, zorder=3)

        if mark_peaks and idx == len(target_cycles) - 1:
            pk = find_cv_peaks(c_pot, c_cur)
            analysis_res["peaks"] = pk
            if pk["e_pa"] is not None and pk["i_pa"] is not None:
                ax.scatter([pk["e_pa"]], [pk["i_pa"]], color="#D55E00", s=32, marker="^", edgecolor="#222222", zorder=5, label=f"Epa = {pk['e_pa']:.3f} V")
            if pk["e_pc"] is not None and pk["i_pc"] is not None:
                ax.scatter([pk["e_pc"]], [pk["i_pc"]], color="#0072B2", s=32, marker="v", edgecolor="#222222", zorder=5, label=f"Epc = {pk['e_pc']:.3f} V")

    ax.axhline(0.0, color="#666666", linestyle="--", linewidth=0.7, zorder=1)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    if title:
        ax.set_title(title)
    elif data.scan_rate_V_s is not None:
        sr_mv = data.scan_rate_V_s * 1000.0
        ax.set_title(f"CV ({sr_mv:.1f} mV/s)")

    if len(target_cycles) > 1 or mark_peaks:
        ax.legend(loc="best", framealpha=0.85, fontsize=6.8)

    return ax, analysis_res


def plot_cv_scan_rate_series(
    datasets: Sequence[StandardizedElectrochemicalData],
    labels: Optional[Sequence[str]] = None,
    current_unit: str = "mA/cm2",
    potential_scale: str = "RHE",
    reference_electrode: Optional[str] = None,
    ph: Optional[float] = None,
    ax: Optional[plt.Axes] = None,
    palette: str = "okabe_ito",
    title: Optional[str] = "Scan Rate Series CV",
) -> Tuple[plt.Axes, Dict[str, Any]]:
    """
    Overlay CV curves across multiple scan rates to visualize kinetics and non-Faradaic charging.
    """
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    colors = get_palette(palette, max(len(datasets), 2))
    rates_used = []

    for idx, d in enumerate(datasets):
        if potential_scale.upper() == "RHE":
            pot = d.get_potential_vs_rhe(reference=reference_electrode, ph=ph)
        else:
            pot = d.potential_V

        if "cm" in current_unit.lower():
            cur = d.get_current_density(unit=current_unit)
        else:
            cur = d.current_A * (1000.0 if "ma" in current_unit.lower() else 1.0)

        # Use last cycle if multi-cycle
        if d.cycle_number is not None:
            last_c = d.unique_cycles[-1]
            mask = (d.cycle_number == last_c)
            pot = pot[mask]
            cur = cur[mask]

        if labels and idx < len(labels):
            lbl = labels[idx]
        elif d.scan_rate_V_s is not None:
            lbl = f"{d.scan_rate_V_s * 1000.0:.0f} mV/s"
            rates_used.append(d.scan_rate_V_s)
        else:
            lbl = f"Curve {idx+1}"

        ax.plot(pot, cur, color=colors[idx % len(colors)], linewidth=1.2, label=lbl, zorder=3)

    ax.axhline(0.0, color="#666666", linestyle="--", linewidth=0.7, zorder=1)
    ax.set_xlabel(r"Potential $E$ (V vs RHE)" if potential_scale.upper() == "RHE" else "Potential $E$ (V)")
    ax.set_ylabel(f"Current Density $j$ ({current_unit})")
    if title:
        ax.set_title(title)
    ax.legend(loc="best", framealpha=0.85, fontsize=6.8)

    return ax, {"scan_rates_V_s": rates_used, "n_series": len(datasets)}

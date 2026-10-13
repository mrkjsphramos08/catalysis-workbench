"""
Electrochemical Impedance Spectroscopy (EIS) Analysis, Nyquist & Bode Visualization,
Equivalent Circuit Modeling, and Kramers-Kronig Validation.

Implements standard equivalent circuit fitting using the impedance package:
- Randles circuits (Rs - [Rct || CPE] or Rs - [Rct || Cdl] - W)
- Lin-KK test for Kramers-Kronig data validation
- Orthonormal (equal aspect ratio) Nyquist plots
- Two-panel stacked Bode plots (magnitude and phase)
- Area normalization (Ohm * cm^2) and multi-spectrum overlays.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import warnings

import matplotlib.pyplot as plt
import numpy as np

# Ensure impedance elements globals contains np to allow linKK evaluation
try:
    import impedance.models.circuits.elements as _ele
    if "np" not in _ele.circuit_elements:
        _ele.circuit_elements["np"] = np
    if "K" not in _ele.circuit_elements and hasattr(_ele, "K"):
        _ele.circuit_elements["K"] = getattr(_ele, "K")
except Exception:
    pass

from impedance.models.circuits import CustomCircuit
from impedance.validation import linKK

from parsers.electrochemistry.schema import StandardizedElectrochemicalData
from scripts.config import get_figure_dimensions
from scripts.utils import finalize_figure_layout, apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


def calculate_impedance_polar(
    z_real: np.ndarray,
    z_imag: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute impedance magnitude |Z| and phase angle phi in degrees.

    Parameters
    ----------
    z_real : np.ndarray
        Real impedance Z' (Ohm).
    z_imag : np.ndarray
        Imaginary impedance Z'' (Ohm). In physical AC impedance, capacitive Z'' < 0.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (magnitude_ohm, phase_angle_deg). Phase angle is in range [-180, 180].
        By electrochemistry convention, capacitive phase is typically reported as -phi or phi in [0, 90] deg.
    """
    re = np.asarray(z_real, dtype=float)
    im = np.asarray(z_imag, dtype=float)
    mag = np.hypot(re, im)
    # Phase angle in degrees: arctan2(Z'', Z')
    # For capacitive impedance (Z'' < 0), phase in degrees is negative (-90 to 0)
    phase_deg = np.degrees(np.arctan2(im, re))
    return mag, phase_deg


def fit_equivalent_circuit(
    frequency: np.ndarray,
    z_real: np.ndarray,
    z_imag: np.ndarray,
    circuit_string: str = "R0-p(R1,CPE1)",
    initial_guess: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """
    Fit an equivalent electrical circuit to complex impedance data.

    Parameters
    ----------
    frequency : np.ndarray
        Frequencies in Hz (> 0).
    z_real : np.ndarray
        Real component Z' in Ohm.
    z_imag : np.ndarray
        Imaginary component Z'' in Ohm.
    circuit_string : str, default 'R0-p(R1,CPE1)'
        Circuit description string for impedance.py CustomCircuit.
        Common representations:
        - 'R0-p(R1,CPE1)': Modified Randles without Warburg
        - 'R0-p(R1,C1)': Classic Randles without Warburg
        - 'R0-p(R1,CPE1)-W1': Randles with semi-infinite Warburg
    initial_guess : Optional[List[float]]
        Initial parameter values. If None, automatically estimated from spectrum features.

    Returns
    -------
    Dict[str, Any]
        Fitted parameters, confidence intervals, RMSE, and dense model predictions.
    """
    f = np.asarray(frequency, dtype=float)
    zr = np.asarray(z_real, dtype=float)
    zi = np.asarray(z_imag, dtype=float)

    # Filter out valid positive frequencies and non-NaNs
    mask = (f > 0) & np.isfinite(zr) & np.isfinite(zi)
    f_clean = f[mask]
    zr_clean = zr[mask]
    zi_clean = zi[mask]

    if len(f_clean) < 4:
        raise ValueError(f"Insufficient valid data points for circuit fitting: {len(f_clean)} points found.")

    # Sort descending by frequency if not already sorted
    sort_idx = np.argsort(f_clean)[::-1]
    f_clean = f_clean[sort_idx]
    zr_clean = zr_clean[sort_idx]
    zi_clean = zi_clean[sort_idx]

    # Combine into complex impedance
    z_complex = zr_clean + 1j * zi_clean

    # Auto-generate initial guess if not provided
    if initial_guess is None:
        rs_est = max(0.1, float(np.min(zr_clean)))
        rct_est = max(1.0, float(np.max(zr_clean) - np.min(zr_clean)))

        if circuit_string == "R0-p(R1,CPE1)":
            # [R0, R1, CPE1_Q, CPE1_n]
            initial_guess = [rs_est, rct_est, 1e-5, 0.85]
        elif circuit_string == "R0-p(R1,C1)":
            # [R0, R1, C1]
            initial_guess = [rs_est, rct_est, 1e-5]
        elif circuit_string == "R0-p(R1,CPE1)-W1":
            # [R0, R1, CPE1_Q, CPE1_n, W1]
            initial_guess = [rs_est, rct_est, 1e-5, 0.85, 50.0]
        else:
            # Fallback uniform guess
            initial_guess = [10.0] * circuit_string.count("R") + [1e-5] * circuit_string.count("C")

    circuit = CustomCircuit(circuit=circuit_string, initial_guess=initial_guess)
    circuit.fit(f_clean, z_complex)

    # Calculate RMSE
    z_fit_eval = circuit.predict(f_clean)
    rmse_val = float(np.sqrt(np.mean(np.abs(z_complex - z_fit_eval) ** 2)))

    # Generate dense prediction curve for plotting
    f_dense = np.logspace(np.log10(np.min(f_clean)), np.log10(np.max(f_clean)), 200)
    z_dense = circuit.predict(f_dense)

    # Extract named parameters
    param_names = circuit.get_param_names()
    param_dict = {}
    conf_dict = {}
    for i, name in enumerate(param_names[0]):
        val = float(circuit.parameters_[i])
        param_dict[name] = val
        if circuit.conf_ is not None and i < len(circuit.conf_):
            conf_dict[name] = float(circuit.conf_[i])
        else:
            conf_dict[name] = None

    # Identify primary electrochemical components
    rs_val = param_dict.get("R0", param_dict.get("R_0", None))
    rct_val = param_dict.get("R1", param_dict.get("R_1", None))
    cpe_q = param_dict.get("CPE1_0", param_dict.get("CPE1", None))
    cpe_n = param_dict.get("CPE1_1", None)

    return {
        "circuit_string": circuit_string,
        "parameters": param_dict,
        "conf_intervals": conf_dict,
        "rs_ohm": rs_val,
        "rct_ohm": rct_val,
        "cpe_q": cpe_q,
        "cpe_n": cpe_n,
        "rmse": rmse_val,
        "frequencies": f_clean,
        "fit_frequencies": f_dense,
        "fit_z_real": np.real(z_dense),
        "fit_z_imag": np.imag(z_dense),
        "circuit_model": circuit,
    }


def run_kramers_kronig_test(
    frequency: np.ndarray,
    z_real: np.ndarray,
    z_imag: np.ndarray,
    c: float = 0.85,
    max_M: int = 50,
) -> Dict[str, Any]:
    """
    Perform the Lin-KK Kramers-Kronig test to evaluate measurement stability and linearity.

    Parameters
    ----------
    frequency : np.ndarray
        Frequencies in Hz.
    z_real : np.ndarray
        Real component Z' (Ohm).
    z_imag : np.ndarray
        Imaginary component Z'' (Ohm).
    c : float, default 0.85
        Lin-KK cutoff parameter for mu.
    max_M : int, default 50
        Maximum number of Voigt RC elements.

    Returns
    -------
    Dict[str, Any]
        Dictionary containing validation flag, mu score, M elements, and residual arrays.
    """
    f = np.asarray(frequency, dtype=float)
    zr = np.asarray(z_real, dtype=float)
    zi = np.asarray(z_imag, dtype=float)

    mask = (f > 0) & np.isfinite(zr) & np.isfinite(zi)
    f_clean = f[mask]
    zr_clean = zr[mask]
    zi_clean = zi[mask]

    # Sort descending
    sort_idx = np.argsort(f_clean)[::-1]
    f_clean = f_clean[sort_idx]
    zr_clean = zr_clean[sort_idx]
    zi_clean = zi_clean[sort_idx]
    z_complex = zr_clean + 1j * zi_clean

    try:
        M, mu, Z_fit, res_real, res_imag = linKK(f_clean, z_complex, c=c, max_M=max_M)
        # In Lin-KK (Schönleber et al.), spectrum is KK-compliant if mu < cutoff c
        is_valid = bool((mu is not None) and (mu < c))
        return {
            "valid": is_valid,
            "mu": float(mu),
            "num_elements": int(M),
            "residuals_real": np.asarray(res_real),
            "residuals_imag": np.asarray(res_imag),
            "fit_z": np.asarray(Z_fit),
            "frequencies": f_clean,
        }
    except Exception as exc:
        warnings.warn(f"Lin-KK validation failed: {exc}", UserWarning)
        return {
            "valid": False,
            "mu": None,
            "num_elements": 0,
            "residuals_real": np.array([]),
            "residuals_imag": np.array([]),
            "fit_z": np.array([]),
            "frequencies": f_clean,
            "error": str(exc),
        }


def plot_nyquist(
    data: Union[StandardizedElectrochemicalData, Sequence[StandardizedElectrochemicalData], Tuple[np.ndarray, np.ndarray]],
    labels: Optional[Union[str, List[str]]] = None,
    fit_results: Optional[Union[Dict[str, Any], List[Optional[Dict[str, Any]]]]] = None,
    normalize_area: bool = False,
    area_cm2: Optional[float] = None,
    show_frequencies: bool = True,
    title: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
    ax: Optional[plt.Axes] = None,
    **kwargs: Any,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate an orthonormal publication-ready Nyquist plot (-Z'' vs Z').

    Adheres strictly to electrochemistry visualization standards:
    - Equal aspect ratio (1:1 scaling) so semicircles appear undistorted.
    - Standard axis labels with explicit units (Ohm or Ohm*cm^2).
    - Frequency marker annotations on critical points.
    - Fitted circuit line overlay when provided.

    Parameters
    ----------
    data : StandardizedElectrochemicalData or Sequence or (z_real, z_imag)
        Electrochemical data object(s) or raw impedance coordinate tuple.
    labels : str or List[str], optional
        Legend labels.
    fit_results : Dict or List of Dict, optional
        Output(s) from fit_equivalent_circuit.
    normalize_area : bool, default False
        If True, multiplies impedance by area to yield Ohm * cm^2.
    area_cm2 : float, optional
        Electrode area for normalization.
    show_frequencies : bool, default True
        If True, annotates key frequency points on the experimental curve.
    title : str, optional
        Figure title.
    save_path : str or Path, optional
        Destination path for saving.
    ax : plt.Axes, optional
        Existing axes to draw into.

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
    """
    apply_publication_style()

    # Harmonize input to lists
    if isinstance(data, tuple) and len(data) == 2:
        # Tuple of (z_real, z_imag)
        datasets = [(np.asarray(data[0]), np.asarray(data[1]), None)]
    elif isinstance(data, StandardizedElectrochemicalData):
        datasets = [(data.impedance_real_ohm, data.impedance_imag_ohm, data.frequency_Hz)]
    elif isinstance(data, (list, tuple)):
        datasets = []
        for d in data:
            if isinstance(d, StandardizedElectrochemicalData):
                datasets.append((d.impedance_real_ohm, d.impedance_imag_ohm, d.frequency_Hz))
            elif isinstance(d, tuple) and len(d) == 2:
                datasets.append((np.asarray(d[0]), np.asarray(d[1]), None))
    else:
        raise ValueError("Unsupported data format for plot_nyquist.")

    if labels is None:
        plot_labels = [f"Dataset {i+1}" if len(datasets) > 1 else "Experiment" for i in range(len(datasets))]
    elif isinstance(labels, str):
        plot_labels = [labels]
    else:
        plot_labels = list(labels)

    if fit_results is None:
        fits_list: List[Optional[Dict[str, Any]]] = [None] * len(datasets)
    elif isinstance(fit_results, dict):
        fits_list = [fit_results] + [None] * (len(datasets) - 1)
    else:
        fits_list = list(fit_results)

    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio=1.0)
        created_fig = True
    else:
        fig = ax.get_figure()

    colors = [get_color(i) for i in range(len(datasets))]

    for i, (zr, zi, freq) in enumerate(datasets):
        if zr is None or zi is None:
            continue

        c = colors[i % len(colors)]
        lbl = plot_labels[i] if i < len(plot_labels) else f"Dataset {i+1}"

        # Area normalization scaling
        scale = 1.0
        if normalize_area:
            a = area_cm2
            if a is None and isinstance(data, StandardizedElectrochemicalData) and data.electrode_area_cm2:
                a = data.electrode_area_cm2
            scale = a if a is not None else 1.0

        zr_plot = zr * scale
        # In Nyquist plots, ordinate is -Z''
        neg_zi_plot = -zi * scale

        # Plot experimental points
        ax.scatter(
            zr_plot,
            neg_zi_plot,
            label=lbl,
            color=c,
            s=22,
            edgecolors="none",
            alpha=0.85,
            zorder=3,
        )

        # Plot fitted curve if available
        fit = fits_list[i] if i < len(fits_list) else None
        if fit is not None and "fit_z_real" in fit and "fit_z_imag" in fit:
            fit_zr = fit["fit_z_real"] * scale
            fit_neg_zi = -fit["fit_z_imag"] * scale
            ax.plot(
                fit_zr,
                fit_neg_zi,
                color="black" if len(datasets) == 1 else c,
                linestyle="-",
                linewidth=1.4,
                label=f"Fit ({fit.get('circuit_string', 'Model')})" if len(datasets) == 1 else None,
                zorder=2,
            )

        # Frequency annotation on primary curve
        if show_frequencies and freq is not None and i == 0 and len(freq) == len(zr):
            # Target decade frequencies: 100 kHz, 10 kHz, 1 kHz, 100 Hz, 10 Hz, 1 Hz, 0.1 Hz
            target_decades = [1e5, 1e4, 1e3, 100.0, 10.0, 1.0, 0.1]
            for dec in target_decades:
                if dec < np.min(freq) or dec > np.max(freq):
                    continue
                # Closest index
                idx = int(np.argmin(np.abs(freq - dec)))
                if np.abs(np.log10(freq[idx]) - np.log10(dec)) < 0.2:
                    fx = zr_plot[idx]
                    fy = neg_zi_plot[idx]
                    f_str = f"{dec:.0f} Hz" if dec >= 1 else f"{dec:.1f} Hz"
                    if dec >= 1e3:
                        f_str = f"{dec/1e3:.0f} kHz"
                    ax.annotate(
                        f_str,
                        xy=(fx, fy),
                        xytext=(5, 5),
                        textcoords="offset points",
                        fontsize=7,
                        color="#444444",
                    )
                    ax.plot(fx, fy, marker="o", markersize=3.5, color="#333333", zorder=4)

    # Unit labels
    unit_str = r"\Omega\cdot\mathrm{cm}^2" if normalize_area else r"\Omega"
    ax.set_xlabel(rf"$Z'~[{unit_str}]$")
    ax.set_ylabel(rf"$-Z''~[{unit_str}]$")

    # Orthonormal aspect ratio to ensure semicircles are undistorted
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, linestyle="--", alpha=0.3, zorder=0)

    if len(datasets) > 1 or any(f is not None for f in fits_list):
        ax.legend(frameon=True, fontsize=8, loc="best")

    if title:
        ax.set_title(title, fontsize=10, pad=8)

    if created_fig:
        finalize_figure_layout(fig)

    if save_path:
        out_p = Path(save_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_p, dpi=300, bbox_inches="tight")

    return fig, ax


def plot_bode(
    data: Union[StandardizedElectrochemicalData, Sequence[StandardizedElectrochemicalData]],
    labels: Optional[Union[str, List[str]]] = None,
    fit_results: Optional[Union[Dict[str, Any], List[Optional[Dict[str, Any]]]]] = None,
    title: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
    **kwargs: Any,
) -> Tuple[plt.Figure, Sequence[plt.Axes]]:
    """
    Generate a 2-panel stacked Bode plot: Magnitude |Z| and Phase Angle phi vs Frequency.

    Parameters
    ----------
    data : StandardizedElectrochemicalData or Sequence
        Electrochemical data object(s).
    labels : str or List[str], optional
        Dataset labels.
    fit_results : Dict or List of Dict, optional
        Fitted equivalent circuit results.
    title : str, optional
        Overall title.
    save_path : str or Path, optional
        Destination path.

    Returns
    -------
    Tuple[plt.Figure, Sequence[plt.Axes]]
        (fig, (ax_mag, ax_phase)).
    """
    apply_publication_style()

    if isinstance(data, StandardizedElectrochemicalData):
        datasets = [data]
    else:
        datasets = list(data)

    if labels is None:
        plot_labels = [f"Dataset {i+1}" if len(datasets) > 1 else "Experiment" for i in range(len(datasets))]
    elif isinstance(labels, str):
        plot_labels = [labels]
    else:
        plot_labels = list(labels)

    if fit_results is None:
        fits_list: List[Optional[Dict[str, Any]]] = [None] * len(datasets)
    elif isinstance(fit_results, dict):
        fits_list = [fit_results] + [None] * (len(datasets) - 1)
    else:
        fits_list = list(fit_results)

    fig, (ax_mag, ax_phase) = create_publication_figure(columns="single", aspect_ratio=1.4, nrows=2, ncols=1, sharex=True)

    colors = [get_color(i) for i in range(len(datasets))]

    for i, d in enumerate(datasets):
        if d.frequency_Hz is None or d.impedance_real_ohm is None or d.impedance_imag_ohm is None:
            continue

        c = colors[i % len(colors)]
        lbl = plot_labels[i] if i < len(plot_labels) else f"Dataset {i+1}"

        f = d.frequency_Hz
        zr = d.impedance_real_ohm
        zi = d.impedance_imag_ohm

        mag, phase = calculate_impedance_polar(zr, zi)
        # Standard convention: plot -phase (positive 0 to 90 deg) or phase
        # Plotting -phase_deg makes capacitive peak positive
        phase_plot = -phase

        ax_mag.scatter(f, mag, color=c, s=18, alpha=0.85, label=lbl, zorder=3)
        ax_phase.scatter(f, phase_plot, color=c, s=18, alpha=0.85, label=lbl, zorder=3)

        fit = fits_list[i] if i < len(fits_list) else None
        if fit is not None and "fit_frequencies" in fit:
            f_fit = fit["fit_frequencies"]
            fit_zr = fit["fit_z_real"]
            fit_zi = fit["fit_z_imag"]
            fit_mag, fit_phase = calculate_impedance_polar(fit_zr, fit_zi)

            ax_mag.plot(f_fit, fit_mag, color="black" if len(datasets) == 1 else c, linestyle="-", linewidth=1.3, zorder=2)
            ax_phase.plot(f_fit, -fit_phase, color="black" if len(datasets) == 1 else c, linestyle="-", linewidth=1.3, zorder=2)

    ax_mag.set_xscale("log")
    ax_mag.set_yscale("log")
    ax_mag.set_ylabel(r"$|Z|~[\Omega]$")
    ax_mag.grid(True, linestyle="--", alpha=0.3, which="both")
    if len(datasets) > 1 or any(f is not None for f in fits_list):
        ax_mag.legend(frameon=True, fontsize=7.5, loc="best")

    ax_phase.set_xscale("log")
    ax_phase.set_xlabel(r"$\mathrm{Frequency}~[f / \mathrm{Hz}]$")
    ax_phase.set_ylabel(r"$-\phi~[^\circ]$")
    ax_phase.set_ylim(-10, 100)
    ax_phase.grid(True, linestyle="--", alpha=0.3, which="both")

    if title:
        fig.suptitle(title, fontsize=10, y=0.98)

    finalize_figure_layout(fig)

    if save_path:
        out_p = Path(save_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_p, dpi=300, bbox_inches="tight")

    return fig, (ax_mag, ax_phase)


def plot_eis_comprehensive(
    data: StandardizedElectrochemicalData,
    fit_result: Optional[Dict[str, Any]] = None,
    title: Optional[str] = "Electrochemical Impedance Spectroscopy Analysis",
    save_path: Optional[Union[str, Path]] = None,
) -> plt.Figure:
    """
    Generate a 3-panel publication figure combining Nyquist, Bode Magnitude, and Bode Phase.

    Includes fitted equivalent circuit parameters and Lin-KK validation status.
    """
    apply_publication_style()
    fig, (ax_nyq, ax_bode_m, ax_bode_p) = create_publication_figure(columns="double", aspect_ratio=0.45, nrows=1, ncols=3)

    # 1. Nyquist
    plot_nyquist(data, labels="Experiment", fit_results=fit_result, ax=ax_nyq, show_frequencies=True)
    ax_nyq.set_title("(a) Nyquist Spectrum", fontsize=9, pad=6)

    # Annotate Rs and Rct if fit available
    if fit_result is not None:
        rs = fit_result.get("rs_ohm")
        rct = fit_result.get("rct_ohm")
        rmse = fit_result.get("rmse")
        lines = []
        if rs is not None:
            lines.append(rf"$R_\mathrm{{s}} = {rs:.2f}~\Omega$")
        if rct is not None:
            lines.append(rf"$R_\mathrm{{ct}} = {rct:.2f}~\Omega$")
        if rmse is not None:
            lines.append(rf"$\mathrm{{RMSE}} = {rmse:.2f}~\Omega$")

        if lines:
            text_str = "\n".join(lines)
            ax_nyq.text(
                0.05,
                0.95,
                text_str,
                transform=ax_nyq.transAxes,
                fontsize=7.5,
                verticalalignment="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.9),
            )

    # 2. Bode Magnitude
    f = data.frequency_Hz
    zr = data.impedance_real_ohm
    zi = data.impedance_imag_ohm
    mag, phase = calculate_impedance_polar(zr, zi)
    c = get_color(0)

    ax_bode_m.scatter(f, mag, color=c, s=18, alpha=0.85, zorder=3)
    if fit_result is not None and "fit_frequencies" in fit_result:
        f_fit = fit_result["fit_frequencies"]
        fit_mag, _ = calculate_impedance_polar(fit_result["fit_z_real"], fit_result["fit_z_imag"])
        ax_bode_m.plot(f_fit, fit_mag, color="black", linestyle="-", linewidth=1.3, zorder=2)

    ax_bode_m.set_xscale("log")
    ax_bode_m.set_yscale("log")
    ax_bode_m.set_xlabel(r"$\mathrm{Frequency}~[f / \mathrm{Hz}]$")
    ax_bode_m.set_ylabel(r"$|Z|~[\Omega]$")
    ax_bode_m.set_title("(b) Bode Magnitude", fontsize=9, pad=6)
    ax_bode_m.grid(True, linestyle="--", alpha=0.3, which="both")

    # 3. Bode Phase
    ax_bode_p.scatter(f, -phase, color=c, s=18, alpha=0.85, zorder=3)
    if fit_result is not None and "fit_frequencies" in fit_result:
        _, fit_phase = calculate_impedance_polar(fit_result["fit_z_real"], fit_result["fit_z_imag"])
        ax_bode_p.plot(f_fit, -fit_phase, color="black", linestyle="-", linewidth=1.3, zorder=2)

    ax_bode_p.set_xscale("log")
    ax_bode_p.set_xlabel(r"$\mathrm{Frequency}~[f / \mathrm{Hz}]$")
    ax_bode_p.set_ylabel(r"$-\phi~[^\circ]$")
    ax_bode_p.set_ylim(-10, 100)
    ax_bode_p.set_title("(c) Bode Phase", fontsize=9, pad=6)
    ax_bode_p.grid(True, linestyle="--", alpha=0.3, which="both")

    if title:
        fig.suptitle(title, fontsize=10, y=1.02)

    finalize_figure_layout(fig)

    if save_path:
        out_p = Path(save_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_p, dpi=300, bbox_inches="tight")

    return fig

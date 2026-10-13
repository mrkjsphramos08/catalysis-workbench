"""
X-Ray Diffraction (XRD) Quantitative Analysis and Publication Visualization.

Implements:
- Bragg d-spacing calculation from 2theta and wavelength
- Scherrer crystallite size estimation with instrumental broadening correction
- Multi-peak deconvolution (Gaussian, Lorentzian, pseudo-Voigt) via lmfit
- Background subtraction (AsLS, SNIP, ModPoly)
- Reference diffraction stick patterns via pymatgen
- Publication-quality plotting (individual, stacked offsets, fit & residuals, reference comparison).
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import warnings
import matplotlib.pyplot as plt
import numpy as np

from parsers.characterization.common import CharacterizationData
from modules.characterization.common_peaks import (
    detect_peaks,
    estimate_baseline,
    fit_peaks,
    subtract_baseline,
)
from scripts.config import get_figure_dimensions
from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


class ScherrerBroadeningWarning(UserWarning):
    """Raised when peak broadening is smaller than or equal to instrumental resolution."""
    pass


def calculate_d_spacing(
    two_theta_deg: Union[float, np.ndarray],
    wavelength_angstrom: float = 1.54060,
) -> Union[float, np.ndarray]:
    """
    Calculate interplanar spacing (d-spacing) using Bragg's Law:
    n * lambda = 2 * d * sin(theta)  =>  d = lambda / (2 * sin(theta))

    Parameters
    ----------
    two_theta_deg : float or np.ndarray
        Diffraction angle 2theta in degrees.
    wavelength_angstrom : float, default 1.54060
        Incident X-ray radiation wavelength in Angstroms (default Cu K-alpha).

    Returns
    -------
    float or np.ndarray
        Interplanar spacing d in Angstroms.
    """
    two_theta = np.asarray(two_theta_deg, dtype=float)
    # Convert theta = (2theta / 2) to radians
    theta_rad = np.radians(two_theta / 2.0)
    sin_theta = np.sin(theta_rad)

    # Guard against non-physical division by zero
    with np.errstate(divide="ignore", invalid="ignore"):
        d_val = wavelength_angstrom / (2.0 * sin_theta)

    if np.ndim(two_theta_deg) == 0:
        return float(d_val)
    return d_val


def calculate_crystallite_size_scherrer(
    two_theta_deg: float,
    fwhm_deg: float,
    wavelength_angstrom: float = 1.54060,
    k_shape_factor: float = 0.9,
    instrumental_fwhm_deg: float = 0.05,
) -> Dict[str, Any]:
    """
    Estimate volume-weighted mean crystallite domain size using the Scherrer equation:
    D = (K * lambda) / (beta * cos(theta))

    Notes
    -----
    1. D reflects coherently scattering crystallite domain size, NOT physical particle or grain size.
    2. Peak broadening beta must be corrected for instrumental broadening.
    3. Assumes Gaussian/pseudo-Voigt instrumental broadening correction: beta^2 = beta_obs^2 - beta_inst^2.

    Parameters
    ----------
    two_theta_deg : float
        Bragg peak position 2theta in degrees.
    fwhm_deg : float
        Observed Full Width at Half Maximum (FWHM) in degrees.
    wavelength_angstrom : float, default 1.54060
        X-ray wavelength in Angstroms.
    k_shape_factor : float, default 0.9
        Scherrer shape factor K (typically 0.89 to 0.94 for spherical/cubic crystallites).
    instrumental_fwhm_deg : float, default 0.05
        Instrumental broadening contribution at this 2theta in degrees.

    Returns
    -------
    Dict[str, Any]
        - 'crystallite_size_nm': Apparent crystallite domain size in nanometers (None if unresolvable).
        - 'crystallite_size_angstrom': Size in Angstroms.
        - 'beta_corrected_rad': Corrected physical broadening beta in radians.
        - 'fwhm_obs_deg': Input observed FWHM in degrees.
        - 'instrumental_fwhm_deg': Instrumental broadening in degrees.
        - 'two_theta_deg': Peak position in degrees.
        - 'k_shape_factor': Shape factor used.
        - 'notes': Methodological caveats.
    """
    if fwhm_deg <= instrumental_fwhm_deg:
        warnings.warn(
            f"Observed FWHM ({fwhm_deg:.4f} deg) <= instrumental broadening ({instrumental_fwhm_deg:.4f} deg). "
            "Broadening is instrument-limited; crystallite size exceeds reliable Scherrer resolution limit (~100-200 nm).",
            ScherrerBroadeningWarning,
        )
        return {
            "crystallite_size_nm": None,
            "crystallite_size_angstrom": None,
            "beta_corrected_rad": None,
            "fwhm_obs_deg": fwhm_deg,
            "instrumental_fwhm_deg": instrumental_fwhm_deg,
            "two_theta_deg": two_theta_deg,
            "k_shape_factor": k_shape_factor,
            "notes": "Instrument-limited broadening; crystallite size exceeds upper Scherrer limit.",
        }

    # Convert observed and instrumental broadening to radians
    beta_obs_rad = np.radians(fwhm_deg)
    beta_inst_rad = np.radians(instrumental_fwhm_deg)

    # Correct for instrumental broadening (root-sum-square subtraction)
    beta_corr_rad = np.sqrt(beta_obs_rad**2 - beta_inst_rad**2)

    theta_rad = np.radians(two_theta_deg / 2.0)
    cos_theta = np.cos(theta_rad)

    # D in Angstroms
    d_angstrom = (k_shape_factor * wavelength_angstrom) / (beta_corr_rad * cos_theta)
    d_nm = d_angstrom / 10.0

    return {
        "crystallite_size_nm": float(d_nm),
        "crystallite_size_angstrom": float(d_angstrom),
        "beta_corrected_rad": float(beta_corr_rad),
        "fwhm_obs_deg": float(fwhm_deg),
        "instrumental_fwhm_deg": float(instrumental_fwhm_deg),
        "two_theta_deg": float(two_theta_deg),
        "k_shape_factor": float(k_shape_factor),
        "notes": "Scherrer domain size reflects coherently diffracting domain dimension, not physical grain/particle size.",
    }


def analyze_xrd_peaks(
    data: CharacterizationData,
    prominence: Optional[float] = None,
    distance: Optional[int] = 15,
    instrumental_fwhm_deg: float = 0.05,
    k_shape_factor: float = 0.9,
) -> List[Dict[str, Any]]:
    """
    Detect major diffraction peaks and compute Bragg d-spacing and Scherrer crystallite size.

    Parameters
    ----------
    data : CharacterizationData
        Standardized XRD dataset.
    prominence : float, optional
        Peak finding prominence threshold.
    distance : int, default 15
        Minimum index separation between peaks.
    instrumental_fwhm_deg : float, default 0.05
        Instrumental broadening.
    k_shape_factor : float, default 0.9
        Scherrer constant.

    Returns
    -------
    List[Dict[str, Any]]
        List of dictionaries with peak parameters (2theta, d-spacing, FWHM, crystallite size).
    """
    wl = data.metadata.get("wavelength_angstrom", 1.54060)
    pk_info = detect_peaks(data.x, data.y, prominence=prominence, distance=distance)

    results = []
    x_arr = data.x
    y_arr = data.y

    for idx, pk_2th in enumerate(pk_info["x_peaks"]):
        pk_idx = pk_info["indices"][idx]
        pk_int = pk_info["y_peaks"][idx]
        d_space = calculate_d_spacing(pk_2th, wavelength_angstrom=wl)

        # Estimate local FWHM around peak
        half_max = pk_int / 2.0
        # Search left and right for half-max
        left_idx = pk_idx
        while left_idx > 0 and y_arr[left_idx] > half_max:
            left_idx -= 1
        right_idx = pk_idx
        while right_idx < len(y_arr) - 1 and y_arr[right_idx] > half_max:
            right_idx += 1

        fwhm_approx = abs(x_arr[right_idx] - x_arr[left_idx]) if right_idx > left_idx else 0.2

        scherrer_res = calculate_crystallite_size_scherrer(
            two_theta_deg=pk_2th,
            fwhm_deg=fwhm_approx,
            wavelength_angstrom=wl,
            k_shape_factor=k_shape_factor,
            instrumental_fwhm_deg=instrumental_fwhm_deg,
        )

        results.append({
            "peak_index": idx + 1,
            "two_theta_deg": float(pk_2th),
            "intensity": float(pk_int),
            "d_spacing_angstrom": float(d_space),
            "fwhm_approx_deg": float(fwhm_approx),
            "crystallite_size_nm": scherrer_res["crystallite_size_nm"],
        })

    return results


def plot_xrd_pattern(
    data: CharacterizationData,
    annotate_peaks: bool = True,
    show_d_spacing: bool = True,
    subtract_bg: bool = False,
    bg_method: str = "asls",
    title: Optional[str] = None,
    color: str = "#0072B2",
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot single XRD diffractogram with optional peak annotations and d-spacing tags.
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        created_fig = True
    else:
        fig = ax.get_figure()

    x_plot = data.x
    y_plot = data.y

    if subtract_bg:
        y_plot, bg = subtract_baseline(x_plot, y_plot, method=bg_method)

    ax.plot(x_plot, y_plot, color=color, linewidth=1.2, label=data.sample_id, zorder=3)

    if annotate_peaks:
        peaks = analyze_xrd_peaks(data)
        # Annotate top 5 most intense peaks
        sorted_peaks = sorted(peaks, key=lambda p: p["intensity"], reverse=True)[:5]
        for p in sorted_peaks:
            th = p["two_theta_deg"]
            intensity = p["intensity"]
            d_tag = f"\n({p['d_spacing_angstrom']:.2f} Å)" if show_d_spacing else ""
            label_text = f"{th:.1f}°{d_tag}"
            ax.annotate(
                label_text,
                xy=(th, intensity),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
                fontsize=7,
                rotation=0,
                arrowprops=dict(arrowstyle="-", color="#555555", lw=0.6),
            )

    ax.set_xlabel(r"Diffraction Angle $2\theta~[^\circ]$")
    ax.set_ylabel(r"Intensity $[\mathrm{a.u.}]$")
    if title:
        ax.set_title(title, fontsize=9, pad=6)

    ax.set_xlim(float(np.min(x_plot)), float(np.max(x_plot)))
    y_max = float(np.max(y_plot))
    ax.set_ylim(0, y_max * 1.15)

    if created_fig:
        fig.tight_layout()
    return fig, ax


def plot_xrd_stacked(
    datasets: Sequence[CharacterizationData],
    labels: Optional[Sequence[str]] = None,
    normalize: bool = True,
    vertical_offset_fraction: float = 0.35,
    title: Optional[str] = None,
    palette: str = "okabe_ito",
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate stacked XRD diffractogram comparison with constant vertical offset.
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio=1.2)
        created_fig = True
    else:
        fig = ax.get_figure()

    colors = get_palette(palette, max(len(datasets), 2))
    names = labels if labels is not None else [d.sample_id for d in datasets]

    current_offset = 0.0
    for i, d in enumerate(datasets):
        c = colors[i % len(colors)]
        lbl = names[i]

        y_sig = d.y.copy()
        if normalize:
            y_min = float(np.min(y_sig))
            y_max = float(np.max(y_sig))
            y_norm = (y_sig - y_min) / (y_max - y_min) if (y_max - y_min) > 0 else y_sig
            step = vertical_offset_fraction
        else:
            y_norm = y_sig
            step = float(np.max(y_sig)) * vertical_offset_fraction

        ax.plot(d.x, y_norm + current_offset, color=c, linewidth=1.2, label=lbl, zorder=3)
        current_offset += step

    ax.set_xlabel(r"Diffraction Angle $2\theta~[^\circ]$")
    ax.set_ylabel(r"Intensity $[\mathrm{a.u.}]$")
    if title:
        ax.set_title(title, fontsize=9, pad=6)

    ax.legend(loc="upper right", frameon=False, fontsize=8)
    if created_fig:
        fig.tight_layout()
    return fig, ax


def plot_xrd_fit(
    data: CharacterizationData,
    fit_result: Dict[str, Any],
    region: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    """
    Generate publication 2-panel XRD deconvolution plot (data, model components, and residuals).
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

    x = data.x
    y = data.y
    fit_y = fit_result["fit_y"]
    res = fit_result["residuals"]

    # Top panel: Data & fit components
    ax_top.scatter(x, y, color="#444444", s=10, alpha=0.6, label="Observed", zorder=2)
    ax_top.plot(x, fit_y, color="#D55E00", linewidth=1.5, label="Total Fit", zorder=4)

    # Component profiles
    comps = fit_result.get("components", {})
    comp_colors = get_palette("tol_muted", max(len(comps), 2))
    c_idx = 0
    for name, curve in comps.items():
        if name.startswith("bg_") or name == "linear":
            ax_top.plot(x, curve, color="#999999", linestyle="--", linewidth=1.0, label="Background", zorder=3)
        else:
            c = comp_colors[c_idx % len(comp_colors)]
            c_idx += 1
            ax_top.plot(x, curve, color=c, linestyle="-", linewidth=1.1, alpha=0.85, label=f"Peak {c_idx}", zorder=3)

    r2 = fit_result.get("r_squared", 1.0)
    ax_top.text(
        0.04,
        0.90,
        rf"$R^2 = {r2:.4f}$",
        transform=ax_top.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.9),
    )

    ax_top.set_ylabel(r"Intensity $[\mathrm{a.u.}]$")
    if title:
        ax_top.set_title(title, fontsize=9, pad=6)
    ax_top.legend(loc="upper right", frameon=False, fontsize=7)

    # Bottom panel: Residuals
    ax_bot.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax_bot.scatter(x, res, color="#0072B2", s=8, alpha=0.7, zorder=2)
    ax_bot.set_xlabel(r"Diffraction Angle $2\theta~[^\circ]$")
    ax_bot.set_ylabel(r"$\Delta I$")

    if region is not None:
        ax_top.set_xlim(region)
        ax_bot.set_xlim(region)

    fig.tight_layout()
    return fig, (ax_top, ax_bot)

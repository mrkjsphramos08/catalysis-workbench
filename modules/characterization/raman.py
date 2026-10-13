"""
Raman Spectroscopy Quantitative Deconvolution and Publication Visualization.

Implements:
- Cosmic-ray despiking and baseline subtraction (AsLS, AirPLS, ModPoly)
- Multi-peak deconvolution using lmfit (Lorentzian, Gaussian, pseudo-Voigt)
- Rigorous distinction between peak height ratios (I1 / I2) and integrated area ratios (A1 / A2)
- Publication-quality spectral plots with residual diagnostic panels.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np

from parsers.characterization.common import CharacterizationData
from modules.characterization.common_peaks import (
    detect_peaks,
    fit_peaks,
    remove_cosmic_rays,
    smooth_spectrum,
    subtract_baseline,
)
from scripts.utils import finalize_figure_layout, apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


def process_raman_spectrum(
    data: CharacterizationData,
    despike: bool = True,
    subtract_bg: bool = True,
    bg_method: str = "asls",
    smooth: bool = False,
    smooth_window: int = 7,
) -> CharacterizationData:
    """
    Standardize, despike, baseline-correct, and optionally smooth a Raman spectrum.
    """
    out = data.copy()
    y_proc = out.y.copy()

    if despike:
        y_proc, spikes = remove_cosmic_rays(y_proc)
        if len(spikes) > 0:
            out.add_processing_step("despike", {"n_spikes": len(spikes)})

    if subtract_bg:
        y_proc, base = subtract_baseline(out.x, y_proc, method=bg_method)
        out.metadata["estimated_baseline"] = base
        out.add_processing_step("subtract_baseline", {"method": bg_method})

    if smooth:
        y_proc = smooth_spectrum(y_proc, window_length=smooth_window)
        out.add_processing_step("smooth", {"window": smooth_window})

    out.y = y_proc
    return out


def calculate_peak_ratios(
    fit_result: Dict[str, Any],
    peak_1_index: int = 1,
    peak_2_index: int = 2,
    label_1: str = "Peak 1",
    label_2: str = "Peak 2",
) -> Dict[str, Any]:
    """
    Calculate peak height and integrated area ratios from a deconvolution result.

    Notes
    -----
    Strictly distinguishes peak height ratios (I1 / I2) from integrated area ratios (A1 / A2).
    """
    peaks = fit_result.get("peaks", [])
    if len(peaks) < max(peak_1_index, peak_2_index):
        raise ValueError(
            f"Fit result contains only {len(peaks)} peaks, cannot extract index {max(peak_1_index, peak_2_index)}."
        )

    p1 = peaks[peak_1_index - 1]
    p2 = peaks[peak_2_index - 1]

    h1, h2 = p1["height"], p2["height"]
    a1, a2 = p1["amplitude"], p2["amplitude"]

    height_ratio = h1 / h2 if h2 > 0 else np.nan
    area_ratio = a1 / a2 if a2 > 0 else np.nan

    return {
        "label_1": label_1,
        "label_2": label_2,
        "peak_1_center": p1["center"],
        "peak_2_center": p2["center"],
        "height_1": float(h1),
        "height_2": float(h2),
        "height_ratio": float(height_ratio),
        "area_1": float(a1),
        "area_2": float(a2),
        "area_ratio": float(area_ratio),
        "disclaimer": "Fitted Raman peak ratios reflect relative vibrational mode intensities and must not be used to infer oxidation states or defect concentrations without calibrated reference standards.",
    }


def plot_raman_spectrum(
    data: CharacterizationData,
    show_baseline: bool = False,
    title: Optional[str] = None,
    color: str = "#0072B2",
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot single Raman spectrum with optional baseline curve.
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        created_fig = True
    else:
        fig = ax.get_figure()

    ax.plot(data.x, data.y, color=color, linewidth=1.2, label=data.sample_id, zorder=3)

    if show_baseline and "estimated_baseline" in data.metadata:
        base = data.metadata["estimated_baseline"]
        ax.plot(data.x, base, color="#999999", linestyle="--", linewidth=1.0, label="Baseline", zorder=2)

    ax.set_xlabel(r"Raman Shift $[\mathrm{cm}^{-1}]$")
    ax.set_ylabel(r"Intensity $[\mathrm{counts}]$")
    if title:
        ax.set_title(title, fontsize=9, pad=6)

    ax.set_xlim(float(np.min(data.x)), float(np.max(data.x)))
    ax.set_ylim(bottom=0)

    if created_fig:
        finalize_figure_layout(fig)
    return fig, ax


def plot_raman_stacked(
    datasets: Sequence[CharacterizationData],
    labels: Optional[Sequence[str]] = None,
    normalize: bool = True,
    vertical_offset_fraction: float = 0.4,
    title: Optional[str] = None,
    palette: str = "okabe_ito",
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot stacked Raman spectra with constant vertical offset.
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio=1.25)
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
            step = float(np.ptp(y_sig)) * vertical_offset_fraction

        ax.plot(d.x, y_norm + current_offset, color=c, linewidth=1.2, label=lbl, zorder=3)
        current_offset += step

    ax.set_xlabel(r"Raman Shift $[\mathrm{cm}^{-1}]$")
    ax.set_ylabel(r"Normalized Intensity / Offset $[\mathrm{a.u.}]$")
    if title:
        ax.set_title(title, fontsize=9, pad=6)

    ax.legend(loc="upper right", frameon=False, fontsize=8)
    if created_fig:
        finalize_figure_layout(fig)
    return fig, ax


def plot_raman_deconvolution(
    data: CharacterizationData,
    fit_result: Dict[str, Any],
    peak_labels: Optional[List[str]] = None,
    title: Optional[str] = None,
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    """
    Generate publication 2-panel Raman deconvolution figure with components and residuals.
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

    # Top panel
    ax_top.scatter(x, y, color="#444444", s=10, alpha=0.5, label="Raw Spectrum", zorder=2)
    ax_top.plot(x, fit_y, color="#D55E00", linewidth=1.5, label="Total Fit", zorder=4)

    comps = fit_result.get("components", {})
    comp_colors = get_palette("tol_muted", max(len(comps), 2))
    c_idx = 0
    for name, curve in comps.items():
        if name.startswith("bg_") or name == "linear":
            ax_top.plot(x, curve, color="#999999", linestyle="--", linewidth=1.0, label="Background", zorder=3)
        else:
            c = comp_colors[c_idx % len(comp_colors)]
            c_lbl = peak_labels[c_idx] if peak_labels and c_idx < len(peak_labels) else f"Peak {c_idx+1}"
            c_idx += 1
            ax_top.plot(x, curve, color=c, linestyle="-", linewidth=1.2, alpha=0.85, label=c_lbl, zorder=3)

    r2 = fit_result.get("r_squared", 1.0)
    ax_top.text(
        0.04,
        0.90,
        rf"$R^2 = {r2:.4f}$",
        transform=ax_top.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.9),
    )

    ax_top.set_ylabel(r"Intensity $[\mathrm{counts}]$")
    if title:
        ax_top.set_title(title, fontsize=9, pad=6)
    ax_top.legend(loc="upper right", frameon=False, fontsize=7)

    # Bottom panel: Residuals
    ax_bot.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax_bot.scatter(x, res, color="#0072B2", s=8, alpha=0.7, zorder=2)
    ax_bot.set_xlabel(r"Raman Shift $[\mathrm{cm}^{-1}]$")
    ax_bot.set_ylabel(r"$\Delta I$")

    finalize_figure_layout(fig)
    return fig, (ax_top, ax_bot)

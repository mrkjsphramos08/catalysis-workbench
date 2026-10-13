"""
Fourier Transform Infrared Spectroscopy (FTIR) Analysis and Publication Visualization.

Implements:
- Descending wavenumber convention (4000 to 400 cm^-1)
- Baseline subtraction (AsLS, ModPoly, SNIP)
- Band integration and multi-peak deconvolution
- Reference functional group annotations (with scientific disclaimers)
- Publication plotting (single, stacked, shaded functional regions, peak fit & residuals).
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np

from parsers.characterization.common import CharacterizationData
from parsers.characterization.ftir import convert_transmittance_to_absorbance
from modules.characterization.common_peaks import (
    detect_peaks,
    fit_peaks,
    integrate_band,
    smooth_spectrum,
    subtract_baseline,
)
from scripts.utils import finalize_figure_layout, apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette

# Standard reference functional group absorption windows (heuristic reference guide)
COMMON_IR_BANDS: List[Dict[str, Any]] = [
    {"name": "O-H / N-H stretch", "range": (3200, 3600), "color": "#56B4E9"},
    {"name": "C-H aliphatic stretch", "range": (2850, 3000), "color": "#009E73"},
    {"name": "C=O carbonyl stretch", "range": (1680, 1750), "color": "#D55E00"},
    {"name": "C=C aromatic / alkene", "range": (1580, 1650), "color": "#CC79A7"},
    {"name": "Fingerprint / M-O", "range": (500, 1000), "color": "#E69F00"},
]


def analyze_ftir_bands(
    data: CharacterizationData,
    prominence: Optional[float] = None,
    distance: int = 15,
) -> List[Dict[str, Any]]:
    """
    Detect major vibrational bands and report wavenumber, intensity, and tentative functional group overlap.

    Notes
    -----
    Tentative functional group matches are reference-based heuristic guides and
    must not be treated as automated structural identification without independent verification.
    """
    # Ensure working in absorbance space for peak finding
    if data.metadata.get("signal_type", "").lower() == "transmittance":
        d_abs = convert_transmittance_to_absorbance(data, clamp_nonpositive=True)
    else:
        d_abs = data

    pk_info = detect_peaks(d_abs.x, d_abs.y, prominence=prominence, distance=distance)
    results = []

    for idx, wn in enumerate(pk_info["x_peaks"]):
        intensity = float(pk_info["y_peaks"][idx])

        # Match with reference ranges
        matched_groups = []
        for g in COMMON_IR_BANDS:
            low, high = g["range"]
            if low <= wn <= high:
                matched_groups.append(g["name"])

        results.append({
            "peak_index": idx + 1,
            "wavenumber_cm1": float(wn),
            "absorbance": intensity,
            "tentative_assignments": matched_groups,
            "disclaimer": "Reference-based heuristic match. Requires researcher validation.",
        })

    return results


def plot_ftir_spectrum(
    data: CharacterizationData,
    annotate_peaks: bool = True,
    show_regions: bool = False,
    title: Optional[str] = None,
    color: str = "#0072B2",
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot single FTIR spectrum with conventional descending wavenumber x-axis (4000 to 400 cm^-1).
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        created_fig = True
    else:
        fig = ax.get_figure()

    ax.plot(data.x, data.y, color=color, linewidth=1.2, label=data.sample_id, zorder=3)

    if show_regions:
        for reg in COMMON_IR_BANDS:
            low, high = reg["range"]
            ax.axvspan(low, high, color=reg["color"], alpha=0.15, label=f"Ref: {reg['name']}")

    if annotate_peaks:
        bands = analyze_ftir_bands(data)
        sorted_bands = sorted(bands, key=lambda b: b["absorbance"], reverse=True)[:5]
        for b in sorted_bands:
            wn = b["wavenumber_cm1"]
            # Find corresponding y on plotted spectrum
            near_idx = np.argmin(np.abs(data.x - wn))
            y_val = float(data.y[near_idx])
            ax.annotate(
                f"{wn:.0f}",
                xy=(wn, y_val),
                xytext=(0, 6 if data.metadata.get("signal_type") == "absorbance" else -12),
                textcoords="offset points",
                ha="center",
                fontsize=7,
                arrowprops=dict(arrowstyle="-", color="#555555", lw=0.6),
            )

    ax.set_xlabel(r"Wavenumber $[\mathrm{cm}^{-1}]$")
    sig_name = data.y_name or "Transmittance"
    unit_name = f" [{data.y_unit}]" if data.y_unit else ""
    ax.set_ylabel(f"{sig_name}{unit_name}")

    if title:
        ax.set_title(title, fontsize=9, pad=6)

    # Invert x-axis (conventional descending FTIR axis: 4000 -> 400 cm^-1)
    x_min = float(np.min(data.x))
    x_max = float(np.max(data.x))
    ax.set_xlim(x_max, x_min)

    if created_fig:
        finalize_figure_layout(fig)
    return fig, ax


def plot_ftir_stacked(
    datasets: Sequence[CharacterizationData],
    labels: Optional[Sequence[str]] = None,
    normalize: bool = True,
    vertical_offset_fraction: float = 0.4,
    title: Optional[str] = None,
    palette: str = "okabe_ito",
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate stacked FTIR spectra with vertical offset on descending wavenumber axis.
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

    ax.set_xlabel(r"Wavenumber $[\mathrm{cm}^{-1}]$")
    ax.set_ylabel(r"Normalized Intensity / Offset $[\mathrm{a.u.}]$")
    if title:
        ax.set_title(title, fontsize=9, pad=6)

    # Invert x-axis
    all_x = np.concatenate([d.x for d in datasets])
    ax.set_xlim(float(np.max(all_x)), float(np.min(all_x)))
    ax.legend(loc="upper right", frameon=False, fontsize=8)

    if created_fig:
        finalize_figure_layout(fig)
    return fig, ax

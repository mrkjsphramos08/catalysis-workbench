"""
DFT Electronic Structure Visualization: Total DOS, Projected DOS (PDOS), and Band Structures.

Provides publication-grade visualizations for total/projected density of states and
electronic band dispersions along high-symmetry k-paths, with energy Fermi-level alignment,
spin polarization channels (spin-up > 0, spin-down < 0), and orbital-projected coloring.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.config import get_figure_dimensions
from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_orbital_color, get_palette


def align_energy(
    energies: Union[np.ndarray, Sequence[float], pd.Series],
    e_fermi: Optional[float] = None,
) -> np.ndarray:
    """
    Align energy values relative to Fermi level (E - E_F).

    Parameters
    ----------
    energies : array-like
        Raw energy values in eV.
    e_fermi : float, optional
        Fermi energy in eV. If None, assumes energies are already referenced to E_F.

    Returns
    -------
    np.ndarray
        Aligned energies relative to E_F.
    """
    arr = np.asarray(energies, dtype=float)
    if e_fermi is not None:
        return arr - float(e_fermi)
    return arr


def plot_dos(
    energies: Union[np.ndarray, Sequence[float], pd.Series],
    dos_up: Union[np.ndarray, Sequence[float], pd.Series],
    dos_down: Optional[Union[np.ndarray, Sequence[float], pd.Series]] = None,
    e_fermi: Optional[float] = None,
    ax: Optional[plt.Axes] = None,
    energy_range: Optional[Tuple[float, float]] = (-6.0, 6.0),
    fill: bool = True,
    color_up: str = "#0072B2",
    color_down: str = "#D55E00",
    label_up: str = "Spin Up",
    label_down: str = "Spin Down",
    orientation: str = "horizontal",
    show_fermi_line: bool = True,
    fermi_line_style: str = "--",
    fermi_line_color: str = "#333333",
    energy_label: str = r"$E - E_{\mathrm{F}}$ (eV)",
    dos_label: str = "DOS (states/eV)",
) -> plt.Axes:
    """
    Plot publication-quality Total Density of States (TDOS).

    Supports spin-unpolarized and spin-polarized channels (spin-up > 0, spin-down < 0).

    Parameters
    ----------
    energies : array-like
        Energy array in eV.
    dos_up : array-like
        Total DOS or spin-up DOS.
    dos_down : array-like, optional
        Spin-down DOS. If provided, plotted as negative DOS values.
    e_fermi : float, optional
        Fermi level in eV.
    ax : plt.Axes, optional
        Existing Axes to plot on. If None, creates a new publication figure.
    energy_range : tuple of (float, float), optional
        Energy window to display.
    fill : bool, default True
        Whether to shade the area under the curve.
    color_up, color_down : str
        Hex color codes for spin channels.
    label_up, label_down : str
        Legend labels.
    orientation : str, default 'horizontal'
        'horizontal' puts energy on x-axis; 'vertical' puts energy on y-axis.
    show_fermi_line : bool, default True
        Whether to draw reference dashed line at E - E_F = 0.
    fermi_line_style : str, default '--'
        Line style for Fermi level.
    fermi_line_color : str, default '#333333'
        Color of Fermi level indicator.
    energy_label, dos_label : str
        Axis labels with physical units.

    Returns
    -------
    plt.Axes
        Matplotlib Axes containing the plot.
    """
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    e_aligned = align_energy(energies, e_fermi)
    d_up = np.asarray(dos_up, dtype=float)

    if orientation == "horizontal":
        ax.plot(e_aligned, d_up, color=color_up, label=label_up if dos_down is not None else "Total DOS")
        if fill:
            # Fill valence band states (E <= 0)
            val_mask = e_aligned <= 0.0
            ax.fill_between(e_aligned[val_mask], 0, d_up[val_mask], color=color_up, alpha=0.35)

        if dos_down is not None:
            d_down = np.asarray(dos_down, dtype=float)
            # Invert spin down
            ax.plot(e_aligned, -d_down, color=color_down, label=label_down)
            if fill:
                ax.fill_between(e_aligned[val_mask], 0, -d_down[val_mask], color=color_down, alpha=0.35)
            ax.axhline(0, color="#888888", linewidth=0.6, linestyle="-")

        if show_fermi_line:
            ax.axvline(0.0, color=fermi_line_color, linestyle=fermi_line_style, linewidth=0.8, label=r"$E_{\mathrm{F}}$")

        ax.set_xlabel(energy_label)
        ax.set_ylabel(dos_label)
        if energy_range:
            ax.set_xlim(energy_range)

    elif orientation == "vertical":
        ax.plot(d_up, e_aligned, color=color_up, label=label_up if dos_down is not None else "Total DOS")
        if fill:
            val_mask = e_aligned <= 0.0
            ax.fill_betweenx(e_aligned[val_mask], 0, d_up[val_mask], color=color_up, alpha=0.35)

        if dos_down is not None:
            d_down = np.asarray(dos_down, dtype=float)
            ax.plot(-d_down, e_aligned, color=color_down, label=label_down)
            if fill:
                ax.fill_betweenx(e_aligned[val_mask], 0, -d_down[val_mask], color=color_down, alpha=0.35)
            ax.axvline(0, color="#888888", linewidth=0.6, linestyle="-")

        if show_fermi_line:
            ax.axhline(0.0, color=fermi_line_color, linestyle=fermi_line_style, linewidth=0.8, label=r"$E_{\mathrm{F}}$")

        ax.set_xlabel(dos_label)
        ax.set_ylabel(energy_label)
        if energy_range:
            ax.set_ylim(energy_range)
    else:
        raise ValueError(f"Unknown orientation '{orientation}'. Choose 'horizontal' or 'vertical'.")

    ax.legend(loc="upper right", framealpha=0.85)
    return ax


def plot_pdos(
    energies: Union[np.ndarray, Sequence[float], pd.Series],
    projections: Dict[str, Union[np.ndarray, Sequence[float], pd.Series, Dict[str, Any]]],
    total_dos: Optional[Union[np.ndarray, Sequence[float], pd.Series]] = None,
    e_fermi: Optional[float] = None,
    ax: Optional[plt.Axes] = None,
    energy_range: Optional[Tuple[float, float]] = (-6.0, 6.0),
    colors: Optional[Dict[str, str]] = None,
    fill_orbitals: bool = True,
    show_fermi_line: bool = True,
    energy_label: str = r"$E - E_{\mathrm{F}}$ (eV)",
    dos_label: str = "PDOS (states/eV)",
) -> plt.Axes:
    """
    Plot Projected Density of States (PDOS) for multiple elements and orbitals.

    Parameters
    ----------
    energies : array-like
        Energy array in eV.
    projections : dict
        Mapping of label -> DOS data or dict with 'up' and 'down' channels.
        Example: {'Pt-d': array, 'O-p': array} or {'Fe-d': {'up': array1, 'down': array2}}
    total_dos : array-like, optional
        Background Total DOS curve for reference.
    e_fermi : float, optional
        Fermi energy in eV.
    ax : plt.Axes, optional
        Existing Axes to plot on.
    energy_range : tuple of (float, float), optional
        Energy window to display.
    colors : dict, optional
        Custom color mapping for projections.
    fill_orbitals : bool, default True
        Whether to shade occupied valence states.
    show_fermi_line : bool, default True
        Whether to draw Fermi level dashed line at 0 eV.
    energy_label, dos_label : str
        Axis labels with units.

    Returns
    -------
    plt.Axes
    """
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    e_aligned = align_energy(energies, e_fermi)
    val_mask = e_aligned <= 0.0

    # Optional background TDOS
    if total_dos is not None:
        t_dos = np.asarray(total_dos, dtype=float)
        ax.plot(e_aligned, t_dos, color="#999999", linewidth=0.8, linestyle=":", label="Total DOS")
        ax.fill_between(e_aligned, 0, t_dos, color="#E0E0E0", alpha=0.3)

    proj_keys = list(projections.keys())
    has_spin_down = False

    for idx, key in enumerate(proj_keys):
        data = projections[key]
        if colors and key in colors:
            col = colors[key]
        else:
            # Check orbital heuristic
            last_part = key.split("-")[-1].lower() if "-" in key else key.lower()
            col = get_orbital_color(last_part) if last_part in ("s", "p", "d", "f", "t2g", "eg") else get_color(idx)

        if isinstance(data, dict) and "up" in data:
            up_series = np.asarray(data["up"], dtype=float)
            ax.plot(e_aligned, up_series, color=col, linewidth=1.1, label=f"{key}")
            if fill_orbitals:
                ax.fill_between(e_aligned[val_mask], 0, up_series[val_mask], color=col, alpha=0.25)

            if "down" in data and data["down"] is not None:
                has_spin_down = True
                down_series = np.asarray(data["down"], dtype=float)
                ax.plot(e_aligned, -down_series, color=col, linewidth=1.1)
                if fill_orbitals:
                    ax.fill_between(e_aligned[val_mask], 0, -down_series[val_mask], color=col, alpha=0.25)
        else:
            series = np.asarray(data, dtype=float)
            ax.plot(e_aligned, series, color=col, linewidth=1.1, label=key)
            if fill_orbitals:
                ax.fill_between(e_aligned[val_mask], 0, series[val_mask], color=col, alpha=0.25)

    if has_spin_down:
        ax.axhline(0, color="#888888", linewidth=0.5, linestyle="-")

    if show_fermi_line:
        ax.axvline(0.0, color="#333333", linestyle="--", linewidth=0.8, label=r"$E_{\mathrm{F}}$")

    ax.set_xlabel(energy_label)
    ax.set_ylabel(dos_label)
    if energy_range:
        ax.set_xlim(energy_range)

    ax.legend(loc="upper right", framealpha=0.85)
    return ax


def plot_band_structure(
    k_distances: Union[np.ndarray, Sequence[float], pd.Series],
    bands: Union[np.ndarray, Sequence[Sequence[float]]],
    high_symmetry_points: Optional[List[Tuple[float, str]]] = None,
    e_fermi: Optional[float] = None,
    bands_down: Optional[Union[np.ndarray, Sequence[Sequence[float]]]] = None,
    ax: Optional[plt.Axes] = None,
    energy_range: Optional[Tuple[float, float]] = (-5.0, 5.0),
    color_up: str = "#0072B2",
    color_down: str = "#D55E00",
    show_fermi_line: bool = True,
    energy_label: str = r"$E - E_{\mathrm{F}}$ (eV)",
    k_label: str = "Wavevector",
) -> plt.Axes:
    """
    Plot electronic band structure along reciprocal-space k-path.

    Parameters
    ----------
    k_distances : array-like
        Cumulative reciprocal distance along path (1D array).
    bands : 2D array-like of shape (n_bands, n_kpoints) or (n_kpoints, n_bands)
        Band energies in eV for spin-up (or unpolarized).
    high_symmetry_points : list of (distance, label), optional
        High-symmetry point positions and LaTeX labels (e.g. [(0.0, r'$\\Gamma$'), (1.0, 'X')]).
    e_fermi : float, optional
        Fermi energy in eV.
    bands_down : 2D array-like, optional
        Band energies for spin-down channel.
    ax : plt.Axes, optional
        Matplotlib Axes to plot on.
    energy_range : tuple of (float, float), optional
        Energy window in eV.
    color_up, color_down : str
        Hex color codes.
    show_fermi_line : bool, default True
        Whether to draw dashed line at E = 0.
    energy_label, k_label : str
        Axis labels.

    Returns
    -------
    plt.Axes
    """
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    k_arr = np.asarray(k_distances, dtype=float)
    b_up = np.asarray(bands, dtype=float)

    # Ensure shape is (n_bands, n_kpoints)
    if b_up.shape[0] == len(k_arr) and b_up.shape[1] != len(k_arr):
        b_up = b_up.T

    if e_fermi is not None:
        b_up = b_up - float(e_fermi)

    for i in range(b_up.shape[0]):
        ax.plot(
            k_arr,
            b_up[i, :],
            color=color_up,
            linewidth=0.9,
            label="Spin Up" if (i == 0 and bands_down is not None) else None,
        )

    if bands_down is not None:
        b_down = np.asarray(bands_down, dtype=float)
        if b_down.shape[0] == len(k_arr) and b_down.shape[1] != len(k_arr):
            b_down = b_down.T
        if e_fermi is not None:
            b_down = b_down - float(e_fermi)

        for i in range(b_down.shape[0]):
            ax.plot(
                k_arr,
                b_down[i, :],
                color=color_down,
                linewidth=0.9,
                linestyle="--",
                label="Spin Down" if i == 0 else None,
            )

    if show_fermi_line:
        ax.axhline(0.0, color="#333333", linestyle="--", linewidth=0.8, label=r"$E_{\mathrm{F}}$")

    # High symmetry ticks and vertical partition lines
    if high_symmetry_points:
        tick_pos = [pt[0] for pt in high_symmetry_points]
        tick_labels = [pt[1] for pt in high_symmetry_points]
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(tick_labels)
        for pos in tick_pos:
            ax.axvline(pos, color="#B0B0B0", linewidth=0.6, linestyle="-")

    ax.set_xlim(float(k_arr.min()), float(k_arr.max()))
    if energy_range:
        ax.set_ylim(energy_range)

    ax.set_ylabel(energy_label)
    ax.set_xlabel(k_label)
    if bands_down is not None:
        ax.legend(loc="upper right", framealpha=0.85)

    return ax

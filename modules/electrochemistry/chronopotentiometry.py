"""
Chronopotentiometry (CP) Galvanostatic Stability Analysis and Publication Plotting.

Processes chronopotentiometric potential-time transient curves:
- Galvanostatic stability benchmarks (e.g. holding at 10 mA/cm^2 for 10-100 h)
- Potential drift rate extraction (mV/h)
- RHE and iR compensation integration
- Overpotential evolution under prolonged Faradaic operation.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np

from modules.electrochemistry.reference_electrodes import convert_potential_to_rhe
from parsers.electrochemistry.schema import StandardizedElectrochemicalData
from scripts.config import get_figure_dimensions
from scripts.utils import finalize_figure_layout, apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


def calculate_cp_stability_metrics(
    time_s: np.ndarray,
    potential_V: np.ndarray,
    initial_settling_time_s: float = 10.0,
    average_window_s: float = 30.0,
) -> Dict[str, Any]:
    """
    Compute galvanostatic stability metrics, potential drift, and drift rate.

    Parameters
    ----------
    time_s : np.ndarray
        Elapsed time in seconds.
    potential_V : np.ndarray
        Potential array in Volts.
    initial_settling_time_s : float, default 10.0
        Initial settling period bypassed before initial potential evaluation.
    average_window_s : float, default 30.0
        Time window over which start and end potentials are averaged.

    Returns
    -------
    Dict[str, Any]
        Dictionary containing initial potential, final potential, total drift (mV),
        drift rate (mV/h), and duration.
    """
    t = np.asarray(time_s, dtype=float)
    pot = np.asarray(potential_V, dtype=float)

    if len(t) < 5 or t[-1] <= initial_settling_time_s:
        raise ValueError(
            f"Insufficient time data for CP metrics. Duration {t[-1]:.1f} s <= settling time {initial_settling_time_s:.1f} s."
        )

    # Initial potential
    init_mask = (t >= initial_settling_time_s) & (t <= initial_settling_time_s + average_window_s)
    if not np.any(init_mask):
        init_mask = t >= initial_settling_time_s
        e_initial = float(pot[init_mask][0]) if np.any(init_mask) else float(pot[0])
    else:
        e_initial = float(np.mean(pot[init_mask]))

    # Final potential
    final_mask = t >= (t[-1] - average_window_s)
    if not np.any(final_mask):
        e_final = float(pot[-1])
    else:
        e_final = float(np.mean(pot[final_mask]))

    # Duration
    duration_s = float(t[-1] - t[0])
    duration_h = duration_s / 3600.0

    # Total drift in mV
    drift_mV = (e_final - e_initial) * 1000.0

    # Drift rate in mV/h
    drift_rate_mV_per_h = (drift_mV / duration_h) if duration_h > 0.01 else 0.0

    return {
        "e_initial_V": e_initial,
        "e_final_V": e_final,
        "potential_drift_mV": drift_mV,
        "drift_rate_mV_per_h": drift_rate_mV_per_h,
        "duration_seconds": duration_s,
        "duration_hours": duration_h,
    }


def plot_chronopotentiometry(
    data: Union[StandardizedElectrochemicalData, Sequence[StandardizedElectrochemicalData]],
    labels: Optional[Union[str, List[str]]] = None,
    time_unit: str = "h",
    convert_to_rhe: bool = True,
    show_drift_rate: bool = True,
    current_label: Optional[str] = r"$j = 10~\mathrm{mA\cdot cm^{-2}}$",
    ph: Optional[float] = None,
    title: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
    ax: Optional[plt.Axes] = None,
    **kwargs: Any,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate publication-ready Chronopotentiometry (potential vs time) stability curves.

    Parameters
    ----------
    data : StandardizedElectrochemicalData or Sequence
        Electrochemical data object(s).
    labels : str or List[str], optional
        Dataset labels.
    time_unit : str, default 'h'
        'h', 'min', or 's'.
    convert_to_rhe : bool, default True
        If True, converts potential to V vs RHE.
    show_drift_rate : bool, default True
        If True, annotates drift rate (mV/h) on single curve.
    current_label : str, optional
        Operating current density descriptor.
    title : str, optional
        Figure title.
    save_path : str or Path, optional
        File destination path.
    ax : plt.Axes, optional
        Existing axes to plot onto.

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
    """
    apply_publication_style()

    if isinstance(data, StandardizedElectrochemicalData):
        datasets = [data]
    else:
        datasets = list(data)

    if labels is None:
        plot_labels = [f"Catalyst {i+1}" if len(datasets) > 1 else "Stability" for i in range(len(datasets))]
    elif isinstance(labels, str):
        plot_labels = [labels]
    else:
        plot_labels = list(labels)

    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio=1.33)
        created_fig = True
    else:
        fig = ax.get_figure()

    colors = [get_color(i) for i in range(len(datasets))]

    time_scales = {"s": 1.0, "min": 60.0, "h": 3600.0}
    t_div = time_scales.get(time_unit.lower(), 3600.0)
    t_label = time_unit.lower()

    for i, d in enumerate(datasets):
        if d.time_s is None:
            continue

        c = colors[i % len(colors)]
        lbl = plot_labels[i] if i < len(plot_labels) else f"Series {i+1}"

        t_plot = (d.time_s - d.time_s[0]) / t_div

        if convert_to_rhe:
            sol_ph = ph if ph is not None else d.ph
            pot_plot = d.get_potential_vs_rhe(ph=sol_ph)
            y_label = r"$E~[\mathrm{V~vs.~RHE}]$"
        else:
            pot_plot = d.potential_V
            ref = d.potential_reference_raw or "Ref"
            y_label = rf"$E~[\mathrm{{V~vs.~{ref}}}]$"

        ax.plot(t_plot, pot_plot, label=lbl, color=c, linewidth=1.4, zorder=3)

        if show_drift_rate and len(datasets) == 1:
            try:
                metrics = calculate_cp_stability_metrics(d.time_s, pot_plot)
                drift_rate = metrics["drift_rate_mV_per_h"]
                drift_text = f"Drift: {drift_rate:+.2f} mV/h"
                ax.annotate(
                    drift_text,
                    xy=(t_plot[-1], pot_plot[-1]),
                    xytext=(-15, 12),
                    textcoords="offset points",
                    ha="right",
                    fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#cccccc", alpha=0.85),
                )
            except Exception:
                pass

    ax.set_xlabel(rf"$\mathrm{{Time}}~[t / \mathrm{{{t_label}}}]$")
    ax.set_ylabel(y_label)
    ax.grid(True, linestyle="--", alpha=0.3, zorder=0)

    if current_label:
        ax.text(
            0.05,
            0.92,
            current_label,
            transform=ax.transAxes,
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#cccccc", alpha=0.85),
        )

    if len(datasets) > 1:
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

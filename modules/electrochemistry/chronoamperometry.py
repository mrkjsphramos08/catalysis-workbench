"""
Chronoamperometry (CA) Stability Analysis, Current Retention, and Publication Plotting.

Processes potentiostatic current-time transient curves:
- Long-term electrocatalytic stability metrics (retention %, degradation rate per hour)
- Initial capacitive current spike isolation
- Multi-catalyst durability comparisons
- Standard IUPAC current conventions (+ anodic, - cathodic).
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np

from parsers.electrochemistry.schema import StandardizedElectrochemicalData
from scripts.config import get_figure_dimensions
from scripts.utils import finalize_figure_layout, apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


def calculate_ca_stability_metrics(
    time_s: np.ndarray,
    current_density: np.ndarray,
    initial_settling_time_s: float = 10.0,
    average_window_s: float = 30.0,
) -> Dict[str, Any]:
    """
    Compute long-term potentiostatic stability metrics from chronoamperometry data.

    Parameters
    ----------
    time_s : np.ndarray
        Elapsed time in seconds.
    current_density : np.ndarray
        Current density array (mA/cm^2 or A).
    initial_settling_time_s : float, default 10.0
        Settling time to bypass initial non-Faradaic double-layer charging spike.
    average_window_s : float, default 30.0
        Duration over which initial and final current levels are averaged.

    Returns
    -------
    Dict[str, Any]
        Dictionary with initial current, final current, retention percentage,
        degradation rate (%/h), and total test duration.
    """
    t = np.asarray(time_s, dtype=float)
    j = np.asarray(current_density, dtype=float)

    if len(t) < 5 or t[-1] <= initial_settling_time_s:
        raise ValueError(
            f"Insufficient time data for CA stability metrics. Duration {t[-1]:.1f} s <= settling time {initial_settling_time_s:.1f} s."
        )

    # Initial current window: [initial_settling_time_s, initial_settling_time_s + average_window_s]
    init_mask = (t >= initial_settling_time_s) & (t <= initial_settling_time_s + average_window_s)
    if not np.any(init_mask):
        # Fallback to points right after settling time
        init_mask = t >= initial_settling_time_s
        j_initial = float(j[init_mask][0]) if np.any(init_mask) else float(j[0])
    else:
        j_initial = float(np.mean(j[init_mask]))

    # Final current window: [t_end - average_window_s, t_end]
    final_mask = t >= (t[-1] - average_window_s)
    if not np.any(final_mask):
        j_final = float(j[-1])
    else:
        j_final = float(np.mean(j[final_mask]))

    # Duration
    duration_s = float(t[-1] - t[0])
    duration_h = duration_s / 3600.0

    # Retention percentage: (j_final / j_initial) * 100
    # Handles both anodic (+j) and cathodic (-j) magnitudes
    if abs(j_initial) > 1e-12:
        retention_pct = float((abs(j_final) / abs(j_initial)) * 100.0)
    else:
        retention_pct = 0.0

    # Linear degradation rate (% loss per hour)
    if duration_h > 0.01:
        degradation_rate_pct_per_h = float((100.0 - retention_pct) / duration_h)
    else:
        degradation_rate_pct_per_h = 0.0

    return {
        "j_initial": j_initial,
        "j_final": j_final,
        "retention_percent": retention_pct,
        "degradation_rate_pct_per_h": degradation_rate_pct_per_h,
        "duration_seconds": duration_s,
        "duration_hours": duration_h,
    }


def plot_chronoamperometry(
    data: Union[StandardizedElectrochemicalData, Sequence[StandardizedElectrochemicalData]],
    labels: Optional[Union[str, List[str]]] = None,
    time_unit: str = "h",
    normalize_area: bool = True,
    show_retention: bool = True,
    potential_label: Optional[str] = None,
    title: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
    ax: Optional[plt.Axes] = None,
    **kwargs: Any,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate publication-ready Chronoamperometry (current vs time) stability curves.

    Parameters
    ----------
    data : StandardizedElectrochemicalData or Sequence
        Electrochemical data object(s).
    labels : str or List[str], optional
        Dataset labels.
    time_unit : str, default 'h'
        'h' (hours), 'min' (minutes), or 's' (seconds).
    normalize_area : bool, default True
        If True, plots current density j (mA/cm^2). If False, plots raw current (mA or A).
    show_retention : bool, default True
        Whether to annotate current retention percentage on the figure.
    potential_label : str, optional
        Label for holding potential (e.g., 'E = 1.60 V vs RHE').
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

        if normalize_area:
            j_plot = d.get_current_density(unit="mA/cm^2")
            y_label = r"$j~[\mathrm{mA\cdot cm^{-2}}]$"
        else:
            j_plot = d.current_A * 1e3
            y_label = r"$I~[\mathrm{mA}]$"

        ax.plot(t_plot, j_plot, label=lbl, color=c, linewidth=1.4, zorder=3)

        # Annotate retention if single dataset or requested
        if show_retention and len(datasets) == 1:
            try:
                metrics = calculate_ca_stability_metrics(d.time_s, j_plot)
                ret_pct = metrics["retention_percent"]
                dur_h = metrics["duration_hours"]
                ret_text = f"Retention: {ret_pct:.1f}% ({dur_h:.1f} h)"
                ax.annotate(
                    ret_text,
                    xy=(t_plot[-1], j_plot[-1]),
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

    # Optional holding potential annotation
    if potential_label:
        ax.text(
            0.05,
            0.92,
            potential_label,
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

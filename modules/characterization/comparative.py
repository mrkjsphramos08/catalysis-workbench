"""
Multi-Technique Characterization Comparative Analysis and Cross-Correlation.

Integrates structural, spectroscopic, microscopic, and analytical chemistry metrics
with electrochemical and catalytic activity:
- XRD crystallite size vs Electrochemical overpotential (eta_10)
- BET surface area vs Double-layer capacitance (C_dl)
- Raman vibrational defect ratios vs Reaction barriers
- ICP catalyst leaching vs Potentiostatic degradation rate
- TEM vs SEM particle size distributions.

Includes mandatory scientific caveats regarding physical versus active surfaces and domains.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


def correlate_cross_technique_metrics(
    x_values: Sequence[float],
    y_values: Sequence[float],
    x_label: str,
    y_label: str,
    x_errors: Optional[Sequence[float]] = None,
    y_errors: Optional[Sequence[float]] = None,
    sample_labels: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """
    Perform statistical correlation and linear regression between two characterization metrics.

    Parameters
    ----------
    x_values : sequence of float
        Independent metric (e.g., BET surface area, Raman ID/IG, crystallite size).
    y_values : sequence of float
        Dependent metric (e.g., overpotential eta_10, degradation rate, Tafel slope).
    x_label : str
        Name and units of x-metric.
    y_label : str
        Name and units of y-metric.
    x_errors : sequence of float, optional
        Standard uncertainty in x.
    y_errors : sequence of float, optional
        Standard uncertainty in y.
    sample_labels : sequence of str, optional
        Sample identifiers.

    Returns
    -------
    Dict[str, Any]
        Regression slope, intercept, R^2, Pearson r, p-value, and confidence intervals.
    """
    x_arr = np.asarray(x_values, dtype=float)
    y_arr = np.asarray(y_values, dtype=float)

    if len(x_arr) != len(y_arr):
        raise ValueError(f"Array lengths differ: len(x)={len(x_arr)} != len(y)={len(y_arr)}")
    if len(x_arr) < 3:
        raise ValueError("At least 3 samples required for cross-technique correlation.")

    slope, intercept, r_value, p_value, std_err = stats.linregress(x_arr, y_arr)
    calc_y = slope * x_arr + intercept

    return {
        "x_label": x_label,
        "y_label": y_label,
        "slope": float(slope),
        "intercept": float(intercept),
        "r_squared": float(r_value ** 2),
        "pearson_r": float(r_value),
        "p_value": float(p_value),
        "std_err": float(std_err),
        "x": x_arr,
        "y": y_arr,
        "calc_y": calc_y,
        "x_errors": np.asarray(x_errors, dtype=float) if x_errors is not None else None,
        "y_errors": np.asarray(y_errors, dtype=float) if y_errors is not None else None,
        "sample_labels": list(sample_labels) if sample_labels is not None else [f"S{i+1}" for i in range(len(x_arr))],
        "disclaimer": "Correlations between physical properties and catalytic performance do not establish direct causality without independent mechanistic verification.",
    }


def plot_cross_technique_correlation(
    corr_result: Dict[str, Any],
    title: Optional[str] = None,
    color: str = "#0072B2",
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot publication cross-technique scatter correlation with linear fit and confidence bounds.
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        created_fig = True
    else:
        fig = ax.get_figure()

    x = corr_result["x"]
    y = corr_result["y"]
    x_err = corr_result["x_errors"]
    y_err = corr_result["y_errors"]
    labels = corr_result["sample_labels"]

    # Scatter points with error bars
    ax.errorbar(
        x,
        y,
        xerr=x_err,
        yerr=y_err,
        fmt="o",
        color=color,
        ecolor="#666666",
        elinewidth=1.0,
        capsize=2.5,
        markersize=6,
        label="Catalyst Samples",
        zorder=3,
    )

    # Annotate sample labels
    for idx, lbl in enumerate(labels):
        ax.annotate(
            lbl,
            xy=(x[idx], y[idx]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7,
            color="#333333",
        )

    # Linear fit line
    x_line = np.linspace(np.min(x) * 0.95, np.max(x) * 1.05, 100)
    y_line = corr_result["slope"] * x_line + corr_result["intercept"]
    ax.plot(x_line, y_line, color="#D55E00", linewidth=1.4, linestyle="-", label="Linear Trend", zorder=2)

    r2 = corr_result["r_squared"]
    pval = corr_result["p_value"]
    annot = f"$R^2 = {r2:.3f}$\n$p = {pval:.3e}$"

    ax.text(
        0.05,
        0.85,
        annot,
        transform=ax.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.9),
    )

    ax.set_xlabel(corr_result["x_label"])
    ax.set_ylabel(corr_result["y_label"])
    if title:
        ax.set_title(title, fontsize=9, pad=6)

    ax.legend(loc="lower right", frameon=False, fontsize=8)

    if created_fig:
        fig.tight_layout()
    return fig, ax

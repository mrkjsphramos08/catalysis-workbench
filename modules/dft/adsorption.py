"""
Adsorption Energy Analysis and Catalyst Descriptor Scaling.

Supports rigorous calculation of adsorption and binding energies,
site comparisons, surface comparisons, and linear scaling relations
(e.g., d-band center or OH* vs O* scaling) with configurable reference conventions.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from scripts.config import get_figure_dimensions
from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


def calculate_adsorption_energy(
    e_total: float,
    e_clean_surface: float,
    e_adsorbate_ref: float,
    n_adsorbates: int = 1,
    sign_convention: str = "negative_favorable",
) -> float:
    """
    Calculate adsorption energy per adsorbate molecule.

    Conventional thermodynamics formula:
    E_ads = (E_total - E_clean_surface - sum(E_ref)) / n_adsorbates

    Parameters
    ----------
    e_total : float
        DFT total electronic energy of slab with adsorbate(s) in eV.
    e_clean_surface : float
        DFT total electronic energy of pristine relaxed slab in eV.
    e_adsorbate_ref : float
        Reference energy of the adsorbate in gas phase or stoichiometric reservoir in eV.
    n_adsorbates : int, default 1
        Number of adsorbed molecules on the surface cell (coverage normalization).
    sign_convention : str, default 'negative_favorable'
        'negative_favorable': E_ads < 0 indicates exothermic adsorption.
        'positive_favorable' (Binding Energy convention): E_bind = -E_ads > 0 indicates stable binding.

    Returns
    -------
    float
        Adsorption energy per adsorbate molecule in eV.
    """
    if n_adsorbates <= 0:
        raise ValueError(f"Number of adsorbates must be positive integer, got {n_adsorbates}.")

    e_ads = (e_total - e_clean_surface - e_adsorbate_ref) / float(n_adsorbates)

    if sign_convention == "negative_favorable":
        return e_ads
    elif sign_convention == "positive_favorable":
        return -e_ads
    else:
        raise ValueError(
            f"Unknown sign convention: '{sign_convention}'. "
            "Choose 'negative_favorable' (exothermic < 0) or 'positive_favorable' (binding > 0)."
        )


def plot_adsorption_energies(
    data: pd.DataFrame,
    surface_col: str = "surface",
    adsorption_col: str = "e_ads",
    site_col: Optional[str] = "site",
    error_col: Optional[str] = "error",
    ax: Optional[plt.Axes] = None,
    palette: str = "okabe_ito",
    y_label: str = r"$\Delta E_{\mathrm{ads}}$ (eV)",
    x_label: str = "Catalyst Surface",
    title: Optional[str] = None,
    reference_note: Optional[str] = None,
) -> plt.Axes:
    """
    Plot grouped bar chart comparing adsorption energies across surfaces and sites.

    Parameters
    ----------
    data : pd.DataFrame
        Table containing surface names, energies, and optional sites/errors.
    surface_col : str, default 'surface'
        Column name for catalyst surfaces (x-axis categories).
    adsorption_col : str, default 'e_ads'
        Column name for adsorption energies in eV.
    site_col : str, optional
        Column name for adsorption sites (fcc, hcp, top, bridge).
    error_col : str, optional
        Column name containing standard errors / uncertainty.
    ax : plt.Axes, optional
        Matplotlib Axes to plot on.
    palette : str, default 'okabe_ito'
        Color palette identifier.
    y_label, x_label : str
        Axis labels with units.
    title : str, optional
        Figure title.
    reference_note : str, optional
        Explanatory text annotation for reference zero convention.

    Returns
    -------
    plt.Axes
    """
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    surfaces = data[surface_col].unique()
    n_surfaces = len(surfaces)
    x_indices = np.arange(n_surfaces)

    if site_col and site_col in data.columns and data[site_col].nunique() > 1:
        sites = data[site_col].unique()
        n_sites = len(sites)
        bar_width = 0.8 / n_sites
        colors = get_palette(palette, n_sites)

        for i, site in enumerate(sites):
            sub_df = data[data[site_col] == site].set_index(surface_col).reindex(surfaces)
            y_vals = sub_df[adsorption_col].values
            y_err = sub_df[error_col].values if (error_col and error_col in sub_df.columns) else None

            offsets = x_indices + (i - (n_sites - 1) / 2.0) * bar_width
            ax.bar(
                offsets,
                y_vals,
                width=bar_width * 0.92,
                yerr=y_err,
                capsize=2.5,
                color=colors[i],
                label=str(site),
                edgecolor="#333333",
                linewidth=0.6,
            )
        ax.legend(title="Site", loc="best", framealpha=0.85)
    else:
        # Single site per surface
        y_vals = data.groupby(surface_col)[adsorption_col].mean().reindex(surfaces).values
        y_err = (
            data.groupby(surface_col)[error_col].mean().reindex(surfaces).values
            if (error_col and error_col in data.columns)
            else None
        )
        bar_color = get_color(0, palette)
        ax.bar(
            x_indices,
            y_vals,
            width=0.55,
            yerr=y_err,
            capsize=3.0,
            color=bar_color,
            edgecolor="#333333",
            linewidth=0.6,
        )

    ax.set_xticks(x_indices)
    ax.set_xticklabels(surfaces, rotation=25 if n_surfaces > 4 else 0, ha="right" if n_surfaces > 4 else "center")
    ax.axhline(0, color="#666666", linewidth=0.7, linestyle="--")

    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    if title:
        ax.set_title(title)

    if reference_note:
        ax.annotate(
            reference_note,
            xy=(0.03, 0.05),
            xycoords="axes fraction",
            fontsize=6.8,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="#F8F8F8", edgecolor="#CCCCCC", alpha=0.9),
        )

    return ax


def plot_adsorption_scaling(
    x_descriptor: Union[np.ndarray, Sequence[float], pd.Series],
    y_adsorption: Union[np.ndarray, Sequence[float], pd.Series],
    labels: Optional[Sequence[str]] = None,
    ax: Optional[plt.Axes] = None,
    x_label: str = r"Descriptor $\Delta E_{\mathrm{ref}}$ (eV)",
    y_label: str = r"$\Delta E_{\mathrm{ads}}$ (eV)",
    marker_color: str = "#0072B2",
    line_color: str = "#D55E00",
    show_fit: bool = True,
    show_ci: bool = True,
    interval_type: str = "confidence",
    annotate_points: bool = True,
) -> Tuple[plt.Axes, Dict[str, Any]]:
    """
    Plot linear scaling relation between two adsorption energies or a physical descriptor.

    Fits y = slope * x + intercept with confidence/prediction bands and regression diagnostics.

    Parameters
    ----------
    x_descriptor : array-like
        Independent descriptor variable (e.g., d-band center, Delta E_O*).
    y_adsorption : array-like
        Dependent adsorption energy (e.g., Delta E_OH*).
    labels : sequence of str, optional
        Catalyst surface labels for each data point.
    ax : plt.Axes, optional
        Matplotlib Axes.
    x_label, y_label : str
        Axis labels with units.
    show_fit : bool, default True
        Whether to draw ordinary least squares fit line.
    show_ci : bool, default True
        Whether to show uncertainty band.
    interval_type : str, default 'confidence'
        'confidence' (95% CI for the mean regression line E(y|x)) or
        'prediction' (95% PI for individual future observation y_new).
    annotate_points : bool, default True
        Whether to label individual data points.

    Returns
    -------
    Tuple[plt.Axes, Dict[str, Any]]
        Axes and fit statistics dictionary (slope, intercept, r_squared, adj_r_squared,
        p_value, stderr, intercept_stderr, residual_std_error, f_statistic, interval_type).
    """
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="square")

    x = np.asarray(x_descriptor, dtype=float)
    y = np.asarray(y_adsorption, dtype=float)

    ax.scatter(x, y, color=marker_color, s=28, edgecolor="#222222", linewidth=0.6, zorder=4)

    fit_stats: Dict[str, Any] = {}

    if show_fit and len(x) >= 2:
        res = stats.linregress(x, y)
        x_dense = np.linspace(x.min() - 0.1 * np.ptp(x), x.max() + 0.1 * np.ptp(x), 100)
        y_fit = res.slope * x_dense + res.intercept
        r2 = res.rvalue ** 2

        ax.plot(
            x_dense,
            y_fit,
            color=line_color,
            linewidth=1.2,
            linestyle="-",
            label=f"Fit: $y = {res.slope:.2f}x + {res.intercept:.2f}$\n$R^2 = {r2:.3f}$",
            zorder=3,
        )

        dof = len(x) - 2
        residuals = y - (res.slope * x + res.intercept)
        s_err = float(np.sqrt(np.sum(residuals ** 2) / dof)) if dof > 0 else 0.0
        x_bar = float(np.mean(x))
        ss_x = float(np.sum((x - x_bar) ** 2))

        if show_ci and len(x) > 2 and ss_x > 0:
            t_val = float(stats.t.ppf(0.975, dof))
            int_type_clean = interval_type.strip().lower()

            if int_type_clean == "prediction":
                # 95% Prediction Interval for individual future observation
                band = t_val * s_err * np.sqrt(1.0 + (1.0 / len(x)) + ((x_dense - x_bar) ** 2) / ss_x)
            else:
                # 95% Confidence Interval for mean regression line
                band = t_val * s_err * np.sqrt((1.0 / len(x)) + ((x_dense - x_bar) ** 2) / ss_x)

            ax.fill_between(x_dense, y_fit - band, y_fit + band, color=line_color, alpha=0.15, zorder=2)

        adj_r2 = float(1.0 - (1.0 - r2) * (len(x) - 1) / dof) if dof > 0 else float(r2)
        f_stat = float((r2 / 1.0) / ((1.0 - r2) / dof)) if (1.0 - r2) > 1e-12 and dof > 0 else float("inf")
        intercept_se = float(s_err * np.sqrt((1.0 / len(x)) + (x_bar ** 2) / ss_x)) if ss_x > 0 else 0.0

        fit_stats = {
            "slope": float(res.slope),
            "intercept": float(res.intercept),
            "r_squared": float(r2),
            "adj_r_squared": float(adj_r2),
            "p_value": float(res.pvalue),
            "stderr": float(res.stderr),
            "intercept_stderr": float(intercept_se),
            "residual_std_error": float(s_err),
            "f_statistic": float(f_stat),
            "interval_type": interval_type,
        }
        ax.legend(loc="upper left", framealpha=0.85, fontsize=7.0)

    if annotate_points and labels is not None:
        for xi, yi, lab in zip(x, y, labels):
            ax.annotate(
                str(lab),
                xy=(xi, yi),
                xytext=(3.5, 3.5),
                textcoords="offset points",
                fontsize=6.8,
                color="#222222",
            )

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    return ax, fit_stats

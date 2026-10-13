"""
Chemical engineering reaction kinetics, reactor performance, and diagnostic plotting.
"""

from modules.chemical_engineering.kinetics import (
    fit_arrhenius,
    plot_arrhenius,
    plot_conversion_selectivity,
    plot_reactor_profiles,
    plot_parity_and_residuals,
)

__all__ = [
    "fit_arrhenius",
    "plot_arrhenius",
    "plot_conversion_selectivity",
    "plot_reactor_profiles",
    "plot_parity_and_residuals",
]

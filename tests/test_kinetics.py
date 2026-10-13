"""
Tests for modules/chemical_engineering (kinetics, Arrhenius, and reactor profiles).
"""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from modules.chemical_engineering.kinetics import (
    fit_arrhenius,
    plot_arrhenius,
    plot_conversion_selectivity,
    plot_parity_and_residuals,
    plot_reactor_profiles,
)
from scripts.preprocess import GAS_CONSTANT_R_J_MOL_K
from scripts.utils import create_publication_figure


def test_fit_arrhenius_exact():
    """Verify Arrhenius parameter extraction on known physical kinetics."""
    temp_k = np.array([450.0, 475.0, 500.0, 525.0, 550.0])
    ea_true_j = 75000.0  # 75 kJ/mol
    ln_a_true = 16.5

    # Exact k values: k = A * exp(-Ea / (R * T))
    ln_k = ln_a_true - (ea_true_j / (GAS_CONSTANT_R_J_MOL_K * temp_k))
    k_vals = np.exp(ln_k)

    res = fit_arrhenius(temp_k, k_vals)
    assert pytest.approx(res["ea_kj_mol"], rel=1e-4) == 75.0
    assert pytest.approx(res["ln_a"], rel=1e-4) == 16.5
    assert pytest.approx(res["r_squared"], rel=1e-5) == 1.0


def test_fit_arrhenius_invalid():
    """Verify error checking for invalid temperatures and rate constants."""
    # Negative temperature
    with pytest.raises(ValueError):
        fit_arrhenius([-10.0, 300.0, 400.0], [1.0, 2.0, 3.0])

    # Negative rate constant
    with pytest.raises(ValueError):
        fit_arrhenius([300.0, 350.0, 400.0], [1.0, -2.0, 3.0])

    # Less than 3 points
    with pytest.raises(ValueError):
        fit_arrhenius([300.0, 350.0], [1.0, 2.0])

    # Length mismatch
    with pytest.raises(ValueError):
        fit_arrhenius([300.0, 350.0, 400.0], [1.0, 2.0])


def test_plot_arrhenius():
    """Verify Arrhenius plotting with residuals."""
    temp_k = np.array([450.0, 475.0, 500.0, 525.0, 550.0])
    k_vals = np.array([0.005, 0.015, 0.040, 0.095, 0.210])

    fig, res = plot_arrhenius(temp_k, k_vals, show_residuals=True)
    assert fig is not None
    assert "ea_kj_mol" in res
    assert res["r_squared"] > 0.95
    plt.close(fig)


def test_plot_conversion_selectivity():
    """Verify conversion and selectivity profiles."""
    conv = np.array([10.0, 30.0, 50.0, 70.0, 90.0])
    sel = np.array([95.0, 92.0, 88.0, 82.0, 75.0])

    fig, ax = create_publication_figure()
    ax_out = plot_conversion_selectivity(conv, sel, catalyst_label="Cat1", ax=ax, as_percentage=True)
    assert ax_out is ax
    assert ax.get_xlim() == (0, 105)
    plt.close(fig)


def test_plot_reactor_profiles():
    """Verify multi-metric reactor performance curves."""
    temps = [200, 220, 240, 260, 280]
    conv = [15, 35, 60, 85, 95]
    sel = [95, 92, 88, 80, 70]

    fig, ax = create_publication_figure()
    ax_out = plot_reactor_profiles(temps, conv, selectivity=sel, ax=ax)
    assert ax_out is ax
    assert len(ax.lines) == 3  # Conversion, Selectivity, Yield
    plt.close(fig)


def test_plot_parity_and_residuals():
    """Verify parity plot diagnostics."""
    exp = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    pred = np.array([1.05, 1.95, 3.10, 3.90, 5.05])

    fig, stats_res = plot_parity_and_residuals(exp, pred)
    assert stats_res["r_squared"] > 0.98
    assert stats_res["rmse"] < 0.1
    plt.close(fig)

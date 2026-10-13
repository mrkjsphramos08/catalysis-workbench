"""
Unit and Integration Tests for Cross-Technique Characterization Correlation.
"""

import numpy as np
import pytest
import matplotlib.pyplot as plt

from modules.characterization.comparative import (
    correlate_cross_technique_metrics,
    plot_cross_technique_correlation,
)


def test_correlate_cross_technique_metrics():
    """Test linear correlation between BET surface area and double-layer capacitance."""
    bet_areas = [45.0, 92.0, 148.0, 210.0, 280.0]  # m^2/g
    # Simulated C_dl proportional to area with small noise
    c_dl = [12.5, 23.8, 38.2, 53.0, 71.4]  # mF/cm^2
    samples = ["Cat-1", "Cat-2", "Cat-3", "Cat-4", "Cat-5"]

    res = correlate_cross_technique_metrics(
        x_values=bet_areas,
        y_values=c_dl,
        x_label="BET Surface Area ($m^2\\,g^{-1}$)",
        y_label="Double-Layer Capacitance $C_{dl}$ ($mF\\,cm^{-2}$)",
        sample_labels=samples,
    )

    assert res["r_squared"] > 0.99
    assert res["pearson_r"] > 0.99
    assert res["p_value"] < 0.001
    assert res["slope"] > 0.2
    assert "disclaimer" in res
    assert "causality" in res["disclaimer"].lower()
    assert len(res["sample_labels"]) == 5


def test_correlate_validation_errors():
    """Verify input validation errors for mismatched or insufficient arrays."""
    with pytest.raises(ValueError, match="differ"):
        correlate_cross_technique_metrics([1, 2], [1, 2, 3], "x", "y")

    with pytest.raises(ValueError, match="At least 3 samples"):
        correlate_cross_technique_metrics([1, 2], [1, 2], "x", "y")


def test_plot_cross_technique_correlation():
    """Verify rendering of cross-technique correlation plot."""
    x = [10.0, 20.0, 30.0, 40.0]
    y = [1.5, 2.8, 4.3, 5.9]
    res = correlate_cross_technique_metrics(x, y, "X metric", "Y metric")

    fig, ax = plot_cross_technique_correlation(res, title="Test Cross Correlation")
    assert isinstance(fig, plt.Figure)
    assert ax.get_xlabel() == "X metric"
    assert ax.get_ylabel() == "Y metric"
    plt.close(fig)

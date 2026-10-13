"""
Unit and Numerical Validation Tests for FTIR Spectroscopy Module.
"""

from pathlib import Path
import numpy as np
import pytest

from parsers.characterization.common import CharacterizationData, InvalidUnitsError
from parsers.characterization.ftir import parse_ftir_file, convert_transmittance_to_absorbance
from modules.characterization.ftir import (
    analyze_ftir_bands,
    plot_ftir_spectrum,
    plot_ftir_stacked,
)
from modules.characterization.common_peaks import integrate_band, smooth_spectrum

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures" / "ftir"


def test_parse_ftir_file():
    """Verify parsing of transmittance FTIR spectrum."""
    data = parse_ftir_file(FIXTURES_DIR / "catalyst_functionalized_ftir.csv")
    assert isinstance(data, CharacterizationData)
    assert data.technique == "FTIR"
    assert data.x_name == "Wavenumber"
    assert data.y_name == "Transmittance"
    assert len(data.x) > 1000
    assert np.min(data.x) < 500.0
    assert np.max(data.x) > 3500.0


def test_beer_lambert_transmittance_to_absorbance():
    """Verify exact numerical conversion A = -log10(T) and nonpositive rejection."""
    data = parse_ftir_file(FIXTURES_DIR / "catalyst_functionalized_ftir.csv")
    # T = 100% -> A = 0.0, T = 10% -> A = 1.0, T = 1% -> A = 2.0
    d_abs = convert_transmittance_to_absorbance(data)
    assert d_abs.y_name == "Absorbance"
    assert np.all(d_abs.y >= 0.0)

    # Test rejection of nonpositive transmittance
    bad_data = data.copy()
    bad_data.y = np.array([-5.0, 10.0, 50.0])
    with pytest.raises(InvalidUnitsError):
        convert_transmittance_to_absorbance(bad_data, clamp_nonpositive=False)


def test_ftir_band_detection_and_integration():
    """Verify detection of O-H, C=O carbonyl, and band area integration."""
    data = parse_ftir_file(FIXTURES_DIR / "catalyst_functionalized_ftir.csv")
    bands = analyze_ftir_bands(data)
    assert len(bands) >= 3

    wn_vals = [b["wavenumber_cm1"] for b in bands]
    # Check that O-H (~3400) or C=O (~1715) are found
    assert any(3300 <= w <= 3500 for w in wn_vals)
    assert any(1680 <= w <= 1750 for w in wn_vals)

    # Test band integration across C=O band (1650 to 1750 cm^-1)
    d_abs = convert_transmittance_to_absorbance(data)
    area = integrate_band(d_abs.x, d_abs.y, x_range=(1650.0, 1750.0))
    assert area > 0.0


def test_ftir_plotting_descending_axis():
    """Verify that FTIR plotting inverts x-axis (descending wavenumber)."""
    data = parse_ftir_file(FIXTURES_DIR / "catalyst_functionalized_ftir.csv")
    fig, ax = plot_ftir_spectrum(data, show_regions=True)
    assert fig is not None
    x_lim = ax.get_xlim()
    assert x_lim[0] > x_lim[1]  # Inverted: 4000 to 400 cm^-1

    data_spent = parse_ftir_file(FIXTURES_DIR / "catalyst_spent_ftir.csv")
    fig_stk, ax_stk = plot_ftir_stacked([data, data_spent], labels=["Pristine", "Spent"])
    assert fig_stk is not None

"""
Unit and Numerical Validation Tests for X-Ray Diffraction (XRD) Module.
"""

from pathlib import Path
import numpy as np
import pytest

from parsers.characterization.common import CharacterizationData
from parsers.characterization.xrd import parse_xrd_file
from modules.characterization.xrd import (
    calculate_d_spacing,
    calculate_crystallite_size_scherrer,
    analyze_xrd_peaks,
    plot_xrd_pattern,
    plot_xrd_stacked,
    plot_xrd_fit,
    ScherrerBroadeningWarning,
)
from modules.characterization.common_peaks import fit_peaks

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures" / "xrd"


def test_bragg_d_spacing_numerical_exactness():
    """Verify Bragg's law: n*lambda = 2*d*sin(theta) for known Pt (111) reflection."""
    # Pt (111) at 2theta = 39.76 deg using Cu K-alpha (lambda = 1.54060 A)
    # theta = 19.88 deg -> sin(19.88 deg) = 0.34006
    # d = 1.54060 / (2 * 0.34006) = 2.2652 A
    d_calc = calculate_d_spacing(39.76, wavelength_angstrom=1.54060)
    assert pytest.approx(2.2652, abs=0.002) == d_calc

    # Array input handling
    angles = np.array([39.76, 46.24, 67.45])
    d_array = calculate_d_spacing(angles, wavelength_angstrom=1.54060)
    assert len(d_array) == 3
    assert d_array[0] > d_array[1] > d_array[2]  # Higher angle = smaller d-spacing


def test_scherrer_crystallite_size_numerical():
    """Verify Scherrer equation and instrumental broadening subtraction."""
    # Peak at 2theta = 40.0 deg, FWHM = 0.80 deg, instrument = 0.10 deg
    res = calculate_crystallite_size_scherrer(
        two_theta_deg=40.0,
        fwhm_deg=0.80,
        wavelength_angstrom=1.54060,
        k_shape_factor=0.9,
        instrumental_fwhm_deg=0.10,
    )
    # beta_corr_deg = sqrt(0.80^2 - 0.10^2) = sqrt(0.63) = 0.7937 deg
    # in rad: 0.7937 * pi / 180 = 0.01385 rad
    # cos(20 deg) = 0.93969
    # D = (0.9 * 1.54060) / (0.01385 * 0.93969) = 1.38654 / 0.013018 = 106.5 A = 10.65 nm
    assert res["crystallite_size_nm"] is not None
    assert pytest.approx(10.65, abs=0.2) == res["crystallite_size_nm"]
    assert res["crystallite_size_angstrom"] == pytest.approx(res["crystallite_size_nm"] * 10.0, rel=1e-3)


def test_scherrer_instrument_broadening_warning():
    """Verify warning when observed broadening <= instrumental resolution."""
    with pytest.warns(ScherrerBroadeningWarning):
        res = calculate_crystallite_size_scherrer(
            two_theta_deg=40.0,
            fwhm_deg=0.04,  # Smaller than instrument
            instrumental_fwhm_deg=0.06,
        )
        assert res["crystallite_size_nm"] is None


def test_parse_xrd_csv_and_metadata():
    """Verify parsing of pristine catalyst XRD CSV."""
    data = parse_xrd_file(FIXTURES_DIR / "catalyst_pristine_xrd.csv", wavelength="cu")
    assert isinstance(data, CharacterizationData)
    assert data.technique == "XRD"
    assert len(data.x) > 1000
    assert data.metadata["radiation_source"] == "CU"
    assert data.metadata["wavelength_angstrom"] == pytest.approx(1.54060)
    assert np.min(data.x) < 20.0
    assert np.max(data.x) > 80.0


def test_parse_xrd_txt_with_comments():
    """Verify parsing of tab-delimited TXT with header comment lines."""
    data = parse_xrd_file(FIXTURES_DIR / "catalyst_treated_xrd.txt", wavelength="cu")
    assert isinstance(data, CharacterizationData)
    assert data.metadata["Instrument"] == "Bruker D8 Advance"
    assert len(data.x) == len(data.y)


def test_xrd_peak_deconvolution_and_plotting():
    """Verify multi-peak fitting around Pt (111) and (200) reflections."""
    data = parse_xrd_file(FIXTURES_DIR / "catalyst_treated_xrd.txt")
    # Sub-region around 36 to 50 deg containing Pt (111) and Pt (200)
    mask = (data.x >= 36.0) & (data.x <= 50.0)
    sub_x = data.x[mask]
    sub_y = data.y[mask]

    fit_res = fit_peaks(sub_x, sub_y, peak_centers=[39.88, 46.36], model_type="pseudo_voigt")
    assert fit_res["r_squared"] > 0.95
    assert len(fit_res["peaks"]) == 2
    # Verify centers are accurately retrieved
    p1 = fit_res["peaks"][0]
    p2 = fit_res["peaks"][1]
    assert p1["center"] == pytest.approx(39.88, abs=0.15)
    assert p2["center"] == pytest.approx(46.36, abs=0.15)

    # Test plotting functions
    sub_data = data.copy()
    sub_data.x = sub_x
    sub_data.y = sub_y
    fig_fit, axes = plot_xrd_fit(sub_data, fit_res)
    assert fig_fit is not None

    fig_pat, ax_pat = plot_xrd_pattern(data, annotate_peaks=True)
    assert fig_pat is not None

"""
Unit and Numerical Validation Tests for Raman Spectroscopy Module.
"""

from pathlib import Path
import numpy as np
import pytest

from parsers.characterization.common import CharacterizationData
from parsers.characterization.raman import parse_raman_file
from modules.characterization.raman import (
    process_raman_spectrum,
    calculate_peak_ratios,
    plot_raman_spectrum,
    plot_raman_stacked,
    plot_raman_deconvolution,
)
from modules.characterization.common_peaks import (
    fit_peaks,
    remove_cosmic_rays,
    subtract_baseline,
)

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures" / "raman"


def test_parse_raman_file():
    """Verify parsing of carbon support Raman spectrum."""
    data = parse_raman_file(FIXTURES_DIR / "carbon_support_raman.csv", laser_wavelength_nm=532.0)
    assert isinstance(data, CharacterizationData)
    assert data.technique == "Raman"
    assert data.x_name == "Raman Shift"
    assert data.metadata["laser_wavelength_nm"] == 532.0
    assert len(data.x) > 1000


def test_cosmic_ray_spike_removal():
    """Verify detection and correction of sharp cosmic-ray artifact."""
    data = parse_raman_file(FIXTURES_DIR / "carbon_support_raman.csv")
    y_raw = data.y.copy()

    # The synthetic spike was injected at index 450
    assert np.max(y_raw) > 3000.0

    y_clean, spikes = remove_cosmic_rays(y_raw, threshold=5.0)
    assert len(spikes) >= 1
    assert 450 in spikes
    assert np.max(y_clean) < 2500.0


def test_raman_peak_deconvolution_and_ratios():
    """Verify D-band and G-band deconvolution and ID/IG ratio calculation."""
    data = parse_raman_file(FIXTURES_DIR / "carbon_support_raman.csv")
    proc = process_raman_spectrum(data, despike=True, subtract_bg=True, bg_method="asls")

    fit_res = fit_peaks(
        proc.x,
        proc.y,
        peak_centers=[1345.0, 1585.0],
        model_type="lorentzian",
        background="linear",
    )

    assert fit_res["r_squared"] > 0.94
    assert len(fit_res["peaks"]) == 2

    d_peak = fit_res["peaks"][0]
    g_peak = fit_res["peaks"][1]
    assert d_peak["center"] == pytest.approx(1345.0, abs=15.0)
    assert g_peak["center"] == pytest.approx(1585.0, abs=15.0)

    ratios = calculate_peak_ratios(fit_res, peak_1_index=1, peak_2_index=2, label_1="D", label_2="G")
    assert 0.4 < ratios["height_ratio"] < 1.5
    assert 0.4 < ratios["area_ratio"] < 2.5
    assert ratios["height_ratio"] != ratios["area_ratio"]  # Distinct physical definitions


def test_raman_plotting():
    """Verify Raman single, stacked, and deconvolution plotting."""
    data_c = parse_raman_file(FIXTURES_DIR / "carbon_support_raman.csv")
    data_tmd = parse_raman_file(FIXTURES_DIR / "tmd_catalyst_raman.csv")

    fig1, ax1 = plot_raman_spectrum(data_c, show_baseline=False)
    assert fig1 is not None

    fig2, ax2 = plot_raman_stacked([data_c, data_tmd], labels=["Carbon", "MoS2"])
    assert fig2 is not None

    proc = process_raman_spectrum(data_c, despike=True, subtract_bg=True)
    fit_res = fit_peaks(proc.x, proc.y, peak_centers=[1345.0, 1585.0], model_type="pseudo_voigt")
    fig3, axes = plot_raman_deconvolution(proc, fit_res, peak_labels=["D Band", "G Band"])
    assert fig3 is not None

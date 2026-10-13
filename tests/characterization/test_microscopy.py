"""
Unit and Numerical Validation Tests for Microscopy (SEM/TEM/HRTEM) Module.
"""

from pathlib import Path
import numpy as np
import pytest

from parsers.characterization.common import CharacterizationData, InvalidUnitsError
from parsers.characterization.microscopy import parse_microscopy_image
from modules.characterization.microscopy import (
    segment_particles,
    analyze_hrtem_lattice_fringes,
    plot_microscopy_with_scalebar,
    plot_segmentation_overlay,
    plot_particle_size_distribution,
    plot_hrtem_fft_panel,
)

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures" / "microscopy"


def test_parse_microscopy_sem():
    """Verify parsing and normalization of 2D SEM image."""
    data = parse_microscopy_image(FIXTURES_DIR / "sem_nanoparticles.png", pixel_size_nm=0.50, technique="SEM")
    assert isinstance(data, CharacterizationData)
    assert data.technique == "SEM"
    assert data.image is not None
    assert data.image.shape == (512, 512)
    assert data.pixel_size_nm == 0.50
    assert 0.0 <= np.min(data.image) <= np.max(data.image) <= 1.0


def test_sem_particle_segmentation():
    """Verify particle detection, size distribution, and morphological metrics."""
    data = parse_microscopy_image(FIXTURES_DIR / "sem_nanoparticles.png", pixel_size_nm=0.50)
    res = segment_particles(data, min_area_px=20, dark_particles=True)

    stats = res["statistics"]
    assert stats["particle_count"] > 20
    assert stats["calibrated"] is True
    assert 5.0 < stats["mean_diameter_nm"] < 25.0
    assert 5.0 < stats["median_diameter_nm"] < 25.0
    assert stats["std_diameter_nm"] > 0


def test_hrtem_lattice_fringe_analysis():
    """Verify HRTEM 2D FFT and real-space fringe d-spacing recovery."""
    data = parse_microscopy_image(FIXTURES_DIR / "hrtem_lattice_crop.png", pixel_size_nm=0.025, technique="HRTEM")
    res = analyze_hrtem_lattice_fringes(data, fringe_direction="vertical")

    assert res["d_spacing_nm"] is not None
    assert res["d_spacing_angstrom"] is not None
    # Generated with d = 2.26 A = 0.226 nm
    assert res["d_spacing_angstrom"] == pytest.approx(2.26, abs=0.15)
    assert res["d_spacing_nm"] == pytest.approx(0.226, abs=0.015)
    assert res["fft_power"].shape == (128, 128)


def test_hrtem_missing_calibration_raises_error():
    """Verify that analyzing lattice spacing without physical calibration raises InvalidUnitsError."""
    data_uncal = parse_microscopy_image(FIXTURES_DIR / "hrtem_lattice_crop.png", pixel_size_nm=None, technique="HRTEM")
    with pytest.raises(InvalidUnitsError):
        analyze_hrtem_lattice_fringes(data_uncal)


def test_microscopy_plotting():
    """Verify rendering of scale bar, segmentation overlay, and FFT panel."""
    data_sem = parse_microscopy_image(FIXTURES_DIR / "sem_nanoparticles.png", pixel_size_nm=0.50)
    fig1, ax1 = plot_microscopy_with_scalebar(data_sem, target_scale_nm=50.0)
    assert fig1 is not None

    seg_res = segment_particles(data_sem, min_area_px=20)
    fig2, ax2 = plot_segmentation_overlay(data_sem, seg_res)
    assert fig2 is not None

    fig3, ax3 = plot_particle_size_distribution(seg_res)
    assert fig3 is not None

    data_hrtem = parse_microscopy_image(FIXTURES_DIR / "hrtem_lattice_crop.png", pixel_size_nm=0.025, technique="HRTEM")
    hrtem_res = analyze_hrtem_lattice_fringes(data_hrtem)
    fig4, axes = plot_hrtem_fft_panel(data_hrtem, hrtem_res)
    assert fig4 is not None

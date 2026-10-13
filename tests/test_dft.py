"""
Tests for modules/dft (electronic_structure and adsorption).
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from modules.dft.adsorption import (
    calculate_adsorption_energy,
    plot_adsorption_energies,
    plot_adsorption_scaling,
)
from modules.dft.electronic_structure import (
    align_energy,
    plot_band_structure,
    plot_dos,
    plot_pdos,
)
from scripts.utils import create_publication_figure


def test_align_energy():
    """Verify Fermi-level energy shifting."""
    raw = np.array([-2.0, 0.0, 3.5])
    e_f = -1.5
    aligned = align_energy(raw, e_fermi=e_f)
    assert pytest.approx(aligned[0]) == -0.5
    assert pytest.approx(aligned[1]) == 1.5
    assert pytest.approx(aligned[2]) == 5.0

    # Without e_fermi
    unaligned = align_energy(raw, e_fermi=None)
    assert np.allclose(unaligned, raw)


def test_plot_dos_spin_polarized():
    """Verify spin-polarized DOS plotting."""
    energies = np.linspace(-5.0, 5.0, 50)
    dos_up = np.abs(np.sin(energies)) + 0.1
    dos_down = np.abs(np.cos(energies)) + 0.1

    fig, ax = create_publication_figure()
    ax_out = plot_dos(
        energies=energies,
        dos_up=dos_up,
        dos_down=dos_down,
        e_fermi=0.0,
        ax=ax,
        orientation="horizontal",
    )
    assert ax_out is ax
    assert len(ax.lines) >= 2  # At least spin up and spin down lines
    plt.close(fig)


def test_plot_dos_vertical_orientation():
    """Verify vertical DOS orientation (energy on y-axis)."""
    energies = np.linspace(-4.0, 4.0, 40)
    dos = np.ones_like(energies)

    fig, ax = create_publication_figure()
    plot_dos(energies=energies, dos_up=dos, ax=ax, orientation="vertical")
    assert ax.get_ylabel() == r"$E - E_{\mathrm{F}}$ (eV)"
    plt.close(fig)


def test_plot_pdos():
    """Verify projected DOS with multiple orbitals."""
    energies = np.linspace(-6.0, 6.0, 60)
    projections = {
        "Pt-d": {"up": np.ones(60), "down": np.ones(60) * 0.8},
        "O-p": np.ones(60) * 0.5,
    }
    fig, ax = create_publication_figure()
    ax_out = plot_pdos(energies=energies, projections=projections, total_dos=np.ones(60) * 2.0, ax=ax)
    assert ax_out is ax
    plt.close(fig)


def test_plot_band_structure():
    """Verify band structure along k-path."""
    k_pts = np.linspace(0.0, 2.0, 50)
    bands = np.array([
        -3.0 + np.cos(k_pts * np.pi),
        0.5 - np.cos(k_pts * np.pi),
    ])
    hs_pts = [(0.0, r"$\Gamma$"), (1.0, "$X$"), (2.0, r"$\Gamma$")]

    fig, ax = create_publication_figure()
    ax_out = plot_band_structure(
        k_distances=k_pts,
        bands=bands,
        high_symmetry_points=hs_pts,
        e_fermi=0.0,
        ax=ax,
    )
    assert ax_out is ax
    assert len(ax.get_xticks()) == 3
    plt.close(fig)


def test_calculate_adsorption_energy():
    """Verify adsorption energy calculations with both sign conventions."""
    e_tot = -105.0
    e_clean = -100.0
    e_ref = -3.5  # Gas molecule
    # E_ads = -105 - (-100) - (-3.5) = -1.5 eV (exothermic)

    e_ads_neg = calculate_adsorption_energy(e_tot, e_clean, e_ref, n_adsorbates=1, sign_convention="negative_favorable")
    assert pytest.approx(e_ads_neg) == -1.5

    e_bind_pos = calculate_adsorption_energy(e_tot, e_clean, e_ref, n_adsorbates=1, sign_convention="positive_favorable")
    assert pytest.approx(e_bind_pos) == 1.5

    # Coverage normalization (2 adsorbates)
    e_ads_cov = calculate_adsorption_energy(e_tot, e_clean, 2 * e_ref, n_adsorbates=2, sign_convention="negative_favorable")
    assert pytest.approx(e_ads_cov) == (-105.0 - (-100.0) - (-7.0)) / 2.0

    with pytest.raises(ValueError):
        calculate_adsorption_energy(e_tot, e_clean, e_ref, n_adsorbates=0)
    with pytest.raises(ValueError):
        calculate_adsorption_energy(e_tot, e_clean, e_ref, sign_convention="invalid_convention")


def test_plot_adsorption_scaling():
    """Verify linear scaling fit between descriptors."""
    x_desc = np.array([0.5, 1.0, 1.5, 2.0, 2.5])
    # Exact linear relation: y = 0.75 * x - 0.20
    y_ads = 0.75 * x_desc - 0.20
    labels = ["M1", "M2", "M3", "M4", "M5"]

    fig, ax = create_publication_figure()
    ax_out, stats_res = plot_adsorption_scaling(
        x_descriptor=x_desc,
        y_adsorption=y_ads,
        labels=labels,
        ax=ax,
        show_fit=True,
    )
    assert pytest.approx(stats_res["slope"], rel=1e-4) == 0.75
    assert pytest.approx(stats_res["intercept"], rel=1e-4) == -0.20
    assert pytest.approx(stats_res["r_squared"], rel=1e-4) == 1.0
    assert "adj_r_squared" in stats_res
    assert "residual_std_error" in stats_res
    assert "intercept_stderr" in stats_res
    assert stats_res["interval_type"] == "confidence"
    plt.close(fig)

    # Test prediction interval
    fig_pi, ax_pi = create_publication_figure()
    _, stats_pi = plot_adsorption_scaling(
        x_descriptor=x_desc,
        y_adsorption=y_ads + np.array([0.02, -0.01, 0.03, -0.02, 0.01]),
        interval_type="prediction",
        ax=ax_pi,
    )
    assert stats_pi["interval_type"] == "prediction"
    assert stats_pi["residual_std_error"] > 0.0
    assert stats_pi["adj_r_squared"] > 0.95
    plt.close(fig_pi)


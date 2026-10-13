"""
Tests for modules/catalysis (reaction_energy and volcano).
"""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from modules.catalysis.reaction_energy import (
    apply_che_correction,
    calculate_gibbs_free_energy,
    convert_reference_potential,
    plot_reaction_pathway,
)
from modules.catalysis.volcano import (
    evaluate_her_volcano,
    evaluate_orr_volcano,
    evaluate_volcano_model,
    plot_her_volcano,
    plot_orr_volcano,
    plot_volcano,
)
from scripts.preprocess import KB_EV_PER_K
from scripts.utils import create_publication_figure


def test_calculate_gibbs_free_energy():
    """Verify Gibbs free energy formulation: G = E_DFT + ZPE - TS + corrections."""
    e_dft = -25.50
    zpe = 0.45
    ts = 0.30
    solv = -0.15

    g = calculate_gibbs_free_energy(e_dft, zpe=zpe, ts=ts, corrections=solv)
    expected = -25.50 + 0.45 - 0.30 - 0.15
    assert pytest.approx(g) == expected


def test_apply_che_correction():
    """Verify Computational Hydrogen Electrode (CHE) potential and pH shifts with signed stoichiometry."""
    delta_g0 = 1.23  # eV at U = 0 V vs RHE, pH = 0

    # 1. Oxidation at U = 1.23 V vs RHE, nu_H = -1: Delta G should be 0.0 eV
    g_u_ox = apply_che_correction(delta_g0, u_potential=1.23, ph=0.0, step_type="oxidation")
    assert pytest.approx(g_u_ox, abs=1e-5) == 0.0

    # 2. Reduction at U = 0.50 V vs RHE, nu_H = +1, Delta G0 = 1.0 eV -> Delta G = 1.50 eV
    g_u_red = apply_che_correction(1.0, u_potential=0.50, nu_H=+1.0)
    assert pytest.approx(g_u_red, abs=1e-5) == 1.50

    # 3. Explicit signed stoichiometry nu_H = -1 (oxidation) vs SHE at pH = 14:
    # Delta G = Delta G0 - [U_SHE + kB*T*ln(10)*pH] -> decreases by Nernstian shift
    expected_nernst = KB_EV_PER_K * 298.15 * np.log(10.0) * 14.0
    g_she_ox = apply_che_correction(0.0, u_potential=0.0, ph=14.0, nu_H=-1.0, reference_electrode="SHE", temperature_k=298.15)
    assert pytest.approx(g_she_ox, rel=1e-4) == -expected_nernst

    # 4. Explicit signed stoichiometry nu_H = +1 (reduction) vs SHE at pH = 14:
    # Delta G = Delta G0 + [U_SHE + kB*T*ln(10)*pH] -> increases by Nernstian shift
    g_she_red = apply_che_correction(0.0, u_potential=0.0, ph=14.0, nu_H=+1.0, reference_electrode="SHE", temperature_k=298.15)
    assert pytest.approx(g_she_red, rel=1e-4) == expected_nernst

    # 5. Invalid step_type or reference_electrode
    with pytest.raises(ValueError):
        apply_che_correction(0.0, step_type="unknown")
    with pytest.raises(ValueError):
        apply_che_correction(0.0, reference_electrode="NHE_UNSUPPORTED")


def test_convert_reference_potential():
    """Verify Nernstian potential conversions between SHE and RHE."""
    # At pH 0, U_RHE = U_SHE
    assert pytest.approx(convert_reference_potential(0.5, "SHE", "RHE", ph=0.0)) == 0.5

    # At pH 14, 298.15 K: U_RHE = U_SHE + 0.05916 * 14 = U_SHE + 0.82824 V
    nernst_val = KB_EV_PER_K * 298.15 * np.log(10.0) * 14.0
    u_rhe = convert_reference_potential(0.0, "SHE", "RHE", ph=14.0, temperature_k=298.15)
    assert pytest.approx(u_rhe, rel=1e-4) == nernst_val

    # Reverse conversion: RHE -> SHE
    u_she = convert_reference_potential(u_rhe, "RHE", "SHE", ph=14.0, temperature_k=298.15)
    assert pytest.approx(u_she, abs=1e-5) == 0.0

    with pytest.raises(ValueError):
        convert_reference_potential(0.0, "SHE", "INVALID")


def test_plot_reaction_pathway():
    """Verify reaction energy diagram plotting and validation."""
    steps = ["R", "Int1", "TS", "Int2", "P"]
    pathways = {
        "Catalyst_A": [0.0, -0.4, 0.8, -0.2, -1.2],
        "Catalyst_B": [0.0, 0.2, 1.4, 0.5, -1.2],
    }
    ts_flags = {
        "Catalyst_A": [False, False, True, False, False],
        "Catalyst_B": [False, False, True, False, False],
    }

    fig, ax = create_publication_figure(columns="double")
    ax_out = plot_reaction_pathway(
        pathways=pathways,
        step_labels=steps,
        transition_states=ts_flags,
        relative_to_initial=True,
        ax=ax,
    )
    assert ax_out is ax
    assert len(ax.get_xticks()) == len(steps)
    plt.close(fig)

    # Dimension mismatch check
    with pytest.raises(ValueError):
        plot_reaction_pathway(
            pathways={"Failing": [0.0, 1.0]},  # Only 2 steps instead of 5
            step_labels=steps,
        )


def test_evaluate_volcano_model():
    """Verify Sabatier two-leg volcano curve and apex location."""
    x = np.linspace(0.0, 2.0, 100)
    # Left leg: y = 1.0 * x; Right leg: y = -1.0 * x + 2.0
    # Intersection at x = 1.0, y = 1.0
    y_curve, (apex_x, apex_y) = evaluate_volcano_model(x, 1.0, 0.0, -1.0, 2.0)
    assert pytest.approx(apex_x) == 1.0
    assert pytest.approx(apex_y) == 1.0
    assert np.max(y_curve) <= 1.0

    with pytest.raises(ValueError):
        evaluate_volcano_model(x, 1.0, 0.0, 1.0, 2.0)  # Parallel slopes


def test_plot_volcano():
    """Verify volcano plot with candidates and model."""
    descriptors = [0.4, 0.7, 1.0, 1.3, 1.6]
    activities = [0.4, 0.7, 0.95, 0.65, 0.35]
    names = ["Ni", "Pd", "Pt", "Cu", "Au"]

    fig, ax = create_publication_figure()
    ax_out, apex = plot_volcano(
        descriptors=descriptors,
        activities=activities,
        catalyst_names=names,
        theoretical_model=(1.0, 0.0, -1.0, 2.0),
        highlight_catalysts=["Pt"],
        ax=ax,
    )
    assert ax_out is ax
    assert apex is not None
    assert pytest.approx(apex[0]) == 1.0
    plt.close(fig)


def test_evaluate_orr_volcano():
    """Verify 4-electron associative ORR volcano model and Sabatier apex."""
    dG_oh = np.array([0.4, 0.7, 0.86, 1.1, 1.4])
    u_lim, eta, apex = evaluate_orr_volcano(dG_oh)

    # Apex should be at dG_OH = (4.92 - 3.20) / 2.0 = 0.86 eV, U_L = 0.86 V
    assert pytest.approx(apex[0], abs=1e-3) == 0.86
    assert pytest.approx(apex[1], abs=1e-3) == 0.86

    # At apex, theoretical overpotential eta = 1.23 - 0.86 = 0.37 V
    assert pytest.approx(1.23 - apex[1], abs=1e-3) == 0.37

    # At dG_OH = 0.86 eV, U_L should be maximized
    idx_apex = np.argmin(np.abs(dG_oh - 0.86))
    assert u_lim[idx_apex] == np.max(u_lim)


def test_evaluate_her_volcano():
    """Verify 2-electron HER volcano model with optimum at Delta G_H* = 0 eV."""
    dG_h = np.array([-0.4, -0.1, 0.0, 0.1, 0.4])
    log_j0, eta, apex = evaluate_her_volcano(dG_h, temperature_k=298.15, log_j0_max=0.0)

    # Apex should be at Delta G_H* = 0.0 eV, log_j0 = 0.0
    assert pytest.approx(apex[0]) == 0.0
    assert pytest.approx(apex[1]) == 0.0

    # Overpotentials are |Delta G_H*|
    assert pytest.approx(eta[0]) == 0.4
    assert pytest.approx(eta[2]) == 0.0

    # log_j0 is maximized at Delta G_H* = 0.0
    assert log_j0[2] == np.max(log_j0)


def test_plot_orr_and_her_volcanoes():
    """Verify specialized ORR and HER volcano plotting routines."""
    fig, ax = create_publication_figure()
    ax_orr, apex_orr = plot_orr_volcano(
        descriptors_dG_OH=[0.5, 0.86, 1.2],
        catalyst_names=["Fe", "Pt", "Au"],
        highlight_catalysts=["Pt"],
        ax=ax,
    )
    assert ax_orr is ax
    assert pytest.approx(apex_orr[0], abs=1e-3) == 0.86
    plt.close(fig)

    fig2, ax2 = create_publication_figure()
    ax_her, apex_her = plot_her_volcano(
        descriptors_dG_H=[-0.3, -0.05, 0.25],
        catalyst_names=["Ni", "Pt", "Au"],
        highlight_catalysts=["Pt"],
        ax=ax2,
    )
    assert ax_her is ax2
    assert pytest.approx(apex_her[0]) == 0.0
    plt.close(fig2)


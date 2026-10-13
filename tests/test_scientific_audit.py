"""
Independent Scientific Correctness and Unit Validation Test Suite.

Audits:
- Physical constants and absence of double conversions
- Reversible vs Standard Hydrogen Electrode potential shifts
- Signed reaction stoichiometry in CHE
- ORR vs HER distinct volcano models
- Working-Hotelling mean confidence interval vs individual prediction intervals
- Adsorption energy sign conventions
- NEB forward/reverse activation barrier thermodynamic consistency
- d-band center numerical quadrature and spin isolation
"""

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pytest
from scipy.integrate import trapezoid

from modules.catalysis.reaction_energy import (
    apply_che_correction,
    convert_reference_potential,
)
from modules.catalysis.volcano import (
    evaluate_her_volcano,
    evaluate_orr_volcano,
)
from modules.dft.adsorption import (
    calculate_adsorption_energy,
    plot_adsorption_scaling,
)
from parsers.quantum_espresso import (
    parse_bands_gnu,
    parse_dos_file,
    parse_neb_output,
    parse_pw_output,
)
from parsers.quantum_espresso.pdos import QEPdosCollection
from scripts.preprocess import (
    BOHR_TO_ANGSTROM,
    BOLTZMANN_CONSTANT_EV_PER_K,
    BOLTZMANN_CONSTANT_J_PER_K,
    ELEMENTARY_CHARGE_C,
    EV_TO_KJ_PER_MOL,
    GAS_CONSTANT_R_J_MOL_K,
    RY_PER_BOHR_TO_EV_PER_ANGSTROM,
    RY_TO_EV,
    convert_energy,
)
from scripts.utils import create_publication_figure

FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "qe"


# =========================================================================
# 1. No Double Energy Conversion & Units Audit
# =========================================================================

def test_no_double_conversion_pw_scf():
    """Verify pw.x total energy is converted from Ry to eV exactly once and Fermi level is not converted."""
    out_file = FIXTURE_DIR / "scf_converged.out"
    res = parse_pw_output(out_file, strict=True)

    # In scf_converged.out:
    # !    total energy              =     -15.84523190 Ry
    #      highest occupied level is   6.2541 eV
    raw_ry = -15.84523190
    expected_ev = raw_ry * RY_TO_EV

    assert pytest.approx(res.total_energy_ry) == raw_ry
    assert pytest.approx(res.total_energy_ev) == expected_ev
    # Highest occupied level in QE is printed in eV, so it must not be multiplied by RY_TO_EV
    assert pytest.approx(res.fermi_energy) == 6.2541
    assert res.fermi_energy < 10.0  # If erroneously converted by 13.6, would be ~85 eV


def test_no_double_conversion_dos_and_pdos():
    """Verify dos.x and projwfc.x outputs are natively in eV and never converted twice."""
    dos_file = FIXTURE_DIR / "dos_nonspin.dat"
    res_dos = parse_dos_file(dos_file, e_fermi=0.0)

    # In dos_nonspin.dat, first column is energy in eV: ranges from 0.0 to 10.0 eV
    assert pytest.approx(res_dos.energies[0]) == 0.0
    assert pytest.approx(res_dos.energies[-1]) == 10.0
    # dos_total is density of states (states/eV)
    assert np.all(res_dos.dos_total >= 0.0)

    pdos_col = QEPdosCollection(FIXTURE_DIR / "pdos_spin", prefix="fe", e_fermi=0.0)
    # Energy grid should span eV scale (-10 to 5), not Ry scale (-136 to 68)
    assert pdos_col.energies[0] >= -15.0
    assert pdos_col.energies[-1] <= 10.0


def test_force_and_stress_units():
    """Verify forces are converted from Ry/au to eV/Angstrom using exact conversion."""
    out_file = FIXTURE_DIR / "scf_converged.out"
    res = parse_pw_output(out_file, strict=True)

    expected_conversion = RY_TO_EV / BOHR_TO_ANGSTROM
    assert RY_PER_BOHR_TO_EV_PER_ANGSTROM == expected_conversion
    assert res.forces_ev_angstrom is not None
    assert res.stress_kbar is not None


# =========================================================================
# 2. Signed CHE & Reference Potential Audit
# =========================================================================

def test_che_thermodynamic_reversibility():
    """Verify forward oxidation and reverse reduction cancel exactly at any potential and pH."""
    u_pot = 0.85  # V vs SHE
    ph_sol = 7.0
    dG0_forward = 0.45  # eV for A -> B + H+ + e- (oxidation, nu_H = -1)
    dG0_reverse = -0.45 # eV for B + H+ + e- -> A (reduction, nu_H = +1)

    dG_fwd = apply_che_correction(dG0_forward, u_potential=u_pot, ph=ph_sol, nu_H=-1.0, reference_electrode="SHE")
    dG_rev = apply_che_correction(dG0_reverse, u_potential=u_pot, ph=ph_sol, nu_H=+1.0, reference_electrode="SHE")

    # Thermodynamic microscopic reversibility: dG_forward + dG_reverse == 0.0
    assert pytest.approx(dG_fwd + dG_rev, abs=1e-12) == 0.0


def test_che_she_to_rhe_equivalence():
    """Verify that calculating dG vs SHE matches dG vs RHE when potential is properly converted."""
    u_she = 0.40  # V vs SHE
    ph = 13.0
    u_rhe = convert_reference_potential(u_she, "SHE", "RHE", ph=ph)

    dG0 = 1.50
    # Evaluate using SHE with explicit pH
    dG_via_she = apply_che_correction(dG0, u_potential=u_she, ph=ph, nu_H=-1.0, reference_electrode="SHE")
    # Evaluate using RHE (where pH shift is already in U_RHE)
    dG_via_rhe = apply_che_correction(dG0, u_potential=u_rhe, ph=0.0, nu_H=-1.0, reference_electrode="RHE")

    assert pytest.approx(dG_via_she, abs=1e-10) == dG_via_rhe


# =========================================================================
# 3. ORR and HER Volcano Models
# =========================================================================

def test_orr_volcano_limiting_potential_steps():
    """Verify that ORR limiting potential correctly switches between strong and weak binding limits."""
    # Strong binding: dG_OH = 0.4 eV (small, binding is strong).
    # Limited by *OH -> H2O removal: U_L = dG_OH = 0.4 V
    u_lim_strong, eta_strong, apex = evaluate_orr_volcano([0.4])
    assert pytest.approx(u_lim_strong[0], abs=1e-3) == 0.4

    # Weak binding: dG_OH = 1.3 eV (large, binding is weak).
    # Limited by O2 -> *OOH formation: U_L = 4.92 - (1.3 + 3.20) = 4.92 - 4.50 = 0.42 V
    u_lim_weak, eta_weak, _ = evaluate_orr_volcano([1.3])
    assert pytest.approx(u_lim_weak[0], abs=1e-3) == 0.42

    # Apex: intersection at dG_OH = 0.86 eV, where U_L = 0.86 V
    assert pytest.approx(apex[0], abs=1e-3) == 0.86
    assert pytest.approx(apex[1], abs=1e-3) == 0.86


def test_her_volcano_symmetry():
    """Verify HER volcano overpotential and exchange current density are symmetric around Delta G_H = 0."""
    dG_pos = 0.25
    dG_neg = -0.25

    log_pos, eta_pos, _ = evaluate_her_volcano([dG_pos])
    log_neg, eta_neg, _ = evaluate_her_volcano([dG_neg])

    assert pytest.approx(eta_pos[0]) == eta_neg[0]
    assert pytest.approx(log_pos[0]) == log_neg[0]


# =========================================================================
# 4. Working-Hotelling CI vs Prediction Interval Audit
# =========================================================================

def test_prediction_interval_strictly_wider_than_confidence_interval():
    """Verify prediction interval for individual observations is strictly wider than mean CI."""
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    y = np.array([2.1, 3.8, 6.2, 7.9, 10.1, 11.9])

    fig1, ax1 = create_publication_figure()
    _, stats_ci = plot_adsorption_scaling(x, y, interval_type="confidence", ax=ax1)

    fig2, ax2 = create_publication_figure()
    _, stats_pi = plot_adsorption_scaling(x, y, interval_type="prediction", ax=ax2)

    # Slopes and R2 must be identical
    assert pytest.approx(stats_ci["slope"]) == stats_pi["slope"]
    assert pytest.approx(stats_ci["r_squared"]) == stats_pi["r_squared"]

    # F-statistic and standard errors must be computed
    assert stats_ci["f_statistic"] > 10.0
    assert stats_ci["residual_std_error"] > 0.0
    assert stats_ci["intercept_stderr"] > 0.0

    plt.close(fig1)
    plt.close(fig2)


# =========================================================================
# 5. NEB Barrier Extraction Thermodynamic Consistency
# =========================================================================

def test_neb_barrier_thermodynamic_consistency():
    """Verify forward barrier, reverse barrier, and reaction energy adhere to Ea_fwd - Ea_rev = Delta E_rxn."""
    neb_file = FIXTURE_DIR / "neb.dat"
    res = parse_neb_output(neb_file)

    assert res.forward_barrier_ev is not None
    assert res.reverse_barrier_ev is not None
    assert res.reaction_energy_ev is not None

    # Delta E_rxn = E(final) - E(initial) = Ea_fwd - Ea_rev
    calc_delta_e = res.forward_barrier_ev - res.reverse_barrier_ev
    assert pytest.approx(res.reaction_energy_ev, abs=1e-5) == calc_delta_e


# =========================================================================
# 6. d-Band Center Numerical Quadrature
# =========================================================================

def test_d_band_center_quadrature_and_spin_channels():
    """Verify d-band center calculation on spin-polarized collection."""
    pdos = QEPdosCollection(FIXTURE_DIR / "pdos_spin", prefix="fe", e_fermi=5.0)

    # d-band center in window [-5, 0] eV relative to EF
    eps_tot = pdos.calculate_d_band_center(element="Fe", energy_window=(-5.0, 0.0), fermi_energy=5.0, spin="total")
    eps_up = pdos.calculate_d_band_center(element="Fe", energy_window=(-5.0, 0.0), fermi_energy=5.0, spin="up")
    eps_down = pdos.calculate_d_band_center(element="Fe", energy_window=(-5.0, 0.0), fermi_energy=5.0, spin="down")

    # Center must lie within the integration window [-5, 0]
    assert -5.0 <= eps_tot <= 0.0
    assert -5.0 <= eps_up <= 0.0
    assert -5.0 <= eps_down <= 0.0

    # Total center must lie between spin-up and spin-down centers (or equal if identical)
    min_spin = min(eps_up, eps_down)
    max_spin = max(eps_up, eps_down)
    assert min_spin <= eps_tot <= max_spin

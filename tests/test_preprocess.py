"""
Tests for scripts/preprocess.py unit conversions and preprocessing workflows.
"""

import numpy as np
import pytest

from scripts.preprocess import (
    AVOGADRO_CONSTANT_MOL_INV,
    BOHR_TO_ANGSTROM,
    BOLTZMANN_CONSTANT_EV_PER_K,
    BOLTZMANN_CONSTANT_J_PER_K,
    ELEMENTARY_CHARGE_C,
    EV_TO_J,
    EV_TO_KCAL_PER_MOL,
    EV_TO_KJ_PER_MOL,
    GAS_CONSTANT_R_J_MOL_K,
    GAS_CONSTANT_R_KJ_MOL_K,
    KCAL_PER_MOL_TO_EV,
    RY_PER_BOHR_TO_EV_PER_ANGSTROM,
    RY_TO_EV,
    RY_TO_EV_CODATA2018,
    RY_TO_EV_QE,
    convert_energy,
    convert_temperature,
)


def test_fundamental_constants_si2019_consistency():
    """Verify fundamental constants adhere to exact 2019 SI definitions and CODATA 2018."""
    # Elementary charge and Boltzmann constant
    assert ELEMENTARY_CHARGE_C == 1.602176634e-19
    assert AVOGADRO_CONSTANT_MOL_INV == 6.02214076e23
    assert BOLTZMANN_CONSTANT_J_PER_K == 1.380649e-23

    # Gas constant R = NA * kB exact
    expected_r = AVOGADRO_CONSTANT_MOL_INV * BOLTZMANN_CONSTANT_J_PER_K
    assert GAS_CONSTANT_R_J_MOL_K == expected_r
    assert GAS_CONSTANT_R_KJ_MOL_K == expected_r / 1000.0

    # kB in eV/K: kB(J/K) / e
    expected_kb_ev = BOLTZMANN_CONSTANT_J_PER_K / ELEMENTARY_CHARGE_C
    assert BOLTZMANN_CONSTANT_EV_PER_K == expected_kb_ev

    # 1 eV in kJ/mol: (e * NA) / 1000
    expected_kj_mol = (ELEMENTARY_CHARGE_C * AVOGADRO_CONSTANT_MOL_INV) / 1000.0
    assert pytest.approx(EV_TO_KJ_PER_MOL, rel=1e-12) == expected_kj_mol

    # 1 kcal/mol: 4184 J/mol / (e * NA) in eV
    expected_kcal_to_ev = 4184.0 / (ELEMENTARY_CHARGE_C * AVOGADRO_CONSTANT_MOL_INV)
    assert pytest.approx(KCAL_PER_MOL_TO_EV, rel=1e-12) == expected_kcal_to_ev
    assert pytest.approx(EV_TO_KCAL_PER_MOL * KCAL_PER_MOL_TO_EV, rel=1e-12) == 1.0

    # Bohr and Rydberg conversions
    assert BOHR_TO_ANGSTROM == 0.529177210903
    assert RY_TO_EV_QE == 13.605693009
    assert pytest.approx(RY_TO_EV_CODATA2018, rel=1e-6) == 13.60569312
    assert RY_PER_BOHR_TO_EV_PER_ANGSTROM == RY_TO_EV / BOHR_TO_ANGSTROM


def test_convert_energy_exact():
    """Verify energy conversions against CODATA constants and round-trips."""
    assert pytest.approx(convert_energy(1.0, "eV", "kJ/mol")) == EV_TO_KJ_PER_MOL
    assert pytest.approx(convert_energy(EV_TO_KJ_PER_MOL, "kJ/mol", "eV")) == 1.0

    # Ry to eV round trip
    assert pytest.approx(convert_energy(1.0, "Ry", "eV")) == RY_TO_EV
    assert pytest.approx(convert_energy(RY_TO_EV, "eV", "Ry")) == 1.0

    # kcal/mol round trip
    val_kcal = convert_energy(1.0, "eV", "kcal/mol")
    assert pytest.approx(convert_energy(val_kcal, "kcal/mol", "eV")) == 1.0

    # Array input
    arr = np.array([0.0, 1.0, 2.0])
    conv = convert_energy(arr, "eV", "kJ/mol")
    assert pytest.approx(conv[1]) == EV_TO_KJ_PER_MOL
    assert pytest.approx(conv[2]) == 2.0 * EV_TO_KJ_PER_MOL


def test_convert_energy_invalid():
    """Verify error on unsupported unit."""
    with pytest.raises(ValueError):
        convert_energy(1.0, "calories", "eV")
    with pytest.raises(ValueError):
        convert_energy(1.0, "eV", "hartree")


def test_convert_temperature():
    """Verify temperature conversions."""
    assert pytest.approx(convert_temperature(0.0, "C", "K")) == 273.15
    assert pytest.approx(convert_temperature(298.15, "K", "C")) == 25.0

    with pytest.raises(ValueError):
        convert_temperature(100.0, "F", "C")


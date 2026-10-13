"""
Comprehensive unit tests for experimental electrochemistry analysis modules.
"""

from pathlib import Path
import pytest
import numpy as np

from modules.electrochemistry.benchmarking import compare_instrument_lsv
from modules.electrochemistry.chronoamperometry import calculate_ca_stability_metrics
from modules.electrochemistry.chronopotentiometry import calculate_cp_stability_metrics
from modules.electrochemistry.cv import calculate_voltammetric_charge, extract_redox_peaks
from modules.electrochemistry.ecsa import (
    CdlInterpretationWarning,
    calculate_double_layer_capacitance,
    estimate_ecsa,
)
from modules.electrochemistry.eis import (
    calculate_impedance_polar,
    fit_equivalent_circuit,
    run_kramers_kronig_test,
)
from modules.electrochemistry.lsv import extract_overpotential_at_current
from modules.electrochemistry.reference_electrodes import (
    apply_ir_compensation,
    calculate_overpotential,
    convert_potential_to_rhe,
    get_standard_electrode_potential,
)
from modules.electrochemistry.tafel import (
    NonSteadyStateWarning,
    calculate_tafel_slope,
)
from parsers.electrochemistry import load_electrochemistry_data


FIXTURES_DIR = Path(__file__).parent / "fixtures" / "electrochemistry"


# =========================================================================
# 1. Reference Electrodes and Thermodynamic Conversions
# =========================================================================

def test_reference_electrode_potentials_298k():
    """Verify standard potentials vs SHE at 298.15 K."""
    assert get_standard_electrode_potential("she") == pytest.approx(0.000)
    assert get_standard_electrode_potential("sce") == pytest.approx(0.2412)
    assert get_standard_electrode_potential("ag/agcl_sat") == pytest.approx(0.1970)
    assert get_standard_electrode_potential("ag/agcl_3m") == pytest.approx(0.2100)
    assert get_standard_electrode_potential("hg/hgo_1m") == pytest.approx(0.0980)


def test_convert_potential_to_rhe():
    """
    Test Nernstian conversion to RHE:
    E_RHE = E_meas + E_ref^0 + 0.05916 * pH
    For Ag/AgCl (sat KCl, E0=0.1970) at pH 14, 298.15 K:
    E_RHE = 0.500 + 0.1970 + 0.059155 * 14 = 0.500 + 0.1970 + 0.82817 = 1.52517 V
    """
    v_meas = np.array([0.0, 0.500, 1.000])
    v_rhe = convert_potential_to_rhe(v_meas, reference_electrode="ag/agcl", ph=14.0)
    expected_0_5 = 0.500 + 0.1970 + (8.314462618 * 298.15 * np.log(10) / 96485.33212) * 14.0
    assert v_rhe[1] == pytest.approx(expected_0_5, rel=1e-4)


def test_temperature_correction_reference_electrode():
    """Test dE/dT temperature correction for SCE."""
    # At T = 308.15 K (dT = +10 K), SCE dE/dT = -0.65 mV/K -> -6.5 mV
    e_std_298 = get_standard_electrode_potential("sce", temperature_k=298.15)
    e_std_308 = get_standard_electrode_potential("sce", temperature_k=308.15)
    assert e_std_308 == pytest.approx(e_std_298 - 0.0065, rel=1e-3)


def test_ir_compensation_and_double_compensation_guard():
    """Test iR compensation formula and safety guard against duplicate correction."""
    v_raw = np.array([1.50, 1.60])
    i_raw = np.array([0.010, 0.020])  # 10 mA, 20 mA (+ anodic)
    ru = 5.0  # 5 Ohm

    # E_corr = E - (I * Ru * fraction) = 1.60 - (0.020 * 5.0 * 1.0) = 1.60 - 0.10 = 1.50 V
    v_corr = apply_ir_compensation(v_raw, i_raw, ru_ohm=ru, compensation_fraction=1.0)
    assert v_corr[1] == pytest.approx(1.50)

    # Test guard on schema object
    data = load_electrochemistry_data(FIXTURES_DIR / "metrohm" / "metrohm_oer_lsv.txt")
    data.ir_compensated = True
    with pytest.raises(ValueError, match="already marked as iR-compensated"):
        data.get_ir_corrected_potential(ru_ohm=4.0)


# =========================================================================
# 2. Cyclic Voltammetry (CV)
# =========================================================================

def test_extract_redox_peaks_and_charge():
    """Test oxidation and reduction peak extraction and voltammetric charge."""
    data = load_electrochemistry_data(FIXTURES_DIR / "metrohm" / "metrohm_cv_multicycle.csv")
    cycle_1 = data.get_cycle(1)

    peaks = extract_redox_peaks(cycle_1.potential_V, cycle_1.current_A)
    assert peaks["anodic_peak_potential_V"] is not None
    assert peaks["cathodic_peak_potential_V"] is not None
    # Oxidation peak is around 0.48 V, reduction around 0.38 V
    assert peaks["anodic_peak_potential_V"] == pytest.approx(0.48, abs=0.03)
    assert peaks["cathodic_peak_potential_V"] == pytest.approx(0.38, abs=0.03)
    assert peaks["delta_ep_mV"] == pytest.approx(100.0, abs=40.0)

    # Voltammetric charge
    q = calculate_voltammetric_charge(cycle_1.time_s, cycle_1.current_A)
    assert q["total_charge_C"] > 0
    assert q["anodic_charge_C"] > 0
    assert q["cathodic_charge_C"] > 0


# =========================================================================
# 3. Linear Sweep Voltammetry (LSV) & Overpotentials
# =========================================================================

def test_extract_overpotential_at_current():
    """Test eta_10 extraction for OER and HER."""
    # OER dataset
    data_oer = load_electrochemistry_data(FIXTURES_DIR / "metrohm" / "metrohm_oer_lsv.txt")
    pot_rhe = data_oer.get_potential_vs_rhe(reference="ag/agcl", ph=14.0)
    j = data_oer.get_current_density(unit="mA/cm^2")

    eta_10 = extract_overpotential_at_current(pot_rhe, j, target_j=10.0, reaction="OER")
    assert eta_10 is not None
    # In alkaline OER on NiFe-LDH, eta_10 is typically 200 - 550 mV
    assert 200.0 < eta_10 < 600.0

    # HER dataset
    data_her = load_electrochemistry_data(FIXTURES_DIR / "biologic" / "biologic_her_lsv.mpt")
    pot_her_rhe = data_her.get_potential_vs_rhe(reference="ag/agcl", ph=0.3)
    j_her = data_her.get_current_density(unit="mA/cm^2")
    eta_10_her = extract_overpotential_at_current(pot_her_rhe, j_her, target_j=10.0, reaction="HER")
    assert eta_10_her is not None
    assert 30.0 < eta_10_her < 250.0


# =========================================================================
# 4. Tafel Analysis & Kinetics
# =========================================================================

def test_calculate_tafel_slope_and_warning():
    """Test Tafel slope calculation and non-steady-state warning at high scan rates."""
    data = load_electrochemistry_data(FIXTURES_DIR / "metrohm" / "metrohm_oer_lsv.txt")
    pot_rhe = data.get_potential_vs_rhe(reference="ag/agcl", ph=14.0)
    j = data_oer = data.get_current_density(unit="mA/cm^2")

    # Fit in the kinetic region (j between 1.0 and 10.0 mA/cm^2)
    res = calculate_tafel_slope(
        pot_rhe,
        j,
        current_range=(1.0, 10.0),
        reaction="OER",
        scan_rate_V_s=0.005,  # 5 mV/s (quasi-steady-state)
    )

    assert 30.0 < res["tafel_slope_mV_dec"] < 160.0
    assert res["r_squared"] > 0.95
    assert res["exchange_current_density_mA_cm2"] is not None

    # Test that scan rate > 5 mV/s triggers NonSteadyStateWarning
    with pytest.warns(NonSteadyStateWarning):
        calculate_tafel_slope(
            pot_rhe,
            j,
            current_range=(1.0, 10.0),
            scan_rate_V_s=0.050,  # 50 mV/s
        )


# =========================================================================
# 5. Double Layer Capacitance (C_dl) & ECSA
# =========================================================================

def test_calculate_cdl_and_ecsa():
    """Test double-layer capacitance linear regression and ECSA."""
    # Synthetic scan rates (mV/s) and non-Faradaic delta_j/2 (mA/cm^2)
    # C_dl = 0.50 mF/cm^2 = 500 uF/cm^2
    scan_rates = np.array([10.0, 20.0, 40.0, 60.0, 80.0, 100.0])  # mV/s
    cdl_true = 0.50  # mF/cm^2
    charging_currents = cdl_true * (scan_rates * 1e-3) * 1000.0  # in uA/cm^2

    res = calculate_double_layer_capacitance(scan_rates, charging_currents, current_unit="uA/cm2")
    assert res["c_dl_uF_cm2"] == pytest.approx(500.0, rel=1e-3)
    assert res["r_squared"] > 0.99

    # ECSA estimate with Cs = 0.040 mF/cm^2 = 40 uF/cm^2
    # Roughness factor = 500 / 40 = 12.5
    ecsa_res = estimate_ecsa(res["c_dl_uF_cm2"], cs_uF_cm2=40.0, geometric_area_cm2=0.196)
    assert ecsa_res["roughness_factor"] == pytest.approx(12.5, rel=1e-3)
    assert ecsa_res["ecsa_cm2"] == pytest.approx(12.5 * 0.196, rel=1e-3)

    # Warning when computing active site density from Cdl
    with pytest.warns(CdlInterpretationWarning):
        estimate_ecsa(res["c_dl_uF_cm2"], cs_uF_cm2=40.0, warn_active_sites=True)


# =========================================================================
# 6. Electrochemical Impedance Spectroscopy (EIS)
# =========================================================================

def test_eis_polar_conversion_and_circuit_fitting():
    """Test polar conversion, Randles circuit fitting, and lin-KK validation."""
    data = load_electrochemistry_data(FIXTURES_DIR / "metrohm" / "metrohm_fra_eis.txt")
    f = data.frequency_Hz
    zr = data.impedance_real_ohm
    zi = data.impedance_imag_ohm

    mag, phase = calculate_impedance_polar(zr, zi)
    assert len(mag) == len(f)
    assert len(phase) == len(f)
    assert np.all(mag > 0)

    # Equivalent circuit fit: Rs - [Rct || CPE]
    fit = fit_equivalent_circuit(f, zr, zi, circuit_string="R0-p(R1,CPE1)")
    assert fit["rs_ohm"] is not None
    assert fit["rct_ohm"] is not None
    # Generated with Rs=5.2 Ohm, Rct=45.0 Ohm
    assert fit["rs_ohm"] == pytest.approx(5.2, abs=1.0)
    assert fit["rct_ohm"] == pytest.approx(45.0, abs=5.0)
    assert fit["rmse"] < 5.0

    # Lin-KK validation test
    kk_res = run_kramers_kronig_test(f, zr, zi)
    assert kk_res["valid"] is True
    assert kk_res["num_elements"] > 0


# =========================================================================
# 7. Chronoamperometry & Chronopotentiometry Durability
# =========================================================================

def test_ca_stability_metrics():
    """Test CA retention percentage and degradation rate."""
    t = np.linspace(0, 3600, 361)
    # Starts at 20 mA/cm^2, drifts to 19.2 mA/cm^2 (96% retention)
    j = 20.0 * (1.0 - 0.04 * (t / 3600.0))
    metrics = calculate_ca_stability_metrics(t, j)

    assert metrics["retention_percent"] == pytest.approx(96.0, rel=1e-2)
    assert metrics["degradation_rate_pct_per_h"] == pytest.approx(4.0, abs=0.1)
    assert metrics["duration_hours"] == pytest.approx(1.0)


def test_cp_stability_metrics():
    """Test CP potential drift rate."""
    t = np.linspace(0, 3600, 361)
    # Starts at 1.500 V, drifts to 1.505 V (+5 mV over 1 hour)
    v = 1.500 + 0.005 * (t / 3600.0)
    metrics = calculate_cp_stability_metrics(t, v)

    assert metrics["potential_drift_mV"] == pytest.approx(5.0, abs=0.2)
    assert metrics["drift_rate_mV_per_h"] == pytest.approx(5.0, abs=0.2)


# =========================================================================
# 8. Cross-Instrument Benchmarking (Metrohm vs BioLogic)
# =========================================================================

def test_cross_instrument_benchmarking():
    """Test comparison of matched OER LSV from Metrohm Autolab and BioLogic EC-Lab."""
    d_metrohm = load_electrochemistry_data(FIXTURES_DIR / "metrohm" / "metrohm_oer_lsv.txt")
    d_biologic = load_electrochemistry_data(FIXTURES_DIR / "biologic" / "biologic_oer_lsv.mpt")

    comp = compare_instrument_lsv(d_metrohm, d_biologic, reaction="OER", target_j=10.0)

    assert comp["conditions_match"] is True
    assert comp["eta10_metrohm_mV"] is not None
    assert comp["eta10_biologic_mV"] is not None
    # Both instruments should yield closely matched eta_10 (< 10 mV discrepancy)
    assert abs(comp["delta_eta10_mV"]) < 10.0
    assert comp["rmse_current_density"] < 2.0

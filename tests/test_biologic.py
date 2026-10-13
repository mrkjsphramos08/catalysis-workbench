"""
Unit tests for BioLogic (EC-Lab) electrochemical data parser.
"""

from pathlib import Path
import pytest
import numpy as np

from parsers.electrochemistry import load_electrochemistry_data
from parsers.electrochemistry.biologic import (
    BioLogicBinaryMprError,
    BioLogicParseError,
    parse_biologic_mpr,
    parse_biologic_mps,
    parse_biologic_mpt,
)
from parsers.electrochemistry.schema import StandardizedElectrochemicalData


FIXTURES_DIR = Path(__file__).parent / "fixtures" / "electrochemistry" / "biologic"


def test_parse_biologic_oer_lsv_mpt():
    """Test reading BioLogic EC-Lab ASCII .mpt export for OER."""
    filepath = FIXTURES_DIR / "biologic_oer_lsv.mpt"
    data = parse_biologic_mpt(filepath)

    assert isinstance(data, StandardizedElectrochemicalData)
    assert data.instrument == "biologic"
    assert data.instrument_software == "EC-Lab"
    assert data.technique == "LSV"
    assert len(data.potential_V) == 260
    assert len(data.current_A) == 260
    # BioLogic I/mA was converted to A
    assert np.max(data.current_A) < 0.1  # In Amperes
    assert data.electrode_area_cm2 == pytest.approx(0.196, rel=1e-3)
    assert data.potential_reference_raw == "Ag/AgCl (3 M KCl)"

    # Current density in mA/cm^2
    j = data.get_current_density(unit="mA/cm^2")
    assert len(j) == 260
    assert np.max(j) > 10.0


def test_parse_biologic_her_lsv_mpt():
    """Test reading BioLogic HER cathodic LSV."""
    filepath = FIXTURES_DIR / "biologic_her_lsv.mpt"
    data = parse_biologic_mpt(filepath)

    assert data.technique == "LSV"
    assert len(data.potential_V) == 225
    # Cathodic currents should be negative under IUPAC convention
    assert np.min(data.current_A) < 0.0


def test_parse_biologic_peis_mpt():
    """Test reading BioLogic PEIS impedance spectrum."""
    filepath = FIXTURES_DIR / "biologic_peis.mpt"
    data = parse_biologic_mpt(filepath)

    assert data.technique == "EIS"
    assert data.frequency_Hz is not None
    assert data.impedance_real_ohm is not None
    assert data.impedance_imag_ohm is not None
    assert len(data.frequency_Hz) == 60
    assert np.all(data.frequency_Hz > 0)
    assert np.all(data.impedance_real_ohm > 0)
    assert np.all(data.impedance_imag_ohm < 0)


def test_parse_biologic_ca_and_cp():
    """Test reading BioLogic CA and CP durability datasets."""
    ca_path = FIXTURES_DIR / "biologic_stability_ca.mpt"
    data_ca = parse_biologic_mpt(ca_path)
    assert data_ca.technique == "CA"
    assert len(data_ca.time_s) == 361

    cp_path = FIXTURES_DIR / "biologic_stability_cp.mpt"
    data_cp = parse_biologic_mpt(cp_path)
    assert data_cp.technique == "CP"
    assert len(data_cp.time_s) == 361


def test_parse_biologic_mps_settings():
    """Test parsing EC-Lab .mps settings file."""
    mps_path = FIXTURES_DIR / "biologic_settings.mps"
    settings = parse_biologic_mps(mps_path)

    assert isinstance(settings, dict)
    assert settings.get("Technique") == "LSV"
    assert float(settings.get("Electrode surface area (cm2)", 0)) == pytest.approx(0.196)


def test_biologic_mpr_handling_nonexistent_or_invalid():
    """Test that invalid or missing .mpr files raise expected exceptions."""
    with pytest.raises(FileNotFoundError):
        parse_biologic_mpr("nonexistent_file.mpr")

    # If an empty/corrupted file is passed as .mpr
    bad_mpr = FIXTURES_DIR / "corrupt.mpr"
    bad_mpr.write_bytes(b"NOT_A_VALID_MPR_FILE")
    try:
        with pytest.raises((BioLogicBinaryMprError, BioLogicParseError, Exception)):
            parse_biologic_mpr(bad_mpr)
    finally:
        if bad_mpr.exists():
            bad_mpr.unlink()


def test_unified_dispatcher_loads_biologic():
    """Test load_electrochemistry_data auto-detects BioLogic files."""
    filepath = FIXTURES_DIR / "biologic_oer_lsv.mpt"
    data = load_electrochemistry_data(filepath)
    assert data.instrument == "biologic"
    assert data.technique == "LSV"

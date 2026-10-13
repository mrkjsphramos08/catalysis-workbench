"""
Unit tests for Metrohm Autolab (NOVA) electrochemical data parser.
"""

from pathlib import Path
import pytest
import numpy as np

from parsers.electrochemistry import load_electrochemistry_data
from parsers.electrochemistry.metrohm import (
    MetrohmBinaryNoxError,
    MetrohmParseError,
    parse_metrohm_export,
)
from parsers.electrochemistry.schema import MissingChannelError, StandardizedElectrochemicalData


FIXTURES_DIR = Path(__file__).parent / "fixtures" / "electrochemistry" / "metrohm"


def test_parse_metrohm_oer_lsv_txt():
    """Test reading NOVA ASCII/TXT export for OER LSV."""
    filepath = FIXTURES_DIR / "metrohm_oer_lsv.txt"
    data = parse_metrohm_export(filepath)

    assert isinstance(data, StandardizedElectrochemicalData)
    assert data.instrument == "metrohm"
    assert data.instrument_software == "NOVA"
    assert data.technique == "LSV"
    assert len(data.potential_V) == 260
    assert len(data.current_A) == 260
    assert data.electrode_area_cm2 == pytest.approx(0.196, rel=1e-3)
    assert "Ag/AgCl" in data.potential_reference_raw

    # Current density normalization
    j = data.get_current_density(unit="mA/cm^2")
    assert len(j) == 260
    assert np.max(j) > 10.0  # OER onset reached


def test_parse_metrohm_cv_multicycle_csv():
    """Test reading NOVA CSV export with multiple cycles."""
    filepath = FIXTURES_DIR / "metrohm_cv_multicycle.csv"
    data = parse_metrohm_export(filepath)

    assert data.instrument == "metrohm"
    assert data.technique == "CV"
    assert data.cycle_number is not None
    assert data.unique_cycles == [1, 2, 3]

    # Test cycle slicing
    cycle_2 = data.get_cycle(2)
    assert cycle_2.technique == "CV"
    assert len(cycle_2.potential_V) == 200
    assert np.all(cycle_2.cycle_number == 2)


def test_parse_metrohm_fra_eis():
    """Test reading NOVA FRA EIS export."""
    filepath = FIXTURES_DIR / "metrohm_fra_eis.txt"
    data = parse_metrohm_export(filepath)

    assert data.technique == "EIS"
    assert data.frequency_Hz is not None
    assert data.impedance_real_ohm is not None
    assert data.impedance_imag_ohm is not None
    assert len(data.frequency_Hz) == 60
    assert np.all(data.frequency_Hz > 0)
    # Z' must be positive for standard electrochemical cell
    assert np.all(data.impedance_real_ohm > 0)
    # Physical Z'' is negative for capacitive impedance
    assert np.all(data.impedance_imag_ohm < 0)


def test_parse_metrohm_stability_ca_and_cp():
    """Test reading CA and CP files with constant potential/current."""
    ca_path = FIXTURES_DIR / "metrohm_stability_ca.txt"
    data_ca = parse_metrohm_export(ca_path)
    assert data_ca.technique == "CA"
    assert len(data_ca.time_s) == 361

    cp_path = FIXTURES_DIR / "metrohm_stability_cp.txt"
    data_cp = parse_metrohm_export(cp_path)
    assert data_cp.technique == "CP"
    assert len(data_cp.time_s) == 361


def test_reject_binary_nox_file():
    """Test that proprietary binary .nox project files are rejected with explicit guidance."""
    nox_path = FIXTURES_DIR / "mock_procedure.nox"
    with pytest.raises(MetrohmBinaryNoxError) as exc_info:
        parse_metrohm_export(nox_path)

    err_msg = str(exc_info.value)
    assert "proprietary binary" in err_msg.lower() or ".nox" in err_msg.lower()
    assert "export" in err_msg.lower()


def test_custom_column_mapping(tmp_path):
    """Test custom column mapping override for non-standard NOVA headers."""
    custom_file = tmp_path / "custom_nova.txt"
    custom_file.write_text(
        "# Custom header\n"
        "MyVoltage\tMyCurrent\tSecs\n"
        "0.100\t0.001\t0.0\n"
        "0.200\t0.005\t1.0\n",
        encoding="utf-8",
    )

    data = parse_metrohm_export(
        custom_file,
        custom_mapping={
            "MyVoltage": "potential",
            "MyCurrent": "current",
            "Secs": "time",
        },
    )
    assert len(data.potential_V) == 2
    assert len(data.current_A) == 2
    assert len(data.time_s) == 2


def test_missing_required_channels_raises_error(tmp_path):
    """Test that missing required physical channels raise MissingChannelError."""
    bad_file = tmp_path / "invalid_channels.txt"
    bad_file.write_text(
        "Temperature\tHumidity\n"
        "25.0\t50.0\n",
        encoding="utf-8",
    )

    with pytest.raises(MissingChannelError):
        parse_metrohm_export(bad_file)

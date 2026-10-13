"""
Unit and Numerical Validation Tests for ICP-OES / ICP-MS Module.
"""

from pathlib import Path
import numpy as np
import pytest

from parsers.characterization.common import CharacterizationData
from parsers.characterization.icp import parse_icp_file
from modules.characterization.icp import (
    fit_icp_calibration,
    calculate_catalyst_leaching,
    plot_icp_calibration,
    plot_icp_leaching_comparison,
)

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures" / "icp"


def test_parse_icp_standards_and_leaching():
    """Verify parsing of ICP calibration standards and leaching study files."""
    data_cal = parse_icp_file(FIXTURES_DIR / "icp_calibration_standards.csv")
    assert isinstance(data_cal, CharacterizationData)
    assert data_cal.technique == "ICP-OES"
    assert "Pt" in data_cal.metadata["elements_detected"]
    assert "intensity" in data_cal.dataframe.columns

    data_leach = parse_icp_file(FIXTURES_DIR / "catalyst_leaching_study.csv")
    assert len(data_leach.dataframe) == 8


def test_icp_calibration_regression_and_lod():
    """Verify multi-point linear regression and IUPAC LOD estimation."""
    data_cal = parse_icp_file(FIXTURES_DIR / "icp_calibration_standards.csv")
    df_pt = data_cal.dataframe[data_cal.dataframe["element"] == "Pt"]
    # Separate standards from blanks
    df_stds = df_pt[df_pt["sample_id"].str.startswith("STD_")]
    df_blanks = df_pt[df_pt["sample_id"].str.startswith("BLANK_")]

    cal_res = fit_icp_calibration(
        standards_concentration=df_stds["concentration"].values,
        standards_intensity=df_stds["intensity"].values,
        blank_intensities=df_blanks["intensity"].values,
    )

    assert cal_res["r_squared"] > 0.999
    assert cal_res["slope"] == pytest.approx(12500.0, abs=100.0)
    assert cal_res["lod_conc"] is not None
    assert 0.001 < cal_res["lod_conc"] < 0.10


def test_catalyst_leaching_calculation():
    """Verify dissolution mass, normalized loss, and leached percentage."""
    # 2.0 mg of 20 wt% metal catalyst in 50 mL electrolyte
    # Initial metal inventory = 2.0 mg * 0.20 = 0.40 mg = 400 ug
    # Measured concentration = 0.40 mg/L, dilution factor = 2.0 -> net = 0.80 mg/L
    # Dissolved mass = 0.80 mg/L * 0.050 L = 0.040 mg = 40 ug
    # Percent leached = 40 ug / 400 ug = 10.0%
    leach_res = calculate_catalyst_leaching(
        measured_concentration_mg_L=0.40,
        electrolyte_volume_mL=50.0,
        catalyst_mass_mg=2.0,
        initial_metal_wt_fraction=0.20,
        dilution_factor=2.0,
    )

    assert leach_res["corrected_concentration_mg_L"] == pytest.approx(0.80)
    assert leach_res["dissolved_metal_mass_ug"] == pytest.approx(40.0)
    assert leach_res["percentage_metal_leached"] == pytest.approx(10.0)
    assert leach_res["mass_normalized_leaching_ug_per_mg"] == pytest.approx(20.0)


def test_icp_plotting():
    """Verify ICP calibration and leaching comparison plots."""
    concs = [0.0, 0.5, 1.0, 5.0, 10.0]
    intensities = [50.0, 6300.0, 12600.0, 62500.0, 125000.0]
    cal_res = fit_icp_calibration(concs, intensities, blank_intensities=[48.0, 52.0, 50.0])

    fig1, axes1 = plot_icp_calibration(cal_res, element="Platinum")
    assert fig1 is not None

    fig2, ax2 = plot_icp_leaching_comparison(
        samples=["Pristine", "10h OER", "50h OER", "100h OER"],
        leached_percentages=[0.1, 2.5, 8.2, 18.5],
        element="Nickel",
    )
    assert fig2 is not None

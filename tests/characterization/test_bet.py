"""
Unit and Numerical Validation Tests for BET and BJH Gas Adsorption Module.
"""

from pathlib import Path
import numpy as np
import pytest

from parsers.characterization.common import CharacterizationData
from parsers.characterization.adsorption_isotherm import parse_isotherm_file
from modules.characterization.bet import (
    calculate_bet_surface_area,
    calculate_bjh_pore_distribution,
    plot_isotherm,
    plot_bet_linear,
    plot_bjh_pore_distribution,
)

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures" / "bet"


def test_parse_isotherm_file():
    """Verify parsing of gas adsorption isotherm with branches."""
    data = parse_isotherm_file(FIXTURES_DIR / "mesoporous_catalyst_n2_isotherm.csv")
    assert isinstance(data, CharacterizationData)
    assert data.technique == "BET"
    assert "adsorption" in data.dataframe["branch"].values
    assert "desorption" in data.dataframe["branch"].values
    assert data.metadata["adsorbate"] == "N2"


def test_bet_surface_area_and_rouquerol_criteria():
    """Verify linear BET regression, specific surface area, and Rouquerol checks."""
    data = parse_isotherm_file(FIXTURES_DIR / "mesoporous_catalyst_n2_isotherm.csv")
    df_ads = data.dataframe[data.dataframe["branch"] == "adsorption"]

    bet_res = calculate_bet_surface_area(
        p_p0=df_ads["p_p0"].values,
        quantity_adsorbed=df_ads["quantity_adsorbed"].values,
        p_range=(0.05, 0.30),
    )

    # Generated with monolayer capacity ~ 35 cm^3/g STP -> S_BET ~ 152 m^2/g
    assert bet_res["r_squared"] > 0.99
    assert bet_res["s_bet_m2_g"] == pytest.approx(152.0, abs=15.0)
    assert bet_res["c_constant"] > 50.0  # Positive and characteristic of N2 physisorption

    # Rouquerol criteria
    rouq = bet_res["rouquerol_criteria"]
    assert rouq["c_positive"] is True
    assert rouq["monotonic_increasing"] is True
    assert rouq["overall_valid"] is True


def test_bjh_pore_size_distribution():
    """Verify BJH mesopore distribution calculation on desorption branch."""
    data = parse_isotherm_file(FIXTURES_DIR / "mesoporous_catalyst_n2_isotherm.csv")
    df_des = data.dataframe[data.dataframe["branch"] == "desorption"]

    bjh_res = calculate_bjh_pore_distribution(
        p_p0=df_des["p_p0"].values,
        quantity_adsorbed=df_des["quantity_adsorbed"].values,
        branch="desorption",
    )

    assert len(bjh_res["pore_diameter_nm"]) > 5
    assert np.all(bjh_res["pore_diameter_nm"] > 1.0)
    assert bjh_res["total_pore_volume_cm3_g"] > 0.05
    # Capillary condensation occurs in mesopore range ~ 3-20 nm
    assert 2.0 < bjh_res["mean_pore_diameter_nm"] < 25.0


def test_bet_and_bjh_plotting():
    """Verify rendering of isotherm, BET linear fit, and BJH pore distribution."""
    data = parse_isotherm_file(FIXTURES_DIR / "mesoporous_catalyst_n2_isotherm.csv")
    fig1, ax1 = plot_isotherm(data, title="N2 Physisorption Isotherm (77 K)")
    assert fig1 is not None

    df_ads = data.dataframe[data.dataframe["branch"] == "adsorption"]
    bet_res = calculate_bet_surface_area(df_ads["p_p0"].values, df_ads["quantity_adsorbed"].values)
    fig2, ax2 = plot_bet_linear(bet_res)
    assert fig2 is not None

    df_des = data.dataframe[data.dataframe["branch"] == "desorption"]
    bjh_res = calculate_bjh_pore_distribution(df_des["p_p0"].values, df_des["quantity_adsorbed"].values)
    fig3, ax3 = plot_bjh_pore_distribution(bjh_res)
    assert fig3 is not None

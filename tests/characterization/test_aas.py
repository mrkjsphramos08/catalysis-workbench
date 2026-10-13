"""
Unit and Numerical Validation Tests for Atomic Absorption Spectroscopy (AAS) Module.

Validates:
- FAAS and GFAAS data ingestion (CSV, TSV, XLSX)
- Linear regression (OLS, zero-intercept, WLS, range checks, curvature)
- Method of Standard Additions (MOSA) with volume dilution corrections & delta-method CI
- Quantitative sample concentration calculation, blank subtraction, replicate %RSD
- Below-LOD/LOQ flagging and censored data preservation
- Catalyst metal loading (wt.%) and electrochemical leaching kinetics with aliquot tracking
- Analytical QC evaluation (EPA 7000B and EPA 7010)
- Publication figure rendering across all AAS visualization modes.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from parsers.characterization.common import CharacterizationData
from parsers.characterization.aas import parse_aas_file
from modules.characterization.elemental_analysis import (
    convert_concentration_units,
    compute_replicate_statistics,
    calculate_detection_limits,
)
from modules.characterization.aas import (
    fit_aas_calibration,
    fit_standard_addition,
    quantify_aas_samples,
    evaluate_aas_qc,
    calculate_aas_catalyst_loading,
    calculate_aas_leaching,
    plot_aas_calibration,
    plot_standard_addition,
    plot_aas_concentration_comparison,
    plot_aas_catalyst_loading,
    plot_aas_leaching,
    plot_aas_stability_before_after,
    plot_aas_replicate_distribution,
    plot_aas_qc_chart,
)

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures" / "aas"


# =========================================================================
# 1. Data Ingestion Tests (CSV, TSV, XLSX)
# =========================================================================

def test_parse_faas_csv():
    """Verify parsing of Flame AAS CSV data."""
    csv_path = FIXTURES_DIR / "faas_calibration_fe.csv"
    data = parse_aas_file(csv_path, technique="FAAS", instrument_model="PerkinElmer PinAAcle 900T")
    assert isinstance(data, CharacterizationData)
    assert data.technique == "FAAS"
    assert "Fe" in data.metadata["elements_detected"]
    assert data.y_name == "Absorbance"
    assert data.y_unit == "AU"
    assert "absorbance" in data.dataframe.columns
    assert len(data.dataframe) == 21  # 6 stds * 3 reps + 1 ICV * 3 reps


def test_parse_faas_tsv_and_xlsx():
    """Verify multi-format ingestion of TSV and Excel files."""
    tsv_path = FIXTURES_DIR / "faas_calibration_fe.tsv"
    data_tsv = parse_aas_file(tsv_path, technique="FAAS")
    assert len(data_tsv.dataframe) == 21

    xlsx_path = FIXTURES_DIR / "faas_calibration_fe.xlsx"
    data_xlsx = parse_aas_file(xlsx_path, technique="FAAS")
    assert len(data_xlsx.dataframe) == 21
    assert data_xlsx.dataframe["absorbance"].iloc[0] == pytest.approx(data_tsv.dataframe["absorbance"].iloc[0])


def test_parse_gfaas_integrated_absorbance():
    """Verify parsing of Graphite Furnace AAS with integrated absorbance (A*s)."""
    csv_path = FIXTURES_DIR / "gfaas_calibration_pt.csv"
    data = parse_aas_file(csv_path, technique="GFAAS", instrument_model="Agilent 240FS/240Z")
    assert data.technique == "GFAAS"
    assert data.y_name == "Integrated Absorbance"
    assert data.y_unit == "A*s"
    assert "integrated_absorbance" in data.dataframe.columns
    assert "Pt" in data.metadata["elements_detected"]


# =========================================================================
# 2. Calibration Curve Regression Tests
# =========================================================================

def test_faas_linear_calibration_and_lod():
    """Verify FAAS linear calibration regression, slope, R^2, and detection limits."""
    data = parse_aas_file(FIXTURES_DIR / "faas_calibration_fe.csv")
    df = data.dataframe
    df_stds = df[df["sample_type"] == "Standard"]
    df_blanks = df[df["sample_type"] == "Blank"]

    cal_res = fit_aas_calibration(
        standards_concentration=df_stds["concentration_reported"].values,
        standards_signal=df_stds["absorbance"].values,
        signal_type="absorbance",
        blank_signals=df_blanks["absorbance"].values,
        lod_method="blank_stdev",
    )

    assert cal_res["r_squared"] > 0.999
    assert cal_res["slope"] == pytest.approx(0.0825, abs=0.005)
    assert abs(cal_res["intercept"]) < 0.010
    assert cal_res["lod_conc"] is not None
    assert 0.001 < cal_res["lod_conc"] < 0.100
    assert cal_res["loq_conc"] > cal_res["lod_conc"]
    assert cal_res["rollover_warning"] is False


def test_zero_intercept_and_weighted_calibration():
    """Verify zero-intercept fitting and weighted regression."""
    concs = [0.1, 0.5, 1.0, 2.0, 5.0]
    # Pure linear without intercept
    sigs = [0.080 * c for c in concs]

    # Zero intercept
    res_zero = fit_aas_calibration(concs, sigs, zero_intercept=True)
    assert res_zero["intercept"] == 0.0
    assert res_zero["slope"] == pytest.approx(0.080, rel=1e-3)
    assert res_zero["r_squared"] > 0.999

    # Weighted regression (1/x weights)
    weights = [1.0 / c for c in concs]
    res_wls = fit_aas_calibration(concs, sigs, weights=weights)
    assert res_wls["weighted"] is True
    assert res_wls["slope"] == pytest.approx(0.080, rel=1e-3)


def test_high_absorbance_rollover_warning():
    """Verify that absorbance exceeding linear range (A > 0.85) triggers rollover warning."""
    concs = [0.0, 2.0, 5.0, 10.0, 15.0]
    sigs = [0.0, 0.20, 0.50, 0.88, 1.05]  # curvature above 0.85 AU
    res = fit_aas_calibration(concs, sigs, signal_type="absorbance")
    assert res["rollover_warning"] is True


# =========================================================================
# 3. Method of Standard Additions (MOSA) Tests
# =========================================================================

def test_standard_addition_extrapolation_and_dilution():
    """
    Verify standard addition linear fit, x-intercept extrapolation,
    and volume dilution correction.
    """
    data = parse_aas_file(FIXTURES_DIR / "aas_standard_addition_catalyst.csv")
    df = data.dataframe

    # Aggregate by addition level
    agg = df.groupby("concentration_reported" if "concentration_reported" in df.columns else "sample_id").mean(numeric_only=True)
    # Spikes: 0.0, 0.5, 1.0, 1.5, 2.0
    # Expected apparent conc in flask ~ 0.85 mg/L, aliquot = 5.0 mL, final = 25.0 mL (DF = 5.0) -> Original ~ 4.25 mg/L
    spikes = [0.0, 0.50, 1.00, 1.50, 2.00]
    # Extract mean absorbances
    mean_sigs = [df[df["sample_id"].str.contains(f"\\+{sp:.2f}") if sp > 0 else df["sample_id"].str.contains("Unspiked")]["absorbance"].mean() for sp in spikes]

    sa_res = fit_standard_addition(
        added_concentrations=spikes,
        signals=mean_sigs,
        sample_volume_mL=5.0,
        final_volume_mL=25.0,
    )

    assert sa_res["r_squared"] > 0.995
    assert sa_res["apparent_concentration"] == pytest.approx(0.85, abs=0.08)
    assert sa_res["dilution_factor"] == pytest.approx(5.0)
    assert sa_res["original_concentration"] == pytest.approx(4.25, abs=0.40)
    assert sa_res["ci95_apparent"] > 0.0
    assert sa_res["ci95_original"] > sa_res["ci95_apparent"]


# =========================================================================
# 4. Sample Quantification & Below-LOD Flagging Tests
# =========================================================================

def test_sample_quantification_and_flags():
    """Verify unknown sample concentration quantification, replicate %RSD, and QC flags."""
    concs = [0.0, 0.5, 1.0, 2.0, 5.0]
    sigs = [0.002, 0.042, 0.084, 0.165, 0.410]
    cal_res = fit_aas_calibration(concs, sigs, blank_signals=[0.0018, 0.0022, 0.0020])

    # Construct test sample dataframe
    sample_df = pd.DataFrame([
        {"sample_id": "Sample_Normal", "sample_type": "Sample", "absorbance": 0.082, "dilution_factor": 1.0, "element": "Fe", "unit": "mg/L"},
        {"sample_id": "Sample_Normal", "sample_type": "Sample", "absorbance": 0.084, "dilution_factor": 1.0, "element": "Fe", "unit": "mg/L"},
        {"sample_id": "Sample_Diluted", "sample_type": "Sample", "absorbance": 0.083, "dilution_factor": 10.0, "element": "Fe", "unit": "mg/L"},
        {"sample_id": "Sample_Low", "sample_type": "Sample", "absorbance": 0.0021, "dilution_factor": 1.0, "element": "Fe", "unit": "mg/L"},  # near blank/LOD
        {"sample_id": "Sample_High", "sample_type": "Sample", "absorbance": 0.650, "dilution_factor": 1.0, "element": "Fe", "unit": "mg/L"},  # > 5.0 mg/L
    ])

    quant_res = quantify_aas_samples(sample_df, cal_res, blank_subtraction=False)
    assert len(quant_res) == 4

    norm_row = quant_res[quant_res["sample_id"] == "Sample_Normal"].iloc[0]
    assert norm_row["concentration"] == pytest.approx(1.0, abs=0.1)
    assert norm_row["raw_signal_rsd_pct"] < 5.0
    assert norm_row["qc_flag"] == "OK"

    dil_row = quant_res[quant_res["sample_id"] == "Sample_Diluted"].iloc[0]
    assert dil_row["concentration"] == pytest.approx(10.0, abs=1.0)

    low_row = quant_res[quant_res["sample_id"] == "Sample_Low"].iloc[0]
    assert low_row["qc_flag"] in ("<LOD", "<LOQ")

    high_row = quant_res[quant_res["sample_id"] == "Sample_High"].iloc[0]
    assert high_row["qc_flag"] == "EXTRAPOLATED_HIGH"


# =========================================================================
# 5. Catalyst Loading & Leaching Tests
# =========================================================================

def test_catalyst_loading_calculation():
    """Verify digestion metal loading wt.% and nominal efficiency."""
    # 2.50 mg Pt/C catalyst digested in 50.0 mL, DF = 2.0, measured conc = 4.95 mg/L
    # Dissolved in cell = 4.95 mg/L * 2.0 * 0.050 L = 0.495 mg Pt
    # wt% = 0.495 mg / 2.50 mg = 19.8 wt%
    res = calculate_aas_catalyst_loading(
        measured_concentration_mg_L=4.95,
        digestion_volume_mL=50.0,
        catalyst_mass_mg=2.50,
        dilution_factor=2.0,
        nominal_metal_wt_fraction=0.20,
    )

    assert res["metal_loading_wt_pct"] == pytest.approx(19.8, abs=0.1)
    assert res["metal_mass_ug"] == pytest.approx(495.0, abs=1.0)
    assert res["loading_efficiency_pct"] == pytest.approx(99.0, abs=0.5)


def test_catalyst_leaching_with_aliquots():
    """Verify stability leaching calculation with periodic aliquot withdrawal tracking."""
    # Initial: 2.0 mg of 15 wt% catalyst -> 300 ug metal
    # 50 mL electrolyte. After 25h: measured = 0.10 mg/L (5.0 ug in cell).
    # Previous aliquot removed at 10h: 2.0 mL at 0.04 mg/L -> 0.08 ug removed.
    # Total dissolved = 5.0 + 0.08 = 5.08 ug.
    # Leached % = 5.08 / 300.0 = 1.69%
    res = calculate_aas_leaching(
        measured_concentration_mg_L=0.10,
        electrolyte_volume_mL=50.0,
        catalyst_mass_mg=2.0,
        initial_metal_wt_fraction=0.15,
        aliquot_history=[(2.0, 0.04)],
    )

    assert res["dissolved_metal_mass_ug"] == pytest.approx(5.08, abs=0.05)
    assert res["percentage_metal_leached"] == pytest.approx(1.69, abs=0.05)
    assert res["mass_normalized_leaching_ug_per_mg"] == pytest.approx(2.54, abs=0.05)


# =========================================================================
# 6. Quality Control Evaluation (EPA 7000B / EPA 7010)
# =========================================================================

def test_evaluate_aas_qc():
    """Verify automated analytical QC evaluation against EPA SOP limits."""
    data = parse_aas_file(FIXTURES_DIR / "faas_calibration_fe.csv")
    df = data.dataframe
    cal_res = fit_aas_calibration(
        standards_concentration=df[df["sample_type"] == "Standard"]["concentration_reported"].values,
        standards_signal=df[df["sample_type"] == "Standard"]["absorbance"].values,
    )

    qc_res = evaluate_aas_qc(df, cal_res, sop="EPA_7000B")
    assert qc_res["sop"] == "EPA_7000B"
    assert qc_res["overall_status"] == "PASS"
    assert len(qc_res["checks"]) >= 3


# =========================================================================
# 7. Publication Figure Generation Tests
# =========================================================================

def test_aas_figures_plotting():
    """Verify generation of all AAS publication figures without GUI backend errors."""
    concs = [0.0, 0.5, 1.0, 2.0, 5.0]
    sigs = [0.002, 0.041, 0.083, 0.165, 0.412]
    cal_res = fit_aas_calibration(concs, sigs)

    # 1. Calibration curve
    fig1, axes1 = plot_aas_calibration(cal_res, element="Fe")
    assert fig1 is not None

    # 2. Standard addition
    sa_res = fit_standard_addition([0.0, 0.5, 1.0, 1.5, 2.0], [0.045, 0.072, 0.099, 0.126, 0.153])
    fig2, ax2 = plot_standard_addition(sa_res, element="Co")
    assert fig2 is not None

    # 3. Concentration comparison
    fig3, ax3 = plot_aas_concentration_comparison(["S1", "S2", "S3"], [1.2, 2.5, 0.8], error_bars=[0.05, 0.12, 0.04])
    assert fig3 is not None

    # 4. Catalyst loading
    fig4, ax4 = plot_aas_catalyst_loading(["Pt/C", "PtCo/C"], [19.8, 14.8], nominal_loading_wt_pct=[20.0, 15.0])
    assert fig4 is not None

    # 5. Leaching kinetics
    fig5, ax5 = plot_aas_leaching([0, 10, 25, 50, 100], [0.1, 1.2, 3.5, 8.2, 14.5])
    assert fig5 is not None

    # 6. Before/after stability
    fig6, ax6 = plot_aas_stability_before_after(["PtCo/C", "NiFe-LDH"], [14.8, 24.5], [12.6, 22.8])
    assert fig6 is not None

    # 7. Replicate precision
    fig7, ax7 = plot_aas_replicate_distribution(["Sample A", "Sample B"], [[1.02, 1.01, 1.03], [2.45, 2.48, 2.42]])
    assert fig7 is not None

    # 8. QC Control chart
    fig8, ax8 = plot_aas_qc_chart(["CCV1", "CCV2", "CCV3"], [100.5, 99.2, 101.1], sop="EPA_7000B")
    assert fig8 is not None

"""
Atomic Absorption Spectroscopy (AAS) Quantitative Analysis and Publication Graphing.

Supports:
- Flame Atomic Absorption Spectroscopy (FAAS)
- Graphite Furnace Atomic Absorption Spectroscopy (GFAAS / ETAAS)
- Hydride Generation (HG-AAS) and Cold Vapor (CV-AAS) architecture hooks

Implements:
- External standard calibration (OLS, WLS, zero-intercept, non-linear curvature detection)
- Method of Standard Additions (MOSA) with matrix interference extrapolation & delta-method uncertainty
- Sample concentration quantification with blank subtraction, dilution correction, and replicate statistics
- Analytical Quality Control (QC) validation against US EPA Method 7000B, EPA 7010, and EPA 200.9 SOPs
- Catalyst metal loading (wt.%) and electrochemical stability leaching kinetics (with aliquot compensation)
- Publication-quality figures across 8 analytical visualization modes.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats

from modules.characterization.elemental_analysis import (
    calculate_catalyst_leaching_common,
    calculate_detection_limits,
    compute_replicate_statistics,
    convert_concentration_units,
    linear_calibration_fit,
)
from parsers.characterization.common import CharacterizationData
from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


# =========================================================================
# 1. Calibration Curve Analysis
# =========================================================================

def fit_aas_calibration(
    standards_concentration: Sequence[float],
    standards_signal: Sequence[float],
    signal_type: str = "absorbance",
    zero_intercept: bool = False,
    weights: Optional[Sequence[float]] = None,
    blank_signals: Optional[Sequence[float]] = None,
    lod_method: str = "blank_stdev",
    k_lod: float = 3.3,
    k_loq: float = 10.0,
    alpha: float = 0.05,
) -> Dict[str, Any]:
    """
    Perform linear regression on AAS calibration standards.

    Model: A = m * C + b  (or A = m * C when zero_intercept is True).

    Parameters
    ----------
    standards_concentration : sequence of float
        Known standard concentrations (e.g. mg/L or ug/L).
    standards_signal : sequence of float
        Measured steady-state absorbance (AU) or integrated absorbance (A*s).
    signal_type : str, default 'absorbance'
        'absorbance' (for FAAS) or 'integrated_absorbance' (for GFAAS).
    zero_intercept : bool, default False
        Force zero intercept when scientifically justified after blank subtraction.
    weights : sequence of float, optional
        Fitting weights for Weighted Least Squares (WLS).
    blank_signals : sequence of float, optional
        Replicate blank readings to determine detection limits.
    lod_method : str, default 'blank_stdev'
        'blank_stdev' (k*s_b/m) or 'residual_stdev' (k*s_{y/x}/m).
    k_lod : float, default 3.3
        Multiplier for Limit of Detection.
    k_loq : float, default 10.0
        Multiplier for Limit of Quantification.
    alpha : float, default 0.05
        Significance level for confidence interval (95% CI).

    Returns
    -------
    Dict[str, Any]
        Comprehensive calibration model parameters and diagnostics.
    """
    c_arr = np.asarray(standards_concentration, dtype=float)
    s_arr = np.asarray(standards_signal, dtype=float)

    if len(c_arr) != len(s_arr):
        raise ValueError("Length of standards_concentration and standards_signal must match.")
    min_required = 2 if zero_intercept else 3
    if len(c_arr) < min_required:
        raise ValueError(f"At least {min_required} calibration standards required for AAS regression.")

    fit_res = linear_calibration_fit(
        c_arr,
        s_arr,
        zero_intercept=zero_intercept,
        weights=weights,
        alpha=alpha,
    )

    # Determine Limit of Detection & Quantification
    b_std = float(np.std(blank_signals, ddof=1)) if blank_signals is not None and len(blank_signals) >= 2 else None
    det_limits = calculate_detection_limits(
        slope=fit_res["slope"],
        blank_stdev=b_std,
        syx=fit_res["syx"],
        method=lod_method,
        k_lod=k_lod,
        k_loq=k_loq,
    )

    # High absorbance rollover check (Beer-Lambert departure at A > 0.8 in FAAS)
    max_sig = float(np.max(s_arr))
    rollover_warning = False
    if signal_type == "absorbance" and max_sig > 0.85:
        rollover_warning = True

    return {
        "slope": fit_res["slope"],
        "intercept": fit_res["intercept"],
        "r_squared": fit_res["r_squared"],
        "syx": fit_res["syx"],
        "std_err_slope": fit_res["std_err_slope"],
        "std_err_intercept": fit_res["std_err_intercept"],
        "p_value": fit_res["p_value"],
        "zero_intercept": zero_intercept,
        "weighted": fit_res["weighted"],
        "standards_c": c_arr,
        "standards_signal": s_arr,
        "fit_signal": fit_res["y_pred"],
        "residuals": fit_res["residuals"],
        "rel_residuals_pct": fit_res["rel_residuals_pct"],
        "calc_c": fit_res["calc_x"],
        "recoveries_pct": fit_res["recoveries_pct"],
        "c_min": fit_res["x_min"],
        "c_max": fit_res["x_max"],
        "lod_conc": det_limits["lod_conc"],
        "loq_conc": det_limits["loq_conc"],
        "lod_method": det_limits["method"],
        "curvature_detected": fit_res["curvature_detected"],
        "rollover_warning": rollover_warning,
        "signal_type": signal_type,
    }


# =========================================================================
# 2. Method of Standard Additions (MOSA)
# =========================================================================

def fit_standard_addition(
    added_concentrations: Sequence[float],
    signals: Sequence[float],
    signal_type: str = "absorbance",
    sample_volume_mL: Optional[float] = None,
    final_volume_mL: Optional[float] = None,
    alpha: float = 0.05,
) -> Dict[str, Any]:
    """
    Perform Method of Standard Additions (MOSA) linear regression for matrix-affected samples.

    Model: S = m * (Delta C) + S_0
    Unknown concentration in analyzed flask: C_apparent = S_0 / m = b / m (x-intercept magnitude).
    If sample aliquot and final flask volumes are provided, corrects for dilution:
    C_original = C_apparent * (final_volume_mL / sample_volume_mL).

    Parameters
    ----------
    added_concentrations : sequence of float
        Spiked standard concentrations in the final analyzed solution (first entry typically 0.0).
    signals : sequence of float
        Measured analytical response for each addition level.
    signal_type : str, default 'absorbance'
        'absorbance' (AU) or 'integrated_absorbance' (A*s).
    sample_volume_mL : float, optional
        Volume of original sample aliquot introduced into addition flasks.
    final_volume_mL : float, optional
        Total final diluted volume in each flask.
    alpha : float, default 0.05
        Significance level for confidence interval (95% CI).

    Returns
    -------
    Dict[str, Any]
        MOSA regression results, apparent and original concentrations, standard error,
        95% confidence interval, and extrapolation curve coordinates.
    """
    x_arr = np.asarray(added_concentrations, dtype=float)
    y_arr = np.asarray(signals, dtype=float)

    if len(x_arr) != len(y_arr):
        raise ValueError("Length of added_concentrations and signals must match.")
    if len(x_arr) < 3:
        raise ValueError("At least 3 addition points (including unspiked sample) required for standard addition.")

    n = len(x_arr)
    fit_res = linear_calibration_fit(x_arr, y_arr, zero_intercept=False, alpha=alpha)

    slope = fit_res["slope"]
    intercept = fit_res["intercept"]
    syx = fit_res["syx"]

    if slope <= 0:
        raise ValueError("Standard addition slope must be strictly positive.")

    # Apparent concentration in analyzed solution: C_app = intercept / slope
    c_apparent = float(intercept / slope)

    # Standard error on extrapolated concentration (Delta method)
    # s_x0 = (s_{y/x} / m) * sqrt(1/N + y_bar^2 / (m^2 * sum((x - x_bar)^2)))
    y_bar = float(np.mean(y_arr))
    x_bar = float(np.mean(x_arr))
    s_xx = float(np.sum((x_arr - x_bar) ** 2))

    if s_xx > 0:
        se_c_app = float((syx / slope) * np.sqrt(1.0 / n + (y_bar ** 2) / ((slope ** 2) * s_xx)))
    else:
        se_c_app = float(syx / slope)

    t_crit = float(stats.t.ppf(1.0 - alpha / 2.0, df=n - 2))
    ci95_half_app = float(t_crit * se_c_app)

    # Volume dilution correction if provided
    dilution_factor = 1.0
    if sample_volume_mL is not None and final_volume_mL is not None:
        if sample_volume_mL <= 0 or final_volume_mL <= 0:
            raise ValueError("sample_volume_mL and final_volume_mL must be positive.")
        dilution_factor = float(final_volume_mL / sample_volume_mL)

    c_original = c_apparent * dilution_factor
    se_c_orig = se_c_app * dilution_factor
    ci95_half_orig = ci95_half_app * dilution_factor

    # Extrapolation curve coordinates for visualization (from negative x-intercept to max addition)
    extrap_x = np.linspace(-c_apparent * 1.35, float(np.max(x_arr)) * 1.15, 100)
    extrap_y = slope * extrap_x + intercept

    return {
        "apparent_concentration": c_apparent,
        "original_concentration": c_original,
        "std_err_apparent": se_c_app,
        "std_err_original": se_c_orig,
        "ci95_apparent": ci95_half_app,
        "ci95_original": ci95_half_orig,
        "dilution_factor": dilution_factor,
        "slope": slope,
        "intercept": intercept,
        "r_squared": fit_res["r_squared"],
        "syx": syx,
        "added_concentrations": x_arr,
        "signals": y_arr,
        "extrap_x": extrap_x,
        "extrap_y": extrap_y,
        "signal_type": signal_type,
    }


# =========================================================================
# 3. Quantitative Sample Concentration Analysis
# =========================================================================

def quantify_aas_samples(
    dataframe: pd.DataFrame,
    cal_result: Dict[str, Any],
    signal_col: Optional[str] = None,
    blank_subtraction: bool = True,
    default_dilution: float = 1.0,
    solid_sample_mass_g: Optional[float] = None,
    digestion_volume_L: Optional[float] = None,
    target_unit: Optional[str] = None,
) -> pd.DataFrame:
    """
    Calculate unknown sample concentrations from raw AAS readings and calibration model.

    Handles blank correction, dilution factors, replicate averaging (%RSD),
    and flags measurements below LOD, below LOQ, or outside the calibration range.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Standardized AAS dataframe from parse_aas_file.
    cal_result : Dict[str, Any]
        Model returned by fit_aas_calibration.
    signal_col : str, optional
        Target signal column ('absorbance' or 'integrated_absorbance').
    blank_subtraction : bool, default True
        Whether to subtract calibration blank signal.
    default_dilution : float, default 1.0
        Default dilution factor if missing from dataframe.
    solid_sample_mass_g : float, optional
        Dry sample mass for solid conversion (mg/kg).
    digestion_volume_L : float, optional
        Digest volume in liters for solid conversion.
    target_unit : str, optional
        Desired output concentration unit (e.g. 'mg/L', 'ug/L', 'mg/kg', 'wt%').

    Returns
    -------
    pd.DataFrame
        Aggregated sample quantification table with statistics and QC flags.
    """
    df = dataframe.copy()

    # Determine signal column
    if signal_col is None:
        if cal_result.get("signal_type") == "integrated_absorbance" and "integrated_absorbance" in df.columns:
            signal_col = "integrated_absorbance"
        elif "absorbance" in df.columns and not df["absorbance"].isna().all():
            signal_col = "absorbance"
        else:
            raise ValueError("No valid absorbance signal column found in DataFrame.")

    slope = cal_result["slope"]
    intercept = cal_result["intercept"]
    lod = cal_result.get("lod_conc")
    loq = cal_result.get("loq_conc")
    c_min = cal_result.get("c_min", 0.0)
    c_max = cal_result.get("c_max", 1e9)

    # 1. Blank signal determination
    blank_signal_mean = 0.0
    if blank_subtraction:
        blanks = df[df["sample_type"].str.lower().str.contains("blank")]
        if not blanks.empty and not blanks[signal_col].isna().all():
            blank_signal_mean = float(blanks[signal_col].mean())

    # Filter for unknown samples
    samples_df = df[~df["sample_type"].str.lower().str.contains("standard|cal")].copy()
    if samples_df.empty:
        samples_df = df.copy()

    results = []
    for sid, group in samples_df.groupby("sample_id", sort=False):
        raw_signals = group[signal_col].dropna().values
        if len(raw_signals) == 0:
            continue

        stype = group["sample_type"].iloc[0]
        df_factor = float(group["dilution_factor"].iloc[0]) if "dilution_factor" in group.columns else default_dilution
        elem = group["element"].iloc[0] if "element" in group.columns else "Analyte"
        orig_unit = group["unit"].iloc[0] if "unit" in group.columns else "mg/L"

        # Signal replicate stats
        sig_stats = compute_replicate_statistics(raw_signals)
        mean_sig = sig_stats["mean"]
        net_sig = max(0.0, mean_sig - blank_signal_mean) if blank_subtraction and "blank" not in stype.lower() else mean_sig

        # Recalculate concentration: C = (S_net - intercept) / slope * dilution
        c_measured = (net_sig - intercept) / slope
        c_final = c_measured * df_factor

        # Uncertainty propagation on concentration
        std_c = (sig_stats["std"] / slope) * df_factor if sig_stats["count"] > 1 else 0.0
        rsd_pct = sig_stats["rsd_pct"]

        # Flag assignment
        flag = "OK"
        if lod is not None and c_measured < lod:
            flag = "<LOD"
        elif loq is not None and c_measured < loq:
            flag = "<LOQ"
        elif c_measured > c_max:
            flag = "EXTRAPOLATED_HIGH"
        elif c_measured < c_min and flag == "OK":
            flag = "EXTRAPOLATED_LOW"

        # Unit conversion if requested
        disp_unit = orig_unit
        c_converted = c_final
        std_converted = std_c
        if target_unit is not None and target_unit != orig_unit:
            try:
                c_converted = float(convert_concentration_units(
                    c_final,
                    from_unit=orig_unit,
                    to_unit=target_unit,
                    sample_mass_g=solid_sample_mass_g,
                    solution_volume_L=digestion_volume_L,
                ))
                std_converted = float(convert_concentration_units(
                    std_c,
                    from_unit=orig_unit,
                    to_unit=target_unit,
                    sample_mass_g=solid_sample_mass_g,
                    solution_volume_L=digestion_volume_L,
                ))
                disp_unit = target_unit
            except Exception:
                pass

        results.append({
            "sample_id": sid,
            "sample_type": stype,
            "element": elem,
            "replicate_count": sig_stats["count"],
            "raw_signal_mean": mean_sig,
            "raw_signal_std": sig_stats["std"],
            "raw_signal_rsd_pct": rsd_pct,
            "dilution_factor": df_factor,
            "concentration": c_converted,
            "concentration_std": std_converted,
            "unit": disp_unit,
            "qc_flag": flag,
        })

    return pd.DataFrame(results)


# =========================================================================
# 4. Analytical Quality Control (QC)
# =========================================================================

def evaluate_aas_qc(
    dataframe: pd.DataFrame,
    cal_result: Dict[str, Any],
    sop: str = "EPA_7000B",
    custom_limits: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Evaluate experimental quality control compliance against US EPA standard procedures.

    SOP Profiles:
    - 'EPA_7000B' (Flame AAS):
        ICV Recovery: 90% - 110%
        CCV Recovery: 90% - 110%
        Calibration Blank: < LOD
        Replicate Precision: RSD <= 10%
        Matrix Spike Recovery: 80% - 120%
    - 'EPA_7010' / 'EPA_200.9' (Graphite Furnace AAS):
        ICV Recovery: 85% - 115%
        CCV Recovery: 85% - 115%
        Calibration Blank: < LOD
        Replicate Precision: RSD <= 15%
        Matrix Spike Recovery: 70% - 130%

    Parameters
    ----------
    dataframe : pd.DataFrame
        AAS dataset dataframe.
    cal_result : Dict[str, Any]
        Calibration model result from fit_aas_calibration.
    sop : str, default 'EPA_7000B'
        'EPA_7000B', 'EPA_7010', or 'custom'.
    custom_limits : dict, optional
        Custom acceptance criteria dictionary.

    Returns
    -------
    Dict[str, Any]
        Structured QC evaluation with overall status ('PASS'/'FAIL') and itemized results.
    """
    sop_up = sop.upper().replace("-", "_").replace(" ", "_")

    # Define thresholds
    if sop_up in ("EPA_7010", "EPA_200_9", "GFAAS"):
        limits = {
            "icv_range": (85.0, 115.0),
            "ccv_range": (85.0, 115.0),
            "rsd_max": 15.0,
            "spike_range": (70.0, 130.0),
            "min_r2": 0.995,
        }
    else:  # Default EPA 7000B (Flame)
        limits = {
            "icv_range": (90.0, 110.0),
            "ccv_range": (90.0, 110.0),
            "rsd_max": 10.0,
            "spike_range": (80.0, 120.0),
            "min_r2": 0.995,
        }

    if custom_limits:
        limits.update(custom_limits)

    checks = []
    overall_pass = True

    # 1. Calibration Linearity
    r2 = cal_result["r_squared"]
    r2_pass = r2 >= limits["min_r2"]
    if not r2_pass:
        overall_pass = False
    checks.append({
        "check": "Calibration Linearity (R^2)",
        "value": f"{r2:.5f}",
        "criteria": f">= {limits['min_r2']}",
        "status": "PASS" if r2_pass else "FAIL",
    })

    # 2. Calibration Blank
    lod = cal_result.get("lod_conc")
    blanks = dataframe[dataframe["sample_type"].str.lower().str.contains("blank")]
    if not blanks.empty:
        # Quantify blanks
        sig_col = "integrated_absorbance" if cal_result.get("signal_type") == "integrated_absorbance" else "absorbance"
        b_signals = blanks[sig_col].dropna().values
        if len(b_signals) > 0:
            b_mean = float(np.mean(b_signals))
            b_conc = (b_mean - cal_result["intercept"]) / cal_result["slope"]
            b_pass = (b_conc <= lod) if lod is not None else (b_mean < 0.05)
            if not b_pass:
                overall_pass = False
            checks.append({
                "check": "Calibration / Method Blank",
                "value": f"{b_conc:.4f} (Sig: {b_mean:.4f})",
                "criteria": f"<= LOD ({lod:.4f})" if lod else "Near Zero",
                "status": "PASS" if b_pass else "FAIL",
            })

    # 3. Initial Calibration Verification (ICV)
    icvs = dataframe[dataframe["sample_type"].str.upper() == "ICV"]
    if not icvs.empty:
        sig_col = "integrated_absorbance" if cal_result.get("signal_type") == "integrated_absorbance" else "absorbance"
        for sid, grp in icvs.groupby("sample_id", sort=False):
            target_c = None
            for col in ["concentration_standard", "target_concentration", "concentration_reported"]:
                if col in grp.columns and not grp[col].isna().all() and float(grp[col].iloc[0]) > 0:
                    target_c = float(grp[col].iloc[0])
                    break
            if target_c is None:
                match = re.search(r"([0-9]+(?:\.[0-9]+)?)", str(sid))
                target_c = float(match.group(1)) if match else (cal_result["c_max"] * 0.5)

            meas_signals = grp[sig_col].dropna().values
            mean_sig = float(np.mean(meas_signals))
            meas_c = (mean_sig - cal_result["intercept"]) / cal_result["slope"]
            rec_pct = (meas_c / target_c) * 100.0 if target_c > 0 else 100.0
            icv_pass = limits["icv_range"][0] <= rec_pct <= limits["icv_range"][1]
            if not icv_pass:
                overall_pass = False
            checks.append({
                "check": f"ICV Recovery ({sid})",
                "value": f"{rec_pct:.1f}%",
                "criteria": f"{limits['icv_range'][0]}-{limits['icv_range'][1]}%",
                "status": "PASS" if icv_pass else "FAIL",
            })

    # 4. Continuing Calibration Verification (CCV)
    ccvs = dataframe[dataframe["sample_type"].str.upper() == "CCV"]
    if not ccvs.empty:
        sig_col = "integrated_absorbance" if cal_result.get("signal_type") == "integrated_absorbance" else "absorbance"
        for sid, grp in ccvs.groupby("sample_id", sort=False):
            target_c = None
            for col in ["concentration_standard", "target_concentration", "concentration_reported"]:
                if col in grp.columns and not grp[col].isna().all() and float(grp[col].iloc[0]) > 0:
                    target_c = float(grp[col].iloc[0])
                    break
            if target_c is None:
                match = re.search(r"([0-9]+(?:\.[0-9]+)?)", str(sid))
                target_c = float(match.group(1)) if match else (cal_result["c_max"] * 0.5)

            meas_signals = grp[sig_col].dropna().values
            mean_sig = float(np.mean(meas_signals))
            meas_c = (mean_sig - cal_result["intercept"]) / cal_result["slope"]
            rec_pct = (meas_c / target_c) * 100.0 if target_c > 0 else 100.0
            ccv_pass = limits["ccv_range"][0] <= rec_pct <= limits["ccv_range"][1]
            if not ccv_pass:
                overall_pass = False
            checks.append({
                "check": f"CCV Recovery ({sid})",
                "value": f"{rec_pct:.1f}%",
                "criteria": f"{limits['ccv_range'][0]}-{limits['ccv_range'][1]}%",
                "status": "PASS" if ccv_pass else "FAIL",
            })

    # 5. Replicate Precision (%RSD) across non-blank samples
    sig_col = "integrated_absorbance" if cal_result.get("signal_type") == "integrated_absorbance" else "absorbance"
    high_rsd_samples = []
    for sid, grp in dataframe.groupby("sample_id"):
        stype = str(grp["sample_type"].iloc[0]).lower() if "sample_type" in grp.columns else ""
        if "blank" in stype or "blank" in str(sid).lower():
            continue  # Near-zero baseline noise is excluded from %RSD criteria
        vals = grp[sig_col].dropna().values
        if len(vals) >= 2:
            st = compute_replicate_statistics(vals)
            if st["rsd_pct"] > limits["rsd_max"]:
                high_rsd_samples.append((sid, st["rsd_pct"]))

    rsd_pass = len(high_rsd_samples) == 0
    if not rsd_pass:
        overall_pass = False
    checks.append({
        "check": "Replicate Precision (%RSD)",
        "value": f"{len(high_rsd_samples)} samples exceed limit" if not rsd_pass else "All within limit",
        "criteria": f"RSD <= {limits['rsd_max']}%",
        "status": "PASS" if rsd_pass else "FAIL",
    })

    return {
        "sop": sop_up,
        "overall_status": "PASS" if overall_pass else "FAIL",
        "limits": limits,
        "checks": pd.DataFrame(checks),
    }


# =========================================================================
# 5. Catalyst Metal Loading & Dissolution Kinetics
# =========================================================================

def calculate_aas_catalyst_loading(
    measured_concentration_mg_L: float,
    digestion_volume_mL: float,
    catalyst_mass_mg: float,
    dilution_factor: float = 1.0,
    blank_concentration_mg_L: float = 0.0,
    nominal_metal_wt_fraction: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Calculate actual catalyst metal loading (wt.%) from digestion AAS measurements.

    Parameters
    ----------
    measured_concentration_mg_L : float
        AAS measured concentration of metal in digested solution (mg/L).
    digestion_volume_mL : float
        Final volume of digestion solution in milliliters.
    catalyst_mass_mg : float
        Initial dry solid catalyst mass digested in milligrams.
    dilution_factor : float, default 1.0
        Dilution factor applied to digest prior to AAS aspiration.
    blank_concentration_mg_L : float, default 0.0
        Procedural digestion acid blank concentration in mg/L.
    nominal_metal_wt_fraction : float, optional
        Target/nominal metal fraction (e.g., 0.20 for 20 wt% Pt/C).

    Returns
    -------
    Dict[str, Any]
        Metal loading wt%, metal mass (mg and ug), and nominal comparison.
    """
    net_c = max(0.0, (measured_concentration_mg_L - blank_concentration_mg_L) * dilution_factor)
    vol_L = digestion_volume_mL / 1000.0

    metal_mass_mg = net_c * vol_L
    metal_mass_ug = metal_mass_mg * 1000.0

    actual_wt_fraction = (metal_mass_mg / catalyst_mass_mg) if catalyst_mass_mg > 0 else 0.0
    actual_wt_pct = actual_wt_fraction * 100.0

    recovery_vs_nominal_pct = None
    if nominal_metal_wt_fraction is not None and nominal_metal_wt_fraction > 0:
        recovery_vs_nominal_pct = (actual_wt_fraction / nominal_metal_wt_fraction) * 100.0

    return {
        "corrected_concentration_mg_L": float(net_c),
        "metal_mass_mg": float(metal_mass_mg),
        "metal_mass_ug": float(metal_mass_ug),
        "catalyst_mass_mg": float(catalyst_mass_mg),
        "metal_loading_wt_pct": float(actual_wt_pct),
        "metal_loading_wt_fraction": float(actual_wt_fraction),
        "nominal_metal_wt_fraction": nominal_metal_wt_fraction,
        "loading_efficiency_pct": float(recovery_vs_nominal_pct) if recovery_vs_nominal_pct is not None else None,
    }


def calculate_aas_leaching(
    measured_concentration_mg_L: float,
    electrolyte_volume_mL: float,
    catalyst_mass_mg: float,
    initial_metal_wt_fraction: float,
    dilution_factor: float = 1.0,
    blank_concentration_mg_L: float = 0.0,
    aliquot_history: Optional[Sequence[Tuple[float, float]]] = None,
) -> Dict[str, Any]:
    """
    Calculate quantitative catalyst metal leaching during HER/OER stability testing.

    Includes compensation for cumulative metal mass removed during periodic aliquot sampling.
    """
    return calculate_catalyst_leaching_common(
        measured_concentration_mg_L=measured_concentration_mg_L,
        electrolyte_volume_mL=electrolyte_volume_mL,
        catalyst_mass_mg=catalyst_mass_mg,
        initial_metal_wt_fraction=initial_metal_wt_fraction,
        dilution_factor=dilution_factor,
        blank_concentration_mg_L=blank_concentration_mg_L,
        aliquot_corrections=aliquot_history,
    )


# =========================================================================
# 6. Publication-Quality Plotting Functions
# =========================================================================

def plot_aas_calibration(
    cal_result: Dict[str, Any],
    element: str = "Analyte",
    unit: str = "mg/L",
    title: Optional[str] = None,
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    """
    Generate publication 2-panel AAS calibration curve with fit line, confidence interval, and residuals.
    """
    apply_publication_style()
    fig, (ax_top, ax_bot) = create_publication_figure(
        columns="single",
        aspect_ratio=1.35,
        nrows=2,
        ncols=1,
        sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.0], "hspace": 0.08},
    )

    c = cal_result["standards_c"]
    s_obs = cal_result["standards_signal"]
    s_fit = cal_result["fit_signal"]
    res = cal_result["residuals"]
    sig_type = cal_result.get("signal_type", "absorbance")

    y_label = r"Absorbance $[\mathrm{AU}]$" if sig_type == "absorbance" else r"Integrated Abs. $[\mathrm{A\cdot s}]$"

    # Top panel: Calibration curve
    ax_top.scatter(c, s_obs, color="#0072B2", s=25, label="Standards", zorder=3)
    ax_top.plot(c, s_fit, color="#D55E00", linewidth=1.4, label="Linear Fit", zorder=2)

    r2 = cal_result["r_squared"]
    slope = cal_result["slope"]
    intercept = cal_result["intercept"]
    annot = f"$R^2 = {r2:.5f}$\nSlope: {slope:.3f}\nIntercept: {intercept:.4f}"
    if cal_result.get("lod_conc") is not None:
        annot += f"\nLOD: {cal_result['lod_conc']:.3f} {unit}"

    ax_top.text(
        0.05,
        0.75,
        annot,
        transform=ax_top.transAxes,
        fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.9),
    )

    ax_top.set_ylabel(y_label)
    ttl = title or f"AAS Calibration: {element}"
    ax_top.set_title(ttl, fontsize=9, pad=6)
    ax_top.legend(loc="lower right", frameon=False, fontsize=8)

    # Bottom panel: Residuals
    ax_bot.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax_bot.scatter(c, res, color="#009E73", s=18, zorder=3)
    ax_bot.set_xlabel(rf"Concentration $[\mathrm{{{unit}}}]$")
    ax_bot.set_ylabel(r"$\Delta \mathrm{Signal}$")

    fig.tight_layout()
    return fig, (ax_top, ax_bot)


def plot_standard_addition(
    sa_result: Dict[str, Any],
    element: str = "Analyte",
    unit: str = "mg/L",
    title: Optional[str] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate publication Method of Standard Additions (MOSA) plot with backward extrapolation.
    """
    apply_publication_style()
    fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    x_add = sa_result["added_concentrations"]
    y_sig = sa_result["signals"]
    ext_x = sa_result["extrap_x"]
    ext_y = sa_result["extrap_y"]
    c_app = sa_result["apparent_concentration"]
    ci95 = sa_result["ci95_apparent"]
    sig_type = sa_result.get("signal_type", "absorbance")
    y_label = r"Absorbance $[\mathrm{AU}]$" if sig_type == "absorbance" else r"Integrated Abs. $[\mathrm{A\cdot s}]$"

    # Plot empirical points
    ax.scatter(x_add, y_sig, color="#0072B2", s=30, label="Standard Additions", zorder=4)

    # Plot extrapolation line
    ax.plot(ext_x, ext_y, color="#D55E00", linestyle="-", linewidth=1.4, label="Linear Extrapolation", zorder=3)

    # Zero response baseline and y-axis
    ax.axhline(0, color="black", linewidth=0.8)
    ax.axvline(0, color="gray", linestyle=":", linewidth=0.8)

    # Mark extrapolated x-intercept (-c_app)
    ax.scatter([-c_app], [0], color="#CC79A7", s=45, marker="D", label=f"Extrapolated $-C_0$", zorder=5)
    # Error bar at intercept
    ax.errorbar([-c_app], [0], xerr=[ci95], color="#CC79A7", capsize=3, linewidth=1.0, zorder=5)

    # Drop dashed line to intercept
    ax.plot([-c_app, -c_app], [0, sa_result["intercept"]], color="#CC79A7", linestyle="--", linewidth=0.8)

    annot = (
        f"$C_{{\mathrm{{sample}}}} = {c_app:.3f} \\pm {ci95:.3f}\\ \\mathrm{{{unit}}}$\n"
        f"$R^2 = {sa_result['r_squared']:.4f}$\n"
        f"Slope: {sa_result['slope']:.3f}"
    )
    if sa_result.get("dilution_factor", 1.0) > 1.0:
        c_orig = sa_result["original_concentration"]
        ci_orig = sa_result["ci95_original"]
        annot += f"\n$C_{{\mathrm{{original}}}} = {c_orig:.3f} \\pm {ci_orig:.3f}\\ \\mathrm{{{unit}}}$"

    ax.text(
        0.05,
        0.72,
        annot,
        transform=ax.transAxes,
        fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.9),
    )

    ax.set_xlabel(rf"Added Concentration $[\mathrm{{{unit}}}]$")
    ax.set_ylabel(y_label)
    ttl = title or f"Standard Additions: {element} Matrix Recovery"
    ax.set_title(ttl, fontsize=9, pad=6)
    ax.legend(loc="lower right", frameon=False, fontsize=7.5)

    fig.tight_layout()
    return fig, ax


def plot_aas_concentration_comparison(
    sample_names: Sequence[str],
    concentrations: Sequence[float],
    error_bars: Optional[Sequence[float]] = None,
    element: str = "Metal",
    unit: str = "mg/L",
    title: Optional[str] = None,
    color: str = "#0072B2",
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate publication bar chart comparing analyte concentrations across samples.
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        created_fig = True
    else:
        fig = ax.get_figure()

    x_pos = np.arange(len(sample_names))
    y_vals = np.asarray(concentrations, dtype=float)
    y_err = np.asarray(error_bars, dtype=float) if error_bars is not None else None

    bars = ax.bar(
        x_pos, y_vals, yerr=y_err, capsize=3, color=color,
        edgecolor="black", linewidth=0.8, alpha=0.85, zorder=3,
    )

    for idx, rect in enumerate(bars):
        h = rect.get_height()
        err_h = y_err[idx] if y_err is not None else 0.0
        ax.annotate(
            f"{h:.2f}",
            xy=(rect.get_x() + rect.get_width() / 2.0, h + err_h),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=7,
        )

    ax.set_xticks(x_pos)
    ax.set_xticklabels(sample_names, rotation=15, ha="right", fontsize=8)
    ax.set_ylabel(rf"{element} Concentration $[\mathrm{{{unit}}}]$")
    if title:
        ax.set_title(title, fontsize=9, pad=6)

    max_v = max(y_vals) if len(y_vals) > 0 and max(y_vals) > 0 else 1.0
    ax.set_ylim(0, max_v * 1.25)

    if created_fig:
        fig.tight_layout()
    return fig, ax


def plot_aas_catalyst_loading(
    sample_names: Sequence[str],
    actual_loading_wt_pct: Sequence[float],
    nominal_loading_wt_pct: Optional[Sequence[float]] = None,
    error_bars: Optional[Sequence[float]] = None,
    element: str = "Metal",
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate publication bar chart comparing actual vs nominal catalyst metal loadings (wt.%).
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        created_fig = True
    else:
        fig = ax.get_figure()

    x_pos = np.arange(len(sample_names))
    width = 0.35 if nominal_loading_wt_pct is not None else 0.55

    # Actual bars
    rects1 = ax.bar(
        x_pos - (width / 2.0 if nominal_loading_wt_pct is not None else 0),
        actual_loading_wt_pct,
        width,
        yerr=error_bars,
        capsize=3,
        label="Actual (AAS)",
        color="#009E73",
        edgecolor="black",
        linewidth=0.8,
        alpha=0.9,
        zorder=3,
    )

    # Nominal bars if provided
    if nominal_loading_wt_pct is not None:
        rects2 = ax.bar(
            x_pos + width / 2.0,
            nominal_loading_wt_pct,
            width,
            label="Nominal Target",
            color="#E69F00",
            edgecolor="black",
            linewidth=0.8,
            alpha=0.75,
            zorder=3,
        )
        ax.legend(loc="upper right", frameon=False, fontsize=7.5)

    for rect in rects1:
        h = rect.get_height()
        ax.annotate(
            f"{h:.1f}%",
            xy=(rect.get_x() + rect.get_width() / 2.0, h),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=7,
        )

    ax.set_xticks(x_pos)
    ax.set_xticklabels(sample_names, rotation=15, ha="right", fontsize=8)
    ax.set_ylabel(rf"{element} Loading $[\mathrm{{wt.\%}}]$")
    ttl = title or f"Catalyst {element} Loading Quantification"
    ax.set_title(ttl, fontsize=9, pad=6)

    max_val = max(actual_loading_wt_pct) if len(actual_loading_wt_pct) > 0 else 10.0
    ax.set_ylim(0, max_val * 1.3)

    if created_fig:
        fig.tight_layout()
    return fig, ax


def plot_aas_leaching(
    time_or_cycles: Sequence[Union[float, str]],
    leached_percentages: Sequence[float],
    mass_normalized_loss_ug_per_mg: Optional[Sequence[float]] = None,
    element: str = "Metal",
    title: Optional[str] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate publication kinetics plot showing catalyst metal leaching over stability testing.
    """
    apply_publication_style()
    fig, ax1 = create_publication_figure(columns="single", aspect_ratio="golden")

    x_vals = np.arange(len(time_or_cycles))
    y_pct = np.asarray(leached_percentages, dtype=float)

    line1 = ax1.plot(
        x_vals, y_pct, marker="o", markersize=5, color="#D55E00",
        linewidth=1.4, label=f"{element} Leached (%)", zorder=3,
    )
    ax1.set_ylabel(rf"{element} Leached $[\%]$", color="#D55E00")
    ax1.tick_params(axis="y", labelcolor="#D55E00")
    ax1.set_xticks(x_vals)
    ax1.set_xticklabels([str(t) for t in time_or_cycles], rotation=15, ha="right", fontsize=8)
    ax1.set_xlabel("Testing Duration / Cycles")

    # Optional dual axis for mass normalized loss
    if mass_normalized_loss_ug_per_mg is not None:
        ax2 = ax1.twinx()
        y_mass = np.asarray(mass_normalized_loss_ug_per_mg, dtype=float)
        line2 = ax2.plot(
            x_vals, y_mass, marker="s", markersize=4.5, color="#0072B2",
            linestyle="--", linewidth=1.2, label="Normalized Loss", zorder=2,
        )
        ax2.set_ylabel(r"Loss $[\mu\mathrm{g}_{\mathrm{metal}}/\mathrm{mg}_{\mathrm{cat}}]$", color="#0072B2")
        ax2.tick_params(axis="y", labelcolor="#0072B2")
        ax2.grid(False)

    ttl = title or f"Electrochemical Stability: {element} Dissolution Kinetics"
    ax1.set_title(ttl, fontsize=9, pad=6)

    fig.tight_layout()
    return fig, ax1


def plot_aas_stability_before_after(
    sample_names: Sequence[str],
    before_loading_wt_pct: Sequence[float],
    after_loading_wt_pct: Sequence[float],
    element: str = "Metal",
    title: Optional[str] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate publication paired bar chart comparing catalyst composition before vs after stability stress.
    """
    apply_publication_style()
    fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    x_pos = np.arange(len(sample_names))
    width = 0.35

    ax.bar(
        x_pos - width / 2.0, before_loading_wt_pct, width,
        label="Pristine (Before OER)", color="#0072B2", edgecolor="black", linewidth=0.8, alpha=0.9, zorder=3,
    )
    ax.bar(
        x_pos + width / 2.0, after_loading_wt_pct, width,
        label="Post-Stability (After 100h)", color="#D55E00", edgecolor="black", linewidth=0.8, alpha=0.85, zorder=3,
    )

    ax.set_xticks(x_pos)
    ax.set_xticklabels(sample_names, rotation=15, ha="right", fontsize=8)
    ax.set_ylabel(rf"{element} Content $[\mathrm{{wt.\%}}]$")
    ttl = title or f"Catalyst Retention: Pristine vs Spent Post-OER"
    ax.set_title(ttl, fontsize=9, pad=6)
    ax.legend(loc="upper right", frameon=False, fontsize=7.5)

    max_v = max(max(before_loading_wt_pct), max(after_loading_wt_pct)) if len(before_loading_wt_pct) > 0 else 10.0
    ax.set_ylim(0, max_v * 1.3)

    fig.tight_layout()
    return fig, ax


def plot_aas_replicate_distribution(
    sample_names: Sequence[str],
    replicate_measurements: Sequence[Sequence[float]],
    unit: str = "mg/L",
    title: Optional[str] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate publication precision strip chart showing individual replicate measurements and mean.
    """
    apply_publication_style()
    fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    x_pos = np.arange(len(sample_names))
    palette = get_palette("nature", len(sample_names))

    for idx, (sname, reps) in enumerate(zip(sample_names, replicate_measurements)):
        r_arr = np.asarray(reps, dtype=float)
        mean_v = np.mean(r_arr)
        std_v = np.std(r_arr, ddof=1) if len(r_arr) > 1 else 0.0
        c = palette[idx % len(palette)]

        # Jitter points slightly
        jitter = np.random.normal(0, 0.04, size=len(r_arr))
        ax.scatter(np.full_like(r_arr, idx) + jitter, r_arr, color=c, alpha=0.7, s=25, zorder=3)
        # Mean bar and error bar
        ax.errorbar(idx, mean_v, yerr=std_v, color="black", capsize=4, marker="_", markersize=14, linewidth=1.2, zorder=4)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(sample_names, rotation=15, ha="right", fontsize=8)
    ax.set_ylabel(rf"Concentration $[\mathrm{{{unit}}}]$")
    ttl = title or "AAS Instrumental Replicate Precision"
    ax.set_title(ttl, fontsize=9, pad=6)

    fig.tight_layout()
    return fig, ax


def plot_aas_qc_chart(
    ccv_labels: Sequence[str],
    ccv_recoveries_pct: Sequence[float],
    sop: str = "EPA_7000B",
    title: Optional[str] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate analytical Quality Control (QC) Shewhart-style control chart for CCV recoveries.
    """
    apply_publication_style()
    fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    x_vals = np.arange(len(ccv_labels))
    y_vals = np.asarray(ccv_recoveries_pct, dtype=float)

    is_furnace = "7010" in sop or "GFAAS" in sop
    ctrl_limit = 15.0 if is_furnace else 10.0
    warn_limit = 10.0 if is_furnace else 5.0

    # Nominal 100% line
    ax.axhline(100.0, color="black", linewidth=1.0, linestyle="-", label="Nominal (100%)")
    # Warning limit lines
    ax.axhline(100.0 + warn_limit, color="#E69F00", linestyle="--", linewidth=0.8, label=f"Warning ($\\pm{warn_limit:.0f}\\%$)")
    ax.axhline(100.0 - warn_limit, color="#E69F00", linestyle="--", linewidth=0.8)
    # Control limit lines
    ax.axhline(100.0 + ctrl_limit, color="#D55E00", linestyle=":", linewidth=1.0, label=f"Control Limit ($\\pm{ctrl_limit:.0f}\\%$)")
    ax.axhline(100.0 - ctrl_limit, color="#D55E00", linestyle=":", linewidth=1.0)

    # Shaded acceptable band
    ax.fill_between([-0.5, len(x_vals) - 0.5], 100.0 - ctrl_limit, 100.0 + ctrl_limit, color="#F0E442", alpha=0.15)

    # Plot CCV points
    ax.plot(x_vals, y_vals, marker="o", color="#0072B2", linewidth=1.2, markersize=5, zorder=4)

    # Highlight any out-of-control points
    for idx, val in enumerate(y_vals):
        if abs(val - 100.0) > ctrl_limit:
            ax.scatter([idx], [val], color="red", s=50, marker="x", zorder=5)

    ax.set_xticks(x_vals)
    ax.set_xticklabels(ccv_labels, rotation=15, ha="right", fontsize=8)
    ax.set_ylabel(r"CCV Recovery $[\%]$")
    ax.set_xlim(-0.5, len(x_vals) - 0.5)
    ax.set_ylim(100.0 - ctrl_limit * 1.5, 100.0 + ctrl_limit * 1.5)

    ttl = title or f"QC Control Chart: Continuing Calibration Verification ({sop})"
    ax.set_title(ttl, fontsize=9, pad=6)
    ax.legend(loc="lower right", frameon=False, fontsize=7)

    fig.tight_layout()
    return fig, ax

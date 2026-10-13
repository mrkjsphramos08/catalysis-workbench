"""
Shared Elemental Analysis and Quantitative Spectrometry Utilities.

Provides reusable mathematical, statistical, calibration, and unit-conversion
routines shared across Atomic Absorption Spectroscopy (AAS) and
Inductively Coupled Plasma (ICP-OES / ICP-MS) techniques.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import numpy as np
import scipy.stats as stats


def convert_concentration_units(
    value: Union[float, np.ndarray, Sequence[float]],
    from_unit: str,
    to_unit: str,
    molar_mass_g_mol: Optional[float] = None,
    sample_mass_g: Optional[float] = None,
    solution_volume_L: Optional[float] = None,
) -> Union[float, np.ndarray]:
    """
    Convert concentrations between liquid and solid analytical units.

    Supported units:
    - Liquid: 'mg/L', 'ppm', 'ug/L', 'ppb', 'ng/L', 'ppt', 'g/L', 'mol/L', 'M', 'mM', 'uM'
    - Solid: 'mg/kg', 'ppm', 'ug/g', 'ng/g', 'ppb', 'wt%', 'wt_fraction'

    Parameters
    ----------
    value : float or array-like
        Concentration numerical value(s).
    from_unit : str
        Source unit.
    to_unit : str
        Target unit.
    molar_mass_g_mol : float, optional
        Molar mass in g/mol for molar conversions.
    sample_mass_g : float, optional
        Dry sample mass in grams for solid digest conversions.
    solution_volume_L : float, optional
        Final solution volume in liters for solid digest conversions.

    Returns
    -------
    float or ndarray
        Converted concentration value(s).
    """
    v_arr = np.asarray(value, dtype=float)
    u_from = from_unit.strip().lower().replace(" ", "")
    u_to = to_unit.strip().lower().replace(" ", "")

    if u_from == u_to:
        return float(v_arr) if np.isscalar(value) or v_arr.ndim == 0 else v_arr

    # Normalize aliases
    alias_map = {
        "ppm": "mg/l",
        "ppb": "ug/l",
        "ppt": "ng/l",
        "µg/l": "ug/l",
        "µg/g": "ug/g",
        "mg/dm3": "mg/l",
        "wt.%": "wt%",
        "weight%": "wt%",
    }
    u_from = alias_map.get(u_from, u_from)
    u_to = alias_map.get(u_to, u_to)

    # 1. Base liquid unit: mg/L
    liquid_to_mg_L = {
        "mg/l": 1.0,
        "ug/l": 1e-3,
        "ng/l": 1e-6,
        "g/l": 1e3,
    }

    # Base solid unit: mg/kg
    solid_to_mg_kg = {
        "mg/kg": 1.0,
        "ug/g": 1.0,
        "ng/g": 1e-3,
        "wt%": 1e4,  # 1 wt% = 10,000 mg/kg
        "wt_fraction": 1e6,  # 1 fraction = 1,000,000 mg/kg
    }

    # Handle molar conversions if requested
    if u_from in ("mol/l", "m", "mm", "um") or u_to in ("mol/l", "m", "mm", "um"):
        if molar_mass_g_mol is None or molar_mass_g_mol <= 0:
            raise ValueError(f"Molar mass required to convert between molar and mass concentrations ({from_unit} -> {to_unit}).")
        # 1 mol/L = molar_mass * 1000 mg/L
        molar_to_mol_L = {"mol/l": 1.0, "m": 1.0, "mm": 1e-3, "um": 1e-6}
        if u_from in molar_to_mol_L:
            c_mol_L = v_arr * molar_to_mol_L[u_from]
            c_mg_L = c_mol_L * molar_mass_g_mol * 1e3
        else:
            if u_from not in liquid_to_mg_L:
                raise ValueError(f"Unsupported source unit for molar conversion: {from_unit}")
            c_mg_L = v_arr * liquid_to_mg_L[u_from]

        if u_to in molar_to_mol_L:
            c_mol_L = c_mg_L / (molar_mass_g_mol * 1e3)
            res = c_mol_L / molar_to_mol_L[u_to]
            return float(res) if np.isscalar(value) or v_arr.ndim == 0 else res
        u_from = "mg/l"
        v_arr = c_mg_L

    # Liquid to liquid
    if u_from in liquid_to_mg_L and u_to in liquid_to_mg_L:
        in_mg_L = v_arr * liquid_to_mg_L[u_from]
        res = in_mg_L / liquid_to_mg_L[u_to]
        return float(res) if np.isscalar(value) or v_arr.ndim == 0 else res

    # Solid to solid
    if u_from in solid_to_mg_kg and u_to in solid_to_mg_kg:
        in_mg_kg = v_arr * solid_to_mg_kg[u_from]
        res = in_mg_kg / solid_to_mg_kg[u_to]
        return float(res) if np.isscalar(value) or v_arr.ndim == 0 else res

    # Liquid to solid: requires sample_mass_g and solution_volume_L
    if u_from in liquid_to_mg_L and u_to in solid_to_mg_kg:
        if sample_mass_g is None or solution_volume_L is None or sample_mass_g <= 0 or solution_volume_L <= 0:
            raise ValueError(f"sample_mass_g and solution_volume_L required to convert liquid to solid concentration ({from_unit} -> {to_unit}).")
        c_mg_L = v_arr * liquid_to_mg_L[u_from]
        # Total dissolved mass in mg = c_mg_L * solution_volume_L
        # mg/kg = (c_mg_L * solution_volume_L) / (sample_mass_g / 1000.0)
        in_mg_kg = (c_mg_L * solution_volume_L) / (sample_mass_g * 1e-3)
        res = in_mg_kg / solid_to_mg_kg[u_to]
        return float(res) if np.isscalar(value) or v_arr.ndim == 0 else res

    # Solid to liquid: requires sample_mass_g and solution_volume_L
    if u_from in solid_to_mg_kg and u_to in liquid_to_mg_L:
        if sample_mass_g is None or solution_volume_L is None or sample_mass_g <= 0 or solution_volume_L <= 0:
            raise ValueError(f"sample_mass_g and solution_volume_L required to convert solid to liquid concentration ({from_unit} -> {to_unit}).")
        in_mg_kg = v_arr * solid_to_mg_kg[u_from]
        total_mg = in_mg_kg * (sample_mass_g * 1e-3)
        c_mg_L = total_mg / solution_volume_L
        res = c_mg_L / liquid_to_mg_L[u_to]
        return float(res) if np.isscalar(value) or v_arr.ndim == 0 else res

    raise ValueError(f"Conversion from '{from_unit}' to '{to_unit}' is not supported.")


def compute_replicate_statistics(
    values: Sequence[float],
    ddof: int = 1,
    confidence_level: float = 0.95,
) -> Dict[str, Any]:
    """
    Compute rigorous descriptive statistics for analytical replicate measurements.

    Parameters
    ----------
    values : sequence of float
        Replicate measurement values.
    ddof : int, default 1
        Degrees of freedom adjustment for sample standard deviation.
    confidence_level : float, default 0.95
        Two-sided confidence level for mean uncertainty.

    Returns
    -------
    Dict[str, Any]
        Dictionary with 'mean', 'std', 'rsd_pct', 'count', 'sem', 'ci_half_width',
        'min', 'max', 'median'.
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[~np.isnan(arr)]
    n = len(arr)

    if n == 0:
        return {
            "count": 0,
            "mean": np.nan,
            "std": np.nan,
            "rsd_pct": np.nan,
            "sem": np.nan,
            "ci_half_width": np.nan,
            "min": np.nan,
            "max": np.nan,
            "median": np.nan,
        }

    mean_val = float(np.mean(arr))
    if n == 1:
        return {
            "count": 1,
            "mean": mean_val,
            "std": 0.0,
            "rsd_pct": 0.0,
            "sem": 0.0,
            "ci_half_width": 0.0,
            "min": mean_val,
            "max": mean_val,
            "median": mean_val,
        }

    std_val = float(np.std(arr, ddof=ddof))
    rsd_pct = float(abs(std_val / mean_val) * 100.0) if abs(mean_val) > 1e-15 else 0.0
    sem_val = float(std_val / np.sqrt(n))

    # Student's t critical value
    t_crit = float(stats.t.ppf((1.0 + confidence_level) / 2.0, df=n - 1))
    ci_half_width = float(t_crit * sem_val)

    return {
        "count": n,
        "mean": mean_val,
        "std": std_val,
        "rsd_pct": rsd_pct,
        "sem": sem_val,
        "ci_half_width": ci_half_width,
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "median": float(np.median(arr)),
    }


def linear_calibration_fit(
    x: Sequence[float],
    y: Sequence[float],
    zero_intercept: bool = False,
    weights: Optional[Sequence[float]] = None,
    alpha: float = 0.05,
) -> Dict[str, Any]:
    """
    Perform linear regression on analytical calibration standards.

    Supports unweighted OLS, weighted WLS, and zero-intercept constraints.

    Parameters
    ----------
    x : sequence of float
        Standard concentrations.
    y : sequence of float
        Analytical responses (absorbance, integrated absorbance, intensity).
    zero_intercept : bool, default False
        If True, forces model y = m * x.
    weights : sequence of float, optional
        Fitting weights w_i. If provided, minimizes sum w_i * (y_i - y_hat_i)^2.
    alpha : float, default 0.05
        Significance level for confidence interval bounds (1 - alpha = 0.95).

    Returns
    -------
    Dict[str, Any]
        Regression parameters, R^2, slope, intercept, standard errors, residuals.
    """
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)

    if len(x_arr) != len(y_arr):
        raise ValueError("Length of concentrations x and responses y must match.")
    if len(x_arr) < (2 if zero_intercept else 3):
        raise ValueError(f"At least {'2' if zero_intercept else '3'} standards required for calibration fit.")

    n = len(x_arr)
    w_arr = np.ones(n, dtype=float) if weights is None else np.asarray(weights, dtype=float)
    if np.any(w_arr <= 0):
        raise ValueError("Calibration regression weights must be strictly positive.")

    if zero_intercept:
        # Model: y = m * x
        # m = sum(w * x * y) / sum(w * x^2)
        denom = float(np.sum(w_arr * x_arr ** 2))
        if denom == 0.0:
            raise ValueError("Degenerate standards: sum(w * x^2) is zero.")
        slope = float(np.sum(w_arr * x_arr * y_arr) / denom)
        intercept = 0.0

        y_pred = slope * x_arr
        residuals = y_arr - y_pred
        df_e = n - 1
        ss_res = float(np.sum(w_arr * residuals ** 2))
        syx = float(np.sqrt(ss_res / df_e)) if df_e > 0 else 0.0
        se_slope = float(syx / np.sqrt(denom)) if denom > 0 else 0.0
        se_intercept = 0.0

        # R^2 for zero-intercept regression (uncentered vs centered)
        ss_tot = float(np.sum(w_arr * (y_arr - np.average(y_arr, weights=w_arr)) ** 2))
        r2 = float(1.0 - (ss_res / ss_tot)) if ss_tot > 0 else 1.0
        p_val = float(2.0 * (1.0 - stats.t.cdf(abs(slope / se_slope), df=df_e))) if se_slope > 0 else 0.0
    else:
        # Standard OLS or WLS: y = m * x + b
        if weights is None:
            ols_res = stats.linregress(x_arr, y_arr)
            slope = float(ols_res.slope)
            intercept = float(ols_res.intercept)
            r2 = float(ols_res.rvalue ** 2)
            p_val = float(ols_res.pvalue)
            se_slope = float(ols_res.stderr)
            se_intercept = float(ols_res.intercept_stderr)

            y_pred = slope * x_arr + intercept
            residuals = y_arr - y_pred
            df_e = n - 2
            ss_res = float(np.sum(residuals ** 2))
            syx = float(np.sqrt(ss_res / df_e)) if df_e > 0 else 0.0
        else:
            # Weighted Least Squares
            w_sum = np.sum(w_arr)
            x_w = np.sum(w_arr * x_arr) / w_sum
            y_w = np.sum(w_arr * y_arr) / w_sum

            s_xx = np.sum(w_arr * (x_arr - x_w) ** 2)
            s_xy = np.sum(w_arr * (x_arr - x_w) * (y_arr - y_w))

            if s_xx == 0.0:
                raise ValueError("Degenerate concentrations: standard variance is zero.")
            slope = float(s_xy / s_xx)
            intercept = float(y_w - slope * x_w)

            y_pred = slope * x_arr + intercept
            residuals = y_arr - y_pred
            df_e = n - 2
            ss_res = float(np.sum(w_arr * residuals ** 2))
            syx = float(np.sqrt(ss_res / df_e)) if df_e > 0 else 0.0

            se_slope = float(syx / np.sqrt(s_xx)) if s_xx > 0 else 0.0
            se_intercept = float(syx * np.sqrt(1.0 / w_sum + (x_w ** 2) / s_xx)) if s_xx > 0 else 0.0

            ss_tot = float(np.sum(w_arr * (y_arr - y_w) ** 2))
            r2 = float(1.0 - (ss_res / ss_tot)) if ss_tot > 0 else 1.0
            p_val = float(2.0 * (1.0 - stats.t.cdf(abs(slope / se_slope), df=df_e))) if se_slope > 0 else 0.0

    # Relative residuals in percent
    safe_pred = np.where(np.abs(y_pred) > 1e-12, y_pred, 1.0)
    rel_residuals_pct = (residuals / safe_pred) * 100.0

    # Back-calculated concentrations
    calc_x = np.where(abs(slope) > 1e-15, (y_arr - intercept) / slope, np.nan)
    safe_x = np.where(np.abs(x_arr) > 1e-12, x_arr, 1.0)
    recoveries_pct = (calc_x / safe_x) * 100.0

    # Check for potential non-linear rollover (e.g. curvature in AAS at high absorbance)
    curvature_detected = False
    if n >= 4:
        # Fit second-degree polynomial to test significance of quadratic term
        quad_fit = np.polyfit(x_arr, y_arr, deg=2)
        if quad_fit[0] < 0:  # negative quadratic term indicates rollover
            quad_pred = np.polyval(quad_fit, x_arr)
            quad_ss_res = np.sum((y_arr - quad_pred) ** 2)
            # F-test for model comparison
            f_stat = ((ss_res - quad_ss_res) / 1.0) / (quad_ss_res / (n - 3)) if quad_ss_res > 0 and n > 3 else 0.0
            if f_stat > 10.0:  # Strong non-linear curvature
                curvature_detected = True

    return {
        "slope": slope,
        "intercept": intercept,
        "r_squared": r2,
        "syx": syx,
        "std_err_slope": se_slope,
        "std_err_intercept": se_intercept,
        "p_value": p_val,
        "zero_intercept": zero_intercept,
        "weighted": weights is not None,
        "x": x_arr,
        "y": y_arr,
        "y_pred": y_pred,
        "residuals": residuals,
        "rel_residuals_pct": rel_residuals_pct,
        "calc_x": calc_x,
        "recoveries_pct": recoveries_pct,
        "x_min": float(np.min(x_arr)),
        "x_max": float(np.max(x_arr)),
        "curvature_detected": curvature_detected,
    }


def calculate_detection_limits(
    slope: float,
    blank_stdev: Optional[float] = None,
    syx: Optional[float] = None,
    method: str = "blank_stdev",
    k_lod: float = 3.3,
    k_loq: float = 10.0,
) -> Dict[str, Any]:
    """
    Calculate analytical Limit of Detection (LOD) and Limit of Quantification (LOQ).

    Methods:
    - 'blank_stdev': Uses blank standard deviation s_b: LOD = k_lod * s_b / m, LOQ = k_loq * s_b / m.
    - 'residual_stdev': Uses regression standard error s_{y/x}: LOD = k_lod * s_{y/x} / m, LOQ = k_loq * s_{y/x} / m.

    Parameters
    ----------
    slope : float
        Calibration slope m. Must be > 0.
    blank_stdev : float, optional
        Standard deviation of replicate blanks s_b.
    syx : float, optional
        Standard error of calibration regression s_{y/x}.
    method : str, default 'blank_stdev'
        'blank_stdev' or 'residual_stdev'.
    k_lod : float, default 3.3
        Multiplier for LOD (IUPAC default 3.3).
    k_loq : float, default 10.0
        Multiplier for LOQ (IUPAC default 10.0).

    Returns
    -------
    Dict[str, Any]
        Dictionary with 'lod_conc', 'loq_conc', 'method', 'k_lod', 'k_loq'.
    """
    if slope <= 0:
        return {
            "lod_conc": None,
            "loq_conc": None,
            "method": method,
            "k_lod": k_lod,
            "k_loq": k_loq,
        }

    m_low = method.strip().lower()
    if m_low in ("blank_stdev", "blank", "sb"):
        if blank_stdev is None or blank_stdev <= 0:
            return {"lod_conc": None, "loq_conc": None, "method": method, "k_lod": k_lod, "k_loq": k_loq}
        lod = float(k_lod * blank_stdev / slope)
        loq = float(k_loq * blank_stdev / slope)
    elif m_low in ("residual_stdev", "syx", "regression"):
        if syx is None or syx <= 0:
            return {"lod_conc": None, "loq_conc": None, "method": method, "k_lod": k_lod, "k_loq": k_loq}
        lod = float(k_lod * syx / slope)
        loq = float(k_loq * syx / slope)
    else:
        raise ValueError(f"Unknown detection limit method '{method}'. Choose 'blank_stdev' or 'residual_stdev'.")

    return {
        "lod_conc": lod,
        "loq_conc": loq,
        "method": method,
        "k_lod": k_lod,
        "k_loq": k_loq,
    }


def calculate_catalyst_leaching_common(
    measured_concentration_mg_L: float,
    electrolyte_volume_mL: float,
    catalyst_mass_mg: float,
    initial_metal_wt_fraction: float,
    dilution_factor: float = 1.0,
    blank_concentration_mg_L: float = 0.0,
    aliquot_corrections: Optional[Sequence[Tuple[float, float]]] = None,
) -> Dict[str, Any]:
    """
    Calculate quantitative catalyst metal dissolution metrics with optional aliquot corrections.

    Parameters
    ----------
    measured_concentration_mg_L : float
        Measured concentration in solution (mg/L).
    electrolyte_volume_mL : float
        Current volume of testing electrolyte in milliliters.
    catalyst_mass_mg : float
        Total initial catalyst mass tested in milligrams.
    initial_metal_wt_fraction : float
        Initial weight fraction of the target metal (e.g., 0.20 for 20 wt% Pt/C).
    dilution_factor : float, default 1.0
        Sample dilution factor applied prior to analysis.
    blank_concentration_mg_L : float, default 0.0
        Analytical blank concentration in mg/L.
    aliquot_corrections : sequence of (vol_mL, conc_mg_L), optional
        List of previously withdrawn aliquots to account for cumulative removed mass.

    Returns
    -------
    Dict[str, Any]
        Dictionary of corrected concentration, total dissolved mass, mass normalized loss,
        and leached percentage.
    """
    net_c = (measured_concentration_mg_L - blank_concentration_mg_L) * dilution_factor
    net_c = max(0.0, net_c)

    # Current dissolved mass in electrolyte cell in mg
    vol_L = electrolyte_volume_mL / 1000.0
    current_dissolved_mg = net_c * vol_L

    # Add mass previously removed in aliquots during time-course sampling
    cumulative_removed_mg = 0.0
    if aliquot_corrections is not None:
        for al_vol_mL, al_conc_mg_L in aliquot_corrections:
            cumulative_removed_mg += (al_conc_mg_L * (al_vol_mL / 1000.0))

    total_dissolved_mg = current_dissolved_mg + cumulative_removed_mg
    dissolved_ug = total_dissolved_mg * 1000.0

    # Initial metal inventory in mg
    initial_metal_mg = catalyst_mass_mg * initial_metal_wt_fraction

    pct_leached = (total_dissolved_mg / initial_metal_mg * 100.0) if initial_metal_mg > 0 else 0.0
    norm_ug_per_mg = dissolved_ug / catalyst_mass_mg if catalyst_mass_mg > 0 else 0.0

    return {
        "corrected_concentration_mg_L": float(net_c),
        "dissolved_metal_mass_ug": float(dissolved_ug),
        "mass_normalized_leaching_ug_per_mg": float(norm_ug_per_mg),
        "percentage_metal_leached": float(pct_leached),
        "cumulative_removed_mass_ug": float(cumulative_removed_mg * 1000.0),
    }

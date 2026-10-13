"""
Shared Peak Fitting, Baseline Correction, and Signal Processing Engine.

Provides unified numerical algorithms using lmfit, pybaselines, and scipy:
- Multi-peak deconvolution (Gaussian, Lorentzian, pseudo-Voigt)
- Baseline estimation and background subtraction (AsLS, SNIP, ModPoly, AirPLS)
- Automated peak finding with prominence and distance constraints
- Analytical and numerical peak integration
- Confidence bounds and regression diagnostics (R^2, reduced chi-square, residuals)
- Savitzky-Golay smoothing and cosmic-ray spike detection.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import numpy as np
import scipy.signal as sp_signal
from scipy.integrate import simpson

import lmfit
from lmfit.models import GaussianModel, LorentzianModel, PseudoVoigtModel, ConstantModel, LinearModel

try:
    import pybaselines
    from pybaselines import Baseline
    HAS_PYBASELINES = True
except ImportError:
    HAS_PYBASELINES = False


def detect_peaks(
    x: np.ndarray,
    y: np.ndarray,
    prominence: Optional[float] = None,
    distance: Optional[int] = None,
    width: Optional[float] = None,
    height: Optional[float] = None,
) -> Dict[str, np.ndarray]:
    """
    Identify significant local maxima in a 1D spectral coordinate series.

    Parameters
    ----------
    x : np.ndarray
        Abscissa array (e.g., 2theta, wavenumber, Raman shift).
    y : np.ndarray
        Signal intensity array.
    prominence : float, optional
        Minimum relative peak prominence. Defaults to 5% of signal ptp.
    distance : int, optional
        Minimum separation in index steps between adjacent peaks.
    width : float, optional
        Minimum peak width in index steps.
    height : float, optional
        Minimum absolute intensity required.

    Returns
    -------
    Dict[str, np.ndarray]
        - 'indices': Index locations of detected peaks.
        - 'x_peaks': Coordinates along x of detected peaks.
        - 'y_peaks': Intensities at detected peaks.
        - 'prominences': Peak prominence values.
    """
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)

    if prominence is None:
        prominence = 0.05 * np.ptp(y_arr) if np.ptp(y_arr) > 0 else 0.01

    indices, props = sp_signal.find_peaks(
        y_arr,
        prominence=prominence,
        distance=distance,
        width=width,
        height=height,
    )

    return {
        "indices": indices,
        "x_peaks": x_arr[indices],
        "y_peaks": y_arr[indices],
        "prominences": props.get("prominences", np.array([])),
    }


def estimate_baseline(
    x: np.ndarray,
    y: np.ndarray,
    method: str = "asls",
    lam: float = 1e5,
    p: float = 0.01,
    poly_order: int = 3,
    **kwargs: Any,
) -> np.ndarray:
    """
    Estimate spectral baseline using pybaselines algorithms with analytical polynomial fallback.

    Parameters
    ----------
    x : np.ndarray
        Abscissa array.
    y : np.ndarray
        Signal intensity array.
    method : str, default 'asls'
        Algorithm: 'asls' (Asymmetric Least Squares), 'snip' (Statistics-sensitive Non-linear
        Iterative Peak-clipping), 'modpoly' (Modified Polynomial), 'airpls', or 'polynomial'.
    lam : float, default 1e5
        Smoothness parameter for spline/whittaker methods.
    p : float, default 0.01
        Asymmetry weighting parameter.
    poly_order : int, default 3
        Order for polynomial-based baselines.

    Returns
    -------
    np.ndarray
        Estimated baseline array with identical shape as y.
    """
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)

    if HAS_PYBASELINES:
        baseline_fitter = Baseline(x_data=x_arr)
        m = method.lower()
        if m == "asls":
            b, _ = baseline_fitter.asls(y_arr, lam=lam, p=p)
            return b
        elif m == "snip":
            max_half_window = kwargs.get("max_half_window", 15)
            b, _ = baseline_fitter.snip(y_arr, max_half_window=max_half_window)
            return b
        elif m in ("modpoly", "polynomial"):
            b, _ = baseline_fitter.modpoly(y_arr, poly_order=poly_order)
            return b
        elif m == "airpls":
            b, _ = baseline_fitter.airpls(y_arr, lam=lam)
            return b
        else:
            b, _ = baseline_fitter.asls(y_arr, lam=lam, p=p)
            return b
    else:
        # Graceful fallback: iterative polynomial thresholding
        coefs = np.polyfit(x_arr, y_arr, deg=poly_order)
        return np.polyval(coefs, x_arr)


def subtract_baseline(
    x: np.ndarray,
    y: np.ndarray,
    method: str = "asls",
    **kwargs: Any,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Subtract estimated baseline from signal, preserving non-negativity where appropriate.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (y_corrected, baseline)
    """
    base = estimate_baseline(x, y, method=method, **kwargs)
    corrected = np.asarray(y, dtype=float) - base
    return corrected, base


def remove_cosmic_rays(
    y: np.ndarray,
    threshold: float = 5.0,
    kernel_size: int = 5,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Detect and correct sharp cosmic-ray spikes using median filtering residual z-scores.

    Parameters
    ----------
    y : np.ndarray
        Raw intensity series.
    threshold : float, default 5.0
        Z-score threshold on the difference between raw and median-filtered signals.
    kernel_size : int, default 5
        Window size for median filter (must be odd).

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (y_despiked, spike_mask_indices)
    """
    y_arr = np.asarray(y, dtype=float)
    if kernel_size % 2 == 0:
        kernel_size += 1

    med = sp_signal.medfilt(y_arr, kernel_size=kernel_size)
    diff = y_arr - med
    mad = np.median(np.abs(diff - np.median(diff)))
    if mad == 0:
        mad = np.std(diff) + 1e-12

    z_score = diff / (1.4826 * mad)
    spike_idx = np.where(z_score > threshold)[0]

    y_clean = y_arr.copy()
    y_clean[spike_idx] = med[spike_idx]
    return y_clean, spike_idx


def smooth_spectrum(
    y: np.ndarray,
    window_length: int = 9,
    polyorder: int = 2,
) -> np.ndarray:
    """Apply Savitzky-Golay smoothing to a 1D spectral series."""
    y_arr = np.asarray(y, dtype=float)
    if window_length >= len(y_arr):
        window_length = len(y_arr) - 1 if (len(y_arr) - 1) % 2 != 0 else len(y_arr) - 2
    if window_length <= polyorder:
        window_length = polyorder + 2 if (polyorder + 2) % 2 != 0 else polyorder + 3
    if window_length % 2 == 0:
        window_length += 1
    return sp_signal.savgol_filter(y_arr, window_length=window_length, polyorder=polyorder)


def fit_peaks(
    x: np.ndarray,
    y: np.ndarray,
    peak_centers: Sequence[float],
    model_type: str = "pseudo_voigt",
    background: str = "linear",
    max_fwhm: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Perform rigorous multi-peak nonlinear least-squares deconvolution using lmfit.

    Parameters
    ----------
    x : np.ndarray
        Abscissa data.
    y : np.ndarray
        Intensity data.
    peak_centers : sequence of float
        Initial center estimates for each peak component.
    model_type : str, default 'pseudo_voigt'
        'gaussian', 'lorentzian', or 'pseudo_voigt'.
    background : str, default 'linear'
        'constant', 'linear', or 'none'.
    max_fwhm : float, optional
        Upper bound constraint for peak full-width at half-maximum.

    Returns
    -------
    Dict[str, Any]
        Dictionary of fit summary, including:
        - 'fit_y': Modeled y-values across input x.
        - 'residuals': Difference (y - fit_y).
        - 'r_squared': Coefficient of determination.
        - 'reduced_chi_sq': Reduced chi-squared metric.
        - 'peaks': List of Dicts with parameters for each component:
          ('center', 'amplitude', 'fWHM', 'height', 'area', 'fraction' if pseudo-Voigt).
        - 'components': Dict mapping component names to individual curve arrays.
        - 'lmfit_result': Full lmfit ModelResult instance.
    """
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)

    if len(peak_centers) == 0:
        raise ValueError("At least one initial peak center must be specified for fit_peaks.")

    # Initialize composite model
    composite_model = None
    params = lmfit.Parameters()

    # Add background component
    if background == "linear":
        bg_model = LinearModel(prefix="bg_")
        params.update(bg_model.make_params(slope=0.0, intercept=float(np.min(y_arr))))
        composite_model = bg_model
    elif background == "constant":
        bg_model = ConstantModel(prefix="bg_")
        params.update(bg_model.make_params(c=float(np.min(y_arr))))
        composite_model = bg_model

    x_range = np.ptp(x_arr)
    initial_sigma = (x_range / (len(peak_centers) * 5.0)) if x_range > 0 else 1.0

    for idx, c_init in enumerate(peak_centers):
        prefix = f"p{idx+1}_"
        m_lower = model_type.lower()
        if m_lower == "gaussian":
            pk_model = GaussianModel(prefix=prefix)
        elif m_lower == "lorentzian":
            pk_model = LorentzianModel(prefix=prefix)
        elif m_lower in ("pseudo_voigt", "pseudovoigt", "voigt"):
            pk_model = PseudoVoigtModel(prefix=prefix)
        else:
            raise ValueError(f"Unsupported model_type: {model_type}. Choose 'gaussian', 'lorentzian', or 'pseudo_voigt'.")

        # Estimate amplitude from local height
        near_idx = np.argmin(np.abs(x_arr - c_init))
        h_est = max(float(y_arr[near_idx]) - float(np.min(y_arr)), 1e-3)
        amp_est = h_est * initial_sigma * np.sqrt(2 * np.pi)

        pk_params = pk_model.make_params(
            center=float(c_init),
            sigma=float(initial_sigma),
            amplitude=float(amp_est),
        )

        # Apply physically meaningful parameter bounds
        pk_params[f"{prefix}center"].set(min=float(c_init - x_range * 0.1), max=float(c_init + x_range * 0.1))
        pk_params[f"{prefix}sigma"].set(min=1e-5)
        if max_fwhm is not None:
            pk_params[f"{prefix}sigma"].set(max=float(max_fwhm / 2.0))
        pk_params[f"{prefix}amplitude"].set(min=0.0)

        if m_lower in ("pseudo_voigt", "pseudovoigt", "voigt"):
            pk_params[f"{prefix}fraction"].set(value=0.5, min=0.0, max=1.0)

        params.update(pk_params)
        if composite_model is None:
            composite_model = pk_model
        else:
            composite_model = composite_model + pk_model

    # Execute fit
    result = composite_model.fit(y_arr, params, x=x_arr)

    # Compute goodness-of-fit
    ss_res = np.sum(result.residual ** 2)
    ss_tot = np.sum((y_arr - np.mean(y_arr)) ** 2)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

    # Extract individual component curves
    comp_dict = result.eval_components(x=x_arr)

    # Extract per-peak parameter summaries
    peaks_summary = []
    for idx in range(len(peak_centers)):
        prefix = f"p{idx+1}_"
        center_val = result.params[f"{prefix}center"].value
        center_err = result.params[f"{prefix}center"].stderr
        amp_val = result.params[f"{prefix}amplitude"].value
        amp_err = result.params[f"{prefix}amplitude"].stderr
        sigma_val = result.params[f"{prefix}sigma"].value
        sigma_err = result.params[f"{prefix}sigma"].stderr

        # FWHM depends on model
        m_lower = model_type.lower()
        if m_lower == "gaussian":
            fwhm_val = sigma_val * 2.35482
        elif m_lower == "lorentzian":
            fwhm_val = sigma_val * 2.0
        else:
            # Pseudo-Voigt approximation
            eta = result.params[f"{prefix}fraction"].value
            fwhm_g = sigma_val * 2.35482
            fwhm_l = sigma_val * 2.0
            fwhm_val = eta * fwhm_l + (1.0 - eta) * fwhm_g

        # Maximum height of this component
        curve = comp_dict.get(prefix[:-1], None)
        height_val = float(np.max(curve)) if curve is not None else float(amp_val / (sigma_val * np.sqrt(2 * np.pi)))

        peaks_summary.append({
            "index": idx + 1,
            "center": float(center_val),
            "center_stderr": float(center_err) if center_err is not None else None,
            "amplitude": float(amp_val),
            "amplitude_stderr": float(amp_err) if amp_err is not None else None,
            "fwhm": float(fwhm_val),
            "height": float(height_val),
            "fraction": float(result.params[f"{prefix}fraction"].value) if f"{prefix}fraction" in result.params else None,
        })

    return {
        "fit_y": np.asarray(result.best_fit),
        "residuals": np.asarray(result.residual),
        "r_squared": float(r2),
        "reduced_chi_sq": float(result.redchi),
        "peaks": peaks_summary,
        "components": comp_dict,
        "lmfit_result": result,
    }


def integrate_band(
    x: np.ndarray,
    y: np.ndarray,
    x_range: Tuple[float, float],
    baseline_subtracted: bool = False,
) -> float:
    """
    Numerically integrate a spectral band over [x_min, x_max] using Simpson's rule.

    Parameters
    ----------
    x : np.ndarray
        Abscissa array.
    y : np.ndarray
        Intensity/Absorbance array.
    x_range : tuple of (float, float)
        Integration bounds (order-independent).
    baseline_subtracted : bool, default False
        If False, subtracts a linear trapezoidal baseline connecting (x1, y1) and (x2, y2).

    Returns
    -------
    float
        Integrated area under the band.
    """
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)

    x1, x2 = sorted(x_range)
    mask = (x_arr >= x1) & (x_arr <= x2)
    if np.sum(mask) < 3:
        raise ValueError(f"Insufficient coordinate points within integration window [{x1}, {x2}].")

    x_sub = x_arr[mask]
    y_sub = y_arr[mask]

    # Ensure strictly increasing order for integration
    sort_idx = np.argsort(x_sub)
    x_sub = x_sub[sort_idx]
    y_sub = y_sub[sort_idx]

    if not baseline_subtracted:
        # Subtract straight baseline connecting boundary points
        base_line = y_sub[0] + (y_sub[-1] - y_sub[0]) * ((x_sub - x_sub[0]) / (x_sub[-1] - x_sub[0]))
        y_int = y_sub - base_line
    else:
        y_int = y_sub

    return float(simpson(y=y_int, x=x_sub))

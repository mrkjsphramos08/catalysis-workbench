"""
X-Ray Diffraction (XRD) Experimental Pattern Parsers.

Supports CSV, TXT, XY, DAT, and XLSX formats from common laboratory diffractometers
(e.g., Bruker, PANalytical, Rigaku, Malvern).
Extracts 2theta, Intensity, optional uncertainty, and radiation source metadata
(Cu, Co, Mo, Cr, or custom wavelength).
"""

from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np
import pandas as pd

from parsers.characterization.common import (
    CharacterizationData,
    CharacterizationParseError,
    compute_file_sha256,
)

# Standard emission line wavelengths in Angstroms
XRAY_WAVELENGTHS_ANGSTROM: Dict[str, float] = {
    "cu": 1.54060,      # Cu K-alpha (weighted average or K-alpha1: 1.54056)
    "cu_kalpha": 1.54060,
    "cu_kalpha1": 1.54056,
    "co": 1.78901,      # Co K-alpha
    "mo": 0.70930,      # Mo K-alpha
    "cr": 2.28970,      # Cr K-alpha
    "fe": 1.93604,      # Fe K-alpha
}


def parse_xrd_file(
    filepath: Union[str, Path],
    wavelength: Optional[Union[float, str]] = "cu",
    two_theta_col: Optional[Union[int, str]] = None,
    intensity_col: Optional[Union[int, str]] = None,
    error_col: Optional[Union[int, str]] = None,
    sample_id: Optional[str] = None,
    delimiter: Optional[str] = None,
) -> CharacterizationData:
    """
    Parse an exported XRD pattern into standardized CharacterizationData.

    Parameters
    ----------
    filepath : str or Path
        Target XRD file.
    wavelength : float or str, default 'cu'
        X-ray wavelength in Angstroms or source identifier ('cu', 'co', 'mo', 'cr', 'fe').
    two_theta_col : int or str, optional
        Column name or index for 2theta. If None, auto-detected.
    intensity_col : int or str, optional
        Column name or index for intensity. If None, auto-detected.
    error_col : int or str, optional
        Column name or index for standard uncertainty.
    sample_id : str, optional
        Sample name identifier. Defaults to filename stem.
    delimiter : str, optional
        Explicit delimiter. If None, auto-detected.

    Returns
    -------
    CharacterizationData
        Standardized XRD dataset with verified units and metadata.
    """
    p = Path(filepath)
    if not p.exists():
        raise FileNotFoundError(f"XRD file not found: {p.resolve()}")

    file_hash = compute_file_sha256(p)
    sid = sample_id or p.stem

    # Resolve wavelength in Angstroms
    if isinstance(wavelength, str):
        wl_key = wavelength.strip().lower()
        if wl_key in XRAY_WAVELENGTHS_ANGSTROM:
            wl_val = XRAY_WAVELENGTHS_ANGSTROM[wl_key]
            radiation_src = wavelength.upper()
        else:
            try:
                wl_val = float(wavelength)
                radiation_src = "Custom"
            except ValueError:
                raise ValueError(
                    f"Unknown X-ray radiation source '{wavelength}'. Available: {list(XRAY_WAVELENGTHS_ANGSTROM.keys())}"
                )
    else:
        wl_val = float(wavelength)
        radiation_src = "Custom"

    # Ingest file based on extension
    ext = p.suffix.lower()
    raw_metadata: Dict[str, Any] = {
        "radiation_source": radiation_src,
        "wavelength_angstrom": wl_val,
        "wavelength_nm": wl_val / 10.0,
    }

    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(p)
    else:
        # Determine header rows, comments, and delimiter for text/csv/dat/xy
        with open(p, "r", encoding="utf-8", errors="replace") as f:
            lines = [f.readline() for _ in range(50)]

        # Extract header metadata comments if present (e.g. #, //, *, [Header])
        skip_rows = 0
        for idx, line in enumerate(lines):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(("#", "//", ";", "*", "%")):
                skip_rows = idx + 1
                # Try to extract key-value pairs from comment
                if ":" in stripped:
                    parts = stripped.lstrip("#/;*% ").split(":", 1)
                    raw_metadata[parts[0].strip()] = parts[1].strip()
            elif re.search(r"[a-zA-Z]", stripped) and not re.search(r"^[0-9.eE+\-\s,;\t]+$", stripped):
                # Probably a header row with column labels
                skip_rows = idx
                break
            elif re.search(r"^[0-9.eE+\-\s,;\t]+$", stripped):
                # Pure numerical line
                skip_rows = idx
                break

        # Read into dataframe
        sep = delimiter if delimiter is not None else r"[\s,;\t]+"
        try:
            df = pd.read_csv(
                p,
                skiprows=skip_rows,
                sep=sep,
                engine="python" if delimiter is None else "c",
                comment="#",
            )
        except Exception:
            # Fallback: whitespace delimiter
            df = pd.read_csv(p, skiprows=skip_rows, delim_whitespace=True)

    # Clean empty columns
    df = df.dropna(how="all", axis=1)

    # Identify 2theta and Intensity columns
    col_2theta = None
    col_int = None
    col_err = None

    if two_theta_col is not None:
        col_2theta = df.columns[two_theta_col] if isinstance(two_theta_col, int) else two_theta_col
    if intensity_col is not None:
        col_int = df.columns[intensity_col] if isinstance(intensity_col, int) else intensity_col
    if error_col is not None:
        col_err = df.columns[error_col] if isinstance(error_col, int) else error_col

    # Auto-detection heuristics if columns not specified
    if col_2theta is None or col_int is None:
        cols = list(df.columns)
        for c in cols:
            c_str = str(c).strip().lower()
            if col_2theta is None and any(k in c_str for k in ["2theta", "2-theta", "2_theta", "two_theta", "theta", "angle", "2th", "pos"]):
                col_2theta = c
            elif col_int is None and any(k in c_str for k in ["intensity", "counts", "cps", "int", "y", "obs", "i"]):
                col_int = c
            elif col_err is None and any(k in c_str for k in ["err", "esd", "sigma", "std"]):
                col_err = c

        # Fallback to positional: column 0 = 2theta, column 1 = intensity
        if col_2theta is None and len(cols) >= 2:
            col_2theta = cols[0]
        if col_int is None and len(cols) >= 2:
            col_int = cols[1]

    if col_2theta is None or col_int is None:
        raise CharacterizationParseError(
            f"Failed to identify 2theta and Intensity columns in {p.name}. Found columns: {list(df.columns)}"
        )

    two_theta_arr = np.asarray(pd.to_numeric(df[col_2theta], errors="coerce"))
    int_arr = np.asarray(pd.to_numeric(df[col_int], errors="coerce"))

    # Mask valid coordinates
    mask = np.isfinite(two_theta_arr) & np.isfinite(int_arr)
    two_theta_clean = two_theta_arr[mask]
    int_clean = int_arr[mask]

    if len(two_theta_clean) < 10:
        raise CharacterizationParseError(f"XRD pattern in {p.name} contains fewer than 10 valid numerical points.")

    # Sort monotonically by 2theta
    sort_idx = np.argsort(two_theta_clean)
    two_theta_clean = two_theta_clean[sort_idx]
    int_clean = int_clean[sort_idx]

    err_clean = None
    if col_err is not None and col_err in df.columns:
        e_arr = np.asarray(pd.to_numeric(df[col_err], errors="coerce"))[mask][sort_idx]
        err_clean = e_arr

    # Populate scan metadata
    raw_metadata["scan_start_2theta"] = float(two_theta_clean[0])
    raw_metadata["scan_end_2theta"] = float(two_theta_clean[-1])
    diffs = np.diff(two_theta_clean)
    raw_metadata["step_size_2theta"] = float(np.median(diffs)) if len(diffs) > 0 else None

    return CharacterizationData(
        technique="XRD",
        sample_id=sid,
        source_file=p,
        file_sha256=file_hash,
        x_name="2Theta",
        x_unit="degrees",
        y_name="Intensity",
        y_unit="counts",
        x=two_theta_clean,
        y=int_clean,
        y_err=err_clean,
        dataframe=df,
        metadata=raw_metadata,
    )

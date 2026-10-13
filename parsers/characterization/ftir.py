"""
Fourier Transform Infrared Spectroscopy (FTIR) Data Parser.

Supports CSV, TXT, DAT, and XLSX exports from standard spectrophotometers
(e.g., Thermo Nicolet, Bruker, PerkinElmer, Shimadzu).
Handles Wavenumber, Transmittance (% or fraction), Absorbance, and Reflectance channels.
Converts Transmittance to Absorbance strictly via A = -log10(T_fraction)
and rejects nonpositive transmittance measurements.
"""

from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np
import pandas as pd

from parsers.characterization.common import (
    CharacterizationData,
    CharacterizationParseError,
    InvalidUnitsError,
    compute_file_sha256,
)


def parse_ftir_file(
    filepath: Union[str, Path],
    wavenumber_col: Optional[Union[int, str]] = None,
    signal_col: Optional[Union[int, str]] = None,
    signal_type: Optional[str] = None,  # 'transmittance', 'absorbance', 'reflectance'
    sample_id: Optional[str] = None,
    delimiter: Optional[str] = None,
) -> CharacterizationData:
    """
    Parse an exported FTIR spectrum into standardized CharacterizationData.

    Parameters
    ----------
    filepath : str or Path
        Target FTIR file path.
    wavenumber_col : int or str, optional
        Wavenumber column name or index. If None, auto-detected.
    signal_col : int or str, optional
        Signal column name or index. If None, auto-detected.
    signal_type : str, optional
        'transmittance', 'absorbance', or 'reflectance'. If None, inferred from headers.
    sample_id : str, optional
        Sample name identifier. Defaults to filename stem.
    delimiter : str, optional
        Explicit delimiter. If None, auto-detected.

    Returns
    -------
    CharacterizationData
        Standardized FTIR dataset.
    """
    p = Path(filepath)
    if not p.exists():
        raise FileNotFoundError(f"FTIR file not found: {p.resolve()}")

    file_hash = compute_file_sha256(p)
    sid = sample_id or p.stem

    ext = p.suffix.lower()
    raw_metadata: Dict[str, Any] = {}

    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(p)
    else:
        with open(p, "r", encoding="utf-8", errors="replace") as f:
            lines = [f.readline() for _ in range(50)]

        skip_rows = 0
        for idx, line in enumerate(lines):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(("#", "//", ";", "*", "%")):
                skip_rows = idx + 1
                if ":" in stripped:
                    parts = stripped.lstrip("#/;*% ").split(":", 1)
                    raw_metadata[parts[0].strip()] = parts[1].strip()
            elif re.search(r"[a-zA-Z]", stripped) and not re.search(r"^[0-9.eE+\-\s,;\t]+$", stripped):
                skip_rows = idx
                break
            elif re.search(r"^[0-9.eE+\-\s,;\t]+$", stripped):
                skip_rows = idx
                break

        sep = delimiter if delimiter is not None else r"[\s,;\t]+"
        try:
            df = pd.read_csv(p, skiprows=skip_rows, sep=sep, engine="python" if delimiter is None else "c", comment="#")
        except Exception:
            df = pd.read_csv(p, skiprows=skip_rows, delim_whitespace=True)

    df = df.dropna(how="all", axis=1)
    cols = list(df.columns)

    col_wn = None
    col_sig = None

    if wavenumber_col is not None:
        col_wn = df.columns[wavenumber_col] if isinstance(wavenumber_col, int) else wavenumber_col
    if signal_col is not None:
        col_sig = df.columns[signal_col] if isinstance(signal_col, int) else signal_col

    # Auto-detection heuristics
    if col_wn is None or col_sig is None:
        for c in cols:
            c_str = str(c).strip().lower()
            if col_wn is None and any(k in c_str for k in ["wavenumber", "cm-1", "cm^-1", "cm_1", "wn", "x", "energy"]):
                col_wn = c
            elif col_sig is None and any(k in c_str for k in ["transmittance", "%t", "t%", "absorbance", "abs", "a", "reflectance", "%r", "y", "intensity"]):
                col_sig = c

        if col_wn is None and len(cols) >= 2:
            col_wn = cols[0]
        if col_sig is None and len(cols) >= 2:
            col_sig = cols[1]

    if col_wn is None or col_sig is None:
        raise CharacterizationParseError(
            f"Failed to identify Wavenumber and Signal columns in {p.name}. Found: {cols}"
        )

    wn_arr = np.asarray(pd.to_numeric(df[col_wn], errors="coerce"))
    sig_arr = np.asarray(pd.to_numeric(df[col_sig], errors="coerce"))

    mask = np.isfinite(wn_arr) & np.isfinite(sig_arr)
    wn_clean = wn_arr[mask]
    sig_clean = sig_arr[mask]

    if len(wn_clean) < 10:
        raise CharacterizationParseError(f"FTIR spectrum in {p.name} contains fewer than 10 valid points.")

    # Determine signal type
    sig_type = signal_type.lower() if signal_type else None
    if sig_type is None:
        sig_header = str(col_sig).lower()
        if any(k in sig_header for k in ["transmittance", "%t", "t%"]):
            sig_type = "transmittance"
        elif any(k in sig_header for k in ["absorbance", "abs"]):
            sig_type = "absorbance"
        elif any(k in sig_header for k in ["reflectance", "%r"]):
            sig_type = "reflectance"
        else:
            # Infer from numerical value range
            # Transmittance is typically 0-100 or 0-1.0
            if np.max(sig_clean) > 2.5:
                sig_type = "transmittance"  # %T
            else:
                sig_type = "absorbance"

    # Sort descending or ascending? We store sorted ascending internally for interpolation/integration,
    # but publication plotting convention displays descending on x-axis.
    sort_idx = np.argsort(wn_clean)
    wn_clean = wn_clean[sort_idx]
    sig_clean = sig_clean[sort_idx]

    unit_str = "%" if sig_type == "transmittance" and np.max(sig_clean) > 1.5 else "fraction" if sig_type == "transmittance" else "a.u."
    raw_metadata["signal_type"] = sig_type
    raw_metadata["sampling_mode"] = raw_metadata.get("mode", "transmission/ATR")

    return CharacterizationData(
        technique="FTIR",
        sample_id=sid,
        source_file=p,
        file_sha256=file_hash,
        x_name="Wavenumber",
        x_unit="cm^-1",
        y_name=sig_type.capitalize(),
        y_unit=unit_str,
        x=wn_clean,
        y=sig_clean,
        dataframe=df,
        metadata=raw_metadata,
    )


def convert_transmittance_to_absorbance(
    data: CharacterizationData,
    clamp_nonpositive: bool = False,
) -> CharacterizationData:
    """
    Convert Transmittance to Absorbance strictly adhering to the Beer-Lambert relation:
    A = -log10(T_fraction)

    Parameters
    ----------
    data : CharacterizationData
        Input FTIR dataset. Must be transmittance.
    clamp_nonpositive : bool, default False
        If False, raises InvalidUnitsError upon encountering nonpositive transmittance values (T <= 0).
        If True, clamps nonpositive values to 1e-4 before log transform.

    Returns
    -------
    CharacterizationData
        New CharacterizationData instance with y converted to Absorbance.
    """
    sig_type = data.metadata.get("signal_type", "transmittance").lower()
    if sig_type != "transmittance":
        # Already absorbance
        return data.copy()

    y_t = data.y.copy()
    # If percentage, scale to fraction
    if np.max(y_t) > 1.5:
        t_frac = y_t / 100.0
    else:
        t_frac = y_t

    nonpos = t_frac <= 0.0
    if np.any(nonpos):
        if not clamp_nonpositive:
            min_val = float(np.min(t_frac))
            raise InvalidUnitsError(
                f"Invalid nonpositive transmittance encountered (min={min_val:.4f}). "
                "Beer-Lambert absorbance -log10(T) is mathematically undefined for T <= 0."
            )
        else:
            t_frac = np.clip(t_frac, 1e-4, 1.0)

    abs_arr = -np.log10(t_frac)

    out = data.copy()
    out.y = abs_arr
    out.y_name = "Absorbance"
    out.y_unit = "a.u."
    out.metadata["signal_type"] = "absorbance"
    out.add_processing_step("convert_to_absorbance", {"formula": "A = -log10(T)"})
    return out

"""
Raman Spectroscopy Data Parser.

Supports CSV, TXT, DAT, and XLSX exports from standard Raman spectrometers
(e.g., Renishaw, Horiba, WITec, Thermo Fisher).
Parses Raman shift (cm^-1), Intensity (counts / a.u.), and laser excitation metadata (532 nm, 633 nm, 785 nm).
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


def parse_raman_file(
    filepath: Union[str, Path],
    raman_shift_col: Optional[Union[int, str]] = None,
    intensity_col: Optional[Union[int, str]] = None,
    laser_wavelength_nm: Optional[float] = 532.0,
    sample_id: Optional[str] = None,
    delimiter: Optional[str] = None,
) -> CharacterizationData:
    """
    Parse an exported Raman spectrum into standardized CharacterizationData.

    Parameters
    ----------
    filepath : str or Path
        Target Raman file path.
    raman_shift_col : int or str, optional
        Raman shift column name or index. If None, auto-detected.
    intensity_col : int or str, optional
        Intensity column name or index. If None, auto-detected.
    laser_wavelength_nm : float, optional
        Excitation laser wavelength in nm (default 532 nm).
    sample_id : str, optional
        Sample name identifier. Defaults to filename stem.
    delimiter : str, optional
        Explicit delimiter. If None, auto-detected.

    Returns
    -------
    CharacterizationData
        Standardized Raman dataset.
    """
    p = Path(filepath)
    if not p.exists():
        raise FileNotFoundError(f"Raman file not found: {p.resolve()}")

    file_hash = compute_file_sha256(p)
    sid = sample_id or p.stem

    ext = p.suffix.lower()
    raw_metadata: Dict[str, Any] = {
        "laser_wavelength_nm": laser_wavelength_nm,
    }

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

    col_shift = None
    col_int = None

    if raman_shift_col is not None:
        col_shift = df.columns[raman_shift_col] if isinstance(raman_shift_col, int) else raman_shift_col
    if intensity_col is not None:
        col_int = df.columns[intensity_col] if isinstance(intensity_col, int) else intensity_col

    if col_shift is None or col_int is None:
        for c in cols:
            c_str = str(c).strip().lower()
            if col_shift is None and any(k in c_str for k in ["shift", "raman", "cm-1", "cm^-1", "wavenumber", "x"]):
                col_shift = c
            elif col_int is None and any(k in c_str for k in ["intensity", "counts", "cps", "ccd", "y", "i"]):
                col_int = c

        if col_shift is None and len(cols) >= 2:
            col_shift = cols[0]
        if col_int is None and len(cols) >= 2:
            col_int = cols[1]

    if col_shift is None or col_int is None:
        raise CharacterizationParseError(
            f"Failed to identify Raman shift and Intensity columns in {p.name}. Found: {cols}"
        )

    shift_arr = np.asarray(pd.to_numeric(df[col_shift], errors="coerce"))
    int_arr = np.asarray(pd.to_numeric(df[col_int], errors="coerce"))

    mask = np.isfinite(shift_arr) & np.isfinite(int_arr)
    shift_clean = shift_arr[mask]
    int_clean = int_arr[mask]

    if len(shift_clean) < 10:
        raise CharacterizationParseError(f"Raman spectrum in {p.name} contains fewer than 10 valid points.")

    sort_idx = np.argsort(shift_clean)
    shift_clean = shift_clean[sort_idx]
    int_clean = int_clean[sort_idx]

    return CharacterizationData(
        technique="Raman",
        sample_id=sid,
        source_file=p,
        file_sha256=file_hash,
        x_name="Raman Shift",
        x_unit="cm^-1",
        y_name="Intensity",
        y_unit="counts",
        x=shift_clean,
        y=int_clean,
        dataframe=df,
        metadata=raw_metadata,
    )

"""
Inductively Coupled Plasma (ICP-OES and ICP-MS) Analytical Data Parser.

Supports CSV, TXT, and XLSX exports from analytical spectrometers
(e.g., Agilent, PerkinElmer, Thermo Fisher).
Parses multi-element concentrations, raw emission/mass intensities,
calibration standards, blanks, replicate measurements, and dilution factors.
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


def parse_icp_file(
    filepath: Union[str, Path],
    element_col: Optional[Union[int, str]] = None,
    concentration_col: Optional[Union[int, str]] = None,
    sample_col: Optional[Union[int, str]] = None,
    intensity_col: Optional[Union[int, str]] = None,
    unit_col: Optional[Union[int, str]] = None,
    dilution_col: Optional[Union[int, str]] = None,
    default_unit: str = "mg/L",
    technique: str = "ICP-OES",
) -> CharacterizationData:
    """
    Parse an ICP-OES or ICP-MS dataset into standardized CharacterizationData.

    Parameters
    ----------
    filepath : str or Path
        Target ICP data file.
    element_col : int or str, optional
        Element/line column. Auto-detected if None.
    concentration_col : int or str, optional
        Measured concentration column. Auto-detected if None.
    sample_col : int or str, optional
        Sample name column. Auto-detected if None.
    intensity_col : int or str, optional
        Raw intensity / counts column. Auto-detected if None.
    unit_col : int or str, optional
        Unit column.
    dilution_col : int or str, optional
        Dilution factor column.
    default_unit : str, default 'mg/L'
        Fallback concentration unit if missing.
    technique : str, default 'ICP-OES'
        'ICP-OES' or 'ICP-MS'.

    Returns
    -------
    CharacterizationData
        Standardized ICP dataset containing clean dataframe and metadata.
    """
    p = Path(filepath)
    if not p.exists():
        raise FileNotFoundError(f"ICP file not found: {p.resolve()}")

    file_hash = compute_file_sha256(p)
    ext = p.suffix.lower()

    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(p)
    else:
        # Check header
        with open(p, "r", encoding="utf-8", errors="replace") as f:
            lines = [f.readline() for _ in range(50)]

        skip_rows = 0
        for idx, line in enumerate(lines):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(("#", "//", ";", "*")):
                skip_rows = idx + 1
            elif re.search(r"[a-zA-Z]", stripped) and not re.search(r"^[0-9.eE+\-\s,;\t]+$", stripped):
                skip_rows = idx
                break

        try:
            df = pd.read_csv(p, skiprows=skip_rows, sep=None, engine="python", comment="#")
        except Exception:
            df = pd.read_csv(p, skiprows=skip_rows)

    df = df.dropna(how="all", axis=1)
    cols = list(df.columns)

    # Column mappings
    c_elem = None
    c_conc = None
    c_samp = None
    c_int = None
    c_df = None

    if element_col is not None:
        c_elem = df.columns[element_col] if isinstance(element_col, int) else element_col
    if concentration_col is not None:
        c_conc = df.columns[concentration_col] if isinstance(concentration_col, int) else concentration_col
    if sample_col is not None:
        c_samp = df.columns[sample_col] if isinstance(sample_col, int) else sample_col
    if intensity_col is not None:
        c_int = df.columns[intensity_col] if isinstance(intensity_col, int) else intensity_col

    for c in cols:
        cl = str(c).strip().lower()
        if c_elem is None and any(k in cl for k in ["element", "analyte", "line", "isotope", "metal"]):
            c_elem = c
        elif c_conc is None and any(k in cl for k in ["conc", "concentration", "ppm", "ppb", "mg/l", "ug/l", "val"]):
            c_conc = c
        elif c_samp is None and any(k in cl for k in ["sample", "id", "name", "label", "solution"]):
            c_samp = c
        elif c_int is None and any(k in cl for k in ["intensity", "counts", "cps", "signal"]):
            c_int = c
        elif c_df is None and any(k in cl for k in ["dilution", "df", "factor"]):
            c_df = c

    # Fallback mappings
    if c_samp is None and len(cols) >= 1:
        c_samp = cols[0]
    if c_elem is None and len(cols) >= 2:
        c_elem = cols[1]
    if c_conc is None and len(cols) >= 3:
        c_conc = cols[2]

    if c_conc is None:
        raise CharacterizationParseError(f"Failed to identify concentration column in {p.name}. Found: {cols}")

    # Standardize column names
    clean_df = pd.DataFrame()
    clean_df["sample_id"] = df[c_samp].astype(str).str.strip() if c_samp else "Sample"
    clean_df["element"] = df[c_elem].astype(str).str.strip() if c_elem else "Element"
    clean_df["concentration"] = pd.to_numeric(df[c_conc], errors="coerce")

    if c_int and c_int in df.columns:
        clean_df["intensity"] = pd.to_numeric(df[c_int], errors="coerce")
    if c_df and c_df in df.columns:
        clean_df["dilution_factor"] = pd.to_numeric(df[c_df], errors="coerce").fillna(1.0)
    else:
        clean_df["dilution_factor"] = 1.0

    clean_df["unit"] = default_unit

    raw_metadata = {
        "elements_detected": list(clean_df["element"].unique()),
        "sample_count": len(clean_df["sample_id"].unique()),
        "default_unit": default_unit,
    }

    return CharacterizationData(
        technique=technique.upper(),
        sample_id=p.stem,
        source_file=p,
        file_sha256=file_hash,
        x_name="Element",
        x_unit="",
        y_name="Concentration",
        y_unit=default_unit,
        dataframe=clean_df,
        metadata=raw_metadata,
    )

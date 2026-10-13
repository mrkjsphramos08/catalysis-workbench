"""
Atomic Absorption Spectroscopy (AAS) Data Parser.

Supports:
- Flame Atomic Absorption Spectroscopy (FAAS)
- Graphite Furnace Atomic Absorption Spectroscopy (GFAAS / ETAAS)
- Architecture hooks for Hydride Generation AAS (HG-AAS) and Cold Vapor AAS (CV-AAS)

Reads exported CSV, TXT, TSV, and XLSX files from major instrument manufacturers
(PerkinElmer, Shimadzu, Thermo Scientific, Agilent, Analytik Jena).
Captures analyte elements, wavelengths, steady-state absorbance, integrated absorbance,
replicate readings, blanks, dilution factors, and preparation metadata.
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


def parse_aas_file(
    filepath: Union[str, Path],
    element_col: Optional[Union[int, str]] = None,
    absorbance_col: Optional[Union[int, str]] = None,
    integrated_abs_col: Optional[Union[int, str]] = None,
    concentration_col: Optional[Union[int, str]] = None,
    sample_col: Optional[Union[int, str]] = None,
    sample_type_col: Optional[Union[int, str]] = None,
    replicate_col: Optional[Union[int, str]] = None,
    unit_col: Optional[Union[int, str]] = None,
    dilution_col: Optional[Union[int, str]] = None,
    wavelength_col: Optional[Union[int, str]] = None,
    default_unit: str = "mg/L",
    technique: str = "FAAS",
    instrument_model: Optional[str] = None,
    background_correction: Optional[str] = None,
) -> CharacterizationData:
    """
    Parse an AAS experimental or calibration export into standardized CharacterizationData.

    Parameters
    ----------
    filepath : str or Path
        Target AAS data file (CSV, TXT, TSV, or XLSX).
    element_col : int or str, optional
        Element/analyte column. Auto-detected if None.
    absorbance_col : int or str, optional
        Steady-state absorbance column (AU). Auto-detected if None.
    integrated_abs_col : int or str, optional
        Integrated absorbance column (A*s) for GFAAS. Auto-detected if None.
    concentration_col : int or str, optional
        Instrument-reported concentration column. Auto-detected if None.
    sample_col : int or str, optional
        Sample name / identifier column. Auto-detected if None.
    sample_type_col : int or str, optional
        Sample type column ('Standard', 'Sample', 'Blank', 'QC', etc.).
    replicate_col : int or str, optional
        Replicate index column.
    unit_col : int or str, optional
        Concentration unit column.
    dilution_col : int or str, optional
        Dilution factor column.
    wavelength_col : int or str, optional
        Analytical wavelength (nm) column.
    default_unit : str, default 'mg/L'
        Fallback concentration unit if missing.
    technique : str, default 'FAAS'
        AAS measurement technique: 'FAAS', 'GFAAS', 'HG-AAS', 'CV-AAS'.
    instrument_model : str, optional
        Spectrometer instrument model.
    background_correction : str, optional
        Background correction status ('D2', 'Zeeman', 'Smith-Hieftje', 'None').

    Returns
    -------
    CharacterizationData
        Standardized AAS dataset containing clean dataframe and metadata.
    """
    p = Path(filepath)
    if not p.exists():
        raise FileNotFoundError(f"AAS file not found: {p.resolve()}")

    file_hash = compute_file_sha256(p)
    ext = p.suffix.lower()

    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(p)
    else:
        # Inspect header lines
        with open(p, "r", encoding="utf-8", errors="replace") as f:
            lines = [f.readline() for _ in range(50)]

        skip_rows = 0
        detected_sep = None
        for idx, line in enumerate(lines):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(("#", "//", ";", "*")):
                skip_rows = idx + 1
            elif re.search(r"[a-zA-Z]", stripped) and not re.search(r"^[0-9.eE+\-\s,;\t]+$", stripped):
                skip_rows = idx
                if "\t" in stripped:
                    detected_sep = "\t"
                elif ";" in stripped and "," not in stripped:
                    detected_sep = ";"
                elif "," in stripped:
                    detected_sep = ","
                break

        try:
            df = pd.read_csv(p, skiprows=skip_rows, sep=detected_sep, engine="python", comment="#")
        except Exception:
            df = pd.read_csv(p, skiprows=skip_rows)

    df = df.dropna(how="all", axis=1)
    cols = list(df.columns)

    # Column mappings
    c_elem = None
    c_abs = None
    c_int_abs = None
    c_conc = None
    c_samp = None
    c_stype = None
    c_rep = None
    c_df = None
    c_wl = None
    c_unit = None

    if element_col is not None:
        c_elem = df.columns[element_col] if isinstance(element_col, int) else element_col
    if absorbance_col is not None:
        c_abs = df.columns[absorbance_col] if isinstance(absorbance_col, int) else absorbance_col
    if integrated_abs_col is not None:
        c_int_abs = df.columns[integrated_abs_col] if isinstance(integrated_abs_col, int) else integrated_abs_col
    if concentration_col is not None:
        c_conc = df.columns[concentration_col] if isinstance(concentration_col, int) else concentration_col
    if sample_col is not None:
        c_samp = df.columns[sample_col] if isinstance(sample_col, int) else sample_col
    if sample_type_col is not None:
        c_stype = df.columns[sample_type_col] if isinstance(sample_type_col, int) else sample_type_col
    if replicate_col is not None:
        c_rep = df.columns[replicate_col] if isinstance(replicate_col, int) else replicate_col
    if unit_col is not None:
        c_unit = df.columns[unit_col] if isinstance(unit_col, int) else unit_col
    if dilution_col is not None:
        c_df = df.columns[dilution_col] if isinstance(dilution_col, int) else dilution_col
    if wavelength_col is not None:
        c_wl = df.columns[wavelength_col] if isinstance(wavelength_col, int) else wavelength_col

    for c in cols:
        cl = str(c).strip().lower()
        if c_elem is None and any(k in cl for k in ["element", "analyte", "metal"]):
            c_elem = c
        elif c_int_abs is None and any(k in cl for k in ["integrated_abs", "int_abs", "peak_area", "area_abs", "a*s", "a.s", "abs_sec"]):
            c_int_abs = c
        elif c_abs is None and any(k in cl for k in ["absorbance", "abs", "mean_abs", "sig", "signal", "peak_height"]) and "int" not in cl and "area" not in cl:
            c_abs = c
        elif c_conc is None and any(k in cl for k in ["conc", "concentration", "ppm", "ppb", "mg/l", "ug/l", "calc_conc"]):
            c_conc = c
        elif c_samp is None and any(k in cl for k in ["sample", "id", "name", "label", "solution"]):
            c_samp = c
        elif c_stype is None and any(k in cl for k in ["type", "sample_type", "role"]):
            c_stype = c
        elif c_rep is None and any(k in cl for k in ["rep", "replicate", "run"]):
            c_rep = c
        elif c_df is None and any(k in cl for k in ["dilution", "df", "factor"]):
            c_df = c
        elif c_wl is None and any(k in cl for k in ["wavelength", "lambda", "nm"]):
            c_wl = c
        elif c_unit is None and any(k in cl for k in ["unit", "conc_unit"]):
            c_unit = c

    # Fallback heuristics
    if c_samp is None and len(cols) >= 1:
        c_samp = cols[0]

    # Verify that at least one response column (absorbance, integrated absorbance, or concentration) exists
    if c_abs is None and c_int_abs is None and c_conc is None:
        raise CharacterizationParseError(f"Failed to identify signal or concentration columns in AAS file '{p.name}'. Found columns: {cols}")

    # Standardize DataFrame
    clean_df = pd.DataFrame()
    clean_df["sample_id"] = df[c_samp].astype(str).str.strip() if c_samp else "Sample"

    if c_elem and c_elem in df.columns:
        clean_df["element"] = df[c_elem].astype(str).str.strip()
    else:
        clean_df["element"] = "Analyte"

    if c_stype and c_stype in df.columns:
        clean_df["sample_type"] = df[c_stype].astype(str).str.strip()
    else:
        # Infer sample type from sample_id name
        clean_df["sample_type"] = clean_df["sample_id"].apply(_infer_sample_type)

    if c_rep and c_rep in df.columns:
        clean_df["replicate"] = pd.to_numeric(df[c_rep], errors="coerce").fillna(1).astype(int)
    else:
        clean_df["replicate"] = 1

    # Absorbance response
    if c_abs and c_abs in df.columns:
        clean_df["absorbance"] = pd.to_numeric(df[c_abs], errors="coerce")
    else:
        clean_df["absorbance"] = np.nan

    # Integrated absorbance for GFAAS
    if c_int_abs and c_int_abs in df.columns:
        clean_df["integrated_absorbance"] = pd.to_numeric(df[c_int_abs], errors="coerce")
    else:
        clean_df["integrated_absorbance"] = np.nan

    # Instrument reported concentration
    if c_conc and c_conc in df.columns:
        clean_df["concentration_reported"] = pd.to_numeric(df[c_conc], errors="coerce")
    else:
        clean_df["concentration_reported"] = np.nan

    # Dilution factor
    if c_df and c_df in df.columns:
        clean_df["dilution_factor"] = pd.to_numeric(df[c_df], errors="coerce").fillna(1.0)
    else:
        clean_df["dilution_factor"] = 1.0

    # Wavelength
    if c_wl and c_wl in df.columns:
        clean_df["wavelength_nm"] = pd.to_numeric(df[c_wl], errors="coerce")
    else:
        clean_df["wavelength_nm"] = np.nan

    # Units
    if c_unit and c_unit in df.columns:
        clean_df["unit"] = df[c_unit].astype(str).str.strip()
    else:
        clean_df["unit"] = default_unit

    # Check primary analytical signal based on technique
    tech_upper = technique.upper()
    if tech_upper == "GFAAS" and not clean_df["integrated_absorbance"].isna().all():
        primary_y_name = "Integrated Absorbance"
        primary_y_unit = "A*s"
        primary_signal = clean_df["integrated_absorbance"].values
    elif not clean_df["absorbance"].isna().all():
        primary_y_name = "Absorbance"
        primary_y_unit = "AU"
        primary_signal = clean_df["absorbance"].values
    else:
        primary_y_name = "Concentration"
        primary_y_unit = default_unit
        primary_signal = clean_df["concentration_reported"].values

    raw_metadata = {
        "technique": tech_upper,
        "instrument_model": instrument_model or "Unknown AAS Spectrometer",
        "background_correction": background_correction or "Unspecified",
        "elements_detected": list(clean_df["element"].unique()),
        "sample_count": len(clean_df["sample_id"].unique()),
        "total_readings": len(clean_df),
        "default_unit": default_unit,
        "wavelengths": [float(w) for w in clean_df["wavelength_nm"].dropna().unique()] if not clean_df["wavelength_nm"].isna().all() else [],
    }

    return CharacterizationData(
        technique=tech_upper,
        sample_id=p.stem,
        source_file=p,
        file_sha256=file_hash,
        x_name="Concentration",
        x_unit=default_unit,
        y_name=primary_y_name,
        y_unit=primary_y_unit,
        y=primary_signal,
        dataframe=clean_df,
        metadata=raw_metadata,
        parser_version="1.0.0-aas",
    )


def _infer_sample_type(sample_id: str) -> str:
    """Heuristically categorize sample role from its ID string."""
    sid = str(sample_id).strip().upper()
    if sid.startswith("STD") or "CAL" in sid or "STANDARD" in sid:
        return "Standard"
    elif "BLANK" in sid or sid.startswith("BLK"):
        return "Blank"
    elif "ICV" in sid:
        return "ICV"
    elif "CCV" in sid:
        return "CCV"
    elif "SPIKE" in sid or "MS" in sid:
        return "Spike"
    elif "CRM" in sid or "REF" in sid:
        return "ReferenceMaterial"
    return "Sample"

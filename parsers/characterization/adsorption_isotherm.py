"""
Gas Adsorption Isotherm Data Parser.

Supports CSV, TXT, and XLSX files from standard volumetric adsorption analyzers
(e.g., Micromeritics 3Flex / ASAP, Quantachrome Autosorb, BELSORP).
Extracts Relative Pressure (p/p0), Adsorbed Quantity (cm^3/g STP or mmol/g),
Adsorption vs Desorption branches, adsorbate gas (N2, Ar, CO2), and temperature (77 K).
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


def parse_isotherm_file(
    filepath: Union[str, Path],
    pressure_col: Optional[Union[int, str]] = None,
    loading_col: Optional[Union[int, str]] = None,
    branch_col: Optional[Union[int, str]] = None,
    adsorbate: str = "N2",
    temperature_K: float = 77.35,
    sample_id: Optional[str] = None,
    delimiter: Optional[str] = None,
) -> CharacterizationData:
    """
    Parse a gas physisorption isotherm dataset into standardized CharacterizationData.

    Parameters
    ----------
    filepath : str or Path
        Target isotherm data file.
    pressure_col : int or str, optional
        Relative pressure (p/p0) column. Auto-detected if None.
    loading_col : int or str, optional
        Adsorbed amount (quantity adsorbed) column. Auto-detected if None.
    branch_col : int or str, optional
        Explicit adsorption/desorption branch indicator column.
    adsorbate : str, default 'N2'
        Adsorbate gas molecule ('N2', 'Ar', 'CO2').
    temperature_K : float, default 77.35
        Measurement temperature in Kelvin.
    sample_id : str, optional
        Sample name identifier. Defaults to filename stem.
    delimiter : str, optional
        Explicit delimiter. If None, auto-detected.

    Returns
    -------
    CharacterizationData
        Standardized physisorption dataset containing adsorption and desorption series.
    """
    p = Path(filepath)
    if not p.exists():
        raise FileNotFoundError(f"Isotherm file not found: {p.resolve()}")

    file_hash = compute_file_sha256(p)
    sid = sample_id or p.stem
    ext = p.suffix.lower()

    raw_metadata: Dict[str, Any] = {
        "adsorbate": adsorbate.upper(),
        "temperature_K": temperature_K,
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
            if stripped.startswith(("#", "//", ";", "*")):
                skip_rows = idx + 1
            elif re.search(r"[a-zA-Z]", stripped) and not re.search(r"^[0-9.eE+\-\s,;\t]+$", stripped):
                skip_rows = idx
                break

        sep = delimiter if delimiter is not None else r"[\s,;\t]+"
        try:
            df = pd.read_csv(p, skiprows=skip_rows, sep=sep, engine="python" if delimiter is None else "c", comment="#")
        except Exception:
            df = pd.read_csv(p, skiprows=skip_rows, delim_whitespace=True)

    df = df.dropna(how="all", axis=1)
    cols = list(df.columns)

    c_p = None
    c_q = None
    c_br = None

    if pressure_col is not None:
        c_p = df.columns[pressure_col] if isinstance(pressure_col, int) else pressure_col
    if loading_col is not None:
        c_q = df.columns[loading_col] if isinstance(loading_col, int) else loading_col
    if branch_col is not None:
        c_br = df.columns[branch_col] if isinstance(branch_col, int) else branch_col

    for c in cols:
        cl = str(c).strip().lower()
        if c_p is None and any(k in cl for k in ["p/p0", "p/po", "relative_pressure", "rel_pressure", "p0", "pe"]):
            c_p = c
        elif c_q is None and any(k in cl for k in ["quantity", "volume", "adsorbed", "loading", "q", "v_ads", "cm3/g", "mmol/g"]):
            c_q = c
        elif c_br is None and any(k in cl for k in ["branch", "type", "adsorp", "desorp"]):
            c_br = c

    if c_p is None and len(cols) >= 2:
        c_p = cols[0]
    if c_q is None and len(cols) >= 2:
        c_q = cols[1]

    if c_p is None or c_q is None:
        raise CharacterizationParseError(f"Failed to identify relative pressure and loading columns in {p.name}. Found: {cols}")

    p_arr = np.asarray(pd.to_numeric(df[c_p], errors="coerce"))
    q_arr = np.asarray(pd.to_numeric(df[c_q], errors="coerce"))

    mask = np.isfinite(p_arr) & np.isfinite(q_arr) & (p_arr >= 0.0) & (p_arr <= 1.05)
    p_clean = p_arr[mask]
    q_clean = q_arr[mask]

    # Detect adsorption vs desorption branches
    # If explicit branch column is present:
    if c_br and c_br in df.columns:
        branch_series = df[c_br][mask].astype(str).str.lower()
        is_ads = branch_series.str.contains("ads").values
    else:
        # Heuristic: Find maximum pressure index.
        # Points before peak are adsorption; points after peak with descending pressure are desorption.
        peak_idx = int(np.argmax(p_clean))
        is_ads = np.zeros(len(p_clean), dtype=bool)
        is_ads[:peak_idx + 1] = True

    # Assemble clean dataframe
    clean_df = pd.DataFrame({
        "p_p0": p_clean,
        "quantity_adsorbed": q_clean,
        "branch": np.where(is_ads, "adsorption", "desorption"),
    })

    raw_metadata["adsorption_points"] = int(np.sum(is_ads))
    raw_metadata["desorption_points"] = int(np.sum(~is_ads))
    raw_metadata["max_p_p0"] = float(np.max(p_clean))

    return CharacterizationData(
        technique="BET",
        sample_id=sid,
        source_file=p,
        file_sha256=file_hash,
        x_name="Relative Pressure (p/p0)",
        x_unit="p/p0",
        y_name="Quantity Adsorbed",
        y_unit="cm^3/g STP",
        x=p_clean,
        y=q_clean,
        dataframe=clean_df,
        metadata=raw_metadata,
    )

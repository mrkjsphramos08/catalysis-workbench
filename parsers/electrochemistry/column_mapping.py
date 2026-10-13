"""
Column Recognition, Fuzzy Aliasing, and Unit Normalization for Electrochemistry Instruments.

Handles customized signal arrangements from Metrohm Autolab (NOVA)
and BioLogic (EC-Lab) exports, accounting for variations in column headers,
units (A vs mA vs uA), and instrument naming conventions.
"""

import re
from typing import Dict, List, Optional, Sequence, Tuple, Union
import pandas as pd


# Canonical physical channel names in StandardizedElectrochemicalData
CANONICAL_CHANNELS = {
    "potential": "potential_V",
    "current": "current_A",
    "time": "time_s",
    "cycle": "cycle_number",
    "segment": "segment_id",
    "frequency": "frequency_Hz",
    "z_real": "impedance_real_ohm",
    "z_imag": "impedance_imag_ohm",
    "z_mag": "impedance_mag_ohm",
    "phase": "phase_angle_deg",
}

# Regex patterns for matching Metrohm NOVA and BioLogic EC-Lab column headers
HEADER_PATTERNS = {
    "potential": [
        r"^(?:we\(?1?\)?\.?)?\s*potential(?:\s*applied)?\s*(?:\(v\)|/v)?$",
        r"^ewe(?:\s*-\s*ece)?(?:/v|\s*\(v\))?$",
        r"^<ewe>(?:/v|\s*\(v\))?$",
        r"^potential(?:/v|\s*\(v\))?$",
        r"^voltage(?:/v|\s*\(v\))?$",
        r"^e\s*(?:\(v\)|/v)$",
        r"^u\s*(?:\(v\)|/v)$",
    ],
    "current": [
        r"^(?:we\(?1?\)?\.?)?\s*current\s*(?:\([a-zµuA-Z]+\)|/[a-zµuA-Z]+)?$",
        r"^<i>?(?:/[a-zµuA-Z]+|\s*\([a-zµuA-Z]+\))?$",
        r"^<i\s*-\s*ic>(?:/[a-zµuA-Z]+)?$",
        r"^i(?:/a|/ma|/ua|/na|/µa|\s*\(a\)|\s*\(ma\)|\s*\(ua\))?$",
    ],
    "time": [
        r"^(?:elapsed\s*)?time\s*(?:\(s\)|\(min\)|\(h\)|/s|/min|/h)?$",
        r"^t\s*(?:\(s\)|/s)$",
    ],
    "cycle": [
        r"^(?:cycle(?:\s*number)?|scan|index|segment|loop)$",
        r"^ox/red$",
    ],
    "segment": [
        r"^(?:segment(?:\s*id)?|step|sweep)$",
    ],
    "frequency": [
        r"^(?:ac\s*)?frequency\s*(?:\(hz\)|/hz)?$",
        r"^freq(?:/hz|\s*\(hz\))?$",
        r"^f\s*(?:\(hz\)|/hz)$",
    ],
    "z_real": [
        r"^(?:re\(z\)|z'|z_real|real\s*z)\s*(?:\([a-zA-Z\sΩ]+\)|/[a-zA-Z\sΩ]+)?$",
        r"^z'\s*\((?:omega|ohm|Ω)\)$",
    ],
    "z_imag": [
        r"^(?:-?im\(z\)|-?z''|-?z\"|z_imag|imag\s*z)\s*(?:\([a-zA-Z\sΩ]+\)|/[a-zA-Z\sΩ]+)?$",
        r"^-?z\"\s*\((?:omega|ohm|Ω)\)$",
    ],
    "z_mag": [
        r"^(?:\|z\||abs\(z\)|mag\(z\)|z_mag)\s*(?:\([a-zA-Z\sΩ]+\)|/[a-zA-Z\sΩ]+)?$",
    ],
    "phase": [
        r"^(?:phase(?:\s*angle)?(?:\(z\))?|-?phase)\s*(?:\(deg\)|/deg|deg)?$",
    ],
}


def _clean_header(col_name: str) -> str:
    """Normalize column header for regex comparison."""
    cleaned = col_name.strip().lower()
    cleaned = re.sub(r"\s+", " ", cleaned)
    return cleaned


def _detect_current_multiplier(col_name: str) -> float:
    """Extract multiplier to convert current into Amperes (A)."""
    clean = col_name.strip().lower()
    if any(u in clean for u in ("/ma", "(ma)", "ma", "milliamp")):
        return 1.0e-3
    elif any(u in clean for u in ("/ua", "(ua)", "ua", "/µa", "(µa)", "µa", "microamp")):
        return 1.0e-6
    elif any(u in clean for u in ("/na", "(na)", "na", "nanoamp")):
        return 1.0e-9
    elif any(u in clean for u in ("/ka", "(ka)", "ka")):
        return 1.0e3
    # Default is Amperes (A)
    return 1.0


def _detect_potential_multiplier(col_name: str) -> float:
    """Extract multiplier to convert potential into Volts (V)."""
    clean = col_name.strip().lower()
    if any(u in clean for u in ("/mv", "(mv)", "millivolt")):
        return 1.0e-3
    return 1.0


def _detect_time_multiplier(col_name: str) -> float:
    """Extract multiplier to convert time into seconds (s)."""
    clean = col_name.strip().lower()
    if any(u in clean for u in ("/min", "(min)", "minutes")):
        return 60.0
    elif any(u in clean for u in ("/h", "(h)", "hours", "hrs")):
        return 3600.0
    elif any(u in clean for u in ("/ms", "(ms)", "millisec")):
        return 1.0e-3
    return 1.0


def _detect_impedance_multiplier(col_name: str) -> float:
    """Extract multiplier to convert impedance into Ohms."""
    clean = col_name.strip().lower()
    if any(u in clean for u in ("/kohm", "(kohm)", "komega", "kΩ")):
        return 1.0e3
    elif any(u in clean for u in ("/mohm", "(mohm)", "momega", "mΩ")):
        return 1.0e6
    return 1.0


def identify_columns(
    columns: Sequence[str],
    custom_mapping: Optional[Dict[str, str]] = None,
) -> Tuple[Dict[str, str], Dict[str, float]]:
    """
    Automatically identify column roles and unit scaling factors.

    Parameters
    ----------
    columns : sequence of str
        Original column headers from raw file.
    custom_mapping : dict, optional
        User-provided override mapping {original_header: canonical_key}.

    Returns
    -------
    Tuple[Dict[str, str], Dict[str, float]]
        - mapping: {original_header: canonical_role}
        - scale_factors: {original_header: multiplier_to_si}
    """
    mapping: Dict[str, str] = {}
    scales: Dict[str, float] = {}
    used_roles = set()

    # Apply user overrides first if provided
    if custom_mapping:
        for orig_col, role in custom_mapping.items():
            if orig_col in columns and role in CANONICAL_CHANNELS:
                mapping[orig_col] = CANONICAL_CHANNELS[role]
                used_roles.add(role)
                if role == "current":
                    scales[orig_col] = _detect_current_multiplier(orig_col)
                elif role == "potential":
                    scales[orig_col] = _detect_potential_multiplier(orig_col)
                elif role == "time":
                    scales[orig_col] = _detect_time_multiplier(orig_col)
                elif role in ("z_real", "z_imag", "z_mag"):
                    scales[orig_col] = _detect_impedance_multiplier(orig_col)
                else:
                    scales[orig_col] = 1.0

    # Auto-detect remaining columns
    for col in columns:
        if col in mapping:
            continue
        cleaned = _clean_header(col)

        matched_role = None
        for role, patterns in HEADER_PATTERNS.items():
            if role in used_roles and role not in ("cycle", "segment"):
                continue  # Prefer first primary match
            for pat in patterns:
                if re.match(pat, cleaned):
                    matched_role = role
                    break
            if matched_role:
                break

        if matched_role:
            canonical_name = CANONICAL_CHANNELS[matched_role]
            mapping[col] = canonical_name
            used_roles.add(matched_role)

            if matched_role == "current":
                scales[col] = _detect_current_multiplier(col)
            elif matched_role == "potential":
                scales[col] = _detect_potential_multiplier(col)
            elif matched_role == "time":
                scales[col] = _detect_time_multiplier(col)
            elif matched_role in ("z_real", "z_imag", "z_mag"):
                scales[col] = _detect_impedance_multiplier(col)
            else:
                scales[col] = 1.0

    return mapping, scales


def standardize_dataframe(
    df: pd.DataFrame,
    mapping: Dict[str, str],
    scale_factors: Dict[str, float],
) -> pd.DataFrame:
    """
    Transform raw dataframe into standardized physical SI units.

    Parameters
    ----------
    df : pd.DataFrame
        Raw dataframe.
    mapping : dict
        {raw_column: canonical_column}.
    scale_factors : dict
        {raw_column: scale_multiplier}.

    Returns
    -------
    pd.DataFrame
        Cleaned dataframe with canonical column names and SI units.
    """
    clean_df = pd.DataFrame()
    for raw_col, canon_col in mapping.items():
        if raw_col in df.columns:
            series = pd.to_numeric(df[raw_col], errors="coerce")
            factor = scale_factors.get(raw_col, 1.0)
            clean_df[canon_col] = series * factor

    # Handle negative imaginary impedance if exported as -Im(Z) vs Im(Z)
    if "impedance_imag_ohm" in clean_df.columns:
        # Check original column name
        orig_imag_cols = [c for c, canon in mapping.items() if canon == "impedance_imag_ohm"]
        if orig_imag_cols:
            orig_col = orig_imag_cols[0].lower()
            # If the column was literally labeled "-Im(Z)" or "-Z''", it is already negative
            # In scientific EIS, Z = Z' + j*Z''. Usually capacitive Z'' is negative.
            # BioLogic often exports -Im(Z)/Ohm as positive numbers.
            # We want standard physical Z'' (which is negative for capacitors).
            if "-im" in orig_col or "-z\"" in orig_col or "-z'" in orig_col:
                # If exported values are positive, negate them to represent true physical Z''
                if (clean_df["impedance_imag_ohm"] > 0).mean() > 0.5:
                    clean_df["impedance_imag_ohm"] = -clean_df["impedance_imag_ohm"]

    return clean_df

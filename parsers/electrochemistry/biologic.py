"""
BioLogic EC-Lab Data Parser.

Parses .mpt (ASCII), .mps (settings), and .mpr (binary via galvani) files produced by BioLogic EC-Lab.
Supports all standard techniques: CV, LSV, PEIS/GEIS, CA, CP, OCV.
Extracts technique headers, scan rates, cycle indices, and formats into StandardizedElectrochemicalData.
"""

from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from parsers.electrochemistry.column_mapping import (
    identify_columns,
    standardize_dataframe,
)
from parsers.electrochemistry.schema import (
    ElectrochemistryParseError,
    MissingChannelError,
    StandardizedElectrochemicalData,
    compute_file_sha256,
    resolve_experimental_conditions,
)

try:
    import galvani
    GALVANI_AVAILABLE = True
except ImportError:
    GALVANI_AVAILABLE = False


class BioLogicParseError(ElectrochemistryParseError):
    """Raised when BioLogic EC-Lab file parsing encounters errors."""
    pass


class BioLogicBinaryMprError(BioLogicParseError):
    """Raised when BioLogic binary .mpr file fails to parse."""
    pass


def parse_biologic_mpt_header(lines: Sequence[str]) -> Tuple[Dict[str, Any], int]:
    """
    Parse header lines of a BioLogic .mpt ASCII export.
    BioLogic files typically begin with:
        EC-Lab ASCII FILE
        Nb header lines : XX
    """
    metadata: Dict[str, Any] = {}
    n_header_lines = 0

    if not lines:
        return metadata, 0

    # Look for "Nb header lines : XX"
    for idx, line in enumerate(lines[:10]):
        m = re.search(r"nb\s+header\s+lines\s*:\s*(\d+)", line, re.IGNORECASE)
        if m:
            n_header_lines = int(m.group(1))
            break

    # If not explicitly stated, search for data header row containing "Ewe" or "time"
    if n_header_lines == 0:
        for idx, line in enumerate(lines[:100]):
            lower = line.lower()
            if any(term in lower for term in ("ewe/v", "time/s", "freq/hz", "i/ma", "i/a", "re(z)", "control/v", "re(z)/ohm")):
                n_header_lines = idx + 1
                break
    else:
        # Verify if line (n_header_lines - 1) is indeed the column row, or adjust if off-by-one
        check_idx = n_header_lines - 1
        if 0 <= check_idx < len(lines):
            curr_lower = lines[check_idx].lower()
            if not any(term in curr_lower for term in ("ewe", "time", "freq", "i/", "re(z)", "control")):
                # Search nearby +/- 3 lines for the actual column row
                for offset in (-1, 1, -2, 2):
                    adj_idx = check_idx + offset
                    if 0 <= adj_idx < len(lines):
                        adj_lower = lines[adj_idx].lower()
                        if any(term in adj_lower for term in ("ewe", "time", "freq", "i/", "re(z)", "control")):
                            n_header_lines = adj_idx + 1
                            break

    # Extract metadata from header lines
    for line in lines[:n_header_lines]:
        line_s = line.strip()
        if ":" in line_s:
            parts = line_s.split(":", 1)
            k = parts[0].strip()
            v = parts[1].strip()
            metadata[k] = v

    return metadata, n_header_lines


def parse_biologic_mps(filepath: Union[str, Path]) -> Dict[str, Any]:
    """
    Parse a BioLogic .mps settings / technique parameters file.
    """
    p = Path(filepath)
    if not p.is_file():
        raise FileNotFoundError(f"MPS settings file not found: {p.resolve()}")

    settings: Dict[str, Any] = {}
    with open(p, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line_s = line.strip()
            if ":" in line_s:
                parts = line_s.split(":", 1)
                settings[parts[0].strip()] = parts[1].strip()
    return settings


def parse_biologic_mpt(
    filepath: Union[str, Path],
    technique: Optional[str] = None,
    custom_mapping: Optional[Dict[str, str]] = None,
    electrode_area_cm2: Optional[float] = None,
    reference_electrode: Optional[str] = None,
    ph: Optional[float] = None,
    catalyst_loading_mg_cm2: Optional[float] = None,
    ru_ohm: Optional[float] = None,
    encoding: str = "latin1",
    temperature_k: Optional[float] = None,
) -> StandardizedElectrochemicalData:
    """
    Parse BioLogic EC-Lab .mpt ASCII data file.

    Parameters
    ----------
    filepath : str or Path
        Path to .mpt file.
    technique : str, optional
        Electrochemical technique ('CV', 'LSV', 'EIS', 'CA', 'CP').
    custom_mapping : dict, optional
        User override for signal mapping.
    electrode_area_cm2 : float, optional
        Geometric surface area in cm^2.
    reference_electrode : str, optional
        Reference electrode used.
    ph : float, optional
        Electrolyte pH.
    catalyst_loading_mg_cm2 : float, optional
        Catalyst loading in mg/cm^2.
    ru_ohm : float, optional
        Uncompensated resistance in Ohms.
    encoding : str, default 'latin1'
        File encoding (EC-Lab exports often use latin1 or cp1252).

    Returns
    -------
    StandardizedElectrochemicalData
    """
    p = Path(filepath)
    if not p.is_file():
        raise FileNotFoundError(f"BioLogic .mpt file not found: {p.resolve()}")

    file_hash = compute_file_sha256(p)

    with open(p, "r", encoding=encoding, errors="replace") as f:
        head_sample = [f.readline() for _ in range(120)]

    metadata, n_header = parse_biologic_mpt_header(head_sample)

    try:
        # EC-Lab .mpt columns are tab-separated
        raw_df = pd.read_csv(
            p,
            skiprows=max(0, n_header - 1),
            sep="\t",
            encoding=encoding,
            engine="python",
            on_bad_lines="skip",
        )
    except Exception as e:
        raise BioLogicParseError(f"Failed to parse tabular data from {p.name}: {e}")

    if raw_df.empty:
        raise BioLogicParseError(f"BioLogic file {p.name} contains no valid data rows.")

    raw_df.columns = [str(c).strip() for c in raw_df.columns]

    # Map columns to standardized physical roles
    mapping, scales = identify_columns(raw_df.columns, custom_mapping=custom_mapping)
    std_df = standardize_dataframe(raw_df, mapping, scales)

    # Technique detection from EC-Lab header metadata or columns
    inferred_technique = technique
    if inferred_technique is None:
        tech_str = ""
        for k, v in metadata.items():
            if any(term in k.lower() for term in ("technique", "method", "loop", "experiment")):
                tech_str += " " + str(v).upper()

        if "CV" in tech_str or "CYCLIC" in tech_str:
            inferred_technique = "CV"
        elif "LSV" in tech_str or "LINEAR" in tech_str:
            inferred_technique = "LSV"
        elif "PEIS" in tech_str or "GEIS" in tech_str or "IMPEDANCE" in tech_str:
            inferred_technique = "EIS"
        elif "CA" in tech_str or "CHRONOAMPEROMETRY" in tech_str:
            inferred_technique = "CA"
        elif "CP" in tech_str or "CHRONOPOTENTIOMETRY" in tech_str:
            inferred_technique = "CP"
        else:
            if "frequency_Hz" in std_df.columns:
                inferred_technique = "EIS"
            elif "potential_V" in std_df.columns and "current_A" in std_df.columns:
                pot = std_df["potential_V"].dropna().values
                if len(pot) > 10 and (np.diff(pot) > 0).any() and (np.diff(pot) < 0).any():
                    inferred_technique = "CV"
                else:
                    inferred_technique = "LSV"
            else:
                inferred_technique = "Unknown"

    # Scan rate detection from EC-Lab metadata (e.g. dE/dt = 20.000 mV/s)
    scan_rate: Optional[float] = None
    for k, v in metadata.items():
        if "de/dt" in k.lower() or "scan rate" in k.lower():
            m = re.search(r"([-+]?\d+\.?\d*)", str(v))
            if m:
                val = float(m.group(1))
                if "mv/s" in str(v).lower():
                    scan_rate = val * 1e-3
                else:
                    scan_rate = val
                break

    if scan_rate is None and "time_s" in std_df.columns and "potential_V" in std_df.columns:
        dt = np.diff(std_df["time_s"].values)
        de = np.diff(std_df["potential_V"].values)
        valid = dt > 0
        if np.any(valid):
            rates = np.abs(de[valid] / dt[valid])
            med_rate = float(np.median(rates))
            if 1e-5 < med_rate < 100.0:
                scan_rate = med_rate

    pot_arr = std_df["potential_V"].values if "potential_V" in std_df.columns else None
    cur_arr = std_df["current_A"].values if "current_A" in std_df.columns else None
    time_arr = std_df["time_s"].values if "time_s" in std_df.columns else None
    cycle_arr = std_df["cycle_number"].values if "cycle_number" in std_df.columns else None

    # BioLogic sometimes records cycle number 0 for initial cycle; offset to 1 if starts at 0
    if cycle_arr is not None and len(cycle_arr) > 0 and np.min(cycle_arr) == 0:
        cycle_arr = cycle_arr + 1
        std_df["cycle_number"] = cycle_arr

    # Metadata extraction for physical parameters
    if electrode_area_cm2 is None:
        for k, v in metadata.items():
            if any(term in k.lower() for term in ("area", "surface area", "electrode surface area")):
                m = re.search(r"([-+]?\d+\.?\d*)", str(v))
                if m:
                    electrode_area_cm2 = float(m.group(1))
                    break

    if ph is None:
        for k, v in metadata.items():
            if "ph" in k.lower() or "electrolyte" in k.lower():
                m = re.search(r"ph\s*[:=]?\s*(\d+\.?\d*)", str(v), re.IGNORECASE)
                if m:
                    ph = float(m.group(1))
                    break

    if ru_ohm is None:
        for k, v in metadata.items():
            if any(term in k.lower() for term in ("ru", "uncompensated", "resistance")):
                m = re.search(r"([-+]?\d+\.?\d*)", str(v))
                if m:
                    ru_ohm = float(m.group(1))
                    break


    if catalyst_loading_mg_cm2 is None:
        for k, v in metadata.items():
            if "characteristic mass" in k.lower() or "loading" in k.lower():
                m = re.search(r"([-+]?\d+\.?\d*)", str(v))
                if m:
                    catalyst_loading_mg_cm2 = float(m.group(1))
                    break

    # Guard required channels
    if inferred_technique != "EIS" and pot_arr is None and cur_arr is None:
        raise MissingChannelError(
            f"BioLogic file {p.name} lacks recognizable potential and current columns. "
            f"Found columns: {list(raw_df.columns)}. Try passing custom_mapping."
        )

    reference_electrode, temperature_k = resolve_experimental_conditions(metadata, reference_electrode, temperature_k)

    return StandardizedElectrochemicalData(
        filepath=p.resolve(),
        file_hash=file_hash,
        instrument="biologic",
        instrument_software="EC-Lab",
        technique=inferred_technique,
        potential_V=pot_arr,
        current_A=cur_arr,
        time_s=time_arr,
        cycle_number=cycle_arr,
        segment_id=std_df["segment_id"].values if "segment_id" in std_df.columns else None,
        frequency_Hz=std_df["frequency_Hz"].values if "frequency_Hz" in std_df.columns else None,
        impedance_real_ohm=std_df["impedance_real_ohm"].values if "impedance_real_ohm" in std_df.columns else None,
        impedance_imag_ohm=std_df["impedance_imag_ohm"].values if "impedance_imag_ohm" in std_df.columns else None,
        impedance_mag_ohm=std_df["impedance_mag_ohm"].values if "impedance_mag_ohm" in std_df.columns else None,
        phase_angle_deg=std_df["phase_angle_deg"].values if "phase_angle_deg" in std_df.columns else None,
        scan_rate_V_s=scan_rate,
        potential_reference_raw=reference_electrode,
        temperature_K=temperature_k,
        ph=ph,
        electrode_area_cm2=electrode_area_cm2,
        catalyst_loading_mg_cm2=catalyst_loading_mg_cm2,
        ru_ohm=ru_ohm,
        raw_metadata=metadata,
        dataframe=std_df,
    )


def parse_biologic_mpr(
    filepath: Union[str, Path],
    technique: Optional[str] = None,
    custom_mapping: Optional[Dict[str, str]] = None,
    electrode_area_cm2: Optional[float] = None,
    reference_electrode: Optional[str] = None,
    ph: Optional[float] = None,
    catalyst_loading_mg_cm2: Optional[float] = None,
    ru_ohm: Optional[float] = None,
    temperature_k: Optional[float] = None,
) -> StandardizedElectrochemicalData:
    """
    Parse BioLogic binary .mpr file using galvani if available.

    If galvani is unavailable or fails on an unsupported binary version,
    raises BioLogicBinaryMprError explaining how to export .mpt from EC-Lab.
    """
    p = Path(filepath)
    if not p.is_file():
        raise FileNotFoundError(f"BioLogic .mpr file not found: {p.resolve()}")

    if not GALVANI_AVAILABLE:
        raise BioLogicBinaryMprError(
            f"Cannot parse binary .mpr file '{p.name}': 'galvani' library is not installed. "
            "Please export the data as an ASCII (.mpt) file from EC-Lab (File -> Export to text), "
            "or install galvani."
        )

    file_hash = compute_file_sha256(p)

    try:
        mpr_obj = galvani.BioLogic.MPRfile(str(p))
        df = pd.DataFrame(mpr_obj.data)
    except Exception as e:
        raise BioLogicBinaryMprError(
            f"Failed to parse binary .mpr file '{p.name}' using galvani: {e}. "
            "Specific binary EC-Lab versions or encrypted formats may not be supported. "
            "Please export the data as .mpt ASCII from EC-Lab."
        )

    df.columns = [str(c).strip() for c in df.columns]

    mapping, scales = identify_columns(df.columns, custom_mapping=custom_mapping)
    std_df = standardize_dataframe(df, mapping, scales)

    inferred_technique = technique or "Unknown"
    if "frequency_Hz" in std_df.columns:
        inferred_technique = "EIS"
    elif "potential_V" in std_df.columns and "current_A" in std_df.columns:
        pot = std_df["potential_V"].dropna().values
        if len(pot) > 10 and (np.diff(pot) > 0).any() and (np.diff(pot) < 0).any():
            inferred_technique = "CV"
        else:
            inferred_technique = "LSV"

    metadata = {}
    reference_electrode, temperature_k = resolve_experimental_conditions(metadata, reference_electrode, temperature_k)

    return StandardizedElectrochemicalData(
        filepath=p.resolve(),
        file_hash=file_hash,
        instrument="biologic",
        instrument_software="EC-Lab",
        technique=inferred_technique,
        potential_V=std_df["potential_V"].values if "potential_V" in std_df.columns else None,
        current_A=std_df["current_A"].values if "current_A" in std_df.columns else None,
        time_s=std_df["time_s"].values if "time_s" in std_df.columns else None,
        cycle_number=std_df["cycle_number"].values if "cycle_number" in std_df.columns else None,
        segment_id=std_df["segment_id"].values if "segment_id" in std_df.columns else None,
        frequency_Hz=std_df["frequency_Hz"].values if "frequency_Hz" in std_df.columns else None,
        impedance_real_ohm=std_df["impedance_real_ohm"].values if "impedance_real_ohm" in std_df.columns else None,
        impedance_imag_ohm=std_df["impedance_imag_ohm"].values if "impedance_imag_ohm" in std_df.columns else None,
        impedance_mag_ohm=std_df["impedance_mag_ohm"].values if "impedance_mag_ohm" in std_df.columns else None,
        phase_angle_deg=std_df["phase_angle_deg"].values if "phase_angle_deg" in std_df.columns else None,
        potential_reference_raw=reference_electrode,
        temperature_K=temperature_k,
        ph=ph,
        electrode_area_cm2=electrode_area_cm2,
        catalyst_loading_mg_cm2=catalyst_loading_mg_cm2,
        ru_ohm=ru_ohm,
        raw_metadata=metadata,
        dataframe=std_df,
    )


def load_biologic(
    filepath: Union[str, Path],
    **kwargs: Any,
) -> StandardizedElectrochemicalData:
    """
    Convenience loader auto-dispatching between .mpt (ASCII) and .mpr (binary).
    """
    p = Path(filepath)
    ext = p.suffix.lower()
    if ext == ".mpt":
        return parse_biologic_mpt(p, **kwargs)
    elif ext == ".mpr":
        return parse_biologic_mpr(p, **kwargs)
    else:
        # Default try as text MPT
        return parse_biologic_mpt(p, **kwargs)

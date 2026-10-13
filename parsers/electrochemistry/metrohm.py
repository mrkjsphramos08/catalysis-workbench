"""
Metrohm Autolab NOVA Data Parser.

Parses ASCII/TXT, CSV, and Excel (XLSX) exports produced by Metrohm NOVA software.
Supports customized signal arrangements, multi-segment/multi-cycle datasets,
metadata headers, and converts raw measurements into StandardizedElectrochemicalData.
Rejects proprietary binary .nox projects with clear user instructions.
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
)


class MetrohmParseError(ElectrochemistryParseError):
    """Raised when Metrohm NOVA file parsing encounters errors."""
    pass


class MetrohmBinaryNoxError(MetrohmParseError):
    """Raised when proprietary binary .nox project files are supplied."""
    pass


def _detect_delimiter(sample_lines: Sequence[str]) -> str:
    """Detect delimiter among tab, comma, semicolon, whitespace."""
    counts = {"\t": 0, ",": 0, ";": 0}
    for line in sample_lines:
        line_s = line.strip()
        if not line_s or line_s.startswith("#"):
            continue
        for delim in counts:
            counts[delim] += line_s.count(delim)

    best_delim = max(counts, key=counts.get)
    if counts[best_delim] > 0:
        return best_delim
    return r"\s+"


def _parse_nova_metadata_lines(lines: Sequence[str]) -> Tuple[Dict[str, Any], int]:
    """
    Parse metadata headers from NOVA text export and find start of numerical table.

    Returns
    -------
    Tuple[Dict[str, Any], int]
        - metadata dictionary
        - line index where column header starts
    """
    metadata: Dict[str, Any] = {}
    header_line_idx = 0

    for idx, line in enumerate(lines):
        line_str = line.strip()
        if not line_str:
            continue

        # If line starts with '#', check if it is metadata or commented header
        if line_str.startswith("#"):
            if ":" in line_str:
                parts = line_str.split(":", 1)
                k = parts[0].strip().replace("#", "").strip()
                v = parts[1].strip()
                metadata[k] = v
                continue
            # If line is '# col1  col2 ...', check if it's the header row
            cleaned = line_str.lstrip("#").strip().lower()
            if any(term in cleaned for term in ("potential", "current", "frequency", "time", "we(1)", "z'")):
                header_line_idx = idx
                break
            continue

        # Non-comment line
        lower_line = line_str.lower()
        if any(term in lower_line for term in ("potential", "current", "ewe", "time", "frequency", "we(1)", "z'")):
            header_line_idx = idx
            break

    return metadata, header_line_idx


def parse_metrohm_export(
    filepath: Union[str, Path],
    technique: Optional[str] = None,
    custom_mapping: Optional[Dict[str, str]] = None,
    electrode_area_cm2: Optional[float] = None,
    reference_electrode: str = "Ag/AgCl",
    ph: Optional[float] = None,
    catalyst_loading_mg_cm2: Optional[float] = None,
    ru_ohm: Optional[float] = None,
    delimiter: Optional[str] = None,
    encoding: str = "utf-8",
) -> StandardizedElectrochemicalData:
    """
    Parse Metrohm Autolab NOVA export (ASCII/TXT, CSV, XLSX) into StandardizedElectrochemicalData.

    Parameters
    ----------
    filepath : str or Path
        Path to exported Metrohm file.
    technique : str, optional
        Electrochemical technique ('CV', 'LSV', 'EIS', 'CA', 'CP').
        If None, automatically inferred from headers or column signals.
    custom_mapping : dict, optional
        User override for signal mapping {raw_header: canonical_role}.
    electrode_area_cm2 : float, optional
        Geometric surface area in cm^2.
    reference_electrode : str, default 'Ag/AgCl'
        Reference electrode used during measurement.
    ph : float, optional
        Electrolyte pH.
    catalyst_loading_mg_cm2 : float, optional
        Catalyst loading in mg/cm^2.
    ru_ohm : float, optional
        Uncompensated solution resistance in Ohms.
    delimiter : str, optional
        Explicit delimiter override. If None, auto-detected.
    encoding : str, default 'utf-8'
        File encoding.

    Returns
    -------
    StandardizedElectrochemicalData
    """
    p = Path(filepath)
    if not p.is_file():
        raise FileNotFoundError(f"Metrohm file not found: {p.resolve()}")

    file_hash = compute_file_sha256(p)
    ext = p.suffix.lower()

    # Reject binary .nox files
    if ext == ".nox":
        raise MetrohmBinaryNoxError(
            f"File '{p.name}' is a proprietary Metrohm NOVA binary project (.nox). "
            "Please open the project in NOVA and export the measurement data as ASCII (.txt, .csv) "
            "or Excel (.xlsx) before importing."
        )

    raw_metadata: Dict[str, Any] = {}
    raw_df: pd.DataFrame

    if ext in (".xlsx", ".xls"):
        # Read Excel export
        try:
            raw_df = pd.read_excel(p, sheet_name=0)
        except Exception as e:
            raise MetrohmParseError(f"Failed to read Metrohm Excel file {p.name}: {e}")
    else:
        # Read text/CSV/ASCII export
        try:
            with open(p, "r", encoding=encoding, errors="replace") as f:
                head_lines = [f.readline() for _ in range(50)]
        except Exception as e:
            raise MetrohmParseError(f"Failed to open Metrohm file {p.name}: {e}")

        raw_metadata, header_line_idx = _parse_nova_metadata_lines(head_lines)

        if delimiter is None:
            sep = _detect_delimiter(head_lines[header_line_idx:])
        else:
            sep = delimiter

        try:
            raw_df = pd.read_csv(
                p,
                skiprows=header_line_idx,
                sep=sep,
                encoding=encoding,
                engine="python",
                on_bad_lines="skip",
            )
        except Exception as e:
            raise MetrohmParseError(f"Failed to parse tabular data from {p.name}: {e}")

    if raw_df.empty:
        raise MetrohmParseError(f"Metrohm file {p.name} contains no valid data rows.")

    # Strip column names
    raw_df.columns = [str(c).strip() for c in raw_df.columns]

    # Column identification and normalization
    mapping, scales = identify_columns(raw_df.columns, custom_mapping=custom_mapping)
    std_df = standardize_dataframe(raw_df, mapping, scales)

    # Technique detection
    inferred_technique = technique
    if inferred_technique is None:
        tech_meta = str(raw_metadata.get("Procedure", raw_metadata.get("Method", ""))).upper()
        if not tech_meta:
            meta_str = " ".join([f"{k} {v}" for k, v in raw_metadata.items()]).upper()
            if "CHRONOPOTENTIOMETRY" in meta_str or " CP " in meta_str:
                tech_meta = "CP"
            elif "CHRONOAMPEROMETRY" in meta_str or " CA " in meta_str:
                tech_meta = "CA"
            elif "IMPEDANCE" in meta_str or "FRA" in meta_str or "EIS" in meta_str:
                tech_meta = "EIS"
            elif "CYCLIC" in meta_str or " CV " in meta_str:
                tech_meta = "CV"
            elif "LINEAR" in meta_str or " LSV " in meta_str:
                tech_meta = "LSV"

        if "CYCLIC" in tech_meta or tech_meta == "CV":
            inferred_technique = "CV"
        elif "LINEAR" in tech_meta or tech_meta == "LSV":
            inferred_technique = "LSV"
        elif "IMPEDANCE" in tech_meta or "EIS" in tech_meta or "FRA" in tech_meta:
            inferred_technique = "EIS"
        elif "POTENTIOMETRY" in tech_meta or tech_meta == "CP":
            inferred_technique = "CP"
        elif "AMPEROMETRY" in tech_meta or tech_meta == "CA":
            inferred_technique = "CA"
        else:
            if "frequency_Hz" in std_df.columns:
                inferred_technique = "EIS"
            elif "potential_V" in std_df.columns and "current_A" in std_df.columns:
                pot = std_df["potential_V"].dropna().values
                cur = std_df["current_A"].dropna().values
                if len(cur) > 10 and (np.ptp(cur) < 1e-6 or (np.ptp(cur) / max(1e-9, np.mean(np.abs(cur)))) < 0.001):
                    inferred_technique = "CP"
                elif len(pot) > 10 and (np.ptp(pot) < 0.005 or (np.ptp(pot) / max(1e-9, np.mean(np.abs(pot)))) < 0.001):
                    inferred_technique = "CA"
                elif len(pot) > 10 and (np.diff(pot) > 0).any() and (np.diff(pot) < 0).any():
                    inferred_technique = "CV"
                else:
                    inferred_technique = "LSV"
            else:
                inferred_technique = "Unknown"

    # Metadata extraction for physical parameters
    if electrode_area_cm2 is None:
        for k, v in raw_metadata.items():
            if any(term in k.lower() for term in ("area", "surface area")):
                m = re.search(r"([-+]?\d+\.?\d*)", str(v))
                if m:
                    electrode_area_cm2 = float(m.group(1))
                    break

    if ph is None:
        for k, v in raw_metadata.items():
            if "ph" in k.lower() or "electrolyte" in k.lower():
                m = re.search(r"ph\s*[:=]?\s*(\d+\.?\d*)", str(v), re.IGNORECASE)
                if m:
                    ph = float(m.group(1))
                    break

    if ru_ohm is None:
        for k, v in raw_metadata.items():
            if any(term in k.lower() for term in ("ru", "uncompensated", "r_u", "resistance")):
                m = re.search(r"([-+]?\d+\.?\d*)", str(v))
                if m:
                    ru_ohm = float(m.group(1))
                    break

    if reference_electrode == "Ag/AgCl" or not reference_electrode:
        for k, v in raw_metadata.items():
            if "reference" in k.lower():
                reference_electrode = str(v).strip()
                break

    if catalyst_loading_mg_cm2 is None:
        for k, v in raw_metadata.items():
            if "characteristic mass" in k.lower() or "loading" in k.lower():
                m = re.search(r"([-+]?\d+\.?\d*)", str(v))
                if m:
                    catalyst_loading_mg_cm2 = float(m.group(1))
                    break

    # Scan rate extraction
    scan_rate: Optional[float] = None
    if "Scan rate" in raw_metadata:
        sr_match = re.search(r"([-+]?\d+\.?\d*)", str(raw_metadata["Scan rate"]))
        if sr_match:
            val = float(sr_match.group(1))
            if "mv/s" in str(raw_metadata["Scan rate"]).lower():
                scan_rate = val * 1e-3
            else:
                scan_rate = val

    # If scan rate not in metadata, estimate from time and potential in linear sweep
    if scan_rate is None and "time_s" in std_df.columns and "potential_V" in std_df.columns:
        dt = np.diff(std_df["time_s"].values)
        de = np.diff(std_df["potential_V"].values)
        valid = dt > 0
        if np.any(valid):
            rates = np.abs(de[valid] / dt[valid])
            med_rate = float(np.median(rates))
            if 1e-5 < med_rate < 100.0:
                scan_rate = med_rate

    # Multi-cycle / segment extraction
    cycle_arr = std_df["cycle_number"].values if "cycle_number" in std_df.columns else None
    if cycle_arr is None and inferred_technique == "CV" and "potential_V" in std_df.columns:
        # Heuristic cycle assignment for continuous CV based on extrema
        pot = std_df["potential_V"].values
        diffs = np.diff(pot)
        sign_changes = np.where(np.diff(np.sign(diffs)))[0]
        # Every 2 sign changes is approximately 1 full cycle
        cycles = np.ones(len(pot), dtype=int)
        if len(sign_changes) >= 2:
            cycle_idx = 1
            for sc_i in range(1, len(sign_changes), 2):
                split_pt = sign_changes[sc_i] + 1
                cycles[split_pt:] = cycle_idx + 1
                cycle_idx += 1
        cycle_arr = cycles
        std_df["cycle_number"] = cycle_arr

    pot_arr = std_df["potential_V"].values if "potential_V" in std_df.columns else None
    cur_arr = std_df["current_A"].values if "current_A" in std_df.columns else None
    time_arr = std_df["time_s"].values if "time_s" in std_df.columns else None

    # Guard required channels for non-EIS
    if inferred_technique != "EIS" and pot_arr is None and cur_arr is None:
        raise MissingChannelError(
            f"Metrohm file {p.name} lacks recognizable potential and current columns. "
            f"Found columns: {list(raw_df.columns)}. Try passing custom_mapping."
        )

    return StandardizedElectrochemicalData(
        filepath=p.resolve(),
        file_hash=file_hash,
        instrument="metrohm",
        instrument_software="NOVA",
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
        electrolyte=raw_metadata.get("Electrolyte"),
        ph=ph,
        electrode_area_cm2=electrode_area_cm2,
        catalyst_loading_mg_cm2=catalyst_loading_mg_cm2,
        ru_ohm=ru_ohm,
        raw_metadata=raw_metadata,
        dataframe=std_df,
    )

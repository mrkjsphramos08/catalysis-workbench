"""
Experimental Electrochemistry Parsers for Metrohm Autolab and BioLogic Workstations.
"""

from pathlib import Path
from typing import Any, Dict, Optional, Union

from parsers.electrochemistry.biologic import (
    BioLogicBinaryMprError,
    BioLogicParseError,
    load_biologic,
    parse_biologic_mpr,
    parse_biologic_mps,
    parse_biologic_mpt,
)
from parsers.electrochemistry.column_mapping import (
    CANONICAL_CHANNELS,
    identify_columns,
    standardize_dataframe,
)
from parsers.electrochemistry.metrohm import (
    MetrohmBinaryNoxError,
    MetrohmParseError,
    parse_metrohm_export,
)
from parsers.electrochemistry.schema import (
    ElectrochemistryParseError,
    MissingChannelError,
    StandardizedElectrochemicalData,
    UnsupportedInstrumentError,
    compute_file_sha256,
)


def load_electrochemistry_data(
    filepath: Union[str, Path],
    instrument: Optional[str] = None,
    **kwargs: Any,
) -> StandardizedElectrochemicalData:
    """
    Unified entrypoint for loading electrochemical datasets from either Metrohm or BioLogic.

    Auto-detects instrument if not explicitly supplied:
    - Files with .mpt, .mpr, .mps -> BioLogic EC-Lab
    - Files containing 'NOVA', 'WE(1)', or other Metrohm headers -> Metrohm Autolab
    """
    p = Path(filepath)
    ext = p.suffix.lower()

    inst = instrument.strip().lower() if instrument else None

    if inst == "biologic" or (inst is None and ext in (".mpt", ".mpr")):
        return load_biologic(p, **kwargs)
    elif inst == "metrohm":
        return parse_metrohm_export(p, **kwargs)
    else:
        # Heuristic detection from contents
        try:
            with open(p, "r", encoding="utf-8", errors="replace") as f:
                head = f.read(2048)
            if "ec-lab" in head.lower() or "nb header lines" in head.lower():
                return load_biologic(p, **kwargs)
            else:
                return parse_metrohm_export(p, **kwargs)
        except Exception:
            return parse_metrohm_export(p, **kwargs)


__all__ = [
    "StandardizedElectrochemicalData",
    "ElectrochemistryParseError",
    "UnsupportedInstrumentError",
    "MissingChannelError",
    "MetrohmParseError",
    "MetrohmBinaryNoxError",
    "BioLogicParseError",
    "BioLogicBinaryMprError",
    "parse_metrohm_export",
    "parse_biologic_mpt",
    "parse_biologic_mpr",
    "parse_biologic_mps",
    "load_biologic",
    "load_electrochemistry_data",
    "identify_columns",
    "standardize_dataframe",
    "compute_file_sha256",
    "CANONICAL_CHANNELS",
]

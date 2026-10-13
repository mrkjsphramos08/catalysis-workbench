"""
Unified Characterization Data Schema and Provenance Tracking.

Provides CharacterizationData and base parser utilities across XRD, FTIR, Raman,
Microscopy (SEM/TEM), ICP-OES/MS, and BET/BJH techniques.
Preserves raw experimental measurements, provenance hashes, unit definitions,
and explicit processing histories.
"""

from dataclasses import dataclass, field
import hashlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import numpy as np
import pandas as pd


def compute_file_sha256(filepath: Union[str, Path]) -> str:
    """Calculate SHA-256 cryptographic digest of a data file."""
    p = Path(filepath)
    if not p.exists() or not p.is_file():
        return "nonexistent"
    hasher = hashlib.sha256()
    with open(p, "rb") as f:
        while chunk := f.read(65536):
            hasher.update(chunk)
    return hasher.hexdigest()


class CharacterizationParseError(Exception):
    """Raised when an instrument export fails to parse."""
    pass


class InvalidUnitsError(CharacterizationParseError):
    """Raised when units are ambiguous, nonphysical, or missing."""
    pass


@dataclass
class CharacterizationData:
    """
    Standardized internal representation of materials characterization datasets.

    Stores continuous 1D signals (XRD, FTIR, Raman, Isotherms, ICP calibrations)
    or tabular multi-channel measurements with complete metadata and history.
    """

    technique: str  # 'XRD', 'FTIR', 'Raman', 'SEM', 'TEM', 'HRTEM', 'ICP', 'BET'
    sample_id: str
    source_file: Path
    file_sha256: str
    x_name: str = "x"
    x_unit: str = "a.u."
    y_name: str = "y"
    y_unit: str = "a.u."
    x: Optional[np.ndarray] = None
    y: Optional[np.ndarray] = None
    y_err: Optional[np.ndarray] = None
    dataframe: Optional[pd.DataFrame] = None
    image: Optional[np.ndarray] = None  # 2D/3D array for microscopy
    pixel_size_nm: Optional[float] = None  # Physical calibration for microscopy (nm/pixel)
    measurement_id: Optional[str] = None
    instrument_manufacturer: Optional[str] = None
    instrument_model: Optional[str] = None
    measurement_date: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    processing_history: List[Dict[str, Any]] = field(default_factory=list)
    parser_version: str = "1.0.0"

    def __post_init__(self) -> None:
        """Harmonize arrays and validate shape consistency."""
        if self.x is not None:
            self.x = np.asarray(self.x, dtype=float)
        if self.y is not None:
            self.y = np.asarray(self.y, dtype=float)
        if self.y_err is not None:
            self.y_err = np.asarray(self.y_err, dtype=float)

        if self.x is not None and self.y is not None:
            if len(self.x) != len(self.y):
                raise ValueError(
                    f"Array length mismatch: len(x)={len(self.x)} != len(y)={len(self.y)}"
                )

    @property
    def n_points(self) -> int:
        """Number of discrete coordinate points recorded."""
        if self.x is not None:
            return len(self.x)
        if self.dataframe is not None:
            return len(self.dataframe)
        return 0

    def add_processing_step(self, action: str, parameters: Dict[str, Any], notes: str = "") -> None:
        """Record an explicit data modification step in provenance history."""
        self.processing_history.append({
            "step": len(self.processing_history) + 1,
            "action": action,
            "parameters": parameters,
            "notes": notes,
        })

    def copy(self) -> "CharacterizationData":
        """Produce a deep copy of the dataset."""
        import copy
        return CharacterizationData(
            technique=self.technique,
            sample_id=self.sample_id,
            source_file=Path(self.source_file),
            file_sha256=self.file_sha256,
            x_name=self.x_name,
            x_unit=self.x_unit,
            y_name=self.y_name,
            y_unit=self.y_unit,
            x=self.x.copy() if self.x is not None else None,
            y=self.y.copy() if self.y is not None else None,
            y_err=self.y_err.copy() if self.y_err is not None else None,
            dataframe=self.dataframe.copy() if self.dataframe is not None else None,
            image=self.image.copy() if self.image is not None else None,
            pixel_size_nm=self.pixel_size_nm,
            measurement_id=self.measurement_id,
            instrument_manufacturer=self.instrument_manufacturer,
            instrument_model=self.instrument_model,
            measurement_date=self.measurement_date,
            metadata=copy.deepcopy(self.metadata),
            processing_history=copy.deepcopy(self.processing_history),
            parser_version=self.parser_version,
        )

    def to_provenance_dict(self) -> Dict[str, Any]:
        """Export complete provenance record as a JSON-serializable dictionary."""
        return {
            "technique": self.technique,
            "sample_id": self.sample_id,
            "measurement_id": self.measurement_id,
            "source_file": str(self.source_file),
            "file_sha256": self.file_sha256,
            "instrument": f"{self.instrument_manufacturer or 'Unknown'} {self.instrument_model or ''}".strip(),
            "measurement_date": self.measurement_date,
            "x_axis": {"name": self.x_name, "unit": self.x_unit},
            "y_axis": {"name": self.y_name, "unit": self.y_unit},
            "n_points": self.n_points,
            "pixel_size_nm": self.pixel_size_nm,
            "metadata": self.metadata,
            "processing_history": self.processing_history,
            "parser_version": self.parser_version,
        }


def record_characterization_provenance(
    figure_name: str,
    input_files: Sequence[Union[str, Path]],
    output_dirs: Sequence[Union[str, Path]],
    technique: str,
    metadata: Optional[Dict[str, Any]] = None,
    *,
    overwrite: bool = False,
) -> Path:
    """Record provenance metadata JSON alongside saved characterization figures."""
    from datetime import datetime, timezone
    import json

    file_records = []
    for f in input_files:
        p = Path(f)
        if p.is_file():
            file_records.append({
                "path": str(p.resolve()),
                "filename": p.name,
                "sha256": compute_file_sha256(p),
                "size_bytes": p.stat().st_size,
            })
        else:
            raise FileNotFoundError(f"Provenance input not found: {p}")

    payload: Dict[str, Any] = {
        "figure_name": figure_name,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_files": file_records,
        "technique": technique,
        "parser_version": "1.0.0-char",
        "analysis_metadata": metadata or {},
    }

    from scripts.utils import write_provenance_json
    return write_provenance_json(payload, figure_name, output_dirs, overwrite=overwrite)

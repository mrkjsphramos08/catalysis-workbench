"""
Unified Data Schema and Storage Representation for Experimental Electrochemistry.

Defines StandardizedElectrochemicalData, shared across Metrohm Autolab (NOVA)
and BioLogic (EC-Lab) instruments.
Preserves raw experimental measurements while providing normalized current densities,
reference electrode potential conversions, and iR compensations.
"""

from dataclasses import dataclass, field
import hashlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd


class ElectrochemistryParseError(Exception):
    """Base exception for electrochemical file parsing failures."""
    pass


class UnsupportedInstrumentError(ElectrochemistryParseError):
    """Raised when an unrecognized instrument file is provided."""
    pass


class MissingChannelError(ElectrochemistryParseError):
    """Raised when required physical channels (potential, current, etc.) are absent."""
    pass


@dataclass
class StandardizedElectrochemicalData:
    """
    Standardized internal representation of experimental electrochemical datasets.

    Shared between Metrohm Autolab (NOVA) and BioLogic (EC-Lab) exports.
    Adheres strictly to IUPAC current conventions:
    - Anodic / oxidation currents are positive (+I)
    - Cathodic / reduction currents are negative (-I)
    """

    filepath: Path
    file_hash: str
    instrument: str  # 'metrohm' or 'biologic'
    instrument_software: str  # 'NOVA' or 'EC-Lab'
    technique: str  # 'CV', 'LSV', 'EIS', 'CA', 'CP', 'OCV'
    
    # Primary continuous channels (numpy arrays)
    potential_V: Optional[np.ndarray] = None  # Raw measured working electrode potential (V)
    current_A: Optional[np.ndarray] = None  # Raw measured current in Amperes (A)
    time_s: Optional[np.ndarray] = None  # Elapsed time in seconds (s)

    # Discretization and sequencing
    cycle_number: Optional[np.ndarray] = None  # 1-indexed or 0-indexed cycle indices
    segment_id: Optional[np.ndarray] = None  # Segment/scan/step index

    # EIS specific channels
    frequency_Hz: Optional[np.ndarray] = None  # AC frequency (Hz)
    impedance_real_ohm: Optional[np.ndarray] = None  # Z' (Ohm)
    impedance_imag_ohm: Optional[np.ndarray] = None  # Z'' (Ohm)
    impedance_mag_ohm: Optional[np.ndarray] = None  # |Z| (Ohm)
    phase_angle_deg: Optional[np.ndarray] = None  # Phase angle in degrees

    # Experimental physical parameters
    scan_rate_V_s: Optional[float] = None  # Scan rate in V/s
    temperature_K: Optional[float] = None  # Cell temperature in K
    potential_reference_raw: str = "unknown"  # 'Ag/AgCl', 'SCE', 'Hg/HgO', 'RHE', etc.
    current_sign_convention: str = "IUPAC"  # IUPAC (+ anodic, - cathodic)
    electrolyte: Optional[str] = None  # e.g. '0.1 M KOH', '0.5 M H2SO4'
    ph: Optional[float] = None  # Electrolyte pH
    electrode_area_cm2: Optional[float] = None  # Geometric surface area (cm^2)
    catalyst_loading_mg_cm2: Optional[float] = None  # Mass loading (mg/cm^2)
    ru_ohm: Optional[float] = None  # Uncompensated solution resistance (Ohm)
    ir_compensated: bool = False  # True if hardware/software compensation already applied

    # Identifiers & Raw metadata
    sample_id: Optional[str] = None
    measurement_id: Optional[str] = None
    raw_metadata: Dict[str, Any] = field(default_factory=dict)
    dataframe: Optional[pd.DataFrame] = None

    def __post_init__(self) -> None:
        """Validate and harmonize array lengths."""
        if self.potential_V is not None:
            self.potential_V = np.asarray(self.potential_V, dtype=float)
        if self.current_A is not None:
            self.current_A = np.asarray(self.current_A, dtype=float)
        if self.time_s is not None:
            self.time_s = np.asarray(self.time_s, dtype=float)
        if self.frequency_Hz is not None:
            self.frequency_Hz = np.asarray(self.frequency_Hz, dtype=float)
        if self.impedance_real_ohm is not None:
            self.impedance_real_ohm = np.asarray(self.impedance_real_ohm, dtype=float)
        if self.impedance_imag_ohm is not None:
            self.impedance_imag_ohm = np.asarray(self.impedance_imag_ohm, dtype=float)

    @property
    def n_points(self) -> int:
        """Total number of continuous data points recorded."""
        if self.potential_V is not None:
            return len(self.potential_V)
        if self.frequency_Hz is not None:
            return len(self.frequency_Hz)
        if self.time_s is not None:
            return len(self.time_s)
        return 0

    @property
    def unique_cycles(self) -> List[int]:
        """List of unique cycle numbers present in the dataset."""
        if self.cycle_number is None:
            return [1]
        unique_vals = np.unique(self.cycle_number)
        return sorted([int(c) for c in unique_vals])

    def get_current_density(
        self,
        area_cm2: Optional[float] = None,
        unit: str = "mA/cm2",
    ) -> np.ndarray:
        """
        Calculate normalized current density.

        Parameters
        ----------
        area_cm2 : float, optional
            Electrode geometric area in cm^2. If None, uses self.electrode_area_cm2.
        unit : str, default 'mA/cm2'
            Target current density unit: 'mA/cm2', 'A/cm2', or 'uA/cm2'.

        Returns
        -------
        np.ndarray
            Current density array.
        """
        if self.current_A is None:
            raise MissingChannelError("Dataset contains no measured current channel.")

        area = area_cm2 if area_cm2 is not None else self.electrode_area_cm2
        if area is None or area <= 0.0:
            raise ValueError(
                "Electrode geometric area (cm^2) must be specified and > 0 to normalize current density."
            )

        j_a_cm2 = self.current_A / float(area)
        u = unit.strip().lower()
        if u in ("ma/cm2", "ma_cm2", "ma/cm^2"):
            return j_a_cm2 * 1000.0
        elif u in ("a/cm2", "a_cm2", "a/cm^2"):
            return j_a_cm2
        elif u in ("ua/cm2", "ua_cm2", "ua/cm^2", "µa/cm2"):
            return j_a_cm2 * 1.0e6
        else:
            raise ValueError(f"Unsupported current density unit: '{unit}'. Use 'mA/cm2', 'A/cm2', or 'uA/cm2'.")

    def get_mass_activity(
        self,
        loading_mg_cm2: Optional[float] = None,
        area_cm2: Optional[float] = None,
        unit: str = "A/mg",
    ) -> np.ndarray:
        """
        Calculate catalyst mass-normalized activity.

        Parameters
        ----------
        loading_mg_cm2 : float, optional
            Catalyst mass loading in mg/cm^2.
        area_cm2 : float, optional
            Geometric area in cm^2.
        unit : str, default 'A/mg'
            'A/mg' or 'mA/mg'.

        Returns
        -------
        np.ndarray
            Mass activity array.
        """
        if self.current_A is None:
            raise MissingChannelError("Dataset contains no measured current channel.")

        load = loading_mg_cm2 if loading_mg_cm2 is not None else self.catalyst_loading_mg_cm2
        area = area_cm2 if area_cm2 is not None else self.electrode_area_cm2

        if load is None or load <= 0.0:
            raise ValueError("Catalyst loading (mg/cm^2) must be specified and > 0.")
        if area is None or area <= 0.0:
            raise ValueError("Electrode geometric area (cm^2) must be specified and > 0.")

        total_mass_mg = float(load) * float(area)
        act = self.current_A / total_mass_mg
        if unit.strip().lower() in ("ma/mg", "ma_mg"):
            return act * 1000.0
        return act

    def get_potential_vs_rhe(
        self,
        reference: Optional[str] = None,
        ph: Optional[float] = None,
        temperature_k: Optional[float] = None,
    ) -> np.ndarray:
        """
        Convert measured potential to Reversible Hydrogen Electrode (RHE) scale.

        Parameters
        ----------
        reference : str, optional
            Reference electrode type (e.g. 'Ag/AgCl', 'SCE', 'Hg/HgO').
            If None, uses self.potential_reference_raw.
        ph : float, optional
            Solution pH. If None, uses self.ph.
        temperature_k : float, optional
            Temperature in Kelvin. If None, uses self.temperature_K.

        Returns
        -------
        np.ndarray
            Potential array in Volts vs RHE.
        """
        from modules.electrochemistry.reference_electrodes import convert_potential_to_rhe

        ref = reference if reference is not None else self.potential_reference_raw
        sol_ph = ph if ph is not None else self.ph
        t_k = temperature_k if temperature_k is not None else self.temperature_K

        if ref is None or ref.strip().lower() in ("unknown", "none", ""):
            raise ValueError(
                "Reference electrode type is unknown. Please provide reference='Ag/AgCl' (or similar) and pH."
            )

        if self.potential_V is None:
            raise MissingChannelError("Dataset contains no measured potential channel.")

        return convert_potential_to_rhe(
            potential_measured=self.potential_V,
            reference_electrode=ref,
            ph=sol_ph,
            temperature_k=t_k,
        )

    def get_ir_corrected_potential(
        self,
        ru_ohm: Optional[float] = None,
        compensation_fraction: float = 1.0,
    ) -> np.ndarray:
        """
        Apply software uncompensated solution resistance (iR) correction:
        E_corrected = E_measured - (I * Ru * compensation_fraction).

        Strict safety rule: If self.ir_compensated is True, raises ValueError
        to prevent duplicate double-compensation.

        Parameters
        ----------
        ru_ohm : float, optional
            Uncompensated solution resistance in Ohms. If None, uses self.ru_ohm.
        compensation_fraction : float, default 1.0
            Fraction of Ru to compensate (typically 0.85 to 1.0).

        Returns
        -------
        np.ndarray
            iR-corrected working electrode potential (V).
        """
        if self.potential_V is None:
            raise MissingChannelError("Dataset contains no measured potential channel.")
        if self.current_A is None:
            raise MissingChannelError("Dataset contains no measured current channel for iR correction.")

        if self.ir_compensated:
            raise ValueError(
                "Dataset is already marked as iR-compensated. Duplicate software correction refused."
            )

        r_uncomp = ru_ohm if ru_ohm is not None else self.ru_ohm
        if r_uncomp is None or r_uncomp < 0.0:
            raise ValueError("Uncompensated resistance Ru (Ohm) must be specified and non-negative.")

        # IUPAC convention: anodic current is positive -> potential drops across solution resistance
        # E_actual = E_measured - I * Ru
        from modules.electrochemistry.reference_electrodes import apply_ir_compensation
        return apply_ir_compensation(self.potential_V, self.current_A, r_uncomp, compensation_fraction)

    def get_cycle(self, cycle_idx: int) -> "StandardizedElectrochemicalData":
        """
        Extract data slice corresponding to a single cycle number.

        Parameters
        ----------
        cycle_idx : int
            Cycle number to isolate.

        Returns
        -------
        StandardizedElectrochemicalData
            Filtered copy containing only the selected cycle.
        """
        if self.cycle_number is None:
            return self

        mask = (self.cycle_number == cycle_idx)
        if not np.any(mask):
            raise ValueError(f"Cycle {cycle_idx} not found. Available cycles: {self.unique_cycles}")

        sub_df = self.dataframe.loc[mask].reset_index(drop=True) if self.dataframe is not None else None

        return StandardizedElectrochemicalData(
            filepath=self.filepath,
            file_hash=self.file_hash,
            instrument=self.instrument,
            instrument_software=self.instrument_software,
            technique=self.technique,
            potential_V=self.potential_V[mask] if self.potential_V is not None else None,
            current_A=self.current_A[mask] if self.current_A is not None else None,
            time_s=self.time_s[mask] if self.time_s is not None else None,
            cycle_number=self.cycle_number[mask] if self.cycle_number is not None else None,
            segment_id=self.segment_id[mask] if self.segment_id is not None else None,
            frequency_Hz=self.frequency_Hz[mask] if self.frequency_Hz is not None else None,
            impedance_real_ohm=self.impedance_real_ohm[mask] if self.impedance_real_ohm is not None else None,
            impedance_imag_ohm=self.impedance_imag_ohm[mask] if self.impedance_imag_ohm is not None else None,
            impedance_mag_ohm=self.impedance_mag_ohm[mask] if self.impedance_mag_ohm is not None else None,
            phase_angle_deg=self.phase_angle_deg[mask] if self.phase_angle_deg is not None else None,
            scan_rate_V_s=self.scan_rate_V_s,
            temperature_K=self.temperature_K,
            potential_reference_raw=self.potential_reference_raw,
            current_sign_convention=self.current_sign_convention,
            electrolyte=self.electrolyte,
            ph=self.ph,
            electrode_area_cm2=self.electrode_area_cm2,
            catalyst_loading_mg_cm2=self.catalyst_loading_mg_cm2,
            ru_ohm=self.ru_ohm,
            ir_compensated=self.ir_compensated,
            sample_id=self.sample_id,
            measurement_id=self.measurement_id,
            raw_metadata=dict(self.raw_metadata),
            dataframe=sub_df,
        )


def compute_file_sha256(filepath: Union[str, Path]) -> str:
    """Compute SHA-256 hash for provenance registration."""
    p = Path(filepath)
    hasher = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def resolve_experimental_conditions(metadata, reference_electrode=None, temperature_k=None):
    """Resolve only explicit reference/temperature information, preserving its origin.

    Temperature headers must state K or C; unlabeled numbers remain unknown.
    Caller values override file headers and are recorded as caller supplied.
    """
    import re
    reference = reference_electrode
    ref_source = "caller" if reference is not None else None
    if reference is None:
        for key, value in metadata.items():
            if "reference" in key.lower() and "electrode" in key.lower():
                reference, ref_source = str(value).strip(), str(key)
                break
    temp_source = "caller" if temperature_k is not None else None
    if temperature_k is None:
        for key, value in metadata.items():
            if "temperature" not in key.lower():
                continue
            match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)", str(value))
            if match is None:
                continue
            unit_text = (str(key) + " " + str(value)).lower()
            if re.search(r"\b(k|kelvin)\b", unit_text):
                temperature_k = float(match.group())
            elif re.search(r"\b(c|celsius)\b", unit_text):
                temperature_k = float(match.group()) + 273.15
            else:
                continue
            temp_source = str(key)
            break
    if temperature_k is not None and (not np.isfinite(temperature_k) or temperature_k <= 0):
        raise ValueError("Temperature must be finite and positive in kelvin.")
    reference = reference or "unknown"
    metadata["resolved_conditions"] = {
        "reference_electrode": reference, "reference_source": ref_source,
        "temperature_k": temperature_k, "temperature_source": temp_source,
    }
    return reference, temperature_k

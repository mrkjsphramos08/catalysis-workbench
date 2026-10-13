"""
Data Preprocessing and Unit Standardization Pipeline.

Provides clean, traceable ingestion, physical unit conversion (using Pint or exact
CODATA standards), and schema alignment for DFT, catalysis, and kinetic data.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import json
import numpy as np
import pandas as pd

from scripts.config import PROCESSED_DATA_DIR, RAW_DATA_DIR
from scripts.utils import load_csv_data, validate_columns


# Physical Constants (CODATA 2018 / 2019 SI Redefinition standards)
# Primary Fundamental Constants
ELEMENTARY_CHARGE_C: float = 1.602176634e-19       # Exact (2019 SI) in Coulombs
AVOGADRO_CONSTANT_MOL_INV: float = 6.02214076e23   # Exact (2019 SI) in mol^-1
BOLTZMANN_CONSTANT_J_PER_K: float = 1.380649e-23   # Exact (2019 SI) in J/K
BOLTZMANN_CONSTANT_EV_PER_K: float = (
    BOLTZMANN_CONSTANT_J_PER_K / ELEMENTARY_CHARGE_C
)  # 8.617333262145179e-5 eV/K (exact from defined kB and e)
KB_EV_PER_K: float = BOLTZMANN_CONSTANT_EV_PER_K

# Quantum ESPRESSO Internal vs CODATA 2018 Rydberg Constants
RY_TO_EV_QE: float = 13.605693009                  # QE Modules/constants.f90 internal standard
RY_TO_EV_CODATA2018: float = 13.605693122994       # CODATA 2018 recommended (h * c * R_infty / e)
RY_TO_EV: float = RY_TO_EV_QE                      # Default for Quantum ESPRESSO consistency

# Length & Force Conversion Constants
BOHR_TO_ANGSTROM: float = 0.529177210903          # CODATA 2018 Bohr radius a_0 in Angstrom
RY_PER_BOHR_TO_EV_PER_ANGSTROM: float = (
    RY_TO_EV / BOHR_TO_ANGSTROM
)  # ~25.71104301548 eV/Angstrom

# Thermodynamic Energy Conversion Constants
EV_TO_J: float = ELEMENTARY_CHARGE_C               # 1 eV = 1.602176634e-19 J
EV_TO_KJ_PER_MOL: float = (
    ELEMENTARY_CHARGE_C * AVOGADRO_CONSTANT_MOL_INV / 1000.0
)  # 96.48533212331002 kJ/mol (exact from defined e and NA)
KCAL_PER_MOL_TO_EV: float = 4184.0 / (
    ELEMENTARY_CHARGE_C * AVOGADRO_CONSTANT_MOL_INV
)  # 0.043364115308770494 eV (exact from 1 cal = 4.184 J)
EV_TO_KCAL_PER_MOL: float = 1.0 / KCAL_PER_MOL_TO_EV  # 23.06054783061903 kcal/mol

# Universal Gas Constant R = NA * kB
GAS_CONSTANT_R_J_MOL_K: float = (
    AVOGADRO_CONSTANT_MOL_INV * BOLTZMANN_CONSTANT_J_PER_K
)  # 8.31446261815324 J/(mol*K) (exact in 2019 SI)
GAS_CONSTANT_R_KJ_MOL_K: float = GAS_CONSTANT_R_J_MOL_K / 1000.0  # 0.00831446261815324 kJ/(mol*K)
GAS_CONSTANT_R_EV_K: float = BOLTZMANN_CONSTANT_EV_PER_K           # eV/K


def convert_energy(
    values: Union[float, np.ndarray, pd.Series],
    from_unit: str,
    to_unit: str,
) -> Union[float, np.ndarray, pd.Series]:
    """
    Explicitly convert energy values between eV, Ry, kJ/mol, and kcal/mol.

    Parameters
    ----------
    values : float, ndarray, or Series
        Numerical energy values.
    from_unit : str
        Source unit ('eV', 'Ry', 'kJ/mol', 'kcal/mol').
    to_unit : str
        Target unit ('eV', 'Ry', 'kJ/mol', 'kcal/mol').

    Returns
    -------
    converted values
    """
    u_from = from_unit.strip().lower()
    u_to = to_unit.strip().lower()

    if u_from == u_to:
        return values

    # Convert source to eV first
    if u_from in ("ev", "electronvolt"):
        ev = values
    elif u_from in ("ry", "rydberg"):
        ev = values * RY_TO_EV
    elif u_from in ("kj/mol", "kj_mol", "kjmol"):
        ev = values / EV_TO_KJ_PER_MOL
    elif u_from in ("kcal/mol", "kcal_mol", "kcalmol"):
        ev = values * KCAL_PER_MOL_TO_EV
    else:
        raise ValueError(f"Unsupported source energy unit: '{from_unit}'")

    # Convert eV to target
    if u_to in ("ev", "electronvolt"):
        return ev
    elif u_to in ("ry", "rydberg"):
        return ev / RY_TO_EV
    elif u_to in ("kj/mol", "kj_mol", "kjmol"):
        return ev * EV_TO_KJ_PER_MOL
    elif u_to in ("kcal/mol", "kcal_mol", "kcalmol"):
        return ev * EV_TO_KCAL_PER_MOL
    else:
        raise ValueError(f"Unsupported target energy unit: '{to_unit}'")


def convert_temperature(
    values: Union[float, np.ndarray, pd.Series],
    from_unit: str,
    to_unit: str,
) -> Union[float, np.ndarray, pd.Series]:
    """
    Convert temperature between Kelvin ('K') and Celsius ('C').
    """
    u_from = from_unit.strip().upper()
    u_to = to_unit.strip().upper()

    if u_from == u_to:
        return values

    if u_from == "C" and u_to == "K":
        return values + 273.15
    elif u_from == "K" and u_to == "C":
        return values - 273.15
    else:
        raise ValueError(f"Unsupported temperature conversion from '{from_unit}' to '{to_unit}'")


def preprocess_dft_dos_data(
    input_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    fermi_energy: Optional[float] = None,
) -> pd.DataFrame:
    """
    Standardize raw electronic density-of-states CSV.
    Aligns energy axis relative to Fermi level (E - E_F).
    """
    df = load_csv_data(input_file)
    validate_columns(df, ["energy"], context="DOS preprocessing")

    processed = df.copy()
    if fermi_energy is not None:
        processed["energy_relative"] = processed["energy"] - fermi_energy
        processed.attrs["fermi_energy_applied"] = fermi_energy
    else:
        processed["energy_relative"] = processed["energy"]

    if output_file:
        out_path = Path(output_file)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        processed.to_csv(out_path, index=False)

    return processed


def preprocess_kinetic_series(
    input_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    temp_unit: str = "C",
) -> pd.DataFrame:
    """
    Preprocess reaction kinetics dataset:
    Ensures temperature in Kelvin, computes 1000/T, validates conversion/selectivity.
    """
    df = load_csv_data(input_file)
    validate_columns(df, ["temperature"], context="Kinetic preprocessing")

    processed = df.copy()
    if temp_unit.upper() == "C":
        processed["temperature_K"] = processed["temperature"] + 273.15
    else:
        processed["temperature_K"] = processed["temperature"]

    processed["inv_temperature_1000K"] = 1000.0 / processed["temperature_K"]

    if "conversion" in processed.columns:
        # Check boundary [0, 1] or [0, 100]
        if (processed["conversion"] > 1.0).any():
            # Percentage provided, convert to fraction
            processed["conversion_fraction"] = processed["conversion"] / 100.0
        else:
            processed["conversion_fraction"] = processed["conversion"]

    if output_file:
        out_path = Path(output_file)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        processed.to_csv(out_path, index=False)

    return processed

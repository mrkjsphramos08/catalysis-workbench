"""
Electrochemical Reference Electrodes, Thermodynamic Potentials, and iR Compensation.

Provides certified reference potentials vs SHE/RHE for standard laboratory electrodes:
- Ag/AgCl (Sat. KCl, 3.0 M, 3.5 M, 1.0 M)
- SCE (Saturated Calomel Electrode)
- Hg/HgO (1.0 M KOH, 0.1 M KOH)
- Hg/Hg2SO4 (Sat. K2SO4)
- RHE (Reversible Hydrogen Electrode)
Includes temperature corrections and uncompensated resistance (iR) corrections.
"""

from typing import Dict, Optional, Sequence, Tuple, Union
import numpy as np

from scripts.preprocess import KB_EV_PER_K


# Reference Electrode potentials vs SHE at T = 298.15 K (25 °C)
# Values from Bard & Faulkner, "Electrochemical Methods", and IUPAC standard conventions
REFERENCE_ELECTRODES_298K: Dict[str, Dict[str, float]] = {
    "she": {
        "e0_vs_she": 0.000,
        "dE_dT": 0.000,  # mV/K
        "description": "Standard Hydrogen Electrode (SHE)",
    },
    "nhe": {
        "e0_vs_she": 0.000,
        "dE_dT": 0.000,
        "description": "Normal Hydrogen Electrode (NHE)",
    },
    "sce": {
        "e0_vs_she": 0.2412,
        "dE_dT": -0.65e-3,  # -0.65 mV/K
        "description": "Saturated Calomel Electrode (SCE, Sat. KCl)",
    },
    "ag/agcl_sat": {
        "e0_vs_she": 0.1970,
        "dE_dT": -1.01e-3,  # -1.01 mV/K
        "description": "Silver/Silver Chloride (Ag/AgCl, Sat. KCl)",
    },
    "ag/agcl": {  # Default alias maps to saturated KCl
        "e0_vs_she": 0.1970,
        "dE_dT": -1.01e-3,
        "description": "Silver/Silver Chloride (Ag/AgCl, Sat. KCl)",
    },
    "ag/agcl_3m": {
        "e0_vs_she": 0.2100,
        "dE_dT": -0.70e-3,
        "description": "Silver/Silver Chloride (Ag/AgCl, 3.0 M KCl)",
    },
    "ag/agcl_3.5m": {
        "e0_vs_she": 0.2050,
        "dE_dT": -0.73e-3,
        "description": "Silver/Silver Chloride (Ag/AgCl, 3.5 M KCl)",
    },
    "ag/agcl_1m": {
        "e0_vs_she": 0.2350,
        "dE_dT": -0.62e-3,
        "description": "Silver/Silver Chloride (Ag/AgCl, 1.0 M KCl)",
    },
    "hg/hgo_1m": {
        "e0_vs_she": 0.0980,  # at 1.0 M KOH (~pH 14)
        "dE_dT": -0.28e-3,
        "description": "Mercury/Mercuric Oxide (Hg/HgO, 1.0 M KOH)",
    },
    "hg/hgo": {  # Default alias maps to 1.0 M KOH
        "e0_vs_she": 0.0980,
        "dE_dT": -0.28e-3,
        "description": "Mercury/Mercuric Oxide (Hg/HgO, 1.0 M KOH)",
    },
    "hg/hgo_0.1m": {
        "e0_vs_she": 0.1650,  # at 0.1 M KOH (~pH 13)
        "dE_dT": -0.28e-3,
        "description": "Mercury/Mercuric Oxide (Hg/HgO, 0.1 M KOH)",
    },
    "hg/hg2so4_sat": {
        "e0_vs_she": 0.6400,
        "dE_dT": -0.80e-3,
        "description": "Mercury/Mercurous Sulfate (Hg/Hg2SO4, Sat. K2SO4)",
    },
    "hg/hg2so4": {
        "e0_vs_she": 0.6400,
        "dE_dT": -0.80e-3,
        "description": "Mercury/Mercurous Sulfate (Hg/Hg2SO4, Sat. K2SO4)",
    },
}


def _normalize_ref_key(ref_name: str) -> str:
    """Normalize reference electrode string for case-insensitive lookup."""
    clean = ref_name.strip().lower()
    clean = clean.replace(" ", "")
    clean = clean.replace("(3mkcl)", "_3m").replace("(3.0mkcl)", "_3m")
    clean = clean.replace("(3.5mkcl)", "_3.5m")
    clean = clean.replace("(1mkcl)", "_1m").replace("(1.0mkcl)", "_1m")
    clean = clean.replace("(sat.kcl)", "_sat").replace("(satkcl)", "_sat").replace("(saturated)", "_sat")
    clean = clean.replace("_sat.kcl", "_sat").replace("_saturated", "_sat")
    clean = clean.replace("(1mkoh)", "_1m").replace("(1.0mkoh)", "_1m")
    clean = clean.replace("(0.1mkoh)", "_0.1m")
    clean = clean.replace("(sat.k2so4)", "_sat").replace("(satk2so4)", "_sat")
    return clean


def get_reference_potential_vs_she(
    reference_name: str,
    temperature_k: float = 298.15,
) -> float:
    """
    Get reference electrode potential vs SHE at specified temperature.

    Parameters
    ----------
    reference_name : str
        Name of reference electrode ('Ag/AgCl', 'SCE', 'Hg/HgO', etc.).
    temperature_k : float, default 298.15
        Cell temperature in Kelvin.

    Returns
    -------
    float
        Potential vs SHE in Volts.
    """
    key = _normalize_ref_key(reference_name)
    if key in ("rhe", "reversiblehydrogenelectrode"):
        raise ValueError("RHE is pH-dependent; use convert_potential_to_rhe instead.")

    if key not in REFERENCE_ELECTRODES_298K:
        supported = list(REFERENCE_ELECTRODES_298K.keys())
        raise ValueError(
            f"Unsupported reference electrode: '{reference_name}'. Supported electrodes: {supported}"
        )

    ref_info = REFERENCE_ELECTRODES_298K[key]
    e0_298 = ref_info["e0_vs_she"]
    dE_dT = ref_info.get("dE_dT", 0.0)

    delta_t = temperature_k - 298.15
    return float(e0_298 + dE_dT * delta_t)


# Alias for API consistency
get_standard_electrode_potential = get_reference_potential_vs_she


def convert_potential_to_rhe(
    potential_measured: Union[float, np.ndarray, Sequence[float]],
    reference_electrode: str,
    ph: float,
    temperature_k: float = 298.15,
) -> Union[float, np.ndarray]:
    """
    Convert measured potential to Reversible Hydrogen Electrode (RHE) scale.

    Thermodynamic equation:
    E_RHE = E_measured + E_ref^0(vs SHE) + (k_B * T * ln(10) / e) * pH
    At 298.15 K, the Nernstian slope is ~0.05916 V / pH.

    Parameters
    ----------
    potential_measured : float or array-like
        Measured potential vs reference electrode (V).
    reference_electrode : str
        Reference electrode type ('Ag/AgCl', 'SCE', 'Hg/HgO', 'RHE', etc.).
    ph : float
        Solution electrolyte pH.
    temperature_k : float, default 298.15
        Cell temperature in Kelvin.

    Returns
    -------
    float or np.ndarray
        Potential on RHE scale (V vs RHE).
    """
    key = _normalize_ref_key(reference_electrode)
    arr = np.asarray(potential_measured, dtype=float)

    if key in ("rhe", "reversiblehydrogenelectrode"):
        # Already on RHE scale
        return arr if arr.ndim > 0 else float(arr)

    e_ref_she = get_reference_potential_vs_she(reference_electrode, temperature_k)
    nernst_slope = KB_EV_PER_K * float(temperature_k) * np.log(10.0)
    e_rhe = arr + e_ref_she + (nernst_slope * float(ph))

    return e_rhe if arr.ndim > 0 else float(e_rhe)


def convert_potential_to_she(
    potential_measured: Union[float, np.ndarray, Sequence[float]],
    reference_electrode: str,
    temperature_k: float = 298.15,
) -> Union[float, np.ndarray]:
    """
    Convert measured potential to Standard Hydrogen Electrode (SHE) scale:
    E_SHE = E_measured + E_ref^0(vs SHE).
    """
    arr = np.asarray(potential_measured, dtype=float)
    e_ref_she = get_reference_potential_vs_she(reference_electrode, temperature_k)
    e_she = arr + e_ref_she
    return e_she if arr.ndim > 0 else float(e_she)


def calculate_overpotential(
    potential_vs_rhe: Union[float, np.ndarray, Sequence[float]],
    reaction: str = "OER",
    temperature_k: float = 298.15,
) -> Union[float, np.ndarray]:
    """
    Calculate theoretical overpotential (eta) for HER or OER.

    Thermodynamic standard potentials vs RHE:
    - HER: E_eq = 0.000 V vs RHE
      eta_HER = E_RHE - 0.000 V (cathodic, typically negative; overpotential magnitude is |E_RHE|)
    - OER: E_eq = 1.229 V vs RHE (at 298.15 K)
      eta_OER = E_RHE - 1.229 V (anodic, positive above equilibrium)

    Parameters
    ----------
    potential_vs_rhe : float or array-like
        Working electrode potential vs RHE (V).
    reaction : str, default 'OER'
        Target catalytic reaction: 'HER' or 'OER'.
    temperature_k : float, default 298.15
        Cell temperature in Kelvin.

    Returns
    -------
    float or np.ndarray
        Overpotential in Volts.
    """
    rxn = reaction.strip().upper()
    arr = np.asarray(potential_vs_rhe, dtype=float)

    if rxn == "HER":
        # Equilibrium potential is 0.000 V vs RHE
        eta = arr - 0.000
    elif rxn == "OER":
        # Equilibrium potential is 1.229 V vs RHE at 298.15 K
        # Temperature dependence of water oxidation E0(T) ~ 1.229 - 0.000846 * (T - 298.15)
        e0_oer = 1.229 - 8.46e-4 * (float(temperature_k) - 298.15)
        eta = arr - e0_oer
    else:
        raise ValueError(f"Unsupported reaction: '{reaction}'. Choose 'HER' or 'OER'.")

    return eta if arr.ndim > 0 else float(eta)


def apply_ir_compensation(
    potential_measured: Union[float, np.ndarray, Sequence[float]],
    current_A: Union[float, np.ndarray, Sequence[float]],
    ru_ohm: float,
    compensation_fraction: float = 1.0,
) -> Union[float, np.ndarray]:
    """
    Apply ohmic drop (iR) correction for uncompensated solution resistance Ru:
    E_corrected = E_measured - (I * Ru * compensation_fraction).

    Parameters
    ----------
    potential_measured : float or array-like
        Measured potential in Volts.
    current_A : float or array-like
        Measured cell current in Amperes (A) (IUPAC: positive anodic, negative cathodic).
    ru_ohm : float
        Uncompensated solution resistance in Ohms.
    compensation_fraction : float, default 1.0
        Fraction of Ru to compensate (typically 0.85 to 1.0).

    Returns
    -------
    float or np.ndarray
        iR-compensated potential in Volts.
    """
    if ru_ohm < 0:
        raise ValueError(f"Uncompensated resistance Ru must be non-negative, got {ru_ohm} Ohm.")
    if not (0.0 <= compensation_fraction <= 1.0):
        raise ValueError(
            f"Compensation fraction must be between 0.0 and 1.0, got {compensation_fraction}."
        )

    e_arr = np.asarray(potential_measured, dtype=float)
    i_arr = np.asarray(current_A, dtype=float)

    e_corr = e_arr - (i_arr * float(ru_ohm) * float(compensation_fraction))
    return e_corr if e_arr.ndim > 0 else float(e_corr)

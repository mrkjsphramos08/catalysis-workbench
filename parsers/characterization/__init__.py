"""
Materials Characterization Parsers Package.

Provides unified data ingestion for:
- X-Ray Diffraction (XRD)
- Fourier Transform Infrared Spectroscopy (FTIR)
- Raman Spectroscopy
- Scanning Electron Microscopy (SEM / FESEM)
- Transmission Electron Microscopy (TEM / HRTEM)
- Inductively Coupled Plasma (ICP-OES / ICP-MS)
- Gas Adsorption Physisorption (BET / BJH).
"""

from parsers.characterization.common import (
    CharacterizationData,
    CharacterizationParseError,
    InvalidUnitsError,
    compute_file_sha256,
    record_characterization_provenance,
)
from parsers.characterization.xrd import parse_xrd_file
from parsers.characterization.ftir import parse_ftir_file, convert_transmittance_to_absorbance
from parsers.characterization.raman import parse_raman_file
from parsers.characterization.microscopy import parse_microscopy_image
from parsers.characterization.icp import parse_icp_file
from parsers.characterization.aas import parse_aas_file
from parsers.characterization.adsorption_isotherm import parse_isotherm_file

__all__ = [
    "CharacterizationData",
    "CharacterizationParseError",
    "InvalidUnitsError",
    "compute_file_sha256",
    "record_characterization_provenance",
    "parse_xrd_file",
    "parse_ftir_file",
    "convert_transmittance_to_absorbance",
    "parse_raman_file",
    "parse_microscopy_image",
    "parse_icp_file",
    "parse_aas_file",
    "parse_isotherm_file",
]

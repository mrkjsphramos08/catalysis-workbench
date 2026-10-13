"""
Quantum ESPRESSO bands.x Output Parser.

Parses Gnuplot-formatted band structure files (*.bands.dat.gnu or *.gnu)
produced by Quantum ESPRESSO's bands.x post-processing utility.
Extracts multiple band paths, verifies k-point coordinate integrity,
aligns energies to Fermi level, and integrates user-provided high-symmetry points.
"""

from dataclasses import dataclass, field
import hashlib
from pathlib import Path
import re
from typing import List, Optional, Tuple, Union

import numpy as np

from parsers.quantum_espresso.pw import QEParseError


@dataclass
class QEBandsOutput:
    """Structured container for parsed band structure dispersion curves."""
    filepath: Path
    file_hash: str
    k_distances: np.ndarray
    bands: np.ndarray  # shape: (n_bands, n_kpoints) in eV
    bands_down: Optional[np.ndarray] = None  # shape: (n_bands, n_kpoints)
    is_spin_polarized: bool = False
    e_fermi: Optional[float] = None
    high_symmetry_points: List[Tuple[float, str]] = field(default_factory=list)

    @property
    def n_bands(self) -> int:
        return self.bands.shape[0]

    @property
    def n_kpoints(self) -> int:
        return len(self.k_distances)

    def get_aligned_bands(self, fermi_energy: Optional[float] = None) -> np.ndarray:
        """Return band energies relative to Fermi level (E - E_F)."""
        ef = fermi_energy if fermi_energy is not None else self.e_fermi
        if ef is not None:
            return self.bands - float(ef)
        return self.bands.copy()

    def get_aligned_bands_down(self, fermi_energy: Optional[float] = None) -> Optional[np.ndarray]:
        """Return spin-down band energies relative to Fermi level."""
        if self.bands_down is None:
            return None
        ef = fermi_energy if fermi_energy is not None else self.e_fermi
        if ef is not None:
            return self.bands_down - float(ef)
        return self.bands_down.copy()


def parse_bands_gnu(
    filepath: Union[str, Path],
    spin_down_filepath: Optional[Union[str, Path]] = None,
    e_fermi: Optional[float] = None,
    high_symmetry_points: Optional[List[Tuple[float, str]]] = None,
) -> QEBandsOutput:
    """
    Parse Quantum ESPRESSO bands.x .gnu output file.

    Parameters
    ----------
    filepath : str or Path
        Path to the primary (or spin-up) .gnu file.
    spin_down_filepath : str or Path, optional
        Path to the companion spin-down .gnu file if spin-polarized calculation.
    e_fermi : float, optional
        Fermi energy in eV.
    high_symmetry_points : list of (float, str), optional
        Explicit user-provided high-symmetry points: [(0.0, r'$\\Gamma$'), (0.866, 'K'), ...].
        Do not fabricate labels if unknown.

    Returns
    -------
    QEBandsOutput
        Parsed band structure object.
    """
    path = Path(filepath)
    if not path.is_file():
        raise FileNotFoundError(f"Band structure file not found: {path.resolve()}")

    hasher = hashlib.sha256()
    with open(path, "rb") as f:
        hasher.update(f.read())
    file_hash = hasher.hexdigest()

    k_coords, bands_matrix = _read_gnu_blocks(path)

    bands_down_matrix = None
    is_spin = False
    if spin_down_filepath is not None:
        p_down = Path(spin_down_filepath)
        if not p_down.is_file():
            raise FileNotFoundError(f"Spin-down band file not found: {p_down.resolve()}")
        k_down, bands_down_matrix = _read_gnu_blocks(p_down)
        if len(k_down) != len(k_coords) or not np.allclose(k_down, k_coords, atol=1e-5):
            raise QEParseError("K-point grids do not match between spin-up and spin-down band files.")
        is_spin = True

    return QEBandsOutput(
        filepath=path.resolve(),
        file_hash=file_hash,
        k_distances=k_coords,
        bands=bands_matrix,
        bands_down=bands_down_matrix,
        is_spin_polarized=is_spin,
        e_fermi=e_fermi,
        high_symmetry_points=high_symmetry_points or [],
    )


def _read_gnu_blocks(filepath: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Read blank-line separated blocks of k-points and band energies from a .gnu file."""
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        content = f.read().strip()

    if not content:
        raise QEParseError(f"Band file {filepath.name} is empty.")

    # In Gnuplot data files, bands are separated by double newlines or blank lines
    blocks = [b.strip() for b in re.split(r"\n\s*\n+", content) if b.strip()]

    if not blocks:
        raise QEParseError(f"Could not parse band blocks from {filepath.name}.")

    bands_list: List[List[float]] = []
    ref_k: Optional[np.ndarray] = None

    for b_idx, block in enumerate(blocks):
        lines = block.split("\n")
        k_vals: List[float] = []
        e_vals: List[float] = []

        for line in lines:
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue
            tokens = line_str.split()
            if len(tokens) >= 2:
                try:
                    k_vals.append(float(tokens[0]))
                    e_vals.append(float(tokens[1]))
                except ValueError:
                    continue

        if not k_vals:
            continue

        k_arr = np.array(k_vals, dtype=float)
        if ref_k is None:
            ref_k = k_arr
        elif len(k_arr) != len(ref_k):
            raise QEParseError(
                f"Inconsistent k-point count in band #{b_idx + 1} ({len(k_arr)}) vs first band ({len(ref_k)})."
            )

        bands_list.append(e_vals)

    if ref_k is None or not bands_list:
        raise QEParseError(f"No valid numerical band data could be extracted from {filepath.name}.")

    bands_matrix = np.array(bands_list, dtype=float)  # shape: (n_bands, n_kpoints)
    return ref_k, bands_matrix

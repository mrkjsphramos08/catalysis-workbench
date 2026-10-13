"""
Quantum ESPRESSO neb.x Nudged Elastic Band Parser.

Parses NEB reaction profile files: prefix.dat (discrete images),
prefix.int (interpolated splined pathway), and prefix.path / output logs
(forward and reverse activation barriers, convergence status).
Does not accept input configuration files as output profiles.
"""

from dataclasses import dataclass, field
import hashlib
from pathlib import Path
import re
from typing import List, Optional, Tuple, Union

import numpy as np

from parsers.quantum_espresso.pw import QEParseError


@dataclass
class QENebOutput:
    """Structured container for Nudged Elastic Band reaction pathway results."""
    filepath: Path
    file_hash: str
    image_indices: np.ndarray
    image_coordinates: np.ndarray
    image_energies_ev: np.ndarray
    raw_energies_ev: np.ndarray
    interpolated_coordinates: Optional[np.ndarray] = None
    interpolated_energies_ev: Optional[np.ndarray] = None
    forward_barrier_ev: Optional[float] = None
    reverse_barrier_ev: Optional[float] = None
    reaction_energy_ev: Optional[float] = None
    is_converged: bool = False
    climbing_image_index: Optional[int] = None
    raw_metadata: dict = field(default_factory=dict)


def parse_neb_output(
    dat_filepath: Union[str, Path],
    int_filepath: Optional[Union[str, Path]] = None,
    path_filepath: Optional[Union[str, Path]] = None,
) -> QENebOutput:
    """
    Parse Quantum ESPRESSO neb.x reaction profile output.

    Parameters
    ----------
    dat_filepath : str or Path
        Path to prefix.dat file containing discrete image energies.
    int_filepath : str or Path, optional
        Path to prefix.int containing splined interpolated profile.
    path_filepath : str or Path, optional
        Path to prefix.path or neb log file with barriers and convergence.

    Returns
    -------
    QENebOutput
        Parsed NEB pathway and barrier statistics.
    """
    p_dat = Path(dat_filepath)
    if not p_dat.is_file():
        raise FileNotFoundError(f"NEB dat file not found: {p_dat.resolve()}")

    with open(p_dat, "r", encoding="utf-8", errors="replace") as f:
        dat_content = f.read()

    # Reject input namelist files
    if "&path" in dat_content.lower() or "&input" in dat_content.lower():
        raise QEParseError(
            f"File {p_dat.name} appears to be a NEB input namelist, not a calculated energy profile."
        )

    hasher = hashlib.sha256()
    hasher.update(dat_content.encode("utf-8"))
    file_hash = hasher.hexdigest()

    # Parse prefix.dat
    # Lines typically: image_idx  reaction_coordinate  energy_eV  [error/force]
    # or: coord  energy_eV
    dat_lines = [l.strip() for l in dat_content.split("\n") if l.strip() and not l.strip().startswith("#")]
    if not dat_lines:
        raise QEParseError(f"No numerical lines found in NEB data file {p_dat.name}.")

    coords: List[float] = []
    energies: List[float] = []
    indices: List[int] = []

    for idx, line in enumerate(dat_lines):
        tokens = [float(t) for t in line.split()]
        if len(tokens) >= 3:
            # Format: Image, Coordinate, Energy
            indices.append(int(tokens[0]))
            coords.append(tokens[1])
            energies.append(tokens[2])
        elif len(tokens) == 2:
            indices.append(idx)
            coords.append(tokens[0])
            energies.append(tokens[1])
        else:
            raise QEParseError(f"Unexpected token count in NEB line: '{line}'")

    img_coords = np.array(coords, dtype=float)
    raw_energies = np.array(energies, dtype=float)
    img_indices = np.array(indices, dtype=int)

    # Relative energy with respect to initial state (Image 0 = 0.0 eV)
    rel_energies = raw_energies - raw_energies[0]
    reaction_energy = float(rel_energies[-1])

    # Parse interpolated file if present or auto-discover
    interpolated_coords = None
    interpolated_energies = None

    if int_filepath is None:
        auto_int = p_dat.with_suffix(".int")
        if auto_int.is_file():
            int_filepath = auto_int

    if int_filepath is not None:
        p_int = Path(int_filepath)
        if p_int.is_file():
            try:
                int_data = np.loadtxt(p_int, comments="#")
                if int_data.ndim == 2 and int_data.shape[1] >= 2:
                    interpolated_coords = int_data[:, 0]
                    # Convert to relative energy if needed
                    int_raw_e = int_data[:, 1]
                    interpolated_energies = int_raw_e - raw_energies[0]
            except Exception:
                pass

    # Forward and Reverse barriers
    fwd_barrier = None
    rev_barrier = None
    is_converged = False
    climbing_idx = None

    if path_filepath is None:
        auto_path = p_dat.with_suffix(".path")
        if auto_path.is_file():
            path_filepath = auto_path

    if path_filepath is not None:
        p_path = Path(path_filepath)
        if p_path.is_file():
            with open(p_path, "r", encoding="utf-8", errors="replace") as f:
                path_content = f.read()

            fwd_match = re.search(r"activation energy \(forward\)\s*=\s*([-+]?\d+\.\d+)\s*ev", path_content, re.IGNORECASE)
            if fwd_match:
                fwd_barrier = float(fwd_match.group(1))

            rev_match = re.search(r"activation energy \(reverse\)\s*=\s*([-+]?\d+\.\d+)\s*ev", path_content, re.IGNORECASE)
            if rev_match:
                rev_barrier = float(rev_match.group(1))

            if "converged" in path_content.lower() or "convergence has been achieved" in path_content.lower():
                is_converged = True

            c_match = re.search(r"climbing image\s*=\s*(\d+)", path_content, re.IGNORECASE)
            if c_match:
                climbing_idx = int(c_match.group(1))

    # Analytical fallback if not in .path file: barrier from peak discrete image
    if fwd_barrier is None:
        max_idx = int(np.argmax(rel_energies))
        fwd_barrier = float(rel_energies[max_idx])
        rev_barrier = float(rel_energies[max_idx] - rel_energies[-1])
        climbing_idx = max_idx

    return QENebOutput(
        filepath=p_dat.resolve(),
        file_hash=file_hash,
        image_indices=img_indices,
        image_coordinates=img_coords,
        image_energies_ev=rel_energies,
        raw_energies_ev=raw_energies,
        interpolated_coordinates=interpolated_coords,
        interpolated_energies_ev=interpolated_energies,
        forward_barrier_ev=fwd_barrier,
        reverse_barrier_ev=rev_barrier,
        reaction_energy_ev=reaction_energy,
        is_converged=is_converged,
        climbing_image_index=climbing_idx,
    )

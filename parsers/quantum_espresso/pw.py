"""
Quantum ESPRESSO pw.x Output Parser.

Parses SCF, NSCF, relax, and vc-relax calculation outputs (.out or .pwo).
Extracts total energies, ionic step history, forces, stress, Fermi energy,
convergence diagnostics, and final atomic structures using ASE and native regex.
Never interprets an incomplete or unconverged calculation as successful.
"""

from dataclasses import dataclass, field
import hashlib
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

try:
    from ase import Atoms
    from ase.io.espresso import read_espresso_out
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False

from scripts.preprocess import (
    BOHR_TO_ANGSTROM,
    RY_PER_BOHR_TO_EV_PER_ANGSTROM,
    RY_TO_EV,
)


class QEParseError(Exception):
    """Base exception for Quantum ESPRESSO output parsing failures."""
    pass


class QEIncompleteCalculationError(QEParseError):
    """Raised when a Quantum ESPRESSO calculation terminated abnormally or did not finish."""
    pass


class QEUnconvergedCalculationError(QEParseError):
    """Raised when an electronic SCF or ionic relaxation calculation failed to converge."""
    pass


@dataclass
class QEPwOutput:
    """Structured container for parsed pw.x calculation results."""
    filepath: Path
    file_hash: str
    calculation_type: str
    is_job_done: bool
    is_scf_converged: bool
    is_relaxation_converged: bool
    total_energy_ry: Optional[float]
    total_energy_ev: Optional[float]
    energy_history_ry: List[float] = field(default_factory=list)
    energy_history_ev: List[float] = field(default_factory=list)
    fermi_energy: Optional[float] = None
    highest_occupied: Optional[float] = None
    is_metal: bool = False
    natoms: int = 0
    symbols: List[str] = field(default_factory=list)
    forces_ev_angstrom: Optional[np.ndarray] = None
    forces_ry_bohr: Optional[np.ndarray] = None
    stress_kbar: Optional[np.ndarray] = None
    atomic_positions: Optional[np.ndarray] = None
    cell: Optional[np.ndarray] = None
    ase_atoms: Optional[Any] = None
    raw_metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_valid_final_result(self) -> bool:
        """Returns True only if the job completed and all convergence criteria were met."""
        if not self.is_job_done or not self.is_scf_converged:
            return False
        if self.calculation_type in ("relax", "vc-relax") and not self.is_relaxation_converged:
            return False
        return self.total_energy_ev is not None


def _compute_sha256(filepath: Path) -> str:
    """Compute SHA-256 checksum of a file without loading entire file in memory at once."""
    hasher = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def parse_pw_output(
    filepath: Union[str, Path],
    strict: bool = True,
    use_ase: bool = True,
) -> QEPwOutput:
    """
    Parse Quantum ESPRESSO pw.x output file.

    Parameters
    ----------
    filepath : str or Path
        Path to the pw.x .out or .pwo file.
    strict : bool, default True
        If True, raises QEIncompleteCalculationError if calculation did not finish (JOB DONE.),
        or QEUnconvergedCalculationError if SCF or BFGS relaxation did not converge.
        If False, returns QEPwOutput with convergence flags set for user inspection.
    use_ase : bool, default True
        Whether to attempt loading structure using ASE's read_espresso_out.

    Returns
    -------
    QEPwOutput
        Structured calculation results and convergence diagnostics.

    Raises
    ------
    FileNotFoundError
        If file does not exist.
    QEIncompleteCalculationError
        If calculation terminated prematurely and strict=True.
    QEUnconvergedCalculationError
        If calculation did not achieve convergence and strict=True.
    QEParseError
        If file format cannot be recognized as a valid QE output.
    """
    path = Path(filepath)
    if not path.is_file():
        raise FileNotFoundError(f"Quantum ESPRESSO output file not found: {path.resolve()}")

    file_hash = _compute_sha256(path)

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()

    if "Program PWSCF" not in content and "PWSCF" not in content and "Quantum ESPRESSO" not in content:
        raise QEParseError(f"File {path.name} does not appear to be a Quantum ESPRESSO pw.x output file.")

    # 1. Job completion status
    is_job_done = bool(re.search(r"JOB DONE\.", content, re.IGNORECASE))

    # 2. Calculation type detection
    calc_type = "scf"
    calc_match = re.search(r"calculation\s*=\s*['\"]?([a-zA-Z0-9_\-]+)['\"]?", content, re.IGNORECASE)
    if calc_match:
        calc_type = calc_match.group(1).lower()
    else:
        # Heuristic detection
        if "End of BFGS Geometry Optimization" in content or "Starting BFGS" in content:
            calc_type = "relax"

    # 3. SCF Convergence status
    # Check for failure flags
    has_unconverged_scf = bool(
        re.search(r"convergence NOT achieved", content, re.IGNORECASE)
        or re.search(r"c_bands:.*not converged", content, re.IGNORECASE)
    )
    has_converged_scf = bool(re.search(r"convergence has been achieved", content, re.IGNORECASE))
    is_scf_converged = has_converged_scf and not has_unconverged_scf

    # 4. Relaxation Convergence status
    is_relax_converged = True
    if calc_type in ("relax", "vc-relax"):
        is_relax_converged = bool(
            re.search(r"End of BFGS Geometry Optimization", content)
            or re.search(r"bfgs converged in", content, re.IGNORECASE)
        )
        if re.search(r"The maximum number of steps has been reached", content, re.IGNORECASE):
            is_relax_converged = False

    # 5. Extract Total Energies and Ionic History
    # QE outputs energy lines: '!    total energy              =     -120.45678901 Ry'
    energy_matches = re.findall(r"!\s+total energy\s+=\s+([-+]?\d+\.\d+)\s+Ry", content)
    energy_history_ry = [float(e) for e in energy_matches]
    energy_history_ev = [float(e) * RY_TO_EV for e in energy_history_ry]

    total_energy_ry = energy_history_ry[-1] if energy_history_ry else None
    total_energy_ev = energy_history_ev[-1] if energy_history_ev else None

    # 6. Extract Fermi energy / Band edge
    fermi_energy: Optional[float] = None
    is_metal = False
    fermi_match = re.search(r"the Fermi energy is\s+([-+]?\d+\.\d+)\s+ev", content, re.IGNORECASE)
    if fermi_match:
        fermi_energy = float(fermi_match.group(1))
        is_metal = True

    highest_occupied: Optional[float] = None
    ho_match = re.search(
        r"highest occupied level is\s+([-+]?\d+\.\d+)\s+ev",
        content,
        re.IGNORECASE,
    )
    if ho_match:
        highest_occupied = float(ho_match.group(1))
    else:
        ho_lu_match = re.search(
            r"highest occupied, lowest unoccupied level \(ev\):\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)",
            content,
            re.IGNORECASE,
        )
        if ho_lu_match:
            highest_occupied = float(ho_lu_match.group(1))

    if fermi_energy is None and highest_occupied is not None:
        fermi_energy = highest_occupied

    # 7. Extract Forces (last occurrence)
    forces_ev_angstrom: Optional[np.ndarray] = None
    forces_ry_bohr: Optional[np.ndarray] = None
    forces_block = re.findall(
        r"Forces acting on atoms \(cartesian axes, Ry/au\):\s*\n((?:\s*atom\s+\d+\s+type\s+\d+\s+force\s*=\s*[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)",
        content,
    )
    if forces_block:
        last_block = forces_block[-1]
        force_lines = re.findall(r"force\s*=\s*([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)", last_block)
        if force_lines:
            f_ry = np.array(force_lines, dtype=float)
            forces_ry_bohr = f_ry
            forces_ev_angstrom = f_ry * RY_PER_BOHR_TO_EV_PER_ANGSTROM

    # 8. Extract Stress (last occurrence)
    stress_kbar: Optional[np.ndarray] = None
    stress_block = re.findall(
        r"total\s+stress\s+\(Ry/bohr\*\*3\)\s+\(kbar\)\s+P=\s*[-+]?\d+\.\d+\s*\n((?:\s*[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n){3})",
        content,
    )
    if stress_block:
        last_stress = stress_block[-1]
        stress_vals = []
        for line in last_stress.strip().split("\n"):
            tokens = [float(t) for t in line.split()]
            if len(tokens) >= 6:
                # tokens 3, 4, 5 are the kbar values
                stress_vals.append(tokens[3:6])
        if len(stress_vals) == 3:
            stress_kbar = np.array(stress_vals, dtype=float)

    # 9. Extract Structure using ASE if available
    ase_atoms = None
    symbols: List[str] = []
    positions: Optional[np.ndarray] = None
    cell: Optional[np.ndarray] = None
    natoms = 0

    if use_ase and ASE_AVAILABLE:
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                atoms_gen = read_espresso_out(f, index=-1)
                if isinstance(atoms_gen, Atoms):
                    ase_atoms = atoms_gen
                    symbols = ase_atoms.get_chemical_symbols()
                    positions = ase_atoms.get_positions()
                    cell = np.array(ase_atoms.get_cell())
                    natoms = len(ase_atoms)
        except Exception:
            # Fall back to regex parsing of atoms if ASE encounters unexpected formatting
            ase_atoms = None

    if not symbols:
        # Fallback regex for number of atoms and positions
        nat_match = re.search(r"number of atoms/cell\s*=\s*(\d+)", content, re.IGNORECASE)
        if nat_match:
            natoms = int(nat_match.group(1))

    # 10. Strict Verification Guard
    if strict:
        if not is_job_done:
            raise QEIncompleteCalculationError(
                f"Calculation in {path.name} did not finish cleanly. 'JOB DONE.' not found. "
                "Never interpret incomplete calculations as valid scientific results."
            )
        if not is_scf_converged:
            raise QEUnconvergedCalculationError(
                f"Electronic SCF in {path.name} failed to reach convergence. "
                "Calculated energies and electronic structures are unphysical."
            )
        if calc_type in ("relax", "vc-relax") and not is_relax_converged:
            raise QEUnconvergedCalculationError(
                f"Ionic relaxation in {path.name} did not converge within the maximum iterations. "
                "Final geometry and forces have not satisfied tolerance criteria."
            )

    return QEPwOutput(
        filepath=path.resolve(),
        file_hash=file_hash,
        calculation_type=calc_type,
        is_job_done=is_job_done,
        is_scf_converged=is_scf_converged,
        is_relaxation_converged=is_relax_converged,
        total_energy_ry=total_energy_ry,
        total_energy_ev=total_energy_ev,
        energy_history_ry=energy_history_ry,
        energy_history_ev=energy_history_ev,
        fermi_energy=fermi_energy,
        highest_occupied=highest_occupied,
        is_metal=is_metal,
        natoms=natoms,
        symbols=symbols,
        forces_ev_angstrom=forces_ev_angstrom,
        forces_ry_bohr=forces_ry_bohr,
        stress_kbar=stress_kbar,
        atomic_positions=positions,
        cell=cell,
        ase_atoms=ase_atoms,
        raw_metadata={
            "qe_version": re.search(r"Program PWSCF\s+v\.?([0-9a-zA-Z\.\-]+)", content)
            and re.search(r"Program PWSCF\s+v\.?([0-9a-zA-Z\.\-]+)", content).group(1),  # type: ignore
        },
    )

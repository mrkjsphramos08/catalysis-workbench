"""
Quantum ESPRESSO neb.x Nudged Elastic Band Parser.

Parses NEB reaction profile files: prefix.dat (discrete images),
prefix.int (interpolated splined pathway), and stdout output logs
(forward and reverse activation barriers, convergence status).
Does not accept input configuration files as output profiles.
"""

from dataclasses import dataclass, field
import hashlib
from pathlib import Path
import re
from typing import Optional, Union

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


def _read_profile(path: Path) -> np.ndarray:
    """Read a finite, rectangular QE numeric table, including Fortran exponents."""
    if path.suffix.lower() == ".path":
        raise QEParseError("QE .path restart files are unsupported; supply .dat/.int profiles or a .out log.")
    try:
        text = path.read_text(encoding="utf-8")
        rows = [line.strip() for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
        values = np.array([[float(token.replace("D", "E").replace("d", "e")) for token in row.split()] for row in rows])
    except (ValueError, UnicodeError) as exc:
        raise QEParseError(f"Invalid numeric NEB profile {path.name}: {exc}") from exc
    if values.ndim != 2 or len(values) < 2 or not np.isfinite(values).all():
        raise QEParseError(f"NEB profile {path.name} needs at least two finite rows.")
    return values


def parse_neb_output(
    dat_filepath: Union[str, Path],
    int_filepath: Optional[Union[str, Path]] = None,
    path_filepath: Optional[Union[str, Path]] = None,
    *,
    log_filepath: Optional[Union[str, Path]] = None,
    data_format: str = "auto",
) -> QENebOutput:
    """Parse QE coordinate/relative-energy/error tables and optional stdout logs.

    Native three-column .dat files contain coordinate, relative energy (eV),
    and residual error. Two-column coordinate/energy profiles remain supported.
    Legacy indexed tables require data_format="legacy_indexed" or an explicit
    '# Image Coord Energy...' header. No numeric layout guessing is performed.

    log_filepath accepts a NEB stdout .out log; path_filepath is retained as an
    alias for that argument, but .path restart files are explicitly rejected.
    Logs are never auto-discovered. Convergence requires a final positive NEB
    convergence marker and JOB DONE in the last run. Without a log, barriers
    are discrete-image estimates and no climbing image/convergence is inferred.
    """
    p_dat = Path(dat_filepath)
    table = _read_profile(p_dat)
    header = p_dat.read_text(encoding="utf-8")
    if data_format not in ("auto", "qe", "legacy_indexed"):
        raise ValueError("data_format must be auto, qe, or legacy_indexed")
    legacy_header = bool(re.search(r"(?im)^\s*#\s*image\s+coord\S*\s+energy", header))
    legacy = data_format == "legacy_indexed" or (data_format == "auto" and legacy_header)
    if legacy:
        if table.shape[1] not in (3, 4) or not np.equal(table[:, 0], np.floor(table[:, 0])).all():
            raise QEParseError("Legacy NEB tables require integer index, coordinate, energy, and optional error.")
        indices = table[:, 0].astype(int)
        coords, energies = table[:, 1], table[:, 2]
    else:
        if table.shape[1] not in (2, 3):
            raise QEParseError("QE NEB profiles require coordinate, energy, and optional residual error columns.")
        indices = np.arange(1, len(table) + 1)
        coords, energies = table[:, 0], table[:, 1]
    if np.any(np.diff(coords) <= 0) or np.any(np.diff(indices) <= 0):
        raise QEParseError("NEB coordinates and image indices must strictly increase.")
    rel = energies - energies[0]
    interp_x = interp_e = None
    p_int = Path(int_filepath) if int_filepath is not None else p_dat.with_suffix(".int")
    if int_filepath is not None or p_int.is_file():
        interpolation = _read_profile(p_int)
        if interpolation.shape[1] != 2 or np.any(np.diff(interpolation[:, 0]) <= 0):
            raise QEParseError("NEB interpolation requires increasing coordinate/energy pairs.")
        interp_x = interpolation[:, 0]
        interp_e = interpolation[:, 1] - energies[0]
    if path_filepath is not None and log_filepath is not None:
        raise ValueError("Supply log_filepath or path_filepath, not both.")
    log_path = log_filepath if log_filepath is not None else path_filepath
    fwd, rev = float(rel.max()), float(rel.max() - rel[-1])
    converged, climbing = False, None
    metadata = {"data_format": "legacy_indexed" if legacy else "qe", "coordinate_unit": "arbitrary",
                "energy_unit": "eV", "barrier_source": "discrete_images", "source_files": [str(p_dat.resolve())]}
    if interp_x is not None:
        metadata["source_files"].append(str(p_int.resolve()))
    if log_path is not None:
        p_log = Path(log_path)
        if p_log.suffix.lower() == ".path":
            raise QEParseError("QE .path restart files are unsupported; use log_filepath with a NEB stdout .out log.")
        log = p_log.read_text(encoding="utf-8")
        if "RESTART INFORMATION" in log.upper():
            raise QEParseError("NEB restart content is unsupported, regardless of filename.")
        starts = list(re.finditer(r"Program\s+NEB\b", log, re.IGNORECASE))
        if starts:
            log = log[starts[-1].start():]
        else:
            raise QEParseError("Expected a NEB stdout log beginning with Program NEB.")
        number = r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?)"
        for direction, is_forward in [(r"(?:->|forward)", True), (r"(?:<-|reverse)", False)]:
            matches = re.findall(r"activation energy\s*\(" + direction + r"\)\s*=\s*" + number + r"\s*eV", log, re.IGNORECASE)
            if matches:
                value = float(matches[-1].replace("D", "E").replace("d", "e"))
                if not np.isfinite(value):
                    raise QEParseError("Nonfinite activation energy in NEB log.")
                if is_forward:
                    fwd = value
                else:
                    rev = value
                metadata["barrier_source"] = "stdout_log"
        positive = list(re.finditer(r"(?im)^\s*(?:neb:\s*convergence achieved\b.*|convergence has been achieved\s*)$", log))
        negative = list(re.finditer(r"(?im)^.*(?:not converged|not achieved|convergence NOT achieved|stopping|Error in routine).*$", log))
        done = log.rfind("JOB DONE.")
        converged = bool(positive and done > positive[-1].end() and (not negative or negative[-1].start() < positive[-1].start()))
        climbs = re.findall(r"climbing image\s*=\s*(\d+)", log, re.IGNORECASE)
        if climbs:
            climbing = int(climbs[-1])
            if climbing not in indices:
                raise QEParseError("Climbing image from log is outside the profile image indices.")
        metadata["source_files"].append(str(p_log.resolve()))
    return QENebOutput(
        filepath=p_dat.resolve(), file_hash=hashlib.sha256(p_dat.read_bytes()).hexdigest(),
        image_indices=indices, image_coordinates=coords, image_energies_ev=rel,
        raw_energies_ev=energies, interpolated_coordinates=interp_x, interpolated_energies_ev=interp_e,
        forward_barrier_ev=fwd, reverse_barrier_ev=rev, reaction_energy_ev=float(rel[-1]),
        is_converged=converged, climbing_image_index=climbing, raw_metadata=metadata,
    )

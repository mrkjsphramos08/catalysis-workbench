"""
Quantum ESPRESSO dos.x Output Parser.

Parses native QE DOS files (.dos or .dos.dat).
Supports non-spin-polarized and collinear spin-polarized states,
integrated density of states (IDOS), and Fermi-level alignment.
"""

from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Optional, Tuple, Union

import numpy as np
import pandas as pd

from parsers.quantum_espresso.pw import QEParseError


@dataclass
class QEDosOutput:
    """Structured container for parsed dos.x results."""
    filepath: Path
    file_hash: str
    is_spin_polarized: bool
    energies: np.ndarray
    dos_total: np.ndarray
    dos_up: Optional[np.ndarray] = None
    dos_down: Optional[np.ndarray] = None
    idos: Optional[np.ndarray] = None
    e_fermi: Optional[float] = None
    energies_relative: Optional[np.ndarray] = None
    dataframe: Optional[pd.DataFrame] = None

    def get_aligned_energies(self, fermi_energy: Optional[float] = None) -> np.ndarray:
        """Return energy array relative to Fermi level (E - E_F)."""
        ef = fermi_energy if fermi_energy is not None else self.e_fermi
        if ef is not None:
            return self.energies - float(ef)
        return self.energies.copy()

    def filter_energy_range(
        self,
        energy_min: float,
        energy_max: float,
        relative_to_fermi: bool = True,
    ) -> "QEDosOutput":
        """Filter dataset within energy limits."""
        ref_e = self.get_aligned_energies() if relative_to_fermi else self.energies
        mask = (ref_e >= energy_min) & (ref_e <= energy_max)
        
        filtered_df = self.dataframe.loc[mask].reset_index(drop=True) if self.dataframe is not None else None
        return QEDosOutput(
            filepath=self.filepath,
            file_hash=self.file_hash,
            is_spin_polarized=self.is_spin_polarized,
            energies=self.energies[mask],
            dos_total=self.dos_total[mask],
            dos_up=self.dos_up[mask] if self.dos_up is not None else None,
            dos_down=self.dos_down[mask] if self.dos_down is not None else None,
            idos=self.idos[mask] if self.idos is not None else None,
            e_fermi=self.e_fermi,
            energies_relative=self.energies_relative[mask] if self.energies_relative is not None else None,
            dataframe=filtered_df,
        )


def parse_dos_file(
    filepath: Union[str, Path],
    e_fermi: Optional[float] = None,
) -> QEDosOutput:
    """
    Parse Quantum ESPRESSO dos.x output file.

    Parameters
    ----------
    filepath : str or Path
        Path to native QE DOS file.
    e_fermi : float, optional
        Fermi energy in eV from pw.x calculation.

    Returns
    -------
    QEDosOutput
        Parsed DOS object with total DOS, spin channels, and DataFrame.
    """
    path = Path(filepath)
    if not path.is_file():
        raise FileNotFoundError(f"DOS file not found: {path.resolve()}")

    hasher = hashlib.sha256()
    with open(path, "rb") as f:
        hasher.update(f.read())
    file_hash = hasher.hexdigest()

    # Read data skipping commented header lines starting with #
    try:
        data = np.loadtxt(path, comments="#")
    except Exception as e:
        raise QEParseError(f"Failed to parse numerical data in {path.name}: {e}")

    if data.ndim != 2 or data.shape[0] == 0:
        raise QEParseError(f"DOS file {path.name} contains empty or malformed table.")

    ncols = data.shape[1]
    energies = data[:, 0]

    # Inspect first comment line for header confirmation if present
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        first_line = f.readline().strip()

    is_spin = False
    if ncols == 3:
        # Non-spin-polarized: E (eV), dos(E), Idos(E)
        is_spin = False
        dos_total = data[:, 1]
        dos_up = None
        dos_down = None
        idos = data[:, 2]
        df_cols = ["energy_ev", "dos_total", "idos"]
        df_data = [energies, dos_total, idos]
    elif ncols == 4:
        # Spin-polarized: E (eV), dosup(E), dosdw(E), Idos(E)
        is_spin = True
        dos_up = data[:, 1]
        dos_down = data[:, 2]
        dos_total = dos_up + dos_down
        idos = data[:, 3]
        df_cols = ["energy_ev", "dos_up", "dos_down", "dos_total", "idos"]
        df_data = [energies, dos_up, dos_down, dos_total, idos]
    else:
        raise QEParseError(
            f"Unsupported DOS column count ({ncols}) in {path.name}. "
            "Expected 3 (non-spin) or 4 (collinear spin-polarized) columns."
        )

    # Relative energy calculation
    rel_energies = energies - float(e_fermi) if e_fermi is not None else None
    if rel_energies is not None:
        df_cols.insert(1, "energy_rel_ev")
        df_data.insert(1, rel_energies)

    df = pd.DataFrame(dict(zip(df_cols, df_data)))

    return QEDosOutput(
        filepath=path.resolve(),
        file_hash=file_hash,
        is_spin_polarized=is_spin,
        energies=energies,
        dos_total=dos_total,
        dos_up=dos_up,
        dos_down=dos_down,
        idos=idos,
        e_fermi=e_fermi,
        energies_relative=rel_energies,
        dataframe=df,
    )

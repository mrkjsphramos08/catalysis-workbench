"""
Quantum ESPRESSO projwfc.x Projected DOS (PDOS) Parser.

Discovers and parses prefix.pdos_tot and prefix.pdos_atm#N(X)_wfc#M(orb) files.
Supports element- and atom-resolved projections, orbital decomposition (s, p, d, f),
component selections (pz, px, py; dz2, dzx, dzy, dx2-y2, dxy; eg, t2g), collinear spin,
and prevents double counting of LDOS and sub-orbital components.
Explicitly rejects noncollinear and spin-orbit coupling formats.
"""

from dataclasses import dataclass, field
import glob
import hashlib
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid

from parsers.quantum_espresso.pw import QEParseError


class QENoncollinearUnsupportedError(QEParseError):
    """Raised when noncollinear or spin-orbit coupling PDOS outputs are encountered."""
    pass


# Quantum ESPRESSO orbital component naming in exact column order
COMPONENTS_ORDER_P = ["pz", "px", "py"]
COMPONENTS_ORDER_D = ["dz2", "dzx", "dzy", "dx2-y2", "dxy"]
COMPONENTS_ORDER_F = ["fz3", "fxz2", "fyz2", "fzx2-zy2", "fxyz", "fx3-3xy2", "f3yx2-y3"]


@dataclass
class QEWfcPdos:
    """Parsed projected DOS for a single atomic wavefunction."""
    filepath: Path
    atom_index: int
    element: str
    wfc_index: int
    orbital: str
    is_spin_polarized: bool
    energies: np.ndarray
    # LDOS total for this wfc
    ldos_up: np.ndarray
    ldos_down: Optional[np.ndarray] = None
    # Individual sub-orbital components
    components: Dict[str, Dict[str, np.ndarray]] = field(default_factory=dict)


class QEPdosCollection:
    """
    Collection of all PDOS files belonging to a Quantum ESPRESSO calculation prefix.
    Provides verified aggregation by element, atom, orbital, and component without double counting.
    """

    def __init__(
        self,
        directory: Union[str, Path],
        prefix: Optional[str] = None,
        e_fermi: Optional[float] = None,
    ):
        self.directory = Path(directory).resolve()
        if not self.directory.is_dir():
            raise NotADirectoryError(f"Directory not found: {self.directory}")

        self.prefix = prefix
        self.e_fermi = e_fermi
        self.wfc_list: List[QEWfcPdos] = []
        self.total_dos_up: Optional[np.ndarray] = None
        self.total_dos_down: Optional[np.ndarray] = None
        self.energies: Optional[np.ndarray] = None
        self.is_spin_polarized: bool = False

        self._discover_and_parse()

    def _discover_and_parse(self) -> None:
        """Scan directory and parse all matching PDOS files."""
        # Find prefix if not supplied
        all_pdos_files = list(self.directory.glob("*.pdos_atm*"))
        if not all_pdos_files:
            raise FileNotFoundError(
                f"No Quantum ESPRESSO PDOS files (*.pdos_atm*) found in {self.directory}."
            )

        if self.prefix is None:
            # Infer prefix from first matching file (e.g. prefix.pdos_atm#1(Pt)_wfc#1(d))
            first_name = all_pdos_files[0].name
            self.prefix = first_name.split(".pdos_atm")[0]

        # Check total pdos if present (prefix.pdos_tot)
        tot_file = self.directory / f"{self.prefix}.pdos_tot"
        if tot_file.is_file():
            self._parse_pdos_tot(tot_file)

        # Pattern: prefix.pdos_atm#<atom>(<elem>)_wfc#<wfc>(<orbital>)
        pattern = re.compile(
            rf"^{re.escape(self.prefix)}\.pdos_atm#(\d+)\(([a-zA-Z]+)\)_wfc#(\d+)\(([spdf])\)$"
        )

        for p_file in sorted(self.directory.glob(f"{self.prefix}.pdos_atm*")):
            m = pattern.match(p_file.name)
            if m:
                atom_idx = int(m.group(1))
                elem = m.group(2)
                wfc_idx = int(m.group(3))
                orbital = m.group(4).lower()

                parsed_wfc = self._parse_single_wfc(p_file, atom_idx, elem, wfc_idx, orbital)
                self.wfc_list.append(parsed_wfc)

        if not self.wfc_list:
            raise QEParseError(f"Found no valid wavefunction PDOS files matching prefix '{self.prefix}'.")

        # Establish common energy grid and spin status
        self.energies = self.wfc_list[0].energies
        self.is_spin_polarized = self.wfc_list[0].is_spin_polarized

    def _parse_pdos_tot(self, filepath: Path) -> None:
        """Parse prefix.pdos_tot file."""
        try:
            data = np.loadtxt(filepath, comments="#")
        except Exception:
            return

        if data.ndim == 2:
            ncols = data.shape[1]
            if ncols == 3:
                # E, dos(E), pdos(E)
                self.total_dos_up = data[:, 2]
                self.is_spin_polarized = False
            elif ncols == 5:
                # E, dosup, dosdw, pdosup, pdosdw
                self.total_dos_up = data[:, 3]
                self.total_dos_down = data[:, 4]
                self.is_spin_polarized = True

    def _parse_single_wfc(
        self,
        filepath: Path,
        atom_idx: int,
        elem: str,
        wfc_idx: int,
        orbital: str,
    ) -> QEWfcPdos:
        """Parse an individual wavefunction PDOS file with column verification."""
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            header = f.readline()

        # Reject noncollinear / spin-orbit calculations explicitly
        if "j=" in header.lower() or "tot" in header.lower() or "j_" in header.lower():
            raise QENoncollinearUnsupportedError(
                f"File {filepath.name} appears to be from a noncollinear or spin-orbit coupling calculation. "
                "Only collinear and non-spin-polarized QE PDOS formats are supported."
            )

        try:
            data = np.loadtxt(filepath, comments="#")
        except Exception as e:
            raise QEParseError(f"Error reading numerical data in {filepath.name}: {e}")

        if data.ndim != 2:
            raise QEParseError(f"File {filepath.name} contains invalid dimensional data.")

        energies = data[:, 0]
        ncols = data.shape[1]

        # Determine column layout based on orbital and ncols
        components: Dict[str, Dict[str, np.ndarray]] = {}

        if orbital == "s":
            if ncols == 3:
                # E, ldos, s_pdos
                is_spin = False
                ldos_up = data[:, 1]
                ldos_down = None
                components["s"] = {"up": data[:, 2], "down": None}
            elif ncols == 5:
                # E, ldosup, ldosdw, s_pdosup, s_pdosdw
                is_spin = True
                ldos_up = data[:, 1]
                ldos_down = data[:, 2]
                components["s"] = {"up": data[:, 3], "down": data[:, 4]}
            else:
                raise QEParseError(f"Unexpected column count {ncols} for s orbital in {filepath.name}")

        elif orbital == "p":
            if ncols == 5:
                # E, ldos, pz, px, py
                is_spin = False
                ldos_up = data[:, 1]
                ldos_down = None
                for i, comp in enumerate(COMPONENTS_ORDER_P):
                    components[comp] = {"up": data[:, 2 + i], "down": None}
            elif ncols == 9:
                # E, ldosup, ldosdw, pz_up, pz_dw, px_up, px_dw, py_up, py_dw
                is_spin = True
                ldos_up = data[:, 1]
                ldos_down = data[:, 2]
                for i, comp in enumerate(COMPONENTS_ORDER_P):
                    col_up = 3 + 2 * i
                    col_dw = 4 + 2 * i
                    components[comp] = {"up": data[:, col_up], "down": data[:, col_dw]}
            else:
                raise QEParseError(f"Unexpected column count {ncols} for p orbital in {filepath.name}")

        elif orbital == "d":
            if ncols == 7:
                # E, ldos, dz2, dzx, dzy, dx2-y2, dxy
                is_spin = False
                ldos_up = data[:, 1]
                ldos_down = None
                for i, comp in enumerate(COMPONENTS_ORDER_D):
                    components[comp] = {"up": data[:, 2 + i], "down": None}
            elif ncols == 13:
                # E, ldosup, ldosdw, dz2_up, dz2_dw, ...
                is_spin = True
                ldos_up = data[:, 1]
                ldos_down = data[:, 2]
                for i, comp in enumerate(COMPONENTS_ORDER_D):
                    col_up = 3 + 2 * i
                    col_dw = 4 + 2 * i
                    components[comp] = {"up": data[:, col_up], "down": data[:, col_dw]}
            else:
                raise QEParseError(f"Unexpected column count {ncols} for d orbital in {filepath.name}")

        elif orbital == "f":
            if ncols == 9:
                # Non-spin f (7 components)
                is_spin = False
                ldos_up = data[:, 1]
                ldos_down = None
                for i, comp in enumerate(COMPONENTS_ORDER_F):
                    components[comp] = {"up": data[:, 2 + i], "down": None}
            elif ncols == 17:
                # Spin-polarized f
                is_spin = True
                ldos_up = data[:, 1]
                ldos_down = data[:, 2]
                for i, comp in enumerate(COMPONENTS_ORDER_F):
                    col_up = 3 + 2 * i
                    col_dw = 4 + 2 * i
                    components[comp] = {"up": data[:, col_up], "down": data[:, col_dw]}
            else:
                raise QEParseError(f"Unexpected column count {ncols} for f orbital in {filepath.name}")
        else:
            raise QEParseError(f"Unrecognized orbital label '{orbital}' in {filepath.name}")

        return QEWfcPdos(
            filepath=filepath,
            atom_index=atom_idx,
            element=elem,
            wfc_index=wfc_idx,
            orbital=orbital,
            is_spin_polarized=is_spin,
            energies=energies,
            ldos_up=ldos_up,
            ldos_down=ldos_down,
            components=components,
        )

    def get_projection(
        self,
        element: Optional[str] = None,
        atom_index: Optional[Union[int, Sequence[int]]] = None,
        orbital: Optional[str] = None,
        component: Optional[Union[str, Sequence[str]]] = None,
    ) -> Dict[str, Optional[np.ndarray]]:
        """
        Aggregate PDOS with rigorous exclusion of LDOS double counting.

        Parameters
        ----------
        element : str, optional
            Filter by element symbol (e.g. 'Pt').
        atom_index : int or sequence of int, optional
            Filter by 1-based atom index.
        orbital : str, optional
            Filter by orbital ('s', 'p', 'd', 'f').
        component : str or sequence of str, optional
            Specific sub-orbital components:
            - 'pz', 'px', 'py'
            - 'dz2', 'dzx', 'dzy', 'dx2-y2', 'dxy'
            - 'eg' (sum of dz2 and dx2-y2)
            - 't2g' (sum of dxy, dzx, dzy)

        Returns
        -------
        Dict[str, Optional[np.ndarray]]
            {'up': sum_array_up, 'down': sum_array_down (or None if unpolarized)}
        """
        target_atoms: Optional[Set[int]] = None
        if atom_index is not None:
            target_atoms = {atom_index} if isinstance(atom_index, int) else set(atom_index)

        target_elem = element.strip() if element else None
        target_orb = orbital.strip().lower() if orbital else None

        # Resolve components
        target_comps: Optional[Set[str]] = None
        if component is not None:
            if isinstance(component, str):
                c_low = component.strip().lower()
                if c_low == "eg":
                    target_comps = {"dz2", "dx2-y2"}
                elif c_low == "t2g":
                    target_comps = {"dxy", "dzx", "dzy"}
                else:
                    target_comps = {c_low}
            else:
                target_comps = set()
                for c in component:
                    c_low = c.strip().lower()
                    if c_low == "eg":
                        target_comps.update(["dz2", "dx2-y2"])
                    elif c_low == "t2g":
                        target_comps.update(["dxy", "dzx", "dzy"])
                    else:
                        target_comps.add(c_low)

        sum_up = np.zeros_like(self.energies)
        sum_down = np.zeros_like(self.energies) if self.is_spin_polarized else None
        matched_any = False

        for wfc in self.wfc_list:
            if target_elem and wfc.element.lower() != target_elem.lower():
                continue
            if target_atoms and wfc.atom_index not in target_atoms:
                continue
            if target_orb and wfc.orbital != target_orb:
                continue

            matched_any = True

            if target_comps:
                # Sum only requested sub-orbital components (NO LDOS)
                for comp_name in target_comps:
                    if comp_name in wfc.components:
                        comp_data = wfc.components[comp_name]
                        sum_up += comp_data["up"]
                        if self.is_spin_polarized and sum_down is not None and comp_data["down"] is not None:
                            sum_down += comp_data["down"]
            else:
                # No sub-component specified: sum the components of this wfc
                # To guarantee NO double counting, sum the sub-orbital components!
                for comp_data in wfc.components.values():
                    sum_up += comp_data["up"]
                    if self.is_spin_polarized and sum_down is not None and comp_data["down"] is not None:
                        sum_down += comp_data["down"]

        if not matched_any:
            raise ValueError(
                f"No PDOS wavefunction matched criteria (element={element}, "
                f"atom={atom_index}, orbital={orbital}, component={component})."
            )

        return {"up": sum_up, "down": sum_down}

    def get_projections_for_plotting(
        self,
        specs: Dict[str, Dict[str, Any]],
    ) -> Dict[str, Any]:
        """
        Format projections dictionary directly compatible with modules.dft.plot_pdos.

        Parameters
        ----------
        specs : dict
            Mapping of label -> kwargs for get_projection().
            Example: {
                r'Pt $5d$': {'element': 'Pt', 'orbital': 'd'},
                r'O $2p$': {'element': 'O', 'orbital': 'p'},
            }

        Returns
        -------
        Dict[str, Any]
            Dict formatted for plot_pdos(projections=...).
        """
        out: Dict[str, Any] = {}
        for label, kwargs in specs.items():
            proj = self.get_projection(**kwargs)
            if self.is_spin_polarized:
                out[label] = {"up": proj["up"], "down": proj["down"]}
            else:
                out[label] = proj["up"]
        return out

    def calculate_d_band_center(
        self,
        element: str,
        energy_window: Tuple[float, float] = (-10.0, 0.0),
        fermi_energy: Optional[float] = None,
        spin: str = "total",
    ) -> float:
        """
        Calculate d-band center (epsilon_d) as standard catalytic descriptor:
        epsilon_d = integral( (E - E_F) * rho_d(E) dE ) / integral( rho_d(E) dE )

        Parameters
        ----------
        element : str
            Element symbol (e.g. 'Pt').
        energy_window : tuple of (float, float), default (-10.0, 0.0)
            Integration energy window in eV relative to Fermi level.
        fermi_energy : float, optional
            Fermi energy in eV. If None, uses self.e_fermi.
        spin : str, default 'total'
            'up', 'down', or 'total'.

        Returns
        -------
        float
            d-band center in eV relative to Fermi level.
        """
        ef = fermi_energy if fermi_energy is not None else self.e_fermi
        if ef is None:
            raise ValueError("Fermi energy must be specified to calculate d-band center relative to E_F.")

        proj = self.get_projection(element=element, orbital="d")
        e_rel = self.energies - float(ef)

        mask = (e_rel >= energy_window[0]) & (e_rel <= energy_window[1])
        e_sub = e_rel[mask]

        if spin == "up":
            dos_sub = proj["up"][mask]
        elif spin == "down":
            if proj["down"] is None:
                raise ValueError("Calculation is non-spin-polarized; cannot isolate 'down' channel.")
            dos_sub = proj["down"][mask]
        else:
            dos_sub = proj["up"][mask] + (proj["down"][mask] if proj["down"] is not None else 0.0)

        denom = trapezoid(dos_sub, e_sub)
        if denom == 0.0:
            raise ValueError(f"Integrated d-band DOS is zero in window {energy_window}.")

        numer = trapezoid(e_sub * dos_sub, e_sub)
        return float(numer / denom)

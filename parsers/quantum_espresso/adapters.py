"""
Lightweight Adapters Connecting Quantum ESPRESSO Parsers to Existing Publication Plotting APIs.

Provides direct bridging between QEPwOutput, QEDosOutput, QEPdosCollection,
QEBandsOutput, and QENebOutput objects and the modular visualization functions:
- plot_dos()
- plot_pdos()
- plot_band_structure()
- calculate_adsorption_energy()
- plot_reaction_pathway()
Also handles research data provenance tracking and metadata JSON generation.
"""

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np

from modules.catalysis.reaction_energy import plot_reaction_pathway
from modules.dft.adsorption import calculate_adsorption_energy
from modules.dft.electronic_structure import (
    plot_band_structure,
    plot_dos,
    plot_pdos,
)
from parsers.quantum_espresso.bands import QEBandsOutput
from parsers.quantum_espresso.dos import QEDosOutput
from parsers.quantum_espresso.neb import QENebOutput
from parsers.quantum_espresso.pdos import QEPdosCollection
from parsers.quantum_espresso.pw import QEPwOutput
from scripts.config import DEFAULT_DPI, OUTPUTS_RASTER_DIR, OUTPUTS_VECTOR_DIR
from scripts.utils import create_publication_figure, save_publication_figure


def plot_qe_dos(
    dos_data: QEDosOutput,
    ax: Optional[plt.Axes] = None,
    e_fermi: Optional[float] = None,
    energy_range: Optional[Tuple[float, float]] = (-6.0, 6.0),
    **kwargs: Any,
) -> plt.Axes:
    """
    Directly render publication DOS from parsed QEDosOutput.
    Reuses modules.dft.electronic_structure.plot_dos.
    """
    fermi = e_fermi if e_fermi is not None else dos_data.e_fermi
    return plot_dos(
        energies=dos_data.energies,
        dos_up=dos_data.dos_up if dos_data.is_spin_polarized else dos_data.dos_total,
        dos_down=dos_data.dos_down if dos_data.is_spin_polarized else None,
        e_fermi=fermi,
        ax=ax,
        energy_range=energy_range,
        **kwargs,
    )


def plot_qe_pdos(
    pdos_collection: QEPdosCollection,
    specs: Optional[Dict[str, Dict[str, Any]]] = None,
    ax: Optional[plt.Axes] = None,
    e_fermi: Optional[float] = None,
    energy_range: Optional[Tuple[float, float]] = (-6.0, 6.0),
    **kwargs: Any,
) -> plt.Axes:
    """
    Directly render publication PDOS from parsed QEPdosCollection.
    Reuses modules.dft.electronic_structure.plot_pdos.
    """
    fermi = e_fermi if e_fermi is not None else pdos_collection.e_fermi

    # If no explicit projection specs provided, create default element orbital projections
    if specs is None:
        specs = {}
        # Discover elements and orbitals
        seen: Set[Tuple[str, str]] = set()
        for wfc in pdos_collection.wfc_list:
            pair = (wfc.element, wfc.orbital)
            if pair not in seen:
                seen.add(pair)
                label = f"{wfc.element}-{wfc.orbital}"
                specs[label] = {"element": wfc.element, "orbital": wfc.orbital}

    projections = pdos_collection.get_projections_for_plotting(specs)
    total_dos = pdos_collection.total_dos_up if pdos_collection.total_dos_up is not None else None

    return plot_pdos(
        energies=pdos_collection.energies,
        projections=projections,
        total_dos=total_dos,
        e_fermi=fermi,
        ax=ax,
        energy_range=energy_range,
        **kwargs,
    )


def plot_qe_bands(
    bands_data: QEBandsOutput,
    ax: Optional[plt.Axes] = None,
    e_fermi: Optional[float] = None,
    energy_range: Optional[Tuple[float, float]] = (-5.0, 5.0),
    **kwargs: Any,
) -> plt.Axes:
    """
    Directly render publication band structure from parsed QEBandsOutput.
    Reuses modules.dft.electronic_structure.plot_band_structure.
    """
    fermi = e_fermi if e_fermi is not None else bands_data.e_fermi
    return plot_band_structure(
        k_distances=bands_data.k_distances,
        bands=bands_data.bands,
        high_symmetry_points=bands_data.high_symmetry_points,
        e_fermi=fermi,
        bands_down=bands_data.bands_down,
        ax=ax,
        energy_range=energy_range,
        **kwargs,
    )


def calculate_adsorption_energy_from_qe(
    slab_output: QEPwOutput,
    slab_adsorbate_output: QEPwOutput,
    adsorbate_ref_ev: float,
    n_adsorbates: int = 1,
    sign_convention: str = "negative_favorable",
) -> float:
    """
    Calculate adsorption energy directly from parsed clean slab and adsorbate pw.x outputs.
    Guarantees both calculations are converged and completed.
    """
    if not slab_output.is_valid_final_result:
        raise ValueError(f"Clean slab calculation {slab_output.filepath.name} is incomplete or unconverged.")
    if not slab_adsorbate_output.is_valid_final_result:
        raise ValueError(
            f"Slab+adsorbate calculation {slab_adsorbate_output.filepath.name} is incomplete or unconverged."
        )

    e_clean = slab_output.total_energy_ev
    e_tot = slab_adsorbate_output.total_energy_ev

    return calculate_adsorption_energy(
        e_total=e_tot,
        e_clean_surface=e_clean,
        e_adsorbate_ref=adsorbate_ref_ev,
        n_adsorbates=n_adsorbates,
        sign_convention=sign_convention,
    )


def plot_qe_neb_profile(
    neb_data: QENebOutput,
    ax: Optional[plt.Axes] = None,
    marker_color: str = "#0072B2",
    spline_color: str = "#D55E00",
    title: Optional[str] = None,
) -> plt.Axes:
    """
    Render continuous and discrete reaction profile from parsed QENebOutput.
    """
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")

    # Plot discrete images
    ax.scatter(
        neb_data.image_coordinates,
        neb_data.image_energies_ev,
        color=marker_color,
        s=32,
        edgecolor="#222222",
        linewidth=0.7,
        label="Images",
        zorder=4,
    )

    # Plot interpolated splined path if present
    if (
        neb_data.interpolated_coordinates is not None
        and neb_data.interpolated_energies_ev is not None
    ):
        ax.plot(
            neb_data.interpolated_coordinates,
            neb_data.interpolated_energies_ev,
            color=spline_color,
            linewidth=1.2,
            label="Interpolated Path",
            zorder=3,
        )

    # Label forward and reverse barriers
    if neb_data.forward_barrier_ev is not None:
        ax.annotate(
            rf"$E_{{\mathrm{{a,fwd}}}} = {neb_data.forward_barrier_ev:.2f}\,$eV"
            + (rf", $E_{{\mathrm{{a,rev}}}} = {neb_data.reverse_barrier_ev:.2f}\,$eV" if neb_data.reverse_barrier_ev is not None else ""),
            xy=(0.04, 0.90),
            xycoords="axes fraction",
            fontsize=7.0,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#F8F8F8", edgecolor="#CCCCCC", alpha=0.9),
        )

    ax.set_xlabel(r"Reaction Coordinate ($\mathrm{\AA}$)")
    ax.set_ylabel(r"Relative Energy $\Delta E$ (eV)")
    ax.axhline(0, color="#888888", linestyle="--", linewidth=0.6)

    if title:
        ax.set_title(title)

    ax.legend(loc="lower right", framealpha=0.85, fontsize=7.0)
    return ax


def record_figure_provenance(
    figure_name: str,
    input_files: Sequence[Union[str, Path]],
    output_dirs: Sequence[Union[str, Path]],
    metadata: Optional[Dict[str, Any]] = None,
) -> Path:
    """
    Write a companion provenance JSON metadata file alongside generated publication figures.

    Parameters
    ----------
    figure_name : str
        Base figure name without extension.
    input_files : sequence of paths
        Calculation and experimental input files used to produce figure.
    output_dirs : sequence of paths
        Directories where figure was saved.
    metadata : dict, optional
        Additional analysis parameters, Fermi levels, or physical assumptions.

    Returns
    -------
    Path
        Path to the primary saved metadata JSON file.
    """
    file_records = []
    for f in input_files:
        p = Path(f)
        if p.is_file():
            hasher = hashlib.sha256()
            with open(p, "rb") as fp:
                hasher.update(fp.read())
            file_records.append({
                "path": str(p.resolve()),
                "filename": p.name,
                "sha256": hasher.hexdigest(),
                "size_bytes": p.stat().st_size,
            })
        else:
            file_records.append({"path": str(f), "filename": p.name, "status": "not_found"})

    payload: Dict[str, Any] = {
        "figure_name": figure_name,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_files": file_records,
        "parser_version": "1.2.0-qe",
        "units": {
            "energy_original": "Ry / eV (Quantum ESPRESSO)",
            "energy_plotted": "eV",
            "fermi_aligned": True,
        },
        "analysis_metadata": metadata or {},
    }

    primary_meta_path = None
    for d in output_dirs:
        p_dir = Path(d)
        p_dir.mkdir(parents=True, exist_ok=True)
        meta_file = p_dir / f"{figure_name}_metadata.json"
        with open(meta_file, "w", encoding="utf-8") as fp:
            json.dump(payload, fp, indent=2)
        if primary_meta_path is None:
            primary_meta_path = meta_file

    return primary_meta_path or Path(f"{figure_name}_metadata.json")

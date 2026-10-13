"""
Comprehensive Test Suite for Quantum ESPRESSO Parsers, Scientific Adapters, and Provenance Tracking.
"""

from pathlib import Path
import tempfile
import matplotlib.pyplot as plt
import numpy as np
import pytest

from parsers.quantum_espresso import (
    QEIncompleteCalculationError,
    QENoncollinearUnsupportedError,
    QEParseError,
    QEUnconvergedCalculationError,
    calculate_adsorption_energy_from_qe,
    parse_bands_gnu,
    parse_dos_file,
    parse_neb_output,
    parse_pw_output,
    plot_qe_bands,
    plot_qe_dos,
    plot_qe_neb_profile,
    plot_qe_pdos,
    record_figure_provenance,
)
from parsers.quantum_espresso.pdos import QEPdosCollection
from scripts.preprocess import RY_TO_EV
from scripts.utils import create_publication_figure

FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "qe"


# =========================================================================
# 1. pw.x SCF & Relaxation Tests
# =========================================================================

def test_parse_pw_scf_converged():
    """Verify nonmagnetic converged SCF output parsing and Ry-to-eV conversion."""
    out_file = FIXTURE_DIR / "scf_converged.out"
    res = parse_pw_output(out_file, strict=True)

    assert res.is_job_done is True
    assert res.is_scf_converged is True
    assert res.is_valid_final_result is True
    assert pytest.approx(res.total_energy_ry) == -15.84523190
    assert pytest.approx(res.total_energy_ev) == -15.84523190 * RY_TO_EV

    # Highest occupied level
    assert pytest.approx(res.highest_occupied) == 6.2541
    assert pytest.approx(res.fermi_energy) == 6.2541

    # Forces and stress
    assert res.forces_ev_angstrom is not None
    assert res.forces_ev_angstrom.shape == (2, 3)
    assert res.stress_kbar is not None
    assert res.stress_kbar.shape == (3, 3)
    assert pytest.approx(res.stress_kbar[0, 0]) == -1.25


def test_parse_pw_spin_converged():
    """Verify spin-polarized metallic SCF output parsing."""
    out_file = FIXTURE_DIR / "scf_spin_converged.out"
    res = parse_pw_output(out_file, strict=True)

    assert res.is_job_done is True
    assert res.is_metal is True
    assert pytest.approx(res.fermi_energy) == 12.3456
    assert pytest.approx(res.total_energy_ry) == -62.15243100
    assert len(res.file_hash) == 64


def test_parse_pw_unconverged_scf():
    """Verify that unconverged SCF calculation is rejected in strict mode."""
    out_file = FIXTURE_DIR / "scf_unconverged.out"
    with pytest.raises(QEUnconvergedCalculationError):
        parse_pw_output(out_file, strict=True)

    # In strict=False mode, allows inspection with is_scf_converged=False
    res = parse_pw_output(out_file, strict=False)
    assert res.is_scf_converged is False
    assert res.is_valid_final_result is False


def test_parse_pw_incomplete():
    """Verify that incomplete (crashed/killed) job is rejected."""
    out_file = FIXTURE_DIR / "scf_incomplete.out"
    with pytest.raises(QEIncompleteCalculationError):
        parse_pw_output(out_file, strict=True)

    res = parse_pw_output(out_file, strict=False)
    assert res.is_job_done is False
    assert res.is_valid_final_result is False


def test_parse_pw_relax_converged():
    """Verify relaxation output with ionic step energy history."""
    out_file = FIXTURE_DIR / "relax_converged.out"
    res = parse_pw_output(out_file, strict=True)

    assert res.calculation_type == "relax"
    assert res.is_relaxation_converged is True
    assert len(res.energy_history_ry) == 3
    # Energies should decrease in BFGS optimization
    assert res.energy_history_ry[0] > res.energy_history_ry[1] > res.energy_history_ry[2]
    assert pytest.approx(res.total_energy_ry) == -120.46100000


# =========================================================================
# 2. dos.x Total DOS Tests
# =========================================================================

def test_parse_dos_nonspin():
    """Verify non-spin-polarized dos.x output parsing and Fermi alignment."""
    dos_file = FIXTURE_DIR / "dos_nonspin.dat"
    res = parse_dos_file(dos_file, e_fermi=5.0)

    assert res.is_spin_polarized is False
    assert len(res.energies) == 101
    assert res.dos_up is None
    assert res.dos_total is not None
    assert res.idos is not None

    # Check Fermi alignment (E - E_F)
    aligned_e = res.get_aligned_energies()
    assert pytest.approx(aligned_e[50]) == 0.0  # Center was at 5.0 eV


def test_parse_dos_spin():
    """Verify collinear spin-polarized dos.x parsing and filtering."""
    dos_file = FIXTURE_DIR / "dos_spin.dat"
    res = parse_dos_file(dos_file, e_fermi=5.0)

    assert res.is_spin_polarized is True
    assert res.dos_up is not None
    assert res.dos_down is not None
    assert len(res.dos_up) == len(res.dos_down)

    # Filter energy range
    filtered = res.filter_energy_range(energy_min=-2.0, energy_max=2.0, relative_to_fermi=True)
    assert len(filtered.energies) < len(res.energies)
    aligned_f = filtered.get_aligned_energies()
    assert aligned_f.min() >= -2.0
    assert aligned_f.max() <= 2.0


def test_parse_dos_malformed():
    """Verify rejection of malformed DOS file."""
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tf:
        tf.write("not_a_number not_a_number\n")
        temp_name = tf.name

    try:
        with pytest.raises(QEParseError):
            parse_dos_file(temp_name)
    finally:
        Path(temp_name).unlink(missing_ok=True)


# =========================================================================
# 3. projwfc.x Projected DOS Tests
# =========================================================================

def test_parse_pdos_collection_spin():
    """Verify discovery and projection aggregation for spin-polarized system."""
    coll = QEPdosCollection(FIXTURE_DIR / "pdos_spin", prefix="fe", e_fermi=5.0)
    assert coll.is_spin_polarized is True
    assert len(coll.wfc_list) == 3

    # Projection for Fe d-orbitals
    fe_d = coll.get_projection(element="Fe", orbital="d")
    assert fe_d["up"] is not None
    assert fe_d["down"] is not None
    assert len(fe_d["up"]) == len(coll.energies)

    # Projection for O p-orbitals
    o_p = coll.get_projection(element="O", orbital="p")
    assert o_p["up"] is not None
    assert o_p["down"] is not None


def test_parse_pdos_suborbitals_no_double_counting():
    """Verify individual suborbital component extraction (eg, t2g, pz) without double counting."""
    coll = QEPdosCollection(FIXTURE_DIR / "pdos_spin", prefix="fe", e_fermi=5.0)

    # Total d
    d_total = coll.get_projection(element="Fe", orbital="d")
    # eg = dz2 + dx2-y2
    d_eg = coll.get_projection(element="Fe", orbital="d", component="eg")
    # t2g = dxy + dzx + dzy
    d_t2g = coll.get_projection(element="Fe", orbital="d", component="t2g")

    # Sum of eg and t2g MUST equal total d
    assert np.allclose(d_eg["up"] + d_t2g["up"], d_total["up"], atol=1e-5)
    assert np.allclose(d_eg["down"] + d_t2g["down"], d_total["down"], atol=1e-5)


def test_parse_pdos_nonspin():
    """Verify non-spin PDOS parsing (Si system)."""
    coll = QEPdosCollection(FIXTURE_DIR / "pdos_nonspin", prefix="si", e_fermi=5.0)
    assert coll.is_spin_polarized is False
    si_p = coll.get_projection(element="Si", orbital="p")
    assert si_p["up"] is not None
    assert si_p["down"] is None


def test_parse_pdos_spin_orbit_rejection():
    """Verify rejection of noncollinear / spin-orbit coupling PDOS formats."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        bad_file = tmp_path / "prefix.pdos_atm#1(Pt)_wfc#1(d)"
        bad_file.write_text("# E (eV)  dos(E)  pdos(E)  j=1.5  j=2.5\n0.0 1.0 1.0 0.5 0.5\n", encoding="utf-8")

        with pytest.raises(QENoncollinearUnsupportedError):
            QEPdosCollection(tmp_path, prefix="prefix")


def test_calculate_d_band_center():
    """Verify d-band center descriptor calculation."""
    coll = QEPdosCollection(FIXTURE_DIR / "pdos_spin", prefix="fe", e_fermi=5.0)
    dbc = coll.calculate_d_band_center(element="Fe", energy_window=(-5.0, 0.0), fermi_energy=5.0)
    # Peak of Fe d is around E = 4.5 -> relative to EF (5.0) is -0.5 eV
    assert -1.5 <= dbc <= 0.0


# =========================================================================
# 4. bands.x Band Structure Tests
# =========================================================================

def test_parse_bands_gnu():
    """Verify electronic band structure parsing from .gnu file."""
    bands_file = FIXTURE_DIR / "bands.gnu"
    hs_pts = [(0.0, r"$\Gamma$"), (0.75, "X"), (1.5, "M")]
    res = parse_bands_gnu(bands_file, e_fermi=0.5, high_symmetry_points=hs_pts)

    assert res.n_bands == 3
    assert res.n_kpoints == 31
    assert len(res.high_symmetry_points) == 3

    aligned_b = res.get_aligned_bands()
    assert aligned_b.shape == (3, 31)
    assert pytest.approx(aligned_b[0, 0]) == res.bands[0, 0] - 0.5


# =========================================================================
# 5. neb.x Reaction Pathway Tests
# =========================================================================

def test_parse_neb_output():
    """Verify NEB output parsing, barriers, and splined path."""
    dat_file = FIXTURE_DIR / "neb.dat"
    res = parse_neb_output(dat_file)

    assert len(res.image_indices) == 5
    assert pytest.approx(res.image_energies_ev[0]) == 0.0  # Relative to initial
    assert pytest.approx(res.forward_barrier_ev) == 0.85
    assert pytest.approx(res.reverse_barrier_ev) == 1.50
    assert res.is_converged is True
    assert res.interpolated_coordinates is not None


def test_parse_neb_reject_input():
    """Verify rejection of NEB input file."""
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tf:
        tf.write("&PATH\n string_method = 'neb'\n/\n")
        temp_name = tf.name

    try:
        with pytest.raises(QEParseError):
            parse_neb_output(temp_name)
    finally:
        Path(temp_name).unlink(missing_ok=True)


# =========================================================================
# 6. Adapters & Provenance Tracking Tests
# =========================================================================

def test_adapter_plot_qe_dos():
    """Verify plot_qe_dos adapter connects parsed DOS to existing plotting function."""
    dos_data = parse_dos_file(FIXTURE_DIR / "dos_spin.dat", e_fermi=5.0)
    fig, ax = create_publication_figure()
    ax_out = plot_qe_dos(dos_data, ax=ax, energy_range=(-4.0, 4.0))
    assert ax_out is ax
    assert len(ax.lines) >= 2
    plt.close(fig)


def test_adapter_plot_qe_pdos():
    """Verify plot_qe_pdos adapter connects parsed PDOS collection to plot_pdos."""
    coll = QEPdosCollection(FIXTURE_DIR / "pdos_spin", prefix="fe", e_fermi=5.0)
    fig, ax = create_publication_figure()
    specs = {
        r"Fe $3d$": {"element": "Fe", "orbital": "d"},
        r"O $2p$": {"element": "O", "orbital": "p"},
    }
    ax_out = plot_qe_pdos(coll, specs=specs, ax=ax)
    assert ax_out is ax
    plt.close(fig)


def test_adapter_plot_qe_bands():
    """Verify plot_qe_bands adapter connects parsed bands to plot_band_structure."""
    bands_data = parse_bands_gnu(FIXTURE_DIR / "bands.gnu", e_fermi=0.0)
    fig, ax = create_publication_figure()
    ax_out = plot_qe_bands(bands_data, ax=ax)
    assert ax_out is ax
    plt.close(fig)


def test_adapter_plot_qe_neb():
    """Verify plot_qe_neb_profile adapter plots continuous path and discrete images."""
    neb_data = parse_neb_output(FIXTURE_DIR / "neb.dat")
    fig, ax = create_publication_figure()
    ax_out = plot_qe_neb_profile(neb_data, ax=ax)
    assert ax_out is ax
    plt.close(fig)


def test_adapter_adsorption_energy_from_qe():
    """Verify adsorption energy calculation directly from two pw.x outputs."""
    slab = parse_pw_output(FIXTURE_DIR / "scf_converged.out", strict=True)
    slab_adsorbate = parse_pw_output(FIXTURE_DIR / "scf_spin_converged.out", strict=True)

    # Calculate adsorption energy
    # E_ads = E(tot) - E(slab) - E(ref)
    ref_ev = -10.0
    e_ads = calculate_adsorption_energy_from_qe(slab, slab_adsorbate, adsorbate_ref_ev=ref_ev, n_adsorbates=1)
    expected = slab_adsorbate.total_energy_ev - slab.total_energy_ev - ref_ev
    assert pytest.approx(e_ads) == expected


def test_adapter_provenance_metadata():
    """Verify JSON provenance metadata recording."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        meta_file = record_figure_provenance(
            figure_name="fig_test_qe",
            input_files=[FIXTURE_DIR / "scf_converged.out"],
            output_dirs=[tmp_path],
            metadata={"fermi_energy_ev": 6.2541, "code": "Quantum ESPRESSO v7.2"},
        )
        assert meta_file.is_file()
        import json
        with open(meta_file, "r") as f:
            data = json.load(f)
        assert data["figure_name"] == "fig_test_qe"
        assert len(data["input_files"]) == 1
        assert len(data["input_files"][0]["sha256"]) == 64
        assert data["analysis_metadata"]["fermi_energy_ev"] == 6.2541

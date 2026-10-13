"""Independent regressions for confirmed scientific and export defects."""

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest
from PIL import Image

from parsers.quantum_espresso import parse_neb_output
from parsers.quantum_espresso.pw import QEParseError
from parsers.quantum_espresso.adapters import record_figure_provenance
from parsers.characterization.common import CharacterizationData, InvalidUnitsError
from parsers.characterization.ftir import parse_ftir_file, convert_transmittance_to_absorbance
from parsers.electrochemistry import load_electrochemistry_data
from modules.electrochemistry.reference_electrodes import convert_potential_to_rhe
from scripts.utils import create_publication_figure, save_publication_figure

UPSTREAM = Path(__file__).parent / "fixtures/qe/upstream_neb"


def test_authentic_qe_neb_columns_and_last_iteration():
    result = parse_neb_output(UPSTREAM / "H2+H.dat", log_filepath=UPSTREAM / "H2+H.out")
    assert result.image_coordinates[1] == pytest.approx(0.1939576458, abs=1e-10)
    assert result.image_energies_ev[3] == pytest.approx(0.2042674385, abs=1e-10)
    assert result.forward_barrier_ev == pytest.approx(0.204267, abs=1e-7)
    assert result.reverse_barrier_ev == pytest.approx(0.204268, abs=1e-7)
    assert result.is_converged
    assert result.climbing_image_index == 4  # QE labels images from one.
    assert result.file_hash == hashlib.sha256((UPSTREAM / "H2+H.dat").read_bytes()).hexdigest()
    assert result.interpolated_energies_ev[0] == pytest.approx(0, abs=1e-8)


def test_neb_without_log_does_not_invent_convergence_or_climbing_image():
    result = parse_neb_output(UPSTREAM / "H2+H.dat")
    assert not result.is_converged
    assert result.climbing_image_index is None
    assert result.forward_barrier_ev == pytest.approx(0.2042674385)


def test_neb_rejects_restart_file():
    with pytest.raises(QEParseError, match="restart"):
        parse_neb_output(UPSTREAM / "H2+H.dat", path_filepath=UPSTREAM / "H2+H.path")


@pytest.mark.parametrize("rows", ["0 0 0\n1 nan 0", "0 0 0\n0 1 0", "0 0\n1 1 2", "0 0 0 0\n1 1 1 1"])
def test_neb_rejects_invalid_or_ambiguous_tables(tmp_path, rows):
    path = tmp_path / "bad.dat"
    path.write_text(rows)
    with pytest.raises(QEParseError):
        parse_neb_output(path)


def test_neb_negative_or_truncated_log_not_converged(tmp_path):
    path = tmp_path / "run.out"
    for text in ["Program NEB\nnot converged\nJOB DONE.", "Program NEB\nneb: convergence achieved in 5 iterations"]:
        path.write_text(text)
        assert not parse_neb_output(UPSTREAM / "H2+H.dat", log_filepath=path).is_converged


@pytest.mark.parametrize("instrument,extension,contents", [
    ("metrohm", ".txt", "Time (s)\tWE(1).Potential (V)\tWE(1).Current (A)\n0\t0.1\t0.01\n1\t0.2\t0.02\n"),
    ("biologic", ".mpt", "EC-Lab ASCII FILE\nNb header lines : 3\ntime/s\tEwe/V\tI/mA\n0\t0.1\t10\n1\t0.2\t20\n"),
])
def test_missing_electrochemical_conditions_remain_unknown(tmp_path, instrument, extension, contents):
    path = tmp_path / ("export" + extension)
    path.write_text(contents)
    data = load_electrochemistry_data(path, instrument=instrument)
    assert data.potential_reference_raw == "unknown"
    assert data.temperature_K is None
    with pytest.raises(ValueError, match="[Rr]eference"):
        data.get_potential_vs_rhe(ph=14)
    with pytest.raises(ValueError, match="[Tt]emperature"):
        data.get_potential_vs_rhe(reference="Ag/AgCl_3m", ph=14)
    before = data.potential_V.copy()
    corrected = data.get_potential_vs_rhe(reference="Ag/AgCl_3m", ph=14, temperature_k=298.15)
    np.testing.assert_allclose(corrected - before, 0.210 + 14 * 0.05915934968478, atol=1e-10)
    np.testing.assert_array_equal(data.potential_V, before)
    assert data.temperature_K is None


def test_rhe_identity_needs_no_ph_or_temperature():
    assert convert_potential_to_rhe(0.4, "RHE", ph=None) == pytest.approx(0.4)


@pytest.mark.parametrize("temperature", [None, 0, -1, np.nan, np.inf])
def test_reference_conversion_requires_finite_positive_temperature(temperature):
    with pytest.raises(ValueError, match="[Tt]emperature"):
        convert_potential_to_rhe(0.4, "Ag/AgCl_3m", 14, temperature_k=temperature)


@pytest.mark.parametrize("unit,values,expected", [
    ("percent", [1.0, 0.1, 0.01], [2, 3, 4]),
    ("fraction", [1.0, 0.1, 0.01], [0, 1, 2]),
])
def test_ftir_low_transmittance_known_values(unit, values, expected):
    data = CharacterizationData("FTIR", "synthetic", Path("unused"), "synthetic", x=np.arange(3), y=values,
                                y_unit=unit, metadata={"signal_type": "transmittance"})
    result = convert_transmittance_to_absorbance(data)
    np.testing.assert_allclose(result.y, expected, atol=1e-14)
    np.testing.assert_array_equal(data.y, values)
    assert result.processing_history[-1]["parameters"]["transmittance_unit"] == unit


def test_ftir_ambiguous_units_and_signal_are_rejected(tmp_path):
    path = tmp_path / "spectrum.csv"
    path.write_text("Wavenumber,Transmittance\n" + "\n".join(f"{400+i},0.1" for i in range(10)))
    with pytest.raises(InvalidUnitsError, match="transmittance_unit"):
        parse_ftir_file(path)
    result = parse_ftir_file(path, transmittance_unit="percent")
    np.testing.assert_allclose(convert_transmittance_to_absorbance(result).y, 3)
    path.write_text(path.read_text().replace("Transmittance", "Signal"))
    with pytest.raises(InvalidUnitsError, match="signal_type"):
        parse_ftir_file(path, transmittance_unit="percent")


def test_provenance_has_no_invented_qe_units(tmp_path):
    source = tmp_path / "sample.csv"
    source.write_bytes(b"potential,current\r\n1,2\r\n")
    output = record_figure_provenance("lsv", [source], [tmp_path], {"technique": "LSV"})
    record = json.loads(output.read_text())
    assert "units" not in record or record["units"] == {}
    assert "qe" not in record["parser_version"].lower()
    assert record["input_files"][0]["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    with pytest.raises(FileNotFoundError):
        record_figure_provenance("missing", [tmp_path / "missing.csv"], [tmp_path])


def test_export_preflights_all_collisions_and_preserves_existing(tmp_path):
    fig, _ = create_publication_figure()
    existing = tmp_path / "figure.svg"
    existing.write_bytes(b"protected-original")
    with pytest.raises(FileExistsError):
        save_publication_figure(fig, "figure", [tmp_path], formats=("png", "svg"), include_default_dirs=False)
    assert existing.read_bytes() == b"protected-original"
    assert not (tmp_path / "figure.png").exists()


def test_export_dimensions_deduplication_and_close(tmp_path, monkeypatch):
    fig, _ = create_publication_figure(columns=25.4, height_mm=25.4)
    calls = []
    original = fig.savefig
    def save(*args, **kwargs):
        calls.append(kwargs["format"])
        return original(*args, **kwargs)
    monkeypatch.setattr(fig, "savefig", save)
    paths = save_publication_figure(fig, "dimensions", [tmp_path, tmp_path / "."], formats=("png", "png"),
                                    dpi=100, include_default_dirs=False, close=True)
    assert calls == ["png"]
    with Image.open(paths["png"]) as img:
        assert img.size == (100, 100)
    assert not plt.fignum_exists(fig.number)


def test_export_render_failure_preserves_existing(tmp_path, monkeypatch):
    fig, _ = create_publication_figure()
    target = tmp_path / "figure.png"
    target.write_bytes(b"protected-original")
    def fail(*args, **kwargs):
        raise RuntimeError("render failed")
    monkeypatch.setattr(fig, "savefig", fail)
    with pytest.raises(RuntimeError, match="render failed"):
        save_publication_figure(fig, "figure", [tmp_path], formats=("png",), overwrite=True,
                                include_default_dirs=False, close=True)
    assert target.read_bytes() == b"protected-original"
    assert not plt.fignum_exists(fig.number)


def test_her_interpolation_known_value_and_unbracketed_target():
    from modules.electrochemistry.lsv import extract_overpotential_at_current
    assert extract_overpotential_at_current(np.array([-0.2, -0.1, 0]), np.array([-20, -5, 0]),
                                           reaction="HER") == pytest.approx(400 / 3)
    assert extract_overpotential_at_current(np.array([-0.2, -0.1]), np.array([-20, -15]),
                                           reaction="HER") is None


def test_conditions_header_units_and_explicit_override():
    from parsers.electrochemistry.schema import resolve_experimental_conditions
    header = {"Reference electrode": "Ag/AgCl_3m", "Temperature (C)": "25"}
    assert resolve_experimental_conditions(header) == ("Ag/AgCl_3m", 298.15)
    assert header["resolved_conditions"]["temperature_source"] == "Temperature (C)"
    assert resolve_experimental_conditions(header, "RHE", 310) == ("RHE", 310)
    assert header["resolved_conditions"]["reference_source"] == "caller"
    assert resolve_experimental_conditions({"Temperature": "25"}) == ("unknown", None)


@pytest.mark.parametrize("reference", ["Ag/AgCl", "Hg/HgO"])
def test_reference_conversion_requires_filling_solution(reference):
    with pytest.raises(ValueError, match="filling"):
        convert_potential_to_rhe(0.0, reference, 14, temperature_k=298.15)


def test_cp_does_not_silently_skip_requested_reference_conversion():
    from parsers.electrochemistry.schema import StandardizedElectrochemicalData
    from modules.electrochemistry.chronopotentiometry import plot_chronopotentiometry
    data = StandardizedElectrochemicalData(Path("synthetic"), "synthetic", "synthetic", "synthetic", "CP",
        potential_V=np.arange(10) / 100, time_s=np.arange(10) * 10,
        potential_reference_raw="Ag/AgCl_3m", temperature_K=298.15)
    with pytest.raises(ValueError, match="pH"):
        plot_chronopotentiometry(data, convert_to_rhe=True)


def test_provenance_collision_preflight(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text("actual source")
    first, second = tmp_path / "first", tmp_path / "second"
    second.mkdir()
    target = second / "figure_metadata.json"
    target.write_bytes(b"protected")
    with pytest.raises(FileExistsError):
        record_figure_provenance("figure", [source], [first, second])
    assert target.read_bytes() == b"protected"
    assert not first.exists()


def test_cli_rejects_uncalibrated_microscopy_before_export(tmp_path, monkeypatch):
    import sys
    from scripts.main_plots import main
    source = tmp_path / "uncalibrated.png"
    Image.fromarray(np.zeros((20, 20), dtype=np.uint8)).save(source)
    output = tmp_path / "derived"
    monkeypatch.setattr(sys, "argv", ["main_plots", "--type", "sem", "--input", str(source),
                                     "--output-root", str(output)])
    with pytest.raises(ValueError, match="calibration"):
        main()
    assert not list(output.rglob("*.png"))


def test_cli_ftir_explicit_units_exports_history(tmp_path, monkeypatch):
    import sys
    from scripts.main_plots import main
    source = tmp_path / "spectrum.csv"
    source.write_text("Wavenumber,Transmittance\n" + "\n".join(f"{400+i},1.0" for i in range(10)))
    output = tmp_path / "derived"
    monkeypatch.setattr(sys, "argv", ["main_plots", "--type", "ftir", "--input", str(source),
        "--transmittance-unit", "percent", "--formats", "png", "--output-root", str(output)])
    main()
    assert list(output.rglob("*.png"))
    metadata = json.loads(next(output.rglob("*_metadata.json")).read_text())
    dataset = metadata["analysis_metadata"]["processed_datasets"][0]
    assert dataset["processing_history"][-1]["parameters"]["transmittance_unit"] == "percent"


def test_cli_rejects_concentration_only_aas_as_calibration(tmp_path, monkeypatch):
    import sys
    from scripts.main_plots import main
    source = Path(__file__).parents[1] / "data/examples/characterization/aas/catalyst_digestion_loading.csv"
    monkeypatch.setattr(sys, "argv", ["main_plots", "--type", "aas", "--input", str(source),
                                     "--output-root", str(tmp_path)])
    with pytest.raises(ValueError, match="measured absorbance"):
        main()
    assert not list(tmp_path.rglob("*.png"))

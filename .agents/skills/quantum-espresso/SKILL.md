---
name: quantum-espresso
description: Guides Quantum ESPRESSO pw.x, DOS/PDOS, bands, NEB, and adsorption-energy work. Use when implementing, fixing, reviewing, or testing QE parsers, spin projections, Fermi alignment, or energy conversions.
---

# Quantum ESPRESSO

Use [scientific-validation](../scientific-validation/SKILL.md) for parser, unit, and quantitative changes. Read [QE file formats](references/qe-file-formats.md) for the affected executable and [energy conventions](references/energy-conventions.md) for thermodynamic work.

- Inspect the existing `parsers/quantum_espresso/` implementation and its tests before extending it. Reuse `parse_pw_output`, `parse_dos_file`, `QEPdosCollection`, `parse_bands_gnu`, `parse_neb_output`, and `adapters.py` where their supported format matches the input.
- Match the producing QE version and official file specification. Filename extensions and synthetic examples alone do not identify a valid format. Preserve source text and reject unsupported layouts explicitly.
- Keep `parse_pw_output(strict=True)` for final-result analysis. Check job completion, electronic convergence, and ionic convergence for relaxation separately. An energy line or `JOB DONE.` alone does not establish convergence. Verify final structures/cell units and ASE extraction status; missing structures are not valid geometries.
- Preserve native spin channels and orbital component order. Mirroring spin-down DOS is a plotting convention, not a negative physical density. Do not add LDOS to its own orbital components. Noncollinear/spin-orbit PDOS is explicitly unsupported.
- Require the provenance of the Fermi energy, energy zero, high-symmetry labels, reference calculations, and missing physical parameters. A fallback to unshifted energies must not be labeled `E - E_F`.
- NEB native columns are covered by pinned upstream QE 6.0 outputs. Read the reference for the tested scope; do not generalize this to every QE version. Reject `.path` restart files and require an explicit stdout log for convergence evidence.
- Classify fixtures using the [QE manifest](../../../tests/fixtures/qe/FIXTURE_MANIFEST.md). Its current examples are synthetic or defensive, not authentic QE outputs; even its format claims need comparison with official docs.

Run relevant cases in `tests/test_quantum_espresso.py`, `tests/test_dft.py`, and `tests/test_scientific_audit.py`. For newly supported layouts add independently sourced, redistributable real output when available, with version, settings, license, and hashes; report the genuine-data gap otherwise.

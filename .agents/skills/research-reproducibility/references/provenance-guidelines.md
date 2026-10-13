# Provenance and regeneration

## Minimum record for a research result

| Record | Contents |
| --- | --- |
| Inputs | Dataset class, source, sample/calculation ID, relative path or private identifier, byte SHA-256, source license/permission, actual instrument or simulation version |
| Conditions | Relevant geometry, method, calibration, electrolyte, temperature, pressure, reference electrode, wavelength/pixel scale, or other domain conditions |
| Processing | Ordered operations, parameters, units/reference conversions, masks/crops, fit windows, exclusions, and intermediate artifacts |
| Environment | Python and relevant installed package versions, platform when relevant, code revision plus dirty state, dependency/environment snapshot |
| Results | Quantities/units, uncertainty, model assumptions, warnings, validation evidence, and output paths/hashes |
| Regeneration | Exact root-relative command or notebook entry point, configuration, seeds and generator version, output destinations |

Store only relevant fields; unknown values stay unknown. Hashes identify bytes, not authenticity, scientific validity, or permission to share.

## Existing facilities and limits

- `parsers/characterization/common.py` provides byte hashing, `CharacterizationData` history and `to_provenance_dict`, and `record_characterization_provenance(figure_name, input_files, output_dirs, technique, metadata=None)`.
- `parsers/quantum_espresso/adapters.py` provides `record_figure_provenance`; `parsers/electrochemistry/schema.py` also provides byte hashing. Reuse the appropriate public helper rather than multiplying implementations. Both figure recorders now reject missing inputs and protect existing sidecars by default (`overwrite=False`).
- NEB now hashes exact source bytes, consistent with the provenance helpers. Verify hash semantics in other parsers when extending them.
- Figure recorders collect inputs and caller-supplied metadata; they do not automatically capture the entire environment or processing history. Pass the needed configuration/history explicitly. Verify every input was found and every intended sidecar exists.
- `record_figure_provenance` no longer invents QE units or Fermi alignment. Record verified domain units/reference in `analysis_metadata`; absent values remain absent. Older sidecars may retain misleading defaults and must not be retroactively treated as corrected.
- `save_publication_figure` writes images only. Calling it alone does not create a provenance sidecar. Sidecars may contain absolute paths; redact private identifiers in a separate shareable copy without altering the protected original record.

## Repeatability decisions

For fixture work, inspect `scripts/generate_synthetic_data.py`, `generate_electrochem_fixtures.py`, `generate_characterization_fixtures.py`, `generate_aas_fixtures.py`, or `tests/fixtures/qe/create_fixtures.py` before executing. Preserve their declared synthetic identity. Some use global NumPy random state, so call order and package versions are part of reproducibility. Do not rerun generators just to get a clean test report.

For figures, the CLI `python -m scripts.main_plots --category <domain>` selects demonstration workflows; inspect `--help` and the implementation for direct data options. Save commands/configurations with the result. Existence of a CLI flag does not prove every intended configuration field is implemented.

For sharing, review the exact staged paths and notebook outputs. `.gitignore` currently does not exclude all research data or generated results, and historical test runs wrote tracked figure artifacts; current export tests use isolated temporary destinations. Stage deliberately; never use a broad cleanup/reset that could erase user work. Keep unrelated historical records intact and clearly date new validation evidence.

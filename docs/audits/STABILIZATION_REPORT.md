# Scientific repository stabilization

Date: 2026-10-13. Scope: confirmed correctness, export, dependency, warning and skill-guidance defects; no new scientific techniques or frameworks. This supersedes the current-state limitations in the earlier [skills audit](../AGENT_SKILLS_VALIDATION.md), which remains a historical record.

## Baseline and preservation

- Interpreter: existing `research-graphing/Scripts/python.exe`, Python 3.11.0 on Windows. No packages were installed, removed, upgraded or downgraded; no environment was recreated.
- Initial Git status: modified `README.md`; untracked `.agents/`, `AGENTS.md`, `docs/AGENT_SKILLS_VALIDATION.md` from the preceding task. Those additions were retained and relevant guidance updated. Scientific implementation files were initially clean. No commit or push was made.
- Baseline: `python -m pytest -q` returned **132 passed, 34 warnings**, 8.08 seconds. The [complete warning summary](BASELINE_WARNINGS.txt) is retained with machine-specific root paths removed.
- Before testing, SHA-256 values were recorded for **730 original files** under `data/`, original fixture directories, `notebooks/`, `figures/` and `outputs/`. Final comparison: **730 identical, zero changed/missing**. New upstream QE fixtures are additions, not replacements.
- The baseline export test rewrote two clean tracked artifacts, `outputs/vector/test_fig.pdf` and `.svg`. Only those test-generated changes were restored from Git. The corrected test now writes exclusively under pytest temporary directories.
- All **157 installed distribution names/versions** match the initial snapshot. [TESTED_ENVIRONMENT.json](TESTED_ENVIRONMENT.json) records them. This does not establish that binaries or external system libraries were independently audited.

## Confirmed scientific defects and corrections

### Quantum ESPRESSO NEB

The old parser interpreted native three-column `.dat` files as image index / coordinate / energy, shifting both physical columns. It also treated an invented synthetic `.path` summary as if it represented official restart output, and selected early log matches rather than final-iteration values.

The [official output-file guide](https://www.quantum-espresso.org/Doc/neb_user_guide/node7.html) defines `.dat` as coordinate / relative energy / residual error and describes `.path` as an undocumented restart file. The corrected parser uses the first two columns, preserves arbitrary coordinate units, rejects restart files, and accepts an explicitly supplied NEB stdout log through `log_filepath`. Two-column profiles remain supported. Legacy indexed tables require an explicit header or `data_format="legacy_indexed"`; the positional `path_filepath` argument remains a log alias, but `.path` restart input is rejected.

Convergence requires a positive marker followed by `JOB DONE.` in the final logged run, without a later failure marker. The last reported forward/reverse barriers and climbing image are used. With no log, discrete-image barriers are identified as estimates; `is_converged=False` means convergence was not established, not proof of failure, and no climbing image is invented. Invalid/nonfinite/ragged tables, repeated coordinates and unsupported shapes raise errors. Hashes use original bytes.

Authentic evidence: unchanged H2+H reference files from the [official QEF/q-e repository at pinned commit 61569eb](https://github.com/QEF/q-e/tree/61569eb480231b649c47f631df9c5c83537df461/NEB/examples/example01/reference). Stdout identifies QE 6.0 and 8 February 2017. The [fixture manifest](../../tests/fixtures/qe/upstream_neb/README.md), exact source URLs/hashes and upstream GPLv2 license accompany the files. No local QE calculation was run.

| Independent expected quantity | Regression result |
| --- | --- |
| Second coordinate | 0.1939576458, tolerance 1e-10 |
| Fourth-image relative energy | 0.2042674385 eV, tolerance 1e-10 |
| Final stdout forward / reverse barriers | 0.204267 / 0.204268 eV, tolerance 1e-7 |
| Climbing image / convergence | Image 4 / positive final evidence |
| Actual restart, negative/truncated log | Restart rejected; no convergence invented |

The figure now labels electronic energy and arbitrary reaction coordinate, rather than claiming free energy or angstroms. Barriers are not Gibbs activation energies without additional thermodynamic corrections. Broader QE release/system compatibility remains unvalidated; callers must associate matching profile and log files. The residual-error unit is not inferred: guide/version conventions require further verification before exposing it quantitatively.

### Electrochemical conditions and derived quantities

BioLogic's implicit Ag/AgCl and the schema's 298.15 K fallback could alter RHE/OER results without experimental evidence. Both instrument families now retain missing reference/temperature as `unknown`/`None`. Explicit caller values override clear headers, temperature headers require K/C units, and `resolved_conditions` records origins. Parse errors propagate instead of triggering an unrelated instrument fallback; `.mps` settings are not measurement data.

RHE conversion requires confirmed filling-solution variant, finite pH, and finite positive temperature when a temperature-dependent correction is performed. Unqualified Ag/AgCl and Hg/HgO are rejected in conversions. RHE-to-RHE identity needs no pH/temperature. The explicitly tabulated standard-reference lookup retains its documented 298.15 K default and historical aliases; these are lookup conventions, not measured cell conditions. OER overpotential/benchmarks use the supplied temperature; HER's RHE equilibrium remains zero.

Physical rationale: `E_RHE = E_measured + E_reference,SHE(T) + (RT/F) ln(10) pH`. A filling-solution choice changes the reference offset; temperature changes the Nernst correction. At explicitly supplied 298.15 K and pH 14, the existing 3 M KCl table value gives an independently checked shift `0.210 + 14 * 0.05915934968478` V (absolute tolerance 1e-10). Tests preserve the raw potential array and confirm absent metadata stays absent after a getter conversion. Tabulated temperature coefficients and ideal-activity assumptions were not recertified against laboratory calibration.

Additional confirmed downstream defects: HER interpolation supplied descending current values to `np.interp`; corrected ascending/bracketed interpolation returns 133.333333 mV for the independent points (-0.2 V, -20 mA/cm2), (-0.1 V, -5 mA/cm2), (0 V, 0 mA/cm2) at -10 mA/cm2. Raw-reference plots no longer report RHE overpotentials. Benchmarks use mA/cm2 independently of display units. Equilibrium lines use the requested scale/temperature. LSV/Tafel iR paths use the schema guard against missing resistance/double compensation. CP no longer silently skips a requested RHE conversion when pH is missing.

Compatibility: existing positional arguments remain; callers that relied on unstated conditions must provide them. Demonstration fixtures use explicitly recorded synthetic 298.15 K; the ambiguous HER demo also records a saturated-KCl assumption. This assumption is not asserted as a measurement. Original notebooks/data were not rewritten; older notebook calls need explicit conditions when rerun.

### FTIR units

The prior maximum-value heuristic interpreted very low percent transmittance as a fraction, changing absorbance by two units. Parsing and conversion now require `percent` or `fraction` from explicit input or an unambiguous column header. Conflicting units and ambiguous signal types raise errors. CSV delimiters preserve headers with spaces. The existing positional `clamp_nonpositive` parameter remains; optional clamping is explicit and recorded, not silently applied. The unsupported `detected_bands` CLI plotting argument was also corrected.

Independent Beer-Lambert checks: 1%, 0.1%, 0.01% give absorbances 2, 3, 4; fractions 1, 0.1, 0.01 give 0, 1, 2, using `A = -log10(T_fraction)` with 1e-14 tolerance. Tests verify unchanged input values and recorded conversion units. Reflectance is not silently converted into absorbance. Sampling mode is no longer fabricated.

### Truthful provenance and missing characterization conditions

The generic recorder emitted QE units and `fermi_aligned=True` even for electrochemistry. Those defaults are removed. Both recorders require real input files, hash exact bytes, deduplicate destinations and reject existing sidecars by default. Callers supply actual analysis information; a missing field is not an inferred fact. Direct CLI NEB includes all consumed profile/log sources; comparisons include both sources; FTIR includes actual processed dataset histories and units; electrochemistry includes resolved condition origins.

Inspection also found direct CLI microscopy scales and ICP/AAS loading/leaching sample conditions hardcoded without measurement evidence. Microscopy now requires calibrated pixels, and techniques are labeled correctly. ICP/AAS convenience routes lacking a way to specify necessary volume/mass/dilution/metal fraction now reject with guidance to use the existing Python functions and explicit experimental conditions. Concentration-only AAS input is also rejected as calibration data: after parsing discarded loading fields, it previously produced an empty plot from missing absorbance values. A CLI regression covers this failure. No new analysis API was introduced. Historical sidecars/figures remain unchanged and must not be treated as retroactively corrected. Recorders still rely on callers for full processing, environment and uncertainty metadata.

## Export reliability and warnings

`save_publication_figure` keeps existing positional arguments and additive default routing, with these deliberate safety changes:

- `overwrite=False`; every target is checked before rendering/writing. Exclusive creation protects against a concurrent writer. Opt-in replacement is atomic per file after all formats render successfully.
- Paths/formats are deduplicated. `include_default_dirs=False` limits exports to explicit destinations; default behavior still includes `outputs/vector` or `outputs/raster`. Returned paths retain the last-destination convention.
- `bbox_inches=None` preserves requested canvas size even when style configuration requests tight cropping. Explicit `bbox_inches="tight"` opts into changed bounds.
- `close=True` closes the supplied figure on success or failure; default library ownership stays with the caller. The CLI closes its figures and supports `--output-root` and `--overwrite`.
- Figure and sidecar sets are not one filesystem transaction. Mid-write I/O failure can leave a partial new set; no guarantee of all-file rollback is made. Default collisions protect existing files.

Regression tests cover preflight collisions across formats/destinations, existing-byte preservation, render failures, duplicate destinations/formats, cleanup, and a 25.4 mm canvas at 100 DPI producing exactly 100 x 100 pixels. Sidecar collisions are tested separately.

| Baseline warning category | Count | Root cause and action |
| --- | ---: | --- |
| Figure layout | 29 | `tight_layout` replaced active constrained layouts, sometimes with incompatible axes. Shared finalization now preserves the active layout engine. |
| Deprecation | 4 | Eagerly evaluated nested `getattr` touched deprecated scikit-image axis properties. Use supported properties directly. |
| Resource management | 1 | Tests retained over 20 figures. Autouse teardown closes only each test's new figures; CLI cleanup covers success/failure. |
| Scientific correctness risk | 0 emitted | Missing-condition/format defects were silent and required separate investigation. |
| Other third-party library warnings | 0 | The layout/deprecation messages above originate in libraries but are triggered by repository usage. |

Removed the pre-existing global `DeprecationWarning` filter. Source compilation additionally found invalid escapes in AAS/ECSA labels and the notebook generator; escaping was repaired without running the generator or changing notebooks. Compiling all 85 Python source/test/style files with warnings enabled produces no warnings. The final full pytest run has no warnings. The electrochemistry demonstration retains its intentional `CdlInterpretationWarning`: Cdl/ECSA alone does not measure active-site density. This scientific warning is not suppressed.

## Dependencies

One tested configuration is Python **3.11.0**, NumPy **2.4.6**, pandas **3.0.6**, SciPy **1.17.1**, Matplotlib **3.11.2**, scikit-image **0.26.0**, pytest **9.1.1**. Exact direct dependencies are in [requirements-tested.txt](../../requirements-tested.txt); the full observed environment is linked above.

`requirements.txt` now permits the installed NumPy 2.x line and declares PyYAML, which was already installed and is needed by the skill validator. `environment.yml` specifies Python 3.11, Conda name `catalysis-workbench`, and the same pip snapshot. The existing venv directory remains `research-graphing`; the names refer to different environment mechanisms. README activation now matches the Conda name and no longer claims Python 3.12 was tested.

`python -m pip check`: **No broken requirements found**. All installed direct versions satisfy both manifests. A clean pip or Conda solve/install was **not performed**. The snapshot is not a transitive lock, package-index availability guarantee, or proof of cross-platform/Python 3.12 compatibility. No installation commands were executed.

## Verification

Commands use `research-graphing/Scripts/python.exe` from the repository root; replace `python` below with that executable on this checkout.

| Check | Result |
| --- | --- |
| `python -m pytest -q` | **164 passed, zero warnings** (7.70 s); original 132 retained, 32 added parameterized regressions |
| `python .agents/skills/scientific-validation/scripts/validate_skills.py` | Seven skills; frontmatter, structure and local references validated |
| `python -m pytest -q .agents/skills/scientific-validation/scripts/test_validate_skills.py` | Five validator tests and seven subtests passed |
| System skill-creator `quick_validate.py`, each existing skill | Seven validated |
| Source compile with warnings enabled | 85 files; zero warnings |
| CLI `--help`, direct exports and negative cases | Verified as described below |
| `git diff --check` | Clean after restoring the two baseline test artifacts |
| Protected original files / installed package versions | 730 unchanged / 157 unchanged |

Existing scientific assertions were retained; old synthetic calls now declare their required units/conditions. Authentic NEB and independent FTIR/reference/HER values supplement synthetic tests. Tests were not relaxed to accept wrong numerical behavior.

### CLI and actual figure inspection

All outputs were generated under a new temporary audit root, outside the repository. The direct checks used all four formats (`--formats pdf svg png tiff`):

```powershell
python -m scripts.main_plots --type neb --input tests/fixtures/qe/upstream_neb/H2+H.dat --neb-log tests/fixtures/qe/upstream_neb/H2+H.out --output-root "$env:TEMP/catalysis-check-neb"
python -m scripts.main_plots --type ftir --input data/examples/characterization/ftir/catalyst_functionalized_ftir.csv --transmittance-unit percent --output-root "$env:TEMP/catalysis-check-ftir"
python -m scripts.main_plots --type aas --input data/examples/characterization/aas/faas_calibration_fe.csv --output-root "$env:TEMP/catalysis-check-aas"
python -m scripts.main_plots --type lsv --input data/examples/electrochemistry/biologic/biologic_her_lsv.mpt --reaction HER --reference-electrode Ag/AgCl_sat --temperature-k 298.15 --output-root "$env:TEMP/catalysis-check-lsv"
python -m scripts.main_plots --type sem --input data/examples/characterization/microscopy/sem_nanoparticles.png --pixel-size 1.0 --output-root "$env:TEMP/catalysis-check-sem"
```

Choose fresh output directories; repeat export requires explicit `--overwrite`. The latter two commands use stated assumptions for synthetic examples only. The complete electrochemistry and characterization demonstrations also passed with `--category <name> --formats png --output-root <isolated-root>`, with only the intentional Cdl interpretation warning described above. No generator was rerun to obtain passing tests.

Negative smoke cases reject missing FTIR units, missing microscopy calibration, unspecified AAS digestion conditions and existing output collisions; existing NEB export hashes remain unchanged. Bare Ag/AgCl was separately observed to fail rather than choose a filling solution. Non-QE sidecars contain no fabricated QE/Fermi fields.

PNG renderings for all five direct workflows were visually inspected. NEB annotation/legend collisions and the AAS annotation overlapping its title were corrected. Reviewed axis labels, scales, legends, residual panels and SEM scale bar. No raster clipping was observed in these examples; automated segmentation remains a candidate analysis needing scientific review for merged particles.

[EXPORT_CHECKS.json](EXPORT_CHECKS.json) records PDF MediaBox, SVG dimensions/text, raster pixel sizes/DPI and PDF TrueType font-stream checks. All canvases are 85 mm wide: NEB/FTIR/LSV 1003 x 620 pixels, AAS 1003 x 743, SEM 1003 x 1003 at 300 DPI (PNG metadata rounds to 299.9994). TIFF dimensions match PNG; vector widths are 240.944882 pt. SVG text is preserved but fonts are not embedded. Font presence checks and PNG inspection do not prove rendering in every journal/Illustrator/Inkscape environment.

## Agent Skills discovery and current guidance

Exactly seven skills remain in `.agents/skills`; no duplicate `.codex/skills` or `.agent/skills` collection was introduced. The root `AGENTS.md`, README and relevant skill references now match corrected signatures and behavior. The [Codex documentation](https://developers.openai.com/codex/skills) and [Antigravity documentation](https://antigravity.google/docs/skills) identify repository `.agents/skills` discovery; their locations were checked against the repository.

Observed runtime evidence: all seven names/descriptions/paths appeared in the current Codex session's available-skills catalog and their instructions were read and applied. This is continuing-session catalog visibility and explicit application, **not an independent fresh-session automatic-selection test**. No Antigravity runtime was available. Fresh-session activation on either platform remains **unverified**.

Manual verification, using a separate fresh conversation for each row:

1. Open this repository as the workspace in a new Codex/Antigravity session after the files are present. Confirm the workspace root. Inspect the skill list/picker (Antigravity: Customizations) for these seven names and their `.agents/skills` paths.
2. Send the exact prompt below without naming skills, then inspect tool/file-read activity for the expected `SKILL.md` files. Ask the agent to stop after a read-only plan; no data edits are needed for this test.
3. Record client/version, time, prompt, observed catalog entries and actual reads. A skill merely naming itself is insufficient evidence. If selection fails, test an explicit `$skill-name` request in Codex or explicit skill name/path in Antigravity and record that separately as explicit invocation, not automatic routing.

| Exact prompt (append "Read only; stop after a plan.") | Expected minimum reads |
| --- | --- |
| Fix a Quantum ESPRESSO PDOS parser. | quantum-espresso; scientific-validation |
| Generate an HER LSV figure from BioLogic data. | electrochemistry; publication-figures |
| Analyze Raman spectra and fit peaks. | materials-characterization; scientific-validation |
| Audit adsorption energy unit conversions. | quantum-espresso; scientific-validation |
| Create journal-ready figures from SEM images. | materials-characterization; publication-figures |
| Ensure original research datasets are preserved. | research-reproducibility |
| Maintain the Python APIs and dependency declarations in this scientific toolkit. | scientific-python |

Additional relevant skills may load; exact exclusive selection is not required. Derived quantitative figures also warrant scientific-validation. Structural validators cannot prove activation or scientific validity.

## Changed files and outstanding limits

| Files | Changes |
| --- | --- |
| `parsers/quantum_espresso/neb.py`, `adapters.py` | Native NEB interpretation, explicit logs/restart rejection, byte hashes, truthful provenance, plot labels/annotations |
| `parsers/electrochemistry/{__init__,biologic,metrohm,schema}.py` | Missing conditions, origin tracking, dispatcher error propagation, guarded corrections |
| `modules/electrochemistry/{reference_electrodes,lsv,tafel,benchmarking,chronopotentiometry}.py` | Explicit corrections, HER interpolation, scale/unit/temperature consistency and downstream propagation |
| `parsers/characterization/{ftir,common}.py` | Explicit units/history and safe truthful sidecars |
| `scripts/utils.py`, `scripts/main_plots.py` | Safe exports/layout, isolated CLI destinations, cleanup, explicit inputs/provenance and fail-fast unsupported convenience routes |
| `modules/characterization/{aas,bet,comparative,ftir,icp,microscopy,raman,xrd}.py`; `modules/electrochemistry/{chronoamperometry,eis}.py` | Layout lifecycle, microscopy property deprecations, annotation placement |
| `modules/electrochemistry/ecsa.py`, `scripts/build_notebooks.py` | Invalid escape repairs only; notebooks not generated |
| `tests/test_stabilization.py`, existing FTIR/electrochemistry/QE/export tests, `tests/conftest.py`, `pytest.ini` | Numerical/protection regressions, explicit fixture conditions, cleanup and warning visibility |
| `tests/fixtures/qe/upstream_neb/` | New authentic reference files, hashes, attribution/license; original fixtures unchanged |
| `requirements.txt`, `requirements-tested.txt`, `environment.yml`, `README.md` | Consistent tested Python/NumPy declarations and migration guidance |
| `AGENTS.md`, existing skill instructions/references, `docs/AGENT_SKILLS_VALIDATION.md`, `docs/audits/` | Current workflows, historical follow-up and audit evidence |

Remaining scientific limits: authentic instrument/binary imports, Raman/peak model selection, BET/BJH applicability, quantitative microscopy calibration, Tafel kinetic-window selection and laboratory reference-electrode calibration are not validated by generated fixtures. Nonmonotonic LSV branches need scientific selection before interpolation. The full repository is not certified publication-ready. Historical notebooks were preserved and not executed. A clean dependency installation, broader QE compatibility and genuine fresh-session skill activation remain unverified. These limits are retained explicitly rather than filled with invented physical metadata.

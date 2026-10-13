# Agent Skills implementation and validation

Audit date: 2026-10-13. This record concerns repository guidance and support tooling, not a new validation of all scientific methods.

## Scope and preserved files

The initial Git working tree was clean. No repository `AGENTS.md`, skill collection, or Codex/Antigravity rules were found outside the ignored Python environment. Existing `CONTEXT.md`, historical scientific/handoff reports, notebooks, parsers, analysis modules, styles, requirements, and research files were retained. README changes add a new section only. No scientific equations or parsers were changed and no dependencies were installed.

Added:

- [Root guidance](../AGENTS.md).
- Seven `SKILL.md` files in [the shared collection](../.agents/skills/): scientific-python, quantum-espresso, electrochemistry, materials-characterization, scientific-validation, publication-figures, research-reproducibility.
- Nine focused references: architecture; QE file formats and energy conventions; instrument imports and electrochemical conventions; characterization methods; numerical validation checklist; figure standards; provenance guidelines.
- A [read-only validator](../.agents/skills/scientific-validation/scripts/validate_skills.py) and [isolated validator tests](../.agents/skills/scientific-validation/scripts/test_validate_skills.py).
- This audit/maintenance record.

## Repository evidence

Inspected README/CONTEXT and existing reports, the Python dependency manifests and active installed packages, `parsers/`, `modules/`, shared script/style utilities, CLI options, fixture generators/manifests, and the pytest suite. Verified callable names and important behaviors against source rather than copying historical architecture claims.

Key boundaries reflected in the skills:

- Characterization has `CharacterizationData` in `parsers/characterization/common.py`, not the specialized schema classes listed in the handoff. Isotherm import is in `adsorption_isotherm.py`.
- QE fixtures are synthetic/defensive. The official NEB guide's three-column output differs from `neb.py`'s current interpretation; official `.path` files are undocumented restart files, not the repository's synthetic summary fixture. See [QE format limitations](../.agents/skills/quantum-espresso/references/qe-file-formats.md).
- BioLogic text import has an Ag/AgCl default and conversion has a temperature fallback; these are not measured conditions. `.mps` parses settings, not measurement data.
- The BET routine evaluates a supplied window; it is not an automatic optimal-window search. Existing peak/FFT functions do not imply Rietveld refinement or automatic lattice indexing.
- `save_publication_figure` adds destinations while retaining default output writes, uses tight bounding boxes, defaults to overwriting, and does not itself generate provenance. SVG text preservation is not font embedding.
- The generic figure provenance recorder emits QE unit/Fermi defaults even for non-QE callers; FTIR percent/fraction conversion uses a magnitude heuristic. Guidance identifies these limits instead of claiming metadata/units are automatically validated.
- Historical "production ready" and universal publication-validation claims were not adopted as new scientific evidence.

## Environment used

Used the existing Windows `research-graphing/Scripts/python.exe`, Python 3.11.0. Relevant installed versions:

| Package | Version |
| --- | --- |
| NumPy / pandas / SciPy | 2.4.6 / 3.0.6 / 1.17.1 |
| Matplotlib / pytest / PyYAML | 3.11.2 / 9.1.1 / 6.0.3 |
| lmfit / statsmodels / pybaselines | 1.3.4 / 0.15.0 / 1.2.1 |
| impedance / galvani / pygaps | 1.7.1 / 0.5.0 / 4.6.1 |
| scikit-image / Pillow / tifffile | 0.26.0 / 12.3.0 / 2026.3.3 |
| Pint / ASE / pymatgen | 0.25.3 / 3.29.0 / 2026.9.24 |
| openpyxl / SciencePlots / seaborn | 3.1.5 / 2.2.2 / 0.13.2 |

NumPy is outside the declared `>=1.26.0,<2.0.0` range. The Conda manifest specifies Python 3.12 and environment name `catalysis-workbench`, while the existing venv uses Python 3.11.0. The default system Python is a different environment; it was not used for regression testing. Passing in this venv does not establish a clean-install result from the manifests. Resolve this separately without silently rebuilding the user's environment.

## Verification

Commands use the project interpreter from the repository root:

```text
python .agents/skills/scientific-validation/scripts/validate_skills.py
python .agents/skills/scientific-validation/scripts/test_validate_skills.py
python -m pytest -q
python -m scripts.main_plots --help
git diff --check
```

Results:

- All seven skills passed the bundled Codex skill-creator `quick_validate.py` check. The available validator was used without installing `skills-ref`.
- The repository validator passed YAML/type/name checks and inline local file links across all skill references and guidance documents, including README. Five isolated validator tests passed, covering valid/invalid YAML, mismatched names, duplicate/unsupported fields, empty/missing manifests, broken/escaping links, and an absent skill collection.
- Manual review checked task descriptions, the six routing examples below, scientific safeguards, and API names against source. A source-definition scan checked 61 scientific API identifiers mentioned in inline code; the two remaining matches were script filenames, not claimed APIs. No full behavioral skill execution was performed.
- The complete existing Python suite passed **132 tests with 34 warnings** both before and after the additions. Final run: `python -m pytest -q --disable-warnings`, 6.92 seconds. Warnings include layout/tight-layout conflicts, more than 20 open figures, and deprecated scikit-image property aliases; these existing issues were not suppressed in code or fixed in this task.
- `python -m scripts.main_plots --help` completed and confirmed the documented category commands. No figures or scientific numerical behavior were intentionally modified, so demonstration batches and new visual inspection were not required.
- The existing tests wrote tracked `outputs/vector/test_fig.pdf` and `.svg`; only artifacts produced by this run from an initially clean state were restored. Final scope review and `git diff --check` exclude those generated changes.

External link verification checked 21 unique URLs on the audit date. Nineteen returned HTTP 200, including the specification, both agent platforms, all QE documentation, instrument documentation, library documentation, and EPA methods. The CHE DOI resolved to the ACS publisher but returned HTTP 403; the IUPAC DOI resolved to the publisher with HTTP 202 and its full HTML was not verified. The IUPAC report's identity was corroborated by the publisher's indexed PDF. These publisher restrictions are verification limits, not a claim of broken DOIs or verified full-text access. Recheck external sources for the producing software/method version during scientific work.

## Discovery review

The [open specification](https://agentskills.io/specification), [Codex local-skill documentation](https://developers.openai.com/codex/skills), and [Antigravity skills documentation](https://antigravity.google/docs/skills) were consulted live. Both platforms document `.agents/skills` as a repository/workspace discovery location. Frontmatter uses only `name` and `description`; no duplicate vendor collection or vendor-specific metadata was added.

The following is a manual review of descriptions and routing, not an observed automatic-activation test:

| Representative task | Expected minimum skills | Description/routing evidence |
| --- | --- | --- |
| Fix a Quantum ESPRESSO PDOS parser | quantum-espresso; scientific-validation | PDOS/parsers and parser changes are explicit triggers |
| Generate an HER LSV figure from BioLogic data | electrochemistry; publication-figures | BioLogic/HER/LSV and figure work are explicit triggers |
| Analyze Raman spectra and fit peaks | materials-characterization; scientific-validation | Raman and peak fitting/quantitative analysis are explicit triggers |
| Audit adsorption energy unit conversions | quantum-espresso; scientific-validation | Adsorption/energy conversions and unit audits are explicit triggers |
| Create journal-ready figures from SEM images | materials-characterization; publication-figures | SEM images and microscopy panels/figures are explicit triggers |
| Ensure original research datasets are preserved | research-reproducibility | Original datasets/raw-data preservation are explicit triggers |

Scientific-validation also applies when a figure task derives new quantities; scientific-python may apply when implementation is needed. No exact exclusive skill set is promised.

Newly created skills were not available through a fresh platform discovery/activation session during this task. To test real discovery, open the repository in a new Codex/Antigravity conversation, inspect the skill list/picker (Antigravity IDE: Customizations), issue each prompt, and record which `SKILL.md` files are actually loaded. Explicitly request a skill by name or file path if discovery is unavailable. Do not count this manual table or the structural validator as automatic activation evidence.

## Follow-up opportunities

Reconcile dependency constraints and the tested environment; validate authentic QE/instrument files with permission and traceable provenance, especially NEB; address existing plotting/resource/deprecation warnings separately; and run the fresh-session activation checks above. These are recommendations, not additions to this guidance-only task.


## Stabilization follow-up (2026-10-13)

This file records the earlier guidance-only audit. Its defect, dependency and warning statements describe that baseline, not the corrected implementation. The [stabilization report](audits/STABILIZATION_REPORT.md) records subsequent fixes, authentic QE reference coverage, tests, data preservation, and current skill verification. All seven skills were visible in the Codex tool catalog during the continuing session and explicitly read; independent fresh-session automatic selection and Antigravity runtime activation remain unverified.

# Working in catalysis-workbench

This repository is a scientific Python post-processing and visualization toolkit for Quantum ESPRESSO, catalysis, electrochemistry, materials characterization, and chemical engineering.

## Architecture and environment

- Read existing functions and their callers before adding replacements. File ingestion and format validation belong in `parsers/`; scientific analysis and domain plots belong in `modules/`. Shared paths, conversions, plotting, and CLI orchestration live in `scripts/`; styles and palettes live in `styles/`. Notebooks demonstrate these APIs.
- Characterization uses `CharacterizationData` in `parsers/characterization/common.py`; electrochemistry uses `StandardizedElectrochemicalData` in `parsers/electrochemistry/schema.py`. Do not infer APIs from historical reports.
- Reuse the existing environment. On this Windows checkout use `.\research-graphing\Scripts\python.exe`; otherwise use the project's activated interpreter. Run from the repository root. `environment.yml` names the Conda environment `catalysis-workbench`.
- Dependencies are declared in `requirements.txt`, the tested direct-dependency snapshot `requirements-tested.txt`, and `environment.yml` (Python 3.11). Check installed versions before using unfamiliar APIs; report drift rather than silently rebuilding the environment. Match existing type hints, NumPy-style docstrings, pathlib paths, and pytest conventions.

## Verification and demonstrations

```powershell
.\research-graphing\Scripts\python.exe -m pytest -q
.\research-graphing\Scripts\python.exe -m scripts.main_plots --help
.\research-graphing\Scripts\python.exe -m scripts.main_plots --category electrochemistry
.\research-graphing\Scripts\python.exe -m scripts.main_plots --category characterization
```

With an activated environment, replace the executable with `python`. Demo commands write figures and provenance to `figures/` and `outputs/` by default. Pass `--output-root <new-directory>` for isolated exports. Collisions fail unless `--overwrite` is explicit. Generators may overwrite synthetic examples and fixtures; inspect destinations before running them. Tests now isolate exports and close their figures. Supply confirmed electrochemical reference/pH/temperature and FTIR transmittance units; see `docs/audits/STABILIZATION_REPORT.md` for migration and validation limits.

## Scientific and data requirements

- For parser or quantitative changes, apply scientific-validation alongside the domain skill. Require meaningful tests before accepting scientifically significant or safety-critical changes. Check dimensions, units, signs, reference states, boundary conditions, and independent numerical examples.
- Treat historical reports and `CONTEXT.md` as background, not proof of scientific validity. Report execution, test results, trusted-reference agreement, assumption suitability, and genuine-data validation separately. Never equate synthetic fixtures with authentic measurements or simulations.
- Preserve original experimental data and simulation outputs, including restart files. Work on derived copies with hashes, parameters, processing history, and software versions. Never silently rewrite measurements or commit confidential/unpublished datasets without explicit permission.
- Reuse `styles/publication.mplstyle`, `styles/palettes.py`, and `scripts.utils.create_publication_figure` / `save_publication_figure`. Inspect rendered figures when graphical output changes. Journal requirements override the generic presets.

## Skill routing

Discover the focused [skills](.agents/skills/) by task and read the relevant `SKILL.md` and only needed references:

| Task | Skill |
| --- | --- |
| Python implementation, APIs, dependencies | [scientific-python](.agents/skills/scientific-python/SKILL.md) |
| QE outputs, DOS/PDOS, bands, NEB, adsorption energies | [quantum-espresso](.agents/skills/quantum-espresso/SKILL.md) |
| Autolab/BioLogic, CV/LSV/EIS, HER/OER, CHE | [electrochemistry](.agents/skills/electrochemistry/SKILL.md) |
| XRD, FTIR, Raman, microscopy, ICP/AAS, BET/BJH | [materials-characterization](.agents/skills/materials-characterization/SKILL.md) |
| Parsers, equations, quantitative analysis, scientific review | [scientific-validation](.agents/skills/scientific-validation/SKILL.md) |
| Figure creation, edits, exports | [publication-figures](.agents/skills/publication-figures/SKILL.md) |
| Datasets, provenance, regeneration, research privacy | [research-reproducibility](.agents/skills/research-reproducibility/SKILL.md) |

Final implementation reports should state changed behavior/files, commands and results, numerical evidence, assumptions and uncertainties, data provenance, figure inspection where applicable, and remaining limitations. Report observed evidence rather than blanket publication-readiness claims.

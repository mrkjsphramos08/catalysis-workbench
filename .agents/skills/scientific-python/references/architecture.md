# Actual architecture and reuse points

Paths below are relative to the repository root. Recheck signatures when implementing; this is a navigation map, not a second API specification.

| Concern | Existing implementation |
| --- | --- |
| QE text output | `parsers/quantum_espresso/{pw,dos,pdos,bands,neb}.py`; public exports in `__init__.py`, figure/adsorption adapters in `adapters.py` |
| DFT | `modules/dft/electronic_structure.py`: alignment, DOS/PDOS/band plots; `adsorption.py`: adsorption energies/scaling |
| Catalysis | `modules/catalysis/reaction_energy.py`: Gibbs energies, CHE, SHE/RHE conversion; `volcano.py`: generic, ORR, HER models |
| Chemical engineering | `modules/chemical_engineering/kinetics.py`: `fit_arrhenius`, conversion/selectivity, reactor-profile, parity/residual plots; no general reactor solver should be inferred |
| Electrochemistry | `parsers/electrochemistry/`: dispatcher, column mapping, instrument parsers, `schema.py`; `modules/electrochemistry/`: technique analysis and reference conversions |
| Characterization | `parsers/characterization/common.py`: `CharacterizationData`, errors, hashes/history; technique parsers, including `adsorption_isotherm.py`; there is no characterization `schema.py` or parser `bet.py` |
| Shared characterization numerics | `modules/characterization/common_peaks.py`: baseline, smoothing, peaks, band integration; `elemental_analysis.py`: concentration conversion, OLS/WLS calibration, detection limits, replicates, leaching |
| Plot infrastructure | `scripts/config.py`, `scripts/utils.py`, `styles/publication.mplstyle`, `styles/palettes.py` |
| Units | `scripts/preprocess.py`: `convert_energy`, `convert_temperature`, documented constants |
| CLI/examples | `scripts/main_plots.py`; `scripts/generate_synthetic_data.py`, `generate_electrochem_fixtures.py`, `generate_characterization_fixtures.py`, `generate_aas_fixtures.py` |

## Environment and tests

This is a root-import workspace, with no `src/` package layout or `pyproject.toml`. `pytest.ini` adds `.` to the Python path and selects `tests/`; `tests/conftest.py` selects Matplotlib Agg. Use `python -m scripts.main_plots` and `python -m pytest` from the root with the project interpreter.

`requirements.txt` specifies compatible ranges including NumPy 2.x. `requirements-tested.txt` pins the direct dependencies observed with Python 3.11.0; this is not a transitive lock or clean-install proof. `environment.yml` uses Python 3.11, Conda name `catalysis-workbench`, and that pip snapshot as one package source. The existing Windows venv directory `research-graphing` is intentionally unchanged. Reuse it; do not rebuild or install packages without authorization. See `docs/audits/STABILIZATION_REPORT.md` for the tested snapshot and limits.

Existing tests include domain test files, `tests/characterization/`, and `tests/test_scientific_audit.py`. Tests use `pytest.approx`, NumPy assertions, synthetic fixtures, and domain error checks. Passing these tests does not establish authentic-instrument compatibility. No repository formatter/linter configuration was found in the initial audit; follow surrounding style instead of introducing a tooling stack.

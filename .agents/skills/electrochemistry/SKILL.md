---
name: electrochemistry
description: Guides Metrohm Autolab NOVA and BioLogic EC-Lab imports and CV, LSV, EIS, CA/CP, Tafel, HER/OER, CHE, Cdl/ECSA, and catalyst benchmarking. Use for electrochemical analysis, parser changes, or instrument-derived figures.
---

# Electrochemistry

Use [instrument imports](references/instrument-imports.md) for raw files and [electrochemical conventions](references/electrochemical-conventions.md) before calculating derived quantities. Apply [scientific-validation](../scientific-validation/SKILL.md) for quantitative/parser work and [publication-figures](../publication-figures/SKILL.md) for figures.

1. Preserve the raw export. Reuse `parsers.electrochemistry.load_electrochemistry_data` and `StandardizedElectrochemicalData`; inspect channel mapping and instrument-specific defaults. Require confirmed reference electrode and filling solution, pH, temperature, area/loading, current sign, and prior iR compensation as needed by the calculation. Missing metadata must not become inferred experimental facts.
2. Choose the existing technique module in `modules/electrochemistry/`: `cv`, `lsv`, `eis`, `chronoamperometry`, `chronopotentiometry`, `tafel`, `ecsa`, `reference_electrodes`, or `benchmarking`. Computational CHE belongs in `modules/catalysis/reaction_energy.py`.
3. Keep potential reference, current normalization, units, cycle/sweep selection, electrolyte/gas, catalyst loading, scan rate, sampling, and processing history with results. Preserve sweep direction and elapsed-time ordering; do not sort a CV into an artificial monotonic curve.
4. Never double-apply iR correction. Verify amperes versus current density, compensation fraction, and whether correction happened in hardware/export. Derived corrected arrays need their own recorded state; the getter does not update the source dataset's compensation flag.
5. Inspect fit intervals, residuals, exclusions, and uncertainty. LSV does not automatically establish steady-state kinetics; a Tafel slope alone does not identify an HER/OER mechanism. A fitted EIS circuit does not uniquely identify a physical mechanism. Cdl-based ECSA requires a justified specific capacitance and is not a direct active-site count.
6. Compare catalysts only with compatible conditions, normalization, reference scales, resistance correction, and durability windows. Keep raw observations and assumptions distinct from activity/mechanistic interpretation.

Validate against `tests/test_metrohm.py`, `tests/test_biologic.py`, `tests/test_electrochemistry.py`, `tests/test_catalysis.py`, and the CHE cases in `tests/test_scientific_audit.py`. Synthetic exports establish functional regression only; report actual instrument/version coverage separately.

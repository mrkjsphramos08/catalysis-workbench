---
name: materials-characterization
description: Guides XRD, FTIR, Raman peak analysis, SEM/FESEM and TEM/HRTEM images, ICP-OES/MS, AAS, and BET/BJH work. Use for characterization imports, quantitative analysis, scientific review, or microscopy figures.
---

# Materials characterization

Start with the technique branch in [characterization methods](references/characterization-methods.md). Apply [scientific-validation](../scientific-validation/SKILL.md) for calculations and [publication-figures](../publication-figures/SKILL.md) for graphical output.

- Reuse parsers in `parsers/characterization/` and analysis in `modules/characterization/`. The shared container is `CharacterizationData` in `common.py`; gas isotherm import is `adsorption_isotherm.parse_isotherm_file`. Do not invent technique-specific schema classes from historical handoff descriptions.
- Preserve source arrays/images and hashes. Use `CharacterizationData.copy()` before processing and `add_processing_step()` for baseline, smoothing, cropping, normalization, calibration, and exclusions. Keep source metadata and processing history in result provenance.
- Reuse `common_peaks.py` for baseline estimation, cosmic-ray removal, smoothing, lmfit peaks, and band integration. Reuse `elemental_analysis.py` for concentration units, calibration, replicate statistics, detection limits, and leaching rather than duplicating ICP/AAS mathematics.
- Select processing based on measurement conditions and the scientific question. Inspect residuals and sensitivity to baseline/window/segmentation choices; good-looking fits do not establish a unique decomposition.
- Separate observed peak positions, intensities, sizes, and concentrations from assignments. Do not automatically assign crystal phases, chemical species, lattice planes, or molecular structures without justified reference data and stated uncertainty. Cross-technique correlation is not mechanistic causation.
- State unsupported operations plainly. Existing XRD peak fits are not Rietveld refinement; current microscopy FFT/fringe tools do not establish automatic crystallographic indexing; simplified BET/BJH routines are not a universal physisorption model.

Use `tests/characterization/test_<technique>.py` and the shared numerical tests affected by the change. Synthetic microscope images and spectra test algorithms, not instrument calibration or sample representativeness. Require calibration and genuine-data checks before claiming validated physical measurements.

---
name: scientific-python
description: Guides implementation and maintenance of this scientific Python toolkit. Use for Python APIs, numerical methods, NumPy or pandas processing, fitting, dependencies, and pytest work.
---

# Scientific Python

Read the [repository guidance](../../../AGENTS.md) and inspect the target function, imports, callers, and tests before creating a replacement.

1. Keep parsing in `parsers/`, analysis and domain plots in `modules/`, shared helpers in `scripts/`, and presentation defaults in `styles/`. Use the [architecture and API map](references/architecture.md) when locating reusable functionality.
2. Preserve public signatures, dataclass fields, array shapes, DataFrame columns, and result dictionary keys unless the task calls for a migration. Copy before modifying raw arrays/tables. Validate finite values, lengths, units, and required channels; do not silently drop observations.
3. Follow existing NumPy-style parameter/return docstrings and type hints; document physical units and reference scales. Use `pathlib.Path` and project paths from `scripts.config`, not machine-specific absolute paths.
4. Use NumPy for array operations and pandas for labeled tabular data. Make dtype conversion, index alignment, missing-value policy, and sorting explicit. Preserve acquisition order for hysteresis, CV cycles, and time series.
5. Reuse SciPy integration, optimization, and statistics as already used locally. Check grid order, sampling, constraints, convergence, and residuals. For peak fits reuse `common_peaks.fit_peaks` (lmfit); preserve bounds and report covariance/identifiability. statsmodels is declared and installed, but is not an existing project regression abstraction; use it only when a requested method needs it, after inspecting the installed API.
6. Verify unfamiliar third-party APIs with `inspect.signature`, installed docstrings/source and `importlib.metadata.version`, or official version-matched documentation. A declared dependency does not prove an API exists. Do not add frameworks or dependencies for functionality already available; keep both dependency manifests consistent when dependency changes are authorized.
7. Reuse energy conversions/constants in `scripts.preprocess`; unit-bearing names and explicit conversions remain essential even though Pint is available. Never conflate eV per entity with molar energies.
8. Preserve informative domain exceptions. Unsupported formats and missing physical parameters must remain visible; do not hide them behind broad exception handlers or plausible defaults.
9. Use pytest with `tmp_path`, numerical tolerances justified by the method, and exception/warning assertions. Pair scientific changes with [scientific-validation](../scientific-validation/SKILL.md); run focused tests then the existing suite. No new tests are needed merely to restate documentation.

For APIs consult [SciPy integration](https://docs.scipy.org/doc/scipy/reference/integrate.html), [lmfit models](https://lmfit.github.io/lmfit-py/builtin_models.html), and [statsmodels regression](https://www.statsmodels.org/stable/regression.html), checking their version against the active interpreter.

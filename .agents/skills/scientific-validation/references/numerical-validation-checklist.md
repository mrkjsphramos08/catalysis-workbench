# Numerical review checklist

Select cases relevant to the change; do not turn an unrelated documentation edit into a scientific rewrite.

| Area | Useful independent check | Failure/boundary case |
| --- | --- | --- |
| Units | One known conversion derived from documented constants; inverse round trip as an additional check | Unsupported units, mixed per-particle/per-mole normalization, double conversion |
| QE | Authentic output layout and known energies; PDOS component sum without LDOS | Incomplete/unconverged output, unsupported spin, mismatched energy grid, absent Fermi level |
| Adsorption/CHE | Hand-calculated stoichiometric energy difference; equivalent SHE/RHE potentials | Zero adsorbates, mismatched reference states, reversed reaction sign |
| DOS/charge integration | Analytic constant or polynomial on a known grid | Nonuniform/reversed grid, zero integral, empty window, NaN/Inf |
| Tafel/Arrhenius | Recover slope from a known synthetic relation with explicit log base and units | Nonpositive log inputs, too few points, transport/transient region |
| Calibration | Hand-computed slope/intercept and dilution example; independent uncertainty calculation | Rank-deficient standards, blank uncertainty, out-of-range unknowns, correlated replicates |
| XRD/spectra | Known peak position, area, FWHM and Bragg geometry | Unresolved instrumental width, overlapping peaks, sensitivity to baseline/smoothing |
| Microscopy | Calibrated known-size geometry and FFT frequency | Unknown pixel scale, resampling, merged particles, edge exclusions |
| BET/BJH | Known-model isotherm and cross-check with a justified reference method | Invalid pressure/volume, failed consistency checks, inappropriate adsorbate or pore regime |

## Existing evidence and gaps

- `tests/test_scientific_audit.py` checks constants, SHE/RHE and CHE signs, adsorption normalization, volcano cases, NEB barrier identities, d-band integration, and scaling diagnostics. Read its expected-value construction; a shared constant or algorithm can hide correlated errors.
- `tests/test_quantum_espresso.py` and `tests/fixtures/qe/FIXTURE_MANIFEST.md` provide synthetic/negative coverage. `tests/test_stabilization.py` additionally checks authentic upstream QE 6.0 NEB files, final-iteration values, explicit conditions, low-percent FTIR, provenance and export protection.
- `tests/test_electrochemistry.py`, instrument tests, and `tests/characterization/` exercise analysis and generated fixtures. Record instrument/version/calibration evidence separately.
- `pytest.ini` no longer filters `DeprecationWarning`. Report warnings without global suppression; compile source separately when checking invalid escapes because cached bytecode can hide compile-time warnings.

## Evidence record

For a material result record: source equation/method/version; explicit inputs and units; reference states; independent expected value; actual value; error and justified tolerance; fixture classification; uncertainty and applicability; commands/tests; and unresolved gaps. Use physical references for equations and installed library documentation for APIs. When they disagree with code, explain the discrepancy and scope a fix rather than quietly changing the reference.

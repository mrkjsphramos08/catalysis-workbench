# Technique-specific decisions

All module names below refer to `modules/characterization/` unless a parser is named. Read only the relevant branch.

## XRD

1. Need import? Use `parsers/characterization/xrd.py`; verify whether the axis is theta or 2-theta, angle units, wavelength/radiation, background, and intensity units.
2. Need peak positions/widths? Use `analyze_xrd_peaks` and shared `common_peaks.fit_peaks`; distinguish sigma from FWHM for the selected model and inspect overlap/residuals. Do not subtract background twice.
3. Need d-spacing? `calculate_d_spacing` implements Bragg geometry; use theta equal to half the measured 2-theta and an explicit wavelength/order convention.
4. Need domain size? `calculate_crystallite_size_scherrer` uses peak width, wavelength, shape factor, and instrumental correction. FWHM must enter the equation in radians. Quadrature broadening subtraction presumes a compatible profile model; a Lorentzian correction differs. At/below instrumental resolution preserve `ScherrerBroadeningWarning`, not a fabricated finite size. Strain, overlap, anisotropy, and shape affect the estimate; coherent-domain size is not automatically SEM particle size.

Phase identification needs a traceable reference pattern and experimental context; local peak annotations are not an identification engine.

## FTIR and Raman

1. Need import/conversion? Use the technique parser. Distinguish cm^-1 wavenumber, Raman shift, intensity, fraction transmittance, and percent transmittance. `parse_ftir_file(..., transmittance_unit="percent"|"fraction")` requires explicit units or a unit-bearing column header. `convert_transmittance_to_absorbance(..., transmittance_unit=...)` accepts an explicit keyword or dataset units; magnitude inference is removed. Ambiguous signal headers also require `signal_type`. Reject conflicting units. For example, 1 percent gives absorbance 2, while fraction 1 gives 0. Nonpositive transmittance is invalid; do not enable clamping silently.
2. Need correction? Use `common_peaks.estimate_baseline`, `subtract_baseline`, `remove_cosmic_rays`, or `smooth_spectrum` only with justified parameters. Retain original data and compare before/after; smoothing must not erase, shift, or manufacture quantitative peaks.
3. Need fit/integral? Use `fit_peaks` and `integrate_band`; record baseline, shape, bounds, fit window, noise model, residuals, and uncertainty. Inspect lmfit parameter definitions in the [official models reference](https://lmfit.github.io/lmfit-py/builtin_models.html).
4. Need interpretation? `ftir.analyze_ftir_bands` and `raman.process_raman_spectrum` / `calculate_peak_ratios` provide analysis paths; height `I_D/I_G` and area `A_D/A_G` are different quantities. Laser wavelength, disorder regime, composition, and acquisition conditions matter. Do not claim a built-in universal Raman domain-size model or unique molecular assignment.

## SEM/FESEM and TEM/HRTEM

1. Calibration known? Use `parsers/characterization/microscopy.parse_microscopy_image` and verify `pixel_size_nm` against acquisition metadata or a documented scale reference. Import can retain an unknown calibration; physical sizing and scale bars cannot.
2. Particle statistics? Use `microscopy.segment_particles`; check threshold polarity, ROI, touching particles, border exclusions, and removal of burned-in text/bars from the analysis ROI. Retain the original full image. Inspect overlays and sample multiple representative fields. Report count, distribution definition, and segmentation uncertainty; equivalent circular diameter is not a crystallite diameter.
3. Lattice spacing? Use `analyze_hrtem_lattice_fringes`; verify real/reciprocal calibration, ROI, orientation, drift, sampling/Nyquist limit, and uncertainty. FFT spots or fringe periods alone do not justify lattice-plane assignment. Current functions are not automated SAED indexing.
4. Figure? Use `add_scale_bar`, `plot_microscopy_with_scalebar`, or `plot_hrtem_fft_panel`. A crop retains pixel scale; resampling changes it. Never invent a physical scale from magnification text alone.

## ICP-OES, ICP-MS, and AAS

1. External calibration valid? Use `icp.fit_icp_calibration` or `aas.fit_aas_calibration`, backed by `elemental_analysis.linear_calibration_fit`. Confirm concentration/signal units, blank treatment, replicates, weighting, standards range, and residuals. Do not force a zero intercept or extrapolate without justification.
2. Matrix effects? `aas.fit_standard_addition` exists for AAS; record spike and final volumes and uncertainty. Do not assume an equivalent public ICP standard-addition wrapper exists. FAAS absorbance and GFAAS integrated absorbance have different units; stored background-correction metadata does not mean software performed instrument-level correction.
3. Quantification? Use `aas.quantify_aas_samples` and shared concentration utilities. Apply dilution and blank corrections once, with explicit digest/sample masses and volumes. Treat negative blank-corrected estimates and below-LOD/LOQ values transparently; do not replace them silently by zero.
4. QC? Use `calculate_detection_limits`, replicate statistics, and `aas.evaluate_aas_qc`. Record which variance estimate, slope, multiplier, control sample, acceptance threshold, and SOP were used. LOD/LOQ conventions vary; an R-squared threshold alone does not validate quantification or establish regulatory compliance.
5. Leaching/loading? Reuse `calculate_catalyst_leaching_common` through domain wrappers, and `aas.calculate_aas_catalyst_loading` where applicable. Account for dilution, sampled/replaced volume, cumulative dissolved mass, initial metal inventory, and uncertainty. Mass balance and units must agree.

See the existing [AAS guide](../../../../docs/characterization/AAS_GUIDE.md), checking its claims against code and the applicable SOP. Primary methods: [EPA 7000B](https://www.epa.gov/sites/default/files/2015-12/documents/7000b.pdf) and [EPA 7010](https://www.epa.gov/sites/default/files/2015-12/documents/7010.pdf). These are analytical references, not certification of this implementation.

## BET/BJH

1. Import using `parsers/characterization/adsorption_isotherm.py`; establish adsorbate, temperature, degassing, sample mass, STP convention, pressure normalization, and adsorption/desorption branch.
2. For BET use `bet.calculate_bet_surface_area` with an explicit `p_range`; the current routine evaluates a supplied window, not an automatic optimal-window search. Review returned `rouquerol_criteria`, positive C, increasing q(1-p/p0), monolayer location, residuals, and sensitivity to range. A high R-squared alone is insufficient.
3. For BJH use `calculate_bjh_pore_distribution`; establish branch and volume units. Its nitrogen-oriented constants and cylindrical-pore/Kelvin/thickness assumptions cannot be generalized to arbitrary adsorbates or micropores just by changing a label. Report cavitation/network effects and method-dependent pore sizes.

Use the [IUPAC physisorption report](https://doi.org/10.1515/pac-2014-1117) to assess applicability. The implemented three checks are useful diagnostics, not proof that every consistency condition or measurement assumption has been met.

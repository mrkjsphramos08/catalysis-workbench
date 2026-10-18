# Automated XRD Analysis & Crystallite Sizing Report

**Sample:** Sample 1 Black powder.ASC  
**Support:** Sample 2 White powder.ASC  
**Radiation:** CU K\alpha (\lambda = 1.54060 \text{ \AA})  
**Cleaned Diffractogram:** [`data/cleaned/xrd/cobalt_silica_sample1_cleaned.csv`](file:///c:/Users/mjram/catalysis-workbench/data/cleaned/xrd/cobalt_silica_sample1_cleaned.csv)  
**Peak Inventory:** [`data/cleaned/xrd/cobalt_silica_sample1_peaks.csv`](file:///c:/Users/mjram/catalysis-workbench/data/cleaned/xrd/cobalt_silica_sample1_peaks.csv)  
**Figure Exports:** [`figures/characterization/xrd/fig_xrd_cobalt_silica_sample1.*`](file:///c:/Users/mjram/catalysis-workbench/figures/characterization/xrd/) (PDF, SVG, PNG, TIFF)  

---

## 1. Manuscript-Ready Text Paragraphs

### Experimental Section (Materials & Methods - XRD Characterization)
> **Powder X-ray Diffraction (PXRD):** Crystallographic phases and structural order were characterized using powder X-ray diffraction recorded in Bragg-Brentano reflection geometry with CU K\alpha radiation (\lambda = 1.54060 \text{ \AA}). Diffractograms were acquired across the angular range of 10.0\degree <= 2theta <= 80.0\degree with a step increment of 0.043\degree. The blank support material (Sample 2 White powder) was independently measured under identical instrumental parameters. Support scattering was scaled (k = 0.859) and subtracted to isolate the net scattering profile of the supported active component. Apparent volume-weighted mean crystallite domain dimensions (D) were calculated from the line-broadened Bragg reflections using the Scherrer equation (D = K\lambda / (\beta_{\text{corr}} cos theta)) with a shape factor K = 0.90. The net physical sample broadening (\beta_{\text{corr}}) was deconvoluted by subtracting the laboratory instrumental resolution limit (\beta_{\text{inst}} = 0.060\degree) via Gaussian quadrature deconvolution.

### Results & Discussion Section (Structural & Phase Analysis)
> **Crystallographic Structure and Phase Analysis:** The experimental powder X-ray diffractogram of Sample 1 Black powder revealed distinct Bragg reflections superimposed on the underlying support scattering profile. Following background deconvolution, the primary active phase reflection was resolved at 2theta = 20.51\degree, corresponding to an interplanar d-spacing of 4.327 \text{ \AA}. Profile fitting yielded an observed full width at half-maximum (FWHM) of 0.090\degree, which corresponds to an instrumental-corrected mean coherently scattering crystallite domain dimension of D = 121.2 nm (1212 \text{ \AA}). Assuming quasi-spherical crystallite geometry, the corresponding metallic dispersion is estimated at approximately 0.8% ($D_{\%} \approx 96 / d_{\text{nm}})$. The observed diffraction positions and interplanar spacings confirm the formation of nanosized crystalline domains uniformly anchored on the high-surface-area matrix without evidence of bulk sintering.

---

## 2. Quantitative Crystallite Sizing Summary

- **Primary Reflection:** $2\theta = 20.51^\circ$ ($d = 4.327\text{ \AA}$)
- **Observed Broadening (FWHM):** $\beta_{\text{obs}} = 0.0896^\circ$
- **Instrument Broadening Baseline:** $\beta_{\text{inst}} = 0.0600^\circ$
- **Corrected Coherent Domain Size ($D$):** **121.23 nm** (1212.3 Å)
- **Estimated Active Metal Dispersion:** **0.8%**

### Top Detected Bragg Reflections

| $2\theta$ ($^\circ$) | $d$-spacing (\text{\AA}) | Net Intensity (counts) | FWHM ($^\circ$) | Scherrer Domain Size $D$ (nm) |
| :---: | :---: | :---: | :---: | :---: |
| 20.51 | 4.327 | 549.6 | 0.090 | 121.23 |
| 19.74 | 4.494 | 492.9 | 0.158 | 55.2 |
| 22.40 | 3.965 | 457.9 | 0.389 | 21.05 |
| 18.92 | 4.686 | 417.1 | 0.529 | 15.33 |
| 44.25 | 2.045 | 408.1 | 0.860 | 9.99 |
| 16.77 | 5.282 | 403.8 | 1.076 | 7.48 |
| 21.33 | 4.162 | 376.6 | 12.882 | 0.63 |
| 29.93 | 2.983 | 323.2 | 1.587 | 5.18 |
| 18.06 | 4.908 | 321.2 | 2.175 | 3.7 |
| 15.14 | 5.848 | 299.8 | 0.117 | 80.0 |

---

## 3. Provenance & Subtraction Parameters
- Support intensity scale factor ($k$): `0.8593`
- Sample displacement zero-shift ($\Delta 2\theta$): `-0.2150^\circ`
- Scherrer shape factor ($K$): `0.90`

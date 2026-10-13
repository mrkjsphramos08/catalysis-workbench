# python_graphing Context & Domain Vocabulary

## Domain Overview
`python_graphing` is a modular, publication-quality scientific graphing workspace designed for computational chemistry, density functional theory (DFT), heterogeneous catalysis, and chemical engineering. It is a post-processing and visualization toolkit that ingests raw/processed calculation outputs and experimental data, enforces thermodynamic and kinetic invariants, and exports vector/raster figures compliant with high-impact journal standards (e.g., Nature, JACS, ACS Catalysis, Journal of Catalysis, AIChE Journal).

## Core Modules & Vocabulary
- **Quantum ESPRESSO Parsers (`parsers/quantum_espresso/`)**:
  - `pw.x`: SCF and structural relaxation (`relax`, `vc-relax`). Extracts total energy in Ry and eV ($1\text{ Ry} = 13.605693009\text{ eV}$), ionic step history, forces ($\text{Ry/Bohr} \to \text{eV/\AA}$), stresses, Fermi energy, SCF and BFGS convergence status, and job completion (`JOB DONE.`).
  - `dos.x`: Total DOS parsing (`prefix.dos`), non-spin-polarized and collinear spin-polarized channels, integrated DOS, Fermi alignment ($E - E_F$).
  - `projwfc.x`: Projected DOS parsing from `prefix.pdos_tot` and `prefix.pdos_atm#N(Element)_wfc#M(Orbital)`. Correct component ordering: $p_z, p_x, p_y$ and $d_{z^2}, d_{zx}, d_{zy}, d_{x^2-y^2}, d_{xy}$. Explicit exclusion of LDOS from orbital sum to prevent double counting. Rejection of unsupported noncollinear/spin-orbit formats.
  - `bands.x`: Electronic band structure parsing from Gnuplot `.gnu` outputs (`prefix.bands.dat.gnu`).
  - `neb.x`: Nudged Elastic Band parsing (`prefix.dat`, `prefix.int`, `prefix.path`) with image coordinates, relative free energies, and forward/reverse activation barriers.
- **DFT & Electronic Structure (`modules/dft`)**:
  - `E_F` / Fermi level alignment: $E - E_F = 0\text{ eV}$.
  - Total DOS (TDOS) and Projected DOS (PDOS): Orbital projections ($s, p, d$, $t_{2g}, e_g$, etc.), spin-polarized channels (spin-up $>0$, spin-down $<0$).
  - Band Structure: Energy vs reciprocal-space k-path ($\Gamma, K, M, X, \dots$), Fermi reference.
  - Adsorption Energy ($E_{\text{ads}}$): $E_{\text{ads}} = \frac{E(\text{slab} + \text{adsorbate}) - E(\text{slab}) - E(\text{ref})}{n_{\text{ads}}}$. Configurable sign conventions and references.
- **Heterogeneous Catalysis (`modules/catalysis`)**:
  - Reaction Energy Profile: Elementary steps, reaction free energy ($\Delta G$), activation barriers ($\Delta G^\ddagger$ / TS), Computational Hydrogen Electrode (CHE) model: $\Delta G(U, \text{pH}) = \Delta G_0 - eU + k_B T \ln(10) \cdot \text{pH}$.
  - Sabatier Volcano Plots: Distinct physical descriptors for Oxygen Reduction Reaction (ORR, descriptor $\Delta G_{\mathrm{OH}^*}$ or $\Delta G_{\mathrm{O}^*}$) vs Hydrogen Evolution Reaction (HER, descriptor $\Delta G_{\mathrm{H}^*}$).
- **Chemical Engineering & Kinetics (`modules/chemical_engineering`)**:
  - Arrhenius Analysis: $\ln(k) = \ln(A) - \frac{E_a}{RT}$, rate constants vs rates, linear regression diagnostics ($R^2$, residuals, confidence intervals).
  - Reactor Performance: Conversion ($X_A$), selectivity ($S_P$), and yield ($Y_P$) vs temperature and residence time / space velocity.
- **Plotting Infrastructure (`scripts/config.py`, `scripts/utils.py`, `styles/publication.mplstyle`, `styles/palettes.py`)**:
  - Single-column (85 mm ~ 3.35 in) and double-column (175 mm ~ 6.89 in) presets.
  - 300+ DPI raster (PNG, TIFF) and vector export (PDF, SVG).
  - Colorblind-safe palettes (Tol, Okabe-Ito, ColorBrewer).
  - Provenance metadata tracking (`*_metadata.json` with SHA-256 hashes and physical units).

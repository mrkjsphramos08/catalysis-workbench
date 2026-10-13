# catalysis-workbench: Computational Chemistry, DFT & Catalysis Research Workspace

`catalysis-workbench` is a modular, publication-quality scientific graphing and post-processing workspace engineered for computational chemistry, density functional theory (DFT), heterogeneous catalysis, and chemical engineering.

Its primary purpose is to process research datasets, analyze computational and experimental results, and generate scientifically rigorous, journal-ready figures adhering to international publishing standards (e.g., Nature Portfolio, ACS, RSC, Elsevier, Wiley, and AIChE).

---

## 1. Project Architecture

```text
catalysis-workbench/
├── data/
│   ├── raw/                             # Pristine raw calculations and experimental datasets (never modified)
│   ├── processed/                       # Standardized, unit-converted datasets with traceable provenance
│   └── examples/                        # Deterministic, documented synthetic benchmark datasets
│       ├── electrochemistry/
│       │   ├── metrohm/                 # Metrohm Autolab NOVA exports (ASCII, CSV, FRA EIS, CA/CP)
│       │   └── biologic/                # BioLogic EC-Lab exports (.mpt, .mps, .mpr)
│
├── parsers/                             # Simulation and experimental instrument parsers
│   ├── __init__.py
│   ├── quantum_espresso/                # Quantum ESPRESSO parsers (pw.x, dos.x, projwfc.x, bands.x, neb.x)
│   └── electrochemistry/                # Potentiostat data parsers
│       ├── __init__.py                  # Unified load_electrochemistry_data dispatcher
│       ├── schema.py                    # StandardizedElectrochemicalData with IUPAC sign conventions
│       ├── column_mapping.py            # Regex channel mapping and unit conversions
│       ├── metrohm.py                   # Metrohm Autolab NOVA ASCII/CSV parser (.nox rejection)
│       └── biologic.py                  # BioLogic EC-Lab .mpt/.mps/.mpr parser
│
├── scripts/
│   ├── __init__.py
│   ├── config.py                        # Centralized paths, journal dimensions, DPI, and export presets
│   ├── utils.py                         # Style loading, figure creation, export, and schema validation
│   ├── preprocess.py                    # Physical unit conversion (eV, Ry, kJ/mol, kcal/mol, K, C) and pipeline
│   ├── generate_synthetic_data.py       # Deterministic generator for benchmark research data
│   ├── generate_electrochem_fixtures.py # Generator for synthetic potentiostat benchmark data
│   ├── build_notebooks.py               # Builder for verified interactive research notebooks
│   └── main_plots.py                    # CLI runner for batch benchmarks, QE, and electrochem figures
│
├── modules/
│   ├── __init__.py
│   ├── dft/                             # Density Functional Theory post-processing (DOS, PDOS, bands, scaling)
│   ├── catalysis/                       # Free energy pathways, CHE models, and Sabatier volcano curves
│   ├── chemical_engineering/            # Arrhenius regression, reaction kinetics, and PFR/CSTR profiles
│   └── electrochemistry/                # Experimental electrochemistry analysis & kinetics
│       ├── __init__.py
│       ├── reference_electrodes.py      # Reference database (SHE/RHE/SCE/Ag-AgCl/Hg-HgO), dE/dT, iR drop
│       ├── cv.py                        # Multi-cycle CV, redox peak detection (E_pa, E_pc, Delta E_p), charge Q
│       ├── lsv.py                       # HER/OER polarization, overpotential extraction (eta_10, eta_50, eta_100)
│       ├── tafel.py                     # Tafel slope regression (mV/dec), j_0, NonSteadyStateWarning
│       ├── ecsa.py                      # Double-layer capacitance C_dl, ECSA, roughness factor, CdlInterpretationWarning
│       ├── eis.py                       # Randles circuit fitting (Rs - [Rct || CPE]), Lin-KK test, Nyquist/Bode
│       ├── chronoamperometry.py         # Potentiostatic stability, current retention %, degradation rate %/h
│       ├── chronopotentiometry.py       # Galvanostatic durability, potential drift rate mV/h
│       └── benchmarking.py              # Cross-instrument benchmarking (Metrohm vs BioLogic), delta j(E)
│
├── figures/                             # Domain-categorized publication figure outputs
│   ├── dft/
│   ├── catalysis/
│   ├── chemical_engineering/
│   └── electrochemistry/
│
├── styles/
│   ├── publication.mplstyle             # Master journal style sheet (sans-serif, inward ticks, 300+ DPI)
│   └── palettes.py                      # Colorblind-safe palettes (Okabe-Ito, Tol Bright/Muted, Orbitals)
│
├── notebooks/
│   ├── exploratory.ipynb                # Interactive Jupyter walkthrough of workflows and figures
│   ├── metrohm_analysis.ipynb           # Metrohm NOVA OER, CV, EIS Nyquist, and CA durability
│   ├── biologic_analysis.ipynb          # BioLogic EC-Lab HER, PEIS Bode, Lin-KK, and CP drift
│   └── electrochemical_comparison.ipynb # Inter-instrument benchmark and residual comparison
│
├── outputs/                             # Centralized multi-format publication exports
│   ├── vector/                          # Editable vector graphics (PDF, SVG)
│   └── raster/                          # High-resolution raster graphics (PNG, TIFF)
│
├── tests/                               # Comprehensive pytest test suite (90 passing tests)
│   ├── fixtures/qe/                     # Authentic redistributable Quantum ESPRESSO test fixtures
│   └── fixtures/electrochemistry/       # Metrohm and BioLogic benchmark fixtures
├── environment.yml                      # Conda / Mamba environment specification (Python 3.12)
├── requirements.txt                     # Pip dependency specification
├── CONTEXT.md                           # Project domain vocabulary and scientific invariants
├── .gitignore
└── README.md
```

---

## 2. Environment Setup & Installation

### Option A: Conda / Miniforge (Recommended)
```bash
conda env create -f environment.yml
conda activate research-graphing
```

### Option B: Python Virtual Environment (venv / pip)
```bash
python -m venv research-graphing
# On Windows PowerShell:
.\research-graphing\Scripts\Activate.ps1
# On Linux / macOS:
source research-graphing/bin/activate

pip install --upgrade pip
pip install -r requirements.txt
```

### Dependency Notes & Limitations
- **Python**: Engineered and verified for Python 3.11 and 3.12.
- **Cantera**: Optional dependency for detailed microkinetic reactor network integration; not required for core graphing and kinetic regression.
- **Quantum ESPRESSO**: The workspace ingests outputs from Quantum ESPRESSO v5.x, v6.x, and v7.x (`pw.x`, `dos.x`, `projwfc.x`, `bands.x`, `neb.x`) and does not require local QE compiler/MPI binaries.

---

## 3. Publication Standards & Design Principles

All generated figures adhere to international journal production guidelines:
- **Dimensions**:
  - Single column: 85 mm (3.35 in)
  - 1.5 column: 115 mm (4.53 in)
  - Double column: 175 mm (6.89 in)
- **Typography**: Clear sans-serif hierarchy (Arial / DejaVu Sans), 8–9 pt labels, 7.5 pt legends and tick labels.
- **Color Accessibility**: Default palettes use Okabe-Ito and Paul Tol colorblind-safe color cycles.
- **Tick & Frame Styling**: Inward-facing ticks (`xtick.direction: in`), top/right spines enabled, clean white background.
- **Vector Preservation**: TrueType fonts embedded as Type 42 (`pdf.fonttype = 42`, `ps.fonttype = 42`, `svg.fonttype = none`), ensuring fully editable text and lines in Adobe Illustrator or Inkscape.
- **Raster Resolution**: 300 DPI standard (configurable up to 600+ DPI).
- **Data Provenance**: Automatic companion `*_metadata.json` files storing SHA-256 hashes of original calculation files, physical units, Fermi alignments, and timestamps.

---

## 4. Quantum ESPRESSO Direct Plotting (CLI)

Plot real Quantum ESPRESSO results directly from the command line:

```bash
# 1. Total DOS from dos.x:
python -m scripts.main_plots --type dos --input path/to/prefix.dos --fermi 6.54

# 2. Projected DOS from projwfc.x (auto-discover all wavefunction files):
python -m scripts.main_plots --type pdos --input path/to/pdos_dir --prefix prefix --fermi 6.54

# 3. Projected DOS with selective element and orbital filtering:
python -m scripts.main_plots --type pdos --input path/to/pdos_dir --prefix prefix --fermi 6.54 \
    --elements Pt O --orbitals d p --energy-range -7.0 5.0 --output my_pt_o_pdos

# 4. Band structure from bands.x:
python -m scripts.main_plots --type bands --input path/to/prefix.bands.dat.gnu --fermi 6.54

# 5. NEB reaction free energy profile from neb.x:
python -m scripts.main_plots --type neb --input path/to/prefix.dat
```

---

## 5. Quantum ESPRESSO Python API & Worked Examples

### Example 1: Plotting Total DOS from `dos.x`
```python
from parsers.quantum_espresso import parse_dos_file, plot_qe_dos
from scripts.utils import create_publication_figure, save_publication_figure

# Parse native QE DOS file (supports non-spin and spin-polarized collinear)
dos_data = parse_dos_file("data/raw/prefix.dos", e_fermi=6.54)

fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
plot_qe_dos(dos_data, ax=ax, energy_range=(-6.0, 6.0))

save_publication_figure(fig, "qe_total_dos", formats=("pdf", "svg", "png", "tiff"))
```

### Example 2: Plotting Metal d-Orbital PDOS and d-Band Center
```python
from parsers.quantum_espresso import QEPdosCollection, plot_qe_pdos
from scripts.utils import create_publication_figure, save_publication_figure

# Discover all prefix.pdos_atm#... files in directory
coll = QEPdosCollection("data/raw/pdos/", prefix="prefix", e_fermi=6.54)

# Calculate d-band center descriptor (epsilon_d)
eps_d = coll.calculate_d_band_center(element="Pt", energy_window=(-8.0, 0.0), fermi_energy=6.54)
print(f"Pt d-band center: {eps_d:.2f} eV")

# Define target orbital components (avoids LDOS double-counting)
specs = {
    r"Pt $5d$ ($e_g$)": {"element": "Pt", "orbital": "d", "component": "eg"},
    r"Pt $5d$ ($t_{2g}$)": {"element": "Pt", "orbital": "d", "component": "t2g"},
}

fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
plot_qe_pdos(coll, specs=specs, ax=ax, energy_range=(-6.0, 4.0))
save_publication_figure(fig, "pt_pdos_components")
```

### Example 3: Comparing Surface and Adsorbate Electronic Properties
```python
from parsers.quantum_espresso import QEPdosCollection, plot_qe_pdos
from scripts.utils import create_publication_figure

coll = QEPdosCollection("data/raw/co_pt111/", prefix="co_pt", e_fermi=7.10)

# Compare substrate d-band and adsorbate pi/sigma states
specs = {
    r"Pt $5d$": {"element": "Pt", "orbital": "d"},
    r"C $2p$":  {"element": "C",  "orbital": "p"},
    r"O $2p$":  {"element": "O",  "orbital": "p"},
}

fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
plot_qe_pdos(coll, specs=specs, ax=ax, energy_range=(-8.0, 4.0))
```

### Example 4: Plotting Band Structures from `bands.x`
```python
from parsers.quantum_espresso import parse_bands_gnu, plot_qe_bands
from scripts.utils import create_publication_figure, save_publication_figure

# High-symmetry k-points along path (user-provided coordinates and labels)
hs_points = [(0.0, r"$\Gamma$"), (0.866, "$K$"), (1.0, "$M$"), (1.577, r"$\Gamma$")]

bands_data = parse_bands_gnu(
    "data/raw/prefix.bands.dat.gnu",
    e_fermi=6.54,
    high_symmetry_points=hs_points,
)

fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
plot_qe_bands(bands_data, ax=ax, energy_range=(-4.0, 4.0))
save_publication_figure(fig, "qe_band_structure")
```

### Example 5: Extracting Adsorption Energies from `pw.x` Outputs
```python
from parsers.quantum_espresso import parse_pw_output, calculate_adsorption_energy_from_qe

# Parse clean surface and surface+adsorbate SCF/relaxation outputs
slab_out = parse_pw_output("data/raw/clean_pt111.out", strict=True)
slab_ads_out = parse_pw_output("data/raw/co_pt111.out", strict=True)

# Gas-phase CO reference energy in eV (from isolated molecule calculation)
e_co_gas = -14.28

# Calculate adsorption energy: E_ads = E(tot) - E(slab) - E(ref)
e_ads = calculate_adsorption_energy_from_qe(
    slab_output=slab_out,
    slab_adsorbate_output=slab_ads_out,
    adsorbate_ref_ev=e_co_gas,
    n_adsorbates=1,
    sign_convention="negative_favorable",
)
print(f"CO adsorption energy on Pt(111): {e_ads:.2f} eV")
```

### Example 6: Plotting NEB Activation Barriers from `neb.x`
```python
from parsers.quantum_espresso import parse_neb_output, plot_qe_neb_profile
from scripts.utils import create_publication_figure, save_publication_figure

# Parses prefix.dat, prefix.int, and prefix.path automatically
neb_data = parse_neb_output("data/raw/co_oxidation.dat")

print(f"Forward barrier: {neb_data.forward_barrier_ev:.2f} eV")
print(f"Reverse barrier: {neb_data.reverse_barrier_ev:.2f} eV")

fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
plot_qe_neb_profile(neb_data, ax=ax, title="CO Oxidation Barrier on Pt(111)")
save_publication_figure(fig, "co_oxidation_neb")
```

---

## 6. Quantum ESPRESSO File Retention Guide

Researchers should retain the following primary outputs from their Quantum ESPRESSO workflows:

| Engine | Required File | Description | Ingestion Function |
| :--- | :--- | :--- | :--- |
| `pw.x` | `prefix.out` | Electronic SCF or relaxation output log | `parse_pw_output()` |
| `dos.x` | `prefix.dos` | Total density of states | `parse_dos_file()` |
| `projwfc.x` | `prefix.pdos_tot` | Total DOS from atomic projection | `QEPdosCollection()` |
| `projwfc.x` | `prefix.pdos_atm#*` | Atomic and orbital projected wavefunctions | `QEPdosCollection()` |
| `bands.x` | `prefix.bands.dat.gnu` | Dispersion curve table formatted for Gnuplot | `parse_bands_gnu()` |
| `neb.x` | `prefix.dat` | Discrete image reaction coordinates and energies | `parse_neb_output()` |
| `neb.x` | `prefix.int` | Splined interpolated reaction pathway | `parse_neb_output()` |
| `neb.x` | `prefix.path` | Summary barrier heights and convergence status | `parse_neb_output()` |

*Note: Raw wavefunction binaries (`.wfc`), charge densities (`.charge-density`), and restart files (`.save/`) are large scratch files and do not need to be retained for graphing.*

---

## 7. Official Quantum ESPRESSO Documentation References

- **pw.x Input Reference**: https://www.quantum-espresso.org/Doc/INPUT_PW.html
- **dos.x Input Reference**: https://www.quantum-espresso.org/Doc/INPUT_DOS.html
- **projwfc.x Input Reference**: https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html
- **bands.x Input Reference**: https://www.quantum-espresso.org/Doc/INPUT_BANDS.html
- **neb.x Input Reference**: https://www.quantum-espresso.org/Doc/INPUT_NEB.html

---

## 8. Experimental Electrochemistry System (Metrohm Autolab & BioLogic)

The electrochemistry suite provides an end-to-end framework for ingesting, validating, standardizing, modeling, and publishing experimental data acquired on **Metrohm Autolab** (NOVA) and **BioLogic** (EC-Lab) potentiostats.

### 8.1 Instrument Compatibility & File Formats
- **Metrohm Autolab (NOVA)**:
  - Supported: Plain ASCII (`.txt`), delimited comma/semicolon/tab tables (`.csv`), and Excel spreadsheets (`.xlsx`).
  - Binary `.nox` handling: Safely rejected with clear instructions guiding users to export to ASCII or CSV from NOVA.
  - Channels parsed: Potential ($E_{\mathrm{we}}$), Current ($I$), Time ($t$), Frequency ($f$), $Z'$, $Z''$, Scan/Cycle numbers, electrode geometric area, uncompensated resistance ($R_u$), reference electrode type, and pH.
- **BioLogic (EC-Lab)**:
  - Supported: ASCII tabular data (`.mpt`) with dynamic header detection (`Nb header lines`), experiment settings (`.mps`), and binary raw data (`.mpr`) via `galvani`.
  - Channels parsed: Potential (`Ewe/V`), Current (`I/mA`), Time (`time/s`), Frequency (`freq/Hz`), `Re(Z)/Ohm`, `-Im(Z)/Ohm`, `|Z|/Ohm`, `Phase(Z)/deg`, Cycle (`cycle number`).

### 8.2 Standardized Schema & IUPAC Sign Convention
All instrument data parse into `StandardizedElectrochemicalData`, enforcing:
- **IUPAC Current Convention**: Anodic (oxidation) currents are strictly positive ($+I$), and cathodic (reduction) currents are strictly negative ($-I$).
- **Normalized Densities**: Flexible extraction via `.get_current_density(unit="mA/cm^2")` or geometric mass loading.
- **Reference Conversions**: Built-in conversion to RHE via `.get_potential_vs_rhe(reference, ph)` incorporating standard potentials and temperature coefficients ($dE/dT$).
- **Software $iR$ Compensation**: `.get_ir_corrected_potential(ru_ohm, fraction)` with double-compensation guard to prevent overcorrection when hardware compensation was active.

### 8.3 Core Analysis Capabilities
1. **Voltammetry & Overpotential Extraction**:
   - `modules/electrochemistry/cv.py`: Redox peak detection ($E_{pa}, E_{pc}, \Delta E_p, E_{1/2}$) and voltammetric charge integration ($Q = \int I dt$).
   - `modules/electrochemistry/lsv.py`: HER and OER polarization curves, overpotential extraction at target current densities ($\eta_{10}, \eta_{50}, \eta_{100}$).
2. **Kinetics & Active Surface Area**:
   - `modules/electrochemistry/tafel.py`: Quasi-steady-state Tafel slopes ($\mathrm{mV/dec}$) and exchange current densities ($j_0$). Emits `NonSteadyStateWarning` if scan rate $\nu > 5\,\mathrm{mV/s}$.
   - `modules/electrochemistry/ecsa.py`: Double-layer capacitance ($C_{dl}$) from scan-rate series, ECSA, and roughness factors. Emits `CdlInterpretationWarning` reminding users that $C_{dl}$ reflects electrostatic capacitance rather than active-site density.
3. **Impedance Spectroscopy & Durability**:
   - `modules/electrochemistry/eis.py`: Complex impedance polar conversion, Randles equivalent circuit fitting ($R_s - [R_{ct} \parallel \mathrm{CPE}]$), orthonormal Nyquist plots, 2-panel stacked Bode plots, and Lin-KK Kramers-Kronig data validation using Schönleber's criterion.
   - `modules/electrochemistry/chronoamperometry.py`: Potentiostatic stability, current retention percentage, and hourly degradation rate ($\%/h$).
   - `modules/electrochemistry/chronopotentiometry.py`: Galvanostatic potential drift rate ($\mathrm{mV/h}$).
4. **Inter-Instrument Benchmarking**:
   - `modules/electrochemistry/benchmarking.py`: Aligns and compares matched catalysts across Metrohm and BioLogic, extracting point-by-point residual current discrepancies $\Delta j(E) = j_{\mathrm{Metrohm}} - j_{\mathrm{BioLogic}}$ and RMSE.

### 8.4 Command-Line Figure Generation
Generate all 10 experimental electrochemistry publication figures in vector (PDF, SVG) and raster (PNG, TIFF) formats:

```bash
# Generate only electrochemistry figures
python -m scripts.main_plots --category electrochemistry

# Generate all workspace figures (DFT, Catalysis, Chemical Engineering, Electrochemistry)
python -m scripts.main_plots --all
```

### 8.5 Interactive Research Notebooks
- [`notebooks/metrohm_analysis.ipynb`](notebooks/metrohm_analysis.ipynb): Metrohm NOVA OER polarization, multi-cycle CV, FRA EIS Nyquist modeling, and CA durability.
- [`notebooks/biologic_analysis.ipynb`](notebooks/biologic_analysis.ipynb): BioLogic EC-Lab HER polarization, PEIS Bode analysis, Lin-KK validation, and CP stability.
- [`notebooks/electrochemical_comparison.ipynb`](notebooks/electrochemical_comparison.ipynb): Cross-instrument benchmark between Metrohm and BioLogic with residual analysis and kinetic concordance tables.

---

## 9. Materials Characterization System (Phase 4)

The workspace provides an autonomous, publication-grade analytical and visualization framework for materials characterization across 7 core techniques:

### 9.1 Supported Techniques & Analytical Capabilities
1. **X-Ray Diffraction (XRD)**:
   - Ingestion: Multi-column CSV and space/comment-delimited diffractograms (`.csv`, `.txt`, `.xy`, `.dat`).
   - Bragg $d$-Spacing: $\lambda = 2d\sin\theta$ for Cu K$\alpha$ ($\lambda = 1.54060\text{ \AA}$) and arbitrary anodes.
   - Scherrer Crystallite Domain Size: Root-sum-square instrumental broadening correction $\beta_{\mathrm{corr}} = \sqrt{\beta_{\mathrm{obs}}^2 - \beta_{\mathrm{inst}}^2}$, with `ScherrerBroadeningWarning` when $\beta_{\mathrm{obs}} \le \beta_{\mathrm{inst}}$.
   - Multi-Peak Deconvolution: Pseudo-Voigt, Gaussian, and Lorentzian composite profile fitting via `lmfit` with residual analysis.
   - Diffractogram Visualizations: Single patterns with peak indexing, stacked comparisons with vertical offsets, and reference stick overlays.
2. **Fourier Transform Infrared Spectroscopy (FTIR)**:
   - Transmittance to Absorbance: Beer-Lambert transform $A = -\log_{10}(T)$ with strict nonpositive transmittance rejection (`InvalidUnitsError`).
   - Axis Inversion: Publication-standard descending wavenumber convention ($4000 \to 400\text{ cm}^{-1}$).
   - Functional Group Mapping: Shaded reference regions with explicit disclaimers regarding preliminary assignments.
3. **Raman Spectroscopy**:
   - Signal Conditioning: Cosmic-ray spike filtering and AsLS / SNIP baseline subtraction.
   - Deconvolution: Multi-component Lorentzian fitting of carbon D-band and G-band, and TMD $E_{2g}^1$ / $A_{1g}$ modes.
   - Structural Defect Metrics: Rigorous physical distinction between peak-height ratio ($I_D/I_G$) and integrated peak-area ratio ($A_D/A_G$).
4. **Electron Microscopy (SEM / TEM / HRTEM)**:
   - Physical Calibration: Strict requirement of `pixel_size_nm` (no uncalibrated assumptions).
   - Scale Bars: Standardized rounded lengths ($1, 2, 5 \times 10^n\text{ nm}$) rendered with high-contrast text.
   - Particle Segmentation: Otsu thresholding, morphological opening, connected components, and border clearance.
   - Size Distributions: Equivalent circular diameter ($d_{\mathrm{eq}} = 2\sqrt{A/\pi}$), aspect ratios, $d_{10}, d_{50}, d_{90}$, and probability density histograms with log-normal curves.
   - HRTEM Lattice Resolution: 2D FFT power spectrum, masked reflection indexing, and spatial fringe intensity profiles for atomic $d$-spacing.
5. **Inductively Coupled Plasma (ICP-OES / ICP-MS)**:
   - Calibration Regression: Multi-point linear regression with standard error and correlation coefficient ($R^2 > 0.999$).
   - Analytical Detection Limits: IUPAC Limit of Detection ($\mathrm{LOD} = 3.3\sigma_{\mathrm{blank}}/S$) and Limit of Quantification ($\mathrm{LOQ} = 10\sigma_{\mathrm{blank}}/S$).
   - Catalyst Metal Leaching: Quantifies dissolved mass ($\mu\text{g}$), normalized loss ($\mu\text{g/mg}_{\mathrm{cat}}$), and dissolution percentage ($\%$) during electrochemical durability testing.
6. **Gas Physisorption (BET & BJH)**:
   - Isotherm Processing: Type IV $\mathrm{N}_2$ (77.35 K) adsorption-desorption branch separation with H1 hysteresis loop.
   - Multipoint BET: Specific surface area ($S_{\mathrm{BET}}$) and $C$ constant regression.
   - IUPAC Rouquerol Consistency: Automates verification of all 3 criteria ($q(1-p/p_0)$ monotonic increase, $C > 0$, monolayer pressure inside range).
   - BJH Mesopore Size Distribution: Evaluates $dV/dD$ vs $D_p$ on desorption branch using Halsey thickness and Kelvin capillary condensation.
7. **Atomic Absorption Spectroscopy (AAS)**:
   - Ingestion: Multi-column CSV, TSV, and XLSX exports (PerkinElmer, Shimadzu, Agilent, Thermo Scientific).
   - Flame AAS (FAAS): Steady-state absorbance ($A$, AU), non-linear rollover detection ($A > 0.85$), and EPA Method 7000B compliance.
   - Graphite Furnace AAS (GFAAS): Integrated absorbance ($A \cdot s$) following the Stabilized Temperature Platform Furnace (STPF) concept and EPA Method 7010 / 200.9.
   - Method of Standard Additions (MOSA): Solves matrix suppression in high-salt electrolytes via linear regression extrapolation ($C_{\mathrm{apparent}} = b/m$) with volume dilution scaling and delta-method confidence intervals.
   - Analytical Quality Control: Automated checks for calibration blanks ($< \mathrm{LOD}$), ICV ($\pm 10\text{--}15\%$), CCV ($\pm 10\text{--}15\%$), replicate precision ($\% \mathrm{RSD} \le 10\text{--}15\%$), and Shewhart control charts.
   - Catalyst Metal Loading & Dissolution: Acid digestion bulk loadings (wt.%) and electrochemical stability leaching kinetics with cumulative aliquot withdrawal compensation.
8. **Cross-Technique Comparative Analysis**:
   - Statistical correlation between physical characterization (BET area, crystallite size, Raman defect ratio) and catalytic activity ($\eta_{10}$, $C_{\mathrm{dl}}$, leaching rate).
   - Linear regression, confidence intervals, and mandatory scientific disclaimers separating correlation from mechanistic causality.

### 9.2 Command-Line Interface (CLI)
Generate all 21 characterization publication figures across vector (PDF, SVG) and raster (PNG, TIFF) formats:

```bash
# Generate all materials characterization and AAS figures:
python -m scripts.main_plots --category characterization

# Direct technique plotting:
python -m scripts.main_plots --type xrd --input data/examples/characterization/xrd/catalyst_pristine_xrd.csv
python -m scripts.main_plots --type ftir --input data/examples/characterization/ftir/catalyst_functionalized_ftir.csv
python -m scripts.main_plots --type raman --input data/examples/characterization/raman/carbon_support_raman.csv
python -m scripts.main_plots --type sem --input data/examples/characterization/microscopy/sem_nanoparticles.png --pixel-size 0.50
python -m scripts.main_plots --type hrtem --input data/examples/characterization/microscopy/hrtem_lattice_crop.png --pixel-size 0.018
python -m scripts.main_plots --type icp --input data/examples/characterization/icp/icp_calibration_standards.csv
python -m scripts.main_plots --type aas --input data/examples/characterization/aas/faas_calibration_fe.csv
python -m scripts.main_plots --type aas --input data/examples/characterization/aas/aas_standard_addition_catalyst.csv --calibration-method standard_addition
python -m scripts.main_plots --type bet --input data/examples/characterization/bet/mesoporous_catalyst_n2_isotherm.csv
```

### 9.3 Interactive Jupyter Notebooks
Located in `notebooks/characterization/`:
- [`xrd_analysis.ipynb`](notebooks/characterization/xrd_analysis.ipynb): Powder XRD, Bragg's law, Scherrer domain size, peak deconvolution, and stacked phase transformation.
- [`ftir_raman_analysis.ipynb`](notebooks/characterization/ftir_raman_analysis.ipynb): FTIR Beer-Lambert conversion, functional group mapping, Raman despiking, D/G height and area ratios.
- [`microscopy_analysis.ipynb`](notebooks/characterization/microscopy_analysis.ipynb): FESEM nanoparticle segmentation, size histograms ($d_{50}$), HRTEM 2D FFT, and lattice fringe spacing.
- [`icp_analysis.ipynb`](notebooks/characterization/icp_analysis.ipynb): ICP-OES calibration, IUPAC LOD/LOQ, blank corrections, and catalyst metal dissolution leaching.
- [`bet_analysis.ipynb`](notebooks/characterization/bet_analysis.ipynb): $\mathrm{N}_2$ Type IV isotherms, Rouquerol-compliant BET surface area, BJH mesopore distribution, and cross-technique correlation ($S_{\mathrm{BET}}$ vs $C_{\mathrm{dl}}$).
- [`aas_analysis.ipynb`](notebooks/characterization/aas_analysis.ipynb): FAAS/GFAAS calibration, EPA 7000B QC evaluation, Method of Standard Additions matrix suppression, catalyst loading (wt.%), and stability leaching kinetics.

---

## 10. Running Automated Tests

Run the complete test suite:

```bash
pytest -v
```

All **132 unit and integration tests** pass across 6 domains:
- **Quantum ESPRESSO (DFT & Catalysis)**: SCF, DOS, PDOS, bands, NEB parsing, Fermi alignment, Ry-to-eV conversions, and Sabatier volcano models.
- **Electrochemistry (Metrohm & BioLogic)**: File formats, IUPAC sign conventions, $iR$ compensation, CV redox integration, HER/OER overpotentials, quasi-steady-state Tafel slopes, and Randles EIS fitting.
- **Materials Characterization & Spectroscopy**:
  - Bragg $d$-spacing and Scherrer equation with root-sum-square instrumental broadening subtraction.
  - Beer-Lambert $A = -\log_{10}(T)$ conversion and nonpositive transmittance rejection.
  - Cosmic-ray spike removal, multi-peak Lorentzian deconvolution, and Raman height vs area ratio distinction.
  - Scale bar calibration, Otsu particle segmentation, and HRTEM 2D FFT lattice fringe line profiles.
  - ICP multi-point linear regression, IUPAC LOD/LOQ, and catalyst leaching dissolution metrics.
  - AAS linear calibration, zero-intercept/WLS regression, Method of Standard Additions (MOSA), EPA 7000B QC validation, and aliquot-tracked leaching kinetics.
  - BET linear transform, Rouquerol consistency criteria, and BJH mesopore distributions ($dV/dD$).
  - Cross-technique statistical correlation and causal disclaimers.
- **Chemical Engineering & Visualization**: Arrhenius kinetics, reactor profile parity, journal dimension constraints, TrueType vector font embedding, and colorblind-safe palettes.


"""
Generate and Validate Interactive Jupyter Notebooks for Materials Characterization.

Creates:
1. notebooks/characterization/xrd_analysis.ipynb
2. notebooks/characterization/ftir_raman_analysis.ipynb
3. notebooks/characterization/microscopy_analysis.ipynb
4. notebooks/characterization/icp_analysis.ipynb
5. notebooks/characterization/bet_analysis.ipynb

Executes all 5 notebooks to verify zero runtime errors.
"""

from pathlib import Path
import nbformat as nbf
from nbclient import NotebookClient

PROJECT_ROOT = Path(__file__).resolve().parent.parent
NOTEBOOKS_CHAR_DIR = PROJECT_ROOT / "notebooks" / "characterization"
NOTEBOOKS_CHAR_DIR.mkdir(parents=True, exist_ok=True)


def build_xrd_notebook() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell(
            "# Quantitative X-Ray Diffraction (XRD) Analysis\n\n"
            "This notebook demonstrates quantitative analysis of powder X-ray diffraction patterns:\n"
            "- Ingestion of raw XRD data (CSV and comment-delimited text)\n"
            "- Bragg's law $d$-spacing calculation: $\\lambda = 2d\\sin\\theta$\n"
            "- Scherrer crystallite domain size calculation with instrumental broadening correction:\n"
            "  $$\\beta_{\\mathrm{corr}} = \\sqrt{\\beta_{\\mathrm{obs}}^2 - \\beta_{\\mathrm{inst}}^2}, \\quad D = \\frac{K \\lambda}{\\beta_{\\mathrm{corr}} \\cos\\theta}$$\n"
            "- Multi-peak pseudo-Voigt deconvolution using `lmfit` with residual analysis\n"
            "- Stacked comparative diffractograms with calibrated vertical offsets."
        ),
        nbf.v4.new_code_cell(
            "import sys\n"
            "from pathlib import Path\n"
            "# Ensure project root is in sys.path\n"
            "root_path = Path.cwd().parent.parent if Path.cwd().name == 'characterization' else Path.cwd()\n"
            "if str(root_path) not in sys.path:\n"
            "    sys.path.insert(0, str(root_path))\n\n"
            "import matplotlib.pyplot as plt\n"
            "import numpy as np\n"
            "from parsers.characterization import parse_xrd_file\n"
            "from modules.characterization.xrd import (\n"
            "    calculate_d_spacing,\n"
            "    calculate_crystallite_size_scherrer,\n"
            "    analyze_xrd_peaks,\n"
            "    plot_xrd_pattern,\n"
            "    plot_xrd_stacked,\n"
            "    plot_xrd_fit,\n"
            ")\n"
            "from modules.characterization.common_peaks import fit_peaks\n"
            "from scripts.utils import apply_publication_style\n\n"
            "apply_publication_style()\n"
            "print('XRD analysis modules loaded successfully.')"
        ),
        nbf.v4.new_markdown_cell(
            "## 1. Single Diffractogram & Peak Analysis\n"
            "Load pristine Pt/C catalyst diffractogram collected with Cu K$\\alpha$ radiation ($\\lambda = 1.54060\\text{ \\AA}$)."
        ),
        nbf.v4.new_code_cell(
            "data_path = root_path / 'data' / 'examples' / 'characterization' / 'xrd' / 'catalyst_pristine_xrd.csv'\n"
            "xrd_data = parse_xrd_file(data_path)\n"
            "print(f'Loaded: {xrd_data.sample_id} ({xrd_data.n_points} points, 2theta: {xrd_data.x.min():.1f} to {xrd_data.x.max():.1f} deg)')\n\n"
            "peaks = analyze_xrd_peaks(xrd_data, prominence=50.0)\n"
            "for p in peaks:\n"
            "    print(f'Peak {p[\"peak_index\"]}: 2theta = {p[\"two_theta_deg\"]:.2f} deg, d = {p[\"d_spacing_angstrom\"]:.3f} A, Size = {p[\"crystallite_size_nm\"]:.1f} nm')\n\n"
            "fig, ax = plot_xrd_pattern(xrd_data, annotate_peaks=True, show_d_spacing=True, title='Pt/C Catalyst Powder XRD')\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 2. Stacked Diffractograms (Pristine vs Treated Catalyst)\n"
            "Compare crystal phase transformations upon thermal treatment."
        ),
        nbf.v4.new_code_cell(
            "treated_path = root_path / 'data' / 'examples' / 'characterization' / 'xrd' / 'catalyst_treated_xrd.txt'\n"
            "treated_data = parse_xrd_file(treated_path)\n\n"
            "fig, ax = plot_xrd_stacked(\n"
            "    [xrd_data, treated_data],\n"
            "    labels=['Pristine Pt/C', 'Thermally Treated Pt/C'],\n"
            "    vertical_offset_fraction=0.35,\n"
            "    title='Phase Evolution: Pristine vs Thermally Activated'\n"
            ")\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 3. High-Resolution Peak Deconvolution & Scherrer Size\n"
            "Deconvolute the Pt(111) reflection around $2\\theta = 39.8^\\circ$ using a pseudo-Voigt profile."
        ),
        nbf.v4.new_code_cell(
            "mask = (xrd_data.x >= 36.0) & (xrd_data.x <= 44.0)\n"
            "sub_x, sub_y = xrd_data.x[mask], xrd_data.y[mask]\n"
            "fit_res = fit_peaks(sub_x, sub_y, peak_centers=[39.76], model_type='pseudo_voigt', background='linear')\n\n"
            "ctr = fit_res['peaks'][0]['center']\n"
            "fwhm = fit_res['peaks'][0]['fwhm']\n"
            "scherrer = calculate_crystallite_size_scherrer(ctr, fwhm, instrumental_fwhm_deg=0.08)\n"
            "print(f'Pt(111) Reflection: Center = {ctr:.2f} deg, FWHM = {fwhm:.3f} deg')\n"
            "print(f'Scherrer Domain Size: {scherrer[\"crystallite_size_nm\"]:.2f} nm ({scherrer[\"crystallite_size_angstrom\"]:.1f} A)')\n\n"
            "sub_data = xrd_data.copy()\n"
            "sub_data.x = sub_x\n"
            "sub_data.y = sub_y\n"
            "fig, axes = plot_xrd_fit(sub_data, fit_res, title=f'Pt(111) Reflection Deconvolution (Domain Size = {scherrer[\"crystallite_size_nm\"]:.1f} nm)')\n"
            "plt.show()"
        ),
    ]
    return nb


def build_ftir_raman_notebook() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell(
            "# Vibrational Spectroscopy: FTIR & Raman Analysis\n\n"
            "This notebook demonstrates quantitative vibrational spectroscopy for catalyst characterization:\n"
            "- FTIR transmittance to absorbance conversion via Beer-Lambert: $A = -\\log_{10}(T)$\n"
            "- Inverted descending wavenumber axis ($4000 \\to 400\\text{ cm}^{-1}$)\n"
            "- Characteristic functional group region mapping\n"
            "- Raman cosmic-ray spike rejection & AsLS baseline subtraction\n"
            "- Raman D-band and G-band deconvolution with explicit distinction between:\n"
            "  - Peak-height ratio ($I_D / I_G$)\n"
            "  - Integrated peak-area ratio ($A_D / A_G$)."
        ),
        nbf.v4.new_code_cell(
            "import sys\n"
            "from pathlib import Path\n"
            "root_path = Path.cwd().parent.parent if Path.cwd().name == 'characterization' else Path.cwd()\n"
            "if str(root_path) not in sys.path:\n"
            "    sys.path.insert(0, str(root_path))\n\n"
            "import matplotlib.pyplot as plt\n"
            "import numpy as np\n"
            "from parsers.characterization import parse_ftir_file, convert_transmittance_to_absorbance, parse_raman_file\n"
            "from modules.characterization.ftir import plot_ftir_spectrum, plot_ftir_stacked, analyze_ftir_bands\n"
            "from modules.characterization.raman import (\n"
            "    process_raman_spectrum,\n"
            "    calculate_peak_ratios,\n"
            "    plot_raman_spectrum,\n"
            "    plot_raman_stacked,\n"
            "    plot_raman_deconvolution,\n"
            ")\n"
            "from modules.characterization.common_peaks import fit_peaks\n"
            "from scripts.utils import apply_publication_style\n\n"
            "apply_publication_style()\n"
            "print('Spectroscopy modules loaded.')"
        ),
        nbf.v4.new_markdown_cell(
            "## 1. FTIR Surface Functionalization & Band Identification\n"
            "Load catalyst FTIR transmittance and convert to absorbance."
        ),
        nbf.v4.new_code_cell(
            "ftir_path = root_path / 'data' / 'examples' / 'characterization' / 'ftir' / 'catalyst_functionalized_ftir.csv'\n"
            "raw_ftir = parse_ftir_file(ftir_path)\n"
            "abs_ftir = convert_transmittance_to_absorbance(raw_ftir)\n\n"
            "fig, ax = plot_ftir_spectrum(abs_ftir, annotate_peaks=True, show_regions=True, title='Functionalized Catalyst Surface FTIR')\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 2. FTIR Degradation Comparison: Pristine vs Spent\n"
            "Examine loss of surface ligands following electrochemical cycling."
        ),
        nbf.v4.new_code_cell(
            "spent_path = root_path / 'data' / 'examples' / 'characterization' / 'ftir' / 'catalyst_spent_ftir.csv'\n"
            "spent_ftir = convert_transmittance_to_absorbance(parse_ftir_file(spent_path))\n\n"
            "fig, ax = plot_ftir_stacked([abs_ftir, spent_ftir], labels=['Pristine', 'Spent Post-OER'], vertical_offset_fraction=0.4, title='FTIR Surface Degradation')\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 3. Raman Spectroscopy: Defect Density & D/G Ratios\n"
            "Deconvolute carbon black support D-band (~1345 cm^-1) and G-band (~1585 cm^-1)."
        ),
        nbf.v4.new_code_cell(
            "raman_path = root_path / 'data' / 'examples' / 'characterization' / 'raman' / 'carbon_support_raman.csv'\n"
            "raw_raman = parse_raman_file(raman_path)\n"
            "proc_raman = process_raman_spectrum(raw_raman, despike=True, subtract_bg=True, bg_method='asls')\n\n"
            "fit_res = fit_peaks(proc_raman.x, proc_raman.y, peak_centers=[1345.0, 1585.0], model_type='lorentzian', background='linear')\n"
            "ratios = calculate_peak_ratios(fit_res, peak_1_index=1, peak_2_index=2, label_1='D', label_2='G')\n"
            "print(f'Peak Height Ratio (I_D / I_G): {ratios[\"height_ratio\"]:.3f}')\n"
            "print(f'Integrated Area Ratio (A_D / A_G): {ratios[\"area_ratio\"]:.3f}')\n\n"
            "fig, axes = plot_raman_deconvolution(proc_raman, fit_res, title='Carbon Support D-Band & G-Band Deconvolution')\n"
            "plt.show()"
        ),
    ]
    return nb


def build_microscopy_notebook() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell(
            "# Electron Microscopy Image Analysis (SEM / TEM / HRTEM)\n\n"
            "This notebook demonstrates quantitative microscopy analysis:\n"
            "- Mathematical physical scale bar rendering with true physical units\n"
            "- Otsu segmentation and morphological filtering of nanoparticles\n"
            "- Equivalent circular diameter statistics ($d_{10}, d_{50}, d_{90}$)\n"
            "- Log-normal particle size distribution fitting\n"
            "- HRTEM 2D FFT power spectrum and atomic lattice fringe $d$-spacing profiling."
        ),
        nbf.v4.new_code_cell(
            "import sys\n"
            "from pathlib import Path\n"
            "root_path = Path.cwd().parent.parent if Path.cwd().name == 'characterization' else Path.cwd()\n"
            "if str(root_path) not in sys.path:\n"
            "    sys.path.insert(0, str(root_path))\n\n"
            "import matplotlib.pyplot as plt\n"
            "import numpy as np\n"
            "from parsers.characterization import parse_microscopy_image\n"
            "from modules.characterization.microscopy import (\n"
            "    segment_particles,\n"
            "    analyze_hrtem_lattice_fringes,\n"
            "    plot_microscopy_with_scalebar,\n"
            "    plot_segmentation_overlay,\n"
            "    plot_particle_size_distribution,\n"
            "    plot_hrtem_fft_panel,\n"
            ")\n"
            "from scripts.utils import apply_publication_style\n\n"
            "apply_publication_style()\n"
            "print('Microscopy analysis modules loaded.')"
        ),
        nbf.v4.new_markdown_cell(
            "## 1. SEM Nanoparticle Segmentation & Size Distribution\n"
            "Load FESEM catalyst image with physical calibration of 0.50 nm/pixel."
        ),
        nbf.v4.new_code_cell(
            "sem_path = root_path / 'data' / 'examples' / 'characterization' / 'microscopy' / 'sem_nanoparticles.png'\n"
            "sem_data = parse_microscopy_image(sem_path, pixel_size_nm=0.50)\n"
            "seg_res = segment_particles(sem_data, min_area_px=15, clear_borders=True)\n"
            "stats = seg_res['statistics']\n"
            "print(f'Segmented {stats[\"particle_count\"]} particles')\n"
            "print(f'Mean Diameter: {stats[\"mean_diameter_nm\"]:.1f} +/- {stats[\"std_diameter_nm\"]:.1f} nm')\n"
            "print(f'Median Diameter (d50): {stats[\"median_diameter_nm\"]:.1f} nm (d10={stats[\"p10_diameter_nm\"]:.1f}, d90={stats[\"p90_diameter_nm\"]:.1f} nm)')\n\n"
            "fig1, ax1 = plot_segmentation_overlay(sem_data, seg_res, title='FESEM Catalyst Nanoparticles & Segmentation')\n"
            "plt.show()\n\n"
            "fig2, ax2 = plot_particle_size_distribution(seg_res, bins=15, title='Nanoparticle Size Distribution')\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 2. HRTEM 2D FFT Lattice Fringe Analysis\n"
            "Load high-resolution TEM crop calibrated at 0.018 nm/pixel (0.18 A/pixel)."
        ),
        nbf.v4.new_code_cell(
            "hrtem_path = root_path / 'data' / 'examples' / 'characterization' / 'microscopy' / 'hrtem_lattice_crop.png'\n"
            "hrtem_data = parse_microscopy_image(hrtem_path, pixel_size_nm=0.018)\n"
            "fringe_res = analyze_hrtem_lattice_fringes(hrtem_data, fringe_direction='horizontal')\n"
            "print(f'Lattice fringe d-spacing: {fringe_res[\"d_spacing_angstrom\"]:.2f} A ({fringe_res[\"d_spacing_nm\"]:.3f} nm)')\n\n"
            "fig, axes = plot_hrtem_fft_panel(hrtem_data, fringe_res, title=f'HRTEM Atomic Lattice Resolution (d = {fringe_res[\"d_spacing_angstrom\"]:.2f} A)')\n"
            "plt.show()"
        ),
    ]
    return nb


def build_icp_notebook() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell(
            "# ICP-OES / ICP-MS Elemental Analysis & Catalyst Leaching\n\n"
            "This notebook demonstrates quantitative analytical chemistry processing:\n"
            "- Multi-point standard calibration with weighted linear regression\n"
            "- IUPAC Limit of Detection ($\\mathrm{LOD} = 3.3\\,\\sigma_{\\mathrm{blank}} / S$) and Limit of Quantification ($\\mathrm{LOQ} = 10\\,\\sigma_{\\mathrm{blank}} / S$)\n"
            "- Solvent blank correction and sample dilution adjustment\n"
            "- Quantitative catalyst metal leaching dissolution loss (% and $\\mu\\text{g/mg}_{\\text{cat}}$)."
        ),
        nbf.v4.new_code_cell(
            "import sys\n"
            "from pathlib import Path\n"
            "root_path = Path.cwd().parent.parent if Path.cwd().name == 'characterization' else Path.cwd()\n"
            "if str(root_path) not in sys.path:\n"
            "    sys.path.insert(0, str(root_path))\n\n"
            "import matplotlib.pyplot as plt\n"
            "import numpy as np\n"
            "import pandas as pd\n"
            "from parsers.characterization import parse_icp_file\n"
            "from modules.characterization.icp import (\n"
            "    fit_icp_calibration,\n"
            "    calculate_catalyst_leaching,\n"
            "    plot_icp_calibration,\n"
            "    plot_icp_leaching_comparison,\n"
            ")\n"
            "from scripts.utils import apply_publication_style\n\n"
            "apply_publication_style()\n"
            "print('ICP modules loaded.')"
        ),
        nbf.v4.new_markdown_cell(
            "## 1. Multi-Point ICP Calibration & LOD Determination\n"
            "Evaluate calibration standards for Pt, Ni, Fe."
        ),
        nbf.v4.new_code_cell(
            "cal_path = root_path / 'data' / 'examples' / 'characterization' / 'icp' / 'icp_calibration_standards.csv'\n"
            "cal_data = parse_icp_file(cal_path)\n"
            "df_cal = cal_data.dataframe\n\n"
            "for elem in ['Pt', 'Ni', 'Fe']:\n"
            "    sub = df_cal[df_cal['element'] == elem]\n"
            "    stds = sub[sub['sample_id'].str.startswith('STD')]\n"
            "    blanks = sub[sub['sample_id'].str.startswith('BLANK')]\n"
            "    res = fit_icp_calibration(stds['concentration'].values, stds['intensity'].values, blank_intensities=blanks['intensity'].values)\n"
            "    print(f'{elem}: Slope = {res[\"slope\"]:.1f} cts/(mg/L), R^2 = {res[\"r_squared\"]:.5f}, LOD = {res[\"lod_conc\"]:.3f} mg/L, LOQ = {res[\"loq_conc\"]:.3f} mg/L')\n\n"
            "# Plot Ni calibration curve\n"
            "sub_ni = df_cal[df_cal['element'] == 'Ni']\n"
            "stds_ni = sub_ni[sub_ni['sample_id'].str.startswith('STD')]\n"
            "blanks_ni = sub_ni[sub_ni['sample_id'].str.startswith('BLANK')]\n"
            "cal_ni = fit_icp_calibration(stds_ni['concentration'].values, stds_ni['intensity'].values, blank_intensities=blanks_ni['intensity'].values)\n"
            "fig, axes = plot_icp_calibration(cal_ni, element='Nickel')\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 2. Catalyst Leaching Durability Study\n"
            "Assess Ni and Fe metal dissolution from 2.0 mg catalyst in 50 mL 1.0 M KOH electrolyte during 100 h OER."
        ),
        nbf.v4.new_code_cell(
            "leach_path = root_path / 'data' / 'examples' / 'characterization' / 'icp' / 'catalyst_leaching_study.csv'\n"
            "leach_data = parse_icp_file(leach_path)\n"
            "df_l = leach_data.dataframe\n"
            "df_l_ni = df_l[df_l['element'] == 'Ni']\n\n"
            "pct_list = []\n"
            "samples = df_l_ni['sample_id'].tolist()\n"
            "for _, row in df_l_ni.iterrows():\n"
            "    l_res = calculate_catalyst_leaching(\n"
            "        measured_concentration_mg_L=float(row['concentration']),\n"
            "        electrolyte_volume_mL=50.0,\n"
            "        catalyst_mass_mg=2.0,\n"
            "        initial_metal_wt_fraction=0.20,\n"
            "        dilution_factor=float(row.get('dilution_factor', 1.0)),\n"
            "    )\n"
            "    pct_list.append(l_res['percentage_metal_leached'])\n"
            "    print(f'{row[\"sample_id\"]}: {l_res[\"percentage_metal_leached\"]:.2f}% leached ({l_res[\"mass_normalized_leaching_ug_per_mg\"]:.1f} ug/mg_cat)')\n\n"
            "fig, ax = plot_icp_leaching_comparison(samples, pct_list, element='Nickel', title='Catalytic Ni Leaching During OER')\n"
            "plt.show()"
        ),
    ]
    return nb


def build_bet_notebook() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell(
            "# Gas Physisorption: BET Surface Area & BJH Pore Distribution\n\n"
            "This notebook demonstrates quantitative gas sorption analysis:\n"
            "- $\\text{N}_2$ 77.35 K adsorption-desorption isotherm parsing & branch separation\n"
            "- Multi-point linear BET regression transform:\n"
            "  $$\\frac{p/p_0}{q(1 - p/p_0)} = \\frac{1}{q_m C} + \\frac{C-1}{q_m C} \\left(\\frac{p}{p_0}\\right)$$\n"
            "- Rigorous IUPAC Rouquerol consistency criteria verification\n"
            "- BJH mesopore size distribution ($dV/dD$) via Halsey thickness and Kelvin capillary condensation equations\n"
            "- Multi-technique cross-correlation: BET specific surface area vs double-layer capacitance ($C_{\\text{dl}}$)."
        ),
        nbf.v4.new_code_cell(
            "import sys\n"
            "from pathlib import Path\n"
            "root_path = Path.cwd().parent.parent if Path.cwd().name == 'characterization' else Path.cwd()\n"
            "if str(root_path) not in sys.path:\n"
            "    sys.path.insert(0, str(root_path))\n\n"
            "import matplotlib.pyplot as plt\n"
            "import numpy as np\n"
            "from parsers.characterization import parse_isotherm_file\n"
            "from modules.characterization.bet import (\n"
            "    calculate_bet_surface_area,\n"
            "    calculate_bjh_pore_distribution,\n"
            "    plot_isotherm,\n"
            "    plot_bet_linear,\n"
            "    plot_bjh_pore_distribution,\n"
            ")\n"
            "from modules.characterization.comparative import (\n"
            "    correlate_cross_technique_metrics,\n"
            "    plot_cross_technique_correlation,\n"
            ")\n"
            "from scripts.utils import apply_publication_style\n\n"
            "apply_publication_style()\n"
            "print('Physisorption modules loaded.')"
        ),
        nbf.v4.new_markdown_cell(
            "## 1. Type IV Isotherm & BET Surface Area with Rouquerol Checks\n"
            "Inspect $\\text{N}_2$ isotherm with H1 hysteresis loop."
        ),
        nbf.v4.new_code_cell(
            "iso_path = root_path / 'data' / 'examples' / 'characterization' / 'bet' / 'mesoporous_catalyst_n2_isotherm.csv'\n"
            "iso_data = parse_isotherm_file(iso_path)\n"
            "df_ads = iso_data.dataframe[iso_data.dataframe['branch'] == 'adsorption']\n"
            "df_des = iso_data.dataframe[iso_data.dataframe['branch'] == 'desorption']\n\n"
            "fig, ax = plot_isotherm(iso_data, title='Type IV N2 Adsorption-Desorption Isotherm (77.3 K)')\n"
            "plt.show()\n\n"
            "bet_res = calculate_bet_surface_area(df_ads['p_p0'].values, df_ads['quantity_adsorbed'].values, p_range=(0.05, 0.30))\n"
            "print(f'S_BET: {bet_res[\"s_bet_m2_g\"]:.1f} m^2/g, C constant = {bet_res[\"c_constant\"]:.1f}, R^2 = {bet_res[\"r_squared\"]:.5f}')\n"
            "rouq = bet_res['rouquerol_criteria']\n"
            "print(f'Rouquerol Check: Monotonic={rouq[\"monotonic_increasing\"]}, C>0={rouq[\"c_positive\"]}, Overall Valid={rouq[\"overall_valid\"]}')\n\n"
            "fig2, ax2 = plot_bet_linear(bet_res, title='BET Linear Transform')\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 2. BJH Mesopore Size Distribution (dV/dD)\n"
            "Calculate incremental and differential pore size distribution on desorption branch."
        ),
        nbf.v4.new_code_cell(
            "bjh_res = calculate_bjh_pore_distribution(df_des['p_p0'].values, df_des['quantity_adsorbed'].values, branch='desorption')\n"
            "print(f'Total Mesopore Volume: {bjh_res[\"total_pore_volume_cm3_g\"]:.3f} cm^3/g')\n"
            "print(f'Mean Pore Diameter: {bjh_res[\"mean_pore_diameter_nm\"]:.1f} nm')\n\n"
            "fig, ax = plot_bjh_pore_distribution(bjh_res, title='BJH Mesopore Size Distribution (Desorption Branch)')\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 3. Multi-Technique Correlation: BET Surface Area vs Double-Layer Capacitance\n"
            "Correlate gas physisorption physical area with electrochemical active surface metrics across 5 catalyst synthesis batches."
        ),
        nbf.v4.new_code_cell(
            "bet_areas = [45.0, 92.0, 148.0, 210.0, 280.0]\n"
            "c_dl = [12.5, 23.8, 38.2, 53.0, 71.4]\n"
            "corr = correlate_cross_technique_metrics(\n"
            "    x_values=bet_areas,\n"
            "    y_values=c_dl,\n"
            "    x_label='BET Specific Surface Area $S_{\\\\mathrm{BET}}~[\\\\mathrm{m^2\\\\cdot g^{-1}}]$',\n"
            "    y_label='Double-Layer Capacitance $C_{\\\\mathrm{dl}}~[\\\\mathrm{mF\\\\cdot cm^{-2}}]$',\n"
            "    sample_labels=['Cat-1', 'Cat-2', 'Cat-3', 'Cat-4', 'Cat-5'],\n"
            ")\n"
            "print(f'Correlation R^2: {corr[\"r_squared\"]:.4f}, p-value: {corr[\"p_value\"]:.3e}')\n"
            "print(f'Disclaimer: {corr[\"disclaimer\"]}')\n\n"
            "fig, ax = plot_cross_technique_correlation(corr, title='Cross-Technique Validation: BET Area vs Double-Layer Capacitance')\n"
            "plt.show()"
        ),
    ]
    return nb


def build_aas_notebook() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell(
            "# Atomic Absorption Spectroscopy (AAS) Quantitative Analysis & Quality Control\n\n"
            "This interactive notebook provides end-to-end data processing for Atomic Absorption Spectroscopy (FAAS & GFAAS):\n"
            "1. **Raw Data Ingestion & Metadata Inspection** (CSV, TSV, XLSX)\n"
            "2. **Calibration Curve Modeling** (Linear OLS, WLS, Zero-Intercept, Curvature Warnings)\n"
            "3. **Analytical Quality Control (QC)** Compliance Verification (US EPA Method 7000B & 7010)\n"
            "4. **Quantitative Sample Concentration & Below-LOD/LOQ Flagging**\n"
            "5. **Method of Standard Additions (MOSA)** for Matrix-Suppressed Samples with Delta-Method Uncertainty\n"
            "6. **Catalyst Metal Loading Quantification (wt.%)** from Acid Digestions\n"
            "7. **Electrochemical Dissolution & Leaching Kinetics** with Cumulative Aliquot Tracking\n"
            "8. **Journal-Quality Multi-Format Figure Exports** with Cryptographic Provenance Tracking."
        ),
        nbf.v4.new_code_cell(
            "import sys\n"
            "from pathlib import Path\n"
            "root_path = Path.cwd().parent.parent if Path.cwd().name == 'characterization' else Path.cwd()\n"
            "if str(root_path) not in sys.path:\n"
            "    sys.path.insert(0, str(root_path))\n\n"
            "import matplotlib.pyplot as plt\n"
            "import numpy as np\n"
            "import pandas as pd\n"
            "from parsers.characterization.aas import parse_aas_file\n"
            "from modules.characterization.aas import (\n"
            "    fit_aas_calibration,\n"
            "    fit_standard_addition,\n"
            "    quantify_aas_samples,\n"
            "    evaluate_aas_qc,\n"
            "    calculate_aas_catalyst_loading,\n"
            "    calculate_aas_leaching,\n"
            "    plot_aas_calibration,\n"
            "    plot_standard_addition,\n"
            "    plot_aas_concentration_comparison,\n"
            "    plot_aas_catalyst_loading,\n"
            "    plot_aas_leaching,\n"
            "    plot_aas_stability_before_after,\n"
            "    plot_aas_replicate_distribution,\n"
            "    plot_aas_qc_chart,\n"
            ")\n"
            "from parsers.characterization.common import record_characterization_provenance\n"
            "from scripts.utils import apply_publication_style\n\n"
            "apply_publication_style()\n"
            "print('AAS processing and visualization modules loaded successfully.')"
        ),
        nbf.v4.new_markdown_cell(
            "## 1. Raw Data Ingestion & Metadata Inspection\n\n"
            "Load Flame AAS (FAAS) Fe calibration standards and inspect instrument metadata."
        ),
        nbf.v4.new_code_cell(
            "base_data = root_path / 'data' / 'examples' / 'characterization' / 'aas'\n"
            "faas_data = parse_aas_file(base_data / 'faas_calibration_fe.csv', technique='FAAS')\n"
            "print('Technique:', faas_data.technique)\n"
            "print('Instrument Model:', faas_data.metadata.get('instrument_model'))\n"
            "print('Background Correction:', faas_data.metadata.get('background_correction'))\n"
            "print('Elements Detected:', faas_data.metadata.get('elements_detected'))\n"
            "faas_df = faas_data.dataframe\n"
            "faas_df.head(10)"
        ),
        nbf.v4.new_markdown_cell(
            "## 2. Calibration Curve Modeling & Detection Limits\n\n"
            "Perform linear regression ($A = mC + b$), derive sensitivity slope, $R^2$, and IUPAC detection limit ($3.3 s_b / m$)."
        ),
        nbf.v4.new_code_cell(
            "stds_df = faas_df[faas_df['sample_type'] == 'Standard']\n"
            "blanks_df = faas_df[faas_df['sample_type'] == 'Blank']\n"
            "cal_fe = fit_aas_calibration(\n"
            "    standards_concentration=stds_df['concentration_reported'].values,\n"
            "    standards_signal=stds_df['absorbance'].values,\n"
            "    signal_type='absorbance',\n"
            "    blank_signals=blanks_df['absorbance'].values,\n"
            "    lod_method='blank_stdev',\n"
            ")\n"
            "print(f\"Slope: {cal_fe['slope']:.4f} AU/(mg/L)\")\n"
            "print(f\"Intercept: {cal_fe['intercept']:.4f} AU\")\n"
            "print(f\"R^2: {cal_fe['r_squared']:.5f}\")\n"
            "print(f\"LOD: {cal_fe['lod_conc']:.4f} mg/L, LOQ: {cal_fe['loq_conc']:.4f} mg/L\")\n\n"
            "fig_cal, (ax_top, ax_bot) = plot_aas_calibration(cal_fe, element='Iron (Fe 248.3 nm)', unit='mg/L')\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 3. Analytical Quality Control (QC) Compliance (US EPA 7000B)\n\n"
            "Verify calibration linearity, blank cleanliness, ICV recovery, and replicate precision against EPA Method 7000B."
        ),
        nbf.v4.new_code_cell(
            "qc_report = evaluate_aas_qc(faas_df, cal_fe, sop='EPA_7000B')\n"
            "print('Overall QC Status:', qc_report['overall_status'])\n"
            "display(qc_report['checks'])\n\n"
            "qc_run_df = pd.read_csv(base_data / 'aas_qc_run_sequence.csv')\n"
            "fig_qc, ax_qc = plot_aas_qc_chart(qc_run_df['QC_Check'], qc_run_df['Recovery_pct'], sop='EPA_7000B')\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 4. Method of Standard Additions (MOSA) for Matrix Suppression\n\n"
            "When analyzing complex high-salt electrolyte samples, matrix interference suppresses analyte sensitivity.\n"
            "Standard addition fits $S = m \\Delta C + S_0$ and extrapolates backward to find sample concentration."
        ),
        nbf.v4.new_code_cell(
            "mosa_data = parse_aas_file(base_data / 'aas_standard_addition_catalyst.csv')\n"
            "mosa_df = mosa_data.dataframe\n"
            "spikes = [0.0, 0.50, 1.00, 1.50, 2.00]\n"
            "sigs = [mosa_df[mosa_df['sample_id'].str.contains(f'\\\\+{sp:.2f}') if sp > 0 else mosa_df['sample_id'].str.contains('Unspiked')]['absorbance'].mean() for sp in spikes]\n\n"
            "sa_res = fit_standard_addition(spikes, sigs, sample_volume_mL=5.0, final_volume_mL=25.0)\n"
            "print(f\"Apparent Flask Conc: {sa_res['apparent_concentration']:.3f} +/- {sa_res['ci95_apparent']:.3f} mg/L\")\n"
            "print(f\"True Original Sample Conc: {sa_res['original_concentration']:.3f} +/- {sa_res['ci95_original']:.3f} mg/L (DF = {sa_res['dilution_factor']:.1f})\")\n\n"
            "fig_sa, ax_sa = plot_standard_addition(sa_res, element='Cobalt (Co in 0.1 M KOH)', unit='mg/L')\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 5. Catalyst Composition & Metal Loading (wt.%)\n\n"
            "Quantify metal loading from acid-digested solid catalyst samples and compare with nominal targets."
        ),
        nbf.v4.new_code_cell(
            "dig_df = pd.read_csv(base_data / 'catalyst_digestion_loading.csv')\n"
            "samples = dig_df[dig_df['Sample_Type'] == 'Sample']\n"
            "pt_samples = samples[samples['Element'] == 'Pt']\n"
            "loadings = [calculate_aas_catalyst_loading(r['Measured_Conc_mg_L'], r['Digest_Vol_mL'], r['Mass_mg'], r['Dilution_Factor'])['metal_loading_wt_pct'] for _, r in pt_samples.iterrows()]\n"
            "print(f\"Pt/C Actual Loading: {np.mean(loadings[:3]):.2f} +/- {np.std(loadings[:3], ddof=1):.2f} wt% (Nominal: 20.0 wt%)\")\n\n"
            "fig_load, ax_load = plot_aas_catalyst_loading(\n"
            "    ['Pt/C', 'PtCo/C (Pt)', 'PtCo/C (Co)', 'NiFe (Ni)', 'NiFe (Fe)'],\n"
            "    [19.8, 14.8, 5.1, 24.5, 12.2],\n"
            "    nominal_loading_wt_pct=[20.0, 15.0, 5.0, 25.0, 12.0],\n"
            "    error_bars=[0.2, 0.2, 0.1, 0.3, 0.2],\n"
            ")\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 6. Electrocatalytic Dissolution & Leaching Kinetics\n\n"
            "Track metal dissolution during a 100-hour alkaline OER stability test with aliquot sampling corrections."
        ),
        nbf.v4.new_code_cell(
            "leach_df = pd.read_csv(base_data / 'catalyst_stability_leaching.csv')\n"
            "co_df = leach_df[leach_df['Element'] == 'Co']\n"
            "time_pts = co_df['Time_hours'].tolist()\n"
            "pcts = []\n"
            "losses = []\n"
            "aliquots = []\n"
            "for _, r in co_df.iterrows():\n"
            "    res = calculate_aas_leaching(\n"
            "        r['Measured_Conc_mg_L'], r['Electrolyte_Volume_mL'], r['Initial_Catalyst_Mass_mg'],\n"
            "        r['Initial_Metal_wt_pct'] / 100.0, aliquot_history=aliquots if aliquots else None\n"
            "    )\n"
            "    pcts.append(res['percentage_metal_leached'])\n"
            "    losses.append(res['mass_normalized_leaching_ug_per_mg'])\n"
            "    if r['Aliquot_Withdrawn_mL'] > 0:\n"
            "        aliquots.append((r['Aliquot_Withdrawn_mL'], r['Measured_Conc_mg_L']))\n\n"
            "print(f\"Total Co Leached after 100h: {pcts[-1]:.2f}%\")\n"
            "print(f\"Mass-normalized metal loss: {losses[-1]:.2f} ug/mg_cat\")\n\n"
            "fig_leach, ax_leach = plot_aas_leaching([f'{t}h' for t in time_pts], pcts, mass_normalized_loss_ug_per_mg=losses, element='Cobalt')\n"
            "plt.show()"
        ),
    ]
    return nb


def generate_and_execute_all() -> None:
    notebooks = {
        "xrd_analysis.ipynb": build_xrd_notebook(),
        "ftir_raman_analysis.ipynb": build_ftir_raman_notebook(),
        "microscopy_analysis.ipynb": build_microscopy_notebook(),
        "icp_analysis.ipynb": build_icp_notebook(),
        "bet_analysis.ipynb": build_bet_notebook(),
        "aas_analysis.ipynb": build_aas_notebook(),
    }

    for fname, nb in notebooks.items():
        out_path = NOTEBOOKS_CHAR_DIR / fname
        print(f"Writing {fname} to {out_path}...")
        with open(out_path, "w", encoding="utf-8") as f:
            nbf.write(nb, f)

        print(f"Executing {fname} to verify zero runtime errors...")
        client = NotebookClient(nb, timeout=60, kernel_name="python3", resources={"metadata": {"path": str(NOTEBOOKS_CHAR_DIR)}})
        client.execute()
        # Save executed notebook with outputs
        with open(out_path, "w", encoding="utf-8") as f:
            nbf.write(nb, f)
        print(f"Successfully executed and validated {fname}.")

    print("\nAll 6 characterization notebooks generated and executed with zero errors!")


if __name__ == "__main__":
    generate_and_execute_all()

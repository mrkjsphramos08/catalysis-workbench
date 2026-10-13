"""
Script to generate the three specialized electrochemistry Jupyter notebooks with verified API calls.
"""

from pathlib import Path
import nbformat as nbf

NOTEBOOKS_DIR = Path(__file__).resolve().parent.parent / "notebooks"
NOTEBOOKS_DIR.mkdir(exist_ok=True)


def create_metrohm_notebook():
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell("""# Metrohm Autolab (NOVA) Experimental Electrochemistry Workflow

This interactive research notebook provides an end-to-end scientific analysis pipeline for experimental electrochemistry datasets acquired from **Metrohm Autolab** potentiostats using **NOVA** software.

### Analysis Modules Covered:
1. **Oxygen Evolution Reaction (OER) LSV**: Calibration of reference potentials ($E$ vs RHE), software uncompensated resistance ($iR$) compensation, and overpotential extraction ($\eta_{10}, \eta_{50}, \eta_{100}$).
2. **Multi-Cycle Cyclic Voltammetry (CV)**: Peak detection ($E_{pa}, E_{pc}, \Delta E_p$) and voltammetric charge integration ($Q = \\int I dt$).
3. **Electrochemical Impedance Spectroscopy (FRA EIS)**: Complex impedance polar conversion, publication Nyquist plots, and Randles circuit fitting ($R_s - [R_{ct} \\parallel \\mathrm{CPE}]$).
4. **Potentiostatic Chronoamperometry (CA)**: Long-term durability analysis, current retention percentage, and degradation rate (\\%/h).
"""),
        nbf.v4.new_code_cell("""import sys
from pathlib import Path

# Add project root to sys.path
project_root = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

import numpy as np
import matplotlib.pyplot as plt

from scripts.utils import apply_publication_style, create_publication_figure
from parsers.electrochemistry import load_electrochemistry_data
from modules.electrochemistry import (
    plot_lsv,
    extract_overpotential_at_current,
    plot_cv,
    extract_redox_peaks,
    calculate_voltammetric_charge,
    plot_nyquist,
    fit_equivalent_circuit,
    calculate_impedance_polar,
    plot_chronoamperometry,
    calculate_ca_stability_metrics,
)

apply_publication_style()
data_dir = project_root / "data" / "examples" / "electrochemistry"
print("Environment successfully initialized with publication styles.")
"""),
        nbf.v4.new_markdown_cell("""## 1. Oxygen Evolution Reaction (OER) Linear Sweep Voltammetry

Load raw Metrohm Autolab ASCII export (`metrohm_oer_lsv.txt`), inspect automatically parsed metadata, and extract the overpotential delivering $10\\,\\mathrm{mA/cm^2}$ geometric current density ($\eta_{10}$).
"""),
        nbf.v4.new_code_cell("""lsv_file = data_dir / "metrohm" / "metrohm_oer_lsv.txt"
data_oer = load_electrochemistry_data(lsv_file, instrument="metrohm", technique="LSV")

print(f"Instrument / Software: {data_oer.instrument_software}")
print(f"Working Electrode Area: {data_oer.electrode_area_cm2} cm^2")
print(f"Electrolyte pH: {data_oer.ph}")
print(f"Reference Electrode: {data_oer.potential_reference_raw}")
print(f"Uncompensated Resistance (Ru): {data_oer.ru_ohm} Ohm")

pot_rhe = data_oer.get_potential_vs_rhe()
j = data_oer.get_current_density(unit="mA/cm^2")

eta_10 = extract_overpotential_at_current(pot_rhe, j, target_j=10.0, reaction="OER")
eta_20 = extract_overpotential_at_current(pot_rhe, j, target_j=20.0, reaction="OER")
print(f"\\n[Metrics] Overpotential eta_10: {eta_10:.1f} mV")
print(f"[Metrics] Overpotential eta_20: {eta_20:.1f} mV" if eta_20 is not None else "[Metrics] Overpotential eta_20: N/A")

ax, bm = plot_lsv(
    data_oer,
    reaction="OER",
    ir_compensate=False,
    target_currents=(10.0, 20.0),
    title="Metrohm Autolab - OER Polarization (1.0 M KOH)",
)
plt.show()
"""),
        nbf.v4.new_markdown_cell("""## 2. Multi-Cycle Cyclic Voltammetry (CV)

Analyze a 3-cycle CV scan (`metrohm_cv_multicycle.csv`), identify redox peak potentials ($E_{pa}, E_{pc}$), peak separation ($\Delta E_p$), and calculate voltammetric charge passed during cycling.
"""),
        nbf.v4.new_code_cell("""cv_file = data_dir / "metrohm" / "metrohm_cv_multicycle.csv"
data_cv = load_electrochemistry_data(cv_file, instrument="metrohm", technique="CV")

cycle_3 = data_cv.get_cycle(3)
peaks = extract_redox_peaks(cycle_3.potential_V, cycle_3.current_A)
q_res = calculate_voltammetric_charge(cycle_3.time_s, cycle_3.current_A)

print(f"Anodic Peak Potential (E_pa): {peaks['anodic_peak_potential_V']:.3f} V | Current: {peaks['anodic_peak_current_A']:.4e} A")
print(f"Cathodic Peak Potential (E_pc): {peaks['cathodic_peak_potential_V']:.3f} V | Current: {peaks['cathodic_peak_current_A']:.4e} A")
print(f"Peak Potential Separation (Delta E_p): {peaks['delta_ep_mV']:.1f} mV")
print(f"Anodic Charge: {q_res['anodic_charge_C']:.4f} C | Cathodic Charge: {q_res['cathodic_charge_C']:.4f} C")

fig, ax = plot_cv(
    data_cv,
    cycles=[1, 2, 3],
    title="Metrohm Autolab - Multi-Cycle CV",
)
plt.show()
"""),
        nbf.v4.new_markdown_cell("""## 3. FRA Electrochemical Impedance Spectroscopy (EIS)

Fit high-frequency to low-frequency AC impedance data to a Randles equivalent circuit ($R_s - [R_{ct} \\parallel \\mathrm{CPE}]$) and display an orthonormal Nyquist representation.
"""),
        nbf.v4.new_code_cell("""eis_file = data_dir / "metrohm" / "metrohm_fra_eis.txt"
data_eis = load_electrochemistry_data(eis_file, instrument="metrohm", technique="EIS")

f = data_eis.frequency_Hz
zr = data_eis.impedance_real_ohm
zi = data_eis.impedance_imag_ohm

fit_res = fit_equivalent_circuit(f, zr, zi)

print(f"Solution Resistance (R_s): {fit_res['rs_ohm']:.2f} Ohm")
print(f"Charge-Transfer Resistance (R_ct): {fit_res['rct_ohm']:.2f} Ohm")
print(f"Fit RMSE: {fit_res['rmse']:.3f} Ohm")

fig, ax = plot_nyquist(
    data_eis,
    fit_results=fit_res,
    title="Metrohm Autolab - FRA Nyquist Spectrum & Randles Fit",
)
plt.show()
"""),
        nbf.v4.new_markdown_cell("""## 4. Potentiostatic Chronoamperometry (CA) Stability

Evaluate catalyst durability under constant bias ($E = 1.60\\,\\mathrm{V}$ vs RHE). Compute current retention percentage and degradation rate (\\%/h).
"""),
        nbf.v4.new_code_cell("""ca_file = data_dir / "metrohm" / "metrohm_stability_ca.txt"
data_ca = load_electrochemistry_data(ca_file, instrument="metrohm", technique="CA")

t = data_ca.time_s
j = data_ca.get_current_density(unit="mA/cm^2")

metrics_ca = calculate_ca_stability_metrics(t, j)

print(f"Duration: {metrics_ca['duration_hours']:.2f} hours")
print(f"Current Retention: {metrics_ca['retention_percent']:.1f} %")
print(f"Degradation Rate: {metrics_ca['degradation_rate_pct_per_h']:.2f} %/h")

fig, ax = plot_chronoamperometry(
    data_ca,
    potential_label="E = 1.60 V vs RHE",
    show_retention=True,
    title="Metrohm Autolab - 1-Hour OER Durability",
)
plt.show()
""")
    ]

    out_path = NOTEBOOKS_DIR / "metrohm_analysis.ipynb"
    with open(out_path, "w", encoding="utf-8") as f:
        nbf.write(nb, f)
    print(f"Created: {out_path}")


def create_biologic_notebook():
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell("""# BioLogic EC-Lab Experimental Electrochemistry Workflow

This interactive research notebook provides an end-to-end scientific processing pipeline for experimental electrochemistry data acquired from **BioLogic** potentiostats using **EC-Lab** software.

### Analysis Modules Covered:
1. **Hydrogen Evolution Reaction (HER) LSV**: Automatic `.mpt` ASCII parsing, RHE conversion, and cathodic overpotential extraction ($\eta_{-10}, \eta_{-50}$).
2. **Potentio Electrochemical Impedance Spectroscopy (PEIS)**: Dual-panel stacked Bode plots (magnitude $|Z|$ and phase angle $\\phi$ vs frequency).
3. **Lin-KK Kramers-Kronig Consistency Validation**: Rigorous statistical verification of stationarity, linearity, and causality using Schönleber's criterion.
4. **Galvanostatic Chronopotentiometry (CP)**: Durability monitoring at $-10\\,\\mathrm{mA/cm^2}$ and potential drift rate estimation ($\\mathrm{mV/h}$).
"""),
        nbf.v4.new_code_cell("""import sys
from pathlib import Path

# Add project root to sys.path
project_root = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

import numpy as np
import matplotlib.pyplot as plt

from scripts.utils import apply_publication_style, create_publication_figure
from parsers.electrochemistry import load_electrochemistry_data
from modules.electrochemistry import (
    plot_lsv,
    extract_overpotential_at_current,
    plot_bode,
    run_kramers_kronig_test,
    plot_chronopotentiometry,
    calculate_cp_stability_metrics,
)

apply_publication_style()
data_dir = project_root / "data" / "examples" / "electrochemistry"
print("BioLogic processing pipeline loaded.")
"""),
        nbf.v4.new_markdown_cell("""## 1. Hydrogen Evolution Reaction (HER) Polarization Curves

Load raw BioLogic `.mpt` export (`biologic_her_lsv.mpt`) and extract cathodic overpotentials ($\eta_{-10}$ and $\eta_{-50}$).
"""),
        nbf.v4.new_code_cell("""her_file = data_dir / "biologic" / "biologic_her_lsv.mpt"
data_her = load_electrochemistry_data(her_file, instrument="biologic", technique="LSV")

print(f"Software: {data_her.instrument_software}")
print(f"Electrode Area: {data_her.electrode_area_cm2} cm^2")
print(f"Reference: {data_her.potential_reference_raw} | pH: {data_her.ph}")

pot_her_rhe = data_her.get_potential_vs_rhe()
j_her = data_her.get_current_density(unit="mA/cm^2")

eta_10_her = extract_overpotential_at_current(pot_her_rhe, j_her, target_j=10.0, reaction="HER")
eta_20_her = extract_overpotential_at_current(pot_her_rhe, j_her, target_j=20.0, reaction="HER")

print(f"\\n[Metrics] Overpotential eta_-10: {eta_10_her:.1f} mV")
print(f"[Metrics] Overpotential eta_-20: {eta_20_her:.1f} mV" if eta_20_her is not None else "[Metrics] Overpotential eta_-20: N/A")

ax, bm = plot_lsv(
    data_her,
    reaction="HER",
    ir_compensate=False,
    target_currents=(10.0, 20.0),
    title="BioLogic EC-Lab - HER Polarization (0.5 M H2SO4)",
)
plt.show()
"""),
        nbf.v4.new_markdown_cell("""## 2. PEIS Stacked Bode Analysis & Lin-KK Validation

Load BioLogic `.mpt` impedance spectrum (`biologic_peis.mpt`), visualize Bode frequency responses, and validate system stability via Lin-KK Kramers-Kronig transform.
"""),
        nbf.v4.new_code_cell("""peis_file = data_dir / "biologic" / "biologic_peis.mpt"
data_peis = load_electrochemistry_data(peis_file, instrument="biologic", technique="EIS")

f = data_peis.frequency_Hz
zr = data_peis.impedance_real_ohm
zi = data_peis.impedance_imag_ohm

kk_res = run_kramers_kronig_test(f, zr, zi)

print(f"Lin-KK Valid: {kk_res['valid']}")
print(f"Lin-KK mu score: {kk_res['mu']:.4f}")
print(f"Number of RC Elements: {kk_res['num_elements']}")

fig, (ax_mag, ax_phase) = plot_bode(
    data_peis,
    title="BioLogic EC-Lab - PEIS Bode Representation",
)
plt.show()
"""),
        nbf.v4.new_markdown_cell("""## 3. Galvanostatic Chronopotentiometry (CP) Stability

Measure potential stability during continuous HER operation at $j = -10\\,\\mathrm{mA/cm^2}$ over 1 hour. Extract potential drift rate in $\\mathrm{mV/h}$.
"""),
        nbf.v4.new_code_cell("""cp_file = data_dir / "biologic" / "biologic_stability_cp.mpt"
data_cp = load_electrochemistry_data(cp_file, instrument="biologic", technique="CP")

t = data_cp.time_s
v = data_cp.get_potential_vs_rhe(reference="ag/agcl", ph=0.3)

metrics_cp = calculate_cp_stability_metrics(t, v)

print(f"Potential Drift: {metrics_cp['potential_drift_mV']:.2f} mV")
print(f"Drift Rate: {metrics_cp['drift_rate_mV_per_h']:.2f} mV/h")

fig, ax = plot_chronopotentiometry(
    data_cp,
    ph=0.3,
    current_label="j = -10 mA/cm^2",
    show_drift_rate=True,
    title="BioLogic EC-Lab - 1-Hour HER Galvanostatic Durability",
)
plt.show()
""")
    ]

    out_path = NOTEBOOKS_DIR / "biologic_analysis.ipynb"
    with open(out_path, "w", encoding="utf-8") as f:
        nbf.write(nb, f)
    print(f"Created: {out_path}")


def create_comparison_notebook():
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell("""# Cross-Instrument Electrochemical Benchmarking
## Metrohm Autolab (NOVA) vs. BioLogic (EC-Lab)

This notebook executes a rigorous inter-instrument comparison between **Metrohm Autolab** and **BioLogic** potentiostats for the identical catalyst system under standardized experimental conditions (1.0 M KOH, scan rate $\\nu = 5\\,\\mathrm{mV/s}$, $R_u = 4.2\\,\\Omega$, $A = 0.196\\,\\mathrm{cm}^2$).

### Objectives:
1. Demonstrate standardized data harmonization across disparate manufacturer formats.
2. Quantify kinetic agreement ($\eta_{10}, \eta_{50}$) and inter-instrument residual current densities $\Delta j(E)$.
3. Compare quasi-steady-state Tafel slopes and assess potential measurement bias.
"""),
        nbf.v4.new_code_cell("""import sys
from pathlib import Path

# Add project root to sys.path
project_root = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scripts.utils import apply_publication_style, create_publication_figure
from parsers.electrochemistry import load_electrochemistry_data
from modules.electrochemistry import (
    compare_instrument_lsv,
    plot_cross_instrument_benchmark,
    calculate_tafel_slope,
    extract_overpotential_at_current,
)

apply_publication_style()
data_dir = project_root / "data" / "examples" / "electrochemistry"
print("Cross-instrument benchmarking environment ready.")
"""),
        nbf.v4.new_markdown_cell("""## 1. Import and Align Instrument Datasets

Load matched OER polarization scans from both workstations.
"""),
        nbf.v4.new_code_cell("""metrohm_file = data_dir / "metrohm" / "metrohm_oer_lsv.txt"
biologic_file = data_dir / "biologic" / "biologic_oer_lsv.mpt"

data_m = load_electrochemistry_data(metrohm_file, instrument="metrohm", technique="LSV")
data_b = load_electrochemistry_data(biologic_file, instrument="biologic", technique="LSV")

print(f"Metrohm Data Points: {len(data_m.potential_V)} | Area: {data_m.electrode_area_cm2} cm^2")
print(f"BioLogic Data Points: {len(data_b.potential_V)} | Area: {data_b.electrode_area_cm2} cm^2")
"""),
        nbf.v4.new_markdown_cell("""## 2. Cross-Instrument Polarization Benchmark & Residuals

Evaluate the residual difference $\Delta j(E) = j_{\\mathrm{Metrohm}}(E) - j_{\\mathrm{BioLogic}}(E)$ across the kinetic polarization window.
"""),
        nbf.v4.new_code_cell("""comp = compare_instrument_lsv(data_m, data_b, reaction="OER", target_j=10.0)

print("--- Inter-Instrument Benchmark Statistics ---")
print(f"Conditions Match: {comp['conditions_match']}")
print(f"Metrohm eta_10:  {comp['eta10_metrohm_mV']:.1f} mV")
print(f"BioLogic eta_10: {comp['eta10_biologic_mV']:.1f} mV")
print(f"Delta eta_10:    {comp['delta_eta10_mV']:.2f} mV")
print(f"Current Density RMSE: {comp['rmse_current_density']:.3f} mA/cm^2")

fig, (ax_top, ax_bot) = plot_cross_instrument_benchmark(
    data_m,
    data_b,
    reaction="OER",
    title="Cross-Instrument OER Benchmark: Metrohm Autolab vs BioLogic",
)
plt.show()
"""),
        nbf.v4.new_markdown_cell("""## 3. Kinetic Parameters Concordance Table

Compare extracted overpotentials ($\eta_{10}, \eta_{50}$) and Tafel slopes between instruments.
"""),
        nbf.v4.new_code_cell("""# Extract overpotentials
pot_m = data_m.get_potential_vs_rhe()
j_m = data_m.get_current_density(unit="mA/cm^2")
pot_b = data_b.get_potential_vs_rhe()
j_b = data_b.get_current_density(unit="mA/cm^2")

eta_m_10 = extract_overpotential_at_current(pot_m, j_m, 10.0, "OER")
eta_b_10 = extract_overpotential_at_current(pot_b, j_b, 10.0, "OER")
eta_m_20 = extract_overpotential_at_current(pot_m, j_m, 20.0, "OER")
eta_b_20 = extract_overpotential_at_current(pot_b, j_b, 20.0, "OER")

# Tafel Slopes in kinetic region (1.0 to 10.0 mA/cm^2)
tafel_m = calculate_tafel_slope(pot_m, j_m, current_range=(1.0, 10.0), reaction="OER", scan_rate_V_s=0.005)
tafel_b = calculate_tafel_slope(pot_b, j_b, current_range=(1.0, 10.0), reaction="OER", scan_rate_V_s=0.005)

summary_df = pd.DataFrame([
    {
        "Parameter": "eta_10 (mV)",
        "Metrohm Autolab": f"{eta_m_10:.1f}",
        "BioLogic EC-Lab": f"{eta_b_10:.1f}",
        "Absolute Delta": f"{abs(eta_m_10 - eta_b_10):.1f} mV"
    },
    {
        "Parameter": "eta_20 (mV)",
        "Metrohm Autolab": f"{eta_m_20:.1f}",
        "BioLogic EC-Lab": f"{eta_b_20:.1f}",
        "Absolute Delta": f"{abs(eta_m_20 - eta_b_20):.1f} mV"
    },
    {
        "Parameter": "Tafel Slope (mV/dec)",
        "Metrohm Autolab": f"{tafel_m['tafel_slope_mV_dec']:.1f}",
        "BioLogic EC-Lab": f"{tafel_b['tafel_slope_mV_dec']:.1f}",
        "Absolute Delta": f"{abs(tafel_m['tafel_slope_mV_dec'] - tafel_b['tafel_slope_mV_dec']):.1f} mV/dec"
    },
    {
        "Parameter": "Exchange Current j_0 (mA/cm^2)",
        "Metrohm Autolab": f"{tafel_m['exchange_current_density_mA_cm2']:.2e}",
        "BioLogic EC-Lab": f"{tafel_b['exchange_current_density_mA_cm2']:.2e}",
        "Absolute Delta": "-"
    }
])

print(summary_df.to_string(index=False))
""")
    ]

    out_path = NOTEBOOKS_DIR / "electrochemical_comparison.ipynb"
    with open(out_path, "w", encoding="utf-8") as f:
        nbf.write(nb, f)
    print(f"Created: {out_path}")


if __name__ == "__main__":
    create_metrohm_notebook()
    create_biologic_notebook()
    create_comparison_notebook()

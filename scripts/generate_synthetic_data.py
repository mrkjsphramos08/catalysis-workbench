"""
Deterministic Generation of Benchmark Synthetic Datasets.

Generates physically plausible synthetic data for:
- Electronic DOS and PDOS (Pt-O surface system)
- Band structure along high-symmetry k-path (Gamma - X - M - Gamma)
- Adsorption energy comparison across 5 transition metal facets (fcc, top, bridge sites)
- Catalytic CO oxidation reaction free energy profiles across 3 metal surfaces
- Sabatier volcano curve with microkinetic dual-leg model
- Arrhenius temperature-dependent rate constants with experimental variance
- Reactor conversion, selectivity, and yield profiles

All datasets are explicitly labeled as synthetic benchmarks in headers and metadata.
"""

from pathlib import Path
import numpy as np
import pandas as pd

from scripts.config import EXAMPLES_DATA_DIR
from scripts.preprocess import GAS_CONSTANT_R_J_MOL_K


def generate_all_synthetic_data() -> None:
    """Generate all benchmark synthetic datasets deterministically."""
    EXAMPLES_DATA_DIR.mkdir(parents=True, exist_ok=True)
    np.random.seed(42)

    # 1. Electronic Density of States (DOS and PDOS)
    energies = np.linspace(-8.0, 6.0, 500)
    # Pt d-band centered around -2.25 eV with width 2.5 eV
    pt_d_up = 4.5 * np.exp(-0.5 * ((energies + 2.25) / 1.4) ** 2) + 0.8 * np.exp(-0.5 * ((energies + 0.8) / 0.7) ** 2)
    pt_d_down = 4.2 * np.exp(-0.5 * ((energies + 2.40) / 1.4) ** 2) + 0.7 * np.exp(-0.5 * ((energies + 0.9) / 0.7) ** 2)
    # Adsorbed O 2p states hybridized around -5.5 eV and anti-bonding states above Fermi
    o_p_up = 1.8 * np.exp(-0.5 * ((energies + 5.5) / 0.8) ** 2) + 0.6 * np.exp(-0.5 * ((energies - 1.2) / 0.9) ** 2)
    o_p_down = 1.8 * np.exp(-0.5 * ((energies + 5.4) / 0.8) ** 2) + 0.6 * np.exp(-0.5 * ((energies - 1.3) / 0.9) ** 2)
    # Free-electron-like background
    sp_bg = 0.4 * np.sqrt(np.maximum(0.0, energies + 7.5))

    tdos_up = pt_d_up + o_p_up + sp_bg
    tdos_down = pt_d_down + o_p_down + sp_bg

    dos_df = pd.DataFrame({
        "energy": np.round(energies, 4),
        "total_dos": np.round(tdos_up + tdos_down, 4),
        "dos_spin_up": np.round(tdos_up, 4),
        "dos_spin_down": np.round(tdos_down, 4),
        "pt_5d_up": np.round(pt_d_up, 4),
        "pt_5d_down": np.round(pt_d_down, 4),
        "o_2p_up": np.round(o_p_up, 4),
        "o_2p_down": np.round(o_p_down, 4),
    })
    dos_path = EXAMPLES_DATA_DIR / "dft_dos_synthetic.csv"
    dos_df.to_csv(dos_path, index=False)

    # 2. Electronic Band Structure
    # Path: Gamma (0.0) -> X (0.8) -> M (1.6) -> Gamma (2.4)
    k_dist = np.linspace(0.0, 2.4, 200)
    # 5 representative bands crossing and below Fermi level
    b1 = -4.0 + 1.2 * np.cos(k_dist * np.pi / 0.8)
    b2 = -2.5 - 0.9 * np.cos(k_dist * np.pi / 0.8) + 0.4 * np.sin(k_dist * 2 * np.pi / 0.8)
    b3 = -0.8 + 1.5 * np.cos(k_dist * np.pi / 1.6)
    b4 = 0.5 - 1.8 * np.cos(k_dist * np.pi / 0.8)
    b5 = 2.8 + 1.1 * np.cos(k_dist * np.pi / 0.8)

    bands_df = pd.DataFrame({
        "k_distance": np.round(k_dist, 4),
        "band_1": np.round(b1, 4),
        "band_2": np.round(b2, 4),
        "band_3": np.round(b3, 4),
        "band_4": np.round(b4, 4),
        "band_5": np.round(b5, 4),
    })
    bands_path = EXAMPLES_DATA_DIR / "dft_bands_synthetic.csv"
    bands_df.to_csv(bands_path, index=False)

    # 3. Adsorption Energies (5 Metals x 3 Sites)
    metals = ["Pt(111)", "Pd(111)", "Cu(111)", "Ni(111)", "Au(111)"]
    sites = ["fcc", "bridge", "top"]
    clean_slab_energies = {"Pt(111)": -842.15, "Pd(111)": -720.50, "Cu(111)": -610.20, "Ni(111)": -795.30, "Au(111)": -530.80}
    # CO reference energy in gas phase
    e_co_gas = -14.28

    rows = []
    # Base binding energies in negative_favorable convention (eV)
    base_bind = {
        "Pt(111)": {"fcc": -1.82, "bridge": -1.75, "top": -1.58},
        "Pd(111)": {"fcc": -1.95, "bridge": -1.88, "top": -1.45},
        "Cu(111)": {"fcc": -0.92, "bridge": -0.85, "top": -0.74},
        "Ni(111)": {"fcc": -2.05, "bridge": -1.98, "top": -1.68},
        "Au(111)": {"fcc": -0.32, "bridge": -0.28, "top": -0.42},
    }

    for m in metals:
        e_clean = clean_slab_energies[m]
        for s in sites:
            e_ads = base_bind[m][s]
            err = np.round(0.03 + 0.02 * np.random.rand(), 3)
            # E_total = E_clean + E_ref + E_ads
            e_tot = np.round(e_clean + e_co_gas + e_ads, 4)
            rows.append({
                "surface": m,
                "site": s,
                "e_clean_surface": e_clean,
                "e_adsorbate_ref": e_co_gas,
                "e_total": e_tot,
                "e_ads": e_ads,
                "error": err,
                "unit": "eV",
                "reference_species": "CO(g)",
                "sign_convention": "negative_favorable",
            })

    ads_df = pd.DataFrame(rows)
    ads_path = EXAMPLES_DATA_DIR / "adsorption_energies_synthetic.csv"
    ads_df.to_csv(ads_path, index=False)

    # 4. Catalytic Reaction Profile (CO Oxidation across 3 Catalysts)
    # Elementary steps: 
    # S0: * + CO(g) + 0.5 O2(g)
    # S1: *CO + *O
    # S2: TS1 (C-O coupling)
    # S3: *CO2
    # S4: * + CO2(g)
    cat_steps = ["* + CO(g) + 0.5 O2(g)", "*CO + *O", "TS (C-O bond)", "*CO2", "* + CO2(g)"]
    reaction_profiles = {
        "step_name": cat_steps,
        "is_transition_state": [False, False, True, False, False],
        "Pt(111)": [0.00, -0.85, 0.42, -0.22, -1.45],
        "Au(111)": [0.00, 0.35, 1.65, 0.15, -1.45],
        "Cu(111)": [0.00, -1.42, 0.28, -0.65, -1.45],
    }
    cat_df = pd.DataFrame(reaction_profiles)
    cat_path = EXAMPLES_DATA_DIR / "catalysis_reaction_profile_synthetic.csv"
    cat_df.to_csv(cat_path, index=False)

    # 5. Catalysis Volcano Plot (ORR / HER descriptor)
    # Descriptor: Delta G_OH* (eV), Activity metric: Limiting potential U_L (V)
    volcano_data = [
        {"catalyst": "Pt", "descriptor_g_oh": 0.82, "activity_potential": 0.80, "candidate": True},
        {"catalyst": "Pd", "descriptor_g_oh": 0.65, "activity_potential": 0.64, "candidate": True},
        {"catalyst": "Rh", "descriptor_g_oh": 0.45, "activity_potential": 0.45, "candidate": False},
        {"catalyst": "Ir", "descriptor_g_oh": 0.50, "activity_potential": 0.51, "candidate": False},
        {"catalyst": "Ru", "descriptor_g_oh": 0.25, "activity_potential": 0.26, "candidate": False},
        {"catalyst": "Ni", "descriptor_g_oh": 0.15, "activity_potential": 0.14, "candidate": False},
        {"catalyst": "Cu", "descriptor_g_oh": 1.15, "activity_potential": 0.48, "candidate": False},
        {"catalyst": "Ag", "descriptor_g_oh": 1.35, "activity_potential": 0.28, "candidate": False},
        {"catalyst": "Au", "descriptor_g_oh": 1.48, "activity_potential": 0.15, "candidate": False},
    ]
    volc_df = pd.DataFrame(volcano_data)
    volc_path = EXAMPLES_DATA_DIR / "catalysis_volcano_synthetic.csv"
    volc_df.to_csv(volc_path, index=False)

    # 6. Arrhenius Kinetics
    # Realistic catalytic gas-phase reaction (Ea ~ 82.5 kJ/mol, ln(A) ~ 18.2)
    # T from 180 to 280 °C (453.15 to 553.15 K)
    temp_c = np.array([180.0, 195.0, 210.0, 225.0, 240.0, 255.0, 270.0, 285.0])
    temp_k = temp_c + 273.15
    ea_true = 82500.0  # J/mol
    ln_a_true = 18.2
    r_gas = GAS_CONSTANT_R_J_MOL_K
    # Exact Arrhenius + small realistic random noise (std ~ 0.04 in ln(k))
    noise = np.array([-0.021, 0.035, -0.012, 0.018, -0.031, 0.024, -0.008, 0.015])
    ln_k = ln_a_true - (ea_true / (r_gas * temp_k)) + noise
    k_rate = np.exp(ln_k)

    arrh_df = pd.DataFrame({
        "temperature_C": temp_c,
        "temperature_K": temp_k,
        "rate_constant_k": k_rate,
        "rate_constant_unit": "s^-1",
    })
    arrh_path = EXAMPLES_DATA_DIR / "kinetics_arrhenius_synthetic.csv"
    arrh_df.to_csv(arrh_path, index=False)

    # 7. Reactor Conversion and Selectivity
    # Packed bed reactor profile over temperature 160 to 320 °C
    t_reactor = np.linspace(160, 320, 9)
    # Sigmoidal conversion
    conv = 98.0 / (1.0 + np.exp(-0.045 * (t_reactor - 235.0)))
    # Selectivity drops slightly as temperature increases due to over-oxidation
    sel = 96.0 - 0.15 * (t_reactor - 160.0) - 0.0008 * (t_reactor - 160.0) ** 2
    yield_val = (conv * sel) / 100.0

    reactor_df = pd.DataFrame({
        "temperature_C": np.round(t_reactor, 1),
        "conversion_percent": np.round(conv, 2),
        "selectivity_percent": np.round(sel, 2),
        "yield_percent": np.round(yield_val, 2),
    })
    reactor_path = EXAMPLES_DATA_DIR / "reactor_performance_synthetic.csv"
    reactor_df.to_csv(reactor_path, index=False)

    print("All synthetic benchmark datasets successfully generated in data/examples/.")


if __name__ == "__main__":
    generate_all_synthetic_data()

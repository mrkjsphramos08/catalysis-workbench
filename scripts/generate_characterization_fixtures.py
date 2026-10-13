"""
Deterministic Generator for Scientific Materials Characterization Benchmark Datasets.

Generates realistic, physically grounded synthetic datasets for:
1. XRD: Pristine and thermally treated catalytic patterns with Cu K-alpha
2. FTIR: Transmittance spectra with functional group bands (OH, CH, C=O, MO)
3. Raman: Carbon support with D and G bands (including a cosmic spike) and TMD catalyst
4. Microscopy: SEM dispersed nanoparticles and HRTEM periodic atomic lattice fringes
5. ICP-OES/MS: Multi-element calibration curves and catalyst dissolution leaching series
6. BET/BJH: IUPAC Type IV N2 physisorption isotherm with H1 hysteresis loop.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image

PROJECT_ROOT = Path(__file__).resolve().parent.parent
EXAMPLES_DIR = PROJECT_ROOT / "data" / "examples" / "characterization"
FIXTURES_DIR = PROJECT_ROOT / "tests" / "characterization" / "fixtures"


def generate_all_fixtures():
    np.random.seed(42)

    # Ensure target directories exist
    for base in [EXAMPLES_DIR, FIXTURES_DIR]:
        for sub in ["xrd", "ftir", "raman", "microscopy", "icp", "bet"]:
            (base / sub).mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # 1. XRD Patterns (Cu K-alpha: 1.54060 A)
    # =========================================================================
    print("Generating XRD benchmark datasets...")
    two_theta = np.linspace(15.0, 85.0, 1401)  # step = 0.05 deg

    # Amorphous halo from carbon support (centered at 24 deg)
    halo = 350.0 * np.exp(-0.5 * ((two_theta - 24.0) / 4.5) ** 2)
    bg = 120.0 + 0.5 * two_theta + halo

    # Pt fcc reflections: (111) at 39.76 deg, (200) at 46.24 deg, (220) at 67.45 deg, (311) at 81.28 deg
    def pseudo_voigt(x, x0, fwhm, height, eta=0.5):
        sigma = fwhm / 2.35482
        gamma = fwhm / 2.0
        g = np.exp(-0.5 * ((x - x0) / sigma) ** 2)
        l = 1.0 / (1.0 + ((x - x0) / gamma) ** 2)
        return height * (eta * l + (1.0 - eta) * g)

    # Pristine Pt/C (small crystallites ~4.5 nm, FWHM ~ 1.8 deg)
    pt_111 = pseudo_voigt(two_theta, 39.76, 1.85, 1800.0)
    pt_200 = pseudo_voigt(two_theta, 46.24, 1.95, 950.0)
    pt_220 = pseudo_voigt(two_theta, 67.45, 2.10, 680.0)
    pt_311 = pseudo_voigt(two_theta, 81.28, 2.30, 520.0)
    noise_1 = np.random.normal(0, 18.0, len(two_theta))

    int_pristine = np.maximum(0, bg + pt_111 + pt_200 + pt_220 + pt_311 + noise_1)

    df_xrd_p = pd.DataFrame({"2Theta": two_theta, "Intensity": np.round(int_pristine, 1)})
    for d in [EXAMPLES_DIR / "xrd", FIXTURES_DIR / "xrd"]:
        df_xrd_p.to_csv(d / "catalyst_pristine_xrd.csv", index=False)

    # Thermally treated Pt/C (crystallite growth ~14 nm, FWHM ~ 0.65 deg, slight lattice compression +0.12 deg)
    pt_111_t = pseudo_voigt(two_theta, 39.88, 0.65, 4200.0)
    pt_200_t = pseudo_voigt(two_theta, 46.36, 0.70, 2200.0)
    pt_220_t = pseudo_voigt(two_theta, 67.57, 0.78, 1500.0)
    pt_311_t = pseudo_voigt(two_theta, 81.40, 0.85, 1150.0)
    noise_2 = np.random.normal(0, 22.0, len(two_theta))

    int_treated = np.maximum(0, bg + pt_111_t + pt_200_t + pt_220_t + pt_311_t + noise_2)
    # Save as space-delimited .txt with metadata comments
    content_txt = "# Instrument: Bruker D8 Advance\n# Radiation: Cu K-alpha (1.54060 A)\n# Voltage: 40 kV, Current: 40 mA\n2Theta\tIntensity\n"
    for th, val in zip(two_theta, int_treated):
        content_txt += f"{th:.2f}\t{val:.1f}\n"

    for d in [EXAMPLES_DIR / "xrd", FIXTURES_DIR / "xrd"]:
        with open(d / "catalyst_treated_xrd.txt", "w", encoding="utf-8") as f:
            f.write(content_txt)

    # =========================================================================
    # 2. FTIR Spectra (4000 to 400 cm^-1)
    # =========================================================================
    print("Generating FTIR benchmark datasets...")
    wn = np.linspace(400.0, 4000.0, 1801)  # step = 2 cm^-1

    # Base baseline ~ 96% Transmittance
    base_t = 96.0 - 5.0 * ((wn - 400.0) / 3600.0)

    # Absorption bands as dips in transmittance (%T)
    # 1. Broad O-H stretch centered at 3420 cm^-1
    b_oh = 32.0 * np.exp(-0.5 * ((wn - 3420.0) / 160.0) ** 2)
    # 2. C-H aliphatic stretch doublet at 2924 and 2853 cm^-1
    b_ch1 = 18.0 * np.exp(-0.5 * ((wn - 2924.0) / 25.0) ** 2)
    b_ch2 = 14.0 * np.exp(-0.5 * ((wn - 2853.0) / 22.0) ** 2)
    # 3. C=O carbonyl at 1715 cm^-1
    b_co = 45.0 * np.exp(-0.5 * ((wn - 1715.0) / 32.0) ** 2)
    # 4. C=C aromatic at 1620 cm^-1
    b_cc = 24.0 * np.exp(-0.5 * ((wn - 1620.0) / 30.0) ** 2)
    # 5. M-O lattice vibration at 580 cm^-1
    b_mo = 38.0 * np.exp(-0.5 * ((wn - 580.0) / 45.0) ** 2)

    noise_ftir = np.random.normal(0, 0.4, len(wn))
    t_pristine = np.clip(base_t - (b_oh + b_ch1 + b_ch2 + b_co + b_cc + b_mo) + noise_ftir, 5.0, 100.0)

    df_ftir_p = pd.DataFrame({"Wavenumber": wn, "Transmittance": np.round(t_pristine, 2)})
    for d in [EXAMPLES_DIR / "ftir", FIXTURES_DIR / "ftir"]:
        df_ftir_p.to_csv(d / "catalyst_functionalized_ftir.csv", index=False)

    # Spent catalyst (decreased C=O, enhanced M-O and surface carbonate ~1420 cm^-1)
    b_carb = 30.0 * np.exp(-0.5 * ((wn - 1420.0) / 40.0) ** 2)
    t_spent = np.clip(base_t - (b_oh * 1.2 + b_ch1 * 0.4 + b_ch2 * 0.4 + b_co * 0.3 + b_cc + b_mo * 1.3 + b_carb) + noise_ftir, 5.0, 100.0)

    df_ftir_s = pd.DataFrame({"Wavenumber": wn, "Transmittance": np.round(t_spent, 2)})
    for d in [EXAMPLES_DIR / "ftir", FIXTURES_DIR / "ftir"]:
        df_ftir_s.to_csv(d / "catalyst_spent_ftir.csv", index=False)

    # =========================================================================
    # 3. Raman Spectra (800 to 2000 cm^-1)
    # =========================================================================
    print("Generating Raman benchmark datasets...")
    shift = np.linspace(800.0, 2000.0, 1201)  # step = 1 cm^-1

    # Fluorescent background (parabolic curvature)
    fl_bg = 200.0 + 0.15 * (shift - 800.0) + 0.0003 * ((shift - 1400.0) ** 2)

    # Carbon D-band (disorder) at 1345 cm^-1, G-band (graphitic) at 1585 cm^-1
    d_band = 850.0 / (1.0 + ((shift - 1345.0) / 48.0) ** 2)
    g_band = 1050.0 / (1.0 + ((shift - 1585.0) / 36.0) ** 2)
    noise_raman = np.random.normal(0, 12.0, len(shift))

    raman_int = fl_bg + d_band + g_band + noise_raman

    # Inject a realistic sharp cosmic-ray spike at index 450 (1250 cm^-1)
    raman_int[450] += 3200.0

    df_raman_c = pd.DataFrame({"Raman_Shift": shift, "Intensity": np.round(raman_int, 1)})
    for d in [EXAMPLES_DIR / "raman", FIXTURES_DIR / "raman"]:
        df_raman_c.to_csv(d / "carbon_support_raman.csv", index=False)

    # TMD MoS2 catalyst (300 to 500 cm^-1)
    shift_tmd = np.linspace(320.0, 480.0, 801)
    e2g = 650.0 / (1.0 + ((shift_tmd - 383.0) / 3.8) ** 2)
    a1g = 920.0 / (1.0 + ((shift_tmd - 408.0) / 4.2) ** 2)
    tmd_bg = 80.0 + np.random.normal(0, 6.0, len(shift_tmd))
    tmd_int = tmd_bg + e2g + a1g

    df_tmd = pd.DataFrame({"Raman_Shift": shift_tmd, "Intensity": np.round(tmd_int, 1)})
    for d in [EXAMPLES_DIR / "raman", FIXTURES_DIR / "raman"]:
        df_tmd.to_csv(d / "tmd_catalyst_raman.csv", index=False)

    # =========================================================================
    # 4. Microscopy Images (SEM and HRTEM)
    # =========================================================================
    print("Generating Microscopy images...")
    # (A) SEM Image: 512x512 pixels, pixel_size_nm = 0.50 nm/px (Field of view: 256 nm)
    sem_canvas = np.random.normal(0.65, 0.05, (512, 512))  # Carbon support texture
    # Add 45 dispersed nanoparticles
    np.random.seed(101)
    n_particles = 45
    x_coords = np.random.randint(40, 472, n_particles)
    y_coords = np.random.randint(40, 472, n_particles)
    radii_nm = np.random.lognormal(mean=1.8, sigma=0.25, size=n_particles)  # ~6-10 nm radius
    radii_px = radii_nm / 0.50

    yy, xx = np.mgrid[:512, :512]
    for xc, yc, r in zip(x_coords, y_coords, radii_px):
        dist = np.sqrt((xx - xc) ** 2 + (yy - yc) ** 2)
        # Dark particle on lighter carbon background
        mask = dist <= r
        sem_canvas[mask] = np.clip(sem_canvas[mask] - 0.45 * (1.0 - (dist[mask] / r) ** 2), 0.05, 1.0)

    sem_uint8 = np.clip(sem_canvas * 255.0, 0, 255).astype(np.uint8)
    sem_img = Image.fromarray(sem_uint8)

    for d in [EXAMPLES_DIR / "microscopy", FIXTURES_DIR / "microscopy"]:
        sem_img.save(d / "sem_nanoparticles.png")

    # (B) HRTEM Lattice Image: 256x256 pixels, pixel_size_nm = 0.025 nm/px (25 pm/px)
    # Target Pt (111) lattice fringes: d = 0.226 nm = 2.26 A
    # In pixels: period = 0.226 / 0.025 = 9.04 pixels
    fringe_period_px = 0.226 / 0.025
    yy_h, xx_h = np.mgrid[:256, :256]
    # Sine modulation with 10% gaussian noise and envelope
    fringes = 0.50 + 0.35 * np.cos(2.0 * np.pi * xx_h / fringe_period_px)
    hrtem_noise = np.random.normal(0, 0.06, (256, 256))
    hrtem_canvas = np.clip(fringes + hrtem_noise, 0.0, 1.0)
    hrtem_uint8 = (hrtem_canvas * 255.0).astype(np.uint8)
    hrtem_img = Image.fromarray(hrtem_uint8)

    for d in [EXAMPLES_DIR / "microscopy", FIXTURES_DIR / "microscopy"]:
        hrtem_img.save(d / "hrtem_lattice_crop.png")

    # =========================================================================
    # 5. ICP-OES Datasets (Calibration & Leaching)
    # =========================================================================
    print("Generating ICP analytical datasets...")
    # Multi-point calibration standards for Pt, Ni, Fe (0.0 to 10.0 mg/L)
    cal_conc = [0.0, 0.1, 0.5, 1.0, 2.5, 5.0, 10.0]
    cal_rows = []
    # Calibration slopes: Pt ~ 12500 cts/(mg/L), Ni ~ 28400 cts/(mg/L), Fe ~ 34200 cts/(mg/L)
    slopes = {"Pt": 12500.0, "Ni": 28400.0, "Fe": 34200.0}
    blanks = {"Pt": 45.0, "Ni": 82.0, "Fe": 115.0}

    for elem, m in slopes.items():
        b0 = blanks[elem]
        for c in cal_conc:
            int_meas = b0 + m * c + np.random.normal(0, 25.0)
            cal_rows.append({
                "sample_id": f"STD_{c}mgL",
                "element": elem,
                "concentration": c,
                "intensity": round(max(10.0, int_meas), 1),
                "unit": "mg/L",
            })
        # Add 3 blank replicates for LOD evaluation
        for bi in range(1, 4):
            cal_rows.append({
                "sample_id": f"BLANK_{bi}",
                "element": elem,
                "concentration": 0.0,
                "intensity": round(b0 + np.random.normal(0, 8.0), 1),
                "unit": "mg/L",
            })

    df_icp_cal = pd.DataFrame(cal_rows)
    for d in [EXAMPLES_DIR / "icp", FIXTURES_DIR / "icp"]:
        df_icp_cal.to_csv(d / "icp_calibration_standards.csv", index=False)

    # Catalyst Leaching stability study (OER testing in 1.0 M KOH, V = 50 mL, m_cat = 2.0 mg, 20 wt% metal)
    # Samples: Catalyst_A (pristine, 10h, 50h, 100h)
    leaching_rows = [
        {"sample_id": "Pristine_0h", "element": "Ni", "concentration": 0.012, "dilution_factor": 1.0, "unit": "mg/L"},
        {"sample_id": "Pristine_0h", "element": "Fe", "concentration": 0.008, "dilution_factor": 1.0, "unit": "mg/L"},
        {"sample_id": "OER_10h", "element": "Ni", "concentration": 0.145, "dilution_factor": 2.0, "unit": "mg/L"},
        {"sample_id": "OER_10h", "element": "Fe", "concentration": 0.082, "dilution_factor": 2.0, "unit": "mg/L"},
        {"sample_id": "OER_50h", "element": "Ni", "concentration": 0.420, "dilution_factor": 2.0, "unit": "mg/L"},
        {"sample_id": "OER_50h", "element": "Fe", "concentration": 0.265, "dilution_factor": 2.0, "unit": "mg/L"},
        {"sample_id": "OER_100h", "element": "Ni", "concentration": 0.890, "dilution_factor": 2.0, "unit": "mg/L"},
        {"sample_id": "OER_100h", "element": "Fe", "concentration": 0.580, "dilution_factor": 2.0, "unit": "mg/L"},
    ]
    df_leach = pd.DataFrame(leaching_rows)
    for d in [EXAMPLES_DIR / "icp", FIXTURES_DIR / "icp"]:
        df_leach.to_csv(d / "catalyst_leaching_study.csv", index=False)

    # =========================================================================
    # 6. BET / BJH N2 Adsorption Isotherm (77.35 K)
    # =========================================================================
    print("Generating BET / BJH isotherm benchmark datasets...")
    # Realistic Type IV isotherm with H1 hysteresis loop for mesoporous catalyst
    p_ads = np.array([
        0.01, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
        0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.98
    ])
    # Monolayer capacity ~35 cm^3/g STP -> S_BET ~ 152 m^2/g, C ~ 115
    q_m_true = 35.0
    c_true = 115.0

    def calc_base_q(p_vals):
        # BET isotherm capped physically at p=0.35 to prevent divergent multilayer infinity
        p_clamped = np.minimum(p_vals, 0.35)
        bet_q = (q_m_true * c_true * p_clamped) / ((1.0 - p_clamped) * (1.0 + (c_true - 1.0) * p_clamped))
        # Add small external surface multilayer contribution beyond 0.35
        ext_q = 20.0 * np.maximum(0.0, p_vals - 0.35)
        return bet_q + ext_q

    # Add mesopore capillary condensation step between p/p0 = 0.50 and 0.85
    meso_cap = 85.0 / (1.0 + np.exp(-18.0 * (p_ads - 0.68)))
    q_ads_tot = calc_base_q(p_ads) + meso_cap

    # Desorption branch with hysteresis loop (capillary evaporation centered around 0.58)
    p_des = np.array([
        0.98, 0.94, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35
    ])
    meso_des = 85.0 / (1.0 + np.exp(-20.0 * (p_des - 0.58)))
    q_des_tot = calc_base_q(p_des) + meso_des

    np.random.seed(42)
    iso_p = np.concatenate([p_ads, p_des])
    iso_q = np.concatenate([q_ads_tot, q_des_tot]) + np.random.normal(0, 0.02, len(iso_p))
    iso_branches = ["adsorption"] * len(p_ads) + ["desorption"] * len(p_des)

    df_isotherm = pd.DataFrame({
        "p_p0": np.round(iso_p, 4),
        "quantity_adsorbed": np.round(iso_q, 2),
        "branch": iso_branches,
    })

    for d in [EXAMPLES_DIR / "bet", FIXTURES_DIR / "bet"]:
        df_isotherm.to_csv(d / "mesoporous_catalyst_n2_isotherm.csv", index=False)

    print("Successfully generated all characterization benchmark fixtures.")

if __name__ == "__main__":
    generate_all_fixtures()

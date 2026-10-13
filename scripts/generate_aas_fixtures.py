"""
Deterministic Generator for Atomic Absorption Spectroscopy (AAS) Benchmark Datasets.

Generates realistic, physically validated synthetic datasets for:
1. Flame AAS (FAAS) Fe calibration standards with triplicates, blanks, and QC checks.
2. Graphite Furnace AAS (GFAAS) Pt integrated absorbance (A*s) calibration standards.
3. Method of Standard Additions (MOSA) for Co in high-salt catalyst matrix.
4. Catalyst metal loading acid digestion (Pt, Co, Ni) with nominal targets.
5. Electrochemical OER/HER stability metal leaching time-series with aliquot tracking.
6. Analytical QC run sequence with continuing calibration verifications (CCV).
"""

from pathlib import Path
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
EXAMPLES_AAS = PROJECT_ROOT / "data" / "examples" / "characterization" / "aas"
FIXTURES_AAS = PROJECT_ROOT / "tests" / "characterization" / "fixtures" / "aas"


def generate_aas_fixtures():
    np.random.seed(42)

    for target_dir in [EXAMPLES_AAS, FIXTURES_AAS]:
        target_dir.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # 1. Flame AAS (FAAS) Calibration: Iron (Fe) at 248.3 nm
    # =========================================================================
    # Standards: 0.0 (Blank), 0.20, 0.50, 1.00, 2.00, 5.00 mg/L in 0.1 M HNO3
    # Model: A = 0.0820 * C + 0.0015 + Gaussian noise
    stds_c = [0.0, 0.20, 0.50, 1.00, 2.00, 5.00]
    faas_rows = []
    slope_faas = 0.0825
    intercept_faas = 0.0018

    # Standards (triplicates)
    for c in stds_c:
        label = "BLANK_CAL" if c == 0.0 else f"STD_{c:.2f}"
        stype = "Blank" if c == 0.0 else "Standard"
        for rep in range(1, 4):
            noise = np.random.normal(0, 0.0008)
            abs_val = max(0.0001, slope_faas * c + intercept_faas + noise)
            faas_rows.append({
                "Sample_ID": label,
                "Sample_Type": stype,
                "Element": "Fe",
                "Wavelength_nm": 248.3,
                "Replicate": rep,
                "Concentration_Standard": c,
                "Absorbance": round(abs_val, 4),
                "Concentration_Reported": round(c + noise / slope_faas, 3) if c > 0 else 0.0,
                "Unit": "mg/L",
                "Dilution_Factor": 1.0,
                "Technique": "FAAS",
                "Background_Correction": "D2",
            })

    # Add ICV and CCV
    for rep in range(1, 4):
        noise = np.random.normal(0, 0.0007)
        c_icv = 1.00
        abs_icv = slope_faas * c_icv + intercept_faas + noise
        faas_rows.append({
            "Sample_ID": "ICV_1.0mgL",
            "Sample_Type": "ICV",
            "Element": "Fe",
            "Wavelength_nm": 248.3,
            "Replicate": rep,
            "Concentration_Standard": 1.00,
            "Absorbance": round(abs_icv, 4),
            "Concentration_Reported": round(c_icv + noise / slope_faas, 3),
            "Unit": "mg/L",
            "Dilution_Factor": 1.0,
            "Technique": "FAAS",
            "Background_Correction": "D2",
        })

    df_faas = pd.DataFrame(faas_rows)

    # =========================================================================
    # 2. Graphite Furnace AAS (GFAAS) Calibration: Platinum (Pt) at 265.9 nm
    # =========================================================================
    # Standards: 0.0, 5.0, 10.0, 25.0, 50.0, 100.0 ug/L
    # Signal: Integrated Absorbance (A*s)
    # Model: Int_Abs = 0.00235 * C + 0.0012 + noise
    gfaas_stds = [0.0, 5.0, 10.0, 25.0, 50.0, 100.0]
    gfaas_rows = []
    slope_gfaas = 0.00238
    int_gfaas = 0.0011

    for c in gfaas_stds:
        label = "BLANK_GFAAS" if c == 0.0 else f"STD_Pt_{c:.0f}ppb"
        stype = "Blank" if c == 0.0 else "Standard"
        for rep in range(1, 3):
            noise = np.random.normal(0, 0.0009)
            int_abs = max(0.0001, slope_gfaas * c + int_gfaas + noise)
            gfaas_rows.append({
                "Sample_ID": label,
                "Sample_Type": stype,
                "Element": "Pt",
                "Wavelength_nm": 265.9,
                "Replicate": rep,
                "Concentration_Standard": c,
                "Integrated_Absorbance": round(int_abs, 4),
                "Concentration_Reported": round(c + noise / slope_gfaas, 2) if c > 0 else 0.0,
                "Unit": "ug/L",
                "Dilution_Factor": 1.0,
                "Technique": "GFAAS",
                "Background_Correction": "Zeeman",
            })

    df_gfaas = pd.DataFrame(gfaas_rows)

    # =========================================================================
    # 3. Method of Standard Additions (MOSA): Cobalt in Alkaline Electrolyte
    # =========================================================================
    # Sample has true concentration in analyzed solution C_app = 0.85 mg/L
    # 5.0 mL sample in 25.0 mL flask -> DF = 5.0 -> True original = 4.25 mg/L
    # Spikes added to final flask: [0.0, 0.50, 1.00, 1.50, 2.00] mg/L
    # Matrix sensitivity: slope = 0.0540 (suppressed vs pure water 0.0710)
    # S = 0.0540 * (C_app + Delta C) = 0.0540 * Delta C + 0.0459
    c_true_app = 0.85
    slope_mosa = 0.0542
    intercept_mosa = slope_mosa * c_true_app  # ~ 0.0461
    mosa_spikes = [0.0, 0.50, 1.00, 1.50, 2.00]
    mosa_rows = []

    for sp in mosa_spikes:
        label = "Sample_Unspiked" if sp == 0.0 else f"Spike_+{sp:.2f}mgL"
        for rep in range(1, 4):
            noise = np.random.normal(0, 0.0006)
            sig = max(0.001, slope_mosa * sp + intercept_mosa + noise)
            mosa_rows.append({
                "Sample_ID": label,
                "Sample_Type": "StandardAddition",
                "Element": "Co",
                "Wavelength_nm": 240.7,
                "Replicate": rep,
                "Added_Concentration": sp,
                "Absorbance": round(sig, 4),
                "Sample_Aliquot_mL": 5.0,
                "Final_Volume_mL": 25.0,
                "Unit": "mg/L",
                "Technique": "FAAS",
                "Background_Correction": "D2",
            })

    df_mosa = pd.DataFrame(mosa_rows)

    # =========================================================================
    # 4. Catalyst Metal Loading Digestion Dataset
    # =========================================================================
    # Catalysts:
    # 1. 20 wt% Pt/C (nominal): 2.50 mg digested in 50.0 mL, DF = 2.0 -> conc ~ 4.95 mg/L Pt -> 19.8 wt% actual
    # 2. 15 wt% Pt - 5 wt% Co / C: 3.00 mg digested in 50.0 mL, DF = 2.0 -> conc Pt ~ 4.41 mg/L (14.7 wt%), Co ~ 1.53 mg/L (5.1 wt%)
    # 3. NiFe-LDH: 4.00 mg digested in 50.0 mL, DF = 5.0 -> Ni conc ~ 3.92 mg/L (24.5 wt%), Fe conc ~ 1.96 mg/L (12.25 wt%)
    digestion_rows = [
        # Blank
        {"Sample_ID": "PROC_BLANK", "Sample_Type": "Blank", "Element": "Pt", "Mass_mg": 0.0, "Digest_Vol_mL": 50.0, "Dilution_Factor": 1.0, "Measured_Conc_mg_L": 0.001, "Unit": "mg/L", "Nominal_wt_pct": 0.0},
        {"Sample_ID": "PROC_BLANK", "Sample_Type": "Blank", "Element": "Co", "Mass_mg": 0.0, "Digest_Vol_mL": 50.0, "Dilution_Factor": 1.0, "Measured_Conc_mg_L": 0.001, "Unit": "mg/L", "Nominal_wt_pct": 0.0},
        {"Sample_ID": "PROC_BLANK", "Sample_Type": "Blank", "Element": "Ni", "Mass_mg": 0.0, "Digest_Vol_mL": 50.0, "Dilution_Factor": 1.0, "Measured_Conc_mg_L": 0.002, "Unit": "mg/L", "Nominal_wt_pct": 0.0},
        # Pt/C
        {"Sample_ID": "Pt_C_Pristine_Rep1", "Sample_Type": "Sample", "Element": "Pt", "Mass_mg": 2.50, "Digest_Vol_mL": 50.0, "Dilution_Factor": 2.0, "Measured_Conc_mg_L": 4.94, "Unit": "mg/L", "Nominal_wt_pct": 20.0},
        {"Sample_ID": "Pt_C_Pristine_Rep2", "Sample_Type": "Sample", "Element": "Pt", "Mass_mg": 2.48, "Digest_Vol_mL": 50.0, "Dilution_Factor": 2.0, "Measured_Conc_mg_L": 4.91, "Unit": "mg/L", "Nominal_wt_pct": 20.0},
        {"Sample_ID": "Pt_C_Pristine_Rep3", "Sample_Type": "Sample", "Element": "Pt", "Mass_mg": 2.52, "Digest_Vol_mL": 50.0, "Dilution_Factor": 2.0, "Measured_Conc_mg_L": 4.99, "Unit": "mg/L", "Nominal_wt_pct": 20.0},
        # PtCo/C
        {"Sample_ID": "PtCo_C_Rep1", "Sample_Type": "Sample", "Element": "Pt", "Mass_mg": 3.00, "Digest_Vol_mL": 50.0, "Dilution_Factor": 2.0, "Measured_Conc_mg_L": 4.41, "Unit": "mg/L", "Nominal_wt_pct": 15.0},
        {"Sample_ID": "PtCo_C_Rep2", "Sample_Type": "Sample", "Element": "Pt", "Mass_mg": 3.02, "Digest_Vol_mL": 50.0, "Dilution_Factor": 2.0, "Measured_Conc_mg_L": 4.45, "Unit": "mg/L", "Nominal_wt_pct": 15.0},
        {"Sample_ID": "PtCo_C_Rep1", "Sample_Type": "Sample", "Element": "Co", "Mass_mg": 3.00, "Digest_Vol_mL": 50.0, "Dilution_Factor": 2.0, "Measured_Conc_mg_L": 1.53, "Unit": "mg/L", "Nominal_wt_pct": 5.0},
        {"Sample_ID": "PtCo_C_Rep2", "Sample_Type": "Sample", "Element": "Co", "Mass_mg": 3.02, "Digest_Vol_mL": 50.0, "Dilution_Factor": 2.0, "Measured_Conc_mg_L": 1.54, "Unit": "mg/L", "Nominal_wt_pct": 5.0},
        # NiFe-LDH
        {"Sample_ID": "NiFe_LDH_Rep1", "Sample_Type": "Sample", "Element": "Ni", "Mass_mg": 4.00, "Digest_Vol_mL": 50.0, "Dilution_Factor": 5.0, "Measured_Conc_mg_L": 3.92, "Unit": "mg/L", "Nominal_wt_pct": 25.0},
        {"Sample_ID": "NiFe_LDH_Rep2", "Sample_Type": "Sample", "Element": "Ni", "Mass_mg": 3.98, "Digest_Vol_mL": 50.0, "Dilution_Factor": 5.0, "Measured_Conc_mg_L": 3.90, "Unit": "mg/L", "Nominal_wt_pct": 25.0},
        {"Sample_ID": "NiFe_LDH_Rep1", "Sample_Type": "Sample", "Element": "Fe", "Mass_mg": 4.00, "Digest_Vol_mL": 50.0, "Dilution_Factor": 5.0, "Measured_Conc_mg_L": 1.96, "Unit": "mg/L", "Nominal_wt_pct": 12.0},
        {"Sample_ID": "NiFe_LDH_Rep2", "Sample_Type": "Sample", "Element": "Fe", "Mass_mg": 3.98, "Digest_Vol_mL": 50.0, "Dilution_Factor": 5.0, "Measured_Conc_mg_L": 1.95, "Unit": "mg/L", "Nominal_wt_pct": 12.0},
    ]
    df_digestion = pd.DataFrame(digestion_rows)

    # =========================================================================
    # 5. Catalyst Metal Leaching Time-Course (HER/OER Stability)
    # =========================================================================
    # Initial: 2.0 mg of PtCo/C (14.8 wt% Pt, 5.1 wt% Co) -> 296 ug Pt, 102 ug Co
    # In 50.0 mL 0.1 M KOH. 2.0 mL aliquot taken at each interval.
    # Time points: 0h, 10h, 25h, 50h, 100h
    leaching_times = [0, 10, 25, 50, 100]
    # Dissolved concentrations in cell (mg/L):
    # Co leaches more rapidly than Pt in alkaline OER
    pt_concs = [0.002, 0.015, 0.038, 0.075, 0.130]  # mg/L
    co_concs = [0.003, 0.045, 0.110, 0.230, 0.420]  # mg/L

    leach_rows = []
    for t, cpt, cco in zip(leaching_times, pt_concs, co_concs):
        leach_rows.append({
            "Time_hours": t,
            "Sample_ID": f"Electrolyte_{t}h",
            "Element": "Pt",
            "Electrolyte_Volume_mL": 50.0,
            "Aliquot_Withdrawn_mL": 2.0 if t > 0 else 0.0,
            "Measured_Conc_mg_L": cpt,
            "Dilution_Factor": 1.0,
            "Initial_Catalyst_Mass_mg": 2.0,
            "Initial_Metal_wt_pct": 14.8,
            "Technique": "GFAAS",
        })
        leach_rows.append({
            "Time_hours": t,
            "Sample_ID": f"Electrolyte_{t}h",
            "Element": "Co",
            "Electrolyte_Volume_mL": 50.0,
            "Aliquot_Withdrawn_mL": 2.0 if t > 0 else 0.0,
            "Measured_Conc_mg_L": cco,
            "Dilution_Factor": 1.0,
            "Initial_Catalyst_Mass_mg": 2.0,
            "Initial_Metal_wt_pct": 5.1,
            "Technique": "FAAS",
        })

    df_leach = pd.DataFrame(leach_rows)

    # =========================================================================
    # 6. Analytical QC Run Sequence & CCV Drift Tracking
    # =========================================================================
    # EPA Method 7000B QC sequence with 6 CCV points (Target = 1.00 mg/L)
    ccv_rec = [100.8, 98.6, 102.1, 99.4, 96.8, 101.5]
    qc_rows = []
    for idx, rec in enumerate(ccv_rec, start=1):
        qc_rows.append({
            "QC_Check": f"CCV_{idx}",
            "Standard_Target_mg_L": 1.00,
            "Measured_Conc_mg_L": round(rec / 100.0, 4),
            "Recovery_pct": rec,
            "SOP": "EPA_7000B",
            "Acceptance_Min_pct": 90.0,
            "Acceptance_Max_pct": 110.0,
            "Status": "PASS",
        })

    df_qc = pd.DataFrame(qc_rows)

    # =========================================================================
    # Save all datasets to both examples and fixtures
    # =========================================================================
    for target_dir in [EXAMPLES_AAS, FIXTURES_AAS]:
        df_faas.to_csv(target_dir / "faas_calibration_fe.csv", index=False)
        df_gfaas.to_csv(target_dir / "gfaas_calibration_pt.csv", index=False)
        df_mosa.to_csv(target_dir / "aas_standard_addition_catalyst.csv", index=False)
        df_digestion.to_csv(target_dir / "catalyst_digestion_loading.csv", index=False)
        df_leach.to_csv(target_dir / "catalyst_stability_leaching.csv", index=False)
        df_qc.to_csv(target_dir / "aas_qc_run_sequence.csv", index=False)

        # Also save TSV and Excel files to verify multi-format parser support
        df_faas.to_csv(target_dir / "faas_calibration_fe.tsv", sep="\t", index=False)
        df_faas.to_excel(target_dir / "faas_calibration_fe.xlsx", index=False)

    print("All AAS benchmark datasets successfully generated.")


if __name__ == "__main__":
    generate_aas_fixtures()

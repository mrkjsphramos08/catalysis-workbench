"""
Synthetic Electrochemical Data Generator for Metrohm Autolab and BioLogic Formats.

Generates realistic but explicitly synthetic test fixtures and example datasets:
- Metrohm Autolab NOVA format (.txt, .csv, multi-cycle, binary .nox mock)
- BioLogic EC-Lab format (.mpt ASCII, .mps settings, multi-cycle)
- Covers LSV (HER & OER), CV, EIS, Chronoamperometry, Chronopotentiometry.
"""

from pathlib import Path
import numpy as np
import pandas as pd


def generate_all_electrochem_datasets(base_dirs: list[Path]) -> None:
    """Generate identical synthetic datasets across the specified target directories."""
    for base_dir in base_dirs:
        base_dir.mkdir(parents=True, exist_ok=True)
        metrohm_dir = base_dir / "metrohm"
        biologic_dir = base_dir / "biologic"
        metrohm_dir.mkdir(parents=True, exist_ok=True)
        biologic_dir.mkdir(parents=True, exist_ok=True)

        # -------------------------------------------------------------------------
        # 1. METROHM LSV (OER on NiFe-LDH in 1.0 M KOH, pH 14, Ag/AgCl sat)
        # -------------------------------------------------------------------------
        # E vs Ag/AgCl from 0.10 to 0.80 V (1.10 to 1.80 V vs RHE)
        e_agcl = np.linspace(0.10, 0.75, 260)
        t_lsv = np.linspace(0, 130, 260)  # scan rate 5 mV/s
        # Butler-Volmer exponential onset near 0.45 V vs Ag/AgCl (1.45 V vs RHE)
        eta_oer = np.maximum(0.0, e_agcl - 0.45)
        # Area = 0.196 cm^2 -> j from 0 to 60 mA/cm^2 -> I = j * Area
        j_mA_cm2 = 0.05 + 0.1 * np.exp(eta_oer * 18.0) * (eta_oer > 0)
        # Add slight noise
        np.random.seed(42)
        j_mA_cm2 += np.random.normal(0, 0.03, size=len(j_mA_cm2))
        i_A = (j_mA_cm2 * 1e-3) * 0.196

        metrohm_lsv_txt = (
            "# Metrohm Autolab NOVA Export (SYNTHETIC BENCHMARK DATASET)\n"
            "# Instrument: Metrohm Autolab PGSTAT302N\n"
            "# Software: NOVA 2.1.5\n"
            "# Procedure: Linear Sweep Voltammetry (OER)\n"
            "# Working Electrode: NiFe-LDH\n"
            "# Reference Electrode: Ag/AgCl (3 M KCl)\n"
            "# Electrolyte: 1.0 M KOH (pH 14.0)\n"
            "# Area: 0.196 cm^2\n"
            "# Scan Rate: 5 mV/s\n"
            "# Ru: 4.2 Ohm\n"
            "Time (s)\tWE(1).Potential (V)\tWE(1).Current (A)\n"
        )
        for t, v, cur in zip(t_lsv, e_agcl, i_A):
            metrohm_lsv_txt += f"{t:.3f}\t{v:.5f}\t{cur:.8e}\n"
        (metrohm_dir / "metrohm_oer_lsv.txt").write_text(metrohm_lsv_txt, encoding="utf-8")

        # -------------------------------------------------------------------------
        # 2. METROHM CV (Multi-cycle, 3 cycles, 50 mV/s)
        # -------------------------------------------------------------------------
        # Triangle wave 0.0 -> 0.8 -> 0.0 V
        t_half = 16.0  # 0.8 V at 50 mV/s takes 16 s
        cycle_records = []
        for c in range(1, 4):
            # Forward scan
            t_fwd = np.linspace(0, t_half, 100)
            v_fwd = np.linspace(0.0, 0.8, 100)
            # Reverse scan
            t_rev = np.linspace(t_half, 2 * t_half, 100)
            v_rev = np.linspace(0.8, 0.0, 100)
            
            t_cyc = np.concatenate([t_fwd, t_rev]) + (c - 1) * 32.0
            v_cyc = np.concatenate([v_fwd, v_rev])
            # Redox peaks: oxidation peak at 0.48 V, reduction peak at 0.38 V
            i_fwd = 0.05e-3 + 0.4e-3 * np.exp(-((v_fwd - 0.48) ** 2) / (2 * 0.04**2))
            i_rev = -0.05e-3 - 0.38e-3 * np.exp(-((v_rev - 0.38) ** 2) / (2 * 0.04**2))
            i_cyc = np.concatenate([i_fwd, i_rev]) + np.random.normal(0, 5e-6, size=200)

            for t_pt, v_pt, i_pt in zip(t_cyc, v_cyc, i_cyc):
                cycle_records.append((t_pt, v_pt, i_pt, c))

        metrohm_cv_csv = (
            "# Metrohm Autolab NOVA Multi-Cycle CV (SYNTHETIC DATASET)\n"
            "# Scan rate: 0.050 V/s\n"
            "# Reference: Ag/AgCl\n"
            "# Area: 0.196 cm^2\n"
            "Time (s),WE(1).Potential (V),WE(1).Current (A),Scan\n"
        )
        for r in cycle_records:
            metrohm_cv_csv += f"{r[0]:.3f},{r[1]:.5f},{r[2]:.8e},{r[3]}\n"
        (metrohm_dir / "metrohm_cv_multicycle.csv").write_text(metrohm_cv_csv, encoding="utf-8")

        # -------------------------------------------------------------------------
        # 3. METROHM EIS (100 kHz to 0.1 Hz, Rs=5.2 Ohm, Rct=45.0 Ohm, CPE)
        # -------------------------------------------------------------------------
        freq_eis = np.logspace(5, -1, 60)
        rs = 5.2
        rct = 45.0
        q_cpe = 3.5e-5
        n_cpe = 0.88
        omega = 2 * np.pi * freq_eis
        z_cpe = 1.0 / (q_cpe * (1j * omega) ** n_cpe)
        z_parallel = (rct * z_cpe) / (rct + z_cpe)
        z_tot = rs + z_parallel
        # Add tiny realistic noise while strictly keeping physical capacitive impedance negative
        zr = z_tot.real + np.random.normal(0, 0.02, size=len(freq_eis))
        zi = np.minimum(-0.01, z_tot.imag + np.random.normal(0, 0.02, size=len(freq_eis)))
        neg_zi = -zi  # Metrohm exports -Z'' as positive

        metrohm_eis_txt = (
            "# Metrohm Autolab NOVA FRA EIS Export (SYNTHETIC DATASET)\n"
            "# Procedure: FRA Impedance\n"
            "# Amplitude: 10 mV RMS\n"
            "# DC Potential: 0.55 V vs Ag/AgCl\n"
            "# Area: 0.196 cm^2\n"
            "Frequency (Hz)\tZ' (Ohm)\t-Z'' (Ohm)\n"
        )
        for f_val, zr_val, neg_zi_val in zip(freq_eis, zr, neg_zi):
            metrohm_eis_txt += f"{f_val:.4e}\t{zr_val:.5f}\t{neg_zi_val:.5f}\n"
        (metrohm_dir / "metrohm_fra_eis.txt").write_text(metrohm_eis_txt, encoding="utf-8")

        # -------------------------------------------------------------------------
        # 4. METROHM CA (1 hour potentiostatic hold at 0.65 V)
        # -------------------------------------------------------------------------
        t_ca = np.linspace(0, 3600, 361)
        # Initial capacitive spike decaying over first 10 s, followed by slow 4% drift
        j_ca = 22.0 * np.exp(-t_ca / 3.0) + 20.0 * (1.0 - 0.04 * (t_ca / 3600.0))
        j_ca += np.random.normal(0, 0.04, size=len(t_ca))
        i_ca = (j_ca * 1e-3) * 0.196

        metrohm_ca_txt = (
            "# Metrohm Autolab Chronoamperometry (SYNTHETIC DATASET)\n"
            "# Procedure: Chronoamperometry\n"
            "# Potential: 0.65 V vs Ag/AgCl\n"
            "# Area: 0.196 cm^2\n"
            "Time (s)\tWE(1).Potential (V)\tWE(1).Current (A)\n"
        )
        for t_val, cur_val in zip(t_ca, i_ca):
            metrohm_ca_txt += f"{t_val:.1f}\t0.65000\t{cur_val:.8e}\n"
        (metrohm_dir / "metrohm_stability_ca.txt").write_text(metrohm_ca_txt, encoding="utf-8")

        # -------------------------------------------------------------------------
        # 5. METROHM CP (1 hour galvanostatic hold at 10 mA/cm^2)
        # -------------------------------------------------------------------------
        t_cp = np.linspace(0, 3600, 361)
        # Potential begins at 0.52 V vs Ag/AgCl and drifts +1.8 mV over 1 hour
        v_cp = 0.520 + 0.0018 * (t_cp / 3600.0) + np.random.normal(0, 0.0003, size=len(t_cp))
        i_cp_fixed = 0.010 * 0.196  # 1.96 mA

        metrohm_cp_txt = (
            "# Metrohm Autolab Chronopotentiometry (SYNTHETIC DATASET)\n"
            "# Procedure: Chronopotentiometry\n"
            "# Applied Current: 1.96 mA (10 mA/cm^2)\n"
            "# Area: 0.196 cm^2\n"
            "Time (s)\tWE(1).Current (A)\tWE(1).Potential (V)\n"
        )
        for t_val, v_val in zip(t_cp, v_cp):
            metrohm_cp_txt += f"{t_val:.1f}\t{i_cp_fixed:.8e}\t{v_val:.5f}\n"
        (metrohm_dir / "metrohm_stability_cp.txt").write_text(metrohm_cp_txt, encoding="utf-8")

        # -------------------------------------------------------------------------
        # 6. METROHM NOX MOCK (Binary mock to test guard)
        # -------------------------------------------------------------------------
        nox_bytes = b"NOVA_BINARY_PROJECT_DATA\x00\x01\x02\x03\x04"
        (metrohm_dir / "mock_procedure.nox").write_bytes(nox_bytes)

        # -------------------------------------------------------------------------
        # 7. BIOLOGIC LSV (.mpt ASCII format)
        # -------------------------------------------------------------------------
        # Matched NiFe-LDH in 1.0 M KOH, E vs Ag/AgCl from 0.10 to 0.75 V
        # BioLogic EC-Lab format:
        # Line 1: EC-Lab ASCII FILE
        # Line 2: Nb header lines : 11
        # Metadata lines...
        # Column line with tabs: mode	ox/red	error	control changes	time/s	control/V	Ewe/V	I/mA	cycle number
        bio_lsv_header = (
            "EC-Lab ASCII FILE\n"
            "Nb header lines : 11\n"
            "Device : VSP-300\n"
            "Electrode surface area : 0.196 cm2\n"
            "Characteristic mass : 0.200 mg\n"
            "Electrolyte : 1.0 M KOH (pH 14.0)\n"
            "Reference electrode : Ag/AgCl (3 M KCl)\n"
            "Technique : LSV\n"
            "Scan rate : 5.0 mV/s\n"
            "Uncompensated resistance : 4.0 Ohm\n"
            "time/s\tcontrol/V\tEwe/V\tI/mA\tcycle number\n"
        )
        # Current in mA matching Metrohm with tiny instrumentation offset
        j_bio_mA = (0.05 + 0.102 * np.exp(eta_oer * 17.9) * (eta_oer > 0) + np.random.normal(0, 0.02, size=len(eta_oer))) * 0.196
        bio_lsv_body = []
        for t_val, v_val, i_val in zip(t_lsv, e_agcl, j_bio_mA):
            bio_lsv_body.append(f"{t_val:.3f}\t{v_val:.5f}\t{v_val:.5f}\t{i_val:.6f}\t1")
        (biologic_dir / "biologic_oer_lsv.mpt").write_text(bio_lsv_header + "\n".join(bio_lsv_body) + "\n", encoding="utf-8")

        # -------------------------------------------------------------------------
        # 8. BIOLOGIC HER LSV (.mpt ASCII format)
        # -------------------------------------------------------------------------
        # Pt/C in 0.5 M H2SO4 (pH 0.3), E vs Ag/AgCl from 0.0 to -0.45 V (-0.20 to -0.65 vs RHE)
        # HER current is cathodic (-I)
        e_her_agcl = np.linspace(0.0, -0.45, 225)
        t_her = np.linspace(0, 90, 225)
        eta_her = np.maximum(0.0, -e_her_agcl - 0.20)
        j_her_mA = (-0.02 - 0.12 * np.exp(eta_her * 22.0) * (eta_her > 0) + np.random.normal(0, 0.02, size=len(eta_her))) * 0.196
        bio_her_header = (
            "EC-Lab ASCII FILE\n"
            "Nb header lines : 10\n"
            "Device : SP-150\n"
            "Electrode surface area : 0.196 cm2\n"
            "Electrolyte : 0.5 M H2SO4 (pH 0.3)\n"
            "Reference electrode : Ag/AgCl\n"
            "Technique : LSV (HER)\n"
            "Scan rate : 5.0 mV/s\n"
            "Ru : 5.0 Ohm\n"
            "time/s\tcontrol/V\tEwe/V\tI/mA\tcycle number\n"
        )
        bio_her_body = [f"{t:.3f}\t{v:.5f}\t{v:.5f}\t{i:.6f}\t1" for t, v, i in zip(t_her, e_her_agcl, j_her_mA)]
        (biologic_dir / "biologic_her_lsv.mpt").write_text(bio_her_header + "\n".join(bio_her_body) + "\n", encoding="utf-8")

        # -------------------------------------------------------------------------
        # 9. BIOLOGIC EIS (.mpt ASCII format)
        # -------------------------------------------------------------------------
        bio_eis_header = (
            "EC-Lab ASCII FILE\n"
            "Nb header lines : 8\n"
            "Device : VSP-300\n"
            "Electrode surface area : 0.196 cm2\n"
            "Technique : PEIS\n"
            "Ewe = 0.55 V vs Ag/AgCl\n"
            "Va = 10.0 mV\n"
            "freq/Hz\tRe(Z)/Ohm\t-Im(Z)/Ohm\t|Z|/Ohm\tPhase(Z)/deg\n"
        )
        mag_eis = np.hypot(zr, zi)
        phase_eis = np.degrees(np.arctan2(zi, zr))
        bio_eis_body = [
            f"{f:.4e}\t{zr_val:.5f}\t{neg_zi_val:.5f}\t{m_val:.5f}\t{ph_val:.3f}"
            for f, zr_val, neg_zi_val, m_val, ph_val in zip(freq_eis, zr, neg_zi, mag_eis, phase_eis)
        ]
        (biologic_dir / "biologic_peis.mpt").write_text(bio_eis_header + "\n".join(bio_eis_body) + "\n", encoding="utf-8")

        # -------------------------------------------------------------------------
        # 10. BIOLOGIC CA & CP (.mpt ASCII format)
        # -------------------------------------------------------------------------
        bio_ca_header = (
            "EC-Lab ASCII FILE\n"
            "Nb header lines : 7\n"
            "Device : VSP-300\n"
            "Electrode surface area : 0.196 cm2\n"
            "Technique : CA\n"
            "E = 0.65 V vs Ag/AgCl\n"
            "time/s\tcontrol/V\tEwe/V\tI/mA\n"
        )
        bio_ca_body = [f"{t:.1f}\t0.65000\t0.65000\t{cur * 1e3:.6f}" for t, cur in zip(t_ca, i_ca)]
        (biologic_dir / "biologic_stability_ca.mpt").write_text(bio_ca_header + "\n".join(bio_ca_body) + "\n", encoding="utf-8")

        bio_cp_header = (
            "EC-Lab ASCII FILE\n"
            "Nb header lines : 7\n"
            "Device : VSP-300\n"
            "Electrode surface area : 0.196 cm2\n"
            "Technique : CP\n"
            "I = 1.96 mA\n"
            "time/s\tcontrol/mA\tEwe/V\tI/mA\n"
        )
        bio_cp_body = [f"{t:.1f}\t1.96000\t{v:.5f}\t1.96000" for t, v in zip(t_cp, v_cp)]
        (biologic_dir / "biologic_stability_cp.mpt").write_text(bio_cp_header + "\n".join(bio_cp_body) + "\n", encoding="utf-8")

        # -------------------------------------------------------------------------
        # 11. BIOLOGIC MPS SETTINGS
        # -------------------------------------------------------------------------
        bio_mps = (
            "EC-Lab Settings File\n"
            "Technique : LSV\n"
            "Eini (V) : 0.100\n"
            "Efin (V) : 0.750\n"
            "dE/dt (mV/s) : 5.0\n"
            "Electrode surface area (cm2) : 0.196\n"
            "Characteristic mass (mg) : 0.200\n"
        )
        (biologic_dir / "biologic_settings.mps").write_text(bio_mps, encoding="utf-8")


if __name__ == "__main__":
    test_fixtures_dir = Path("tests/fixtures/electrochemistry")
    examples_dir = Path("data/examples/electrochemistry")
    generate_all_electrochem_datasets([test_fixtures_dir, examples_dir])
    print("Generated all synthetic electrochemistry fixtures and example files successfully.")

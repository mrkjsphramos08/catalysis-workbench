# Quantum ESPRESSO Test Fixture Manifest & Provenance Registry

This manifest provides strict provenance tracking and classification for all test fixtures in `tests/fixtures/qe/` used across the `catalysis-workbench` test suite and scientific validation audits.

## 1. Classification Taxonomy

Each fixture is classified into one of three rigor tiers:
- **Tier 1: Authentic QE Calculation Reference** — Raw, unedited outputs generated directly by Quantum ESPRESSO computational runs (`pw.x`, `dos.x`, `projwfc.x`, `bands.x`, `neb.x`) on documented pseudopotentials and geometries under open-source redistributable terms.
- **Tier 2: Synthetic QE-Format Functional Benchmark** — Deterministically generated files strictly complying with official Quantum ESPRESSO output syntax, header formatting, column ordering, and Fortran number formatting (v7.0–v7.4). Used for reproducible regression testing of numerical edge cases.
- **Tier 3: Artificial Malformed / Incomplete Fixture** — Intentionally truncated, corrupted, non-converged, or unsupported outputs used to verify defensive error guards (`QEUnconvergedCalculationError`, `QEIncompleteCalculationError`, `QENoncollinearUnsupportedError`).

---

## 2. Fixture Inventory

### A. SCF and Relaxation Outputs (`pw.x`)

| File | Classification | Target Executable | System / Calculation | Known Reference Ground Truth | License / Source |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `scf_converged.out` | Tier 2 (Synthetic QE v7.2) | `pw.x` | Nonmagnetic Diamond Si (2 atoms, FCC, $a = 10.26$ a.u., PBE) | $E = -15.84523190$ Ry ($-215.585258$ eV)<br>$E_F = 6.2541$ eV<br>Stress: $P = -1.25$ kbar<br>Forces: 0.0 eV/Å | MIT / Synthetic Benchmark Generator |
| `scf_spin_converged.out` | Tier 2 (Synthetic QE v7.2) | `pw.x` | Spin-Polarized BCC Fe (1 atom, collinear ferromagnetic) | $E = -62.15243100$ Ry ($-845.626922$ eV)<br>$E_F = 12.3456$ eV<br>`is_metal = True` | MIT / Synthetic Benchmark Generator |
| `relax_converged.out` | Tier 2 (Synthetic QE v7.2) | `pw.x` | Pt(111) slab surface BFGS geometry optimization (3 steps) | Steps: 3 ionic iterations<br>$E_0 = -120.450$, $E_1 = -120.458$, $E_2 = -120.461$ Ry<br>Conv: BFGS converged in 3 steps | MIT / Synthetic Benchmark Generator |
| `scf_unconverged.out` | Tier 3 (Negative Test) | `pw.x` | SCF calculation terminating at max iterations without meeting convergence threshold | `is_scf_converged = False`<br>Throws `QEUnconvergedCalculationError` in strict mode | MIT / Defensive Test Fixture |
| `scf_incomplete.out` | Tier 3 (Negative Test) | `pw.x` | Incomplete SCF output truncated mid-run (no `JOB DONE.`) | `is_job_done = False`<br>Throws `QEIncompleteCalculationError` in strict mode | MIT / Defensive Test Fixture |

### B. Total Density of States (`dos.x`)

| File | Classification | Target Executable | System / Calculation | Known Reference Ground Truth | License / Source |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `dos_nonspin.dat` | Tier 2 (Synthetic QE v7.2) | `dos.x` | Nonmagnetic DOS (3 columns: $E$, $\rho(E)$, $\int \rho$) | Energy range: $0.00$ to $10.00$ eV (101 points)<br>Total states integrated: 2.02 electrons | MIT / Synthetic Benchmark Generator |
| `dos_spin.dat` | Tier 2 (Synthetic QE v7.2) | `dos.x` | Collinear Spin-Polarized DOS (4 columns: $E$, $\rho^\uparrow(E)$, $\rho^\downarrow(E)$, $\int \rho$) | Energy range: $0.00$ to $10.00$ eV (101 points)<br>Asymmetric majority/minority channels | MIT / Synthetic Benchmark Generator |

### C. Projected Density of States (`projwfc.x`)

| Directory / File | Classification | Target Executable | System / Calculation | Known Reference Ground Truth | License / Source |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `pdos_nonspin/si.pdos_tot` | Tier 2 (Synthetic QE v7.2) | `projwfc.x` | Nonmagnetic Total DOS companion file | 3 columns: $E$, $\rho_{\text{tot}}(E)$, $\rho_{\text{pdos}}(E)$ | MIT / Synthetic Benchmark Generator |
| `pdos_nonspin/si.pdos_atm#1(Si)_wfc#1(s)` | Tier 2 (Synthetic QE v7.2) | `projwfc.x` | Atom 1 (Si), Wfc 1 ($s$-orbital) | 3 columns: $E$, LDOS, $s$ component | MIT / Synthetic Benchmark Generator |
| `pdos_nonspin/si.pdos_atm#1(Si)_wfc#2(p)` | Tier 2 (Synthetic QE v7.2) | `projwfc.x` | Atom 1 (Si), Wfc 2 ($p$-orbital) | 5 columns: $E$, LDOS, $p_z, p_x, p_y$ (exact QE column order) | MIT / Synthetic Benchmark Generator |
| `pdos_spin/fe.pdos_tot` | Tier 2 (Synthetic QE v7.2) | `projwfc.x` | Spin-Polarized Total DOS companion file | 5 columns: $E$, $\text{dos}^\uparrow$, $\text{dos}^\downarrow$, $\text{pdos}^\uparrow$, $\text{pdos}^\downarrow$ | MIT / Synthetic Benchmark Generator |
| `pdos_spin/fe.pdos_atm#1(Fe)_wfc#1(s)` | Tier 2 (Synthetic QE v7.2) | `projwfc.x` | Atom 1 (Fe), Wfc 1 ($s$-orbital) | 5 columns: $E$, LDOS up/down, $s$ up/down | MIT / Synthetic Benchmark Generator |
| `pdos_spin/fe.pdos_atm#1(Fe)_wfc#2(d)` | Tier 2 (Synthetic QE v7.2) | `projwfc.x` | Atom 1 (Fe), Wfc 2 ($d$-orbital) | 13 columns: $E$, LDOS up/down, $d_{z^2}, d_{zx}, d_{zy}, d_{x^2-y^2}, d_{xy}$ up/down | MIT / Synthetic Benchmark Generator |
| `pdos_spin/fe.pdos_atm#2(O)_wfc#1(p)` | Tier 2 (Synthetic QE v7.2) | `projwfc.x` | Atom 2 (O), Wfc 1 ($p$-orbital) | 9 columns: $E$, LDOS up/down, $p_z, p_x, p_y$ up/down | MIT / Synthetic Benchmark Generator |
| `pdos_spin_orbit_unsupported.dat` | Tier 3 (Negative Test) | `projwfc.x` | Noncollinear / spin-orbit calculation output containing $j=1.5, j=2.5$ headers | Throws `QENoncollinearUnsupportedError` | MIT / Defensive Test Fixture |

### D. Electronic Band Structures (`bands.x`)

| File | Classification | Target Executable | System / Calculation | Known Reference Ground Truth | License / Source |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `bands.gnu` | Tier 2 (Synthetic QE v7.2) | `bands.x` | Multi-band dispersion table in Gnuplot format | 3 electronic bands<br>31 k-points along path $\Gamma \to X \to M$<br>Separated by double newlines | MIT / Synthetic Benchmark Generator |

### E. Nudged Elastic Band Reaction Pathways (`neb.x`)

| File | Classification | Target Executable | System / Calculation | Known Reference Ground Truth | License / Source |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `neb.dat` | Tier 2 (Synthetic QE v7.2) | `neb.x` | Discrete reaction pathway images | 7 images (index 0 to 6)<br>Reaction coordinate: $0.00$ to $3.60$ Å<br>Peak discrete barrier: $+0.850$ eV at image 3 | MIT / Synthetic Benchmark Generator |
| `neb.int` | Tier 2 (Synthetic QE v7.2) | `neb.x` | Spline-interpolated continuous reaction pathway | 61 interpolated coordinates and energies | MIT / Synthetic Benchmark Generator |
| `neb.path` | Tier 2 (Synthetic QE v7.2) | `neb.x` | Convergence log and barrier summary | Forward barrier: $0.85$ eV<br>Reverse barrier: $1.15$ eV<br>Reaction energy: $-0.30$ eV<br>Climbing image: 3 | MIT / Synthetic Benchmark Generator |

---

## 3. Ground Truth Verification Summary

1. **Energy Values:** In `pw.x`, total energies are stated in Rydberg units and verified against exact conversion ($1\text{ Ry} = 13.605693009\text{ eV}$). In `dos.x`, `projwfc.x`, `bands.x`, and `neb.x`, energies are verified to be native electronvolts (eV) and are strictly shielded from accidental double conversions.
2. **Double Counting Protection:** In `projwfc.x`, total orbital sums are verified against atom-level and component-level projections to confirm that individual sub-orbitals ($d_{z^2}, d_{zx}, d_{zy}, d_{x^2-y^2}, d_{xy}$) sum to the orbital total and do not double-count the LDOS column.
3. **Format Integrity:** In `bands.x`, Gnuplot block delimiters (`\n\n`) are validated for k-point sequence consistency. In `neb.x`, input namelists (`&PATH`, `&INPUT`) are explicitly recognized and rejected to prevent confusing configuration files with output profiles.

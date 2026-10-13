# Atomic Absorption Spectroscopy (AAS) Research Guide & Analytical Reference

**Workspace:** `python_graphing`  
**Domain:** Quantitative Analytical Chemistry, Materials Characterization & Catalyst Stability  
**Standards Compliance:** US EPA Method 7000B (Flame AAS), US EPA Method 7010 (GFAAS), US EPA Method 200.9 (STPF-GFAAS)  
**Author:** Senior Scientific Python Developer & Analytical Chemistry Specialist  

---

## 1. Overview & Theoretical Principles

Atomic Absorption Spectroscopy (AAS) is an elemental analysis technique based on the absorption of optical radiation (light) by free gas-phase ground-state atoms. According to the Beer-Lambert law:

$$T = \frac{I}{I_0} = 10^{-A} \implies A = -\log_{10}(T) = a \cdot b \cdot C$$

where:
- $A$ is the analytical absorbance (dimensionless, in Absorbance Units, AU)
- $I_0$ is the incident light intensity from a line source (Hollow Cathode Lamp, HCL, or Electrodeless Discharge Lamp, EDL)
- $I$ is the transmitted light intensity through the atomic vapor
- $a$ is the atomic absorption coefficient
- $b$ is the optical path length of the atomizer
- $C$ is the free ground-state atomic concentration in the atomization cell.

In `python_graphing`, AAS is integrated alongside Inductively Coupled Plasma (ICP-OES / ICP-MS), X-Ray Diffraction (XRD), Gas Physisorption (BET/BJH), and Electron Microscopy to provide elemental composition, catalyst loading, and electrochemical dissolution kinetics for heterogeneous catalysts and electrocatalysts (HER / OER / ORR).

---

## 2. Supported Atomization Techniques

### 2.1 Flame Atomic Absorption Spectroscopy (FAAS)
- **Atomization Mechanism:** Pneumatic nebulization followed by desolvation, vaporization, and thermal dissociation in an air-acetylene ($T \approx 2300^\circ\text{C}$) or nitrous oxide-acetylene ($T \approx 2800^\circ\text{C}$) laminar flame burner (optical path $b = 5\text{--}10\text{ cm}$).
- **Signal Characteristic:** Steady-state, time-averaged optical absorbance ($A$, AU).
- **Working Dynamic Range:** Typically $0.1\text{ to }10\text{ mg/L}$ ($\text{ppm}$).
- **Regulatory SOP:** **US EPA Method 7000B** (*Flame Atomic Absorption Spectrophotometry*).
- **Physical Limitations:** High sample aspiration rate ($3\text{--}6\text{ mL/min}$); nebulizer transport efficiency is low ($\sim 5\text{--}10\%$); matrix salts can clog nebulizer and burner slots. High absorbance ($A > 0.85$) departs from linearity due to stray light and unabsorbed lamp emission lines.

### 2.2 Graphite Furnace AAS (GFAAS / ETAAS)
- **Atomization Mechanism:** Electrothermal heating of small liquid aliquots ($10\text{--}50\ \mu\text{L}$) inside a longitudinally or transversely heated pyrolytic graphite tube through sequential programmed thermal stages:
  1. *Drying:* Gentle evaporation of solvent ($100\text{--}130^\circ\text{C}$).
  2. *Pyrolysis / Ashing:* Matrix volatilization without analyte loss ($400\text{--}1400^\circ\text{C}$).
  3. *Atomization:* Rapid thermal ramp ($1800\text{--}2600^\circ\text{C}$) creating transient free ground-state atomic vapor.
  4. *Cleanout / Tube Bake:* Desorption of non-volatile residues ($2700^\circ\text{C}$).
- **Signal Characteristic:** Transient, time-resolved pulse. The quantitative response is **Integrated Absorbance** (peak area in absorbance-seconds, $A \cdot s$), adhering to the **Stabilized Temperature Platform Furnace (STPF)** concept.
- **Working Dynamic Range:** $0.1\text{ to }100\ \mu\text{g/L}$ ($\text{ppb}$).
- **Regulatory SOPs:** **US EPA Method 7010** and **US EPA Method 200.9** (*Determination of Trace Elements by Stabilized Temperature Graphite Furnace Atomic Absorption Spectrometry*).
- **Matrix Modifiers:** Chemical matrix modifiers (e.g. $\text{Pd(NO}_3)_2$, $\text{Mg(NO}_3)_2$, $\text{NH}_4\text{H}_2\text{PO}_4$) are commonly required to thermally stabilize volatile analytes during pyrolysis.

### 2.3 Optional Techniques (HG-AAS & CV-AAS)
- **Hydride Generation AAS (HG-AAS):** Reductive derivatization via sodium borohydride ($\text{NaBH}_4$ in acid) converting hydride-forming metalloids (As, Se, Sb, Bi, Sn, Ge) into volatile gaseous hydrides ($\text{AsH}_3$, $\text{H}_2\text{Se}$) swept into a heated quartz T-cell.
- **Cold Vapor AAS (CV-AAS):** Specific reduction of mercury ($\text{Hg}^{2+} \to \text{Hg}^0$) with $\text{SnCl}_2$ or $\text{NaBH}_4$ at ambient room temperature, swept into an unheated quartz optical cell (EPA Method 7470A/7471B).

---

## 3. Data Import & Ingestion Specifications

The `parsers.characterization.aas.parse_aas_file` function supports CSV, TXT, TSV, and Excel (XLSX) instrument exports from major analytical manufacturers (PerkinElmer Syngistix/WinLab32, Agilent SpectrAA, Thermo Scientific SOLAAR, Shimadzu Wizard, Analytik Jena Aspect LS).

### 3.1 Field Schema & Capture
| Standard Field | Internal Dataclass | Units | Description |
| :--- | :--- | :--- | :--- |
| `Sample_ID` | `sample_id` | String | Unique specimen label or alphanumeric run tag |
| `Sample_Type` | `sample_type` | Category | `Standard`, `Sample`, `Blank`, `ICV`, `CCV`, `Spike`, `ReferenceMaterial` |
| `Element` | `element` | String | Analyte elemental symbol (e.g., Fe, Co, Ni, Pt, Ir) |
| `Wavelength_nm`| `wavelength_nm` | nm | Resonant hollow cathode spectral line |
| `Absorbance` | `absorbance` | AU | Steady-state peak height absorbance (FAAS) |
| `Integrated_Absorbance`| `integrated_absorbance` | $A \cdot s$ | Peak area under transient absorbance profile (GFAAS) |
| `Concentration_Standard`| `concentration` | $\text{mg/L}$ or $\mu\text{g/L}$ | Nominal concentration of reference calibration standard |
| `Concentration_Reported`| `concentration_reported` | $\text{mg/L}$ or $\mu\text{g/L}$ | Raw instrument-recalculated value (preserved separately) |
| `Replicate` | `replicate` | Integer | Repetition index (1, 2, 3...) for statistical pooling |
| `Dilution_Factor` | `dilution_factor` | Float | Volumetric pre-dilution ratio ($V_{\text{final}} / V_{\text{aliquot}}$) |
| `Background_Correction` | `background_correction` | String | `D2` (deuterium arc), `Zeeman`, `Smith-Hieftje`, or `None` |

---

## 4. Calibration Curve Analysis

### 4.1 Linear External Calibration Model
For external aqueous standards:

$$A = m \cdot C + b$$

where $m$ is the analytical calibration sensitivity ($\text{AU}/(\text{mg/L})$), $C$ is analyte concentration, and $b$ is the non-zero regression intercept.

In `modules.characterization.aas.fit_aas_calibration`:
1. **Zero-Intercept Control:** Fitting through zero ($A = m \cdot C$) is permitted only when `zero_intercept=True` is explicitly specified and scientifically justified (e.g. after subtraction of certified reagents and solvent blanks).
2. **Weighted Least Squares (WLS):** Heteroscedastic instrumental variance at elevated concentrations can be compensated using empirical weights $w_i = 1/s_i^2$ or $w_i = 1/C_i$.
3. **Range Validation:** Concentrations evaluated above $C_{\max}$ are flagged with `"EXTRAPOLATED_HIGH"` to prevent inaccurate reporting beyond the verified linear dynamic range.
4. **Beer-Lambert Rollover Warning:** High steady-state absorbance ($A > 0.85\text{ AU}$) triggers a rollover warning due to stray radiation and resonance line self-reversal in the flame.

---

## 5. Method of Standard Additions (MOSA)

### 5.1 Matrix Suppression & Rotational Interference
Electrochemical testing electrolytes (e.g., $1.0\ \text{M}\ \text{KOH}$, $0.5\ \text{M}\ \text{H}_2\text{SO}_4$, saturated $\text{KPF}_6$) alter physical aspiration rates, surface tension, droplet size distributions, and ionization equilibrium, causing matrix suppression (rotational interference).

Under standard additions, known increments of standard standard ($\Delta C_i$) are spiked into identical aliquots of the sample matrix:

$$S_i = m \cdot \Delta C_i + S_0$$

where $S_0 = m \cdot C_{\text{apparent}}$ is the unspiked sample signal.

### 5.2 Mathematical Formulation & Extrapolation Uncertainty
The apparent unknown concentration is derived from the horizontal $x$-intercept:

$$C_{\text{apparent}} = \frac{S_0}{m} = \frac{b}{m}$$

Accounting for sample aliquot volume $V_{\text{sample}}$ and final flask volume $V_{\text{final}}$:

$$C_{\text{original}} = C_{\text{apparent}} \times \left(\frac{V_{\text{final}}}{V_{\text{sample}}}\right)$$

The standard error on $C_{\text{apparent}}$ is calculated via the delta method:

$$s_{C} = \frac{s_{y/x}}{m} \sqrt{\frac{1}{N} + \frac{\bar{y}^2}{m^2 \sum (\Delta C_i - \overline{\Delta C})^2}}$$

Confidence interval at $(1 - \alpha)$:

$$\text{CI}_{1-\alpha} = C_{\text{apparent}} \pm t_{\alpha/2, N-2} \cdot s_C$$

> [!WARNING]
> Standard addition compensates exclusively for **matrix-dependent proportional sensitivity changes** (rotational interferences). It **does not** eliminate spectral line overlaps or additive non-specific background absorption. Background absorption must be corrected optically via Deuterium ($D_2$) continuum lamps or Zeeman effect splitting.

---

## 6. Analytical Quality Control (QC) Protocols

Quality control rules are implemented in `modules.characterization.aas.evaluate_aas_qc` based on standard regulatory methods:

| QC Parameter | US EPA Method 7000B (Flame AAS) | US EPA Method 7010 (GFAAS) | Action on Failure |
| :--- | :--- | :--- | :--- |
| **Calibration Linearity** | $R^2 \ge 0.995$ | $R^2 \ge 0.995$ | Recalibrate instrument |
| **Calibration Blank** | $< \text{LOD}$ | $< \text{LOD}$ | Investigate solvent contamination |
| **Initial Cal. Verif. (ICV)** | $90\% \le \text{Recovery} \le 110\%$ | $85\% \le \text{Recovery} \le 115\%$ | Remake second-source stock standard |
| **Continuing Cal. Verif. (CCV)**| $90\% \le \text{Recovery} \le 110\%$ | $85\% \le \text{Recovery} \le 115\%$ | Re-run bracketed sample block |
| **Replicate Precision (%RSD)** | $\% \text{RSD} \le 10\%$ (samples $\ge \text{LOQ}$) | $\% \text{RSD} \le 15\%$ (samples $\ge \text{LOQ}$) | Check nebulizer uptake / furnace injection |
| **Matrix Spike Recovery** | $80\% \le \text{Recovery} \le 120\%$ | $70\% \le \text{Recovery} \le 130\%$ | Switch to Method of Standard Additions |

---

## 7. Electrocatalyst Loading & Stability Leaching

### 7.1 Bulk Catalyst Metal Loading (wt.%)
Quantification of target metal loading (e.g. Pt, Ir, Ru, Ni, Fe, Co) from acid digestion:

$$\text{Metal Loading (wt.\%)} = \frac{(C_{\text{sample}} - C_{\text{blank}}) \times (V_{\text{digest}} / 1000) \times \text{DF}}{m_{\text{catalyst}}} \times 100\%$$

where:
- $C$ is in $\text{mg/L}$
- $V_{\text{digest}}$ is in $\text{mL}$
- $m_{\text{catalyst}}$ is initial dry catalyst mass in $\text{mg}$
- $\text{DF}$ is volumetric dilution factor.

### 7.2 Electrochemical Leaching with Aliquot Compensation
During continuous chronopotentiometric or chronoamperometric stress tests (e.g., 100h OER in $1.0\ \text{M}\ \text{KOH}$), periodic withdrawal of analytical aliquots removes dissolved metal from the cell inventory. Total cumulative dissolved metal at time-step $k$ is:

$$m_{\text{dissolved, total}}(t_k) = C_k \cdot V_k + \sum_{j=1}^{k-1} C_j \cdot V_{\text{aliquot}, j}$$

$$\text{Percentage Leached } \% = \frac{m_{\text{dissolved, total}}}{m_{\text{initial\_metal\_inventory}}} \times 100\%$$

where $m_{\text{initial\_metal\_inventory}} = m_{\text{catalyst}} \times w_{\text{metal, actual}}$.

---

## 8. CLI Runner Commands

The unified workspace CLI supports direct AAS operations:

```bash
# 1. Plot FAAS external calibration curve
python -m scripts.main_plots --type aas --input data/examples/characterization/aas/faas_calibration_fe.csv

# 2. Method of Standard Additions
python -m scripts.main_plots --type aas --input data/examples/characterization/aas/aas_standard_addition_catalyst.csv --calibration-method standard_addition

# 3. Batch generate all 21 publication figures across characterization (including 8 AAS figures)
python -m scripts.main_plots --category characterization
```

---

## 9. References & Regulatory Standards
1. **US EPA Method 7000B** (2007), *Flame Atomic Absorption Spectrophotometry*, Test Methods for Evaluating Solid Waste, Physical/Chemical Methods (SW-846).
2. **US EPA Method 7010** (2007), *Graphite Furnace Atomic Absorption Spectrophotometry*, SW-846.
3. **US EPA Method 200.9** (1994), *Determination of Trace Elements by Stabilized Temperature Graphite Furnace Atomic Absorption Spectrometry*, Rev. 2.2.
4. **Welz, B. & Sperling, M.** (1999), *Atomic Absorption Spectrometry*, 3rd ed., Wiley-VCH.
5. **Slavin, W., Manning, D.C., Carnrick, G.R.** (1981), *The Stabilized Temperature Platform Furnace*, *Atomic Spectroscopy* 2, 137–145.

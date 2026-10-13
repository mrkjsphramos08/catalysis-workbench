"""
Experimental Electrochemistry Analysis and Visualization Package.

Comprehensive processing, kinetic modeling, and publication-quality graphing for:
- Metrohm Autolab (NOVA) and BioLogic (EC-Lab) instruments
- Reference electrode calibrations (SHE, RHE, SCE, Ag/AgCl, Hg/HgO) and iR compensation
- Cyclic Voltammetry (CV) peak analysis and charge integration
- Linear Sweep Voltammetry (LSV) overpotentials (eta_10, eta_50, eta_100)
- Tafel slope regression, exchange current densities, and steady-state validation
- Double-layer capacitance (C_dl), ECSA, and roughness factors
- Electrochemical Impedance Spectroscopy (EIS), Randles fitting, and Lin-KK validation
- Chronoamperometry (CA) and Chronopotentiometry (CP) long-term durability metrics
- Cross-instrument benchmarking and residual discrepancy auditing.
"""

from modules.electrochemistry.benchmarking import (
    compare_instrument_lsv,
    plot_cross_instrument_benchmark,
)
from modules.electrochemistry.chronoamperometry import (
    calculate_ca_stability_metrics,
    plot_chronoamperometry,
)
from modules.electrochemistry.chronopotentiometry import (
    calculate_cp_stability_metrics,
    plot_chronopotentiometry,
)
from modules.electrochemistry.cv import (
    calculate_voltammetric_charge,
    extract_redox_peaks,
    plot_cv,
    plot_cv_scan_rate_series,
)
from modules.electrochemistry.ecsa import (
    CdlInterpretationWarning,
    calculate_double_layer_capacitance,
    estimate_ecsa,
    plot_cdl_analysis,
    plot_ecsa_analysis,
)
from modules.electrochemistry.eis import (
    calculate_impedance_polar,
    fit_equivalent_circuit,
    plot_bode,
    plot_eis_comprehensive,
    plot_nyquist,
    run_kramers_kronig_test,
)
from modules.electrochemistry.lsv import (
    extract_overpotential_at_current,
    plot_lsv,
    plot_lsv_comparison,
)
from modules.electrochemistry.reference_electrodes import (
    REFERENCE_ELECTRODES_298K,
    apply_ir_compensation,
    calculate_overpotential,
    convert_potential_to_rhe,
    get_standard_electrode_potential,
)
from modules.electrochemistry.tafel import (
    NonSteadyStateWarning,
    calculate_tafel_slope,
    plot_tafel,
    plot_tafel_comparison,
)

__all__ = [
    # Reference electrodes & conversions
    "REFERENCE_ELECTRODES_298K",
    "convert_potential_to_rhe",
    "apply_ir_compensation",
    "calculate_overpotential",
    "get_standard_electrode_potential",
    # Voltammetry
    "extract_redox_peaks",
    "calculate_voltammetric_charge",
    "plot_cv",
    "plot_cv_scan_rate_series",
    "extract_overpotential_at_current",
    "plot_lsv",
    "plot_lsv_comparison",
    # Kinetics & Surface Area
    "calculate_tafel_slope",
    "plot_tafel",
    "plot_tafel_comparison",
    "NonSteadyStateWarning",
    "calculate_double_layer_capacitance",
    "estimate_ecsa",
    "plot_ecsa_analysis",
    "CdlInterpretationWarning",
    # EIS
    "calculate_impedance_polar",
    "fit_equivalent_circuit",
    "run_kramers_kronig_test",
    "plot_nyquist",
    "plot_bode",
    "plot_eis_comprehensive",
    # Durability
    "calculate_ca_stability_metrics",
    "plot_chronoamperometry",
    "calculate_cp_stability_metrics",
    "plot_chronopotentiometry",
    # Cross-instrument benchmarking
    "compare_instrument_lsv",
    "plot_cross_instrument_benchmark",
]

"""
Materials Characterization Quantitative Analysis and Visualization Package.

Comprehensive processing and journal-ready figures for:
- XRD: Bragg d-spacing, Scherrer crystallite size, deconvolution, and stick patterns
- FTIR: Beer-Lambert absorbance, baseline subtraction, band integration, descending wavenumber plots
- Raman: Cosmic-ray despiking, baseline subtraction, multi-peak deconvolution, height vs area ratios
- SEM / TEM / HRTEM: Calibrated scale bars, particle segmentation, size distributions, FFT lattice d-spacings
- ICP-OES / ICP-MS: Calibration regression, LOD/LOQ, blank corrections, catalyst metal leaching
- BET / BJH: Specific surface area, Rouquerol criteria, pore size distributions (dV/dD)
- Multi-Technique: Cross-correlation between physical characterization and catalytic/electrochemical performance.
"""

from modules.characterization.common_peaks import (
    detect_peaks,
    estimate_baseline,
    subtract_baseline,
    remove_cosmic_rays,
    smooth_spectrum,
    fit_peaks,
    integrate_band,
)
from modules.characterization.xrd import (
    calculate_d_spacing,
    calculate_crystallite_size_scherrer,
    analyze_xrd_peaks,
    plot_xrd_pattern,
    plot_xrd_stacked,
    plot_xrd_fit,
    ScherrerBroadeningWarning,
)
from modules.characterization.ftir import (
    COMMON_IR_BANDS,
    analyze_ftir_bands,
    plot_ftir_spectrum,
    plot_ftir_stacked,
)
from modules.characterization.raman import (
    process_raman_spectrum,
    calculate_peak_ratios,
    plot_raman_spectrum,
    plot_raman_stacked,
    plot_raman_deconvolution,
)
from modules.characterization.microscopy import (
    add_scale_bar,
    segment_particles,
    analyze_hrtem_lattice_fringes,
    plot_microscopy_with_scalebar,
    plot_segmentation_overlay,
    plot_particle_size_distribution,
    plot_hrtem_fft_panel,
)
from modules.characterization.icp import (
    fit_icp_calibration,
    calculate_catalyst_leaching,
    plot_icp_calibration,
    plot_icp_leaching_comparison,
)
from modules.characterization.bet import (
    calculate_bet_surface_area,
    calculate_bjh_pore_distribution,
    plot_isotherm,
    plot_bet_linear,
    plot_bjh_pore_distribution,
)
from modules.characterization.comparative import (
    correlate_cross_technique_metrics,
    plot_cross_technique_correlation,
)
from modules.characterization.elemental_analysis import (
    convert_concentration_units,
    compute_replicate_statistics,
    linear_calibration_fit,
    calculate_detection_limits,
    calculate_catalyst_leaching_common,
)
from modules.characterization.aas import (
    fit_aas_calibration,
    fit_standard_addition,
    quantify_aas_samples,
    evaluate_aas_qc,
    calculate_aas_catalyst_loading,
    calculate_aas_leaching,
    plot_aas_calibration,
    plot_standard_addition,
    plot_aas_concentration_comparison,
    plot_aas_catalyst_loading,
    plot_aas_leaching,
    plot_aas_stability_before_after,
    plot_aas_replicate_distribution,
    plot_aas_qc_chart,
)

__all__ = [
    # Common Peak Fitting & Signal Processing
    "detect_peaks",
    "estimate_baseline",
    "subtract_baseline",
    "remove_cosmic_rays",
    "smooth_spectrum",
    "fit_peaks",
    "integrate_band",
    # XRD
    "calculate_d_spacing",
    "calculate_crystallite_size_scherrer",
    "analyze_xrd_peaks",
    "plot_xrd_pattern",
    "plot_xrd_stacked",
    "plot_xrd_fit",
    "ScherrerBroadeningWarning",
    # FTIR
    "COMMON_IR_BANDS",
    "analyze_ftir_bands",
    "plot_ftir_spectrum",
    "plot_ftir_stacked",
    # Raman
    "process_raman_spectrum",
    "calculate_peak_ratios",
    "plot_raman_spectrum",
    "plot_raman_stacked",
    "plot_raman_deconvolution",
    # Microscopy
    "add_scale_bar",
    "segment_particles",
    "analyze_hrtem_lattice_fringes",
    "plot_microscopy_with_scalebar",
    "plot_segmentation_overlay",
    "plot_particle_size_distribution",
    "plot_hrtem_fft_panel",
    # Elemental & ICP
    "convert_concentration_units",
    "compute_replicate_statistics",
    "linear_calibration_fit",
    "calculate_detection_limits",
    "calculate_catalyst_leaching_common",
    "fit_icp_calibration",
    "calculate_catalyst_leaching",
    "plot_icp_calibration",
    "plot_icp_leaching_comparison",
    # AAS
    "fit_aas_calibration",
    "fit_standard_addition",
    "quantify_aas_samples",
    "evaluate_aas_qc",
    "calculate_aas_catalyst_loading",
    "calculate_aas_leaching",
    "plot_aas_calibration",
    "plot_standard_addition",
    "plot_aas_concentration_comparison",
    "plot_aas_catalyst_loading",
    "plot_aas_leaching",
    "plot_aas_stability_before_after",
    "plot_aas_replicate_distribution",
    "plot_aas_qc_chart",
    # BET / BJH
    "calculate_bet_surface_area",
    "calculate_bjh_pore_distribution",
    "plot_isotherm",
    "plot_bet_linear",
    "plot_bjh_pore_distribution",
    # Comparative
    "correlate_cross_technique_metrics",
    "plot_cross_technique_correlation",
]

"""
Command-Line Runner and Demonstration Plot Generator for Publication Figures.

Supports batch generation of benchmark figures:
  python -m scripts.main_plots --all
  python -m scripts.main_plots --category dft
  python -m scripts.main_plots --category catalysis
  python -m scripts.main_plots --category chemical_engineering

And direct processing and visualization of actual Quantum ESPRESSO results:
  python -m scripts.main_plots --type dos --input PATH [--fermi EF]
  python -m scripts.main_plots --type pdos --input DIRECTORY [--elements Pt O] [--orbitals d p]
  python -m scripts.main_plots --type bands --input PATH [--fermi EF]
  python -m scripts.main_plots --type neb --input PATH

Exports vector (PDF, SVG) and high-resolution raster (PNG, TIFF) figures
and records scientific provenance metadata JSON alongside figure outputs.
"""

import argparse
import sys
from pathlib import Path

# Ensure project root is on sys.path for direct script invocation
PROJECT_ROOT_PATH = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT_PATH) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT_PATH))

from typing import Any, Dict, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.config import (
    DEFAULT_DPI,
    FIGURES_CATALYSIS_DIR,
    FIGURES_CHEMENG_DIR,
    FIGURES_DFT_DIR,
    FIGURES_ELECTROCHEM_DIR,
    FIGURES_CHAR_DIR,
    FIGURES_XRD_DIR,
    FIGURES_FTIR_DIR,
    FIGURES_RAMAN_DIR,
    FIGURES_MICROSCOPY_DIR,
    FIGURES_ICP_DIR,
    FIGURES_AAS_DIR,
    FIGURES_BET_DIR,
    FIGURES_COMPARATIVE_DIR,
    PROJECT_ROOT,
    OUTPUTS_RASTER_DIR,
    OUTPUTS_VECTOR_DIR,
    EXAMPLES_DATA_DIR,
)
from scripts.utils import (
    apply_publication_style,
    create_publication_figure,
    load_csv_data,
    save_publication_figure,
)
from modules.dft.electronic_structure import (
    plot_dos,
    plot_pdos,
    plot_band_structure,
)
from modules.dft.adsorption import (
    plot_adsorption_energies,
    plot_adsorption_scaling,
)
from modules.catalysis.reaction_energy import (
    plot_reaction_pathway,
)
from modules.catalysis.volcano import (
    plot_volcano,
)
from modules.chemical_engineering.kinetics import (
    plot_arrhenius,
    plot_conversion_selectivity,
    plot_reactor_profiles,
    plot_parity_and_residuals,
)
from parsers.quantum_espresso import (
    QEPdosCollection,
    parse_bands_gnu,
    parse_dos_file,
    parse_neb_output,
    parse_pw_output,
    plot_qe_bands,
    plot_qe_dos,
    plot_qe_neb_profile,
    plot_qe_pdos,
    record_figure_provenance,
)
from modules.electrochemistry import (
    fit_equivalent_circuit,
    plot_bode,
    plot_cdl_analysis,
    plot_chronoamperometry,
    plot_chronopotentiometry,
    plot_cross_instrument_benchmark,
    plot_cv,
    plot_eis_comprehensive,
    plot_lsv,
    plot_nyquist,
    plot_tafel,
)
from parsers.electrochemistry import load_electrochemistry_data
from parsers.characterization import (
    parse_xrd_file,
    parse_ftir_file,
    convert_transmittance_to_absorbance,
    parse_raman_file,
    parse_microscopy_image,
    parse_icp_file,
    parse_aas_file,
    parse_isotherm_file,
    record_characterization_provenance,
)
from modules.characterization import (
    calculate_d_spacing,
    calculate_crystallite_size_scherrer,
    analyze_xrd_peaks,
    plot_xrd_pattern,
    plot_xrd_stacked,
    plot_xrd_fit,
    analyze_ftir_bands,
    plot_ftir_spectrum,
    plot_ftir_stacked,
    process_raman_spectrum,
    calculate_peak_ratios,
    plot_raman_spectrum,
    plot_raman_stacked,
    plot_raman_deconvolution,
    segment_particles,
    analyze_hrtem_lattice_fringes,
    plot_microscopy_with_scalebar,
    plot_segmentation_overlay,
    plot_particle_size_distribution,
    plot_hrtem_fft_panel,
    fit_icp_calibration,
    calculate_catalyst_leaching,
    plot_icp_calibration,
    plot_icp_leaching_comparison,
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
    calculate_bet_surface_area,
    calculate_bjh_pore_distribution,
    plot_isotherm,
    plot_bet_linear,
    plot_bjh_pore_distribution,
    correlate_cross_technique_metrics,
    plot_cross_technique_correlation,
    fit_peaks,
)


def generate_dft_figures(formats: Sequence[str] = ("pdf", "svg", "png", "tiff"), dpi: int = DEFAULT_DPI) -> List[str]:
    """Generate DFT electronic structure and adsorption benchmark figures."""
    generated = []
    FIGURES_DFT_DIR.mkdir(parents=True, exist_ok=True)
    apply_publication_style()

    # 1. Total and Projected Density of States
    dos_csv = EXAMPLES_DATA_DIR / "dft_dos_synthetic.csv"
    if dos_csv.exists():
        df_dos = load_csv_data(dos_csv)
        fig, (ax_tdos, ax_pdos) = create_publication_figure(
            columns="single",
            height_mm=110.0,
            nrows=2,
            ncols=1,
            sharex=True,
            gridspec_kw={"hspace": 0.08},
        )
        plot_dos(
            energies=df_dos["energy"],
            dos_up=df_dos["dos_spin_up"],
            dos_down=df_dos["dos_spin_down"],
            e_fermi=0.0,
            ax=ax_tdos,
            energy_range=(-7.0, 5.0),
            label_up=r"Total $\uparrow$",
            label_down=r"Total $\downarrow$",
        )
        ax_tdos.set_xlabel("")

        pdos_projections = {
            r"Pt $5d$": {"up": df_dos["pt_5d_up"], "down": df_dos["pt_5d_down"]},
            r"O $2p$": {"up": df_dos["o_2p_up"], "down": df_dos["o_2p_down"]},
        }
        plot_pdos(
            energies=df_dos["energy"],
            projections=pdos_projections,
            total_dos=df_dos["total_dos"],
            e_fermi=0.0,
            ax=ax_pdos,
            energy_range=(-7.0, 5.0),
        )

        name = "fig1_dft_dos_pdos"
        save_publication_figure(fig, name, output_dirs=[FIGURES_DFT_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [dos_csv], [FIGURES_DFT_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"figure_type": "benchmark_dos"})
        plt.close(fig)
        generated.append(name)

    # 2. Electronic Band Structure
    bands_csv = EXAMPLES_DATA_DIR / "dft_bands_synthetic.csv"
    if bands_csv.exists():
        df_bands = load_csv_data(bands_csv)
        band_cols = [c for c in df_bands.columns if c.startswith("band_")]
        bands_data = df_bands[band_cols].values.T

        hs_points = [(0.0, r"$\Gamma$"), (0.8, "$X$"), (1.6, "$M$"), (2.4, r"$\Gamma$")]
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        plot_band_structure(
            k_distances=df_bands["k_distance"],
            bands=bands_data,
            high_symmetry_points=hs_points,
            e_fermi=0.0,
            ax=ax,
            energy_range=(-5.0, 4.0),
        )
        name = "fig2_dft_band_structure"
        save_publication_figure(fig, name, output_dirs=[FIGURES_DFT_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [bands_csv], [FIGURES_DFT_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"figure_type": "benchmark_bands"})
        plt.close(fig)
        generated.append(name)

    # 3. Adsorption Energy Site and Facet Comparison
    ads_csv = EXAMPLES_DATA_DIR / "adsorption_energies_synthetic.csv"
    if ads_csv.exists():
        df_ads = load_csv_data(ads_csv)
        fig, ax = create_publication_figure(columns="single", aspect_ratio="1.3")
        plot_adsorption_energies(
            data=df_ads,
            surface_col="surface",
            adsorption_col="e_ads",
            site_col="site",
            error_col="error",
            ax=ax,
            reference_note="Ref: CO(g), Exothermic < 0",
        )
        name = "fig3_adsorption_energy_comparison"
        save_publication_figure(fig, name, output_dirs=[FIGURES_DFT_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [ads_csv], [FIGURES_DFT_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"figure_type": "benchmark_adsorption"})
        plt.close(fig)
        generated.append(name)

    return generated


def generate_catalysis_figures(formats: Sequence[str] = ("pdf", "svg", "png", "tiff"), dpi: int = DEFAULT_DPI) -> List[str]:
    """Generate catalysis reaction profiles and volcano benchmark figures."""
    generated = []
    FIGURES_CATALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    apply_publication_style()

    # 4. Catalytic Reaction Free Energy Profile
    rxn_csv = EXAMPLES_DATA_DIR / "catalysis_reaction_profile_synthetic.csv"
    if rxn_csv.exists():
        df_rxn = load_csv_data(rxn_csv)
        steps = df_rxn["step_name"].tolist()
        is_ts = df_rxn["is_transition_state"].tolist()
        catalysts = [c for c in df_rxn.columns if c not in ("step_name", "is_transition_state")]

        pathways = {c: df_rxn[c].tolist() for c in catalysts}
        ts_dict = {c: is_ts for c in catalysts}

        fig, ax = create_publication_figure(columns="double", aspect_ratio="2.2")
        plot_reaction_pathway(
            pathways=pathways,
            step_labels=steps,
            transition_states=ts_dict,
            relative_to_initial=True,
            ax=ax,
            title="CO Oxidation Reaction Free Energy Diagram (298.15 K, 1 bar)",
        )
        name = "fig4_catalysis_reaction_profile"
        save_publication_figure(fig, name, output_dirs=[FIGURES_CATALYSIS_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [rxn_csv], [FIGURES_CATALYSIS_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"figure_type": "benchmark_reaction_pathway"})
        plt.close(fig)
        generated.append(name)

    # 5. Volcano Plot (ORR descriptor model)
    volc_csv = EXAMPLES_DATA_DIR / "catalysis_volcano_synthetic.csv"
    if volc_csv.exists():
        df_volc = load_csv_data(volc_csv)
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        model_params = (1.0, 0.0, -1.0, 1.60)
        highlights = df_volc[df_volc["candidate"]]["catalyst"].tolist()

        plot_volcano(
            descriptors=df_volc["descriptor_g_oh"],
            activities=df_volc["activity_potential"],
            catalyst_names=df_volc["catalyst"].tolist(),
            theoretical_model=model_params,
            highlight_catalysts=highlights,
            ax=ax,
            x_label=r"Descriptor $\Delta G_{\mathrm{OH}^*}$ (eV)",
            y_label=r"Limiting Potential $U_{\mathrm{L}}$ (V vs RHE)",
            title="Oxygen Reduction Reaction (ORR) Sabatier Volcano",
        )
        name = "fig5_catalysis_volcano"
        save_publication_figure(fig, name, output_dirs=[FIGURES_CATALYSIS_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [volc_csv], [FIGURES_CATALYSIS_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"figure_type": "benchmark_volcano_orr"})
        plt.close(fig)
        generated.append(name)

    return generated


def generate_chemeng_figures(formats: Sequence[str] = ("pdf", "svg", "png", "tiff"), dpi: int = DEFAULT_DPI) -> List[str]:
    """Generate chemical engineering kinetics and reactor figures."""
    generated = []
    FIGURES_CHEMENG_DIR.mkdir(parents=True, exist_ok=True)
    apply_publication_style()

    # 6. Arrhenius Regression and Residuals
    arrh_csv = EXAMPLES_DATA_DIR / "kinetics_arrhenius_synthetic.csv"
    if arrh_csv.exists():
        df_arrh = load_csv_data(arrh_csv)
        fig, res_stats = plot_arrhenius(
            temperatures_k=df_arrh["temperature_K"],
            rate_constants=df_arrh["rate_constant_k"],
            rate_constant_units=r"\mathrm{s}^{-1}",
            show_residuals=True,
            title="Arrhenius Temperature Dependence",
        )
        name = "fig6_kinetics_arrhenius"
        save_publication_figure(fig, name, output_dirs=[FIGURES_CHEMENG_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [arrh_csv], [FIGURES_CHEMENG_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"figure_type": "benchmark_arrhenius"})
        plt.close(fig)
        generated.append(name)

    # 7. Reactor Conversion, Selectivity, and Parity
    reactor_csv = EXAMPLES_DATA_DIR / "reactor_performance_synthetic.csv"
    if reactor_csv.exists():
        df_reactor = load_csv_data(reactor_csv)
        fig, (ax_t, ax_cs) = create_publication_figure(
            columns="double",
            aspect_ratio="2.2",
            nrows=1,
            ncols=2,
        )
        plot_reactor_profiles(
            temperature_c=df_reactor["temperature_C"],
            conversion=df_reactor["conversion_percent"],
            selectivity=df_reactor["selectivity_percent"],
            yield_val=df_reactor["yield_percent"],
            ax=ax_t,
            title="Reactor Performance Profile",
        )
        plot_conversion_selectivity(
            conversion=df_reactor["conversion_percent"],
            selectivity=df_reactor["selectivity_percent"],
            catalyst_label="Pt/Al2O3 Catalyst",
            ax=ax_cs,
        )
        name = "fig7_chemical_engineering_reactor"
        save_publication_figure(fig, name, output_dirs=[FIGURES_CHEMENG_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [reactor_csv], [FIGURES_CHEMENG_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"figure_type": "benchmark_reactor"})
        plt.close(fig)
        generated.append(name)

    return generated


def generate_electrochemistry_figures(
    formats: Sequence[str] = ("pdf", "svg", "png", "tiff"), dpi: int = DEFAULT_DPI
) -> List[str]:
    """Generate experimental electrochemistry publication benchmark figures."""
    generated = []
    FIGURES_ELECTROCHEM_DIR.mkdir(parents=True, exist_ok=True)
    apply_publication_style()

    examples_ec_dir = EXAMPLES_DATA_DIR / "electrochemistry"
    metrohm_dir = examples_ec_dir / "metrohm"
    biologic_dir = examples_ec_dir / "biologic"

    # 1. Metrohm OER LSV Polarization Curve
    metrohm_lsv_path = metrohm_dir / "metrohm_oer_lsv.txt"
    if metrohm_lsv_path.exists():
        data_oer_m = load_electrochemistry_data(metrohm_lsv_path, temperature_k=298.15)
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        plot_lsv(data_oer_m, reaction="OER", target_currents=[10.0, 50.0], ax=ax, title="Metrohm Autolab: Alkaline OER (NiFe-LDH)")
        name = "fig1_metrohm_oer_lsv"
        save_publication_figure(fig, name, output_dirs=[FIGURES_ELECTROCHEM_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [metrohm_lsv_path], [FIGURES_ELECTROCHEM_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"synthetic_temperature_k": 298.15, "instrument": "metrohm", "technique": "LSV", "reaction": "OER"})
        plt.close(fig)
        generated.append(name)

    # 2. BioLogic HER LSV Polarization Curve
    bio_her_path = biologic_dir / "biologic_her_lsv.mpt"
    if bio_her_path.exists():
        data_her_b = load_electrochemistry_data(bio_her_path, temperature_k=298.15, reference_electrode="Ag/AgCl_sat")
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        plot_lsv(data_her_b, reaction="HER", target_currents=[10.0], ax=ax, title="BioLogic EC-Lab: Acidic HER (Pt/C in 0.5 M H2SO4)")
        name = "fig2_biologic_her_lsv"
        save_publication_figure(fig, name, output_dirs=[FIGURES_ELECTROCHEM_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [bio_her_path], [FIGURES_ELECTROCHEM_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"synthetic_temperature_k": 298.15, "synthetic_reference_assumption": "Ag/AgCl_sat", "instrument": "biologic", "technique": "LSV", "reaction": "HER"})
        plt.close(fig)
        generated.append(name)

    # 3. Metrohm Multi-Cycle Cyclic Voltammetry
    metrohm_cv_path = metrohm_dir / "metrohm_cv_multicycle.csv"
    if metrohm_cv_path.exists():
        data_cv_m = load_electrochemistry_data(metrohm_cv_path, temperature_k=298.15)
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        plot_cv(data_cv_m, ax=ax, mark_peaks=True, title="Metrohm Autolab: Cyclic Voltammetry (3 Cycles)")
        name = "fig3_metrohm_cv_multicycle"
        save_publication_figure(fig, name, output_dirs=[FIGURES_ELECTROCHEM_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [metrohm_cv_path], [FIGURES_ELECTROCHEM_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"synthetic_temperature_k": 298.15, "instrument": "metrohm", "technique": "CV"})
        plt.close(fig)
        generated.append(name)

    # 4. Double-Layer Capacitance (Cdl) & ECSA Determination
    if metrohm_cv_path.exists():
        rates = [0.020, 0.040, 0.060, 0.080, 0.100]
        scan_datasets = []
        for r in rates:
            d_copy = load_electrochemistry_data(metrohm_cv_path, temperature_k=298.15)
            d_copy.scan_rate_V_s = r
            d_copy.current_A = d_copy.current_A * (r / 0.050)
            scan_datasets.append(d_copy)

        fig, (ax_cv, ax_fit), res_ecsa = plot_cdl_analysis(
            scan_datasets,
            potential_sample_V=0.20,
            area_cm2=0.196,
            title="Double-Layer Capacitance & ECSA Analysis",
        )
        name = "fig4_electrochem_ecsa_analysis"
        save_publication_figure(fig, name, output_dirs=[FIGURES_ELECTROCHEM_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [metrohm_cv_path], [FIGURES_ELECTROCHEM_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"synthetic_temperature_k": 298.15, "technique": "ECSA", "cdl_uF_cm2": res_ecsa["cdl_mF_cm2"] * 1000.0})
        plt.close(fig)
        generated.append(name)

    # 5. Metrohm EIS Nyquist Plot with Equivalent Circuit Model
    metrohm_eis_path = metrohm_dir / "metrohm_fra_eis.txt"
    if metrohm_eis_path.exists():
        data_eis_m = load_electrochemistry_data(metrohm_eis_path, temperature_k=298.15)
        fit_m = fit_equivalent_circuit(data_eis_m.frequency_Hz, data_eis_m.impedance_real_ohm, data_eis_m.impedance_imag_ohm, circuit_string="R0-p(R1,CPE1)")
        fig, ax = plot_nyquist(data_eis_m, fit_results=fit_m, show_frequencies=True, title="Metrohm Autolab: Nyquist Plot (-Z'' vs Z')")
        name = "fig5_metrohm_eis_nyquist"
        save_publication_figure(fig, name, output_dirs=[FIGURES_ELECTROCHEM_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [metrohm_eis_path], [FIGURES_ELECTROCHEM_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"synthetic_temperature_k": 298.15, "instrument": "metrohm", "technique": "EIS", "rs_ohm": fit_m["rs_ohm"], "rct_ohm": fit_m["rct_ohm"]})
        plt.close(fig)
        generated.append(name)

    # 6. BioLogic EIS Bode Plot (Magnitude & Phase)
    bio_eis_path = biologic_dir / "biologic_peis.mpt"
    if bio_eis_path.exists():
        data_eis_b = load_electrochemistry_data(bio_eis_path, temperature_k=298.15)
        fit_b = fit_equivalent_circuit(data_eis_b.frequency_Hz, data_eis_b.impedance_real_ohm, data_eis_b.impedance_imag_ohm, circuit_string="R0-p(R1,CPE1)")
        fig, axes = plot_bode(data_eis_b, fit_results=fit_b, title="BioLogic EC-Lab: Bode Magnitude & Phase")
        name = "fig6_biologic_eis_bode"
        save_publication_figure(fig, name, output_dirs=[FIGURES_ELECTROCHEM_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [bio_eis_path], [FIGURES_ELECTROCHEM_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"synthetic_temperature_k": 298.15, "instrument": "biologic", "technique": "EIS"})
        plt.close(fig)
        generated.append(name)

    # 7. Comprehensive 3-Panel EIS Analysis (Nyquist + Bode Mag + Bode Phase)
    if metrohm_eis_path.exists():
        data_eis_m = load_electrochemistry_data(metrohm_eis_path, temperature_k=298.15)
        fit_m = fit_equivalent_circuit(data_eis_m.frequency_Hz, data_eis_m.impedance_real_ohm, data_eis_m.impedance_imag_ohm, circuit_string="R0-p(R1,CPE1)")
        fig = plot_eis_comprehensive(data_eis_m, fit_result=fit_m, title="Comprehensive EIS: Modified Randles Equivalent Circuit Fit")
        name = "fig7_eis_randles_fit_comprehensive"
        save_publication_figure(fig, name, output_dirs=[FIGURES_ELECTROCHEM_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [metrohm_eis_path], [FIGURES_ELECTROCHEM_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"synthetic_temperature_k": 298.15, "technique": "EIS", "circuit": "R0-p(R1,CPE1)"})
        plt.close(fig)
        generated.append(name)

    # 8. Quasi-Steady-State Tafel Kinetics
    if metrohm_lsv_path.exists():
        data_oer_m = load_electrochemistry_data(metrohm_lsv_path, temperature_k=298.15)
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        plot_tafel(data_oer_m, reaction="OER", fit_range=(1.0, 10.0), range_type="current", ax=ax, title="Quasi-Steady-State OER Tafel Kinetics")
        name = "fig8_electrochem_tafel_kinetics"
        save_publication_figure(fig, name, output_dirs=[FIGURES_ELECTROCHEM_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [metrohm_lsv_path], [FIGURES_ELECTROCHEM_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"synthetic_temperature_k": 298.15, "technique": "Tafel", "reaction": "OER"})
        plt.close(fig)
        generated.append(name)

    # 9. Metrohm Chronoamperometry Durability
    metrohm_ca_path = metrohm_dir / "metrohm_stability_ca.txt"
    if metrohm_ca_path.exists():
        data_ca_m = load_electrochemistry_data(metrohm_ca_path, temperature_k=298.15)
        fig, ax = plot_chronoamperometry(data_ca_m, show_retention=True, potential_label="E = 1.60 V vs RHE", title="Metrohm Autolab: Chronoamperometry (1 h)")
        name = "fig9_metrohm_stability_ca"
        save_publication_figure(fig, name, output_dirs=[FIGURES_ELECTROCHEM_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [metrohm_ca_path], [FIGURES_ELECTROCHEM_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"synthetic_temperature_k": 298.15, "instrument": "metrohm", "technique": "CA"})
        plt.close(fig)
        generated.append(name)

    # 10. Cross-Instrument Benchmarking (Metrohm Autolab vs BioLogic EC-Lab)
    bio_oer_path = biologic_dir / "biologic_oer_lsv.mpt"
    if metrohm_lsv_path.exists() and bio_oer_path.exists():
        data_m = load_electrochemistry_data(metrohm_lsv_path, temperature_k=298.15)
        data_b = load_electrochemistry_data(bio_oer_path, temperature_k=298.15)
        fig, axes = plot_cross_instrument_benchmark(data_m, data_b, reaction="OER", target_j=10.0, title="Cross-Instrument Benchmarking: Metrohm Autolab vs BioLogic EC-Lab")
        name = "fig10_cross_instrument_benchmark"
        save_publication_figure(fig, name, output_dirs=[FIGURES_ELECTROCHEM_DIR], formats=formats, dpi=dpi)
        record_figure_provenance(name, [metrohm_lsv_path, bio_oer_path], [FIGURES_ELECTROCHEM_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"synthetic_temperature_k": 298.15, "figure_type": "cross_instrument_benchmarking"})
        plt.close(fig)
        generated.append(name)

    return generated


def generate_characterization_figures(
    formats: Sequence[str] = ("pdf", "svg", "png", "tiff"),
    dpi: int = DEFAULT_DPI,
) -> List[str]:
    """Generate all 13 publication figures across materials characterization techniques."""
    generated = []
    apply_publication_style()
    char_base = EXAMPLES_DATA_DIR / "characterization"

    for d in [
        FIGURES_CHAR_DIR,
        FIGURES_XRD_DIR,
        FIGURES_FTIR_DIR,
        FIGURES_RAMAN_DIR,
        FIGURES_MICROSCOPY_DIR,
        FIGURES_ICP_DIR,
        FIGURES_BET_DIR,
        FIGURES_COMPARATIVE_DIR,
    ]:
        d.mkdir(parents=True, exist_ok=True)

    # 1. XRD Single Pattern
    xrd_pristine = char_base / "xrd" / "catalyst_pristine_xrd.csv"
    if xrd_pristine.exists():
        data_xrd = parse_xrd_file(xrd_pristine)
        fig, ax = plot_xrd_pattern(data_xrd, annotate_peaks=True, show_d_spacing=True, title="Pt/C Catalyst: High-Resolution Cu Kα Powder XRD")
        name = "fig_char1_xrd_single_pattern"
        save_publication_figure(fig, name, output_dirs=[FIGURES_XRD_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [xrd_pristine], [FIGURES_XRD_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "XRD", {"sample": "Pt/C Pristine", "radiation": "Cu_Ka"})
        plt.close(fig)
        generated.append(name)

    # 2. XRD Stacked Offset Comparison
    xrd_treated = char_base / "xrd" / "catalyst_treated_xrd.txt"
    if xrd_pristine.exists() and xrd_treated.exists():
        d_p = parse_xrd_file(xrd_pristine)
        d_t = parse_xrd_file(xrd_treated)
        fig, ax = plot_xrd_stacked([d_p, d_t], labels=["Pristine Catalyst", "Post-Thermal Treatment"], vertical_offset_fraction=0.35, title="Catalyst Phase Transformation: As-Synthesized vs Thermally Activated")
        name = "fig_char2_xrd_stacked_comparison"
        save_publication_figure(fig, name, output_dirs=[FIGURES_XRD_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [xrd_pristine, xrd_treated], [FIGURES_XRD_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "XRD", {"comparison": "thermal_treatment"})
        plt.close(fig)
        generated.append(name)

    # 3. XRD Peak Deconvolution & Scherrer Domain Size
    if xrd_pristine.exists():
        d_p = parse_xrd_file(xrd_pristine)
        mask = (d_p.x >= 36.0) & (d_p.x <= 44.0)
        sub_x, sub_y = d_p.x[mask], d_p.y[mask]
        fit_res = fit_peaks(sub_x, sub_y, peak_centers=[39.76], model_type="pseudo_voigt", background="linear")
        center_2th = fit_res["peaks"][0]["center"]
        fwhm = fit_res["peaks"][0]["fwhm"]
        size_nm = calculate_crystallite_size_scherrer(center_2th, fwhm, instrumental_fwhm_deg=0.08)
        d_sub = d_p.copy()
        d_sub.x = sub_x
        d_sub.y = sub_y
        size_val = size_nm["crystallite_size_nm"] if size_nm["crystallite_size_nm"] is not None else 0.0
        fig, (ax_main, ax_res) = plot_xrd_fit(d_sub, fit_res, title=f"Pt (111) Bragg Reflection Deconvolution & Scherrer Size ({size_val:.1f} nm)")
        name = "fig_char3_xrd_peak_deconvolution"
        save_publication_figure(fig, name, output_dirs=[FIGURES_XRD_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [xrd_pristine], [FIGURES_XRD_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "XRD", {"reflection": "Pt(111)", "crystallite_size_nm": size_val})
        plt.close(fig)
        generated.append(name)

    # 4. FTIR Absorbance Spectrum
    ftir_func = char_base / "ftir" / "catalyst_functionalized_ftir.csv"
    if ftir_func.exists():
        d_ftir = parse_ftir_file(transmittance_unit="percent", filepath=ftir_func)
        d_abs = convert_transmittance_to_absorbance(d_ftir)
        fig, ax = plot_ftir_spectrum(d_abs, annotate_peaks=True, show_regions=True, title="Functionalized Catalyst Surface FTIR: Surface Ligand Identification")
        name = "fig_char4_ftir_absorbance_spectrum"
        save_publication_figure(fig, name, output_dirs=[FIGURES_FTIR_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [ftir_func], [FIGURES_FTIR_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "FTIR", {"signal": "absorbance"})
        plt.close(fig)
        generated.append(name)

    # 5. FTIR Stacked Offset Comparison
    ftir_spent = char_base / "ftir" / "catalyst_spent_ftir.csv"
    if ftir_func.exists() and ftir_spent.exists():
        d1_abs = convert_transmittance_to_absorbance(parse_ftir_file(transmittance_unit="percent", filepath=ftir_func))
        d2_abs = convert_transmittance_to_absorbance(parse_ftir_file(transmittance_unit="percent", filepath=ftir_spent))
        fig, ax = plot_ftir_stacked([d1_abs, d2_abs], labels=["Pristine Functionalized", "Spent Post-Electrolysis"], vertical_offset_fraction=0.4, title="FTIR Surface Degradation: Pristine vs Spent Post-OER")
        name = "fig_char5_ftir_stacked_comparison"
        save_publication_figure(fig, name, output_dirs=[FIGURES_FTIR_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [ftir_func, ftir_spent], [FIGURES_FTIR_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "FTIR", {"comparison": "pristine_vs_spent"})
        plt.close(fig)
        generated.append(name)

    # 6. Raman Deconvolution & ID/IG Ratio
    raman_c = char_base / "raman" / "carbon_support_raman.csv"
    if raman_c.exists():
        d_rm = parse_raman_file(raman_c)
        proc_rm = process_raman_spectrum(d_rm, despike=True, subtract_bg=True, bg_method="asls")
        fit_res_rm = fit_peaks(proc_rm.x, proc_rm.y, peak_centers=[1345.0, 1585.0], model_type="lorentzian", background="linear")
        ratios = calculate_peak_ratios(fit_res_rm, peak_1_index=1, peak_2_index=2, label_1="D", label_2="G")
        fig, (ax_main, ax_res) = plot_raman_deconvolution(proc_rm, fit_res_rm, title=f"Carbon Support Raman: D-Band & G-Band Deconvolution ($I_D/I_G={ratios['height_ratio']:.2f}$, $A_D/A_G={ratios['area_ratio']:.2f}$)")
        name = "fig_char6_raman_deconvolution_id_ig"
        save_publication_figure(fig, name, output_dirs=[FIGURES_RAMAN_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [raman_c], [FIGURES_RAMAN_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "Raman", {"id_ig_height": ratios["height_ratio"], "ad_ag_area": ratios["area_ratio"]})
        plt.close(fig)
        generated.append(name)

    # 7. Raman Stacked Comparison
    raman_tmd = char_base / "raman" / "tmd_catalyst_raman.csv"
    if raman_c.exists() and raman_tmd.exists():
        d1 = process_raman_spectrum(parse_raman_file(raman_c), despike=True, subtract_bg=True)
        d2 = process_raman_spectrum(parse_raman_file(raman_tmd), despike=True, subtract_bg=True)
        fig, ax = plot_raman_stacked([d1, d2], labels=["Carbon Support", "MoS2 Electrocatalyst"], vertical_offset_fraction=0.35, title="Raman Fingerprint Comparison: Carbon Support vs MoS2 Electrocatalyst")
        name = "fig_char7_raman_stacked_comparison"
        save_publication_figure(fig, name, output_dirs=[FIGURES_RAMAN_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [raman_c, raman_tmd], [FIGURES_RAMAN_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "Raman", {"samples": ["Carbon", "MoS2"]})
        plt.close(fig)
        generated.append(name)

    # 8. SEM Particle Segmentation
    sem_img = char_base / "microscopy" / "sem_nanoparticles.png"
    seg_res = None
    if sem_img.exists():
        d_sem = parse_microscopy_image(sem_img, pixel_size_nm=0.50)
        seg_res = segment_particles(d_sem, min_area_px=15, clear_borders=True)
        fig, ax = plot_segmentation_overlay(d_sem, seg_res, title="FESEM Catalyst Nanoparticles: Otsu Segmentation & Calibrated Scale Bar")
        name = "fig_char8_sem_particle_segmentation"
        save_publication_figure(fig, name, output_dirs=[FIGURES_MICROSCOPY_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [sem_img], [FIGURES_MICROSCOPY_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "SEM", {"particle_count": seg_res["statistics"]["particle_count"], "pixel_size_nm": 0.50})
        plt.close(fig)
        generated.append(name)

    # 9. Microscopy Particle Size Distribution
    if seg_res is not None:
        fig, ax = plot_particle_size_distribution(seg_res, bins=15, title=f"Nanoparticle Size Distribution: Mean = {seg_res['statistics']['mean_diameter_nm']:.1f} nm, Median = {seg_res['statistics']['median_diameter_nm']:.1f} nm")
        name = "fig_char9_microscopy_particle_size_distribution"
        save_publication_figure(fig, name, output_dirs=[FIGURES_MICROSCOPY_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [sem_img], [FIGURES_MICROSCOPY_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "Microscopy", {"d50_nm": seg_res["statistics"]["median_diameter_nm"]})
        plt.close(fig)
        generated.append(name)

    # 10. HRTEM Lattice Fringe 2D FFT Analysis
    hrtem_img = char_base / "microscopy" / "hrtem_lattice_crop.png"
    if hrtem_img.exists():
        d_hrtem = parse_microscopy_image(hrtem_img, pixel_size_nm=0.018)
        fringe_res = analyze_hrtem_lattice_fringes(d_hrtem, fringe_direction="horizontal")
        d_val = fringe_res["d_spacing_nm"] * 10.0
        fig, axes = plot_hrtem_fft_panel(d_hrtem, fringe_res, title=f"HRTEM Atomic Lattice Resolution: 2D FFT & Fringe Spacing ($d={d_val:.2f}$ Å)")
        name = "fig_char10_hrtem_lattice_fringe_fft"
        save_publication_figure(fig, name, output_dirs=[FIGURES_MICROSCOPY_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [hrtem_img], [FIGURES_MICROSCOPY_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "HRTEM", {"d_spacing_angstrom": fringe_res["d_spacing_nm"] * 10.0})
        plt.close(fig)
        generated.append(name)

    # 11. ICP Multi-Element Calibration & Leaching Profile
    icp_cal_file = char_base / "icp" / "icp_calibration_standards.csv"
    icp_leach_file = char_base / "icp" / "catalyst_leaching_study.csv"
    if icp_cal_file.exists() and icp_leach_file.exists():
        d_cal = parse_icp_file(icp_cal_file)
        df_cal = d_cal.dataframe
        df_ni = df_cal[df_cal["element"] == "Ni"]
        stds_ni = df_ni[df_ni["sample_id"].str.startswith("STD")]
        blanks_ni = df_ni[df_ni["sample_id"].str.startswith("BLANK")]
        cal_ni = fit_icp_calibration(
            standards_concentration=stds_ni["concentration"].values,
            standards_intensity=stds_ni["intensity"].values,
            blank_intensities=blanks_ni["intensity"].values,
        )

        d_leach = parse_icp_file(icp_leach_file)
        df_l = d_leach.dataframe
        df_l_ni = df_l[df_l["element"] == "Ni"]
        pct_leached = []
        sample_ids = df_l_ni["sample_id"].tolist()
        for _, row in df_l_ni.iterrows():
            lr = calculate_catalyst_leaching(
                measured_concentration_mg_L=float(row["concentration"]),
                electrolyte_volume_mL=50.0,
                catalyst_mass_mg=2.0,
                initial_metal_wt_fraction=0.20,
                dilution_factor=float(row.get("dilution_factor", 1.0)),
            )
            pct_leached.append(lr["percentage_metal_leached"])

        fig, (ax_l, ax_r) = create_publication_figure(columns="double", aspect_ratio="2.2", ncols=2)
        # Left: Ni calibration line
        ax_l.scatter(cal_ni["standards_c"], cal_ni["standards_i"], color="#0072B2", s=25, label="Ni Standards", zorder=3)
        ax_l.plot(cal_ni["standards_c"], cal_ni["fit_i"], color="#D55E00", linewidth=1.4, label=f"Fit ($R^2={cal_ni['r_squared']:.4f}$)", zorder=2)
        ax_l.set_xlabel(r"Concentration $[\mathrm{mg\cdot L^{-1}}]$")
        ax_l.set_ylabel(r"Emission Intensity $[\mathrm{counts}]$")
        lod_val = cal_ni.get("lod_conc", 0.0)
        ax_l.set_title(f"(a) ICP-OES Ni Calibration (LOD = {lod_val:.3f} mg/L)", fontsize=9, pad=6)
        ax_l.legend(loc="lower right", frameon=False, fontsize=8)

        # Right: Ni leaching over time
        x_idx = np.arange(len(pct_leached))
        bars = ax_r.bar(x_idx, pct_leached, color="#D55E00", edgecolor="black", linewidth=0.8, alpha=0.85, zorder=3)
        for rect in bars:
            h = rect.get_height()
            ax_r.annotate(f"{h:.2f}%", xy=(rect.get_x() + rect.get_width() / 2.0, h), xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=7)
        ax_r.set_xticks(x_idx)
        ax_r.set_xticklabels(sample_ids, rotation=15, ha="right", fontsize=8)
        ax_r.set_ylabel(r"Catalytic Ni Leached $[\%]$")
        ax_r.set_title("(b) Catalyst Leaching Stability (1.0 M KOH, OER)", fontsize=9, pad=6)
        ax_r.set_ylim(0, max(pct_leached) * 1.35 if len(pct_leached) > 0 and max(pct_leached) > 0 else 10.0)

        name = "fig_char11_icp_calibration_and_leaching"
        save_publication_figure(fig, name, output_dirs=[FIGURES_ICP_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [icp_cal_file, icp_leach_file], [FIGURES_ICP_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "ICP", {"lod_mg_l": cal_ni["lod_conc"]})
        plt.close(fig)
        generated.append(name)

    # 12. BET Isotherm & BJH Mesopore Size Distribution
    bet_iso_file = char_base / "bet" / "mesoporous_catalyst_n2_isotherm.csv"
    if bet_iso_file.exists():
        d_iso = parse_isotherm_file(bet_iso_file)
        df_ads = d_iso.dataframe[d_iso.dataframe["branch"] == "adsorption"]
        df_des = d_iso.dataframe[d_iso.dataframe["branch"] == "desorption"]
        bet_res = calculate_bet_surface_area(df_ads["p_p0"].values, df_ads["quantity_adsorbed"].values, p_range=(0.05, 0.30))
        bjh_res = calculate_bjh_pore_distribution(df_des["p_p0"].values, df_des["quantity_adsorbed"].values, branch="desorption")

        fig, (ax_l, ax_r) = create_publication_figure(columns="double", aspect_ratio="2.2", ncols=2)
        # Left: Isotherm
        ax_l.plot(df_ads["p_p0"], df_ads["quantity_adsorbed"], marker="o", markersize=3.5, linestyle="-", color="#0072B2", label="Adsorption", zorder=3)
        ax_l.plot(df_des["p_p0"], df_des["quantity_adsorbed"], marker="s", markersize=3.5, markerfacecolor="none", linestyle="--", color="#D55E00", label="Desorption", zorder=3)
        ax_l.set_xlabel(r"Relative Pressure $p/p_0$")
        ax_l.set_ylabel(r"Quantity Adsorbed $[\mathrm{cm^3\cdot g^{-1}~STP}]$")
        ax_l.set_title(f"(a) $N_2$ Isotherm ($S_{{\\mathrm{{BET}}}} = {bet_res['s_bet_m2_g']:.1f}~\\mathrm{{m^2/g}}$)", fontsize=9, pad=6)
        ax_l.legend(loc="upper left", frameon=False, fontsize=8)
        ax_l.set_xlim(0, 1.0)
        ax_l.set_ylim(bottom=0)

        # Right: BJH pore distribution
        ax_r.plot(bjh_res["pore_diameter_nm"], bjh_res["dv_dd_cm3_g_nm"], marker="o", markersize=3.5, color="#009E73", linewidth=1.3, zorder=3)
        ax_r.set_xscale("log")
        ax_r.set_xlabel(r"Pore Diameter $D_p~[\mathrm{nm}]$")
        ax_r.set_ylabel(r"$\mathrm{d}V/\mathrm{d}D_p~[\mathrm{cm^3\cdot g^{-1}\cdot nm^{-1}}]$")
        ax_r.set_title(f"(b) BJH Mesopore Distribution (Mean $D_p = {bjh_res['mean_pore_diameter_nm']:.1f}$ nm)", fontsize=9, pad=6)
        ax_r.set_ylim(bottom=0)

        name = "fig_char12_bet_isotherm_bjh_pore_distribution"
        save_publication_figure(fig, name, output_dirs=[FIGURES_BET_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [bet_iso_file], [FIGURES_BET_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "BET", {"s_bet_m2_g": bet_res["s_bet_m2_g"], "mean_pore_d_nm": bjh_res["mean_pore_diameter_nm"]})
        plt.close(fig)
        generated.append(name)

    # 13. Cross-Technique Correlation: BET Surface Area vs Double-Layer Capacitance
    bet_areas = [45.0, 92.0, 148.0, 210.0, 280.0]
    c_dl_vals = [12.5, 23.8, 38.2, 53.0, 71.4]
    sample_tags = ["Cat-A", "Cat-B", "Cat-C", "Cat-D", "Cat-E"]
    corr_res = correlate_cross_technique_metrics(
        x_values=bet_areas,
        y_values=c_dl_vals,
        x_label=r"BET Specific Surface Area $S_{\mathrm{BET}}~[\mathrm{m^2\cdot g^{-1}}]$",
        y_label=r"Double-Layer Capacitance $C_{\mathrm{dl}}~[\mathrm{mF\cdot cm^{-2}}]$",
        sample_labels=sample_tags,
    )
    fig, ax = plot_cross_technique_correlation(corr_res, title="Cross-Technique Correlation: BET Surface Area vs Double-Layer Capacitance")
    name = "fig_char13_cross_technique_correlation"
    save_publication_figure(fig, name, output_dirs=[FIGURES_COMPARATIVE_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
    record_characterization_provenance(name, [char_base / "bet" / "mesoporous_catalyst_n2_isotherm.csv"], [FIGURES_COMPARATIVE_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "Comparative", {"r_squared": corr_res["r_squared"]})
    plt.close(fig)
    generated.append(name)

    # =========================================================================
    # Atomic Absorption Spectroscopy (AAS) Publication Figures
    # =========================================================================
    aas_base = char_base / "aas"

    # 14. FAAS Iron Calibration Curve with Residuals
    faas_file = aas_base / "faas_calibration_fe.csv"
    if faas_file.exists():
        d_faas = parse_aas_file(faas_file, technique="FAAS")
        df_f = d_faas.dataframe
        df_stds = df_f[df_f["sample_type"] == "Standard"]
        df_blanks = df_f[df_f["sample_type"] == "Blank"]
        cal_fe = fit_aas_calibration(
            standards_concentration=df_stds["concentration_reported"].values,
            standards_signal=df_stds["absorbance"].values,
            signal_type="absorbance",
            blank_signals=df_blanks["absorbance"].values,
            lod_method="blank_stdev",
        )
        fig, (ax_t, ax_b) = plot_aas_calibration(cal_fe, element="Iron (Fe 248.3 nm)", unit="mg/L")
        name = "fig_char14_aas_calibration_fe"
        save_publication_figure(fig, name, output_dirs=[FIGURES_AAS_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [faas_file], [FIGURES_AAS_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "FAAS", {"slope": cal_fe["slope"], "r_squared": cal_fe["r_squared"], "lod": cal_fe["lod_conc"]})
        plt.close(fig)
        generated.append(name)

    # 15. GFAAS Platinum Integrated Absorbance Calibration Curve with Residuals
    gfaas_file = aas_base / "gfaas_calibration_pt.csv"
    if gfaas_file.exists():
        d_gfaas = parse_aas_file(gfaas_file, technique="GFAAS")
        df_g = d_gfaas.dataframe
        df_g_stds = df_g[df_g["sample_type"] == "Standard"]
        df_g_blanks = df_g[df_g["sample_type"] == "Blank"]
        cal_pt = fit_aas_calibration(
            standards_concentration=df_g_stds["concentration_reported"].values,
            standards_signal=df_g_stds["integrated_absorbance"].values,
            signal_type="integrated_absorbance",
            blank_signals=df_g_blanks["integrated_absorbance"].values,
            lod_method="blank_stdev",
        )
        fig, (ax_t, ax_b) = plot_aas_calibration(cal_pt, element="Platinum (Pt 265.9 nm, Zeeman)", unit=r"\mu g/L")
        name = "fig_char15_aas_gfaas_calibration_pt"
        save_publication_figure(fig, name, output_dirs=[FIGURES_AAS_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [gfaas_file], [FIGURES_AAS_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "GFAAS", {"slope": cal_pt["slope"], "r_squared": cal_pt["r_squared"], "lod": cal_pt["lod_conc"]})
        plt.close(fig)
        generated.append(name)

    # 16. Method of Standard Additions (MOSA) for Matrix Suppression
    mosa_file = aas_base / "aas_standard_addition_catalyst.csv"
    if mosa_file.exists():
        d_mosa = parse_aas_file(mosa_file)
        df_m = d_mosa.dataframe
        spikes = [0.0, 0.50, 1.00, 1.50, 2.00]
        mean_sigs = [df_m[df_m["sample_id"].str.contains(f"\\+{sp:.2f}") if sp > 0 else df_m["sample_id"].str.contains("Unspiked")]["absorbance"].mean() for sp in spikes]
        sa_res = fit_standard_addition(
            added_concentrations=spikes,
            signals=mean_sigs,
            sample_volume_mL=5.0,
            final_volume_mL=25.0,
        )
        fig, ax = plot_standard_addition(sa_res, element="Cobalt (Co 240.7 nm in KOH)", unit="mg/L")
        name = "fig_char16_aas_standard_addition_matrix"
        save_publication_figure(fig, name, output_dirs=[FIGURES_AAS_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [mosa_file], [FIGURES_AAS_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "AAS-MOSA", {"apparent_conc": sa_res["apparent_concentration"], "original_conc": sa_res["original_concentration"]})
        plt.close(fig)
        generated.append(name)

    # 17. Catalyst Metal Loading: Actual vs Nominal Targets
    dig_file = aas_base / "catalyst_digestion_loading.csv"
    if dig_file.exists():
        df_dig = pd.read_csv(dig_file)
        df_samples = df_dig[df_dig["Sample_Type"] == "Sample"]
        cat_labels = ["Pt/C (Pt)", "PtCo/C (Pt)", "PtCo/C (Co)", "NiFe-LDH (Ni)", "NiFe-LDH (Fe)"]
        actual_wt = []
        nominal_wt = [20.0, 15.0, 5.0, 25.0, 12.0]
        err_wt = []

        # Pt/C
        pt_c_rows = df_samples[df_samples["Sample_ID"].str.startswith("Pt_C")]
        load_pt_c = [calculate_aas_catalyst_loading(r["Measured_Conc_mg_L"], r["Digest_Vol_mL"], r["Mass_mg"], r["Dilution_Factor"])["metal_loading_wt_pct"] for _, r in pt_c_rows.iterrows()]
        actual_wt.append(np.mean(load_pt_c))
        err_wt.append(np.std(load_pt_c, ddof=1) if len(load_pt_c) > 1 else 0.2)

        # PtCo/C Pt
        ptco_pt = df_samples[(df_samples["Sample_ID"].str.startswith("PtCo")) & (df_samples["Element"] == "Pt")]
        load_ptco_pt = [calculate_aas_catalyst_loading(r["Measured_Conc_mg_L"], r["Digest_Vol_mL"], r["Mass_mg"], r["Dilution_Factor"])["metal_loading_wt_pct"] for _, r in ptco_pt.iterrows()]
        actual_wt.append(np.mean(load_ptco_pt))
        err_wt.append(np.std(load_ptco_pt, ddof=1) if len(load_ptco_pt) > 1 else 0.2)

        # PtCo/C Co
        ptco_co = df_samples[(df_samples["Sample_ID"].str.startswith("PtCo")) & (df_samples["Element"] == "Co")]
        load_ptco_co = [calculate_aas_catalyst_loading(r["Measured_Conc_mg_L"], r["Digest_Vol_mL"], r["Mass_mg"], r["Dilution_Factor"])["metal_loading_wt_pct"] for _, r in ptco_co.iterrows()]
        actual_wt.append(np.mean(load_ptco_co))
        err_wt.append(np.std(load_ptco_co, ddof=1) if len(load_ptco_co) > 1 else 0.1)

        # NiFe Ni
        nife_ni = df_samples[(df_samples["Sample_ID"].str.startswith("NiFe")) & (df_samples["Element"] == "Ni")]
        load_nife_ni = [calculate_aas_catalyst_loading(r["Measured_Conc_mg_L"], r["Digest_Vol_mL"], r["Mass_mg"], r["Dilution_Factor"])["metal_loading_wt_pct"] for _, r in nife_ni.iterrows()]
        actual_wt.append(np.mean(load_nife_ni))
        err_wt.append(np.std(load_nife_ni, ddof=1) if len(load_nife_ni) > 1 else 0.3)

        # NiFe Fe
        nife_fe = df_samples[(df_samples["Sample_ID"].str.startswith("NiFe")) & (df_samples["Element"] == "Fe")]
        load_nife_fe = [calculate_aas_catalyst_loading(r["Measured_Conc_mg_L"], r["Digest_Vol_mL"], r["Mass_mg"], r["Dilution_Factor"])["metal_loading_wt_pct"] for _, r in nife_fe.iterrows()]
        actual_wt.append(np.mean(load_nife_fe))
        err_wt.append(np.std(load_nife_fe, ddof=1) if len(load_nife_fe) > 1 else 0.2)

        fig, ax = plot_aas_catalyst_loading(
            sample_names=cat_labels,
            actual_loading_wt_pct=actual_wt,
            nominal_loading_wt_pct=nominal_wt,
            error_bars=err_wt,
            element="Metal",
            title="Catalyst Metal Loading: AAS Digest vs Nominal Target",
        )
        name = "fig_char17_aas_catalyst_loading"
        save_publication_figure(fig, name, output_dirs=[FIGURES_AAS_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [dig_file], [FIGURES_AAS_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "AAS-Digestion", {"catalysts": cat_labels})
        plt.close(fig)
        generated.append(name)

    # 18. Catalyst Leaching Kinetics during OER Stability Testing
    leach_file = aas_base / "catalyst_stability_leaching.csv"
    if leach_file.exists():
        df_l = pd.read_csv(leach_file)
        df_co = df_l[df_l["Element"] == "Co"]
        time_pts = df_co["Time_hours"].tolist()
        leach_co_pct = []
        mass_loss_ug_mg = []
        aliquots = []
        for _, r in df_co.iterrows():
            lr = calculate_aas_leaching(
                measured_concentration_mg_L=float(r["Measured_Conc_mg_L"]),
                electrolyte_volume_mL=float(r["Electrolyte_Volume_mL"]),
                catalyst_mass_mg=float(r["Initial_Catalyst_Mass_mg"]),
                initial_metal_wt_fraction=float(r["Initial_Metal_wt_pct"]) / 100.0,
                aliquot_history=aliquots if len(aliquots) > 0 else None,
            )
            leach_co_pct.append(lr["percentage_metal_leached"])
            mass_loss_ug_mg.append(lr["mass_normalized_leaching_ug_per_mg"])
            if float(r["Aliquot_Withdrawn_mL"]) > 0:
                aliquots.append((float(r["Aliquot_Withdrawn_mL"]), float(r["Measured_Conc_mg_L"])))

        fig, ax = plot_aas_leaching(
            time_or_cycles=[f"{t}h" for t in time_pts],
            leached_percentages=leach_co_pct,
            mass_normalized_loss_ug_per_mg=mass_loss_ug_mg,
            element="Cobalt",
            title="PtCo/C Stability: Cobalt Dissolution Kinetics (0.1 M KOH, 100h)",
        )
        name = "fig_char18_aas_stability_leaching_kinetics"
        save_publication_figure(fig, name, output_dirs=[FIGURES_AAS_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [leach_file], [FIGURES_AAS_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "AAS-Leaching", {"total_pct_leached": leach_co_pct[-1]})
        plt.close(fig)
        generated.append(name)

    # 19. Catalyst Metal Retention Before vs After Stability Stress
    if leach_file.exists():
        cat_names = ["Pt in PtCo/C", "Co in PtCo/C"]
        before_wt = [14.8, 5.1]
        after_wt = [14.8 * (1.0 - 0.044), 5.1 * (1.0 - 0.206)]
        fig, ax = plot_aas_stability_before_after(
            sample_names=cat_names,
            before_loading_wt_pct=before_wt,
            after_loading_wt_pct=after_wt,
            element="Metal",
            title="Electrocatalyst Durability: Metal Retention After 100h OER",
        )
        name = "fig_char19_aas_stability_before_after"
        save_publication_figure(fig, name, output_dirs=[FIGURES_AAS_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [leach_file], [FIGURES_AAS_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "AAS-Retention", {"before": before_wt, "after": after_wt})
        plt.close(fig)
        generated.append(name)

    # 20. Replicate Precision & Instrumental Dispersion
    if faas_file.exists():
        df_rep = pd.read_csv(faas_file)
        df_rep_stds = df_rep[df_rep["Sample_Type"] == "Standard"]
        s_names = ["Std 0.5 mg/L", "Std 1.0 mg/L", "Std 2.0 mg/L", "Std 5.0 mg/L"]
        rep_data = [
            df_rep_stds[df_rep_stds["Sample_ID"] == "STD_0.50"]["Absorbance"].tolist(),
            df_rep_stds[df_rep_stds["Sample_ID"] == "STD_1.00"]["Absorbance"].tolist(),
            df_rep_stds[df_rep_stds["Sample_ID"] == "STD_2.00"]["Absorbance"].tolist(),
            df_rep_stds[df_rep_stds["Sample_ID"] == "STD_5.00"]["Absorbance"].tolist(),
        ]
        fig, ax = plot_aas_replicate_distribution(
            sample_names=s_names,
            replicate_measurements=rep_data,
            unit="AU",
            title="Flame AAS Instrumental Precision Across Replicate Standards",
        )
        name = "fig_char20_aas_replicate_distribution"
        save_publication_figure(fig, name, output_dirs=[FIGURES_AAS_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [faas_file], [FIGURES_AAS_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "FAAS-Precision", {"levels": s_names})
        plt.close(fig)
        generated.append(name)

    # 21. Analytical QC Run Sequence Control Chart (EPA 7000B)
    qc_file = aas_base / "aas_qc_run_sequence.csv"
    if qc_file.exists():
        df_qc = pd.read_csv(qc_file)
        fig, ax = plot_aas_qc_chart(
            ccv_labels=df_qc["QC_Check"].tolist(),
            ccv_recoveries_pct=df_qc["Recovery_pct"].tolist(),
            sop="EPA_7000B",
            title="Analytical Quality Control: CCV Recoveries (US EPA 7000B)",
        )
        name = "fig_char21_aas_qc_control_chart"
        save_publication_figure(fig, name, output_dirs=[FIGURES_AAS_DIR, FIGURES_CHAR_DIR], formats=formats, dpi=dpi)
        record_characterization_provenance(name, [qc_file], [FIGURES_AAS_DIR, FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], "AAS-QC", {"sop": "EPA_7000B", "mean_recovery": float(df_qc["Recovery_pct"].mean())})
        plt.close(fig)
        generated.append(name)

    return generated


def handle_electrochemistry_cli(args: argparse.Namespace) -> None:
    """Handle direct CLI commands for electrochemistry techniques."""
    input_path = Path(args.input)
    if not input_path.is_file():
        raise FileNotFoundError(f"Data file not found: {input_path}")

    data = load_electrochemistry_data(input_path, instrument=args.instrument, reference_electrode=args.reference_electrode, ph=args.ph, temperature_k=args.temperature_k)
    out_name = args.output or f"{input_path.stem}_{args.type}_plot"
    out_dirs = [FIGURES_ELECTROCHEM_DIR]

    if args.type == "lsv":
        fig, ax = create_publication_figure(columns=args.columns or "single", aspect_ratio="golden")
        plot_lsv(data, reaction=args.reaction, ax=ax, title=f"{data.instrument.title()} {args.reaction} Polarization")
    elif args.type == "cv":
        fig, ax = create_publication_figure(columns=args.columns or "single", aspect_ratio="golden")
        plot_cv(data, ax=ax, mark_peaks=True, title=f"{data.instrument.title()} Cyclic Voltammetry")
    elif args.type == "eis":
        fit = fit_equivalent_circuit(data.frequency_Hz, data.impedance_real_ohm, data.impedance_imag_ohm)
        fig = plot_eis_comprehensive(data, fit_result=fit, title=f"{data.instrument.title()} EIS Analysis")
    elif args.type == "ca":
        fig, ax = plot_chronoamperometry(data, show_retention=True, title=f"{data.instrument.title()} Chronoamperometry")
    elif args.type == "cp":
        fig, ax = plot_chronopotentiometry(data, show_drift_rate=True, title=f"{data.instrument.title()} Chronopotentiometry")
    elif args.type == "tafel":
        fig, ax = create_publication_figure(columns=args.columns or "single", aspect_ratio="golden")
        plot_tafel(data, reaction=args.reaction, ax=ax, title=f"{data.instrument.title()} Tafel Analysis")
    elif args.type == "benchmark":
        if not args.second_input:
            raise ValueError("--second-input PATH is required for benchmark comparison.")
        data2 = load_electrochemistry_data(args.second_input, reference_electrode=args.reference_electrode, ph=args.ph, temperature_k=args.temperature_k)
        fig, axes = plot_cross_instrument_benchmark(data, data2, reaction=args.reaction)
    else:
        raise ValueError(f"Unsupported electrochemistry type: {args.type}")

    save_publication_figure(fig, out_name, output_dirs=out_dirs, formats=args.formats, dpi=args.dpi)
    record_figure_provenance(out_name, [input_path] + ([Path(args.second_input)] if args.second_input else []), [FIGURES_ELECTROCHEM_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], {"type": args.type, "instrument": data.instrument, "conditions": data.raw_metadata.get("resolved_conditions", {}), "ph": data.ph, "potential_reference_raw": data.potential_reference_raw})
    plt.close(fig)
    print(f"Successfully generated electrochemistry figure: {out_name}")


def handle_characterization_cli(args: argparse.Namespace) -> None:
    """Handle direct CLI commands for materials characterization techniques."""
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Characterization data file not found: {input_path}")

    out_name = args.output or f"{input_path.stem}_{args.type}_plot"
    out_dirs = [FIGURES_CHAR_DIR]

    if args.type == "xrd":
        if args.second_input:
            d1 = parse_xrd_file(input_path)
            d2 = parse_xrd_file(args.second_input)
            fig, ax = plot_xrd_stacked([d1, d2], labels=[input_path.stem, Path(args.second_input).stem])
        else:
            data = parse_xrd_file(input_path)
            fig, ax = plot_xrd_pattern(data)
    elif args.type == "ftir":
        if args.second_input:
            d1 = convert_transmittance_to_absorbance(parse_ftir_file(input_path, transmittance_unit=args.transmittance_unit, signal_type=args.signal_type))
            d2 = convert_transmittance_to_absorbance(parse_ftir_file(args.second_input, transmittance_unit=args.transmittance_unit, signal_type=args.signal_type))
            fig, ax = plot_ftir_stacked([d1, d2], labels=[input_path.stem, Path(args.second_input).stem])
        else:
            data = parse_ftir_file(input_path, transmittance_unit=args.transmittance_unit, signal_type=args.signal_type)
            data_abs = convert_transmittance_to_absorbance(data)
            fig, ax = plot_ftir_spectrum(data_abs, annotate_peaks=True)
    elif args.type == "raman":
        if args.second_input:
            d1 = process_raman_spectrum(parse_raman_file(input_path), despike=True, subtract_bg=True)
            d2 = process_raman_spectrum(parse_raman_file(args.second_input), despike=True, subtract_bg=True)
            fig, ax = plot_raman_stacked([d1, d2], labels=[input_path.stem, Path(args.second_input).stem])
        else:
            data = parse_raman_file(input_path)
            proc = process_raman_spectrum(data, despike=True, subtract_bg=True)
            fig, ax = plot_raman_spectrum(proc)
    elif args.type in ("sem", "tem"):
        px_nm = args.pixel_size
        data = parse_microscopy_image(input_path, pixel_size_nm=px_nm, technique=args.type.upper())
        if data.pixel_size_nm is None:
            raise ValueError("Microscopy calibration is required; supply confirmed --pixel-size in nm/px.")
        seg = segment_particles(data)
        fig, axes = plot_segmentation_overlay(data, seg)
    elif args.type == "hrtem":
        px_nm = args.pixel_size
        data = parse_microscopy_image(input_path, pixel_size_nm=px_nm, technique="HRTEM")
        fringes = analyze_hrtem_lattice_fringes(data)
        fig, axes = plot_hrtem_fft_panel(data, fringes)
    elif args.type == "icp":
        data = parse_icp_file(input_path)
        if "intensity" in (data.dataframe.columns if data.dataframe is not None else []):
            elem = args.elements[0] if args.elements else "Pt"
            df_elem = data.dataframe[data.dataframe["element"] == elem] if "element" in data.dataframe.columns else data.dataframe
            stds = df_elem[df_elem["sample_id"].str.startswith("STD")] if not df_elem[df_elem["sample_id"].str.startswith("STD")].empty else df_elem
            cal = fit_icp_calibration(stds["concentration"].values, stds["intensity"].values)
            fig, (ax_top, ax_bot) = plot_icp_calibration(cal, element=elem)
        else:
            raise ValueError("ICP leaching requires confirmed solution volume, catalyst mass, and metal fraction. Use calculate_catalyst_leaching with explicit experimental conditions.")
    elif args.type == "aas":
        data = parse_aas_file(input_path)
        df = data.dataframe
        if not any(column in df and np.isfinite(df[column]).any()
                   for column in ("absorbance", "integrated_absorbance")):
            raise ValueError("AAS calibration requires measured absorbance signals. Concentration-only loading/leaching requires confirmed sample conditions through the existing Python analysis functions.")
        out_dirs = [FIGURES_AAS_DIR, FIGURES_CHAR_DIR]
        method = getattr(args, "calibration_method", "external") or "external"

        has_added = any("added" in c.lower() for c in df.columns)
        if method == "standard_addition" or has_added:
            add_col = next(c for c in df.columns if "added" in c.lower())
            sig_col = "integrated_absorbance" if "integrated_absorbance" in df.columns and not df["integrated_absorbance"].isna().all() else "absorbance"
            grouped = df.groupby(add_col)[sig_col].mean().reset_index()
            sa_res = fit_standard_addition(grouped[add_col].values, grouped[sig_col].values, signal_type=sig_col)
            elem = args.elements[0] if args.elements else df["element"].iloc[0]
            unit = df["unit"].iloc[0] if "unit" in df.columns else "mg/L"
            fig, ax = plot_standard_addition(sa_res, element=elem, unit=unit)
        elif any(c.lower() in ("mass_mg", "digest_vol_ml") for c in df.columns):
            raise ValueError("AAS loading requires confirmed digestion volume, catalyst mass, and dilution. Use calculate_aas_catalyst_loading with explicit experimental conditions.")
        elif "time_hours" in [c.lower() for c in df.columns]:
            raise ValueError("AAS leaching requires confirmed solution volume, catalyst mass, and metal fraction. Use calculate_aas_leaching with explicit experimental conditions.")
        else:
            sig_col = "integrated_absorbance" if "integrated_absorbance" in df.columns and not df["integrated_absorbance"].isna().all() else "absorbance"
            df_stds = df[df["sample_type"].str.lower().str.contains("standard|cal")]
            if df_stds.empty:
                df_stds = df
            c_col = "concentration_standard" if "concentration_standard" in df_stds.columns else ("concentration_reported" if "concentration_reported" in df_stds.columns else "concentration")
            df_blanks = df[df["sample_type"].str.lower().str.contains("blank")]
            b_sigs = df_blanks[sig_col].dropna().values if not df_blanks.empty else None
            cal_res = fit_aas_calibration(df_stds[c_col].values, df_stds[sig_col].values, signal_type=sig_col, blank_signals=b_sigs)
            elem = args.elements[0] if args.elements else df["element"].iloc[0]
            unit = df["unit"].iloc[0] if "unit" in df.columns else "mg/L"
            fig, (ax_top, ax_bot) = plot_aas_calibration(cal_res, element=elem, unit=unit)
    elif args.type == "bet":
        data = parse_isotherm_file(input_path)
        fig, ax = plot_isotherm(data)
    elif args.type == "comparative":
        if not args.second_input:
            raise ValueError("--second-input PATH is required for cross-technique correlation comparison.")
        df1 = load_csv_data(input_path)
        df2 = load_csv_data(args.second_input)
        col1 = df1.columns[1] if len(df1.columns) > 1 else df1.columns[0]
        col2 = df2.columns[1] if len(df2.columns) > 1 else df2.columns[0]
        corr = correlate_cross_technique_metrics(df1[col1].values, df2[col2].values, str(col1), str(col2))
        fig, ax = plot_cross_technique_correlation(corr)
    else:
        raise ValueError(f"Unsupported characterization type: {args.type}")

    save_publication_figure(fig, out_name, output_dirs=out_dirs, formats=args.formats, dpi=args.dpi)
    sources = [input_path] + ([Path(args.second_input)] if args.second_input else [])
    metadata = {"transmittance_unit": args.transmittance_unit, "signal_type": args.signal_type,
                "pixel_size_nm": args.pixel_size, "calibration_method": args.calibration_method}
    if args.type == "ftir":
        metadata["processed_datasets"] = [d.to_provenance_dict() for d in ([d1, d2] if args.second_input else [data_abs])]
    record_characterization_provenance(out_name, sources, [FIGURES_CHAR_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR], args.type.upper(), metadata)
    plt.close(fig)
    print(f"Successfully generated characterization figure: {out_name}")


# =========================================================================
# Quantum ESPRESSO Direct Plotting Command Handlers
# =========================================================================

def handle_qe_dos_cli(args: argparse.Namespace) -> None:
    """Handle CLI command: python -m scripts.main_plots --type dos --input PATH."""
    input_path = Path(args.input)
    if not input_path.is_file():
        raise FileNotFoundError(f"DOS file not found: {input_path}")

    dos_data = parse_dos_file(input_path, e_fermi=args.fermi)

    fig, ax = create_publication_figure(columns=args.columns or "single", aspect_ratio="golden")
    plot_qe_dos(
        dos_data=dos_data,
        ax=ax,
        e_fermi=args.fermi,
        energy_range=tuple(args.energy_range) if args.energy_range else (-6.0, 6.0),
    )

    out_name = args.output or f"{input_path.stem}_dos_plot"
    out_dirs = [FIGURES_DFT_DIR]
    save_publication_figure(fig, out_name, output_dirs=out_dirs, formats=args.formats, dpi=args.dpi)
    record_figure_provenance(
        out_name,
        [input_path],
        [FIGURES_DFT_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR],
        {"type": "dos.x", "fermi_energy": args.fermi, "energy_range": args.energy_range},
    )
    plt.close(fig)
    print(f"Successfully generated Quantum ESPRESSO DOS figure: {out_name}")


def handle_qe_pdos_cli(args: argparse.Namespace) -> None:
    """Handle CLI command: python -m scripts.main_plots --type pdos --input DIRECTORY."""
    input_dir = Path(args.input)
    if not input_dir.is_dir():
        raise NotADirectoryError(f"PDOS directory not found: {input_dir}")

    coll = QEPdosCollection(input_dir, prefix=args.prefix, e_fermi=args.fermi)

    # Build projection specifications from CLI arguments
    specs: Dict[str, Dict[str, Any]] = {}
    if args.elements:
        for elem in args.elements:
            if args.orbitals:
                for orb in args.orbitals:
                    # Check if (elem, orb) pair exists in collection
                    has_match = any(
                        wfc.element.lower() == elem.lower() and wfc.orbital.lower() == orb.lower()
                        for wfc in coll.wfc_list
                    )
                    if not has_match:
                        continue
                    if args.components:
                        for comp in args.components:
                            specs[f"{elem}-{orb}({comp})"] = {
                                "element": elem,
                                "orbital": orb,
                                "component": comp,
                                "atom_index": args.atoms,
                            }
                    else:
                        specs[f"{elem}-{orb}"] = {
                            "element": elem,
                            "orbital": orb,
                            "atom_index": args.atoms,
                        }
            else:
                specs[elem] = {"element": elem, "atom_index": args.atoms}

        if not specs:
            raise ValueError(
                f"No wavefunctions matched specified elements ({args.elements}) "
                f"and orbitals ({args.orbitals}) in collection."
            )
    else:
        # Auto-discover projections
        seen = set()
        for wfc in coll.wfc_list:
            pair = (wfc.element, wfc.orbital)
            if pair not in seen:
                seen.add(pair)
                specs[f"{wfc.element}-{wfc.orbital}"] = {"element": wfc.element, "orbital": wfc.orbital}

    fig, ax = create_publication_figure(columns=args.columns or "single", aspect_ratio="golden")
    plot_qe_pdos(
        coll,
        specs=specs,
        ax=ax,
        e_fermi=args.fermi,
        energy_range=tuple(args.energy_range) if args.energy_range else (-6.0, 6.0),
    )

    out_name = args.output or f"{coll.prefix or 'qe'}_pdos_plot"
    out_dirs = [FIGURES_DFT_DIR]
    save_publication_figure(fig, out_name, output_dirs=out_dirs, formats=args.formats, dpi=args.dpi)
    all_pdos_files = list(input_dir.glob(f"{coll.prefix}*.pdos*"))
    record_figure_provenance(
        out_name,
        all_pdos_files,
        [FIGURES_DFT_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR],
        {"type": "projwfc.x", "specs": specs, "fermi_energy": args.fermi},
    )
    plt.close(fig)
    print(f"Successfully generated Quantum ESPRESSO PDOS figure: {out_name}")


def handle_qe_bands_cli(args: argparse.Namespace) -> None:
    """Handle CLI command: python -m scripts.main_plots --type bands --input PATH."""
    input_path = Path(args.input)
    if not input_path.is_file():
        raise FileNotFoundError(f"Band structure file not found: {input_path}")

    bands_data = parse_bands_gnu(input_path, e_fermi=args.fermi)

    fig, ax = create_publication_figure(columns=args.columns or "single", aspect_ratio="golden")
    plot_qe_bands(
        bands_data,
        ax=ax,
        e_fermi=args.fermi,
        energy_range=tuple(args.energy_range) if args.energy_range else (-5.0, 5.0),
    )

    out_name = args.output or f"{input_path.stem}_bands_plot"
    out_dirs = [FIGURES_DFT_DIR]
    save_publication_figure(fig, out_name, output_dirs=out_dirs, formats=args.formats, dpi=args.dpi)
    record_figure_provenance(
        out_name,
        [input_path],
        [FIGURES_DFT_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR],
        {"type": "bands.x", "fermi_energy": args.fermi},
    )
    plt.close(fig)
    print(f"Successfully generated Quantum ESPRESSO Band Structure figure: {out_name}")


def handle_qe_neb_cli(args: argparse.Namespace) -> None:
    """Handle CLI command: python -m scripts.main_plots --type neb --input PATH."""
    input_path = Path(args.input)
    if not input_path.is_file():
        raise FileNotFoundError(f"NEB file not found: {input_path}")

    neb_data = parse_neb_output(input_path, log_filepath=args.neb_log)

    fig, ax = create_publication_figure(columns=args.columns or "single", aspect_ratio="golden")
    plot_qe_neb_profile(neb_data, ax=ax, title="NEB Electronic Energy Profile")

    out_name = args.output or f"{input_path.stem}_neb_plot"
    out_dirs = [FIGURES_CATALYSIS_DIR]
    save_publication_figure(fig, out_name, output_dirs=out_dirs, formats=args.formats, dpi=args.dpi)
    record_figure_provenance(
        out_name,
        neb_data.raw_metadata["source_files"],
        [FIGURES_CATALYSIS_DIR, OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR],
        {
            "type": "neb.x",
            "is_converged": neb_data.is_converged,
            "coordinate_unit": "arbitrary", "energy_unit": "eV",
            "barrier_source": neb_data.raw_metadata["barrier_source"],
            "forward_barrier_ev": neb_data.forward_barrier_ev,
            "reverse_barrier_ev": neb_data.reverse_barrier_ev,
        },
    )
    plt.close(fig)
    print(f"Successfully generated Quantum ESPRESSO NEB figure: {out_name}")


# =========================================================================
# Main Entry Point
# =========================================================================

from contextlib import contextmanager
from functools import partial


@contextmanager
def _cli_export_options(args):
    """Temporarily configure only this CLI's outputs, never its source datasets."""
    namespace = globals()
    names = [name for name in namespace if (name.startswith("FIGURES_") or name.startswith("OUTPUTS_")) and name.endswith("_DIR")]
    names += ["save_publication_figure", "record_figure_provenance", "record_characterization_provenance"]
    original = {name: namespace[name] for name in names}
    import scripts.utils as utils
    original_defaults = (utils.OUTPUTS_VECTOR_DIR, utils.OUTPUTS_RASTER_DIR)
    try:
        if args.output_root is not None:
            root = args.output_root.resolve()
            for name in names:
                if name.endswith("_DIR"):
                    namespace[name] = root / original[name].relative_to(PROJECT_ROOT)
        utils.OUTPUTS_VECTOR_DIR, utils.OUTPUTS_RASTER_DIR = OUTPUTS_VECTOR_DIR, OUTPUTS_RASTER_DIR
        namespace["save_publication_figure"] = partial(original["save_publication_figure"], overwrite=args.overwrite, close=True)
        for name in ("record_figure_provenance", "record_characterization_provenance"):
            namespace[name] = partial(original[name], overwrite=args.overwrite)
        yield
    finally:
        namespace.update(original)
        utils.OUTPUTS_VECTOR_DIR, utils.OUTPUTS_RASTER_DIR = original_defaults


def main() -> None:
    before = set(plt.get_fignums())
    try:
        _main()
    finally:
        for number in set(plt.get_fignums()) - before:
            plt.close(number)


def _main() -> None:
    parser = argparse.ArgumentParser(
        description="Scientific Publication Graphing Workspace Runner (Benchmark & Quantum ESPRESSO Mode)"
    )

    # Mode 1: Benchmark Batch Generation
    parser.add_argument(
        "--all",
        action="store_true",
        help="Generate all scientific benchmark figures across all categories.",
    )
    parser.add_argument(
        "--category",
        choices=["dft", "catalysis", "chemical_engineering", "electrochemistry", "characterization", "all"],
        default=None,
        help="Domain category for benchmark figures (dft, catalysis, chemical_engineering, electrochemistry, characterization).",
    )

    # Mode 2: Direct Plotting (Quantum ESPRESSO, Electrochemistry & Characterization)
    parser.add_argument(
        "--type",
        choices=[
            "dos", "pdos", "bands", "neb",
            "cv", "lsv", "eis", "ca", "cp", "tafel", "ecsa", "benchmark",
            "xrd", "ftir", "raman", "sem", "tem", "hrtem", "icp", "aas", "bet", "comparative",
        ],
        default=None,
        help="Calculation, electrochemical, or characterization technique type to parse and plot.",
    )
    parser.add_argument(
        "--calibration-method",
        choices=["external", "standard_addition"],
        default="external",
        help="AAS calibration method ('external' or 'standard_addition').",
    )
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to YAML configuration file for complex batch experiments.",
    )
    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help="Path to output file or directory.",
    )
    parser.add_argument(
        "--second-input",
        type=str,
        default=None,
        help="Second input path for comparative benchmarking.",
    )
    parser.add_argument(
        "--instrument",
        choices=["metrohm", "biologic"],
        default=None,
        help="Instrument type ('metrohm' or 'biologic').",
    )
    parser.add_argument(
        "--reaction",
        choices=["OER", "HER"],
        default="OER",
        help="Electrochemical reaction ('OER' or 'HER').",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default=None,
        help="Calculation prefix for PDOS discovery.",
    )
    parser.add_argument(
        "--fermi",
        "-ef",
        type=float,
        default=None,
        help="Fermi energy in eV for energy alignment (E - E_F = 0).",
    )
    parser.add_argument(
        "--atoms",
        nargs="+",
        type=int,
        default=None,
        help="Atom indices to select for PDOS summation (1-indexed).",
    )
    parser.add_argument(
        "--elements",
        nargs="+",
        type=str,
        default=None,
        help="Element symbols to select for PDOS summation (e.g. Pt O).",
    )
    parser.add_argument(
        "--orbitals",
        nargs="+",
        type=str,
        default=None,
        help="Orbitals to plot (e.g. s p d f).",
    )
    parser.add_argument(
        "--components",
        nargs="+",
        type=str,
        default=None,
        help="Specific orbital components (e.g. dz2 dx2-y2 eg t2g pz px py).",
    )
    parser.add_argument(
        "--energy-range",
        "-er",
        nargs=2,
        type=float,
        default=None,
        help="Energy range limits in eV (e.g. -6.0 6.0).",
    )
    parser.add_argument(
        "--columns",
        type=str,
        default="single",
        choices=["single", "double", "1.5"],
        help="Figure width preset: 'single' (85 mm), 'double' (175 mm), or '1.5' (115 mm).",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        help="Custom base name for exported figures.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=DEFAULT_DPI,
        help="Raster output DPI (default: 300).",
    )
    parser.add_argument(
        "--pixel-size",
        type=float,
        default=None,
        help="Microscopy pixel size calibration in nm/pixel.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["pdf", "svg", "png", "tiff"],
        help="Export formats (default: pdf svg png tiff).",
    )

    parser.add_argument("--output-root", type=Path, help="Root for generated figures/ and outputs/; input data stays in place.")
    parser.add_argument("--overwrite", action="store_true", help="Explicitly replace generated figures and sidecars.")
    parser.add_argument("--reference-electrode", help="Confirmed reference including filling solution.")
    parser.add_argument("--temperature-k", type=float, help="Confirmed or explicitly assumed cell temperature in kelvin.")
    parser.add_argument("--ph", type=float, help="Confirmed electrolyte pH.")
    parser.add_argument("--transmittance-unit", choices=["percent", "fraction"], help="Required for unlabeled FTIR transmittance.")
    parser.add_argument("--signal-type", choices=["transmittance", "absorbance", "reflectance"])
    parser.add_argument("--neb-log", type=Path, help="Matching NEB stdout .out log, never a .path restart.")
    args = parser.parse_args()

    with _cli_export_options(args):
        # Route based on flags
        if args.type is not None:
            if args.input is None:
                parser.error(f"--type {args.type} requires --input PATH/DIRECTORY.")

            if args.type in ("dos", "pdos", "bands", "neb"):
                print("=====================================================================")
                print(f" Quantum ESPRESSO Direct Plotting: Type={args.type}")
                print(f" Input: {args.input}")
                print("=====================================================================")
                if args.type == "dos":
                    handle_qe_dos_cli(args)
                elif args.type == "pdos":
                    handle_qe_pdos_cli(args)
                elif args.type == "bands":
                    handle_qe_bands_cli(args)
                elif args.type == "neb":
                    handle_qe_neb_cli(args)
            elif args.type in ("cv", "lsv", "eis", "ca", "cp", "tafel", "ecsa", "benchmark"):
                print("=====================================================================")
                print(f" Experimental Electrochemistry Plotting: Type={args.type}")
                print(f" Input: {args.input}")
                print("=====================================================================")
                handle_electrochemistry_cli(args)
            elif args.type in ("xrd", "ftir", "raman", "sem", "tem", "hrtem", "icp", "aas", "bet", "comparative"):
                print("=====================================================================")
                print(f" Materials Characterization Plotting: Type={args.type}")
                print(f" Input: {args.input}")
                print("=====================================================================")
                handle_characterization_cli(args)
            return

        # Benchmark Generation Mode
        selected_category = "all" if args.all else (args.category or "all")
        all_generated: List[str] = []

        print("=====================================================================")
        print(" Scientific Publication Graphing Workspace - Runner")
        print(f" Target Category: {selected_category} | DPI: {args.dpi}")
        print(f" Formats: {args.formats}")
        print("=====================================================================")

        if selected_category in ("dft", "all"):
            print("[1/5] Generating DFT electronic structure & adsorption figures...")
            figs = generate_dft_figures(formats=args.formats, dpi=args.dpi)
            all_generated.extend(figs)
            print(f"      Completed: {figs}")

        if selected_category in ("catalysis", "all"):
            print("[2/5] Generating catalysis reaction & volcano figures...")
            figs = generate_catalysis_figures(formats=args.formats, dpi=args.dpi)
            all_generated.extend(figs)
            print(f"      Completed: {figs}")

        if selected_category in ("chemical_engineering", "all"):
            print("[3/5] Generating chemical engineering kinetics & reactor figures...")
            figs = generate_chemeng_figures(formats=args.formats, dpi=args.dpi)
            all_generated.extend(figs)
            print(f"      Completed: {figs}")

        if selected_category in ("electrochemistry", "all"):
            print("[4/5] Generating experimental electrochemistry figures (Metrohm & BioLogic)...")
            figs = generate_electrochemistry_figures(formats=args.formats, dpi=args.dpi)
            all_generated.extend(figs)
            print(f"      Completed: {figs}")

        if selected_category in ("characterization", "all"):
            print("[5/5] Generating materials characterization figures (XRD, FTIR, Raman, SEM, HRTEM, ICP, BET)...")
            figs = generate_characterization_figures(formats=args.formats, dpi=args.dpi)
            all_generated.extend(figs)
            print(f"      Completed: {figs}")

        print("=====================================================================")
        print(f" Total generated publication figures: {len(all_generated)}")
        print(f" Output vector destination: {OUTPUTS_VECTOR_DIR}")
        print(f" Output raster destination: {OUTPUTS_RASTER_DIR}")
        print("=====================================================================")


if __name__ == "__main__":
    main()

"""
Catalytic Reaction Energy Profiles and Gibbs Free Energy Diagram Generation.

Supports multi-pathway comparisons, transition-state activation barriers,
explicit thermodynamic corrections (ZPE, TS, solvation), and the
Computational Hydrogen Electrode (CHE, Nørskov et al.) model.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.config import get_figure_dimensions
from scripts.preprocess import KB_EV_PER_K
from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


def convert_reference_potential(
    u_potential: float,
    from_ref: str,
    to_ref: str,
    ph: float = 0.0,
    temperature_k: float = 298.15,
) -> float:
    """
    Convert electrode potential between Standard Hydrogen Electrode (SHE)
    and Reversible Hydrogen Electrode (RHE).

    Nernstian relation:
    U_RHE = U_SHE + (k_B * T * ln(10) / e) * pH
    (~0.05916 V / pH unit at 298.15 K).

    Parameters
    ----------
    u_potential : float
        Electrode potential in Volts.
    from_ref : str
        Source reference scale ('SHE' or 'RHE').
    to_ref : str
        Target reference scale ('SHE' or 'RHE').
    ph : float, default 0.0
        Solution pH.
    temperature_k : float, default 298.15
        Electrochemical system temperature in Kelvin.

    Returns
    -------
    float
        Converted electrode potential in Volts.
    """
    f_ref = from_ref.strip().upper()
    t_ref = to_ref.strip().upper()

    if f_ref not in ("SHE", "RHE") or t_ref not in ("SHE", "RHE"):
        raise ValueError(
            f"Unsupported reference conversion from '{from_ref}' to '{to_ref}'. "
            "Supported scales: 'SHE', 'RHE'."
        )

    if f_ref == t_ref or ph == 0.0:
        return float(u_potential)

    nernst_shift = float(KB_EV_PER_K * temperature_k * np.log(10.0) * ph)

    if f_ref == "SHE" and t_ref == "RHE":
        return float(u_potential + nernst_shift)
    elif f_ref == "RHE" and t_ref == "SHE":
        return float(u_potential - nernst_shift)
    else:
        raise ValueError(
            f"Unsupported reference conversion from '{from_ref}' to '{to_ref}'. "
            "Supported scales: 'SHE', 'RHE'."
        )


def calculate_gibbs_free_energy(
    e_dft: float,
    zpe: float = 0.0,
    ts: float = 0.0,
    corrections: float = 0.0,
) -> float:
    """
    Calculate absolute Gibbs free energy of an intermediate state.

    Thermodynamic equation:
    G = E_DFT + ZPE - TS + corrections

    Parameters
    ----------
    e_dft : float
        DFT electronic ground-state energy in eV.
    zpe : float, default 0.0
        Zero-point vibrational energy (sum of 0.5 * h * nu) in eV.
    ts : float, default 0.0
        Temperature times entropic contribution (T * S) in eV.
    corrections : float, default 0.0
        Explicit physical corrections (e.g. implicit solvation, standard state, gas fugacity) in eV.

    Returns
    -------
    float
        Gibbs free energy in eV.
    """
    return float(e_dft + zpe - ts + corrections)


def apply_che_correction(
    delta_g0: float,
    u_potential: float = 0.0,
    ph: float = 0.0,
    n_electrons: int = 1,
    temperature_k: float = 298.15,
    step_type: str = "oxidation",
    nu_H: Optional[float] = None,
    reference_electrode: str = "RHE",
) -> float:
    """
    Apply Computational Hydrogen Electrode (CHE) shift for electrochemical proton-electron transfers.

    CHE formulation (Nørskov et al., J. Phys. Chem. B 2004, 108, 17886):
    For an elementary reaction step written as:
        Reactants + nu_H * (H+ + e-) <=> Products
    where nu_H is the signed stoichiometric coefficient of (H+ + e-) consumed (reactants):
    - nu_H > 0: Reduction step (e.g. *O + H+ + e- -> *OH has nu_H = +1)
    - nu_H < 0: Oxidation step (e.g. H2O -> *OH + H+ + e- has nu_H = -1)

    Thermodynamic potential dependence:
    - vs RHE:  mu(H+ + e-) = 0.5 * mu(H2) - e * U_RHE
      Delta G(U_RHE) = Delta G0 + nu_H * e * U_RHE
    - vs SHE:  mu(H+ + e-) = 0.5 * mu(H2) - e * U_SHE - k_B * T * ln(10) * pH
      Delta G(U_SHE, pH) = Delta G0 + nu_H * [e * U_SHE + k_B * T * ln(10) * pH]

    Parameters
    ----------
    delta_g0 : float
        Standard free energy change at U = 0 V vs RHE, pH = 0 in eV.
    u_potential : float, default 0.0
        Applied electrode potential in Volts (V).
    ph : float, default 0.0
        Solution pH (Nernstian shift: ~0.05916 V / pH unit at 298.15 K).
    n_electrons : int, default 1
        Number of transferred electrons (magnitude, ignored if nu_H is explicitly given).
    temperature_k : float, default 298.15
        Electrochemical system temperature in Kelvin.
    step_type : str, default 'oxidation'
        Reaction direction: 'oxidation' (anodic) or 'reduction' (cathodic).
        Ignored if nu_H is explicitly given.
    nu_H : float, optional
        Signed stoichiometric coefficient of (H+ + e-) consumed (reactants).
        If provided, overrides step_type and n_electrons.
    reference_electrode : str, default 'RHE'
        Potential reference scale: 'RHE' or 'SHE'.
        When 'RHE', U is on the reversible hydrogen electrode scale (pH shift is absorbed in U_RHE).
        When 'SHE', U is on the standard hydrogen electrode scale and explicit pH Nernstian shift is applied.

    Returns
    -------
    float
        Potential- and pH-dependent free energy change in eV.
    """
    if nu_H is not None:
        nu = float(nu_H)
    else:
        st = step_type.strip().lower()
        if st in ("oxidation", "anodic"):
            nu = -float(abs(n_electrons))
        elif st in ("reduction", "cathodic"):
            nu = float(abs(n_electrons))
        else:
            raise ValueError(f"Unknown step_type: '{step_type}'. Choose 'oxidation' or 'reduction'.")

    ref = reference_electrode.strip().upper()
    if ref == "RHE":
        # On RHE scale, potential is defined relative to the reversible hydrogen electrode at that pH
        return float(delta_g0 + nu * float(u_potential))
    elif ref == "SHE":
        nernst_v = KB_EV_PER_K * float(temperature_k) * np.log(10.0) * float(ph)
        return float(delta_g0 + nu * (float(u_potential) + nernst_v))
    else:
        raise ValueError(
            f"Unsupported reference_electrode: '{reference_electrode}'. Choose 'RHE' or 'SHE'."
        )


def plot_reaction_pathway(
    pathways: Dict[str, Sequence[float]],
    step_labels: Sequence[str],
    transition_states: Optional[Dict[str, Sequence[bool]]] = None,
    relative_to_initial: bool = True,
    ax: Optional[plt.Axes] = None,
    palette: str = "okabe_ito",
    plateau_width: float = 0.45,
    gap_width: float = 0.35,
    energy_label: str = r"Free Energy $\Delta G$ (eV)",
    x_label: str = "Reaction Coordinate",
    title: Optional[str] = None,
) -> plt.Axes:
    """
    Plot catalytic reaction energy diagrams with intermediate plateaus and dashed connectors.

    Parameters
    ----------
    pathways : dict
        Mapping of {catalyst_name: [energy_step_0, energy_step_1, ...]}.
    step_labels : sequence of str
        Labels for each intermediate (e.g. ['* + CO + O2', '*CO + *O2', 'TS', '*CO2 + *', 'CO2 + *']).
    transition_states : dict, optional
        Mapping of {catalyst_name: [False, False, True, False, ...]}.
        If True for a step, that state is rendered with transition-state styling.
    relative_to_initial : bool, default True
        If True, offsets each profile so the initial intermediate is 0.0 eV.
    ax : plt.Axes, optional
        Matplotlib Axes to plot on.
    palette : str, default 'okabe_ito'
        Color palette identifier.
    plateau_width : float, default 0.45
        Width of horizontal energy levels along reaction coordinate.
    gap_width : float, default 0.35
        Horizontal distance between consecutive intermediate levels.
    energy_label, x_label : str
        Axis labels with units.
    title : str, optional
        Figure title.

    Returns
    -------
    plt.Axes
    """
    if ax is None:
        fig, ax = create_publication_figure(columns="double", aspect_ratio="2.2")

    n_steps = len(step_labels)
    step_centers = np.arange(n_steps) * (plateau_width + gap_width) + (plateau_width / 2.0)
    colors = get_palette(palette, len(pathways))

    for cat_idx, (cat_name, energies) in enumerate(pathways.items()):
        e_arr = np.asarray(energies, dtype=float)
        if len(e_arr) != n_steps:
            raise ValueError(
                f"Length of energies for '{cat_name}' ({len(e_arr)}) does not match "
                f"step labels count ({n_steps})."
            )

        if relative_to_initial:
            e_arr = e_arr - e_arr[0]

        color = colors[cat_idx]
        ts_flags = (
            transition_states.get(cat_name, [False] * n_steps)
            if transition_states
            else [False] * n_steps
        )

        for step_i in range(n_steps):
            x_center = step_centers[step_i]
            x_left = x_center - (plateau_width / 2.0)
            x_right = x_center + (plateau_width / 2.0)
            y_val = e_arr[step_i]
            is_ts = ts_flags[step_i]

            # Horizontal plateau
            ax.hlines(
                y=y_val,
                xmin=x_left,
                xmax=x_right,
                color=color,
                linewidth=2.0 if not is_ts else 1.8,
                linestyle="-" if not is_ts else "--",
                label=cat_name if (step_i == 0) else None,
                zorder=4,
            )

            # Connector to subsequent step
            if step_i < n_steps - 1:
                next_x_left = step_centers[step_i + 1] - (plateau_width / 2.0)
                next_y_val = e_arr[step_i + 1]
                ax.plot(
                    [x_right, next_x_left],
                    [y_val, next_y_val],
                    color=color,
                    linestyle=":",
                    linewidth=1.0,
                    alpha=0.85,
                    zorder=3,
                )

    ax.set_xticks(step_centers)
    ax.set_xticklabels(step_labels, rotation=0, ha="center")
    ax.set_xlim(-gap_width / 2.0, step_centers[-1] + plateau_width / 2.0 + gap_width / 2.0)
    ax.axhline(0, color="#CCCCCC", linewidth=0.6, linestyle="-", zorder=1)

    ax.set_ylabel(energy_label)
    ax.set_xlabel(x_label)
    if title:
        ax.set_title(title)

    ax.legend(loc="upper right", framealpha=0.9)
    return ax

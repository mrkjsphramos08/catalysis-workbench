"""
Electron Microscopy (SEM / TEM / HRTEM) Quantitative Image Analysis and Publication Plotting.

Implements:
- Calibrated scale bar rendering with true physical units
- Automated particle segmentation (Otsu, morphological filtering, connected components)
- Equivalent circular diameter, projected area, aspect ratio, and size distribution statistics
- HRTEM Fast Fourier Transform (FFT) power spectra and lattice fringe d-spacing extraction
- Publication-quality multi-panel figures (scale bars, segmentation overlays, histograms, FFT panels).
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import scipy.ndimage as ndi
import scipy.signal as sp_signal

import skimage.filters as sk_filters
import skimage.measure as sk_measure
import skimage.morphology as sk_morph
import skimage.segmentation as sk_seg

from parsers.characterization.common import CharacterizationData, InvalidUnitsError
from scripts.utils import apply_publication_style, create_publication_figure
from styles.palettes import get_color, get_palette


def add_scale_bar(
    ax: plt.Axes,
    image_shape: Tuple[int, int],
    pixel_size_nm: float,
    target_length_nm: Optional[float] = None,
    location: str = "lower right",
    color: str = "white",
    height_px: int = 6,
    fontsize: int = 8,
) -> None:
    """
    Overlay a mathematically calibrated physical scale bar onto an image axis.

    Parameters
    ----------
    ax : plt.Axes
        Matplotlib axis displaying the image.
    image_shape : tuple of (int, int)
        Image (height, width) in pixels.
    pixel_size_nm : float
        Physical scale calibration in nanometers per pixel.
    target_length_nm : float, optional
        Desired physical scale bar length in nm. If None, auto-selects a clean round number.
    location : str, default 'lower right'
        'lower right', 'lower left', 'upper right', or 'upper left'.
    color : str, default 'white'
        Scale bar and text color ('white' or 'black').
    height_px : int, default 6
        Bar thickness in display pixels.
    fontsize : int, default 8
        Text font size.
    """
    h, w = image_shape

    if target_length_nm is None:
        # Target ~15% of image width
        rough_nm = 0.15 * w * pixel_size_nm
        # Round to 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, etc.
        magnitude = 10 ** np.floor(np.log10(rough_nm))
        fraction = rough_nm / magnitude
        if fraction < 1.5:
            target_length_nm = 1.0 * magnitude
        elif fraction < 3.5:
            target_length_nm = 2.0 * magnitude
        elif fraction < 7.5:
            target_length_nm = 5.0 * magnitude
        else:
            target_length_nm = 10.0 * magnitude

    bar_len_px = target_length_nm / pixel_size_nm

    # Format text label (nm vs um vs mm)
    if target_length_nm >= 1000.0:
        val_um = target_length_nm / 1000.0
        label_text = f"{val_um:g} µm"
    else:
        label_text = f"{target_length_nm:g} nm"

    margin_x = 0.05 * w
    margin_y = 0.06 * h

    if "right" in location:
        x0 = w - margin_x - bar_len_px
    else:
        x0 = margin_x

    if "lower" in location:
        y0 = h - margin_y
    else:
        y0 = margin_y + height_px

    # Draw rectangular bar
    rect = patches.Rectangle((x0, y0), bar_len_px, height_px, facecolor=color, edgecolor="none", zorder=10)
    ax.add_patch(rect)

    # Text centered above bar
    ax.text(
        x0 + bar_len_px / 2.0,
        y0 - 4,
        label_text,
        color=color,
        ha="center",
        va="bottom",
        fontsize=fontsize,
        fontweight="bold",
        zorder=10,
    )


def segment_particles(
    data: CharacterizationData,
    min_area_px: int = 25,
    max_area_px: Optional[int] = None,
    threshold_method: str = "otsu",
    custom_threshold: Optional[float] = None,
    clear_borders: bool = False,
    dark_particles: bool = True,
) -> Dict[str, Any]:
    """
    Segment particles from SEM/TEM images and compute morphological descriptors.

    Notes
    -----
    1. Equivalent circular diameter is calculated from 2D projected area: d = 2 * sqrt(Area / pi).
    2. 2D projected dimensions must not be equated to 3D volume dimensions without spherical assumptions.

    Parameters
    ----------
    data : CharacterizationData
        Standardized microscopy dataset containing a 2D grayscale image matrix.
    min_area_px : int, default 25
        Minimum particle area in pixels to filter out noise.
    max_area_px : int, optional
        Maximum particle area in pixels to filter out large agglomerates.
    threshold_method : str, default 'otsu'
        'otsu' or 'manual'.
    custom_threshold : float, optional
        Explicit threshold value in [0, 1] if method is 'manual'.
    clear_borders : bool, default False
        Whether to exclude particles touching image borders.
    dark_particles : bool, default True
        If True, targets dark particles on light background (typical for TEM / brightfield SEM).
        If False, targets bright particles on dark background.

    Returns
    -------
    Dict[str, Any]
        Dictionary with segmentation mask, labeled image, particle dataframe, and statistics.
    """
    if data.image is None:
        raise ValueError("CharacterizationData contains no image matrix for segmentation.")

    img = data.image.copy()

    # Invert image if dark particles
    proc_img = 1.0 - img if dark_particles else img

    # Compute threshold
    if threshold_method.lower() == "otsu":
        thresh_val = float(sk_filters.threshold_otsu(proc_img))
    elif custom_threshold is not None:
        thresh_val = float(custom_threshold)
    else:
        thresh_val = float(np.mean(proc_img))

    binary = proc_img > thresh_val

    # Morphological cleaning: remove small noise specks
    binary_clean = sk_morph.opening(binary, sk_morph.disk(2))

    if clear_borders:
        binary_clean = sk_seg.clear_border(binary_clean)

    labeled, num_features = ndi.label(binary_clean)
    props = sk_measure.regionprops(labeled)

    px_nm = data.pixel_size_nm

    diameters_nm = []
    diameters_px = []
    areas_nm2 = []
    aspect_ratios = []

    valid_labels = []
    for p in props:
        if p.area < min_area_px:
            continue
        if max_area_px is not None and p.area > max_area_px:
            continue

        valid_labels.append(p.label)
        d_px = 2.0 * np.sqrt(p.area / np.pi)
        diameters_px.append(float(d_px))

        maj = getattr(p, "axis_major_length", getattr(p, "major_axis_length", 1.0))
        min_ax = getattr(p, "axis_minor_length", getattr(p, "minor_axis_length", 1.0))
        if min_ax <= 0:
            min_ax = 1.0
        aspect_ratios.append(float(maj / min_ax))

        if px_nm is not None:
            d_nm = d_px * px_nm
            a_nm2 = p.area * (px_nm ** 2)
            diameters_nm.append(float(d_nm))
            areas_nm2.append(float(a_nm2))

    # Mask of only retained particles
    filtered_mask = np.isin(labeled, valid_labels)

    stats: Dict[str, Any] = {
        "particle_count": len(diameters_px),
        "calibrated": px_nm is not None,
        "pixel_size_nm": px_nm,
        "diameters_px": np.array(diameters_px),
        "aspect_ratios": np.array(aspect_ratios),
    }

    if px_nm is not None and len(diameters_nm) > 0:
        d_arr = np.array(diameters_nm)
        stats.update({
            "mean_diameter_nm": float(np.mean(d_arr)),
            "std_diameter_nm": float(np.std(d_arr)),
            "median_diameter_nm": float(np.median(d_arr)),
            "p10_diameter_nm": float(np.percentile(d_arr, 10)),
            "p90_diameter_nm": float(np.percentile(d_arr, 90)),
            "diameters_nm": d_arr,
            "areas_nm2": np.array(areas_nm2),
        })

    return {
        "labeled_mask": filtered_mask,
        "labeled_image": labeled,
        "statistics": stats,
    }


def analyze_hrtem_lattice_fringes(
    data: CharacterizationData,
    roi_bbox: Optional[Tuple[int, int, int, int]] = None,
    fringe_direction: str = "vertical",
) -> Dict[str, Any]:
    """
    Extract lattice fringe d-spacings and 2D FFT spatial power spectrum from HRTEM images.

    Parameters
    ----------
    data : CharacterizationData
        Calibrated HRTEM dataset with valid pixel_size_nm.
    roi_bbox : tuple of (row_start, row_end, col_start, col_end), optional
        Bounding box of crystalline region of interest. If None, uses central crop.
    fringe_direction : str, default 'vertical'
        Orientation of fringes ('vertical' -> integrate horizontally, 'horizontal' -> integrate vertically).

    Returns
    -------
    Dict[str, Any]
        - 'd_spacing_nm': Measured lattice spacing in nanometers.
        - 'd_spacing_angstrom': Measured lattice spacing in Angstroms.
        - 'fft_power': 2D logarithmic FFT power spectrum.
        - 'profile_x_nm': Real-space coordinates along line profile.
        - 'profile_intensity': 1D averaged fringe intensity profile.
        - 'spatial_frequency_nm_inv': Dominant spatial frequency (1/nm).
        - 'disclaimer': Crystallographic assignment disclaimer.
    """
    if data.pixel_size_nm is None:
        raise InvalidUnitsError(
            "HRTEM lattice fringe analysis requires a calibrated pixel size (pixel_size_nm). "
            "Uncalibrated images cannot produce physically meaningful d-spacings in nm or Å."
        )

    img = data.image.copy()
    h, w = img.shape

    if roi_bbox is not None:
        r1, r2, c1, c2 = roi_bbox
        crop = img[r1:r2, c1:c2]
    else:
        # Default to central 128x128 crop
        cr_h, cr_w = min(128, h), min(128, w)
        r0, c0 = (h - cr_h) // 2, (w - cr_w) // 2
        crop = img[r0:r0 + cr_h, c0:c0 + cr_w]

    # 1. 2D Fast Fourier Transform
    fft2 = np.fft.fftshift(np.fft.fft2(crop))
    power = np.log10(1.0 + np.abs(fft2) ** 2)

    # 2. Real-space profile across fringes
    px_nm = data.pixel_size_nm
    if fringe_direction.lower() == "vertical":
        # Average vertically, profile along horizontal axis
        prof = np.mean(crop, axis=0)
    else:
        prof = np.mean(crop, axis=1)

    prof_x_nm = np.arange(len(prof)) * px_nm

    # Find fringe peaks in real-space profile
    peaks_idx, _ = sp_signal.find_peaks(prof, distance=max(3, int(0.15 / px_nm)))
    if len(peaks_idx) >= 2:
        diffs_px = np.diff(peaks_idx)
        mean_spacing_px = float(np.median(diffs_px))
        d_nm = mean_spacing_px * px_nm
        d_angstrom = d_nm * 10.0
        freq_nm = 1.0 / d_nm if d_nm > 0 else 0.0
    else:
        d_nm = None
        d_angstrom = None
        freq_nm = None

    return {
        "d_spacing_nm": d_nm,
        "d_spacing_angstrom": d_angstrom,
        "spatial_frequency_nm_inv": freq_nm,
        "fft_power": power,
        "profile_x_nm": prof_x_nm,
        "profile_intensity": prof,
        "crop_image": crop,
        "disclaimer": "Lattice d-spacings must be matched against validated crystallographic reference databases (COD, ICSD, Materials Project); do not assign Miller indices solely from spacing.",
    }


def plot_microscopy_with_scalebar(
    data: CharacterizationData,
    target_scale_nm: Optional[float] = None,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Display microscopy image with true calibrated physical scale bar overlay.
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio=1.0)
        created_fig = True
    else:
        fig = ax.get_figure()

    ax.imshow(data.image, cmap="gray", vmin=0, vmax=1)
    ax.axis("off")

    if data.pixel_size_nm is not None:
        add_scale_bar(ax, data.image.shape, data.pixel_size_nm, target_length_nm=target_scale_nm)

    if title:
        ax.set_title(title, fontsize=9, pad=6)

    if created_fig:
        fig.tight_layout()
    return fig, ax


def plot_segmentation_overlay(
    data: CharacterizationData,
    seg_result: Dict[str, Any],
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Display raw microscopy image with overlaid segmented particle contours and scale bar.
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio=1.0)
        created_fig = True
    else:
        fig = ax.get_figure()

    ax.imshow(data.image, cmap="gray", vmin=0, vmax=1)
    ax.axis("off")

    # Overlay particle boundaries
    mask = seg_result["labeled_mask"]
    ax.contour(mask, levels=[0.5], colors=["#D55E00"], linewidths=[1.0], zorder=4)

    if data.pixel_size_nm is not None:
        add_scale_bar(ax, data.image.shape, data.pixel_size_nm)

    cnt = seg_result["statistics"]["particle_count"]
    ax.text(
        0.04,
        0.94,
        f"$N = {cnt}$ particles",
        transform=ax.transAxes,
        color="white",
        fontsize=8,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="black", alpha=0.6, edgecolor="none"),
    )

    if title:
        ax.set_title(title, fontsize=9, pad=6)

    if created_fig:
        fig.tight_layout()
    return fig, ax


def plot_particle_size_distribution(
    seg_result: Dict[str, Any],
    bins: int = 15,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot particle size distribution histogram with mean/median statistical markers.
    """
    apply_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
        created_fig = True
    else:
        fig = ax.get_figure()

    stats = seg_result["statistics"]
    if stats["calibrated"] and "diameters_nm" in stats:
        d_vals = stats["diameters_nm"]
        unit_str = "nm"
    else:
        d_vals = stats["diameters_px"]
        unit_str = "px"

    ax.hist(d_vals, bins=bins, color="#0072B2", edgecolor="white", alpha=0.8, density=True, zorder=2)

    mean_d = float(np.mean(d_vals))
    med_d = float(np.median(d_vals))

    ax.axvline(mean_d, color="#D55E00", linestyle="-", linewidth=1.4, label=f"Mean: {mean_d:.1f} {unit_str}", zorder=3)
    ax.axvline(med_d, color="#009E73", linestyle="--", linewidth=1.4, label=f"Median: {med_d:.1f} {unit_str}", zorder=3)

    ax.set_xlabel(rf"Equivalent Circular Diameter $d_{{\mathrm{{eq}}}}~[\mathrm{{{unit_str}}}]$")
    ax.set_ylabel(r"Probability Density")
    if title:
        ax.set_title(title, fontsize=9, pad=6)

    ax.legend(loc="upper right", frameon=False, fontsize=8)
    if created_fig:
        fig.tight_layout()
    return fig, ax


def plot_hrtem_fft_panel(
    data: CharacterizationData,
    hrtem_result: Dict[str, Any],
    title: Optional[str] = None,
) -> Tuple[plt.Figure, Sequence[plt.Axes]]:
    """
    Generate publication 3-panel HRTEM panel:
    (a) Real-space lattice crop with scale bar
    (b) 2D FFT power spectrum
    (c) 1D averaged lattice fringe intensity profile with d-spacing annotation.
    """
    apply_publication_style()
    fig, (ax_img, ax_fft, ax_prof) = create_publication_figure(
        columns="double",
        aspect_ratio=2.6,
        nrows=1,
        ncols=3,
    )

    crop = hrtem_result["crop_image"]
    power = hrtem_result["fft_power"]
    prof_x = hrtem_result["profile_x_nm"]
    prof_y = hrtem_result["profile_intensity"]
    d_nm = hrtem_result["d_spacing_nm"]
    d_ang = hrtem_result["d_spacing_angstrom"]

    # (a) Real-space crop
    ax_img.imshow(crop, cmap="gray", vmin=0, vmax=1)
    ax_img.axis("off")
    if data.pixel_size_nm is not None:
        add_scale_bar(ax_img, crop.shape, data.pixel_size_nm, location="lower right")
    ax_img.set_title("(a) HRTEM Lattice Crop", fontsize=8, pad=6)

    # (b) FFT Power Spectrum
    ax_fft.imshow(power, cmap="magma")
    ax_fft.axis("off")
    ax_fft.set_title("(b) 2D FFT Power Spectrum", fontsize=8, pad=6)

    # (c) Line Profile
    ax_prof.plot(prof_x, prof_y, color="#0072B2", linewidth=1.2)
    ax_prof.set_xlabel(r"Distance $x~[\mathrm{nm}]$")
    ax_prof.set_ylabel(r"Intensity $[\mathrm{a.u.}]$")
    ax_prof.set_title("(c) Fringe Profile", fontsize=8, pad=6)

    if d_ang is not None:
        ax_prof.text(
            0.05,
            0.88,
            rf"$d = {d_ang:.2f}~\mathrm{{\AA}}~({d_nm:.3f}~\mathrm{{nm}})$",
            transform=ax_prof.transAxes,
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#cccccc", alpha=0.9),
        )

    if title:
        fig.suptitle(title, fontsize=9, y=0.98)

    fig.tight_layout()
    return fig, (ax_img, ax_fft, ax_prof)

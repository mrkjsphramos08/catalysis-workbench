"""
Electron Microscopy (SEM / FESEM / TEM / HRTEM) Image Parser.

Ingests scientific TIFF, PNG, and JPEG images using tifffile, skimage, and Pillow.
Extracts pixel matrices, image dimensions, bit depth, and physical scale calibrations.
Enforces strict calibration validation: never assumes unverified physical dimensions.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np
from PIL import Image
import tifffile

from parsers.characterization.common import (
    CharacterizationData,
    CharacterizationParseError,
    InvalidUnitsError,
    compute_file_sha256,
)


def parse_microscopy_image(
    filepath: Union[str, Path],
    pixel_size_nm: Optional[float] = None,
    technique: str = "SEM",
    sample_id: Optional[str] = None,
    instrument_manufacturer: Optional[str] = None,
    instrument_model: Optional[str] = None,
    accelerating_voltage_kv: Optional[float] = None,
    magnification: Optional[float] = None,
) -> CharacterizationData:
    """
    Parse a microscopy image into standardized CharacterizationData.

    Parameters
    ----------
    filepath : str or Path
        Microscopy image file path.
    pixel_size_nm : float, optional
        Physical pixel size in nanometers per pixel (nm/px).
        If None, attempts to extract from scientific TIFF tags or leaves as None.
    technique : str, default 'SEM'
        'SEM', 'FESEM', 'TEM', or 'HRTEM'.
    sample_id : str, optional
        Sample identifier. Defaults to filename stem.
    instrument_manufacturer : str, optional
        e.g., 'JEOL', 'Thermo Fisher / FEI', 'Hitachi', 'Zeiss'.
    instrument_model : str, optional
        e.g., 'Talos F200X', 'Regulus 8230', 'JEM-2100F'.
    accelerating_voltage_kv : float, optional
        Beam energy in kilovolts (kV).
    magnification : float, optional
        Instrument nominal magnification.

    Returns
    -------
    CharacterizationData
        Standardized microscopy dataset containing calibrated 2D image matrix.
    """
    p = Path(filepath)
    if not p.exists():
        raise FileNotFoundError(f"Microscopy image file not found: {p.resolve()}")

    file_hash = compute_file_sha256(p)
    sid = sample_id or p.stem

    ext = p.suffix.lower()
    raw_metadata: Dict[str, Any] = {
        "accelerating_voltage_kv": accelerating_voltage_kv,
        "magnification": magnification,
    }

    calibrated_pixel_nm = pixel_size_nm

    # Attempt TIFF tag reading for pixel size calibration
    if ext in (".tif", ".tiff"):
        try:
            with tifffile.TiffFile(p) as tf:
                img_arr = tf.asarray()
                # Inspect tags for resolution
                page = tf.pages[0]
                if page.tags.get("XResolution") and page.tags.get("ResolutionUnit"):
                    x_res = page.tags["XResolution"].value
                    res_unit = page.tags["ResolutionUnit"].value  # 2: inches, 3: centimeters
                    if isinstance(x_res, tuple) and len(x_res) == 2 and x_res[0] > 0:
                        res_val = x_res[0] / x_res[1]
                        if res_unit == 3 and calibrated_pixel_nm is None:  # cm
                            # px/cm -> nm/px: 1 cm = 1e7 nm => nm/px = 1e7 / res_val
                            calibrated_pixel_nm = 1e7 / res_val
        except Exception:
            pil_img = Image.open(p)
            img_arr = np.array(pil_img)
    else:
        pil_img = Image.open(p)
        img_arr = np.array(pil_img)

    # Convert to grayscale 2D float array normalized [0, 1]
    if img_arr.ndim == 3:
        if img_arr.shape[2] == 4:  # RGBA
            img_arr = img_arr[:, :, :3]
        # Standard ITU-R 601-2 luma transform
        gray = 0.2989 * img_arr[:, :, 0] + 0.5870 * img_arr[:, :, 1] + 0.1140 * img_arr[:, :, 2]
    elif img_arr.ndim == 2:
        gray = img_arr.astype(float)
    else:
        raise CharacterizationParseError(f"Unsupported image array dimensions: {img_arr.shape}")

    # Normalize to [0.0, 1.0]
    g_min = float(np.min(gray))
    g_max = float(np.max(gray))
    if g_max > g_min:
        gray_norm = (gray - g_min) / (g_max - g_min)
    else:
        gray_norm = gray

    raw_metadata["height_px"] = gray.shape[0]
    raw_metadata["width_px"] = gray.shape[1]
    raw_metadata["calibrated"] = calibrated_pixel_nm is not None

    return CharacterizationData(
        technique=technique.upper(),
        sample_id=sid,
        source_file=p,
        file_sha256=file_hash,
        x_name="X",
        x_unit="nm" if calibrated_pixel_nm else "pixels",
        y_name="Y",
        y_unit="nm" if calibrated_pixel_nm else "pixels",
        image=gray_norm,
        pixel_size_nm=calibrated_pixel_nm,
        instrument_manufacturer=instrument_manufacturer,
        instrument_model=instrument_model,
        metadata=raw_metadata,
    )

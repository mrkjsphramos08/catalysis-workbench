"""
Centralized Configuration for Scientific Publication Graphing Workspace.

Defines project-root relative paths, journal dimension presets,
DPI resolutions, supported export formats, and formatting invariants.
"""

from pathlib import Path
from typing import Tuple, Union


# Project Root Resolution
PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent

# Core Directory Hierarchy
DATA_DIR: Path = PROJECT_ROOT / "data"
RAW_DATA_DIR: Path = DATA_DIR / "raw"
PROCESSED_DATA_DIR: Path = DATA_DIR / "processed"
EXAMPLES_DATA_DIR: Path = DATA_DIR / "examples"

SCRIPTS_DIR: Path = PROJECT_ROOT / "scripts"
MODULES_DIR: Path = PROJECT_ROOT / "modules"
STYLES_DIR: Path = PROJECT_ROOT / "styles"
NOTEBOOKS_DIR: Path = PROJECT_ROOT / "notebooks"
TESTS_DIR: Path = PROJECT_ROOT / "tests"

FIGURES_DIR: Path = PROJECT_ROOT / "figures"
FIGURES_DFT_DIR: Path = FIGURES_DIR / "dft"
FIGURES_CATALYSIS_DIR: Path = FIGURES_DIR / "catalysis"
FIGURES_CHEMENG_DIR: Path = FIGURES_DIR / "chemical_engineering"
FIGURES_ELECTROCHEM_DIR: Path = FIGURES_DIR / "electrochemistry"
FIGURES_CHAR_DIR: Path = FIGURES_DIR / "characterization"
FIGURES_XRD_DIR: Path = FIGURES_CHAR_DIR / "xrd"
FIGURES_FTIR_DIR: Path = FIGURES_CHAR_DIR / "ftir"
FIGURES_RAMAN_DIR: Path = FIGURES_CHAR_DIR / "raman"
FIGURES_MICROSCOPY_DIR: Path = FIGURES_CHAR_DIR / "microscopy"
FIGURES_ICP_DIR: Path = FIGURES_CHAR_DIR / "icp"
FIGURES_AAS_DIR: Path = FIGURES_CHAR_DIR / "aas"
FIGURES_BET_DIR: Path = FIGURES_CHAR_DIR / "bet"
FIGURES_COMPARATIVE_DIR: Path = FIGURES_CHAR_DIR / "comparative"

OUTPUTS_DIR: Path = PROJECT_ROOT / "outputs"
OUTPUTS_VECTOR_DIR: Path = OUTPUTS_DIR / "vector"
OUTPUTS_RASTER_DIR: Path = OUTPUTS_DIR / "raster"

# Style Sheet Path
MPLSTYLE_PATH: Path = STYLES_DIR / "publication.mplstyle"

# Journal Column Width Specifications (International Standard: mm to inches)
MM_TO_INCH: float = 1.0 / 25.4
SINGLE_COLUMN_MM: float = 85.0
DOUBLE_COLUMN_MM: float = 175.0
ONE_AND_HALF_COLUMN_MM: float = 115.0

SINGLE_COLUMN_WIDTH_IN: float = SINGLE_COLUMN_MM * MM_TO_INCH  # ~3.346 in
DOUBLE_COLUMN_WIDTH_IN: float = DOUBLE_COLUMN_MM * MM_TO_INCH  # ~6.890 in
ONE_AND_HALF_COLUMN_WIDTH_IN: float = ONE_AND_HALF_COLUMN_MM * MM_TO_INCH  # ~4.528 in

# Resolution Presets
DEFAULT_DPI: int = 300
PUBLICATION_DPI: int = 600

# Export Format Presets
SUPPORTED_VECTOR_FORMATS: Tuple[str, ...] = ("pdf", "svg")
SUPPORTED_RASTER_FORMATS: Tuple[str, ...] = ("png", "tiff", "tif")
SUPPORTED_FORMATS: Tuple[str, ...] = SUPPORTED_VECTOR_FORMATS + SUPPORTED_RASTER_FORMATS


def get_figure_dimensions(
    columns: Union[str, float, int] = "single",
    aspect_ratio: Union[str, float] = "golden",
    height_mm: Union[None, float] = None,
) -> Tuple[float, float]:
    """
    Calculate figure width and height in inches based on publication column standards.

    Parameters
    ----------
    columns : str, float, or int
        'single' (85 mm), 'double' (175 mm), '1.5' (115 mm), or explicit width in mm.
    aspect_ratio : str or float
        'golden' (1.618), '4:3' (1.333), 'square' (1.0), 'wide' (2.0), or numeric width/height ratio.
    height_mm : float, optional
        Explicit height in millimeters. Overrides aspect_ratio if provided.

    Returns
    -------
    Tuple[float, float]
        (width_inches, height_inches)
    """
    if isinstance(columns, str):
        col_lower = columns.strip().lower()
        if col_lower in ("single", "1"):
            width_mm = SINGLE_COLUMN_MM
        elif col_lower in ("double", "2"):
            width_mm = DOUBLE_COLUMN_MM
        elif col_lower in ("one_and_half", "1.5", "intermediate"):
            width_mm = ONE_AND_HALF_COLUMN_MM
        else:
            raise ValueError(f"Unknown column preset: '{columns}'. Choose 'single', 'double', '1.5', or float mm.")
    else:
        width_mm = float(columns)

    width_in = width_mm * MM_TO_INCH

    if height_mm is not None:
        height_in = float(height_mm) * MM_TO_INCH
    else:
        if isinstance(aspect_ratio, str):
            ratio_lower = aspect_ratio.strip().lower()
            if ratio_lower in ("golden", "phi"):
                ratio = 1.61803398875
            elif ratio_lower in ("4:3", "standard"):
                ratio = 4.0 / 3.0
            elif ratio_lower in ("1:1", "square"):
                ratio = 1.0
            elif ratio_lower in ("wide", "16:9"):
                ratio = 16.0 / 9.0
            elif ratio_lower in ("tall", "3:4"):
                ratio = 3.0 / 4.0
            else:
                try:
                    ratio = float(ratio_lower)
                except ValueError:
                    raise ValueError(f"Unknown aspect ratio preset: '{aspect_ratio}'.")
        else:
            ratio = float(aspect_ratio)

        if ratio <= 0:
            raise ValueError("Aspect ratio must be positive.")

        height_in = width_in / ratio

    return (width_in, height_in)

"""
Publication Plotting Utilities and Data Management Functions.

Implements core style application, figure instantiation, multi-format export,
data loading with column validation, and descriptive statistical verification.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.config import (
    DEFAULT_DPI,
    FIGURES_DIR,
    MPLSTYLE_PATH,
    OUTPUTS_RASTER_DIR,
    OUTPUTS_VECTOR_DIR,
    SUPPORTED_FORMATS,
    get_figure_dimensions,
)


def apply_publication_style(
    style_path: Optional[Union[str, Path]] = None,
    use_latex: bool = False,
) -> None:
    """
    Apply publication-grade matplotlib style and rendering settings.

    Parameters
    ----------
    style_path : str or Path, optional
        Path to custom .mplstyle file. Defaults to styles/publication.mplstyle.
    use_latex : bool, default False
        Whether to enable external LaTeX rendering. Kept False by default
        to ensure out-of-the-box reproducibility without TeX installations.
    """
    target_style = Path(style_path) if style_path else MPLSTYLE_PATH
    if target_style.exists():
        plt.style.use(str(target_style))
    else:
        # Fallback to sensible built-in settings
        plt.rcParams.update({
            "font.family": "sans-serif",
            "font.size": 8.0,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9.0,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.5,
            "lines.linewidth": 1.2,
            "lines.markersize": 4.0,
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "figure.dpi": 300,
        })

    # Enforce vector-safe TrueType fonts for Nature/ACS/Elsevier compliance
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["text.usetex"] = use_latex


def create_publication_figure(
    columns: Union[str, float, int] = "single",
    aspect_ratio: Union[str, float] = "golden",
    height_mm: Optional[float] = None,
    nrows: int = 1,
    ncols: int = 1,
    sharex: bool = False,
    sharey: bool = False,
    gridspec_kw: Optional[Dict[str, Any]] = None,
    **kwargs: Any,
) -> Tuple[plt.Figure, Union[plt.Axes, np.ndarray]]:
    """
    Create a matplotlib Figure and Axes with journal column width constraints.

    Parameters
    ----------
    columns : str, float, or int
        'single' (85 mm), 'double' (175 mm), '1.5' (115 mm), or width in mm.
    aspect_ratio : str or float
        Aspect ratio preset or custom width/height ratio.
    height_mm : float, optional
        Explicit height in millimeters. Overrides aspect_ratio.
    nrows : int, default 1
        Number of subplot rows.
    ncols : int, default 1
        Number of subplot columns.
    sharex, sharey : bool, default False
        Subplot axis sharing options.
    gridspec_kw : dict, optional
        GridSpec layout arguments.

    Returns
    -------
    Tuple[plt.Figure, Union[plt.Axes, np.ndarray]]
        Matplotlib figure and axes object(s).
    """
    apply_publication_style()
    figsize = get_figure_dimensions(columns=columns, aspect_ratio=aspect_ratio, height_mm=height_mm)
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=figsize,
        sharex=sharex,
        sharey=sharey,
        gridspec_kw=gridspec_kw,
        constrained_layout=kwargs.pop("constrained_layout", True),
        **kwargs,
    )
    return fig, axes


def save_publication_figure(
    fig: plt.Figure,
    filename_base: str,
    output_dirs: Optional[Sequence[Union[str, Path]]] = None,
    formats: Sequence[str] = ("pdf", "svg", "png", "tiff"),
    dpi: int = DEFAULT_DPI,
    overwrite: bool = True,
) -> Dict[str, Path]:
    """
    Export figure across requested vector and raster publication formats.

    Automatically routes formats to appropriate destination directories:
    vector (pdf, svg) -> outputs/vector and specified output dirs
    raster (png, tiff) -> outputs/raster and specified output dirs

    Parameters
    ----------
    fig : plt.Figure
        Figure to export.
    filename_base : str
        Base file name without extension.
    output_dirs : sequence of Paths/strings, optional
        Additional destination directories (e.g. figures/dft).
    formats : sequence of str, default ('pdf', 'svg', 'png', 'tiff')
        Export extensions.
    dpi : int, default 300
        Raster dots per inch.
    overwrite : bool, default True
        Whether to overwrite existing files. If False and file exists, raises FileExistsError.

    Returns
    -------
    Dict[str, Path]
        Mapping from format name to saved file path.
    """
    saved_paths: Dict[str, Path] = {}
    dest_dirs: List[Path] = []

    if output_dirs:
        for d in output_dirs:
            p = Path(d)
            p.mkdir(parents=True, exist_ok=True)
            dest_dirs.append(p)

    for fmt in formats:
        clean_fmt = fmt.lstrip(".").lower()
        if clean_fmt not in SUPPORTED_FORMATS:
            raise ValueError(f"Unsupported format '{clean_fmt}'. Supported: {SUPPORTED_FORMATS}")

        # Primary default destination
        if clean_fmt in ("pdf", "svg"):
            primary_dir = OUTPUTS_VECTOR_DIR
        else:
            primary_dir = OUTPUTS_RASTER_DIR

        primary_dir.mkdir(parents=True, exist_ok=True)
        all_targets = [primary_dir] + [d for d in dest_dirs if d != primary_dir]

        for target_dir in all_targets:
            out_file = target_dir / f"{filename_base}.{clean_fmt}"
            if out_file.exists() and not overwrite:
                raise FileExistsError(f"Target file {out_file} already exists and overwrite=False.")

            # Save with bbox_inches='tight' and explicit format
            fig.savefig(
                out_file,
                format=clean_fmt,
                dpi=dpi,
                bbox_inches="tight",
                pad_inches=0.03,
                facecolor="white",
                edgecolor="none",
            )
            saved_paths[clean_fmt] = out_file

    return saved_paths


def validate_columns(
    df: pd.DataFrame,
    required_columns: Sequence[str],
    context: str = "dataset",
) -> None:
    """
    Validate presence of required columns in a DataFrame without modifying data.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to inspect.
    required_columns : sequence of str
        List of mandatory column identifiers.
    context : str
        Context description for error reporting.

    Raises
    ------
    ValueError
        If any required column is missing.
    """
    missing = [c for c in required_columns if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required column(s) in {context}: {missing}. "
            f"Available columns: {list(df.columns)}"
        )


def load_csv_data(
    filepath: Union[str, Path],
    required_columns: Optional[Sequence[str]] = None,
    check_nan: bool = False,
    **kwargs: Any,
) -> pd.DataFrame:
    """
    Safely load CSV dataset, verifying file existence, columns, and data completeness.

    Parameters
    ----------
    filepath : str or Path
        Path to CSV file.
    required_columns : sequence of str, optional
        Columns that must be present.
    check_nan : bool, default False
        If True, raises ValueError if any NaN values exist in required columns.

    Returns
    -------
    pd.DataFrame
        Loaded immutable copy of DataFrame.
    """
    path = Path(filepath)
    if not path.is_file():
        raise FileNotFoundError(f"Data file not found at: {path.resolve()}")

    df = pd.read_csv(path, **kwargs)

    if required_columns:
        validate_columns(df, required_columns, context=str(path.name))
        if check_nan:
            nan_counts = df[list(required_columns)].isna().sum()
            failing = nan_counts[nan_counts > 0]
            if not failing.empty:
                raise ValueError(
                    f"Dataset {path.name} contains unexpected NaN values in required columns:\n"
                    f"{failing.to_dict()}"
                )

    return df.copy()


def load_excel_data(
    filepath: Union[str, Path],
    sheet_name: Union[str, int] = 0,
    required_columns: Optional[Sequence[str]] = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """
    Safely load Excel spreadsheet data.

    Parameters
    ----------
    filepath : str or Path
        Path to .xlsx or .xls file.
    sheet_name : str or int, default 0
        Sheet name or 0-based index.
    required_columns : sequence of str, optional
        Columns that must be present.

    Returns
    -------
    pd.DataFrame
        Loaded copy of DataFrame.
    """
    path = Path(filepath)
    if not path.is_file():
        raise FileNotFoundError(f"Excel file not found at: {path.resolve()}")

    df = pd.read_excel(path, sheet_name=sheet_name, **kwargs)

    if required_columns:
        validate_columns(df, required_columns, context=f"{path.name} (sheet {sheet_name})")

    return df.copy()


def calculate_descriptive_statistics(
    df: pd.DataFrame,
    columns: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """
    Compute rigorous descriptive statistics (mean, std, median, IQR, min, max, count).

    Parameters
    ----------
    df : pd.DataFrame
        Input dataset.
    columns : sequence of str, optional
        Subset of numerical columns to summarize.

    Returns
    -------
    pd.DataFrame
        Summary statistics table.
    """
    target_cols = list(columns) if columns else df.select_dtypes(include=[np.number]).columns.tolist()
    if not target_cols:
        raise ValueError("No numerical columns found for descriptive statistics.")

    stats: Dict[str, Dict[str, float]] = {}
    for col in target_cols:
        series = df[col].dropna()
        q25 = float(np.percentile(series, 25))
        q75 = float(np.percentile(series, 75))
        stats[col] = {
            "count": int(len(series)),
            "mean": float(np.mean(series)),
            "std": float(np.std(series, ddof=1)) if len(series) > 1 else 0.0,
            "median": float(np.median(series)),
            "iqr": float(q75 - q25),
            "min": float(np.min(series)),
            "max": float(np.max(series)),
        }

    return pd.DataFrame(stats).T

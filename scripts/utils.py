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


def finalize_figure_layout(fig: plt.Figure) -> None:
    """Keep an existing layout engine; use tight layout only for unmanaged figures."""
    if fig.get_layout_engine() is None:
        fig.tight_layout()


def save_publication_figure(
    fig: plt.Figure,
    filename_base: str,
    output_dirs: Optional[Sequence[Union[str, Path]]] = None,
    formats: Sequence[str] = ("pdf", "svg", "png", "tiff"),
    dpi: int = DEFAULT_DPI,
    overwrite: bool = False,
    *,
    include_default_dirs: bool = True,
    bbox_inches: Optional[str] = None,
    close: bool = False,
) -> Dict[str, Path]:
    """Export figures with collision protection and stable physical dimensions.

    Existing positional arguments are preserved. output_dirs are additional
    destinations unless include_default_dirs=False. Equivalent paths/formats
    are deduplicated. Default exports preserve the canvas size; bbox_inches=
    'tight' explicitly opts into cropping. overwrite=True is required to replace
    files. All collisions/formats are checked and all formats rendered before
    writing; replacements use atomic per-file publication (not a batch transaction).
    Returned paths retain the last-destination convention. close=True closes the
    supplied figure even on failure; otherwise ownership remains with the caller.
    """
    from io import BytesIO
    import os
    import tempfile

    try:
        if not filename_base or Path(filename_base).name != filename_base or any(c in filename_base for c in '/\\:') or filename_base in ('.', '..'):
            raise ValueError("filename_base must be a plain filename without directories.")
        if not np.isfinite(dpi) or dpi <= 0:
            raise ValueError("dpi must be finite and positive.")
        if bbox_inches not in (None, "tight"):
            raise ValueError("bbox_inches must be None or 'tight'.")
        clean_formats = list(dict.fromkeys(fmt.lstrip(".").lower() for fmt in formats))
        if not clean_formats or any(fmt not in SUPPORTED_FORMATS for fmt in clean_formats):
            raise ValueError(f"Choose at least one supported format: {SUPPORTED_FORMATS}")
        targets = {}
        for fmt in clean_formats:
            default = OUTPUTS_VECTOR_DIR if fmt in ("pdf", "svg") else OUTPUTS_RASTER_DIR
            dirs = ([default] if include_default_dirs else []) + list(output_dirs or [])
            destinations = list(dict.fromkeys(Path(d).resolve() for d in dirs))
            if not destinations:
                raise ValueError("Supply output_dirs when include_default_dirs=False.")
            targets[fmt] = [d / f"{filename_base}.{fmt}" for d in destinations]
        for paths in targets.values():
            for path in paths:
                if path.exists() and (not overwrite or not path.is_file()):
                    raise FileExistsError(f"Output exists: {path}; choose a new name or explicitly set overwrite=True.")
        rendered = {}
        # rc override is necessary because savefig(None) otherwise inherits
        # savefig.bbox='tight' from the publication style.
        with mpl.rc_context({"savefig.bbox": None}):
            for fmt in clean_formats:
                buffer = BytesIO()
                fig.savefig(buffer, format=fmt, dpi=dpi, bbox_inches=bbox_inches,
                            pad_inches=0.03, facecolor="white", edgecolor="none")
                rendered[fmt] = buffer.getvalue()
        saved = {}
        for fmt, paths in targets.items():
            for path in paths:
                path.parent.mkdir(parents=True, exist_ok=True)
                if not overwrite:
                    # Exclusive creation also protects against a concurrent writer.
                    with path.open("xb") as stream:
                        stream.write(rendered[fmt])
                else:
                    temp_name = None
                    try:
                        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
                            temp_name = stream.name
                            stream.write(rendered[fmt])
                        os.replace(temp_name, path)
                    finally:
                        if temp_name is not None:
                            Path(temp_name).unlink(missing_ok=True)
                saved[fmt] = path
        return saved
    finally:
        if close:
            plt.close(fig)


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


def write_provenance_json(payload: Dict[str, Any], figure_name: str,
                          output_dirs: Sequence[Union[str, Path]], *, overwrite: bool = False) -> Path:
    """Write a truthful caller-built record; preflight all destinations before writing."""
    import json
    if not figure_name or Path(figure_name).name != figure_name or any(c in figure_name for c in '/\\:'):
        raise ValueError("figure_name must be a plain filename.")
    targets = list(dict.fromkeys(Path(d).resolve() / f"{figure_name}_metadata.json" for d in output_dirs))
    if not targets:
        raise ValueError("At least one provenance output directory is required.")
    content = json.dumps(payload, indent=2, allow_nan=False) + "\n"
    for target in targets:
        if target.exists() and (not overwrite or not target.is_file()):
            raise FileExistsError(f"Provenance output exists: {target}")
    for target in targets:
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w" if overwrite else "x", encoding="utf-8") as stream:
            stream.write(content)
    return targets[0]

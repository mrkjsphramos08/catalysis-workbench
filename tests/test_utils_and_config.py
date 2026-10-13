"""
Tests for scripts/config.py and scripts/utils.py.
"""

from pathlib import Path
import tempfile
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from scripts.config import (
    DOUBLE_COLUMN_MM,
    DOUBLE_COLUMN_WIDTH_IN,
    MM_TO_INCH,
    SINGLE_COLUMN_MM,
    SINGLE_COLUMN_WIDTH_IN,
    get_figure_dimensions,
)
from scripts.utils import (
    apply_publication_style,
    calculate_descriptive_statistics,
    create_publication_figure,
    load_csv_data,
    save_publication_figure,
    validate_columns,
)


def test_figure_dimensions_presets():
    """Verify single and double column dimension calculations."""
    w_single, h_single = get_figure_dimensions(columns="single", aspect_ratio="golden")
    assert pytest.approx(w_single, rel=1e-3) == SINGLE_COLUMN_WIDTH_IN
    assert pytest.approx(h_single, rel=1e-3) == SINGLE_COLUMN_WIDTH_IN / 1.61803398875

    w_double, h_double = get_figure_dimensions(columns="double", aspect_ratio="4:3")
    assert pytest.approx(w_double, rel=1e-3) == DOUBLE_COLUMN_WIDTH_IN
    assert pytest.approx(h_double, rel=1e-3) == DOUBLE_COLUMN_WIDTH_IN / (4.0 / 3.0)

    # Custom height
    w_cust, h_cust = get_figure_dimensions(columns=100.0, height_mm=50.0)
    assert pytest.approx(w_cust, rel=1e-3) == 100.0 * MM_TO_INCH
    assert pytest.approx(h_cust, rel=1e-3) == 50.0 * MM_TO_INCH


def test_figure_dimensions_invalid():
    """Verify invalid presets raise ValueError."""
    with pytest.raises(ValueError):
        get_figure_dimensions(columns="nonexistent")
    with pytest.raises(ValueError):
        get_figure_dimensions(columns="single", aspect_ratio="invalid_ratio")
    with pytest.raises(ValueError):
        get_figure_dimensions(columns="single", aspect_ratio=-1.5)


def test_create_and_save_publication_figure():
    """Verify figure creation and multi-format export (PDF, SVG, PNG, TIFF)."""
    fig, ax = create_publication_figure(columns="single", aspect_ratio="golden")
    ax.plot([0, 1, 2], [0, 1, 4], label="Test Curve")
    ax.set_xlabel("X (units)")
    ax.set_ylabel("Y (units)")

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        saved = save_publication_figure(
            fig=fig,
            filename_base="test_fig",
            output_dirs=[tmp_path],
            formats=("pdf", "svg", "png", "tiff"),
            dpi=150,
            overwrite=True,
        )

        assert "pdf" in saved and saved["pdf"].exists()
        assert "svg" in saved and saved["svg"].exists()
        assert "png" in saved and saved["png"].exists()
        assert "tiff" in saved and saved["tiff"].exists()

        # Check overwrite = False behavior
        with pytest.raises(FileExistsError):
            save_publication_figure(
                fig=fig,
                filename_base="test_fig",
                output_dirs=[tmp_path],
                formats=("png",),
                overwrite=False,
            )

    plt.close(fig)


def test_validate_columns():
    """Verify column presence check."""
    df = pd.DataFrame({"col_a": [1, 2], "col_b": [3, 4]})
    validate_columns(df, ["col_a", "col_b"])

    with pytest.raises(ValueError) as excinfo:
        validate_columns(df, ["col_a", "col_c"], context="test_df")
    assert "col_c" in str(excinfo.value)


def test_load_csv_data_nan_check():
    """Verify load_csv_data handles missing files and NaN detection."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        csv_file = Path(tmp_dir) / "test_data.csv"
        df_src = pd.DataFrame({"a": [1.0, 2.0, np.nan], "b": [4.0, 5.0, 6.0]})
        df_src.to_csv(csv_file, index=False)

        # Loading without check_nan succeeds
        loaded = load_csv_data(csv_file, required_columns=["a", "b"])
        assert len(loaded) == 3

        # Loading with check_nan raises ValueError
        with pytest.raises(ValueError):
            load_csv_data(csv_file, required_columns=["a", "b"], check_nan=True)

    with pytest.raises(FileNotFoundError):
        load_csv_data("non_existent_file.csv")


def test_calculate_descriptive_statistics():
    """Verify computation of summary statistics."""
    df = pd.DataFrame({"x": [1.0, 2.0, 3.0, 4.0, 5.0]})
    stats = calculate_descriptive_statistics(df, columns=["x"])
    assert stats.loc["x", "count"] == 5
    assert pytest.approx(stats.loc["x", "mean"]) == 3.0
    assert pytest.approx(stats.loc["x", "median"]) == 3.0
    assert pytest.approx(stats.loc["x", "min"]) == 1.0
    assert pytest.approx(stats.loc["x", "max"]) == 5.0

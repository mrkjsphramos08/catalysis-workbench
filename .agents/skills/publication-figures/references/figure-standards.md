# Figure implementation and review

## Existing API

```python
from scripts.utils import create_publication_figure, save_publication_figure
from styles.palettes import get_palette

fig, ax = create_publication_figure(columns="single")
# Add the scientifically validated data and labels using the existing domain API.
paths = save_publication_figure(
    fig, "unique_analysis_name", formats=("pdf", "svg", "png", "tiff"),
    dpi=300, overwrite=False,
)
```

`create_publication_figure` applies the shared style and accepts subplot settings, `height_mm`, and `aspect_ratio`. It returns an axes object for one panel or an array for multiple panels. Close figures after use, especially in batches/tests.

`save_publication_figure` defaults to `overwrite=False` and preflights all collisions before rendering/writing. Existing positional parameters are preserved. `output_dirs` adds to primary PDF/SVG `outputs/vector/` and raster `outputs/raster/` destinations; set `include_default_dirs=False` to write only explicit destinations. Equivalent destinations/formats are deduplicated. Returned paths retain the last-destination convention. Replacements with `overwrite=True` are atomic per file, not one transaction across all outputs. `close=True` releases the figure even on failure; otherwise the caller owns cleanup. The CLI closes its figures and offers `--output-root` and `--overwrite`. The helper does not generate provenance.

## Current defaults, not journal guarantees

| Property | Repository setting |
| --- | --- |
| Single / 1.5 / double column | 85 / 115 / 175 mm in `scripts/config.py` |
| Labels / ticks / legend | 8.5 / 7.5 / 7.5 pt in `publication.mplstyle` |
| Fonts | Arial with DejaVu Sans and other sans-serif fallbacks |
| Raster | 300 DPI default; 600 DPI publication preset available in config |
| Vector fonts | `pdf.fonttype=42`, `ps.fonttype=42`, `svg.fonttype="none"` applied by style utility |
| Colors | Existing Okabe-Ito/Tol palettes through `styles/palettes.py` |

PDF Type 42 preserves TrueType embedding where supported; SVG `fonttype="none"` preserves text but does not embed the fonts. Check recipient font availability or the journal's outlining requirements. See [Matplotlib fonts](https://matplotlib.org/stable/users/explain/text/fonts.html). LaTeX is optional and off by default.

Exports default to `bbox_inches=None` and override style cropping to preserve requested canvas dimensions. Explicit `bbox_inches="tight"` opts into changed physical bounds. Keep an existing constrained layout engine; use `finalize_figure_layout` instead of calling `tight_layout` on constrained figures. Inspect exported page size, embedded fonts/text, raster pixel dimensions and DPI, and all requested formats; do not infer compliance from the figure object's size alone. For line art, microscopy, and mixed panels use the journal's resolution rules.

## Scientific content checks

- Define error bars (SD/SE/CI), replicate count, fit interval and exclusions. Show fit residuals and distinguish observations from model curves. Avoid unsupported significant digits.
- Format scientific symbols and units consistently; verify minus signs, subscripts, Greek letters, and logarithmic-axis labels. Align legends and panel letters without reducing text below readable size.
- Use `modules.characterization.microscopy.add_scale_bar` with confirmed pixel calibration; a burned-in bar is not a substitute for validated sizing. Recheck scale after resizing. Disclose crops and contrast adjustments; retain raw images.
- Reuse `record_figure_provenance` in QE `adapters.py` or `record_characterization_provenance` in characterization `common.py` explicitly. See [provenance guidance](../../research-reproducibility/references/provenance-guidelines.md) for fields they do not collect automatically.

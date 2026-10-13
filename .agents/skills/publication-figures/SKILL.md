---
name: publication-figures
description: Creates and reviews publication figures and exports using the repository's Matplotlib utilities. Use for journal-ready plots, HER/LSV figures, fit diagnostics, SEM/TEM panels, annotations, or PDF/SVG/PNG/TIFF output.
---

# Publication figures

Reuse [figure standards](references/figure-standards.md) and the existing domain plot function before building a new figure.

1. Identify the scientific message, journal dimensions, required formats, physical units, reference scales, uncertainty, and intended panel layout. Apply the relevant domain skill; apply scientific-validation if deriving or changing quantities.
2. Reuse `styles/publication.mplstyle`, `styles/palettes.py`, `scripts.utils.create_publication_figure`, and `save_publication_figure`. Preserve shared styling; avoid a parallel theme/export framework.
3. Select `columns="single"`, `"1.5"`, `"double"`, or a numeric width in millimeters. Use journal-specific requirements when supplied; these project presets do not establish universal journal compliance.
4. Label quantity and units, potential/energy reference, normalization, spin channels, and uncertainty definitions. Use accessible palettes plus distinguishable line/marker styles, consistent legends, residuals for fitted data, and annotations that do not obscure measurements.
5. Preserve original data. Display offsets, normalization, clipping, smoothing, background subtraction, or contrast changes must be scientifically justified and disclosed. Never modify measurements solely to improve appearance.
6. Export the requested vector/raster formats and record provenance with the existing domain recorders. Inspect actual rendered figures at intended print size and enlarged view: labels, symbols, font substitution, clipping, contrast, error bars, legends, residuals, microscopy scale bars, and panel alignment. A completed script or nonempty file is not visual validation.
7. Report inspected files and remaining export/rendering limitations. If rendering is unavailable, explicitly mark visual review incomplete.

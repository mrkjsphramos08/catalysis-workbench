---
name: research-reproducibility
description: Protects original research datasets and makes analyses repeatable. Use for raw-data preservation, SHA-256 provenance, processing history, environment capture, synthetic fixtures, figure regeneration, or research-data sharing and Git review.
---

# Research reproducibility

Use [provenance guidelines](references/provenance-guidelines.md) when recording or reviewing an analysis.

1. Classify inputs as original experimental data, original simulation outputs, processed data, synthetic demonstrations, or reference datasets; generated figures are outputs, not measurements. Record origin, permission/license, and validation status independently of filenames.
2. Treat original datasets as immutable, including simulation restart files and instrument exports. Hash their exact bytes with SHA-256 before processing. Work in derived copies; never silently rewrite measurements or delete/move research files as cleanup.
3. Reuse existing parsers, `CharacterizationData.copy`, `add_processing_step`, `to_provenance_dict`, and domain figure provenance functions. Verify which fields and source files are actually recorded; complete missing processing/environment metadata explicitly.
4. Record input hashes, analysis parameters, scientific assumptions, unit/reference conversions, selection/exclusions, processing sequence, random seeds, software versions, code revision, and reproduction commands. Keep results and figures tied to that configuration.
5. Keep `data/raw`, `data/processed`, `data/examples`, `tests/fixtures`, `figures`, and `outputs` semantically distinct. Generators can overwrite examples/fixtures and exporters can overwrite figures; inspect targets and preserve originals before rerunning them.
6. Use deterministic synthetic data with recorded seeds and generator versions, label them as synthetic, and test expected scientific properties. A fixed seed alone does not guarantee identical bytes across software versions or different call orders.
7. Review `git status` and the scoped diff before delivery. Never commit confidential or unpublished research datasets without explicit permission. Check notebooks, metadata, figures, and embedded paths for private data as well as raw files. Do not upload data to external validation services without authorization.
8. Report the environment used and unresolved dependency drift, provenance gaps, and genuine-data coverage. Preserve an environment snapshot for a research result without dumping a machine-wide package inventory into shared project dependency files.

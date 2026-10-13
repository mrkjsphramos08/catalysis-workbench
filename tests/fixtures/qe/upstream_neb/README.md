# Authentic upstream Quantum ESPRESSO NEB reference

These unchanged output files come from QEF/q-e `NEB/examples/example01/reference` at commit `61569eb480231b649c47f631df9c5c83537df461`. The stdout identifies QE 6.0, dated 8 February 2017. No local QE calculation was run. This is authentic upstream simulation reference data, not a measurement and not a new synthetic fixture.

[provenance.json](provenance.json) records exact URLs and byte SHA-256 hashes. [QE-LICENSE.txt](QE-LICENSE.txt) copies the upstream GPL version 2 license; upstream copyright and attribution apply.

- `.dat`: seven images, coordinate / relative energy in eV / residual error. Image 4 has 0.2042674385 eV relative to image 1.
- `.int`: interpolated coordinate / energy pairs.
- `.out`: 13 optimization iterations, final forward/reverse barriers 0.204267 / 0.204268 eV; climbing image 4; convergence marker and JOB DONE.
- `.path`: authentic restart file used solely to verify explicit rejection.

Regression coverage does not certify every QE release, restart variant, or scientific suitability of arbitrary pathways. The older synthetic fixtures in the parent directory are preserved unchanged.

Repository `.gitattributes` disables line-ending conversion for the original reference files. Native stdout trailing spaces are retained and exempted only for that file from Git whitespace checks.

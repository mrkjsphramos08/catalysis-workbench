# QE formats and current support

| Executable | Existing API and checks | Official reference |
| --- | --- | --- |
| `pw.x` | `pw.parse_pw_output(filepath, strict=True, use_ase=True)` returns `QEPwOutput`; separate Ry/eV energies, completion/convergence, optional ASE geometry | [PW inputs and units](https://www.quantum-espresso.org/Doc/INPUT_PW.html) |
| `dos.x` | `dos.parse_dos_file`; non-spin energy/DOS/integrated DOS and collinear spin energy/up/down/integrated DOS | [DOS output notes](https://www.quantum-espresso.org/Doc/INPUT_DOS.html) |
| `projwfc.x` | `pdos.QEPdosCollection(directory, prefix=None, e_fermi=None)`; `get_projection`, `get_projections_for_plotting`, `calculate_d_band_center` | [Projection file layout and orbital order](https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html) |
| `bands.x` | `bands.parse_bands_gnu`; repeated path/energy blocks, optional spin-down file, explicit symmetry points; not arbitrary binary/XML band data | [Bands output options](https://www.quantum-espresso.org/Doc/INPUT_BANDS.html) |
| `neb.x` | `neb.parse_neb_output(dat_filepath, int_filepath=None, path_filepath=None, *, log_filepath=None, data_format="auto")`; native profiles and explicit stdout logs | [NEB inputs](https://www.quantum-espresso.org/Doc/INPUT_NEB.html), [actual output files](https://www.quantum-espresso.org/Doc/neb_user_guide/node7.html) |

## Projection checks

For non-spin `s/p/d`, expect energy, LDOS, then 1/3/5 components. For collinear spin, LDOS and each component have up/down pairs. QE's real-harmonic order is `pz, px, py` and `dz2, dzx, dzy, dx2-y2, dxy`. Match atom, element, wavefunction, and prefix from filenames; verify energy grids before combining files. A sum of projections is basis-dependent and is not automatically identical to total DOS. For d-band centers record the integration window, reference energy, spin selection, and zero-density behavior.

## NEB format and validation scope

Native `.dat` columns are reaction coordinate, relative electronic energy (eV), and residual error; `.int` contains coordinate/energy pairs. Coordinates are not assumed to be angstroms. Two-column profiles are supported. Legacy indexed tables need `data_format="legacy_indexed"` or an explicit `# Image Coord Energy...` header. No numeric column guessing is performed.

A sibling `.int` may be read automatically. Pass `log_filepath="run.out"` explicitly for final-iteration barriers, climbing image and convergence. Convergence requires a positive marker plus `JOB DONE.` in the final run. Without a log, barriers are discrete-image estimates and convergence/climbing-image status is unknown (`False`/`None`). Never confuse this with proof of nonconvergence. The old `path_filepath` argument remains a log alias, but actual `.path` restart files are rejected; they are not a documented analysis format.

Authentic coverage is the pinned upstream QE 6.0 H2+H example in `tests/fixtures/qe/upstream_neb/`, including byte hashes, upstream license and restart rejection. Other NEB versions and arbitrary systems remain unvalidated. Existing synthetic fixtures remain unchanged for legacy regression. See the [stabilization report](../../../../docs/audits/STABILIZATION_REPORT.md). Electronic-energy barriers are not free energies without thermodynamic corrections.

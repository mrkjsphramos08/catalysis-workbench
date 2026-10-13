# Instrument import decisions

| Input | Reuse | Required checks |
| --- | --- | --- |
| NOVA ASCII/CSV/Excel | `metrohm.parse_metrohm_export` | Header/delimiter/decimal separator, column units, current convention, exported settings; `.nox` is explicitly rejected with `MetrohmBinaryNoxError` |
| EC-Lab `.mpt` | `biologic.parse_biologic_mpt` | `Nb header lines`, encoding, technique/cycle columns, mA to A, and sign of exported `-Im(Z)` |
| EC-Lab `.mps` | `biologic.parse_biologic_mps` | Settings dictionary, not measurement arrays; do not route it as a time-series dataset |
| EC-Lab `.mpr` | `biologic.parse_biologic_mpr` via galvani | Installed galvani version, actual binary layout, and real-file evidence; text fixtures do not validate binary support |

Paths are inside `parsers/electrochemistry/`. The dispatcher uses extensions/headers when no instrument is supplied; parser errors propagate without retrying a different instrument. Use an explicit instrument when known. Reuse `column_mapping.py` and `custom_mapping` rather than introducing duplicate unit conversion rules.

Metrohm/BioLogic imports leave missing reference and temperature as `unknown` and `None`. Explicit caller values override unambiguous headers; temperature headers need K/C units. `raw_metadata["resolved_conditions"]` records origins. Supply `reference_electrode`, `temperature_k` and pH from confirmed conditions for RHE corrections. Unqualified Ag/AgCl and Hg/HgO are rejected for conversions; use the filling-solution variant. RHE identity does not require pH or temperature. Standard-potential table lookup retains its documented 298.15 K default; it is not a measured cell temperature.

Authoritative resources: [Metrohm NOVA export capabilities](https://www.metrohm-autolab.com/products/software/nova.html), [NOVA manual](https://www.metrohm-autolab.com/images/documentation/nova-user-manual.pdf), and [BioLogic technical note on MPT text export](https://my.biologic.net/documents/ec-lab-technical-notes-05-importing-an-ec-lab-text-file-into-excel-on-line/). Use the manual matching the instrument/software release; manufacturer export support does not prove this parser handles every option.

---
name: scientific-validation
description: Checks physical and numerical correctness. Use alongside domain skills for parser changes, equations, quantitative analysis or peak fitting, unit and reference-state audits, numerical regression, publication results, assumptions, and result discrepancies.
---

# Scientific validation

Treat this guidance as a priority for scientific changes, alongside the relevant domain skill. Use the [numerical validation checklist](references/numerical-validation-checklist.md) for detailed cases.

1. State the quantity, equation, units, sign convention, normalization, reference state, applicability, and expected behavior before changing code. Verify unfamiliar methods/constants against an authoritative source and the producing software version.
2. Inspect implementation and fixtures independently of historical claims. A parsed field, plausible plot, converged optimizer, or passing test does not establish physical correctness.
3. Check dimensional consistency end to end: input units, conversion count, numerical result, uncertainty, exported labels, and provenance. Compare compatible thermodynamic states and experimental conditions.
4. Construct an independent numerical example, analytical limit, or trusted external reference. Avoid deriving expected values by calling the same implementation or copying the same faulty equation. Justify tolerances from precision, conditioning, discretization, or measurement uncertainty.
5. Test known values, boundaries, malformed/missing data, and informative error/warning behavior. For scientifically significant or safety-critical implementation changes, require tests before accepting the implementation. Run the affected domain tests and full regression suite; investigate failures rather than adjusting expectations to make them pass.
6. Evaluate parameter/input uncertainty and sensitivity where they matter. Distinguish standard deviation, standard error, confidence intervals, prediction intervals, and systematic calibration/model uncertainty. Report unidentifiable parameters and unvalidated assumptions.
7. For figures inspect rendered output in addition to numerical tests. For datasets establish origin and genuine-data coverage with [research-reproducibility](../research-reproducibility/SKILL.md).

Report these five levels separately, marking gaps explicitly:

| Evidence level | Required statement |
| --- | --- |
| Execution | What command completed |
| Automated regression | Which tests passed/failed and warnings |
| Trusted numerical agreement | Independent value/source and error/tolerance |
| Scientific applicability | Assumptions checked and remaining uncertainty |
| Genuine-data validation | Actual experimental/simulation source and scope; otherwise not validated |

To maintain the skill collection itself, run `python .agents/skills/scientific-validation/scripts/validate_skills.py` from the repository root. This checks metadata and local links, not scientific truth or runtime skill activation. See [maintenance and discovery checks](../../../docs/AGENT_SKILLS_VALIDATION.md).

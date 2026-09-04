# Analysis code organisation

The question folders contain the primary cohort GLM, contrast, and group-summary
entry points. The final topography renderer is shared across all questions and is
kept as the single script `python/plot_q1_q10_topographies_mne.py`.

## Shared baseline implementation

`matlab/run_q1_procedure_cohort_estimation.m` is the shared Q1--Q4 implementation.
Its arguments select the phase and question, so it also generates Q2/Q3 Procedure
and DOM models:

```matlab
run_q1_procedure_cohort_estimation("procedure", "q2")
run_q1_procedure_cohort_estimation("dom", "q2")
run_q1_procedure_cohort_estimation("procedure", "q3")
run_q1_procedure_cohort_estimation("dom", "q3")
```

The convenience entries are provided in `matlab/q2/` and `matlab/q3/`.
Q2 and Q3 AOE use `matlab/run_q2_aoe_rebuild.m` and
`matlab/run_q3_aoe_rebuild.m`, respectively, because those models required the
corrected AOE event-label handling.

## Question folders

- `matlab/q1_aoe`, `q1_dom`, `q1_procedure`: Q1 cohort and group entries.
- `matlab/q2`, `q3`: Q2/Q3 Procedure/DOM convenience entries and AOE rebuilds.
- `matlab/q4_aoe`, `q4_dom`, `q4_procedure`: Q4 cohort and group entries.
- `q5_aoe`, `q5_dom`, `q5_procedure`: Q5 primary cohort/summary scripts.
- `q6_aoe`, `q6_dom`, `q6_procedure`: Q6 primary cohort/summary scripts.
- `q7_aoe`, `q7_dom`, `q7_procedure`: Q7 primary cohort/summary scripts.
- `q8_aoe`, `q8_dom`, `q8_procedure`: Q8 primary cohort/summary scripts.
- `matlab/q9`: corrected timeline-defined HILoad cohort and group analysis.
- `matlab/`: Q10 cohort estimation, contrast, and summary scripts.

All scripts require the private processed inputs and the local NIRS-SPM/SPM
toolboxes described in the root README; those data and toolboxes are not stored
in this public repository.

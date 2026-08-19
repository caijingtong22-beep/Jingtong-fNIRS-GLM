# Continuation guide: Q5–Q8 fNIRS GLMs

This document is the hand-off specification for continuing the current Q10
workflow with Q5–Q8.  It deliberately uses only direct condition comparisons;
it does **not** depend on resolving the meaning or coding of Felipe's
`global_baseline` regressor.

## Current repository and analysis state

- Repository root: `C:\aliceproject\Jingtong-fNIRS-GLM`
- Data and NIRS-SPM input files remain in the excluded source project:
  `C:\aliceproject\GLM`.
- Q10 is complete through first-level specification, estimation, contrast,
  individual/channel/ROI summaries, and figures.
- Q10 is `HILoadHE - HILoadLE` during the specified bleeding-related phase.
- Q10 first-level outputs use all 61 channels.  The current ROI summary uses
  all 45 physical long channels and retains QC-flagged channels.
- Short-channel regression has **not** been used.  Therefore these are
  reproducible exploratory results, not final neurophysiological inference.
- Outstanding Q10 work: NIRS-SPM spatial-map output and final reporting/
  interpretation.  Do not represent the present HbO/HbR contrasts as a clean
  canonical neural activation: their ROI-level contrast betas were frequently
  same-signed.

## Scientific conditions and coding

The 2 × 2 cells are:

| Cell label | Intrinsic load | Extraneous load |
|---|---:|---:|
| `LILE` | low | low |
| `HILE` | high | low |
| `LIHE` | low | high |
| `HIHE` | high | high |

The requested phases are analysed separately: **Procedure**, **DOM** and
**AOE**.  A phase-specific event name must be constructed only from event
labels that actually exist in the full event table.  Inspect and assert the
available labels before model creation; never guess onsets or reuse a phase
from another question.

## Governing modelling rule

Felipe's meeting recommendation was one small, purpose-built GLM per scientific
question/phase, rather than one large all-purpose design.  For a direct
condition comparison, use:

`X = [condition_A, condition_B, constant]`

where both event regressors are HRF-convolved by `specification_batch`; the
constant is added by NIRS-SPM.  The comparison is a contrast of the two task
columns, for example `[-1, +1, 0]` in the implemented order
`[condition_A, condition_B, constant]`.

No separate `global` regressor is required for these direct comparisons.  The
intercept/constant is not HRF-convolved.  Do **not** add a column representing
all remaining samples: together with an intercept and a target-condition
indicator it is perfectly collinear.

Always recover regressor indices from `SPM_nirs.xX.name` before constructing a
contrast.  Never assume column positions merely because they were expected.

## Q5–Q8 planned questions

The project-level mapping used in prior discussion is:

| Question | Factorial question | Preferred model type | Phase count |
|---|---|---|---:|
| Q5 | intrinsic-load main effect | 4-cell factorial contrast | 3 |
| Q6 | extraneous-load main effect | 4-cell factorial contrast | 3 |
| Q7 | intrinsic × extraneous interaction | 4-cell factorial contrast | 3 |
| Q8 | combined high/high versus low/low (`HIHE - LILE`) | 2-cell direct contrast | 3 |

Thus there are **12 first-level scientific models**: Q5–Q8 × Procedure/DOM/AOE,
each fitted separately for HbO and HbR.  “Model” here means a per-question,
per-phase specification per chromophore; it does not mean a different model
for every channel, since NIRS-SPM estimates channels together.

### Q5: intrinsic-load main effect

For a genuine main effect averaged over extraneous load, retain all four cells
for the relevant phase:

`X = [LILE, HILE, LIHE, HIHE, constant]`

Use:

`c_intrinsic = [-0.5, +0.5, -0.5, +0.5, 0]`

Positive beta means high intrinsic load > low intrinsic load, averaged equally
over low and high extraneous load.  Its reverse is the negative vector.

### Q6: extraneous-load main effect

Use the same four-cell model in the same implemented order:

`c_extraneous = [-0.5, -0.5, +0.5, +0.5, 0]`

Positive beta means high extraneous load > low extraneous load, averaged equally
over intrinsic load.

### Q7: interaction

Use the four-cell model:

`c_interaction = [+1, -1, -1, +1, 0]`

This is `(HIHE - HILE) - (LIHE - LILE)`: whether the extraneous-load effect is
larger at high intrinsic load.  Scaling by a common constant is immaterial to
the t statistic, but retain this stated vector consistently across all output
files and figure labels.

### Q8: combined high/high versus low/low

Use the smallest direct comparison:

`X = [LILE, HIHE, constant]`

`c_Q8 = [-1, +1, 0]`

Positive beta means `HIHE > LILE`.

## Important distinction: why Q5–Q7 have five columns

Felipe's "three regressors" rule applies directly to a **two-condition**
comparison, as implemented in Q10 and Q8.  A main effect averaged across the
other factor and an interaction cannot be recovered from only two cells without
changing their scientific meaning.  Q5–Q7 therefore need the four condition
regressors plus a constant.  This is still a small, explicit factorial model;
it avoids any unresolved global baseline definition.

If the supervisory team instead wants simple effects rather than factorial main
effects, make separate Q5/Q6 stratum-specific two-cell models and name them
explicitly (for example `HILE - LILE` at low extraneous load).  Do not label a
simple effect as a factorial main effect.

## Input and subject policy (copy Q10 unless a phase requires otherwise)

Use Q10 as the tested input policy:

- `subject001`, `subject002`, and `subject010`: V2 trial-1/full-timeline input
  under `first_level\nirs_spm_inputs_v2_fulltimeline_trial1`, with its
  corresponding selected events.  This avoids duplicate/repeated condition
  occurrences.
- Other subjects: `first_level\nirs_spm_inputs\subjectXXX` and the full event
  table at `amended_snirf\event_tables_full_annotations\event_table_all_subjects_full_annotations.csv`.
- Use raw-input channel metadata for all subjects, including V2 subjects.
- Before specification, assert the required cell × phase event occurs exactly
  once after the subject-specific selection policy.  If it does not, record the
  subject/phase as a controlled error; do not silently choose an occurrence.

## Implementation sequence for each question × phase

1. Copy `run_q10_cohort_specification.m` and rename it clearly, for example
   `run_q5_procedure_cohort_specification.m`.
2. Change `q10_conditions` into a generic event selector for the required
   phase and condition cells; write the selected-event CSV and
   `multiple_conditions.mat` per subject.
3. Set a question-specific output root such as
   `felipe_small_models/q5_intrinsic_procedure/cohort_specification`.
4. Keep the NIRS-SPM specification settings used in Q10 unless a documented
   methodological decision changes them: `hbo`/`hbr`, `DCT, 128`, `none`, and
   the same `specification_batch` call.
5. Update QA assertions:
   - Q5–Q7: expect 5 total columns and one constant.
   - Q8: expect 3 total columns and one constant.
   - calculate rank, task-regressor correlations and VIF; stop if rank is
     deficient, absolute correlation is at least 0.95, or VIF exceeds 10.
6. Run specification for a single subject first; inspect its design matrix and
   selected events.  Only then run the cohort.
7. Copy/adapt `run_q10_cohort_estimation.m`; update expected beta rows to 5
   for Q5–Q7 or 3 for Q8.  Estimate HbO and HbR separately.
8. Copy/adapt `run_q10_cohort_compute_contrast.m`.  Locate the named regressors
   in `SPM_nirs.xX.name`, build the contrast by name, and retain its existing
   channel-level beta/SE/t/p/BH-FDR output logic.
9. Copy/adapt the Q10 group and ROI summary scripts.  Preserve both result
   levels:
   - all 61 channels for the primary current first-level/group summary;
   - all 45 long channels, QC retained, for the anatomical ROI summaries.
10. Produce figures only after checking that the contrast direction in each
    title, legend and file name matches the implemented vector.

## Required output and quality checks

For every question/phase retain:

- event-selection CSV and condition MAT;
- design QA (columns, rank, correlations, VIF);
- estimation log and assertion of finite beta dimensions;
- individual channel contrast CSV (HbO/HbR) with beta, SE, t, p and BH-FDR;
- group channel summary (unweighted one-sample t test on subject-level betas;
  state clearly this is not a formal NIRS-SPM second-level model);
- long-channel ROI summaries, with 45 long channels and no QC exclusion;
- explicit contrast definition CSV and a README note of direction.

For every HbO/HbR interpretation, inspect direction jointly.  A typical
canonical increase would be HbO positive and HbR negative for the same
comparison; same-signed HbO/HbR values warrant cautious interpretation,
especially because short-channel regression is absent.

## Files to reuse

- `code/matlab/run_q10_cohort_specification.m`
- `code/matlab/run_q10_cohort_estimation.m`
- `code/matlab/run_q10_cohort_compute_contrast.m`
- `code/matlab/summarize_q10_group_contrasts.m`
- `code/matlab/summarize_q10_individual_contrasts.m`
- `code/matlab/summarize_q10_roi_long45_allqc.m`
- `code/matlab/summarize_q10_roi_cortical_system_long45_allqc.m`
- `code/matlab/plot_q10_roi_long45_allqc.m`
- `code/matlab/plot_q10_cortical_system_long45_allqc.m`

The new conversation should begin by reading this document, then inspecting
the original phase-specific event labels and running one subject/one phase
through specification QA before running any cohort estimation.

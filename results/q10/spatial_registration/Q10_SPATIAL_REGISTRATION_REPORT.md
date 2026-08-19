# Q10 spatial registration report

## Purpose

This registration is for NIRS-SPM visualisation of Q10 and later small-model contrasts. It does not define the GLM or its ROI statistics.

## Inputs

- Original allocation: `C:\aliceproject\GLM\channel_allocation.csv`
- Subject001 SNIRF measurement list: `C:\aliceproject\GLM\amended_snirf\Subj001_Sess101_run001.snirf`
- Subject001 SNIRF-derived probe geometry: `C:\aliceproject\GLM\first_level\nirs_spm_inputs\subject001\subject001_probe_geometry.mat`
- SPM template: `C:\aliceproject\GLM\toolboxes\spm\canonical\avg152T1.nii`

## Validated mapping

- Total channels: 61
- Long channels spatially rendered: 45
- Short channels not rendered: 16
- Unique SNIRF source-detector pairs: 61
- Validated long-channel pairs present in SNIRF: 45 / 45

## Channel 40 audit correction

The original allocation has `ch40:Src11-Det22` in both source and VC fields. That pair is absent from the subject001 SNIRF measurement list. The allocation’s established conversion pattern maps raw `Det22` to VC/SNIRF `Det06`; SNIRF contains `Src11-Det06`. This registration therefore uses `ch40:Src11-Det06`, without altering the original allocation file.

## Coordinate method

For each validated long channel, the coordinate is the midpoint of source and detector `sourcePos3D` / `detectorPos3D` from the subject001 SNIRF-derived probe geometry. Coordinates are projected to the cortical surface with `projection_CS`: 1. NIRS-SPM then renders the channel positions on the SPM canonical template.

## Limits

- This is MNI-compatible montage registration, not an atlas-based cortical label assignment.
- Short channels remain available to the GLM but are intentionally omitted from spatial rendering.
- Q10 maps made with this file are subject-level unless a separately specified group-level model is fitted.

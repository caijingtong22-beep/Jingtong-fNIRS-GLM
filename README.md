# Jingtong fNIRS GLM

Reproducible MATLAB/NIRS-SPM scripts and Q10 (HILoad HE minus HILoad LE) outputs for the Jingtong fNIRS project.

## Repository contents

- `code/matlab/`: Q10 first-level specification, estimation, contrast, summary, and plotting scripts.
- `metadata/channel_allocation.csv`: channel type and project ROI allocation.
- `results/q10/`: all Q10 CSV summaries generated in the current analysis, final figures, and interactive HbO-HbR quadrant plot.

## Deliberately excluded

Raw SNIRF files, participant-level NIRS-SPM MAT/SPM outputs, time-series data, recordings/transcripts, thesis files, and third-party toolboxes are not included.

## Current analytical scope

Q10 contrast is defined as HILoadHE minus HILoadLE. The present results retain all 45 physical long channels, including QC-flagged channels, and do not use short-channel regression. They should therefore be interpreted as current exploratory/reproducible results rather than a final neurophysiological inference.

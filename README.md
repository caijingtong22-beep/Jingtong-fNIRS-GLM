# Jingtong fNIRS GLM

Reproducible primary first-level GLM, group-analysis, and unified topography scripts for the Jingtong fNIRS project.

## Repository contents

- `code/matlab/`: Q1--Q10 primary first-level GLM and group-analysis scripts.
- `code/python/plot_q1_q10_topographies_mne.py`: the unified MNE-Python topography script used for the final Q1--Q10 maps.
- `metadata/channel_allocation.csv`: channel type and project ROI allocation.

## Deliberately excluded

Raw SNIRF files, participant-level NIRS-SPM MAT/SPM outputs, participant-level CSV files, group result tables, generated figures, time-series data, recordings/transcripts, thesis files, and third-party toolboxes are not included.

## Running the unified plotting script

The plotting script expects the locally retained group-level CSV files and the validated channel-coordinate CSV. Set the environment variable `FNIRS_REVIEW_ROOT` to the private review-data directory before running it. These data are intentionally not stored in this public repository.

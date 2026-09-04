"""Q1-Q10 primary-result topomaps using Jess's MNE topomap framework.

Only the plotting/interpolation mechanism is borrowed.  Channel geometry,
ROI allocation, beta values, and p-values are all from this project.
"""
import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
import mne
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


THIS_FILE = Path(__file__).resolve()
# The data and registration roots are configurable so that this public script
# does not depend on the author's private local directory layout.  Set
# FNIRS_REVIEW_ROOT to the local review-data directory before running.
REVIEW_ROOT = Path(os.environ.get("FNIRS_REVIEW_ROOT", THIS_FILE.parents[3]))
QUESTION_SPECS = {
    1: {"label": "LILE - ExperimentalBaseline", "folder": "q1_lile_vs_experimentalbaseline", "stem": "q1_lile_minus_experimentalbaseline", "n": 29},
    2: {"label": "HILE - ExperimentalBaseline", "folder": "q2_hile_vs_experimentalbaseline", "stem": "q2_hile_minus_experimentalbaseline", "n": 29},
    3: {"label": "LIHE - ExperimentalBaseline", "folder": "q3_lihe_vs_experimentalbaseline", "stem": "q3_lihe_minus_experimentalbaseline", "n": 29},
    4: {"label": "HIHE - ExperimentalBaseline", "folder": "q4_hihe_vs_experimentalbaseline", "stem": "q4_hihe_minus_experimentalbaseline", "n": 29},
    5: {"label": "High intrinsic - low intrinsic", "folder": "q5_intrinsic", "stem": "q5_high_intrinsic_minus_low_intrinsic", "n": 30},
    6: {"label": "High extraneous - low extraneous", "folder": "q6_extraneous", "stem": "q6_high_extraneous_minus_low_extraneous", "n": 30},
    7: {"label": "Interaction", "folder": "q7_interaction", "stem": "q7_interaction", "n": 30},
    8: {"label": "HIHE - LILE", "folder": "q8_hihe_vs_lile", "stem": "q8_HIHE_minus_LILE", "n": 30},
}
STAGES = ("procedure", "dom", "aoe")
# Q9 is the corrected single timeline-defined HILoad contrast; Q10 is also
# a single direct contrast.  Both are mapped with the same 45 registered
# long channels as Q1-Q8.
SPECIAL_INPUTS = {
    9: {
        "label": "HILoad - ExperimentalBaseline",
        "files": {
            "": REVIEW_ROOT / "02_group_level_results" /
            "q9_hiload_vs_experimentalbaseline_group_channel_summary_n29.csv",
        },
    },
    10: {
        "label": "HE - LE",
        "files": {
            "": REVIEW_ROOT / "02_group_level_results" /
            "q10__q10_HE_minus_LE_group_channel_summary_n30.csv",
        },
    },
}
REGISTRATION_FILE = REVIEW_ROOT / "spatial mapping" / "q10_validated_channel_mni_coordinates.csv"
OUT_DIR = REVIEW_ROOT / "04_figures_main" / "mne_framework_trial"

ROI_COLOURS = {
    "Medial Frontopolar PFC": "#0059B3",
    "Lateral Frontopolar PFC": "#D95614",
    "DL PFC": "#7326A8",
    "Dorsolateral PFC": "#7326A8",
    "VL PFC": "#1A8C75",
    "Ventrolateral PFC": "#1A8C75",
    "Medial PFC": "#BE3158",
    "PMC": "#1496CC",
    "Premotor Cortex": "#1496CC",
    "SMA": "#947000",
    "Supplementary Motor Area": "#947000",
    "M1": "#4D4D4D",
    "Primary Motor Cortex": "#4D4D4D",
}
ROI_SHORT = {
    "Medial Frontopolar PFC": "mFPC", "Lateral Frontopolar PFC": "lFPFC",
    "DL PFC": "DL PFC", "Dorsolateral PFC": "DL PFC",
    "VL PFC": "VL PFC", "Ventrolateral PFC": "VL PFC",
    "Medial PFC": "mPFC", "PMC": "PMC", "Premotor Cortex": "PMC",
    "SMA": "SMA", "Supplementary Motor Area": "SMA",
    "M1": "M1", "Primary Motor Cortex": "M1",
}
PROJECT_CMAP = LinearSegmentedColormap.from_list(
    "project_blue_white_red", ["#1440FF", "#FFFFFF", "#FF1A1A"], N=256
)


def main(questions: tuple[int, ...] | None = None) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if questions is None:
        questions = tuple(QUESTION_SPECS) + tuple(SPECIAL_INPUTS)
    known_questions = set(QUESTION_SPECS) | set(SPECIAL_INPUTS)
    unknown = set(questions) - known_questions
    if unknown:
        raise ValueError(f"Unknown question number(s): {sorted(unknown)}")
    registration = pd.read_csv(REGISTRATION_FILE)
    registration = registration.loc[
        registration["channel_type"].eq("Long") & registration["rendered_geometry"].eq(1),
        ["channel_id", "roi", "midpoint_2d_x", "midpoint_2d_y"],
    ].copy()
    for question in questions:
        if question in QUESTION_SPECS:
            spec = QUESTION_SPECS[question]
            stage_files = {stage: group_input_file(spec, stage) for stage in STAGES}
        else:
            spec = SPECIAL_INPUTS[question]
            stage_files = spec["files"]
        for stage, beta_file in stage_files.items():
            if not beta_file.is_file():
                raise FileNotFoundError(beta_file)
            stage_beta = pd.read_csv(beta_file)
            stage_beta = stage_beta.loc[stage_beta["channel_type"].eq("Long")].copy()
            stage_limit = max(1.0, np.ceil(stage_beta["mean_contrast_beta"].abs().max() * 2) / 2)
            for hemoglobin, display_name in (("HbO", "HbO2"), ("HbR", "HbR")):
                beta = stage_beta.loc[stage_beta["hemoglobin"].eq(hemoglobin)].copy()
                plot_one(beta, registration, question, spec["label"], stage, display_name, stage_limit)


def group_input_file(spec: dict, stage: str) -> Path:
    folder = f"{spec['folder']}_{stage}__cohort_specification__"
    stem = f"{spec['stem']}_{stage}_group_channel_summary_n{spec['n']}.csv"
    file = REVIEW_ROOT / "02_group_level_results" / f"{folder}{stem}"
    if not file.is_file():
        raise FileNotFoundError(file)
    return file


def plot_one(beta: pd.DataFrame, registration: pd.DataFrame, question: int, contrast_label: str,
             stage: str, display_name: str, limit: float) -> None:
    data = beta.merge(registration, on="channel_id", how="inner", validate="one_to_one",
                      suffixes=("_beta", ""))
    if len(data) != 45:
        raise ValueError(f"Expected 45 registered long channels, found {len(data)}")

    # Midpoints are the project's own physical probe coordinates.  They are
    # converted from centimetre-scale geometry to metres for MNE only; their
    # relative geometry is not changed or reordered.
    pos = data[["midpoint_2d_x", "midpoint_2d_y"]].to_numpy(dtype=float) / 100.0
    values = data["mean_contrast_beta"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(10, 10), constrained_layout=True)
    image, _ = mne.viz.plot_topomap(
        values, pos, axes=ax, show=False, sensors=False, contours=0,
        image_interp="cubic", extrapolate="head", border="mean",
        cmap=PROJECT_CMAP, vlim=(-limit, limit),
    )
    xx, yy, alpha = fade_unsupported_interpolation(image, pos)
    draw_supported_contours(ax, image, xx, yy, alpha, limit)
    fig.colorbar(image, ax=ax, shrink=0.78, pad=0.04, label="Mean contrast beta (µM)")

    # Project-specific overlays, rather than Jess's channel/ROI tables.
    for roi, subset in data.groupby("roi", sort=False):
        colour = ROI_COLOURS.get(roi, "#333333")
        xy = subset[["midpoint_2d_x", "midpoint_2d_y"]].to_numpy(dtype=float) / 100.0
        ax.scatter(xy[:, 0], xy[:, 1], s=185, facecolors="none", edgecolors=colour,
                   linewidths=2.1, zorder=5)
        centre = xy.mean(axis=0)
        ax.text(*centre, ROI_SHORT.get(roi, roi), ha="center", va="center", color=colour,
                fontsize=10, fontweight="bold", zorder=7,
                path_effects=[patheffects.withStroke(linewidth=3, foreground="white")])

    ax.scatter(pos[:, 0], pos[:, 1], s=145, c=values, cmap=PROJECT_CMAP, vmin=-limit, vmax=limit,
               edgecolors="black", linewidths=0.85, zorder=6)
    p_column = next((column for column in ("p_two_sided", "group_p_two_sided", "p_value")
                     if column in data.columns), None)
    if p_column is None:
        raise ValueError("No two-sided group p-value column was found.")
    significant = data[p_column].to_numpy(dtype=float) < 0.05
    ax.scatter(pos[significant, 0], pos[significant, 1], marker="*", s=260, c="black",
               linewidths=0.6, zorder=8)
    for xy, channel_id in zip(pos, data["channel_id"].astype(int)):
        ax.annotate(str(channel_id), xy, xytext=(0, -16), textcoords="offset points",
                    ha="center", va="top", fontsize=9, fontweight="bold", color="black", zorder=9)

    title_prefix = f"Q{question} ({stage.upper()})" if stage else f"Q{question}"
    ax.set_title(f"{title_prefix}: {contrast_label}, {display_name}",
                 fontsize=13, fontweight="bold", pad=15)
    ax.set_axis_off()
    stage_stem = f"_{stage}" if stage else ""
    stem = f"q{question}{stage_stem}_{display_name.lower()}_mne_own_registration_fade_to_white"
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_DIR / f"{stem}.pdf",
                bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fade_unsupported_interpolation(image, pos: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fade MNE's head extrapolation to white beyond measured-channel support.

    The head outline and channel placement remain MNE's standard topomap view,
    so M1 is not pushed downward by a tight custom outline.  Only the colour
    opacity is attenuated where the nearest physical channel is too distant.
    """
    values = np.asarray(image.get_array())
    height, width = values.shape[-2:]
    x0, x1, y0, y1 = image.get_extent()
    xx, yy = np.meshgrid(np.linspace(x0, x1, width), np.linspace(y0, y1, height))
    distances = np.sqrt((xx[..., None] - pos[:, 0]) ** 2 + (yy[..., None] - pos[:, 1]) ** 2)
    nearest = distances.min(axis=2)
    pairwise = np.sqrt(((pos[:, None, :] - pos[None, :, :]) ** 2).sum(axis=2))
    pairwise[pairwise == 0] = np.inf
    local_spacing = np.median(pairwise.min(axis=1))
    support = np.exp(-(nearest ** 2) / (2 * (2.0 * local_spacing) ** 2))
    alpha = np.clip((support - 0.04) / 0.76, 0, 1)
    image.set_alpha(alpha)
    return xx, yy, alpha


def draw_supported_contours(ax, image, xx: np.ndarray, yy: np.ndarray,
                            alpha: np.ndarray, limit: float) -> None:
    """Add contour lines only where the interpolated field has support."""
    values = np.asarray(image.get_array(), dtype=float)
    # MNE's default head outline is the 9.5-cm sphere centred at (0, 0).
    # Masking this explicitly prevents contour lines from appearing outside it.
    inside_head = xx ** 2 + yy ** 2 <= 0.095 ** 2
    masked = np.ma.masked_where((alpha < 0.14) | (~inside_head), values)
    ax.contour(xx, yy, masked, levels=np.linspace(-limit, limit, 21),
               colors="#202020", linewidths=0.45, alpha=0.65, zorder=3)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render Q1-Q10 MNE topomaps.")
    parser.add_argument("--questions", nargs="+", type=int, choices=tuple(range(1, 11)),
                        default=tuple(range(1, 11)), help="Question numbers to render (default: 1 through 10).")
    args = parser.parse_args()
    main(tuple(args.questions))

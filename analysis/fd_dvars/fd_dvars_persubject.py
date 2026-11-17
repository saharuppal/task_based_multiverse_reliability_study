#!/bin/bash
# """
# ********** per subjet
# Compute censor-excluded FD-DVARS Spearman's rho for ALL TWELVE
# denoising pipelines.
#
# Spike columns are identified from confoundevs.txt where present.
# For non-spike-reg pipelines, no one-hot spike column exists, so zero
# frames get censored and the "corrected" rho equals the rho over all
# available (non-NaN) frames 
#
# Expected layout:
# further_analysis/fd_dvars_dof/
#     <task>/
#         sub-0XX/
#             ses-X/
#                 fd_raw.txt
#                 <pipeline>/
#                     dvars_raw.txt
#
# Outputs ONE file:
#   fd_dvars_per_subject.csv   — one row per subject x task x pipeline,
#                                 sessions combined via a PLAIN AVERAGE
#                                 of the two session-level rho values
#                                 (not pooled raw FD/DVARS values).
#                                 Session-level results are computed as
#                                 an in-memory intermediary only and are
#                                 not written.
# """

module load FSL/6.0.7.17
source $FSLDIR/etc/fslconf/fsl.sh
module load SciPy-bundle/2025.06-gfbf-2025a

python3 - <<'PY'

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE_DIR = Path("further_analysis/fd_dvars_dof")
OUT_SUBJECT = BASE_DIR / "fd_dvars_per_subject.csv"
FEAT_ROOT = Path("first_levels")

def build_design_path(task: str, subject: str, session: str, pipeline: str) -> Path:
    feat_dir = (
        FEAT_ROOT / task / subject / session
        / f"{subject}_{session}_{task}_{pipeline}.feat"
    )
    return feat_dir / "confoundevs.txt"   # was: design.mat

# ----------------------------------------------------
# TASKS AND PIPELINES
# ----------------------------------------------------

TASKS = ["nback", "gonogo", "mid", "oasis"]

PIPELINES = [
    "6HMP", "24HMP",
    "6HMP_WM_CSF", "24HMP_WM_CSF",
    "6HMP_aCompCor", "24HMP_aCompCor",
    "6HMP_spike_reg", "24HMP_spike_reg",
    "6HMP_WM_CSF_spike_reg", "24HMP_WM_CSF_spike_reg",
    "6HMP_aCompCor_spike_reg", "24HMP_aCompCor_spike_reg",
]


# ----------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------

def load_1d(path: Path):
    """Load a whitespace/newline-delimited 1D text file as a float array.
    'n/a' and 'nan' both become NaN. Values keep their original temporal
    position."""
    arr = np.genfromtxt(
        path, dtype=float,
        missing_values=["n/a", "nan"], filling_values=np.nan,
    )
    return np.atleast_1d(arr).astype(float)


def read_confound_matrix(path: Path) -> np.ndarray:
    """Parse the plain whitespace-delimited confoundevs.txt matrix."""
    return np.loadtxt(path)


def get_spike_censored_indices(design_path: Path):
    """
    Identify spike/outlier regressors in confoundevs.txt.

    Spike columns here are raw, un-demeaned one-hot columns: exactly
    two distinct values (0.0 and 1.0), with 1.0 occurring exactly once.
    Continuous confound columns (motion params, WM, CSF, aCompCor) have
    many distinct values and are ignored -- for non-spike
    pipelines this returns an empty set (not an error).
    """
    if not design_path.exists():
        return None

    try:
        design = read_confound_matrix(design_path)
    except Exception as e:
        warnings.warn(f"Could not parse confoundevs.txt: {design_path}\nError: {e}")
        return None

    censored = set()

    for col in range(design.shape[1]):
        values, counts = np.unique(design[:, col], return_counts=True)

        if len(values) == 2 and set(values) == {0.0, 1.0} and np.sum(counts == 1) == 1:
            singleton_idx = np.flatnonzero(counts == 1)[0]
            spike_value = values[singleton_idx]
            spike_row = np.flatnonzero(design[:, col] == spike_value)[0]
            censored.add(int(spike_row))

    return censored


def align_fd_dvars(fd: np.ndarray, dvars: np.ndarray):
    """
    Align FD and DVARS timepoint-for-timepoint. Returns aligned arrays
    PLUS the original (full-run, 0-indexed) timepoint each aligned
    position corresponds to -- needed so censored-frame exclusion
    (identified in original design.mat indexing) can be applied
    correctly even after any first-frame trimming here.

    If lengths differ by more than one frame, alignment is not
    attempted because this could indicate actual volume removal
    somewhere within run.
    """
    n_fd = len(fd)
    n_dv = len(dvars)

    if n_fd == n_dv:
        fd_al, dv_al = fd, dvars
        orig_idx = np.arange(n_fd)
        note = "matched"

    elif abs(n_fd - n_dv) == 1:
        if n_fd > n_dv:
            fd_al, dv_al = fd[1:], dvars
            orig_idx = np.arange(1, n_fd)
            note = "trimmed_fd_first_frame"
        else:
            fd_al, dv_al = fd, dvars[1:]
            orig_idx = np.arange(1, n_dv)
            note = "trimmed_dvars_first_frame"

    else:
        return None, None, None, f"LENGTH_MISMATCH_fd{n_fd}_dvars{n_dv}", None

    mask = np.isfinite(fd_al) & np.isfinite(dv_al)
    n_dropped_nan = int((~mask).sum())
    fd_al, dv_al, orig_idx = fd_al[mask], dv_al[mask], orig_idx[mask]

    if n_dropped_nan:
        note += f"+dropped_{n_dropped_nan}_nan_pairs"

    return fd_al, dv_al, orig_idx, note, len(fd_al)


def compute_rho(fd_al, dv_al, note_suffix=""):
    if fd_al is None or len(fd_al) < 2:
        return np.nan, np.nan, (note_suffix + "+too_few_observations" if fd_al is not None else note_suffix)
    rho, pval = spearmanr(fd_al, dv_al)
    return rho, pval, note_suffix


# ----------------------------------------------------
# MAIN ANALYSIS
# ----------------------------------------------------

def main():
    rows = []
    missing_fd = []
    missing_dvars = []
    missing_design = []

    if not BASE_DIR.is_dir():
        raise FileNotFoundError(f"BASE_DIR does not exist:\n{BASE_DIR.resolve()}")

    print("=" * 70)
    print("FD-DVARS SPEARMAN CORRELATION (with censored-frame correction)")
    print("=" * 70)
    print(f"Base directory:      {BASE_DIR.resolve()}")
    print(f"Per-subject output:  {OUT_SUBJECT.resolve()}")
    print()

    for task in TASKS:
        task_dir = BASE_DIR / task
        if not task_dir.is_dir():
            warnings.warn(f"Task directory not found, skipping: {task_dir}")
            continue

        print(f"Processing task: {task}")

        for sub_dir in sorted(p for p in task_dir.glob("sub-*") if p.is_dir()):
            subject = sub_dir.name

            for ses_dir in sorted(p for p in sub_dir.glob("ses-*") if p.is_dir()):
                session = ses_dir.name

                fd_path = ses_dir / "fd_raw.txt"
                if not fd_path.exists():
                    missing_fd.append(str(fd_path))
                    continue
                try:
                    fd = load_1d(fd_path)
                except Exception as e:
                    warnings.warn(f"Could not load FD file: {fd_path}\nError: {e}")
                    continue

                for pipeline in PIPELINES:  # all 12 pipelines
                    dvars_path = ses_dir / pipeline / "dvars_raw.txt"
                    if not dvars_path.exists():
                        missing_dvars.append(str(dvars_path))
                        continue
                    try:
                        dvars = load_1d(dvars_path)
                    except Exception as e:
                        warnings.warn(f"Could not load DVARS file: {dvars_path}\nError: {e}")
                        continue

                    fd_al, dv_al, orig_idx, note, n_used = align_fd_dvars(fd, dvars)

                    row = {
                        "task": task, "subject": subject, "session": session,
                        "pipeline": pipeline,
                        "n_fd": len(fd), "n_dvars": len(dvars),
                        "alignment_note": note,
                        "n_censored_frames_removed": np.nan,
                        "n_used_censor_excluded": n_used,
                        "spearman_rho_censor_excluded": np.nan,
                        "spearman_pval_censor_excluded": np.nan,
                        "censor_correction_note": "",
                    }

                    if fd_al is not None:
                        design_path = build_design_path(task, subject, session, pipeline)
                        censored = get_spike_censored_indices(design_path)

                        if censored is None:
                            missing_design.append(str(design_path))
                            row["censor_correction_note"] = "design_mat_unavailable_NOT_COMPUTED"
                        else:
                            # exclude both the censored frame (t) and the
                            # following frame (t+1), since DVARS(t+1) is
                            # also built from the forced-zero residual(t)
                            exclude = set()
                            for t in censored:
                                exclude.add(t)
                                exclude.add(t + 1)

                            keep_mask = np.array([idx not in exclude for idx in orig_idx])
                            fd_c = fd_al[keep_mask]
                            dv_c = dv_al[keep_mask]
                            n_removed = int((~keep_mask).sum())

                            rho_c, pval_c, note_c = compute_rho(
                                fd_c, dv_c, note_suffix=f"n_spikes_detected={len(censored)}"
                            )

                            row["n_censored_frames_removed"] = n_removed
                            row["n_used_censor_excluded"] = len(fd_c)
                            row["spearman_rho_censor_excluded"] = rho_c
                            row["spearman_pval_censor_excluded"] = pval_c
                            row["censor_correction_note"] = note_c if note_c else "corrected"
                    else:
                        row["censor_correction_note"] = "alignment_failed"

                    rows.append(row)

    df = pd.DataFrame(rows)   # in-memory session-level intermediary — not written to disk

# ---------------------------------------------------------------
# Subject-level summary: PLAIN AVERAGE of the two session-level
# rho values e.g. not pooling RAW FD/DVARS values across sessions
    subj_rows = []
    if len(df):
        grp_cols = ["task", "subject", "pipeline"]
        for keys, g in df.groupby(grp_cols):
            rho_vals = g["spearman_rho_censor_excluded"].dropna()
            rho_subj = rho_vals.mean() if len(rho_vals) else np.nan
            subj_rows.append({
                "task": keys[0], "subject": keys[1], "pipeline": keys[2],
                "n_sessions_available": len(g),
                "n_sessions_used_in_average": len(rho_vals),
                "spearman_rho_subject_level": rho_subj,
            })

    subj_df = pd.DataFrame(subj_rows)
    subj_df.to_csv(OUT_SUBJECT, index=False)

    print()
    print("=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"Wrote {len(subj_df)} subject-level rows to: {OUT_SUBJECT.resolve()}")
    print(f"(session-level results used as an in-memory intermediary — {len(df)} rows, not saved to disk)")

    if len(df):
        n_mismatch = df["alignment_note"].fillna("").str.startswith("LENGTH_MISMATCH").sum()
        n_nan_rho_corr = df["spearman_rho_censor_excluded"].isna().sum()
        n_design_missing = df["censor_correction_note"].fillna("").str.startswith("design_mat_unavailable").sum()

        print("\nSUMMARY")
        print("-" * 70)
        print(f"Total run x pipeline results:                  {len(df)}")
        print(f"Unresolved FD/DVARS length mismatches:         {n_mismatch}")
        print(f"NaN censor-excluded rho:                       {n_nan_rho_corr}")
        print(f"Rows where design.mat missing (NOT COMPUTED):  {n_design_missing}")

        if n_design_missing:
            print("\n*** design.mat could not be found for some rows -- rho was NOT      ***")
            print("*** computed for these (left as NaN), not silently substituted.      ***")
            print("*** Verify DESIGN_MAT_ROOT / build_design_path() and rerun.          ***")

    if missing_fd:
        print(f"\n{len(missing_fd)} missing fd_raw.txt file(s). Examples:")
        for p in missing_fd[:5]:
            print(f"  {p}")

    if missing_dvars:
        print(f"\n{len(missing_dvars)} missing dvars_raw.txt file(s). Examples:")
        for p in missing_dvars[:5]:
            print(f"  {p}")

    if missing_design:
        print(f"\n{len(set(missing_design))} unique missing design.mat file(s). Examples:")
        for p in sorted(set(missing_design))[:5]:
            print(f"  {p}")

    print()
    print("=" * 70)


if __name__ == "__main__":
    main()

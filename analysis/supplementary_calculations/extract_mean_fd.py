#!/usr/bin/env python3
"""
Extract mean framewise displacement (FD) per subject/session/task
from fmriprep confounds TSV files.

Output: mean_fd.csv with columns:
    subject, session, task, n_volumes, n_fd_missing, mean_fd, median_fd, max_fd, pct_fd_gt_0.5

Usage:
    python extract_mean_fd.py
    python extract_mean_fd.py --fmriprep_dir /path/to/fmriprep --output mean_fd.csv
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd

FMRIPREP_DIR = "$basedir/fmriprep_output"
TASKS = ["nback", "gonogo", "mid", "oasis"]
SESSIONS = ["ses-1", "ses-2"]

def extract_fd(tsv_path):
    """Read a confounds TSV and return FD summary stats."""
    df = pd.read_csv(tsv_path, sep="\t")
    
    if "framewise_displacement" not in df.columns:
        return None
    
    fd = df["framewise_displacement"]
    n_total = len(fd)
    n_missing = fd.isna().sum()
    
    # First row is always NaN for FD — expected
    fd_clean = fd.dropna().values
    
    if len(fd_clean) == 0:
        return None
    
    return {
        "n_volumes": n_total,
        "n_fd_missing": int(n_missing),
        "mean_fd": float(np.mean(fd_clean)),
        "median_fd": float(np.median(fd_clean)),
        "max_fd": float(np.max(fd_clean)),
        "pct_fd_gt_0.2": float(np.mean(fd_clean > 0.2) * 100),
        "pct_fd_gt_0.5": float(np.mean(fd_clean > 0.5) * 100),
    }

def main(fmriprep_dir, output_path):
    records = []

    sub_dirs = sorted(glob.glob(os.path.join(fmriprep_dir, "sub-*")))
    subjects = [os.path.basename(d) for d in sub_dirs if os.path.isdir(d)]
    
    print(f"Found {len(subjects)} subjects")
    
    for sub in subjects:
        for ses in SESSIONS:
            func_dir = os.path.join(fmriprep_dir, sub, ses, "func")
            if not os.path.isdir(func_dir):
                continue
            
            for task in TASKS:
                pattern = os.path.join(
                    func_dir,
                    f"{sub}_{ses}_task-{task}_desc-confounds_timeseries.tsv"
                )
                matches = glob.glob(pattern)
                
                if not matches:
                    # Try with run entity in case some have run-01 etc.
                    pattern_run = os.path.join(
                        func_dir,
                        f"{sub}_{ses}_task-{task}_run-*_desc-confounds_timeseries.tsv"
                    )
                    matches = sorted(glob.glob(pattern_run))
                
                for tsv_path in matches:
                    stats = extract_fd(tsv_path)
                    if stats is None:
                        print(f"  WARNING: no FD column in {tsv_path}")
                        continue
                    
                    record = {
                        "subject": sub,
                        "session": ses,
                        "task": task,
                        "tsv_file": os.path.basename(tsv_path),
                        **stats,
                    }
                    records.append(record)
                    print(f"  {sub} {ses} {task}: mean_fd={stats['mean_fd']:.4f}, "
                          f"n_vol={stats['n_volumes']}, pct>0.5mm={stats['pct_fd_gt_0.5']:.1f}%")
    
    if not records:
        print("No records found — check FMRIPREP_DIR path.")
        return
    
    out_df = pd.DataFrame(records)
    
    # Clean up subject/session labels for easier merging later
    out_df["subject_id"] = out_df["subject"].str.replace("sub-", "", regex=False)
    out_df["session_id"] = out_df["session"].str.replace("ses-", "", regex=False)
    
    # Reorder columns
    cols = ["subject", "subject_id", "session", "session_id", "task",
            "mean_fd", "median_fd", "max_fd",
            "pct_fd_gt_0.2", "pct_fd_gt_0.5",
            "n_volumes", "n_fd_missing", "tsv_file"]
    out_df = out_df[cols]
    
    out_df.to_csv(output_path, index=False)
    print(f"\nSaved {len(out_df)} rows to: {output_path}")
    print("\nSummary:")
    print(out_df.groupby("task")["mean_fd"].describe().round(4))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fmriprep_dir", default=FMRIPREP_DIR)
    parser.add_argument("--output", default="mean_fd.csv")
    args = parser.parse_args()
    
    main(args.fmriprep_dir, args.output)

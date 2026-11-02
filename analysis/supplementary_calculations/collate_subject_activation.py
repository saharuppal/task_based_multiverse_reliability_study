#!/usr/bin/env python3
"""
Collate subject-level activation txt files into three CSVs:
    - mean_z.csv
    - n_suprathreshold.csv
    - mean_top1pct_z.csv

Each CSV has columns: subject, session, task, pipeline, cope, <measure>

Usage:
    python collate_subject_activation.py
"""

import os
import glob

BASE = "/path/to/project"
INPUT_DIR = os.path.join(BASE, "further_analysis", "subject_activation")
OUT_DIR = os.path.join(BASE, "further_analysis")

def parse_txt(path):
    record = {}
    with open(path) as f:
        for line in f:
            if ": " in line:
                key, val = line.strip().split(": ", 1)
                record[key.strip()] = val.strip()
    return record

def main():
    pattern = os.path.join(INPUT_DIR, "*", "*", "*", "*", "*_activation.txt")
    files = sorted(glob.glob(pattern))
    print(f"Found {len(files)} result files")

    mean_z_rows = []
    nvox_rows = []
    top1pct_rows = []

    for f in files:
        rec = parse_txt(f)
        base = {
            "subject":  rec.get("subject", "NA"),
            "session":  rec.get("session", "NA"),
            "task":     rec.get("task", "NA"),
            "pipeline": rec.get("pipeline", "NA"),
            "cope":     rec.get("cope", "NA"),
        }
        mean_z_rows.append({**base, "mean_z":          rec.get("mean_z", "NA")})
        nvox_rows.append(  {**base, "n_suprathreshold": rec.get("n_suprathreshold", "NA")})
        top1pct_rows.append({**base, "mean_top1pct_z":  rec.get("mean_top1pct_z", "NA")})

    import csv

    def write_csv(rows, filename, value_col):
        path = os.path.join(OUT_DIR, filename)
        fieldnames = ["subject", "session", "task", "pipeline", "cope", value_col]
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"Saved {len(rows)} rows to: {path}")

    write_csv(mean_z_rows,   "mean_z.csv",          "mean_z")
    write_csv(nvox_rows,     "n_suprathreshold.csv", "n_suprathreshold")
    write_csv(top1pct_rows,  "mean_top1pct_z.csv",   "mean_top1pct_z")

if __name__ == "__main__":
    main()
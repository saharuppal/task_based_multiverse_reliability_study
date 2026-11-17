############################
# FOR RUNNING IN IMPERIAL RDS
#####################

python3 - <<'PY'
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(".")

tasks = {
    "gonogo": 2,
    "nback": 2,
    "mid": 2,
    "oasis": 3,
}

# Use one pipeline purely to obtain the task EVs.
# We have already established that these are identical across pipelines.
pipeline = "24HMP_aCompCor_spike_reg"

def read_design_matrix(path):
    with open(path) as f:
        lines = f.readlines()

    start = next(
        i for i, line in enumerate(lines)
        if line.strip() == "/Matrix"
    )

    return np.array([
        [float(x) for x in line.split()]
        for line in lines[start + 1:]
        if line.strip()
    ])


results = []

# Loop through tasks, subjects, and sessions
for task, n_evs in tasks.items():

    firstlevel_dir = ROOT / "first_levels" / task
    fmriprep_dir = ROOT / "fmriprep_MW_new"

    for subject_dir in sorted(firstlevel_dir.glob("sub-*")):

        subject = subject_dir.name

        for session_dir in sorted(subject_dir.glob("ses-*")):

            session = session_dir.name

            feat_dir = (
                session_dir /
                f"{subject}_{session}_{task}_{pipeline}.feat"
            )

            design_file = feat_dir / "design.mat"

            confounds_file = (
                fmriprep_dir /
                subject /
                session /
                "func" /
                f"{subject}_{session}_task-{task}_desc-confounds_timeseries.tsv"
            )

            if not design_file.exists():
                continue

            # Check if confounds file exists for FD
            if not confounds_file.exists():
                continue

            design = read_design_matrix(design_file)

            confounds = pd.read_csv(confounds_file, sep="\t")

            # Extract FD as a NumPy array, non-numeric values to NaN
            fd = pd.to_numeric(
                confounds["framewise_displacement"],
                errors="coerce"
            ).to_numpy()

            # Checks for length mismatch between FD and design matrix
            if len(fd) != len(design):
                print(
                    f"LENGTH MISMATCH: "
                    f"{subject} {session} {task}: "
                    f"FD={len(fd)}, design={len(design)}"
                )
                continue

            valid = ~np.isnan(fd)

            # Calc number of valid FD values and total number of volumes
            row = {
                "subject": subject,
                "session": session,
                "task": task,
                "n_volumes": len(fd),
                "n_valid_fd": valid.sum(),
            }

            # Compute correlations between FD and each task EV.
            # Pearson's r
            for i in range(n_evs):

                ev = design[:, i]

                row[f"EV{i+1}_FD_r"] = np.corrcoef(
                    fd[valid],
                    ev[valid]
                )[0, 1]

            results.append(row)


out = pd.DataFrame(results)

output = ROOT / "further_analysis" / "fd_task_correlations.csv"

out.to_csv(output, index=False)

print()
print("Saved:", output)
print()
print(out.to_string(index=False))
print()
print("Number of runs:", len(out))
PY
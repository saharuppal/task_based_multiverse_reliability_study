############################
# FOR RUNNING IN IMPERIAL RDS
#####################

python3 - <<'PY'
import pandas as pd
import numpy as np

df = pd.read_csv("fd_task_correlations.csv")

# Convert EV columns into long format
long = df.melt(
    id_vars=["subject", "session", "task", "n_volumes", "n_valid_fd"],
    value_vars=["EV1_FD_r", "EV2_FD_r"],
    var_name="EV",
    value_name="r"
)

# Remove missing correlations
long = long.dropna(subset=["r"])

# Fisher z transform
long["fisher_z"] = np.arctanh(long["r"])

# Summarise by task and EV
summary = (
    long
    .groupby(["task", "EV"])
    .agg(
        n=("r", "count"),
        mean_r=("r", "mean"),
        median_r=("r", "median"),
        sd_r=("r", "std"),
        min_r=("r", "min"),
        max_r=("r", "max"),
        mean_fisher_z=("fisher_z", "mean")
    )
    .reset_index()
)

# Back-transform mean Fisher z
summary["fisher_mean_r"] = np.tanh(summary["mean_fisher_z"])

# Add proportions above useful descriptive magnitudes
for threshold in [0.2, 0.3]:
    prop = (
        long.assign(abs_r=long["r"].abs() > threshold)
        .groupby(["task", "EV"])["abs_r"]
        .mean()
        .reset_index(name="proportion_abs_r_above")
    )
    prop["threshold"] = threshold

    summary = summary.merge(
        prop,
        on=["task", "EV"],
        how="left"
    )

summary.to_csv("fd_task_correlations_summary.csv", index=False)

print(summary.to_string(index=False))
PY
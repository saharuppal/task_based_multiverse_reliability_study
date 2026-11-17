"""
Group-level analysis of FD-DVARS Spearman correlations across a denoising
pipeline multiverse (12 pipelines x 4 tasks x N subjects; sessions already
combined into one rho per subject/task/pipeline).

WHAT THIS SCRIPT DOES, STEP BY STEP
------------------------------------
Step 1  Sanity check          -- confirm the data is complete/well-formed before trusting anything downstream.
Step 2  Pipeline comparison   -- median/IQR of FD-DVARS rho per pipeline, overall and per task.
Step 3  Within-run ranking    -- rank the 12 pipelines *within each (task, subject)* (this is the core
                                  comparison, because it controls for subject/task differences by
                                  construction: every pipeline is applied to the exact same subject/task,
                                  so ranking within-run isolates the pipeline effect from everything else
                                  that varies between runs).
Step 4  Factor decomposition  -- (optional) descriptive paired differences for 24 vs 6 HMP, tissue
                                  regressors, and spike regression, stratified by task.
Step 5  Motion-bias check     -- (optional, OFF by default) does run-level motion severity relate to
                                  FD-DVARS rho, and does that relationship differ by pipeline? Not required
                                  to validate Steps 2/3 (within-run ranking already controls for run-level
                                  motion differences) -- only relevant if you want to ask whether a
                                  pipeline's benefit is motion-dependent. Requires fd_raw.txt files
                                  locally if enabled -- leave RUN_STEP5 = False unless you have them.

Inputs
------
fd_dvars_per_subject.csv  -- per-(task, subject, pipeline) Spearman rho between FD and DVARS,
                              sessions already combined via plain average. One row per
                              (task, subject, pipeline). Also carries n_fd_total, n_dvars_total,
                              n_used_censor_excluded_total (summed across sessions).
tdof_summary.tsv          -- REAL per-SESSION nuisance-regressor-column counts, read directly from
                              each run's design.mat. Since the rho file no longer has a session
                              column, this is aggregated to (task, subject, pipeline) here by
                              averaging each subject's sessions -- nothing about pipeline
                              composition is hard-coded/assumed.

Output
------
A single labeled Excel workbook (fd_dvars_group_results.xlsx) with one sheet per step,
human-readable column headers, and a README sheet explaining each tab.
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

try:
    from statsmodels.stats.multitest import multipletests
    HAVE_STATSMODELS = True
except ImportError:
    HAVE_STATSMODELS = False

# ============================================================
# CONFIG
# ============================================================
# EDIT THESE to your actual local paths before running.
RHO_FILE = Path("/Users/Liz/Documents/Work/Multiverse_pipelines/xlsx/fd_dvars_per_subject.csv")
TDOF_FILE = Path("/Users/Liz/Documents/Work/Multiverse_pipelines/xlsx/tdof_summary.tsv")
FD_RAW_BASE = Path("/Users/Liz/Documents/Work/Multiverse_pipelines")  # only used if RUN_STEP5=True

OUT_DIR = Path("/Users/Liz/Documents/Work/Multiverse_pipelines/fd_dvars_group_analysis_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_XLSX = OUT_DIR / "new_fd_dvars_results_spikes_censored_persubject_factor_decomp.xlsx"

RUN_STEP4 = True
RUN_STEP4_STATS = True  # requires statsmodels -- see docstring for what this does
STEP4_STRATIFY_BY_TASK = True
RUN_STEP5 = False   # see docstring -- optional, needs fd_raw.txt locally if enabled

PIPELINES = [
    "6HMP", "24HMP",
    "6HMP_WM_CSF", "24HMP_WM_CSF",
    "6HMP_aCompCor", "24HMP_aCompCor",
    "6HMP_spike_reg", "24HMP_spike_reg",
    "6HMP_WM_CSF_spike_reg", "24HMP_WM_CSF_spike_reg",
    "6HMP_aCompCor_spike_reg", "24HMP_aCompCor_spike_reg",
]
TASKS = ["nback", "gonogo", "mid", "oasis"]

# a place to collect (sheet_name -> (dataframe, description)) for the final workbook
SHEETS = {}


# ============================================================
# HELPERS
# ============================================================

def fisher_mean(x):
    """Mean of a correlation-like quantity, computed properly via Fisher z:
    transform each value to z, average in z-space, back-transform to r-space.
    A naive arithmetic mean of raw correlation coefficients is a biased
    estimator because r's sampling distribution is skewed away from 0."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan
    x = np.clip(x, -0.999999, 0.999999)
    return np.tanh(np.arctanh(x).mean())


def load_tdof_subject_level(path: Path) -> pd.DataFrame:
    """
    tdof_summary.tsv is per-SESSION (task, subject, session, pipeline,
    n_nuisance_cols, npts, pct_tdof_loss). The rho file is now per-subject
    (sessions already combined), so this aggregates tdof to the same
    (task, subject, pipeline) grain by averaging across that subject's
    sessions -- pct_tdof_loss is averaged directly (not recomputed from
    summed n_nuisance_cols/npts), since it's already the correct per-session
    ratio and averaging ratios directly is more correct than a ratio of sums
    here (npts can differ across sessions, e.g. a truncated acquisition).
    """
    if path.suffix.lower() in (".xlsx", ".xls"):
        t = pd.read_excel(path)
    else:
        t = pd.read_csv(path, sep="\t")
    t = t.rename(columns={"sub": "subject", "ses": "session"})
    required = {"task", "subject", "session", "pipeline", "n_nuisance_cols", "npts", "pct_tdof_loss"}
    missing = required - set(t.columns)
    if missing:
        raise ValueError(
            f"tdof_summary is missing expected columns: {missing}. "
            f"Found: {list(t.columns)}. Update the rename/required list to match your file."
        )
    agg = (
        t.groupby(["task", "subject", "pipeline"])
        .agg(
            n_nuisance_cols=("n_nuisance_cols", "mean"),
            npts=("npts", "mean"),
            pct_tdof_loss=("pct_tdof_loss", "mean"),
            n_sessions_in_tdof=("session", "nunique"),
        )
        .reset_index()
    )
    return agg


def load_rho_data(path: Path) -> pd.DataFrame:
    """
    Loads the subject-level (sessions already combined) rho file produced
    by fd_dvars_all.sh: fd_dvars_per_subject.csv. Renames
    spearman_rho_subject_level -> spearman_rho for consistency with the
    rest of this script.
    """
    if path.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)

    if "spearman_rho" not in df.columns:
        if "spearman_rho_subject_level" in df.columns:
            df = df.rename(columns={"spearman_rho_subject_level": "spearman_rho"})
        else:
            raise ValueError(
                f"Could not find 'spearman_rho' or 'spearman_rho_subject_level' in {path}. "
                f"Found columns: {list(df.columns)}"
            )
    required = {"task", "subject", "pipeline", "spearman_rho"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"rho file is missing expected columns: {missing}. Found: {list(df.columns)}")
    return df


def load_mean_fd_per_subject(base_dir, tasks, run_index):
    """
    Only used if RUN_STEP5=True. Averages mean FD across a subject's
    sessions (to match the now session-combined rho file) rather than
    returning one row per session.
    """
    rows = []
    for task in tasks:
        task_dir = base_dir / task
        if not task_dir.is_dir():
            continue
        for sub_dir in sorted(task_dir.glob("sub-*")):
            for ses_dir in sorted(sub_dir.glob("ses-*")):
                fd_path = ses_dir / "fd_raw.txt"
                if not fd_path.exists():
                    continue
                fd = np.genfromtxt(fd_path, dtype=float,
                                    missing_values=["n/a", "nan"], filling_values=np.nan)
                fd = np.atleast_1d(fd)
                rows.append({
                    "task": task, "subject": sub_dir.name, "session": ses_dir.name,
                    "mean_fd": np.nanmean(fd),
                })
    fd_df_session = pd.DataFrame(rows)
    fd_df = (
        fd_df_session.groupby(["task", "subject"])["mean_fd"]
        .mean()
        .reset_index()
    )
    check = run_index.merge(fd_df, on=["task", "subject"], how="left")
    n_missing = check["mean_fd"].isna().sum()
    if n_missing:
        warnings.warn(f"{n_missing} subject(s) have no matching mean FD -- check FD_RAW_BASE path.")
    return fd_df


def paired_descriptive(df_piv, pipe_a, pipe_b, label):
    diff = (df_piv[pipe_b] - df_piv[pipe_a]).dropna()
    return {
        "Comparison": label, "Baseline_pipeline": pipe_a, "Compared_pipeline": pipe_b,
        "N_runs": len(diff),
        "Median_difference": diff.median(),
        "IQR_of_difference": diff.quantile(0.75) - diff.quantile(0.25),
        "Mean_difference": diff.mean(),
        "Interpretation": (
            "Positive = Compared pipeline has HIGHER FD-DVARS rho than Baseline "
            "(i.e. Baseline shows better motion/signal separation)."
        ),
    }


def paired_wilcoxon(df_piv, pipe_a, pipe_b, label):
    diff = (df_piv[pipe_b] - df_piv[pipe_a]).dropna()
    if len(diff) < 2 or np.allclose(diff, 0):
        return {"Comparison": label, "Baseline_pipeline": pipe_a, "Compared_pipeline": pipe_b,
                "N_runs": len(diff), "Median_difference": np.nan, "W": np.nan, "z": np.nan,
                "p_value": np.nan, "effect_size_r": np.nan}
    res = stats.wilcoxon(diff, zero_method="pratt", mode="approx")
    z_stat = getattr(res, "zstatistic", np.nan)
    effect_r = abs(z_stat) / np.sqrt(len(diff)) if np.isfinite(z_stat) else np.nan
    return {"Comparison": label, "Baseline_pipeline": pipe_a, "Compared_pipeline": pipe_b,
            "N_runs": len(diff), "Median_difference": diff.median(),
            "W": res.statistic, "z": z_stat, "p_value": res.pvalue, "effect_size_r": effect_r}


# ============================================================
# LOAD + MERGE
# ============================================================

df = load_rho_data(RHO_FILE)
tdof_subj = load_tdof_subject_level(TDOF_FILE)
df = df.merge(
    tdof_subj[["task", "subject", "pipeline", "n_nuisance_cols", "npts", "pct_tdof_loss", "n_sessions_in_tdof"]],
    on=["task", "subject", "pipeline"], how="left",
)
n_unmatched = df["n_nuisance_cols"].isna().sum()
df["pct_nuisance"] = df["pct_tdof_loss"]  # already the correct per-subject-averaged ratio
df["run_id"] = df["task"] + "_" + df["subject"]  # no session component -- sessions already combined


# ============================================================
# STEP 1: SANITY CHECK
# ============================================================
runs_per_task = df[["task", "subject"]].drop_duplicates().groupby("task").size()
rho_all = df["spearman_rho"].dropna()

sanity = pd.DataFrame([
    {"Check": "Total rows", "Value": len(df)},
    {"Check": "Unique subjects", "Value": df["subject"].nunique()},
    {"Check": "Unique tasks", "Value": ", ".join(sorted(df["task"].unique()))},
    {"Check": "Pipelines found", "Value": len(df["pipeline"].unique())},
    {"Check": "Pipelines expected but missing", "Value": str(set(PIPELINES) - set(df["pipeline"].unique())) or "None"},
    {"Check": "Missing/NaN spearman_rho", "Value": df["spearman_rho"].isna().sum()},
    {"Check": "Rows with no tdof_summary match (regressor count missing)", "Value": n_unmatched},
    {"Check": "Rows built from only 1 session (not 2)", "Value": int((df.get("n_sessions_in_tdof", pd.Series(dtype=float)) == 1).sum())},
    {"Check": "Overall median rho", "Value": round(rho_all.median(), 4)},
    {"Check": "Overall IQR rho", "Value": f"{rho_all.quantile(0.25):.4f} - {rho_all.quantile(0.75):.4f}"},
    {"Check": "Overall Fisher-z mean rho", "Value": round(fisher_mean(rho_all), 4)},
    {"Check": "Skewness of rho distribution", "Value": round(stats.skew(rho_all), 4)},
])
runs_per_task_df = runs_per_task.reset_index().rename(columns={0: "N_subjects"})
runs_per_task_df.columns = ["Task", "N_subjects"]

SHEETS["01_Sanity_Check"] = (sanity, "Data completeness checks -- confirm these look right before trusting any other sheet.")
SHEETS["01_Subjects_per_Task"] = (runs_per_task_df, "Number of unique subjects found per task (sessions already combined).")


# ============================================================
# STEP 2: PIPELINE COMPARISON
# ============================================================
pipe_summary = (
    df.groupby("pipeline")["spearman_rho"]
    .agg(N_observations="count", Median_rho="median",
         Q25=lambda x: x.quantile(0.25), Q75=lambda x: x.quantile(0.75),
         Fisher_z_mean_rho=lambda x: fisher_mean(x))
    .reindex(PIPELINES).reset_index()
)
pipe_summary["IQR"] = pipe_summary["Q75"] - pipe_summary["Q25"]
pipe_summary = pipe_summary.rename(columns={"pipeline": "Pipeline"}).sort_values("Median_rho")

task_pipe_median = (
    df.pivot_table(index="pipeline", columns="task", values="spearman_rho", aggfunc="median")
    .reindex(index=PIPELINES).reset_index().rename(columns={"pipeline": "Pipeline"})
)
task_pipe_iqr = (
    df.pivot_table(index="pipeline", columns="task", values="spearman_rho",
                    aggfunc=lambda x: x.quantile(0.75) - x.quantile(0.25))
    .reindex(index=PIPELINES).reset_index().rename(columns={"pipeline": "Pipeline"})
)

SHEETS["02_Pipeline_Overall"] = (
    pipe_summary,
    "Median/IQR FD-DVARS rho per pipeline, pooled across all tasks/subjects (one value per subject, "
    "sessions already combined). Lower median rho = better motion/signal separation for that pipeline. "
    "Sorted best-first.",
)
SHEETS["02_Median_by_Task"] = (
    task_pipe_median,
    "Median FD-DVARS rho per pipeline, broken out by task -- shows whether the best pipeline differs by task.",
)
SHEETS["02_IQR_by_Task"] = (
    task_pipe_iqr,
    "Spread (IQR) of FD-DVARS rho per pipeline within each task -- a pipeline with a low median but huge IQR is inconsistent.",
)


# ============================================================
# STEP 3: WITHIN-RUN RANKING
# ============================================================
df["rank"] = df.groupby("run_id")["spearman_rho"].rank(method="average", ascending=True)

rank_summary = (
    df.groupby("pipeline")["rank"]
    .agg(Median_rank="median", Q25_rank=lambda x: x.quantile(0.25), Q75_rank=lambda x: x.quantile(0.75),
         Pct_runs_ranked_1st=lambda x: (x == 1).mean() * 100,
         Pct_runs_top3=lambda x: (x <= 3).mean() * 100,
         Pct_runs_bottom3=lambda x: (x >= 10).mean() * 100)
    .reindex(PIPELINES).reset_index().rename(columns={"pipeline": "Pipeline"})
)
rank_summary["IQR_rank"] = rank_summary["Q75_rank"] - rank_summary["Q25_rank"]
rank_summary = rank_summary.sort_values("Median_rank")

task_rank_median = (
    df.pivot_table(index="pipeline", columns="task", values="rank", aggfunc="median")
    .reindex(index=PIPELINES).reset_index().rename(columns={"pipeline": "Pipeline"})
)

SHEETS["03_WithinRun_Rank"] = (
    rank_summary,
    "Each (task, subject)'s 12 pipelines are ranked 1 (lowest/best rho) to 12 (highest/worst rho) against "
    "EACH OTHER, then those ranks are summarised per pipeline. This is the core comparison: because every "
    "pipeline is applied to the identical subject/task, ranking within-run cancels out subject/task "
    "differences and isolates the pipeline effect specifically.",
)
SHEETS["03_Rank_by_Task"] = (
    task_rank_median,
    "Median within-run rank per pipeline, broken out by task -- the direct answer to 'does the best "
    "pipeline depend on task'.",
)


# ============================================================
# STEP 4: FACTOR DECOMPOSITION
# ============================================================
if RUN_STEP4:
    pairs_24_vs_6 = [
        ("6HMP", "24HMP"), ("6HMP_WM_CSF", "24HMP_WM_CSF"),
        ("6HMP_aCompCor", "24HMP_aCompCor"), ("6HMP_spike_reg", "24HMP_spike_reg"),
        ("6HMP_WM_CSF_spike_reg", "24HMP_WM_CSF_spike_reg"),
        ("6HMP_aCompCor_spike_reg", "24HMP_aCompCor_spike_reg"),
    ]
    pairs_tissue = [
        ("6HMP", "6HMP_WM_CSF", "WM_CSF vs None (6HMP)"),
        ("24HMP", "24HMP_WM_CSF", "WM_CSF vs None (24HMP)"),
        ("6HMP", "6HMP_aCompCor", "aCompCor vs None (6HMP)"),
        ("24HMP", "24HMP_aCompCor", "aCompCor vs None (24HMP)"),
        ("6HMP_WM_CSF", "6HMP_aCompCor", "aCompCor vs WM_CSF (6HMP)"),
        ("24HMP_WM_CSF", "24HMP_aCompCor", "aCompCor vs WM_CSF (24HMP)"),
    ]
    pairs_spike = [
        ("6HMP", "6HMP_spike_reg", "Spike vs None (6HMP)"),
        ("24HMP", "24HMP_spike_reg", "Spike vs None (24HMP)"),
        ("6HMP_WM_CSF", "6HMP_WM_CSF_spike_reg", "Spike vs None (6HMP_WM_CSF)"),
        ("24HMP_WM_CSF", "24HMP_WM_CSF_spike_reg", "Spike vs None (24HMP_WM_CSF)"),
        ("6HMP_aCompCor", "6HMP_aCompCor_spike_reg", "Spike vs None (6HMP_aCompCor)"),
        ("24HMP_aCompCor", "24HMP_aCompCor_spike_reg", "Spike vs None (24HMP_aCompCor)"),
    ]

    def run_family(df_sub, pairs, three_col=False):
        results = []
        for p in pairs:
            pa, pb, label = (p[0], p[1], p[2]) if three_col else (p[0], p[1], f"{p[1]} vs {p[0]}")
            piv = df_sub.pivot(index="run_id", columns="pipeline", values="spearman_rho")
            if pa not in piv.columns or pb not in piv.columns:
                continue
            results.append(paired_wilcoxon(piv, pa, pb, label) if RUN_STEP4_STATS
                            else paired_descriptive(piv, pa, pb, label))
        out = pd.DataFrame(results)
        if RUN_STEP4_STATS and HAVE_STATSMODELS and len(out) and out["p_value"].notna().any():
            valid = out["p_value"].notna()
            _, p_adj, _, _ = multipletests(out.loc[valid, "p_value"], method="fdr_bh")
            out.loc[valid, "p_value_FDR_corrected"] = p_adj
        return out

    task_scopes = TASKS if STEP4_STRATIFY_BY_TASK else [None]
    all_step4 = []
    for task in task_scopes:
        df_sub = df if task is None else df[df["task"] == task]
        scope_label = "ALL_TASKS_POOLED" if task is None else task
        for pairs, fam, three_col in [(pairs_24_vs_6, "Motion_model_24_vs_6HMP", False),
                                       (pairs_tissue, "Tissue_regressors", True),
                                       (pairs_spike, "Spike_regression", True)]:
            r = run_family(df_sub, pairs, three_col=three_col)
            if len(r):
                r["Factor_family"] = fam
                r["Task"] = scope_label
                all_step4.append(r)

    step4_out = pd.concat(all_step4, ignore_index=True) if all_step4 else pd.DataFrame()
    SHEETS["04_Factor_Decomposition"] = (
        step4_out,
        "Descriptive paired differences (Compared - Baseline pipeline) isolating one factor at a time "
        "(24 vs 6 motion params; adding WM/CSF or aCompCor; adding spike regression), stratified by task. "
        "Median_difference > 0 means the 'Compared_pipeline' has HIGHER rho (worse motion/signal "
        "separation) than 'Baseline_pipeline'.",
    )


# ============================================================
# STEP 5 (OPTIONAL, OFF BY DEFAULT)
# ============================================================
if RUN_STEP5:
    run_index = df[["task", "subject"]].drop_duplicates()
    fd_df = load_mean_fd_per_subject(FD_RAW_BASE, TASKS, run_index)
    df = df.merge(fd_df, on=["task", "subject"], how="left")

    motion_bias = []
    for pipe in PIPELINES:
        sub = df[df["pipeline"] == pipe].dropna(subset=["mean_fd", "spearman_rho"])
        if len(sub) < 3:
            motion_bias.append({"Pipeline": pipe, "N": len(sub), "Spearman_r": np.nan, "p_value": np.nan})
            continue
        r_sp, p_sp = stats.spearmanr(sub["mean_fd"], sub["spearman_rho"])
        motion_bias.append({"Pipeline": pipe, "N": len(sub), "Spearman_r": r_sp, "p_value": p_sp})
    df_mb = pd.DataFrame(motion_bias)
    SHEETS["05_Motion_Bias_OPTIONAL"] = (
        df_mb,
        "OPTIONAL/exploratory: correlation between a subject's real mean FD (averaged across sessions) "
        "and that subject's FD-DVARS rho, per pipeline. NOT required to validate Steps 2/3. Only relevant "
        "if asking whether a pipeline's apparent benefit depends on how much a subject moved.",
    )


# ============================================================
# WRITE LABELED WORKBOOK
# ============================================================
readme_rows = [{"Sheet": name, "Description": desc} for name, (_, desc) in SHEETS.items()]
readme_df = pd.DataFrame(readme_rows)

with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
    readme_df.to_excel(writer, sheet_name="00_README", index=False)
    for name, (data, _) in SHEETS.items():
        data.to_excel(writer, sheet_name=name[:31], index=False)  # Excel sheet name limit = 31 chars

print(f"Done. Workbook written to: {OUT_XLSX.resolve()}")
print(f"Sheets: {list(SHEETS.keys())}")
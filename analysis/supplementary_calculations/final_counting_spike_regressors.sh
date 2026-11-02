#!/bin/bash

# ── Config ────────────────────────────────────────────────────────────────
base_dir="/path/to/project"
out_file="spike_regressor_counts.csv"

tasks=("gonogo" "oasis" "nback" "mid")
sessions=(1 2)

# ── Diagnostic: confirm base_dir exists and is readable ───────────────────
echo "=== DIAGNOSTIC: checking base_dir ==="
echo "base_dir = '${base_dir}'"
if [[ -d "$base_dir" ]]; then
  echo "base_dir exists: YES"
else
  echo "base_dir exists: NO -- STOPPING, fix base_dir path"
  exit 1
fi
#!/bin/bash

# ── Config ────────────────────────────────────────────────────────────────
base_dir="/path/to/project"
out_file="spike_regressor_counts.csv"

tasks=("gonogo" "oasis" "nback" "mid")
sessions=(1 2)

# ── Write header ──────────────────────────────────────────────────────────
echo "subject,session,task,n_motion_outliers,file_found" > "$out_file"

n_found=0
n_missing=0

# ── Loop over actual subject directories (skips .html files automatically) ──
for sub_dir in "${base_dir}"/sub-*/; do
  sub=$(basename "$sub_dir")   # e.g. sub-001

  for ses in "${sessions[@]}"; do
    ses_label="ses-${ses}"

    for task in "${tasks[@]}"; do
      tsv_file="${base_dir}/${sub}/${ses_label}/func/${sub}_${ses_label}_task-${task}_desc-confounds_timeseries.tsv"

      if [[ -f "$tsv_file" ]]; then
        n_outliers=$(head -n 1 "$tsv_file" | tr '\t' '\n' | grep -c "motion_outlier")
        found="yes"
        n_found=$((n_found + 1))
      else
        n_outliers="NA"
        found="no"
        n_missing=$((n_missing + 1))
      fi

      echo "${sub},${ses_label},${task},${n_outliers},${found}" >> "$out_file"
    done
  done
done

echo
echo "=== SUMMARY ==="
echo "Found:   ${n_found}"
echo "Missing: ${n_missing}"
echo "Results written to ${out_file}"
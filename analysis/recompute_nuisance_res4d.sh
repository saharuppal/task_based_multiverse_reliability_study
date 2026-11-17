#!/bin/bash
# recompute_nuisance_res4d.sh
# Usage: recompute_nuisance_res4d.sh <SUB> <SES>
#   SUB e.g. 017   SES e.g. 1
set -uo pipefail   # NOTE: no -e — a failure in one pipeline must not kill the whole job
module load FSL/6.0.7.17 2>/dev/null || true
source "$FSLDIR/etc/fslconf/fsl.sh" 2>/dev/null || true

SUB=$1
SES=$2

PROJECT_ROOT=/rds/general/project/task_reliability_study/live/MultiverseReliability
FMRIPREP_DIR=$PROJECT_ROOT/fmriprep_MW_new
FIRSTLEVEL_DIR=$PROJECT_ROOT/first_levels
OUT_DIR=$PROJECT_ROOT/further_analysis/fd_dvars_dof
TDOF_SUMMARY=$OUT_DIR/tdof_summary.tsv
mkdir -p "$OUT_DIR"
[ -f "$TDOF_SUMMARY" ] || echo -e "sub\tses\ttask\tpipeline\tn_nuisance_cols\tnpts\tpct_tdof_loss" > "$TDOF_SUMMARY"

# ASSUMPTION: BIDS task labels below — verify against your actual filenames
TASKS=(gonogo nback mid oasis)
declare -A NTASKEV=( [gonogo]=2 [nback]=2 [mid]=2 [oasis]=3 )

PIPELINES=(6HMP 24HMP 6HMP_WM_CSF 24HMP_WM_CSF 6HMP_aCompCor 24HMP_aCompCor \
6HMP_spike_reg 24HMP_spike_reg 6HMP_WM_CSF_spike_reg 24HMP_WM_CSF_spike_reg \
6HMP_aCompCor_spike_reg 24HMP_aCompCor_spike_reg)

# All per-pipeline work lives in this function. Its internal commands run under
# their own set -e (function body), so ONE failing command aborts only THIS
# pipeline (function returns non-zero) — the caller catches that and continues.
process_pipeline() {
  local TASK=$1 PIPE=$2 NTASK=$3 TASK_OUT=$4
  set -e

  local FEATDIR=$FIRSTLEVEL_DIR/$TASK/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_${TASK}_${PIPE}.feat
  local DESIGN=$FEATDIR/design.mat
  local MASK=$FEATDIR/mask.nii.gz
  local FUNC=$FEATDIR/filtered_func_data.nii.gz
  local OUTSUB=$TASK_OUT/$PIPE

  if [ ! -d "$FEATDIR" ]; then
    echo "WARN: missing FEAT dir $FEATDIR" >&2; return 1
  fi
  for f in "$DESIGN" "$MASK" "$FUNC"; do
    if [ ! -s "$f" ]; then
      echo "WARN: missing/empty required file $f (FEAT run likely incomplete)" >&2
      return 1
    fi
  done

  mkdir -p "$OUTSUB"

  local NEED_GLM=1 NEED_DVARS=1 NEED_TDOF=1
  if [ -s "$OUTSUB/res4d_nuisance.nii.gz" ]; then
    NEED_GLM=0
    NEED_TDOF=0   # already logged to tdof_summary.tsv when first computed
  fi
  if [ -s "$OUTSUB/dvars_raw.txt" ] && [ "$(wc -l < "$OUTSUB/dvars_raw.txt")" -gt 1 ]; then
    NEED_DVARS=0  # >1 line rules out the old fsl_motion_outliers failure mode (single "0" line)
  fi

  if [ "$NEED_GLM" -eq 0 ] && [ "$NEED_DVARS" -eq 0 ]; then
    echo "$SUB $SES $TASK $PIPE already done. Skipping."
    return 0
  fi

  local NPTS TM1
  NPTS=$(fslval "$FUNC" dim4 | tr -d ' ')
  TM1=$((NPTS - 1))

  if [ "$NEED_GLM" -eq 1 ]; then
    # ---- Strip task EV columns from the ALREADY-FILTERED design matrix ----
    Vest2Text "$DESIGN" "$OUTSUB/design_full.txt"
    awk -v n="$NTASK" '{s=""; for(i=n+1;i<=NF;i++) s=s $i" "; print s}' \
      "$OUTSUB/design_full.txt" > "$OUTSUB/design_nuisance.txt"
    Text2Vest "$OUTSUB/design_nuisance.txt" "$OUTSUB/design_nuisance.mat"

    # ---- Nuisance-only GLM residuals ----
    fsl_glm -i "$FUNC" -d "$OUTSUB/design_nuisance.mat" -m "$MASK" \
      -o "$OUTSUB/betas_nuisance" --out_res="$OUTSUB/res4d_nuisance.nii.gz" --demean
  fi

  if [ "$NEED_DVARS" -eq 1 ]; then
    # ---- Raw DVARS on nuisance-only residuals ----
    # Pure FSL: DVARS_t = sqrt(mean_over_mask((Y_t - Y_(t-1))^2)), computed on the
    # whole run at once via a lagged-pair subtraction + fslmeants, not a per-timepoint loop.
    fslroi "$OUTSUB/res4d_nuisance.nii.gz" "$OUTSUB/_lag_t1" 0 "$TM1"   # volumes 1..T-1
    fslroi "$OUTSUB/res4d_nuisance.nii.gz" "$OUTSUB/_lag_t2" 1 "$TM1"   # volumes 2..T
    fslmaths "$OUTSUB/_lag_t2" -sub "$OUTSUB/_lag_t1" -sqr "$OUTSUB/_diffsq"
    fslmeants -i "$OUTSUB/_diffsq" -m "$MASK" -o "$OUTSUB/_meansq.txt"
    { echo "nan"; awk '{print sqrt($1)}' "$OUTSUB/_meansq.txt"; } > "$OUTSUB/dvars_raw.txt"
    rm -f "$OUTSUB/_lag_t1.nii.gz" "$OUTSUB/_lag_t2.nii.gz" "$OUTSUB/_diffsq.nii.gz" "$OUTSUB/_meansq.txt"
  fi

  if [ "$NEED_TDOF" -eq 1 ]; then
    # ---- tDOF loss ----
    local NCOLS_NUISANCE PCT_LOSS
    NCOLS_NUISANCE=$(awk 'NR==1{print NF}' "$OUTSUB/design_nuisance.txt")
    PCT_LOSS=$(echo "scale=3; 100*$NCOLS_NUISANCE/$NPTS" | bc)
    echo -e "sub-${SUB}\tses-${SES}\t${TASK}\t${PIPE}\t${NCOLS_NUISANCE}\t${NPTS}\t${PCT_LOSS}" >> "$TDOF_SUMMARY"
  fi
}

for TASK in "${TASKS[@]}"; do
  TASK_OUT=$OUT_DIR/$TASK/sub-${SUB}/ses-${SES}
  mkdir -p "$TASK_OUT"

  # ---- FD (pipeline-independent, extracted once per task) ----
  CONF=$FMRIPREP_DIR/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}_desc-confounds_timeseries.tsv
  if [ -f "$CONF" ]; then
    FD_COL=$(head -1 "$CONF" | tr '\t' '\n' | grep -n '^framewise_displacement$' | cut -d: -f1)
    awk -F'\t' -v c="$FD_COL" 'NR>1{print $c}' "$CONF" > "$TASK_OUT/fd_raw.txt"
  else
    echo "WARN: missing confounds file $CONF" >&2
  fi

  NTASK=${NTASKEV[$TASK]}

  for PIPE in "${PIPELINES[@]}"; do
    if ! process_pipeline "$TASK" "$PIPE" "$NTASK" "$TASK_OUT"; then
      echo "SKIPPED: sub-${SUB} ses-${SES} ${TASK} ${PIPE} (see WARN above)" >&2
    fi
  done
done

echo "Done: sub-${SUB} ses-${SES}"

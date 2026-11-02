#!/bin/bash
# top1pct_zstat_template.sh
# Calculates the 99th percentile Z threshold and mean of top 1% Z-values
# from the unthresholded zstat image for a given task/pipeline/cope.
# Expects TASK, PIPELINE, COPE to be exported by the calling script.

set -euo pipefail

# ── Input unthresholded zstat image ───────────────────────────────────────────
ZSTAT="/path/to/project/group_levels/${TASK}/${PIPELINE}.gfeat/${COPE}.feat/stats/zstat1.nii.gz"

if [ ! -f "$ZSTAT" ]; then
    echo "ERROR: zstat1 not found: ${ZSTAT}"
    exit 1
fi

# ── Output directory ───────────────────────────────────────────────────────────
OUT_DIR="/path/to/project/group_levels/top1pct_zstat/${TASK}/${PIPELINE}/${COPE}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# ── Step 1: Get 99th percentile Z value (lower bound of top 1%) ───────────────
P99=$(fslstats "$ZSTAT" -P 99 | awk '{print $1}')
echo "99th percentile Z : ${P99}"

# ── Step 2: Threshold image to keep only voxels above 99th percentile ─────────
fslmaths "$ZSTAT" -thr "$P99" top1pct_mask.nii.gz

# ── Step 3: Count voxels in top 1% and get their mean ─────────────────────────
NVOX=$(fslstats top1pct_mask.nii.gz -V | awk '{print $1}')
MEAN_TOP1=$(fslstats top1pct_mask.nii.gz -M)

echo "N voxels in top 1%  : ${NVOX}"
echo "Mean Z in top 1%    : ${MEAN_TOP1}"

# ── Step 4: Write results ──────────────────────────────────────────────────────
cat > top1pct_result.txt << EOF
P99_Z: ${P99}
Mean_top1pct_Z: ${MEAN_TOP1}
N_voxels_top1pct: ${NVOX}
EOF

# Single-line summary for easy parsing — mirrors other output styles
echo "Mean_top1pct_Z: ${MEAN_TOP1}" > ${TASK}_${COPE}_${PIPELINE}_top1pct_zstat.txt

echo "Done: ${TASK} / ${PIPELINE} / ${COPE} — P99 = ${P99}, Mean top 1% = ${MEAN_TOP1}"

# ── Cleanup ────────────────────────────────────────────────────────────────────
rm -f top1pct_mask.nii.gz

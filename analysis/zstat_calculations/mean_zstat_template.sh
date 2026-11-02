#!/bin/bash
# mean_zstat_template.sh
# Calculates mean Z-stat of all non-zero voxels in thresh_zstat1
# for a given task/pipeline/cope.
# Expects TASK, PIPELINE, COPE to be exported by the calling script.

set -euo pipefail

# ── Input thresholded zstat image ─────────────────────────────────────────────
THRESH_ZSTAT="/path/to/project/group_levels/${TASK}/${PIPELINE}.gfeat/${COPE}.feat/thresh_zstat1.nii.gz"

if [ ! -f "$THRESH_ZSTAT" ]; then
    echo "ERROR: thresh_zstat1 not found: ${THRESH_ZSTAT}"
    exit 1
fi

# ── Output directory ───────────────────────────────────────────────────────────
OUT_DIR="/path/to/project/group_levels/mean_zstat/${TASK}/${PIPELINE}/${COPE}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# ── Calculate mean of non-zero voxels ─────────────────────────────────────────
# fslstats -M returns mean of all non-zero voxels
# fslstats -V returns number of non-zero voxels and volume
MEAN_Z=$(fslstats "$THRESH_ZSTAT" -M)
NVOX=$(fslstats "$THRESH_ZSTAT" -V | awk '{print $1}')

# Guard against empty image (no suprathreshold voxels)
if [ "$NVOX" -eq 0 ]; then
    echo "WARNING: No suprathreshold voxels for ${TASK}/${PIPELINE}/${COPE} — writing NA"
    echo "Mean_Z: NA" > ${TASK}_${COPE}_${PIPELINE}_mean_zstat.txt
    echo "Mean_Z: NA"
    exit 0
fi

echo "Mean Z-stat : ${MEAN_Z}"
echo "N voxels    : ${NVOX}"

# Full results file
cat > mean_zstat_result.txt << EOF
Mean_Z: ${MEAN_Z}
N_voxels: ${NVOX}
EOF

# Single-line summary — consistent with Dice/Spearman/ICC output style
echo "Mean_Z: ${MEAN_Z}" > ${TASK}_${COPE}_${PIPELINE}_mean_zstat.txt

echo "Done: ${TASK} / ${PIPELINE} / ${COPE} — Mean Z = ${MEAN_Z}"

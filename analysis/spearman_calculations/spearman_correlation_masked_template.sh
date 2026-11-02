#!/bin/bash
# spearman_correlation_masked_template.sh
# Computes voxel-wise Spearman correlation between ses-1 and ses-2
# thresholded zstat images within a task-specific ROI mask.
# Expects TASK, PIPELINE, COPE to be exported by the calling script.

set -euo pipefail

# ── Task-specific masks ────────────────────────────────────────────────────────
declare -A TASK_MASKS
TASK_MASKS["nback"]="/path/to/project/group_levels/templates/working_memory.nii.gz"
TASK_MASKS["oasis"]="/path/to/project/task_reliability_study/live/MultiverseReliability/group_levels/templates/emotional.nii.gz"
TASK_MASKS["gonogo"]="/path/to/project/task_reliability_study/live/MultiverseReliability/group_levels/templates/inhibition.nii.gz"
TASK_MASKS["mid"]="/path/to/project/task_reliability_study/live/MultiverseReliability/group_levels/templates/monetary_incentive.nii.gz"

MASK="${TASK_MASKS[$TASK]}"

if [ ! -f "$MASK" ]; then
    echo "ERROR: Mask not found for task ${TASK}: ${MASK}"
    exit 1
fi

# ── Input zstat images (unthresholded — use zstat1 for continuous values) ──────
# Using thresh_zstat1 keeps it consistent with the Dice calculation.

SES1_ZSTAT="/path/to/project/group_levels/means/${TASK}/ses-1/${PIPELINE}.gfeat/${COPE}.feat/stats/zstat1.nii.gz"
SES2_ZSTAT="/path/to/project/task_reliability_study/live/MultiverseReliability/group_levels/means/${TASK}/ses-2/${PIPELINE}.gfeat/${COPE}.feat/stats/zstat1.nii.gz"

if [ ! -f "$SES1_ZSTAT" ]; then
    echo "ERROR: ses-1 zstat not found: ${SES1_ZSTAT}"
    exit 1
fi
if [ ! -f "$SES2_ZSTAT" ]; then
    echo "ERROR: ses-2 zstat not found: ${SES2_ZSTAT}"
    exit 1
fi

# ── Output directory ───────────────────────────────────────────────────────────
OUT_DIR="/path/to/project/group_levels/spearman_masked/${TASK}/${PIPELINE}/${COPE}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# ── Step 1: Extract masked voxel values for each session ──────────────────────
# fslmeants with -m dumps one value per masked voxel to stdout
fslmeants \
    -i "$SES1_ZSTAT" \
    -m "$MASK" \
    --showall \
    -o ses1_voxels_raw.txt

fslmeants \
    -i "$SES2_ZSTAT" \
    -m "$MASK" \
    --showall \
    -o ses2_voxels_raw.txt

# fslmeants --showall outputs: x y z val per row when used on a 3D image.
# Extract just the value column (column 4):
awk '{print $4}' ses1_voxels_raw.txt > ses1_voxels.txt
awk '{print $4}' ses2_voxels_raw.txt > ses2_voxels.txt

# ── Step 2: Compute Spearman correlation in Python ────────────────────────────
python3 - <<'PYEOF'
import sys
import numpy as np
from scipy.stats import spearmanr

ses1 = np.loadtxt("ses1_voxels.txt")
ses2 = np.loadtxt("ses2_voxels.txt")

if len(ses1) == 0 or len(ses2) == 0:
    print("WARNING: No voxels found within mask — writing NA")
    with open("spearman_result.txt", "w") as f:
        f.write("Spearman_r: NA\nSpearman_p: NA\nN_voxels: 0\n")
    sys.exit(0)

if len(ses1) != len(ses2):
    print(f"ERROR: Voxel count mismatch: ses1={len(ses1)}, ses2={len(ses2)}")
    sys.exit(1)

rho, pval = spearmanr(ses1, ses2)
n = len(ses1)

print(f"Spearman r : {rho:.6f}")
print(f"p-value    : {pval:.6e}")
print(f"N voxels   : {n}")

with open("spearman_result.txt", "w") as f:
    f.write(f"Spearman_r: {rho:.6f}\n")
    f.write(f"Spearman_p: {pval:.6e}\n")
    f.write(f"N_voxels: {n}\n")
PYEOF

# ── Step 3: Write a clean single-value file (mirrors dice output style) ────────
RHO=$(grep "Spearman_r" spearman_result.txt | awk '{print $2}')
echo "Spearman_r: $RHO" > ${TASK}_${COPE}_${PIPELINE}_spearman.txt

echo "Done: ${TASK} / ${PIPELINE} / ${COPE} — r = ${RHO}"
#!/bin/bash
# icc_masked_extract_template.sh
# Thresholds neurosynth mask at 0, then extracts ICC(3,1) values
# within that mask for a given task/pipeline/cope.
# Expects TASK, PIPELINE, COPE to be exported by the calling script.

set -euo pipefail

# ── Task-specific masks ────────────────────────────────────────────────────────
declare -A TASK_MASKS
TASK_MASKS["nback"]="/rds/general/project/task_reliability_study/live/MultiverseReliability/group_levels/templates/working_memory.nii.gz"
TASK_MASKS["oasis"]="/rds/general/project/task_reliability_study/live/MultiverseReliability/group_levels/templates/emotional.nii.gz"
TASK_MASKS["gonogo"]="/rds/general/project/task_reliability_study/live/MultiverseReliability/group_levels/templates/inhibition.nii.gz"
TASK_MASKS["mid"]="/rds/general/project/task_reliability_study/live/MultiverseReliability/group_levels/templates/monetary_incentive.nii.gz"

MASK="${TASK_MASKS[$TASK]}"

if [ ! -f "$MASK" ]; then
    echo "ERROR: Mask not found for task ${TASK}: ${MASK}"
    exit 1
fi

# ── Input ICC image ────────────────────────────────────────────────────────────
ICC_DIR="/rds/general/project/task_reliability_study/live/MultiverseReliability/group_levels/export_filtered_func_data/analysis_scripts/icc_outputs/${TASK}/${PIPELINE}/${COPE}"

if [ ! -d "$ICC_DIR" ]; then
    echo "ERROR: ICC directory not found: ${ICC_DIR}"
    exit 1
fi

# Find the ICC NIfTI — try common names first, then fall back to any .nii.gz
ICC_IMG=""
for candidate in icc_map.nii.gz icc.nii.gz ICC_map.nii.gz icc_output.nii.gz; do
    if [ -f "${ICC_DIR}/${candidate}" ]; then
        ICC_IMG="${ICC_DIR}/${candidate}"
        break
    fi
done

# Replace the candidate search block with a direct path
ICC_IMG="${ICC_DIR}/icc3_est.nii.gz"

if [ ! -f "$ICC_IMG" ]; then
    echo "ERROR: ICC image not found: ${ICC_IMG}"
    exit 1
fi

echo "Using ICC image: ${ICC_IMG}"

# ── Output directory ───────────────────────────────────────────────────────────
OUT_DIR="/rds/general/project/task_reliability_study/live/MultiverseReliability/group_levels/icc_masked/${TASK}/${PIPELINE}/${COPE}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# ── Step 1: Threshold mask at 0 (binarise — keep voxels > 0) ──────────────────
fslmaths "$MASK" -thr 0.0001 -bin mask_bin.nii.gz

# ── Step 2: Apply binarised mask to ICC image ──────────────────────────────────
fslmaths "$ICC_IMG" -mas mask_bin.nii.gz icc_masked.nii.gz

# ── Step 3: Extract voxel values and compute summary stats in Python ───────────
python3 << PYEOF
import sys
import numpy as np
import nibabel as nib

icc_img  = nib.load("icc_masked.nii.gz")
mask_img = nib.load("mask_bin.nii.gz")

icc_data  = icc_img.get_fdata().flatten()
mask_data = mask_img.get_fdata().flatten()

# Keep only voxels inside the mask
mask_bool = mask_data > 0
icc_vals  = icc_data[mask_bool]

n = int(mask_bool.sum())
print(f"N voxels in mask: {n}")

if n == 0:
    print("WARNING: No voxels found within mask — writing NA")
    with open("icc_masked_summary.txt", "w") as f:
        f.write("Mean_ICC: NA\nMedian_ICC: NA\nSD_ICC: NA\nMin_ICC: NA\nMax_ICC: NA\nN_voxels: 0\n")
    with open("${TASK}_${COPE}_${PIPELINE}_icc_masked.txt", "w") as f:
        f.write("Mean_ICC: NA\n")
    sys.exit(0)

mean_icc   = float(np.mean(icc_vals))
median_icc = float(np.median(icc_vals))
sd_icc     = float(np.std(icc_vals))
min_icc    = float(np.min(icc_vals))
max_icc    = float(np.max(icc_vals))

print(f"Mean ICC   : {mean_icc:.6f}")
print(f"Median ICC : {median_icc:.6f}")
print(f"SD ICC     : {sd_icc:.6f}")
print(f"Min ICC    : {min_icc:.6f}")
print(f"Max ICC    : {max_icc:.6f}")
print(f"N voxels   : {n}")

# Full summary file
with open("icc_masked_summary.txt", "w") as f:
    f.write(f"Mean_ICC: {mean_icc:.6f}\n")
    f.write(f"Median_ICC: {median_icc:.6f}\n")
    f.write(f"SD_ICC: {sd_icc:.6f}\n")
    f.write(f"Min_ICC: {min_icc:.6f}\n")
    f.write(f"Max_ICC: {max_icc:.6f}\n")
    f.write(f"N_voxels: {n}\n")

# Single-line summary — mirrors Dice/Spearman output style
with open("${TASK}_${COPE}_${PIPELINE}_icc_masked.txt", "w") as f:
    f.write(f"Mean_ICC: {mean_icc:.6f}\n")

# Save all voxel values for distribution plots later
np.savetxt("icc_voxel_values.txt", icc_vals, fmt="%.6f")
PYEOF

echo "Done: ${TASK} / ${PIPELINE} / ${COPE} — Mean ICC = $(grep Mean_ICC icc_masked_summary.txt | awk '{print $2}')"

# ── Step 4: Clean up intermediate files ───────────────────────────────────────
# rm -f mask_bin.nii.gz icc_masked.nii.gz

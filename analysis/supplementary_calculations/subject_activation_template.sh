#!/bin/bash
#PBS -N subject_activation
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=00:30:00
#PBS -J 0-2399
#PBS -o /path/to/project/group_levels/subject_activation_logs/subject_activation.o
#PBS -e /path/to/project/group_levels/subject_activation_logs/subject_activation.e

set -euo pipefail

mkdir -p /path/to/project/group_levels/subject_activation_logs

module purge
module load tools/prod
module load FSL/6.0.7.17
source $FSLDIR/etc/fslconf/fsl.sh

BASE="/path/to/project"

PIPELINES=(
    "6HMP"
    "24HMP"
    "6HMP_aCompCor"
    "6HMP_spike_reg"
    "24HMP_aCompCor"
    "6HMP_WM_CSF"
    "24HMP_WM_CSF"
    "24HMP_spike_reg"
    "6HMP_aCompCor_spike_reg"
    "24HMP_aCompCor_spike_reg"
    "6HMP_WM_CSF_spike_reg"
    "24HMP_WM_CSF_spike_reg"
)

TASKS=(nback mid gonogo oasis)
SESSIONS=(ses-1 ses-2)
SUBJECTS=($(ls ${BASE}/first_levels/nback/ | grep "^sub-"))

# Build flat COMBOS array: 4 tasks x 12 pipelines x 25 subjects x 2 sessions = 2400
COMBOS=()
for TASK in "${TASKS[@]}"; do
    for PIPELINE in "${PIPELINES[@]}"; do
        for SUB in "${SUBJECTS[@]}"; do
            for SES in "${SESSIONS[@]}"; do
                COMBOS+=("${TASK}:${PIPELINE}:${SUB}:${SES}")
            done
        done
    done
done

COMBO="${COMBOS[$PBS_ARRAY_INDEX]}"
TASK=$(echo "$COMBO"     | cut -d: -f1)
PIPELINE=$(echo "$COMBO" | cut -d: -f2)
SUB=$(echo "$COMBO"      | cut -d: -f3)
SES=$(echo "$COMBO"      | cut -d: -f4)

echo "================================================"
echo "Job index : $PBS_ARRAY_INDEX"
echo "Task      : $TASK"
echo "Pipeline  : $PIPELINE"
echo "Subject   : $SUB"
echo "Session   : $SES"
echo "================================================"

# ── Task-specific zstat numbers ────────────────────────────────────────────────
case "$TASK" in
    nback|mid|gonogo)
        ZSTAT_NUMS=(3)
        COPE_LABELS=("cope3")
        ;;
    oasis)
        ZSTAT_NUMS=(4 5)
        COPE_LABELS=("cope4" "cope5")
        ;;
    *)
        echo "ERROR: Unknown task: ${TASK}"; exit 1 ;;
esac

# ── Feat directory ─────────────────────────────────────────────────────────────
FEAT_DIR="${BASE}/first_levels/${TASK}/${SUB}/${SES}/${SUB}_${SES}_${TASK}_${PIPELINE}.feat"

if [ ! -d "$FEAT_DIR" ]; then
    echo "WARNING: feat dir not found, skipping: ${FEAT_DIR}"
    exit 0
fi

OUT_DIR="${BASE}/group_levels/subject_activation/${TASK}/${PIPELINE}/${SUB}/${SES}"
mkdir -p "$OUT_DIR"

# ── Extract activation measures ────────────────────────────────────────────────
for i in "${!COPE_LABELS[@]}"; do
    COPE_LABEL="${COPE_LABELS[$i]}"
    N="${ZSTAT_NUMS[$i]}"
    ZSTAT="${FEAT_DIR}/stats/zstat${N}.nii.gz"
    OUT_FILE="${OUT_DIR}/${SUB}_${SES}_${TASK}_${PIPELINE}_${COPE_LABEL}_activation.txt"

    if [ ! -f "$ZSTAT" ]; then
        echo "WARNING: zstat${N} not found, writing NA: ${ZSTAT}"
        cat > "$OUT_FILE" << EOF
subject: ${SUB}
session: ${SES}
task: ${TASK}
pipeline: ${PIPELINE}
cope: ${COPE_LABEL}
mean_z: NA
n_suprathreshold: NA
mean_top1pct_z: NA
p99_z: NA
EOF
        continue
    fi

    # Mean Z and N suprathreshold voxels — threshold at Z>3.1
    TMP_THRESH=$(mktemp /tmp/thresh_XXXXXX)
    fslmaths "$ZSTAT" -thr 3.1 "$TMP_THRESH"
    MEAN_Z=$(fslstats "${TMP_THRESH}.nii.gz" -M)
    NVOX=$(fslstats "${TMP_THRESH}.nii.gz" -V | awk '{print $1}')
    rm -f ${TMP_THRESH}*

    # Top 1% Z — from unthresholded image
    P99=$(fslstats "$ZSTAT" -P 99 | awk '{print $1}')
    TMP_TOP=$(mktemp /tmp/top1pct_XXXXXX)
    fslmaths "$ZSTAT" -thr "$P99" "$TMP_TOP"
    MEAN_TOP1=$(fslstats "${TMP_TOP}.nii.gz" -M)
    rm -f ${TMP_TOP}*

    cat > "$OUT_FILE" << EOF
subject: ${SUB}
session: ${SES}
task: ${TASK}
pipeline: ${PIPELINE}
cope: ${COPE_LABEL}
mean_z: ${MEAN_Z}
n_suprathreshold: ${NVOX}
mean_top1pct_z: ${MEAN_TOP1}
p99_z: ${P99}
EOF

    echo "Done: ${SUB} ${SES} ${TASK} ${PIPELINE} ${COPE_LABEL} | mean_z=${MEAN_Z} nvox=${NVOX} top1pct=${MEAN_TOP1}"
done
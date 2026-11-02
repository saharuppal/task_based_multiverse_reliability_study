# Set output directories 
OUT_DIR="/path/to/project/group_levels/dice_coefficients/${TASK}/${PIPELINE}/${COPE}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# First binarise the thresholded zstat images:
fslmaths \
/path/to/project/group_levels/means/${TASK}/ses-1/${PIPELINE}.gfeat/${COPE}.feat/thresh_zstat1.nii.gz \
-bin \
thresh_zstat1_ses-1_bin.nii.gz

fslmaths \
/path/to/project/group_levels/means/${TASK}/ses-2/${PIPELINE}.gfeat/${COPE}.feat/thresh_zstat1.nii.gz \
-bin \
thresh_zstat1_ses-2_bin.nii.gz

# Multiply one by the other:
fslmaths thresh_zstat1_ses-1_bin.nii.gz \
-mul thresh_zstat1_ses-2_bin.nii.gz \
${TASK}_${COPE}_${PIPELINE}_overlap.nii.gz

# Count number of voxels in visit 1 and 2, plus overlap image
fslstats thresh_zstat1_ses-1_bin.nii.gz -V > ${TASK}_${COPE}_${PIPELINE}_ses-1_volume.txt
fslstats thresh_zstat1_ses-2_bin.nii.gz -V > ${TASK}_${COPE}_${PIPELINE}_ses-2_volume.txt
fslstats ${TASK}_${COPE}_${PIPELINE}_overlap.nii.gz -V > ${TASK}_${COPE}_${PIPELINE}_overlap_volume.txt

# Then calculate dice by 2*overlap/sum of visit 1 and visit 2
V1=$(awk '{print $1}' ${TASK}_${COPE}_${PIPELINE}_ses-1_volume.txt)
V2=$(awk '{print $1}' ${TASK}_${COPE}_${PIPELINE}_ses-2_volume.txt)
OV=$(awk '{print $1}' ${TASK}_${COPE}_${PIPELINE}_overlap_volume.txt)
DICE=$(echo "scale=6; (2*$OV)/($V1+$V2)" | bc -l)

echo "Dice: $DICE"
echo "Dice: $DICE" > ${TASK}_${COPE}_${PIPELINE}_dice_coefficient.txt
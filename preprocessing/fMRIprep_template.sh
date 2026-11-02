#!/bin/bash

# ---- EDIT THESE ----
BASEDIR=/path/to/project                       
SIF=/path/to/fmriprep-24.1.1.simg              
LICENSE_DIR=/path/to/freesurfer_license_dir    
# --------------------

SUB=$1

export SINGULARITYENV_FS_LICENSE=/license/license.txt

if [ ! -f ${BASEDIR}/${SUB}.html ]; then

singularity run -B ${BASEDIR}:/data -B $TMPDIR:/tmp \
-B ${BASEDIR}:/home/fmriprep --home /home/fmriprep \
-B ${LICENSE_DIR}:/license \
--cleanenv ${SIF} \
/data/bids /data/fmriprep_MW_new -w /data/Work participant \
--participant-label $SUB \
--write-graph \
--output-spaces func MNI152NLin6Asym:res-2 \
--skip-bids-validation \
--ignore slicetiming \
--fs-license-file $SINGULARITYENV_FS_LICENSE

fi
#!/bin/bash

module load FSL/6.0.7.17
source "$FSLDIR/etc/fslconf/fsl.sh"

datadir=/path/to/first_levels_output

for TASK in gonogo nback oasis mid
do

for SUB in sub-001 sub-002 sub-003 sub-004 sub-005 sub-006 sub-007 sub-008 sub-009 sub-011 sub-012 sub-013 sub-014 sub-015 sub-016 sub-017 sub-018 sub-019 sub-020 sub-021 sub-022 sub-024 sub-025 sub-026 sub-027
do

for SESSION in ses-1 ses-2
do

for PIPELINE in 6HMP 24HMP 6HMP_aCompCor 6HMP_spike_reg 24HMP_aCompCor  6HMP_WM_CSF 24HMP_WM_CSF \
 24HMP_spike_reg 6HMP_aCompCor_spike_reg 24HMP_aCompCor_spike_reg 6HMP_WM_CSF_spike_reg 24HMP_WM_CSF_spike_reg
do

mkdir $datadir/$TASK/$SUB/$SESSION/${SUB}_${SESSION}_${TASK}_${PIPELINE}.feat/reg

cp $FSLDIR/etc/flirtsch/ident.mat $datadir/$TASK/$SUB/$SESSION/${SUB}_${SESSION}_${TASK}_${PIPELINE}.feat/reg/example_func2standard.mat

cp $datadir/$TASK/$SUB/$SESSION/${SUB}_${SESSION}_${TASK}_${PIPELINE}.feat/mean_func.nii.gz $datadir/$TASK/$SUB/$SESSION/${SUB}_${SESSION}_${TASK}_${PIPELINE}.feat/reg/standard.nii.gz

done 

done 

done

done

exit
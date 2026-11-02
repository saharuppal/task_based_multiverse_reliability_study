#!/bin/bash
set -euo pipefail
shopt -s nullglob

datadir=/path/to/fmriprep_output

SESSIONS=(ses-1 ses-2)
TASKS=(gonogo nback oasis mid rest)

PIPELINES=(
  6HMP 24HMP 6HMP_aCompCor 6HMP_spike_reg 24HMP_aCompCor 6HMP_WM_CSF 24HMP_WM_CSF
  24HMP_spike_reg 6HMP_aCompCor_spike_reg 24HMP_aCompCor_spike_reg 6HMP_WM_CSF_spike_reg 24HMP_WM_CSF_spike_reg
)

# Regressor building blocks
HMP6=(trans_x trans_y trans_z rot_x rot_y rot_z)

HMP24=(
  trans_x trans_y trans_z rot_x rot_y rot_z
  trans_x_derivative1 trans_y_derivative1 trans_z_derivative1 rot_x_derivative1 rot_y_derivative1 rot_z_derivative1
  trans_x_power2 trans_y_power2 trans_z_power2 rot_x_power2 rot_y_power2 rot_z_power2
  trans_x_derivative1_power2 trans_y_derivative1_power2 trans_z_derivative1_power2
  rot_x_derivative1_power2 rot_y_derivative1_power2 rot_z_derivative1_power2
)

ACOMP=(
  w_comp_cor_00 w_comp_cor_01 w_comp_cor_02 w_comp_cor_03 w_comp_cor_04
  c_comp_cor_00 c_comp_cor_01 c_comp_cor_02 c_comp_cor_03 c_comp_cor_04
)

WMCSF=(
  csf white_matter 
)

for SUBDIR in "$datadir"/sub-*; do
  [[ -d "$SUBDIR" ]] || continue
  SUB=$(basename "$SUBDIR")

  for SESSION in "${SESSIONS[@]}"; do
    OUT_DIR="$datadir/$SUB/$SESSION/func"
    [[ -d "$OUT_DIR" ]] || continue

    for TASK in "${TASKS[@]}"; do
      # All matching confounds files (handles run-1/run-2, acq-*, etc.)
      confound_files=( "$OUT_DIR/${SUB}_${SESSION}_task-${TASK}"*_desc-confounds_timeseries.tsv )

      if (( ${#confound_files[@]} == 0 )); then
        echo "SKIP: No confounds found for $SUB $SESSION task-$TASK"
        continue
      fi

      for CONFOUND_FILE in "${confound_files[@]}"; do
        FILE_BASE=$(basename "$CONFOUND_FILE" "_desc-confounds_timeseries.tsv")

        # Read header once and build a name->column map
        IFS=$'\t' read -r -a headers < <(head -n 1 "$CONFOUND_FILE")
        declare -A col=()
        for i in "${!headers[@]}"; do
          col["${headers[$i]}"]=$((i+1))
        done

        for PIPELINE in "${PIPELINES[@]}"; do
          FINAL_OUT="${OUT_DIR}/${FILE_BASE}_${PIPELINE}_no_header.txt"
          WITH_HEADER="${OUT_DIR}/${FILE_BASE}_${PIPELINE}.txt"

          # NO OVERWRITE: skip if final already exists
          if [[ -f "$FINAL_OUT" ]]; then
            echo "EXISTS: $FINAL_OUT"
            continue
          fi

          # Build regressor list for this pipeline
          unset pipeline_regressors temp_array temp_files
          declare -a pipeline_regressors=()
          declare -a temp_array=()
          declare -a temp_files=()

          case "$PIPELINE" in
            6HMP|6HMP_spike_reg)
              pipeline_regressors=("${HMP6[@]}")
              ;;
            24HMP|24HMP_spike_reg)
              pipeline_regressors=("${HMP24[@]}")
              ;;
            6HMP_aCompCor|6HMP_aCompCor_spike_reg)
              pipeline_regressors=("${HMP6[@]}" "${ACOMP[@]}")
              ;;
            24HMP_aCompCor|24HMP_aCompCor_spike_reg)
              pipeline_regressors=("${HMP24[@]}" "${ACOMP[@]}")
              ;;
            6HMP_WM_CSF|6HMP_WM_CSF_spike_reg)
              pipeline_regressors=("${HMP6[@]}" "${WMCSF[@]}")
              ;;
            24HMP_WM_CSF|24HMP_WM_CSF_spike_reg)
              pipeline_regressors=("${HMP24[@]}" "${WMCSF[@]}")
              ;;
            *)
              echo "WARNING: Unknown pipeline $PIPELINE (skipping)" >&2
              continue
              ;;
          esac

          # Add motion outliers for spike_reg pipelines
          if [[ "$PIPELINE" == *"spike_reg"* ]]; then
            for h in "${headers[@]}"; do
              [[ "$h" == motion_outlier* ]] && temp_array+=("$h")
            done
            pipeline_regressors+=("${temp_array[@]}")
          fi

          echo "MAKE: $FILE_BASE | $PIPELINE"

          # Extract each regressor column to a temp file
          for regressor in "${pipeline_regressors[@]}"; do
            idx="${col[$regressor]:-}"
            if [[ -z "$idx" ]]; then
              echo "  WARNING: $regressor not found in $FILE_BASE" >&2
              continue
            fi

            tmp="${OUT_DIR}/tmp_${FILE_BASE}_${PIPELINE}_${regressor}.txt"
            cut -f "$idx" "$CONFOUND_FILE" > "$tmp"
            temp_files+=("$tmp")
          done

          if (( ${#temp_files[@]} == 0 )); then
            echo "  ERROR: No columns extracted for $FILE_BASE $PIPELINE (skipping)" >&2
            continue
          fi

          # Combine, remove header row, replace n/a, cleanup temps
          paste "${temp_files[@]}" > "$WITH_HEADER"
          tail -n +2 "$WITH_HEADER" > "$FINAL_OUT"
          sed -i 's/n\/a/0/g' "$FINAL_OUT"
          rm -f "${temp_files[@]}"
        done
      done
    done
  done
done

echo "Done."
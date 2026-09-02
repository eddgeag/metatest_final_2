#!/usr/bin/env bash
# Backend nativo: usa dependencias locales existentes, no instala nada.
set -Eeuo pipefail
REPO="$METATEST_REPO_ROOT"
OUT="$METATEST_RESULTS_DIR"
REPORT="$METATEST_RUN_DIR"
ARTIFACTS="$REPORT/job_artifacts"
mkdir -p "$ARTIFACTS" "$REPORT/job_report"
start="$(date +%s)"
final=NOT_FINISHED
finish() {
  rc=$?
  trap - EXIT
  {
    printf 'METATEST_PORTABLE_RUN=%s\n' "$([[ $rc == 0 ]] && echo PASS || echo FAIL)"
    printf 'EXIT_STATUS=%s\nFINAL_STATUS=%s\nSCENARIO=%s\nSTAN_O=%s\nRESUME=%s\nWALLTIME_SECONDS=%s\nBACKEND=local\n'       "$rc" "$final" "$METATEST_SCENARIO_OVERRIDE" "$METATEST_STAN_O" "$METATEST_RESUME" "$(( $(date +%s)-start ))"
    printf 'NATIVE_ENVIRONMENT_NOT_CERTIFIED=TRUE\n'
    printf 'NATIVE_RUNTIME_FINGERPRINT=%s\n' "${METATEST_NATIVE_FINGERPRINT:-UNKNOWN}"
    printf 'CODE_COMMIT=%s\n' "${METATEST_CODE_COMMIT:-UNKNOWN}"
    printf 'RUNTIME_BASE_COMMIT=%s\n' "${METATEST_RUNTIME_BASE_COMMIT:-UNKNOWN}"
  } > "$REPORT/job_report/99_FINAL_STATUS.env"
  exit "$rc"
}
trap finish EXIT
[[ -x /usr/bin/time ]]
source "$REPO/scripts/portable/native_env.sh"
metatest_native_prepare "$REPO" "$ARTIFACTS"
if metatest_native_smoke "$REPO" "$ARTIFACTS" > "$ARTIFACTS/01_native_identity.txt" 2>&1; then
  cat "$ARTIFACTS/01_native_identity.txt"
else
  cat "$ARTIFACTS/01_native_identity.txt" >&2
  exit 1
fi
metatest_native_guard_resume "$OUT/.portable_state/$METATEST_SCENARIO_OVERRIDE" "$METATEST_RESUME"
printf 'NATIVE_RUNTIME_FINGERPRINT=%s\n' "$METATEST_NATIVE_FINGERPRINT" > "$ARTIFACTS/native_runtime.env"
export DATA_ROOT="$METATEST_DATA_DIR" RESULTS_ROOT="$OUT"
export COMPARISON_ROOT="$OUT/comparacion_modelamientos"
export ONLY_ANALYSIS_NAME="$METATEST_SCENARIO_OVERRIDE"
export PIPELINE_R_CORES=2 BRMS_INTERNAL_CORES=1 BLAS_THREADS=1
export CMDSTAN_OPTIMIZATION_LEVEL="$METATEST_STAN_O"
export METATEST_PHASE_TIMINGS="$ARTIFACTS/phase_timings.tsv"
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
printf 'scenario\tmodel\tstage\tstart_epoch\tend_epoch\tseconds\texit_code\n' > "$METATEST_PHASE_TIMINGS"
printf 'PIPELINE_START_UTC=%s\n' "$(date -u +%FT%TZ)" > "$ARTIFACTS/04_timing.env"
set +e
/usr/bin/time -v -o "$ARTIFACTS/05_time_and_stderr.txt"   bash "$REPO/scripts/portable/orchestrator.sh" "$REPO/$METATEST_TEMPLATE_REL"   2>&1 | tee "$REPORT/job_report/03_pipeline_console.log"
rc="${PIPESTATUS[0]}"
set -e
printf 'PIPELINE_END_UTC=%s\nPIPELINE_EXIT_STATUS=%s\n' "$(date -u +%FT%TZ)" "$rc" >> "$ARTIFACTS/04_timing.env"
((rc == 0)) || exit "$rc"
source "$REPO/scripts/portable/validate.sh"
for model in modelamiento2 modelamiento3; do
  dir="$COMPARISON_ROOT/$model/$METATEST_SCENARIO_ROOT_OVERRIDE/$METATEST_SCENARIO_OVERRIDE"
  validate_framework_dir "$dir"
  [[ -f "$dir/.FRAMEWORK_BIDIRECCIONAL_COMPLETE" ]]
done
printf 'PILOT_FUNCTIONAL_TEST=PASS\nSELECTED_QUERY_COLUMN=PASS\n' > "$ARTIFACTS/06_pilot_validation.txt"
printf 'FINAL_VALIDATION=PASS\n' > "$ARTIFACTS/07_final_validation.txt"
final=PORTABLE_FINAL_VALIDATION_PASS

#!/usr/bin/env bash
set -Eeuo pipefail
umask 022

# Runner SLURM. Las herramientas científicas (incluido Rscript) se invocan
# exclusivamente dentro del contenedor creado desde la SQSH sellada.

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

need_host() {
  command -v "$1" >/dev/null 2>&1 || die "falta comando del nodo: $1"
}

sha_file() {
  sha256sum "$1" | awk '{print $1}'
}

# Compatible con el Git antiguo del nodo login de CEDIA, que no admite git -C.
git_at() {
  local directory="$1"
  shift
  (cd "$directory" && git "$@")
}

[[ -n "${SLURM_JOB_ID:-}" ]] || die "este script debe ejecutarse mediante SLURM"
for required_var in METATEST_MODE METATEST_REPO_ROOT METATEST_ASSET_DIR METATEST_DATA_DIR METATEST_RESULTS_DIR METATEST_RUN_DIR; do
  [[ -n "${!required_var:-}" ]] || die "falta variable: $required_var"
done
case "$METATEST_MODE" in q00|all) ;; *) die "modo inválido: $METATEST_MODE" ;; esac

for command_name in enroot unsquashfs sha256sum awk grep git tar stat date mkdir mktemp find ln tee df chmod; do
  need_host "$command_name"
done

REPO_ROOT="$METATEST_REPO_ROOT"
ASSET_DIR="$METATEST_ASSET_DIR"
DATA_DIR="$METATEST_DATA_DIR"
RESULTS_DIR="$METATEST_RESULTS_DIR"
RUN_DIR="$METATEST_RUN_DIR"
CONFIG="${REPO_ROOT}/config/cedia_release.env"
[[ -s "$CONFIG" ]] || die "falta configuración: $CONFIG"
# shellcheck disable=SC1090
source "$CONFIG"

IMAGE="${ASSET_DIR}/${METATEST_IMAGE_FILENAME}"
BACKEND_ARCHIVE="${ASSET_DIR}/${METATEST_BACKEND_ARCHIVE_FILENAME}"
INPUT_MANIFEST="${REPO_ROOT}/reproducibilidad/release/inputs.sha256"
PIPELINE_RESULTS="${RUN_DIR}/pipeline_results"
ARTIFACTS="${PIPELINE_RESULTS}/job_artifacts"
JOB_REPORT="${RUN_DIR}/job_report"
TMP_DIR=""
ENROOT_NAME="metatest-release-${SLURM_JOB_ID}"
ENROOT_DATA_ROOT="${METATEST_ENROOT_DATA_ROOT:-${HOME}/.local/share/enroot}"
ENROOT_ROOTFS="${ENROOT_DATA_ROOT}/${ENROOT_NAME}"
START_EPOCH="$(date -u +%s)"
FINAL_STATUS="NOT_FINISHED"
FAILED_LINE=""
FAILED_COMMAND=""

mkdir -p "$PIPELINE_RESULTS" "$ARTIFACTS" "$JOB_REPORT"

on_err() {
  local status="$1" line="$2" command="$3"
  if [[ -z "$FAILED_LINE" ]]; then
    FAILED_LINE="$line"
    FAILED_COMMAND="$command"
  fi
  return "$status"
}

cleanup() {
  local cleanup_status=0
  if enroot list 2>/dev/null | grep -Fxq "$ENROOT_NAME"; then
    enroot remove --force "$ENROOT_NAME" > "${JOB_REPORT}/98_enroot_remove.log" 2>&1 || cleanup_status=1
  fi
  if [[ -n "$TMP_DIR" && -d "$TMP_DIR" ]]; then
    case "$TMP_DIR" in
      "${RUN_DIR}/.metatest_runtime_${SLURM_JOB_ID}."*)
        chmod -R u+rwX "$TMP_DIR" 2>/dev/null || true
        find "$TMP_DIR" -depth -delete || cleanup_status=1
        ;;
      *) printf 'REFUSING_UNSAFE_TMP_CLEANUP=%s\n' "$TMP_DIR" >&2; cleanup_status=1 ;;
    esac
  fi
  return "$cleanup_status"
}

finish() {
  local status=$?
  trap - ERR EXIT HUP INT TERM
  set +e
  cleanup
  cleanup_status=$?
  (( status == 0 && cleanup_status != 0 )) && status=70
  END_EPOCH="$(date -u +%s)"
  command -v sstat >/dev/null 2>&1 && \
    sstat -j "${SLURM_JOB_ID}.batch" --format=JobID,MaxRSS,AveRSS,MaxVMSize,AveCPU -P \
      > "${JOB_REPORT}/96_sstat.txt" 2>&1 || true
  command -v sacct >/dev/null 2>&1 && \
    sacct -j "$SLURM_JOB_ID" --format=JobID,JobName,State,ExitCode,Elapsed,Timelimit,MaxRSS,AllocCPUS,NodeList -P \
      > "${JOB_REPORT}/97_sacct.txt" 2>&1 || true
  {
    printf 'METATEST_RELEASE_RUN=%s\n' "$([[ $status -eq 0 ]] && printf PASS || printf FAIL)"
    printf 'EXIT_STATUS=%s\n' "$status"
    printf 'FINAL_STATUS=%s\n' "$FINAL_STATUS"
    printf 'FAILED_LINE=%s\n' "$FAILED_LINE"
    printf 'FAILED_COMMAND=%s\n' "$FAILED_COMMAND"
    printf 'JOB_ID=%s\n' "$SLURM_JOB_ID"
    printf 'MODE=%s\n' "$METATEST_MODE"
    printf 'RESUME=%s\n' "${METATEST_RESUME:-0}"
    printf 'RUN_DIR=%s\n' "$RUN_DIR"
    printf 'START_UTC=%s\n' "$(date -u -d "@${START_EPOCH}" +%Y-%m-%dT%H:%M:%SZ)"
    printf 'END_UTC=%s\n' "$(date -u -d "@${END_EPOCH}" +%Y-%m-%dT%H:%M:%SZ)"
    printf 'WALLTIME_SECONDS=%s\n' "$((END_EPOCH - START_EPOCH))"
    printf 'IMAGE_SHA256=%s\n' "$METATEST_IMAGE_SHA256"
    printf 'BACKEND_ARCHIVE_SHA256=%s\n' "$METATEST_BACKEND_ARCHIVE_SHA256"
    printf 'HOST_RSCRIPT_REQUIRED=FALSE\n'
    printf 'IMAGE_LAUNCH_MODE=ENROOT_CREATE_THEN_START\n'
    printf 'DATA_MOUNT=READ_ONLY\n'
    printf 'BACKEND_PAYLOAD_MOUNT=READ_ONLY\n'
    printf 'BASILISK_STATE=WRITABLE_EPHEMERAL\n'
  } > "${JOB_REPORT}/99_FINAL_STATUS.env"
  exit "$status"
}

trap 'on_err "$?" "$LINENO" "$BASH_COMMAND"' ERR
trap finish EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

printf '============================================================\n'
printf 'METATEST RELEASE — %s\n' "$METATEST_MODE"
printf '============================================================\n'
printf 'JOB_ID=%s\n' "$SLURM_JOB_ID"
printf 'UTC_START=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'RELEASE_ID=%s\n' "$METATEST_RELEASE_ID"
printf 'VALIDATED_Q00_REFERENCE_JOB=%s\n' "$METATEST_VALIDATED_Q00_JOB"
printf 'RUN_DIR=%s\n' "$RUN_DIR"
printf 'HOST_RSCRIPT_REQUIRED=FALSE\n'

[[ -d "$REPO_ROOT" && -d "$DATA_DIR" && -d "$RESULTS_DIR" ]] || die "faltan rutas de ejecución"
[[ -s "$IMAGE" && -s "$BACKEND_ARCHIVE" ]] || die "faltan artefactos sellados"
[[ "$(stat -c '%s' "$IMAGE")" == "$METATEST_IMAGE_BYTES" ]] || die "tamaño de imagen inesperado"
[[ "$(stat -c '%s' "$BACKEND_ARCHIVE")" == "$METATEST_BACKEND_ARCHIVE_BYTES" ]] || die "tamaño de backend inesperado"
[[ "$(sha_file "$IMAGE")" == "$METATEST_IMAGE_SHA256" ]] || die "SHA256 de imagen inesperado"
[[ "$(sha_file "$BACKEND_ARCHIVE")" == "$METATEST_BACKEND_ARCHIVE_SHA256" ]] || die "SHA256 de backend inesperado"
unsquashfs -s "$IMAGE" > "${JOB_REPORT}/00_unsquashfs_superblock.txt" 2>&1 || die "imagen SquashFS ilegible"

[[ "$(sha_file "${REPO_ROOT}/renv.lock")" == "$METATEST_RENV_LOCK_SHA256" ]] || die "renv.lock inesperado"
[[ "$(sha_file "${REPO_ROOT}/renv/settings.json")" == "$METATEST_RENV_SETTINGS_SHA256" ]] || die "renv/settings.json inesperado"
[[ "$(sha_file "${REPO_ROOT}/reproducibilidad/release/renv.lock.validated")" == "$METATEST_RENV_LOCK_SHA256" ]] || die "copia validada del lock inesperada"
[[ "$(sha_file "${REPO_ROOT}/reproducibilidad/release/renv.settings.validated.json")" == "$METATEST_RENV_SETTINGS_SHA256" ]] || die "copia validada de settings inesperada"
(
  cd "$DATA_DIR"
  sha256sum -c "$INPUT_MANIFEST"
) > "${JOB_REPORT}/01_inputs_sha256_check.txt" 2>&1 || die "hash de input inesperado"

HEAD_COMMIT="$(git_at "$REPO_ROOT" rev-parse HEAD)"
[[ -z "$(git_at "$REPO_ROOT" status --porcelain --untracked-files=all)" ]] || die "worktree de release sucio al iniciar"
[[ "$(git_at "$REPO_ROOT" rev-parse "${METATEST_RELEASE_TAG}^{commit}")" == "$HEAD_COMMIT" ]] || die "tag de release no apunta a HEAD"
[[ "$(git_at "$REPO_ROOT" merge-base "$METATEST_VALIDATED_BASE_COMMIT" "$HEAD_COMMIT")" == "$METATEST_VALIDATED_BASE_COMMIT" ]] \
  || die "base validada ausente"

FREE_GB="$(df -Pk "$RUN_DIR" | awk 'END {print int($4/1024/1024)}')"
[[ "$FREE_GB" =~ ^[0-9]+$ ]] || die "no se pudo medir disco"
(( FREE_GB >= METATEST_DEFAULT_MIN_FREE_GB )) || die "disco libre insuficiente: ${FREE_GB} GiB"

TMP_DIR="$(mktemp -d "${RUN_DIR}/.metatest_runtime_${SLURM_JOB_ID}.XXXXXX")"
printf '[01] Extrayendo backend sellado...\n'
tar -C "$TMP_DIR" -xf "$BACKEND_ARCHIVE"
BACKEND_PAYLOAD="${TMP_DIR}/${METATEST_BACKEND_DIRECTORY_NAME}"
[[ -x "${BACKEND_PAYLOAD}/1.18.0/MOFA2/1.16.0/mofa_env/bin/python" ]] || die "backend extraído incompleto"

BASILISK_STATE="${TMP_DIR}/basilisk_external_state"
mkdir -p "${BASILISK_STATE}/1.18.0/MOFA2/1.16.0"
ln -s /workspace/basilisk_payload/1.18.0/0 "${BASILISK_STATE}/1.18.0/0"
ln -s /workspace/basilisk_payload/1.18.0/MOFA2/1.16.0/mofa_env \
  "${BASILISK_STATE}/1.18.0/MOFA2/1.16.0/mofa_env"

enroot list | grep -Fxq "$ENROOT_NAME" && die "contenedor temporal ya existente: $ENROOT_NAME"
printf '[02] Creando Enroot desde SQSH (sin squashfuse)...\n'
ENROOT_MAX_PROCESSORS="${SLURM_CPUS_PER_TASK:-$METATEST_DEFAULT_SLURM_CPUS}" \
  enroot create --name "$ENROOT_NAME" "$IMAGE" > "${JOB_REPORT}/02_enroot_create.log" 2>&1
enroot list | grep -Fxq "$ENROOT_NAME" || die "enroot create no registró el contenedor"
[[ -d "$ENROOT_ROOTFS" ]] || die "no se localizó rootfs: $ENROOT_ROOTFS"
mkdir -p \
  "${ENROOT_ROOTFS}/workspace/metatest" \
  "${ENROOT_ROOTFS}/workspace/data" \
  "${ENROOT_ROOTFS}/workspace/results" \
  "${ENROOT_ROOTFS}/workspace/basilisk_payload" \
  "${ENROOT_ROOTFS}/workspace/basilisk_external"

printf '[03] Ejecutando smoke, preflight y pipeline dentro de Enroot...\n'
set +e
ENROOT_MOUNT_HOME=no enroot start --root \
  --mount "${REPO_ROOT}:/workspace/metatest:none:bind,rw" \
  --mount "${DATA_DIR}:/workspace/data:none:bind,ro" \
  --mount "${PIPELINE_RESULTS}:/workspace/results:none:bind,rw" \
  --mount "${BACKEND_PAYLOAD}:/workspace/basilisk_payload:none:bind,ro" \
  --mount "${BASILISK_STATE}:/workspace/basilisk_external:none:bind,rw" \
  "$ENROOT_NAME" \
  /usr/bin/env -i \
  HOME=/root USER=root LOGNAME=root \
  PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  LANG=C.UTF-8 LC_ALL=C.UTF-8 TZ=Etc/UTC \
  R_LIBS=/opt/metatest/renv-library \
  R_LIBS_USER=/opt/metatest/renv-library \
  R_ENVIRON_USER=/dev/null \
  R_PROFILE_USER=/dev/null \
  RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
  BASILISK_EXTERNAL_DIR=/workspace/basilisk_external \
  BASILISK_EXTERNAL_CONDA=/workspace/basilisk_external/1.18.0/0 \
  BASILISK_USE_SYSTEM_DIR=0 \
  BASILISK_NO_DESTROY=1 \
  BASILISK_NO_FALLBACK_R=1 \
  RETICULATE_PYTHON=/workspace/basilisk_external/1.18.0/MOFA2/1.16.0/mofa_env/bin/python \
  RETICULATE_USE_MANAGED_VENV=no \
  PYTHONDONTWRITEBYTECODE=1 \
  PIP_NO_INDEX=1 PIP_DISABLE_PIP_VERSION_CHECK=1 \
  CONDA_OFFLINE=true UV_OFFLINE=1 \
  http_proxy=http://127.0.0.1:9 https_proxy=http://127.0.0.1:9 \
  HTTP_PROXY=http://127.0.0.1:9 HTTPS_PROXY=http://127.0.0.1:9 \
  ALL_PROXY=http://127.0.0.1:9 NO_PROXY=localhost,127.0.0.1 \
  CMDSTAN="$METATEST_CMDSTAN_PATH" CMDSTAN_PATH="$METATEST_CMDSTAN_PATH" \
  METATEST_MODE="$METATEST_MODE" \
  METATEST_SCENARIO="$METATEST_Q00_SCENARIO" \
  METATEST_HEAD_COMMIT="$HEAD_COMMIT" \
  EXPECTED_LOCK="$METATEST_RENV_LOCK_SHA256" \
  EXPECTED_SETTINGS="$METATEST_RENV_SETTINGS_SHA256" \
  EXPECTED_R="$METATEST_R_VERSION" \
  EXPECTED_RENV="$METATEST_RENV_VERSION" \
  EXPECTED_BIOC="$METATEST_BIOCONDUCTOR_VERSION" \
  EXPECTED_MOFA2="$METATEST_MOFA2_VERSION" \
  EXPECTED_BASILISK="$METATEST_BASILISK_VERSION" \
  EXPECTED_RETICULATE="$METATEST_RETICULATE_VERSION" \
  EXPECTED_MOFAPY2="$METATEST_MOFAPY2_VERSION" \
  EXPECTED_CMDSTAN="$METATEST_CMDSTAN_VERSION" \
  EXPECTED_CMDSTANR="$METATEST_CMDSTANR_VERSION" \
  EXPECTED_CMDSTANR_SHA="$METATEST_CMDSTANR_REMOTE_SHA" \
  EXPECTED_GGPLOT2="$METATEST_GGPLOT2_VERSION" \
  PIPELINE_R_CORES="$METATEST_DEFAULT_PIPELINE_R_CORES" \
  BRMS_INTERNAL_CORES="$METATEST_DEFAULT_BRMS_INTERNAL_CORES" \
  BLAS_THREADS="$METATEST_DEFAULT_BLAS_THREADS" \
  MIN_FREE_DISK_GB="$METATEST_DEFAULT_MIN_FREE_GB" \
  /bin/bash --noprofile --norc -s <<'INNER' \
  2>&1 | tee "${JOB_REPORT}/03_enroot_pipeline_console.log"
set -Eeuo pipefail

REPO=/workspace/metatest
DATA_ROOT=/workspace/data
RESULTS_ROOT=/workspace/results
COMPARISON_ROOT="${RESULTS_ROOT}/comparacion_modelamientos"
ARTIFACTS="${RESULTS_ROOT}/job_artifacts"
mkdir -p "$ARTIFACTS"
cd "$REPO"

[[ "$(sha256sum renv.lock | awk '{print $1}')" == "$EXPECTED_LOCK" ]]
[[ "$(sha256sum renv/settings.json | awk '{print $1}')" == "$EXPECTED_SETTINGS" ]]
[[ "$(sha256sum /opt/metatest/project/renv.lock | awk '{print $1}')" == "$EXPECTED_LOCK" ]]
[[ "$(sha256sum /opt/metatest/project/renv/settings.json | awk '{print $1}')" == "$EXPECTED_SETTINGS" ]]
[[ -x "${CMDSTAN}/bin/stanc" ]]
"${CMDSTAN}/bin/stanc" --version

if { : > /workspace/data/.metatest_write_probe; } 2>/dev/null; then
  rm -f /workspace/data/.metatest_write_probe
  echo 'ERROR: data mount escribible' >&2
  exit 21
fi
if { : > /workspace/basilisk_payload/.metatest_write_probe; } 2>/dev/null; then
  rm -f /workspace/basilisk_payload/.metatest_write_probe
  echo 'ERROR: backend payload escribible' >&2
  exit 22
fi
: > /workspace/basilisk_external/.metatest_state_write_probe
rm -f /workspace/basilisk_external/.metatest_state_write_probe
: > "${ARTIFACTS}/.results_write_probe"
rm -f "${ARTIFACTS}/.results_write_probe"
printf 'MOUNT_POLICY=PASS\n' | tee "${ARTIFACTS}/00_mount_policy.txt"

Rscript --vanilla - <<'RS_IDENTITY' | tee "${ARTIFACTS}/01_runtime_identity.txt"
.libPaths(unique(c("/opt/metatest/renv-library", .Library)))
stopifnot(as.character(getRversion()) == Sys.getenv("EXPECTED_R"))
stopifnot(as.character(packageVersion("renv")) == Sys.getenv("EXPECTED_RENV"))
stopifnot(as.character(BiocManager::version()) == Sys.getenv("EXPECTED_BIOC"))
stopifnot(as.character(packageVersion("MOFA2")) == Sys.getenv("EXPECTED_MOFA2"))
stopifnot(as.character(packageVersion("basilisk")) == Sys.getenv("EXPECTED_BASILISK"))
stopifnot(as.character(packageVersion("reticulate")) == Sys.getenv("EXPECTED_RETICULATE"))
stopifnot(as.character(packageVersion("cmdstanr")) == Sys.getenv("EXPECTED_CMDSTANR"))
stopifnot(as.character(packageVersion("ggplot2")) == Sys.getenv("EXPECTED_GGPLOT2"))
desc <- utils::packageDescription("cmdstanr")
stopifnot(identical(unname(desc[["RemoteSha"]]), Sys.getenv("EXPECTED_CMDSTANR_SHA")))
suppressPackageStartupMessages(library(cmdstanr))
set_cmdstan_path(Sys.getenv("CMDSTAN"))
stopifnot(identical(as.character(cmdstan_version()), Sys.getenv("EXPECTED_CMDSTAN")))
suppressPackageStartupMessages(library(MOFA2))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(basilisk))
expected_python <- "/workspace/basilisk_external/1.18.0/MOFA2/1.16.0/mofa_env/bin/python"
stopifnot(identical(normalizePath(py_config()$python), normalizePath(expected_python)))
py_run_string("import mofapy2; assert mofapy2.version.__version__ == '0.7.0'")
env <- get("mofa_env", envir = asNamespace("MOFA2"))
expected_env <- "/workspace/basilisk_external/1.18.0/MOFA2/1.16.0/mofa_env"
stopifnot(identical(normalizePath(basilisk::obtainEnvironmentPath(env)), normalizePath(expected_env)))
proc <- basilisk::basiliskStart(env)
on.exit(basilisk::basiliskStop(proc), add = TRUE)
observed <- basilisk::basiliskRun(proc, function() {
  m <- reticulate::import("mofapy2", convert = FALSE)
  reticulate::py_to_r(m$version$`__version__`)
})
stopifnot(identical(observed, Sys.getenv("EXPECTED_MOFAPY2")))
cat("SEALED_RUNTIME_IDENTITY=PASS\n")
RS_IDENTITY

export REPO_ROOT="$REPO"
export WORKSPACE_ROOT=/workspace
export DATA_ROOT RESULTS_ROOT COMPARISON_ROOT
export PAUSE_FILE="${RESULTS_ROOT}/.server_runs/PAUSAR_PIPELINE"
if [[ "$METATEST_MODE" == q00 ]]; then
  export ONLY_ANALYSIS_NAME="$METATEST_SCENARIO"
else
  unset ONLY_ANALYSIS_NAME
fi

export OMP_NUM_THREADS="$BLAS_THREADS"
export OPENBLAS_NUM_THREADS="$BLAS_THREADS"
export MKL_NUM_THREADS="$BLAS_THREADS"
export BLIS_NUM_THREADS="$BLAS_THREADS"
export VECLIB_MAXIMUM_THREADS="$BLAS_THREADS"
export NUMEXPR_NUM_THREADS="$BLAS_THREADS"
export STAN_NUM_THREADS="$BLAS_THREADS"
export TBB_NUM_THREADS="$BLAS_THREADS"
export MALLOC_ARENA_MAX=2
export MAKEFLAGS=-j1
export CMDSTAN_OPTIMIZATION_LEVEL=0
export O=0
export MC_CORES="$BRMS_INTERNAL_CORES"
export MIN_AVAILABLE_MEM_GB=12
export MEMORY_WAIT_SECONDS=60
export MEMORY_WAIT_ATTEMPTS=0

bash scripts/reproducibility/01_preflight_servidor.sh \
  2>&1 | tee "${ARTIFACTS}/02_preflight.txt"
grep -Fxq 'PRE-FLIGHT BÁSICO OK' "${ARTIFACTS}/02_preflight.txt"
printf 'PREFLIGHT=PASS\n' | tee "${ARTIFACTS}/03_preflight_pass.txt"

printf 'PIPELINE_START_UTC=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  | tee "${ARTIFACTS}/04_timing.env"
printf 'MODE=%s\n' "$METATEST_MODE" | tee -a "${ARTIFACTS}/04_timing.env"
printf 'GIT_HEAD=%s\n' "$METATEST_HEAD_COMMIT" | tee -a "${ARTIFACTS}/04_timing.env"

set +e
if [[ -x /usr/bin/time ]]; then
  /usr/bin/time -v bash scripts/run_comparacion_modelamiento2_vs_3.sh \
    2> >(tee "${ARTIFACTS}/05_time_and_stderr.txt" >&2)
  pipeline_status=$?
else
  bash scripts/run_comparacion_modelamiento2_vs_3.sh
  pipeline_status=$?
fi
set -e
printf 'PIPELINE_END_UTC=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "${ARTIFACTS}/04_timing.env"
printf 'PIPELINE_EXIT_STATUS=%s\n' "$pipeline_status" | tee -a "${ARTIFACTS}/04_timing.env"
(( pipeline_status == 0 )) || exit "$pipeline_status"

validate_heatmap_signal() {
  local dir="$1" signal="$2"
  local manifest="${dir}/diagnostico_bayes_features/final_annotated_bidirectional_module_heatmap_${signal}/tables/00_final_annotated_bidirectional_module_heatmap_manifest_${signal}.csv"
  local signature="${dir}/diagnostico_bayes_features/framework_bidireccional_${signal}/tables/01_bidirectional_feature_signature.csv"
  if [[ -s "$manifest" ]]; then
    return 0
  fi
  [[ -s "$signature" ]] || return 1
  Rscript --vanilla - "$signature" <<'RS_HEAT'
args <- commandArgs(trailingOnly = TRUE)
x <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)
required <- c("evidence_grade", "bidirectional_pattern")
if (length(setdiff(required, colnames(x))) > 0L) quit(status = 1L)
eligible <- x$evidence_grade %in% c(
  "A_strong", "B_moderate", "C_exploratory", "D_ambiguous_dependent"
) & !is.na(x$bidirectional_pattern) & x$bidirectional_pattern != "uncertain"
eligible[is.na(eligible)] <- FALSE
quit(status = if (any(eligible)) 1L else 0L)
RS_HEAT
}

validate_framework_dir() {
  local dir="$1" signal relative id_resolution
  for signal in reconstruido observado; do
    for relative in \
      "framework_bidireccional_${signal}/tables/01_bidirectional_feature_signature.csv" \
      "framework_bidireccional_${signal}/tables/04_bidirectional_final_table_corrected.csv" \
      "framework_bidireccional_${signal}/plots_interpretables_bidireccional/00_plot_manifest.csv"
    do
      [[ -s "${dir}/diagnostico_bayes_features/${relative}" ]] || return 1
    done
    validate_heatmap_signal "$dir" "$signal" || return 1
  done
  [[ -s "${dir}/diagnostico_bayes_features/mofa_dual_feature_differential_ranking/20_observed_all_feature_differential_ranking.csv" ]] || return 1
  [[ -s "${dir}/diagnostico_bayes_features/mofa_dual_feature_differential_ranking/21_reconstructed_all_feature_differential_ranking.csv" ]] || return 1
  [[ -s "${dir}/diagnostico_bayes_features/database_biological_annotation/tables/00_database_annotation_manifest.csv" ]] || return 1
  id_resolution="${dir}/diagnostico_bayes_features/database_biological_annotation/tables/31_tx_pr_identifier_resolution.csv"
  [[ -s "$id_resolution" ]] || return 1
  Rscript --vanilla - "$id_resolution" <<'RS_ID'
args <- commandArgs(trailingOnly = TRUE)
x <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)
quit(status = if ("selected_query" %in% colnames(x)) 0L else 1L)
RS_ID
}

if [[ "$METATEST_MODE" == q00 ]]; then
  for model in modelamiento2 modelamiento3; do
    scenario_dir="${COMPARISON_ROOT}/${model}/escenarios/${METATEST_SCENARIO}"
    validate_framework_dir "$scenario_dir" || {
      echo "ERROR: validación de framework falló: $model" >&2
      exit 51
    }
    [[ -f "${scenario_dir}/.FRAMEWORK_BIDIRECCIONAL_COMPLETE" ]] || exit 52
    printf 'FRAMEWORK_COMPLETE_%s=PASS\n' "$model"
  done
  [[ -f "${COMPARISON_ROOT}/resumen_corridas.tsv" ]] || exit 53
  grep -Fq "$METATEST_SCENARIO" "${COMPARISON_ROOT}/resumen_corridas.tsv" || exit 54
  printf 'HEATMAP_VALIDATION=PASS_OR_JUSTIFIED_ABSENCE\n' | tee "${ARTIFACTS}/06_q00_validation.txt"
  printf 'SELECTED_QUERY_COLUMN=PASS\n' | tee -a "${ARTIFACTS}/06_q00_validation.txt"
  printf 'Q00_FUNCTIONAL_TEST=PASS\n' | tee -a "${ARTIFACTS}/06_q00_validation.txt"
  printf 'FINAL_VALIDATION=PASS\n' | tee "${ARTIFACTS}/07_final_validation.txt"
else
  [[ -f "${COMPARISON_ROOT}/resumen_corridas.tsv" ]] || exit 55
  printf 'ALL_ACTIVE_SCENARIOS_PIPELINE=PASS\n' | tee "${ARTIFACTS}/06_all_validation.txt"
fi
INNER
ENROOT_STATUS="${PIPESTATUS[0]}"
set -e
(( ENROOT_STATUS == 0 )) || die "pipeline dentro de Enroot falló: status=${ENROOT_STATUS}"

[[ -z "$(git_at "$REPO_ROOT" status --porcelain --untracked-files=all)" ]] \
  || die "el pipeline dejó cambios Git no ignorados en el worktree"
FINAL_STATUS="$([[ "$METATEST_MODE" == q00 ]] && printf Q00_FINAL_VALIDATION_PASS || printf ALL_PIPELINE_PASS)"
printf '============================================================\n'
printf 'METATEST_RELEASE_%s=PASS\n' "${METATEST_MODE^^}"
printf 'RUN_DIR=%s\n' "$RUN_DIR"

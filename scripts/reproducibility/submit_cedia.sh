#!/usr/bin/env bash
set -Eeuo pipefail
umask 022

# Interfaz única de envío para el release sellado de METATEST.
# Este script se ejecuta en login y deliberadamente no requiere R/Rscript.

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

need() {
  command -v "$1" >/dev/null 2>&1 || die "falta comando de login: $1"
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

MODE="${1:-}"
case "$MODE" in
  q00|all) ;;
  *)
    printf 'Uso: bash scripts/reproducibility/submit_cedia.sh q00|all\n' >&2
    exit 2
    ;;
esac
[[ $# -eq 1 ]] || die "se acepta exactamente un modo: q00 o all"

for command_name in git sbatch sha256sum awk stat date mkdir find grep; do
  need "$command_name"
done

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." >/dev/null 2>&1 && pwd -P)"
CONFIG="${REPO_ROOT}/config/cedia_release.env"
[[ -s "$CONFIG" ]] || die "falta configuración de release: $CONFIG"
# shellcheck disable=SC1090
source "$CONFIG"

ASSET_DIR="${METATEST_ASSET_DIR:-${REPO_ROOT}/../metatest_release_assets_private/${METATEST_RELEASE_ID}}"
DATA_DIR="${METATEST_DATA_DIR:-${REPO_ROOT}/../data}"
RESULTS_DIR="${METATEST_RESULTS_DIR:-${REPO_ROOT}/../results}"
IMAGE="${ASSET_DIR}/${METATEST_IMAGE_FILENAME}"
BACKEND_ARCHIVE="${ASSET_DIR}/${METATEST_BACKEND_ARCHIVE_FILENAME}"
INPUT_MANIFEST="${REPO_ROOT}/reproducibilidad/release/inputs.sha256"
RUNNER="${REPO_ROOT}/scripts/reproducibility/run_cedia_release.sh"

[[ -d "${REPO_ROOT}/.git" || -f "${REPO_ROOT}/.git" ]] || die "no es un worktree Git: $REPO_ROOT"
[[ -d "$ASSET_DIR" ]] || die "falta METATEST_ASSET_DIR: $ASSET_DIR"
[[ -d "$DATA_DIR" ]] || die "falta METATEST_DATA_DIR: $DATA_DIR"
[[ -d "$RESULTS_DIR" ]] || die "falta METATEST_RESULTS_DIR: $RESULTS_DIR"
[[ -s "$IMAGE" ]] || die "falta imagen sellada: $IMAGE"
[[ -s "$BACKEND_ARCHIVE" ]] || die "falta backend sellado: $BACKEND_ARCHIVE"
[[ -s "$INPUT_MANIFEST" ]] || die "falta manifest de inputs: $INPUT_MANIFEST"
[[ -x "$RUNNER" ]] || die "runner no ejecutable: $RUNNER"

HEAD_COMMIT="$(git_at "$REPO_ROOT" rev-parse HEAD)"
[[ -z "$(git_at "$REPO_ROOT" status --porcelain --untracked-files=all)" ]] \
  || die "el worktree de release no está limpio"
[[ "$(git_at "$REPO_ROOT" rev-parse "${METATEST_RELEASE_TAG}^{commit}" 2>/dev/null || true)" == "$HEAD_COMMIT" ]] \
  || die "el tag ${METATEST_RELEASE_TAG} no apunta a HEAD"
[[ "$(git_at "$REPO_ROOT" merge-base "$METATEST_VALIDATED_BASE_COMMIT" "$HEAD_COMMIT")" == "$METATEST_VALIDATED_BASE_COMMIT" ]] \
  || die "HEAD no desciende del commit científico validado"

[[ "$(sha_file "${REPO_ROOT}/renv.lock")" == "$METATEST_RENV_LOCK_SHA256" ]] \
  || die "renv.lock no coincide con el lock efectivo sellado"
[[ "$(sha_file "${REPO_ROOT}/renv/settings.json")" == "$METATEST_RENV_SETTINGS_SHA256" ]] \
  || die "renv/settings.json no coincide con el runtime sellado"
[[ "$(stat -c '%s' "$IMAGE")" == "$METATEST_IMAGE_BYTES" ]] \
  || die "tamaño inesperado de la imagen sellada"
[[ "$(stat -c '%s' "$BACKEND_ARCHIVE")" == "$METATEST_BACKEND_ARCHIVE_BYTES" ]] \
  || die "tamaño inesperado del backend sellado"

while IFS= read -r manifest_line; do
  [[ -n "$manifest_line" ]] || continue
  input_name="${manifest_line#*  }"
  [[ -s "${DATA_DIR}/${input_name}" ]] || die "falta input declarado: ${DATA_DIR}/${input_name}"
done < "$INPUT_MANIFEST"

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
RESUME="${METATEST_RESUME:-0}"
[[ "$RESUME" == 0 || "$RESUME" == 1 ]] || die "METATEST_RESUME debe ser 0 o 1"

if [[ -n "${METATEST_RUN_DIR:-}" ]]; then
  RUN_DIR="$METATEST_RUN_DIR"
  [[ "$RESUME" == 1 ]] || die "METATEST_RUN_DIR requiere METATEST_RESUME=1"
  [[ -d "$RUN_DIR" ]] || die "no existe METATEST_RUN_DIR: $RUN_DIR"
else
  [[ "$RESUME" == 0 ]] || die "METATEST_RESUME=1 requiere METATEST_RUN_DIR"
  RUN_DIR="${RESULTS_DIR}/reproducibilidad_release/${METATEST_RELEASE_ID}/${STAMP}_${MODE}"
  mkdir -p "$RUN_DIR"
fi

PIPELINE_RESULTS="${RUN_DIR}/pipeline_results"
mkdir -p "${PIPELINE_RESULTS}" "${RUN_DIR}/logs" "${RUN_DIR}/submission"
if [[ "$RESUME" == 0 ]] && find "$PIPELINE_RESULTS" -mindepth 1 -print -quit | grep -q .; then
  die "una ejecución nueva exige pipeline_results vacío: $PIPELINE_RESULTS"
fi

PARTITION="${METATEST_SLURM_PARTITION:-$METATEST_DEFAULT_SLURM_PARTITION}"
CPUS="${METATEST_SLURM_CPUS:-$METATEST_DEFAULT_SLURM_CPUS}"
MEM="${METATEST_SLURM_MEM:-$METATEST_DEFAULT_SLURM_MEM}"
TIME_LIMIT="${METATEST_SLURM_TIME:-$METATEST_DEFAULT_SLURM_TIME}"
[[ "$CPUS" =~ ^[1-9][0-9]*$ ]] || die "METATEST_SLURM_CPUS inválido: $CPUS"

{
  printf 'SUBMIT_UTC=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'MODE=%s\n' "$MODE"
  printf 'RESUME=%s\n' "$RESUME"
  printf 'RELEASE_ID=%s\n' "$METATEST_RELEASE_ID"
  printf 'RELEASE_TAG=%s\n' "$METATEST_RELEASE_TAG"
  printf 'GIT_HEAD=%s\n' "$HEAD_COMMIT"
  printf 'VALIDATED_Q00_JOB=%s\n' "$METATEST_VALIDATED_Q00_JOB"
  printf 'ASSET_DIR=%s\n' "$ASSET_DIR"
  printf 'DATA_DIR=%s\n' "$DATA_DIR"
  printf 'RESULTS_DIR=%s\n' "$RESULTS_DIR"
  printf 'RUN_DIR=%s\n' "$RUN_DIR"
  printf 'PARTITION=%s\n' "$PARTITION"
  printf 'CPUS=%s\n' "$CPUS"
  printf 'MEM=%s\n' "$MEM"
  printf 'TIME=%s\n' "$TIME_LIMIT"
  printf 'HOST_RSCRIPT_REQUIRED=FALSE\n'
} > "${RUN_DIR}/submission/00_request.env"
sha256sum "$CONFIG" "$RUNNER" "$INPUT_MANIFEST" \
  > "${RUN_DIR}/submission/01_release_code.sha256"

EXPORTS="ALL,METATEST_MODE=${MODE},METATEST_RESUME=${RESUME},METATEST_REPO_ROOT=${REPO_ROOT},METATEST_ASSET_DIR=${ASSET_DIR},METATEST_DATA_DIR=${DATA_DIR},METATEST_RESULTS_DIR=${RESULTS_DIR},METATEST_RUN_DIR=${RUN_DIR}"
SBATCH_REPLY="$(
  sbatch --parsable \
    --job-name="metatest_${MODE}" \
    --partition="$PARTITION" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$CPUS" \
    --mem="$MEM" \
    --time="$TIME_LIMIT" \
    --no-requeue \
    --chdir="$REPO_ROOT" \
    --output="${RUN_DIR}/logs/slurm_%j.out" \
    --error="${RUN_DIR}/logs/slurm_%j.err" \
    --export="$EXPORTS" \
    "$RUNNER"
)"
JOB_ID="${SBATCH_REPLY%%;*}"
[[ "$JOB_ID" =~ ^[0-9]+$ ]] || die "respuesta sbatch inesperada: $SBATCH_REPLY"
printf '%s\n' "$JOB_ID" > "${RUN_DIR}/submission/02_job_id.txt"

printf 'JOB_ID=%s\n' "$JOB_ID"
printf 'RUN_DIR=%s\n' "$RUN_DIR"
printf 'STDOUT=%s\n' "${RUN_DIR}/logs/slurm_${JOB_ID}.out"
printf 'STDERR=%s\n' "${RUN_DIR}/logs/slurm_${JOB_ID}.err"
printf '# Ctrl+C detiene tail, no el job.\n'
printf "tail -F '%s' '%s'\n" \
  "${RUN_DIR}/logs/slurm_${JOB_ID}.out" \
  "${RUN_DIR}/logs/slurm_${JOB_ID}.err"
printf "squeue -j '%s' -o '%s'\n" "$JOB_ID" '%.18i %.9P %.24j %.2t %.10M %.10l %.6D %R'

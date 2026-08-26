#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." >/dev/null 2>&1 && pwd -P)"
WORKSPACE_ROOT="$(cd -- "${REPO_ROOT}/.." >/dev/null 2>&1 && pwd -P)"
DATA_ROOT="${DATA_ROOT:-${WORKSPACE_ROOT}/data}"
RESULTS_ROOT="${RESULTS_ROOT:-${WORKSPACE_ROOT}/results}"
COMPARISON_ROOT="${COMPARISON_ROOT:-${RESULTS_ROOT}/comparacion_modelamientos}"
export REPO_ROOT WORKSPACE_ROOT DATA_ROOT RESULTS_ROOT COMPARISON_ROOT
cd "${REPO_ROOT}"

ACTION="${1:-}"
RUN_BASE="${RESULTS_ROOT}/.server_runs"
CURRENT="${RUN_BASE}/current"
PAUSE_FILE="${RUN_BASE}/PAUSAR_PIPELINE"

case "${ACTION}" in
  pausa)
    touch "${PAUSE_FILE}"
    echo "Pausa solicitada. El orquestador terminará al llegar al siguiente límite seguro."
    echo "No libera inmediatamente la RAM del paso R que ya está corriendo."
    ;;
  detener)
    RUN_DIR="$(readlink -f "${CURRENT}")"
    PGID="$(cat "${RUN_DIR}/pipeline.pgid" 2>/dev/null || true)"
    [[ "${PGID}" =~ ^[0-9]+$ ]] || { echo "No se encontró PGID válido."; exit 1; }
    echo "Enviando TERM al grupo ${PGID}..."
    kill -TERM -- "-${PGID}"
    ;;
  reanudar)
    rm -f "${PAUSE_FILE}"
    exec "${SCRIPT_DIR}/03_lanzar_pipeline_servidor.sh" --background
    ;;
  *)
    echo "Uso: $0 {pausa|detener|reanudar}"
    exit 2
    ;;
esac

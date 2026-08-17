#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ACTION="${1:-}"
CURRENT="${ROOT}/.server_runs/current"

case "${ACTION}" in
  pausa)
    touch "${ROOT}/PAUSAR_PIPELINE"
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
    rm -f "${ROOT}/PAUSAR_PIPELINE"
    exec "${ROOT}/03_lanzar_pipeline_servidor.sh" --background
    ;;
  *)
    echo "Uso: $0 {pausa|detener|reanudar}"
    exit 2
    ;;
esac


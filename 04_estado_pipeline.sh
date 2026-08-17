#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CURRENT="${ROOT}/.server_runs/current"

if [[ ! -L "${CURRENT}" && ! -d "${CURRENT}" ]]; then
  echo "No hay ejecución registrada."
  exit 1
fi

RUN_DIR="$(readlink -f "${CURRENT}")"
PID="$(cat "${RUN_DIR}/pipeline.pid" 2>/dev/null || true)"

echo "RUN: ${RUN_DIR}"
echo "PID: ${PID:-NO_REGISTRADO}"
if [[ -n "${PID}" ]] && kill -0 "${PID}" 2>/dev/null; then
  echo "ESTADO: EN EJECUCIÓN"
  ps -o pid,ppid,pgid,%cpu,%mem,etime,rss,cmd -p "${PID}" || true
else
  echo "ESTADO: DETENIDO O TERMINADO"
fi

echo
free -h
echo
df -h "${ROOT}"
echo
echo "Ramas completas:"
find "${ROOT}/comparacion_modelamientos" -name .PIPELINE_COMPLETE -type f 2>/dev/null | wc -l
echo "Frameworks completos:"
find "${ROOT}/comparacion_modelamientos" -name .BIDIR_COMPLETE -type f 2>/dev/null | wc -l
echo
echo "Últimas 50 líneas:"
tail -n 50 "${RUN_DIR}/pipeline_master.log" 2>/dev/null || true


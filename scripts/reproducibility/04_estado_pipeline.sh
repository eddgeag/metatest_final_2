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

RUN_BASE="${RESULTS_ROOT}/.server_runs"
CURRENT="${RUN_BASE}/current"

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
df -h "${RESULTS_ROOT}"
echo
echo "Ramas completas:"
find "${COMPARISON_ROOT}" -name .PIPELINE_COMPLETE -type f 2>/dev/null | wc -l
echo "Frameworks completos:"
find "${COMPARISON_ROOT}" -name .FRAMEWORK_BIDIRECCIONAL_COMPLETE -type f 2>/dev/null | wc -l
echo
echo "Últimas 50 líneas:"
tail -n 50 "${RUN_DIR}/pipeline_master.log" 2>/dev/null || true

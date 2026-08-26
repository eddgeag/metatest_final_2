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

CONFIG_FILE="${CONFIG_FILE:-${REPO_ROOT}/config/recursos_servidor.env}"
[[ -f "${CONFIG_FILE}" ]] && source "${CONFIG_FILE}"

REPORT_DIR="${RESULTS_ROOT}/.server_runs/preflight_$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "${REPORT_DIR}"
REPORT="${REPORT_DIR}/preflight.txt"
FAIL=0
exec > >(tee "${REPORT}") 2>&1

need_cmd () {
  local cmd="$1"
  if command -v "${cmd}" >/dev/null 2>&1; then
    echo "OK command ${cmd}: $(command -v "${cmd}")"
  else
    echo "FALTA command: ${cmd}"
    FAIL=1
  fi
}

echo "PRE-FLIGHT METATEST"
  echo "UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "REPO_ROOT: ${REPO_ROOT}"
  echo "DATA_ROOT: ${DATA_ROOT}"
  echo "RESULTS_ROOT: ${RESULTS_ROOT}"
  echo "COMPARISON_ROOT: ${COMPARISON_ROOT}"
  echo
  echo "== Sistema =="
  uname -a
  [[ -r /etc/os-release ]] && cat /etc/os-release
  echo
  echo "== CPU =="
  nproc
  lscpu 2>/dev/null | grep -E '^(CPU\(s\)|Thread|Core|Socket|Model name)' || true
  echo
  echo "== Memoria =="
  free -h
  swapon --show 2>/dev/null || true
  echo
  echo "== Disco =="
  df -h "${RESULTS_ROOT}"
  echo
  echo "== Scheduler/cgroup =="
  echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
  echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-}"
  echo "SLURM_MEM_PER_NODE=${SLURM_MEM_PER_NODE:-}"
  [[ -r /sys/fs/cgroup/cpu.max ]] && { printf 'cpu.max='; cat /sys/fs/cgroup/cpu.max; }
  [[ -r /sys/fs/cgroup/memory.max ]] && { printf 'memory.max='; cat /sys/fs/cgroup/memory.max; }
  echo
  echo "== Herramientas =="
  for c in bash R Rscript git make gcc g++ python3 awk find sha256sum flock nohup setsid; do
    need_cmd "${c}"
  done
  echo
  echo "== Versiones =="
  R --version | head -n 1 || true
  gcc --version 2>/dev/null | head -n 1 || true
  g++ --version 2>/dev/null | head -n 1 || true
  make --version 2>/dev/null | head -n 1 || true
  git --version 2>/dev/null || true
  echo
  echo "== Git =="
  git status --short 2>/dev/null || true
git describe --always --dirty --tags 2>/dev/null || true

REQUIRED=(
  renv.lock
  renv/activate.R
  scripts/pipeline/preprocess_data.R
  scripts/pipeline/INTEGRACION.R
  scripts/pipeline/modelamiento2.R
  scripts/pipeline/modelamiento3.R
  scripts/pipeline/panel_final_multimodelo.R
  scripts/pipeline/downstream.R
  scripts/pipeline/FRAMEWORK_BIDIRECCIONAL.R
  scripts/run_comparacion_modelamiento2_vs_3.sh
  templates/current/plantilla_1_todas_clinicas_barrido.tsv
  templates/current/plantilla_2_con_FLI_barrido.tsv
  templates/current/plantilla_3_sin_FLI_barrido.tsv
  config/recursos_servidor.env
  reproducibilidad/entorno_R_y_cmdstan.txt
  scripts/utilities/crear_mapa_clinico_maestro.R
  scripts/utilities/crear_mapa_proteico_maestro.R
)

for f in "${REQUIRED[@]}"; do
  if [[ ! -s "${f}" ]]; then
    echo "FALTA archivo requerido: ${f}"
    FAIL=1
  fi
done

for TEMPLATE_FILE in \
  templates/current/plantilla_1_todas_clinicas_barrido.tsv \
  templates/current/plantilla_2_con_FLI_barrido.tsv \
  templates/current/plantilla_3_sin_FLI_barrido.tsv
do
  if ! awk -F '\t' '
    /^#/ || NF==0 {next}
    $1=="ANALYSIS_NAME" {if(NF!=11) exit 2; header++; next}
    {if(NF!=11 || seen[$1]++) exit 3; n++}
    END {if(header!=1 || n!=67) exit 4}
  ' "${TEMPLATE_FILE}"; then
    echo "Plantilla activa inválida (se esperan 67 filas, 11 columnas e IDs únicos): ${TEMPLATE_FILE}"
    FAIL=1
  else
    echo "OK plantilla activa: ${TEMPLATE_FILE}"
  fi
done

if [[ ! -d "${DATA_ROOT}" ]]; then
  echo "FALTA DATA_ROOT: ${DATA_ROOT}"
  FAIL=1
else
  DATA_FILES=(
    "transcriptomica.xlsx"
    "citokina.csv"
    "miokina.csv"
    "5_03. METAHEALTH_POS_3 grupos_filtered_normalized.xlsx"
    "5_04. METAHEALTH_NEG_3 grupos_filtered_normalized.xlsx"
    "clinicos.xlsx"
  )
  for f in "${DATA_FILES[@]}"; do
    if [[ ! -s "${DATA_ROOT}/${f}" ]]; then
      echo "FALTA input requerido: ${DATA_ROOT}/${f}"
      FAIL=1
    fi
  done
fi

if [[ -f reproducibilidad/entorno_R_y_cmdstan.txt ]]; then
  LOCAL_R="$(sed -n 's/^R_VERSION=R version \([^ ]*\).*/\1/p' reproducibilidad/entorno_R_y_cmdstan.txt | head -n1)"
  SERVER_R="$(Rscript --vanilla -e 'cat(as.character(getRversion()))' 2>/dev/null || true)"
  if [[ -n "${LOCAL_R}" && "${LOCAL_R}" != "${SERVER_R}" ]]; then
    echo "ERROR: R local congelado=${LOCAL_R}; R servidor=${SERVER_R}"
    FAIL=1
  fi
fi

FREE_GB="$(df -Pk "${RESULTS_ROOT}" | awk 'NR==2 {print int($4/1024/1024)}')"
MIN_FREE_DISK_GB="${MIN_FREE_DISK_GB:-100}"
if (( FREE_GB < MIN_FREE_DISK_GB )); then
  echo "ERROR: solo ${FREE_GB} GiB libres; mínimo solicitado ${MIN_FREE_DISK_GB} GiB."
  FAIL=1
fi

if (( FAIL != 0 )); then
  echo "PRE-FLIGHT FALLIDO. Corrige lo indicado antes de restaurar o ejecutar."
  exit 1
fi

echo "PRE-FLIGHT BÁSICO OK"
echo "Reporte: ${REPORT}"

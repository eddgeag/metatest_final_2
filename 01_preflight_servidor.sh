#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT}"

REPORT_DIR="${ROOT}/.server_runs/preflight_$(date -u +%Y%m%dT%H%M%SZ)"
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
  echo "ROOT: ${ROOT}"
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
  df -h "${ROOT}"
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
  preprocess_data.R INTEGRACION.R modelamiento2.R modelamiento3.R
  panel_final_multimodelo.R downstream.R FRAMEWORK_BIDIRECCIONAL.R
  run_comparacion_modelamiento2_vs_3_memoria_servidor.sh
  plantillas/plantilla_1_todas_clinicas_barrido.tsv
  plantillas/plantilla_2_con_FLI_barrido.tsv
  plantillas/plantilla_3_sin_FLI_barrido.tsv
)

for f in "${REQUIRED[@]}"; do
  if [[ ! -s "${f}" ]]; then
    echo "FALTA archivo requerido: ${f}"
    FAIL=1
  fi
done

for f in plantillas/plantilla_{1_todas_clinicas,2_con_FLI,3_sin_FLI}_barrido.tsv; do
  [[ -f "${f}" ]] || continue
  if ! awk -F '\t' '
    /^#/ || NF==0 {next}
    $1=="ANALYSIS_NAME" {if(NF!=11) exit 2; header++; next}
    {if(NF!=11 || seen[$1]++) exit 3; n++}
    END {if(header!=1 || n!=67) exit 4}
  ' "${f}"; then
    echo "Plantilla inválida (se esperan 67 escenarios, 11 columnas, IDs únicos): ${f}"
    FAIL=1
  else
    echo "OK plantilla: ${f}"
  fi
done

if [[ -f reproducibilidad/entorno_R_y_cmdstan.txt ]]; then
  LOCAL_R="$(sed -n 's/^R_VERSION=R version \([^ ]*\).*/\1/p' reproducibilidad/entorno_R_y_cmdstan.txt | head -n1)"
  SERVER_R="$(Rscript --vanilla -e 'cat(as.character(getRversion()))' 2>/dev/null || true)"
  if [[ -n "${LOCAL_R}" && "${LOCAL_R}" != "${SERVER_R}" ]]; then
    echo "ERROR: R local congelado=${LOCAL_R}; R servidor=${SERVER_R}"
    FAIL=1
  fi
fi

FREE_GB="$(df -Pk "${ROOT}" | awk 'NR==2 {print int($4/1024/1024)}')"
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

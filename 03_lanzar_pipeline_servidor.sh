#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT}"

CONFIG_FILE="${CONFIG_FILE:-${ROOT}/recursos_servidor.env}"
[[ -f "${CONFIG_FILE}" ]] && source "${CONFIG_FILE}"

MODE="${1:---background}"
RUN_BASE="${ROOT}/.server_runs"
CURRENT="${RUN_BASE}/current"
LOCK_FILE="${RUN_BASE}/pipeline.lock"
mkdir -p "${RUN_BASE}"

detect_cpu_limit () {
  local n quota period cg
  n="$(nproc)"
  if [[ -n "${SLURM_CPUS_PER_TASK:-}" && "${SLURM_CPUS_PER_TASK}" =~ ^[0-9]+$ ]]; then
    n="${SLURM_CPUS_PER_TASK}"
  elif [[ -r /sys/fs/cgroup/cpu.max ]]; then
    read -r quota period < /sys/fs/cgroup/cpu.max
    if [[ "${quota}" != "max" && "${quota}" =~ ^[0-9]+$ && "${period}" =~ ^[0-9]+$ ]]; then
      cg=$(( (quota + period - 1) / period ))
      (( cg < n )) && n="${cg}"
    fi
  fi
  (( n < 1 )) && n=1
  echo "${n}"
}

detect_memory_gb () {
  local kb bytes host_gb cg_gb
  kb="$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)"
  host_gb=$((kb / 1024 / 1024))
  if [[ -n "${SLURM_MEM_PER_NODE:-}" && "${SLURM_MEM_PER_NODE}" =~ ^[0-9]+$ ]]; then
    echo $((SLURM_MEM_PER_NODE / 1024))
    return
  fi
  if [[ -r /sys/fs/cgroup/memory.max ]]; then
    bytes="$(cat /sys/fs/cgroup/memory.max)"
    if [[ "${bytes}" != "max" && "${bytes}" =~ ^[0-9]+$ ]]; then
      cg_gb=$((bytes / 1024 / 1024 / 1024))
      if (( cg_gb > 0 && cg_gb < host_gb )); then
        echo "${cg_gb}"
        return
      fi
    fi
  fi
  echo "${host_gb}"
}

if [[ "${MODE}" == "--foreground" ]]; then
  RUN_ID="run_$(date -u +%Y%m%dT%H%M%SZ)"
  RUN_DIR="${RUN_BASE}/${RUN_ID}"
  mkdir -p "${RUN_DIR}"
  ln -sfn "${RUN_DIR}" "${CURRENT}"
  echo "$$" > "${RUN_DIR}/pipeline.pid"
  echo "$(ps -o pgid= -p $$ | tr -d ' ')" > "${RUN_DIR}/pipeline.pgid"
  ./01_preflight_servidor.sh > "${RUN_DIR}/preflight_lanzamiento.log" 2>&1
fi

if [[ "${MODE}" == "--worker" || "${MODE}" == "--foreground" ]]; then
  exec 9>"${LOCK_FILE}"
  if ! flock -n 9; then
    echo "ERROR: ya existe otra ejecución del pipeline."
    exit 1
  fi

  CPU_VISIBLE="$(detect_cpu_limit)"
  MEM_VISIBLE_GB="$(detect_memory_gb)"
  AUTO_CPU_PERCENT="${AUTO_CPU_PERCENT:-50}"
  AUTO_MAX_PIPELINE_R_CORES="${AUTO_MAX_PIPELINE_R_CORES:-2}"
  PIPELINE_R_CORES="${PIPELINE_R_CORES:-$((CPU_VISIBLE * AUTO_CPU_PERCENT / 100))}"
  (( PIPELINE_R_CORES < 1 )) && PIPELINE_R_CORES=1
  (( PIPELINE_R_CORES > AUTO_MAX_PIPELINE_R_CORES )) && PIPELINE_R_CORES="${AUTO_MAX_PIPELINE_R_CORES}"

  BRMS_INTERNAL_CORES="${BRMS_INTERNAL_CORES:-1}"
  BLAS_THREADS="${BLAS_THREADS:-1}"
  MIN_AVAILABLE_MEM_GB="${MIN_AVAILABLE_MEM_GB:-12}"
  MEMORY_WAIT_SECONDS="${MEMORY_WAIT_SECONDS:-60}"
  MEMORY_WAIT_ATTEMPTS="${MEMORY_WAIT_ATTEMPTS:-0}"
  MIN_FREE_DISK_GB="${MIN_FREE_DISK_GB:-100}"
  CMDSTAN_OPTIMIZATION_LEVEL="${CMDSTAN_OPTIMIZATION_LEVEL:-0}"
  ORCHESTRATOR="${ORCHESTRATOR:-run_comparacion_modelamiento2_vs_3_memoria_servidor.sh}"

  if (( MEM_VISIBLE_GB < MIN_AVAILABLE_MEM_GB + 2 )); then
    echo "ERROR: memoria visible=${MEM_VISIBLE_GB} GiB; el umbral de inicio es ${MIN_AVAILABLE_MEM_GB} GiB."
    echo "Reduce MIN_AVAILABLE_MEM_GB solo después de medir un escenario q00."
    exit 1
  fi

  FREE_GB="$(df -Pk "${ROOT}" | awk 'NR==2 {print int($4/1024/1024)}')"
  if (( FREE_GB < MIN_FREE_DISK_GB )); then
    echo "ERROR: espacio libre=${FREE_GB} GiB; mínimo=${MIN_FREE_DISK_GB} GiB."
    exit 1
  fi

  RUN_DIR="$(readlink -f "${CURRENT}")"
  {
    echo "utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "git=$(git describe --always --dirty --tags 2>/dev/null || true)"
    echo "cpu_visible=${CPU_VISIBLE}"
    echo "memory_visible_gb=${MEM_VISIBLE_GB}"
    echo "PIPELINE_R_CORES=${PIPELINE_R_CORES}"
    echo "BRMS_INTERNAL_CORES=${BRMS_INTERNAL_CORES}"
    echo "BLAS_THREADS=${BLAS_THREADS}"
    echo "MIN_AVAILABLE_MEM_GB=${MIN_AVAILABLE_MEM_GB}"
    echo "CMDSTAN_OPTIMIZATION_LEVEL=${CMDSTAN_OPTIMIZATION_LEVEL}"
  } > "${RUN_DIR}/recursos_efectivos.env"

  export PIPELINE_R_CORES BRMS_INTERNAL_CORES BLAS_THREADS
  export MIN_AVAILABLE_MEM_GB MEMORY_WAIT_SECONDS MEMORY_WAIT_ATTEMPTS
  export CMDSTAN_OPTIMIZATION_LEVEL
  export PAUSE_FILE="${ROOT}/PAUSAR_PIPELINE"

  rm -f "${PAUSE_FILE}"
  exec bash "${ORCHESTRATOR}"
fi

if [[ "${MODE}" != "--background" ]]; then
  echo "Uso: $0 [--background|--foreground]"
  exit 2
fi

RUN_ID="run_$(date -u +%Y%m%dT%H%M%SZ)"
RUN_DIR="${RUN_BASE}/${RUN_ID}"
mkdir -p "${RUN_DIR}"
ln -sfn "${RUN_DIR}" "${CURRENT}"

./01_preflight_servidor.sh > "${RUN_DIR}/preflight_lanzamiento.log" 2>&1

sha256sum renv.lock ./*.R ./*.sh plantillas/*.tsv \
  > "${RUN_DIR}/SHA256SUMS_lanzamiento.txt" 2>/dev/null || true
git status --short > "${RUN_DIR}/git_status.txt" 2>/dev/null || true
git rev-parse HEAD > "${RUN_DIR}/git_commit.txt" 2>/dev/null || true

RUNNER=("${ROOT}/03_lanzar_pipeline_servidor.sh" --worker)
NICE_LEVEL="${NICE_LEVEL:-10}"
if command -v ionice >/dev/null 2>&1; then
  LAUNCH=(ionice -c "${IONICE_CLASS:-2}" -n "${IONICE_LEVEL:-7}" nice -n "${NICE_LEVEL}")
else
  LAUNCH=(nice -n "${NICE_LEVEL}")
fi

nohup setsid "${LAUNCH[@]}" "${RUNNER[@]}" \
  > "${RUN_DIR}/pipeline_master.log" 2>&1 < /dev/null &
PID=$!
echo "${PID}" > "${RUN_DIR}/pipeline.pid"
echo "${PID}" > "${RUN_DIR}/pipeline.pgid"

sleep 2
if ! kill -0 "${PID}" 2>/dev/null; then
  echo "ERROR: el proceso terminó al arrancar."
  tail -n 80 "${RUN_DIR}/pipeline_master.log"
  exit 1
fi

echo "PIPELINE INICIADO"
echo "PID/PGID: ${PID}"
echo "RUN: ${RUN_DIR}"
echo "LOG: ${RUN_DIR}/pipeline_master.log"
echo "Estado: ./04_estado_pipeline.sh"
echo "Pausa segura: touch PAUSAR_PIPELINE"

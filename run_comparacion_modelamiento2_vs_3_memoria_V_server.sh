#!/usr/bin/env bash
set -euo pipefail

## ============================================================
## COMPARAR modelamiento2.R vs modelamiento3.R
##
## - Usa las 3 plantillas.
## - Ejecuta preprocess_data.R e INTEGRACION.R una sola vez por
##   escenario y copia ese mismo upstream a las dos variantes.
## - Guarda todo separado en:
##     ./comparacion_modelamientos/modelamiento2/
##     ./comparacion_modelamientos/modelamiento3/
## - Al final ejecuta FRAMEWORK_BIDIRECCIONAL.R en las 6 raíces.
## ============================================================

COMPARISON_ROOT="./comparacion_modelamientos"
SHARED_ROOT="${COMPARISON_ROOT}/upstream_comun"
STATUS_FILE="${COMPARISON_ROOT}/resumen_corridas.tsv"
# 
# TEMPLATES=(
#   "plantilla_1_todas_clinicas_barrido.tsv"
#   "plantilla_2_con_FLI_barrido.tsv"
#   "plantilla_3_sin_FLI_barrido.tsv"
# )
TEMPLATES=(
  "plantillas_v2/plantilla_1_todas_clinicas_barrido.tsv"
  "plantillas_v2/plantilla_2_con_FLI_barrido.tsv"
  "plantillas_v2/plantilla_3_sin_FLI_barrido.tsv"
)

## Las raíces se obtienen de OUTDIR_BASE en cada fila de la plantilla.
## Se guardan aquí para ejecutar después el framework bidireccional.
ROOT_NAMES=()
declare -A ROOT_NAMES_SEEN=()

MODEL_LABELS=("modelamiento2" "modelamiento3")
MODEL_SCRIPTS=("modelamiento2.R" "modelamiento3.R")

PREPROCESS_SCRIPT="preprocess_data.R"
INTEGRATION_SCRIPT="INTEGRACION.R"
PANEL_SCRIPT="panel_final_multimodelo.R"
DOWNSTREAM_SCRIPT="downstream.R"
FRAMEWORK_SCRIPT="FRAMEWORK_BIDIRECCIONAL.R"

FRAMEWORK_DIAG="diagnostico_bayes_features"
FRAMEWORK_MODE="both"

## ------------------------------------------------------------
## Control de memoria y paralelismo
##
## No modifica la plantilla ni los parámetros estadísticos.
## Solo limita la cantidad de procesos/hilos que pueden coexistir.
##
## Con PIPELINE_R_CORES=2:
## - modelamiento2.R/modelamiento3.R calculan cores_cv=1;
## - panel_final_multimodelo.R calcula N_CORES=1.
##
## Los valores se pueden cambiar al lanzar el script, por ejemplo:
##   PIPELINE_R_CORES=3 BLAS_THREADS=1 bash este_script.sh
## ------------------------------------------------------------
PIPELINE_R_CORES="${PIPELINE_R_CORES:-1}"
BRMS_INTERNAL_CORES="${BRMS_INTERNAL_CORES:-1}"
BLAS_THREADS="${BLAS_THREADS:-1}"
MIN_AVAILABLE_MEM_GB="${MIN_AVAILABLE_MEM_GB:-12}"
MEMORY_WAIT_SECONDS="${MEMORY_WAIT_SECONDS:-15}"
MEMORY_WAIT_ATTEMPTS="${MEMORY_WAIT_ATTEMPTS:-8}"

CMDSTAN_OPTIMIZATION_LEVEL="${CMDSTAN_OPTIMIZATION_LEVEL:-0}"

for VALUE_NAME in   PIPELINE_R_CORES   BRMS_INTERNAL_CORES   BLAS_THREADS   MIN_AVAILABLE_MEM_GB   MEMORY_WAIT_SECONDS   MEMORY_WAIT_ATTEMPTS
do
  VALUE="${!VALUE_NAME}"
  if [[ ! "${VALUE}" =~ ^[0-9]+$ ]] || (( VALUE < 1 )); then
    echo "ERROR: ${VALUE_NAME} debe ser un entero >= 1; recibido: ${VALUE}"
    exit 1
  fi
done

if [[ ! "${CMDSTAN_OPTIMIZATION_LEVEL}" =~ ^[0-3]$ ]]; then
  echo "ERROR: CMDSTAN_OPTIMIZATION_LEVEL debe estar entre 0 y 3."
  echo "Valor recibido: ${CMDSTAN_OPTIMIZATION_LEVEL}"
  exit 1
fi
export PIPELINE_R_CORES
export BRMS_INTERNAL_CORES
export CMDSTAN_OPTIMIZATION_LEVEL

## CmdStan lee O y genera -O0, -O1, -O2 o -O3.
export O="${CMDSTAN_OPTIMIZATION_LEVEL}"

export OMP_NUM_THREADS="${BLAS_THREADS}"
export OPENBLAS_NUM_THREADS="${BLAS_THREADS}"
export MKL_NUM_THREADS="${BLAS_THREADS}"
export BLIS_NUM_THREADS="${BLAS_THREADS}"
export VECLIB_MAXIMUM_THREADS="${BLAS_THREADS}"
export NUMEXPR_NUM_THREADS="${BLAS_THREADS}"
export STAN_NUM_THREADS="${BLAS_THREADS}"
export TBB_NUM_THREADS="${BLAS_THREADS}"
export MALLOC_ARENA_MAX="${MALLOC_ARENA_MAX:-2}"

## La compilación de CmdStan siempre será serial.
export MAKEFLAGS="-j1"

export MC_CORES="${BRMS_INTERNAL_CORES}"

R_RUNTIME_DIR="${COMPARISON_ROOT}/.runtime_memoria"
R_PROFILE_LIMITED="${R_RUNTIME_DIR}/Rprofile_pipeline_memoria.R"
R_MAKEVARS_LIMITED="${R_RUNTIME_DIR}/Makevars_pipeline_memoria"

## ------------------------------------------------------------
## Validación
## ------------------------------------------------------------

REQUIRED_FILES=(
  "${TEMPLATES[@]}"
  "${PREPROCESS_SCRIPT}"
  "${INTEGRATION_SCRIPT}"
  "${MODEL_SCRIPTS[@]}"
  "${PANEL_SCRIPT}"
  "${DOWNSTREAM_SCRIPT}"
  "${FRAMEWORK_SCRIPT}"
)

for FILE in "${REQUIRED_FILES[@]}"; do
  if [[ ! -f "${FILE}" ]]; then
    echo "ERROR: no encuentro ${FILE}"
    exit 1
  fi
done

## Evita mezclar una corrida nueva con resultados anteriores.
## Permitir reanudar una ejecución previa.
mkdir -p \
  "${SHARED_ROOT}" \
  "${COMPARISON_ROOT}/logs_bidir" \
  "${R_RUNTIME_DIR}"

## Perfil temporal de R usado solamente durante esta comparación.
## Conserva ~/.Rprofile si existe, limita parallel::detectCores()
## y deja los procesos internos de BRMS en un solo core.
cat > "${R_PROFILE_LIMITED}" <<'RS_PROFILE'
.user_profile <- path.expand("~/.Rprofile")
if (file.exists(.user_profile)) {
  try(sys.source(.user_profile, envir = .GlobalEnv), silent = TRUE)
}

.pipeline_cores <- suppressWarnings(
  as.integer(Sys.getenv("PIPELINE_R_CORES", unset = "2"))
)
.brms_internal_cores <- suppressWarnings(
  as.integer(Sys.getenv("BRMS_INTERNAL_CORES", unset = "1"))
)

if (is.na(.pipeline_cores) || .pipeline_cores < 1L) {
  .pipeline_cores <- 2L
}
if (is.na(.brms_internal_cores) || .brms_internal_cores < 1L) {
  .brms_internal_cores <- 1L
}

options(mc.cores = .brms_internal_cores)

## Los scripts calculan sus workers a partir de parallel::detectCores().
## Se limita únicamente durante esta ejecución, sin editar los scripts R.
.parallel_ns <- asNamespace("parallel")
try(unlockBinding("detectCores", .parallel_ns), silent = TRUE)
assign(
  "detectCores",
  local({
    n <- .pipeline_cores
    function(all.tests = FALSE, logical = TRUE) n
  }),
  envir = .parallel_ns
)
try(lockBinding("detectCores", .parallel_ns), silent = TRUE)

rm(
  .user_profile,
  .pipeline_cores,
  .brms_internal_cores,
  .parallel_ns
)
RS_PROFILE

export R_PROFILE_USER="${R_PROFILE_LIMITED}"

cat > "${R_MAKEVARS_LIMITED}" <<'MAKEVARS'
MAKEFLAGS = -j1

CXXFLAGS = -O0 -g0
CXX11FLAGS = -O0 -g0
CXX14FLAGS = -O0 -g0
CXX17FLAGS = -O0 -g0
MAKEVARS

export R_MAKEVARS_USER="${R_MAKEVARS_LIMITED}"

## ------------------------------------------------------------
## Verificación real de la configuración de CmdStan
## ------------------------------------------------------------

CMDSTAN_PATH="$(
  Rscript -e 'cat(cmdstanr::cmdstan_path())'
)"

if [[ ! -d "${CMDSTAN_PATH}" ]]; then
  echo "ERROR: no existe la instalación de CmdStan:"
  echo "${CMDSTAN_PATH}"
  exit 1
fi

CMDSTAN_O_ACTIVE="$(
  make \
    --no-print-directory \
    -s \
    -C "${CMDSTAN_PATH}" \
    print-compiler-flags |
    awk '
      /O \(Optimization Level\)/ && !found {
        value = $NF
        found = 1
      }
      END {
        if (found) print value
      }
    '
)"
if [[ "${CMDSTAN_O_ACTIVE}" != "${CMDSTAN_OPTIMIZATION_LEVEL}" ]]; then
  echo "ERROR: CmdStan no recibió el nivel de optimización esperado."
  echo "Esperado: ${CMDSTAN_OPTIMIZATION_LEVEL}"
  echo "Detectado: ${CMDSTAN_O_ACTIVE:-NO_DETECTADO}"
  echo "Ruta: ${CMDSTAN_PATH}"
  exit 1
fi

export CMDSTAN_PATH

echo "CmdStan path          : ${CMDSTAN_PATH}"
echo "CmdStan optimization  : -O${CMDSTAN_O_ACTIVE}"
## No sobrescribir el resumen anterior al reanudar.
if [[ ! -f "${STATUS_FILE}" ]]; then
  printf \
    'modelamiento\tplantilla\traiz\tanalisis\tQ_TX\tQ_PR\tQ_ME\tQ_CL\truta\n' \
    > "${STATUS_FILE}"
fi

## ------------------------------------------------------------
## Funciones
## ------------------------------------------------------------

available_memory_kb () {
  awk '/^MemAvailable:/ {print $2; exit}' /proc/meminfo
}

wait_for_memory () {
  local REQUIRED_KB=$((MIN_AVAILABLE_MEM_GB * 1024 * 1024))
  local ATTEMPT
  local AVAILABLE_KB
  local AVAILABLE_GB

  for ((ATTEMPT = 1; ATTEMPT <= MEMORY_WAIT_ATTEMPTS; ATTEMPT++)); do
    AVAILABLE_KB="$(available_memory_kb)"
    AVAILABLE_KB="${AVAILABLE_KB:-0}"
    AVAILABLE_GB="$(awk -v kb="${AVAILABLE_KB}" 'BEGIN {printf "%.2f", kb/1024/1024}')"

    if (( AVAILABLE_KB >= REQUIRED_KB )); then
      echo "RAM disponible antes del paso: ${AVAILABLE_GB} GiB"
      return 0
    fi

    echo "AVISO: solo hay ${AVAILABLE_GB} GiB disponibles; se requieren ${MIN_AVAILABLE_MEM_GB} GiB."
    echo "Esperando ${MEMORY_WAIT_SECONDS} s para que el sistema libere memoria (${ATTEMPT}/${MEMORY_WAIT_ATTEMPTS})..."
    sleep "${MEMORY_WAIT_SECONDS}"
  done

  echo "ERROR: no se alcanzaron ${MIN_AVAILABLE_MEM_GB} GiB de RAM disponible."
  echo "Cierra otros procesos o aumenta swap antes de reanudar."
  return 1
}

remember_root_name () {
  local ROOT="$1"

  if [[ -z "${ROOT_NAMES_SEEN[${ROOT}]+x}" ]]; then
    ROOT_NAMES+=("${ROOT}")
    ROOT_NAMES_SEEN["${ROOT}"]=1
  fi
}

resolve_root_name () {
  local TEMPLATE_PATH="$1"
  local REQUESTED_ROOT
  local BASE
  local ENTRY
  local -a MATCHES=()
  local MATCH

  TEMPLATE_PATH="${TEMPLATE_PATH%/}"
  REQUESTED_ROOT="$(basename "${TEMPLATE_PATH}")"

  if [[ -z "${REQUESTED_ROOT}" || "${REQUESTED_ROOT}" == "." || "${REQUESTED_ROOT}" == "/" ]]; then
    echo "ERROR: OUTDIR_BASE inválido en la plantilla: ${TEMPLATE_PATH}" >&2
    return 1
  fi

  ## Si ya existe una carpeta con el mismo nombre y capitalización, usarla.
  for BASE in     "${COMPARISON_ROOT}/modelamiento2"     "${COMPARISON_ROOT}/modelamiento3"     "${SHARED_ROOT}"
  do
    if [[ -d "${BASE}/${REQUESTED_ROOT}" ]]; then
      printf '%s\n' "${REQUESTED_ROOT}"
      return 0
    fi
  done

  ## Evitar duplicar una raíz que ya existe con otra capitalización,
  ## por ejemplo escenarios_sin_FLI frente a escenarios_SIN_FLI.
  for BASE in     "${COMPARISON_ROOT}/modelamiento2"     "${COMPARISON_ROOT}/modelamiento3"     "${SHARED_ROOT}"
  do
    [[ -d "${BASE}" ]] || continue

    while IFS= read -r ENTRY; do
      [[ -n "${ENTRY}" ]] || continue
      if [[ "${ENTRY,,}" == "${REQUESTED_ROOT,,}" ]]; then
        MATCHES+=("${ENTRY}")
      fi
    done < <(find "${BASE}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n')
  done

  if ((${#MATCHES[@]} > 0)); then
    MATCH="$(printf '%s\n' "${MATCHES[@]}" | sort -u | head -n 1)"
    echo "AVISO: la plantilla indica '${REQUESTED_ROOT}', pero la estructura existente usa '${MATCH}'." >&2
    echo "Se reutilizará '${MATCH}' sin modificar la plantilla." >&2
    printf '%s\n' "${MATCH}"
    return 0
  fi

  printf '%s\n' "${REQUESTED_ROOT}"
}

run_step () {
  local STEP_NAME="$1"
  local SCRIPT_NAME="$2"
  local LOG_FILE="${LOG_DIR}/${STEP_NAME}.log"
  local R_PORT

  mkdir -p "${LOG_DIR}"

  echo
  echo "============================================================"
  echo "RUNNING : ${SCRIPT_NAME}"
  echo "VARIANTE: ${MODEL_VARIANT}"
  echo "ANALISIS: ${ANALYSIS_NAME}"
  echo "LOG     : ${LOG_FILE}"
  echo "============================================================"

  ## Obtener un puerto local libre para los clústeres PSOCK de R.
  R_PORT="$(
    python3 - <<'PY'
import socket

with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
    sock.bind(("127.0.0.1", 0))
    print(sock.getsockname()[1])
PY
  )"

  if [[ -z "${R_PORT}" ]]; then
    echo "ERROR: no se pudo obtener un puerto libre para R."
    exit 1
  fi

  export R_PARALLEL_PORT="${R_PORT}"

  echo "R_PARALLEL_PORT    : ${R_PARALLEL_PORT}"
  echo "PIPELINE_R_CORES   : ${PIPELINE_R_CORES}"
  echo "BRMS_INTERNAL_CORES: ${BRMS_INTERNAL_CORES}"
  echo "BLAS_THREADS       : ${BLAS_THREADS}"
  echo "CMDSTAN_PATH       : ${CMDSTAN_PATH}"
  echo "CMDSTAN_OPTIMIZATION: -O${CMDSTAN_O_ACTIVE}"
  echo "MAKEFLAGS          : ${MAKEFLAGS}"

  wait_for_memory

  Rscript "${SCRIPT_NAME}" 2>&1 | tee "${LOG_FILE}"
}


export_config () {
  export ANALYSIS_NAME
  export OUTDIR_BASE
  export Q_TX
  export Q_PR
  export Q_ME
  export Q_CL
  export PANEL_CHOICE
  export CLINICAL_INCLUDE
  export CLINICAL_EXCLUDE
  export MODEL_SOURCE
  export MODEL_VARIANT
  export DIAG_TARGET
}

write_environment () {
  local FILE="$1"

  {
    echo "ANALYSIS_NAME=${ANALYSIS_NAME}"
    echo "OUTDIR_BASE=${OUTDIR_BASE}"
    echo "ANALYSIS_DIR=${ANALYSIS_DIR}"
    echo "Q_TX=${Q_TX}"
    echo "Q_PR=${Q_PR}"
    echo "Q_ME=${Q_ME}"
    echo "Q_CL=${Q_CL}"
    echo "PANEL_CHOICE=${PANEL_CHOICE}"
    echo "CLINICAL_INCLUDE=${CLINICAL_INCLUDE}"
    echo "CLINICAL_EXCLUDE=${CLINICAL_EXCLUDE}"
    echo "MODEL_SOURCE=${MODEL_SOURCE}"
    echo "MODEL_VARIANT=${MODEL_VARIANT}"
    echo "DIAG_TARGET=${DIAG_TARGET}"
    echo "TEMPLATE_FILE=${TEMPLATE_FILE}"
    echo "PIPELINE_R_CORES=${PIPELINE_R_CORES}"
    echo "BRMS_INTERNAL_CORES=${BRMS_INTERNAL_CORES}"
    echo "BLAS_THREADS=${BLAS_THREADS}"
    echo "CMDSTAN_PATH=${CMDSTAN_PATH}"
    echo "CMDSTAN_OPTIMIZATION_LEVEL=${CMDSTAN_OPTIMIZATION_LEVEL}"
    echo "CMDSTAN_O_ACTIVE=${CMDSTAN_O_ACTIVE}"
    echo "MAKEFLAGS=${MAKEFLAGS}"
    echo "MIN_AVAILABLE_MEM_GB=${MIN_AVAILABLE_MEM_GB}"
    echo "R_PROFILE_USER=${R_PROFILE_USER}"
    echo "R_MAKEVARS_USER=${R_MAKEVARS_USER}"
  } > "${FILE}"
}

normalize_model_dir () {
  ## panel_final_multimodelo.R y downstream.R trabajan con
  ## MODEL_SOURCE=NEW y normalmente esperan:
  ##
  ##   ANALYSIS_DIR/modelamiento/
  ##
  ## Si modelamiento2.R o modelamiento3.R guarda en una carpeta
  ## con su propio nombre, se crea un enlace simbólico compatible.

  if [[ -d "${ANALYSIS_DIR}/modelamiento" ]]; then
    return 0
  fi

  if [[ -d "${ANALYSIS_DIR}/${MODEL_VARIANT}" ]]; then
    ln -s \
      "${MODEL_VARIANT}" \
      "${ANALYSIS_DIR}/modelamiento"

    return 0
  fi

  echo "ERROR: ${MODEL_VARIANT} no generó:"
  echo "  ${ANALYSIS_DIR}/modelamiento"
  echo "ni"
  echo "  ${ANALYSIS_DIR}/${MODEL_VARIANT}"
  exit 1
}

run_downstream () {
  local TARGET="$1"
  local ORIGINAL_TARGET="${DIAG_TARGET}"
  local DONE_FILE

  case "${TARGET}" in
    BAYES_INPUT)
      DONE_FILE="${ANALYSIS_DIR}/diagnostico_bayes_features/rds/diagnostic_objects.rds"
      ;;

    PANEL_FINAL)
      DONE_FILE="${ANALYSIS_DIR}/diagnostico_panel_final/rds/diagnostic_objects.rds"
      ;;

    *)
      echo "ERROR: TARGET downstream no reconocido: ${TARGET}"
      exit 1
      ;;
  esac

  if [[ -s "${DONE_FILE}" ]]; then
    echo
    echo "SKIP downstream ${TARGET}: ya está completo."
    echo "Archivo detectado: ${DONE_FILE}"
    return 0
  fi

  DIAG_TARGET="${TARGET}"
  export DIAG_TARGET

  run_step \
    "05_downstream_${TARGET}" \
    "${DOWNSTREAM_SCRIPT}"

  DIAG_TARGET="${ORIGINAL_TARGET}"
  export DIAG_TARGET
}
## ------------------------------------------------------------
## Plantilla(s) × escenarios definidos en cada TSV × 2 modelamientos
## ------------------------------------------------------------

RUN_ID=0

for TEMPLATE_INDEX in "${!TEMPLATES[@]}"; do

  TEMPLATE_FILE="${TEMPLATES[${TEMPLATE_INDEX}]}"

while IFS= read -r TSV_LINE || [[ -n "${TSV_LINE}" ]]
do

  ## Quitar retorno de carro si el TSV tiene finales de línea Windows.
  TSV_LINE="${TSV_LINE%$'\r'}"

  ## Convertir TAB a un separador que Bash no considere espacio.
  ## Así se conservan correctamente los campos TSV vacíos.
  IFS=$'\x1f' read -r \
    ANALYSIS_NAME \
    TEMPLATE_OUTDIR_BASE \
    Q_TX \
    Q_PR \
    Q_ME \
    Q_CL \
    PANEL_CHOICE \
    CLINICAL_INCLUDE \
    CLINICAL_EXCLUDE \
    TEMPLATE_MODEL_SOURCE \
    DIAG_TARGET \
    EXTRA_COLS \
    <<< "${TSV_LINE//$'\t'/$'\x1f'}"
    ## Saltar líneas vacías.
    [[ -z "${ANALYSIS_NAME// }" ]] && continue

    ## Saltar comentarios.
    [[ "${ANALYSIS_NAME}" =~ ^# ]] && continue

    ## Saltar cabecera.
    [[ "${ANALYSIS_NAME}" == "ANALYSIS_NAME" ]] && continue

    if [[ -z "${TEMPLATE_OUTDIR_BASE}" ]]; then
      echo "ERROR: OUTDIR_BASE vacío en ${ANALYSIS_NAME}"
      exit 1
    fi

    ## Cada fila decide su familia: escenarios, con FLI o sin FLI.
    ROOT_NAME="$(resolve_root_name "${TEMPLATE_OUTDIR_BASE}")"
    remember_root_name "${ROOT_NAME}"

    RUN_ID=$((RUN_ID + 1))

    ## Defaults.
    Q_TX="${Q_TX:-0.97}"
    Q_PR="${Q_PR:-0.40}"
    Q_ME="${Q_ME:-0.97}"
    Q_CL="${Q_CL:-0.70}"

    PANEL_CHOICE="${PANEL_CHOICE:-FINAL}"

    CLINICAL_INCLUDE="${CLINICAL_INCLUDE:-}"
    CLINICAL_EXCLUDE="${CLINICAL_EXCLUDE:-}"

    DIAG_TARGET="${DIAG_TARGET:-BOTH}"

    PANEL_CHOICE="${PANEL_CHOICE^^}"
    DIAG_TARGET="${DIAG_TARGET^^}"

    if [[ -z "${CLINICAL_INCLUDE}" ]]; then
      echo "ERROR: CLINICAL_INCLUDE vacío en ${ANALYSIS_NAME}"
      exit 1
    fi

    case "${DIAG_TARGET}" in
      BAYES_INPUT|PANEL_FINAL|BOTH)
        ;;
      *)
        echo "ERROR: DIAG_TARGET inválido: ${DIAG_TARGET}"
        echo "ANALYSIS_NAME: ${ANALYSIS_NAME}"
        exit 1
        ;;
    esac

    echo
    echo "################################################################"
    echo "ESCENARIO ${RUN_ID}"
    echo "PLANTILLA: ${TEMPLATE_FILE}"
    echo "RAÍZ     : ${ROOT_NAME}"
    echo "ANÁLISIS : ${ANALYSIS_NAME}"
    echo "Q_TX     : ${Q_TX}"
    echo "Q_PR     : ${Q_PR}"
    echo "Q_ME     : ${Q_ME}"
    echo "Q_CL     : ${Q_CL}"
    echo "################################################################"

    ## ==========================================================
    ## A. Upstream común
    ##
    ## El preprocesamiento y MOFA se ejecutan una sola vez.
    ## Después se copia exactamente el mismo resultado a las dos
    ## ramas de modelamiento.
    ## ==========================================================

    OUTDIR_BASE="${SHARED_ROOT}/${ROOT_NAME}"
    ANALYSIS_DIR="${OUTDIR_BASE}/${ANALYSIS_NAME}"
    LOG_DIR="${ANALYSIS_DIR}/logs"
    
    MODEL_SOURCE="${TEMPLATE_MODEL_SOURCE:-NEW}"
    MODEL_VARIANT="upstream_comun"
    
    if [[ -f "${ANALYSIS_DIR}/.UPSTREAM_COMPLETE" ]]; then
    
      echo
      echo "SKIP upstream: preprocesamiento e integración ya completos."
      echo "ANÁLISIS: ${ANALYSIS_NAME}"
    
    else
    
      mkdir -p "${LOG_DIR}"
    
      export_config
    
      write_environment \
        "${LOG_DIR}/run_environment_upstream.txt"
    
      echo "${TEMPLATE_FILE}" \
        > "${LOG_DIR}/template_file_used.txt"
    
      run_step \
        "01_preprocess_data" \
        "${PREPROCESS_SCRIPT}"
    
      run_step \
        "02_integracion_mofa" \
        "${INTEGRATION_SCRIPT}"
    
      touch "${ANALYSIS_DIR}/.UPSTREAM_COMPLETE"
    
    fi
    
    SHARED_ANALYSIS_DIR="${ANALYSIS_DIR}"

    ## ==========================================================
    ## B. Ejecutar modelamiento2 y modelamiento3
    ## ==========================================================

    for MODEL_INDEX in "${!MODEL_LABELS[@]}"; do

      MODEL_VARIANT="${MODEL_LABELS[${MODEL_INDEX}]}"
      MODEL_SCRIPT="${MODEL_SCRIPTS[${MODEL_INDEX}]}"

      OUTDIR_BASE="${COMPARISON_ROOT}/${MODEL_VARIANT}/${ROOT_NAME}"
      ANALYSIS_DIR="${OUTDIR_BASE}/${ANALYSIS_NAME}"
      LOG_DIR="${ANALYSIS_DIR}/logs"

      ## Respetar MODEL_SOURCE indicado en la plantilla.
      MODEL_SOURCE="${TEMPLATE_MODEL_SOURCE:-NEW}"

      mkdir -p "${OUTDIR_BASE}"
      
      ## Si toda esta variante ya terminó, no repetir nada.
      if [[ -f "${ANALYSIS_DIR}/.PIPELINE_COMPLETE" ]]; then
        echo
        echo "SKIP variante completa:"
        echo "MODELO   : ${MODEL_VARIANT}"
        echo "ANÁLISIS : ${ANALYSIS_NAME}"
        continue
      fi
      
      ## Copiar el upstream únicamente cuando la rama todavía no existe.
      ## No volver a copiar sobre una rama parcial porque podría crear
      ## directorios anidados o sobrescribir resultados.
      if [[ ! -d "${ANALYSIS_DIR}" ]]; then
      
        cp -a \
          "${SHARED_ANALYSIS_DIR}" \
          "${ANALYSIS_DIR}"
      
      else
      
        echo
        echo "REANUDANDO rama existente:"
        echo "${ANALYSIS_DIR}"
      
      fi
      
      LOG_DIR="${ANALYSIS_DIR}/logs"
      mkdir -p "${LOG_DIR}"
      
      export_config
      
      write_environment \
        "${LOG_DIR}/run_environment_${MODEL_VARIANT}.txt"
      
      echo "${SHARED_ANALYSIS_DIR}" \
        > "${LOG_DIR}/upstream_source.txt"
      
      ## --------------------------------------------------------
      ## Modelamiento específico
      ## --------------------------------------------------------
      
      MODEL_METRICS_FILE="${ANALYSIS_DIR}/modelamiento/tables/estadisticas_modelo_global.csv"
      MODEL_FEATURES_FILE="${ANALYSIS_DIR}/modelamiento/rds/features_finales_integradas_mofa_brms.rds"
      
      if [[ -s "${MODEL_METRICS_FILE}" &&
            -s "${MODEL_FEATURES_FILE}" ]]; then
      
        echo
        echo "SKIP ${MODEL_VARIANT}: modelamiento ya completo."
        echo "Métricas: ${MODEL_METRICS_FILE}"
        echo "Features: ${MODEL_FEATURES_FILE}"
      
      else
      
        run_step \
          "03_${MODEL_VARIANT}" \
          "${MODEL_SCRIPT}"
      
      fi
      
      normalize_model_dir
      
      ## --------------------------------------------------------
      ## Panel final multimodelo
      ## --------------------------------------------------------
      
      PANEL_OBJECT_FILE="${ANALYSIS_DIR}/panel_final_multimodelo/rds/panel_final_multimodel_object.rds"
      PANEL_BRMS_FILE="${ANALYSIS_DIR}/panel_final_multimodelo/rds/brms_confirmatory_object.rds"
      
      if [[ -s "${PANEL_OBJECT_FILE}" &&
            -s "${PANEL_BRMS_FILE}" ]]; then
      
        echo
        echo "SKIP panel final: ya está completo."
        echo "Panel: ${PANEL_OBJECT_FILE}"
      
      else
      
        run_step \
          "04_panel_multimodelo_${MODEL_VARIANT}" \
          "${PANEL_SCRIPT}"
      
      fi
      
      ## --------------------------------------------------------
      ## Downstream
      ## --------------------------------------------------------
      
      case "${DIAG_TARGET}" in
      
        BAYES_INPUT)
          run_downstream "BAYES_INPUT"
          ;;
      
        PANEL_FINAL)
          run_downstream "PANEL_FINAL"
          ;;
      
        BOTH)
          run_downstream "BAYES_INPUT"
          run_downstream "PANEL_FINAL"
          ;;
      
      esac
      
      touch "${ANALYSIS_DIR}/.PIPELINE_COMPLETE"

      ## Registrar la corrida completa.
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${MODEL_VARIANT}" \
        "${TEMPLATE_FILE}" \
        "${ROOT_NAME}" \
        "${ANALYSIS_NAME}" \
        "${Q_TX}" \
        "${Q_PR}" \
        "${Q_ME}" \
        "${Q_CL}" \
        "${ANALYSIS_DIR}" \
        >> "${STATUS_FILE}"

      echo
      echo "============================================================"
      echo "CORRIDA TERMINADA"
      echo "MODELO    : ${MODEL_VARIANT}"
      echo "ANÁLISIS  : ${ANALYSIS_NAME}"
      echo "RESULTADOS: ${ANALYSIS_DIR}"
      echo "============================================================"

    done

  done < "${TEMPLATE_FILE}"

done

## ------------------------------------------------------------
## Ejecutar FRAMEWORK_BIDIRECCIONAL.R en las seis raíces
## ------------------------------------------------------------

for MODEL_VARIANT in "${MODEL_LABELS[@]}"; do

  for ROOT_NAME in "${ROOT_NAMES[@]}"; do

    ROOT_DIR="${COMPARISON_ROOT}/${MODEL_VARIANT}/${ROOT_NAME}"

    BIDIR_LOG_DIR="${COMPARISON_ROOT}/logs_bidir/${MODEL_VARIANT}"
    BIDIR_LOG="${BIDIR_LOG_DIR}/bidir_${ROOT_NAME}.log"

    mkdir -p "${BIDIR_LOG_DIR}"

    echo
    echo "============================================================"
    echo "FRAMEWORK BIDIRECCIONAL"
    echo "MODELO: ${MODEL_VARIANT}"
    echo "ROOT  : ${ROOT_DIR}"
    echo "LOG   : ${BIDIR_LOG}"
    echo "============================================================"

    wait_for_memory

    Rscript "${FRAMEWORK_SCRIPT}" \
      "${ROOT_DIR}" \
      "${FRAMEWORK_DIAG}" \
      "${FRAMEWORK_MODE}" \
      2>&1 | tee "${BIDIR_LOG}"

  done

done

## ------------------------------------------------------------
## Crear índices para localizar las estadísticas comparables
## ------------------------------------------------------------

find \
  "${COMPARISON_ROOT}/modelamiento2" \
  "${COMPARISON_ROOT}/modelamiento3" \
  -type f \
  -name 'estadisticas_modelo_global.csv' \
  -print \
  | sort \
  > "${COMPARISON_ROOT}/archivos_estadisticas_modelo_global.txt"

find \
  "${COMPARISON_ROOT}/modelamiento2" \
  "${COMPARISON_ROOT}/modelamiento3" \
  -type f \
  -name 'estadisticas_modelo_por_clase.csv' \
  -print \
  | sort \
  > "${COMPARISON_ROOT}/archivos_estadisticas_modelo_por_clase.txt"

echo
echo "============================================================"
echo "COMPARACIÓN TERMINADA"
echo "============================================================"
echo "Control de memoria:"
echo "  PIPELINE_R_CORES    = ${PIPELINE_R_CORES}"
echo "  BRMS_INTERNAL_CORES = ${BRMS_INTERNAL_CORES}"
echo "  BLAS_THREADS        = ${BLAS_THREADS}"
echo
echo "Escenarios upstream ejecutados : ${RUN_ID}"
echo "Corridas de modelamiento       : $((RUN_ID * 2))"
echo
echo "Resultados modelamiento2:"
echo "  ${COMPARISON_ROOT}/modelamiento2"
echo
echo "Resultados modelamiento3:"
echo "  ${COMPARISON_ROOT}/modelamiento3"
echo
echo "Upstream común:"
echo "  ${COMPARISON_ROOT}/upstream_comun"
echo
echo "Resumen de corridas:"
echo "  ${STATUS_FILE}"
echo
echo "Logs del framework:"
echo "  ${COMPARISON_ROOT}/logs_bidir"
echo
echo "Índice de estadísticas globales:"
echo "  ${COMPARISON_ROOT}/archivos_estadisticas_modelo_global.txt"
echo
echo "Índice de estadísticas por clase:"
echo "  ${COMPARISON_ROOT}/archivos_estadisticas_modelo_por_clase.txt"
echo "============================================================"
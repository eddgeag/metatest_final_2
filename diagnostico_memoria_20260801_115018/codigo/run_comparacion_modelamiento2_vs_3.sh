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
    "plantilla_base.tsv"
)

ROOT_NAMES=(
  "escenarios"
  "escenarios_con_FLI"
  "escenarios_SIN_FLI"
)

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
  "${COMPARISON_ROOT}/logs_bidir"

## No sobrescribir el resumen anterior al reanudar.
if [[ ! -f "${STATUS_FILE}" ]]; then
  printf \
    'modelamiento\tplantilla\traiz\tanalisis\tQ_TX\tQ_PR\tQ_ME\tQ_CL\truta\n' \
    > "${STATUS_FILE}"
fi

## ------------------------------------------------------------
## Funciones
## ------------------------------------------------------------

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

  echo "R_PARALLEL_PORT: ${R_PARALLEL_PORT}"

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
## 3 plantillas × 66 escenarios × 2 modelamientos
## ------------------------------------------------------------

RUN_ID=0

for TEMPLATE_INDEX in "${!TEMPLATES[@]}"; do

  TEMPLATE_FILE="${TEMPLATES[${TEMPLATE_INDEX}]}"
  ROOT_NAME="${ROOT_NAMES[${TEMPLATE_INDEX}]}"

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
    
    MODEL_SOURCE="NEW"
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

      ## Los scripts panel_final_multimodelo.R y downstream.R
      ## deben leer la estructura correspondiente a NEW.
      MODEL_SOURCE="NEW"

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
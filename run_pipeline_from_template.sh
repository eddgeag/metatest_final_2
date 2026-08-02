#!/usr/bin/env bash
set -euo pipefail

## ============================================================
## MASTER PIPELINE MULTIÓMICO PCOS DESDE PLANTILLA TXT/TSV
##
## Uso:
##   chmod +x run_pipeline_from_template.sh
##   ./run_pipeline_from_template.sh plantilla_pipeline.tsv
##
## Plantilla TAB-separada con columnas mínimas:
## ANALYSIS_NAME  OUTDIR_BASE  Q_TX  Q_PR  Q_ME  Q_CL  PANEL_CHOICE  CLINICAL_INCLUDE  CLINICAL_EXCLUDE
##
## Columnas opcionales:
## MODEL_SOURCE  DIAG_TARGET
##
## MODEL_SOURCE:
##   NEW  = usa modelamiento/
##   OLD  = usa modelado_anterior/
##   BOTH = ejecuta ambos desde cero: modelado_anterior/ y modelamiento/
##
## DIAG_TARGET:
##   BAYES_INPUT  = diagnóstico de features finales BRMS/MOFA
##   PANEL_FINAL  = diagnóstico del panel final multimodelo
##   BOTH         = ambos diagnósticos
##
## Ejecuta por cada fila:
## 1) preprocess_data.R
## 2) INTEGRACION.R
## 3) modelamiento_anterior.R si MODEL_SOURCE=OLD/BOTH
## 4) modelamiento2.R si MODEL_SOURCE=NEW/BOTH
## 5) panel_final_multimodelo.R para cada fuente disponible
## 6) downstream.R para diagnósticos válidos
## ============================================================

TEMPLATE_FILE="${1:-}"

if [[ -z "${TEMPLATE_FILE}" ]]; then
  echo "ERROR: Debes pasar una plantilla como argumento."
  echo "Uso:"
  echo "  ./run_pipeline_from_template.sh plantilla_pipeline.tsv"
  exit 1
fi

if [[ ! -f "${TEMPLATE_FILE}" ]]; then
  echo "ERROR: No existe la plantilla: ${TEMPLATE_FILE}"
  exit 1
fi

## -----------------------------
## Scripts requeridos
## -----------------------------

REQUIRED_SCRIPTS=(
  "preprocess_data.R"
  "INTEGRACION.R"
  "modelamiento_anterior.R"
  "modelamiento3.R"
  "panel_final_multimodelo.R"
  "downstream.R"
)

for s in "${REQUIRED_SCRIPTS[@]}"; do
  if [[ ! -f "${s}" ]]; then
    echo "ERROR: No encuentro el script requerido: ${s}"
    exit 1
  fi
done

## -----------------------------
## Función para ejecutar pasos
## -----------------------------

run_step () {
  local step_name="$1"
  local script_name="$2"
  local log_file="${LOG_DIR}/${step_name}.log"

  echo
  echo "============================================================"
  echo "RUNNING: ${script_name}"
  echo "LOG    : ${log_file}"
  echo "============================================================"

  Rscript "${script_name}" 2>&1 | tee "${log_file}"

  echo
  echo "OK: ${script_name}"
}

run_with_model_source () {
  local src="$1"
  local step_name="$2"
  local script_name="$3"

  local old_model_source="${MODEL_SOURCE}"

  MODEL_SOURCE="${src}"
  export MODEL_SOURCE

  run_step "${step_name}" "${script_name}"

  MODEL_SOURCE="${old_model_source}"
  export MODEL_SOURCE
}

run_downstream_single () {
  local src="$1"
  local tgt="$2"

  local old_model_source="${MODEL_SOURCE}"
  local old_diag_target="${DIAG_TARGET}"

  MODEL_SOURCE="${src}"
  DIAG_TARGET="${tgt}"
  export MODEL_SOURCE
  export DIAG_TARGET

  run_step "05_downstream_${src}_${tgt}" "downstream.R"

  MODEL_SOURCE="${old_model_source}"
  DIAG_TARGET="${old_diag_target}"
  export MODEL_SOURCE
  export DIAG_TARGET
}

## -----------------------------
## Función para auditar configuración
## -----------------------------

print_config () {
  echo "============================================================"
  echo "PIPELINE MULTIÓMICO"
  echo "============================================================"
  echo "RUN_ID             = ${RUN_ID}"
  echo "ANALYSIS_NAME      = ${ANALYSIS_NAME}"
  echo "OUTDIR_BASE        = ${OUTDIR_BASE}"
  echo "ANALYSIS_DIR       = ${ANALYSIS_DIR}"
  echo "Q_TX               = ${Q_TX}"
  echo "Q_PR               = ${Q_PR}"
  echo "Q_ME               = ${Q_ME}"
  echo "Q_CL               = ${Q_CL}"
  echo "PANEL_CHOICE       = ${PANEL_CHOICE}"
  echo "MODEL_SOURCE       = ${MODEL_SOURCE}"
  echo "DIAG_TARGET        = ${DIAG_TARGET}"
  echo "CLINICAL_INCLUDE   = ${CLINICAL_INCLUDE}"
  echo "CLINICAL_EXCLUDE   = ${CLINICAL_EXCLUDE}"
  echo "LOG_DIR            = ${LOG_DIR}"
  echo "============================================================"
}

## -----------------------------
## Leer plantilla línea a línea
## -----------------------------

RUN_ID=0
MODEL_SOURCE_ALLOWED=(NEW OLD BOTH)
DIAG_TARGET_ALLOWED=(BAYES_INPUT PANEL_FINAL BOTH)

while IFS=$'\t' read -r \
  ANALYSIS_NAME \
  OUTDIR_BASE \
  Q_TX \
  Q_PR \
  Q_ME \
  Q_CL \
  PANEL_CHOICE \
  CLINICAL_INCLUDE \
  CLINICAL_EXCLUDE \
  MODEL_SOURCE \
  DIAG_TARGET \
  EXTRA_COLS
do

  ## Saltar líneas vacías
  [[ -z "${ANALYSIS_NAME// }" ]] && continue

  ## Saltar comentarios
  [[ "${ANALYSIS_NAME}" =~ ^# ]] && continue

  ## Saltar cabecera
  if [[ "${ANALYSIS_NAME}" == "ANALYSIS_NAME" ]]; then
    continue
  fi

  RUN_ID=$((RUN_ID + 1))

  ## Defaults seguros
  OUTDIR_BASE="${OUTDIR_BASE:-.}"
  Q_TX="${Q_TX:-0.97}"
  Q_PR="${Q_PR:-0.40}"
  Q_ME="${Q_ME:-0.97}"
  Q_CL="${Q_CL:-0.70}"

  ## En esta versión, FINAL es más claro que RELAXED.
  ## Si dejas RELAXED, el downstream lo recodifica a FINAL.
  PANEL_CHOICE="${PANEL_CHOICE:-FINAL}"

  CLINICAL_INCLUDE="${CLINICAL_INCLUDE:-}"
  CLINICAL_EXCLUDE="${CLINICAL_EXCLUDE:-}"

  ## Para corrida completa fresca por defecto, usar NEW.
  ## BOTH ahora ejecuta OLD y NEW desde cero.
  MODEL_SOURCE="${MODEL_SOURCE:-NEW}"

  ## Por defecto ejecuta los dos diagnósticos:
  ## 1) features BRMS/MOFA
  ## 2) panel final multimodelo
  DIAG_TARGET="${DIAG_TARGET:-BOTH}"

  MODEL_SOURCE="${MODEL_SOURCE^^}"
  DIAG_TARGET="${DIAG_TARGET^^}"
  PANEL_CHOICE="${PANEL_CHOICE^^}"

  if [[ ! " ${MODEL_SOURCE_ALLOWED[*]} " =~ " ${MODEL_SOURCE} " ]]; then
    echo "ERROR: MODEL_SOURCE debe ser NEW, OLD o BOTH. Valor actual: ${MODEL_SOURCE}"
    exit 1
  fi

  if [[ ! " ${DIAG_TARGET_ALLOWED[*]} " =~ " ${DIAG_TARGET} " ]]; then
    echo "ERROR: DIAG_TARGET debe ser BAYES_INPUT, PANEL_FINAL o BOTH. Valor actual: ${DIAG_TARGET}"
    exit 1
  fi

  PIPELINE_MODEL_SOURCE="${MODEL_SOURCE}"
  PIPELINE_DIAG_TARGET="${DIAG_TARGET}"

  if [[ -z "${ANALYSIS_NAME}" ]]; then
    echo "ERROR: ANALYSIS_NAME vacío en corrida ${RUN_ID}"
    exit 1
  fi

  if [[ -z "${CLINICAL_INCLUDE}" ]]; then
    echo "ERROR: CLINICAL_INCLUDE vacío en corrida ${RUN_ID} (${ANALYSIS_NAME})"
    exit 1
  fi

  ANALYSIS_DIR="${OUTDIR_BASE}/${ANALYSIS_NAME}"
  LOG_DIR="${ANALYSIS_DIR}/logs"

  mkdir -p "${LOG_DIR}"

  ## Exportar variables para R
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
  export DIAG_TARGET

  print_config

  env | grep -E "ANALYSIS_NAME|OUTDIR_BASE|Q_TX|Q_PR|Q_ME|Q_CL|CLINICAL_INCLUDE|CLINICAL_EXCLUDE|PANEL_CHOICE|MODEL_SOURCE|DIAG_TARGET" \
    > "${LOG_DIR}/run_environment.txt"

  echo "${TEMPLATE_FILE}" > "${LOG_DIR}/template_file_used.txt"

  ## ============================================================
  ## Ejecutar pipeline completo
  ## ============================================================

  run_step "01_preprocess_data"      "preprocess_data.R"
  run_step "02_integracion_mofa"     "INTEGRACION.R"

  ## ------------------------------------------------------------
  ## Modelamiento
  ## ------------------------------------------------------------
  ## NEW  -> crea ANALYSIS_DIR/modelamiento/
  ## OLD  -> crea ANALYSIS_DIR/modelado_anterior/
  ## BOTH -> crea ambos desde cero, en el orden: OLD y luego NEW

  case "${PIPELINE_MODEL_SOURCE}" in
    NEW)
      run_step "03_modelamiento_NEW" "modelamiento2.R"
      ;;

    OLD)
      run_step "03_modelamiento_OLD" "modelamiento_anterior.R"
      ;;

    BOTH)
      run_step "03a_modelamiento_OLD" "modelamiento_anterior.R"
      run_step "03b_modelamiento_NEW" "modelamiento2.R"
      ;;
  esac

  ## ------------------------------------------------------------
  ## Panel final multimodelo
  ## ------------------------------------------------------------
  ## Se ejecuta por fuente explícita para que los logs no se mezclen.

  case "${PIPELINE_MODEL_SOURCE}" in
    NEW)
      run_with_model_source "NEW" "04_panel_multimodelo_NEW" "panel_final_multimodelo.R"
      ;;

    OLD)
      run_with_model_source "OLD" "04_panel_multimodelo_OLD" "panel_final_multimodelo.R"
      ;;

    BOTH)
      run_with_model_source "OLD" "04a_panel_multimodelo_OLD" "panel_final_multimodelo.R"
      run_with_model_source "NEW" "04b_panel_multimodelo_NEW" "panel_final_multimodelo.R"
      ;;
  esac

  ## ------------------------------------------------------------
  ## Downstream
  ## ------------------------------------------------------------
  ## BAYES_INPUT existe para NEW porque modelamiento2.R guarda
  ## features_finales_integradas_mofa_brms.*
  ## Para OLD no existe esa tabla; por eso OLD se diagnostica como PANEL_FINAL.

  case "${PIPELINE_MODEL_SOURCE}:${PIPELINE_DIAG_TARGET}" in
    NEW:BAYES_INPUT)
      run_downstream_single "NEW" "BAYES_INPUT"
      ;;

    NEW:PANEL_FINAL)
      run_downstream_single "NEW" "PANEL_FINAL"
      ;;

    NEW:BOTH)
      run_downstream_single "NEW" "BAYES_INPUT"
      run_downstream_single "NEW" "PANEL_FINAL"
      ;;

    OLD:BAYES_INPUT)
      echo "AVISO: OLD no genera features_finales_integradas_mofa_brms.*; se omite BAYES_INPUT para OLD."
      ;;

    OLD:PANEL_FINAL|OLD:BOTH)
      run_downstream_single "OLD" "PANEL_FINAL"
      ;;

    BOTH:BAYES_INPUT)
      run_downstream_single "NEW" "BAYES_INPUT"
      echo "AVISO: OLD no genera features_finales_integradas_mofa_brms.*; se omite BAYES_INPUT para OLD."
      ;;

    BOTH:PANEL_FINAL)
      run_downstream_single "OLD" "PANEL_FINAL"
      run_downstream_single "NEW" "PANEL_FINAL"
      ;;

    BOTH:BOTH)
      run_downstream_single "OLD" "PANEL_FINAL"
      run_downstream_single "NEW" "BAYES_INPUT"
      run_downstream_single "NEW" "PANEL_FINAL"
      ;;
  esac

  MODEL_SOURCE="${PIPELINE_MODEL_SOURCE}"
  DIAG_TARGET="${PIPELINE_DIAG_TARGET}"
  export MODEL_SOURCE
  export DIAG_TARGET

  echo
  echo "============================================================"
  echo "CORRIDA TERMINADA"
  echo "ANALYSIS_NAME = ${ANALYSIS_NAME}"
  echo "Resultados en = ${ANALYSIS_DIR}"
  echo "Logs en       = ${LOG_DIR}"
  echo "============================================================"

done < "${TEMPLATE_FILE}"

echo
echo "============================================================"
echo "TODAS LAS CORRIDAS TERMINARON"
echo "Número de corridas ejecutadas: ${RUN_ID}"
echo "============================================================"
#!/usr/bin/env bash
set -Eeuo pipefail
# Orquestador de ejecución. No modifica algoritmos ni resultados ya creados.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
usage() {
  cat <<'EOF'
USO:
  bash scripts/orchestration/run_metatest_vigente.sh [--phase biomarkers|upstream|selection|all] [OPCIONES]

Opciones:
  --results-root DIR  Carpeta de resultados externa al repositorio (default ../results)
  --data-root DIR     Carpeta de datos externa al repositorio (default ../data)
  --upstream DIR      Raíz upstream_comun (default RESULTS_ROOT/comparacion_modelamientos/upstream_comun)
  --selection-run DIR Selección congelada; si falta, resolver última selección válida
  --from N           Primera etapa de biomarcadores (01..10,16..18)
  --to N             Última etapa de biomarcadores (01..10,16..18)
  --with-console     Añadir exploración etapa 19 (solo consola, no guarda resultados)
  --dry-run          Mostrar plan y requisitos sin ejecutar R
  --no-check         No verificar paquetes R (no recomendado)
  -h, --help         Ayuda

IMPORTANTE: --phase all ejecuta el upstream completo definido en sus plantillas (potencialmente 201 escenarios).
Las etapas 20-22 están diseñadas para consola gráfica de RStudio: no se ejecutan automáticamente.
EOF
}
PHASE=biomarkers; DRY=0; CHECK=1; WITH_CONSOLE=0; FROM=01; TO=18
RESULTS_ROOT="${RESULTS_ROOT:-$REPO/../results}"
DATA_ROOT="${DATA_ROOT:-$REPO/../data}"
UPSTREAM_DIR="${UPSTREAM_DIR:-}"; SELECTION_RUN="${SELECTION_RUN:-}"
while (($#)); do
  case "$1" in
    --phase|--results-root|--data-root|--upstream|--selection-run|--from|--to)
      (($#>=2)) || { echo "Falta valor para $1" >&2; exit 2; }
      case "$1" in
        --phase) PHASE="$2";; --results-root) RESULTS_ROOT="$2";;
        --data-root) DATA_ROOT="$2";; --upstream) UPSTREAM_DIR="$2";;
        --selection-run) SELECTION_RUN="$2";; --from) FROM="$2";; --to) TO="$2";;
      esac; shift 2;;
    --with-console) WITH_CONSOLE=1; shift;;
    --dry-run) DRY=1; shift;;
    --no-check) CHECK=0; shift;;
    --help|-h) usage; exit 0;;
    *) echo "Opción desconocida: $1" >&2; usage; exit 2;;
  esac
done
[[ "$PHASE" =~ ^(all|upstream|selection|biomarkers)$ ]] || { echo 'Phase inválida' >&2; exit 2; }
case "$FROM" in 1) FROM=01;; 2) FROM=02;; 3) FROM=03;; 4) FROM=04;; 5) FROM=05;; 6) FROM=06;; 7) FROM=07;; 8) FROM=08;; 9) FROM=09;; esac
case "$TO" in 1) TO=01;; 2) TO=02;; 3) TO=03;; 4) TO=04;; 5) TO=05;; 6) TO=06;; 7) TO=07;; 8) TO=08;; 9) TO=09;; esac
STEPS=(01 02 03 04 05 06 07 08 09 10 16 17 18)
FILES=(
  01_run_bidirectional_framework_selected_scenarios.R
  02_run_latent_supervised_audit.R
  03_run_mofa_native_enrichment_part1.R
  04_run_class_pathway_bridge.R
  05_run_crossomic_bridge.R
  06_run_observed_train_corroboration.R
  07_run_test_transportability_audit.R
  08_define_biological_programs.R
  09_build_biological_synthesis_figures.R
  10_build_final_predictive_panel.R
  16_sparse_group_lasso_hierarchical.R
  17_evidenceA_plsda_pca_validation_SCATTER_INTERSECTIONS_LOADINGS_FINAL.R
  18_build_final_consensus_panel.R
)
seenfrom=0; seenTo=0
for s in "${STEPS[@]}"; do [[ "$s" == "$FROM" ]] && seenfrom=1; [[ "$s" == "$TO" ]] && seenTo=1; done
((seenfrom && seenTo)) || { echo "--from/--to debe ser una etapa de: ${STEPS[*]}" >&2; exit 2; }
((10#$FROM <= 10#$TO)) || { echo '--from mayor que --to' >&2; exit 2; }
absolute() { mkdir -p "$1"; (cd "$1" && pwd -P); }
# El dry-run no crea carpetas de salida si no existen.
resolve() { if ((DRY)); then realpath -m "$1"; else absolute "$1"; fi; }
RESULTS_ROOT="$(resolve "$RESULTS_ROOT")"
DATA_ROOT="$(realpath -m "$DATA_ROOT")"
if [[ -z "$UPSTREAM_DIR" ]]; then UPSTREAM_DIR="$RESULTS_ROOT/comparacion_modelamientos/upstream_comun"; fi
UPSTREAM_DIR="$(realpath -m "$UPSTREAM_DIR")"
COMPARISON_ROOT="$(dirname "$UPSTREAM_DIR")"
if [[ -n "$SELECTION_RUN" ]]; then SELECTION_RUN="$(realpath -m "$SELECTION_RUN")"; fi
export RESULTS_ROOT DATA_ROOT UPSTREAM_DIR COMPARISON_ROOT SELECTION_RUN
export PIPELINE_R_CORES="${PIPELINE_R_CORES:-2}" BRMS_INTERNAL_CORES="${BRMS_INTERNAL_CORES:-1}"
export BLAS_THREADS="${BLAS_THREADS:-1}" CMDSTAN_OPTIMIZATION_LEVEL="${CMDSTAN_OPTIMIZATION_LEVEL:-0}"
export OMP_NUM_THREADS="$BLAS_THREADS" OPENBLAS_NUM_THREADS="$BLAS_THREADS"
export MKL_NUM_THREADS="$BLAS_THREADS" VECLIB_MAXIMUM_THREADS="$BLAS_THREADS"
export NUMEXPR_NUM_THREADS="$BLAS_THREADS" STAN_NUM_THREADS="$BLAS_THREADS"
export METATEST_REPO_ROOT="$REPO" R_PROFILE_USER="$REPO/scripts/reproducibility/activate_renv.Rprofile"
# Usar entorno del proyecto también con Rscript no interactivo.
LOG_ROOT="$RESULTS_ROOT/logs/metatest_vigente_$(date -u +%Y%m%dT%H%M%SZ)_$$"
if ((DRY)); then
  echo "PLAN | fase=$PHASE | repo=$REPO | data=$DATA_ROOT | results=$RESULTS_ROOT"
  echo "UPSTREAM_DIR=$UPSTREAM_DIR"
  echo "SELECTION_RUN=${SELECTION_RUN:-[se buscará último válido]}"
  case "$PHASE" in upstream|all) echo "Upstream: scripts/run_upstream_only_local.sh (usa plantillas actuales; NO ejecuta modelos posteriores)";; esac
  case "$PHASE" in selection|all) echo "Selector local: scripts/selection/select_best_upstream_scenario_downstream.R (debe existir en repo local)";; esac
  if [[ "$PHASE" == all || "$PHASE" == biomarkers ]]; then
    for i in "${!STEPS[@]}"; do s="${STEPS[$i]}"; if ((10#$s>=10#$FROM && 10#$s<=10#$TO)); then echo "$s: scripts/biomarker_selection/${FILES[$i]}"; fi; done
    ((WITH_CONSOLE)) && echo '19: scripts/biomarker_selection/19_probe_final_small_panel_console.R'
  fi
  echo "Logs en: $LOG_ROOT"; exit 0
fi
mkdir -p "$LOG_ROOT"
RUNLOG="$LOG_ROOT/estado.tsv"
printf 'utc_start\tstep\tscript\tsha256\tstatus\n' > "$RUNLOG"
{
  date -u --iso-8601=seconds; printf 'REPO=%s\n' "$REPO"
  git rev-parse HEAD 2>/dev/null || true
  git status --short 2>/dev/null || true
  printf 'DATA_ROOT=%s\nRESULTS_ROOT=%s\nUPSTREAM_DIR=%s\n' "$DATA_ROOT" "$RESULTS_ROOT" "$UPSTREAM_DIR"
  printf 'Rscript=%s\n' "$(command -v Rscript || echo MISSING)"
} > "$LOG_ROOT/ejecucion.txt"
command -v Rscript >/dev/null || { echo 'ERROR: No se encontró Rscript' >&2; exit 1; }
run_step() {
  local step="$1" script="$2" ts sum safe
  [[ -f "$REPO/$script" ]] || { echo "ERROR: falta $script" >&2; exit 1; }
  safe="${step}_$(basename "$script" .R)"
  ts="$(date -u +%FT%TZ)"; sum="$(sha256sum "$REPO/$script" | cut -d' ' -f1)"
  echo "===== $ts STEP $step $script ====="
  if (if [[ "$script" == *.sh ]]; then bash "$REPO/$script"; else Rscript "$REPO/$script"; fi) > >(tee "$LOG_ROOT/$safe.log") 2> >(tee "$LOG_ROOT/$safe.err" >&2); then
    printf '%s\t%s\t%s\t%s\tPASS\n' "$ts" "$step" "$script" "$sum" >> "$RUNLOG"
  else
    printf '%s\t%s\t%s\t%s\tFAIL\n' "$ts" "$step" "$script" "$sum" >> "$RUNLOG"
    echo "ERROR: falló $step. Logs: $LOG_ROOT" >&2; exit 1
  fi
}
if ((CHECK)); then
  Rscript "$REPO/scripts/reproducibility/check_renv_vigente.R" | tee "$LOG_ROOT/preflight.log"
fi
if [[ "$PHASE" == upstream || "$PHASE" == all ]]; then
  [[ -d "$DATA_ROOT" ]] || { echo "Falta DATA_ROOT=$DATA_ROOT" >&2; exit 1; }
  # El upstream histórico requiere inputs y recursos ya preparados.
  run_step upstream scripts/run_upstream_only_local.sh
fi
if [[ "$PHASE" == selection || "$PHASE" == all ]]; then
  [[ -f "$REPO/scripts/selection/select_best_upstream_scenario_downstream.R" ]] || {
    echo 'ERROR: falta el selector local scripts/selection/select_best_upstream_scenario_downstream.R (no estaba en scripts.zip).' >&2; exit 1;
  }
  run_step selection scripts/selection/select_best_upstream_scenario_downstream.R
fi
if [[ "$PHASE" == all || "$PHASE" == biomarkers ]]; then
  if [[ -z "$SELECTION_RUN" ]]; then
    select_root="$UPSTREAM_DIR/downstream/seleccion_mejor_escenario_train"
    [[ -d "$select_root" ]] || { echo "No existe directorio de selección: $select_root" >&2; exit 1; }
    newest=""
    while IFS= read -r -d '' cand; do
      [[ -s "$cand/SELECTED_SCENARIOS_BY_ANALYSIS.csv" ]] || continue
      if [[ -z "$newest" || "$cand" -nt "$newest" ]]; then newest="$cand"; fi
    done < <(find "$select_root" -mindepth 1 -maxdepth 1 -type d -print0)
    [[ -n "$newest" ]] || { echo "No existe selección completa en $select_root" >&2; exit 1; }
    SELECTION_RUN="$(realpath "$newest")"
  fi
  [[ -s "$SELECTION_RUN/SELECTED_SCENARIOS_BY_ANALYSIS.csv" ]] || {
    echo "Selección incompleta: $SELECTION_RUN" >&2; exit 1;
  }
  export SELECTION_RUN
  printf 'SELECTION_RUN=%s\n' "$SELECTION_RUN" >> "$LOG_ROOT/ejecucion.txt"
  for i in "${!STEPS[@]}"; do
    s="${STEPS[$i]}"
    if ((10#$s>=10#$FROM && 10#$s<=10#$TO)); then
      run_step "$s" "scripts/biomarker_selection/${FILES[$i]}"
    fi
  done
  if ((WITH_CONSOLE)); then
    run_step 19 scripts/biomarker_selection/19_probe_final_small_panel_console.R
  fi
fi
echo "PIPELINE COMPLETADO. Resumen: $RUNLOG"

#!/usr/bin/env bash
# Recolección sólo de evidencia/logs y resúmenes; no datos brutos ni RDS.
set -Eeuo pipefail
OUT="$(realpath -e "${1:?Uso: collect.sh OUTDIR [BASELINE_PILOT_ROOT]}")"
BASELINE="${2:-}"
[[ -d "$OUT/run_reports" ]] || { echo 'No hay run_reports' >&2; exit 1; }
exec 8>"$OUT/.portable_run.lock"
flock -n 8 || { echo 'Todavía se está ejecutando un escenario' >&2; exit 1; }
for state in "$OUT"/run_reports/submission_*/array_job_id.txt; do
  [[ -s "$state" ]] || continue
  if command -v squeue >/dev/null 2>&1; then
    job="$(cat "$state")"
    queue="$(squeue -h -u "${USER:?}" -o '%i')" || { echo 'No se pudo consultar Slurm; reintenta la recolección' >&2; exit 3; }
    active="$(printf '%s\n' "$queue" | awk -v j="$job" '$0==j || index($0,j "_")==1')"
    [[ -z "$active" ]] || { echo "Array $job todavía activo"; exit 3; }
  fi
done
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
TMP="$(mktemp -d "$OUT/evidence_${STAMP}.XXXXXX")"
mkdir -p "$TMP/reports" "$TMP/summary_tables"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
if [[ -f "$REPO/PORTABLE_SHA256SUMS.txt" ]]; then cp -p "$REPO/PORTABLE_SHA256SUMS.txt" "$TMP/reports/"; fi
# No empaquetar un rootfs/backend temporal si una limpieza quedó interrumpida.
while IFS= read -r file; do
  rel="${file#"$OUT/run_reports/"}"
  mkdir -p "$TMP/reports/$(dirname "$rel")"
  cp -p "$file" "$TMP/reports/$rel"
done < <(find "$OUT/run_reports" -type d -name '.metatest_runtime_*' -prune -o -type f -print)
printf 'report\tscenario\tstan_o\tresume\tstatus\ttotal_seconds\tpipeline_time\tmaxrss_kb\n' > "$TMP/timing_summary.tsv"
while IFS= read -r final; do
  report="$(dirname "$(dirname "$final")")"
  get() { awk -F= -v k="$1" '$1==k {print substr($0,length(k)+2);exit}' "$final"; }
  timefile="$report/job_artifacts/05_time_and_stderr.txt"
  wall="$(awk '/Elapsed \(wall clock\) time/ {sub(/^.*\):[[:space:]]*/,""); print; exit}' "$timefile" 2>/dev/null || true)"
  rss="$(awk -F: '/Maximum resident set size/ {gsub(/[[:space:]]/,"",$2);print $2;exit}' "$timefile" 2>/dev/null || true)"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'     "$(basename "$report")" "$(get SCENARIO)" "$(get STAN_O)" "$(get RESUME)"     "$(get METATEST_PORTABLE_RUN)" "$(get WALLTIME_SECONDS)" "$wall" "$rss" >> "$TMP/timing_summary.tsv"
done < <(find "$OUT/run_reports" -name 99_FINAL_STATUS.env -type f | sort)
gather_tables() {
  local source="$1" dest="$2" file rel
  [[ -d "$source" ]] || return 0
  while IFS= read -r file; do
    rel="${file#"$source/"}"
    mkdir -p "$dest/$(dirname "$rel")"
    cp -p "$file" "$dest/$rel"
  done < <(find "$source" -type f \(       -name estadisticas_modelo_global.csv -o -name run_config.csv       -o -name cv_best_stable_grid.csv -o -name cv_stability_summary_candidates.csv       -o -name cv_repeated_summary_candidates.csv -o -name test_all_grids_summary.csv       -o -name resumen_corridas.tsv \))
}
gather_tables "$OUT/comparacion_modelamientos" "$TMP/summary_tables"
if [[ -n "$BASELINE" ]]; then
  BASELINE="$(realpath -e "$BASELINE")"
  mkdir -p "$TMP/baseline"
  for task in "$BASELINE"/tasks/task_*; do
    [[ -d "$task" ]] || continue
    dest="$TMP/baseline/$(basename "$task")"
    mkdir -p "$dest"
    cp -a "$task/job_report" "$task/pipeline_results/job_artifacts" "$dest/"
    gather_tables "$task/pipeline_results/comparacion_modelamientos" "$dest/summary_tables"
  done
fi
(
  cd "$TMP"
  sha256sum timing_summary.tsv
  find reports summary_tables -type f -print | sort | while IFS= read -r f; do sha256sum "$f"; done
  if [[ -d baseline ]]; then
    find baseline -type f -print | sort | while IFS= read -r f; do sha256sum "$f"; done
  fi
) > "$TMP/SHA256SUMS_EVIDENCE.txt"
PACKAGE="$OUT/PARA_CHAT_PORTABLE_${STAMP}.tar.gz"
tar -czf "$PACKAGE" -C "$TMP" .
cat "$TMP/timing_summary.tsv"
printf 'PAQUETE_PARA_CHAT=%s\n' "$PACKAGE"
printf 'NOTA=un PASS funcional no certifica equivalencia estadistica O0/O1\n'

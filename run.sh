#!/usr/bin/env bash
set -Eeuo pipefail
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
REPO="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
usage() {
  cat <<'HELP'
METATEST portable (Linux/Bash 4+)
./run.sh --backend enroot --plantillas templates/current/plantilla_1_todas_clinicas_barrido.tsv \
  --scenario P1_TODAS_qtx90_qpr45_qme90_qcl45 --data ../data --outdir ../results/bench_O1 \
  --assets ../assets --stan-o 1 --slurm
--backend local|enroot      local usa el R ya instalado; enroot usa los activos sellados.
--plantillas FILE           Repetible. Omisión: las tres plantillas current.
--scenario ID              Repetible. Omisión: todos los análisis de las plantillas elegidas.
--fila N                   Fila de análisis (1=q00), excluye comentarios/cabecera; una plantilla.
--pilot q90                Los tres q90 balanceados (no 201 escenarios).
--data DIR                 Default: ../data respecto al repositorio.
--outdir DIR               Default: ../results respecto al repositorio.
--assets DIR               Default: ../assets; necesario sólo para enroot.
--stan-o 0|1               Default 0. O1 experimental: rebuild sólo del CmdStan temporal.
--slurm                    Array secuencial, una tarea por escenario (MaxJobs=1).
--partition NAME           Default cpu-dev. --cpus N (8), --mem VALUE (64G), --time HH:MM:SS (8h).
--resume                   Sólo mismo perfil, código e inputs; no sirve como benchmark limpio.
--list                     Lista el plan sin ejecutar ni crear archivos.
--dry-run                  Muestra rutas y plan sin ejecutar ni crear archivos.
--check                    Verifica Git, codigo y rutas sin ejecutar ni enviar jobs.
Rutas relativas explícitas: respecto al directorio desde donde invocas el comando.
No instala R, paquetes, Enroot ni descarga datos. No publica ni cambia ningún tag.
HELP
}
BACKEND=local DATA="$REPO/../data" OUT="$REPO/../results" ASSETS="$REPO/../assets"
OPT=0 SUBMIT=0 LIST=0 DRY=0 CHECK=0 RESUME=0 ROW="" PILOT=""
PARTITION=cpu-dev CPUS=8 MEM=64G LIMIT=08:00:00
TEMPLATES=() SELECTED=()
while (($#)); do
  case "$1" in
    --help|-h) usage; exit 0;;
    --slurm) SUBMIT=1; shift;;
    --list) LIST=1; shift;;
    --dry-run) DRY=1; shift;;
    --check) CHECK=1; shift;;
    --resume) RESUME=1; shift;;
    --backend|--data|--outdir|--assets|--stan-o|--partition|--cpus|--mem|--time|--fila|--pilot|--scenario|--plantillas|--plantilla)
      (($# >= 2)) || die "falta valor: $1"
      case "$1" in
        --backend) BACKEND="$2";; --data) DATA="$2";; --outdir) OUT="$2";;
        --assets) ASSETS="$2";; --stan-o) OPT="$2";; --partition) PARTITION="$2";;
        --cpus) CPUS="$2";; --mem) MEM="$2";; --time) LIMIT="$2";;
        --fila) ROW="$2";; --pilot) PILOT="$2";; --scenario) SELECTED+=("$2");;
        --plantillas|--plantilla) TEMPLATES+=("$2");;
      esac
      shift 2;;
    *) die "opción desconocida: $1";;
  esac
done
[[ "$BACKEND" == local || "$BACKEND" == enroot ]] || die "backend inválido"
[[ "$OPT" == 0 || "$OPT" == 1 ]] || die "stan-o debe ser 0 o 1"
[[ "$CPUS" =~ ^[1-9][0-9]*$ ]] || die "CPUS inválido"
((CPUS >= 2)) || die "este perfil usa PIPELINE_R_CORES=2; pide al menos 2 CPUs"
if [[ "$BACKEND" == local && "$OPT" != 0 ]]; then
  die "O1 automático sólo en contenedor efímero; no se modifica el CmdStan local"
fi
if (("${#TEMPLATES[@]}" == 0)); then
  TEMPLATES=("$REPO"/templates/current/*.tsv)
fi
for i in "${!TEMPLATES[@]}"; do
  TEMPLATES[$i]="$(realpath -e "${TEMPLATES[$i]}")"
  [[ "${TEMPLATES[$i]}" == "$REPO/"* && -f "${TEMPLATES[$i]}" ]] ||
    die "la plantilla debe estar dentro del repo: ${TEMPLATES[$i]}"
done
[[ -z "$ROW" || "${#TEMPLATES[@]}" == 1 ]] || die "--fila requiere una sola plantilla"
[[ -z "$ROW" || "$ROW" =~ ^[1-9][0-9]*$ ]] || die "--fila inválida"
[[ -z "$ROW" || ("${#SELECTED[@]}" == 0 && -z "$PILOT") ]] || die "no mezcles fila y scenario/pilot"
[[ -z "$PILOT" || "$PILOT" == q90 ]] || die "piloto desconocido"
if [[ "$PILOT" == q90 ]]; then
  (("${#SELECTED[@]}" == 0)) || die "no mezcles --pilot y --scenario"
  SELECTED=(P1_TODAS_qtx90_qpr45_qme90_qcl45 P2_CON_FLI_qtx90_qpr45_qme90_qcl45 P3_SIN_FLI_qtx90_qpr45_qme90_qcl45)
fi
PLAN_TEMPLATE=() PLAN_SCENARIO=() PLAN_ROOT=()
declare -A SEEN=() FOUND=()
for template in "${TEMPLATES[@]}"; do
  awk -F '\t' '
    {sub(/\r$/, "")}
    /^#/ || NF==0 {next}
    $1=="ANALYSIS_NAME" {if(NF!=11) exit 2; h++; next}
    {if(NF!=11 || seen[$1]++) exit 3; n++}
    END {if(h!=1 || n<1) exit 4}
  ' "$template" || die "plantilla inválida: $template"
  row=0
  while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%$'\r'}"
    [[ -z "$line" || "$line" == \#* || "$line" == ANALYSIS_NAME$'\t'* ]] && continue
    row=$((row+1))
    IFS=$'\x1f' read -r scenario root rest <<< "${line//$'\t'/$'\x1f'}"
    root="${root#./}"
    [[ "$scenario" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && "$root" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] ||
      die "ID/raíz inseguros"
    [[ -z "$ROW" || "$row" == "$ROW" ]] || continue
    if (("${#SELECTED[@]}" > 0)); then
      include=0
      for wanted in "${SELECTED[@]}"; do [[ "$scenario" == "$wanted" ]] && include=1; done
      ((include)) || continue
    fi
    [[ -z "${SEEN[$scenario]:-}" ]] || die "escenario duplicado: $scenario"
    SEEN[$scenario]=1 FOUND[$scenario]=1
    PLAN_TEMPLATE+=("${template#"$REPO/"}")
    PLAN_SCENARIO+=("$scenario") PLAN_ROOT+=("$root")
  done < "$template"
done
(("${#PLAN_SCENARIO[@]}" > 0)) || die "selección vacía"
for wanted in "${SELECTED[@]}"; do [[ -n "${FOUND[$wanted]:-}" ]] || die "no existe escenario: $wanted"; done
DATA="$(realpath -m "$DATA")" OUT="$(realpath -m "$OUT")" ASSETS="$(realpath -m "$ASSETS")"
for path in "$DATA" "$OUT" "$ASSETS" "$REPO"; do
  [[ "$path" != *$'\n'* && "$path" != *$'\t'* ]] || die "ruta con TAB/newline"
done
print_plan() {
  printf 'template\tscenario\troot\n'
  for i in "${!PLAN_SCENARIO[@]}"; do
    printf '%s\t%s\t%s\n' "${PLAN_TEMPLATE[$i]}" "${PLAN_SCENARIO[$i]}" "${PLAN_ROOT[$i]}"
  done
}
if ((LIST)); then print_plan; exit 0; fi
printf 'REPO=%s\nDATA=%s\nOUTDIR=%s\nASSETS=%s\nBACKEND=%s\nSTAN_O=%s\nSCENARIOS=%s\nRESUME=%s\n'   "$REPO" "$DATA" "$OUT" "$ASSETS" "$BACKEND" "$OPT" "${#PLAN_SCENARIO[@]}" "$RESUME"
if ((DRY)); then print_plan; exit 0; fi
[[ -d "$DATA" ]] || die "no existe --data: $DATA"
[[ "$OUT" != / && "$OUT" != "$REPO" && "$OUT" != "$DATA" && "$OUT" != "$ASSETS" && "$OUT" != "$REPO/"* && "$OUT" != "$DATA/"* && "$OUT" != "$ASSETS/"* ]] ||
  die "outdir debe estar fuera del repositorio/datos/assets"
source "$REPO/scripts/portable/provenance.sh"
metatest_verify_checkout "$REPO" || die "checkout no reproducible para este perfil"
HEAD_COMMIT="$METATEST_CODE_COMMIT"
metatest_print_provenance
if [[ "$BACKEND" == enroot ]]; then
  # Chequeo barato antes de enviar; los hashes completos se verifican en cómputo.
  [[ -s "$ASSETS/metatest-runtime-metatest-reproducible-20260830.sqsh" && -s "$ASSETS/metatest-mofa-backend-metatest-reproducible-20260830.tar" ]] || die "faltan los dos activos sellados en --assets"
fi
if ((CHECK)); then
  printf 'CHECK=PASS\nMODELS_EXECUTED=FALSE\nSLURM_SUBMITTED=FALSE\n'
  printf 'NOTA=check estatico; hashes completos de datos/assets y smoke se verifican dentro del job\n'
  exit 0
fi
if ((SUBMIT)); then
  command -v sbatch >/dev/null || die "falta sbatch"
  [[ -z "${SLURM_JOB_ID:-}" ]] || die "no enviar otro array desde un job"
  mkdir -p "$OUT/run_reports"
  CONTROL="$(mktemp -d "$OUT/run_reports/submission_$(date -u +%Y%m%dT%H%M%SZ).XXXXXX")"
  print_plan > "$CONTROL/plan.tsv"
  COMMON=(--backend "$BACKEND" --data "$DATA" --outdir "$OUT" --assets "$ASSETS" --stan-o "$OPT" --cpus "$CPUS")
  ((RESUME)) && COMMON+=(--resume)
  {
    printf 'REPO=%q\n' "$REPO"
    printf 'PLAN=%q\n' "$CONTROL/plan.tsv"
    printf 'COMMON=('; printf '%q ' "${COMMON[@]}"; printf ')\n'
  } > "$CONTROL/request.sh"
  (
    cd "$REPO"
    sha256sum PORTABLE_SHA256SUMS.txt
    git ls-files -z | xargs -0 sha256sum
    for template in "${TEMPLATES[@]}"; do sha256sum "${template#"$REPO/"}"; done
  ) > "$CONTROL/code.sha256"
  reply="$(sbatch --parsable --job-name=metatest_perf --partition="$PARTITION"     --nodes=1 --ntasks=1 --cpus-per-task="$CPUS" --mem="$MEM" --time="$LIMIT"     --no-requeue --array="0-$((${#PLAN_SCENARIO[@]}-1))%1" --chdir="$REPO"     --output="$CONTROL/slurm_%A_%a.out" --error="$CONTROL/slurm_%A_%a.err"     "$REPO/run.sbatch" --plan-state "$CONTROL/request.sh")"
  job="${reply%%;*}"
  [[ "$job" =~ ^[0-9]+$ ]] || die "respuesta sbatch inesperada: $reply"
  printf '%s\n' "$job" > "$CONTROL/array_job_id.txt"
  printf 'SUBMITTED=PASS\nARRAY_JOB_ID=%s\nCONTROL=%s\n' "$job" "$CONTROL"
  printf 'squeue -j %q\n' "$job"
  printf 'tail -F %q %q\n' "$CONTROL/slurm_${job}_0.out" "$CONTROL/slurm_${job}_0.err"
  printf 'bash %q %q\n' "$REPO/scripts/portable/collect.sh" "$OUT"
  exit 0
fi
if [[ -z "${SLURM_JOB_ID:-}" && "$(hostname -s)" == login* ]]; then
  die "nodo login: añade --slurm; no ejecutar cómputo aquí"
fi
mkdir -p "$OUT/run_reports" "$OUT/.portable_state"
exec 9>"$OUT/.portable_run.lock"
flock -n 9 || die "otro proceso está usando este outdir"
(cd "$DATA" && sha256sum -c "$REPO/reproducibilidad/release/inputs.sha256") || die "inputs diferentes al dataset validado"
SOURCE_SHA="$(cd "$REPO" && sha256sum PORTABLE_SHA256SUMS.txt | awk '{print $1}')"
for i in "${!PLAN_SCENARIO[@]}"; do
  scenario="${PLAN_SCENARIO[$i]}" root="${PLAN_ROOT[$i]}"
  STATE="$OUT/.portable_state/$scenario"
  TEMPLATE_SHA="$(sha256sum "$REPO/${PLAN_TEMPLATE[$i]}" | awk '{print $1}')"
  PROFILE="$(printf '%s\n' "$HEAD_COMMIT" "$SOURCE_SHA" "$TEMPLATE_SHA" "$DATA" "$BACKEND" "$OPT" "$scenario" "${PLAN_TEMPLATE[$i]}" | sha256sum | awk '{print $1}')"
  exists=0
  for variant in upstream_comun modelamiento2 modelamiento3; do
    [[ -e "$OUT/comparacion_modelamientos/$variant/$root/$scenario" ]] && exists=1
  done
  if [[ -d "$STATE" ]] || ((exists)); then
    ((RESUME)) || die "escenario ya intentado/existente; usa otro outdir o --resume: $scenario"
    [[ -f "$STATE/profile.sha256" && "$(cat "$STATE/profile.sha256")" == "$PROFILE" ]] ||
      die "perfil previo ausente o diferente; no se reanuda"
  else
    mkdir "$STATE"
    printf '%s\n' "$PROFILE" > "$STATE/profile.sha256"
  fi
  REPORT="$(mktemp -d "$OUT/run_reports/$(date -u +%Y%m%dT%H%M%SZ)_${scenario}.XXXXXX")"
  metatest_print_provenance > "$REPORT/code_provenance.env"
  export METATEST_MODE=pilot METATEST_REPO_ROOT="$REPO" METATEST_DATA_DIR="$DATA"
  export METATEST_RESULTS_DIR="$OUT" METATEST_ASSET_DIR="$ASSETS" METATEST_RUN_DIR="$REPORT"
  export METATEST_SCENARIO_OVERRIDE="$scenario" METATEST_SCENARIO_ROOT_OVERRIDE="$root"
  export METATEST_TEMPLATE_REL="${PLAN_TEMPLATE[$i]}" METATEST_REPORT_REL="${REPORT#"$OUT/"}"
  export METATEST_STAN_O="$OPT" METATEST_RESUME="$RESUME"
  export METATEST_HOST_CPUS="$CPUS"
  {
    printf 'SCENARIO=%s\nROOT=%s\nTEMPLATE=%s\nBACKEND=%s\nSTAN_O=%s\nRESUME=%s\nDATA=%s\nOUTDIR=%s\nBASE_COMMIT=%s\n'       "$scenario" "$root" "${PLAN_TEMPLATE[$i]}" "$BACKEND" "$OPT" "$RESUME" "$DATA" "$OUT" "$HEAD_COMMIT"
  } > "$REPORT/request.env"
  bash "$REPO/scripts/portable/${BACKEND}.sh"
done

#!/usr/bin/env bash
# Smoke de dependencias nativas: no ejecuta escenarios ni compila Stan.
set -Eeuo pipefail
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
REPO="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)"
DATA="$REPO/../data" OUT="$REPO/../results/prueba_nativa"
while (($#)); do
  case "$1" in
    --data|--outdir)
      (($# >= 2)) || die "falta valor: $1"
      case "$1" in --data) DATA="$2";; --outdir) OUT="$2";; esac
      shift 2;;
    --help|-h)
      printf 'Uso: bash scripts/portable/check_native.sh [--data ../data] [--outdir ../results/prueba_nativa]\n'
      printf 'No instala, no compila, no ejecuta el pipeline ni envia jobs. Escribe un informe fuera del repo.\n'
      exit 0;;
    *) die "opcion desconocida: $1";;
  esac
done
[[ "$(hostname -s)" != login* || -n "${SLURM_JOB_ID:-}" ]] || die 'no ejecutar R en un nodo login'
DATA="$(realpath -e "$DATA")" OUT="$(realpath -m "$OUT")"
for path in "$DATA" "$OUT" "$REPO"; do
  [[ "$path" != *$'\n'* && "$path" != *$'\t'* ]] || die 'ruta con TAB/newline'
done
[[ -d "$DATA" && "$OUT" != / && "$OUT" != "$REPO" && "$OUT" != "$REPO/"* && "$OUT" != "$DATA" && "$OUT" != "$DATA/"* ]] || die 'outdir debe estar fuera del repo y los datos'
source "$REPO/scripts/portable/provenance.sh"
metatest_verify_checkout "$REPO" || die 'checkout no verificable'
mkdir -p "$OUT/run_reports"
exec 9>"$OUT/.portable_run.lock"
flock -n 9 || die 'outdir ocupado por otra ejecucion'
REPORT="$(mktemp -d "$OUT/run_reports/native_check_$(date -u +%Y%m%dT%H%M%SZ).XXXXXX")"
finish() {
  local rc=$?
  trap - EXIT
  printf 'NATIVE_CHECK_EXIT_STATUS=%s\nPIPELINE_EXECUTED=FALSE\nSLURM_SUBMITTED=FALSE\n' "$rc" > "$REPORT/status.env"
  printf 'NATIVE_CHECK_REPORT=%s\n' "$REPORT"
  exit "$rc"
}
trap finish EXIT
metatest_print_provenance > "$REPORT/code_provenance.env"
(cd "$DATA" && sha256sum -c "$REPO/reproducibilidad/release/inputs.sha256") | tee "$REPORT/inputs_check.txt"
source "$REPO/scripts/portable/native_env.sh"
# Redireccionar sin tuberia para conservar las variables exportadas por prepare.
if metatest_native_prepare "$REPO" "$REPORT" > "$REPORT/native_prepare.log" 2>&1; then
  cat "$REPORT/native_prepare.log"
else
  cat "$REPORT/native_prepare.log" >&2
  exit 1
fi
if metatest_native_smoke "$REPO" "$REPORT" > "$REPORT/native_smoke.log" 2>&1; then
  cat "$REPORT/native_smoke.log"
else
  cat "$REPORT/native_smoke.log" >&2
  exit 1
fi
[[ -x /usr/bin/time ]] || die 'falta /usr/bin/time, requerido por el runner nativo; no se instala automaticamente'
printf 'NATIVE_CHECK=PASS\nSTAN_COMPILATION_TEST=NOT_RUN\nPIPELINE_EXECUTED=FALSE\n'

#!/usr/bin/env bash
# Incrementa el barrido inverso ya completado, sin reemplazar plantillas originales.
set -euo pipefail
PKG="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
REPO="${1:-$HOME/Documents/PhD/metatest/metatest}"
MODE="${2:-}"
[[ -d "$REPO" ]] || { echo "Repositorio inexistente: $REPO" >&2; exit 1; }
REPO="$(cd -- "$REPO" && pwd -P)"
[[ -z "$MODE" || "$MODE" == "--dry-run" ]] || { echo "Opcion desconocida: $MODE" >&2; exit 2; }
ORIGINAL="$REPO/scripts/run_upstream_only_local.sh"
[[ -f "$ORIGINAL" ]] || { echo "Falta runner: $ORIGINAL" >&2; exit 1; }
DATA_ROOT="${DATA_ROOT:-$(dirname "$REPO")/data}"
RESULTS_ROOT="${RESULTS_ROOT:-$(dirname "$REPO")/results_barrido_txme_inverso_20261021}"
COMPARISON_ROOT="${COMPARISON_ROOT:-$RESULTS_ROOT/comparacion_modelamientos}"
UP="$COMPARISON_ROOT/upstream_comun"
[[ -d "$DATA_ROOT" ]] || { echo "Falta DATA_ROOT: $DATA_ROOT" >&2; exit 1; }
TS=("$PKG/plantillas/plantilla_1_SOLO_NUEVOS.tsv"
    "$PKG/plantillas/plantilla_2_SOLO_NUEVOS.tsv"
    "$PKG/plantillas/plantilla_3_SOLO_NUEVOS.tsv")
python3 - "${TS[@]}" "$UP" "$MODE" <<'PY'
import csv,sys
from pathlib import Path
prefixes=('P1_TODAS','P2_CON_FLI','P3_SIN_FLI')
rootdirs=('escenarios','escenarios_con_FLI','escenarios_SIN_FLI')
expected_cols=('ANALYSIS_NAME','OUTDIR_BASE','Q_TX','Q_PR','Q_ME','Q_CL','PANEL_CHOICE','CLINICAL_INCLUDE','CLINICAL_EXCLUDE','MODEL_SOURCE','DIAG_TARGET')
expected_new={(80,80),(85,85),(90,90),(93,93),(94,94),(96,96),(97,97),(98,98)}
expected_old={(80,97),(85,95),(90,93),(92,92),(93,90),(95,85),(97,80),(95,95)}
for family,(path,prefix,rootdir) in enumerate(zip(sys.argv[1:4],prefixes,rootdirs),1):
  p=Path(path)
  assert p.is_file(),f'Falta plantilla: {p}'
  rows=list(csv.DictReader((line for line in p.open(encoding='utf-8',newline='') if line.strip() and not line.startswith('#')),delimiter='\t'))
  assert rows and tuple(rows[0])==expected_cols and len(rows)==8, ('Columnas/filas invalidas',p)
  pairs=set();clin=None
  for r in rows:
    tx=round(100*float(r['Q_TX']));me=round(100*float(r['Q_ME']))
    assert (tx,me) in expected_new and (tx,me) not in pairs
    pairs.add((tx,me))
    name=f'{prefix}_qtx{tx:02d}_qpr00_qme{me:02d}_qcl00'
    assert r['ANALYSIS_NAME']==name and r['OUTDIR_BASE']==f'./{rootdir}'
    assert abs(float(r['Q_PR']))<1e-10 and abs(float(r['Q_CL']))<1e-10
    assert (r['PANEL_CHOICE'],r['MODEL_SOURCE'],r['DIAG_TARGET'])==('FINAL','NEW','BOTH')
    fields=(r['CLINICAL_INCLUDE'],r['CLINICAL_EXCLUDE'])
    if clin is None: clin=fields
    assert fields==clin and fields[0], ('Variables clinicas inconsistentes',name)
  assert pairs==expected_new
  existing=Path(sys.argv[4])/rootdir
  if sys.argv[5]!='--dry-run':
    for tx,me in expected_old:
      name=f'{prefix}_qtx{tx:02d}_qpr00_qme{me:02d}_qcl00'
      assert (existing/name/'.UPSTREAM_COMPLETE').is_file(), f'No esta completo el barrido anterior: {existing/name}'
    for r in rows:
      name=r['ANALYSIS_NAME']
      scen=existing/name
      assert not scen.exists() or (scen/'.UPSTREAM_COMPLETE').is_file(),f'Escenario nuevo parcial detectado, revisar antes de reiniciar: {scen}'
  print(f'P{family} OK: 8 nuevos, 8 inversos anteriores. PR y CL sin filtro.')
print('24 escenarios nuevos; 48 para comparacion final (24 anteriores +24 nuevos).')
PY
if [[ "$MODE" == "--dry-run" ]]; then
  echo "PREVALIDACION COMPLETA: ningun modelo ejecutado. Salida prevista: $UP"
  exit 0
fi
# Runner temporal en scripts para preservar determinacion de REPO; no se modifica codigo persistente.
TMP="$(mktemp "$REPO/scripts/.txme_diagonal_XXXXXXXX.sh")"
cleanup(){ rm -f -- "$TMP"; }
trap cleanup EXIT
python3 - "$ORIGINAL" "$TMP" "${TS[@]}" <<'PY'
from pathlib import Path
import sys,re
src,dest=map(Path,sys.argv[1:3]);ts=[str(Path(x).resolve()) for x in sys.argv[3:]]
text=src.read_text(encoding='utf-8')
assert 'ONLY_ANALYSIS_NAME' in text and 'UPSTREAM_ONLY: preprocess_data + INTEGRACION/MOFA completos.' in text
assert re.search(r'(?m)^exit 0\s*$',text), 'No se reconoce el corte upstream-only; abortar'
pattern=r'(?ms)^TEMPLATES=\(\n.*?^\)'
repl='TEMPLATES=(\n'+''.join('  "'+t+'"\n' for t in ts)+')'
new,n=re.subn(pattern,lambda m:repl,text,count=1)
assert n==1,'No se localizo bloque TEMPLATES; no ejecutar'
dest.write_text(new,encoding='utf-8')
print('Runner temporal preparado, el original permanece intacto.')
PY
bash -n "$TMP"
mkdir -p "$RESULTS_ROOT/logs_lanzamiento"
cd "$REPO"
unset ONLY_ANALYSIS_NAME
export DATA_ROOT RESULTS_ROOT COMPARISON_ROOT
export PIPELINE_R_CORES="${PIPELINE_R_CORES:-2}"
export BRMS_INTERNAL_CORES="${BRMS_INTERNAL_CORES:-1}"
export BLAS_THREADS="${BLAS_THREADS:-1}"
export CMDSTAN_OPTIMIZATION_LEVEL="${CMDSTAN_OPTIMIZATION_LEVEL:-0}"
echo "Iniciando 24 integraciones diagonales nuevas; 720 ajustes MOFA previstos."
echo "Se escriben junto a las 24 integraciones inversas existentes en: $UP"
bash "$TMP" 2>&1 | tee -a "$RESULTS_ROOT/logs_lanzamiento/barrido_TXME_diagonal_general.log"
echo "BARRIDO DIAGONAL TERMINADO: $UP"

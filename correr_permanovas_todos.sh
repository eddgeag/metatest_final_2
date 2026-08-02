#!/usr/bin/env bash
set -euo pipefail

SCRIPT="comparar_escenarios_permanova.R"
N_PERM=9999
OBS_TOL=0.95

BASE_DIRS=(
  "./escenarios"
  "./escenarios_con_FLI"
  "./escenarios_SIN_FLI"
)

for BASE in "${BASE_DIRS[@]}"; do
  echo "============================================================"
  echo "Corriendo PERMANOVA para: ${BASE}"
  echo "============================================================"

  if [ ! -d "$BASE" ]; then
    echo "ERROR: No existe la carpeta: $BASE"
    exit 1
  fi

  Rscript "$SCRIPT" "$BASE" ALL AUTO "$N_PERM" "$OBS_TOL"

  echo ""
  echo "Terminado: ${BASE}"
  echo ""
done

echo "============================================================"
echo "TODAS LAS COMPARACIONES TERMINARON CORRECTAMENTE"
echo "============================================================"
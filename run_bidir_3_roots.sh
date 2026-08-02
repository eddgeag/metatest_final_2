#!/usr/bin/env bash
set -euo pipefail

SCRIPT="FRAMEWORK_BIDIRECCIONAL.R"
DIAG="diagnostico_bayes_features"
MODE="both"

ROOTS=(
  "./escenarios"
  "./escenarios_con_FLI"
  "./escenarios_SIN_FLI"
)

mkdir -p logs_bidir

for ROOT in "${ROOTS[@]}"; do
TAG=$(basename "$ROOT")
LOG="logs_bidir/bidir_${TAG}.log"

echo "===================================================="
echo "ROOT: $ROOT"
echo "LOG : $LOG"
echo "===================================================="

Rscript "$SCRIPT" "$ROOT" "$DIAG" "$MODE" 2>&1 | tee "$LOG"
done

echo "OK: corridas terminadas."
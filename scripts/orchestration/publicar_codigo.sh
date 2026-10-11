#!/usr/bin/env bash
set -Eeuo pipefail
# Uso: bash scripts/orchestration/publicar_codigo.sh "Actualizo figuras..."
[[ $# -eq 1 && -n "$1" ]] || { echo 'Uso: publicar_codigo.sh "mensaje"' >&2; exit 2; }
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"; cd "$REPO"
[[ "$(git branch --show-current)" == main ]] || { echo 'Se publica solo desde main' >&2; exit 1; }
# Nunca versionar datos ni resultados.
git add -A -- scripts/plots_presentation scripts/biomarker_selection scripts/orchestration scripts/reproducibility
[[ -z "$(git status --porcelain --untracked-files=no)" ]] && { echo 'Sin cambios preparados'; exit 0; }
# No publicar cambios del pipeline sin una revisión explícita.
git diff --cached --name-only | grep -Ev '^(scripts/plots_presentation/|scripts/biomarker_selection/|scripts/orchestration/|scripts/reproducibility/)' && {
  echo 'Hay cambios preparados fuera de las carpetas permitidas; revisar manualmente.' >&2; exit 1;
} || true
git commit -m "$1"
git pull --rebase origin main
git push origin main

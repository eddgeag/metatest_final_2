#!/usr/bin/env bash
#SBATCH --job-name=metatest_portable
#SBATCH --partition=cpu-dev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
set -Eeuo pipefail
if [[ "${1:-}" == --plan-state ]]; then
  source "$2"
  CONTROL="$(dirname "$2")"
  (cd "$REPO" && sha256sum -c "$CONTROL/code.sha256") || exit 70
  line="$(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "$PLAN")"
  IFS=$'\t' read -r template scenario root <<< "$line"
  [[ -n "$template" && -n "$scenario" && -n "$root" ]] || exit 71
  cd "$REPO"
  exec bash "$REPO/run.sh" "${COMMON[@]}" --plantillas "$REPO/$template" --scenario "$scenario"
else
  # Uso directo: sbatch run.sbatch --backend enroot --scenario ID ...
  # SLURM copia este archivo a spool; por eso usar el directorio de envío.
  REPO="${SLURM_SUBMIT_DIR:?}"
  exec bash "$REPO/run.sh" "$@"
fi

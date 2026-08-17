#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT}"

OUT_DIR="reproducibilidad"
mkdir -p "${OUT_DIR}"

{
  echo "generated_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "git_commit=$(git rev-parse HEAD 2>/dev/null || echo NO_GIT)"
  echo "git_describe=$(git describe --always --dirty --tags 2>/dev/null || echo NO_GIT)"
  echo "kernel=$(uname -srmo)"
  if [[ -r /etc/os-release ]]; then
    sed 's/^/os_/' /etc/os-release
  fi
  echo "R=$(R --version | head -n 1)"
  echo "gcc=$(gcc --version 2>/dev/null | head -n 1 || echo NO_GCC)"
  echo "gxx=$(g++ --version 2>/dev/null | head -n 1 || echo NO_GXX)"
  echo "make=$(make --version 2>/dev/null | head -n 1 || echo NO_MAKE)"
  echo "python=$(python3 --version 2>&1 || echo NO_PYTHON3)"
} > "${OUT_DIR}/entorno_externo_local.txt"

Rscript - <<'RS' > "reproducibilidad/entorno_R_y_cmdstan.txt" 2>&1
cat("R_VERSION=", R.version.string, "\n", sep = "")
cat("PLATFORM=", R.version$platform, "\n", sep = "")
cat("RENVDIR=", Sys.getenv("RENV_PROJECT", unset = getwd()), "\n", sep = "")
if (requireNamespace("renv", quietly = TRUE)) {
  cat("RENV_VERSION=", as.character(utils::packageVersion("renv")), "\n", sep = "")
  print(renv::status())
}
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  cat("CMDSTANR_VERSION=", as.character(utils::packageVersion("cmdstanr")), "\n", sep = "")
  path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
  version <- tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) "")
  cat("CMDSTAN_PATH=", path, "\n", sep = "")
  cat("CMDSTAN_VERSION=", version, "\n", sep = "")
  try(cmdstanr::check_cmdstan_toolchain(fix = FALSE), silent = FALSE)
} else {
  cat("CMDSTANR_NOT_INSTALLED\n")
}
sessionInfo()
RS

Rscript - <<'RS' > "${OUT_DIR}/requisitos_sistema_sugeridos.txt" 2>&1 || true
if (requireNamespace("renv", quietly = TRUE) &&
    "sysreqs" %in% getNamespaceExports("renv")) {
  renv::sysreqs(check = TRUE, report = TRUE, collapse = TRUE)
} else {
  cat("renv::sysreqs() no está disponible en esta versión de renv.\n")
}
RS

find . -maxdepth 2 -type f \
  \( -name '*.R' -o -name '*.sh' -o -name '*.tsv' -o -name 'renv.lock' \) \
  -not -path './comparacion_modelamientos/*' \
  -not -path './renv/library/*' \
  -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > "${OUT_DIR}/SHA256SUMS_codigo.txt"

echo "Manifiesto externo creado en: ${ROOT}/${OUT_DIR}"
echo "Añade reproducibilidad/ a Git y crea un nuevo tag antes de copiar al servidor."


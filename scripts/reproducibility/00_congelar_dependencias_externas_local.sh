#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." >/dev/null 2>&1 && pwd -P)"
WORKSPACE_ROOT="$(cd -- "${REPO_ROOT}/.." >/dev/null 2>&1 && pwd -P)"
DATA_ROOT="${DATA_ROOT:-${WORKSPACE_ROOT}/data}"
RESULTS_ROOT="${RESULTS_ROOT:-${WORKSPACE_ROOT}/results}"
COMPARISON_ROOT="${COMPARISON_ROOT:-${RESULTS_ROOT}/comparacion_modelamientos}"
export REPO_ROOT WORKSPACE_ROOT DATA_ROOT RESULTS_ROOT COMPARISON_ROOT
cd "${REPO_ROOT}"

OUT_DIR="${REPO_ROOT}/reproducibilidad"
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
  echo "make=$(make --version 2>/dev/null | sed -n '1p' || echo NO_MAKE)"
  echo "python=$(python3 --version 2>&1 || echo NO_PYTHON3)"
} > "${OUT_DIR}/entorno_externo_local.txt"

Rscript --vanilla - <<'RS' > "${OUT_DIR}/entorno_R_y_cmdstan.txt" 2>&1
source("renv/activate.R")
cat("R_VERSION=", R.version.string, "\n", sep = "")
cat("PLATFORM=", R.version$platform, "\n", sep = "")
cat("RENVDIR=", Sys.getenv("RENV_PROJECT", unset = getwd()), "\n", sep = "")
if (requireNamespace("renv", quietly = TRUE)) {
  cat("RENV_VERSION=", as.character(utils::packageVersion("renv")), "\n", sep = "")
  renv::status()
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

Rscript --vanilla - <<'RS' > "${OUT_DIR}/entorno_mofa_basilisk.txt" 2>&1
source("renv/activate.R")

pkgs <- c("MOFA2", "basilisk", "basilisk.utils", "reticulate")

for (p in pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(
      toupper(gsub("\\.", "_", p)),
      "_VERSION=",
      as.character(packageVersion(p)),
      "\n",
      sep = ""
    )
  } else {
    cat(
      toupper(gsub("\\.", "_", p)),
      "_VERSION=NO_INSTALADO\n",
      sep = ""
    )
  }
}

if (requireNamespace("MOFA2", quietly = TRUE)) {
  ns <- asNamespace("MOFA2")

  if (exists("mofa_env", envir = ns, inherits = FALSE)) {
    env <- get("mofa_env", envir = ns)

    cat("BASILISK_ENV_CLASS=", class(env)[1], "\n", sep = "")
    cat("BASILISK_ENVNAME=", slot(env, "envname"), "\n", sep = "")
    cat("BASILISK_PKGNAME=", slot(env, "pkgname"), "\n", sep = "")
    cat(
      "BASILISK_PACKAGES=",
      paste(slot(env, "packages"), collapse = "|"),
      "\n",
      sep = ""
    )
    cat(
      "BASILISK_CHANNELS=",
      paste(slot(env, "channels"), collapse = "|"),
      "\n",
      sep = ""
    )
    cat(
      "BASILISK_PIP=",
      paste(slot(env, "pip"), collapse = "|"),
      "\n",
      sep = ""
    )
  } else {
    cat("BASILISK_ENVIRONMENT=NO_DECLARADO\n")
  }
}
RS

Rscript --vanilla - <<'RS' > "${OUT_DIR}/requisitos_sistema_sugeridos.txt" 2>&1 || true
source("renv/activate.R")
if (requireNamespace("renv", quietly = TRUE) &&
    "sysreqs" %in% getNamespaceExports("renv")) {
  renv::sysreqs(check = TRUE, report = TRUE, collapse = TRUE)
} else {
  cat("renv::sysreqs() no está disponible en esta versión de renv.\n")
}
RS

find . -maxdepth 4 -type f \
  \( -name '*.R' -o -name '*.sh' -o -name '*.tsv' -o -name '*.csv' -o -name '*.env' -o -name '*.json' -o -name 'renv.lock' \) \
  -not -path './renv/library/*' \
  -print0 \
  | while IFS= read -r -d '' f; do
      if ! git check-ignore -q -- "$f"; then
        printf '%s\0' "$f"
      fi
    done \
  | sort -z \
  | xargs -0 -r sha256sum \
  > "${OUT_DIR}/SHA256SUMS_codigo.txt"

# Normalizar whitespace de los manifiestos de texto para mantener Git limpio.
for f in \
  "${OUT_DIR}/entorno_R_y_cmdstan.txt" \
  "${OUT_DIR}/entorno_externo_local.txt" \
  "${OUT_DIR}/entorno_mofa_basilisk.txt" \
  "${OUT_DIR}/requisitos_sistema_sugeridos.txt"
do
  if [[ -f "$f" ]]; then
    sed -i 's/[[:space:]]\+$//' "$f"
  fi
done

echo "Manifiesto externo creado en: ${OUT_DIR}"
echo "Revisa reproducibilidad/ antes de versionar o copiar al servidor."

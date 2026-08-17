#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT}"

export RENV_CONFIG_INSTALL_JOBS="${RENV_CONFIG_INSTALL_JOBS:-1}"
export MAKEFLAGS="${MAKEFLAGS:--j1}"
export AUTO_INSTALL_CMDSTAN="${AUTO_INSTALL_CMDSTAN:-1}"

./01_preflight_servidor.sh

LOG_DIR="${ROOT}/.server_runs/restore_$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "${LOG_DIR}"

Rscript -e '
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}
if ("sysreqs" %in% getNamespaceExports("renv")) {
  renv::sysreqs(check = TRUE, report = TRUE, collapse = TRUE)
}
' 2>&1 | tee "${LOG_DIR}/01_requisitos_sistema.txt"

echo
echo "Si el informe anterior muestra paquetes Linux faltantes, instálalos con el administrador"
echo "y vuelve a ejecutar este script. Continuando con restore()..."

Rscript -e '
options(renv.config.install.jobs = 1L)
renv::restore(prompt = FALSE, clean = FALSE, retry = FALSE)
renv::status()
' 2>&1 | tee "${LOG_DIR}/02_renv_restore.txt"

Rscript - <<'RS' 2>&1 | tee "${LOG_DIR}/03_validacion_paquetes_cmdstan.txt"
critical <- c("brms", "cmdstanr", "MOFA2", "mixOmics", "caret", "clusterProfiler")
missing <- critical[!vapply(critical, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Faltan paquetes críticos: ", paste(missing, collapse = ", "))

cat("Paquetes críticos OK\n")
cmdstanr::check_cmdstan_toolchain(fix = FALSE)

expected_file <- "reproducibilidad/entorno_R_y_cmdstan.txt"
expected <- ""
if (file.exists(expected_file)) {
  lines <- readLines(expected_file, warn = FALSE)
  hit <- grep("^CMDSTAN_VERSION=", lines, value = TRUE)
  if (length(hit)) expected <- sub("^CMDSTAN_VERSION=", "", hit[[1]])
}

installed <- tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) "")
needs_install <- !nzchar(installed) || (nzchar(expected) && expected != installed)
if (needs_install && identical(Sys.getenv("AUTO_INSTALL_CMDSTAN", "1"), "1")) {
  if (!nzchar(expected)) {
    stop("No se registró CMDSTAN_VERSION en el manifiesto local; no se instalará una versión arbitraria.")
  }
  cat("Instalando CmdStan ", expected, " de forma serial...\n", sep = "")
  cmdstanr::install_cmdstan(version = expected, cores = 1, overwrite = FALSE, quiet = FALSE)
  installed <- tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) "")
}
if (!nzchar(installed)) {
  stop("CmdStan no está disponible después de la restauración.")
}
if (nzchar(expected) && expected != installed) {
  stop("Versión CmdStan distinta después de restaurar: esperada=", expected, "; activa=", installed)
}
cat("CmdStan version: ", installed, "\n", sep = "")
cat("CmdStan path: ", cmdstanr::cmdstan_path(), "\n", sep = "")
cat("ENTORNO R Y CMDSTAN OK\n")
RS

echo "ENTORNO RESTAURADO Y VALIDADO"
echo "Logs: ${LOG_DIR}"

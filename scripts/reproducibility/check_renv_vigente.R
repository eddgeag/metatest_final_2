#!/usr/bin/env Rscript
# Fail-fast runtime/package gate: this does not train or change packages.
required <- c("renv", "MOFA2", "brms", "cmdstanr", "posterior", "msigdbr",
              "msgl", "sglOptim", "glmnet", "mixOmics", "igraph", "ggplot2",
              "ggrepel", "dplyr", "tidyr", "patchwork", "ragg")
expected_r <- "4.4.2"
if (as.character(getRversion()) != expected_r) {
  stop("Requiere R ", expected_r, "; versión actual: ", as.character(getRversion()))
}
if (!requireNamespace("renv", quietly = TRUE)) stop("Instalar renv y restaurar renv.lock")
root <- Sys.getenv("METATEST_REPO_ROOT", unset = getwd())
lock_path <- file.path(root, "renv.lock")
if (!file.exists(lock_path)) stop("Falta renv.lock")
lock <- renv::lockfile_read(lock_path)
missing_in_lock <- setdiff(required, names(lock$Packages))
if (length(missing_in_lock)) {
  stop("Actualizar renv.lock antes de ejecutar (faltan): ", paste(missing_in_lock, collapse = ", "),
       "\nVer scripts/reproducibility/actualizar_renv_vigente.R")
}
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Paquetes no instalados: ", paste(missing, collapse = ", "))
not_matching <- required[vapply(required, function(p) {
  utils::packageDescription(p, fields = "Version") != lock$Packages[[p]]$Version
}, logical(1))]
if (length(not_matching)) stop("Versiones distintas a renv.lock: ", paste(not_matching, collapse = ", "))
message("R y paquetes clave: OK (", length(required), " paquetes)")
if (!nzchar(Sys.getenv("CMDSTAN", unset = "")) &&
    !dir.exists(path.expand("~/.cmdstan/cmdstan-2.36.0"))) {
  warning("No se ha detectado CMDSTAN; comprobar antes de ejecutar 02")
}

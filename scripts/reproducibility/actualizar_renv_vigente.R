#!/usr/bin/env Rscript
# EJECUTAR UNA VEZ EN LOCAL; revisar el diff de renv.lock antes del commit.
# Instala paquetes ausentes sin hacer renv::snapshot() global silenciosamente.
stopifnot(as.character(getRversion()) == "4.4.2")
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
root <- normalizePath(getwd())
if (!file.exists(file.path(root, "renv.lock"))) stop("Ejecutar desde la raíz de metatest")
renv::restore(project = root, prompt = FALSE)
renv::load(project = root)
# Los dos paquetes históricos de sparse group lasso se instalan desde el
# archivo de CRAN si no están presentes en la biblioteca del proyecto.
# Es obligatorio congelar su procedencia en el lockfile final.
if (!requireNamespace("sglOptim", quietly = TRUE)) {
  renv::install("sglOptim@1.3.8", prompt = FALSE)
}
if (!requireNamespace("msgl", quietly = TRUE)) {
  renv::install("msgl@2.3.9", prompt = FALSE)
}
if (!requireNamespace("msigdbr", quietly = TRUE) ||
    as.character(utils::packageVersion("msigdbr")) != "26.1.1") {
  renv::install("msigdbr@26.1.1", prompt = FALSE)
}
for (p in c("sglOptim", "msgl", "msigdbr")) {
  if (!requireNamespace(p, quietly=TRUE)) stop("Falta paquete: ", p)
}
# Como el .gitignore es restrictivo y se usan scripts autónomos,
# registrar explícitamente estos paquetes en el snapshot.
# El snapshot del conjunto original debe mantenerse también.
renv::snapshot(project = root, type = "all", prompt = FALSE)
message("LISTO: inspeccionar git diff -- renv.lock; ejecutar preflight y luego commit")

# Verificacion del runtime nativo, sin restaurar paquetes ni entrenar modelos.
fail <- function(...) stop(..., call. = FALSE)
must <- function(condition, ...) if (!isTRUE(condition)) fail(...)
norm <- function(path) normalizePath(path, mustWork = TRUE)
sha <- function(path) {
  must(file.exists(path), "No existe para calcular SHA256: ", path)
  value <- system2("sha256sum", shQuote(path), stdout = TRUE)
  must(is.null(attr(value, "status")), "Fallo sha256sum: ", path)
  sub(" .*", "", value[[1]])
}

verify_library <- function(repo, lib) {
  must(identical(as.character(getRversion()), "4.4.2"), "Se requiere R 4.4.2")
  # Comprobar la herencia real; no corregir .libPaths() dentro de esta prueba.
  must(identical(norm(.libPaths()), unique(norm(c(lib, .Library)))),
       "R no heredo exclusivamente la biblioteca elegida y la biblioteca base")
  lock <- jsonlite::fromJSON(file.path(repo, "renv.lock"), simplifyVector = FALSE)
  must(identical(lock$R$Version, "4.4.2") && length(lock$Packages) == 302L,
       "Lockfile ajeno al perfil nativo auditado")
  ip <- installed.packages(lib.loc = lib, fields = c("RemoteSha", "Built"), noCache = TRUE)
  expected <- vapply(lock$Packages, function(x) x$Version, character(1))
  observed <- ip[match(names(expected), ip[, "Package"]), "Version"]
  bad <- is.na(observed) | observed != expected
  if (any(bad)) fail("Versiones faltantes/distintas: ", paste(names(expected)[bad], collapse = ", "))
  for (p in names(lock$Packages)) {
    remote <- lock$Packages[[p]]$RemoteSha
    if (!is.null(remote)) {
      must(identical(unname(ip[p, "RemoteSha"]), remote), "RemoteSha distinto: ", p)
    }
  }
  ip[order(ip[, "Package"]), c("Package", "Version", "Built", "RemoteSha"), drop = FALSE]
}

resolve_backend <- function() {
  original <- path.expand(basilisk.utils::getExternalDir())
  cache <- file.path(norm(dirname(original)), basename(original))
  must(dir.exists(cache), "No existe la cache Basilisk: ", cache)
  root <- dirname(cache)
  conda <- norm(basilisk.utils::getCondaDir())
  # No normalizePath() del Python: no perder el prefijo del entorno por un symlink.
  python <- file.path(cache, "MOFA2", "1.16.0", "mofa_env", "bin", "python")
  must(file.access(python, 1L) == 0L, "No existe Python MOFA: ", python,
       ". No se instalara automaticamente")
  must(file.access(file.path(conda, "bin", "conda"), 1L) == 0L,
       "No existe Conda base: ", conda, ". No se instalara automaticamente")
  must(file.access(cache, 2L) == 0L && file.access(dirname(dirname(python)), 2L) == 0L,
       "Basilisk necesita permisos de escritura para sus locks en la cache")
  result <- c(BASILISK_EXTERNAL_DIR = root, BASILISK_EXTERNAL_CONDA = conda,
              RETICULATE_PYTHON = python)
  must(!any(grepl("[\t\r\n:]", result)), "Ruta de backend no representable")
  result
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  must(length(args) >= 1L && args[[1]] %in% c("resolve", "check"), "Modo nativo invalido")
  repo <- norm(Sys.getenv("METATEST_REPO_ROOT"))
  lib <- norm(Sys.getenv("METATEST_R_LIBRARY"))
  inventory <- verify_library(repo, lib)
  backend <- resolve_backend()
  if (identical(args[[1]], "resolve")) {
    cat(paste(names(backend), backend, sep = "\t"), sep = "\n")
    return(invisible(NULL))
  }
  must(length(args) == 2L, "Falta directorio de informe")
  report <- norm(args[[2]])
  must(identical(norm(Sys.getenv("RETICULATE_PYTHON")), norm(backend[["RETICULATE_PYTHON"]])),
       "Python no coincide con el entorno Basilisk esperado")
  for (p in c("renv", "ggplot2", "cmdstanr", "brms", "MOFA2", "basilisk",
              "basilisk.utils", "reticulate", "igraph", "readxl", "openxlsx")) {
    ns <- suppressPackageStartupMessages(loadNamespace(p, lib.loc = lib))
    must(identical(norm(getNamespaceInfo(ns, "path")), norm(file.path(lib, p))),
         "Paquete cargado fuera de la biblioteca elegida: ", p)
    cat("CARGA_OK:", p, as.character(packageVersion(p)), "\n")
  }
  must(igraph::vcount(igraph::make_ring(3L)) == 3L, "Prueba basica igraph fallida")
  cat("IGRAPH_LOAD=PASS\nR_LOCK_VERSIONES=PASS\nR_LIBRARY_INHERITANCE=PASS\n")
  cmdstan <- norm(Sys.getenv("CMDSTAN"))
  must(identical(norm(cmdstanr::cmdstan_path()), cmdstan), "Ruta CmdStan distinta")
  must(identical(as.character(cmdstanr::cmdstan_version()), "2.36.0"), "CmdStan debe ser 2.36.0")
  config <- reticulate::py_config()
  must(identical(norm(config$python), norm(backend[["RETICULATE_PYTHON"]])), "Reticulate uso otro Python")
  reticulate::py_run_string(paste(
    "import sys, json, importlib.metadata as metadata, mofapy2",
    "from mofapy2.run.entry_point import entry_point",
    "assert mofapy2.__version__ == '0.7.0'",
    "metatest_native_python = json.dumps({'version': sys.version, 'executable': sys.executable, 'prefix': sys.prefix, 'packages': sorted([(d.metadata.get('Name', ''), d.version) for d in metadata.distributions()])}, sort_keys=True)",
    sep = "\n"
  ))
  python_identity <- jsonlite::fromJSON(reticulate::py_eval("metatest_native_python"), simplifyVector = FALSE)
  cat("RETICULATE_PYTHON:", config$python, "\nRETICULATE_MOFAPY2_IMPORT=PASS\n")

  # Si la resolucion de Basilisk pidiera crear un entorno, abortar antes de hacerlo.
  suppressMessages(trace("setupBasiliskEnv", where = asNamespace("basilisk"),
    tracer = quote(stop("La prueba nativa no permite crear entornos", call. = FALSE)), print = FALSE))
  on.exit(suppressMessages(untrace("setupBasiliskEnv", where = asNamespace("basilisk"))), add = TRUE)
  env <- get("mofa_env", envir = asNamespace("MOFA2"))
  expected_env <- dirname(dirname(backend[["RETICULATE_PYTHON"]]))
  must(identical(norm(basilisk::obtainEnvironmentPath(env)), norm(expected_env)), "Basilisk resolvio otro entorno")
  proc <- basilisk::basiliskStart(env)
  on.exit(basilisk::basiliskStop(proc), add = TRUE, after = FALSE)
  observed <- basilisk::basiliskRun(proc, function() {
    mod <- reticulate::import("mofapy2", convert = FALSE)
    list(version = reticulate::py_to_r(mod$version$`__version__`),
         python = reticulate::py_config()$python,
         mofa = as.character(utils::packageVersion("MOFA2")))
  })
  must(identical(observed$version, "0.7.0") && identical(observed$mofa, "1.16.0") &&
       identical(norm(observed$python), norm(config$python)), "Basilisk uso otro runtime")
  cat("BASILISK_MOFAPY2_BRIDGE=PASS\n")

  # Huella de metadatos + binarios criticos, no una certificacion bit a bit del SO.
  binaries <- c(Rscript = Sys.getenv("METATEST_RSCRIPT"),
                CXX = unname(Sys.which("g++")),
                stanc = file.path(cmdstan, "bin", "stanc"),
                cmdstan_makefile = file.path(cmdstan, "makefile"),
                python = backend[["RETICULATE_PYTHON"]])
  if (file.exists(file.path(R.home("lib"), "libR.so")))
    binaries <- c(binaries, libR = file.path(R.home("lib"), "libR.so"))
  if (file.exists(file.path(cmdstan, "make", "local")))
    binaries <- c(binaries, cmdstan_local = file.path(cmdstan, "make", "local"))
  if (!is.null(config$libpython) && nzchar(config$libpython) && file.exists(config$libpython))
    binaries <- c(binaries, libpython = config$libpython)
  shared <- unlist(lapply(inventory[, "Package"], function(p) {
    paths <- list.files(file.path(lib, p, "libs"), pattern = "\\.so$", full.names = TRUE, recursive = TRUE)
    if (length(paths)) setNames(paths, paste0("R/", p, "/", substring(paths, nchar(file.path(lib, p, "libs")) + 2L))) else character()
  }), use.names = TRUE)
  binaries <- c(binaries, shared)
  hashes <- vapply(binaries, sha, character(1))
  identity <- list(schema = 1L, R = R.version.string, platform = R.version$platform,
    os = as.list(Sys.info()[c("sysname", "release", "machine")]),
    library = lib, R_packages = lapply(seq_len(nrow(inventory)), function(i) as.list(inventory[i, ])),
    cmdstan = cmdstan, cmdstan_version = "2.36.0", stan_optimization = 0L,
    startup = "PROJECT_LIBRARY_AND_BASE_NO_USER_PROFILE", python = python_identity,
    basilisk = as.list(backend), binaries_sha256 = as.list(hashes))
  jsonlite::write_json(identity, file.path(report, "native_identity.json"),
                       auto_unbox = TRUE, pretty = TRUE, na = "null")
  writeLines(capture.output(sessionInfo()), file.path(report, "sessionInfo.txt"))
  cat("NATIVE_DEPENDENCIES=PASS\nNO_INSTALLATIONS_REQUESTED=TRUE\n")
}

main()

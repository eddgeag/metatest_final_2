launcher_args <- base::commandArgs(trailingOnly = TRUE)

if (length(launcher_args) != 4L) {
  stop(
    "Uso interno: launcher FRAMEWORK ROOT DIAG MODE"
  )
}

framework_file <- normalizePath(
  path.expand(launcher_args[[1]]),
  mustWork = TRUE
)

framework_args <- launcher_args[2:4]

requested_root <- normalizePath(
  path.expand(framework_args[[1]]),
  mustWork = TRUE
)

exact_scenario <- normalizePath(
  path.expand(Sys.getenv("BIDIR_EXACT_SCENARIO", unset = "")),
  mustWork = TRUE
)

if (!identical(dirname(exact_scenario), requested_root)) {
  stop(
    "El escenario exacto no pertenece a la raíz solicitada:\n",
    exact_scenario,
    "\nRaíz:\n",
    requested_root
  )
}

base_command_args <- base::commandArgs
base_list_dirs <- base::list.dirs

commandArgs <- function(trailingOnly = FALSE) {
  if (isTRUE(trailingOnly)) {
    return(framework_args)
  }

  base_command_args(trailingOnly = FALSE)
}

list.dirs <- function(path = ".", full.names = TRUE, recursive = TRUE) {
  normalized_path <- normalizePath(
    path.expand(path),
    mustWork = TRUE
  )

  if (identical(normalized_path, requested_root) && identical(recursive, FALSE)) {
    if (isTRUE(full.names)) {
      return(exact_scenario)
    }

    return(basename(exact_scenario))
  }

  base_list_dirs(
    path = path,
    full.names = full.names,
    recursive = recursive
  )
}

sys.source(framework_file, envir = .GlobalEnv)

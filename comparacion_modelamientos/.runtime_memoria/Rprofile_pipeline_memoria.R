.user_profile <- path.expand("~/.Rprofile")
if (file.exists(.user_profile)) {
  try(sys.source(.user_profile, envir = .GlobalEnv), silent = TRUE)
}

.pipeline_cores <- suppressWarnings(
  as.integer(Sys.getenv("PIPELINE_R_CORES", unset = "2"))
)
.brms_internal_cores <- suppressWarnings(
  as.integer(Sys.getenv("BRMS_INTERNAL_CORES", unset = "1"))
)

if (is.na(.pipeline_cores) || .pipeline_cores < 1L) {
  .pipeline_cores <- 2L
}
if (is.na(.brms_internal_cores) || .brms_internal_cores < 1L) {
  .brms_internal_cores <- 1L
}

options(mc.cores = .brms_internal_cores)

## Los scripts calculan sus workers a partir de parallel::detectCores().
## Se limita únicamente durante esta ejecución, sin editar los scripts R.
.parallel_ns <- asNamespace("parallel")
try(unlockBinding("detectCores", .parallel_ns), silent = TRUE)
assign(
  "detectCores",
  local({
    n <- .pipeline_cores
    function(all.tests = FALSE, logical = TRUE) n
  }),
  envir = .parallel_ns
)
try(lockBinding("detectCores", .parallel_ns), silent = TRUE)

rm(
  .user_profile,
  .pipeline_cores,
  .brms_internal_cores,
  .parallel_ns
)

#!/usr/bin/env Rscript

## ============================================================================
## SELECCION DE LOS SEIS MEJORES ESPACIOS LATENTES
##
## Competiciones independientes:
##   TODAS     x modelamiento2 / modelamiento3
##   CON_FLI   x modelamiento2 / modelamiento3
##   SIN_FLI   x modelamiento2 / modelamiento3
##
## Principios:
##   * La seleccion usa exclusivamente TRAIN.
##   * El MOFA original de cada escenario se evalua con todas las variables que
##     sobrevivieron a su prefiltrado.
##   * Las features Bayes son exactamente las usadas por
##     diagnostico_bayes_features.
##   * Las vistas se combinan con peso igual mediante kernels normalizados.
##   * Rango insuficiente y mala condicion son penalizaciones continuas; no son
##     filtros duros.
##   * Solo faltantes, no finitos o imposibilidad tecnica de construir la
##     representacion invalidan un candidato.
##   * q00 compite y se informa como referencia de sensibilidad, pero la
##     similitud consigo mismo no participa en el ranking.
##   * El numero de features no es un objetivo. Solo constituye el ultimo
##     desempate determinista y, por decision del proyecto, favorece primero
##     el panel util/mapeado mayor y despues el panel listado mayor si todo lo
##     demas resulta exactamente igual.
##
## Este script NO reentrena MOFA ni los modelos Bayes.
## ============================================================================

options(stringsAsFactors = FALSE, width = 220)

SCRIPT_VERSION <- "1.0.0"
CANONICAL_VIEWS <- c("transcriptomica", "proteomica", "metabolomica", "clinical")
VIEW_CODE <- c(
  transcriptomica = "tx",
  proteomica = "pr",
  metabolomica = "me",
  clinical = "cl"
)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

usage <- function(status = 0L) {
  cat(
    "Uso:\n",
    "  Rscript seleccionar_mejores_espacios_latentes.R [opciones]\n\n",
    "Opciones:\n",
    "  --project-root RUTA       Raiz del proyecto (default: .)\n",
    "  --comparison-root RUTA    Raiz de comparacion (default: comparacion_modelamientos)\n",
    "  --output-dir RUTA         Carpeta nueva de salida (default: timestamp)\n",
    "  --families LISTA          TODAS,CON_FLI,SIN_FLI (default: las tres)\n",
    "  --models LISTA            modelamiento2,modelamiento3 (default: ambos)\n",
    "  --ridge-grid LISTA        Regularizacion relativa SVD\n",
    "                            (default: 1e-8,1e-6,1e-4)\n",
    "  --no-test                 No proyectar TEST tras congelar ganadores\n",
    "  --no-copy-models          No copiar los RDS MOFA ganadores a la salida\n",
    "  --help                    Mostrar ayuda\n\n",
    "Ejemplo:\n",
    "  Rscript seleccionar_mejores_espacios_latentes.R --project-root .\n",
    sep = ""
  )
  quit(save = "no", status = status)
}

parse_csv_arg <- function(x) {
  trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
}

parse_args <- function(args) {
  out <- list(
    project_root = ".",
    comparison_root = "comparacion_modelamientos",
    output_dir = NULL,
    families = c("TODAS", "CON_FLI", "SIN_FLI"),
    models = c("modelamiento2", "modelamiento3"),
    ridge_grid = c(1e-8, 1e-6, 1e-4),
    do_test = TRUE,
    copy_models = TRUE
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--help") usage(0L)
    if (key == "--no-test") {
      out$do_test <- FALSE
      i <- i + 1L
      next
    }
    if (key == "--no-copy-models") {
      out$copy_models <- FALSE
      i <- i + 1L
      next
    }
    if (!key %in% c(
      "--project-root", "--comparison-root", "--output-dir",
      "--families", "--models", "--ridge-grid"
    )) {
      stop("Argumento no reconocido: ", key)
    }
    if (i == length(args)) stop("Falta valor para: ", key)
    value <- args[[i + 1L]]
    if (key == "--project-root") out$project_root <- value
    if (key == "--comparison-root") out$comparison_root <- value
    if (key == "--output-dir") out$output_dir <- value
    if (key == "--families") out$families <- toupper(parse_csv_arg(value))
    if (key == "--models") out$models <- tolower(parse_csv_arg(value))
    if (key == "--ridge-grid") {
      out$ridge_grid <- suppressWarnings(as.numeric(parse_csv_arg(value)))
    }
    i <- i + 2L
  }

  allowed_families <- c("TODAS", "CON_FLI", "SIN_FLI")
  allowed_models <- c("modelamiento2", "modelamiento3")
  if (!length(out$families) || any(!out$families %in% allowed_families)) {
    stop("--families admite: ", paste(allowed_families, collapse = ","))
  }
  if (!length(out$models) || any(!out$models %in% allowed_models)) {
    stop("--models admite: ", paste(allowed_models, collapse = ","))
  }
  if (!length(out$ridge_grid) || any(!is.finite(out$ridge_grid)) || any(out$ridge_grid <= 0)) {
    stop("--ridge-grid debe contener numeros positivos y finitos.")
  }
  out$ridge_grid <- sort(unique(out$ridge_grid))
  out
}

normalize_existing <- function(path, label) {
  if (!file.exists(path) && !dir.exists(path)) stop("No existe ", label, ": ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

normalize_maybe <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

path_inside <- function(path, root) {
  p <- normalize_maybe(path)
  r <- sub("/+$", "", normalize_maybe(root))
  identical(p, r) || startsWith(p, paste0(r, "/"))
}

relative_to <- function(path, root) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  p <- normalize_maybe(path)
  r <- sub("/+$", "", normalize_maybe(root))
  if (!path_inside(p, r)) return(p)
  if (identical(p, r)) return(".")
  substring(p, nchar(r) + 2L)
}

actual_child_dir <- function(parent, requested) {
  exact <- file.path(parent, requested)
  if (dir.exists(exact)) return(normalizePath(exact, winslash = "/"))
  if (!dir.exists(parent)) return(NA_character_)
  children <- list.dirs(parent, recursive = FALSE, full.names = TRUE)
  hit <- children[tolower(basename(children)) == tolower(requested)]
  if (!length(hit)) return(NA_character_)
  normalizePath(sort(hit)[[1L]], winslash = "/")
}

first_existing <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  hit <- paths[file.exists(paths)]
  if (length(hit)) normalizePath(hit[[1L]], winslash = "/") else NA_character_
}

find_one <- function(dir, pattern) {
  if (!dir.exists(dir)) return(NA_character_)
  x <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = FALSE)
  x <- sort(x[file.exists(x)])
  if (length(x)) normalizePath(x[[1L]], winslash = "/") else NA_character_
}

safe_read_rds <- function(path) {
  tryCatch(
    list(ok = TRUE, value = readRDS(path), error = NA_character_),
    error = function(e) list(ok = FALSE, value = NULL, error = conditionMessage(e))
  )
}

safe_read_csv <- function(path) {
  tryCatch(
    list(
      ok = TRUE,
      value = read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
      error = NA_character_
    ),
    error = function(e) list(ok = FALSE, value = NULL, error = conditionMessage(e))
  )
}

write_tsv <- function(x, path) {
  if (is.null(x) || !ncol(as.data.frame(x))) {
    writeLines(character(), path, useBytes = TRUE)
    return(invisible(path))
  }
  write.table(
    x, path, sep = "\t", row.names = FALSE, col.names = TRUE,
    quote = TRUE, na = "NA", fileEncoding = "UTF-8"
  )
}

bind_rows_fill <- function(xs) {
  xs <- xs[!vapply(xs, is.null, logical(1))]
  if (!length(xs)) return(data.frame())
  cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    miss <- setdiff(cols, names(x))
    for (nm in miss) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, xs)
}

canonical_view <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))
  out[y %in% c("tx", "transcriptomica", "transcriptomics", "transcriptomic")] <- "transcriptomica"
  out[y %in% c("pr", "proteomica", "proteomics", "proteomic")] <- "proteomica"
  out[y %in% c("me", "metabolomica", "metabolomics", "metabolomic")] <- "metabolomica"
  out[y %in% c("cl", "clinical", "clinica", "clinico", "clinicos")] <- "clinical"
  out
}

canonical_prefix_ok <- function(root_name, analysis) {
  if (root_name == "escenarios") return(grepl("^P1_(ALL|TODAS)_", analysis))
  if (root_name == "escenarios_con_FLI") return(grepl("^P2_CON_FLI_", analysis))
  if (root_name == "escenarios_SIN_FLI") return(grepl("^P3_SIN_FLI_", analysis))
  FALSE
}

family_specifications <- function() {
  list(
    TODAS = list(root = "escenarios", prefix = "P1_"),
    CON_FLI = list(root = "escenarios_con_FLI", prefix = "P2_"),
    SIN_FLI = list(root = "escenarios_SIN_FLI", prefix = "P3_")
  )
}

is_q00 <- function(analysis) {
  grepl("_qtx00_qpr00_qme00_qcl00$", analysis)
}

scenario_names <- function(root) {
  if (is.na(root) || !dir.exists(root)) return(character())
  sort(basename(list.dirs(root, recursive = FALSE, full.names = TRUE)))
}

mofa_file_for <- function(shared_dir, analysis) {
  first_existing(c(
    file.path(shared_dir, "integracion", "rds", paste0("modelo_mofa_final_", analysis, ".rds")),
    find_one(file.path(shared_dir, "integracion", "rds"), "^modelo_mofa_final_.*\\.rds$")
  ))
}

ready_file_for <- function(shared_dir, analysis) {
  first_existing(c(
    file.path(shared_dir, paste0("ready_for_modeling_", analysis, ".rds")),
    find_one(shared_dir, "^ready_for_modeling_.*\\.rds$")
  ))
}

panel_csv_for <- function(branch_dir) {
  file.path(
    branch_dir, "diagnostico_bayes_features", "tables",
    "features_usadas_para_diagnostico_bayes_input.csv"
  )
}

diagnostic_objects_for <- function(branch_dir) {
  file.path(branch_dir, "diagnostico_bayes_features", "rds", "diagnostic_objects.rds")
}

build_registry <- function(comparison_root, families, models) {
  specs <- family_specifications()
  rows <- list()
  i <- 0L

  for (family in families) {
    spec <- specs[[family]]
    upstream_root <- actual_child_dir(file.path(comparison_root, "upstream_comun"), spec$root)
    upstream_names <- scenario_names(upstream_root)
    if (!length(upstream_names)) {
      warning("No se detectaron escenarios upstream para ", family, ".")
      next
    }

    for (analysis in upstream_names) {
      shared_dir <- file.path(upstream_root, analysis)
      for (model_variant in models) {
        model_root <- actual_child_dir(file.path(comparison_root, model_variant), spec$root)
        branch_dir <- if (!is.na(model_root)) file.path(model_root, analysis) else NA_character_
        i <- i + 1L
        rows[[i]] <- data.frame(
          family = family,
          family_root = spec$root,
          modelamiento = model_variant,
          analysis = analysis,
          q00 = is_q00(analysis),
          canonical_prefix_ok = canonical_prefix_ok(spec$root, analysis),
          shared_dir = shared_dir,
          branch_dir = branch_dir,
          shared_exists = dir.exists(shared_dir),
          branch_exists = !is.na(branch_dir) && dir.exists(branch_dir),
          ready_file = ready_file_for(shared_dir, analysis),
          mofa_file = mofa_file_for(shared_dir, analysis),
          panel_csv = if (!is.na(branch_dir)) panel_csv_for(branch_dir) else NA_character_,
          diagnostic_objects = if (!is.na(branch_dir)) diagnostic_objects_for(branch_dir) else NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  registry <- bind_rows_fill(rows)
  if (!nrow(registry)) stop("No se construyo ningun candidato.")
  registry
}

assert_numeric_matrix <- function(x, label) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (length(dim(x)) != 2L || any(dim(x) == 0L)) stop(label, " esta vacia.")
  if (any(!is.finite(x))) stop(label, " contiene NA/NaN/Inf.")
  x
}

effective_number <- function(x) {
  x <- pmax(as.numeric(x), 0)
  x <- x[is.finite(x)]
  if (!length(x) || sum(x) <= 0) return(0)
  p <- x / sum(x)
  p <- p[p > 0]
  exp(-sum(p * log(p)))
}

svd_diagnostics <- function(x, expected_rank = min(dim(x))) {
  x <- as.matrix(x)
  if (!length(x) || any(!is.finite(x))) {
    return(list(rank = 0L, condition = Inf, min_singular = 0, singular_values = numeric()))
  }
  d <- tryCatch(svd(x, nu = 0L, nv = 0L)$d, error = function(e) numeric())
  if (!length(d) || !all(is.finite(d)) || max(d) <= 0) {
    return(list(rank = 0L, condition = Inf, min_singular = 0, singular_values = d))
  }
  tol <- max(dim(x)) * max(d) * .Machine$double.eps
  rank_eff <- sum(d > tol)
  cond <- if (rank_eff < expected_rank || rank_eff == 0L) Inf else max(d) / min(d[d > tol])
  list(
    rank = as.integer(rank_eff),
    condition = cond,
    min_singular = if (rank_eff) min(d[d > tol]) else 0,
    singular_values = d
  )
}

numeric_stability <- function(rank, expected_rank, condition) {
  if (!is.finite(rank) || expected_rank <= 0) return(0)
  rank_fraction <- min(max(rank / expected_rank, 0), 1)
  if (!is.finite(condition) || condition <= 0) return(0)
  rank_fraction / (1 + log10(max(condition, 1)))
}

center_kernel <- function(K) {
  K <- assert_numeric_matrix(K, "Kernel")
  if (nrow(K) != ncol(K)) stop("El kernel no es cuadrado.")
  n <- nrow(K)
  H <- diag(n) - matrix(1 / n, nrow = n, ncol = n)
  out <- H %*% K %*% H
  dimnames(out) <- dimnames(K)
  out
}

row_standardize <- function(X) {
  X <- assert_numeric_matrix(X, "Matriz de features")
  ctr <- rowMeans(X)
  Xc <- X - ctr
  s <- apply(Xc, 1L, stats::sd)
  keep <- is.finite(s) & s > sqrt(.Machine$double.eps)
  if (!any(keep)) return(matrix(0, nrow = 0L, ncol = ncol(X), dimnames = list(NULL, colnames(X))))
  Xc <- Xc[keep, , drop = FALSE]
  sweep(Xc, 1L, s[keep], "/")
}

kernel_from_features <- function(X, standardize_rows = TRUE) {
  X <- assert_numeric_matrix(X, "Matriz para kernel")
  if (standardize_rows) {
    Xc <- row_standardize(X)
  } else {
    Xc <- X - rowMeans(X)
  }
  if (!nrow(Xc) || sum(Xc * Xc) <= .Machine$double.eps) {
    K <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    dimnames(K) <- list(colnames(X), colnames(X))
    return(K)
  }
  K <- crossprod(Xc)
  center_kernel(K)
}

kernel_from_latent <- function(Z, balance_factors = FALSE) {
  Z <- assert_numeric_matrix(Z, "Espacio latente")
  Zc <- scale(Z, center = TRUE, scale = FALSE)
  Zc <- as.matrix(Zc)
  if (balance_factors) {
    s <- apply(Zc, 2L, stats::sd)
    keep <- is.finite(s) & s > sqrt(.Machine$double.eps)
    if (!any(keep)) {
      K <- matrix(0, nrow = nrow(Z), ncol = nrow(Z), dimnames = list(rownames(Z), rownames(Z)))
      return(K)
    }
    Zc <- sweep(Zc[, keep, drop = FALSE], 2L, s[keep], "/")
  }
  K <- tcrossprod(Zc)
  dimnames(K) <- list(rownames(Z), rownames(Z))
  center_kernel(K)
}

normalize_kernel <- function(K) {
  K <- center_kernel(K)
  nrm <- sqrt(sum(K * K))
  if (!is.finite(nrm) || nrm <= .Machine$double.eps) return(K * 0)
  K / nrm
}

combine_view_kernels <- function(kernels, expected_views = CANONICAL_VIEWS) {
  if (!all(expected_views %in% names(kernels))) {
    stop("Faltan kernels de vistas: ", paste(setdiff(expected_views, names(kernels)), collapse = ", "))
  }
  kernels <- kernels[expected_views]
  dims <- lapply(kernels, dim)
  if (!all(vapply(dims[-1L], identical, logical(1), dims[[1L]]))) {
    stop("Los kernels de vistas no tienen dimensiones identicas.")
  }
  normalized <- lapply(kernels, normalize_kernel)
  out <- Reduce("+", normalized) / length(normalized)
  dimnames(out) <- dimnames(normalized[[1L]])
  center_kernel(out)
}

linear_cka <- function(K, L) {
  K <- center_kernel(K)
  L <- center_kernel(L)
  if (!identical(dim(K), dim(L))) stop("CKA requiere kernels de igual dimension.")
  den <- sqrt(sum(K * K) * sum(L * L))
  if (!is.finite(den) || den <= .Machine$double.eps) return(0)
  val <- sum(K * L) / den
  min(max(as.numeric(val), 0), 1)
}

kernel_distances <- function(K) {
  K <- center_kernel(K)
  d2 <- outer(diag(K), diag(K), "+") - 2 * K
  d2[d2 < 0 & d2 > -1e-10] <- 0
  d2 <- pmax(d2, 0)
  sqrt(d2)
}

kernel_distance_correlation <- function(K, L) {
  if (!identical(dim(K), dim(L))) return(NA_real_)
  d1 <- kernel_distances(K)
  d2 <- kernel_distances(L)
  idx <- upper.tri(d1)
  x <- d1[idx]
  y <- d2[idx]
  if (length(x) < 3L || stats::sd(x) == 0 || stats::sd(y) == 0) return(0)
  as.numeric(stats::cor(x, y, method = "spearman", use = "complete.obs"))
}

loo_cka_q10 <- function(K, L) {
  n <- nrow(K)
  if (n < 6L || !identical(dim(K), dim(L))) return(NA_real_)
  vals <- vapply(seq_len(n), function(i) {
    keep <- setdiff(seq_len(n), i)
    linear_cka(K[keep, keep, drop = FALSE], L[keep, keep, drop = FALSE])
  }, numeric(1))
  as.numeric(stats::quantile(vals, probs = 0.10, na.rm = TRUE, names = FALSE, type = 8))
}

procrustes_similarity <- function(A, B) {
  A <- assert_numeric_matrix(A, "A Procrustes")
  B <- assert_numeric_matrix(B, "B Procrustes")
  if (nrow(A) != nrow(B)) return(NA_real_)
  k <- min(ncol(A), ncol(B))
  if (k < 1L) return(0)
  A <- scale(A[, seq_len(k), drop = FALSE], center = TRUE, scale = TRUE)
  B <- scale(B[, seq_len(k), drop = FALSE], center = TRUE, scale = TRUE)
  A[!is.finite(A)] <- 0
  B[!is.finite(B)] <- 0
  na <- sqrt(sum(A * A))
  nb <- sqrt(sum(B * B))
  if (na <= 0 || nb <= 0) return(0)
  A <- A / na
  B <- B / nb
  s <- svd(crossprod(B, A), nu = k, nv = k)
  min(max(sum(s$d), 0), 1)
}

factor_recovery_r2 <- function(Z, Zhat) {
  Z <- assert_numeric_matrix(Z, "Z referencia")
  Zhat <- assert_numeric_matrix(Zhat, "Z recuperado")
  k <- min(ncol(Z), ncol(Zhat))
  out <- rep(0, k)
  for (j in seq_len(k)) {
    if (stats::sd(Z[, j]) > 0 && stats::sd(Zhat[, j]) > 0) {
      r <- suppressWarnings(stats::cor(Z[, j], Zhat[, j], use = "complete.obs"))
      out[[j]] <- if (is.finite(r)) r^2 else 0
    }
  }
  names(out) <- colnames(Z)[seq_len(k)] %||% paste0("Factor", seq_len(k))
  out
}

project_latent_equal_view <- function(X_list, W_list, ridge_rel) {
  views <- intersect(CANONICAL_VIEWS, intersect(names(X_list), names(W_list)))
  if (!length(views)) stop("No hay vistas para proyectar al espacio latente.")
  sample_ids <- colnames(X_list[[views[[1L]]]])
  K <- ncol(W_list[[views[[1L]]]])
  X_blocks <- list()
  W_blocks <- list()

  for (v in views) {
    Xv <- assert_numeric_matrix(X_list[[v]], paste0("X ", v))
    Wv <- assert_numeric_matrix(W_list[[v]], paste0("W ", v))
    if (nrow(Xv) != nrow(Wv)) stop("X y W no coinciden en features para ", v, ".")
    if (ncol(Wv) != K) stop("Numero de factores inconsistente en ", v, ".")
    if (!setequal(colnames(Xv), sample_ids)) stop("IDs de muestra distintos entre vistas.")
    Xv <- Xv[, sample_ids, drop = FALSE]
    block_scale <- sqrt(max(nrow(Wv), 1L))
    X_blocks[[v]] <- Xv / block_scale
    W_blocks[[v]] <- Wv / block_scale
  }

  X_stack <- do.call(rbind, X_blocks)
  W_stack <- do.call(rbind, W_blocks)
  s <- svd(W_stack, nu = min(dim(W_stack)), nv = min(dim(W_stack)))
  if (!length(s$d) || max(s$d) <= 0) stop("W apilado no contiene señal proyectable.")
  lambda <- ridge_rel * max(s$d)^2
  ux <- crossprod(s$u, X_stack)
  shrink <- s$d / (s$d^2 + lambda)
  coefficients <- s$v %*% sweep(ux, 1L, shrink, "*")
  Zhat <- t(coefficients)
  rownames(Zhat) <- sample_ids
  colnames(Zhat) <- colnames(W_list[[views[[1L]]]]) %||% paste0("Factor", seq_len(K))
  diag <- svd_diagnostics(W_stack, expected_rank = K)
  list(
    Z = Zhat,
    ridge_rel = ridge_rel,
    lambda = lambda,
    rank = diag$rank,
    condition = diag$condition,
    min_singular = diag$min_singular,
    stability = numeric_stability(diag$rank, K, diag$condition)
  )
}

extract_z <- function(model) {
  z <- MOFA2::get_factors(model, as.data.frame = FALSE)
  while (is.list(z) && length(z) == 1L && !is.matrix(z)) z <- z[[1L]]
  if (is.list(z) && !is.matrix(z)) {
    if (length(z) != 1L) stop("Se esperaban un solo grupo MOFA para Z.")
    z <- z[[1L]]
  }
  z <- assert_numeric_matrix(z, "Z de MOFA")
  if (!is.null(rownames(z)) && anyDuplicated(rownames(z))) {
    stop("Z de MOFA contiene IDs de muestra duplicados.")
  }
  z
}

extract_w <- function(model) {
  w <- MOFA2::get_expectations(model, variable = "W")
  while (
    is.list(w) && length(w) == 1L && is.list(w[[1L]]) &&
      !is.matrix(w[[1L]]) && !is.data.frame(w[[1L]])
  ) {
    w <- w[[1L]]
  }
  if (!is.list(w) || !length(w)) stop("No se pudo extraer W como lista de vistas.")

  out <- list()
  for (i in seq_along(w)) {
    item <- w[[i]]
    while (is.list(item) && length(item) == 1L && !is.matrix(item)) item <- item[[1L]]
    nm <- if (is.null(names(w)) || length(names(w)) < i) "" else names(w)[[i]]
    cv <- canonical_view(nm)
    if (is.na(cv)) stop("Nombre de vista MOFA no reconocido: ", nm)
    if (cv %in% names(out)) stop("Vista MOFA duplicada tras normalizacion: ", cv)
    item <- assert_numeric_matrix(item, paste0("W ", cv))
    if (!is.null(rownames(item)) && anyDuplicated(rownames(item))) {
      stop("W de ", cv, " contiene features duplicadas.")
    }
    out[[cv]] <- item
  }
  if (!all(CANONICAL_VIEWS %in% names(out))) {
    stop("W no contiene las cuatro vistas canonicas.")
  }
  out[CANONICAL_VIEWS]
}

align_z_to_samples <- function(Z, sample_ids) {
  if (nrow(Z) != length(sample_ids) && ncol(Z) == length(sample_ids)) Z <- t(Z)
  if (nrow(Z) != length(sample_ids)) stop("Z no coincide con el numero de muestras TRAIN.")
  if (is.null(rownames(Z))) {
    rownames(Z) <- sample_ids
  } else if (setequal(rownames(Z), sample_ids)) {
    Z <- Z[sample_ids, , drop = FALSE]
  } else {
    stop("Los IDs de Z no coinciden con TRAIN.")
  }
  if (is.null(colnames(Z))) colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))
  Z
}

extract_rfm_views <- function(rfm, split = c("train", "test")) {
  split <- match.arg(split)
  out <- list()
  for (v in CANONICAL_VIEWS) {
    code <- VIEW_CODE[[v]]
    obj <- rfm[[code]]
    if (is.null(obj) || is.null(obj[[split]])) {
      if (split == "train") stop("Falta rfm$", code, "$", split, ".")
      return(NULL)
    }
    X <- assert_numeric_matrix(obj[[split]], paste0("rfm$", code, "$", split))
    if (is.null(rownames(X)) || is.null(colnames(X))) {
      stop("La matriz ", code, "/", split, " necesita rownames de features y colnames de muestras.")
    }
    if (anyDuplicated(rownames(X))) stop("Features duplicadas en ", code, "/", split, ".")
    if (anyDuplicated(colnames(X))) stop("IDs de muestra duplicados en ", code, "/", split, ".")
    out[[v]] <- X
  }
  out
}

align_x_and_w <- function(X_list, W_list, Z = NULL) {
  out_x <- list()
  out_w <- list()
  first_ids <- NULL
  K <- if (!is.null(Z)) ncol(Z) else ncol(W_list[[1L]])

  for (v in CANONICAL_VIEWS) {
    X <- X_list[[v]]
    W <- W_list[[v]]
    if (ncol(W) != K && nrow(W) == K) W <- t(W)
    if (ncol(W) != K) stop("W de ", v, " no coincide con K factores.")
    if (is.null(rownames(W))) stop("W de ", v, " no tiene nombres de features.")
    if (!all(rownames(W) %in% rownames(X))) {
      missing <- setdiff(rownames(W), rownames(X))
      stop("X de ", v, " no contiene ", length(missing), " features usadas por W.")
    }
    X <- X[rownames(W), , drop = FALSE]
    if (is.null(first_ids)) first_ids <- colnames(X)
    if (!setequal(colnames(X), first_ids)) stop("Las vistas no contienen las mismas muestras.")
    X <- X[, first_ids, drop = FALSE]
    if (!is.null(colnames(W)) && !is.null(Z) && setequal(colnames(W), colnames(Z))) {
      W <- W[, colnames(Z), drop = FALSE]
    } else if (!is.null(Z)) {
      colnames(W) <- colnames(Z)
    }
    out_x[[v]] <- X
    out_w[[v]] <- W
  }
  list(X = out_x, W = out_w, sample_ids = first_ids)
}

extract_groups <- function(rfm, split, ids) {
  obj <- rfm[[paste0("grupo_", split)]]
  value <- NULL
  if (is.data.frame(obj) || is.matrix(obj)) {
    obj <- as.data.frame(obj, stringsAsFactors = FALSE)
    preferred <- intersect(c("grupo", "Group", "group", "y"), names(obj))
    col <- if (length(preferred)) preferred[[1L]] else if (ncol(obj)) names(obj)[[1L]] else NULL
    if (!is.null(col)) {
      if (!is.null(rownames(obj)) && all(ids %in% rownames(obj))) {
        value <- obj[ids, col]
      } else if (nrow(obj) == length(ids)) {
        value <- obj[[col]]
      }
    }
  } else if (!is.null(obj) && length(obj) == length(ids)) {
    value <- obj
    if (!is.null(names(obj)) && all(ids %in% names(obj))) value <- obj[ids]
  }
  if (is.null(value)) value <- sub("_.*", "", ids)
  as.character(value)
}

compute_reconstruction_r2 <- function(X_list, W_list, Z) {
  rows <- list()
  i <- 0L
  for (v in CANONICAL_VIEWS) {
    X <- X_list[[v]]
    W <- W_list[[v]]
    obs_centered <- X - rowMeans(X)
    tss <- rowSums(obs_centered^2)
    valid <- is.finite(tss) & tss > sqrt(.Machine$double.eps)
    pred_all <- W %*% t(Z)
    sse <- rowSums((X - pred_all)^2)
    r2_feature <- rep(NA_real_, nrow(X))
    r2_feature[valid] <- 1 - sse[valid] / tss[valid]
    total_view <- 100 * mean(r2_feature[valid], na.rm = TRUE)
    obs_var <- apply(X, 1L, stats::var)

    for (k in seq_len(ncol(Z))) {
      pred_var <- (W[, k]^2) * stats::var(Z[, k])
      ratio <- pred_var[valid] / obs_var[valid]
      ratio <- ratio[is.finite(ratio)]
      i <- i + 1L
      rows[[i]] <- data.frame(
        view = v,
        factor = colnames(Z)[[k]],
        r2_factor_percent = if (length(ratio)) 100 * mean(ratio) else 0,
        r2_total_view_percent = total_view,
        n_features_full = nrow(X),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows_fill(rows)
}

extract_mofa_variance_table <- function(api_variance, Z, X_list) {
  tryCatch({
    if (!is.list(api_variance)) stop("Objeto de varianza MOFA no reconocido.")
    total <- api_variance$r2_total
    per_factor <- api_variance$r2_per_factor
    while (is.list(total) && length(total) == 1L) total <- total[[1L]]
    while (is.list(per_factor) && length(per_factor) == 1L) per_factor <- per_factor[[1L]]
    if (is.null(total) || is.null(per_factor)) stop("Faltan r2_total/r2_per_factor.")

    total_names <- names(total)
    total <- as.numeric(total)
    if (length(total) != length(CANONICAL_VIEWS)) stop("r2_total no tiene cuatro vistas.")
    mapped_total_names <- canonical_view(total_names)
    if (length(mapped_total_names) == length(total) && !anyNA(mapped_total_names) &&
        setequal(mapped_total_names, CANONICAL_VIEWS)) {
      names(total) <- mapped_total_names
      total <- total[CANONICAL_VIEWS]
    } else {
      names(total) <- CANONICAL_VIEWS
    }

    per_factor <- as.matrix(per_factor)
    storage.mode(per_factor) <- "double"
    row_views <- canonical_view(rownames(per_factor))
    col_views <- canonical_view(colnames(per_factor))
    rows_are_views <- length(row_views) == nrow(per_factor) && !anyNA(row_views) &&
      setequal(row_views, CANONICAL_VIEWS)
    cols_are_views <- length(col_views) == ncol(per_factor) && !anyNA(col_views) &&
      setequal(col_views, CANONICAL_VIEWS)
    if (cols_are_views && !rows_are_views) {
      per_factor <- t(per_factor)
      row_views <- canonical_view(rownames(per_factor))
      rows_are_views <- TRUE
    }
    if (rows_are_views) {
      rownames(per_factor) <- row_views
      per_factor <- per_factor[CANONICAL_VIEWS, , drop = FALSE]
    } else if (nrow(per_factor) == length(CANONICAL_VIEWS)) {
      rownames(per_factor) <- CANONICAL_VIEWS
    } else if (ncol(per_factor) == length(CANONICAL_VIEWS)) {
      per_factor <- t(per_factor)
      rownames(per_factor) <- CANONICAL_VIEWS
    } else {
      stop("No se pudo orientar r2_per_factor por vistas.")
    }

    factor_names <- colnames(Z)
    if (ncol(per_factor) != length(factor_names)) stop("r2_per_factor no coincide con Z.")
    if (!is.null(colnames(per_factor)) && all(factor_names %in% colnames(per_factor))) {
      per_factor <- per_factor[, factor_names, drop = FALSE]
    } else {
      colnames(per_factor) <- factor_names
    }
    if (any(!is.finite(total)) || any(!is.finite(per_factor))) {
      stop("La varianza explicada de MOFA contiene valores no finitos.")
    }

    rows <- vector("list", length(CANONICAL_VIEWS) * length(factor_names))
    i <- 0L
    for (v in CANONICAL_VIEWS) {
      for (f in factor_names) {
        i <- i + 1L
        rows[[i]] <- data.frame(
          view = v,
          factor = f,
          r2_factor_percent = per_factor[v, f],
          r2_total_view_percent = total[[v]],
          n_features_full = nrow(X_list[[v]]),
          variance_source = "MOFA2::get_variance_explained",
          stringsAsFactors = FALSE
        )
      }
    }
    bind_rows_fill(rows)
  }, error = function(e) NULL)
}

summarize_full_quality <- function(r2_table, Z) {
  view_r2 <- tapply(r2_table$r2_total_view_percent, r2_table$view, function(x) x[[1L]])
  view_r2 <- view_r2[CANONICAL_VIEWS]
  view_features <- tapply(r2_table$n_features_full, r2_table$view, function(x) x[[1L]])
  view_features <- view_features[CANONICAL_VIEWS]
  r2_matrix <- xtabs(r2_factor_percent ~ view + factor, data = r2_table)
  r2_matrix <- r2_matrix[CANONICAL_VIEWS, colnames(Z), drop = FALSE]
  r2_pos <- pmax(r2_matrix, 0)
  factor_strength <- colSums(r2_pos)
  sharedness <- vapply(seq_len(ncol(r2_pos)), function(k) {
    effective_number(r2_pos[, k]) / nrow(r2_pos)
  }, numeric(1))
  zdiag <- svd_diagnostics(scale(Z, center = TRUE, scale = FALSE), expected_rank = ncol(Z))
  z_energy <- zdiag$singular_values^2

  list(
    metrics = list(
      n_factors = ncol(Z),
      n_full_features_total = sum(view_features),
      n_full_features_transcriptomica = view_features[["transcriptomica"]],
      n_full_features_proteomica = view_features[["proteomica"]],
      n_full_features_metabolomica = view_features[["metabolomica"]],
      n_full_features_clinical = view_features[["clinical"]],
      full_r2_view_mean = mean(view_r2, na.rm = TRUE),
      full_r2_view_min = min(view_r2, na.rm = TRUE),
      full_view_balance = effective_number(pmax(view_r2, 0)) / length(view_r2),
      full_factor_r2_mean = mean(factor_strength, na.rm = TRUE),
      full_factor_r2_min = min(factor_strength, na.rm = TRUE),
      full_factor_balance = effective_number(factor_strength) / max(length(factor_strength), 1L),
      full_factor_sharedness_mean = mean(sharedness, na.rm = TRUE),
      full_factor_sharedness_min = min(sharedness, na.rm = TRUE),
      full_z_effective_rank_fraction = effective_number(z_energy) / max(ncol(Z), 1L),
      full_z_rank_fraction = zdiag$rank / max(ncol(Z), 1L),
      full_z_condition = zdiag$condition,
      full_z_numeric_stability = numeric_stability(zdiag$rank, ncol(Z), zdiag$condition)
    ),
    factor_strength = factor_strength,
    factor_sharedness = sharedness,
    r2_matrix = r2_matrix,
    view_r2 = view_r2
  )
}

build_full_context <- function(shared_dir, analysis, ridge_grid, do_test = TRUE) {
  ready_file <- ready_file_for(shared_dir, analysis)
  mofa_file <- mofa_file_for(shared_dir, analysis)
  if (is.na(ready_file)) stop("No se encontro ready_for_modeling para ", analysis, ".")
  if (is.na(mofa_file)) stop("No se encontro modelo_mofa_final para ", analysis, ".")

  rfm_read <- safe_read_rds(ready_file)
  if (!rfm_read$ok) stop("No se pudo leer RFM: ", rfm_read$error)
  model_read <- safe_read_rds(mofa_file)
  if (!model_read$ok) stop("No se pudo leer MOFA: ", model_read$error)
  rfm <- rfm_read$value
  model <- model_read$value

  W <- extract_w(model)
  X_train_raw <- extract_rfm_views(rfm, "train")
  train_ids <- colnames(X_train_raw[[1L]])
  if (length(train_ids) < 6L) stop("Se requieren al menos seis muestras TRAIN para la estabilidad LOO.")
  Z <- align_z_to_samples(extract_z(model), train_ids)
  aligned <- align_x_and_w(X_train_raw, W, Z)
  X_train <- aligned$X
  W <- aligned$W
  Z <- Z[aligned$sample_ids, , drop = FALSE]

  X_test <- NULL
  test_error <- NA_character_
  if (do_test) {
    X_test <- tryCatch({
      raw <- extract_rfm_views(rfm, "test")
      aligned_test <- align_x_and_w(raw, W, Z = NULL)
      overlap <- intersect(aligned$sample_ids, aligned_test$sample_ids)
      if (length(overlap)) {
        stop(
          "TRAIN y TEST comparten ", length(overlap),
          " IDs; se cancela exclusivamente la proyeccion TEST."
        )
      }
      aligned_test$X
    }, error = function(e) {
      test_error <<- conditionMessage(e)
      NULL
    })
  }

  observed_view_kernels <- lapply(X_train, kernel_from_features, standardize_rows = TRUE)
  reconstructed <- lapply(CANONICAL_VIEWS, function(v) W[[v]] %*% t(Z))
  names(reconstructed) <- CANONICAL_VIEWS
  reconstructed_view_kernels <- lapply(reconstructed, kernel_from_features, standardize_rows = FALSE)
  K_observed_full <- combine_view_kernels(observed_view_kernels)
  K_reconstructed_full <- combine_view_kernels(reconstructed_view_kernels)
  K_z_native <- kernel_from_latent(Z, balance_factors = FALSE)
  K_z_balanced <- kernel_from_latent(Z, balance_factors = TRUE)

  api_variance <- tryCatch(MOFA2::get_variance_explained(model), error = function(e) NULL)
  r2_table_recomputed <- compute_reconstruction_r2(X_train, W, Z)
  r2_table_recomputed$variance_source <- "recalculada_desde_X_W_Z"
  r2_table_api <- extract_mofa_variance_table(api_variance, Z, X_train)
  r2_table <- if (!is.null(r2_table_api) && nrow(r2_table_api)) {
    r2_table_api
  } else {
    r2_table_recomputed
  }
  variance_source <- unique(r2_table$variance_source)[[1L]]
  full_quality <- summarize_full_quality(r2_table, Z)
  ridge_central <- ridge_grid[[ceiling(length(ridge_grid) / 2)]]
  full_projection_train <- project_latent_equal_view(X_train, W, ridge_central)
  full_projection_test <- NULL
  if (!is.null(X_test)) {
    full_projection_test <- tryCatch(
      project_latent_equal_view(X_test, W, ridge_central),
      error = function(e) {
        test_error <<- conditionMessage(e)
        NULL
      }
    )
  }

  list(
    analysis = analysis,
    shared_dir = shared_dir,
    ready_file = ready_file,
    mofa_file = mofa_file,
    model = model,
    rfm = rfm,
    Z = Z,
    W = W,
    X_train = X_train,
    X_test = X_test,
    train_ids = rownames(Z),
    test_ids = if (!is.null(X_test)) colnames(X_test[[1L]]) else character(),
    group_train = extract_groups(rfm, "train", rownames(Z)),
    group_test = if (!is.null(X_test)) extract_groups(rfm, "test", colnames(X_test[[1L]])) else character(),
    observed_view_kernels = observed_view_kernels,
    reconstructed_view_kernels = reconstructed_view_kernels,
    K_observed_full = K_observed_full,
    K_reconstructed_full = K_reconstructed_full,
    K_z_native = K_z_native,
    K_z_balanced = K_z_balanced,
    r2_table = r2_table,
    r2_table_api = r2_table_api,
    r2_table_recomputed = r2_table_recomputed,
    variance_source = variance_source,
    api_variance = api_variance,
    full_quality = full_quality,
    full_projection_train = full_projection_train,
    full_projection_test = full_projection_test,
    ridge_central = ridge_central,
    test_error = test_error
  )
}

candidate_mapping_rows <- function(panel_tbl, diagnostic_obj) {
  panel_features <- unique(as.character(panel_tbl$feature_model))
  rows <- list()
  i <- 0L

  weights <- NULL
  if (is.list(diagnostic_obj) && !is.null(diagnostic_obj$df_weights)) {
    weights <- tryCatch(as.data.frame(diagnostic_obj$df_weights), error = function(e) NULL)
  }

  if (!is.null(weights) && all(c("feature_model", "feature_id") %in% names(weights))) {
    view_col <- intersect(c("view_original", "view", "View", "view_code"), names(weights))
    if (length(view_col)) {
      sub <- weights[as.character(weights$feature_model) %in% panel_features, , drop = FALSE]
      sub$view_canonical <- canonical_view(sub[[view_col[[1L]]]])
      sub <- sub[!is.na(sub$view_canonical), , drop = FALSE]
      if (nrow(sub)) {
        out <- unique(data.frame(
          feature_model = as.character(sub$feature_model),
          view = sub$view_canonical,
          feature_id = as.character(sub$feature_id),
          mapping_source = "diagnostic_objects$df_weights",
          stringsAsFactors = FALSE
        ))
        out <- out[nzchar(out$feature_id) & !is.na(out$feature_id), , drop = FALSE]
        if (nrow(out)) return(out)
      }
    }
  }

  view_col <- intersect(c("view", "view_code", "View"), names(panel_tbl))
  if (!length(view_col)) stop("La tabla Bayes no contiene una columna de vista.")
  candidate_cols <- intersect(
    c("feature_original", "feature_id", "feature_label", "label_pretty"),
    names(panel_tbl)
  )

  for (r in seq_len(nrow(panel_tbl))) {
    fm <- as.character(panel_tbl$feature_model[[r]])
    view <- canonical_view(panel_tbl[[view_col[[1L]]]][[r]])
    candidates <- as.character(unlist(panel_tbl[r, candidate_cols, drop = FALSE], use.names = FALSE))
    candidates <- c(candidates, sub("^[a-zA-Z]{2}_", "", fm))
    candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
    for (id in candidates) {
      i <- i + 1L
      rows[[i]] <- data.frame(
        feature_model = fm,
        view = view,
        feature_id = id,
        mapping_source = "fallback_panel_columns",
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows_fill(rows)
}

load_panel_for_context <- function(branch_dir, ctx) {
  panel_file <- panel_csv_for(branch_dir)
  diag_file <- diagnostic_objects_for(branch_dir)
  if (!file.exists(panel_file)) stop("No existe el CSV diagnostico Bayes: ", panel_file)
  if (!file.exists(diag_file)) stop("No existe diagnostic_objects.rds: ", diag_file)

  panel_read <- safe_read_csv(panel_file)
  if (!panel_read$ok) stop("No se pudo leer el panel Bayes: ", panel_read$error)
  panel_tbl <- as.data.frame(panel_read$value, stringsAsFactors = FALSE)
  if (!"feature_model" %in% names(panel_tbl)) stop("El panel no contiene feature_model.")
  panel_tbl <- panel_tbl[!duplicated(as.character(panel_tbl$feature_model)), , drop = FALSE]
  panel_tbl <- panel_tbl[!is.na(panel_tbl$feature_model) & nzchar(panel_tbl$feature_model), , drop = FALSE]
  if (!nrow(panel_tbl)) stop("El panel Bayes esta vacio.")
  numeric_cols <- vapply(panel_tbl, is.numeric, logical(1))
  if (any(numeric_cols) && any(!is.finite(as.matrix(panel_tbl[, numeric_cols, drop = FALSE])))) {
    stop("La tabla Bayes contiene valores numericos no finitos.")
  }

  diag_read <- safe_read_rds(diag_file)
  if (!diag_read$ok) stop("No se pudo leer diagnostic_objects: ", diag_read$error)
  if (!is.list(diag_read$value)) stop("diagnostic_objects no es una lista.")
  diagnostic_panel <- diag_read$value$panel_tbl
  if (!is.null(diagnostic_panel)) {
    diagnostic_panel <- tryCatch(as.data.frame(diagnostic_panel), error = function(e) NULL)
    if (is.null(diagnostic_panel) || !"feature_model" %in% names(diagnostic_panel)) {
      stop("diagnostic_objects$panel_tbl no es auditable mediante feature_model.")
    }
    csv_set <- sort(unique(as.character(panel_tbl$feature_model)))
    rds_set <- sort(unique(as.character(diagnostic_panel$feature_model)))
    csv_set <- csv_set[!is.na(csv_set) & nzchar(csv_set)]
    rds_set <- rds_set[!is.na(rds_set) & nzchar(rds_set)]
    if (!identical(csv_set, rds_set)) {
      stop("El set Bayes del CSV no coincide con diagnostic_objects$panel_tbl.")
    }
  }
  mapping <- candidate_mapping_rows(panel_tbl, diag_read$value)
  if (!nrow(mapping)) stop("No se pudo construir el mapeo feature_model -> feature_id.")

  selected_ids <- list()
  selected_models <- list()
  mapping$matched <- FALSE
  for (v in CANONICAL_VIEWS) {
    idx <- which(mapping$view == v)
    ids <- if (length(idx)) unique(mapping$feature_id[idx]) else character()
    matched_ids <- intersect(ids, intersect(rownames(ctx$W[[v]]), rownames(ctx$X_train[[v]])))
    selected_ids[[v]] <- matched_ids
    if (length(idx) && length(matched_ids)) {
      mapping$matched[idx] <- mapping$feature_id[idx] %in% matched_ids
    }
    selected_models[[v]] <- unique(mapping$feature_model[idx][mapping$feature_id[idx] %in% matched_ids])
    if (!length(matched_ids)) {
      stop("La vista ", v, " no tiene features Bayes mapeadas simultaneamente a X y W.")
    }
  }

  matched_models <- unique(mapping$feature_model[mapping$matched])
  panel_models <- unique(as.character(panel_tbl$feature_model))
  mapping_fraction <- length(intersect(panel_models, matched_models)) / length(panel_models)
  if (!is.finite(mapping_fraction) || mapping_fraction <= 0) stop("Ninguna feature Bayes pudo mapearse.")

  X_train <- lapply(CANONICAL_VIEWS, function(v) {
    ctx$X_train[[v]][selected_ids[[v]], , drop = FALSE]
  })
  W <- lapply(CANONICAL_VIEWS, function(v) {
    ctx$W[[v]][selected_ids[[v]], , drop = FALSE]
  })
  names(X_train) <- names(W) <- CANONICAL_VIEWS

  X_test <- NULL
  if (!is.null(ctx$X_test)) {
    X_test <- lapply(CANONICAL_VIEWS, function(v) {
      ctx$X_test[[v]][selected_ids[[v]], , drop = FALSE]
    })
    names(X_test) <- CANONICAL_VIEWS
  }

  list(
    panel_file = panel_file,
    diagnostic_file = diag_file,
    panel_table = panel_tbl,
    mapping = mapping,
    selected_ids = selected_ids,
    selected_models = selected_models,
    matched_models = matched_models,
    mapping_fraction = mapping_fraction,
    X_train = X_train,
    X_test = X_test,
    W = W
  )
}

loading_coverage_tables <- function(ctx, panel) {
  rows <- list()
  i <- 0L
  for (v in CANONICAL_VIEWS) {
    Wfull <- ctx$W[[v]]
    Wpanel <- panel$W[[v]]
    for (k in seq_len(ncol(Wfull))) {
      denom <- sum(Wfull[, k]^2)
      num <- sum(Wpanel[, k]^2)
      coverage <- if (denom > 0) num / denom else 1
      i <- i + 1L
      rows[[i]] <- data.frame(
        view = v,
        factor = colnames(ctx$Z)[[k]],
        loading_energy_full = denom,
        loading_energy_panel = num,
        loading_energy_coverage = min(max(coverage, 0), 1),
        stringsAsFactors = FALSE
      )
    }
  }
  tbl <- bind_rows_fill(rows)

  factor_cov <- vapply(colnames(ctx$Z), function(f) {
    z <- tbl[tbl$factor == f, , drop = FALSE]
    denom <- sum(z$loading_energy_full)
    if (denom > 0) sum(z$loading_energy_panel) / denom else 1
  }, numeric(1))
  view_cov <- vapply(CANONICAL_VIEWS, function(v) {
    z <- tbl[tbl$view == v, , drop = FALSE]
    denom <- sum(z$loading_energy_full)
    if (denom > 0) sum(z$loading_energy_panel) / denom else 1
  }, numeric(1))

  list(
    table = tbl,
    factor_coverage = factor_cov,
    view_coverage = view_cov,
    factor_min = min(factor_cov),
    factor_mean = mean(factor_cov),
    view_min = min(view_cov),
    view_mean = mean(view_cov)
  )
}

evaluate_panel <- function(ctx, panel, ridge_grid, do_test = TRUE, keep_details = FALSE) {
  observed_panel_view_kernels <- lapply(panel$X_train, kernel_from_features, standardize_rows = TRUE)
  reconstructed_panel <- lapply(CANONICAL_VIEWS, function(v) panel$W[[v]] %*% t(ctx$Z))
  names(reconstructed_panel) <- CANONICAL_VIEWS
  reconstructed_panel_view_kernels <- lapply(
    reconstructed_panel, kernel_from_features, standardize_rows = FALSE
  )
  K_observed_panel <- combine_view_kernels(observed_panel_view_kernels)
  K_reconstructed_panel <- combine_view_kernels(reconstructed_panel_view_kernels)

  ridge_results <- lapply(ridge_grid, function(ridge) {
    projection <- project_latent_equal_view(panel$X_train, panel$W, ridge)
    Kzhat_native <- kernel_from_latent(projection$Z, balance_factors = FALSE)
    Kzhat_balanced <- kernel_from_latent(projection$Z, balance_factors = TRUE)
    factor_r2 <- factor_recovery_r2(ctx$Z, projection$Z)
    data.frame(
      ridge_rel = ridge,
      zhat_cka_z_native = linear_cka(Kzhat_native, ctx$K_z_native),
      zhat_cka_z_balanced = linear_cka(Kzhat_balanced, ctx$K_z_balanced),
      zhat_cka_full_projection = linear_cka(
        Kzhat_balanced,
        kernel_from_latent(ctx$full_projection_train$Z, balance_factors = TRUE)
      ),
      zhat_procrustes = procrustes_similarity(ctx$Z, projection$Z),
      factor_recovery_min = min(factor_r2),
      factor_recovery_mean = mean(factor_r2),
      projection_rank = projection$rank,
      projection_condition = projection$condition,
      projection_stability = projection$stability,
      stringsAsFactors = FALSE
    )
  })
  ridge_table <- bind_rows_fill(ridge_results)
  central_idx <- which.min(abs(ridge_table$ridge_rel - ctx$ridge_central))
  central_projection <- project_latent_equal_view(
    panel$X_train, panel$W, ridge_table$ridge_rel[[central_idx]]
  )
  K_zhat_balanced <- kernel_from_latent(central_projection$Z, balance_factors = TRUE)
  central_factor_r2 <- factor_recovery_r2(ctx$Z, central_projection$Z)

  coverage <- loading_coverage_tables(ctx, panel)
  view_rows <- list()
  for (v in CANONICAL_VIEWS) {
    diag_v <- svd_diagnostics(panel$W[[v]], expected_rank = ncol(ctx$Z))
    view_rows[[v]] <- data.frame(
      view = v,
      panel_features_matched = nrow(panel$W[[v]]),
      observed_cka_full_view = linear_cka(
        observed_panel_view_kernels[[v]], ctx$observed_view_kernels[[v]]
      ),
      observed_cka_z_balanced = linear_cka(
        observed_panel_view_kernels[[v]], ctx$K_z_balanced
      ),
      reconstructed_cka_full_view = linear_cka(
        reconstructed_panel_view_kernels[[v]], ctx$reconstructed_view_kernels[[v]]
      ),
      reconstructed_cka_z_balanced = linear_cka(
        reconstructed_panel_view_kernels[[v]], ctx$K_z_balanced
      ),
      panel_W_rank = diag_v$rank,
      panel_W_rank_fraction = diag_v$rank / ncol(ctx$Z),
      panel_W_condition = diag_v$condition,
      panel_W_numeric_stability = numeric_stability(
        diag_v$rank, ncol(ctx$Z), diag_v$condition
      ),
      loading_energy_coverage_view = coverage$view_coverage[[v]],
      stringsAsFactors = FALSE
    )
  }
  view_table <- bind_rows_fill(view_rows)

  factor_strength <- ctx$full_quality$factor_strength
  factor_weights <- if (sum(pmax(factor_strength, 0)) > 0) {
    pmax(factor_strength, 0) / sum(pmax(factor_strength, 0))
  } else {
    rep(1 / length(factor_strength), length(factor_strength))
  }
  factor_weighted <- sum(central_factor_r2[names(factor_strength)] * factor_weights, na.rm = TRUE)

  test_valid <- FALSE
  test_metrics <- list(
    test_cka_panel_vs_full_projection = NA_real_,
    test_distance_panel_vs_full_projection = NA_real_,
    test_factor_recovery_min = NA_real_,
    test_factor_recovery_mean = NA_real_
  )
  panel_projection_test <- NULL
  if (do_test && !is.null(panel$X_test) && !is.null(ctx$full_projection_test)) {
    panel_projection_test <- tryCatch(
      project_latent_equal_view(panel$X_test, panel$W, ctx$ridge_central),
      error = function(e) NULL
    )
    if (!is.null(panel_projection_test)) {
      K_test_panel <- kernel_from_latent(panel_projection_test$Z, balance_factors = TRUE)
      K_test_full <- kernel_from_latent(ctx$full_projection_test$Z, balance_factors = TRUE)
      test_factor <- factor_recovery_r2(ctx$full_projection_test$Z, panel_projection_test$Z)
      test_metrics <- list(
        test_cka_panel_vs_full_projection = linear_cka(K_test_panel, K_test_full),
        test_distance_panel_vs_full_projection = kernel_distance_correlation(K_test_panel, K_test_full),
        test_factor_recovery_min = min(test_factor),
        test_factor_recovery_mean = mean(test_factor)
      )
      test_valid <- all(is.finite(unlist(test_metrics)))
    }
  }

  metrics <- c(
    list(
      n_panel_features = nrow(panel$panel_table),
      n_panel_features_matched = length(panel$matched_models),
      mapping_fraction = panel$mapping_fraction,
      observed_cka_full = linear_cka(K_observed_panel, ctx$K_observed_full),
      observed_cka_z_native = linear_cka(K_observed_panel, ctx$K_z_native),
      observed_cka_z_balanced = linear_cka(K_observed_panel, ctx$K_z_balanced),
      observed_distance_full = kernel_distance_correlation(K_observed_panel, ctx$K_observed_full),
      observed_view_min = min(view_table$observed_cka_full_view),
      observed_view_mean = mean(view_table$observed_cka_full_view),
      observed_loo_q10 = loo_cka_q10(K_observed_panel, ctx$K_observed_full),
      reconstructed_cka_full = linear_cka(K_reconstructed_panel, ctx$K_reconstructed_full),
      reconstructed_cka_z_native = linear_cka(K_reconstructed_panel, ctx$K_z_native),
      reconstructed_cka_z_balanced = linear_cka(K_reconstructed_panel, ctx$K_z_balanced),
      reconstructed_distance_full = kernel_distance_correlation(
        K_reconstructed_panel, ctx$K_reconstructed_full
      ),
      reconstructed_view_min = min(view_table$reconstructed_cka_full_view),
      reconstructed_view_mean = mean(view_table$reconstructed_cka_full_view),
      reconstructed_loo_q10 = loo_cka_q10(K_reconstructed_panel, ctx$K_reconstructed_full),
      zhat_cka_z_native = ridge_table$zhat_cka_z_native[[central_idx]],
      zhat_cka_z_balanced = ridge_table$zhat_cka_z_balanced[[central_idx]],
      zhat_cka_full_projection = ridge_table$zhat_cka_full_projection[[central_idx]],
      zhat_procrustes = ridge_table$zhat_procrustes[[central_idx]],
      zhat_ridge_cka_min = min(ridge_table$zhat_cka_z_balanced),
      zhat_ridge_cka_median = stats::median(ridge_table$zhat_cka_z_balanced),
      zhat_loo_q10 = loo_cka_q10(K_zhat_balanced, ctx$K_z_balanced),
      factor_recovery_min = min(central_factor_r2),
      factor_recovery_mean = mean(central_factor_r2),
      factor_recovery_weighted = factor_weighted,
      loading_coverage_factor_min = coverage$factor_min,
      loading_coverage_factor_mean = coverage$factor_mean,
      loading_coverage_view_min = coverage$view_min,
      loading_coverage_view_mean = coverage$view_mean,
      panel_joint_rank = central_projection$rank,
      panel_joint_rank_fraction = central_projection$rank / ncol(ctx$Z),
      panel_joint_condition = central_projection$condition,
      panel_numeric_stability = central_projection$stability,
      panel_view_rank_min = min(view_table$panel_W_rank_fraction),
      panel_view_numeric_stability_min = min(view_table$panel_W_numeric_stability),
      test_projection_valid = test_valid
    ),
    test_metrics
  )

  details <- NULL
  if (keep_details) {
    details <- list(
      panel_table = panel$panel_table,
      panel_mapping = panel$mapping,
      selected_feature_ids = panel$selected_ids,
      W_panel = panel$W,
      Z_full_train = ctx$Z,
      Z_full_projection_train = ctx$full_projection_train$Z,
      Z_panel_recovered_train = central_projection$Z,
      Z_full_projection_test = if (!is.null(ctx$full_projection_test)) ctx$full_projection_test$Z else NULL,
      Z_panel_recovered_test = if (!is.null(panel_projection_test)) panel_projection_test$Z else NULL,
      K_observed_full = ctx$K_observed_full,
      K_observed_panel = K_observed_panel,
      K_reconstructed_full = ctx$K_reconstructed_full,
      K_reconstructed_panel = K_reconstructed_panel,
      K_z_full = ctx$K_z_balanced,
      K_z_panel_recovered = K_zhat_balanced,
      ridge_table = ridge_table,
      factor_recovery = central_factor_r2,
      coverage = coverage,
      view_table = view_table
    )
  }

  list(
    metrics = metrics,
    view_table = view_table,
    coverage_table = coverage$table,
    ridge_table = ridge_table,
    details = details
  )
}

compare_to_q00 <- function(ctx, qctx) {
  common <- intersect(ctx$train_ids, qctx$train_ids)
  if (length(common) < 6L) {
    return(list(
      q00_common_train = length(common),
      q00_cka_z = NA_real_,
      q00_cka_observed_full = NA_real_,
      q00_cka_reconstructed_full = NA_real_,
      q00_distance_z = NA_real_
    ))
  }
  idx <- match(common, ctx$train_ids)
  idx0 <- match(common, qctx$train_ids)
  Kz <- ctx$K_z_balanced[idx, idx, drop = FALSE]
  Kz0 <- qctx$K_z_balanced[idx0, idx0, drop = FALSE]
  Ko <- ctx$K_observed_full[idx, idx, drop = FALSE]
  Ko0 <- qctx$K_observed_full[idx0, idx0, drop = FALSE]
  Kr <- ctx$K_reconstructed_full[idx, idx, drop = FALSE]
  Kr0 <- qctx$K_reconstructed_full[idx0, idx0, drop = FALSE]
  list(
    q00_common_train = length(common),
    q00_cka_z = linear_cka(Kz, Kz0),
    q00_cka_observed_full = linear_cka(Ko, Ko0),
    q00_cka_reconstructed_full = linear_cka(Kr, Kr0),
    q00_distance_z = kernel_distance_correlation(Kz, Kz0)
  )
}

flatten_named_list <- function(x) {
  out <- unlist(x, recursive = TRUE, use.names = TRUE)
  as.list(out)
}

empty_result <- function(registry_row, error_message) {
  data.frame(
    family = registry_row$family,
    family_root = registry_row$family_root,
    modelamiento = registry_row$modelamiento,
    analysis = registry_row$analysis,
    q00 = registry_row$q00,
    technical_valid = FALSE,
    technical_error = error_message,
    stringsAsFactors = FALSE
  )
}

percentile_high <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  n <- sum(ok)
  if (!n) return(out)
  if (n == 1L || length(unique(x[ok])) == 1L) {
    out[ok] <- 0.5
  } else {
    out[ok] <- (rank(x[ok], ties.method = "average") - 1) / (n - 1)
  }
  out
}

row_median_complete <- function(df, cols) {
  mat <- as.matrix(df[, cols, drop = FALSE])
  apply(mat, 1L, function(x) if (all(is.finite(x))) stats::median(x) else NA_real_)
}

pareto_fronts <- function(scores) {
  scores <- as.matrix(scores)
  n <- nrow(scores)
  fronts <- rep(NA_integer_, n)
  remaining <- which(apply(scores, 1L, function(x) all(is.finite(x))))
  front <- 1L
  while (length(remaining)) {
    nondominated <- vapply(remaining, function(i) {
      !any(vapply(setdiff(remaining, i), function(j) {
        all(scores[j, ] >= scores[i, ]) && any(scores[j, ] > scores[i, ])
      }, logical(1)))
    }, logical(1))
    current <- remaining[nondominated]
    fronts[current] <- front
    remaining <- setdiff(remaining, current)
    front <- front + 1L
  }
  fronts
}

rank_competition <- function(df) {
  quality_metrics <- c(
    "full_r2_view_mean", "full_r2_view_min", "full_view_balance",
    "full_factor_balance", "full_z_effective_rank_fraction",
    "full_z_rank_fraction", "full_z_numeric_stability"
  )
  observed_metrics <- c(
    "observed_cka_full", "observed_cka_z_balanced", "observed_distance_full",
    "observed_view_min", "observed_loo_q10"
  )
  integrated_metrics <- c(
    "reconstructed_cka_full", "reconstructed_cka_z_balanced",
    "reconstructed_distance_full", "reconstructed_view_min",
    "reconstructed_loo_q10", "zhat_cka_z_balanced",
    "zhat_cka_full_projection", "zhat_ridge_cka_min", "zhat_loo_q10",
    "factor_recovery_min", "factor_recovery_weighted"
  )
  coverage_metrics <- c(
    "mapping_fraction", "loading_coverage_factor_min",
    "loading_coverage_view_min", "panel_joint_rank_fraction",
    "panel_numeric_stability", "panel_view_rank_min",
    "panel_view_numeric_stability_min"
  )
  all_metrics <- unique(c(quality_metrics, observed_metrics, integrated_metrics, coverage_metrics))
  missing <- setdiff(all_metrics, names(df))
  for (nm in missing) df[[nm]] <- NA_real_

  valid <- as.logical(df$technical_valid)
  for (nm in all_metrics) {
    values <- suppressWarnings(as.numeric(df[[nm]]))
    values[!valid] <- NA_real_
    df[[paste0("pct_", nm)]] <- percentile_high(values)
  }

  df$domain_full_mofa <- row_median_complete(df, paste0("pct_", quality_metrics))
  df$domain_observed_recovery <- row_median_complete(df, paste0("pct_", observed_metrics))
  df$domain_integrated_recovery <- row_median_complete(df, paste0("pct_", integrated_metrics))
  df$domain_biological_coverage <- row_median_complete(df, paste0("pct_", coverage_metrics))
  domain_cols <- c(
    "domain_full_mofa", "domain_observed_recovery",
    "domain_integrated_recovery", "domain_biological_coverage"
  )
  domain_matrix <- as.matrix(df[, domain_cols, drop = FALSE])
  df$rankable <- valid & apply(domain_matrix, 1L, function(x) all(is.finite(x)))
  df$domain_minimax <- apply(domain_matrix, 1L, function(x) {
    if (all(is.finite(x))) min(x) else NA_real_
  })
  df$domain_median <- apply(domain_matrix, 1L, function(x) {
    if (all(is.finite(x))) stats::median(x) else NA_real_
  })
  df$stability_min <- apply(
    as.matrix(df[, c("observed_loo_q10", "reconstructed_loo_q10", "zhat_loo_q10"), drop = FALSE]),
    1L,
    function(x) if (all(is.finite(x))) min(x) else NA_real_
  )
  df$pareto_front <- NA_integer_
  if (any(df$rankable)) {
    df$pareto_front[df$rankable] <- pareto_fronts(domain_matrix[df$rankable, , drop = FALSE])
  }

  condition_sort <- suppressWarnings(as.numeric(df$panel_joint_condition))
  condition_sort[!is.finite(condition_sort)] <- Inf
  n_features_matched_sort <- suppressWarnings(as.numeric(df$n_panel_features_matched))
  n_features_matched_sort[!is.finite(n_features_matched_sort)] <- -Inf
  n_features_sort <- suppressWarnings(as.numeric(df$n_panel_features))
  n_features_sort[!is.finite(n_features_sort)] <- -Inf
  bad <- !df$rankable
  front_sort <- df$pareto_front
  front_sort[is.na(front_sort)] <- Inf
  order_idx <- order(
    bad,
    front_sort,
    -df$domain_minimax,
    -df$domain_median,
    -df$stability_min,
    -df$loading_coverage_factor_min,
    -df$loading_coverage_view_min,
    -df$panel_numeric_stability,
    condition_sort,
    -n_features_matched_sort,
    -n_features_sort,
    df$analysis,
    na.last = TRUE
  )
  df <- df[order_idx, , drop = FALSE]
  df$rank_selection <- seq_len(nrow(df))
  df$winner <- FALSE
  if (any(df$rankable)) {
    first_rankable <- which(df$rankable)[[1L]]
    df$winner[[first_rankable]] <- TRUE
  }
  df$tiebreak_rule <- paste(
    "Pareto -> maximin de dominios -> mediana de dominios -> estabilidad LOO ->",
    "cobertura de factores -> cobertura de vistas -> estabilidad numerica ->",
    "menor condicion -> mayor numero de features mapeadas y luego listadas solo como ultimo desempate"
  )
  df
}

rank_all_competitions <- function(raw) {
  keys <- unique(raw[, c("family", "modelamiento"), drop = FALSE])
  out <- list()
  for (i in seq_len(nrow(keys))) {
    idx <- raw$family == keys$family[[i]] & raw$modelamiento == keys$modelamiento[[i]]
    out[[i]] <- rank_competition(raw[idx, , drop = FALSE])
  }
  bind_rows_fill(out)
}

kernel_coordinates <- function(K, n_axes = 2L) {
  K <- center_kernel(K)
  e <- eigen((K + t(K)) / 2, symmetric = TRUE)
  keep <- which(e$values > sqrt(.Machine$double.eps))
  if (!length(keep)) {
    out <- matrix(0, nrow = nrow(K), ncol = n_axes)
  } else {
    keep <- keep[seq_len(min(length(keep), n_axes))]
    out <- sweep(e$vectors[, keep, drop = FALSE], 2L, sqrt(e$values[keep]), "*")
    if (ncol(out) < n_axes) out <- cbind(out, matrix(0, nrow(out), n_axes - ncol(out)))
  }
  rownames(out) <- rownames(K)
  colnames(out) <- paste0("Axis", seq_len(n_axes))
  out
}

prefix_matrix_columns <- function(x, prefix) {
  if (is.null(x)) return(NULL)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  names(x) <- paste0(prefix, "__", names(x))
  x
}

scores_table <- function(ids, groups, matrices) {
  out <- data.frame(sample_id = ids, group = groups, stringsAsFactors = FALSE)
  for (nm in names(matrices)) {
    mat <- matrices[[nm]]
    if (is.null(mat)) next
    mat <- as.matrix(mat)
    if (!is.null(rownames(mat)) && setequal(rownames(mat), ids)) mat <- mat[ids, , drop = FALSE]
    out <- cbind(out, prefix_matrix_columns(mat, nm))
  }
  out
}

full_weights_table <- function(W_list, selected_ids = NULL) {
  rows <- lapply(CANONICAL_VIEWS, function(v) {
    W <- as.matrix(W_list[[v]])
    out <- data.frame(
      view = v,
      feature_id = rownames(W),
      selected_by_bayes = if (is.null(selected_ids)) FALSE else {
        rownames(W) %in% selected_ids[[v]]
      },
      stringsAsFactors = FALSE
    )
    cbind(out, as.data.frame(W, stringsAsFactors = FALSE))
  })
  bind_rows_fill(rows)
}

plot_one_space <- function(M, title, groups, ids) {
  M <- as.matrix(M)
  if (ncol(M) < 2L) M <- cbind(M, 0)
  axis_names <- colnames(M)
  if (is.null(axis_names) || length(axis_names) < 2L) axis_names <- c("Eje 1", "Eje 2")
  groups <- as.character(groups)
  groups[is.na(groups) | !nzchar(groups)] <- "NA"
  groups <- as.factor(groups)
  cols <- grDevices::hcl.colors(max(nlevels(groups), 1L), palette = "Dark 3")
  point_cols <- cols[as.integer(groups)]
  graphics::plot(
    M[, 1], M[, 2], col = point_cols, pch = 19, cex = 1.1,
    xlab = axis_names[[1L]],
    ylab = axis_names[[2L]],
    main = title
  )
  graphics::text(M[, 1], M[, 2], labels = ids, pos = 3, cex = 0.55)
  graphics::legend("topright", legend = levels(groups), col = cols, pch = 19, cex = 0.7, bty = "n")
}

save_winner_plot <- function(ctx, details, path) {
  tryCatch({
    grDevices::png(path, width = 2600, height = 1700, res = 210)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    graphics::par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

    plot_one_space(
      ctx$Z[, seq_len(min(2L, ncol(ctx$Z))), drop = FALSE],
      "MOFA original completo · TRAIN", ctx$group_train, ctx$train_ids
    )
    plot_one_space(
      kernel_coordinates(details$K_observed_panel),
      "Panel Bayes observado · TRAIN", ctx$group_train, ctx$train_ids
    )
    plot_one_space(
      kernel_coordinates(details$K_reconstructed_panel),
      "Panel Bayes reconstruido por MOFA · TRAIN", ctx$group_train, ctx$train_ids
    )
    plot_one_space(
      details$Z_panel_recovered_train[, seq_len(min(2L, ncol(details$Z_panel_recovered_train))), drop = FALSE],
      "Z recuperado con panel Bayes · TRAIN", ctx$group_train, ctx$train_ids
    )
    if (!is.null(details$Z_full_projection_test)) {
      plot_one_space(
        details$Z_full_projection_test[, seq_len(min(2L, ncol(details$Z_full_projection_test))), drop = FALSE],
        "Proyeccion completa · TEST", ctx$group_test, ctx$test_ids
      )
    } else {
      graphics::plot.new(); graphics::title("TEST completo no disponible")
    }
    if (!is.null(details$Z_panel_recovered_test)) {
      plot_one_space(
        details$Z_panel_recovered_test[, seq_len(min(2L, ncol(details$Z_panel_recovered_test))), drop = FALSE],
        "Z recuperado con panel · TEST", ctx$group_test, ctx$test_ids
      )
    } else {
      graphics::plot.new(); graphics::title("TEST panel no disponible")
    }
    TRUE
  }, error = function(e) {
    warning("No se pudo crear el plot del ganador: ", conditionMessage(e))
    FALSE
  })
}

save_winner_bundle <- function(
    winner_row, comparison_root, output_dir, ridge_grid,
    do_test = TRUE, copy_models = TRUE
) {
  family <- winner_row$family[[1L]]
  root_name <- winner_row$family_root[[1L]]
  model_variant <- winner_row$modelamiento[[1L]]
  analysis <- winner_row$analysis[[1L]]
  shared_dir <- file.path(comparison_root, "upstream_comun", root_name, analysis)
  branch_dir <- file.path(comparison_root, model_variant, root_name, analysis)
  ctx <- build_full_context(shared_dir, analysis, ridge_grid, do_test = do_test)
  panel <- load_panel_for_context(branch_dir, ctx)
  evaluated <- evaluate_panel(ctx, panel, ridge_grid, do_test = do_test, keep_details = TRUE)

  winner_dir <- file.path(output_dir, "ganadores", family, model_variant)
  dir.create(winner_dir, recursive = TRUE, showWarnings = FALSE)
  unique_model_dir <- file.path(output_dir, "modelos_mofa_originales")
  dir.create(unique_model_dir, recursive = TRUE, showWarnings = FALSE)
  copied_model <- file.path(unique_model_dir, paste0(analysis, ".rds"))
  if (copy_models && !file.exists(copied_model)) {
    ok <- file.copy(ctx$mofa_file, copied_model, overwrite = FALSE, copy.mode = TRUE)
    if (!ok) warning("No se pudo copiar el MOFA original: ", ctx$mofa_file)
  }

  write_tsv(winner_row, file.path(winner_dir, "01_metricas_y_ranking_ganador.tsv"))
  write_tsv(panel$panel_table, file.path(winner_dir, "02_features_bayes.tsv"))
  write_tsv(panel$mapping, file.path(winner_dir, "03_mapeo_features_bayes_a_mofa.tsv"))
  write_tsv(ctx$r2_table, file.path(winner_dir, "04a_varianza_mofa_original_usada.tsv"))
  if (!is.null(ctx$r2_table_api)) {
    write_tsv(ctx$r2_table_api, file.path(winner_dir, "04b_varianza_mofa_api.tsv"))
  }
  write_tsv(
    ctx$r2_table_recomputed,
    file.path(winner_dir, "04c_varianza_mofa_recalculada_X_W_Z.tsv")
  )
  write_tsv(
    full_weights_table(ctx$W, panel$selected_ids),
    file.path(winner_dir, "04d_loadings_mofa_original_todas_features.tsv")
  )
  write_tsv(evaluated$view_table, file.path(winner_dir, "05_metricas_por_vista.tsv"))
  write_tsv(evaluated$coverage_table, file.path(winner_dir, "06_cobertura_loading_factor_vista.tsv"))
  write_tsv(evaluated$ridge_table, file.path(winner_dir, "07_sensibilidad_regularizacion.tsv"))

  train_scores <- scores_table(
    ctx$train_ids,
    ctx$group_train,
    list(
      MOFA_original = ctx$Z,
      MOFA_proyeccion_completa = evaluated$details$Z_full_projection_train,
      Bayes_Z_recuperado = evaluated$details$Z_panel_recovered_train,
      Bayes_observado_geometria = kernel_coordinates(evaluated$details$K_observed_panel),
      Bayes_reconstruido_geometria = kernel_coordinates(evaluated$details$K_reconstructed_panel)
    )
  )
  write_tsv(train_scores, file.path(winner_dir, "08_espacios_latentes_train.tsv"))

  if (!is.null(evaluated$details$Z_full_projection_test)) {
    test_scores <- scores_table(
      ctx$test_ids,
      ctx$group_test,
      list(
        MOFA_proyeccion_completa = evaluated$details$Z_full_projection_test,
        Bayes_Z_recuperado = evaluated$details$Z_panel_recovered_test
      )
    )
    write_tsv(test_scores, file.path(winner_dir, "09_espacios_latentes_test.tsv"))
  }

  test_metric_names <- grep("^test_", names(evaluated$metrics), value = TRUE)
  test_metric_values <- evaluated$metrics[test_metric_names]
  test_metric_table <- as.data.frame(test_metric_values, stringsAsFactors = FALSE)
  test_metric_table$test_projection_error <- ctx$test_error
  write_tsv(
    test_metric_table,
    file.path(winner_dir, "09b_metricas_validacion_test.tsv")
  )

  bundle <- list(
    script_version = SCRIPT_VERSION,
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    family = family,
    modelamiento = model_variant,
    analysis = analysis,
    source_mofa_file = ctx$mofa_file,
    copied_mofa_file = if (copy_models) copied_model else NA_character_,
    source_ready_for_modeling = ctx$ready_file,
    source_panel_file = panel$panel_file,
    test_projection_error = ctx$test_error,
    ranking_row = winner_row,
    full_quality = ctx$full_quality,
    variance_source = ctx$variance_source,
    api_variance_explained = ctx$api_variance,
    r2_used = ctx$r2_table,
    r2_api = ctx$r2_table_api,
    r2_recomputed = ctx$r2_table_recomputed,
    panel_evaluation = evaluated
  )
  saveRDS(bundle, file.path(winner_dir, "10_resultado_completo_ganador.rds"), compress = "xz")
  writeLines(
    c(
      paste0("MOFA original fuente: ", ctx$mofa_file),
      paste0("MOFA original copiado: ", if (copy_models) copied_model else "NO SOLICITADO"),
      paste0("ready_for_modeling fuente: ", ctx$ready_file),
      paste0("panel Bayes fuente: ", panel$panel_file),
      paste0(
        "estado proyeccion TEST: ",
        if (is.na(ctx$test_error)) "sin error tecnico" else ctx$test_error
      )
    ),
    file.path(winner_dir, "11_procedencia_archivos.txt"),
    useBytes = TRUE
  )
  save_winner_plot(ctx, evaluated$details, file.path(winner_dir, "12_comparacion_espacios_latentes.png"))

  invisible(list(ctx = ctx, panel = panel, evaluated = evaluated, winner_dir = winner_dir))
}

methodology_text <- function(args, project_root, comparison_root, output_dir) {
  c(
    "SELECCION DE ESPACIOS LATENTES · METODO Y CONFIGURACION",
    paste0("Version del script: ", SCRIPT_VERSION),
    paste0("Fecha: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("PROJECT_ROOT: ", project_root),
    paste0("COMPARISON_ROOT: ", comparison_root),
    paste0("OUTPUT_DIR: ", output_dir),
    paste0("Familias: ", paste(args$families, collapse = ", ")),
    paste0("Modelamientos: ", paste(args$models, collapse = ", ")),
    paste0("Ridge grid: ", paste(format(args$ridge_grid, scientific = TRUE), collapse = ", ")),
    paste0("Validacion TEST posterior: ", args$do_test),
    "",
    "FUENTE DE CANDIDATOS",
    "Se escanean directamente las carpetas canonicas bajo upstream_comun y cada modelamiento.",
    "Las plantillas y resumen_corridas.tsv no definen ni duplican candidatos.",
    "",
    "SELECCION",
    "La seleccion usa solo TRAIN. TEST se proyecta despues de congelar los seis ganadores.",
    "q00 compite; su similitud con q00 se informa pero no entra al ranking.",
    "Las cuatro vistas reciben el mismo peso mediante kernels centrados y normalizados.",
    "Se usan todos los factores presentes en el MOFA seleccionado de cada escenario.",
    "No se exige que cada factor sea compartido: los factores privados se conservan si tienen soporte y son recuperables.",
    "La varianza explicada oficial de MOFA2 tiene prioridad; X/W/Z se recalcula como auditoria y respaldo.",
    "Rango y condicion son continuos. Solo faltantes, no finitos o una construccion imposible invalidan.",
    "CKA lineal compara matrices de similitud entre individuos y no exige alinear Factor1 con PC1.",
    "Z del panel se recupera por inversion SVD-Tikhonov de W restringido, con perdida media igual por vista.",
    "",
    "DOMINIOS",
    "1) Calidad del MOFA original completo.",
    "2) Recuperacion observada del panel Bayes.",
    "3) Recuperacion integrada y de Z mediante W restringido al panel.",
    "4) Cobertura biologica efectiva de factores y vistas.",
    "",
    "ORDEN",
    "Pareto -> maximin -> mediana -> estabilidad LOO -> cobertura factor/vista ->",
    "estabilidad numerica -> condicion -> mayor panel mapeado/listado solo si todo lo anterior empata.",
    "El numero de features no forma parte de ningun dominio."
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (!requireNamespace("MOFA2", quietly = TRUE)) {
    stop("Se requiere MOFA2 en el entorno R que genero los modelos.")
  }

  project_root <- normalize_existing(args$project_root, "PROJECT_ROOT")
  comparison_root <- args$comparison_root
  if (!grepl("^/", comparison_root)) comparison_root <- file.path(project_root, comparison_root)
  comparison_root <- normalize_existing(comparison_root, "COMPARISON_ROOT")
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- args$output_dir %||% file.path(
    project_root, paste0("seleccion_espacios_latentes_", timestamp)
  )
  if (!grepl("^/", output_dir)) output_dir <- file.path(project_root, output_dir)
  output_dir <- normalize_maybe(output_dir)
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = FALSE))) {
    stop("OUTPUT_DIR ya existe y no esta vacio: ", output_dir)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("PROJECT_ROOT    : ", project_root)
  message("COMPARISON_ROOT : ", comparison_root)
  message("OUTPUT_DIR      : ", output_dir)
  writeLines(
    methodology_text(args, project_root, comparison_root, output_dir),
    file.path(output_dir, "00_metodologia_y_configuracion.txt"),
    useBytes = TRUE
  )

  registry <- build_registry(comparison_root, args$families, args$models)
  registry_out <- registry
  path_cols <- c("shared_dir", "branch_dir", "ready_file", "mofa_file", "panel_csv", "diagnostic_objects")
  for (nm in intersect(path_cols, names(registry_out))) {
    registry_out[[nm]] <- vapply(
      registry_out[[nm]],
      function(x) if (is.na(x)) NA_character_ else relative_to(x, project_root),
      character(1)
    )
  }
  write_tsv(registry_out, file.path(output_dir, "01_inventario_candidatos.tsv"))

  result_rows <- list()
  view_rows <- list()
  coverage_rows <- list()
  ridge_rows <- list()
  q00_rows <- list()
  controls <- list()
  result_i <- view_i <- coverage_i <- ridge_i <- q_i <- control_i <- 0L

  add_control <- function(family, modelamiento, analysis, severity, message_text) {
    control_i <<- control_i + 1L
    controls[[control_i]] <<- data.frame(
      family = family,
      modelamiento = modelamiento,
      analysis = analysis,
      severity = severity,
      message = message_text,
      stringsAsFactors = FALSE
    )
  }

  for (family in args$families) {
    family_registry <- registry[registry$family == family, , drop = FALSE]
    analyses <- unique(family_registry$analysis)
    q00_names <- analyses[is_q00(analyses)]
    if (length(q00_names) != 1L) {
      stop("La familia ", family, " debe contener exactamente un q00; encontrados: ", length(q00_names))
    }
    q00_name <- q00_names[[1L]]
    q00_reg <- family_registry[family_registry$analysis == q00_name, , drop = FALSE][1L, ]
    message("\n[", family, "] Cargando referencia q00: ", q00_name)
    qctx <- tryCatch(
      build_full_context(q00_reg$shared_dir, q00_name, args$ridge_grid, do_test = FALSE),
      error = function(e) e
    )
    if (inherits(qctx, "error")) stop("No se pudo construir q00 para ", family, ": ", conditionMessage(qctx))

    for (j in seq_along(analyses)) {
      analysis <- analyses[[j]]
      message("[", family, " ", j, "/", length(analyses), "] ", analysis)
      rows_analysis <- family_registry[family_registry$analysis == analysis, , drop = FALSE]
      if (analysis == q00_name) {
        ctx <- qctx
      } else {
        ctx <- tryCatch(
          build_full_context(rows_analysis$shared_dir[[1L]], analysis, args$ridge_grid, do_test = FALSE),
          error = function(e) e
        )
      }

      if (inherits(ctx, "error")) {
        msg <- conditionMessage(ctx)
        for (r in seq_len(nrow(rows_analysis))) {
          result_i <- result_i + 1L
          result_rows[[result_i]] <- empty_result(rows_analysis[r, ], msg)
          add_control(family, rows_analysis$modelamiento[[r]], analysis, "CRITICAL", msg)
        }
        next
      }

      if (!setequal(ctx$train_ids, qctx$train_ids)) {
        msg <- "El conjunto de IDs TRAIN no coincide con el q00 de su familia."
        for (r in seq_len(nrow(rows_analysis))) {
          result_i <- result_i + 1L
          result_rows[[result_i]] <- empty_result(rows_analysis[r, ], msg)
          add_control(family, rows_analysis$modelamiento[[r]], analysis, "CRITICAL", msg)
        }
        if (analysis != q00_name) {
          rm(ctx)
          invisible(gc(verbose = FALSE))
        }
        next
      }

      qcmp <- compare_to_q00(ctx, qctx)
      q_i <- q_i + 1L
      q00_rows[[q_i]] <- data.frame(
        family = family,
        family_root = rows_analysis$family_root[[1L]],
        analysis = analysis,
        q00 = is_q00(analysis),
        as.data.frame(flatten_named_list(qcmp), stringsAsFactors = FALSE),
        stringsAsFactors = FALSE
      )

      for (r in seq_len(nrow(rows_analysis))) {
        reg <- rows_analysis[r, ]
        model_variant <- reg$modelamiento[[1L]]
        if (!isTRUE(reg$canonical_prefix_ok[[1L]])) {
          msg <- "Prefijo de escenario incompatible con la familia canonica."
          result_i <- result_i + 1L
          result_rows[[result_i]] <- empty_result(reg, msg)
          add_control(family, model_variant, analysis, "CRITICAL", msg)
          next
        }
        if (!isTRUE(reg$branch_exists[[1L]])) {
          msg <- "No existe la rama del modelamiento."
          result_i <- result_i + 1L
          result_rows[[result_i]] <- empty_result(reg, msg)
          add_control(family, model_variant, analysis, "CRITICAL", msg)
          next
        }

        evaluated <- tryCatch({
          panel <- load_panel_for_context(reg$branch_dir[[1L]], ctx)
          ev <- evaluate_panel(ctx, panel, args$ridge_grid, do_test = FALSE, keep_details = FALSE)
          list(panel = panel, evaluation = ev)
        }, error = function(e) e)

        if (inherits(evaluated, "error")) {
          msg <- conditionMessage(evaluated)
          result_i <- result_i + 1L
          result_rows[[result_i]] <- empty_result(reg, msg)
          add_control(family, model_variant, analysis, "CRITICAL", msg)
          next
        }

        base_row <- list(
          family = family,
          family_root = reg$family_root[[1L]],
          modelamiento = model_variant,
          analysis = analysis,
          q00 = is_q00(analysis),
          technical_valid = TRUE,
          technical_error = NA_character_,
          variance_source = ctx$variance_source,
          mofa_file = relative_to(reg$mofa_file[[1L]], project_root),
          ready_file = relative_to(reg$ready_file[[1L]], project_root),
          panel_file = relative_to(reg$panel_csv[[1L]], project_root)
        )
        full_metrics <- ctx$full_quality$metrics
        panel_metrics <- evaluated$evaluation$metrics
        result_i <- result_i + 1L
        result_rows[[result_i]] <- as.data.frame(
          c(base_row, flatten_named_list(full_metrics), flatten_named_list(panel_metrics)),
          stringsAsFactors = FALSE
        )

        vtbl <- evaluated$evaluation$view_table
        vtbl$family <- family
        vtbl$modelamiento <- model_variant
        vtbl$analysis <- analysis
        view_i <- view_i + 1L
        view_rows[[view_i]] <- vtbl[, c("family", "modelamiento", "analysis", setdiff(names(vtbl), c("family", "modelamiento", "analysis"))), drop = FALSE]

        ctbl <- evaluated$evaluation$coverage_table
        ctbl$family <- family
        ctbl$modelamiento <- model_variant
        ctbl$analysis <- analysis
        coverage_i <- coverage_i + 1L
        coverage_rows[[coverage_i]] <- ctbl[, c("family", "modelamiento", "analysis", setdiff(names(ctbl), c("family", "modelamiento", "analysis"))), drop = FALSE]

        rtbl <- evaluated$evaluation$ridge_table
        rtbl$family <- family
        rtbl$modelamiento <- model_variant
        rtbl$analysis <- analysis
        ridge_i <- ridge_i + 1L
        ridge_rows[[ridge_i]] <- rtbl[, c("family", "modelamiento", "analysis", setdiff(names(rtbl), c("family", "modelamiento", "analysis"))), drop = FALSE]

        if (evaluated$panel$mapping_fraction < 1) {
          add_control(
            family, model_variant, analysis, "MEDIUM",
            paste0("Mapeo parcial de features Bayes: ", signif(evaluated$panel$mapping_fraction, 4))
          )
        }
        if (!identical(ctx$variance_source, "MOFA2::get_variance_explained")) {
          add_control(
            family, model_variant, analysis, "LOW",
            "No se pudo normalizar la salida oficial de varianza MOFA2; se uso el calculo reproducible X/W/Z."
          )
        }
        if (ctx$full_quality$metrics$full_z_rank_fraction < 1) {
          add_control(
            family, model_variant, analysis, "MEDIUM",
            "Z original no alcanza rango numerico completo; se conserva como penalizacion continua."
          )
        }
        panel_metric_now <- evaluated$evaluation$metrics
        if (panel_metric_now$panel_joint_rank_fraction < 1 ||
            panel_metric_now$panel_view_rank_min < 1 ||
            !is.finite(panel_metric_now$panel_joint_condition)) {
          add_control(
            family, model_variant, analysis, "MEDIUM",
            paste0(
              "Panel con representacion/rango incompleto o condicion no finita: rango conjunto=",
              signif(panel_metric_now$panel_joint_rank_fraction, 4),
              "; rango minimo por vista=", signif(panel_metric_now$panel_view_rank_min, 4),
              "; condicion=", format(panel_metric_now$panel_joint_condition, digits = 5),
              ". Se mantiene y se penaliza numericamente."
            )
          )
        }
      }

      if (analysis != q00_name) {
        rm(ctx)
        invisible(gc(verbose = FALSE))
      }
    }
    rm(qctx)
    invisible(gc(verbose = FALSE))
  }

  raw_metrics <- bind_rows_fill(result_rows)
  numeric_candidates <- setdiff(
    names(raw_metrics),
    c(
      "family", "family_root", "modelamiento", "analysis", "q00",
      "technical_valid", "technical_error", "variance_source",
      "mofa_file", "ready_file", "panel_file"
    )
  )
  for (nm in numeric_candidates) raw_metrics[[nm]] <- suppressWarnings(as.numeric(raw_metrics[[nm]]))
  raw_metrics$q00 <- as.logical(raw_metrics$q00)
  raw_metrics$technical_valid <- as.logical(raw_metrics$technical_valid)

  write_tsv(raw_metrics, file.path(output_dir, "02_metricas_escenario_modelo.tsv"))
  write_tsv(bind_rows_fill(view_rows), file.path(output_dir, "03_metricas_por_vista.tsv"))
  write_tsv(bind_rows_fill(coverage_rows), file.path(output_dir, "04_cobertura_loading_factor_vista.tsv"))
  write_tsv(bind_rows_fill(ridge_rows), file.path(output_dir, "05_sensibilidad_regularizacion.tsv"))
  write_tsv(bind_rows_fill(q00_rows), file.path(output_dir, "06_sensibilidad_q00.tsv"))
  write_tsv(bind_rows_fill(controls), file.path(output_dir, "07_controles_tecnicos.tsv"))

  ranked <- rank_all_competitions(raw_metrics)
  write_tsv(ranked, file.path(output_dir, "08_rankings_completos.tsv"))
  winners <- ranked[ranked$winner %in% TRUE, , drop = FALSE]
  winners <- winners[order(winners$family, winners$modelamiento), , drop = FALSE]
  write_tsv(winners, file.path(output_dir, "09_seis_ganadores.tsv"))

  expected_winners <- length(args$families) * length(args$models)
  if (nrow(winners) != expected_winners) {
    stop(
      "Se esperaban ", expected_winners, " ganadores y se obtuvieron ", nrow(winners),
      ". Revise 07_controles_tecnicos.tsv y 08_rankings_completos.tsv."
    )
  }

  message("\nGuardando objetos y espacios de los ganadores...")
  for (i in seq_len(nrow(winners))) {
    message(
      "  [", i, "/", nrow(winners), "] ", winners$family[[i]], " / ",
      winners$modelamiento[[i]], " / ", winners$analysis[[i]]
    )
    save_winner_bundle(
      winners[i, , drop = FALSE], comparison_root, output_dir,
      args$ridge_grid, do_test = args$do_test, copy_models = args$copy_models
    )
  }

  writeLines(capture.output(sessionInfo()), file.path(output_dir, "10_session_info.txt"), useBytes = TRUE)
  summary_lines <- c(
    "SELECCION DE ESPACIOS LATENTES COMPLETADA",
    paste0("Version: ", SCRIPT_VERSION),
    paste0("Fecha: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("Candidatos evaluados: ", nrow(raw_metrics)),
    paste0("Candidatos tecnicamente validos: ", sum(raw_metrics$technical_valid %in% TRUE)),
    paste0("Ganadores: ", nrow(winners)),
    "",
    paste0(
      winners$family, " / ", winners$modelamiento, " -> ", winners$analysis,
      " | Pareto=", winners$pareto_front,
      " | maximin=", format(winners$domain_minimax, digits = 4),
      " | mediana=", format(winners$domain_median, digits = 4)
    ),
    "",
    "TEST no participo en el ranking. Sus proyecciones aparecen solo dentro de los paquetes ganadores.",
    "El numero de features mapeadas/listadas solo se uso como ultimo desempate determinista."
  )
  writeLines(summary_lines, file.path(output_dir, "RESUMEN_SELECCION.txt"), useBytes = TRUE)
  cat(paste(summary_lines, collapse = "\n"), "\n")
  cat("\nSalida completa:\n", output_dir, "\n", sep = "")
  invisible(list(registry = registry, raw = raw_metrics, ranking = ranked, winners = winners))
}

if (sys.nframe() == 0L) main()

#!/usr/bin/env Rscript

## ============================================================================
## PCA DE ESPACIOS OBSERVADO, RECONSTRUIDO Y LATENTE
##
## Para CADA escenario x modelamiento:
##   1) Recupera exactamente las features Bayes usadas por el selector.
##   2) PCA observado:
##        - usa solo las features seleccionadas
##        - ajusta centrado/escalado SOLO en TRAIN
##        - balancea vistas dividiendo cada bloque por sqrt(numero de features)
##        - ajusta PCA en TRAIN
##        - proyecta TEST con los parametros congelados de TRAIN
##   3) PCA reconstruido:
##        - TRAIN = W_panel %*% t(Z_MOFA_train)
##        - TEST  = W_panel %*% t(Z_panel_test_recuperado)
##        - el Z de TEST se recupera exclusivamente desde las features
##          seleccionadas del panel mediante la misma proyeccion SVD-Tikhonov
##          usada por seleccionar_mejores_espacios_latentes.R
##        - PCA ajustado en TRAIN y TEST proyectado posteriormente
##   4) PCA de factores:
##        - TRAIN = Z original del MOFA
##        - TEST  = Z recuperado desde el panel seleccionado
##        - PCA ajustado exclusivamente sobre Z TRAIN
##        - TEST proyectado con predict()
##
## IMPORTANTE:
##   * TEST nunca se usa para centrar, escalar, elegir variables ni ajustar PCA.
##   * Este script NO reentrena MOFA ni modelos Bayes.
##   * Reutiliza las funciones auditadas del script de seleccion.
##
## Uso:
##   Rscript pca_espacios_train_project_test.R \
##     --project-root . \
##     --selection-script seleccionar_mejores_espacios_latentes.R
##
## Opciones:
##   --project-root RUTA
##   --comparison-root RUTA
##   --selection-script RUTA
##   --output-dir RUTA
##   --families TODAS,CON_FLI,SIN_FLI
##   --models modelamiento2,modelamiento3
##   --ridge-rel 1e-6
## ============================================================================

options(stringsAsFactors = FALSE, width = 220)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

parse_csv_arg <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1L]])

parse_args_pca <- function(args) {
  out <- list(
    project_root = ".",
    comparison_root = "comparacion_modelamientos",
    selection_script = "seleccionar_mejores_espacios_latentes.R",
    output_dir = NULL,
    families = c("TODAS", "CON_FLI", "SIN_FLI"),
    models = c("modelamiento2", "modelamiento3"),
    ridge_rel = 1e-6
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--help") {
      cat(
        "Uso:\n",
        "  Rscript pca_espacios_train_project_test.R [opciones]\n\n",
        "Opciones:\n",
        "  --project-root RUTA\n",
        "  --comparison-root RUTA\n",
        "  --selection-script RUTA\n",
        "  --output-dir RUTA\n",
        "  --families TODAS,CON_FLI,SIN_FLI\n",
        "  --models modelamiento2,modelamiento3\n",
        "  --ridge-rel 1e-6\n",
        sep = ""
      )
      quit(save = "no", status = 0L)
    }

    allowed <- c(
      "--project-root", "--comparison-root", "--selection-script",
      "--output-dir", "--families", "--models", "--ridge-rel"
    )
    if (!key %in% allowed) stop("Argumento no reconocido: ", key)
    if (i == length(args)) stop("Falta valor para: ", key)
    value <- args[[i + 1L]]

    if (key == "--project-root") out$project_root <- value
    if (key == "--comparison-root") out$comparison_root <- value
    if (key == "--selection-script") out$selection_script <- value
    if (key == "--output-dir") out$output_dir <- value
    if (key == "--families") out$families <- toupper(parse_csv_arg(value))
    if (key == "--models") out$models <- tolower(parse_csv_arg(value))
    if (key == "--ridge-rel") out$ridge_rel <- suppressWarnings(as.numeric(value))

    i <- i + 2L
  }

  if (!all(out$families %in% c("TODAS", "CON_FLI", "SIN_FLI"))) {
    stop("--families admite TODAS,CON_FLI,SIN_FLI")
  }
  if (!all(out$models %in% c("modelamiento2", "modelamiento3"))) {
    stop("--models admite modelamiento2,modelamiento3")
  }
  if (length(out$ridge_rel) != 1L || !is.finite(out$ridge_rel) || out$ridge_rel <= 0) {
    stop("--ridge-rel debe ser un numero positivo y finito.")
  }

  out
}

resolve_path <- function(path, root, must_work = TRUE) {
  if (!grepl("^/", path)) path <- file.path(root, path)
  normalizePath(path, winslash = "/", mustWork = must_work)
}

## ----------------------------------------------------------------------------
## Estandarizacion TRAIN -> TEST
## ----------------------------------------------------------------------------

fit_train_feature_transform <- function(X_train, X_test, view_name) {
  ## Entrada: features x samples
  X_train <- as.matrix(X_train)
  X_test <- as.matrix(X_test)
  storage.mode(X_train) <- "double"
  storage.mode(X_test) <- "double"

  if (!identical(rownames(X_train), rownames(X_test))) {
    if (!setequal(rownames(X_train), rownames(X_test))) {
      stop("TRAIN/TEST no contienen las mismas features en vista ", view_name)
    }
    X_test <- X_test[rownames(X_train), , drop = FALSE]
  }

  mu <- rowMeans(X_train)
  sdv <- apply(X_train, 1L, stats::sd)
  keep <- is.finite(mu) & is.finite(sdv) & sdv > sqrt(.Machine$double.eps)

  if (!any(keep)) {
    stop("No quedan features con varianza TRAIN positiva en vista ", view_name)
  }

  Xtr <- X_train[keep, , drop = FALSE]
  Xte <- X_test[keep, , drop = FALSE]
  mu <- mu[keep]
  sdv <- sdv[keep]

  Xtr <- sweep(Xtr, 1L, mu, "-")
  Xtr <- sweep(Xtr, 1L, sdv, "/")
  Xte <- sweep(Xte, 1L, mu, "-")
  Xte <- sweep(Xte, 1L, sdv, "/")

  ## Peso aproximadamente igual por vista:
  ## despues de z-score, Frobenius esperado ~ p_v * n.
  ## Dividir por sqrt(p_v) evita que la vista con mas variables domine.
  block_scale <- sqrt(nrow(Xtr))
  Xtr <- Xtr / block_scale
  Xte <- Xte / block_scale

  list(
    train = Xtr,
    test = Xte,
    kept_features = rownames(Xtr),
    center = mu,
    scale = sdv,
    block_scale = block_scale
  )
}

stack_views_for_pca <- function(train_list, test_list) {
  if (!all(CANONICAL_VIEWS %in% names(train_list)) ||
      !all(CANONICAL_VIEWS %in% names(test_list))) {
    stop("Faltan vistas canonicas.")
  }

  transforms <- list()
  train_blocks <- list()
  test_blocks <- list()

  for (v in CANONICAL_VIEWS) {
    tr <- fit_train_feature_transform(train_list[[v]], test_list[[v]], v)
    transforms[[v]] <- tr

    A <- tr$train
    B <- tr$test

    rownames(A) <- paste0(v, "::", rownames(A))
    rownames(B) <- paste0(v, "::", rownames(B))

    train_blocks[[v]] <- A
    test_blocks[[v]] <- B
  }

  ## PCA requiere samples x features
  train_mat <- t(do.call(rbind, train_blocks))
  test_mat <- t(do.call(rbind, test_blocks))

  if (any(!is.finite(train_mat)) || any(!is.finite(test_mat))) {
    stop("Matriz PCA contiene valores no finitos tras transformacion TRAIN.")
  }

  list(
    train = train_mat,
    test = test_mat,
    transforms = transforms
  )
}

fit_pca_train_project_test <- function(train_mat, test_mat, ncomp = 2L) {
  train_mat <- as.matrix(train_mat)
  test_mat <- as.matrix(test_mat)

  if (!identical(colnames(train_mat), colnames(test_mat))) {
    if (!setequal(colnames(train_mat), colnames(test_mat))) {
      stop("TRAIN y TEST no tienen las mismas variables para PCA.")
    }
    test_mat <- test_mat[, colnames(train_mat), drop = FALSE]
  }

  max_rank <- min(nrow(train_mat) - 1L, ncol(train_mat))
  if (max_rank < 1L) stop("Rango insuficiente para PCA.")

  ## Ya se centro/escalo cada feature con TRAIN; no se repite aqui.
  pca <- stats::prcomp(
    train_mat,
    center = FALSE,
    scale. = FALSE,
    rank. = min(max(ncomp, 2L), max_rank)
  )

  tr_scores <- pca$x
  te_scores <- stats::predict(pca, newdata = test_mat)

  if (ncol(tr_scores) < 2L) {
    tr_scores <- cbind(tr_scores, PC2 = 0)
    te_scores <- cbind(te_scores, PC2 = 0)
  }

  var_explained <- pca$sdev^2 / sum(pca$sdev^2)

  list(
    pca = pca,
    train_scores = tr_scores,
    test_scores = te_scores,
    variance_explained = var_explained
  )
}

fit_pca_latent_train_project_test <- function(Z_train, Z_test, ncomp = 2L) {
  Z_train <- as.matrix(Z_train)
  Z_test <- as.matrix(Z_test)
  storage.mode(Z_train) <- "double"
  storage.mode(Z_test) <- "double"

  if (!identical(colnames(Z_train), colnames(Z_test))) {
    if (!setequal(colnames(Z_train), colnames(Z_test))) {
      stop("Z TRAIN y Z TEST no tienen los mismos factores.")
    }
    Z_test <- Z_test[, colnames(Z_train), drop = FALSE]
  }

  if (any(!is.finite(Z_train)) || any(!is.finite(Z_test))) {
    stop("Z contiene valores no finitos.")
  }

  ## Aqui SI se centra/escala con TRAIN porque los factores pueden tener
  ## magnitudes distintas. predict() aplica esos mismos parametros a TEST.
  pca <- stats::prcomp(
    Z_train,
    center = TRUE,
    scale. = TRUE,
    rank. = min(max(ncomp, 2L), min(nrow(Z_train) - 1L, ncol(Z_train)))
  )

  tr_scores <- pca$x
  te_scores <- stats::predict(pca, newdata = Z_test)

  if (ncol(tr_scores) < 2L) {
    tr_scores <- cbind(tr_scores, PC2 = 0)
    te_scores <- cbind(te_scores, PC2 = 0)
  }

  var_explained <- pca$sdev^2 / sum(pca$sdev^2)

  list(
    pca = pca,
    train_scores = tr_scores,
    test_scores = te_scores,
    variance_explained = var_explained
  )
}

## ----------------------------------------------------------------------------
## Reconstruccion de features seleccionadas
## ----------------------------------------------------------------------------

reconstruct_selected <- function(W_panel, Z_samples) {
  out <- lapply(CANONICAL_VIEWS, function(v) {
    Wv <- as.matrix(W_panel[[v]])
    if (!is.null(colnames(Wv)) && setequal(colnames(Wv), colnames(Z_samples))) {
      Wv <- Wv[, colnames(Z_samples), drop = FALSE]
    } else {
      if (ncol(Wv) != ncol(Z_samples)) {
        stop("Numero de factores incompatible para reconstruccion en ", v)
      }
      colnames(Wv) <- colnames(Z_samples)
    }
    Xhat <- Wv %*% t(Z_samples)  # features x samples
    colnames(Xhat) <- rownames(Z_samples)
    rownames(Xhat) <- rownames(Wv)
    Xhat
  })
  names(out) <- CANONICAL_VIEWS
  out
}

## ----------------------------------------------------------------------------
## Tablas y plots
## ----------------------------------------------------------------------------

scores_to_df <- function(scores, ids, groups, split, space) {
  scores <- as.matrix(scores)
  if (!is.null(rownames(scores)) && setequal(rownames(scores), ids)) {
    scores <- scores[ids, , drop = FALSE]
  }
  data.frame(
    sample_id = ids,
    group = as.character(groups),
    split = split,
    space = space,
    PC1 = scores[, 1L],
    PC2 = scores[, 2L],
    stringsAsFactors = FALSE
  )
}

variance_row <- function(obj, space) {
  ve <- obj$variance_explained
  data.frame(
    space = space,
    pc1_variance_fraction = if (length(ve) >= 1L) ve[[1L]] else NA_real_,
    pc2_variance_fraction = if (length(ve) >= 2L) ve[[2L]] else NA_real_,
    stringsAsFactors = FALSE
  )
}

plot_space_base <- function(df, title, xlab, ylab, color_map) {
  groups <- as.character(df$group)
  cols <- unname(color_map[groups])
  cols[is.na(cols)] <- "grey50"
  pch <- ifelse(df$split == "TRAIN", 19, 17)

  graphics::plot(
    df$PC1, df$PC2,
    col = cols, pch = pch, cex = 1.1,
    xlab = xlab, ylab = ylab,
    main = title
  )
  graphics::abline(h = 0, v = 0, lty = 3, col = "grey80")
}

save_three_panel_plot <- function(all_scores, variance_tbl, out_png, subtitle) {
  groups <- sort(unique(as.character(all_scores$group)))
  groups[is.na(groups) | !nzchar(groups)] <- "NA"
  palette <- grDevices::hcl.colors(max(length(groups), 1L), palette = "Dark 3")
  color_map <- setNames(palette, groups)

  grDevices::png(out_png, width = 3300, height = 1200, res = 220)
  on.exit(grDevices::dev.off(), add = TRUE)

  graphics::par(mfrow = c(1, 3), mar = c(4.5, 4.5, 4.2, 1.2), oma = c(0, 0, 2.5, 0))

  spaces <- c("OBSERVED_SELECTED", "RECONSTRUCTED_SELECTED", "LATENT_FACTORS")
  titles <- c(
    "Observado · features seleccionadas",
    "Reconstruido · features seleccionadas",
    "Factores latentes"
  )

  for (i in seq_along(spaces)) {
    sp <- spaces[[i]]
    d <- all_scores[all_scores$space == sp, , drop = FALSE]
    v <- variance_tbl[variance_tbl$space == sp, , drop = FALSE]

    pc1 <- if (nrow(v)) 100 * v$pc1_variance_fraction[[1L]] else NA_real_
    pc2 <- if (nrow(v)) 100 * v$pc2_variance_fraction[[1L]] else NA_real_

    xlab <- if (is.finite(pc1)) sprintf("PC1 (%.1f%%)", pc1) else "PC1"
    ylab <- if (is.finite(pc2)) sprintf("PC2 (%.1f%%)", pc2) else "PC2"

    plot_space_base(d, titles[[i]], xlab, ylab, color_map)

    if (i == 1L) {
      graphics::legend(
        "topright",
        legend = names(color_map),
        col = unname(color_map),
        pch = 19,
        cex = 0.72,
        bty = "n",
        title = "Grupo"
      )
      graphics::legend(
        "bottomright",
        legend = c("TRAIN", "TEST"),
        col = "black",
        pch = c(19, 17),
        cex = 0.72,
        bty = "n",
        title = "Split"
      )
    }
  }

  graphics::mtext(subtitle, outer = TRUE, side = 3, line = 0.6, cex = 1.05, font = 2)
}

safe_dir_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

## ----------------------------------------------------------------------------
## Analisis de un candidato
## ----------------------------------------------------------------------------

analyse_candidate_pca <- function(reg, comparison_root, output_dir, ridge_rel) {
  family <- reg$family[[1L]]
  model_variant <- reg$modelamiento[[1L]]
  analysis <- reg$analysis[[1L]]

  shared_dir <- reg$shared_dir[[1L]]
  branch_dir <- reg$branch_dir[[1L]]

  if (!isTRUE(reg$branch_exists[[1L]])) {
    stop("No existe rama de modelamiento: ", branch_dir)
  }

  ctx <- build_full_context(
    shared_dir = shared_dir,
    analysis = analysis,
    ridge_grid = ridge_rel,
    do_test = TRUE
  )

  if (is.null(ctx$X_test)) {
    stop("TEST no disponible/proyectable: ", ctx$test_error %||% "sin detalle")
  }

  panel <- load_panel_for_context(branch_dir, ctx)
  if (is.null(panel$X_test)) stop("Panel TEST no disponible.")

  ## --------------------------------------------------------------------------
  ## 1) OBSERVADO
  ## --------------------------------------------------------------------------
  observed_stack <- stack_views_for_pca(panel$X_train, panel$X_test)
  pca_observed <- fit_pca_train_project_test(
    observed_stack$train,
    observed_stack$test,
    ncomp = 2L
  )

  ## --------------------------------------------------------------------------
  ## Recuperacion de Z TEST exclusivamente desde las features seleccionadas.
  ## Tambien calculamos Z recuperado TRAIN como diagnostico, aunque el PCA
  ## latente solicitado se referencia contra el Z original de TRAIN.
  ## --------------------------------------------------------------------------
  Z_panel_train <- project_latent_equal_view(
    panel$X_train, panel$W, ridge_rel
  )$Z

  Z_panel_test <- project_latent_equal_view(
    panel$X_test, panel$W, ridge_rel
  )$Z

  ## Alinear nombres de factores con Z original
  if (ncol(Z_panel_train) == ncol(ctx$Z)) colnames(Z_panel_train) <- colnames(ctx$Z)
  if (ncol(Z_panel_test) == ncol(ctx$Z)) colnames(Z_panel_test) <- colnames(ctx$Z)

  ## --------------------------------------------------------------------------
  ## 2) RECONSTRUIDO
  ##
  ## TRAIN: reconstruccion de las features seleccionadas usando Z original.
  ## TEST:  reconstruccion usando Z recuperado del panel TEST.
  ## --------------------------------------------------------------------------
  reconstructed_train <- reconstruct_selected(panel$W, ctx$Z)
  reconstructed_test <- reconstruct_selected(panel$W, Z_panel_test)

  reconstructed_stack <- stack_views_for_pca(
    reconstructed_train,
    reconstructed_test
  )

  pca_reconstructed <- fit_pca_train_project_test(
    reconstructed_stack$train,
    reconstructed_stack$test,
    ncomp = 2L
  )

  ## --------------------------------------------------------------------------
  ## 3) FACTORES LATENTES
  ##
  ## TRAIN: Z original de MOFA.
  ## TEST: Z recuperado desde las features seleccionadas.
  ## --------------------------------------------------------------------------
  pca_latent <- fit_pca_latent_train_project_test(
    ctx$Z,
    Z_panel_test,
    ncomp = 2L
  )

  train_ids <- ctx$train_ids
  test_ids <- ctx$test_ids
  group_train <- ctx$group_train
  group_test <- ctx$group_test

  scores <- rbind(
    scores_to_df(
      pca_observed$train_scores, train_ids, group_train,
      "TRAIN", "OBSERVED_SELECTED"
    ),
    scores_to_df(
      pca_observed$test_scores, test_ids, group_test,
      "TEST", "OBSERVED_SELECTED"
    ),
    scores_to_df(
      pca_reconstructed$train_scores, train_ids, group_train,
      "TRAIN", "RECONSTRUCTED_SELECTED"
    ),
    scores_to_df(
      pca_reconstructed$test_scores, test_ids, group_test,
      "TEST", "RECONSTRUCTED_SELECTED"
    ),
    scores_to_df(
      pca_latent$train_scores, train_ids, group_train,
      "TRAIN", "LATENT_FACTORS"
    ),
    scores_to_df(
      pca_latent$test_scores, test_ids, group_test,
      "TEST", "LATENT_FACTORS"
    )
  )

  variance_tbl <- rbind(
    variance_row(pca_observed, "OBSERVED_SELECTED"),
    variance_row(pca_reconstructed, "RECONSTRUCTED_SELECTED"),
    variance_row(pca_latent, "LATENT_FACTORS")
  )

  n_features_tbl <- data.frame(
    view = CANONICAL_VIEWS,
    n_selected_features = vapply(panel$X_train[CANONICAL_VIEWS], nrow, integer(1)),
    n_observed_features_used_pca = vapply(
      observed_stack$transforms[CANONICAL_VIEWS],
      function(x) length(x$kept_features),
      integer(1)
    ),
    n_reconstructed_features_used_pca = vapply(
      reconstructed_stack$transforms[CANONICAL_VIEWS],
      function(x) length(x$kept_features),
      integer(1)
    ),
    stringsAsFactors = FALSE
  )

  candidate_dir <- file.path(
    output_dir,
    family,
    model_variant,
    safe_dir_name(analysis)
  )
  dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)

  write_tsv(scores, file.path(candidate_dir, "01_pca_scores_train_test.tsv"))
  write_tsv(variance_tbl, file.path(candidate_dir, "02_pca_variance_explained.tsv"))
  write_tsv(n_features_tbl, file.path(candidate_dir, "03_features_used_by_view.tsv"))

  latent_raw <- rbind(
    data.frame(
      sample_id = train_ids,
      group = group_train,
      split = "TRAIN",
      as.data.frame(ctx$Z, stringsAsFactors = FALSE),
      stringsAsFactors = FALSE
    ),
    data.frame(
      sample_id = test_ids,
      group = group_test,
      split = "TEST",
      as.data.frame(Z_panel_test, stringsAsFactors = FALSE),
      stringsAsFactors = FALSE
    )
  )
  write_tsv(latent_raw, file.path(candidate_dir, "04_latent_scores_reference_train_projected_test.tsv"))

  latent_recovery_train <- data.frame(
    sample_id = train_ids,
    group = group_train,
    as.data.frame(ctx$Z, stringsAsFactors = FALSE),
    stringsAsFactors = FALSE
  )
  colnames(latent_recovery_train)[-(1:2)] <- paste0(
    "Z_original__", colnames(ctx$Z)
  )
  zrec <- as.data.frame(Z_panel_train[train_ids, , drop = FALSE], stringsAsFactors = FALSE)
  colnames(zrec) <- paste0("Z_panel_recovered__", colnames(Z_panel_train))
  latent_recovery_train <- cbind(latent_recovery_train, zrec)
  write_tsv(
    latent_recovery_train,
    file.path(candidate_dir, "05_latent_train_original_vs_panel_recovered.tsv")
  )

  saveRDS(
    list(
      family = family,
      modelamiento = model_variant,
      analysis = analysis,
      ridge_rel = ridge_rel,
      selected_feature_ids = panel$selected_ids,
      pca_observed = pca_observed$pca,
      pca_reconstructed = pca_reconstructed$pca,
      pca_latent = pca_latent$pca,
      observed_transform = observed_stack$transforms,
      reconstructed_transform = reconstructed_stack$transforms,
      Z_train_original = ctx$Z,
      Z_train_panel_recovered = Z_panel_train,
      Z_test_panel_recovered = Z_panel_test,
      scores = scores,
      variance_explained = variance_tbl
    ),
    file.path(candidate_dir, "06_pca_objects_and_transforms.rds"),
    compress = "xz"
  )

  save_three_panel_plot(
    all_scores = scores,
    variance_tbl = variance_tbl,
    out_png = file.path(candidate_dir, "07_pca_observed_reconstructed_latent_train_test.png"),
    subtitle = paste(family, model_variant, analysis, sep = " · ")
  )

  data.frame(
    family = family,
    modelamiento = model_variant,
    analysis = analysis,
    q00 = reg$q00[[1L]],
    technical_valid = TRUE,
    technical_error = NA_character_,
    n_train = length(train_ids),
    n_test = length(test_ids),
    n_selected_features = sum(vapply(panel$X_train, nrow, integer(1))),
    mapping_fraction = panel$mapping_fraction,
    ridge_rel = ridge_rel,
    output_dir = candidate_dir,
    stringsAsFactors = FALSE
  )
}

## ----------------------------------------------------------------------------
## MAIN
## ----------------------------------------------------------------------------

main_pca <- function() {
  args <- parse_args_pca(commandArgs(trailingOnly = TRUE))

  project_root <- normalizePath(
    args$project_root,
    winslash = "/",
    mustWork = TRUE
  )

  comparison_root <- resolve_path(
    args$comparison_root,
    project_root,
    must_work = TRUE
  )

  selection_script <- resolve_path(
    args$selection_script,
    project_root,
    must_work = TRUE
  )

  ## Carga las funciones del selector. Su main() no se ejecuta al hacer source().
  source(selection_script, local = .GlobalEnv)

  if (!exists("build_registry", mode = "function") ||
      !exists("build_full_context", mode = "function") ||
      !exists("load_panel_for_context", mode = "function") ||
      !exists("project_latent_equal_view", mode = "function")) {
    stop("El script de seleccion no expone las funciones esperadas.")
  }

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- args$output_dir %||%
    file.path(project_root, paste0("pca_espacios_train_test_", timestamp))
  output_dir <- resolve_path(output_dir, project_root, must_work = FALSE)

  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = FALSE))) {
    stop("OUTPUT_DIR ya existe y no esta vacio: ", output_dir)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  registry <- build_registry(
    comparison_root = comparison_root,
    families = args$families,
    models = args$models
  )

  registry <- registry[
    registry$family %in% args$families &
      registry$modelamiento %in% args$models,
    ,
    drop = FALSE
  ]

  if (!nrow(registry)) stop("No hay candidatos para analizar.")

  message("Candidatos PCA: ", nrow(registry))
  message("Salida: ", output_dir)

  manifest_rows <- vector("list", nrow(registry))

  for (i in seq_len(nrow(registry))) {
    reg <- registry[i, , drop = FALSE]

    message(
      "[", i, "/", nrow(registry), "] ",
      reg$family[[1L]], " / ",
      reg$modelamiento[[1L]], " / ",
      reg$analysis[[1L]]
    )

    result <- tryCatch(
      analyse_candidate_pca(
        reg = reg,
        comparison_root = comparison_root,
        output_dir = output_dir,
        ridge_rel = args$ridge_rel
      ),
      error = function(e) {
        data.frame(
          family = reg$family[[1L]],
          modelamiento = reg$modelamiento[[1L]],
          analysis = reg$analysis[[1L]],
          q00 = reg$q00[[1L]],
          technical_valid = FALSE,
          technical_error = conditionMessage(e),
          n_train = NA_integer_,
          n_test = NA_integer_,
          n_selected_features = NA_integer_,
          mapping_fraction = NA_real_,
          ridge_rel = args$ridge_rel,
          output_dir = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    )

    manifest_rows[[i]] <- result
    invisible(gc(verbose = FALSE))
  }

  manifest <- do.call(rbind, manifest_rows)
  write_tsv(manifest, file.path(output_dir, "00_manifest_pca_escenarios.tsv"))

  writeLines(
    c(
      "PCA OBSERVADO / RECONSTRUIDO / LATENTE",
      paste0("Fecha: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
      paste0("Selection script: ", selection_script),
      paste0("Comparison root: ", comparison_root),
      paste0("Ridge relative: ", args$ridge_rel),
      paste0("Candidatos: ", nrow(manifest)),
      paste0("Validos: ", sum(manifest$technical_valid %in% TRUE)),
      "",
      "REGLA DE VALIDACION:",
      "Todos los PCA se ajustan exclusivamente en TRAIN.",
      "TEST se transforma con medias/desviaciones/rotaciones aprendidas en TRAIN.",
      "En observado/reconstruido, cada vista se balancea por sqrt(numero de features).",
      "Z TEST se recupera solo desde las features seleccionadas del panel.",
      "TEST no participa en seleccion de features, ajuste MOFA ni ajuste PCA."
    ),
    file.path(output_dir, "METODOLOGIA_PCA.txt"),
    useBytes = TRUE
  )

  writeLines(
    capture.output(sessionInfo()),
    file.path(output_dir, "sessionInfo.txt"),
    useBytes = TRUE
  )

  cat(
    "\nCompletado.\n",
    "Salida: ", output_dir, "\n",
    "Validos: ", sum(manifest$technical_valid %in% TRUE),
    " / ", nrow(manifest), "\n",
    sep = ""
  )

  invisible(manifest)
}

if (sys.nframe() == 0L) main_pca()

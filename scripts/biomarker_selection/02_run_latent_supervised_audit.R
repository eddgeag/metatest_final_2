#!/usr/bin/env Rscript
# =============================================================================
# METATEST — Auditoría supervisada bayesiana del espacio latente MOFA
# Version 0.2.0
#
# OBJETIVO CIENTÍFICO
#   Relacionar, a posteriori, los factores latentes aprendidos de forma no
#   supervisada por MOFA con las clases NP / FD / SO.
#
#   Este script NO:
#     - reentrena MOFA;
#     - vuelve a seleccionar escenarios;
#     - selecciona biomarcadores;
#     - ejecuta FGSEA;
#     - usa TEST para elegir hiperparámetros;
#     - convierte la reconstrucción en un nuevo predictor.
#
#   Este script SÍ:
#     1) reutiliza Z_TRAIN y Z_TEST proyectado del framework bidireccional;
#     2) orienta cada factor de forma determinista a partir de W;
#     3) estandariza Z usando exclusivamente TRAIN;
#     4) ajusta por factor:
#
#          Z_f ~ Student-t(mu_g, sigma, nu),  g in {NP, FD, SO}
#
#        mediante:
#
#          Z_f = beta_NP I(NP) + beta_FD I(FD) + beta_SO I(SO) + error_t
#
#     5) estima los contrastes posteriores:
#
#          FD - NP
#          SO - NP
#          FD - SO
#
#     6) usa TEST exclusivamente como auditoría direccional/transportabilidad;
#     7) exporta una tabla que puede enlazarse después con:
#          - GSEA nativo de MOFA: factor -> pathways
#          - framework bidireccional: feature -> clase
#
# NOTA DE IDENTIFICACIÓN
#   El signo de un factor MOFA es arbitrario: (Z_f, W_f) y (-Z_f, -W_f)
#   representan la misma reconstrucción. Para que la dirección sea reproducible,
#   cada factor se orienta haciendo positivo el loading de mayor valor absoluto
#   entre todas las vistas. Esta orientación NO usa las etiquetas NP/FD/SO.
# =============================================================================

options(stringsAsFactors = FALSE, digits = 17, warn = 1)

# ============================= CONFIGURACIÓN ================================
UPSTREAM_DIR <- Sys.getenv(
  "UPSTREAM_DIR",
  unset = paste0(
    "~/Documents/PhD/metatest/results_upstream_local_20260908/",
    "comparacion_modelamientos/upstream_comun"
  )
)

SELECTION_RUN <- Sys.getenv("SELECTION_RUN", unset = "")
FAMILIES <- c("P1", "P2", "P3")
GROUPS <- c("NP", "FD", "SO")

SEED <- as.integer(Sys.getenv("LATENT_BAYES_SEED", unset = "20260917"))

CHAINS <- as.integer(Sys.getenv("LATENT_BAYES_CHAINS", unset = "4"))
ITER <- as.integer(Sys.getenv("LATENT_BAYES_ITER", unset = "4000"))
WARMUP <- as.integer(Sys.getenv("LATENT_BAYES_WARMUP", unset = "2000"))
CORES <- as.integer(Sys.getenv("LATENT_BAYES_CORES", unset = "2"))
ADAPT_DELTA <- as.numeric(Sys.getenv("LATENT_BAYES_ADAPT_DELTA", unset = "0.99"))
MAX_TREEDEPTH <- as.integer(Sys.getenv("LATENT_BAYES_MAX_TREEDEPTH", unset = "12"))

MAKE_PLOTS <- TRUE
SAVE_FITS <- TRUE
NUM_TOL <- 1e-12
# ============================================================================

assert <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x) || !ncol(x)) x <- data.frame(status = character())
  utils::write.table(
    x, path,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    na = "NA",
    qmethod = "double",
    fileEncoding = "UTF-8"
  )
}

canonical_path <- function(p, must_work = TRUE) {
  normalizePath(path.expand(p), winslash = "/", mustWork = must_work)
}

get_slot_safe <- function(object, name) {
  assert(isS4(object), "El objeto MOFA no es S4")
  assert(name %in% names(attributes(object)), paste("Falta slot", name))
  attr(object, name, exact = TRUE)
}

find_latest_selection_run <- function(upstream) {
  base <- file.path(upstream, "downstream", "seleccion_mejor_escenario_train")
  assert(dir.exists(base), paste("No existe:", base))
  dirs <- list.dirs(base, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "SELECTED_SCENARIOS_BY_ANALYSIS.csv"))]
  assert(length(dirs) > 0L, paste("No hay corridas de selección completas en", base))
  info <- file.info(dirs)
  dirs[order(info$mtime, basename(dirs), decreasing = TRUE)][1L]
}

find_bidirectional_run <- function(selection_run, family) {
  root <- file.path(selection_run, family, "BIDIRECTIONAL_FRAMEWORK")
  assert(dir.exists(root), paste("No existe BIDIRECTIONAL_FRAMEWORK para", family))

  latest_file <- file.path(root, "LATEST_RUN.txt")
  if (file.exists(latest_file)) {
    p <- trimws(readLines(latest_file, warn = FALSE)[1L])
    if (nzchar(p) && dir.exists(p) &&
        file.exists(file.path(p, "rds", "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"))) {
      return(canonical_path(p, TRUE))
    }
  }

  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "rds", "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"))]
  assert(length(dirs) > 0L, paste("No hay corrida bidireccional completa para", family))
  info <- file.info(dirs)
  canonical_path(dirs[order(info$mtime, basename(dirs), decreasing = TRUE)][1L], TRUE)
}

extract_W <- function(model, factor_names) {
  ex <- get_slot_safe(model, "expectations")
  assert(is.list(ex) && "W" %in% names(ex), "expectations no contiene W")

  w <- ex$W
  while (is.list(w) && length(w) == 1L && is.list(w[[1L]]) && !is.matrix(w[[1L]])) {
    w <- w[[1L]]
  }
  assert(is.list(w), "W no es una lista por vista")

  out <- vector("list", length(w))
  names(out) <- names(w)

  for (v in names(w)) {
    m <- as.matrix(w[[v]])
    assert(!is.null(rownames(m)), paste("W sin nombres de features en", v))

    if (is.null(colnames(m))) {
      assert(ncol(m) == length(factor_names),
             paste("No se pueden identificar factores en W para", v))
      colnames(m) <- factor_names
    }

    assert(all(factor_names %in% colnames(m)),
           paste("W no contiene todos los factores de Z en", v))

    out[[v]] <- m[, factor_names, drop = FALSE]
  }

  out
}

build_orientation <- function(W, factor_names) {
  rows <- vector("list", length(factor_names))
  names(rows) <- factor_names

  for (f in factor_names) {
    candidates <- do.call(
      rbind,
      lapply(names(W), function(v) {
        m <- W[[v]]
        data.frame(
          view = v,
          feature = rownames(m),
          loading_native = as.numeric(m[, f]),
          stringsAsFactors = FALSE
        )
      })
    )

    candidates$abs_loading <- abs(candidates$loading_native)
    idx <- which.max(candidates$abs_loading)

    anchor <- candidates[idx, , drop = FALSE]
    sgn <- ifelse(anchor$loading_native < 0, -1, 1)

    rows[[f]] <- data.frame(
      factor = f,
      orientation_sign = sgn,
      anchor_view = anchor$view,
      anchor_feature = anchor$feature,
      anchor_loading_native = anchor$loading_native,
      anchor_loading_oriented = anchor$loading_native * sgn,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

orient_scores <- function(Z, orientation) {
  Z <- as.matrix(Z)
  signs <- orientation$orientation_sign
  names(signs) <- orientation$factor
  assert(all(colnames(Z) %in% names(signs)), "Faltan signos de orientación")
  sweep(Z, 2, signs[colnames(Z)], "*")
}

fit_train_scaler <- function(Z) {
  center <- colMeans(Z)
  scale <- apply(Z, 2, stats::sd)

  bad <- !is.finite(scale) | scale <= NUM_TOL
  assert(!any(bad),
         paste("Hay factores con varianza nula/no finita:", paste(names(scale)[bad], collapse = ", ")))

  list(center = center, scale = scale)
}

apply_scaler <- function(Z, scaler) {
  Z <- as.matrix(Z)
  X <- sweep(Z, 2, scaler$center[colnames(Z)], "-")
  sweep(X, 2, scaler$scale[colnames(Z)], "/")
}

posterior_summary_vec <- function(x) {
  x <- as.numeric(x)
  c(
    posterior_mean = mean(x),
    posterior_sd = stats::sd(x),
    q025 = unname(stats::quantile(x, 0.025)),
    q50 = unname(stats::quantile(x, 0.50)),
    q975 = unname(stats::quantile(x, 0.975)),
    p_gt_0 = mean(x > 0),
    p_lt_0 = mean(x < 0),
    p_direction = max(mean(x > 0), mean(x < 0))
  )
}

safe_group_summary <- function(x, y) {
  do.call(
    rbind,
    lapply(GROUPS, function(g) {
      z <- x[y == g]
      data.frame(
        group = g,
        n = length(z),
        mean = mean(z),
        median = stats::median(z),
        sd = if (length(z) > 1L) stats::sd(z) else NA_real_,
        stringsAsFactors = FALSE
      )
    })
  )
}

fit_factor_bayes <- function(score, group, factor_name, seed) {
  dat <- data.frame(
    score = as.numeric(score),
    grupo = factor(as.character(group), levels = GROUPS)
  )

  priors <- c(
    brms::prior(normal(0, 1), class = "b"),
    brms::prior(exponential(1), class = "sigma"),
    brms::prior(gamma(2, 0.1), class = "nu")
  )

  backend <- if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"

  fit <- brms::brm(
    formula = score ~ 0 + grupo,
    data = dat,
    family = brms::student(),
    prior = priors,
    chains = CHAINS,
    iter = ITER,
    warmup = WARMUP,
    cores = CORES,
    seed = seed,
    backend = backend,
    refresh = 0,
    control = list(
      adapt_delta = ADAPT_DELTA,
      max_treedepth = MAX_TREEDEPTH
    )
  )

  draws <- posterior::as_draws_df(fit)

  coef_names <- paste0("b_grupo", GROUPS)
  assert(all(coef_names %in% names(draws)),
         paste("No se encontraron coeficientes de grupo esperados para", factor_name))

  mu_NP <- draws[["b_grupoNP"]]
  mu_FD <- draws[["b_grupoFD"]]
  mu_SO <- draws[["b_grupoSO"]]

  group_draws <- list(NP = mu_NP, FD = mu_FD, SO = mu_SO)
  contrasts <- list(
    FD_minus_NP = mu_FD - mu_NP,
    SO_minus_NP = mu_SO - mu_NP,
    FD_minus_SO = mu_FD - mu_SO
  )

  group_rows <- do.call(
    rbind,
    lapply(names(group_draws), function(g) {
      s <- posterior_summary_vec(group_draws[[g]])
      data.frame(
        factor = factor_name,
        group = g,
        t(s),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })
  )

  contrast_rows <- do.call(
    rbind,
    lapply(names(contrasts), function(k) {
      s <- posterior_summary_vec(contrasts[[k]])
      data.frame(
        factor = factor_name,
        contrast = k,
        t(s),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })
  )

  diag <- posterior::summarise_draws(posterior::as_draws(fit))
  rhat_max <- suppressWarnings(max(diag$rhat, na.rm = TRUE))
  ess_bulk_min <- suppressWarnings(min(diag$ess_bulk, na.rm = TRUE))
  ess_tail_min <- suppressWarnings(min(diag$ess_tail, na.rm = TRUE))

  np <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
  divergences <- if (is.null(np)) {
    NA_integer_
  } else {
    sum(np$Parameter == "divergent__" & np$Value > 0)
  }

  diagnostics <- data.frame(
    factor = factor_name,
    backend = backend,
    chains = CHAINS,
    iter = ITER,
    warmup = WARMUP,
    rhat_max = rhat_max,
    ess_bulk_min = ess_bulk_min,
    ess_tail_min = ess_tail_min,
    divergences = divergences,
    stringsAsFactors = FALSE
  )

  list(
    fit = fit,
    group_posterior = group_rows,
    contrasts = contrast_rows,
    diagnostics = diagnostics
  )
}

make_test_direction_audit <- function(
  Z_test_scaled,
  y_test,
  train_contrasts
) {
  if (is.null(Z_test_scaled) || !nrow(Z_test_scaled)) return(data.frame())

  y_test <- factor(y_test, levels = GROUPS)

  rows <- list()
  rr <- 1L

  for (f in colnames(Z_test_scaled)) {
    sm <- safe_group_summary(Z_test_scaled[, f], y_test)
    means <- setNames(sm$mean, sm$group)

    # unname() evita que R concatene el nombre interno del escalar
    # (p.ej. "FD_minus_NP.FD"), lo que impedía emparejarlo con
    # los nombres exactos de los contrastes de TRAIN.
    test_delta <- c(
      FD_minus_NP = unname(means["FD"] - means["NP"]),
      SO_minus_NP = unname(means["SO"] - means["NP"]),
      FD_minus_SO = unname(means["FD"] - means["SO"])
    )

    tc <- train_contrasts[train_contrasts$factor == f, , drop = FALSE]

    for (contrast in names(test_delta)) {
      tr <- tc[tc$contrast == contrast, , drop = FALSE]
      if (!nrow(tr)) next

      train_med <- tr$q50[1]
      test_val <- unname(test_delta[contrast])

      rows[[rr]] <- data.frame(
        factor = f,
        contrast = contrast,
        train_posterior_median = train_med,
        train_p_direction = tr$p_direction[1],
        test_mean_contrast = test_val,
        direction_concordant =
          is.finite(train_med) &&
          is.finite(test_val) &&
          sign(train_med) != 0 &&
          sign(train_med) == sign(test_val),
        stringsAsFactors = FALSE
      )
      rr <- rr + 1L
    }
  }

  if (!length(rows)) data.frame() else do.call(rbind, rows)
}


latent_lda_test_audit <- function(
  Z_train_scaled,
  y_train,
  Z_test_scaled,
  y_test
) {
  Z_train_scaled <- as.matrix(Z_train_scaled)
  Z_test_scaled <- as.matrix(Z_test_scaled)
  y_train <- factor(y_train, levels = GROUPS)
  y_test <- factor(y_test, levels = GROUPS)

  assert(ncol(Z_train_scaled) >= 1L, "No hay factores para LDA")
  assert(identical(colnames(Z_train_scaled), colnames(Z_test_scaled)),
         "Factores TRAIN/TEST distintos en auditoría LDA")
  assert(all(table(y_train)[GROUPS] > 0), "Alguna clase falta en TRAIN")
  assert(all(table(y_test)[GROUPS] > 0), "Alguna clase falta en TEST")

  # Auditoría mínima y sin tuning:
  # - solo factores latentes;
  # - prior igual por clase;
  # - covarianza común estimada exclusivamente en TRAIN.
  fit <- MASS::lda(
    x = Z_train_scaled,
    grouping = y_train,
    prior = rep(1 / length(GROUPS), length(GROUPS))
  )

  pr <- predict(fit, newdata = Z_test_scaled)
  pred <- factor(as.character(pr$class), levels = GROUPS)

  P <- as.matrix(pr$posterior)
  assert(all(GROUPS %in% colnames(P)), "LDA no devolvió probabilidades para todas las clases")
  P <- P[, GROUPS, drop = FALSE]

  cm <- table(
    reference = y_test,
    predicted = pred
  )

  safe_ratio <- function(num, den) {
    if (!is.finite(den) || den <= 0) return(NA_real_)
    num / den
  }

  class_metrics <- do.call(
    rbind,
    lapply(GROUPS, function(g) {
      TP <- as.numeric(cm[g, g])
      FN <- as.numeric(sum(cm[g, ]) - TP)
      FP <- as.numeric(sum(cm[, g]) - TP)
      TN <- as.numeric(sum(cm) - TP - FN - FP)

      recall <- safe_ratio(TP, TP + FN)
      precision <- safe_ratio(TP, TP + FP)
      specificity <- safe_ratio(TN, TN + FP)

      # Si la clase tiene soporte real y no hay TP, F1 debe ser 0,
      # no NA excluido del macro-promedio.
      f1 <- if ((TP + FN) > 0) {
        den <- 2 * TP + FP + FN
        if (den > 0) 2 * TP / den else 0
      } else {
        NA_real_
      }

      data.frame(
        class = g,
        support = TP + FN,
        TP = TP,
        FN = FN,
        FP = FP,
        TN = TN,
        recall = recall,
        sensitivity = recall,
        precision = precision,
        specificity = specificity,
        F1 = f1,
        stringsAsFactors = FALSE
      )
    })
  )

  idx_true <- cbind(
    seq_len(nrow(P)),
    match(as.character(y_test), GROUPS)
  )
  p_true <- P[idx_true]

  Y <- matrix(
    0,
    nrow = length(y_test),
    ncol = length(GROUPS),
    dimnames = list(NULL, GROUPS)
  )
  Y[idx_true] <- 1

  accuracy <- mean(pred == y_test)
  balanced_accuracy <- mean(class_metrics$recall)
  macro_F1 <- mean(class_metrics$F1)

  summary <- data.frame(
    method = "LDA_equal_priors_on_oriented_scaled_latent_factors",
    role = "exploratory_TEST_audit_not_model_selection",
    n_factors = ncol(Z_train_scaled),
    n_train = nrow(Z_train_scaled),
    n_test = nrow(Z_test_scaled),
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    macro_F1 = macro_F1,
    logloss = -mean(log(pmax(p_true, 1e-15))),
    brier = mean(rowSums((P - Y)^2)),
    stringsAsFactors = FALSE
  )

  predictions <- data.frame(
    sample = rownames(Z_test_scaled),
    true_group = as.character(y_test),
    predicted_group = as.character(pred),
    p_NP = P[, "NP"],
    p_FD = P[, "FD"],
    p_SO = P[, "SO"],
    p_true = p_true,
    stringsAsFactors = FALSE
  )

  confusion <- as.data.frame(cm, stringsAsFactors = FALSE)

  cm_norm <- matrix(
    NA_real_,
    nrow = nrow(cm),
    ncol = ncol(cm),
    dimnames = dimnames(cm)
  )

  for (i in seq_len(nrow(cm))) {
    den <- sum(cm[i, ])
    if (den > 0) cm_norm[i, ] <- cm[i, ] / den
  }

  confusion_normalized <- as.data.frame(
    as.table(cm_norm),
    stringsAsFactors = FALSE
  )
  names(confusion_normalized) <- c("reference", "predicted", "proportion")

  list(
    fit = fit,
    summary = summary,
    predictions = predictions,
    confusion = confusion,
    confusion_normalized = confusion_normalized,
    class_metrics = class_metrics
  )
}

make_score_table <- function(
  Z_train_native,
  Z_test_native,
  Z_train_oriented_scaled,
  Z_test_oriented_scaled,
  y_train,
  y_test
) {
  tr <- data.frame(
    sample = rownames(Z_train_native),
    split = "TRAIN",
    group = as.character(y_train[rownames(Z_train_native)]),
    stringsAsFactors = FALSE
  )

  for (f in colnames(Z_train_native)) {
    tr[[paste0(f, "__native")]] <- Z_train_native[, f]
    tr[[paste0(f, "__oriented_scaled")]] <- Z_train_oriented_scaled[, f]
  }

  if (is.null(Z_test_native) || !nrow(Z_test_native)) return(tr)

  te <- data.frame(
    sample = rownames(Z_test_native),
    split = "TEST_PROJECTED",
    group = as.character(y_test[rownames(Z_test_native)]),
    stringsAsFactors = FALSE
  )

  for (f in colnames(Z_test_native)) {
    te[[paste0(f, "__native")]] <- Z_test_native[, f]
    te[[paste0(f, "__oriented_scaled")]] <- Z_test_oriented_scaled[, f]
  }

  common <- union(names(tr), names(te))
  for (cc in setdiff(common, names(tr))) tr[[cc]] <- NA
  for (cc in setdiff(common, names(te))) te[[cc]] <- NA

  rbind(tr[, common, drop = FALSE], te[, common, drop = FALSE])
}

make_plots <- function(
  scores,
  contrasts,
  outdir,
  family
) {
  if (!MAKE_PLOTS || !requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(NULL))
  }

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Forest plot de contrastes por factor.
  if (nrow(contrasts)) {
    p <- ggplot2::ggplot(
      contrasts,
      ggplot2::aes(
        x = q50,
        y = factor,
        xmin = q025,
        xmax = q975
      )
    ) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::geom_errorbarh(height = 0.15) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~ contrast, scales = "free_y") +
      ggplot2::labs(
        title = paste(family, "| Contrastes bayesianos en el espacio latente"),
        x = "Diferencia posterior en unidades SD de TRAIN",
        y = "Factor"
      ) +
      ggplot2::theme_bw()

    ggplot2::ggsave(
      filename = file.path(outdir, "01_BAYESIAN_LATENT_CONTRASTS.pdf"),
      plot = p, width = 10, height = 6
    )
    ggplot2::ggsave(
      filename = file.path(outdir, "01_BAYESIAN_LATENT_CONTRASTS.png"),
      plot = p, width = 10, height = 6, dpi = 160
    )
  }

  # Scatter Factor1/Factor2 si existen al menos dos factores.
  oriented_cols <- grep("__oriented_scaled$", names(scores), value = TRUE)

  if (length(oriented_cols) >= 2L) {
    xcol <- oriented_cols[1]
    ycol <- oriented_cols[2]

    p2 <- ggplot2::ggplot(
      scores,
      ggplot2::aes(
        x = .data[[xcol]],
        y = .data[[ycol]],
        color = group,
        shape = split
      )
    ) +
      ggplot2::geom_point(size = 3.2, stroke = 1) +
      ggplot2::scale_shape_manual(
        values = c(TRAIN = 17, TEST_PROJECTED = 16)
      ) +
      ggplot2::labs(
        title = paste(family, "| Espacio latente MOFA"),
        x = sub("__oriented_scaled$", "", xcol),
        y = sub("__oriented_scaled$", "", ycol)
      ) +
      ggplot2::theme_bw()

    ggplot2::ggsave(
      filename = file.path(outdir, "02_LATENT_SPACE_FIRST_TWO_FACTORS.pdf"),
      plot = p2, width = 8, height = 6
    )
    ggplot2::ggsave(
      filename = file.path(outdir, "02_LATENT_SPACE_FIRST_TWO_FACTORS.png"),
      plot = p2, width = 8, height = 6, dpi = 160
    )
  }

  invisible(NULL)
}

run_family <- function(selection_run, family, stamp) {
  message("\n============================================================")
  message("LATENT SUPERVISED AUDIT | ", family)
  message("============================================================")

  bidir_run <- find_bidirectional_run(selection_run, family)
  object_path <- file.path(
    bidir_run,
    "rds",
    "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"
  )

  obj <- readRDS(object_path)

  needed <- c(
    "Z_train",
    "Z_test_projected",
    "model_path",
    "selected_scenario_row"
  )
  assert(all(needed %in% names(obj)),
         paste("BIDIRECTIONAL_FRAMEWORK_OBJECT incompleto en", family))

  Z_train_native <- as.matrix(obj$Z_train)
  Z_test_native <- as.matrix(obj$Z_test_projected)

  assert(!is.null(rownames(Z_train_native)), "Z_train sin IDs")
  assert(!is.null(colnames(Z_train_native)), "Z_train sin factores")
  assert(identical(colnames(Z_train_native), colnames(Z_test_native)),
         "Factores TRAIN/TEST no coinciden")

  model_path <- canonical_path(obj$model_path, TRUE)
  model <- readRDS(model_path)

  W <- extract_W(model, colnames(Z_train_native))
  orientation <- build_orientation(W, colnames(Z_train_native))

  Z_train_oriented <- orient_scores(Z_train_native, orientation)
  Z_test_oriented <- orient_scores(Z_test_native, orientation)

  scaler <- fit_train_scaler(Z_train_oriented)
  Z_train_scaled <- apply_scaler(Z_train_oriented, scaler)
  Z_test_scaled <- apply_scaler(Z_test_oriented, scaler)

  # Etiquetas: se toman del ready_for_modeling usado por el mismo framework.
  rfm_path <- canonical_path(obj$rfm_path, TRUE)
  rfm <- readRDS(rfm_path)

  assert(all(c("grupo_train", "grupo_test") %in% names(rfm)),
         "ready_for_modeling no contiene grupo_train/grupo_test")
  assert(all(c("id", "grupo") %in% names(rfm$grupo_train)),
         "grupo_train inválido")
  assert(all(c("id", "grupo") %in% names(rfm$grupo_test)),
         "grupo_test inválido")

  y_train <- factor(
    as.character(rfm$grupo_train$grupo),
    levels = GROUPS
  )
  names(y_train) <- as.character(rfm$grupo_train$id)

  y_test <- factor(
    as.character(rfm$grupo_test$grupo),
    levels = GROUPS
  )
  names(y_test) <- as.character(rfm$grupo_test$id)

  assert(all(rownames(Z_train_scaled) %in% names(y_train)),
         "Faltan etiquetas TRAIN")
  assert(all(rownames(Z_test_scaled) %in% names(y_test)),
         "Faltan etiquetas TEST")

  y_train <- y_train[rownames(Z_train_scaled)]
  y_test <- y_test[rownames(Z_test_scaled)]

  tab_train <- table(y_train)
  assert(all(tab_train[GROUPS] > 0), "Alguna clase falta en TRAIN")

  out <- file.path(
    bidir_run,
    "LATENT_SUPERVISED_AUDIT",
    stamp
  )
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "rds"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "plots"), recursive = TRUE, showWarnings = FALSE)

  # Auditoría descriptiva de TRAIN antes del ajuste.
  train_descriptive <- do.call(
    rbind,
    lapply(colnames(Z_train_scaled), function(f) {
      x <- safe_group_summary(Z_train_scaled[, f], y_train)
      x$factor <- f
      x[, c("factor", "group", "n", "mean", "median", "sd")]
    })
  )

  fits <- list()
  group_post <- list()
  contrasts <- list()
  diagnostics <- list()

  for (ii in seq_along(colnames(Z_train_scaled))) {
    f <- colnames(Z_train_scaled)[ii]
    message("Ajustando ", family, " | ", f, " ...")

    ans <- fit_factor_bayes(
      score = Z_train_scaled[, f],
      group = y_train,
      factor_name = f,
      seed = SEED + match(family, FAMILIES) * 1000L + ii
    )

    fits[[f]] <- ans$fit
    group_post[[f]] <- ans$group_posterior
    contrasts[[f]] <- ans$contrasts
    diagnostics[[f]] <- ans$diagnostics

    if (SAVE_FITS) {
      saveRDS(
        ans$fit,
        file.path(out, "rds", paste0("BRMS_", f, ".rds")),
        version = 3
      )
    }
  }

  group_post <- do.call(rbind, group_post)
  contrasts <- do.call(rbind, contrasts)
  diagnostics <- do.call(rbind, diagnostics)

  test_audit <- make_test_direction_audit(
    Z_test_scaled = Z_test_scaled,
    y_test = y_test,
    train_contrasts = contrasts
  )

  # Clasificación secundaria del espacio latente:
  # referencia simple, sin tuning y entrenada exclusivamente con TRAIN.
  latent_lda_audit <- latent_lda_test_audit(
    Z_train_scaled = Z_train_scaled,
    y_train = y_train,
    Z_test_scaled = Z_test_scaled,
    y_test = y_test
  )

  scores <- make_score_table(
    Z_train_native = Z_train_native,
    Z_test_native = Z_test_native,
    Z_train_oriented_scaled = Z_train_scaled,
    Z_test_oriented_scaled = Z_test_scaled,
    y_train = y_train,
    y_test = y_test
  )

  # Parámetros de orientación y escalado para trazabilidad.
  orientation$train_center_after_orientation <- scaler$center[orientation$factor]
  orientation$train_scale_after_orientation <- scaler$scale[orientation$factor]

  write_csv(
    orientation,
    file.path(out, "tables", "01_FACTOR_ORIENTATION_AND_SCALING.csv")
  )
  write_csv(
    scores,
    file.path(out, "tables", "02_LATENT_SCORES_TRAIN_AND_PROJECTED_TEST.csv")
  )
  write_csv(
    train_descriptive,
    file.path(out, "tables", "03_TRAIN_LATENT_DESCRIPTIVE_BY_GROUP.csv")
  )
  write_csv(
    group_post,
    file.path(out, "tables", "04_TRAIN_BAYESIAN_GROUP_LOCATIONS.csv")
  )
  write_csv(
    contrasts,
    file.path(out, "tables", "05_TRAIN_BAYESIAN_FACTOR_CONTRASTS.csv")
  )
  write_csv(
    diagnostics,
    file.path(out, "tables", "06_BAYESIAN_DIAGNOSTICS.csv")
  )
  write_csv(
    test_audit,
    file.path(out, "tables", "07_TEST_PROJECTED_DIRECTIONAL_AUDIT.csv")
  )
  write_csv(
    latent_lda_audit$summary,
    file.path(out, "tables", "08_TEST_LATENT_LDA_AUDIT_SUMMARY.csv")
  )
  write_csv(
    latent_lda_audit$predictions,
    file.path(out, "tables", "09_TEST_LATENT_LDA_AUDIT_PREDICTIONS.csv")
  )
  write_csv(
    latent_lda_audit$confusion,
    file.path(out, "tables", "10_TEST_LATENT_LDA_AUDIT_CONFUSION.csv")
  )
  write_csv(
    latent_lda_audit$confusion_normalized,
    file.path(out, "tables", "10B_TEST_LATENT_LDA_AUDIT_CONFUSION_NORMALIZED.csv")
  )
  write_csv(
    latent_lda_audit$class_metrics,
    file.path(out, "tables", "10C_TEST_LATENT_LDA_AUDIT_CLASS_METRICS.csv")
  )

  saveRDS(
    list(
      family = family,
      scenario = obj$selected_scenario_row,
      bidirectional_run = bidir_run,
      model_path = model_path,
      rfm_path = rfm_path,
      orientation = orientation,
      scaler = scaler,
      Z_train_native = Z_train_native,
      Z_test_native = Z_test_native,
      Z_train_oriented_scaled = Z_train_scaled,
      Z_test_oriented_scaled = Z_test_scaled,
      y_train = y_train,
      y_test = y_test,
      train_descriptive = train_descriptive,
      group_posterior = group_post,
      contrasts = contrasts,
      diagnostics = diagnostics,
      test_directional_audit = test_audit,
      latent_LDA_TEST_audit = list(
        summary = latent_lda_audit$summary,
        predictions = latent_lda_audit$predictions,
        confusion = latent_lda_audit$confusion,
        confusion_normalized = latent_lda_audit$confusion_normalized,
        class_metrics = latent_lda_audit$class_metrics
      )
    ),
    file.path(out, "rds", "LATENT_SUPERVISED_AUDIT_OBJECT.rds"),
    version = 3
  )

  make_plots(
    scores = scores,
    contrasts = contrasts,
    outdir = file.path(out, "plots"),
    family = family
  )

  readme <- c(
    paste("METATEST |", family, "| AUDITORÍA SUPERVISADA BAYESIANA DEL ESPACIO LATENTE"),
    "",
    paste("BIDIRECTIONAL_FRAMEWORK:", bidir_run),
    paste("MOFA model:", model_path),
    paste("TRAIN:", paste(names(tab_train), as.integer(tab_train), collapse = ", ")),
    paste("TEST projected:", paste(names(table(y_test)), as.integer(table(y_test)), collapse = ", ")),
    paste("Factores:", ncol(Z_train_scaled)),
    "",
    "ROL DEL ANÁLISIS:",
    "- MOFA permanece congelado y no supervisado.",
    "- La supervisión se aplica post-hoc únicamente a los scores latentes Z.",
    "- El objetivo es caracterizar qué factores se desplazan entre NP, FD y SO.",
    "- Los contrastes están expresados en unidades de SD de TRAIN.",
    "",
    "ORIENTACIÓN:",
    "- El signo de cada factor MOFA es arbitrario.",
    "- Para hacer reproducible la dirección, cada factor se orienta usando el loading",
    "  de mayor valor absoluto entre todas las vistas, forzándolo a ser positivo.",
    "- Esta orientación no usa las etiquetas de clase.",
    "",
    "INFERENCIA:",
    "- Modelo por factor: Student-t con una localización específica para NP, FD y SO.",
    "- Contrastes principales: FD-NP, SO-NP y FD-SO.",
    "- p_gt_0 / p_lt_0 son probabilidades posteriores de signo, no p-values.",
    "- No se aplica un umbral automático para declarar 'significativo'.",
    "",
    "TEST:",
    "- TEST no participa en el ajuste bayesiano.",
    "- Se usa para comprobar si la dirección media proyectada coincide con TRAIN.",
    "- Se añade una LDA de referencia usando únicamente los factores latentes:",
    "  se ajusta en TRAIN, usa priors iguales, no hace tuning y se congela antes de TEST.",
    "- Sus métricas (accuracy, balanced accuracy, macro-F1, log-loss y Brier) son",
    "  una auditoría exploratoria de cuánto contenido de clase conserva el espacio latente.",
    "- Dado que TEST ya fue inspeccionado durante el desarrollo, estas métricas NO son",
    "  validación confirmatoria ni se utilizan para cambiar el modelo.",
    "",
    "PUENTE INTERPRETATIVO:",
    "- factor <-> clase: este script.",
    "- factor <-> pathway: GSEA nativo de MOFA.",
    "- feature <-> clase: framework bidireccional.",
    "- La siguiente etapa debe integrar estas tres capas, no volver a entrenar un clasificador de biomarcadores."
  )
  writeLines(readme, file.path(out, "READ_ME_FIRST.txt"), useBytes = TRUE)

  writeLines(out, file.path(bidir_run, "LATENT_SUPERVISED_AUDIT", "LATEST_RUN.txt"))

  data.frame(
    family = family,
    n_factors = ncol(Z_train_scaled),
    n_train = nrow(Z_train_scaled),
    n_test = nrow(Z_test_scaled),
    max_rhat = max(diagnostics$rhat_max, na.rm = TRUE),
    total_divergences = sum(diagnostics$divergences, na.rm = TRUE),
    output = out,
    stringsAsFactors = FALSE
  )
}

main <- function() {
  assert(requireNamespace("brms", quietly = TRUE),
         "Se requiere el paquete brms")
  assert(requireNamespace("posterior", quietly = TRUE),
         "Se requiere el paquete posterior")
  assert(requireNamespace("MOFA2", quietly = TRUE),
         "Se requiere MOFA2 para restaurar correctamente el objeto MOFA")
  assert(requireNamespace("MASS", quietly = TRUE),
         "Se requiere MASS para la auditoría LDA del espacio latente")

  assert(CHAINS >= 2L, "CHAINS debe ser >= 2")
  assert(ITER > WARMUP, "ITER debe ser mayor que WARMUP")
  assert(CORES >= 1L, "CORES debe ser >= 1")

  upstream <- canonical_path(UPSTREAM_DIR, TRUE)

  selection_run <- if (nzchar(SELECTION_RUN)) {
    canonical_path(SELECTION_RUN, TRUE)
  } else {
    find_latest_selection_run(upstream)
  }

  stamp <- paste0(
    format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"),
    "_pid",
    Sys.getpid()
  )

  statuses <- list()

  for (family in FAMILIES) {
    statuses[[family]] <- tryCatch(
      run_family(selection_run, family, stamp),
      error = function(e) {
        data.frame(
          family = family,
          n_factors = NA_integer_,
          n_train = NA_integer_,
          n_test = NA_integer_,
          max_rhat = NA_real_,
          total_divergences = NA_integer_,
          output = NA_character_,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  }

  status <- do.call(rbind, statuses)

  root_out <- file.path(
    selection_run,
    paste0("LATENT_SUPERVISED_AUDIT_RUN_", stamp, ".csv")
  )
  write_csv(status, root_out)

  print(status, row.names = FALSE)

  if ("error" %in% names(status) &&
      any(!is.na(status$error) & nzchar(status$error))) {
    stop("Alguna familia falló. Revisar la tabla de estado.", call. = FALSE)
  }

  invisible(status)
}

if (sys.nframe() == 0L) main()

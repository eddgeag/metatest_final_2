#!/usr/bin/env Rscript
# =============================================================================
# METATEST — Framework bidireccional multi-omico para escenarios MOFA ya elegidos
# Version 1.2.0-frozen
#
# Principios:
#   1) P1, P2 y P3 se analizan INDEPENDIENTEMENTE.
#   2) El escenario ya viene elegido por select_best_upstream_scenario_downstream.R.
#   3) La seleccion/gradacion de biomarcadores usa EXCLUSIVAMENTE TRAIN.
#   4) TEST se proyecta al MOFA de TRAIN congelado y se usa solo para validacion.
#   5) No se reajustan mu/sigma/prior/patrones/evidence_grade con TEST.
#   6) T1/T2/T3 constituyen resultados inferenciales/candidatos del framework,
#      no una firma predictiva final.
#   7) El ensemble T1/T2 en TEST se conserva exclusivamente como auditoria
#      exploratoria de transportabilidad.
#   8) Se exportan pools TRAIN-only de alta confianza para downstream:
#      A_reconstructed + A_observed + mismo patron.
#   9) La prediccion BRMS y FGSEA se realizan en scripts downstream separados.
# Modelo interpretable por feature:
#     X_j | Y=k ~ Student-t(df_t, mu_jk, sigma_jk)
#     P(Y=k | X_j=x) ∝ pi_k f_jk(x)
#
# Dependencia global de baja N:
#   - separacion robusta sobre rangos con peso igual para NP/FD/SO.
#   - permutacion libre de etiquetas exclusivamente en TRAIN.
#   - p/q global no depende de aproximaciones asintoticas.
#
# Direccion de baja N:
#   - medianas robustas por grupo.
#   - jackknife leave-one-out para estabilidad del signo y del grupo superior.
#   - contrastes pairwise con permutacion exacta de rangos cuando es posible.
#
# Salida por familia:
#   <SELECTION_RUN>/P1/BIDIRECTIONAL_FRAMEWORK/<stamp>/
#   <SELECTION_RUN>/P2/BIDIRECTIONAL_FRAMEWORK/<stamp>/
#   <SELECTION_RUN>/P3/BIDIRECTIONAL_FRAMEWORK/<stamp>/
# =============================================================================

options(stringsAsFactors = FALSE, digits = 17, warn = 1)

# ============================= CONFIGURACION ================================
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
VIEW_MAP <- c(
  transcriptomica = "tx",
  proteomica = "pr",
  metabolomica = "me",
  clinical = "cl"
)

SEED <- 20260911L
N_PERM_GLOBAL <- 9999L
N_GRID_DENSITY <- 512L
DF_T <- 4
SCALE_METHOD <- "robust"       # robust | zscore | none
PRIOR_MODE <- "equal"          # equal | train_prevalence
RIDGE_REL <- 1e-6              # proyeccion TEST: ridge relativo a diag(W'W)
SCALE_SHRINKAGE_K <- 4         # regularizacion de sigma grupo -> sigma global (baja N)
CHUNK_SIZE <- 250L
MAKE_PLOTS <- TRUE

# Umbrales OPERATIVOS y editables. Se guardan en PARAMETERS.rds.
# No son equivalencia estadistica ni umbrales biologicos universales.
Q_A <- 0.10
Q_B <- 0.20
P_C <- 0.05
MI_A <- 0.05
MI_DEPENDENCY <- 0.20
JK_A <- 0.875                  # >=14/16 LOO, aproximadamente
JK_B <- 0.75
PAIR_Q_SUPPORT <- 0.20
DIRECTION_CONCORDANCE_MIN <- 2L  # de 3 contrastes obs vs reconstruido

NUM_TOL <- 1e-12
SIGMA_FLOOR <- 1e-4
# ============================================================================

assert <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x) || !ncol(x)) x <- data.frame(status = character())
  utils::write.table(
    x, path, sep = ",", row.names = FALSE, col.names = TRUE,
    na = "NA", qmethod = "double", fileEncoding = "UTF-8"
  )
}

pkg_version <- function(p) {
  if (requireNamespace(p, quietly = TRUE)) {
    as.character(utils::packageVersion(p))
  } else {
    "not_installed"
  }
}

canonical_path <- function(p, must_work = TRUE) {
  normalizePath(path.expand(p), winslash = "/", mustWork = must_work)
}

get_slot <- function(object, name) {
  assert(isS4(object), "El objeto MOFA no es S4")
  assert(name %in% names(attributes(object)), paste("Falta slot", name))
  attr(object, name, exact = TRUE)
}

fmt_num <- function(x, digits = 4) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
}

safe_scale_value <- function(x, fallback = 1) {
  x <- as.numeric(x)
  x[!is.finite(x) | x <= NUM_TOL] <- fallback
  x
}

argmax_name <- function(x) {
  if (!length(x) || all(!is.finite(x))) return(NA_character_)
  names(x)[which.max(x)]
}

row_argmax <- function(m) {
  m <- as.matrix(m)
  colnames(m)[max.col(m, ties.method = "first")]
}

balanced_accuracy <- function(truth, pred, levels = GROUPS) {
  truth <- factor(truth, levels = levels)
  pred <- factor(pred, levels = levels)
  recalls <- vapply(levels, function(g) {
    idx <- which(truth == g)
    if (!length(idx)) return(NA_real_)
    mean(pred[idx] == g)
  }, numeric(1))
  mean(recalls, na.rm = TRUE)
}

binary_auc_rank <- function(y01, score) {
  y01 <- as.integer(y01)
  ok <- is.finite(score) & !is.na(y01)
  y01 <- y01[ok]; score <- score[ok]
  n1 <- sum(y01 == 1L); n0 <- sum(y01 == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y01 == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

macro_auc_ovr <- function(truth, prob_mat, levels = GROUPS) {
  truth <- factor(truth, levels = levels)
  P <- as.matrix(prob_mat)[, levels, drop = FALSE]
  aucs <- vapply(levels, function(g) {
    binary_auc_rank(as.integer(truth == g), P[, g])
  }, numeric(1))
  if (all(!is.finite(aucs))) NA_real_ else mean(aucs[is.finite(aucs)])
}

# ------------------------------ CLI -----------------------------------------
parse_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key %in% c("--help", "-h")) {
      cat(
        "Uso:\n",
        "  Rscript run_bidirectional_framework_selected_scenarios.R \\\n",
        "    --selection-run /ruta/al/202609..._pid... \\\n",
        "    --families P1,P2,P3\n\n",
        "Opciones:\n",
        "  --upstream DIR\n",
        "  --selection-run DIR\n",
        "  --families P1,P2,P3\n",
        "  --seed N\n",
        "  --n-perm N\n",
        "  --n-grid N\n",
        "  --df-t X\n",
        "  --scale robust|zscore|none\n",
        "  --prior equal|train_prevalence\n",
        "  --ridge X\n",
        "  --scale-shrink-k X\n",
        "  --chunk-size N\n",
        "  --no-plots\n",
        sep = ""
      )
      return(FALSE)
    }
    if (key == "--no-plots") {
      MAKE_PLOTS <<- FALSE; i <- i + 1L; next
    }

    allowed <- c(
      "--upstream", "--selection-run", "--families", "--seed", "--n-perm",
      "--n-grid", "--df-t", "--scale", "--prior", "--ridge", "--scale-shrink-k", "--chunk-size"
    )
    assert(key %in% allowed && i < length(args), paste("Argumento invalido/incompleto:", key))
    val <- args[i + 1L]
    switch(
      key,
      `--upstream` = { UPSTREAM_DIR <<- val },
      `--selection-run` = { SELECTION_RUN <<- val },
      `--families` = { FAMILIES <<- strsplit(val, ",", fixed = TRUE)[[1L]] },
      `--seed` = { SEED <<- as.integer(val) },
      `--n-perm` = { N_PERM_GLOBAL <<- as.integer(val) },
      `--n-grid` = { N_GRID_DENSITY <<- as.integer(val) },
      `--df-t` = { DF_T <<- as.numeric(val) },
      `--scale` = { SCALE_METHOD <<- tolower(val) },
      `--prior` = { PRIOR_MODE <<- tolower(val) },
      `--ridge` = { RIDGE_REL <<- as.numeric(val) },
      `--scale-shrink-k` = { SCALE_SHRINKAGE_K <<- as.numeric(val) },
      `--chunk-size` = { CHUNK_SIZE <<- as.integer(val) }
    )
    i <- i + 2L
  }
  TRUE
}

find_latest_selection_run <- function(upstream) {
  base <- file.path(upstream, "downstream", "seleccion_mejor_escenario_train")
  assert(dir.exists(base), paste("No existe el directorio de seleccion:", base))
  dirs <- list.dirs(base, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "SELECTED_SCENARIOS_BY_ANALYSIS.csv"))]
  assert(length(dirs) > 0L, paste("No hay corridas de seleccion completas en", base))
  info <- file.info(dirs)
  dirs[order(info$mtime, basename(dirs), decreasing = TRUE)][1L]
}

# -------------------------- Escenario seleccionado --------------------------
read_selected_row <- function(selection_run, family) {
  family_dir <- file.path(selection_run, family)
  assert(dir.exists(family_dir), paste("No existe carpeta de familia:", family_dir))

  f <- file.path(family_dir, "SELECTED_SCENARIO.csv")
  if (file.exists(f)) {
    x <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    root_f <- file.path(selection_run, "SELECTED_SCENARIOS_BY_ANALYSIS.csv")
    assert(file.exists(root_f), paste("No existe", f, "ni", root_f))
    all <- utils::read.csv(root_f, stringsAsFactors = FALSE, check.names = FALSE)
    x <- all[all$family == family, , drop = FALSE]
  }

  assert(nrow(x) == 1L, paste("Se esperaba exactamente un escenario seleccionado para", family))
  assert(all(c("family", "scenario", "model_path") %in% names(x)), "Tabla de seleccion incompleta")
  assert(identical(as.character(x$family[1]), family), "Familia de seleccion inconsistente")
  x
}

scenario_dir_from_model_path <- function(model_path) {
  model_path <- canonical_path(model_path, must_work = TRUE)
  # <scenario>/integracion/rds/modelo_mofa_final_*.rds
  canonical_path(dirname(dirname(dirname(model_path))), must_work = TRUE)
}

# ------------------------------ MOFA / RFM ----------------------------------
extract_Z_W <- function(model, train_ids) {
  assert(isS4(model) && "MOFA" %in% class(model), "El modelo final no es un S4 MOFA")
  ex <- get_slot(model, "expectations")
  assert(is.list(ex) && all(c("Z", "W") %in% names(ex)), "expectations no contiene Z/W")

  z_obj <- ex$Z
  while (is.list(z_obj) && length(z_obj) == 1L && !is.matrix(z_obj)) z_obj <- z_obj[[1L]]
  z <- as.matrix(z_obj)
  if (!all(train_ids %in% rownames(z)) && all(train_ids %in% colnames(z))) z <- t(z)
  assert(all(train_ids %in% rownames(z)), "Z no contiene todos los IDs de TRAIN")
  z <- z[train_ids, , drop = FALSE]
  assert(ncol(z) >= 1L && all(is.finite(z)), "Z TRAIN invalido")
  if (is.null(colnames(z))) colnames(z) <- paste0("Factor", seq_len(ncol(z)))

  w <- ex$W
  while (is.list(w) && length(w) == 1L && is.list(w[[1L]]) && !is.matrix(w[[1L]])) {
    w <- w[[1L]]
  }
  assert(is.list(w), "W no es lista por vista")
  assert(all(names(VIEW_MAP) %in% names(w)), "W no contiene las cuatro vistas esperadas")

  for (v in names(VIEW_MAP)) {
    wv <- as.matrix(w[[v]])
    assert(!is.null(rownames(wv)), paste("W sin rownames en", v))
    if (is.null(colnames(wv))) colnames(wv) <- colnames(z)
    assert(setequal(colnames(wv), colnames(z)), paste("Factores W/Z no coinciden en", v))
    w[[v]] <- wv[, colnames(z), drop = FALSE]
  }
  list(Z = z, W = w)
}

validate_rfm <- function(rfm, selection_run) {
  needed <- c(unname(VIEW_MAP), "grupo_train", "grupo_test")
  assert(all(needed %in% names(rfm)), "ready_for_modeling no contiene vistas/grupos esperados")
  assert(all(c("id", "grupo") %in% names(rfm$grupo_train)), "grupo_train sin id/grupo")
  assert(all(c("id", "grupo") %in% names(rfm$grupo_test)), "grupo_test sin id/grupo")

  tr <- rfm$grupo_train[, c("id", "grupo"), drop = FALSE]
  te <- rfm$grupo_test[, c("id", "grupo"), drop = FALSE]
  tr[] <- lapply(tr, as.character); te[] <- lapply(te, as.character)
  assert(!anyDuplicated(tr$id) && !anyDuplicated(te$id), "IDs duplicados en split")
  assert(!length(intersect(tr$id, te$id)), "TRAIN y TEST comparten IDs")
  assert(all(tr$grupo %in% GROUPS) && all(te$grupo %in% GROUPS), "Grupos inesperados")

  mf <- file.path(selection_run, "TRAIN_MANIFEST_USED.csv")
  if (file.exists(mf)) {
    manifest <- utils::read.csv(mf, stringsAsFactors = FALSE, check.names = FALSE)
    assert(all(c("sample", "grupo") %in% names(manifest)), "TRAIN_MANIFEST_USED.csv invalido")
    manifest$sample <- as.character(manifest$sample)
    manifest$grupo <- as.character(manifest$grupo)
    assert(setequal(tr$id, manifest$sample), "TRAIN de ready_for_modeling != TRAIN usado en seleccion")
    idx <- match(manifest$sample, tr$id)
    assert(identical(tr$grupo[idx], manifest$grupo), "Etiquetas TRAIN != manifiesto de seleccion")
    tr <- tr[idx, , drop = FALSE]
  }

  for (v in names(VIEW_MAP)) {
    code <- VIEW_MAP[[v]]
    obj <- rfm[[code]]
    assert(is.list(obj) && all(c("train", "test") %in% names(obj)), paste("Vista", code, "sin train/test"))
    xtr <- as.matrix(obj$train); xte <- as.matrix(obj$test)
    assert(all(tr$id %in% colnames(xtr)), paste("TRAIN IDs faltan en", code))
    assert(all(te$id %in% colnames(xte)), paste("TEST IDs faltan en", code))
    assert(!is.null(rownames(xtr)) && !is.null(rownames(xte)), paste("Features sin nombre en", code))
    assert(setequal(rownames(xtr), rownames(xte)), paste("Features TRAIN/TEST distintas en", code))
    assert(all(is.finite(xtr)) && all(is.finite(xte)), paste("NA/Inf en", code))
  }
  list(train = tr, test = te)
}

ridge_project_view <- function(W, X, ridge_rel = RIDGE_REL) {
  # W: features x factors; X: features x samples
  G <- crossprod(W)
  diag_mean <- mean(diag(G))
  lambda <- ridge_rel * max(diag_mean, NUM_TOL)
  A <- G + diag(lambda, nrow(G))
  rhs <- crossprod(W, X)
  sol <- tryCatch(solve(A, rhs), error = function(e) qr.solve(A, rhs, tol = 1e-10))
  z <- t(sol)
  colnames(z) <- colnames(W)
  rownames(z) <- colnames(X)
  attr(z, "ridge_lambda") <- lambda
  z
}

assemble_family_data <- function(selection_run, family, selected_row) {
  model_path <- canonical_path(selected_row$model_path[1], TRUE)
  scenario_dir <- scenario_dir_from_model_path(model_path)
  scenario <- basename(scenario_dir)
  assert(identical(as.character(selected_row$scenario[1]), scenario), "scenario != basename(model_path)")

  rfm_path <- file.path(scenario_dir, paste0("ready_for_modeling_", scenario, ".rds"))
  assert(file.exists(rfm_path), paste("No existe", rfm_path))

  rfm <- readRDS(rfm_path)
  split <- validate_rfm(rfm, selection_run)
  train_ids <- split$train$id
  test_ids <- split$test$id
  y_train <- factor(split$train$grupo, levels = GROUPS)
  y_test <- factor(split$test$grupo, levels = GROUPS)
  names(y_train) <- train_ids; names(y_test) <- test_ids

  model <- readRDS(model_path)
  zw <- extract_Z_W(model, train_ids)
  Z_train <- zw$Z
  W <- zw$W

  feature_map_parts <- list()
  obs_train_parts <- list(); obs_test_parts <- list()
  rec_train_parts <- list(); rec_test_parts <- list()
  projection_parts <- list(); projection_audit <- list()

  for (v in names(VIEW_MAP)) {
    code <- VIEW_MAP[[v]]
    Xtr0 <- as.matrix(rfm[[code]]$train)
    Xte0 <- as.matrix(rfm[[code]]$test)
    Xtr0 <- Xtr0[, train_ids, drop = FALSE]
    Xte0 <- Xte0[, test_ids, drop = FALSE]

    Wv0 <- W[[v]]
    common <- intersect(rownames(Wv0), rownames(Xtr0))
    common <- common[common %in% rownames(Xte0)]
    assert(length(common) >= ncol(Z_train), paste("Muy pocas features comunes para proyectar", v))

    Wv <- Wv0[common, colnames(Z_train), drop = FALSE]
    Xtr <- Xtr0[common, , drop = FALSE]
    Xte <- Xte0[common, , drop = FALSE]

    Zte_v <- ridge_project_view(Wv, Xte, RIDGE_REL)
    Ztr_proj_v <- ridge_project_view(Wv, Xtr, RIDGE_REL)
    projection_parts[[v]] <- Zte_v

    projection_audit[[v]] <- data.frame(
      family = family,
      scenario = scenario,
      view = v,
      view_code = code,
      n_features_model = nrow(Wv0),
      n_features_common = length(common),
      n_factors = ncol(Wv),
      ridge_lambda = attr(Zte_v, "ridge_lambda"),
      train_projection_rmse_vs_mofa_Z = sqrt(mean((Ztr_proj_v - Z_train)^2)),
      stringsAsFactors = FALSE
    )

    feature_model <- paste0(code, "__", common)
    fmap <- data.frame(
      feature_model = feature_model,
      feature_original = common,
      view = v,
      view_code = code,
      stringsAsFactors = FALSE
    )
    feature_map_parts[[v]] <- fmap

    # Observado: samples x features.
    Otr <- t(Xtr); Ote <- t(Xte)
    colnames(Otr) <- feature_model; colnames(Ote) <- feature_model

    # Reconstruido TRAIN usando Z posterior de MOFA.
    Rtr <- Z_train %*% t(Wv)
    colnames(Rtr) <- feature_model

    obs_train_parts[[v]] <- Otr
    obs_test_parts[[v]] <- Ote
    rec_train_parts[[v]] <- Rtr
  }

  # TEST: cada vista proyecta Z con el mismo W de TRAIN; las vistas pesan igual.
  Z_test_array <- simplify2array(projection_parts)
  if (length(dim(Z_test_array)) == 2L) {
    Z_test <- Z_test_array
  } else {
    Z_test <- apply(Z_test_array, c(1, 2), mean)
  }
  rownames(Z_test) <- test_ids
  colnames(Z_test) <- colnames(Z_train)

  # Reconstruccion TEST con Z_test comun.
  for (v in names(VIEW_MAP)) {
    code <- VIEW_MAP[[v]]
    fmap <- feature_map_parts[[v]]
    common <- fmap$feature_original
    Wv <- W[[v]][common, colnames(Z_train), drop = FALSE]
    Rte <- Z_test %*% t(Wv)
    colnames(Rte) <- fmap$feature_model
    rec_test_parts[[v]] <- Rte
  }

  feature_map <- do.call(rbind, feature_map_parts)
  rownames(feature_map) <- NULL
  assert(!anyDuplicated(feature_map$feature_model), "feature_model duplicado tras prefijo de vista")

  obs_train <- do.call(cbind, obs_train_parts)
  obs_test <- do.call(cbind, obs_test_parts)
  rec_train <- do.call(cbind, rec_train_parts)
  rec_test <- do.call(cbind, rec_test_parts)
  common_cols <- feature_map$feature_model
  obs_train <- obs_train[train_ids, common_cols, drop = FALSE]
  obs_test <- obs_test[test_ids, common_cols, drop = FALSE]
  rec_train <- rec_train[train_ids, common_cols, drop = FALSE]
  rec_test <- rec_test[test_ids, common_cols, drop = FALSE]

  list(
    family = family,
    scenario = scenario,
    scenario_dir = scenario_dir,
    model_path = model_path,
    rfm_path = rfm_path,
    rfm = rfm,
    feature_map = feature_map,
    y_train = y_train,
    y_test = y_test,
    Z_train = Z_train,
    Z_test = Z_test,
    observed_train = obs_train,
    observed_test = obs_test,
    reconstructed_train = rec_train,
    reconstructed_test = rec_test,
    projection_audit = do.call(rbind, projection_audit)
  )
}

# --------------------------- Escalado TRAIN-only ----------------------------
fit_scaler <- function(X, method = SCALE_METHOD) {
  X <- as.matrix(X)
  if (method == "none") {
    center <- rep(0, ncol(X)); scale <- rep(1, ncol(X))
  } else if (method == "zscore") {
    center <- colMeans(X)
    scale <- apply(X, 2, stats::sd)
  } else if (method == "robust") {
    center <- apply(X, 2, stats::median)
    scale <- apply(X, 2, stats::mad)
    bad <- !is.finite(scale) | scale <= NUM_TOL
    if (any(bad)) scale[bad] <- apply(X[, bad, drop = FALSE], 2, stats::sd)
  } else {
    stop("SCALE_METHOD invalido", call. = FALSE)
  }
  scale <- safe_scale_value(scale, 1)
  names(center) <- colnames(X); names(scale) <- colnames(X)
  list(center = center, scale = scale, method = method)
}

apply_scaler <- function(X, scaler) {
  X <- as.matrix(X)
  X <- sweep(X, 2, scaler$center[colnames(X)], "-")
  sweep(X, 2, scaler$scale[colnames(X)], "/")
}

# -------------------------- Banco de permutaciones --------------------------
make_global_bank <- function(y, n_perm, seed) {
  set.seed(seed)
  y <- factor(y, levels = GROUPS)
  n <- length(y)
  idx <- t(replicate(n_perm, sample.int(n)))
  memberships <- lapply(GROUPS, function(g) matrix(as.numeric(y[idx] == g), nrow = n_perm))
  names(memberships) <- GROUPS
  list(indices = idx, memberships = memberships, n_perm = n_perm, seed = seed)
}

rank_matrix <- function(X) {
  X <- as.matrix(X)
  R <- apply(X, 2, rank, ties.method = "average", na.last = "keep")
  if (is.null(dim(R))) R <- matrix(R, ncol = 1L)
  colnames(R) <- colnames(X); rownames(R) <- rownames(X)
  R
}

global_rank_equal_prior_permutation <- function(X, y, bank, chunk_size = CHUNK_SIZE) {
  # Estadistica global robusta y balanceada entre los tres grupos.
  # Se usa una separacion sobre rangos con peso IGUAL para NP/FD/SO:
  #   B = Var_equal(mu_NP, mu_FD, mu_SO)
  #   W = mean_equal(var_NP, var_FD, var_SO)
  #   E = B / (B + W)
  # donde las varianzas intragrupo usan denominador n_g - 1.
  # Esto evita que FD=3 quede infraponderado por prevalencia muestral.
  X <- as.matrix(X)
  y <- factor(y, levels = GROUPS)
  R <- rank_matrix(X)
  K <- length(GROUPS)

  group_means <- sapply(GROUPS, function(g) colMeans(R[y == g, , drop = FALSE]))
  group_vars <- sapply(GROUPS, function(g) {
    apply(R[y == g, , drop = FALSE], 2, stats::var)
  })
  if (is.null(dim(group_means))) group_means <- matrix(group_means, ncol = K)
  if (is.null(dim(group_vars))) group_vars <- matrix(group_vars, ncol = K)
  colnames(group_means) <- GROUPS
  colnames(group_vars) <- GROUPS

  mu_equal <- rowMeans(group_means)
  between <- rowMeans((group_means - mu_equal)^2)
  within <- rowMeans(group_vars)
  observed <- ifelse(
    between + within > NUM_TOL,
    pmin(1, pmax(0, between / (between + within))),
    0
  )
  p <- rep(1, ncol(R))

  starts <- seq.int(1L, ncol(R), by = chunk_size)
  for (s in starts) {
    idx <- s:min(ncol(R), s + chunk_size - 1L)
    r <- R[, idx, drop = FALSE]
    r2 <- r^2
    B <- bank$n_perm
    sum_mu <- matrix(0, nrow = B, ncol = length(idx))
    sum_mu2 <- matrix(0, nrow = B, ncol = length(idx))
    sum_var <- matrix(0, nrow = B, ncol = length(idx))

    for (g in GROUPS) {
      ng <- sum(y == g)
      M <- bank$memberships[[g]]
      sums <- M %*% r
      sums2 <- M %*% r2
      means <- sums / ng
      vars <- pmax(0, sums2 - sums^2 / ng) / (ng - 1)
      sum_mu <- sum_mu + means
      sum_mu2 <- sum_mu2 + means^2
      sum_var <- sum_var + vars
    }

    mu_bar <- sum_mu / K
    null_between <- pmax(0, sum_mu2 / K - mu_bar^2)
    null_within <- sum_var / K
    denom <- null_between + null_within
    null_stat <- ifelse(denom > NUM_TOL, null_between / denom, 0)
    obs <- observed[idx]
    p[idx] <- vapply(seq_along(idx), function(j) {
      (1 + sum(null_stat[, j] >= obs[j] - NUM_TOL)) / (bank$n_perm + 1)
    }, numeric(1))
  }

  data.frame(
    feature_model = colnames(X),
    rank_equal_group_separation = observed,
    p_global = p,
    q_global = stats::p.adjust(p, method = "BH"),
    global_permutation_stat = "rank_equal_prior_between_over_between_plus_equal_within",
    stringsAsFactors = FALSE
  )
}

pairwise_exact_rank_permutation <- function(X, y, chunk_size = CHUNK_SIZE) {
  X <- as.matrix(X)
  y <- factor(y, levels = GROUPS)
  pairs <- list(c("FD", "NP"), c("SO", "NP"), c("SO", "FD"))
  out <- data.frame(feature_model = colnames(X), stringsAsFactors = FALSE)

  for (pair in pairs) {
    Aname <- pair[1]; Bname <- pair[2]
    rows <- which(y %in% pair)
    yy <- droplevels(factor(y[rows], levels = pair))
    Xp <- X[rows, , drop = FALSE]
    Rp <- rank_matrix(Xp)
    na <- sum(yy == Aname); nb <- sum(yy == Bname); n <- length(rows)
    combos <- utils::combn(n, na)
    A <- matrix(0, nrow = ncol(combos), ncol = n)
    for (i in seq_len(nrow(A))) A[i, combos[, i]] <- 1
    nalloc <- nrow(A)

    obs_rank_delta <- colMeans(Rp[yy == Aname, , drop = FALSE]) -
      colMeans(Rp[yy == Bname, , drop = FALSE])
    pvals <- numeric(ncol(Xp))

    starts <- seq.int(1L, ncol(Xp), by = chunk_size)
    for (s in starts) {
      idx <- s:min(ncol(Xp), s + chunk_size - 1L)
      r <- Rp[, idx, drop = FALSE]
      sumA <- A %*% r
      meanA <- sumA / na
      total <- matrix(colSums(r), nrow = nalloc, ncol = length(idx), byrow = TRUE)
      meanB <- (total - sumA) / nb
      nul <- abs(meanA - meanB)
      obs <- abs(obs_rank_delta[idx])
      pvals[idx] <- vapply(seq_along(idx), function(j) {
        sum(nul[, j] >= obs[j] - NUM_TOL) / nalloc
      }, numeric(1))
    }

    tag <- paste0(Aname, "_minus_", Bname)
    out[[paste0("rank_delta__", tag)]] <- obs_rank_delta
    out[[paste0("p_pair__", tag)]] <- pvals
    out[[paste0("q_pair__", tag)]] <- stats::p.adjust(pvals, method = "BH")
    out[[paste0("n_alloc__", tag)]] <- nalloc
  }
  out
}

# ------------------------- Student-t / MI / JS ------------------------------
make_priors <- function(y, mode = PRIOR_MODE) {
  if (mode == "equal") {
    p <- rep(1 / length(GROUPS), length(GROUPS)); names(p) <- GROUPS
    return(p)
  }
  if (mode == "train_prevalence") {
    tab <- table(factor(y, levels = GROUPS))
    p <- as.numeric(tab) / sum(tab); names(p) <- GROUPS
    return(p)
  }
  stop("PRIOR_MODE invalido", call. = FALSE)
}

fit_group_params <- function(x, y, shrinkage_k = SCALE_SHRINKAGE_K) {
  y <- factor(y, levels = GROUPS)
  mu <- sigma_raw <- sigma <- n <- setNames(numeric(length(GROUPS)), GROUPS)

  global_scale <- stats::mad(x)
  if (!is.finite(global_scale) || global_scale <= NUM_TOL) global_scale <- stats::sd(x)
  if (!is.finite(global_scale) || global_scale <= NUM_TOL) global_scale <- 1
  global_scale <- max(global_scale, SIGMA_FLOOR)

  for (g in GROUPS) {
    z <- x[y == g]
    n[g] <- length(z)
    mu[g] <- stats::median(z)
    s <- stats::mad(z)
    if (!is.finite(s) || s <= NUM_TOL) s <- stats::sd(z)
    if (!is.finite(s) || s <= NUM_TOL) s <- global_scale
    s <- max(s, SIGMA_FLOOR)
    sigma_raw[g] <- s

    # Contraccion de varianza hacia la escala global. Con N pequeno, especialmente
    # FD=3, evita densidades excesivamente estrechas por azar. Al crecer n_g,
    # aumenta el peso de la escala especifica del grupo.
    w <- n[g] / (n[g] + max(shrinkage_k, 0))
    sigma[g] <- sqrt(w * s^2 + (1 - w) * global_scale^2)
    sigma[g] <- max(sigma[g], SIGMA_FLOOR)
  }
  list(mu = mu, sigma_raw = sigma_raw, sigma = sigma, n = n, global_scale = global_scale)
}

t_density <- function(x, mu, sigma, df_t = DF_T) {
  sigma <- max(sigma, SIGMA_FLOOR)
  stats::dt((x - mu) / sigma, df = df_t) / sigma
}

posterior_at_x <- function(x, mu, sigma, priors, df_t = DF_T) {
  dens <- vapply(GROUPS, function(g) t_density(x, mu[g], sigma[g], df_t), numeric(1))
  num <- priors[GROUPS] * dens
  den <- sum(num)
  if (!is.finite(den) || den <= 0) return(setNames(rep(NA_real_, length(GROUPS)), GROUPS))
  out <- num / den
  names(out) <- GROUPS
  out
}

mi_js_feature <- function(mu, sigma, priors, n_grid = N_GRID_DENSITY, df_t = DF_T) {
  xmin <- min(mu - 8 * sigma)
  xmax <- max(mu + 8 * sigma)
  if (!is.finite(xmin) || !is.finite(xmax) || xmax <= xmin + NUM_TOL) {
    return(c(MI_bits = 0, MI_norm = 0, JS_FD_NP = 0, JS_SO_NP = 0, JS_SO_FD = 0))
  }
  grid <- seq(xmin, xmax, length.out = n_grid)
  dx <- mean(diff(grid))
  dens <- sapply(GROUPS, function(g) t_density(grid, mu[g], sigma[g], df_t))
  colnames(dens) <- GROUPS
  mix <- as.numeric(dens %*% priors[GROUPS])
  eps <- 1e-15
  mi_nats <- sum(vapply(GROUPS, function(g) {
    sum(priors[g] * dens[, g] * log((dens[, g] + eps) / (mix + eps))) * dx
  }, numeric(1)))
  hy <- -sum(priors * log(priors + eps))

  js_pair <- function(a, b) {
    fa <- dens[, a]; fb <- dens[, b]; m <- 0.5 * (fa + fb)
    val <- 0.5 * sum(fa * log((fa + eps) / (m + eps))) * dx +
      0.5 * sum(fb * log((fb + eps) / (m + eps))) * dx
    val / log(2)
  }
  c(
    MI_bits = mi_nats / log(2),
    MI_norm = if (hy > 0) mi_nats / hy else 0,
    JS_FD_NP = js_pair("FD", "NP"),
    JS_SO_NP = js_pair("SO", "NP"),
    JS_SO_FD = js_pair("SO", "FD")
  )
}

student_t_descriptors <- function(X, y, priors, n_grid = N_GRID_DENSITY, df_t = DF_T) {
  X <- as.matrix(X)
  rows <- vector("list", ncol(X))
  for (j in seq_len(ncol(X))) {
    x <- X[, j]
    par <- fit_group_params(x, y)
    info <- mi_js_feature(par$mu, par$sigma, priors, n_grid, df_t)
    q10 <- as.numeric(stats::quantile(x, 0.10, names = FALSE))
    q50 <- as.numeric(stats::quantile(x, 0.50, names = FALSE))
    q90 <- as.numeric(stats::quantile(x, 0.90, names = FALSE))
    plo <- posterior_at_x(q10, par$mu, par$sigma, priors, df_t)
    pmid <- posterior_at_x(q50, par$mu, par$sigma, priors, df_t)
    phi <- posterior_at_x(q90, par$mu, par$sigma, priors, df_t)
    top <- argmax_name(par$mu)
    high <- argmax_name(phi)
    enrich <- if (is.na(high)) NA_real_ else log(pmax(phi[high], 1e-15) / pmax(priors[high], 1e-15))

    rows[[j]] <- data.frame(
      feature_model = colnames(X)[j],
      n_used = length(x),
      feature_sd = stats::sd(x),
      constant_feature = !is.finite(stats::sd(x)) || stats::sd(x) <= NUM_TOL,
      MI_bits = unname(info["MI_bits"]),
      MI_norm = unname(info["MI_norm"]),
      JS_FD_NP = unname(info["JS_FD_NP"]),
      JS_SO_NP = unname(info["JS_SO_NP"]),
      JS_SO_FD = unname(info["JS_SO_FD"]),
      mu_NP = par$mu["NP"], mu_FD = par$mu["FD"], mu_SO = par$mu["SO"],
      sigma_raw_NP = par$sigma_raw["NP"], sigma_raw_FD = par$sigma_raw["FD"], sigma_raw_SO = par$sigma_raw["SO"],
      sigma_NP = par$sigma["NP"], sigma_FD = par$sigma["FD"], sigma_SO = par$sigma["SO"],
      sigma_global = par$global_scale,
      prior_NP = priors["NP"], prior_FD = priors["FD"], prior_SO = priors["SO"],
      q10 = q10, q50 = q50, q90 = q90,
      post_low_NP = plo["NP"], post_low_FD = plo["FD"], post_low_SO = plo["SO"],
      post_mid_NP = pmid["NP"], post_mid_FD = pmid["FD"], post_mid_SO = pmid["SO"],
      post_high_NP = phi["NP"], post_high_FD = phi["FD"], post_high_SO = phi["SO"],
      top_location_group = top,
      posterior_high_group = high,
      concordant_high = !is.na(top) && !is.na(high) && identical(top, high),
      posterior_enrichment_high = enrich,
      delta_FD_minus_NP = par$mu["FD"] - par$mu["NP"],
      delta_SO_minus_NP = par$mu["SO"] - par$mu["NP"],
      delta_SO_minus_FD = par$mu["SO"] - par$mu["FD"],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# ------------------------------- Jackknife ----------------------------------
jackknife_direction <- function(X, y, full_desc) {
  X <- as.matrix(X)
  y <- factor(y, levels = GROUPS)
  n <- nrow(X); p <- ncol(X)
  full_top <- full_desc$top_location_group
  full_d <- cbind(
    FD_minus_NP = full_desc$delta_FD_minus_NP,
    SO_minus_NP = full_desc$delta_SO_minus_NP,
    SO_minus_FD = full_desc$delta_SO_minus_FD
  )
  top_same <- matrix(FALSE, nrow = n, ncol = p)
  delta_arr <- array(NA_real_, dim = c(n, p, 3L),
                     dimnames = list(NULL, colnames(X), colnames(full_d)))

  for (i in seq_len(n)) {
    keep <- setdiff(seq_len(n), i)
    med <- sapply(GROUPS, function(g) {
      apply(X[keep[y[keep] == g], , drop = FALSE], 2, stats::median)
    })
    if (is.null(dim(med))) med <- matrix(med, ncol = length(GROUPS))
    colnames(med) <- GROUPS
    top_i <- GROUPS[max.col(med, ties.method = "first")]
    top_same[i, ] <- top_i == full_top
    delta_arr[i, , "FD_minus_NP"] <- med[, "FD"] - med[, "NP"]
    delta_arr[i, , "SO_minus_NP"] <- med[, "SO"] - med[, "NP"]
    delta_arr[i, , "SO_minus_FD"] <- med[, "SO"] - med[, "FD"]
  }

  out <- data.frame(
    feature_model = colnames(X),
    jackknife_top_stability = colMeans(top_same, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  for (k in colnames(full_d)) {
    d <- delta_arr[, , k, drop = FALSE][, , 1L]
    full_sign <- sign(full_d[, k])
    same <- sweep(sign(d), 2, full_sign, "==")
    same[, full_sign == 0] <- FALSE
    out[[paste0("jk_sign_stability__", k)]] <- colMeans(same, na.rm = TRUE)
    out[[paste0("jk_prob_gt0__", k)]] <- colMeans(d > 0, na.rm = TRUE)
    out[[paste0("jk_prob_lt0__", k)]] <- colMeans(d < 0, na.rm = TRUE)
  }
  out
}

classify_pattern <- function(df) {
  p_fd_np <- df$jk_prob_gt0__FD_minus_NP
  p_so_np <- df$jk_prob_gt0__SO_minus_NP
  p_so_fd <- df$jk_prob_gt0__SO_minus_FD
  top <- df$top_location_group
  high <- df$posterior_high_group
  stab <- df$jackknife_top_stability

  pattern <- rep("uncertain", nrow(df))
  active <- !df$constant_feature
  shared <- active & p_fd_np >= JK_B & p_so_np >= JK_B & p_so_fd > 0.25 & p_so_fd < 0.75
  pattern[shared] <- "shared_FD_SO_high"

  for (g in GROUPS) {
    idx <- active & pattern == "uncertain" & top == g & high == g & stab >= JK_B
    pattern[idx] <- paste0(g, "_like_high")
  }
  idx <- active & pattern == "uncertain" & p_so_fd >= JK_A
  pattern[idx] <- "SO_gt_FD_splitter"
  idx <- active & pattern == "uncertain" & p_so_fd <= (1 - JK_A)
  pattern[idx] <- "FD_gt_SO_splitter"
  pattern
}

evidence_grade <- function(df) {
  min_pair_q <- pmin(
    df$q_pair__FD_minus_NP,
    df$q_pair__SO_minus_NP,
    df$q_pair__SO_minus_FD,
    na.rm = TRUE
  )
  pattern_ok <- df$bidirectional_pattern != "uncertain"
  strong <- is.finite(df$q_global) & df$q_global < Q_A &
    is.finite(df$MI_norm) & df$MI_norm >= MI_A &
    is.finite(df$jackknife_top_stability) & df$jackknife_top_stability >= JK_A &
    df$concordant_high %in% TRUE &
    is.finite(df$posterior_enrichment_high) & df$posterior_enrichment_high > 0 & pattern_ok
  moderate <- !strong & is.finite(df$q_global) & df$q_global < Q_B &
    is.finite(df$jackknife_top_stability) & df$jackknife_top_stability >= JK_B & pattern_ok
  exploratory <- !strong & !moderate & pattern_ok &
    ((is.finite(df$p_global) & df$p_global < P_C) | (is.finite(min_pair_q) & min_pair_q < PAIR_Q_SUPPORT))
  dependent <- !strong & !moderate & !exploratory &
    ((is.finite(df$q_global) & df$q_global < Q_B) |
       (is.finite(df$p_global) & df$p_global < P_C) |
       (is.finite(df$MI_norm) & df$MI_norm >= MI_DEPENDENCY))
  out <- rep("E_no_evidence", nrow(df))
  out[dependent] <- "D_ambiguous_dependent"
  out[exploratory] <- "C_exploratory"
  out[moderate] <- "B_moderate"
  out[strong] <- "A_strong"
  out
}

run_train_framework <- function(X, y, bank, signal_label) {
  priors <- make_priors(y, PRIOR_MODE)
  glob <- global_rank_equal_prior_permutation(X, y, bank)
  pair <- pairwise_exact_rank_permutation(X, y)
  desc <- student_t_descriptors(X, y, priors)
  jk <- jackknife_direction(X, y, desc)

  out <- Reduce(function(a, b) merge(a, b, by = "feature_model", all = TRUE, sort = FALSE),
                list(desc, glob, pair, jk))
  out$bidirectional_pattern <- classify_pattern(out)
  out$evidence_grade <- evidence_grade(out)
  out$signal <- signal_label
  grade_order <- c("A_strong", "B_moderate", "C_exploratory", "D_ambiguous_dependent", "E_no_evidence")
  out$grade_rank <- match(out$evidence_grade, grade_order)
  out <- out[order(out$grade_rank, out$q_global, -out$MI_norm, -pmax(out$JS_FD_NP, out$JS_SO_NP, out$JS_SO_FD)), ]
  rownames(out) <- NULL
  out
}

# ---------------------------- TRAIN obs vs rec ------------------------------
prefix_except <- function(df, prefix, except = "feature_model") {
  n <- names(df)
  names(df)[!(n %in% except)] <- paste0(prefix, n[!(n %in% except)])
  df
}

build_final_train_table <- function(feature_map, rec, obs) {
  r <- prefix_except(rec, "rec_")
  o <- prefix_except(obs, "obs_")
  x <- merge(feature_map, r, by = "feature_model", all.x = TRUE, sort = FALSE)
  x <- merge(x, o, by = "feature_model", all.x = TRUE, sort = FALSE)

  contrast_cols <- c("FD_minus_NP", "SO_minus_NP", "SO_minus_FD")
  same <- sapply(contrast_cols, function(k) {
    sign(x[[paste0("rec_delta_", k)]]) == sign(x[[paste0("obs_delta_", k)]]) &
      sign(x[[paste0("rec_delta_", k)]]) != 0
  })
  if (is.null(dim(same))) same <- matrix(same, ncol = length(contrast_cols))
  x$direction_concordance_n <- rowSums(same, na.rm = TRUE)
  x$direction_concordance_prop <- x$direction_concordance_n / length(contrast_cols)

  x$cross_signal_support <- ifelse(
    x$direction_concordance_n < DIRECTION_CONCORDANCE_MIN,
    "discordant_or_weak_direction",
    ifelse(
      x$obs_evidence_grade %in% c("A_strong", "B_moderate", "C_exploratory"),
      "concordant_with_observed_statistical_support",
      "concordant_descriptive_observed_support"
    )
  )

  x$priority_tier <- "T4_not_selected"
  x$priority_tier[
    x$rec_evidence_grade %in% c("A_strong", "B_moderate") &
      x$direction_concordance_n >= DIRECTION_CONCORDANCE_MIN &
      x$obs_evidence_grade %in% c("A_strong", "B_moderate", "C_exploratory")
  ] <- "T1_integrated_concordant"
  x$priority_tier[
    x$priority_tier == "T4_not_selected" &
      x$rec_evidence_grade %in% c("A_strong", "B_moderate") &
      x$direction_concordance_n >= DIRECTION_CONCORDANCE_MIN
  ] <- "T2_reconstructed_strong_direction_concordant"
  x$priority_tier[
    x$priority_tier == "T4_not_selected" &
      x$rec_evidence_grade == "C_exploratory" &
      x$direction_concordance_n >= DIRECTION_CONCORDANCE_MIN
  ] <- "T3_exploratory_concordant"
  x$biomarker_candidate_primary <- x$priority_tier %in% c(
    "T1_integrated_concordant", "T2_reconstructed_strong_direction_concordant"
  )
  x$biomarker_candidate_exploratory <- x$priority_tier %in% c(
    "T1_integrated_concordant", "T2_reconstructed_strong_direction_concordant",
    "T3_exploratory_concordant"
  )
  # Alias conservado para downstreams simples: por defecto = conjunto primario.
  x$biomarker_candidate_default <- x$biomarker_candidate_primary
  # -------------------------------------------------------------------------
  # Pools de alta confianza para downstream.
  #
  # No modifican T1/T2/T3 ni evidence_grade.
  # Se derivan exclusivamente de TRAIN.
  # -------------------------------------------------------------------------
  
  x$rec_obs_same_pattern <-
    !is.na(x$rec_bidirectional_pattern) &
    !is.na(x$obs_bidirectional_pattern) &
    x$rec_bidirectional_pattern == x$obs_bidirectional_pattern &
    x$rec_bidirectional_pattern != "uncertain"
  
  x$high_confidence_AA_same_pattern <-
    x$rec_evidence_grade == "A_strong" &
    x$obs_evidence_grade == "A_strong" &
    x$rec_obs_same_pattern &
    x$direction_concordance_n >= DIRECTION_CONCORDANCE_MIN
  
  class_like_patterns <- c(
    "NP_like_high",
    "FD_like_high",
    "SO_like_high"
  )
  
  x$high_confidence_class_signature <-
    x$high_confidence_AA_same_pattern &
    x$rec_bidirectional_pattern %in% class_like_patterns
  
  x$high_confidence_class <- NA_character_
  
  x$high_confidence_class[
    x$high_confidence_class_signature
  ] <- sub(
    "_like_high$",
    "",
    x$rec_bidirectional_pattern[
      x$high_confidence_class_signature
    ]
  )
  # FDR por vista como auditoria adicional. La gradacion primaria conserva FDR global
  # entre todas las features del analisis para no cambiar silenciosamente la familia
  # de hipotesis.
  for (prefix in c("rec_", "obs_")) {
    pcols <- c(
      "p_global",
      "p_pair__FD_minus_NP",
      "p_pair__SO_minus_NP",
      "p_pair__SO_minus_FD"
    )
    for (pc in pcols) {
      src <- paste0(prefix, pc)
      if (!src %in% names(x)) next
      dst <- paste0(prefix, sub("^p_", "q_within_view_", pc))
      x[[dst]] <- ave(
        x[[src]], x$view_code,
        FUN = function(z) stats::p.adjust(z, method = "BH")
      )
    }
  }

  tier_order <- c(
    "T1_integrated_concordant", "T2_reconstructed_strong_direction_concordant",
    "T3_exploratory_concordant", "T4_not_selected"
  )
  x <- x[order(match(x$priority_tier, tier_order), x$rec_grade_rank, x$rec_q_global, -x$rec_MI_norm), ]
  rownames(x) <- NULL
  x
}

# ---------------------------- TEST congelado --------------------------------
posterior_matrix_frozen <- function(x, row, df_t = DF_T) {
  mu <- c(NP = row$mu_NP, FD = row$mu_FD, SO = row$mu_SO)
  sigma <- c(NP = row$sigma_NP, FD = row$sigma_FD, SO = row$sigma_SO)
  priors <- c(NP = row$prior_NP, FD = row$prior_FD, SO = row$prior_SO)
  out <- t(vapply(x, function(x0) posterior_at_x(x0, mu, sigma, priors, df_t), numeric(length(GROUPS))))
  colnames(out) <- GROUPS
  out
}

validate_test_frozen <- function(train_core, X_test, y_test, signal_label) {
  X_test <- as.matrix(X_test)
  y_test <- factor(y_test, levels = GROUPS)
  feature_summaries <- vector("list", ncol(X_test))
  pred_rows <- vector("list", ncol(X_test))

  core_idx <- match(colnames(X_test), train_core$feature_model)
  assert(!anyNA(core_idx), "TEST contiene features sin parametros TRAIN")

  for (j in seq_len(ncol(X_test))) {
    feat <- colnames(X_test)[j]
    row <- train_core[core_idx[j], , drop = FALSE]
    P <- posterior_matrix_frozen(X_test[, j], row)
    pred <- factor(colnames(P)[max.col(P, ties.method = "first")], levels = GROUPS)
    true_idx <- cbind(seq_len(nrow(P)), match(as.character(y_test), GROUPS))
    p_true <- P[true_idx]
    onehot <- model.matrix(~ y_test - 1)
    colnames(onehot) <- sub("^y_test", "", colnames(onehot))
    onehot <- onehot[, GROUPS, drop = FALSE]
    logloss <- -mean(log(pmax(p_true, 1e-15)))
    brier <- mean(rowSums((P - onehot)^2))
    prior_true <- vapply(as.character(y_test), function(g) row[[paste0("prior_", g)]], numeric(1))
    enrich_true <- log(pmax(p_true, 1e-15) / pmax(prior_true, 1e-15))

    med_test <- setNames(vapply(GROUPS, function(g) stats::median(X_test[y_test == g, j]), numeric(1)), GROUPS)
    d_test <- c(
      FD_minus_NP = unname(med_test["FD"] - med_test["NP"]),
      SO_minus_NP = unname(med_test["SO"] - med_test["NP"]),
      SO_minus_FD = unname(med_test["SO"] - med_test["FD"])
    )
    d_train <- c(
      FD_minus_NP = row$delta_FD_minus_NP,
      SO_minus_NP = row$delta_SO_minus_NP,
      SO_minus_FD = row$delta_SO_minus_FD
    )
    rep_dir <- sign(d_test) == sign(d_train) & sign(d_train) != 0

    feature_summaries[[j]] <- data.frame(
      feature_model = feat,
      signal = signal_label,
      test_mean_p_true = mean(p_true),
      test_median_p_true = stats::median(p_true),
      test_mean_log_posterior_enrichment_true = mean(enrich_true),
      test_logloss = logloss,
      test_brier = brier,
      test_accuracy = mean(pred == y_test),
      test_balanced_accuracy = balanced_accuracy(y_test, pred),
      test_macro_auc_ovr = macro_auc_ovr(y_test, P),
      test_delta_FD_minus_NP = d_test["FD_minus_NP"],
      test_delta_SO_minus_NP = d_test["SO_minus_NP"],
      test_delta_SO_minus_FD = d_test["SO_minus_FD"],
      direction_replicated_FD_minus_NP = rep_dir["FD_minus_NP"],
      direction_replicated_SO_minus_NP = rep_dir["SO_minus_NP"],
      direction_replicated_SO_minus_FD = rep_dir["SO_minus_FD"],
      direction_replicated_n = sum(rep_dir, na.rm = TRUE),
      direction_replicated_prop = mean(rep_dir, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    pred_rows[[j]] <- data.frame(
      sample = rownames(X_test),
      true_group = as.character(y_test),
      feature_model = feat,
      signal = signal_label,
      x_scaled_using_train = X_test[, j],
      p_NP = P[, "NP"], p_FD = P[, "FD"], p_SO = P[, "SO"],
      predicted_group = as.character(pred),
      p_true = p_true,
      log_posterior_enrichment_true = enrich_true,
      stringsAsFactors = FALSE
    )
  }

  list(
    feature_summary = do.call(rbind, feature_summaries),
    predictions = do.call(rbind, pred_rows)
  )
}

ensemble_from_frozen_predictions <- function(predictions, selected_features) {
  
  empty_result <- function() {
    list(
      summary = data.frame(),
      predictions = data.frame(),
      confusion = data.frame(),
      confusion_normalized = data.frame(),
      class_metrics = data.frame()
    )
  }
  
  if (!length(selected_features)) {
    return(empty_result())
  }
  
  x <- predictions[
    predictions$feature_model %in% selected_features,
    ,
    drop = FALSE
  ]
  
  if (!nrow(x)) {
    return(empty_result())
  }
  
  # -------------------------------------------------------------------------
  # Pooling geometrico de probabilidades de las features seleccionadas en TRAIN.
  # TEST no reajusta pesos ni parametros.
  # -------------------------------------------------------------------------
  samples <- unique(x$sample)
  
  rows <- lapply(samples, function(id) {
    
    d <- x[x$sample == id, , drop = FALSE]
    
    probs <- as.matrix(
      d[, c("p_NP", "p_FD", "p_SO"), drop = FALSE]
    )
    colnames(probs) <- GROUPS
    
    logp <- colMeans(log(pmax(probs, 1e-15)))
    
    p <- exp(logp - max(logp))
    p <- p / sum(p)
    
    truth <- d$true_group[1]
    pred <- GROUPS[which.max(p)]
    
    data.frame(
      sample = id,
      true_group = truth,
      p_NP = p["NP"],
      p_FD = p["FD"],
      p_SO = p["SO"],
      predicted_group = pred,
      n_features = nrow(d),
      stringsAsFactors = FALSE
    )
  })
  
  pred <- do.call(rbind, rows)
  
  P <- as.matrix(
    pred[, c("p_NP", "p_FD", "p_SO"), drop = FALSE]
  )
  colnames(P) <- GROUPS
  
  truth <- factor(pred$true_group, levels = GROUPS)
  pp <- factor(pred$predicted_group, levels = GROUPS)
  
  ptrue <- P[
    cbind(
      seq_len(nrow(P)),
      match(pred$true_group, GROUPS)
    )
  ]
  
  Y <- model.matrix(~ truth - 1)
  colnames(Y) <- sub("^truth", "", colnames(Y))
  Y <- Y[, GROUPS, drop = FALSE]
  
  # -------------------------------------------------------------------------
  # Matriz de confusion:
  # filas    = clase REAL
  # columnas = clase PREDICHA
  # -------------------------------------------------------------------------
  cm <- table(
    reference = truth,
    predicted = pp
  )
  
  confusion <- as.data.frame(
    cm,
    stringsAsFactors = FALSE
  )
  
  # -------------------------------------------------------------------------
  # Matriz normalizada POR CLASE REAL.
  # Cada fila suma 1.
  # Permite ver, por ejemplo, que proporcion de FD reales se predice
  # correctamente como FD, NP o SO.
  # -------------------------------------------------------------------------
  cm_norm <- matrix(
    NA_real_,
    nrow = nrow(cm),
    ncol = ncol(cm),
    dimnames = dimnames(cm)
  )
  
  for (i in seq_len(nrow(cm))) {
    den <- sum(cm[i, ])
    
    if (den > 0) {
      cm_norm[i, ] <- cm[i, ] / den
    }
  }
  
  confusion_normalized <- as.data.frame(
    as.table(cm_norm),
    stringsAsFactors = FALSE
  )
  
  names(confusion_normalized) <- c(
    "reference",
    "predicted",
    "proportion"
  )
  
  # -------------------------------------------------------------------------
  # Metricas one-vs-rest por clase.
  # -------------------------------------------------------------------------
  safe_ratio <- function(num, den) {
    if (!is.finite(den) || den <= 0) {
      return(NA_real_)
    }
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
      
      f1 <- if (
        is.finite(precision) &&
        is.finite(recall) &&
        (precision + recall) > 0
      ) {
        2 * precision * recall / (precision + recall)
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
  
  # -------------------------------------------------------------------------
  # Metricas globales de la firma congelada.
  # -------------------------------------------------------------------------
  summary <- data.frame(
    method = "equal_weight_log_pool_of_T1_T2_train_selected_features",
    analysis_role = "exploratory_transportability_audit_not_final_predictor",
    n_features = length(selected_features),
    n_test = nrow(pred),
    
    accuracy = mean(pp == truth),
    
    balanced_accuracy = balanced_accuracy(
      truth,
      pp
    ),
    
    macro_auc_ovr = macro_auc_ovr(
      truth,
      P
    ),
    
    macro_precision = mean(
      class_metrics$precision,
      na.rm = TRUE
    ),
    
    macro_recall = mean(
      class_metrics$recall,
      na.rm = TRUE
    ),
    
    macro_F1 = mean(
      class_metrics$F1,
      na.rm = TRUE
    ),
    
    logloss = -mean(
      log(pmax(ptrue, 1e-15))
    ),
    
    brier = mean(
      rowSums((P - Y)^2)
    ),
    
    mean_p_true = mean(ptrue),
    
    stringsAsFactors = FALSE
  )
  
  list(
    summary = summary,
    predictions = pred,
    confusion = confusion,
    confusion_normalized = confusion_normalized,
    class_metrics = class_metrics
  )
}

annotate_transcriptomics <- function(feature_map, rfm) {
  out <- feature_map
  out$GeneSymbol <- NA_character_
  out$EntrezGeneID <- NA_character_
  out$EnsemblID <- NA_character_
  out$Description <- NA_character_
  if (is.null(rfm$features_metadata) || !is.data.frame(rfm$features_metadata)) return(out)
  md <- rfm$features_metadata
  tx_idx <- which(out$view_code == "tx")
  if (!length(tx_idx)) return(out)
  candidate_cols <- intersect(
    c("GeneSymbol", "EntrezGeneID", "ProbeID", "EnsemblID", "PrimaryAccession", "RefSeqAccession"),
    names(md)
  )
  if (!length(candidate_cols)) return(out)

  for (ii in tx_idx) {
    key <- as.character(out$feature_original[ii])
    hit <- NA_integer_
    for (cc in candidate_cols) {
      mm <- match(key, as.character(md[[cc]]))
      if (!is.na(mm)) { hit <- mm; break }
    }
    if (!is.na(hit)) {
      if ("GeneSymbol" %in% names(md)) out$GeneSymbol[ii] <- as.character(md$GeneSymbol[hit])
      if ("EntrezGeneID" %in% names(md)) out$EntrezGeneID[ii] <- as.character(md$EntrezGeneID[hit])
      if ("EnsemblID" %in% names(md)) out$EnsemblID[ii] <- as.character(md$EnsemblID[hit])
      if ("Description" %in% names(md)) out$Description[ii] <- as.character(md$Description[hit])
    }
  }
  out
}


# ---------------------------- Diagnosticos rec -------------------------------
reconstruction_diagnostics <- function(data) {
  fmap <- data$feature_map
  rows <- lapply(names(VIEW_MAP), function(v) {
    code <- VIEW_MAP[[v]]
    feats <- fmap$feature_model[fmap$view_code == code]
    Otr <- data$observed_train[, feats, drop = FALSE]
    Rtr <- data$reconstructed_train[, feats, drop = FALSE]
    Ote <- data$observed_test[, feats, drop = FALSE]
    Rte <- data$reconstructed_test[, feats, drop = FALSE]
    calc <- function(O, R) {
      sse <- sum((O - R)^2)
      sst <- sum(sweep(O, 2, colMeans(O), "-")^2)
      c(RMSE = sqrt(mean((O - R)^2)), R2 = if (sst > NUM_TOL) 1 - sse / sst else NA_real_)
    }
    tr <- calc(Otr, Rtr); te <- calc(Ote, Rte)
    data.frame(
      view = v, view_code = code, n_features = length(feats),
      train_reconstruction_RMSE = tr["RMSE"], train_reconstruction_R2 = tr["R2"],
      test_reconstruction_RMSE = te["RMSE"], test_reconstruction_R2 = te["R2"],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# ------------------------------- Graficos -----------------------------------
open_plot_device <- function(path, type = c("pdf", "png"), width = 10, height = 7) {
  type <- match.arg(type)
  if (type == "pdf") {
    grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  } else {
    grDevices::png(path, width = width, height = height, units = "in", res = 150,
                   type = if (capabilities("cairo")) "cairo" else getOption("bitmapType", "Xlib"))
  }
}

make_family_plots <- function(final_train, test_rec, outdir, family) {
  if (!MAKE_PLOTS) return(invisible(NULL))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  for (type in c("pdf", "png")) {
    ext <- type
    # 1) Evidence grade by view.
    path <- file.path(outdir, paste0("01_RECONSTRUCTED_EVIDENCE_BY_VIEW.", ext))
    open_plot_device(path, type)
    tab <- table(final_train$view_code, final_train$rec_evidence_grade)
    graphics::barplot(t(tab), beside = FALSE, las = 1,
                      main = paste(family, "| Evidencia TRAIN reconstruida por vista"),
                      xlab = "Vista", ylab = "Numero de features")
    graphics::legend("topright", legend = colnames(tab), fill = seq_len(ncol(tab)), bty = "n", cex = 0.75)
    grDevices::dev.off()

    # 2) MI TRAIN vs probabilidad clase verdadera TEST.
    path <- file.path(outdir, paste0("02_TRAIN_MI_VS_TEST_PTRUE.", ext))
    open_plot_device(path, type)
    idx <- match(final_train$feature_model, test_rec$feature_model)
    y <- test_rec$test_mean_p_true[idx]
    col <- ifelse(final_train$biomarker_candidate_default, "black", "grey70")
    graphics::plot(final_train$rec_MI_norm, y, pch = 16, col = col,
                   xlab = "MI normalizada TRAIN (reconstruido)",
                   ylab = "P media de clase verdadera en TEST",
                   main = paste(family, "| Transporte congelado TRAIN -> TEST"))
    graphics::abline(h = 1 / 3, lty = 2)
    graphics::legend("bottomright", c("Candidato TRAIN T1/T2", "Resto", "Prior igual 1/3"),
                     col = c("black", "grey70", "black"), pch = c(16, 16, NA), lty = c(NA, NA, 2), bty = "n")
    grDevices::dev.off()

    # 3) Concordancia obs/reconstruido en TRAIN.
    path <- file.path(outdir, paste0("03_OBSERVED_RECONSTRUCTED_DIRECTION_CONCORDANCE.", ext))
    open_plot_device(path, type)
    tb <- table(final_train$direction_concordance_n)
    graphics::barplot(tb, main = paste(family, "| Concordancia de direccion observado vs reconstruido"),
                      xlab = "Numero de contrastes concordantes (de 3)", ylab = "Features")
    grDevices::dev.off()
  }
}

# ----------------------------- Ejecucion familia ----------------------------
run_family <- function(selection_run, family, stamp, parameters) {
  message("\n============================================================")
  message("FAMILIA ", family)
  message("============================================================")
  selected <- read_selected_row(selection_run, family)
  data <- assemble_family_data(selection_run, family, selected)
  message("Escenario: ", data$scenario)
  message("Features comunes MOFA/RFM: ", nrow(data$feature_map))
  message("TRAIN: ", length(data$y_train), " | TEST: ", length(data$y_test))

  base_out <- file.path(selection_run, family, "BIDIRECTIONAL_FRAMEWORK")
  out <- file.path(base_out, stamp)
  assert(!dir.exists(out), paste("La salida ya existe:", out))
  dirs <- c("tables", "rds", "plots")
  for (d in dirs) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)

  # Anotacion opcional para TX.
  data$feature_map <- annotate_transcriptomics(data$feature_map, data$rfm)

  # Escalado se aprende solo en TRAIN; TEST usa esos parametros congelados.
  scaler_rec <- fit_scaler(data$reconstructed_train, SCALE_METHOD)
  scaler_obs <- fit_scaler(data$observed_train, SCALE_METHOD)
  Xrec_tr <- apply_scaler(data$reconstructed_train, scaler_rec)
  Xrec_te <- apply_scaler(data$reconstructed_test, scaler_rec)
  Xobs_tr <- apply_scaler(data$observed_train, scaler_obs)
  Xobs_te <- apply_scaler(data$observed_test, scaler_obs)

  bank <- make_global_bank(data$y_train, N_PERM_GLOBAL, SEED + match(family, FAMILIES) * 1000L)

  message("Framework TRAIN reconstruido...")
  core_rec <- run_train_framework(Xrec_tr, data$y_train, bank, "reconstructed_train")
  message("Framework TRAIN observado...")
  core_obs <- run_train_framework(Xobs_tr, data$y_train, bank, "observed_train")

  final_train <- build_final_train_table(data$feature_map, core_rec, core_obs)

  message("Validacion TEST congelada: reconstruido...")
  test_rec <- validate_test_frozen(core_rec, Xrec_te, data$y_test, "reconstructed_test_frozen")
  message("Validacion TEST congelada: observado...")
  test_obs <- validate_test_frozen(core_obs, Xobs_te, data$y_test, "observed_test_frozen")

  # TEST nunca decide estas features: candidate_default se fijo antes.
  selected_features <- final_train$feature_model[final_train$biomarker_candidate_default]
  ensemble <- ensemble_from_frozen_predictions(test_rec$predictions, selected_features)


  # Diagnosticos de proyeccion/reconstruccion.
  rec_diag <- reconstruction_diagnostics(data)

  # Resumen por vista.
  view_summary <- do.call(rbind, lapply(unique(final_train$view_code), function(v) {
    z <- final_train[final_train$view_code == v, , drop = FALSE]
    data.frame(
      view_code = v,
      view = unique(z$view)[1],
      n_features = nrow(z),
      n_A_rec = sum(z$rec_evidence_grade == "A_strong", na.rm = TRUE),
      n_B_rec = sum(z$rec_evidence_grade == "B_moderate", na.rm = TRUE),
      n_C_rec = sum(z$rec_evidence_grade == "C_exploratory", na.rm = TRUE),
      n_A_obs = sum(
        z$obs_evidence_grade == "A_strong",
        na.rm = TRUE
      ),
      
      n_high_confidence_AA_same_pattern = sum(
        z$high_confidence_AA_same_pattern,
        na.rm = TRUE
      ),
      
      n_high_confidence_class_signature = sum(
        z$high_confidence_class_signature,
        na.rm = TRUE
      ),
      n_T1 = sum(z$priority_tier == "T1_integrated_concordant", na.rm = TRUE),
      n_T2 = sum(z$priority_tier == "T2_reconstructed_strong_direction_concordant", na.rm = TRUE),
      n_default_biomarkers = sum(z$biomarker_candidate_default, na.rm = TRUE),
      median_MI_rec = stats::median(z$rec_MI_norm, na.rm = TRUE),
      median_direction_concordance = stats::median(z$direction_concordance_prop, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  # Unir test summary al listado final sin alterar selection.
  final_with_test <- final_train
  test_rec_pref <- prefix_except(test_rec$feature_summary, "testrec_")
  test_obs_pref <- prefix_except(test_obs$feature_summary, "testobs_")
  final_with_test <- merge(final_with_test, test_rec_pref, by = "feature_model", all.x = TRUE, sort = FALSE)
  final_with_test <- merge(final_with_test, test_obs_pref, by = "feature_model", all.x = TRUE, sort = FALSE)
  tier_order <- c("T1_integrated_concordant", "T2_reconstructed_strong_direction_concordant", "T3_exploratory_concordant", "T4_not_selected")
  final_with_test <- final_with_test[order(match(final_with_test$priority_tier, tier_order), final_with_test$rec_q_global, -final_with_test$rec_MI_norm), ]

  # Escribir tablas.
  tdir <- file.path(out, "tables")
  write_csv(data$feature_map, file.path(tdir, "00_FEATURE_MAP.csv"))
  write_csv(data$projection_audit, file.path(tdir, "01_TEST_LATENT_PROJECTION_AUDIT.csv"))
  write_csv(rec_diag, file.path(tdir, "02_RECONSTRUCTION_DIAGNOSTICS.csv"))
  write_csv(core_rec, file.path(tdir, "03_TRAIN_RECONSTRUCTED_FRAMEWORK_ALL_FEATURES.csv"))
  write_csv(core_obs, file.path(tdir, "04_TRAIN_OBSERVED_FRAMEWORK_ALL_FEATURES.csv"))
  write_csv(final_train, file.path(tdir, "05_TRAIN_BIDIRECTIONAL_INTEGRATED_SELECTION.csv"))
  write_csv(final_train[final_train$biomarker_candidate_primary, , drop = FALSE],
            file.path(tdir, "06_BIOMARKER_CANDIDATES_PRIMARY_TRAIN_ONLY.csv"))
  write_csv(final_train[final_train$biomarker_candidate_exploratory, , drop = FALSE],
            file.path(tdir, "06B_BIOMARKER_CANDIDATES_EXPLORATORY_TRAIN_ONLY.csv"))
  
  write_csv(
    final_train[
      final_train$high_confidence_AA_same_pattern,
      ,
      drop = FALSE
    ],
    file.path(
      tdir,
      "06C_HIGH_CONFIDENCE_AA_SAME_PATTERN_TRAIN_ONLY.csv"
    )
  )
  
  write_csv(
    final_train[
      final_train$high_confidence_class_signature,
      ,
      drop = FALSE
    ],
    file.path(
      tdir,
      "06D_HIGH_CONFIDENCE_CLASS_SIGNATURE_TRAIN_ONLY.csv"
    )
  )
  transcriptomic_universe <- final_train[
    final_train$view_code == "tx",
    ,
    drop = FALSE
  ]
  
  write_csv(
    transcriptomic_universe,
    file.path(
      tdir,
      "06E_TRANSCRIPTOMIC_UNIVERSE_TRAIN_ONLY.csv"
    )
  )
  write_csv(
    final_train[
      final_train$view_code == "tx" &
        final_train$high_confidence_class_signature,
      ,
      drop = FALSE
    ],
    file.path(
      tdir,
      "06F_TRANSCRIPTOMIC_HIGH_CONFIDENCE_CLASS_SIGNATURE.csv"
    )
  )
  write_csv(view_summary, file.path(tdir, "07_VIEW_SUMMARY_TRAIN.csv"))
  write_csv(test_rec$feature_summary, file.path(tdir, "08_TEST_RECONSTRUCTED_FROZEN_FEATURE_METRICS.csv"))
  write_csv(test_obs$feature_summary, file.path(tdir, "09_TEST_OBSERVED_FROZEN_FEATURE_METRICS.csv"))
  write_csv(final_with_test, file.path(tdir, "10_FINAL_FEATURE_TABLE_WITH_TEST_VALIDATION.csv"))
  write_csv(test_rec$predictions, file.path(tdir, "11_TEST_RECONSTRUCTED_FROZEN_POSTERIORS_LONG.csv"))
  write_csv(test_obs$predictions, file.path(tdir, "12_TEST_OBSERVED_FROZEN_POSTERIORS_LONG.csv"))
  write_csv(ensemble$summary, file.path(tdir, "13_TEST_BIOMARKER_SIGNATURE_ENSEMBLE_METRICS.csv"))
  write_csv(ensemble$predictions, file.path(tdir, "14_TEST_BIOMARKER_SIGNATURE_ENSEMBLE_PREDICTIONS.csv"))
  write_csv(ensemble$confusion, file.path(tdir, "15_TEST_BIOMARKER_SIGNATURE_ENSEMBLE_CONFUSION.csv"))
  write_csv(
    ensemble$confusion_normalized,
    file.path(
      tdir,
      "15B_TEST_BIOMARKER_SIGNATURE_ENSEMBLE_CONFUSION_NORMALIZED.csv"
    )
  )
  
  write_csv(
    ensemble$class_metrics,
    file.path(
      tdir,
      "15C_TEST_BIOMARKER_SIGNATURE_ENSEMBLE_CLASS_METRICS.csv"
    )
  )


  # RDS reproducibles.
  saveRDS(
    list(reconstructed = scaler_rec, observed = scaler_obs),
    file.path(out, "rds", "TRAIN_SCALERS_FROZEN.rds"), version = 3
  )
  saveRDS(bank, file.path(out, "rds", "TRAIN_PERMUTATION_BANK.rds"), version = 3)
  saveRDS(
    list(
      core_reconstructed_train = core_rec,
      core_observed_train = core_obs,
      final_train = final_train,
      high_confidence_all_pattern = final_train[
        final_train$high_confidence_AA_same_pattern,
        ,
        drop = FALSE
      ],
      
      high_confidence_class_signature = final_train[
        final_train$high_confidence_class_signature,
        ,
        drop = FALSE
      ],
      
      transcriptomic_universe = final_train[
        final_train$view_code == "tx",
        ,
        drop = FALSE
      ],
      test_reconstructed = test_rec,
      test_observed = test_obs,
      ensemble = ensemble,
      feature_map = data$feature_map,
      projection_audit = data$projection_audit,
      reconstruction_diagnostics = rec_diag,
      Z_train = data$Z_train,
      Z_test_projected = data$Z_test,
      selected_scenario_row = selected,
      model_path = data$model_path,
      rfm_path = data$rfm_path
    ),
    file.path(out, "rds", "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"), version = 3
  )

  saveRDS(parameters, file.path(out, "rds", "PARAMETERS.rds"), version = 3)
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "SESSION_INFO.txt"))

  script_arg <- commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))]
  if (length(script_arg)) {
    sp <- sub("^--file=", "", script_arg[1])
    if (file.exists(sp)) file.copy(sp, file.path(out, "SCRIPT_USED.R"), overwrite = FALSE)
  }

  make_family_plots(final_train, test_rec$feature_summary, file.path(out, "plots"), family)

  # Resumen legible.
  n_bio <- sum(final_train$biomarker_candidate_default, na.rm = TRUE)
  summary_lines <- c(
    paste("METATEST |", family, "| FRAMEWORK BIDIRECCIONAL TRAIN-ONLY + TEST CONGELADO"),
    paste("Escenario seleccionado:", data$scenario),
    paste("Modelo MOFA:", data$model_path),
    paste("ready_for_modeling:", data$rfm_path),
    paste("TRAIN:", length(data$y_train), "muestras |", paste(names(table(data$y_train)), table(data$y_train), collapse = ", ")),
    paste("TEST:", length(data$y_test), "muestras |", paste(names(table(data$y_test)), table(data$y_test), collapse = ", ")),
    paste("Features evaluadas:", nrow(final_train)),
    paste("Candidatos primarios T1/T2 TRAIN-only:", n_bio),
    paste("Candidatos exploratorios T1/T2/T3 TRAIN-only:", sum(final_train$biomarker_candidate_exploratory, na.rm = TRUE)),
    "",
    "SELECCION:",
    "- Scenario selection: heredada del downstream de seleccion; no se repite aqui.",
    "- Biomarker discovery/grade/pattern: SOLO TRAIN.",
    "- Evidencia principal: reconstruccion MOFA + Student-t con escala regularizada para baja N + MI/JS + separacion robusta de rangos con peso igual NP/FD/SO + permutacion TRAIN + jackknife.",
    "- Observado TRAIN: soporte de concordancia, no sustituye la reconstruccion integrada.",
    "- T1/T2/T3 se conservan como resultados inferenciales del framework y NO se redefinen.",
    "- Se exporta un pool de alta confianza A_reconstructed + A_observed + mismo patron.",
    "- La firma class-specific exige ademas NP_like_high, FD_like_high o SO_like_high.",
    "- Estos pools se definen exclusivamente con TRAIN.",
    "",
    "TEST:",
    "- TEST se proyecta con W del MOFA TRAIN congelado y todos sus factores.",
    "- Los parametros Student-t, priors, patrones, grades y seleccion T1/T2 NO se reajustan en TEST.",
    "- El ensemble T1/T2 en TEST se conserva como auditoria exploratoria de transportabilidad; NO constituye el predictor final.",
    "- Las probabilidades de las features se combinan mediante pooling geometrico de peso igual.",
    "- TEST reporta accuracy, balanced accuracy, macro-AUC OVR, macro-F1, log-loss, Brier y probabilidad de clase verdadera.",
    "- Se reportan matriz de confusion en conteos, matriz normalizada por clase real y precision/recall/F1/specificity por NP, FD y SO.",
    "LIMITACIONES:",
    "",
    "- N=16 TRAIN y FD=3: inferencia exploratoria; q-values no eliminan la incertidumbre estructural.",
    "- La proyeccion TEST es lineal ridge sobre W y combina vistas con peso igual.",
    "- El ensemble TEST usa pooling geometrico igualitario y se considera diagnostico de firma, no modelo causal.",
    "- Los umbrales de evidence_grade son operativos/editables y deben predeclararse antes de interpretar resultados.",
    "",
    "DOWNSTREAM:",
    "- Prediccion final: Modelamiento 3 sobre valores OBSERVADOS usando como pool candidato la firma class-specific A/A concordante.",
    "- FGSEA: se realizara en script separado usando TODOS los transcriptos elegibles y rankings direccionales observados TRAIN.",
    "- El universo transcriptomico completo se exporta sin filtrar por significacion.",
    "- TEST no se utilizara para seleccionar features, tamaños de firma ni configuraciones downstream.",
    paste("Salida:", out)
  )
  writeLines(summary_lines, file.path(out, "READ_ME_FIRST.txt"), useBytes = TRUE)
  writeLines(out, file.path(base_out, "LATEST_RUN.txt"), useBytes = TRUE)
  n_high_conf <- sum(
    final_train$high_confidence_AA_same_pattern,
    na.rm = TRUE
  )
  
  n_high_conf_class <- sum(
    final_train$high_confidence_class_signature,
    na.rm = TRUE
  )
  status <- data.frame(
    family = family,
    scenario = data$scenario,
    n_train = length(data$y_train),
    n_test = length(data$y_test),
    n_features = nrow(final_train),
    n_biomarker_candidates = n_bio,
    n_high_confidence_AA_same_pattern = n_high_conf,
    n_high_confidence_class_signature = n_high_conf_class,
    n_biomarker_candidates_exploratory = sum(final_train$biomarker_candidate_exploratory, na.rm = TRUE),
    n_T1 = sum(final_train$priority_tier == "T1_integrated_concordant", na.rm = TRUE),
    n_T2 = sum(final_train$priority_tier == "T2_reconstructed_strong_direction_concordant", na.rm = TRUE),
    output = out,
    error = "",
    stringsAsFactors = FALSE
  )
  write_csv(status, file.path(out, "RUN_STATUS.csv"))
  status
}

# -------------------------------- main ---------------------------------------
main <- function() {
  if (!parse_cli()) return(invisible(NULL))
  assert(getRversion() >= "4.1", "Se requiere R >= 4.1")
  assert(length(FAMILIES) > 0L && !anyDuplicated(FAMILIES) && all(FAMILIES %in% c("P1", "P2", "P3")), "FAMILIES invalidas")
  assert(is.finite(SEED) && SEED >= 0, "SEED invalida")
  assert(is.finite(N_PERM_GLOBAL) && N_PERM_GLOBAL >= 99, "N_PERM_GLOBAL demasiado pequeno")
  assert(is.finite(N_GRID_DENSITY) && N_GRID_DENSITY >= 128, "N_GRID_DENSITY demasiado pequeno")
  assert(is.finite(DF_T) && DF_T > 0, "DF_T invalido")
  assert(SCALE_METHOD %in% c("robust", "zscore", "none"), "SCALE_METHOD invalido")
  assert(PRIOR_MODE %in% c("equal", "train_prevalence"), "PRIOR_MODE invalido")
  assert(is.finite(RIDGE_REL) && RIDGE_REL >= 0, "RIDGE_REL invalido")
  assert(is.finite(SCALE_SHRINKAGE_K) && SCALE_SHRINKAGE_K >= 0, "SCALE_SHRINKAGE_K invalido")
  assert(requireNamespace("MOFA2", quietly = TRUE), "Se requiere MOFA2 instalado para restaurar de forma segura el objeto MOFA")

  upstream <- canonical_path(UPSTREAM_DIR, TRUE)
  selection_run <- if (nzchar(SELECTION_RUN)) {
    canonical_path(SELECTION_RUN, TRUE)
  } else {
    find_latest_selection_run(upstream)
  }
  assert(file.exists(file.path(selection_run, "SELECTED_SCENARIOS_BY_ANALYSIS.csv")),
         "SELECTION_RUN no parece una corrida de seleccion completa")

  stamp <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "_pid", Sys.getpid())
  packages <- c("base", "stats", "utils", "graphics", "grDevices", "methods", "MOFA2")
  parameters <- list(
    script_version = "1.2.0-frozen",
    framework_status = "FROZEN",
    
    framework_role =
      "TRAIN_only_inference_and_candidate_characterization",
    
    downstream_prediction_pool =
      "A_reconstructed_AND_A_observed_AND_same_class_like_pattern",
    
    downstream_prediction_data =
      "observed_features",
    
    downstream_fgsea_universe =
      "all_transcriptomic_features_entering_selected_MOFA_and_framework",
    
    downstream_fgsea_primary_signal =
      "observed_TRAIN_directional_contrasts",
    
    test_ensemble_role =
      "exploratory_transportability_audit_not_final_predictor",
    date_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    R_version = R.version.string,
    package_versions = setNames(lapply(packages, pkg_version), packages),
    upstream = upstream,
    selection_run = selection_run,
    families = FAMILIES,
    groups = GROUPS,
    view_map = VIEW_MAP,
    seed = SEED,
    n_perm_global = N_PERM_GLOBAL,
    n_grid_density = N_GRID_DENSITY,
    df_t = DF_T,
    scale_method = SCALE_METHOD,
    prior_mode = PRIOR_MODE,
    ridge_rel = RIDGE_REL,
    scale_shrinkage_k = SCALE_SHRINKAGE_K,
    chunk_size = CHUNK_SIZE,
    
    thresholds = list(
      Q_A = Q_A,
      Q_B = Q_B,
      P_C = P_C,
      MI_A = MI_A,
      MI_DEPENDENCY = MI_DEPENDENCY,
      JK_A = JK_A,
      JK_B = JK_B,
      PAIR_Q_SUPPORT = PAIR_Q_SUPPORT,
      DIRECTION_CONCORDANCE_MIN = DIRECTION_CONCORDANCE_MIN
    ),
    
    selection_scope = "TRAIN_only",
    test_scope = "frozen_external_validation",
    test_signature_method = "equal_weight_log_pool_of_T1_T2_train_selected_features"
  )

  statuses <- list()
  for (family in FAMILIES) {
    statuses[[family]] <- tryCatch(
      run_family(selection_run, family, stamp, parameters),
      error = function(e) {
        family_root <- file.path(selection_run, family, "BIDIRECTIONAL_FRAMEWORK", stamp)
        dir.create(family_root, recursive = TRUE, showWarnings = FALSE)
        writeLines(conditionMessage(e), file.path(family_root, "FAMILY_ERROR.txt"))
        data.frame(
          family = family, scenario = NA_character_, n_train = NA_integer_, n_test = NA_integer_,
          n_features = NA_integer_, n_biomarker_candidates = NA_integer_,
          n_biomarker_candidates_exploratory = NA_integer_, n_T1 = NA_integer_, n_T2 = NA_integer_,
          output = family_root, error = conditionMessage(e), stringsAsFactors = FALSE
        )
      }
    )
  }
  status <- do.call(rbind, statuses)
  root_status <- file.path(selection_run, paste0("BIDIRECTIONAL_FRAMEWORK_RUN_", stamp, ".csv"))
  write_csv(status, root_status)

  cat("\nFRAMEWORK BIDIRECCIONAL TERMINADO\n")
  cat("Selection run: ", selection_run, "\n", sep = "")
  print(status, row.names = FALSE)
  if ("error" %in% names(status) && any(!is.na(status$error) & nzchar(status$error))) {
    stop("Alguna familia fallo. Revisar FAMILY_ERROR.txt en su salida.", call. = FALSE)
  }
  invisible(selection_run)
}

if (sys.nframe() == 0L) main()

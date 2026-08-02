rm(list = ls())

## ============================================================
## PANEL FINAL MULTIMODELO DESDE FEATURES FINALES
## ------------------------------------------------------------
## Entrada esperada:
##   - df_train_final.rds
##   - df_test_final.rds
##   - features_finales_integradas_mofa_brms.rds o .csv
##
## Objetivo:
##   Usar las features finales ya seleccionadas por MOFA/BRMS y
##   construir un panel final por consenso multimodelo.
##
## Modelos selectores:
##   1) sPLSDA
##   2) LASSO multinomial
##   3) PLSDA
##
## Modelo base de auditoría:
##   - Multinomial clásico
##
## Criterio final:
##   MOFA + al menos 2 de 3 modelos supervisados selectores/ranking
##
## BRMS confirmatorio:
##   1) Panel final completo
##
## TEST se usa solo como auditoría final.
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(caret)
  library(pROC)
  library(glmnet)
  library(mixOmics)
  library(nnet)
  library(doParallel)
  library(foreach)
  library(brms)
  library(posterior)
})
## Evitar que MOFA2 enmascare predict()
if ("package:MOFA2" %in% search()) {
  detach("package:MOFA2", unload = FALSE, character.only = TRUE)
}
set.seed(123)

## ============================================================
## 0) CONFIGURACIÓN PARAMETRIZADA
## Lee parámetros desde run_pipeline_from_template.sh
## ============================================================

get_env_chr <- function(name, default = NULL, required = FALSE) {
  val <- Sys.getenv(name, unset = NA_character_)
  
  if (is.na(val) || !nzchar(val)) {
    if (required) {
      stop("Variable de entorno obligatoria no definida: ", name)
    }
    return(default)
  }
  
  val
}

ANALYSIS_NAME <- get_env_chr("ANALYSIS_NAME", required = TRUE)
OUTDIR_BASE   <- get_env_chr("OUTDIR_BASE", ".")

MODEL_SOURCE <- toupper(get_env_chr("MODEL_SOURCE", "NEW"))

if (!MODEL_SOURCE %in% c("NEW", "OLD", "BOTH")) {
  stop(
    "MODEL_SOURCE debe ser NEW, OLD o BOTH. Valor actual: ",
    MODEL_SOURCE
  )
}

ANALYSIS_DIR <- file.path(OUTDIR_BASE, ANALYSIS_NAME)

## ============================================================
## Selección de entrada:
##   NEW  -> ANALYSIS_DIR/modelamiento
##   OLD  -> ANALYSIS_DIR/modelado_anterior
##   BOTH -> relanza el script para NEW y OLD
## ============================================================

if (MODEL_SOURCE == "BOTH") {
  
  cmd_args <- commandArgs(FALSE)
  script_file <- sub("^--file=", "", grep("^--file=", cmd_args, value = TRUE)[1])
  
  if (is.na(script_file) || !nzchar(script_file)) {
    stop(
      "MODEL_SOURCE=BOTH requiere ejecutar con Rscript panel_final_multimodelo.R. ",
      "Si estás usando source(), corre MODEL_SOURCE=NEW y MODEL_SOURCE=OLD por separado."
    )
  }
  
  for (src in c("NEW", "OLD")) {
    
    cat("\n============================================================\n")
    cat("RELANZANDO PANEL FINAL PARA MODEL_SOURCE =", src, "\n")
    cat("============================================================\n")
    
    status <- system2(
      command = file.path(R.home("bin"), "Rscript"),
      args = shQuote(script_file),
      env = c(
        paste0("ANALYSIS_NAME=", ANALYSIS_NAME),
        paste0("OUTDIR_BASE=", OUTDIR_BASE),
        paste0("MODEL_SOURCE=", src)
      )
    )
    
    if (!identical(status, 0L)) {
      stop("Falló panel final para MODEL_SOURCE = ", src)
    }
  }
  
  quit(save = "no", status = 0)
}

MODEL_DIR_NAME <- dplyr::case_when(
  MODEL_SOURCE == "NEW" ~ "modelamiento",
  MODEL_SOURCE == "OLD" ~ "modelado_anterior",
  TRUE ~ NA_character_
)

PANEL_DIR_NAME <- dplyr::case_when(
  MODEL_SOURCE == "NEW" ~ "panel_final_multimodelo",
  MODEL_SOURCE == "OLD" ~ "panel_final_multimodelo_anterior",
  TRUE ~ NA_character_
)

## Entradas desde el script de modelamiento BRMS/MOFA
INPUT_MODELAMIENTO_DIR <- file.path(ANALYSIS_DIR, MODEL_DIR_NAME)
INPUT_RDS_DIR          <- file.path(INPUT_MODELAMIENTO_DIR, "rds")
INPUT_TABLES_DIR       <- file.path(INPUT_MODELAMIENTO_DIR, "tables")

## Salidas propias de este script
PANEL_DIR  <- file.path(ANALYSIS_DIR, PANEL_DIR_NAME)
TABLES_DIR <- file.path(PANEL_DIR, "tables")
RDS_DIR    <- file.path(PANEL_DIR, "rds")
PLOTS_DIR  <- file.path(PANEL_DIR, "plots")

dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PANEL_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR,      recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR,    recursive = TRUE, showWarnings = FALSE)

out_file <- function(...) {
  file.path(PANEL_DIR, ...)
}

out_table <- function(...) {
  file.path(TABLES_DIR, ...)
}

out_rds <- function(...) {
  file.path(RDS_DIR, ...)
}

out_plot <- function(...) {
  file.path(PLOTS_DIR, ...)
}
REF_CLASS <- "NP"

cat("\n============================================================\n")
cat("PANEL FINAL MULTIMODELO\n")
cat("============================================================\n")
cat("ANALYSIS_NAME:", ANALYSIS_NAME, "\n")
cat("OUTDIR_BASE:", OUTDIR_BASE, "\n")
cat("ANALYSIS_DIR:", ANALYSIS_DIR, "\n")
cat("MODEL_SOURCE:", MODEL_SOURCE, "\n")
cat("MODEL_DIR_NAME:", MODEL_DIR_NAME, "\n")
cat("PANEL_DIR_NAME:", PANEL_DIR_NAME, "\n")
cat("INPUT_MODELAMIENTO_DIR:", INPUT_MODELAMIENTO_DIR, "\n")
cat("INPUT_RDS_DIR:", INPUT_RDS_DIR, "\n")
cat("INPUT_TABLES_DIR:", INPUT_TABLES_DIR, "\n")
cat("PANEL_DIR:", PANEL_DIR, "\n")
cat("TABLES_DIR:", TABLES_DIR, "\n")
cat("RDS_DIR:", RDS_DIR, "\n")
cat("PLOTS_DIR:", PLOTS_DIR, "\n")
## CV
K_CV <- 3
N_REPEATS_CV <- 10

## Semillas para estabilidad
SEEDS_SPLSDA <- 1:30
SEEDS_GLMNET <- 1:30

## Umbral de frecuencia para que una feature pase cada modelo
FREQ_THR_MODEL <- 0.60

## Umbral MOFA dentro de vista.
## 0.50 = mitad superior del score MOFA dentro de su vista.
MOFA_PERCENTILE_THR <- 0.50

## Modelos densos: PLSDA no es sparse estricto.
## Se toma el top 35% de features por importancia.
TOP_PROP_DENSE_MODELS <- 0.35

## sPLSDA
NCOMP_SPLSDA <- 2
KEEPX_GRID <- seq(5, 70, 5)

## PLSDA
NCOMP_PLSDA <- 2

## Elastic/Ridge
LAMBDA_GRID <- 10^seq(-3, 1, length.out = 30)
ALPHA_GRID_LASSO <- c(1.00)

## BRMS confirmatorio del panel final
DO_BRMS_CONFIRMATORY <- TRUE

BRMS_CONFIRM_MAX_FEATURES <- 25
BRMS_CONFIRM_MIN_FEATURES <- 3

## Segundo BRMS confirmatorio: top-N del panel final
BRMS_CONFIRM_TOP_N <- NA_integer_

## Ponderación del consenso.
## LASSO se mantiene, pero pesa menos por su mayor dispersión.
W_SPLSDA <- 1.0
W_LASSO  <- 0.5
W_PLSDA  <- 1.0



BRMS_CONFIRM_ITER <- 4000
BRMS_CONFIRM_CHAINS <- 4
BRMS_CONFIRM_WARMUP <- floor(BRMS_CONFIRM_ITER / 2)

BRMS_CONFIRM_PRIOR_SD_B <- 0.3
BRMS_CONFIRM_PRIOR_SD_INT <- 3


## Paralelización
ncores_avail <- parallel::detectCores()

if (is.na(ncores_avail)) {
  ncores_avail <- 2
}

N_CORES <- max(1, ncores_avail - 2)

cat("N_CORES:", N_CORES, "\n")
cat("============================================================\n")

## ============================================================
## 1) HELPERS GENERALES
## ============================================================

find_existing_file <- function(fname) {
  
  candidates <- unique(c(
    file.path(INPUT_RDS_DIR, fname),
    file.path(INPUT_TABLES_DIR, fname),
    file.path(INPUT_MODELAMIENTO_DIR, fname),
    file.path(ANALYSIS_DIR, fname),
    fname
  ))
  
  hit <- candidates[file.exists(candidates)]
  
  if (length(hit) == 0) {
    stop(
      "No encuentro el archivo: ", fname, "\n",
      "Busqué en:\n",
      paste(candidates, collapse = "\n")
    )
  }
  
  hit[1]
}

sanitize_names <- function(v) {
  v2 <- gsub("[^A-Za-z0-9]+", "_", v)
  v2 <- gsub("__+", "_", v2)
  v2 <- gsub("^_+|_+$", "", v2)
  v2 <- ifelse(grepl("^[0-9]", v2), paste0("X", v2), v2)
  v2 <- make.names(v2, unique = FALSE)
  v2 <- gsub("\\.+", "_", v2)
  v2 <- gsub("^_+|_+$", "", v2)
  v2 <- gsub("__+", "_", v2)
  make.unique(v2, sep = "_")
}

.print_tbl <- function(x, n = Inf) print(tibble::as_tibble(x), n = n)

safe_balacc <- function(pred, obs) {
  obs <- factor(obs)
  pred <- factor(pred, levels = levels(obs))
  out <- sapply(levels(obs), function(cl) {
    idx <- obs == cl
    if (sum(idx) == 0) return(NA_real_)
    mean(pred[idx] == cl)
  })
  mean(out, na.rm = TRUE)
}

safe_logloss <- function(prob_mat, y_true) {
  y_true <- factor(y_true)
  lev <- levels(y_true)
  prob_mat <- as.matrix(prob_mat)
  prob_mat <- prob_mat[, lev, drop = FALSE]
  Y <- model.matrix(~ y_true - 1)
  colnames(Y) <- lev
  eps <- 1e-15
  -mean(rowSums(Y * log(pmax(pmin(prob_mat, 1 - eps), eps))))
}

safe_multiclass_auc <- function(prob_mat, y_true) {
  y_true <- factor(y_true)
  lev <- levels(y_true)
  prob_mat <- as.matrix(prob_mat)
  prob_mat <- prob_mat[, lev, drop = FALSE]
  tryCatch(
    as.numeric(pROC::multiclass.roc(y_true, prob_mat)$auc),
    error = function(e) NA_real_
  )
}

class_metrics_from_class <- function(pred, obs) {
  obs <- factor(obs)
  pred <- factor(pred, levels = levels(obs))
  lev <- levels(obs)
  
  dplyr::bind_rows(lapply(lev, function(cl) {
    TP <- sum(pred == cl & obs == cl, na.rm = TRUE)
    FN <- sum(pred != cl & obs == cl, na.rm = TRUE)
    FP <- sum(pred == cl & obs != cl, na.rm = TRUE)
    TN <- sum(pred != cl & obs != cl, na.rm = TRUE)
    
    tibble::tibble(
      Class = cl,
      TP = TP,
      FN = FN,
      FP = FP,
      TN = TN,
      Sensitivity = ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_),
      Specificity = ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_),
      Precision   = ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_),
      F1 = ifelse(
        (2 * TP + FP + FN) > 0,
        2 * TP / (2 * TP + FP + FN),
        NA_real_
      )
    )
  }))
}

macro_metrics_from_class <- function(pred, obs) {
  cm <- class_metrics_from_class(pred, obs)
  
  tibble::tibble(
    MacroSensitivity = mean(cm$Sensitivity, na.rm = TRUE),
    MacroSpecificity = mean(cm$Specificity, na.rm = TRUE),
    MacroPrecision   = mean(cm$Precision, na.rm = TRUE),
    MacroF1          = mean(cm$F1, na.rm = TRUE)
  )
}

metrics_from_class <- function(pred, obs) {
  obs <- factor(obs)
  pred <- factor(pred, levels = levels(obs))
  bal <- safe_balacc(pred, obs)
  
  tibble::tibble(
    Accuracy = mean(pred == obs),
    BalAcc = bal,
    BER = 1 - bal
  ) %>%
    dplyr::bind_cols(
      macro_metrics_from_class(pred, obs)
    )
}

metrics_from_probs <- function(prob_mat, obs) {
  obs <- factor(obs)
  lev <- levels(obs)
  prob_mat <- as.matrix(prob_mat)
  prob_mat <- prob_mat[, lev, drop = FALSE]
  
  pred <- factor(
    colnames(prob_mat)[max.col(prob_mat, ties.method = "first")],
    levels = lev
  )
  
  metrics_from_class(pred, obs) %>%
    dplyr::mutate(
      AUC = safe_multiclass_auc(prob_mat, obs),
      LogLoss = safe_logloss(prob_mat, obs)
    )
}

safe_mixomics_auc <- function(fit, X_new, y_new, ncomp_use) {
  
  auc_obj <- tryCatch(
    mixOmics::auroc(
      object = fit,
      newdata = as.matrix(X_new),
      outcome.test = factor(y_new),
      roc.comp = ncomp_use,
      plot = FALSE,
      print = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(auc_obj)) return(NA_real_)
  
  comp_name <- paste0("Comp", ncomp_use)
  comp_auc <- auc_obj[[comp_name]]
  
  if (is.null(comp_auc)) {
    comp_name2 <- grep(
      paste0("^Comp", ncomp_use, "$|^comp", ncomp_use, "$"),
      names(auc_obj),
      value = TRUE
    )[1]
    comp_auc <- auc_obj[[comp_name2]]
  }
  
  if (is.null(comp_auc)) return(NA_real_)
  
  if (is.data.frame(comp_auc) || is.matrix(comp_auc)) {
    auc_col <- grep("auc", colnames(comp_auc), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(auc_col)) {
      vals <- suppressWarnings(as.numeric(comp_auc[, auc_col]))
      vals <- vals[is.finite(vals)]
      if (length(vals) > 0) return(mean(vals, na.rm = TRUE))
    }
  }
  
  flat <- unlist(comp_auc, recursive = TRUE, use.names = TRUE)
  nm <- names(flat)
  
  if (!is.null(nm)) {
    vals <- suppressWarnings(as.numeric(flat[grepl("auc", nm, ignore.case = TRUE)]))
  } else {
    vals <- suppressWarnings(as.numeric(flat))
  }
  
  vals <- vals[is.finite(vals) & vals >= 0 & vals <= 1]
  
  if (length(vals) == 0) return(NA_real_)
  
  mean(vals, na.rm = TRUE)
}


## ============================================================
## Helpers BRMS confirmatorio
## ============================================================

.detect_backend_safe <- function() {
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    v <- try(cmdstanr::cmdstan_version(), silent = TRUE)
    if (!inherits(v, "try-error")) return("cmdstanr")
  }
  "rstan"
}

.make_ctrl_safe <- function(backend) {
  if (identical(backend, "cmdstanr")) {
    list(adapt_delta = 0.95, max_treedepth = 12)
  } else {
    list(adapt_delta = 0.95, max_treedepth = 12, init_r = 0.1)
  }
}

make_multinom_priors <- function(y_levels,
                                 ref = REF_CLASS,
                                 sd_b = 0.3,
                                 sd_int = 3) {
  
  levs <- setdiff(y_levels, ref)
  
  do.call(c, lapply(levs, function(lev) {
    
    dpar <- paste0("mu", lev)
    
    c(
      brms::set_prior(
        paste0("student_t(3,0,", sd_b, ")"),
        class = "b",
        dpar = dpar
      ),
      brms::set_prior(
        paste0("student_t(3,0,", sd_int, ")"),
        class = "Intercept",
        dpar = dpar
      )
    )
  }))
}

brms_probs <- function(fit, newdata) {
  
  pp <- brms::posterior_epred(fit, newdata = newdata)
  M <- apply(pp, c(2, 3), mean)
  
  cls <- dimnames(pp)[[3]]
  
  if (is.null(cls)) {
    cls <- fit$family$names
  }
  
  colnames(M) <- cls
  M
}

extract_brms_confirmatory_posteriors <- function(fit, panel_tbl) {
  
  fx <- as.data.frame(brms::fixef(fit, robust = TRUE))
  fx$term <- rownames(fx)
  
  fx_tbl <- fx %>%
    tibble::as_tibble() %>%
    dplyr::filter(!grepl("Intercept", term)) %>%
    dplyr::mutate(
      class = sub("^mu([^_]+)_.*$", "\\1", term),
      feature_model = sub("^mu[^_]+_", "", term),
      abs_estimate = abs(Estimate),
      ci_excludes_zero = (`Q2.5` > 0) | (`Q97.5` < 0),
      posterior_direction_ci = dplyr::case_when(
        `Q2.5` > 0 ~ "positive",
        `Q97.5` < 0 ~ "negative",
        TRUE ~ "uncertain"
      )
    )
  
  draws_df <- posterior::as_draws_df(fit)
  draws_raw <- as.data.frame(draws_df)
  
  coef_cols <- grep("^b_mu", colnames(draws_raw), value = TRUE)
  
  prob_dir_tbl <- dplyr::bind_rows(lapply(coef_cols, function(cc) {
    
    x <- as.numeric(draws_raw[[cc]])
    
    term <- sub("^b_", "", cc)
    class <- sub("^mu([^_]+)_.*$", "\\1", term)
    feature_model <- sub("^mu[^_]+_", "", term)
    
    tibble::tibble(
      term = term,
      class = class,
      feature_model = feature_model,
      prob_positive = mean(x > 0, na.rm = TRUE),
      prob_negative = mean(x < 0, na.rm = TRUE),
      posterior_direction_prob = dplyr::case_when(
        mean(x > 0, na.rm = TRUE) >= 0.95 ~ "positive",
        mean(x < 0, na.rm = TRUE) >= 0.95 ~ "negative",
        TRUE ~ "uncertain"
      )
    )
  }))
  
  fx_tbl %>%
    dplyr::left_join(
      prob_dir_tbl,
      by = c("term", "class", "feature_model")
    ) %>%
    dplyr::left_join(
      panel_tbl %>%
        dplyr::select(
          dplyr::any_of(c(
            "feature_model",
            "view",
            "view_code",
            "feature_label",
            "rank_mofa_for_selection",
            "score_mofa_for_selection",
            "score_mofa_percentile_view",
            "mofa_pass",
            "splsda_freq",
            "splsda_pass",
            "lasso_freq",
            "lasso_pass",
            "plsda_freq",
            "plsda_pass",
            "n_supervised_pass",
            "weighted_supervised_score",
            "candidate_panel_final",
            "consensus_score"
          ))
        ),
      by = "feature_model"
    ) %>%
    dplyr::arrange(
      dplyr::desc(ci_excludes_zero),
      dplyr::desc(abs_estimate),
      class,
      feature_model
    )
}

select_panel_final_for_brms <- function(panel_final,
                                        max_features = Inf,
                                        min_features = 3,
                                        panel_name = "panel_final") {
  
  panel_use <- panel_final
  
  if (nrow(panel_use) == 0 || nrow(panel_use) < min_features) {
    return(list(
      panel_name = NA_character_,
      panel_tbl = tibble::tibble(),
      features = character(0)
    ))
  }
  
  panel_use <- panel_use %>%
    dplyr::arrange(
      dplyr::desc(n_supervised_pass),
      dplyr::desc(weighted_supervised_score),
      dplyr::desc(score_mofa_percentile_view),
      dplyr::desc(consensus_score),
      rank_mofa_for_selection
    )
  
  if (is.finite(max_features) && nrow(panel_use) > max_features) {
    panel_use <- panel_use %>%
      dplyr::slice_head(n = max_features)
  }
  
  list(
    panel_name = panel_name,
    panel_tbl = panel_use,
    features = panel_use$feature_model
  )
}

make_brms_panel_variants <- function(panel_final,
                                     top_n = NULL,
                                     min_features = 3) {
  
  list(
    all = select_panel_final_for_brms(
      panel_final = panel_final,
      max_features = Inf,
      min_features = min_features,
      panel_name = "panel_final_all"
    )
  )
}

make_train_indices <- function(y, k = 3, repeats = 10, seed = 123) {
  set.seed(seed)
  y <- factor(y)
  min_class <- min(table(y))
  k <- min(k, min_class)
  if (k < 2) stop("No se puede hacer CV: alguna clase tiene menos de 2 muestras.")
  caret::createMultiFolds(y, k = k, times = repeats)
}

make_valid_indices_from_train <- function(train_indices, n) {
  lapply(train_indices, function(idx_train) setdiff(seq_len(n), idx_train))
}

get_feature_view <- function(feature_model) {
  view_code <- sub("_.*$", "", feature_model)
  dplyr::case_when(
    view_code == "tx" ~ "transcriptomica",
    view_code == "pr" ~ "proteomica",
    view_code == "me" ~ "metabolomica",
    view_code == "cl" ~ "clinical",
    TRUE ~ view_code
  )
}

## ============================================================
## 2) CARGAR DATOS FINALES
## ============================================================

df_train_path <- find_existing_file("df_train_final.rds")
df_test_path  <- find_existing_file("df_test_final.rds")

df_train <- readRDS(df_train_path)
df_test  <- readRDS(df_test_path)

df_train$y <- factor(df_train$y)
df_train$y <- stats::relevel(df_train$y, ref = REF_CLASS)
df_test$y <- factor(df_test$y, levels = levels(df_train$y))

NON_BIOMARKER_COLS <- c("y", "Factor1", "Factor2")

common_features <- setdiff(
  intersect(colnames(df_train), colnames(df_test)),
  NON_BIOMARKER_COLS
)

if (any(c("Factor1", "Factor2") %in% colnames(df_train))) {
  cat("\n--- Columnas latentes excluidas del panel final ---\n")
  print(intersect(c("Factor1", "Factor2"), colnames(df_train)))
}

new_names <- make.unique(sanitize_names(common_features), sep = "_")

colnames(df_train)[match(common_features, colnames(df_train))] <- new_names
colnames(df_test)[match(common_features, colnames(df_test))] <- new_names

df_train <- df_train[, c("y", new_names), drop = FALSE]
df_test  <- df_test[,  c("y", new_names), drop = FALSE]

is_num <- vapply(df_train[, -1, drop = FALSE], is.numeric, logical(1))
if (!all(is_num)) {
  cat("\n--- Variables no numéricas eliminadas ---\n")
  print(names(is_num)[!is_num])
}

df_train <- df_train[, c("y", names(is_num)[is_num]), drop = FALSE]
df_test  <- df_test[,  c("y", names(is_num)[is_num]), drop = FALSE]
final_features <- setdiff(colnames(df_train), "y")

cat("\n--- Dimensiones df_train / df_test ---\n")
print(dim(df_train))
print(dim(df_test))
cat("\n--- Distribución train ---\n")
print(table(df_train$y))
cat("\n--- Distribución test ---\n")
print(table(df_test$y))
cat("\n--- Número de features finales de entrada ---\n")
print(length(final_features))

## ============================================================
## 3) CARGAR SCORE MOFA DE FEATURES FINALES
## ============================================================
load_mofa_feature_table <- function() {
  
  rds_candidates <- unique(c(
    file.path(INPUT_RDS_DIR, "features_finales_integradas_mofa_brms.rds"),
    file.path(INPUT_MODELAMIENTO_DIR, "features_finales_integradas_mofa_brms.rds"),
    file.path(ANALYSIS_DIR, "features_finales_integradas_mofa_brms.rds"),
    "features_finales_integradas_mofa_brms.rds"
  ))
  
  csv_candidates <- unique(c(
    file.path(INPUT_TABLES_DIR, "features_finales_integradas_mofa_brms.csv"),
    file.path(INPUT_MODELAMIENTO_DIR, "features_finales_integradas_mofa_brms.csv"),
    file.path(ANALYSIS_DIR, "features_finales_integradas_mofa_brms.csv"),
    "features_finales_integradas_mofa_brms.csv"
  ))
  
  rds_hit <- rds_candidates[file.exists(rds_candidates)]
  csv_hit <- csv_candidates[file.exists(csv_candidates)]
  
  if (length(rds_hit) > 0) {
    cat("\nTabla MOFA/BRMS cargada desde RDS:\n")
    cat(rds_hit[1], "\n")
    return(readRDS(rds_hit[1]))
  }
  
  if (length(csv_hit) > 0) {
    cat("\nTabla MOFA/BRMS cargada desde CSV:\n")
    cat(csv_hit[1], "\n")
    return(read.csv(csv_hit[1], check.names = FALSE))
  }
  
  warning(
    "No encontré features_finales_integradas_mofa_brms.rds/csv. ",
    "Se creará tabla mínima sin score MOFA."
  )
  
  tibble::tibble(
    feature_model = final_features,
    view_code = sub("_.*$", "", final_features),
    view = get_feature_view(final_features),
    feature_label = sub("^[^_]+_", "", final_features),
    rank_mofa_for_selection = NA_integer_,
    score_mofa_for_selection = NA_real_,
    score_mofa_percentile_view = NA_real_
  )
}

mofa_tbl <- load_mofa_feature_table() %>% tibble::as_tibble()
if (!"feature_model" %in% colnames(mofa_tbl)) stop("La tabla MOFA no tiene columna feature_model.")

mofa_tbl <- mofa_tbl %>%
  dplyr::filter(feature_model %in% final_features) %>%
  dplyr::distinct(feature_model, .keep_all = TRUE)

missing_mofa <- setdiff(final_features, mofa_tbl$feature_model)
if (length(missing_mofa) > 0) {
  mofa_tbl <- dplyr::bind_rows(
    mofa_tbl,
    tibble::tibble(
      feature_model = missing_mofa,
      view_code = sub("_.*$", "", missing_mofa),
      view = get_feature_view(missing_mofa),
      feature_label = sub("^[^_]+_", "", missing_mofa),
      rank_mofa_for_selection = NA_integer_,
      score_mofa_for_selection = NA_real_,
      score_mofa_percentile_view = NA_real_
    )
  )
}

if (!"view_code" %in% colnames(mofa_tbl)) mofa_tbl$view_code <- sub("_.*$", "", mofa_tbl$feature_model)
if (!"view" %in% colnames(mofa_tbl)) mofa_tbl$view <- get_feature_view(mofa_tbl$feature_model)
if (!"feature_label" %in% colnames(mofa_tbl)) mofa_tbl$feature_label <- sub("^[^_]+_", "", mofa_tbl$feature_model)
if (!"score_mofa_percentile_view" %in% colnames(mofa_tbl)) mofa_tbl$score_mofa_percentile_view <- NA_real_
if (!"score_mofa_for_selection" %in% colnames(mofa_tbl)) mofa_tbl$score_mofa_for_selection <- NA_real_

mofa_pass_tbl <- mofa_tbl %>%
  dplyr::mutate(
    mofa_pass = dplyr::case_when(
      is.na(score_mofa_percentile_view) ~ TRUE,
      score_mofa_percentile_view >= MOFA_PERCENTILE_THR ~ TRUE,
      TRUE ~ FALSE
    )
  )

cat("\n--- Auditoría tabla MOFA cargada ---\n")
print(dim(mofa_pass_tbl))
cat("\nFeatures finales sin score MOFA:\n")
print(setdiff(final_features, mofa_pass_tbl$feature_model))

## ============================================================
## 4) FOLDS CV COMUNES
## ============================================================

train_indices <- make_train_indices(df_train$y, k = K_CV, repeats = N_REPEATS_CV, seed = 123)
valid_indices <- make_valid_indices_from_train(train_indices, nrow(df_train))

cat("\n--- Número de folds/repeticiones CV ---\n")
print(length(train_indices))

fold_audit <- dplyr::bind_rows(lapply(seq_along(valid_indices), function(i) {
  tb <- table(df_train$y[valid_indices[[i]]])
  tibble::tibble(FoldID = names(valid_indices)[i], Class = names(tb), n = as.integer(tb))
}))
cat("\n--- Auditoría folds de validación ---\n")
.print_tbl(fold_audit, n = Inf)

## ============================================================
## 5) MODELO 1: sPLSDA
## ============================================================

extract_splsda_features <- function(model) {
  L <- model$loadings$X
  unique(unlist(lapply(seq_len(ncol(L)), function(j) rownames(L)[abs(L[, j]) > 0])))
}

run_splsda_one_seed <- function(seed, df_train, ncomp = 2, keepX_grid = seq(5, 70, 5), folds = 3, nrepeat = 5) {
  set.seed(seed)
  X <- df_train[, -1, drop = FALSE]
  Y <- df_train$y
  keepX_grid <- keepX_grid[keepX_grid <= ncol(X)]
  keepX_grid <- keepX_grid[keepX_grid >= 1]
  if (length(keepX_grid) == 0) keepX_grid <- seq_len(min(10, ncol(X)))
  tune_obj <- tryCatch(
    mixOmics::tune.splsda(
      X = X,
      Y = Y,
      ncomp = ncomp,
      validation = "Mfold",
      folds = min(folds, min(table(Y))),
      nrepeat = nrepeat,
      test.keepX = keepX_grid,
      progressBar = FALSE,
      dist = "max.dist",
      measure = "BER"
    ),
    error = function(e) e
  )
  if (inherits(tune_obj, "error")) {
    return(list(seed = seed, ok = FALSE, error = conditionMessage(tune_obj), keepX = NA_character_, features = character(0)))
  }
  keepX <- as.integer(tune_obj$choice.keepX)
  mdl <- tryCatch(
    mixOmics::splsda(X = X, Y = Y, ncomp = ncomp, keepX = keepX, scale = TRUE),
    error = function(e) e
  )
  if (inherits(mdl, "error")) {
    return(list(seed = seed, ok = FALSE, error = conditionMessage(mdl), keepX = paste(keepX, collapse = ";"), features = character(0)))
  }
  list(seed = seed, ok = TRUE, error = NA_character_, keepX = paste(keepX, collapse = ";"), features = extract_splsda_features(mdl), model = mdl)
}

cat("\n============================================================\n")
cat("MODELO 1: sPLSDA estabilidad por semillas\n")
cat("============================================================\n")

splsda_seed_results <- lapply(SEEDS_SPLSDA, function(s) {
  run_splsda_one_seed(s, df_train, NCOMP_SPLSDA, KEEPX_GRID, K_CV, nrepeat = 5)
})

splsda_seed_audit <- tibble::tibble(
  seed = SEEDS_SPLSDA,
  ok = vapply(splsda_seed_results, `[[`, logical(1), "ok"),
  error = vapply(splsda_seed_results, function(x) if (is.null(x$error)) NA_character_ else x$error, character(1)),
  keepX = vapply(splsda_seed_results, function(x) if (is.null(x$keepX)) NA_character_ else as.character(x$keepX), character(1)),
  n_features = vapply(splsda_seed_results, function(x) length(x$features), integer(1))
)
.print_tbl(splsda_seed_audit, n = Inf)

splsda_features_long <- dplyr::bind_rows(lapply(splsda_seed_results, function(x) {
  if (!isTRUE(x$ok) || length(x$features) == 0) return(NULL)
  tibble::tibble(seed = x$seed, feature_model = x$features)
}))

n_splsda_ok <- sum(splsda_seed_audit$ok)
if (n_splsda_ok == 0) {
  warning("sPLSDA falló en todas las semillas.")
  splsda_feature_tbl <- tibble::tibble(feature_model = final_features, splsda_n_selected = 0L, splsda_freq = 0, splsda_pass = FALSE)
} else {
  splsda_feature_tbl <- splsda_features_long %>%
    dplyr::count(feature_model, name = "splsda_n_selected") %>%
    dplyr::mutate(splsda_freq = splsda_n_selected / n_splsda_ok, splsda_pass = splsda_freq >= FREQ_THR_MODEL)
  splsda_feature_tbl <- tibble::tibble(feature_model = final_features) %>%
    dplyr::left_join(splsda_feature_tbl, by = "feature_model") %>%
    dplyr::mutate(
      splsda_n_selected = dplyr::coalesce(splsda_n_selected, 0L),
      splsda_freq = dplyr::coalesce(splsda_freq, 0),
      splsda_pass = dplyr::coalesce(splsda_pass, FALSE)
    )
}

cat("\n--- Top sPLSDA por frecuencia ---\n")
.print_tbl(splsda_feature_tbl %>% dplyr::arrange(dplyr::desc(splsda_freq)), n = 30)

## CV sPLSDA con keepX mediano
keepX_ok <- splsda_seed_audit$keepX[splsda_seed_audit$ok]
if (length(keepX_ok) > 0) {
  keepX_mat <- do.call(rbind, strsplit(keepX_ok, ";"))
  keepX_mat <- apply(keepX_mat, 2, as.integer)
  keepX_final_vec <- as.integer(round(apply(keepX_mat, 2, median, na.rm = TRUE)))
} else {
  keepX_final_vec <- rep(min(10, ncol(df_train) - 1), NCOMP_SPLSDA)
}

splsda_cv_pred <- dplyr::bind_rows(lapply(seq_along(train_indices), function(i) {
  idx_tr <- train_indices[[i]]
  idx_va <- valid_indices[[i]]
  dat_tr <- df_train[idx_tr, , drop = FALSE]
  dat_va <- df_train[idx_va, , drop = FALSE]
  fit <- tryCatch(
    mixOmics::splsda(X = dat_tr[, -1, drop = FALSE], Y = dat_tr$y, ncomp = NCOMP_SPLSDA, keepX = keepX_final_vec, scale = TRUE),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(NULL)
  pred <- stats::predict(fit, as.matrix(dat_va[, -1, drop = FALSE]))
  y_pred <- pred$MajorityVote$max.dist[, paste0("comp", NCOMP_SPLSDA)]
  
  auc_fold <- safe_mixomics_auc(
    fit = fit,
    X_new = dat_va[, -1, drop = FALSE],
    y_new = dat_va$y,
    ncomp_use = NCOMP_SPLSDA
  )
  
  tibble::tibble(
    FoldID = names(train_indices)[i],
    id = rownames(dat_va),
    obs = dat_va$y,
    pred = factor(y_pred, levels = levels(df_train$y)),
    AUC = auc_fold
  )
  
  }))

splsda_cv_metrics <- splsda_cv_pred %>%
  dplyr::group_by(FoldID) %>%
  dplyr::group_modify(~ {
    metrics_from_class(
      pred = factor(.x$pred, levels = levels(df_train$y)),
      obs  = factor(.x$obs,  levels = levels(df_train$y))
    ) %>%
      dplyr::mutate(
        AUC = dplyr::first(.x$AUC),
        LogLoss = NA_real_
      )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Model = "sPLSDA") %>%
  dplyr::select(
    Model,
    FoldID,
    Accuracy,
    BalAcc,
    BER,
    MacroSensitivity,
    MacroSpecificity,
    MacroPrecision,
    MacroF1,
    AUC,
    LogLoss
  )
## ============================================================
## 6) MODELO 2 Y 3: GLMNET MULTINOMIAL
## ============================================================

extract_glmnet_selected <- function(cfit, lambda, top_prop = NULL, nonzero = TRUE) {
  coefs <- coef(cfit$finalModel, s = lambda)
  if (!is.list(coefs)) coefs <- list(coefs)
  all_tab <- dplyr::bind_rows(lapply(names(coefs), function(cl) {
    m <- as.matrix(coefs[[cl]])
    tibble::tibble(class = cl, feature_model = rownames(m), coef = as.numeric(m[, 1]))
  })) %>%
    dplyr::filter(feature_model != "(Intercept)") %>%
    dplyr::group_by(feature_model) %>%
    dplyr::summarise(
      coef_max_abs = max(abs(coef), na.rm = TRUE),
      coef_mean_abs = mean(abs(coef), na.rm = TRUE),
      coef_classes_nonzero = sum(abs(coef) > 0, na.rm = TRUE),
      .groups = "drop"
    )
  if (nonzero) {
    all_tab <- all_tab %>% dplyr::mutate(selected = coef_max_abs > 0)
  } else {
    n_top <- max(1L, ceiling(nrow(all_tab) * top_prop))
    top_features <- all_tab %>%
      dplyr::arrange(dplyr::desc(coef_max_abs)) %>%
      dplyr::slice_head(n = n_top) %>%
      dplyr::pull(feature_model)
    all_tab <- all_tab %>% dplyr::mutate(selected = feature_model %in% top_features)
  }
  all_tab
}

run_glmnet_seed <- function(seed, df_train, alpha_grid, lambda_grid, train_indices, metric = "logLoss", maximize = FALSE, weights_balanced = TRUE) {
  set.seed(seed)
  y <- df_train$y
  weights <- if (weights_balanced) as.numeric(1 / table(y)[y]) else rep(1, length(y))
  ctrl <- caret::trainControl(
    method = "repeatedcv",
    number = K_CV,
    repeats = N_REPEATS_CV,
    index = train_indices,
    classProbs = TRUE,
    summaryFunction = caret::multiClassSummary,
    savePredictions = "final",
    allowParallel = TRUE
  )
  tryCatch(
    caret::train(
      y ~ .,
      data = df_train,
      method = "glmnet",
      family = "multinomial",
      trControl = ctrl,
      tuneGrid = expand.grid(alpha = alpha_grid, lambda = lambda_grid),
      metric = metric,
      maximize = maximize,
      weights = weights
    ),
    error = function(e) e
  )
}

summarise_glmnet_feature_stability <- function(results, final_features, prefix, freq_thr) {
  audit <- tibble::tibble(
    seed = vapply(results, `[[`, numeric(1), "seed"),
    ok = vapply(results, `[[`, logical(1), "ok"),
    error = vapply(results, function(x) if (is.null(x$error)) NA_character_ else x$error, character(1))
  )
  n_ok <- sum(audit$ok)
  if (n_ok == 0) {
    out <- tibble::tibble(feature_model = final_features)
    out[[paste0(prefix, "_n_selected")]] <- 0L
    out[[paste0(prefix, "_freq")]] <- 0
    out[[paste0(prefix, "_coef_max_abs_median")]] <- NA_real_
    out[[paste0(prefix, "_coef_max_abs_mean")]] <- NA_real_
    out[[paste0(prefix, "_pass")]] <- FALSE
    return(list(audit = audit, features = out))
  }
  long <- dplyr::bind_rows(lapply(results, function(x) {
    if (!isTRUE(x$ok)) return(NULL)
    x$selected %>% dplyr::filter(selected) %>% dplyr::mutate(seed = x$seed)
  }))
  out <- long %>%
    dplyr::group_by(feature_model) %>%
    dplyr::summarise(
      n_selected = dplyr::n_distinct(seed),
      coef_max_abs_median = median(coef_max_abs, na.rm = TRUE),
      coef_max_abs_mean = mean(coef_max_abs, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(freq = n_selected / n_ok, pass = freq >= freq_thr)
  out <- tibble::tibble(feature_model = final_features) %>%
    dplyr::left_join(out, by = "feature_model") %>%
    dplyr::mutate(
      n_selected = dplyr::coalesce(n_selected, 0L),
      freq = dplyr::coalesce(freq, 0),
      pass = dplyr::coalesce(pass, FALSE)
    )
  names(out) <- gsub("^n_selected$", paste0(prefix, "_n_selected"), names(out))
  names(out) <- gsub("^freq$", paste0(prefix, "_freq"), names(out))
  names(out) <- gsub("^coef_max_abs_median$", paste0(prefix, "_coef_max_abs_median"), names(out))
  names(out) <- gsub("^coef_max_abs_mean$", paste0(prefix, "_coef_max_abs_mean"), names(out))
  names(out) <- gsub("^pass$", paste0(prefix, "_pass"), names(out))
  list(audit = audit, features = out)
}

summarise_glmnet_cv_metrics <- function(results, model_name, df_train) {
  
  ok_results <- results[vapply(results, function(x) isTRUE(x$ok), logical(1))]
  
  if (length(ok_results) == 0) {
    return(tibble::tibble())
  }
  
  ## Para métricas CV usa un solo ajuste, porque con folds fijos las semillas duplican lo mismo.
  x <- ok_results[[1]]
  
  pred <- x$fit$pred
  bt <- x$fit$bestTune
  
  pred <- pred %>%
    dplyr::filter(alpha == bt$alpha, lambda == bt$lambda)
  
  prob_cols <- intersect(levels(df_train$y), colnames(pred))
  
  pred %>%
    dplyr::group_by(Resample) %>%
    dplyr::group_modify(~ {
      metrics_from_probs(
        prob_mat = as.matrix(.x[, prob_cols, drop = FALSE]),
        obs = factor(.x$obs, levels = levels(df_train$y))
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Model = model_name,
      seed = x$seed,
      alpha = bt$alpha,
      lambda = bt$lambda,
      FoldID = Resample
    ) %>%
    dplyr::select(
      Model,
      FoldID,
      Accuracy,
      BalAcc,
      BER,
      MacroSensitivity,
      MacroSpecificity,
      MacroPrecision,
      MacroF1,
      AUC,
      LogLoss,
      seed,
      alpha,
      lambda
    )
  }



run_multinom_cv <- function(df_train, train_indices, valid_indices) {
  
  dplyr::bind_rows(lapply(seq_along(train_indices), function(i) {
    
    idx_tr <- train_indices[[i]]
    idx_va <- valid_indices[[i]]
    
    dat_tr <- df_train[idx_tr, , drop = FALSE]
    dat_va <- df_train[idx_va, , drop = FALSE]
    
    fit <- tryCatch(
      nnet::multinom(
        y ~ .,
        data = dat_tr,
        trace = FALSE,
        MaxNWts = 10000,
        maxit = 1000
      ),
      error = function(e) e
    )
    
    if (inherits(fit, "error")) return(NULL)
    
    prob <- tryCatch(
      predict(fit, newdata = dat_va, type = "probs"),
      error = function(e) NULL
    )
    
    if (is.null(prob)) return(NULL)
    
    prob <- as.matrix(prob)
    
    if (is.null(colnames(prob))) {
      colnames(prob) <- levels(df_train$y)
    }
    
    prob <- prob[, levels(df_train$y), drop = FALSE]
    
    pred <- factor(
      colnames(prob)[max.col(prob, ties.method = "first")],
      levels = levels(df_train$y)
    )
    
    prob_df <- as.data.frame(prob, check.names = FALSE)
    colnames(prob_df) <- paste0("prob_", colnames(prob_df))
    
    dplyr::bind_cols(
      tibble::tibble(
        FoldID = names(train_indices)[i],
        id = rownames(dat_va),
        obs = dat_va$y,
        pred = pred
      ),
      prob_df
    )
  }))
}


cat("\n============================================================\n")
cat("MODELO 2: LASSO multinomial\n")
cat("============================================================\n")

cl <- parallel::makeCluster(N_CORES)
doParallel::registerDoParallel(cl)
lasso_results <- lapply(SEEDS_GLMNET, function(s) {
  fit <- run_glmnet_seed(s, df_train, ALPHA_GRID_LASSO, LAMBDA_GRID, train_indices)
  if (inherits(fit, "error")) return(list(seed = s, ok = FALSE, error = conditionMessage(fit), fit = NULL))
  sel <- extract_glmnet_selected(fit, fit$bestTune$lambda, nonzero = TRUE)
  list(seed = s, ok = TRUE, error = NA_character_, fit = fit, selected = sel)
})
parallel::stopCluster(cl)
foreach::registerDoSEQ()

lasso_sum <- summarise_glmnet_feature_stability(lasso_results, final_features, "lasso", FREQ_THR_MODEL)
lasso_audit <- lasso_sum$audit
lasso_feature_tbl <- lasso_sum$features
.print_tbl(lasso_audit, n = Inf)

cat("\n--- Top LASSO por frecuencia ---\n")
.print_tbl(
  lasso_feature_tbl %>%
    dplyr::arrange(dplyr::desc(lasso_freq), dplyr::desc(lasso_coef_max_abs_mean)),
  n = 30
)

lasso_cv_metrics <- summarise_glmnet_cv_metrics(lasso_results, "LASSO", df_train)





## ============================================================
## 7) MODELO 3: PLSDA
## ============================================================

extract_plsda_top_features <- function(model, top_prop = 0.35) {
  L <- model$loadings$X[, seq_len(min(NCOMP_PLSDA, ncol(model$loadings$X))), drop = FALSE]
  imp <- apply(abs(L), 1, max, na.rm = TRUE)
  n_top <- max(1L, ceiling(length(imp) * top_prop))
  tibble::tibble(feature_model = names(imp), plsda_importance = as.numeric(imp)) %>%
    dplyr::arrange(dplyr::desc(plsda_importance)) %>%
    dplyr::mutate(plsda_rank = dplyr::row_number(), selected = plsda_rank <= n_top)
}

cat("\n============================================================\n")
cat("MODELO 3: PLSDA\n")
cat("============================================================\n")

plsda_fold_results <- lapply(seq_along(train_indices), function(i) {
  idx_tr <- train_indices[[i]]
  idx_va <- valid_indices[[i]]
  dat_tr <- df_train[idx_tr, , drop = FALSE]
  dat_va <- df_train[idx_va, , drop = FALSE]
  fit <- tryCatch(
    mixOmics::plsda(X = dat_tr[, -1, drop = FALSE], Y = dat_tr$y, ncomp = NCOMP_PLSDA, scale = TRUE),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(list(FoldID = names(train_indices)[i], ok = FALSE, error = conditionMessage(fit), selected = NULL, pred = NULL))
  selected <- extract_plsda_top_features(fit, top_prop = TOP_PROP_DENSE_MODELS)
  pred <- stats::predict(fit, as.matrix(dat_va[, -1, drop = FALSE]))
  y_pred <- pred$MajorityVote$max.dist[, paste0("comp", NCOMP_PLSDA)]
  auc_fold <- safe_mixomics_auc(
    fit = fit,
    X_new = dat_va[, -1, drop = FALSE],
    y_new = dat_va$y,
    ncomp_use = NCOMP_PLSDA
  )
  
  pred_tbl <- tibble::tibble(
    FoldID = names(train_indices)[i],
    id = rownames(dat_va),
    obs = dat_va$y,
    pred = factor(y_pred, levels = levels(df_train$y)),
    AUC = auc_fold
  )
  list(FoldID = names(train_indices)[i], ok = TRUE, error = NA_character_, selected = selected, pred = pred_tbl)
})

plsda_audit <- tibble::tibble(
  FoldID = names(train_indices),
  ok = vapply(plsda_fold_results, `[[`, logical(1), "ok"),
  error = vapply(plsda_fold_results, function(x) if (is.null(x$error)) NA_character_ else x$error, character(1))
)
.print_tbl(plsda_audit, n = Inf)

plsda_selected_long <- dplyr::bind_rows(lapply(plsda_fold_results, function(x) {
  if (!isTRUE(x$ok)) return(NULL)
  x$selected %>% dplyr::filter(selected) %>% dplyr::mutate(FoldID = x$FoldID)
}))
n_plsda_ok <- sum(plsda_audit$ok)

if (n_plsda_ok == 0) {
  plsda_feature_tbl <- tibble::tibble(
    feature_model = final_features,
    plsda_n_selected = 0L,
    plsda_importance_median = NA_real_,
    plsda_importance_mean = NA_real_,
    plsda_freq = 0,
    plsda_pass = FALSE
  )
} else {
  plsda_feature_tbl <- plsda_selected_long %>%
    dplyr::group_by(feature_model) %>%
    dplyr::summarise(
      plsda_n_selected = dplyr::n_distinct(FoldID),
      plsda_importance_median = median(plsda_importance, na.rm = TRUE),
      plsda_importance_mean = mean(plsda_importance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(plsda_freq = plsda_n_selected / n_plsda_ok, plsda_pass = plsda_freq >= FREQ_THR_MODEL)
  plsda_feature_tbl <- tibble::tibble(feature_model = final_features) %>%
    dplyr::left_join(plsda_feature_tbl, by = "feature_model") %>%
    dplyr::mutate(
      plsda_n_selected = dplyr::coalesce(plsda_n_selected, 0L),
      plsda_freq = dplyr::coalesce(plsda_freq, 0),
      plsda_pass = dplyr::coalesce(plsda_pass, FALSE)
    )
}

cat("\n--- Top PLSDA por frecuencia/importancia ---\n")
.print_tbl(plsda_feature_tbl %>% dplyr::arrange(dplyr::desc(plsda_freq), dplyr::desc(plsda_importance_mean)), n = 30)

plsda_cv_pred <- dplyr::bind_rows(lapply(plsda_fold_results, function(x) {
  if (!isTRUE(x$ok)) return(NULL)
  x$pred
}))
plsda_cv_metrics <- plsda_cv_pred %>%
  dplyr::group_by(FoldID) %>%
  dplyr::group_modify(~ {
    metrics_from_class(
      pred = factor(.x$pred, levels = levels(df_train$y)),
      obs  = factor(.x$obs,  levels = levels(df_train$y))
    ) %>%
      dplyr::mutate(
        AUC = dplyr::first(.x$AUC),
        LogLoss = NA_real_
      )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Model = "PLSDA") %>%
  dplyr::select(
    Model,
    FoldID,
    Accuracy,
    BalAcc,
    BER,
    MacroSensitivity,
    MacroSpecificity,
    MacroPrecision,
    MacroF1,
    AUC,
    LogLoss
  )

cat("\n============================================================\n")
cat("MODELO BASE: Multinomial clásico\n")
cat("============================================================\n")

multinom_cv_pred <- run_multinom_cv(
  df_train = df_train,
  train_indices = train_indices,
  valid_indices = valid_indices
)

multinom_cv_metrics <- multinom_cv_pred %>%
  dplyr::group_by(FoldID) %>%
  dplyr::group_modify(~ {
    
    prob_cols <- grep("^prob_", colnames(.x), value = TRUE)
    prob_mat <- as.matrix(.x[, prob_cols, drop = FALSE])
    colnames(prob_mat) <- sub("^prob_", "", colnames(prob_mat))
    
    metrics_from_probs(
      prob_mat = prob_mat,
      obs = factor(.x$obs, levels = levels(df_train$y))
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Model = "Multinom_classic") %>%
  dplyr::select(Model, FoldID, Accuracy, BalAcc, BER,
                MacroSensitivity, MacroSpecificity, MacroPrecision, MacroF1,
                AUC, LogLoss)
## ============================================================
## 8) CONSOLIDAR MÉTRICAS CV
## ============================================================

cv_metrics_all <- dplyr::bind_rows(
  splsda_cv_metrics,
  lasso_cv_metrics %>%
    dplyr::select(
      dplyr::any_of(c(
        "Model",
        "FoldID",
        "Accuracy",
        "BalAcc",
        "BER",
        "MacroSensitivity",
        "MacroSpecificity",
        "MacroPrecision",
        "MacroF1",
        "AUC",
        "LogLoss",
        "seed"
      ))
    ),
  plsda_cv_metrics,
  multinom_cv_metrics
)

cv_metrics_summary <- cv_metrics_all %>%
  dplyr::group_by(Model) %>%
  dplyr::summarise(
    n_eval = dplyr::n(),
    mean_Accuracy = mean(Accuracy, na.rm = TRUE),
    sd_Accuracy = sd(Accuracy, na.rm = TRUE),
    mean_BalAcc = mean(BalAcc, na.rm = TRUE),
    sd_BalAcc = sd(BalAcc, na.rm = TRUE),
    mean_BER = mean(BER, na.rm = TRUE),
    sd_BER = sd(BER, na.rm = TRUE),
    mean_MacroSensitivity = mean(MacroSensitivity, na.rm = TRUE),
    sd_MacroSensitivity = sd(MacroSensitivity, na.rm = TRUE),
    mean_MacroSpecificity = mean(MacroSpecificity, na.rm = TRUE),
    sd_MacroSpecificity = sd(MacroSpecificity, na.rm = TRUE),
    mean_MacroPrecision = mean(MacroPrecision, na.rm = TRUE),
    sd_MacroPrecision = sd(MacroPrecision, na.rm = TRUE),
    mean_MacroF1 = mean(MacroF1, na.rm = TRUE),
    sd_MacroF1 = sd(MacroF1, na.rm = TRUE),
    mean_AUC = mean(AUC, na.rm = TRUE),
    sd_AUC = sd(AUC, na.rm = TRUE),
    mean_LogLoss = mean(LogLoss, na.rm = TRUE),
    sd_LogLoss = sd(LogLoss, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(mean_BalAcc), mean_LogLoss)

cat("\n============================================================\n")
cat("RESUMEN CV POR MODELO\n")
cat("============================================================\n")
.print_tbl(cv_metrics_summary, n = Inf)

## ============================================================
## 9) PANEL FINAL: CONSENSO DE FEATURES
## ============================================================

panel_feature_tbl <- tibble::tibble(feature_model = final_features) %>%
  dplyr::left_join(mofa_pass_tbl, by = "feature_model") %>%
  dplyr::left_join(splsda_feature_tbl, by = "feature_model") %>%
  dplyr::left_join(lasso_feature_tbl, by = "feature_model") %>%
  dplyr::left_join(plsda_feature_tbl, by = "feature_model") %>%
  dplyr::mutate(
    view = dplyr::coalesce(view, get_feature_view(feature_model)),
    view_code = dplyr::coalesce(view_code, sub("_.*$", "", feature_model)),
    feature_label = dplyr::coalesce(feature_label, sub("^[^_]+_", "", feature_model)),
    mofa_pass = dplyr::coalesce(mofa_pass, TRUE),
    splsda_pass = dplyr::coalesce(splsda_pass, FALSE),
    lasso_pass = dplyr::coalesce(lasso_pass, FALSE),
    plsda_pass = dplyr::coalesce(plsda_pass, FALSE),
    
    n_supervised_pass =
      as.integer(splsda_pass) +
      as.integer(lasso_pass) +
      as.integer(plsda_pass),
    
    n_core_pass =
      as.integer(splsda_pass) +
      as.integer(plsda_pass),
    
    weighted_supervised_score =
      W_SPLSDA * as.numeric(splsda_pass) +
      W_LASSO  * as.numeric(lasso_pass) +
      W_PLSDA  * as.numeric(plsda_pass),
    
    ## Panel final:
    ## - MOFA debe pasar
    ## - al menos 2/3 modelos supervisados
    ## - al menos uno de los modelos estructurales debe apoyar: sPLSDA o PLSDA
    candidate_panel_final =
      mofa_pass &
      n_supervised_pass >= 2 &
      n_core_pass >= 1,
    
    consensus_score =
      as.numeric(mofa_pass) +
      dplyr::coalesce(score_mofa_percentile_view, 0) +
      weighted_supervised_score
  ) %>%
  dplyr::arrange(
    dplyr::desc(candidate_panel_final),
    dplyr::desc(n_supervised_pass),
    dplyr::desc(weighted_supervised_score),
    dplyr::desc(score_mofa_percentile_view),
    dplyr::desc(consensus_score),
    rank_mofa_for_selection
  )

panel_final <- panel_feature_tbl %>%
  dplyr::filter(candidate_panel_final)

cat("\n============================================================\n")
cat("PANEL FINAL: MOFA + >=2/3 MODELOS SUPERVISADOS\n")
cat("============================================================\n")
.print_tbl(panel_final, n = Inf)

cat("\nNúmero de features panel final:\n")
print(nrow(panel_final))

cat("\n============================================================\n")
cat("TOP CONSENSO GLOBAL\n")
cat("============================================================\n")
.print_tbl(
  panel_feature_tbl %>%
    dplyr::select(
      feature_model,
      view,
      feature_label,
      rank_mofa_for_selection,
      score_mofa_for_selection,
      score_mofa_percentile_view,
      mofa_pass,
      splsda_freq,
      splsda_pass,
      lasso_freq,
      lasso_pass,
      plsda_freq,
      plsda_pass,
      n_supervised_pass,
      n_core_pass,
      weighted_supervised_score,
      candidate_panel_final,
      consensus_score
    ) %>%
    dplyr::slice_head(n = 50),
  n = 50
)

## ============================================================
## 10) AUDITORÍA TEST DEL PANEL FINAL
## ============================================================

fit_eval_models_on_panel <- function(panel_features, panel_name, df_train, df_test) {
  
  panel_features <- intersect(
    panel_features,
    setdiff(colnames(df_train), NON_BIOMARKER_COLS)
  )
  
  if (length(panel_features) < 2) {
    warning("Panel ", panel_name, " tiene menos de 2 features. No se evalúa.")
    return(NULL)
  }
  
  train_sub <- df_train[, c("y", panel_features), drop = FALSE]
  test_sub  <- df_test[,  c("y", panel_features), drop = FALSE]
  
  y_test <- test_sub$y
  
  X_train <- as.matrix(train_sub[, -1, drop = FALSE])
  X_test  <- as.matrix(test_sub[, -1, drop = FALSE])
  
  out <- list()
  
  ## ------------------------------------------------------------
  ## TEST audit 1: LASSO multinomial
  ## ------------------------------------------------------------
  fit_lasso <- tryCatch(
    glmnet::cv.glmnet(
      x = X_train,
      y = train_sub$y,
      family = "multinomial",
      alpha = 1,
      type.measure = "class",
      nfolds = min(3, min(table(train_sub$y)))
    ),
    error = function(e) e
  )
  
  if (!inherits(fit_lasso, "error")) {
    prob_lasso <- predict(
      fit_lasso,
      newx = X_test,
      s = "lambda.min",
      type = "response"
    )
    
    prob_lasso <- as.matrix(prob_lasso[, , 1])
    prob_lasso <- prob_lasso[, levels(y_test), drop = FALSE]
    
    out[["LASSO"]] <- metrics_from_probs(prob_lasso, y_test) %>%
      dplyr::mutate(
        Panel = panel_name,
        Model = "LASSO",
        n_features = length(panel_features)
      )
  }
  
  ## ------------------------------------------------------------
  ## TEST audit 2: Multinomial clásico
  ## ------------------------------------------------------------
  fit_multinom <- tryCatch(
    nnet::multinom(
      y ~ .,
      data = train_sub,
      trace = FALSE,
      MaxNWts = 10000,
      maxit = 1000
    ),
    error = function(e) e
  )
  
  if (!inherits(fit_multinom, "error")) {
    
    prob_multinom <- predict(
      fit_multinom,
      newdata = test_sub,
      type = "probs"
    )
    
    prob_multinom <- as.matrix(prob_multinom)
    
    if (is.null(colnames(prob_multinom))) {
      colnames(prob_multinom) <- levels(y_test)
    }
    
    prob_multinom <- prob_multinom[, levels(y_test), drop = FALSE]
    
    out[["Multinom_classic"]] <- metrics_from_probs(prob_multinom, y_test) %>%
      dplyr::mutate(
        Panel = panel_name,
        Model = "Multinom_classic",
        n_features = length(panel_features)
      )
  }
  dplyr::bind_rows(out)
}
test_panel_metrics <- dplyr::bind_rows(
  fit_eval_models_on_panel(
    panel_final$feature_model,
    "panel_final",
    df_train,
    df_test
  )
)


## ------------------------------------------------------------
## Sensibilidad/especificidad por clase para modelos TEST
## ------------------------------------------------------------

test_panel_class_metrics <- tibble::tibble()

## Recalcular predicciones por modelo de auditoría TEST para métricas por clase
get_test_class_metrics_on_panel <- function(panel_features, panel_name, df_train, df_test) {
  
  panel_features <- intersect(
    panel_features,
    setdiff(colnames(df_train), NON_BIOMARKER_COLS)
  )
  
  if (length(panel_features) < 2) return(tibble::tibble())
  
  train_sub <- df_train[, c("y", panel_features), drop = FALSE]
  test_sub  <- df_test[,  c("y", panel_features), drop = FALSE]
  y_test <- test_sub$y
  
  X_train <- as.matrix(train_sub[, -1, drop = FALSE])
  X_test  <- as.matrix(test_sub[, -1, drop = FALSE])
  
  out <- list()
  
  fit_lasso <- tryCatch(
    glmnet::cv.glmnet(
      x = X_train,
      y = train_sub$y,
      family = "multinomial",
      alpha = 1,
      type.measure = "class",
      nfolds = min(3, min(table(train_sub$y)))
    ),
    error = function(e) e
  )
  
  if (!inherits(fit_lasso, "error")) {
    prob_lasso <- predict(
      fit_lasso,
      newx = X_test,
      s = "lambda.min",
      type = "response"
    )
    prob_lasso <- as.matrix(prob_lasso[, , 1])
    prob_lasso <- prob_lasso[, levels(y_test), drop = FALSE]
    pred_lasso <- factor(colnames(prob_lasso)[max.col(prob_lasso)], levels = levels(y_test))
    
    out[["LASSO"]] <- class_metrics_from_class(pred_lasso, y_test) %>%
      dplyr::mutate(Panel = panel_name, Model = "LASSO")
  }
  
  fit_multinom <- tryCatch(
    nnet::multinom(
      y ~ .,
      data = train_sub,
      trace = FALSE,
      MaxNWts = 10000,
      maxit = 1000
    ),
    error = function(e) e
  )
  
  if (!inherits(fit_multinom, "error")) {
    prob_multinom <- predict(fit_multinom, newdata = test_sub, type = "probs")
    prob_multinom <- as.matrix(prob_multinom)
    if (is.null(colnames(prob_multinom))) colnames(prob_multinom) <- levels(y_test)
    prob_multinom <- prob_multinom[, levels(y_test), drop = FALSE]
    pred_multinom <- factor(colnames(prob_multinom)[max.col(prob_multinom)], levels = levels(y_test))
    
    out[["Multinom_classic"]] <- class_metrics_from_class(pred_multinom, y_test) %>%
      dplyr::mutate(Panel = panel_name, Model = "Multinom_classic")
  }
  
  dplyr::bind_rows(out)
}

test_panel_class_metrics <- get_test_class_metrics_on_panel(
  panel_features = panel_final$feature_model,
  panel_name = "panel_final",
  df_train = df_train,
  df_test = df_test
)
cat("\n============================================================\n")
cat("AUDITORÍA TEST DEL PANEL FINAL\n")
cat("============================================================\n")
.print_tbl(test_panel_metrics, n = Inf)


## ============================================================
## 10B) BRMS CONFIRMATORIO DEL PANEL FINAL
## ------------------------------------------------------------
## Un modelo:
##   1) panel_final_all
## No participa en selección.
## ============================================================

brms_confirmatory <- list()
brms_confirmatory_metrics <- tibble::tibble()
brms_confirmatory_posteriors <- tibble::tibble()
brms_confirmatory_predictions <- tibble::tibble()
brms_confirmatory_confusion <- tibble::tibble()
brms_confirmatory_class_metrics <- tibble::tibble()
brms_confirmatory_signature <- tibble::tibble()
brms_confirmatory_contrasts <- tibble::tibble()

brms_confirmatory_plot_manifest <- tibble::tibble(
  plot = character(),
  description = character()
)

.safe_file_tag <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

if (isTRUE(DO_BRMS_CONFIRMATORY)) {
  
  cat("\n============================================================\n")
  cat("BRMS CONFIRMATORIO DEL PANEL FINAL\n")
  cat("============================================================\n")
  
  brms_variants <- make_brms_panel_variants(
    panel_final = panel_final,
    top_n = BRMS_CONFIRM_TOP_N,
    min_features = BRMS_CONFIRM_MIN_FEATURES
  )
  
  for (variant_name in names(brms_variants)) {
    
    brms_panel_obj <- brms_variants[[variant_name]]
    
    if (length(brms_panel_obj$features) < BRMS_CONFIRM_MIN_FEATURES) {
      warning("No se ajusta BRMS para ", variant_name, ": pocas features.")
      next
    }
    
    file_tag <- .safe_file_tag(brms_panel_obj$panel_name)
    
    brms_panel_features <- intersect(
      brms_panel_obj$features,
      setdiff(colnames(df_train), NON_BIOMARKER_COLS)
    )
    
    cat("\n------------------------------------------------------------\n")
    cat("BRMS variante:", brms_panel_obj$panel_name, "\n")
    cat("Número inicial de features:", length(brms_panel_features), "\n")
    cat("------------------------------------------------------------\n")
    
    brms_train <- df_train[, c("y", brms_panel_features), drop = FALSE]
    brms_test  <- df_test[,  c("y", brms_panel_features), drop = FALSE]
    
    sds_brms <- vapply(brms_train[, -1, drop = FALSE], sd, numeric(1), na.rm = TRUE)
    keep_brms <- is.finite(sds_brms) & sds_brms > 0
    
    if (!all(keep_brms)) {
      cat("\n--- Variables eliminadas del BRMS por sd=0 ---\n")
      print(names(keep_brms)[!keep_brms])
    }
    
    brms_panel_features <- names(keep_brms)[keep_brms]
    
    if (length(brms_panel_features) < BRMS_CONFIRM_MIN_FEATURES) {
      warning("No se ajusta BRMS para ", variant_name, ": pocas features tras sd=0.")
      next
    }
    
    brms_train <- brms_train[, c("y", brms_panel_features), drop = FALSE]
    brms_test  <- brms_test[,  c("y", brms_panel_features), drop = FALSE]
    
    brms_train$y <- factor(brms_train$y)
    brms_train$y <- stats::relevel(brms_train$y, ref = REF_CLASS)
    brms_test$y <- factor(brms_test$y, levels = levels(brms_train$y))
    
    brms_panel_tbl_final <- brms_panel_obj$panel_tbl %>%
      dplyr::filter(feature_model %in% brms_panel_features)
    
    cat("\n--- Features usadas finalmente en BRMS ---\n")
    .print_tbl(
      brms_panel_tbl_final %>%
        dplyr::select(
          feature_model,
          view,
          feature_label,
          score_mofa_percentile_view,
          splsda_freq,
          lasso_freq,
          plsda_freq,
          n_supervised_pass,
          weighted_supervised_score,
          consensus_score
        ),
      n = Inf
    )
    
    backend_brms <- .detect_backend_safe()
    ctrl_brms <- .make_ctrl_safe(backend_brms)
    
    priors_brms <- make_multinom_priors(
      y_levels = levels(brms_train$y),
      ref = REF_CLASS,
      sd_b = BRMS_CONFIRM_PRIOR_SD_B,
      sd_int = BRMS_CONFIRM_PRIOR_SD_INT
    )
    
    fit_brms_confirmatory <- tryCatch(
      brms::brm(
        y ~ .,
        data = brms_train,
        family = brms::categorical(),
        prior = priors_brms,
        chains = BRMS_CONFIRM_CHAINS,
        iter = BRMS_CONFIRM_ITER,
        warmup = BRMS_CONFIRM_WARMUP,
        refresh = 0,
        backend = backend_brms,
        control = ctrl_brms,
        cores = min(BRMS_CONFIRM_CHAINS, max(1, N_CORES)),
        seed = 123,
        save_pars = brms::save_pars(all = TRUE)
      ),
      error = function(e) e
    )
    
    if (inherits(fit_brms_confirmatory, "error")) {
      warning("Falló BRMS ", brms_panel_obj$panel_name, ": ",
              conditionMessage(fit_brms_confirmatory))
      next
    }
    
    P_brms_confirm <- brms_probs(
      fit_brms_confirmatory,
      brms_test
    )[, levels(brms_test$y), drop = FALSE]
    
    pred_brms_confirm <- factor(
      colnames(P_brms_confirm)[max.col(P_brms_confirm, ties.method = "first")],
      levels = levels(brms_test$y)
    )
    
    brms_metrics_one <- metrics_from_probs(
      P_brms_confirm,
      brms_test$y
    ) %>%
      dplyr::mutate(
        Panel = brms_panel_obj$panel_name,
        Model = paste0("BRMS_", variant_name),
        n_features = length(brms_panel_features),
        backend = backend_brms
      )
    
    prob_brms_confirm_df <- as.data.frame(P_brms_confirm, check.names = FALSE)
    colnames(prob_brms_confirm_df) <- paste0("prob_", colnames(prob_brms_confirm_df))
    
    brms_predictions_one <- dplyr::bind_cols(
      tibble::tibble(
        id = rownames(brms_test),
        real = brms_test$y,
        predicho = pred_brms_confirm,
        correcto = brms_test$y == pred_brms_confirm,
        prob_max = apply(P_brms_confirm, 1, max),
        Panel = brms_panel_obj$panel_name,
        Model = paste0("BRMS_", variant_name)
      ),
      prob_brms_confirm_df
    )
    
    brms_confusion_one <- as.data.frame(
      table(
        Prediction = pred_brms_confirm,
        Reference = brms_test$y
      )
    ) %>%
      dplyr::mutate(
        Panel = brms_panel_obj$panel_name,
        Model = paste0("BRMS_", variant_name)
      )
    
    brms_class_one <- class_metrics_from_class(
      pred = pred_brms_confirm,
      obs = brms_test$y
    ) %>%
      dplyr::mutate(
        Panel = brms_panel_obj$panel_name,
        Model = paste0("BRMS_", variant_name)
      )
    
    smry_brms <- summary(fit_brms_confirmatory)
    
    diag_brms <- list(
      rhat = max(smry_brms$fixed[, "Rhat"], na.rm = TRUE),
      ess_bulk = min(smry_brms$fixed[, "Bulk_ESS"], na.rm = TRUE),
      ess_tail = min(smry_brms$fixed[, "Tail_ESS"], na.rm = TRUE)
    )
    
    loo_brms_confirmatory <- tryCatch(
      loo::loo(fit_brms_confirmatory, moment_match = FALSE),
      error = function(e) NULL
    )
    
    waic_brms_confirmatory <- tryCatch(
      loo::waic(fit_brms_confirmatory),
      error = function(e) NULL
    )
    
    if (!is.null(loo_brms_confirmatory)) {
      pareto_k <- loo::pareto_k_values(loo_brms_confirmatory)
    } else {
      pareto_k <- rep(NA_real_, nrow(brms_train))
    }
    
    pareto_diag <- tibble::tibble(
      pareto_k_max = ifelse(all(is.na(pareto_k)), NA_real_, max(pareto_k, na.rm = TRUE)),
      pareto_k_mean = ifelse(all(is.na(pareto_k)), NA_real_, mean(pareto_k, na.rm = TRUE)),
      pareto_k_n_gt_0_7 = sum(pareto_k > 0.7, na.rm = TRUE),
      pareto_k_n_gt_1 = sum(pareto_k > 1, na.rm = TRUE)
    )
    
    brms_metrics_one <- brms_metrics_one %>%
      dplyr::mutate(
        rhat = diag_brms$rhat,
        ess_bulk = diag_brms$ess_bulk,
        ess_tail = diag_brms$ess_tail,
        elpd_loo = if (!is.null(loo_brms_confirmatory)) loo_brms_confirmatory$estimates["elpd_loo", "Estimate"] else NA_real_,
        p_loo = if (!is.null(loo_brms_confirmatory)) loo_brms_confirmatory$estimates["p_loo", "Estimate"] else NA_real_,
        looic = if (!is.null(loo_brms_confirmatory)) loo_brms_confirmatory$estimates["looic", "Estimate"] else NA_real_,
        pareto_k_max = pareto_diag$pareto_k_max,
        pareto_k_mean = pareto_diag$pareto_k_mean,
        pareto_k_n_gt_0_7 = pareto_diag$pareto_k_n_gt_0_7,
        pareto_k_n_gt_1 = pareto_diag$pareto_k_n_gt_1,
        elpd_waic = if (!is.null(waic_brms_confirmatory)) waic_brms_confirmatory$estimates["elpd_waic", "Estimate"] else NA_real_,
        p_waic = if (!is.null(waic_brms_confirmatory)) waic_brms_confirmatory$estimates["p_waic", "Estimate"] else NA_real_,
        waic = if (!is.null(waic_brms_confirmatory)) waic_brms_confirmatory$estimates["waic", "Estimate"] else NA_real_
      )
    
    brms_posteriors_one <- extract_brms_confirmatory_posteriors(
      fit = fit_brms_confirmatory,
      panel_tbl = brms_panel_tbl_final
    ) %>%
      dplyr::mutate(
        Panel = brms_panel_obj$panel_name,
        Model = paste0("BRMS_", variant_name)
      )
    
    ## ------------------------------------------------------------
    ## Contraste posterior directo SO - FD
    ## Positivo: mayor log-odds hacia SO que hacia FD
    ## Negativo: mayor log-odds hacia FD que hacia SO
    ## ------------------------------------------------------------
    
    draws_raw_contrast <- as.data.frame(posterior::as_draws_df(fit_brms_confirmatory))
    
    brms_contrast_so_fd_one <- dplyr::bind_rows(lapply(brms_panel_features, function(feat) {
      
      col_so <- paste0("b_muSO_", feat)
      col_fd <- paste0("b_muFD_", feat)
      
      if (!all(c(col_so, col_fd) %in% colnames(draws_raw_contrast))) {
        return(NULL)
      }
      
      delta <- as.numeric(draws_raw_contrast[[col_so]]) -
        as.numeric(draws_raw_contrast[[col_fd]])
      
      tibble::tibble(
        Panel = brms_panel_obj$panel_name,
        Model = paste0("BRMS_", variant_name),
        contrast = "SO_vs_FD",
        feature_model = feat,
        Estimate = median(delta, na.rm = TRUE),
        Q2.5 = as.numeric(stats::quantile(delta, 0.025, na.rm = TRUE)),
        Q97.5 = as.numeric(stats::quantile(delta, 0.975, na.rm = TRUE)),
        prob_SO_gt_FD = mean(delta > 0, na.rm = TRUE),
        prob_FD_gt_SO = mean(delta < 0, na.rm = TRUE),
        ci_excludes_zero = (
          stats::quantile(delta, 0.025, na.rm = TRUE) > 0 |
            stats::quantile(delta, 0.975, na.rm = TRUE) < 0
        ),
        posterior_direction = dplyr::case_when(
          mean(delta > 0, na.rm = TRUE) >= 0.95 ~ "SO > FD",
          mean(delta < 0, na.rm = TRUE) >= 0.95 ~ "FD > SO",
          TRUE ~ "uncertain"
        )
      )
    })) %>%
      dplyr::left_join(
        brms_panel_tbl_final %>%
          dplyr::select(
            dplyr::any_of(c(
              "feature_model",
              "view",
              "view_code",
              "feature_label",
              "score_mofa_percentile_view",
              "splsda_freq",
              "lasso_freq",
              "plsda_freq",
              "n_supervised_pass",
              "weighted_supervised_score",
              "consensus_score"
            ))
          ),
        by = "feature_model"
      )
    
    ## Firma final por variable
    brms_signature_one <- brms_posteriors_one %>%
      dplyr::select(
        Panel,
        Model,
        feature_model,
        view,
        view_code,
        feature_label,
        class,
        Estimate,
        prob_positive,
        prob_negative
      ) %>%
      tidyr::pivot_wider(
        names_from = class,
        values_from = c(Estimate, prob_positive, prob_negative),
        names_sep = "_"
      ) %>%
      dplyr::mutate(
        beta_FD = dplyr::coalesce(Estimate_FD, NA_real_),
        beta_SO = dplyr::coalesce(Estimate_SO, NA_real_),
        signature_group = dplyr::case_when(
          is.finite(beta_FD) & is.finite(beta_SO) & beta_FD < 0 & beta_SO < 0 ~ "NP-like",
          is.finite(beta_FD) & is.finite(beta_SO) & beta_FD > 0 & beta_SO > 0 ~ "non-NP-like",
          is.finite(beta_FD) & beta_FD > 0 & (is.na(beta_SO) | abs(beta_FD) >= abs(beta_SO)) ~ "FD-like",
          is.finite(beta_SO) & beta_SO > 0 & (is.na(beta_FD) | abs(beta_SO) > abs(beta_FD)) ~ "SO-like",
          TRUE ~ "uncertain"
        )
      )
    
    ## Acumular tablas
    brms_confirmatory_metrics <- dplyr::bind_rows(brms_confirmatory_metrics, brms_metrics_one)
    brms_confirmatory_predictions <- dplyr::bind_rows(brms_confirmatory_predictions, brms_predictions_one)
    brms_confirmatory_confusion <- dplyr::bind_rows(brms_confirmatory_confusion, brms_confusion_one)
    brms_confirmatory_class_metrics <- dplyr::bind_rows(brms_confirmatory_class_metrics, brms_class_one)
    brms_confirmatory_posteriors <- dplyr::bind_rows(brms_confirmatory_posteriors, brms_posteriors_one)
    brms_confirmatory_signature <- dplyr::bind_rows(brms_confirmatory_signature, brms_signature_one)
    brms_confirmatory_contrasts <- dplyr::bind_rows(brms_confirmatory_contrasts, brms_contrast_so_fd_one)
    ## Añadir BRMS a auditoría TEST global
    test_panel_metrics <- dplyr::bind_rows(
      test_panel_metrics,
      brms_metrics_one %>%
        dplyr::select(dplyr::any_of(colnames(test_panel_metrics)))
    )
    
    test_panel_class_metrics <- dplyr::bind_rows(
      test_panel_class_metrics,
      brms_class_one
    )
    
    ## Guardar RDS específicos
    saveRDS(
      fit_brms_confirmatory,
      out_rds(paste0("fit_brms_confirmatory_", file_tag, ".rds"))
    )
    
    saveRDS(
      loo_brms_confirmatory,
      out_rds(paste0("loo_brms_confirmatory_", file_tag, ".rds"))
    )
    
    saveRDS(
      waic_brms_confirmatory,
      out_rds(paste0("waic_brms_confirmatory_", file_tag, ".rds"))
    )
    
    brms_confirmatory[[variant_name]] <- list(
      panel_name = brms_panel_obj$panel_name,
      panel_features = brms_panel_features,
      panel_tbl = brms_panel_tbl_final,
      fit = fit_brms_confirmatory,
      metrics = brms_metrics_one,
      posteriors = brms_posteriors_one,
      predictions = brms_predictions_one,
      confusion = brms_confusion_one,
      class_metrics = brms_class_one,
      signature = brms_signature_one,
      contrasts = brms_contrast_so_fd_one,
      loo = loo_brms_confirmatory,
      waic = waic_brms_confirmatory
    )
  }
  
  ## Guardar tablas agregadas BRMS
  write.csv(brms_confirmatory_metrics, out_table("12_brms_confirmatory_metrics.csv"), row.names = FALSE)
  write.csv(brms_confirmatory_posteriors, out_table("13_brms_confirmatory_posteriors.csv"), row.names = FALSE)
  write.csv(brms_confirmatory_signature, out_table("13b_brms_confirmatory_signature.csv"), row.names = FALSE)
  write.csv(brms_confirmatory_contrasts, out_table("13c_brms_confirmatory_contrasts_SO_FD.csv"), row.names = FALSE)
  write.csv(brms_confirmatory_predictions, out_table("14_brms_confirmatory_predictions.csv"), row.names = FALSE)
  write.csv(brms_confirmatory_confusion, out_table("15_brms_confirmatory_confusion.csv"), row.names = FALSE)
  write.csv(brms_confirmatory_class_metrics, out_table("15b_brms_confirmatory_class_metrics.csv"), row.names = FALSE)
  
  saveRDS(
    brms_confirmatory,
    out_rds("brms_confirmatory_object.rds")
  )
}
## ============================================================
## 11) EXPORTAR RESULTADOS
## ============================================================


write.csv(
  test_panel_class_metrics,
  out_table("11b_test_panel_class_metrics.csv"),
  row.names = FALSE
)


write.csv(fold_audit, out_table("00_audit_cv_folds.csv"), row.names = FALSE)
write.csv(splsda_seed_audit,
          out_table("00_audit_splsda_seeds.csv"),
          row.names = FALSE)
write.csv(lasso_audit,
          out_table("00_audit_LASSO_seeds.csv"),
          row.names = FALSE)
write.csv(plsda_audit,
          out_table("00_audit_plsda_folds.csv"),
          row.names = FALSE)

write.csv(mofa_pass_tbl,
          out_table("01_mofa_features_input.csv"),
          row.names = FALSE)
write.csv(splsda_feature_tbl,
          out_table("02_splsda_feature_stability.csv"),
          row.names = FALSE)
write.csv(lasso_feature_tbl,
          out_table("03_LASSO_feature_stability.csv"),
          row.names = FALSE)
write.csv(plsda_feature_tbl,
          out_table("05_plsda_feature_stability.csv"),
          row.names = FALSE)

write.csv(panel_feature_tbl,
          out_table("06_panel_feature_consensus_all.csv"),
          row.names = FALSE)
write.csv(panel_final,
          out_table("07_PANEL_FINAL_mofa_2of3.csv"),
          row.names = FALSE)
write.csv(cv_metrics_all,
          out_table("09_cv_metrics_all_models_long.csv"),
          row.names = FALSE)
write.csv(cv_metrics_summary,
          out_table("10_cv_metrics_summary_by_model.csv"),
          row.names = FALSE)
write.csv(test_panel_metrics,
          out_table("11_test_panel_metrics_audit.csv"),
          row.names = FALSE)
## ============================================================
## 12) PLOTS DE AUDITORÍA DEL PANEL FINAL MULTIMODELO
## No reajusta modelos. Usa objetos ya calculados.
## ============================================================

cat("\n============================================================\n")
cat("GENERANDO PLOTS DE AUDITORÍA PANEL FINAL MULTIMODELO\n")
cat("============================================================\n")

METRIC_PANEL_PLOTS_DIR <- out_plot("metricas_panel")
dir.create(METRIC_PANEL_PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

out_panel_plot <- function(filename) {
  file.path(METRIC_PANEL_PLOTS_DIR, filename)
}

.save_panel_plot <- function(plot, filename, width = 8, height = 6) {
  ggsave(
    filename = out_panel_plot(filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

panel_plot_manifest <- tibble::tibble(
  plot = character(),
  description = character()
)

.add_panel_manifest <- function(filename, description) {
  panel_plot_manifest <<- dplyr::bind_rows(
    panel_plot_manifest,
    tibble::tibble(
      plot = filename,
      description = description
    )
  )
}

## -----------------------------
## 12.1 Auditoría folds CV
## -----------------------------

if (exists("fold_audit") && nrow(fold_audit) > 0) {
  
  p_folds <- fold_audit %>%
    ggplot(aes(x = FoldID, y = Class, fill = n)) +
    geom_tile() +
    geom_text(aes(label = n), size = 3) +
    labs(
      title = "Auditoría de folds CV",
      subtitle = "Número de muestras por clase en cada validación",
      x = "Fold / repetición",
      y = "Clase",
      fill = "n"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  print(p_folds)
  .save_panel_plot(p_folds, "01_auditoria_folds_cv.png", 12, 5)
  .add_panel_manifest(
    "01_auditoria_folds_cv.png",
    "Distribución de clases en los folds de validación."
  )
}

## -----------------------------
## 12.2 Éxito/fallo por modelo
## -----------------------------

audit_models_tbl <- dplyr::bind_rows(if (exists("splsda_seed_audit")) {
  splsda_seed_audit %>%
    dplyr::transmute(
      Model = "sPLSDA",
      unit = as.character(seed),
      ok = ok,
      error = error
    )
}, if (exists("lasso_audit")) {
  lasso_audit %>%
    dplyr::transmute(
      Model = "LASSO",
      unit = as.character(seed),
      ok = ok,
      error = error
    )
},  if (exists("plsda_audit")) {
  plsda_audit %>%
    dplyr::transmute(
      Model = "PLSDA",
      unit = as.character(FoldID),
      ok = ok,
      error = error
    )
})

if (nrow(audit_models_tbl) > 0) {
  
  audit_success_tbl <- audit_models_tbl %>%
    dplyr::group_by(Model) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_ok = sum(ok, na.rm = TRUE),
      n_fail = sum(!ok, na.rm = TRUE),
      prop_ok = n_ok / n_total,
      .groups = "drop"
    )
  
  p_audit_success <- audit_success_tbl %>%
    ggplot(aes(x = Model, y = prop_ok)) +
    geom_col() +
    geom_text(aes(label = paste0(n_ok, "/", n_total)), vjust = -0.4, size = 4) +
    coord_cartesian(ylim = c(0, 1.05)) +
    labs(
      title = "Éxito de ajuste por modelo",
      subtitle = "Proporción de semillas/folds que corrieron correctamente",
      x = "Modelo",
      y = "Proporción OK"
    ) +
    theme_minimal()
  
  print(p_audit_success)
  .save_panel_plot(p_audit_success, "02_auditoria_exito_modelos.png", 7, 5)
  .add_panel_manifest(
    "02_auditoria_exito_modelos.png",
    "Proporción de ejecuciones correctas por modelo."
  )
}

## -----------------------------
## 12.3 Métricas CV promedio por modelo
## -----------------------------

if (exists("cv_metrics_summary") && nrow(cv_metrics_summary) > 0) {
  
  cv_metric_cols <- intersect(
    c(
      "mean_Accuracy",
      "mean_BalAcc",
      "mean_BER",
      "mean_MacroSensitivity",
      "mean_MacroSpecificity",
      "mean_MacroF1",
      "mean_AUC",
      "mean_LogLoss"
    ),
    colnames(cv_metrics_summary)
  )
  
  p_cv_summary <- cv_metrics_summary %>%
    dplyr::select(Model, dplyr::all_of(cv_metric_cols)) %>%
    tidyr::pivot_longer(
      cols = -Model,
      names_to = "metric",
      values_to = "value"
    ) %>%
    ggplot(aes(x = Model, y = value)) +
    geom_col() +
    facet_wrap(~ metric, scales = "free_y") +
    labs(
      title = "Resumen CV por modelo",
      subtitle = "Métricas promedio sobre CV/repeticiones",
      x = "Modelo",
      y = "Valor"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_cv_summary)
  .save_panel_plot(p_cv_summary, "03_cv_metricas_promedio_por_modelo.png", 10, 7)
  .add_panel_manifest(
    "03_cv_metricas_promedio_por_modelo.png",
    "Métricas CV promedio por modelo."
  )
}

## -----------------------------
## 12.4 Distribución de métricas CV
## -----------------------------

if (exists("cv_metrics_all") && nrow(cv_metrics_all) > 0) {
  
  cv_long_metric_cols <- intersect(
    c(
      "Accuracy",
      "BalAcc",
      "BER",
      "MacroSensitivity",
      "MacroSpecificity",
      "MacroF1",
      "AUC",
      "LogLoss"
    ),
    colnames(cv_metrics_all)
  )
  
  p_cv_dist <- cv_metrics_all %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(cv_long_metric_cols),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::filter(is.finite(value)) %>%
    ggplot(aes(x = Model, y = value)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.8) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(
      title = "Distribución de métricas CV",
      subtitle = "Cada punto representa una semilla, fold o repetición disponible",
      x = "Modelo",
      y = "Valor"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_cv_dist)
  .save_panel_plot(p_cv_dist, "04_cv_metricas_distribucion_por_modelo.png", 11, 7)
  .add_panel_manifest(
    "04_cv_metricas_distribucion_por_modelo.png",
    "Distribución de métricas CV por modelo."
  )
}

## -----------------------------
## 12.5 Frecuencia de selección por modelo
## -----------------------------

if (exists("panel_feature_tbl") && nrow(panel_feature_tbl) > 0) {
  
  freq_cols <- intersect(
    c("splsda_freq", "lasso_freq", "plsda_freq"),
    colnames(panel_feature_tbl)
  )
  if (length(freq_cols) > 0) {
    
    panel_feature_tbl_for_freq <- panel_feature_tbl %>%
      dplyr::mutate(
        mean_supervised_freq = rowMeans(
          dplyr::select(., dplyr::all_of(freq_cols)),
          na.rm = TRUE
        )
      )
    
    top_features_freq <- panel_feature_tbl_for_freq %>%
      dplyr::arrange(
        dplyr::desc(candidate_panel_final),
        dplyr::desc(n_supervised_pass),
        dplyr::desc(weighted_supervised_score),
        dplyr::desc(mean_supervised_freq),
        dplyr::desc(score_mofa_percentile_view)
      ) %>%
      dplyr::slice_head(n = 40) %>%
      dplyr::pull(feature_model)
    
    feature_freq_long <- panel_feature_tbl_for_freq %>%
      dplyr::select(
        feature_model,
        view,
        feature_label,
        score_mofa_percentile_view,
        consensus_score,
        n_supervised_pass,
        dplyr::all_of(freq_cols)
      ) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(freq_cols),
        names_to = "model",
        values_to = "freq"
      ) %>%
      dplyr::mutate(
        model = dplyr::recode(
          model,
          splsda_freq = "sPLSDA",
          lasso_freq = "LASSO",
          plsda_freq = "PLSDA"
        )
      )
    
    p_freq_heat <- feature_freq_long %>%
      dplyr::filter(feature_model %in% top_features_freq) %>%
      dplyr::mutate(feature_model = factor(feature_model, levels = rev(top_features_freq))) %>%
      ggplot(aes(x = model, y = feature_model, fill = freq)) +
      geom_tile() +
      geom_text(aes(label = sprintf("%.2f", freq)), size = 2.6) +
      labs(
        title = "Frecuencia de selección por modelo supervisado",
        subtitle = "Top 40 features por consenso multimodelo",
        x = "Modelo",
        y = "Feature",
        fill = "Frecuencia"
      ) +
      theme_minimal()
    
    print(p_freq_heat)
    .save_panel_plot(p_freq_heat, "05_frecuencia_seleccion_top40_heatmap.png", 9, 10)
    .add_panel_manifest(
      "05_frecuencia_seleccion_top40_heatmap.png",
      "Heatmap de frecuencia de selección por modelo supervisado."
    )
  }
}

## -----------------------------
## 12.6 Número de modelos que apoyan cada feature
## -----------------------------

if (exists("panel_feature_tbl") && nrow(panel_feature_tbl) > 0 &&
    "n_supervised_pass" %in% colnames(panel_feature_tbl)) {
  
  p_npass <- panel_feature_tbl %>%
    dplyr::count(n_supervised_pass, name = "n_features") %>%
    ggplot(aes(x = factor(n_supervised_pass), y = n_features)) +
    geom_col() +
    geom_text(aes(label = n_features), vjust = -0.4, size = 4) +
    labs(
      title = "Distribución del soporte supervisado",
      subtitle = "Número de modelos supervisados que seleccionan cada feature",
      x = "Número de modelos supervisados que pasan",
      y = "Número de features"
    ) +
    theme_minimal()
  
  print(p_npass)
  .save_panel_plot(p_npass, "06_distribucion_n_modelos_supervisados_pass.png", 7, 5)
  .add_panel_manifest(
    "06_distribucion_n_modelos_supervisados_pass.png",
    "Distribución de features según número de modelos supervisados que las apoyan."
  )
}

## -----------------------------
## 12.7 Consenso binario MOFA + modelos
## -----------------------------

if (exists("panel_feature_tbl") && nrow(panel_feature_tbl) > 0) {
  
  pass_cols <- intersect(
    c("mofa_pass", "splsda_pass", "lasso_pass", "plsda_pass"),
    colnames(panel_feature_tbl)
  )
  if (length(pass_cols) > 0) {
    
    top_consensus_features <- panel_feature_tbl %>%
      dplyr::arrange(
        dplyr::desc(candidate_panel_final),
        dplyr::desc(n_supervised_pass),
        dplyr::desc(weighted_supervised_score),
        dplyr::desc(score_mofa_percentile_view),
        rank_mofa_for_selection
      ) %>%
      dplyr::slice_head(n = 50) %>%
      dplyr::pull(feature_model)
    
    pass_long <- panel_feature_tbl %>%
      dplyr::filter(feature_model %in% top_consensus_features) %>%
      dplyr::select(feature_model, dplyr::all_of(pass_cols)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(pass_cols),
        names_to = "source",
        values_to = "pass"
      ) %>%
      dplyr::mutate(
        source = dplyr::recode(
          source,
          mofa_pass = "MOFA",
          splsda_pass = "sPLSDA",
          lasso_pass = "LASSO",
          plsda_pass = "PLSDA"
        ),
        pass = as.integer(pass),
        feature_model = factor(feature_model, levels = rev(top_consensus_features))
      )
    
    p_pass_heat <- pass_long %>%
      ggplot(aes(x = source, y = feature_model, fill = pass)) +
      geom_tile() +
      geom_text(aes(label = pass), size = 2.8) +
      labs(
        title = "Consenso binario MOFA + modelos supervisados",
        subtitle = "Top 50 features ordenadas por consenso",
        x = "Fuente",
        y = "Feature",
        fill = "Pass"
      ) +
      theme_minimal()
    
    print(p_pass_heat)
    .save_panel_plot(p_pass_heat, "07_consenso_binario_top50_heatmap.png", 9, 11)
    .add_panel_manifest(
      "07_consenso_binario_top50_heatmap.png",
      "Heatmap binario de paso/no paso para MOFA y modelos supervisados."
    )
  }
}

## -----------------------------
## 12.8 Conteo de features por vista y panel
## -----------------------------

if (exists("panel_feature_tbl") && nrow(panel_feature_tbl) > 0) {
  
  panel_count_by_view <- dplyr::bind_rows(
    panel_feature_tbl %>%
      dplyr::mutate(panel_type = "all_input_features"),
    panel_feature_tbl %>%
      dplyr::filter(candidate_panel_final) %>%
      dplyr::mutate(panel_type = "panel_final_mofa_2of3")
  ) %>%
    dplyr::count(panel_type, view, name = "n_features")
  
  p_panel_view <- panel_count_by_view %>%
    ggplot(aes(x = view, y = n_features)) +
    geom_col() +
    geom_text(aes(label = n_features), vjust = -0.4, size = 3.5) +
    facet_wrap(~ panel_type, scales = "free_y") +
    labs(
      title = "Composición del panel por vista",
      x = "Vista",
      y = "Número de features"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_panel_view)
  .save_panel_plot(p_panel_view, "08_panel_conteo_por_vista.png", 10, 6)
  .add_panel_manifest(
    "08_panel_conteo_por_vista.png",
    "Conteo de features por vista para entrada total, panel estricto y panel relajado."
  )
  
  write.csv(
    panel_count_by_view,
    out_table("12_panel_count_by_view_for_plots.csv"),
    row.names = FALSE
  )
}

## -----------------------------
## 12.9 MOFA percentil vs consenso supervisado
## -----------------------------

if (exists("panel_feature_tbl") && nrow(panel_feature_tbl) > 0 &&
    all(c("score_mofa_percentile_view", "n_supervised_pass", "candidate_panel_final") %in% colnames(panel_feature_tbl))) {
  
  p_mofa_consensus <- panel_feature_tbl %>%
    ggplot(aes(
      x = score_mofa_percentile_view,
      y = n_supervised_pass,
      shape = candidate_panel_final
    )) +
    geom_jitter(width = 0.01, height = 0.08, size = 2.8, alpha = 0.8) +
    geom_vline(xintercept = MOFA_PERCENTILE_THR, linetype = 2) +
    scale_y_continuous(breaks = 0:3, limits = c(-0.2, 3.2)) +
    labs(
      title = "Relación entre score MOFA y consenso supervisado",
      subtitle = "Línea vertical = umbral MOFA dentro de vista",
      x = "Percentil MOFA dentro de vista",
      y = "Número de modelos supervisados que pasan",
      shape = "Panel final"
    ) +
    theme_minimal()
  
  print(p_mofa_consensus)
  .save_panel_plot(p_mofa_consensus, "09_mofa_percentil_vs_consenso_supervisado.png", 8, 6)
  .add_panel_manifest(
    "09_mofa_percentil_vs_consenso_supervisado.png",
    "Relación entre percentil MOFA y número de modelos supervisados que apoyan la feature."
  )
}

## -----------------------------
## 12.10 Métricas TEST del panel final
## -----------------------------

if (exists("test_panel_metrics") && !is.null(test_panel_metrics) && nrow(test_panel_metrics) > 0) {
  
  test_metric_cols <- intersect(
    c(
      "Accuracy",
      "BalAcc",
      "BER",
      "MacroSensitivity",
      "MacroSpecificity",
      "MacroF1",
      "AUC",
      "LogLoss"
    ),
    colnames(test_panel_metrics)
  )
  
  p_test_panel <- test_panel_metrics %>%
    dplyr::select(Panel, Model, n_features, dplyr::all_of(test_metric_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(test_metric_cols),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::filter(is.finite(value)) %>%
    ggplot(aes(x = interaction(Panel, Model, sep = " | "), y = value)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.3, size = 3) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(
      title = "Auditoría TEST del panel final",
      subtitle = "TEST solo como auditoría; no usar como selección del panel",
      x = "Panel | Modelo",
      y = "Valor"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_test_panel)
  .save_panel_plot(p_test_panel, "10_test_panel_metricas_auditoria.png", 11, 7)
  .add_panel_manifest(
    "10_test_panel_metricas_auditoria.png",
    "Métricas en TEST para panel estricto y relajado."
  )
  
  ## -----------------------------
  ## 12.10B Sensibilidad/especificidad TEST por clase
  ## -----------------------------
  
  if (exists("test_panel_class_metrics") &&
      !is.null(test_panel_class_metrics) &&
      nrow(test_panel_class_metrics) > 0) {
    
    p_test_class_metrics <- test_panel_class_metrics %>%
      dplyr::select(Panel, Model, Class, Sensitivity, Specificity, Precision, F1) %>%
      tidyr::pivot_longer(
        cols = c(Sensitivity, Specificity, Precision, F1),
        names_to = "metric",
        values_to = "value"
      ) %>%
      ggplot(aes(x = Class, y = value)) +
      geom_col() +
      geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.3, size = 3) +
      facet_grid(metric ~ Model) +
      coord_cartesian(ylim = c(0, 1.05)) +
      labs(
        title = "Auditoría TEST: métricas por clase",
        subtitle = "Sensibilidad, especificidad, precisión y F1 por clase",
        x = "Clase",
        y = "Valor"
      ) +
      theme_minimal()
    
    print(p_test_class_metrics)
    
    .save_panel_plot(
      p_test_class_metrics,
      "10b_test_panel_metricas_por_clase.png",
      12,
      8
    )
    
    .add_panel_manifest(
      "10b_test_panel_metricas_por_clase.png",
      "Sensibilidad, especificidad, precisión y F1 por clase en TEST."
    )
  }
}

## -----------------------------
## 12.11 Ranking panel final
## -----------------------------

if (exists("panel_final") && nrow(panel_final) > 0) {
  
  freq_cols_final <- intersect(
    c("splsda_freq", "lasso_freq", "plsda_freq"),
    colnames(panel_final)
  )
  
  panel_final_rank_plot <- panel_final %>%
    dplyr::mutate(
      mean_supervised_freq = if (length(freq_cols_final) > 0) {
        rowMeans(dplyr::select(., dplyr::all_of(freq_cols_final)), na.rm = TRUE)
      } else {
        NA_real_
      }
    ) %>%
    dplyr::arrange(
      dplyr::desc(n_supervised_pass),
      dplyr::desc(weighted_supervised_score),
      dplyr::desc(score_mofa_percentile_view),
      dplyr::desc(mean_supervised_freq),
      rank_mofa_for_selection
    ) %>%
    dplyr::slice_head(n = 30) %>%
    dplyr::mutate(feature_model = factor(feature_model, levels = rev(feature_model)))
  
  p_panel_final_rank <- panel_final_rank_plot %>%
    ggplot(aes(x = consensus_score, y = feature_model)) +
    geom_col() +
    geom_text(
      aes(label = paste0("sup=", n_supervised_pass, " | w=", sprintf("%.1f", weighted_supervised_score))),
      hjust = -0.1,
      size = 3
    ) +
    labs(
      title = "Top panel final por consenso",
      subtitle = "MOFA + al menos 2/3 modelos supervisados; LASSO ponderado",
      x = "Consensus score",
      y = "Feature"
    ) +
    theme_minimal()
  
  print(p_panel_final_rank)
  
  .save_panel_plot(
    p_panel_final_rank,
    "11_top_panel_final_consenso.png",
    9,
    8
  )
  
  .add_panel_manifest(
    "11_top_panel_final_consenso.png",
    "Ranking visual de las principales features del panel final."
  )
}


## -----------------------------
## 12.12 PCA final del panel_final_all
## -----------------------------

if (exists("panel_final") && nrow(panel_final) > 1) {
  
  pca_features <- intersect(
    panel_final$feature_model,
    setdiff(colnames(df_train), NON_BIOMARKER_COLS)
  )
  
  if (length(pca_features) >= 2) {
    
    X_train_pca <- as.matrix(df_train[, pca_features, drop = FALSE])
    X_test_pca  <- as.matrix(df_test[,  pca_features, drop = FALSE])
    
    ## PCA ajustado SOLO con train; test se proyecta después.
    ## Esto evita fuga de información desde test al PCA.
    pca_panel_all <- prcomp(
      X_train_pca,
      center = TRUE,
      scale. = TRUE
    )
    
    scores_train <- as.data.frame(pca_panel_all$x[, 1:2, drop = FALSE])
    scores_train$id <- rownames(df_train)
    scores_train$y <- df_train$y
    scores_train$Split <- "Train"
    
    scores_test <- as.data.frame(
      predict(pca_panel_all, newdata = X_test_pca)[, 1:2, drop = FALSE]
    )
    scores_test$id <- rownames(df_test)
    scores_test$y <- df_test$y
    scores_test$Split <- "Test"
    
    pca_scores_panel_all <- dplyr::bind_rows(scores_train, scores_test)
    
    var_exp <- summary(pca_panel_all)$importance["Proportion of Variance", 1:2] * 100
    
    p_pca_panel_all <- pca_scores_panel_all %>%
      ggplot(aes(x = PC1, y = PC2, shape = Split, label = id)) +
      geom_point(aes(color = y), size = 3, alpha = 0.9) +
      scale_color_manual(
        values = c(
          "NP" = "black",
          "SO" = "blue",
          "FD" = "red"
        )
      ) +
      geom_text(
        vjust = -0.8,
        size = 2.7,
        show.legend = FALSE
      ) +
      labs(
        title = "PCA final del panel BRMS all",
        subtitle = "PCA ajustado solo en train; test proyectado. Panel = panel_final_all",
        x = paste0("PC1 (", sprintf("%.1f", var_exp[1]), "%)"),
        y = paste0("PC2 (", sprintf("%.1f", var_exp[2]), "%)"),
        color = "Clase",
        shape = "Split"
      ) +
      theme_minimal()
    
    print(p_pca_panel_all)
    
    .save_panel_plot(
      p_pca_panel_all,
      "12_pca_panel_final_all_train_test.png",
      8,
      6
    )
    
    .add_panel_manifest(
      "12_pca_panel_final_all_train_test.png",
      "PCA del panel_final_all ajustado solo con train y proyección de test."
    )
    
    write.csv(
      pca_scores_panel_all,
      out_table("12b_pca_panel_final_all_scores.csv"),
      row.names = FALSE
    )
    
    saveRDS(
      pca_panel_all,
      out_rds("pca_panel_final_all_train_only.rds")
    )
  }
}


## -----------------------------
## 12.13 Métricas BRMS confirmatorio
## -----------------------------

if (exists("brms_confirmatory_metrics") &&
    nrow(brms_confirmatory_metrics) > 0) {
  
  brms_metrics_01 <- brms_confirmatory_metrics %>%
    dplyr::select(
      Panel,
      Model,
      n_features,
      Accuracy,
      BalAcc,
      BER,
      MacroSensitivity,
      MacroSpecificity,
      MacroF1,
      AUC
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        Accuracy,
        BalAcc,
        BER,
        MacroSensitivity,
        MacroSpecificity,
        MacroF1,
        AUC
      ),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::filter(is.finite(value))
  
  p_brms_metrics_01 <- brms_metrics_01 %>%
    ggplot(aes(x = metric, y = value)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.3, size = 3) +
    facet_wrap(~ Model) +
    coord_cartesian(ylim = c(0, 1.05)) +
    labs(
      title = "BRMS confirmatorio: métricas finales 0-1",
      subtitle = "Accuracy, BalAcc, BER, sensibilidad, especificidad, F1 y AUC",
      x = "Métrica",
      y = "Valor"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_brms_metrics_01)
  
  .save_panel_plot(
    p_brms_metrics_01,
    "13_brms_confirmatory_metricas_0_1.png",
    12,
    6
  )
  
  .add_panel_manifest(
    "13_brms_confirmatory_metricas_0_1.png",
    "Métricas BRMS finales en escala 0-1."
  )
  
  p_brms_logloss <- brms_confirmatory_metrics %>%
    ggplot(aes(x = Model, y = LogLoss)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.3f", LogLoss)), vjust = -0.3, size = 3.5) +
    labs(
      title = "BRMS confirmatorio: LogLoss en TEST",
      subtitle = "Menor LogLoss indica mejor calibración probabilística",
      x = "Modelo",
      y = "LogLoss"
    ) +
    theme_minimal()
  
  print(p_brms_logloss)
  
  .save_panel_plot(
    p_brms_logloss,
    "14_brms_confirmatory_logloss.png",
    7,
    5
  )
  
  .add_panel_manifest(
    "14_brms_confirmatory_logloss.png",
    "LogLoss BRMS en TEST."
  )
}


## -----------------------------
## 12.13B Diagnóstico bayesiano BRMS
## -----------------------------

if (exists("brms_confirmatory_metrics") &&
    nrow(brms_confirmatory_metrics) > 0) {
  
  brms_diag_cols <- intersect(
    c(
      "rhat",
      "ess_bulk",
      "ess_tail",
      "pareto_k_max",
      "pareto_k_mean",
      "pareto_k_n_gt_0_7",
      "pareto_k_n_gt_1",
      "elpd_loo",
      "looic",
      "waic"
    ),
    colnames(brms_confirmatory_metrics)
  )
  
  p_brms_diag <- brms_confirmatory_metrics %>%
    dplyr::select(Panel, Model, dplyr::all_of(brms_diag_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(brms_diag_cols),
      names_to = "diagnostic",
      values_to = "value"
    ) %>%
    dplyr::filter(is.finite(value)) %>%
    ggplot(aes(x = Model, y = value)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.3, size = 3) +
    facet_wrap(~ diagnostic, scales = "free_y") +
    labs(
      title = "BRMS confirmatorio: diagnóstico bayesiano",
      subtitle = "Rhat, ESS, Pareto-k, LOO y WAIC del modelo confirmatorio",
      x = "Modelo",
      y = "Valor"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_brms_diag)
  
  .save_panel_plot(
    p_brms_diag,
    "14b_brms_confirmatory_diagnostico_bayesiano.png",
    12,
    7
  )
  
  .add_panel_manifest(
    "14b_brms_confirmatory_diagnostico_bayesiano.png",
    "Diagnóstico bayesiano BRMS: Rhat, ESS, Pareto-k, LOO y WAIC."
  )
}
## -----------------------------
## 12.14 BRMS sensibilidad/especificidad por clase
## -----------------------------

if (exists("brms_confirmatory_class_metrics") &&
    nrow(brms_confirmatory_class_metrics) > 0) {
  
  p_brms_class_metrics <- brms_confirmatory_class_metrics %>%
    dplyr::select(Panel, Model, Class, Sensitivity, Specificity, Precision, F1) %>%
    tidyr::pivot_longer(
      cols = c(Sensitivity, Specificity, Precision, F1),
      names_to = "metric",
      values_to = "value"
    ) %>%
    ggplot(aes(x = Class, y = value)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.3, size = 3) +
    facet_grid(metric ~ Model) +
    coord_cartesian(ylim = c(0, 1.05)) +
    labs(
      title = "BRMS confirmatorio: métricas por clase",
      subtitle = "Sensibilidad, especificidad, precisión y F1 por clase",
      x = "Clase",
      y = "Valor"
    ) +
    theme_minimal()
  
  print(p_brms_class_metrics)
  
  .save_panel_plot(
    p_brms_class_metrics,
    "15_brms_confirmatory_metricas_por_clase.png",
    12,
    8
  )
  
  .add_panel_manifest(
    "15_brms_confirmatory_metricas_por_clase.png",
    "Métricas BRMS por clase."
  )
}
## -----------------------------
## 12.15 Forest plot BRMS con interpretación de firma
## -----------------------------

if (exists("brms_confirmatory_posteriors") &&
    nrow(brms_confirmatory_posteriors) > 0) {
  
  brms_top_features_plot <- brms_confirmatory_posteriors %>%
    dplyr::group_by(Model, feature_model) %>%
    dplyr::summarise(
      max_abs_estimate = max(abs(Estimate), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(Model) %>%
    dplyr::arrange(dplyr::desc(max_abs_estimate), .by_group = TRUE) %>%
    dplyr::slice_head(n = 25) %>%
    dplyr::ungroup()
  
  brms_forest_tbl <- brms_confirmatory_posteriors %>%
    dplyr::inner_join(
      brms_top_features_plot %>% dplyr::select(Model, feature_model),
      by = c("Model", "feature_model")
    ) %>%
    dplyr::mutate(
      base_label = ifelse(
        !is.na(feature_label) & feature_label != "",
        paste0(view_code, "_", feature_label),
        feature_model
      ),
      contrast_label = paste0(class, " vs NP"),
      direction_prob = pmax(prob_positive, prob_negative, na.rm = TRUE),
      signature_text = dplyr::case_when(
        class == "FD" & Estimate > 0 ~ "FD > NP",
        class == "FD" & Estimate < 0 ~ "NP > FD",
        class == "SO" & Estimate > 0 ~ "SO > NP",
        class == "SO" & Estimate < 0 ~ "NP > SO",
        TRUE ~ "incierto"
      ),
      evidence_text = dplyr::case_when(
        ci_excludes_zero ~ "IC95 excluye 0",
        TRUE ~ "IC95 cruza 0"
      ),
      label_prob = paste0(
        signature_text,
        "\nPr=", sprintf("%.2f", direction_prob)
      ),
      feature_plot = paste0(
        base_label,
        " | ",
        signature_text,
        " | ",
        evidence_text
      ),
      feature_plot = factor(feature_plot, levels = rev(unique(feature_plot)))
    )
  
  p_brms_forest_signature <- brms_forest_tbl %>%
    ggplot(aes(x = feature_plot, y = Estimate)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_pointrange(
      aes(
        ymin = `Q2.5`,
        ymax = `Q97.5`,
        shape = ci_excludes_zero
      ),
      size = 0.35
    ) +
    geom_text(
      aes(label = label_prob),
      hjust = ifelse(brms_forest_tbl$Estimate >= 0, -0.05, 1.05),
      size = 2.4,
      show.legend = FALSE
    ) +
    coord_flip() +
    facet_grid(class ~ Model, scales = "free_y") +
    labs(
      title = "BRMS confirmatorio: forest plot e interpretación de firma",
      subtitle = "Referencia = NP. Positivo favorece la clase del facet frente a NP; negativo favorece NP.",
      x = "Biomarcador | dirección | evidencia",
      y = "Estimación posterior",
      shape = "IC 95% excluye 0"
    ) +
    theme_minimal()
  
  print(p_brms_forest_signature)
  
  .save_panel_plot(
    p_brms_forest_signature,
    "16_brms_confirmatory_forestplot_firma.png",
    15,
    10
  )
  
  .add_panel_manifest(
    "16_brms_confirmatory_forestplot_firma.png",
    "Forest plot BRMS con dirección, probabilidad posterior e interpretación por clase."
  )
}

## -----------------------------
## 12.15B Forest plot del contraste posterior SO - FD
## -----------------------------

if (exists("brms_confirmatory_contrasts") &&
    nrow(brms_confirmatory_contrasts) > 0) {
  
  contrast_top_features <- brms_confirmatory_contrasts %>%
    dplyr::group_by(Model, feature_model) %>%
    dplyr::summarise(
      max_abs_estimate = max(abs(Estimate), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(Model) %>%
    dplyr::arrange(dplyr::desc(max_abs_estimate), .by_group = TRUE) %>%
    dplyr::slice_head(n = 25) %>%
    dplyr::ungroup()
  
  contrast_plot_tbl <- brms_confirmatory_contrasts %>%
    dplyr::inner_join(
      contrast_top_features %>% dplyr::select(Model, feature_model),
      by = c("Model", "feature_model")
    ) %>%
    dplyr::mutate(
      base_label = ifelse(
        !is.na(feature_label) & feature_label != "",
        paste0(view_code, "_", feature_label),
        feature_model
      ),
      direction_prob = pmax(prob_SO_gt_FD, prob_FD_gt_SO, na.rm = TRUE),
      signature_text = dplyr::case_when(
        Estimate > 0 ~ "SO > FD",
        Estimate < 0 ~ "FD > SO",
        TRUE ~ "incierto"
      ),
      evidence_text = dplyr::case_when(
        ci_excludes_zero ~ "IC95 excluye 0",
        TRUE ~ "IC95 cruza 0"
      ),
      label_prob = paste0(
        signature_text,
        "\nPr=", sprintf("%.2f", direction_prob)
      ),
      feature_plot = paste0(
        base_label,
        " | ",
        signature_text,
        " | ",
        evidence_text
      ),
      feature_plot = factor(feature_plot, levels = rev(unique(feature_plot)))
    )
  
  p_brms_contrast_so_fd <- contrast_plot_tbl %>%
    ggplot(aes(x = feature_plot, y = Estimate)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_pointrange(
      aes(
        ymin = Q2.5,
        ymax = Q97.5,
        shape = ci_excludes_zero
      ),
      size = 0.35
    ) +
    geom_text(
      aes(label = label_prob),
      hjust = ifelse(contrast_plot_tbl$Estimate >= 0, -0.05, 1.05),
      size = 2.4,
      show.legend = FALSE
    ) +
    coord_flip() +
    facet_wrap(~ Model, scales = "free_y") +
    labs(
      title = "BRMS confirmatorio: contraste posterior SO - FD",
      subtitle = "Positivo favorece SO frente a FD; negativo favorece FD frente a SO.",
      x = "Biomarcador | dirección | evidencia",
      y = "Estimación posterior SO - FD",
      shape = "IC 95% excluye 0"
    ) +
    theme_minimal()
  
  print(p_brms_contrast_so_fd)
  
  .save_panel_plot(
    p_brms_contrast_so_fd,
    "16b_brms_confirmatory_contraste_SO_FD.png",
    14,
    9
  )
  
  .add_panel_manifest(
    "16b_brms_confirmatory_contraste_SO_FD.png",
    "Forest plot del contraste posterior SO - FD."
  )
}


## -----------------------------
## 12.16 Firma final BRMS por biomarcador y contraste
## -----------------------------

if (exists("brms_confirmatory_posteriors") &&
    nrow(brms_confirmatory_posteriors) > 0) {
  
  ## Contrastes contra referencia NP: FD vs NP y SO vs NP
  signature_ref_tbl <- brms_confirmatory_posteriors %>%
    dplyr::mutate(
      contrast = paste0(class, "_vs_NP"),
      contrast_label = paste0(class, " vs NP"),
      direction_group = dplyr::case_when(
        class == "FD" & Estimate > 0 ~ "FD > NP",
        class == "FD" & Estimate < 0 ~ "NP > FD",
        class == "SO" & Estimate > 0 ~ "SO > NP",
        class == "SO" & Estimate < 0 ~ "NP > SO",
        TRUE ~ "uncertain"
      ),
      direction_prob = pmax(prob_positive, prob_negative, na.rm = TRUE)
    ) %>%
    dplyr::select(
      Panel,
      Model,
      feature_model,
      view,
      view_code,
      feature_label,
      contrast,
      contrast_label,
      Estimate,
      Q2.5,
      Q97.5,
      ci_excludes_zero,
      direction_group,
      direction_prob
    )
  
  ## Contraste directo SO vs FD
  signature_so_fd_tbl <- tibble::tibble()
  
  if (exists("brms_confirmatory_contrasts") &&
      nrow(brms_confirmatory_contrasts) > 0) {
    
    signature_so_fd_tbl <- brms_confirmatory_contrasts %>%
      dplyr::mutate(
        contrast_label = "SO vs FD",
        direction_group = dplyr::case_when(
          Estimate > 0 ~ "SO > FD",
          Estimate < 0 ~ "FD > SO",
          TRUE ~ "uncertain"
        ),
        direction_prob = pmax(prob_SO_gt_FD, prob_FD_gt_SO, na.rm = TRUE)
      ) %>%
      dplyr::select(
        Panel,
        Model,
        feature_model,
        view,
        view_code,
        feature_label,
        contrast,
        contrast_label,
        Estimate,
        Q2.5,
        Q97.5,
        ci_excludes_zero,
        direction_group,
        direction_prob
      )
  }
  
  signature_contrast_plot_tbl <- dplyr::bind_rows(
    signature_ref_tbl,
    signature_so_fd_tbl
  ) %>%
    dplyr::mutate(
      feature_plot = ifelse(
        !is.na(feature_label) & feature_label != "",
        paste0(view_code, "_", feature_label),
        feature_model
      ),
      evidence_label = dplyr::case_when(
        ci_excludes_zero ~ "*",
        TRUE ~ ""
      ),
      label_tile = paste0(
        direction_group,
        "\nβ=", sprintf("%.2f", Estimate),
        "\nPr=", sprintf("%.2f", direction_prob),
        evidence_label
      )
    ) %>%
    dplyr::group_by(Model, feature_plot) %>%
    dplyr::mutate(max_abs_feature = max(abs(Estimate), na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Model, dplyr::desc(max_abs_feature), feature_plot) %>%
    dplyr::mutate(
      feature_plot = factor(feature_plot, levels = rev(unique(feature_plot))),
      contrast_label = factor(
        contrast_label,
        levels = c("FD vs NP", "SO vs NP", "SO vs FD")
      )
    )
  
  p_brms_signature <- signature_contrast_plot_tbl %>%
    ggplot(aes(x = contrast_label, y = feature_plot, fill = direction_group)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label_tile), size = 2.4, color = "black") +
    facet_wrap(~ Model, scales = "free_y") +
    labs(
      title = "BRMS confirmatorio: firma final por biomarcador",
      subtitle = "NP es la referencia. * indica IC95% que excluye 0. Se añade contraste posterior SO vs FD.",
      x = "Contraste",
      y = "Biomarcador",
      fill = "Dirección"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  print(p_brms_signature)
  
  .save_panel_plot(
    p_brms_signature,
    "17_brms_confirmatory_firma_final_por_biomarcador.png",
    13,
    10
  )
  
  .add_panel_manifest(
    "17_brms_confirmatory_firma_final_por_biomarcador.png",
    "Heatmap de firma final por biomarcador usando contrastes FD vs NP, SO vs NP y SO vs FD."
  )
}

table_plot_map <- tibble::tribble(
  ~table, ~main_plot,
  
  "00_audit_cv_folds.csv", "01_auditoria_folds_cv.png",
  "00_audit_splsda_seeds.csv", "02_auditoria_exito_modelos.png",
  "00_audit_LASSO_seeds.csv", "02_auditoria_exito_modelos.png",
  "00_audit_plsda_folds.csv", "02_auditoria_exito_modelos.png",
  
  "01_mofa_features_input.csv", "09_mofa_percentil_vs_consenso_supervisado.png",
  "02_splsda_feature_stability.csv", "05_frecuencia_seleccion_top40_heatmap.png",
  "03_LASSO_feature_stability.csv", "05_frecuencia_seleccion_top40_heatmap.png",
  "05_plsda_feature_stability.csv", "05_frecuencia_seleccion_top40_heatmap.png",
  
  "06_panel_feature_consensus_all.csv", "07_consenso_binario_top50_heatmap.png",
  "07_PANEL_FINAL_mofa_2of3.csv", "11_top_panel_final_consenso.png",
  
  "09_cv_metrics_all_models_long.csv", "04_cv_metricas_distribucion_por_modelo.png",
  "10_cv_metrics_summary_by_model.csv", "03_cv_metricas_promedio_por_modelo.png",
  
  "11_test_panel_metrics_audit.csv", "10_test_panel_metricas_auditoria.png",
  "11b_test_panel_class_metrics.csv", "10b_test_panel_metricas_por_clase.png",
  
  "12_brms_confirmatory_metrics.csv", "13_brms_confirmatory_metricas_0_1.png",
  "13_brms_confirmatory_posteriors.csv", "16_brms_confirmatory_forestplot_firma.png",
  "13b_brms_confirmatory_signature.csv", "17_brms_confirmatory_firma_final_por_biomarcador.png",
  "13c_brms_confirmatory_contrasts_SO_FD.csv", "16b_brms_confirmatory_contraste_SO_FD.png",
  
  "14_brms_confirmatory_predictions.csv", "13_brms_confirmatory_metricas_0_1.png",
  "15_brms_confirmatory_confusion.csv", "15_brms_confirmatory_metricas_por_clase.png",
  "15b_brms_confirmatory_class_metrics.csv", "15_brms_confirmatory_metricas_por_clase.png"
)
write.csv(
  table_plot_map,
  out_table("00_table_plot_map.csv"),
  row.names = FALSE
)
write.csv(
  panel_plot_manifest,
  out_table("13_manifest_panel_metric_plots.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("PLOTS PANEL FINAL MULTIMODELO TERMINADOS\n")
cat("Plots guardados en:", METRIC_PANEL_PLOTS_DIR, "\n")
cat("Manifest guardado en:", out_table("13_manifest_panel_metric_plots.csv"), "\n")
cat("============================================================\n")

panel_config <- list(
  analysis_name = ANALYSIS_NAME,
  outdir_base = OUTDIR_BASE,
  analysis_dir = ANALYSIS_DIR,
  model_source = MODEL_SOURCE,
  model_dir_name = MODEL_DIR_NAME,
  panel_dir_name = PANEL_DIR_NAME,
  input_modelamiento_dir = INPUT_MODELAMIENTO_DIR,
  input_rds_dir = INPUT_RDS_DIR,
  input_tables_dir = INPUT_TABLES_DIR,
  panel_dir = PANEL_DIR,
  tables_dir = TABLES_DIR,
  rds_dir = RDS_DIR,
  plots_dir = PLOTS_DIR,
  ref_class = REF_CLASS,
  freq_thr_model = FREQ_THR_MODEL,
  mofa_percentile_thr = MOFA_PERCENTILE_THR,
  top_prop_plsda = TOP_PROP_DENSE_MODELS,
  ncomp_splsda = NCOMP_SPLSDA,
  ncomp_plsda = NCOMP_PLSDA,
  keepx_grid = KEEPX_GRID,
  lambda_grid = LAMBDA_GRID,
  alpha_grid_lasso = ALPHA_GRID_LASSO,
  do_brms_confirmatory = DO_BRMS_CONFIRMATORY,
  brms_confirm_max_features = BRMS_CONFIRM_MAX_FEATURES,
  brms_confirm_min_features = BRMS_CONFIRM_MIN_FEATURES,
  brms_confirm_iter = BRMS_CONFIRM_ITER,
  brms_confirm_chains = BRMS_CONFIRM_CHAINS,
  brms_confirm_warmup = BRMS_CONFIRM_WARMUP,
  brms_confirm_prior_sd_b = BRMS_CONFIRM_PRIOR_SD_B,
  brms_confirm_prior_sd_int = BRMS_CONFIRM_PRIOR_SD_INT,
  seeds_splsda = SEEDS_SPLSDA,
  seeds_glmnet = SEEDS_GLMNET,
  k_cv = K_CV,
  n_repeats_cv = N_REPEATS_CV,
  n_cores = N_CORES,
  brms_confirm_top_n = NA_integer_,
  w_splsda = W_SPLSDA,
  w_lasso = W_LASSO,
  w_plsda = W_PLSDA
)

saveRDS(
  list(
    config = panel_config,
    fold_audit = fold_audit,
    splsda_seed_audit = splsda_seed_audit,
    lasso_audit = lasso_audit,
    plsda_audit = plsda_audit,
    cv_metrics_all = cv_metrics_all,
    cv_metrics_summary = cv_metrics_summary,
    multinom_cv_pred = multinom_cv_pred,
    multinom_cv_metrics = multinom_cv_metrics,
    panel_feature_tbl = panel_feature_tbl,
    panel_final = panel_final,
    test_panel_metrics = test_panel_metrics,
    test_panel_class_metrics = test_panel_class_metrics,
    brms_confirmatory = brms_confirmatory,
    brms_confirmatory_metrics = brms_confirmatory_metrics,
    brms_confirmatory_posteriors = brms_confirmatory_posteriors,
    brms_confirmatory_predictions = brms_confirmatory_predictions,
    brms_confirmatory_confusion = brms_confirmatory_confusion,
    brms_confirmatory_class_metrics = brms_confirmatory_class_metrics,
    brms_confirmatory_signature = brms_confirmatory_signature,
    brms_confirmatory_contrasts = brms_confirmatory_contrasts,
    mofa_pass_tbl = mofa_pass_tbl,
    splsda_feature_tbl = splsda_feature_tbl,
    lasso_feature_tbl = lasso_feature_tbl,
    plsda_feature_tbl = plsda_feature_tbl
  ),
  out_rds("panel_final_multimodel_object.rds")
)

saveRDS(panel_config, out_rds("panel_config.rds"))
cat("\n============================================================\n")
cat("SCRIPT TERMINADO\n")
cat("Resultados guardados en:\n")
cat("Panel final:", PANEL_DIR, "\n")
cat("Tablas     :", TABLES_DIR, "\n")
cat("RDS        :", RDS_DIR, "\n")
cat("Plots      :", PLOTS_DIR, "\n")
cat("Archivo principal panel final:\n")
cat(out_table("07_PANEL_FINAL_mofa_2of3.csv"), "\n")
cat("Objeto RDS principal:\n")
cat(out_rds("panel_final_multimodel_object.rds"), "\n")
cat("BRMS confirmatorio:\n")
cat(out_rds("brms_confirmatory_object.rds"), "\n")
cat(out_table("12_brms_confirmatory_metrics.csv"), "\n")
cat(out_table("13_brms_confirmatory_posteriors.csv"), "\n")
cat(out_table("13b_brms_confirmatory_signature.csv"), "\n")
cat(out_table("13c_brms_confirmatory_contrasts_SO_FD.csv"), "\n")
cat(out_table("15b_brms_confirmatory_class_metrics.csv"), "\n")

cat("============================================================\n")
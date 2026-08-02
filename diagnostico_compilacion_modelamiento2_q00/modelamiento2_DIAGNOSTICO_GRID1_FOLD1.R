

rm(list = ls())
# gc()

##===librerias =====
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(glmnet); library(nnet); library(MASS)
  library(caret); library(pROC); library(FNN)
  library(ggplot2); library(ggridges)
  library(MOFA2); library(AnnotationDbi);library(brms)
})

set.seed(123)

## ============================================================
## CONFIGURACIÓN mRMR / REDUNDANCIA
## ------------------------------------------------------------
## No cambia el score MOFA.
## Solo evita seleccionar features casi duplicadas dentro de cada vista.
## ============================================================

MRMR_RHO_CUTOFF <- 0.95
MRMR_COR_METHOD <- "spearman"   # más robusto que pearson con n pequeño

## ============================================================
## CONFIGURACIÓN PRIORS BRMS
## ------------------------------------------------------------
## Interceptos: centrados en prevalencia suavizada de TRAIN.
## Coeficientes: prior común débilmente regularizante.
## No se usa MOFA ni vista para definir el prior de coeficientes.
## ============================================================

PRIOR_PREVALENCE_ALPHA <- 0.5
PRIOR_INTERCEPT_SD     <- 1.0
PRIOR_COEF_SD          <- 0.5
## ============================================================
## CONFIGURACIÓN PARAMETRIZABLE DEL ANÁLISIS
## Debe usar el mismo ANALYSIS_NAME que preprocess_data.R e INTEGRACION.R
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

ANALYSIS_DIR <- file.path(OUTDIR_BASE, ANALYSIS_NAME)

RFM_FILE <- file.path(
  ANALYSIS_DIR,
  paste0("ready_for_modeling_", ANALYSIS_NAME, ".rds")
)

MOFA_MODEL_FILE <- file.path(
  ANALYSIS_DIR,
  "integracion",
  "rds",
  paste0("modelo_mofa_final_", ANALYSIS_NAME, ".rds")
)
PREPROCESS_CONFIG_FILE <- file.path(
  ANALYSIS_DIR,
  "run_config_preprocess.rds"
)
OUTDIR <- file.path(
  ANALYSIS_DIR,
  "modelamiento"
)

PLOTS_DIR <- file.path(
  OUTDIR,
  "plots"
)

TABLES_DIR <- file.path(
  OUTDIR,
  "tables"
)

RDS_DIR <- file.path(
  OUTDIR,
  "rds"
)

dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTDIR,       recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR,      recursive = TRUE, showWarnings = FALSE)

out_file <- function(...) {
  file.path(OUTDIR, ...)
}

out_plot <- function(...) {
  file.path(PLOTS_DIR, ...)
}

out_table <- function(...) {
  file.path(TABLES_DIR, ...)
}

out_rds <- function(...) {
  file.path(RDS_DIR, ...)
}

cat("\n============================================================\n")
cat("ANÁLISIS:", ANALYSIS_NAME, "\n")
cat("OUTDIR_BASE:", OUTDIR_BASE, "\n")
cat("ANALYSIS_DIR:", ANALYSIS_DIR, "\n")
cat("RFM_FILE:", RFM_FILE, "\n")
cat("MOFA_MODEL_FILE:", MOFA_MODEL_FILE, "\n")
cat("PREPROCESS_CONFIG_FILE:", PREPROCESS_CONFIG_FILE, "\n")
cat("OUTDIR modelamiento:", OUTDIR, "\n")
cat("PLOTS_DIR:", PLOTS_DIR, "\n")
cat("TABLES_DIR:", TABLES_DIR, "\n")
cat("RDS_DIR:", RDS_DIR, "\n")
cat("============================================================\n")
##===Funciones =====

refactor_levels <- function(x) factor(as.character(x))

sanitize_names <- function(v){
  v2 <- gsub("[^A-Za-z0-9]+", "_", v)
  v2 <- gsub("__+", "_", v2)
  v2 <- gsub("^_+|_+$", "", v2)
  v2 <- ifelse(grepl("^[0-9]", v2), paste0("X", v2), v2)
  v2 <- gsub("_+$", "", v2)
  v2 <- make.names(v2, unique = FALSE)
  v2 <- gsub("\\.+", "_", v2)
  v2 <- gsub("^_+|_+$", "", v2)
  v2 <- gsub("__+", "_", v2)
  make.unique(v2, sep = "_")
}

std_and_prune <- function(df){
  y <- df$y
  X <- df[, setdiff(colnames(df),"y"), drop=FALSE]
  nzv <- caret::nearZeroVar(X)
  if(length(nzv)) X <- X[, -nzv, drop=FALSE]
  if(ncol(X) > 2){
    cc <- stats::cor(X, use="pairwise.complete.obs")
    bad <- suppressWarnings(caret::findCorrelation(cc, cutoff=0.95, names=FALSE))
    if(length(bad)) X <- X[, -bad, drop=FALSE]
  }
  mu <- vapply(X, mean, 0.0, na.rm=TRUE)
  sd <- vapply(X,   sd, 0.0, na.rm=TRUE)
  sd[sd == 0] <- 1
  Xs <- sweep(sweep(X, 2, mu, "-"), 2, sd, "/")
  data.frame(y=y, Xs, check.names=FALSE)
}

select_features_view <- function(mod,
                                 view,
                                 data_list,
                                 metadata = NULL,
                                 q_thr = NULL,
                                 n_keep = NULL,
                                 factors = 1:2,
                                 scale_weights = FALSE) {
  
  aux <- MOFA2::plot_top_weights(
    mod,
    view      = view,
    factors   = factors,
    abs       = TRUE,
    scale     = scale_weights,
    nfeatures = Inf
  )
  
  weight_tbl <- aux$data %>%
    as.data.frame() %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(
      importance = max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(importance))
  
  if (!is.null(n_keep)) {
    
    n_keep <- as.integer(n_keep)
    n_keep <- max(1L, min(n_keep, nrow(weight_tbl)))
    
    feats <- weight_tbl %>%
      dplyr::slice_head(n = n_keep) %>%
      dplyr::pull(feature)
    
  } else {
    
    if (is.null(q_thr)) {
      stop("Debes dar q_thr o n_keep.")
    }
    
    thr <- stats::quantile(
      weight_tbl$importance,
      probs = q_thr,
      na.rm = TRUE
    )
    
    feats <- weight_tbl %>%
      dplyr::filter(importance >= thr) %>%
      dplyr::pull(feature)
  }
  
  feats <- unique(as.character(feats))
  
  sub_train <- data_list$train[feats, , drop = FALSE]
  sub_test  <- data_list$test[feats,  , drop = FALSE]
  
  if (!is.null(metadata) && all(c("EntrezGeneID", "GeneSymbol") %in% colnames(metadata))) {
    
    map_dic <- metadata %>%
      dplyr::select(EntrezGeneID, GeneSymbol) %>%
      dplyr::distinct()
    
    idx <- rownames(sub_train)
    
    ids_map <- data.frame(EntrezGeneID = idx) %>%
      dplyr::left_join(map_dic, by = "EntrezGeneID")
    
    new_ids <- ifelse(
      !is.na(ids_map$GeneSymbol),
      ids_map$GeneSymbol,
      ids_map$EntrezGeneID
    )
    
    rownames(sub_train) <- new_ids
    rownames(sub_test)  <- new_ids
  }
  
  list(
    train = sub_train,
    test  = sub_test,
    selected_features = feats,
    weight_table = weight_tbl
  )
}


project_all_views <- function(X_list, W_list){
  W_all <- NULL; X_all <- NULL
  for(v in names(W_list)){
    Xv <- X_list[[v]]; if(is.list(Xv) && "test" %in% names(Xv)) Xv <- Xv$test
    common_feats <- intersect(rownames(W_list[[v]]), rownames(Xv))
    if(length(common_feats) == 0) next
    W_sub <- as.matrix(W_list[[v]][common_feats, , drop=FALSE])
    X_sub <- as.matrix(Xv[common_feats, , drop=FALSE])
    W_all <- rbind(W_all, W_sub); X_all <- rbind(X_all, X_sub)
  }
  stopifnot(nrow(W_all) == nrow(X_all))
  Z_est <- solve(t(W_all) %*% W_all) %*% t(W_all) %*% X_all
  t(Z_est)
}

build_sets_from_grid <- function(thr_grid,
                                 tst_Data,
                                 mod,
                                 view_map,
                                 metadata_by_code = NULL,
                                 factors = 1:2,
                                 scale_weights = FALSE) {
  
  out_train <- vector("list", nrow(thr_grid))
  out_test  <- vector("list", nrow(thr_grid))
  out_info  <- vector("list", nrow(thr_grid))
  
  y_train <- factor(tst_Data$grupo_train$grupo)
  names(y_train) <- tst_Data$grupo_train$id
  
  y_test <- factor(tst_Data$grupo_test$grupo, levels = levels(y_train))
  names(y_test) <- tst_Data$grupo_test$id
  
  view_names <- names(view_map)
  view_codes <- unname(view_map)
  
  for (i in seq_len(nrow(thr_grid))) {
    
    train_list_i <- list()
    test_list_i  <- list()
    info_i       <- list()
    
    for (j in seq_along(view_map)) {
      
      view_name <- view_names[j]
      view_code <- view_codes[j]
      
      n_col <- paste0("n_", view_code)
      q_col <- paste0("q_", view_code)
      
      n_keep_i <- NULL
      q_thr_i  <- NULL
      
      if (n_col %in% colnames(thr_grid)) {
        n_keep_i <- thr_grid[[n_col]][i]
        
        if (!is.null(n_keep_i) && (is.na(n_keep_i) || !is.finite(n_keep_i))) {
          stop(
            "n_keep_i inválido para GridID=", i,
            ", vista=", view_name,
            ", columna=", n_col,
            ", valor=", n_keep_i
          )
        }
        
      } else if (q_col %in% colnames(thr_grid)) {
        q_thr_i <- thr_grid[[q_col]][i]
      } else {
        stop(
          "No encuentro columna ", n_col, " ni ", q_col,
          " para la vista ", view_name, " / código ", view_code
        )
      }
      
      metadata_i <- NULL
      
      if (!is.null(metadata_by_code) && view_code %in% names(metadata_by_code)) {
        metadata_i <- metadata_by_code[[view_code]]
      }
      
      sel <- select_features_view(
        mod           = mod,
        view          = j,
        data_list     = tst_Data[[view_code]],
        metadata      = metadata_i,
        q_thr         = q_thr_i,
        n_keep        = n_keep_i,
        factors       = factors,
        scale_weights = scale_weights
      )
      
      rownames(sel$train) <- paste0(view_code, "__", rownames(sel$train))
      rownames(sel$test)  <- paste0(view_code, "__", rownames(sel$test))
      
      train_list_i[[view_code]] <- sel$train
      test_list_i[[view_code]]  <- sel$test
      
      info_i[[view_code]] <- data.frame(
        view_name = view_name,
        view_code = view_code,
        n_selected = nrow(sel$train),
        stringsAsFactors = FALSE
      )
    }
    
    train_feats <- do.call(rbind, train_list_i)
    test_feats  <- do.call(rbind, test_list_i)
    
    train_feats <- train_feats[, tst_Data$grupo_train$id, drop = FALSE]
    test_feats  <- test_feats[,  tst_Data$grupo_test$id,  drop = FALSE]
    
    out_train[[i]] <- as.data.frame(t(train_feats), check.names = FALSE)
    out_test[[i]]  <- as.data.frame(t(test_feats),  check.names = FALSE)
    
    out_info[[i]] <- dplyr::bind_rows(info_i) %>%
      dplyr::mutate(GridID = i)
  }
  
  list(
    train   = out_train,
    test    = out_test,
    y_train = y_train[tst_Data$grupo_train$id],
    y_test  = y_test[tst_Data$grupo_test$id],
    info    = dplyr::bind_rows(out_info)
  )
}

metrics_from_probs <- function(prob_mat, y_true){
  lev <- levels(y_true); prob_mat <- prob_mat[, lev, drop=FALSE]
  pred_labels <- apply(prob_mat, 1, function(p) lev[which.max(p)])
  pred_class  <- factor(pred_labels, levels=lev)
  cm <- caret::confusionMatrix(pred_class, y_true)
  auc_obj <- tryCatch(pROC::multiclass.roc(y_true, prob_mat), error=function(e) NULL)
  auc_num <- if(is.null(auc_obj)) NA_real_ else as.numeric(auc_obj$auc)
  Y <- model.matrix(~ y_true - 1); colnames(Y) <- lev
  eps <- 1e-15
  logloss <- -mean(rowSums(Y * log(pmax(pmin(prob_mat, 1-eps), eps))))
  tibble(Accuracy=as.numeric(cm$overall["Accuracy"]),
         BalAcc  =mean(cm$byClass[,"Balanced Accuracy"], na.rm=TRUE),
         Kappa   =as.numeric(cm$overall["Kappa"]),
         AUC=auc_num, LogLoss=logloss)
}

.get_global_or_null <- function(name) {
  if (exists(name, inherits = TRUE)) {
    get(name, inherits = TRUE)
  } else {
    NULL
  }
}

.infer_feature_view_tbl <- function(feature_names) {
  
  tibble::tibble(
    feature_model = feature_names,
    view_code = sub("_.*$", "", feature_names),
    view = dplyr::case_when(
      view_code == "tx" ~ "transcriptomica",
      view_code == "pr" ~ "proteomica",
      view_code == "me" ~ "metabolomica",
      view_code == "cl" ~ "clinical",
      TRUE ~ view_code
    )
  )
}

.make_class_prevalence_tbl <- function(y,
                                       y_levels,
                                       ref = "NP",
                                       alpha = PRIOR_PREVALENCE_ALPHA) {
  
  y_levels <- as.character(y_levels)
  
  if (!ref %in% y_levels) {
    stop("La clase de referencia ref = ", ref, " no está en y_levels.")
  }
  
  if (is.null(y)) {
    n_by_class <- rep(1L, length(y_levels))
    names(n_by_class) <- y_levels
  } else {
    y <- factor(y, levels = y_levels)
    n_by_class <- table(y)
  }
  
  n_total <- sum(n_by_class)
  K <- length(y_levels)
  
  pi_smooth <- (as.numeric(n_by_class) + alpha) /
    (n_total + K * alpha)
  
  names(pi_smooth) <- y_levels
  
  pi_ref <- pi_smooth[[ref]]
  
  tibble::tibble(
    class = y_levels,
    n_class = as.integer(n_by_class[y_levels]),
    prevalence_smooth = as.numeric(pi_smooth[y_levels]),
    ref_class = ref,
    prevalence_ref_smooth = as.numeric(pi_ref),
    intercept_prior_mean = log(prevalence_smooth / prevalence_ref_smooth)
  )
}

make_multinom_prior_audit_tbl <- function(y = NULL,
                                          y_levels = NULL,
                                          feature_names = NULL,
                                          ref = "NP",
                                          view_map = NULL,
                                          r2_score = NULL,
                                          prevalence_alpha = PRIOR_PREVALENCE_ALPHA,
                                          intercept_sd = PRIOR_INTERCEPT_SD,
                                          coef_sd = PRIOR_COEF_SD,
                                          ...) {
  
  if (is.null(y_levels)) {
    if (!is.null(y)) {
      y_levels <- levels(factor(y))
    } else {
      stop("Debes proporcionar y o y_levels.")
    }
  }
  
  y_levels <- as.character(y_levels)
  levs <- setdiff(y_levels, ref)
  
  prev_tbl <- .make_class_prevalence_tbl(
    y = y,
    y_levels = y_levels,
    ref = ref,
    alpha = prevalence_alpha
  )
  
  intercept_audit <- prev_tbl %>%
    dplyr::filter(class %in% levs) %>%
    dplyr::mutate(
      prior_type = "Intercept",
      dpar = paste0("mu", class),
      prior_family = "normal",
      prior_mean = intercept_prior_mean,
      prior_sd = intercept_sd,
      feature_model = NA_character_,
      view = NA_character_,
      view_code = NA_character_,
      prior_source = "train_prevalence_smoothed"
    ) %>%
    dplyr::select(
      prior_type,
      dpar,
      class,
      feature_model,
      view,
      view_code,
      n_class,
      prevalence_smooth,
      ref_class,
      prevalence_ref_smooth,
      prior_family,
      prior_mean,
      prior_sd,
      prior_source
    )
  
  if (is.null(feature_names)) {
    return(intercept_audit)
  }
  
  feature_view_tbl <- .infer_feature_view_tbl(feature_names)
  
  coef_audit <- tidyr::expand_grid(
    class = levs,
    feature_view_tbl
  ) %>%
    dplyr::mutate(
      prior_type = "b",
      dpar = paste0("mu", class),
      n_class = NA_integer_,
      prevalence_smooth = NA_real_,
      ref_class = ref,
      prevalence_ref_smooth = NA_real_,
      prior_family = "normal",
      prior_mean = 0,
      prior_sd = coef_sd,
      prior_source = "common_weakly_regularizing"
    ) %>%
    dplyr::select(
      prior_type,
      dpar,
      class,
      feature_model,
      view,
      view_code,
      n_class,
      prevalence_smooth,
      ref_class,
      prevalence_ref_smooth,
      prior_family,
      prior_mean,
      prior_sd,
      prior_source
    )
  
  dplyr::bind_rows(intercept_audit, coef_audit)
}

make_multinom_priors <- function(y_levels = NULL,
                                 ref = "NP",
                                 sd_b = NULL,
                                 sd_int = NULL,
                                 y = NULL,
                                 feature_names = NULL,
                                 view_map = NULL,
                                 r2_score = NULL,
                                 prevalence_alpha = PRIOR_PREVALENCE_ALPHA,
                                 intercept_sd = PRIOR_INTERCEPT_SD,
                                 coef_sd = PRIOR_COEF_SD,
                                 ...) {
  
  ## Compatibilidad con llamadas antiguas
  if (!is.null(sd_b)) {
    coef_sd <- sd_b
  }
  
  if (!is.null(sd_int)) {
    intercept_sd <- sd_int
  }
  
  if (is.null(y_levels)) {
    if (!is.null(y)) {
      y_levels <- levels(factor(y))
    } else {
      stop("Debes proporcionar y_levels o y.")
    }
  }
  
  y_levels <- as.character(y_levels)
  levs <- setdiff(y_levels, ref)
  
  prior_audit <- make_multinom_prior_audit_tbl(
    y = y,
    y_levels = y_levels,
    feature_names = feature_names,
    ref = ref,
    prevalence_alpha = prevalence_alpha,
    intercept_sd = intercept_sd,
    coef_sd = coef_sd
  )
  
  prior_list <- list()
  
  for (lev in levs) {
    
    dpar <- paste0("mu", lev)
    
    int_row <- prior_audit %>%
      dplyr::filter(
        prior_type == "Intercept",
        class == lev
      ) %>%
      dplyr::slice(1)
    
    prior_list[[length(prior_list) + 1L]] <- brms::set_prior(
      paste0(
        "normal(",
        signif(int_row$prior_mean, 6),
        ",",
        signif(int_row$prior_sd, 6),
        ")"
      ),
      class = "Intercept",
      dpar = dpar
    )
    
    prior_list[[length(prior_list) + 1L]] <- brms::set_prior(
      paste0("normal(0,", signif(coef_sd, 6), ")"),
      class = "b",
      dpar = dpar
    )
  }
  
  do.call(c, prior_list)
}

brms_probs <- function(fit, newdata){
  pp <- posterior_epred(fit, newdata=newdata)
  M  <- apply(pp, c(2,3), mean)
  cls <- dimnames(pp)[[3]]; if(is.null(cls)) cls <- fit$family$names
  colnames(M) <- cls; M
}

count_coefs_nnet <- function(fit){
  cf <- coef(fit); if(is.list(cf)) cf <- do.call(rbind, cf)
  cf <- as.matrix(cf)
  if(ncol(cf) > 0) cf_noint <- cf[, -1, drop=FALSE] else cf_noint <- cf
  sum(abs(cf_noint) > 0)
}

count_coefs_brms <- function(fit){
  fe <- tryCatch(as.data.frame(brms::fixef(fit)), error=function(e) NULL)
  if(is.null(fe)) return(NA_integer_)
  keep <- !grepl("Intercept", rownames(fe))
  sum(abs(fe$Estimate[keep]) > 0)
}

extract_biomarkers_brms <- function(fit, top_k=50){
  fx <- as.data.frame(brms::fixef(fit, robust=TRUE)); fx$term <- rownames(fx)
  fx <- fx[!grepl("Intercept", fx$term), ]
  fx$feature <- sub("^b_", "", sub("^.*?b_", "b_", fx$term))
  fx$class   <- sub("^mu([^_]+).*", "\\1", fx$term)
  fx$nonzero <- (fx$Q2.5 > 0) | (fx$Q97.5 < 0)
  agg <- fx %>% group_by(feature) %>% summarise(any_nonzero=any(nonzero),
                                                max_abs_est=max(abs(Estimate), na.rm=TRUE),
                                                n_classes=sum(nonzero), .groups="drop") %>%
    arrange(desc(any_nonzero), desc(n_classes), desc(max_abs_est))
  list(per_class = fx[fx$nonzero, c("class","feature","Estimate","Q2.5","Q97.5")],
       ranking   = head(agg, top_k))
}


evaluate_all <- function(return.list, thr_grid, ref="NP",
                         do_brms=FALSE, brms_iter=2000, brms_chains=2,
                         brms_algorithm=c("mcmc","vb")){
  options(contrasts = c("contr.treatment","contr.poly"))
  brms_algorithm <- match.arg(brms_algorithm)
  res_all <- list()
  
  detect_backend <- function(){
    if (requireNamespace("cmdstanr", quietly=TRUE)) {
      v <- try(cmdstanr::cmdstan_version(), silent=TRUE)
      if (!inherits(v, "try-error")) return("cmdstanr")
    }
    "rstan"
  }
  backend_choice <- detect_backend()
  brms_errs <- list()
  
  rank_one <- function(df){
    df %>% mutate(r1 = rank(-BalAcc, ties.method="min"),
                  r2 = rank(LogLoss,  ties.method="min"),
                  r3 = rank(-AUC,    ties.method="min"),
                  RankSum = r1 + r2 + r3) %>%
      arrange(RankSum, desc(BalAcc), AUC, LogLoss)
  }
  
  # construye control según backend
  make_ctrl <- function(backend){
    if(identical(backend, "cmdstanr")){
      list(adapt_delta=0.95, max_treedepth=12)                 # sin init_r
    } else {
      list(adapt_delta=0.95, max_treedepth=12, init_r=0.1)     # rstan permite init_r
    }
  }
  
  fit_brms_with_fallback <- function(df_train, ref_lbl, algo=brms_algorithm,
                                     iter=brms_iter, chains=brms_chains, backend=backend_choice){
    priors <- make_multinom_priors(
      y_levels = levels(df_train$y),
      y = df_train$y,
      feature_names = setdiff(colnames(df_train), "y"),
      ref = ref_lbl,
      view_map = .get_global_or_null("view_map"),
      r2_score = .get_global_or_null("r2_score")
    )
    
    family_cat <- brms::categorical()
    vb_err <- NULL; mcmc_err <- NULL
    
    if(algo %in% c("vb","meanfield","fullrank")){
      fit_vb <- tryCatch(
        brm(y ~ ., data=df_train, family=family_cat, prior=priors,
            algorithm="fullrank", iter=max(4000, iter), refresh=0, backend=backend,
            control=list(tol_rel_obj=1e-4, eval_elbo=200, adapt_engaged=TRUE)),
        error=function(e){ vb_err <<- conditionMessage(e); NULL }
      )
      if(!is.null(fit_vb)){
        kk <- tryCatch({ ll <- brms::loo(fit_vb); max(ll$diagnostics$pareto_k) }, error=function(e) Inf)
        if(is.finite(kk) && kk <= 0.7) return(list(fit=fit_vb, label="BRMS_VB", err=NULL))
        vb_err <- sprintf("VB Pareto-k=%.2f", kk)
      }
      if(!is.null(vb_err)) message("[VB] ", vb_err)
    }
    
    fit_mcmc <- tryCatch(
      brm(y ~ ., data=df_train, family=family_cat, prior=priors,
          chains=max(2, chains), iter=max(2000, iter), warmup=floor(max(2000, iter)/2),
          refresh=0, backend=backend, control=make_ctrl(backend)),
      error=function(e){ mcmc_err <<- conditionMessage(e); NULL }
    )
    if(!is.null(fit_mcmc)) return(list(fit=fit_mcmc, label="BRMS_MCMC", err=NULL))
    
    errs <- c(vb_err, mcmc_err); errs <- errs[!vapply(errs, is.null, logical(1))]
    if(length(errs)==0) errs <- NA_character_
    list(fit=NULL, label=NA_character_, err=paste(unlist(errs), collapse=" | "))
  }
  
  for(i in seq_along(return.list$train)){
    X_train <- return.list$train[[i]]
    X_test  <- return.list$test[[i]]
    y_train <- stats::relevel(refactor_levels(return.list$y_train), ref=ref)
    y_test  <- stats::relevel(refactor_levels(return.list$y_test),  ref=ref)
    
    df_train <- data.frame(y=y_train, X_train, check.names=FALSE)
    df_test  <- data.frame(y=y_test,  X_test,  check.names=FALSE)
    
    df_train$y <- droplevels(df_train$y)
    df_test$y  <- factor(df_test$y, levels=levels(df_train$y))
    
    common_cols <- setdiff(intersect(colnames(df_train), colnames(df_test)), "y")
    new_names <- sanitize_names(common_cols); new_names <- make.unique(new_names, sep="_")
    colnames(df_train)[match(common_cols, colnames(df_train))] <- new_names
    colnames(df_test)[match(common_cols,  colnames(df_test))]  <- new_names
    
    df_train <- df_train[, c("y", new_names), drop=FALSE]
    df_test  <- df_test[,  c("y", new_names), drop=FALSE]
    
    p <- ncol(df_train) - 1L
    rows <- list()
    
    fit_nnet <- nnet::multinom(y ~ ., data=df_train, trace=FALSE, MaxNWts=200000)
    probs_nnet <- stats::predict(fit_nnet, newdata=df_test, type="probs")
    if(is.null(colnames(probs_nnet))) colnames(probs_nnet) <- levels(df_test$y)
    m_nnet <- metrics_from_probs(probs_nnet, df_test$y)
    rows[["NNET"]] <- m_nnet %>% mutate(Model="NNET", GridID=i, p=p, n_coefs=count_coefs_nnet(fit_nnet))
    
    if(do_brms){
      fr <- fit_brms_with_fallback(df_train, ref_lbl=ref, algo=brms_algorithm,
                                   iter=brms_iter, chains=brms_chains, backend=backend_choice)
      if(!is.null(fr$fit)){
        P <- brms_probs(fr$fit, df_test)[, levels(df_test$y), drop=FALSE]
        m <- metrics_from_probs(P, df_test$y)
        rows[[fr$label]] <- m %>% mutate(Model=fr$label, GridID=i, p=p, n_coefs=count_coefs_brms(fr$fit))
        attr(rows[[fr$label]], "biomarkers") <- extract_biomarkers_brms(fr$fit, top_k=100)
      } else if(!is.null(fr$err)){
        brms_errs[[length(brms_errs)+1]] <- sprintf("Grid %d: %s", i, fr$err)
      }
    }
    
    res_all[[i]] <- dplyr::bind_rows(rows)
  }
  
  res_tbl <- dplyr::bind_rows(res_all) %>%
    dplyr::mutate(dplyr::across(c(Accuracy,BalAcc,Kappa,AUC,LogLoss), as.numeric)) %>%
    dplyr::left_join(thr_grid %>% mutate(GridID = dplyr::row_number()), by="GridID")
  
  rank_one <- rank_one
  best_by_grid <- res_tbl %>% group_by(GridID) %>% group_modify(~ rank_one(.x) %>% slice_head(n=1)) %>% ungroup()
  best_overall <- rank_one(res_tbl) %>% slice_head(n=1)
  
  list(results=res_tbl, best_by_grid=best_by_grid, best_overall=best_overall,
       brms_errors=unique(unlist(brms_errs)), backend=backend_choice)
}
## ============================================================
## Funciones CV interna estratificada para BRMS_MCMC
## ============================================================

.make_stratified_folds <- function(y, K = 3, seed = 123) {
  set.seed(seed)
  
  y <- factor(y)
  tab <- table(y)
  
  if (K > min(tab)) {
    stop(
      "K = ", K, " no es válido. Clase minoritaria: ",
      paste(names(tab), tab, sep = "=", collapse = ", "),
      ". Usa K <= ", min(tab)
    )
  }
  
  fold_id <- integer(length(y))
  names(fold_id) <- names(y)
  
  for (cl in levels(y)) {
    idx <- which(y == cl)
    idx <- sample(idx)
    fold_id[idx] <- rep(seq_len(K), length.out = length(idx))
  }
  
  folds <- lapply(seq_len(K), function(k) which(fold_id == k))
  names(folds) <- paste0("Fold", seq_len(K))
  
  folds
}


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


.prepare_cv_data <- function(return.list, grid_id, train_idx, valid_idx, ref = "NP") {
  
  X_all <- return.list$train[[grid_id]]
  y_all <- stats::relevel(refactor_levels(return.list$y_train), ref = ref)
  
  if (is.null(rownames(X_all))) {
    rownames(X_all) <- names(y_all)
  }
  
  X_train <- X_all[train_idx, , drop = FALSE]
  X_valid <- X_all[valid_idx, , drop = FALSE]
  
  y_train <- y_all[train_idx]
  y_valid <- y_all[valid_idx]
  
  df_train <- data.frame(y = y_train, X_train, check.names = FALSE)
  df_valid <- data.frame(y = y_valid, X_valid, check.names = FALSE)
  
  df_train$y <- droplevels(df_train$y)
  df_valid$y <- factor(df_valid$y, levels = levels(df_train$y))
  
  common_cols <- setdiff(intersect(colnames(df_train), colnames(df_valid)), "y")
  new_names <- make.unique(sanitize_names(common_cols), sep = "_")
  
  colnames(df_train)[match(common_cols, colnames(df_train))] <- new_names
  colnames(df_valid)[match(common_cols, colnames(df_valid))] <- new_names
  
  df_train <- df_train[, c("y", new_names), drop = FALSE]
  df_valid <- df_valid[, c("y", new_names), drop = FALSE]
  
  list(
    df_train = df_train,
    df_valid = df_valid
  )
}


.fit_one_brms_cv <- function(task_id,
                             tasks,
                             return.list,
                             folds,
                             ref = "NP",
                             brms_iter = 2000,
                             brms_chains = 2,
                             seed = 123) {
  
  grid_id <- tasks$GridID[task_id]
  fold_id <- tasks$Fold[task_id]
  
  valid_idx <- folds[[fold_id]]
  train_idx <- setdiff(seq_along(return.list$y_train), valid_idx)
  
  dat <- .prepare_cv_data(
    return.list = return.list,
    grid_id     = grid_id,
    train_idx   = train_idx,
    valid_idx   = valid_idx,
    ref         = ref
  )
  
  df_train <- dat$df_train
  df_valid <- dat$df_valid
  
  p <- ncol(df_train) - 1L
  
  if (length(levels(droplevels(df_train$y))) < length(levels(refactor_levels(return.list$y_train)))) {
    return(tibble::tibble(
      GridID = grid_id,
      Fold = fold_id,
      p = p,
      Model = "BRMS_MCMC",
      Accuracy = NA_real_,
      BalAcc = NA_real_,
      Kappa = NA_real_,
      AUC = NA_real_,
      LogLoss = NA_real_,
      error = "ETAPA_PREPARACION_FOLD: Fold sin todas las clases en entrenamiento"
    ))
  }
  
  backend <- .detect_backend_safe()
  ctrl <- .make_ctrl_safe(backend)
  
  priors <- make_multinom_priors(
    y_levels = levels(df_train$y),
    y = df_train$y,
    feature_names = setdiff(colnames(df_train), "y"),
    ref = ref,
    view_map = .get_global_or_null("view_map"),
    r2_score = .get_global_or_null("r2_score")
  )
  
  fit <- tryCatch(
    brms::brm(
      y ~ .,
      data    = df_train,
      family  = brms::categorical(),
      prior   = priors,
      chains  = brms_chains,
      iter    = brms_iter,
      warmup  = floor(brms_iter / 2),
      refresh = 0,
      backend = backend,
      control = ctrl,
      cores   = 1,
      seed    = seed + 1000L * grid_id + fold_id
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(tibble::tibble(
      GridID = grid_id,
      Fold = fold_id,
      p = p,
      Model = "BRMS_MCMC",
      Accuracy = NA_real_,
      BalAcc = NA_real_,
      Kappa = NA_real_,
      AUC = NA_real_,
      LogLoss = NA_real_,
      error = paste0(
        "ETAPA_AJUSTE_BRMS: ",
        conditionMessage(fit)
      )
    ))
  }
  
  P <- tryCatch(
    brms_probs(fit, df_valid)[, levels(df_valid$y), drop = FALSE],
    error = function(e) e
  )
  
  if (inherits(P, "error")) {
    return(tibble::tibble(
      GridID = grid_id,
      Fold = fold_id,
      p = p,
      Model = "BRMS_MCMC",
      Accuracy = NA_real_,
      BalAcc = NA_real_,
      Kappa = NA_real_,
      AUC = NA_real_,
      LogLoss = NA_real_,
      error = paste0(
        "ETAPA_PREDICCION_POSTERIOR_EPRED: ",
        conditionMessage(P)
      )
    ))
  }
  
  met <- metrics_from_probs(P, df_valid$y)
  
  dplyr::bind_cols(
    tibble::tibble(
      GridID = grid_id,
      Fold = fold_id,
      p = p,
      Model = "BRMS_MCMC"
    ),
    met,
    tibble::tibble(error = NA_character_)
  )
}
## ============================================================
## Regla matemática de selección de grid
## ------------------------------------------------------------
## Para cada grid g:
##
##   LogLoss_g = promedio del LogLoss en los K folds
##   SE_g      = SD(LogLoss_g,fold) / sqrt(K_g)
##
## Primero se identifica:
##
##   g* = argmin_g LogLoss_g
##
## Luego se define el conjunto 1-SE:
##
##   C_1SE = { g : LogLoss_g <= LogLoss_g* + SE_g* }
##
## Dentro de C_1SE se aplica una regla jerárquica:
##
##   1) menor número de variables p_g
##   2) menor variabilidad de LogLoss en CV repetida
##   3) mayor frecuencia de selección en CV repetida
##   4) mayor Balanced Accuracy media
##   5) mayor AUC media
##   6) menor LogLoss medio
##
## No se usa una función ponderada tipo:
##
##   J = LogLoss + lambda_1 * SD + lambda_2 * p/n
##
## Por tanto, no se introducen pesos arbitrarios.
## El conjunto TEST se usa solo como auditoría externa.
## ============================================================

evaluate_brms_cv_grid <- function(return.list,
                                  thr_grid,
                                  ref = "NP",
                                  K = 3,
                                  p_min = 35,
                                  p_max = 85,
                                  brms_iter = 2000,
                                  brms_chains = 2,
                                  cores = 2,
                                  seed = 123) {
  
  y_train <- stats::relevel(refactor_levels(return.list$y_train), ref = ref)
  
  cat("\n--- Distribución global en train ---\n")
  print(table(y_train))
  
  folds <- .make_stratified_folds(y_train, K = K, seed = seed)
  
  cat("\n--- Auditoría de folds de validación ---\n")
  fold_audit <- dplyr::bind_rows(lapply(seq_along(folds), function(k) {
    tb <- table(y_train[folds[[k]]])
    data.frame(
      Fold = k,
      Clase = names(tb),
      n = as.integer(tb)
    )
  }))
  print(fold_audit)
  
  tasks <- expand.grid(
    GridID = seq_len(nrow(thr_grid)),
    Fold   = seq_len(K)
  )
  
  if (cores > 1) {
    
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterSetRNGStream(cl, seed)
    
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(dplyr)
        library(tibble)
        library(brms)
        library(caret)
        library(pROC)
      })
      NULL
    })
    parallel::clusterExport(
      cl,
      varlist = c(
        "tasks",
        "return.list",
        "folds",
        "ref",
        "brms_iter",
        "brms_chains",
        "seed",
        
        "refactor_levels", "sanitize_names",
        "metrics_from_probs", "brms_probs",
        
        "PRIOR_PREVALENCE_ALPHA",
        "PRIOR_INTERCEPT_SD",
        "PRIOR_COEF_SD",
        
        ".get_global_or_null",
        ".infer_feature_view_tbl",
        ".make_class_prevalence_tbl",
        "make_multinom_prior_audit_tbl",
        "make_multinom_priors",
        
        ".detect_backend_safe", ".make_ctrl_safe",
        ".prepare_cv_data", ".fit_one_brms_cv"
      ),
      envir = environment()
    )
    
    cv_list <- parallel::parLapply(
      cl,
      seq_len(nrow(tasks)),
      function(ii) {
        .fit_one_brms_cv(
          task_id     = ii,
          tasks       = tasks,
          return.list = return.list,
          folds       = folds,
          ref         = ref,
          brms_iter   = brms_iter,
          brms_chains = brms_chains,
          seed        = seed
        )
      }
    )
    
  } else {
    
    cv_list <- lapply(seq_len(nrow(tasks)), function(ii) {
      .fit_one_brms_cv(
        task_id     = ii,
        tasks       = tasks,
        return.list = return.list,
        folds       = folds,
        ref         = ref,
        brms_iter   = brms_iter,
        brms_chains = brms_chains,
        seed        = seed
      )
    })
  }
  
  cv_results <- dplyr::bind_rows(cv_list) %>%
    dplyr::left_join(
      thr_grid %>% dplyr::mutate(GridID = dplyr::row_number()),
      by = "GridID"
    ) %>%
    dplyr::mutate(
      diagnostic_status = dplyr::case_when(
        !is.na(error) & nzchar(error) ~ "ERROR_CAPTURADO",
        is.na(LogLoss) ~ "LOGLOSS_NA_SIN_ERROR_CAPTURADO",
        TRUE ~ "OK"
      )
    )
  
  ## ============================================================
  ## DIAGNÓSTICO CV
  ## No cambia resultados, selección ni matemática.
  ## Guarda los resultados antes de aplicar filtros o detenerse.
  ## ============================================================
  
  debug_file_all <- out_table(
    paste0(
      "DEBUG_cv_brms_todos_los_folds_seed_",
      seed,
      ".csv"
    )
  )
  
  debug_file_failed <- out_table(
    paste0(
      "DEBUG_cv_brms_folds_fallidos_seed_",
      seed,
      ".csv"
    )
  )
  
  write.csv(
    cv_results,
    debug_file_all,
    row.names = FALSE
  )
  
  cv_failed_debug <- cv_results %>%
    dplyr::filter(
      diagnostic_status != "OK" |
        is.na(LogLoss)
    ) %>%
    dplyr::arrange(
      GridID,
      Fold
    )
  
  write.csv(
    cv_failed_debug,
    debug_file_failed,
    row.names = FALSE
  )
  
  cat("\n============================================================\n")
  cat("DIAGNÓSTICO DE LA CV BRMS\n")
  cat("============================================================\n")
  cat("Archivo con todos los folds:\n")
  cat(debug_file_all, "\n\n")
  cat("Archivo con folds fallidos:\n")
  cat(debug_file_failed, "\n")
  
  cat("\n--- Estado por GridID y Fold ---\n")
  
  print(
    cv_results %>%
      dplyr::select(
        GridID,
        Fold,
        p,
        diagnostic_status,
        Accuracy,
        BalAcc,
        AUC,
        LogLoss,
        error
      ) %>%
      dplyr::arrange(
        GridID,
        Fold
      ),
    n = Inf,
    width = Inf
  )
  
  cat("\n--- Número de folds por estado ---\n")
  
  print(
    cv_results %>%
      dplyr::count(
        diagnostic_status,
        name = "n_folds"
      ),
    n = Inf,
    width = Inf
  )
  
  cat("\n--- Errores únicos detectados ---\n")
  
  print(
    cv_results %>%
      dplyr::filter(
        !is.na(error),
        nzchar(error)
      ) %>%
      dplyr::distinct(
        diagnostic_status,
        error
      ),
    n = Inf,
    width = Inf
  )
  
  n_train_total <- length(return.list$y_train)
  
  cv_summary <- cv_results %>%
    dplyr::group_by(GridID) %>%
    dplyr::summarise(
      Model = "BRMS_MCMC",
      n_folds_ok     = sum(!is.na(LogLoss)),
      n_folds_failed = sum(is.na(LogLoss)),
      p = dplyr::first(p),
      mean_Accuracy = mean(Accuracy, na.rm = TRUE),
      mean_BalAcc   = mean(BalAcc,   na.rm = TRUE),
      mean_Kappa    = mean(Kappa,    na.rm = TRUE),
      mean_AUC      = mean(AUC,      na.rm = TRUE),
      mean_LogLoss  = mean(LogLoss,  na.rm = TRUE),
      sd_LogLoss    = sd(LogLoss,    na.rm = TRUE),
      q_tx = dplyr::first(q_tx),
      q_pr = dplyr::first(q_pr),
      q_me = dplyr::first(q_me),
      q_cl = dplyr::first(q_cl),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      se_LogLoss = sd_LogLoss / sqrt(n_folds_ok),
      perfect_cv = mean_Accuracy >= 0.999 & mean_BalAcc >= 0.999,
      p_per_sample = p / n_train_total
    ) %>%
    dplyr::filter(n_folds_ok == K) %>%
    dplyr::filter(p >= p_min, p <= p_max)
  
  if (nrow(cv_summary) == 0) {
    
    cat("\n============================================================\n")
    cat("NINGUNA GRID COMPLETÓ LOS", K, "FOLDS\n")
    cat("============================================================\n")
    
    cat("\n--- Folds con LogLoss ausente ---\n")
    
    print(
      cv_results %>%
        dplyr::filter(
          is.na(LogLoss)
        ) %>%
        dplyr::select(
          GridID,
          Fold,
          p,
          diagnostic_status,
          error
        ) %>%
        dplyr::arrange(
          GridID,
          Fold
        ),
      n = Inf,
      width = Inf
    )
    
    cat("\n--- Resumen de folds válidos por grid antes del filtro ---\n")
    
    print(
      cv_results %>%
        dplyr::group_by(GridID) %>%
        dplyr::summarise(
          p = dplyr::first(p),
          n_folds_total = dplyr::n(),
          n_folds_ok = sum(!is.na(LogLoss)),
          n_folds_failed = sum(is.na(LogLoss)),
          errores = paste(
            unique(error[!is.na(error) & nzchar(error)]),
            collapse = " | "
          ),
          .groups = "drop"
        ) %>%
        dplyr::arrange(GridID),
      n = Inf,
      width = Inf
    )
    
    stop(
      paste0(
        "No quedó ningún grid válido tras CV. ",
        "Los errores exactos fueron impresos y guardados en: ",
        debug_file_failed
      ),
      call. = FALSE
    )
  }
  
  # Evita seleccionar grids perfectas si existen alternativas no perfectas
  cv_candidates <- cv_summary %>%
    dplyr::filter(!perfect_cv)
  
  if (nrow(cv_candidates) == 0) {
    message("Todas las grids válidas tienen CV perfecta; se usará la regla 1-SE sobre todas.")
    cv_candidates <- cv_summary
  }
  ## Regla jerárquica sin pesos arbitrarios:
  ## 1) identifica el menor LogLoss medio en CV
  ## 2) define el conjunto 1-SE: mean_LogLoss <= min_LogLoss + SE_del_mejor
  ## 3) dentro de 1-SE, elige menor p
  ## 4) desempata por menor sd_LogLoss
  ## 5) desempata por mayor Balanced Accuracy
  ## 6) desempate final técnico por mayor AUC, menor LogLoss y GridID
  
  best_idx <- which.min(cv_candidates$mean_LogLoss)
  best_logloss <- cv_candidates$mean_LogLoss[best_idx]
  best_se <- cv_candidates$se_LogLoss[best_idx]
  
  if (!is.finite(best_se)) {
    best_se <- 0
  }
  
  logloss_cutoff <- best_logloss + best_se
  
  one_se_candidates <- cv_candidates %>%
    dplyr::filter(
      is.finite(mean_LogLoss),
      mean_LogLoss <= logloss_cutoff
    )
  
  if (nrow(one_se_candidates) == 0) {
    warning("No hubo grids dentro de 1-SE; se usa la grid de menor LogLoss como fallback.")
    one_se_candidates <- cv_candidates %>%
      dplyr::slice(best_idx)
  }
  
  best_cv <- one_se_candidates %>%
    dplyr::arrange(
      p,
      sd_LogLoss,
      dplyr::desc(mean_BalAcc),
      dplyr::desc(mean_AUC),
      mean_LogLoss,
      GridID
    ) %>%
    dplyr::slice(1)
  
  candidate_grid_ids <- cv_candidates$GridID
  one_se_grid_ids <- one_se_candidates$GridID
  
  cv_summary <- cv_summary %>%
    dplyr::mutate(
      used_for_selection = GridID %in% candidate_grid_ids,
      best_logloss_cv = best_logloss,
      best_se_cv = best_se,
      logloss_1se_cutoff = logloss_cutoff,
      within_1se = GridID %in% one_se_grid_ids,
      selected_by_rule = GridID %in% best_cv$GridID,
      selection_rule = dplyr::case_when(
        selected_by_rule ~ "hierarchical_1se_min_p_tie_sd_balacc",
        within_1se ~ "within_1se_candidate",
        used_for_selection ~ "valid_nonperfect_candidate",
        TRUE ~ ""
      )
    ) %>%
    dplyr::arrange(
      dplyr::desc(selected_by_rule),
      dplyr::desc(within_1se),
      p,
      sd_LogLoss,
      dplyr::desc(mean_BalAcc),
      mean_LogLoss
    )
  list(
    folds      = folds,
    cv_results = cv_results,
    cv_summary = cv_summary,
    best_cv    = best_cv
  )
}

## ============================================================
## Auditoría en TEST para TODAS las grids BRMS_MCMC
## Entrena cada grid con todo train y evalúa en test
## No usar como selección principal; usar como sensibilidad
## ============================================================

.prepare_full_test_data <- function(return.list, grid_id, ref = "NP") {
  
  X_train <- return.list$train[[grid_id]]
  X_test  <- return.list$test[[grid_id]]
  
  y_train <- stats::relevel(refactor_levels(return.list$y_train), ref = ref)
  y_test  <- stats::relevel(refactor_levels(return.list$y_test),  ref = ref)
  
  if (is.null(rownames(X_train))) rownames(X_train) <- names(y_train)
  if (is.null(rownames(X_test)))  rownames(X_test)  <- names(y_test)
  
  df_train <- data.frame(y = y_train, X_train, check.names = FALSE)
  df_test  <- data.frame(y = y_test,  X_test,  check.names = FALSE)
  
  df_train$y <- droplevels(df_train$y)
  df_test$y  <- factor(df_test$y, levels = levels(df_train$y))
  
  common_cols <- setdiff(intersect(colnames(df_train), colnames(df_test)), "y")
  new_names <- make.unique(sanitize_names(common_cols), sep = "_")
  
  colnames(df_train)[match(common_cols, colnames(df_train))] <- new_names
  colnames(df_test)[match(common_cols,  colnames(df_test))]  <- new_names
  
  df_train <- df_train[, c("y", new_names), drop = FALSE]
  df_test  <- df_test[,  c("y", new_names), drop = FALSE]
  
  list(
    df_train = df_train,
    df_test  = df_test
  )
}


.fit_one_brms_test_grid <- function(task_id,
                                    tasks,
                                    return.list,
                                    thr_grid,
                                    ref = "NP",
                                    brms_iter = 4000,
                                    brms_chains = 2,
                                    seed = 123) {
  
  grid_id <- tasks$GridID[task_id]
  
  dat <- .prepare_full_test_data(
    return.list = return.list,
    grid_id     = grid_id,
    ref         = ref
  )
  
  df_train <- dat$df_train
  df_test  <- dat$df_test
  
  p <- ncol(df_train) - 1L
  
  backend <- .detect_backend_safe()
  ctrl <- .make_ctrl_safe(backend)
  
  priors <- make_multinom_priors(
    y_levels = levels(df_train$y),
    y = df_train$y,
    feature_names = setdiff(colnames(df_train), "y"),
    ref = ref,
    view_map = .get_global_or_null("view_map"),
    r2_score = .get_global_or_null("r2_score")
  )
  
  fit <- tryCatch(
    brms::brm(
      y ~ .,
      data      = df_train,
      family    = brms::categorical(),
      prior     = priors,
      chains    = brms_chains,
      iter      = brms_iter,
      warmup    = floor(brms_iter / 2),
      refresh   = 0,
      backend   = backend,
      control   = ctrl,
      cores     = 1,
      seed      = seed + 1000L * grid_id,
      save_pars = brms::save_pars(all = TRUE)
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(list(
      summary = tibble::tibble(
        GridID = grid_id,
        Model = "BRMS_MCMC",
        p = p,
        Accuracy = NA_real_,
        BalAcc = NA_real_,
        Kappa = NA_real_,
        AUC = NA_real_,
        LogLoss = NA_real_,
        error = conditionMessage(fit)
      ),
      predictions = NULL,
      confusion = NULL
    ))
  }
  
  P <- tryCatch(
    brms_probs(fit, df_test)[, levels(df_test$y), drop = FALSE],
    error = function(e) e
  )
  
  if (inherits(P, "error")) {
    return(list(
      summary = tibble::tibble(
        GridID = grid_id,
        Model = "BRMS_MCMC",
        p = p,
        Accuracy = NA_real_,
        BalAcc = NA_real_,
        Kappa = NA_real_,
        AUC = NA_real_,
        LogLoss = NA_real_,
        error = conditionMessage(P)
      ),
      predictions = NULL,
      confusion = NULL
    ))
  }
  
  met <- metrics_from_probs(P, df_test$y)
  
  pred_lab <- factor(
    colnames(P)[max.col(P, ties.method = "first")],
    levels = levels(df_test$y)
  )
  
  cm <- caret::confusionMatrix(pred_lab, df_test$y)
  
  prob_df <- as.data.frame(P, check.names = FALSE)
  colnames(prob_df) <- paste0("prob_", colnames(prob_df))
  
  pred_df <- dplyr::bind_cols(
    tibble::tibble(
      GridID   = grid_id,
      p        = p,
      id       = rownames(df_test),
      real     = df_test$y,
      predicho = pred_lab,
      correcto = df_test$y == pred_lab,
      prob_max = apply(P, 1, max)
    ),
    prob_df
  )
  
  conf_df <- as.data.frame(cm$table)
  colnames(conf_df) <- c("Prediction", "Reference", "n")
  conf_df$GridID <- grid_id
  conf_df$p <- p
  
  summary_df <- dplyr::bind_cols(
    tibble::tibble(
      GridID = grid_id,
      Model = "BRMS_MCMC",
      p = p
    ),
    met,
    tibble::tibble(error = NA_character_)
  )
  
  list(
    summary = summary_df,
    predictions = pred_df,
    confusion = conf_df
  )
}

## ============================================================
## Estabilidad CV repetida sobre grids candidatas
## No usa test para seleccionar
## ============================================================

evaluate_cv_stability_candidates <- function(CV_BRM,
                                             return.list,
                                             thr_grid,
                                             ref = "NP",
                                             K = 3,
                                             p_min = 35,
                                             p_max = 50,
                                             n_candidates = 12,
                                             cv_seeds = c(123, 321, 456, 789, 1001),
                                             brms_iter = 2000,
                                             brms_chains = 2,
                                             cores = 2) {
  
  ## 1) Selecciona candidatas desde la primera CV usando la regla 1-SE.
  ## La CV repetida NO crea un score nuevo; solo se usa como desempate.
  
  candidate_pool <- CV_BRM$cv_summary %>%
    dplyr::filter(
      used_for_selection %in% TRUE,
      within_1se %in% TRUE
    ) %>%
    dplyr::arrange(
      p,
      sd_LogLoss,
      dplyr::desc(mean_BalAcc),
      dplyr::desc(mean_AUC),
      mean_LogLoss,
      GridID
    )
  
  if (nrow(candidate_pool) == 0) {
    
    warning(
      "No hay candidatas within_1se en CV_BRM$cv_summary. ",
      "Se usará la grid marcada como selected_by_rule."
    )
    
    candidate_pool <- CV_BRM$cv_summary %>%
      dplyr::filter(selected_by_rule %in% TRUE)
  }
  
  if (nrow(candidate_pool) == 0) {
    stop("No hay grids candidatas para estabilidad.")
  }
  
  candidate_ids <- candidate_pool %>%
    dplyr::slice_head(n = n_candidates) %>%
    dplyr::pull(GridID) %>%
    unique()
  
  cat("\n--- GridIDs candidatas para estabilidad ---\n")
  print(candidate_ids)
  
  ## 2) Subset de thr_grid y return.list
  thr_grid_stab <- thr_grid[candidate_ids, , drop = FALSE] %>%
    dplyr::mutate(OrigGridID = candidate_ids)
  
  return_stab <- list(
    train   = return.list$train[candidate_ids],
    test    = return.list$test[candidate_ids],
    y_train = return.list$y_train,
    y_test  = return.list$y_test
  )
  
  ## 3) Repite CV con distintas semillas
  cv_rep_list <- lapply(seq_along(cv_seeds), function(ii) {
    
    s <- cv_seeds[ii]
    
    cat("\n============================================================\n")
    cat("CV repetida", ii, "de", length(cv_seeds), "| seed =", s, "\n")
    cat("============================================================\n")
    
    cv_i <- evaluate_brms_cv_grid(
      return.list  = return_stab,
      thr_grid     = thr_grid_stab,
      ref          = ref,
      K            = K,
      p_min        = p_min,
      p_max        = p_max,
      brms_iter    = brms_iter,
      brms_chains  = brms_chains,
      cores        = cores,
      seed         = s
    )
    
    best_local_id <- as.integer(cv_i$best_cv$GridID[1])
    best_orig_id  <- thr_grid_stab$OrigGridID[best_local_id]
    
    cv_i$cv_summary %>%
      dplyr::mutate(
        OrigGridID = thr_grid_stab$OrigGridID[GridID],
        repeat_id = ii,
        cv_seed = s,
        selected_in_repeat = OrigGridID == best_orig_id
      )
  })
  
  cv_repeated_summary <- dplyr::bind_rows(cv_rep_list)
  
  ## 4) Resumen de estabilidad por grid original
  stability_summary <- cv_repeated_summary %>%
    dplyr::group_by(OrigGridID) %>%
    dplyr::summarise(
      n_repeats = dplyr::n(),
      n_selected = sum(selected_in_repeat, na.rm = TRUE),
      selection_freq = n_selected / n_repeats,
      
      p = median(p, na.rm = TRUE),
      
      mean_cv_LogLoss = mean(mean_LogLoss, na.rm = TRUE),
      sd_cv_LogLoss   = sd(mean_LogLoss, na.rm = TRUE),
      min_cv_LogLoss  = min(mean_LogLoss, na.rm = TRUE),
      max_cv_LogLoss  = max(mean_LogLoss, na.rm = TRUE),
      
      mean_cv_BalAcc = mean(mean_BalAcc, na.rm = TRUE),
      sd_cv_BalAcc   = sd(mean_BalAcc, na.rm = TRUE),
      
      mean_cv_AUC = mean(mean_AUC, na.rm = TRUE),
      sd_cv_AUC   = sd(mean_AUC, na.rm = TRUE),
      
      prop_perfect_cv = mean(perfect_cv, na.rm = TRUE),
      
      q_tx = dplyr::first(q_tx),
      q_pr = dplyr::first(q_pr),
      q_me = dplyr::first(q_me),
      q_cl = dplyr::first(q_cl),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      candidate_pool %>%
        dplyr::select(
          OrigGridID = GridID,
          initial_mean_LogLoss = mean_LogLoss,
          initial_se_LogLoss = se_LogLoss,
          initial_sd_LogLoss = sd_LogLoss,
          initial_mean_BalAcc = mean_BalAcc,
          initial_mean_AUC = mean_AUC,
          initial_p = p,
          logloss_1se_cutoff,
          within_1se,
          selected_by_rule
        ),
      by = "OrigGridID"
    ) %>%
    dplyr::mutate(
      final_selection_rule = "hierarchical_1se_min_p_tie_repeated_sd_balacc"
    ) %>%
    dplyr::arrange(
      p,
      sd_cv_LogLoss,
      dplyr::desc(selection_freq),
      dplyr::desc(mean_cv_BalAcc),
      mean_cv_LogLoss,
      dplyr::desc(mean_cv_AUC),
      OrigGridID
    )
  
  ## 5) Selección final jerárquica.
  ## No se excluyen grids mediante un umbral adicional de prop_perfect_cv.
  ## prop_perfect_cv queda solo como auditoría.
  
  best_stable <- stability_summary %>%
    dplyr::arrange(
      p,
      sd_cv_LogLoss,
      dplyr::desc(selection_freq),
      dplyr::desc(mean_cv_BalAcc),
      mean_cv_LogLoss,
      dplyr::desc(mean_cv_AUC),
      OrigGridID
    ) %>%
    dplyr::slice(1)
  
  list(
    candidate_ids = candidate_ids,
    cv_repeated_summary = cv_repeated_summary,
    stability_summary = stability_summary,
    best_stable = best_stable
  )
}
evaluate_test_all_grids_brms <- function(return.list,
                                         thr_grid,
                                         ref = "NP",
                                         brms_iter = 4000,
                                         brms_chains = 2,
                                         cores = 2,
                                         seed = 123) {
  
  tasks <- data.frame(GridID = seq_len(nrow(thr_grid)))
  
  if (cores > 1) {
    
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterSetRNGStream(cl, seed)
    
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(dplyr)
        library(tibble)
        library(brms)
        library(caret)
        library(pROC)
      })
      NULL
    })
    
    parallel::clusterExport(
      cl,
      varlist = c(
        ## Objetos locales necesarios en los workers
        "tasks",
        "return.list",
        "thr_grid",
        "ref",
        "brms_iter",
        "brms_chains",
        "seed",
        
        ## Funciones auxiliares
        "refactor_levels",
        "sanitize_names",
        "metrics_from_probs",
        "brms_probs",
        
        ## Priors
        "PRIOR_PREVALENCE_ALPHA",
        "PRIOR_INTERCEPT_SD",
        "PRIOR_COEF_SD",
        
        ".get_global_or_null",
        ".infer_feature_view_tbl",
        ".make_class_prevalence_tbl",
        "make_multinom_prior_audit_tbl",
        "make_multinom_priors",
        
        ## BRMS/test
        ".detect_backend_safe",
        ".make_ctrl_safe",
        ".prepare_full_test_data",
        ".fit_one_brms_test_grid"
      ),
      envir = environment()
    )
    
    out_list <- parallel::parLapply(
      cl,
      seq_len(nrow(tasks)),
      function(ii) {
        .fit_one_brms_test_grid(
          task_id     = ii,
          tasks       = tasks,
          return.list = return.list,
          thr_grid    = thr_grid,
          ref         = ref,
          brms_iter   = brms_iter,
          brms_chains = brms_chains,
          seed        = seed
        )
      }
    )
    
  } else {
    
    out_list <- lapply(seq_len(nrow(tasks)), function(ii) {
      .fit_one_brms_test_grid(
        task_id     = ii,
        tasks       = tasks,
        return.list = return.list,
        thr_grid    = thr_grid,
        ref         = ref,
        brms_iter   = brms_iter,
        brms_chains = brms_chains,
        seed        = seed
      )
    })
  }
  
  test_summary <- dplyr::bind_rows(lapply(out_list, `[[`, "summary")) %>%
    dplyr::left_join(
      thr_grid %>% dplyr::mutate(GridID = dplyr::row_number()),
      by = "GridID"
    ) %>%
    dplyr::arrange(LogLoss, dplyr::desc(BalAcc), dplyr::desc(AUC), p)
  
  test_predictions <- dplyr::bind_rows(lapply(out_list, `[[`, "predictions")) %>%
    dplyr::left_join(
      thr_grid %>% dplyr::mutate(GridID = dplyr::row_number()),
      by = "GridID"
    )
  
  test_confusion <- dplyr::bind_rows(lapply(out_list, `[[`, "confusion")) %>%
    dplyr::left_join(
      thr_grid %>% dplyr::mutate(GridID = dplyr::row_number()),
      by = "GridID"
    )
  
  list(
    test_summary = test_summary,
    test_predictions = test_predictions,
    test_confusion = test_confusion
  )
}
##===ejecucion=====

##=== ejecucion sin datos_para_modelar.rds =====

if (!file.exists(RFM_FILE)) {
  stop(
    "No existe RFM_FILE: ", RFM_FILE,
    "\nPrimero corre el preprocesamiento para ANALYSIS_NAME = ", ANALYSIS_NAME
  )
}

if (!file.exists(MOFA_MODEL_FILE)) {
  stop(
    "No existe MOFA_MODEL_FILE: ", MOFA_MODEL_FILE,
    "\nPrimero corre la integración MOFA para ANALYSIS_NAME = ", ANALYSIS_NAME
  )
}

rfm <- readRDS(RFM_FILE)
tst_Data <- rfm

mod <- readRDS(MOFA_MODEL_FILE)

preprocess_config <- NULL

if (file.exists(PREPROCESS_CONFIG_FILE)) {
  preprocess_config <- readRDS(PREPROCESS_CONFIG_FILE)
  
  cat("\n--- Configuración de preprocesamiento detectada ---\n")
  print(preprocess_config)
} else {
  warning("No se encontró run_config_preprocess.rds en: ", PREPROCESS_CONFIG_FILE)
}

## -----------------------------
## 1) Auditoría básica del split
## -----------------------------

view_map <- c(
  transcriptomica = "tx",
  proteomica      = "pr",
  metabolomica    = "me",
  clinical        = "cl"
)

ids_train <- rfm$grupo_train$id
ids_test  <- rfm$grupo_test$id

audit_ids <- lapply(names(view_map), function(v) {
  nm <- view_map[[v]]
  data.frame(
    vista = v,
    train_ok = identical(colnames(rfm[[nm]]$train), ids_train),
    test_ok  = identical(colnames(rfm[[nm]]$test),  ids_test),
    n_train  = ncol(rfm[[nm]]$train),
    n_test   = ncol(rfm[[nm]]$test)
  )
}) |> dplyr::bind_rows()

print(audit_ids)

if (!all(audit_ids$train_ok)) {
  stop("El orden de muestras train no coincide entre rfm$grupo_train y alguna vista.")
}

if (!all(audit_ids$test_ok)) {
  stop("El orden de muestras test no coincide entre rfm$grupo_test y alguna vista.")
}

# ## ---------------------------------------------
# ## 2) Extraer factores reales desde el modelo MOFA
# ## ---------------------------------------------
# 
# extract_mofa_factors <- function(mod, ids_train, factors = 1:2) {
#   
#   F_raw <- MOFA2::get_factors(
#     mod,
#     factors = factors,
#     as.data.frame = FALSE
#   )
#   
#   if (is.list(F_raw)) {
#     F <- F_raw[[1]]
#   } else {
#     F <- F_raw
#   }
#   
#   F <- as.matrix(F)
#   
#   ## Asegurar orientación: muestras x factores
#   if (!all(ids_train %in% rownames(F)) && all(ids_train %in% colnames(F))) {
#     F <- t(F)
#   }
#   
#   missing_ids <- setdiff(ids_train, rownames(F))
#   if (length(missing_ids) > 0) {
#     stop(
#       "Estos ids_train no están en los factores MOFA: ",
#       paste(missing_ids, collapse = ", ")
#     )
#   }
#   
#   F <- F[ids_train, , drop = FALSE]
#   F <- F[, seq_along(factors), drop = FALSE]
#   colnames(F) <- paste0("Factor", factors)
#   
#   F
# }
# 
# Z_seed <- extract_mofa_factors(
#   mod = mod,
#   ids_train = ids_train,
#   factors = 1:2
# )

## ---------------------------------------------
## 3) Extraer pesos/loadings desde el modelo MOFA
## ---------------------------------------------

extract_mofa_weights <- function(mod, rfm, view_map, factors = 1:2) {
  
  W_raw <- MOFA2::get_weights(
    mod,
    views = "all",
    factors = factors,
    as.data.frame = FALSE
  )
  
  ## Algunos objetos MOFA devuelven una lista anidada por grupo
  if (
    is.list(W_raw) &&
    length(W_raw) == 1 &&
    is.list(W_raw[[1]]) &&
    !is.matrix(W_raw[[1]])
  ) {
    W_raw <- W_raw[[1]]
  }
  
  if (is.null(names(W_raw))) {
    names(W_raw) <- names(view_map)[seq_along(W_raw)]
  }
  
  out <- list()
  
  for (v in names(view_map)) {
    
    nm <- view_map[[v]]
    
    if (!v %in% names(W_raw)) {
      stop("No encuentro la vista ", v, " en los pesos del modelo MOFA.")
    }
    
    W <- as.matrix(W_raw[[v]])
    feats_rfm <- rownames(rfm[[nm]]$train)
    
    ## Asegurar orientación: features x factores
    common_rows <- intersect(rownames(W), feats_rfm)
    common_cols <- intersect(colnames(W), feats_rfm)
    
    if (length(common_rows) == 0 && length(common_cols) > 0) {
      W <- t(W)
    }
    
    common_feats <- intersect(rownames(W), feats_rfm)
    
    if (length(common_feats) == 0) {
      stop("No hay features comunes entre pesos MOFA y rfm para la vista: ", v)
    }
    
    W <- W[common_feats, seq_along(factors), drop = FALSE]
    colnames(W) <- paste0("Factor", factors)
    
    out[[v]] <- W
  }
  
  out
}

pesos_list <- extract_mofa_weights(
  mod = mod,
  rfm = rfm,
  view_map = view_map,
  factors = 1:2
)

## ---------------------------------------------
## 4) Respuesta y orden correcto de entrenamiento
## ---------------------------------------------

# y_seed <- rfm$grupo_train$grupo[match(rownames(Z_seed), rfm$grupo_train$id)]
# y_seed <- factor(y_seed)
# names(y_seed) <- rownames(Z_seed)

## ---------------------------------------------
## 5) Auditoría de dimensiones
## ---------------------------------------------

# cat("\n--- Auditoría Z_seed / y_seed ---\n")
# print(dim(Z_seed))
# print(table(y_seed))
# print(all(names(y_seed) == rownames(Z_seed)))

cat("\n--- Auditoría pesos_list ---\n")
print(lapply(pesos_list, dim))

cat("\n--- Features comunes pesos vs rfm ---\n")
audit_weights <- lapply(names(view_map), function(v) {
  nm <- view_map[[v]]
  data.frame(
    vista = v,
    n_pesos = nrow(pesos_list[[v]]),
    n_rfm_train = nrow(rfm[[nm]]$train),
    n_common_train = length(intersect(rownames(pesos_list[[v]]), rownames(rfm[[nm]]$train))),
    n_common_test  = length(intersect(rownames(pesos_list[[v]]), rownames(rfm[[nm]]$test)))
  )
}) |> dplyr::bind_rows()

print(audit_weights)

##===grilla de thresholds=====

.get_mofa_view_weights <- function(mod,
                                   view_map,
                                   factors = 1:2,
                                   weight_source = c("r2_total", "r2_per_factor_sum", "equal"),
                                   group_mofa = NULL,
                                   weight_power = 1) {
  
  weight_source <- match.arg(weight_source)
  view_names <- names(view_map)
  
  if (identical(weight_source, "equal")) {
    w <- rep(1, length(view_names))
    names(w) <- view_names
    return(w / sum(w))
  }
  
  ve <- MOFA2::get_variance_explained(mod)
  
  if (is.null(group_mofa)) {
    if (is.list(ve$r2_total)) {
      group_mofa <- names(ve$r2_total)[1]
    } else {
      group_mofa <- "group1"
    }
  }
  
  if (identical(weight_source, "r2_total")) {
    
    r2 <- if (is.list(ve$r2_total)) {
      ve$r2_total[[group_mofa]]
    } else {
      ve$r2_total
    }
    
    r2 <- r2[view_names]
    w <- as.numeric(r2)
    names(w) <- view_names
  }
  
  if (identical(weight_source, "r2_per_factor_sum")) {
    
    r2pf <- if (is.list(ve$r2_per_factor)) {
      ve$r2_per_factor[[group_mofa]]
    } else {
      ve$r2_per_factor
    }
    
    factor_names <- paste0("Factor", factors)
    factor_names <- intersect(factor_names, rownames(r2pf))
    
    if (length(factor_names) == 0) {
      stop("No encuentro los factores solicitados en r2_per_factor.")
    }
    
    r2pf <- r2pf[factor_names, view_names, drop = FALSE]
    w <- colSums(r2pf, na.rm = TRUE)
  }
  
  w[is.na(w)] <- 0
  w[w < 0] <- 0
  
  if (sum(w) == 0) {
    w <- rep(1, length(view_names))
    names(w) <- view_names
  }
  
  w <- w ^ weight_power
  w / sum(w)
}


.allocate_budget_by_weight <- function(P,
                                       weights,
                                       n_available,
                                       min_keep = NULL,
                                       max_keep = NULL) {
  
  view_names <- names(weights)
  
  if (is.null(min_keep)) {
    min_keep <- rep(1L, length(view_names))
    names(min_keep) <- view_names
  }
  
  if (is.null(max_keep)) {
    max_keep <- n_available
  }
  
  min_keep <- min_keep[view_names]
  max_keep <- max_keep[view_names]
  n_available <- n_available[view_names]
  
  min_keep[is.na(min_keep)] <- 1L
  max_keep[is.na(max_keep)] <- n_available[is.na(max_keep)]
  
  min_keep <- pmax(0L, as.integer(min_keep))
  max_keep <- pmin(as.integer(max_keep), as.integer(n_available))
  min_keep <- pmin(min_keep, max_keep)
  
  if (P < sum(min_keep)) {
    stop(
      "p_total = ", P,
      " es menor que la suma de min_keep = ", sum(min_keep)
    )
  }
  
  if (P > sum(max_keep)) {
    stop(
      "p_total = ", P,
      " es mayor que la suma de max_keep = ", sum(max_keep)
    )
  }
  
  raw <- P * weights
  n_keep <- floor(raw)
  names(n_keep) <- view_names
  
  n_keep <- pmax(n_keep, min_keep)
  n_keep <- pmin(n_keep, max_keep)
  
  remainder <- raw - floor(raw)
  names(remainder) <- view_names
  
  while (sum(n_keep) < P) {
    candidates <- view_names[n_keep < max_keep]
    if (length(candidates) == 0) break
    
    add_to <- candidates[which.max(remainder[candidates])]
    n_keep[add_to] <- n_keep[add_to] + 1L
    remainder[add_to] <- -Inf
  }
  
  while (sum(n_keep) > P) {
    candidates <- view_names[n_keep > min_keep]
    if (length(candidates) == 0) break
    
    remove_from <- candidates[which.min(remainder[candidates])]
    n_keep[remove_from] <- n_keep[remove_from] - 1L
    remainder[remove_from] <- Inf
  }
  out <- as.integer(n_keep)
  names(out) <- view_names
  out
}


make_mofa_weighted_grid <- function(mod,
                                    rfm,
                                    view_map,
                                    p_total_grid =  c(25, 35, 45, 55, 65, 75, 85),
                                    factors = 1:2,
                                    weight_source = c("r2_total", "r2_per_factor_sum", "equal"),
                                    group_mofa = NULL,
                                    weight_power = 1,
                                    min_keep = NULL,
                                    max_keep = NULL) {
  
  weight_source <- match.arg(weight_source)
  
  view_names <- names(view_map)
  view_codes <- unname(view_map)
  
  weights <- .get_mofa_view_weights(
    mod           = mod,
    view_map      = view_map,
    factors       = factors,
    weight_source = weight_source,
    group_mofa    = group_mofa,
    weight_power  = weight_power
  )
  
  n_available <- vapply(view_names, function(v) {
    code <- view_map[[v]]
    nrow(rfm[[code]]$train)
  }, integer(1))
  
  if (is.null(min_keep)) {
    min_keep <- rep(1L, length(view_names))
    names(min_keep) <- view_names
  }
  
  if (is.null(max_keep)) {
    max_keep <- n_available
  }
  
  out <- lapply(p_total_grid, function(P) {
    
    n_keep <- .allocate_budget_by_weight(
      P           = P,
      weights     = weights,
      n_available = n_available,
      min_keep    = min_keep,
      max_keep    = max_keep
    )
    
    n_keep <- n_keep[view_names]
    
    if (any(is.na(n_keep))) {
      stop(
        "n_keep tiene NA. Revisa names(n_keep). names actuales: ",
        paste(names(n_keep), collapse = ", ")
      )
    }
    
    q_keep <- 1 - (n_keep / n_available[view_names])
    names(q_keep) <- view_names
    q_keep <- pmax(pmin(q_keep, 0.999), 0)
    
    df <- data.frame(
      p_target = P,
      p_actual = sum(n_keep),
      weight_source = weight_source,
      weight_power = weight_power,
      stringsAsFactors = FALSE
    )
    
    for (v in view_names) {
      code <- view_map[[v]]
      df[[paste0("weight_", code)]] <- weights[v]
      df[[paste0("n_", code)]] <- as.integer(n_keep[v])
      df[[paste0("q_", code)]] <- as.numeric(q_keep[v])
    }
    
    df
  })
  
  dplyr::bind_rows(out)
}


thr_grid <- make_mofa_weighted_grid(
  mod           = mod,
  rfm           = rfm,
  view_map      = view_map,
  p_total_grid  =c(25, 35, 45, 55, 65, 75, 85),
  factors       = 1:2,
  weight_source = "r2_total",
  weight_power  = 1,
  min_keep      = NULL,
  max_keep      = NULL
)
cat("\n--- Grid ponderada por varianza explicada de MOFA ---\n")
print(thr_grid)

cat("\n--- Auditoría de thr_grid ---\n")
print(str(thr_grid))

n_cols_grid <- paste0("n_", unname(view_map))
q_cols_grid <- paste0("q_", unname(view_map))

cat("\nColumnas n esperadas:\n")
print(n_cols_grid)

cat("\nColumnas q esperadas:\n")
print(q_cols_grid)

cat("\nValores n por grid:\n")
print(thr_grid[, n_cols_grid, drop = FALSE])

cat("\nValores q por grid:\n")
print(thr_grid[, q_cols_grid, drop = FALSE])

if (any(is.na(thr_grid[, n_cols_grid, drop = FALSE]))) {
  stop("ERROR: hay NA en columnas n_* de thr_grid.")
}

if (any(is.na(thr_grid[, q_cols_grid, drop = FALSE]))) {
  stop("ERROR: hay NA en columnas q_* de thr_grid.")
}

cat("\nSuma de n_* por grid:\n")
print(rowSums(thr_grid[, n_cols_grid, drop = FALSE]))

write.csv(
  thr_grid,
  out_table("thr_grid_mofa_weighted.csv"),
  row.names = FALSE
)

saveRDS(
  thr_grid,
  out_rds("thr_grid_mofa_weighted.rds")
)


metadata_by_code <- list()

if (!is.null(rfm$features_metadata)) {
  metadata_by_code[[view_map[["transcriptomica"]]]] <- rfm$features_metadata
}

## ============================================================
## SCORE MOFA PARA TODAS LAS FEATURES CANDIDATAS DE TODAS LAS GRIDS
## No modifica modelos. No guarda archivos. Solo imprime y deja objetos en memoria.
## Requiere:
##   mod, rfm, view_map, thr_grid
## Opcional:
##   grid_id, df_train
## ============================================================

cat("\n============================================================\n")
cat("SCORE MOFA DE BIOMARCADORES CANDIDATOS\n")
cat("============================================================\n")

## -----------------------------
## 0) Fallback seguro para sanitize_names
## -----------------------------
if (!exists("sanitize_names")) {
  sanitize_names <- function(v){
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
}

.print_tbl <- function(x, n = Inf) {
  print(tibble::as_tibble(x), n = n)
}

## -----------------------------
## 1) Chequeos mínimos
## -----------------------------
stopifnot(exists("mod"))
stopifnot(exists("rfm"))
stopifnot(exists("view_map"))
stopifnot(exists("thr_grid"))

view_names <- names(view_map)
view_codes <- unname(view_map)

if (is.null(view_names) || any(view_names == "")) {
  stop("view_map debe tener nombres de vistas MOFA.")
}

if (any(is.na(view_map)) || any(view_map == "")) {
  stop("view_map debe mapear cada vista MOFA a un código de rfm.")
}

for (v in view_names) {
  code <- view_map[[v]]
  if (!code %in% names(rfm)) {
    stop("No encuentro rfm[['", code, "']] para la vista ", v)
  }
  if (is.null(rfm[[code]]$train)) {
    stop("rfm[['", code, "']] no tiene matriz train.")
  }
}

## Usar todos los factores disponibles en el modelo.
## Si quieres restringirlo, cambia esta línea a: factors_score <- 1:2
factors_score <- 1:2
factor_names_score <- paste0("Factor", factors_score)

cat("\n--- Factores usados para score ---\n")
print(factor_names_score)

## ============================================================
## 2) Extraer pesos MOFA de forma robusta
## ============================================================

.extract_weights_for_score <- function(mod, rfm, view_map, factors_score) {
  
  W_raw <- MOFA2::get_weights(
    mod,
    views = "all",
    factors = factors_score,
    as.data.frame = FALSE
  )
  
  ## Algunos MOFA devuelven lista anidada por grupo
  if (
    is.list(W_raw) &&
    length(W_raw) == 1 &&
    is.list(W_raw[[1]]) &&
    !is.matrix(W_raw[[1]])
  ) {
    W_raw <- W_raw[[1]]
  }
  
  if (is.null(names(W_raw))) {
    names(W_raw) <- names(view_map)[seq_along(W_raw)]
  }
  
  out <- list()
  
  for (v in names(view_map)) {
    
    code <- view_map[[v]]
    
    if (!v %in% names(W_raw)) {
      stop("No encuentro la vista '", v, "' en los pesos MOFA.")
    }
    
    W <- as.matrix(W_raw[[v]])
    feats_rfm <- rownames(rfm[[code]]$train)
    
    if (is.null(rownames(W)) || is.null(colnames(W))) {
      stop("La matriz de pesos MOFA para ", v, " no tiene dimnames.")
    }
    
    ## Asegurar orientación: features x factores
    common_rows <- intersect(rownames(W), feats_rfm)
    common_cols <- intersect(colnames(W), feats_rfm)
    
    if (length(common_rows) == 0 && length(common_cols) > 0) {
      W <- t(W)
    }
    
    common_feats <- intersect(rownames(W), feats_rfm)
    
    if (length(common_feats) == 0) {
      stop("No hay features comunes entre pesos MOFA y rfm para vista: ", v)
    }
    
    factor_cols <- paste0("Factor", factors_score)
    factor_cols <- intersect(factor_cols, colnames(W))
    
    if (length(factor_cols) == 0) {
      stop("No encuentro factores solicitados en pesos MOFA para vista: ", v)
    }
    
    W <- W[common_feats, factor_cols, drop = FALSE]
    
    out[[v]] <- W
  }
  
  out
}

weights_by_view <- .extract_weights_for_score(
  mod = mod,
  rfm = rfm,
  view_map = view_map,
  factors_score = factors_score
)

cat("\n--- Dimensiones pesos MOFA usados para score ---\n")
print(lapply(weights_by_view, dim))

## ============================================================
## 3) Extraer varianza explicada MOFA por vista y factor
## ============================================================

.extract_r2_for_score <- function(mod, view_map, factors_score, group_mofa = NULL) {
  
  ve <- MOFA2::get_variance_explained(mod)
  
  if (is.null(group_mofa)) {
    if (is.list(ve$r2_total)) {
      group_mofa <- names(ve$r2_total)[1]
    } else {
      group_mofa <- "group1"
    }
  }
  
  view_names <- names(view_map)
  factor_names <- paste0("Factor", factors_score)
  
  r2_total <- if (is.list(ve$r2_total)) {
    ve$r2_total[[group_mofa]]
  } else {
    ve$r2_total
  }
  
  r2_total <- as.numeric(r2_total[view_names])
  names(r2_total) <- view_names
  
  r2_per_factor <- if (is.list(ve$r2_per_factor)) {
    ve$r2_per_factor[[group_mofa]]
  } else {
    ve$r2_per_factor
  }
  
  factor_names_ok <- intersect(factor_names, rownames(r2_per_factor))
  view_names_ok <- intersect(view_names, colnames(r2_per_factor))
  
  r2_per_factor <- r2_per_factor[factor_names_ok, view_names_ok, drop = FALSE]
  
  list(
    group_mofa = group_mofa,
    r2_total = r2_total,
    r2_per_factor = r2_per_factor
  )
}

r2_score <- .extract_r2_for_score(
  mod = mod,
  view_map = view_map,
  factors_score = factors_score
)

cat("\n--- Grupo MOFA usado para varianza explicada ---\n")
print(r2_score$group_mofa)

cat("\n--- R2 total por vista ---\n")
print(r2_score$r2_total)

cat("\n--- R2 por factor y vista ---\n")
print(r2_score$r2_per_factor)

## ============================================================
## 4) Tabla base de importancia MOFA por feature
## ============================================================
.make_feature_score_base <- function(weights_by_view,
                                     r2_score,
                                     view_map,
                                     metadata_by_code = NULL) {
  
  view_names <- names(view_map)
  
  out <- lapply(view_names, function(v) {
    
    code <- view_map[[v]]
    W <- weights_by_view[[v]]
    
    absW <- abs(W)
    
    ## R2 de cada factor en esta vista
    factor_r2_vec <- rep(0, ncol(absW))
    names(factor_r2_vec) <- colnames(absW)
    
    common_factors <- intersect(colnames(absW), rownames(r2_score$r2_per_factor))
    
    if (length(common_factors) > 0 && v %in% colnames(r2_score$r2_per_factor)) {
      factor_r2_vec[common_factors] <- as.numeric(
        r2_score$r2_per_factor[common_factors, v]
      )
    }
    
    factor_r2_vec[is.na(factor_r2_vec)] <- 0
    
    ## Si por alguna razón todos los R2 son cero, usar pesos iguales
    if (sum(factor_r2_vec, na.rm = TRUE) == 0) {
      factor_r2_vec[] <- 1
    }
    
    ## Score directo MOFA:
    ## suma_k |loading_feature,k| * R2_factor,k_en_vista
    weighted_absW <- sweep(absW, 2, factor_r2_vec, "*")
    
    score_mofa_for_selection <- rowSums(weighted_absW, na.rm = TRUE)
    
    ## Factor dominante según contribución ponderada, no solo loading absoluto
    dominant_factor <- colnames(weighted_absW)[
      max.col(weighted_absW, ties.method = "first")
    ]
    
    loading_abs_max <- apply(absW, 1, max, na.rm = TRUE)
    
    dominant_factor_loading <- absW[
      cbind(seq_len(nrow(absW)), match(dominant_factor, colnames(absW)))
    ]
    
    dominant_factor_contribution <- weighted_absW[
      cbind(seq_len(nrow(weighted_absW)), match(dominant_factor, colnames(weighted_absW)))
    ]
    
    factor_r2 <- factor_r2_vec[dominant_factor]
    
    ## Percentil del loading dentro de vista: solo descriptivo
    loading_rank <- rank(loading_abs_max, ties.method = "average")
    loading_percentile_view <- loading_rank / length(loading_rank)
    
    ## Percentil del score MOFA ponderado dentro de vista
    score_rank <- rank(score_mofa_for_selection, ties.method = "average")
    score_mofa_percentile_view <- score_rank / length(score_rank)
    
    ## Escalado min-max del score dentro de vista
    rng_score <- range(score_mofa_for_selection, na.rm = TRUE)
    
    if (diff(rng_score) == 0) {
      score_mofa_scaled_view <- rep(1, length(score_mofa_for_selection))
    } else {
      score_mofa_scaled_view <- 
        (score_mofa_for_selection - rng_score[1]) / diff(rng_score)
    }
    
    ## Escalado min-max del loading, solo descriptivo
    rng_loading <- range(loading_abs_max, na.rm = TRUE)
    
    if (diff(rng_loading) == 0) {
      loading_scaled_view <- rep(1, length(loading_abs_max))
    } else {
      loading_scaled_view <- 
        (loading_abs_max - rng_loading[1]) / diff(rng_loading)
    }
    
    view_r2_total <- as.numeric(r2_score$r2_total[v])
    
    feature_original <- rownames(W)
    feature_label <- feature_original
    
    ## Mapeo opcional EntrezGeneID -> GeneSymbol
    if (!is.null(metadata_by_code) && code %in% names(metadata_by_code)) {
      
      metadata_i <- metadata_by_code[[code]]
      
      if (!is.null(metadata_i) && all(c("EntrezGeneID", "GeneSymbol") %in% colnames(metadata_i))) {
        
        map_dic <- metadata_i %>%
          dplyr::transmute(
            feature_original = as.character(EntrezGeneID),
            feature_label_meta = as.character(GeneSymbol)
          ) %>%
          dplyr::filter(!is.na(feature_original), feature_original != "") %>%
          dplyr::group_by(feature_original) %>%
          dplyr::summarise(
            feature_label_meta = dplyr::first(
              feature_label_meta[!is.na(feature_label_meta) & feature_label_meta != ""],
              default = NA_character_
            ),
            .groups = "drop"
          )
        
        label_tbl <- tibble::tibble(
          feature_original = feature_original
        ) %>%
          dplyr::left_join(map_dic, by = "feature_original")
        
        feature_label <- ifelse(
          !is.na(label_tbl$feature_label_meta) & label_tbl$feature_label_meta != "",
          label_tbl$feature_label_meta,
          label_tbl$feature_original
        )
      }
    }
    
    tibble::tibble(
      view = v,
      view_code = code,
      feature_original = feature_original,
      feature_label = feature_label,
      feature_prefixed_original = paste0(code, "__", feature_original),
      feature_prefixed_label = paste0(code, "__", feature_label),
      feature_model = sanitize_names(paste0(code, "__", feature_label)),
      
      dominant_factor = dominant_factor,
      loading_abs_max = as.numeric(loading_abs_max),
      dominant_factor_loading = as.numeric(dominant_factor_loading),
      factor_r2 = as.numeric(factor_r2),
      dominant_factor_contribution = as.numeric(dominant_factor_contribution),
      
      view_r2_total = view_r2_total,
      
      loading_percentile_view = as.numeric(loading_percentile_view),
      loading_scaled_view = as.numeric(loading_scaled_view),
      
      score_mofa_for_selection = as.numeric(score_mofa_for_selection),
      score_mofa_percentile_view = as.numeric(score_mofa_percentile_view),
      score_mofa_scaled_view = as.numeric(score_mofa_scaled_view)
    )
  })
  
  dplyr::bind_rows(out)
}

feature_score_base <- .make_feature_score_base(
  weights_by_view  = weights_by_view,
  r2_score         = r2_score,
  view_map         = view_map,
  metadata_by_code = metadata_by_code
)


cat("\n--- Base de score MOFA por feature, top 20 por loading ---\n")
.print_tbl(
  feature_score_base %>%
    dplyr::arrange(dplyr::desc(loading_abs_max)) %>%
    dplyr::select(
      view, view_code, feature_original, feature_label, feature_model,
      dominant_factor, loading_abs_max,
      loading_percentile_view, view_r2_total, factor_r2
    ) %>%
    dplyr::slice_head(n = 20),
  n = 20
)

## ============================================================
## 5) Reconstruir selección de features por TODAS las grids
## ============================================================

.reconstruct_grid_selections <- function(feature_score_base, thr_grid, view_map) {
  
  view_names <- names(view_map)
  out <- list()
  
  for (g in seq_len(nrow(thr_grid))) {
    
    for (v in view_names) {
      
      code <- view_map[[v]]
      n_col <- paste0("n_", code)
      q_col <- paste0("q_", code)
      
      feat_v <- feature_score_base %>%
        dplyr::filter(view == v) %>%
        dplyr::arrange(
          dplyr::desc(loading_abs_max),
          dplyr::desc(loading_percentile_view),
          feature_original
        )
      
      if (n_col %in% colnames(thr_grid)) {
        
        n_keep <- as.integer(thr_grid[[n_col]][g])
        n_keep <- max(1L, min(n_keep, nrow(feat_v)))
        
        selected_v <- feat_v %>%
          dplyr::slice_head(n = n_keep)
        
      } else if (q_col %in% colnames(thr_grid)) {
        
        q_thr <- thr_grid[[q_col]][g]
        thr <- stats::quantile(feat_v$loading_abs_max, probs = q_thr, na.rm = TRUE)
        
        selected_v <- feat_v %>%
          dplyr::filter(loading_abs_max >= thr)
        
      } else {
        stop("No encuentro columna ", n_col, " ni ", q_col, " para vista ", v)
      }
      
      selected_v <- selected_v %>%
        dplyr::mutate(
          GridID = g,
          p_target = if ("p_target" %in% colnames(thr_grid)) thr_grid$p_target[g] else NA_real_,
          p_actual = if ("p_actual" %in% colnames(thr_grid)) thr_grid$p_actual[g] else NA_real_,
          n_selected_view = nrow(selected_v),
          selected_in_grid = TRUE
        )
      
      out[[paste(g, v, sep = "__")]] <- selected_v
    }
  }
  
  dplyr::bind_rows(out)
}

grid_feature_selection_tbl <- .reconstruct_grid_selections(
  feature_score_base = feature_score_base,
  thr_grid = thr_grid,
  view_map = view_map
)

cat("\n--- Auditoría selección reconstruida por grid y vista ---\n")
.print_tbl(
  grid_feature_selection_tbl %>%
    dplyr::count(GridID, p_target, p_actual, view, view_code, name = "n_selected") %>%
    dplyr::arrange(GridID, view),
  n = Inf
)

cat("\n--- Suma seleccionada por grid reconstruida ---\n")
.print_tbl(
  grid_feature_selection_tbl %>%
    dplyr::count(GridID, p_target, p_actual, name = "n_selected_total") %>%
    dplyr::arrange(GridID),
  n = Inf
)

## ============================================================
## 6) Frecuencia de selección y score final
## ============================================================

n_grids_total <- nrow(thr_grid)

selection_stability_tbl <- grid_feature_selection_tbl %>%
  dplyr::group_by(
    view, view_code,
    feature_original,
    feature_prefixed_original,
    feature_model
  ) %>%
  dplyr::summarise(
    n_grids_selected = dplyr::n_distinct(GridID),
    selection_freq_grid = n_grids_selected / n_grids_total,
    first_grid_selected = min(GridID, na.rm = TRUE),
    min_p_selected = min(p_actual, na.rm = TRUE),
    max_p_selected = max(p_actual, na.rm = TRUE),
    .groups = "drop"
  )

mofa_candidate_score_tbl <- feature_score_base %>%
  dplyr::left_join(
    selection_stability_tbl,
    by = c(
      "view", "view_code",
      "feature_original",
      "feature_prefixed_original",
      "feature_model"
    )
  ) %>%
  dplyr::mutate(
    n_grids_selected = dplyr::coalesce(n_grids_selected, 0L),
    selection_freq_grid = dplyr::coalesce(selection_freq_grid, 0),
    first_grid_selected = ifelse(is.na(first_grid_selected), NA_integer_, first_grid_selected),
    min_p_selected = ifelse(is.na(min_p_selected), NA_real_, min_p_selected),
    max_p_selected = ifelse(is.na(max_p_selected), NA_real_, max_p_selected)
  )

## Score MOFA puro:
## - 45% peso relativo dentro de vista
## - 25% importancia del factor dentro de vista
## - 20% importancia global de la vista
## - 10% estabilidad en grids
mofa_candidate_score_tbl <- mofa_candidate_score_tbl %>%
  dplyr::mutate(
    ## Score principal:
    ## ya viene calculado como sum(abs(loading) * R2_factor)
    score_mofa_raw = score_mofa_for_selection,
    
    ## Score secundario descriptivo:
    ## pondera el score MOFA por frecuencia de aparición en grids antiguas por loading.
    ## No se usa para seleccionar el modelo, solo para interpretación.
    score_mofa_stability = score_mofa_for_selection * selection_freq_grid,
    
    rank_mofa_for_selection = rank(-score_mofa_for_selection, ties.method = "first"),
    rank_mofa_raw = rank(-score_mofa_raw, ties.method = "first"),
    rank_mofa_stability = rank(-score_mofa_stability, ties.method = "first")
  ) %>%
  dplyr::arrange(rank_mofa_for_selection)
## ============================================================
## 7) Marcar features presentes en grid final y en BRMS final
## ============================================================

if (exists("grid_id")) {
  final_grid_features <- score_grid_feature_selection_tbl %>%
    dplyr::filter(GridID == grid_id) %>%
    dplyr::pull(feature_model) %>%
    unique()
  
  mofa_candidate_score_tbl <- mofa_candidate_score_tbl %>%
    dplyr::mutate(
      selected_in_final_grid = feature_model %in% final_grid_features
    )
} else {
  mofa_candidate_score_tbl <- mofa_candidate_score_tbl %>%
    dplyr::mutate(
      selected_in_final_grid = NA
    )
}

if (exists("df_train")) {
  brms_final_features <- setdiff(colnames(df_train), "y")
  
  mofa_candidate_score_tbl <- mofa_candidate_score_tbl %>%
    dplyr::mutate(
      present_in_brms_final = feature_model %in% brms_final_features
    )
} else {
  mofa_candidate_score_tbl <- mofa_candidate_score_tbl %>%
    dplyr::mutate(
      present_in_brms_final = NA
    )
}

## ============================================================
## 8) Salidas por pantalla
## ============================================================

cat("\n============================================================\n")
cat("TOP 50 BIOMARCADORES CANDIDATOS POR SCORE MOFA PONDERADO\n")
cat("============================================================\n")

.print_tbl(
  mofa_candidate_score_tbl %>%
    dplyr::select(
      rank_mofa_for_selection,
      rank_mofa_stability,
      view, view_code,
      feature_original, feature_label, feature_model,
      dominant_factor,
      loading_abs_max,
      factor_r2,
      dominant_factor_contribution,
      score_mofa_for_selection,
      score_mofa_percentile_view,
      selection_freq_grid,
      score_mofa_stability,
      selected_in_final_grid,
      present_in_brms_final
    ) %>%
    dplyr::slice_head(n = 50),
  n = 50
)

cat("\n============================================================\n")
cat("TOP POR VISTA\n")
cat("============================================================\n")

.print_tbl(
  mofa_candidate_score_tbl %>%
    dplyr::group_by(view) %>%
    dplyr::slice_min(order_by = rank_mofa_for_selection, n = 15, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(view, rank_mofa_for_selection) %>%
    dplyr::select(
      view,
      rank_mofa_for_selection,
      rank_mofa_stability,
      feature_original,
      feature_label,
      feature_model,
      dominant_factor,
      loading_abs_max,
      factor_r2,
      dominant_factor_contribution,
      score_mofa_for_selection,
      score_mofa_percentile_view,
      selection_freq_grid,
      score_mofa_stability,
      selected_in_final_grid,
      present_in_brms_final
    ),
  n = Inf
)



if (exists("df_train")) {
  cat("\n============================================================\n")
  cat("FEATURES PRESENTES EN BRMS FINAL ORDENADAS POR SCORE MOFA\n")
  cat("============================================================\n")
  
  .print_tbl(
    mofa_candidate_score_tbl %>%
      dplyr::filter(present_in_brms_final %in% TRUE) %>%
      dplyr::arrange(rank_mofa_stability) %>%
      dplyr::select(
        rank_mofa_for_selection,
        rank_mofa_stability,
        view,
        feature_original,
        feature_label,
        feature_model,
        dominant_factor,
        loading_abs_max,
        loading_percentile_view,
        factor_r2,
        score_mofa_for_selection,
        score_mofa_stability,
        selected_in_final_grid
      ),
    n = Inf
  )
}

cat("\n============================================================\n")
cat("RESUMEN DEL SCORE\n")
cat("============================================================\n")

cat("\nNúmero total de features evaluadas:\n")
print(nrow(mofa_candidate_score_tbl))

cat("\nNúmero de features seleccionadas al menos en una grid:\n")
print(sum(mofa_candidate_score_tbl$n_grids_selected > 0))

cat("\nNúmero de features seleccionadas en todas las grids:\n")
print(sum(mofa_candidate_score_tbl$selection_freq_grid == 1))

cat("\nDistribución por vista, candidatas seleccionadas al menos una vez:\n")
.print_tbl(
  mofa_candidate_score_tbl %>%
    dplyr::filter(n_grids_selected > 0) %>%
    dplyr::count(view, view_code, name = "n_features_selected_any_grid") %>%
    dplyr::arrange(view),
  n = Inf
)

cat("\nDistribución por vista, top 50 global:\n")
.print_tbl(
  mofa_candidate_score_tbl %>%
    dplyr::slice_head(n = 50) %>%
    dplyr::count(view, view_code, name = "n_top50") %>%
    dplyr::arrange(view),
  n = Inf
)

cat("\n============================================================\n")
cat("Objetos creados en memoria:\n")
cat(" - feature_score_base\n")
cat(" - grid_feature_selection_tbl\n")
cat(" - selection_stability_tbl\n")
cat(" - mofa_candidate_score_tbl\n")
cat("============================================================\n")
## ============================================================
## Construir return.list usando score MOFA como criterio de selección
## Sustituye la selección previa basada solo en loading_abs_max
## ============================================================
## ============================================================
## Selección mRMR-like:
## - Relevancia = score MOFA ya calculado
## - Redundancia = correlación absoluta en TRAIN
## - No crea score nuevo
## ============================================================

.select_mofa_mrmr_view <- function(score_v,
                                   X_train_view,
                                   n_keep,
                                   score_col = "score_mofa_for_selection",
                                   rho_cutoff = 0.95,
                                   cor_method = "spearman") {
  
  if (nrow(score_v) == 0) {
    stop("score_v vacío en .select_mofa_mrmr_view().")
  }
  
  n_keep <- max(1L, min(as.integer(n_keep), nrow(score_v)))
  
  score_v <- score_v %>%
    dplyr::arrange(
      dplyr::desc(.data[[score_col]]),
      dplyr::desc(score_mofa_percentile_view),
      dplyr::desc(loading_abs_max),
      feature_model
    )
  
  selected_idx <- integer(0)
  skipped_list <- list()
  
  for (ii in seq_len(nrow(score_v))) {
    
    candidate_feature <- score_v$feature_original[ii]
    
    if (!candidate_feature %in% rownames(X_train_view)) {
      next
    }
    
    if (length(selected_idx) == 0) {
      selected_idx <- c(selected_idx, ii)
    } else {
      
      selected_features <- score_v$feature_original[selected_idx]
      selected_features <- selected_features[selected_features %in% rownames(X_train_view)]
      
      x_candidate <- as.numeric(X_train_view[candidate_feature, ])
      X_selected  <- t(X_train_view[selected_features, , drop = FALSE])
      
      rho_vec <- suppressWarnings(
        stats::cor(
          x_candidate,
          X_selected,
          use = "pairwise.complete.obs",
          method = cor_method
        )
      )
      
      rho_vec <- as.numeric(rho_vec)
      names(rho_vec) <- selected_features
      
      if (length(rho_vec) == 0 || all(is.na(rho_vec))) {
        max_abs_rho <- NA_real_
        redundant_with <- NA_character_
      } else {
        max_abs_rho <- max(abs(rho_vec), na.rm = TRUE)
        redundant_with <- names(rho_vec)[which.max(abs(rho_vec))]
      }
      
      if (is.finite(max_abs_rho) && max_abs_rho > rho_cutoff) {
        
        skipped_list[[length(skipped_list) + 1L]] <- score_v[ii, , drop = FALSE] %>%
          dplyr::mutate(
            mrmr_status = "skipped_redundant",
            max_abs_rho_to_selected = max_abs_rho,
            redundant_with_original = redundant_with,
            rho_cutoff = rho_cutoff,
            cor_method = cor_method
          )
        
      } else {
        selected_idx <- c(selected_idx, ii)
      }
    }
    
    if (length(selected_idx) >= n_keep) {
      break
    }
  }
  
  selected_core <- score_v[selected_idx, , drop = FALSE] %>%
    dplyr::mutate(
      mrmr_status = "selected_nonredundant",
      max_abs_rho_to_selected = NA_real_,
      redundant_with_original = NA_character_,
      rho_cutoff = rho_cutoff,
      cor_method = cor_method
    )
  
  ## Si el filtro por redundancia deja menos features que n_keep,
  ## se rellena por ranking MOFA para no romper tamaños de la grilla.
  ## Versión estricta:
  ## Si una feature fue descartada por redundancia, no se reintroduce como relleno.
  ## Por tanto, una vista puede quedar con menos features que n_keep.
  selected_final <- selected_core
  
  
  skipped_tbl <- if (length(skipped_list) > 0) {
    dplyr::bind_rows(skipped_list)
  } else {
    score_v[0, , drop = FALSE] %>%
      dplyr::mutate(
        mrmr_status = character(),
        max_abs_rho_to_selected = numeric(),
        redundant_with_original = character(),
        rho_cutoff = numeric(),
        cor_method = character()
      )
  }
  
  list(
    selected = selected_final,
    skipped = skipped_tbl
  )
}
build_sets_from_score_grid <- function(score_tbl,
                                       thr_grid,
                                       tst_Data,
                                       view_map,
                                       score_col = "score_mofa_for_selection",
                                       refactor_output = TRUE) {
  
  if (!score_col %in% colnames(score_tbl)) {
    stop("No existe score_col = ", score_col, " en score_tbl.")
  }
  
  out_train <- vector("list", nrow(thr_grid))
  out_test  <- vector("list", nrow(thr_grid))
  out_info  <- vector("list", nrow(thr_grid))
  out_selection <- vector("list", nrow(thr_grid))
  out_skipped <- vector("list", nrow(thr_grid))
  
  y_train <- factor(tst_Data$grupo_train$grupo)
  names(y_train) <- tst_Data$grupo_train$id
  
  y_test <- factor(tst_Data$grupo_test$grupo, levels = levels(y_train))
  names(y_test) <- tst_Data$grupo_test$id
  
  view_names <- names(view_map)
  view_codes <- unname(view_map)
  
  for (i in seq_len(nrow(thr_grid))) {
    
    train_list_i <- list()
    test_list_i  <- list()
    info_i       <- list()
    sel_i        <- list()
    skipped_i    <- list()
    
    for (j in seq_along(view_map)) {
      
      view_name <- view_names[j]
      view_code <- view_codes[j]
      
      n_col <- paste0("n_", view_code)
      
      if (!n_col %in% colnames(thr_grid)) {
        stop("No encuentro columna ", n_col, " en thr_grid.")
      }
      
      n_keep_i <- as.integer(thr_grid[[n_col]][i])
      
      if (is.na(n_keep_i) || !is.finite(n_keep_i)) {
        stop(
          "n_keep_i inválido para GridID=", i,
          ", vista=", view_name,
          ", columna=", n_col,
          ", valor=", n_keep_i
        )
      }
      
      score_v <- score_tbl %>%
        dplyr::filter(
          .data$view == .env$view_name,
          .data$view_code == .env$view_code
        ) %>%
        dplyr::arrange(
          dplyr::desc(.data[[score_col]]),
          dplyr::desc(score_mofa_percentile_view),
          dplyr::desc(loading_abs_max),
          feature_model
        )
      
      if (nrow(score_v) == 0) {
        stop("No hay features en score_tbl para vista ", view_name, " / ", view_code)
      }
      
      n_keep_i <- max(1L, min(n_keep_i, nrow(score_v)))
      
      mrmr_i <- .select_mofa_mrmr_view(
        score_v       = score_v,
        X_train_view  = tst_Data[[view_code]]$train,
        n_keep        = n_keep_i,
        score_col     = score_col,
        rho_cutoff    = MRMR_RHO_CUTOFF,
        cor_method    = MRMR_COR_METHOD
      )
      
      selected_v <- mrmr_i$selected %>%
        dplyr::mutate(
          GridID = i,
          p_target = if ("p_target" %in% colnames(thr_grid)) thr_grid$p_target[i] else NA_real_,
          p_actual = if ("p_actual" %in% colnames(thr_grid)) thr_grid$p_actual[i] else NA_real_,
          n_requested_view = n_keep_i,
          n_selected_view = nrow(mrmr_i$selected),
          selected_by_score = TRUE,
          selected_by_mrmr = TRUE,
          mrmr_strict = TRUE,
          score_selection_column = score_col
        )
      
      skipped_v <- mrmr_i$skipped %>%
        dplyr::mutate(
          GridID = i,
          p_target = if ("p_target" %in% colnames(thr_grid)) thr_grid$p_target[i] else NA_real_,
          p_actual = if ("p_actual" %in% colnames(thr_grid)) thr_grid$p_actual[i] else NA_real_,
          n_requested_view = n_keep_i,
          n_selected_view = nrow(mrmr_i$selected),
          selected_by_score = FALSE,
          selected_by_mrmr = FALSE,
          mrmr_strict = TRUE,
          score_selection_column = score_col
        )
      
      feats_original <- selected_v$feature_original
      feats_model    <- selected_v$feature_model
      
      missing_train <- setdiff(feats_original, rownames(tst_Data[[view_code]]$train))
      missing_test  <- setdiff(feats_original, rownames(tst_Data[[view_code]]$test))
      
      if (length(missing_train) > 0) {
        stop(
          "Features no encontradas en train para vista ", view_name, ": ",
          paste(missing_train, collapse = ", ")
        )
      }
      
      if (length(missing_test) > 0) {
        stop(
          "Features no encontradas en test para vista ", view_name, ": ",
          paste(missing_test, collapse = ", ")
        )
      }
      
      sub_train <- tst_Data[[view_code]]$train[feats_original, , drop = FALSE]
      sub_test  <- tst_Data[[view_code]]$test[feats_original,  , drop = FALSE]
      
      ## Usar directamente el nombre sanitizado del modelo.
      ## Así el BRMS y el score quedan alineados.
      rownames(sub_train) <- feats_model
      rownames(sub_test)  <- feats_model
      
      train_list_i[[view_code]] <- sub_train
      test_list_i[[view_code]]  <- sub_test
      
      info_i[[view_code]] <- data.frame(
        view_name = view_name,
        view_code = view_code,
        n_selected = nrow(sub_train),
        score_col = score_col,
        stringsAsFactors = FALSE
      )
      
      sel_i[[view_code]] <- selected_v
      skipped_i[[view_code]] <- skipped_v
      
    }
    
    train_feats <- do.call(rbind, train_list_i)
    test_feats  <- do.call(rbind, test_list_i)
    
    if (anyDuplicated(rownames(train_feats))) {
      stop(
        "Hay nombres duplicados tras selección por score en GridID=", i,
        ": ",
        paste(unique(rownames(train_feats)[duplicated(rownames(train_feats))]), collapse = ", ")
      )
    }
    
    train_feats <- train_feats[, tst_Data$grupo_train$id, drop = FALSE]
    test_feats  <- test_feats[,  tst_Data$grupo_test$id,  drop = FALSE]
    
    out_train[[i]] <- as.data.frame(t(train_feats), check.names = FALSE)
    out_test[[i]]  <- as.data.frame(t(test_feats),  check.names = FALSE)
    
    out_info[[i]] <- dplyr::bind_rows(info_i) %>%
      dplyr::mutate(GridID = i)
    
    out_selection[[i]] <- dplyr::bind_rows(sel_i)
    out_skipped[[i]] <- dplyr::bind_rows(skipped_i)
    
  }
  
  list(
    train = out_train,
    test = out_test,
    y_train = y_train[tst_Data$grupo_train$id],
    y_test  = y_test[tst_Data$grupo_test$id],
    info = dplyr::bind_rows(out_info),
    score_selection = dplyr::bind_rows(out_selection),
    mrmr_skipped = dplyr::bind_rows(out_skipped)
  )
}

## ------------------------------------------------------------
## Construir return.list final usando score MOFA ponderado
## ------------------------------------------------------------

return.list_score <- build_sets_from_score_grid(
  score_tbl = mofa_candidate_score_tbl,
  thr_grid  = thr_grid,
  tst_Data  = tst_Data,
  view_map  = view_map,
  score_col = "score_mofa_for_selection"
)
# grid_size_after_mrmr_tbl <- score_grid_feature_selection_tbl %>%
#   dplyr::count(GridID, name = "p_after_mrmr") %>%
#   dplyr::left_join(
#     thr_grid %>% dplyr::mutate(GridID = dplyr::row_number()),
#     by = "GridID"
#   ) %>%
#   dplyr::select(GridID, p_target, p_actual_before_mrmr = p_actual, p_after_mrmr)
# 
# print(grid_size_after_mrmr_tbl)
# 
# write.csv(
#   grid_size_after_mrmr_tbl,
#   out_table("grid_size_after_mrmr_tbl.csv"),
#   row.names = FALSE
# )
return.list <- return.list_score

score_grid_feature_selection_tbl <- return.list_score$score_selection
mrmr_skipped_features_tbl <- return.list_score$mrmr_skipped

## ============================================================
## AUDITORÍA TAMAÑO DE GRID ANTES/DESPUÉS DE mRMR
## Guarda tabla y plot en las carpetas del análisis
## ============================================================

grid_size_after_mrmr_tbl <- score_grid_feature_selection_tbl %>%
  dplyr::count(GridID, name = "p_after_mrmr") %>%
  dplyr::left_join(
    thr_grid %>%
      dplyr::mutate(GridID = dplyr::row_number()),
    by = "GridID"
  ) %>%
  dplyr::mutate(
    p_lost_by_mrmr = p_actual - p_after_mrmr,
    prop_kept_after_mrmr = p_after_mrmr / p_actual
  ) %>%
  dplyr::select(
    GridID,
    p_target,
    p_actual_before_mrmr = p_actual,
    p_after_mrmr,
    p_lost_by_mrmr,
    prop_kept_after_mrmr
  )

cat("\n--- Tamaño de grid antes/después de mRMR ---\n")
print(grid_size_after_mrmr_tbl)

write.csv(
  grid_size_after_mrmr_tbl,
  out_table("grid_size_after_mrmr_tbl.csv"),
  row.names = FALSE
)

## Carpeta de plots mRMR, igual que el resto de auditorías mRMR
MRMR_PLOTS_DIR <- out_plot("mrmr")
dir.create(MRMR_PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

grid_size_plot_long <- grid_size_after_mrmr_tbl %>%
  tidyr::pivot_longer(
    cols = c(p_actual_before_mrmr, p_after_mrmr),
    names_to = "stage",
    values_to = "n_features"
  ) %>%
  dplyr::mutate(
    stage = dplyr::recode(
      stage,
      p_actual_before_mrmr = "Antes de mRMR",
      p_after_mrmr = "Después de mRMR"
    )
  )

gg_grid_size_after_mrmr <- ggplot(
  grid_size_plot_long,
  aes(x = factor(GridID), y = n_features, fill = stage)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  geom_text(
    aes(label = n_features),
    position = position_dodge(width = 0.75),
    vjust = -0.4,
    size = 3
  ) +
  geom_text(
    data = grid_size_after_mrmr_tbl,
    aes(
      x = factor(GridID),
      y = p_actual_before_mrmr + 4,
      label = paste0("-", p_lost_by_mrmr, " vars")
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  labs(
    title = "Tamaño de las grids antes y después de mRMR",
    subtitle = paste0(
      "mRMR estricto: no se rellenan variables eliminadas por redundancia | rho > ",
      MRMR_RHO_CUTOFF
    ),
    x = "GridID",
    y = "Número de variables",
    fill = NULL
  ) +
  theme_minimal()

print(gg_grid_size_after_mrmr)

ggsave(
  filename = file.path(MRMR_PLOTS_DIR, "00_grid_size_after_mrmr.png"),
  plot = gg_grid_size_after_mrmr,
  width = 9,
  height = 6,
  dpi = 300
)

## Opcional: guardar también PDF
ggsave(
  filename = file.path(MRMR_PLOTS_DIR, "00_grid_size_after_mrmr.pdf"),
  plot = gg_grid_size_after_mrmr,
  width = 9,
  height = 6
)



mrmr_selection_audit_tbl <- score_grid_feature_selection_tbl %>%
  dplyr::count(
    GridID,
    p_target,
    p_actual,
    view,
    view_code,
    mrmr_status,
    name = "n_features"
  ) %>%
  dplyr::arrange(GridID, view, mrmr_status)

## ------------------------------------------------------------
## Auditoría robusta: redundancia mRMR por grid y vista
## ------------------------------------------------------------
## Si mRMR no descarta ninguna feature, mrmr_skipped_features_tbl
## queda con 0 filas. En ese caso, se fuerza una tabla completa
## GridID × vista con n_skipped_redundant = 0.
##
## Esto NO cambia la selección, NO cambia el score MOFA,
## NO cambia return.list y NO cambia BRMS.
## Solo hace robusta la tabla/figura de auditoría.
## ------------------------------------------------------------

mrmr_redundancy_summary_raw <- mrmr_skipped_features_tbl %>%
  dplyr::count(
    GridID,
    p_target,
    p_actual,
    view,
    view_code,
    name = "n_skipped_redundant"
  )

mrmr_grid_view_template <- tidyr::expand_grid(
  GridID = seq_len(nrow(thr_grid)),
  view   = names(view_map)
) %>%
  dplyr::mutate(
    view_code = unname(view_map[view])
  ) %>%
  dplyr::left_join(
    thr_grid %>%
      dplyr::mutate(GridID = dplyr::row_number()) %>%
      dplyr::select(GridID, p_target, p_actual),
    by = "GridID"
  ) %>%
  dplyr::select(
    GridID,
    p_target,
    p_actual,
    view,
    view_code
  )

mrmr_redundancy_summary_tbl <- mrmr_grid_view_template %>%
  dplyr::left_join(
    mrmr_redundancy_summary_raw,
    by = c("GridID", "p_target", "p_actual", "view", "view_code")
  ) %>%
  dplyr::mutate(
    n_skipped_redundant = dplyr::coalesce(n_skipped_redundant, 0L)
  ) %>%
  dplyr::arrange(GridID, view)

mrmr_strict_size_audit_tbl <- score_grid_feature_selection_tbl %>%
  dplyr::group_by(
    GridID,
    p_target,
    p_actual,
    view,
    view_code
  ) %>%
  dplyr::summarise(
    n_requested_view = dplyr::first(n_requested_view),
    n_selected_view = dplyr::n(),
    n_lost_by_mrmr = n_requested_view - n_selected_view,
    .groups = "drop"
  ) %>%
  dplyr::arrange(GridID, view)


cat("\n============================================================\n")
cat("AUDITORÍA return.list BASADO EN SCORE MOFA PONDERADO\n")
cat("============================================================\n")

cat("\n--- Conteo por grid y vista usando score MOFA ponderado ---\n")

cat("\n--- Auditoría mRMR: seleccionadas por estado ---\n")
.print_tbl(mrmr_selection_audit_tbl, n = Inf)

cat("\n--- Auditoría mRMR: features descartadas por redundancia ---\n")
.print_tbl(mrmr_redundancy_summary_tbl, n = Inf)

cat("\n--- Auditoría mRMR estricto: variables solicitadas vs seleccionadas ---\n")
.print_tbl(mrmr_strict_size_audit_tbl, n = Inf)



cat("\n--- Top features descartadas por redundancia mRMR ---\n")
.print_tbl(
  mrmr_skipped_features_tbl %>%
    dplyr::arrange(
      GridID,
      view,
      dplyr::desc(max_abs_rho_to_selected),
      dplyr::desc(score_mofa_for_selection)
    ) %>%
    dplyr::select(
      GridID,
      view,
      view_code,
      feature_original,
      feature_label,
      feature_model,
      score_mofa_for_selection,
      max_abs_rho_to_selected,
      redundant_with_original,
      rho_cutoff,
      cor_method
    ) %>%
    dplyr::slice_head(n = 50),
  n = 50
)


.print_tbl(
  score_grid_feature_selection_tbl %>%
    dplyr::count(GridID, p_target, p_actual, view, view_code, name = "n_selected") %>%
    dplyr::arrange(GridID, view),
  n = Inf
)

cat("\n--- Suma por grid usando score MOFA ponderado ---\n")
.print_tbl(
  score_grid_feature_selection_tbl %>%
    dplyr::count(GridID, p_target, p_actual, name = "n_selected_total") %>%
    dplyr::arrange(GridID),
  n = Inf
)

cat("\n--- Top 30 features seleccionadas por score MOFA ponderado en GridID 1 ---\n")
.print_tbl(
  score_grid_feature_selection_tbl %>%
    dplyr::filter(GridID == 1) %>%
    dplyr::arrange(dplyr::desc(score_mofa_for_selection)) %>%
    dplyr::select(
      GridID, view, feature_original, feature_label, feature_model,
      dominant_factor,
      loading_abs_max,
      factor_r2,
      dominant_factor_contribution,
      score_mofa_for_selection,
      score_mofa_percentile_view
    ) %>%
    dplyr::slice_head(n = 30),
  n = 30
)

cat("\n--- Diferencias contra selección previa por loading_abs_max ---\n")

diff_score_vs_loading <- dplyr::bind_rows(lapply(seq_len(nrow(thr_grid)), function(g) {
  
  old_features <- grid_feature_selection_tbl %>%
    dplyr::filter(GridID == g) %>%
    dplyr::pull(feature_model) %>%
    unique()
  
  new_features <- score_grid_feature_selection_tbl %>%
    dplyr::filter(GridID == g) %>%
    dplyr::pull(feature_model) %>%
    unique()
  
  data.frame(
    GridID = g,
    p_actual = thr_grid$p_actual[g],
    n_old = length(old_features),
    n_new = length(new_features),
    n_common = length(intersect(old_features, new_features)),
    n_only_old = length(setdiff(old_features, new_features)),
    n_only_new = length(setdiff(new_features, old_features)),
    prop_common = length(intersect(old_features, new_features)) / length(new_features)
  )
}))

.print_tbl(diff_score_vs_loading, n = Inf)

write.csv(
  score_grid_feature_selection_tbl,
  out_table("score_grid_feature_selection_tbl.csv"),
  row.names = FALSE
)
write.csv(
  mrmr_selection_audit_tbl,
  out_table("mrmr_selection_audit_tbl.csv"),
  row.names = FALSE
)

write.csv(
  mrmr_skipped_features_tbl,
  out_table("mrmr_skipped_features_tbl.csv"),
  row.names = FALSE
)

write.csv(
  mrmr_redundancy_summary_tbl,
  out_table("mrmr_redundancy_summary_tbl.csv"),
  row.names = FALSE
)

write.csv(
  mrmr_strict_size_audit_tbl,
  out_table("mrmr_strict_size_audit_tbl.csv"),
  row.names = FALSE
)


write.csv(
  diff_score_vs_loading,
  out_table("diff_score_vs_loading_grid.csv"),
  row.names = FALSE
)
write.csv(
  feature_score_base,
  out_table("feature_score_base.csv"),
  row.names = FALSE
)

write.csv(
  grid_feature_selection_tbl,
  out_table("grid_feature_selection_tbl_loading_based.csv"),
  row.names = FALSE
)

write.csv(
  selection_stability_tbl,
  out_table("selection_stability_tbl.csv"),
  row.names = FALSE
)

write.csv(
  mofa_candidate_score_tbl,
  out_table("mofa_candidate_score_tbl.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    feature_score_base = feature_score_base,
    grid_feature_selection_tbl = grid_feature_selection_tbl,
    selection_stability_tbl = selection_stability_tbl,
    mofa_candidate_score_tbl = mofa_candidate_score_tbl,
    score_grid_feature_selection_tbl = score_grid_feature_selection_tbl,
    mrmr_selection_audit_tbl = mrmr_selection_audit_tbl,
    mrmr_skipped_features_tbl = mrmr_skipped_features_tbl,
    mrmr_redundancy_summary_tbl = mrmr_redundancy_summary_tbl,
    mrmr_strict_size_audit_tbl = mrmr_strict_size_audit_tbl,
    diff_score_vs_loading = diff_score_vs_loading
  ),
  out_rds("mofa_score_selection_objects.rds")
)
MRMR_PLOTS_DIR <- out_plot("mrmr")
dir.create(MRMR_PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

gg_mrmr_status <- mrmr_selection_audit_tbl %>%
  ggplot(aes(x = factor(GridID), y = n_features, fill = mrmr_status)) +
  geom_col(position = "stack") +
  facet_wrap(~ view, scales = "free_y") +
  labs(
    title = "Auditoría mRMR: features seleccionadas por estado",
    subtitle = paste0("Redundancia definida como |rho| > ", MRMR_RHO_CUTOFF),
    x = "GridID",
    y = "Número de features",
    fill = "Estado mRMR"
  ) +
  theme_minimal()

print(gg_mrmr_status)

ggsave(
  filename = file.path(MRMR_PLOTS_DIR, "01_mrmr_estado_seleccion_por_grid_vista.png"),
  plot = gg_mrmr_status,
  width = 10,
  height = 7,
  dpi = 300
)

gg_mrmr_skipped <- mrmr_redundancy_summary_tbl %>%
  ggplot(aes(x = factor(GridID), y = n_skipped_redundant)) +
  geom_col() +
  facet_wrap(~ view, scales = "free_y") +
  labs(
    title = "Auditoría mRMR: features descartadas por redundancia",
    subtitle = paste0("Descartadas cuando |rho| > ", MRMR_RHO_CUTOFF),
    x = "GridID",
    y = "Features descartadas"
  ) +
  theme_minimal()

print(gg_mrmr_skipped)

ggsave(
  filename = file.path(MRMR_PLOTS_DIR, "02_mrmr_descartadas_por_grid_vista.png"),
  plot = gg_mrmr_skipped,
  width = 10,
  height = 7,
  dpi = 300
)


gg_mrmr_lost <- mrmr_strict_size_audit_tbl %>%
  ggplot(aes(x = factor(GridID), y = n_lost_by_mrmr)) +
  geom_col() +
  facet_wrap(~ view, scales = "free_y") +
  labs(
    title = "Auditoría mRMR estricto: variables eliminadas por redundancia",
    subtitle = paste0("No se rellena después de aplicar |rho| > ", MRMR_RHO_CUTOFF),
    x = "GridID",
    y = "Variables eliminadas"
  ) +
  theme_minimal()

print(gg_mrmr_lost)

ggsave(
  filename = file.path(MRMR_PLOTS_DIR, "03_mrmr_variables_perdidas_por_grid_vista.png"),
  plot = gg_mrmr_lost,
  width = 10,
  height = 7,
  dpi = 300
)






if (nrow(mrmr_skipped_features_tbl) > 0) {
  
  gg_mrmr_rho <- mrmr_skipped_features_tbl %>%
    ggplot(aes(x = max_abs_rho_to_selected)) +
    geom_histogram(bins = 20) +
    facet_wrap(~ view, scales = "free_y") +
    geom_vline(xintercept = MRMR_RHO_CUTOFF, linetype = 2) +
    labs(
      title = "Distribución de |rho| en features descartadas por mRMR",
      x = "|rho| frente a feature ya seleccionada",
      y = "Número de features"
    ) +
    theme_minimal()
  
  print(gg_mrmr_rho)
  
  ggsave(
    filename = file.path(MRMR_PLOTS_DIR, "04_mrmr_distribucion_rho_descartadas.png"),
    plot = gg_mrmr_rho,
    width = 9,
    height = 6,
    dpi = 300
  )
}

detach("package:MOFA2", unload = TRUE)
## ============================================================
## Selección del grid por CV interna estratificada
## Train: NP=6, SO=7, FD=3 -> K = 3
## Test queda intacto para evaluación final
## ============================================================

ref <- "NP"
backend <- .detect_backend_safe()

K_cv <- 3
ncores_avail <- parallel::detectCores()

if (is.na(ncores_avail)) {
  ncores_avail <- 2
}

cores_cv <- min(3, max(1, ncores_avail - 1))



## ============================================================
## DIAGNÓSTICO SERIAL DE COMPILACIÓN
## GridID 1 / Fold 1
##
## No ejecuta la selección de grid.
## No ejecuta CV repetida.
## No ejecuta TEST.
## No modifica el modelo original.
## ============================================================

options(
  warn = 1,
  width = 200
)

DIAG_DIR <- Sys.getenv(
  "DIAG_DIR",
  unset = "diagnostico_compilacion_modelamiento2_q00"
)

dir.create(
  DIAG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n")
cat("============================================================\n")
cat("DIAGNÓSTICO SERIAL GRID 1 / FOLD 1\n")
cat("============================================================\n")

cat("\nObjetos disponibles:\n")
cat("return.list:", exists("return.list"), "\n")
cat("thr_grid:", exists("thr_grid"), "\n")
cat(".prepare_cv_data:", exists(".prepare_cv_data"), "\n")
cat("make_multinom_priors:", exists("make_multinom_priors"), "\n")

stopifnot(
  exists("return.list"),
  exists("thr_grid"),
  exists(".prepare_cv_data"),
  exists("make_multinom_priors")
)

ref_diag <- "NP"
grid_diag <- 1L
fold_diag <- 1L
K_diag <- 3L
seed_diag <- 123L

y_train_diag <- stats::relevel(
  refactor_levels(return.list$y_train),
  ref = ref_diag
)

folds_diag <- .make_stratified_folds(
  y = y_train_diag,
  K = K_diag,
  seed = seed_diag
)

valid_idx_diag <- folds_diag[[fold_diag]]

train_idx_diag <- setdiff(
  seq_along(return.list$y_train),
  valid_idx_diag
)

cat("\n--- Distribución completa ---\n")
print(table(y_train_diag))

cat("\n--- Índices de train del fold ---\n")
print(train_idx_diag)

cat("\n--- Índices de validación del fold ---\n")
print(valid_idx_diag)

cat("\n--- Distribución del fold de entrenamiento ---\n")
print(table(y_train_diag[train_idx_diag]))

cat("\n--- Distribución del fold de validación ---\n")
print(table(y_train_diag[valid_idx_diag]))

dat_diag <- .prepare_cv_data(
  return.list = return.list,
  grid_id = grid_diag,
  train_idx = train_idx_diag,
  valid_idx = valid_idx_diag,
  ref = ref_diag
)

df_train_diag <- dat_diag$df_train
df_valid_diag <- dat_diag$df_valid

p_diag <- ncol(df_train_diag) - 1L

cat("\n--- Dimensiones ---\n")
cat("GridID:", grid_diag, "\n")
cat("Fold:", fold_diag, "\n")
cat("n train:", nrow(df_train_diag), "\n")
cat("n validación:", nrow(df_valid_diag), "\n")
cat("p:", p_diag, "\n")

cat("\n--- Clases en train ---\n")
print(table(df_train_diag$y))

cat("\n--- Clases en validación ---\n")
print(table(df_valid_diag$y))

feature_names_diag <- setdiff(
  colnames(df_train_diag),
  "y"
)

feature_audit_diag <- data.frame(
  feature = feature_names_diag,
  class = vapply(
    df_train_diag[, feature_names_diag, drop = FALSE],
    function(x) paste(class(x), collapse = "/"),
    character(1)
  ),
  n_na_train = vapply(
    df_train_diag[, feature_names_diag, drop = FALSE],
    function(x) sum(is.na(x)),
    integer(1)
  ),
  n_nonfinite_train = vapply(
    df_train_diag[, feature_names_diag, drop = FALSE],
    function(x) {
      if (!is.numeric(x)) return(NA_integer_)
      sum(!is.finite(x))
    },
    integer(1)
  ),
  sd_train = vapply(
    df_train_diag[, feature_names_diag, drop = FALSE],
    function(x) {
      if (!is.numeric(x)) return(NA_real_)
      stats::sd(x, na.rm = TRUE)
    },
    numeric(1)
  ),
  n_unique_train = vapply(
    df_train_diag[, feature_names_diag, drop = FALSE],
    function(x) length(unique(x)),
    integer(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  feature_audit_diag,
  file.path(
    DIAG_DIR,
    "auditoria_features_grid1_fold1.csv"
  ),
  row.names = FALSE
)

cat("\n--- Auditoría de features ---\n")
print(feature_audit_diag, row.names = FALSE)

cat("\nFeatures no numéricas:\n")
print(
  feature_audit_diag$feature[
    !feature_audit_diag$class %in% c(
      "numeric",
      "integer"
    )
  ]
)

cat("\nFeatures con NA:\n")
print(
  feature_audit_diag$feature[
    feature_audit_diag$n_na_train > 0
  ]
)

cat("\nFeatures con valores no finitos:\n")
print(
  feature_audit_diag$feature[
    !is.na(feature_audit_diag$n_nonfinite_train) &
      feature_audit_diag$n_nonfinite_train > 0
  ]
)

cat("\nFeatures constantes o sin desviación estándar finita:\n")
print(
  feature_audit_diag$feature[
    !is.finite(feature_audit_diag$sd_train) |
      feature_audit_diag$sd_train == 0
  ]
)

backend_diag <- .detect_backend_safe()
ctrl_diag <- .make_ctrl_safe(backend_diag)

cat("\n============================================================\n")
cat("BACKEND DETECTADO:", backend_diag, "\n")
cat("============================================================\n")

if (identical(backend_diag, "cmdstanr")) {
  cat("\nCmdStan path:\n")
  print(cmdstanr::cmdstan_path())

  cat("\nCmdStan version:\n")
  print(cmdstanr::cmdstan_version())
}

priors_diag <- make_multinom_priors(
  y_levels = levels(df_train_diag$y),
  y = df_train_diag$y,
  feature_names = feature_names_diag,
  ref = ref_diag,
  view_map = .get_global_or_null("view_map"),
  r2_score = .get_global_or_null("r2_score")
)

cat("\n--- Priors ---\n")
print(priors_diag)

capture.output(
  priors_diag,
  file = file.path(
    DIAG_DIR,
    "priors_grid1_fold1.txt"
  )
)

## ------------------------------------------------------------
## Generar código Stan y standata antes de compilar
## ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENERANDO CÓDIGO STAN\n")
cat("============================================================\n")

stancode_diag <- brms::make_stancode(
  formula = y ~ .,
  data = df_train_diag,
  family = brms::categorical(),
  prior = priors_diag
)

writeLines(
  stancode_diag,
  con = file.path(
    DIAG_DIR,
    "modelo_grid1_fold1.stan"
  )
)

standata_diag <- brms::make_standata(
  formula = y ~ .,
  data = df_train_diag,
  family = brms::categorical(),
  prior = priors_diag
)

saveRDS(
  standata_diag,
  file.path(
    DIAG_DIR,
    "standata_grid1_fold1.rds"
  )
)

cat("\nCódigo Stan guardado correctamente.\n")

## ------------------------------------------------------------
## Prueba mínima del compilador
## Determina si falla cualquier modelo Stan o solo el BRMS.
## ------------------------------------------------------------

cat("\n============================================================\n")
cat("PRUEBA MÍNIMA DEL COMPILADOR STAN\n")
cat("============================================================\n")

smoke_stan_file <- file.path(
  DIAG_DIR,
  "stan_smoke_test.stan"
)

writeLines(
  c(
    "parameters {",
    "  real y;",
    "}",
    "model {",
    "  y ~ normal(0, 1);",
    "}"
  ),
  smoke_stan_file
)

smoke_error <- NULL

smoke_ok <- tryCatch(
  {
    if (identical(backend_diag, "cmdstanr")) {
      smoke_model <- cmdstanr::cmdstan_model(
        stan_file = smoke_stan_file,
        quiet = FALSE,
        force_recompile = TRUE
      )
    } else {
      smoke_model <- rstan::stan_model(
        file = smoke_stan_file,
        verbose = TRUE,
        auto_write = FALSE
      )
    }

    TRUE
  },
  error = function(e) {
    smoke_error <<- conditionMessage(e)
    FALSE
  }
)

cat("\nResultado prueba mínima:\n")
cat("smoke_ok =", smoke_ok, "\n")

if (!is.null(smoke_error)) {
  cat("smoke_error =", smoke_error, "\n")

  writeLines(
    smoke_error,
    file.path(
      DIAG_DIR,
      "ERROR_stan_smoke_test.txt"
    )
  )
}

## ------------------------------------------------------------
## Compilación directa del modelo BRMS real
##
## Se utiliza el mismo:
## - GridID
## - Fold
## - backend
## - fórmula
## - familia
## - priors
## - control
## - semilla efectiva
##
## Se reducen las iteraciones porque este script solo diagnostica
## la compilación. El código Stan compilado no depende del número
## de iteraciones.
## ------------------------------------------------------------

cat("\n============================================================\n")
cat("INICIO DE COMPILACIÓN BRMS REAL\n")
cat("GRID 1 / FOLD 1 / PROCESO PRINCIPAL / SERIAL\n")
cat("============================================================\n")

flush.console()

fit_diag <- tryCatch(
  brms::brm(
    y ~ .,
    data = df_train_diag,
    family = brms::categorical(),
    prior = priors_diag,
    chains = 1,
    iter = 100,
    warmup = 50,
    refresh = 1,
    backend = backend_diag,
    control = ctrl_diag,
    cores = 1,
    seed = seed_diag + 1000L * grid_diag + fold_diag
  ),
  error = function(e) {
    cat("\n")
    cat("============================================================\n")
    cat("ERROR DIRECTO DE BRMS\n")
    cat("============================================================\n")
    cat(conditionMessage(e), "\n")

    writeLines(
      conditionMessage(e),
      file.path(
        DIAG_DIR,
        "ERROR_brms_grid1_fold1.txt"
      )
    )

    NULL
  }
)

if (is.null(fit_diag)) {
  cat("\nDIAGNÓSTICO TERMINADO CON ERROR DE COMPILACIÓN/AJUSTE.\n")
  quit(
    save = "no",
    status = 42L,
    runLast = FALSE
  )
}

saveRDS(
  fit_diag,
  file.path(
    DIAG_DIR,
    "fit_diagnostico_grid1_fold1.rds"
  )
)

cat("\n")
cat("============================================================\n")
cat("COMPILACIÓN BRMS COMPLETADA CORRECTAMENTE EN SERIAL\n")
cat("============================================================\n")
cat(
  "Esto indicaría que el fallo original depende probablemente ",
  "de la ejecución paralela de varias compilaciones.\n",
  sep = ""
)

quit(
  save = "no",
  status = 0L,
  runLast = FALSE
)


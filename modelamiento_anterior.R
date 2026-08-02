

rm(list = ls())

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
## CONFIGURACIÓN PIPELINE ACTUAL - MODELO ANTERIOR
## Guarda todo en ANALYSIS_DIR/modelado_anterior
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

OUTDIR <- file.path(
  ANALYSIS_DIR,
  "modelado_anterior"
)

PLOTS_DIR <- file.path(OUTDIR, "plots")
TABLES_DIR <- file.path(OUTDIR, "tables")
RDS_DIR <- file.path(OUTDIR, "rds")

dir.create(OUTDIR,     recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR,    recursive = TRUE, showWarnings = FALSE)

out_file <- function(...) file.path(OUTDIR, ...)
out_plot <- function(...) file.path(PLOTS_DIR, ...)
out_table <- function(...) file.path(TABLES_DIR, ...)
out_rds <- function(...) file.path(RDS_DIR, ...)

cat("\n============================================================\n")
cat("MODELO ANTERIOR ADAPTADO AL PIPELINE ACTUAL\n")
cat("ANÁLISIS:", ANALYSIS_NAME, "\n")
cat("ANALYSIS_DIR:", ANALYSIS_DIR, "\n")
cat("RFM_FILE:", RFM_FILE, "\n")
cat("MOFA_MODEL_FILE:", MOFA_MODEL_FILE, "\n")
cat("OUTDIR anterior:", OUTDIR, "\n")
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


select_features_view <- function(mod, view, q_thr, data_list, metadata = NULL){
  aux <- plot_top_weights(mod, view=view, factors=1:2, abs=TRUE, scale=TRUE, nfeatures=Inf)
  thr <- stats::quantile(aux$data$value, probs=q_thr, na.rm=TRUE)
  feats <- aux$data$feature[aux$data$value >= thr]
  sub_train <- data_list$train[as.character(feats), , drop=FALSE]
  sub_test  <- data_list$test[as.character(feats),  , drop=FALSE]
  if(!is.null(metadata) && all(c("EntrezGeneID","GeneSymbol") %in% colnames(metadata))){
    map_dic <- metadata %>% dplyr::select(EntrezGeneID, GeneSymbol) %>% distinct()
    idx <- rownames(sub_train)
    ids_map <- data.frame(EntrezGeneID = idx) %>% dplyr::left_join(map_dic, by = "EntrezGeneID")
    new_ids <- ifelse(!is.na(ids_map$GeneSymbol), ids_map$GeneSymbol, ids_map$EntrezGeneID)
    rownames(sub_train) <- new_ids; rownames(sub_test) <- new_ids
  }
  list(train=sub_train, test=sub_test)
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

build_sets_from_grid <- function(thr_grid, tst_Data, pesos_list, Z_seed, y_seed, mod){
  
  out_train <- vector("list", nrow(thr_grid))
  out_test  <- vector("list", nrow(thr_grid))
  
  ids_train <- tst_Data$grupo_train$id
  ids_test  <- tst_Data$grupo_test$id
  
  X_list <- list(
    transcriptomica = tst_Data$tx$test,
    proteomica      = tst_Data$pr$test,
    metabolomica    = tst_Data$me$test,
    clinical        = tst_Data$cl$test
  )
  
  factores_Test <- project_all_views(X_list, pesos_list)
  
  if (!all(ids_test %in% rownames(factores_Test))) {
    stop("No todos los ids_test están en factores_Test proyectado.")
  }
  
  factores_Test <- factores_Test[ids_test, , drop = FALSE]
  colnames(factores_Test) <- paste0("Factor", seq_len(ncol(factores_Test)))
  
  if (!all(ids_train %in% rownames(Z_seed))) {
    stop("No todos los ids_train están en Z_seed.")
  }
  
  Z_seed <- Z_seed[ids_train, , drop = FALSE]
  colnames(Z_seed) <- paste0("Factor", seq_len(ncol(Z_seed)))
  
  for(i in seq_len(nrow(thr_grid))){
    
    tx <- select_features_view(
      mod       = mod,
      view      = "transcriptomica",
      q_thr     = thr_grid$q_tx[i],
      data_list = tst_Data$tx,
      metadata  = tst_Data$features_metadata
    )
    
    pr <- select_features_view(
      mod       = mod,
      view      = "proteomica",
      q_thr     = thr_grid$q_pr[i],
      data_list = tst_Data$pr,
      metadata  = NULL
    )
    
    me <- select_features_view(
      mod       = mod,
      view      = "metabolomica",
      q_thr     = thr_grid$q_me[i],
      data_list = tst_Data$me,
      metadata  = NULL
    )
    
    cl <- select_features_view(
      mod       = mod,
      view      = "clinical",
      q_thr     = thr_grid$q_cl[i],
      data_list = tst_Data$cl,
      metadata  = NULL
    )
    ## ------------------------------------------------------------
    ## Compatibilidad downstream:
    ## añadir prefijo de vista antes de rbind
    ## Luego sanitize_names() convertirá tx__GENE -> tx_GENE
    ## ------------------------------------------------------------
    
    rownames(tx$train) <- paste0("tx__", rownames(tx$train))
    rownames(tx$test)  <- paste0("tx__", rownames(tx$test))
    
    rownames(pr$train) <- paste0("pr__", rownames(pr$train))
    rownames(pr$test)  <- paste0("pr__", rownames(pr$test))
    
    rownames(me$train) <- paste0("me__", rownames(me$train))
    rownames(me$test)  <- paste0("me__", rownames(me$test))
    
    rownames(cl$train) <- paste0("cl__", rownames(cl$train))
    rownames(cl$test)  <- paste0("cl__", rownames(cl$test))
    train_feats <- rbind(tx$train, pr$train, me$train, cl$train)
    test_feats  <- rbind(tx$test,  pr$test,  me$test,  cl$test)
    
    train_feats <- train_feats[, ids_train, drop = FALSE]
    test_feats  <- test_feats[,  ids_test,  drop = FALSE]
    
    out_train[[i]] <- cbind(Z_seed, t(train_feats)) %>%
      as.data.frame(check.names = FALSE)
    
    out_test[[i]] <- cbind(factores_Test, t(test_feats)) %>%
      as.data.frame(check.names = FALSE)
  }
  
  y_test <- factor(tst_Data$grupo_test$grupo, levels = levels(y_seed))
  names(y_test) <- ids_test
  
  list(
    train   = out_train,
    test    = out_test,
    y_train = y_seed,
    y_test  = y_test
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

make_multinom_priors <- function(y_levels, ref="NP", sd_b=0.3, sd_int=3){
  levs <- setdiff(y_levels, ref)
  do.call(c, lapply(levs, function(lev){
    dpar <- paste0("mu", lev)
    c(set_prior(paste0("student_t(3,0,", sd_b, ")"), class="b", dpar=dpar),
      set_prior(paste0("student_t(3,0,", sd_int,")"), class="Intercept", dpar=dpar))
  }))
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
    priors <- make_multinom_priors(levels(df_train$y), ref=ref_lbl, sd_b=0.3, sd_int=3)
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
extract_mofa_factors_current <- function(mod, ids_train, factors = 1:2) {
  
  F_raw <- MOFA2::get_factors(
    mod,
    factors = factors,
    as.data.frame = FALSE
  )
  
  if (is.list(F_raw)) {
    F <- F_raw[[1]]
  } else {
    F <- F_raw
  }
  
  F <- as.matrix(F)
  
  ## Asegurar orientación: muestras x factores
  if (!all(ids_train %in% rownames(F)) && all(ids_train %in% colnames(F))) {
    F <- t(F)
  }
  
  missing_ids <- setdiff(ids_train, rownames(F))
  
  if (length(missing_ids) > 0) {
    stop(
      "Estos ids_train no están en los factores MOFA: ",
      paste(missing_ids, collapse = ", ")
    )
  }
  
  F <- F[ids_train, , drop = FALSE]
  F <- F[, seq_along(factors), drop = FALSE]
  colnames(F) <- paste0("Factor", factors)
  
  F
}


extract_mofa_weights_current <- function(mod, rfm, view_map, factors = 1:2) {
  
  W_raw <- MOFA2::get_weights(
    mod,
    views = "all",
    factors = factors,
    as.data.frame = FALSE
  )
  
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
      stop("No encuentro la vista ", v, " en los pesos del modelo MOFA.")
    }
    
    W <- as.matrix(W_raw[[v]])
    feats_rfm <- rownames(rfm[[code]]$train)
    
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
##===ejecucion=====
##===ejecucion adaptada al pipeline actual =====

if (!file.exists(RFM_FILE)) {
  stop("No existe RFM_FILE: ", RFM_FILE)
}

if (!file.exists(MOFA_MODEL_FILE)) {
  stop("No existe MOFA_MODEL_FILE: ", MOFA_MODEL_FILE)
}

tst_Data <- readRDS(RFM_FILE)
rfm <- tst_Data
mod <- readRDS(MOFA_MODEL_FILE)

view_map <- c(
  transcriptomica = "tx",
  proteomica      = "pr",
  metabolomica    = "me",
  clinical        = "cl"
)

ids_train <- rfm$grupo_train$id
ids_test  <- rfm$grupo_test$id

cat("\n--- Auditoría de orden de muestras ---\n")

audit_ids <- lapply(names(view_map), function(v) {
  code <- view_map[[v]]
  data.frame(
    vista = v,
    train_ok = identical(colnames(rfm[[code]]$train), ids_train),
    test_ok  = identical(colnames(rfm[[code]]$test),  ids_test),
    n_train  = ncol(rfm[[code]]$train),
    n_test   = ncol(rfm[[code]]$test)
  )
}) %>%
  dplyr::bind_rows()

print(audit_ids)

if (!all(audit_ids$train_ok)) {
  stop("El orden de muestras train no coincide entre grupo_train y alguna vista.")
}

if (!all(audit_ids$test_ok)) {
  stop("El orden de muestras test no coincide entre grupo_test y alguna vista.")
}

## ------------------------------------------------------------
## Factores y pesos MOFA como en el script anterior
## ------------------------------------------------------------

Z_seed <- extract_mofa_factors_current(
  mod       = mod,
  ids_train = ids_train,
  factors   = 1:2
)

y_seed <- factor(rfm$grupo_train$grupo[match(rownames(Z_seed), rfm$grupo_train$id)])
names(y_seed) <- rownames(Z_seed)

pesos_list <- extract_mofa_weights_current(
  mod      = mod,
  rfm      = rfm,
  view_map = view_map,
  factors  = 1:2
)

cat("\n--- Auditoría Z_seed / y_seed ---\n")
print(dim(Z_seed))
print(table(y_seed))
print(all(names(y_seed) == rownames(Z_seed)))

cat("\n--- Auditoría pesos_list ---\n")
print(lapply(pesos_list, dim))

## ------------------------------------------------------------
## Grid antiguo por cuantiles
## ------------------------------------------------------------

thr_grid <- expand.grid(
  q_tx = c(0.90, 0.95, 0.97),
  q_pr = c(0.70, 0.85),
  q_me = c(0.70, 0.85, 0.95, 0.97),
  q_cl = c(0.65, 0.70)
)

cat("\n--- Grid antiguo por cuantiles MOFA ---\n")
print(thr_grid)

write.csv(
  thr_grid,
  out_table("thr_grid_anterior_quantiles.csv"),
  row.names = FALSE
)

saveRDS(
  thr_grid,
  out_rds("thr_grid_anterior_quantiles.rds")
)

## ------------------------------------------------------------
## Construcción de datasets antiguos
## Incluye Factor1/Factor2 + features seleccionadas por loading
## ------------------------------------------------------------

return.list <- build_sets_from_grid(
  thr_grid   = thr_grid,
  tst_Data   = tst_Data,
  pesos_list = pesos_list,
  Z_seed     = Z_seed,
  y_seed     = y_seed,
  mod        = mod
)

saveRDS(
  return.list,
  out_rds("return_list_modelo_anterior.rds")
)

detach("package:MOFA2", unload = TRUE)

## ------------------------------------------------------------
## Evaluación antigua sobre todas las grids
## ------------------------------------------------------------

EVAL <- evaluate_all(
  return.list    = return.list,
  thr_grid       = thr_grid,
  ref            = "NP",
  do_brms        = TRUE,
  brms_algorithm = "mcmc",
  brms_iter      = 4000,
  brms_chains    = 4
)

saveRDS(
  EVAL,
  out_rds("EVAL_modelo_anterior.rds")
)

write.csv(
  EVAL$results,
  out_table("eval_modelo_anterior_results.csv"),
  row.names = FALSE
)

write.csv(
  EVAL$best_by_grid,
  out_table("eval_modelo_anterior_best_by_grid.csv"),
  row.names = FALSE
)

write.csv(
  EVAL$best_overall,
  out_table("eval_modelo_anterior_best_overall.csv"),
  row.names = FALSE
)

evaluacion <- EVAL$best_by_grid
evaluacion <- evaluacion[evaluacion$Accuracy < 1, ]

grid_id_auto <- evaluacion$GridID[which.max(evaluacion$BalAcc)]

cat("\n--- Grid automático antiguo por mayor BalAcc excluyendo Accuracy=1 ---\n")
print(grid_id_auto)

## ------------------------------------------------------------
## Grid final anterior
## Se conserva el comportamiento del script viejo: grid 33
## ------------------------------------------------------------

grid_id <- 33L

if (!grid_id %in% seq_len(nrow(thr_grid))) {
  stop("grid_id = 33L no existe en thr_grid.")
}

cat("\nGrid anterior usado para ajuste final:", grid_id, "\n")
cat("\nThresholds del grid anterior final:\n")
print(
  thr_grid %>%
    dplyr::mutate(GridID = dplyr::row_number()) %>%
    dplyr::filter(GridID == grid_id)
)

## ------------------------------------------------------------
## Datos consistentes con evaluate_all
## ------------------------------------------------------------

ref <- "NP"
backend <- EVAL$backend

X_train <- return.list$train[[grid_id]]
X_test  <- return.list$test[[grid_id]]

y_train <- stats::relevel(refactor_levels(return.list$y_train), ref = ref)
y_test  <- stats::relevel(refactor_levels(return.list$y_test),  ref = ref)

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

cat("\n--- Dimensiones modelo anterior final ---\n")
cat("TRAIN:", nrow(df_train), "x", ncol(df_train), "\n")
cat("TEST :", nrow(df_test),  "x", ncol(df_test),  "\n")
cat("p final anterior:", ncol(df_train) - 1L, "\n")

saveRDS(df_train, out_rds("df_train_final.rds"))
saveRDS(df_test,  out_rds("df_test_final.rds"))

## ------------------------------------------------------------
## BRMS final anterior
## ------------------------------------------------------------

priors <- make_multinom_priors(
  levels(df_train$y),
  ref    = ref,
  sd_b   = 0.3,
  sd_int = 3
)

ctrl <- if (identical(backend, "cmdstanr")) {
  list(adapt_delta = 0.95, max_treedepth = 12)
} else {
  list(adapt_delta = 0.95, max_treedepth = 12, init_r = 0.1)
}

fit_brm <- brms::brm(
  y ~ .,
  data      = df_train,
  family    = brms::categorical(),
  prior     = priors,
  chains    = 4,
  iter      = 8000,
  warmup    = 4000,
  refresh   = 0,
  backend   = backend,
  control   = ctrl,
  save_pars = brms::save_pars(all = TRUE)
)

saveRDS(
  fit_brm,
  out_rds("fit_brm_final_modelo_anterior.rds")
)

## ------------------------------------------------------------
## Probabilidades, métricas y matriz de confusión
## ------------------------------------------------------------

P <- brms_probs(fit_brm, df_test)[, levels(df_test$y), drop = FALSE]

metrics_grid_anterior <- metrics_from_probs(P, df_test$y)

pred_lab <- factor(
  colnames(P)[max.col(P, ties.method = "first")],
  levels = levels(df_test$y)
)

cm_grid_anterior <- caret::confusionMatrix(pred_lab, df_test$y)

summary(fit_brm)

diag <- list(
  rhat     = max(summary(fit_brm)$fixed[, "Rhat"], na.rm = TRUE),
  ess_bulk = min(summary(fit_brm)$fixed[, "Bulk_ESS"], na.rm = TRUE),
  ess_tail = min(summary(fit_brm)$fixed[, "Tail_ESS"], na.rm = TRUE)
)

print(diag)
print(metrics_grid_anterior)
print(cm_grid_anterior)

saveRDS(P, out_rds("probabilidades_test_modelo_anterior.rds"))
saveRDS(cm_grid_anterior, out_rds("confusion_matrix_modelo_anterior.rds"))

## ------------------------------------------------------------
## LOO / WAIC como en el anterior, pero protegido contra fallo
## ------------------------------------------------------------

loo_brm <- tryCatch(
  loo::loo(fit_brm, moment_match = TRUE),
  error = function(e) {
    warning("LOO con moment_match falló: ", conditionMessage(e))
    NULL
  }
)

waic_brm <- tryCatch(
  loo::waic(fit_brm),
  error = function(e) {
    warning("WAIC falló: ", conditionMessage(e))
    NULL
  }
)

print(loo_brm)
print(waic_brm)

saveRDS(loo_brm, out_rds("loo_brm_modelo_anterior.rds"))
saveRDS(waic_brm, out_rds("waic_brm_modelo_anterior.rds"))

## ------------------------------------------------------------
## ROC multiclase
## ------------------------------------------------------------

lev <- levels(df_test$y)

roc_list <- lapply(
  lev,
  function(lv) {
    pROC::roc(
      as.integer(df_test$y == lv),
      as.numeric(P[, lv]),
      quiet = TRUE
    )
  }
)

names(roc_list) <- lev

mroc <- pROC::multiclass.roc(df_test$y, P)

roc_df <- do.call(
  rbind,
  lapply(names(roc_list), function(lv) {
    r <- roc_list[[lv]]
    data.frame(
      specificity = rev(r$specificities),
      sensitivity = rev(r$sensitivities),
      class = lv
    )
  })
)

gg_roc <- ggplot(
  roc_df,
  aes(x = 1 - specificity, y = sensitivity, color = class)
) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_equal() +
  labs(
    title = sprintf("Modelo anterior: ROC one-vs-all (macro AUC = %.3f)", as.numeric(mroc$auc)),
    x = "1 - Especificidad",
    y = "Sensibilidad"
  ) +
  theme_minimal()

print(gg_roc)

ggsave(
  filename = out_plot("roc_one_vs_all_modelo_anterior.png"),
  plot = gg_roc,
  width = 7,
  height = 6,
  dpi = 300
)

## ------------------------------------------------------------
## Tablas finales
## ------------------------------------------------------------

df_metrics <- metrics_grid_anterior %>%
  dplyr::mutate(
    rhat      = diag$rhat,
    ess_bulk  = diag$ess_bulk,
    ess_tail  = diag$ess_tail,
    
    elpd_loo = if (!is.null(loo_brm)) loo_brm$estimates["elpd_loo", "Estimate"] else NA_real_,
    p_loo    = if (!is.null(loo_brm)) loo_brm$estimates["p_loo", "Estimate"] else NA_real_,
    looic    = if (!is.null(loo_brm)) loo_brm$estimates["looic", "Estimate"] else NA_real_,
    
    elpd_waic = if (!is.null(waic_brm)) waic_brm$estimates["elpd_waic", "Estimate"] else NA_real_,
    p_waic    = if (!is.null(waic_brm)) waic_brm$estimates["p_waic", "Estimate"] else NA_real_,
    waic      = if (!is.null(waic_brm)) waic_brm$estimates["waic", "Estimate"] else NA_real_
  )

stats_class <- as.data.frame(cm_grid_anterior$byClass) %>%
  tibble::rownames_to_column("class")

pred_individual_modelo_anterior <- data.frame(
  id       = rownames(df_test),
  real     = df_test$y,
  predicho = pred_lab,
  correcto = df_test$y == pred_lab,
  prob_max = apply(P, 1, max),
  P,
  check.names = FALSE
)

write.csv(
  df_metrics,
  out_table("estadisticas_modelo_global.csv"),
  row.names = FALSE
)

write.csv(
  stats_class,
  out_table("estadisticas_modelo_por_clase.csv"),
  row.names = FALSE
)

write.csv(
  pred_individual_modelo_anterior,
  out_table("predicciones_individuales_modelo_anterior.csv"),
  row.names = FALSE
)

## También guarda nombres antiguos por compatibilidad
write.csv(
  stats_class,
  out_table("estadisticas_modelo.csv"),
  row.names = FALSE
)

write.csv(
  df_metrics,
  out_table("estadisticas_modelo2.csv"),
  row.names = FALSE
)

run_config_anterior <- data.frame(
  analysis_name = ANALYSIS_NAME,
  analysis_dir = ANALYSIS_DIR,
  output_dir = OUTDIR,
  model_type = "modelo_anterior",
  selection_logic = "quantile_loading_mofa_plus_latent_factors",
  includes_latent_factors = TRUE,
  selected_grid_id_fixed = grid_id,
  selected_grid_id_auto_balacc = grid_id_auto,
  ref_class = ref,
  backend = backend,
  n_train = nrow(df_train),
  n_test = nrow(df_test),
  p_final = ncol(df_train) - 1L,
  brms_iter_eval = 4000,
  brms_chains_eval = 4,
  brms_iter_final = 8000,
  brms_chains_final = 4,
  stringsAsFactors = FALSE
)

write.csv(
  run_config_anterior,
  out_table("run_config_modelo_anterior.csv"),
  row.names = FALSE
)

saveRDS(
  run_config_anterior,
  out_rds("run_config_modelo_anterior.rds")
)

## ============================================================
## PANEL COMPATIBLE PARA DOWNSTREAM - MODELO ANTERIOR
## No es consenso multimodelo; es el panel final usado por el BRMS anterior
## Se guarda con los mismos nombres esperados por el downstream
## ============================================================

PANEL_OLD_DIR <- file.path(
  ANALYSIS_DIR,
  "panel_final_multimodelo_anterior",
  "tables"
)

dir.create(PANEL_OLD_DIR, recursive = TRUE, showWarnings = FALSE)

old_panel_tbl <- tibble::tibble(
  feature_model = setdiff(colnames(df_train), c("y", "Factor1", "Factor2")),
  Candidate = "Yes",
  model_source = "OLD",
  panel_origin = "df_train_final_modelado_anterior_grid33"
) %>%
  dplyr::distinct(feature_model, .keep_all = TRUE)

write.csv(
  old_panel_tbl,
  file.path(PANEL_OLD_DIR, "06_panel_feature_consensus_all.csv"),
  row.names = FALSE
)

write.csv(
  old_panel_tbl,
  file.path(PANEL_OLD_DIR, "07_PANEL_FINAL_STRICT_all_models.csv"),
  row.names = FALSE
)

write.csv(
  old_panel_tbl,
  file.path(PANEL_OLD_DIR, "08_PANEL_FINAL_RELAXED_mofa_3of4.csv"),
  row.names = FALSE
)

cat("\n--- Panel compatible modelo anterior guardado en ---\n")
cat(PANEL_OLD_DIR, "\n")
cat("n features panel anterior:", nrow(old_panel_tbl), "\n")


## ============================================================
## PLOTS DE AUDITORÍA MÉTRICAS - MODELO ANTERIOR
## Post-hoc: no reajusta modelos; usa tablas/RDS ya generados
## ============================================================

cat("\n============================================================\n")
cat("GENERANDO PLOTS DE AUDITORÍA - MODELO ANTERIOR\n")
cat("============================================================\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

METRIC_PLOTS_DIR <- out_plot("metricas_modelo_anterior")
dir.create(METRIC_PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

out_metric_plot <- function(filename) {
  file.path(METRIC_PLOTS_DIR, filename)
}

.read_csv_safe <- function(filename) {
  f <- out_table(filename)
  if (!file.exists(f)) {
    warning("No existe: ", f)
    return(NULL)
  }
  read.csv(f, check.names = FALSE)
}

.save_plot <- function(plot, filename, width = 8, height = 6) {
  ggsave(
    filename = out_metric_plot(filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

manifest <- tibble::tibble(
  plot = character(),
  description = character()
)

.add_manifest <- function(filename, description) {
  manifest <<- dplyr::bind_rows(
    manifest,
    tibble::tibble(
      plot = filename,
      description = description
    )
  )
}

## -----------------------------
## 1) Cargar tablas/RDS principales
## -----------------------------

eval_results <- .read_csv_safe("eval_modelo_anterior_results.csv")
eval_best    <- .read_csv_safe("eval_modelo_anterior_best_by_grid.csv")
pred_best    <- .read_csv_safe("predicciones_individuales_modelo_anterior.csv")
global_met   <- .read_csv_safe("estadisticas_modelo_global.csv")
class_met    <- .read_csv_safe("estadisticas_modelo_por_clase.csv")
run_cfg      <- .read_csv_safe("run_config_modelo_anterior.csv")

selected_grid_id <- NA_integer_

if (exists("grid_id")) {
  selected_grid_id <- as.integer(grid_id)
} else if (!is.null(run_cfg) && "selected_grid_id_fixed" %in% colnames(run_cfg)) {
  selected_grid_id <- as.integer(run_cfg$selected_grid_id_fixed[1])
}

cat("\nGrid final modelo anterior:", selected_grid_id, "\n")

## -----------------------------
## 2) Métricas por grid: NNET vs BRMS
## -----------------------------

if (!is.null(eval_results)) {
  
  metric_cols <- intersect(
    c("Accuracy", "BalAcc", "Kappa", "AUC", "LogLoss"),
    colnames(eval_results)
  )
  
  p_eval_all <- eval_results %>%
    dplyr::mutate(
      selected_final = GridID == selected_grid_id,
      GridID = factor(GridID)
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(metric_cols),
      names_to = "metric",
      values_to = "value"
    ) %>%
    ggplot(aes(x = p, y = value)) +
    geom_point(aes(shape = selected_final), size = 2.7, alpha = 0.85) +
    geom_line(aes(group = Model), alpha = 0.5) +
    facet_grid(metric ~ Model, scales = "free_y") +
    labs(
      title = "Modelo anterior: métricas por grid",
      subtitle = "Comparación NNET y BRMS_MCMC usando selección antigua por cuantiles MOFA.",
      x = "Número de predictores",
      y = "Valor",
      shape = "Grid final"
    ) +
    theme_minimal()
  
  print(p_eval_all)
  .save_plot(p_eval_all, "01_metricas_por_grid_nnet_brms.png", 11, 8)
  .add_manifest(
    "01_metricas_por_grid_nnet_brms.png",
    "Accuracy, Balanced Accuracy, Kappa, AUC y LogLoss por grid para NNET y BRMS_MCMC."
  )
}

## -----------------------------
## 3) BRMS: LogLoss por grid
## -----------------------------

if (!is.null(eval_results) && "Model" %in% colnames(eval_results)) {
  
  eval_brm <- eval_results %>%
    dplyr::filter(Model == "BRMS_MCMC") %>%
    dplyr::mutate(
      selected_final = GridID == selected_grid_id,
      GridID_factor = factor(GridID, levels = sort(unique(GridID)))
    )
  
  if (nrow(eval_brm) > 0) {
    
    p_brm_logloss <- eval_brm %>%
      ggplot(aes(x = GridID_factor, y = LogLoss)) +
      geom_col(alpha = 0.8) +
      geom_point(
        aes(shape = selected_final),
        size = 3
      ) +
      labs(
        title = "Modelo anterior: LogLoss BRMS por grid",
        subtitle = "El punto marca el grid final antiguo usado para el ajuste final.",
        x = "GridID",
        y = "LogLoss",
        shape = "Grid final"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    
    print(p_brm_logloss)
    .save_plot(p_brm_logloss, "02_brms_logloss_por_grid.png", 10, 6)
    .add_manifest(
      "02_brms_logloss_por_grid.png",
      "LogLoss del BRMS_MCMC para cada grid del modelo anterior."
    )
  }
}

## -----------------------------
## 4) BRMS: métricas predictivas por grid
## -----------------------------

if (!is.null(eval_results) && "Model" %in% colnames(eval_results)) {
  
  eval_brm <- eval_results %>%
    dplyr::filter(Model == "BRMS_MCMC") %>%
    dplyr::mutate(
      selected_final = GridID == selected_grid_id
    )
  
  brm_metric_cols <- intersect(
    c("Accuracy", "BalAcc", "AUC", "LogLoss"),
    colnames(eval_brm)
  )
  
  if (nrow(eval_brm) > 0 && length(brm_metric_cols) > 0) {
    
    p_brm_metrics <- eval_brm %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(brm_metric_cols),
        names_to = "metric",
        values_to = "value"
      ) %>%
      ggplot(aes(x = p, y = value)) +
      geom_point(aes(shape = selected_final), size = 3, alpha = 0.85) +
      geom_line(aes(group = metric), alpha = 0.5) +
      facet_wrap(~ metric, scales = "free_y") +
      labs(
        title = "Modelo anterior: BRMS métricas por tamaño de panel",
        x = "Número de predictores",
        y = "Valor",
        shape = "Grid final"
      ) +
      theme_minimal()
    
    print(p_brm_metrics)
    .save_plot(p_brm_metrics, "03_brms_metricas_por_p.png", 10, 7)
    .add_manifest(
      "03_brms_metricas_por_p.png",
      "Métricas BRMS_MCMC frente al número de predictores."
    )
  }
}

## -----------------------------
## 5) Matriz de confusión modelo final anterior
## -----------------------------

if (!is.null(pred_best) && all(c("real", "predicho") %in% colnames(pred_best))) {
  
  conf_best <- pred_best %>%
    dplyr::count(real, predicho, name = "n")
  
  p_conf <- conf_best %>%
    ggplot(aes(x = real, y = predicho, fill = n)) +
    geom_tile() +
    geom_text(aes(label = n), size = 5) +
    labs(
      title = "Modelo anterior: matriz de confusión BRMS final",
      x = "Clase real",
      y = "Clase predicha",
      fill = "n"
    ) +
    theme_minimal()
  
  print(p_conf)
  .save_plot(p_conf, "04_confusion_matrix_modelo_anterior.png", 6, 5)
  .add_manifest(
    "04_confusion_matrix_modelo_anterior.png",
    "Matriz de confusión del modelo anterior final."
  )
}

## -----------------------------
## 6) Probabilidades predictivas individuales
## -----------------------------

if (!is.null(pred_best)) {
  
  prob_cols <- setdiff(
    colnames(pred_best),
    c("id", "real", "predicho", "correcto", "prob_max")
  )
  
  prob_cols <- prob_cols[vapply(pred_best[prob_cols], is.numeric, logical(1))]
  
  if (length(prob_cols) > 0 && "id" %in% colnames(pred_best)) {
    
    prob_long <- pred_best %>%
      dplyr::select(id, real, predicho, correcto, dplyr::all_of(prob_cols)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(prob_cols),
        names_to = "class",
        values_to = "prob"
      ) %>%
      dplyr::mutate(
        id = factor(id, levels = pred_best$id)
      )
    
    p_prob_heat <- prob_long %>%
      ggplot(aes(x = class, y = id, fill = prob)) +
      geom_tile() +
      geom_text(aes(label = sprintf("%.2f", prob)), size = 3) +
      labs(
        title = "Modelo anterior: probabilidades predictivas por individuo",
        subtitle = "Filas = sujetos test; columnas = clases.",
        x = "Clase",
        y = "Individuo",
        fill = "Probabilidad"
      ) +
      theme_minimal()
    
    print(p_prob_heat)
    .save_plot(p_prob_heat, "05_probabilidades_individuales_modelo_anterior.png", 7, 6)
    .add_manifest(
      "05_probabilidades_individuales_modelo_anterior.png",
      "Heatmap de probabilidades predictivas individuales del modelo anterior."
    )
  }
}

## -----------------------------
## 7) Métricas globales modelo final
## -----------------------------

if (!is.null(global_met)) {
  
  global_metric_cols <- intersect(
    c(
      "Accuracy", "BalAcc", "Kappa", "AUC", "LogLoss",
      "rhat", "ess_bulk", "ess_tail",
      "elpd_loo", "p_loo", "looic",
      "elpd_waic", "p_waic", "waic"
    ),
    colnames(global_met)
  )
  
  if (length(global_metric_cols) > 0) {
    
    p_global <- global_met %>%
      dplyr::select(dplyr::all_of(global_metric_cols)) %>%
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = "metric",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        metric = factor(metric, levels = metric)
      ) %>%
      ggplot(aes(x = metric, y = value)) +
      geom_col() +
      facet_wrap(~ metric, scales = "free_y") +
      labs(
        title = "Modelo anterior: resumen global BRMS final",
        x = NULL,
        y = "Valor"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_blank())
    
    print(p_global)
    .save_plot(p_global, "06_metricas_globales_modelo_anterior.png", 12, 8)
    .add_manifest(
      "06_metricas_globales_modelo_anterior.png",
      "Resumen visual de métricas predictivas, convergencia, LOO y WAIC."
    )
  }
}

## -----------------------------
## 8) Métricas por clase
## -----------------------------

if (!is.null(class_met)) {
  
  if (!"class" %in% names(class_met)) {
    candidate_class_col <- intersect(
      c("", "X", "...1", "Class", "class"),
      names(class_met)
    )
    
    if (length(candidate_class_col) > 0) {
      class_met <- class_met %>%
        dplyr::rename(class = dplyr::all_of(candidate_class_col[1]))
    } else {
      class_met <- class_met %>%
        dplyr::mutate(class = paste0("class_", dplyr::row_number()))
    }
  }
  
  class_metric_cols <- intersect(
    c(
      "Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value",
      "Precision", "Recall", "F1", "Balanced Accuracy"
    ),
    names(class_met)
  )
  
  if (length(class_metric_cols) > 0) {
    
    p_class <- class_met %>%
      dplyr::select(class, dplyr::all_of(class_metric_cols)) %>%
      tidyr::pivot_longer(
        cols = -class,
        names_to = "metric",
        values_to = "value"
      ) %>%
      ggplot(aes(x = class, y = value)) +
      geom_col() +
      facet_wrap(~ metric, scales = "free_y") +
      labs(
        title = "Modelo anterior: métricas por clase",
        x = "Clase",
        y = "Valor"
      ) +
      theme_minimal()
    
    print(p_class)
    .save_plot(p_class, "07_metricas_por_clase_modelo_anterior.png", 10, 7)
    .add_manifest(
      "07_metricas_por_clase_modelo_anterior.png",
      "Sensitivity, specificity, precision, recall, F1 y Balanced Accuracy por clase."
    )
  }
}

## -----------------------------
## 9) Pareto-k desde LOO
## -----------------------------

loo_file <- out_rds("loo_brm_modelo_anterior.rds")

if (file.exists(loo_file)) {
  
  loo_brm_saved <- readRDS(loo_file)
  
  if (!is.null(loo_brm_saved)) {
    
    pareto_vals <- tryCatch(
      loo::pareto_k_values(loo_brm_saved),
      error = function(e) NULL
    )
    
    if (!is.null(pareto_vals)) {
      
      pareto_k_tbl <- tibble::tibble(
        obs = seq_along(pareto_vals),
        pareto_k = as.numeric(pareto_vals)
      )
      
      p_pareto <- pareto_k_tbl %>%
        ggplot(aes(x = obs, y = pareto_k)) +
        geom_point(size = 3) +
        geom_hline(yintercept = 0.7, linetype = 2) +
        geom_hline(yintercept = 1.0, linetype = 2) +
        labs(
          title = "Modelo anterior: diagnóstico Pareto-k de LOO",
          subtitle = "Líneas de referencia: 0.7 y 1.0.",
          x = "Observación",
          y = "Pareto-k"
        ) +
        theme_minimal()
      
      print(p_pareto)
      .save_plot(p_pareto, "08_pareto_k_loo_modelo_anterior.png", 7, 5)
      .add_manifest(
        "08_pareto_k_loo_modelo_anterior.png",
        "Valores Pareto-k del LOO para detectar observaciones influyentes."
      )
      
      write.csv(
        pareto_k_tbl,
        out_table("pareto_k_loo_modelo_anterior.csv"),
        row.names = FALSE
      )
    }
  }
}

## -----------------------------
## 10) Coeficientes BRMS con intervalo creíble
## -----------------------------

fit_file <- out_rds("fit_brm_final_modelo_anterior.rds")

if (file.exists(fit_file)) {
  
  fit_saved <- readRDS(fit_file)
  
  coef_tbl <- tryCatch({
    fx <- as.data.frame(brms::fixef(fit_saved, robust = TRUE))
    fx$term <- rownames(fx)
    
    fx %>%
      tibble::as_tibble() %>%
      dplyr::filter(!grepl("Intercept", term)) %>%
      dplyr::mutate(
        class = sub("^mu([^_]+)_.*$", "\\1", term),
        feature_model = sub("^mu[^_]+_", "", term),
        ci_excludes_zero = (`Q2.5` > 0) | (`Q97.5` < 0),
        abs_estimate = abs(Estimate)
      ) %>%
      dplyr::select(
        class,
        feature_model,
        Estimate,
        Est.Error,
        Q2.5,
        Q97.5,
        abs_estimate,
        ci_excludes_zero,
        term
      )
  }, error = function(e) {
    warning("No se pudieron extraer coeficientes BRMS: ", conditionMessage(e))
    NULL
  })
  
  if (!is.null(coef_tbl) && nrow(coef_tbl) > 0) {
    
    write.csv(
      coef_tbl,
      out_table("coeficientes_brms_modelo_anterior.csv"),
      row.names = FALSE
    )
    
    top_features_coef <- coef_tbl %>%
      dplyr::group_by(feature_model) %>%
      dplyr::summarise(
        max_abs = max(abs_estimate, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(max_abs)) %>%
      dplyr::slice_head(n = 25) %>%
      dplyr::pull(feature_model)
    
    coef_plot_tbl <- coef_tbl %>%
      dplyr::filter(feature_model %in% top_features_coef) %>%
      dplyr::mutate(
        feature_model = factor(feature_model, levels = rev(top_features_coef))
      )
    
    p_coef <- coef_plot_tbl %>%
      ggplot(aes(x = Estimate, y = feature_model)) +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_errorbarh(
        aes(xmin = `Q2.5`, xmax = `Q97.5`),
        height = 0.2
      ) +
      geom_point(aes(shape = ci_excludes_zero), size = 2.5) +
      facet_wrap(~ class, scales = "free_x") +
      labs(
        title = "Modelo anterior: top coeficientes BRMS con intervalo creíble",
        subtitle = "Top 25 features por máximo |Estimate| entre clases.",
        x = "Estimate posterior",
        y = "Feature",
        shape = "IC 95% excluye 0"
      ) +
      theme_minimal()
    
    print(p_coef)
    .save_plot(p_coef, "09_coeficientes_brms_top25_modelo_anterior.png", 12, 9)
    .add_manifest(
      "09_coeficientes_brms_top25_modelo_anterior.png",
      "Coeficientes BRMS finales con intervalos creíbles por clase."
    )
  }
}

## -----------------------------
## 11) Guardar índice de plots
## -----------------------------

write.csv(
  manifest,
  out_table("manifest_metricas_modelo_anterior_plots.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("PLOTS DE AUDITORÍA - MODELO ANTERIOR TERMINADOS\n")
cat("Plots guardados en:", METRIC_PLOTS_DIR, "\n")
cat("Manifest guardado en:", out_table("manifest_metricas_modelo_anterior_plots.csv"), "\n")
cat("============================================================\n")
cat("\n============================================================\n")
cat("MODELO ANTERIOR TERMINADO\n")
cat("Guardado en:", OUTDIR, "\n")
cat("Tablas:", TABLES_DIR, "\n")
cat("Plots:", PLOTS_DIR, "\n")
cat("RDS:", RDS_DIR, "\n")
cat("============================================================\n")
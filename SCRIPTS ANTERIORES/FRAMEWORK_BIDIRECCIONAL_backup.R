## ============================================================
## FRAMEWORK EXPLICATIVO BIDIRECCIONAL MULTI-ÓMICO
## ------------------------------------------------------------
## Entrada:
##   X_S : matriz/data.frame muestras x features seleccionadas
##   y   : grupo clínico NP / FD / SO
##   v   : vista de cada feature: tx / pr / me / cl
##
## Modelo:
##   Y ~ Categorical(pi)
##   X_j | Y = k ~ Student-t(df_t, mu_jk, sigma_jk)
##
## Salidas:
##   1) feature_signature
##   2) view_summary
##   3) permutation_distribution
##   4) plots opcionales
##
## No usa MOFA, BRMS, mRMR, train/test ni objetos intermedios.
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

## ============================================================
## 1) HELPERS
## ============================================================
.save_plot_publication <- function(filename_base,
                                   plot,
                                   width = 12,
                                   height = 8,
                                   dpi = 600) {
  
  png_file <- paste0(filename_base, ".png")
  pdf_file <- paste0(filename_base, ".pdf")
  
  ## PNG de alta resolución con antialiasing si ragg está disponible
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(
      filename = png_file,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      device = ragg::agg_png,
      bg = "white"
    )
  } else {
    ggplot2::ggsave(
      filename = png_file,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }
  
  ## PDF vectorial: este es el importante para tesis/paper
  ggplot2::ggsave(
    filename = pdf_file,
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )
  
  invisible(
    list(
      png = png_file,
      pdf = pdf_file
    )
  )
}
.safe_scale <- function(x, method = c("zscore", "robust", "none")) {
  method <- match.arg(method)
  x <- as.numeric(x)
  
  if (method == "none") return(x)
  
  if (method == "zscore") {
    mu <- mean(x, na.rm = TRUE)
    ss <- sd(x, na.rm = TRUE)
    if (!is.finite(ss) || ss == 0) ss <- 1
    return((x - mu) / ss)
  }
  
  if (method == "robust") {
    mu <- median(x, na.rm = TRUE)
    ss <- stats::mad(x, na.rm = TRUE)
    if (!is.finite(ss) || ss == 0) ss <- sd(x, na.rm = TRUE)
    if (!is.finite(ss) || ss == 0) ss <- 1
    return((x - mu) / ss)
  }
}

.prepare_views <- function(v, features) {
  
  if (is.data.frame(v)) {
    
    if (!"feature_model" %in% colnames(v)) {
      stop("Si v es data.frame debe tener columna feature_model.")
    }
    
    view_col <- intersect(c("view_code", "view"), colnames(v))[1]
    
    if (is.na(view_col)) {
      stop("Si v es data.frame debe tener columna view_code o view.")
    }
    
    out <- v[[view_col]]
    names(out) <- v$feature_model
    out <- out[features]
    
  } else {
    
    out <- as.character(v)
    
    if (is.null(names(out))) {
      if (length(out) != length(features)) {
        stop("Si v no tiene names, length(v) debe ser igual a ncol(X_S).")
      }
      names(out) <- features
    }
    
    out <- out[features]
  }
  
  if (any(is.na(out) | out == "")) {
    bad <- names(out)[is.na(out) | out == ""]
    stop(
      "Hay features sin vista en v: ",
      paste(bad, collapse = ", ")
    )
  }
  
  out
}

.dens_t <- function(x, mu, sigma, df_t = 4) {
  sigma <- pmax(sigma, 1e-8)
  stats::dt((x - mu) / sigma, df = df_t) / sigma
}

.fit_group_t_params <- function(x,
                                y,
                                y_levels,
                                scale_floor = 1e-4) {
  
  out <- lapply(y_levels, function(k) {
    
    xk <- x[y == k]
    xk <- xk[is.finite(xk)]
    
    n_k <- length(xk)
    
    if (n_k == 0) {
      return(tibble(
        group = k,
        n = 0L,
        mu = NA_real_,
        sigma = NA_real_
      ))
    }
    
    mu_k <- median(xk, na.rm = TRUE)
    
    sigma_k <- stats::mad(xk, na.rm = TRUE)
    
    if (!is.finite(sigma_k) || sigma_k <= 0) {
      sigma_k <- sd(xk, na.rm = TRUE)
    }
    
    if (!is.finite(sigma_k) || sigma_k <= 0) {
      sigma_k <- stats::mad(x, na.rm = TRUE)
    }
    
    if (!is.finite(sigma_k) || sigma_k <= 0) {
      sigma_k <- sd(x, na.rm = TRUE)
    }
    
    if (!is.finite(sigma_k) || sigma_k <= 0) {
      sigma_k <- 1
    }
    
    tibble(
      group = k,
      n = n_k,
      mu = mu_k,
      sigma = max(sigma_k, scale_floor)
    )
  })
  
  bind_rows(out)
}

.make_priors <- function(y,
                         y_levels,
                         alpha = 0.5) {
  
  tab <- table(factor(y, levels = y_levels))
  
  pi_k <- (as.numeric(tab) + alpha) /
    (sum(tab) + length(y_levels) * alpha)
  
  names(pi_k) <- y_levels
  
  pi_k
}

.posterior_group_at_x <- function(x0,
                                  params,
                                  pi_k,
                                  y_levels,
                                  df_t = 4) {
  
  dens <- sapply(y_levels, function(k) {
    row <- params[params$group == k, , drop = FALSE]
    .dens_t(x0, row$mu, row$sigma, df_t = df_t)
  })
  
  num <- as.numeric(pi_k[y_levels]) * dens
  den <- sum(num, na.rm = TRUE)
  
  if (!is.finite(den) || den <= 0) {
    out <- rep(NA_real_, length(y_levels))
    names(out) <- y_levels
    return(out)
  }
  
  out <- num / den
  names(out) <- y_levels
  out
}

.mi_js_from_params <- function(params,
                               pi_k,
                               y_levels,
                               df_t = 4,
                               n_grid = 2048,
                               eps = 1e-12) {
  
  mus <- params$mu
  sig <- params$sigma
  
  xmin <- min(mus - 8 * sig, na.rm = TRUE)
  xmax <- max(mus + 8 * sig, na.rm = TRUE)
  
  if (!is.finite(xmin) || !is.finite(xmax) || xmin == xmax) {
    return(list(
      MI_nats = NA_real_,
      MI_bits = NA_real_,
      MI_norm = NA_real_,
      JS_pairs = tibble()
    ))
  }
  
  grid <- seq(xmin, xmax, length.out = n_grid)
  dx <- mean(diff(grid))
  
  dens_mat <- sapply(y_levels, function(k) {
    row <- params[params$group == k, , drop = FALSE]
    .dens_t(grid, row$mu, row$sigma, df_t = df_t)
  })
  
  dens_mat <- as.matrix(dens_mat)
  colnames(dens_mat) <- y_levels
  
  pi_vec <- as.numeric(pi_k[y_levels])
  names(pi_vec) <- y_levels
  
  mix <- as.numeric(dens_mat %*% pi_vec)
  
  ## Mutual information:
  ## I(X;Y) = sum_k pi_k ∫ f_k(x) log[f_k(x)/f_mix(x)] dx
  mi_k <- sapply(y_levels, function(k) {
    fk <- dens_mat[, k]
    pk <- pi_vec[k]
    sum(pk * fk * log((fk + eps) / (mix + eps))) * dx
  })
  
  MI_nats <- sum(mi_k, na.rm = TRUE)
  MI_bits <- MI_nats / log(2)
  
  H_y <- -sum(pi_vec * log(pi_vec + eps))
  MI_norm <- ifelse(H_y > 0, MI_nats / H_y, NA_real_)
  
  ## Jensen-Shannon pairwise entre grupos
  pairs <- combn(y_levels, 2, simplify = FALSE)
  
  JS_pairs <- bind_rows(lapply(pairs, function(pp) {
    
    a <- pp[1]
    b <- pp[2]
    
    fa <- dens_mat[, a]
    fb <- dens_mat[, b]
    m <- 0.5 * fa + 0.5 * fb
    
    js <- 0.5 * sum(fa * log((fa + eps) / (m + eps))) * dx +
      0.5 * sum(fb * log((fb + eps) / (m + eps))) * dx
    
    tibble(
      contrast = paste0(a, "_vs_", b),
      JS_nats = js,
      JS_bits = js / log(2)
    )
  }))
  
  list(
    MI_nats = MI_nats,
    MI_bits = MI_bits,
    MI_norm = MI_norm,
    JS_pairs = JS_pairs
  )
}

.bootstrap_direction <- function(x,
                                 y,
                                 y_levels,
                                 B_boot = 499,
                                 seed = 123) {
  
  set.seed(seed)
  
  idx_by_group <- split(seq_along(y), y)
  idx_by_group <- idx_by_group[y_levels]
  
  boot_mu <- matrix(
    NA_real_,
    nrow = B_boot,
    ncol = length(y_levels),
    dimnames = list(NULL, y_levels)
  )
  
  for (b in seq_len(B_boot)) {
    for (k in y_levels) {
      idx <- idx_by_group[[k]]
      if (length(idx) == 0) next
      
      idx_b <- sample(idx, size = length(idx), replace = TRUE)
      boot_mu[b, k] <- median(x[idx_b], na.rm = TRUE)
    }
  }
  
  top_boot <- apply(boot_mu, 1, function(z) {
    if (all(!is.finite(z))) return(NA_character_)
    names(z)[which.max(z)]
  })
  
  top_prob_tbl <- tibble(
    group = y_levels,
    prob_top_location = sapply(y_levels, function(k) {
      mean(top_boot == k, na.rm = TRUE)
    })
  )
  
  get_delta <- function(a, b) {
    d <- boot_mu[, a] - boot_mu[, b]
    
    tibble(
      contrast = paste0(a, "_minus_", b),
      estimate = median(d, na.rm = TRUE),
      ci_low = as.numeric(quantile(d, 0.025, na.rm = TRUE)),
      ci_high = as.numeric(quantile(d, 0.975, na.rm = TRUE)),
      prob_gt0 = mean(d > 0, na.rm = TRUE),
      prob_lt0 = mean(d < 0, na.rm = TRUE)
    )
  }
  
  delta_tbl <- bind_rows(
    get_delta("FD", "NP"),
    get_delta("SO", "NP"),
    get_delta("SO", "FD")
  )
  
  list(
    top_prob_tbl = top_prob_tbl,
    delta_tbl = delta_tbl
  )
}

.classify_pattern <- function(params,
                              post_high,
                              post_low,
                              pi_k,
                              boot_dir,
                              y_levels,
                              prob_thr_strong = 0.90,
                              prob_thr_moderate = 0.80) {
  
  loc <- params$mu
  names(loc) <- params$group
  
  top_group <- names(loc)[which.max(loc)]
  low_group <- names(loc)[which.min(loc)]
  
  high_group <- names(post_high)[which.max(post_high)]
  low_post_group <- names(post_low)[which.max(post_low)]
  
  top_prob <- boot_dir$top_prob_tbl %>%
    filter(group == top_group) %>%
    pull(prob_top_location)
  
  if (length(top_prob) == 0) top_prob <- NA_real_
  
  concordant_high <- isTRUE(top_group == high_group)
  
  delta_tbl <- boot_dir$delta_tbl
  
  p_FD_gt_NP <- delta_tbl %>%
    filter(contrast == "FD_minus_NP") %>%
    pull(prob_gt0)
  
  p_SO_gt_NP <- delta_tbl %>%
    filter(contrast == "SO_minus_NP") %>%
    pull(prob_gt0)
  
  p_SO_gt_FD <- delta_tbl %>%
    filter(contrast == "SO_minus_FD") %>%
    pull(prob_gt0)
  
  if (length(p_FD_gt_NP) == 0) p_FD_gt_NP <- NA_real_
  if (length(p_SO_gt_NP) == 0) p_SO_gt_NP <- NA_real_
  if (length(p_SO_gt_FD) == 0) p_SO_gt_FD <- NA_real_
  
  posterior_enrichment_high <- log(
    pmax(post_high[high_group], 1e-12) /
      pmax(pi_k[high_group], 1e-12)
  )
  
  pattern <- case_when(
    is.finite(p_FD_gt_NP) &&
      is.finite(p_SO_gt_NP) &&
      p_FD_gt_NP >= prob_thr_moderate &&
      p_SO_gt_NP >= prob_thr_moderate &&
      p_SO_gt_FD > 0.20 &&
      p_SO_gt_FD < 0.80 ~ "shared_FD_SO_high",
    
    top_group == "FD" &&
      concordant_high &&
      is.finite(top_prob) &&
      top_prob >= prob_thr_moderate ~ "FD_like_high",
    
    top_group == "SO" &&
      concordant_high &&
      is.finite(top_prob) &&
      top_prob >= prob_thr_moderate ~ "SO_like_high",
    
    top_group == "NP" &&
      concordant_high &&
      is.finite(top_prob) &&
      top_prob >= prob_thr_moderate ~ "NP_like_high",
    
    is.finite(p_SO_gt_FD) &&
      p_SO_gt_FD >= prob_thr_strong ~ "SO_gt_FD_splitter",
    
    is.finite(p_SO_gt_FD) &&
      p_SO_gt_FD <= (1 - prob_thr_strong) ~ "FD_gt_SO_splitter",
    
    TRUE ~ "uncertain"
  )
  
  tibble(
    top_location_group = top_group,
    low_location_group = low_group,
    posterior_high_group = high_group,
    posterior_low_group = low_post_group,
    top_location_prob_boot = as.numeric(top_prob),
    concordant_high = concordant_high,
    posterior_enrichment_high = as.numeric(posterior_enrichment_high),
    bidirectional_pattern = pattern
  )
}

.evidence_grade <- function(p_perm,
                            q_perm,
                            MI_norm,
                            top_location_prob_boot,
                            concordant_high,
                            posterior_enrichment_high,
                            pattern,
                            q_A = 0.10,
                            q_B = 0.20,
                            p_C = 0.05,
                            q_D = 0.20,
                            p_D = 0.05,
                            MI_D = 0.20) {
  
  concordant_high <- dplyr::coalesce(as.logical(concordant_high), FALSE)
  pattern <- dplyr::coalesce(as.character(pattern), "uncertain")
  
  has_dependency <- (
    (is.finite(q_perm) & q_perm < q_D) |
      (is.finite(p_perm) & p_perm < p_D) |
      (is.finite(MI_norm) & MI_norm >= MI_D)
  )
  
  dplyr::case_when(
    ## ------------------------------------------------------------
    ## A: evidencia fuerte
    ## Dependencia FDR significativa, dirección estable,
    ## concordancia bidireccional y patrón interpretable.
    ## ------------------------------------------------------------
    is.finite(q_perm) &
      q_perm < q_A &
      is.finite(MI_norm) &
      MI_norm >= 0.05 &
      is.finite(top_location_prob_boot) &
      top_location_prob_boot >= 0.90 &
      concordant_high &
      is.finite(posterior_enrichment_high) &
      posterior_enrichment_high > 0 &
      pattern != "uncertain" ~ "A_strong",
    
    ## ------------------------------------------------------------
    ## B: evidencia moderada
    ## Dependencia aceptable y patrón interpretable,
    ## pero menor estabilidad que A.
    ## ------------------------------------------------------------
    is.finite(q_perm) &
      q_perm < q_B &
      is.finite(top_location_prob_boot) &
      top_location_prob_boot >= 0.80 &
      pattern != "uncertain" ~ "B_moderate",
    
    ## ------------------------------------------------------------
    ## C: exploratoria
    ## Hay señal nominal y patrón direccional,
    ## pero no alcanza A/B.
    ## ------------------------------------------------------------
    is.finite(p_perm) &
      p_perm < p_C &
      pattern != "uncertain" ~ "C_exploratory",
    
    ## ------------------------------------------------------------
    ## D: dependiente pero ambigua
    ## Hay dependencia feature-grupo, pero el patrón no es claro.
    ## Esto corrige el problema de D_uncertain.
    ## ------------------------------------------------------------
    has_dependency &
      pattern == "uncertain" ~ "D_ambiguous_dependent",
    
    ## ------------------------------------------------------------
    ## E: sin evidencia suficiente
    ## No hay dependencia clara ni patrón interpretable.
    ## ------------------------------------------------------------
    TRUE ~ "E_no_evidence"
  )
}



.make_final_bidirectional_table <- function(feature_signature,
                                            y_levels = c("NP", "FD", "SO"),
                                            remove_prefix = TRUE) {
  
  fs <- tibble::as_tibble(feature_signature)
  
  get_col <- function(df, nm, default = NA_real_) {
    if (nm %in% colnames(df)) {
      df[[nm]]
    } else {
      rep(default, nrow(df))
    }
  }
  
  ## ------------------------------------------------------------
  ## Contraste principal: el de mayor magnitud absoluta.
  ## ------------------------------------------------------------
  
  est_mat <- cbind(
    FD_minus_NP = get_col(fs, "estimate__FD_minus_NP"),
    SO_minus_NP = get_col(fs, "estimate__SO_minus_NP"),
    SO_minus_FD = get_col(fs, "estimate__SO_minus_FD")
  )
  
  idx_main <- apply(abs(est_mat), 1, function(z) {
    if (all(!is.finite(z))) return(NA_integer_)
    which.max(z)
  })
  
  main_contrast <- ifelse(
    is.na(idx_main),
    NA_character_,
    colnames(est_mat)[idx_main]
  )
  
  main_estimate <- rep(NA_real_, nrow(fs))
  
  ok_idx <- which(!is.na(idx_main))
  
  if (length(ok_idx) > 0) {
    main_estimate[ok_idx] <- est_mat[cbind(ok_idx, idx_main[ok_idx])]
  }
  
  main_ci_low <- dplyr::case_when(
    main_contrast == "FD_minus_NP" ~ get_col(fs, "ci_low__FD_minus_NP"),
    main_contrast == "SO_minus_NP" ~ get_col(fs, "ci_low__SO_minus_NP"),
    main_contrast == "SO_minus_FD" ~ get_col(fs, "ci_low__SO_minus_FD"),
    TRUE ~ NA_real_
  )
  
  main_ci_high <- dplyr::case_when(
    main_contrast == "FD_minus_NP" ~ get_col(fs, "ci_high__FD_minus_NP"),
    main_contrast == "SO_minus_NP" ~ get_col(fs, "ci_high__SO_minus_NP"),
    main_contrast == "SO_minus_FD" ~ get_col(fs, "ci_high__SO_minus_FD"),
    TRUE ~ NA_real_
  )
  
  main_prob_gt0 <- dplyr::case_when(
    main_contrast == "FD_minus_NP" ~ get_col(fs, "prob_gt0__FD_minus_NP"),
    main_contrast == "SO_minus_NP" ~ get_col(fs, "prob_gt0__SO_minus_NP"),
    main_contrast == "SO_minus_FD" ~ get_col(fs, "prob_gt0__SO_minus_FD"),
    TRUE ~ NA_real_
  )
  
  main_prob_lt0 <- dplyr::case_when(
    main_contrast == "FD_minus_NP" ~ get_col(fs, "prob_lt0__FD_minus_NP"),
    main_contrast == "SO_minus_NP" ~ get_col(fs, "prob_lt0__SO_minus_NP"),
    main_contrast == "SO_minus_FD" ~ get_col(fs, "prob_lt0__SO_minus_FD"),
    TRUE ~ NA_real_
  )
  
  main_prob_direction <- ifelse(
    is.finite(main_estimate) & main_estimate >= 0,
    main_prob_gt0,
    main_prob_lt0
  )
  
  fs %>%
    mutate(
      feature_label = if (remove_prefix) {
        sub("^[^_]+_", "", feature_model)
      } else {
        feature_model
      },
      
      biological_axis = dplyr::case_when(
        bidirectional_pattern == "SO_like_high" ~ "SO_like_axis",
        bidirectional_pattern == "FD_like_high" ~ "FD_like_axis",
        bidirectional_pattern == "NP_like_high" ~ "NP_like_axis",
        bidirectional_pattern %in% c("SO_gt_FD_splitter", "FD_gt_SO_splitter") ~ "FD_SO_splitter_axis",
        bidirectional_pattern == "shared_FD_SO_high" ~ "shared_FD_SO_axis",
        evidence_grade == "D_ambiguous_dependent" ~ "ambiguous_dependent_axis",
        evidence_grade == "E_no_evidence" ~ "no_evidence_axis",
        TRUE ~ "unclassified_axis"
      ),
      
      dependency_status = dplyr::case_when(
        evidence_grade %in% c("A_strong", "B_moderate", "C_exploratory") ~
          "dependent_interpretable",
        evidence_grade == "D_ambiguous_dependent" ~
          "dependent_but_directionally_ambiguous",
        evidence_grade == "E_no_evidence" ~
          "no_sufficient_evidence",
        TRUE ~ "unclassified"
      ),
      
      pattern_interpretation = dplyr::case_when(
        bidirectional_pattern == "SO_like_high" ~
          "Feature alta característica de SO; feature alta favorece posterior SO.",
        
        bidirectional_pattern == "FD_like_high" ~
          "Feature alta característica de FD; feature alta favorece posterior FD.",
        
        bidirectional_pattern == "NP_like_high" ~
          "Feature alta característica de NP; feature alta favorece posterior NP.",
        
        bidirectional_pattern == "SO_gt_FD_splitter" ~
          "Feature separadora FD/SO, con valores mayores en SO que en FD.",
        
        bidirectional_pattern == "FD_gt_SO_splitter" ~
          "Feature separadora FD/SO, con valores mayores en FD que en SO.",
        
        bidirectional_pattern == "shared_FD_SO_high" ~
          "Feature alta en FD y SO frente a NP; patrón compartido no-NP.",
        
        evidence_grade == "D_ambiguous_dependent" ~
          "Existe dependencia feature-grupo, pero la dirección bidireccional no es limpia.",
        
        evidence_grade == "E_no_evidence" ~
          "No hay evidencia suficiente de dependencia feature-grupo ni patrón direccional claro.",
        
        TRUE ~ "Patrón no clasificado."
      ),
      
      main_contrast = main_contrast,
      main_contrast_estimate = main_estimate,
      main_contrast_ci_low = main_ci_low,
      main_contrast_ci_high = main_ci_high,
      main_contrast_prob_direction = main_prob_direction,
      
      main_reading = paste0(
        "Máx feature|grupo = ", top_location_group,
        "; máx grupo|feature alta = ", posterior_high_group,
        "; patrón = ", bidirectional_pattern,
        "; evidencia = ", evidence_grade
      )
    ) %>%
    select(
      feature_model,
      feature_label,
      view_code,
      evidence_grade,
      bidirectional_pattern,
      biological_axis,
      dependency_status,
      pattern_interpretation,
      
      n_used,
      MI_norm,
      MI_bits,
      max_JS_bits,
      p_perm,
      q_perm,
      
      top_location_group,
      posterior_high_group,
      concordant_high,
      top_location_prob_boot,
      posterior_enrichment_high,
      
      any_of(c("mu_NP", "mu_FD", "mu_SO")),
      any_of(c("post_high_NP", "post_high_FD", "post_high_SO")),
      
      main_contrast,
      main_contrast_estimate,
      main_contrast_ci_low,
      main_contrast_ci_high,
      main_contrast_prob_direction,
      main_reading
    ) %>%
    arrange(
      factor(
        evidence_grade,
        levels = c(
          "A_strong",
          "B_moderate",
          "C_exploratory",
          "D_ambiguous_dependent",
          "E_no_evidence"
        )
      ),
      view_code,
      desc(MI_norm),
      desc(max_JS_bits)
    )
}

.compute_one_feature <- function(feature,
                                 x,
                                 y,
                                 v_j,
                                 y_levels = c("NP", "FD", "SO"),
                                 df_t = 4,
                                 n_grid = 2048,
                                 B_perm = 999,
                                 B_boot = 499,
                                 scale_method = "zscore",
                                 seed = 123,
                                 perm_stat = c("MI", "maxJS")) {
  
  perm_stat <- match.arg(perm_stat)
  
  keep <- is.finite(x) & !is.na(y)
  x <- as.numeric(x[keep])
  y <- factor(y[keep], levels = y_levels)
  
  x <- .safe_scale(x, method = scale_method)
  
  if (length(unique(y)) < 2 || length(unique(x)) < 2) {
    return(list(
      feature_row = tibble(
        feature_model = feature,
        view_code = v_j,
        n_used = length(x),
        error = "feature_sin_variacion_o_grupos_insuficientes"
      ),
      perm_tbl = tibble()
    ))
  }
  
  pi_k <- .make_priors(y, y_levels)
  
  params <- .fit_group_t_params(
    x = x,
    y = y,
    y_levels = y_levels
  )
  
  stat_obj <- .mi_js_from_params(
    params = params,
    pi_k = pi_k,
    y_levels = y_levels,
    df_t = df_t,
    n_grid = n_grid
  )
  
  JS_wide <- stat_obj$JS_pairs %>%
    select(contrast, JS_bits) %>%
    pivot_wider(
      names_from = contrast,
      values_from = JS_bits,
      names_prefix = "JS_bits_"
    )
  
  max_JS_bits <- max(stat_obj$JS_pairs$JS_bits, na.rm = TRUE)
  
  T_obs <- if (perm_stat == "MI") {
    stat_obj$MI_bits
  } else {
    max_JS_bits
  }
  
  ## ------------------------------------------------------------
  ## Permutación de etiquetas:
  ## H0: independencia entre feature y grupo.
  ## ------------------------------------------------------------
  
  set.seed(seed)
  
  T_perm <- numeric(B_perm)
  
  for (b in seq_len(B_perm)) {
    
    y_perm <- sample(y)
    pi_perm <- .make_priors(y_perm, y_levels)
    
    params_perm <- .fit_group_t_params(
      x = x,
      y = y_perm,
      y_levels = y_levels
    )
    
    stat_perm <- .mi_js_from_params(
      params = params_perm,
      pi_k = pi_perm,
      y_levels = y_levels,
      df_t = df_t,
      n_grid = n_grid
    )
    
    T_perm[b] <- if (perm_stat == "MI") {
      stat_perm$MI_bits
    } else {
      max(stat_perm$JS_pairs$JS_bits, na.rm = TRUE)
    }
  }
  
  p_perm <- (1 + sum(T_perm >= T_obs, na.rm = TRUE)) / (B_perm + 1)
  
  perm_tbl <- tibble(
    feature_model = feature,
    view_code = v_j,
    perm_id = seq_len(B_perm),
    T_perm = T_perm,
    T_obs = T_obs,
    perm_stat = perm_stat
  )
  
  ## ------------------------------------------------------------
  ## Bootstrap direccional.
  ## ------------------------------------------------------------
  
  boot_dir <- .bootstrap_direction(
    x = x,
    y = y,
    y_levels = y_levels,
    B_boot = B_boot,
    seed = seed + 999
  )
  
  delta_wide <- boot_dir$delta_tbl %>%
    pivot_wider(
      names_from = contrast,
      values_from = c(estimate, ci_low, ci_high, prob_gt0, prob_lt0),
      names_sep = "__"
    )
  
  top_prob_wide <- boot_dir$top_prob_tbl %>%
    pivot_wider(
      names_from = group,
      values_from = prob_top_location,
      names_prefix = "prob_top_location_"
    )
  
  ## ------------------------------------------------------------
  ## Posteriores P(grupo | feature baja/media/alta).
  ## ------------------------------------------------------------
  
  q10 <- as.numeric(quantile(x, 0.10, na.rm = TRUE))
  q50 <- as.numeric(quantile(x, 0.50, na.rm = TRUE))
  q90 <- as.numeric(quantile(x, 0.90, na.rm = TRUE))
  
  post_low <- .posterior_group_at_x(
    x0 = q10,
    params = params,
    pi_k = pi_k,
    y_levels = y_levels,
    df_t = df_t
  )
  
  post_mid <- .posterior_group_at_x(
    x0 = q50,
    params = params,
    pi_k = pi_k,
    y_levels = y_levels,
    df_t = df_t
  )
  
  post_high <- .posterior_group_at_x(
    x0 = q90,
    params = params,
    pi_k = pi_k,
    y_levels = y_levels,
    df_t = df_t
  )
  
  post_tbl <- tibble(
    group = y_levels,
    post_low = as.numeric(post_low[y_levels]),
    post_mid = as.numeric(post_mid[y_levels]),
    post_high = as.numeric(post_high[y_levels]),
    prior = as.numeric(pi_k[y_levels])
  ) %>%
    pivot_wider(
      names_from = group,
      values_from = c(post_low, post_mid, post_high, prior),
      names_sep = "_"
    )
  
  pattern_tbl <- .classify_pattern(
    params = params,
    post_high = post_high,
    post_low = post_low,
    pi_k = pi_k,
    boot_dir = boot_dir,
    y_levels = y_levels
  )
  
  params_wide <- params %>%
    select(group, n, mu, sigma) %>%
    pivot_wider(
      names_from = group,
      values_from = c(n, mu, sigma),
      names_sep = "_"
    )
  
  feature_row <- bind_cols(
    tibble(
      feature_model = feature,
      view_code = v_j,
      n_used = length(x),
      MI_bits = stat_obj$MI_bits,
      MI_nats = stat_obj$MI_nats,
      MI_norm = stat_obj$MI_norm,
      max_JS_bits = max_JS_bits,
      T_obs = T_obs,
      perm_stat = perm_stat,
      p_perm = p_perm,
      q10 = q10,
      q50 = q50,
      q90 = q90,
      error = NA_character_
    ),
    params_wide,
    JS_wide,
    delta_wide,
    top_prob_wide,
    post_tbl,
    pattern_tbl
  )
  
  list(
    feature_row = feature_row,
    perm_tbl = perm_tbl
  )
}

## ============================================================
## 2) FUNCIÓN PRINCIPAL
## ============================================================

run_bidirectional_explainability <- function(X_S,
                                             y,
                                             v,
                                             y_levels = c("NP", "FD", "SO"),
                                             df_t = 4,
                                             n_grid = 2048,
                                             B_perm = 999,
                                             B_boot = 499,
                                             scale_method = c("zscore", "robust", "none"),
                                             perm_stat = c("MI", "maxJS"),
                                             outdir = NULL,
                                             seed = 123) {
  
  scale_method <- match.arg(scale_method)
  perm_stat <- match.arg(perm_stat)
  
  X_S <- as.data.frame(X_S, check.names = FALSE)
  
  if (is.null(colnames(X_S))) {
    stop("X_S debe tener nombres de columnas.")
  }
  
  y <- factor(y, levels = y_levels)
  
  if (length(y) != nrow(X_S)) {
    stop("length(y) debe ser igual a nrow(X_S).")
  }
  
  if (any(is.na(y))) {
    stop("y contiene valores fuera de y_levels o NA.")
  }
  
  feature_names <- colnames(X_S)
  v <- .prepare_views(v, feature_names)
  
  cat("\n============================================================\n")
  cat("FRAMEWORK EXPLICATIVO BIDIRECCIONAL MULTI-ÓMICO\n")
  cat("============================================================\n")
  cat("Muestras:", nrow(X_S), "\n")
  cat("Features:", ncol(X_S), "\n")
  cat("Grupos:\n")
  print(table(y))
  cat("Vistas:\n")
  print(table(v))
  cat("B_perm:", B_perm, "| B_boot:", B_boot, "\n")
  cat("perm_stat:", perm_stat, "\n")
  cat("scale_method:", scale_method, "\n")
  cat("============================================================\n")
  
  results <- vector("list", length(feature_names))
  names(results) <- feature_names
  
  for (jj in seq_along(feature_names)) {
    
    feat <- feature_names[jj]
    
    cat(
      sprintf(
        "[%d/%d] %s\n",
        jj,
        length(feature_names),
        feat
      )
    )
    
    results[[jj]] <- .compute_one_feature(
      feature = feat,
      x = X_S[[feat]],
      y = y,
      v_j = v[[feat]],
      y_levels = y_levels,
      df_t = df_t,
      n_grid = n_grid,
      B_perm = B_perm,
      B_boot = B_boot,
      scale_method = scale_method,
      seed = seed + jj,
      perm_stat = perm_stat
    )
  }
  
  feature_signature <- bind_rows(lapply(results, `[[`, "feature_row")) %>%
    mutate(
      q_perm = p.adjust(p_perm, method = "BH"),
      evidence_grade = .evidence_grade(
        p_perm = p_perm,
        q_perm = q_perm,
        MI_norm = MI_norm,
        top_location_prob_boot = top_location_prob_boot,
        concordant_high = concordant_high,
        posterior_enrichment_high = posterior_enrichment_high,
        pattern = bidirectional_pattern
      )
    ) %>%
    arrange(
      factor(
        evidence_grade,
        levels = c(
          "A_strong",
          "B_moderate",
          "C_exploratory",
          "D_ambiguous_dependent",
          "E_no_evidence"
        )
      ),
      q_perm,
      desc(MI_norm),
      desc(max_JS_bits)
    )
  
  final_bidirectional_table <- .make_final_bidirectional_table(
    feature_signature = feature_signature,
    y_levels = y_levels,
    remove_prefix = TRUE
  )
  
  
  
  permutation_distribution <- bind_rows(lapply(results, `[[`, "perm_tbl"))
  
  view_summary <- feature_signature %>%
    group_by(view_code) %>%
    summarise(
      n_features = n(),
      n_A = sum(evidence_grade == "A_strong", na.rm = TRUE),
      n_B = sum(evidence_grade == "B_moderate", na.rm = TRUE),
      n_C = sum(evidence_grade == "C_exploratory", na.rm = TRUE),
      n_D_ambiguous = sum(evidence_grade == "D_ambiguous_dependent", na.rm = TRUE),
      n_E_no_evidence = sum(evidence_grade == "E_no_evidence", na.rm = TRUE),
      
      mean_MI_bits = mean(MI_bits, na.rm = TRUE),
      median_MI_bits = median(MI_bits, na.rm = TRUE),
      sum_MI_bits = sum(MI_bits, na.rm = TRUE),
      mean_MI_norm = mean(MI_norm, na.rm = TRUE),
      median_max_JS_bits = median(max_JS_bits, na.rm = TRUE),
      n_q_lt_0_10 = sum(q_perm < 0.10, na.rm = TRUE),
      n_q_lt_0_20 = sum(q_perm < 0.20, na.rm = TRUE),
      prop_concordant = mean(concordant_high, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(
      desc(n_A),
      desc(n_B),
      desc(n_C),
      desc(mean_MI_norm),
      desc(sum_MI_bits)
    )
  
  if (!is.null(outdir)) {
    
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
    
    write.csv(
      feature_signature,
      file.path(outdir, "tables", "01_bidirectional_feature_signature.csv"),
      row.names = FALSE
    )
    
    write.csv(
      view_summary,
      file.path(outdir, "tables", "02_bidirectional_view_summary.csv"),
      row.names = FALSE
    )
    
    write.csv(
      permutation_distribution,
      file.path(outdir, "tables", "03_permutation_distribution.csv"),
      row.names = FALSE
    )
    write.csv(
      final_bidirectional_table,
      file.path(outdir, "tables", "04_bidirectional_final_table_corrected.csv"),
      row.names = FALSE
    )
    ## Plot 1: MI por feature
    p_mi <- feature_signature %>%
      mutate(
        feature_model = factor(feature_model, levels = rev(feature_model)),
        evidence_grade = factor(
          evidence_grade,
          levels = c(
            "A_strong",
            "B_moderate",
            "C_exploratory",
            "D_ambiguous_dependent",
            "E_no_evidence"
          )
        )
      ) %>%
      ggplot(aes(x = feature_model, y = MI_bits, fill = evidence_grade)) +
      geom_col() +
      coord_flip() +
      facet_grid(view_code ~ ., scales = "free_y", space = "free_y") +
      labs(
        title = "Información mutua feature–grupo",
        subtitle = "I(X_j;Y). Mayor valor = mayor dependencia entre feature y grupo",
        x = "Feature",
        y = "MI bits",
        fill = "Evidencia"
      ) +
      theme_minimal()
    
    ggsave(
      file.path(outdir, "plots", "01_MI_feature_group.png"),
      p_mi,
      width = 10,
      height = max(6, 0.18 * nrow(feature_signature)),
      dpi = 300
    )
    
    ## Plot 2: JS por contraste
    js_cols <- grep("^JS_bits_", colnames(feature_signature), value = TRUE)
    
    if (length(js_cols) > 0) {
      p_js <- feature_signature %>%
        select(feature_model, view_code, evidence_grade, all_of(js_cols)) %>%
        pivot_longer(
          cols = all_of(js_cols),
          names_to = "contrast",
          values_to = "JS_bits"
        ) %>%
        mutate(
          contrast = sub("^JS_bits_", "", contrast),
          feature_model = factor(feature_model, levels = rev(unique(feature_signature$feature_model)))
        ) %>%
        ggplot(aes(x = feature_model, y = JS_bits, fill = contrast)) +
        geom_col(position = "dodge") +
        coord_flip() +
        facet_grid(view_code ~ ., scales = "free_y", space = "free_y") +
        labs(
          title = "Distancia Jensen–Shannon entre grupos",
          subtitle = "Mayor valor = mayor separación de distribuciones entre pares de grupos",
          x = "Feature",
          y = "JS bits",
          fill = "Contraste"
        ) +
        theme_minimal()
      
      ggsave(
        file.path(outdir, "plots", "02_JS_pairwise_group_distances.png"),
        p_js,
        width = 11,
        height = max(6, 0.18 * nrow(feature_signature)),
        dpi = 300
      )
    }
    
    ## Plot 3: posterior del grupo para feature alta
    post_high_cols <- grep("^post_high_", colnames(feature_signature), value = TRUE)
    
    if (length(post_high_cols) > 0) {
      p_post_high <- feature_signature %>%
        select(feature_model, view_code, evidence_grade, all_of(post_high_cols)) %>%
        pivot_longer(
          cols = all_of(post_high_cols),
          names_to = "group",
          values_to = "posterior"
        ) %>%
        mutate(
          group = sub("^post_high_", "", group),
          feature_model = factor(feature_model, levels = rev(unique(feature_signature$feature_model)))
        ) %>%
        ggplot(aes(x = feature_model, y = posterior, fill = group)) +
        geom_col(position = "dodge") +
        coord_flip() +
        facet_grid(view_code ~ ., scales = "free_y", space = "free_y") +
        labs(
          title = "P(grupo | feature alta)",
          subtitle = "Posterior calculado en el percentil 90 de cada feature",
          x = "Feature",
          y = "Probabilidad posterior",
          fill = "Grupo"
        ) +
        theme_minimal()
      
      ggsave(
        file.path(outdir, "plots", "03_posterior_group_given_high_feature.png"),
        p_post_high,
        width = 11,
        height = max(6, 0.18 * nrow(feature_signature)),
        dpi = 300
      )
    }
    
    ## Plot 4: resumen por vista
    p_view <- view_summary %>%
      pivot_longer(
        cols = c(
          n_A,
          n_B,
          n_C,
          n_D_ambiguous,
          n_E_no_evidence
        ),
        names_to = "grade",
        values_to = "n"
      ) %>%
      mutate(
        grade = recode(
          grade,
          n_A = "A_strong",
          n_B = "B_moderate",
          n_C = "C_exploratory",
          n_D_ambiguous = "D_ambiguous_dependent",
          n_E_no_evidence = "E_no_evidence"
        ),
        grade = factor(
          grade,
          levels = c(
            "A_strong",
            "B_moderate",
            "C_exploratory",
            "D_ambiguous_dependent",
            "E_no_evidence"
          )
        )
      ) %>%
      ggplot(aes(x = view_code, y = n, fill = grade)) +
      geom_col() +
      labs(
        title = "Resumen multi-ómico por vista",
        subtitle = "Número de features por grado de evidencia",
        x = "Vista",
        y = "Número de features",
        fill = "Evidencia"
      ) +
      theme_minimal()
    
    ggsave(
      file.path(outdir, "plots", "04_view_summary_evidence_grade.png"),
      p_view,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
  
  list(
    feature_signature = feature_signature,
    final_bidirectional_table = final_bidirectional_table,
    view_summary = view_summary,
    permutation_distribution = permutation_distribution
  )
}


## ============================================================
## INPUT RECONSTRUIDO MOFA2 PARA FRAMEWORK BIDIRECCIONAL
## ============================================================

make_mofa_reconstructed_bidir_input <- function(scenario_dir,
                                                diag_target = "diagnostico_bayes_features",
                                                keep_factors = FALSE,
                                                save_reconstructed = TRUE) {
  
  scenario_dir <- normalizePath(path.expand(scenario_dir), mustWork = TRUE)
  
  diagnostic_rds <- file.path(
    scenario_dir,
    diag_target,
    "rds",
    "diagnostic_objects.rds"
  )
  
  diagnostic_rds <- normalizePath(path.expand(diagnostic_rds), mustWork = TRUE)
  
  diagnostic_objects <- readRDS(diagnostic_rds)
  
  ## ------------------------------------------------------------
  ## Helpers internos
  ## ------------------------------------------------------------
  
  resolve_path_from_diagnostic <- function(path_raw, diagnostic_rds) {
    
    path_raw <- path.expand(path_raw)
    
    if (file.exists(path_raw)) {
      return(normalizePath(path_raw, mustWork = TRUE))
    }
    
    repo_guess <- normalizePath(
      file.path(dirname(diagnostic_rds), "..", "..", "..", ".."),
      mustWork = FALSE
    )
    
    path_alt <- file.path(repo_guess, sub("^\\./", "", path_raw))
    
    if (file.exists(path_alt)) {
      return(normalizePath(path_alt, mustWork = TRUE))
    }
    
    stop(
      "No encuentro el archivo MOFA2:\n",
      path_raw,
      "\nRuta alternativa probada:\n",
      path_alt
    )
  }
  
  load_mofa_model_robust <- function(path) {
    
    path <- path.expand(path)
    
    if (!file.exists(path)) {
      stop("No existe el archivo del modelo MOFA2:\n", path)
    }
    
    ext <- tolower(tools::file_ext(path))
    
    if (ext %in% c("hdf5", "h5")) {
      return(MOFA2::load_model(path))
    }
    
    obj <- tryCatch(
      readRDS(path),
      error = function(e) NULL
    )
    
    if (!is.null(obj)) {
      return(obj)
    }
    
    MOFA2::load_model(path)
  }
  
  standardize_factor <- function(x) {
    x <- as.character(x)
    ifelse(
      grepl("^Factor", x),
      x,
      paste0("Factor", gsub("[^0-9]", "", x))
    )
  }
  
  standardize_view <- function(x) {
    
    x <- tolower(as.character(x))
    
    dplyr::case_when(
      grepl("trans", x) ~ "tx",
      grepl("prot",  x) ~ "pr",
      grepl("metab", x) ~ "me",
      grepl("clin",  x) ~ "cl",
      x %in% c("tx", "pr", "me", "cl") ~ x,
      TRUE ~ x
    )
  }
  
  ## ------------------------------------------------------------
  ## Extraer train/test ya proyectados
  ## ------------------------------------------------------------
  
  if (is.null(diagnostic_objects$df_train)) {
    stop("No encuentro diagnostic_objects$df_train")
  }
  
  if (is.null(diagnostic_objects$df_test)) {
    stop("No encuentro diagnostic_objects$df_test")
  }
  
  if (is.null(diagnostic_objects$df_weights)) {
    stop("No encuentro diagnostic_objects$df_weights")
  }
  
  df_train <- diagnostic_objects$df_train
  df_test  <- diagnostic_objects$df_test
  
  factor_cols <- grep("^Factor[0-9]+$", names(df_train), value = TRUE)
  factor_cols <- factor_cols[order(as.integer(gsub("^Factor", "", factor_cols)))]
  
  if (length(factor_cols) == 0) {
    stop("No encuentro columnas Factor1, Factor2, ... en df_train")
  }
  
  if (!all(factor_cols %in% names(df_test))) {
    stop("df_test no contiene los mismos factores que df_train")
  }
  
  feature_cols <- setdiff(names(df_train), c("y", factor_cols))
  
  if (!all(feature_cols %in% names(df_test))) {
    stop("df_test no contiene las mismas features que df_train")
  }
  
  ## ------------------------------------------------------------
  ## Cargar modelo MOFA2
  ## ------------------------------------------------------------
  
  if (is.null(diagnostic_objects$diagnostic_config$mofa_model_file)) {
    stop("No encuentro diagnostic_objects$diagnostic_config$mofa_model_file")
  }
  
  mofa_model_file <- resolve_path_from_diagnostic(
    diagnostic_objects$diagnostic_config$mofa_model_file,
    diagnostic_rds
  )
  
  mofa_model <- load_mofa_model_robust(mofa_model_file)
  
  ## ------------------------------------------------------------
  ## Mapa feature_model -> feature interna MOFA
  ## ------------------------------------------------------------
  
  needed_cols <- c("feature_model", "feature", "view_code")
  missing_cols <- setdiff(needed_cols, names(diagnostic_objects$df_weights))
  
  if (length(missing_cols) > 0) {
    stop(
      "Faltan columnas en diagnostic_objects$df_weights: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  feature_map <- diagnostic_objects$df_weights %>%
    dplyr::transmute(
      feature_model = as.character(feature_model),
      mofa_feature  = as.character(feature),
      view_code     = as.character(view_code)
    ) %>%
    dplyr::distinct() %>%
    dplyr::filter(feature_model %in% feature_cols)
  
  missing_map <- setdiff(feature_cols, feature_map$feature_model)
  
  if (length(missing_map) > 0) {
    stop(
      "Features sin mapa en df_weights:\n",
      paste(missing_map, collapse = ", ")
    )
  }
  
  ## ------------------------------------------------------------
  ## Pesos MOFA
  ## ------------------------------------------------------------
  
  weights_long <- MOFA2::get_weights(
    mofa_model,
    views = "all",
    factors = "all",
    as.data.frame = TRUE
  ) %>%
    dplyr::mutate(
      mofa_feature = as.character(feature),
      factor       = standardize_factor(factor),
      view_code    = standardize_view(view),
      weight       = as.numeric(value)
    ) %>%
    dplyr::filter(factor %in% factor_cols)
  
  ## ------------------------------------------------------------
  ## Reconstrucción: Xhat = Z %*% t(W)
  ## ------------------------------------------------------------
  
  reconstruct_dataset <- function(df, split_name) {
    
    Z <- as.matrix(df[, factor_cols, drop = FALSE])
    
    reconstruct_view <- function(vc) {
      
      map_v <- feature_map %>%
        dplyr::filter(view_code == vc)
      
      if (nrow(map_v) == 0) {
        return(NULL)
      }
      
      W_v <- weights_long %>%
        dplyr::filter(
          view_code == vc,
          mofa_feature %in% map_v$mofa_feature
        ) %>%
        dplyr::select(mofa_feature, factor, weight) %>%
        tidyr::pivot_wider(
          names_from = factor,
          values_from = weight,
          values_fill = 0
        )
      
      missing_weights <- setdiff(map_v$mofa_feature, W_v$mofa_feature)
      
      if (length(missing_weights) > 0) {
        stop(
          "Faltan pesos MOFA para vista ", vc, ":\n",
          paste(missing_weights, collapse = ", ")
        )
      }
      
      W_v <- W_v %>%
        dplyr::arrange(match(mofa_feature, map_v$mofa_feature))
      
      W_mat <- W_v %>%
        dplyr::select(dplyr::all_of(factor_cols)) %>%
        as.matrix()
      
      Xhat <- Z %*% t(W_mat)
      colnames(Xhat) <- map_v$feature_model
      
      Xhat
    }
    
    Xhat_list <- lapply(unique(feature_map$view_code), reconstruct_view)
    Xhat_list <- Filter(Negate(is.null), Xhat_list)
    
    Xhat_mat <- do.call(cbind, Xhat_list)
    Xhat_mat <- Xhat_mat[, feature_cols, drop = FALSE]
    
    sample_id <- rownames(df)
    
    if (
      is.null(sample_id) ||
      anyDuplicated(sample_id) ||
      all(sample_id == as.character(seq_len(nrow(df))))
    ) {
      sample_id <- paste0(split_name, "_", seq_len(nrow(df)))
    }
    
    if (keep_factors) {
      data.frame(
        sample_id = sample_id,
        split = split_name,
        y = df$y,
        df[, factor_cols, drop = FALSE],
        Xhat_mat,
        check.names = FALSE
      )
    } else {
      data.frame(
        sample_id = sample_id,
        split = split_name,
        y = df$y,
        Xhat_mat,
        check.names = FALSE
      )
    }
  }
  
  df_all_reconstructed <- dplyr::bind_rows(
    reconstruct_dataset(df_train, "train"),
    reconstruct_dataset(df_test,  "test")
  )
  
  ## ------------------------------------------------------------
  ## Guardar reconstrucción
  ## ------------------------------------------------------------
  
  outdir_recon <- file.path(
    scenario_dir,
    diag_target,
    "mofa_train_test_reconstructed"
  )
  
  if (save_reconstructed) {
    
    dir.create(outdir_recon, recursive = TRUE, showWarnings = FALSE)
    
    readr::write_csv(
      df_all_reconstructed,
      file.path(outdir_recon, "train_test_mofa_reconstructed_data.csv")
    )
    
    saveRDS(
      df_all_reconstructed,
      file.path(outdir_recon, "train_test_mofa_reconstructed_data.rds")
    )
  }
  
  ## ------------------------------------------------------------
  ## Devolver entrada directa para framework
  ## ------------------------------------------------------------
  
  X_S <- df_all_reconstructed %>%
    dplyr::select(dplyr::all_of(feature_cols))
  
  y <- df_all_reconstructed$y
  
  v <- feature_map$view_code
  names(v) <- feature_map$feature_model
  v <- v[colnames(X_S)]
  
  list(
    X_S = X_S,
    y = y,
    v = v,
    df_all_reconstructed = df_all_reconstructed,
    feature_map = feature_map,
    factor_cols = factor_cols,
    feature_cols = feature_cols,
    outdir_reconstructed = outdir_recon
  )
}



## ============================================================
## INPUT OBSERVADO PARA FRAMEWORK BIDIRECCIONAL
## ============================================================

make_observed_bidir_input <- function(scenario_dir,
                                      diag_target = "diagnostico_bayes_features",
                                      save_observed = TRUE) {
  
  scenario_dir <- normalizePath(path.expand(scenario_dir), mustWork = TRUE)
  
  diagnostic_rds <- file.path(
    scenario_dir,
    diag_target,
    "rds",
    "diagnostic_objects.rds"
  )
  
  diagnostic_rds <- normalizePath(path.expand(diagnostic_rds), mustWork = TRUE)
  
  diagnostic_objects <- readRDS(diagnostic_rds)
  
  if (is.null(diagnostic_objects$df_train)) {
    stop("No encuentro diagnostic_objects$df_train")
  }
  
  if (is.null(diagnostic_objects$df_test)) {
    stop("No encuentro diagnostic_objects$df_test")
  }
  
  df_train <- diagnostic_objects$df_train
  df_test  <- diagnostic_objects$df_test
  
  factor_cols <- grep("^Factor[0-9]+$", names(df_train), value = TRUE)
  factor_cols <- factor_cols[order(as.integer(gsub("^Factor", "", factor_cols)))]
  
  feature_cols <- setdiff(names(df_train), c("y", factor_cols))
  
  if (!all(feature_cols %in% names(df_test))) {
    stop("df_test no contiene las mismas features observadas que df_train")
  }
  
  make_split <- function(df, split_name) {
    
    sample_id <- rownames(df)
    
    if (
      is.null(sample_id) ||
      anyDuplicated(sample_id) ||
      all(sample_id == as.character(seq_len(nrow(df))))
    ) {
      sample_id <- paste0(split_name, "_", seq_len(nrow(df)))
    }
    
    data.frame(
      sample_id = sample_id,
      split = split_name,
      y = df$y,
      df[, feature_cols, drop = FALSE],
      check.names = FALSE
    )
  }
  
  df_all_observed <- dplyr::bind_rows(
    make_split(df_train, "train"),
    make_split(df_test,  "test")
  )
  
  outdir_obs <- file.path(
    scenario_dir,
    diag_target,
    "observed_train_test_signal"
  )
  
  if (save_observed) {
    
    dir.create(outdir_obs, recursive = TRUE, showWarnings = FALSE)
    
    readr::write_csv(
      df_all_observed,
      file.path(outdir_obs, "train_test_observed_data.csv")
    )
    
    saveRDS(
      df_all_observed,
      file.path(outdir_obs, "train_test_observed_data.rds")
    )
  }
  
  X_S <- df_all_observed %>%
    dplyr::select(dplyr::all_of(feature_cols))
  
  y <- df_all_observed$y
  
  ## Preferir mapa real desde df_weights; si no existe, usar prefijo tx/pr/me/cl
  if (
    !is.null(diagnostic_objects$df_weights) &&
    all(c("feature_model", "view_code") %in% names(diagnostic_objects$df_weights))
  ) {
    
    feature_map <- diagnostic_objects$df_weights %>%
      dplyr::transmute(
        feature_model = as.character(feature_model),
        view_code = as.character(view_code)
      ) %>%
      dplyr::distinct() %>%
      dplyr::filter(feature_model %in% feature_cols)
    
  } else {
    
    feature_map <- tibble::tibble(
      feature_model = feature_cols,
      view_code = sub("_.*$", "", feature_cols)
    )
  }
  
  missing_view <- setdiff(feature_cols, feature_map$feature_model)
  
  if (length(missing_view) > 0) {
    stop(
      "Features observadas sin vista asignada:\n",
      paste(missing_view, collapse = ", ")
    )
  }
  
  v <- feature_map$view_code
  names(v) <- feature_map$feature_model
  v <- v[colnames(X_S)]
  
  list(
    X_S = X_S,
    y = y,
    v = v,
    df_all_observed = df_all_observed,
    feature_map = feature_map,
    factor_cols = factor_cols,
    feature_cols = feature_cols,
    outdir_observed = outdir_obs
  )
}

plot_bidirectional_explainability <- function(bidir,
                                              X_S = NULL,
                                              outdir = "framework_bidireccional_firma",
                                              top_n = Inf,
                                              y_levels = c("NP", "FD", "SO"),
                                              remove_prefix = TRUE,
                                              network_rho_thr = 0.60,
                                              network_max_edges = 80,
                                              show_opposite_edges = FALSE,
                                              network_split_by_group = TRUE,
                                              network_group_max_edges = 25,
                                              network_group_rho_thr = 0.70,
                                              network_label_width = 28) {
  if (is.null(bidir$feature_signature)) {
    stop("bidir debe tener bidir$feature_signature.")
  }
  
  feature_signature <- tibble::as_tibble(bidir$feature_signature)
  
  PLOTS_DIR <- file.path(outdir, "plots_interpretables_bidireccional")
  TABLES_DIR <- file.path(outdir, "tables")
  
  dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
  
  out_plot <- function(x) file.path(PLOTS_DIR, x)
  out_table <- function(x) file.path(TABLES_DIR, x)
  
  grade_levels <- c(
    "A_strong",
    "B_moderate",
    "C_exploratory",
    "D_ambiguous_dependent",
    "E_no_evidence"
  )
  feature_signature <- feature_signature %>%
    mutate(
      evidence_grade = factor(
        evidence_grade,
        levels = grade_levels
      ),
      feature_label = if (remove_prefix) {
        sub("^[^_]+_", "", feature_model)
      } else {
        feature_model
      },
      feature_label = make.unique(feature_label),
      feature_plot = paste0(feature_label, " [", view_code, "]")
    ) %>%
    arrange(
      evidence_grade,
      q_perm,
      desc(MI_norm),
      desc(MI_bits),
      desc(max_JS_bits)
    )
  
  if (is.finite(top_n)) {
    feature_signature <- feature_signature %>%
      slice_head(n = top_n)
  }
  
  feature_levels <- rev(feature_signature$feature_plot)
  
  feature_signature <- feature_signature %>%
    mutate(
      feature_plot = factor(feature_plot, levels = feature_levels)
    )
  
  plot_manifest <- tibble(
    plot = character(),
    description = character()
  )
  
  add_manifest <- function(filename, description) {
    plot_manifest <<- bind_rows(
      plot_manifest,
      tibble(
        plot = filename,
        description = description
      )
    )
  }
  make_short_network_label <- function(x, width = network_label_width) {
    
    x <- as.character(x)
    x <- trimws(x)
    x[is.na(x) | x == "" | x == "NA"] <- NA_character_
    
    x <- ifelse(
      is.na(x),
      NA_character_,
      ifelse(
        nchar(x) > width,
        paste0(substr(x, 1, max(width - 3, 1)), "..."),
        x
      )
    )
    
    x
  }
  
  make_network_label_from_df <- function(df,
                                         fallback_col = "feature_plot",
                                         width = network_label_width) {
    
    preferred_cols <- c(
      "manual_biological_annotation",
      "database_semantic_category",
      "ontology_label",
      "term_name",
      "feature_label",
      fallback_col
    )
    
    preferred_cols <- intersect(preferred_cols, colnames(df))
    
    if (length(preferred_cols) == 0) {
      return(rep(NA_character_, nrow(df)))
    }
    
    out <- rep(NA_character_, nrow(df))
    
    for (cc in preferred_cols) {
      
      val <- as.character(df[[cc]])
      val <- trimws(val)
      val[is.na(val) | val == "" | val == "NA"] <- NA_character_
      
      idx <- is.na(out) & !is.na(val)
      out[idx] <- val[idx]
    }
    
    out <- make_short_network_label(out, width = width)
    out[is.na(out)] <- as.character(df[[fallback_col]])[is.na(out)]
    make.unique(out)
  }
  
  edge_type_label <- function(x) {
    dplyr::recode(
      as.character(x),
      posterior_high_group = "Directa",
      shared_FD_SO = "Compartida FD/SO",
      FD_SO_splitter = "Separadora FD/SO",
      ambiguous_posterior = "Ambigua",
      .default = as.character(x)
    )
  }
  
  feature_belongs_to_group <- function(df, group_interest) {
    
    pattern <- as.character(df$bidirectional_pattern)
    posterior <- as.character(df$posterior_high_group)
    
    if (group_interest == "NP") {
      return(
        posterior == "NP" |
          pattern == "NP_like_high"
      )
    }
    
    if (group_interest == "FD") {
      return(
        posterior == "FD" |
          pattern %in% c(
            "FD_like_high",
            "shared_FD_SO_high",
            "FD_gt_SO_splitter",
            "SO_gt_FD_splitter"
          )
      )
    }
    
    if (group_interest == "SO") {
      return(
        posterior == "SO" |
          pattern %in% c(
            "SO_like_high",
            "shared_FD_SO_high",
            "FD_gt_SO_splitter",
            "SO_gt_FD_splitter"
          )
      )
    }
    
    rep(FALSE, nrow(df))
  }
  
  phenotype_modules_for_group <- function(group_interest) {
    
    if (group_interest == "NP") {
      return(c("NP-like"))
    }
    
    if (group_interest == "FD") {
      return(c("FD-like", "FD/SO splitter", "shared FD/SO"))
    }
    
    if (group_interest == "SO") {
      return(c("SO-like", "FD/SO splitter", "shared FD/SO"))
    }
    
    character()
  }
  ## ============================================================
  ## PLOT 1. FEATURE | GRUPO
  ## Heatmap de localización robusta por grupo
  ## ============================================================
  
  mu_cols <- intersect(paste0("mu_", y_levels), colnames(feature_signature))
  
  if (length(mu_cols) > 0) {
    
    mu_long <- feature_signature %>%
      select(
        feature_model,
        feature_plot,
        view_code,
        evidence_grade,
        bidirectional_pattern,
        all_of(mu_cols)
      ) %>%
      pivot_longer(
        cols = all_of(mu_cols),
        names_to = "group",
        values_to = "mu"
      ) %>%
      mutate(
        group = sub("^mu_", "", group),
        group = factor(group, levels = y_levels)
      )
    
    p_mu <- ggplot(
      mu_long,
      aes(x = group, y = feature_plot, fill = mu)
    ) +
      geom_tile() +
      geom_text(
        aes(label = sprintf("%.2f", mu)),
        size = 3
      ) +
      facet_grid(
        view_code ~ .,
        scales = "free_y",
        space = "free_y"
      ) +
      labs(
        title = "Dirección grupo → feature",
        subtitle = "Localización robusta de cada feature por grupo. Valores en escala z.",
        x = "Grupo clínico",
        y = "Feature",
        fill = "Media/mediana\nrobusta"
      ) +
      theme_minimal() +
      theme(
        strip.text.y = element_text(angle = 0),
        axis.text.y = element_text(size = 8)
      )
    
    ggsave(
      out_plot("01_feature_given_group_heatmap.png"),
      p_mu,
      width = 8,
      height = max(6, 0.25 * nrow(feature_signature)),
      dpi = 300
    )
    
    add_manifest(
      "01_feature_given_group_heatmap.png",
      "Muestra en qué grupo está alta o baja cada feature: lectura feature | grupo."
    )
  }
  
  ## ============================================================
  ## PLOT 2. CONTRASTES DE GRUPO
  ## FD-NP, SO-NP, SO-FD con IC bootstrap
  ## ============================================================
  
  delta_cols <- grep(
    "^(estimate|ci_low|ci_high|prob_gt0|prob_lt0)__",
    colnames(feature_signature),
    value = TRUE
  )
  
  if (length(delta_cols) > 0) {
    
    delta_long <- feature_signature %>%
      select(
        feature_model,
        feature_plot,
        view_code,
        evidence_grade,
        bidirectional_pattern,
        all_of(delta_cols)
      ) %>%
      pivot_longer(
        cols = all_of(delta_cols),
        names_to = c("metric", "contrast"),
        names_sep = "__",
        values_to = "value"
      ) %>%
      pivot_wider(
        names_from = metric,
        values_from = value
      ) %>%
      mutate(
        contrast = recode(
          contrast,
          FD_minus_NP = "FD - NP",
          SO_minus_NP = "SO - NP",
          SO_minus_FD = "SO - FD"
        ),
        contrast = factor(
          contrast,
          levels = c("FD - NP", "SO - NP", "SO - FD")
        ),
        prob_label = paste0("P>0=", sprintf("%.2f", prob_gt0))
      )
    
    p_delta <- ggplot(
      delta_long,
      aes(
        x = estimate,
        y = feature_plot
      )
    ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      geom_errorbarh(
        aes(xmin = ci_low, xmax = ci_high),
        height = 0.15
      ) +
      geom_point(
        aes(shape = evidence_grade),
        size = 2
      ) +
      facet_grid(
        view_code ~ contrast,
        scales = "free_y",
        space = "free_y"
      ) +
      labs(
        title = "Contrastes grupo → feature",
        subtitle = "Diferencias robustas de localización con IC bootstrap. Escala z.",
        x = "Diferencia entre grupos",
        y = "Feature",
        shape = "Evidencia"
      ) +
      theme_minimal() +
      theme(
        strip.text.y = element_text(angle = 0),
        axis.text.y = element_text(size = 7)
      )
    
    ggsave(
      out_plot("02_group_contrasts_bootstrap_CI.png"),
      p_delta,
      width = 13,
      height = max(6, 0.25 * nrow(feature_signature)),
      dpi = 300
    )
    
    add_manifest(
      "02_group_contrasts_bootstrap_CI.png",
      "Contrastes FD-NP, SO-NP y SO-FD con IC bootstrap para explicar dirección diferencial."
    )
  }
  
  ## ============================================================
  ## PLOT 3. DEPENDENCIA FEATURE-GRUPO
  ## Información mutua y q-value por permutación
  ## ============================================================
  
  p_mi_q <- feature_signature %>%
    mutate(
      minus_log10_q = -log10(pmax(q_perm, 1e-300)),
      q_threshold_010 = -log10(0.10),
      q_threshold_020 = -log10(0.20)
    ) %>%
    ggplot(
      aes(
        x = MI_norm,
        y = minus_log10_q
      )
    ) +
    geom_hline(
      yintercept = -log10(0.10),
      linetype = "dashed"
    ) +
    geom_point(
      aes(
        shape = evidence_grade,
        size = max_JS_bits
      )
    ) +
    geom_text(
      aes(label = feature_label),
      size = 2.6,
      vjust = -0.6,
      check_overlap = TRUE
    ) +
    facet_wrap(
      ~ view_code,
      scales = "free"
    ) +
    labs(
      title = "Dependencia bidireccional feature–grupo",
      subtitle = "Información mutua normalizada vs significación por permutación. Línea: q = 0.10.",
      x = "Información mutua normalizada I(X;Y)/H(Y)",
      y = "-log10(q permutation)",
      shape = "Evidencia",
      size = "Máx. JS"
    ) +
    theme_minimal()
  
  ggsave(
    out_plot("03_mutual_information_vs_qvalue.png"),
    p_mi_q,
    width = 11,
    height = 7,
    dpi = 300
  )
  
  add_manifest(
    "03_mutual_information_vs_qvalue.png",
    "Resume la fuerza estadística de dependencia feature–grupo usando información mutua y q-value por permutación."
  )
  
  ## ============================================================
  ## PLOT 4. DISTANCIAS JENSEN-SHANNON POR CONTRASTE
  ## ============================================================
  
  js_cols <- grep("^JS_bits_", colnames(feature_signature), value = TRUE)
  
  if (length(js_cols) > 0) {
    
    js_long <- feature_signature %>%
      select(
        feature_model,
        feature_plot,
        view_code,
        evidence_grade,
        bidirectional_pattern,
        all_of(js_cols)
      ) %>%
      pivot_longer(
        cols = all_of(js_cols),
        names_to = "contrast",
        values_to = "JS_bits"
      ) %>%
      mutate(
        contrast = sub("^JS_bits_", "", contrast),
        contrast = recode(
          contrast,
          NP_vs_FD = "NP vs FD",
          NP_vs_SO = "NP vs SO",
          FD_vs_SO = "FD vs SO"
        )
      )
    
    p_js_heat <- ggplot(
      js_long,
      aes(
        x = contrast,
        y = feature_plot,
        fill = JS_bits
      )
    ) +
      geom_tile() +
      geom_text(
        aes(label = sprintf("%.3f", JS_bits)),
        size = 2.7
      ) +
      facet_grid(
        view_code ~ .,
        scales = "free_y",
        space = "free_y"
      ) +
      labs(
        title = "Distancia entre distribuciones de grupo",
        subtitle = "Jensen–Shannon por contraste. Mayor valor = mayor separación distribucional.",
        x = "Contraste",
        y = "Feature",
        fill = "JS bits"
      ) +
      theme_minimal() +
      theme(
        strip.text.y = element_text(angle = 0),
        axis.text.y = element_text(size = 8)
      )
    
    ggsave(
      out_plot("04_JS_distance_heatmap.png"),
      p_js_heat,
      width = 9,
      height = max(6, 0.25 * nrow(feature_signature)),
      dpi = 300
    )
    
    add_manifest(
      "04_JS_distance_heatmap.png",
      "Indica qué par de grupos separa mejor cada feature mediante distancia Jensen–Shannon."
    )
  }
  
  ## ============================================================
  ## PLOT 5. GRUPO | FEATURE ALTA
  ## P(NP/FD/SO | feature en percentil 90)
  ## ============================================================
  
  post_high_cols <- grep("^post_high_", colnames(feature_signature), value = TRUE)
  
  if (length(post_high_cols) > 0) {
    
    post_high_long <- feature_signature %>%
      select(
        feature_model,
        feature_plot,
        view_code,
        evidence_grade,
        bidirectional_pattern,
        all_of(post_high_cols)
      ) %>%
      pivot_longer(
        cols = all_of(post_high_cols),
        names_to = "group",
        values_to = "posterior"
      ) %>%
      mutate(
        group = sub("^post_high_", "", group),
        group = factor(group, levels = y_levels)
      )
    
    p_post_high <- ggplot(
      post_high_long,
      aes(
        x = group,
        y = feature_plot,
        fill = posterior
      )
    ) +
      geom_tile() +
      geom_text(
        aes(label = sprintf("%.2f", posterior)),
        size = 3
      ) +
      facet_grid(
        view_code ~ .,
        scales = "free_y",
        space = "free_y"
      ) +
      labs(
        title = "Dirección feature → grupo",
        subtitle = "Probabilidad posterior del grupo cuando la feature está alta, evaluada en el percentil 90.",
        x = "Grupo clínico",
        y = "Feature",
        fill = "P(grupo | feature alta)"
      ) +
      theme_minimal() +
      theme(
        strip.text.y = element_text(angle = 0),
        axis.text.y = element_text(size = 8)
      )
    
    ggsave(
      out_plot("05_posterior_group_given_high_feature_heatmap.png"),
      p_post_high,
      width = 8.5,
      height = max(6, 0.25 * nrow(feature_signature)),
      dpi = 300
    )
    
    add_manifest(
      "05_posterior_group_given_high_feature_heatmap.png",
      "Lectura inversa grupo | feature: muestra hacia qué grupo empuja una feature alta."
    )
  }
  
  
  ## ============================================================
  ## PLOT 5B. HEATMAP BIDIRECCIONAL DOBLE
  ## feature | grupo  +  grupo | feature alta
  ## ------------------------------------------------------------
  ## Panel izquierdo:
  ##   mu_NP, mu_FD, mu_SO
  ##
  ## Panel derecho:
  ##   P(NP | feature alta), P(FD | feature alta), P(SO | feature alta)
  ##
  ## Nota:
  ##   Para poder usar una misma escala visual:
  ##   - mu se reescala 0-1 dentro de cada feature.
  ##   - posterior ya está entre 0-1.
  ##   Los números escritos en las celdas son los valores reales.
  ## ============================================================
  
  mu_cols <- intersect(paste0("mu_", y_levels), colnames(feature_signature))
  post_high_cols <- grep("^post_high_", colnames(feature_signature), value = TRUE)
  
  if (length(mu_cols) > 0 && length(post_high_cols) > 0) {
    
    mu_bidir_long <- feature_signature %>%
      select(
        feature_model,
        feature_plot,
        view_code,
        evidence_grade,
        bidirectional_pattern,
        all_of(mu_cols)
      ) %>%
      pivot_longer(
        cols = all_of(mu_cols),
        names_to = "group",
        values_to = "raw_value"
      ) %>%
      mutate(
        group = sub("^mu_", "", group),
        group = factor(group, levels = y_levels),
        bidir_panel = "feature | grupo\nμ robusta",
        label_value = sprintf("%.2f", raw_value)
      ) %>%
      group_by(feature_plot) %>%
      mutate(
        fill_value = ifelse(
          max(raw_value, na.rm = TRUE) > min(raw_value, na.rm = TRUE),
          (raw_value - min(raw_value, na.rm = TRUE)) /
            (max(raw_value, na.rm = TRUE) - min(raw_value, na.rm = TRUE)),
          0.5
        )
      ) %>%
      ungroup()
    
    post_bidir_long <- feature_signature %>%
      select(
        feature_model,
        feature_plot,
        view_code,
        evidence_grade,
        bidirectional_pattern,
        all_of(post_high_cols)
      ) %>%
      pivot_longer(
        cols = all_of(post_high_cols),
        names_to = "group",
        values_to = "raw_value"
      ) %>%
      mutate(
        group = sub("^post_high_", "", group),
        group = factor(group, levels = y_levels),
        bidir_panel = "grupo | feature alta\nP(grupo | q90)",
        fill_value = raw_value,
        label_value = sprintf("%.2f", raw_value)
      )
    
    bidir_heat_long <- bind_rows(
      mu_bidir_long,
      post_bidir_long
    ) %>%
      mutate(
        bidir_panel = factor(
          bidir_panel,
          levels = c(
            "feature | grupo\nμ robusta",
            "grupo | feature alta\nP(grupo | q90)"
          )
        )
      )
    
    p_bidir_double <- ggplot(
      bidir_heat_long,
      aes(
        x = group,
        y = feature_plot,
        fill = fill_value
      )
    ) +
      geom_tile() +
      geom_text(
        aes(label = label_value),
        size = 3
      ) +
      facet_grid(
        view_code ~ bidir_panel,
        scales = "free_y",
        space = "free_y"
      ) +
      labs(
        title = "Mapa bidireccional de la firma multi-ómica",
        subtitle = "Izquierda: distribución feature | grupo. Derecha: posterior grupo | feature alta.",
        x = "Grupo clínico",
        y = "Feature",
        fill = "Intensidad\nrelativa"
      ) +
      theme_minimal() +
      theme(
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8)
      )
    
    ggsave(
      out_plot("05b_bidirectional_double_heatmap_feature_group_and_group_feature.png"),
      p_bidir_double,
      width = 11,
      height = max(6, 0.28 * nrow(feature_signature)),
      dpi = 300
    )
    
    add_manifest(
      "05b_bidirectional_double_heatmap_feature_group_and_group_feature.png",
      "Heatmap doble que muestra simultáneamente feature | grupo y grupo | feature alta."
    )
  }
  
  
  ## ============================================================
  ## PLOT 6. BAJA / MEDIA / ALTA FEATURE → POSTERIOR DEL GRUPO
  ## Este es el más didáctico para bidireccionalidad.
  ## ============================================================
  
  post_cols <- grep(
    "^(post_low|post_mid|post_high)_",
    colnames(feature_signature),
    value = TRUE
  )
  
  if (length(post_cols) > 0) {
    
    post_long <- feature_signature %>%
      select(
        feature_model,
        feature_plot,
        view_code,
        evidence_grade,
        bidirectional_pattern,
        all_of(post_cols)
      ) %>%
      pivot_longer(
        cols = all_of(post_cols),
        names_to = c("feature_state", "group"),
        names_pattern = "^(post_low|post_mid|post_high)_(.*)$",
        values_to = "posterior"
      ) %>%
      mutate(
        feature_state = recode(
          feature_state,
          post_low = "Feature baja",
          post_mid = "Feature media",
          post_high = "Feature alta"
        ),
        feature_state = factor(
          feature_state,
          levels = c("Feature baja", "Feature media", "Feature alta")
        ),
        group = factor(group, levels = y_levels)
      )
    
    p_post_states <- ggplot(
      post_long,
      aes(
        x = feature_state,
        y = posterior,
        group = group,
        linetype = group,
        shape = group
      )
    ) +
      geom_line() +
      geom_point(size = 2) +
      facet_wrap(
        ~ feature_plot,
        scales = "free_y"
      ) +
      labs(
        title = "Trayectoria posterior grupo | feature",
        subtitle = "Cómo cambia P(grupo) cuando la feature pasa de baja a media y alta.",
        x = NULL,
        y = "P(grupo | valor de feature)",
        linetype = "Grupo",
        shape = "Grupo"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 7)
      )
    
    ggsave(
      out_plot("06_posterior_trajectories_low_mid_high.png"),
      p_post_states,
      width = 13,
      height = max(8, 0.45 * ceiling(nrow(feature_signature) / 3)),
      dpi = 300
    )
    
    add_manifest(
      "06_posterior_trajectories_low_mid_high.png",
      "Plot didáctico de bidireccionalidad: muestra cómo cambia P(NP/FD/SO) al aumentar cada feature."
    )
  }
  
  ## ============================================================
  ## PLOT 7. CONCORDANCIA BIDIRECCIONAL
  ## Compara:
  ##   - grupo con mayor localización feature | grupo
  ##   - grupo favorecido por feature alta grupo | feature
  ## ============================================================
  
  needed_conc <- c(
    "top_location_group",
    "posterior_high_group",
    "concordant_high",
    "bidirectional_pattern"
  )
  
  if (all(needed_conc %in% colnames(feature_signature))) {
    
    conc_long <- feature_signature %>%
      select(
        feature_model,
        feature_plot,
        view_code,
        evidence_grade,
        bidirectional_pattern,
        top_location_group,
        posterior_high_group,
        concordant_high
      ) %>%
      pivot_longer(
        cols = c(top_location_group, posterior_high_group),
        names_to = "direction",
        values_to = "group"
      ) %>%
      mutate(
        direction = recode(
          direction,
          top_location_group = "Máx. feature | grupo",
          posterior_high_group = "Máx. grupo | feature alta"
        ),
        direction = factor(
          direction,
          levels = c(
            "Máx. feature | grupo",
            "Máx. grupo | feature alta"
          )
        ),
        group = factor(group, levels = y_levels),
        concordance_label = ifelse(concordant_high, "concordante", "discordante")
      )
    
    p_conc <- ggplot(
      conc_long,
      aes(
        x = direction,
        y = feature_plot,
        fill = group
      )
    ) +
      geom_tile() +
      geom_text(
        aes(label = as.character(group)),
        size = 3
      ) +
      facet_grid(
        view_code ~ .,
        scales = "free_y",
        space = "free_y"
      ) +
      labs(
        title = "Concordancia bidireccional",
        subtitle = "Compara el grupo donde la feature está más alta con el grupo favorecido por una feature alta.",
        x = NULL,
        y = "Feature",
        fill = "Grupo"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 20, hjust = 1),
        strip.text.y = element_text(angle = 0),
        axis.text.y = element_text(size = 8)
      )
    
    ggsave(
      out_plot("07_bidirectional_concordance.png"),
      p_conc,
      width = 9,
      height = max(6, 0.25 * nrow(feature_signature)),
      dpi = 300
    )
    
    add_manifest(
      "07_bidirectional_concordance.png",
      "Plot central del framework: resume si feature | grupo y grupo | feature alta apuntan al mismo grupo."
    )
  }
  
  ## ============================================================
  ## PLOT 8. PATRONES BIDIRECCIONALES
  ## Conteo de FD-like, SO-like, NP-like, splitter, etc.
  ## ============================================================
  
  if ("bidirectional_pattern" %in% colnames(feature_signature)) {
    
    p_patterns <- feature_signature %>%
      count(
        view_code,
        bidirectional_pattern,
        evidence_grade,
        name = "n"
      ) %>%
      ggplot(
        aes(
          x = bidirectional_pattern,
          y = n,
          fill = evidence_grade
        )
      ) +
      geom_col() +
      facet_wrap(~ view_code, scales = "free_y") +
      labs(
        title = "Patrones bidireccionales por vista",
        subtitle = "Clasificación interpretable de las features seleccionadas.",
        x = "Patrón bidireccional",
        y = "Número de features",
        fill = "Evidencia"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    ggsave(
      out_plot("08_bidirectional_patterns_by_view.png"),
      p_patterns,
      width = 11,
      height = 6,
      dpi = 300
    )
    
    add_manifest(
      "08_bidirectional_patterns_by_view.png",
      "Resume cuántas features FD-like, SO-like, NP-like, compartidas o inciertas hay por vista ómica."
    )
  }
  
  ## ============================================================
  ## PLOT 9. RESUMEN MULTI-ÓMICO POR VISTA
  ## ============================================================
  
  if (!is.null(bidir$view_summary)) {
    
    view_summary <- tibble::as_tibble(bidir$view_summary)
    
    grade_cols <- intersect(
      c(
        "n_A",
        "n_B",
        "n_C",
        "n_D_ambiguous",
        "n_E_no_evidence"
      ),
      colnames(view_summary)
    )
    
    if (length(grade_cols) > 0) {
      
      p_view_grade <- view_summary %>%
        select(view_code, all_of(grade_cols)) %>%
        pivot_longer(
          cols = all_of(grade_cols),
          names_to = "grade",
          values_to = "n"
        ) %>%
        mutate(
          grade = recode(
            grade,
            n_A = "A_strong",
            n_B = "B_moderate",
            n_C = "C_exploratory",
            n_D_ambiguous = "D_ambiguous_dependent",
            n_E_no_evidence = "E_no_evidence"
          ),
          grade = factor(grade, levels = grade_levels)
        ) %>%
        ggplot(
          aes(
            x = view_code,
            y = n,
            fill = grade
          )
        ) +
        geom_col() +
        labs(
          title = "Resumen multi-ómico de evidencia",
          subtitle = "Número de features por grado de evidencia en cada vista.",
          x = "Vista ómica",
          y = "Número de features",
          fill = "Evidencia"
        ) +
        theme_minimal()
      
      ggsave(
        out_plot("09_view_summary_evidence_grade.png"),
        p_view_grade,
        width = 8,
        height = 6,
        dpi = 300
      )
      
      add_manifest(
        "09_view_summary_evidence_grade.png",
        "Resume la contribución explicativa por vista ómica."
      )
    }
    
    if (all(c("view_code", "mean_MI_norm", "sum_MI_bits") %in% colnames(view_summary))) {
      
      p_view_mi <- view_summary %>%
        ggplot(
          aes(
            x = reorder(view_code, mean_MI_norm),
            y = mean_MI_norm
          )
        ) +
        geom_col() +
        coord_flip() +
        labs(
          title = "Intensidad media de dependencia por vista",
          subtitle = "Promedio de información mutua normalizada por feature.",
          x = "Vista ómica",
          y = "Media I(X;Y)/H(Y)"
        ) +
        theme_minimal()
      
      ggsave(
        out_plot("10_view_mean_mutual_information.png"),
        p_view_mi,
        width = 8,
        height = 5,
        dpi = 300
      )
      
      add_manifest(
        "10_view_mean_mutual_information.png",
        "Compara qué vista concentra mayor dependencia media feature–grupo."
      )
    }
  }
  
  ## ============================================================
  ## PLOT 11. FAMILIA FEATURE-GRUPO POR CONTEXTO
  ## ------------------------------------------------------------
  ## Objetivo:
  ##   Generar figuras separadas para:
  ##
  ##   Contextos unarios:
  ##     11_feature_group_NP
  ##     11_feature_group_FD
  ##     11_feature_group_SO
  ##
  ##   Contextos binarios:
  ##     11_feature_group_NP_vs_FD
  ##     11_feature_group_NP_vs_SO
  ##     11_feature_group_FD_vs_SO
  ##
  ## Lógica:
  ##   - No se grafican correlaciones feature-feature.
  ##   - No se mezclan todos los grupos en una sola red.
  ##   - Cada figura responde una pregunta concreta.
  ##   - No se escriben tablas nuevas para esta familia.
  ##
  ## Nodos:
  ##   - biomarcadores
  ##   - grupo(s) ancla del contexto
  ##
  ## Aristas:
  ##   - biomarcador -> grupo favorecido
  ##
  ## Tamaño nodo:
  ##   - MI_norm
  ##
  ## Color nodo:
  ##   - view_code
  ##
  ## Tipo de línea:
  ##   - grupo-like
  ##   - shared FD/SO
  ##   - splitter FD/SO
  ##   - ambigua, solo si se decide incluir D
  ## ============================================================
  
  required_network_cols <- c(
    "feature_model",
    "feature_plot",
    "feature_label",
    "view_code",
    "evidence_grade",
    "bidirectional_pattern",
    "posterior_high_group",
    "posterior_enrichment_high",
    "MI_norm",
    "max_JS_bits"
  )
  
  if (all(required_network_cols %in% colnames(feature_signature))) {
    
    ## ------------------------------------------------------------
    ## Parámetros internos de visualización.
    ## No cambian el framework estadístico.
    ## ------------------------------------------------------------
    
    plot11_grade_keep <- c(
      "A_strong",
      "B_moderate",
      "C_exploratory"
    )
    
    plot11_include_ambiguous_D <- FALSE
    
    if (isTRUE(plot11_include_ambiguous_D)) {
      plot11_grade_keep <- c(
        plot11_grade_keep,
        "D_ambiguous_dependent"
      )
    }
    
    plot11_max_unary  <- min(12L, as.integer(network_group_max_edges))
    plot11_max_binary <- min(16L, as.integer(network_group_max_edges))
    
    ## ------------------------------------------------------------
    ## Helpers internos de PLOT 11.
    ## ------------------------------------------------------------
    
    plot11_safe_num <- function(df, nm, default = NA_real_) {
      
      if (
        length(nm) != 1 ||
        is.na(nm) ||
        !nzchar(nm) ||
        !nm %in% colnames(df)
      ) {
        return(rep(default, nrow(df)))
      }
      
      suppressWarnings(as.numeric(df[[nm]]))
    }
    
    plot11_pattern_short <- function(x) {
      
      dplyr::recode(
        as.character(x),
        NP_like_high = "NP-like",
        FD_like_high = "FD-like",
        SO_like_high = "SO-like",
        shared_FD_SO_high = "shared FD/SO",
        FD_gt_SO_splitter = "FD>SO splitter",
        SO_gt_FD_splitter = "SO>FD splitter",
        uncertain = "ambiguous",
        .default = as.character(x)
      )
    }
    
    plot11_edge_relation <- function(x) {
      
      dplyr::case_when(
        x %in% c(
          "NP_like_high",
          "FD_like_high",
          "SO_like_high"
        ) ~ "Grupo-like",
        
        x == "shared_FD_SO_high" ~ "Shared FD/SO",
        
        x %in% c(
          "FD_gt_SO_splitter",
          "SO_gt_FD_splitter"
        ) ~ "Splitter FD/SO",
        
        TRUE ~ "Ambigua"
      )
    }
    
    plot11_axis_label <- function(pattern, context_id) {
      
      dplyr::case_when(
        pattern == "NP_like_high" ~ "NP-like axis",
        pattern == "FD_like_high" ~ "FD-like axis",
        pattern == "SO_like_high" ~ "SO-like axis",
        
        pattern == "shared_FD_SO_high" &
          context_id %in% c("NP_vs_FD", "NP_vs_SO") ~
          "shared FD/SO\nnon-NP axis",
        
        pattern == "shared_FD_SO_high" ~
          "shared FD/SO axis",
        
        pattern == "FD_gt_SO_splitter" ~
          "FD>SO splitter axis",
        
        pattern == "SO_gt_FD_splitter" ~
          "SO>FD splitter axis",
        
        TRUE ~ "ambiguous axis"
      )
    }
    
    plot11_axis_id <- function(pattern, side_group, context_id) {
      
      paste(
        "axis",
        context_id,
        side_group,
        pattern,
        sep = "__"
      )
    }
    
    
    
    plot11_context_defs <- tibble::tibble(
      context_id = c(
        "NP",
        "FD",
        "SO",
        "NP_vs_FD",
        "NP_vs_SO",
        "FD_vs_SO"
      ),
      context_type = c(
        "unary",
        "unary",
        "unary",
        "binary",
        "binary",
        "binary"
      ),
      left_group = c(
        "NP",
        "FD",
        "SO",
        "NP",
        "NP",
        "FD"
      ),
      right_group = c(
        NA_character_,
        NA_character_,
        NA_character_,
        "FD",
        "SO",
        "SO"
      ),
      contrast_estimate_col = c(
        NA_character_,
        NA_character_,
        NA_character_,
        "estimate__FD_minus_NP",
        "estimate__SO_minus_NP",
        "estimate__SO_minus_FD"
      ),
      contrast_prob_gt_col = c(
        NA_character_,
        NA_character_,
        NA_character_,
        "prob_gt0__FD_minus_NP",
        "prob_gt0__SO_minus_NP",
        "prob_gt0__SO_minus_FD"
      ),
      contrast_prob_lt_col = c(
        NA_character_,
        NA_character_,
        NA_character_,
        "prob_lt0__FD_minus_NP",
        "prob_lt0__SO_minus_NP",
        "prob_lt0__SO_minus_FD"
      ),
      contrast_js_col = c(
        NA_character_,
        NA_character_,
        NA_character_,
        "JS_bits_NP_vs_FD",
        "JS_bits_NP_vs_SO",
        "JS_bits_FD_vs_SO"
      ),
      patterns_keep = list(
        c("NP_like_high"),
        c("FD_like_high"),
        c("SO_like_high"),
        c(
          "NP_like_high",
          "FD_like_high",
          "shared_FD_SO_high"
        ),
        c(
          "NP_like_high",
          "SO_like_high",
          "shared_FD_SO_high"
        ),
        c(
          "FD_like_high",
          "SO_like_high",
          "FD_gt_SO_splitter",
          "SO_gt_FD_splitter"
        )
      )
    )
    
    plot11_base <- feature_signature %>%
      mutate(
        evidence_grade_chr = as.character(evidence_grade),
        bidirectional_pattern_chr = as.character(bidirectional_pattern),
        posterior_high_group_chr = as.character(posterior_high_group),
        feature_label_chr = as.character(feature_label),
        feature_label_chr = ifelse(
          is.na(feature_label_chr) | feature_label_chr == "",
          as.character(feature_model),
          feature_label_chr
        ),
        feature_label_clean = paste0(
          feature_label_chr,
          " [",
          view_code,
          "]"
        ),
        feature_label_clean = make_short_network_label(
          feature_label_clean,
          width = network_label_width
        ),
        pattern_short = plot11_pattern_short(
          bidirectional_pattern_chr
        ),
        node_label = feature_label_clean
      ) %>%
      filter(
        evidence_grade_chr %in% plot11_grade_keep,
        bidirectional_pattern_chr != "uncertain",
        evidence_grade_chr != "E_no_evidence"
      )
    
    ## ------------------------------------------------------------
    ## Función interna para construir cada contexto de familia 11.
    ## ------------------------------------------------------------
    plot11_make_context <- function(context_id,
                                    context_type,
                                    left_group,
                                    right_group,
                                    patterns_keep,
                                    contrast_estimate_col = NA_character_,
                                    contrast_prob_gt_col = NA_character_,
                                    contrast_prob_lt_col = NA_character_,
                                    contrast_js_col = NA_character_) {
      
      df_ctx <- plot11_base %>%
        filter(
          bidirectional_pattern_chr %in% patterns_keep
        )
      
      if (nrow(df_ctx) == 0) {
        cat("\n[PLOT 11 DAG] Sin biomarcadores para contexto: ", context_id, "\n")
        return(invisible(NULL))
      }
      
      df_ctx$contrast_estimate <- plot11_safe_num(
        df_ctx,
        contrast_estimate_col
      )
      
      df_ctx$contrast_prob_gt <- plot11_safe_num(
        df_ctx,
        contrast_prob_gt_col
      )
      
      df_ctx$contrast_prob_lt <- plot11_safe_num(
        df_ctx,
        contrast_prob_lt_col
      )
      
      df_ctx$contrast_js <- plot11_safe_num(
        df_ctx,
        contrast_js_col
      )
      
      df_ctx <- df_ctx %>%
        mutate(
          side_group = dplyr::case_when(
            context_type == "unary" ~ left_group,
            
            bidirectional_pattern_chr == paste0(left_group, "_like_high") ~
              left_group,
            
            !is.na(right_group) &
              bidirectional_pattern_chr == paste0(right_group, "_like_high") ~
              right_group,
            
            context_id %in% c("NP_vs_FD", "NP_vs_SO") &
              bidirectional_pattern_chr == "shared_FD_SO_high" ~
              right_group,
            
            context_id == "FD_vs_SO" &
              bidirectional_pattern_chr == "FD_gt_SO_splitter" ~
              "FD",
            
            context_id == "FD_vs_SO" &
              bidirectional_pattern_chr == "SO_gt_FD_splitter" ~
              "SO",
            
            TRUE ~ posterior_high_group_chr
          ),
          
          edge_relation = plot11_edge_relation(
            bidirectional_pattern_chr
          ),
          
          axis_label = plot11_axis_label(
            bidirectional_pattern_chr,
            context_id
          ),
          
          direction_prob = dplyr::case_when(
            is.finite(contrast_estimate) &
              contrast_estimate >= 0 ~ contrast_prob_gt,
            
            is.finite(contrast_estimate) &
              contrast_estimate < 0 ~ contrast_prob_lt,
            
            TRUE ~ top_location_prob_boot
          ),
          
          context_score = dplyr::coalesce(
            contrast_js,
            abs(contrast_estimate),
            posterior_enrichment_high,
            MI_norm,
            max_JS_bits,
            0
          )
        ) %>%
        mutate(
          axis_id = plot11_axis_id(
            bidirectional_pattern_chr,
            side_group,
            context_id
          )
        )
      
      valid_groups <- unique(c(left_group, right_group))
      valid_groups <- valid_groups[
        !is.na(valid_groups) &
          nzchar(valid_groups)
      ]
      
      df_ctx <- df_ctx %>%
        filter(side_group %in% valid_groups)
      
      if (nrow(df_ctx) == 0) {
        cat("\n[PLOT 11 DAG] Sin biomarcadores válidos para contexto: ", context_id, "\n")
        return(invisible(NULL))
      }
      
      ## ------------------------------------------------------------
      ## Ranking y recorte.
      ## En binarios se intenta balancear ambos lados.
      ## ------------------------------------------------------------
      
      if (context_type == "unary") {
        
        df_ctx <- df_ctx %>%
          arrange(
            factor(
              evidence_grade_chr,
              levels = grade_levels
            ),
            desc(posterior_enrichment_high),
            desc(MI_norm),
            desc(top_location_prob_boot),
            feature_label_chr
          ) %>%
          slice_head(n = plot11_max_unary)
        
      } else {
        
        n_per_side <- ceiling(plot11_max_binary / 2)
        
        df_ctx <- df_ctx %>%
          arrange(
            factor(
              evidence_grade_chr,
              levels = grade_levels
            ),
            desc(context_score),
            desc(direction_prob),
            desc(MI_norm),
            feature_label_chr
          ) %>%
          group_by(side_group) %>%
          slice_head(n = n_per_side) %>%
          ungroup() %>%
          arrange(
            factor(
              side_group,
              levels = valid_groups
            ),
            factor(
              evidence_grade_chr,
              levels = grade_levels
            ),
            desc(context_score),
            desc(MI_norm),
            feature_label_chr
          )
      }
      
      if (nrow(df_ctx) == 0) {
        return(invisible(NULL))
      }
      
      ## ------------------------------------------------------------
      ## Layout DAG manual:
      ##   biomarcador -> eje/patrón -> grupo
      ## ------------------------------------------------------------
      
      axis_order <- df_ctx %>%
        group_by(
          side_group,
          axis_id,
          axis_label,
          edge_relation
        ) %>%
        summarise(
          n_features = n(),
          mean_score = mean(context_score, na.rm = TRUE),
          mean_MI_norm = mean(MI_norm, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(
          factor(side_group, levels = valid_groups),
          factor(
            edge_relation,
            levels = c(
              "Grupo-like",
              "Shared FD/SO",
              "Splitter FD/SO",
              "Ambigua"
            )
          ),
          desc(mean_score),
          desc(mean_MI_norm),
          axis_label
        )
      
      df_ctx <- df_ctx %>%
        mutate(
          axis_id = factor(
            axis_id,
            levels = axis_order$axis_id
          )
        ) %>%
        arrange(
          axis_id,
          factor(
            evidence_grade_chr,
            levels = grade_levels
          ),
          desc(context_score),
          desc(MI_norm),
          feature_label_chr
        ) %>%
        mutate(
          y = rev(seq_len(n())),
          feature_x = 0.18,
          axis_x = 0.58,
          group_x = 0.92,
          feature_label_x = 0.145,
          feature_label_hjust = 1
        )
      
      rng_w <- range(df_ctx$context_score, na.rm = TRUE)
      
      if (all(is.finite(rng_w)) && diff(rng_w) > 0) {
        
        df_ctx <- df_ctx %>%
          mutate(
            edge_width = 0.45 + 2.0 *
              (context_score - rng_w[1]) /
              (rng_w[2] - rng_w[1])
          )
        
      } else {
        
        df_ctx <- df_ctx %>%
          mutate(edge_width = 0.9)
      }
      
      df_ctx <- df_ctx %>%
        mutate(
          edge_alpha = dplyr::case_when(
            edge_relation == "Grupo-like" ~ 0.90,
            edge_relation == "Shared FD/SO" ~ 0.70,
            edge_relation == "Splitter FD/SO" ~ 0.85,
            edge_relation == "Ambigua" ~ 0.35,
            TRUE ~ 0.60
          )
        )
      
      axis_nodes <- df_ctx %>%
        group_by(
          axis_id,
          axis_label,
          side_group,
          edge_relation
        ) %>%
        summarise(
          x = first(axis_x),
          y = median(y, na.rm = TRUE),
          n_features = n(),
          mean_context_score = mean(context_score, na.rm = TRUE),
          mean_MI_norm = mean(MI_norm, na.rm = TRUE),
          axis_edge_width = mean(edge_width, na.rm = TRUE),
          axis_edge_alpha = mean(edge_alpha, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          axis_label_plot = paste0(
            axis_label,
            "\n",
            "n = ",
            n_features
          )
        )
      
      group_nodes <- axis_nodes %>%
        group_by(side_group) %>%
        summarise(
          x = 0.92,
          y = median(y, na.rm = TRUE),
          n_axes = n(),
          .groups = "drop"
        ) %>%
        rename(group = side_group)
      
      axis_group_edges <- axis_nodes %>%
        left_join(
          group_nodes %>%
            rename(
              group_x = x,
              group_y = y
            ),
          by = c("side_group" = "group")
        ) %>%
        mutate(
          xend = group_x,
          yend = group_y
        )
      
      context_title <- dplyr::case_when(
        context_type == "unary" ~ paste0(
          "DAG explicativo de evidencia bidireccional: ",
          left_group
        ),
        
        TRUE ~ paste0(
          "DAG explicativo de evidencia bidireccional: ",
          left_group,
          " vs ",
          right_group
        )
      )
      
      context_subtitle <- dplyr::case_when(
        context_type == "unary" ~ paste0(
          "Biomarcador -> eje ",
          left_group,
          "-like -> grupo. Flechas = evidencia estadística, no causalidad."
        ),
        
        context_id == "NP_vs_FD" ~
          "Biomarcador -> eje bidireccional -> grupo. Incluye NP-like, FD-like y shared FD/SO como eje no-NP. No causal.",
        
        context_id == "NP_vs_SO" ~
          "Biomarcador -> eje bidireccional -> grupo. Incluye NP-like, SO-like y shared FD/SO como eje no-NP. No causal.",
        
        context_id == "FD_vs_SO" ~
          "Biomarcador -> eje bidireccional -> grupo. Incluye FD-like, SO-like y splitters FD/SO. No causal.",
        
        TRUE ~ "Biomarcador -> eje bidireccional -> grupo. No causal."
      )
      
      p11_ctx <- ggplot() +
        
        ## --------------------------------------------------------
      ## Aristas biomarcador -> eje.
      ## --------------------------------------------------------
      geom_segment(
        data = df_ctx,
        aes(
          x = feature_x,
          y = y,
          xend = axis_x - 0.055,
          yend = y,
          linewidth = edge_width,
          linetype = edge_relation,
          alpha = edge_alpha
        ),
        color = "grey30",
        arrow = grid::arrow(
          length = grid::unit(0.085, "inches"),
          type = "closed"
        )
      ) +
        
        ## --------------------------------------------------------
      ## Aristas eje -> grupo.
      ## --------------------------------------------------------
      geom_segment(
        data = axis_group_edges,
        aes(
          x = x + 0.075,
          y = y,
          xend = xend - 0.055,
          yend = yend,
          linewidth = axis_edge_width,
          linetype = edge_relation,
          alpha = axis_edge_alpha
        ),
        color = "grey20",
        arrow = grid::arrow(
          length = grid::unit(0.095, "inches"),
          type = "closed"
        )
      ) +
        
        ## --------------------------------------------------------
      ## Nodos eje/patrón.
      ## --------------------------------------------------------
      geom_label(
        data = axis_nodes,
        aes(
          x = x,
          y = y,
          label = axis_label_plot
        ),
        size = 3.35,
        fontface = "bold",
        fill = "grey93",
        color = "black",
        label.size = 0.25,
        lineheight = 0.92
      ) +
        
        ## --------------------------------------------------------
      ## Nodos grupo.
      ## --------------------------------------------------------
      geom_label(
        data = group_nodes,
        aes(
          x = x,
          y = y,
          label = group
        ),
        size = 5,
        fontface = "bold",
        fill = "grey88",
        color = "black",
        label.size = 0.30
      ) +
        
        ## --------------------------------------------------------
      ## Nodos biomarcador.
      ## --------------------------------------------------------
      geom_point(
        data = df_ctx,
        aes(
          x = feature_x,
          y = y,
          fill = view_code,
          size = MI_norm
        ),
        shape = 21,
        color = "black",
        alpha = 0.95
      ) +
        
        geom_text(
          data = df_ctx,
          aes(
            x = feature_label_x,
            y = y,
            label = node_label,
            hjust = feature_label_hjust
          ),
          size = 3.05
        ) +
        
        ## --------------------------------------------------------
      ## Encabezados de capas.
      ## --------------------------------------------------------
      annotate(
        "text",
        x = 0.18,
        y = max(df_ctx$y, na.rm = TRUE) + 1.05,
        label = "Biomarcadores",
        fontface = "bold",
        size = 4
      ) +
        annotate(
          "text",
          x = 0.58,
          y = max(df_ctx$y, na.rm = TRUE) + 1.05,
          label = "Eje bidireccional",
          fontface = "bold",
          size = 4
        ) +
        annotate(
          "text",
          x = 0.92,
          y = max(df_ctx$y, na.rm = TRUE) + 1.05,
          label = "Grupo",
          fontface = "bold",
          size = 4
        ) +
        
        scale_linewidth_identity(
          guide = "none"
        ) +
        
        scale_alpha_identity(
          guide = "none"
        ) +
        
        scale_linetype_manual(
          values = c(
            "Grupo-like" = "solid",
            "Shared FD/SO" = "dashed",
            "Splitter FD/SO" = "longdash",
            "Ambigua" = "dotted"
          ),
          drop = FALSE
        ) +
        
        scale_size_continuous(
          range = c(3.2, 7.5),
          limits = c(
            0,
            max(df_ctx$MI_norm, na.rm = TRUE)
          )
        ) +
        
        coord_cartesian(
          xlim = c(-0.10, 1.04),
          ylim = c(
            min(df_ctx$y, na.rm = TRUE) - 0.8,
            max(df_ctx$y, na.rm = TRUE) + 1.45
          ),
          clip = "off"
        ) +
        
        labs(
          title = context_title,
          subtitle = context_subtitle,
          x = NULL,
          y = NULL,
          fill = "Vista ómica",
          size = "MI norm",
          linetype = "Tipo de eje"
        ) +
        
        guides(
          linetype = guide_legend(
            override.aes = list(
              linewidth = 0.8,
              alpha = 1
            )
          )
        ) +
        
        theme_void() +
        theme(
          plot.title = element_text(
            face = "bold",
            size = 13
          ),
          plot.subtitle = element_text(
            size = 10
          ),
          plot.margin = margin(22, 90, 22, 135),
          legend.position = "right"
        )
      
      filename_id <- paste0(
        "11_feature_group_",
        context_id
      )
      
      .save_plot_publication(
        filename_base = file.path(
          PLOTS_DIR,
          filename_id
        ),
        plot = p11_ctx,
        width = ifelse(context_type == "unary", 11.5, 12.8),
        height = max(
          5.8,
          0.42 * nrow(df_ctx) + 2.3
        ),
        dpi = 600
      )
      
      add_manifest(
        paste0(
          filename_id,
          ".png / ",
          filename_id,
          ".pdf"
        ),
        paste0(
          "Familia 11 como DAG explicativo por contexto: ",
          context_id,
          ". Capas biomarcador -> eje bidireccional -> grupo; no representa causalidad."
        )
      )
      
      invisible(p11_ctx)
    }
    
    ## ------------------------------------------------------------
    ## Ejecutar los 6 contextos de la familia 11.
    ## ------------------------------------------------------------
    
    for (ii in seq_len(nrow(plot11_context_defs))) {
      
      plot11_make_context(
        context_id = plot11_context_defs$context_id[ii],
        context_type = plot11_context_defs$context_type[ii],
        left_group = plot11_context_defs$left_group[ii],
        right_group = plot11_context_defs$right_group[ii],
        patterns_keep = plot11_context_defs$patterns_keep[[ii]],
        contrast_estimate_col = plot11_context_defs$contrast_estimate_col[ii],
        contrast_prob_gt_col = plot11_context_defs$contrast_prob_gt_col[ii],
        contrast_prob_lt_col = plot11_context_defs$contrast_prob_lt_col[ii],
        contrast_js_col = plot11_context_defs$contrast_js_col[ii]
      )
    }
  }

  ## ============================================================
  ## PLOT 12. RED FEATURE-FEATURE HÍBRIDA MULTI-ÓMICA
  ## ------------------------------------------------------------
  ## Objetivo:
  ##   detectar módulos de features coordinadas entre muestras.
  ##
  ## Arista:
  ##   Spearman rho entre features, filtrado por |rho| >= network_rho_thr.
  ##
  ## Capa bidireccional:
  ##   cada arista se clasifica según compatibilidad de eje:
  ##   - same_axis
  ##   - compatible_axis
  ##   - ambiguous_axis
  ##   - opposite_or_mixed_axis
  ##
  ## Nodo:
  ##   color = vista ómica
  ##   tamaño = MI_norm
  ##
  ## Se excluyen:
  ##   E_no_evidence
  ## ============================================================
  
  if (!is.null(X_S)) {
    
    X_net <- as.data.frame(X_S, check.names = FALSE)
    
    required_ff_cols <- c(
      "feature_model",
      "feature_plot",
      "view_code",
      "evidence_grade",
      "bidirectional_pattern",
      "MI_norm",
      "max_JS_bits",
      "q_perm"
    )
    
    if (all(required_ff_cols %in% colnames(feature_signature))) {
      
      ff_features <- feature_signature %>%
        filter(
          evidence_grade %in% c(
            "A_strong",
            "B_moderate",
            "C_exploratory"
          ),
          bidirectional_pattern != "uncertain",
          feature_model %in% colnames(X_net)
        ) %>%
        mutate(
          biological_axis = dplyr::case_when(
            bidirectional_pattern == "SO_like_high" ~ "SO_like_axis",
            bidirectional_pattern == "FD_like_high" ~ "FD_like_axis",
            bidirectional_pattern == "NP_like_high" ~ "NP_like_axis",
            bidirectional_pattern %in% c(
              "SO_gt_FD_splitter",
              "FD_gt_SO_splitter"
            ) ~ "FD_SO_splitter_axis",
            bidirectional_pattern == "shared_FD_SO_high" ~ "shared_FD_SO_axis",
            evidence_grade == "D_ambiguous_dependent" ~ "ambiguous_dependent_axis",
            TRUE ~ "unclassified_axis"
          ),
          feature_label_network = make_short_network_label(as.character(feature_plot))
        )
      
      if (nrow(ff_features) >= 3) {
        
        feature_ids <- ff_features$feature_model
        
        X_cor <- X_net[, feature_ids, drop = FALSE]
        
        X_cor[] <- lapply(X_cor, function(z) {
          suppressWarnings(as.numeric(z))
        })
        
        cor_mat <- suppressWarnings(
          stats::cor(
            X_cor,
            method = "spearman",
            use = "pairwise.complete.obs"
          )
        )
        
        cor_mat[!is.finite(cor_mat)] <- 0
        diag(cor_mat) <- 1
        
        edge_idx <- which(
          upper.tri(cor_mat) & abs(cor_mat) >= network_rho_thr,
          arr.ind = TRUE
        )
        
        if (nrow(edge_idx) > 0) {
          
          ff_edges <- tibble(
            from = colnames(cor_mat)[edge_idx[, 1]],
            to = colnames(cor_mat)[edge_idx[, 2]],
            rho = cor_mat[edge_idx],
            abs_rho = abs(cor_mat[edge_idx])
          ) %>%
            left_join(
              ff_features %>%
                select(
                  feature_model,
                  from_label = feature_label_network,
                  from_view = view_code,
                  from_grade = evidence_grade,
                  from_pattern = bidirectional_pattern,
                  from_axis = biological_axis
                ),
              by = c("from" = "feature_model")
            ) %>%
            left_join(
              ff_features %>%
                select(
                  feature_model,
                  to_label = feature_label_network,
                  to_view = view_code,
                  to_grade = evidence_grade,
                  to_pattern = bidirectional_pattern,
                  to_axis = biological_axis
                ),
              by = c("to" = "feature_model")
            ) %>%
            mutate(
              rho_sign = ifelse(rho >= 0, "positive", "negative"),
              
              edge_axis_relation = dplyr::case_when(
                from_axis == "ambiguous_dependent_axis" |
                  to_axis == "ambiguous_dependent_axis" ~
                  "ambiguous_axis",
                
                from_axis == to_axis ~
                  "same_axis",
                
                from_axis == "shared_FD_SO_axis" &
                  to_axis %in% c(
                    "SO_like_axis",
                    "FD_like_axis",
                    "FD_SO_splitter_axis"
                  ) ~
                  "compatible_axis",
                
                to_axis == "shared_FD_SO_axis" &
                  from_axis %in% c(
                    "SO_like_axis",
                    "FD_like_axis",
                    "FD_SO_splitter_axis"
                  ) ~
                  "compatible_axis",
                
                from_axis == "FD_SO_splitter_axis" &
                  to_axis %in% c(
                    "SO_like_axis",
                    "FD_like_axis",
                    "shared_FD_SO_axis"
                  ) ~
                  "compatible_axis",
                
                to_axis == "FD_SO_splitter_axis" &
                  from_axis %in% c(
                    "SO_like_axis",
                    "FD_like_axis",
                    "shared_FD_SO_axis"
                  ) ~
                  "compatible_axis",
                
                TRUE ~
                  "opposite_or_mixed_axis"
              ),
              edge_axis_relation = factor(
                edge_axis_relation,
                levels = c(
                  "same_axis",
                  "compatible_axis",
                  "ambiguous_axis",
                  "opposite_or_mixed_axis"
                )
              ),
              
              edge_priority = dplyr::case_when(
                edge_axis_relation == "same_axis" ~ 3,
                edge_axis_relation == "compatible_axis" ~ 2,
                edge_axis_relation == "ambiguous_axis" ~ 1,
                TRUE ~ 0
              ),
              
              edge_score = abs_rho * edge_priority
            )
          
          if (!show_opposite_edges) {
            ff_edges <- ff_edges %>%
              filter(edge_axis_relation != "opposite_or_mixed_axis")
          }
          ff_edges <- ff_edges %>%
            arrange(
              desc(edge_score),
              desc(abs_rho)
            )
          
          if (nrow(ff_edges) > network_max_edges) {
            ff_edges <- ff_edges %>%
              slice_head(n = network_max_edges)
          }
          
          if (nrow(ff_edges) > 0) {
            
            node_ids <- sort(unique(c(ff_edges$from, ff_edges$to)))
            node_ids <- node_ids[node_ids %in% colnames(cor_mat)]
            
            cor_sub <- cor_mat[node_ids, node_ids, drop = FALSE]
            cor_sub <- as.matrix(cor_sub)
            
            rownames(cor_sub) <- node_ids
            colnames(cor_sub) <- node_ids
            
            ## ------------------------------------------------------------
            ## Distancia para layout:
            ## 1 - |rho| mantiene cerca features muy correlacionadas,
            ## tanto positivas como negativas.
            ## Importante: no usar pmax() directamente sobre la matriz,
            ## porque puede perder dimnames y romper cmdscale().
            ## ------------------------------------------------------------
            
            dist_mat <- 1 - abs(cor_sub)
            dist_mat[!is.finite(dist_mat)] <- 1
            dist_mat[dist_mat < 0] <- 0
            dist_mat <- (dist_mat + t(dist_mat)) / 2
            diag(dist_mat) <- 0
            
            rownames(dist_mat) <- node_ids
            colnames(dist_mat) <- node_ids
            
            dist_sub <- stats::as.dist(dist_mat)
            
            coords <- try(
              stats::cmdscale(dist_sub, k = 2),
              silent = TRUE
            )
            
            bad_coords <- inherits(coords, "try-error") ||
              is.null(coords) ||
              length(coords) == 0
            
            if (!bad_coords) {
              coords <- as.matrix(coords)
              bad_coords <- ncol(coords) < 2 ||
                nrow(coords) != length(node_ids) ||
                any(!is.finite(coords))
            }
            
            if (bad_coords) {
              
              theta <- seq(
                0,
                2 * pi,
                length.out = length(node_ids) + 1
              )[-1]
              
              coords <- cbind(
                cos(theta),
                sin(theta)
              )
            }
            
            rownames(coords) <- node_ids
            
            coords_tbl <- tibble::tibble(
              feature_model = node_ids,
              x = as.numeric(coords[node_ids, 1]),
              y = as.numeric(coords[node_ids, 2])
            )
            ff_nodes <- ff_features %>%
              filter(feature_model %in% node_ids) %>%
              left_join(coords_tbl, by = "feature_model") %>%
              mutate(
                node_id = feature_model,
                node_label = feature_label_network,
                
                phenotype_module = dplyr::recode(
                  biological_axis,
                  SO_like_axis = "SO-like",
                  FD_like_axis = "FD-like",
                  NP_like_axis = "NP-like",
                  FD_SO_splitter_axis = "FD/SO splitter",
                  shared_FD_SO_axis = "shared FD/SO",
                  ambiguous_dependent_axis = "ambiguous dependent",
                  unclassified_axis = "unclassified",
                  .default = biological_axis
                ),
                
                phenotype_module = factor(
                  phenotype_module,
                  levels = c(
                    "NP-like",
                    "FD-like",
                    "SO-like",
                    "FD/SO splitter",
                    "shared FD/SO",
                    "ambiguous dependent",
                    "unclassified"
                  )
                )
              ) %>%
              select(
                node_id,
                node_label,
                x,
                y,
                view_code,
                evidence_grade,
                bidirectional_pattern,
                biological_axis,
                phenotype_module,
                MI_norm,
                max_JS_bits,
                q_perm
              )
            
            ## ------------------------------------------------------------
            ## Zonas fenotípicas visibles:
            ## correlación define la red; direccionalidad asigna el eje.
            ## Las envolventes ayudan a ver los grupos de interés.
            ## ------------------------------------------------------------
            
            hull_tbl <- ff_nodes %>%
              filter(
                is.finite(x),
                is.finite(y),
                !is.na(phenotype_module)
              ) %>%
              group_by(phenotype_module) %>%
              filter(n() >= 3) %>%
              slice(chull(x, y)) %>%
              ungroup()
            
            phenotype_label_tbl <- ff_nodes %>%
              filter(
                is.finite(x),
                is.finite(y),
                !is.na(phenotype_module)
              ) %>%
              group_by(phenotype_module) %>%
              summarise(
                x = median(x, na.rm = TRUE),
                y = max(y, na.rm = TRUE) + 0.08,
                n_features = n(),
                .groups = "drop"
              ) %>%
              mutate(
                label = paste0(
                  phenotype_module,
                  "\n",
                  "n = ",
                  n_features
                )
              )
            node_degree_tbl <- ff_edges %>%
              select(from, to, abs_rho, edge_axis_relation, edge_score) %>%
              tidyr::pivot_longer(
                cols = c(from, to),
                names_to = "endpoint",
                values_to = "node_id"
              ) %>%
              group_by(node_id) %>%
              summarise(
                degree = n(),
                weighted_degree = sum(abs_rho, na.rm = TRUE),
                same_axis_degree = sum(edge_axis_relation == "same_axis", na.rm = TRUE),
                compatible_axis_degree = sum(edge_axis_relation == "compatible_axis", na.rm = TRUE),
                mean_edge_score = mean(edge_score, na.rm = TRUE),
                .groups = "drop"
              )
            
            ff_nodes <- ff_nodes %>%
              left_join(
                node_degree_tbl,
                by = "node_id"
              ) %>%
              mutate(
                degree = dplyr::coalesce(degree, 0L),
                weighted_degree = dplyr::coalesce(weighted_degree, 0),
                same_axis_degree = dplyr::coalesce(same_axis_degree, 0L),
                compatible_axis_degree = dplyr::coalesce(compatible_axis_degree, 0L),
                mean_edge_score = dplyr::coalesce(mean_edge_score, 0),
                
                label_show = dplyr::case_when(
                  degree >= 2 ~ TRUE,
                  MI_norm >= stats::quantile(MI_norm, 0.75, na.rm = TRUE) ~ TRUE,
                  evidence_grade == "A_strong" ~ TRUE,
                  TRUE ~ FALSE
                )
              )
            ff_edges_plot <- ff_edges %>%
              left_join(
                ff_nodes %>%
                  select(node_id, x, y) %>%
                  rename(
                    from_x = x,
                    from_y = y
                  ),
                by = c("from" = "node_id")
              ) %>%
              left_join(
                ff_nodes %>%
                  select(node_id, x, y) %>%
                  rename(
                    to_x = x,
                    to_y = y
                  ),
                by = c("to" = "node_id")
              ) %>%
              mutate(
                same_endpoint = is.finite(from_x) &
                  is.finite(from_y) &
                  is.finite(to_x) &
                  is.finite(to_y) &
                  abs(from_x - to_x) < 1e-10 &
                  abs(from_y - to_y) < 1e-10
              ) %>%
              filter(
                is.finite(from_x),
                is.finite(from_y),
                is.finite(to_x),
                is.finite(to_y),
                !same_endpoint
              )
            if (nrow(ff_edges_plot) == 0) {
              cat("\nNo quedan aristas válidas para graficar la red feature-feature tras filtrar endpoints idénticos.\n")
            } else {
              
              
            write.csv(
              ff_nodes,
              out_table("07_feature_feature_hybrid_network_nodes.csv"),
              row.names = FALSE
            )
            
            write.csv(
              ff_edges,
              out_table("08_feature_feature_hybrid_network_edges.csv"),
              row.names = FALSE
            )
            
            write.csv(
              node_degree_tbl,
              out_table("09b_feature_feature_node_degree_summary.csv"),
              row.names = FALSE
            )
            
            phenotype_module_summary <- ff_nodes %>%
              count(
                phenotype_module,
                view_code,
                evidence_grade,
                name = "n_features"
              ) %>%
              arrange(
                phenotype_module,
                view_code,
                evidence_grade
              )
            
            write.csv(
              phenotype_module_summary,
              out_table("09_feature_feature_phenotype_module_summary.csv"),
              row.names = FALSE
            )
            
            edge_module_summary <- ff_edges %>%
              mutate(
                from_module = dplyr::recode(
                  from_axis,
                  SO_like_axis = "SO-like",
                  FD_like_axis = "FD-like",
                  NP_like_axis = "NP-like",
                  FD_SO_splitter_axis = "FD/SO splitter",
                  shared_FD_SO_axis = "shared FD/SO",
                  ambiguous_dependent_axis = "ambiguous dependent",
                  unclassified_axis = "unclassified",
                  .default = from_axis
                ),
                to_module = dplyr::recode(
                  to_axis,
                  SO_like_axis = "SO-like",
                  FD_like_axis = "FD-like",
                  NP_like_axis = "NP-like",
                  FD_SO_splitter_axis = "FD/SO splitter",
                  shared_FD_SO_axis = "shared FD/SO",
                  ambiguous_dependent_axis = "ambiguous dependent",
                  unclassified_axis = "unclassified",
                  .default = to_axis
                ),
                module_pair = ifelse(
                  from_module <= to_module,
                  paste(from_module, to_module, sep = " -- "),
                  paste(to_module, from_module, sep = " -- ")
                )
              ) %>%
              group_by(
                module_pair,
                edge_axis_relation,
                rho_sign
              ) %>%
              summarise(
                n_edges = n(),
                mean_abs_rho = mean(abs_rho, na.rm = TRUE),
                median_abs_rho = median(abs_rho, na.rm = TRUE),
                max_abs_rho = max(abs_rho, na.rm = TRUE),
                mean_edge_score = mean(edge_score, na.rm = TRUE),
                .groups = "drop"
              ) %>%
              arrange(
                desc(n_edges),
                desc(mean_abs_rho),
                desc(mean_edge_score)
              )
            
            write.csv(
              edge_module_summary,
              out_table("09c_feature_feature_edge_module_summary.csv"),
              row.names = FALSE
            )
            
            
            p_ff_network <- ggplot() +
              
              ## ------------------------------------------------------------
            ## Fondo: grupos fenotípicos de interés
            ## ------------------------------------------------------------
            geom_polygon(
              data = hull_tbl,
              aes(
                x = x,
                y = y,
                group = phenotype_module
              ),
              fill = "grey92",
              color = "grey45",
              linewidth = 0.5,
              linetype = "dashed",
              alpha = 0.45,
              show.legend = FALSE
            ) +
              
              geom_label(
                data = phenotype_label_tbl,
                aes(
                  x = x,
                  y = y,
                  label = label
                ),
                size = 3.2,
                fontface = "bold",
                label.size = 0.25,
                fill = "white",
                alpha = 0.90,
                show.legend = FALSE
              ) +
              
              ## ------------------------------------------------------------
            ## Aristas: correlación entre features
            ## ------------------------------------------------------------
            geom_curve(
              data = ff_edges_plot,
              aes(
                x = from_x,
                y = from_y,
                xend = to_x,
                yend = to_y,
                linewidth = abs_rho,
                linetype = edge_axis_relation,
                color = rho_sign
              ),
              curvature = 0.12,
              alpha = 0.60
            ) +
              
              ## ------------------------------------------------------------
            ## Nodos: features
            ## ------------------------------------------------------------
            geom_point(
              data = ff_nodes,
              aes(
                x = x,
                y = y,
                fill = view_code,
                size = MI_norm
              ),
              shape = 21,
              color = "black",
              alpha = 0.95
            ) +
              
              geom_label(
                data = ff_nodes %>%
                  filter(label_show),
                aes(
                  x = x,
                  y = y,
                  label = node_label
                ),
                size = 2.8,
                vjust = -0.9,
                label.size = 0.15,
                fill = "white",
                alpha = 0.85,
                check_overlap = TRUE
              ) +
              
              scale_linewidth_continuous(
                range = c(0.3, 2.2)
              ) +
              
              scale_linetype_manual(
                values = c(
                  same_axis = "solid",
                  compatible_axis = "longdash",
                  ambiguous_axis = "dotted",
                  opposite_or_mixed_axis = "twodash"
                ),
                drop = FALSE
              ) +
              
              coord_equal(clip = "off") +
              
              labs(
                title = "Red de co-variación feature–feature por módulos fenotípicos",
                subtitle = paste0(
                  "Aristas: Spearman |rho| ≥ ",
                  network_rho_thr,
                  ". Se priorizan aristas fuertes y compatibles con el eje bidireccional. No representa causalidad."
                ),
                x = NULL,
                y = NULL,
                fill = "Vista",
                size = "MI norm",
                linewidth = "|Spearman rho|",
                linetype = "Relación direccional",
                color = "Signo rho"
              ) +
              
              theme_void() +
              theme(
                plot.margin = margin(30, 90, 30, 90),
                legend.position = "right"
              )
            .save_plot_publication(
              filename_base = file.path(
                PLOTS_DIR,
                "12_feature_feature_phenotype_module_network"
              ),
              plot = p_ff_network,
              width = 14,
              height = 11,
              dpi = 600
            )
            
            add_manifest(
              "12_feature_feature_phenotype_module_network.png / 12_feature_feature_phenotype_module_network.pdf",
              paste0(
                "Red feature–feature por módulos fenotípicos: aristas por Spearman |rho| >= ",
                network_rho_thr,
                "; las zonas visibles corresponden al eje fenotípico asignado por direccionalidad bidireccional."
              )
            )
            ## ------------------------------------------------------------
            ## PLOT 12B. RED FEATURE-FEATURE SEPARADA POR GRUPO
            ## ------------------------------------------------------------
            
            if (isTRUE(network_split_by_group)) {
              
              for (group_interest in y_levels) {
                
                modules_g <- phenotype_modules_for_group(group_interest)
                
                ff_nodes_g <- ff_nodes %>%
                  filter(
                    as.character(phenotype_module) %in% modules_g
                  )
                
                ff_edges_plot_g <- ff_edges_plot %>%
                  filter(
                    from %in% ff_nodes_g$node_id,
                    to %in% ff_nodes_g$node_id,
                    abs_rho >= network_group_rho_thr
                  ) %>%
                  arrange(
                    factor(
                      edge_axis_relation,
                      levels = c(
                        "same_axis",
                        "compatible_axis",
                        "ambiguous_axis",
                        "opposite_or_mixed_axis"
                      )
                    ),
                    desc(abs_rho)
                  ) %>%
                  slice_head(n = network_group_max_edges)
                
                if (nrow(ff_edges_plot_g) == 0) {
                  next
                }
                
                node_ids_g <- sort(unique(c(ff_edges_plot_g$from, ff_edges_plot_g$to)))
                
                ff_nodes_g <- ff_nodes_g %>%
                  filter(node_id %in% node_ids_g) %>%
                  mutate(
                    node_label = make_short_network_label(node_label),
                    label_show = dplyr::case_when(
                      degree >= 2 ~ TRUE,
                      MI_norm >= stats::quantile(MI_norm, 0.75, na.rm = TRUE) ~ TRUE,
                      evidence_grade == "A_strong" ~ TRUE,
                      TRUE ~ FALSE
                    )
                  )
                
                if (nrow(ff_nodes_g) < 2) {
                  next
                }
                
                hull_tbl_g <- ff_nodes_g %>%
                  filter(
                    is.finite(x),
                    is.finite(y),
                    !is.na(phenotype_module)
                  ) %>%
                  group_by(phenotype_module) %>%
                  filter(n() >= 3) %>%
                  slice(chull(x, y)) %>%
                  ungroup()
                
                phenotype_label_tbl_g <- ff_nodes_g %>%
                  filter(
                    is.finite(x),
                    is.finite(y),
                    !is.na(phenotype_module)
                  ) %>%
                  group_by(phenotype_module) %>%
                  summarise(
                    x = median(x, na.rm = TRUE),
                    y = max(y, na.rm = TRUE) + 0.08,
                    n_features = n(),
                    .groups = "drop"
                  ) %>%
                  mutate(
                    label = paste0(
                      phenotype_module,
                      "\n",
                      "n = ",
                      n_features
                    )
                  )
                
                write.csv(
                  ff_nodes_g,
                  out_table(paste0("07_feature_feature_hybrid_network_nodes_", group_interest, ".csv")),
                  row.names = FALSE
                )
                
                write.csv(
                  ff_edges_plot_g,
                  out_table(paste0("08_feature_feature_hybrid_network_edges_", group_interest, ".csv")),
                  row.names = FALSE
                )
                
                p_ff_network_g <- ggplot() +
                  geom_polygon(
                    data = hull_tbl_g,
                    aes(
                      x = x,
                      y = y,
                      group = phenotype_module
                    ),
                    fill = "grey92",
                    color = "grey45",
                    linewidth = 0.5,
                    linetype = "dashed",
                    alpha = 0.45,
                    show.legend = FALSE
                  ) +
                  geom_label(
                    data = phenotype_label_tbl_g,
                    aes(
                      x = x,
                      y = y,
                      label = label
                    ),
                    size = 3.2,
                    fontface = "bold",
                    label.size = 0.25,
                    fill = "white",
                    alpha = 0.90,
                    show.legend = FALSE
                  ) +
                  geom_curve(
                    data = ff_edges_plot_g,
                    aes(
                      x = from_x,
                      y = from_y,
                      xend = to_x,
                      yend = to_y,
                      linewidth = abs_rho,
                      linetype = edge_axis_relation,
                      color = rho_sign
                    ),
                    curvature = 0.12,
                    alpha = 0.60
                  ) +
                  geom_point(
                    data = ff_nodes_g,
                    aes(
                      x = x,
                      y = y,
                      fill = view_code,
                      size = MI_norm
                    ),
                    shape = 21,
                    color = "black",
                    alpha = 0.95
                  ) +
                  geom_label(
                    data = ff_nodes_g %>%
                      filter(label_show),
                    aes(
                      x = x,
                      y = y,
                      label = node_label
                    ),
                    size = 2.8,
                    vjust = -0.9,
                    label.size = 0.15,
                    fill = "white",
                    alpha = 0.85,
                    check_overlap = TRUE
                  ) +
                  scale_linewidth_continuous(
                    range = c(0.3, 2.2)
                  ) +
                  scale_linetype_manual(
                    values = c(
                      same_axis = "solid",
                      compatible_axis = "longdash",
                      ambiguous_axis = "dotted",
                      opposite_or_mixed_axis = "twodash"
                    ),
                    drop = FALSE
                  ) +
                  coord_equal(clip = "off") +
                  labs(
                    title = paste0("Red de co-variación feature–feature: eje ", group_interest),
                    subtitle = paste0(
                      "Aristas: Spearman |rho| ≥ ",
                      network_group_rho_thr,
                      ". Se muestran conexiones fuertes dentro de módulos relevantes para ",
                      group_interest,
                      ". No causal."
                    ),
                    x = NULL,
                    y = NULL,
                    fill = "Vista",
                    size = "MI norm",
                    linewidth = "|Spearman rho|",
                    linetype = "Relación direccional",
                    color = "Signo rho"
                  ) +
                  theme_void() +
                  theme(
                    plot.margin = margin(30, 90, 30, 90),
                    legend.position = "right"
                  )
                
                .save_plot_publication(
                  filename_base = file.path(
                    PLOTS_DIR,
                    paste0("12_feature_feature_phenotype_module_network_", group_interest)
                  ),
                  plot = p_ff_network_g,
                  width = 13,
                  height = 10,
                  dpi = 600
                )
                
                add_manifest(
                  paste0(
                    "12_feature_feature_phenotype_module_network_",
                    group_interest,
                    ".png / 12_feature_feature_phenotype_module_network_",
                    group_interest,
                    ".pdf"
                  ),
                  paste0(
                    "Red feature–feature separada para ",
                    group_interest,
                    ": filtra módulos fenotípicos relevantes y reduce densidad de aristas."
                  )
                )
              }
            }
            
            }
            
          }
        }
      }
    }
  }
  ## ============================================================
  ## Guardar manifest
  ## ============================================================
  
  write.csv(
    plot_manifest,
    file.path(PLOTS_DIR, "00_plot_manifest.csv"),
    row.names = FALSE
  )
  
  cat("\n============================================================\n")
  cat("PLOTS INTERPRETABLES GUARDADOS EN:\n")
  cat(PLOTS_DIR, "\n")
  cat("============================================================\n")
  print(plot_manifest, n = Inf)
  
  invisible(
    list(
      feature_signature_plot = feature_signature,
      plot_manifest = plot_manifest,
      plots_dir = PLOTS_DIR
    )
  )
}

summarize_feature_modules_from_network <- function(outdir = "framework_bidireccional_firma",
                                                   min_module_size = 2) {
  
  tables_dir <- file.path(outdir, "tables")
  
  nodes_file <- file.path(
    tables_dir,
    "07_feature_feature_hybrid_network_nodes.csv"
  )
  
  edges_file <- file.path(
    tables_dir,
    "08_feature_feature_hybrid_network_edges.csv"
  )
  
  if (!file.exists(nodes_file)) {
    stop("No existe el archivo de nodos: ", nodes_file)
  }
  
  if (!file.exists(edges_file)) {
    stop("No existe el archivo de aristas: ", edges_file)
  }
  
  nodes <- read.csv(nodes_file, check.names = FALSE)
  edges <- read.csv(edges_file, check.names = FALSE)
  
  if (!all(c("node_id", "node_label") %in% colnames(nodes))) {
    stop("La tabla de nodos debe tener node_id y node_label.")
  }
  
  if (!all(c("from", "to") %in% colnames(edges))) {
    stop("La tabla de aristas debe tener from y to.")
  }
  
  node_ids <- unique(nodes$node_id)
  
  ## ------------------------------------------------------------
  ## Componentes conectados sin depender de igraph.
  ## Cada componente conectado = módulo empírico.
  ## ------------------------------------------------------------
  
  adj <- setNames(vector("list", length(node_ids)), node_ids)
  
  for (ii in seq_len(nrow(edges))) {
    a <- as.character(edges$from[ii])
    b <- as.character(edges$to[ii])
    
    if (a %in% node_ids && b %in% node_ids) {
      adj[[a]] <- unique(c(adj[[a]], b))
      adj[[b]] <- unique(c(adj[[b]], a))
    }
  }
  
  visited <- setNames(rep(FALSE, length(node_ids)), node_ids)
  module_list <- list()
  module_counter <- 0L
  
  for (id in node_ids) {
    
    if (visited[[id]]) next
    
    module_counter <- module_counter + 1L
    
    queue <- id
    component <- character()
    visited[[id]] <- TRUE
    
    while (length(queue) > 0) {
      
      current <- queue[1]
      queue <- queue[-1]
      component <- c(component, current)
      
      neigh <- adj[[current]]
      neigh <- neigh[!visited[neigh]]
      
      if (length(neigh) > 0) {
        visited[neigh] <- TRUE
        queue <- c(queue, neigh)
      }
    }
    
    module_list[[module_counter]] <- component
  }
  
  module_members_raw <- bind_rows(
    lapply(seq_along(module_list), function(i) {
      tibble(
        module_id = sprintf("M%02d", i),
        node_id = module_list[[i]]
      )
    })
  )
  
  module_members <- module_members_raw %>%
    left_join(nodes, by = "node_id") %>%
    group_by(module_id) %>%
    mutate(module_size = n()) %>%
    ungroup() %>%
    filter(module_size >= min_module_size)
  
  ## ------------------------------------------------------------
  ## Aristas anotadas por módulo.
  ## ------------------------------------------------------------
  
  edge_modules <- edges %>%
    left_join(
      module_members %>%
        select(module_id, from = node_id),
      by = "from"
    ) %>%
    left_join(
      module_members %>%
        select(module_id_to = module_id, to = node_id),
      by = "to"
    ) %>%
    filter(
      !is.na(module_id),
      !is.na(module_id_to),
      module_id == module_id_to
    )
  
  ## ------------------------------------------------------------
  ## Resumen de módulos.
  ## ------------------------------------------------------------
  
  mode_chr <- function(x) {
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) return(NA_character_)
    names(sort(table(x), decreasing = TRUE))[1]
  }
  
  collapse_unique <- function(x) {
    x <- unique(x[!is.na(x) & x != ""])
    if (length(x) == 0) return(NA_character_)
    paste(x, collapse = " | ")
  }
  
  module_top_features <- module_members %>%
    group_by(module_id) %>%
    summarise(
      top_features_by_MI = {
        ord <- order(MI_norm, decreasing = TRUE, na.last = TRUE)
        labs <- node_label[ord]
        labs <- labs[!is.na(labs) & labs != ""]
        paste(head(labs, 6), collapse = " | ")
      },
      .groups = "drop"
    )
  module_top_features <- module_members %>%
    group_by(module_id) %>%
    summarise(
      top_features_by_MI = {
        ord <- order(MI_norm, decreasing = TRUE, na.last = TRUE)
        labs <- node_label[ord]
        labs <- labs[!is.na(labs) & labs != ""]
        paste(head(labs, 6), collapse = " | ")
      },
      .groups = "drop"
    )
  
  module_summary <- module_members %>%
    group_by(module_id) %>%
    summarise(
      n_features = n(),
      views_present = collapse_unique(view_code),
      dominant_view = mode_chr(view_code),
      phenotype_modules_present = collapse_unique(phenotype_module),
      dominant_phenotype_module = mode_chr(phenotype_module),
      dominant_evidence_grade = mode_chr(evidence_grade),
      mean_MI_norm = mean(MI_norm, na.rm = TRUE),
      median_MI_norm = median(MI_norm, na.rm = TRUE),
      min_q_perm = min(q_perm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      module_top_features,
      by = "module_id"
    ) %>%
    left_join(
      edge_modules %>%
        group_by(module_id) %>%
        summarise(
          n_edges = n(),
          mean_abs_rho = mean(abs_rho, na.rm = TRUE),
          median_abs_rho = median(abs_rho, na.rm = TRUE),
          max_abs_rho = max(abs_rho, na.rm = TRUE),
          edge_axis_relations = collapse_unique(edge_axis_relation),
          .groups = "drop"
        ),
      by = "module_id"
    ) %>%
    mutate(
      n_edges = ifelse(is.na(n_edges), 0L, n_edges),
      suggested_neutral_name = paste0(
        module_id,
        "_",
        gsub("[^A-Za-z0-9]+", "_", dominant_phenotype_module),
        "_module"
      )
    ) %>%
    arrange(
      desc(n_features),
      desc(mean_MI_norm),
      desc(median_abs_rho)
    )
  
  ## ------------------------------------------------------------
  ## Tabla de anotación manual por módulo.
  ## No se fuerza adiposidad/inflamación/metabolismo.
  ## ------------------------------------------------------------
  
  module_annotation_template <- module_summary %>%
    transmute(
      module_id,
      suggested_neutral_name,
      n_features,
      views_present,
      dominant_view,
      dominant_phenotype_module,
      phenotype_modules_present,
      dominant_evidence_grade,
      mean_MI_norm,
      median_abs_rho,
      top_features_by_MI,
      
      manual_biological_annotation = "",
      annotation_basis = "",
      annotation_confidence = "",
      candidate_pathways = "",
      notes = ""
    )
  
  write.csv(
    module_members,
    file.path(tables_dir, "10_feature_feature_neutral_module_members.csv"),
    row.names = FALSE
  )
  
  write.csv(
    module_summary,
    file.path(tables_dir, "11_feature_feature_neutral_module_summary.csv"),
    row.names = FALSE
  )
  
  write.csv(
    module_annotation_template,
    file.path(tables_dir, "12_feature_feature_module_annotation_template.csv"),
    row.names = FALSE
  )
  
  list(
    module_members = module_members,
    module_summary = module_summary,
    module_annotation_template = module_annotation_template
  )
}


## ============================================================
## RANKING DIFERENCIAL DUAL:
## OBSERVADO/SIN INTEGRAR vs RECONSTRUIDO/INTEGRADO
## ------------------------------------------------------------
## Objetivo:
##   1) Extraer todas las features que entran en MOFA.
##   2) Crear dos matrices comparables:
##        A) observed     = señal original usada por MOFA
##        B) reconstructed = señal integrada Xhat = Z %*% t(W)
##   3) Hacer el mismo ranking diferencial en ambos flujos.
##   4) Marcar biomarcadores finales.
##   5) Generar plots de auditoría visual.
##
## Metabolitos:
##   - se rankean estadísticamente
##   - NO se anotan biológicamente si no hay HMDB/KEGG/ChEBI
## ============================================================

.standardize_factor_global <- function(x) {
  x <- as.character(x)
  ifelse(
    grepl("^Factor", x),
    x,
    paste0("Factor", gsub("[^0-9]", "", x))
  )
}

.standardize_view_global <- function(x) {
  
  x <- tolower(as.character(x))
  
  dplyr::case_when(
    grepl("trans", x) ~ "tx",
    grepl("prot",  x) ~ "pr",
    grepl("metab", x) ~ "me",
    grepl("clin",  x) ~ "cl",
    x %in% c("tx", "pr", "me", "cl") ~ x,
    TRUE ~ x
  )
}

.get_sample_ids_safe <- function(df, split_name) {
  
  sample_id <- rownames(df)
  
  if (
    is.null(sample_id) ||
    anyDuplicated(sample_id) ||
    all(sample_id == as.character(seq_len(nrow(df))))
  ) {
    sample_id <- paste0(split_name, "_", seq_len(nrow(df)))
  }
  
  sample_id
}

.resolve_path_from_diagnostic_global <- function(path_raw, diagnostic_rds) {
  
  path_raw <- path.expand(path_raw)
  
  if (file.exists(path_raw)) {
    return(normalizePath(path_raw, mustWork = TRUE))
  }
  
  repo_guess <- normalizePath(
    file.path(dirname(diagnostic_rds), "..", "..", "..", ".."),
    mustWork = FALSE
  )
  
  path_alt <- file.path(repo_guess, sub("^\\./", "", path_raw))
  
  if (file.exists(path_alt)) {
    return(normalizePath(path_alt, mustWork = TRUE))
  }
  
  stop(
    "No encuentro el archivo MOFA2:\n",
    path_raw,
    "\nRuta alternativa probada:\n",
    path_alt
  )
}

.load_mofa_model_robust_global <- function(path) {
  
  path <- path.expand(path)
  
  if (!file.exists(path)) {
    stop("No existe el archivo del modelo MOFA2:\n", path)
  }
  
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c("hdf5", "h5")) {
    return(MOFA2::load_model(path))
  }
  
  obj <- tryCatch(
    readRDS(path),
    error = function(e) NULL
  )
  
  if (!is.null(obj)) {
    return(obj)
  }
  
  MOFA2::load_model(path)
}

.load_diagnostic_context_for_dual_ranking <- function(scenario_dir,
                                                      diag_target = "diagnostico_bayes_features") {
  
  scenario_dir <- normalizePath(path.expand(scenario_dir), mustWork = TRUE)
  
  diagnostic_rds <- file.path(
    scenario_dir,
    diag_target,
    "rds",
    "diagnostic_objects.rds"
  )
  
  diagnostic_rds <- normalizePath(path.expand(diagnostic_rds), mustWork = TRUE)
  
  diagnostic_objects <- readRDS(diagnostic_rds)
  
  if (is.null(diagnostic_objects$df_train)) {
    stop("No encuentro diagnostic_objects$df_train")
  }
  
  if (is.null(diagnostic_objects$df_test)) {
    stop("No encuentro diagnostic_objects$df_test")
  }
  
  if (is.null(diagnostic_objects$df_weights)) {
    stop("No encuentro diagnostic_objects$df_weights")
  }
  
  if (is.null(diagnostic_objects$diagnostic_config$mofa_model_file)) {
    stop("No encuentro diagnostic_objects$diagnostic_config$mofa_model_file")
  }
  
  mofa_model_file <- .resolve_path_from_diagnostic_global(
    diagnostic_objects$diagnostic_config$mofa_model_file,
    diagnostic_rds
  )
  
  mofa_model <- .load_mofa_model_robust_global(mofa_model_file)
  
  df_train <- diagnostic_objects$df_train
  df_test  <- diagnostic_objects$df_test
  
  factor_cols <- grep("^Factor[0-9]+$", names(df_train), value = TRUE)
  factor_cols <- factor_cols[order(as.integer(gsub("^Factor", "", factor_cols)))]
  
  final_biomarkers <- setdiff(
    names(df_train),
    c("y", factor_cols)
  )
  
  y_tbl <- dplyr::bind_rows(
    tibble::tibble(
      sample_id = .get_sample_ids_safe(df_train, "train"),
      split = "train",
      y = as.character(df_train$y)
    ),
    tibble::tibble(
      sample_id = .get_sample_ids_safe(df_test, "test"),
      split = "test",
      y = as.character(df_test$y)
    )
  )
  
  list(
    scenario_dir = scenario_dir,
    diag_target = diag_target,
    diagnostic_rds = diagnostic_rds,
    diagnostic_objects = diagnostic_objects,
    df_train = df_train,
    df_test = df_test,
    factor_cols = factor_cols,
    final_biomarkers = final_biomarkers,
    y_tbl = y_tbl,
    mofa_model = mofa_model,
    mofa_model_file = mofa_model_file
  )
}

.get_mofa_weights_all_features <- function(mofa_model,
                                           factor_cols) {
  
  MOFA2::get_weights(
    mofa_model,
    views = "all",
    factors = "all",
    as.data.frame = TRUE
  ) %>%
    dplyr::mutate(
      view_code = .standardize_view_global(view),
      mofa_feature = as.character(feature),
      factor = .standardize_factor_global(factor),
      weight = as.numeric(value)
    ) %>%
    dplyr::filter(factor %in% factor_cols)
}

.make_feature_map_all_mofa <- function(diagnostic_objects,
                                       weights_long,
                                       final_biomarkers) {
  
  mofa_features <- weights_long %>%
    dplyr::distinct(
      view_code,
      mofa_feature
    )
  
  map_diag <- diagnostic_objects$df_weights %>%
    dplyr::transmute(
      view_code = .standardize_view_global(view_code),
      mofa_feature = as.character(feature),
      feature_model = as.character(feature_model)
    ) %>%
    dplyr::distinct()
  
  mofa_features %>%
    dplyr::left_join(
      map_diag,
      by = c("view_code", "mofa_feature")
    ) %>%
    dplyr::mutate(
      feature_model = dplyr::if_else(
        is.na(feature_model) | feature_model == "",
        paste0(view_code, "_", mofa_feature),
        feature_model
      ),
      feature_model = make.unique(feature_model),
      feature_label = sub("^[^_]+_", "", feature_model),
      final_biomarker = feature_model %in% final_biomarkers,
      annotation_eligible = view_code %in% c("tx", "pr", "cl"),
      annotation_route = dplyr::case_when(
        view_code == "tx" ~ "gene_level_GO_Reactome_MSigDB",
        view_code == "pr" ~ "protein_level_UniProt_GO_Reactome",
        view_code == "cl" ~ "clinical_ontology_EFO_HPO_MeSH_OLS",
        view_code == "me" ~ "not_annotated_metabolite_id_missing",
        TRUE ~ "unknown"
      )
    ) %>%
    dplyr::distinct(
      view_code,
      mofa_feature,
      .keep_all = TRUE
    )
}

.flatten_mofa_data_matrix <- function(mat,
                                      view_name,
                                      mofa_group_name = NA_character_) {
  
  mat <- as.matrix(mat)
  
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    return(NULL)
  }
  
  as.data.frame(mat, check.names = FALSE) %>%
    tibble::rownames_to_column("mofa_feature") %>%
    tidyr::pivot_longer(
      cols = -mofa_feature,
      names_to = "sample_id",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      view = as.character(view_name),
      view_code = .standardize_view_global(view),
      mofa_group = as.character(mofa_group_name),
      mofa_feature = as.character(mofa_feature),
      sample_id = as.character(sample_id),
      value = as.numeric(value)
    ) %>%
    dplyr::select(
      sample_id,
      mofa_group,
      view,
      view_code,
      mofa_feature,
      value
    )
}

.get_mofa_observed_data_long <- function(mofa_model) {
  
  dat_df <- tryCatch(
    MOFA2::get_data(
      mofa_model,
      views = "all",
      as.data.frame = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.data.frame(dat_df)) {
    
    cn <- colnames(dat_df)
    
    sample_col <- intersect(
      c("sample", "samples", "sample_id", "Sample", "SampleID"),
      cn
    )[1]
    
    feature_col <- intersect(
      c("feature", "features", "mofa_feature", "Feature"),
      cn
    )[1]
    
    view_col <- intersect(
      c("view", "views", "View"),
      cn
    )[1]
    
    value_col <- intersect(
      c("value", "values", "data", "x"),
      cn
    )[1]
    
    if (
      !is.na(sample_col) &&
      !is.na(feature_col) &&
      !is.na(view_col) &&
      !is.na(value_col)
    ) {
      
      return(
        tibble::tibble(
          sample_id = as.character(dat_df[[sample_col]]),
          view = as.character(dat_df[[view_col]]),
          view_code = .standardize_view_global(dat_df[[view_col]]),
          mofa_feature = as.character(dat_df[[feature_col]]),
          value = as.numeric(dat_df[[value_col]])
        ) %>%
          dplyr::mutate(mofa_group = NA_character_) %>%
          dplyr::select(
            sample_id,
            mofa_group,
            view,
            view_code,
            mofa_feature,
            value
          )
      )
    }
  }
  
  dat_list <- tryCatch(
    MOFA2::get_data(
      mofa_model,
      views = "all",
      as.data.frame = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(dat_list) || !is.list(dat_list)) {
    stop("No pude extraer datos observados desde MOFA2::get_data().")
  }
  
  out <- list()
  counter <- 0L
  
  for (lvl1 in names(dat_list)) {
    
    obj1 <- dat_list[[lvl1]]
    
    if (is.matrix(obj1) || is.data.frame(obj1)) {
      
      counter <- counter + 1L
      
      out[[counter]] <- .flatten_mofa_data_matrix(
        mat = obj1,
        view_name = lvl1,
        mofa_group_name = NA_character_
      )
      
    } else if (is.list(obj1)) {
      
      for (lvl2 in names(obj1)) {
        
        obj2 <- obj1[[lvl2]]
        
        if (!(is.matrix(obj2) || is.data.frame(obj2))) {
          next
        }
        
        lvl1_view <- .standardize_view_global(lvl1)
        lvl2_view <- .standardize_view_global(lvl2)
        
        if (lvl1_view %in% c("tx", "pr", "me", "cl")) {
          view_name <- lvl1
          group_name <- lvl2
        } else if (lvl2_view %in% c("tx", "pr", "me", "cl")) {
          view_name <- lvl2
          group_name <- lvl1
        } else {
          view_name <- lvl2
          group_name <- lvl1
        }
        
        counter <- counter + 1L
        
        out[[counter]] <- .flatten_mofa_data_matrix(
          mat = obj2,
          view_name = view_name,
          mofa_group_name = group_name
        )
      }
    }
  }
  
  out <- Filter(Negate(is.null), out)
  
  if (length(out) == 0) {
    stop("MOFA2::get_data() devolvió una estructura no convertible a tabla larga.")
  }
  
  dplyr::bind_rows(out)
}

.join_y_to_observed_mofa <- function(mofa_data_long,
                                     y_tbl) {
  
  joined <- mofa_data_long %>%
    dplyr::left_join(
      y_tbl,
      by = "sample_id"
    )
  
  if (!any(is.na(joined$y))) {
    return(joined)
  }
  
  mofa_samples <- unique(mofa_data_long$sample_id)
  
  y_train <- y_tbl %>%
    dplyr::filter(split == "train")
  
  if (length(mofa_samples) == nrow(y_train)) {
    
    warning(
      "No hubo match perfecto sample_id entre MOFA y diagnostic_objects. ",
      "Asigno y usando orden TRAIN."
    )
    
    fallback_y <- tibble::tibble(
      sample_id = mofa_samples,
      split = "train",
      y = y_train$y
    )
    
  } else if (length(mofa_samples) == nrow(y_tbl)) {
    
    warning(
      "No hubo match perfecto sample_id entre MOFA y diagnostic_objects. ",
      "Asigno y usando orden TRAIN+TEST."
    )
    
    fallback_y <- tibble::tibble(
      sample_id = mofa_samples,
      split = y_tbl$split,
      y = y_tbl$y
    )
    
  } else {
    
    stop(
      "No pude asignar grupos a muestras MOFA.\n",
      "Muestras MOFA: ", length(mofa_samples), "\n",
      "Filas y_tbl: ", nrow(y_tbl), "\n"
    )
  }
  
  mofa_data_long %>%
    dplyr::left_join(
      fallback_y,
      by = "sample_id"
    )
}

.make_observed_mofa_wide <- function(ctx,
                                     feature_map) {
  
  obs_long <- .get_mofa_observed_data_long(ctx$mofa_model)
  
  obs_long <- .join_y_to_observed_mofa(
    mofa_data_long = obs_long,
    y_tbl = ctx$y_tbl
  )
  
  obs_long %>%
    dplyr::left_join(
      feature_map %>%
        dplyr::select(
          view_code,
          mofa_feature,
          feature_model
        ),
      by = c("view_code", "mofa_feature")
    ) %>%
    dplyr::filter(!is.na(feature_model)) %>%
    dplyr::select(
      sample_id,
      split,
      y,
      feature_model,
      value
    ) %>%
    tidyr::pivot_wider(
      names_from = feature_model,
      values_from = value,
      values_fn = mean
    ) %>%
    dplyr::arrange(
      factor(split, levels = c("train", "test")),
      sample_id
    )
}

.make_reconstructed_mofa_wide <- function(ctx,
                                          feature_map,
                                          weights_long) {
  
  df_z <- dplyr::bind_rows(
    data.frame(
      sample_id = .get_sample_ids_safe(ctx$df_train, "train"),
      split = "train",
      y = ctx$df_train$y,
      ctx$df_train[, ctx$factor_cols, drop = FALSE],
      check.names = FALSE
    ),
    data.frame(
      sample_id = .get_sample_ids_safe(ctx$df_test, "test"),
      split = "test",
      y = ctx$df_test$y,
      ctx$df_test[, ctx$factor_cols, drop = FALSE],
      check.names = FALSE
    )
  )
  
  Z <- as.matrix(df_z[, ctx$factor_cols, drop = FALSE])
  
  reconstruct_view <- function(vc) {
    
    map_v <- feature_map %>%
      dplyr::filter(view_code == vc)
    
    if (nrow(map_v) == 0) {
      return(NULL)
    }
    
    W_v <- weights_long %>%
      dplyr::filter(
        view_code == vc,
        mofa_feature %in% map_v$mofa_feature
      ) %>%
      dplyr::select(
        mofa_feature,
        factor,
        weight
      ) %>%
      tidyr::pivot_wider(
        names_from = factor,
        values_from = weight,
        values_fill = 0
      )
    
    missing_weights <- setdiff(map_v$mofa_feature, W_v$mofa_feature)
    
    if (length(missing_weights) > 0) {
      warning(
        "Faltan pesos MOFA para vista ", vc, ": ",
        paste(missing_weights, collapse = ", ")
      )
    }
    
    map_v <- map_v %>%
      dplyr::filter(mofa_feature %in% W_v$mofa_feature)
    
    W_v <- W_v %>%
      dplyr::arrange(match(mofa_feature, map_v$mofa_feature))
    
    W_mat <- W_v %>%
      dplyr::select(dplyr::all_of(ctx$factor_cols)) %>%
      as.matrix()
    
    Xhat <- Z %*% t(W_mat)
    colnames(Xhat) <- map_v$feature_model
    
    Xhat
  }
  
  Xhat_list <- lapply(
    unique(feature_map$view_code),
    reconstruct_view
  )
  
  Xhat_list <- Filter(Negate(is.null), Xhat_list)
  
  if (length(Xhat_list) == 0) {
    stop("No se pudo reconstruir ninguna vista desde pesos MOFA.")
  }
  
  Xhat_mat <- do.call(cbind, Xhat_list)
  
  data.frame(
    sample_id = df_z$sample_id,
    split = df_z$split,
    y = df_z$y,
    Xhat_mat,
    check.names = FALSE
  )
}

.align_observed_and_reconstructed_wide <- function(df_observed,
                                                   df_reconstructed) {
  
  sample_common <- intersect(
    df_observed$sample_id,
    df_reconstructed$sample_id
  )
  
  if (length(sample_common) < 6) {
    stop(
      "Muy pocas muestras comunes entre observed y reconstructed: ",
      length(sample_common)
    )
  }
  
  feat_obs <- setdiff(
    colnames(df_observed),
    c("sample_id", "split", "y")
  )
  
  feat_rec <- setdiff(
    colnames(df_reconstructed),
    c("sample_id", "split", "y")
  )
  
  feature_common <- intersect(feat_obs, feat_rec)
  
  if (length(feature_common) < 3) {
    stop(
      "Muy pocas features comunes entre observed y reconstructed: ",
      length(feature_common)
    )
  }
  
  df_observed2 <- df_observed %>%
    dplyr::filter(sample_id %in% sample_common) %>%
    dplyr::arrange(match(sample_id, sample_common)) %>%
    dplyr::select(
      sample_id,
      split,
      y,
      dplyr::all_of(feature_common)
    )
  
  df_reconstructed2 <- df_reconstructed %>%
    dplyr::filter(sample_id %in% sample_common) %>%
    dplyr::arrange(match(sample_id, sample_common)) %>%
    dplyr::select(
      sample_id,
      split,
      y,
      dplyr::all_of(feature_common)
    )
  
  list(
    observed = df_observed2,
    reconstructed = df_reconstructed2,
    sample_common = sample_common,
    feature_common = feature_common
  )
}

make_mofa_dual_feature_ranking_input <- function(scenario_dir,
                                                 diag_target = "diagnostico_bayes_features") {
  
  ctx <- .load_diagnostic_context_for_dual_ranking(
    scenario_dir = scenario_dir,
    diag_target = diag_target
  )
  
  weights_long <- .get_mofa_weights_all_features(
    mofa_model = ctx$mofa_model,
    factor_cols = ctx$factor_cols
  )
  
  feature_map <- .make_feature_map_all_mofa(
    diagnostic_objects = ctx$diagnostic_objects,
    weights_long = weights_long,
    final_biomarkers = ctx$final_biomarkers
  )
  
  df_observed <- .make_observed_mofa_wide(
    ctx = ctx,
    feature_map = feature_map
  )
  
  df_reconstructed <- .make_reconstructed_mofa_wide(
    ctx = ctx,
    feature_map = feature_map,
    weights_long = weights_long
  )
  
  aligned <- .align_observed_and_reconstructed_wide(
    df_observed = df_observed,
    df_reconstructed = df_reconstructed
  )
  
  feature_map <- feature_map %>%
    dplyr::filter(feature_model %in% aligned$feature_common)
  
  list(
    ctx = ctx,
    feature_map = feature_map,
    df_observed = aligned$observed,
    df_reconstructed = aligned$reconstructed,
    sample_common = aligned$sample_common,
    feature_common = aligned$feature_common
  )
}

## ============================================================
## TEST BAYESIANO/WELCH POR CONTRASTE
## Misma matemática que definiste.
## ============================================================

contrast_t_test <- function(datas, factor_, response, contrast_expr) {
  
  if (!all(c(factor_, response) %in% colnames(datas))) {
    stop("El factor y la variable de respuesta deben estar en las columnas del dataframe.")
  }
  
  datas <- datas %>%
    dplyr::select(dplyr::all_of(c(factor_, response))) %>%
    dplyr::filter(!is.na(.data[[factor_]]), !is.na(.data[[response]]))
  
  datas[[factor_]] <- factor(datas[[factor_]])
  
  contrast_clean <- gsub(" ", "", contrast_expr)
  groups <- unlist(strsplit(contrast_clean, "-"))
  
  if (length(groups) != 2) {
    warning("Contraste no reconocido: ", contrast_expr)
    return(list(statistic = NA_real_, p.value = NA_real_))
  }
  
  g1 <- groups[1]
  g2 <- groups[2]
  
  if (!all(c(g1, g2) %in% levels(datas[[factor_]]))) {
    warning("Los grupos del contraste no están en el factor: ", contrast_expr)
    return(list(statistic = NA_real_, p.value = NA_real_))
  }
  
  x <- datas[[response]][datas[[factor_]] == g1]
  y <- datas[[response]][datas[[factor_]] == g2]
  
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  if (length(x) < 2 || length(y) < 2) {
    warning("Muy pocas observaciones para ", response, " en contraste ", contrast_expr)
    return(list(statistic = NA_real_, p.value = NA_real_))
  }
  
  vx <- stats::var(x)
  vy <- stats::var(y)
  
  if (!is.finite(vx) || !is.finite(vy) || vx == 0 || vy == 0) {
    warning("Varianza cero/no finita para ", response, " en contraste ", contrast_expr)
    
    tt <- tryCatch(
      stats::t.test(x, y, var.equal = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(tt)) {
      return(list(statistic = NA_real_, p.value = NA_real_))
    }
    
    return(list(
      statistic = unname(tt$statistic),
      p.value = unname(tt$p.value)
    ))
  }
  
  m_x <- mean(x)
  m_y <- mean(y)
  
  m <- c(m_x, m_y)
  n0 <- c(1 / vx, 1 / vy)
  sig.med <- median(c(sd(x), sd(y)))
  kappa <- 1
  
  res <- tryCatch(
    Bolstad::bayes.t.test(
      x,
      y,
      var.equal = FALSE,
      prior = "joint.conj",
      m = m,
      n0 = n0,
      sig.med = sig.med,
      kappa = kappa
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(
      statistic = unname(res$statistic),
      p.value = unname(res$p.value)
    ))
  }
  
  tt <- tryCatch(
    stats::t.test(x, y, var.equal = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(tt)) {
    return(list(statistic = NA_real_, p.value = NA_real_))
  }
  
  list(
    statistic = unname(tt$statistic),
    p.value = unname(tt$p.value)
  )
}

.compute_contrast_descriptives <- function(datas,
                                           factor_,
                                           response,
                                           contrast_expr) {
  
  contrast_clean <- gsub(" ", "", contrast_expr)
  groups <- unlist(strsplit(contrast_clean, "-"))
  
  if (length(groups) != 2) {
    return(tibble::tibble())
  }
  
  g1 <- groups[1]
  g2 <- groups[2]
  
  d <- datas %>%
    dplyr::select(dplyr::all_of(c(factor_, response))) %>%
    dplyr::filter(
      !is.na(.data[[factor_]]),
      !is.na(.data[[response]])
    )
  
  x <- d[[response]][d[[factor_]] == g1]
  y <- d[[response]][d[[factor_]] == g2]
  
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  mean_g1 <- mean(x, na.rm = TRUE)
  mean_g2 <- mean(y, na.rm = TRUE)
  
  median_g1 <- median(x, na.rm = TRUE)
  median_g2 <- median(y, na.rm = TRUE)
  
  sd_g1 <- sd(x, na.rm = TRUE)
  sd_g2 <- sd(y, na.rm = TRUE)
  
  pooled_sd <- sqrt(
    ((length(x) - 1) * sd_g1^2 + (length(y) - 1) * sd_g2^2) /
      pmax(length(x) + length(y) - 2, 1)
  )
  
  estimate_mean <- mean_g1 - mean_g2
  estimate_median <- median_g1 - median_g2
  
  cohen_d <- ifelse(
    is.finite(pooled_sd) && pooled_sd > 0,
    estimate_mean / pooled_sd,
    NA_real_
  )
  
  tibble::tibble(
    group_1 = g1,
    group_2 = g2,
    n_group_1 = length(x),
    n_group_2 = length(y),
    mean_group_1 = mean_g1,
    mean_group_2 = mean_g2,
    median_group_1 = median_g1,
    median_group_2 = median_g2,
    sd_group_1 = sd_g1,
    sd_group_2 = sd_g2,
    estimate_mean = estimate_mean,
    estimate_median = estimate_median,
    cohen_d = cohen_d,
    direction = dplyr::case_when(
      is.finite(estimate_mean) & estimate_mean > 0 ~ paste0("higher_in_", g1),
      is.finite(estimate_mean) & estimate_mean < 0 ~ paste0("higher_in_", g2),
      TRUE ~ "no_clear_direction"
    )
  )
}

rank_features_by_contrast <- function(df_wide,
                                      feature_map,
                                      signal_mode,
                                      y_levels = c("NP", "FD", "SO")) {
  
  df_rank <- df_wide %>%
    dplyr::mutate(
      y = factor(y, levels = y_levels),
      y_nonNP = dplyr::case_when(
        y == "NP" ~ "NP",
        y %in% c("FD", "SO") ~ "NONNP",
        TRUE ~ NA_character_
      ),
      y_nonNP = factor(y_nonNP, levels = c("NONNP", "NP"))
    )
  
  feature_cols <- setdiff(
    colnames(df_rank),
    c("sample_id", "split", "y", "y_nonNP")
  )
  
  contrast_defs <- list(
    FD_vs_NP = list(
      factor_col = "y",
      contrast_expr = "FD-NP"
    ),
    SO_vs_NP = list(
      factor_col = "y",
      contrast_expr = "SO-NP"
    ),
    SO_vs_FD = list(
      factor_col = "y",
      contrast_expr = "SO-FD"
    ),
    NONNP_vs_NP = list(
      factor_col = "y_nonNP",
      contrast_expr = "NONNP-NP"
    )
  )
  
  out <- list()
  counter <- 0L
  
  for (feat in feature_cols) {
    
    for (contrast_name in names(contrast_defs)) {
      
      factor_col <- contrast_defs[[contrast_name]]$factor_col
      contrast_expr <- contrast_defs[[contrast_name]]$contrast_expr
      
      dat_feat <- df_rank %>%
        dplyr::select(
          dplyr::all_of(c(factor_col, feat))
        )
      
      test_res <- suppressWarnings(
        contrast_t_test(
          datas = dat_feat,
          factor_ = factor_col,
          response = feat,
          contrast_expr = contrast_expr
        )
      )
      
      desc_res <- .compute_contrast_descriptives(
        datas = dat_feat,
        factor_ = factor_col,
        response = feat,
        contrast_expr = contrast_expr
      )
      
      counter <- counter + 1L
      
      out[[counter]] <- tibble::tibble(
        signal_mode = signal_mode,
        feature_model = feat,
        contrast = contrast_name,
        factor_used = factor_col,
        contrast_expr = contrast_expr,
        statistic = as.numeric(test_res$statistic),
        p_value = as.numeric(test_res$p.value)
      ) %>%
        dplyr::bind_cols(desc_res)
    }
  }
  
  dplyr::bind_rows(out) %>%
    dplyr::group_by(signal_mode, contrast) %>%
    dplyr::mutate(
      q_value = p.adjust(p_value, method = "BH")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      rank_score = dplyr::case_when(
        is.finite(statistic) ~ statistic,
        is.finite(p_value) & is.finite(estimate_mean) ~
          sign(estimate_mean) * (-log10(pmax(p_value, 1e-300))),
        TRUE ~ NA_real_
      ),
      abs_rank_score = abs(rank_score),
      minus_log10_p = -log10(pmax(p_value, 1e-300)),
      minus_log10_q = -log10(pmax(q_value, 1e-300)),
      nominal_p_lt_0_05 = is.finite(p_value) & p_value < 0.05,
      q_lt_0_20 = is.finite(q_value) & q_value < 0.20,
      q_lt_0_10 = is.finite(q_value) & q_value < 0.10
    ) %>%
    dplyr::left_join(
      feature_map,
      by = "feature_model"
    ) %>%
    dplyr::arrange(
      signal_mode,
      contrast,
      dplyr::desc(abs_rank_score),
      p_value
    ) %>%
    dplyr::select(
      signal_mode,
      feature_model,
      feature_label,
      view_code,
      mofa_feature,
      final_biomarker,
      annotation_eligible,
      annotation_route,
      contrast,
      factor_used,
      contrast_expr,
      group_1,
      group_2,
      n_group_1,
      n_group_2,
      mean_group_1,
      mean_group_2,
      median_group_1,
      median_group_2,
      sd_group_1,
      sd_group_2,
      estimate_mean,
      estimate_median,
      cohen_d,
      direction,
      statistic,
      rank_score,
      abs_rank_score,
      p_value,
      q_value,
      minus_log10_p,
      minus_log10_q,
      nominal_p_lt_0_05,
      q_lt_0_20,
      q_lt_0_10
    )
}

plot_dual_ranking_audit <- function(ranking_observed,
                                    ranking_reconstructed,
                                    outdir) {
  
  plots_dir <- file.path(outdir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  ranking_all <- dplyr::bind_rows(
    ranking_observed,
    ranking_reconstructed
  ) %>%
    dplyr::mutate(
      feature_plot = paste0(feature_label, " [", view_code, "]"),
      signal_mode = factor(
        signal_mode,
        levels = c("observed", "reconstructed")
      )
    )
  
  ## ------------------------------------------------------------
  ## Plot 1. Volcano/ranking por señal, vista y contraste.
  ## ------------------------------------------------------------
  
  p_volcano <- ranking_all %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = rank_score,
        y = minus_log10_p
      )
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        shape = final_biomarker,
        size = abs_rank_score
      ),
      alpha = 0.75
    ) +
    ggplot2::facet_grid(
      signal_mode + view_code ~ contrast,
      scales = "free"
    ) +
    ggplot2::labs(
      title = "Auditoría visual del ranking diferencial",
      subtitle = "Observed = sin integrar; reconstructed = señal integrada MOFA. Biomarcadores finales resaltados por shape.",
      x = "Rank score firmado",
      y = "-log10(p)",
      shape = "Biomarcador final",
      size = "|rank score|"
    ) +
    ggplot2::theme_minimal()
  
  .save_plot_publication(
    filename_base = file.path(
      plots_dir,
      "01_dual_ranking_volcano_all_features"
    ),
    plot = p_volcano,
    width = 16,
    height = 11,
    dpi = 600
  )
  
  ## ------------------------------------------------------------
  ## Plot 2. Comparación observed vs reconstructed en biomarcadores.
  ## ------------------------------------------------------------
  
  final_wide <- ranking_all %>%
    dplyr::filter(final_biomarker) %>%
    dplyr::select(
      feature_model,
      feature_plot,
      view_code,
      contrast,
      signal_mode,
      rank_score,
      p_value,
      q_value
    ) %>%
    tidyr::pivot_wider(
      names_from = signal_mode,
      values_from = c(rank_score, p_value, q_value)
    ) %>%
    dplyr::mutate(
      concordant_direction = sign(rank_score_observed) == sign(rank_score_reconstructed),
      min_p = pmin(p_value_observed, p_value_reconstructed, na.rm = TRUE),
      minus_log10_min_p = -log10(pmax(min_p, 1e-300))
    )
  
  if (nrow(final_wide) > 0) {
    
    p_compare <- final_wide %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = rank_score_observed,
          y = rank_score_reconstructed
        )
      ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed"
      ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      ggplot2::geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dotted"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          shape = view_code,
          size = minus_log10_min_p
        ),
        alpha = 0.85
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = feature_plot),
        size = 2.7,
        vjust = -0.7,
        check_overlap = TRUE
      ) +
      ggplot2::facet_wrap(
        ~ contrast,
        scales = "free"
      ) +
      ggplot2::labs(
        title = "Concordancia del ranking: observed vs reconstructed",
        subtitle = "Solo biomarcadores finales. Cuadrantes concordantes indican misma dirección diferencial.",
        x = "Rank score observed / sin integrar",
        y = "Rank score reconstructed / integrado",
        shape = "Vista",
        size = "-log10(min p)"
      ) +
      ggplot2::theme_minimal()
    
    .save_plot_publication(
      filename_base = file.path(
        plots_dir,
        "02_observed_vs_reconstructed_final_biomarkers"
      ),
      plot = p_compare,
      width = 14,
      height = 9,
      dpi = 600
    )
  }
  
  ## ------------------------------------------------------------
  ## Plot 3. Heatmap de rank score para biomarcadores finales.
  ## ------------------------------------------------------------
  
  p_heat <- ranking_all %>%
    dplyr::filter(final_biomarker) %>%
    dplyr::mutate(
      signal_contrast = paste(signal_mode, contrast, sep = " | "),
      feature_plot = factor(
        feature_plot,
        levels = rev(unique(feature_plot[order(view_code, feature_plot)]))
      )
    ) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = signal_contrast,
        y = feature_plot,
        fill = rank_score
      )
    ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", rank_score)),
      size = 2.7
    ) +
    ggplot2::facet_grid(
      view_code ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::labs(
      title = "Rank score por biomarcador final",
      subtitle = "Comparación visual de señal sin integrar vs señal integrada.",
      x = "Señal | contraste",
      y = "Biomarcador final",
      fill = "Rank score"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text.y = ggplot2::element_text(angle = 0)
    )
  
  .save_plot_publication(
    filename_base = file.path(
      plots_dir,
      "03_final_biomarker_rankscore_heatmap_dual_signal"
    ),
    plot = p_heat,
    width = 16,
    height = max(7, 0.35 * length(unique(ranking_all$feature_model[ranking_all$final_biomarker]))),
    dpi = 600
  )
  
  ## ------------------------------------------------------------
  ## Plot 4. Top biomarcadores finales por ranking.
  ## ------------------------------------------------------------
  
  p_top <- ranking_all %>%
    dplyr::filter(final_biomarker) %>%
    dplyr::group_by(signal_mode, contrast) %>%
    dplyr::slice_max(
      order_by = abs_rank_score,
      n = 15,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      feature_plot = paste0(feature_label, " [", view_code, "]"),
      feature_plot = stats::reorder(feature_plot, abs_rank_score)
    ) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = feature_plot,
        y = abs_rank_score
      )
    ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_grid(
      signal_mode ~ contrast,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::labs(
      title = "Top biomarcadores finales por magnitud del ranking",
      subtitle = "Auditoría de qué biomarcadores dominan cada flujo y contraste.",
      x = "Biomarcador final",
      y = "|rank score|"
    ) +
    ggplot2::theme_minimal()
  
  .save_plot_publication(
    filename_base = file.path(
      plots_dir,
      "04_top_final_biomarkers_by_dual_ranking"
    ),
    plot = p_top,
    width = 16,
    height = 10,
    dpi = 600
  )
  
  invisible(
    list(
      ranking_all = ranking_all,
      final_wide = final_wide,
      plots_dir = plots_dir
    )
  )
}


.make_dual_ranking_input_from_bidir_results <- function(all_bidir_results) {
  
  if (is.null(all_bidir_results$observed)) {
    stop("No existe all_bidir_results$observed. Ejecuta SIGNAL_MODE=both.")
  }
  
  if (is.null(all_bidir_results$reconstructed)) {
    stop("No existe all_bidir_results$reconstructed. Ejecuta SIGNAL_MODE=both.")
  }
  
  inp_obs <- all_bidir_results$observed$input
  inp_rec <- all_bidir_results$reconstructed$input
  
  if (is.null(inp_obs$df_all_observed)) {
    stop("No encuentro inp_observed$df_all_observed.")
  }
  
  if (is.null(inp_rec$df_all_reconstructed)) {
    stop("No encuentro inp_reconstructed$df_all_reconstructed.")
  }
  
  df_observed <- inp_obs$df_all_observed
  df_reconstructed <- inp_rec$df_all_reconstructed
  
  feature_obs <- setdiff(
    colnames(df_observed),
    c("sample_id", "split", "y")
  )
  
  feature_rec <- setdiff(
    colnames(df_reconstructed),
    c("sample_id", "split", "y")
  )
  
  feature_common <- intersect(feature_obs, feature_rec)
  
  if (length(feature_common) < 3) {
    stop(
      "Muy pocas features comunes entre observed y reconstructed: ",
      length(feature_common)
    )
  }
  
  sample_common <- intersect(
    df_observed$sample_id,
    df_reconstructed$sample_id
  )
  
  if (length(sample_common) < 6) {
    
    if (nrow(df_observed) == nrow(df_reconstructed)) {
      
      warning(
        "No hubo match suficiente por sample_id. ",
        "Alineo observed/reconstructed por orden de filas."
      )
      
      df_observed$sample_id <- paste0("sample_", seq_len(nrow(df_observed)))
      df_reconstructed$sample_id <- paste0("sample_", seq_len(nrow(df_reconstructed)))
      
      sample_common <- df_observed$sample_id
      
    } else {
      
      stop(
        "Muy pocas muestras comunes entre observed y reconstructed: ",
        length(sample_common)
      )
    }
  }
  
  df_observed2 <- df_observed %>%
    dplyr::filter(sample_id %in% sample_common) %>%
    dplyr::arrange(match(sample_id, sample_common)) %>%
    dplyr::select(
      sample_id,
      split,
      y,
      dplyr::all_of(feature_common)
    )
  
  df_reconstructed2 <- df_reconstructed %>%
    dplyr::filter(sample_id %in% sample_common) %>%
    dplyr::arrange(match(sample_id, sample_common)) %>%
    dplyr::select(
      sample_id,
      split,
      y,
      dplyr::all_of(feature_common)
    )
  
  if (!identical(as.character(df_observed2$y), as.character(df_reconstructed2$y))) {
    stop("Los grupos y no coinciden entre observed y reconstructed tras alinear muestras.")
  }
  
  ## ------------------------------------------------------------
  ## Feature map.
  ## Preferimos el mapa del flujo observado.
  ## ------------------------------------------------------------
  
  feature_map <- inp_obs$feature_map %>%
    dplyr::as_tibble() %>%
    dplyr::filter(feature_model %in% feature_common)
  
  missing_features <- setdiff(feature_common, feature_map$feature_model)
  
  if (length(missing_features) > 0) {
    
    fallback_map <- tibble::tibble(
      feature_model = missing_features,
      view_code = sub("_.*$", "", missing_features)
    )
    
    feature_map <- dplyr::bind_rows(
      feature_map,
      fallback_map
    )
  }
  
  ## ------------------------------------------------------------
  ## Asegurar columna mofa_feature.
  ## inp_obs$feature_map normalmente NO tiene mofa_feature.
  ## inp_rec$feature_map sí suele tenerla.
  ## No usar dplyr::if_else() aquí porque evalúa ambas ramas.
  ## ------------------------------------------------------------
  
  if (!"mofa_feature" %in% colnames(feature_map)) {
    feature_map$mofa_feature <- feature_map$feature_model
  }
  
  feature_map <- feature_map %>%
    dplyr::mutate(
      feature_model = as.character(feature_model),
      view_code = .standardize_view_global(view_code),
      mofa_feature = as.character(mofa_feature),
      feature_label = sub("^[^_]+_", "", feature_model),
      
      ## En este framework, las features comunes son las biomarcadoras
      ## finales que se están examinando.
      final_biomarker = TRUE,
      
      annotation_eligible = view_code %in% c("tx", "pr", "cl"),
      annotation_route = dplyr::case_when(
        view_code == "tx" ~ "gene_level_GO_Reactome_MSigDB",
        view_code == "pr" ~ "protein_level_UniProt_GO_Reactome",
        view_code == "cl" ~ "clinical_ontology_EFO_HPO_MeSH_OLS",
        view_code == "me" ~ "not_annotated_metabolite_id_missing",
        TRUE ~ "unknown"
      )
    ) %>%
    dplyr::distinct(
      feature_model,
      .keep_all = TRUE
    ) %>%
    dplyr::select(
      feature_model,
      feature_label,
      view_code,
      mofa_feature,
      final_biomarker,
      annotation_eligible,
      annotation_route
    )
  
  list(
    df_observed = df_observed2,
    df_reconstructed = df_reconstructed2,
    feature_map = feature_map,
    feature_common = feature_common,
    sample_common = sample_common
  )
}


run_mofa_dual_feature_differential_ranking <- function(scenario_dir,
                                                       diag_target = "diagnostico_bayes_features",
                                                       all_bidir_results,
                                                       outdir = NULL) {
  
  scenario_dir <- normalizePath(path.expand(scenario_dir), mustWork = TRUE)
  
  if (is.null(outdir)) {
    outdir <- file.path(
      scenario_dir,
      diag_target,
      "mofa_dual_feature_differential_ranking"
    )
  }
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\n============================================================\n")
  cat("RANKING DIFERENCIAL DUAL: OBSERVED vs RECONSTRUCTED\n")
  cat("Fuente observed     : all_bidir_results$observed$input$df_all_observed\n")
  cat("Fuente reconstructed: all_bidir_results$reconstructed$input$df_all_reconstructed\n")
  cat("Outdir:\n")
  cat(outdir, "\n")
  cat("============================================================\n")
  
  inp <- .make_dual_ranking_input_from_bidir_results(
    all_bidir_results = all_bidir_results
  )
  
  ranking_observed <- rank_features_by_contrast(
    df_wide = inp$df_observed,
    feature_map = inp$feature_map,
    signal_mode = "observed"
  )
  
  ranking_reconstructed <- rank_features_by_contrast(
    df_wide = inp$df_reconstructed,
    feature_map = inp$feature_map,
    signal_mode = "reconstructed"
  )
  
  ranking_observed_final <- ranking_observed %>%
    dplyr::filter(final_biomarker)
  
  ranking_reconstructed_final <- ranking_reconstructed %>%
    dplyr::filter(final_biomarker)
  
  write.csv(
    ranking_observed,
    file.path(outdir, "20_observed_all_feature_differential_ranking.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ranking_reconstructed,
    file.path(outdir, "21_reconstructed_all_feature_differential_ranking.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ranking_observed_final,
    file.path(outdir, "22_observed_final_biomarker_differential_ranking.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ranking_reconstructed_final,
    file.path(outdir, "23_reconstructed_final_biomarker_differential_ranking.csv"),
    row.names = FALSE
  )
  
  write.csv(
    inp$feature_map,
    file.path(outdir, "24_feature_map_for_annotation.csv"),
    row.names = FALSE
  )
  
  plot_obj <- plot_dual_ranking_audit(
    ranking_observed = ranking_observed,
    ranking_reconstructed = ranking_reconstructed,
    outdir = outdir
  )
  
  cat("\nRanking diferencial dual guardado:\n")
  cat(" - 20_observed_all_feature_differential_ranking.csv\n")
  cat(" - 21_reconstructed_all_feature_differential_ranking.csv\n")
  cat(" - 22_observed_final_biomarker_differential_ranking.csv\n")
  cat(" - 23_reconstructed_final_biomarker_differential_ranking.csv\n")
  cat(" - 24_feature_map_for_annotation.csv\n")
  cat(" - plots/01_dual_ranking_volcano_all_features.png/pdf\n")
  cat(" - plots/02_observed_vs_reconstructed_final_biomarkers.png/pdf\n")
  cat(" - plots/03_final_biomarker_rankscore_heatmap_dual_signal.png/pdf\n")
  cat(" - plots/04_top_final_biomarkers_by_dual_ranking.png/pdf\n")
  
  invisible(
    list(
      input = inp,
      ranking_observed = ranking_observed,
      ranking_reconstructed = ranking_reconstructed,
      ranking_observed_final = ranking_observed_final,
      ranking_reconstructed_final = ranking_reconstructed_final,
      plots = plot_obj,
      outdir = outdir
    )
  )
}

## ============================================================
## ANOTACIÓN BIOLÓGICA DINÁMICA CON BASES DE DATOS
## ------------------------------------------------------------
## No hardcodea features concretas.
##
## Entradas principales:
##   - 24_feature_map_for_annotation.csv
##   - 20_observed_all_feature_differential_ranking.csv
##   - 21_reconstructed_all_feature_differential_ranking.csv
##   - módulos observed/reconstructed si existen
##
## Reglas:
##   tx/pr:
##     - usa feature_label y feature_model como queries
##     - intenta resolver IDs con gprofiler2::gconvert()
##     - ejecuta enrichment con gprofiler2::gost()
##
##   cl:
##     - usa query derivada automáticamente del feature_label
##     - busca en OLS
##
##   me:
##     - no se anota biológicamente sin ID químico
##     - queda explícitamente como unidentified_metabolomic_signal
##
## Archivos opcionales para mejorar mapeo sin hardcodear:
##   database_biological_annotation/manual_maps/gene_manual_map.csv
##   database_biological_annotation/manual_maps/protein_manual_map.csv
##   database_biological_annotation/manual_maps/clinical_manual_map.csv
##   database_biological_annotation/manual_maps/metabolite_manual_map.csv
## ============================================================

.collapse_unique_annotation <- function(x,
                                        sep = " | ",
                                        max_n = 8) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & x != "" & x != "NA"]
  if (length(x) == 0) return(NA_character_)
  paste(head(x, max_n), collapse = sep)
}

.safe_min_annotation <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x, na.rm = TRUE)
}

.safe_max_annotation <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x, na.rm = TRUE)
}

.safe_read_csv_annotation <- function(path) {
  
  if (!file.exists(path)) return(NULL)
  
  x <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
  
  tibble::as_tibble(x)
}
## ============================================================
## REFERENCE MAPS GLOBALES
## ------------------------------------------------------------
## Los mapas en reference_maps/ son maestros del proyecto.
## Los mapas del escenario en manual_maps/ solo actúan como override.
## ============================================================

.find_reference_maps_dir <- function(start_dir = getwd()) {
  
  start_dir <- normalizePath(path.expand(start_dir), mustWork = FALSE)
  
  candidates <- unique(c(
    file.path(getwd(), "reference_maps"),
    file.path(dirname(getwd()), "reference_maps"),
    file.path(start_dir, "reference_maps")
  ))
  
  cur <- start_dir
  
  for (i in seq_len(12)) {
    candidates <- unique(c(candidates, file.path(cur, "reference_maps")))
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  
  candidates <- candidates[file.exists(candidates) & dir.exists(candidates)]
  
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  
  normalizePath(candidates[1], mustWork = TRUE)
}

.row_has_any_nonempty <- function(df, cols) {
  
  cols <- intersect(cols, colnames(df))
  
  if (length(cols) == 0) {
    return(rep(FALSE, nrow(df)))
  }
  
  mat <- df[, cols, drop = FALSE]
  
  apply(
    mat,
    1,
    function(z) {
      z <- as.character(z)
      any(!is.na(z) & trimws(z) != "")
    }
  )
}

.align_columns_for_bind <- function(x, y) {
  
  x <- tibble::as_tibble(x)
  y <- tibble::as_tibble(y)
  
  all_cols <- union(colnames(x), colnames(y))
  
  for (cc in setdiff(all_cols, colnames(x))) {
    x[[cc]] <- NA_character_
  }
  
  for (cc in setdiff(all_cols, colnames(y))) {
    y[[cc]] <- NA_character_
  }
  
  x <- x[, all_cols, drop = FALSE]
  y <- y[, all_cols, drop = FALSE]
  
  ## Los mapas manuales son metadatos textuales.
  ## Evita errores tipo:
  ## logical vs character en database_id, ontology_id, notes, etc.
  x <- x %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        ~ as.character(.x)
      )
    )
  
  y <- y %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        ~ as.character(.x)
      )
    )
  
  list(
    x = x,
    y = y
  )
}

.merge_reference_master_csv <- function(master_file,
                                        target_file,
                                        key_col = "feature_model",
                                        override_cols = character()) {
  
  if (!file.exists(master_file)) {
    return(invisible(NULL))
  }
  
  master <- .safe_read_csv_annotation(master_file)
  
  if (is.null(master) || nrow(master) == 0) {
    return(invisible(NULL))
  }
  
  master <- tibble::as_tibble(master)
  
  if (!key_col %in% colnames(master)) {
    warning("Mapa maestro sin columna ", key_col, ": ", master_file)
    return(invisible(NULL))
  }
  
  if (!file.exists(target_file)) {
    
    .safe_write_csv_annotation(
      master,
      target_file
    )
    
    return(invisible(target_file))
  }
  
  local <- .safe_read_csv_annotation(target_file)
  
  if (is.null(local) || nrow(local) == 0) {
    
    .safe_write_csv_annotation(
      master,
      target_file
    )
    
    return(invisible(target_file))
  }
  
  local <- tibble::as_tibble(local)
  
  if (!key_col %in% colnames(local)) {
    warning("Mapa local sin columna ", key_col, ": ", target_file)
    return(invisible(NULL))
  }
  
  aligned <- .align_columns_for_bind(local, master)
  local  <- aligned$x
  master <- aligned$y
  
  ## Solo las filas locales con información real deben sobreescribir
  ## al mapa maestro. Las plantillas vacías no deben bloquear el mapa global.
  local_override <- local[
    .row_has_any_nonempty(local, override_cols),
    ,
    drop = FALSE
  ]
  
  local_empty <- local[
    !.row_has_any_nonempty(local, override_cols),
    ,
    drop = FALSE
  ]
  
  merged <- dplyr::bind_rows(
    local_override,
    master,
    local_empty
  ) %>%
    dplyr::mutate(
      "{key_col}" := as.character(.data[[key_col]])
    ) %>%
    dplyr::filter(
      !is.na(.data[[key_col]]),
      .data[[key_col]] != ""
    ) %>%
    dplyr::distinct(
      .data[[key_col]],
      .keep_all = TRUE
    )
  
  .safe_write_csv_annotation(
    merged,
    target_file
  )
  
  invisible(target_file)
}

.stage_reference_maps_into_manual_dir <- function(manual_dir) {
  
  ref_dir <- .find_reference_maps_dir(start_dir = dirname(manual_dir))
  
  if (is.na(ref_dir) || !nzchar(ref_dir)) {
    
    cat("\n[REFERENCE_MAPS] No encontré carpeta reference_maps/. Se usarán solo manual_maps locales.\n")
    return(invisible(NULL))
  }
  
  cat("\n[REFERENCE_MAPS] Usando mapas maestros desde:\n")
  cat("  ", ref_dir, "\n")
  
  dir.create(manual_dir, recursive = TRUE, showWarnings = FALSE)
  
  protein_master <- file.path(ref_dir, "protein_manual_map_master.csv")
  clinical_master <- file.path(ref_dir, "clinical_manual_map_master.csv")
  
  protein_target <- file.path(manual_dir, "protein_manual_map.csv")
  clinical_target <- file.path(manual_dir, "clinical_manual_map.csv")
  
  .merge_reference_master_csv(
    master_file = protein_master,
    target_file = protein_target,
    key_col = "feature_model",
    override_cols = c(
      "mapped_gene_symbol",
      "mapped_uniprot_id",
      "database_id",
      "mapping_confidence",
      "notes"
    )
  )
  
  .merge_reference_master_csv(
    master_file = clinical_master,
    target_file = clinical_target,
    key_col = "feature_model",
    override_cols = c(
      "clinical_query",
      "ontology_source",
      "ontology_id",
      "ontology_label",
      "semantic_category",
      "biological_axis",
      "mapping_confidence",
      "notes"
    )
  )
  
  ## Copia trazable de los mapas maestros usados.
  used_dir <- file.path(manual_dir, "reference_maps_used")
  dir.create(used_dir, recursive = TRUE, showWarnings = FALSE)
  
  files_to_copy <- c(
    "clinical_manual_map_master.csv",
    "clinical_master_alias_map.tsv",
    "clinical_master_map.tsv",
    "protein_manual_map_master.csv",
    "protein_manual_map_by_label.csv",
    "protein_master_alias_map.tsv",
    "protein_master_map.tsv",
    "protein_enrichment_gene_sets.tsv"
  )
  
  for (ff in files_to_copy) {
    src <- file.path(ref_dir, ff)
    if (file.exists(src)) {
      file.copy(
        from = src,
        to = file.path(used_dir, ff),
        overwrite = TRUE
      )
    }
  }
  
  invisible(ref_dir)
}



.safe_write_csv_annotation <- function(x, path) {
  
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  cat("\n[WRITE_CSV] Intentando guardar:\n")
  cat("  Archivo:", path, "\n")
  cat("  Clase objeto:", paste(class(x), collapse = " / "), "\n")
  
  if (is.null(x)) {
    cat("  Estado: objeto NULL. Se guarda tabla vacía.\n")
    x <- tibble::tibble()
  }
  
  x <- tibble::as_tibble(x)
  
  cat("  Dimensión:", nrow(x), "x", ncol(x), "\n")
  
  list_cols <- names(x)[vapply(x, is.list, logical(1))]
  
  if (length(list_cols) > 0) {
    
    cat("  COLUMNAS TIPO LIST DETECTADAS:\n")
    cat("  -", paste(list_cols, collapse = "\n  - "), "\n")
    
    ## Guardar diagnóstico crudo antes de convertir
    debug_rds <- paste0(path, ".raw_debug.rds")
    saveRDS(x, debug_rds)
    cat("  Objeto crudo guardado en:\n")
    cat("  ", debug_rds, "\n")
    
    ## Convertir columnas list a texto para que write.csv no falle
    x <- x %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(list_cols),
          function(z) {
            vapply(
              z,
              function(elem) {
                elem <- unlist(elem)
                elem <- elem[!is.na(elem)]
                if (length(elem) == 0) return(NA_character_)
                paste(as.character(elem), collapse = " | ")
              },
              character(1)
            )
          }
        )
      )
    
    cat("  Columnas list convertidas a texto.\n")
  }
  
  tryCatch(
    {
      write.csv(
        x,
        path,
        row.names = FALSE,
        na = ""
      )
      cat("  Estado: OK\n")
    },
    error = function(e) {
      cat("\n[ERROR WRITE_CSV]\n")
      cat("  Archivo:", path, "\n")
      cat("  Mensaje:", conditionMessage(e), "\n")
      cat("  Columnas:\n")
      print(str(x))
      stop(e)
    }
  )
  
  invisible(path)
}
.clean_query_text_annotation <- function(x) {
  
  x <- as.character(x)
  x <- trimws(x)
  
  x <- gsub("^tx_", "", x)
  x <- gsub("^pr_", "", x)
  x <- gsub("^cl_", "", x)
  x <- gsub("^me_", "", x)
  
  x <- gsub("_0M$", "", x)
  x <- gsub("_M0$", "", x)
  x <- gsub("_baseline$", "", x, ignore.case = TRUE)
  
  x <- gsub("[._]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

.make_identifier_candidates <- function(feature_model,
                                        feature_label,
                                        view_code) {
  
  fm <- as.character(feature_model)
  fl <- as.character(feature_label)
  
  candidates <- unique(c(
    fl,
    fm,
    gsub("^tx_", "", fm),
    gsub("^pr_", "", fm),
    gsub("^cl_", "", fm),
    gsub("^me_", "", fm),
    gsub("_", "-", fl),
    gsub("_", "", fl),
    gsub("-", "", fl),
    gsub("[^A-Za-z0-9-]", "", fl),
    .clean_query_text_annotation(fl)
  ))
  
  candidates <- trimws(candidates)
  candidates <- candidates[!is.na(candidates) & candidates != ""]
  
  tibble::tibble(
    feature_model = fm,
    feature_label = fl,
    view_code = view_code,
    query_candidate = candidates,
    candidate_rank = seq_along(candidates)
  )
}

.read_optional_manual_map <- function(path,
                                      expected_cols) {
  
  x <- .safe_read_csv_annotation(path)
  
  if (is.null(x)) {
    return(NULL)
  }
  
  missing_cols <- setdiff(expected_cols, colnames(x))
  
  if (length(missing_cols) > 0) {
    warning(
      "El archivo manual existe pero le faltan columnas: ",
      basename(path),
      " -> ",
      paste(missing_cols, collapse = ", ")
    )
    return(NULL)
  }
  
  tibble::as_tibble(x)
}

.write_manual_mapping_templates <- function(feature_map,
                                            manual_dir) {
  
  dir.create(manual_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## Primero cargar mapas maestros globales.
  ## Después, las plantillas locales se crean solo para lo que falte.
  .stage_reference_maps_into_manual_dir(manual_dir)
  
  
  fm <- tibble::as_tibble(feature_map) %>%
    dplyr::mutate(
      feature_model = as.character(feature_model),
      feature_label = as.character(feature_label),
      view_code = as.character(view_code)
    )
  
  gene_template <- fm %>%
    dplyr::filter(view_code == "tx") %>%
    dplyr::transmute(
      feature_model,
      feature_label,
      database_id = "",
      database_id_type = "",
      mapped_gene_symbol = "",
      mapping_confidence = "",
      notes = ""
    )
  
  protein_template <- fm %>%
    dplyr::filter(view_code == "pr") %>%
    dplyr::transmute(
      feature_model,
      feature_label,
      database_id = "",
      database_id_type = "",
      mapped_gene_symbol = "",
      mapped_uniprot_id = "",
      mapping_confidence = "",
      notes = ""
    )
  
  clinical_template <- fm %>%
    dplyr::filter(view_code == "cl") %>%
    dplyr::transmute(
      feature_model,
      feature_label,
      clinical_query = "",
      ontology_source = "",
      ontology_id = "",
      ontology_label = "",
      semantic_category = "",
      biological_axis = "",
      mapping_confidence = "",
      notes = ""
    )
  
  metabolite_template <- fm %>%
    dplyr::filter(view_code == "me") %>%
    dplyr::transmute(
      feature_model,
      feature_label,
      metabolite_name = "",
      HMDB_ID = "",
      KEGG_ID = "",
      ChEBI_ID = "",
      PubChem_CID = "",
      identification_level = "",
      mapping_confidence = "",
      notes = ""
    )
  
  gene_file <- file.path(manual_dir, "gene_manual_map.csv")
  protein_file <- file.path(manual_dir, "protein_manual_map.csv")
  clinical_file <- file.path(manual_dir, "clinical_manual_map.csv")
  metabolite_file <- file.path(manual_dir, "metabolite_manual_map.csv")
  
  if (!file.exists(gene_file)) {
    .safe_write_csv_annotation(gene_template, gene_file)
  }
  
  if (!file.exists(protein_file)) {
    .safe_write_csv_annotation(protein_template, protein_file)
  }
  
  if (!file.exists(clinical_file)) {
    .safe_write_csv_annotation(clinical_template, clinical_file)
  }
  
  if (!file.exists(metabolite_file)) {
    .safe_write_csv_annotation(metabolite_template, metabolite_file)
  }
}

.apply_optional_manual_maps <- function(feature_map,
                                        manual_dir) {
  
  fm <- tibble::as_tibble(feature_map) %>%
    dplyr::mutate(
      feature_model = as.character(feature_model),
      feature_label = as.character(feature_label),
      view_code = as.character(view_code)
    )
  
  gene_map <- .read_optional_manual_map(
    file.path(manual_dir, "gene_manual_map.csv"),
    c("feature_model", "mapped_gene_symbol")
  )
  
  protein_map <- .read_optional_manual_map(
    file.path(manual_dir, "protein_manual_map.csv"),
    c("feature_model", "mapped_gene_symbol")
  )
  
  clinical_map <- .read_optional_manual_map(
    file.path(manual_dir, "clinical_manual_map.csv"),
    c("feature_model", "clinical_query")
  )
  
  metabolite_map <- .read_optional_manual_map(
    file.path(manual_dir, "metabolite_manual_map.csv"),
    c("feature_model", "HMDB_ID", "KEGG_ID", "ChEBI_ID")
  )
  
  if (!is.null(gene_map)) {
    fm <- fm %>%
      dplyr::left_join(
        gene_map %>%
          dplyr::select(
            feature_model,
            manual_mapped_gene_symbol_tx = mapped_gene_symbol,
            manual_database_id_tx = dplyr::any_of("database_id"),
            manual_database_id_type_tx = dplyr::any_of("database_id_type"),
            manual_mapping_confidence_tx = dplyr::any_of("mapping_confidence"),
            manual_notes_tx = dplyr::any_of("notes")
          ),
        by = "feature_model"
      )
  }
  
  if (!is.null(protein_map)) {
    fm <- fm %>%
      dplyr::left_join(
        protein_map %>%
          dplyr::select(
            feature_model,
            manual_mapped_gene_symbol_pr = mapped_gene_symbol,
            manual_mapped_uniprot_id = dplyr::any_of("mapped_uniprot_id"),
            manual_database_id_pr = dplyr::any_of("database_id"),
            manual_database_id_type_pr = dplyr::any_of("database_id_type"),
            manual_mapping_confidence_pr = dplyr::any_of("mapping_confidence"),
            manual_notes_pr = dplyr::any_of("notes")
          ),
        by = "feature_model"
      )
  }
  
  if (!is.null(clinical_map)) {
    fm <- fm %>%
      dplyr::left_join(
        clinical_map %>%
          dplyr::select(
            feature_model,
            manual_clinical_query = clinical_query,
            manual_ontology_source = dplyr::any_of("ontology_source"),
            manual_ontology_id = dplyr::any_of("ontology_id"),
            manual_ontology_label = dplyr::any_of("ontology_label"),
            manual_semantic_category = dplyr::any_of("semantic_category"),
            manual_biological_axis = dplyr::any_of("biological_axis"),
            manual_clinical_confidence = dplyr::any_of("mapping_confidence"),
            manual_clinical_notes = dplyr::any_of("notes")
          ),
        by = "feature_model"
      )
  }
  
  if (!is.null(metabolite_map)) {
    fm <- fm %>%
      dplyr::left_join(
        metabolite_map %>%
          dplyr::select(
            feature_model,
            metabolite_name = dplyr::any_of("metabolite_name"),
            HMDB_ID,
            KEGG_ID,
            ChEBI_ID,
            PubChem_CID = dplyr::any_of("PubChem_CID"),
            metabolite_identification_level = dplyr::any_of("identification_level"),
            metabolite_mapping_confidence = dplyr::any_of("mapping_confidence")
          ),
        by = "feature_model"
      )
  }
  
  fm
}

.make_feature_annotation_base_dynamic <- function(feature_map,
                                                  manual_dir) {
  
  .write_manual_mapping_templates(
    feature_map = feature_map,
    manual_dir = manual_dir
  )
  
  fm <- .apply_optional_manual_maps(
    feature_map = feature_map,
    manual_dir = manual_dir
  )
  
  fm %>%
    dplyr::mutate(
      feature_model = as.character(feature_model),
      feature_label = as.character(feature_label),
      view_code = as.character(view_code),
      
      query_text_clean = .clean_query_text_annotation(feature_label),
      
      manual_gene_symbol = dplyr::case_when(
        view_code == "tx" &
          "manual_mapped_gene_symbol_tx" %in% colnames(.) ~
          as.character(manual_mapped_gene_symbol_tx),
        
        view_code == "pr" &
          "manual_mapped_gene_symbol_pr" %in% colnames(.) ~
          as.character(manual_mapped_gene_symbol_pr),
        
        TRUE ~ NA_character_
      ),
      
      manual_gene_symbol = dplyr::na_if(trimws(manual_gene_symbol), ""),
      
      database_query_type = dplyr::case_when(
        view_code == "tx" ~ "gene_or_transcript_identifier",
        view_code == "pr" ~ "protein_identifier_or_protein_label",
        view_code == "cl" ~ "clinical_variable_label",
        view_code == "me" ~ "metabolomic_feature_identifier",
        TRUE ~ "unknown"
      ),
      
      clinical_query = dplyr::case_when(
        view_code == "cl" &
          "manual_clinical_query" %in% colnames(.) &
          !is.na(manual_clinical_query) &
          manual_clinical_query != "" ~ as.character(manual_clinical_query),
        
        view_code == "cl" ~ query_text_clean,
        
        TRUE ~ NA_character_
      ),
      
      database_semantic_category = dplyr::case_when(
        view_code == "cl" &
          "manual_semantic_category" %in% colnames(.) &
          !is.na(manual_semantic_category) &
          manual_semantic_category != "" ~ as.character(manual_semantic_category),
        
        TRUE ~ NA_character_
      ),
      
      database_biological_axis = dplyr::case_when(
        view_code == "cl" &
          "manual_biological_axis" %in% colnames(.) &
          !is.na(manual_biological_axis) &
          manual_biological_axis != "" ~ as.character(manual_biological_axis),
        
        TRUE ~ NA_character_
      ),
      
      database_mapping_confidence = dplyr::case_when(
        view_code == "cl" &
          "manual_clinical_confidence" %in% colnames(.) &
          !is.na(manual_clinical_confidence) &
          manual_clinical_confidence != "" ~ as.character(manual_clinical_confidence),
        
        view_code == "pr" &
          "manual_mapping_confidence_pr" %in% colnames(.) &
          !is.na(manual_mapping_confidence_pr) &
          manual_mapping_confidence_pr != "" ~ as.character(manual_mapping_confidence_pr),
        
        view_code == "tx" &
          "manual_mapping_confidence_tx" %in% colnames(.) &
          !is.na(manual_mapping_confidence_tx) &
          manual_mapping_confidence_tx != "" ~ as.character(manual_mapping_confidence_tx),
        
        TRUE ~ NA_character_
      ),
      
      clinical_curated_no_ols = dplyr::case_when(
        view_code == "cl" &
          "manual_ontology_source" %in% colnames(.) &
          !is.na(manual_ontology_source) &
          manual_ontology_source == "curated_clinical_concept" ~ TRUE,
        
        TRUE ~ FALSE
      ),
      
      
      
      metabolite_has_database_id = dplyr::case_when(
        view_code != "me" ~ FALSE,
        "HMDB_ID" %in% colnames(.) &
          !is.na(HMDB_ID) &
          HMDB_ID != "" ~ TRUE,
        "KEGG_ID" %in% colnames(.) &
          !is.na(KEGG_ID) &
          KEGG_ID != "" ~ TRUE,
        "ChEBI_ID" %in% colnames(.) &
          !is.na(ChEBI_ID) &
          ChEBI_ID != "" ~ TRUE,
        TRUE ~ FALSE
      ),
      
      annotation_status_initial = dplyr::case_when(
        view_code %in% c("tx", "pr") ~ "ready_for_identifier_resolution",
        view_code == "cl" ~ "ready_for_ols",
        view_code == "me" & metabolite_has_database_id ~ "metabolite_id_available_manual",
        view_code == "me" ~ "unidentified_metabolomic_signal",
        TRUE ~ "unknown_or_unmapped"
      )
    )
}

.resolve_gene_protein_identifiers_gprofiler <- function(feature_annotation_base,
                                                        organism = "hsapiens") {
  
  gp_features <- feature_annotation_base %>%
    dplyr::filter(view_code %in% c("tx", "pr"))
  
  if (nrow(gp_features) == 0) {
    return(tibble::tibble())
  }
  
  candidates <- purrr::map_dfr(
    seq_len(nrow(gp_features)),
    function(i) {
      .make_identifier_candidates(
        feature_model = gp_features$feature_model[i],
        feature_label = gp_features$feature_label[i],
        view_code = gp_features$view_code[i]
      )
    }
  )
  
  manual_candidates <- gp_features %>%
    dplyr::filter(
      !is.na(manual_gene_symbol),
      manual_gene_symbol != ""
    ) %>%
    dplyr::transmute(
      feature_model,
      feature_label,
      view_code,
      query_candidate = manual_gene_symbol,
      candidate_rank = 0L
    )
  
  candidates <- dplyr::bind_rows(
    manual_candidates,
    candidates
  ) %>%
    dplyr::distinct(
      feature_model,
      query_candidate,
      .keep_all = TRUE
    ) %>%
    dplyr::arrange(feature_model, candidate_rank)
  
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    
    return(
      candidates %>%
        dplyr::group_by(feature_model) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          resolved_id = NA_character_,
          resolved_name = NA_character_,
          resolved_namespace = NA_character_,
          resolution_status = "gprofiler2_not_installed"
        )
    )
  }
  
  all_queries <- unique(candidates$query_candidate)
  all_queries <- all_queries[!is.na(all_queries) & all_queries != ""]
  
  conv <- tryCatch(
    gprofiler2::gconvert(
      query = all_queries,
      organism = organism,
      target = "ENSG",
      mthreshold = Inf,
      filter_na = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(conv) || nrow(conv) == 0) {
    
    return(
      candidates %>%
        dplyr::group_by(feature_model) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          resolved_id = NA_character_,
          resolved_name = NA_character_,
          resolved_namespace = NA_character_,
          resolution_status = "no_gprofiler_conversion"
        )
    )
  }
  
  conv <- tibble::as_tibble(conv)
  
  query_col <- intersect(c("input", "incoming", "query"), colnames(conv))[1]
  target_col <- intersect(c("target", "converted", "ensg"), colnames(conv))[1]
  name_col <- intersect(c("name", "names", "description"), colnames(conv))[1]
  ns_col <- intersect(c("namespace", "source"), colnames(conv))[1]
  
  if (is.na(query_col) || is.na(target_col)) {
    
    return(
      candidates %>%
        dplyr::group_by(feature_model) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          resolved_id = NA_character_,
          resolved_name = NA_character_,
          resolved_namespace = NA_character_,
          resolution_status = "unexpected_gconvert_columns"
        )
    )
  }
  
  conv2 <- conv %>%
    dplyr::transmute(
      query_candidate = as.character(.data[[query_col]]),
      resolved_id = as.character(.data[[target_col]]),
      resolved_name = if (!is.na(name_col)) {
        as.character(.data[[name_col]])
      } else {
        NA_character_
      },
      resolved_namespace = if (!is.na(ns_col)) {
        as.character(.data[[ns_col]])
      } else {
        NA_character_
      }
    ) %>%
    dplyr::filter(
      !is.na(resolved_id),
      resolved_id != ""
    ) %>%
    dplyr::distinct()
  
  candidates %>%
    dplyr::left_join(
      conv2,
      by = "query_candidate"
    ) %>%
    dplyr::mutate(
      has_resolution = !is.na(resolved_id) & resolved_id != ""
    ) %>%
    dplyr::arrange(
      feature_model,
      dplyr::desc(has_resolution),
      candidate_rank
    ) %>%
    dplyr::group_by(feature_model) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      resolution_status = dplyr::if_else(
        has_resolution,
        "resolved_by_gprofiler_gconvert",
        "not_resolved_by_gprofiler_gconvert"
      )
    ) %>%
    dplyr::select(
      feature_model,
      feature_label,
      view_code,
      selected_query = query_candidate,
      candidate_rank,
      resolved_id,
      resolved_name,
      resolved_namespace,
      resolution_status
    )
}

.ols_search_one_dynamic <- function(query,
                                    feature_model,
                                    feature_label,
                                    rows = 5) {
  
  empty_ols_result <- function(status) {
    tibble::tibble(
      feature_model = as.character(feature_model),
      feature_label = as.character(feature_label),
      clinical_query = as.character(query)[1],
      ontology_source = NA_character_,
      ontology_id = NA_character_,
      ontology_label = NA_character_,
      ontology_iri = NA_character_,
      ontology_description = NA_character_,
      ols_rank = NA_integer_,
      ols_status = status
    )
  }
  
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    return(empty_ols_result("jsonlite_not_installed"))
  }
  
  query <- as.character(query)
  query <- query[1]
  query <- trimws(query)
  
  if (length(query) == 0 || is.na(query) || !nzchar(query)) {
    return(empty_ols_result("empty_query"))
  }
  
  url <- paste0(
    "https://www.ebi.ac.uk/ols4/api/search?q=",
    utils::URLencode(query, reserved = TRUE),
    "&rows=",
    rows * 10
  )
  
  obj <- tryCatch(
    jsonlite::fromJSON(url),
    error = function(e) NULL
  )
  
  docs_raw <- tryCatch(
    obj$response$docs,
    error = function(e) NULL
  )
  
  if (is.null(docs_raw)) {
    return(empty_ols_result("no_ols_response_docs"))
  }
  
  docs <- tryCatch(
    tibble::as_tibble(docs_raw),
    error = function(e) tibble::tibble()
  )
  
  if (!is.data.frame(docs) || nrow(docs) < 1) {
    return(empty_ols_result("no_ols_hit"))
  }
  
  preferred_ontologies <- c(
    "efo",
    "hp",
    "mesh",
    "ncit",
    "mondo",
    "oba",
    "pato"
  )
  
  if ("ontology_name" %in% colnames(docs)) {
    docs <- docs %>%
      dplyr::mutate(
        ontology_name_chr = as.character(ontology_name),
        ontology_priority = dplyr::case_when(
          tolower(ontology_name_chr) %in% preferred_ontologies ~
            match(tolower(ontology_name_chr), preferred_ontologies),
          TRUE ~ 999L
        )
      ) %>%
      dplyr::arrange(ontology_priority)
  }
  
  get_col <- function(df, nm, default = NA_character_) {
    
    if (!nm %in% colnames(df)) {
      return(rep(default, nrow(df)))
    }
    
    out <- df[[nm]]
    
    if (is.list(out)) {
      out <- vapply(
        out,
        function(z) {
          z <- unlist(z)
          z <- z[!is.na(z) & z != ""]
          if (length(z) == 0) return(NA_character_)
          as.character(z[1])
        },
        character(1)
      )
    } else {
      out <- as.character(out)
    }
    
    out
  }
  
  get_description <- function(df) {
    
    if (!"description" %in% colnames(df)) {
      return(rep(NA_character_, nrow(df)))
    }
    
    out <- df$description
    
    if (is.list(out)) {
      return(
        vapply(
          out,
          function(z) {
            z <- unlist(z)
            z <- z[!is.na(z) & z != ""]
            if (length(z) == 0) return(NA_character_)
            as.character(z[1])
          },
          character(1)
        )
      )
    }
    
    as.character(out)
  }
  
  ontology_id_vec <- dplyr::coalesce(
    get_col(docs, "obo_id"),
    get_col(docs, "short_form"),
    get_col(docs, "iri")
  )
  
  out <- tibble::tibble(
    feature_model = as.character(feature_model),
    feature_label = as.character(feature_label),
    clinical_query = query,
    ontology_source = get_col(docs, "ontology_name"),
    ontology_id = ontology_id_vec,
    ontology_label = get_col(docs, "label"),
    ontology_iri = get_col(docs, "iri"),
    ontology_description = get_description(docs),
    ols_rank = seq_len(nrow(docs)),
    ols_status = "ols_hit"
  ) %>%
    dplyr::slice_head(n = rows)
  
  out
}

.run_clinical_ols_annotation_dynamic <- function(feature_annotation_base,
                                                 rows = 5) {
  
  cl_tbl <- feature_annotation_base %>%
    dplyr::filter(view_code == "cl") %>%
    dplyr::mutate(
      clinical_query = dplyr::case_when(
        !is.na(clinical_query) & clinical_query != "" ~ clinical_query,
        TRUE ~ query_text_clean
      ),
      clinical_query = as.character(clinical_query)
    )
  
  if (nrow(cl_tbl) == 0) {
    return(tibble::tibble())
  }
  
  get_optional_chr <- function(df, col, default = NA_character_) {
    if (col %in% colnames(df)) {
      out <- as.character(df[[col]])
      out[is.na(out)] <- default
      return(out)
    }
    rep(default, nrow(df))
  }
  
  ## ------------------------------------------------------------
  ## 1) Clínicas curadas que NO deben pasar por OLS.
  ## Ejemplos:
  ##   FLI      -> fatty liver index
  ##   CUN_BAE  -> body adiposity estimator
  ##   HOMA-IR  -> insulin resistance index
  ##   QUICKI   -> insulin sensitivity index
  ## ------------------------------------------------------------
  
  curated_no_ols <- cl_tbl %>%
    dplyr::mutate(
      manual_ontology_source_chr = get_optional_chr(., "manual_ontology_source", ""),
      manual_ontology_id_chr = get_optional_chr(., "manual_ontology_id", ""),
      manual_ontology_label_chr = get_optional_chr(., "manual_ontology_label", ""),
      manual_clinical_notes_chr = get_optional_chr(., "manual_clinical_notes", "")
    ) %>%
    dplyr::filter(
      manual_ontology_source_chr == "curated_clinical_concept"
    )
  
  curated_no_ols_out <- if (nrow(curated_no_ols) > 0) {
    
    curated_no_ols %>%
      dplyr::transmute(
        feature_model = as.character(feature_model),
        feature_label = as.character(feature_label),
        clinical_query = as.character(clinical_query),
        ontology_source = manual_ontology_source_chr,
        ontology_id = manual_ontology_id_chr,
        ontology_label = dplyr::case_when(
          !is.na(manual_ontology_label_chr) &
            manual_ontology_label_chr != "" ~ manual_ontology_label_chr,
          TRUE ~ clinical_query
        ),
        ontology_iri = NA_character_,
        ontology_description = dplyr::case_when(
          !is.na(manual_clinical_notes_chr) &
            manual_clinical_notes_chr != "" ~ manual_clinical_notes_chr,
          TRUE ~ "Curated clinical concept from reference_maps; OLS was intentionally skipped."
        ),
        ols_rank = 1L,
        ols_status = "curated_master_no_ols"
      )
    
  } else {
    
    tibble::tibble()
  }
  
  ## ------------------------------------------------------------
  ## 2) Resto de clínicas:
  ##    usar OLS, pero ya con query canónica en inglés.
  ##    Ejemplo:
  ##      IMC  -> body mass index
  ##      Peso -> body weight
  ## ------------------------------------------------------------
  
  cl_tbl_ols <- cl_tbl %>%
    dplyr::mutate(
      manual_ontology_source_chr = get_optional_chr(., "manual_ontology_source", "")
    ) %>%
    dplyr::filter(
      manual_ontology_source_chr != "curated_clinical_concept"
    )
  
  ols_out <- if (nrow(cl_tbl_ols) > 0) {
    
    purrr::map_dfr(
      seq_len(nrow(cl_tbl_ols)),
      function(i) {
        
        tryCatch(
          {
            .ols_search_one_dynamic(
              query = cl_tbl_ols$clinical_query[i],
              feature_model = cl_tbl_ols$feature_model[i],
              feature_label = cl_tbl_ols$feature_label[i],
              rows = rows
            )
          },
          error = function(e) {
            tibble::tibble(
              feature_model = as.character(cl_tbl_ols$feature_model[i]),
              feature_label = as.character(cl_tbl_ols$feature_label[i]),
              clinical_query = as.character(cl_tbl_ols$clinical_query[i]),
              ontology_source = NA_character_,
              ontology_id = NA_character_,
              ontology_label = NA_character_,
              ontology_iri = NA_character_,
              ontology_description = NA_character_,
              ols_rank = NA_integer_,
              ols_status = paste0("ols_error: ", conditionMessage(e))
            )
          }
        )
      }
    )
    
  } else {
    
    tibble::tibble()
  }
  
  dplyr::bind_rows(
    curated_no_ols_out,
    ols_out
  )
}

.make_gprofiler_query_sets_dynamic <- function(feature_annotation_base,
                                               id_resolution,
                                               module_members = NULL) {
  
  resolved_features <- feature_annotation_base %>%
    dplyr::filter(view_code %in% c("tx", "pr")) %>%
    dplyr::left_join(
      id_resolution %>%
        dplyr::select(
          feature_model,
          resolved_id,
          selected_query,
          resolution_status
        ),
      by = "feature_model"
    ) %>%
    dplyr::mutate(
      enrichment_id = dplyr::case_when(
        !is.na(resolved_id) & resolved_id != "" ~ resolved_id,
        !is.na(manual_gene_symbol) & manual_gene_symbol != "" ~ manual_gene_symbol,
        !is.na(selected_query) & selected_query != "" ~ selected_query,
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(
      !is.na(enrichment_id),
      enrichment_id != ""
    )
  
  query_sets <- list()
  
  all_ids <- unique(resolved_features$enrichment_id)
  all_ids <- all_ids[!is.na(all_ids) & all_ids != ""]
  
  if (length(all_ids) >= 2) {
    query_sets[["all_tx_pr"]] <- all_ids
  }
  
  for (vc in c("tx", "pr")) {
    
    ids_v <- resolved_features %>%
      dplyr::filter(view_code == vc) %>%
      dplyr::pull(enrichment_id) %>%
      unique()
    
    ids_v <- ids_v[!is.na(ids_v) & ids_v != ""]
    
    if (length(ids_v) >= 2) {
      query_sets[[paste0("view_", vc)]] <- ids_v
    }
  }
  
  if (!is.null(module_members) && nrow(module_members) > 0) {
    
    module_features <- module_members %>%
      dplyr::left_join(
        resolved_features %>%
          dplyr::select(
            feature_model,
            enrichment_id
          ),
        by = c("node_id" = "feature_model")
      ) %>%
      dplyr::filter(
        !is.na(enrichment_id),
        enrichment_id != ""
      )
    
    module_sets <- module_features %>%
      dplyr::group_by(signal_label, module_id) %>%
      dplyr::summarise(
        ids = list(unique(enrichment_id)),
        n_ids = length(unique(enrichment_id)),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_ids >= 2)
    
    if (nrow(module_sets) > 0) {
      for (i in seq_len(nrow(module_sets))) {
        
        nm <- paste0(
          "module_",
          module_sets$signal_label[i],
          "_",
          module_sets$module_id[i]
        )
        
        query_sets[[nm]] <- module_sets$ids[[i]]
      }
    }
  }
  
  query_sets
}

.run_gprofiler_annotation_dynamic <- function(query_sets,
                                              organism = "hsapiens",
                                              user_threshold = 0.20) {
  
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    return(
      list(
        result = tibble::tibble(),
        status = "gprofiler2_not_installed"
      )
    )
  }
  
  query_sets <- query_sets[
    vapply(
      query_sets,
      function(z) length(unique(z[!is.na(z) & z != ""])) >= 2,
      logical(1)
    )
  ]
  
  if (length(query_sets) == 0) {
    return(
      list(
        result = tibble::tibble(),
        status = "no_query_set_with_at_least_2_ids"
      )
    )
  }
  
  gost_res <- tryCatch(
    gprofiler2::gost(
      query = query_sets,
      organism = organism,
      ordered_query = FALSE,
      multi_query = FALSE,
      significant = TRUE,
      user_threshold = user_threshold,
      correction_method = "fdr",
      sources = c(
        "GO:BP",
        "GO:MF",
        "GO:CC",
        "REAC",
        "KEGG",
        "WP"
      ),
      evcodes = TRUE
    ),
    error = function(e) {
      structure(
        list(error_message = conditionMessage(e)),
        class = "gprofiler_error"
      )
    }
  )
  
  if (inherits(gost_res, "gprofiler_error")) {
    return(
      list(
        result = tibble::tibble(),
        status = paste0("gprofiler_error: ", gost_res$error_message)
      )
    )
  }
  
  if (
    is.null(gost_res) ||
    is.null(gost_res$result) ||
    nrow(gost_res$result) == 0
  ) {
    return(
      list(
        result = tibble::tibble(),
        status = "no_significant_gprofiler_terms"
      )
    )
  }
  
  list(
    result = tibble::as_tibble(gost_res$result),
    status = "ok"
  )
}

.extract_gprofiler_feature_terms_dynamic <- function(gprofiler_terms,
                                                     feature_annotation_base,
                                                     id_resolution) {
  
  if (is.null(gprofiler_terms) || nrow(gprofiler_terms) == 0) {
    return(tibble::tibble())
  }
  
  gp <- tibble::as_tibble(gprofiler_terms)
  
  if (!"intersection" %in% colnames(gp)) {
    return(tibble::tibble())
  }
  
  resolved_features <- feature_annotation_base %>%
    dplyr::filter(view_code %in% c("tx", "pr")) %>%
    dplyr::left_join(
      id_resolution %>%
        dplyr::select(
          feature_model,
          resolved_id,
          selected_query
        ),
      by = "feature_model"
    ) %>%
    dplyr::mutate(
      enrichment_id = dplyr::case_when(
        !is.na(resolved_id) & resolved_id != "" ~ resolved_id,
        !is.na(manual_gene_symbol) & manual_gene_symbol != "" ~ manual_gene_symbol,
        !is.na(selected_query) & selected_query != "" ~ selected_query,
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(
      !is.na(enrichment_id),
      enrichment_id != ""
    )
  
  gp <- gp %>%
    dplyr::mutate(
      intersection_chr = vapply(
        intersection,
        function(z) {
          z <- unlist(z)
          z <- z[!is.na(z) & z != ""]
          paste(z, collapse = ",")
        },
        character(1)
      )
    )
  
  gp_long <- gp %>%
    dplyr::select(
      dplyr::any_of(c(
        "query",
        "source",
        "term_id",
        "term_name",
        "p_value",
        "intersection_size",
        "query_size",
        "term_size",
        "intersection_chr"
      ))
    ) %>%
    tidyr::separate_rows(
      intersection_chr,
      sep = ","
    ) %>%
    dplyr::rename(enrichment_id = intersection_chr) %>%
    dplyr::filter(
      !is.na(enrichment_id),
      enrichment_id != ""
    )
  
  gp_long %>%
    dplyr::left_join(
      resolved_features %>%
        dplyr::select(
          feature_model,
          feature_label,
          view_code,
          enrichment_id
        ),
      by = "enrichment_id"
    ) %>%
    dplyr::filter(!is.na(feature_model)) %>%
    dplyr::arrange(
      feature_model,
      p_value
    )
}

.read_module_members_for_annotation_dynamic <- function(scenario_dir,
                                                        diag_target) {
  
  module_paths <- tibble::tibble(
    signal_label = c("reconstructed", "observed"),
    module_file = file.path(
      scenario_dir,
      diag_target,
      c(
        "framework_bidireccional_reconstruido",
        "framework_bidireccional_observado"
      ),
      "tables",
      "10_feature_feature_neutral_module_members.csv"
    )
  )
  
  purrr::map_dfr(
    seq_len(nrow(module_paths)),
    function(i) {
      
      mm <- .safe_read_csv_annotation(module_paths$module_file[i])
      
      if (is.null(mm)) {
        return(tibble::tibble())
      }
      
      tibble::as_tibble(mm) %>%
        dplyr::mutate(
          signal_label = module_paths$signal_label[i],
          module_file = module_paths$module_file[i]
        )
    }
  )
}

.summarise_ranking_support_for_annotation_dynamic <- function(ranking_observed,
                                                              ranking_reconstructed) {
  
  ranking_all <- dplyr::bind_rows(
    ranking_observed,
    ranking_reconstructed
  )
  
  if (nrow(ranking_all) == 0) {
    return(tibble::tibble())
  }
  
  ranking_all <- ranking_all %>%
    dplyr::mutate(
      signal_mode = as.character(signal_mode),
      signal_mode = dplyr::case_when(
        signal_mode %in% c("observed", "observado", "raw", "original", "sin_integrar") ~
          "observed",
        signal_mode %in% c("reconstructed", "reconstruido", "integrado", "mofa_reconstructed") ~
          "reconstructed",
        TRUE ~ signal_mode
      )
    )
  
  ranking_wide <- ranking_all %>%
    dplyr::group_by(feature_model, signal_mode) %>%
    dplyr::slice_max(
      order_by = abs_rank_score,
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      feature_model,
      signal_mode,
      best_contrast = contrast,
      best_direction = direction,
      best_rank_score = rank_score,
      best_abs_rank_score = abs_rank_score,
      best_p_value = p_value,
      best_q_value = q_value
    ) %>%
    tidyr::pivot_wider(
      names_from = signal_mode,
      values_from = c(
        best_contrast,
        best_direction,
        best_rank_score,
        best_abs_rank_score,
        best_p_value,
        best_q_value
      ),
      names_sep = "_"
    )
  
  required_cols <- c(
    "best_contrast_observed",
    "best_contrast_reconstructed",
    "best_direction_observed",
    "best_direction_reconstructed",
    "best_rank_score_observed",
    "best_rank_score_reconstructed",
    "best_abs_rank_score_observed",
    "best_abs_rank_score_reconstructed",
    "best_p_value_observed",
    "best_p_value_reconstructed",
    "best_q_value_observed",
    "best_q_value_reconstructed"
  )
  
  for (cc in required_cols) {
    if (!cc %in% colnames(ranking_wide)) {
      ranking_wide[[cc]] <- NA
    }
  }
  
  ranking_wide %>%
    dplyr::mutate(
      best_rank_score_observed = as.numeric(best_rank_score_observed),
      best_rank_score_reconstructed = as.numeric(best_rank_score_reconstructed),
      best_abs_rank_score_observed = as.numeric(best_abs_rank_score_observed),
      best_abs_rank_score_reconstructed = as.numeric(best_abs_rank_score_reconstructed),
      best_p_value_observed = as.numeric(best_p_value_observed),
      best_p_value_reconstructed = as.numeric(best_p_value_reconstructed),
      
      dual_direction_concordant = dplyr::case_when(
        is.finite(best_rank_score_observed) &
          is.finite(best_rank_score_reconstructed) ~
          sign(best_rank_score_observed) == sign(best_rank_score_reconstructed),
        TRUE ~ NA
      ),
      
      best_dual_abs_rank_score = dplyr::case_when(
        is.finite(best_abs_rank_score_observed) &
          is.finite(best_abs_rank_score_reconstructed) ~
          pmax(best_abs_rank_score_observed, best_abs_rank_score_reconstructed),
        
        is.finite(best_abs_rank_score_observed) ~
          best_abs_rank_score_observed,
        
        is.finite(best_abs_rank_score_reconstructed) ~
          best_abs_rank_score_reconstructed,
        
        TRUE ~ NA_real_
      ),
      
      best_dual_p_value = dplyr::case_when(
        is.finite(best_p_value_observed) &
          is.finite(best_p_value_reconstructed) ~
          pmin(best_p_value_observed, best_p_value_reconstructed),
        
        is.finite(best_p_value_observed) ~
          best_p_value_observed,
        
        is.finite(best_p_value_reconstructed) ~
          best_p_value_reconstructed,
        
        TRUE ~ NA_real_
      )
    )
}

.build_feature_database_annotation_dynamic <- function(feature_annotation_base,
                                                       id_resolution,
                                                       clinical_ols,
                                                       gprofiler_feature_terms,
                                                       ranking_support) {
  
  gp_feature_summary <- if (
    !is.null(gprofiler_feature_terms) &&
    nrow(gprofiler_feature_terms) > 0
  ) {
    
    gprofiler_feature_terms %>%
      dplyr::group_by(feature_model) %>%
      dplyr::summarise(
        n_database_terms = dplyr::n_distinct(term_id),
        top_database_terms = .collapse_unique_annotation(
          paste0(source, ":", term_name),
          max_n = 6
        ),
        top_term_ids = .collapse_unique_annotation(term_id, max_n = 6),
        min_gprofiler_p_value = .safe_min_annotation(p_value),
        gprofiler_sources = .collapse_unique_annotation(source, max_n = 6),
        .groups = "drop"
      )
    
  } else {
    
    tibble::tibble(
      feature_model = character(),
      n_database_terms = integer(),
      top_database_terms = character(),
      top_term_ids = character(),
      min_gprofiler_p_value = numeric(),
      gprofiler_sources = character()
    )
  }
  
  clinical_top <- if (
    !is.null(clinical_ols) &&
    nrow(clinical_ols) > 0
  ) {
    
    clinical_ols %>%
      dplyr::group_by(feature_model) %>%
      dplyr::arrange(ols_rank, .by_group = TRUE) %>%
      dplyr::summarise(
        clinical_query = dplyr::first(clinical_query),
        ontology_source = dplyr::first(ontology_source),
        ontology_id = dplyr::first(ontology_id),
        ontology_label = dplyr::first(ontology_label),
        ontology_iri = dplyr::first(ontology_iri),
        ols_status = dplyr::first(ols_status),
        .groups = "drop"
      )
    
  } else {
    
    tibble::tibble(
      feature_model = character(),
      clinical_query = character(),
      ontology_source = character(),
      ontology_id = character(),
      ontology_label = character(),
      ontology_iri = character(),
      ols_status = character()
    )
  }
  
  feature_annotation_base %>%
    dplyr::left_join(
      id_resolution,
      by = c("feature_model", "feature_label", "view_code")
    ) %>%
    dplyr::left_join(
      gp_feature_summary,
      by = "feature_model"
    ) %>%
    dplyr::left_join(
      clinical_top,
      by = "feature_model",
      suffix = c("", "_ols")
    ) %>%
    dplyr::left_join(
      ranking_support,
      by = "feature_model"
    ) %>%
    dplyr::mutate(
      n_database_terms = dplyr::coalesce(n_database_terms, 0L),
      
      metabolite_database_ids = dplyr::case_when(
        view_code == "me" ~ .collapse_unique_annotation(
          c(
            if ("HMDB_ID" %in% colnames(.)) HMDB_ID else NA_character_,
            if ("KEGG_ID" %in% colnames(.)) KEGG_ID else NA_character_,
            if ("ChEBI_ID" %in% colnames(.)) ChEBI_ID else NA_character_,
            if ("PubChem_CID" %in% colnames(.)) PubChem_CID else NA_character_
          ),
          max_n = 10
        ),
        TRUE ~ NA_character_
      ),
      
      annotation_status = dplyr::case_when(
        view_code == "me" &
          !is.na(metabolite_database_ids) &
          metabolite_database_ids != "" ~
          "metabolite_id_available_no_pathway_enrichment",
        
        view_code == "me" ~
          "unidentified_metabolomic_signal",
        
        view_code == "cl" &
          !is.na(ontology_label) &
          ontology_label != "" ~
          "database_supported_clinical_ontology_mapping",
        
        view_code %in% c("tx", "pr") &
          n_database_terms > 0 ~
          "database_supported_functional_enrichment",
        
        view_code %in% c("tx", "pr") &
          !is.na(resolved_id) &
          resolved_id != "" ~
          "identifier_resolved_no_enriched_term",
        
        view_code %in% c("tx", "pr") ~
          "identifier_not_resolved",
        
        view_code == "cl" ~
          "clinical_ontology_mapping_not_found",
        
        TRUE ~
          "unmapped"
      ),
      
      database_supported_annotation = dplyr::case_when(
        view_code == "me" &
          !is.na(metabolite_database_ids) &
          metabolite_database_ids != "" ~
          metabolite_database_ids,
        
        view_code == "me" ~
          NA_character_,
        
        view_code == "cl" &
          !is.na(ontology_label) &
          ontology_label != "" ~
          paste0(
            ontology_source,
            ":",
            ontology_id,
            " - ",
            ontology_label
          ),
        
        view_code %in% c("tx", "pr") &
          !is.na(top_database_terms) &
          top_database_terms != "" ~
          top_database_terms,
        
        view_code %in% c("tx", "pr") &
          !is.na(resolved_id) &
          resolved_id != "" ~
          paste0("Resolved identifier: ", resolved_id),
        
        TRUE ~
          NA_character_
      ),
      
      evidence_type = dplyr::case_when(
        view_code == "tx" & n_database_terms > 0 ~
          "gene_identifier_resolved_then_gprofiler_enrichment",
        
        view_code == "pr" & n_database_terms > 0 ~
          "protein_label_or_id_resolved_then_gprofiler_enrichment",
        
        view_code %in% c("tx", "pr") &
          !is.na(resolved_id) &
          resolved_id != "" ~
          "identifier_resolution_only",
        
        view_code == "cl" &
          !is.na(ontology_label) &
          ontology_label != "" ~
          "ols_ontology_search_top_hit",
        
        view_code == "me" &
          !is.na(metabolite_database_ids) &
          metabolite_database_ids != "" ~
          "manual_metabolite_identifier_available",
        
        view_code == "me" ~
          "metabolite_identity_missing_no_pathway_annotation",
        
        TRUE ~
          "no_database_hit"
      ),
      
      phenotype_axis_interpretation = dplyr::case_when(
        view_code == "me" &
          annotation_status == "unidentified_metabolomic_signal" ~
          "Metabolomic signal retained statistically but not biologically annotated because metabolite identity is unavailable.",
        
        view_code == "me" ~
          "Metabolite identifier available, but pathway annotation is not attempted unless a curated metabolite-to-pathway map is supplied.",
        
        view_code == "cl" ~
          "Clinical variable mapped through ontology search; this supports semantic interpretation, not pathway enrichment.",
        
        view_code %in% c("tx", "pr") &
          n_database_terms > 0 ~
          "Database-supported gene/protein functional signal.",
        
        view_code %in% c("tx", "pr") &
          !is.na(resolved_id) &
          resolved_id != "" ~
          "Identifier resolved, but no significant enriched database term was detected.",
        
        view_code %in% c("tx", "pr") ~
          "Identifier could not be resolved automatically; use optional manual mapping template if needed.",
        
        TRUE ~
          "No interpretation assigned."
      )
    ) %>%
    dplyr::arrange(
      factor(
        annotation_status,
        levels = c(
          "database_supported_functional_enrichment",
          "database_supported_clinical_ontology_mapping",
          "identifier_resolved_no_enriched_term",
          "identifier_not_resolved",
          "clinical_ontology_mapping_not_found",
          "metabolite_id_available_no_pathway_enrichment",
          "unidentified_metabolomic_signal",
          "unmapped"
        )
      ),
      view_code,
      feature_model
    )
}

.build_module_database_annotation_dynamic <- function(module_members,
                                                      feature_database_annotation,
                                                      gprofiler_feature_terms) {
  
  if (is.null(module_members) || nrow(module_members) == 0) {
    return(
      list(
        module_members_annotated = tibble::tibble(),
        module_summary = tibble::tibble()
      )
    )
  }
  
  module_members_annotated <- module_members %>%
    dplyr::left_join(
      feature_database_annotation %>%
        dplyr::select(
          feature_model,
          feature_label,
          view_code_annotation = view_code,
          database_supported_annotation,
          annotation_status,
          evidence_type,
          phenotype_axis_interpretation,
          n_database_terms,
          top_database_terms,
          resolved_id,
          selected_query,
          clinical_query,
          ontology_label
        ),
      by = c("node_id" = "feature_model")
    )
  
  module_term_summary <- if (
    !is.null(gprofiler_feature_terms) &&
    nrow(gprofiler_feature_terms) > 0
  ) {
    
    module_members %>%
      dplyr::select(signal_label, module_id, node_id) %>%
      dplyr::left_join(
        gprofiler_feature_terms,
        by = c("node_id" = "feature_model")
      ) %>%
      dplyr::filter(!is.na(term_id)) %>%
      dplyr::group_by(signal_label, module_id) %>%
      dplyr::summarise(
        module_top_gene_protein_terms = .collapse_unique_annotation(
          paste0(source, ":", term_name),
          max_n = 8
        ),
        module_top_term_ids = .collapse_unique_annotation(term_id, max_n = 8),
        module_min_gprofiler_p_value = .safe_min_annotation(p_value),
        .groups = "drop"
      )
    
  } else {
    
    tibble::tibble(
      signal_label = character(),
      module_id = character(),
      module_top_gene_protein_terms = character(),
      module_top_term_ids = character(),
      module_min_gprofiler_p_value = numeric()
    )
  }
  
  module_summary <- module_members_annotated %>%
    dplyr::group_by(signal_label, module_id) %>%
    dplyr::summarise(
      n_features = dplyr::n(),
      views_present = .collapse_unique_annotation(view_code),
      annotation_statuses = .collapse_unique_annotation(annotation_status),
      n_database_supported = sum(
        annotation_status %in% c(
          "database_supported_functional_enrichment",
          "database_supported_clinical_ontology_mapping"
        ),
        na.rm = TRUE
      ),
      n_identifier_resolved = sum(
        annotation_status %in% c(
          "database_supported_functional_enrichment",
          "identifier_resolved_no_enriched_term"
        ),
        na.rm = TRUE
      ),
      n_unidentified_metabolites = sum(
        annotation_status == "unidentified_metabolomic_signal",
        na.rm = TRUE
      ),
      dominant_phenotype_module = .collapse_unique_annotation(
        phenotype_module,
        max_n = 3
      ),
      top_feature_annotations = .collapse_unique_annotation(
        database_supported_annotation,
        max_n = 8
      ),
      module_interpretation = dplyr::case_when(
        n_database_supported > 0 &
          n_unidentified_metabolites > 0 ~
          "Mixed module: database-supported tx/pr/cl features plus unidentified metabolomic signals.",
        
        n_database_supported > 0 ~
          "Module with database-supported annotation.",
        
        n_identifier_resolved > 0 ~
          "Module with resolved tx/pr identifiers but no significant functional enrichment.",
        
        n_unidentified_metabolites > 0 ~
          "Module contains unidentified metabolomic signals; no metabolite pathway annotation attempted.",
        
        TRUE ~
          "Module without sufficient database-supported annotation."
      ),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      module_term_summary,
      by = c("signal_label", "module_id")
    ) %>%
    dplyr::arrange(
      signal_label,
      dplyr::desc(n_database_supported),
      module_id
    )
  
  list(
    module_members_annotated = module_members_annotated,
    module_summary = module_summary
  )
}

.plot_database_annotation_audit_dynamic <- function(feature_database_annotation,
                                                    module_summary,
                                                    ranking_support,
                                                    outdir) {
  
  plots_dir <- file.path(outdir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  p_status <- feature_database_annotation %>%
    dplyr::count(
      view_code,
      annotation_status,
      name = "n"
    ) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = view_code,
        y = n,
        fill = annotation_status
      )
    ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = "Estado de anotación biológica por vista",
      subtitle = "Proceso dinámico: no usa nombres de features hardcodeados.",
      x = "Vista",
      y = "Número de features",
      fill = "Estado"
    ) +
    ggplot2::theme_minimal()
  
  .save_plot_publication(
    filename_base = file.path(plots_dir, "01_annotation_status_by_view"),
    plot = p_status,
    width = 10,
    height = 6,
    dpi = 600
  )
  
  p_feature_terms <- feature_database_annotation %>%
    dplyr::mutate(
      n_database_terms_plot = dplyr::coalesce(
        as.numeric(n_database_terms),
        0
      ),
      feature_plot = paste0(feature_label, " [", view_code, "]"),
      feature_plot = factor(
        feature_plot,
        levels = rev(feature_plot[order(view_code, n_database_terms_plot)])
      )
    ) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = n_database_terms_plot,
        y = feature_plot
      )
    ) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(
      view_code ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::labs(
      title = "Número de términos funcionales por feature",
      subtitle = "tx/pr: g:Profiler; cl: OLS; me: solo si hay IDs químicos externos.",
      x = "Número de términos de base de datos",
      y = "Feature"
    ) +
    ggplot2::theme_minimal()
  
  .save_plot_publication(
    filename_base = file.path(plots_dir, "02_database_term_count_by_feature"),
    plot = p_feature_terms,
    width = 10,
    height = max(7, 0.32 * nrow(feature_database_annotation)),
    dpi = 600
  )
  
  if (!is.null(module_summary) && nrow(module_summary) > 0) {
    
    p_module <- module_summary %>%
      dplyr::mutate(
        module_plot = paste(signal_label, module_id, sep = " | ")
      ) %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = module_plot,
          y = n_database_supported
        )
      ) +
      ggplot2::geom_col() +
      ggplot2::geom_text(
        ggplot2::aes(
          label = paste0(
            "DB=", n_database_supported,
            "\nresolved=", n_identifier_resolved,
            "\nME_unid=", n_unidentified_metabolites
          )
        ),
        vjust = -0.2,
        size = 3
      ) +
      ggplot2::labs(
        title = "Soporte de base de datos por módulo",
        subtitle = "DB = anotación soportada; resolved = ID tx/pr resuelto; ME_unid = metabolitos sin identidad.",
        x = "Flujo | módulo",
        y = "Features con soporte de base de datos"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
    .save_plot_publication(
      filename_base = file.path(plots_dir, "03_database_support_by_module"),
      plot = p_module,
      width = 11,
      height = 6,
      dpi = 600
    )
  }
  
  if (!is.null(ranking_support) && nrow(ranking_support) > 0) {
    
    ranking_plot_tbl <- feature_database_annotation %>%
      dplyr::left_join(
        ranking_support,
        by = "feature_model"
      )
    
    needed_rank_cols <- c(
      "best_rank_score_observed",
      "best_rank_score_reconstructed",
      "best_abs_rank_score_observed",
      "best_abs_rank_score_reconstructed"
    )
    
    for (cc in needed_rank_cols) {
      if (!cc %in% colnames(ranking_plot_tbl)) {
        ranking_plot_tbl[[cc]] <- NA_real_
      }
    }
    
    ranking_plot_tbl <- ranking_plot_tbl %>%
      dplyr::mutate(
        best_rank_score_observed = as.numeric(best_rank_score_observed),
        best_rank_score_reconstructed = as.numeric(best_rank_score_reconstructed),
        n_database_terms_plot = dplyr::coalesce(
          as.numeric(n_database_terms),
          0
        ),
        feature_plot = paste0(feature_label, " [", view_code, "]")
      ) %>%
      dplyr::filter(
        is.finite(best_rank_score_observed),
        is.finite(best_rank_score_reconstructed)
      )
    
    if (nrow(ranking_plot_tbl) > 0) {
      
      p_rank_annot <- ranking_plot_tbl %>%
        ggplot2::ggplot(
          ggplot2::aes(
            x = best_rank_score_observed,
            y = best_rank_score_reconstructed
          )
        ) +
        ggplot2::geom_hline(
          yintercept = 0,
          linetype = "dashed"
        ) +
        ggplot2::geom_vline(
          xintercept = 0,
          linetype = "dashed"
        ) +
        ggplot2::geom_abline(
          slope = 1,
          intercept = 0,
          linetype = "dotted"
        ) +
        ggplot2::geom_point(
          ggplot2::aes(
            shape = view_code,
            size = n_database_terms_plot
          ),
          alpha = 0.85
        ) +
        ggplot2::geom_text(
          ggplot2::aes(label = feature_plot),
          size = 2.7,
          vjust = -0.7,
          check_overlap = TRUE
        ) +
        ggplot2::labs(
          title = "Ranking dual con anotación biológica dinámica",
          subtitle = "Observed vs reconstructed. Tamaño = número de términos funcionales tx/pr.",
          x = "Mejor rank score observado",
          y = "Mejor rank score reconstruido",
          shape = "Vista",
          size = "N términos DB"
        ) +
        ggplot2::theme_minimal()
      
      .save_plot_publication(
        filename_base = file.path(
          plots_dir,
          "04_dual_ranking_with_database_annotation"
        ),
        plot = p_rank_annot,
        width = 11,
        height = 8,
        dpi = 600
      )
      
    } else {
      
      cat("\n[BIO-PLOT-04] No se generó plot 04 porque no hay columnas/rank scores observed y reconstructed finitos.\n")
      cat("[BIO-PLOT-04] La anotación queda guardada igual.\n")
    }
  }
  
  
  invisible(plots_dir)
}

run_database_biological_annotation_dynamic <- function(scenario_dir,
                                                       diag_target = "diagnostico_bayes_features",
                                                       outdir = NULL,
                                                       organism = "hsapiens",
                                                       gprofiler_threshold = 0.20) {
  
  scenario_dir <- normalizePath(path.expand(scenario_dir), mustWork = TRUE)
  
  ranking_dir <- file.path(
    scenario_dir,
    diag_target,
    "mofa_dual_feature_differential_ranking"
  )
  
  if (is.null(outdir)) {
    outdir <- file.path(
      scenario_dir,
      diag_target,
      "database_biological_annotation"
    )
  }
  
  tables_dir <- file.path(outdir, "tables")
  plots_dir  <- file.path(outdir, "plots")
  manual_dir <- file.path(outdir, "manual_maps")
  
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(manual_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\n============================================================\n")
  cat("ANOTACIÓN BIOLÓGICA DINÁMICA CON BASES DE DATOS\n")
  cat("No hardcodea features concretas.\n")
  cat("Outdir:\n")
  cat(outdir, "\n")
  cat("============================================================\n")
  
  feature_map_file <- file.path(
    ranking_dir,
    "24_feature_map_for_annotation.csv"
  )
  
  observed_file <- file.path(
    ranking_dir,
    "20_observed_all_feature_differential_ranking.csv"
  )
  
  reconstructed_file <- file.path(
    ranking_dir,
    "21_reconstructed_all_feature_differential_ranking.csv"
  )
  
  feature_map <- .safe_read_csv_annotation(feature_map_file)
  ranking_observed <- .safe_read_csv_annotation(observed_file)
  ranking_reconstructed <- .safe_read_csv_annotation(reconstructed_file)
  
  if (is.null(feature_map)) {
    stop("No existe: ", feature_map_file)
  }
  
  if (is.null(ranking_observed)) {
    stop("No existe: ", observed_file)
  }
  
  if (is.null(ranking_reconstructed)) {
    stop("No existe: ", reconstructed_file)
  }
  
  module_members <- .read_module_members_for_annotation_dynamic(
    scenario_dir = scenario_dir,
    diag_target = diag_target
  )
  
  cat("\n[BIO-01] Construyendo feature_annotation_base...\n")
  
  feature_annotation_base <- .make_feature_annotation_base_dynamic(
    feature_map = feature_map,
    manual_dir = manual_dir
  )
  
  cat("[BIO-01] OK:", nrow(feature_annotation_base), "features\n")
  
  cat("\n[BIO-02] Resolviendo IDs tx/pr con gprofiler2::gconvert...\n")
  
  id_resolution <- .resolve_gene_protein_identifiers_gprofiler(
    feature_annotation_base = feature_annotation_base,
    organism = organism
  )
  
  cat("[BIO-02] OK:", nrow(id_resolution), "filas\n")
  
  cat("\n[BIO-03] Resumiendo soporte de ranking dual...\n")
  
  ranking_support <- .summarise_ranking_support_for_annotation_dynamic(
    ranking_observed = ranking_observed,
    ranking_reconstructed = ranking_reconstructed
  )
  
  cat("[BIO-03] OK:", nrow(ranking_support), "features\n")
  
  cat("\n[BIO-04] Ejecutando OLS para clínicas...\n")
  
  clinical_ols <- .run_clinical_ols_annotation_dynamic(
    feature_annotation_base = feature_annotation_base,
    rows = 5
  )
  
  cat("[BIO-04] OK:", nrow(clinical_ols), "filas OLS\n")
  
  cat("\n[BIO-05] Construyendo query sets para enrichment...\n")
  
  query_sets <- .make_gprofiler_query_sets_dynamic(
    feature_annotation_base = feature_annotation_base,
    id_resolution = id_resolution,
    module_members = module_members
  )
  
  cat("[BIO-05] OK:", length(query_sets), "query sets\n")
  
  cat("\n[BIO-06] Ejecutando gprofiler2::gost...\n")
  
  gprofiler_obj <- .run_gprofiler_annotation_dynamic(
    query_sets = query_sets,
    organism = organism,
    user_threshold = gprofiler_threshold
  )
  
  gprofiler_terms <- gprofiler_obj$result
  
  cat("[BIO-06] Estado:", gprofiler_obj$status, "\n")
  cat("[BIO-06] Términos:", nrow(gprofiler_terms), "\n")
  
  gprofiler_feature_terms <- .extract_gprofiler_feature_terms_dynamic(
    gprofiler_terms = gprofiler_terms,
    feature_annotation_base = feature_annotation_base,
    id_resolution = id_resolution
  )
  
  feature_database_annotation <- .build_feature_database_annotation_dynamic(
    feature_annotation_base = feature_annotation_base,
    id_resolution = id_resolution,
    clinical_ols = clinical_ols,
    gprofiler_feature_terms = gprofiler_feature_terms,
    ranking_support = ranking_support
  )
  
  module_annotation <- .build_module_database_annotation_dynamic(
    module_members = module_members,
    feature_database_annotation = feature_database_annotation,
    gprofiler_feature_terms = gprofiler_feature_terms
  )
  
  module_members_annotated <- module_annotation$module_members_annotated
  module_summary <- module_annotation$module_summary
  
  annotation_run_manifest <- tibble::tibble(
    item = c(
      "n_features",
      "n_tx",
      "n_pr",
      "n_cl",
      "n_me",
      "n_id_resolved_tx_pr",
      "gprofiler_status",
      "n_gprofiler_terms",
      "n_gprofiler_feature_term_links",
      "n_clinical_ols_rows",
      "n_module_members",
      "n_module_summaries",
      "manual_maps_dir"
    ),
    value = c(
      nrow(feature_annotation_base),
      sum(feature_annotation_base$view_code == "tx", na.rm = TRUE),
      sum(feature_annotation_base$view_code == "pr", na.rm = TRUE),
      sum(feature_annotation_base$view_code == "cl", na.rm = TRUE),
      sum(feature_annotation_base$view_code == "me", na.rm = TRUE),
      sum(!is.na(id_resolution$resolved_id) & id_resolution$resolved_id != "", na.rm = TRUE),
      gprofiler_obj$status,
      nrow(gprofiler_terms),
      nrow(gprofiler_feature_terms),
      nrow(clinical_ols),
      nrow(module_members_annotated),
      nrow(module_summary),
      manual_dir
    )
  )
  
  .safe_write_csv_annotation(
    feature_annotation_base,
    file.path(tables_dir, "30_feature_annotation_base.csv")
  )
  
  .safe_write_csv_annotation(
    id_resolution,
    file.path(tables_dir, "31_tx_pr_identifier_resolution.csv")
  )
  
  .safe_write_csv_annotation(
    clinical_ols,
    file.path(tables_dir, "32_clinical_ols_mapping.csv")
  )
  
  .safe_write_csv_annotation(
    gprofiler_terms,
    file.path(tables_dir, "33_gprofiler_enrichment_terms.csv")
  )
  
  .safe_write_csv_annotation(
    gprofiler_feature_terms,
    file.path(tables_dir, "34_gprofiler_feature_term_links.csv")
  )
  
  .safe_write_csv_annotation(
    feature_database_annotation,
    file.path(tables_dir, "35_feature_database_annotation.csv")
  )
  
  .safe_write_csv_annotation(
    module_members_annotated,
    file.path(tables_dir, "36_module_members_database_annotation.csv")
  )
  
  .safe_write_csv_annotation(
    module_summary,
    file.path(tables_dir, "37_module_database_annotation_summary.csv")
  )
  
  .safe_write_csv_annotation(
    ranking_support,
    file.path(tables_dir, "38_feature_dual_ranking_support_summary.csv")
  )
  
  .safe_write_csv_annotation(
    annotation_run_manifest,
    file.path(tables_dir, "00_database_annotation_manifest.csv")
  )
  
  .plot_database_annotation_audit_dynamic(
    feature_database_annotation = feature_database_annotation,
    module_summary = module_summary,
    ranking_support = ranking_support,
    outdir = outdir
  )
  
  cat("\nAnotación biológica dinámica guardada:\n")
  cat(" - tables/30_feature_annotation_base.csv\n")
  cat(" - tables/31_tx_pr_identifier_resolution.csv\n")
  cat(" - tables/32_clinical_ols_mapping.csv\n")
  cat(" - tables/33_gprofiler_enrichment_terms.csv\n")
  cat(" - tables/34_gprofiler_feature_term_links.csv\n")
  cat(" - tables/35_feature_database_annotation.csv\n")
  cat(" - tables/36_module_members_database_annotation.csv\n")
  cat(" - tables/37_module_database_annotation_summary.csv\n")
  cat(" - tables/38_feature_dual_ranking_support_summary.csv\n")
  cat(" - manual_maps/*.csv\n")
  cat(" - plots/01_annotation_status_by_view.png/pdf\n")
  cat(" - plots/02_database_term_count_by_feature.png/pdf\n")
  cat(" - plots/03_database_support_by_module.png/pdf\n")
  cat(" - plots/04_dual_ranking_with_database_annotation.png/pdf\n")
  
  invisible(
    list(
      feature_annotation_base = feature_annotation_base,
      id_resolution = id_resolution,
      clinical_ols = clinical_ols,
      gprofiler_terms = gprofiler_terms,
      gprofiler_feature_terms = gprofiler_feature_terms,
      feature_database_annotation = feature_database_annotation,
      module_members_annotated = module_members_annotated,
      module_summary = module_summary,
      ranking_support = ranking_support,
      manifest = annotation_run_manifest,
      outdir = outdir,
      manual_dir = manual_dir
    )
  )
}



## ============================================================
## RED FINAL BIOLÓGICA ANOTADA POR MÓDULOS
## ------------------------------------------------------------
## Última figura interpretativa.
## Usa tablas ya calculadas:
##   framework_bidireccional_<observado/reconstruido>/tables/
##   database_biological_annotation/tables/
## No recalcula bidireccionalidad.
## ============================================================

plot_final_annotated_biological_module_network <- function(scenario_dir,
                                                           diag_target = "diagnostico_bayes_features",
                                                           signal_label = c("observado", "reconstruido"),
                                                           network_max_edges = 80,
                                                           min_module_size = 2,
                                                           top_terms_per_module = 4,
                                                           save_tables = TRUE) {
  
  keep_valid_edges_for_curve <- function(df) {
    
    df %>%
      dplyr::mutate(
        same_endpoint = is.finite(from_x) &
          is.finite(from_y) &
          is.finite(to_x) &
          is.finite(to_y) &
          abs(from_x - to_x) < 1e-10 &
          abs(from_y - to_y) < 1e-10
      ) %>%
      dplyr::filter(
        is.finite(from_x),
        is.finite(from_y),
        is.finite(to_x),
        is.finite(to_y),
        !same_endpoint
      ) %>%
      dplyr::select(-same_endpoint)
  }
  
  signal_label <- match.arg(signal_label)
  
  scenario_dir <- normalizePath(
    path.expand(scenario_dir),
    mustWork = TRUE
  )
  
  signal_key <- dplyr::case_when(
    signal_label == "observado" ~ "observed",
    signal_label == "reconstruido" ~ "reconstructed",
    TRUE ~ signal_label
  )
  
  signal_keys <- unique(c(signal_label, signal_key))
  
  bidir_outdir <- file.path(
    scenario_dir,
    diag_target,
    paste0("framework_bidireccional_", signal_label)
  )
  
  annotation_outdir <- file.path(
    scenario_dir,
    diag_target,
    "database_biological_annotation"
  )
  
  bidir_tables_dir <- file.path(bidir_outdir, "tables")
  annotation_tables_dir <- file.path(annotation_outdir, "tables")
  annotation_plots_dir <- file.path(annotation_outdir, "plots")
  
  dir.create(annotation_plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(annotation_tables_dir, recursive = TRUE, showWarnings = FALSE)
  
  read_tbl <- function(path, required = TRUE) {
    
    if (!file.exists(path)) {
      if (required) {
        stop("No existe archivo requerido:\n", path)
      }
      return(tibble::tibble())
    }
    
    tibble::as_tibble(
      read.csv(
        path,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    )
  }
  
  write_tbl <- function(x, path) {
    
    if (exists(".safe_write_csv_annotation", mode = "function")) {
      .safe_write_csv_annotation(x, path)
    } else {
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      write.csv(x, path, row.names = FALSE, na = "")
    }
    
    invisible(path)
  }
  
  get_chr_col <- function(df, col, default = NA_character_) {
    
    if (col %in% colnames(df)) {
      out <- as.character(df[[col]])
      out[is.na(out)] <- default
      return(out)
    }
    
    rep(default, nrow(df))
  }
  
  get_num_col <- function(df, col, default = NA_real_) {
    
    if (col %in% colnames(df)) {
      return(suppressWarnings(as.numeric(df[[col]])))
    }
    
    rep(default, nrow(df))
  }
  
  collapse_unique_local <- function(x, sep = " | ", max_n = 4) {
    
    x <- unique(as.character(x))
    x <- x[!is.na(x) & x != "" & x != "NA"]
    
    if (length(x) == 0) {
      return(NA_character_)
    }
    
    paste(head(x, max_n), collapse = sep)
  }
  
  mode_chr_local <- function(x) {
    
    x <- as.character(x)
    x <- x[!is.na(x) & x != ""]
    
    if (length(x) == 0) {
      return(NA_character_)
    }
    
    names(sort(table(x), decreasing = TRUE))[1]
  }
  
  trunc_chr <- function(x, n = 45) {
    
    x <- as.character(x)
    x[is.na(x)] <- ""
    
    ifelse(
      nchar(x) > n,
      paste0(substr(x, 1, n - 3), "..."),
      x
    )
  }
  
  phenotype_axis_from_pattern <- function(pattern, evidence_grade = NA_character_) {
    
    dplyr::case_when(
      pattern == "SO_like_high" ~ "SO-like",
      pattern == "FD_like_high" ~ "FD-like",
      pattern == "NP_like_high" ~ "NP-like",
      pattern == "SO_gt_FD_splitter" ~ "SO>FD splitter",
      pattern == "FD_gt_SO_splitter" ~ "FD>SO splitter",
      pattern == "shared_FD_SO_high" ~ "FD/SO shared",
      evidence_grade == "D_ambiguous_dependent" ~ "Ambiguous dependent",
      evidence_grade == "E_no_evidence" ~ "No evidence",
      TRUE ~ "Unclassified"
    )
  }
  
  filter_signal <- function(df) {
    
    if (!"signal_label" %in% colnames(df)) {
      return(df)
    }
    
    df %>%
      dplyr::filter(
        tolower(as.character(signal_label)) %in% tolower(signal_keys)
      )
  }
  
  ## ------------------------------------------------------------
  ## 1) Leer entradas
  ## ------------------------------------------------------------
  
  nodes_raw <- read_tbl(
    file.path(
      bidir_tables_dir,
      "07_feature_feature_hybrid_network_nodes.csv"
    ),
    required = TRUE
  )
  
  edges_raw <- read_tbl(
    file.path(
      bidir_tables_dir,
      "08_feature_feature_hybrid_network_edges.csv"
    ),
    required = TRUE
  )
  
  feature_signature <- read_tbl(
    file.path(
      bidir_tables_dir,
      "01_bidirectional_feature_signature.csv"
    ),
    required = TRUE
  )
  
  feature_annotation <- read_tbl(
    file.path(
      annotation_tables_dir,
      "35_feature_database_annotation.csv"
    ),
    required = TRUE
  )
  
  module_members_annotated <- read_tbl(
    file.path(
      annotation_tables_dir,
      "36_module_members_database_annotation.csv"
    ),
    required = FALSE
  ) %>%
    filter_signal()
  
  module_summary <- read_tbl(
    file.path(
      annotation_tables_dir,
      "37_module_database_annotation_summary.csv"
    ),
    required = FALSE
  ) %>%
    filter_signal()
  
  if (!all(c("node_id", "x", "y") %in% colnames(nodes_raw))) {
    stop("La tabla de nodos debe contener node_id, x, y.")
  }
  
  if (!all(c("from", "to") %in% colnames(edges_raw))) {
    stop("La tabla de aristas debe contener from y to.")
  }
  
  ## ------------------------------------------------------------
  ## 2) Preparar nodos evitando errores por columnas ausentes
  ## ------------------------------------------------------------
  
  feature_signature_compact <- feature_signature %>%
    dplyr::mutate(
      feature_model = as.character(feature_model)
    ) %>%
    dplyr::distinct(
      feature_model,
      .keep_all = TRUE
    )
  
  feature_annotation_compact <- feature_annotation %>%
    dplyr::mutate(
      feature_model = as.character(feature_model)
    ) %>%
    dplyr::distinct(
      feature_model,
      .keep_all = TRUE
    )
  
  sig_nodes <- tibble::tibble(
    node_id = as.character(feature_signature_compact$feature_model),
    
    sig_feature_label = dplyr::case_when(
      "feature_label" %in% colnames(feature_signature_compact) ~
        get_chr_col(feature_signature_compact, "feature_label"),
      TRUE ~ sub("^[^_]+_", "", as.character(feature_signature_compact$feature_model))
    ),
    
    sig_view_code = get_chr_col(feature_signature_compact, "view_code"),
    sig_evidence_grade = get_chr_col(feature_signature_compact, "evidence_grade"),
    sig_bidirectional_pattern = get_chr_col(feature_signature_compact, "bidirectional_pattern"),
    sig_posterior_high_group = get_chr_col(feature_signature_compact, "posterior_high_group"),
    sig_MI_norm = get_num_col(feature_signature_compact, "MI_norm"),
    sig_max_JS_bits = get_num_col(feature_signature_compact, "max_JS_bits"),
    sig_q_perm = get_num_col(feature_signature_compact, "q_perm"),
    sig_posterior_enrichment_high = get_num_col(feature_signature_compact, "posterior_enrichment_high")
  )
  
  ann_nodes <- tibble::tibble(
    node_id = as.character(feature_annotation_compact$feature_model),
    
    ann_database_supported_annotation = get_chr_col(
      feature_annotation_compact,
      "database_supported_annotation"
    ),
    
    ann_annotation_status = get_chr_col(
      feature_annotation_compact,
      "annotation_status"
    ),
    
    ann_evidence_type = get_chr_col(
      feature_annotation_compact,
      "evidence_type"
    ),
    
    ann_phenotype_axis_interpretation = get_chr_col(
      feature_annotation_compact,
      "phenotype_axis_interpretation"
    ),
    
    ann_n_database_terms = get_num_col(
      feature_annotation_compact,
      "n_database_terms"
    ),
    
    ann_top_database_terms = get_chr_col(
      feature_annotation_compact,
      "top_database_terms"
    ),
    
    ann_resolved_id = get_chr_col(
      feature_annotation_compact,
      "resolved_id"
    ),
    
    ann_selected_query = get_chr_col(
      feature_annotation_compact,
      "selected_query"
    ),
    
    ann_clinical_query = get_chr_col(
      feature_annotation_compact,
      "clinical_query"
    ),
    
    ann_ontology_label = get_chr_col(
      feature_annotation_compact,
      "ontology_label"
    )
  )
  
  if (
    nrow(module_members_annotated) > 0 &&
    all(c("node_id", "module_id") %in% colnames(module_members_annotated))
  ) {
    
    mod_nodes <- module_members_annotated %>%
      dplyr::mutate(
        node_id = as.character(node_id),
        mod_module_id = as.character(module_id),
        mod_phenotype_module = get_chr_col(., "phenotype_module")
      ) %>%
      dplyr::distinct(
        node_id,
        .keep_all = TRUE
      ) %>%
      dplyr::select(
        node_id,
        mod_module_id,
        mod_phenotype_module
      )
    
  } else if ("module_id" %in% colnames(nodes_raw)) {
    
    mod_nodes <- nodes_raw %>%
      dplyr::mutate(
        node_id = as.character(node_id),
        mod_module_id = as.character(module_id),
        mod_phenotype_module = get_chr_col(., "phenotype_module")
      ) %>%
      dplyr::distinct(
        node_id,
        .keep_all = TRUE
      ) %>%
      dplyr::select(
        node_id,
        mod_module_id,
        mod_phenotype_module
      )
    
  } else {
    
    mod_nodes <- tibble::tibble(
      node_id = as.character(nodes_raw$node_id),
      mod_module_id = "M01",
      mod_phenotype_module = NA_character_
    )
  }
  
  nodes <- nodes_raw %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      node_id = as.character(node_id),
      x = suppressWarnings(as.numeric(x)),
      y = suppressWarnings(as.numeric(y)),
      node_label_raw = get_chr_col(., "node_label")
    ) %>%
    dplyr::left_join(sig_nodes, by = "node_id") %>%
    dplyr::left_join(ann_nodes, by = "node_id") %>%
    dplyr::left_join(mod_nodes, by = "node_id") %>%
    dplyr::mutate(
      feature_model = node_id,
      
      feature_label = dplyr::case_when(
        !is.na(sig_feature_label) &
          sig_feature_label != "" ~ sig_feature_label,
        
        !is.na(node_label_raw) &
          node_label_raw != "" ~ sub("\\s*\\[.*$", "", node_label_raw),
        
        TRUE ~ sub("^[^_]+_", "", node_id)
      ),
      
      view_code = dplyr::case_when(
        !is.na(sig_view_code) &
          sig_view_code != "" ~ sig_view_code,
        
        TRUE ~ sub("_.*$", "", node_id)
      ),
      
      evidence_grade = dplyr::case_when(
        !is.na(sig_evidence_grade) &
          sig_evidence_grade != "" ~ sig_evidence_grade,
        
        TRUE ~ "unclassified"
      ),
      
      bidirectional_pattern = dplyr::case_when(
        !is.na(sig_bidirectional_pattern) &
          sig_bidirectional_pattern != "" ~ sig_bidirectional_pattern,
        
        TRUE ~ "unclassified"
      ),
      
      posterior_high_group = sig_posterior_high_group,
      posterior_enrichment_high = sig_posterior_enrichment_high,
      MI_norm = sig_MI_norm,
      max_JS_bits = sig_max_JS_bits,
      q_perm = sig_q_perm,
      
      phenotype_axis = phenotype_axis_from_pattern(
        bidirectional_pattern,
        evidence_grade
      ),
      
      annotation_label = dplyr::case_when(
        !is.na(ann_database_supported_annotation) &
          ann_database_supported_annotation != "" ~ ann_database_supported_annotation,
        
        !is.na(ann_top_database_terms) &
          ann_top_database_terms != "" ~ ann_top_database_terms,
        
        !is.na(ann_ontology_label) &
          ann_ontology_label != "" ~ ann_ontology_label,
        
        !is.na(ann_selected_query) &
          ann_selected_query != "" ~ ann_selected_query,
        
        !is.na(ann_clinical_query) &
          ann_clinical_query != "" ~ ann_clinical_query,
        
        TRUE ~ "no_database_annotation"
      ),
      
      annotation_label_short = trunc_chr(annotation_label, 38),
      feature_label_short = trunc_chr(feature_label, 24),
      
      node_plot_label = paste0(
        feature_label_short,
        "\n",
        annotation_label_short
      ),
      
      module_id = dplyr::case_when(
        !is.na(mod_module_id) &
          mod_module_id != "" ~ mod_module_id,
        
        TRUE ~ "M_unassigned"
      ),
      
      phenotype_module = mod_phenotype_module
    ) %>%
    dplyr::filter(
      is.finite(x),
      is.finite(y)
    )
  
  if (nrow(nodes) == 0) {
    stop("No hay nodos válidos para graficar.")
  }
  
  ## ------------------------------------------------------------
  ## 3) Aristas feature-feature
  ## ------------------------------------------------------------
  
  edges_feature_feature <- tibble::tibble(
    from = get_chr_col(edges_raw, "from"),
    to = get_chr_col(edges_raw, "to"),
    rho = get_num_col(edges_raw, "rho"),
    abs_rho = get_num_col(edges_raw, "abs_rho"),
    edge_axis_relation = get_chr_col(edges_raw, "edge_axis_relation")
  ) %>%
    dplyr::mutate(
      abs_rho = dplyr::case_when(
        is.finite(abs_rho) ~ abs_rho,
        is.finite(rho) ~ abs(rho),
        TRUE ~ NA_real_
      ),
      
      rho_sign = dplyr::case_when(
        is.finite(rho) & rho >= 0 ~ "positive",
        is.finite(rho) & rho < 0 ~ "negative",
        TRUE ~ "unknown"
      ),
      
      edge_relation = dplyr::case_when(
        !is.na(edge_axis_relation) &
          edge_axis_relation != "" ~ edge_axis_relation,
        
        TRUE ~ rho_sign
      )
    ) %>%
    dplyr::filter(
      from %in% nodes$node_id,
      to %in% nodes$node_id,
      is.finite(abs_rho)
    ) %>%
    dplyr::arrange(
      dplyr::desc(abs_rho)
    )
  
  if (nrow(edges_feature_feature) > network_max_edges) {
    edges_feature_feature <- edges_feature_feature %>%
      dplyr::slice_head(n = network_max_edges)
  }
  
  edges_feature_feature_plot <- edges_feature_feature %>%
    dplyr::left_join(
      nodes %>%
        dplyr::select(node_id, from_x = x, from_y = y),
      by = c("from" = "node_id")
    ) %>%
    dplyr::left_join(
      nodes %>%
        dplyr::select(node_id, to_x = x, to_y = y),
      by = c("to" = "node_id")
    ) %>%
    dplyr::mutate(
      edge_class = paste0("feature_feature_", rho_sign),
      edge_weight = dplyr::case_when(
        is.finite(abs_rho) ~ abs_rho,
        TRUE ~ 0.60
      )
    ) %>%
    keep_valid_edges_for_curve()
  ## ------------------------------------------------------------
  ## 4) Aristas feature -> grupo
  ## ------------------------------------------------------------
  
  group_levels <- c("NP", "FD", "SO")
  
  edges_feature_group_main <- nodes %>%
    dplyr::filter(
      !bidirectional_pattern %in% c(
        "shared_FD_SO_high",
        "SO_gt_FD_splitter",
        "FD_gt_SO_splitter"
      ),
      posterior_high_group %in% group_levels
    ) %>%
    dplyr::transmute(
      from = node_id,
      to_group = as.character(posterior_high_group),
      edge_class = "feature_to_posterior_group",
      edge_weight = pmax(posterior_enrichment_high, 0, na.rm = TRUE),
      edge_label = bidirectional_pattern
    )
  
  edges_feature_group_shared <- nodes %>%
    dplyr::filter(
      bidirectional_pattern == "shared_FD_SO_high"
    ) %>%
    tidyr::crossing(to_group = c("FD", "SO")) %>%
    dplyr::transmute(
      from = node_id,
      to_group = to_group,
      edge_class = "shared_FD_SO",
      edge_weight = pmax(posterior_enrichment_high, 0, na.rm = TRUE),
      edge_label = bidirectional_pattern
    )
  
  edges_feature_group_splitter <- nodes %>%
    dplyr::filter(
      bidirectional_pattern %in% c(
        "SO_gt_FD_splitter",
        "FD_gt_SO_splitter"
      )
    ) %>%
    tidyr::crossing(to_group = c("FD", "SO")) %>%
    dplyr::transmute(
      from = node_id,
      to_group = to_group,
      edge_class = "FD_SO_splitter",
      edge_weight = dplyr::case_when(
        is.finite(max_JS_bits) ~ max_JS_bits,
        TRUE ~ 0.50
      ),
      edge_label = bidirectional_pattern
    )
  
  edges_feature_group <- dplyr::bind_rows(
    edges_feature_group_main,
    edges_feature_group_shared,
    edges_feature_group_splitter
  ) %>%
    dplyr::mutate(
      edge_weight = suppressWarnings(as.numeric(edge_weight)),
      edge_weight = dplyr::case_when(
        is.finite(edge_weight) & edge_weight > 0 ~ edge_weight,
        TRUE ~ 0.25
      )
    )
  
  x_rng <- range(nodes$x, na.rm = TRUE)
  y_rng <- range(nodes$y, na.rm = TRUE)
  
  x_span <- diff(x_rng)
  y_span <- diff(y_rng)
  
  if (!is.finite(x_span) || x_span == 0) x_span <- 1
  if (!is.finite(y_span) || y_span == 0) y_span <- 1
  
  group_nodes <- tibble::tibble(
    group = group_levels,
    x = x_rng[1] - 0.55 * x_span,
    y = seq(
      from = y_rng[2],
      to = y_rng[1],
      length.out = length(group_levels)
    )
  )
  
  edges_feature_group_plot <- edges_feature_group %>%
    dplyr::left_join(
      nodes %>%
        dplyr::select(node_id, from_x = x, from_y = y),
      by = c("from" = "node_id")
    ) %>%
    dplyr::left_join(
      group_nodes %>%
        dplyr::select(to_group = group, to_x = x, to_y = y),
      by = "to_group"
    ) %>%
    keep_valid_edges_for_curve()
  
  ## ------------------------------------------------------------
  ## 5) Resumen de módulos
  ## ------------------------------------------------------------
  
  module_centers <- nodes %>%
    dplyr::group_by(module_id) %>%
    dplyr::summarise(
      n_features = dplyr::n(),
      x = median(x, na.rm = TRUE),
      y = max(y, na.rm = TRUE) + 0.10 * y_span,
      dominant_axis = mode_chr_local(phenotype_axis),
      dominant_view = mode_chr_local(view_code),
      module_features = collapse_unique_local(feature_label, max_n = 6),
      module_annotations_from_features = collapse_unique_local(
        annotation_label,
        max_n = top_terms_per_module
      ),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_features >= min_module_size)
  
  if (nrow(module_summary) > 0) {
    
    module_centers <- module_centers %>%
      dplyr::left_join(
        module_summary %>%
          dplyr::mutate(
            module_id = as.character(module_id)
          ) %>%
          dplyr::select(
            module_id,
            dplyr::any_of(c(
              "dominant_phenotype_module",
              "top_feature_annotations",
              "module_top_gene_protein_terms",
              "module_interpretation",
              "n_database_supported",
              "n_unidentified_metabolites"
            ))
          ) %>%
          dplyr::distinct(module_id, .keep_all = TRUE),
        by = "module_id"
      )
    
  } else {
    
    module_centers <- module_centers %>%
      dplyr::mutate(
        dominant_phenotype_module = NA_character_,
        top_feature_annotations = NA_character_,
        module_top_gene_protein_terms = NA_character_,
        module_interpretation = NA_character_,
        n_database_supported = NA_integer_,
        n_unidentified_metabolites = NA_integer_
      )
  }
  
  module_centers <- module_centers %>%
    dplyr::mutate(
      module_axis_label = dplyr::case_when(
        !is.na(dominant_phenotype_module) &
          dominant_phenotype_module != "" ~ as.character(dominant_phenotype_module),
        TRUE ~ dominant_axis
      ),
      
      module_annotation_label = dplyr::case_when(
        !is.na(top_feature_annotations) &
          top_feature_annotations != "" ~ as.character(top_feature_annotations),
        
        !is.na(module_top_gene_protein_terms) &
          module_top_gene_protein_terms != "" ~ as.character(module_top_gene_protein_terms),
        
        !is.na(module_annotations_from_features) &
          module_annotations_from_features != "" ~ as.character(module_annotations_from_features),
        
        TRUE ~ "sin anotación fuerte"
      ),
      
      module_plot_label = paste0(
        module_id,
        " | ",
        trunc_chr(module_axis_label, 34),
        "\n",
        trunc_chr(module_annotation_label, 70),
        "\n",
        "n=",
        n_features
      )
    )
  
  hull_tbl <- nodes %>%
    dplyr::filter(
      module_id %in% module_centers$module_id
    ) %>%
    dplyr::group_by(module_id) %>%
    dplyr::filter(dplyr::n() >= 3) %>%
    dplyr::slice(chull(x, y)) %>%
    dplyr::ungroup()
  
  ## ------------------------------------------------------------
  ## 6) Guardar tablas
  ## ------------------------------------------------------------
  
  if (save_tables) {
    
    write_tbl(
      nodes,
      file.path(
        annotation_tables_dir,
        paste0("39_final_annotated_network_nodes_", signal_label, ".csv")
      )
    )
    
    write_tbl(
      edges_feature_feature_plot,
      file.path(
        annotation_tables_dir,
        paste0("40_final_annotated_network_feature_feature_edges_", signal_label, ".csv")
      )
    )
    
    write_tbl(
      edges_feature_group_plot,
      file.path(
        annotation_tables_dir,
        paste0("41_final_annotated_network_feature_group_edges_", signal_label, ".csv")
      )
    )
    
    write_tbl(
      module_centers,
      file.path(
        annotation_tables_dir,
        paste0("42_final_annotated_network_modules_", signal_label, ".csv")
      )
    )
  }
  
  ## ------------------------------------------------------------
  ## 7) Plot final
  ## ------------------------------------------------------------
  
  p <- ggplot2::ggplot() +
    
    ggplot2::geom_polygon(
      data = hull_tbl,
      ggplot2::aes(
        x = x,
        y = y,
        group = module_id
      ),
      fill = "grey95",
      color = "grey55",
      linetype = "dashed",
      linewidth = 0.4,
      alpha = 0.45,
      show.legend = FALSE
    ) +
    
    ggplot2::geom_curve(
      data = edges_feature_feature_plot,
      ggplot2::aes(
        x = from_x,
        y = from_y,
        xend = to_x,
        yend = to_y,
        linewidth = edge_weight,
        linetype = edge_class
      ),
      curvature = 0.08,
      alpha = 0.35
    ) +
    
    ggplot2::geom_curve(
      data = edges_feature_group_plot,
      ggplot2::aes(
        x = from_x,
        y = from_y,
        xend = to_x,
        yend = to_y,
        linewidth = edge_weight,
        linetype = edge_class
      ),
      curvature = -0.18,
      alpha = 0.55,
      arrow = grid::arrow(
        length = grid::unit(0.10, "inches"),
        type = "closed"
      )
    ) +
    
    ggplot2::geom_point(
      data = group_nodes,
      ggplot2::aes(
        x = x,
        y = y
      ),
      shape = 22,
      size = 9,
      fill = "grey85",
      color = "black"
    ) +
    
    ggplot2::geom_text(
      data = group_nodes,
      ggplot2::aes(
        x = x,
        y = y,
        label = group
      ),
      fontface = "bold",
      size = 5
    ) +
    
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(
        x = x,
        y = y,
        fill = view_code,
        size = MI_norm,
        shape = phenotype_axis
      ),
      color = "black",
      alpha = 0.95
    ) +
    
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(
        x = x,
        y = y,
        label = feature_label_short
      ),
      size = 3,
      vjust = -0.9,
      check_overlap = TRUE
    ) +
    
    ggplot2::geom_label(
      data = module_centers,
      ggplot2::aes(
        x = x,
        y = y,
        label = module_plot_label
      ),
      size = 3.1,
      linewidth = 0.25,
      fill = "white",
      alpha = 0.92,
      fontface = "bold",
      show.legend = FALSE
    ) +
    
    ggplot2::scale_shape_manual(
      values = c(
        "NP-like" = 21,
        "FD-like" = 22,
        "SO-like" = 24,
        "SO>FD splitter" = 25,
        "FD>SO splitter" = 23,
        "FD/SO shared" = 21,
        "Ambiguous dependent" = 4,
        "No evidence" = 1,
        "Unclassified" = 1
      ),
      drop = FALSE
    ) +
    
    ggplot2::scale_linetype_manual(
      values = c(
        "feature_feature_positive" = "solid",
        "feature_feature_negative" = "dotted",
        "feature_feature_unknown" = "dashed",
        "feature_to_posterior_group" = "solid",
        "shared_FD_SO" = "longdash",
        "FD_SO_splitter" = "twodash"
      ),
      drop = FALSE
    ) +
    
    ggplot2::scale_linewidth_continuous(
      range = c(0.25, 1.8)
    ) +
    
    ggplot2::coord_equal(
      clip = "off"
    ) +
    
    ggplot2::labs(
      title = paste0(
        "Red biológica anotada de módulos — ",
        signal_label
      ),
      subtitle = paste0(
        "Features conectadas por correlación; flechas hacia NP/FD/SO según patrón bidireccional. ",
        "Las etiquetas de módulo resumen la anotación funcional/ontológica."
      ),
      x = NULL,
      y = NULL,
      fill = "Vista ómica",
      shape = "Patrón bidireccional",
      linewidth = "Fuerza",
      linetype = "Tipo de relación"
    ) +
    
    ggplot2::theme_void() +
    
    ggplot2::theme(
      legend.position = "right",
      plot.margin = ggplot2::margin(30, 120, 30, 120),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10)
    )
  
  filename_base <- file.path(
    annotation_plots_dir,
    paste0("05_final_annotated_biological_module_network_", signal_label)
  )
  
  .save_plot_publication(
    filename_base = filename_base,
    plot = p,
    width = 16,
    height = 12,
    dpi = 600
  )
  
  cat("\nRed biológica anotada guardada:\n")
  cat(" - ", paste0(filename_base, ".png"), "\n")
  cat(" - ", paste0(filename_base, ".pdf"), "\n")
  
  invisible(
    list(
      signal_label = signal_label,
      plot = p,
      nodes = nodes,
      feature_feature_edges = edges_feature_feature_plot,
      feature_group_edges = edges_feature_group_plot,
      modules = module_centers,
      plot_file_png = paste0(filename_base, ".png"),
      plot_file_pdf = paste0(filename_base, ".pdf")
    )
  )
}


## ============================================================
## HEATMAPS BIDIRECCIONALES ANOTADOS POR CONTEXTO
## ------------------------------------------------------------
## No reemplaza ni elimina la red anotada.
## Genera figuras por contexto:
##   NP, FD, SO, NP_vs_FD, NP_vs_SO, FD_vs_SO
##
## Paneles:
##   A) resumen estadístico: evidencia, Top->High, patrón, MI, Enrich
##   B) grupo -> feature : mu_NP / mu_FD / mu_SO
##   C) feature alta -> grupo : post_high_NP / post_high_FD / post_high_SO
##
## Nota:
##   La anotación biológica no se muestra en el plot.
##   Se conserva en las tablas de anotación para interpretación escrita.
## Usa:
##   - 01_bidirectional_feature_signature.csv
##   - 35_feature_database_annotation.csv
##   - 36_module_members_database_annotation.csv
##   - 37_module_database_annotation_summary.csv
##   - 42_final_annotated_network_modules_<signal>.csv   (si existe)
## ============================================================

plot_final_annotated_bidirectional_module_heatmap <- function(
    scenario_dir,
    diag_target = "diagnostico_bayes_features",
    signal_label = c("observado", "reconstruido"),
    y_levels = c("NP", "FD", "SO"),
    save_tables = TRUE
) {
  
  signal_label <- match.arg(signal_label)
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Necesitas instalar patchwork: install.packages('patchwork')")
  }
  
  scenario_dir <- normalizePath(path.expand(scenario_dir), mustWork = TRUE)
  
  bidir_tables_dir <- file.path(
    scenario_dir,
    diag_target,
    paste0("framework_bidireccional_", signal_label),
    "tables"
  )
  
  annotation_tables_dir <- file.path(
    scenario_dir,
    diag_target,
    "database_biological_annotation",
    "tables"
  )
  
  outdir <- file.path(
    scenario_dir,
    diag_target,
    paste0("final_annotated_bidirectional_module_heatmap_", signal_label)
  )
  
  plots_dir  <- file.path(outdir, "plots")
  tables_dir <- file.path(outdir, "tables")
  
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## ----------------------------------------------------------
  ## helpers
  ## ----------------------------------------------------------
  
  read_tbl_local <- function(path, required = TRUE) {
    if (!file.exists(path)) {
      if (required) stop("No existe archivo requerido:\n", path)
      return(tibble::tibble())
    }
    tibble::as_tibble(
      read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    )
  }
  
  write_tbl_local <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (exists(".safe_write_csv_annotation", mode = "function")) {
      .safe_write_csv_annotation(x, path)
    } else {
      write.csv(x, path, row.names = FALSE, na = "")
    }
    invisible(path)
  }
  
  save_plot_local <- function(filename_base, plot, width = 18, height = 10, dpi = 600) {
    if (exists(".save_plot_publication", mode = "function")) {
      .save_plot_publication(
        filename_base = filename_base,
        plot = plot,
        width = width,
        height = height,
        dpi = dpi
      )
    } else {
      ggplot2::ggsave(
        filename = paste0(filename_base, ".png"),
        plot = plot,
        width = width,
        height = height,
        dpi = dpi,
        bg = "white"
      )
      ggplot2::ggsave(
        filename = paste0(filename_base, ".pdf"),
        plot = plot,
        width = width,
        height = height,
        device = cairo_pdf,
        bg = "white"
      )
    }
    invisible(NULL)
  }
  
  filter_signal <- function(df) {
    if (nrow(df) == 0) return(df)
    
    if ("signal_label" %in% colnames(df)) {
      return(
        df %>%
          dplyr::filter(
            is.na(signal_label) |
              signal_label == "" |
              as.character(signal_label) == !!signal_label
          )
      )
    }
    
    if ("signal_mode" %in% colnames(df)) {
      expected_mode <- dplyr::recode(
        signal_label,
        observado = "observed",
        reconstruido = "reconstructed",
        .default = signal_label
      )
      
      return(
        df %>%
          dplyr::filter(
            is.na(signal_mode) |
              signal_mode == "" |
              as.character(signal_mode) %in% c(expected_mode, signal_label)
          )
      )
    }
    
    df
  }
  
  coalesce_text_cols <- function(df, cols, default = NA_character_) {
    cols <- intersect(cols, colnames(df))
    if (length(cols) == 0) return(rep(default, nrow(df)))
    
    out <- rep(NA_character_, nrow(df))
    for (cc in cols) {
      val <- as.character(df[[cc]])
      val <- trimws(val)
      val[val == "" | val == "NA"] <- NA_character_
      idx <- is.na(out) & !is.na(val)
      out[idx] <- val[idx]
    }
    out[is.na(out)] <- default
    out
  }
  
  get_num_col <- function(df, nm, default = NA_real_) {
    if (nm %in% colnames(df)) suppressWarnings(as.numeric(df[[nm]])) else rep(default, nrow(df))
  }
  
  wrap_two_lines <- function(x,
                             width = 36,
                             max_lines = 4,
                             add_ellipsis = FALSE) {
    
    x <- as.character(x)
    x <- trimws(x)
    x[is.na(x) | x == "" | x == "NA"] <- "sin anotación"
    
    ## Limpiar prefijos técnicos largos que no ayudan visualmente
    x <- gsub("curated_clinical_concept:\\s*-\\s*", "", x, ignore.case = TRUE)
    x <- gsub("Resolved identifier:\\s*", "", x, ignore.case = TRUE)
    x <- gsub("efo:EFO:[0-9]+\\s*-\\s*", "", x, ignore.case = TRUE)
    x <- gsub("mesh:mesh:[A-Za-z0-9]+\\s*-\\s*", "", x, ignore.case = TRUE)
    x <- gsub("\\[\\.\\.\\.\\]", "", x)
    x <- gsub("\\.\\.\\.", "", x)
    x <- gsub("\\s*\\|\\s*", " | ", x)
    x <- gsub("\\s+", " ", x)
    x <- trimws(x)
    
    wrap_one <- function(xx) {
      
      lines <- strwrap(xx, width = width)
      
      if (length(lines) == 0) {
        lines <- xx
      }
      
      if (!is.null(max_lines) && is.finite(max_lines) && length(lines) > max_lines) {
        
        lines <- lines[seq_len(max_lines)]
        
        if (isTRUE(add_ellipsis)) {
          lines[max_lines] <- paste0(
            substr(lines[max_lines], 1, max(1, width - 3)),
            "..."
          )
        }
      }
      
      paste(lines, collapse = "\n")
    }
    
    vapply(x, wrap_one, character(1))
  }
  
  pattern_short <- function(x) {
    dplyr::recode(
      as.character(x),
      NP_like_high      = "NP-like",
      FD_like_high      = "FD-like",
      SO_like_high      = "SO-like",
      shared_FD_SO_high = "shared FD/SO",
      FD_gt_SO_splitter = "FD > SO split",
      SO_gt_FD_splitter = "SO > FD split",
      uncertain         = "uncertain",
      .default = as.character(x)
    )
  }
  
  evidence_short <- function(x) {
    dplyr::recode(
      as.character(x),
      A_strong                = "A",
      B_moderate              = "B",
      C_exploratory           = "C",
      D_ambiguous_dependent   = "D",
      E_no_evidence           = "E",
      .default = as.character(x)
    )
  }
  
  ## ----------------------------------------------------------
  ## leer tablas
  ## ----------------------------------------------------------
  
  feature_signature <- read_tbl_local(
    file.path(bidir_tables_dir, "01_bidirectional_feature_signature.csv"),
    required = TRUE
  )
  
  feature_annotation <- read_tbl_local(
    file.path(annotation_tables_dir, "35_feature_database_annotation.csv"),
    required = FALSE
  ) %>% filter_signal()
  
  module_members_annotated <- read_tbl_local(
    file.path(annotation_tables_dir, "36_module_members_database_annotation.csv"),
    required = FALSE
  ) %>% filter_signal()
  
  module_summary_37 <- read_tbl_local(
    file.path(annotation_tables_dir, "37_module_database_annotation_summary.csv"),
    required = FALSE
  ) %>% filter_signal()
  
  module_plot_42 <- read_tbl_local(
    file.path(
      annotation_tables_dir,
      paste0("42_final_annotated_network_modules_", signal_label, ".csv")
    ),
    required = FALSE
  )
  
  ## ----------------------------------------------------------
  ## preparar feature_signature con fallbacks
  ## ----------------------------------------------------------
  
  feature_signature$feature_model <- as.character(feature_signature$feature_model)
  
  fallback_chr <- function(df, nm, default = NA_character_) {
    if (!nm %in% colnames(df)) df[[nm]] <- default
    df
  }
  
  fallback_num <- function(df, nm, default = NA_real_) {
    if (!nm %in% colnames(df)) df[[nm]] <- default
    df
  }
  
  if (!"feature_label" %in% colnames(feature_signature)) {
    feature_signature$feature_label <- sub("^[^_]+_", "", feature_signature$feature_model)
  }
  
  feature_signature <- feature_signature %>%
    { fallback_chr(., "posterior_high_group") } %>%
    { fallback_chr(., "top_location_group") } %>%
    { fallback_chr(., "bidirectional_pattern") } %>%
    { fallback_chr(., "evidence_grade") } %>%
    { fallback_chr(., "view_code") } %>%
    { fallback_chr(., "concordant_high", NA) } %>%
    { fallback_num(., "MI_norm") } %>%
    { fallback_num(., "posterior_enrichment_high") } %>%
    { fallback_num(., "q_perm") }
  
  mu_cols <- intersect(paste0("mu_", y_levels), colnames(feature_signature))
  post_cols <- intersect(paste0("post_high_", y_levels), colnames(feature_signature))
  
  if (length(mu_cols) != 3 || length(post_cols) != 3) {
    stop("No encuentro todas las columnas mu_* y post_high_* esperadas.")
  }
  
  ## ----------------------------------------------------------
  ## preparar módulos
  ## ----------------------------------------------------------
  
  if (nrow(module_members_annotated) > 0) {
    feature_key <- if ("feature_model" %in% colnames(module_members_annotated)) {
      "feature_model"
    } else if ("node_id" %in% colnames(module_members_annotated)) {
      "node_id"
    } else {
      NA_character_
    }
    
    if (!is.na(feature_key)) {
      module_members_compact <- module_members_annotated %>%
        dplyr::transmute(
          feature_model = as.character(.data[[feature_key]]),
          module_id = if ("module_id" %in% colnames(module_members_annotated)) as.character(module_id) else NA_character_
        ) %>%
        dplyr::filter(!is.na(feature_model), feature_model != "") %>%
        dplyr::distinct(feature_model, .keep_all = TRUE)
    } else {
      module_members_compact <- tibble::tibble(feature_model = character(), module_id = character())
    }
  } else {
    module_members_compact <- tibble::tibble(feature_model = character(), module_id = character())
  }
  
  if (nrow(module_plot_42) > 0 && "module_id" %in% colnames(module_plot_42)) {
    module_compact <- module_plot_42 %>%
      dplyr::mutate(module_id = as.character(module_id)) %>%
      dplyr::distinct(module_id, .keep_all = TRUE) %>%
      dplyr::transmute(
        module_id,
        module_axis = coalesce_text_cols(
          .,
          c("module_axis_label", "dominant_axis", "dominant_phenotype_module"),
          default = NA_character_
        ),
        module_label_raw = coalesce_text_cols(
          .,
          c(
            "module_axis_label",
            "module_annotation_label",
            "module_top_gene_protein_terms",
            "top_feature_annotations",
            "module_interpretation",
            "dominant_phenotype_module"
          ),
          default = NA_character_
        ),
        module_n_database_supported = get_num_col(., "n_database_supported"),
        module_n_unidentified = get_num_col(., "n_unidentified_metabolites")
      )
  } else if (nrow(module_summary_37) > 0 && "module_id" %in% colnames(module_summary_37)) {
    module_compact <- module_summary_37 %>%
      dplyr::mutate(module_id = as.character(module_id)) %>%
      dplyr::distinct(module_id, .keep_all = TRUE) %>%
      dplyr::transmute(
        module_id,
        module_axis = coalesce_text_cols(
          .,
          c("dominant_phenotype_module"),
          default = NA_character_
        ),
        module_label_raw = coalesce_text_cols(
          .,
          c("module_interpretation", "top_feature_annotations", "dominant_phenotype_module"),
          default = NA_character_
        ),
        module_n_database_supported = get_num_col(., "n_database_supported"),
        module_n_unidentified = get_num_col(., "n_unidentified_metabolites")
      )
  } else {
    module_compact <- tibble::tibble(
      module_id = character(),
      module_axis = character(),
      module_label_raw = character(),
      module_n_database_supported = numeric(),
      module_n_unidentified = numeric()
    )
  }
  
  if (nrow(feature_annotation) > 0 && "feature_model" %in% colnames(feature_annotation)) {
    feature_annotation_compact <- feature_annotation %>%
      dplyr::mutate(feature_model = as.character(feature_model)) %>%
      dplyr::distinct(feature_model, .keep_all = TRUE) %>%
      dplyr::transmute(
        feature_model,
        database_supported_annotation = coalesce_text_cols(
          .,
          c("database_supported_annotation", "phenotype_axis_interpretation",
            "ontology_label", "manual_ontology_label", "metabolite_name",
            "annotation_status"),
          default = NA_character_
        ),
        feature_annotation_status = coalesce_text_cols(
          .,
          c("annotation_status", "evidence_type"),
          default = NA_character_
        )
      )
  } else {
    feature_annotation_compact <- tibble::tibble(
      feature_model = character(),
      database_supported_annotation = character(),
      feature_annotation_status = character()
    )
  }
  
  ## ----------------------------------------------------------
  ## tabla base final
  ## ----------------------------------------------------------
  
  heatmap_base <- feature_signature %>%
    dplyr::mutate(
      feature_model = as.character(feature_model),
      feature_label = as.character(feature_label),
      feature_label = trimws(feature_label),
      feature_label = dplyr::na_if(feature_label, ""),
      feature_label = dplyr::coalesce(feature_label, feature_model),
      
      view_code = as.character(view_code),
      evidence_grade = as.character(evidence_grade),
      bidirectional_pattern = as.character(bidirectional_pattern),
      posterior_high_group = as.character(posterior_high_group),
      top_location_group = as.character(top_location_group),
      
      MI_norm = suppressWarnings(as.numeric(MI_norm)),
      posterior_enrichment_high = suppressWarnings(as.numeric(posterior_enrichment_high)),
      q_perm = suppressWarnings(as.numeric(q_perm)),
      
      concordant_high = as.character(concordant_high)
    ) %>%
    dplyr::left_join(module_members_compact, by = "feature_model") %>%
    dplyr::left_join(module_compact, by = "module_id") %>%
    dplyr::left_join(feature_annotation_compact, by = "feature_model") %>%
    dplyr::mutate(
      module_id = dplyr::if_else(
        is.na(module_id) | module_id == "",
        "UNASSIGNED",
        module_id
      ),
      module_axis = dplyr::coalesce(module_axis, pattern_short(bidirectional_pattern), "unassigned"),
      module_label_raw = dplyr::coalesce(
        module_label_raw,
        database_supported_annotation,
        module_axis,
        "sin anotación"
      ),
      module_label = wrap_two_lines(
        module_label_raw,
        width = 42,
        max_lines = 4,
        add_ellipsis = FALSE
      ),
      
      feature_short = wrap_two_lines(
        feature_label,
        width = 24,
        max_lines = 2,
        add_ellipsis = FALSE
      ),
      pattern_label = pattern_short(bidirectional_pattern),
      evidence_label = evidence_short(evidence_grade),
      top_to_high = paste0(
        ifelse(is.na(top_location_group) | top_location_group == "", "?", top_location_group),
        " \u2192 ",
        ifelse(is.na(posterior_high_group) | posterior_high_group == "", "?", posterior_high_group)
      ),
      concordance_label = dplyr::case_when(
        tolower(concordant_high) %in% c("true", "t", "1") ~ "yes",
        tolower(concordant_high) %in% c("false", "f", "0") ~ "no",
        TRUE ~ "?"
      )
    ) %>%
    dplyr::filter(
      evidence_grade %in% c("A_strong", "B_moderate", "C_exploratory", "D_ambiguous_dependent"),
      bidirectional_pattern != "uncertain"
    )
  
  if (nrow(heatmap_base) == 0) {
    cat("\n[HEATMAP BIDIR CONTEXTUAL] No hay features elegibles.\n")
    return(invisible(NULL))
  }
  
  ## ----------------------------------------------------------
  ## contextos: SIN GLOBAL y SIN top-N
  ## ----------------------------------------------------------
  
  context_defs <- tibble::tibble(
    context_id = c("NP", "FD", "SO", "NP_vs_FD", "NP_vs_SO", "FD_vs_SO"),
    context_title = c("NP", "FD", "SO", "NP vs FD", "NP vs SO", "FD vs SO")
  )
  
  select_context_df <- function(df, context_id) {
    if (context_id == "NP") {
      return(df %>% dplyr::filter(bidirectional_pattern == "NP_like_high"))
    }
    
    if (context_id == "FD") {
      return(df %>% dplyr::filter(bidirectional_pattern == "FD_like_high"))
    }
    
    if (context_id == "SO") {
      return(df %>% dplyr::filter(bidirectional_pattern == "SO_like_high"))
    }
    
    if (context_id == "NP_vs_FD") {
      return(
        df %>%
          dplyr::filter(bidirectional_pattern %in% c("NP_like_high", "FD_like_high"))
      )
    }
    
    if (context_id == "NP_vs_SO") {
      return(
        df %>%
          dplyr::filter(bidirectional_pattern %in% c("NP_like_high", "SO_like_high"))
      )
    }
    
    if (context_id == "FD_vs_SO") {
      return(
        df %>%
          dplyr::filter(
            bidirectional_pattern %in% c(
              "FD_like_high",
              "SO_like_high",
              "shared_FD_SO_high",
              "FD_gt_SO_splitter",
              "SO_gt_FD_splitter"
            )
          )
      )
    }
    
    df
  }
  
  ## ----------------------------------------------------------
  ## graficar un contexto
  ## ----------------------------------------------------------
  
  make_one_context_plot <- function(df_ctx, context_id, context_title) {
    
    if (nrow(df_ctx) == 0) {
      cat("\n[HEATMAP BIDIR CONTEXTUAL] Sin features para ", context_id, "\n")
      return(NULL)
    }
    
    evidence_order <- c(
      "A_strong",
      "B_moderate",
      "C_exploratory",
      "D_ambiguous_dependent"
    )
    
    pattern_order <- c(
      "NP_like_high",
      "FD_like_high",
      "SO_like_high",
      "shared_FD_SO_high",
      "FD_gt_SO_splitter",
      "SO_gt_FD_splitter"
    )
    
    evidence_order <- c(
      "A_strong",
      "B_moderate",
      "C_exploratory",
      "D_ambiguous_dependent"
    )
    
    pattern_order <- c(
      "NP_like_high",
      "FD_like_high",
      "SO_like_high",
      "shared_FD_SO_high",
      "FD_gt_SO_splitter",
      "SO_gt_FD_splitter"
    )
    
    ## ----------------------------------------------------------
    ## Orden interno del contexto
    ## ----------------------------------------------------------
    ## Se ordena por:
    ##   1) evidencia A > B > C > D
    ##   2) patrón bidireccional
    ##   3) q_perm
    ##   4) MI_norm
    ##   5) Enrich
    ## ----------------------------------------------------------
    
    df_ctx <- df_ctx %>%
      dplyr::mutate(
        evidence_grade = as.character(evidence_grade),
        bidirectional_pattern = as.character(bidirectional_pattern),
        
        evidence_rank = match(evidence_grade, evidence_order),
        evidence_rank = ifelse(
          is.na(evidence_rank),
          length(evidence_order) + 1L,
          evidence_rank
        ),
        
        pattern_rank = match(bidirectional_pattern, pattern_order),
        pattern_rank = ifelse(
          is.na(pattern_rank),
          length(pattern_order) + 1L,
          pattern_rank
        ),
        
        q_perm_sort = suppressWarnings(as.numeric(q_perm)),
        q_perm_sort = ifelse(
          is.finite(q_perm_sort),
          q_perm_sort,
          Inf
        ),
        
        MI_norm_sort = suppressWarnings(as.numeric(MI_norm)),
        MI_norm_sort = ifelse(
          is.finite(MI_norm_sort),
          MI_norm_sort,
          -Inf
        ),
        
        enrich_sort = suppressWarnings(as.numeric(posterior_enrichment_high)),
        enrich_sort = ifelse(
          is.finite(enrich_sort),
          enrich_sort,
          -Inf
        )
      ) %>%
      dplyr::arrange(
        evidence_rank,
        pattern_rank,
        q_perm_sort,
        dplyr::desc(MI_norm_sort),
        dplyr::desc(enrich_sort),
        feature_label
      ) %>%
      dplyr::mutate(
        row_id = dplyr::row_number(),
        y_pos = rev(row_id),
        
        ## Sin vista en el eje Y.
        ## Si hay nombres repetidos, make.unique evita solapamientos internos.
        feature_plot = make.unique(as.character(feature_short))
      )
    
    ## ----------------------------------------------------------
    ## Separadores horizontales por evidencia
    ## ----------------------------------------------------------
    ## Este objeto es el que usan p_info, p_mu y p_post.
    ## Antes faltaba y por eso fallaba:
    ##   object 'evidence_bounds' not found
    ## ----------------------------------------------------------
    
    evidence_bounds <- df_ctx %>%
      dplyr::mutate(
        evidence_label = as.character(evidence_label),
        evidence_label = dplyr::case_when(
          is.na(evidence_label) | evidence_label == "" ~ evidence_grade,
          TRUE ~ evidence_label
        )
      ) %>%
      dplyr::group_by(evidence_label) %>%
      dplyr::summarise(
        ymin = min(y_pos, na.rm = TRUE),
        ymax = max(y_pos, na.rm = TRUE),
        ymid = mean(c(ymin, ymax), na.rm = TRUE),
        .groups = "drop"
      )
    
    evidence_bounds <- df_ctx %>%
      dplyr::group_by(evidence_label) %>%
      dplyr::summarise(
        ymin = min(y_pos),
        ymax = max(y_pos),
        ymid = mean(c(min(y_pos), max(y_pos))),
        .groups = "drop"
      )
    
    ## ---------------------------
    ## panel info / anotación
    ## ---------------------------
    
    info_df <- df_ctx %>%
      dplyr::transmute(
        y_pos,
        feature_plot,
        col_module  = module_label,
        col_view    = view_code,
        col_evid    = evidence_label,
        col_tophigh = top_to_high,
        col_pattern = pattern_label,
        col_mi      = sprintf("%.2f", MI_norm),
        col_enrich  = sprintf("%.2f", posterior_enrichment_high)
      )
    
    info_long <- info_df %>%
      tidyr::pivot_longer(
        cols = c(
          col_module,
          col_view,
          col_evid,
          col_tophigh,
          col_pattern,
          col_mi,
          col_enrich
        ),
        names_to = "colname",
        values_to = "label"
      ) %>%
      dplyr::mutate(
        x_pos = dplyr::case_when(
          colname == "col_module"  ~ 1.00,
          colname == "col_view"    ~ 3.25,
          colname == "col_evid"    ~ 3.85,
          colname == "col_tophigh" ~ 4.65,
          colname == "col_pattern" ~ 5.65,
          colname == "col_mi"      ~ 6.55,
          colname == "col_enrich"  ~ 7.35,
          TRUE ~ NA_real_
        ),
        hjust_val = dplyr::case_when(
          colname == "col_module" ~ 0,
          TRUE ~ 0.5
        ),
        text_size = dplyr::case_when(
          colname == "col_module" ~ 2.35,
          TRUE ~ 2.55
        )
      )
    
    p_info <- ggplot2::ggplot(info_long, ggplot2::aes(x = x_pos, y = y_pos)) +
      ggplot2::geom_text(
        ggplot2::aes(
          label = label,
          hjust = hjust_val,
          size = text_size
        ),
        lineheight = 0.96
      ) +
      ggplot2::scale_size_identity() +
      ggplot2::geom_hline(
        data = evidence_bounds,
        ggplot2::aes(yintercept = ymin - 0.5),
        linewidth = 0.25,
        color = "grey60"
      ) +
      ggplot2::scale_x_continuous(
        breaks = c(1.00, 3.25, 3.85, 4.65, 5.65, 6.55, 7.35),
        labels = c(
          "Módulo",
          "Vista",
          "Evid.",
          "Top→High",
          "Patrón",
          "MI",
          "Enrich"
        ),
        position = "top",
        expand = ggplot2::expansion(mult = c(0.01, 0.01))
      ) +
      ggplot2::scale_y_continuous(
        breaks = df_ctx$y_pos,
        labels = df_ctx$feature_plot,
        expand = ggplot2::expansion(mult = c(0.01, 0.01))
      ) +
      ggplot2::coord_cartesian(xlim = c(0.95, 7.65), clip = "off") +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(face = "bold"),
        axis.text.y = ggplot2::element_text(size = 7.8),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
    
    ## ---------------------------
    ## panel mu: grupo -> feature
    ## ---------------------------
    
    mu_long <- df_ctx %>%
      dplyr::select(feature_model, y_pos, dplyr::all_of(mu_cols)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(mu_cols),
        names_to = "metric",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        group = sub("^mu_", "", metric),
        x_pos = match(group, y_levels),
        label = sprintf("%.2f", value)
      )
    
    mu_lim <- max(abs(mu_long$value), na.rm = TRUE)
    if (!is.finite(mu_lim) || mu_lim <= 0) mu_lim <- 1
    
    p_mu <- ggplot2::ggplot(mu_long, ggplot2::aes(x = x_pos, y = y_pos, fill = value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.25) +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 2.7) +
      ggplot2::geom_hline(
        data = evidence_bounds,
        ggplot2::aes(yintercept = ymin - 0.5),
        linewidth = 0.25,
        color = "grey60"
      ) +
      ggplot2::scale_x_continuous(
        breaks = 1:3,
        labels = paste0("\u03bc ", y_levels),
        position = "top",
        expand = ggplot2::expansion(mult = c(0.02, 0.02))
      ) +
      ggplot2::scale_y_continuous(
        breaks = df_ctx$y_pos,
        labels = rep("", nrow(df_ctx)),
        expand = ggplot2::expansion(mult = c(0.01, 0.01))
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#b2182b",
        mid = "white",
        high = "#2166ac",
        midpoint = 0,
        limits = c(-mu_lim, mu_lim),
        name = expression(mu["grupo"])
      ) +
      ggplot2::labs(
        title = "grupo \u2192 feature",
        x = NULL, y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(face = "bold"),
        axis.text.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
    
    ## ---------------------------
    ## panel posterior: feature alta -> grupo
    ## ---------------------------
    
    post_long <- df_ctx %>%
      dplyr::select(feature_model, y_pos, dplyr::all_of(post_cols)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(post_cols),
        names_to = "metric",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        group = sub("^post_high_", "", metric),
        x_pos = match(group, y_levels),
        label = sprintf("%.2f", value)
      )
    
    p_post <- ggplot2::ggplot(post_long, ggplot2::aes(x = x_pos, y = y_pos, fill = value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.25) +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 2.7) +
      ggplot2::geom_hline(
        data = evidence_bounds,
        ggplot2::aes(yintercept = ymin - 0.5),
        linewidth = 0.25,
        color = "grey60"
      ) +
      ggplot2::scale_x_continuous(
        breaks = 1:3,
        labels = paste0("P ", y_levels),
        position = "top",
        expand = ggplot2::expansion(mult = c(0.02, 0.02))
      ) +
      ggplot2::scale_y_continuous(
        breaks = df_ctx$y_pos,
        labels = rep("", nrow(df_ctx)),
        expand = ggplot2::expansion(mult = c(0.01, 0.01))
      ) +
      ggplot2::scale_fill_gradient(
        low = "white",
        high = "#2166ac",
        limits = c(0, 1),
        name = "P(grupo | alta)"
      ) +
      ggplot2::labs(
        title = "feature alta \u2192 grupo",
        x = NULL, y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(face = "bold"),
        axis.text.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
    
    combined_plot <-
      combined_plot <-
      p_info + p_mu + p_post +
      patchwork::plot_layout(widths = c(11.5, 2.6, 2.6)) +
      patchwork::plot_annotation(
        title = paste0(
          "Mapa bidireccional anotado por módulos biológicos: ",
          context_title, " / ", signal_label
        ),
        subtitle = paste0(
          "Panel izquierdo: anotación y resumen por feature. ",
          "Centro: medias robustas por grupo (\u03bc_NP, \u03bc_FD, \u03bc_SO). ",
          "Derecha: posterior P(grupo | feature alta). ",
          "Top→High = top_location_group \u2192 posterior_high_group; ",
          "Enrich = posterior_enrichment_high."
        ),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 13),
          plot.subtitle = ggplot2::element_text(size = 9)
        )
      )
    
    if (isTRUE(save_tables)) {
      write_tbl_local(
        df_ctx,
        file.path(
          tables_dir,
          paste0(
            "05_final_annotated_bidirectional_module_heatmap_input_",
            signal_label, "_", context_id, ".csv"
          )
        )
      )
    }
    
    filename_base <- file.path(
      plots_dir,
      paste0(
        "05_final_annotated_bidirectional_module_heatmap_",
        signal_label, "_", context_id
      )
    )
    
    save_plot_local(
      filename_base = filename_base,
      plot = combined_plot,
      width = 28,
      height = max(8.0, 0.56 * nrow(df_ctx) + 3.4),
      dpi = 600
    )
    invisible(list(
      context_id = context_id,
      input = df_ctx,
      plot = combined_plot,
      filename_base = filename_base
    ))
  }
  
  ## ----------------------------------------------------------
  ## ejecutar contextos
  ## ----------------------------------------------------------
  
  outputs <- list()
  
  for (ii in seq_len(nrow(context_defs))) {
    ctx_id <- context_defs$context_id[ii]
    ctx_title <- context_defs$context_title[ii]
    
    df_ctx <- select_context_df(heatmap_base, ctx_id)
    
    outputs[[ctx_id]] <- make_one_context_plot(
      df_ctx = df_ctx,
      context_id = ctx_id,
      context_title = ctx_title
    )
  }
  
  manifest <- tibble::tibble(
    signal_label = signal_label,
    context_id = names(outputs),
    generated = vapply(outputs, function(x) !is.null(x), logical(1)),
    plot_png = file.path(
      plots_dir,
      paste0("05_final_annotated_bidirectional_module_heatmap_", signal_label, "_", names(outputs), ".png")
    ),
    plot_pdf = file.path(
      plots_dir,
      paste0("05_final_annotated_bidirectional_module_heatmap_", signal_label, "_", names(outputs), ".pdf")
    )
  )
  
  write_tbl_local(
    manifest,
    file.path(
      tables_dir,
      paste0("00_final_annotated_bidirectional_module_heatmap_manifest_", signal_label, ".csv")
    )
  )
  
  cat("\n============================================================\n")
  cat("HEATMAPS BIDIRECCIONALES ANOTADOS GENERADOS:", signal_label, "\n")
  cat("Outdir:\n")
  cat(outdir, "\n")
  print(manifest)
  cat("============================================================\n")
  
  invisible(list(
    heatmap_base = heatmap_base,
    outputs = outputs,
    manifest = manifest,
    outdir = outdir
  ))
}


## ============================================================

## ============================================================
## HEATMAPS BIDIRECCIONALES ANOTADOS POR CONTEXTO
## ------------------------------------------------------------
## No reemplaza ni elimina la red anotada.
## Genera figuras por contexto:
##   NP, FD, SO, NP_vs_FD, NP_vs_SO, FD_vs_SO
##
## Paneles:
##   A) tabla/anotación (módulo, vista, evidencia, top->high, patrón, MI, enrich)
##   B) grupo -> feature : mu_NP / mu_FD / mu_SO
##   C) feature alta -> grupo : post_high_NP / post_high_FD / post_high_SO
##
## Usa:
##   - 01_bidirectional_feature_signature.csv
##   - 35_feature_database_annotation.csv
##   - 36_module_members_database_annotation.csv
##   - 37_module_database_annotation_summary.csv
##   - 42_final_annotated_network_modules_<signal>.csv   (si existe)
## ============================================================

plot_final_annotated_bidirectional_module_heatmap <- function(
    scenario_dir,
    diag_target = "diagnostico_bayes_features",
    signal_label = c("observado", "reconstruido"),
    y_levels = c("NP", "FD", "SO"),
    save_tables = TRUE
) {
  
  signal_label <- match.arg(signal_label)
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Necesitas instalar patchwork: install.packages('patchwork')")
  }
  
  scenario_dir <- normalizePath(path.expand(scenario_dir), mustWork = TRUE)
  
  bidir_tables_dir <- file.path(
    scenario_dir,
    diag_target,
    paste0("framework_bidireccional_", signal_label),
    "tables"
  )
  
  annotation_tables_dir <- file.path(
    scenario_dir,
    diag_target,
    "database_biological_annotation",
    "tables"
  )
  
  outdir <- file.path(
    scenario_dir,
    diag_target,
    paste0("final_annotated_bidirectional_module_heatmap_", signal_label)
  )
  
  plots_dir  <- file.path(outdir, "plots")
  tables_dir <- file.path(outdir, "tables")
  
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## ----------------------------------------------------------
  ## helpers
  ## ----------------------------------------------------------
  
  read_tbl_local <- function(path, required = TRUE) {
    if (!file.exists(path)) {
      if (required) stop("No existe archivo requerido:\n", path)
      return(tibble::tibble())
    }
    tibble::as_tibble(
      read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    )
  }
  
  write_tbl_local <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (exists(".safe_write_csv_annotation", mode = "function")) {
      .safe_write_csv_annotation(x, path)
    } else {
      write.csv(x, path, row.names = FALSE, na = "")
    }
    invisible(path)
  }
  
  save_plot_local <- function(filename_base, plot, width = 18, height = 10, dpi = 600) {
    if (exists(".save_plot_publication", mode = "function")) {
      .save_plot_publication(
        filename_base = filename_base,
        plot = plot,
        width = width,
        height = height,
        dpi = dpi
      )
    } else {
      ggplot2::ggsave(
        filename = paste0(filename_base, ".png"),
        plot = plot,
        width = width,
        height = height,
        dpi = dpi,
        bg = "white"
      )
      ggplot2::ggsave(
        filename = paste0(filename_base, ".pdf"),
        plot = plot,
        width = width,
        height = height,
        device = cairo_pdf,
        bg = "white"
      )
    }
    invisible(NULL)
  }
  
  filter_signal <- function(df) {
    if (nrow(df) == 0) return(df)
    
    if ("signal_label" %in% colnames(df)) {
      return(
        df %>%
          dplyr::filter(
            is.na(signal_label) |
              signal_label == "" |
              as.character(signal_label) == !!signal_label
          )
      )
    }
    
    if ("signal_mode" %in% colnames(df)) {
      expected_mode <- dplyr::recode(
        signal_label,
        observado = "observed",
        reconstruido = "reconstructed",
        .default = signal_label
      )
      
      return(
        df %>%
          dplyr::filter(
            is.na(signal_mode) |
              signal_mode == "" |
              as.character(signal_mode) %in% c(expected_mode, signal_label)
          )
      )
    }
    
    df
  }
  
  coalesce_text_cols <- function(df, cols, default = NA_character_) {
    cols <- intersect(cols, colnames(df))
    if (length(cols) == 0) return(rep(default, nrow(df)))
    
    out <- rep(NA_character_, nrow(df))
    for (cc in cols) {
      val <- as.character(df[[cc]])
      val <- trimws(val)
      val[val == "" | val == "NA"] <- NA_character_
      idx <- is.na(out) & !is.na(val)
      out[idx] <- val[idx]
    }
    out[is.na(out)] <- default
    out
  }
  
  get_num_col <- function(df, nm, default = NA_real_) {
    if (nm %in% colnames(df)) suppressWarnings(as.numeric(df[[nm]])) else rep(default, nrow(df))
  }
  
  wrap_two_lines <- function(x, width = 28, max_lines = 2) {
    x <- as.character(x)
    x <- trimws(x)
    x[is.na(x) | x == "" | x == "NA"] <- "sin anotación"
    
    x <- gsub("curated_clinical_concept:\\s*-\\s*", "", x, ignore.case = TRUE)
    x <- gsub("Resolved identifier:\\s*", "", x, ignore.case = TRUE)
    x <- gsub("efo:EFO:[0-9]+\\s*-\\s*", "", x, ignore.case = TRUE)
    x <- gsub("mesh:mesh:[A-Za-z0-9]+\\s*-\\s*", "", x, ignore.case = TRUE)
    x <- gsub("\\s*\\|\\s*", " | ", x)
    x <- gsub("\\s+", " ", x)
    
    wrap_one <- function(xx) {
      lines <- strwrap(xx, width = width)
      if (length(lines) == 0) lines <- xx
      if (length(lines) > max_lines) {
        kept <- lines[seq_len(max_lines)]
        kept[max_lines] <- paste0(
          substr(kept[max_lines], 1, max(1, width - 3)),
          "..."
        )
        lines <- kept
      }
      paste(lines, collapse = "\n")
    }
    
    vapply(x, wrap_one, character(1))
  }
  
  pattern_short <- function(x) {
    dplyr::recode(
      as.character(x),
      NP_like_high      = "NP-like",
      FD_like_high      = "FD-like",
      SO_like_high      = "SO-like",
      shared_FD_SO_high = "shared FD/SO",
      FD_gt_SO_splitter = "FD > SO split",
      SO_gt_FD_splitter = "SO > FD split",
      uncertain         = "uncertain",
      .default = as.character(x)
    )
  }
  
  evidence_short <- function(x) {
    dplyr::recode(
      as.character(x),
      A_strong                = "A",
      B_moderate              = "B",
      C_exploratory           = "C",
      D_ambiguous_dependent   = "D",
      E_no_evidence           = "E",
      .default = as.character(x)
    )
  }
  
  ## ----------------------------------------------------------
  ## leer tablas
  ## ----------------------------------------------------------
  
  feature_signature <- read_tbl_local(
    file.path(bidir_tables_dir, "01_bidirectional_feature_signature.csv"),
    required = TRUE
  )
  
  feature_annotation <- read_tbl_local(
    file.path(annotation_tables_dir, "35_feature_database_annotation.csv"),
    required = FALSE
  ) %>% filter_signal()
  
  module_members_annotated <- read_tbl_local(
    file.path(annotation_tables_dir, "36_module_members_database_annotation.csv"),
    required = FALSE
  ) %>% filter_signal()
  
  module_summary_37 <- read_tbl_local(
    file.path(annotation_tables_dir, "37_module_database_annotation_summary.csv"),
    required = FALSE
  ) %>% filter_signal()
  
  module_plot_42 <- read_tbl_local(
    file.path(
      annotation_tables_dir,
      paste0("42_final_annotated_network_modules_", signal_label, ".csv")
    ),
    required = FALSE
  )
  
  ## ----------------------------------------------------------
  ## preparar feature_signature con fallbacks
  ## ----------------------------------------------------------
  
  feature_signature$feature_model <- as.character(feature_signature$feature_model)
  
  fallback_chr <- function(df, nm, default = NA_character_) {
    if (!nm %in% colnames(df)) df[[nm]] <- default
    df
  }
  
  fallback_num <- function(df, nm, default = NA_real_) {
    if (!nm %in% colnames(df)) df[[nm]] <- default
    df
  }
  
  if (!"feature_label" %in% colnames(feature_signature)) {
    feature_signature$feature_label <- sub("^[^_]+_", "", feature_signature$feature_model)
  }
  
  feature_signature <- feature_signature %>%
    { fallback_chr(., "posterior_high_group") } %>%
    { fallback_chr(., "top_location_group") } %>%
    { fallback_chr(., "bidirectional_pattern") } %>%
    { fallback_chr(., "evidence_grade") } %>%
    { fallback_chr(., "view_code") } %>%
    { fallback_chr(., "concordant_high", NA) } %>%
    { fallback_num(., "MI_norm") } %>%
    { fallback_num(., "posterior_enrichment_high") } %>%
    { fallback_num(., "q_perm") }
  
  mu_cols <- intersect(paste0("mu_", y_levels), colnames(feature_signature))
  post_cols <- intersect(paste0("post_high_", y_levels), colnames(feature_signature))
  
  if (length(mu_cols) != 3 || length(post_cols) != 3) {
    stop("No encuentro todas las columnas mu_* y post_high_* esperadas.")
  }
  
  ## ----------------------------------------------------------
  ## preparar módulos
  ## ----------------------------------------------------------
  
  if (nrow(module_members_annotated) > 0) {
    feature_key <- if ("feature_model" %in% colnames(module_members_annotated)) {
      "feature_model"
    } else if ("node_id" %in% colnames(module_members_annotated)) {
      "node_id"
    } else {
      NA_character_
    }
    
    if (!is.na(feature_key)) {
      module_members_compact <- module_members_annotated %>%
        dplyr::transmute(
          feature_model = as.character(.data[[feature_key]]),
          module_id = if ("module_id" %in% colnames(module_members_annotated)) as.character(module_id) else NA_character_
        ) %>%
        dplyr::filter(!is.na(feature_model), feature_model != "") %>%
        dplyr::distinct(feature_model, .keep_all = TRUE)
    } else {
      module_members_compact <- tibble::tibble(feature_model = character(), module_id = character())
    }
  } else {
    module_members_compact <- tibble::tibble(feature_model = character(), module_id = character())
  }
  
  if (nrow(module_plot_42) > 0 && "module_id" %in% colnames(module_plot_42)) {
    module_compact <- module_plot_42 %>%
      dplyr::mutate(module_id = as.character(module_id)) %>%
      dplyr::distinct(module_id, .keep_all = TRUE) %>%
      dplyr::transmute(
        module_id,
        module_axis = coalesce_text_cols(
          .,
          c("module_axis_label", "dominant_axis", "dominant_phenotype_module"),
          default = NA_character_
        ),
        module_label_raw = coalesce_text_cols(
          .,
          c("module_plot_label", "module_annotation_label", "module_interpretation",
            "top_feature_annotations", "dominant_phenotype_module"),
          default = NA_character_
        ),
        module_n_database_supported = get_num_col(., "n_database_supported"),
        module_n_unidentified = get_num_col(., "n_unidentified_metabolites")
      )
  } else if (nrow(module_summary_37) > 0 && "module_id" %in% colnames(module_summary_37)) {
    module_compact <- module_summary_37 %>%
      dplyr::mutate(module_id = as.character(module_id)) %>%
      dplyr::distinct(module_id, .keep_all = TRUE) %>%
      dplyr::transmute(
        module_id,
        module_axis = coalesce_text_cols(
          .,
          c("dominant_phenotype_module"),
          default = NA_character_
        ),
        module_label_raw = coalesce_text_cols(
          .,
          c("module_interpretation", "top_feature_annotations", "dominant_phenotype_module"),
          default = NA_character_
        ),
        module_n_database_supported = get_num_col(., "n_database_supported"),
        module_n_unidentified = get_num_col(., "n_unidentified_metabolites")
      )
  } else {
    module_compact <- tibble::tibble(
      module_id = character(),
      module_axis = character(),
      module_label_raw = character(),
      module_n_database_supported = numeric(),
      module_n_unidentified = numeric()
    )
  }
  
  if (nrow(feature_annotation) > 0 && "feature_model" %in% colnames(feature_annotation)) {
    feature_annotation_compact <- feature_annotation %>%
      dplyr::mutate(feature_model = as.character(feature_model)) %>%
      dplyr::distinct(feature_model, .keep_all = TRUE) %>%
      dplyr::transmute(
        feature_model,
        database_supported_annotation = coalesce_text_cols(
          .,
          c("database_supported_annotation", "phenotype_axis_interpretation",
            "ontology_label", "manual_ontology_label", "metabolite_name",
            "annotation_status"),
          default = NA_character_
        ),
        feature_annotation_status = coalesce_text_cols(
          .,
          c("annotation_status", "evidence_type"),
          default = NA_character_
        )
      )
  } else {
    feature_annotation_compact <- tibble::tibble(
      feature_model = character(),
      database_supported_annotation = character(),
      feature_annotation_status = character()
    )
  }
  
  ## ----------------------------------------------------------
  ## tabla base final
  ## ----------------------------------------------------------
  
  heatmap_base <- feature_signature %>%
    dplyr::mutate(
      feature_model = as.character(feature_model),
      feature_label = as.character(feature_label),
      feature_label = trimws(feature_label),
      feature_label = dplyr::na_if(feature_label, ""),
      feature_label = dplyr::coalesce(feature_label, feature_model),
      
      view_code = as.character(view_code),
      evidence_grade = as.character(evidence_grade),
      bidirectional_pattern = as.character(bidirectional_pattern),
      posterior_high_group = as.character(posterior_high_group),
      top_location_group = as.character(top_location_group),
      
      MI_norm = suppressWarnings(as.numeric(MI_norm)),
      posterior_enrichment_high = suppressWarnings(as.numeric(posterior_enrichment_high)),
      q_perm = suppressWarnings(as.numeric(q_perm)),
      
      concordant_high = as.character(concordant_high)
    ) %>%
    dplyr::left_join(module_members_compact, by = "feature_model") %>%
    dplyr::left_join(module_compact, by = "module_id") %>%
    dplyr::left_join(feature_annotation_compact, by = "feature_model") %>%
    dplyr::mutate(
      module_id = dplyr::if_else(
        is.na(module_id) | module_id == "",
        "UNASSIGNED",
        module_id
      ),
      module_axis = dplyr::coalesce(module_axis, pattern_short(bidirectional_pattern), "unassigned"),
      module_label_raw = dplyr::coalesce(
        module_label_raw,
        database_supported_annotation,
        module_axis,
        "sin anotación"
      ),
      module_label = wrap_two_lines(module_label_raw, width = 28, max_lines = 2),
      
      feature_short = wrap_two_lines(feature_label, width = 22, max_lines = 2),
      pattern_label = pattern_short(bidirectional_pattern),
      evidence_label = evidence_short(evidence_grade),
      top_to_high = paste0(
        ifelse(is.na(top_location_group) | top_location_group == "", "?", top_location_group),
        " \u2192 ",
        ifelse(is.na(posterior_high_group) | posterior_high_group == "", "?", posterior_high_group)
      ),
      concordance_label = dplyr::case_when(
        tolower(concordant_high) %in% c("true", "t", "1") ~ "yes",
        tolower(concordant_high) %in% c("false", "f", "0") ~ "no",
        TRUE ~ "?"
      )
    ) %>%
    dplyr::filter(
      evidence_grade %in% c("A_strong", "B_moderate", "C_exploratory", "D_ambiguous_dependent"),
      bidirectional_pattern != "uncertain"
    )
  
  if (nrow(heatmap_base) == 0) {
    cat("\n[HEATMAP BIDIR CONTEXTUAL] No hay features elegibles.\n")
    return(invisible(NULL))
  }
  
  ## ----------------------------------------------------------
  ## contextos: SIN GLOBAL y SIN top-N
  ## ----------------------------------------------------------
  
  context_defs <- tibble::tibble(
    context_id = c("NP", "FD", "SO", "NP_vs_FD", "NP_vs_SO", "FD_vs_SO"),
    context_title = c("NP", "FD", "SO", "NP vs FD", "NP vs SO", "FD vs SO")
  )
  
  select_context_df <- function(df, context_id) {
    if (context_id == "NP") {
      return(df %>% dplyr::filter(bidirectional_pattern == "NP_like_high"))
    }
    
    if (context_id == "FD") {
      return(df %>% dplyr::filter(bidirectional_pattern == "FD_like_high"))
    }
    
    if (context_id == "SO") {
      return(df %>% dplyr::filter(bidirectional_pattern == "SO_like_high"))
    }
    
    if (context_id == "NP_vs_FD") {
      return(
        df %>%
          dplyr::filter(bidirectional_pattern %in% c("NP_like_high", "FD_like_high"))
      )
    }
    
    if (context_id == "NP_vs_SO") {
      return(
        df %>%
          dplyr::filter(bidirectional_pattern %in% c("NP_like_high", "SO_like_high"))
      )
    }
    
    if (context_id == "FD_vs_SO") {
      return(
        df %>%
          dplyr::filter(
            bidirectional_pattern %in% c(
              "FD_like_high",
              "SO_like_high",
              "shared_FD_SO_high",
              "FD_gt_SO_splitter",
              "SO_gt_FD_splitter"
            )
          )
      )
    }
    
    df
  }
  
  ## ----------------------------------------------------------
  ## graficar un contexto
  ## ----------------------------------------------------------
  
  make_one_context_plot <- function(df_ctx, context_id, context_title) {
    
    if (nrow(df_ctx) == 0) {
      cat("\n[HEATMAP BIDIR CONTEXTUAL] Sin features para ", context_id, "\n")
      return(NULL)
    }
    
    df_ctx <- df_ctx %>%
      dplyr::arrange(
        module_axis,
        module_id,
        factor(evidence_grade,
               levels = c("A_strong", "B_moderate", "C_exploratory", "D_ambiguous_dependent")),
        bidirectional_pattern,
        dplyr::desc(MI_norm),
        feature_label
      ) %>%
      dplyr::mutate(
        row_id = dplyr::row_number(),
        y_pos = rev(row_id),
        feature_plot = paste0(feature_short, " [", view_code, "]")
      )
    
    module_bounds <- df_ctx %>%
      dplyr::group_by(module_id, module_label) %>%
      dplyr::summarise(
        ymin = min(y_pos),
        ymax = max(y_pos),
        ymid = mean(c(min(y_pos), max(y_pos))),
        .groups = "drop"
      )
    
    evidence_bounds <- df_ctx %>%
      dplyr::mutate(
        evidence_label = as.character(evidence_label),
        evidence_label = dplyr::case_when(
          is.na(evidence_label) | evidence_label == "" ~ as.character(evidence_grade),
          TRUE ~ evidence_label
        )
      ) %>%
      dplyr::group_by(evidence_label) %>%
      dplyr::summarise(
        ymin = min(y_pos, na.rm = TRUE),
        ymax = max(y_pos, na.rm = TRUE),
        ymid = mean(c(ymin, ymax), na.rm = TRUE),
        .groups = "drop"
      )
    
    ## ---------------------------
    ## panel info / anotación
    ## ---------------------------
    
    info_df <- df_ctx %>%
      dplyr::transmute(
        y_pos,
        feature_plot,
        col_evid    = evidence_label,
        col_tophigh = top_to_high,
        col_pattern = pattern_label,
        col_mi      = sprintf("%.2f", MI_norm),
        col_enrich  = sprintf("%.2f", posterior_enrichment_high)
      )
    
    info_long <- info_df %>%
      tidyr::pivot_longer(
        cols = c(
          col_evid,
          col_tophigh,
          col_pattern,
          col_mi,
          col_enrich
        ),
        names_to = "colname",
        values_to = "label"
      ) %>%
      dplyr::mutate(
        x_pos = dplyr::case_when(
          colname == "col_evid"    ~ 1,
          colname == "col_tophigh" ~ 2,
          colname == "col_pattern" ~ 3,
          colname == "col_mi"      ~ 4,
          colname == "col_enrich"  ~ 5,
          TRUE ~ NA_real_
        )
      )
    
    p_info <- ggplot2::ggplot(info_long, ggplot2::aes(x = x_pos, y = y_pos)) +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        size = 2.75,
        lineheight = 0.94
      ) +
      ggplot2::geom_hline(
        data = evidence_bounds,
        ggplot2::aes(yintercept = ymin - 0.5),
        linewidth = 0.25,
        color = "grey65"
      ) +
      ggplot2::scale_x_continuous(
        breaks = 1:5,
        labels = c("Evid.", "Top→High", "Patrón", "MI", "Enrich"),
        position = "top",
        expand = ggplot2::expansion(mult = c(0.03, 0.03))
      ) +
      ggplot2::scale_y_continuous(
        breaks = df_ctx$y_pos,
        labels = df_ctx$feature_plot,
        expand = ggplot2::expansion(mult = c(0.01, 0.01))
      ) +
      ggplot2::coord_cartesian(xlim = c(0.5, 5.5), clip = "off") +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(face = "bold"),
        axis.text.y = ggplot2::element_text(size = 7.8),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
    
    ## ---------------------------
    ## panel mu: grupo -> feature
    ## ---------------------------
    
    mu_long <- df_ctx %>%
      dplyr::select(feature_model, y_pos, dplyr::all_of(mu_cols)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(mu_cols),
        names_to = "metric",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        group = sub("^mu_", "", metric),
        x_pos = match(group, y_levels),
        label = sprintf("%.2f", value)
      )
    
    mu_lim <- max(abs(mu_long$value), na.rm = TRUE)
    if (!is.finite(mu_lim) || mu_lim <= 0) mu_lim <- 1
    
    p_mu <- ggplot2::ggplot(mu_long, ggplot2::aes(x = x_pos, y = y_pos, fill = value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.25) +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 2.7) +
      ggplot2::geom_hline(
        data = evidence_bounds,
        ggplot2::aes(yintercept = ymin - 0.5),
        linewidth = 0.25,
        color = "grey65"
      ) +
      ggplot2::scale_x_continuous(
        breaks = 1:3,
        labels = paste0("\u03bc ", y_levels),
        position = "top",
        expand = ggplot2::expansion(mult = c(0.02, 0.02))
      ) +
      ggplot2::scale_y_continuous(
        breaks = df_ctx$y_pos,
        labels = rep("", nrow(df_ctx)),
        expand = ggplot2::expansion(mult = c(0.01, 0.01))
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#b2182b",
        mid = "white",
        high = "#2166ac",
        midpoint = 0,
        limits = c(-mu_lim, mu_lim),
        name = expression(mu["grupo"])
      ) +
      ggplot2::labs(
        title = "grupo \u2192 feature",
        x = NULL, y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(face = "bold"),
        axis.text.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
    
    ## ---------------------------
    ## panel posterior: feature alta -> grupo
    ## ---------------------------
    
    post_long <- df_ctx %>%
      dplyr::select(feature_model, y_pos, dplyr::all_of(post_cols)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(post_cols),
        names_to = "metric",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        group = sub("^post_high_", "", metric),
        x_pos = match(group, y_levels),
        label = sprintf("%.2f", value)
      )
    
    p_post <- ggplot2::ggplot(post_long, ggplot2::aes(x = x_pos, y = y_pos, fill = value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.25) +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 2.7) +
      ggplot2::geom_hline(
        data = evidence_bounds,
        ggplot2::aes(yintercept = ymin - 0.5),
        linewidth = 0.25,
        color = "grey65"
      ) +
      ggplot2::scale_x_continuous(
        breaks = 1:3,
        labels = paste0("P ", y_levels),
        position = "top",
        expand = ggplot2::expansion(mult = c(0.02, 0.02))
      ) +
      ggplot2::scale_y_continuous(
        breaks = df_ctx$y_pos,
        labels = rep("", nrow(df_ctx)),
        expand = ggplot2::expansion(mult = c(0.01, 0.01))
      ) +
      ggplot2::scale_fill_gradient(
        low = "white",
        high = "#2166ac",
        limits = c(0, 1),
        name = "P(grupo | alta)"
      ) +
      ggplot2::labs(
        title = "feature alta \u2192 grupo",
        x = NULL, y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(face = "bold"),
        axis.text.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
    
    combined_plot <-
      p_info + p_mu + p_post +
      patchwork::plot_layout(widths = c(4.8, 2.8, 2.8)) +
      patchwork::plot_annotation(
        title = paste0(
          "Mapa bidireccional por contexto: ",
          context_title, " / ", signal_label
        ),
        subtitle = paste0(
          "Izquierda: evidencia, concordancia Top→High, patrón, MI y Enrich. ",
          "Centro: dirección grupo \u2192 feature mediante \u03bc_NP, \u03bc_FD, \u03bc_SO. ",
          "Derecha: dirección feature alta \u2192 grupo mediante P(grupo | feature alta)."
        ),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 13),
          plot.subtitle = ggplot2::element_text(size = 9)
        )
      )
    
    if (isTRUE(save_tables)) {
      write_tbl_local(
        df_ctx,
        file.path(
          tables_dir,
          paste0(
            "05_final_annotated_bidirectional_module_heatmap_input_",
            signal_label, "_", context_id, ".csv"
          )
        )
      )
    }
    
    filename_base <- file.path(
      plots_dir,
      paste0(
        "05_final_annotated_bidirectional_module_heatmap_",
        signal_label, "_", context_id
      )
    )
    
    save_plot_local(
      filename_base = filename_base,
      plot = combined_plot,
      width = 18,
      height = max(6.5, 0.34 * nrow(df_ctx) + 2.8),
      dpi = 600
    )
    invisible(list(
      context_id = context_id,
      input = df_ctx,
      plot = combined_plot,
      filename_base = filename_base
    ))
  }
  
  ## ----------------------------------------------------------
  ## ejecutar contextos
  ## ----------------------------------------------------------
  
  outputs <- list()
  
  for (ii in seq_len(nrow(context_defs))) {
    ctx_id <- context_defs$context_id[ii]
    ctx_title <- context_defs$context_title[ii]
    
    df_ctx <- select_context_df(heatmap_base, ctx_id)
    
    outputs[[ctx_id]] <- make_one_context_plot(
      df_ctx = df_ctx,
      context_id = ctx_id,
      context_title = ctx_title
    )
  }
  
  manifest <- tibble::tibble(
    signal_label = signal_label,
    context_id = names(outputs),
    generated = vapply(outputs, function(x) !is.null(x), logical(1)),
    plot_png = file.path(
      plots_dir,
      paste0("05_final_annotated_bidirectional_module_heatmap_", signal_label, "_", names(outputs), ".png")
    ),
    plot_pdf = file.path(
      plots_dir,
      paste0("05_final_annotated_bidirectional_module_heatmap_", signal_label, "_", names(outputs), ".pdf")
    )
  )
  
  write_tbl_local(
    manifest,
    file.path(
      tables_dir,
      paste0("00_final_annotated_bidirectional_module_heatmap_manifest_", signal_label, ".csv")
    )
  )
  
  cat("\n============================================================\n")
  cat("HEATMAPS BIDIRECCIONALES ANOTADOS GENERADOS:", signal_label, "\n")
  cat("Outdir:\n")
  cat(outdir, "\n")
  print(manifest)
  cat("============================================================\n")
  
  invisible(list(
    heatmap_base = heatmap_base,
    outputs = outputs,
    manifest = manifest,
    outdir = outdir
  ))
}
## ============================================================
## 4) RUNNER GENERAL
## ============================================================
.load_existing_bidir_result <- function(outdir) {
  
  tables_dir <- file.path(outdir, "tables")
  
  feature_file <- file.path(
    tables_dir,
    "01_bidirectional_feature_signature.csv"
  )
  
  view_file <- file.path(
    tables_dir,
    "02_bidirectional_view_summary.csv"
  )
  
  perm_file <- file.path(
    tables_dir,
    "03_permutation_distribution.csv"
  )
  
  final_file <- file.path(
    tables_dir,
    "04_bidirectional_final_table_corrected.csv"
  )
  
  required_files <- c(
    feature_file,
    view_file,
    perm_file,
    final_file
  )
  
  if (!all(file.exists(required_files))) {
    return(NULL)
  }
  
  feature_signature <- tryCatch(
    read.csv(feature_file, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  
  view_summary <- tryCatch(
    read.csv(view_file, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  
  permutation_distribution <- tryCatch(
    read.csv(perm_file, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  
  final_bidirectional_table <- tryCatch(
    read.csv(final_file, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  
  if (
    is.null(feature_signature) ||
    is.null(view_summary) ||
    is.null(permutation_distribution) ||
    is.null(final_bidirectional_table)
  ) {
    return(NULL)
  }
  
  required_cols <- c(
    "feature_model",
    "view_code",
    "MI_bits",
    "MI_norm",
    "p_perm",
    "q_perm",
    "evidence_grade",
    "bidirectional_pattern"
  )
  
  missing_cols <- setdiff(required_cols, colnames(feature_signature))
  
  if (length(missing_cols) > 0) {
    cat("\n[CACHE BIDIR] Existe tabla previa, pero faltan columnas:\n")
    cat(paste(missing_cols, collapse = ", "), "\n")
    cat("Se recalculará el bidireccional.\n")
    return(NULL)
  }
  
  list(
    feature_signature = tibble::as_tibble(feature_signature),
    final_bidirectional_table = tibble::as_tibble(final_bidirectional_table),
    view_summary = tibble::as_tibble(view_summary),
    permutation_distribution = tibble::as_tibble(permutation_distribution),
    loaded_from_cache = TRUE,
    cache_files = required_files
  )
}
run_bidirectional_pipeline <- function(inp,
                                       signal_label,
                                       scenario_dir,
                                       diag_target,
                                       B_perm = 999,
                                       B_boot = 499,
                                       network_rho_thr = 0.60,
                                       network_max_edges = 80,
                                       seed = 123,
                                       reuse_if_exists = TRUE) {
  
  outdir <- file.path(
    scenario_dir,
    diag_target,
    paste0("framework_bidireccional_", signal_label)
  )
  
  cat("\n============================================================\n")
  cat("EJECUTANDO FRAMEWORK BIDIRECCIONAL:", signal_label, "\n")
  cat("Outdir:\n")
  cat(outdir, "\n")
  cat("============================================================\n")
  
  X_S <- inp$X_S
  y   <- inp$y
  v   <- inp$v
  
  force_recalc <- Sys.getenv("BIDIR_FORCE_RECALC", unset = "0")
  force_recalc <- force_recalc %in% c("1", "TRUE", "true", "yes", "YES", "si", "SI")
  
  res_bidir <- NULL
  bidir_reused <- FALSE
  
  if (reuse_if_exists && !force_recalc) {
    
    res_bidir <- .load_existing_bidir_result(outdir)
    
    if (!is.null(res_bidir)) {
      
      bidir_reused <- TRUE
      
      cat("\n[CACHE BIDIR] Resultado bidireccional ya existe. No se recalcula:\n")
      cat(" - tables/01_bidirectional_feature_signature.csv\n")
      cat(" - tables/02_bidirectional_view_summary.csv\n")
      cat(" - tables/03_permutation_distribution.csv\n")
      cat(" - tables/04_bidirectional_final_table_corrected.csv\n")
    }
  }
  
  if (is.null(res_bidir)) {
    
    if (force_recalc) {
      cat("\n[CACHE BIDIR] BIDIR_FORCE_RECALC=1. Se fuerza recálculo.\n")
    } else {
      cat("\n[CACHE BIDIR] No hay resultado completo previo. Se calcula bidireccional.\n")
    }
    
    res_bidir <- run_bidirectional_explainability(
      X_S = X_S,
      y = y,
      v = v,
      y_levels = c("NP", "FD", "SO"),
      B_perm = B_perm,
      B_boot = B_boot,
      outdir = outdir,
      seed = seed
    )
  }
  
  ## ------------------------------------------------------------
  ## Plots interpretables:
  ## Si el bidireccional se reutiliza y los plots ya existen, no se regeneran.
  ## Si faltan plots, se reconstruyen usando las tablas cargadas.
  ## ------------------------------------------------------------
  
  plot_manifest_file <- file.path(
    outdir,
    "plots_interpretables_bidireccional",
    "00_plot_manifest.csv"
  )
  
  plots_bidir <- NULL
  
  if (bidir_reused && file.exists(plot_manifest_file) && !force_recalc) {
    
    cat("\n[CACHE PLOTS] Plots interpretables ya existen. No se regeneran:\n")
    cat(plot_manifest_file, "\n")
    
    plots_bidir <- list(
      plots_dir = dirname(plot_manifest_file),
      reused_from_cache = TRUE
    )
    
  } else {
    
    plots_bidir <- tryCatch(
      {
        plot_bidirectional_explainability(
          bidir = res_bidir,
          X_S = X_S,
          outdir = outdir,
          top_n = Inf,
          y_levels = c("NP", "FD", "SO"),
          remove_prefix = TRUE,
          network_rho_thr = network_rho_thr,
          network_max_edges = network_max_edges,
          show_opposite_edges = FALSE
        )
      },
      error = function(e) {
        cat("\nERROR EN PLOTS INTERPRETABLES PARA:", signal_label, "\n")
        cat(conditionMessage(e), "\n")
        cat("El pipeline continúa con tablas y siguiente modo de señal.\n")
        NULL
      }
    )
  }
  
  ## ------------------------------------------------------------
  ## Módulos:
  ## Si ya existen los summaries, no repetir.
  ## Si solo existen nodos/aristas, resumir.
  ## ------------------------------------------------------------
  
  nodes_file <- file.path(
    outdir,
    "tables",
    "07_feature_feature_hybrid_network_nodes.csv"
  )
  
  edges_file <- file.path(
    outdir,
    "tables",
    "08_feature_feature_hybrid_network_edges.csv"
  )
  
  module_members_file <- file.path(
    outdir,
    "tables",
    "10_feature_feature_neutral_module_members.csv"
  )
  
  module_summary_file <- file.path(
    outdir,
    "tables",
    "11_feature_feature_neutral_module_summary.csv"
  )
  
  module_template_file <- file.path(
    outdir,
    "tables",
    "12_feature_feature_module_annotation_template.csv"
  )
  
  module_bidir <- NULL
  
  if (
    bidir_reused &&
    file.exists(module_members_file) &&
    file.exists(module_summary_file) &&
    file.exists(module_template_file) &&
    !force_recalc
  ) {
    
    cat("\n[CACHE MODULES] Resumen de módulos ya existe. No se recalcula.\n")
    
    module_bidir <- list(
      module_members = tibble::as_tibble(
        read.csv(module_members_file, check.names = FALSE, stringsAsFactors = FALSE)
      ),
      module_summary = tibble::as_tibble(
        read.csv(module_summary_file, check.names = FALSE, stringsAsFactors = FALSE)
      ),
      module_annotation_template = tibble::as_tibble(
        read.csv(module_template_file, check.names = FALSE, stringsAsFactors = FALSE)
      ),
      reused_from_cache = TRUE
    )
    
  } else if (file.exists(nodes_file) && file.exists(edges_file)) {
    
    module_bidir <- summarize_feature_modules_from_network(
      outdir = outdir,
      min_module_size = 2
    )
    
  } else {
    
    cat("\nNo se creó red feature-feature suficiente para resumir módulos en:\n")
    cat(outdir, "\n")
  }
  
  invisible(
    list(
      signal_label = signal_label,
      outdir = outdir,
      input = inp,
      res_bidir = res_bidir,
      plots_bidir = plots_bidir,
      module_bidir = module_bidir,
      reused_bidir_from_cache = bidir_reused
    )
  )
}


## ============================================================
## 3) DEFINIR ESCENARIO Y MODO DE SEÑAL
## ============================================================

args <- commandArgs(trailingOnly = TRUE)

SCENARIO_DIR <- if (length(args) >= 1) {
  args[1]
} else {
  stop(
    "Debes pasar el directorio del escenario.\n",
    "Ejemplo:\n",
    "Rscript framework_bidireccional.R ./escenarios/MIESCENARIO"
  )
}

DIAG_TARGET <- if (length(args) >= 2) {
  args[2]
} else {
  Sys.getenv("DIAG_TARGET", unset = "diagnostico_bayes_features")
}

SIGNAL_MODE <- if (length(args) >= 3) {
  args[3]
} else {
  Sys.getenv("BIDIR_SIGNAL_MODE", unset = "both")
}

SCENARIO_DIR <- normalizePath(path.expand(SCENARIO_DIR), mustWork = TRUE)

SIGNAL_MODE <- tolower(SIGNAL_MODE)

SIGNAL_MODE <- dplyr::recode(
  SIGNAL_MODE,
  reconstruido = "reconstructed",
  reconstruida = "reconstructed",
  reconstructed = "reconstructed",
  observado = "observed",
  observada = "observed",
  observed = "observed",
  ambos = "both",
  both = "both",
  .default = SIGNAL_MODE
)

if (!SIGNAL_MODE %in% c("reconstructed", "observed", "both")) {
  stop(
    "SIGNAL_MODE inválido: ", SIGNAL_MODE, "\n",
    "Usa: reconstructed, observed o both"
  )
}

signal_modes_to_run <- if (SIGNAL_MODE == "both") {
  c("reconstructed", "observed")
} else {
  SIGNAL_MODE
}

cat("\n============================================================\n")
cat("FRAMEWORK BIDIRECCIONAL\n")
cat("Escenario   :", SCENARIO_DIR, "\n")
cat("Diagnóstico :", DIAG_TARGET, "\n")
cat("Modo señal  :", paste(signal_modes_to_run, collapse = " + "), "\n")
cat("============================================================\n")


## ============================================================
## 5) CREAR INPUTS Y EJECUTAR
## ============================================================

all_bidir_results <- list()

if ("reconstructed" %in% signal_modes_to_run) {
  
  inp_reconstructed <- make_mofa_reconstructed_bidir_input(
    scenario_dir = SCENARIO_DIR,
    diag_target = DIAG_TARGET,
    keep_factors = FALSE,
    save_reconstructed = TRUE
  )
  
  all_bidir_results$reconstructed <- run_bidirectional_pipeline(
    inp = inp_reconstructed,
    signal_label = "reconstruido",
    scenario_dir = SCENARIO_DIR,
    diag_target = DIAG_TARGET,
    B_perm = 999,
    B_boot = 499,
    network_rho_thr = 0.60,
    network_max_edges = 80,
    seed = 123
  )
}

if ("observed" %in% signal_modes_to_run) {
  
  inp_observed <- make_observed_bidir_input(
    scenario_dir = SCENARIO_DIR,
    diag_target = DIAG_TARGET,
    save_observed = TRUE
  )
  
  all_bidir_results$observed <- run_bidirectional_pipeline(
    inp = inp_observed,
    signal_label = "observado",
    scenario_dir = SCENARIO_DIR,
    diag_target = DIAG_TARGET,
    B_perm = 999,
    B_boot = 499,
    network_rho_thr = 0.60,
    network_max_edges = 80,
    seed = 456
  )
}





## ============================================================
## 6) RANKING DIFERENCIAL DUAL SOBRE TODAS LAS FEATURES MOFA
##    observed/sin integrar vs reconstructed/integrado
## ============================================================

ranking_dual_mofa <- tryCatch(
  {
    run_mofa_dual_feature_differential_ranking(
      scenario_dir = SCENARIO_DIR,
      diag_target = DIAG_TARGET,
      all_bidir_results = all_bidir_results
    )
  },
  error = function(e) {
    cat("\nERROR EN RANKING DIFERENCIAL DUAL MOFA\n")
    cat(conditionMessage(e), "\n")
    cat("El framework bidireccional queda completo, pero sin ranking dual.\n")
    NULL
  }
)


## ============================================================
## 7) ANOTACIÓN BIOLÓGICA DINÁMICA CON BASES DE DATOS
##    Sin hardcodear features concretas
## ============================================================

annotation_bio <- tryCatch(
  {
    run_database_biological_annotation_dynamic(
      scenario_dir = SCENARIO_DIR,
      diag_target = DIAG_TARGET,
      organism = "hsapiens",
      gprofiler_threshold = 0.20
    )
  },
  error = function(e) {
    
    cat("\n============================================================\n")
    cat("ERROR EN ANOTACIÓN BIOLÓGICA DINÁMICA\n")
    cat("============================================================\n")
    cat("Mensaje:\n")
    cat(conditionMessage(e), "\n\n")
    
    cat("Llamada asociada al error:\n")
    print(conditionCall(e))
    
    cat("\nTraceback base R:\n")
    traceback()
    
    if (requireNamespace("rlang", quietly = TRUE)) {
      cat("\nTraceback rlang:\n")
      print(rlang::trace_back())
    }
    
    debug_dir <- file.path(
      SCENARIO_DIR,
      DIAG_TARGET,
      "database_biological_annotation",
      "debug"
    )
    
    dir.create(debug_dir, recursive = TRUE, showWarnings = FALSE)
    
    saveRDS(
      e,
      file.path(debug_dir, "annotation_error_condition.rds")
    )
    
    cat("\nError guardado en:\n")
    cat(file.path(debug_dir, "annotation_error_condition.rds"), "\n")
    
    cat("\nEl framework queda completo, pero sin anotación biológica automática.\n")
    NULL
  }
)

## ============================================================
## 8) RED FINAL BIOLÓGICA ANOTADA POR MÓDULOS
## ============================================================

final_annotated_network <- tryCatch(
  {
    if (is.null(annotation_bio)) {
      
      cat("\nNo se genera red biológica anotada porque annotation_bio es NULL.\n")
      NULL
      
    } else {
      
      signal_labels_for_network <- dplyr::recode(
        signal_modes_to_run,
        reconstructed = "reconstruido",
        observed = "observado"
      )
      
      names(signal_labels_for_network) <- signal_labels_for_network
      
      purrr::map(
        signal_labels_for_network,
        function(sig_lab) {
          
          plot_final_annotated_biological_module_network(
            scenario_dir = SCENARIO_DIR,
            diag_target = DIAG_TARGET,
            signal_label = sig_lab,
            network_max_edges = 80,
            min_module_size = 2,
            top_terms_per_module = 4,
            save_tables = TRUE
          )
        }
      )
    }
  },
  error = function(e) {
    
    cat("\n============================================================\n")
    cat("ERROR EN RED FINAL BIOLÓGICA ANOTADA\n")
    cat("============================================================\n")
    cat("Mensaje:\n")
    cat(conditionMessage(e), "\n")
    
    cat("\nEl framework queda completo, pero sin la red biológica anotada final.\n")
    
    NULL
 
     }
)

## ============================================================
## 8B) HEATMAPS FINALES BIDIRECCIONALES ANOTADOS POR CONTEXTO
##     Figura adicional, no reemplaza la red.
## ============================================================

## ============================================================
## 8B) HEATMAPS FINALES BIDIRECCIONALES ANOTADOS POR CONTEXTO
##     Figura adicional, no reemplaza la red.
## ============================================================

final_annotated_bidirectional_heatmap <- tryCatch(
  {
    if (is.null(annotation_bio)) {
      
      cat("\nNo se generan heatmaps bidireccionales anotados porque annotation_bio es NULL.\n")
      NULL
      
    } else {
      
      signal_labels_for_heatmap <- dplyr::recode(
        signal_modes_to_run,
        reconstructed = "reconstruido",
        observed = "observado"
      )
      
      purrr::map(
        signal_labels_for_heatmap,
        function(sig_lab) {
          plot_final_annotated_bidirectional_module_heatmap(
            scenario_dir = SCENARIO_DIR,
            diag_target = DIAG_TARGET,
            signal_label = sig_lab,
            y_levels = c("NP", "FD", "SO"),
            save_tables = TRUE
          )
        }
      )
    }
  },
  error = function(e) {
    
    cat("\n============================================================\n")
    cat("ERROR EN HEATMAPS BIDIRECCIONALES ANOTADOS POR CONTEXTO\n")
    cat("============================================================\n")
    cat("Mensaje:\n")
    cat(conditionMessage(e), "\n")
    
    cat("\nEl framework queda completo, pero sin los heatmaps bidireccionales anotados.\n")
    
    NULL
  }
)

cat("\n============================================================\n")
cat("FRAMEWORK BIDIRECCIONAL COMPLETADO\n")
cat("Resultados guardados en:\n")

for (nm in names(all_bidir_results)) {
  cat(" - ", nm, ": ", all_bidir_results[[nm]]$outdir, "\n", sep = "")
}

if (!is.null(ranking_dual_mofa)) {
  cat(" - ranking_dual_mofa: ", ranking_dual_mofa$outdir, "\n", sep = "")
}

if (!is.null(annotation_bio)) {
  cat(" - database_annotation: ", annotation_bio$outdir, "\n", sep = "")
}

cat("============================================================\n")
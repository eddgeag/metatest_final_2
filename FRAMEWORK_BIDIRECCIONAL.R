## ============================================================
## FRAMEWORK DE EVIDENCIA BIDIRECCIONAL MULTI-ÓMICA
## ============================================================
##
## Objetivo:
##   Evaluar e interpretar biomarcadores multi-ómicos mediante dos
##   direcciones complementarias:
##
##   1) Grupo -> feature:
##        X_j | Y = k
##        Distribución robusta de cada feature dentro de NP, FD y SO.
##
##   2) Feature -> grupo:
##        P(Y = k | X_j alto)
##        Posterior del grupo clínico cuando una feature se encuentra alta.
##
## Modelo estadístico por feature:
##   Y ~ Categorical(pi)
##   X_j | Y = k ~ Student-t(df_t, mu_jk, sigma_jk)
##
## Evidencia estimada:
##   - Información mutua feature-grupo.
##   - Distancia Jensen-Shannon entre grupos.
##   - Permutación de etiquetas para p/q-value.
##   - Bootstrap para estabilidad direccional.
##   - Clasificación interpretable del patrón bidireccional.
##
## Entradas posibles:
##   - Señal observada: datos originales usados en el modelo.
##   - Señal reconstruida: señal integrada reconstruida desde MOFA2.
##
## Salidas principales:
##   - Tabla bidireccional final corregida.
##   - Ranking dual observado/reconstruido.
##   - Anotación biológica dinámica.
##   - Heatmaps finales anotados por contexto.
##
## Nota metodológica:
##   Este script estima asociación, direccionalidad estadística y
##   enriquecimiento posterior. No implementa inferencia causal.
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
                                              outdir = "framework_bidireccional_firma",
                                              top_n = Inf,
                                              y_levels = c("NP", "FD", "SO"),
                                              remove_prefix = TRUE) {
  if (is.null(bidir$feature_signature)) {
    stop("bidir debe tener bidir$feature_signature.")
  }
  
  feature_signature <- tibble::as_tibble(bidir$feature_signature)
  
  PLOTS_DIR <- file.path(outdir, "plots_interpretables_bidireccional")

  
  dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

  
  out_plot <- function(x) file.path(PLOTS_DIR, x)
  
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
  
  ## ============================================================
  ## PLOT 05B. HEATMAP BIDIRECCIONAL DOBLE
  ## feature | grupo  +  grupo | feature alta
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
  ## PLOT 07. CONCORDANCIA BIDIRECCIONAL
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
  ## PLOT 09. RESUMEN MULTI-ÓMICO POR VISTA
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
  
  if (nrow(x) == 0 && ncol(x) == 0) {
    x <- tibble::tibble(.empty_table = character())
    cat("  Tabla 0 x 0 convertida a tabla vacía con cabecera.\n")
  }
  
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
                                               id_resolution) {
  
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
  
  ## Los módulos derivados de red fueron desactivados.
  ## La anotación biológica se realiza a nivel de feature y se resume
  ## posteriormente en el heatmap final anotado.
  
  tibble::tibble()
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
  
  ## Los módulos derivados de red están desactivados.
  ## Se devuelven tablas vacías para conservar compatibilidad
  ## con el resto del pipeline.
  
  list(
    module_members_annotated = tibble::tibble(),
    module_summary = tibble::tibble()
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
    id_resolution = id_resolution
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
## FIGURA FINAL: HEATMAP BIDIRECCIONAL ANOTADO POR CONTEXTO
## ------------------------------------------------------------
## Figura principal final del framework.
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
##   - 36_module_members_database_annotation.csv          (opcional)
##   - 37_module_database_annotation_summary.csv          (opcional)
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
    
    file_size <- file.info(path)$size
    
    if (!is.finite(file_size) || file_size == 0) {
      if (required) {
        stop("Archivo requerido vacío:\n", path)
      }
      return(tibble::tibble())
    }
    
    out <- tryCatch(
      {
        read.csv(
          path,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      },
      error = function(e) {
        if (required) {
          stop(
            "No se pudo leer archivo requerido:\n",
            path,
            "\nMensaje:\n",
            conditionMessage(e)
          )
        }
        return(NULL)
      }
    )
    
    if (is.null(out)) {
      return(tibble::tibble())
    }
    
    tibble::as_tibble(out)
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
  
  if (nrow(module_summary_37) > 0 && "module_id" %in% colnames(module_summary_37)) {
    
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
          outdir = outdir,
          top_n = Inf,
          y_levels = c("NP", "FD", "SO"),
          remove_prefix = TRUE
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
  
  module_bidir <- NULL
  
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

SCENARIOS_ROOT <- if (length(args) >= 1) args[1] else "./escenarios"
DIAG_TARGET    <- if (length(args) >= 2) args[2] else "diagnostico_bayes_features"
SIGNAL_MODE    <- if (length(args) >= 3) args[3] else "both"

SCENARIOS_ROOT <- normalizePath(path.expand(SCENARIOS_ROOT), mustWork = TRUE)

scenario_dirs <- list.dirs(SCENARIOS_ROOT, recursive = FALSE, full.names = TRUE)

scenario_dirs <- scenario_dirs[
  file.exists(file.path(scenario_dirs, DIAG_TARGET, "rds", "diagnostic_objects.rds"))
]

if (length(scenario_dirs) == 0) {
  stop("No encontré escenarios con: ", file.path(DIAG_TARGET, "rds", "diagnostic_objects.rds"))
}

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
  stop("SIGNAL_MODE inválido: ", SIGNAL_MODE)
}

signal_modes_to_run <- if (SIGNAL_MODE == "both") {
  c("reconstructed", "observed")
} else {
  SIGNAL_MODE
}

for (SCENARIO_DIR in scenario_dirs) {
  
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
## 8) FIGURA FINAL: HEATMAPS BIDIRECCIONALES ANOTADOS
## ======================================================
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

}
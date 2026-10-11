#!/usr/bin/env Rscript
# =============================================================================
# METATEST — PASO 5
# Puente cross-ómico clase <-> factor <-> pathway <-> feature
# Version 0.1.0
#
# OBJETIVO
#   Construir, SIN reentrenar MOFA ni utilizar TEST:
#
#       clase
#         <->
#       factor/cara
#         <->
#       pathway transcriptómico
#         <->
#       features TX / PR / ME / CL
#
#   y recuperar en paralelo el alineamiento latente clase-feature:
#
#       A_jc = sum_f L_jf * beta_fc
#
#   junto con su descomposición:
#
#       A_jfc = L_jf * beta_fc
#
#   donde:
#
#       L_jf = orientation_sign_f *
#              train_scale_after_orientation_f *
#              W_jf
#
#   Esta transformación pone W en la misma métrica que los beta de la auditoría
#   latente, porque beta fue estimado sobre Z orientado y estandarizado.
#
# DISEÑO
#   1) Módulos factor-cara:
#      una feature pertenece direccionalmente a la cara del factor si su L_jf
#      tiene el signo de esa cara. No se exige representación mínima por vista.
#
#   2) Importancia dentro de cada vista:
#      se reportan loading, rango y percentil dentro de view x factor.
#      NO se comparan ingenuamente magnitudes crudas entre ómicas.
#
#   3) Pathway -> feature cross-ómica:
#      para cada pathway significativo del bridge clase-pathway se construye un
#      perfil latente transcriptómico P_q a partir de TODOS los genes miembros
#      que apuntan hacia la cara enriquecida, ponderados por su carga en el
#      factor ancla.
#
#          P_q = sum_g a_g * L_g / sum_g a_g
#
#      con:
#
#          a_g = max(side_sign * L_g,f_anchor, 0)
#
#      Después:
#
#          cosine(j,q) = cos(L_j, P_q)
#
#      Esto NO afirma pertenencia bioquímica de PR/ME/CL al pathway. Afirma
#      coherencia en el mismo programa latente multi-ómico.
#
#   4) Clase -> feature:
#      la incertidumbre de beta se propaga desde los fits BRMS ya existentes.
#      Los factores fueron ajustados por separado en la auditoría; por ello la
#      propagación conjunta usa el producto de los posteriores marginales
#      factor-wise (aproximación de independencia).
#
#   5) Biomarcadores:
#      este paso GENERA/PRIORIZA candidatos y aporta anotación mecanística.
#      NO valida biomarcadores. El framework bidireccional TRAIN se añade como
#      evidencia ortogonal, nunca como filtro obligatorio.
#
#   6) TEST:
#      no se usa.
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

# El pathway bridge ya decidió qué evidencia es primaria.
PRIMARY_PATHWAY_LABEL <- "PRIMARY_TWO_SIDE_FDR"

# Solo para exportación compacta de enlaces de mayor coherencia.
# No selecciona biomarcadores ni modifica la tabla ALL.
TOP_LINKS_PER_PATHWAY_VIEW <- as.integer(
  Sys.getenv("CROSSOMIC_TOP_LINKS_PER_PATHWAY_VIEW", unset = "10")
)

NUM_TOL <- 1e-12
# ============================================================================

assert <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
}

canonical_path <- function(p, must_work = TRUE) {
  normalizePath(path.expand(p), winslash = "/", mustWork = must_work)
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (is.null(x) || !ncol(x)) {
    x <- data.frame(status = character())
  }

  utils::write.table(
    x,
    path,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    na = "NA",
    qmethod = "double",
    fileEncoding = "UTF-8"
  )
}

safe_quantiles <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(q025 = NA_real_, q50 = NA_real_, q975 = NA_real_))
  }

  q <- stats::quantile(
    x,
    probs = c(0.025, 0.50, 0.975),
    names = FALSE,
    na.rm = TRUE
  )

  c(q025 = q[1L], q50 = q[2L], q975 = q[3L])
}

empirical_percentile <- function(x) {
  n <- length(x)

  if (n == 0L) return(numeric())
  if (n == 1L) return(1)

  (rank(x, ties.method = "average") - 1) / (n - 1)
}

side_sign_value <- function(side) {
  ifelse(
    side == "positive",
    1,
    ifelse(side == "negative", -1, NA_real_)
  )
}

cosine_one_to_many <- function(feature_matrix, pathway_vector) {
  P <- as.numeric(pathway_vector)
  denom_p <- sqrt(sum(P^2))

  if (!is.finite(denom_p) || denom_p <= NUM_TOL) {
    return(rep(NA_real_, nrow(feature_matrix)))
  }

  denom_x <- sqrt(rowSums(feature_matrix^2))

  num <- as.numeric(feature_matrix %*% P)

  out <- num / (denom_x * denom_p)
  out[!is.finite(out)] <- NA_real_

  # Protección numérica.
  out <- pmax(-1, pmin(1, out))
  out
}

# --------------------------- localizar corridas -----------------------------

find_latest_selection_run <- function(upstream) {
  base <- file.path(
    upstream,
    "downstream",
    "seleccion_mejor_escenario_train"
  )

  assert(
    dir.exists(base),
    paste("No existe:", base)
  )

  dirs <- list.dirs(
    base,
    recursive = FALSE,
    full.names = TRUE
  )

  dirs <- dirs[file.exists(file.path(
    dirs,
    "SELECTED_SCENARIOS_BY_ANALYSIS.csv"
  ))]

  assert(
    length(dirs) > 0L,
    paste("No hay corridas de selección completas en", base)
  )

  info <- file.info(dirs)

  dirs[
    order(
      info$mtime,
      basename(dirs),
      decreasing = TRUE
    )
  ][1L]
}

find_bidirectional_run <- function(selection_run, family) {
  root <- file.path(
    selection_run,
    family,
    "BIDIRECTIONAL_FRAMEWORK"
  )

  assert(
    dir.exists(root),
    paste("No existe BIDIRECTIONAL_FRAMEWORK para", family)
  )

  latest_file <- file.path(root, "LATEST_RUN.txt")

  if (file.exists(latest_file)) {
    p <- trimws(readLines(latest_file, warn = FALSE)[1L])

    if (
      nzchar(p) &&
      dir.exists(p) &&
      file.exists(file.path(
        p,
        "rds",
        "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"
      ))
    ) {
      return(canonical_path(p, TRUE))
    }
  }

  dirs <- list.dirs(
    root,
    recursive = FALSE,
    full.names = TRUE
  )

  dirs <- dirs[file.exists(file.path(
    dirs,
    "rds",
    "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"
  ))]

  assert(
    length(dirs) > 0L,
    paste("No hay corrida bidireccional completa para", family)
  )

  info <- file.info(dirs)

  canonical_path(
    dirs[
      order(
        info$mtime,
        basename(dirs),
        decreasing = TRUE
      )
    ][1L],
    TRUE
  )
}

find_latest_subrun <- function(
  parent,
  folder,
  required_file
) {
  root <- file.path(parent, folder)

  assert(
    dir.exists(root),
    paste("No existe:", root)
  )

  latest_file <- file.path(root, "LATEST_RUN.txt")

  if (file.exists(latest_file)) {
    p <- trimws(readLines(latest_file, warn = FALSE)[1L])

    if (
      nzchar(p) &&
      dir.exists(p) &&
      file.exists(file.path(p, required_file))
    ) {
      return(canonical_path(p, TRUE))
    }
  }

  dirs <- list.dirs(
    root,
    recursive = FALSE,
    full.names = TRUE
  )

  dirs <- dirs[file.exists(file.path(
    dirs,
    required_file
  ))]

  assert(
    length(dirs) > 0L,
    paste(
      "No hay corrida completa de",
      folder,
      "con",
      required_file
    )
  )

  info <- file.info(dirs)

  canonical_path(
    dirs[
      order(
        info$mtime,
        basename(dirs),
        decreasing = TRUE
      )
    ][1L],
    TRUE
  )
}

find_latent_audit_run <- function(bidir_run) {
  find_latest_subrun(
    parent = bidir_run,
    folder = "LATENT_SUPERVISED_AUDIT",
    required_file = file.path(
      "tables",
      "01_FACTOR_ORIENTATION_AND_SCALING.csv"
    )
  )
}

find_enrichment_run <- function(bidir_run) {
  find_latest_subrun(
    parent = bidir_run,
    folder = "MOFA_NATIVE_ENRICHMENT",
    required_file = file.path(
      "tables",
      "03_MOFA_NATIVE_ENRICHMENT_ALL.csv"
    )
  )
}

find_class_pathway_run <- function(bidir_run) {
  find_latest_subrun(
    parent = bidir_run,
    folder = "CLASS_PATHWAY_BRIDGE",
    required_file = file.path(
      "tables",
      "03_CLASS_PATHWAY_BRIDGE_ALL.csv"
    )
  )
}

# ----------------------- loadings orientados/escalados ----------------------

extract_all_view_scaled_loadings <- function(
  model,
  orientation,
  feature_map
) {
  factors <- as.character(orientation$factor)

  signs <- orientation$orientation_sign
  names(signs) <- orientation$factor

  scales <- orientation$train_scale_after_orientation
  names(scales) <- orientation$factor

  assert(
    all(is.finite(scales)) &&
      all(scales > NUM_TOL),
    "Escalas latentes inválidas"
  )

  view_map <- unique(
    feature_map[, c("view", "view_code"), drop = FALSE]
  )

  assert(
    !anyDuplicated(view_map$view),
    "Una vista tiene más de un view_code"
  )

  rows <- list()
  mats <- list()

  rr <- 1L

  for (v in MOFA2::views_names(model)) {
    W <- MOFA2::get_weights(
      model,
      views = v,
      factors = factors,
      scale = FALSE,
      as.data.frame = FALSE
    )[[1L]]

    W <- as.matrix(W)

    assert(
      all(factors %in% colnames(W)),
      paste("Faltan factores en W para", v)
    )

    W <- W[, factors, drop = FALSE]

    L <- sweep(
      W,
      2,
      signs[factors],
      "*"
    )

    L <- sweep(
      L,
      2,
      scales[factors],
      "*"
    )

    mats[[v]] <- L

    code_hit <- view_map$view_code[
      match(v, view_map$view)
    ]

    code <- if (
      length(code_hit) == 1L &&
      !is.na(code_hit) &&
      nzchar(code_hit)
    ) {
      as.character(code_hit)
    } else {
      make.names(v)
    }

    fmap_v <- feature_map[
      feature_map$view == v,
      ,
      drop = FALSE
    ]

    for (f in factors) {
      feature_original <- rownames(L)

      idx <- match(
        feature_original,
        fmap_v$feature_original
      )

      feature_model <- rep(
        NA_character_,
        length(feature_original)
      )

      if (nrow(fmap_v)) {
        feature_model <- as.character(
          fmap_v$feature_model[idx]
        )
      }

      missing_key <- is.na(feature_model) |
        !nzchar(feature_model)

      feature_model[missing_key] <- paste0(
        code,
        "__",
        feature_original[missing_key]
      )

      tmp <- data.frame(
        view = v,
        view_code = code,
        feature_original = feature_original,
        feature_model = feature_model,
        bidirectional_universe =
          !missing_key,
        factor = f,
        loading_native =
          as.numeric(W[, f]),
        loading_oriented_scaled =
          as.numeric(L[, f]),
        abs_loading_oriented_scaled =
          abs(as.numeric(L[, f])),
        stringsAsFactors = FALSE
      )

      tmp$abs_loading_rank_within_view_factor <-
        rank(
          -tmp$abs_loading_oriented_scaled,
          ties.method = "min"
        )

      tmp$abs_loading_percentile_within_view_factor <-
        empirical_percentile(
          tmp$abs_loading_oriented_scaled
        )

      rows[[rr]] <- tmp
      rr <- rr + 1L
    }
  }

  list(
    long = do.call(rbind, rows),
    matrices = mats
  )
}

add_factor_view_r2 <- function(
  loading_long,
  model
) {
  r2 <- tryCatch(
    MOFA2::get_variance_explained(
      model,
      as.data.frame = TRUE
    )$r2_per_factor,
    error = function(e) NULL
  )

  loading_long$r2_factor_view_raw <- NA_real_

  if (
    is.null(r2) ||
    !is.data.frame(r2) ||
    !all(c("view", "factor", "value") %in% names(r2))
  ) {
    return(loading_long)
  }

  # Si hubiera más de un grupo, promediamos solo como anotación descriptiva.
  r2a <- stats::aggregate(
    value ~ view + factor,
    data = r2,
    FUN = mean,
    na.rm = TRUE
  )

  key_load <- paste(
    loading_long$view,
    loading_long$factor,
    sep = "\r"
  )

  key_r2 <- paste(
    r2a$view,
    r2a$factor,
    sep = "\r"
  )

  loading_long$r2_factor_view_raw <-
    r2a$value[
      match(key_load, key_r2)
    ]

  loading_long
}

# -------------------------- módulos factor-cara -----------------------------

build_factor_side_modules <- function(
  loading_long,
  side_class
) {
  assert(
    all(c(
      "factor",
      "oriented_side",
      "assigned_class_nominal",
      "p_assigned_class"
    ) %in% names(side_class)),
    "Tabla factor-side-clase incompleta"
  )

  rows <- list()
  rr <- 1L

  for (i in seq_len(nrow(side_class))) {
    sc <- side_class[i, , drop = FALSE]

    f <- as.character(sc$factor)
    side <- as.character(sc$oriented_side)
    sgn <- side_sign_value(side)

    z <- loading_long[
      loading_long$factor == f,
      ,
      drop = FALSE
    ]

    z$oriented_side <- side
    z$side_sign <- sgn
    z$side_alignment_loading <-
      sgn * z$loading_oriented_scaled

    z$same_factor_side <-
      is.finite(z$side_alignment_loading) &
      z$side_alignment_loading > NUM_TOL

    # Percentil SOLO entre features que apuntan a esta cara,
    # calculado dentro de cada view x factor.
    z$side_loading_rank_within_view_factor <- NA_integer_
    z$side_loading_percentile_within_view_factor <- NA_real_

    split_idx <- split(
      seq_len(nrow(z)),
      z$view
    )

    for (idx in split_idx) {
      hit <- idx[z$same_factor_side[idx]]

      if (!length(hit)) next

      vals <- z$side_alignment_loading[hit]

      z$side_loading_rank_within_view_factor[hit] <-
        rank(
          -vals,
          ties.method = "min"
        )

      z$side_loading_percentile_within_view_factor[hit] <-
        empirical_percentile(vals)
    }

    # Añadir todas las probabilidades de clase de esa cara.
    for (nm in setdiff(names(sc), c("factor", "oriented_side"))) {
      z[[nm]] <- sc[[nm]][1L]
    }

    # Solo las features que realmente apuntan hacia la cara
    # forman el módulo direccional. No hay mínimo ni cuota por vista.
    z <- z[
      z$same_factor_side,
      ,
      drop = FALSE
    ]

    rows[[rr]] <- z
    rr <- rr + 1L
  }

  if (!length(rows)) {
    return(data.frame())
  }

  do.call(rbind, rows)
}

# --------------------------- posterior beta ---------------------------------

find_factor_fit <- function(
  latent_run,
  factor_name
) {
  candidates <- c(
    file.path(
      latent_run,
      "rds",
      paste0("BRMS_", factor_name, ".rds")
    ),
    file.path(
      latent_run,
      "rds",
      paste0(
        "BRMS_",
        make.names(factor_name),
        ".rds"
      )
    )
  )

  hit <- candidates[file.exists(candidates)]

  assert(
    length(hit) >= 1L,
    paste(
      "No se encontró fit BRMS para",
      factor_name
    )
  )

  hit[1L]
}

extract_group_draws <- function(
  latent_run,
  factors
) {
  out <- list()

  for (f in factors) {
    fit_path <- find_factor_fit(
      latent_run,
      f
    )

    fit <- readRDS(fit_path)

    dr <- posterior::as_draws_df(fit)

    coef_names <- paste0(
      "b_grupo",
      GROUPS
    )

    assert(
      all(coef_names %in% names(dr)),
      paste(
        "Faltan coeficientes de grupo en",
        basename(fit_path)
      )
    )

    m <- as.matrix(
      dr[, coef_names, drop = FALSE]
    )

    colnames(m) <- GROUPS
    out[[f]] <- m
  }

  # Los modelos factor-wise son independientes.
  # Usamos el mismo número de draws por factor y los emparejamos por índice;
  # esto produce una muestra Monte Carlo del producto de marginales bajo la
  # aproximación de independencia ya implícita en el ajuste por factor.
  n_common <- min(
    vapply(
      out,
      nrow,
      integer(1)
    )
  )

  assert(
    n_common >= 100L,
    "Muy pocos draws posteriores comunes"
  )

  out <- lapply(
    out,
    function(x) {
      x[
        seq_len(n_common),
        ,
        drop = FALSE
      ]
    }
  )

  attr(out, "n_common") <- n_common
  out
}

build_beta_draws <- function(
  group_draws,
  factors
) {
  contrast_defs <- list(
    FD_minus_NP = c("FD", "NP"),
    SO_minus_NP = c("SO", "NP"),
    SO_minus_FD = c("SO", "FD")
  )

  out <- list()

  for (cn in names(contrast_defs)) {
    cc <- contrast_defs[[cn]]

    m <- sapply(
      factors,
      function(f) {
        group_draws[[f]][, cc[1L]] -
          group_draws[[f]][, cc[2L]]
      }
    )

    m <- as.matrix(m)
    colnames(m) <- factors
    out[[cn]] <- m
  }

  out
}

# ------------------------ A_jc y descomposición -----------------------------

summarize_class_feature_alignment <- function(
  L_all,
  feature_meta,
  beta_draws
) {
  contrast_rows <- list()
  contribution_rows <- list()
  sign_store <- list()

  rr <- 1L
  cc <- 1L

  factors <- colnames(L_all)

  for (contrast_name in names(beta_draws)) {
    B <- beta_draws[[contrast_name]]

    assert(
      identical(colnames(B), factors),
      "Orden de factores beta != L"
    )

    A <- B %*% t(L_all)

    # Guardar solo signo para soporte conjunto por clase.
    sign_store[[contrast_name]] <- A > 0

    means <- colMeans(A)
    sds <- apply(A, 2, stats::sd)
    qs <- apply(
      A,
      2,
      safe_quantiles
    )

    pgt <- colMeans(A > 0)
    plt <- colMeans(A < 0)

    contrast_rows[[rr]] <- data.frame(
      feature_meta,
      contrast = contrast_name,
      posterior_mean = means,
      posterior_sd = sds,
      q025 = qs["q025", ],
      q50 = qs["q50", ],
      q975 = qs["q975", ],
      p_gt_0 = pgt,
      p_lt_0 = plt,
      p_direction = pmax(pgt, plt),
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L

    # Descomposición A_jfc.
    mean_contrib_matrix <- matrix(
      NA_real_,
      nrow = nrow(L_all),
      ncol = length(factors),
      dimnames = list(
        rownames(L_all),
        factors
      )
    )

    tmp_factor_rows <- list()

    for (f in factors) {
      D <- outer(
        B[, f],
        L_all[, f]
      )

      dm <- colMeans(D)
      ds <- apply(D, 2, stats::sd)
      dq <- apply(
        D,
        2,
        safe_quantiles
      )
      dpgt <- colMeans(D > 0)
      dplt <- colMeans(D < 0)

      mean_contrib_matrix[, f] <- dm

      tmp_factor_rows[[f]] <- data.frame(
        feature_meta,
        contrast = contrast_name,
        factor = f,
        contribution_posterior_mean = dm,
        contribution_posterior_sd = ds,
        contribution_q025 =
          dq["q025", ],
        contribution_q50 =
          dq["q50", ],
        contribution_q975 =
          dq["q975", ],
        contribution_p_gt_0 = dpgt,
        contribution_p_lt_0 = dplt,
        contribution_p_direction =
          pmax(dpgt, dplt),
        stringsAsFactors = FALSE
      )
    }

    denom <- rowSums(
      abs(mean_contrib_matrix)
    )

    for (f in factors) {
      z <- tmp_factor_rows[[f]]

      z$contribution_abs_share <- ifelse(
        denom > NUM_TOL,
        abs(
          mean_contrib_matrix[, f]
        ) / denom,
        NA_real_
      )

      dominant_factor <- apply(
        abs(mean_contrib_matrix),
        1,
        function(x) {
          mx <- max(x, na.rm = TRUE)
          hit <- which(
            abs(x - mx) <= sqrt(.Machine$double.eps)
          )
          paste(
            factors[hit],
            collapse = ";"
          )
        }
      )

      z$dominant_factor_for_contrast <-
        dominant_factor

      contribution_rows[[cc]] <- z
      cc <- cc + 1L
    }

    rm(A)
    gc(FALSE)
  }

  contrasts <- do.call(
    rbind,
    contrast_rows
  )

  contributions <- do.call(
    rbind,
    contribution_rows
  )

  # Probabilidad posterior de una firma "alta en clase" respecto a AMBAS otras.
  # NP high: FD-NP < 0 AND SO-NP < 0
  # FD high: FD-NP > 0 AND SO-FD < 0
  # SO high: SO-NP > 0 AND SO-FD > 0
  assert(
    all(c(
      "FD_minus_NP",
      "SO_minus_NP",
      "SO_minus_FD"
    ) %in% names(sign_store)),
    "Faltan contrastes para soporte class-specific"
  )

  p_np <- colMeans(
    !sign_store$FD_minus_NP &
      !sign_store$SO_minus_NP
  )

  p_fd <- colMeans(
    sign_store$FD_minus_NP &
      !sign_store$SO_minus_FD
  )

  p_so <- colMeans(
    sign_store$SO_minus_NP &
      sign_store$SO_minus_FD
  )

  pmat <- cbind(
    NP = p_np,
    FD = p_fd,
    SO = p_so
  )

  max_idx <- apply(
    pmat,
    1,
    function(x) {
      mx <- max(x, na.rm = TRUE)
      hit <- which(abs(x - mx) <= sqrt(.Machine$double.eps))
      if (length(hit) == 1L) hit else NA_integer_
    }
  )

  sorted <- t(
    apply(
      pmat,
      1,
      sort,
      decreasing = TRUE
    )
  )

  nominal_class <- rep(
    NA_character_,
    nrow(pmat)
  )

  nominal_prob <- rep(
    NA_real_,
    nrow(pmat)
  )

  unique_hit <- !is.na(max_idx)

  nominal_class[unique_hit] <-
    GROUPS[max_idx[unique_hit]]

  nominal_prob[unique_hit] <-
    pmat[
      cbind(
        which(unique_hit),
        max_idx[unique_hit]
      )
    ]

  class_specific <- data.frame(
    feature_meta,
    p_latent_NP_high_vs_both = p_np,
    p_latent_FD_high_vs_both = p_fd,
    p_latent_SO_high_vs_both = p_so,
    latent_class_nominal =
      nominal_class,
    p_latent_class_nominal =
      nominal_prob,
    latent_class_probability_margin =
      sorted[, 1L] - sorted[, 2L],
    latent_class_assignment_tied =
      is.na(max_idx),
    stringsAsFactors = FALSE
  )

  list(
    contrasts = contrasts,
    contributions = contributions,
    class_specific = class_specific
  )
}

# --------------------------- matrices L globales ----------------------------

make_feature_factor_matrix <- function(
  loading_long,
  factors
) {
  keys <- unique(
    loading_long[, c(
      "view",
      "view_code",
      "feature_original",
      "feature_model",
      "bidirectional_universe"
    )]
  )

  feature_key <- paste(
    keys$view,
    keys$feature_original,
    sep = "\r"
  )

  L <- matrix(
    NA_real_,
    nrow = nrow(keys),
    ncol = length(factors),
    dimnames = list(
      feature_key,
      factors
    )
  )

  for (f in factors) {
    z <- loading_long[
      loading_long$factor == f,
      ,
      drop = FALSE
    ]

    key_z <- paste(
      z$view,
      z$feature_original,
      sep = "\r"
    )

    idx <- match(
      feature_key,
      key_z
    )

    assert(
      all(!is.na(idx)),
      paste("Faltan loadings para", f)
    )

    L[, f] <- z$loading_oriented_scaled[idx]
  }

  list(
    matrix = L,
    meta = keys
  )
}

# -------------------------- perfiles de pathways ----------------------------

load_feature_set_cache <- function(
  enrichment_run,
  libraries
) {
  out <- list()

  for (lib in unique(libraries)) {
    path <- file.path(
      enrichment_run,
      "rds",
      paste0(
        "NATIVE_ENRICHMENT_",
        lib,
        ".rds"
      )
    )

    if (!file.exists(path)) {
      out[[lib]] <- NULL
      next
    }

    out[[lib]] <- readRDS(path)
  }

  out
}

build_pathway_profiles <- function(
  primary_pathways,
  fs_cache,
  L_tx,
  factors
) {
  if (!nrow(primary_pathways)) {
    return(list(
      profiles = data.frame(),
      vectors = list(),
      members = list()
    ))
  }

  rows <- list()
  vectors <- list()
  members_out <- list()

  rr <- 1L

  for (i in seq_len(nrow(primary_pathways))) {
    p <- primary_pathways[i, , drop = FALSE]

    lib <- as.character(p$library)
    pathway <- as.character(p$pathway)
    f_anchor <- as.character(p$factor)
    side <- as.character(p$oriented_side)
    sgn <- side_sign_value(side)

    obj <- fs_cache[[lib]]

    if (
      is.null(obj) ||
      !"feature_set_matrix" %in% names(obj)
    ) {
      next
    }

    membership <- as.matrix(
      obj$feature_set_matrix
    )

    if (!pathway %in% rownames(membership)) {
      next
    }

    members <- colnames(membership)[
      membership[pathway, ] == 1
    ]

    members <- intersect(
      members,
      rownames(L_tx)
    )

    if (!length(members)) next

    anchor_values <-
      sgn * L_tx[members, f_anchor]

    side_drivers <- members[
      is.finite(anchor_values) &
        anchor_values > NUM_TOL
    ]

    if (!length(side_drivers)) next

    weights <- sgn *
      L_tx[side_drivers, f_anchor]

    weights <- pmax(
      as.numeric(weights),
      0
    )

    if (
      !all(is.finite(weights)) ||
      sum(weights) <= NUM_TOL
    ) {
      next
    }

    P <- colSums(
      sweep(
        L_tx[
          side_drivers,
          factors,
          drop = FALSE
        ],
        1,
        weights,
        "*"
      )
    ) / sum(weights)

    vectors[[rr]] <- P
    names(vectors)[rr] <- paste(
      lib,
      pathway,
      f_anchor,
      side,
      sep = "||"
    )

    members_out[[rr]] <- data.frame(
      library = lib,
      pathway = pathway,
      factor = f_anchor,
      oriented_side = side,
      feature_original =
        side_drivers,
      anchor_weight =
        weights,
      stringsAsFactors = FALSE
    )

    row <- data.frame(
      family = if ("family" %in% names(p)) p$family else NA_character_,
      library = lib,
      pathway = pathway,
      factor = f_anchor,
      oriented_side = side,
      assigned_class_nominal =
        if ("assigned_class_nominal" %in% names(p)) p$assigned_class_nominal else NA_character_,
      p_assigned_class =
        if ("p_assigned_class" %in% names(p)) p$p_assigned_class else NA_real_,
      q_both_sides_within_library_factor =
        if ("q_both_sides_within_library_factor" %in% names(p)) {
          p$q_both_sides_within_library_factor
        } else {
          NA_real_
        },
      n_pathway_members_in_model =
        length(members),
      n_side_driver_genes =
        length(side_drivers),
      stringsAsFactors = FALSE
    )

    for (f in factors) {
      row[[paste0(
        "pathway_profile_",
        f
      )]] <- as.numeric(P[f])
    }

    rows[[rr]] <- row
    rr <- rr + 1L
  }

  if (!length(rows)) {
    return(list(
      profiles = data.frame(),
      vectors = list(),
      members = list()
    ))
  }

  # Reindex names consistently.
  vectors <- vectors[
    seq_along(rows)
  ]
  members_out <- members_out[
    seq_along(rows)
  ]

  list(
    profiles = do.call(rbind, rows),
    vectors = vectors,
    members = members_out
  )
}

build_pathway_feature_links <- function(
  pathway_profiles,
  pathway_vectors,
  pathway_members,
  module_features,
  L_all,
  feature_meta,
  factors
) {
  if (
    !nrow(pathway_profiles) ||
    !nrow(module_features)
  ) {
    return(data.frame())
  }

  feature_key <- paste(
    feature_meta$view,
    feature_meta$feature_original,
    sep = "\r"
  )

  rows <- list()
  rr <- 1L

  for (i in seq_len(nrow(pathway_profiles))) {
    p <- pathway_profiles[i, , drop = FALSE]

    f_anchor <- as.character(p$factor)
    side <- as.character(p$oriented_side)

    mf <- module_features[
      module_features$factor == f_anchor &
        module_features$oriented_side == side,
      ,
      drop = FALSE
    ]

    if (!nrow(mf)) next

    mf_key <- paste(
      mf$view,
      mf$feature_original,
      sep = "\r"
    )

    idx <- match(
      mf_key,
      feature_key
    )

    keep <- !is.na(idx)

    if (!any(keep)) next

    mf <- mf[keep, , drop = FALSE]
    idx <- idx[keep]

    P <- pathway_vectors[[i]]

    cosv <- cosine_one_to_many(
      L_all[
        idx,
        factors,
        drop = FALSE
      ],
      P[factors]
    )

    drivers <- pathway_members[[i]]$feature_original

    z <- data.frame(
      family =
        if ("family" %in% names(p)) p$family else NA_character_,
      library = p$library,
      pathway = p$pathway,
      factor = f_anchor,
      oriented_side = side,
      assigned_class_nominal =
        p$assigned_class_nominal,
      p_assigned_class =
        p$p_assigned_class,
      q_both_sides_within_library_factor =
        p$q_both_sides_within_library_factor,
      view = mf$view,
      view_code = mf$view_code,
      feature_original =
        mf$feature_original,
      feature_model =
        mf$feature_model,
      pathway_feature_cosine = cosv,
      feature_is_transcriptomic_pathway_driver =
        mf$view_code == "tx" &
        mf$feature_original %in% drivers,
      stringsAsFactors = FALSE
    )

    rows[[rr]] <- z
    rr <- rr + 1L
  }

  if (!length(rows)) {
    return(data.frame())
  }

  do.call(rbind, rows)
}

top_links_per_pathway_view <- function(
  links
) {
  if (!nrow(links)) return(data.frame())

  key <- interaction(
    links$library,
    links$pathway,
    links$factor,
    links$oriented_side,
    links$view,
    drop = TRUE,
    lex.order = TRUE
  )

  idx <- unlist(
    lapply(
      split(seq_len(nrow(links)), key),
      function(ii) {
        z <- links$pathway_feature_cosine[ii]

        ord <- order(
          z,
          decreasing = TRUE,
          na.last = NA
        )

        if (!length(ord)) return(integer())

        ii[
          ord[
            seq_len(
              min(
                TOP_LINKS_PER_PATHWAY_VIEW,
                length(ord)
              )
            )
          ]
        ]
      }
    ),
    use.names = FALSE
  )

  links[idx, , drop = FALSE]
}

# --------------------- soporte bidireccional TRAIN --------------------------

extract_bidirectional_support <- function(
  bidir_object
) {
  assert(
    "final_train" %in% names(bidir_object),
    "BIDIRECTIONAL_FRAMEWORK_OBJECT no contiene final_train"
  )

  x <- bidir_object$final_train

  wanted <- c(
    "feature_model",
    "feature_original",
    "view",
    "view_code",
    "GeneSymbol",
    "EntrezGeneID",
    "EnsemblID",
    "rec_evidence_grade",
    "obs_evidence_grade",
    "rec_bidirectional_pattern",
    "obs_bidirectional_pattern",
    "direction_concordance_n",
    "direction_concordance_prop",
    "cross_signal_support",
    "priority_tier",
    "biomarker_candidate_primary",
    "biomarker_candidate_exploratory",
    "high_confidence_AA_same_pattern",
    "high_confidence_class_signature",
    "high_confidence_class",
    "rec_q_global",
    "obs_q_global",
    "rec_MI_norm",
    "obs_MI_norm"
  )

  wanted <- intersect(
    wanted,
    names(x)
  )

  x[, wanted, drop = FALSE]
}

# ---------------------- matriz final biomarcadores --------------------------

attach_top_pathway_to_modules <- function(
  module_features,
  links
) {
  module_features$top_primary_pathway <- NA_character_
  module_features$top_primary_pathway_library <- NA_character_
  module_features$top_primary_pathway_cosine <- NA_real_
  module_features$n_primary_pathways_in_module <- 0L

  if (!nrow(links)) {
    return(module_features)
  }

  key_mod <- paste(
    module_features$factor,
    module_features$oriented_side,
    module_features$feature_model,
    sep = "\r"
  )

  key_link <- paste(
    links$factor,
    links$oriented_side,
    links$feature_model,
    sep = "\r"
  )

  split_links <- split(
    seq_len(nrow(links)),
    key_link
  )

  for (i in seq_len(nrow(module_features))) {
    k <- key_mod[i]
    idx <- split_links[[k]]

    if (is.null(idx) || !length(idx)) next

    z <- links[idx, , drop = FALSE]

    module_features$n_primary_pathways_in_module[i] <-
      length(unique(
        paste(
          z$library,
          z$pathway,
          sep = "||"
        )
      ))

    valid <- is.finite(
      z$pathway_feature_cosine
    )

    if (!any(valid)) next

    zv <- z[valid, , drop = FALSE]

    j <- which.max(
      zv$pathway_feature_cosine
    )

    module_features$top_primary_pathway[i] <-
      zv$pathway[j]

    module_features$top_primary_pathway_library[i] <-
      zv$library[j]

    module_features$top_primary_pathway_cosine[i] <-
      zv$pathway_feature_cosine[j]
  }

  module_features
}

attach_class_specific_support <- function(
  module_features,
  class_specific
) {
  idx <- match(
    module_features$feature_model,
    class_specific$feature_model
  )

  support_cols <- setdiff(
    names(class_specific),
    c(
      "view",
      "view_code",
      "feature_original",
      "feature_model",
      "bidirectional_universe"
    )
  )

  for (cc in support_cols) {
    module_features[[cc]] <-
      class_specific[[cc]][idx]
  }

  module_features$latent_support_for_module_class <- NA_real_

  module_features$latent_support_for_module_class[
    module_features$assigned_class_nominal == "NP"
  ] <- module_features$p_latent_NP_high_vs_both[
    module_features$assigned_class_nominal == "NP"
  ]

  module_features$latent_support_for_module_class[
    module_features$assigned_class_nominal == "FD"
  ] <- module_features$p_latent_FD_high_vs_both[
    module_features$assigned_class_nominal == "FD"
  ]

  module_features$latent_support_for_module_class[
    module_features$assigned_class_nominal == "SO"
  ] <- module_features$p_latent_SO_high_vs_both[
    module_features$assigned_class_nominal == "SO"
  ]

  module_features$module_class_matches_latent_nominal <-
    !is.na(module_features$assigned_class_nominal) &
    !is.na(module_features$latent_class_nominal) &
    module_features$assigned_class_nominal ==
      module_features$latent_class_nominal

  module_features
}

attach_bidirectional_support <- function(
  modules,
  bidir_support
) {
  idx <- match(
    modules$feature_model,
    bidir_support$feature_model
  )

  add_cols <- setdiff(
    names(bidir_support),
    c(
      "feature_model",
      "feature_original",
      "view",
      "view_code"
    )
  )

  for (cc in add_cols) {
    modules[[cc]] <-
      bidir_support[[cc]][idx]
  }

  expected_pattern <- paste0(
    modules$assigned_class_nominal,
    "_like_high"
  )

  modules$rec_pattern_matches_module_class <-
    !is.na(modules$rec_bidirectional_pattern) &
    modules$rec_bidirectional_pattern ==
      expected_pattern

  modules$obs_pattern_matches_module_class <-
    !is.na(modules$obs_bidirectional_pattern) &
    modules$obs_bidirectional_pattern ==
      expected_pattern

  modules
}

# ----------------------------- una familia ----------------------------------

run_family <- function(
  selection_run,
  family,
  stamp
) {
  message("\n============================================================")
  message("CROSS-OMIC BRIDGE | ", family)
  message("============================================================")

  bidir_run <- find_bidirectional_run(
    selection_run,
    family
  )

  latent_run <- find_latent_audit_run(
    bidir_run
  )

  enrichment_run <- find_enrichment_run(
    bidir_run
  )

  class_pathway_run <- find_class_pathway_run(
    bidir_run
  )

  bidir_object_path <- file.path(
    bidir_run,
    "rds",
    "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"
  )

  bidir_object <- readRDS(
    bidir_object_path
  )

  assert(
    all(c(
      "model_path",
      "feature_map",
      "final_train"
    ) %in% names(bidir_object)),
    "BIDIRECTIONAL_FRAMEWORK_OBJECT incompleto"
  )

  model <- readRDS(
    canonical_path(
      bidir_object$model_path,
      TRUE
    )
  )

  orientation_file <- file.path(
    latent_run,
    "tables",
    "01_FACTOR_ORIENTATION_AND_SCALING.csv"
  )

  orientation <- utils::read.csv(
    orientation_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_orientation <- c(
    "factor",
    "orientation_sign",
    "train_scale_after_orientation"
  )

  assert(
    all(required_orientation %in% names(orientation)),
    paste(
      "Faltan columnas en orientación:",
      paste(
        setdiff(
          required_orientation,
          names(orientation)
        ),
        collapse = ", "
      )
    )
  )

  factors <- as.character(
    orientation$factor
  )

  # -------------------------------------------------------------------------
  # 1) Loadings L_jf en la misma métrica de los beta.
  # -------------------------------------------------------------------------
  loadings <- extract_all_view_scaled_loadings(
    model = model,
    orientation = orientation,
    feature_map = bidir_object$feature_map
  )

  loading_long <- add_factor_view_r2(
    loadings$long,
    model
  )

  # -------------------------------------------------------------------------
  # 2) Cara del factor -> clase.
  # -------------------------------------------------------------------------
  side_class_file <- file.path(
    class_pathway_run,
    "tables",
    "01_FACTOR_SIDE_CLASS_POSTERIOR.csv"
  )

  side_class <- utils::read.csv(
    side_class_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  module_features <- build_factor_side_modules(
    loading_long = loading_long,
    side_class = side_class
  )

  # -------------------------------------------------------------------------
  # 3) Matriz feature x factor L.
  # -------------------------------------------------------------------------
  ff <- make_feature_factor_matrix(
    loading_long = loading_long,
    factors = factors
  )

  L_all <- ff$matrix
  feature_meta <- ff$meta

  # TX específicamente para construir perfiles de pathways.
  tx_view <- unique(
    loading_long$view[
      loading_long$view_code == "tx"
    ]
  )

  assert(
    length(tx_view) == 1L,
    paste(
      "Se esperaba una única vista tx; encontradas:",
      paste(tx_view, collapse = ", ")
    )
  )

  tx_keys <- paste(
    tx_view,
    rownames(loadings$matrices[[tx_view]]),
    sep = "\r"
  )

  idx_tx <- match(
    tx_keys,
    rownames(L_all)
  )

  assert(
    all(!is.na(idx_tx)),
    "No se pudo reconstruir L transcriptómico"
  )

  L_tx <- L_all[
    idx_tx,
    factors,
    drop = FALSE
  ]

  rownames(L_tx) <-
    rownames(
      loadings$matrices[[tx_view]]
    )

  # -------------------------------------------------------------------------
  # 4) Pathways primarios -> perfiles latentes.
  # -------------------------------------------------------------------------
  bridge_file <- file.path(
    class_pathway_run,
    "tables",
    "03_CLASS_PATHWAY_BRIDGE_ALL.csv"
  )

  class_pathway <- utils::read.csv(
    bridge_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  assert(
    "pathway_evidence" %in% names(class_pathway),
    "El class-pathway bridge no contiene pathway_evidence"
  )

  primary_pathways <- class_pathway[
    class_pathway$pathway_evidence ==
      PRIMARY_PATHWAY_LABEL,
    ,
    drop = FALSE
  ]

  fs_cache <- load_feature_set_cache(
    enrichment_run = enrichment_run,
    libraries = primary_pathways$library
  )

  pathway_obj <- build_pathway_profiles(
    primary_pathways = primary_pathways,
    fs_cache = fs_cache,
    L_tx = L_tx,
    factors = factors
  )

  pathway_profiles <- pathway_obj$profiles

  pathway_links <- build_pathway_feature_links(
    pathway_profiles =
      pathway_profiles,
    pathway_vectors =
      pathway_obj$vectors,
    pathway_members =
      pathway_obj$members,
    module_features =
      module_features,
    L_all = L_all,
    feature_meta =
      feature_meta,
    factors = factors
  )

  pathway_links_top <-
    top_links_per_pathway_view(
      pathway_links
    )

  # -------------------------------------------------------------------------
  # 5) Clase -> feature: A_jc y A_jfc.
  # -------------------------------------------------------------------------
  group_draws <- extract_group_draws(
    latent_run = latent_run,
    factors = factors
  )

  beta_draws <- build_beta_draws(
    group_draws = group_draws,
    factors = factors
  )

  alignment <- summarize_class_feature_alignment(
    L_all = L_all,
    feature_meta = feature_meta,
    beta_draws = beta_draws
  )

  # -------------------------------------------------------------------------
  # 6) Matriz integrada para descubrimiento/priorización de biomarcadores.
  # -------------------------------------------------------------------------
  biomarker_matrix <- attach_top_pathway_to_modules(
    module_features =
      module_features,
    links = pathway_links
  )

  biomarker_matrix <- attach_class_specific_support(
    module_features =
      biomarker_matrix,
    class_specific =
      alignment$class_specific
  )

  bidir_support <- extract_bidirectional_support(
    bidir_object
  )

  biomarker_matrix <- attach_bidirectional_support(
    modules = biomarker_matrix,
    bidir_support = bidir_support
  )

  # Orden solo para lectura; NO constituye score compuesto ni selección.
  ord <- order(
    biomarker_matrix$assigned_class_nominal,
    biomarker_matrix$factor,
    biomarker_matrix$view_code,
    -biomarker_matrix$latent_support_for_module_class,
    -biomarker_matrix$side_loading_percentile_within_view_factor,
    -biomarker_matrix$top_primary_pathway_cosine,
    na.last = TRUE
  )

  biomarker_matrix <-
    biomarker_matrix[
      ord,
      ,
      drop = FALSE
    ]

  # Resumen de composición por módulo/vista.
  if (nrow(module_features)) {
    module_view_summary <- stats::aggregate(
      feature_model ~
        factor +
        oriented_side +
        assigned_class_nominal +
        view +
        view_code,
      data = module_features,
      FUN = length
    )

    names(module_view_summary)[
      names(module_view_summary) ==
        "feature_model"
    ] <- "n_same_side_features"
  } else {
    module_view_summary <- data.frame()
  }

  # -------------------------------------------------------------------------
  # Salidas.
  # -------------------------------------------------------------------------
  out <- file.path(
    bidir_run,
    "CROSSOMIC_BRIDGE",
    stamp
  )

  dir.create(
    file.path(out, "tables"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    file.path(out, "rds"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  input_audit <- data.frame(
    family = family,
    bidirectional_run = bidir_run,
    latent_audit_run = latent_run,
    enrichment_run = enrichment_run,
    class_pathway_run =
      class_pathway_run,
    model_path =
      canonical_path(
        bidir_object$model_path,
        TRUE
      ),
    n_factors = length(factors),
    n_features_all_views =
      nrow(feature_meta),
    n_primary_pathways =
      nrow(primary_pathways),
    n_posterior_draws_product_marginals =
      attr(
        group_draws,
        "n_common"
      ),
    posterior_factor_independence_assumed =
      TRUE,
    TEST_used = FALSE,
    new_model_fitted = FALSE,
    hard_minimum_features_per_view =
      FALSE,
    stringsAsFactors = FALSE
  )

  write_csv(
    input_audit,
    file.path(
      out,
      "tables",
      "00_INPUT_AUDIT.csv"
    )
  )

  write_csv(
    loading_long,
    file.path(
      out,
      "tables",
      "01_ORIENTED_SCALED_LOADINGS_ALL_VIEWS.csv"
    )
  )

  write_csv(
    module_features,
    file.path(
      out,
      "tables",
      "02_FACTOR_SIDE_CROSSOMIC_MODULE_FEATURES.csv"
    )
  )

  write_csv(
    module_view_summary,
    file.path(
      out,
      "tables",
      "02B_FACTOR_SIDE_MODULE_VIEW_SUMMARY.csv"
    )
  )

  write_csv(
    pathway_profiles,
    file.path(
      out,
      "tables",
      "03_PRIMARY_PATHWAY_LATENT_PROFILES.csv"
    )
  )

  pathway_driver_table <- if (
    length(pathway_obj$members)
  ) {
    do.call(
      rbind,
      pathway_obj$members
    )
  } else {
    data.frame()
  }

  write_csv(
    pathway_driver_table,
    file.path(
      out,
      "tables",
      "03B_PRIMARY_PATHWAY_SIDE_DRIVER_GENES.csv"
    )
  )

  write_csv(
    pathway_links,
    file.path(
      out,
      "tables",
      "04_PATHWAY_FEATURE_CROSSOMIC_LINKS_ALL.csv"
    )
  )

  write_csv(
    pathway_links_top,
    file.path(
      out,
      "tables",
      "04B_PATHWAY_FEATURE_TOP_LINKS_PER_VIEW.csv"
    )
  )

  write_csv(
    alignment$contrasts,
    file.path(
      out,
      "tables",
      "05_CLASS_FEATURE_LATENT_CONTRASTS.csv"
    )
  )

  write_csv(
    alignment$contributions,
    file.path(
      out,
      "tables",
      "06_CLASS_FEATURE_FACTOR_CONTRIBUTIONS.csv"
    )
  )

  write_csv(
    alignment$class_specific,
    file.path(
      out,
      "tables",
      "07_CLASS_SPECIFIC_LATENT_FEATURE_SUPPORT.csv"
    )
  )

  write_csv(
    biomarker_matrix,
    file.path(
      out,
      "tables",
      "08_BIOMARKER_EVIDENCE_MATRIX.csv"
    )
  )

  saveRDS(
    list(
      family = family,
      input_audit = input_audit,
      orientation = orientation,
      loadings_long = loading_long,
      L_all = L_all,
      feature_meta = feature_meta,
      factor_side_modules = module_features,
      module_view_summary =
        module_view_summary,
      primary_pathways =
        primary_pathways,
      pathway_profiles =
        pathway_profiles,
      pathway_driver_genes =
        pathway_driver_table,
      pathway_feature_links =
        pathway_links,
      class_feature_contrasts =
        alignment$contrasts,
      factor_contributions =
        alignment$contributions,
      class_specific_latent_support =
        alignment$class_specific,
      biomarker_evidence_matrix =
        biomarker_matrix,
      bidirectional_support =
        bidir_support,
      parameters = list(
        primary_pathway_label =
          PRIMARY_PATHWAY_LABEL,
        top_links_per_pathway_view =
          TOP_LINKS_PER_PATHWAY_VIEW,
        loading_transform =
          "L_jf = orientation_sign_f * train_scale_after_orientation_f * W_jf",
        test_used = FALSE
      )
    ),
    file.path(
      out,
      "rds",
      "CROSSOMIC_BRIDGE_OBJECT.rds"
    ),
    version = 3
  )

  readme <- c(
    paste(
      "METATEST |",
      family,
      "| CROSS-OMIC BRIDGE"
    ),
    "",
    "OBJETIVO:",
    "- Generar y priorizar biomarcadores candidatos dentro de la estructura",
    "  multi-ómica aprendida por MOFA, conservando interpretación biológica.",
    "- NO constituye validación final de biomarcadores.",
    "",
    "ESCALA LATENTE:",
    "- La auditoría bayesiana modeló Z orientado y estandarizado.",
    "- Por consistencia, los loadings usados aquí son:",
    "    L_jf = orientation_sign_f * train_scale_after_orientation_f * W_jf",
    "- El centro de Z solo induce una constante/intercepto en la reconstrucción",
    "  y no entra en los contrastes de clase.",
    "",
    "MÓDULOS FACTOR-CARA:",
    "- Una feature se vincula a una cara si side_sign * L_jf > 0.",
    "- No hay cuota ni mínimo de features por vista.",
    "- Se reportan rangos/percentiles dentro de cada view x factor, no se",
    "  comparan directamente magnitudes crudas entre ómicas.",
    "- r2_factor_view_raw se conserva como contexto del peso del factor en la vista.",
    "",
    "PATHWAY <-> FEATURE CROSS-ÓMICA:",
    "- Los pathways proceden SOLO de PRIMARY_TWO_SIDE_FDR del class-pathway bridge.",
    "- Cada pathway se representa mediante un vector latente transcriptómico",
    "  construido con todos sus genes miembros que apuntan hacia la cara",
    "  enriquecida y ponderados por su carga en el factor ancla.",
    "- cosine(feature, pathway) mide semejanza del perfil latente.",
    "- Para PR/ME/CL esto NO significa pertenencia bioquímica al pathway;",
    "  significa coherencia con el mismo programa latente multi-ómico.",
    "",
    "CLASE <-> FEATURE:",
    "- A_jc = sum_f L_jf * beta_fc.",
    "- A_jfc = L_jf * beta_fc permite identificar qué factor conduce el contraste.",
    "- La incertidumbre se propaga desde los fits BRMS ya existentes.",
    "- Como los factores se ajustaron por separado, la combinación usa el producto",
    "  de posteriores marginales: aproximación de independencia factor-wise.",
    "",
    "SOPORTE CLASS-SPECIFIC:",
    "- p_latent_NP_high_vs_both = P(FD-NP < 0 Y SO-NP < 0).",
    "- p_latent_FD_high_vs_both = P(FD-NP > 0 Y SO-FD < 0).",
    "- p_latent_SO_high_vs_both = P(SO-NP > 0 Y SO-FD > 0).",
    "- Son probabilidades de soporte latente; no son validación observada.",
    "",
    "FRAMEWORK BIDIRECCIONAL:",
    "- Se añade final_train como evidencia ortogonal TRAIN-only.",
    "- evidence_grade, pattern, priority_tier y high-confidence NO son filtros",
    "  obligatorios para pertenecer al módulo.",
    "- Esto evita repetir el cuello de botella A/A.",
    "",
    "TEST:",
    "- No participa.",
    "",
    "SALIDAS PRINCIPALES:",
    "- 01_ORIENTED_SCALED_LOADINGS_ALL_VIEWS.csv",
    "- 02_FACTOR_SIDE_CROSSOMIC_MODULE_FEATURES.csv",
    "- 03_PRIMARY_PATHWAY_LATENT_PROFILES.csv",
    "- 04_PATHWAY_FEATURE_CROSSOMIC_LINKS_ALL.csv",
    "- 05_CLASS_FEATURE_LATENT_CONTRASTS.csv",
    "- 06_CLASS_FEATURE_FACTOR_CONTRIBUTIONS.csv",
    "- 07_CLASS_SPECIFIC_LATENT_FEATURE_SUPPORT.csv",
    "- 08_BIOMARKER_EVIDENCE_MATRIX.csv",
    "",
    "INTERPRETACIÓN DE 08_BIOMARKER_EVIDENCE_MATRIX:",
    "- module class: qué clase ocupa la cara del factor.",
    "- side loading percentile: cuánto participa la feature en esa cara dentro",
    "  de su propia vista.",
    "- latent_support_for_module_class: soporte posterior de que la feature sea",
    "  alta en esa clase frente a las otras dos dentro del modelo latente.",
    "- top pathway cosine: coherencia de la feature con el programa transcriptómico",
    "  de un pathway de esa misma cara.",
    "- bidirectional fields: corroboración independiente de la firma clase-feature.",
    "",
    "SIGUIENTE PASO:",
    "- Corroboración observada TRAIN: score transcriptómico de pathway vs PR/ME/CL",
    "  mediante asociaciones robustas/permutación, sin usarla para construir",
    "  retrospectivamente los módulos.",
    "- TEST, si se usa después, solo como auditoría exploratoria de dirección."
  )

  writeLines(
    readme,
    file.path(
      out,
      "READ_ME_FIRST.txt"
    ),
    useBytes = TRUE
  )

  root_latest <- file.path(
    bidir_run,
    "CROSSOMIC_BRIDGE"
  )

  dir.create(
    root_latest,
    recursive = TRUE,
    showWarnings = FALSE
  )

  writeLines(
    out,
    file.path(
      root_latest,
      "LATEST_RUN.txt"
    )
  )

  data.frame(
    family = family,
    n_features_all_views =
      nrow(feature_meta),
    n_module_feature_rows =
      nrow(module_features),
    n_primary_pathways =
      nrow(primary_pathways),
    n_pathway_feature_links =
      nrow(pathway_links),
    n_biomarker_evidence_rows =
      nrow(biomarker_matrix),
    output = out,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------- MAIN ------------------------------------

main <- function() {
  assert(
    requireNamespace("MOFA2", quietly = TRUE),
    "Se requiere MOFA2"
  )

  assert(
    requireNamespace("posterior", quietly = TRUE),
    "Se requiere posterior"
  )

  assert(
    requireNamespace("brms", quietly = TRUE),
    "Se requiere brms para leer correctamente los fits brmsfit"
  )

  upstream <- canonical_path(
    UPSTREAM_DIR,
    TRUE
  )

  selection_run <- if (
    nzchar(SELECTION_RUN)
  ) {
    canonical_path(
      SELECTION_RUN,
      TRUE
    )
  } else {
    find_latest_selection_run(
      upstream
    )
  }

  stamp <- paste0(
    format(
      Sys.time(),
      "%Y%m%dT%H%M%SZ",
      tz = "UTC"
    ),
    "_pid",
    Sys.getpid()
  )

  statuses <- list()

  for (family in FAMILIES) {
    statuses[[family]] <- tryCatch(
      run_family(
        selection_run = selection_run,
        family = family,
        stamp = stamp
      ),
      error = function(e) {
        data.frame(
          family = family,
          n_features_all_views =
            NA_integer_,
          n_module_feature_rows =
            NA_integer_,
          n_primary_pathways =
            NA_integer_,
          n_pathway_feature_links =
            NA_integer_,
          n_biomarker_evidence_rows =
            NA_integer_,
          output = NA_character_,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  }

  status <- do.call(
    rbind,
    statuses
  )

  status_file <- file.path(
    selection_run,
    paste0(
      "CROSSOMIC_BRIDGE_RUN_",
      stamp,
      ".csv"
    )
  )

  write_csv(
    status,
    status_file
  )

  print(
    status,
    row.names = FALSE
  )

  if (
    any(
      !is.na(status$error) &
        nzchar(status$error)
    )
  ) {
    stop(
      "Alguna familia falló. Revisar la tabla de estado.",
      call. = FALSE
    )
  }

  invisible(status)
}

if (sys.nframe() == 0L) main()

#!/usr/bin/env Rscript
# =============================================================================
# METATEST — PASO 6A
# Corroboración observada TRAIN de módulos cross-ómicos
# Version 0.1.1
#
# OBJETIVO
#   Corroborar en datos OBSERVADOS TRAIN, sin redefinir ni seleccionar módulos:
#
#       A) pathway transcriptómico observado <-> clase
#       B) pathway transcriptómico observado <-> PR / ME / CL observados
#
#   Todo lo que se prueba aquí fue definido previamente por:
#
#       LATENT_SUPERVISED_AUDIT
#       -> MOFA_NATIVE_ENRICHMENT
#       -> CLASS_PATHWAY_BRIDGE
#       -> CROSSOMIC_BRIDGE
#
# PRINCIPIOS
#   - TRAIN solamente.
#   - TEST NO se abre ni se usa.
#   - No se reentrena MOFA.
#   - No se recalculan pathways ni módulos.
#   - No se usan asociaciones observadas para rescatar/eliminar retrospectivamente
#     features del módulo.
#   - La evidencia observada se añade como una capa adicional.
#
# SCORE OBSERVADO DE PATHWAY
#   Para un pathway q, factor f y cara s ya congelados:
#
#       Xstd_ig = (X_ig - mean_TRAIN_g) / sd_TRAIN_g
#
#       a_g = anchor_weight_g > 0
#       w_g = a_g / sum_h a_h
#
#       Score_iq = sum_g w_g * Xstd_ig
#
#   Los genes y anchor_weight proceden del PASO 5.
#   No se vuelven a elegir usando las clases.
#
# CORROBORACIÓN PATHWAY <-> CLASE
#   El PASO 4 asignó una clase a la cara del factor.
#   Se pregunta si el score observado del pathway es mayor en esa clase que en
#   las otras dos. Estadístico principal:
#
#       T_q = median(score | clase asignada)
#             - max[median(score | cada otra clase)]
#
#   p por permutación TRAIN, preservando los tamaños NP/FD/SO.
#
# CORROBORACIÓN PATHWAY <-> FEATURE
#   Solo PR / ME / CL. TX no se usa como "corroboración cross-ómica" porque el
#   score del pathway se construye a partir de transcriptómica.
#
#   Primario descriptivo:
#       rho_total = Spearman(score_pathway, feature_observada)
#
#   Dirección esperada:
#       sign_expected = sign(cosine_latente_pathway_feature)
#
#   Se reporta concordancia de signo y p de permutación direccional.
#
#   Secundario mecanístico:
#       rho_within_class
#
#   calculado tras centrar pathway y feature por la mediana de su clase, con
#   permutaciones restringidas dentro de NP/FD/SO. Esto pregunta si existe
#   covariación adicional más allá del desplazamiento entre clases.
#
# ROBUSTEZ
#   - Spearman
#   - permutación
#   - jackknife leave-one-out
#   - BH como auditoría de multiplicidad; NO como filtro duro.
#
# SALIDA PARA PASO 6B
#   Se congela explícitamente:
#   - genes usados por pathway
#   - media/SD TRAIN por gen
#   - pesos normalizados
#   - relaciones pathway-feature ya preespecificadas
#   para proyectarlas después a TEST sin reselección.
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

NPERM <- as.integer(
  Sys.getenv("OBSERVED_TRAIN_NPERM", unset = "5000")
)

SEED <- as.integer(
  Sys.getenv("OBSERVED_TRAIN_SEED", unset = "20260922")
)

PADJ_METHOD <- Sys.getenv(
  "OBSERVED_TRAIN_PADJ",
  unset = "BH"
)

ALPHA_AUDIT <- as.numeric(
  Sys.getenv("OBSERVED_TRAIN_ALPHA", unset = "0.10")
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

pathway_key <- function(
  library,
  pathway,
  factor,
  oriented_side
) {
  paste(
    library,
    pathway,
    factor,
    oriented_side,
    sep = "||"
  )
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)

  if (sum(ok) < 4L) return(NA_real_)

  x <- x[ok]
  y <- y[ok]

  if (
    stats::sd(x) <= NUM_TOL ||
    stats::sd(y) <= NUM_TOL
  ) {
    return(NA_real_)
  }

  suppressWarnings(
    stats::cor(
      x,
      y,
      method = "spearman",
      use = "complete.obs"
    )
  )
}

group_median_center <- function(x, group) {
  out <- rep(NA_real_, length(x))

  for (g in GROUPS) {
    idx <- which(group == g & is.finite(x))

    if (!length(idx)) next

    out[idx] <- x[idx] - stats::median(
      x[idx],
      na.rm = TRUE
    )
  }

  out
}

jackknife_spearman <- function(
  x,
  y,
  expected_sign = NA_real_
) {
  ok <- is.finite(x) & is.finite(y)
  ids <- which(ok)

  if (length(ids) < 5L) {
    return(list(
      values = numeric(),
      median = NA_real_,
      min = NA_real_,
      max = NA_real_,
      sign_stability_vs_full = NA_real_,
      expected_sign_stability = NA_real_
    ))
  }

  full <- safe_spearman(
    x[ids],
    y[ids]
  )

  vals <- vapply(
    ids,
    function(drop_id) {
      keep <- setdiff(ids, drop_id)

      safe_spearman(
        x[keep],
        y[keep]
      )
    },
    numeric(1)
  )

  finite_vals <- vals[is.finite(vals)]

  if (!length(finite_vals)) {
    return(list(
      values = vals,
      median = NA_real_,
      min = NA_real_,
      max = NA_real_,
      sign_stability_vs_full = NA_real_,
      expected_sign_stability = NA_real_
    ))
  }

  sign_stability <- if (
    is.finite(full) &&
    abs(full) > NUM_TOL
  ) {
    mean(
      sign(finite_vals) == sign(full)
    )
  } else {
    NA_real_
  }

  expected_stability <- if (
    is.finite(expected_sign) &&
    abs(expected_sign) == 1
  ) {
    mean(
      sign(finite_vals) == expected_sign
    )
  } else {
    NA_real_
  }

  list(
    values = vals,
    median = stats::median(
      finite_vals
    ),
    min = min(finite_vals),
    max = max(finite_vals),
    sign_stability_vs_full =
      sign_stability,
    expected_sign_stability =
      expected_stability
  )
}

permute_indices_unrestricted <- function(
  n,
  nperm
) {
  t(
    replicate(
      nperm,
      sample.int(n),
      simplify = "matrix"
    )
  )
}

permute_indices_within_group <- function(
  group,
  nperm
) {
  n <- length(group)
  base <- seq_len(n)

  out <- matrix(
    NA_integer_,
    nrow = nperm,
    ncol = n
  )

  split_idx <- split(
    base,
    group
  )

  for (b in seq_len(nperm)) {
    p <- base

    for (idx in split_idx) {
      p[idx] <- sample(
        idx,
        length(idx),
        replace = FALSE
      )
    }

    out[b, ] <- p
  }

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
    paste(
      "No hay corridas de selección completas en",
      base
    )
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

find_bidirectional_run <- function(
  selection_run,
  family
) {
  root <- file.path(
    selection_run,
    family,
    "BIDIRECTIONAL_FRAMEWORK"
  )

  assert(
    dir.exists(root),
    paste(
      "No existe BIDIRECTIONAL_FRAMEWORK para",
      family
    )
  )

  latest_file <- file.path(
    root,
    "LATEST_RUN.txt"
  )

  if (file.exists(latest_file)) {
    p <- trimws(
      readLines(
        latest_file,
        warn = FALSE
      )[1L]
    )

    if (
      nzchar(p) &&
      dir.exists(p) &&
      file.exists(file.path(
        p,
        "rds",
        "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"
      ))
    ) {
      return(
        canonical_path(
          p,
          TRUE
        )
      )
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
    paste(
      "No hay corrida bidireccional completa para",
      family
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

find_latest_subrun <- function(
  parent,
  folder,
  required_file
) {
  root <- file.path(
    parent,
    folder
  )

  assert(
    dir.exists(root),
    paste(
      "No existe:",
      root
    )
  )

  latest_file <- file.path(
    root,
    "LATEST_RUN.txt"
  )

  if (file.exists(latest_file)) {
    p <- trimws(
      readLines(
        latest_file,
        warn = FALSE
      )[1L]
    )

    if (
      nzchar(p) &&
      dir.exists(p) &&
      file.exists(
        file.path(
          p,
          required_file
        )
      )
    ) {
      return(
        canonical_path(
          p,
          TRUE
        )
      )
    }
  }

  dirs <- list.dirs(
    root,
    recursive = FALSE,
    full.names = TRUE
  )

  dirs <- dirs[file.exists(
    file.path(
      dirs,
      required_file
    )
  )]

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

find_crossomic_run <- function(
  bidir_run
) {
  find_latest_subrun(
    parent = bidir_run,
    folder = "CROSSOMIC_BRIDGE",
    required_file = file.path(
      "rds",
      "CROSSOMIC_BRIDGE_OBJECT.rds"
    )
  )
}

# ---------------------- ready_for_modeling TRAIN ----------------------------

get_train_matrix <- function(
  rfm,
  view_code
) {
  assert(
    view_code %in% names(rfm),
    paste(
      "ready_for_modeling no contiene vista",
      view_code
    )
  )

  assert(
    is.list(rfm[[view_code]]) &&
      "train" %in% names(
        rfm[[view_code]]
      ),
    paste(
      "ready_for_modeling no contiene",
      view_code,
      "$train"
    )
  )

  X <- as.matrix(
    rfm[[view_code]]$train
  )

  assert(
    !is.null(rownames(X)) &&
      !is.null(colnames(X)),
    paste(
      view_code,
      "$train necesita rownames=features y colnames=samples"
    )
  )

  X
}

extract_train_labels <- function(
  rfm,
  sample_ids
) {
  assert(
    "grupo_train" %in% names(rfm),
    "ready_for_modeling no contiene grupo_train"
  )

  g <- rfm$grupo_train

  if (is.data.frame(g)) {
    assert(
      all(c("id", "grupo") %in% names(g)),
      "grupo_train data.frame debe contener id/grupo"
    )

    idx <- match(
      sample_ids,
      as.character(g$id)
    )

    assert(
      all(!is.na(idx)),
      "Hay muestras TRAIN sin etiqueta en grupo_train"
    )

    y <- as.character(
      g$grupo[idx]
    )

  } else if (
    is.atomic(g) &&
    !is.null(names(g))
  ) {
    idx <- match(
      sample_ids,
      names(g)
    )

    assert(
      all(!is.na(idx)),
      "Hay muestras TRAIN sin etiqueta en grupo_train"
    )

    y <- as.character(
      g[idx]
    )

  } else {
    stop(
      "Estructura de grupo_train no reconocida",
      call. = FALSE
    )
  }

  assert(
    all(y %in% GROUPS),
    paste(
      "Etiquetas inesperadas:",
      paste(
        unique(
          y[!y %in% GROUPS]
        ),
        collapse = ", "
      )
    )
  )

  factor(
    y,
    levels = GROUPS
  )
}

validate_train_views <- function(
  rfm
) {
  mats <- lapply(
    c("tx", "pr", "me", "cl"),
    function(v) {
      get_train_matrix(
        rfm,
        v
      )
    }
  )

  names(mats) <- c(
    "tx",
    "pr",
    "me",
    "cl"
  )

  ids <- lapply(
    mats,
    colnames
  )

  ref <- ids$tx

  for (v in names(ids)) {
    assert(
      identical(
        ids[[v]],
        ref
      ),
      paste(
        "TRAIN IDs no coinciden entre tx y",
        v
      )
    )
  }

  mats
}

# ------------------------- scores de pathways -------------------------------

build_train_pathway_scores <- function(
  drivers,
  profiles,
  Xtx,
  sample_ids,
  labels
) {
  needed_drivers <- c(
    "library",
    "pathway",
    "factor",
    "oriented_side",
    "feature_original",
    "anchor_weight"
  )

  assert(
    all(
      needed_drivers %in%
        names(drivers)
    ),
    paste(
      "Faltan columnas drivers:",
      paste(
        setdiff(
          needed_drivers,
          names(drivers)
        ),
        collapse = ", "
      )
    )
  )

  needed_profiles <- c(
    "library",
    "pathway",
    "factor",
    "oriented_side",
    "assigned_class_nominal",
    "p_assigned_class"
  )

  assert(
    all(
      needed_profiles %in%
        names(profiles)
    ),
    paste(
      "Faltan columnas profiles:",
      paste(
        setdiff(
          needed_profiles,
          names(profiles)
        ),
        collapse = ", "
      )
    )
  )

  drivers$key <- pathway_key(
    drivers$library,
    drivers$pathway,
    drivers$factor,
    drivers$oriented_side
  )

  profiles$key <- pathway_key(
    profiles$library,
    profiles$pathway,
    profiles$factor,
    profiles$oriented_side
  )

  keys <- unique(
    profiles$key
  )

  score_rows <- list()
  spec_rows <- list()
  summary_rows <- list()

  rr <- 1L
  ss <- 1L
  uu <- 1L

  for (k in keys) {
    p <- profiles[
      profiles$key == k,
      ,
      drop = FALSE
    ]

    assert(
      nrow(p) == 1L,
      paste(
        "Perfil duplicado para",
        k
      )
    )

    d <- drivers[
      drivers$key == k,
      ,
      drop = FALSE
    ]

    if (!nrow(d)) next

    in_tx <- d$feature_original %in%
      rownames(Xtx)

    d$present_in_train_tx <- in_tx

    mu <- rep(
      NA_real_,
      nrow(d)
    )

    sdv <- rep(
      NA_real_,
      nrow(d)
    )

    for (j in seq_len(nrow(d))) {
      if (!in_tx[j]) next

      x <- as.numeric(
        Xtx[
          d$feature_original[j],
          sample_ids
        ]
      )

      mu[j] <- mean(
        x,
        na.rm = TRUE
      )

      sdv[j] <- stats::sd(
        x,
        na.rm = TRUE
      )
    }

    d$train_mean <- mu
    d$train_sd <- sdv

    d$usable_for_score <-
      d$present_in_train_tx &
      is.finite(d$train_mean) &
      is.finite(d$train_sd) &
      d$train_sd > NUM_TOL &
      is.finite(d$anchor_weight) &
      d$anchor_weight > 0

    usable <- which(
      d$usable_for_score
    )

    if (!length(usable)) next

    raw_w <- as.numeric(
      d$anchor_weight[usable]
    )

    norm_w <- raw_w / sum(raw_w)

    d$normalized_anchor_weight <- 0
    d$normalized_anchor_weight[usable] <-
      norm_w

    Z <- matrix(
      NA_real_,
      nrow = length(usable),
      ncol = length(sample_ids)
    )

    for (jj in seq_along(usable)) {
      j <- usable[jj]

      x <- as.numeric(
        Xtx[
          d$feature_original[j],
          sample_ids
        ]
      )

      Z[jj, ] <-
        (
          x -
            d$train_mean[j]
        ) /
        d$train_sd[j]
    }

    score <- as.numeric(
      crossprod(
        norm_w,
        Z
      )
    )

    score_rows[[rr]] <- data.frame(
      family =
        if ("family" %in% names(p)) {
          p$family[1L]
        } else {
          NA_character_
        },
      key = k,
      library = p$library[1L],
      pathway = p$pathway[1L],
      factor = p$factor[1L],
      oriented_side =
        p$oriented_side[1L],
      assigned_class_nominal =
        p$assigned_class_nominal[1L],
      p_assigned_class =
        p$p_assigned_class[1L],
      q_both_sides_within_library_factor =
        if (
          "q_both_sides_within_library_factor" %in%
            names(p)
        ) {
          p$q_both_sides_within_library_factor[1L]
        } else {
          NA_real_
        },
      sample_id = sample_ids,
      grupo = as.character(labels),
      pathway_score_observed_train =
        score,
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L

    d$family <-
      if ("family" %in% names(p)) {
        p$family[1L]
      } else {
        NA_character_
      }

    d$assigned_class_nominal <-
      p$assigned_class_nominal[1L]

    d$p_assigned_class <-
      p$p_assigned_class[1L]

    spec_rows[[ss]] <- d
    ss <- ss + 1L

    summary_rows[[uu]] <- data.frame(
      family =
        if ("family" %in% names(p)) {
          p$family[1L]
        } else {
          NA_character_
        },
      key = k,
      library = p$library[1L],
      pathway = p$pathway[1L],
      factor = p$factor[1L],
      oriented_side =
        p$oriented_side[1L],
      assigned_class_nominal =
        p$assigned_class_nominal[1L],
      p_assigned_class =
        p$p_assigned_class[1L],
      n_driver_genes_defined =
        nrow(d),
      n_driver_genes_present_train =
        sum(d$present_in_train_tx),
      n_driver_genes_used_score =
        sum(d$usable_for_score),
      stringsAsFactors = FALSE
    )

    uu <- uu + 1L
  }

  assert(
    length(score_rows) > 0L,
    "No se pudo construir ningún score observado de pathway"
  )

  list(
    scores = do.call(
      rbind,
      score_rows
    ),
    scoring_spec = do.call(
      rbind,
      spec_rows
    ),
    pathway_qc = do.call(
      rbind,
      summary_rows
    )
  )
}

# ---------------- pathway observado <-> clase -------------------------------

pathway_class_audit <- function(
  pathway_scores,
  nperm
) {
  keys <- unique(
    pathway_scores$key
  )

  rows <- list()
  rr <- 1L

  for (k in keys) {
    z <- pathway_scores[
      pathway_scores$key == k,
      ,
      drop = FALSE
    ]

    assigned <- unique(
      z$assigned_class_nominal
    )

    assert(
      length(assigned) == 1L &&
        assigned %in% GROUPS,
      paste(
        "Clase asignada inválida para",
        k
      )
    )

    score <- z$pathway_score_observed_train
    group <- factor(
      z$grupo,
      levels = GROUPS
    )

    medians <- vapply(
      GROUPS,
      function(g) {
        stats::median(
          score[group == g],
          na.rm = TRUE
        )
      },
      numeric(1)
    )

    other <- setdiff(
      GROUPS,
      assigned
    )

    stat_obs <-
      medians[assigned] -
      max(
        medians[other],
        na.rm = TRUE
      )

    # Estadístico preespecificado: la clase asignada debe tener el score
    # más alto que AMBAS otras clases.
    null_stats <- numeric(
      nperm
    )

    for (b in seq_len(nperm)) {
      gp <- sample(
        group,
        length(group),
        replace = FALSE
      )

      pm <- vapply(
        GROUPS,
        function(g) {
          stats::median(
            score[gp == g],
            na.rm = TRUE
          )
        },
        numeric(1)
      )

      null_stats[b] <-
        pm[assigned] -
        max(
          pm[other],
          na.rm = TRUE
        )
    }

    p_perm <- (
      1 +
        sum(
          null_stats >= stat_obs,
          na.rm = TRUE
        )
    ) / (
      nperm + 1
    )

    ord <- order(
      medians,
      decreasing = TRUE
    )

    observed_top <- GROUPS[
      ord[1L]
    ]

    rows[[rr]] <- data.frame(
      family = z$family[1L],
      key = k,
      library = z$library[1L],
      pathway = z$pathway[1L],
      factor = z$factor[1L],
      oriented_side =
        z$oriented_side[1L],
      assigned_class_nominal =
        assigned,
      p_assigned_class =
        z$p_assigned_class[1L],
      median_NP =
        unname(medians["NP"]),
      median_FD =
        unname(medians["FD"]),
      median_SO =
        unname(medians["SO"]),
      observed_top_class =
        observed_top,
      assigned_class_is_observed_top =
        observed_top == assigned,
      assigned_vs_best_other_median_difference =
        unname(stat_obs),
      p_perm_assigned_high_vs_both =
        p_perm,
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  out <- do.call(
    rbind,
    rows
  )

  out$q_perm_assigned_high_within_library <- NA_real_

  split_idx <- split(
    seq_len(nrow(out)),
    out$library
  )

  for (idx in split_idx) {
    out$q_perm_assigned_high_within_library[idx] <-
      stats::p.adjust(
        out$p_perm_assigned_high_vs_both[idx],
        method = PADJ_METHOD
      )
  }

  out$observed_class_support_q10 <-
    is.finite(
      out$q_perm_assigned_high_within_library
    ) &
    out$q_perm_assigned_high_within_library <=
      ALPHA_AUDIT &
    out$assigned_class_is_observed_top

  out
}

# ---------------- pathway observado <-> PR/ME/CL ----------------------------

collect_observed_link_data <- function(
  links,
  pathway_scores,
  train_mats
) {
  # TX se excluye deliberadamente de esta corroboración.
  links <- links[
    links$view_code %in%
      c("pr", "me", "cl"),
    ,
    drop = FALSE
  ]

  if (!nrow(links)) {
    return(list())
  }

  links$key <- pathway_key(
    links$library,
    links$pathway,
    links$factor,
    links$oriented_side
  )

  split_links <- split(
    seq_len(nrow(links)),
    links$key
  )

  out <- list()

  for (k in names(split_links)) {
    idx_links <- split_links[[k]]

    score_df <- pathway_scores[
      pathway_scores$key == k,
      ,
      drop = FALSE
    ]

    if (!nrow(score_df)) next

    sample_ids <- as.character(
      score_df$sample_id
    )

    score <- as.numeric(
      score_df$pathway_score_observed_train
    )

    group <- factor(
      score_df$grupo,
      levels = GROUPS
    )

    link_rows <- links[
      idx_links,
      ,
      drop = FALSE
    ]

    entries <- list()
    ee <- 1L

    for (i in seq_len(nrow(link_rows))) {
      lr <- link_rows[i, , drop = FALSE]

      vc <- as.character(
        lr$view_code
      )

      X <- train_mats[[vc]]

      feature <- as.character(
        lr$feature_original
      )

      if (!feature %in% rownames(X)) {
        next
      }

      assert(
        all(
          sample_ids %in%
            colnames(X)
        ),
        paste(
          "Faltan muestras TRAIN en",
          vc
        )
      )

      x <- as.numeric(
        X[
          feature,
          sample_ids
        ]
      )

      entries[[ee]] <- list(
        link_row = lr,
        feature_values = x
      )

      ee <- ee + 1L
    }

    if (length(entries)) {
      out[[k]] <- list(
        score_df = score_df,
        score = score,
        group = group,
        entries = entries
      )
    }
  }

  out
}

audit_pathway_feature_links <- function(
  observed_link_data,
  nperm
) {
  if (!length(observed_link_data)) {
    return(data.frame())
  }

  rank_standardize_vector <- function(v) {
    r <- rank(
      v,
      ties.method = "average"
    )

    r <- r - mean(r)
    denom <- sqrt(sum(r^2))

    if (
      !is.finite(denom) ||
      denom <= NUM_TOL
    ) {
      return(
        rep(
          NA_real_,
          length(v)
        )
      )
    }

    r / denom
  }

  rank_standardize_rows <- function(M) {
    R <- t(
      apply(
        M,
        1,
        rank,
        ties.method = "average"
      )
    )

    R <- sweep(
      R,
      1,
      rowMeans(R),
      "-"
    )

    denom <- sqrt(
      rowSums(
        R^2
      )
    )

    good <- is.finite(denom) &
      denom > NUM_TOL

    R[good, ] <- R[good, , drop = FALSE] /
      denom[good]

    R[!good, ] <- NA_real_

    R
  }

  matrix_spearman_complete <- function(
    X,
    y
  ) {
    Y <- rank_standardize_vector(
      y
    )

    XR <- rank_standardize_rows(
      X
    )

    as.numeric(
      XR %*% Y
    )
  }

  matrix_jackknife <- function(
    X,
    y,
    expected_sign
  ) {
    n <- ncol(X)
    m <- nrow(X)

    vals <- matrix(
      NA_real_,
      nrow = m,
      ncol = n
    )

    for (drop in seq_len(n)) {
      keep <- setdiff(
        seq_len(n),
        drop
      )

      vals[, drop] <-
        matrix_spearman_complete(
          X[, keep, drop = FALSE],
          y[keep]
        )
    }

    full <- matrix_spearman_complete(
      X,
      y
    )

    median_jk <- apply(
      vals,
      1,
      stats::median,
      na.rm = TRUE
    )

    min_jk <- apply(
      vals,
      1,
      min,
      na.rm = TRUE
    )

    max_jk <- apply(
      vals,
      1,
      max,
      na.rm = TRUE
    )

    stability_full <- rep(
      NA_real_,
      m
    )

    stability_expected <- rep(
      NA_real_,
      m
    )

    for (j in seq_len(m)) {
      v <- vals[j, ]
      v <- v[is.finite(v)]

      if (!length(v)) next

      if (
        is.finite(full[j]) &&
        abs(full[j]) > NUM_TOL
      ) {
        stability_full[j] <-
          mean(
            sign(v) ==
              sign(full[j])
          )
      }

      if (
        is.finite(
          expected_sign[j]
        ) &&
        abs(
          expected_sign[j]
        ) == 1
      ) {
        stability_expected[j] <-
          mean(
            sign(v) ==
              expected_sign[j]
          )
      }
    }

    list(
      full = full,
      median = median_jk,
      min = min_jk,
      max = max_jk,
      sign_stability_vs_full =
        stability_full,
      expected_sign_stability =
        stability_expected
    )
  }

  rows <- list()
  rr <- 1L

  for (k in names(observed_link_data)) {
    obj <- observed_link_data[[k]]

    score <- as.numeric(
      obj$score
    )

    group <- factor(
      obj$group,
      levels = GROUPS
    )

    n <- length(score)
    m <- length(
      obj$entries
    )

    assert(
      all(is.finite(score)),
      paste(
        "Score de pathway contiene NA/Inf:",
        k
      )
    )

    # Matriz features x samples.
    X <- do.call(
      rbind,
      lapply(
        obj$entries,
        function(e) {
          as.numeric(
            e$feature_values
          )
        }
      )
    )

    assert(
      ncol(X) == n,
      "Dimensión X/pathway incompatible"
    )

    complete_feature <- apply(
      X,
      1,
      function(v) {
        all(is.finite(v)) &&
          stats::sd(v) >
            NUM_TOL
      }
    )

    meta <- do.call(
      rbind,
      lapply(
        obj$entries,
        function(e) {
          e$link_row
        }
      )
    )

    expected_sign <- ifelse(
      is.finite(
        meta$pathway_feature_cosine
      ) &
        abs(
          meta$pathway_feature_cosine
        ) > NUM_TOL,
      sign(
        meta$pathway_feature_cosine
      ),
      NA_real_
    )

    # Inicializar resultados para todas las features.
    rho <- rep(
      NA_real_,
      m
    )

    p_two <- rep(
      NA_real_,
      m
    )

    p_dir <- rep(
      NA_real_,
      m
    )

    rho_within <- rep(
      NA_real_,
      m
    )

    p_within_two <- rep(
      NA_real_,
      m
    )

    p_within_dir <- rep(
      NA_real_,
      m
    )

    jk_median <- rep(
      NA_real_,
      m
    )

    jk_min <- rep(
      NA_real_,
      m
    )

    jk_max <- rep(
      NA_real_,
      m
    )

    jk_stability_full <- rep(
      NA_real_,
      m
    )

    jk_stability_expected <- rep(
      NA_real_,
      m
    )

    idx <- which(
      complete_feature
    )

    if (length(idx)) {
      Xc <- X[
        idx,
        ,
        drop = FALSE
      ]

      expected_c <-
        expected_sign[idx]

      # -----------------------
      # Spearman global.
      # -----------------------
      Yrank <- rank_standardize_vector(
        score
      )

      Xrank <- rank_standardize_rows(
        Xc
      )

      rho_c <- as.numeric(
        Xrank %*% Yrank
      )

      rho[idx] <- rho_c

      # Una sola familia de permutaciones por pathway,
      # reutilizada para todas sus features.
      perm_global <-
        permute_indices_unrestricted(
          n,
          nperm
        )

      Yperm <- t(
        vapply(
          seq_len(nperm),
          function(b) {
            Yrank[
              perm_global[b, ]
            ]
          },
          numeric(n)
        )
      )

      rho_perm <- Xrank %*%
        t(Yperm)

      p_two_c <- (
        1 +
          rowSums(
            sweep(
              abs(
                rho_perm
              ),
              1,
              abs(
                rho_c
              ),
              ">="
            ),
            na.rm = TRUE
          )
      ) / (
        nperm + 1
      )

      p_two[idx] <- p_two_c

      valid_expected <-
        is.finite(
          expected_c
        )

      if (any(valid_expected)) {
        aligned_perm <- sweep(
          rho_perm[
            valid_expected,
            ,
            drop = FALSE
          ],
          1,
          expected_c[
            valid_expected
          ],
          "*"
        )

        aligned_obs <-
          expected_c[
            valid_expected
          ] *
          rho_c[
            valid_expected
          ]

        p_dir_c <- (
          1 +
            rowSums(
              sweep(
                aligned_perm,
                1,
                aligned_obs,
                ">="
              ),
              na.rm = TRUE
            )
        ) / (
          nperm + 1
        )

        p_dir[
          idx[
            valid_expected
          ]
        ] <- p_dir_c
      }

      # -----------------------
      # Within-class centered.
      # -----------------------
      score_within <-
        group_median_center(
          score,
          group
        )

      X_within <- Xc

      for (g in GROUPS) {
        gi <- which(
          group == g
        )

        if (!length(gi)) next

        med <- apply(
          X_within[
            ,
            gi,
            drop = FALSE
          ],
          1,
          stats::median,
          na.rm = TRUE
        )

        X_within[
          ,
          gi
        ] <- sweep(
          X_within[
            ,
            gi,
            drop = FALSE
          ],
          1,
          med,
          "-"
        )
      }

      valid_within <-
        apply(
          X_within,
          1,
          function(v) {
            all(is.finite(v)) &&
              stats::sd(v) >
                NUM_TOL
          }
        ) &
        all(
          is.finite(
            score_within
          )
        ) &
        stats::sd(
          score_within
        ) > NUM_TOL

      if (any(valid_within)) {
        Xw <- X_within[
          valid_within,
          ,
          drop = FALSE
        ]

        expected_w <-
          expected_c[
            valid_within
          ]

        Ywrank <-
          rank_standardize_vector(
            score_within
          )

        Xwrank <-
          rank_standardize_rows(
            Xw
          )

        rho_w <- as.numeric(
          Xwrank %*%
            Ywrank
        )

        target_idx <-
          idx[
            valid_within
          ]

        rho_within[
          target_idx
        ] <- rho_w

        perm_within <-
          permute_indices_within_group(
            group,
            nperm
          )

        Ywperm <- t(
          vapply(
            seq_len(nperm),
            function(b) {
              Ywrank[
                perm_within[b, ]
              ]
            },
            numeric(n)
          )
        )

        rho_w_perm <-
          Xwrank %*%
          t(Ywperm)

        p_w_two <- (
          1 +
            rowSums(
              sweep(
                abs(
                  rho_w_perm
                ),
                1,
                abs(
                  rho_w
                ),
                ">="
              ),
              na.rm = TRUE
            )
        ) / (
          nperm + 1
        )

        p_within_two[
          target_idx
        ] <- p_w_two

        valid_expected_w <-
          is.finite(
            expected_w
          )

        if (
          any(
            valid_expected_w
          )
        ) {
          aligned_w_perm <-
            sweep(
              rho_w_perm[
                valid_expected_w,
                ,
                drop = FALSE
              ],
              1,
              expected_w[
                valid_expected_w
              ],
              "*"
            )

          aligned_w_obs <-
            expected_w[
              valid_expected_w
            ] *
            rho_w[
              valid_expected_w
            ]

          p_w_dir <- (
            1 +
              rowSums(
                sweep(
                  aligned_w_perm,
                  1,
                  aligned_w_obs,
                  ">="
                ),
                na.rm = TRUE
              )
          ) / (
            nperm + 1
          )

          p_within_dir[
            target_idx[
              valid_expected_w
            ]
          ] <- p_w_dir
        }
      }

      # -----------------------
      # Jackknife vectorizado.
      # -----------------------
      jk <- matrix_jackknife(
        Xc,
        score,
        expected_c
      )

      jk_median[idx] <-
        jk$median

      jk_min[idx] <-
        jk$min

      jk_max[idx] <-
        jk$max

      jk_stability_full[idx] <-
        jk$sign_stability_vs_full

      jk_stability_expected[idx] <-
        jk$expected_sign_stability
    }

    z <- data.frame(
      family =
        if (
          "family" %in%
            names(meta)
        ) {
          meta$family
        } else {
          NA_character_
        },
      key = k,
      library =
        meta$library,
      pathway =
        meta$pathway,
      factor =
        meta$factor,
      oriented_side =
        meta$oriented_side,
      assigned_class_nominal =
        meta$assigned_class_nominal,
      p_assigned_class =
        meta$p_assigned_class,
      view =
        meta$view,
      view_code =
        meta$view_code,
      feature_original =
        meta$feature_original,
      feature_model =
        meta$feature_model,
      latent_pathway_feature_cosine =
        meta$pathway_feature_cosine,
      expected_observed_sign =
        expected_sign,
      rho_spearman_train =
        rho,
      expected_aligned_rho_train =
        expected_sign * rho,
      direction_concordant_latent_train =
        ifelse(
          is.finite(rho) &
            abs(rho) >
              NUM_TOL &
            is.finite(
              expected_sign
            ),
          sign(rho) ==
            expected_sign,
          NA
        ),
      p_perm_directional_train =
        p_dir,
      p_perm_two_sided_train =
        p_two,
      rho_within_class_centered_train =
        rho_within,
      expected_aligned_rho_within_class =
        expected_sign *
        rho_within,
      direction_concordant_within_class =
        ifelse(
          is.finite(
            rho_within
          ) &
            abs(
              rho_within
            ) >
              NUM_TOL &
            is.finite(
              expected_sign
            ),
          sign(
            rho_within
          ) ==
            expected_sign,
          NA
        ),
      p_perm_directional_within_class =
        p_within_dir,
      p_perm_two_sided_within_class =
        p_within_two,
      jackknife_rho_median =
        jk_median,
      jackknife_rho_min =
        jk_min,
      jackknife_rho_max =
        jk_max,
      jackknife_sign_stability_vs_full =
        jk_stability_full,
      jackknife_expected_sign_stability =
        jk_stability_expected,
      n_complete =
        rowSums(
          is.finite(X)
        ),
      stringsAsFactors = FALSE
    )

    rows[[rr]] <- z
    rr <- rr + 1L
  }

  if (!length(rows)) {
    return(data.frame())
  }

  out <- do.call(
    rbind,
    rows
  )

  out$q_perm_directional_within_library_view <-
    NA_real_

  out$q_perm_two_sided_within_library_view <-
    NA_real_

  out$q_perm_directional_within_class_library_view <-
    NA_real_

  key_bh <- interaction(
    out$library,
    out$view_code,
    drop = TRUE,
    lex.order = TRUE
  )

  split_idx <- split(
    seq_len(nrow(out)),
    key_bh
  )

  for (idx in split_idx) {
    out$q_perm_directional_within_library_view[idx] <-
      stats::p.adjust(
        out$p_perm_directional_train[idx],
        method = PADJ_METHOD
      )

    out$q_perm_two_sided_within_library_view[idx] <-
      stats::p.adjust(
        out$p_perm_two_sided_train[idx],
        method = PADJ_METHOD
      )

    out$q_perm_directional_within_class_library_view[idx] <-
      stats::p.adjust(
        out$p_perm_directional_within_class[idx],
        method = PADJ_METHOD
      )
  }

  out$observed_train_directional_support_q10 <-
    is.finite(
      out$q_perm_directional_within_library_view
    ) &
    out$q_perm_directional_within_library_view <=
      ALPHA_AUDIT &
    out$direction_concordant_latent_train

  out
}

# ---------------------- resumen por biomarcador -----------------------------

summarize_observed_links_by_candidate <- function(
  assoc
) {
  if (!nrow(assoc)) {
    return(data.frame())
  }

  key <- paste(
    assoc$factor,
    assoc$oriented_side,
    assoc$feature_model,
    sep = "\r"
  )

  split_idx <- split(
    seq_len(nrow(assoc)),
    key
  )

  rows <- list()
  rr <- 1L

  for (idx in split_idx) {
    z <- assoc[
      idx,
      ,
      drop = FALSE
    ]

    valid_dir <-
      !is.na(
        z$direction_concordant_latent_train
      )

    rows[[rr]] <- data.frame(
      factor = z$factor[1L],
      oriented_side =
        z$oriented_side[1L],
      feature_model =
        z$feature_model[1L],
      observed_train_n_pathway_links_tested =
        nrow(z),
      observed_train_n_direction_concordant =
        sum(
          z$direction_concordant_latent_train,
          na.rm = TRUE
        ),
      observed_train_prop_direction_concordant =
        if (any(valid_dir)) {
          mean(
            z$direction_concordant_latent_train[
              valid_dir
            ]
          )
        } else {
          NA_real_
        },
      observed_train_median_expected_aligned_rho =
        stats::median(
          z$expected_aligned_rho_train,
          na.rm = TRUE
        ),
      observed_train_max_expected_aligned_rho =
        max(
          z$expected_aligned_rho_train,
          na.rm = TRUE
        ),
      observed_train_best_q_directional =
        min(
          z$q_perm_directional_within_library_view,
          na.rm = TRUE
        ),
      observed_train_n_q10_directional =
        sum(
          z$observed_train_directional_support_q10,
          na.rm = TRUE
        ),
      observed_train_median_jackknife_expected_sign_stability =
        stats::median(
          z$jackknife_expected_sign_stability,
          na.rm = TRUE
        ),
      observed_train_median_within_class_aligned_rho =
        stats::median(
          z$expected_aligned_rho_within_class,
          na.rm = TRUE
        ),
      stringsAsFactors = FALSE
    )

    # Evitar Inf cuando todo es NA.
    numeric_names <- c(
      "observed_train_median_expected_aligned_rho",
      "observed_train_max_expected_aligned_rho",
      "observed_train_best_q_directional",
      "observed_train_median_jackknife_expected_sign_stability",
      "observed_train_median_within_class_aligned_rho"
    )

    for (nm in numeric_names) {
      if (
        !is.finite(
          rows[[rr]][[nm]]
        )
      ) {
        rows[[rr]][[nm]] <- NA_real_
      }
    }

    rr <- rr + 1L
  }

  do.call(
    rbind,
    rows
  )
}

attach_observed_train_to_biomarker_matrix <- function(
  biomarker_matrix,
  observed_summary
) {
  if (!nrow(observed_summary)) {
    return(
      biomarker_matrix
    )
  }

  key_b <- paste(
    biomarker_matrix$factor,
    biomarker_matrix$oriented_side,
    biomarker_matrix$feature_model,
    sep = "\r"
  )

  key_o <- paste(
    observed_summary$factor,
    observed_summary$oriented_side,
    observed_summary$feature_model,
    sep = "\r"
  )

  idx <- match(
    key_b,
    key_o
  )

  add_cols <- setdiff(
    names(observed_summary),
    c(
      "factor",
      "oriented_side",
      "feature_model"
    )
  )

  for (cc in add_cols) {
    biomarker_matrix[[cc]] <-
      observed_summary[[cc]][idx]
  }

  biomarker_matrix
}

# ----------------------------- una familia ----------------------------------

run_family <- function(
  selection_run,
  family,
  stamp
) {
  message("\n============================================================")
  message("OBSERVED TRAIN CORROBORATION | ", family)
  message("============================================================")

  bidir_run <- find_bidirectional_run(
    selection_run,
    family
  )

  crossomic_run <- find_crossomic_run(
    bidir_run
  )

  crossomic_path <- file.path(
    crossomic_run,
    "rds",
    "CROSSOMIC_BRIDGE_OBJECT.rds"
  )

  crossomic <- readRDS(
    crossomic_path
  )

  required_crossomic <- c(
    "pathway_profiles",
    "pathway_driver_genes",
    "pathway_feature_links",
    "biomarker_evidence_matrix"
  )

  assert(
    all(
      required_crossomic %in%
        names(crossomic)
    ),
    paste(
      "CROSSOMIC_BRIDGE_OBJECT incompleto:",
      paste(
        setdiff(
          required_crossomic,
          names(crossomic)
        ),
        collapse = ", "
      )
    )
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
    "rfm_path" %in%
      names(bidir_object),
    "BIDIRECTIONAL_FRAMEWORK_OBJECT no contiene rfm_path"
  )

  rfm_path <- canonical_path(
    bidir_object$rfm_path,
    TRUE
  )

  rfm <- readRDS(
    rfm_path
  )

  train_mats <- validate_train_views(
    rfm
  )

  sample_ids <- colnames(
    train_mats$tx
  )

  labels <- extract_train_labels(
    rfm,
    sample_ids
  )

  assert(
    identical(
      names(
        table(labels)
      ),
      GROUPS
    ),
    "No están representadas NP/FD/SO en TRAIN"
  )

  message(
    "TRAIN: ",
    paste(
      paste0(
        names(table(labels)),
        "=",
        as.integer(table(labels))
      ),
      collapse = " | "
    )
  )

  # -------------------------------------------------------------------------
  # 1) Score observado TRAIN de cada pathway predefinido.
  # -------------------------------------------------------------------------
  pathway_score_obj <-
    build_train_pathway_scores(
      drivers =
        crossomic$pathway_driver_genes,
      profiles =
        crossomic$pathway_profiles,
      Xtx =
        train_mats$tx,
      sample_ids =
        sample_ids,
      labels =
        labels
    )

  pathway_scores <-
    pathway_score_obj$scores

  pathway_scoring_spec <-
    pathway_score_obj$scoring_spec

  pathway_qc <-
    pathway_score_obj$pathway_qc

  # -------------------------------------------------------------------------
  # 2) ¿El pathway observado conserva la clase esperada?
  # -------------------------------------------------------------------------
  class_audit <-
    pathway_class_audit(
      pathway_scores =
        pathway_scores,
      nperm =
        NPERM
    )

  # -------------------------------------------------------------------------
  # 3) ¿Las asociaciones pathway -> PR/ME/CL aparecen en observado TRAIN?
  # -------------------------------------------------------------------------
  observed_link_data <-
    collect_observed_link_data(
      links =
        crossomic$pathway_feature_links,
      pathway_scores =
        pathway_scores,
      train_mats =
        train_mats
    )

  assoc <-
    audit_pathway_feature_links(
      observed_link_data =
        observed_link_data,
      nperm =
        NPERM
    )

  observed_summary <-
    summarize_observed_links_by_candidate(
      assoc
    )

  biomarker_matrix_train <-
    attach_observed_train_to_biomarker_matrix(
      biomarker_matrix =
        crossomic$biomarker_evidence_matrix,
      observed_summary =
        observed_summary
    )

  # Añadir soporte observado del pathway de la clase a cada biomarcador,
  # según su factor/cara.
  module_class_summary <- if (
    nrow(class_audit)
  ) {
    keys <- unique(
      class_audit[, c(
        "factor",
        "oriented_side"
      )]
    )

    rows <- vector(
      "list",
      nrow(keys)
    )

    for (i in seq_len(nrow(keys))) {
      k <- keys[i, , drop = FALSE]

      z <- class_audit[
        class_audit$factor ==
          k$factor &
          class_audit$oriented_side ==
            k$oriented_side,
        ,
        drop = FALSE
      ]

      rows[[i]] <- data.frame(
        factor = k$factor,
        oriented_side =
          k$oriented_side,
        n_pathways_class_audited =
          nrow(z),
        prop_pathways_assigned_class_observed_top =
          mean(
            z$assigned_class_is_observed_top,
            na.rm = TRUE
          ),
        n_pathways_class_support_q10 =
          sum(
            z$observed_class_support_q10,
            na.rm = TRUE
          ),
        stringsAsFactors = FALSE
      )
    }

    do.call(
      rbind,
      rows
    )
  } else {
    data.frame()
  }

  if (
    nrow(module_class_summary)
  ) {
    key_b <- paste(
      biomarker_matrix_train$factor,
      biomarker_matrix_train$oriented_side,
      sep = "\r"
    )

    key_m <- paste(
      module_class_summary$factor,
      module_class_summary$oriented_side,
      sep = "\r"
    )

    idx <- match(
      key_b,
      key_m
    )

    for (
      cc in setdiff(
        names(module_class_summary),
        c(
          "factor",
          "oriented_side"
        )
      )
    ) {
      biomarker_matrix_train[[cc]] <-
        module_class_summary[[cc]][idx]
    }
  }

  # -------------------------------------------------------------------------
  # Salidas.
  # -------------------------------------------------------------------------
  out <- file.path(
    bidir_run,
    "OBSERVED_TRAIN_CORROBORATION",
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
    bidirectional_run =
      bidir_run,
    crossomic_run =
      crossomic_run,
    ready_for_modeling =
      rfm_path,
    n_train =
      length(sample_ids),
    n_NP =
      sum(labels == "NP"),
    n_FD =
      sum(labels == "FD"),
    n_SO =
      sum(labels == "SO"),
    n_primary_pathways_scored =
      length(
        unique(
          pathway_scores$key
        )
      ),
    n_pathway_feature_links_observed_tested =
      nrow(assoc),
    nperm =
      NPERM,
    seed =
      SEED,
    padj_method =
      PADJ_METHOD,
    alpha_audit =
      ALPHA_AUDIT,
    TEST_used =
      FALSE,
    new_model_fitted =
      FALSE,
    modules_reselected =
      FALSE,
    pathways_reselected =
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
    pathway_qc,
    file.path(
      out,
      "tables",
      "01_PATHWAY_SCORE_QC.csv"
    )
  )

  write_csv(
    pathway_scoring_spec,
    file.path(
      out,
      "tables",
      "01B_PATHWAY_SCORING_SPEC_TRAIN_FROZEN.csv"
    )
  )

  write_csv(
    pathway_scores,
    file.path(
      out,
      "tables",
      "02_PATHWAY_SCORES_OBSERVED_TRAIN.csv"
    )
  )

  write_csv(
    class_audit,
    file.path(
      out,
      "tables",
      "03_PATHWAY_CLASS_DIRECTION_OBSERVED_TRAIN.csv"
    )
  )

  write_csv(
    assoc,
    file.path(
      out,
      "tables",
      "04_PATHWAY_FEATURE_ASSOCIATIONS_OBSERVED_TRAIN.csv"
    )
  )

  write_csv(
    observed_summary,
    file.path(
      out,
      "tables",
      "05_CANDIDATE_OBSERVED_TRAIN_SUPPORT_SUMMARY.csv"
    )
  )

  write_csv(
    module_class_summary,
    file.path(
      out,
      "tables",
      "05B_MODULE_PATHWAY_CLASS_SUPPORT_SUMMARY.csv"
    )
  )

  write_csv(
    biomarker_matrix_train,
    file.path(
      out,
      "tables",
      "06_BIOMARKER_EVIDENCE_MATRIX_WITH_OBSERVED_TRAIN.csv"
    )
  )

  # Congelar inputs exactos que debe consumir 6B.
  test_projection_contract <- list(
    family = family,
    created_from_train_only = TRUE,
    TEST_used_in_creation = FALSE,
    sample_ids_train =
      sample_ids,
    groups_train =
      as.character(labels),
    pathway_scoring_spec =
      pathway_scoring_spec,
    pathway_profiles =
      crossomic$pathway_profiles,
    pathway_feature_links =
      crossomic$pathway_feature_links,
    pathway_class_audit_train =
      class_audit,
    pathway_feature_associations_train =
      assoc,
    biomarker_evidence_matrix_train =
      biomarker_matrix_train,
    parameters = list(
      nperm = NPERM,
      seed = SEED,
      padj_method =
        PADJ_METHOD,
      alpha_audit =
        ALPHA_AUDIT
    )
  )

  saveRDS(
    list(
      family = family,
      input_audit =
        input_audit,
      pathway_qc =
        pathway_qc,
      pathway_scoring_spec =
        pathway_scoring_spec,
      pathway_scores_train =
        pathway_scores,
      pathway_class_audit_train =
        class_audit,
      pathway_feature_associations_train =
        assoc,
      candidate_observed_train_summary =
        observed_summary,
      module_class_summary =
        module_class_summary,
      biomarker_evidence_matrix_train =
        biomarker_matrix_train,
      TEST_PROJECTION_CONTRACT =
        test_projection_contract
    ),
    file.path(
      out,
      "rds",
      "OBSERVED_TRAIN_CORROBORATION_OBJECT.rds"
    ),
    version = 3
  )

  saveRDS(
    test_projection_contract,
    file.path(
      out,
      "rds",
      "TEST_PROJECTION_CONTRACT_FROZEN.rds"
    ),
    version = 3
  )

  readme <- c(
    paste(
      "METATEST |",
      family,
      "| PASO 6A OBSERVED TRAIN CORROBORATION"
    ),
    "",
    "OBJETIVO:",
    "- Corroborar en datos observados TRAIN relaciones ya descubiertas.",
    "- No construir nuevas relaciones usando el observado.",
    "",
    "PATHWAY SCORE:",
    "- Genes y anchor_weight provienen del CROSSOMIC_BRIDGE.",
    "- Cada gen se estandariza con media y SD calculadas SOLO en TRAIN.",
    "- Los anchor_weight positivos se normalizan para sumar 1.",
    "- No hay selección por clase durante la construcción del score.",
    "",
    "PATHWAY -> CLASE:",
    "- Se espera que la clase asignada por el bridge tenga el score más alto.",
    "- Estadístico: mediana(clase asignada) - mejor mediana de las otras dos.",
    "- Permutación de etiquetas dentro de TRAIN preservando tamaños de grupo.",
    "- BH por librería se reporta como auditoría, no como filtro.",
    "",
    "PATHWAY -> PR/ME/CL:",
    "- Spearman global evalúa la relación fenotípica observada.",
    "- expected_sign procede del cosine latente del PASO 5.",
    "- Se reporta concordancia del signo observado con el signo latente.",
    "- p_perm_directional_train usa la dirección preespecificada por el bridge.",
    "- También se reporta p de permutación bilateral.",
    "- Jackknife leave-one-out evalúa estabilidad del signo.",
    "",
    "WITHIN-CLASS:",
    "- Pathway y feature se centran por la mediana de cada clase.",
    "- Después se calcula Spearman.",
    "- Las permutaciones se restringen dentro de NP/FD/SO.",
    "- Es evidencia secundaria: pregunta por covariación adicional no explicada",
    "  únicamente por el desplazamiento entre clases.",
    "",
    "BIOMARCADORES:",
    "- 06_BIOMARKER_EVIDENCE_MATRIX_WITH_OBSERVED_TRAIN.csv añade esta nueva",
    "  capa de evidencia a los candidatos del PASO 5.",
    "- No se genera un score compuesto ni una selección binaria.",
    "- Un resultado no significativo no elimina automáticamente un candidato.",
    "",
    "TEST:",
    "- TEST NO se usa en este script.",
    "- TEST_PROJECTION_CONTRACT_FROZEN.rds congela genes, medias/SD TRAIN,",
    "  pesos y relaciones a auditar en PASO 6B.",
    "",
    "SIGUIENTE PASO:",
    "- PASO 6B: aplicar exactamente este contrato a TEST sin reselección.",
    "- Con n=10 TEST se priorizará concordancia direccional/transportabilidad",
    "  sobre inferencia de significación."
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
    "OBSERVED_TRAIN_CORROBORATION"
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
    n_pathways_scored =
      length(
        unique(
          pathway_scores$key
        )
      ),
    n_pathways_class_direction_concordant =
      sum(
        class_audit$assigned_class_is_observed_top,
        na.rm = TRUE
      ),
    n_crossomic_links_tested =
      nrow(assoc),
    n_crossomic_links_direction_concordant =
      if (nrow(assoc)) {
        sum(
          assoc$direction_concordant_latent_train,
          na.rm = TRUE
        )
      } else {
        0L
      },
    output = out,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------- MAIN ------------------------------------

main <- function() {
  assert(
    NPERM >= 1000L,
    "Usar al menos 1000 permutaciones"
  )

  assert(
    is.finite(ALPHA_AUDIT) &&
      ALPHA_AUDIT > 0 &&
      ALPHA_AUDIT < 1,
    "OBSERVED_TRAIN_ALPHA inválido"
  )

  assert(
    PADJ_METHOD %in%
      stats::p.adjust.methods,
    "OBSERVED_TRAIN_PADJ inválido"
  )

  set.seed(
    SEED
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
        selection_run =
          selection_run,
        family =
          family,
        stamp =
          stamp
      ),
      error = function(e) {
        data.frame(
          family = family,
          n_pathways_scored =
            NA_integer_,
          n_pathways_class_direction_concordant =
            NA_integer_,
          n_crossomic_links_tested =
            NA_integer_,
          n_crossomic_links_direction_concordant =
            NA_integer_,
          output =
            NA_character_,
          error =
            conditionMessage(e),
          stringsAsFactors =
            FALSE
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
      "OBSERVED_TRAIN_CORROBORATION_RUN_",
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

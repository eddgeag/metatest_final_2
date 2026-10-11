#!/usr/bin/env Rscript
# =============================================================================
# METATEST — PASO 6B
# Auditoría TEST congelada + auditoría pooled descriptiva
# Version 0.1.0
#
# OBJETIVO
#   Aplicar a TEST, sin ninguna reselección, el contrato congelado por PASO 6A:
#
#       TEST_PROJECTION_CONTRACT_FROZEN.rds
#
#   y evaluar transportabilidad de:
#
#       1) pathway observado <-> clase esperada
#       2) pathway observado <-> PR / ME / CL
#       3) biomarcador candidato <-> clase esperada
#
#   Después construir una auditoría descriptiva TRAIN + TEST (n=26), siempre
#   identificando dataset y SIN tratar el pooled como validación independiente.
#
# REGLAS
#   - TEST no selecciona pathways, genes, pesos, módulos ni candidatos.
#   - Scores TEST usan EXCLUSIVAMENTE mean/SD TRAIN y pesos congelados en 6A.
#   - Las parejas pathway-feature son exactamente las auditadas en TRAIN.
#   - Con TEST pequeño, la evidencia principal es DIRECCIONAL / TRANSPORTABILIDAD.
#   - No se usa significación TEST para eliminar biomarcadores.
#   - POOLED es DESCRIPTIVO; no genera evidencia confirmatoria.
#
# SALIDA FINAL
#   BIOMARKER_EVIDENCE_MATRIX_FINAL.csv integra:
#     latente + pathway + bidireccional + observado TRAIN + auditoría TEST.
#
#   Después de este paso no se añade otro modelo. La siguiente fase es síntesis
#   biológica y visualización.
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

pathway_key <- function(library, pathway, factor, oriented_side) {
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

jackknife_spearman <- function(x, y, expected_sign = NA_real_) {
  ok <- is.finite(x) & is.finite(y)
  ids <- which(ok)

  if (length(ids) < 5L) {
    return(list(
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

  vals <- vals[is.finite(vals)]

  if (!length(vals)) {
    return(list(
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
    mean(sign(vals) == sign(full))
  } else {
    NA_real_
  }

  expected_stability <- if (
    is.finite(expected_sign) &&
    abs(expected_sign) == 1
  ) {
    mean(sign(vals) == expected_sign)
  } else {
    NA_real_
  }

  list(
    median = stats::median(vals),
    min = min(vals),
    max = max(vals),
    sign_stability_vs_full = sign_stability,
    expected_sign_stability = expected_stability
  )
}

median_center_by_dataset <- function(x, dataset) {
  out <- rep(NA_real_, length(x))

  for (d in unique(dataset)) {
    idx <- which(dataset == d & is.finite(x))

    if (!length(idx)) next

    out[idx] <- x[idx] - stats::median(
      x[idx],
      na.rm = TRUE
    )
  }

  out
}

transport_status <- function(train_value, test_value) {
  if (
    !is.finite(train_value) ||
    !is.finite(test_value) ||
    abs(train_value) <= NUM_TOL ||
    abs(test_value) <= NUM_TOL
  ) {
    return("UNRESOLVED")
  }

  if (sign(train_value) == sign(test_value)) {
    "DIRECTION_PRESERVED"
  } else {
    "DIRECTION_REVERSED"
  }
}

# --------------------------- localizar corridas -----------------------------

find_latest_selection_run <- function(upstream) {
  base <- file.path(
    upstream,
    "downstream",
    "seleccion_mejor_escenario_train"
  )

  assert(dir.exists(base), paste("No existe:", base))

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

find_bidirectional_run <- function(selection_run, family) {
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

find_latest_subrun <- function(parent, folder, required_file) {
  root <- file.path(
    parent,
    folder
  )

  assert(
    dir.exists(root),
    paste("No existe:", root)
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
      return(canonical_path(p, TRUE))
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

find_observed_train_run <- function(bidir_run) {
  find_latest_subrun(
    parent = bidir_run,
    folder = "OBSERVED_TRAIN_CORROBORATION",
    required_file = file.path(
      "rds",
      "TEST_PROJECTION_CONTRACT_FROZEN.rds"
    )
  )
}

# ------------------------ ready_for_modeling --------------------------------

get_dataset_matrix <- function(rfm, view_code, dataset) {
  assert(
    view_code %in% names(rfm),
    paste(
      "ready_for_modeling no contiene vista",
      view_code
    )
  )

  assert(
    is.list(rfm[[view_code]]) &&
      dataset %in% names(rfm[[view_code]]),
    paste(
      "ready_for_modeling no contiene",
      view_code,
      "$",
      dataset
    )
  )

  X <- as.matrix(
    rfm[[view_code]][[dataset]]
  )

  assert(
    !is.null(rownames(X)) &&
      !is.null(colnames(X)),
    paste(
      view_code,
      "$",
      dataset,
      " necesita rownames=features y colnames=samples"
    )
  )

  assert(
    all(is.finite(X)),
    paste(
      "NA/Inf en",
      view_code,
      dataset
    )
  )

  X
}

get_split_table <- function(rfm, dataset) {
  nm <- paste0(
    "grupo_",
    dataset
  )

  assert(
    nm %in% names(rfm),
    paste(
      "ready_for_modeling no contiene",
      nm
    )
  )

  x <- rfm[[nm]]

  assert(
    is.data.frame(x) &&
      all(c("id", "grupo") %in% names(x)),
    paste(
      nm,
      "debe ser data.frame con id/grupo"
    )
  )

  out <- x[, c("id", "grupo"), drop = FALSE]
  out[] <- lapply(out, as.character)

  assert(
    !anyDuplicated(out$id),
    paste("IDs duplicados en", nm)
  )

  assert(
    all(out$grupo %in% GROUPS),
    paste("Grupos inesperados en", nm)
  )

  out
}

load_all_views <- function(rfm, dataset, split_table) {
  views <- c(
    tx = "tx",
    pr = "pr",
    me = "me",
    cl = "cl"
  )

  mats <- lapply(
    views,
    function(v) {
      X <- get_dataset_matrix(
        rfm,
        v,
        dataset
      )

      assert(
        all(
          split_table$id %in%
            colnames(X)
        ),
        paste(
          "Faltan IDs",
          dataset,
          "en vista",
          v
        )
      )

      X[
        ,
        split_table$id,
        drop = FALSE
      ]
    }
  )

  names(mats) <- names(views)

  ref <- split_table$id

  for (v in names(mats)) {
    assert(
      identical(
        colnames(mats[[v]]),
        ref
      ),
      paste(
        "Orden de IDs incorrecto en",
        dataset,
        v
      )
    )
  }

  mats
}

validate_train_contract <- function(contract, train_split) {
  assert(
    isTRUE(
      contract$created_from_train_only
    ),
    "El contrato no declara created_from_train_only=TRUE"
  )

  assert(
    identical(
      contract$TEST_used_in_creation,
      FALSE
    ),
    "El contrato indica uso de TEST durante su creación"
  )

  assert(
    setequal(
      as.character(
        contract$sample_ids_train
      ),
      train_split$id
    ),
    "TRAIN del contrato != TRAIN de ready_for_modeling"
  )

  idx <- match(
    contract$sample_ids_train,
    train_split$id
  )

  assert(
    identical(
      as.character(
        contract$groups_train
      ),
      train_split$grupo[idx]
    ),
    "Etiquetas TRAIN del contrato != ready_for_modeling"
  )

  invisible(TRUE)
}

# ------------------------ scores congelados ---------------------------------

project_pathway_scores <- function(
  scoring_spec,
  profiles,
  Xtx,
  sample_ids,
  labels,
  dataset_label
) {
  needed <- c(
    "library",
    "pathway",
    "factor",
    "oriented_side",
    "feature_original",
    "train_mean",
    "train_sd",
    "usable_for_score",
    "normalized_anchor_weight"
  )

  assert(
    all(needed %in% names(scoring_spec)),
    paste(
      "Scoring spec incompleto:",
      paste(
        setdiff(
          needed,
          names(scoring_spec)
        ),
        collapse = ", "
      )
    )
  )

  spec <- scoring_spec

  if (!"key" %in% names(spec)) {
    spec$key <- pathway_key(
      spec$library,
      spec$pathway,
      spec$factor,
      spec$oriented_side
    )
  }

  prof <- profiles

  if (!"key" %in% names(prof)) {
    prof$key <- pathway_key(
      prof$library,
      prof$pathway,
      prof$factor,
      prof$oriented_side
    )
  }

  keys <- unique(
    prof$key
  )

  rows <- list()
  qc <- list()
  rr <- 1L
  qq <- 1L

  for (k in keys) {
    p <- prof[
      prof$key == k,
      ,
      drop = FALSE
    ]

    assert(
      nrow(p) == 1L,
      paste(
        "Perfil pathway duplicado:",
        k
      )
    )

    d <- spec[
      spec$key == k &
        spec$usable_for_score,
      ,
      drop = FALSE
    ]

    assert(
      nrow(d) > 0L,
      paste(
        "Pathway sin genes congelados utilizables:",
        k
      )
    )

    assert(
      all(
        d$feature_original %in%
          rownames(Xtx)
      ),
      paste(
        "TEST/TRAIN no contiene todos los genes congelados de",
        k
      )
    )

    assert(
      all(
        is.finite(
          d$train_mean
        )
      ) &&
        all(
          is.finite(
            d$train_sd
          )
        ) &&
        all(
          d$train_sd >
            NUM_TOL
        ),
      paste(
        "mean/SD TRAIN inválidos en",
        k
      )
    )

    w <- as.numeric(
      d$normalized_anchor_weight
    )

    assert(
      all(is.finite(w)) &&
        all(w >= 0) &&
        abs(sum(w) - 1) < 1e-8,
      paste(
        "Pesos congelados no suman 1 en",
        k
      )
    )

    X <- Xtx[
      d$feature_original,
      sample_ids,
      drop = FALSE
    ]

    Z <- sweep(
      X,
      1,
      as.numeric(
        d$train_mean
      ),
      "-"
    )

    Z <- sweep(
      Z,
      1,
      as.numeric(
        d$train_sd
      ),
      "/"
    )

    assert(
      all(is.finite(Z)),
      paste(
        "Score estandarizado contiene NA/Inf:",
        k
      )
    )

    score <- as.numeric(
      crossprod(
        w,
        Z
      )
    )

    rows[[rr]] <- data.frame(
      family =
        if ("family" %in% names(p)) {
          p$family[1L]
        } else {
          NA_character_
        },
      dataset = dataset_label,
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
      sample_id = sample_ids,
      grupo = as.character(labels),
      pathway_score =
        score,
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L

    qc[[qq]] <- data.frame(
      family =
        if ("family" %in% names(p)) {
          p$family[1L]
        } else {
          NA_character_
        },
      dataset = dataset_label,
      key = k,
      library = p$library[1L],
      pathway = p$pathway[1L],
      factor = p$factor[1L],
      oriented_side =
        p$oriented_side[1L],
      n_genes_frozen =
        nrow(d),
      weight_sum =
        sum(w),
      train_scaling_reused =
        TRUE,
      stringsAsFactors = FALSE
    )

    qq <- qq + 1L
  }

  list(
    scores = do.call(
      rbind,
      rows
    ),
    qc = do.call(
      rbind,
      qc
    )
  )
}

# ------------------------- pathway <-> clase --------------------------------

audit_pathway_class_direction <- function(
  scores,
  train_class_audit = NULL
) {
  keys <- unique(
    scores$key
  )

  rows <- list()
  rr <- 1L

  for (k in keys) {
    z <- scores[
      scores$key == k,
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
        "Clase asignada inválida:",
        k
      )
    )

    medians <- vapply(
      GROUPS,
      function(g) {
        stats::median(
          z$pathway_score[
            z$grupo == g
          ],
          na.rm = TRUE
        )
      },
      numeric(1)
    )

    other <- setdiff(
      GROUPS,
      assigned
    )

    diff_expected <- unname(
      medians[assigned] -
        max(
          medians[other],
          na.rm = TRUE
        )
    )

    observed_top <- GROUPS[
      order(
        medians,
        decreasing = TRUE
      )[1L]
    ]

    train_row <- NULL

    if (
      !is.null(train_class_audit) &&
      nrow(train_class_audit)
    ) {
      train_row <- train_class_audit[
        train_class_audit$key == k,
        ,
        drop = FALSE
      ]
    }

    train_top_ok <- if (
      !is.null(train_row) &&
      nrow(train_row) == 1L &&
      "assigned_class_is_observed_top" %in%
        names(train_row)
    ) {
      as.logical(
        train_row$assigned_class_is_observed_top
      )
    } else {
      NA
    }

    train_diff <- if (
      !is.null(train_row) &&
      nrow(train_row) == 1L &&
      "assigned_vs_best_other_median_difference" %in%
        names(train_row)
    ) {
      as.numeric(
        train_row$assigned_vs_best_other_median_difference
      )
    } else {
      NA_real_
    }

    rows[[rr]] <- data.frame(
      family = z$family[1L],
      dataset = unique(z$dataset)[1L],
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
        diff_expected,
      train_assigned_class_is_observed_top =
        train_top_ok,
      train_assigned_vs_best_other_difference =
        train_diff,
      train_test_assigned_class_top_both =
        if (
          !is.na(train_top_ok)
        ) {
          train_top_ok &&
            observed_top == assigned
        } else {
          NA
        },
      class_direction_transport_status =
        if (
          is.na(train_top_ok)
        ) {
          "TRAIN_NOT_AVAILABLE"
        } else if (
          train_top_ok &&
          observed_top == assigned
        ) {
          "PRESERVED"
        } else if (
          train_top_ok &&
          observed_top != assigned
        ) {
          "TEST_DISCORDANT"
        } else if (
          !train_top_ok &&
          observed_top == assigned
        ) {
          "TEST_SUPPORTS"
        } else {
          "NEITHER"
        },
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  do.call(
    rbind,
    rows
  )
}

# ---------------- pathway <-> feature TEST ----------------------------------

audit_pathway_feature_test <- function(
  train_assoc,
  test_scores,
  test_mats
) {
  if (!nrow(train_assoc)) {
    return(data.frame())
  }

  needed <- c(
    "key",
    "library",
    "pathway",
    "factor",
    "oriented_side",
    "assigned_class_nominal",
    "view_code",
    "feature_original",
    "feature_model",
    "latent_pathway_feature_cosine",
    "expected_observed_sign",
    "rho_spearman_train"
  )

  assert(
    all(needed %in% names(train_assoc)),
    paste(
      "Asociaciones TRAIN incompletas:",
      paste(
        setdiff(
          needed,
          names(train_assoc)
        ),
        collapse = ", "
      )
    )
  )

  rows <- list()
  rr <- 1L

  for (i in seq_len(nrow(train_assoc))) {
    tr <- train_assoc[i, , drop = FALSE]

    score_df <- test_scores[
      test_scores$key == tr$key,
      ,
      drop = FALSE
    ]

    if (!nrow(score_df)) next

    vc <- as.character(
      tr$view_code
    )

    if (!vc %in% names(test_mats)) next

    X <- test_mats[[vc]]

    feature <- as.character(
      tr$feature_original
    )

    present <- feature %in%
      rownames(X)

    rho_test <- NA_real_
    jk <- list(
      median = NA_real_,
      min = NA_real_,
      max = NA_real_,
      sign_stability_vs_full = NA_real_,
      expected_sign_stability = NA_real_
    )

    if (present) {
      x <- as.numeric(
        X[
          feature,
          score_df$sample_id
        ]
      )

      y <- as.numeric(
        score_df$pathway_score
      )

      rho_test <- safe_spearman(
        y,
        x
      )

      jk <- jackknife_spearman(
        y,
        x,
        expected_sign =
          as.numeric(
            tr$expected_observed_sign
          )
      )
    }

    expected_sign <- as.numeric(
      tr$expected_observed_sign
    )

    rho_train <- as.numeric(
      tr$rho_spearman_train
    )

    test_latent_concordance <- if (
      is.finite(rho_test) &&
      abs(rho_test) > NUM_TOL &&
      is.finite(expected_sign)
    ) {
      sign(rho_test) ==
        expected_sign
    } else {
      NA
    }

    train_test_same <- if (
      is.finite(rho_train) &&
      abs(rho_train) > NUM_TOL &&
      is.finite(rho_test) &&
      abs(rho_test) > NUM_TOL
    ) {
      sign(rho_train) ==
        sign(rho_test)
    } else {
      NA
    }

    triple <- if (
      is.finite(expected_sign) &&
      is.finite(rho_train) &&
      abs(rho_train) > NUM_TOL &&
      is.finite(rho_test) &&
      abs(rho_test) > NUM_TOL
    ) {
      sign(rho_train) ==
        expected_sign &&
        sign(rho_test) ==
          expected_sign
    } else {
      NA
    }

    rows[[rr]] <- data.frame(
      family =
        if ("family" %in% names(tr)) {
          tr$family
        } else {
          NA_character_
        },
      key = tr$key,
      library = tr$library,
      pathway = tr$pathway,
      factor = tr$factor,
      oriented_side =
        tr$oriented_side,
      assigned_class_nominal =
        tr$assigned_class_nominal,
      p_assigned_class =
        if (
          "p_assigned_class" %in%
            names(tr)
        ) {
          tr$p_assigned_class
        } else {
          NA_real_
        },
      view =
        if ("view" %in% names(tr)) {
          tr$view
        } else {
          NA_character_
        },
      view_code = vc,
      feature_original =
        feature,
      feature_model =
        tr$feature_model,
      feature_present_in_test =
        present,
      latent_pathway_feature_cosine =
        tr$latent_pathway_feature_cosine,
      expected_observed_sign =
        expected_sign,
      rho_train =
        rho_train,
      rho_test =
        rho_test,
      delta_rho_test_minus_train =
        if (
          is.finite(rho_test) &&
          is.finite(rho_train)
        ) {
          rho_test - rho_train
        } else {
          NA_real_
        },
      expected_aligned_rho_train =
        if (
          is.finite(expected_sign) &&
          is.finite(rho_train)
        ) {
          expected_sign *
            rho_train
        } else {
          NA_real_
        },
      expected_aligned_rho_test =
        if (
          is.finite(expected_sign) &&
          is.finite(rho_test)
        ) {
          expected_sign *
            rho_test
        } else {
          NA_real_
        },
      latent_test_direction_concordant =
        test_latent_concordance,
      train_test_same_sign =
        train_test_same,
      latent_train_test_same_sign =
        triple,
      transport_status =
        transport_status(
          rho_train,
          rho_test
        ),
      jackknife_test_rho_median =
        jk$median,
      jackknife_test_rho_min =
        jk$min,
      jackknife_test_rho_max =
        jk$max,
      jackknife_test_sign_stability_vs_full =
        jk$sign_stability_vs_full,
      jackknife_test_expected_sign_stability =
        jk$expected_sign_stability,
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  if (!length(rows)) {
    return(data.frame())
  }

  do.call(
    rbind,
    rows
  )
}

# ------------------ biomarcador observado por clase -------------------------

candidate_feature_class_audit <- function(
  biomarker_matrix,
  train_mats,
  test_mats,
  train_split,
  test_split
) {
  needed <- c(
    "factor",
    "oriented_side",
    "assigned_class_nominal",
    "view_code",
    "feature_original",
    "feature_model"
  )

  assert(
    all(needed %in% names(biomarker_matrix)),
    paste(
      "Biomarker matrix incompleta:",
      paste(
        setdiff(
          needed,
          names(biomarker_matrix)
        ),
        collapse = ", "
      )
    )
  )

  unique_candidates <- unique(
    biomarker_matrix[
      ,
      needed,
      drop = FALSE
    ]
  )

  rows <- list()
  long_rows <- list()
  rr <- 1L
  ll <- 1L

  for (i in seq_len(nrow(unique_candidates))) {
    cnd <- unique_candidates[i, , drop = FALSE]

    vc <- as.character(
      cnd$view_code
    )

    feature <- as.character(
      cnd$feature_original
    )

    if (
      !vc %in% names(train_mats) ||
      !vc %in% names(test_mats)
    ) {
      next
    }

    present_train <- feature %in%
      rownames(train_mats[[vc]])

    present_test <- feature %in%
      rownames(test_mats[[vc]])

    if (
      !present_train ||
      !present_test
    ) {
      next
    }

    xtr <- as.numeric(
      train_mats[[vc]][
        feature,
        train_split$id
      ]
    )

    xte <- as.numeric(
      test_mats[[vc]][
        feature,
        test_split$id
      ]
    )

    assigned <- as.character(
      cnd$assigned_class_nominal
    )

    med_tr <- vapply(
      GROUPS,
      function(g) {
        stats::median(
          xtr[
            train_split$grupo == g
          ],
          na.rm = TRUE
        )
      },
      numeric(1)
    )

    med_te <- vapply(
      GROUPS,
      function(g) {
        stats::median(
          xte[
            test_split$grupo == g
          ],
          na.rm = TRUE
        )
      },
      numeric(1)
    )

    other <- setdiff(
      GROUPS,
      assigned
    )

    diff_tr <- unname(
      med_tr[assigned] -
        max(
          med_tr[other],
          na.rm = TRUE
        )
    )

    diff_te <- unname(
      med_te[assigned] -
        max(
          med_te[other],
          na.rm = TRUE
        )
    )

    top_tr <- GROUPS[
      order(
        med_tr,
        decreasing = TRUE
      )[1L]
    ]

    top_te <- GROUPS[
      order(
        med_te,
        decreasing = TRUE
      )[1L]
    ]

    rows[[rr]] <- data.frame(
      factor =
        cnd$factor,
      oriented_side =
        cnd$oriented_side,
      assigned_class_nominal =
        assigned,
      view_code =
        vc,
      feature_original =
        feature,
      feature_model =
        cnd$feature_model,
      train_median_NP =
        unname(med_tr["NP"]),
      train_median_FD =
        unname(med_tr["FD"]),
      train_median_SO =
        unname(med_tr["SO"]),
      test_median_NP =
        unname(med_te["NP"]),
      test_median_FD =
        unname(med_te["FD"]),
      test_median_SO =
        unname(med_te["SO"]),
      train_observed_top_class =
        top_tr,
      test_observed_top_class =
        top_te,
      train_assigned_class_is_top =
        top_tr == assigned,
      test_assigned_class_is_top =
        top_te == assigned,
      train_assigned_vs_best_other_difference =
        diff_tr,
      test_assigned_vs_best_other_difference =
        diff_te,
      train_test_candidate_direction_preserved =
        diff_tr > 0 &&
        diff_te > 0,
      candidate_transport_status =
        if (
          diff_tr > 0 &&
          diff_te > 0
        ) {
          "PRESERVED"
        } else if (
          diff_tr > 0 &&
          diff_te <= 0
        ) {
          "TEST_DISCORDANT"
        } else if (
          diff_tr <= 0 &&
          diff_te > 0
        ) {
          "TEST_SUPPORTS"
        } else {
          "NEITHER"
        },
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L

    long_rows[[ll]] <- rbind(
      data.frame(
        factor =
          cnd$factor,
        oriented_side =
          cnd$oriented_side,
        assigned_class_nominal =
          assigned,
        view_code =
          vc,
        feature_original =
          feature,
        feature_model =
          cnd$feature_model,
        dataset = "TRAIN",
        sample_id =
          train_split$id,
        grupo =
          train_split$grupo,
        observed_value =
          xtr,
        stringsAsFactors = FALSE
      ),
      data.frame(
        factor =
          cnd$factor,
        oriented_side =
          cnd$oriented_side,
        assigned_class_nominal =
          assigned,
        view_code =
          vc,
        feature_original =
          feature,
        feature_model =
          cnd$feature_model,
        dataset = "TEST",
        sample_id =
          test_split$id,
        grupo =
          test_split$grupo,
        observed_value =
          xte,
        stringsAsFactors = FALSE
      )
    )

    ll <- ll + 1L
  }

  list(
    summary = if (
      length(rows)
    ) {
      do.call(
        rbind,
        rows
      )
    } else {
      data.frame()
    },
    long = if (
      length(long_rows)
    ) {
      do.call(
        rbind,
        long_rows
      )
    } else {
      data.frame()
    }
  )
}

# ------------------------ pooled descriptivo --------------------------------

pooled_pathway_feature_audit <- function(
  train_assoc,
  test_assoc,
  train_scores,
  test_scores,
  train_mats,
  test_mats
) {
  if (
    !nrow(train_assoc) ||
    !nrow(test_assoc)
  ) {
    return(data.frame())
  }

  key_pair_train <- paste(
    train_assoc$key,
    train_assoc$feature_model,
    sep = "\r"
  )

  key_pair_test <- paste(
    test_assoc$key,
    test_assoc$feature_model,
    sep = "\r"
  )

  common_pairs <- intersect(
    key_pair_train,
    key_pair_test
  )

  rows <- list()
  rr <- 1L

  for (pair in common_pairs) {
    tr <- train_assoc[
      key_pair_train == pair,
      ,
      drop = FALSE
    ]

    te <- test_assoc[
      key_pair_test == pair,
      ,
      drop = FALSE
    ]

    if (
      nrow(tr) != 1L ||
      nrow(te) != 1L
    ) {
      next
    }

    vc <- as.character(
      tr$view_code
    )

    feature <- as.character(
      tr$feature_original
    )

    score_tr <- train_scores[
      train_scores$key == tr$key,
      ,
      drop = FALSE
    ]

    score_te <- test_scores[
      test_scores$key == te$key,
      ,
      drop = FALSE
    ]

    if (
      !nrow(score_tr) ||
      !nrow(score_te)
    ) {
      next
    }

    if (
      !feature %in%
        rownames(
          train_mats[[vc]]
        ) ||
      !feature %in%
        rownames(
          test_mats[[vc]]
        )
    ) {
      next
    }

    xtr <- as.numeric(
      train_mats[[vc]][
        feature,
        score_tr$sample_id
      ]
    )

    xte <- as.numeric(
      test_mats[[vc]][
        feature,
        score_te$sample_id
      ]
    )

    ytr <- as.numeric(
      score_tr$pathway_score
    )

    yte <- as.numeric(
      score_te$pathway_score
    )

    xall <- c(
      xtr,
      xte
    )

    yall <- c(
      ytr,
      yte
    )

    dataset <- c(
      rep(
        "TRAIN",
        length(xtr)
      ),
      rep(
        "TEST",
        length(xte)
      )
    )

    rho_all <- safe_spearman(
      yall,
      xall
    )

    x_centered <- median_center_by_dataset(
      xall,
      dataset
    )

    y_centered <- median_center_by_dataset(
      yall,
      dataset
    )

    rho_dataset_centered <-
      safe_spearman(
        y_centered,
        x_centered
      )

    rows[[rr]] <- data.frame(
      family =
        if ("family" %in% names(tr)) {
          tr$family
        } else {
          NA_character_
        },
      key = tr$key,
      library = tr$library,
      pathway = tr$pathway,
      factor = tr$factor,
      oriented_side =
        tr$oriented_side,
      assigned_class_nominal =
        tr$assigned_class_nominal,
      view_code =
        vc,
      feature_original =
        feature,
      feature_model =
        tr$feature_model,
      latent_pathway_feature_cosine =
        tr$latent_pathway_feature_cosine,
      expected_observed_sign =
        tr$expected_observed_sign,
      rho_train =
        tr$rho_spearman_train,
      rho_test =
        te$rho_test,
      rho_pooled_descriptive =
        rho_all,
      rho_pooled_dataset_centered_descriptive =
        rho_dataset_centered,
      train_test_same_sign =
        te$train_test_same_sign,
      latent_train_test_same_sign =
        te$latent_train_test_same_sign,
      pooled_is_independent_validation =
        FALSE,
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  if (!length(rows)) {
    return(data.frame())
  }

  do.call(
    rbind,
    rows
  )
}

pooled_pathway_class_audit <- function(scores_all) {
  keys <- unique(
    scores_all$key
  )

  rows <- list()
  rr <- 1L

  for (k in keys) {
    z <- scores_all[
      scores_all$key == k,
      ,
      drop = FALSE
    ]

    assigned <- unique(
      z$assigned_class_nominal
    )

    if (
      length(assigned) != 1L ||
      !assigned %in% GROUPS
    ) {
      next
    }

    med <- vapply(
      GROUPS,
      function(g) {
        stats::median(
          z$pathway_score[
            z$grupo == g
          ],
          na.rm = TRUE
        )
      },
      numeric(1)
    )

    other <- setdiff(
      GROUPS,
      assigned
    )

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
      pooled_median_NP =
        unname(med["NP"]),
      pooled_median_FD =
        unname(med["FD"]),
      pooled_median_SO =
        unname(med["SO"]),
      pooled_observed_top_class =
        GROUPS[
          order(
            med,
            decreasing = TRUE
          )[1L]
        ],
      pooled_assigned_vs_best_other_difference =
        unname(
          med[assigned] -
            max(
              med[other],
              na.rm = TRUE
            )
        ),
      pooled_is_independent_validation =
        FALSE,
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  if (!length(rows)) {
    return(data.frame())
  }

  do.call(
    rbind,
    rows
  )
}

# --------------------- resumen para matriz final ----------------------------

summarize_test_links_by_candidate <- function(test_assoc) {
  if (!nrow(test_assoc)) {
    return(data.frame())
  }

  key <- paste(
    test_assoc$factor,
    test_assoc$oriented_side,
    test_assoc$feature_model,
    sep = "\r"
  )

  split_idx <- split(
    seq_len(nrow(test_assoc)),
    key
  )

  rows <- list()
  rr <- 1L

  for (idx in split_idx) {
    z <- test_assoc[
      idx,
      ,
      drop = FALSE
    ]

    rows[[rr]] <- data.frame(
      factor = z$factor[1L],
      oriented_side =
        z$oriented_side[1L],
      feature_model =
        z$feature_model[1L],
      test_n_pathway_links_audited =
        nrow(z),
      test_n_latent_direction_concordant =
        sum(
          z$latent_test_direction_concordant,
          na.rm = TRUE
        ),
      test_prop_latent_direction_concordant =
        if (
          any(
            !is.na(
              z$latent_test_direction_concordant
            )
          )
        ) {
          mean(
            z$latent_test_direction_concordant,
            na.rm = TRUE
          )
        } else {
          NA_real_
        },
      test_n_train_test_same_sign =
        sum(
          z$train_test_same_sign,
          na.rm = TRUE
        ),
      test_prop_train_test_same_sign =
        if (
          any(
            !is.na(
              z$train_test_same_sign
            )
          )
        ) {
          mean(
            z$train_test_same_sign,
            na.rm = TRUE
          )
        } else {
          NA_real_
        },
      test_n_latent_train_test_same_sign =
        sum(
          z$latent_train_test_same_sign,
          na.rm = TRUE
        ),
      test_prop_latent_train_test_same_sign =
        if (
          any(
            !is.na(
              z$latent_train_test_same_sign
            )
          )
        ) {
          mean(
            z$latent_train_test_same_sign,
            na.rm = TRUE
          )
        } else {
          NA_real_
        },
      test_median_expected_aligned_rho =
        if (
          any(
            is.finite(
              z$expected_aligned_rho_test
            )
          )
        ) {
          stats::median(
            z$expected_aligned_rho_test,
            na.rm = TRUE
          )
        } else {
          NA_real_
        },
      test_median_jackknife_expected_sign_stability =
        if (
          any(
            is.finite(
              z$jackknife_test_expected_sign_stability
            )
          )
        ) {
          stats::median(
            z$jackknife_test_expected_sign_stability,
            na.rm = TRUE
          )
        } else {
          NA_real_
        },
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  do.call(
    rbind,
    rows
  )
}

summarize_test_pathway_class_by_module <- function(test_class) {
  if (!nrow(test_class)) {
    return(data.frame())
  }

  key <- paste(
    test_class$factor,
    test_class$oriented_side,
    sep = "\r"
  )

  split_idx <- split(
    seq_len(nrow(test_class)),
    key
  )

  rows <- list()
  rr <- 1L

  for (idx in split_idx) {
    z <- test_class[
      idx,
      ,
      drop = FALSE
    ]

    rows[[rr]] <- data.frame(
      factor = z$factor[1L],
      oriented_side =
        z$oriented_side[1L],
      test_n_pathways_class_audited =
        nrow(z),
      test_prop_pathways_assigned_class_observed_top =
        mean(
          z$assigned_class_is_observed_top,
          na.rm = TRUE
        ),
      train_test_n_pathways_class_top_both =
        sum(
          z$train_test_assigned_class_top_both,
          na.rm = TRUE
        ),
      train_test_prop_pathways_class_top_both =
        mean(
          z$train_test_assigned_class_top_both,
          na.rm = TRUE
        ),
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  do.call(
    rbind,
    rows
  )
}

attach_by_key <- function(base, extra, key_cols) {
  if (!nrow(extra)) {
    return(base)
  }

  key_base <- do.call(
    paste,
    c(
      base[key_cols],
      sep = "\r"
    )
  )

  key_extra <- do.call(
    paste,
    c(
      extra[key_cols],
      sep = "\r"
    )
  )

  idx <- match(
    key_base,
    key_extra
  )

  add_cols <- setdiff(
    names(extra),
    key_cols
  )

  for (cc in add_cols) {
    base[[cc]] <-
      extra[[cc]][idx]
  }

  base
}

# ----------------------------- una familia ----------------------------------

run_family <- function(selection_run, family, stamp) {
  message("\n============================================================")
  message("TEST TRANSPORTABILITY + POOLED AUDIT | ", family)
  message("============================================================")

  bidir_run <- find_bidirectional_run(
    selection_run,
    family
  )

  train_run <- find_observed_train_run(
    bidir_run
  )

  contract_path <- file.path(
    train_run,
    "rds",
    "TEST_PROJECTION_CONTRACT_FROZEN.rds"
  )

  contract <- readRDS(
    contract_path
  )

  assert(
    identical(
      as.character(
        contract$family
      ),
      family
    ),
    "Familia del contrato != familia actual"
  )

  required_contract <- c(
    "pathway_scoring_spec",
    "pathway_profiles",
    "pathway_feature_links",
    "pathway_class_audit_train",
    "pathway_feature_associations_train",
    "biomarker_evidence_matrix_train"
  )

  assert(
    all(
      required_contract %in%
        names(contract)
    ),
    paste(
      "Contrato 6A incompleto:",
      paste(
        setdiff(
          required_contract,
          names(contract)
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

  train_split <- get_split_table(
    rfm,
    "train"
  )

  test_split <- get_split_table(
    rfm,
    "test"
  )

  assert(
    !length(
      intersect(
        train_split$id,
        test_split$id
      )
    ),
    "TRAIN y TEST comparten IDs"
  )

  validate_train_contract(
    contract,
    train_split
  )

  train_mats <- load_all_views(
    rfm,
    "train",
    train_split
  )

  test_mats <- load_all_views(
    rfm,
    "test",
    test_split
  )

  # Features TRAIN/TEST deben ser las mismas por vista.
  for (v in names(train_mats)) {
    assert(
      setequal(
        rownames(
          train_mats[[v]]
        ),
        rownames(
          test_mats[[v]]
        )
      ),
      paste(
        "Features TRAIN/TEST distintas en",
        v
      )
    )
  }

  message(
    "TEST: ",
    paste(
      paste0(
        names(
          table(
            factor(
              test_split$grupo,
              levels = GROUPS
            )
          )
        ),
        "=",
        as.integer(
          table(
            factor(
              test_split$grupo,
              levels = GROUPS
            )
          )
        )
      ),
      collapse = " | "
    )
  )

  # -------------------------------------------------------------------------
  # 1) Proyección congelada de pathway scores a TRAIN y TEST.
  #    TRAIN se recalcula con el mismo contrato solo para poder construir
  #    la auditoría pooled bajo exactamente la misma función.
  # -------------------------------------------------------------------------
  train_scores_obj <- project_pathway_scores(
    scoring_spec =
      contract$pathway_scoring_spec,
    profiles =
      contract$pathway_profiles,
    Xtx =
      train_mats$tx,
    sample_ids =
      train_split$id,
    labels =
      factor(
        train_split$grupo,
        levels = GROUPS
      ),
    dataset_label =
      "TRAIN"
  )

  test_scores_obj <- project_pathway_scores(
    scoring_spec =
      contract$pathway_scoring_spec,
    profiles =
      contract$pathway_profiles,
    Xtx =
      test_mats$tx,
    sample_ids =
      test_split$id,
    labels =
      factor(
        test_split$grupo,
        levels = GROUPS
      ),
    dataset_label =
      "TEST"
  )

  train_scores <-
    train_scores_obj$scores

  test_scores <-
    test_scores_obj$scores

  # -------------------------------------------------------------------------
  # 2) Pathway observado <-> clase en TEST.
  # -------------------------------------------------------------------------
  test_pathway_class <-
    audit_pathway_class_direction(
      scores =
        test_scores,
      train_class_audit =
        contract$pathway_class_audit_train
    )

  # -------------------------------------------------------------------------
  # 3) Pathway observado <-> PR/ME/CL en TEST.
  #    Parejas EXACTAS congeladas por 6A.
  # -------------------------------------------------------------------------
  test_assoc <-
    audit_pathway_feature_test(
      train_assoc =
        contract$pathway_feature_associations_train,
      test_scores =
        test_scores,
      test_mats =
        test_mats
    )

  # -------------------------------------------------------------------------
  # 4) Biomarcadores candidatos: dirección observada por clase en TEST.
  # -------------------------------------------------------------------------
  candidate_audit <-
    candidate_feature_class_audit(
      biomarker_matrix =
        contract$biomarker_evidence_matrix_train,
      train_mats =
        train_mats,
      test_mats =
        test_mats,
      train_split =
        train_split,
      test_split =
        test_split
    )

  candidate_class_summary <-
    candidate_audit$summary

  candidate_long <-
    candidate_audit$long

  # -------------------------------------------------------------------------
  # 5) POOLED descriptivo.
  # -------------------------------------------------------------------------
  pooled_scores <- rbind(
    train_scores,
    test_scores
  )

  pooled_assoc <-
    pooled_pathway_feature_audit(
      train_assoc =
        contract$pathway_feature_associations_train,
      test_assoc =
        test_assoc,
      train_scores =
        train_scores,
      test_scores =
        test_scores,
      train_mats =
        train_mats,
      test_mats =
        test_mats
    )

  pooled_class <-
    pooled_pathway_class_audit(
      pooled_scores
    )

  # -------------------------------------------------------------------------
  # 6) Integrar TEST a la matriz maestra de biomarcadores.
  # -------------------------------------------------------------------------
  final_matrix <-
    contract$biomarker_evidence_matrix_train

  test_link_summary <-
    summarize_test_links_by_candidate(
      test_assoc
    )

  test_module_summary <-
    summarize_test_pathway_class_by_module(
      test_pathway_class
    )

  final_matrix <- attach_by_key(
    base = final_matrix,
    extra =
      candidate_class_summary,
    key_cols = c(
      "factor",
      "oriented_side",
      "feature_model"
    )
  )

  final_matrix <- attach_by_key(
    base = final_matrix,
    extra =
      test_link_summary,
    key_cols = c(
      "factor",
      "oriented_side",
      "feature_model"
    )
  )

  final_matrix <- attach_by_key(
    base = final_matrix,
    extra =
      test_module_summary,
    key_cols = c(
      "factor",
      "oriented_side"
    )
  )

  # -------------------------------------------------------------------------
  # Salidas.
  # -------------------------------------------------------------------------
  out <- file.path(
    bidir_run,
    "TEST_TRANSPORTABILITY_AUDIT",
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
    observed_train_run =
      train_run,
    frozen_contract =
      contract_path,
    ready_for_modeling =
      rfm_path,
    n_train =
      nrow(train_split),
    n_test =
      nrow(test_split),
    test_NP =
      sum(
        test_split$grupo == "NP"
      ),
    test_FD =
      sum(
        test_split$grupo == "FD"
      ),
    test_SO =
      sum(
        test_split$grupo == "SO"
      ),
    n_frozen_pathways =
      length(
        unique(
          test_scores$key
        )
      ),
    n_frozen_pathway_feature_pairs =
      nrow(
        contract$pathway_feature_associations_train
      ),
    n_test_pathway_feature_pairs_evaluated =
      nrow(test_assoc),
    test_used_for_selection =
      FALSE,
    pathways_reselected_in_test =
      FALSE,
    genes_reselected_in_test =
      FALSE,
    weights_refit_in_test =
      FALSE,
    train_scaling_reused_for_test =
      TRUE,
    pooled_is_independent_validation =
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
    test_scores_obj$qc,
    file.path(
      out,
      "tables",
      "01_TEST_PATHWAY_SCORE_PROJECTION_QC.csv"
    )
  )

  write_csv(
    test_scores,
    file.path(
      out,
      "tables",
      "02_PATHWAY_SCORES_TEST.csv"
    )
  )

  write_csv(
    test_pathway_class,
    file.path(
      out,
      "tables",
      "03_PATHWAY_CLASS_DIRECTION_TEST.csv"
    )
  )

  write_csv(
    test_assoc,
    file.path(
      out,
      "tables",
      "04_PATHWAY_FEATURE_ASSOCIATIONS_TEST.csv"
    )
  )

  write_csv(
    candidate_class_summary,
    file.path(
      out,
      "tables",
      "05_BIOMARKER_CLASS_DIRECTION_TRAIN_TEST.csv"
    )
  )

  write_csv(
    test_link_summary,
    file.path(
      out,
      "tables",
      "06_BIOMARKER_TEST_PATHWAY_LINK_SUMMARY.csv"
    )
  )

  write_csv(
    test_module_summary,
    file.path(
      out,
      "tables",
      "06B_MODULE_PATHWAY_CLASS_TEST_SUMMARY.csv"
    )
  )

  write_csv(
    final_matrix,
    file.path(
      out,
      "tables",
      "07_BIOMARKER_EVIDENCE_MATRIX_FINAL.csv"
    )
  )

  # Auditoría pooled, explícitamente separada.
  write_csv(
    pooled_scores,
    file.path(
      out,
      "tables",
      "08_POOLED_PATHWAY_SCORES_DESCRIPTIVE.csv"
    )
  )

  write_csv(
    candidate_long,
    file.path(
      out,
      "tables",
      "09_POOLED_BIOMARKER_VALUES_LONG_DESCRIPTIVE.csv"
    )
  )

  write_csv(
    pooled_assoc,
    file.path(
      out,
      "tables",
      "10_POOLED_PATHWAY_FEATURE_ASSOCIATIONS_DESCRIPTIVE.csv"
    )
  )

  write_csv(
    pooled_class,
    file.path(
      out,
      "tables",
      "11_POOLED_PATHWAY_CLASS_DESCRIPTIVE.csv"
    )
  )

  saveRDS(
    list(
      family = family,
      input_audit =
        input_audit,
      frozen_contract =
        contract,
      pathway_scores_train_frozen_function =
        train_scores,
      pathway_scores_test =
        test_scores,
      pathway_class_test =
        test_pathway_class,
      pathway_feature_test =
        test_assoc,
      biomarker_class_train_test =
        candidate_class_summary,
      biomarker_values_pooled_long =
        candidate_long,
      biomarker_test_pathway_summary =
        test_link_summary,
      module_pathway_test_summary =
        test_module_summary,
      biomarker_evidence_matrix_final =
        final_matrix,
      pooled = list(
        pathway_scores =
          pooled_scores,
        pathway_feature_associations =
          pooled_assoc,
        pathway_class =
          pooled_class,
        independent_validation =
          FALSE
      )
    ),
    file.path(
      out,
      "rds",
      "TEST_TRANSPORTABILITY_AUDIT_OBJECT.rds"
    ),
    version = 3
  )

  readme <- c(
    paste(
      "METATEST |",
      family,
      "| PASO 6B TEST TRANSPORTABILITY AUDIT"
    ),
    "",
    "OBJETIVO:",
    "- Auditar en TEST relaciones completamente congeladas antes de abrir TEST.",
    "- Añadir una auditoría pooled TRAIN+TEST estrictamente descriptiva.",
    "",
    "CONTRATO:",
    paste(
      "-",
      contract_path
    ),
    "- Pathways, genes, weights, relaciones pathway-feature y candidatos",
    "  proceden del PASO 6A.",
    "- TEST no modifica ninguna de esas decisiones.",
    "",
    "PATHWAY SCORE TEST:",
    "- Usa los mismos genes congelados.",
    "- Usa exactamente mean/SD TRAIN.",
    "- Usa exactamente normalized_anchor_weight TRAIN.",
    "- No hay normalización usando TEST.",
    "",
    "PATHWAY -> CLASE TEST:",
    "- Se comparan medianas NP/FD/SO.",
    "- La salida principal es si la clase asignada por el bridge sigue siendo",
    "  la clase con mayor score observado.",
    "- No se usa p-value TEST como criterio de selección.",
    "",
    "PATHWAY -> FEATURE TEST:",
    "- Solo se auditan las parejas PR/ME/CL ya probadas en TRAIN.",
    "- Se reporta rho TEST, delta respecto a TRAIN y tres concordancias:",
    "    latent -> TEST",
    "    TRAIN -> TEST",
    "    latent -> TRAIN -> TEST",
    "- Jackknife TEST se reporta como estabilidad descriptiva.",
    "",
    "BIOMARCADOR -> CLASE TEST:",
    "- Para cada candidato ya congelado se reportan medianas por clase en",
    "  TRAIN y TEST.",
    "- Se indica si la clase asignada al módulo sigue ocupando el extremo alto.",
    "- Un resultado discordante TEST no elimina automáticamente el candidato.",
    "",
    "POOLED TRAIN+TEST:",
    "- Es únicamente una auditoría descriptiva n=26.",
    "- Mantiene dataset=TRAIN/TEST en las tablas long.",
    "- Todos los pathway scores siguen usando escalado TRAIN.",
    "- rho_pooled_descriptive NO es validación independiente.",
    "- rho_pooled_dataset_centered_descriptive reduce desplazamientos globales",
    "  TRAIN/TEST y se reporta solo como descriptor.",
    "",
    "TABLA FINAL:",
    "- 07_BIOMARKER_EVIDENCE_MATRIX_FINAL.csv",
    "  integra evidencia latente, pathways, bidireccional, observado TRAIN",
    "  y transportabilidad TEST.",
    "- No se calcula un score compuesto obligatorio ni un ranking automático.",
    "",
    "FIN DEL DESARROLLO ANALÍTICO:",
    "- Después de PASO 6B no se añade otro modelo.",
    "- La siguiente etapa es síntesis biológica, priorización narrativa y plots:",
    "  biomarcadores, mecanismos clase-factor-pathway-ómica, concordancia",
    "  TRAIN/TEST y comparación descriptiva P1/P2/P3."
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
    "TEST_TRANSPORTABILITY_AUDIT"
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
    n_test_pathways =
      length(
        unique(
          test_scores$key
        )
      ),
    n_test_pathways_assigned_class_top =
      sum(
        test_pathway_class$assigned_class_is_observed_top,
        na.rm = TRUE
      ),
    n_test_crossomic_links =
      nrow(test_assoc),
    n_test_crossomic_links_latent_concordant =
      if (nrow(test_assoc)) {
        sum(
          test_assoc$latent_test_direction_concordant,
          na.rm = TRUE
        )
      } else {
        0L
      },
    n_test_crossomic_links_train_test_same_sign =
      if (nrow(test_assoc)) {
        sum(
          test_assoc$train_test_same_sign,
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
          n_test_pathways =
            NA_integer_,
          n_test_pathways_assigned_class_top =
            NA_integer_,
          n_test_crossomic_links =
            NA_integer_,
          n_test_crossomic_links_latent_concordant =
            NA_integer_,
          n_test_crossomic_links_train_test_same_sign =
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
      "TEST_TRANSPORTABILITY_AUDIT_RUN_",
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

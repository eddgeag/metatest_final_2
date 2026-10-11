#!/usr/bin/env Rscript

# ==============================================================================
# METATEST — 19_probe_final_small_panel_console.R
#
# PRUEBA SOLAMENTE:
#   * NO escribe CSV
#   * NO guarda figuras
#   * NO modifica ninguna selección
#   * imprime candidatos/resúmenes PCA en terminal
#   * ABRE los PCA PC1-PC2 en ventanas gráficas en pantalla
#
# Working rule for a SMALL exploratory final panel:
#   1) CORE: feature strict Evidence-A presente en P1 + P2 + P3
#   2) VIP >= 1 en los 3 PLS-DA
#   3) conservar el cuartil superior (Q75) de
#      median_abs_loading_across_models
#
# SGL se imprime como corroboración, pero NO filtra.
#
# PCA:
#   * se ajusta SOLO con TRAIN
#   * TEST se proyecta con transformación TRAIN congelada
#   * se muestran las 26 observaciones en el resumen descriptivo
# ==============================================================================

options(stringsAsFactors = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.", call. = FALSE)
}

UPSTREAM_DIR <- Sys.getenv(
  "UPSTREAM_DIR",
  unset = path.expand(
    "~/Documents/PhD/metatest/results_upstream_local_20260908/comparacion_modelamientos/upstream_comun"
  )
)

FINAL_CONSENSUS_RUN <- Sys.getenv("FINAL_CONSENSUS_RUN", unset = "")
SMALL_QUANTILE <- as.numeric(
  Sys.getenv("PANEL19_SMALL_QUANTILE", unset = "0.75")
)

CLASS_LEVELS <- c("NP", "FD", "SO")
VIEW_LEVELS <- c("tx", "pr", "me", "cl")

stopifnot(
  is.finite(SMALL_QUANTILE),
  SMALL_QUANTILE > 0,
  SMALL_QUANTILE < 1
)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

assert <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
}

canonical <- function(path, must_work = TRUE) {
  normalizePath(
    path.expand(path),
    winslash = "/",
    mustWork = must_work
  )
}

read_csv_checked <- function(path, required = character()) {
  assert(file.exists(path), paste("Missing:", path))

  x <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  miss <- setdiff(required, names(x))
  assert(
    !length(miss),
    paste(
      "Missing columns in",
      basename(path),
      ":",
      paste(miss, collapse = ", ")
    )
  )

  x
}

find_latest_selection <- function(upstream) {
  root <- file.path(
    upstream,
    "downstream",
    "seleccion_mejor_escenario_train"
  )

  assert(dir.exists(root), paste("Missing selection root:", root))

  dirs <- list.dirs(
    root,
    recursive = FALSE,
    full.names = TRUE
  )

  dirs <- dirs[
    file.exists(
      file.path(
        dirs,
        "SELECTED_SCENARIOS_BY_ANALYSIS.csv"
      )
    )
  ]

  assert(length(dirs) > 0L, "No selection run found.")
  info <- file.info(dirs)

  canonical(
    dirs[
      order(
        info$mtime,
        basename(dirs),
        decreasing = TRUE
      )
    ][1L]
  )
}

find_latest_bio <- function(selection) {
  root <- file.path(selection, "BIOLOGICAL_PROGRAMS")
  assert(dir.exists(root), paste("Missing BIOLOGICAL_PROGRAMS:", root))

  latest <- file.path(root, "LATEST_RUN.txt")

  if (file.exists(latest)) {
    p <- trimws(readLines(latest, warn = FALSE)[1L])
    if (nzchar(p) && dir.exists(p)) return(canonical(p))
  }

  dirs <- list.dirs(
    root,
    recursive = FALSE,
    full.names = TRUE
  )

  assert(length(dirs) > 0L, "No BIOLOGICAL_PROGRAMS run found.")
  info <- file.info(dirs)

  canonical(
    dirs[
      order(
        info$mtime,
        basename(dirs),
        decreasing = TRUE
      )
    ][1L]
  )
}

find_consensus <- function(bio) {
  if (nzchar(FINAL_CONSENSUS_RUN)) {
    return(canonical(FINAL_CONSENSUS_RUN))
  }

  root <- file.path(bio, "FINAL_CONSENSUS_PANEL")
  assert(dir.exists(root), paste("Missing FINAL_CONSENSUS_PANEL:", root))

  latest <- file.path(root, "LATEST_RUN.txt")

  if (file.exists(latest)) {
    p <- trimws(readLines(latest, warn = FALSE)[1L])

    if (
      nzchar(p) &&
      file.exists(
        file.path(
          p,
          "tables",
          "00_FINAL_PANEL_MASTER.csv"
        )
      )
    ) {
      return(canonical(p))
    }
  }

  dirs <- list.dirs(
    root,
    recursive = FALSE,
    full.names = TRUE
  )

  dirs <- dirs[
    file.exists(
      file.path(
        dirs,
        "tables",
        "00_FINAL_PANEL_MASTER.csv"
      )
    )
  ]

  assert(length(dirs) > 0L, "No complete FINAL_CONSENSUS_PANEL run found.")
  info <- file.info(dirs)

  canonical(
    dirs[
      order(
        info$mtime,
        basename(dirs),
        decreasing = TRUE
      )
    ][1L]
  )
}

extract_group_table <- function(g) {
  if (is.data.frame(g)) {
    assert(
      all(c("id", "grupo") %in% names(g)),
      "Group table must contain id and grupo."
    )

    return(
      data.frame(
        sample_id = as.character(g$id),
        grupo = factor(
          as.character(g$grupo),
          levels = CLASS_LEVELS
        ),
        stringsAsFactors = FALSE
      )
    )
  }

  if (is.atomic(g) && !is.null(names(g))) {
    return(
      data.frame(
        sample_id = names(g),
        grupo = factor(
          as.character(g),
          levels = CLASS_LEVELS
        ),
        stringsAsFactors = FALSE
      )
    )
  }

  stop("Unrecognized group structure.", call. = FALSE)
}

eta2_oneway <- function(score, group) {
  ok <- is.finite(score) & !is.na(group)
  score <- score[ok]
  group <- droplevels(factor(group[ok]))

  if (length(score) < 3L || nlevels(group) < 2L) {
    return(NA_real_)
  }

  grand <- mean(score)
  ss_total <- sum((score - grand)^2)

  if (!is.finite(ss_total) || ss_total <= 0) {
    return(NA_real_)
  }

  ss_between <- sum(
    vapply(
      levels(group),
      function(g) {
        x <- score[group == g]
        length(x) * (mean(x) - grand)^2
      },
      numeric(1)
    )
  )

  ss_between / ss_total
}

mean_silhouette_2d <- function(scores, group) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    return(NA_real_)
  }

  scores <- as.matrix(scores)
  ok <- stats::complete.cases(scores) & !is.na(group)

  scores <- scores[ok, , drop = FALSE]
  group <- droplevels(factor(group[ok]))

  if (
    nrow(scores) < 4L ||
    nlevels(group) < 2L ||
    any(table(group) < 2L)
  ) {
    return(NA_real_)
  }

  sil <- cluster::silhouette(
    as.integer(group),
    stats::dist(scores)
  )

  mean(sil[, "sil_width"])
}

build_matrix_from_p1 <- function(rfm, features, train_meta, test_meta) {
  train_cols <- list()
  test_cols <- list()
  skipped <- character()

  for (i in seq_len(nrow(features))) {
    view <- tolower(as.character(features$view_code[i]))
    feature <- as.character(features$feature_original[i])

    key <- paste(
      toupper(view),
      features$feature_label[i],
      sep = " | "
    )

    if (
      !view %in% names(rfm) ||
      is.null(rfm[[view]]$train) ||
      is.null(rfm[[view]]$test)
    ) {
      skipped <- c(skipped, paste0(key, " [view absent]"))
      next
    }

    tr <- as.matrix(rfm[[view]]$train)
    te <- as.matrix(rfm[[view]]$test)

    if (!feature %in% rownames(tr)) {
      skipped <- c(skipped, paste0(key, " [feature absent TRAIN]"))
      next
    }

    a <- as.numeric(
      tr[
        feature,
        train_meta$sample_id
      ]
    )

    if (any(!is.finite(a))) {
      skipped <- c(skipped, paste0(key, " [nonfinite TRAIN]"))
      next
    }

    if (!is.finite(stats::sd(a)) || stats::sd(a) <= 0) {
      skipped <- c(skipped, paste0(key, " [zero variance TRAIN]"))
      next
    }

    b <- if (feature %in% rownames(te)) {
      as.numeric(
        te[
          feature,
          test_meta$sample_id
        ]
      )
    } else {
      rep(
        NA_real_,
        nrow(test_meta)
      )
    }

    if (any(!is.finite(b))) {
      b[
        !is.finite(b)
      ] <- stats::median(a)
    }

    machine_key <- paste(
      view,
      feature,
      sep = "::"
    )

    train_cols[[machine_key]] <- a
    test_cols[[machine_key]] <- b
  }

  if (!length(train_cols)) {
    return(
      list(
        train = matrix(
          numeric(),
          nrow = nrow(train_meta),
          ncol = 0
        ),
        test = matrix(
          numeric(),
          nrow = nrow(test_meta),
          ncol = 0
        ),
        skipped = skipped
      )
    )
  }

  xtr <- do.call(cbind, train_cols)
  xte <- do.call(cbind, test_cols)

  rownames(xtr) <- train_meta$sample_id
  rownames(xte) <- test_meta$sample_id

  list(
    train = xtr,
    test = xte,
    skipped = skipped
  )
}

open_screen_device <- function(title = "METATEST PCA") {
  sys <- Sys.info()[["sysname"]]

  if (identical(sys, "Linux")) {
    if (!nzchar(Sys.getenv("DISPLAY"))) {
      stop(
        paste0(
          "No graphical DISPLAY is available. ",
          "Run this script from your graphical Ubuntu session, not through a headless SSH shell."
        ),
        call. = FALSE
      )
    }

    grDevices::x11(
      width = 7.5,
      height = 6.5,
      title = title
    )
  } else if (identical(sys, "Windows")) {
    grDevices::windows(
      width = 7.5,
      height = 6.5,
      title = title
    )
  } else if (identical(sys, "Darwin")) {
    grDevices::quartz(
      width = 7.5,
      height = 6.5,
      title = title
    )
  } else {
    grDevices::dev.new(
      width = 7.5,
      height = 6.5
    )
  }
}

print_pca_summary <- function(name, feature_df, rfm, train_meta, test_meta) {
  cat("\n------------------------------------------------------------\n")
  cat("PCA SET:", name, "\n")
  cat("Requested features:", nrow(feature_df), "\n")

  if (nrow(feature_df) < 2L) {
    cat("SKIP: fewer than 2 features; PCA is not informative.\n")
    return(invisible(NULL))
  }

  mat <- build_matrix_from_p1(
    rfm = rfm,
    features = feature_df,
    train_meta = train_meta,
    test_meta = test_meta
  )

  cat("Usable features in P1 observed matrix:", ncol(mat$train), "\n")

  if (length(mat$skipped)) {
    cat("Skipped:", length(mat$skipped), "\n")
    cat(
      paste0("  - ", mat$skipped),
      sep = "\n"
    )
    cat("\n")
  }

  if (ncol(mat$train) < 2L) {
    cat("SKIP: fewer than 2 usable observed features.\n")
    return(invisible(NULL))
  }

  pca <- stats::prcomp(
    mat$train,
    center = TRUE,
    scale. = TRUE
  )

  tr_score <- as.data.frame(pca$x)

  te_score <- as.data.frame(
    stats::predict(
      pca,
      newdata = mat$test
    )
  )

  tr_score$sample_id <- rownames(mat$train)
  te_score$sample_id <- rownames(mat$test)

  tr_score$grupo <- train_meta$grupo[
    match(
      tr_score$sample_id,
      train_meta$sample_id
    )
  ]

  te_score$grupo <- test_meta$grupo[
    match(
      te_score$sample_id,
      test_meta$sample_id
    )
  ]

  tr_score$split <- "TRAIN"
  te_score$split <- "TEST"

  all_score <- rbind(
    tr_score,
    te_score
  )

  var <- pca$sdev^2
  var <- var / sum(var)

  nshow <- min(3L, length(var))

  cat(
    "TRAIN-fit variance explained:",
    paste(
      paste0(
        "PC",
        seq_len(nshow),
        "=",
        sprintf("%.1f%%", 100 * var[seq_len(nshow)])
      ),
      collapse = " | "
    ),
    "\n"
  )

  eta1_all <- eta2_oneway(
    all_score$PC1,
    all_score$grupo
  )

  eta2_all <- if ("PC2" %in% names(all_score)) {
    eta2_oneway(
      all_score$PC2,
      all_score$grupo
    )
  } else {
    NA_real_
  }

  eta1_train <- eta2_oneway(
    tr_score$PC1,
    tr_score$grupo
  )

  eta2_train <- if ("PC2" %in% names(tr_score)) {
    eta2_oneway(
      tr_score$PC2,
      tr_score$grupo
    )
  } else {
    NA_real_
  }

  sil_all <- if ("PC2" %in% names(all_score)) {
    mean_silhouette_2d(
      all_score[, c("PC1", "PC2"), drop = FALSE],
      all_score$grupo
    )
  } else {
    NA_real_
  }

  sil_train <- if ("PC2" %in% names(tr_score)) {
    mean_silhouette_2d(
      tr_score[, c("PC1", "PC2"), drop = FALSE],
      tr_score$grupo
    )
  } else {
    NA_real_
  }

  cat(
    "TRAIN eta^2: PC1=",
    sprintf("%.3f", eta1_train),
    " | PC2=",
    sprintf("%.3f", eta2_train),
    "\n",
    sep = ""
  )

  cat(
    "ALL-26 descriptive eta^2: PC1=",
    sprintf("%.3f", eta1_all),
    " | PC2=",
    sprintf("%.3f", eta2_all),
    "\n",
    sep = ""
  )

  cat(
    "Mean silhouette PC1-PC2: TRAIN=",
    sprintf("%.3f", sil_train),
    " | ALL-26 descriptive=",
    sprintf("%.3f", sil_all),
    "\n",
    sep = ""
  )

  if ("PC2" %in% names(all_score)) {
    plot_df <- all_score
    plot_df$grupo <- factor(
      plot_df$grupo,
      levels = CLASS_LEVELS
    )

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = PC1,
        y = PC2,
        colour = grupo,
        shape = split
      )
    ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linewidth = 0.3,
        linetype = 3
      ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linewidth = 0.3,
        linetype = 3
      ) +
      ggplot2::geom_point(
        size = 3.2,
        alpha = 0.90
      ) +
      ggplot2::labs(
        title = paste0(
          "PCA — ",
          name
        ),
        subtitle = paste0(
          "TRAIN-fit; TEST projected | ",
          "eta2 ALL: PC1=",
          sprintf("%.3f", eta1_all),
          ", PC2=",
          sprintf("%.3f", eta2_all),
          " | silhouette ALL=",
          sprintf("%.3f", sil_all)
        ),
        x = paste0(
          "PC1 (",
          sprintf("%.1f%%", 100 * var[1L]),
          ")"
        ),
        y = paste0(
          "PC2 (",
          sprintf("%.1f%%", 100 * var[2L]),
          ")"
        ),
        colour = "Phenotype",
        shape = "Split"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      )

    open_screen_device(
      paste0(
        "METATEST PCA - ",
        name
      )
    )

    print(p)
  }

  invisible(
    list(
      pca = pca,
      train_scores = tr_score,
      test_scores = te_score
    )
  )
}

# ------------------------------------------------------------------------------
# Locate latest final consensus
# ------------------------------------------------------------------------------

SELECTION <- find_latest_selection(UPSTREAM_DIR)
BIO <- find_latest_bio(SELECTION)
CONSENSUS <- find_consensus(BIO)

master <- read_csv_checked(
  file.path(
    CONSENSUS,
    "tables",
    "00_FINAL_PANEL_MASTER.csv"
  ),
  required = c(
    "view_code",
    "feature_original",
    "feature_label",
    "consensus_tier",
    "n_VIP_ge_threshold",
    "VIP_median_across_present",
    "median_abs_loading_across_models",
    "n_SGL_selected",
    "families"
  )
)

master$view_code <- tolower(as.character(master$view_code))
master$consensus_tier <- as.character(master$consensus_tier)

cat("============================================================\n")
cat("METATEST — SMALL FINAL PANEL PROBE (CONSOLE ONLY)\n")
cat("Consensus run:", CONSENSUS, "\n")
cat("============================================================\n")

# ------------------------------------------------------------------------------
# Current final-panel structure
# ------------------------------------------------------------------------------

cat("\nCURRENT CONSENSUS\n")
cat("CORE      :", sum(master$consensus_tier == "CORE"), "\n")
cat("SHARED_2  :", sum(master$consensus_tier == "SHARED_2"), "\n")
cat("SPECIFIC  :", sum(master$consensus_tier == "SPECIFIC"), "\n")

core_vip3 <- master[
  master$consensus_tier == "CORE" &
    master$n_VIP_ge_threshold == 3L,
  ,
  drop = FALSE
]

cat("\nCORE + VIP>=1 in all 3 PLS-DA:", nrow(core_vip3), "\n")

# ------------------------------------------------------------------------------
# Sensitivity: contribution threshold
# ------------------------------------------------------------------------------

q_grid <- c(0.67, 0.75, 0.80, 0.85)

cat("\nSMALL-PANEL SENSITIVITY\n")
cat(
  "Eligibility is fixed: CORE + VIP>=1 in P1/P2/P3.\n",
  "Only the loading-contribution quantile changes.\n",
  sep = ""
)

for (q in q_grid) {
  thr <- stats::quantile(
    core_vip3$median_abs_loading_across_models,
    probs = q,
    na.rm = TRUE,
    names = FALSE
  )

  z <- core_vip3[
    core_vip3$median_abs_loading_across_models >= thr,
    ,
    drop = FALSE
  ]

  views <- table(
    factor(
      toupper(z$view_code),
      levels = toupper(VIEW_LEVELS)
    )
  )

  cat(
    sprintf(
      "Q%.0f -> n=%d | TX=%d PR=%d ME=%d CL=%d | threshold=%.5f\n",
      100 * q,
      nrow(z),
      views["TX"],
      views["PR"],
      views["ME"],
      views["CL"],
      thr
    )
  )
}

# ------------------------------------------------------------------------------
# Working small panel: Q75
# ------------------------------------------------------------------------------

small_threshold <- stats::quantile(
  core_vip3$median_abs_loading_across_models,
  probs = SMALL_QUANTILE,
  na.rm = TRUE,
  names = FALSE
)

small_panel <- core_vip3[
  core_vip3$median_abs_loading_across_models >=
    small_threshold,
  ,
  drop = FALSE
]

small_panel <- small_panel[
  order(
    -small_panel$median_abs_loading_across_models,
    -small_panel$VIP_median_across_present,
    match(
      small_panel$view_code,
      VIEW_LEVELS
    ),
    small_panel$feature_label
  ),
  ,
  drop = FALSE
]

cat("\n============================================================\n")
cat("WORKING SMALL PANEL\n")
cat(
  "Rule: CORE + VIP>=1 in all 3 + upper ",
  sprintf("%.0f%%", 100 * (1 - SMALL_QUANTILE)),
  " of median |loading| across models/components\n",
  sep = ""
)
cat("Loading threshold:", sprintf("%.6f", small_threshold), "\n")
cat("n =", nrow(small_panel), "\n")

view_small <- table(
  factor(
    toupper(small_panel$view_code),
    levels = toupper(VIEW_LEVELS)
  )
)

cat(
  "Views: TX=", view_small["TX"],
  " PR=", view_small["PR"],
  " ME=", view_small["ME"],
  " CL=", view_small["CL"],
  "\n",
  sep = ""
)

cat("\nFEATURES\n")

print(
  small_panel[
    ,
    c(
      "view_code",
      "feature_label",
      "feature_original",
      "VIP_median_across_present",
      "median_abs_loading_across_models",
      "n_SGL_selected"
    )
  ],
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# SHARED_2 + SPECIFIC context
# ------------------------------------------------------------------------------

shared2 <- master[
  master$consensus_tier == "SHARED_2",
  ,
  drop = FALSE
]

specific <- master[
  master$consensus_tier == "SPECIFIC",
  ,
  drop = FALSE
]

cat("\n============================================================\n")
cat("SHARED_2 CONTEXT\n")
cat("n =", nrow(shared2), "\n")

print(
  shared2[
    order(
      -shared2$n_VIP_ge_threshold,
      -shared2$VIP_median_across_present
    ),
    c(
      "families",
      "view_code",
      "feature_label",
      "n_VIP_ge_threshold",
      "VIP_median_across_present",
      "median_abs_loading_across_models",
      "n_SGL_selected"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nSPECIFIC COUNTS BY ANALYSIS\n")

for (fam in c("P1", "P2", "P3")) {
  z <- specific[
    specific$families == fam,
    ,
    drop = FALSE
  ]

  cat(
    fam,
    ": n=",
    nrow(z),
    " | ",
    paste(
      names(
        table(
          toupper(z$view_code)
        )
      ),
      as.integer(
        table(
          toupper(z$view_code)
        )
      ),
      sep = "=",
      collapse = " "
    ),
    "\n",
    sep = ""
  )
}

# ------------------------------------------------------------------------------
# Load P1 observed data and print PCA summaries
# ------------------------------------------------------------------------------

cons_manifest <- read_csv_checked(
  file.path(
    CONSENSUS,
    "RUN_MANIFEST.csv"
  ),
  required = c("item", "value")
)

panel17 <- as.character(
  cons_manifest$value[
    match(
      "PANEL17_RUN",
      cons_manifest$item
    )
  ]
)

assert(
  length(panel17) == 1L &&
    !is.na(panel17) &&
    file.exists(
      file.path(
        panel17,
        "RUN_MANIFEST.csv"
      )
    ),
  "Could not resolve PANEL17_RUN from consensus manifest."
)

p17_manifest <- read_csv_checked(
  file.path(
    panel17,
    "RUN_MANIFEST.csv"
  ),
  required = c("item", "value")
)

panel10 <- as.character(
  p17_manifest$value[
    match(
      "PANEL10_RUN",
      p17_manifest$item
    )
  ]
)

assert(
  length(panel10) == 1L &&
    !is.na(panel10),
  "Could not resolve PANEL10_RUN."
)

input_audit <- read_csv_checked(
  file.path(
    panel10,
    "tables",
    "00_INPUT_AUDIT.csv"
  ),
  required = c(
    "family",
    "ready_for_modeling"
  )
)

p1_row <- input_audit[
  input_audit$family == "P1",
  ,
  drop = FALSE
]

assert(
  nrow(p1_row) == 1L,
  "Expected exactly one P1 input row."
)

rfm_path <- as.character(
  p1_row$ready_for_modeling[1L]
)

assert(
  file.exists(rfm_path),
  paste("Missing P1 ready_for_modeling:", rfm_path)
)

rfm <- readRDS(rfm_path)

train_meta <- extract_group_table(
  rfm$grupo_train
)

test_meta <- extract_group_table(
  rfm$grupo_test
)

assert(
  !anyNA(train_meta$grupo),
  "Bad TRAIN phenotype labels."
)

assert(
  !anyNA(test_meta$grupo),
  "Bad TEST phenotype labels."
)

cat("\n============================================================\n")
cat("PCA DESCRIPTIVE AUDIT\n")
cat(
  "PCA is fitted on TRAIN only; TEST is projected. ",
  "ALL-26 statistics are descriptive only.\n",
  sep = ""
)
cat("============================================================\n")

# Recommended primary PCA
print_pca_summary(
  "SMALL_CORE_TRIAL",
  small_panel,
  rfm,
  train_meta,
  test_meta
)

# Reference / sensitivity PCA
core_all <- master[
  master$consensus_tier == "CORE",
  ,
  drop = FALSE
]

print_pca_summary(
  "CORE_ALL_3OF3",
  core_all,
  rfm,
  train_meta,
  test_meta
)

print_pca_summary(
  "SHARED_2_ONLY",
  shared2,
  rfm,
  train_meta,
  test_meta
)

specific_p1 <- specific[
  specific$families == "P1",
  ,
  drop = FALSE
]

print_pca_summary(
  "SPECIFIC_P1_ONLY",
  specific_p1,
  rfm,
  train_meta,
  test_meta
)

specific_p2 <- specific[
  specific$families == "P2",
  ,
  drop = FALSE
]

print_pca_summary(
  "SPECIFIC_P2_ONLY",
  specific_p2,
  rfm,
  train_meta,
  test_meta
)

specific_p3 <- specific[
  specific$families == "P3",
  ,
  drop = FALSE
]

print_pca_summary(
  "SPECIFIC_P3_ONLY",
  specific_p3,
  rfm,
  train_meta,
  test_meta
)

cat("\n============================================================\n")
cat("RECOMMENDED NEXT STEP\n")
cat(
  "Primary figure: PCA PC1-PC2 of SMALL_CORE_TRIAL, TRAIN-fit + TEST projected.\n",
  "Secondary sensitivity: CORE_ALL and SHARED_2_ONLY.\n",
  "SPECIFIC_P1 can be supplementary; P2/P3 specific PCA should be skipped if <2 features.\n",
  "Do not use these PCA statistics to re-select the panel.\n",
  sep = ""
)
cat("============================================================\n")

cat(
  "\nPLOTS OPENED ON SCREEN.\n",
  "No PNG/PDF files were written.\n",
  "Press [ENTER] when you are finished inspecting the PCA windows to close the script...\n",
  sep = ""
)

invisible(
  readLines(
    con = file("stdin"),
    n = 1L
  )
)


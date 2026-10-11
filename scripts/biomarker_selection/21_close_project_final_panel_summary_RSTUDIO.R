# ==============================================================================
# METATEST — 20_probe_small_final_panel_RSTUDIO.R
#
# USO:
#   Ejecutar desde la CONSOLA de RStudio, NO desde la pestaña Terminal:
#
#   source("20_probe_small_final_panel_RSTUDIO.R")
#
# NO guarda archivos. Todos los gráficos se imprimen en el panel "Plots" de
# RStudio y quedan disponibles en el historial de gráficos.
#
# Selección exploratoria del panel pequeño:
#   CORE (3/3) + VIP>=1 en P1/P2/P3 + Q75 de median(|loading|)
#
# SGL se mantiene como anotación, NO como filtro.
# PCA:
#   ajuste en TRAIN;
#   TEST proyectado con centro/escala/rotación de TRAIN.
# ==============================================================================

options(stringsAsFactors = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Instala ggplot2 antes de ejecutar este script.", call. = FALSE)
}

UPSTREAM_DIR <- Sys.getenv(
  "UPSTREAM_DIR",
  unset = path.expand(
    "~/Documents/PhD/metatest/results_upstream_local_20260908/comparacion_modelamientos/upstream_comun"
  )
)

FINAL_CONSENSUS_RUN <- Sys.getenv("FINAL_CONSENSUS_RUN", unset = "")
SMALL_QUANTILE <- as.numeric(
  Sys.getenv("PANEL20_SMALL_QUANTILE", unset = "0.75")
)

CLASS_LEVELS <- c("NP", "FD", "SO")
VIEW_LEVELS <- c("tx", "pr", "me", "cl")
TOP_LOADING_LABELS <- as.integer(
  Sys.getenv("PANEL20_TOP_LOADING_LABELS", unset = "27")
)

assert <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
}

canonical <- function(path, must_work = TRUE) {
  normalizePath(path.expand(path), winslash = "/", mustWork = must_work)
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

  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
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

    if (nzchar(p) && dir.exists(p)) {
      return(canonical(p))
    }
  }

  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)

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

  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)

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

build_observed_matrix <- function(rfm, features, train_meta, test_meta) {
  train_cols <- list()
  test_cols <- list()
  metadata <- list()
  rr <- 1L

  for (i in seq_len(nrow(features))) {
    view <- tolower(as.character(features$view_code[i]))
    feature <- as.character(features$feature_original[i])

    if (
      !view %in% names(rfm) ||
      is.null(rfm[[view]]$train) ||
      is.null(rfm[[view]]$test)
    ) {
      next
    }

    tr <- as.matrix(rfm[[view]]$train)
    te <- as.matrix(rfm[[view]]$test)

    if (!feature %in% rownames(tr)) {
      next
    }

    a <- as.numeric(
      tr[
        feature,
        train_meta$sample_id
      ]
    )

    if (
      any(!is.finite(a)) ||
      !is.finite(stats::sd(a)) ||
      stats::sd(a) <= 0
    ) {
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

    metadata[[rr]] <- data.frame(
      machine_key = machine_key,
      view_code = view,
      feature_original = feature,
      feature_label = as.character(features$feature_label[i]),
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  assert(length(train_cols) >= 2L, "Fewer than 2 usable features.")

  xtr <- do.call(cbind, train_cols)
  xte <- do.call(cbind, test_cols)

  rownames(xtr) <- train_meta$sample_id
  rownames(xte) <- test_meta$sample_id

  list(
    train = xtr,
    test = xte,
    metadata = do.call(rbind, metadata)
  )
}

# ==============================================================================
# Locate latest consensus
# ==============================================================================

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
    "n_SGL_selected"
  )
)

master$view_code <- tolower(as.character(master$view_code))
master$consensus_tier <- as.character(master$consensus_tier)

core_vip3 <- master[
  master$consensus_tier == "CORE" &
    master$n_VIP_ge_threshold == 3L,
  ,
  drop = FALSE
]

assert(
  nrow(core_vip3) >= 2L,
  "No CORE+VIP3 pool available."
)

small_threshold <- stats::quantile(
  core_vip3$median_abs_loading_across_models,
  probs = SMALL_QUANTILE,
  na.rm = TRUE,
  names = FALSE
)

small_panel <- core_vip3[
  core_vip3$median_abs_loading_across_models >= small_threshold,
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

cat("============================================================\n")
cat("SMALL FINAL PANEL — RStudio probe\n")
cat("Consensus run:", CONSENSUS, "\n")
cat("Rule: CORE + VIP>=1 in all 3 + Q75 median |loading|\n")
cat("n =", nrow(small_panel), "\n")

views <- table(
  factor(
    toupper(small_panel$view_code),
    levels = toupper(VIEW_LEVELS)
  )
)

cat(
  "TX=", views["TX"],
  " PR=", views["PR"],
  " ME=", views["ME"],
  " CL=", views["CL"],
  "\n",
  sep = ""
)

cat("\nFINAL MARKERS\n")

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

# ==============================================================================
# Resolve observed P1 matrix
# ==============================================================================

cons_manifest <- read_csv_checked(
  file.path(
    CONSENSUS,
    "RUN_MANIFEST.csv"
  ),
  required = c(
    "item",
    "value"
  )
)

panel17 <- as.character(
  cons_manifest$value[
    match(
      "PANEL17_RUN",
      cons_manifest$item
    )
  ]
)

p17_manifest <- read_csv_checked(
  file.path(
    panel17,
    "RUN_MANIFEST.csv"
  ),
  required = c(
    "item",
    "value"
  )
)

panel10 <- as.character(
  p17_manifest$value[
    match(
      "PANEL10_RUN",
      p17_manifest$item
    )
  ]
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

p1 <- input_audit[
  input_audit$family == "P1",
  ,
  drop = FALSE
]

assert(nrow(p1) == 1L, "Expected one P1 input row.")

rfm_path <- as.character(
  p1$ready_for_modeling[1L]
)

assert(
  file.exists(rfm_path),
  paste("Missing ready_for_modeling:", rfm_path)
)

rfm <- readRDS(rfm_path)

train_meta <- extract_group_table(
  rfm$grupo_train
)

test_meta <- extract_group_table(
  rfm$grupo_test
)

mat <- build_observed_matrix(
  rfm = rfm,
  features = small_panel,
  train_meta = train_meta,
  test_meta = test_meta
)

cat(
  "\nUsable final markers in observed matrix:",
  ncol(mat$train),
  "\n"
)

# ==============================================================================
# PCA TRAIN-fit; TEST projection
# ==============================================================================

pca <- stats::prcomp(
  mat$train,
  center = TRUE,
  scale. = TRUE
)

tr_scores <- as.data.frame(
  pca$x
)

te_scores <- as.data.frame(
  stats::predict(
    pca,
    newdata = mat$test
  )
)

tr_scores$sample_id <- rownames(mat$train)
te_scores$sample_id <- rownames(mat$test)

tr_scores$grupo <- train_meta$grupo[
  match(
    tr_scores$sample_id,
    train_meta$sample_id
  )
]

te_scores$grupo <- test_meta$grupo[
  match(
    te_scores$sample_id,
    test_meta$sample_id
  )
]

tr_scores$split <- "TRAIN"
te_scores$split <- "TEST_EXPLORATORY_REUSED"

scores <- rbind(
  tr_scores,
  te_scores
)

var <- pca$sdev^2
var <- var / sum(var)

eta_pc1 <- eta2_oneway(
  scores$PC1,
  scores$grupo
)

eta_pc2 <- eta2_oneway(
  scores$PC2,
  scores$grupo
)

sil <- mean_silhouette_2d(
  scores[, c("PC1", "PC2"), drop = FALSE],
  scores$grupo
)

cat(
  "PCA variance: PC1=",
  sprintf("%.1f%%", 100 * var[1L]),
  " | PC2=",
  sprintf("%.1f%%", 100 * var[2L]),
  "\n",
  sep = ""
)

cat(
  "ALL-26 descriptive eta2: PC1=",
  sprintf("%.3f", eta_pc1),
  " | PC2=",
  sprintf("%.3f", eta_pc2),
  " | silhouette=",
  sprintf("%.3f", sil),
  "\n",
  sep = ""
)

# ==============================================================================
# PCA loadings for FINAL MARKERS
# ==============================================================================

loading <- as.data.frame(
  pca$rotation[
    ,
    c("PC1", "PC2"),
    drop = FALSE
  ]
)

loading$machine_key <- rownames(loading)

loading <- merge(
  loading,
  mat$metadata,
  by = "machine_key",
  all.x = TRUE,
  sort = FALSE
)

loading$loading_strength <- sqrt(
  loading$PC1^2 +
    loading$PC2^2
)

loading <- loading[
  order(
    -loading$loading_strength,
    loading$feature_label
  ),
  ,
  drop = FALSE
]

loading$rank <- seq_len(
  nrow(loading)
)

loading$display_label <- paste0(
  toupper(loading$view_code),
  " | ",
  loading$feature_label
)

cat("\nFINAL MARKER PCA LOADINGS\n")

print(
  loading[
    ,
    c(
      "rank",
      "view_code",
      "feature_label",
      "PC1",
      "PC2",
      "loading_strength"
    )
  ],
  row.names = FALSE
)

# ==============================================================================
# PLOT 1 — PCA scores
# ==============================================================================

p_scores <- ggplot2::ggplot(
  scores,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    colour = grupo,
    shape = split
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 3,
    linewidth = 0.3
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = 3,
    linewidth = 0.3
  ) +
  ggplot2::geom_point(
    size = 3.4,
    alpha = 0.90
  ) +
  ggplot2::labs(
    title =
      "Small final panel — PCA scores",
    subtitle =
      paste0(
        "TRAIN-fit; TEST projected | eta2 PC1=",
        sprintf("%.3f", eta_pc1),
        ", PC2=",
        sprintf("%.3f", eta_pc2),
        " | silhouette=",
        sprintf("%.3f", sil)
      ),
    x =
      paste0(
        "PC1 (",
        sprintf("%.1f%%", 100 * var[1L]),
        ")"
      ),
    y =
      paste0(
        "PC2 (",
        sprintf("%.1f%%", 100 * var[2L]),
        ")"
      ),
    colour = "Phenotype",
    shape = "Split"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank(),
    plot.title =
      ggplot2::element_text(face = "bold"),
    legend.position =
      "bottom"
  )

print(p_scores)

# ==============================================================================
# PLOT 2 — FINAL MARKERS loading map
# ==============================================================================

label_n <- min(
  TOP_LOADING_LABELS,
  nrow(loading)
)

loading_plot <- loading[
  seq_len(label_n),
  ,
  drop = FALSE
]

p_load <- ggplot2::ggplot(
  loading,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    colour = view_code
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 3,
    linewidth = 0.3
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = 3,
    linewidth = 0.3
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      y = 0,
      xend = PC1,
      yend = PC2,
      colour = view_code
    ),
    arrow = grid::arrow(
      length = grid::unit(0.11, "inches")
    ),
    linewidth = 0.45,
    alpha = 0.75
  ) +
  ggplot2::geom_point(
    size = 2.2
  )

if (
  requireNamespace(
    "ggrepel",
    quietly = TRUE
  )
) {
  p_load <- p_load +
    ggrepel::geom_text_repel(
      data = loading_plot,
      ggplot2::aes(
        label = display_label
      ),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.15,
      segment.alpha = 0.4
    )
} else {
  p_load <- p_load +
    ggplot2::geom_text(
      data = loading_plot,
      ggplot2::aes(
        label = display_label
      ),
      size = 3,
      check_overlap = TRUE,
      hjust = -0.05,
      vjust = 0.4
    )
}

p_load <- p_load +
  ggplot2::labs(
    title =
      "Small final panel — PCA loading map",
    subtitle =
      paste0(
        nrow(loading),
        " final markers; labels ordered by sqrt(PC1^2 + PC2^2)."
      ),
    x = "PC1 loading",
    y = "PC2 loading",
    colour = "Omic"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank(),
    plot.title =
      ggplot2::element_text(face = "bold"),
    legend.position =
      "bottom"
  )

print(p_load)

# ==============================================================================
# PLOT 3 — classical-style biplot: samples + final marker arrows
# ==============================================================================

# Scale loading arrows only for visualization so they occupy the score plane.
score_radius <- max(
  sqrt(
    scores$PC1^2 +
      scores$PC2^2
  ),
  na.rm = TRUE
)

loading_radius <- max(
  sqrt(
    loading$PC1^2 +
      loading$PC2^2
  ),
  na.rm = TRUE
)

arrow_multiplier <- if (
  is.finite(score_radius) &&
  is.finite(loading_radius) &&
  loading_radius > 0
) {
  0.72 *
    score_radius /
    loading_radius
} else {
  1
}

loading$biplot_x <- loading$PC1 *
  arrow_multiplier

loading$biplot_y <- loading$PC2 *
  arrow_multiplier

loading_biplot_labels <- loading[
  seq_len(label_n),
  ,
  drop = FALSE
]

p_biplot <- ggplot2::ggplot() +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 3,
    linewidth = 0.3
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = 3,
    linewidth = 0.3
  ) +
  ggplot2::geom_point(
    data = scores,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      colour = grupo,
      shape = split
    ),
    size = 3.2,
    alpha = 0.88
  ) +
  ggplot2::geom_segment(
    data = loading,
    ggplot2::aes(
      x = 0,
      y = 0,
      xend = biplot_x,
      yend = biplot_y
    ),
    arrow = grid::arrow(
      length = grid::unit(0.10, "inches")
    ),
    linewidth = 0.35,
    alpha = 0.45
  )

if (
  requireNamespace(
    "ggrepel",
    quietly = TRUE
  )
) {
  p_biplot <- p_biplot +
    ggrepel::geom_text_repel(
      data = loading_biplot_labels,
      ggplot2::aes(
        x = biplot_x,
        y = biplot_y,
        label = display_label
      ),
      size = 2.8,
      max.overlaps = Inf,
      box.padding = 0.20,
      point.padding = 0.10,
      segment.alpha = 0.3
    )
} else {
  p_biplot <- p_biplot +
    ggplot2::geom_text(
      data = loading_biplot_labels,
      ggplot2::aes(
        x = biplot_x,
        y = biplot_y,
        label = display_label
      ),
      size = 2.8,
      check_overlap = TRUE
    )
}

p_biplot <- p_biplot +
  ggplot2::labs(
    title =
      "Small final panel — PCA biplot",
    subtitle =
      "Samples: PCA scores. Arrows: final biomarkers (loadings; visual scaling only).",
    x =
      paste0(
        "PC1 (",
        sprintf("%.1f%%", 100 * var[1L]),
        ")"
      ),
    y =
      paste0(
        "PC2 (",
        sprintf("%.1f%%", 100 * var[2L]),
        ")"
      ),
    colour = "Phenotype",
    shape = "Split"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank(),
    plot.title =
      ggplot2::element_text(face = "bold"),
    legend.position =
      "bottom"
  )

print(p_biplot)

# ==============================================================================
# PLOT 4/5 — signed contribution bars for PCA PC1 / PC2
# ==============================================================================

make_signed_bar <- function(component) {
  z <- loading

  z$loading_value <- z[[component]]
  z$abs_loading_value <- abs(z$loading_value)

  z <- z[
    order(
      -z$abs_loading_value,
      z$feature_label
    ),
    ,
    drop = FALSE
  ]

  z$direction <- ifelse(
    z$loading_value >= 0,
    "Positive",
    "Negative"
  )

  z <- z[
    order(z$loading_value),
    ,
    drop = FALSE
  ]

  z$display_label <- factor(
    z$display_label,
    levels = z$display_label
  )

  ggplot2::ggplot(
    z,
    ggplot2::aes(
      x = display_label,
      y = loading_value,
      fill = direction
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = 3
    ) +
    ggplot2::labs(
      title =
        paste0(
          "Small final panel — ",
          component,
          " marker loadings"
        ),
      subtitle =
        "All final biomarkers; sign = PCA direction, magnitude = contribution.",
      x = NULL,
      y =
        paste0(
          component,
          " loading"
        ),
      fill = "Direction"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor =
        ggplot2::element_blank(),
      plot.title =
        ggplot2::element_text(face = "bold"),
      legend.position =
        "bottom"
    )
}

print(
  make_signed_bar("PC1")
)

print(
  make_signed_bar("PC2")
)

cat("\n============================================================\n")
cat("Five plots were sent to the RStudio Plots pane:\n")
cat("  1. PCA scores (samples)\n")
cat("  2. PCA loading map (final biomarkers)\n")
cat("  3. PCA biplot (samples + final biomarkers)\n")
cat("  4. PC1 signed loading bars\n")
cat("  5. PC2 signed loading bars\n")
cat("\nUse the Back/Forward arrows in the RStudio Plots pane to navigate.\n")
cat("NO plot files were written.\n")
cat("============================================================\n")


# ==============================================================================
# 6. FINAL CLOSING SUMMARY — NO RECALCULATION, NO FILES WRITTEN
# ==============================================================================
#
# This section only joins already-frozen outputs:
#   * strict Evidence-A direction from framework
#   * existing PLS-DA scores/loadings
#   * frozen Q5 pathway-feature catalog
#
# It does NOT refit PLS-DA, rerun enrichment, redefine Evidence A, or use TEST
# to select biomarkers.
# ==============================================================================

cat("\n\n============================================================\n")
cat("FINAL CLOSING SUMMARY — FROZEN OUTPUTS ONLY\n")
cat("============================================================\n")

# ------------------------------------------------------------------------------
# 6A. Read frozen Evidence-A annotations from script 17
# ------------------------------------------------------------------------------

evidence_all <- read_csv_checked(
  file.path(
    panel17,
    "tables",
    "02_EVIDENCE_A_FEATURES_FOR_PLSDA_ALL_FAMILIES.csv"
  ),
  required = c(
    "family",
    "view_code",
    "feature_original",
    "feature_label",
    "framework_rec_bidirectional_pattern",
    "framework_obs_bidirectional_pattern",
    "framework_rec_evidence_grade",
    "framework_obs_evidence_grade",
    "framework_direction_concordance_prop"
  )
)

evidence_all$family <- as.character(evidence_all$family)
evidence_all$view_code <- tolower(as.character(evidence_all$view_code))
evidence_all$feature_original <- as.character(evidence_all$feature_original)
evidence_all$feature_label <- as.character(evidence_all$feature_label)

evidence_all$feature_key_final <- paste(
  evidence_all$view_code,
  evidence_all$feature_original,
  sep = "\r"
)

small_panel$feature_key_final <- paste(
  tolower(as.character(small_panel$view_code)),
  as.character(small_panel$feature_original),
  sep = "\r"
)

master$feature_key_final <- paste(
  tolower(as.character(master$view_code)),
  as.character(master$feature_original),
  sep = "\r"
)

pattern_to_class <- function(x) {
  out <- rep(NA_character_, length(x))

  out[x == "NP_like_high"] <- "NP"
  out[x == "FD_like_high"] <- "FD"
  out[x == "SO_like_high"] <- "SO"

  out
}

# ------------------------------------------------------------------------------
# 6B. WHY ARE P1 CLINICAL FEATURES NOT IN THE SMALL FINAL PANEL?
# ------------------------------------------------------------------------------

p1_clin <- evidence_all[
  evidence_all$family == "P1" &
    evidence_all$view_code == "cl",
  ,
  drop = FALSE
]

if (nrow(p1_clin)) {
  clinical_master <- master[
    master$view_code == "cl",
    ,
    drop = FALSE
  ]

  keep_master_cols <- c(
    "feature_key_final",
    "n_families",
    "families",
    "consensus_tier",
    "VIP_P1",
    "VIP_P2",
    "VIP_P3",
    "n_VIP_ge_threshold",
    "median_abs_loading_across_models",
    "n_SGL_selected"
  )

  keep_master_cols <- keep_master_cols[
    keep_master_cols %in% names(clinical_master)
  ]

  clinical_audit <- merge(
    p1_clin,
    clinical_master[
      ,
      keep_master_cols,
      drop = FALSE
    ],
    by = "feature_key_final",
    all.x = TRUE,
    sort = FALSE
  )

  clinical_audit$in_small_final_panel <-
    clinical_audit$feature_key_final %in%
      small_panel$feature_key_final

  clinical_audit$small_panel_reason <- vapply(
    seq_len(nrow(clinical_audit)),
    function(i) {
      if (isTRUE(clinical_audit$in_small_final_panel[i])) {
        return("IN_SMALL_PANEL")
      }

      nf <- suppressWarnings(
        as.integer(
          clinical_audit$n_families[i]
        )
      )

      nv <- suppressWarnings(
        as.integer(
          clinical_audit$n_VIP_ge_threshold[i]
        )
      )

      ml <- suppressWarnings(
        as.numeric(
          clinical_audit$median_abs_loading_across_models[i]
        )
      )

      if (!is.finite(nf) || nf < 3L) {
        return(
          paste0(
            "NOT_CORE_3OF3 (present in ",
            ifelse(is.finite(nf), nf, 0L),
            "/3)"
          )
        )
      }

      if (!is.finite(nv) || nv < 3L) {
        return(
          paste0(
            "CORE_BUT_VIP_NOT_3OF3 (",
            ifelse(is.finite(nv), nv, 0L),
            "/3)"
          )
        )
      }

      if (
        !is.finite(ml) ||
        ml < small_threshold
      ) {
        return(
          paste0(
            "CORE+VIP3_BUT_BELOW_Q75_LOADING (",
            sprintf("%.4f", ml),
            " < ",
            sprintf("%.4f", small_threshold),
            ")"
          )
        )
      }

      "OTHER"
    },
    character(1)
  )

  clinical_audit <- clinical_audit[
    order(
      -clinical_audit$n_families,
      -clinical_audit$VIP_P1,
      clinical_audit$feature_label
    ),
    ,
    drop = FALSE
  ]

  cat("\nP1 CLINICAL STRICT-A FEATURES — WHY THEY DO/DO NOT ENTER SMALL PANEL\n")

  clinical_print_cols <- c(
    "feature_label",
    "framework_rec_bidirectional_pattern",
    "framework_obs_bidirectional_pattern",
    "n_families",
    "families",
    "VIP_P1",
    "VIP_P2",
    "VIP_P3",
    "median_abs_loading_across_models",
    "n_SGL_selected",
    "small_panel_reason"
  )

  clinical_print_cols <- clinical_print_cols[
    clinical_print_cols %in% names(clinical_audit)
  ]

  print(
    clinical_audit[
      ,
      clinical_print_cols,
      drop = FALSE
    ],
    row.names = FALSE
  )
} else {
  cat("\nNo P1 clinical strict-A features were found.\n")
}

# ------------------------------------------------------------------------------
# 6C. FRAMEWORK DIRECTION FOR EACH FINAL SMALL-PANEL FEATURE
# ------------------------------------------------------------------------------

ev_small <- evidence_all[
  evidence_all$feature_key_final %in%
    small_panel$feature_key_final,
  ,
  drop = FALSE
]

ev_split <- split(
  ev_small,
  ev_small$feature_key_final
)

direction_summary <- do.call(
  rbind,
  lapply(
    ev_split,
    function(z) {
      rec_patterns <- sort(
        unique(
          as.character(
            z$framework_rec_bidirectional_pattern
          )
        )
      )

      obs_patterns <- sort(
        unique(
          as.character(
            z$framework_obs_bidirectional_pattern
          )
        )
      )

      all_patterns <- sort(
        unique(
          c(
            rec_patterns,
            obs_patterns
          )
        )
      )

      pattern_consensus <- if (
        length(all_patterns) == 1L
      ) {
        all_patterns[1L]
      } else {
        paste(
          all_patterns,
          collapse = ";"
        )
      }

      direction_class <- if (
        length(all_patterns) == 1L
      ) {
        pattern_to_class(
          all_patterns[1L]
        )
      } else {
        NA_character_
      }

      data.frame(
        feature_key_final =
          z$feature_key_final[1L],
        view_code =
          z$view_code[1L],
        feature_original =
          z$feature_original[1L],
        feature_label =
          z$feature_label[1L],
        families_framework =
          paste(
            sort(
              unique(
                z$family
              )
            ),
            collapse = ";"
          ),
        framework_pattern =
          pattern_consensus,
        framework_direction_class =
          direction_class,
        rec_obs_same_direction_all =
          length(all_patterns) == 1L,
        mean_direction_concordance_prop =
          mean(
            suppressWarnings(
              as.numeric(
                z$framework_direction_concordance_prop
              )
            ),
            na.rm = TRUE
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(direction_summary) <- NULL

# ------------------------------------------------------------------------------
# 6D. PLS-DA DIRECTIONAL CORROBORATION FROM EXISTING SCORES + LOADINGS
# ------------------------------------------------------------------------------

model_status <- read_csv_checked(
  file.path(
    panel17,
    "tables",
    "14_MODEL_STATUS_ALL_FAMILIES.csv"
  ),
  required = c(
    "family",
    "chosen_ncomp"
  )
)

pls_direction_rows <- list()
rr_pls <- 1L

for (fam in CLASS_LEVELS[0]) {
  # never entered; placeholder avoids accidental partial-evaluation issues
}

for (fam in c("P1", "P2", "P3")) {
  chosen <- suppressWarnings(
    as.integer(
      model_status$chosen_ncomp[
        match(
          fam,
          model_status$family
        )
      ]
    )
  )

  if (!is.finite(chosen) || chosen < 1L) {
    next
  }

  score_path <- file.path(
    panel17,
    "tables",
    fam,
    paste0(
      fam,
      "_08_PLSDA_SCORES.csv"
    )
  )

  load_path <- file.path(
    panel17,
    "tables",
    fam,
    paste0(
      fam,
      "_07_PLSDA_LOADINGS.csv"
    )
  )

  scores_f <- read_csv_checked(
    score_path,
    required = c(
      "sample_id",
      "grupo",
      "split"
    )
  )

  loads_f <- read_csv_checked(
    load_path,
    required = c(
      "component",
      "loading",
      "view_code",
      "feature_original",
      "feature_label"
    )
  )

  score_cols <- paste0(
    "PLS",
    seq_len(chosen)
  )

  assert(
    all(score_cols %in% names(scores_f)),
    paste(
      fam,
      "score table lacks selected PLS components."
    )
  )

  train_scores_f <- scores_f[
    scores_f$split == "TRAIN",
    ,
    drop = FALSE
  ]

  class_means <- do.call(
    rbind,
    lapply(
      CLASS_LEVELS,
      function(cls) {
        z <- train_scores_f[
          train_scores_f$grupo == cls,
          score_cols,
          drop = FALSE
        ]

        out <- colMeans(
          z,
          na.rm = TRUE
        )

        data.frame(
          grupo = cls,
          t(out),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  loads_f$view_code <- tolower(
    as.character(
      loads_f$view_code
    )
  )

  loads_f$feature_original <- as.character(
    loads_f$feature_original
  )

  loads_f$feature_key_final <- paste(
    loads_f$view_code,
    loads_f$feature_original,
    sep = "\r"
  )

  loads_f <- loads_f[
    loads_f$feature_key_final %in%
      small_panel$feature_key_final &
      loads_f$component <= chosen,
    ,
    drop = FALSE
  ]

  for (fk in unique(loads_f$feature_key_final)) {
    z <- loads_f[
      loads_f$feature_key_final == fk,
      ,
      drop = FALSE
    ]

    loading_vec <- rep(
      0,
      chosen
    )

    for (cc in seq_len(chosen)) {
      hit <- z$loading[
        z$component == cc
      ]

      if (length(hit)) {
        loading_vec[cc] <- suppressWarnings(
          as.numeric(
            hit[1L]
          )
        )
      }
    }

    align <- vapply(
      CLASS_LEVELS,
      function(cls) {
        means_cls <- as.numeric(
          class_means[
            class_means$grupo == cls,
            score_cols,
            drop = FALSE
          ]
        )

        sum(
          loading_vec *
            means_cls
        )
      },
      numeric(1)
    )

    ord <- order(
      align,
      decreasing = TRUE
    )

    pls_high_class <- names(align)[
      ord[1L]
    ]

    margin <- if (
      length(ord) >= 2L
    ) {
      align[
        ord[1L]
      ] -
        align[
          ord[2L]
        ]
    } else {
      NA_real_
    }

    pls_direction_rows[[rr_pls]] <- data.frame(
      feature_key_final = fk,
      family = fam,
      PLS_high_class =
        pls_high_class,
      PLS_alignment_NP =
        unname(
          align["NP"]
        ),
      PLS_alignment_FD =
        unname(
          align["FD"]
        ),
      PLS_alignment_SO =
        unname(
          align["SO"]
        ),
      PLS_alignment_margin =
        unname(
          margin
        ),
      stringsAsFactors = FALSE
    )

    rr_pls <- rr_pls + 1L
  }
}

pls_direction_long <- if (
  length(pls_direction_rows)
) {
  do.call(
    rbind,
    pls_direction_rows
  )
} else {
  data.frame()
}

if (nrow(pls_direction_long)) {
  pls_direction_long <- merge(
    pls_direction_long,
    direction_summary[
      ,
      c(
        "feature_key_final",
        "framework_direction_class"
      ),
      drop = FALSE
    ],
    by = "feature_key_final",
    all.x = TRUE,
    sort = FALSE
  )

  pls_direction_long$PLS_framework_concordant <-
    !is.na(
      pls_direction_long$framework_direction_class
    ) &
    pls_direction_long$PLS_high_class ==
      pls_direction_long$framework_direction_class

  pls_split <- split(
    pls_direction_long,
    pls_direction_long$feature_key_final
  )

  pls_consensus <- do.call(
    rbind,
    lapply(
      pls_split,
      function(z) {
        data.frame(
          feature_key_final =
            z$feature_key_final[1L],
          PLS_high_classes =
            paste(
              paste0(
                z$family,
                ":",
                z$PLS_high_class
              ),
              collapse = ";"
            ),
          PLS_framework_concordant_n =
            sum(
              z$PLS_framework_concordant,
              na.rm = TRUE
            ),
          PLS_framework_evaluable_n =
            sum(
              !is.na(
                z$framework_direction_class
              )
            ),
          min_PLS_alignment_margin =
            min(
              z$PLS_alignment_margin,
              na.rm = TRUE
            ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  direction_summary <- merge(
    direction_summary,
    pls_consensus,
    by = "feature_key_final",
    all.x = TRUE,
    sort = FALSE
  )
}

direction_summary <- merge(
  direction_summary,
  small_panel[
    ,
    c(
      "feature_key_final",
      "VIP_median_across_present",
      "median_abs_loading_across_models",
      "n_SGL_selected"
    ),
    drop = FALSE
  ],
  by = "feature_key_final",
  all.x = TRUE,
  sort = FALSE
)

direction_summary <- direction_summary[
  order(
    match(
      direction_summary$framework_direction_class,
      CLASS_LEVELS
    ),
    match(
      direction_summary$view_code,
      VIEW_LEVELS
    ),
    -direction_summary$median_abs_loading_across_models,
    direction_summary$feature_label
  ),
  ,
  drop = FALSE
]

cat("\nFINAL SMALL PANEL — FRAMEWORK DIRECTION + PLS-DA CORROBORATION\n")

direction_print_cols <- c(
  "framework_direction_class",
  "framework_pattern",
  "view_code",
  "feature_label",
  "PLS_high_classes",
  "PLS_framework_concordant_n",
  "PLS_framework_evaluable_n",
  "min_PLS_alignment_margin",
  "VIP_median_across_present",
  "median_abs_loading_across_models",
  "n_SGL_selected"
)

direction_print_cols <- direction_print_cols[
  direction_print_cols %in%
    names(
      direction_summary
    )
]

print(
  direction_summary[
    ,
    direction_print_cols,
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nDIRECTION COUNTS\n")

print(
  table(
    Direction =
      direction_summary$framework_direction_class,
    Omic =
      toupper(
        direction_summary$view_code
      )
  )
)

if (
  "PLS_framework_concordant_n" %in%
    names(direction_summary)
) {
  cat(
    "\nFeatures concordant framework <-> PLS-DA in all evaluable P1/P2/P3 models: ",
    sum(
      direction_summary$PLS_framework_concordant_n ==
        direction_summary$PLS_framework_evaluable_n,
      na.rm = TRUE
    ),
    "/",
    sum(
      is.finite(
        direction_summary$PLS_framework_evaluable_n
      )
    ),
    "\n",
    sep = ""
  )
}

# ------------------------------------------------------------------------------
# 6E. FROZEN PATHWAY LINKS MATCHED TO EACH FEATURE'S PHENOTYPIC DIRECTION
# ------------------------------------------------------------------------------

panel10_manifest <- read_csv_checked(
  file.path(
    panel10,
    "RUN_MANIFEST.csv"
  ),
  required = c(
    "item",
    "value"
  )
)

bio_synthesis_run <- as.character(
  panel10_manifest$value[
    match(
      "BIO_SYNTHESIS_RUN",
      panel10_manifest$item
    )
  ]
)

catalog_path <- file.path(
  bio_synthesis_run,
  "plot_data",
  "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
  "00_INDEX",
  "FIG5_MASTER_PATHWAY_FEATURE_CATALOG.csv"
)

if (file.exists(catalog_path)) {
  pathway_catalog <- read_csv_checked(
    catalog_path,
    required = c(
      "family",
      "program_id",
      "program_short",
      "library",
      "pathway",
      "pathway_label",
      "phenotypic_side",
      "view_code",
      "feature_original",
      "feature_label",
      "pathway_feature_cosine"
    )
  )

  pathway_catalog$family <- as.character(
    pathway_catalog$family
  )

  pathway_catalog$view_code <- tolower(
    as.character(
      pathway_catalog$view_code
    )
  )

  pathway_catalog$feature_original <- as.character(
    pathway_catalog$feature_original
  )

  pathway_catalog$feature_key_final <- paste(
    pathway_catalog$view_code,
    pathway_catalog$feature_original,
    sep = "\r"
  )

  pathway_catalog$pathway_feature_cosine <- suppressWarnings(
    as.numeric(
      pathway_catalog$pathway_feature_cosine
    )
  )

  path_final <- merge(
    pathway_catalog[
      pathway_catalog$feature_key_final %in%
        small_panel$feature_key_final,
      ,
      drop = FALSE
    ],
    direction_summary[
      ,
      c(
        "feature_key_final",
        "framework_direction_class",
        "framework_pattern"
      ),
      drop = FALSE
    ],
    by = "feature_key_final",
    all.x = TRUE,
    sort = FALSE
  )

  path_final$target_phenotypic_side <- ifelse(
    is.na(
      path_final$framework_direction_class
    ),
    NA_character_,
    paste0(
      path_final$framework_direction_class,
      "_side"
    )
  )

  # Exact side match: do not force a pathway to a phenotype if Q5 did not
  # annotate that pathway on the same phenotypic side as the final feature.
  path_final_side <- path_final[
    !is.na(
      path_final$target_phenotypic_side
    ) &
    path_final$phenotypic_side ==
      path_final$target_phenotypic_side,
    ,
    drop = FALSE
  ]

  # ---------------------------------------------------------------------------
  # Feature-level top phenotype-aligned pathways
  # ---------------------------------------------------------------------------

  feature_path_rows <- lapply(
    split(
      direction_summary,
      direction_summary$feature_key_final
    ),
    function(zdir) {
      fk <- zdir$feature_key_final[1L]

      zz <- path_final_side[
        path_final_side$feature_key_final == fk,
        ,
        drop = FALSE
      ]

      if (!nrow(zz)) {
        return(
          data.frame(
            feature_key_final = fk,
            n_matching_side_pathways = 0L,
            top_matching_side_pathways =
              "NO_EXACT_SIDE_PATHWAY_LINK",
            top_programs =
              "NO_EXACT_SIDE_PROGRAM_LINK",
            stringsAsFactors = FALSE
          )
        )
      }

      zz$rank_strength <- abs(
        zz$pathway_feature_cosine
      )

      zz$rank_strength[
        !is.finite(
          zz$rank_strength
        )
      ] <- -Inf

      zz <- zz[
        order(
          -zz$rank_strength,
          zz$library,
          zz$pathway_label
        ),
        ,
        drop = FALSE
      ]

      unique_path_key <- paste(
        zz$family,
        zz$library,
        zz$pathway,
        sep = "\r"
      )

      zz_unique <- zz[
        !duplicated(
          unique_path_key
        ),
        ,
        drop = FALSE
      ]

      topn <- zz_unique[
        seq_len(
          min(
            3L,
            nrow(
              zz_unique
            )
          )
        ),
        ,
        drop = FALSE
      ]

      data.frame(
        feature_key_final = fk,
        n_matching_side_pathways =
          length(
            unique(
              paste(
                zz$family,
                zz$library,
                zz$pathway,
                sep = "\r"
              )
            )
          ),
        top_matching_side_pathways =
          paste(
            paste0(
              topn$library,
              "::",
              topn$pathway_label
            ),
            collapse = " | "
          ),
        top_programs =
          paste(
            unique(
              topn$program_short
            ),
            collapse = " | "
          ),
        stringsAsFactors = FALSE
      )
    }
  )

  feature_path_summary <- do.call(
    rbind,
    feature_path_rows
  )

  final_feature_summary <- merge(
    direction_summary,
    feature_path_summary,
    by = "feature_key_final",
    all.x = TRUE,
    sort = FALSE
  )

  final_feature_summary <- final_feature_summary[
    order(
      match(
        final_feature_summary$framework_direction_class,
        CLASS_LEVELS
      ),
      match(
        final_feature_summary$view_code,
        VIEW_LEVELS
      ),
      -final_feature_summary$median_abs_loading_across_models
    ),
    ,
    drop = FALSE
  ]

  cat("\nFINAL PANEL — DIRECTION + TOP PHENOTYPE-ALIGNED PATHWAYS\n")

  final_feature_print_cols <- c(
    "framework_direction_class",
    "view_code",
    "feature_label",
    "PLS_framework_concordant_n",
    "PLS_framework_evaluable_n",
    "n_matching_side_pathways",
    "top_matching_side_pathways",
    "top_programs"
  )

  final_feature_print_cols <- final_feature_print_cols[
    final_feature_print_cols %in%
      names(
        final_feature_summary
      )
  ]

  print(
    final_feature_summary[
      ,
      final_feature_print_cols,
      drop = FALSE
    ],
    row.names = FALSE
  )

  # ---------------------------------------------------------------------------
  # Pathway-level closing summary
  # ---------------------------------------------------------------------------

  if (nrow(path_final_side)) {
    path_final_side$abs_cosine <- abs(
      path_final_side$pathway_feature_cosine
    )

    pathway_key <- paste(
      path_final_side$framework_direction_class,
      path_final_side$library,
      path_final_side$pathway,
      path_final_side$program_id,
      sep = "\r"
    )

    pathway_split <- split(
      path_final_side,
      pathway_key
    )

    pathway_summary <- do.call(
      rbind,
      lapply(
        pathway_split,
        function(z) {
          data.frame(
            framework_direction_class =
              z$framework_direction_class[1L],
            library =
              z$library[1L],
            pathway =
              z$pathway[1L],
            pathway_label =
              z$pathway_label[1L],
            program_short =
              z$program_short[1L],
            n_final_features =
              length(
                unique(
                  z$feature_key_final
                )
              ),
            n_views =
              length(
                unique(
                  z$view_code
                )
              ),
            n_families =
              length(
                unique(
                  z$family
                )
              ),
            max_abs_pathway_feature_cosine =
              if (
                any(
                  is.finite(
                    z$abs_cosine
                  )
                )
              ) {
                max(
                  z$abs_cosine,
                  na.rm = TRUE
                )
              } else {
                NA_real_
              },
            features =
              paste(
                sort(
                  unique(
                    z$feature_label
                  )
                ),
                collapse = ";"
              ),
            stringsAsFactors = FALSE
          )
        }
      )
    )

    pathway_summary <- pathway_summary[
      order(
        match(
          pathway_summary$framework_direction_class,
          CLASS_LEVELS
        ),
        -pathway_summary$n_final_features,
        -pathway_summary$n_views,
        -pathway_summary$n_families,
        -pathway_summary$max_abs_pathway_feature_cosine,
        pathway_summary$pathway_label
      ),
      ,
      drop = FALSE
    ]

    cat(
      "\nPHENOTYPE-ALIGNED PATHWAY SUMMARY\n",
      "NOTE: for ME features these are frozen cross-omic pathway associations, ",
      "NOT direct chemical identification of the untargeted metabolite.\n",
      sep = ""
    )

    for (cls in CLASS_LEVELS) {
      z <- pathway_summary[
        pathway_summary$framework_direction_class == cls,
        ,
        drop = FALSE
      ]

      cat(
        "\n--- ",
        cls,
        " direction: top pathways ---\n",
        sep = ""
      )

      if (!nrow(z)) {
        cat("No exact-side pathway links.\n")
        next
      }

      print(
        z[
          seq_len(
            min(
              10L,
              nrow(z)
            )
          ),
          c(
            "library",
            "pathway_label",
            "program_short",
            "n_final_features",
            "n_views",
            "n_families",
            "max_abs_pathway_feature_cosine",
            "features"
          ),
          drop = FALSE
        ],
        row.names = FALSE
      )
    }

    # Metabolomic final features only.
    me_pathway_summary <- path_final_side[
      path_final_side$view_code == "me",
      ,
      drop = FALSE
    ]

    cat(
      "\nMETABOLOMIC FINAL FEATURES — PHENOTYPE-ALIGNED CROSS-OMIC PATHWAYS\n",
      "(These are pathway associations, not direct metabolite IDs.)\n",
      sep = ""
    )

    for (cls in CLASS_LEVELS) {
      zz <- me_pathway_summary[
        me_pathway_summary$framework_direction_class == cls,
        ,
        drop = FALSE
      ]

      if (!nrow(zz)) {
        cat(
          cls,
          ": no exact-side pathway links.\n",
          sep = ""
        )
        next
      }

      split_me <- split(
        zz,
        paste(
          zz$library,
          zz$pathway,
          zz$program_id,
          sep = "\r"
        )
      )

      me_sum <- do.call(
        rbind,
        lapply(
          split_me,
          function(a) {
            data.frame(
              library =
                a$library[1L],
              pathway_label =
                a$pathway_label[1L],
              program_short =
                a$program_short[1L],
              n_ME_final_features =
                length(
                  unique(
                    a$feature_key_final
                  )
                ),
              n_families =
                length(
                  unique(
                    a$family
                  )
                ),
              max_abs_cosine =
                if (
                  any(
                    is.finite(
                      a$pathway_feature_cosine
                    )
                  )
                ) {
                  max(
                    abs(
                      a$pathway_feature_cosine
                    ),
                    na.rm = TRUE
                  )
                } else {
                  NA_real_
                },
              ME_features =
                paste(
                  sort(
                    unique(
                      a$feature_label
                    )
                  ),
                  collapse = ";"
                ),
              stringsAsFactors = FALSE
            )
          }
        )
      )

      me_sum <- me_sum[
        order(
          -me_sum$n_ME_final_features,
          -me_sum$n_families,
          -me_sum$max_abs_cosine,
          me_sum$pathway_label
        ),
        ,
        drop = FALSE
      ]

      cat(
        "\n",
        cls,
        " metabolomic-side summary:\n",
        sep = ""
      )

      print(
        me_sum[
          seq_len(
            min(
              10L,
              nrow(
                me_sum
              )
            )
          ),
          ,
          drop = FALSE
        ],
        row.names = FALSE
      )
    }

    # -------------------------------------------------------------------------
    # Two final RStudio plots; nothing is saved.
    # -------------------------------------------------------------------------

    direction_plot_df <- as.data.frame(
      table(
        factor(
          direction_summary$framework_direction_class,
          levels = CLASS_LEVELS
        ),
        factor(
          toupper(
            direction_summary$view_code
          ),
          levels = toupper(
            VIEW_LEVELS
          )
        )
      ),
      stringsAsFactors = FALSE
    )

    names(direction_plot_df) <- c(
      "Direction",
      "Omic",
      "n_features"
    )

    direction_plot_df <- direction_plot_df[
      direction_plot_df$n_features > 0L,
      ,
      drop = FALSE
    ]

    p_final_direction <- ggplot2::ggplot(
      direction_plot_df,
      ggplot2::aes(
        x = Direction,
        y = n_features,
        fill = Omic
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::labs(
        title =
          "Final 27-feature panel — phenotype direction",
        subtitle =
          "Direction comes from frozen bidirectional Evidence-A framework.",
        x = "Framework phenotype-high direction",
        y = "Final features",
        fill = "Omic"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor =
          ggplot2::element_blank(),
        plot.title =
          ggplot2::element_text(
            face = "bold"
          ),
        legend.position =
          "bottom"
      )

    print(
      p_final_direction
    )

    top_path_plot <- do.call(
      rbind,
      lapply(
        CLASS_LEVELS,
        function(cls) {
          z <- pathway_summary[
            pathway_summary$framework_direction_class == cls,
            ,
            drop = FALSE
          ]

          if (!nrow(z)) {
            return(NULL)
          }

          z[
            seq_len(
              min(
                8L,
                nrow(z)
              )
            ),
            ,
            drop = FALSE
          ]
        }
      )
    )

    if (
      !is.null(
        top_path_plot
      ) &&
      nrow(
        top_path_plot
      )
    ) {
      top_path_plot$display_pathway <- paste0(
        top_path_plot$framework_direction_class,
        " | ",
        top_path_plot$pathway_label
      )

      top_path_plot$display_pathway <- factor(
        top_path_plot$display_pathway,
        levels = rev(
          unique(
            top_path_plot$display_pathway
          )
        )
      )

      p_final_pathways <- ggplot2::ggplot(
        top_path_plot,
        ggplot2::aes(
          x = display_pathway,
          y = n_final_features,
          fill = library
        )
      ) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::facet_wrap(
          ~ framework_direction_class,
          scales = "free_y"
        ) +
        ggplot2::labs(
          title =
            "Final panel — phenotype-aligned pathway summary",
          subtitle =
            "Existing frozen Q5 feature-pathway links; no enrichment was recalculated.",
          x = NULL,
          y = "Number of final features linked",
          fill = "Library"
        ) +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
          panel.grid.minor =
            ggplot2::element_blank(),
          plot.title =
            ggplot2::element_text(
              face = "bold"
            ),
          legend.position =
            "bottom"
        )

      print(
        p_final_pathways
      )
    }
  } else {
    cat(
      "\nNo exact phenotypic-side pathway links were found for the final panel.\n"
    )
  }
} else {
  cat(
    "\nWARNING: frozen Q5 pathway-feature catalog was not found:\n",
    catalog_path,
    "\nPathway summary skipped; no new enrichment was attempted.\n",
    sep = ""
  )
}

cat("\n============================================================\n")
cat("PROJECT-CLOSING INTERPRETATION OBJECTS CREATED IN MEMORY ONLY\n")
cat("  clinical_audit        -> P1 clinical exclusion audit\n")
cat("  direction_summary     -> final feature direction + PLS corroboration\n")
cat("  final_feature_summary -> final feature + phenotype-aligned pathways (if Q5 catalog found)\n")
cat("  pathway_summary       -> pathway-level final summary (if Q5 catalog found)\n")
cat("NO FILES WERE WRITTEN BY THIS CLOSING SECTION.\n")
cat("============================================================\n")

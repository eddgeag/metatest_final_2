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

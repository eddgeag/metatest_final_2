#!/usr/bin/env Rscript

# ==============================================================================
# METATEST — 18_build_final_consensus_panel.R
#
# PURPOSE
# -------
# Standalone final integration step. It does NOT rerun PLS-DA or redefine
# Evidence A. It starts from the completed script-17 output and builds:
#
#   CORE      = strict Evidence-A feature present in P1 + P2 + P3
#   SHARED_2  = strict Evidence-A feature present in exactly 2 analyses
#   SPECIFIC  = strict Evidence-A feature present in exactly 1 analysis
#
# Final biological signatures exported:
#   CORE panel      = 3/3 analyses
#   EXTENDED panel  = CORE + SHARED_2 (>=2/3 analyses)
#
# PLS-DA VIP/loadings and Sparse Group Lasso membership are annotations only.
# They do NOT redefine CORE / SHARED_2 / SPECIFIC.
#
# IMPORTANT
# ---------
# Loading signs are interpretable within each fitted PLS-DA model/component.
# They are NOT averaged across P1/P2/P3 because latent-component signs can flip
# between separately fitted models.
# ==============================================================================

options(stringsAsFactors = FALSE)

# ==============================================================================
# 0. CONFIGURATION
# ==============================================================================

UPSTREAM_DIR <- Sys.getenv(
  "UPSTREAM_DIR",
  unset = path.expand(
    "~/Documents/PhD/metatest/results_upstream_local_20260908/comparacion_modelamientos/upstream_comun"
  )
)

PANEL17_RUN <- Sys.getenv("PANEL17_RUN", unset = "")
BIO_PROGRAM_RUN <- Sys.getenv("BIO_PROGRAM_RUN", unset = "")
SGL_RUN <- Sys.getenv("SGL_RUN", unset = "")
OUTPUT_ROOT_ENV <- Sys.getenv("PANEL18_OUTPUT_ROOT", unset = "")

FAMILIES <- c("P1", "P2", "P3")
VIEW_LEVELS <- c("tx", "pr", "me", "cl")
VIP_THRESHOLD <- as.numeric(Sys.getenv("PANEL18_VIP_THRESHOLD", unset = "1"))
LOADING_BAR_TOPN <- as.integer(Sys.getenv("PANEL18_LOADING_BAR_TOPN", unset = "25"))
LOADING_MAP_TOPN <- as.integer(Sys.getenv("PANEL18_LOADING_MAP_TOPN", unset = "30"))
DPI <- as.integer(Sys.getenv("PANEL18_DPI", unset = "180"))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.", call. = FALSE)
}

# ==============================================================================
# 1. HELPERS
# ==============================================================================

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
  assert(file.exists(path), paste("Missing required file:", path))

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

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

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

first_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x) & x != "NA"]
  if (!length(x)) return(NA_character_)
  x[1L]
}

theme_panel <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

find_latest_selection_run <- function(upstream) {
  root <- file.path(
    upstream,
    "downstream",
    "seleccion_mejor_escenario_train"
  )

  assert(dir.exists(root), paste("Selection root does not exist:", root))

  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[
    file.exists(
      file.path(
        dirs,
        "SELECTED_SCENARIOS_BY_ANALYSIS.csv"
      )
    )
  ]

  assert(length(dirs) > 0L, "No complete selection run found.")
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

find_latest_bio_run <- function(selection_run) {
  if (nzchar(BIO_PROGRAM_RUN)) {
    return(canonical(BIO_PROGRAM_RUN))
  }

  root <- file.path(selection_run, "BIOLOGICAL_PROGRAMS")
  assert(dir.exists(root), paste("BIOLOGICAL_PROGRAMS does not exist:", root))

  latest <- file.path(root, "LATEST_RUN.txt")

  if (file.exists(latest)) {
    p <- trimws(readLines(latest, warn = FALSE)[1L])

    if (
      nzchar(p) &&
      dir.exists(p) &&
      file.exists(file.path(p, "BIOLOGICAL_PROGRAMS_OBJECT.rds"))
    ) {
      return(canonical(p))
    }
  }

  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[
    file.exists(
      file.path(
        dirs,
        "BIOLOGICAL_PROGRAMS_OBJECT.rds"
      )
    )
  ]

  assert(length(dirs) > 0L, "No complete BIOLOGICAL_PROGRAMS run found.")
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

is_complete_panel17_run <- function(p) {
  if (!nzchar(p) || !dir.exists(p)) return(FALSE)

  req <- c(
    file.path("tables", "02_EVIDENCE_A_FEATURES_FOR_PLSDA_ALL_FAMILIES.csv"),
    file.path("tables", "10_PLSDA_VIP_PRIORITY_ALL_FAMILIES.csv"),
    file.path("tables", "11_PLSDA_LOADINGS_ALL_FAMILIES.csv")
  )

  all(file.exists(file.path(p, req)))
}

find_panel17_run <- function(bio_run) {
  if (nzchar(PANEL17_RUN)) {
    p <- canonical(PANEL17_RUN)
    assert(is_complete_panel17_run(p), "PANEL17_RUN is incomplete.")
    return(p)
  }

  root <- file.path(bio_run, "FINAL_EVIDENCE_A_PLSDA_PCA")
  assert(dir.exists(root), paste("Missing FINAL_EVIDENCE_A_PLSDA_PCA:", root))

  latest <- file.path(root, "LATEST_RUN.txt")

  if (file.exists(latest)) {
    p <- trimws(readLines(latest, warn = FALSE)[1L])
    if (is_complete_panel17_run(p)) return(canonical(p))
  }

  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[
    vapply(
      dirs,
      is_complete_panel17_run,
      logical(1)
    )
  ]

  assert(length(dirs) > 0L, "No complete script-17 run found.")
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

find_optional_sgl_run <- function(bio_run) {
  if (nzchar(SGL_RUN)) {
    p <- canonical(SGL_RUN)
    f <- file.path(
      p,
      "tables",
      "08_FINAL_SPARSE_GROUP_PANEL_ALL_FAMILIES.csv"
    )

    if (file.exists(f)) return(p)

    warning("SGL_RUN supplied but final SGL panel table is absent; continuing without SGL.")
    return(NA_character_)
  }

  root <- file.path(
    bio_run,
    "FINAL_SPARSE_GROUP_LASSO_PANEL"
  )

  if (!dir.exists(root)) return(NA_character_)

  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)

  dirs <- dirs[
    file.exists(
      file.path(
        dirs,
        "tables",
        "08_FINAL_SPARSE_GROUP_PANEL_ALL_FAMILIES.csv"
      )
    )
  ]

  if (!length(dirs)) return(NA_character_)

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

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

# ==============================================================================
# 2. LOCATE INPUTS
# ==============================================================================

SELECTION <- find_latest_selection_run(UPSTREAM_DIR)
BIO_RUN <- find_latest_bio_run(SELECTION)
PANEL17 <- find_panel17_run(BIO_RUN)
SGL <- find_optional_sgl_run(BIO_RUN)

if (nzchar(OUTPUT_ROOT_ENV)) {
  OUT_ROOT <- path.expand(OUTPUT_ROOT_ENV)
} else {
  OUT_ROOT <- file.path(
    BIO_RUN,
    "FINAL_CONSENSUS_PANEL",
    paste0(
      format(Sys.time(), "%Y%m%dT%H%M%S"),
      "_pid",
      Sys.getpid()
    )
  )
}

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
OUT_ROOT <- canonical(OUT_ROOT)

DIR_TABLES <- file.path(OUT_ROOT, "tables")
DIR_FIGURES <- file.path(OUT_ROOT, "figures")
dir.create(DIR_TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIGURES, recursive = TRUE, showWarnings = FALSE)

message("============================================================")
message("METATEST — FINAL CONSENSUS PANEL")
message("BIOLOGICAL_PROGRAMS run: ", BIO_RUN)
message("Script-17 run: ", PANEL17)
message("Sparse Group Lasso run: ", ifelse(is.na(SGL), "NOT FOUND / OPTIONAL", SGL))
message("Output: ", OUT_ROOT)
message("============================================================")

# ==============================================================================
# 3. READ SCRIPT-17 TABLES
# ==============================================================================

selected <- read_csv_checked(
  file.path(
    PANEL17,
    "tables",
    "02_EVIDENCE_A_FEATURES_FOR_PLSDA_ALL_FAMILIES.csv"
  ),
  required = c(
    "family",
    "predictor_key",
    "view_code",
    "feature_original",
    "feature_label"
  )
)

vip <- read_csv_checked(
  file.path(
    PANEL17,
    "tables",
    "10_PLSDA_VIP_PRIORITY_ALL_FAMILIES.csv"
  ),
  required = c(
    "family",
    "predictor_key",
    "view_code",
    "feature_original",
    "feature_label",
    "VIP_max_selected_components"
  )
)

loadings <- read_csv_checked(
  file.path(
    PANEL17,
    "tables",
    "11_PLSDA_LOADINGS_ALL_FAMILIES.csv"
  ),
  required = c(
    "family",
    "predictor_key",
    "component",
    "loading",
    "abs_loading",
    "view_code",
    "feature_original",
    "feature_label"
  )
)

selected$family <- as.character(selected$family)
selected$view_code <- tolower(as.character(selected$view_code))
selected$feature_original <- as.character(selected$feature_original)
selected$feature_label <- as.character(selected$feature_label)

vip$family <- as.character(vip$family)
vip$view_code <- tolower(as.character(vip$view_code))
vip$feature_original <- as.character(vip$feature_original)
vip$VIP_max_selected_components <- safe_num(vip$VIP_max_selected_components)

loadings$family <- as.character(loadings$family)
loadings$view_code <- tolower(as.character(loadings$view_code))
loadings$feature_original <- as.character(loadings$feature_original)
loadings$component <- as.integer(loadings$component)
loadings$loading <- safe_num(loadings$loading)
loadings$abs_loading <- safe_num(loadings$abs_loading)

selected <- selected[
  selected$family %in% FAMILIES &
  selected$view_code %in% VIEW_LEVELS,
  ,
  drop = FALSE
]

# One row per family x measurable feature.
selected$key <- paste(
  selected$view_code,
  selected$feature_original,
  sep = "\r"
)

selected_family_unique <- selected[
  !duplicated(
    paste(
      selected$family,
      selected$key,
      sep = "\r"
    )
  ),
  ,
  drop = FALSE
]

# ==============================================================================
# 4. CONSENSUS: CORE / SHARED_2 / SPECIFIC
# ==============================================================================

split_features <- split(
  selected_family_unique,
  selected_family_unique$key
)

master_rows <- lapply(
  split_features,
  function(z) {
    fam <- sort(unique(as.character(z$family)))

    in_P1 <- "P1" %in% fam
    in_P2 <- "P2" %in% fam
    in_P3 <- "P3" %in% fam
    n_fam <- length(fam)

    exact_intersection <- if (
      in_P1 && in_P2 && in_P3
    ) {
      "P1&P2&P3"
    } else if (
      in_P1 && in_P2
    ) {
      "P1&P2"
    } else if (
      in_P1 && in_P3
    ) {
      "P1&P3"
    } else if (
      in_P2 && in_P3
    ) {
      "P2&P3"
    } else if (in_P1) {
      "P1 only"
    } else if (in_P2) {
      "P2 only"
    } else {
      "P3 only"
    }

    tier <- if (
      n_fam == 3L
    ) {
      "CORE"
    } else if (
      n_fam == 2L
    ) {
      "SHARED_2"
    } else {
      "SPECIFIC"
    }

    data.frame(
      feature_key_consensus = z$key[1L],
      view_code = first_nonempty(z$view_code),
      feature_original = first_nonempty(z$feature_original),
      feature_label = first_nonempty(z$feature_label),
      n_families = n_fam,
      families = paste(fam, collapse = ";"),
      in_P1 = in_P1,
      in_P2 = in_P2,
      in_P3 = in_P3,
      exact_intersection = exact_intersection,
      consensus_tier = tier,
      in_CORE_panel = n_fam == 3L,
      in_EXTENDED_panel = n_fam >= 2L,
      stringsAsFactors = FALSE
    )
  }
)

panel_master <- do.call(rbind, master_rows)
rownames(panel_master) <- NULL

# ==============================================================================
# 5. VIP ANNOTATION — CONTRIBUTION ONLY, NOT MEMBERSHIP
# ==============================================================================

vip$key <- paste(
  vip$view_code,
  vip$feature_original,
  sep = "\r"
)

vip$VIP_ge_threshold <- vip$VIP_max_selected_components >= VIP_THRESHOLD

for (family in FAMILIES) {
  z <- vip[
    vip$family == family,
    ,
    drop = FALSE
  ]

  # Defensive aggregation if duplicates ever occur.
  if (nrow(z)) {
    z_split <- split(z, z$key)

    z_agg <- do.call(
      rbind,
      lapply(
        z_split,
        function(a) {
          data.frame(
            feature_key_consensus = a$key[1L],
            vip = max(
              a$VIP_max_selected_components,
              na.rm = TRUE
            ),
            vip_ge = any(
              a$VIP_ge_threshold,
              na.rm = TRUE
            ),
            stringsAsFactors = FALSE
          )
        }
      )
    )

    names(z_agg)[names(z_agg) == "vip"] <- paste0("VIP_", family)
    names(z_agg)[names(z_agg) == "vip_ge"] <- paste0("VIP_GE_", family)

    panel_master <- merge(
      panel_master,
      z_agg,
      by = "feature_key_consensus",
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    panel_master[[paste0("VIP_", family)]] <- NA_real_
    panel_master[[paste0("VIP_GE_", family)]] <- FALSE
  }
}

for (family in FAMILIES) {
  nm <- paste0("VIP_GE_", family)
  panel_master[[nm]][is.na(panel_master[[nm]])] <- FALSE
}

panel_master$n_VIP_ge_threshold <- rowSums(
  panel_master[
    ,
    paste0("VIP_GE_", FAMILIES),
    drop = FALSE
  ]
)

vip_cols <- paste0("VIP_", FAMILIES)

panel_master$VIP_max_across_present <- apply(
  panel_master[, vip_cols, drop = FALSE],
  1L,
  function(x) {
    x <- safe_num(x)
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    max(x)
  }
)

panel_master$VIP_median_across_present <- apply(
  panel_master[, vip_cols, drop = FALSE],
  1L,
  function(x) {
    x <- safe_num(x)
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    stats::median(x)
  }
)

# Secondary descriptive subset only.
panel_master$high_contribution_CORE_VIP_all3 <-
  panel_master$in_CORE_panel &
  panel_master$n_VIP_ge_threshold == 3L

# ==============================================================================
# 6. LOADINGS ANNOTATION
# ==============================================================================
#
# Signs are kept family/component-specific. We do not calculate a cross-analysis
# signed mean because each PLS component may be sign-reflected independently.
# ==============================================================================

loadings$key <- paste(
  loadings$view_code,
  loadings$feature_original,
  sep = "\r"
)

available_components <- sort(
  unique(
    loadings$component[
      is.finite(loadings$component)
    ]
  )
)

for (family in FAMILIES) {
  for (cc in available_components) {
    z <- loadings[
      loadings$family == family &
      loadings$component == cc,
      ,
      drop = FALSE
    ]

    if (!nrow(z)) next

    z_split <- split(z, z$key)

    z_agg <- do.call(
      rbind,
      lapply(
        z_split,
        function(a) {
          ii <- which.max(a$abs_loading)

          data.frame(
            feature_key_consensus = a$key[1L],
            loading = a$loading[ii],
            abs_loading = a$abs_loading[ii],
            stringsAsFactors = FALSE
          )
        }
      )
    )

    names(z_agg)[names(z_agg) == "loading"] <-
      paste0("loading_", family, "_C", cc)

    names(z_agg)[names(z_agg) == "abs_loading"] <-
      paste0("abs_loading_", family, "_C", cc)

    panel_master <- merge(
      panel_master,
      z_agg,
      by = "feature_key_consensus",
      all.x = TRUE,
      sort = FALSE
    )
  }
}

abs_loading_cols <- grep(
  "^abs_loading_P[123]_C[0-9]+$",
  names(panel_master),
  value = TRUE
)

panel_master$max_abs_loading_across_models <- apply(
  panel_master[, abs_loading_cols, drop = FALSE],
  1L,
  function(x) {
    x <- safe_num(x)
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    max(x)
  }
)

panel_master$median_abs_loading_across_models <- apply(
  panel_master[, abs_loading_cols, drop = FALSE],
  1L,
  function(x) {
    x <- safe_num(x)
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    stats::median(x)
  }
)

# ==============================================================================
# 7. OPTIONAL SPARSE GROUP LASSO CORROBORATION
# ==============================================================================

panel_master$SGL_P1 <- FALSE
panel_master$SGL_P2 <- FALSE
panel_master$SGL_P3 <- FALSE

if (!is.na(SGL)) {
  sgl <- read_csv_checked(
    file.path(
      SGL,
      "tables",
      "08_FINAL_SPARSE_GROUP_PANEL_ALL_FAMILIES.csv"
    ),
    required = c(
      "family",
      "view_code",
      "feature_original"
    )
  )

  sgl$family <- as.character(sgl$family)
  sgl$view_code <- tolower(as.character(sgl$view_code))
  sgl$feature_original <- as.character(sgl$feature_original)

  sgl$key <- paste(
    sgl$view_code,
    sgl$feature_original,
    sep = "\r"
  )

  for (family in FAMILIES) {
    keys_here <- unique(
      sgl$key[
        sgl$family == family
      ]
    )

    panel_master[[paste0("SGL_", family)]] <-
      panel_master$feature_key_consensus %in% keys_here
  }
}

panel_master$n_SGL_selected <- rowSums(
  panel_master[
    ,
    paste0("SGL_", FAMILIES),
    drop = FALSE
  ]
)

# ==============================================================================
# 8. PANEL FILES
# ==============================================================================

tier_order <- c("CORE", "SHARED_2", "SPECIFIC")

panel_master$consensus_tier <- factor(
  panel_master$consensus_tier,
  levels = tier_order
)

panel_master <- panel_master[
  order(
    panel_master$consensus_tier,
    match(
      panel_master$view_code,
      VIEW_LEVELS
    ),
    -panel_master$n_VIP_ge_threshold,
    -panel_master$VIP_max_across_present,
    panel_master$feature_label
  ),
  ,
  drop = FALSE
]

panel_master$consensus_tier <- as.character(
  panel_master$consensus_tier
)

core_panel <- panel_master[
  panel_master$in_CORE_panel,
  ,
  drop = FALSE
]

extended_panel <- panel_master[
  panel_master$in_EXTENDED_panel,
  ,
  drop = FALSE
]

specific_panel <- panel_master[
  panel_master$consensus_tier == "SPECIFIC",
  ,
  drop = FALSE
]

core_vip_all3 <- panel_master[
  panel_master$high_contribution_CORE_VIP_all3,
  ,
  drop = FALSE
]

write_csv(
  panel_master,
  file.path(
    DIR_TABLES,
    "00_FINAL_PANEL_MASTER.csv"
  )
)

write_csv(
  core_panel,
  file.path(
    DIR_TABLES,
    "01_FINAL_PANEL_CORE_3OF3.csv"
  )
)

write_csv(
  extended_panel,
  file.path(
    DIR_TABLES,
    "02_FINAL_PANEL_EXTENDED_GE2OF3.csv"
  )
)

write_csv(
  specific_panel,
  file.path(
    DIR_TABLES,
    "03_CONTEXT_SPECIFIC_1OF3.csv"
  )
)

write_csv(
  core_vip_all3,
  file.path(
    DIR_TABLES,
    "04_CORE_HIGH_CONTRIBUTION_VIP_ALL3_DESCRIPTIVE.csv"
  )
)

# Long loadings restricted to the extended final signature.
extended_keys <- extended_panel$feature_key_consensus

final_loadings_long <- loadings[
  loadings$key %in% extended_keys,
  ,
  drop = FALSE
]

final_loadings_long$consensus_tier <-
  panel_master$consensus_tier[
    match(
      final_loadings_long$key,
      panel_master$feature_key_consensus
    )
  ]

write_csv(
  final_loadings_long,
  file.path(
    DIR_TABLES,
    "05_FINAL_EXTENDED_PANEL_LOADINGS_LONG.csv"
  )
)

# ==============================================================================
# 9. COUNTS
# ==============================================================================

count_summary <- data.frame(
  panel = c(
    "CORE_3OF3",
    "EXTENDED_GE2OF3",
    "SPECIFIC_1OF3",
    "CORE_VIP_ALL3_DESCRIPTIVE"
  ),
  n_features = c(
    nrow(core_panel),
    nrow(extended_panel),
    nrow(specific_panel),
    nrow(core_vip_all3)
  ),
  stringsAsFactors = FALSE
)

write_csv(
  count_summary,
  file.path(
    DIR_TABLES,
    "06_FINAL_PANEL_COUNTS.csv"
  )
)

count_by_view <- as.data.frame(
  table(
    factor(
      panel_master$consensus_tier,
      levels = tier_order
    ),
    factor(
      toupper(panel_master$view_code),
      levels = toupper(VIEW_LEVELS)
    )
  ),
  stringsAsFactors = FALSE
)

names(count_by_view) <- c(
  "consensus_tier",
  "view_code",
  "n_features"
)

count_by_view <- count_by_view[
  count_by_view$n_features > 0L,
  ,
  drop = FALSE
]

write_csv(
  count_by_view,
  file.path(
    DIR_TABLES,
    "07_FINAL_PANEL_COUNTS_BY_OMIC.csv"
  )
)

# ==============================================================================
# 10. FIGURE — CONSENSUS PRESENCE HEATMAP
# ==============================================================================

heat <- do.call(
  rbind,
  lapply(
    seq_len(nrow(panel_master)),
    function(i) {
      data.frame(
        feature_key_consensus =
          panel_master$feature_key_consensus[i],
        display_label =
          paste0(
            toupper(panel_master$view_code[i]),
            " | ",
            panel_master$feature_label[i]
          ),
        consensus_tier =
          panel_master$consensus_tier[i],
        family = FAMILIES,
        present = c(
          panel_master$in_P1[i],
          panel_master$in_P2[i],
          panel_master$in_P3[i]
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

# Put CORE first, then SHARED_2, then SPECIFIC, preserving omic ordering.
ordered_labels <- unique(
  paste0(
    toupper(panel_master$view_code),
    " | ",
    panel_master$feature_label
  )
)

heat$display_label <- factor(
  heat$display_label,
  levels = rev(ordered_labels)
)

heat$present_label <- ifelse(
  heat$present,
  "Strict Evidence-A",
  "Absent"
)

p_heat <- ggplot2::ggplot(
  heat,
  ggplot2::aes(
    x = family,
    y = display_label,
    fill = present_label
  )
) +
  ggplot2::geom_tile(
    linewidth = 0.15
  ) +
  ggplot2::facet_grid(
    consensus_tier ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::labs(
    title =
      "Final strict Evidence-A consensus across P1, P2 and P3",
    subtitle =
      paste0(
        "CORE=3/3; SHARED_2=2/3; SPECIFIC=1/3. ",
        "Membership is defined before considering VIP or Sparse Group Lasso."
      ),
    x = "Analysis",
    y = NULL,
    fill = "Membership"
  ) +
  theme_panel(8) +
  ggplot2::theme(
    axis.text.y =
      ggplot2::element_text(size = 5.5),
    panel.spacing.y =
      grid::unit(0.6, "lines")
  )

ggplot2::ggsave(
  file.path(
    DIR_FIGURES,
    "01_FINAL_PANEL_CONSENSUS_HEATMAP.png"
  ),
  p_heat,
  width = 8,
  height = max(
    10,
    0.12 * length(ordered_labels) + 4
  ),
  limitsize = FALSE,
  dpi = DPI
)

# ==============================================================================
# 11. FIGURE — OMIC COMPOSITION BY CONSENSUS TIER
# ==============================================================================

plot_counts <- count_by_view
plot_counts$consensus_tier <- factor(
  plot_counts$consensus_tier,
  levels = tier_order
)

p_comp <- ggplot2::ggplot(
  plot_counts,
  ggplot2::aes(
    x = consensus_tier,
    y = n_features,
    fill = view_code
  )
) +
  ggplot2::geom_col() +
  ggplot2::labs(
    title =
      "Final consensus signature — omic composition",
    subtitle =
      "CORE and SHARED_2 form the extended final transversal signature.",
    x = "Consensus tier",
    y = "Number of features",
    fill = "Omic"
  ) +
  theme_panel(11)

ggplot2::ggsave(
  file.path(
    DIR_FIGURES,
    "02_FINAL_PANEL_COMPOSITION_BY_OMIC.png"
  ),
  p_comp,
  width = 8,
  height = 6,
  dpi = DPI
)

# ==============================================================================
# 12. FIGURE — VIP RECURRENCE IN FINAL EXTENDED PANEL
# ==============================================================================

vip_recurrence <- as.data.frame(
  table(
    factor(
      extended_panel$n_VIP_ge_threshold,
      levels = 0:3
    ),
    factor(
      toupper(extended_panel$view_code),
      levels = toupper(VIEW_LEVELS)
    )
  ),
  stringsAsFactors = FALSE
)

names(vip_recurrence) <- c(
  "n_analyses_VIP_ge_threshold",
  "view_code",
  "n_features"
)

vip_recurrence <- vip_recurrence[
  vip_recurrence$n_features > 0L,
  ,
  drop = FALSE
]

p_vip_rec <- ggplot2::ggplot(
  vip_recurrence,
  ggplot2::aes(
    x = n_analyses_VIP_ge_threshold,
    y = n_features,
    fill = view_code
  )
) +
  ggplot2::geom_col() +
  ggplot2::labs(
    title =
      "PLS-DA contribution recurrence within the final extended signature",
    subtitle =
      paste0(
        "Number of P1/P2/P3 models in which VIP >= ",
        VIP_THRESHOLD,
        ". VIP is descriptive and does not redefine final-panel membership."
      ),
    x = paste0("Number of analyses with VIP >= ", VIP_THRESHOLD),
    y = "Number of final features",
    fill = "Omic"
  ) +
  theme_panel(11)

ggplot2::ggsave(
  file.path(
    DIR_FIGURES,
    "03_FINAL_PANEL_VIP_RECURRENCE_BY_OMIC.png"
  ),
  p_vip_rec,
  width = 8,
  height = 6,
  dpi = DPI
)

# ==============================================================================
# 13. SIGNED LOADING BARS — FINAL EXTENDED PANEL
# ==============================================================================

DIR_LOADINGS <- file.path(
  DIR_FIGURES,
  "LOADINGS_FINAL_EXTENDED"
)

dir.create(
  DIR_LOADINGS,
  recursive = TRUE,
  showWarnings = FALSE
)

for (family in FAMILIES) {
  fam_dir <- file.path(
    DIR_LOADINGS,
    family
  )

  dir.create(
    fam_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  comps_here <- sort(
    unique(
      final_loadings_long$component[
        final_loadings_long$family == family
      ]
    )
  )

  comps_here <- comps_here[
    is.finite(comps_here)
  ]

  for (cc in comps_here) {
    z <- final_loadings_long[
      final_loadings_long$family == family &
      final_loadings_long$component == cc,
      ,
      drop = FALSE
    ]

    if (!nrow(z)) next

    z <- z[
      order(
        -z$abs_loading,
        z$feature_label
      ),
      ,
      drop = FALSE
    ]

    z <- z[
      seq_len(
        min(
          LOADING_BAR_TOPN,
          nrow(z)
        )
      ),
      ,
      drop = FALSE
    ]

    z$direction <- ifelse(
      z$loading >= 0,
      "Positive",
      "Negative"
    )

    z$display_label <- paste0(
      toupper(z$view_code),
      " | ",
      z$feature_label
    )

    z <- z[
      order(z$loading),
      ,
      drop = FALSE
    ]

    z$display_label <- factor(
      z$display_label,
      levels = z$display_label
    )

    p_bar <- ggplot2::ggplot(
      z,
      ggplot2::aes(
        x = display_label,
        y = loading,
        fill = direction
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = 2
      ) +
      ggplot2::labs(
        title =
          paste0(
            family,
            " — final extended panel loadings: component ",
            cc
          ),
        subtitle =
          paste0(
            "Top ",
            nrow(z),
            " CORE+SHARED_2 features. ",
            "Sign gives within-model direction; magnitude gives loading contribution."
          ),
        x = NULL,
        y = "PLS loading",
        fill = "Direction"
      ) +
      theme_panel(9)

    ggplot2::ggsave(
      file.path(
        fam_dir,
        paste0(
          family,
          "_FINAL_EXTENDED_LOADINGS_COMP",
          cc,
          "_TOP",
          nrow(z),
          ".png"
        )
      ),
      p_bar,
      width = 9,
      height = max(
        6,
        0.24 * nrow(z) + 2
      ),
      dpi = DPI
    )
  }
}

# ==============================================================================
# 14. LOADING MAPS — BIPLOT-LIKE FEATURE MAP
# ==============================================================================
#
# Feature loadings only, not sample scores. This avoids an overloaded classical
# biplot while showing which final features jointly contribute to two components.
# ==============================================================================

make_loading_map <- function(family, comp_x, comp_y) {
  dx <- final_loadings_long[
    final_loadings_long$family == family &
    final_loadings_long$component == comp_x,
    ,
    drop = FALSE
  ]

  dy <- final_loadings_long[
    final_loadings_long$family == family &
    final_loadings_long$component == comp_y,
    ,
    drop = FALSE
  ]

  if (!nrow(dx) || !nrow(dy)) return(invisible(FALSE))

  m <- merge(
    dx[
      ,
      c(
        "key",
        "view_code",
        "feature_label",
        "loading",
        "abs_loading"
      )
    ],
    dy[
      ,
      c(
        "key",
        "loading",
        "abs_loading"
      )
    ],
    by = "key",
    suffixes = c("_x", "_y"),
    all = FALSE,
    sort = FALSE
  )

  if (!nrow(m)) return(invisible(FALSE))

  m$importance <- pmax(
    m$abs_loading_x,
    m$abs_loading_y
  )

  m <- m[
    order(
      -m$importance,
      m$feature_label
    ),
    ,
    drop = FALSE
  ]

  m <- m[
    seq_len(
      min(
        LOADING_MAP_TOPN,
        nrow(m)
      )
    ),
    ,
    drop = FALSE
  ]

  p_map <- ggplot2::ggplot(
    m,
    ggplot2::aes(
      x = loading_x,
      y = loading_y,
      colour = view_code
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    ggplot2::geom_point(
      size = 2.7,
      alpha = 0.90
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = feature_label
      ),
      size = 3,
      hjust = -0.05,
      vjust = 0.4,
      check_overlap = TRUE
    ) +
    ggplot2::labs(
      title =
        paste0(
          family,
          " — final extended loading map: comp ",
          comp_x,
          " vs comp ",
          comp_y
        ),
      subtitle =
        paste0(
          "Top ",
          nrow(m),
          " CORE+SHARED_2 features by maximum absolute loading across both components."
        ),
      x = paste0("PLS loading component ", comp_x),
      y = paste0("PLS loading component ", comp_y),
      colour = "Omic"
    ) +
    theme_panel(9)

  fam_dir <- file.path(
    DIR_LOADINGS,
    family
  )

  ggplot2::ggsave(
    file.path(
      fam_dir,
      paste0(
        family,
        "_FINAL_EXTENDED_LOADING_MAP_COMP",
        comp_x,
        "_",
        comp_y,
        ".png"
      )
    ),
    p_map,
    width = 8,
    height = 6,
    dpi = DPI
  )

  invisible(TRUE)
}

for (family in FAMILIES) {
  comps_here <- sort(
    unique(
      final_loadings_long$component[
        final_loadings_long$family == family
      ]
    )
  )

  comps_here <- comps_here[
    is.finite(comps_here)
  ]

  if (all(c(1L, 2L) %in% comps_here)) {
    make_loading_map(family, 1L, 2L)
  }

  if (all(c(2L, 3L) %in% comps_here)) {
    make_loading_map(family, 2L, 3L)
  }

  if (all(c(1L, 3L) %in% comps_here)) {
    make_loading_map(family, 1L, 3L)
  }
}

# ==============================================================================
# 15. RUN MANIFEST + README
# ==============================================================================

manifest <- data.frame(
  item = c(
    "BIO_PROGRAM_RUN",
    "PANEL17_RUN",
    "SGL_RUN",
    "OUTPUT_ROOT",
    "CORE_RULE",
    "EXTENDED_RULE",
    "VIP_ROLE",
    "LOADINGS_ROLE",
    "SGL_ROLE",
    "TEST_ROLE"
  ),
  value = c(
    BIO_RUN,
    PANEL17,
    ifelse(is.na(SGL), "NOT_AVAILABLE", SGL),
    OUT_ROOT,
    "Strict Evidence-A feature present in P1 + P2 + P3.",
    "Strict Evidence-A feature present in at least two of P1/P2/P3.",
    paste0(
      "Descriptive annotation only; VIP>=",
      VIP_THRESHOLD,
      " does not define CORE/SHARED_2/SPECIFIC."
    ),
    paste0(
      "Interpretation only. Signs are retained within family/component; ",
      "no cross-analysis signed mean is computed because latent signs may flip."
    ),
    "Corroborative annotation only; Sparse Group Lasso does not filter the consensus panel.",
    "No TEST-based selection is performed in this script."
  ),
  stringsAsFactors = FALSE
)

write_csv(
  manifest,
  file.path(
    OUT_ROOT,
    "RUN_MANIFEST.csv"
  )
)

readme <- c(
  "METATEST — FINAL CONSENSUS PANEL",
  "",
  "PRIMARY PANEL DEFINITIONS",
  "  CORE      = strict Evidence-A feature present in P1 + P2 + P3.",
  "  SHARED_2  = strict Evidence-A feature present in exactly two analyses.",
  "  SPECIFIC  = strict Evidence-A feature present in exactly one analysis.",
  "",
  "  CORE PANEL     = CORE only.",
  "  EXTENDED PANEL = CORE + SHARED_2.",
  "",
  "Neither VIP nor Sparse Group Lasso redefines those panel memberships.",
  "",
  "MAIN TABLES",
  "  tables/00_FINAL_PANEL_MASTER.csv",
  "  tables/01_FINAL_PANEL_CORE_3OF3.csv",
  "  tables/02_FINAL_PANEL_EXTENDED_GE2OF3.csv",
  "  tables/03_CONTEXT_SPECIFIC_1OF3.csv",
  "  tables/04_CORE_HIGH_CONTRIBUTION_VIP_ALL3_DESCRIPTIVE.csv",
  "  tables/05_FINAL_EXTENDED_PANEL_LOADINGS_LONG.csv",
  "  tables/06_FINAL_PANEL_COUNTS.csv",
  "  tables/07_FINAL_PANEL_COUNTS_BY_OMIC.csv",
  "",
  "MAIN FIGURES",
  "  figures/01_FINAL_PANEL_CONSENSUS_HEATMAP.png",
  "  figures/02_FINAL_PANEL_COMPOSITION_BY_OMIC.png",
  "  figures/03_FINAL_PANEL_VIP_RECURRENCE_BY_OMIC.png",
  "  figures/LOADINGS_FINAL_EXTENDED/<P1|P2|P3>/...",
  "",
  "INTERPRETATION",
  "  The consensus tier reflects reproducibility across the three scientific analyses.",
  "  VIP and loadings describe how the retained features contribute to each PLS-DA model.",
  "  Sparse Group Lasso is recorded only as an independent competitive corroboration.",
  "",
  "IMPORTANT",
  "  PLS loading signs are not averaged across P1/P2/P3 because independently fitted",
  "  latent components can be reflected without changing the model."
)

writeLines(
  readme,
  file.path(
    OUT_ROOT,
    "READ_ME_FIRST.txt"
  ),
  useBytes = TRUE
)

writeLines(
  paste(
    "COMPLETE",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  file.path(
    OUT_ROOT,
    "RUN_COMPLETE.flag"
  )
)

latest_root <- file.path(
  BIO_RUN,
  "FINAL_CONSENSUS_PANEL"
)

dir.create(
  latest_root,
  recursive = TRUE,
  showWarnings = FALSE
)

writeLines(
  OUT_ROOT,
  file.path(
    latest_root,
    "LATEST_RUN.txt"
  )
)

message("")
message("============================================================")
message("FINAL CONSENSUS PANEL COMPLETE")
message("CORE (3/3): ", nrow(core_panel))
message("EXTENDED (>=2/3): ", nrow(extended_panel))
message("SPECIFIC (1/3): ", nrow(specific_panel))
message(
  "CORE + VIP>=",
  VIP_THRESHOLD,
  " in all 3 (descriptive): ",
  nrow(core_vip_all3)
)
message("Output: ", OUT_ROOT)
message("============================================================")

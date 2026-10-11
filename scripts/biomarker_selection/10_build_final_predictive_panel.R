#!/usr/bin/env Rscript

# ==============================================================================
# METATEST — 10_build_final_predictive_panel.R
#
# FINAL STAGE
#
# Objective:
#   1. Build a biologically anchored feature shortlist from the frozen outputs
#      of 09_build_biological_synthesis_figures.R.
#   2. Fit a sparse multinomial Elastic Net panel for NP / FD / SO using TRAIN
#      only.
#   3. Tune alpha/lambda exclusively inside TRAIN using repeated stratified
#      3-fold CV.
#   4. Audit selection stability across the repeated TRAIN resamples.
#   5. Fit one final model on all TRAIN.
#   6. Evaluate TEST exactly once. TEST never selects candidates, alpha, lambda
#      or panel features.
#   7. Keep P1/P2/P3 as distinct analyses and report a cross-analysis consensus
#      table without declaring one analysis the winner.
#
# IMPORTANT:
#   - No minimum number of features per omic is imposed.
#   - No quota by TX / PR / ME / CL is imposed.
#   - Q7 ALL / TEST separation metrics are annotation only; they do not enter
#     candidate eligibility or model tuning.
#   - TX eligibility = frozen TX drivers.
#   - PR / ME / CL eligibility = frozen cross-omic links already defined
#     upstream for the frozen biological pathways.
#   - The model itself decides the final sparse panel inside TRAIN.
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

SELECTION_RUN <- Sys.getenv("SELECTION_RUN", unset = "")
BIO_PROGRAM_RUN <- Sys.getenv("BIO_PROGRAM_RUN", unset = "")
BIO_SYNTHESIS_RUN <- Sys.getenv("BIO_SYNTHESIS_RUN", unset = "")

FAMILIES <- trimws(
  strsplit(
    Sys.getenv("PANEL_FAMILIES", unset = "P1,P2,P3"),
    ",",
    fixed = TRUE
  )[[1]]
)
FAMILIES <- FAMILIES[nzchar(FAMILIES)]

PANEL_SEED <- as.integer(
  Sys.getenv("PANEL_SEED", unset = "20261001")
)

PANEL_CV_REPEATS <- as.integer(
  Sys.getenv("PANEL_CV_REPEATS", unset = "30")
)

PANEL_CV_FOLDS <- as.integer(
  Sys.getenv("PANEL_CV_FOLDS", unset = "3")
)

PANEL_ALPHA_GRID <- as.numeric(
  trimws(
    strsplit(
      Sys.getenv(
        "PANEL_ALPHA_GRID",
        unset = "0.25,0.50,0.75,1.00"
      ),
      ",",
      fixed = TRUE
    )[[1]]
  )
)

PANEL_USE_CLASS_WEIGHTS <- tolower(
  Sys.getenv(
    "PANEL_USE_CLASS_WEIGHTS",
    unset = "true"
  )
) %in% c("true", "1", "yes", "y")

PANEL_DPI <- as.integer(
  Sys.getenv("PANEL_DPI", unset = "180")
)

PANEL_OUTPUT_ROOT <- Sys.getenv(
  "PANEL_OUTPUT_ROOT",
  unset = ""
)

CLASS_LEVELS <- c("NP", "FD", "SO")
VIEW_LEVELS <- c("tx", "pr", "me", "cl")

if (!is.finite(PANEL_SEED)) {
  stop("PANEL_SEED must be finite.", call. = FALSE)
}

if (
  !is.finite(PANEL_CV_REPEATS) ||
  PANEL_CV_REPEATS < 2L
) {
  stop("PANEL_CV_REPEATS must be >= 2.", call. = FALSE)
}

if (
  !is.finite(PANEL_CV_FOLDS) ||
  PANEL_CV_FOLDS < 2L
) {
  stop("PANEL_CV_FOLDS must be >= 2.", call. = FALSE)
}

if (
  !length(PANEL_ALPHA_GRID) ||
  any(!is.finite(PANEL_ALPHA_GRID)) ||
  any(PANEL_ALPHA_GRID <= 0) ||
  any(PANEL_ALPHA_GRID > 1)
) {
  stop(
    "PANEL_ALPHA_GRID must contain values in (0, 1].",
    call. = FALSE
  )
}

# ==============================================================================
# 1. DEPENDENCIES
# ==============================================================================

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop(
    paste0(
      "Package 'glmnet' is required. Install it with:\n",
      "install.packages('glmnet')"
    ),
    call. = FALSE
  )
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    paste0(
      "Package 'ggplot2' is required. Install it with:\n",
      "install.packages('ggplot2')"
    ),
    call. = FALSE
  )
}

# ==============================================================================
# 2. GENERIC HELPERS
# ==============================================================================

assert <- function(ok, msg) {
  if (!isTRUE(ok)) {
    stop(msg, call. = FALSE)
  }
}

canonical_path <- function(p, must_work = TRUE) {
  normalizePath(
    path.expand(p),
    winslash = "/",
    mustWork = must_work
  )
}

read_csv_checked <- function(
    path,
    required = character()
) {
  assert(
    file.exists(path),
    paste("Missing required file:", path)
  )

  x <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  missing <- setdiff(
    required,
    names(x)
  )

  assert(
    !length(missing),
    paste(
      "Missing columns in",
      basename(path),
      ":",
      paste(missing, collapse = ", ")
    )
  )

  x
}

write_csv <- function(
    x,
    path
) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

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

as_flag <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  toupper(
    trimws(
      as.character(x)
    )
  ) %in% c(
    "TRUE",
    "T",
    "1",
    "YES",
    "Y"
  )
}

collapse_unique <- function(
    x,
    sep = ";"
) {
  x <- as.character(x)

  x <- x[
    !is.na(x) &
      nzchar(x) &
      x != "NA"
  ]

  x <- sort(
    unique(x)
  )

  if (!length(x)) {
    return(NA_character_)
  }

  paste(
    x,
    collapse = sep
  )
}

max_or_na <- function(x) {
  x <- suppressWarnings(
    as.numeric(x)
  )

  x <- x[
    is.finite(x)
  ]

  if (!length(x)) {
    return(NA_real_)
  }

  max(x)
}

mean_or_na <- function(x) {
  x <- suppressWarnings(
    as.numeric(x)
  )

  x <- x[
    is.finite(x)
  ]

  if (!length(x)) {
    return(NA_real_)
  }

  mean(x)
}

first_nonempty <- function(x) {
  x <- as.character(x)

  x <- x[
    !is.na(x) &
      nzchar(x) &
      x != "NA"
  ]

  if (!length(x)) {
    return(NA_character_)
  }

  x[1L]
}

safe_filename <- function(x) {
  x <- gsub(
    "[^A-Za-z0-9_.-]+",
    "_",
    as.character(x)
  )

  gsub(
    "_+",
    "_",
    x
  )
}

plot_theme <- function(base_size = 11) {
  ggplot2::theme_bw(
    base_size = base_size
  ) +
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
}

# ==============================================================================
# 3. LOCATE FROZEN RUNS
# ==============================================================================

find_latest_selection_run <- function(upstream) {
  base <- file.path(
    upstream,
    "downstream",
    "seleccion_mejor_escenario_train"
  )

  assert(
    dir.exists(base),
    paste(
      "Selection root does not exist:",
      base
    )
  )

  dirs <- list.dirs(
    base,
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

  assert(
    length(dirs) > 0L,
    paste(
      "No complete selection runs found in",
      base
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

find_latest_bio_program_run <- function(
    selection_run
) {
  root <- file.path(
    selection_run,
    "BIOLOGICAL_PROGRAMS"
  )

  assert(
    dir.exists(root),
    paste(
      "BIOLOGICAL_PROGRAMS does not exist:",
      root
    )
  )

  latest <- file.path(
    root,
    "LATEST_RUN.txt"
  )

  if (file.exists(latest)) {
    p <- trimws(
      readLines(
        latest,
        warn = FALSE
      )[1L]
    )

    if (
      nzchar(p) &&
      dir.exists(p) &&
      file.exists(
        file.path(
          p,
          "BIOLOGICAL_PROGRAMS_OBJECT.rds"
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

  dirs <- dirs[
    file.exists(
      file.path(
        dirs,
        "BIOLOGICAL_PROGRAMS_OBJECT.rds"
      )
    )
  ]

  assert(
    length(dirs) > 0L,
    paste(
      "No complete biological-program runs in",
      root
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

resolve_bio_program_run <- function() {
  if (nzchar(BIO_PROGRAM_RUN)) {
    return(
      canonical_path(
        BIO_PROGRAM_RUN,
        TRUE
      )
    )
  }

  selection_run <- if (
    nzchar(SELECTION_RUN)
  ) {
    canonical_path(
      SELECTION_RUN,
      TRUE
    )
  } else {
    find_latest_selection_run(
      UPSTREAM_DIR
    )
  }

  find_latest_bio_program_run(
    selection_run
  )
}

find_latest_synthesis_run <- function(
    bio_program_run
) {
  if (nzchar(BIO_SYNTHESIS_RUN)) {
    out <- canonical_path(
      BIO_SYNTHESIS_RUN,
      TRUE
    )

    return(out)
  }

  root <- file.path(
    bio_program_run,
    "BIOLOGICAL_SYNTHESIS_FIGURES"
  )

  assert(
    dir.exists(root),
    paste(
      "BIOLOGICAL_SYNTHESIS_FIGURES does not exist:",
      root
    )
  )

  latest <- file.path(
    root,
    "LATEST_RUN.txt"
  )

  required_rel <- file.path(
    "plot_data",
    "FIG7_INDIVIDUAL_MANIFESTATION",
    "00_INDEX",
    "Q7_PATHWAY_PHENOTYPE_SEPARATION_MASTER.csv"
  )

  if (file.exists(latest)) {
    p <- trimws(
      readLines(
        latest,
        warn = FALSE
      )[1L]
    )

    if (
      nzchar(p) &&
      dir.exists(p) &&
      file.exists(
        file.path(
          p,
          required_rel
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

  dirs <- dirs[
    file.exists(
      file.path(
        dirs,
        required_rel
      )
    )
  ]

  assert(
    length(dirs) > 0L,
    paste(
      "No completed biological synthesis run with Q7 found in",
      root
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

BIO_RUN <- resolve_bio_program_run()

SYNTHESIS_RUN <- find_latest_synthesis_run(
  BIO_RUN
)

if (nzchar(PANEL_OUTPUT_ROOT)) {
  OUT_ROOT <- path.expand(
    PANEL_OUTPUT_ROOT
  )
} else {
  stamp <- format(
    Sys.time(),
    "%Y%m%dT%H%M%S"
  )

  OUT_ROOT <- file.path(
    BIO_RUN,
    "FINAL_PREDICTIVE_PANEL",
    paste0(
      stamp,
      "_pid",
      Sys.getpid()
    )
  )
}

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

OUT_ROOT <- canonical_path(
  OUT_ROOT,
  TRUE
)

DIR_TABLES <- file.path(
  OUT_ROOT,
  "tables"
)

DIR_FIGURES <- file.path(
  OUT_ROOT,
  "figures"
)

DIR_MODELS <- file.path(
  OUT_ROOT,
  "models"
)

dir.create(
  DIR_TABLES,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  DIR_FIGURES,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  DIR_MODELS,
  recursive = TRUE,
  showWarnings = FALSE
)

message("============================================================")
message("METATEST — FINAL BIOLOGICAL SHORTLIST + PREDICTIVE PANEL")
message("BIOLOGICAL_PROGRAMS run: ", BIO_RUN)
message("BIOLOGICAL_SYNTHESIS run: ", SYNTHESIS_RUN)
message("Output: ", OUT_ROOT)
message("============================================================")

# ==============================================================================
# 4. LOAD FROZEN BIOLOGICAL OUTPUTS
# ==============================================================================

q5_catalog_path <- file.path(
  SYNTHESIS_RUN,
  "plot_data",
  "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
  "00_INDEX",
  "FIG5_MASTER_PATHWAY_FEATURE_CATALOG.csv"
)

q6_pathway_path <- file.path(
  SYNTHESIS_RUN,
  "plot_data",
  "FIG6_COMMON_SPECIFIC_REINFORCED",
  "00_INDEX",
  "Q6_PATHWAY_MASTER.csv"
)

q7_separation_path <- file.path(
  SYNTHESIS_RUN,
  "plot_data",
  "FIG7_INDIVIDUAL_MANIFESTATION",
  "00_INDEX",
  "Q7_PATHWAY_PHENOTYPE_SEPARATION_MASTER.csv"
)

catalog <- read_csv_checked(
  q5_catalog_path,
  required = c(
    "family",
    "program_id",
    "program_short",
    "library",
    "pathway",
    "pathway_label",
    "phenotypic_side",
    "factor_side",
    "view_code",
    "feature_original",
    "feature_label",
    "is_tx_driver",
    "tx_driver_weight",
    "pathway_feature_cosine",
    "available_in_pooled"
  )
)

q6 <- read_csv_checked(
  q6_pathway_path,
  required = c(
    "family",
    "program_id",
    "program_short",
    "library",
    "pathway",
    "phenotypic_side",
    "pathway_analysis_scope",
    "pathway_phenotype_scope",
    "multimodal_reinforcement_status",
    "n_nonTX_views_linked"
  )
)

q7 <- read_csv_checked(
  q7_separation_path,
  required = c(
    "family",
    "program_id",
    "library",
    "pathway",
    "phenotypic_side",
    "eta2_TRAIN",
    "eta2_ALL"
  )
)

catalog$family <- as.character(
  catalog$family
)

catalog$view_code <- tolower(
  as.character(
    catalog$view_code
  )
)

catalog$feature_original <- as.character(
  catalog$feature_original
)

catalog$feature_label <- as.character(
  catalog$feature_label
)

catalog$is_tx_driver <- as_flag(
  catalog$is_tx_driver
)

catalog$available_in_pooled <- as_flag(
  catalog$available_in_pooled
)

catalog$tx_driver_weight <- suppressWarnings(
  as.numeric(
    catalog$tx_driver_weight
  )
)

catalog$pathway_feature_cosine <- suppressWarnings(
  as.numeric(
    catalog$pathway_feature_cosine
  )
)

catalog$feature_key <- paste(
  catalog$view_code,
  catalog$feature_original,
  sep = "\r"
)

catalog$predictor_key <- paste(
  toupper(
    catalog$view_code
  ),
  catalog$feature_original,
  sep = "::"
)

# ==============================================================================
# 5. BUILD BIOLOGICAL CANDIDATE UNIVERSE
# ==============================================================================

# ------------------------------------------------------------------------------
# Q6 annotations.
# ------------------------------------------------------------------------------

q6_keep <- intersect(
  c(
    "family",
    "program_id",
    "program_short",
    "library",
    "pathway",
    "phenotypic_side",
    "analyses_with_pathway",
    "n_analyses_with_pathway",
    "pathway_analysis_scope",
    "phenotypic_sides_with_pathway",
    "n_phenotypic_sides_with_pathway",
    "pathway_phenotype_scope",
    "program_analysis_scope",
    "program_phenotype_scope",
    "pathway_program_joint_status",
    "strict_threeway_pathway_in_strict_program",
    "multimodal_reinforcement_status",
    "n_nonTX_views_linked"
  ),
  names(q6)
)

q6_ann <- unique(
  q6[
    ,
    q6_keep,
    drop = FALSE
  ]
)

merge_key <- c(
  "family",
  "program_id",
  "library",
  "pathway",
  "phenotypic_side"
)

catalog_q6 <- merge(
  catalog,
  q6_ann,
  by = merge_key,
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("", "_q6")
)

# ------------------------------------------------------------------------------
# Q7 annotations.
#
# IMPORTANT:
#   These are carried for interpretation only.
#   They do NOT determine candidate eligibility.
# ------------------------------------------------------------------------------

q7_keep <- intersect(
  c(
    "family",
    "program_id",
    "library",
    "pathway",
    "phenotypic_side",
    "eta2_TRAIN",
    "eta2_TEST",
    "eta2_ALL",
    "g_NP_FD_TRAIN",
    "g_NP_SO_TRAIN",
    "g_FD_SO_TRAIN",
    "g_NP_FD_ALL",
    "g_NP_SO_ALL",
    "g_FD_SO_ALL",
    "permanova_R2_TRAIN",
    "permanova_R2_TEST",
    "permanova_R2_ALL",
    "separation_index_ALL",
    "relative_separation_class_ALL"
  ),
  names(q7)
)

q7_ann <- unique(
  q7[
    ,
    q7_keep,
    drop = FALSE
  ]
)

candidate_occurrences <- merge(
  catalog_q6,
  q7_ann,
  by = merge_key,
  all.x = TRUE,
  sort = FALSE
)

# ------------------------------------------------------------------------------
# Frozen biological eligibility.
#
# TX:
#   only frozen TX drivers.
#
# PR / ME / CL:
#   frozen cross-omic links with a finite non-zero cosine.
#
# No outcome-derived ALL / TEST threshold is used.
# ------------------------------------------------------------------------------

candidate_occurrences$eligible_tx <-
  candidate_occurrences$view_code == "tx" &
  candidate_occurrences$is_tx_driver

candidate_occurrences$eligible_non_tx <-
  candidate_occurrences$view_code %in%
    c("pr", "me", "cl") &
  is.finite(
    candidate_occurrences$pathway_feature_cosine
  ) &
  abs(
    candidate_occurrences$pathway_feature_cosine
  ) > 0

candidate_occurrences$eligible_biological <-
  candidate_occurrences$available_in_pooled &
  (
    candidate_occurrences$eligible_tx |
      candidate_occurrences$eligible_non_tx
  )

candidate_occurrences$candidate_reason <- ifelse(
  candidate_occurrences$eligible_tx,
  "TX_DRIVER_FROZEN",
  ifelse(
    candidate_occurrences$eligible_non_tx,
    "NON_TX_FROZEN_CROSSOMIC_LINK",
    NA_character_
  )
)

candidate_occurrences <- candidate_occurrences[
  candidate_occurrences$eligible_biological &
    candidate_occurrences$family %in%
      FAMILIES &
    candidate_occurrences$view_code %in%
      VIEW_LEVELS &
    !is.na(
      candidate_occurrences$feature_original
    ) &
    nzchar(
      candidate_occurrences$feature_original
    ),
  ,
  drop = FALSE
]

assert(
  nrow(candidate_occurrences) > 0L,
  "No biologically eligible candidate features were found."
)

write_csv(
  candidate_occurrences,
  file.path(
    DIR_TABLES,
    "01_BIOLOGICAL_CANDIDATE_OCCURRENCES.csv"
  )
)

# ==============================================================================
# 6. AGGREGATE BIOLOGICAL SHORTLIST TO FEATURE LEVEL
# ==============================================================================

split_candidate <- split(
  candidate_occurrences,
  paste(
    candidate_occurrences$family,
    candidate_occurrences$feature_key,
    sep = "\r\r"
  )
)

shortlist_rows <- lapply(
  split_candidate,
  function(z) {
    data.frame(
      family =
        first_nonempty(
          z$family
        ),

      feature_key =
        first_nonempty(
          z$feature_key
        ),

      predictor_key =
        first_nonempty(
          z$predictor_key
        ),

      view_code =
        first_nonempty(
          z$view_code
        ),

      feature_original =
        first_nonempty(
          z$feature_original
        ),

      feature_label =
        first_nonempty(
          z$feature_label
        ),

      candidate_reason =
        collapse_unique(
          z$candidate_reason
        ),

      n_programs =
        length(
          unique(
            z$program_id
          )
        ),

      n_pathways =
        length(
          unique(
            paste(
              z$library,
              z$pathway,
              sep = "::"
            )
          )
        ),

      programs =
        collapse_unique(
          z$program_id
        ),

      pathways =
        collapse_unique(
          paste(
            z$library,
            z$pathway,
            sep = "::"
          )
        ),

      phenotypic_sides =
        collapse_unique(
          z$phenotypic_side
        ),

      pathway_analysis_scopes =
        if (
          "pathway_analysis_scope" %in%
            names(z)
        ) {
          collapse_unique(
            z$pathway_analysis_scope
          )
        } else {
          NA_character_
        },

      pathway_phenotype_scopes =
        if (
          "pathway_phenotype_scope" %in%
            names(z)
        ) {
          collapse_unique(
            z$pathway_phenotype_scope
          )
        } else {
          NA_character_
        },

      pathway_program_joint_status =
        if (
          "pathway_program_joint_status" %in%
            names(z)
        ) {
          collapse_unique(
            z$pathway_program_joint_status
          )
        } else {
          NA_character_
        },

      multimodal_reinforcement_status =
        if (
          "multimodal_reinforcement_status" %in%
            names(z)
        ) {
          collapse_unique(
            z$multimodal_reinforcement_status
          )
        } else {
          NA_character_
        },

      max_tx_driver_weight =
        max_or_na(
          z$tx_driver_weight
        ),

      max_abs_pathway_feature_cosine =
        max_or_na(
          abs(
            z$pathway_feature_cosine
          )
        ),

      max_eta2_TRAIN =
        if (
          "eta2_TRAIN" %in%
            names(z)
        ) {
          max_or_na(
            z$eta2_TRAIN
          )
        } else {
          NA_real_
        },

      max_eta2_ALL_descriptive =
        if (
          "eta2_ALL" %in%
            names(z)
        ) {
          max_or_na(
            z$eta2_ALL
          )
        } else {
          NA_real_
        },

      max_permanova_R2_TRAIN =
        if (
          "permanova_R2_TRAIN" %in%
            names(z)
        ) {
          max_or_na(
            z$permanova_R2_TRAIN
          )
        } else {
          NA_real_
        },

      max_permanova_R2_ALL_descriptive =
        if (
          "permanova_R2_ALL" %in%
            names(z)
        ) {
          max_or_na(
            z$permanova_R2_ALL
          )
        } else {
          NA_real_
        },

      stringsAsFactors = FALSE
    )
  }
)

biological_shortlist <- do.call(
  rbind,
  shortlist_rows
)

rownames(
  biological_shortlist
) <- NULL

biological_shortlist <- biological_shortlist[
  order(
    match(
      biological_shortlist$family,
      FAMILIES
    ),
    match(
      biological_shortlist$view_code,
      VIEW_LEVELS
    ),
    biological_shortlist$feature_label
  ),
  ,
  drop = FALSE
]

write_csv(
  biological_shortlist,
  file.path(
    DIR_TABLES,
    "02_BIOLOGICAL_FEATURE_SHORTLIST.csv"
  )
)

# Candidate-count figure.
candidate_counts <- stats::aggregate(
  predictor_key ~ family + view_code,
  data = biological_shortlist,
  FUN = length
)

names(candidate_counts)[
  names(candidate_counts) ==
    "predictor_key"
] <- "n_candidates"

p_candidates <- ggplot2::ggplot(
  candidate_counts,
  ggplot2::aes(
    x = view_code,
    y = n_candidates
  )
) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(
    ~ family,
    scales = "free_y"
  ) +
  ggplot2::labs(
    title =
      "Biological candidate universe entering predictive modeling",
    subtitle =
      "TX = frozen drivers; PR/ME/CL = frozen cross-omic links. No omic quota.",
    x = "View",
    y = "Candidate features"
  ) +
  plot_theme(11)

ggplot2::ggsave(
  file.path(
    DIR_FIGURES,
    "01_BIOLOGICAL_CANDIDATE_COUNTS_BY_VIEW.png"
  ),
  p_candidates,
  width = 9,
  height = 5.5,
  dpi = PANEL_DPI
)

# ==============================================================================
# 7. LOAD EXACT READY-FOR-MODELING MATRICES PER FAMILY
# ==============================================================================

extract_group_table <- function(g) {
  if (is.data.frame(g)) {
    assert(
      all(
        c(
          "id",
          "grupo"
        ) %in%
          names(g)
      ),
      "grupo_train/grupo_test must contain id and grupo."
    )

    return(
      data.frame(
        sample_id =
          as.character(
            g$id
          ),
        grupo =
          as.character(
            g$grupo
          ),
        stringsAsFactors = FALSE
      )
    )
  }

  if (
    is.atomic(g) &&
    !is.null(
      names(g)
    )
  ) {
    return(
      data.frame(
        sample_id =
          names(g),
        grupo =
          as.character(g),
        stringsAsFactors = FALSE
      )
    )
  }

  stop(
    "Unrecognized group structure in ready_for_modeling.",
    call. = FALSE
  )
}

load_family_ready_data <- function(
    bio_program_run,
    family
) {
  family_dir <- file.path(
    bio_program_run,
    family
  )

  audit <- read_csv_checked(
    file.path(
      family_dir,
      "00_PROGRAM_DEFINITION_AUDIT.csv"
    ),
    required = c(
      "family",
      "observed_train_run"
    )
  )

  observed_train_run <- as.character(
    audit$observed_train_run[1L]
  )

  input_audit <- read_csv_checked(
    file.path(
      observed_train_run,
      "tables",
      "00_INPUT_AUDIT.csv"
    ),
    required =
      "ready_for_modeling"
  )

  ready_for_modeling <- as.character(
    input_audit$ready_for_modeling[1L]
  )

  assert(
    file.exists(
      ready_for_modeling
    ),
    paste(
      "ready_for_modeling does not exist:",
      ready_for_modeling
    )
  )

  rfm <- readRDS(
    ready_for_modeling
  )

  assert(
    all(
      c(
        VIEW_LEVELS,
        "grupo_train",
        "grupo_test"
      ) %in%
        names(rfm)
    ),
    paste(
      "ready_for_modeling is missing expected elements for",
      family
    )
  )

  train_meta <- extract_group_table(
    rfm$grupo_train
  )

  test_meta <- extract_group_table(
    rfm$grupo_test
  )

  train_meta$split <- "TRAIN"
  test_meta$split <- "TEST"

  train_meta$grupo <- factor(
    train_meta$grupo,
    levels = CLASS_LEVELS
  )

  test_meta$grupo <- factor(
    test_meta$grupo,
    levels = CLASS_LEVELS
  )

  assert(
    !anyNA(
      train_meta$grupo
    ),
    paste(
      "Unknown TRAIN group label in",
      family
    )
  )

  assert(
    !anyNA(
      test_meta$grupo
    ),
    paste(
      "Unknown TEST group label in",
      family
    )
  )

  matrices_train <- list()
  matrices_test <- list()

  for (v in VIEW_LEVELS) {
    assert(
      all(
        c(
          "train",
          "test"
        ) %in%
          names(
            rfm[[v]]
          )
      ),
      paste(
        family,
        v,
        "does not contain train/test matrices."
      )
    )

    tr <- as.matrix(
      rfm[[v]]$train
    )

    te <- as.matrix(
      rfm[[v]]$test
    )

    assert(
      !is.null(
        rownames(tr)
      ) &&
        !is.null(
          colnames(tr)
        ),
      paste(
        "Missing TRAIN dimnames:",
        family,
        v
      )
    )

    assert(
      !is.null(
        rownames(te)
      ) &&
        !is.null(
          colnames(te)
        ),
      paste(
        "Missing TEST dimnames:",
        family,
        v
      )
    )

    assert(
      all(
        train_meta$sample_id %in%
          colnames(tr)
      ),
      paste(
        "TRAIN sample IDs do not align:",
        family,
        v
      )
    )

    assert(
      all(
        test_meta$sample_id %in%
          colnames(te)
      ),
      paste(
        "TEST sample IDs do not align:",
        family,
        v
      )
    )

    matrices_train[[v]] <-
      tr[
        ,
        train_meta$sample_id,
        drop = FALSE
      ]

    matrices_test[[v]] <-
      te[
        ,
        test_meta$sample_id,
        drop = FALSE
      ]
  }

  list(
    family =
      family,
    ready_for_modeling =
      canonical_path(
        ready_for_modeling,
        TRUE
      ),
    train_meta =
      train_meta,
    test_meta =
      test_meta,
    train =
      matrices_train,
    test =
      matrices_test
  )
}

family_data <- lapply(
  FAMILIES,
  function(family) {
    message("[INPUT] ", family)

    load_family_ready_data(
      BIO_RUN,
      family
    )
  }
)

names(
  family_data
) <- FAMILIES

input_audit_rows <- lapply(
  family_data,
  function(x) {
    data.frame(
      family =
        x$family,
      ready_for_modeling =
        x$ready_for_modeling,
      n_train =
        nrow(
          x$train_meta
        ),
      n_test =
        nrow(
          x$test_meta
        ),
      train_class_counts =
        paste(
          names(
            table(
              x$train_meta$grupo
            )
          ),
          as.integer(
            table(
              x$train_meta$grupo
            )
          ),
          sep = "=",
          collapse = ";"
        ),
      test_class_counts =
        paste(
          names(
            table(
              x$test_meta$grupo
            )
          ),
          as.integer(
            table(
              x$test_meta$grupo
            )
          ),
          sep = "=",
          collapse = ";"
        ),
      stringsAsFactors = FALSE
    )
  }
)

input_audit <- do.call(
  rbind,
  input_audit_rows
)

write_csv(
  input_audit,
  file.path(
    DIR_TABLES,
    "00_INPUT_AUDIT.csv"
  )
)

# ==============================================================================
# 8. BUILD MODEL MATRICES FROM THE BIOLOGICAL SHORTLIST
# ==============================================================================

build_family_model_matrix <- function(
    family,
    shortlist,
    family_input
) {
  s <- shortlist[
    shortlist$family ==
      family,
    ,
    drop = FALSE
  ]

  assert(
    nrow(s) > 0L,
    paste(
      "No shortlist candidates for",
      family
    )
  )

  train_ids <-
    family_input$train_meta$sample_id

  test_ids <-
    family_input$test_meta$sample_id

  train_cols <- list()
  test_cols <- list()
  technical_rows <- list()
  rr <- 1L

  for (
    i in seq_len(
      nrow(s)
    )
  ) {
    v <- as.character(
      s$view_code[i]
    )

    f <- as.character(
      s$feature_original[i]
    )

    pred <- as.character(
      s$predictor_key[i]
    )

    tr_mat <- family_input$train[[v]]
    te_mat <- family_input$test[[v]]

    present_train <-
      f %in%
        rownames(
          tr_mat
        )

    present_test <-
      f %in%
        rownames(
          te_mat
        )

    if (!present_train) {
      technical_rows[[rr]] <- data.frame(
        family = family,
        predictor_key = pred,
        view_code = v,
        feature_original = f,
        eligible_for_model = FALSE,
        reason = "ABSENT_FROM_TRAIN_MATRIX",
        train_sd = NA_real_,
        n_nonfinite_train = NA_integer_,
        n_nonfinite_test = NA_integer_,
        stringsAsFactors = FALSE
      )

      rr <- rr + 1L
      next
    }

    tr <- as.numeric(
      tr_mat[
        f,
        train_ids
      ]
    )

    te <- if (present_test) {
      as.numeric(
        te_mat[
          f,
          test_ids
        ]
      )
    } else {
      rep(
        NA_real_,
        length(test_ids)
      )
    }

    n_bad_train <- sum(
      !is.finite(tr)
    )

    n_bad_test <- sum(
      !is.finite(te)
    )

    train_sd <- if (
      n_bad_train == 0L
    ) {
      stats::sd(tr)
    } else {
      NA_real_
    }

    if (
      n_bad_train > 0L
    ) {
      technical_rows[[rr]] <- data.frame(
        family = family,
        predictor_key = pred,
        view_code = v,
        feature_original = f,
        eligible_for_model = FALSE,
        reason = "NONFINITE_IN_TRAIN",
        train_sd = train_sd,
        n_nonfinite_train = n_bad_train,
        n_nonfinite_test = n_bad_test,
        stringsAsFactors = FALSE
      )

      rr <- rr + 1L
      next
    }

    if (
      !is.finite(train_sd) ||
      train_sd <= 0
    ) {
      technical_rows[[rr]] <- data.frame(
        family = family,
        predictor_key = pred,
        view_code = v,
        feature_original = f,
        eligible_for_model = FALSE,
        reason = "ZERO_VARIANCE_IN_TRAIN",
        train_sd = train_sd,
        n_nonfinite_train = n_bad_train,
        n_nonfinite_test = n_bad_test,
        stringsAsFactors = FALSE
      )

      rr <- rr + 1L
      next
    }

    # TEST missing values, if any, are filled using the TRAIN median only.
    # This does not use TEST labels or TEST distribution.
    train_median <- stats::median(
      tr
    )

    te[
      !is.finite(te)
    ] <- train_median

    train_cols[[pred]] <- tr
    test_cols[[pred]] <- te

    technical_rows[[rr]] <- data.frame(
      family = family,
      predictor_key = pred,
      view_code = v,
      feature_original = f,
      eligible_for_model = TRUE,
      reason = "ELIGIBLE",
      train_sd = train_sd,
      n_nonfinite_train = n_bad_train,
      n_nonfinite_test = n_bad_test,
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  technical_audit <- do.call(
    rbind,
    technical_rows
  )

  assert(
    length(train_cols) > 0L,
    paste(
      "No technically eligible predictors for",
      family
    )
  )

  x_train <- do.call(
    cbind,
    train_cols
  )

  x_test <- do.call(
    cbind,
    test_cols
  )

  x_train <- as.matrix(
    x_train
  )

  x_test <- as.matrix(
    x_test
  )

  rownames(
    x_train
  ) <- train_ids

  rownames(
    x_test
  ) <- test_ids

  y_train <- factor(
    family_input$train_meta$grupo,
    levels = CLASS_LEVELS
  )

  y_test <- factor(
    family_input$test_meta$grupo,
    levels = CLASS_LEVELS
  )

  assert(
    min(
      table(y_train)
    ) >= PANEL_CV_FOLDS,
    paste(
      "PANEL_CV_FOLDS =",
      PANEL_CV_FOLDS,
      "is too large for the smallest TRAIN class in",
      family,
      ". Counts:",
      paste(
        names(
          table(y_train)
        ),
        table(y_train),
        collapse = ";"
      )
    )
  )

  list(
    x_train = x_train,
    y_train = y_train,
    x_test = x_test,
    y_test = y_test,
    train_meta =
      family_input$train_meta,
    test_meta =
      family_input$test_meta,
    technical_audit =
      technical_audit
  )
}

model_data <- lapply(
  FAMILIES,
  function(family) {
    message("[MATRIX] ", family)

    build_family_model_matrix(
      family =
        family,
      shortlist =
        biological_shortlist,
      family_input =
        family_data[[family]]
    )
  }
)

names(
  model_data
) <- FAMILIES

technical_audit_all <- do.call(
  rbind,
  lapply(
    FAMILIES,
    function(family) {
      model_data[[family]]$technical_audit
    }
  )
)

write_csv(
  technical_audit_all,
  file.path(
    DIR_TABLES,
    "03_TECHNICAL_FEATURE_ELIGIBILITY_AUDIT.csv"
  )
)

# ==============================================================================
# 9. MODELING HELPERS
# ==============================================================================

make_class_weights <- function(y) {
  y <- factor(
    y,
    levels = CLASS_LEVELS
  )

  if (!PANEL_USE_CLASS_WEIGHTS) {
    return(
      rep(
        1,
        length(y)
      )
    )
  }

  tab <- table(y)

  w_class <- length(y) /
    (
      length(CLASS_LEVELS) *
        as.numeric(tab)
    )

  names(
    w_class
  ) <- names(tab)

  as.numeric(
    w_class[
      as.character(y)
    ]
  )
}

make_stratified_foldid <- function(
    y,
    k,
    seed
) {
  set.seed(seed)

  y <- factor(
    y,
    levels = CLASS_LEVELS
  )

  foldid <- integer(
    length(y)
  )

  for (
    cls in CLASS_LEVELS
  ) {
    idx <- which(
      y == cls
    )

    idx <- sample(
      idx,
      length(idx),
      replace = FALSE
    )

    foldid[idx] <- rep(
      seq_len(k),
      length.out =
        length(idx)
    )
  }

  assert(
    all(
      foldid >= 1L &
        foldid <= k
    ),
    "Could not create stratified fold assignment."
  )

  foldid
}

cv_one_alpha <- function(
    x,
    y,
    alpha,
    foldid,
    weights
) {
  fit <- glmnet::cv.glmnet(
    x = x,
    y = y,
    family = "multinomial",
    type.multinomial = "grouped",
    alpha = alpha,
    foldid = foldid,
    nfolds = max(foldid),
    type.measure = "deviance",
    standardize = TRUE,
    intercept = TRUE,
    weights = weights,
    keep = FALSE,
    parallel = FALSE
  )

  idx_1se <- which.min(
    abs(
      fit$lambda -
        fit$lambda.1se
    )
  )

  idx_min <- which.min(
    abs(
      fit$lambda -
        fit$lambda.min
    )
  )

  data.frame(
    alpha = alpha,
    lambda_1se = fit$lambda.1se,
    lambda_min = fit$lambda.min,
    cvm_1se = fit$cvm[idx_1se],
    cvm_min = fit$cvm[idx_min],
    stringsAsFactors = FALSE
  )
}

tune_elastic_net <- function(
    x,
    y,
    family,
    repeats,
    alpha_grid,
    seed
) {
  rows <- list()
  rr <- 1L

  weights <- make_class_weights(
    y
  )

  for (
    r in seq_len(
      repeats
    )
  ) {
    foldid <- make_stratified_foldid(
      y = y,
      k = PANEL_CV_FOLDS,
      seed =
        seed +
        1000L * r
    )

    for (
      a in alpha_grid
    ) {
      ans <- try(
        cv_one_alpha(
          x = x,
          y = y,
          alpha = a,
          foldid = foldid,
          weights = weights
        ),
        silent = TRUE
      )

      if (
        inherits(
          ans,
          "try-error"
        )
      ) {
        warning(
          "Tuning failed for ",
          family,
          " | repeat ",
          r,
          " | alpha ",
          a,
          ": ",
          as.character(ans)
        )

        next
      }

      ans$family <- family
      ans$repeat_id <- r

      rows[[rr]] <- ans
      rr <- rr + 1L
    }
  }

  assert(
    length(rows) > 0L,
    paste(
      "No successful tuning fits for",
      family
    )
  )

  tuning <- do.call(
    rbind,
    rows
  )

  alpha_summary <- do.call(
    rbind,
    lapply(
      sort(
        unique(
          tuning$alpha
        )
      ),
      function(a) {
        z <- tuning[
          tuning$alpha == a,
          ,
          drop = FALSE
        ]

        data.frame(
          family = family,
          alpha = a,
          n_successful_repeats =
            nrow(z),
          mean_cvm_1se =
            mean_or_na(
              z$cvm_1se
            ),
          sd_cvm_1se =
            stats::sd(
              z$cvm_1se,
              na.rm = TRUE
            ),
          mean_cvm_min =
            mean_or_na(
              z$cvm_min
            ),
          median_lambda_1se =
            stats::median(
              z$lambda_1se,
              na.rm = TRUE
            ),
          median_lambda_min =
            stats::median(
              z$lambda_min,
              na.rm = TRUE
            ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  alpha_summary <- alpha_summary[
    order(
      alpha_summary$mean_cvm_1se,
      alpha_summary$alpha
    ),
    ,
    drop = FALSE
  ]

  best_alpha <- alpha_summary$alpha[1L]

  best_rows <- tuning[
    tuning$alpha ==
      best_alpha,
    ,
    drop = FALSE
  ]

  list(
    tuning =
      tuning,
    alpha_summary =
      alpha_summary,
    alpha =
      best_alpha,
    lambda_1se =
      stats::median(
        best_rows$lambda_1se,
        na.rm = TRUE
      ),
    lambda_min =
      stats::median(
        best_rows$lambda_min,
        na.rm = TRUE
      )
  )
}

extract_multinomial_coefficients <- function(
    fit,
    lambda
) {
  cf <- stats::coef(
    fit,
    s = lambda
  )

  rows <- list()
  rr <- 1L

  for (
    cls in names(cf)
  ) {
    m <- as.matrix(
      cf[[cls]]
    )

    rows[[rr]] <- data.frame(
      class =
        cls,
      predictor_key =
        rownames(m),
      coefficient =
        as.numeric(
          m[, 1L]
        ),
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  out <- do.call(
    rbind,
    rows
  )

  out
}

selected_predictors_from_fit <- function(
    fit,
    lambda,
    tol = 1e-12
) {
  cf <- extract_multinomial_coefficients(
    fit,
    lambda
  )

  cf <- cf[
    cf$predictor_key !=
      "(Intercept)",
    ,
    drop = FALSE
  ]

  selected <- tapply(
    abs(
      cf$coefficient
    ) >
      tol,
    cf$predictor_key,
    any
  )

  names(
    selected
  )[
    selected
  ]
}

probability_matrix <- function(
    pred,
    class_levels = CLASS_LEVELS
) {
  d <- dim(pred)

  if (
    length(d) == 3L
  ) {
    pred <- pred[
      ,
      ,
      1L,
      drop = FALSE
    ]

    pred <- matrix(
      pred,
      nrow = d[1L],
      ncol = d[2L],
      dimnames = list(
        dimnames(pred)[[1L]],
        dimnames(pred)[[2L]]
      )
    )
  } else {
    pred <- as.matrix(
      pred
    )
  }

  missing <- setdiff(
    class_levels,
    colnames(pred)
  )

  assert(
    !length(missing),
    paste(
      "Prediction matrix missing classes:",
      paste(
        missing,
        collapse = ", "
      )
    )
  )

  pred <- pred[
    ,
    class_levels,
    drop = FALSE
  ]

  rs <- rowSums(pred)

  assert(
    all(
      is.finite(rs) &
        rs > 0
    ),
    "Invalid probability rows."
  )

  pred / rs
}

fit_final_model <- function(
    x,
    y,
    alpha,
    lambda,
    weights
) {
  glmnet::glmnet(
    x = x,
    y = y,
    family = "multinomial",
    type.multinomial = "grouped",
    alpha = alpha,
    lambda = lambda,
    standardize = TRUE,
    intercept = TRUE,
    weights = weights
  )
}

repeated_oof_and_stability <- function(
    x,
    y,
    alpha,
    lambda,
    repeats,
    seed
) {
  n <- nrow(x)

  prob_sum <- matrix(
    0,
    nrow = n,
    ncol =
      length(
        CLASS_LEVELS
      ),
    dimnames = list(
      rownames(x),
      CLASS_LEVELS
    )
  )

  prob_n <- integer(n)

  selection_count <- setNames(
    rep(
      0L,
      ncol(x)
    ),
    colnames(x)
  )

  n_resample_models <- 0L

  for (
    r in seq_len(
      repeats
    )
  ) {
    foldid <- make_stratified_foldid(
      y = y,
      k = PANEL_CV_FOLDS,
      seed =
        seed +
        50000L +
        r
    )

    for (
      fold in seq_len(
        PANEL_CV_FOLDS
      )
    ) {
      idx_train <- which(
        foldid != fold
      )

      idx_valid <- which(
        foldid == fold
      )

      y_sub <- factor(
        y[idx_train],
        levels = CLASS_LEVELS
      )

      assert(
        all(
          table(y_sub) > 0L
        ),
        "A resampling training subset lost a class."
      )

      w_sub <- make_class_weights(
        y_sub
      )

      fit <- fit_final_model(
        x = x[
          idx_train,
          ,
          drop = FALSE
        ],
        y = y_sub,
        alpha = alpha,
        lambda = lambda,
        weights = w_sub
      )

      pred <- stats::predict(
        fit,
        newx = x[
          idx_valid,
          ,
          drop = FALSE
        ],
        s = lambda,
        type = "response"
      )

      prob <- probability_matrix(
        pred
      )

      prob_sum[
        idx_valid,
        colnames(prob)
      ] <- prob_sum[
        idx_valid,
        colnames(prob)
      ] +
        prob

      prob_n[
        idx_valid
      ] <- prob_n[
        idx_valid
      ] +
        1L

      selected <- selected_predictors_from_fit(
        fit,
        lambda
      )

      if (length(selected)) {
        selection_count[
          selected
        ] <- selection_count[
          selected
        ] +
          1L
      }

      n_resample_models <-
        n_resample_models +
        1L
    }
  }

  assert(
    all(
      prob_n > 0L
    ),
    "Some TRAIN samples did not receive OOF predictions."
  )

  oof_prob <- prob_sum /
    prob_n

  stability <- data.frame(
    predictor_key =
      names(
        selection_count
      ),
    selection_count =
      as.integer(
        selection_count
      ),
    n_resample_models =
      n_resample_models,
    stability_frequency =
      as.numeric(
        selection_count
      ) /
      n_resample_models,
    stringsAsFactors = FALSE
  )

  list(
    oof_prob =
      oof_prob,
    stability =
      stability
  )
}

classification_metrics <- function(
    truth,
    prob,
    split_label
) {
  truth <- factor(
    truth,
    levels = CLASS_LEVELS
  )

  prob <- probability_matrix(
    prob
  )

  pred <- factor(
    CLASS_LEVELS[
      max.col(
        prob,
        ties.method = "first"
      )
    ],
    levels = CLASS_LEVELS
  )

  cm <- table(
    truth = truth,
    predicted = pred
  )

  accuracy <- mean(
    pred == truth
  )

  per_class <- lapply(
    CLASS_LEVELS,
    function(cls) {
      tp <- cm[
        cls,
        cls
      ]

      fn <- sum(
        cm[
          cls,
          ,
          drop = TRUE
        ]
      ) -
        tp

      fp <- sum(
        cm[
          ,
          cls,
          drop = TRUE
        ]
      ) -
        tp

      recall <- if (
        tp + fn > 0
      ) {
        tp /
          (
            tp + fn
          )
      } else {
        NA_real_
      }

      precision <- if (
        tp + fp > 0
      ) {
        tp /
          (
            tp + fp
          )
      } else {
        NA_real_
      }

      f1 <- if (
        is.finite(
          precision
        ) &&
        is.finite(
          recall
        ) &&
        precision + recall > 0
      ) {
        2 *
          precision *
          recall /
          (
            precision + recall
          )
      } else {
        NA_real_
      }

      data.frame(
        split =
          split_label,
        class =
          cls,
        recall =
          recall,
        precision =
          precision,
        f1 =
          f1,
        stringsAsFactors = FALSE
      )
    }
  )

  per_class <- do.call(
    rbind,
    per_class
  )

  balanced_accuracy <- mean(
    per_class$recall,
    na.rm = TRUE
  )

  macro_f1 <- mean(
    per_class$f1,
    na.rm = TRUE
  )

  truth_index <- match(
    as.character(truth),
    CLASS_LEVELS
  )

  eps <- 1e-15

  log_loss <- -mean(
    log(
      pmax(
        prob[
          cbind(
            seq_len(
              nrow(prob)
            ),
            truth_index
          )
        ],
        eps
      )
    )
  )

  Y <- matrix(
    0,
    nrow = nrow(prob),
    ncol = length(CLASS_LEVELS),
    dimnames = list(
      NULL,
      CLASS_LEVELS
    )
  )

  Y[
    cbind(
      seq_len(
        nrow(Y)
      ),
      truth_index
    )
  ] <- 1

  multiclass_brier <- mean(
    rowSums(
      (
        prob -
          Y
      )^2
    )
  )

  summary <- data.frame(
    split =
      split_label,
    n =
      length(truth),
    accuracy =
      accuracy,
    balanced_accuracy =
      balanced_accuracy,
    macro_f1 =
      macro_f1,
    log_loss =
      log_loss,
    multiclass_brier =
      multiclass_brier,
    stringsAsFactors = FALSE
  )

  list(
    summary =
      summary,
    per_class =
      per_class,
    confusion =
      cm,
    predicted =
      pred,
    probabilities =
      prob
  )
}

plot_confusion <- function(
    cm,
    title,
    path
) {
  d <- as.data.frame(
    cm,
    stringsAsFactors = FALSE
  )

  names(d) <- c(
    "truth",
    "predicted",
    "n"
  )

  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = predicted,
      y = truth,
      fill = n
    )
  ) +
    ggplot2::geom_tile(
      color = "white"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = n
      ),
      size = 4
    ) +
    ggplot2::labs(
      title = title,
      x = "Predicted",
      y = "Truth",
      fill = "N"
    ) +
    plot_theme(11)

  ggplot2::ggsave(
    path,
    p,
    width = 5.5,
    height = 4.8,
    dpi = PANEL_DPI
  )
}

# ==============================================================================
# 10. FIT ONE FINAL PANEL PER P1 / P2 / P3
# ==============================================================================

all_tuning <- list()
all_alpha_summary <- list()
all_stability <- list()
all_final_panel <- list()
all_coefficients <- list()
all_metrics <- list()
all_per_class <- list()
all_predictions <- list()

for (
  family in FAMILIES
) {
  message("")
  message("============================================================")
  message("[PANEL] ", family)
  message("============================================================")

  dat <- model_data[[family]]

  x_train <- dat$x_train
  y_train <- dat$y_train
  x_test <- dat$x_test
  y_test <- dat$y_test

  message(
    "TRAIN n=",
    nrow(x_train),
    " | candidates=",
    ncol(x_train)
  )

  # ---------------------------------------------------------------------------
  # 10.1 TRAIN-only repeated CV hyperparameter tuning.
  # ---------------------------------------------------------------------------

  tuning <- tune_elastic_net(
    x = x_train,
    y = y_train,
    family = family,
    repeats = PANEL_CV_REPEATS,
    alpha_grid = PANEL_ALPHA_GRID,
    seed =
      PANEL_SEED +
      match(
        family,
        FAMILIES
      ) *
      100000L
  )

  all_tuning[[family]] <-
    tuning$tuning

  all_alpha_summary[[family]] <-
    tuning$alpha_summary

  # Plot tuning.
  p_tune <- ggplot2::ggplot(
    tuning$alpha_summary,
    ggplot2::aes(
      x = alpha,
      y = mean_cvm_1se
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(
      size = 2.8
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin =
          mean_cvm_1se -
          sd_cvm_1se,
        ymax =
          mean_cvm_1se +
          sd_cvm_1se
      ),
      width = 0.03
    ) +
    ggplot2::labs(
      title = paste0(
        family,
        " — TRAIN-only Elastic Net tuning"
      ),
      subtitle =
        "Repeated stratified 3-fold CV; model selection uses deviance at lambda.1se.",
      x = "Elastic Net alpha",
      y = "Mean CV deviance at lambda.1se"
    ) +
    plot_theme(11)

  ggplot2::ggsave(
    file.path(
      DIR_FIGURES,
      paste0(
        "02_",
        family,
        "_TRAIN_TUNING.png"
      )
    ),
    p_tune,
    width = 7,
    height = 5,
    dpi = PANEL_DPI
  )

  # ---------------------------------------------------------------------------
  # 10.2 Fit full TRAIN using median repeated-CV lambda.1se.
  #
  # A priori fallback:
  #   if lambda.1se produces an empty panel, use median lambda.min.
  # ---------------------------------------------------------------------------

  weights_full <- make_class_weights(
    y_train
  )

  lambda_used <- tuning$lambda_1se
  lambda_rule <- "median_repeated_lambda_1se"

  fit_full <- fit_final_model(
    x = x_train,
    y = y_train,
    alpha = tuning$alpha,
    lambda = lambda_used,
    weights = weights_full
  )

  selected <- selected_predictors_from_fit(
    fit_full,
    lambda_used
  )

  if (!length(selected)) {
    warning(
      family,
      ": lambda.1se selected no predictors; using the pre-specified fallback lambda.min."
    )

    lambda_used <- tuning$lambda_min
    lambda_rule <- "fallback_median_repeated_lambda_min"

    fit_full <- fit_final_model(
      x = x_train,
      y = y_train,
      alpha = tuning$alpha,
      lambda = lambda_used,
      weights = weights_full
    )

    selected <- selected_predictors_from_fit(
      fit_full,
      lambda_used
    )
  }

  assert(
    length(selected) > 0L,
    paste(
      family,
      "produced an empty panel even at lambda.min."
    )
  )

  # ---------------------------------------------------------------------------
  # 10.3 TRAIN repeated OOF audit + selection stability.
  # ---------------------------------------------------------------------------

  resample_audit <- repeated_oof_and_stability(
    x = x_train,
    y = y_train,
    alpha = tuning$alpha,
    lambda = lambda_used,
    repeats = PANEL_CV_REPEATS,
    seed =
      PANEL_SEED +
      match(
        family,
        FAMILIES
      ) *
      200000L
  )

  stability <- resample_audit$stability
  stability$family <- family
  stability$selected_in_final_model <-
    stability$predictor_key %in%
    selected

  all_stability[[family]] <-
    stability

  # ---------------------------------------------------------------------------
  # 10.4 Final TEST prediction — TEST touched only here.
  # ---------------------------------------------------------------------------

  test_pred_raw <- stats::predict(
    fit_full,
    newx = x_test,
    s = lambda_used,
    type = "response"
  )

  test_prob <- probability_matrix(
    test_pred_raw
  )

  train_eval <- classification_metrics(
    truth = y_train,
    prob = resample_audit$oof_prob,
    split_label = "TRAIN_REPEATED_CV"
  )

  test_eval <- classification_metrics(
    truth = y_test,
    prob = test_prob,
    split_label = "TEST_FINAL"
  )

  metrics_family <- rbind(
    train_eval$summary,
    test_eval$summary
  )

  metrics_family$family <- family
  metrics_family$alpha <- tuning$alpha
  metrics_family$lambda <- lambda_used
  metrics_family$lambda_rule <- lambda_rule
  metrics_family$n_candidates_entering_model <-
    ncol(x_train)
  metrics_family$n_final_panel_features <-
    length(selected)

  all_metrics[[family]] <-
    metrics_family

  per_class_family <- rbind(
    train_eval$per_class,
    test_eval$per_class
  )

  per_class_family$family <- family

  all_per_class[[family]] <-
    per_class_family

  # ---------------------------------------------------------------------------
  # 10.5 Final coefficients and panel table.
  # ---------------------------------------------------------------------------

  coef_long <- extract_multinomial_coefficients(
    fit_full,
    lambda_used
  )

  coef_long$family <- family
  coef_long$alpha <- tuning$alpha
  coef_long$lambda <- lambda_used
  coef_long$lambda_rule <- lambda_rule

  all_coefficients[[family]] <-
    coef_long

  coef_feature <- coef_long[
    coef_long$predictor_key %in%
      selected,
    ,
    drop = FALSE
  ]

  coef_wide <- reshape(
    coef_feature[
      ,
      c(
        "predictor_key",
        "class",
        "coefficient"
      ),
      drop = FALSE
    ],
    idvar = "predictor_key",
    timevar = "class",
    direction = "wide"
  )

  names(coef_wide) <- sub(
    "^coefficient\\.",
    "coef_",
    names(coef_wide)
  )

  s_family <- biological_shortlist[
    biological_shortlist$family ==
      family,
    ,
    drop = FALSE
  ]

  panel <- merge(
    data.frame(
      predictor_key = selected,
      stringsAsFactors = FALSE
    ),
    s_family,
    by = "predictor_key",
    all.x = TRUE,
    sort = FALSE
  )

  panel <- merge(
    panel,
    stability[
      ,
      c(
        "predictor_key",
        "selection_count",
        "n_resample_models",
        "stability_frequency"
      ),
      drop = FALSE
    ],
    by = "predictor_key",
    all.x = TRUE,
    sort = FALSE
  )

  panel <- merge(
    panel,
    coef_wide,
    by = "predictor_key",
    all.x = TRUE,
    sort = FALSE
  )

  panel$alpha <-
    tuning$alpha

  panel$lambda <-
    lambda_used

  panel$lambda_rule <-
    lambda_rule

  panel <- panel[
    order(
      -panel$stability_frequency,
      match(
        panel$view_code,
        VIEW_LEVELS
      ),
      panel$feature_label
    ),
    ,
    drop = FALSE
  ]

  all_final_panel[[family]] <-
    panel

  # ---------------------------------------------------------------------------
  # 10.6 Prediction tables.
  # ---------------------------------------------------------------------------

  pred_train <- data.frame(
    family = family,
    sample_id =
      rownames(
        x_train
      ),
    split =
      "TRAIN_REPEATED_CV",
    truth =
      as.character(
        y_train
      ),
    predicted =
      as.character(
        train_eval$predicted
      ),
    resample_audit$oof_prob,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  pred_test <- data.frame(
    family = family,
    sample_id =
      rownames(
        x_test
      ),
    split =
      "TEST_FINAL",
    truth =
      as.character(
        y_test
      ),
    predicted =
      as.character(
        test_eval$predicted
      ),
    test_prob,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  all_predictions[[family]] <- rbind(
    pred_train,
    pred_test
  )

  # ---------------------------------------------------------------------------
  # 10.7 Family-specific files.
  # ---------------------------------------------------------------------------

  fam_table_dir <- file.path(
    DIR_TABLES,
    family
  )

  fam_fig_dir <- file.path(
    DIR_FIGURES,
    family
  )

  fam_model_dir <- file.path(
    DIR_MODELS,
    family
  )

  dir.create(
    fam_table_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    fam_fig_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    fam_model_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  write_csv(
    panel,
    file.path(
      fam_table_dir,
      paste0(
        family,
        "_FINAL_PANEL.csv"
      )
    )
  )

  write_csv(
    stability[
      order(
        -stability$stability_frequency
      ),
      ,
      drop = FALSE
    ],
    file.path(
      fam_table_dir,
      paste0(
        family,
        "_SELECTION_STABILITY.csv"
      )
    )
  )

  write_csv(
    metrics_family,
    file.path(
      fam_table_dir,
      paste0(
        family,
        "_METRICS.csv"
      )
    )
  )

  write_csv(
    per_class_family,
    file.path(
      fam_table_dir,
      paste0(
        family,
        "_PER_CLASS_METRICS.csv"
      )
    )
  )

  train_cm_out <- data.frame(
    truth =
      rownames(
        train_eval$confusion
      ),
    as.data.frame.matrix(
      train_eval$confusion
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  test_cm_out <- data.frame(
    truth =
      rownames(
        test_eval$confusion
      ),
    as.data.frame.matrix(
      test_eval$confusion
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write_csv(
    train_cm_out,
    file.path(
      fam_table_dir,
      paste0(
        family,
        "_TRAIN_REPEATED_CV_CONFUSION.csv"
      )
    )
  )

  write_csv(
    test_cm_out,
    file.path(
      fam_table_dir,
      paste0(
        family,
        "_TEST_FINAL_CONFUSION.csv"
      )
    )
  )

  write_csv(
    all_predictions[[family]],
    file.path(
      fam_table_dir,
      paste0(
        family,
        "_PREDICTIONS.csv"
      )
    )
  )

  write_csv(
    tuning$tuning,
    file.path(
      fam_table_dir,
      paste0(
        family,
        "_TUNING_ALL_REPEATS.csv"
      )
    )
  )

  write_csv(
    tuning$alpha_summary,
    file.path(
      fam_table_dir,
      paste0(
        family,
        "_TUNING_ALPHA_SUMMARY.csv"
      )
    )
  )

  saveRDS(
    list(
      family =
        family,
      fit =
        fit_full,
      alpha =
        tuning$alpha,
      lambda =
        lambda_used,
      lambda_rule =
        lambda_rule,
      selected_predictors =
        selected,
      class_levels =
        CLASS_LEVELS,
      candidate_columns =
        colnames(
          x_train
        ),
      biological_synthesis_run =
        SYNTHESIS_RUN,
      bio_program_run =
        BIO_RUN
    ),
    file.path(
      fam_model_dir,
      paste0(
        family,
        "_FINAL_PANEL_MODEL.rds"
      )
    )
  )

  # Stability figure: final panel only.
  panel_plot <- panel

  panel_plot$plot_label <- paste0(
    toupper(
      panel_plot$view_code
    ),
    " | ",
    panel_plot$feature_label
  )

  panel_plot$plot_label <- factor(
    panel_plot$plot_label,
    levels =
      rev(
        panel_plot$plot_label
      )
  )

  p_stability <- ggplot2::ggplot(
    panel_plot,
    ggplot2::aes(
      x = stability_frequency,
      y = plot_label
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = paste0(
        family,
        " — final predictive panel"
      ),
      subtitle =
        paste0(
          "Final full-TRAIN non-zero predictors; bars = selection frequency across ",
          PANEL_CV_REPEATS * PANEL_CV_FOLDS,
          " TRAIN resample fits."
        ),
      x = "Selection stability",
      y = NULL
    ) +
    ggplot2::coord_cartesian(
      xlim = c(
        0,
        1
      )
    ) +
    plot_theme(10)

  ggplot2::ggsave(
    file.path(
      fam_fig_dir,
      paste0(
        family,
        "_FINAL_PANEL_STABILITY.png"
      )
    ),
    p_stability,
    width = 9,
    height = max(
      5,
      min(
        18,
        3 +
          0.32 *
          nrow(panel_plot)
      )
    ),
    dpi = PANEL_DPI,
    limitsize = FALSE
  )

  # Coefficient heatmap.
  coef_plot <- coef_feature

  coef_plot$predictor_key <- factor(
    coef_plot$predictor_key,
    levels =
      rev(
        panel$predictor_key
      )
  )

  p_coef <- ggplot2::ggplot(
    coef_plot,
    ggplot2::aes(
      x = class,
      y = predictor_key,
      fill = coefficient
    )
  ) +
    ggplot2::geom_tile(
      color = "white"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf(
          "%.2f",
          coefficient
        )
      ),
      size = 2.8
    ) +
    ggplot2::scale_fill_gradient2(
      midpoint = 0,
      low = "#2166AC",
      mid = "white",
      high = "#B2182B"
    ) +
    ggplot2::labs(
      title = paste0(
        family,
        " — final multinomial coefficients"
      ),
      subtitle =
        "Coefficients from the model fitted once on all TRAIN.",
      x = "Class",
      y = NULL,
      fill = "Coefficient"
    ) +
    plot_theme(10)

  ggplot2::ggsave(
    file.path(
      fam_fig_dir,
      paste0(
        family,
        "_FINAL_PANEL_COEFFICIENTS.png"
      )
    ),
    p_coef,
    width = 8,
    height = max(
      5,
      min(
        18,
        3 +
          0.30 *
          length(
            selected
          )
      )
    ),
    dpi = PANEL_DPI,
    limitsize = FALSE
  )

  plot_confusion(
    train_eval$confusion,
    paste0(
      family,
      " — TRAIN repeated-CV confusion"
    ),
    file.path(
      fam_fig_dir,
      paste0(
        family,
        "_TRAIN_REPEATED_CV_CONFUSION.png"
      )
    )
  )

  plot_confusion(
    test_eval$confusion,
    paste0(
      family,
      " — TEST final confusion"
    ),
    file.path(
      fam_fig_dir,
      paste0(
        family,
        "_TEST_FINAL_CONFUSION.png"
      )
    )
  )

  message(
    "[PANEL] ",
    family,
    " | alpha=",
    sprintf(
      "%.2f",
      tuning$alpha
    ),
    " | lambda=",
    signif(
      lambda_used,
      4
    ),
    " | n_panel=",
    length(selected)
  )
}

# ==============================================================================
# 11. CROSS-ANALYSIS CONSENSUS
# ==============================================================================

tuning_all <- do.call(
  rbind,
  all_tuning
)

alpha_summary_all <- do.call(
  rbind,
  all_alpha_summary
)

stability_all <- do.call(
  rbind,
  all_stability
)

final_panel_all <- do.call(
  rbind,
  lapply(
    names(
      all_final_panel
    ),
    function(family) {
      z <- all_final_panel[[family]]
      z$family <- family
      z
    }
  )
)

coefficients_all <- do.call(
  rbind,
  all_coefficients
)

metrics_all <- do.call(
  rbind,
  all_metrics
)

per_class_all <- do.call(
  rbind,
  all_per_class
)

predictions_all <- do.call(
  rbind,
  all_predictions
)

rownames(
  final_panel_all
) <- NULL

# Predictor/biomarker consensus across P1/P2/P3.
consensus_split <- split(
  final_panel_all,
  paste(
    final_panel_all$view_code,
    final_panel_all$feature_original,
    sep = "\r"
  )
)

consensus_rows <- lapply(
  consensus_split,
  function(z) {
    families_selected <- sort(
      unique(
        as.character(
          z$family
        )
      )
    )

    n_families <- length(
      families_selected
    )

    predictive_scope <- if (
      n_families >= 3L
    ) {
      "Selected P1-P2-P3"
    } else if (
      n_families == 2L
    ) {
      "Selected two analyses"
    } else {
      "Selected one analysis"
    }

    data.frame(
      view_code =
        first_nonempty(
          z$view_code
        ),

      feature_original =
        first_nonempty(
          z$feature_original
        ),

      feature_label =
        first_nonempty(
          z$feature_label
        ),

      n_families_selected =
        n_families,

      families_selected =
        paste(
          families_selected,
          collapse = ";"
        ),

      predictive_selection_scope =
        predictive_scope,

      mean_stability_frequency =
        mean_or_na(
          z$stability_frequency
        ),

      min_stability_frequency =
        if (
          any(
            is.finite(
              z$stability_frequency
            )
          )
        ) {
          min(
            z$stability_frequency,
            na.rm = TRUE
          )
        } else {
          NA_real_
        },

      max_stability_frequency =
        max_or_na(
          z$stability_frequency
        ),

      biological_pathway_scopes =
        collapse_unique(
          z$pathway_analysis_scopes
        ),

      biological_phenotype_scopes =
        collapse_unique(
          z$pathway_phenotype_scopes
        ),

      phenotypic_sides =
        collapse_unique(
          z$phenotypic_sides
        ),

      multimodal_reinforcement_status =
        collapse_unique(
          z$multimodal_reinforcement_status
        ),

      stringsAsFactors = FALSE
    )
  }
)

consensus <- do.call(
  rbind,
  consensus_rows
)

rownames(
  consensus
) <- NULL

consensus <- consensus[
  order(
    -consensus$n_families_selected,
    -consensus$mean_stability_frequency,
    match(
      consensus$view_code,
      VIEW_LEVELS
    ),
    consensus$feature_label
  ),
  ,
  drop = FALSE
]

# IMPORTANT:
# This consensus is descriptive. A feature selected in only one analysis is NOT
# automatically inferior to one selected in all three analyses.

write_csv(
  tuning_all,
  file.path(
    DIR_TABLES,
    "04_TUNING_ALL_FAMILIES.csv"
  )
)

write_csv(
  alpha_summary_all,
  file.path(
    DIR_TABLES,
    "05_TUNING_ALPHA_SUMMARY_ALL_FAMILIES.csv"
  )
)

write_csv(
  stability_all,
  file.path(
    DIR_TABLES,
    "06_SELECTION_STABILITY_ALL_FAMILIES.csv"
  )
)

write_csv(
  final_panel_all,
  file.path(
    DIR_TABLES,
    "07_FINAL_PANEL_ALL_FAMILIES.csv"
  )
)

write_csv(
  coefficients_all,
  file.path(
    DIR_TABLES,
    "08_FINAL_COEFFICIENTS_ALL_FAMILIES.csv"
  )
)

write_csv(
  metrics_all,
  file.path(
    DIR_TABLES,
    "09_MODEL_METRICS_ALL_FAMILIES.csv"
  )
)

write_csv(
  per_class_all,
  file.path(
    DIR_TABLES,
    "10_PER_CLASS_METRICS_ALL_FAMILIES.csv"
  )
)

write_csv(
  predictions_all,
  file.path(
    DIR_TABLES,
    "11_PREDICTIONS_ALL_FAMILIES.csv"
  )
)

write_csv(
  consensus,
  file.path(
    DIR_TABLES,
    "12_CROSS_ANALYSIS_PANEL_CONSENSUS.csv"
  )
)

# ==============================================================================
# 12. SUMMARY FIGURES
# ==============================================================================

metrics_long <- metrics_all[
  ,
  c(
    "family",
    "split",
    "accuracy",
    "balanced_accuracy",
    "macro_f1"
  ),
  drop = FALSE
]

metrics_plot <- reshape(
  metrics_long,
  varying = c(
    "accuracy",
    "balanced_accuracy",
    "macro_f1"
  ),
  v.names = "value",
  timevar = "metric",
  times = c(
    "Accuracy",
    "Balanced accuracy",
    "Macro F1"
  ),
  direction = "long"
)

metrics_plot$metric <- factor(
  metrics_plot$metric,
  levels = c(
    "Accuracy",
    "Balanced accuracy",
    "Macro F1"
  )
)

p_metrics <- ggplot2::ggplot(
  metrics_plot,
  ggplot2::aes(
    x = metric,
    y = value,
    shape = split,
    group = split
  )
) +
  ggplot2::geom_point(
    size = 3
  ) +
  ggplot2::facet_wrap(
    ~ family
  ) +
  ggplot2::coord_cartesian(
    ylim = c(
      0,
      1
    )
  ) +
  ggplot2::labs(
    title =
      "Final predictive-panel performance",
    subtitle =
      "TRAIN = repeated-CV development estimate; TEST = single final holdout evaluation.",
    x = NULL,
    y = "Metric",
    shape = "Split"
  ) +
  plot_theme(11)

ggplot2::ggsave(
  file.path(
    DIR_FIGURES,
    "03_FINAL_PANEL_PERFORMANCE_SUMMARY.png"
  ),
  p_metrics,
  width = 10,
  height = 5.8,
  dpi = PANEL_DPI
)

if (nrow(consensus)) {
  cons_plot <- consensus

  cons_plot$plot_label <- paste0(
    toupper(
      cons_plot$view_code
    ),
    " | ",
    cons_plot$feature_label
  )

  cons_plot$plot_label <- factor(
    cons_plot$plot_label,
    levels =
      rev(
        cons_plot$plot_label
      )
  )

  p_consensus <- ggplot2::ggplot(
    cons_plot,
    ggplot2::aes(
      x = n_families_selected,
      y = plot_label
    )
  ) +
    ggplot2::geom_point(
      size = 3
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(
        FAMILIES
      ),
      limits = c(
        0.8,
        length(FAMILIES) + 0.2
      )
    ) +
    ggplot2::labs(
      title =
        "Cross-analysis predictive-panel selection",
      subtitle =
        "Number of P1/P2/P3 models selecting each biomarker. This is an audit, not a quality score.",
      x = "Number of analyses selecting the biomarker",
      y = NULL
    ) +
    plot_theme(10)

  ggplot2::ggsave(
    file.path(
      DIR_FIGURES,
      "04_CROSS_ANALYSIS_PANEL_CONSENSUS.png"
    ),
    p_consensus,
    width = 9,
    height = max(
      5,
      min(
        24,
        3 +
          0.30 *
          nrow(cons_plot)
      )
    ),
    dpi = PANEL_DPI,
    limitsize = FALSE
  )
}

# ==============================================================================
# 13. RUN MANIFEST + README
# ==============================================================================

manifest <- data.frame(
  item = c(
    "BIO_PROGRAM_RUN",
    "BIO_SYNTHESIS_RUN",
    "PANEL_OUTPUT_ROOT",
    "FAMILIES",
    "PANEL_SEED",
    "PANEL_CV_REPEATS",
    "PANEL_CV_FOLDS",
    "PANEL_ALPHA_GRID",
    "PANEL_USE_CLASS_WEIGHTS",
    "candidate_rule_TX",
    "candidate_rule_PR_ME_CL",
    "tuning_rule",
    "lambda_rule",
    "test_rule"
  ),
  value = c(
    BIO_RUN,
    SYNTHESIS_RUN,
    OUT_ROOT,
    paste(
      FAMILIES,
      collapse = ","
    ),
    as.character(
      PANEL_SEED
    ),
    as.character(
      PANEL_CV_REPEATS
    ),
    as.character(
      PANEL_CV_FOLDS
    ),
    paste(
      PANEL_ALPHA_GRID,
      collapse = ","
    ),
    as.character(
      PANEL_USE_CLASS_WEIGHTS
    ),
    "Frozen TX drivers only",
    "Frozen non-TX crossomic links with finite non-zero pathway_feature_cosine",
    "Repeated stratified TRAIN-only CV minimizing mean deviance at lambda.1se across alpha",
    "Median repeated lambda.1se; pre-specified fallback to median lambda.min only if panel is empty",
    "TEST is predicted once after all candidate/model decisions are frozen"
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
  "METATEST — FINAL BIOLOGICAL SHORTLIST + PREDICTIVE PANEL",
  "",
  "PURPOSE",
  "  This stage is intentionally separate from 09_build_biological_synthesis_figures.R.",
  "  Script 09 answers biological questions Q1-Q7.",
  "  Script 10 converts the frozen biology into candidate biomarkers and then",
  "  develops sparse NP/FD/SO predictive panels.",
  "",
  "BIOLOGICAL SHORTLIST",
  "  - TX: only frozen TX drivers are eligible.",
  "  - PR/ME/CL: frozen cross-omic links are eligible.",
  "  - No minimum representation per omic.",
  "  - No omic quota.",
  "  - Shared/specific, phenotype-side and multimodal annotations are preserved.",
  "  - Q7 ALL/TEST separation values are descriptive annotations only.",
  "",
  "PREDICTIVE MODEL",
  "  - Multinomial grouped Elastic Net.",
  "  - Alpha and lambda are selected using TRAIN only.",
  "  - Repeated stratified 3-fold CV is used because the smallest TRAIN class is FD.",
  "  - lambda.1se is preferred to favor a smaller panel.",
  "  - If lambda.1se yields zero predictors, the pre-specified fallback is lambda.min.",
  "  - Stability is audited across repeated TRAIN resample fits.",
  "",
  "TEST",
  "  - TEST does not select candidates.",
  "  - TEST does not tune alpha or lambda.",
  "  - TEST does not choose the final panel.",
  "  - TEST is evaluated once after the full-TRAIN panel is frozen.",
  "",
  "P1 / P2 / P3",
  "  - They remain distinct analyses.",
  "  - No winner is selected among P1/P2/P3.",
  "  - A cross-analysis table reports whether a final biomarker is selected in",
  "    one, two or all three analyses.",
  "  - Selection in one analysis is not considered biologically inferior.",
  "",
  "MAIN OUTPUTS",
  "  tables/02_BIOLOGICAL_FEATURE_SHORTLIST.csv",
  "  tables/07_FINAL_PANEL_ALL_FAMILIES.csv",
  "  tables/09_MODEL_METRICS_ALL_FAMILIES.csv",
  "  tables/11_PREDICTIONS_ALL_FAMILIES.csv",
  "  tables/12_CROSS_ANALYSIS_PANEL_CONSENSUS.csv",
  "  tables/P1/P1_FINAL_PANEL.csv",
  "  tables/P2/P2_FINAL_PANEL.csv",
  "  tables/P3/P3_FINAL_PANEL.csv",
  "  models/P1/P1_FINAL_PANEL_MODEL.rds",
  "  models/P2/P2_FINAL_PANEL_MODEL.rds",
  "  models/P3/P3_FINAL_PANEL_MODEL.rds"
)

writeLines(
  readme,
  con = file.path(
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
  con = file.path(
    OUT_ROOT,
    "RUN_COMPLETE.flag"
  )
)

latest_root <- file.path(
  BIO_RUN,
  "FINAL_PREDICTIVE_PANEL"
)

dir.create(
  latest_root,
  recursive = TRUE,
  showWarnings = FALSE
)

writeLines(
  OUT_ROOT,
  con = file.path(
    latest_root,
    "LATEST_RUN.txt"
  )
)

message("")
message("============================================================")
message("FINAL PREDICTIVE PANEL COMPLETE")
message("Output: ", OUT_ROOT)
message("============================================================")

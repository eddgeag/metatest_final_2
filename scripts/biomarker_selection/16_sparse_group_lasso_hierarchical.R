#!/usr/bin/env Rscript

# ==============================================================================
# METATEST — 16_sparse_group_lasso_hierarchical.R
#
# FINAL STRUCTURED SUPERVISED PANEL ATTEMPT
#
# Scientific idea:
#
#   script-10 frozen biological candidates
#        -> outcome-free PRIMARY PROGRAM assignment
#        -> non-overlapping PROGRAM x OMIC groups
#        -> multinomial SPARSE GROUP LASSO on ORIGINAL FEATURES
#        -> group sparsity + within-group feature sparsity simultaneously
#        -> one direct measurable panel
#
# Penalty:
#
#   loss(beta)
#     + lambda * [
#         (1-alpha) * sum_g w_g ||beta_g||_2
#         + alpha * sum_j |beta_j|
#       ]
#
# alpha = 0 : pure group lasso
# alpha = 1 : pure lasso
# 0 < alpha < 1 : sparse group lasso
#
# We intentionally restrict the alpha grid to (0,1), because this script is
# testing the structured sparse-group hypothesis rather than reproducing the
# flat LASSO or pure group-lasso endpoints.
#
# IMPORTANT:
# * ORIGINAL features are modeled; PCA is NOT a selection bottleneck.
# * No omic quota.
# * No fixed panel size.
# * No feature/block stability threshold decides membership in the final panel.
#   Membership is determined by the final full-TRAIN sparse-group model chosen
#   with a repeated stratified 3-fold 1-SE-like rule.
# * Stability is exported as an audit annotation.
# * TEST is not used for tuning or selection.
# * TEST has already been inspected during method development, therefore it is
#   explicitly labeled TEST_EXPLORATORY_REUSED.
# * With FD=3, fully nested multinomial CV is not feasible; TRAIN performance
#   is a conditional development estimate.
#
# Package:
#   msgl — multinomial sparse group lasso.
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
PANEL10_RUN <- Sys.getenv("PANEL10_RUN", unset = "")

FAMILIES <- trimws(
  strsplit(
    Sys.getenv("PANEL16_FAMILIES", unset = "P1,P2,P3"),
    ",",
    fixed = TRUE
  )[[1]]
)
FAMILIES <- FAMILIES[nzchar(FAMILIES)]

CLASS_LEVELS <- c("NP", "FD", "SO")
VIEW_LEVELS <- c("tx", "pr", "me", "cl")

SEED <- as.integer(
  Sys.getenv("PANEL16_SEED", unset = "20261001")
)

N_FOLDS <- as.integer(
  Sys.getenv("PANEL16_FOLDS", unset = "3")
)

CV_REPEATS <- as.integer(
  Sys.getenv("PANEL16_CV_REPEATS", unset = "20")
)

STABILITY_REPEATS <- as.integer(
  Sys.getenv("PANEL16_STABILITY_REPEATS", unset = "30")
)

ALPHA_GRID <- as.numeric(
  trimws(
    strsplit(
      Sys.getenv(
        "PANEL16_ALPHA_GRID",
        unset = "0.15,0.35,0.50,0.65,0.85"
      ),
      ",",
      fixed = TRUE
    )[[1]]
  )
)

# msgl uses a scalar lambda argument as the minimum FRACTION of lambda.max
# when a regularization path is generated. d controls path length.
LAMBDA_MIN_REL <- as.numeric(
  Sys.getenv("PANEL16_LAMBDA_MIN_REL", unset = "0.05")
)

PATH_LENGTH <- as.integer(
  Sys.getenv("PANEL16_PATH_LENGTH", unset = "50")
)

USE_CLASS_WEIGHTS <- tolower(
  Sys.getenv("PANEL16_USE_CLASS_WEIGHTS", unset = "true")
) %in% c("true", "1", "yes", "y")

STABILITY_AUDIT_THRESHOLD <- as.numeric(
  Sys.getenv("PANEL16_STABILITY_AUDIT_THRESHOLD", unset = "0.50")
)

DPI <- as.integer(
  Sys.getenv("PANEL16_DPI", unset = "180")
)

OUTPUT_ROOT_ENV <- Sys.getenv(
  "PANEL16_OUTPUT_ROOT",
  unset = ""
)

if (!requireNamespace("msgl", quietly = TRUE)) {
  stop(
    paste0(
      "Package 'msgl' is required.\n",
      "Install it with:\n",
      "install.packages('msgl')"
    ),
    call. = FALSE
  )
}

if (!requireNamespace("sglOptim", quietly = TRUE)) {
  stop(
    "Package 'sglOptim' is required by msgl.",
    call. = FALSE
  )
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "Package 'ggplot2' is required.",
    call. = FALSE
  )
}

if (N_FOLDS != 3L) {
  stop(
    "This script is designed for the current 3-fold stratified TRAIN design.",
    call. = FALSE
  )
}

if (
  any(!is.finite(ALPHA_GRID)) ||
  any(ALPHA_GRID <= 0) ||
  any(ALPHA_GRID >= 1)
) {
  stop(
    "ALPHA_GRID must contain only values strictly inside (0,1).",
    call. = FALSE
  )
}

if (
  !is.finite(LAMBDA_MIN_REL) ||
  LAMBDA_MIN_REL <= 0 ||
  LAMBDA_MIN_REL >= 1
) {
  stop(
    "LAMBDA_MIN_REL must be in (0,1).",
    call. = FALSE
  )
}

if (
  !is.finite(PATH_LENGTH) ||
  PATH_LENGTH < 10L
) {
  stop(
    "PATH_LENGTH must be >=10.",
    call. = FALSE
  )
}

# ==============================================================================
# 1. GENERIC HELPERS
# ==============================================================================

assert <- function(ok, msg) {
  if (!isTRUE(ok)) {
    stop(msg, call. = FALSE)
  }
}

canonical <- function(path, must_work = TRUE) {
  normalizePath(
    path.expand(path),
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

rbind_fill <- function(...) {
  xs <- list(...)

  if (
    length(xs) == 1L &&
    is.list(xs[[1L]]) &&
    !is.data.frame(xs[[1L]])
  ) {
    xs <- xs[[1L]]
  }

  xs <- xs[
    vapply(
      xs,
      is.data.frame,
      logical(1)
    )
  ]

  if (!length(xs)) {
    return(data.frame())
  }

  all_names <- Reduce(
    union,
    lapply(xs, names)
  )

  xs <- lapply(
    xs,
    function(x) {
      miss <- setdiff(
        all_names,
        names(x)
      )

      for (nm in miss) {
        x[[nm]] <- NA
      }

      x[
        ,
        all_names,
        drop = FALSE
      ]
    }
  )

  out <- do.call(
    rbind,
    xs
  )

  rownames(out) <- NULL
  out
}

theme_panel <- function(base_size = 11) {
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
# 2. LOCATE FROZEN RUNS
# ==============================================================================

find_latest_selection_run <- function(
    upstream
) {
  root <- file.path(
    upstream,
    "downstream",
    "seleccion_mejor_escenario_train"
  )

  assert(
    dir.exists(root),
    paste(
      "Selection root does not exist:",
      root
    )
  )

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

  assert(
    length(dirs) > 0L,
    "No completed scenario-selection run found."
  )

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

find_latest_bio_run <- function(
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
        canonical(p)
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
    "No completed BIOLOGICAL_PROGRAMS run found."
  )

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

resolve_bio_run <- function() {
  if (nzchar(BIO_PROGRAM_RUN)) {
    return(
      canonical(
        BIO_PROGRAM_RUN
      )
    )
  }

  selection_run <- if (
    nzchar(SELECTION_RUN)
  ) {
    canonical(
      SELECTION_RUN
    )
  } else {
    find_latest_selection_run(
      UPSTREAM_DIR
    )
  }

  find_latest_bio_run(
    selection_run
  )
}

find_panel10_run <- function(
    bio_run
) {
  if (nzchar(PANEL10_RUN)) {
    return(
      canonical(
        PANEL10_RUN
      )
    )
  }

  root <- file.path(
    bio_run,
    "FINAL_PREDICTIVE_PANEL"
  )

  assert(
    dir.exists(root),
    paste(
      "FINAL_PREDICTIVE_PANEL does not exist:",
      root
    )
  )

  required <- c(
    file.path(
      "tables",
      "00_INPUT_AUDIT.csv"
    ),
    file.path(
      "tables",
      "01_BIOLOGICAL_CANDIDATE_OCCURRENCES.csv"
    ),
    file.path(
      "tables",
      "02_BIOLOGICAL_FEATURE_SHORTLIST.csv"
    )
  )

  is_complete <- function(p) {
    nzchar(p) &&
      dir.exists(p) &&
      all(
        file.exists(
          file.path(
            p,
            required
          )
        )
      )
  }

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

    if (is_complete(p)) {
      return(
        canonical(p)
      )
    }
  }

  dirs <- list.dirs(
    root,
    recursive = FALSE,
    full.names = TRUE
  )

  dirs <- dirs[
    vapply(
      dirs,
      is_complete,
      logical(1)
    )
  ]

  assert(
    length(dirs) > 0L,
    "No completed script-10 run found."
  )

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

BIO_RUN <- resolve_bio_run()
PANEL10 <- find_panel10_run(
  BIO_RUN
)

if (nzchar(OUTPUT_ROOT_ENV)) {
  OUT_ROOT <- path.expand(
    OUTPUT_ROOT_ENV
  )
} else {
  OUT_ROOT <- file.path(
    BIO_RUN,
    "FINAL_SPARSE_GROUP_LASSO_PANEL",
    paste0(
      format(
        Sys.time(),
        "%Y%m%dT%H%M%S"
      ),
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

OUT_ROOT <- canonical(
  OUT_ROOT
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
message("METATEST — MULTINOMIAL SPARSE GROUP LASSO")
message("BIOLOGICAL_PROGRAMS run: ", BIO_RUN)
message("Script-10 run: ", PANEL10)
message("Output: ", OUT_ROOT)
message("============================================================")

# ==============================================================================
# 3. READ FROZEN BIOLOGICAL CANDIDATES
# ==============================================================================

input_audit <- read_csv_checked(
  file.path(
    PANEL10,
    "tables",
    "00_INPUT_AUDIT.csv"
  ),
  required = c(
    "family",
    "ready_for_modeling",
    "n_train",
    "n_test"
  )
)

occ <- read_csv_checked(
  file.path(
    PANEL10,
    "tables",
    "01_BIOLOGICAL_CANDIDATE_OCCURRENCES.csv"
  ),
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
    "tx_driver_weight",
    "pathway_feature_cosine",
    "predictor_key",
    "eligible_biological"
  )
)

shortlist <- read_csv_checked(
  file.path(
    PANEL10,
    "tables",
    "02_BIOLOGICAL_FEATURE_SHORTLIST.csv"
  ),
  required = c(
    "family",
    "predictor_key",
    "view_code",
    "feature_original",
    "feature_label",
    "candidate_reason"
  )
)

input_audit$family <- as.character(
  input_audit$family
)

occ$family <- as.character(
  occ$family
)

occ$view_code <- tolower(
  as.character(
    occ$view_code
  )
)

occ$predictor_key <- as.character(
  occ$predictor_key
)

occ$eligible_biological <- as_flag(
  occ$eligible_biological
)

occ <- occ[
  occ$family %in% FAMILIES &
    occ$view_code %in% VIEW_LEVELS &
    !is.na(occ$eligible_biological) &
    occ$eligible_biological,
  ,
  drop = FALSE
]

shortlist$family <- as.character(
  shortlist$family
)

shortlist$view_code <- tolower(
  as.character(
    shortlist$view_code
  )
)

shortlist$predictor_key <- as.character(
  shortlist$predictor_key
)

shortlist <- shortlist[
  shortlist$family %in% FAMILIES &
    shortlist$view_code %in% VIEW_LEVELS,
  ,
  drop = FALSE
]

assert(
  nrow(occ) > 0L,
  "No eligible biological occurrences."
)

assert(
  nrow(shortlist) > 0L,
  "No biological shortlist candidates."
)

# ==============================================================================
# 4. OUTCOME-FREE PRIMARY PROGRAM ASSIGNMENT
# ==============================================================================
#
# Same principle used in the previous hierarchical scripts:
#
# TX:
#   RMS absolute tx_driver_weight within each candidate program.
#
# PR / ME / CL:
#   RMS absolute pathway_feature_cosine within each candidate program.
#
# Ties:
#   more linked pathways -> larger maximal link -> program_id.
#
# Each original feature gets exactly ONE primary program so the sparse-group
# grouping is non-overlapping, as required by this analysis.
# ==============================================================================

occ$tx_driver_weight_num <- suppressWarnings(
  as.numeric(
    occ$tx_driver_weight
  )
)

occ$pathway_feature_cosine_num <- suppressWarnings(
  as.numeric(
    occ$pathway_feature_cosine
  )
)

occ$link_abs <- ifelse(
  occ$view_code == "tx",
  abs(
    occ$tx_driver_weight_num
  ),
  abs(
    occ$pathway_feature_cosine_num
  )
)

occ$link_abs[
  !is.finite(
    occ$link_abs
  )
] <- 0

occ$pathway_key_local <- paste(
  occ$library,
  occ$pathway,
  sep = "::"
)

affinity_split <- split(
  occ,
  paste(
    occ$family,
    occ$predictor_key,
    occ$program_id,
    sep = "\r"
  )
)

affinity_rows <- lapply(
  affinity_split,
  function(z) {
    links <- z$link_abs[
      is.finite(
        z$link_abs
      )
    ]

    if (!length(links)) {
      links <- 0
    }

    data.frame(
      family =
        first_nonempty(
          z$family
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
      program_id =
        first_nonempty(
          z$program_id
        ),
      program_short =
        first_nonempty(
          z$program_short
        ),
      program_affinity_rms =
        sqrt(
          mean(
            links^2
          )
        ),
      program_affinity_max =
        max(
          links
        ),
      n_linked_pathways_in_program =
        length(
          unique(
            z$pathway_key_local
          )
        ),
      phenotypic_sides_in_program =
        collapse_unique(
          z$phenotypic_side
        ),
      stringsAsFactors = FALSE
    )
  }
)

program_affinity <- do.call(
  rbind,
  affinity_rows
)

feature_split <- split(
  program_affinity,
  paste(
    program_affinity$family,
    program_affinity$predictor_key,
    sep = "\r"
  )
)

primary_program <- do.call(
  rbind,
  lapply(
    feature_split,
    function(z) {
      z <- z[
        order(
          -z$program_affinity_rms,
          -z$n_linked_pathways_in_program,
          -z$program_affinity_max,
          z$program_id
        ),
        ,
        drop = FALSE
      ]

      out <- z[
        1L,
        ,
        drop = FALSE
      ]

      out$n_candidate_programs_for_feature <-
        nrow(z)

      out$all_candidate_programs_for_feature <-
        paste(
          sort(
            unique(
              z$program_id
            )
          ),
          collapse = ";"
        )

      out
    }
  )
)

rownames(
  primary_program
) <- NULL

primary_program$program_view_block <- paste(
  primary_program$family,
  primary_program$program_id,
  toupper(
    primary_program$view_code
  ),
  sep = "::"
)

write_csv(
  program_affinity,
  file.path(
    DIR_TABLES,
    "00_FEATURE_PROGRAM_AFFINITY_ALL.csv"
  )
)

write_csv(
  primary_program,
  file.path(
    DIR_TABLES,
    "01_PRIMARY_PROGRAM_ASSIGNMENT.csv"
  )
)

# ==============================================================================
# 5. LOAD EXACT TRAIN / TEST INPUT
# ==============================================================================

extract_group_table <- function(g) {
  if (is.data.frame(g)) {
    assert(
      all(
        c(
          "id",
          "grupo"
        ) %in% names(g)
      ),
      "Group table must contain id and grupo."
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
    "Unrecognized group structure.",
    call. = FALSE
  )
}

load_family_input <- function(
    family
) {
  z <- input_audit[
    input_audit$family == family,
    ,
    drop = FALSE
  ]

  assert(
    nrow(z) == 1L,
    paste(
      "Expected one input row for",
      family
    )
  )

  path <- as.character(
    z$ready_for_modeling[1L]
  )

  assert(
    file.exists(path),
    paste(
      "Missing ready_for_modeling:",
      path
    )
  )

  rfm <- readRDS(
    path
  )

  assert(
    all(
      c(
        VIEW_LEVELS,
        "grupo_train",
        "grupo_test"
      ) %in% names(rfm)
    ),
    paste(
      "Incomplete ready_for_modeling:",
      family
    )
  )

  train_meta <- extract_group_table(
    rfm$grupo_train
  )

  test_meta <- extract_group_table(
    rfm$grupo_test
  )

  train_meta$grupo <- factor(
    train_meta$grupo,
    levels = CLASS_LEVELS
  )

  test_meta$grupo <- factor(
    test_meta$grupo,
    levels = CLASS_LEVELS
  )

  assert(
    !anyNA(train_meta$grupo),
    paste(
      "Bad TRAIN labels:",
      family
    )
  )

  assert(
    !anyNA(test_meta$grupo),
    paste(
      "Bad TEST labels:",
      family
    )
  )

  tr <- list()
  te <- list()

  for (view in VIEW_LEVELS) {
    xtr <- as.matrix(
      rfm[[view]]$train
    )

    xte <- as.matrix(
      rfm[[view]]$test
    )

    assert(
      all(
        train_meta$sample_id %in%
          colnames(xtr)
      ),
      paste(
        "TRAIN sample mismatch:",
        family,
        view
      )
    )

    assert(
      all(
        test_meta$sample_id %in%
          colnames(xte)
      ),
      paste(
        "TEST sample mismatch:",
        family,
        view
      )
    )

    tr[[view]] <- xtr[
      ,
      train_meta$sample_id,
      drop = FALSE
    ]

    te[[view]] <- xte[
      ,
      test_meta$sample_id,
      drop = FALSE
    ]
  }

  list(
    train_meta = train_meta,
    test_meta = test_meta,
    train = tr,
    test = te,
    ready_for_modeling =
      canonical(path)
  )
}

family_input <- lapply(
  FAMILIES,
  load_family_input
)

names(
  family_input
) <- FAMILIES

# ==============================================================================
# 6. ORIGINAL-FEATURE MODEL MATRICES + PROGRAM x OMIC GROUPS
# ==============================================================================

build_model_matrix <- function(
    family,
    inp
) {
  s <- shortlist[
    shortlist$family == family,
    ,
    drop = FALSE
  ]

  p <- primary_program[
    primary_program$family == family,
    ,
    drop = FALSE
  ]

  s <- merge(
    s,
    p[
      ,
      c(
        "family",
        "predictor_key",
        "program_id",
        "program_short",
        "program_affinity_rms",
        "program_affinity_max",
        "n_linked_pathways_in_program",
        "phenotypic_sides_in_program",
        "n_candidate_programs_for_feature",
        "all_candidate_programs_for_feature",
        "program_view_block"
      ),
      drop = FALSE
    ],
    by = c(
      "family",
      "predictor_key"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  original_order <- shortlist$predictor_key[
    shortlist$family == family
  ]

  s <- s[
    match(
      original_order,
      s$predictor_key
    ),
    ,
    drop = FALSE
  ]

  assert(
    !anyNA(
      s$program_view_block
    ),
    paste(
      "Primary-program mapping missing:",
      family
    )
  )

  train_ids <- inp$train_meta$sample_id
  test_ids <- inp$test_meta$sample_id

  tr_cols <- list()
  te_cols <- list()
  audit <- list()
  rr <- 1L

  for (i in seq_len(nrow(s))) {
    view <- as.character(
      s$view_code[i]
    )

    feature <- as.character(
      s$feature_original[i]
    )

    key <- as.character(
      s$predictor_key[i]
    )

    xtr_view <- inp$train[[view]]
    xte_view <- inp$test[[view]]

    if (
      !feature %in%
        rownames(xtr_view)
    ) {
      audit[[rr]] <- data.frame(
        family = family,
        predictor_key = key,
        technically_eligible = FALSE,
        reason = "ABSENT_TRAIN",
        stringsAsFactors = FALSE
      )

      rr <- rr + 1L
      next
    }

    a <- as.numeric(
      xtr_view[
        feature,
        train_ids
      ]
    )

    if (
      any(
        !is.finite(a)
      )
    ) {
      audit[[rr]] <- data.frame(
        family = family,
        predictor_key = key,
        technically_eligible = FALSE,
        reason = "NONFINITE_TRAIN",
        stringsAsFactors = FALSE
      )

      rr <- rr + 1L
      next
    }

    sd_a <- stats::sd(a)

    if (
      !is.finite(sd_a) ||
      sd_a <= 0
    ) {
      audit[[rr]] <- data.frame(
        family = family,
        predictor_key = key,
        technically_eligible = FALSE,
        reason = "ZERO_VARIANCE_TRAIN",
        stringsAsFactors = FALSE
      )

      rr <- rr + 1L
      next
    }

    b <- if (
      feature %in%
        rownames(xte_view)
    ) {
      as.numeric(
        xte_view[
          feature,
          test_ids
        ]
      )
    } else {
      rep(
        NA_real_,
        length(test_ids)
      )
    }

    n_imp <- sum(
      !is.finite(b)
    )

    if (n_imp > 0L) {
      b[
        !is.finite(b)
      ] <- stats::median(
        a
      )
    }

    tr_cols[[key]] <- a
    te_cols[[key]] <- b

    audit[[rr]] <- data.frame(
      family = family,
      predictor_key = key,
      technically_eligible = TRUE,
      reason = "OK",
      test_values_imputed_from_TRAIN_median =
        n_imp,
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  assert(
    length(tr_cols) >= 2L,
    paste(
      "Too few eligible predictors:",
      family
    )
  )

  x_train <- do.call(
    cbind,
    tr_cols
  )

  x_test <- do.call(
    cbind,
    te_cols
  )

  rownames(
    x_train
  ) <- train_ids

  rownames(
    x_test
  ) <- test_ids

  meta <- s[
    match(
      colnames(x_train),
      s$predictor_key
    ),
    ,
    drop = FALSE
  ]

  assert(
    all(
      meta$predictor_key ==
        colnames(x_train)
    ),
    paste(
      "Metadata alignment failed:",
      family
    )
  )

  block_levels <- unique(
    as.character(
      meta$program_view_block
    )
  )

  grouping <- match(
    as.character(
      meta$program_view_block
    ),
    block_levels
  )

  names(
    grouping
  ) <- meta$predictor_key

  block_map <- data.frame(
    family = family,
    group_id =
      seq_along(
        block_levels
      ),
    program_view_block =
      block_levels,
    stringsAsFactors = FALSE
  )

  block_meta <- unique(
    meta[
      ,
      c(
        "program_view_block",
        "program_id",
        "program_short",
        "view_code",
        "phenotypic_sides_in_program"
      ),
      drop = FALSE
    ]
  )

  block_map <- merge(
    block_map,
    block_meta,
    by = "program_view_block",
    all.x = TRUE,
    sort = FALSE
  )

  block_map <- block_map[
    match(
      block_levels,
      block_map$program_view_block
    ),
    ,
    drop = FALSE
  ]

  list(
    x_train = x_train,
    x_test = x_test,
    meta = meta,
    grouping = grouping,
    block_map = block_map,
    audit =
      rbind_fill(
        audit
      )
  )
}

family_matrix <- lapply(
  FAMILIES,
  function(family) {
    message(
      "[FEATURE MATRIX] ",
      family
    )

    build_model_matrix(
      family,
      family_input[[family]]
    )
  }
)

names(
  family_matrix
) <- FAMILIES

technical_audit <- rbind_fill(
  lapply(
    family_matrix,
    function(z) {
      z$audit
    }
  )
)

block_map_master <- rbind_fill(
  lapply(
    family_matrix,
    function(z) {
      z$block_map
    }
  )
)

write_csv(
  technical_audit,
  file.path(
    DIR_TABLES,
    "02_TECHNICAL_FEATURE_ELIGIBILITY.csv"
  )
)

write_csv(
  block_map_master,
  file.path(
    DIR_TABLES,
    "03_PROGRAM_VIEW_GROUP_MAP.csv"
  )
)

# ==============================================================================
# 7. CV + MSGL HELPERS
# ==============================================================================

make_sample_weights <- function(y) {
  n <- length(y)

  if (!USE_CLASS_WEIGHTS) {
    return(
      rep(
        1 / n,
        n
      )
    )
  }

  y <- factor(
    y,
    levels = CLASS_LEVELS
  )

  tab <- table(y)

  assert(
    all(
      tab > 0L
    ),
    "All classes must be present."
  )

  raw <- 1 /
    as.numeric(
      tab[
        as.character(y)
      ]
    )

  raw /
    sum(raw)
}

make_cv_indices <- function(
    y,
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

  for (cls in CLASS_LEVELS) {
    idx <- which(
      y == cls
    )

    assert(
      length(idx) >= N_FOLDS,
      paste(
        "Class",
        cls,
        "has fewer observations than folds."
      )
    )

    idx <- sample(
      idx,
      length(idx),
      replace = FALSE
    )

    foldid[idx] <- rep(
      seq_len(N_FOLDS),
      length.out =
        length(idx)
    )
  }

  lapply(
    seq_len(N_FOLDS),
    function(fold) {
      which(
        foldid == fold
      )
    }
  )
}

feature_count_from_cv <- function(cvfit) {
  fs <- sglOptim::features_stat(
    cvfit
  )

  if (is.null(dim(fs))) {
    return(
      as.numeric(fs)
    )
  }

  colMeans(
    as.matrix(fs)
  )
}

extract_msgl_probability <- function(
    response_item,
    sample_ids = NULL
) {
  p <- as.matrix(
    response_item
  )

  if (
    nrow(p) ==
      length(CLASS_LEVELS)
  ) {
    # msgl documentation: K x N.
    rn <- rownames(p)

    if (
      !is.null(rn) &&
      all(
        CLASS_LEVELS %in%
          rn
      )
    ) {
      p <- p[
        CLASS_LEVELS,
        ,
        drop = FALSE
      ]
    }

    p <- t(p)
  } else if (
    ncol(p) ==
      length(CLASS_LEVELS)
  ) {
    cn <- colnames(p)

    if (
      !is.null(cn) &&
      all(
        CLASS_LEVELS %in%
          cn
      )
    ) {
      p <- p[
        ,
        CLASS_LEVELS,
        drop = FALSE
      ]
    }
  } else {
    stop(
      "Could not identify class dimension in msgl probabilities.",
      call. = FALSE
    )
  }

  colnames(p) <- CLASS_LEVELS

  if (
    !is.null(sample_ids) &&
    nrow(p) ==
      length(sample_ids)
  ) {
    rownames(p) <- sample_ids
  }

  rs <- rowSums(p)

  assert(
    all(
      is.finite(rs) &
        rs > 0
    ),
    "Invalid msgl probability rows."
  )

  p / rs
}

selected_feature_indices <- function(
    fit,
    model_index,
    p,
    feature_names = NULL
) {
  ff <- sglOptim::features(
    fit
  )[[model_index]]

  if (is.logical(ff)) {
    idx <- which(ff)
  } else if (is.numeric(ff)) {
    idx <- as.integer(ff)
  } else if (is.character(ff)) {
    suppressWarnings(
      idx_num <- as.integer(ff)
    )

    if (
      all(
        is.finite(idx_num)
      )
    ) {
      idx <- idx_num
    } else {
      assert(
        !is.null(feature_names),
        "Character feature identifiers returned by sglOptim::features(), but feature_names were not supplied."
      )

      idx <- match(
        ff,
        feature_names
      )

      idx <- idx[
        is.finite(idx)
      ]
    }
  } else {
    idx <- integer()
  }

  idx <- unique(
    idx[
      idx >= 1L &
        idx <= p
    ]
  )

  sort(idx)
}

fit_msgl_path <- function(
    x,
    y,
    grouping,
    alpha
) {
  msgl::fit(
    x = x,
    classes = factor(
      y,
      levels = CLASS_LEVELS
    ),
    sampleWeights =
      make_sample_weights(y),
    grouping = grouping,
    alpha = alpha,
    standardize = TRUE,
    lambda =
      LAMBDA_MIN_REL,
    d =
      PATH_LENGTH,
    intercept = TRUE
  )
}

# ==============================================================================
# 8. REPEATED 3-FOLD TUNING OF ALPHA + PATH INDEX
# ==============================================================================

run_repeated_cv_grid <- function(
    family,
    x,
    y,
    grouping,
    seed
) {
  rows <- list()
  predictions <- list()

  rr <- 1L
  pp <- 1L

  for (alpha in ALPHA_GRID) {
    message(
      "  alpha=",
      alpha
    )

    for (r in seq_len(CV_REPEATS)) {
      cv_indices <- make_cv_indices(
        y,
        seed +
          round(
            alpha *
              10000
          ) +
          1000L *
            r
      )

      cvfit <- try(
        msgl::cv(
          x = x,
          classes =
            factor(
              y,
              levels =
                CLASS_LEVELS
            ),
          sampleWeights =
            make_sample_weights(y),
          grouping =
            grouping,
          alpha =
            alpha,
          standardize =
            TRUE,
          lambda =
            LAMBDA_MIN_REL,
          d =
            PATH_LENGTH,
          fold =
            N_FOLDS,
          cv.indices =
            cv_indices,
          intercept =
            TRUE,
          use_parallel =
            FALSE
        ),
        silent = TRUE
      )

      if (
        inherits(
          cvfit,
          "try-error"
        )
      ) {
        warning(
          "msgl::cv failed: ",
          family,
          " alpha=",
          alpha,
          " repeat=",
          r,
          " | ",
          as.character(
            cvfit
          )
        )

        next
      }

      # Err/features_stat/nmod/features are sglOptim generics.
      # msgl registers S3 methods for class 'msgl' but does not export
      # those generic names from the msgl namespace.
      err <- as.numeric(
        sglOptim::Err(
          cvfit,
          type = "loglike"
        )
      )

      feat_n <- feature_count_from_cv(
        cvfit
      )

      n_models <- min(
        length(err),
        length(feat_n),
        length(
          cvfit$response
        )
      )

      assert(
        n_models > 0L,
        "msgl CV returned no path models."
      )

      err <- err[
        seq_len(
          n_models
        )
      ]

      feat_n <- feat_n[
        seq_len(
          n_models
        )
      ]

      for (j in seq_len(n_models)) {
        rows[[rr]] <- data.frame(
          family = family,
          alpha = alpha,
          repeat_id = r,
          model_index = j,
          cv_loglike_error =
            err[j],
          expected_n_features =
            feat_n[j],
          stringsAsFactors = FALSE
        )

        rr <- rr + 1L
      }

      # Keep OOF probabilities for every model; final aggregation only uses
      # the model selected after the repeated-CV table is summarized.
      predictions[[pp]] <- list(
        alpha = alpha,
        repeat_id = r,
        response =
          cvfit$response
      )

      pp <- pp + 1L
    }
  }

  tab <- rbind_fill(
    rows
  )

  assert(
    nrow(tab) > 0L,
    paste(
      "All msgl CV fits failed for",
      family
    )
  )

  agg_split <- split(
    tab,
    paste(
      tab$alpha,
      tab$model_index,
      sep = "\r"
    )
  )

  summary_rows <- lapply(
    agg_split,
    function(z) {
      n_ok <- sum(
        is.finite(
          z$cv_loglike_error
        )
      )

      data.frame(
        family = family,
        alpha =
          z$alpha[1L],
        model_index =
          z$model_index[1L],
        n_repeats =
          n_ok,
        mean_cv_loglike =
          mean(
            z$cv_loglike_error,
            na.rm = TRUE
          ),
        sd_cv_loglike =
          stats::sd(
            z$cv_loglike_error,
            na.rm = TRUE
          ),
        se_cv_loglike =
          stats::sd(
            z$cv_loglike_error,
            na.rm = TRUE
          ) /
          sqrt(
            n_ok
          ),
        mean_expected_n_features =
          mean(
            z$expected_n_features,
            na.rm = TRUE
          ),
        stringsAsFactors = FALSE
      )
    }
  )

  summary <- do.call(
    rbind,
    summary_rows
  )

  summary <- summary[
    is.finite(
      summary$mean_cv_loglike
    ),
    ,
    drop = FALSE
  ]

  assert(
    nrow(summary) > 0L,
    paste(
      "No valid tuning summary:",
      family
    )
  )

  min_idx <- which.min(
    summary$mean_cv_loglike
  )

  min_error <-
    summary$mean_cv_loglike[
      min_idx
    ]

  min_se <-
    summary$se_cv_loglike[
      min_idx
    ]

  if (!is.finite(min_se)) {
    min_se <- 0
  }

  threshold <- min_error +
    min_se

  candidates <- summary[
    summary$mean_cv_loglike <=
      threshold &
      summary$mean_expected_n_features >
      0,
    ,
    drop = FALSE
  ]

  if (!nrow(candidates)) {
    candidates <- summary[
      summary$mean_expected_n_features >
        0,
      ,
      drop = FALSE
    ]

    candidates <- candidates[
      order(
        candidates$mean_cv_loglike,
        candidates$mean_expected_n_features
      ),
      ,
      drop = FALSE
    ]

    selection_rule <-
      "MIN_LOGLOSS_NONZERO_FALLBACK"
  } else {
    # 1-SE-like parsimony rule:
    # among statistically near-equivalent configurations, prefer the one
    # with the fewest expected selected original features.
    candidates <- candidates[
      order(
        candidates$mean_expected_n_features,
        candidates$mean_cv_loglike,
        abs(
          candidates$alpha -
            0.5
        ),
        candidates$model_index
      ),
      ,
      drop = FALSE
    ]

    selection_rule <-
      "REPEATED_CV_1SE_LIKE_MIN_FEATURES"
  }

  chosen <- candidates[
    1L,
    ,
    drop = FALSE
  ]

  summary$within_1se_like <-
    summary$mean_cv_loglike <=
      threshold

  summary$selected_config <-
    summary$alpha ==
      chosen$alpha &
      summary$model_index ==
        chosen$model_index

  list(
    raw = tab,
    summary = summary,
    predictions = predictions,
    chosen = chosen,
    min_error = min_error,
    min_se = min_se,
    one_se_threshold = threshold,
    selection_rule = selection_rule
  )
}

aggregate_chosen_oof <- function(
    cv_result,
    chosen_alpha,
    chosen_index,
    sample_ids
) {
  selected <- cv_result$predictions[
    vapply(
      cv_result$predictions,
      function(z) {
        isTRUE(
          all.equal(
            z$alpha,
            chosen_alpha
          )
        ) &&
          length(
            z$response
          ) >=
          chosen_index
      },
      logical(1)
    )
  ]

  assert(
    length(selected) > 0L,
    "No OOF predictions available for chosen sparse-group configuration."
  )

  prob_sum <- matrix(
    0,
    nrow =
      length(
        sample_ids
      ),
    ncol =
      length(
        CLASS_LEVELS
      ),
    dimnames = list(
      sample_ids,
      CLASS_LEVELS
    )
  )

  prob_n <- integer(
    length(
      sample_ids
    )
  )

  for (z in selected) {
    p <- extract_msgl_probability(
      z$response[[chosen_index]],
      sample_ids =
        sample_ids
    )

    assert(
      nrow(p) ==
        length(
          sample_ids
        ),
      "OOF response dimension mismatch."
    )

    prob_sum <- prob_sum +
      p

    prob_n <- prob_n + 1L
  }

  prob_sum /
    prob_n
}

# ==============================================================================
# 9. STABILITY AUDIT AT THE CHOSEN STRUCTURED PENALTY
# ==============================================================================

run_stability_audit <- function(
    family,
    x,
    y,
    meta,
    grouping,
    alpha,
    model_index,
    seed
) {
  feature_count <- setNames(
    integer(
      ncol(x)
    ),
    colnames(x)
  )

  blocks <- unique(
    as.character(
      meta$program_view_block
    )
  )

  block_count <- setNames(
    integer(
      length(blocks)
    ),
    blocks
  )

  audit <- list()
  aa <- 1L
  n_models <- 0L

  for (r in seq_len(STABILITY_REPEATS)) {
    cv_indices <- make_cv_indices(
      y,
      seed +
        800000L +
        1000L *
          r
    )

    # cv_indices are validation folds; fit on complementary folds.
    for (fold in seq_len(N_FOLDS)) {
      idx_valid <-
        cv_indices[[fold]]

      idx_train <- setdiff(
        seq_len(
          nrow(x)
        ),
        idx_valid
      )

      y_sub <- factor(
        y[idx_train],
        levels =
          CLASS_LEVELS
      )

      assert(
        all(
          table(y_sub) > 0L
        ),
        "A fold-training subset lost a class."
      )

      fit <- try(
        fit_msgl_path(
          x = x[
            idx_train,
            ,
            drop = FALSE
          ],
          y = y_sub,
          grouping = grouping,
          alpha = alpha
        ),
        silent = TRUE
      )

      if (
        inherits(
          fit,
          "try-error"
        )
      ) {
        warning(
          "Stability fit failed: ",
          family,
          " repeat=",
          r,
          " fold=",
          fold
        )

        next
      }

      actual_index <- min(
        model_index,
        sglOptim::nmod(
          fit
        )
      )

      selected_idx <- selected_feature_indices(
        fit,
        actual_index,
        p =
          ncol(x),
        feature_names =
          colnames(x)
      )

      selected_keys <- colnames(x)[
        selected_idx
      ]

      if (length(selected_keys)) {
        feature_count[
          selected_keys
        ] <- feature_count[
          selected_keys
        ] +
          1L

        selected_blocks <- unique(
          meta$program_view_block[
            match(
              selected_keys,
              meta$predictor_key
            )
          ]
        )

        block_count[
          selected_blocks
        ] <- block_count[
          selected_blocks
        ] +
          1L
      } else {
        selected_blocks <-
          character()
      }

      n_models <- n_models + 1L

      audit[[aa]] <- data.frame(
        family = family,
        repeat_id = r,
        fold = fold,
        alpha = alpha,
        requested_model_index =
          model_index,
        actual_model_index =
          actual_index,
        n_selected_features =
          length(
            selected_keys
          ),
        n_selected_blocks =
          length(
            selected_blocks
          ),
        stringsAsFactors = FALSE
      )

      aa <- aa + 1L
    }
  }

  assert(
    n_models > 0L,
    paste(
      "No successful stability fits:",
      family
    )
  )

  feature_stability <- data.frame(
    family = family,
    predictor_key =
      names(
        feature_count
      ),
    feature_selection_count =
      as.integer(
        feature_count
      ),
    n_fold_models =
      n_models,
    feature_stability =
      as.numeric(
        feature_count
      ) /
      n_models,
    stringsAsFactors = FALSE
  )

  feature_stability <- merge(
    feature_stability,
    meta,
    by = c(
      "family",
      "predictor_key"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  block_stability <- data.frame(
    family = family,
    program_view_block =
      names(
        block_count
      ),
    block_selection_count =
      as.integer(
        block_count
      ),
    n_fold_models =
      n_models,
    block_stability =
      as.numeric(
        block_count
      ) /
      n_models,
    stringsAsFactors = FALSE
  )

  block_meta <- unique(
    meta[
      ,
      c(
        "family",
        "program_view_block",
        "program_id",
        "program_short",
        "view_code",
        "phenotypic_sides_in_program"
      ),
      drop = FALSE
    ]
  )

  block_stability <- merge(
    block_stability,
    block_meta,
    by = c(
      "family",
      "program_view_block"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  feature_stability$audit_stability_tier <-
    ifelse(
      feature_stability$feature_stability >=
        STABILITY_AUDIT_THRESHOLD,
      "HIGH_STABILITY_AUDIT",
      "LOWER_STABILITY_AUDIT"
    )

  block_stability$audit_stability_tier <-
    ifelse(
      block_stability$block_stability >=
        STABILITY_AUDIT_THRESHOLD,
      "HIGH_STABILITY_AUDIT",
      "LOWER_STABILITY_AUDIT"
    )

  list(
    feature = feature_stability,
    block = block_stability,
    audit =
      rbind_fill(
        audit
      )
  )
}

# ==============================================================================
# 10. METRICS
# ==============================================================================

multiclass_logloss <- function(
    truth,
    prob
) {
  truth <- factor(
    truth,
    levels = CLASS_LEVELS
  )

  idx <- match(
    as.character(truth),
    CLASS_LEVELS
  )

  -mean(
    log(
      pmax(
        prob[
          cbind(
            seq_len(
              nrow(prob)
            ),
            idx
          )
        ],
        1e-15
      )
    )
  )
}

multiclass_brier <- function(
    truth,
    prob
) {
  truth <- factor(
    truth,
    levels = CLASS_LEVELS
  )

  Y <- matrix(
    0,
    nrow =
      nrow(prob),
    ncol =
      length(
        CLASS_LEVELS
      ),
    dimnames = list(
      NULL,
      CLASS_LEVELS
    )
  )

  idx <- match(
    as.character(truth),
    CLASS_LEVELS
  )

  Y[
    cbind(
      seq_len(
        nrow(Y)
      ),
      idx
    )
  ] <- 1

  mean(
    rowSums(
      (
        prob -
          Y
      )^2
    )
  )
}

evaluate_model <- function(
    truth,
    prob,
    family,
    split
) {
  truth <- factor(
    truth,
    levels = CLASS_LEVELS
  )

  pred <- factor(
    CLASS_LEVELS[
      max.col(
        prob,
        ties.method = "first"
      )
    ],
    levels =
      CLASS_LEVELS
  )

  cm <- table(
    truth = truth,
    predicted = pred
  )

  per_class <- do.call(
    rbind,
    lapply(
      CLASS_LEVELS,
      function(cls) {
        tp <- cm[
          cls,
          cls
        ]

        fn <- sum(
          cm[
            cls,
            ]
        ) -
          tp

        fp <- sum(
          cm[
            ,
            cls
          ]
        ) -
          tp

        recall <- if (
          tp + fn > 0
        ) {
          tp /
            (
              tp +
                fn
            )
        } else {
          NA_real_
        }

        precision <- if (
          tp + fp > 0
        ) {
          tp /
            (
              tp +
                fp
            )
        } else {
          0
        }

        f1 <- if (
          is.finite(recall) &&
          is.finite(precision) &&
          recall +
            precision >
            0
        ) {
          2 *
            recall *
            precision /
            (
              recall +
                precision
            )
        } else if (
          is.finite(recall) &&
          recall == 0
        ) {
          0
        } else {
          NA_real_
        }

        data.frame(
          family = family,
          split = split,
          class = cls,
          recall = recall,
          precision = precision,
          f1 = f1,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  assert(
    all(
      is.finite(
        per_class$recall
      )
    ),
    "Undefined class recall."
  )

  assert(
    all(
      is.finite(
        per_class$f1
      )
    ),
    "Undefined class F1."
  )

  accuracy <- mean(
    pred ==
      truth
  )

  bt <- stats::binom.test(
    sum(
      pred ==
        truth
    ),
    length(truth)
  )

  brier <- multiclass_brier(
    truth,
    prob
  )

  list(
    summary = data.frame(
      family = family,
      split = split,
      n = length(truth),
      accuracy = accuracy,
      accuracy_ci95_low =
        bt$conf.int[1L],
      accuracy_ci95_high =
        bt$conf.int[2L],
      balanced_accuracy =
        mean(
          per_class$recall
        ),
      macro_f1 =
        mean(
          per_class$f1
        ),
      log_loss =
        multiclass_logloss(
          truth,
          prob
        ),
      multiclass_brier_unscaled =
        brier,
      multiclass_brier_scaled_0_1 =
        brier /
        2,
      mean_max_probability =
        mean(
          apply(
            prob,
            1L,
            max
          )
        ),
      stringsAsFactors = FALSE
    ),
    per_class = per_class,
    confusion = cm,
    predicted = pred,
    prob = prob
  )
}

# ==============================================================================
# 11. RUN P1 / P2 / P3
# ==============================================================================

tuning_raw_all <- list()
tuning_summary_all <- list()
panel_all <- list()
block_selected_all <- list()
feature_stability_all <- list()
block_stability_all <- list()
stability_audit_all <- list()
metrics_all <- list()
per_class_all <- list()
predictions_all <- list()

for (family in FAMILIES) {
  message("")
  message("============================================================")
  message("[SPARSE GROUP LASSO] ", family)
  message("============================================================")

  fm <- family_matrix[[family]]
  inp <- family_input[[family]]

  x_train <- fm$x_train
  x_test <- fm$x_test
  y_train <- inp$train_meta$grupo
  y_test <- inp$test_meta$grupo
  meta <- fm$meta
  grouping <- fm$grouping

  message(
    "TRAIN n=",
    nrow(x_train),
    " | biological features=",
    ncol(x_train),
    " | program×omic groups=",
    length(
      unique(
        grouping
      )
    )
  )

  cv_res <- run_repeated_cv_grid(
    family = family,
    x = x_train,
    y = y_train,
    grouping = grouping,
    seed =
      SEED +
      1000000L *
        match(
          family,
          FAMILIES
        )
  )

  chosen_alpha <- as.numeric(
    cv_res$chosen$alpha[1L]
  )

  chosen_index <- as.integer(
    cv_res$chosen$model_index[1L]
  )

  message(
    "Chosen alpha=",
    chosen_alpha,
    " | path index=",
    chosen_index,
    " | mean expected features=",
    sprintf(
      "%.2f",
      cv_res$chosen$mean_expected_n_features[1L]
    ),
    " | rule=",
    cv_res$selection_rule
  )

  # Full TRAIN model using the exact same path definition.
  final_fit <- fit_msgl_path(
    x = x_train,
    y = y_train,
    grouping = grouping,
    alpha = chosen_alpha
  )

  actual_index <- min(
    chosen_index,
    sglOptim::nmod(
      final_fit
    )
  )

  selected_idx <- selected_feature_indices(
    final_fit,
    actual_index,
    p =
      ncol(x_train),
    feature_names =
      colnames(x_train)
  )

  assert(
    length(selected_idx) > 0L,
    paste(
      "Chosen sparse-group model selected zero features:",
      family
    )
  )

  selected_keys <- colnames(
    x_train
  )[
    selected_idx
  ]

  selected_blocks <- unique(
    meta$program_view_block[
      match(
        selected_keys,
        meta$predictor_key
      )
    ]
  )

  message(
    "Final full-TRAIN panel=",
    length(selected_keys),
    " features | ",
    length(selected_blocks),
    " program×omic groups"
  )

  # OOF probabilities from the repeated CVs used for conditional tuning.
  train_prob <- aggregate_chosen_oof(
    cv_result = cv_res,
    chosen_alpha = chosen_alpha,
    chosen_index = chosen_index,
    sample_ids =
      rownames(
        x_train
      )
  )

  # Final TEST projection.
  pred_test <- stats::predict(
    final_fit,
    x_test
  )

  test_prob <- extract_msgl_probability(
    pred_test$response[[actual_index]],
    sample_ids =
      rownames(
        x_test
      )
  )

  stability <- run_stability_audit(
    family = family,
    x = x_train,
    y = y_train,
    meta = meta,
    grouping = grouping,
    alpha = chosen_alpha,
    model_index = chosen_index,
    seed =
      SEED +
      2000000L *
        match(
          family,
          FAMILIES
        )
  )

  panel <- meta[
    match(
      selected_keys,
      meta$predictor_key
    ),
    ,
    drop = FALSE
  ]

  panel$selected_in_full_train <- TRUE

  panel <- merge(
    panel,
    stability$feature[
      ,
      c(
        "family",
        "predictor_key",
        "feature_selection_count",
        "n_fold_models",
        "feature_stability",
        "audit_stability_tier"
      ),
      drop = FALSE
    ],
    by = c(
      "family",
      "predictor_key"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  panel <- merge(
    panel,
    stability$block[
      ,
      c(
        "family",
        "program_view_block",
        "block_selection_count",
        "block_stability",
        "audit_stability_tier"
      ),
      drop = FALSE
    ],
    by = c(
      "family",
      "program_view_block"
    ),
    all.x = TRUE,
    sort = FALSE,
    suffixes = c(
      "_feature",
      "_block"
    )
  )

  panel$chosen_alpha <- chosen_alpha
  panel$chosen_path_index <- chosen_index
  panel$selection_rule <- cv_res$selection_rule

  panel <- panel[
    order(
      match(
        panel$view_code,
        VIEW_LEVELS
      ),
      -panel$block_stability,
      -panel$feature_stability,
      panel$program_id,
      panel$feature_label
    ),
    ,
    drop = FALSE
  ]

  blocks_selected <- stability$block[
    stability$block$program_view_block %in%
      selected_blocks,
    ,
    drop = FALSE
  ]

  blocks_selected$selected_in_full_train <- TRUE
  blocks_selected$chosen_alpha <- chosen_alpha
  blocks_selected$chosen_path_index <- chosen_index

  ev_train <- evaluate_model(
    y_train,
    train_prob,
    family,
    "TRAIN_REPEATED_CV_TUNING_CONDITIONAL"
  )

  ev_test <- evaluate_model(
    y_test,
    test_prob,
    family,
    "TEST_EXPLORATORY_REUSED"
  )

  metrics <- rbind(
    ev_train$summary,
    ev_test$summary
  )

  metrics$chosen_alpha <- chosen_alpha
  metrics$chosen_path_index <- chosen_index
  metrics$n_panel_features <-
    length(selected_keys)
  metrics$n_selected_program_view_groups <-
    length(selected_blocks)
  metrics$selection_rule <-
    cv_res$selection_rule

  predictions <- rbind(
    data.frame(
      family = family,
      sample_id =
        rownames(
          train_prob
        ),
      split =
        "TRAIN_REPEATED_CV_TUNING_CONDITIONAL",
      truth =
        as.character(
          y_train
        ),
      predicted =
        as.character(
          ev_train$predicted
        ),
      NP =
        train_prob[, "NP"],
      FD =
        train_prob[, "FD"],
      SO =
        train_prob[, "SO"],
      stringsAsFactors = FALSE
    ),
    data.frame(
      family = family,
      sample_id =
        rownames(
          test_prob
        ),
      split =
        "TEST_EXPLORATORY_REUSED",
      truth =
        as.character(
          y_test
        ),
      predicted =
        as.character(
          ev_test$predicted
        ),
      NP =
        test_prob[, "NP"],
      FD =
        test_prob[, "FD"],
      SO =
        test_prob[, "SO"],
      stringsAsFactors = FALSE
    )
  )

  fam_dir <- file.path(
    DIR_TABLES,
    family
  )

  dir.create(
    fam_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  write_csv(
    cv_res$raw,
    file.path(
      fam_dir,
      paste0(
        family,
        "_01_REPEATED_CV_RAW.csv"
      )
    )
  )

  write_csv(
    cv_res$summary,
    file.path(
      fam_dir,
      paste0(
        family,
        "_02_ALPHA_LAMBDA_PATH_TUNING.csv"
      )
    )
  )

  write_csv(
    stability$feature,
    file.path(
      fam_dir,
      paste0(
        family,
        "_03_FEATURE_STABILITY_AUDIT.csv"
      )
    )
  )

  write_csv(
    stability$block,
    file.path(
      fam_dir,
      paste0(
        family,
        "_04_BLOCK_STABILITY_AUDIT.csv"
      )
    )
  )

  write_csv(
    panel,
    file.path(
      fam_dir,
      paste0(
        family,
        "_05_FINAL_SPARSE_GROUP_PANEL.csv"
      )
    )
  )

  write_csv(
    blocks_selected,
    file.path(
      fam_dir,
      paste0(
        family,
        "_06_SELECTED_PROGRAM_VIEW_GROUPS.csv"
      )
    )
  )

  write_csv(
    metrics,
    file.path(
      fam_dir,
      paste0(
        family,
        "_07_METRICS.csv"
      )
    )
  )

  write_csv(
    rbind(
      ev_train$per_class,
      ev_test$per_class
    ),
    file.path(
      fam_dir,
      paste0(
        family,
        "_08_PER_CLASS_METRICS.csv"
      )
    )
  )

  write_csv(
    predictions,
    file.path(
      fam_dir,
      paste0(
        family,
        "_09_PREDICTIONS.csv"
      )
    )
  )

  write_csv(
    as.data.frame.matrix(
      ev_train$confusion
    ),
    file.path(
      fam_dir,
      paste0(
        family,
        "_10_TRAIN_CONFUSION.csv"
      )
    )
  )

  write_csv(
    as.data.frame.matrix(
      ev_test$confusion
    ),
    file.path(
      fam_dir,
      paste0(
        family,
        "_11_TEST_EXPLORATORY_CONFUSION.csv"
      )
    )
  )

  saveRDS(
    list(
      family = family,
      chosen_alpha =
        chosen_alpha,
      chosen_path_index =
        chosen_index,
      actual_final_model_index =
        actual_index,
      selection_rule =
        cv_res$selection_rule,
      selected_features =
        selected_keys,
      selected_blocks =
        selected_blocks,
      panel =
        panel,
      final_fit =
        final_fit,
      source_panel10 =
        PANEL10
    ),
    file.path(
      DIR_MODELS,
      paste0(
        family,
        "_SPARSE_GROUP_LASSO_FINAL.rds"
      )
    )
  )

  tuning_raw_all[[family]] <-
    cv_res$raw

  tuning_summary_all[[family]] <-
    cv_res$summary

  panel_all[[family]] <-
    panel

  block_selected_all[[family]] <-
    blocks_selected

  feature_stability_all[[family]] <-
    stability$feature

  block_stability_all[[family]] <-
    stability$block

  stability_audit_all[[family]] <-
    stability$audit

  metrics_all[[family]] <-
    metrics

  per_class_all[[family]] <-
    rbind(
      ev_train$per_class,
      ev_test$per_class
    )

  predictions_all[[family]] <-
    predictions

  # Panel composition plot.
  counts <- as.data.frame(
    table(
      toupper(
        panel$view_code
      )
    ),
    stringsAsFactors = FALSE
  )

  names(counts) <- c(
    "view",
    "n_features"
  )

  p_comp <- ggplot2::ggplot(
    counts,
    ggplot2::aes(
      x = view,
      y = n_features
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title =
        paste0(
          family,
          " — sparse-group final panel"
        ),
      subtitle =
        paste0(
          "alpha=",
          chosen_alpha,
          "; ",
          length(selected_blocks),
          " selected program×omic groups"
        ),
      x = "Omic view",
      y = "Selected original measurable features"
    ) +
    theme_panel(11)

  ggplot2::ggsave(
    file.path(
      DIR_FIGURES,
      paste0(
        family,
        "_FINAL_PANEL_COMPOSITION.png"
      )
    ),
    p_comp,
    width = 7,
    height = 5,
    dpi = DPI
  )
}

# ==============================================================================
# 12. MASTER OUTPUTS
# ==============================================================================

tuning_raw_master <- rbind_fill(
  tuning_raw_all
)

tuning_summary_master <- rbind_fill(
  tuning_summary_all
)

panel_master <- rbind_fill(
  panel_all
)

selected_blocks_master <- rbind_fill(
  block_selected_all
)

feature_stability_master <- rbind_fill(
  feature_stability_all
)

block_stability_master <- rbind_fill(
  block_stability_all
)

stability_audit_master <- rbind_fill(
  stability_audit_all
)

metrics_master <- rbind_fill(
  metrics_all
)

per_class_master <- rbind_fill(
  per_class_all
)

predictions_master <- rbind_fill(
  predictions_all
)

write_csv(
  tuning_raw_master,
  file.path(
    DIR_TABLES,
    "04_REPEATED_CV_RAW_ALL_FAMILIES.csv"
  )
)

write_csv(
  tuning_summary_master,
  file.path(
    DIR_TABLES,
    "05_ALPHA_LAMBDA_PATH_TUNING_ALL_FAMILIES.csv"
  )
)

write_csv(
  feature_stability_master,
  file.path(
    DIR_TABLES,
    "06_FEATURE_STABILITY_AUDIT_ALL_FAMILIES.csv"
  )
)

write_csv(
  block_stability_master,
  file.path(
    DIR_TABLES,
    "07_BLOCK_STABILITY_AUDIT_ALL_FAMILIES.csv"
  )
)

write_csv(
  panel_master,
  file.path(
    DIR_TABLES,
    "08_FINAL_SPARSE_GROUP_PANEL_ALL_FAMILIES.csv"
  )
)

write_csv(
  selected_blocks_master,
  file.path(
    DIR_TABLES,
    "09_SELECTED_PROGRAM_VIEW_GROUPS_ALL_FAMILIES.csv"
  )
)

write_csv(
  metrics_master,
  file.path(
    DIR_TABLES,
    "10_METRICS_ALL_FAMILIES.csv"
  )
)

write_csv(
  per_class_master,
  file.path(
    DIR_TABLES,
    "11_PER_CLASS_METRICS_ALL_FAMILIES.csv"
  )
)

write_csv(
  predictions_master,
  file.path(
    DIR_TABLES,
    "12_PREDICTIONS_ALL_FAMILIES.csv"
  )
)

write_csv(
  stability_audit_master,
  file.path(
    DIR_TABLES,
    "13_STABILITY_RESAMPLE_AUDIT_ALL_FAMILIES.csv"
  )
)

# ==============================================================================
# 13. CROSS-ANALYSIS CONSENSUS
# ==============================================================================

consensus_split <- split(
  panel_master,
  paste(
    panel_master$view_code,
    panel_master$feature_original,
    sep = "\r"
  )
)

consensus <- do.call(
  rbind,
  lapply(
    consensus_split,
    function(z) {
      families <- sort(
        unique(
          as.character(
            z$family
          )
        )
      )

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
        n_families =
          length(
            families
          ),
        families =
          paste(
            families,
            collapse = ";"
          ),
        n_distinct_program_view_blocks =
          length(
            unique(
              z$program_view_block
            )
          ),
        mean_feature_stability =
          mean_or_na(
            z$feature_stability
          ),
        mean_block_stability =
          mean_or_na(
            z$block_stability
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(
  consensus
) <- NULL

consensus <- consensus[
  order(
    -consensus$n_families,
    -consensus$mean_feature_stability,
    -consensus$mean_block_stability,
    match(
      consensus$view_code,
      VIEW_LEVELS
    ),
    consensus$feature_label
  ),
  ,
  drop = FALSE
]

write_csv(
  consensus,
  file.path(
    DIR_TABLES,
    "14_CROSS_ANALYSIS_SPARSE_GROUP_CONSENSUS.csv"
  )
)

# ==============================================================================
# 14. MANIFEST + README
# ==============================================================================

manifest <- data.frame(
  item = c(
    "BIO_PROGRAM_RUN",
    "PANEL10_RUN",
    "OUTPUT_ROOT",
    "MSGL_VERSION",
    "GROUPING",
    "ALPHA_GRID",
    "LAMBDA_MIN_REL",
    "PATH_LENGTH",
    "CV_REPEATS",
    "STABILITY_REPEATS",
    "SELECTION_RULE",
    "OMIC_QUOTA",
    "PCA_SELECTION",
    "TEST_STATUS",
    "VALIDATION_LIMITATION"
  ),
  value = c(
    BIO_RUN,
    PANEL10,
    OUT_ROOT,
    as.character(
      utils::packageVersion(
        "msgl"
      )
    ),
    "PRIMARY_PROGRAM_X_OMIC",
    paste(
      ALPHA_GRID,
      collapse = ","
    ),
    as.character(
      LAMBDA_MIN_REL
    ),
    as.character(
      PATH_LENGTH
    ),
    as.character(
      CV_REPEATS
    ),
    as.character(
      STABILITY_REPEATS
    ),
    "Repeated stratified 3-fold CV; 1-SE-like near-optimal set; minimum expected number of original selected features.",
    "NONE",
    "NO",
    "TEST_EXPLORATORY_REUSED",
    paste0(
      "TRAIN tuning is conditional on the biological candidate universe and ",
      "primary-program assignment derived from full TRAIN development data. ",
      "FD has only 3 TRAIN samples, preventing a reliable fully nested ",
      "multinomial design. TEST was inspected during prior method development."
    )
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
  "METATEST — MULTINOMIAL SPARSE GROUP LASSO",
  "",
  "QUESTION",
  "Can the biological candidate universe yield a compact supervised panel when",
  "program×omic structure is included directly in the penalty?",
  "",
  "STRUCTURE",
  "  script-10 biological candidate",
  "    -> one outcome-free primary biological program",
  "    -> PROGRAM x OMIC group",
  "    -> sparse group lasso",
  "       * group penalty selects/deselects program×omic blocks",
  "       * L1 penalty selects individual features inside selected blocks",
  "",
  "NO PCA is used for feature selection.",
  "NO omic quota is imposed.",
  "NO fixed final panel size is imposed.",
  "",
  "TUNING",
  paste0(
    "  alpha grid: ",
    paste(
      ALPHA_GRID,
      collapse = ", "
    )
  ),
  paste0(
    "  path length: ",
    PATH_LENGTH
  ),
  paste0(
    "  minimum lambda fraction: ",
    LAMBDA_MIN_REL
  ),
  paste0(
    "  repeated stratified 3-fold CV: ",
    CV_REPEATS,
    " repeats"
  ),
  "  choice: 1-SE-like near-optimal configurations -> fewest expected features.",
  "",
  "STABILITY",
  "  Stability is an audit annotation only.",
  "  It does NOT delete a feature chosen by the final full-TRAIN SGL model.",
  "",
  "MAIN OUTPUTS",
  "  tables/05_ALPHA_LAMBDA_PATH_TUNING_ALL_FAMILIES.csv",
  "  tables/06_FEATURE_STABILITY_AUDIT_ALL_FAMILIES.csv",
  "  tables/07_BLOCK_STABILITY_AUDIT_ALL_FAMILIES.csv",
  "  tables/08_FINAL_SPARSE_GROUP_PANEL_ALL_FAMILIES.csv",
  "  tables/09_SELECTED_PROGRAM_VIEW_GROUPS_ALL_FAMILIES.csv",
  "  tables/10_METRICS_ALL_FAMILIES.csv",
  "  tables/11_PER_CLASS_METRICS_ALL_FAMILIES.csv",
  "  tables/14_CROSS_ANALYSIS_SPARSE_GROUP_CONSENSUS.csv",
  "",
  "VALIDATION NOTE",
  "  TRAIN metrics are conditional development estimates.",
  "  TEST is labeled TEST_EXPLORATORY_REUSED because it has already been inspected",
  "  during previous method-development iterations."
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
  "FINAL_SPARSE_GROUP_LASSO_PANEL"
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
message("SPARSE GROUP LASSO COMPLETE")
message("Output: ", OUT_ROOT)
message("")

for (family in FAMILIES) {
  pan <- panel_all[[family]]

  counts <- table(
    toupper(
      pan$view_code
    )
  )

  message(
    family,
    ": ",
    nrow(pan),
    " features | ",
    length(
      unique(
        pan$program_view_block
      )
    ),
    " groups | views: ",
    paste(
      paste0(
        names(counts),
        "=",
        as.integer(counts)
      ),
      collapse = ", "
    )
  )
}

message("============================================================")

#!/usr/bin/env Rscript

# ==============================================================================
# METATEST — 17_evidenceA_plsda_pca_validation.R
#
# LAST ALTERNATIVE AGREED:
#
#   script-10 biological candidates
#        ∩
#   strict Evidence-A pool from the TRAIN-only bidirectional Bayesian framework
#        -> PLS-DA (supervised representation)
#        -> PCA on EXACTLY the same features (unsupervised audit only)
#
# PRIMARY EVIDENCE-A DEFINITION USED HERE
# ---------------------------------------
# The framework already exports:
#
#   06C_HIGH_CONFIDENCE_AA_SAME_PATTERN_TRAIN_ONLY.csv
#
# This pool requires, upstream:
#   * reconstructed evidence_grade == A_strong
#   * observed evidence_grade      == A_strong
#   * same non-uncertain bidirectional pattern
#   * sufficient reconstructed/observed directional concordance
#
# This script DOES NOT redefine Evidence A.
#
# WHY ALL AA PATTERNS, NOT ONLY CLASS-LIKE?
# -----------------------------------------
# PLS-DA may benefit from:
#   NP_like_high / FD_like_high / SO_like_high
#   FD/SO splitters
#   shared FD/SO signals
#
# Therefore the default pool is 06C, not the narrower class-signature 06D.
#
# PLS-DA
# ------
# * Original OBSERVED features from ready_for_modeling.
# * Feature eligibility is frozen before PLS-DA.
# * Number of PLS components chosen in TRAIN by repeated stratified 3-fold CV.
# * Selection rule: smallest ncomp within a 1-SE-like BER region.
# * No VIP threshold determines model membership.
# * VIP>=1 is exported only as a descriptive prioritization annotation.
#
# PCA
# ---
# * PCA uses exactly the same Evidence-A feature matrix.
# * PCA is fit on TRAIN only.
# * TEST is projected with frozen TRAIN center/scale/rotation.
# * PCA does not select or remove features.
#
# TEST
# ----
# TEST has already been inspected during prior method development.
# It is therefore labelled TEST_EXPLORATORY_REUSED.
#
# With TRAIN n=16 and FD=3, repeated CV is a conditional development audit,
# not an independent fully nested validation.
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
    Sys.getenv("PANEL17_FAMILIES", unset = "P1,P2,P3"),
    ",",
    fixed = TRUE
  )[[1]]
)
FAMILIES <- FAMILIES[nzchar(FAMILIES)]

CLASS_LEVELS <- c("NP", "FD", "SO")
VIEW_LEVELS <- c("tx", "pr", "me", "cl")

SEED <- as.integer(
  Sys.getenv("PANEL17_SEED", unset = "20261001")
)

N_FOLDS <- as.integer(
  Sys.getenv("PANEL17_FOLDS", unset = "3")
)

CV_REPEATS <- as.integer(
  Sys.getenv("PANEL17_CV_REPEATS", unset = "50")
)

MAX_PLS_COMPONENTS <- as.integer(
  Sys.getenv("PANEL17_MAX_PLS_COMPONENTS", unset = "3")
)

VIP_REFERENCE_THRESHOLD <- as.numeric(
  Sys.getenv("PANEL17_VIP_REFERENCE_THRESHOLD", unset = "1")
)

LOADING_BAR_TOPN <- as.integer(
  Sys.getenv("PANEL17_LOADING_BAR_TOPN", unset = "20")
)

LOADING_MAP_TOPN <- as.integer(
  Sys.getenv("PANEL17_LOADING_MAP_TOPN", unset = "25")
)

DPI <- as.integer(
  Sys.getenv("PANEL17_DPI", unset = "180")
)

OUTPUT_ROOT_ENV <- Sys.getenv(
  "PANEL17_OUTPUT_ROOT",
  unset = ""
)

if (!requireNamespace("mixOmics", quietly = TRUE)) {
  stop(
    paste0(
      "Package 'mixOmics' is required.\n",
      "Install it with:\n",
      "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
      "BiocManager::install('mixOmics')"
    ),
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
  !is.finite(CV_REPEATS) ||
  CV_REPEATS < 5L
) {
  stop(
    "CV_REPEATS must be >=5.",
    call. = FALSE
  )
}

if (
  !is.finite(MAX_PLS_COMPONENTS) ||
  MAX_PLS_COMPONENTS < 1L
) {
  stop(
    "MAX_PLS_COMPONENTS must be >=1.",
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

  miss <- setdiff(
    required,
    names(x)
  )

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
# 2. LOCATE SELECTION / BIOLOGICAL / SCRIPT-10 RUNS
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
    "No complete selection run found."
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
    "No complete BIOLOGICAL_PROGRAMS run found."
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
    "No complete script-10 run found."
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

SELECTION <- if (
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

BIO_RUN <- if (
  nzchar(BIO_PROGRAM_RUN)
) {
  canonical(
    BIO_PROGRAM_RUN
  )
} else {
  find_latest_bio_run(
    SELECTION
  )
}

PANEL10 <- find_panel10_run(
  BIO_RUN
)

# ==============================================================================
# 3. LOCATE STRICT FRAMEWORK EVIDENCE-A POOL
# ==============================================================================

find_framework_run <- function(
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
      "Missing BIDIRECTIONAL_FRAMEWORK root:",
      root
    )
  )

  required <- file.path(
    "tables",
    "06C_HIGH_CONFIDENCE_AA_SAME_PATTERN_TRAIN_ONLY.csv"
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
          required
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
        required
      )
    )
  ]

  assert(
    length(dirs) > 0L,
    paste(
      "No complete bidirectional-framework run with strict AA pool for",
      family
    )
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

FRAMEWORK_RUNS <- setNames(
  lapply(
    FAMILIES,
    function(family) {
      find_framework_run(
        SELECTION,
        family
      )
    }
  ),
  FAMILIES
)

# ==============================================================================
# 4. OUTPUT
# ==============================================================================

if (nzchar(OUTPUT_ROOT_ENV)) {
  OUT_ROOT <- path.expand(
    OUTPUT_ROOT_ENV
  )
} else {
  OUT_ROOT <- file.path(
    BIO_RUN,
    "FINAL_EVIDENCE_A_PLSDA_PCA",
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
message("METATEST — STRICT EVIDENCE-A -> PLS-DA + PCA AUDIT")
message("Selection run: ", SELECTION)
message("BIOLOGICAL_PROGRAMS run: ", BIO_RUN)
message("Script-10 run: ", PANEL10)
message("Output: ", OUT_ROOT)
message("============================================================")

# ==============================================================================
# 5. READ SCRIPT-10 BIOLOGICAL CANDIDATES + INPUT AUDIT
# ==============================================================================

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
    "feature_label"
  )
)

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

shortlist$family <- as.character(
  shortlist$family
)

shortlist$view_code <- tolower(
  as.character(
    shortlist$view_code
  )
)

shortlist$feature_original <- as.character(
  shortlist$feature_original
)

shortlist <- shortlist[
  shortlist$family %in% FAMILIES &
    shortlist$view_code %in% VIEW_LEVELS,
  ,
  drop = FALSE
]

# ==============================================================================
# 6. INTERSECTION: SCRIPT-10 CANDIDATES ∩ STRICT AA FRAMEWORK
# ==============================================================================

aa_all <- list()
intersection_all <- list()
intersection_audit <- list()

for (family in FAMILIES) {
  framework_run <- FRAMEWORK_RUNS[[family]]

  aa <- read_csv_checked(
    file.path(
      framework_run,
      "tables",
      "06C_HIGH_CONFIDENCE_AA_SAME_PATTERN_TRAIN_ONLY.csv"
    ),
    required = c(
      "feature_model",
      "feature_original",
      "view_code",
      "rec_evidence_grade",
      "obs_evidence_grade",
      "rec_bidirectional_pattern",
      "obs_bidirectional_pattern",
      "direction_concordance_n",
      "high_confidence_AA_same_pattern"
    )
  )

  aa$family <- family
  aa$view_code <- tolower(
    as.character(
      aa$view_code
    )
  )
  aa$feature_original <- as.character(
    aa$feature_original
  )

  # Defensive confirmation: the exported strict pool must really be A/A.
  assert(
    all(
      aa$rec_evidence_grade ==
        "A_strong"
    ),
    paste(
      "Strict AA file contains non-A reconstructed evidence:",
      family
    )
  )

  assert(
    all(
      aa$obs_evidence_grade ==
        "A_strong"
    ),
    paste(
      "Strict AA file contains non-A observed evidence:",
      family
    )
  )

  aa$key_intersection <- paste(
    aa$family,
    aa$view_code,
    aa$feature_original,
    sep = "\r"
  )

  s <- shortlist[
    shortlist$family == family,
    ,
    drop = FALSE
  ]

  s$key_intersection <- paste(
    s$family,
    s$view_code,
    s$feature_original,
    sep = "\r"
  )

  in_aa <- s$key_intersection %in%
    aa$key_intersection

  audit <- data.frame(
    family = family,
    predictor_key =
      s$predictor_key,
    view_code =
      s$view_code,
    feature_original =
      s$feature_original,
    feature_label =
      s$feature_label,
    in_script10_biological_shortlist =
      TRUE,
    in_strict_AA_framework =
      in_aa,
    selected_for_plsda =
      in_aa,
    stringsAsFactors = FALSE
  )

  selected <- s[
    in_aa,
    ,
    drop = FALSE
  ]

  if (nrow(selected)) {
    aa_sub <- aa[
      match(
        selected$key_intersection,
        aa$key_intersection
      ),
      ,
      drop = FALSE
    ]

    extra_cols <- setdiff(
      names(aa_sub),
      c(
        "family",
        "view_code",
        "feature_original",
        "key_intersection"
      )
    )

    aa_extra <- aa_sub[
      ,
      extra_cols,
      drop = FALSE
    ]

    names(aa_extra) <- paste0(
      "framework_",
      names(aa_extra)
    )

    selected <- cbind(
      selected,
      aa_extra
    )
  }

  aa_all[[family]] <- aa
  intersection_all[[family]] <- selected
  intersection_audit[[family]] <- audit
}

aa_master <- rbind_fill(
  aa_all
)

selected_master <- rbind_fill(
  intersection_all
)

audit_master <- rbind_fill(
  intersection_audit
)

write_csv(
  aa_master,
  file.path(
    DIR_TABLES,
    "00_STRICT_AA_FRAMEWORK_POOL_ALL_FAMILIES.csv"
  )
)

write_csv(
  audit_master,
  file.path(
    DIR_TABLES,
    "01_SCRIPT10_TO_EVIDENCE_A_INTERSECTION_AUDIT.csv"
  )
)

write_csv(
  selected_master,
  file.path(
    DIR_TABLES,
    "02_EVIDENCE_A_FEATURES_FOR_PLSDA_ALL_FAMILIES.csv"
  )
)

# ==============================================================================
# 7. LOAD OBSERVED TRAIN / TEST MATRICES
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
      "Expected one input-audit row for",
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

  list(
    rfm = rfm,
    train_meta = train_meta,
    test_meta = test_meta,
    path = canonical(path)
  )
}

build_evidence_matrix <- function(
    family,
    inp,
    selected
) {
  selected <- selected[
    selected$family == family,
    ,
    drop = FALSE
  ]

  if (!nrow(selected)) {
    return(
      list(
        x_train = matrix(
          numeric(),
          nrow = nrow(
            inp$train_meta
          ),
          ncol = 0
        ),
        x_test = matrix(
          numeric(),
          nrow = nrow(
            inp$test_meta
          ),
          ncol = 0
        ),
        metadata = selected,
        technical_audit = data.frame()
      )
    )
  }

  tr_cols <- list()
  te_cols <- list()
  audit <- list()
  rr <- 1L

  for (i in seq_len(nrow(selected))) {
    view <- as.character(
      selected$view_code[i]
    )

    feature <- as.character(
      selected$feature_original[i]
    )

    key <- as.character(
      selected$predictor_key[i]
    )

    xtr <- as.matrix(
      inp$rfm[[view]]$train
    )

    xte <- as.matrix(
      inp$rfm[[view]]$test
    )

    if (
      !feature %in%
        rownames(xtr)
    ) {
      audit[[rr]] <- data.frame(
        family = family,
        predictor_key = key,
        eligible = FALSE,
        reason = "ABSENT_TRAIN",
        stringsAsFactors = FALSE
      )

      rr <- rr + 1L
      next
    }

    a <- as.numeric(
      xtr[
        feature,
        inp$train_meta$sample_id
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
        eligible = FALSE,
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
        eligible = FALSE,
        reason = "ZERO_VARIANCE_TRAIN",
        stringsAsFactors = FALSE
      )

      rr <- rr + 1L
      next
    }

    b <- if (
      feature %in%
        rownames(xte)
    ) {
      as.numeric(
        xte[
          feature,
          inp$test_meta$sample_id
        ]
      )
    } else {
      rep(
        NA_real_,
        nrow(
          inp$test_meta
        )
      )
    }

    n_imp <- sum(
      !is.finite(b)
    )

    if (n_imp > 0L) {
      b[
        !is.finite(b)
      ] <- stats::median(a)
    }

    tr_cols[[key]] <- a
    te_cols[[key]] <- b

    audit[[rr]] <- data.frame(
      family = family,
      predictor_key = key,
      eligible = TRUE,
      reason = "OK",
      test_values_imputed_from_TRAIN_median =
        n_imp,
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  if (!length(tr_cols)) {
    return(
      list(
        x_train = matrix(
          numeric(),
          nrow = nrow(
            inp$train_meta
          ),
          ncol = 0
        ),
        x_test = matrix(
          numeric(),
          nrow = nrow(
            inp$test_meta
          ),
          ncol = 0
        ),
        metadata = selected[
          0,
          ,
          drop = FALSE
        ],
        technical_audit =
          rbind_fill(
            audit
          )
      )
    )
  }

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
  ) <- inp$train_meta$sample_id

  rownames(
    x_test
  ) <- inp$test_meta$sample_id

  meta <- selected[
    match(
      colnames(x_train),
      selected$predictor_key
    ),
    ,
    drop = FALSE
  ]

  list(
    x_train = x_train,
    x_test = x_test,
    metadata = meta,
    technical_audit =
      rbind_fill(
        audit
      )
  )
}

family_input <- setNames(
  lapply(
    FAMILIES,
    load_family_input
  ),
  FAMILIES
)

family_matrix <- setNames(
  lapply(
    FAMILIES,
    function(family) {
      build_evidence_matrix(
        family = family,
        inp = family_input[[family]],
        selected = selected_master
      )
    }
  ),
  FAMILIES
)

technical_audit <- rbind_fill(
  lapply(
    family_matrix,
    function(z) {
      z$technical_audit
    }
  )
)

write_csv(
  technical_audit,
  file.path(
    DIR_TABLES,
    "03_EVIDENCE_A_TECHNICAL_ELIGIBILITY.csv"
  )
)

# ==============================================================================
# 8. STRATIFIED FOLDS + METRICS
# ==============================================================================

make_foldid <- function(
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
        "has fewer samples than folds."
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

  foldid
}

extract_pls_class <- function(
    pred_obj,
    ncomp,
    n_expected,
    distance = "max.dist"
) {
  cl <- pred_obj$class[[distance]]

  if (is.null(dim(cl))) {
    out <- as.character(cl)

    assert(
      length(out) ==
        n_expected,
      "PLS-DA predicted-class length mismatch."
    )

    return(out)
  }

  cl <- as.matrix(cl)

  if (nrow(cl) == n_expected) {
    idx <- min(
      ncomp,
      ncol(cl)
    )

    return(
      as.character(
        cl[, idx]
      )
    )
  }

  if (ncol(cl) == n_expected) {
    idx <- min(
      ncomp,
      nrow(cl)
    )

    return(
      as.character(
        cl[idx, ]
      )
    )
  }

  stop(
    "Could not resolve PLS-DA predicted-class matrix orientation.",
    call. = FALSE
  )
}

classification_metrics <- function(
    truth,
    pred,
    family,
    split,
    ncomp = NA_integer_,
    repeat_id = NA_integer_
) {
  truth <- factor(
    truth,
    levels = CLASS_LEVELS
  )

  pred <- factor(
    pred,
    levels = CLASS_LEVELS
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
          ncomp = ncomp,
          repeat_id = repeat_id,
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

  balanced_accuracy <- mean(
    per_class$recall
  )

  macro_f1 <- mean(
    per_class$f1
  )

  list(
    summary = data.frame(
      family = family,
      split = split,
      ncomp = ncomp,
      repeat_id = repeat_id,
      n = length(truth),
      accuracy = accuracy,
      balanced_accuracy =
        balanced_accuracy,
      BER =
        1 -
        balanced_accuracy,
      macro_f1 =
        macro_f1,
      stringsAsFactors = FALSE
    ),
    per_class = per_class,
    confusion = cm
  )
}

# ==============================================================================
# 9. REPEATED STRATIFIED 3-FOLD PLS-DA TUNING
# ==============================================================================

tune_plsda_components <- function(
    family,
    x,
    y,
    seed
) {
  max_comp <- min(
    MAX_PLS_COMPONENTS,
    ncol(x),
    nrow(x) - 1L
  )

  assert(
    max_comp >= 1L,
    paste(
      "No valid PLS components:",
      family
    )
  )

  repeat_metrics <- list()
  predictions <- list()
  rr <- 1L
  pp <- 1L

  for (r in seq_len(CV_REPEATS)) {
    foldid <- make_foldid(
      y,
      seed +
        1000L *
          r
    )

    pred_by_comp <- lapply(
      seq_len(max_comp),
      function(k) {
        rep(
          NA_character_,
          length(y)
        )
      }
    )

    for (fold in seq_len(N_FOLDS)) {
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
        "A PLS-DA fold-training subset lost a class."
      )

      fit <- mixOmics::plsda(
        X = x[
          idx_train,
          ,
          drop = FALSE
        ],
        Y = y_sub,
        ncomp = max_comp,
        scale = TRUE
      )

      pr <- predict(
        fit,
        x[
          idx_valid,
          ,
          drop = FALSE
        ],
        dist = "max.dist"
      )

      for (k in seq_len(max_comp)) {
        pred_by_comp[[k]][
          idx_valid
        ] <- extract_pls_class(
          pr,
          ncomp = k,
          n_expected =
            length(
              idx_valid
            ),
          distance =
            "max.dist"
        )
      }
    }

    for (k in seq_len(max_comp)) {
      assert(
        !anyNA(
          pred_by_comp[[k]]
        ),
        paste(
          "Missing CV predictions:",
          family,
          "repeat",
          r,
          "component",
          k
        )
      )

      ev <- classification_metrics(
        truth = y,
        pred =
          pred_by_comp[[k]],
        family = family,
        split =
          "TRAIN_REPEATED_3FOLD_CONDITIONAL",
        ncomp = k,
        repeat_id = r
      )

      repeat_metrics[[rr]] <-
        ev$summary

      rr <- rr + 1L

      predictions[[pp]] <- data.frame(
        family = family,
        repeat_id = r,
        ncomp = k,
        sample_id =
          rownames(x),
        truth =
          as.character(y),
        predicted =
          pred_by_comp[[k]],
        stringsAsFactors = FALSE
      )

      pp <- pp + 1L
    }
  }

  metrics <- rbind_fill(
    repeat_metrics
  )

  pred_long <- rbind_fill(
    predictions
  )

  summary_split <- split(
    metrics,
    metrics$ncomp
  )

  summary_rows <- lapply(
    summary_split,
    function(z) {
      n_ok <- sum(
        is.finite(
          z$BER
        )
      )

      data.frame(
        family = family,
        ncomp =
          z$ncomp[1L],
        n_repeats =
          n_ok,
        mean_accuracy =
          mean(
            z$accuracy,
            na.rm = TRUE
          ),
        mean_balanced_accuracy =
          mean(
            z$balanced_accuracy,
            na.rm = TRUE
          ),
        mean_BER =
          mean(
            z$BER,
            na.rm = TRUE
          ),
        sd_BER =
          stats::sd(
            z$BER,
            na.rm = TRUE
          ),
        se_BER =
          stats::sd(
            z$BER,
            na.rm = TRUE
          ) /
          sqrt(
            n_ok
          ),
        mean_macro_f1 =
          mean(
            z$macro_f1,
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

  min_idx <- which.min(
    summary$mean_BER
  )

  min_ber <- summary$mean_BER[
    min_idx
  ]

  min_se <- summary$se_BER[
    min_idx
  ]

  if (!is.finite(min_se)) {
    min_se <- 0
  }

  threshold <- min_ber +
    min_se

  candidates <- summary[
    summary$mean_BER <=
      threshold,
    ,
    drop = FALSE
  ]

  candidates <- candidates[
    order(
      candidates$ncomp,
      candidates$mean_BER
    ),
    ,
    drop = FALSE
  ]

  chosen_ncomp <- as.integer(
    candidates$ncomp[1L]
  )

  summary$within_1se_like <-
    summary$mean_BER <=
      threshold

  summary$selected_ncomp <-
    summary$ncomp ==
      chosen_ncomp

  list(
    repeat_metrics = metrics,
    prediction_long = pred_long,
    summary = summary,
    chosen_ncomp =
      chosen_ncomp,
    min_BER = min_ber,
    one_se_threshold =
      threshold
  )
}

# ==============================================================================
# 10. FINAL PLS-DA FIT + VIP / LOADINGS
# ==============================================================================

extract_vip_table <- function(
    fit,
    metadata,
    chosen_ncomp,
    family
) {
  v <- mixOmics::vip(
    fit
  )

  v <- as.matrix(v)

  if (is.null(rownames(v))) {
    rownames(v) <- colnames(
      fit$X
    )
  }

  n_use <- min(
    chosen_ncomp,
    ncol(v)
  )

  v_use <- v[
    ,
    seq_len(n_use),
    drop = FALSE
  ]

  vip_max <- apply(
    v_use,
    1L,
    max,
    na.rm = TRUE
  )

  vip_mean <- rowMeans(
    v_use,
    na.rm = TRUE
  )

  out <- data.frame(
    family = family,
    predictor_key =
      rownames(v_use),
    VIP_max_selected_components =
      vip_max,
    VIP_mean_selected_components =
      vip_mean,
    VIP_ge_reference_threshold =
      vip_max >=
      VIP_REFERENCE_THRESHOLD,
    VIP_reference_threshold =
      VIP_REFERENCE_THRESHOLD,
    stringsAsFactors = FALSE
  )

  for (k in seq_len(n_use)) {
    out[[paste0(
      "VIP_comp",
      k
    )]] <- v_use[, k]
  }

  out <- merge(
    out,
    metadata,
    by = c(
      "family",
      "predictor_key"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  out <- out[
    order(
      -out$VIP_max_selected_components,
      match(
        out$view_code,
        VIEW_LEVELS
      ),
      out$feature_label
    ),
    ,
    drop = FALSE
  ]

  out$VIP_rank <-
    seq_len(
      nrow(out)
    )

  out
}

extract_loading_table <- function(
    fit,
    metadata,
    chosen_ncomp,
    family
) {
  L <- as.matrix(
    fit$loadings$X
  )

  n_use <- min(
    chosen_ncomp,
    ncol(L)
  )

  rows <- list()
  rr <- 1L

  for (k in seq_len(n_use)) {
    rows[[rr]] <- data.frame(
      family = family,
      predictor_key =
        rownames(L),
      component = k,
      loading =
        L[, k],
      abs_loading =
        abs(
          L[, k]
        ),
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  out <- rbind_fill(
    rows
  )

  out <- merge(
    out,
    metadata,
    by = c(
      "family",
      "predictor_key"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  out <- out[
    order(
      out$component,
      -out$abs_loading
    ),
    ,
    drop = FALSE
  ]

  out
}

make_signed_loading_barplot <- function(
    dat,
    family,
    component,
    outdir,
    top_n = LOADING_BAR_TOPN,
    prefix = "FINAL"
) {
  sub <- dat[
    dat$family == family &
      dat$component == component,
    ,
    drop = FALSE
  ]

  if (!nrow(sub)) {
    return(invisible(NULL))
  }

  sub <- sub[
    order(-sub$abs_loading, sub$feature_label),
    ,
    drop = FALSE
  ]

  sub <- sub[
    seq_len(min(top_n, nrow(sub))),
    ,
    drop = FALSE
  ]

  sub$direction <- ifelse(
    sub$loading >= 0,
    "Positive",
    "Negative"
  )

  sub$display_label <- paste0(
    toupper(sub$view_code),
    " | ",
    sub$feature_label
  )

  sub <- sub[
    order(sub$loading),
    ,
    drop = FALSE
  ]

  sub$display_label <- factor(
    sub$display_label,
    levels = sub$display_label
  )

  p <- ggplot2::ggplot(
    sub,
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
      title = paste0(
        family,
        " — final transversal loadings: component ",
        component
      ),
      subtitle = paste0(
        "Top ",
        min(top_n, nrow(sub)),
        " features from the final CORE+SHARED_2 transversal set. ",
        "Bar sign = direction; absolute height = contribution magnitude."
      ),
      x = NULL,
      y = "PLS loading"
    ) +
    theme_panel(9)

  ggplot2::ggsave(
    file.path(
      outdir,
      paste0(
        family,
        "_",
        prefix,
        "_LOADINGS_COMP",
        component,
        "_TOP",
        min(top_n, nrow(sub)),
        ".png"
      )
    ),
    p,
    width = 9,
    height = max(6, 0.24 * nrow(sub) + 2),
    dpi = DPI
  )

  invisible(sub)
}

make_loading_scatter_map <- function(
    dat,
    family,
    comp_x,
    comp_y,
    outdir,
    top_n = LOADING_MAP_TOPN,
    prefix = "FINAL"
) {
  dx <- dat[
    dat$family == family &
      dat$component == comp_x,
    ,
    drop = FALSE
  ]

  dy <- dat[
    dat$family == family &
      dat$component == comp_y,
    ,
    drop = FALSE
  ]

  if (!nrow(dx) || !nrow(dy)) {
    return(invisible(NULL))
  }

  m <- merge(
    dx[, c("family", "predictor_key", "loading", "abs_loading", "view_code", "feature_label")],
    dy[, c("family", "predictor_key", "loading", "abs_loading")],
    by = c("family", "predictor_key"),
    suffixes = c("_x", "_y"),
    all = FALSE,
    sort = FALSE
  )

  if (!nrow(m)) {
    return(invisible(NULL))
  }

  m$importance <- pmax(
    m$abs_loading_x,
    m$abs_loading_y
  )

  m <- m[
    order(-m$importance, m$feature_label),
    ,
    drop = FALSE
  ]

  m_top <- m[
    seq_len(min(top_n, nrow(m))),
    ,
    drop = FALSE
  ]

  p <- ggplot2::ggplot(
    m_top,
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
    ggplot2::geom_point(size = 2.5, alpha = 0.90) +
    ggplot2::geom_text(
      ggplot2::aes(label = feature_label),
      size = 3,
      hjust = -0.05,
      vjust = 0.4,
      check_overlap = TRUE
    ) +
    ggplot2::labs(
      title = paste0(
        family,
        " — final transversal loading map: comp ",
        comp_x,
        " vs comp ",
        comp_y
      ),
      subtitle = paste0(
        "Top ",
        min(top_n, nrow(m_top)),
        " final CORE+SHARED_2 features by max absolute loading across both components."
      ),
      x = paste0("PLS loading comp ", comp_x),
      y = paste0("PLS loading comp ", comp_y),
      colour = "Omic"
    ) +
    theme_panel(9)

  ggplot2::ggsave(
    file.path(
      outdir,
      paste0(
        family,
        "_",
        prefix,
        "_LOADING_MAP_COMP",
        comp_x,
        "_",
        comp_y,
        ".png"
      )
    ),
    p,
    width = 8,
    height = 6,
    dpi = DPI
  )

  invisible(m_top)
}

# ==============================================================================
# 11. PCA AUDIT ON EXACTLY THE SAME FEATURES
# ==============================================================================


run_pca_audit <- function(
    family,
    x_train,
    x_test,
    train_meta,
    test_meta
) {
  pca <- stats::prcomp(
    x_train,
    center = TRUE,
    scale. = TRUE
  )

  tr_scores <- as.data.frame(
    pca$x
  )

  te_scores <- as.data.frame(
    stats::predict(
      pca,
      newdata = x_test
    )
  )

  tr_scores$sample_id <-
    rownames(x_train)

  te_scores$sample_id <-
    rownames(x_test)

  tr_scores$grupo <-
    as.character(
      train_meta$grupo[
        match(
          tr_scores$sample_id,
          train_meta$sample_id
        )
      ]
    )

  te_scores$grupo <-
    as.character(
      test_meta$grupo[
        match(
          te_scores$sample_id,
          test_meta$sample_id
        )
      ]
    )

  tr_scores$split <- "TRAIN"
  te_scores$split <- "TEST_EXPLORATORY_REUSED"

  variance <- pca$sdev^2

  variance <- variance /
    sum(variance)

  var_table <- data.frame(
    family = family,
    PC =
      paste0(
        "PC",
        seq_along(
          variance
        )
      ),
    variance_explained =
      variance,
    cumulative_variance =
      cumsum(
        variance
      ),
    stringsAsFactors = FALSE
  )

  list(
    model = pca,
    train_scores = tr_scores,
    test_scores = te_scores,
    variance = var_table
  )
}

# ==============================================================================
# 12. RUN FAMILIES
# ==============================================================================

feature_count_summary <- list()
cv_repeat_all <- list()
cv_summary_all <- list()
pls_metrics_all <- list()
pls_per_class_all <- list()
pls_predictions_all <- list()
vip_all <- list()
loading_all <- list()
pca_variance_all <- list()
pca_scores_all <- list()
model_status_all <- list()

for (family in FAMILIES) {
  message("")
  message("============================================================")
  message("[EVIDENCE A -> PLS-DA + PCA] ", family)
  message("============================================================")

  fm <- family_matrix[[family]]
  inp <- family_input[[family]]

  n_script10 <- sum(
    shortlist$family ==
      family
  )

  n_aa_framework <- sum(
    aa_master$family ==
      family
  )

  n_intersection <- ncol(
    fm$x_train
  )

  view_counts <- if (
    n_intersection > 0L
  ) {
    table(
      toupper(
        fm$metadata$view_code
      )
    )
  } else {
    integer()
  }

  feature_count_summary[[family]] <- data.frame(
    family = family,
    n_script10_biological_candidates =
      n_script10,
    n_strict_AA_framework_all =
      n_aa_framework,
    n_intersection_script10_x_AA =
      n_intersection,
    n_TX =
      if ("TX" %in% names(view_counts)) {
        as.integer(
          view_counts["TX"]
        )
      } else {
        0L
      },
    n_PR =
      if ("PR" %in% names(view_counts)) {
        as.integer(
          view_counts["PR"]
        )
      } else {
        0L
      },
    n_ME =
      if ("ME" %in% names(view_counts)) {
        as.integer(
          view_counts["ME"]
        )
      } else {
        0L
      },
    n_CL =
      if ("CL" %in% names(view_counts)) {
        as.integer(
          view_counts["CL"]
        )
      } else {
        0L
      },
    stringsAsFactors = FALSE
  )

  message(
    "script-10 candidates=",
    n_script10,
    " | strict AA framework=",
    n_aa_framework,
    " | intersection=",
    n_intersection
  )

  if (length(view_counts)) {
    message(
      "Views: ",
      paste(
        paste0(
          names(view_counts),
          "=",
          as.integer(
            view_counts
          )
        ),
        collapse = ", "
      )
    )
  }

  fam_tab <- file.path(
    DIR_TABLES,
    family
  )

  fam_fig <- file.path(
    DIR_FIGURES,
    family
  )

  dir.create(
    fam_tab,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    fam_fig,
    recursive = TRUE,
    showWarnings = FALSE
  )

  write_csv(
    fm$metadata,
    file.path(
      fam_tab,
      paste0(
        family,
        "_01_EVIDENCE_A_MODEL_FEATURES.csv"
      )
    )
  )

  if (n_intersection < 2L) {
    warning(
      family,
      ": fewer than 2 strict Evidence-A biological features. ",
      "PLS-DA/PCA model skipped."
    )

    model_status_all[[family]] <- data.frame(
      family = family,
      status = "SKIPPED",
      reason = "FEWER_THAN_2_STRICT_AA_FEATURES",
      n_features = n_intersection,
      stringsAsFactors = FALSE
    )

    next
  }

  x_train <- fm$x_train
  x_test <- fm$x_test

  y_train <- inp$train_meta$grupo
  y_test <- inp$test_meta$grupo

  # ---------------------------------------------------------------------------
  # PLS-DA tuning
  # ---------------------------------------------------------------------------

  cv <- tune_plsda_components(
    family = family,
    x = x_train,
    y = y_train,
    seed =
      SEED +
      1000000L *
        match(
          family,
          FAMILIES
        )
  )

  chosen_ncomp <- cv$chosen_ncomp

  message(
    "Chosen PLS-DA components=",
    chosen_ncomp
  )

  final_fit_ncomp <- min(
    MAX_PLS_COMPONENTS,
    3L,
    ncol(x_train),
    nrow(x_train) - 1L
  )

  final_pls <- mixOmics::plsda(
    X = x_train,
    Y = y_train,
    ncomp = final_fit_ncomp,
    scale = TRUE
  )

  # Aggregate OOF predictions at chosen component across repeats by majority
  # vote. Per-repeat metrics are the main TRAIN development estimate.
  pred_chosen <- cv$prediction_long[
    cv$prediction_long$ncomp ==
      chosen_ncomp,
    ,
    drop = FALSE
  ]

  vote_split <- split(
    pred_chosen,
    pred_chosen$sample_id
  )

  oof_vote <- vapply(
    vote_split,
    function(z) {
      tab <- table(
        factor(
          z$predicted,
          levels = CLASS_LEVELS
        )
      )

      CLASS_LEVELS[
        which.max(tab)
      ]
    },
    character(1)
  )

  oof_vote <- oof_vote[
    rownames(
      x_train
    )
  ]

  ev_train_vote <- classification_metrics(
    truth = y_train,
    pred = oof_vote,
    family = family,
    split =
      "TRAIN_REPEATED_3FOLD_MAJORITY_VOTE",
    ncomp = chosen_ncomp,
    repeat_id = NA_integer_
  )

  pr_test <- predict(
    final_pls,
    x_test,
    dist = "max.dist"
  )

  test_pred <- extract_pls_class(
    pr_test,
    ncomp =
      chosen_ncomp,
    n_expected =
      nrow(
        x_test
      ),
    distance =
      "max.dist"
  )

  ev_test <- classification_metrics(
    truth = y_test,
    pred = test_pred,
    family = family,
    split =
      "TEST_EXPLORATORY_REUSED",
    ncomp = chosen_ncomp,
    repeat_id = NA_integer_
  )

  # Mean repeated CV metrics for selected ncomp.
  selected_cv_metrics <- cv$repeat_metrics[
    cv$repeat_metrics$ncomp ==
      chosen_ncomp,
    ,
    drop = FALSE
  ]

  train_dev_summary <- data.frame(
    family = family,
    split =
      "TRAIN_REPEATED_3FOLD_CONDITIONAL_MEAN",
    ncomp =
      chosen_ncomp,
    n_repeats =
      nrow(
        selected_cv_metrics
      ),
    accuracy =
      mean(
        selected_cv_metrics$accuracy
      ),
    balanced_accuracy =
      mean(
        selected_cv_metrics$balanced_accuracy
      ),
    BER =
      mean(
        selected_cv_metrics$BER
      ),
    macro_f1 =
      mean(
        selected_cv_metrics$macro_f1
      ),
    stringsAsFactors = FALSE
  )

  test_summary <- ev_test$summary[
    ,
    c(
      "family",
      "split",
      "ncomp",
      "accuracy",
      "balanced_accuracy",
      "BER",
      "macro_f1"
    ),
    drop = FALSE
  ]

  test_summary$n_repeats <- NA_integer_

  test_summary <- test_summary[
    ,
    names(
      train_dev_summary
    ),
    drop = FALSE
  ]

  metric_out <- rbind(
    train_dev_summary,
    test_summary
  )

  # ---------------------------------------------------------------------------
  # VIP / loadings
  # ---------------------------------------------------------------------------

  vip_tab <- extract_vip_table(
    final_pls,
    metadata =
      fm$metadata,
    chosen_ncomp =
      chosen_ncomp,
    family = family
  )

  loading_tab <- extract_loading_table(
    final_pls,
    metadata =
      fm$metadata,
    chosen_ncomp =
      chosen_ncomp,
    family = family
  )

  # ---------------------------------------------------------------------------
  # PLS score coordinates for plotting/audit
  # ---------------------------------------------------------------------------

  n_score_components <- min(
    3L,
    ncol(
      final_pls$variates$X
    )
  )

  tr_pls_scores <- as.data.frame(
    final_pls$variates$X[
      ,
      seq_len(
        n_score_components
      ),
      drop = FALSE
    ]
  )

  colnames(
    tr_pls_scores
  ) <- paste0(
    "PLS",
    seq_len(
      ncol(
        tr_pls_scores
      )
    )
  )

  tr_pls_scores$sample_id <-
    rownames(
      x_train
    )

  tr_pls_scores$grupo <-
    as.character(
      y_train
    )

  tr_pls_scores$split <-
    "TRAIN"

  test_variates <- pr_test$variates

  if (is.list(test_variates) &&
      "X" %in% names(test_variates)) {
    te_score_mat <- as.matrix(
      test_variates$X
    )
  } else {
    te_score_mat <- as.matrix(
      test_variates
    )
  }

  if (nrow(te_score_mat) != nrow(x_test) &&
      ncol(te_score_mat) == nrow(x_test)) {
    te_score_mat <- t(
      te_score_mat
    )
  }

  n_score_use <- min(
    n_score_components,
    ncol(
      te_score_mat
    )
  )

  te_pls_scores <- as.data.frame(
    te_score_mat[
      ,
      seq_len(
        n_score_use
      ),
      drop = FALSE
    ]
  )

  colnames(
    te_pls_scores
  ) <- paste0(
    "PLS",
    seq_len(
      ncol(
        te_pls_scores
      )
    )
  )

  te_pls_scores$sample_id <-
    rownames(
      x_test
    )

  te_pls_scores$grupo <-
    as.character(
      y_test
    )

  te_pls_scores$split <-
    "TEST_EXPLORATORY_REUSED"

  pls_scores <- rbind_fill(
    tr_pls_scores,
    te_pls_scores
  )

  # ---------------------------------------------------------------------------
  # PCA audit on exactly the same features
  # ---------------------------------------------------------------------------

  pca <- run_pca_audit(
    family = family,
    x_train = x_train,
    x_test = x_test,
    train_meta =
      inp$train_meta,
    test_meta =
      inp$test_meta
  )

  pca_scores <- rbind_fill(
    pca$train_scores,
    pca$test_scores
  )

  # ---------------------------------------------------------------------------
  # Tables
  # ---------------------------------------------------------------------------

  write_csv(
    cv$repeat_metrics,
    file.path(
      fam_tab,
      paste0(
        family,
        "_02_PLSDA_REPEATED_CV_BY_REPEAT.csv"
      )
    )
  )

  write_csv(
    cv$summary,
    file.path(
      fam_tab,
      paste0(
        family,
        "_03_PLSDA_COMPONENT_TUNING.csv"
      )
    )
  )

  write_csv(
    metric_out,
    file.path(
      fam_tab,
      paste0(
        family,
        "_04_PLSDA_METRICS.csv"
      )
    )
  )

  write_csv(
    rbind(
      ev_train_vote$per_class,
      ev_test$per_class
    ),
    file.path(
      fam_tab,
      paste0(
        family,
        "_05_PLSDA_PER_CLASS_METRICS.csv"
      )
    )
  )

  write_csv(
    vip_tab,
    file.path(
      fam_tab,
      paste0(
        family,
        "_06_PLSDA_VIP_PRIORITY.csv"
      )
    )
  )

  write_csv(
    loading_tab,
    file.path(
      fam_tab,
      paste0(
        family,
        "_07_PLSDA_LOADINGS.csv"
      )
    )
  )

  write_csv(
    pls_scores,
    file.path(
      fam_tab,
      paste0(
        family,
        "_08_PLSDA_SCORES.csv"
      )
    )
  )

  write_csv(
    pca$variance,
    file.path(
      fam_tab,
      paste0(
        family,
        "_09_PCA_VARIANCE_EXPLAINED.csv"
      )
    )
  )

  write_csv(
    pca_scores,
    file.path(
      fam_tab,
      paste0(
        family,
        "_10_PCA_SCORES.csv"
      )
    )
  )

  write_csv(
    as.data.frame.matrix(
      ev_test$confusion
    ),
    file.path(
      fam_tab,
      paste0(
        family,
        "_11_TEST_EXPLORATORY_CONFUSION.csv"
      )
    )
  )

  saveRDS(
    list(
      family = family,
      evidence_A_features =
        fm$metadata,
      chosen_ncomp =
        chosen_ncomp,
      fitted_ncomp_for_scores =
        final_fit_ncomp,
      plsda =
        final_pls,
      pca =
        pca$model,
      source_framework_run =
        FRAMEWORK_RUNS[[family]],
      source_panel10_run =
        PANEL10
    ),
    file.path(
      DIR_MODELS,
      paste0(
        family,
        "_EVIDENCE_A_PLSDA_PCA.rds"
      )
    )
  )

  # ---------------------------------------------------------------------------
  # Figures
  # ---------------------------------------------------------------------------

  pls_pairs <- list(
    c(1L, 2L),
    c(2L, 3L),
    c(1L, 3L)
  )

  for (pair in pls_pairs) {
    xvar <- paste0(
      "PLS",
      pair[1L]
    )

    yvar <- paste0(
      "PLS",
      pair[2L]
    )

    if (
      all(
        c(
          xvar,
          yvar
        ) %in%
          names(
            pls_scores
          )
      )
    ) {
      plot_df <- pls_scores
      plot_df$score_x <- plot_df[[xvar]]
      plot_df$score_y <- plot_df[[yvar]]

      p_pls <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x = score_x,
          y = score_y,
          colour = grupo,
          shape = split
        )
      ) +
        ggplot2::geom_point(
          size = 3,
          alpha = 0.90
        ) +
        ggplot2::labs(
          title =
            paste0(
              family,
              " — strict Evidence-A PLS-DA: ",
              xvar,
              " vs ",
              yvar
            ),
          subtitle =
            paste0(
              "Predictive model uses ",
              chosen_ncomp,
              " component(s); ",
              "up to 3 components are fitted only to display latent score geometry. ",
              "TRAIN fit; TEST projected."
            ),
          x = paste0(
            "PLS component ",
            pair[1L]
          ),
          y = paste0(
            "PLS component ",
            pair[2L]
          ),
          colour = "Phenotype",
          shape = "Split"
        ) +
        theme_panel(10)

      ggplot2::ggsave(
        file.path(
          fam_fig,
          paste0(
            family,
            "_01_PLSDA_SCORES_COMP",
            pair[1L],
            "_",
            pair[2L],
            ".png"
          )
        ),
        p_pls,
        width = 7,
        height = 6,
        dpi = DPI
      )
    }
  }

  pca_pairs <- list(
    c(1L, 2L),
    c(2L, 3L),
    c(1L, 3L)
  )

  for (pair in pca_pairs) {
    xvar <- paste0(
      "PC",
      pair[1L]
    )

    yvar <- paste0(
      "PC",
      pair[2L]
    )

    if (
      all(
        c(
          xvar,
          yvar
        ) %in%
          names(
            pca_scores
          )
      )
    ) {
      plot_df <- pca_scores
      plot_df$score_x <- plot_df[[xvar]]
      plot_df$score_y <- plot_df[[yvar]]

      x_pct <- 100 *
        pca$variance$variance_explained[
          pair[1L]
        ]

      y_pct <- 100 *
        pca$variance$variance_explained[
          pair[2L]
        ]

      p_pca <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x = score_x,
          y = score_y,
          colour = grupo,
          shape = split
        )
      ) +
        ggplot2::geom_point(
          size = 3,
          alpha = 0.90
        ) +
        ggplot2::labs(
          title =
            paste0(
              family,
              " — PCA audit: ",
              xvar,
              " vs ",
              yvar
            ),
          subtitle =
            paste0(
              "Unsupervised PCA on exactly the same strict Evidence-A features; ",
              "TRAIN fit; TEST projected."
            ),
          x =
            paste0(
              xvar,
              " (",
              sprintf(
                "%.1f%%",
                x_pct
              ),
              ")"
            ),
          y =
            paste0(
              yvar,
              " (",
              sprintf(
                "%.1f%%",
                y_pct
              ),
              ")"
            ),
          colour = "Phenotype",
          shape = "Split"
        ) +
        theme_panel(10)

      ggplot2::ggsave(
        file.path(
          fam_fig,
          paste0(
            family,
            "_02_PCA_SCORES_COMP",
            pair[1L],
            "_",
            pair[2L],
            ".png"
          )
        ),
        p_pca,
        width = 7,
        height = 6,
        dpi = DPI
      )
    }
  }

  top_vip <- vip_tab[
    seq_len(
      min(
        30L,
        nrow(
          vip_tab
        )
      )
    ),
    ,
    drop = FALSE
  ]

  if (nrow(top_vip)) {
    top_vip$display_label <- paste0(
      toupper(
        top_vip$view_code
      ),
      " | ",
      top_vip$feature_label
    )

    top_vip$display_label <- factor(
      top_vip$display_label,
      levels =
        rev(
          unique(
            top_vip$display_label
          )
        )
    )

    p_vip <- ggplot2::ggplot(
      top_vip,
      ggplot2::aes(
        x = display_label,
        y = VIP_max_selected_components
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::geom_hline(
        yintercept =
          VIP_REFERENCE_THRESHOLD,
        linetype = 2
      ) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title =
          paste0(
            family,
            " — PLS-DA VIP prioritization"
          ),
        subtitle =
          "VIP>=1 is descriptive only; it does not redefine the Evidence-A model.",
        x = NULL,
        y = "Maximum VIP across selected PLS components"
      ) +
      theme_panel(9)

    ggplot2::ggsave(
      file.path(
        fam_fig,
        paste0(
          family,
          "_03_PLSDA_VIP_TOP30.png"
        )
      ),
      p_vip,
      width = 9,
      height = max(
        6,
        0.24 *
          nrow(
            top_vip
          ) +
          2
      ),
      dpi = DPI
    )
  }

  cv_repeat_all[[family]] <-
    cv$repeat_metrics

  cv_summary_all[[family]] <-
    cv$summary

  pls_metrics_all[[family]] <-
    metric_out

  pls_per_class_all[[family]] <-
    rbind(
      ev_train_vote$per_class,
      ev_test$per_class
    )

  pls_predictions_all[[family]] <-
    data.frame(
      family = family,
      sample_id =
        c(
          rownames(
            x_train
          ),
          rownames(
            x_test
          )
        ),
      split =
        c(
          rep(
            "TRAIN_REPEATED_3FOLD_MAJORITY_VOTE",
            nrow(
              x_train
            )
          ),
          rep(
            "TEST_EXPLORATORY_REUSED",
            nrow(
              x_test
            )
          )
        ),
      truth =
        c(
          as.character(
            y_train
          ),
          as.character(
            y_test
          )
        ),
      predicted =
        c(
          as.character(
            oof_vote
          ),
          as.character(
            test_pred
          )
        ),
      chosen_ncomp =
        chosen_ncomp,
      stringsAsFactors = FALSE
    )

  vip_all[[family]] <-
    vip_tab

  loading_all[[family]] <-
    loading_tab

  pca_variance_all[[family]] <-
    pca$variance

  pca_scores_all[[family]] <-
    pca_scores

  model_status_all[[family]] <- data.frame(
    family = family,
    status = "COMPLETE",
    reason = "",
    n_features =
      n_intersection,
    chosen_ncomp =
      chosen_ncomp,
    fitted_ncomp_for_scores =
      final_fit_ncomp,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# 13. MASTER OUTPUTS
# ==============================================================================

feature_count_master <- rbind_fill(
  feature_count_summary
)

cv_repeat_master <- rbind_fill(
  cv_repeat_all
)

cv_summary_master <- rbind_fill(
  cv_summary_all
)

pls_metrics_master <- rbind_fill(
  pls_metrics_all
)

pls_per_class_master <- rbind_fill(
  pls_per_class_all
)

pls_predictions_master <- rbind_fill(
  pls_predictions_all
)

vip_master <- rbind_fill(
  vip_all
)

loading_master <- rbind_fill(
  loading_all
)

pca_variance_master <- rbind_fill(
  pca_variance_all
)

pca_scores_master <- rbind_fill(
  pca_scores_all
)

model_status_master <- rbind_fill(
  model_status_all
)

write_csv(
  feature_count_master,
  file.path(
    DIR_TABLES,
    "04_EVIDENCE_A_FEATURE_COUNTS_BY_VIEW.csv"
  )
)

write_csv(
  cv_repeat_master,
  file.path(
    DIR_TABLES,
    "05_PLSDA_REPEATED_CV_BY_REPEAT_ALL_FAMILIES.csv"
  )
)

write_csv(
  cv_summary_master,
  file.path(
    DIR_TABLES,
    "06_PLSDA_COMPONENT_TUNING_ALL_FAMILIES.csv"
  )
)

write_csv(
  pls_metrics_master,
  file.path(
    DIR_TABLES,
    "07_PLSDA_METRICS_ALL_FAMILIES.csv"
  )
)

write_csv(
  pls_per_class_master,
  file.path(
    DIR_TABLES,
    "08_PLSDA_PER_CLASS_METRICS_ALL_FAMILIES.csv"
  )
)

write_csv(
  pls_predictions_master,
  file.path(
    DIR_TABLES,
    "09_PLSDA_PREDICTIONS_ALL_FAMILIES.csv"
  )
)

write_csv(
  vip_master,
  file.path(
    DIR_TABLES,
    "10_PLSDA_VIP_PRIORITY_ALL_FAMILIES.csv"
  )
)

write_csv(
  loading_master,
  file.path(
    DIR_TABLES,
    "11_PLSDA_LOADINGS_ALL_FAMILIES.csv"
  )
)

write_csv(
  pca_variance_master,
  file.path(
    DIR_TABLES,
    "12_PCA_VARIANCE_EXPLAINED_ALL_FAMILIES.csv"
  )
)

write_csv(
  pca_scores_master,
  file.path(
    DIR_TABLES,
    "13_PCA_SCORES_ALL_FAMILIES.csv"
  )
)

write_csv(
  model_status_master,
  file.path(
    DIR_TABLES,
    "14_MODEL_STATUS_ALL_FAMILIES.csv"
  )
)

# ==============================================================================
# 13B. CROSS-ANALYSIS STRICT-A INTERSECTIONS: P1 / P2 / P3
# ==============================================================================
#
# This section does NOT change the upstream Evidence-A definition and does NOT
# use TEST for selection. It only asks how reproducibly the same measurable
# feature survives the already-frozen strict-A biological filter across the
# three scientific analyses.
#
# Tiers:
#   CORE        = present in P1 + P2 + P3
#   SHARED_2    = present in exactly two analyses
#   SPECIFIC    = present in exactly one analysis
#
# The proposed transversal set for a later final classifier is CORE + SHARED_2.
# ANALYSIS-SPECIFIC features remain exported as biological context and are not
# deleted from the per-analysis PLS-DA models above.
# ==============================================================================

if (nrow(selected_master)) {
  consensus_split <- split(
    selected_master,
    paste(
      selected_master$view_code,
      selected_master$feature_original,
      sep = "\r"
    )
  )

  consensus <- do.call(
    rbind,
    lapply(
      consensus_split,
      function(z) {
        fam <- sort(
          unique(
            as.character(
              z$family
            )
          )
        )

        in_P1 <- "P1" %in% fam
        in_P2 <- "P2" %in% fam
        in_P3 <- "P3" %in% fam

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

        n_fam <- length(fam)

        consensus_tier <- if (
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
            n_fam,
          families =
            paste(
              fam,
              collapse = ";"
            ),
          in_P1 =
            in_P1,
          in_P2 =
            in_P2,
          in_P3 =
            in_P3,
          exact_intersection =
            exact_intersection,
          consensus_tier =
            consensus_tier,
          in_final_transversal_core_shared2 =
            n_fam >= 2L,
          framework_patterns =
            paste(
              sort(
                unique(
                  as.character(
                    z$framework_rec_bidirectional_pattern
                  )
                )
              ),
              collapse = ";"
            ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  rownames(
    consensus
  ) <- NULL

  consensus$exact_intersection <- factor(
    consensus$exact_intersection,
    levels = c(
      "P1&P2&P3",
      "P1&P2",
      "P1&P3",
      "P2&P3",
      "P1 only",
      "P2 only",
      "P3 only"
    )
  )

  consensus$consensus_tier <- factor(
    consensus$consensus_tier,
    levels = c(
      "CORE",
      "SHARED_2",
      "SPECIFIC"
    )
  )

  consensus <- consensus[
    order(
      consensus$consensus_tier,
      consensus$exact_intersection,
      match(
        consensus$view_code,
        VIEW_LEVELS
      ),
      consensus$feature_label
    ),
    ,
    drop = FALSE
  ]

  # Restore character columns before CSV export.
  consensus$exact_intersection <-
    as.character(
      consensus$exact_intersection
    )

  consensus$consensus_tier <-
    as.character(
      consensus$consensus_tier
    )

  write_csv(
    consensus,
    file.path(
      DIR_TABLES,
      "15_CROSS_ANALYSIS_STRICT_AA_CONSENSUS.csv"
    )
  )

  final_transversal <- consensus[
    consensus$in_final_transversal_core_shared2,
    ,
    drop = FALSE
  ]

  final_transversal <- final_transversal[
    order(
      -final_transversal$n_families,
      match(
        final_transversal$view_code,
        VIEW_LEVELS
      ),
      final_transversal$feature_label
    ),
    ,
    drop = FALSE
  ]

  write_csv(
    final_transversal,
    file.path(
      DIR_TABLES,
      "16_FINAL_TRANSVERSAL_CORE_SHARED2_SET.csv"
    )
  )

  # ---------------------------------------------------------------------------
  # Exact P1/P2/P3 intersection counts, overall and by omic
  # ---------------------------------------------------------------------------

  intersection_levels <- c(
    "P1&P2&P3",
    "P1&P2",
    "P1&P3",
    "P2&P3",
    "P1 only",
    "P2 only",
    "P3 only"
  )

  intersection_summary <- as.data.frame(
    table(
      factor(
        consensus$exact_intersection,
        levels = intersection_levels
      )
    ),
    stringsAsFactors = FALSE
  )

  names(intersection_summary) <- c(
    "exact_intersection",
    "n_features"
  )

  intersection_summary <- intersection_summary[
    intersection_summary$n_features > 0L,
    ,
    drop = FALSE
  ]

  write_csv(
    intersection_summary,
    file.path(
      DIR_TABLES,
      "17_STRICT_AA_INTERSECTION_COUNTS.csv"
    )
  )

  intersection_by_view <- as.data.frame(
    table(
      factor(
        consensus$exact_intersection,
        levels = intersection_levels
      ),
      factor(
        toupper(
          consensus$view_code
        ),
        levels = toupper(
          VIEW_LEVELS
        )
      )
    ),
    stringsAsFactors = FALSE
  )

  names(intersection_by_view) <- c(
    "exact_intersection",
    "view_code",
    "n_features"
  )

  intersection_by_view <- intersection_by_view[
    intersection_by_view$n_features > 0L,
    ,
    drop = FALSE
  ]

  write_csv(
    intersection_by_view,
    file.path(
      DIR_TABLES,
      "18_STRICT_AA_INTERSECTION_COUNTS_BY_OMIC.csv"
    )
  )

  # ---------------------------------------------------------------------------
  # Tier counts by omic: CORE / SHARED_2 / SPECIFIC
  # ---------------------------------------------------------------------------

  tier_by_view <- as.data.frame(
    table(
      factor(
        consensus$consensus_tier,
        levels = c(
          "CORE",
          "SHARED_2",
          "SPECIFIC"
        )
      ),
      factor(
        toupper(
          consensus$view_code
        ),
        levels = toupper(
          VIEW_LEVELS
        )
      )
    ),
    stringsAsFactors = FALSE
  )

  names(tier_by_view) <- c(
    "consensus_tier",
    "view_code",
    "n_features"
  )

  tier_by_view <- tier_by_view[
    tier_by_view$n_features > 0L,
    ,
    drop = FALSE
  ]

  write_csv(
    tier_by_view,
    file.path(
      DIR_TABLES,
      "19_STRICT_AA_CONSENSUS_TIER_COUNTS_BY_OMIC.csv"
    )
  )

  # ---------------------------------------------------------------------------
  # Figure 1: exact intersections, stacked by omic
  # ---------------------------------------------------------------------------

  plot_intersection <- intersection_by_view
  plot_intersection$exact_intersection <- factor(
    plot_intersection$exact_intersection,
    levels = intersection_levels
  )

  p_intersection <- ggplot2::ggplot(
    plot_intersection,
    ggplot2::aes(
      x = exact_intersection,
      y = n_features,
      fill = view_code
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title =
        "Strict Evidence-A feature intersections across P1, P2 and P3",
      subtitle =
        paste0(
          "Exact set membership; stacked bars show omic composition. ",
          "CORE=P1&P2&P3; SHARED_2=any exact two-analysis overlap."
        ),
      x = "Exact analysis intersection",
      y = "Number of strict Evidence-A features",
      fill = "Omic"
    ) +
    theme_panel(11) +
    ggplot2::theme(
      axis.text.x =
        ggplot2::element_text(
          angle = 30,
          hjust = 1
        )
    )

  ggplot2::ggsave(
    file.path(
      DIR_FIGURES,
      "FINAL_04_STRICT_AA_P1_P2_P3_INTERSECTIONS_BY_OMIC.png"
    ),
    p_intersection,
    width = 9,
    height = 6,
    dpi = DPI
  )

  # ---------------------------------------------------------------------------
  # Figure 2: CORE / SHARED_2 / SPECIFIC composition by omic
  # ---------------------------------------------------------------------------

  plot_tier <- tier_by_view
  plot_tier$consensus_tier <- factor(
    plot_tier$consensus_tier,
    levels = c(
      "CORE",
      "SHARED_2",
      "SPECIFIC"
    )
  )

  p_tier <- ggplot2::ggplot(
    plot_tier,
    ggplot2::aes(
      x = consensus_tier,
      y = n_features,
      fill = view_code
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title =
        "Cross-analysis strict Evidence-A consensus",
      subtitle =
        paste0(
          "CORE=3/3 analyses; SHARED_2=2/3; SPECIFIC=1/3. ",
          "The exported transversal set is CORE + SHARED_2."
        ),
      x = "Cross-analysis consensus tier",
      y = "Number of strict Evidence-A features",
      fill = "Omic"
    ) +
    theme_panel(11)

  ggplot2::ggsave(
    file.path(
      DIR_FIGURES,
      "FINAL_05_STRICT_AA_CONSENSUS_TIERS_BY_OMIC.png"
    ),
    p_tier,
    width = 8,
    height = 6,
    dpi = DPI
  )

  # ---------------------------------------------------------------------------
  # Final transversal loadings: only CORE + SHARED_2 features
  # ---------------------------------------------------------------------------

  final_loading_master <- merge(
    loading_master,
    final_transversal[, c(
      "view_code",
      "feature_original",
      "consensus_tier",
      "exact_intersection"
    )],
    by = c("view_code", "feature_original"),
    all = FALSE,
    sort = FALSE
  )

  if (nrow(final_loading_master)) {
    write_csv(
      final_loading_master,
      file.path(
        DIR_TABLES,
        "20_FINAL_TRANSVERSAL_PLSDA_LOADINGS.csv"
      )
    )

    dir.create(
      file.path(
        DIR_FIGURES,
        "FINAL_TRANSVERSAL_LOADINGS"
      ),
      recursive = TRUE,
      showWarnings = FALSE
    )

    for (family_i in FAMILIES) {
      fam_dat <- final_loading_master[
        final_loading_master$family == family_i,
        ,
        drop = FALSE
      ]

      if (!nrow(fam_dat)) {
        next
      }

      fam_dir <- file.path(
        DIR_FIGURES,
        "FINAL_TRANSVERSAL_LOADINGS",
        family_i
      )

      dir.create(
        fam_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )

      comps_here <- sort(unique(fam_dat$component))

      for (cc in comps_here) {
        make_signed_loading_barplot(
          dat = final_loading_master,
          family = family_i,
          component = cc,
          outdir = fam_dir,
          top_n = LOADING_BAR_TOPN,
          prefix = "FINAL_TRANSVERSAL"
        )
      }

      if (all(c(1, 2) %in% comps_here)) {
        make_loading_scatter_map(
          dat = final_loading_master,
          family = family_i,
          comp_x = 1L,
          comp_y = 2L,
          outdir = fam_dir,
          top_n = LOADING_MAP_TOPN,
          prefix = "FINAL_TRANSVERSAL"
        )
      }

      if (all(c(2, 3) %in% comps_here)) {
        make_loading_scatter_map(
          dat = final_loading_master,
          family = family_i,
          comp_x = 2L,
          comp_y = 3L,
          outdir = fam_dir,
          top_n = LOADING_MAP_TOPN,
          prefix = "FINAL_TRANSVERSAL"
        )
      }

      if (all(c(1, 3) %in% comps_here)) {
        make_loading_scatter_map(
          dat = final_loading_master,
          family = family_i,
          comp_x = 1L,
          comp_y = 3L,
          outdir = fam_dir,
          top_n = LOADING_MAP_TOPN,
          prefix = "FINAL_TRANSVERSAL"
        )
      }
    }
  }

  message("")
  message("------------------------------------------------------------")
  message("CROSS-ANALYSIS STRICT-A CONSENSUS")
  message(
    "CORE (3/3): ",
    sum(
      consensus$n_families == 3L
    )
  )
  message(
    "SHARED_2 (2/3): ",
    sum(
      consensus$n_families == 2L
    )
  )
  message(
    "SPECIFIC (1/3): ",
    sum(
      consensus$n_families == 1L
    )
  )
  message(
    "FINAL TRANSVERSAL CORE+SHARED_2: ",
    nrow(
      final_transversal
    )
  )
  message("------------------------------------------------------------")
}

# ==============================================================================
# 14. MANIFEST + README
# ==============================================================================

manifest <- data.frame(
  item = c(
    "SELECTION_RUN",
    "BIO_PROGRAM_RUN",
    "PANEL10_RUN",
    "OUTPUT_ROOT",
    "EVIDENCE_A_SOURCE",
    "EVIDENCE_A_RULE",
    "PLSDA_PACKAGE",
    "CV_FOLDS",
    "CV_REPEATS",
    "MAX_PLS_COMPONENTS",
    "PLSDA_COMPONENT_SELECTION",
    "PLSDA_SCORE_PLOTS",
    "VIP_ROLE",
    "PCA_ROLE",
    "TEST_STATUS",
    "VALIDATION_LIMITATION"
  ),
  value = c(
    SELECTION,
    BIO_RUN,
    PANEL10,
    OUT_ROOT,
    "06C_HIGH_CONFIDENCE_AA_SAME_PATTERN_TRAIN_ONLY.csv",
    paste0(
      "Intersection of script-10 biological shortlist with upstream strict AA pool: ",
      "A_strong reconstructed + A_strong observed + same non-uncertain pattern + ",
      "required directional concordance. This script does not redefine evidence_grade."
    ),
    as.character(
      utils::packageVersion(
        "mixOmics"
      )
    ),
    as.character(
      N_FOLDS
    ),
    as.character(
      CV_REPEATS
    ),
    as.character(
      MAX_PLS_COMPONENTS
    ),
    "Repeated stratified 3-fold TRAIN CV; choose smallest ncomp within one-SE-like BER region.",
    "Non-faceted scatter plots for component pairs 1-2, 2-3 and 1-3. Classification uses chosen_ncomp; component 3 is visualization/audit only when not selected.",
    paste0(
      "Descriptive prioritization only. VIP>=",
      VIP_REFERENCE_THRESHOLD,
      " does not alter the PLS-DA feature set."
    ),
    "Unsupervised audit on exactly the same strict Evidence-A features; TRAIN fit, TEST projection.",
    "TEST_EXPLORATORY_REUSED",
    paste0(
      "The strict Evidence-A pool and biological shortlist are defined from the full TRAIN ",
      "development set. Therefore repeated PLS-DA CV is conditional, not fully nested. ",
      "FD has only 3 TRAIN samples. TEST has already been inspected in prior method development."
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
  "METATEST — STRICT EVIDENCE-A -> PLS-DA + PCA AUDIT",
  "",
  "FINAL AGREED ALTERNATIVE",
  "",
  "1. Start from script-10 biological candidates.",
  "2. Keep only features that also belong to the strict A/A bidirectional-framework pool:",
  "     reconstructed A_strong",
  "     observed A_strong",
  "     same non-uncertain bidirectional pattern",
  "     required reconstructed/observed direction concordance.",
  "3. Run supervised PLS-DA on OBSERVED TRAIN values.",
  "4. Choose ncomp with repeated stratified 3-fold TRAIN CV.",
  "5. Fit final PLS-DA on all TRAIN and project reused TEST.",
  "6. Run PCA on EXACTLY the same features as an unsupervised audit only.",
  "7. Export non-faceted score scatter plots for component pairs 1-2, 2-3 and 1-3.",
  "   PLS-DA classification still uses the TRAIN-selected chosen_ncomp; component 3 is display/audit only when not selected.",
  "",
  "NO PLS-DA feature preselection is performed after Evidence A.",
  "VIP>=1 is only a descriptive prioritization annotation.",
  "",
  "MAIN TABLES",
  "  tables/01_SCRIPT10_TO_EVIDENCE_A_INTERSECTION_AUDIT.csv",
  "  tables/02_EVIDENCE_A_FEATURES_FOR_PLSDA_ALL_FAMILIES.csv",
  "  tables/04_EVIDENCE_A_FEATURE_COUNTS_BY_VIEW.csv",
  "  tables/06_PLSDA_COMPONENT_TUNING_ALL_FAMILIES.csv",
  "  tables/07_PLSDA_METRICS_ALL_FAMILIES.csv",
  "  tables/08_PLSDA_PER_CLASS_METRICS_ALL_FAMILIES.csv",
  "  tables/10_PLSDA_VIP_PRIORITY_ALL_FAMILIES.csv",
  "  tables/11_PLSDA_LOADINGS_ALL_FAMILIES.csv",
  "  tables/12_PCA_VARIANCE_EXPLAINED_ALL_FAMILIES.csv",
  "  tables/13_PCA_SCORES_ALL_FAMILIES.csv",
  "  tables/15_CROSS_ANALYSIS_STRICT_AA_CONSENSUS.csv",
  "  tables/16_FINAL_TRANSVERSAL_CORE_SHARED2_SET.csv",
  "  tables/17_STRICT_AA_INTERSECTION_COUNTS.csv",
  "  tables/18_STRICT_AA_INTERSECTION_COUNTS_BY_OMIC.csv",
  "  tables/19_STRICT_AA_CONSENSUS_TIER_COUNTS_BY_OMIC.csv",
  "  tables/20_FINAL_TRANSVERSAL_PLSDA_LOADINGS.csv",
  "",
  "CROSS-ANALYSIS CONSENSUS",
  "  CORE      = strict Evidence-A feature present in P1 + P2 + P3.",
  "  SHARED_2  = strict Evidence-A feature present in exactly two analyses.",
  "  SPECIFIC  = strict Evidence-A feature present in exactly one analysis.",
  "  FINAL TRANSVERSAL SET = CORE + SHARED_2; this is exported but is not",
  "  automatically refit as a new classifier in this script.",
  "",
  "INTERPRETATION",
  "If PLS-DA separates NP/FD/SO while PCA on the exact same features does not,",
  "that supports the hypothesis that supervised phenotype-aligned directions exist",
  "but are not dominant directions of total variance.",
  "",
  "If the strict A/A pool itself is still overwhelmingly metabolomic, PLS-DA cannot",
  "manufacture a genuinely multimodal signature; that would be evidence about the data,",
  "not a reason to impose an omic quota.",
  "",
  "TEST",
  "TEST is labeled TEST_EXPLORATORY_REUSED because it has already been inspected",
  "during earlier model-development iterations."
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
  "FINAL_EVIDENCE_A_PLSDA_PCA"
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
message("STRICT EVIDENCE-A -> PLS-DA + PCA COMPLETE")
message("Output: ", OUT_ROOT)
message("")

for (family in FAMILIES) {
  z <- feature_count_master[
    feature_count_master$family ==
      family,
    ,
    drop = FALSE
  ]

  message(
    family,
    ": script10=",
    z$n_script10_biological_candidates,
    " | strict AA=",
    z$n_strict_AA_framework_all,
    " | intersection=",
    z$n_intersection_script10_x_AA,
    " | TX=",
    z$n_TX,
    " PR=",
    z$n_PR,
    " ME=",
    z$n_ME,
    " CL=",
    z$n_CL
  )
}

message("============================================================")

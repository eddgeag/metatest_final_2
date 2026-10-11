# ==============================================================================
# METATEST — 22_final_panels_by_question_RSTUDIO.R
#
# CIERRE DEL PROYECTO — SOLO LECTURA DE RESULTADOS YA CALCULADOS
#
# Ejecutar en la CONSOLA de RStudio:
#
#   source("22_final_panels_by_question_RSTUDIO.R")
#
# NO guarda CSV, PNG, PDF ni RDS.
# NO reajusta PLS-DA.
# NO recalcula Evidence A.
# NO recalcula enriquecimientos.
# NO usa TEST para seleccionar nada.
#
# PRINCIPIO:
#   P1, P2 y P3 son preguntas científicas distintas.
#
# Por tanto:
#   PANEL P1 = strict Evidence-A de P1
#   PANEL P2 = strict Evidence-A de P2
#   PANEL P3 = strict Evidence-A de P3
#
# El solapamiento P1/P2/P3 se muestra DESPUÉS solo como auditoría.
#
# VIP, ||loading|| y concordancia PLS-DA son anotaciones del panel,
# no filtros de pertenencia.
#
# FIGURAS EN RSTUDIO:
#   1-3.  Direccionalidad del panel P1/P2/P3
#   4.    Clínicas por pregunta
#   5-7.  Resumen de vías por pregunta, usando TODAS las features del panel
#   8-10. Feature <-> pathway, versión legible por pregunta
#   11.   Auditoría de solapamiento P1/P2/P3
#
# Todos los puntos de los paneles aparecen en los plots de direccionalidad.
# Para los plots feature<->pathway se etiqueta un subconjunto de visualización
# para evitar un gráfico ilegible; la selección del panel NO cambia.
# ==============================================================================

options(stringsAsFactors = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Instala ggplot2.", call. = FALSE)
}

UPSTREAM_DIR <- Sys.getenv(
  "UPSTREAM_DIR",
  unset = path.expand(
    "~/Documents/PhD/metatest/results_upstream_local_20260908/comparacion_modelamientos/upstream_comun"
  )
)

PANEL17_RUN <- Sys.getenv("PANEL17_RUN", unset = "")

FAMILIES <- c("P1", "P2", "P3")
CLASS_LEVELS <- c("NP", "FD", "SO")
VIEW_LEVELS <- c("tx", "pr", "me", "cl")
VIP_REFERENCE <- 1

# Solo controla cuántas TX/ME se ETIQUETAN en plots densos.
# PR y CL se etiquetan siempre.
TOP_LABEL_PER_OMIC_DIRECTION <- as.integer(
  Sys.getenv(
    "PANEL22_TOP_LABEL_PER_OMIC_DIRECTION",
    unset = "5"
  )
)

# Número de vías mostradas por dirección en los resúmenes.
TOP_PATHWAYS_PER_DIRECTION <- as.integer(
  Sys.getenv(
    "PANEL22_TOP_PATHWAYS_PER_DIRECTION",
    unset = "8"
  )
)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

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

read_csv_checked <- function(path, required = character()) {
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

find_latest_selection <- function(upstream) {
  root <- file.path(
    upstream,
    "downstream",
    "seleccion_mejor_escenario_train"
  )

  assert(
    dir.exists(root),
    paste("Missing selection root:", root)
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

find_latest_bio <- function(selection) {
  root <- file.path(
    selection,
    "BIOLOGICAL_PROGRAMS"
  )

  assert(
    dir.exists(root),
    paste("Missing BIOLOGICAL_PROGRAMS:", root)
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
      dir.exists(p)
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

  assert(
    length(dirs) > 0L,
    "No BIOLOGICAL_PROGRAMS run found."
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

is_complete_panel17 <- function(p) {
  if (
    !nzchar(p) ||
    !dir.exists(p)
  ) {
    return(FALSE)
  }

  req <- c(
    file.path(
      "tables",
      "02_EVIDENCE_A_FEATURES_FOR_PLSDA_ALL_FAMILIES.csv"
    ),
    file.path(
      "tables",
      "10_PLSDA_VIP_PRIORITY_ALL_FAMILIES.csv"
    ),
    file.path(
      "tables",
      "11_PLSDA_LOADINGS_ALL_FAMILIES.csv"
    ),
    file.path(
      "tables",
      "14_MODEL_STATUS_ALL_FAMILIES.csv"
    )
  )

  all(
    file.exists(
      file.path(
        p,
        req
      )
    )
  )
}

find_panel17 <- function(bio) {
  if (nzchar(PANEL17_RUN)) {
    p <- canonical(
      PANEL17_RUN
    )

    assert(
      is_complete_panel17(p),
      "PANEL17_RUN is incomplete."
    )

    return(p)
  }

  root <- file.path(
    bio,
    "FINAL_EVIDENCE_A_PLSDA_PCA"
  )

  assert(
    dir.exists(root),
    paste(
      "Missing FINAL_EVIDENCE_A_PLSDA_PCA:",
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
      is_complete_panel17(p)
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
    vapply(
      dirs,
      is_complete_panel17,
      logical(1)
    )
  ]

  assert(
    length(dirs) > 0L,
    "No complete script-17 run found."
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

pattern_to_class <- function(pattern) {
  out <- rep(
    NA_character_,
    length(pattern)
  )

  out[pattern == "NP_like_high"] <- "NP"
  out[pattern == "FD_like_high"] <- "FD"
  out[pattern == "SO_like_high"] <- "SO"

  out
}

class_to_side <- function(cls) {
  ifelse(
    is.na(cls),
    NA_character_,
    paste0(
      cls,
      "_side"
    )
  )
}

theme_final <- function(base_size = 10) {
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

add_repel_or_text <- function(
    p,
    data,
    mapping,
    size = 3
) {
  if (!nrow(data)) {
    return(p)
  }

  if (
    requireNamespace(
      "ggrepel",
      quietly = TRUE
    )
  ) {
    p +
      ggrepel::geom_text_repel(
        data = data,
        mapping = mapping,
        inherit.aes = FALSE,
        size = size,
        max.overlaps = Inf,
        box.padding = 0.25,
        point.padding = 0.15,
        segment.alpha = 0.35
      )
  } else {
    p +
      ggplot2::geom_text(
        data = data,
        mapping = mapping,
        inherit.aes = FALSE,
        size = size,
        check_overlap = TRUE,
        hjust = -0.05,
        vjust = 0.4
      )
  }
}

# ------------------------------------------------------------------------------
# Locate frozen runs
# ------------------------------------------------------------------------------

SELECTION <- find_latest_selection(
  UPSTREAM_DIR
)

BIO <- find_latest_bio(
  SELECTION
)

PANEL17 <- find_panel17(
  BIO
)

manifest17 <- read_csv_checked(
  file.path(
    PANEL17,
    "RUN_MANIFEST.csv"
  ),
  required = c(
    "item",
    "value"
  )
)

PANEL10 <- as.character(
  manifest17$value[
    match(
      "PANEL10_RUN",
      manifest17$item
    )
  ]
)

assert(
  length(PANEL10) == 1L &&
    !is.na(PANEL10) &&
    dir.exists(PANEL10),
  "Could not resolve PANEL10_RUN."
)

manifest10 <- read_csv_checked(
  file.path(
    PANEL10,
    "RUN_MANIFEST.csv"
  ),
  required = c(
    "item",
    "value"
  )
)

SYNTHESIS_RUN <- as.character(
  manifest10$value[
    match(
      "BIO_SYNTHESIS_RUN",
      manifest10$item
    )
  ]
)

assert(
  length(SYNTHESIS_RUN) == 1L &&
    !is.na(SYNTHESIS_RUN) &&
    dir.exists(SYNTHESIS_RUN),
  "Could not resolve BIO_SYNTHESIS_RUN."
)

Q5_CATALOG <- file.path(
  SYNTHESIS_RUN,
  "plot_data",
  "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
  "00_INDEX",
  "FIG5_MASTER_PATHWAY_FEATURE_CATALOG.csv"
)

assert(
  file.exists(Q5_CATALOG),
  paste(
    "Frozen Q5 catalog was not found:",
    Q5_CATALOG
  )
)

cat("============================================================\n")
cat("METATEST — FINAL PANELS BY SCIENTIFIC QUESTION\n")
cat("P1/P2/P3 are treated independently.\n")
cat("Cross-analysis overlap is audit only.\n")
cat("Panel17:", PANEL17, "\n")
cat("Q5 catalog:", Q5_CATALOG, "\n")
cat("NO FILES WILL BE WRITTEN.\n")
cat("============================================================\n")

# ------------------------------------------------------------------------------
# Frozen tables
# ------------------------------------------------------------------------------

features <- read_csv_checked(
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
    "feature_label",
    "framework_rec_bidirectional_pattern",
    "framework_obs_bidirectional_pattern",
    "framework_rec_evidence_grade",
    "framework_obs_evidence_grade",
    "framework_direction_concordance_prop"
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

model_status <- read_csv_checked(
  file.path(
    PANEL17,
    "tables",
    "14_MODEL_STATUS_ALL_FAMILIES.csv"
  ),
  required = c(
    "family",
    "chosen_ncomp"
  )
)

catalog <- read_csv_checked(
  Q5_CATALOG,
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
    "is_tx_driver",
    "tx_driver_weight",
    "pathway_feature_cosine"
  )
)

# Basic types
features$family <- as.character(
  features$family
)

features$view_code <- tolower(
  as.character(
    features$view_code
  )
)

features$feature_original <- as.character(
  features$feature_original
)

features$feature_label <- as.character(
  features$feature_label
)

vip$family <- as.character(
  vip$family
)

vip$VIP_max_selected_components <- suppressWarnings(
  as.numeric(
    vip$VIP_max_selected_components
  )
)

loadings$family <- as.character(
  loadings$family
)

loadings$component <- suppressWarnings(
  as.integer(
    loadings$component
  )
)

loadings$loading <- suppressWarnings(
  as.numeric(
    loadings$loading
  )
)

loadings$abs_loading <- suppressWarnings(
  as.numeric(
    loadings$abs_loading
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

# ------------------------------------------------------------------------------
# Build one independent panel per question.
#
# PANEL MEMBERSHIP = strict Evidence-A feature already exported by script 17.
# No VIP/loadings threshold removes features here.
# ------------------------------------------------------------------------------

PANEL_BY_QUESTION <- list()

for (fam in FAMILIES) {
  f <- features[
    features$family == fam,
    ,
    drop = FALSE
  ]

  # Because 06C is same-pattern strict A, reconstructed and observed direction
  # should be the same. Preserve upstream truth and assert it.
  assert(
    all(
      f$framework_rec_evidence_grade == "A_strong"
    ),
    paste(
      fam,
      "contains non-A_strong reconstructed rows."
    )
  )

  assert(
    all(
      f$framework_obs_evidence_grade == "A_strong"
    ),
    paste(
      fam,
      "contains non-A_strong observed rows."
    )
  )

  assert(
    all(
      f$framework_rec_bidirectional_pattern ==
        f$framework_obs_bidirectional_pattern
    ),
    paste(
      fam,
      "contains reconstructed/observed pattern disagreement."
    )
  )

  f$framework_pattern <-
    f$framework_rec_bidirectional_pattern

  f$framework_direction <-
    pattern_to_class(
      f$framework_pattern
    )

  f$phenotypic_side <-
    class_to_side(
      f$framework_direction
    )

  # VIP is annotation only.
  v <- vip[
    vip$family == fam,
    c(
      "predictor_key",
      "VIP_max_selected_components"
    ),
    drop = FALSE
  ]

  f <- merge(
    f,
    v,
    by = "predictor_key",
    all.x = TRUE,
    sort = FALSE
  )

  f$VIP_ge_1 <-
    f$VIP_max_selected_components >=
      VIP_REFERENCE

  # PLS loading norm over the TRAIN-selected components.
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

  assert(
    is.finite(chosen) &&
      chosen >= 1L,
    paste(
      "Invalid chosen_ncomp for",
      fam
    )
  )

  lf <- loadings[
    loadings$family == fam &
      loadings$component <= chosen,
    ,
    drop = FALSE
  ]

  load_split <- split(
    lf,
    lf$predictor_key
  )

  load_summary <- do.call(
    rbind,
    lapply(
      load_split,
      function(z) {
        data.frame(
          predictor_key =
            z$predictor_key[1L],
          loading_norm =
            sqrt(
              sum(
                z$loading^2,
                na.rm = TRUE
              )
            ),
          max_abs_loading =
            max(
              z$abs_loading,
              na.rm = TRUE
            ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  f <- merge(
    f,
    load_summary,
    by = "predictor_key",
    all.x = TRUE,
    sort = FALSE
  )

  # ---------------------------------------------------------------------------
  # Existing PLS-DA directional corroboration.
  # No model is refit.
  # ---------------------------------------------------------------------------

  score_path <- file.path(
    PANEL17,
    "tables",
    fam,
    paste0(
      fam,
      "_08_PLSDA_SCORES.csv"
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

  score_cols <- paste0(
    "PLS",
    seq_len(chosen)
  )

  assert(
    all(
      score_cols %in%
        names(
          scores_f
        )
    ),
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

        mm <- colMeans(
          z,
          na.rm = TRUE
        )

        data.frame(
          grupo = cls,
          t(mm),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  pls_rows <- lapply(
    split(
      lf,
      lf$predictor_key
    ),
    function(z) {
      beta <- rep(
        0,
        chosen
      )

      for (cc in seq_len(chosen)) {
        hit <- z$loading[
          z$component == cc
        ]

        if (length(hit)) {
          beta[cc] <- hit[1L]
        }
      }

      align <- vapply(
        CLASS_LEVELS,
        function(cls) {
          mu <- as.numeric(
            class_means[
              class_means$grupo == cls,
              score_cols,
              drop = FALSE
            ]
          )

          sum(
            beta *
              mu
          )
        },
        numeric(1)
      )

      ord <- order(
        align,
        decreasing = TRUE
      )

      data.frame(
        predictor_key =
          z$predictor_key[1L],
        PLS_direction =
          names(
            align
          )[
            ord[1L]
          ],
        PLS_direction_margin =
          if (
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
          },
        stringsAsFactors = FALSE
      )
    }
  )

  pls_direction <- do.call(
    rbind,
    pls_rows
  )

  f <- merge(
    f,
    pls_direction,
    by = "predictor_key",
    all.x = TRUE,
    sort = FALSE
  )

  f$PLS_framework_concordant <-
    !is.na(
      f$framework_direction
    ) &
    f$PLS_direction ==
      f$framework_direction

  f$view_code <- factor(
    f$view_code,
    levels = VIEW_LEVELS
  )

  f$framework_direction <- factor(
    f$framework_direction,
    levels = CLASS_LEVELS
  )

  PANEL_BY_QUESTION[[fam]] <- f
}

# ------------------------------------------------------------------------------
# Console counts only — no long tables.
# ------------------------------------------------------------------------------

cat("\nQUESTION-SPECIFIC PANELS\n")

for (fam in FAMILIES) {
  f <- PANEL_BY_QUESTION[[fam]]

  vv <- table(
    factor(
      toupper(
        as.character(
          f$view_code
        )
      ),
      levels = toupper(
        VIEW_LEVELS
      )
    )
  )

  cat(
    fam,
    ": n=",
    nrow(f),
    " | TX=",
    vv["TX"],
    " PR=",
    vv["PR"],
    " ME=",
    vv["ME"],
    " CL=",
    vv["CL"],
    " | VIP>=1=",
    sum(
      f$VIP_ge_1,
      na.rm = TRUE
    ),
    " | PLS direction concordant=",
    sum(
      f$PLS_framework_concordant,
      na.rm = TRUE
    ),
    "/",
    nrow(f),
    "\n",
    sep = ""
  )
}

# ==============================================================================
# PLOT STORAGE
# ==============================================================================

PLOTS_FINAL <- list()

# ==============================================================================
# 1-3. PANEL DIRECTIONALITY BY QUESTION
# ==============================================================================

for (fam in FAMILIES) {
  f <- PANEL_BY_QUESTION[[fam]]

  # All PR + CL labels are ALWAYS displayed.
  always_label <- f[
    as.character(
      f$view_code
    ) %in%
      c(
        "pr",
        "cl"
      ),
    ,
    drop = FALSE
  ]

  # TX/ME: top contributors only for text labels.
  txme <- f[
    as.character(
      f$view_code
    ) %in%
      c(
        "tx",
        "me"
      ),
    ,
    drop = FALSE
  ]

  txme_split <- split(
    txme,
    paste(
      txme$framework_direction,
      txme$view_code,
      sep = "\r"
    )
  )

  txme_label <- do.call(
    rbind,
    lapply(
      txme_split,
      function(z) {
        z <- z[
          order(
            -z$loading_norm,
            -z$VIP_max_selected_components,
            z$feature_label
          ),
          ,
          drop = FALSE
        ]

        z[
          seq_len(
            min(
              TOP_LABEL_PER_OMIC_DIRECTION,
              nrow(z)
            )
          ),
          ,
          drop = FALSE
        ]
      }
    )
  )

  label_df <- rbind(
    always_label,
    txme_label
  )

  label_df <- label_df[
    !duplicated(
      label_df$predictor_key
    ),
    ,
    drop = FALSE
  ]

  f$VIP_status <- ifelse(
    f$VIP_ge_1,
    "VIP>=1",
    "VIP<1"
  )

  f$PLS_concordance <- ifelse(
    f$PLS_framework_concordant,
    "Framework = PLS",
    "Framework != PLS"
  )

  p <- ggplot2::ggplot(
    f,
    ggplot2::aes(
      x = framework_direction,
      y = loading_norm,
      colour = view_code,
      shape = PLS_concordance,
      alpha = VIP_status,
      size = VIP_max_selected_components
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(
      ~ view_code,
      scales = "free_y",
      ncol = 2
    ) +
    ggplot2::scale_alpha_manual(
      values = c(
        "VIP>=1" = 0.95,
        "VIP<1" = 0.28
      )
    ) +
    ggplot2::labs(
      title =
        paste0(
          fam,
          " — strict Evidence-A panel by phenotype direction"
        ),
      subtitle =
        paste0(
          "All panel features are shown. ",
          "PR/CL are always labelled; top TX/ME are labelled for readability. ",
          "Point size=VIP; y=PLS loading norm; shape=framework/PLS directional concordance."
        ),
      x =
        "Framework phenotype-high direction",
      y =
        "PLS loading norm across selected components",
      colour =
        "Omic",
      shape =
        "Directional corroboration",
      alpha =
        "VIP annotation",
      size =
        "VIP"
    ) +
    theme_final(9)

  p <- add_repel_or_text(
    p = p,
    data = label_df,
    mapping = ggplot2::aes(
      x = framework_direction,
      y = loading_norm,
      label = feature_label
    ),
    size = 2.8
  )

  PLOTS_FINAL[[paste0(
    fam,
    "_direction"
  )]] <- p

  print(p)
}

# ==============================================================================
# 4. CLINICAL FEATURES — ALL QUESTIONS
# ==============================================================================

clinical_all <- do.call(
  rbind,
  lapply(
    FAMILIES,
    function(fam) {
      z <- PANEL_BY_QUESTION[[fam]]

      z <- z[
        as.character(
          z$view_code
        ) == "cl",
        ,
        drop = FALSE
      ]

      if (!nrow(z)) {
        return(NULL)
      }

      z$family <- fam
      z
    }
  )
)

if (
  !is.null(clinical_all) &&
  nrow(clinical_all)
) {
  clinical_all$family <- factor(
    clinical_all$family,
    levels = FAMILIES
  )

  clinical_all$feature_label <- factor(
    clinical_all$feature_label,
    levels = rev(
      unique(
        clinical_all$feature_label[
          order(
            clinical_all$family,
            -clinical_all$loading_norm
          )
        ]
      )
    )
  )

  p_clin <- ggplot2::ggplot(
    clinical_all,
    ggplot2::aes(
      x = family,
      y = feature_label,
      colour = framework_direction,
      size = VIP_max_selected_components,
      shape = PLS_framework_concordant
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.9
    ) +
    ggplot2::labs(
      title =
        "Clinical Evidence-A features retained by each scientific question",
      subtitle =
        "No cross-question intersection filter is applied. This is why P1 recovers anthropometry/adiposity.",
      x =
        "Scientific question",
      y =
        "Clinical feature",
      colour =
        "Framework direction",
      size =
        "VIP",
      shape =
        "Framework = PLS direction"
    ) +
    theme_final(10)

  PLOTS_FINAL$clinical_by_question <-
    p_clin

  print(
    p_clin
  )
}

# ==============================================================================
# Build phenotype-aligned Q5 feature-pathway links for each question
# ==============================================================================

PATHWAY_LINKS <- list()
PATHWAY_SUMMARY <- list()
FEATURE_PATHWAY_DISPLAY <- list()

for (fam in FAMILIES) {
  f <- PANEL_BY_QUESTION[[fam]]

  f$key_feature <- paste(
    tolower(
      as.character(
        f$view_code
      )
    ),
    f$feature_original,
    sep = "\r"
  )

  cata <- catalog[
    catalog$family == fam,
    ,
    drop = FALSE
  ]

  cata$key_feature <- paste(
    cata$view_code,
    cata$feature_original,
    sep = "\r"
  )

  # Add framework direction of the feature.
  ann <- f[
    ,
    c(
      "key_feature",
      "framework_direction",
      "framework_pattern",
      "VIP_max_selected_components",
      "loading_norm",
      "PLS_framework_concordant"
    ),
    drop = FALSE
  ]

  ann$framework_direction <- as.character(
    ann$framework_direction
  )

  links <- merge(
    cata,
    ann,
    by = "key_feature",
    all = FALSE,
    sort = FALSE
  )

  links$expected_side <- class_to_side(
    links$framework_direction
  )

  # Keep only biological pathways on the same phenotype side as the feature.
  links <- links[
    !is.na(
      links$expected_side
    ) &
      links$phenotypic_side ==
        links$expected_side,
    ,
    drop = FALSE
  ]

  # Association-strength annotation:
  # TX membership comes from frozen driver weight;
  # non-TX comes from frozen pathway_feature_cosine.
  links$association_strength <- ifelse(
    links$view_code == "tx" &
      is.finite(
        links$tx_driver_weight
      ),
    abs(
      links$tx_driver_weight
    ),
    abs(
      links$pathway_feature_cosine
    )
  )

  links$association_strength[
    !is.finite(
      links$association_strength
    )
  ] <- NA_real_

  # Unique feature-pathway rows for counts.
  links$key_pathway <- paste(
    links$framework_direction,
    links$library,
    links$pathway,
    links$program_id,
    sep = "\r"
  )

  links$key_feature_pathway <- paste(
    links$key_feature,
    links$key_pathway,
    sep = "\r"
  )

  links_unique <- links[
    !duplicated(
      links$key_feature_pathway
    ),
    ,
    drop = FALSE
  ]

  PATHWAY_LINKS[[fam]] <-
    links_unique

  # ---------------------------------------------------------------------------
  # Pathway summary using ALL panel features.
  # ---------------------------------------------------------------------------

  psplit <- split(
    links_unique,
    links_unique$key_pathway
  )

  psum <- do.call(
    rbind,
    lapply(
      psplit,
      function(z) {
        data.frame(
          family = fam,
          framework_direction =
            z$framework_direction[1L],
          library =
            z$library[1L],
          pathway =
            z$pathway[1L],
          pathway_label =
            z$pathway_label[1L],
          program_short =
            z$program_short[1L],
          n_panel_features =
            length(
              unique(
                z$key_feature
              )
            ),
          n_TX =
            length(
              unique(
                z$key_feature[
                  z$view_code == "tx"
                ]
              )
            ),
          n_PR =
            length(
              unique(
                z$key_feature[
                  z$view_code == "pr"
                ]
              )
            ),
          n_ME =
            length(
              unique(
                z$key_feature[
                  z$view_code == "me"
                ]
              )
            ),
          n_CL =
            length(
              unique(
                z$key_feature[
                  z$view_code == "cl"
                ]
              )
            ),
          max_strength =
            if (
              any(
                is.finite(
                  z$association_strength
                )
              )
            ) {
              max(
                z$association_strength,
                na.rm = TRUE
              )
            } else {
              NA_real_
            },
          stringsAsFactors = FALSE
        )
      }
    )
  )

  psum <- psum[
    order(
      match(
        psum$framework_direction,
        CLASS_LEVELS
      ),
      -psum$n_panel_features,
      -psum$max_strength,
      psum$pathway_label
    ),
    ,
    drop = FALSE
  ]

  PATHWAY_SUMMARY[[fam]] <-
    psum

  # Top pathways PER DIRECTION for plots.
  top_path_rows <- do.call(
    rbind,
    lapply(
      CLASS_LEVELS,
      function(cls) {
        z <- psum[
          psum$framework_direction == cls,
          ,
          drop = FALSE
        ]

        if (!nrow(z)) {
          return(NULL)
        }

        z[
          seq_len(
            min(
              TOP_PATHWAYS_PER_DIRECTION,
              nrow(z)
            )
          ),
          ,
          drop = FALSE
        ]
      }
    )
  )

  top_keys <- paste(
    top_path_rows$framework_direction,
    top_path_rows$library,
    top_path_rows$pathway,
    sep = "\r"
  )

  links_unique$key_pathway_simple <- paste(
    links_unique$framework_direction,
    links_unique$library,
    links_unique$pathway,
    sep = "\r"
  )

  display_links <- links_unique[
    links_unique$key_pathway_simple %in%
      top_keys,
    ,
    drop = FALSE
  ]

  # Plot-display feature subset only:
  # all PR/CL + strongest TX/ME per phenotype direction.
  display_features <- list()
  dd <- 1L

  for (cls in CLASS_LEVELS) {
    zf <- f[
      as.character(
        f$framework_direction
      ) == cls,
      ,
      drop = FALSE
    ]

    keep_pc <- zf[
      as.character(
        zf$view_code
      ) %in%
        c(
          "pr",
          "cl"
        ),
      ,
      drop = FALSE
    ]

    keep_tm <- do.call(
      rbind,
      lapply(
        c(
          "tx",
          "me"
        ),
        function(vw) {
          zz <- zf[
            as.character(
              zf$view_code
            ) == vw,
            ,
            drop = FALSE
          ]

          if (!nrow(zz)) {
            return(NULL)
          }

          zz <- zz[
            order(
              -zz$loading_norm,
              -zz$VIP_max_selected_components,
              zz$feature_label
            ),
            ,
            drop = FALSE
          ]

          zz[
            seq_len(
              min(
                TOP_LABEL_PER_OMIC_DIRECTION,
                nrow(zz)
              )
            ),
            ,
            drop = FALSE
          ]
        }
      )
    )

    display_features[[dd]] <- rbind(
      keep_pc,
      keep_tm
    )

    dd <- dd + 1L
  }

  display_features <- do.call(
    rbind,
    display_features
  )

  display_keys <- paste(
    tolower(
      as.character(
        display_features$view_code
      )
    ),
    display_features$feature_original,
    sep = "\r"
  )

  display_links <- display_links[
    display_links$key_feature %in%
      display_keys,
    ,
    drop = FALSE
  ]

  FEATURE_PATHWAY_DISPLAY[[fam]] <-
    display_links
}

# ==============================================================================
# 5-7. PATHWAY SUMMARY PLOTS BY QUESTION — ALL PANEL FEATURES
# ==============================================================================

for (fam in FAMILIES) {
  psum <- PATHWAY_SUMMARY[[fam]]

  top <- do.call(
    rbind,
    lapply(
      CLASS_LEVELS,
      function(cls) {
        z <- psum[
          psum$framework_direction == cls,
          ,
          drop = FALSE
        ]

        if (!nrow(z)) {
          return(NULL)
        }

        z[
          seq_len(
            min(
              TOP_PATHWAYS_PER_DIRECTION,
              nrow(z)
            )
          ),
          ,
          drop = FALSE
        ]
      }
    )
  )

  long <- do.call(
    rbind,
    lapply(
      c(
        "TX",
        "PR",
        "ME",
        "CL"
      ),
      function(vw) {
        data.frame(
          family =
            top$family,
          framework_direction =
            top$framework_direction,
          library =
            top$library,
          pathway_label =
            top$pathway_label,
          Omic = vw,
          n_features =
            top[[
              paste0(
                "n_",
                vw
              )
            ]],
          stringsAsFactors = FALSE
        )
      }
    )
  )

  long <- long[
    long$n_features > 0L,
    ,
    drop = FALSE
  ]

  long$pathway_display <- paste0(
    long$library,
    " | ",
    long$pathway_label
  )

  # Factor order independently roughly follows ranking table.
  path_levels <- rev(
    unique(
      paste0(
        top$library,
        " | ",
        top$pathway_label
      )
    )
  )

  long$pathway_display <- factor(
    long$pathway_display,
    levels = path_levels
  )

  long$framework_direction <- factor(
    long$framework_direction,
    levels = CLASS_LEVELS
  )

  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = pathway_display,
      y = n_features,
      fill = Omic
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(
      ~ framework_direction,
      scales = "free_y",
      ncol = 1
    ) +
    ggplot2::labs(
      title =
        paste0(
          fam,
          " — phenotype-aligned pathway summary"
        ),
      subtitle =
        paste0(
          "Uses ALL strict Evidence-A panel features. ",
          "Bar segments are the number of TX/PR/ME/CL features linked to each frozen Q5 pathway."
        ),
      x = NULL,
      y =
        "Number of panel features linked",
      fill =
        "Omic"
    ) +
    theme_final(8)

  PLOTS_FINAL[[paste0(
    fam,
    "_pathway_summary"
  )]] <- p

  print(p)
}

# ==============================================================================
# 8-10. FEATURE <-> PATHWAY ASSOCIATION PLOTS BY QUESTION
# ==============================================================================

for (fam in FAMILIES) {
  z <- FEATURE_PATHWAY_DISPLAY[[fam]]

  if (
    is.null(z) ||
    !nrow(z)
  ) {
    next
  }

  z$feature_display <- paste0(
    toupper(
      z$view_code
    ),
    " | ",
    z$feature_label
  )

  z$pathway_display <- paste0(
    z$library,
    " | ",
    z$pathway_label
  )

  # Stable ordering by phenotype direction, omic and feature label.
  feature_order <- unique(
    z$feature_display[
      order(
        match(
          z$framework_direction,
          CLASS_LEVELS
        ),
        match(
          z$view_code,
          VIEW_LEVELS
        ),
        z$feature_label
      )
    ]
  )

  z$feature_display <- factor(
    z$feature_display,
    levels = rev(
      feature_order
    )
  )

  z$framework_direction <- factor(
    z$framework_direction,
    levels = CLASS_LEVELS
  )

  # Strength only controls bubble size; it is not a new selection statistic.
  z$association_strength_plot <-
    z$association_strength

  # If some links have no finite strength, still show them with a minimal point.
  min_pos <- suppressWarnings(
    min(
      z$association_strength_plot[
        is.finite(
          z$association_strength_plot
        ) &
          z$association_strength_plot > 0
      ],
      na.rm = TRUE
    )
  )

  if (!is.finite(min_pos)) {
    min_pos <- 1
  }

  z$association_strength_plot[
    !is.finite(
      z$association_strength_plot
    ) |
      z$association_strength_plot <= 0
  ] <- min_pos * 0.5

  p <- ggplot2::ggplot(
    z,
    ggplot2::aes(
      x = pathway_display,
      y = feature_display,
      size = association_strength_plot,
      colour = view_code
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.82
    ) +
    ggplot2::facet_wrap(
      ~ framework_direction,
      scales = "free",
      ncol = 1
    ) +
    ggplot2::labs(
      title =
        paste0(
          fam,
          " — feature <-> pathway associations"
        ),
      subtitle =
        paste0(
          "Frozen Q5 links on the SAME phenotype side as each feature. ",
          "All PR/CL are shown; strongest TX/ME are shown for readability. ",
          "This display filter does NOT change panel membership."
        ),
      x =
        "Phenotype-aligned pathway",
      y =
        "Panel feature",
      size =
        "|association|",
      colour =
        "Omic"
    ) +
    theme_final(7) +
    ggplot2::theme(
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1
        )
    )

  PLOTS_FINAL[[paste0(
    fam,
    "_feature_pathway"
  )]] <- p

  print(p)
}

# ==============================================================================
# 11. OVERLAP AUDIT — AFTER DEFINING P1/P2/P3 INDEPENDENTLY
# ==============================================================================

membership <- do.call(
  rbind,
  lapply(
    FAMILIES,
    function(fam) {
      f <- PANEL_BY_QUESTION[[fam]]

      data.frame(
        family = fam,
        feature_key = paste(
          tolower(
            as.character(
              f$view_code
            )
          ),
          f$feature_original,
          sep = "\r"
        ),
        view_code =
          tolower(
            as.character(
              f$view_code
            )
          ),
        feature_label =
          f$feature_label,
        stringsAsFactors = FALSE
      )
    }
  )
)

all_keys <- unique(
  membership$feature_key
)

overlap_rows <- lapply(
  all_keys,
  function(k) {
    z <- membership[
      membership$feature_key == k,
      ,
      drop = FALSE
    ]

    fam <- sort(
      unique(
        z$family
      )
    )

    exact <- if (
      identical(
        fam,
        c(
          "P1",
          "P2",
          "P3"
        )
      )
    ) {
      "P1&P2&P3"
    } else if (
      identical(
        fam,
        c(
          "P1",
          "P2"
        )
      )
    ) {
      "P1&P2"
    } else if (
      identical(
        fam,
        c(
          "P1",
          "P3"
        )
      )
    ) {
      "P1&P3"
    } else if (
      identical(
        fam,
        c(
          "P2",
          "P3"
        )
      )
    ) {
      "P2&P3"
    } else if (
      identical(
        fam,
        "P1"
      )
    ) {
      "P1 only"
    } else if (
      identical(
        fam,
        "P2"
      )
    ) {
      "P2 only"
    } else {
      "P3 only"
    }

    data.frame(
      exact_overlap = exact,
      view_code = z$view_code[1L],
      feature_label = z$feature_label[1L],
      stringsAsFactors = FALSE
    )
  }
)

overlap <- do.call(
  rbind,
  overlap_rows
)

overlap_levels <- c(
  "P1&P2&P3",
  "P1&P2",
  "P1&P3",
  "P2&P3",
  "P1 only",
  "P2 only",
  "P3 only"
)

overlap_count <- as.data.frame(
  table(
    factor(
      overlap$exact_overlap,
      levels = overlap_levels
    ),
    factor(
      toupper(
        overlap$view_code
      ),
      levels = toupper(
        VIEW_LEVELS
      )
    )
  ),
  stringsAsFactors = FALSE
)

names(overlap_count) <- c(
  "Exact_overlap",
  "Omic",
  "n_features"
)

overlap_count <- overlap_count[
  overlap_count$n_features > 0L,
  ,
  drop = FALSE
]

p_overlap <- ggplot2::ggplot(
  overlap_count,
  ggplot2::aes(
    x = Exact_overlap,
    y = n_features,
    fill = Omic
  )
) +
  ggplot2::geom_col() +
  ggplot2::labs(
    title =
      "P1/P2/P3 panel overlap — audit only",
    subtitle =
      "Each question-specific panel was defined independently before this comparison.",
    x =
      "Exact panel overlap",
    y =
      "Number of strict Evidence-A features",
    fill =
      "Omic"
  ) +
  theme_final(10) +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 30,
        hjust = 1
      )
  )

PLOTS_FINAL$overlap_audit <-
  p_overlap

print(
  p_overlap
)

cat("\n============================================================\n")
cat("PLOTS SENT TO RSTUDIO PLOTS HISTORY\n")
cat("  P1/P2/P3 directionality\n")
cat("  clinical features by question\n")
cat("  P1/P2/P3 pathway summaries\n")
cat("  P1/P2/P3 feature<->pathway maps\n")
cat("  cross-question overlap audit\n")
cat("\nNamed plot objects are in PLOTS_FINAL.\n")
cat("Example: print(PLOTS_FINAL$P1_direction)\n")
cat("Example: print(PLOTS_FINAL$P1_feature_pathway)\n")
cat("Example: print(PLOTS_FINAL$clinical_by_question)\n")
cat("NO FILES WERE WRITTEN.\n")
cat("============================================================\n")

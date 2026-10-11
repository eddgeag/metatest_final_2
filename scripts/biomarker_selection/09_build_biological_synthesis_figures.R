#!/usr/bin/env Rscript

# ===== METATEST — 09_build_biological_synthesis_figures.R =====
#
# PREGUNTAS BIOLÓGICAS
#
# 1. ¿Qué significan F1 y F2?
# 2. ¿Qué pathways caracterizan NP, SO y FD en cada pregunta?
# 3. ¿Qué pathways representan programas transcriptómicos redundantes comunes?
# 4. ¿Qué genes sostienen esos programas?
# 5. ¿Qué features PR, ME y CL se asocian a los pathways de cada programa TX?
# 6. ¿Qué programas y pathways son compartidos, específicos o
#    multimodalmente reforzados entre P1/P2/P3 y los lados fenotípicos?
# 7. ¿Cómo se manifiestan los pathways/programas y su arquitectura
#    multimodal individuo por individuo?
#
# Este script se construye secuencialmente siguiendo esas siete preguntas.
# TRAIN define e interpreta.
# TEST no redefine factores, pathways ni programas.

options(stringsAsFactors = FALSE)


# ===== 0. CONFIGURACIÓN =====


UPSTREAM_DIR <- Sys.getenv(
  "UPSTREAM_DIR",
  unset = path.expand(
    "~/Documents/PhD/metatest/results_upstream_local_20260908/comparacion_modelamientos/upstream_comun"
  )
)

# Opcionales. Si se dejan vacíos, se localiza automáticamente la última corrida.
SELECTION_RUN   <- Sys.getenv("SELECTION_RUN", unset = "")
BIO_PROGRAM_RUN <- Sys.getenv("BIO_PROGRAM_RUN", unset = "")

FAMILIES <- trimws(strsplit(
  Sys.getenv("BIO_FIG_FAMILIES", unset = "P1,P2,P3"),
  ",", fixed = TRUE
)[[1]])
FAMILIES <- FAMILIES[nzchar(FAMILIES)]

FIG_DPI <- as.integer(Sys.getenv("BIO_FIG_DPI", unset = "180"))
MAKE_PNG <- tolower(Sys.getenv("BIO_FIG_MAKE_PNG", unset = "true")) %in%
  c("true", "1", "yes", "y")

# Mostrar ID de muestra junto a cada punto. Con n=16 es legible y útil para
# detectar visualmente muestras extremas. Se puede desactivar con false.
SHOW_SAMPLE_LABELS <- tolower(
  Sys.getenv("BIO_FIG1_SHOW_SAMPLE_LABELS", unset = "true")
) %in% c("true", "1", "yes", "y")

# ===== CONFIGURACIÓN VISUAL GLOBAL =====

FIG_TEXT_SCALE <- as.numeric(
  Sys.getenv(
    "BIO_FIG_TEXT_SCALE",
    unset = "1.15"
  )
)


# ===== CONFIGURACIÓN GRÁFICA DE LA PREGUNTA 4 =====

# Si hay <= este número de drivers, se muestran todos.
# Si hay más, se reduce usando conjuntamente:
#   - normalized_anchor_weight alto
#   - |loading_oriented_scaled| alto
FIG4_PRESENTATION_MAX_GENES <- as.integer(
  Sys.getenv(
    "BIO_FIG4_PRESENTATION_MAX_GENES",
    unset = "15"
  )
)

# Percentiles gráficos. 0.50 = mitad superior.
FIG4_PRESENTATION_DRIVER_PERCENTILE_MIN <- as.numeric(
  Sys.getenv(
    "BIO_FIG4_PRESENTATION_DRIVER_PERCENTILE_MIN",
    unset = "0.50"
  )
)

FIG4_PRESENTATION_LOADING_PERCENTILE_MIN <- as.numeric(
  Sys.getenv(
    "BIO_FIG4_PRESENTATION_LOADING_PERCENTILE_MIN",
    unset = "0.50"
  )
)

FIG4_PRESENTATION_BASE_SIZE <- as.numeric(
  Sys.getenv(
    "BIO_FIG4_PRESENTATION_BASE_SIZE",
    unset = "12"
  )
)

FIG4_PRESENTATION_AXIS_TEXT_SIZE <- as.numeric(
  Sys.getenv(
    "BIO_FIG4_PRESENTATION_AXIS_TEXT_SIZE",
    unset = "10"
  )
)


# ===== CONFIGURACIÓN GRÁFICA DE LA PREGUNTA 5 =====

FIG5_PRESENTATION_MAX_FEATURES <- as.integer(
  Sys.getenv(
    "BIO_FIG5_PRESENTATION_MAX_FEATURES",
    unset = "25"
  )
)

# Reserva mínima por vista cuando hay que reducir.
# Evita que PR/CL consuman todos los huecos y desaparezcan TX/ME.
FIG5_PRESENTATION_MIN_PER_VIEW <- as.integer(
  Sys.getenv(
    "BIO_FIG5_PRESENTATION_MIN_PER_VIEW",
    unset = "3"
  )
)

FIG5A_BASE_SIZE <- as.numeric(
  Sys.getenv(
    "BIO_FIG5A_BASE_SIZE",
    unset = "12"
  )
)

FIG5A_AXIS_TEXT_SIZE <- as.numeric(
  Sys.getenv(
    "BIO_FIG5A_AXIS_TEXT_SIZE",
    unset = "10"
  )
)

FIG5A_STRIP_TEXT_SIZE <- as.numeric(
  Sys.getenv(
    "BIO_FIG5A_STRIP_TEXT_SIZE",
    unset = "10"
  )
)

FIG5A_LABEL_WRAP <- as.integer(
  Sys.getenv(
    "BIO_FIG5A_LABEL_WRAP",
    unset = "18"
  )
)

FIG5A_INCH_PER_FEATURE <- as.numeric(
  Sys.getenv(
    "BIO_FIG5A_INCH_PER_FEATURE",
    unset = "0.38"
  )
)

# Número máximo de features en las matrices de PRESENTACIÓN.
#
# Si una vía tiene <= este número:
#   -> se muestran TODAS.
#
# Si tiene > este número:
#   -> PR y CL se conservan prioritariamente;
#   -> TX y ME se reducen por variabilidad descriptiva en los 26 individuos;
#   -> para TX se priorizan primero los genes driver;
#   -> la selección es EXCLUSIVAMENTE gráfica.
#
# 01_ALL nunca se filtra.
# ===== CONFIGURACIÓN VISUAL FIG5B =====

FIG5B_CLUSTER_WITHIN_VIEW <- TRUE

FIG5B_CLUSTER_METHOD <- Sys.getenv(
  "BIO_FIG5B_CLUSTER_METHOD",
  unset = "average"
)

FIG5B_BASE_SIZE <- as.numeric(
  Sys.getenv(
    "BIO_FIG5B_BASE_SIZE",
    unset = "14"
  )
)

FIG5B_AXIS_TEXT_SIZE <- as.numeric(
  Sys.getenv(
    "BIO_FIG5B_AXIS_TEXT_SIZE",
    unset = "14"
  )
)

FIG5B_CELL_TEXT_SIZE <- as.numeric(
  Sys.getenv(
    "BIO_FIG5B_CELL_TEXT_SIZE",
    unset = "3.2"
  )
)
# ===== CONFIGURACIÓN GRÁFICA DE LA PREGUNTA 7 =====

# Límite exclusivamente visual para los Z-scores.
# Los CSV conservan los valores completos.
FIG7_Z_LIMIT <- as.numeric(
  Sys.getenv(
    "BIO_FIG7_Z_LIMIT",
    unset = "3"
  )
)


# Mostrar IDs individuales en FIG7C.
FIG7_SHOW_SAMPLE_LABELS <- tolower(
  Sys.getenv(
    "BIO_FIG7_SHOW_SAMPLE_LABELS",
    unset = "true"
  )
) %in% c(
  "true",
  "1",
  "yes",
  "y"
)

# Permutaciones para la separación multimodal descriptiva.
FIG7_PERMANOVA_NPERM <- as.integer(
  Sys.getenv(
    "BIO_FIG7_PERMANOVA_NPERM",
    unset = "999"
  )
)

FIG7_PERMANOVA_SEED <- as.integer(
  Sys.getenv(
    "BIO_FIG7_PERMANOVA_SEED",
    unset = "20261001"
  )
)

if (
  !is.finite(FIG7_PERMANOVA_NPERM) ||
  FIG7_PERMANOVA_NPERM < 99L
) {
  stop(
    "BIO_FIG7_PERMANOVA_NPERM must be >= 99",
    call. = FALSE
  )
}
if (
  !is.finite(FIG7_Z_LIMIT) ||
  FIG7_Z_LIMIT <= 0
) {
  stop(
    "BIO_FIG7_Z_LIMIT must be > 0",
    call. = FALSE
  )
}

if (
  !is.finite(FIG5_PRESENTATION_MAX_FEATURES) ||
  FIG5_PRESENTATION_MAX_FEATURES < 4L
) {
  stop(
    "BIO_FIG5_PRESENTATION_MAX_FEATURES must be >= 4",
    call. = FALSE
  )
}

# ===== CONTROL DE REANUDACIÓN =====

BIO_FIG_RESUME_DIR <- Sys.getenv(
  "BIO_FIG_RESUME_DIR",
  unset = ""
)

BIO_FIG_START_QUESTION <- as.integer(
  Sys.getenv(
    "BIO_FIG_START_QUESTION",
    unset = "1"
  )
)

BIO_FIG_FORCE_RAW <- toupper(
  trimws(
    Sys.getenv(
      "BIO_FIG_FORCE",
      unset = ""
    )
  )
)

BIO_FIG_FORCE <- trimws(
  strsplit(
    BIO_FIG_FORCE_RAW,
    ",",
    fixed = TRUE
  )[[1]]
)

BIO_FIG_FORCE <- BIO_FIG_FORCE[
  nzchar(BIO_FIG_FORCE)
]

BIO_FIG_RESUME <- nzchar(
  BIO_FIG_RESUME_DIR
)

if (
  !BIO_FIG_START_QUESTION %in%
  c(1L, 5L, 6L, 7L)
) {
  stop(
    paste0(
      "BIO_FIG_START_QUESTION currently supports 1, 5, 6 or 7. ",
      "Use 1 for a full run, 5 to resume Question 5, ",
      "6 to build Question 6, or 7 to build Question 7."
    ),
    call. = FALSE
  )
}


if (
  BIO_FIG_START_QUESTION %in%
  c(5L, 6L, 7L) &&
  !BIO_FIG_RESUME
) {
  stop(
    paste0(
      "BIO_FIG_START_QUESTION=5, 6 or 7 requires ",
      "BIO_FIG_RESUME_DIR."
    ),
    call. = FALSE
  )
}

#===== 1. DEPENDENCIAS=====


if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "Package 'ggplot2' is required. Install it before running script 09.",
    call. = FALSE
  )
}

# =====2. HELPERS GENERALES =====


assert <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
}
run_question <- function(q) {
  q >= BIO_FIG_START_QUESTION
}


force_requested <- function(...) {
  
  tags <- toupper(
    as.character(
      c(...)
    )
  )
  
  if (!length(BIO_FIG_FORCE)) {
    return(FALSE)
  }
  
  any(
    BIO_FIG_FORCE %in%
      c(
        "ALL",
        "Q5",
        tags
      )
  )
}


files_complete <- function(paths) {
  
  paths <- as.character(paths)
  
  if (
    !length(paths) ||
    !all(file.exists(paths))
  ) {
    return(FALSE)
  }
  
  info <- file.info(paths)
  
  all(
    is.finite(info$size) &
      info$size > 0
  )
}
canonical_path <- function(p, must_work = TRUE) {
  normalizePath(path.expand(p), winslash = "/", mustWork = must_work)
}

read_csv_checked <- function(path, required = character()) {
  assert(file.exists(path), paste("Missing required file:", path))
  
  x <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  missing <- setdiff(required, names(x))
  
  assert(
    !length(missing),
    paste(
      "Missing columns in", basename(path), ":",
      paste(missing, collapse = ", ")
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

write_text <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, con = path, useBytes = TRUE)
}

factor_number <- function(x) {
  y <- suppressWarnings(as.integer(gsub("[^0-9]+", "", as.character(x))))
  y[is.na(y)] <- seq_len(sum(is.na(y))) + 10000L
  y
}

phenotype_label <- function(x) {
  # Solo mejora la legibilidad gráfica. NO modifica la categoría guardada por 08.
  x <- as.character(x)
  x[x == "NP_side"]        <- "NP-side"
  x[x == "SO_side"]        <- "SO-side"
  x[x == "FD_side"]        <- "FD-side"
  x[x == "nonFD_side"]     <- "non-FD-side"
  x[x == "ambiguous_side"] <- "ambiguous-side"
  x
}

axis_title_from_side_map <- function(side_map, factor_name) {
  z <- side_map[side_map$factor == factor_name, , drop = FALSE]
  assert(nrow(z) == 2L, paste("Expected two sides for", factor_name))

  neg <- z[z$oriented_side == "negative", , drop = FALSE]
  pos <- z[z$oriented_side == "positive", , drop = FALSE]
  assert(nrow(neg) == 1L && nrow(pos) == 1L,
         paste("Could not identify negative/positive sides for", factor_name))

  paste0(
    phenotype_label(neg$phenotypic_side[1L]),
    "  ←  ", factor_name, "  →  ",
    phenotype_label(pos$phenotypic_side[1L])
  )
}

side_caption_from_map <- function(side_map, factor_name) {
  z <- side_map[side_map$factor == factor_name, , drop = FALSE]
  neg <- z[z$oriented_side == "negative", , drop = FALSE]
  pos <- z[z$oriented_side == "positive", , drop = FALSE]

  fmt_prob <- function(v) {
    if (!length(v) || is.na(v)) return("NA")
    sprintf("%.2f", as.numeric(v))
  }

  paste0(
    factor_name, ": − = ", phenotype_label(neg$phenotypic_side[1L]),
    " [nominal ", neg$assigned_class_nominal[1L],
    ", P=", fmt_prob(neg$p_assigned_class[1L]), "]",
    "; + = ", phenotype_label(pos$phenotypic_side[1L]),
    " [nominal ", pos$assigned_class_nominal[1L],
    ", P=", fmt_prob(pos$p_assigned_class[1L]), "]"
  )
}

save_plot <- function(plot, pdf_path, png_path = NULL, width = 8.2, height = 7.0) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    limitsize = FALSE
  )

  if (MAKE_PNG && !is.null(png_path)) {
    ggplot2::ggsave(
      filename = png_path,
      plot = plot,
      width = width,
      height = height,
      dpi = FIG_DPI,
      limitsize = FALSE
    )
  }
}
wrap_plot_text <- function(
    x,
    width = 85L
) {
  paste(
    strwrap(
      as.character(x),
      width = width
    ),
    collapse = "\n"
  )
}
fig_theme <- function(base_size = 11) {
  
  base_size <- base_size * FIG_TEXT_SCALE
  
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(size = base_size),
      plot.caption = ggplot2::element_text(size = base_size - 2, hjust = 0),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}


# =====3. LOCALIZAR LA CORRIDA DE 08 =====

find_latest_selection_run <- function(upstream) {
  base <- file.path(
    upstream,
    "downstream",
    "seleccion_mejor_escenario_train"
  )
  assert(dir.exists(base), paste("Selection root does not exist:", base))

  dirs <- list.dirs(base, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "SELECTED_SCENARIOS_BY_ANALYSIS.csv"))]
  assert(length(dirs) > 0L, paste("No complete selection runs found in", base))

  info <- file.info(dirs)
  canonical_path(
    dirs[order(info$mtime, basename(dirs), decreasing = TRUE)][1L],
    TRUE
  )
}

find_latest_bio_program_run <- function(selection_run) {
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
      return(canonical_path(p, TRUE))
    }
  }

  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "BIOLOGICAL_PROGRAMS_OBJECT.rds"))]
  assert(length(dirs) > 0L, paste("No complete biological-program runs in", root))

  info <- file.info(dirs)
  canonical_path(
    dirs[order(info$mtime, basename(dirs), decreasing = TRUE)][1L],
    TRUE
  )
}

resolve_bio_program_run <- function() {
  if (nzchar(BIO_PROGRAM_RUN)) {
    return(canonical_path(BIO_PROGRAM_RUN, TRUE))
  }

  selection_run <- if (nzchar(SELECTION_RUN)) {
    canonical_path(SELECTION_RUN, TRUE)
  } else {
    find_latest_selection_run(UPSTREAM_DIR)
  }

  find_latest_bio_program_run(selection_run)
}


# ===== PREGUNTA 1. ¿QUÉ SIGNIFICAN F1 Y F2? =====

# ===== 4.1 CARGAR AUDITORÍA CLASE ↔ FACTOR =====
load_fig1_inputs <- function(bio_program_run, family) {
  family_dir <- file.path(bio_program_run, family)
  assert(dir.exists(family_dir), paste("Missing family output from 08:", family_dir))

  audit <- read_csv_checked(
    file.path(family_dir, "00_PROGRAM_DEFINITION_AUDIT.csv"),
    required = c("family", "latent_run")
  )

  side_map <- read_csv_checked(
    file.path(family_dir, "01_FACTOR_PHENOTYPIC_SIDE_MAP.csv"),
    required = c(
      "factor", "oriented_side",
      "assigned_class_nominal", "p_assigned_class",
      "assignment_strong", "factor_role", "phenotypic_side"
    )
  )

  assert(nrow(audit) >= 1L, paste("Empty audit table for", family))
  latent_run <- as.character(audit$latent_run[1L])
  assert(dir.exists(latent_run), paste("Latent run stored by 08 does not exist:", latent_run))

  scores_path <- file.path(
    latent_run,
    "tables",
    "02_LATENT_SCORES_TRAIN_AND_PROJECTED_TEST.csv"
  )

  scores <- read_csv_checked(
    scores_path,
    required = c("sample", "split", "group")
  )

  # ---------------------------------------------------------------------------
  # Normalizar únicamente la etiqueta del split para visualización.
  #
  # Acepta TRAIN, TEST, PROJECTED_TEST, TEST_PROJECTED, etc.
  # No modifica scores ni realiza ninguna nueva proyección.
  # ---------------------------------------------------------------------------
  
  scores$split_original <- as.character(scores$split)
  
  split_upper <- toupper(
    trimws(scores$split_original)
  )
  
  scores$split <- ifelse(
    grepl("TRAIN", split_upper),
    "TRAIN",
    ifelse(
      grepl("TEST", split_upper),
      "TEST",
      NA_character_
    )
  )
  
  
  train <- scores[
    scores$split == "TRAIN",
    ,
    drop = FALSE
  ]
  
  test <- scores[
    scores$split == "TEST",
    ,
    drop = FALSE
  ]
  
  
  assert(
    nrow(train) > 0L,
    paste(
      "No TRAIN rows detected for",
      family,
      "| original split values:",
      paste(
        unique(scores$split_original),
        collapse = ", "
      )
    )
  )
  
  
  assert(
    nrow(test) > 0L,
    paste(
      "No projected TEST rows detected for",
      family,
      "| original split values:",
      paste(
        unique(scores$split_original),
        collapse = ", "
      )
    )
  )
  list(
    family = family,
    latent_run = latent_run,
    side_map = side_map,
    train_scores = train,
    test_scores = test,
    all_scores = scores
  )
}

# ===== 4.2 FIGURA 1 — ESPACIO LATENTE TRAIN =====


make_figure1 <- function(inp, outdir) {
  family <- inp$family
  side_map <- inp$side_map
  train <- inp$train_scores

  # Orden reproducible de factores. Para esta figura necesitamos los dos primeros.
  factors <- unique(as.character(side_map$factor))
  factors <- factors[order(factor_number(factors), factors)]

  assert(
    length(factors) >= 2L,
    paste("Figure 1 requires at least two audited factors in", family)
  )

  f1 <- factors[1L]
  f2 <- factors[2L]

  xcol <- paste0(f1, "__oriented_scaled")
  ycol <- paste0(f2, "__oriented_scaled")

  assert(xcol %in% names(train), paste("Missing score column:", xcol))
  assert(ycol %in% names(train), paste("Missing score column:", ycol))

  dat <- data.frame(
    family = family,
    sample = as.character(train$sample),
    group = as.character(train$group),
    factor_x = f1,
    factor_y = f2,
    score_x = suppressWarnings(as.numeric(train[[xcol]])),
    score_y = suppressWarnings(as.numeric(train[[ycol]])),
    split = "TRAIN",
    stringsAsFactors = FALSE
  )

  assert(all(is.finite(dat$score_x)), paste("Non-finite", f1, "scores in", family))
  assert(all(is.finite(dat$score_y)), paste("Non-finite", f2, "scores in", family))

  dat$group <- factor(dat$group, levels = c("NP", "FD", "SO"))

  x_title <- axis_title_from_side_map(side_map, f1)
  y_title <- axis_title_from_side_map(side_map, f2)

  caption_lines <- c(
    side_caption_from_map(side_map, f1),
    side_caption_from_map(side_map, f2),
    "La asignación fenotípica procede del audit clase↔factor congelado en 08; el signo matemático por sí solo no se interpreta biológicamente."
  )

  # Centroides descriptivos de TRAIN: ayudan a ver dónde se sitúa cada clase,
  # pero NO intervienen en ninguna selección ni redefinición del factor.
  centroids <- stats::aggregate(
    cbind(score_x, score_y) ~ group,
    data = dat,
    FUN = mean
  )

  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(x = score_x, y = score_y, color = group)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.45) +
    ggplot2::geom_point(size = 3.4, alpha = 0.90) +
    ggplot2::geom_point(
      data = centroids,
      ggplot2::aes(x = score_x, y = score_y, color = group),
      shape = 4,
      stroke = 1.5,
      size = 5,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      title = paste0(
        family,
        " — Pregunta 1: interpretación del espacio latente MOFA"
      ),
      subtitle = paste0(
        "TRAIN. Puntos = individuos; × = centroide descriptivo de cada clase. ",
        "F1 y F2 se interpretan mediante sus lados fenotípicos."
      ),
      x = x_title,
      y = y_title,
      color = "Clase",
      caption = paste(caption_lines, collapse = "\n")
    ) +
    fig_theme(11)

  if (SHOW_SAMPLE_LABELS) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = sample),
        size = 2.6,
        vjust = -0.8,
        show.legend = FALSE,
        check_overlap = TRUE
      )
  }

  # Datos exactos de la figura.
  centroids$family <- family
  centroids$record_type <- "class_centroid_TRAIN"
  dat$record_type <- "individual_TRAIN"

  out_dat <- dat[, c(
    "family", "record_type", "sample", "group",
    "factor_x", "factor_y", "score_x", "score_y", "split"
  )]

  write_csv(
    out_dat,
    file.path(outdir, "plot_data", paste0("FIG1_", family, "_LATENT_FACTOR_CLASS_MAP.csv"))
  )

  ggplot2::ggsave(
    filename = file.path(
      outdir,
      "figures",
      paste0("FIG1_", family, "_LATENT_FACTOR_CLASS_MAP.png")
    ),
    plot = p,
    width = 8.5,
    height = 7.2,
    dpi = FIG_DPI,
    limitsize = FALSE
  )
  list(
    plot = p,
    data = dat,
    centroids = centroids,
    factor_x = f1,
    factor_y = f2,
    x_title = x_title,
    y_title = y_title
  )
}
# ===== 4.3 FIGURA 1B — AUDITORÍA DESCRIPTIVA TRAIN + TEST PROYECTADO =====

make_figure1b_test_audit <- function(inp, outdir) {
  
  family <- inp$family
  side_map <- inp$side_map
  scores <- inp$all_scores
  
  factors <- unique(as.character(side_map$factor))
  factors <- factors[order(factor_number(factors), factors)]
  
  assert(
    length(factors) >= 2L,
    paste("Figure 1B requires at least two audited factors in", family)
  )
  
  f1 <- factors[1L]
  f2 <- factors[2L]
  
  xcol <- paste0(f1, "__oriented_scaled")
  ycol <- paste0(f2, "__oriented_scaled")
  
  assert(xcol %in% names(scores), paste("Missing score column:", xcol))
  assert(ycol %in% names(scores), paste("Missing score column:", ycol))
  
  
  dat <- data.frame(
    family = family,
    sample = as.character(scores$sample),
    group = as.character(scores$group),
    split = as.character(scores$split),
    factor_x = f1,
    factor_y = f2,
    score_x = suppressWarnings(as.numeric(scores[[xcol]])),
    score_y = suppressWarnings(as.numeric(scores[[ycol]])),
    stringsAsFactors = FALSE
  )
  
  dat <- dat[
    dat$split %in% c("TRAIN", "TEST"),
    ,
    drop = FALSE
  ]
  
  dat$group <- factor(
    dat$group,
    levels = c("NP", "FD", "SO")
  )
  
  dat$split <- factor(
    dat$split,
    levels = c("TRAIN", "TEST")
  )
  
  
  # Centroides SIEMPRE calculados exclusivamente con TRAIN.
  train_dat <- dat[
    dat$split == "TRAIN",
    ,
    drop = FALSE
  ]
  
  centroids <- stats::aggregate(
    cbind(score_x, score_y) ~ group,
    data = train_dat,
    FUN = mean
  )
  
  
  x_title <- axis_title_from_side_map(
    side_map,
    f1
  )
  
  y_title <- axis_title_from_side_map(
    side_map,
    f2
  )
  
  
  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = score_x,
      y = score_y,
      color = group,
      shape = split
    )
  ) +
    
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = 2,
      linewidth = 0.45
    ) +
    
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = 2,
      linewidth = 0.45
    ) +
    
    ggplot2::geom_point(
      size = 3.5,
      stroke = 1.1
    ) +
    
    ggplot2::geom_point(
      data = centroids,
      ggplot2::aes(
        x = score_x,
        y = score_y,
        color = group
      ),
      shape = 4,
      stroke = 1.5,
      size = 5,
      inherit.aes = FALSE
    ) +
    
    ggplot2::scale_shape_manual(
      values = c(
        "TRAIN" = 16,
        "TEST" = 1
      )
    ) +
    
    ggplot2::labs(
      title = paste0(
        family,
        " — auditoría descriptiva del espacio latente"
      ),
      
      subtitle = paste0(
        "TRAIN = sólido; TEST proyectado = abierto; ",
        "× = centroide calculado exclusivamente en TRAIN."
      ),
      
      x = x_title,
      y = y_title,
      
      color = "Clase",
      shape = "Split",
      
      caption = paste0(
        "TEST se proyecta sobre el espacio latente ya aprendido. ",
        "No se usa para orientar factores, redefinir clases, ",
        "seleccionar pathways ni recalibrar el modelo."
      )
    ) +
    
    fig_theme(11)
  
  
  ggplot2::ggsave(
    filename = file.path(
      outdir,
      "figures",
      paste0(
        "FIG1B_",
        family,
        "_TRAIN_TEST_LATENT_AUDIT.png"
      )
    ),
    plot = p,
    width = 8.5,
    height = 7.2,
    dpi = FIG_DPI,
    limitsize = FALSE
  )
  
  
  write_csv(
    dat,
    file.path(
      outdir,
      "plot_data",
      paste0(
        "FIG1B_",
        family,
        "_TRAIN_TEST_LATENT_AUDIT.csv"
      )
    )
  )
  
  
  list(
    plot = p,
    data = dat,
    centroids_train = centroids,
    factor_x = f1,
    factor_y = f2
  )
}
# ===== 4.2 FIGURA 1 — ESPACIO LATENTE TRAIN =====


BIO_RUN <- resolve_bio_program_run()

if (BIO_FIG_RESUME) {
  
  OUT_ROOT <- canonical_path(
    BIO_FIG_RESUME_DIR,
    must_work = TRUE
  )
  
  STAMP <- basename(
    OUT_ROOT
  )
  
  message(
    "RESUME MODE: ",
    OUT_ROOT
  )
  
} else {
  
  STAMP <- format(
    Sys.time(),
    "%Y%m%dT%H%M%S"
  )
  
  OUT_ROOT <- file.path(
    BIO_RUN,
    "BIOLOGICAL_SYNTHESIS_FIGURES",
    paste0(
      STAMP,
      "_pid",
      Sys.getpid()
    )
  )
}

DIR_FIG <- file.path(OUT_ROOT, "figures")
DIR_DAT <- file.path(OUT_ROOT, "plot_data")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_DAT, recursive = TRUE, showWarnings = FALSE)
latest_dir <- file.path(
  BIO_RUN,
  "BIOLOGICAL_SYNTHESIS_FIGURES"
)

dir.create(
  latest_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

writeLines(
  canonical_path(
    OUT_ROOT,
    must_work = TRUE
  ),
  file.path(
    latest_dir,
    "LATEST_INCOMPLETE_RUN.txt"
  )
)

writeLines(
  paste(
    "RUNNING",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ),
  file.path(
    OUT_ROOT,
    "RUNNING.flag"
  )
)

if (
  file.exists(
    file.path(
      OUT_ROOT,
      "RUN_COMPLETE.flag"
    )
  )
) {
  unlink(
    file.path(
      OUT_ROOT,
      "RUN_COMPLETE.flag"
    )
  )
}


message("============================================================")
message("METATEST — biological synthesis figures | PART 1")
message("BIOLOGICAL_PROGRAMS run: ", BIO_RUN)
message("Output: ", OUT_ROOT)
message("Building Questions 1–7")
message("============================================================")

fig1_results <- list()
side_rows <- list()
manifest_rows <- list()
rr <- 1L

fig1b_results <- list()
if (run_question(1L)) {
for (family in FAMILIES) {
  message("\n[FIG1] ", family)

  inp <- load_fig1_inputs(BIO_RUN, family)
  res <- make_figure1(inp, OUT_ROOT)
  fig1_results[[family]] <- res
  res_audit <- make_figure1b_test_audit(
    inp,
    OUT_ROOT
  )
  
  fig1b_results[[family]] <- res_audit


  sm <- inp$side_map
  sm$family <- family
  side_rows[[family]] <- sm

  manifest_rows[[rr]] <- data.frame(
    figure = "FIG1",
    family = family,
    question = "How do NP, FD and SO occupy the first two latent factors and what phenotype does each factor side represent?",
    data_split = "TRAIN only",
    source_side_map = file.path(BIO_RUN, family, "01_FACTOR_PHENOTYPIC_SIDE_MAP.csv"),
    source_latent_run = inp$latent_run,
    factor_x = res$factor_x,
    factor_y = res$factor_y,
    png = if (MAKE_PNG) file.path(DIR_FIG, paste0("FIG1_", family, "_LATENT_FACTOR_CLASS_MAP.png")) else NA_character_,
    analytical_decision_added = FALSE,
    stringsAsFactors = FALSE
  )
  rr <- rr + 1L
}



# Tabla conjunta que deja explícita la traducción signo matemático -> lado fenotípico.
all_sides <- do.call(rbind, side_rows)
write_csv(
  all_sides,
  file.path(DIR_DAT, "FIG1_FACTOR_SIDE_LABELS.csv")
)

manifest <- do.call(rbind, manifest_rows)
write_csv(manifest, file.path(OUT_ROOT, "00_FIGURE_BUILD_MANIFEST.csv"))


  write_text(
    c(
      "METATEST — biological synthesis figures",
      "",
      "PREGUNTA 1 — ¿Qué significan F1 y F2?",
      "  - F1 x F2 en TRAIN.",
      "  - La identidad biológica procede de la auditoría clase ↔ factor.",
      "  - El signo matemático no define la identidad fenotípica.",
      "",
      "PREGUNTA 2 — ¿Qué pathways caracterizan cada lado fenotípico?",
      "  - Solo pathways TX primarios con q <= 0.10.",
      "  - P1, P2 y P3 se mantienen separados.",
      "  - Se conserva un overview y todos los pathways significativos.",
      "",
      "PREGUNTA 3 — ¿Qué pathways representan programas TX redundantes?",
      "  - Los pathways se agrupan usando la definición congelada en 08.",
      "  - Se evalúa conservación entre P1/P2/P3 y entre lados fenotípicos.",
      "",
      "PREGUNTA 4 — ¿Qué genes sostienen esos programas?",
      "  - Se conservan todos los genes TX observados de los pathways primarios.",
      "  - Los drivers constituyen un subconjunto de los genes observados.",
      "  - La fuerza driver se muestra mediante normalized_anchor_weight congelado en TRAIN.",
      "  - normalized_anchor_weight no es logFC ni un nuevo efecto diferencial.",
      "  - Se auditan genes compartidos entre pathways, programas y análisis.",
      "",
      "TRAIN define la interpretación.",
      "TEST no redefine factores, pathways, programas ni genes."   ,   "",
      "TRAIN define la interpretación.",
      "TEST no interviene en la definición actual."
    ),
    file.path(OUT_ROOT, "READ_ME_FIRST.txt")
  )



message("\n============================================================")
message("FIGURE 1 COMPLETE")
message("Output: ", OUT_ROOT)
message("============================================================")
}
# ===== PREGUNTA 2. ¿QUÉ PATHWAYS CARACTERIZAN NP, SO Y FD EN CADA PREGUNTA? =====

# ===== 5.1 CARGAR PATHWAYS TX PRIMARIOS SIGNIFICATIVOS =====


#
# ESTRUCTURA:
#
#   P1/
#     Factor1_NP_side/
#       GO_BP/
#         ALL_SIGNIFICANT.png
#         TOP_DISPLAY.png
#       REACTOME/
#       ...
#
# PRINCIPIOS:
#   - solo transcriptómica;
#   - solo pathways primarios significativos q <= 0.10;
#   - P1/P2/P3 separados;
#   - factor/lado fenotípico separados;
#   - librerías separadas;
#   - ALL_SIGNIFICANT conserva todos los pathways significativos;
#   - no se seleccionan top pathways en esta pregunta;
#   - la reducción a programas y pathways representativos pertenece a la Pregunta 3.
#
# ==============================================================================


clean_pathway_label <- function(pathway) {
  
  x <- as.character(pathway)
  
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^REACTOME_", "", x)
  x <- sub("^GOBP_", "", x)
  x <- sub("^GOCC_", "", x)
  x <- sub("^GOMF_", "", x)
  x <- sub("^KEGG_MEDICUS_", "", x)
  
  gsub("_", " ", x, fixed = TRUE)
}


short_program_id <- function(program_id) {
  
  x <- as.character(program_id)
  
  out <- sub("^.*__(P[0-9]+)$", "\\1", x)
  
  bad <- out == x
  
  if (any(bad)) {
    out[bad] <- paste0("Prog", seq_len(sum(bad)))
  }
  
  out
}


safe_filename <- function(x) {
  
  x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  
  x
}


load_fig2_inputs <- function(bio_program_run, family) {
  
  family_dir <- file.path(
    bio_program_run,
    family
  )
  
  # ---------------------------------------------------------------------------
  # Pathways primarios: salida ANTERIOR a la definición de programas TX
  # ---------------------------------------------------------------------------
  
  path <- file.path(
    family_dir,
    "02_PRIMARY_PATHWAYS_WITH_PHENOTYPIC_SIDE.csv"
  )
  
  x <- read_csv_checked(
    path,
    required = c(
      "pathway",
      "family",
      "library",
      "factor",
      "oriented_side",
      "set_statistic",
      "q_mofa_within_side",
      "n_features_in_set",
      "q_both_sides_within_library_factor",
      "significant_both_sides_adjusted",
      "factor_role",
      "phenotypic_side",
      "p_assigned_class",
      "assignment_strong"
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Mapa completo de lados fenotípicos.
  #
  # Es necesario porque puede existir un lado fenotípico real
  # aunque tenga CERO pathways significativos.
  #
  # Ejemplo actual:
  #   P2 / Factor1 / NP_side = 0 pathways significativos
  # ---------------------------------------------------------------------------
  
  side_map <- read_csv_checked(
    file.path(
      family_dir,
      "01_FACTOR_PHENOTYPIC_SIDE_MAP.csv"
    ),
    required = c(
      "factor",
      "oriented_side",
      "factor_role",
      "phenotypic_side"
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Criterio primario de significación
  # ---------------------------------------------------------------------------
  
  x <- x[
    !is.na(x$q_both_sides_within_library_factor) &
      x$q_both_sides_within_library_factor <= 0.10 &
      x$significant_both_sides_adjusted,
    ,
    drop = FALSE
  ]
  
  
  x$family <- family
  
  x$minus_log10_q <- -log10(
    pmax(
      as.numeric(x$q_both_sides_within_library_factor),
      .Machine$double.xmin
    )
  )
  
  x$pathway_label <- clean_pathway_label(
    x$pathway
  )
  
  x$factor_side <- paste0(
    x$factor,
    "_",
    x$phenotypic_side
  )
  
  
  side_map$factor_side <- paste0(
    side_map$factor,
    "_",
    side_map$phenotypic_side
  )
  
  
  list(
    pathways = x,
    side_map = side_map
  )
}

# ===== 5.2 FIGURA 2A — OVERVIEW POR LADO FENOTÍPICO Y LIBRERÍA =====

make_fig2_overview <- function(dat, side_map, family, outdir) {
  
  # Librerías que produjeron al menos algún pathway significativo
  # en este análisis.
  libraries <- sort(unique(dat$library))
  
  # Todos los lados fenotípicos definidos por la auditoría,
  # tengan o no pathways significativos.
  sides <- unique(
    side_map[, c(
      "factor",
      "factor_role",
      "phenotypic_side",
      "factor_side"
    )]
  )
  
  
  # Todas las combinaciones lado × librería.
  grid <- merge(
    sides,
    data.frame(
      library = libraries,
      stringsAsFactors = FALSE
    ),
    by = NULL
  )
  
  
  # Conteos observados.
  counts <- stats::aggregate(
    pathway ~ factor_side + library,
    data = dat,
    FUN = length
  )
  
  names(counts)[
    names(counts) == "pathway"
  ] <- "n_pathways"
  
  
  # Añadir explícitamente los ceros.
  summary_dat <- merge(
    grid,
    counts,
    by = c(
      "factor_side",
      "library"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  
  summary_dat$n_pathways[
    is.na(summary_dat$n_pathways)
  ] <- 0L
  
  
  summary_dat$label <- as.character(
    summary_dat$n_pathways
  )
  
  
  p <- ggplot2::ggplot(
    summary_dat,
    ggplot2::aes(
      x = library,
      y = factor_side
    )
  ) +
    
    ggplot2::geom_tile(
      fill = "white",
      color = "grey60",
      linewidth = 0.5
    ) +
    
    ggplot2::geom_text(
      ggplot2::aes(
        label = label
      ),
      size = 4
    ) +
    
    ggplot2::labs(
      title = paste0(
        family,
        " — Pregunta 2: pathways TX significativos"
      ),
      
      subtitle = paste0(
        "Número de pathways con ",
        "q_both_sides_within_library_factor ≤ 0.10."
      ),
      
      x = "Librería",
      
      y = "Factor / lado fenotípico"
    ) +
    
    fig_theme(11) +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      
      legend.position = "none"
    )
  
  
  family_root <- file.path(
    outdir,
    "figures",
    "FIG2_PATHWAYS_TX",
    family
  )
  
  dir.create(
    family_root,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  ggplot2::ggsave(
    filename = file.path(
      family_root,
      paste0(
        "FIG2A_",
        family,
        "_OVERVIEW.png"
      )
    ),
    plot = p,
    width = 9,
    height = max(
      4,
      2.0 +
        0.7 *
        length(unique(summary_dat$factor_side))
    ),
    dpi = FIG_DPI,
    limitsize = FALSE
  )
  
  
  write_csv(
    summary_dat,
    file.path(
      family_root,
      paste0(
        "FIG2A_",
        family,
        "_OVERVIEW.csv"
      )
    )
  )
  
  
  invisible(summary_dat)
}

# ===== 5.3 FIGURA 2B — TODOS LOS PATHWAYS SIGNIFICATIVOS =====

make_fig2_block <- function(
    dat,
    family,
    factor_side,
    library,
    outdir
) {
  
  z <- dat[
    dat$factor_side == factor_side &
      dat$library == library,
    ,
    drop = FALSE
  ]
  
  
  block_dir <- file.path(
    outdir,
    "figures",
    "FIG2_PATHWAYS_TX",
    family,
    safe_filename(factor_side),
    safe_filename(library)
  )
  
  data_dir <- file.path(
    outdir,
    "plot_data",
    "FIG2_PATHWAYS_TX",
    family,
    safe_filename(factor_side),
    safe_filename(library)
  )
  
  dir.create(
    block_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    data_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # ---------------------------------------------------------------------------
  # Caso explícito: cero pathways significativos
  # ---------------------------------------------------------------------------
  
  if (!nrow(z)) {
    
    status <- data.frame(
      family = family,
      factor_side = factor_side,
      library = library,
      n_significant_pathways = 0L,
      criterion =
        "q_both_sides_within_library_factor <= 0.10",
      stringsAsFactors = FALSE
    )
    
    write_csv(
      status,
      file.path(
        data_dir,
        "NO_SIGNIFICANT_PATHWAYS.csv"
      )
    )
    
    return(
      invisible(
        list(
          all_significant = z
        )
      )
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Orden visual exclusivamente por evidencia estadística.
  # No hay programas ni representatives en esta pregunta.
  # ---------------------------------------------------------------------------
  
  z <- z[
    order(
      z$q_both_sides_within_library_factor,
      z$pathway
    ),
    ,
    drop = FALSE
  ]
  
  
  z$pathway_label_plot <- factor(
    z$pathway_label,
    levels = rev(
      z$pathway_label
    )
  )
  
  
  p_all <- ggplot2::ggplot(
    z,
    ggplot2::aes(
      x = minus_log10_q,
      y = pathway_label_plot
    )
  ) +
    
    ggplot2::geom_point(
      size = 3
    ) +
    
    ggplot2::labs(
      title = paste0(
        family,
        " — ",
        factor_side,
        " — ",
        library
      ),
      
      subtitle = paste0(
        "Pregunta 2 — todos los pathways TX ",
        "significativos (q ≤ 0.10)."
      ),
      
      x = expression(-log[10](q)),
      
      y = NULL,
      
      caption = paste0(
        "Cada punto representa un pathway significativo. ",
        "La reducción por redundancia se realiza en la Pregunta 3."
      )
    ) +
    
    fig_theme(10) +
    
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        size = 8
      ),
      
      panel.grid.minor =
        ggplot2::element_blank(),
      
      legend.position = "none"
    )
  
  
  height_all <- max(
    4.5,
    min(
      18,
      2.8 + 0.32 * nrow(z)
    )
  )
  
  
  ggplot2::ggsave(
    filename = file.path(
      block_dir,
      "ALL_SIGNIFICANT.png"
    ),
    plot = p_all,
    width = 10,
    height = height_all,
    dpi = FIG_DPI,
    limitsize = FALSE
  )
  
  
  write_csv(
    z,
    file.path(
      data_dir,
      "ALL_SIGNIFICANT.csv"
    )
  )
  
  
  invisible(
    list(
      all_significant = z
    )
  )
}
# ===== 5.4 EJECUTAR PREGUNTA 2 =====
if (run_question(2L)) {
message("")
message("============================================================")
message("METATEST — FIGURE 2")
message("TX pathways separated by analysis / phenotype-side / library")
message("============================================================")


fig2_results <- list()


for (family in FAMILIES) {
  
  message("\n[FIG2] ", family)
  
  
  inp2 <- load_fig2_inputs(
    bio_program_run = BIO_RUN,
    family = family
  )
  
  dat <- inp2$pathways
  side_map <- inp2$side_map
  
  # ---------------------------------------------------------------------------
  # Overview
  # ---------------------------------------------------------------------------
  
  make_fig2_overview(
    dat = dat,
    side_map = side_map,
    family = family,
    outdir = OUT_ROOT
  )
  
  
  # ---------------------------------------------------------------------------
  # Figuras separadas
  # ---------------------------------------------------------------------------
  
  factor_sides <- unique(
    side_map$factor_side
  )
  
  libraries <- sort(
    unique(dat$library)
  )
  
  blocks <- expand.grid(
    factor_side = factor_sides,
    library = libraries,
    stringsAsFactors = FALSE
  )
  
  family_results <- list()
  
  
  for (i in seq_len(nrow(blocks))) {
    
    factor_side_i <- blocks$factor_side[i]
    library_i <- blocks$library[i]
    
    
    message(
      "  ",
      factor_side_i,
      " | ",
      library_i
    )
    
    
    key <- paste(
      factor_side_i,
      library_i,
      sep = "__"
    )
    
    
    family_results[[key]] <- make_fig2_block(
      dat = dat,
      family = family,
      factor_side = factor_side_i,
      library = library_i,
      outdir = OUT_ROOT
    )
  }
  
  
  fig2_results[[family]] <- family_results
}


message("")
message("FIGURE 2 COMPLETE")
message(
  "Output: ",
  file.path(
    OUT_ROOT,
    "figures",
    "FIG2_PATHWAYS_TX"
  )
)
message("")
}
# ===== PREGUNTA 3. ¿QUÉ PATHWAYS SON REDUNDANTES Y FORMAN PROGRAMAS TX? =====

# ===== 6.1 CARGAR PROGRAMAS TX CONGELADOS POR 08 =====

load_fig3_inputs <- function(bio_program_run, family) {
  
  path <- file.path(
    bio_program_run,
    family,
    "07_TX_PROGRAM_PATHWAY_MEMBERSHIP.csv"
  )
  
  x <- read_csv_checked(
    path,
    required = c(
      "program_id",
      "membership_type",
      "library",
      "pathway",
      "factor",
      "factor_role",
      "phenotypic_side",
      "q_both_sides_within_library_factor",
      "n_features_in_set",
      "representative_library",
      "representative_pathway"
    )
  )
  
  
  # Solo membership principal.
  # SECONDARY/bridge se conservará para una figura específica posterior.
  x <- x[
    x$membership_type == "PRIMARY" &
      !is.na(x$q_both_sides_within_library_factor) &
      x$q_both_sides_within_library_factor <= 0.10,
    ,
    drop = FALSE
  ]
  
  
  assert(
    nrow(x) > 0L,
    paste(
      "No primary TX program memberships found for",
      family
    )
  )
  
  
  x$family <- family
  
  x$program_short <- short_program_id(
    x$program_id
  )
  
  x$pathway_label <- clean_pathway_label(
    x$pathway
  )
  
  x$representative_label <- clean_pathway_label(
    x$representative_pathway
  )
  
  x$is_representative <-
    x$library == x$representative_library &
    x$pathway == x$representative_pathway
  
  x$factor_side <- paste0(
    x$factor,
    "_",
    x$phenotypic_side
  )
  
  
  x
}

# ===== 6.2 RESUMEN DE PROGRAMAS TX =====

summarize_fig3_programs <- function(dat) {
  
  ids <- unique(dat$program_id)
  
  rows <- lapply(
    ids,
    function(pid) {
      
      z <- dat[
        dat$program_id == pid,
        ,
        drop = FALSE
      ]
      
      data.frame(
        family = z$family[1L],
        factor = z$factor[1L],
        factor_role = z$factor_role[1L],
        phenotypic_side = z$phenotypic_side[1L],
        factor_side = z$factor_side[1L],
        
        program_id = pid,
        program_short = z$program_short[1L],
        
        n_primary_pathways =
          length(unique(z$pathway)),
        
        n_libraries =
          length(unique(z$library)),
        
        representative_library =
          z$representative_library[1L],
        
        representative_pathway =
          z$representative_pathway[1L],
        
        representative_label =
          z$representative_label[1L],
        
        best_q =
          min(
            z$q_both_sides_within_library_factor,
            na.rm = TRUE
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  do.call(
    rbind,
    rows
  )
}
# ===== 6.3 FIGURA 3A — PATHWAYS AGRUPADOS EN PROGRAMAS TX =====

make_fig3_program_map <- function(
    dat,
    family,
    factor_side,
    outdir
) {
  
  z <- dat[
    dat$factor_side == factor_side,
    ,
    drop = FALSE
  ]
  
  if (!nrow(z)) {
    return(invisible(NULL))
  }
  
  
  z <- z[
    order(
      z$program_short,
      z$library,
      z$q_both_sides_within_library_factor,
      z$pathway
    ),
    ,
    drop = FALSE
  ]
  
  
  z$row_id <- paste(
    z$program_short,
    z$library,
    z$pathway,
    sep = "|||"
  )
  
  z$row_id <- factor(
    z$row_id,
    levels = rev(unique(z$row_id))
  )
  
  
  label_dictionary <- setNames(
    z$pathway_label,
    as.character(z$row_id)
  )
  
  
  z$pathway_role <- ifelse(
    z$is_representative,
    "Representative medoid",
    "Program member"
  )
  
  z$pathway_role <- factor(
    z$pathway_role,
    levels = c(
      "Program member",
      "Representative medoid"
    )
  )
  
  
  p <- ggplot2::ggplot(
    z,
    ggplot2::aes(
      x = program_short,
      y = row_id
    )
  ) +
    
    ggplot2::geom_point(
      ggplot2::aes(
        fill = library,
        shape = pathway_role
      ),
      size = 3.2,
      color = "black",
      stroke = 0.5
    ) +
    
    ggplot2::scale_shape_manual(
      values = c(
        "Program member" = 21,
        "Representative medoid" = 23
      )
    ) +
    
    ggplot2::scale_y_discrete(
      labels = function(x) {
        unname(
          label_dictionary[x]
        )
      }
    ) +
    
    ggplot2::labs(
      title = paste0(
        family,
        " — ",
        factor_side,
        " — programas transcriptómicos"
      ),
      
      subtitle = paste0(
        "Pathways significativos agrupados por redundancia ",
        "de genes TX observados."
      ),
      
      x = "Programa TX",
      y = NULL,
      
      fill = "Librería",
      shape = "Rol",
      
      caption = paste0(
        "Representative medoid = pathway estructuralmente más central ",
        "del programa. No implica mayor importancia biológica ni ",
        "mayor significación."
      )
    ) +
    
    fig_theme(10) +
    
    ggplot2::theme(
      axis.text.y =
        ggplot2::element_text(
          size = 7.5
        )
    )
  
  
  fig_height <- max(
    5,
    min(
      20,
      3 + 0.28 * nrow(z)
    )
  )
  
  
  dir_fig <- file.path(
    outdir,
    "figures",
    "FIG3_TX_PROGRAMS",
    "01_WITHIN_ANALYSIS",
    family,
    safe_filename(factor_side)
  )
  
  dir_dat <- file.path(
    outdir,
    "plot_data",
    "FIG3_TX_PROGRAMS",
    "01_WITHIN_ANALYSIS",
    family,
    safe_filename(factor_side)
  )
  
  dir.create(
    dir_fig,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    dir_dat,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  ggplot2::ggsave(
    filename = file.path(
      dir_fig,
      paste0(
        "FIG3A_",
        family,
        "_",
        safe_filename(factor_side),
        "_PROGRAM_MEMBERSHIP.png"
      )
    ),
    plot = p,
    width = 11,
    height = fig_height,
    dpi = FIG_DPI,
    limitsize = FALSE
  )
  
  
  write_csv(
    z,
    file.path(
      dir_dat,
      paste0(
        "FIG3A_",
        family,
        "_",
        safe_filename(factor_side),
        "_PROGRAM_MEMBERSHIP.csv"
      )
    )
  )
  
  
  invisible(z)
}
# ===== 6.4 EJECUTAR PREGUNTA 3A — PROGRAMAS DENTRO DE CADA ANÁLISIS =====
if (run_question(3L)) {
message("")
message("============================================================")
message("METATEST — QUESTION 3A")
message("Significant pathways -> frozen TX programs")
message("============================================================")


fig3_results <- list()
fig3_program_summaries <- list()


for (family in FAMILIES) {
  
  message("\n[FIG3] ", family)
  
  
  dat3 <- load_fig3_inputs(
    bio_program_run = BIO_RUN,
    family = family
  )
  
  
  program_summary <- summarize_fig3_programs(
    dat3
  )
  
  
  fig3_program_summaries[[family]] <-
    program_summary
  
  
  write_csv(
    program_summary,
    file.path(
      OUT_ROOT,
      "plot_data",
      "FIG3_TX_PROGRAMS",
      "01_WITHIN_ANALYSIS",
      family,
      paste0(
        "FIG3_",
        family,
        "_PROGRAM_SUMMARY.csv"
      )
    )
  )
  
  
  sides <- unique(
    dat3$factor_side
  )
  
  
  family_res <- list()
  
  
  for (side_i in sides) {
    
    message(
      "  ",
      side_i
    )
    
    family_res[[side_i]] <-
      make_fig3_program_map(
        dat = dat3,
        family = family,
        factor_side = side_i,
        outdir = OUT_ROOT
      )
  }
  
  
  fig3_results[[family]] <-
    family_res
}


all_program_summary <- do.call(
  rbind,
  fig3_program_summaries
)


write_csv(
  all_program_summary,
  file.path(
    OUT_ROOT,
    "plot_data",
    "FIG3_TX_PROGRAMS",
    "00_INDEX",
    "FIG3_ALL_ANALYSES_PROGRAM_SUMMARY.csv"
  )
)

message("")
message("QUESTION 3A COMPLETE")
}
# ===== 6.5 CARGAR ALINEAMIENTO ENTRE ANÁLISIS CONGELADO POR 08 =====

load_fig3_cross_alignment <- function(bio_program_run) {
  
  path <- file.path(
    bio_program_run,
    "cross_family",
    "01_PROGRAM_PAIRWISE_ALIGNMENT_PHENOTYPE_BASED.csv"
  )
  
  x <- read_csv_checked(
    path,
    required = c(
      "family_A",
      "family_B",
      "factor_role",
      "phenotypic_side",
      "program_A",
      "program_B",
      "representative_A",
      "representative_B",
      "n_genes_A",
      "n_genes_B",
      "n_shared_observed_genes",
      "jaccard_observed_genes",
      "overlap_coefficient_observed_genes",
      "n_shared_pathway_names",
      "jaccard_pathway_names",
      "best_for_A",
      "best_for_B",
      "reciprocal_best_match"
    )
  )
  
  as_flag <- function(v) {
    
    if (is.logical(v)) {
      return(v)
    }
    
    toupper(
      trimws(as.character(v))
    ) %in% c(
      "TRUE",
      "T",
      "1",
      "YES",
      "Y"
    )
  }
  
  x$best_for_A <- as_flag(x$best_for_A)
  x$best_for_B <- as_flag(x$best_for_B)
  x$reciprocal_best_match <- as_flag(
    x$reciprocal_best_match
  )
  
  x$reciprocal_positive_alignment <-
    x$reciprocal_best_match &
    x$n_shared_observed_genes > 0
  
  x$program_A_short <- short_program_id(
    x$program_A
  )
  
  x$program_B_short <- short_program_id(
    x$program_B
  )
  
  x$family_pair <- paste0(
    x$family_A,
    "_vs_",
    x$family_B
  )
  
  x
}
# ===== 6.6 FIGURA 3B — ALINEAMIENTO PAIRWISE DE PROGRAMAS TX =====


make_fig3b_alignment_matrix <- function(
    dat,
    phenotypic_side,
    family_A,
    family_B,
    outdir
) {
  
  z <- dat[
    dat$phenotypic_side == phenotypic_side &
      dat$family_A == family_A &
      dat$family_B == family_B,
    ,
    drop = FALSE
  ]
  
  
  if (!nrow(z)) {
    return(invisible(NULL))
  }
  
  
  z$program_A_short <- factor(
    z$program_A_short,
    levels = sort(
      unique(z$program_A_short)
    )
  )
  
  z$program_B_short <- factor(
    z$program_B_short,
    levels = sort(
      unique(z$program_B_short)
    )
  )
  
  
  z$cell_label <- sprintf(
    "%.2f",
    z$jaccard_observed_genes
  )
  
  
  z$cell_label[
    z$reciprocal_positive_alignment
  ] <- paste0(
    z$cell_label[
      z$reciprocal_positive_alignment
    ],
    " *"
  )
  
  
  p <- ggplot2::ggplot(
    z,
    ggplot2::aes(
      x = program_B_short,
      y = program_A_short,
      fill = jaccard_observed_genes
    )
  ) +
    
    ggplot2::geom_tile(
      color = "white",
      linewidth = 0.5
    ) +
    
    ggplot2::geom_text(
      ggplot2::aes(
        label = cell_label
      ),
      size = 3.5
    ) +
    
    ggplot2::labs(
      title = paste0(
        phenotype_label(phenotypic_side),
        " — ",
        family_A,
        " vs ",
        family_B
      ),
      
      subtitle = paste0(
        "Similitud entre programas TX basada en genes observados."
      ),
      
      x = paste0(
        family_B,
        " — programa TX"
      ),
      
      y = paste0(
        family_A,
        " — programa TX"
      ),
      
      fill = "Jaccard\ngenes",
      
      caption = paste0(
        "* = reciprocal best match con al menos un gen observado compartido. ",
        "El Jaccard se mantiene como magnitud continua; no se aplica ",
        "un umbral adicional para declarar similitud."
      )
    ) +
    
    fig_theme(10) +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
  
  
  dir_fig <- file.path(
    outdir,
    "figures",
    "FIG3_TX_PROGRAMS",
    "02_CROSS_ANALYSIS_SAME_PHENOTYPE",
    safe_filename(phenotypic_side)
  )
  
  dir_dat <- file.path(
    outdir,
    "plot_data",
    "FIG3_TX_PROGRAMS",
    "02_CROSS_ANALYSIS_SAME_PHENOTYPE",
    safe_filename(phenotypic_side)
  )
  
  
  dir.create(
    dir_fig,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    dir_dat,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  file_stub <- paste0(
    "FIG3B_",
    safe_filename(phenotypic_side),
    "_",
    family_A,
    "_vs_",
    family_B
  )
  
  
  ggplot2::ggsave(
    filename = file.path(
      dir_fig,
      paste0(
        file_stub,
        ".png"
      )
    ),
    plot = p,
    width = max(
      6,
      2.5 +
        0.8 *
        length(
          unique(z$program_B_short)
        )
    ),
    height = max(
      5,
      2.5 +
        0.7 *
        length(
          unique(z$program_A_short)
        )
    ),
    dpi = FIG_DPI,
    limitsize = FALSE
  )
  
  
  write_csv(
    z,
    file.path(
      dir_dat,
      paste0(
        file_stub,
        ".csv"
      )
    )
  )
  
  
  invisible(z)
}


# ===== 6.7 ESTADO CROSS-ANALYSIS DE CADA PROGRAMA =====


build_fig3_program_cross_status <- function(
    program_summary,
    alignment
) {
  
  programs <- unique(
    program_summary[, c(
      "family",
      "factor_role",
      "phenotypic_side",
      "program_id",
      "program_short",
      "representative_pathway",
      "representative_label"
    )]
  )
  
  
  aligned <- alignment[
    alignment$reciprocal_positive_alignment,
    ,
    drop = FALSE
  ]
  
  
  rows <- lapply(
    seq_len(nrow(programs)),
    function(i) {
      
      p <- programs[i, , drop = FALSE]
      
      pid <- p$program_id[1L]
      fam <- p$family[1L]
      
      
      as_A <- aligned[
        aligned$program_A == pid,
        ,
        drop = FALSE
      ]
      
      as_B <- aligned[
        aligned$program_B == pid,
        ,
        drop = FALSE
      ]
      
      
      matched_families <- c(
        as.character(as_A$family_B),
        as.character(as_B$family_A)
      )
      
      matched_programs <- c(
        as.character(as_A$program_B),
        as.character(as_B$program_A)
      )
      
      
      matched_families <- sort(
        unique(
          matched_families[
            nzchar(matched_families)
          ]
        )
      )
      
      matched_programs <- sort(
        unique(
          matched_programs[
            nzchar(matched_programs)
          ]
        )
      )
      
      
      n_other_analyses <- length(
        matched_families
      )
      
      
      status <- if (
        n_other_analyses == 0L
      ) {
        
        "Analysis-specific"
        
      } else if (
        n_other_analyses == 1L
      ) {
        
        "Aligned with one other analysis"
        
      } else {
        
        "Aligned with both other analyses"
      }
      
      
      jaccards <- c(
        as.numeric(
          as_A$jaccard_observed_genes
        ),
        as.numeric(
          as_B$jaccard_observed_genes
        )
      )
      
      jaccards <- jaccards[
        is.finite(jaccards)
      ]
      
      
      data.frame(
        family = fam,
        factor_role =
          p$factor_role[1L],
        
        phenotypic_side =
          p$phenotypic_side[1L],
        
        program_id = pid,
        program_short =
          p$program_short[1L],
        
        representative_pathway =
          p$representative_pathway[1L],
        
        representative_label =
          p$representative_label[1L],
        
        n_other_analyses_aligned =
          n_other_analyses,
        
        matched_families =
          if (length(matched_families)) {
            paste(
              matched_families,
              collapse = ";"
            )
          } else {
            NA_character_
          },
        
        matched_programs =
          if (length(matched_programs)) {
            paste(
              matched_programs,
              collapse = ";"
            )
          } else {
            NA_character_
          },
        
        max_reciprocal_jaccard =
          if (length(jaccards)) {
            max(jaccards)
          } else {
            NA_real_
          },
        
        cross_analysis_status =
          status,
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  
  do.call(
    rbind,
    rows
  )
}
# ===== 6.8 FIGURA 3C — ESPECIFICIDAD / COMPARTICIÓN POR FENOTIPO =====


make_fig3c_cross_status <- function(
    status_dat,
    outdir
) {
  
  summary_dat <- stats::aggregate(
    program_id ~
      family +
      phenotypic_side +
      cross_analysis_status,
    data = status_dat,
    FUN = length
  )
  
  
  names(summary_dat)[
    names(summary_dat) == "program_id"
  ] <- "n_programs"
  
  
  summary_dat$phenotype_label <-
    phenotype_label(
      summary_dat$phenotypic_side
    )
  
  
  p <- ggplot2::ggplot(
    summary_dat,
    ggplot2::aes(
      x = family,
      y = n_programs,
      fill = cross_analysis_status
    )
  ) +
    
    ggplot2::geom_col() +
    
    ggplot2::facet_wrap(
      ~ phenotype_label,
      scales = "free_y"
    ) +
    
    ggplot2::labs(
      title =
        "Programas TX específicos y alineados entre P1/P2/P3",
      
      subtitle = paste0(
        "Separados por lado fenotípico. ",
        "La correspondencia requiere reciprocal best match ",
        "y al menos un gen observado compartido."
      ),
      
      x = "Análisis",
      y = "Número de programas TX",
      
      fill = "Estado cross-analysis",
      
      caption = paste0(
        "Analysis-specific = sin correspondencia recíproca positiva ",
        "con otro análisis. Esto no implica ausencia de biología relacionada; ",
        "solo ausencia de una correspondencia de programa bajo la regla congelada."
      )
    ) +
    
    fig_theme(10)
  
  
  dir_fig <- file.path(
    outdir,
    "figures",
    "FIG3_TX_PROGRAMS",
    "02_CROSS_ANALYSIS_SAME_PHENOTYPE"
    )
  
  dir_dat <- file.path(
    outdir,
    "plot_data",
    "FIG3_TX_PROGRAMS",
    "02_CROSS_ANALYSIS_SAME_PHENOTYPE"
  )
  
  dir.create(
    dir_fig,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    dir_dat,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  ggplot2::ggsave(
    filename = file.path(
      dir_fig,
      "FIG3C_PROGRAM_CROSS_ANALYSIS_STATUS.png"
    ),
    plot = p,
    width = 10,
    height = 7,
    dpi = FIG_DPI,
    limitsize = FALSE
  )
  
  
  write_csv(
    summary_dat,
    file.path(
      dir_dat,
      "FIG3C_PROGRAM_CROSS_ANALYSIS_STATUS_SUMMARY.csv"
    )
  )
  
  
  invisible(summary_dat)
}

# ===== 6.9 EJECUTAR PREGUNTA 3B — COMPARACIÓN P1/P2/P3 =====

if (run_question(3L)) {
message("")
message("============================================================")
message("METATEST — QUESTION 3B")
message("Cross-analysis TX program alignment")
message("============================================================")


cross_alignment <- load_fig3_cross_alignment(
  BIO_RUN
)


write_csv(
  cross_alignment,
  file.path(
    OUT_ROOT,
    "plot_data",
    "FIG3_TX_PROGRAMS",
    "02_CROSS_ANALYSIS_SAME_PHENOTYPE",
    "FIG3B_ALL_PAIRWISE_PROGRAM_ALIGNMENTS.csv"
  )
)


phenotypic_sides <- unique(
  cross_alignment$phenotypic_side
)


family_pairs <- unique(
  cross_alignment[, c(
    "family_A",
    "family_B"
  )]
)


for (side_i in phenotypic_sides) {
  
  message(
    "\n[FIG3B] ",
    side_i
  )
  
  
  for (
    i in seq_len(
      nrow(family_pairs)
    )
  ) {
    
    make_fig3b_alignment_matrix(
      dat = cross_alignment,
      phenotypic_side = side_i,
      family_A =
        family_pairs$family_A[i],
      family_B =
        family_pairs$family_B[i],
      outdir = OUT_ROOT
    )
  }
}


program_cross_status <-
  build_fig3_program_cross_status(
    program_summary =
      all_program_summary,
    alignment =
      cross_alignment
  )


write_csv(
  program_cross_status,
  file.path(
    OUT_ROOT,
    "plot_data",
    "FIG3_TX_PROGRAMS",
    "02_CROSS_ANALYSIS_SAME_PHENOTYPE",
    "FIG3_PROGRAM_CROSS_ANALYSIS_STATUS.csv"
  )
)


make_fig3c_cross_status(
  status_dat =
    program_cross_status,
  outdir =
    OUT_ROOT
)


message("")
message("QUESTION 3B COMPLETE")

}

# ===== 6.10 PREGUNTA 3C — ESPECIFICIDAD POR LADO FENOTÍPICO =====


load_fig3_program_genes <- function(
    bio_program_run,
    families
) {
  
  rows <- lapply(
    families,
    function(family) {
      
      path <- file.path(
        bio_program_run,
        family,
        "08_TX_PROGRAM_GENE_MEMBERSHIP.csv"
      )
      
      x <- read_csv_checked(
        path,
        required = c(
          "program_id",
          "factor",
          "factor_role",
          "phenotypic_side",
          "feature_original",
          "n_pathways",
          "n_primary_pathways"
        )
      )
      
      # Solo genes sustentados por pathways PRIMARY.
      # Las afiliaciones secondary/bridge no definen esta comparación.
      x <- x[
        !is.na(x$feature_original) &
          nzchar(as.character(x$feature_original)) &
          as.numeric(x$n_primary_pathways) > 0,
        ,
        drop = FALSE
      ]
      
      x$family <- family
      
      x
    }
  )
  
  do.call(
    rbind,
    rows
  )
}

# ===== 6.11 SIMILITUD ENTRE TODOS LOS PROGRAMAS TX =====


build_all_program_pairwise <- function(
    program_summary,
    program_genes
) {
  
  meta <- unique(
    program_summary[, c(
      "family",
      "factor",
      "factor_role",
      "phenotypic_side",
      "program_id",
      "program_short",
      "representative_pathway",
      "representative_label"
    )]
  )
  
  gene_sets <- split(
    as.character(program_genes$feature_original),
    as.character(program_genes$program_id)
  )
  
  gene_sets <- lapply(
    gene_sets,
    unique
  )
  
  
  rows <- list()
  rr <- 1L
  
  
  if (nrow(meta) < 2L) {
    return(data.frame())
  }
  
  
  for (i in seq_len(nrow(meta) - 1L)) {
    
    for (j in (i + 1L):nrow(meta)) {
      
      A <- meta[i, , drop = FALSE]
      B <- meta[j, , drop = FALSE]
      
      ga <- gene_sets[[A$program_id]]
      gb <- gene_sets[[B$program_id]]
      
      if (is.null(ga)) ga <- character()
      if (is.null(gb)) gb <- character()
      
      shared <- intersect(
        ga,
        gb
      )
      
      union_genes <- union(
        ga,
        gb
      )
      
      
      jaccard <- if (
        length(union_genes) > 0L
      ) {
        
        length(shared) /
          length(union_genes)
        
      } else {
        
        NA_real_
      }
      
      
      overlap <- if (
        min(
          length(ga),
          length(gb)
        ) > 0L
      ) {
        
        length(shared) /
          min(
            length(ga),
            length(gb)
          )
        
      } else {
        
        NA_real_
      }
      
      
      same_family <-
        A$family == B$family
      
      same_side <-
        A$phenotypic_side ==
        B$phenotypic_side
      
      
      comparison_scope <- if (
        same_family && same_side
      ) {
        
        "same_analysis_same_phenotype"
        
      } else if (
        same_family && !same_side
      ) {
        
        "same_analysis_cross_phenotype"
        
      } else if (
        !same_family && same_side
      ) {
        
        "cross_analysis_same_phenotype"
        
      } else {
        
        "cross_analysis_cross_phenotype"
      }
      
      
      rows[[rr]] <- data.frame(
        
        family_A =
          A$family,
        
        family_B =
          B$family,
        
        phenotypic_side_A =
          A$phenotypic_side,
        
        phenotypic_side_B =
          B$phenotypic_side,
        
        program_A =
          A$program_id,
        
        program_B =
          B$program_id,
        
        program_A_short =
          A$program_short,
        
        program_B_short =
          B$program_short,
        
        representative_A =
          A$representative_pathway,
        
        representative_B =
          B$representative_pathway,
        
        n_genes_A =
          length(ga),
        
        n_genes_B =
          length(gb),
        
        n_shared_observed_genes =
          length(shared),
        
        jaccard_observed_genes =
          jaccard,
        
        overlap_coefficient_observed_genes =
          overlap,
        
        comparison_scope =
          comparison_scope,
        
        stringsAsFactors = FALSE
      )
      
      
      rr <- rr + 1L
    }
  }
  
  
  out <- do.call(
    rbind,
    rows
  )
  
  out$pair_id <- seq_len(
    nrow(out)
  )
  
  out
}
# ===== 6.12 RECIPROCAL-BEST ENTRE LADOS FENOTÍPICOS =====


mark_cross_phenotype_reciprocal <- function(
    pairwise
) {
  
  pairwise$reciprocal_best_cross_phenotype <-
    FALSE
  
  
  idx <- which(
    pairwise$comparison_scope ==
      "same_analysis_cross_phenotype"
  )
  
  
  if (!length(idx)) {
    return(pairwise)
  }
  
  
  z <- pairwise[
    idx,
    ,
    drop = FALSE
  ]
  
  
  # Representación dirigida para poder preguntar:
  # para cada programa, ¿cuál es su mejor correspondencia
  # entre los OTROS lados fenotípicos del mismo análisis?
  long <- rbind(
    
    data.frame(
      pair_id = z$pair_id,
      source = z$program_A,
      target = z$program_B,
      jaccard =
        z$jaccard_observed_genes,
      stringsAsFactors = FALSE
    ),
    
    data.frame(
      pair_id = z$pair_id,
      source = z$program_B,
      target = z$program_A,
      jaccard =
        z$jaccard_observed_genes,
      stringsAsFactors = FALSE
    )
  )
  
  
  best_by_source <- tapply(
    long$jaccard,
    long$source,
    function(v) {
      max(
        v,
        na.rm = TRUE
      )
    }
  )
  
  
  long$is_best <- mapply(
    function(source, value) {
      
      best <- unname(
        best_by_source[source]
      )
      
      is.finite(best) &&
        is.finite(value) &&
        abs(value - best) < 1e-12
      
    },
    long$source,
    long$jaccard
  )
  
  
  best_count <- stats::aggregate(
    is_best ~ pair_id,
    data = long,
    FUN = sum
  )
  
  
  reciprocal_ids <-
    best_count$pair_id[
      best_count$is_best == 2L
    ]
  
  
  reciprocal_ids <-
    reciprocal_ids[
      pairwise$n_shared_observed_genes[
        match(
          reciprocal_ids,
          pairwise$pair_id
        )
      ] > 0L
    ]
  
  
  pairwise$reciprocal_best_cross_phenotype[
    pairwise$pair_id %in%
      reciprocal_ids
  ] <- TRUE
  
  
  pairwise
}

# ===== 6.13 ESTADO DE ESPECIFICIDAD FENOTÍPICA =====


build_program_phenotype_status <- function(
    program_summary,
    pairwise
) {
  
  programs <- unique(
    program_summary[, c(
      "family",
      "program_id",
      "program_short",
      "phenotypic_side"
    )]
  )
  
  
  aligned <- pairwise[
    pairwise$comparison_scope ==
      "same_analysis_cross_phenotype" &
      pairwise$reciprocal_best_cross_phenotype,
    ,
    drop = FALSE
  ]
  
  
  rows <- lapply(
    seq_len(nrow(programs)),
    function(i) {
      
      p <- programs[
        i,
        ,
        drop = FALSE
      ]
      
      pid <- p$program_id[1L]
      
      
      as_A <- aligned[
        aligned$program_A == pid,
        ,
        drop = FALSE
      ]
      
      as_B <- aligned[
        aligned$program_B == pid,
        ,
        drop = FALSE
      ]
      
      
      matched_sides <- c(
        as.character(
          as_A$phenotypic_side_B
        ),
        as.character(
          as_B$phenotypic_side_A
        )
      )
      
      
      matched_programs <- c(
        as.character(
          as_A$program_B
        ),
        as.character(
          as_B$program_A
        )
      )
      
      
      matched_sides <- sort(
        unique(
          matched_sides[
            !is.na(matched_sides) &
              nzchar(matched_sides)
          ]
        )
      )
      
      
      matched_programs <- sort(
        unique(
          matched_programs[
            !is.na(matched_programs) &
              nzchar(matched_programs)
          ]
        )
      )
      
      
      jaccards <- c(
        as.numeric(
          as_A$jaccard_observed_genes
        ),
        as.numeric(
          as_B$jaccard_observed_genes
        )
      )
      
      jaccards <- jaccards[
        is.finite(jaccards)
      ]
      
      
      n_other_sides <-
        length(matched_sides)
      
      
      phenotype_status <- if (
        n_other_sides == 0L
      ) {
        
        "Phenotype-side-specific"
        
      } else if (
        n_other_sides == 1L
      ) {
        
        "Aligned with one other phenotype side"
        
      } else {
        
        "Aligned with multiple phenotype sides"
      }
      
      
      data.frame(
        
        family =
          p$family[1L],
        
        program_id =
          pid,
        
        phenotypic_side =
          p$phenotypic_side[1L],
        
        n_other_phenotype_sides_aligned =
          n_other_sides,
        
        matched_phenotypic_sides =
          if (length(matched_sides)) {
            paste(
              matched_sides,
              collapse = ";"
            )
          } else {
            NA_character_
          },
        
        matched_cross_side_programs =
          if (length(matched_programs)) {
            paste(
              matched_programs,
              collapse = ";"
            )
          } else {
            NA_character_
          },
        
        max_cross_side_jaccard =
          if (length(jaccards)) {
            max(jaccards)
          } else {
            NA_real_
          },
        
        cross_phenotype_status =
          phenotype_status,
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  
  do.call(
    rbind,
    rows
  )
}

# ===== 6.14 TRÍADAS ESTRICTAS P1 ↔ P2 ↔ P3 =====


build_strict_three_analysis_triads <- function(
    alignment
) {
  
  z <- alignment[
    alignment$reciprocal_positive_alignment,
    ,
    drop = FALSE
  ]
  
  
  normalize_pair <- function(
    dat,
    left_family,
    right_family
  ) {
    
    direct <- dat[
      dat$family_A == left_family &
        dat$family_B == right_family,
      ,
      drop = FALSE
    ]
    
    if (nrow(direct)) {
      
      return(
        data.frame(
          phenotypic_side =
            direct$phenotypic_side,
          program_left =
            direct$program_A,
          program_right =
            direct$program_B,
          representative_left =
            direct$representative_A,
          representative_right =
            direct$representative_B,
          jaccard =
            direct$jaccard_observed_genes,
          stringsAsFactors = FALSE
        )
      )
    }
    
    
    reverse <- dat[
      dat$family_A == right_family &
        dat$family_B == left_family,
      ,
      drop = FALSE
    ]
    
    
    data.frame(
      phenotypic_side =
        reverse$phenotypic_side,
      program_left =
        reverse$program_B,
      program_right =
        reverse$program_A,
      representative_left =
        reverse$representative_B,
      representative_right =
        reverse$representative_A,
      jaccard =
        reverse$jaccard_observed_genes,
      stringsAsFactors = FALSE
    )
  }
  
  
  e12 <- normalize_pair(
    z,
    "P1",
    "P2"
  )
  
  e13 <- normalize_pair(
    z,
    "P1",
    "P3"
  )
  
  e23 <- normalize_pair(
    z,
    "P2",
    "P3"
  )
  
  
  rows <- list()
  rr <- 1L
  
  
  for (side_i in unique(
    c(
      e12$phenotypic_side,
      e13$phenotypic_side,
      e23$phenotypic_side
    )
  )) {
    
    a12 <- e12[
      e12$phenotypic_side == side_i,
      ,
      drop = FALSE
    ]
    
    a13 <- e13[
      e13$phenotypic_side == side_i,
      ,
      drop = FALSE
    ]
    
    a23 <- e23[
      e23$phenotypic_side == side_i,
      ,
      drop = FALSE
    ]
    
    
    if (
      !nrow(a12) ||
      !nrow(a13) ||
      !nrow(a23)
    ) {
      next
    }
    
    
    for (i in seq_len(nrow(a12))) {
      
      p1 <- a12$program_left[i]
      p2 <- a12$program_right[i]
      
      
      cand13 <- a13[
        a13$program_left == p1,
        ,
        drop = FALSE
      ]
      
      
      if (!nrow(cand13)) {
        next
      }
      
      
      for (j in seq_len(nrow(cand13))) {
        
        p3 <- cand13$program_right[j]
        
        
        hit23 <- a23[
          a23$program_left == p2 &
            a23$program_right == p3,
          ,
          drop = FALSE
        ]
        
        
        if (!nrow(hit23)) {
          next
        }
        
        
        rows[[rr]] <- data.frame(
          
          phenotypic_side =
            side_i,
          
          program_P1 =
            p1,
          
          program_P2 =
            p2,
          
          program_P3 =
            p3,
          
          representative_P1 =
            a12$representative_left[i],
          
          representative_P2 =
            a12$representative_right[i],
          
          representative_P3 =
            cand13$representative_right[j],
          
          jaccard_P1_P2 =
            a12$jaccard[i],
          
          jaccard_P1_P3 =
            cand13$jaccard[j],
          
          jaccard_P2_P3 =
            hit23$jaccard[1L],
          
          strict_three_way_reciprocal =
            TRUE,
          
          stringsAsFactors = FALSE
        )
        
        
        rr <- rr + 1L
      }
    }
  }
  
  
  if (!length(rows)) {
    return(data.frame())
  }
  
  
  unique(
    do.call(
      rbind,
      rows
    )
  )
}

# ===== 6.15 MASTER PROGRAM AUDIT =====


build_master_program_audit <- function(
    program_cross_status,
    phenotype_status,
    strict_triads
) {
  
  master <- merge(
    program_cross_status,
    phenotype_status[, c(
      "family",
      "program_id",
      "n_other_phenotype_sides_aligned",
      "matched_phenotypic_sides",
      "matched_cross_side_programs",
      "max_cross_side_jaccard",
      "cross_phenotype_status"
    )],
    by = c(
      "family",
      "program_id"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  
  
  strict_programs <- character()
  
  if (nrow(strict_triads)) {
    
    strict_programs <- unique(
      c(
        strict_triads$program_P1,
        strict_triads$program_P2,
        strict_triads$program_P3
      )
    )
  }
  
  
  master$strict_P1_P2_P3_triad <-
    master$program_id %in%
    strict_programs
  
  
  master$strict_three_analysis_status <-
    ifelse(
      master$strict_P1_P2_P3_triad,
      "Strict P1-P2-P3 triad",
      ifelse(
        master$n_other_analyses_aligned >= 2L,
        "Aligned with both; not strict triad",
        "Not a three-analysis triad"
      )
    )
  
  
  master$combined_scope <- paste(
    master$cross_analysis_status,
    master$cross_phenotype_status,
    sep = " | "
  )
  
  
  master <- master[
    order(
      master$phenotypic_side,
      master$family,
      master$program_short
    ),
    ,
    drop = FALSE
  ]
  
  
  master
}

# ===== 6.16 FIGURA 3D — MAPA MAESTRO DE ESPECIFICIDAD =====


make_fig3_master_audit <- function(
    master,
    outdir
) {
  
  master$program_label <- paste(
    master$family,
    phenotype_label(
      master$phenotypic_side
    ),
    master$program_short,
    sep = " | "
  )
  
  
  plot_dat <- rbind(
    
    data.frame(
      program_label =
        master$program_label,
      dimension =
        "Entre análisis\n(mismo lado)",
      status =
        master$cross_analysis_status,
      stringsAsFactors = FALSE
    ),
    
    data.frame(
      program_label =
        master$program_label,
      dimension =
        "Entre lados\n(mismo análisis)",
      status =
        master$cross_phenotype_status,
      stringsAsFactors = FALSE
    ),
    
    data.frame(
      program_label =
        master$program_label,
      dimension =
        "Tríada estricta\nP1-P2-P3",
      status =
        master$strict_three_analysis_status,
      stringsAsFactors = FALSE
    )
  )
  
  
  program_order <- rev(
    unique(
      master$program_label
    )
  )
  
  
  plot_dat$program_label <- factor(
    plot_dat$program_label,
    levels = program_order
  )
  
  
  p <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(
      x = dimension,
      y = program_label
    )
  ) +
    
    ggplot2::geom_tile(
      fill = "white",
      color = "grey70",
      linewidth = 0.5
    ) +
    
    ggplot2::geom_text(
      ggplot2::aes(
        label = status
      ),
      size = 2.7
    ) +
    
    ggplot2::labs(
      title =
        "Pregunta 3 — mapa maestro de programas TX",
      
      subtitle = paste0(
        "Cada programa se audita por conservación entre análisis, ",
        "especificidad fenotípica y pertenencia a una tríada estricta."
      ),
      
      x = NULL,
      y = NULL,
      
      caption = paste0(
        "La especificidad es relativa a las correspondencias ",
        "reciprocal-best basadas en genes TX observados; ",
        "no implica exclusividad biológica absoluta."
      )
    ) +
    
    fig_theme(10) +
    
    ggplot2::theme(
      panel.grid =
        ggplot2::element_blank(),
      
      axis.text.y =
        ggplot2::element_text(
          size = 7
        )
    )
  
  
  dir_fig <- file.path(
    outdir,
    "figures",
    "FIG3_TX_PROGRAMS",
    "04_GLOBAL_AUDIT"
  )
  
  dir.create(
    dir_fig,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  ggplot2::ggsave(
    filename = file.path(
      dir_fig,
      "FIG3D_MASTER_PROGRAM_AUDIT.png"
    ),
    plot = p,
    width = 13,
    height = max(
      7,
      3 + 0.38 * nrow(master)
    ),
    dpi = FIG_DPI,
    limitsize = FALSE
  )
  
  
  invisible(p)
}
# ===== 6.17 EJECUTAR AUDITORÍA COMPLETA DE ESPECIFICIDAD =====

if (run_question(3L)) {
message("")
message("============================================================")
message("METATEST — QUESTION 3C")
message("Analysis specificity × phenotype-side specificity")
message("============================================================")


program_genes <- load_fig3_program_genes(
  BIO_RUN,
  FAMILIES
)


all_program_pairwise <-
  build_all_program_pairwise(
    program_summary =
      all_program_summary,
    program_genes =
      program_genes
  )


all_program_pairwise <-
  mark_cross_phenotype_reciprocal(
    all_program_pairwise
  )


phenotype_status <-
  build_program_phenotype_status(
    program_summary =
      all_program_summary,
    pairwise =
      all_program_pairwise
  )


strict_triads <-
  build_strict_three_analysis_triads(
    cross_alignment
  )


master_program_audit <-
  build_master_program_audit(
    program_cross_status =
      program_cross_status,
    phenotype_status =
      phenotype_status,
    strict_triads =
      strict_triads
  )

global_dat_dir <- file.path(
  OUT_ROOT,
  "plot_data",
  "FIG3_TX_PROGRAMS",
  "04_GLOBAL_AUDIT"
)

dir.create(
  global_dat_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


write_csv(
  all_program_pairwise,
  file.path(
    global_dat_dir,
    "ALL_PROGRAM_PAIRWISE_SIMILARITY.csv"
  )
)


write_csv(
  strict_triads,
  file.path(
    global_dat_dir,
    "STRICT_P1_P2_P3_TRIADS.csv"
  )
)


write_csv(
  master_program_audit,
  file.path(
    global_dat_dir,
    "FIG3_MASTER_PROGRAM_AUDIT.csv"
  )
)


make_fig3_master_audit(
  master =
    master_program_audit,
  outdir =
    OUT_ROOT
)
for (family in FAMILIES) {
  
  z <- all_program_pairwise[
    all_program_pairwise$comparison_scope ==
      "same_analysis_cross_phenotype" &
      all_program_pairwise$family_A ==
      family &
      all_program_pairwise$family_B ==
      family,
    ,
    drop = FALSE
  ]
  
  
  family_dir <- file.path(
    OUT_ROOT,
    "plot_data",
    "FIG3_TX_PROGRAMS",
    "03_CROSS_PHENOTYPE_WITHIN_ANALYSIS",
    family
  )
  
  
  dir.create(
    family_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  write_csv(
    z,
    file.path(
      family_dir,
      paste0(
        "FIG3_",
        family,
        "_CROSS_PHENOTYPE_PAIRWISE.csv"
      )
    )
  )
}

# ===== 6.18 ÍNDICE DE SALIDAS DE LA PREGUNTA 3 =====


fig3_output_index <- data.frame(
  
  folder = c(
    "01_WITHIN_ANALYSIS",
    "02_CROSS_ANALYSIS_SAME_PHENOTYPE",
    "03_CROSS_PHENOTYPE_WITHIN_ANALYSIS",
    "04_GLOBAL_AUDIT"
  ),
  
  question = c(
    "¿Qué pathways redundantes forman cada programa dentro de cada análisis y lado?",
    "¿Qué programas del mismo lado se conservan entre P1/P2/P3?",
    "¿Un programa asignado a un lado también se parece a programas de otros lados dentro del mismo análisis?",
    "¿Cuál es el estado conjunto análisis × fenotipo de cada programa?"
  ),
  
  interpretation = c(
    "Reducción pathways -> programas",
    "Especificidad respecto a la pregunta científica",
    "Especificidad respecto al lado fenotípico",
    "Resumen maestro y tríadas estrictas P1-P2-P3"
  ),
  
  stringsAsFactors = FALSE
)


index_dir <- file.path(
  OUT_ROOT,
  "plot_data",
  "FIG3_TX_PROGRAMS",
  "00_INDEX"
)

dir.create(
  index_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


write_csv(
  fig3_output_index,
  file.path(
    index_dir,
    "FIG3_OUTPUT_INDEX.csv"
  )
)


write_csv(
  master_program_audit,
  file.path(
    index_dir,
    "FIG3_MASTER_PROGRAM_AUDIT.csv"
  )
)


message("")
message("QUESTION 3 COMPLETE")

}
# ===== PREGUNTA 4. ¿QUÉ GENES SOSTIENEN LOS PROGRAMAS TX? =====


# ===== 7.1 CARGAR GENES, ANOTACIÓN Y PESOS CONGELADOS =====

load_fig4_inputs <- function(
    bio_program_run,
    family
) {
  
  family_dir <- file.path(
    bio_program_run,
    family
  )
  
  # Auditoría: permite recuperar exactamente las corridas usadas por 08.
  audit <- read_csv_checked(
    file.path(
      family_dir,
      "00_PROGRAM_DEFINITION_AUDIT.csv"
    ),
    required = c(
      "family",
      "enrichment_run",
      "observed_train_run",
      "crossomic_run"
    )
  )
  
  assert(
    nrow(audit) >= 1L,
    paste("Empty program audit for", family)
  )
  
  enrichment_run <- as.character(
    audit$enrichment_run[1L]
  )
  
  observed_train_run <- as.character(
    audit$observed_train_run[1L]
  )
  crossomic_run <- as.character(
    audit$crossomic_run[1L]
  )
  # ===== Gen -> programa =====
  
  program_genes <- read_csv_checked(
    file.path(
      family_dir,
      "08_TX_PROGRAM_GENE_MEMBERSHIP.csv"
    ),
    required = c(
      "program_id",
      "factor",
      "factor_role",
      "phenotypic_side",
      "feature_original",
      "n_pathways",
      "pathways",
      "n_primary_pathways",
      "is_driver_in_any_pathway",
      "n_driver_pathways",
      "gene_bridge_across_pathways"
    )
  )
  
  program_genes <- program_genes[
    !is.na(program_genes$feature_original) &
      nzchar(as.character(program_genes$feature_original)) &
      as.numeric(program_genes$n_primary_pathways) > 0,
    ,
    drop = FALSE
  ]
  
  program_genes$feature_original <-
    as.character(program_genes$feature_original)
  
  program_genes$family <- family
  
  program_genes$program_short <-
    short_program_id(
      program_genes$program_id
    )
  
  program_genes$factor_side <- paste0(
    program_genes$factor,
    "_",
    program_genes$phenotypic_side
  )
  
  # ===== Pathway -> programa =====
  
  program_pathways <- read_csv_checked(
    file.path(
      family_dir,
      "07_TX_PROGRAM_PATHWAY_MEMBERSHIP.csv"
    ),
    required = c(
      "program_id",
      "membership_type",
      "library",
      "pathway",
      "factor",
      "factor_role",
      "phenotypic_side"
    )
  )
  
  program_pathways <- program_pathways[
    program_pathways$membership_type == "PRIMARY",
    ,
    drop = FALSE
  ]
  
  program_pathways$family <- family
  
  program_pathways$program_short <-
    short_program_id(
      program_pathways$program_id
    )
  
  # ===== Pathway -> todos los genes observados =====
  
  observed <- read_csv_checked(
    file.path(
      family_dir,
      "03_PRIMARY_PATHWAY_GENE_OBSERVED_MEMBERSHIP_LONG.csv"
    ),
    required = c(
      "family",
      "library",
      "pathway",
      "factor",
      "factor_role",
      "phenotypic_side",
      "feature_original",
      "formal_pathway_membership",
      "observed_in_frozen_tx_universe",
      "is_side_driver"
    )
  )
  
  observed$feature_original <-
    as.character(observed$feature_original)
  
  # ===== Signo orientado correspondiente a cada lado fenotípico =====
  
  side_map <- read_csv_checked(
    file.path(
      family_dir,
      "01_FACTOR_PHENOTYPIC_SIDE_MAP.csv"
    ),
    required = c(
      "factor",
      "oriented_side",
      "phenotypic_side"
    )
  )
  
  side_key <- unique(
    side_map[, c(
      "factor",
      "phenotypic_side",
      "oriented_side"
    )]
  )
  
  # ===== Anotación feature_original -> GeneSymbol =====
  
  enrichment_object_path <- file.path(
    enrichment_run,
    "rds",
    "MOFA_NATIVE_ENRICHMENT_OBJECT.rds"
  )
  
  assert(
    file.exists(enrichment_object_path),
    paste(
      "Missing enrichment object:",
      enrichment_object_path
    )
  )
  
  enrichment_object <- readRDS(
    enrichment_object_path
  )
  
  tx_annotation <- as.data.frame(
    enrichment_object$tx_annotation,
    stringsAsFactors = FALSE
  )
  # ===== Loading latente TX orientado y escalado =====
  
  tx_latent_loadings <- read_csv_checked(
    file.path(
      crossomic_run,
      "tables",
      "01_ORIENTED_SCALED_LOADINGS_ALL_VIEWS.csv"
    ),
    required = c(
      "view_code",
      "feature_original",
      "factor",
      "loading_oriented_scaled"
    )
  )
  
  tx_latent_loadings <- tx_latent_loadings[
    tolower(
      as.character(
        tx_latent_loadings$view_code
      )
    ) == "tx",
    ,
    drop = FALSE
  ]
  
  tx_latent_loadings$feature_original <-
    as.character(
      tx_latent_loadings$feature_original
    )
  
  tx_latent_loadings$loading_oriented_scaled <-
    suppressWarnings(
      as.numeric(
        tx_latent_loadings$loading_oriented_scaled
      )
    )
  
  tx_latent_loadings$abs_loading_oriented_scaled <-
    abs(
      tx_latent_loadings$loading_oriented_scaled
    )
  
  tx_latent_loadings <- unique(
    tx_latent_loadings[
      ,
      c(
        "feature_original",
        "factor",
        "loading_oriented_scaled",
        "abs_loading_oriented_scaled"
      ),
      drop = FALSE
    ]
  )
  assert(
    all(
      c(
        "feature_original",
        "GeneSymbol"
      ) %in% names(tx_annotation)
    ),
    "tx_annotation must contain feature_original and GeneSymbol"
  )
  
  tx_annotation$feature_original <-
    as.character(tx_annotation$feature_original)
  
  tx_annotation <- unique(
    tx_annotation[, c(
      "feature_original",
      "GeneSymbol"
    )]
  )
  
  # ===== Pesos driver ya congelados en OBSERVED TRAIN =====
  #
  # NO se vuelven a calcular.
  # Son los mismos usados para construir el pathway score observado.
  
  observed_object_path <- file.path(
    observed_train_run,
    "rds",
    "OBSERVED_TRAIN_CORROBORATION_OBJECT.rds"
  )
  
  assert(
    file.exists(observed_object_path),
    paste(
      "Missing observed TRAIN object:",
      observed_object_path
    )
  )
  
  observed_object <- readRDS(
    observed_object_path
  )
  
  driver_weights <- as.data.frame(
    observed_object$pathway_scoring_spec,
    stringsAsFactors = FALSE
  )
  
  assert(
    all(
      c(
        "library",
        "pathway",
        "factor",
        "oriented_side",
        "feature_original",
        "anchor_weight",
        "normalized_anchor_weight"
      ) %in% names(driver_weights)
    ),
    paste(
      "pathway_scoring_spec does not contain the expected driver-weight columns"
    )
  )
  
  driver_weights$feature_original <-
    as.character(driver_weights$feature_original)
  
  driver_weights <- unique(
    driver_weights[, c(
      "library",
      "pathway",
      "factor",
      "oriented_side",
      "feature_original",
      "anchor_weight",
      "normalized_anchor_weight"
    )]
  )
  
  # ===== Construir pathway x gen dentro de cada programa =====
  
  pathway_gene <- merge(
    program_pathways[
      ,
      c(
        "family",
        "program_id",
        "program_short",
        "library",
        "pathway",
        "factor",
        "factor_role",
        "phenotypic_side"
      ),
      drop = FALSE
    ],
    observed[
      ,
      c(
        "family",
        "library",
        "pathway",
        "factor",
        "factor_role",
        "phenotypic_side",
        "feature_original",
        "is_side_driver"
      ),
      drop = FALSE
    ],
    by = c(
      "family",
      "library",
      "pathway",
      "factor",
      "factor_role",
      "phenotypic_side"
    ),
    all = FALSE,
    sort = FALSE
  )
  
  # Añadir oriented_side sin asumir que positivo/negativo
  # significa el mismo fenotipo entre modelos.
  
  pathway_gene <- merge(
    pathway_gene,
    side_key,
    by = c(
      "factor",
      "phenotypic_side"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  
  # Añadir pesos congelados únicamente donde corresponden.
  
  pathway_gene <- merge(
    pathway_gene,
    driver_weights,
    by = c(
      "library",
      "pathway",
      "factor",
      "oriented_side",
      "feature_original"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  
  # Añadir símbolo génico.
  
  pathway_gene <- merge(
    pathway_gene,
    tx_annotation,
    by = "feature_original",
    all.x = TRUE,
    sort = FALSE
  )
  
  program_genes <- merge(
    program_genes,
    tx_annotation,
    by = "feature_original",
    all.x = TRUE,
    sort = FALSE
  )
  
  pathway_gene$gene_label <- ifelse(
    !is.na(pathway_gene$GeneSymbol) &
      nzchar(pathway_gene$GeneSymbol),
    pathway_gene$GeneSymbol,
    pathway_gene$feature_original
  )
  
  program_genes$gene_label <- ifelse(
    !is.na(program_genes$GeneSymbol) &
      nzchar(program_genes$GeneSymbol),
    program_genes$GeneSymbol,
    program_genes$feature_original
  )
  
  pathway_gene$factor_side <- paste0(
    pathway_gene$factor,
    "_",
    pathway_gene$phenotypic_side
  )
  
  pathway_gene$pathway_label <-
    clean_pathway_label(
      pathway_gene$pathway
    )
  pathway_gene <- merge(
    pathway_gene,
    tx_latent_loadings,
    by = c(
      "feature_original",
      "factor"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  list(
    program_genes = program_genes,
    pathway_gene = pathway_gene,
    program_pathways = program_pathways,
    tx_annotation = tx_annotation,
    driver_weights = driver_weights,
    enrichment_run = enrichment_run,
    observed_train_run = observed_train_run
  )
}

# ===== 7.2 CONSTRUIR TABLA MAESTRA DE GENES POR PROGRAMA =====


build_fig4_gene_master <- function(
    inputs,
    program_audit = NULL
) {
  
  x <- inputs$program_genes
  
  
  # ---------------------------------------------------------------------------
  # Número de programas distintos en los que aparece cada gen,
  # dentro de toda la familia P1/P2/P3 correspondiente.
  # ---------------------------------------------------------------------------
  
  n_programs_gene <- tapply(
    x$program_id,
    x$feature_original,
    function(z) {
      length(
        unique(z)
      )
    }
  )
  
  
  x$n_programs_in_analysis <- as.integer(
    n_programs_gene[
      x$feature_original
    ]
  )
  
  
  # ---------------------------------------------------------------------------
  # Número de lados fenotípicos distintos en los que aparece.
  #
  # Esto NO implica expresión específica/no específica.
  # Describe únicamente membership en programas significativos.
  # ---------------------------------------------------------------------------
  
  n_sides_gene <- tapply(
    x$phenotypic_side,
    x$feature_original,
    function(z) {
      length(
        unique(z)
      )
    }
  )
  
  
  x$n_phenotypic_sides_in_analysis <-
    as.integer(
      n_sides_gene[
        x$feature_original
      ]
    )
  
  
  x$gene_program_scope <- ifelse(
    x$n_programs_in_analysis > 1L,
    "Multi-program gene",
    "Single-program gene"
  )
  
  
  x$gene_phenotype_scope <- ifelse(
    x$n_phenotypic_sides_in_analysis > 1L,
    "Cross-phenotype-side membership",
    "Single-phenotype-side membership"
  )
  
  
  # ---------------------------------------------------------------------------
  # Añadir la auditoría de Programa 3.
  # Así cada gen hereda el contexto del programa que sostiene.
  # ---------------------------------------------------------------------------
  
  if (
    !is.null(program_audit) &&
    nrow(program_audit)
  ) {
    
    keep <- intersect(
      c(
        "family",
        "program_id",
        "cross_analysis_status",
        "cross_phenotype_status",
        "strict_P1_P2_P3_triad",
        "strict_three_analysis_status"
      ),
      names(program_audit)
    )
    
    
    x <- merge(
      x,
      unique(
        program_audit[
          ,
          keep,
          drop = FALSE
        ]
      ),
      by = c(
        "family",
        "program_id"
      ),
      all.x = TRUE,
      sort = FALSE
    )
  }
  
  
  x <- x[
    order(
      x$family,
      x$phenotypic_side,
      x$program_short,
      -as.integer(
        x$is_driver_in_any_pathway
      ),
      -as.integer(
        x$n_primary_pathways
      ),
      x$feature_original
    ),
    ,
    drop = FALSE
  ]
  
  
  x
}

# ===== 7.3 RESUMIR GENES OBSERVADOS, DRIVERS Y PUENTES =====


summarize_fig4_program_genes <- function(
    gene_master
) {
  
  ids <- unique(
    gene_master$program_id
  )
  
  
  rows <- lapply(
    ids,
    function(pid) {
      
      z <- gene_master[
        gene_master$program_id == pid,
        ,
        drop = FALSE
      ]
      
      
      data.frame(
        
        family =
          z$family[1L],
        
        factor =
          z$factor[1L],
        
        factor_role =
          z$factor_role[1L],
        
        phenotypic_side =
          z$phenotypic_side[1L],
        
        factor_side =
          z$factor_side[1L],
        
        program_id =
          pid,
        
        program_short =
          z$program_short[1L],
        
        n_observed_genes =
          length(
            unique(
              z$feature_original
            )
          ),
        
        n_driver_genes =
          length(
            unique(
              z$feature_original[
                z$is_driver_in_any_pathway %in%
                  TRUE
              ]
            )
          ),
        
        n_multi_pathway_genes =
          length(
            unique(
              z$feature_original[
                as.numeric(
                  z$n_primary_pathways
                ) > 1
              ]
            )
          ),
        
        n_bridge_across_pathways =
          length(
            unique(
              z$feature_original[
                z$gene_bridge_across_pathways %in%
                  TRUE
              ]
            )
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  
  do.call(
    rbind,
    rows
  )
}



# ===== 7.4 FIGURA 4A — PATHWAY × GEN CON FUERZA DRIVER =====

make_fig4_pathway_gene_map <- function(
    pathway_gene,
    gene_master,
    family,
    program_id,
    outdir
) {
  
  z <- pathway_gene[
    pathway_gene$family == family &
      pathway_gene$program_id == program_id,
    ,
    drop = FALSE
  ]
  
  if (!nrow(z)) {
    return(invisible(NULL))
  }
  
  g <- gene_master[
    gene_master$family == family &
      gene_master$program_id == program_id,
    ,
    drop = FALSE
  ]
  
  program_short <-
    short_program_id(
      program_id
    )[1L]
  
  phenotypic_side <-
    as.character(
      z$phenotypic_side[1L]
    )
  
  factor_side <- paste0(
    z$factor[1L],
    "_",
    phenotypic_side
  )
  
  # ===== Fuerza driver por gen =====
  #
  # normalized_anchor_weight ya fue calculado y congelado aguas arriba.
  # NO se recalcula aquí.
  #
  # Para genes observados que NO son drivers, el peso queda NA.
  
  z$driver_weight <- ifelse(
    z$is_side_driver %in% TRUE &
      is.finite(
        as.numeric(
          z$normalized_anchor_weight
        )
      ),
    as.numeric(
      z$normalized_anchor_weight
    ),
    NA_real_
  )
  
  # Fuerza máxima de cada gen dentro del programa.
  # Solo se usa para ordenar visualmente las columnas.
  
  gene_strength <- tapply(
    z$driver_weight,
    z$gene_label,
    function(v) {
      if (all(is.na(v))) {
        0
      } else {
        max(
          v,
          na.rm = TRUE
        )
      }
    }
  )
  
  gene_levels <- names(
    sort(
      gene_strength,
      decreasing = TRUE
    )
  )
  
  pathway_order <- unique(
    z[
      order(
        z$library,
        z$pathway
      ),
      "pathway_label"
    ]
  )
  
  z$gene_plot <- factor(
    z$gene_label,
    levels = gene_levels
  )
  
  z$pathway_plot <- factor(
    z$pathway_label,
    levels = rev(
      pathway_order
    )
  )
  
  # ===== Matriz completa pathway × gen =====
  
  grid <- expand.grid(
    pathway_plot =
      levels(z$pathway_plot),
    
    gene_plot =
      levels(z$gene_plot),
    
    stringsAsFactors = FALSE
  )
  
  grid$membership <- FALSE
  grid$is_side_driver <- FALSE
  grid$driver_weight <- NA_real_
  
  key_grid <- paste(
    grid$pathway_plot,
    grid$gene_plot,
    sep = "\r"
  )
  
  key_obs <- paste(
    as.character(
      z$pathway_plot
    ),
    as.character(
      z$gene_plot
    ),
    sep = "\r"
  )
  
  idx <- match(
    key_obs,
    key_grid
  )
  
  valid_idx <- !is.na(idx)
  
  grid$membership[
    idx[valid_idx]
  ] <- TRUE
  
  grid$is_side_driver[
    idx[
      valid_idx &
        z$is_side_driver %in% TRUE
    ]
  ] <- TRUE
  
  # Si existiera más de una fila para la misma celda,
  # conservar la mayor contribución congelada.
  
  for (j in which(valid_idx)) {
    
    k <- idx[j]
    
    w <- suppressWarnings(
      as.numeric(
        z$driver_weight[j]
      )
    )
    
    if (!is.finite(w)) {
      next
    }
    
    if (
      is.na(grid$driver_weight[k]) ||
      w > grid$driver_weight[k]
    ) {
      grid$driver_weight[k] <- w
    }
  }
  
  grid$pathway_plot <- factor(
    grid$pathway_plot,
    levels = levels(
      z$pathway_plot
    )
  )
  
  grid$gene_plot <- factor(
    grid$gene_plot,
    levels = levels(
      z$gene_plot
    )
  )
  
  # ===== Figura =====
  #
  # Blanco = no pertenece al pathway.
  # Gris claro = gen observado del pathway, pero no driver cuantificado.
  # Intensidad = normalized_anchor_weight del driver.
  # * = driver según enriquecimiento del lado.
  
  p <- ggplot2::ggplot(
    grid,
    ggplot2::aes(
      x = gene_plot,
      y = pathway_plot
    )
  ) +
    
    ggplot2::geom_tile(
      fill = "white",
      color = "grey90",
      linewidth = 0.25
    ) +
    
    ggplot2::geom_tile(
      data = grid[
        grid$membership,
        ,
        drop = FALSE
      ],
      fill = "grey92",
      color = "grey85",
      linewidth = 0.25
    )
  
  weighted_cells <- grid[
    is.finite(grid$driver_weight) &
      grid$driver_weight > 0,
    ,
    drop = FALSE
  ]
  
  if (nrow(weighted_cells)) {
    
    p <- p +
      ggplot2::geom_tile(
        data = weighted_cells,
        ggplot2::aes(
          fill = driver_weight
        ),
        color = "grey75",
        linewidth = 0.25
      ) +
      ggplot2::scale_fill_gradient(
        low = "grey80",
        high = "black",
        name = "Normalized\nanchor weight"
      )
  }
  
  # Mantener explícita la condición de driver incluso si,
  # por alguna razón técnica, no tuvo peso utilizable en el score observado.
  
  p <- p +
    ggplot2::geom_point(
      data = grid[
        grid$is_side_driver,
        ,
        drop = FALSE
      ],
      shape = 8,
      size = 1.8,
      inherit.aes = FALSE,
      ggplot2::aes(
        x = gene_plot,
        y = pathway_plot
      )
    ) +
    
    ggplot2::labs(
      
      title = paste0(
        family,
        " — ",
        phenotype_label(
          phenotypic_side
        ),
        " — ",
        program_short
      ),
      
      subtitle = paste0(
        "Genes TX observados y contribución relativa de los drivers ",
        "a los pathways primarios del programa."
      ),
      
      x = "Gen",
      
      y = "Pathway primario",
      
      caption = paste0(
        "Gris claro = gen observado del pathway. ",
        "Intensidad = normalized_anchor_weight del driver congelado en TRAIN. ",
        "* = driver del lado fenotípico. ",
        "El peso no es logFC ni un nuevo efecto diferencial."
      )
    ) +
    
    fig_theme(9) +
    
    ggplot2::theme(
      axis.text.x =
        ggplot2::element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = 6.5
        ),
      
      axis.text.y =
        ggplot2::element_text(
          size = 7
        ),
      
      panel.grid =
        ggplot2::element_blank(),
      
      legend.position =
        "bottom"
    )
  
  dir_fig <- file.path(
    outdir,
    "figures",
    "FIG4_TX_GENES",
    "01_WITHIN_PROGRAM",
    family,
    safe_filename(
      factor_side
    ),
    safe_filename(
      program_short
    )
  )
  
  dir_dat <- file.path(
    outdir,
    "plot_data",
    "FIG4_TX_GENES",
    "01_WITHIN_PROGRAM",
    family,
    safe_filename(
      factor_side
    ),
    safe_filename(
      program_short
    )
  )
  
  dir.create(
    dir_fig,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    dir_dat,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  ggplot2::ggsave(
    filename = file.path(
      dir_fig,
      "FIG4A_PATHWAY_GENE_DRIVER_WEIGHT.png"
    ),
    plot = p,
    
    width = max(
      8,
      min(
        24,
        4 +
          0.32 *
          length(
            gene_levels
          )
      )
    ),
    
    height = max(
      5,
      min(
        18,
        3 +
          0.38 *
          length(
            pathway_order
          )
      )
    ),
    
    dpi = FIG_DPI,
    limitsize = FALSE
  )
  
  # Tabla exacta usada para dibujar el heatmap.
  write_csv(
    grid,
    file.path(
      dir_dat,
      "FIG4A_PATHWAY_GENE_DRIVER_WEIGHT.csv"
    )
  )
  
  # Mantener las tablas completas.
  write_csv(
    z,
    file.path(
      dir_dat,
      "PATHWAY_GENE_MEMBERSHIP_WITH_DRIVER_WEIGHT.csv"
    )
  )
  
  write_csv(
    g,
    file.path(
      dir_dat,
      "GENE_MASTER_TABLE.csv"
    )
  )
  
  write_csv(
    g[
      g$is_driver_in_any_pathway %in% TRUE,
      ,
      drop = FALSE
    ],
    file.path(
      dir_dat,
      "DRIVER_GENES.csv"
    )
  )
  
  write_csv(
    g[
      as.numeric(
        g$n_primary_pathways
      ) > 1,
      ,
      drop = FALSE
    ],
    file.path(
      dir_dat,
      "SHARED_PATHWAY_GENES.csv"
    )
  )
  
  invisible(
    list(
      plot = p,
      heatmap_data = grid,
      pathway_gene = z,
      gene_master = g
    )
  )
}

# ===== 7.4B FIGURA 4A-P — VERSIÓN REDUCIDA PARA PRESENTACIÓN =====

make_fig4_pathway_gene_presentation <- function(
    pathway_gene,
    family,
    program_id,
    outdir,
    max_genes =
      FIG4_PRESENTATION_MAX_GENES,
    driver_percentile_min =
      FIG4_PRESENTATION_DRIVER_PERCENTILE_MIN,
    loading_percentile_min =
      FIG4_PRESENTATION_LOADING_PERCENTILE_MIN
) {
  
  z <- pathway_gene[
    pathway_gene$family == family &
      pathway_gene$program_id == program_id,
    ,
    drop = FALSE
  ]
  
  if (!nrow(z)) {
    return(invisible(NULL))
  }
  
  
  program_short <-
    short_program_id(
      program_id
    )[1L]
  
  phenotypic_side <-
    as.character(
      z$phenotypic_side[1L]
    )
  
  factor_side <- paste0(
    z$factor[1L],
    "_",
    phenotypic_side
  )
  
  # ---------------------------------------------------------------------------
  # 1. Peso driver congelado
  # ---------------------------------------------------------------------------
  
  z$driver_weight <- ifelse(
    z$is_side_driver %in% TRUE &
      is.finite(
        as.numeric(
          z$normalized_anchor_weight
        )
      ),
    as.numeric(
      z$normalized_anchor_weight
    ),
    NA_real_
  )
  # ---------------------------------------------------------------------------
  # 2. Resumen de fuerza por gen
  # ---------------------------------------------------------------------------
  
  driver_rows <-
    z$is_side_driver %in% TRUE &
    is.finite(z$driver_weight) &
    z$driver_weight > 0 &
    is.finite(
      z$abs_loading_oriented_scaled
    )
  
  
  driver_gene <- z[
    driver_rows,
    c(
      "gene_label",
      "driver_weight",
      "abs_loading_oriented_scaled"
    ),
    drop = FALSE
  ]
  
  
  if (!nrow(driver_gene)) {
    
    message(
      "  No usable driver genes in ",
      family,
      " ",
      program_short
    )
    
    return(
      invisible(NULL)
    )
  }
  
  
  driver_summary <- stats::aggregate(
    cbind(
      driver_weight,
      abs_loading_oriented_scaled
    ) ~ gene_label,
    data = driver_gene,
    FUN = max
  )
  
  
  driver_summary$driver_percentile <-
    rank(
      driver_summary$driver_weight,
      ties.method = "average"
    ) /
    nrow(driver_summary)
  
  
  driver_summary$loading_percentile <-
    rank(
      driver_summary$abs_loading_oriented_scaled,
      ties.method = "average"
    ) /
    nrow(driver_summary)
  
  
  # ---------------------------------------------------------------------------
  # 3. Selección exclusivamente gráfica
  #
  # Si el programa ya tiene pocos drivers:
  #   -> se conservan TODOS.
  #
  # Si es grande:
  #   -> peso driver alto Y loading latente alto.
  # ---------------------------------------------------------------------------
  
  if (
    nrow(driver_summary) <=
    max_genes
  ) {
    
    selected_genes <- driver_summary
    
  } else {
    
    selected_genes <- driver_summary[
      driver_summary$driver_percentile >=
        driver_percentile_min &
        driver_summary$loading_percentile >=
        loading_percentile_min,
      ,
      drop = FALSE
    ]
    
    
    # Fallback técnico únicamente para evitar una figura vacía.
    if (!nrow(selected_genes)) {
      
      warning(
        "No genes passed both FIG4 presentation percentiles in ",
        family,
        " ",
        program_short,
        "; using strongest drivers ordered by latent loading as fallback."
      )
      
      selected_genes <- driver_summary[
        order(
          -driver_summary$driver_weight,
          -driver_summary$abs_loading_oriented_scaled
        ),
        ,
        drop = FALSE
      ]
    }
    
    
    selected_genes <- selected_genes[
      order(
        -selected_genes$driver_weight,
        -selected_genes$abs_loading_oriented_scaled
      ),
      ,
      drop = FALSE
    ]
    
    
    if (
      nrow(selected_genes) >
      max_genes
    ) {
      
      selected_genes <- selected_genes[
        seq_len(max_genes),
        ,
        drop = FALSE
      ]
    }
  }
  
  
  keep_genes <-
    as.character(
      selected_genes$gene_label
    )
  
  
  z <- z[
    z$gene_label %in%
      keep_genes,
    ,
    drop = FALSE
  ]
  
  
  gene_levels <-
    as.character(
      selected_genes$gene_label
    )
  
  pathway_order <- unique(
    z[
      order(
        z$library,
        z$pathway
      ),
      "pathway_label"
    ]
  )
  
  
  z$gene_plot <- factor(
    z$gene_label,
    levels =
      gene_levels
  )
  
  
  z$pathway_plot <- factor(
    z$pathway_label,
    levels =
      rev(
        pathway_order
      )
  )
  
  
  # ---------------------------------------------------------------------------
  # 6. Matriz completa para los genes seleccionados
  # ---------------------------------------------------------------------------
  
  grid <- expand.grid(
    
    pathway_plot =
      levels(
        z$pathway_plot
      ),
    
    gene_plot =
      levels(
        z$gene_plot
      ),
    
    stringsAsFactors =
      FALSE
  )
  
  
  grid$membership <- FALSE
  grid$is_side_driver <- FALSE
  grid$driver_weight <- NA_real_

  
  key_grid <- paste(
    grid$pathway_plot,
    grid$gene_plot,
    sep = "\r"
  )
  
  
  key_obs <- paste(
    as.character(
      z$pathway_plot
    ),
    as.character(
      z$gene_plot
    ),
    sep = "\r"
  )
  
  
  idx <- match(
    key_obs,
    key_grid
  )
  
  
  valid_idx <- !is.na(idx)
  
  
  grid$membership[
    idx[valid_idx]
  ] <- TRUE
  
  
  grid$is_side_driver[
    idx[
      valid_idx &
        z$is_side_driver %in% TRUE
    ]
  ] <- TRUE
  
  
  for (j in which(valid_idx)) {
    
    k <- idx[j]
    
    w <- suppressWarnings(
      as.numeric(
        z$driver_weight[j]
      )
    )
    

    
    if (
      is.finite(w) &&
      (
        is.na(
          grid$driver_weight[k]
        ) ||
        w >
        grid$driver_weight[k]
      )
    ) {
      
      grid$driver_weight[k] <- w
    }
    
    

  }
  
  
  grid$pathway_plot <- factor(
    grid$pathway_plot,
    levels =
      levels(
        z$pathway_plot
      )
  )
  
  
  grid$gene_plot <- factor(
    grid$gene_plot,
    levels =
      levels(
        z$gene_plot
      )
  )
  
  
  # ---------------------------------------------------------------------------
  # 7. Figura
  # ---------------------------------------------------------------------------
  
  p <- ggplot2::ggplot(
    grid,
    ggplot2::aes(
      x = gene_plot,
      y = pathway_plot
    )
  ) +
    
    # No membership
    ggplot2::geom_tile(
      fill = "white",
      color = "grey90",
      linewidth = 0.3
    ) +
    
    # Membership observada
    ggplot2::geom_tile(
      data = grid[
        grid$membership,
        ,
        drop = FALSE
      ],
      fill = "grey92",
      color = "grey85",
      linewidth = 0.3
    )
  
  
  weighted_cells <- grid[
    is.finite(
      grid$driver_weight
    ) &
      grid$driver_weight > 0,
    ,
    drop = FALSE
  ]
  
  
  if (nrow(weighted_cells)) {
    
    p <- p +
      
      ggplot2::geom_tile(
        data =
          weighted_cells,
        
        ggplot2::aes(
          fill =
            driver_weight
        ),
        
        color =
          "grey70",
        
        linewidth =
          0.3
      ) +
      
      ggplot2::scale_fill_gradient(
        low = "#DEEBF7",
        high = "#08519C",
        name = "Anchor weight"
      )
  }
  
  
  p <- p +
    
    ggplot2::labs(
      
      title = paste0(
        family,
        " — ",
        phenotype_label(
          phenotypic_side
        ),
        " — ",
        program_short
      ),
      
      subtitle ="Principales genes driver del programa."
      ,
      
      x =
        "Gen",
      
      y =
        "Pathway primario",
      
      caption = NULL
    ) +
    
    fig_theme(
      FIG4_PRESENTATION_BASE_SIZE
    ) + 
    ggplot2::theme(
      
      axis.text.x =
        ggplot2::element_text(
          angle = 60,
          hjust = 1,
          vjust = 1,
          size = FIG4_PRESENTATION_AXIS_TEXT_SIZE
          ),
      
      axis.text.y =
        ggplot2::element_text(
          size = FIG4_PRESENTATION_AXIS_TEXT_SIZE
          ),
      
      panel.grid =
        ggplot2::element_blank(),
      
      legend.position =
        "bottom"
    )
  
  
  # ---------------------------------------------------------------------------
  # 8. Salidas separadas de la figura completa
  # ---------------------------------------------------------------------------
  
  dir_fig <- file.path(
    outdir,
    "figures",
    "FIG4_TX_GENES",
    "01B_WITHIN_PROGRAM_PRESENTATION",
    family,
    safe_filename(
      factor_side
    ),
    safe_filename(
      program_short
    )
  )
  
  
  dir_dat <- file.path(
    outdir,
    "plot_data",
    "FIG4_TX_GENES",
    "01B_WITHIN_PROGRAM_PRESENTATION",
    family,
    safe_filename(
      factor_side
    ),
    safe_filename(
      program_short
    )
  )
  
  
  dir.create(
    dir_fig,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  dir.create(
    dir_dat,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  ggplot2::ggsave(
    
    filename = file.path(
      dir_fig,
      "FIG4A_PRESENTATION_DRIVER_WEIGHT.png"
    ),
    
    plot = p,
    
    width = max(
      11,
      min(
        20,
        5 +
          0.62 *
          length(
            gene_levels
          )
      )
    ),
    
    height = max(
      5,
      min(
        12,
        3 +
          0.55 *
          length(
            pathway_order
          )
      )
    ),
    
    dpi =
      FIG_DPI,
    
    limitsize =
      FALSE
  )
  
  
  write_csv(
    grid,
    file.path(
      dir_dat,
      "FIG4A_PRESENTATION_DRIVER_WEIGHT.csv"
    )
  )
  
  
  write_csv(
    z,
    file.path(
      dir_dat,
      "FIG4A_PRESENTATION_SELECTED_GENES.csv"
    )
  )
  
  
  invisible(
    list(
      plot = p,
      selected_genes =
        keep_genes,
      data = z
    )
  )
}


# ===== 7.5 GENES COMPARTIDOS ENTRE PROGRAMAS ALINEADOS =====


build_fig4_aligned_program_gene_audit <- function(
    gene_master,
    cross_alignment
) {
  
  align <- cross_alignment[
    cross_alignment$reciprocal_positive_alignment,
    ,
    drop = FALSE
  ]
  
  
  rows <- list()
  rr <- 1L
  
  
  for (
    i in seq_len(
      nrow(align)
    )
  ) {
    
    pa <- align$program_A[i]
    pb <- align$program_B[i]
    
    
    ga <- unique(
      as.character(
        gene_master$feature_original[
          gene_master$program_id ==
            pa
        ]
      )
    )
    
    
    gb <- unique(
      as.character(
        gene_master$feature_original[
          gene_master$program_id ==
            pb
        ]
      )
    )
    
    
    genes <- sort(
      union(
        ga,
        gb
      )
    )
    
    
    if (!length(genes)) {
      next
    }
    
    
    rows[[rr]] <- data.frame(
      
      family_A =
        align$family_A[i],
      
      family_B =
        align$family_B[i],
      
      phenotypic_side =
        align$phenotypic_side[i],
      
      program_A =
        pa,
      
      program_B =
        pb,
      
      feature_original =
        genes,
      
      present_A =
        genes %in% ga,
      
      present_B =
        genes %in% gb,
      
      shared_between_programs =
        genes %in% intersect(
          ga,
          gb
        ),
      
      jaccard_programs =
        align$jaccard_observed_genes[i],
      
      stringsAsFactors = FALSE
    )
    
    
    rr <- rr + 1L
  }
  
  
  if (!length(rows)) {
    return(
      data.frame()
    )
  }
  
  
  do.call(
    rbind,
    rows
  )
}

# ===== 7.6 GENES QUE SOSTIENEN TRÍADAS P1/P2/P3 =====


build_fig4_strict_triad_gene_audit <- function(
    gene_master,
    strict_triads
) {
  
  if (!nrow(strict_triads)) {
    return(
      data.frame()
    )
  }
  
  
  rows <- list()
  rr <- 1L
  
  
  for (
    i in seq_len(
      nrow(strict_triads)
    )
  ) {
    
    p1 <- strict_triads$program_P1[i]
    p2 <- strict_triads$program_P2[i]
    p3 <- strict_triads$program_P3[i]
    
    
    g1 <- unique(
      as.character(
        gene_master$feature_original[
          gene_master$program_id ==
            p1
        ]
      )
    )
    
    
    g2 <- unique(
      as.character(
        gene_master$feature_original[
          gene_master$program_id ==
            p2
        ]
      )
    )
    
    
    g3 <- unique(
      as.character(
        gene_master$feature_original[
          gene_master$program_id ==
            p3
        ]
      )
    )
    
    
    genes <- sort(
      unique(
        c(
          g1,
          g2,
          g3
        )
      )
    )
    
    
    if (!length(genes)) {
      next
    }
    
    
    present_P1 <- genes %in% g1
    present_P2 <- genes %in% g2
    present_P3 <- genes %in% g3
    
    
    support_n <-
      as.integer(
        present_P1
      ) +
      as.integer(
        present_P2
      ) +
      as.integer(
        present_P3
      )
    
    
    rows[[rr]] <- data.frame(
      
      phenotypic_side =
        strict_triads$phenotypic_side[i],
      
      program_P1 =
        p1,
      
      program_P2 =
        p2,
      
      program_P3 =
        p3,
      
      feature_original =
        genes,
      
      present_P1 =
        present_P1,
      
      present_P2 =
        present_P2,
      
      present_P3 =
        present_P3,
      
      n_analyses_supporting_gene =
        support_n,
      
      gene_scope = ifelse(
        support_n == 3L,
        "Shared P1-P2-P3",
        ifelse(
          support_n == 2L,
          "Shared two analyses",
          "Analysis-specific within triad"
        )
      ),
      
      stringsAsFactors = FALSE
    )
    
    
    rr <- rr + 1L
  }
  
  
  if (!length(rows)) {
    return(
      data.frame()
    )
  }
  
  
  do.call(
    rbind,
    rows
  )
}
# ===== 7.7 AUDITORÍA DE GENES ENTRE LADOS FENOTÍPICOS =====


build_fig4_gene_phenotype_scope <- function(
    gene_master
) {
  
  keys <- unique(
    gene_master[
      ,
      c(
        "family",
        "feature_original"
      ),
      drop = FALSE
    ]
  )
  
  
  rows <- lapply(
    seq_len(
      nrow(keys)
    ),
    function(i) {
      
      family_i <-
        keys$family[i]
      
      gene_i <-
        keys$feature_original[i]
      
      
      z <- gene_master[
        gene_master$family ==
          family_i &
          gene_master$feature_original ==
          gene_i,
        ,
        drop = FALSE
      ]
      
      
      sides <- sort(
        unique(
          as.character(
            z$phenotypic_side
          )
        )
      )
      
      
      programs <- sort(
        unique(
          as.character(
            z$program_id
          )
        )
      )
      
      
      driver_sides <- sort(
        unique(
          as.character(
            z$phenotypic_side[
              z$is_driver_in_any_pathway %in%
                TRUE
            ]
          )
        )
      )
      
      
      data.frame(
        
        family =
          family_i,
        
        feature_original =
          gene_i,
        
        n_programs =
          length(
            programs
          ),
        
        n_phenotypic_sides =
          length(
            sides
          ),
        
        phenotypic_sides =
          paste(
            sides,
            collapse = ";"
          ),
        
        programs =
          paste(
            programs,
            collapse = ";"
          ),
        
        driver_in_any_program =
          any(
            z$is_driver_in_any_pathway %in%
              TRUE
          ),
        
        driver_phenotypic_sides =
          if (
            length(
              driver_sides
            )
          ) {
            paste(
              driver_sides,
              collapse = ";"
            )
          } else {
            NA_character_
          },
        
        membership_scope =
          if (
            length(sides) == 1L
          ) {
            "Single-phenotype-side membership"
          } else {
            "Cross-phenotype-side membership"
          },
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  
  do.call(
    rbind,
    rows
  )
}

# ===== 7.8 FIGURA 4B — RESUMEN GÉNICO DE LOS PROGRAMAS =====


make_fig4_gene_summary <- function(
    program_summary,
    outdir
) {
  
  plot_dat <- rbind(
    
    data.frame(
      family =
        program_summary$family,
      
      phenotypic_side =
        program_summary$phenotypic_side,
      
      program_short =
        program_summary$program_short,
      
      gene_type =
        "Observed genes",
      
      n_genes =
        program_summary$n_observed_genes,
      
      stringsAsFactors = FALSE
    ),
    
    data.frame(
      family =
        program_summary$family,
      
      phenotypic_side =
        program_summary$phenotypic_side,
      
      program_short =
        program_summary$program_short,
      
      gene_type =
        "Driver genes",
      
      n_genes =
        program_summary$n_driver_genes,
      
      stringsAsFactors = FALSE
    ),
    
    data.frame(
      family =
        program_summary$family,
      
      phenotypic_side =
        program_summary$phenotypic_side,
      
      program_short =
        program_summary$program_short,
      
      gene_type =
        "Genes in >1 primary pathway",
      
      n_genes =
        program_summary$n_multi_pathway_genes,
      
      stringsAsFactors = FALSE
    )
  )
  
  
  plot_dat$program_label <- paste(
    plot_dat$family,
    phenotype_label(
      plot_dat$phenotypic_side
    ),
    plot_dat$program_short,
    sep = " | "
  )
  
  
  p <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(
      x = gene_type,
      y = program_label
    )
  ) +
    
    ggplot2::geom_tile(
      fill = "white",
      color = "grey70",
      linewidth = 0.5
    ) +
    
    ggplot2::geom_text(
      ggplot2::aes(
        label = n_genes
      ),
      size = 3.2
    ) +
    
    ggplot2::labs(
      
      title =
        "Pregunta 4 — genes que sostienen los programas TX",
      
      subtitle = paste0(
        "Genes observados, drivers y genes compartidos por ",
        "más de un pathway primario."
      ),
      
      x = NULL,
      y = NULL,
      
      caption = paste0(
        "Los drivers son un subconjunto de los genes observados. ",
        "Los conteos son descriptivos y no constituyen un ranking."
      )
    ) +
    
    fig_theme(10) +
    
    ggplot2::theme(
      panel.grid =
        ggplot2::element_blank(),
      
      axis.text.x =
        ggplot2::element_text(
          angle = 25,
          hjust = 1
        )
    )
  
  
  dir_fig <- file.path(
    outdir,
    "figures",
    "FIG4_TX_GENES",
    "00_INDEX"
  )
  
  
  dir.create(
    dir_fig,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  ggplot2::ggsave(
    filename = file.path(
      dir_fig,
      "FIG4B_PROGRAM_GENE_SUMMARY.png"
    ),
    plot = p,
    width = 10,
    height = max(
      7,
      3 +
        0.40 *
        length(
          unique(
            plot_dat$program_label
          )
        )
    ),
    dpi = FIG_DPI,
    limitsize = FALSE
  )
  
  
  invisible(p)
}

# ===== 7.9 EJECUTAR PREGUNTA 4 =====
if (run_question(4L)) {
message("")
message("============================================================")
message("METATEST — QUESTION 4")
message("Genes supporting frozen TX programs")
message("============================================================")


fig4_inputs <- list()
fig4_gene_master_list <- list()
fig4_program_summary_list <- list()


for (family in FAMILIES) {
  
  message(
    "\n[FIG4] ",
    family
  )
  
  
  inp4 <- load_fig4_inputs(
    bio_program_run =
      BIO_RUN,
    
    family =
      family
  )
  
  
  fig4_inputs[[family]] <-
    inp4
  
  
  gene_master <- build_fig4_gene_master(
    inputs =
      inp4,
    
    program_audit =
      master_program_audit
  )
  
  
  fig4_gene_master_list[[family]] <-
    gene_master
  
  
  program_gene_summary <-
    summarize_fig4_program_genes(
      gene_master
    )
  
  
  fig4_program_summary_list[[family]] <-
    program_gene_summary
  
  
  programs <- unique(
    gene_master$program_id
  )
  
  
  for (pid in programs) {
    
    message(
      "  ",
      pid
    )
    
    
    make_fig4_pathway_gene_map(
      pathway_gene =
        inp4$pathway_gene,
      
      gene_master =
        gene_master,
      
      family =
        family,
      
      program_id =
        pid,
      
      outdir =
        OUT_ROOT
    )
    
    make_fig4_pathway_gene_presentation(
      
      pathway_gene =
        inp4$pathway_gene,
      
      family =
        family,
      
      program_id =
        pid,
      
      outdir =
        OUT_ROOT,
      
      max_genes =
        FIG4_PRESENTATION_MAX_GENES,
      
      driver_percentile_min =
        FIG4_PRESENTATION_DRIVER_PERCENTILE_MIN,
      
      loading_percentile_min =
        FIG4_PRESENTATION_LOADING_PERCENTILE_MIN
    )
  }
}


all_gene_master <- do.call(
  rbind,
  fig4_gene_master_list
)


all_program_gene_summary <- do.call(
  rbind,
  fig4_program_summary_list
)


# ---------------------------------------------------------------------------
# Auditorías transversales
# ---------------------------------------------------------------------------

aligned_program_gene_audit <-
  build_fig4_aligned_program_gene_audit(
    gene_master =
      all_gene_master,
    
    cross_alignment =
      cross_alignment
  )


strict_triad_gene_audit <-
  build_fig4_strict_triad_gene_audit(
    gene_master =
      all_gene_master,
    
    strict_triads =
      strict_triads
  )


gene_phenotype_scope <-
  build_fig4_gene_phenotype_scope(
    all_gene_master
  )


# ---------------------------------------------------------------------------
# Directorios maestros
# ---------------------------------------------------------------------------

fig4_index_dir <- file.path(
  OUT_ROOT,
  "plot_data",
  "FIG4_TX_GENES",
  "00_INDEX"
)


fig4_cross_dir <- file.path(
  OUT_ROOT,
  "plot_data",
  "FIG4_TX_GENES",
  "02_CROSS_ANALYSIS"
)


fig4_phenotype_dir <- file.path(
  OUT_ROOT,
  "plot_data",
  "FIG4_TX_GENES",
  "03_CROSS_PHENOTYPE"
)


dir.create(
  fig4_index_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  fig4_cross_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  fig4_phenotype_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ---------------------------------------------------------------------------
# Exportar
# ---------------------------------------------------------------------------

write_csv(
  all_gene_master,
  file.path(
    fig4_index_dir,
    "FIG4_MASTER_PROGRAM_GENE_TABLE.csv"
  )
)


write_csv(
  all_program_gene_summary,
  file.path(
    fig4_index_dir,
    "FIG4_PROGRAM_GENE_SUMMARY.csv"
  )
)


write_csv(
  aligned_program_gene_audit,
  file.path(
    fig4_cross_dir,
    "FIG4_ALIGNED_PROGRAM_GENE_AUDIT.csv"
  )
)


write_csv(
  strict_triad_gene_audit,
  file.path(
    fig4_cross_dir,
    "FIG4_STRICT_TRIAD_GENE_AUDIT.csv"
  )
)


write_csv(
  gene_phenotype_scope,
  file.path(
    fig4_phenotype_dir,
    "FIG4_GENE_PHENOTYPE_MEMBERSHIP_SCOPE.csv"
  )
)


# ---------------------------------------------------------------------------
# Figura resumen
# ---------------------------------------------------------------------------

make_fig4_gene_summary(
  program_summary =
    all_program_gene_summary,
  
  outdir =
    OUT_ROOT
)

fig4_output_index <- data.frame(
  
  folder = c(
    "00_INDEX",
    "01_WITHIN_PROGRAM",
    "02_CROSS_ANALYSIS",
    "03_CROSS_PHENOTYPE"
  ),
  
  question = c(
    "¿Qué genes sostienen globalmente cada programa?",
    "¿Qué genes sostienen cada pathway dentro de cada programa y cuáles son drivers?",
    "¿Qué genes explican la conservación de programas entre P1/P2/P3?",
    "¿Qué genes aparecen únicamente en un lado o atraviesan varios lados fenotípicos?"
  ),
  
  interpretation = c(
    "Resumen maestro",
    "Pathway -> genes observados -> drivers",
    "Conservación génica entre análisis",
    "Especificidad de membership fenotípico"
  ),
  
  stringsAsFactors = FALSE
)


write_csv(
  fig4_output_index,
  file.path(
    fig4_index_dir,
    "FIG4_OUTPUT_INDEX.csv"
  )
)
message("")
message("QUESTION 4 COMPLETE")


message("Output: ", file.path(
  OUT_ROOT,
  "figures",
  "FIG4_TX_GENES"
))
}
# ===== PREGUNTA 5. ¿CUÁL ES LA ARQUITECTURA MULTIMODAL ASOCIADA A LOS PATHWAYS BIOLÓGICOS DE CADA FENOTIPO? =====


# ===== 8.1 HELPERS ESPECÍFICOS DE LA PREGUNTA 5 =====

fig5_pathway_key <- function(
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
    sep = "\r"
  )
}


fig5_as_flag <- function(x) {
  
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


fig5_view_label <- function(x) {
  
  x <- tolower(
    as.character(x)
  )
  
  dictionary <- c(
    tx = "TX",
    pr = "PR",
    me = "ME",
    cl = "CL"
  )
  
  out <- unname(
    dictionary[x]
  )
  
  bad <- is.na(out)
  
  out[bad] <- toupper(
    x[bad]
  )
  
  out
}


fig5_wrap <- function(
    x,
    width = 24L
) {
  
  vapply(
    as.character(x),
    
    function(z) {
      
      paste(
        strwrap(
          z,
          width = width
        ),
        collapse = "\n"
      )
    },
    
    character(1)
  )
}


# ===== 8.2 CARGAR LOS 26 INDIVIDUOS SIN DISTINGUIR TRAIN / TEST =====

fig5_extract_group_table <- function(g) {
  
  if (is.data.frame(g)) {
    
    assert(
      all(
        c(
          "id",
          "grupo"
        ) %in%
          names(g)
      ),
      "grupo_train/grupo_test must contain id and grupo"
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
        
        stringsAsFactors =
          FALSE
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
        
        stringsAsFactors =
          FALSE
      )
    )
  }
  
  
  stop(
    "Unrecognized group structure in ready_for_modeling",
    call. = FALSE
  )
}


fig5_combine_train_test_matrix <- function(
    train,
    test,
    view_code
) {
  
  train <- as.matrix(
    train
  )
  
  test <- as.matrix(
    test
  )
  
  
  assert(
    !is.null(
      rownames(train)
    ) &&
      !is.null(
        colnames(train)
      ),
    paste(
      view_code,
      "TRAIN matrix has no feature/sample names"
    )
  )
  
  
  assert(
    !is.null(
      rownames(test)
    ) &&
      !is.null(
        colnames(test)
      ),
    paste(
      view_code,
      "TEST matrix has no feature/sample names"
    )
  )
  
  
  all_features <- union(
    rownames(train),
    rownames(test)
  )
  
  
  all_samples <- c(
    colnames(train),
    colnames(test)
  )
  
  
  assert(
    !anyDuplicated(
      all_samples
    ),
    paste(
      "Duplicated sample IDs in pooled",
      view_code
    )
  )
  
  
  out <- matrix(
    NA_real_,
    
    nrow =
      length(
        all_features
      ),
    
    ncol =
      length(
        all_samples
      ),
    
    dimnames =
      list(
        all_features,
        all_samples
      )
  )
  
  
  out[
    rownames(train),
    colnames(train)
  ] <- train
  
  
  out[
    rownames(test),
    colnames(test)
  ] <- test
  
  
  out
}


load_fig5_pooled_data <- function(
    ready_for_modeling
) {
  
  rfm <- readRDS(
    ready_for_modeling
  )
  
  
  views <- c(
    "tx",
    "pr",
    "me",
    "cl"
  )
  
  
  matrices <- list()
  
  
  for (v in views) {
    
    assert(
      v %in%
        names(rfm),
      paste(
        "ready_for_modeling missing view",
        v
      )
    )
    
    
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
        "ready_for_modeling",
        v,
        "does not contain train/test"
      )
    )
    
    
    matrices[[v]] <-
      fig5_combine_train_test_matrix(
        train =
          rfm[[v]]$train,
        
        test =
          rfm[[v]]$test,
        
        view_code =
          v
      )
  }
  
  
  assert(
    all(
      c(
        "grupo_train",
        "grupo_test"
      ) %in%
        names(rfm)
    ),
    "ready_for_modeling does not contain grupo_train/grupo_test"
  )
  
  
  sample_meta <- rbind(
    fig5_extract_group_table(
      rfm$grupo_train
    ),
    
    fig5_extract_group_table(
      rfm$grupo_test
    )
  )
  
  
  assert(
    !anyDuplicated(
      sample_meta$sample_id
    ),
    "Duplicated sample IDs in pooled metadata"
  )
  
  
  sample_meta$grupo <- factor(
    sample_meta$grupo,
    levels = c(
      "NP",
      "FD",
      "SO"
    )
  )
  
  
  # Todas las vistas deben usar las mismas 26 muestras.
  reference_ids <- sample_meta$sample_id
  
  
  for (v in views) {
    
    assert(
      all(
        reference_ids %in%
          colnames(
            matrices[[v]]
          )
      ),
      paste(
        "Some pooled samples are missing from view",
        v
      )
    )
    
    
    matrices[[v]] <-
      matrices[[v]][
        ,
        reference_ids,
        drop = FALSE
      ]
  }
  
  
  list(
    matrices =
      matrices,
    
    sample_meta =
      sample_meta
  )
}
# ===== 8.2B RECONSTRUIR COMPONENTE LATENTE MULTIMODAL EN LOS 26 INDIVIDUOS =====

load_fig5_latent_signal <- function(
    crossomic_run,
    latent_run,
    sample_meta
) {
  
  # ---------------------------------------------------------------------------
  # Loadings ya orientados y escalados por el script 05:
  #
  #   L_jf =
  #     orientation_sign_f *
  #     train_scale_after_orientation_f *
  #     W_jf
  #
  # NO se recalculan loadings ni se vuelve a cargar/refitear MOFA.
  # ---------------------------------------------------------------------------
  
  loading_path <- file.path(
    crossomic_run,
    "tables",
    "01_ORIENTED_SCALED_LOADINGS_ALL_VIEWS.csv"
  )
  
  
  loading_long <- read_csv_checked(
    loading_path,
    required = c(
      "view_code",
      "feature_original",
      "factor",
      "loading_oriented_scaled"
    )
  )
  
  
  loading_long$view_code <- tolower(
    as.character(
      loading_long$view_code
    )
  )
  
  
  loading_long$feature_original <- as.character(
    loading_long$feature_original
  )
  
  
  loading_long$factor <- as.character(
    loading_long$factor
  )
  
  
  loading_long$loading_oriented_scaled <- suppressWarnings(
    as.numeric(
      loading_long$loading_oriented_scaled
    )
  )
  
  
  factors <- unique(
    loading_long$factor
  )
  
  
  factors <- factors[
    order(
      factor_number(factors),
      factors
    )
  ]
  
  
  # ---------------------------------------------------------------------------
  # Scores latentes TRAIN + TEST ya proyectado.
  #
  # TEST NO se vuelve a proyectar.
  # ---------------------------------------------------------------------------
  
  score_path <- file.path(
    latent_run,
    "tables",
    "02_LATENT_SCORES_TRAIN_AND_PROJECTED_TEST.csv"
  )
  
  
  scores <- read_csv_checked(
    score_path,
    required = c(
      "sample",
      "split",
      "group"
    )
  )
  
  
  score_cols <- paste0(
    factors,
    "__oriented_scaled"
  )
  
  
  assert(
    all(
      score_cols %in%
        names(scores)
    ),
    paste(
      "Missing oriented-scaled latent score columns:",
      paste(
        setdiff(
          score_cols,
          names(scores)
        ),
        collapse = ", "
      )
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Usar exactamente el mismo orden de los 26 individuos que en OBSERVED.
  # ---------------------------------------------------------------------------
  
  sample_ids <- as.character(
    sample_meta$sample_id
  )
  
  
  idx <- match(
    sample_ids,
    as.character(
      scores$sample
    )
  )
  
  
  assert(
    all(
      !is.na(idx)
    ),
    "Some pooled samples are missing from latent score table"
  )
  
  
  scores <- scores[
    idx,
    ,
    drop = FALSE
  ]
  
  
  assert(
    identical(
      as.character(scores$sample),
      sample_ids
    ),
    "Latent scores could not be aligned to pooled sample order"
  )
  
  
  # ---------------------------------------------------------------------------
  # Z = factores x individuos.
  # ---------------------------------------------------------------------------
  
  Z <- t(
    as.matrix(
      scores[
        ,
        score_cols,
        drop = FALSE
      ]
    )
  )
  
  
  storage.mode(Z) <- "double"
  
  
  rownames(Z) <- factors
  
  colnames(Z) <- sample_ids
  
  
  assert(
    all(
      is.finite(Z)
    ),
    "Non-finite oriented-scaled latent scores"
  )
  
  
  # ---------------------------------------------------------------------------
  # Reconstrucción CENTRADA del componente latente:
  #
  #   X_latent = L %*% Z_oriented_scaled
  #
  # Algebraicamente:
  #
  #   L Z_scaled = W Z_native - constante_feature
  #
  # Esa constante específica de feature NO altera Spearman.
  #
  # Por tanto esta matriz es apropiada para comparar la estructura
  # de covariación LATENTE con la estructura OBSERVED.
  #
  # Importante:
  # representa solamente los factores auditados presentes en "factors".
  # En el análisis actual: Factor1 y Factor2.
  # ---------------------------------------------------------------------------
  
  matrices <- list()
  
  
  for (
    v in c(
      "tx",
      "pr",
      "me",
      "cl"
    )
  ) {
    
    d <- loading_long[
      loading_long$view_code ==
        v,
      ,
      drop = FALSE
    ]
    
    
    assert(
      nrow(d) > 0L,
      paste(
        "No oriented-scaled loadings found for view",
        v
      )
    )
    
    
    features <- unique(
      d$feature_original
    )
    
    
    L <- matrix(
      NA_real_,
      nrow = length(features),
      ncol = length(factors),
      dimnames = list(
        features,
        factors
      )
    )
    
    
    for (f in factors) {
      
      df <- d[
        d$factor ==
          f,
        c(
          "feature_original",
          "loading_oriented_scaled"
        ),
        drop = FALSE
      ]
      
      
      assert(
        !anyDuplicated(
          df$feature_original
        ),
        paste(
          "Duplicated feature loadings in",
          v,
          f
        )
      )
      
      
      ii <- match(
        features,
        df$feature_original
      )
      
      
      assert(
        all(
          !is.na(ii)
        ),
        paste(
          "Missing loadings in",
          v,
          f
        )
      )
      
      
      L[, f] <- as.numeric(
        df$loading_oriented_scaled[ii]
      )
    }
    
    
    assert(
      all(
        is.finite(L)
      ),
      paste(
        "Non-finite latent loadings in view",
        v
      )
    )
    
    
    X_latent <- L %*% Z
    
    
    rownames(X_latent) <- features
    
    colnames(X_latent) <- sample_ids
    
    
    matrices[[v]] <- X_latent
  }
  
  
  list(
    matrices = matrices,
    factors = factors,
    loading_source = loading_path,
    score_source = score_path
  )
}

# ===== 8.3 CARGAR PROGRAMAS, PATHWAYS Y CONTEXTO MULTIMODAL CONGELADO =====

load_fig5_inputs <- function(
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
      "crossomic_run",
      "observed_train_run",
      "latent_run"
    )
  )
  
  
  assert(
    nrow(audit) >= 1L,
    paste(
      "Empty program audit for",
      family
    )
  )
  
  
  crossomic_run <- as.character(
    audit$crossomic_run[1L]
  )
  latent_run <- as.character(
    audit$latent_run[1L]
  )
  
  
  assert(
    dir.exists(
      latent_run
    ),
    paste(
      "Latent run does not exist:",
      latent_run
    )
  )
  
  observed_train_run <- as.character(
    audit$observed_train_run[1L]
  )
  
  
  # ---------------------------------------------------------------------------
  # Ready-for-modeling exacto usado por la corrida.
  # ---------------------------------------------------------------------------
  
  observed_audit <- read_csv_checked(
    file.path(
      observed_train_run,
      "tables",
      "00_INPUT_AUDIT.csv"
    ),
    
    required =
      "ready_for_modeling"
  )
  
  
  ready_for_modeling <-
    as.character(
      observed_audit$ready_for_modeling[1L]
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
  
  
  pooled <- load_fig5_pooled_data(
    ready_for_modeling
  )
  
  latent_pooled <- load_fig5_latent_signal(
    crossomic_run =
      crossomic_run,
    
    latent_run =
      latent_run,
    
    sample_meta =
      pooled$sample_meta
  )
  # ---------------------------------------------------------------------------
  # Programas.
  # ---------------------------------------------------------------------------
  
  programs <- read_csv_checked(
    file.path(
      family_dir,
      "06_TX_PROGRAM_DEFINITION.csv"
    ),
    
    required = c(
      "family",
      "factor",
      "factor_role",
      "phenotypic_side",
      "program_id",
      "representative_library",
      "representative_pathway",
      "n_primary_pathways"
    )
  )
  
  
  programs$family <-
    family
  
  
  programs$program_short <-
    short_program_id(
      programs$program_id
    )
  
  
  programs$factor_side <- paste0(
    programs$factor,
    "_",
    programs$phenotypic_side
  )
  
  
  # ---------------------------------------------------------------------------
  # Pathways PRIMARY.
  # ---------------------------------------------------------------------------
  
  memberships <- read_csv_checked(
    file.path(
      family_dir,
      "07_TX_PROGRAM_PATHWAY_MEMBERSHIP.csv"
    ),
    
    required = c(
      "program_id",
      "membership_type",
      "library",
      "pathway",
      "factor",
      "oriented_side",
      "factor_role",
      "phenotypic_side"
    )
  )
  
  
  memberships <- memberships[
    memberships$membership_type ==
      "PRIMARY",
    ,
    drop = FALSE
  ]
  
  
  memberships$family <-
    family
  
  
  memberships$program_short <-
    short_program_id(
      memberships$program_id
    )
  
  
  memberships$factor_side <- paste0(
    memberships$factor,
    "_",
    memberships$phenotypic_side
  )
  
  
  memberships$pathway_label <-
    clean_pathway_label(
      memberships$pathway
    )
  
  
  memberships$pathway_key <-
    fig5_pathway_key(
      memberships$library,
      memberships$pathway,
      memberships$factor,
      memberships$oriented_side
    )
  
  
  # ---------------------------------------------------------------------------
  # Contexto latente PR / ME / CL ya definido aguas arriba.
  #
  # Aquí NO se usan las correlaciones para decidir qué features pertenecen
  # al contexto del pathway.
  # ---------------------------------------------------------------------------
  
  links <- read_csv_checked(
    file.path(
      crossomic_run,
      "tables",
      "04_PATHWAY_FEATURE_CROSSOMIC_LINKS_ALL.csv"
    ),
    
    required = c(
      "family",
      "library",
      "pathway",
      "factor",
      "oriented_side",
      "view_code",
      "feature_original",
      "feature_model",
      "pathway_feature_cosine"
    )
  )
  
  
  links <- links[
    links$family ==
      family &
      tolower(
        as.character(
          links$view_code
        )
      ) %in%
      c(
        "pr",
        "me",
        "cl"
      ),
    ,
    drop = FALSE
  ]
  
  
  links$view_code <-
    tolower(
      as.character(
        links$view_code
      )
    )
  
  
  links$feature_original <-
    as.character(
      links$feature_original
    )
  
  
  links$feature_model <-
    as.character(
      links$feature_model
    )
  
  
  links$pathway_feature_cosine <-
    suppressWarnings(
      as.numeric(
        links$pathway_feature_cosine
      )
    )
  
  
  links$pathway_key <-
    fig5_pathway_key(
      links$library,
      links$pathway,
      links$factor,
      links$oriented_side
    )
  
  
  # Solo pathways que pertenecen a programas PRIMARY.
  pmap <- memberships[
    ,
    c(
      "program_id",
      "program_short",
      "library",
      "pathway",
      "pathway_label",
      "factor",
      "oriented_side",
      "factor_role",
      "phenotypic_side",
      "factor_side",
      "pathway_key"
    ),
    drop = FALSE
  ]
  
  
  links <- merge(
    links,
    pmap,
    
    by = c(
      "library",
      "pathway",
      "factor",
      "oriented_side",
      "pathway_key"
    ),
    
    all = FALSE,
    sort = FALSE
  )
  
  
  # ---------------------------------------------------------------------------
  # Scoring specification congelada.
  #
  # Se usa únicamente para aplicar exactamente el MISMO score de pathway
  # a los 26 individuos. No se refitea ni renormaliza nada.
  # ---------------------------------------------------------------------------
  
  observed_object_path <- file.path(
    observed_train_run,
    "rds",
    "OBSERVED_TRAIN_CORROBORATION_OBJECT.rds"
  )
  
  
  assert(
    file.exists(
      observed_object_path
    ),
    paste(
      "Missing observed object:",
      observed_object_path
    )
  )
  
  
  observed_object <- readRDS(
    observed_object_path
  )
  
  
  scoring_spec <- as.data.frame(
    observed_object$pathway_scoring_spec,
    stringsAsFactors = FALSE
  )
  
  
  required_spec <- c(
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
    all(
      required_spec %in%
        names(
          scoring_spec
        )
    ),
    paste(
      "Missing pathway scoring columns:",
      paste(
        setdiff(
          required_spec,
          names(
            scoring_spec
          )
        ),
        collapse = ", "
      )
    )
  )
  
  
  scoring_spec$feature_original <-
    as.character(
      scoring_spec$feature_original
    )
  
  
  scoring_spec$pathway_key <-
    fig5_pathway_key(
      scoring_spec$library,
      scoring_spec$pathway,
      scoring_spec$factor,
      scoring_spec$oriented_side
    )
  
  
  list(
    family =
      family,
    
    programs =
      programs,
    
    memberships =
      memberships,
    
    non_tx_links =
      links,
    
    scoring_spec =
      scoring_spec,
    
    matrices =
      pooled$matrices,
    
    latent_matrices =
      latent_pooled$matrices,
    
    latent_factors =
      latent_pooled$factors,
    
    sample_meta =
      pooled$sample_meta,
    
    ready_for_modeling =
      ready_for_modeling,
    
    latent_run =
      latent_run,
    
    crossomic_run =
      crossomic_run
  )
}


# ===== 8.4 APLICAR LOS SCORES DE PATHWAY CONGELADOS A LOS 26 INDIVIDUOS =====

build_fig5_pooled_pathway_scores <- function(
    inputs
) {
  
  spec <- inputs$scoring_spec
  
  Xtx <- inputs$matrices$tx
  
  sample_meta <- inputs$sample_meta
  
  
  keys <- unique(
    inputs$memberships$pathway_key
  )
  
  
  rows <- list()
  
  rr <- 1L
  
  
  for (k in keys) {
    
    d <- spec[
      spec$pathway_key ==
        k,
      ,
      drop = FALSE
    ]
    
    
    if (!nrow(d)) {
      next
    }
    
    
    usable <-
      fig5_as_flag(
        d$usable_for_score
      ) &
      is.finite(
        as.numeric(
          d$train_mean
        )
      ) &
      is.finite(
        as.numeric(
          d$train_sd
        )
      ) &
      as.numeric(
        d$train_sd
      ) > 0 &
      is.finite(
        as.numeric(
          d$normalized_anchor_weight
        )
      ) &
      as.numeric(
        d$normalized_anchor_weight
      ) > 0 &
      d$feature_original %in%
      rownames(Xtx)
    
    
    d <- d[
      usable,
      ,
      drop = FALSE
    ]
    
    
    if (!nrow(d)) {
      next
    }
    
    
    sample_ids <- sample_meta$sample_id
    
    
    Z <- matrix(
      NA_real_,
      
      nrow =
        nrow(d),
      
      ncol =
        length(
          sample_ids
        ),
      
      dimnames =
        list(
          d$feature_original,
          sample_ids
        )
    )
    
    
    for (
      j in seq_len(
        nrow(d)
      )
    ) {
      
      values <- as.numeric(
        Xtx[
          d$feature_original[j],
          sample_ids
        ]
      )
      
      
      Z[j, ] <-
        (
          values -
            as.numeric(
              d$train_mean[j]
            )
        ) /
        as.numeric(
          d$train_sd[j]
        )
    }
    
    
    weights <- as.numeric(
      d$normalized_anchor_weight
    )
    
    
    score <- as.numeric(
      crossprod(
        weights,
        Z
      )
    )
    
    
    m <- inputs$memberships[
      inputs$memberships$pathway_key ==
        k,
      ,
      drop = FALSE
    ]
    
    
    assert(
      nrow(m) >= 1L,
      paste(
        "Could not map pathway score",
        k
      )
    )
    
    
    rows[[rr]] <- data.frame(
      
      family =
        inputs$family,
      
      program_id =
        m$program_id[1L],
      
      program_short =
        m$program_short[1L],
      
      factor =
        m$factor[1L],
      
      factor_role =
        m$factor_role[1L],
      
      phenotypic_side =
        m$phenotypic_side[1L],
      
      factor_side =
        m$factor_side[1L],
      
      library =
        m$library[1L],
      
      pathway =
        m$pathway[1L],
      
      pathway_label =
        m$pathway_label[1L],
      
      pathway_key =
        k,
      
      sample_id =
        sample_ids,
      
      grupo =
        as.character(
          sample_meta$grupo
        ),
      
      pathway_score =
        score,
      
      stringsAsFactors =
        FALSE
    )
    
    
    rr <- rr + 1L
  }
  
  
  if (!length(rows)) {
    
    return(
      data.frame()
    )
  }
  
  
  do.call(
    rbind,
    rows
  )
}


# ===== 8.5 CONSTRUIR CATÁLOGO TX + PR + ME + CL POR PATHWAY =====

build_fig5_feature_catalog <- function(
    inputs,
    tx_pathway_gene
) {
  
  # ---------------------------------------------------------------------------
  # TX
  # ---------------------------------------------------------------------------
  
  tx <- tx_pathway_gene[
    tx_pathway_gene$family ==
      inputs$family,
    ,
    drop = FALSE
  ]
  
  
  tx$pathway_key <-
    fig5_pathway_key(
      tx$library,
      tx$pathway,
      tx$factor,
      tx$oriented_side
    )
  
  
  tx$view_code <-
    "tx"
  
  
  tx$view_label <-
    "TX"
  
  
  tx$feature_label <-
    as.character(
      tx$gene_label
    )
  
  
  tx$feature_model <-
    as.character(
      tx$feature_original
    )
  
  
  tx$is_tx_member <-
    TRUE
  
  
  tx$is_tx_driver <-
    tx$is_side_driver %in%
    TRUE
  
  
  tx$tx_driver_weight <-
    suppressWarnings(
      as.numeric(
        tx$normalized_anchor_weight
      )
    )
  
  
  tx$pathway_feature_cosine <-
    NA_real_
  
  
  tx_catalog <- tx[
    ,
    c(
      "family",
      "program_id",
      "program_short",
      "library",
      "pathway",
      "pathway_label",
      "pathway_key",
      "factor",
      "factor_role",
      "oriented_side",
      "phenotypic_side",
      "factor_side",
      "view_code",
      "view_label",
      "feature_original",
      "feature_model",
      "feature_label",
      "is_tx_member",
      "is_tx_driver",
      "tx_driver_weight",
      "pathway_feature_cosine"
    ),
    drop = FALSE
  ]
  
  
  # ---------------------------------------------------------------------------
  # PR / ME / CL
  # ---------------------------------------------------------------------------
  
  z <- inputs$non_tx_links
  
  
  z$view_label <-
    fig5_view_label(
      z$view_code
    )
  
  
  z$feature_label <-
    as.character(
      z$feature_original
    )
  
  
  bad <-
    is.na(
      z$feature_label
    ) |
    !nzchar(
      z$feature_label
    ) |
    z$feature_label ==
    "NA"
  
  
  z$feature_label[
    bad
  ] <- z$feature_model[
    bad
  ]
  
  
  z$is_tx_member <-
    FALSE
  
  
  z$is_tx_driver <-
    FALSE
  
  
  z$tx_driver_weight <-
    NA_real_
  
  
  non_tx_catalog <- z[
    ,
    c(
      "family",
      "program_id",
      "program_short",
      "library",
      "pathway",
      "pathway_label",
      "pathway_key",
      "factor",
      "factor_role",
      "oriented_side",
      "phenotypic_side",
      "factor_side",
      "view_code",
      "view_label",
      "feature_original",
      "feature_model",
      "feature_label",
      "is_tx_member",
      "is_tx_driver",
      "tx_driver_weight",
      "pathway_feature_cosine"
    ),
    drop = FALSE
  ]
  
  
  catalog <- rbind(
    tx_catalog,
    non_tx_catalog
  )
  
  
  catalog$feature_original <-
    as.character(
      catalog$feature_original
    )
  
  
  catalog$feature_key <- paste(
    catalog$view_code,
    catalog$feature_original,
    sep = "\r"
  )
  
  
  # ---------------------------------------------------------------------------
  # Evitar duplicados exactos pathway-feature.
  # ---------------------------------------------------------------------------
  
  duplicate_key <- paste(
    catalog$pathway_key,
    catalog$feature_key,
    sep = "\r"
  )
  
  
  catalog <- catalog[
    !duplicated(
      duplicate_key
    ),
    ,
    drop = FALSE
  ]
  
  
  # ---------------------------------------------------------------------------
  # ¿Está realmente disponible la feature en los 26?
  # ---------------------------------------------------------------------------
  
  catalog$available_in_pooled <- mapply(
    
    function(v, f) {
      
      v %in%
        names(
          inputs$matrices
        ) &&
        f %in%
        rownames(
          inputs$matrices[[v]]
        )
    },
    
    catalog$view_code,
    catalog$feature_original
  )
  
  
  catalog
}


# ===== 8.6 VARIABILIDAD DESCRIPTIVA DE CADA FEATURE EN LOS 26 =====

add_fig5_pooled_variance <- function(
    catalog,
    matrices
) {
  
  feature_meta <- unique(
    catalog[
      ,
      c(
        "feature_key",
        "view_code",
        "view_label",
        "feature_original",
        "feature_label",
        "available_in_pooled"
      ),
      drop = FALSE
    ]
  )
  
  
  variance <- rep(
    NA_real_,
    nrow(
      feature_meta
    )
  )
  
  
  for (
    i in seq_len(
      nrow(
        feature_meta
      )
    )
  ) {
    
    if (
      !feature_meta$available_in_pooled[i]
    ) {
      next
    }
    
    
    v <- feature_meta$view_code[i]
    
    f <- feature_meta$feature_original[i]
    
    
    values <- as.numeric(
      matrices[[v]][
        f,
        ,
        drop = TRUE
      ]
    )
    
    
    values <- values[
      is.finite(
        values
      )
    ]
    
    
    if (
      length(values) >=
      2L
    ) {
      
      variance[i] <- stats::var(
        values
      )
    }
  }
  
  
  feature_meta$pooled_variance <-
    variance
  
  
  # ---------------------------------------------------------------------------
  # Percentil de variabilidad DENTRO de cada vista.
  #
  # Esto evita comparar directamente la escala numérica de TX, ME, PR y CL.
  # ---------------------------------------------------------------------------
  
  feature_meta$variance_percentile_within_view <-
    NA_real_
  
  
  for (
    v in unique(
      feature_meta$view_code
    )
  ) {
    
    idx <- which(
      feature_meta$view_code ==
        v &
        is.finite(
          feature_meta$pooled_variance
        )
    )
    
    
    if (!length(idx)) {
      next
    }
    
    
    feature_meta$variance_percentile_within_view[
      idx
    ] <-
      rank(
        feature_meta$pooled_variance[
          idx
        ],
        ties.method = "average"
      ) /
      length(idx)
  }
  
  
  catalog <- merge(
    catalog,
    feature_meta[
      ,
      c(
        "feature_key",
        "pooled_variance",
        "variance_percentile_within_view"
      ),
      drop = FALSE
    ],
    
    by =
      "feature_key",
    
    all.x =
      TRUE,
    
    sort =
      FALSE
  )
  
  
  catalog
}


# ===== 8.7 SELECCIÓN ADAPTATIVA EXCLUSIVAMENTE PARA PRESENTACIÓN =====
select_fig5_presentation_features <- function(
    catalog,
    max_features =
      FIG5_PRESENTATION_MAX_FEATURES,
    min_per_view =
      FIG5_PRESENTATION_MIN_PER_VIEW
) {
  
  feature_meta <- unique(
    catalog[
      ,
      c(
        "feature_key",
        "view_code",
        "view_label",
        "feature_original",
        "feature_label",
        "available_in_pooled",
        "pooled_variance",
        "variance_percentile_within_view"
      ),
      drop = FALSE
    ]
  )
  
  
  driver_any <- tapply(
    catalog$is_tx_driver,
    catalog$feature_key,
    function(x) {
      any(
        x %in% TRUE
      )
    }
  )
  
  
  feature_meta$is_tx_driver_any <-
    as.logical(
      driver_any[
        feature_meta$feature_key
      ]
    )
  
  
  feature_meta$is_tx_driver_any[
    is.na(
      feature_meta$is_tx_driver_any
    )
  ] <- FALSE
  
  
  feature_meta <- feature_meta[
    feature_meta$available_in_pooled,
    ,
    drop = FALSE
  ]
  
  
  feature_meta$selected_for_presentation <-
    FALSE
  
  feature_meta$presentation_selection_reason <-
    NA_character_
  
  
  # Si ya es suficientemente pequeña: conservar todo.
  if (
    nrow(feature_meta) <=
    max_features
  ) {
    
    feature_meta$selected_for_presentation <-
      TRUE
    
    feature_meta$presentation_selection_reason <-
      "All features retained: pathway/program already readable"
    
    return(
      feature_meta
    )
  }
  
  
  # -------------------------------------------------------------------------
  # Cuando hay que reducir:
  # reservar representación de cada vista presente.
  # -------------------------------------------------------------------------
  
  view_order <- c(
    "tx",
    "pr",
    "me",
    "cl"
  )
  
  
  views_present <- view_order[
    view_order %in%
      unique(
        feature_meta$view_code
      )
  ]
  
  
  reserve_n <- min(
    min_per_view,
    max(
      1L,
      floor(
        max_features /
          length(
            views_present
          )
      )
    )
  )
  
  
  selected <- feature_meta[
    0,
    ,
    drop = FALSE
  ]
  
  
  for (v in views_present) {
    
    d <- feature_meta[
      feature_meta$view_code ==
        v,
      ,
      drop = FALSE
    ]
    
    
    if (
      v == "tx"
    ) {
      
      # TX: priorizar drivers y después variabilidad.
      d <- d[
        order(
          -as.integer(
            d$is_tx_driver_any
          ),
          -d$variance_percentile_within_view,
          -d$pooled_variance,
          d$feature_label
        ),
        ,
        drop = FALSE
      ]
      
    } else {
      
      # PR / ME / CL: mayor variabilidad dentro de la propia vista.
      d <- d[
        order(
          -d$variance_percentile_within_view,
          -d$pooled_variance,
          d$feature_label
        ),
        ,
        drop = FALSE
      ]
    }
    
    
    take <- min(
      reserve_n,
      nrow(d)
    )
    
    
    if (
      take > 0L
    ) {
      
      d <- d[
        seq_len(take),
        ,
        drop = FALSE
      ]
      
      d$selected_for_presentation <-
        TRUE
      
      d$presentation_selection_reason <-
        paste0(
          toupper(v),
          " retained in per-view presentation reserve"
        )
      
      selected <- rbind(
        selected,
        d
      )
    }
  }
  
  
  # -------------------------------------------------------------------------
  # Rellenar los huecos restantes con las features más variables.
  # -------------------------------------------------------------------------
  
  remaining_slots <-
    max_features -
    nrow(selected)
  
  
  if (
    remaining_slots > 0L
  ) {
    
    rest <- feature_meta[
      !feature_meta$feature_key %in%
        selected$feature_key,
      ,
      drop = FALSE
    ]
    
    
    rest$tx_driver_priority <-
      as.integer(
        rest$view_code == "tx" &
          rest$is_tx_driver_any
      )
    
    
    rest <- rest[
      order(
        -rest$variance_percentile_within_view,
        -rest$tx_driver_priority,
        -rest$pooled_variance,
        match(
          rest$view_code,
          view_order
        ),
        rest$feature_label
      ),
      ,
      drop = FALSE
    ]
    
    
    take <- min(
      remaining_slots,
      nrow(rest)
    )
    
    
    if (
      take > 0L
    ) {
      
      rest <- rest[
        seq_len(take),
        ,
        drop = FALSE
      ]
      
      rest$selected_for_presentation <-
        TRUE
      
      rest$presentation_selection_reason <-
        "Added by high within-view variability after per-view reserve"
      rest$tx_driver_priority <- NULL
      selected <- rbind(
        selected,
        rest
      )
    }
  }
  
  
  feature_meta$selected_for_presentation <-
    feature_meta$feature_key %in%
    selected$feature_key
  
  
  idx <- match(
    feature_meta$feature_key,
    selected$feature_key
  )
  
  
  feature_meta$presentation_selection_reason <-
    selected$presentation_selection_reason[
      idx
    ]
  
  
  feature_meta
}
# ===== 8.8 CALCULAR SPEARMAN PATHWAY ↔ PR / ME / CL EN LOS 26 Y POR FENOTIPO =====

compute_fig5_pathway_feature_correlations <- function(
    catalog,
    pathway_scores,
    matrices,
    sample_meta,
    scope = c(
      "GLOBAL",
      "NP",
      "FD",
      "SO"
    )
) {
  
  scope <- match.arg(
    scope
  )
  
  
  z <- catalog[
    catalog$view_code %in%
      c(
        "pr",
        "me",
        "cl"
      ) &
      catalog$available_in_pooled,
    ,
    drop = FALSE
  ]
  
  
  if (!nrow(z)) {
    return(
      data.frame()
    )
  }
  
  
  rows <- list()
  
  rr <- 1L
  
  
  for (
    i in seq_len(
      nrow(z)
    )
  ) {
    
    row <- z[
      i,
      ,
      drop = FALSE
    ]
    
    
    s <- pathway_scores[
      pathway_scores$pathway_key ==
        row$pathway_key,
      ,
      drop = FALSE
    ]
    
    
    if (!nrow(s)) {
      next
    }
    
    
    if (
      scope !=
      "GLOBAL"
    ) {
      
      s <- s[
        s$grupo ==
          scope,
        ,
        drop = FALSE
      ]
    }
    
    
    if (!nrow(s)) {
      next
    }
    
    
    X <- matrices[[as.character(row$view_code[1L])]]
    
    feature <-
      row$feature_original
    
    
    if (
      !feature %in%
      rownames(X)
    ) {
      next
    }
    
    
    x <- as.numeric(
      X[
        feature,
        s$sample_id
      ]
    )
    
    
    y <- as.numeric(
      s$pathway_score
    )
    
    
    ok <-
      is.finite(x) &
      is.finite(y)
    
    
    n_complete <- sum(
      ok
    )
    
    
    rho <- if (
      n_complete >=
      3L &&
      stats::sd(
        x[ok]
      ) > 0 &&
      stats::sd(
        y[ok]
      ) > 0
    ) {
      
      suppressWarnings(
        stats::cor(
          x[ok],
          y[ok],
          method =
            "spearman"
        )
      )
      
    } else {
      
      NA_real_
    }
    
    
    rows[[rr]] <- data.frame(
      
      family =
        row$family,
      
      program_id =
        row$program_id,
      
      program_short =
        row$program_short,
      
      factor =
        row$factor,
      
      factor_role =
        row$factor_role,
      
      phenotypic_side =
        row$phenotypic_side,
      
      factor_side =
        row$factor_side,
      
      library =
        row$library,
      
      pathway =
        row$pathway,
      
      pathway_label =
        row$pathway_label,
      
      pathway_key =
        row$pathway_key,
      
      view_code =
        row$view_code,
      
      view_label =
        row$view_label,
      
      feature_original =
        row$feature_original,
      
      feature_label =
        row$feature_label,
      
      scope =
        scope,
      
      n_complete =
        n_complete,
      
      rho_spearman =
        rho,
      
      pathway_feature_cosine =
        row$pathway_feature_cosine,
      
      stringsAsFactors =
        FALSE
    )
    
    
    rr <- rr + 1L
  }
  
  
  if (!length(rows)) {
    
    return(
      data.frame()
    )
  }
  
  
  do.call(
    rbind,
    rows
  )
}


# ===== 8.9 FIGURA 5A — MAPA INTEGRADO PATHWAY × TX / PR / ME / CL =====

make_fig5_program_architecture_map <- function(
    catalog,
    pathway_correlations,
    memberships,
    family,
    program_id,
    scope,
    outdir,
    branch
) {
  
  z <- catalog[
    catalog$family ==
      family &
      catalog$program_id ==
      program_id,
    ,
    drop = FALSE
  ]
  
  
  if (!nrow(z)) {
    
    return(
      invisible(NULL)
    )
  }
  
  
  pmap <- memberships[
    memberships$family ==
      family &
      memberships$program_id ==
      program_id,
    ,
    drop = FALSE
  ]
  
  
  if (!nrow(pmap)) {
    
    return(
      invisible(NULL)
    )
  }
  
  
  meta <- unique(
    z[
      ,
      c(
        "family",
        "program_id",
        "program_short",
        "factor_side",
        "phenotypic_side"
      ),
      drop = FALSE
    ]
  )
  
  
  program_short <-
    meta$program_short[1L]
  
  
  factor_side <-
    meta$factor_side[1L]
  
  
  phenotypic_side <-
    meta$phenotypic_side[1L]
  
  branch_tag <- if (
    branch == "02_PRESENTATION"
  ) {
    "PRESENTATION"
  } else {
    "ALL"
  }
  
  
  fig5a_png_target <- file.path(
    outdir,
    "figures",
    "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
    branch,
    family,
    safe_filename(factor_side),
    safe_filename(program_short),
    "01_PROGRAM_PATHWAY_MAP",
    scope,
    paste0(
      "FIG5A_",
      scope,
      "_PROGRAM_PATHWAY_MULTIMODAL_MAP.png"
    )
  )
  
  
  fig5a_csv_target <- file.path(
    outdir,
    "plot_data",
    "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
    branch,
    family,
    safe_filename(factor_side),
    safe_filename(program_short),
    "01_PROGRAM_PATHWAY_MAP",
    scope,
    paste0(
      "FIG5A_",
      scope,
      "_PROGRAM_PATHWAY_MULTIMODAL_MAP.csv"
    )
  )
  
  
  if (
    BIO_FIG_RESUME &&
    !force_requested(
      "FIG5A",
      paste0(
        "FIG5A_",
        branch_tag
      )
    ) &&
    files_complete(
      c(
        fig5a_png_target,
        fig5a_csv_target
      )
    )
  ) {
    
    message(
      "    [SKIP FIG5A] ",
      family,
      " | ",
      program_short,
      " | ",
      branch_tag,
      " | ",
      scope
    )
    
    return(
      invisible(NULL)
    )
  }
  # ---------------------------------------------------------------------------
  # Features únicas del programa.
  # ---------------------------------------------------------------------------
  
  feature_meta <- unique(
    z[
      ,
      c(
        "feature_key",
        "view_code",
        "view_label",
        "feature_original",
        "feature_label",
        "pooled_variance",
        "variance_percentile_within_view"
      ),
      drop = FALSE
    ]
  )
  
  
  # Drivers primero dentro de TX.
  driver_any <- tapply(
    z$is_tx_driver,
    z$feature_key,
    function(x) {
      any(
        x %in%
          TRUE
      )
    }
  )
  
  
  feature_meta$driver_any <-
    as.logical(
      driver_any[
        feature_meta$feature_key
      ]
    )
  
  
  feature_meta$driver_any[
    is.na(
      feature_meta$driver_any
    )
  ] <- FALSE
  
  
  feature_meta <- feature_meta[
    order(
      match(
        feature_meta$view_code,
        c(
          "tx",
          "pr",
          "me",
          "cl"
        )
      ),
      -as.integer(
        feature_meta$driver_any
      ),
      -feature_meta$variance_percentile_within_view,
      feature_meta$feature_label
    ),
    ,
    drop = FALSE
  ]
  
  
  feature_levels <- rev(
    feature_meta$feature_key
  )
  
  
  feature_labels <- setNames(
    fig5_wrap(
      feature_meta$feature_label,
      width = FIG5A_LABEL_WRAP
    ),
    feature_meta$feature_key
  )
  
  
  pmap <- pmap[
    order(
      pmap$library,
      pmap$pathway
    ),
    ,
    drop = FALSE
  ]
  
  
  pathway_levels <- pmap$pathway_key
  
  
  pathway_labels <- setNames(
    paste0(
      fig5_wrap(
        pmap$pathway_label,
        width = FIG5A_LABEL_WRAP
      ),
      "\n[",
      pmap$library,
      "]"
    ),
    pmap$pathway_key
  )
  
  
  grid <- expand.grid(
    
    feature_key =
      feature_levels,
    
    pathway_key =
      pathway_levels,
    
    stringsAsFactors =
      FALSE
  )
  
  
  grid <- merge(
    grid,
    feature_meta[
      ,
      c(
        "feature_key",
        "view_code",
        "view_label",
        "feature_label"
      ),
      drop = FALSE
    ],
    by =
      "feature_key",
    all.x =
      TRUE,
    sort =
      FALSE
  )
  
  
  grid$is_tx_member <-
    FALSE
  
  
  grid$is_tx_driver <-
    FALSE
  
  
  grid$tx_driver_weight <-
    NA_real_
  
  
  grid$non_tx_linked <-
    FALSE
  
  
  grid$rho_spearman <-
    NA_real_
  
  
  # ---------------------------------------------------------------------------
  # TX membership.
  # ---------------------------------------------------------------------------
  
  tx_rel <- z[
    z$view_code ==
      "tx",
    ,
    drop = FALSE
  ]
  
  
  if (
    nrow(tx_rel)
  ) {
    
    key_grid <- paste(
      grid$feature_key,
      grid$pathway_key,
      sep = "\r"
    )
    
    
    key_tx <- paste(
      tx_rel$feature_key,
      tx_rel$pathway_key,
      sep = "\r"
    )
    
    
    idx <- match(
      key_tx,
      key_grid
    )
    
    
    valid <- !is.na(
      idx
    )
    
    
    grid$is_tx_member[
      idx[valid]
    ] <- TRUE
    
    
    grid$is_tx_driver[
      idx[
        valid &
          tx_rel$is_tx_driver
      ]
    ] <- TRUE
    
    
    for (
      i in which(valid)
    ) {
      
      if (
        is.finite(
          tx_rel$tx_driver_weight[i]
        )
      ) {
        
        grid$tx_driver_weight[
          idx[i]
        ] <-
          tx_rel$tx_driver_weight[i]
      }
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # PR / ME / CL correlations.
  # ---------------------------------------------------------------------------
  
  corr <- pathway_correlations[
    pathway_correlations$family ==
      family &
      pathway_correlations$program_id ==
      program_id &
      pathway_correlations$scope ==
      scope,
    ,
    drop = FALSE
  ]
  
  
  if (
    nrow(corr)
  ) {
    
    corr$feature_key <- paste(
      corr$view_code,
      corr$feature_original,
      sep = "\r"
    )
    
    
    key_grid <- paste(
      grid$feature_key,
      grid$pathway_key,
      sep = "\r"
    )
    
    
    key_corr <- paste(
      corr$feature_key,
      corr$pathway_key,
      sep = "\r"
    )
    
    
    idx <- match(
      key_corr,
      key_grid
    )
    
    
    valid <- !is.na(
      idx
    )
    
    
    grid$non_tx_linked[
      idx[valid]
    ] <- TRUE
    
    
    grid$rho_spearman[
      idx[valid]
    ] <-
      corr$rho_spearman[
        valid
      ]
  }
  
  
  grid$feature_plot <- factor(
    grid$feature_key,
    levels =
      feature_levels
  )
  
  
  grid$pathway_plot <- factor(
    grid$pathway_key,
    levels =
      pathway_levels
  )
  
  
  grid$view_label <- factor(
    grid$view_label,
    levels = c(
      "TX",
      "PR",
      "ME",
      "CL"
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Figura.
  #
  # TX:
  #   gris = membership.
  #   punto negro = driver.
  #   tamaño punto = normalized_anchor_weight.
  #
  # PR/ME/CL:
  #   color = Spearman pathway-score ↔ feature.
  # ---------------------------------------------------------------------------
  
  p <- ggplot2::ggplot(
    grid,
    ggplot2::aes(
      x = pathway_plot,
      y = feature_plot
    )
  ) +
    
    ggplot2::geom_tile(
      fill = "white",
      color = "grey92",
      linewidth = 0.25
    ) +
    
    ggplot2::geom_tile(
      data = grid[
        grid$is_tx_member,
        ,
        drop = FALSE
      ],
      fill = "grey82",
      color = "grey70",
      linewidth = 0.25
    ) +
    
    ggplot2::geom_tile(
      data = grid[
        grid$non_tx_linked &
          !is.finite(
            grid$rho_spearman
          ),
        ,
        drop = FALSE
      ],
      fill = "grey94",
      color = "grey85",
      linewidth = 0.25
    )
  
  
  weighted <- grid[
    grid$non_tx_linked &
      is.finite(
        grid$rho_spearman
      ),
    ,
    drop = FALSE
  ]
  
  
  if (
    nrow(weighted)
  ) {
    
    p <- p +
      
      ggplot2::geom_tile(
        data =
          weighted,
        
        ggplot2::aes(
          fill =
            rho_spearman
        ),
        
        color =
          "grey78",
        
        linewidth =
          0.25
      ) +
      
      ggplot2::scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = c(
          -1,
          1
        ),
        name = "Spearman\nrho"
      )
  }
  
  
  drivers <- grid[
    grid$is_tx_driver &
      is.finite(
        grid$tx_driver_weight
      ),
    ,
    drop = FALSE
  ]
  
  
  if (
    nrow(drivers)
  ) {
    
    p <- p +
      
      ggplot2::geom_point(
        data =
          drivers,
        
        ggplot2::aes(
          size =
            tx_driver_weight
        ),
        
        shape = 21,
        fill = "black",
        color = "black",
        inherit.aes = TRUE
      ) +
      
      ggplot2::scale_size_continuous(
        name =
          "TX normalized\nanchor weight",
        range = c(
          1.5,
          5
        )
      )
  }
  
  
  n_scope <- if (
    scope ==
    "GLOBAL"
  ) {
    
    length(
      unique(
        pathway_scores_global_ids <- pathway_correlations$feature_original[0]
      )
    )
    
  } else {
    
    NA_integer_
  }
  
  
  p <- p +
    
    ggplot2::scale_x_discrete(
      labels =
        pathway_labels
    ) +
    
    ggplot2::scale_y_discrete(
      labels =
        feature_labels
    ) +
    
    ggplot2::facet_grid(
      view_label ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    
    ggplot2::labs(
      
      title = paste0(
        family,
        " — ",
        phenotype_label(
          phenotypic_side
        ),
        " — ",
        program_short
      ),
      
      subtitle = paste0(
        "Arquitectura multimodal descriptiva — ",
        scope,
        ". TX = membership/driver; PR/ME/CL = Spearman con el score congelado del pathway."
      ),
      
      x =
        "Pathway TX primario",
      
      y =
        NULL,
      
      caption = paste0(
        "Los 26 individuos se analizan conjuntamente en GLOBAL; ",
        "NP, FD y SO son análisis descriptivos intra-fenotipo. ",
        "No se utiliza ninguna correlación de esta figura para redefinir ",
        "pathways, programas o features."
      )
    ) +
    
    fig_theme(
      FIG5A_BASE_SIZE
    )+    
    ggplot2::theme(
      panel.grid =
        ggplot2::element_blank(),
      
      axis.text.x =
        ggplot2::element_text(
          angle = 55,
          hjust = 1,
          vjust = 1,
          size = FIG5A_AXIS_TEXT_SIZE
        ),
      
      axis.text.y =
        ggplot2::element_text(
          size = FIG5A_AXIS_TEXT_SIZE
        ),
      strip.text.y =
        ggplot2::element_text(
          face = "bold",
          angle = 0,
          size = FIG5A_STRIP_TEXT_SIZE
        ),
      
      legend.position =
        "bottom"
    )
  
  
  dir_fig <- file.path(
    outdir,
    "figures",
    "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
    branch,
    family,
    safe_filename(
      factor_side
    ),
    safe_filename(
      program_short
    ),
    "01_PROGRAM_PATHWAY_MAP",
    scope
  )
  
  
  dir_dat <- file.path(
    outdir,
    "plot_data",
    "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
    branch,
    family,
    safe_filename(
      factor_side
    ),
    safe_filename(
      program_short
    ),
    "01_PROGRAM_PATHWAY_MAP",
    scope
  )
  
  
  dir.create(
    dir_fig,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  dir.create(
    dir_dat,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  ggplot2::ggsave(
    
    filename = file.path(
      dir_fig,
      paste0(
        "FIG5A_",
        scope,
        "_PROGRAM_PATHWAY_MULTIMODAL_MAP.png"
      )
    ),
    
    plot = p,
    
    width = max(
      9,
      min(
        24,
        5 +
          1.0 *
          length(
            pathway_levels
          )
      )
    ),
    
    height = max(
      8,
      min(
        40,
        5 +
          FIG5A_INCH_PER_FEATURE *
          length(
            feature_levels
          )
      )
    ),
    dpi =
      FIG_DPI,
    
    limitsize =
      FALSE
  )
  
  
  write_csv(
    grid,
    file.path(
      dir_dat,
      paste0(
        "FIG5A_",
        scope,
        "_PROGRAM_PATHWAY_MULTIMODAL_MAP.csv"
      )
    )
  )
  
  
  invisible(p)
}


# ===== 8.10 EXTRAER MATRIZ DE FEATURES DE UNA VÍA =====

build_fig5_pathway_feature_matrix <- function(
    catalog,
    matrices,
    sample_meta,
    family,
    pathway_key,
    selected_feature_keys = NULL
) {
  
  z <- catalog[
    catalog$family ==
      family &
      catalog$pathway_key ==
      pathway_key &
      catalog$available_in_pooled,
    ,
    drop = FALSE
  ]
  
  
  if (
    !is.null(
      selected_feature_keys
    )
  ) {
    
    z <- z[
      z$feature_key %in%
        selected_feature_keys,
      ,
      drop = FALSE
    ]
  }
  
  
  if (!nrow(z)) {
    
    return(NULL)
  }
  
  
  feature_meta <- unique(
    z[
      ,
      c(
        "feature_key",
        "view_code",
        "view_label",
        "feature_original",
        "feature_label",
        "is_tx_driver",
        "pooled_variance",
        "variance_percentile_within_view"
      ),
      drop = FALSE
    ]
  )
  
  
  driver_any <- tapply(
    z$is_tx_driver,
    z$feature_key,
    function(x) {
      any(
        x %in%
          TRUE
      )
    }
  )
  
  
  feature_meta$driver_any <-
    as.logical(
      driver_any[
        feature_meta$feature_key
      ]
    )
  
  
  feature_meta$driver_any[
    is.na(
      feature_meta$driver_any
    )
  ] <- FALSE
  
  
  # ---------------------------------------------------------------------------
  # Orden:
  #
  # TX -> PR -> ME -> CL.
  #
  # Dentro de TX:
  #   drivers primero.
  #
  # Dentro de cada vista:
  #   mayor variabilidad primero, solo como orden gráfico.
  # ---------------------------------------------------------------------------
  
  feature_meta <- feature_meta[
    order(
      match(
        feature_meta$view_code,
        c(
          "tx",
          "pr",
          "me",
          "cl"
        )
      ),
      -as.integer(
        feature_meta$driver_any
      ),
      -feature_meta$variance_percentile_within_view,
      feature_meta$feature_label
    ),
    ,
    drop = FALSE
  ]
  
  
  sample_ids <- sample_meta$sample_id
  
  
  M <- matrix(
    NA_real_,
    
    nrow =
      nrow(
        feature_meta
      ),
    
    ncol =
      length(
        sample_ids
      ),
    
    dimnames =
      list(
        feature_meta$feature_key,
        sample_ids
      )
  )
  
  
  for (
    i in seq_len(
      nrow(
        feature_meta
      )
    )
  ) {
    
    v <-
      feature_meta$view_code[i]
    
    f <-
      feature_meta$feature_original[i]
    
    
    M[i, ] <- as.numeric(
      matrices[[v]][
        f,
        sample_ids
      ]
    )
  }
  
  
  list(
    matrix =
      M,
    
    feature_meta =
      feature_meta,
    
    pathway_meta =
      z[1L, , drop = FALSE]
  )
}

# ===== 8.10B ORDEN GLOBAL POR CLUSTERING DENTRO DE CADA VISTA =====

get_fig5_clustered_feature_order <- function(
    pathway_matrix,
    sample_meta
) {
  
  M <- pathway_matrix$matrix
  
  feature_meta <- pathway_matrix$feature_meta
  
  assert(
    identical(
      rownames(M),
      feature_meta$feature_key
    ),
    "Feature metadata and matrix rows are not aligned"
  )
  if (
    !isTRUE(
      FIG5B_CLUSTER_WITHIN_VIEW
    )
  ) {
    
    return(
      feature_meta$feature_key
    )
  }
  
  # ---------------------------------------------------------------------------
  # El clustering SIEMPRE se define usando GLOBAL.
  #
  # Después este mismo orden se reutiliza para GLOBAL / NP / FD / SO.
  #
  # No se mezclan vistas:
  #   TX se clusteriza con TX
  #   PR con PR
  #   ME con ME
  #   CL con CL
  #
  # El dendrograma NO se representa.
  # ---------------------------------------------------------------------------
  
  view_order <- c(
    "tx",
    "pr",
    "me",
    "cl"
  )
  
  
  ordered_keys <- character()
  
  
  for (v in view_order) {
    
    idx <- which(
      feature_meta$view_code == v
    )
    
    
    if (!length(idx)) {
      next
    }
    
    
    keys_v <- feature_meta$feature_key[idx]
    
    
    # Una sola feature: no hay nada que clusterizar.
    if (length(keys_v) == 1L) {
      
      ordered_keys <- c(
        ordered_keys,
        keys_v
      )
      
      next
    }
    
    
    X <- M[
      keys_v,
      ,
      drop = FALSE
    ]
    
    
    # Eliminar únicamente features sin variación suficiente
    # PARA EL CÁLCULO DEL ORDEN.
    # No se eliminan de la figura.
    usable <- apply(
      X,
      1,
      function(z) {
        
        z <- z[
          is.finite(z)
        ]
        
        length(z) >= 3L &&
          is.finite(
            stats::sd(z)
          ) &&
          stats::sd(z) > 0
      }
    )
    
    
    usable_keys <- keys_v[
      usable
    ]
    
    
    nonusable_keys <- keys_v[
      !usable
    ]
    
    
    if (length(usable_keys) >= 2L) {
      
      R <- suppressWarnings(
        stats::cor(
          t(
            X[
              usable_keys,
              ,
              drop = FALSE
            ]
          ),
          use = "pairwise.complete.obs",
          method = "spearman"
        )
      )
      
      
      # -----------------------------------------------------------------------
      # Distancia:
      #
      #   d = 1 - rho
      #
      # rho = +1 -> distancia 0
      # rho =  0 -> distancia 1
      # rho = -1 -> distancia 2
      #
      # Por tanto agrupa patrones que covarían en la MISMA dirección.
      # -----------------------------------------------------------------------
      
      R[
        !is.finite(R)
      ] <- 0
      
      
      diag(R) <- 1
      
      
      D <- stats::as.dist(
        1 - R
      )
      
      
      hc <- stats::hclust(
        D,
        method = FIG5B_CLUSTER_METHOD
      )
      
      
      clustered_keys <- usable_keys[
        hc$order
      ]
      
    } else {
      
      clustered_keys <- usable_keys
    }
    
    
    ordered_keys <- c(
      ordered_keys,
      clustered_keys,
      nonusable_keys
    )
  }
  
  
  # Control de integridad.
  assert(
    setequal(
      ordered_keys,
      rownames(M)
    ),
    "Clustering order lost or duplicated FIG5B features"
  )
  
  
  ordered_keys
}
# ===== 8.11 CALCULAR MATRIZ SPEARMAN FEATURE × FEATURE =====

compute_fig5_feature_correlation_matrix <- function(
    pathway_matrix,
    sample_meta,
    scope = c(
      "GLOBAL",
      "NP",
      "FD",
      "SO"
    )
) {
  
  scope <- match.arg(
    scope
  )
  
  
  M <- pathway_matrix$matrix
  
  
  ids <- colnames(
    M
  )
  
  
  meta <- sample_meta[
    match(
      ids,
      sample_meta$sample_id
    ),
    ,
    drop = FALSE
  ]
  
  
  if (
    scope !=
    "GLOBAL"
  ) {
    
    keep <- meta$grupo ==
      scope
    
    
    M <- M[
      ,
      keep,
      drop = FALSE
    ]
    
    
    meta <- meta[
      keep,
      ,
      drop = FALSE
    ]
  }
  
  
  n_samples <- ncol(
    M
  )
  
  
  if (
    n_samples <
    3L
  ) {
    
    return(
      list(
        rho =
          matrix(
            NA_real_,
            nrow = nrow(M),
            ncol = nrow(M),
            dimnames = list(
              rownames(M),
              rownames(M)
            )
          ),
        
        n =
          matrix(
            n_samples,
            nrow = nrow(M),
            ncol = nrow(M),
            dimnames = list(
              rownames(M),
              rownames(M)
            )
          ),
        
        scope =
          scope,
        
        n_samples =
          n_samples
      )
    )
  }
  
  
  rho <- suppressWarnings(
    stats::cor(
      t(M),
      use =
        "pairwise.complete.obs",
      method =
        "spearman"
    )
  )
  
  
  valid <- is.finite(
    M
  )
  
  
  n_pair <- tcrossprod(
    valid
  )
  
  
  dimnames(
    n_pair
  ) <- list(
    rownames(M),
    rownames(M)
  )
  
  
  list(
    rho =
      rho,
    
    n =
      n_pair,
    
    scope =
      scope,
    
    n_samples =
      n_samples
  )
}


# ===== 8.12 FIGURA 5B — MATRIZ SPEARMAN MULTIMODAL POR PATHWAY =====

# ===== 8.12 FIGURA 5B — MATRIZ SPEARMAN MULTIMODAL POR PATHWAY =====

make_fig5_pathway_correlation_heatmap <- function(
    correlation_object,
    pathway_matrix,
    family,
    outdir,
    branch,
    feature_order = NULL,
    signal_type = "OBSERVED"
) {
  
  R <- correlation_object$rho
  
  N <- correlation_object$n
  
  scope <- correlation_object$scope
  
  feature_meta <-
    pathway_matrix$feature_meta
  
  
  pmeta <-
    pathway_matrix$pathway_meta
  
  
  if (
    !nrow(R)
  ) {
    
    return(
      invisible(NULL)
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Orden.
  #
  # Si recibimos feature_order:
  #   usa exactamente el orden GLOBAL previamente calculado.
  #
  # En caso contrario:
  #   conserva el orden original.
  # ---------------------------------------------------------------------------
  
  if (
    is.null(
      feature_order
    )
  ) {
    
    feature_order <-
      feature_meta$feature_key
  }
  
  
  feature_order <- feature_order[
    feature_order %in%
      rownames(R)
  ]
  
  
  assert(
    length(feature_order) ==
      nrow(R),
    paste(
      "FIG5B feature-order mismatch:",
      family,
      pmeta$pathway[1L],
      scope,
      signal_type
    )
  )
  
  
  R <- R[
    feature_order,
    feature_order,
    drop = FALSE
  ]
  
  
  N <- N[
    feature_order,
    feature_order,
    drop = FALSE
  ]
  
  
  feature_meta <- feature_meta[
    match(
      feature_order,
      feature_meta$feature_key
    ),
    ,
    drop = FALSE
  ]
  
  
  n_features <-
    length(
      feature_order
    )
  
  
  grid <- expand.grid(
    
    row_i =
      seq_len(
        n_features
      ),
    
    col_i =
      seq_len(
        n_features
      ),
    
    stringsAsFactors =
      FALSE
  )
  
  
  grid$rho <- R[
    cbind(
      grid$row_i,
      grid$col_i
    )
  ]
  
  
  grid$n_pair <- N[
    cbind(
      grid$row_i,
      grid$col_i
    )
  ]
  
  labels <- fig5_wrap(
    feature_meta$feature_label,
    width = 16L
  )
  
  
  # ---------------------------------------------------------------------------
  # Separadores TX / PR / ME / CL.
  #
  # El clustering ocurre DENTRO de cada vista, de modo que los cuatro
  # bloques siguen siendo contiguos.
  # ---------------------------------------------------------------------------
  
  block_counts <- table(
    factor(
      feature_meta$view_code,
      levels = c(
        "tx",
        "pr",
        "me",
        "cl"
      )
    )
  )
  
  
  block_counts <- as.integer(
    block_counts
  )
  
  
  boundaries <- cumsum(
    block_counts
  )
  
  
  boundaries <- boundaries[
    boundaries > 0 &
      boundaries <
      n_features
  ] + 0.5
  
  
  p <- ggplot2::ggplot(
    grid,
    ggplot2::aes(
      x = col_i,
      y = row_i,
      fill = rho
    )
  ) +
    
    ggplot2::geom_tile(
      color = "grey92",
      linewidth = 0.20
    ) +
    
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(
        -1,
        1
      ),
      na.value = "grey90",
      name = "Spearman\nrho"
    ) +
    
    ggplot2::scale_x_continuous(
      breaks =
        seq_len(
          n_features
        ),
      labels =
        labels,
      position =
        "top"
    ) +
    
    ggplot2::scale_y_continuous(
      breaks =
        seq_len(
          n_features
        ),
      labels =
        labels,
      trans =
        "reverse"
    )
  
  
  if (
    length(
      boundaries
    )
  ) {
    
    p <- p +
      
      ggplot2::geom_vline(
        xintercept =
          boundaries,
        linewidth =
          0.9
      ) +
      
      ggplot2::geom_hline(
        yintercept =
          boundaries,
        linewidth =
          0.9
      )
  }
  
  
  # ---------------------------------------------------------------------------
  # Valores numéricos:
  #
  # PRESENTATION -> siempre.
  # ALL          -> únicamente si la matriz ya es pequeña.
  # ---------------------------------------------------------------------------
  
  if (
    branch ==
    "02_PRESENTATION" ||
    n_features <=
    FIG5_PRESENTATION_MAX_FEATURES
  ) {
    
    grid$label <- ifelse(
      is.finite(
        grid$rho
      ),
      sprintf(
        "%.2f",
        grid$rho
      ),
      ""
    )
    
    
    p <- p +
      
      ggplot2::geom_text(
        data = grid,
        
        ggplot2::aes(
          label = label
        ),
        
        size =
          if (
            n_features <= 15L
          ) {
            FIG5B_CELL_TEXT_SIZE + 0.5
          } else {
            FIG5B_CELL_TEXT_SIZE
          }
      )
  }
  
  
  signal_text <- if (
    signal_type ==
    "OBSERVED"
  ) {
    
    "features observadas/preprocesadas"
    
  } else {
    
    "componente latente centrada explicada por los factores auditados"
  }
  
  
  p <- p +
    
    ggplot2::labs(
      
      title = paste0(
        family,
        " — ",
        phenotype_label(
          pmeta$phenotypic_side[1L]
        ),
        " — ",
        clean_pathway_label(
          pmeta$pathway[1L]
        )
      ),
      
      subtitle = paste0(
        signal_type,
        " — covariación multimodal descriptiva — ",
        scope,
        " (n=",
        correlation_object$n_samples,
        "). ",
        "Filas/columnas = TX | PR | ME | CL."
      ),
      
      x = NULL,
      y = NULL,
      
      caption = paste0(
        "Cada celda es rho de Spearman entre ",
        signal_text,
        ". Orden = bloques ómicos con clustering jerárquico ",
        "interno definido en GLOBAL; el dendrograma no se muestra. ",
        "El clustering es exclusivamente visual y no redefine ",
        "pathways, programas ni features."
      )
    ) +
    
    fig_theme(
      FIG5B_BASE_SIZE
    ) +
    
    ggplot2::theme(
      
      panel.grid =
        ggplot2::element_blank(),
      
      axis.text.x =
        ggplot2::element_text(
          angle = 60,
          hjust = 0,
          vjust = 0,
          size = FIG5B_AXIS_TEXT_SIZE
        ),
      
      axis.text.y =
        ggplot2::element_text(
          size = FIG5B_AXIS_TEXT_SIZE
        ),
      
      plot.title =
        ggplot2::element_text(
          face = "bold",
          size = FIG5B_BASE_SIZE + 3
        ),
      
      plot.subtitle =
        ggplot2::element_text(
          size = FIG5B_BASE_SIZE
        ),
      
      legend.position =
        "bottom"
    )
  
  
  factor_side <-
    pmeta$factor_side[1L]
  
  
  program_short <-
    pmeta$program_short[1L]
  
  
  pathway_folder <- paste0(
    safe_filename(
      pmeta$library[1L]
    ),
    "__",
    safe_filename(
      clean_pathway_label(
        pmeta$pathway[1L]
      )
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # NUEVO:
  #
  # .../02_PATHWAY_CORRELATIONS/pathway/01_OBSERVED/GLOBAL
  # .../02_PATHWAY_CORRELATIONS/pathway/02_LATENT/GLOBAL
  # ---------------------------------------------------------------------------
  
  signal_folder <- if (
    signal_type ==
    "OBSERVED"
  ) {
    "01_OBSERVED"
  } else {
    "02_LATENT"
  }
  
  
  dir_fig <- file.path(
    outdir,
    "figures",
    "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
    branch,
    family,
    safe_filename(
      factor_side
    ),
    safe_filename(
      program_short
    ),
    "02_PATHWAY_CORRELATIONS",
    pathway_folder,
    signal_folder,
    scope
  )
  
  
  dir_dat <- file.path(
    outdir,
    "plot_data",
    "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
    branch,
    family,
    safe_filename(
      factor_side
    ),
    safe_filename(
      program_short
    ),
    "02_PATHWAY_CORRELATIONS",
    pathway_folder,
    signal_folder,
    scope
  )
  
  
  dir.create(
    dir_fig,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  dir.create(
    dir_dat,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  file_stub <- paste0(
    "FIG5B_",
    signal_type,
    "_",
    scope,
    "_FEATURE_FEATURE_SPEARMAN"
  )
  
  
  # ---------------------------------------------------------------------------
  # PRESENTATION: más espacio por feature para que las etiquetas sean legibles.
  #
  # ALL sigue permitiendo matrices grandes, pero no intenta hacer legibles
  # cientos de etiquetas a tamaño de presentación.
  # ---------------------------------------------------------------------------
  fig_size <- if (
    branch ==
    "02_PRESENTATION"
  ) {
    
    max(
      11,
      min(
        18,
        6 +
          0.45 *
          n_features
      )
    )
    
  } else {
    
    max(
      10,
      min(
        45,
        5 +
          0.42 *
          n_features
      )
    )
  }
  
  
  ggplot2::ggsave(
    
    filename = file.path(
      dir_fig,
      paste0(
        file_stub,
        ".png"
      )
    ),
    
    plot = p,
    
    width =
      fig_size,
    
    height =
      fig_size,
    
    dpi =
      FIG_DPI,
    
    limitsize =
      FALSE
  )
  
  
  rho_df <- data.frame(
    feature =
      rownames(R),
    R,
    check.names =
      FALSE,
    stringsAsFactors =
      FALSE
  )
  
  
  n_df <- data.frame(
    feature =
      rownames(N),
    N,
    check.names =
      FALSE,
    stringsAsFactors =
      FALSE
  )
  
  
  write_csv(
    rho_df,
    file.path(
      dir_dat,
      paste0(
        file_stub,
        "_RHO_MATRIX.csv"
      )
    )
  )
  
  
  write_csv(
    n_df,
    file.path(
      dir_dat,
      paste0(
        file_stub,
        "_N_MATRIX.csv"
      )
    )
  )
  
  
  feature_export <- feature_meta
  
  feature_export$scope <-
    scope
  
  feature_export$signal_type <-
    signal_type
  
  feature_export$n_scope <-
    correlation_object$n_samples
  
  feature_export$display_order <-
    seq_len(
      nrow(
        feature_export
      )
    )
  
  
  write_csv(
    feature_export,
    file.path(
      dir_dat,
      paste0(
        file_stub,
        "_FEATURES_AND_ORDER.csv"
      )
    )
  )
  
  
  invisible(p)
}


fig5b_expected_files <- function(
    pathway_matrix,
    family,
    outdir,
    branch,
    scope,
    signal_type
) {
  
  pmeta <-
    pathway_matrix$pathway_meta
  
  factor_side <-
    pmeta$factor_side[1L]
  
  program_short <-
    pmeta$program_short[1L]
  
  pathway_folder <- paste0(
    safe_filename(
      pmeta$library[1L]
    ),
    "__",
    safe_filename(
      clean_pathway_label(
        pmeta$pathway[1L]
      )
    )
  )
  
  signal_folder <- if (
    signal_type == "OBSERVED"
  ) {
    "01_OBSERVED"
  } else {
    "02_LATENT"
  }
  
  file_stub <- paste0(
    "FIG5B_",
    signal_type,
    "_",
    scope,
    "_FEATURE_FEATURE_SPEARMAN"
  )
  
  dir_fig <- file.path(
    outdir,
    "figures",
    "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
    branch,
    family,
    safe_filename(factor_side),
    safe_filename(program_short),
    "02_PATHWAY_CORRELATIONS",
    pathway_folder,
    signal_folder,
    scope
  )
  
  dir_dat <- file.path(
    outdir,
    "plot_data",
    "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
    branch,
    family,
    safe_filename(factor_side),
    safe_filename(program_short),
    "02_PATHWAY_CORRELATIONS",
    pathway_folder,
    signal_folder,
    scope
  )
  
  c(
    png =
      file.path(
        dir_fig,
        paste0(
          file_stub,
          ".png"
        )
      ),
    
    rho =
      file.path(
        dir_dat,
        paste0(
          file_stub,
          "_RHO_MATRIX.csv"
        )
      ),
    
    n =
      file.path(
        dir_dat,
        paste0(
          file_stub,
          "_N_MATRIX.csv"
        )
      ),
    
    features =
      file.path(
        dir_dat,
        paste0(
          file_stub,
          "_FEATURES_AND_ORDER.csv"
        )
      )
  )
}


fig5b_should_run <- function(
    pathway_matrix,
    family,
    outdir,
    branch,
    scope,
    signal_type
) {
  
  branch_tag <- if (
    branch == "02_PRESENTATION"
  ) {
    "PRESENTATION"
  } else {
    "ALL"
  }
  
  
  forced <- force_requested(
    "FIG5B",
    paste0(
      "FIG5B_",
      branch_tag
    ),
    paste0(
      "FIG5B_",
      branch_tag,
      "_",
      signal_type
    )
  )
  
  
  if (
    !BIO_FIG_RESUME ||
    forced
  ) {
    return(TRUE)
  }
  
  
  !files_complete(
    fig5b_expected_files(
      pathway_matrix =
        pathway_matrix,
      family =
        family,
      outdir =
        outdir,
      branch =
        branch,
      scope =
        scope,
      signal_type =
        signal_type
    )
  )
}


# ===== 8.13 EJECUTAR PREGUNTA 5 =====
# ---------------------------------------------------------------------------
# Si Q1-Q4 fueron omitidas durante RESUME, Q5 solo necesita reconstruir
# las entradas congeladas de Q4. No se regeneran las figuras de Q4.
# ---------------------------------------------------------------------------
if (run_question(5L)) {
if (
  !exists("fig4_inputs") ||
  is.null(fig4_inputs)
) {
  
  fig4_inputs <- list()
  
  for (family in FAMILIES) {
    
    message(
      "[RESUME] Rebuilding frozen Q4 input for ",
      family
    )
    
    fig4_inputs[[family]] <-
      load_fig4_inputs(
        bio_program_run =
          BIO_RUN,
        family =
          family
      )
  }
}
message("")
message("============================================================")
message("METATEST — QUESTION 5")
message("Pathway-centered multimodal architecture | pooled descriptive")
message("============================================================")


fig5_inputs <- list()

fig5_catalog_list <- list()

fig5_pathway_scores_list <- list()

fig5_correlations_list <- list()

fig5_presentation_audit_list <- list()


for (family in FAMILIES) {
  
  message(
    "\n[FIG5] ",
    family
  )
  
  
  inp5 <- load_fig5_inputs(
    bio_program_run =
      BIO_RUN,
    
    family =
      family
  )
  
  
  fig5_inputs[[family]] <-
    inp5
  
  
  # ---------------------------------------------------------------------------
  # Scores congelados aplicados a los 26.
  # ---------------------------------------------------------------------------
  
  pathway_scores <-
    build_fig5_pooled_pathway_scores(
      inp5
    )
  
  
  fig5_pathway_scores_list[[family]] <-
    pathway_scores
  
  
  # ---------------------------------------------------------------------------
  # Catálogo TX + PR + ME + CL.
  #
  # Reutilizamos Q4 para TX.
  # ---------------------------------------------------------------------------
  
  catalog <-
    build_fig5_feature_catalog(
      
      inputs =
        inp5,
      
      tx_pathway_gene =
        fig4_inputs[[family]]$pathway_gene
    )
  
  
  catalog <-
    add_fig5_pooled_variance(
      catalog =
        catalog,
      
      matrices =
        inp5$matrices
    )
  
  
  fig5_catalog_list[[family]] <-
    catalog
  
  
  # ---------------------------------------------------------------------------
  # Correlaciones pathway ↔ feature descriptivas.
  # ---------------------------------------------------------------------------
  
  corr_scope <- list()
  
  
  for (
    scope in c(
      "GLOBAL",
      "NP",
      "FD",
      "SO"
    )
  ) {
    
    corr_scope[[scope]] <-
      compute_fig5_pathway_feature_correlations(
        
        catalog =
          catalog,
        
        pathway_scores =
          pathway_scores,
        
        matrices =
          inp5$matrices,
        
        sample_meta =
          inp5$sample_meta,
        
        scope =
          scope
      )
  }
  
  
  all_corr <- do.call(
    rbind,
    corr_scope
  )
  
  
  fig5_correlations_list[[family]] <-
    all_corr
  
  
  # =======================================================================
  # 01_ALL
  # =======================================================================
  
  programs_i <- unique(
    inp5$programs$program_id
  )
  
  
  for (
    pid in programs_i
  ) {
    
    message(
      "  [ALL] ",
      pid
    )
    
    
    z_program <- catalog[
      catalog$program_id ==
        pid &
        catalog$available_in_pooled,
      ,
      drop = FALSE
    ]
    
    
    if (!nrow(z_program)) {
      next
    }
    
    
    # -----------------------------------------------------------------------
    # 5A — programa completo.
    # -----------------------------------------------------------------------
    
    for (
      scope in c(
        "GLOBAL",
        "NP",
        "FD",
        "SO"
      )
    ) {
      
      make_fig5_program_architecture_map(
        
        catalog =
          z_program,
        
        pathway_correlations =
          all_corr,
        
        memberships =
          inp5$memberships,
        
        family =
          family,
        
        program_id =
          pid,
        
        scope =
          scope,
        
        outdir =
          OUT_ROOT,
        
        branch =
          "01_ALL"
      )
    }
    
    
    # -----------------------------------------------------------------------
    # 5B — una matriz de correlación por pathway.
    # -----------------------------------------------------------------------
    
    pathways_i <- unique(
      z_program$pathway_key
    )
    
    for (
      pk in pathways_i
    ) {
      
      # =======================================================================
      # OBSERVED
      # =======================================================================
      
      path_matrix_observed <-
        build_fig5_pathway_feature_matrix(
          
          catalog =
            catalog,
          
          matrices =
            inp5$matrices,
          
          sample_meta =
            inp5$sample_meta,
          
          family =
            family,
          
          pathway_key =
            pk,
          
          selected_feature_keys =
            NULL
        )
      
      
      if (
        is.null(
          path_matrix_observed
        )
      ) {
        next
      }
      
      
      feature_order_observed <-
        get_fig5_clustered_feature_order(
          
          pathway_matrix =
            path_matrix_observed,
          
          sample_meta =
            inp5$sample_meta
        )
      
      
      # =======================================================================
      # LATENT
      #
      # Exactamente las MISMAS features.
      # Lo único que cambia son sus valores:
      #
      # OBSERVED = datos preprocesados
      # LATENT   = L %*% Z_oriented_scaled
      # =======================================================================
      
      path_matrix_latent <-
        build_fig5_pathway_feature_matrix(
          
          catalog =
            catalog,
          
          matrices =
            inp5$latent_matrices,
          
          sample_meta =
            inp5$sample_meta,
          
          family =
            family,
          
          pathway_key =
            pk,
          
          selected_feature_keys =
            NULL
        )
      
      
      assert(
        !is.null(
          path_matrix_latent
        ),
        paste(
          "Could not construct latent pathway matrix:",
          family,
          pk
        )
      )
      
      
      assert(
        setequal(
          path_matrix_observed$feature_meta$feature_key,
          path_matrix_latent$feature_meta$feature_key
        ),
        paste(
          "OBSERVED/LATENT feature mismatch:",
          family,
          pk
        )
      )
      
      
      feature_order_latent <-
        get_fig5_clustered_feature_order(
          
          pathway_matrix =
            path_matrix_latent,
          
          sample_meta =
            inp5$sample_meta
        )
      
      
      # =======================================================================
      # GLOBAL / NP / FD / SO
      # =======================================================================
      
      for (
        scope in c(
          "GLOBAL",
          "NP",
          "FD",
          "SO"
        )
      ) {
        
        # ---------------------------------------------------------------------
        # OBSERVED
        # ---------------------------------------------------------------------
        
        if (
          fig5b_should_run(
            pathway_matrix =
            path_matrix_observed,
            family =
            family,
            outdir =
            OUT_ROOT,
            branch =
            "01_ALL",
            scope =
            scope,
            signal_type =
            "OBSERVED"
          )
        ) {
          
          cor_observed <-
            compute_fig5_feature_correlation_matrix(
              
              pathway_matrix =
                path_matrix_observed,
              
              sample_meta =
                inp5$sample_meta,
              
              scope =
                scope
            )
          
          
          make_fig5_pathway_correlation_heatmap(
            
            correlation_object =
              cor_observed,
            
            pathway_matrix =
              path_matrix_observed,
            
            family =
              family,
            
            outdir =
              OUT_ROOT,
            
            branch =
              "01_ALL",
            
            feature_order =
              feature_order_observed,
            
            signal_type =
              "OBSERVED"
          )
          
        } else {
          
          message(
            "    [SKIP FIG5B] ALL | OBSERVED | ",
            scope
          )
        }
        
        
        # ---------------------------------------------------------------------
        # LATENT
        # ---------------------------------------------------------------------
        
        if (
          fig5b_should_run(
            pathway_matrix =
            path_matrix_latent,
            family =
            family,
            outdir =
            OUT_ROOT,
            branch =
            "01_ALL",
            scope =
            scope,
            signal_type =
            "LATENT"
          )
        ) {
          
          cor_latent <-
            compute_fig5_feature_correlation_matrix(
              
              pathway_matrix =
                path_matrix_latent,
              
              sample_meta =
                inp5$sample_meta,
              
              scope =
                scope
            )
          
          
          make_fig5_pathway_correlation_heatmap(
            
            correlation_object =
              cor_latent,
            
            pathway_matrix =
              path_matrix_latent,
            
            family =
              family,
            
            outdir =
              OUT_ROOT,
            
            branch =
              "01_ALL",
            
            feature_order =
              feature_order_latent,
            
            signal_type =
              "LATENT"
          )
          
        } else {
          
          message(
            "    [SKIP FIG5B] ALL | LATENT | ",
            scope
          )
        }
      }
    }
    
    # =====================================================================
    # 02_PRESENTATION
    # =====================================================================
    
    presentation_selection <-
      select_fig5_presentation_features(
        
        catalog =
          z_program,
        
        max_features =
          FIG5_PRESENTATION_MAX_FEATURES
      )
    
    
    presentation_selection$family <-
      family
    
    
    presentation_selection$program_id <-
      pid
    
    
    fig5_presentation_audit_list[[paste(family,pid,sep = "__")]] <-
      presentation_selection
    
    
    selected_program_keys <-
      presentation_selection$feature_key[
        presentation_selection$selected_for_presentation
      ]
    
    
    z_program_presentation <- z_program[
      z_program$feature_key %in%
        selected_program_keys,
      ,
      drop = FALSE
    ]
    
    
    message(
      "  [PRESENTATION] ",
      pid,
      " | ",
      length(
        unique(
          z_program$feature_key
        )
      ),
      " -> ",
      length(
        unique(
          z_program_presentation$feature_key
        )
      ),
      " features"
    )
    
    
    # -----------------------------------------------------------------------
    # 5A presentation.
    # -----------------------------------------------------------------------
    
    for (
      scope in c(
        "GLOBAL",
        "NP",
        "FD",
        "SO"
      )
    ) {
      
      make_fig5_program_architecture_map(
        
        catalog =
          z_program_presentation,
        
        pathway_correlations =
          all_corr,
        
        memberships =
          inp5$memberships,
        
        family =
          family,
        
        program_id =
          pid,
        
        scope =
          scope,
        
        outdir =
          OUT_ROOT,
        
        branch =
          "02_PRESENTATION"
      )
    }
    
    
    # -----------------------------------------------------------------------
    # 5B presentation:
    #
    # IMPORTANTE:
    # selección independiente POR PATHWAY.
    #
    # Si la vía ya tiene <= MAX, conserva TODO.
    # -----------------------------------------------------------------------
    
    for (
      pk in pathways_i
    ) {
      
      z_path <- catalog[
        catalog$family ==
          family &
          catalog$pathway_key ==
          pk &
          catalog$available_in_pooled,
        ,
        drop = FALSE
      ]
      
      
      if (!nrow(z_path)) {
        next
      }
      
      
      path_selection <-
        select_fig5_presentation_features(
          
          catalog =
            z_path,
          
          max_features =
            FIG5_PRESENTATION_MAX_FEATURES
        )
      
      
      selected_path_keys <-
        path_selection$feature_key[
          path_selection$selected_for_presentation
        ]
      
      message("")
      message("QUESTION 5 COMPLETE")
      
      message(
        "Output: ",
        file.path(
          OUT_ROOT,
          "figures",
          "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE"
        )
      )
      
    }
  }
}


# ===== 8.14 TABLAS MAESTRAS =====

all_fig5_catalog <- do.call(
  rbind,
  fig5_catalog_list
)


all_fig5_pathway_scores <- do.call(
  rbind,
  fig5_pathway_scores_list
)


all_fig5_correlations <- do.call(
  rbind,
  fig5_correlations_list
)


all_fig5_presentation_audit <- if (
  length(
    fig5_presentation_audit_list
  )
) {
  
  do.call(
    rbind,
    fig5_presentation_audit_list
  )
  
} else {
  
  data.frame()
}


fig5_index_dir <- file.path(
  OUT_ROOT,
  "plot_data",
  "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
  "00_INDEX"
)


dir.create(
  fig5_index_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


write_csv(
  all_fig5_catalog,
  file.path(
    fig5_index_dir,
    "FIG5_MASTER_PATHWAY_FEATURE_CATALOG.csv"
  )
)


write_csv(
  all_fig5_pathway_scores,
  file.path(
    fig5_index_dir,
    "FIG5_POOLED_FROZEN_PATHWAY_SCORES.csv"
  )
)


write_csv(
  all_fig5_correlations,
  file.path(
    fig5_index_dir,
    "FIG5_PATHWAY_FEATURE_SPEARMAN_GLOBAL_AND_PHENOTYPES.csv"
  )
)


write_csv(
  all_fig5_presentation_audit,
  file.path(
    fig5_index_dir,
    "FIG5_PRESENTATION_SELECTION_AUDIT.csv"
  )
)


fig5_output_index <- data.frame(
  
  folder = c(
    "00_INDEX",
    "01_ALL",
    "02_PRESENTATION"
  ),
  
  content = c(
    
    paste0(
      "Catálogo pathway-feature, scores congelados y ",
      "correlaciones descriptivas GLOBAL/NP/FD/SO."
    ),
    
    paste0(
      "Todas las features TX/PR/ME/CL asociadas a cada pathway. ",
      "Sin reducción gráfica."
    ),
    
    paste0(
      "Versión adaptativa. Si n_features <= ",
      FIG5_PRESENTATION_MAX_FEATURES,
      ", conserva todas; si excede el máximo, prioriza PR/CL ",
      "y reduce principalmente TX/ME por variabilidad dentro de vista."
    )
  ),
  
  stringsAsFactors =
    FALSE
)


write_csv(
  fig5_output_index,
  file.path(
    fig5_index_dir,
    "FIG5_OUTPUT_INDEX.csv"
  )
)


# ===== 8.15 ACTUALIZAR README =====

cat(
  
  paste0(
    
    "\n",
    
    paste0(
      "PREGUNTA 5 — ¿Cuál es la arquitectura multimodal asociada ",
      "a los pathways biológicos de cada fenotipo?\n"
    ),
    
    "  - Q5 es exclusivamente descriptiva.\n",
    
    paste0(
      "  - GLOBAL combina los 26 individuos sin distinguir ",
      "TRAIN/TEST en la interpretación.\n"
    ),
    
    "  - También se muestran NP, FD y SO por separado.\n",
    
    paste0(
      "  - Los pathways, programas y asociaciones latentes ",
      "ya estaban definidos previamente; Q5 no los redefine.\n"
    ),
    
    paste0(
      "  - FIG5A integra TX + PR + ME + CL alrededor de los ",
      "pathways de cada programa.\n"
    ),
    
    paste0(
      "  - TX se representa mediante membership, condición de driver ",
      "y normalized_anchor_weight.\n"
    ),
    
    paste0(
      "  - PR/ME/CL se representan mediante rho de Spearman ",
      "con el score congelado del pathway.\n"
    ),
    
    paste0(
      "  - FIG5B representa correlaciones Spearman feature-feature ",
      "dentro del conjunto multimodal asociado a cada pathway.\n"
    ),
    paste0(
      "  - FIG5B se exporta en dos representaciones: ",
      "01_OBSERVED y 02_LATENT.\n"
    ),
    
    paste0(
      "  - OBSERVED usa los valores preprocesados originales; ",
      "LATENT usa L %*% Z_oriented_scaled con los factores auditados.\n"
    ),
    
    paste0(
      "  - OBSERVED y LATENT usan exactamente las mismas features; ",
      "LATENT no introduce ninguna selección adicional.\n"
    ),
    
    paste0(
      "  - El clustering se calcula separadamente para OBSERVED y LATENT ",
      "usando GLOBAL y únicamente dentro de cada vista; ",
      "el mismo orden se reutiliza después en NP/FD/SO.\n"
    ),
    paste0(
      "  - Las matrices feature-feature no crean módulos nuevos ",
      "ni seleccionan pathways.\n"
    ),
    
    paste0(
      "  - 01_ALL conserva todas las features disponibles.\n"
    ),
    
    paste0(
      "  - 02_PRESENTATION aplica únicamente reducción gráfica ",
      "cuando se supera BIO_FIG5_PRESENTATION_MAX_FEATURES=",
      FIG5_PRESENTATION_MAX_FEATURES,
      ".\n"
    ),
    
    paste0(
      "  - La reducción de presentación NO usa rho, p-valores ni ",
      "resultados intra-fenotipo; usa únicamente variabilidad descriptiva.\n"
    ),
    
    paste0(
      "  - Las correlaciones intra-fenotipo son informativas y deben ",
      "interpretarse con cautela debido al tamaño muestral, especialmente FD.\n"
    )
  ),
  
  file = file.path(
    OUT_ROOT,
    "READ_ME_FIRST.txt"
  ),
  
  append = TRUE
)

message("")
message("QUESTION 5 COMPLETE")

}



# =============================================================================
# PREGUNTA 6
#
# ¿QUÉ PROGRAMAS Y PATHWAYS SON COMPARTIDOS, ESPECÍFICOS O
# MULTIMODALMENTE REFORZADOS?
#
# PRINCIPIOS:
#   - no redefine pathways ni programas;
#   - conserva P1/P2/P3 como preguntas científicas distintas;
#   - "compartido" se evalúa por análisis y por lado fenotípico;
#   - los pathways se auditan también DENTRO del contexto de su programa;
#   - "multimodalmente reforzado" = pathway TX con contexto crossómico
#     PR/ME/CL ya definido aguas arriba;
#   - las correlaciones GLOBAL/NP/FD/SO son descriptivas y NO se usan
#     para declarar refuerzo;
#   - las figuras de Q3-Q5 se COPIAN, nunca se mueven.
# =============================================================================


# ===== 9.1 HELPERS =====

q6_max_abs <- function(x) {
  
  x <- suppressWarnings(
    as.numeric(x)
  )
  
  x <- x[
    is.finite(x)
  ]
  
  if (!length(x)) {
    return(NA_real_)
  }
  
  max(
    abs(x)
  )
}


q6_occurrence_key <- function(
    family,
    program_id,
    library,
    pathway,
    phenotypic_side
) {
  
  paste(
    family,
    program_id,
    library,
    pathway,
    phenotypic_side,
    sep = "\r"
  )
}


# ===== 9.2 CARGAR MEMBERSHIP PATHWAY -> PROGRAMA =====

load_q6_pathway_membership <- function(
    bio_program_run,
    families
) {
  
  rows <- lapply(
    families,
    
    function(family) {
      
      x <- read_csv_checked(
        file.path(
          bio_program_run,
          family,
          "07_TX_PROGRAM_PATHWAY_MEMBERSHIP.csv"
        ),
        required = c(
          "program_id",
          "membership_type",
          "library",
          "pathway",
          "factor",
          "factor_role",
          "phenotypic_side"
        )
      )
      
      
      x <- x[
        x$membership_type ==
          "PRIMARY",
        ,
        drop = FALSE
      ]
      
      
      x$family <- family
      
      x$program_short <-
        short_program_id(
          x$program_id
        )
      
      x$factor_side <- paste0(
        x$factor,
        "_",
        x$phenotypic_side
      )
      
      x$pathway_label <-
        clean_pathway_label(
          x$pathway
        )
      
      
      unique(
        x[
          ,
          c(
            "family",
            "program_id",
            "program_short",
            "factor",
            "factor_role",
            "phenotypic_side",
            "factor_side",
            "library",
            "pathway",
            "pathway_label"
          ),
          drop = FALSE
        ]
      )
    }
  )
  
  
  do.call(
    rbind,
    rows
  )
}


# ===== 9.3 PREPARAR MASTER DE PROGRAMAS =====

prepare_q6_program_master <- function(
    program_master,
    pathway_membership
) {
  
  pm <- program_master
  
  
  meta <- unique(
    pathway_membership[
      ,
      c(
        "family",
        "program_id",
        "factor",
        "factor_side"
      ),
      drop = FALSE
    ]
  )
  
  
  pm <- merge(
    pm,
    meta,
    by = c(
      "family",
      "program_id"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  
  
  strict_flag <- if (
    "strict_P1_P2_P3_triad" %in%
    names(pm)
  ) {
    
    fig5_as_flag(
      pm$strict_P1_P2_P3_triad
    )
    
  } else {
    
    rep(
      FALSE,
      nrow(pm)
    )
  }
  
  
  n_other <- suppressWarnings(
    as.integer(
      pm$n_other_analyses_aligned
    )
  )
  
  
  pm$program_analysis_scope <- ifelse(
    strict_flag,
    
    "Shared P1-P2-P3 (strict triad)",
    
    ifelse(
      !is.na(n_other) &
        n_other >= 2L,
      
      "Shared P1-P2-P3",
      
      ifelse(
        !is.na(n_other) &
          n_other == 1L,
        
        "Shared two analyses",
        
        "Analysis-specific"
      )
    )
  )
  
  
  pm$program_phenotype_scope <-
    as.character(
      pm$cross_phenotype_status
    )
  
  
  pm
}


# ===== 9.4 PATHWAYS: CONSERVACIÓN + ESPECIFICIDAD + REFUERZO =====

build_q6_pathway_master <- function(
    pathway_membership,
    program_master,
    fig5_catalog,
    fig5_correlations
) {
  
  x <- pathway_membership
  
  
  # ---------------------------------------------------------------------------
  # Identidad exacta del pathway dentro del MISMO lado fenotípico.
  #
  # Librería + nombre de pathway + lado.
  # ---------------------------------------------------------------------------
  
  x$pathway_side_id <- paste(
    x$phenotypic_side,
    x$library,
    x$pathway,
    sep = "\r"
  )
  
  
  family_map <- tapply(
    x$family,
    x$pathway_side_id,
    function(v) {
      paste(
        sort(
          unique(v)
        ),
        collapse = ";"
      )
    }
  )
  
  
  n_family_map <- tapply(
    x$family,
    x$pathway_side_id,
    function(v) {
      length(
        unique(v)
      )
    }
  )
  
  
  x$analyses_with_pathway <-
    as.character(
      family_map[
        x$pathway_side_id
      ]
    )
  
  
  x$n_analyses_with_pathway <-
    as.integer(
      n_family_map[
        x$pathway_side_id
      ]
    )
  
  
  x$pathway_analysis_scope <- ifelse(
    x$n_analyses_with_pathway >= 3L,
    "Shared P1-P2-P3",
    ifelse(
      x$n_analyses_with_pathway == 2L,
      "Shared two analyses",
      "Analysis-specific"
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Especificidad por lado fenotípico.
  #
  # Aquí quitamos temporalmente phenotypic_side de la identidad para preguntar
  # si la MISMA vía aparece también en otros lados.
  # ---------------------------------------------------------------------------
  
  x$pathway_biology_id <- paste(
    x$library,
    x$pathway,
    sep = "\r"
  )
  
  
  side_map <- tapply(
    x$phenotypic_side,
    x$pathway_biology_id,
    function(v) {
      paste(
        sort(
          unique(v)
        ),
        collapse = ";"
      )
    }
  )
  
  
  n_side_map <- tapply(
    x$phenotypic_side,
    x$pathway_biology_id,
    function(v) {
      length(
        unique(v)
      )
    }
  )
  
  
  x$phenotypic_sides_with_pathway <-
    as.character(
      side_map[
        x$pathway_biology_id
      ]
    )
  
  
  x$n_phenotypic_sides_with_pathway <-
    as.integer(
      n_side_map[
        x$pathway_biology_id
      ]
    )
  
  
  x$pathway_phenotype_scope <- ifelse(
    x$n_phenotypic_sides_with_pathway == 1L,
    "Phenotype-side-specific",
    "Cross-phenotype pathway"
  )
  
  
  # ---------------------------------------------------------------------------
  # Contexto del programa padre.
  # ---------------------------------------------------------------------------
  
  pm_keep <- program_master[
    ,
    c(
      "family",
      "program_id",
      "program_analysis_scope",
      "program_phenotype_scope",
      "n_other_analyses_aligned",
      "strict_P1_P2_P3_triad"
    ),
    drop = FALSE
  ]
  
  
  x <- merge(
    x,
    pm_keep,
    by = c(
      "family",
      "program_id"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  
  
  pathway_shared <-
    x$n_analyses_with_pathway >= 2L
  
  
  program_shared <-
    suppressWarnings(
      as.integer(
        x$n_other_analyses_aligned
      )
    ) >= 1L
  
  
  program_shared[
    is.na(program_shared)
  ] <- FALSE
  
  
  x$pathway_program_joint_status <- ifelse(
    pathway_shared &
      program_shared,
    
    "Shared pathway within shared-program context",
    
    ifelse(
      pathway_shared &
        !program_shared,
      
      "Shared pathway across non-aligned programs",
      
      ifelse(
        !pathway_shared &
          program_shared,
        
        "Shared program; pathway analysis-specific",
        
        "Analysis-specific pathway in analysis-specific program"
      )
    )
  )
  
  
  strict_program <- fig5_as_flag(
    x$strict_P1_P2_P3_triad
  )
  
  
  x$strict_threeway_pathway_in_strict_program <-
    x$n_analyses_with_pathway >= 3L &
    strict_program
  
  
  # ---------------------------------------------------------------------------
  # Refuerzo multimodal.
  #
  # Se basa en el catálogo de links crossómicos PR/ME/CL ya definido
  # aguas arriba. NO se utiliza un umbral de rho.
  # ---------------------------------------------------------------------------
  
  catalog <- fig5_catalog
  
  
  catalog$occurrence_key <- q6_occurrence_key(
    catalog$family,
    catalog$program_id,
    catalog$library,
    catalog$pathway,
    catalog$phenotypic_side
  )
  
  
  if (
    "available_in_pooled" %in%
    names(catalog)
  ) {
    
    available <- fig5_as_flag(
      catalog$available_in_pooled
    )
    
  } else {
    
    available <- rep(
      TRUE,
      nrow(catalog)
    )
  }
  
  
  catalog <- catalog[
    available &
      tolower(
        as.character(
          catalog$view_code
        )
      ) %in%
      c(
        "pr",
        "me",
        "cl"
      ),
    ,
    drop = FALSE
  ]
  
  
  x$occurrence_key <- q6_occurrence_key(
    x$family,
    x$program_id,
    x$library,
    x$pathway,
    x$phenotypic_side
  )
  
  
  for (
    v in c(
      "pr",
      "me",
      "cl"
    )
  ) {
    
    x[[paste0(
      "n_",
      v,
      "_features_linked"
    )]] <- 0L
    
    
    x[[paste0(
      "max_abs_cosine_",
      v
    )]] <- NA_real_
  }
  
  
  for (
    i in seq_len(
      nrow(x)
    )
  ) {
    
    z <- catalog[
      catalog$occurrence_key ==
        x$occurrence_key[i],
      ,
      drop = FALSE
    ]
    
    
    for (
      v in c(
        "pr",
        "me",
        "cl"
      )
    ) {
      
      zv <- z[
        tolower(
          as.character(
            z$view_code
          )
        ) == v,
        ,
        drop = FALSE
      ]
      
      
      x[[paste0(
        "n_",
        v,
        "_features_linked"
      )]][i] <-
        length(
          unique(
            zv$feature_original
          )
        )
      
      
      x[[paste0(
        "max_abs_cosine_",
        v
      )]][i] <-
        q6_max_abs(
          zv$pathway_feature_cosine
        )
    }
  }
  
  
  x$n_nonTX_views_linked <-
    as.integer(
      x$n_pr_features_linked > 0L
    ) +
    as.integer(
      x$n_me_features_linked > 0L
    ) +
    as.integer(
      x$n_cl_features_linked > 0L
    )
  
  
  x$multimodal_reinforcement_status <- ifelse(
    x$n_nonTX_views_linked > 0L,
    "Multimodally reinforced",
    "TX-only in Q5 crossomic context"
  )
  
  
  # ---------------------------------------------------------------------------
  # Resumen descriptivo de correlaciones por grupo.
  #
  # NO cambia la clasificación anterior.
  # ---------------------------------------------------------------------------
  
  corr <- fig5_correlations
  
  
  corr$occurrence_key <- q6_occurrence_key(
    corr$family,
    corr$program_id,
    corr$library,
    corr$pathway,
    corr$phenotypic_side
  )
  
  
  for (
    scope in c(
      "GLOBAL",
      "NP",
      "FD",
      "SO"
    )
  ) {
    
    x[[paste0(
      "n_nonTX_views_with_finite_rho_",
      scope
    )]] <- 0L
    
    
    x[[paste0(
      "max_abs_rho_",
      scope
    )]] <- NA_real_
  }
  
  
  for (
    i in seq_len(
      nrow(x)
    )
  ) {
    
    for (
      scope in c(
        "GLOBAL",
        "NP",
        "FD",
        "SO"
      )
    ) {
      
      z <- corr[
        corr$occurrence_key ==
          x$occurrence_key[i] &
          corr$scope ==
          scope &
          is.finite(
            suppressWarnings(
              as.numeric(
                corr$rho_spearman
              )
            )
          ),
        ,
        drop = FALSE
      ]
      
      
      x[[paste0(
        "n_nonTX_views_with_finite_rho_",
        scope
      )]][i] <-
        length(
          unique(
            tolower(
              as.character(
                z$view_code
              )
            )
          )
        )
      
      
      x[[paste0(
        "max_abs_rho_",
        scope
      )]][i] <-
        q6_max_abs(
          z$rho_spearman
        )
    }
  }
  
  
  x <- x[
    order(
      x$phenotypic_side,
      -x$n_analyses_with_pathway,
      -x$n_nonTX_views_linked,
      x$family,
      x$program_short,
      x$library,
      x$pathway
    ),
    ,
    drop = FALSE
  ]
  
  
  x
}


# ===== 9.5 RESUMEN Q6 A NIVEL DE PROGRAMA =====

build_q6_program_master <- function(
    program_master,
    pathway_master
) {
  
  pm <- program_master
  
  
  pm$n_primary_pathways_q6 <- 0L
  
  pm$n_pathways_shared_P1_P2_P3 <- 0L
  
  pm$n_pathways_shared_two_analyses <- 0L
  
  pm$n_pathways_analysis_specific <- 0L
  
  pm$n_pathways_multimodally_reinforced <- 0L
  
  pm$n_nonTX_views_linked <- 0L
  
  pm$has_PR_support <- FALSE
  
  pm$has_ME_support <- FALSE
  
  pm$has_CL_support <- FALSE
  
  
  for (
    scope in c(
      "GLOBAL",
      "NP",
      "FD",
      "SO"
    )
  ) {
    
    pm[[paste0(
      "max_abs_rho_",
      scope
    )]] <- NA_real_
  }
  
  
  for (
    i in seq_len(
      nrow(pm)
    )
  ) {
    
    z <- pathway_master[
      pathway_master$family ==
        pm$family[i] &
        pathway_master$program_id ==
        pm$program_id[i],
      ,
      drop = FALSE
    ]
    
    
    if (!nrow(z)) {
      next
    }
    
    
    path_ids <- paste(
      z$library,
      z$pathway,
      sep = "\r"
    )
    
    
    pm$n_primary_pathways_q6[i] <-
      length(
        unique(
          path_ids
        )
      )
    
    
    pm$n_pathways_shared_P1_P2_P3[i] <-
      sum(
        z$n_analyses_with_pathway >= 3L
      )
    
    
    pm$n_pathways_shared_two_analyses[i] <-
      sum(
        z$n_analyses_with_pathway == 2L
      )
    
    
    pm$n_pathways_analysis_specific[i] <-
      sum(
        z$n_analyses_with_pathway == 1L
      )
    
    
    pm$n_pathways_multimodally_reinforced[i] <-
      sum(
        z$n_nonTX_views_linked > 0L
      )
    
    
    pm$has_PR_support[i] <-
      any(
        z$n_pr_features_linked > 0L
      )
    
    
    pm$has_ME_support[i] <-
      any(
        z$n_me_features_linked > 0L
      )
    
    
    pm$has_CL_support[i] <-
      any(
        z$n_cl_features_linked > 0L
      )
    
    
    pm$n_nonTX_views_linked[i] <-
      as.integer(
        pm$has_PR_support[i]
      ) +
      as.integer(
        pm$has_ME_support[i]
      ) +
      as.integer(
        pm$has_CL_support[i]
      )
    
    
    for (
      scope in c(
        "GLOBAL",
        "NP",
        "FD",
        "SO"
      )
    ) {
      
      pm[[paste0(
        "max_abs_rho_",
        scope
      )]][i] <-
        q6_max_abs(
          z[[paste0(
            "max_abs_rho_",
            scope
          )]]
        )
    }
  }
  
  
  pm$program_multimodal_reinforcement_status <-
    ifelse(
      pm$n_nonTX_views_linked > 0L,
      "Multimodally reinforced",
      "TX-only in Q5 crossomic context"
    )
  
  
  pm$program_pathway_conservation <- ifelse(
    pm$n_pathways_shared_P1_P2_P3 > 0L,
    
    "Contains pathways shared P1-P2-P3",
    
    ifelse(
      pm$n_pathways_shared_two_analyses > 0L,
      "Contains pathways shared by two analyses",
      "All constituent pathways are analysis-specific"
    )
  )
  
  
  pm
}


# ===== 9.6 MATRIZ RESUMEN PATHWAY: P1/P2/P3 | PR/ME/CL =====

build_q6_pathway_presence_summary <- function(
    pathway_master
) {
  
  ids <- unique(
    pathway_master$pathway_side_id
  )
  
  
  rows <- lapply(
    ids,
    
    function(id) {
      
      z <- pathway_master[
        pathway_master$pathway_side_id ==
          id,
        ,
        drop = FALSE
      ]
      
      
      data.frame(
        
        pathway_side_id =
          id,
        
        phenotypic_side =
          z$phenotypic_side[1L],
        
        library =
          z$library[1L],
        
        pathway =
          z$pathway[1L],
        
        pathway_label =
          z$pathway_label[1L],
        
        P1 =
          "P1" %in%
          z$family,
        
        P2 =
          "P2" %in%
          z$family,
        
        P3 =
          "P3" %in%
          z$family,
        
        PR =
          any(
            z$n_pr_features_linked >
              0L
          ),
        
        ME =
          any(
            z$n_me_features_linked >
              0L
          ),
        
        CL =
          any(
            z$n_cl_features_linked >
              0L
          ),
        
        n_analyses_with_pathway =
          max(
            z$n_analyses_with_pathway
          ),
        
        n_nonTX_views_linked =
          max(
            z$n_nonTX_views_linked
          ),
        
        pathway_analysis_scope =
          z$pathway_analysis_scope[1L],
        
        pathway_phenotype_scope =
          z$pathway_phenotype_scope[1L],
        
        stringsAsFactors =
          FALSE
      )
    }
  )
  
  
  do.call(
    rbind,
    rows
  )
}


make_q6_pathway_presence_figures <- function(
    summary_dat,
    outdir
) {
  
  sides <- unique(
    summary_dat$phenotypic_side
  )
  
  
  for (
    side_i in sides
  ) {
    
    z <- summary_dat[
      summary_dat$phenotypic_side ==
        side_i,
      ,
      drop = FALSE
    ]
    
    
    z <- z[
      order(
        -z$n_analyses_with_pathway,
        -z$n_nonTX_views_linked,
        z$pathway_label
      ),
      ,
      drop = FALSE
    ]
    
    
    z$row_label <- make.unique(
      paste0(
        z$pathway_label,
        " [",
        z$library,
        "]"
      )
    )
    
    
    long <- do.call(
      rbind,
      
      lapply(
        c(
          "P1",
          "P2",
          "P3",
          "PR",
          "ME",
          "CL"
        ),
        
        function(d) {
          
          data.frame(
            row_label =
              z$row_label,
            dimension =
              d,
            present =
              as.logical(
                z[[d]]
              ),
            stringsAsFactors =
              FALSE
          )
        }
      )
    )
    
    
    long$row_label <- factor(
      long$row_label,
      levels =
        rev(
          z$row_label
        )
    )
    
    
    long$dimension <- factor(
      long$dimension,
      levels = c(
        "P1",
        "P2",
        "P3",
        "PR",
        "ME",
        "CL"
      )
    )
    
    
    p <- ggplot2::ggplot(
      long,
      ggplot2::aes(
        x = dimension,
        y = row_label,
        fill = present
      )
    ) +
      
      ggplot2::geom_tile(
        color = "grey80",
        linewidth = 0.35
      ) +
      
      ggplot2::scale_fill_manual(
        values = c(
          "FALSE" = "white",
          "TRUE" = "grey25"
        ),
        guide = "none"
      ) +
      
      ggplot2::labs(
        title = paste0(
          "Pregunta 6 — ",
          phenotype_label(
            side_i
          )
        ),
        subtitle =
          "P1/P2/P3 = presencia TX; PR/ME/CL = contexto crossómico asociado.",
        x = NULL,
        y = NULL
      ) +
      
      fig_theme(11) +
      
      ggplot2::theme(
        panel.grid =
          ggplot2::element_blank(),
        axis.text.y =
          ggplot2::element_text(
            size = 7.5
          ),
        axis.text.x =
          ggplot2::element_text(
            size = 10,
            face = "bold"
          )
      )
    
    
    dir.create(
      file.path(
        outdir,
        "00_INDEX"
      ),
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    
    ggplot2::ggsave(
      filename = file.path(
        outdir,
        "00_INDEX",
        paste0(
          "FIG6A_",
          safe_filename(
            side_i
          ),
          "_PATHWAY_ANALYSIS_MULTIMODAL_MATRIX.png"
        )
      ),
      plot = p,
      width = 9,
      height = max(
        6,
        min(
          30,
          3 +
            0.28 *
            nrow(z)
        )
      ),
      dpi = FIG_DPI,
      limitsize = FALSE
    )
  }
}


# ===== 9.7 COPIAR FIGURAS YA EXISTENTES SIN MOVERLAS =====

q6_copy_manifest_rows <- list()

q6_copy_rr <- 1L


q6_copy_files <- function(
    sources,
    destination_dir,
    category,
    entity_type,
    entity_id
) {
  
  sources <- unique(
    as.character(
      sources
    )
  )
  
  
  sources <- sources[
    file.exists(
      sources
    )
  ]
  
  
  if (!length(sources)) {
    return(
      invisible(NULL)
    )
  }
  
  
  dir.create(
    destination_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  for (
    src in sources
  ) {
    
    dst <- file.path(
      destination_dir,
      basename(src)
    )
    
    
    ok <- file.copy(
      src,
      dst,
      overwrite = TRUE
    )
    
    
    q6_copy_manifest_rows[[
      q6_copy_rr
    ]] <<- data.frame(
      
      category =
        category,
      
      entity_type =
        entity_type,
      
      entity_id =
        entity_id,
      
      source =
        src,
      
      destination =
        dst,
      
      copied =
        isTRUE(ok),
      
      stringsAsFactors =
        FALSE
    )
    
    
    q6_copy_rr <<-
      q6_copy_rr + 1L
  }
  
  
  invisible(NULL)
}


q6_program_categories <- function(
    row
) {
  
  n_other <- suppressWarnings(
    as.integer(
      row$n_other_analyses_aligned
    )
  )
  
  
  categories <- if (
    !is.na(n_other) &&
    n_other >= 2L
  ) {
    
    "01_COMMON/01_P1_P2_P3"
    
  } else if (
    !is.na(n_other) &&
    n_other == 1L
  ) {
    
    "01_COMMON/02_TWO_ANALYSES"
    
  } else {
    
    "02_ANALYSIS_SPECIFIC"
  }
  
  
  if (
    row$program_multimodal_reinforcement_status ==
    "Multimodally reinforced"
  ) {
    
    categories <- c(
      categories,
      "03_MULTIMODALLY_REINFORCED"
    )
  }
  
  
  if (
    row$program_phenotype_scope ==
    "Phenotype-side-specific"
  ) {
    
    categories <- c(
      categories,
      "04_PHENOTYPE_SPECIFIC"
    )
  }
  
  
  unique(
    categories
  )
}


q6_pathway_categories <- function(
    row
) {
  
  categories <- if (
    row$n_analyses_with_pathway >= 3L
  ) {
    
    "01_COMMON/01_P1_P2_P3"
    
  } else if (
    row$n_analyses_with_pathway == 2L
  ) {
    
    "01_COMMON/02_TWO_ANALYSES"
    
  } else {
    
    "02_ANALYSIS_SPECIFIC"
  }
  
  
  if (
    row$multimodal_reinforcement_status ==
    "Multimodally reinforced"
  ) {
    
    categories <- c(
      categories,
      "03_MULTIMODALLY_REINFORCED"
    )
  }
  
  
  if (
    row$pathway_phenotype_scope ==
    "Phenotype-side-specific"
  ) {
    
    categories <- c(
      categories,
      "04_PHENOTYPE_SPECIFIC"
    )
  }
  
  
  unique(
    categories
  )
}


# ===== 9.8 EJECUTAR PREGUNTA 6 =====

if (run_question(6L)) {
  
  message("")
  message("============================================================")
  message("METATEST — QUESTION 6")
  message("Common / specific / multimodally reinforced biology")
  message("============================================================")
  
  
  q6_root_fig <- file.path(
    OUT_ROOT,
    "figures",
    "FIG6_COMMON_SPECIFIC_REINFORCED"
  )
  
  
  q6_root_dat <- file.path(
    OUT_ROOT,
    "plot_data",
    "FIG6_COMMON_SPECIFIC_REINFORCED"
  )
  
  
  dir.create(
    q6_root_fig,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  dir.create(
    q6_root_dat,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # ---------------------------------------------------------------------------
  # Leer resultados congelados de Q3 y Q5.
  # Esto permite arrancar directamente con BIO_FIG_START_QUESTION=6.
  # ---------------------------------------------------------------------------
  
  q6_program_master_path <- file.path(
    OUT_ROOT,
    "plot_data",
    "FIG3_TX_PROGRAMS",
    "00_INDEX",
    "FIG3_MASTER_PROGRAM_AUDIT.csv"
  )
  
  
  q6_catalog_path <- file.path(
    OUT_ROOT,
    "plot_data",
    "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
    "00_INDEX",
    "FIG5_MASTER_PATHWAY_FEATURE_CATALOG.csv"
  )
  
  
  q6_corr_path <- file.path(
    OUT_ROOT,
    "plot_data",
    "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
    "00_INDEX",
    "FIG5_PATHWAY_FEATURE_SPEARMAN_GLOBAL_AND_PHENOTYPES.csv"
  )
  
  
  program_master_q3 <- read_csv_checked(
    q6_program_master_path
  )
  
  
  fig5_catalog_q6 <- read_csv_checked(
    q6_catalog_path
  )
  
  
  fig5_corr_q6 <- read_csv_checked(
    q6_corr_path
  )
  
  
  pathway_membership_q6 <-
    load_q6_pathway_membership(
      bio_program_run =
        BIO_RUN,
      families =
        FAMILIES
    )
  
  
  program_master_q6 <-
    prepare_q6_program_master(
      program_master =
        program_master_q3,
      pathway_membership =
        pathway_membership_q6
    )
  
  
  pathway_master_q6 <-
    build_q6_pathway_master(
      pathway_membership =
        pathway_membership_q6,
      program_master =
        program_master_q6,
      fig5_catalog =
        fig5_catalog_q6,
      fig5_correlations =
        fig5_corr_q6
    )
  
  
  program_master_q6 <-
    build_q6_program_master(
      program_master =
        program_master_q6,
      pathway_master =
        pathway_master_q6
    )
  
  
  pathway_presence_q6 <-
    build_q6_pathway_presence_summary(
      pathway_master_q6
    )
  
  
  # ---------------------------------------------------------------------------
  # Tablas maestras.
  # ---------------------------------------------------------------------------
  
  dir.create(
    file.path(
      q6_root_dat,
      "00_INDEX"
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  write_csv(
    program_master_q6,
    file.path(
      q6_root_dat,
      "00_INDEX",
      "Q6_PROGRAM_MASTER.csv"
    )
  )
  
  
  write_csv(
    pathway_master_q6,
    file.path(
      q6_root_dat,
      "00_INDEX",
      "Q6_PATHWAY_MASTER.csv"
    )
  )
  
  
  write_csv(
    pathway_presence_q6,
    file.path(
      q6_root_dat,
      "00_INDEX",
      "Q6_PATHWAY_ANALYSIS_MULTIMODAL_MATRIX.csv"
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Nueva figura resumen.
  # ---------------------------------------------------------------------------
  
  make_q6_pathway_presence_figures(
    summary_dat =
      pathway_presence_q6,
    outdir =
      q6_root_fig
  )
  
  
  # ---------------------------------------------------------------------------
  # Reorganizar COPIAS de figuras de PROGRAMAS.
  # ---------------------------------------------------------------------------
  
  for (
    i in seq_len(
      nrow(
        program_master_q6
      )
    )
  ) {
    
    r <- program_master_q6[
      i,
      ,
      drop = FALSE
    ]
    
    
    family_i <-
      as.character(
        r$family
      )
    
    side_i <-
      as.character(
        r$phenotypic_side
      )
    
    factor_side_i <-
      as.character(
        r$factor_side
      )
    
    program_short_i <-
      as.character(
        r$program_short
      )
    
    
    src_q3 <- file.path(
      OUT_ROOT,
      "figures",
      "FIG3_TX_PROGRAMS",
      "01_WITHIN_ANALYSIS",
      family_i,
      safe_filename(
        factor_side_i
      ),
      paste0(
        "FIG3A_",
        family_i,
        "_",
        safe_filename(
          factor_side_i
        ),
        "_PROGRAM_MEMBERSHIP.png"
      )
    )
    
    
    src_q4 <- file.path(
      OUT_ROOT,
      "figures",
      "FIG4_TX_GENES",
      "01B_WITHIN_PROGRAM_PRESENTATION",
      family_i,
      safe_filename(
        factor_side_i
      ),
      safe_filename(
        program_short_i
      ),
      "FIG4A_PRESENTATION_DRIVER_WEIGHT.png"
    )
    
    
    src_q5_dir <- file.path(
      OUT_ROOT,
      "figures",
      "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
      "02_PRESENTATION",
      family_i,
      safe_filename(
        factor_side_i
      ),
      safe_filename(
        program_short_i
      ),
      "01_PROGRAM_PATHWAY_MAP"
    )
    
    
    src_q5 <- if (
      dir.exists(
        src_q5_dir
      )
    ) {
      
      list.files(
        src_q5_dir,
        pattern = "\\.png$",
        recursive = TRUE,
        full.names = TRUE
      )
      
    } else {
      
      character()
    }
    
    
    sources <- c(
      src_q3,
      src_q4,
      src_q5
    )
    
    
    categories <-
      q6_program_categories(
        r
      )
    
    
    for (
      category_i in categories
    ) {
      
      destination <- file.path(
        q6_root_fig,
        category_i,
        "PROGRAMS",
        safe_filename(
          side_i
        ),
        family_i,
        safe_filename(
          program_short_i
        )
      )
      
      
      q6_copy_files(
        sources =
          sources,
        destination_dir =
          destination,
        category =
          category_i,
        entity_type =
          "PROGRAM",
        entity_id =
          as.character(
            r$program_id
          )
      )
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # Reorganizar COPIAS de matrices de PATHWAYS.
  # ---------------------------------------------------------------------------
  
  for (
    i in seq_len(
      nrow(
        pathway_master_q6
      )
    )
  ) {
    
    r <- pathway_master_q6[
      i,
      ,
      drop = FALSE
    ]
    
    
    pathway_folder <- paste0(
      safe_filename(
        r$library
      ),
      "__",
      safe_filename(
        clean_pathway_label(
          r$pathway
        )
      )
    )
    
    
    source_dir <- file.path(
      OUT_ROOT,
      "figures",
      "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
      "02_PRESENTATION",
      as.character(
        r$family
      ),
      safe_filename(
        r$factor_side
      ),
      safe_filename(
        r$program_short
      ),
      "02_PATHWAY_CORRELATIONS",
      pathway_folder
    )
    
    
    sources <- if (
      dir.exists(
        source_dir
      )
    ) {
      
      list.files(
        source_dir,
        pattern = "\\.png$",
        recursive = TRUE,
        full.names = TRUE
      )
      
    } else {
      
      character()
    }
    
    
    categories <-
      q6_pathway_categories(
        r
      )
    
    
    for (
      category_i in categories
    ) {
      
      destination <- file.path(
        q6_root_fig,
        category_i,
        "PATHWAYS",
        safe_filename(
          r$phenotypic_side
        ),
        as.character(
          r$family
        ),
        safe_filename(
          r$program_short
        ),
        pathway_folder
      )
      
      
      q6_copy_files(
        sources =
          sources,
        destination_dir =
          destination,
        category =
          category_i,
        entity_type =
          "PATHWAY",
        entity_id =
          paste(
            r$family,
            r$program_id,
            r$library,
            r$pathway,
            sep = " | "
          )
      )
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # Manifest de copias.
  # ---------------------------------------------------------------------------
  
  q6_copy_manifest <- if (
    length(
      q6_copy_manifest_rows
    )
  ) {
    
    do.call(
      rbind,
      q6_copy_manifest_rows
    )
    
  } else {
    
    data.frame(
      category = character(),
      entity_type = character(),
      entity_id = character(),
      source = character(),
      destination = character(),
      copied = logical(),
      stringsAsFactors = FALSE
    )
  }
  
  
  write_csv(
    q6_copy_manifest,
    file.path(
      q6_root_dat,
      "00_INDEX",
      "Q6_COPIED_FIGURES_MANIFEST.csv"
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # README.
  # ---------------------------------------------------------------------------
  
  cat(
    
    paste0(
      "\n",
      "PREGUNTA 6 — programas y pathways compartidos, específicos ",
      "o multimodalmente reforzados.\n",
      "  - La clasificación se realiza a nivel de programa Y de pathway.\n",
      "  - Un pathway puede ser compartido aunque su programa padre no lo sea.\n",
      "  - Un programa compartido puede contener pathways conservados y otros ",
      "específicos de un análisis.\n",
      "  - Multimodally reinforced significa presencia de contexto crossómico ",
      "PR/ME/CL previamente definido; no implica un pathway emergente nuevo.\n",
      "  - GLOBAL/NP/FD/SO rho se conserva como descripción y no determina ",
      "la clasificación de refuerzo.\n",
      "  - Las figuras de Q3-Q5 se COPIAN a FIG6; los originales no se mueven.\n"
    ),
    
    file =
      file.path(
        OUT_ROOT,
        "READ_ME_FIRST.txt"
      ),
    
    append =
      TRUE
  )
  
  
  message("")
  message("QUESTION 6 COMPLETE")
  message(
    "Output: ",
    q6_root_fig
  )
}



# =============================================================================
# PREGUNTA 7
#
# ¿CÓMO SE MANIFIESTAN LOS PATHWAYS / PROGRAMAS Y SU ARQUITECTURA
# MULTIMODAL INDIVIDUO POR INDIVIDUO?
#
# PRINCIPIOS:
#   - no redefine pathways, programas ni features;
#   - reutiliza los pathway scores congelados de Q5;
#   - TRAIN define cualquier centrado/escalado;
#   - TEST recibe exactamente la transformación definida en TRAIN;
#   - PR/ME/CL usan exclusivamente los links crossómicos ya definidos;
#   - las correlaciones de Q5 NO seleccionan features en Q7;
#   - no clasifica individuos como típicos/atípicos;
#   - Q7 es descriptiva.
# =============================================================================


# ===== 10.1 HELPERS =====

q7_normalize_split <- function(x) {
  
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  ifelse(
    grepl(
      "TRAIN",
      x
    ),
    "TRAIN",
    ifelse(
      grepl(
        "TEST",
        x
      ),
      "TEST",
      NA_character_
    )
  )
}


q7_order_samples <- function(meta) {
  
  meta <- unique(
    meta[
      ,
      c(
        "sample_id",
        "grupo",
        "split"
      ),
      drop = FALSE
    ]
  )
  
  
  meta <- meta[
    order(
      match(
        as.character(
          meta$grupo
        ),
        c(
          "NP",
          "FD",
          "SO"
        )
      ),
      match(
        meta$split,
        c(
          "TRAIN",
          "TEST"
        )
      ),
      meta$sample_id
    ),
    ,
    drop = FALSE
  ]
  
  
  meta$sample_label <- paste0(
    meta$sample_id,
    "\n",
    as.character(
      meta$grupo
    ),
    "-",
    ifelse(
      meta$split ==
        "TRAIN",
      "TR",
      "TE"
    )
  )
  
  
  meta
}


q7_train_standardize_long <- function(
    dat,
    value_col,
    key_cols,
    output_col
) {
  
  assert(
    all(
      c(
        value_col,
        key_cols,
        "split"
      ) %in%
        names(dat)
    ),
    "Q7 standardization input is missing required columns"
  )
  
  
  key <- do.call(
    paste,
    c(
      dat[
        ,
        key_cols,
        drop = FALSE
      ],
      sep = "\r"
    )
  )
  
  
  out <- rep(
    NA_real_,
    nrow(dat)
  )
  
  train_mean <- rep(
    NA_real_,
    nrow(dat)
  )
  
  train_sd <- rep(
    NA_real_,
    nrow(dat)
  )
  
  
  for (
    k in unique(key)
  ) {
    
    idx <- which(
      key == k
    )
    
    
    values <- suppressWarnings(
      as.numeric(
        dat[[value_col]][idx]
      )
    )
    
    
    is_train <-
      dat$split[idx] ==
      "TRAIN" &
      is.finite(values)
    
    
    if (
      sum(is_train) <
      2L
    ) {
      next
    }
    
    
    mu <- mean(
      values[is_train],
      na.rm = TRUE
    )
    
    
    s <- stats::sd(
      values[is_train],
      na.rm = TRUE
    )
    
    
    if (
      !is.finite(s) ||
      s <= 0
    ) {
      next
    }
    
    
    out[idx] <-
      (
        values -
          mu
      ) /
      s
    
    
    train_mean[idx] <- mu
    
    train_sd[idx] <- s
  }
  
  
  dat[[output_col]] <-
    out
  
  
  dat[[paste0(
    output_col,
    "_TRAIN_mean"
  )]] <- train_mean
  
  
  dat[[paste0(
    output_col,
    "_TRAIN_sd"
  )]] <- train_sd
  
  
  dat
}


q7_clip_z <- function(x) {
  
  pmax(
    -FIG7_Z_LIMIT,
    pmin(
      FIG7_Z_LIMIT,
      as.numeric(x)
    )
  )
}


# ===== 10.2 METADATA TRAIN / TEST =====

load_q7_sample_meta <- function(
    inp
) {
  
  score_path <- file.path(
    inp$latent_run,
    "tables",
    "02_LATENT_SCORES_TRAIN_AND_PROJECTED_TEST.csv"
  )
  
  
  scores <- read_csv_checked(
    score_path,
    required = c(
      "sample",
      "split",
      "group"
    )
  )
  
  
  scores$split <-
    q7_normalize_split(
      scores$split
    )
  
  
  meta <- unique(
    scores[
      ,
      c(
        "sample",
        "split"
      ),
      drop = FALSE
    ]
  )
  
  
  names(meta)[
    names(meta) ==
      "sample"
  ] <- "sample_id"
  
  
  sample_ids <-
    as.character(
      inp$sample_meta$sample_id
    )
  
  
  idx <- match(
    sample_ids,
    meta$sample_id
  )
  
  
  assert(
    all(
      !is.na(idx)
    ),
    paste(
      "Q7 could not map TRAIN/TEST split for",
      inp$family
    )
  )
  
  
  out <- data.frame(
    
    family =
      inp$family,
    
    sample_id =
      sample_ids,
    
    grupo =
      as.character(
        inp$sample_meta$grupo
      ),
    
    split =
      meta$split[idx],
    
    stringsAsFactors =
      FALSE
  )
  
  
  assert(
    all(
      out$split %in%
        c(
          "TRAIN",
          "TEST"
        )
    ),
    paste(
      "Q7 has undefined TRAIN/TEST samples in",
      inp$family
    )
  )
  
  
  out
}


# ===== 10.3 PATHWAY SCORES INDIVIDUALES =====

prepare_q7_pathway_scores <- function(
    pathway_scores,
    sample_meta
) {
  
  meta <- sample_meta[
    ,
    c(
      "family",
      "sample_id",
      "split"
    ),
    drop = FALSE
  ]
  
  
  x <- merge(
    pathway_scores,
    meta,
    by = c(
      "family",
      "sample_id"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  
  
  assert(
    !any(
      is.na(
        x$split
      )
    ),
    paste0(
      "Q7 pathway scores contain samples without ",
      "TRAIN/TEST assignment"
    )
  )
  
  
  x <-
    q7_train_standardize_long(
      
      dat =
        x,
      
      value_col =
        "pathway_score",
      
      key_cols =
        c(
          "family",
          "pathway_key"
        ),
      
      output_col =
        "pathway_z_train"
    )
  
  
  x$pathway_z_plot <-
    q7_clip_z(
      x$pathway_z_train
    )
  
  
  x
}


# ===== 10.4 FIGURA 7A — FINGERPRINT INDIVIDUAL =====

make_q7_pathway_fingerprint <- function(
    pathway_scores,
    family,
    factor_side,
    outdir
) {
  
  z <- pathway_scores[
    pathway_scores$family ==
      family &
      pathway_scores$factor_side ==
      factor_side,
    ,
    drop = FALSE
  ]
  
  
  if (!nrow(z)) {
    return(
      invisible(NULL)
    )
  }
  
  
  sample_meta <-
    q7_order_samples(
      z
    )
  
  
  z$sample_plot <- factor(
    z$sample_id,
    levels =
      sample_meta$sample_id,
    labels =
      sample_meta$sample_label
  )
  
  
  pathway_meta <- unique(
    z[
      ,
      c(
        "program_short",
        "library",
        "pathway",
        "pathway_label"
      ),
      drop = FALSE
    ]
  )
  
  
  pathway_meta <- pathway_meta[
    order(
      pathway_meta$program_short,
      pathway_meta$library,
      pathway_meta$pathway_label
    ),
    ,
    drop = FALSE
  ]
  
  
  pathway_meta$row_label <- make.unique(
    paste0(
      pathway_meta$program_short,
      " | ",
      pathway_meta$pathway_label,
      " [",
      pathway_meta$library,
      "]"
    )
  )
  
  
  key_meta <- paste(
    pathway_meta$program_short,
    pathway_meta$library,
    pathway_meta$pathway,
    sep = "\r"
  )
  
  
  key_z <- paste(
    z$program_short,
    z$library,
    z$pathway,
    sep = "\r"
  )
  
  
  z$row_label <-
    pathway_meta$row_label[
      match(
        key_z,
        key_meta
      )
    ]
  
  
  z$pathway_plot <- factor(
    z$row_label,
    levels =
      rev(
        pathway_meta$row_label
      )
  )
  
  
  p <- ggplot2::ggplot(
    z,
    ggplot2::aes(
      x = sample_plot,
      y = pathway_plot,
      fill = pathway_z_plot
    )
  ) +
    
    ggplot2::geom_tile(
      color = "grey88",
      linewidth = 0.25
    ) +
    
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(
        -FIG7_Z_LIMIT,
        FIG7_Z_LIMIT
      ),
      name = "Pathway score\nZ(TRAIN)"
    ) +
    
    ggplot2::labs(
      
      title = paste0(
        family,
        " — ",
        phenotype_label(
          z$phenotypic_side[1L]
        ),
        " — manifestación individual"
      ),
      
      subtitle = paste0(
        "Cada columna = individuo; cada fila = pathway dentro de su programa. ",
        "TR = TRAIN; TE = TEST."
      ),
      
      x = NULL,
      y = NULL,
      
      caption = paste0(
        "Pathway score congelado. La estandarización visual utiliza ",
        "exclusivamente media y SD de TRAIN y se aplica sin cambios a TEST. ",
        "Color limitado gráficamente a ±",
        FIG7_Z_LIMIT,
        " SD."
      )
    ) +
    
    fig_theme(10) +
    
    ggplot2::theme(
      
      panel.grid =
        ggplot2::element_blank(),
      
      axis.text.x =
        ggplot2::element_text(
          angle = 60,
          hjust = 1,
          vjust = 1,
          size = 7
        ),
      
      axis.text.y =
        ggplot2::element_text(
          size = 7
        ),
      
      legend.position =
        "bottom"
    )
  
  
  dir_fig <- file.path(
    outdir,
    "figures",
    "FIG7_INDIVIDUAL_MANIFESTATION",
    "01_PATHWAY_FINGERPRINT",
    family,
    safe_filename(
      factor_side
    )
  )
  
  
  dir_dat <- file.path(
    outdir,
    "plot_data",
    "FIG7_INDIVIDUAL_MANIFESTATION",
    "01_PATHWAY_FINGERPRINT",
    family,
    safe_filename(
      factor_side
    )
  )
  
  
  dir.create(
    dir_fig,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  dir.create(
    dir_dat,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  ggplot2::ggsave(
    
    filename = file.path(
      dir_fig,
      "FIG7A_INDIVIDUAL_PATHWAY_FINGERPRINT.png"
    ),
    
    plot = p,
    
    width = 16,
    
    height = max(
      6,
      min(
        28,
        3 +
          0.30 *
          nrow(
            pathway_meta
          )
      )
    ),
    
    dpi = FIG_DPI,
    
    limitsize = FALSE
  )
  
  
  write_csv(
    z,
    file.path(
      dir_dat,
      "FIG7A_INDIVIDUAL_PATHWAY_FINGERPRINT.csv"
    )
  )
  
  
  invisible(p)
}


# ===== 10.5 SCORE MULTIMODAL INDIVIDUAL =====
#
# TX:
#   score congelado del pathway.
#
# PR / ME / CL:
#
#             sum_j w_j z_ij
#   Score = -------------------
#               sum_j |w_j|
#
# donde:
#   z_ij = valor de feature estandarizado usando SOLO TRAIN
#   w_j  = pathway_feature_cosine congelado aguas arriba.
#
# No se utilizan rho GLOBAL/NP/FD/SO para seleccionar features.
# =============================================================================

build_q7_multimodal_manifestation <- function(
    q6_pathways,
    fig5_catalog,
    pathway_scores,
    inputs_list,
    sample_meta
) {
  
  required_q6 <- c(
    "family",
    "program_id",
    "program_short",
    "factor_side",
    "phenotypic_side",
    "library",
    "pathway",
    "pathway_label",
    "n_nonTX_views_linked"
  )
  
  
  assert(
    all(
      required_q6 %in%
        names(
          q6_pathways
        )
    ),
    paste(
      "Q7 requires Q6_PATHWAY_MASTER columns:",
      paste(
        setdiff(
          required_q6,
          names(
            q6_pathways
          )
        ),
        collapse = ", "
      )
    )
  )
  
  
  q <- q6_pathways[
    suppressWarnings(
      as.integer(
        q6_pathways$n_nonTX_views_linked
      )
    ) > 0L,
    ,
    drop = FALSE
  ]
  
  
  if (!nrow(q)) {
    return(
      data.frame()
    )
  }
  
  
  q <- unique(
    q[
      ,
      required_q6,
      drop = FALSE
    ]
  )
  
  
  rows <- list()
  
  rr <- 1L
  
  
  for (
    i in seq_len(
      nrow(q)
    )
  ) {
    
    qi <- q[
      i,
      ,
      drop = FALSE
    ]
    
    
    family_i <-
      as.character(
        qi$family
      )
    
    
    inp <-
      inputs_list[[
        family_i
      ]]
    
    
    meta <- sample_meta[
      sample_meta$family ==
        family_i,
      ,
      drop = FALSE
    ]
    
    
    meta <-
      q7_order_samples(
        meta
      )
    
    
    zcat <- fig5_catalog[
      fig5_catalog$family ==
        family_i &
        fig5_catalog$program_id ==
        qi$program_id &
        fig5_catalog$library ==
        qi$library &
        fig5_catalog$pathway ==
        qi$pathway &
        fig5_catalog$phenotypic_side ==
        qi$phenotypic_side,
      ,
      drop = FALSE
    ]
    
    
    if (!nrow(zcat)) {
      
      warning(
        "Q7 could not map Q6 pathway to Q5 catalog: ",
        family_i,
        " | ",
        qi$library,
        " | ",
        qi$pathway
      )
      
      next
    }
    
    
    pathway_keys <- unique(
      as.character(
        zcat$pathway_key
      )
    )
    
    
    assert(
      length(
        pathway_keys
      ) == 1L,
      paste(
        "Q7 pathway maps to multiple pathway_key values:",
        family_i,
        qi$library,
        qi$pathway
      )
    )
    
    
    pk <- pathway_keys[1L]
    
    
    # -------------------------------------------------------------------------
    # TX
    # -------------------------------------------------------------------------
    
    tx <- pathway_scores[
      pathway_scores$family ==
        family_i &
        pathway_scores$pathway_key ==
        pk,
      ,
      drop = FALSE
    ]
    
    
    if (nrow(tx)) {
      
      rows[[rr]] <- data.frame(
        
        family =
          family_i,
        
        program_id =
          as.character(
            qi$program_id
          ),
        
        program_short =
          as.character(
            qi$program_short
          ),
        
        factor_side =
          as.character(
            qi$factor_side
          ),
        
        phenotypic_side =
          as.character(
            qi$phenotypic_side
          ),
        
        library =
          as.character(
            qi$library
          ),
        
        pathway =
          as.character(
            qi$pathway
          ),
        
        pathway_label =
          as.character(
            qi$pathway_label
          ),
        
        pathway_key =
          pk,
        
        sample_id =
          as.character(
            tx$sample_id
          ),
        
        grupo =
          as.character(
            tx$grupo
          ),
        
        split =
          as.character(
            tx$split
          ),
        
        view_code =
          "tx",
        
        view_label =
          "TX",
        
        raw_score =
          as.numeric(
            tx$pathway_score
          ),
        
        n_features_used =
          NA_integer_,
        
        stringsAsFactors =
          FALSE
      )
      
      
      rr <- rr + 1L
    }
    
    
    # -------------------------------------------------------------------------
    # PR / ME / CL
    # -------------------------------------------------------------------------
    
    for (
      v in c(
        "pr",
        "me",
        "cl"
      )
    ) {
      
      zv <- zcat[
        tolower(
          as.character(
            zcat$view_code
          )
        ) == v,
        ,
        drop = FALSE
      ]
      
      
      if (
        "available_in_pooled" %in%
        names(zv)
      ) {
        
        zv <- zv[
          fig5_as_flag(
            zv$available_in_pooled
          ),
          ,
          drop = FALSE
        ]
      }
      
      
      zv <- zv[
        is.finite(
          suppressWarnings(
            as.numeric(
              zv$pathway_feature_cosine
            )
          )
        ),
        ,
        drop = FALSE
      ]
      
      
      if (!nrow(zv)) {
        next
      }
      
      
      features <- unique(
        as.character(
          zv$feature_original
        )
      )
      
      
      weights <- rep(
        NA_real_,
        length(
          features
        )
      )
      
      
      names(weights) <-
        features
      
      
      for (
        f in features
      ) {
        
        w <- unique(
          suppressWarnings(
            as.numeric(
              zv$pathway_feature_cosine[
                zv$feature_original ==
                  f
              ]
            )
          )
        )
        
        
        w <- w[
          is.finite(w)
        ]
        
        
        assert(
          length(w) <=
            1L,
          paste(
            "Q7 found non-unique crossomic weight for",
            family_i,
            qi$pathway,
            v,
            f
          )
        )
        
        
        if (length(w)) {
          weights[f] <- w
        }
      }
      
      
      X <- inp$matrices[[
        v
      ]]
      
      
      train_ids <- meta$sample_id[
        meta$split ==
          "TRAIN"
      ]
      
      
      sample_ids <-
        meta$sample_id
      
      
      usable <- features[
        features %in%
          rownames(X) &
          is.finite(
            weights[
              features
            ]
          ) &
          abs(
            weights[
              features
            ]
          ) >
          0
      ]
      
      
      if (!length(usable)) {
        next
      }
      
      
      Z <- matrix(
        NA_real_,
        nrow =
          length(usable),
        ncol =
          length(sample_ids),
        dimnames =
          list(
            usable,
            sample_ids
          )
      )
      
      
      keep_feature <- rep(
        FALSE,
        length(usable)
      )
      
      
      for (
        j in seq_along(
          usable
        )
      ) {
        
        f <- usable[j]
        
        
        train_values <- as.numeric(
          X[
            f,
            train_ids
          ]
        )
        
        
        train_values <- train_values[
          is.finite(
            train_values
          )
        ]
        
        
        if (
          length(
            train_values
          ) <
          2L
        ) {
          next
        }
        
        
        mu <- mean(
          train_values
        )
        
        
        s <- stats::sd(
          train_values
        )
        
        
        if (
          !is.finite(s) ||
          s <= 0
        ) {
          next
        }
        
        
        values <- as.numeric(
          X[
            f,
            sample_ids
          ]
        )
        
        
        Z[j, ] <-
          (
            values -
              mu
          ) /
          s
        
        
        keep_feature[j] <-
          TRUE
      }
      
      
      Z <- Z[
        keep_feature,
        ,
        drop = FALSE
      ]
      
      
      if (!nrow(Z)) {
        next
      }
      
      
      w <- weights[
        rownames(Z)
      ]
      
      
      score <- rep(
        NA_real_,
        ncol(Z)
      )
      
      
      for (
        j in seq_len(
          ncol(Z)
        )
      ) {
        
        values <- Z[
          ,
          j
        ]
        
        
        ok <-
          is.finite(values) &
          is.finite(w)
        
        
        denom <- sum(
          abs(
            w[ok]
          )
        )
        
        
        if (
          any(ok) &&
          is.finite(denom) &&
          denom > 0
        ) {
          
          score[j] <-
            sum(
              w[ok] *
                values[ok]
            ) /
            denom
        }
      }
      
      
      rows[[rr]] <- data.frame(
        
        family =
          family_i,
        
        program_id =
          as.character(
            qi$program_id
          ),
        
        program_short =
          as.character(
            qi$program_short
          ),
        
        factor_side =
          as.character(
            qi$factor_side
          ),
        
        phenotypic_side =
          as.character(
            qi$phenotypic_side
          ),
        
        library =
          as.character(
            qi$library
          ),
        
        pathway =
          as.character(
            qi$pathway
          ),
        
        pathway_label =
          as.character(
            qi$pathway_label
          ),
        
        pathway_key =
          pk,
        
        sample_id =
          sample_ids,
        
        grupo =
          as.character(
            meta$grupo
          ),
        
        split =
          meta$split,
        
        view_code =
          v,
        
        view_label =
          fig5_view_label(
            v
          ),
        
        raw_score =
          score,
        
        n_features_used =
          nrow(Z),
        
        stringsAsFactors =
          FALSE
      )
      
      
      rr <- rr + 1L
    }
  }
  
  
  if (!length(rows)) {
    return(
      data.frame()
    )
  }
  
  
  out <- do.call(
    rbind,
    rows
  )
  
  
  out <-
    q7_train_standardize_long(
      
      dat =
        out,
      
      value_col =
        "raw_score",
      
      key_cols =
        c(
          "family",
          "pathway_key",
          "view_code"
        ),
      
      output_col =
        "manifestation_z_train"
    )
  
  
  out$manifestation_z_plot <-
    q7_clip_z(
      out$manifestation_z_train
    )
  
  
  out
}


# ===== 10.6 FIGURA 7B — MANIFESTACIÓN MULTIMODAL =====

make_q7_multimodal_heatmaps <- function(
    manifestation,
    outdir
) {
  
  if (!nrow(manifestation)) {
    return(
      invisible(NULL)
    )
  }
  
  
  program_keys <- unique(
    manifestation[
      ,
      c(
        "family",
        "factor_side",
        "phenotypic_side",
        "program_id",
        "program_short"
      ),
      drop = FALSE
    ]
  )
  
  
  for (
    i in seq_len(
      nrow(
        program_keys
      )
    )
  ) {
    
    k <- program_keys[
      i,
      ,
      drop = FALSE
    ]
    
    
    z <- manifestation[
      manifestation$family ==
        k$family &
        manifestation$program_id ==
        k$program_id,
      ,
      drop = FALSE
    ]
    
    
    if (!nrow(z)) {
      next
    }
    
    
    sample_meta <-
      q7_order_samples(
        z
      )
    
    
    pathway_meta <- unique(
      z[
        ,
        c(
          "pathway_key",
          "library",
          "pathway_label"
        ),
        drop = FALSE
      ]
    )
    
    
    pathway_meta <- pathway_meta[
      order(
        pathway_meta$library,
        pathway_meta$pathway_label
      ),
      ,
      drop = FALSE
    ]
    
    
    pathway_meta$pathway_display <- make.unique(
      paste0(
        pathway_meta$pathway_label,
        "\n[",
        pathway_meta$library,
        "]"
      )
    )
    
    
    grid <- expand.grid(
      
      sample_id =
        sample_meta$sample_id,
      
      pathway_key =
        pathway_meta$pathway_key,
      
      view_code =
        c(
          "tx",
          "pr",
          "me",
          "cl"
        ),
      
      stringsAsFactors =
        FALSE
    )
    
    
    grid <- merge(
      grid,
      z[
        ,
        c(
          "sample_id",
          "pathway_key",
          "view_code",
          "manifestation_z_train",
          "manifestation_z_plot",
          "n_features_used"
        ),
        drop = FALSE
      ],
      by = c(
        "sample_id",
        "pathway_key",
        "view_code"
      ),
      all.x = TRUE,
      sort = FALSE
    )
    
    
    grid <- merge(
      grid,
      sample_meta[
        ,
        c(
          "sample_id",
          "sample_label",
          "grupo",
          "split"
        ),
        drop = FALSE
      ],
      by =
        "sample_id",
      all.x =
        TRUE,
      sort =
        FALSE
    )
    
    
    grid <- merge(
      grid,
      pathway_meta[
        ,
        c(
          "pathway_key",
          "pathway_display"
        ),
        drop = FALSE
      ],
      by =
        "pathway_key",
      all.x =
        TRUE,
      sort =
        FALSE
    )
    
    
    grid$sample_plot <- factor(
      grid$sample_id,
      levels =
        rev(
          sample_meta$sample_id
        ),
      labels =
        rev(
          sample_meta$sample_label
        )
    )
    
    
    grid$view_label <- factor(
      fig5_view_label(
        grid$view_code
      ),
      levels = c(
        "TX",
        "PR",
        "ME",
        "CL"
      )
    )
    
    
    grid$pathway_display <- factor(
      grid$pathway_display,
      levels =
        pathway_meta$pathway_display
    )
    
    
    p <- ggplot2::ggplot(
      grid,
      ggplot2::aes(
        x = view_label,
        y = sample_plot,
        fill = manifestation_z_plot
      )
    ) +
      
      ggplot2::geom_tile(
        color = "grey85",
        linewidth = 0.35
      ) +
      
      ggplot2::scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = c(
          -FIG7_Z_LIMIT,
          FIG7_Z_LIMIT
        ),
        na.value = "grey92",
        name = "Individual\nZ(TRAIN)"
      ) +
      
      ggplot2::facet_grid(
        . ~ pathway_display,
        scales = "free_x",
        space = "free_x"
      ) +
      
      ggplot2::labs(
        
        title = paste0(
          k$family,
          " — ",
          phenotype_label(
            k$phenotypic_side
          ),
          " — ",
          k$program_short
        ),
        
        subtitle = paste0(
          "Manifestación individual de pathways multimodalmente reforzados. ",
          "TX = score congelado; PR/ME/CL = score ponderado crossómico."
        ),
        
        x = NULL,
        y = NULL,
        
        caption = paste0(
          "Todos los escalados utilizan exclusivamente TRAIN. ",
          "Gris = vista sin link utilizable. ",
          "Estos scores no redefinen grupos, pathways ni individuos."
        )
      ) +
      
      fig_theme(10) +
      
      ggplot2::theme(
        
        panel.grid =
          ggplot2::element_blank(),
        
        axis.text.x =
          ggplot2::element_text(
            size = 9,
            face = "bold"
          ),
        
        axis.text.y =
          ggplot2::element_text(
            size = 7
          ),
        
        strip.text.x =
          ggplot2::element_text(
            size = 8
          ),
        
        legend.position =
          "bottom"
      )
    
    
    dir_fig <- file.path(
      outdir,
      "figures",
      "FIG7_INDIVIDUAL_MANIFESTATION",
      "02_MULTIMODAL_MANIFESTATION",
      as.character(
        k$family
      ),
      safe_filename(
        k$factor_side
      ),
      safe_filename(
        k$program_short
      )
    )
    
    
    dir_dat <- file.path(
      outdir,
      "plot_data",
      "FIG7_INDIVIDUAL_MANIFESTATION",
      "02_MULTIMODAL_MANIFESTATION",
      as.character(
        k$family
      ),
      safe_filename(
        k$factor_side
      ),
      safe_filename(
        k$program_short
      )
    )
    
    
    dir.create(
      dir_fig,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    
    dir.create(
      dir_dat,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    
    ggplot2::ggsave(
      
      filename = file.path(
        dir_fig,
        "FIG7B_INDIVIDUAL_MULTIMODAL_MANIFESTATION.png"
      ),
      
      plot = p,
      
      width = max(
        10,
        min(
          24,
          4 +
            2.2 *
            nrow(
              pathway_meta
            )
        )
      ),
      
      height = 10,
      
      dpi = FIG_DPI,
      
      limitsize = FALSE
    )
    
    
    write_csv(
      grid,
      file.path(
        dir_dat,
        "FIG7B_INDIVIDUAL_MULTIMODAL_MANIFESTATION.csv"
      )
    )
  }
  
  
  invisible(NULL)
}


# ===== 10.7 FIGURA 7C — CONCORDANCIA TX ↔ PR/ME/CL =====

make_q7_individual_concordance <- function(
    manifestation,
    outdir
) {
  
  if (!nrow(manifestation)) {
    return(
      invisible(NULL)
    )
  }
  
  
  pathways <- unique(
    manifestation[
      ,
      c(
        "family",
        "factor_side",
        "phenotypic_side",
        "program_id",
        "program_short",
        "library",
        "pathway",
        "pathway_label",
        "pathway_key"
      ),
      drop = FALSE
    ]
  )
  
  
  for (
    i in seq_len(
      nrow(pathways)
    )
  ) {
    
    pmeta <- pathways[
      i,
      ,
      drop = FALSE
    ]
    
    
    z <- manifestation[
      manifestation$family ==
        pmeta$family &
        manifestation$pathway_key ==
        pmeta$pathway_key,
      ,
      drop = FALSE
    ]
    
    
    tx <- z[
      z$view_code ==
        "tx",
      c(
        "sample_id",
        "manifestation_z_train"
      ),
      drop = FALSE
    ]
    
    
    if (!nrow(tx)) {
      next
    }
    
    
    names(tx)[
      names(tx) ==
        "manifestation_z_train"
    ] <- "tx_z"
    
    
    non_tx <- z[
      z$view_code %in%
        c(
          "pr",
          "me",
          "cl"
        ),
      ,
      drop = FALSE
    ]
    
    
    if (!nrow(non_tx)) {
      next
    }
    
    
    pairs <- merge(
      non_tx,
      tx,
      by =
        "sample_id",
      all.x =
        TRUE,
      sort =
        FALSE
    )
    
    
    pairs <- pairs[
      is.finite(
        pairs$tx_z
      ) &
        is.finite(
          pairs$manifestation_z_train
        ),
      ,
      drop = FALSE
    ]
    
    
    if (!nrow(pairs)) {
      next
    }
    
    
    pairs$view_label <- factor(
      pairs$view_label,
      levels = c(
        "PR",
        "ME",
        "CL"
      )
    )
    
    
    pairs$split <- factor(
      pairs$split,
      levels = c(
        "TRAIN",
        "TEST"
      )
    )
    
    
    p <- ggplot2::ggplot(
      pairs,
      ggplot2::aes(
        x = tx_z,
        y = manifestation_z_train,
        color = grupo,
        shape = split
      )
    ) +
      
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = 2,
        linewidth = 0.4
      ) +
      
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = 2,
        linewidth = 0.4
      ) +
      
      ggplot2::geom_point(
        size = 3.2,
        stroke = 1
      ) +
      
      ggplot2::scale_shape_manual(
        values = c(
          "TRAIN" = 16,
          "TEST" = 1
        )
      ) +
      
      ggplot2::facet_wrap(
        ~ view_label
      ) +
      
      ggplot2::labs(
        
        title = paste0(
          pmeta$family,
          " — ",
          phenotype_label(
            pmeta$phenotypic_side
          ),
          " — ",
          pmeta$pathway_label
        ),
        
        subtitle =
          "Cada punto = individuo; sólido = TRAIN; abierto = TEST.",
        
        x =
          "TX pathway score — Z(TRAIN)",
        
        y =
          "Manifestación de la vista — Z(TRAIN)",
        
        color =
          "Grupo",
        
        shape =
          "Split",
        
        caption = paste0(
          "Figura descriptiva. No se ajusta regresión y no se utiliza ",
          "para seleccionar individuos ni redefinir el pathway."
        )
      ) +
      
      fig_theme(10) +
      
      ggplot2::theme(
        legend.position =
          "bottom"
      )
    
    
    if (
      FIG7_SHOW_SAMPLE_LABELS
    ) {
      
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(
            label =
              sample_id
          ),
          size = 2.3,
          vjust = -0.8,
          show.legend = FALSE,
          check_overlap = TRUE
        )
    }
    
    
    pathway_folder <- paste0(
      safe_filename(
        pmeta$library
      ),
      "__",
      safe_filename(
        clean_pathway_label(
          pmeta$pathway
        )
      )
    )
    
    
    dir_fig <- file.path(
      outdir,
      "figures",
      "FIG7_INDIVIDUAL_MANIFESTATION",
      "03_TX_NON_TX_CONCORDANCE",
      as.character(
        pmeta$family
      ),
      safe_filename(
        pmeta$factor_side
      ),
      safe_filename(
        pmeta$program_short
      ),
      pathway_folder
    )
    
    
    dir_dat <- file.path(
      outdir,
      "plot_data",
      "FIG7_INDIVIDUAL_MANIFESTATION",
      "03_TX_NON_TX_CONCORDANCE",
      as.character(
        pmeta$family
      ),
      safe_filename(
        pmeta$factor_side
      ),
      safe_filename(
        pmeta$program_short
      ),
      pathway_folder
    )
    
    
    dir.create(
      dir_fig,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    
    dir.create(
      dir_dat,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    
    ggplot2::ggsave(
      
      filename = file.path(
        dir_fig,
        "FIG7C_INDIVIDUAL_TX_NON_TX_CONCORDANCE.png"
      ),
      
      plot = p,
      
      width = 11,
      height = 6.5,
      
      dpi = FIG_DPI,
      
      limitsize = FALSE
    )
    
    
    write_csv(
      pairs,
      file.path(
        dir_dat,
        "FIG7C_INDIVIDUAL_TX_NON_TX_CONCORDANCE.csv"
      )
    )
  }
  
  
  invisible(NULL)
}

# ===== 10.8 FIGURA 7D — SEPARACIÓN FENOTÍPICA POR PATHWAY =====
#
# Objetivo:
#   resumir pathway por pathway cuánto difieren NP / FD / SO.
#
# Se calculan TRES scopes:
#   ALL   = TRAIN + TEST (26 individuos): descriptivo principal.
#   TRAIN = auditoría interna.
#   TEST  = auditoría descriptiva externa; NO selecciona pathways.
#
# TX:
#   eta2 global + Hedges g pairwise.
#
# MULTIMODAL:
#   pseudo-F y R2 tipo PERMANOVA usando TX/PR/ME/CL disponibles.
#
# Importante:
#   "Higher/Lower phenotypic separation" es una clasificación RELATIVA
#   dentro de cada P1/P2/P3, no una nueva definición biológica.
# =============================================================================


q7_hedges_g <- function(
    x1,
    x2
) {
  
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  
  x1 <- x1[
    is.finite(x1)
  ]
  
  x2 <- x2[
    is.finite(x2)
  ]
  
  
  n1 <- length(x1)
  n2 <- length(x2)
  
  
  if (
    n1 < 2L ||
    n2 < 2L
  ) {
    return(NA_real_)
  }
  
  
  s1 <- stats::var(x1)
  s2 <- stats::var(x2)
  
  df <- n1 + n2 - 2L
  
  
  if (
    !is.finite(s1) ||
    !is.finite(s2) ||
    df <= 0L
  ) {
    return(NA_real_)
  }
  
  
  sp2 <- (
    (n1 - 1L) * s1 +
      (n2 - 1L) * s2
  ) / df
  
  
  if (
    !is.finite(sp2) ||
    sp2 <= 0
  ) {
    return(NA_real_)
  }
  
  
  d <- (
    mean(x1) -
      mean(x2)
  ) / sqrt(sp2)
  
  
  J <- 1 -
    3 / (
      4 * (n1 + n2) -
        9
    )
  
  
  J * d
}


q7_scalar_separation_stats <- function(
    dat,
    scope = "ALL"
) {
  
  z <- dat
  
  
  if (
    scope != "ALL"
  ) {
    
    z <- z[
      z$split ==
        scope,
      ,
      drop = FALSE
    ]
  }
  
  
  z <- z[
    is.finite(
      z$pathway_z_train
    ) &
      !is.na(
        z$grupo
      ),
    ,
    drop = FALSE
  ]
  
  
  groups <- c(
    "NP",
    "FD",
    "SO"
  )
  
  
  out <- list(
    
    n_total =
      nrow(z),
    
    n_NP =
      sum(z$grupo == "NP"),
    
    n_FD =
      sum(z$grupo == "FD"),
    
    n_SO =
      sum(z$grupo == "SO"),
    
    mean_NP =
      NA_real_,
    
    mean_FD =
      NA_real_,
    
    mean_SO =
      NA_real_,
    
    eta2 =
      NA_real_,
    
    anova_F =
      NA_real_,
    
    anova_p =
      NA_real_,
    
    welch_p =
      NA_real_,
    
    g_NP_FD =
      NA_real_,
    
    g_NP_SO =
      NA_real_,
    
    g_FD_SO =
      NA_real_,
    
    max_abs_g =
      NA_real_,
    
    dominant_contrast =
      NA_character_
  )
  
  
  if (!nrow(z)) {
    return(out)
  }
  
  
  for (
    g in groups
  ) {
    
    v <- z$pathway_z_train[
      z$grupo ==
        g
    ]
    
    v <- v[
      is.finite(v)
    ]
    
    
    out[[
      paste0(
        "mean_",
        g
      )
    ]] <- if (
      length(v)
    ) {
      
      mean(v)
      
    } else {
      
      NA_real_
    }
  }
  
  
  present <- groups[
    vapply(
      groups,
      function(g) {
        sum(
          z$grupo ==
            g
        ) > 0L
      },
      logical(1)
    )
  ]
  
  
  if (
    length(present) >=
    2L
  ) {
    
    y <- z$pathway_z_train
    
    grand <- mean(y)
    
    ss_total <- sum(
      (
        y -
          grand
      )^2
    )
    
    
    ss_between <- sum(
      vapply(
        present,
        
        function(g) {
          
          yg <- y[
            z$grupo ==
              g
          ]
          
          length(yg) *
            (
              mean(yg) -
                grand
            )^2
        },
        
        numeric(1)
      )
    )
    
    
    if (
      is.finite(ss_total) &&
      ss_total > 0
    ) {
      
      out$eta2 <-
        ss_between /
        ss_total
      
      
      df_b <-
        length(present) -
        1L
      
      df_w <-
        length(y) -
        length(present)
      
      ss_w <-
        ss_total -
        ss_between
      
      
      if (
        df_b > 0L &&
        df_w > 0L &&
        is.finite(ss_w) &&
        ss_w > 0
      ) {
        
        out$anova_F <-
          (
            ss_between /
              df_b
          ) /
          (
            ss_w /
              df_w
          )
        
        
        out$anova_p <-
          stats::pf(
            out$anova_F,
            df_b,
            df_w,
            lower.tail = FALSE
          )
      }
    }
    
    
    if (
      all(
        vapply(
          present,
          
          function(g) {
            sum(
              z$grupo ==
              g
            ) >= 2L
          },
          
          logical(1)
        )
      )
    ) {
      
      wt <- try(
        stats::oneway.test(
          pathway_z_train ~ grupo,
          data = z,
          var.equal = FALSE
        ),
        silent = TRUE
      )
      
      
      if (
        !inherits(
          wt,
          "try-error"
        )
      ) {
        
        out$welch_p <-
          wt$p.value
      }
    }
  }
  
  
  out$g_NP_FD <- q7_hedges_g(
    z$pathway_z_train[
      z$grupo ==
        "NP"
    ],
    z$pathway_z_train[
      z$grupo ==
        "FD"
    ]
  )
  
  
  out$g_NP_SO <- q7_hedges_g(
    z$pathway_z_train[
      z$grupo ==
        "NP"
    ],
    z$pathway_z_train[
      z$grupo ==
        "SO"
    ]
  )
  
  
  out$g_FD_SO <- q7_hedges_g(
    z$pathway_z_train[
      z$grupo ==
        "FD"
    ],
    z$pathway_z_train[
      z$grupo ==
        "SO"
    ]
  )
  
  
  gv <- c(
    "NP-FD" =
      out$g_NP_FD,
    "NP-SO" =
      out$g_NP_SO,
    "FD-SO" =
      out$g_FD_SO
  )
  
  
  finite_g <- is.finite(
    gv
  )
  
  
  if (
    any(finite_g)
  ) {
    
    out$max_abs_g <-
      max(
        abs(
          gv[
            finite_g
          ]
        )
      )
    
    
    out$dominant_contrast <-
      names(gv)[
        which.max(
          ifelse(
            finite_g,
            abs(gv),
            -Inf
          )
        )
      ]
  }
  
  
  out
}


q7_permanova_core <- function(
    X,
    group
) {
  
  X <- as.matrix(X)
  
  group <- as.character(
    group
  )
  
  
  ok <-
    stats::complete.cases(X) &
    !is.na(group)
  
  
  X <- X[
    ok,
    ,
    drop = FALSE
  ]
  
  group <- group[
    ok
  ]
  
  
  lev <- unique(
    group
  )
  
  
  n <- nrow(X)
  k <- length(lev)
  
  
  if (
    n <= k ||
    k < 2L ||
    ncol(X) < 1L
  ) {
    
    return(
      c(
        R2 = NA_real_,
        F = NA_real_
      )
    )
  }
  
  
  grand <- colMeans(X)
  
  
  ss_total <- sum(
    rowSums(
      (
        X -
          matrix(
            grand,
            nrow = n,
            ncol = ncol(X),
            byrow = TRUE
          )
      )^2
    )
  )
  
  
  if (
    !is.finite(ss_total) ||
    ss_total <= 0
  ) {
    
    return(
      c(
        R2 = NA_real_,
        F = NA_real_
      )
    )
  }
  
  
  ss_between <- 0
  
  
  for (
    g in lev
  ) {
    
    idx <- which(
      group ==
        g
    )
    
    
    centroid <- colMeans(
      X[
        idx,
        ,
        drop = FALSE
      ]
    )
    
    
    ss_between <-
      ss_between +
      length(idx) *
      sum(
        (
          centroid -
            grand
        )^2
      )
  }
  
  
  ss_within <-
    ss_total -
    ss_between
  
  
  df_b <- k - 1L
  df_w <- n - k
  
  
  Fv <- if (
    df_b > 0L &&
    df_w > 0L &&
    ss_within > 0
  ) {
    
    (
      ss_between /
        df_b
    ) /
      (
        ss_within /
          df_w
      )
    
  } else {
    
    NA_real_
  }
  
  
  c(
    R2 =
      ss_between /
      ss_total,
    
    F =
      Fv
  )
}


q7_dispersion_F <- function(
    X,
    group
) {
  
  X <- as.matrix(X)
  
  group <- as.character(
    group
  )
  
  
  ok <-
    stats::complete.cases(X) &
    !is.na(group)
  
  
  X <- X[
    ok,
    ,
    drop = FALSE
  ]
  
  group <- group[
    ok
  ]
  
  
  lev <- unique(
    group
  )
  
  
  if (
    length(lev) < 2L ||
    nrow(X) <=
    length(lev)
  ) {
    
    return(NA_real_)
  }
  
  
  d <- rep(
    NA_real_,
    nrow(X)
  )
  
  
  for (
    g in lev
  ) {
    
    idx <- which(
      group ==
        g
    )
    
    
    centroid <- colMeans(
      X[
        idx,
        ,
        drop = FALSE
      ]
    )
    
    
    d[idx] <- sqrt(
      rowSums(
        (
          X[
            idx,
            ,
            drop = FALSE
          ] -
            matrix(
              centroid,
              nrow =
                length(idx),
              ncol =
                ncol(X),
              byrow =
                TRUE
            )
        )^2
      )
    )
  }
  
  
  grand <- mean(d)
  
  
  ss_total <- sum(
    (
      d -
        grand
    )^2
  )
  
  
  ss_between <- sum(
    vapply(
      lev,
      
      function(g) {
        
        dg <- d[
          group ==
            g
        ]
        
        length(dg) *
          (
            mean(dg) -
              grand
          )^2
      },
      
      numeric(1)
    )
  )
  
  
  ss_within <-
    ss_total -
    ss_between
  
  
  df_b <-
    length(lev) -
    1L
  
  df_w <-
    length(d) -
    length(lev)
  
  
  if (
    df_b <= 0L ||
    df_w <= 0L ||
    !is.finite(ss_within) ||
    ss_within <= 0
  ) {
    
    return(NA_real_)
  }
  
  
  (
    ss_between /
      df_b
  ) /
    (
      ss_within /
        df_w
    )
}


q7_multimodal_separation_stats <- function(
    dat,
    scope = "ALL",
    nperm =
      FIG7_PERMANOVA_NPERM,
    seed =
      FIG7_PERMANOVA_SEED
) {
  
  z <- dat
  
  
  views_order <- c(
    "tx",
    "pr",
    "me",
    "cl"
  )
  
  
  views <- views_order[
    views_order %in%
      unique(
        as.character(
          z$view_code[
            is.finite(
              z$manifestation_z_train
            )
          ]
        )
      )
  ]
  
  
  out <- list(
    
    n_multimodal =
      0L,
    
    n_views =
      length(views),
    
    views_used =
      if (
        length(views)
      ) {
        
        paste(
          toupper(views),
          collapse = ";"
        )
        
      } else {
        
        NA_character_
      },
    
    permanova_R2 =
      NA_real_,
    
    permanova_F =
      NA_real_,
    
    permanova_p =
      NA_real_,
    
    permdisp_F =
      NA_real_,
    
    permdisp_p =
      NA_real_
  )
  
  
  if (
    length(views) < 2L ||
    !"tx" %in% views
  ) {
    
    return(out)
  }
  
  
  meta <- unique(
    z[
      ,
      c(
        "sample_id",
        "grupo",
        "split"
      ),
      drop = FALSE
    ]
  )
  
  
  if (
    scope != "ALL"
  ) {
    
    meta <- meta[
      meta$split ==
        scope,
      ,
      drop = FALSE
    ]
  }
  
  
  if (!nrow(meta)) {
    return(out)
  }
  
  
  X <- matrix(
    NA_real_,
    nrow =
      nrow(meta),
    ncol =
      length(views),
    dimnames =
      list(
        meta$sample_id,
        views
      )
  )
  
  
  for (
    v in views
  ) {
    
    zv <- z[
      z$view_code ==
        v &
        is.finite(
          z$manifestation_z_train
        ),
      c(
        "sample_id",
        "manifestation_z_train"
      ),
      drop = FALSE
    ]
    
    
    if (nrow(zv)) {
      
      vals <- tapply(
        zv$manifestation_z_train,
        zv$sample_id,
        mean,
        na.rm = TRUE
      )
      
      
      idx <- match(
        names(vals),
        rownames(X)
      )
      
      good <- !is.na(idx)
      
      
      X[
        idx[good],
        v
      ] <- as.numeric(
        vals[good]
      )
    }
  }
  
  
  cc <- stats::complete.cases(X)
  
  
  X <- X[
    cc,
    ,
    drop = FALSE
  ]
  
  meta <- meta[
    cc,
    ,
    drop = FALSE
  ]
  
  
  out$n_multimodal <-
    nrow(X)
  
  
  if (
    nrow(X) < 4L ||
    length(
      unique(
        meta$grupo
      )
    ) < 2L
  ) {
    
    return(out)
  }
  
  
  obs <- q7_permanova_core(
    X,
    meta$grupo
  )
  
  
  out$permanova_R2 <-
    unname(
      obs["R2"]
    )
  
  out$permanova_F <-
    unname(
      obs["F"]
    )
  
  
  out$permdisp_F <-
    q7_dispersion_F(
      X,
      meta$grupo
    )
  
  
  if (
    !is.finite(
      out$permanova_F
    ) ||
    nperm < 1L
  ) {
    
    return(out)
  }
  
  
  set.seed(seed)
  
  
  perm_F <- rep(
    NA_real_,
    nperm
  )
  
  perm_D <- rep(
    NA_real_,
    nperm
  )
  
  
  for (
    b in seq_len(
      nperm
    )
  ) {
    
    gp <- sample(
      meta$grupo,
      replace = FALSE
    )
    
    
    perm_F[b] <-
      q7_permanova_core(
        X,
        gp
      )["F"]
    
    
    perm_D[b] <-
      q7_dispersion_F(
        X,
        gp
      )
  }
  
  
  n_valid_F <- sum(
    is.finite(
      perm_F
    )
  )
  
  
  if (
    n_valid_F > 0L
  ) {
    
    out$permanova_p <-
      (
        1 +
          sum(
            perm_F >=
              out$permanova_F,
            na.rm = TRUE
          )
      ) /
      (
        1 +
          n_valid_F
      )
  }
  
  
  n_valid_D <- sum(
    is.finite(
      perm_D
    )
  )
  
  
  if (
    is.finite(
      out$permdisp_F
    ) &&
    n_valid_D > 0L
  ) {
    
    out$permdisp_p <-
      (
        1 +
          sum(
            perm_D >=
              out$permdisp_F,
            na.rm = TRUE
          )
      ) /
      (
        1 +
          n_valid_D
      )
  }
  
  
  out
}


build_q7d_pathway_separation_master <- function(
    pathway_scores,
    manifestation,
    q6_pathways
) {
  
  keys <- unique(
    pathway_scores[
      ,
      c(
        "family",
        "program_id",
        "program_short",
        "factor_side",
        "phenotypic_side",
        "library",
        "pathway",
        "pathway_label",
        "pathway_key"
      ),
      drop = FALSE
    ]
  )
  
  
  rows <- vector(
    "list",
    nrow(keys)
  )
  
  
  scopes <- c(
    "ALL",
    "TRAIN",
    "TEST"
  )
  
  
  for (
    i in seq_len(
      nrow(keys)
    )
  ) {
    
    k <- keys[
      i,
      ,
      drop = FALSE
    ]
    
    
    tx <- pathway_scores[
      pathway_scores$family ==
        k$family &
        pathway_scores$program_id ==
        k$program_id &
        pathway_scores$pathway_key ==
        k$pathway_key,
      ,
      drop = FALSE
    ]
    
    
    mm <- manifestation[
      manifestation$family ==
        k$family &
        manifestation$program_id ==
        k$program_id &
        manifestation$pathway_key ==
        k$pathway_key,
      ,
      drop = FALSE
    ]
    
    
    row <- as.list(
      k[
        1L,
        ,
        drop = FALSE
      ]
    )
    
    
    for (
      scope in scopes
    ) {
      
      s <-
        q7_scalar_separation_stats(
          tx,
          scope
        )
      
      
      for (
        nm in names(s)
      ) {
        
        row[[
          paste0(
            nm,
            "_",
            scope
          )
        ]] <- s[[nm]]
      }
      
      
      m <-
        q7_multimodal_separation_stats(
          
          dat =
            mm,
          
          scope =
            scope,
          
          nperm =
            FIG7_PERMANOVA_NPERM,
          
          seed =
            FIG7_PERMANOVA_SEED +
            i * 17L +
            match(
              scope,
              scopes
            )
        )
      
      
      for (
        nm in names(m)
      ) {
        
        row[[
          paste0(
            nm,
            "_",
            scope
          )
        ]] <- m[[nm]]
      }
    }
    
    
    rows[[i]] <-
      as.data.frame(
        row,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
  }
  
  
  out <- do.call(
    rbind,
    rows
  )
  
  
  # ---------------------------------------------------------------------------
  # Añadir contexto de Q6.
  # ---------------------------------------------------------------------------
  ann_cols <- intersect(
    c(
      "family",
      "program_id",
      "library",
      "pathway",
      "phenotypic_side",
      
      # Conservación de la vía entre P1/P2/P3
      "analyses_with_pathway",
      "n_analyses_with_pathway",
      "pathway_analysis_scope",
      
      # Especificidad respecto al lado fenotípico
      "phenotypic_sides_with_pathway",
      "n_phenotypic_sides_with_pathway",
      "pathway_phenotype_scope",
      
      # Contexto conjunto vía ↔ programa
      "program_analysis_scope",
      "program_phenotype_scope",
      "pathway_program_joint_status",
      "strict_threeway_pathway_in_strict_program",
      
      # Refuerzo multimodal
      "multimodal_reinforcement_status",
      "n_nonTX_views_linked"
    ),
    names(
      q6_pathways
    )
  )
  
  if (
    length(ann_cols) >= 5L
  ) {
    
    ann <- unique(
      q6_pathways[
        ,
        ann_cols,
        drop = FALSE
      ]
    )
    
    
    merge_by <- intersect(
      c(
        "family",
        "program_id",
        "library",
        "pathway",
        "phenotypic_side"
      ),
      names(ann)
    )
    
    
    out <- merge(
      out,
      ann,
      by =
        merge_by,
      all.x =
        TRUE,
      sort =
        FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # FDR: se ajusta SEPARADAMENTE por P1/P2/P3 y por scope.
  # Los p/q son secundarios; el ranking descriptivo usa magnitud de efecto.
  # ---------------------------------------------------------------------------
  
  for (
    scope in scopes
  ) {
    
    out[[
      paste0(
        "welch_q_",
        scope
      )
    ]] <- NA_real_
    
    
    out[[
      paste0(
        "permanova_q_",
        scope
      )
    ]] <- NA_real_
    
    
    for (
      fam in unique(
        out$family
      )
    ) {
      
      idx <- which(
        out$family ==
          fam
      )
      
      
      p1 <- suppressWarnings(
        as.numeric(
          out[[
            paste0(
              "welch_p_",
              scope
            )
          ]][idx]
        )
      )
      
      
      p2 <- suppressWarnings(
        as.numeric(
          out[[
            paste0(
              "permanova_p_",
              scope
            )
          ]][idx]
        )
      )
      
      
      ok1 <- is.finite(p1)
      ok2 <- is.finite(p2)
      
      
      if (any(ok1)) {
        
        out[[
          paste0(
            "welch_q_",
            scope
          )
        ]][
          idx[ok1]
        ] <- stats::p.adjust(
          p1[ok1],
          method = "BH"
        )
      }
      
      
      if (any(ok2)) {
        
        out[[
          paste0(
            "permanova_q_",
            scope
          )
        ]][
          idx[ok2]
        ] <- stats::p.adjust(
          p2[ok2],
          method = "BH"
        )
      }
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # Índice descriptivo principal usando TODAS las muestras.
  #
  # Promedio de:
  #   - eta2 TX en ALL
  #   - R2 PERMANOVA multimodal en ALL
  #
  # Si solo está disponible una de ambas, se conserva esa.
  # ---------------------------------------------------------------------------
  
  out$separation_index_ALL <- rowMeans(
    cbind(
      suppressWarnings(
        as.numeric(
          out$eta2_ALL
        )
      ),
      suppressWarnings(
        as.numeric(
          out$permanova_R2_ALL
        )
      )
    ),
    na.rm = TRUE
  )
  
  
  out$separation_index_ALL[
    !is.finite(
      out$separation_index_ALL
    )
  ] <- NA_real_
  
  
  # ---------------------------------------------------------------------------
  # Clasificación RELATIVA dentro de cada análisis.
  #
  # No se llama "common", porque common/shared ya tiene un significado
  # específico en Q6.
  # ---------------------------------------------------------------------------
  
  out$relative_separation_percentile_ALL <-
    NA_real_
  
  out$relative_separation_class_ALL <-
    NA_character_
  
  
  for (
    fam in unique(
      out$family
    )
  ) {
    
    idx <- which(
      out$family ==
        fam &
        is.finite(
          out$separation_index_ALL
        )
    )
    
    
    if (!length(idx)) {
      next
    }
    
    
    pct <- rank(
      out$separation_index_ALL[idx],
      ties.method = "average"
    ) /
      length(idx)
    
    
    out$relative_separation_percentile_ALL[
      idx
    ] <- pct
    
    
    out$relative_separation_class_ALL[
      idx
    ] <- ifelse(
      pct <= 1 / 3,
      "Lower phenotypic separation",
      ifelse(
        pct <= 2 / 3,
        "Intermediate phenotypic separation",
        "Higher phenotypic separation"
      )
    )
  }
  
  
  out <- out[
    order(
      out$family,
      -out$separation_index_ALL,
      out$phenotypic_side,
      out$program_short,
      out$pathway_label
    ),
    ,
    drop = FALSE
  ]
  
  
  out
}

make_q7d_separation_landscape <- function(
    master,
    outdir
) {
  
  # ---------------------------------------------------------------------------
  # Métricas mostradas.
  #
  # ALL es la descripción principal con los 26 individuos.
  # TRAIN / TEST permanecen como auditoría.
  # ---------------------------------------------------------------------------
  
  metrics <- c(
    
    eta2_ALL =
      "TX eta2\nALL",
    
    permanova_R2_ALL =
      "Multimodal R2\nALL",
    
    eta2_TRAIN =
      "TX eta2\nTRAIN",
    
    permanova_R2_TRAIN =
      "Multimodal R2\nTRAIN",
    
    eta2_TEST =
      "TX eta2\nTEST",
    
    permanova_R2_TEST =
      "Multimodal R2\nTEST"
  )
  
  
  # ---------------------------------------------------------------------------
  # Orden BIOLÓGICO de conservación.
  #
  # Esta clasificación NO se recalcula aquí.
  # Procede directamente de Q6.
  # ---------------------------------------------------------------------------
  
  conservation_levels <- c(
    "Shared P1-P2-P3",
    "Shared two analyses",
    "Analysis-specific"
  )
  
  
  conservation_labels <- c(
    
    "Shared P1-P2-P3" =
      "Compartida P1-P2-P3",
    
    "Shared two analyses" =
      "Compartida en dos análisis",
    
    "Analysis-specific" =
      "Específica del análisis"
  )
  
  
  for (
    fam in unique(
      master$family
    )
  ) {
    
    z <- master[
      master$family ==
        fam,
      ,
      drop = FALSE
    ]
    
    
    if (!nrow(z)) {
      next
    }
    
    
    # -------------------------------------------------------------------------
    # Validar que Q6 haya llegado correctamente a Q7D.
    # -------------------------------------------------------------------------
    
    assert(
      "pathway_analysis_scope" %in%
        names(z),
      paste0(
        "FIG7D requires pathway_analysis_scope from Q6. ",
        "Missing in ",
        fam
      )
    )
    
    
    # -------------------------------------------------------------------------
    # Clasificación de conservación congelada en Q6.
    # -------------------------------------------------------------------------
    
    z$conservation_scope <-
      as.character(
        z$pathway_analysis_scope
      )
    
    
    z$conservation_scope[
      is.na(
        z$conservation_scope
      ) |
        !nzchar(
          z$conservation_scope
        )
    ] <- "Unknown"
    
    
    # Mantener primero las tres categorías esperadas.
    # Cualquier estado inesperado queda al final y no se pierde.
    
    extra_conservation <- setdiff(
      unique(
        z$conservation_scope
      ),
      conservation_levels
    )
    
    
    conservation_order <- c(
      conservation_levels,
      sort(
        extra_conservation
      )
    )
    
    
    z$conservation_scope <- factor(
      z$conservation_scope,
      levels =
        conservation_order
    )
    
    
    # Etiqueta legible.
    
    z$conservation_label <- as.character(
      z$conservation_scope
    )
    
    
    known <- z$conservation_label %in%
      names(
        conservation_labels
      )
    
    
    z$conservation_label[
      known
    ] <- conservation_labels[
      z$conservation_label[
        known
      ]
    ]
    
    
    conservation_label_order <- c(
      unname(
        conservation_labels[
          conservation_levels[
            conservation_levels %in%
              names(
                conservation_labels
              )
          ]
        ]
      ),
      setdiff(
        unique(
          z$conservation_label
        ),
        unname(
          conservation_labels
        )
      )
    )
    
    
    z$conservation_label <- factor(
      z$conservation_label,
      levels =
        conservation_label_order
    )
    
    
    # -------------------------------------------------------------------------
    # Lado fenotípico.
    # -------------------------------------------------------------------------
    
    z$side_label <-
      phenotype_label(
        z$phenotypic_side
      )
    
    
    preferred_side_order <- c(
      "NP-side",
      "FD-side",
      "SO-side",
      "non-FD-side",
      "ambiguous-side"
    )
    
    
    side_levels <- c(
      preferred_side_order[
        preferred_side_order %in%
          unique(
            z$side_label
          )
      ],
      setdiff(
        unique(
          z$side_label
        ),
        preferred_side_order
      )
    )
    
    
    z$side_label <- factor(
      z$side_label,
      levels =
        side_levels
    )
    
    
    # -------------------------------------------------------------------------
    # Orden interno:
    #
    # 1. shared P1/P2/P3
    # 2. shared two analyses
    # 3. analysis-specific
    #
    # Dentro de cada categoría:
    #   - lado fenotípico
    #   - mayor separación ALL primero
    #   - programa / pathway para desempatar
    # -------------------------------------------------------------------------
    
    z <- z[
      order(
        z$conservation_scope,
        z$side_label,
        -z$separation_index_ALL,
        z$program_short,
        z$pathway_label
      ),
      ,
      drop = FALSE
    ]
    
    
    # -------------------------------------------------------------------------
    # Etiqueta de fila.
    #
    # No repetimos shared/specific en cada fila porque ya aparece como bloque.
    # -------------------------------------------------------------------------
    
    z$row_label <- make.unique(
      paste0(
        z$program_short,
        " | ",
        z$pathway_label,
        " [",
        z$library,
        "]"
      )
    )
    
    
    # -------------------------------------------------------------------------
    # Tabla larga para heatmap.
    # -------------------------------------------------------------------------
    
    long <- do.call(
      rbind,
      
      lapply(
        names(metrics),
        
        function(m) {
          
          data.frame(
            
            row_label =
              z$row_label,
            
            conservation_label =
              z$conservation_label,
            
            side_label =
              z$side_label,
            
            metric =
              unname(
                metrics[m]
              ),
            
            value =
              suppressWarnings(
                as.numeric(
                  z[[m]]
                )
              ),
            
            stringsAsFactors =
              FALSE
          )
        }
      )
    )
    
    
    # -------------------------------------------------------------------------
    # Conservar exactamente el orden ya establecido.
    # -------------------------------------------------------------------------
    
    long$row_label <- factor(
      long$row_label,
      levels =
        rev(
          z$row_label
        )
    )
    
    
    long$conservation_label <- factor(
      long$conservation_label,
      levels =
        levels(
          z$conservation_label
        )
    )
    
    
    long$side_label <- factor(
      long$side_label,
      levels =
        levels(
          z$side_label
        )
    )
    
    
    long$metric <- factor(
      long$metric,
      levels =
        unname(
          metrics
        )
    )
    
    
    long$cell_label <- ifelse(
      is.finite(
        long$value
      ),
      sprintf(
        "%.2f",
        long$value
      ),
      ""
    )
    
    
    # -------------------------------------------------------------------------
    # Figura.
    #
    # FILAS:
    #   conservación Q6 × lado fenotípico.
    #
    # Dentro de cada bloque:
    #   pathways ordenados por separation_index_ALL.
    # -------------------------------------------------------------------------
    
    p <- ggplot2::ggplot(
      long,
      ggplot2::aes(
        x = metric,
        y = row_label,
        fill = value
      )
    ) +
      
      ggplot2::geom_tile(
        color = "grey85",
        linewidth = 0.3
      ) +
      
      ggplot2::geom_text(
        ggplot2::aes(
          label =
            cell_label
        ),
        size = 2.7
      ) +
      
      ggplot2::scale_fill_gradient(
        low = "white",
        high = "#08519C",
        limits = c(
          0,
          1
        ),
        na.value = "grey92",
        name = "Separación"
      ) +
      
      ggplot2::facet_grid(
        
        rows =
          ggplot2::vars(
            conservation_label,
            side_label
          ),
        
        scales =
          "free_y",
        
        space =
          "free_y",
        
        switch =
          "y"
      ) +
      
      ggplot2::labs(
        
        title = paste0(
          fam,
          " — separación fenotípica por pathway"
        ),
        
        subtitle = paste0(
          "Vías agrupadas según conservación P1/P2/P3 definida en Q6 ",
          "y ordenadas dentro de cada bloque por separación fenotípica ALL."
        ),
        
        x = NULL,
        y = NULL,
        
        caption = paste0(
          "ALL = TRAIN + TEST (26 individuos) y constituye la descripción principal. ",
          "TX eta2 = proporción descriptiva de variación del pathway score ",
          "asociada a NP/FD/SO. ",
          "Multimodal R2 = PERMANOVA euclídea sobre TX/PR/ME/CL disponibles. ",
          "La clasificación compartida/específica procede exclusivamente de Q6; ",
          "no se redefine con estas métricas. ",
          "TRAIN y TEST se muestran separadamente como auditoría."
        )
      ) +
      
      fig_theme(10) +
      
      ggplot2::theme(
        
        panel.grid =
          ggplot2::element_blank(),
        
        axis.text.x =
          ggplot2::element_text(
            size = 8.5,
            face = "bold"
          ),
        
        axis.text.y =
          ggplot2::element_text(
            size = 7
          ),
        
        strip.placement =
          "outside",
        
        strip.text.y.left =
          ggplot2::element_text(
            size = 8.5,
            face = "bold",
            angle = 0,
            hjust = 1
          ),
        
        legend.position =
          "bottom"
      )
    
    
    # -------------------------------------------------------------------------
    # Salidas.
    # -------------------------------------------------------------------------
    
    dir_fig <- file.path(
      outdir,
      "figures",
      "FIG7_INDIVIDUAL_MANIFESTATION",
      "04_PATHWAY_PHENOTYPE_SEPARATION"
    )
    
    
    dir_dat <- file.path(
      outdir,
      "plot_data",
      "FIG7_INDIVIDUAL_MANIFESTATION",
      "04_PATHWAY_PHENOTYPE_SEPARATION"
    )
    
    
    dir.create(
      dir_fig,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    
    dir.create(
      dir_dat,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    
    ggplot2::ggsave(
      
      filename = file.path(
        dir_fig,
        paste0(
          "FIG7D_",
          fam,
          "_PATHWAY_PHENOTYPE_SEPARATION.png"
        )
      ),
      
      plot = p,
      
      width = 13.5,
      
      height = max(
        9,
        min(
          34,
          4 +
            0.32 *
            nrow(z)
        )
      ),
      
      dpi =
        FIG_DPI,
      
      limitsize =
        FALSE
    )
    
    
    # -------------------------------------------------------------------------
    # CSV específico de la figura.
    #
    # Se ordena exactamente igual:
    # conservación -> lado -> separación.
    # -------------------------------------------------------------------------
    
    write_csv(
      z,
      file.path(
        dir_dat,
        paste0(
          "FIG7D_",
          fam,
          "_PATHWAY_PHENOTYPE_SEPARATION.csv"
        )
      )
    )
  }
  
  
  invisible(NULL)
}
# ===== 10.9 EJECUTAR PREGUNTA 7 =====

if (run_question(7L)) {
  
  message("")
  message("============================================================")
  message("METATEST — QUESTION 7")
  message("Individual biological and multimodal manifestation")
  message("============================================================")
  
  
  q7_root_dat <- file.path(
    OUT_ROOT,
    "plot_data",
    "FIG7_INDIVIDUAL_MANIFESTATION",
    "00_INDEX"
  )
  
  
  dir.create(
    q7_root_dat,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # ---------------------------------------------------------------------------
  # Q6 define qué pathways poseen refuerzo multimodal.
  # ---------------------------------------------------------------------------
  
  q6_pathway_master_path <- file.path(
    OUT_ROOT,
    "plot_data",
    "FIG6_COMMON_SPECIFIC_REINFORCED",
    "00_INDEX",
    "Q6_PATHWAY_MASTER.csv"
  )
  
  
  assert(
    file.exists(
      q6_pathway_master_path
    ),
    paste(
      "Question 7 requires Question 6 output:",
      q6_pathway_master_path
    )
  )
  
  
  q6_pathways <- read_csv_checked(
    q6_pathway_master_path
  )
  
  
  # ---------------------------------------------------------------------------
  # Q5: scores congelados y catálogo crossómico.
  # ---------------------------------------------------------------------------
  
  q5_index_dir <- file.path(
    OUT_ROOT,
    "plot_data",
    "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
    "00_INDEX"
  )
  
  
  pathway_scores_raw <- read_csv_checked(
    file.path(
      q5_index_dir,
      "FIG5_POOLED_FROZEN_PATHWAY_SCORES.csv"
    ),
    required = c(
      "family",
      "program_id",
      "program_short",
      "factor",
      "phenotypic_side",
      "factor_side",
      "library",
      "pathway",
      "pathway_label",
      "pathway_key",
      "sample_id",
      "grupo",
      "pathway_score"
    )
  )
  
  
  fig5_catalog_q7 <- read_csv_checked(
    file.path(
      q5_index_dir,
      "FIG5_MASTER_PATHWAY_FEATURE_CATALOG.csv"
    ),
    required = c(
      "family",
      "program_id",
      "program_short",
      "library",
      "pathway",
      "pathway_label",
      "pathway_key",
      "phenotypic_side",
      "factor_side",
      "view_code",
      "feature_original",
      "pathway_feature_cosine"
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Reconstruir únicamente las entradas necesarias.
  #
  # Esto NO vuelve a ejecutar Q5.
  # ---------------------------------------------------------------------------
  
  q7_inputs <- list()
  
  q7_meta_list <- list()
  
  
  for (
    family in FAMILIES
  ) {
    
    message(
      "[Q7 INPUT] ",
      family
    )
    
    
    q7_inputs[[family]] <-
      load_fig5_inputs(
        bio_program_run =
          BIO_RUN,
        family =
          family
      )
    
    
    q7_meta_list[[family]] <-
      load_q7_sample_meta(
        q7_inputs[[
          family
        ]]
      )
  }
  
  
  q7_sample_meta <- do.call(
    rbind,
    q7_meta_list
  )
  
  
  write_csv(
    q7_sample_meta,
    file.path(
      q7_root_dat,
      "Q7_SAMPLE_METADATA.csv"
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # FIG7A
  # ---------------------------------------------------------------------------
  
  pathway_scores_q7 <-
    prepare_q7_pathway_scores(
      pathway_scores =
        pathway_scores_raw,
      sample_meta =
        q7_sample_meta
    )
  
  
  write_csv(
    pathway_scores_q7,
    file.path(
      q7_root_dat,
      "Q7_INDIVIDUAL_PATHWAY_SCORES.csv"
    )
  )
  
  
  fingerprint_blocks <- unique(
    pathway_scores_q7[
      ,
      c(
        "family",
        "factor_side"
      ),
      drop = FALSE
    ]
  )
  
  
  for (
    i in seq_len(
      nrow(
        fingerprint_blocks
      )
    )
  ) {
    
    make_q7_pathway_fingerprint(
      
      pathway_scores =
        pathway_scores_q7,
      
      family =
        fingerprint_blocks$family[i],
      
      factor_side =
        fingerprint_blocks$factor_side[i],
      
      outdir =
        OUT_ROOT
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # FIG7B + FIG7C
  # ---------------------------------------------------------------------------
  
  manifestation_q7 <-
    build_q7_multimodal_manifestation(
      
      q6_pathways =
        q6_pathways,
      
      fig5_catalog =
        fig5_catalog_q7,
      
      pathway_scores =
        pathway_scores_q7,
      
      inputs_list =
        q7_inputs,
      
      sample_meta =
        q7_sample_meta
    )
  
  
  write_csv(
    manifestation_q7,
    file.path(
      q7_root_dat,
      "Q7_INDIVIDUAL_MULTIMODAL_MANIFESTATION.csv"
    )
  )
  
  
  make_q7_multimodal_heatmaps(
    manifestation =
      manifestation_q7,
    outdir =
      OUT_ROOT
  )
  
  
  make_q7_individual_concordance(
    manifestation =
      manifestation_q7,
    outdir =
      OUT_ROOT
  )
  # ---------------------------------------------------------------------------
  # FIG7D — separación fenotípica pathway por pathway
  #
  # ALL = 26 individuos y es el resumen descriptivo principal.
  # TRAIN / TEST se mantienen como auditorías separadas.
  # ---------------------------------------------------------------------------
  
  q7_separation_master <-
    build_q7d_pathway_separation_master(
      
      pathway_scores =
        pathway_scores_q7,
      
      manifestation =
        manifestation_q7,
      
      q6_pathways =
        q6_pathways
    )
  
  
  write_csv(
    q7_separation_master,
    file.path(
      q7_root_dat,
      "Q7_PATHWAY_PHENOTYPE_SEPARATION_MASTER.csv"
    )
  )
  
  
  make_q7d_separation_landscape(
    master =
      q7_separation_master,
    outdir =
      OUT_ROOT
  )
  
  # ---------------------------------------------------------------------------
  # Índice.
  # ---------------------------------------------------------------------------
  q7_output_index <- data.frame(
    
    folder = c(
      "01_PATHWAY_FINGERPRINT",
      "02_MULTIMODAL_MANIFESTATION",
      "03_TX_NON_TX_CONCORDANCE",
      "04_PATHWAY_PHENOTYPE_SEPARATION"
    ),
    
    question = c(
      
      paste0(
        "¿Cómo se manifiestan los pathway scores congelados ",
        "individuo por individuo?"
      ),
      
      paste0(
        "¿Cómo se manifiesta en cada individuo un pathway ",
        "multimodalmente reforzado a través de TX/PR/ME/CL?"
      ),
      
      paste0(
        "¿Qué individuos explican visualmente la concordancia ",
        "TX ↔ PR/ME/CL?"
      ),
      
      paste0(
        "¿Qué pathways muestran mayor o menor separación fenotípica ",
        "entre NP, FD y SO?"
      )
    ),
    
    interpretation = c(
      
      "Fingerprint individual pathway/programa",
      
      paste0(
        "Score individual por vista usando únicamente ",
        "links crossómicos congelados"
      ),
      
      paste0(
        "Scatter descriptivo individuo por individuo; ",
        "sin regresión ni selección"
      ),
      
      paste0(
        "ALL = resumen descriptivo principal con los 26 individuos; ",
        "TRAIN/TEST = auditoría. eta2 TX, Hedges g pairwise y ",
        "PERMANOVA R2 multimodal."
      )
    ),
    
    stringsAsFactors =
      FALSE
  )  
  write_csv(
    q7_output_index,
    file.path(
      q7_root_dat,
      "Q7_OUTPUT_INDEX.csv"
    )
  )
  
  
  cat(
    
    paste0(
      "\n",
      "PREGUNTA 7 — manifestación individual de la arquitectura biológica.\n",
      "  - Q7 no redefine pathways, programas, links ni grupos.\n",
      "  - FIG7A usa los pathway scores congelados de Q5.\n",
      "  - Los Z individuales se definen exclusivamente con media/SD de TRAIN.\n",
      "  - TEST recibe la transformación definida en TRAIN.\n",
      "  - FIG7B se limita a pathways multimodalmente reforzados en Q6.\n",
      "  - PR/ME/CL usan únicamente features y pathway_feature_cosine ",
      "congelados aguas arriba.\n",
      "  - Las rho descriptivas de Q5 no seleccionan features en Q7.\n",
      "  - FIG7C muestra individuos; no ajusta regresión ni define outliers.\n",
        "  - FIG7D resume la separación NP/FD/SO pathway por pathway.\n",
      "  - ALL combina TRAIN + TEST (26 individuos) y es la descripción principal.\n",
      "  - TRAIN y TEST se reportan además por separado como auditoría.\n",
      "  - TEST no selecciona pathways ni modifica ninguna definición congelada.\n",
      "  - eta2 y PERMANOVA R2 se interpretan prioritariamente como tamaños de efecto.\n",
      "  - Los p/q se exportan como información secundaria, no como criterio de ranking.\n",
      "  - Hedges g conserva el signo para NP-FD, NP-SO y FD-SO.\n",
      "  - Higher/Intermediate/Lower phenotypic separation es una clasificación relativa dentro de P1/P2/P3.\n"
    ),
    
    file =
      file.path(
        OUT_ROOT,
        "READ_ME_FIRST.txt"
      ),
    
    append =
      TRUE
  )
  
  
  message("")
  message("QUESTION 7 COMPLETE")
  
  message(
    "Output: ",
    file.path(
      OUT_ROOT,
      "figures",
      "FIG7_INDIVIDUAL_MANIFESTATION"
    )
  )
}
# ============================================================================
# MARCAR LA CORRIDA COMO COMPLETA
# ============================================================================

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


if (
  file.exists(
    file.path(
      OUT_ROOT,
      "RUNNING.flag"
    )
  )
) {
  
  unlink(
    file.path(
      OUT_ROOT,
      "RUNNING.flag"
    )
  )
}


writeLines(
  canonical_path(
    OUT_ROOT,
    must_work = TRUE
  ),
  file.path(
    latest_dir,
    "LATEST_RUN.txt"
  )
)


if (
  file.exists(
    file.path(
      latest_dir,
      "LATEST_INCOMPLETE_RUN.txt"
    )
  )
) {
  
  unlink(
    file.path(
      latest_dir,
      "LATEST_INCOMPLETE_RUN.txt"
    )
  )
}


message(
  "Output root: ",
  OUT_ROOT
)


#!/usr/bin/env Rscript
# =============================================================================
# METATEST — PASO 4
# Puente clase <-> factor <-> pathway
# Version 0.1.0
#
# OBJETIVO
#   Combinar, SIN reentrenar nada:
#
#     1) Auditoría bayesiana del espacio latente:
#          clase <-> factor
#
#     2) Enriquecimiento nativo MOFA:
#          factor <-> pathway
#
#   para obtener:
#
#          clase <-> factor <-> pathway
#
# PRINCIPIO
#   Para cada factor orientado, la auditoría latente ya estimó las localizaciones
#   posteriores de NP, FD y SO.
#
#   En cada draw posterior:
#
#       cara positiva  = clase con mayor localización latente
#       cara negativa  = clase con menor localización latente
#
#   Esto permite estimar:
#
#       P(cara positiva = NP/FD/SO | datos)
#       P(cara negativa = NP/FD/SO | datos)
#
#   Los pathways enriquecidos en una "oriented_side" del factor heredan esa
#   distribución posterior de clase.
#
#   IMPORTANTE:
#   - No se calcula un nuevo GSEA.
#   - No se recalculan W ni Z.
#   - No se usa TEST.
#   - No se usa el framework bidireccional para seleccionar pathways.
#   - No se construyen todavía módulos cross-ómicos.
#
# EVIDENCIA DE PATHWAY
#   Primaria:
#       q_both_sides_within_library_factor <= 0.10
#
#   Secundaria:
#       q_mofa_within_side <= 0.10
#       pero q_both_sides_within_library_factor > 0.10
#
#   Se conserva TODO en la tabla "ALL"; los thresholds solo etiquetan evidencia.
# =============================================================================

options(stringsAsFactors = FALSE, digits = 17, warn = 1)

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

PATHWAY_FDR <- as.numeric(
  Sys.getenv("CLASS_PATHWAY_FDR", unset = "0.10")
)

CLASS_ASSIGN_STRONG <- as.numeric(
  Sys.getenv("CLASS_PATHWAY_ASSIGN_STRONG", unset = "0.90")
)

TOP_PATHWAYS_PLOT <- as.integer(
  Sys.getenv("CLASS_PATHWAY_TOP_PLOT", unset = "15")
)

MAKE_PLOTS <- TRUE

assert <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
}

canonical_path <- function(p, must_work = TRUE) {
  normalizePath(path.expand(p), winslash = "/", mustWork = must_work)
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x) || !ncol(x)) x <- data.frame(status = character())
  utils::write.table(
    x, path,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    na = "NA",
    qmethod = "double",
    fileEncoding = "UTF-8"
  )
}

entropy_normalized <- function(p) {
  p <- as.numeric(p)
  p <- p[is.finite(p) & p > 0]
  if (!length(p)) return(NA_real_)
  if (length(p) == 1L) return(0)
  -sum(p * log(p)) / log(length(GROUPS))
}

find_latest_selection_run <- function(upstream) {
  base <- file.path(upstream, "downstream", "seleccion_mejor_escenario_train")
  assert(dir.exists(base), paste("No existe:", base))
  dirs <- list.dirs(base, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "SELECTED_SCENARIOS_BY_ANALYSIS.csv"))]
  assert(length(dirs) > 0L, paste("No hay corridas de selección completas en", base))
  info <- file.info(dirs)
  dirs[order(info$mtime, basename(dirs), decreasing = TRUE)][1L]
}

find_bidirectional_run <- function(selection_run, family) {
  root <- file.path(selection_run, family, "BIDIRECTIONAL_FRAMEWORK")
  assert(dir.exists(root), paste("No existe BIDIRECTIONAL_FRAMEWORK para", family))

  latest_file <- file.path(root, "LATEST_RUN.txt")
  if (file.exists(latest_file)) {
    p <- trimws(readLines(latest_file, warn = FALSE)[1L])
    if (
      nzchar(p) &&
      dir.exists(p) &&
      file.exists(file.path(p, "rds", "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"))
    ) {
      return(canonical_path(p, TRUE))
    }
  }

  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "rds", "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"))]
  assert(length(dirs) > 0L, paste("No hay corrida bidireccional completa para", family))
  info <- file.info(dirs)
  canonical_path(dirs[order(info$mtime, basename(dirs), decreasing = TRUE)][1L], TRUE)
}

find_latest_subrun <- function(parent, folder, required_file) {
  root <- file.path(parent, folder)
  assert(dir.exists(root), paste("No existe:", root))

  latest_file <- file.path(root, "LATEST_RUN.txt")
  if (file.exists(latest_file)) {
    p <- trimws(readLines(latest_file, warn = FALSE)[1L])
    if (
      nzchar(p) &&
      dir.exists(p) &&
      file.exists(file.path(p, required_file))
    ) {
      return(canonical_path(p, TRUE))
    }
  }

  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, required_file))]
  assert(
    length(dirs) > 0L,
    paste("No hay corrida completa de", folder, "con", required_file)
  )
  info <- file.info(dirs)
  canonical_path(dirs[order(info$mtime, basename(dirs), decreasing = TRUE)][1L], TRUE)
}

find_latent_audit_run <- function(bidir_run) {
  find_latest_subrun(
    parent = bidir_run,
    folder = "LATENT_SUPERVISED_AUDIT",
    required_file = file.path("tables", "05_TRAIN_BAYESIAN_FACTOR_CONTRASTS.csv")
  )
}

find_enrichment_run <- function(bidir_run) {
  find_latest_subrun(
    parent = bidir_run,
    folder = "MOFA_NATIVE_ENRICHMENT",
    required_file = file.path("tables", "03_MOFA_NATIVE_ENRICHMENT_ALL.csv")
  )
}

find_factor_fit <- function(latent_run, factor_name) {
  candidates <- c(
    file.path(latent_run, "rds", paste0("BRMS_", factor_name, ".rds")),
    file.path(latent_run, "rds", paste0("BRMS_", make.names(factor_name), ".rds"))
  )

  hit <- candidates[file.exists(candidates)]

  assert(
    length(hit) >= 1L,
    paste("No se encontró el fit BRMS para", factor_name, "en", file.path(latent_run, "rds"))
  )

  hit[1L]
}

extract_group_location_draws <- function(fit) {
  draws <- posterior::as_draws_df(fit)

  coef_names <- paste0("b_grupo", GROUPS)

  assert(
    all(coef_names %in% names(draws)),
    paste("El fit BRMS no contiene:", paste(coef_names, collapse = ", "))
  )

  out <- as.matrix(draws[, coef_names, drop = FALSE])
  colnames(out) <- GROUPS
  out
}

make_factor_side_class_posterior <- function(factor_name, group_draws) {
  assert(all(GROUPS %in% colnames(group_draws)), "Faltan grupos en los draws")
  X <- group_draws[, GROUPS, drop = FALSE]
  
  tol <- sqrt(.Machine$double.eps)
  
  max_idx <- apply(X, 1, function(x) {
    mx <- max(x)
    hit <- which(abs(x - mx) <= tol)
    
    if (length(hit) == 1L) hit else NA_integer_
  })
  
  min_idx <- apply(X, 1, function(x) {
    mn <- min(x)
    hit <- which(abs(x - mn) <= tol)
    
    if (length(hit) == 1L) hit else NA_integer_
  })
  probs_positive <- prop.table(
    table(
      factor(
        GROUPS[max_idx[!is.na(max_idx)]],
        levels = GROUPS
      )
    )
  )
  
  probs_negative <- prop.table(
    table(
      factor(
        GROUPS[min_idx[!is.na(min_idx)]],
        levels = GROUPS
      )
    )
  )
  build_side <- function(side, probs) {
    probs <- as.numeric(probs)
    names(probs) <- GROUPS

    ord <- order(probs, decreasing = TRUE)
    top_class <- GROUPS[ord[1L]]
    second_class <- GROUPS[ord[2L]]

    data.frame(
      factor = factor_name,
      oriented_side = side,
      p_side_NP = unname(probs["NP"]),
      p_side_FD = unname(probs["FD"]),
      p_side_SO = unname(probs["SO"]),
      assigned_class_nominal = top_class,
      p_assigned_class = unname(probs[top_class]),
      second_class = second_class,
      p_second_class = unname(probs[second_class]),
      probability_margin = unname(probs[top_class] - probs[second_class]),
      normalized_entropy = entropy_normalized(probs),
      assignment_strong = unname(probs[top_class]) >= CLASS_ASSIGN_STRONG,
      stringsAsFactors = FALSE
    )
  }

  rbind(
    build_side("positive", probs_positive),
    build_side("negative", probs_negative)
  )
}

make_pairwise_side_alignment <- function(factor_name, group_draws) {
  contrast_defs <- list(
    FD_minus_NP = c("FD", "NP"),
    SO_minus_NP = c("SO", "NP"),
    FD_minus_SO = c("FD", "SO")
  )

  rows <- list()
  rr <- 1L

  for (contrast_name in names(contrast_defs)) {
    cc <- contrast_defs[[contrast_name]]
    delta <- group_draws[, cc[1L]] - group_draws[, cc[2L]]

    for (side in c("positive", "negative")) {
      side_sign <- if (side == "positive") 1 else -1
      gamma <- side_sign * delta

      p_first <- mean(gamma > 0)
      p_second <- mean(gamma < 0)

      q <- stats::quantile(
        gamma,
        probs = c(0.025, 0.50, 0.975),
        names = FALSE
      )

      rows[[rr]] <- data.frame(
        factor = factor_name,
        oriented_side = side,
        contrast = contrast_name,
        first_class = cc[1L],
        second_class = cc[2L],
        alignment_posterior_mean = mean(gamma),
        alignment_q025 = q[1L],
        alignment_q50 = q[2L],
        alignment_q975 = q[3L],
        p_side_toward_first_class = p_first,
        p_side_toward_second_class = p_second,
        nominal_pairwise_class = if (p_first >= p_second) cc[1L] else cc[2L],
        p_nominal_pairwise_class = max(p_first, p_second),
        stringsAsFactors = FALSE
      )

      rr <- rr + 1L
    }
  }

  do.call(rbind, rows)
}

classify_pathway_evidence <- function(q_native, q_two_sides) {
  if (is.finite(q_two_sides) && q_two_sides <= PATHWAY_FDR) {
    return("PRIMARY_TWO_SIDE_FDR")
  }

  if (is.finite(q_native) && q_native <= PATHWAY_FDR) {
    return("NATIVE_WITHIN_SIDE_ONLY")
  }

  "NOT_FDR_SUPPORTED"
}

build_class_pathway_bridge <- function(enrichment, side_class) {
  needed_enrichment <- c(
    "factor",
    "oriented_side",
    "library",
    "pathway",
    "p_value",
    "q_mofa_within_side",
    "q_both_sides_within_library_factor"
  )

  assert(
    all(needed_enrichment %in% names(enrichment)),
    paste(
      "03_MOFA_NATIVE_ENRICHMENT_ALL.csv no contiene:",
      paste(setdiff(needed_enrichment, names(enrichment)), collapse = ", ")
    )
  )

  bridge <- merge(
    enrichment,
    side_class,
    by = c("factor", "oriented_side"),
    all.x = TRUE,
    sort = FALSE
  )

  assert(
    all(!is.na(bridge$assigned_class_nominal)),
    "Hay pathways sin asignación factor-side -> clase"
  )

  bridge$pathway_evidence <- mapply(
    classify_pathway_evidence,
    q_native = bridge$q_mofa_within_side,
    q_two_sides = bridge$q_both_sides_within_library_factor,
    USE.NAMES = FALSE
  )

  bridge$class_pathway_primary <-
    bridge$pathway_evidence == "PRIMARY_TWO_SIDE_FDR"

  bridge$class_pathway_native_secondary <-
    bridge$pathway_evidence == "NATIVE_WITHIN_SIDE_ONLY"

  bridge
}

make_class_pathway_summary <- function(bridge) {
  x <- bridge[bridge$class_pathway_primary, , drop = FALSE]

  if (!nrow(x)) return(data.frame())

  keys <- unique(
    x[, c("assigned_class_nominal", "factor", "oriented_side", "library")]
  )

  rows <- vector("list", nrow(keys))

  for (i in seq_len(nrow(keys))) {
    k <- keys[i, , drop = FALSE]

    z <- x[
      x$assigned_class_nominal == k$assigned_class_nominal &
        x$factor == k$factor &
        x$oriented_side == k$oriented_side &
        x$library == k$library,
      ,
      drop = FALSE
    ]

    rows[[i]] <- data.frame(
      assigned_class_nominal = k$assigned_class_nominal,
      factor = k$factor,
      oriented_side = k$oriented_side,
      library = k$library,
      n_primary_pathways = nrow(z),
      p_assigned_class = unique(z$p_assigned_class)[1L],
      assignment_strong = unique(z$assignment_strong)[1L],
      best_q_two_sides = min(
        z$q_both_sides_within_library_factor,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

safe_neglog10 <- function(x) {
  -log10(pmax(as.numeric(x), .Machine$double.xmin))
}

make_bridge_plots <- function(bridge, outdir, family) {
  if (!MAKE_PLOTS || !requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(NULL))
  }

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  x <- bridge[bridge$class_pathway_primary, , drop = FALSE]
  if (!nrow(x)) return(invisible(NULL))

  x$score <- safe_neglog10(x$q_both_sides_within_library_factor)

  key <- interaction(
    x$assigned_class_nominal,
    x$factor,
    x$library,
    drop = TRUE,
    lex.order = TRUE
  )

  idx_keep <- unlist(
    lapply(
      split(seq_len(nrow(x)), key),
      function(idx) {
        ord <- order(
          x$q_both_sides_within_library_factor[idx],
          decreasing = FALSE
        )

        idx[
          ord[
            seq_len(
              min(TOP_PATHWAYS_PLOT, length(ord))
            )
          ]
        ]
      }
    ),
    use.names = FALSE
  )

  xp <- x[idx_keep, , drop = FALSE]

  xp$pathway_label <- gsub(
    "^HALLMARK_|^REACTOME_|^KEGG_MEDICUS_|^GOBP_|^GOCC_|^GOMF_",
    "",
    xp$pathway
  )

  p <- ggplot2::ggplot(
    xp,
    ggplot2::aes(
      x = score,
      y = stats::reorder(pathway_label, score),
      shape = oriented_side
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = p_assigned_class)
    ) +
    ggplot2::facet_grid(
      assigned_class_nominal + factor ~ library,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::labs(
      title = paste(family, "| puente clase–factor–pathway"),
      x = "-log10(q; ambas caras del factor)",
      y = "Pathway",
      size = "P(clase | cara)",
      shape = "Cara orientada"
    ) +
    ggplot2::theme_bw()

  ggplot2::ggsave(
    filename = file.path(outdir, "01_CLASS_FACTOR_PATHWAY_PRIMARY.png"),
    plot = p,
    width = 14,
    height = 10,
    dpi = 170
  )

  ggplot2::ggsave(
    filename = file.path(outdir, "01_CLASS_FACTOR_PATHWAY_PRIMARY.pdf"),
    plot = p,
    width = 14,
    height = 10
  )

  invisible(NULL)
}

run_family <- function(selection_run, family, stamp) {
  message("\n============================================================")
  message("CLASS <-> FACTOR <-> PATHWAY | ", family)
  message("============================================================")

  bidir_run <- find_bidirectional_run(selection_run, family)
  latent_run <- find_latent_audit_run(bidir_run)
  enrichment_run <- find_enrichment_run(bidir_run)

  enrichment_file <- file.path(
    enrichment_run,
    "tables",
    "03_MOFA_NATIVE_ENRICHMENT_ALL.csv"
  )

  enrichment <- utils::read.csv(
    enrichment_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  assert(nrow(enrichment) > 0L, "La tabla de enriquecimiento está vacía")

  factors <- unique(as.character(enrichment$factor))
  assert(length(factors) >= 1L, "No se encontraron factores")

  side_class_rows <- list()
  pairwise_rows <- list()

  for (f in factors) {
    fit_path <- find_factor_fit(latent_run, f)

    message("Leyendo draws existentes: ", basename(fit_path))

    fit <- readRDS(fit_path)
    group_draws <- extract_group_location_draws(fit)

    side_class_rows[[f]] <-
      make_factor_side_class_posterior(
        factor_name = f,
        group_draws = group_draws
      )

    pairwise_rows[[f]] <-
      make_pairwise_side_alignment(
        factor_name = f,
        group_draws = group_draws
      )
  }

  side_class <- do.call(rbind, side_class_rows)
  pairwise <- do.call(rbind, pairwise_rows)

  bridge <- build_class_pathway_bridge(
    enrichment = enrichment,
    side_class = side_class
  )

  if (!"family" %in% names(bridge)) {
    bridge$family <- family
  }

  summary <- make_class_pathway_summary(bridge)

  out <- file.path(
    bidir_run,
    "CLASS_PATHWAY_BRIDGE",
    stamp
  )

  dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "rds"), recursive = TRUE, showWarnings = FALSE)

  input_audit <- data.frame(
    family = family,
    bidirectional_run = bidir_run,
    latent_audit_run = latent_run,
    enrichment_run = enrichment_run,
    n_factors = length(factors),
    n_enrichment_rows = nrow(enrichment),
    pathway_fdr_primary = PATHWAY_FDR,
    class_assignment_strong_threshold = CLASS_ASSIGN_STRONG,
    TEST_used = FALSE,
    new_model_fitted = FALSE,
    stringsAsFactors = FALSE
  )

  write_csv(
    input_audit,
    file.path(out, "tables", "00_INPUT_AUDIT.csv")
  )

  write_csv(
    side_class,
    file.path(out, "tables", "01_FACTOR_SIDE_CLASS_POSTERIOR.csv")
  )

  write_csv(
    pairwise,
    file.path(out, "tables", "02_FACTOR_SIDE_PAIRWISE_ALIGNMENT.csv")
  )

  write_csv(
    bridge,
    file.path(out, "tables", "03_CLASS_PATHWAY_BRIDGE_ALL.csv")
  )

  write_csv(
    bridge[bridge$class_pathway_primary, , drop = FALSE],
    file.path(out, "tables", "04_CLASS_PATHWAY_BRIDGE_PRIMARY_Q10.csv")
  )

  write_csv(
    bridge[bridge$class_pathway_native_secondary, , drop = FALSE],
    file.path(out, "tables", "05_CLASS_PATHWAY_BRIDGE_NATIVE_ONLY_Q10.csv")
  )

  write_csv(
    summary,
    file.path(out, "tables", "06_CLASS_PATHWAY_SUMMARY.csv")
  )

  saveRDS(
    list(
      family = family,
      bidirectional_run = bidir_run,
      latent_audit_run = latent_run,
      enrichment_run = enrichment_run,
      input_audit = input_audit,
      factor_side_class_posterior = side_class,
      factor_side_pairwise_alignment = pairwise,
      class_pathway_bridge = bridge,
      class_pathway_summary = summary,
      parameters = list(
        pathway_fdr = PATHWAY_FDR,
        class_assign_strong = CLASS_ASSIGN_STRONG
      )
    ),
    file.path(out, "rds", "CLASS_PATHWAY_BRIDGE_OBJECT.rds"),
    version = 3
  )

  make_bridge_plots(
    bridge = bridge,
    outdir = file.path(out, "plots"),
    family = family
  )

  readme <- c(
    paste("METATEST |", family, "| PUENTE CLASE <-> FACTOR <-> PATHWAY"),
    "",
    paste("LATENT_SUPERVISED_AUDIT:", latent_run),
    paste("MOFA_NATIVE_ENRICHMENT:", enrichment_run),
    "",
    "OBJETIVO:",
    "- Traducir las caras orientadas de cada factor MOFA a NP/FD/SO.",
    "- Heredar esa asignación a los pathways enriquecidos en cada cara.",
    "",
    "NO SE REENTRENA NADA:",
    "- Los draws BRMS se leen de la auditoría latente ya terminada.",
    "- El enriquecimiento se lee de la Parte 1 ya terminada.",
    "- No se recalculan Z, W ni pathways.",
    "- TEST no participa.",
    "",
    "ASIGNACIÓN FACTOR-SIDE -> CLASE:",
    "- En cada draw posterior, la cara positive se asigna a la clase con",
    "  mayor localización latente y la cara negative a la clase con menor",
    "  localización latente.",
    "- Se reportan P(side=NP), P(side=FD), P(side=SO).",
    "- assigned_class_nominal es simplemente la mayor probabilidad.",
    paste0(
      "- assignment_strong = TRUE si P(clase|cara) >= ",
      CLASS_ASSIGN_STRONG,
      "."
    ),
    "- Este umbral NO elimina pathways; solo etiqueta la claridad de la asignación.",
    "",
    "EVIDENCIA DE PATHWAY:",
    paste0("- PRIMARY_TWO_SIDE_FDR: q_both_sides <= ", PATHWAY_FDR, "."),
    paste0(
      "- NATIVE_WITHIN_SIDE_ONLY: q_mofa_within_side <= ",
      PATHWAY_FDR,
      " pero no supera el ajuste conjunto de ambas caras."
    ),
    "- NOT_FDR_SUPPORTED: resto.",
    "",
    "INTERPRETACIÓN:",
    "- Un pathway asociado a FD significa que está enriquecido en una cara",
    "  del factor que la auditoría latente sitúa hacia FD.",
    "- No implica que el pathway sea exclusivo de FD ni que sea causal.",
    "- La asignación depende del factor y de la cara orientada.",
    "",
    "SALIDAS:",
    "- 01_FACTOR_SIDE_CLASS_POSTERIOR.csv",
    "- 02_FACTOR_SIDE_PAIRWISE_ALIGNMENT.csv",
    "- 03_CLASS_PATHWAY_BRIDGE_ALL.csv",
    "- 04_CLASS_PATHWAY_BRIDGE_PRIMARY_Q10.csv",
    "- 05_CLASS_PATHWAY_BRIDGE_NATIVE_ONLY_Q10.csv",
    "- 06_CLASS_PATHWAY_SUMMARY.csv",
    "",
    "SIGUIENTE PASO:",
    "- Construir el puente cross-ómico clase <-> factor <-> TX/PR/ME/CL.",
    "- El framework bidireccional se incorporará después como soporte adicional."
  )

  writeLines(
    readme,
    file.path(out, "READ_ME_FIRST.txt"),
    useBytes = TRUE
  )

  root_latest <- file.path(bidir_run, "CLASS_PATHWAY_BRIDGE")
  dir.create(root_latest, recursive = TRUE, showWarnings = FALSE)

  writeLines(
    out,
    file.path(root_latest, "LATEST_RUN.txt")
  )

  data.frame(
    family = family,
    n_factors = length(factors),
    n_pathway_rows = nrow(bridge),
    n_primary_pathway_rows =
      sum(bridge$class_pathway_primary, na.rm = TRUE),
    n_strong_class_assignments =
      sum(side_class$assignment_strong, na.rm = TRUE),
    output = out,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

main <- function() {
  assert(
    requireNamespace("posterior", quietly = TRUE),
    "Se requiere el paquete posterior"
  )

  assert(
    is.finite(PATHWAY_FDR) &&
      PATHWAY_FDR > 0 &&
      PATHWAY_FDR < 1,
    "CLASS_PATHWAY_FDR inválido"
  )

  assert(
    is.finite(CLASS_ASSIGN_STRONG) &&
      CLASS_ASSIGN_STRONG > 0.5 &&
      CLASS_ASSIGN_STRONG < 1,
    "CLASS_PATHWAY_ASSIGN_STRONG inválido"
  )

  upstream <- canonical_path(UPSTREAM_DIR, TRUE)

  selection_run <- if (nzchar(SELECTION_RUN)) {
    canonical_path(SELECTION_RUN, TRUE)
  } else {
    find_latest_selection_run(upstream)
  }

  stamp <- paste0(
    format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"),
    "_pid",
    Sys.getpid()
  )

  statuses <- list()

  for (family in FAMILIES) {
    statuses[[family]] <- tryCatch(
      run_family(
        selection_run = selection_run,
        family = family,
        stamp = stamp
      ),
      error = function(e) {
        data.frame(
          family = family,
          n_factors = NA_integer_,
          n_pathway_rows = NA_integer_,
          n_primary_pathway_rows = NA_integer_,
          n_strong_class_assignments = NA_integer_,
          output = NA_character_,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  }

  status <- do.call(rbind, statuses)

  status_file <- file.path(
    selection_run,
    paste0("CLASS_PATHWAY_BRIDGE_RUN_", stamp, ".csv")
  )

  write_csv(status, status_file)
  print(status, row.names = FALSE)

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

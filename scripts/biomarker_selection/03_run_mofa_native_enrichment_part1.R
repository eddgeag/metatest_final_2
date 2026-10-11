#!/usr/bin/env Rscript
# =============================================================================
# METATEST — PARTE 1
# Enriquecimiento nativo de MOFA por factor, con orientación latente congelada
# Version 0.1.4
#
# OBJETIVO
#   Cerrar el enlace:
#
#       FACTOR MOFA  <->  PATHWAY
#
#   usando MOFA2::run_enrichment() sobre los pesos transcriptómicos W.
#
#   El análisis de clase <-> factor ya fue realizado en:
#       LATENT_SUPERVISED_AUDIT
#
#   El framework bidireccional ya cubre:
#       clase <-> feature
#
#   Este script NO:
#     - reentrena MOFA;
#     - modifica Z o W;
#     - vuelve a seleccionar escenarios;
#     - utiliza TEST;
#     - clasifica sujetos;
#     - construye módulos cross-ómicos todavía;
#     - utiliza las grades A/B/C para seleccionar genes del GSEA.
#
# PRINCIPIO DE DIRECCIÓN
#   MOFA2::run_enrichment() utiliza los pesos W como estadístico local.
#   Con sign="positive", los pesos negativos se fijan a 0.
#   Con sign="negative", los pesos positivos se fijan a 0 y los negativos
#   se transforman a magnitud positiva.
#   Por ello se ejecutan ambas caras por separado.
#
#   El signo de un factor MOFA es arbitrario. La auditoría latente ya fijó una
#   orientación reproducible, sin usar las clases. Este script NO altera el
#   objeto MOFA; simplemente reetiqueta la cara nativa positiva/negativa como
#   cara ORIENTADA positiva/negativa según orientation_sign.
#
# FUENTES DE GENE SETS (MSigDB humano, msigdbr):
#   - HALLMARK
#   - REACTOME
#   - KEGG_MEDICUS
#   - GO_BP
#   - GO_CC
#   - GO_MF
#
# El script guarda la versión exacta de MSigDB utilizada.

# =============================================================================

options(stringsAsFactors = FALSE, digits = 17, warn = 1)

# ============================= CONFIGURACIÓN ================================
UPSTREAM_DIR <- Sys.getenv(
  "UPSTREAM_DIR",
  unset = paste0(
    "~/Documents/PhD/metatest/results_upstream_local_20260908/",
    "comparacion_modelamientos/upstream_comun"
  )
)

SELECTION_RUN <- Sys.getenv("SELECTION_RUN", unset = "")
FAMILIES <- c("P1", "P2", "P3")

TRANSCRIPTOMIC_VIEW <- Sys.getenv(
  "MOFA_GSEA_TRANSCRIPTOMIC_VIEW",
  unset = "transcriptomica"
)

VIEW_CODE_MAP <- c(
  transcriptomica = "tx",
  proteomica = "pr",
  metabolomica = "me",
  clinical = "cl"
)

SET_STATISTIC <- Sys.getenv(
  "MOFA_GSEA_SET_STATISTIC",
  unset = "mean.diff"
)

STATISTICAL_TEST <- Sys.getenv(
  "MOFA_GSEA_TEST",
  unset = "permutation"
)

NPERM <- as.integer(
  Sys.getenv("MOFA_GSEA_NPERM", unset = "5000")
)
GSEA_SEED_BASE <- as.integer(
  Sys.getenv(
    "MOFA_GSEA_SEED",
    unset = "20260918"
  )
)
MIN_SET_SIZE <- as.integer(
  Sys.getenv("MOFA_GSEA_MIN_SET_SIZE", unset = "10")
)

ALPHA <- as.numeric(
  Sys.getenv("MOFA_GSEA_ALPHA", unset = "0.10")
)

PADJ_METHOD <- Sys.getenv(
  "MOFA_GSEA_PADJ", unset = "BH"
)

TOP_FEATURES_PER_PATHWAY <- as.integer(
  Sys.getenv("MOFA_GSEA_TOP_FEATURES", unset = "10")
)

MAKE_PLOTS <- TRUE

LIBRARIES <- list(
  HALLMARK = list(collection = "H", subcollection = NULL),
  REACTOME = list(collection = "C2", subcollection = "CP:REACTOME"),
  KEGG_MEDICUS = list(collection = "C2", subcollection = "CP:KEGG_MEDICUS"),## REVISAR
  GO_BP = list(collection = "C5", subcollection = "GO:BP"),
  GO_CC = list(collection = "C5", subcollection = "GO:CC"),
  GO_MF = list(collection = "C5", subcollection = "GO:MF")
)
# ============================================================================

assert <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
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

canonical_path <- function(p, must_work = TRUE) {
  normalizePath(path.expand(p), winslash = "/", mustWork = must_work)
}



flip_side <- function(x) {
  ifelse(
    x == "positive", "negative",
    ifelse(x == "negative", "positive", NA_character_)
  )
}

# --------------------------- localizar corridas -----------------------------

find_latest_selection_run <- function(upstream) {
  base <- file.path(upstream, "downstream", "seleccion_mejor_escenario_train")
  assert(dir.exists(base), paste("No existe:", base))

  dirs <- list.dirs(base, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(
    dirs,
    "SELECTED_SCENARIOS_BY_ANALYSIS.csv"
  ))]

  assert(length(dirs) > 0L, paste(
    "No hay corridas de selección completas en",
    base
  ))

  info <- file.info(dirs)
  dirs[order(info$mtime, basename(dirs), decreasing = TRUE)][1L]
}

find_bidirectional_run <- function(selection_run, family) {
  root <- file.path(
    selection_run,
    family,
    "BIDIRECTIONAL_FRAMEWORK"
  )

  assert(dir.exists(root), paste(
    "No existe BIDIRECTIONAL_FRAMEWORK para",
    family
  ))

  latest_file <- file.path(root, "LATEST_RUN.txt")

  if (file.exists(latest_file)) {
    p <- trimws(readLines(latest_file, warn = FALSE)[1L])

    if (
      nzchar(p) &&
      dir.exists(p) &&
      file.exists(file.path(
        p,
        "rds",
        "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"
      ))
    ) {
      return(canonical_path(p, TRUE))
    }
  }

  dirs <- list.dirs(
    root,
    recursive = FALSE,
    full.names = TRUE
  )

  dirs <- dirs[file.exists(file.path(
    dirs,
    "rds",
    "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"
  ))]

  assert(length(dirs) > 0L, paste(
    "No hay corrida bidireccional completa para",
    family
  ))

  info <- file.info(dirs)

  canonical_path(
    dirs[order(
      info$mtime,
      basename(dirs),
      decreasing = TRUE
    )][1L],
    TRUE
  )
}

find_latent_audit_run <- function(bidirectional_run) {
  root <- file.path(
    bidirectional_run,
    "LATENT_SUPERVISED_AUDIT"
  )

  assert(dir.exists(root), paste(
    "No existe LATENT_SUPERVISED_AUDIT en",
    bidirectional_run
  ))

  latest_file <- file.path(root, "LATEST_RUN.txt")

  if (file.exists(latest_file)) {
    p <- trimws(readLines(latest_file, warn = FALSE)[1L])

    if (
      nzchar(p) &&
      dir.exists(p) &&
      file.exists(file.path(
        p,
        "tables",
        "01_FACTOR_ORIENTATION_AND_SCALING.csv"
      ))
    ) {
      return(canonical_path(p, TRUE))
    }
  }

  dirs <- list.dirs(
    root,
    recursive = FALSE,
    full.names = TRUE
  )

  dirs <- dirs[file.exists(file.path(
    dirs,
    "tables",
    "01_FACTOR_ORIENTATION_AND_SCALING.csv"
  ))]

  assert(length(dirs) > 0L, paste(
    "No hay auditoría latente completa en",
    root
  ))

  info <- file.info(dirs)

  canonical_path(
    dirs[order(
      info$mtime,
      basename(dirs),
      decreasing = TRUE
    )][1L],
    TRUE
  )
}

# ------------------------------- W MOFA -------------------------------------

extract_W_view <- function(model, view, factor_names = NULL) {
  
  factors <- if (is.null(factor_names)) "all" else factor_names
  
  W <- MOFA2::get_weights(
    model,
    views = view,
    factors = factors,
    scale = FALSE,
    as.data.frame = FALSE
  )[[1L]]
  
  W <- as.matrix(W)
  
  assert(
    !is.null(rownames(W)),
    "W transcriptómico no tiene nombres de features"
  )
  
  assert(
    !is.null(colnames(W)),
    "W transcriptómico no tiene nombres de factores"
  )
  
  W
}

# ------------------ compatibilidad get_data()/run_enrichment ----------------
#
# Algunos objetos RDS de METATEST conservan correctamente W/Z, pero la vista
# dentro de model@data puede estar materializada como un objeto S4 en lugar de
# la lista de matrices por grupo que MOFA2::get_data() espera.
#
# run_enrichment() llama siempre a get_data(), incluso cuando usamos
# set.statistic="mean.diff" + permutation. Para no modificar el modelo original
# ni cambiar la metodología, reconstruimos SOLO EN UNA COPIA TEMPORAL la vista
# observada TRAIN usando el ready_for_modeling del mismo escenario.
#
# No usa TEST, no usa etiquetas, no recalcula W/Z y no reentrena MOFA.

materialize_train_view_for_enrichment <- function(
  model,
  rfm,
  view
) {
  assert(
    view %in% names(VIEW_CODE_MAP),
    paste("No hay código interno para la vista", view)
  )

  code <- unname(VIEW_CODE_MAP[[view]])

  assert(
    code %in% names(rfm) &&
      is.list(rfm[[code]]) &&
      "train" %in% names(rfm[[code]]),
    paste("ready_for_modeling no contiene", code, "$train")
  )

  groups <- MOFA2::groups_names(model)

  assert(
    length(groups) == 1L,
    paste(
      "Este adaptador asume un único grupo MOFA; encontrados:",
      paste(groups, collapse = ", ")
    )
  )

  feature_names <- MOFA2::features_names(model)[[view]]
  sample_names <- MOFA2::samples_names(model)[[groups[1L]]]

  X <- as.matrix(rfm[[code]]$train)

  assert(
    !is.null(rownames(X)) && !is.null(colnames(X)),
    paste(code, "$train debe tener rownames=features y colnames=samples")
  )

  assert(
    all(feature_names %in% rownames(X)),
    paste(
      "Faltan features del modelo en ready_for_modeling:",
      sum(!feature_names %in% rownames(X))
    )
  )

  assert(
    all(sample_names %in% colnames(X)),
    paste(
      "Faltan muestras TRAIN del modelo en ready_for_modeling:",
      sum(!sample_names %in% colnames(X))
    )
  )

  X <- X[
    feature_names,
    sample_names,
    drop = FALSE
  ]

  assert(
    all(is.finite(X)),
    "La matriz TRAIN materializada contiene NA/Inf"
  )

  model_enrichment <- model

  # Estructura esperada por MOFA2::get_data():
  # model@data[[view]][[group]] = features x samples
  model_enrichment@data[[view]] <- stats::setNames(
    list(X),
    groups
  )

  # Verificación explícita antes de gastar miles de permutaciones.
  probe <- MOFA2::get_data(
    model_enrichment,
    views = view,
    as.data.frame = FALSE,
    add_intercept = FALSE
  )

  assert(
    is.list(probe) &&
      length(probe) == 1L &&
      is.list(probe[[1L]]) &&
      length(probe[[1L]]) == 1L,
    "La vista temporal sigue sin tener la estructura esperada por MOFA2::get_data()"
  )

  Y <- as.matrix(probe[[1L]][[1L]])

  assert(
    identical(rownames(Y), feature_names),
    "Features materializadas != features MOFA"
  )

  assert(
    identical(colnames(Y), sample_names),
    "Muestras materializadas != muestras MOFA"
  )

  audit <- data.frame(
    view = view,
    view_code = code,
    group = groups[1L],
    n_features = nrow(X),
    n_samples = ncol(X),
    source = "ready_for_modeling TRAIN of the frozen selected scenario",
    model_original_modified = FALSE,
    W_modified = FALSE,
    Z_modified = FALSE,
    TEST_used = FALSE,
    stringsAsFactors = FALSE
  )

  list(
    model = model_enrichment,
    audit = audit
  )
}

# ------------------------- anotación TX -> genes ----------------------------

strip_ensembl_version <- function(x) {
  x <- as.character(x)
  sub("\\..*$", "", x)
}

prepare_tx_annotation <- function(feature_map, model_features) {
  assert(is.data.frame(feature_map), "feature_map no es data.frame")

  needed <- c(
    "feature_original",
    "view"
  )
  assert(
    all(needed %in% names(feature_map)),
    paste("feature_map no contiene:", paste(needed, collapse = ", "))
  )

  fm <- feature_map[
    feature_map$view == TRANSCRIPTOMIC_VIEW |
      (
        "view_code" %in% names(feature_map) &
        feature_map$view_code == "tx"
      ),
    ,
    drop = FALSE
  ]

  assert(nrow(fm) > 0L, "No se encontraron features transcriptómicas en feature_map")

  fm$feature_original <- as.character(fm$feature_original)

  # Una fila por feature original del modelo.
  fm <- fm[!duplicated(fm$feature_original), , drop = FALSE]

  idx <- match(model_features, fm$feature_original)

  out <- data.frame(
    feature_original = model_features,
    GeneSymbol = NA_character_,
    EntrezGeneID = NA_character_,
    EnsemblID = NA_character_,
    stringsAsFactors = FALSE
  )

  for (cc in c("GeneSymbol", "EntrezGeneID", "EnsemblID")) {
    if (cc %in% names(fm)) {
      out[[cc]] <- as.character(fm[[cc]][idx])
    }
  }

  # Sanitización.
  out$GeneSymbol <- toupper(trimws(out$GeneSymbol))
  out$GeneSymbol[out$GeneSymbol %in% c("", "NA", "NAN")] <- NA_character_

  out$EntrezGeneID <- trimws(out$EntrezGeneID)
  out$EntrezGeneID[out$EntrezGeneID %in% c("", "NA", "NAN")] <- NA_character_

  out$EnsemblID <- strip_ensembl_version(trimws(out$EnsemblID))
  out$EnsemblID[out$EnsemblID %in% c("", "NA", "NAN")] <- NA_character_

  out
}

choose_mapping_id <- function(tx_annotation, msig) {
  candidates <- list(
    GeneSymbol = list(
      model = tx_annotation$GeneSymbol,
      msig = toupper(as.character(msig$gene_symbol))
    ),
    EnsemblID = list(
      model = strip_ensembl_version(tx_annotation$EnsemblID),
      msig = strip_ensembl_version(as.character(msig$ensembl_gene))
    ),
    EntrezGeneID = list(
      model = as.character(tx_annotation$EntrezGeneID),
      msig = as.character(msig$ncbi_gene)
    )
  )

  counts <- vapply(
    candidates,
    function(z) {
      sum(
        !is.na(z$model) &
          nzchar(z$model) &
          z$model %in% unique(z$msig)
      )
    },
    numeric(1)
  )

  assert(
    max(counts) > 0,
    "No existe solapamiento entre la anotación transcriptómica y MSigDB"
  )

  # En empate: GeneSymbol > Ensembl > Entrez.
  priority <- c("GeneSymbol", "EnsemblID", "EntrezGeneID")
  best <- priority[priority %in% names(counts)[counts == max(counts)]][1L]

  list(
    id_type = best,
    counts = counts,
    model_id = candidates[[best]]$model,
    msig_id = candidates[[best]]$msig
  )
}

build_feature_set_matrix <- function(
  msig,
  tx_annotation,
  mapping,
  min_size = MIN_SET_SIZE
) {
  msig <- as.data.frame(msig)

  assert(
    all(c("gs_name", "gs_description") %in% names(msig)),
    "msigdbr no devolvió gs_name/gs_description"
  )

  msig$.mapping_id <- mapping$msig_id
  msig <- msig[
    !is.na(msig$.mapping_id) &
      nzchar(msig$.mapping_id),
    ,
    drop = FALSE
  ]

  # ID elegido para cada feature MOFA.
  feature_id <- mapping$model_id
  names(feature_id) <- tx_annotation$feature_original

  mapped_feature <- !is.na(feature_id) &
    nzchar(feature_id) &
    feature_id %in% unique(msig$.mapping_id)

  mapped_features <- names(feature_id)[mapped_feature]

  assert(
    length(mapped_features) >= min_size,
    paste(
      "Muy pocas features transcriptómicas mapeadas a MSigDB:",
      length(mapped_features)
    )
  )

  pathways <- split(
    msig$.mapping_id,
    msig$gs_name
  )

  pathways <- lapply(pathways, unique)

  mat <- matrix(
    0L,
    nrow = length(pathways),
    ncol = length(mapped_features),
    dimnames = list(
      names(pathways),
      mapped_features
    )
  )

  mapped_ids <- feature_id[mapped_features]

  for (i in seq_along(pathways)) {
    hit <- mapped_ids %in% pathways[[i]]
    mat[i, hit] <- 1L
  }
  sizes_before_filter <- rowSums(mat)
  
  n_gene_sets_before_filter <-
    length(sizes_before_filter)
  
  max_gene_set_overlap_before_filter <-
    if (length(sizes_before_filter)) {
      max(sizes_before_filter, na.rm = TRUE)
    } else {
      0L
    }
  
  n_gene_sets_ge_5 <-
    sum(sizes_before_filter >= 5L)
  
  n_gene_sets_ge_10 <-
    sum(sizes_before_filter >= 10L)
  
  keep <- sizes_before_filter >= min_size
  
  mat <- mat[keep, , drop = FALSE]

  # Una colección puede existir en MSigDB pero no ser evaluable con el
  # universo transcriptómico del MOFA. No reducimos MIN_SET_SIZE
  # automáticamente para rescatarla: se devuelve una matriz de 0 filas y
  # la colección se registra como no evaluable.

  # Metadatos pathway.
  meta_cols <- intersect(
    c(
      "gs_name",
      "gs_description",
      "gs_collection",
      "gs_subcollection",
      "gs_collection_name",
      "gs_source_species",
      "db_version"
    ),
    names(msig)
  )

  meta <- unique(msig[, meta_cols, drop = FALSE])
  meta <- meta[!duplicated(meta$gs_name), , drop = FALSE]

  meta$n_features_model <- rowSums(mat)[
    match(meta$gs_name, rownames(mat))
  ]

  meta <- meta[
    meta$gs_name %in% rownames(mat),
    ,
    drop = FALSE
  ]

  list(
    matrix = mat,
    metadata = meta,
    n_model_features = nrow(tx_annotation),
    n_model_features_mapped = length(mapped_features),
    n_unique_gene_ids_mapped = length(unique(mapped_ids)),
    mapped_features = mapped_features,
    mapping_id_type = mapping$id_type,
    mapping_overlap_counts = mapping$counts,
    n_gene_sets_before_filter = n_gene_sets_before_filter,
    max_gene_set_overlap_before_filter =
      max_gene_set_overlap_before_filter,
    n_gene_sets_ge_5 = n_gene_sets_ge_5,
    n_gene_sets_ge_10 = n_gene_sets_ge_10
  )
}

# ---------------------------- MSigDB ----------------------------------------

get_msig_library <- function(spec) {
  args <- list(
    db_species = "HS",
    species = "Homo sapiens",
    collection = spec$collection
  )

  if (!is.null(spec$subcollection)) {
    args$subcollection <- spec$subcollection
  }

  x <- do.call(msigdbr::msigdbr, args)
  as.data.frame(x)
}

# ----------------------- MOFA native enrichment -----------------------------
run_enrichment_mc_corrected <- function(...) {
  
  f <- MOFA2::run_enrichment
  
  old_expr <- quote(
    results$p.values <- Reduce("+", ll) / nperm
  )
  
  new_expr <- quote(
    results$p.values <-
      (Reduce("+", ll) + 1) / (nperm + 1)
  )
  
  n_replaced <- 0L
  
  replace_expr <- function(x) {
    
    if (identical(x, old_expr)) {
      n_replaced <<- n_replaced + 1L
      return(new_expr)
    }
    
    if (
      is.call(x) ||
      is.expression(x) ||
      is.pairlist(x)
    ) {
      for (i in seq_along(x)) {
        x[[i]] <- Recall(x[[i]])
      }
    }
    
    x
  }
  
  body(f) <- replace_expr(body(f))
  
  assert(
    n_replaced == 1L,
    paste0(
      "No se pudo aplicar de forma inequívoca la corrección Monte Carlo ",
      "a MOFA2::run_enrichment(). Encontrados: ",
      n_replaced,
      ". Revisar versión de MOFA2."
    )
  )
  
  environment(f) <- asNamespace("MOFA2")
  
  f(...)
}
run_native_enrichment <- function(
  model,
  view,
  factors,
  feature_sets,
  sign_side,
  log_file
) {
  assert(sign_side %in% c("positive", "negative"), "sign_side inválido")

  result <- NULL

  # MOFA2 imprime cada permutación. Lo redirigimos a un log reproducible.
  invisible(
    capture.output(
      result <- run_enrichment_mc_corrected(
        object = model,
        view = view,
        feature.sets = feature_sets,
        factors = factors,
        set.statistic = SET_STATISTIC,
        statistical.test = STATISTICAL_TEST,
        sign = sign_side,
        min.size = MIN_SET_SIZE,
        nperm = NPERM,
        p.adj.method = PADJ_METHOD,
        alpha = ALPHA,
        verbose = TRUE
      ),
      file = log_file,
      append = TRUE,
      type = "output"
    )
  )

  assert(is.list(result), "run_enrichment no devolvió una lista")
  assert(
    all(c(
      "pval",
      "pval.adj",
      "feature.statistics",
      "set.statistics",
      "feature.sets"
    ) %in% names(result)),
    "Salida de run_enrichment incompleta"
  )

  result
}

enrichment_to_long <- function(
  res,
  family,
  library_name,
  native_side,
  orientation,
  pathway_meta
) {
  p <- as.matrix(res$pval)
  q <- as.matrix(res$pval.adj)
  s <- as.matrix(res$set.statistics)

  assert(identical(dim(p), dim(q)), "pval/pval.adj con dimensiones distintas")
  assert(identical(dim(p), dim(s)), "pval/set.statistics con dimensiones distintas")

  rows <- vector("list", ncol(p))

  for (jj in seq_len(ncol(p))) {
    f <- colnames(p)[jj]

    oo <- orientation[
      orientation$factor == f,
      ,
      drop = FALSE
    ]

    assert(nrow(oo) == 1L, paste(
      "No se encontró orientación única para",
      f
    ))

    oriented_side <- if (
      oo$orientation_sign[1] == 1
    ) {
      native_side
    } else {
      flip_side(native_side)
    }

    z <- data.frame(
      family = family,
      library = library_name,
      factor = f,
      orientation_sign = oo$orientation_sign[1],
      native_side = native_side,
      oriented_side = oriented_side,
      pathway = rownames(p),
      set_statistic = s[, jj],
      p_value = p[, jj],
      q_mofa_within_side = q[, jj],
      stringsAsFactors = FALSE
    )

    z$n_features_in_set <- rowSums(res$feature.sets)[
      match(z$pathway, rownames(res$feature.sets))
    ]

    rows[[jj]] <- z
  }

  out <- do.call(rbind, rows)

  if (!is.null(pathway_meta) && nrow(pathway_meta)) {
    pm <- pathway_meta
    names(pm)[names(pm) == "gs_name"] <- "pathway"

    out <- merge(
      out,
      pm,
      by = "pathway",
      all.x = TRUE,
      sort = FALSE
    )
  }

  out
}

add_two_side_fdr <- function(x) {
  x$q_both_sides_within_library_factor <- NA_real_

  key <- interaction(
    x$family,
    x$library,
    x$factor,
    drop = TRUE,
    lex.order = TRUE
  )

  split_idx <- split(seq_len(nrow(x)), key)

  for (idx in split_idx) {
    x$q_both_sides_within_library_factor[idx] <-
      stats::p.adjust(
        x$p_value[idx],
        method = PADJ_METHOD
      )
  }

  x
}

# -------------------------- drivers / leading features ----------------------

make_top_feature_table <- function(
  res,
  family,
  library_name,
  native_side,
  orientation,
  W_native,
  tx_annotation,
  result_long
) {
  fs <- as.matrix(res$feature.statistics)
  membership <- as.matrix(res$feature.sets)

  assert(
    all(colnames(fs) %in% colnames(W_native)),
    "Factores feature.statistics != W"
  )

  # Solo pathways con soporte nativo MOFA al alpha indicado.
  sig <- result_long[
    result_long$q_mofa_within_side <= ALPHA &
      is.finite(result_long$q_mofa_within_side),
    ,
    drop = FALSE
  ]

  if (!nrow(sig)) return(data.frame())

  rows <- list()
  rr <- 1L

  for (i in seq_len(nrow(sig))) {
    f <- sig$factor[i]
    pathway <- sig$pathway[i]

    if (
      !pathway %in% rownames(membership) ||
      !f %in% colnames(fs)
    ) next

    members <- colnames(membership)[
      membership[pathway, ] == 1
    ]

    members <- intersect(
      members,
      rownames(fs)
    )

    if (!length(members)) next

    stat <- fs[members, f]
    ord <- order(stat, decreasing = TRUE)
    members <- members[ord]
    stat <- stat[ord]

    keep_n <- min(
      TOP_FEATURES_PER_PATHWAY,
      length(members)
    )

    members <- members[seq_len(keep_n)]
    stat <- stat[seq_len(keep_n)]

    oo <- orientation[
      orientation$factor == f,
      ,
      drop = FALSE
    ]

    wnat <- W_native[members, f]
    wor <- wnat * oo$orientation_sign[1]

    ann_idx <- match(
      members,
      tx_annotation$feature_original
    )

    rows[[rr]] <- data.frame(
      family = family,
      library = library_name,
      factor = f,
      native_side = native_side,
      oriented_side = if (
        oo$orientation_sign[1] == 1
      ) native_side else flip_side(native_side),
      pathway = pathway,
      pathway_q_mofa = sig$q_mofa_within_side[i],
      pathway_q_both_sides =
        sig$q_both_sides_within_library_factor[i],
      rank_within_pathway = seq_len(keep_n),
      feature_original = members,
      GeneSymbol = tx_annotation$GeneSymbol[ann_idx],
      EntrezGeneID = tx_annotation$EntrezGeneID[ann_idx],
      EnsemblID = tx_annotation$EnsemblID[ann_idx],
      feature_statistic_abs_filtered_weight = as.numeric(stat),
      loading_native = as.numeric(wnat),
      loading_oriented = as.numeric(wor),
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

# -------------------------------- plots -------------------------------------

save_native_plots <- function(
  res,
  outdir,
  library_name,
  native_side
) {
  if (!MAKE_PLOTS) return(invisible(NULL))

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Heatmap nativo.
  png_path <- file.path(
    outdir,
    paste0(
      library_name,
      "__",
      native_side,
      "__NATIVE_HEATMAP.png"
    )
  )

  grDevices::png(
    png_path,
    width = 1800,
    height = 1300,
    res = 160
  )

  try(
    MOFA2::plot_enrichment_heatmap(
      res,
      alpha = ALPHA
    ),
    silent = TRUE
  )

  grDevices::dev.off()

  # Gráficos factor-wise nativos.
  for (f in colnames(res$pval.adj)) {
    has_sig <- any(
      is.finite(res$pval.adj[, f]) &
        res$pval.adj[, f] <= ALPHA
    )

    if (!has_sig) next

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p1 <- tryCatch(
        MOFA2::plot_enrichment(
          res,
          factor = f,
          alpha = ALPHA,
          max.pathways = 20
        ),
        error = function(e) NULL
      )

      if (!is.null(p1)) {
        ggplot2::ggsave(
          filename = file.path(
            outdir,
            paste0(
              library_name,
              "__",
              native_side,
              "__",
              f,
              "__TOP_PATHWAYS.png"
            )
          ),
          plot = p1,
          width = 10,
          height = 7,
          dpi = 160
        )
      }

      p2 <- tryCatch(
        MOFA2::plot_enrichment_detailed(
          res,
          factor = f,
          alpha = ALPHA,
          max.genes = TOP_FEATURES_PER_PATHWAY,
          max.pathways = 10
        ),
        error = function(e) NULL
      )

      if (!is.null(p2)) {
        ggplot2::ggsave(
          filename = file.path(
            outdir,
            paste0(
              library_name,
              "__",
              native_side,
              "__",
              f,
              "__DETAILED.png"
            )
          ),
          plot = p2,
          width = 12,
          height = 8,
          dpi = 160
        )
      }
    }
  }

  invisible(NULL)
}

# ----------------------------- una familia ----------------------------------

run_family <- function(
  selection_run,
  family,
  stamp,
  msig_cache
) {
  message("\n============================================================")
  message("MOFA NATIVE ENRICHMENT | ", family)
  message("============================================================")

  bidir_run <- find_bidirectional_run(
    selection_run,
    family
  )

  latent_run <- find_latent_audit_run(
    bidir_run
  )

  obj_path <- file.path(
    bidir_run,
    "rds",
    "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"
  )

  obj <- readRDS(obj_path)

  assert(
    all(c(
      "model_path",
      "rfm_path",
      "feature_map"
    ) %in% names(obj)),
    "BIDIRECTIONAL_FRAMEWORK_OBJECT incompleto"
  )

  model_path <- canonical_path(
    obj$model_path,
    TRUE
  )

  model <- readRDS(model_path)

  rfm_path <- canonical_path(
    obj$rfm_path,
    TRUE
  )

  rfm <- readRDS(rfm_path)

  materialized <- materialize_train_view_for_enrichment(
    model = model,
    rfm = rfm,
    view = TRANSCRIPTOMIC_VIEW
  )

  model_enrichment <- materialized$model
  data_compatibility_audit <- materialized$audit

  orientation_file <- file.path(
    latent_run,
    "tables",
    "01_FACTOR_ORIENTATION_AND_SCALING.csv"
  )

  orientation <- utils::read.csv(
    orientation_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  assert(
    all(c(
      "factor",
      "orientation_sign"
    ) %in% names(orientation)),
    "Tabla de orientación incompleta"
  )

  factors <- as.character(orientation$factor)

  W_native <- extract_W_view(
    model = model,
    view = TRANSCRIPTOMIC_VIEW,
    factor_names = factors
  )

  tx_annotation <- prepare_tx_annotation(
    feature_map = obj$feature_map,
    model_features = rownames(W_native)
  )

  out <- file.path(
    bidir_run,
    "MOFA_NATIVE_ENRICHMENT",
    stamp
  )

  dir.create(
    file.path(out, "tables"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    file.path(out, "rds"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    file.path(out, "plots"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    file.path(out, "logs"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  write_csv(
    tx_annotation,
    file.path(
      out,
      "tables",
      "00_TRANSCRIPTOMIC_FEATURE_ANNOTATION.csv"
    )
  )

  write_csv(
    orientation,
    file.path(
      out,
      "tables",
      "01_FACTOR_ORIENTATION_FROM_LATENT_AUDIT.csv"
    )
  )

  write_csv(
    data_compatibility_audit,
    file.path(
      out,
      "tables",
      "01B_MOFA_DATA_COMPATIBILITY_AUDIT.csv"
    )
  )

  all_results_long <- list()
  all_top_features <- list()
  library_audit <- list()
  native_objects <- list()

  rr <- 1L
  tt <- 1L
  aa <- 1L

  for (library_name in names(LIBRARIES)) {
    message("Preparando ", library_name, "...")

    msig <- msig_cache[[library_name]]

    prep <- tryCatch(
      {
        mapping <- choose_mapping_id(
          tx_annotation,
          msig
        )

        fs <- build_feature_set_matrix(
          msig = msig,
          tx_annotation = tx_annotation,
          mapping = mapping,
          min_size = MIN_SET_SIZE
        )

        list(
          ok = TRUE,
          mapping = mapping,
          fs = fs,
          error = NA_character_
        )
      },
      error = function(e) {
        list(
          ok = FALSE,
          mapping = NULL,
          fs = NULL,
          error = conditionMessage(e)
        )
      }
    )

    db_version <- if (
      "db_version" %in% names(msig)
    ) {
      paste(unique(msig$db_version), collapse = ";")
    } else {
      NA_character_
    }

    if (!isTRUE(prep$ok)) {
      library_audit[[aa]] <- data.frame(
        family = family,
        library = library_name,
        db_version = db_version,
        status = "SKIPPED_PREPARATION_ERROR",
        reason = prep$error,
        mapping_id_type = NA_character_,
        overlap_GeneSymbol = NA_real_,
        overlap_EnsemblID = NA_real_,
        overlap_EntrezGeneID = NA_real_,
        n_model_tx_features = nrow(tx_annotation),
        n_model_tx_features_mapped = NA_integer_,
        n_unique_gene_ids_mapped = NA_integer_,
        n_gene_sets_before_filter = NA_integer_,
        max_gene_set_overlap_before_filter = NA_integer_,
        n_gene_sets_ge_5 = NA_integer_,
        n_gene_sets_ge_10 = NA_integer_,
        n_gene_sets_after_min_size = 0L,
        min_set_size = MIN_SET_SIZE,
        stringsAsFactors = FALSE
      )
      
      aa <- aa + 1L
      
      message(
        "  OMITIDA ",
        library_name,
        ": ",
        prep$error
      )
      
      next
    }

    mapping <- prep$mapping
    fs <- prep$fs

    if (nrow(fs$matrix) == 0L) {
      library_audit[[aa]] <- data.frame(
        family = family,
        library = library_name,
        db_version = db_version,
        status = "SKIPPED_INSUFFICIENT_COVERAGE",
        reason = paste0(
          "0 gene sets con >= ", MIN_SET_SIZE,
          " features del universo MOFA después del mapeo"
        ),
        mapping_id_type = fs$mapping_id_type,
        overlap_GeneSymbol = unname(mapping$counts["GeneSymbol"]),
        overlap_EnsemblID = unname(mapping$counts["EnsemblID"]),
        overlap_EntrezGeneID = unname(mapping$counts["EntrezGeneID"]),
        n_model_tx_features = fs$n_model_features,
        n_model_tx_features_mapped = fs$n_model_features_mapped,
        n_unique_gene_ids_mapped = fs$n_unique_gene_ids_mapped,
        n_gene_sets_before_filter =
          fs$n_gene_sets_before_filter,
        
        max_gene_set_overlap_before_filter =
          fs$max_gene_set_overlap_before_filter,
        
        n_gene_sets_ge_5 =
          fs$n_gene_sets_ge_5,
        
        n_gene_sets_ge_10 =
          fs$n_gene_sets_ge_10,
        n_gene_sets_after_min_size = 0L,
        min_set_size = MIN_SET_SIZE,
        stringsAsFactors = FALSE
      )
      aa <- aa + 1L
      message(
        "  OMITIDA ", library_name,
        ": ningún gene set alcanza MIN_SET_SIZE=", MIN_SET_SIZE,
        ". No se modifica el umbral."
      )
      next
    }

    library_audit[[aa]] <- data.frame(
      family = family,
      library = library_name,
      db_version = db_version,
      status = "EVALUATED",
      reason = NA_character_,
      mapping_id_type = fs$mapping_id_type,
      overlap_GeneSymbol = unname(mapping$counts["GeneSymbol"]),
      overlap_EnsemblID = unname(mapping$counts["EnsemblID"]),
      overlap_EntrezGeneID = unname(mapping$counts["EntrezGeneID"]),
      n_model_tx_features = fs$n_model_features,
      n_model_tx_features_mapped = fs$n_model_features_mapped,
      n_unique_gene_ids_mapped = fs$n_unique_gene_ids_mapped,
      
      n_gene_sets_before_filter =
        fs$n_gene_sets_before_filter,
      
      max_gene_set_overlap_before_filter =
        fs$max_gene_set_overlap_before_filter,
      
      n_gene_sets_ge_5 =
        fs$n_gene_sets_ge_5,
      
      n_gene_sets_ge_10 =
        fs$n_gene_sets_ge_10,
      
      n_gene_sets_after_min_size = nrow(fs$matrix),
      min_set_size = MIN_SET_SIZE,
      stringsAsFactors = FALSE
    )
    aa <- aa + 1L

    native_objects[[library_name]] <- list()

    side_long <- list()
    side_res <- list()

    for (native_side in c("positive", "negative")) {
      message(
        "  ",
        library_name,
        " | native ",
        native_side,
        " | ",
        STATISTICAL_TEST,
        " nperm=",
        NPERM
      )

      log_file <- file.path(
        out,
        "logs",
        paste0(
          library_name,
          "__",
          native_side,
          ".log"
        )
      )
      seed_here <-
        GSEA_SEED_BASE +
        10000L * match(family, FAMILIES) +
        100L * match(library_name, names(LIBRARIES)) +
        match(native_side, c("positive", "negative"))
      
      set.seed(seed_here)
      res <- run_native_enrichment(
        model = model_enrichment,
        view = TRANSCRIPTOMIC_VIEW,
        factors = factors,
        feature_sets = fs$matrix,
        sign_side = native_side,
        log_file = log_file
      )

      side_res[[native_side]] <- res

      long <- enrichment_to_long(
        res = res,
        family = family,
        library_name = library_name,
        native_side = native_side,
        orientation = orientation,
        pathway_meta = fs$metadata
      )

      side_long[[native_side]] <- long

      save_native_plots(
        res = res,
        outdir = file.path(
          out,
          "plots",
          library_name
        ),
        library_name = library_name,
        native_side = native_side
      )
    }

    combined <- do.call(
      rbind,
      side_long
    )

    combined <- add_two_side_fdr(
      combined
    )

    all_results_long[[rr]] <- combined
    rr <- rr + 1L

    # Top features después de haber calculado q de ambas caras.
    for (native_side in c("positive", "negative")) {
      side_combined <- combined[
        combined$native_side == native_side,
        ,
        drop = FALSE
      ]

      top <- make_top_feature_table(
        res = side_res[[native_side]],
        family = family,
        library_name = library_name,
        native_side = native_side,
        orientation = orientation,
        W_native = W_native,
        tx_annotation = tx_annotation,
        result_long = side_combined
      )

      if (nrow(top)) {
        all_top_features[[tt]] <- top
        tt <- tt + 1L
      }
    }

    native_objects[[library_name]]$feature_set_matrix <- fs$matrix
    native_objects[[library_name]]$metadata <- fs$metadata
    native_objects[[library_name]]$mapping_id_type <- fs$mapping_id_type
    native_objects[[library_name]]$positive <- side_res$positive
    native_objects[[library_name]]$negative <- side_res$negative

    saveRDS(
      native_objects[[library_name]],
      file.path(
        out,
        "rds",
        paste0(
          "NATIVE_ENRICHMENT_",
          library_name,
          ".rds"
        )
      ),
      version = 3
    )
  }

  assert(
    length(all_results_long) > 0L,
    "Ninguna colección pudo evaluarse con el universo transcriptómico actual"
  )

  results_long <- do.call(
    rbind,
    all_results_long
  )

  results_long$significant_mofa_native <-
    is.finite(results_long$q_mofa_within_side) &
    results_long$q_mofa_within_side <= ALPHA

  results_long$significant_both_sides_adjusted <-
    is.finite(results_long$q_both_sides_within_library_factor) &
    results_long$q_both_sides_within_library_factor <= ALPHA

  # Orden: factor, lado orientado, librería, q.
  results_long <- results_long[
    order(
      results_long$factor,
      results_long$oriented_side,
      results_long$library,
      results_long$q_mofa_within_side,
      -results_long$set_statistic
    ),
    ,
    drop = FALSE
  ]

  top_features <- if (
    length(all_top_features)
  ) {
    do.call(rbind, all_top_features)
  } else {
    data.frame()
  }

  audit <- do.call(
    rbind,
    library_audit
  )

  write_csv(
    audit,
    file.path(
      out,
      "tables",
      "02_GENESET_MAPPING_AUDIT.csv"
    )
  )

  write_csv(
    results_long,
    file.path(
      out,
      "tables",
      "03_MOFA_NATIVE_ENRICHMENT_ALL.csv"
    )
  )

  write_csv(
    results_long[
      results_long$significant_mofa_native,
      ,
      drop = FALSE
    ],
    file.path(
      out,
      "tables",
      "04_MOFA_NATIVE_ENRICHMENT_SIGNIFICANT_Q10.csv"
    )
  )

  write_csv(
    results_long[
      results_long$significant_both_sides_adjusted,
      ,
      drop = FALSE
    ],
    file.path(
      out,
      "tables",
      "04B_MOFA_NATIVE_ENRICHMENT_SIGNIFICANT_TWO_SIDES_Q10.csv"
    )
  )

  write_csv(
    top_features,
    file.path(
      out,
      "tables",
      "05_TOP_TRANSCRIPTOMIC_FEATURES_PER_ENRICHED_PATHWAY.csv"
    )
  )

  # Copia de contexto clase-factor, sin reinterpretarlo todavía.
  latent_contrast_file <- file.path(
    latent_run,
    "tables",
    "05_TRAIN_BAYESIAN_FACTOR_CONTRASTS.csv"
  )

  if (file.exists(latent_contrast_file)) {
    file.copy(
      latent_contrast_file,
      file.path(
        out,
        "tables",
        "06_CLASS_FACTOR_CONTEXT_FROM_LATENT_AUDIT.csv"
      ),
      overwrite = TRUE
    )
  }

  saveRDS(
    list(
      family = family,
      bidirectional_run = bidir_run,
      latent_audit_run = latent_run,
      model_path = model_path,
      rfm_path = rfm_path,
      data_compatibility_audit = data_compatibility_audit,
      transcriptomic_view = TRANSCRIPTOMIC_VIEW,
      factor_orientation = orientation,
      W_transcriptomic_native = W_native,
      tx_annotation = tx_annotation,
      library_audit = audit,
      enrichment_long = results_long,
      top_features = top_features,
      parameters = list(
        set_statistic = SET_STATISTIC,
        statistical_test = STATISTICAL_TEST,
        nperm = NPERM,
        gsea_seed_base = GSEA_SEED_BASE,
        monte_carlo_pvalue_correction =
          "(b + 1) / (nperm + 1)",
        min_set_size = MIN_SET_SIZE,
        alpha = ALPHA,
        p_adj_method = PADJ_METHOD,
        libraries = LIBRARIES
      )
    ),
    file.path(
      out,
      "rds",
      "MOFA_NATIVE_ENRICHMENT_OBJECT.rds"
    ),
    version = 3
  )

  readme <- c(
    paste(
      "METATEST |",
      family,
      "| ENRIQUECIMIENTO NATIVO MOFA"
    ),
    "",
    paste("MOFA model:", model_path),
    paste("BIDIRECTIONAL_FRAMEWORK:", bidir_run),
    paste("LATENT_SUPERVISED_AUDIT:", latent_run),
    paste("Transcriptomic view:", TRANSCRIPTOMIC_VIEW),
    "",
    "PREGUNTA:",
    "- ¿Qué procesos biológicos están representados por cada factor MOFA?",
    "",
    "MÉTODO:",
    "- MOFA2::run_enrichment().",
    "- Compatibilidad técnica: la vista TRAIN observada se materializa en una copia",
    "  temporal del objeto MOFA desde ready_for_modeling del MISMO escenario, porque",
    "  algunos RDS guardan model@data en una representación S4 que get_data() no",
    "  puede indexar como lista de grupos.",
    "- Esto NO modifica W/Z, no reentrena MOFA, no usa TEST ni etiquetas de clase.",
    paste("- set.statistic =", SET_STATISTIC),
    paste("- statistical.test =", STATISTICAL_TEST),
    paste("- nperm =", NPERM),
    paste("- permutation seed base =", GSEA_SEED_BASE),
    paste(
      "- permutation p-value correction =",
      "(b + 1) / (nperm + 1)"
    ),
    paste("- min.size =", MIN_SET_SIZE),
    paste("- FDR =", PADJ_METHOD, "| alpha =", ALPHA),
    "",
    "DIRECCIÓN:",
    "- MOFA analiza por separado las magnitudes de los loadings positivos y negativos",
    "  mediante sign='positive' y sign='negative'.",
    "- Para conservar las dos caras del factor, se ejecuta sign='positive' y",
    "  sign='negative' por separado.",
    "- native_side corresponde al signo almacenado originalmente en el modelo.",
    "- oriented_side aplica la orientación congelada de la auditoría latente.",
    "- NO se modifica el objeto MOFA ni se utilizan etiquetas NP/FD/SO para",
    "  seleccionar pathways.",
    "",
    "FDR:",
    "- q_mofa_within_side es la corrección BH nativa de MOFA dentro de cada",
    "  librería, factor y análisis positive/negative.",
    "- q_both_sides_within_library_factor es una auditoría adicional que corrige",
    "  conjuntamente las dos caras del factor dentro de la misma librería.",
    "- Las colecciones MSigDB se mantienen como familias biológicas separadas;",
    "  no se aplica una única BH mezclando GO, Reactome, KEGG y Hallmark.",
    "",
    "COBERTURA DE COLECCIONES:",
    "- Una colección sin ningún gene set con MIN_SET_SIZE features mapeadas se",
    "  registra como SKIPPED_INSUFFICIENT_COVERAGE y NO aborta la familia.",
    "- MIN_SET_SIZE no se reduce automáticamente para rescatar una colección.",
    "- 02_GENESET_MAPPING_AUDIT.csv documenta qué colecciones fueron evaluables.",
    "",
    "GENE SETS:",
    "- HALLMARK",
    "- REACTOME",
    "- KEGG_MEDICUS",
    "- GO Biological Process",
    "- GO Cellular Component",
    "- GO Molecular Function",
    "- La versión exacta de MSigDB se registra en 02_GENESET_MAPPING_AUDIT.csv.",
    "",
    "IMPORTANTE:",
    "- Este análisis es factor -> pathway.",
    "- La relación clase -> factor procede exclusivamente de la auditoría latente.",
    "- La relación clase -> feature procede del framework bidireccional.",
    "- Todavía NO se construyen módulos cross-ómicos.",
    "- TEST no participa en este análisis.",
    "",
    "SALIDAS PRINCIPALES:",
    "- 03_MOFA_NATIVE_ENRICHMENT_ALL.csv",
    "- 04_MOFA_NATIVE_ENRICHMENT_SIGNIFICANT_Q10.csv",
    "- 04B_MOFA_NATIVE_ENRICHMENT_SIGNIFICANT_TWO_SIDES_Q10.csv",
    "- 05_TOP_TRANSCRIPTOMIC_FEATURES_PER_ENRICHED_PATHWAY.csv",
    "- 06_CLASS_FACTOR_CONTEXT_FROM_LATENT_AUDIT.csv"
  )

  writeLines(
    readme,
    file.path(out, "READ_ME_FIRST.txt"),
    useBytes = TRUE
  )

  dir.create(
    file.path(
      bidir_run,
      "MOFA_NATIVE_ENRICHMENT"
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )

  writeLines(
    out,
    file.path(
      bidir_run,
      "MOFA_NATIVE_ENRICHMENT",
      "LATEST_RUN.txt"
    )
  )

  data.frame(
    family = family,
    n_factors = length(factors),
    n_libraries = length(LIBRARIES),
    n_enrichment_rows = nrow(results_long),
    n_significant_native_q10 =
      sum(results_long$significant_mofa_native, na.rm = TRUE),
    n_significant_two_sides_q10 =
      sum(results_long$significant_both_sides_adjusted, na.rm = TRUE),
    output = out,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------- MAIN ------------------------------------

main <- function() {
  assert(
    requireNamespace("MOFA2", quietly = TRUE),
    "Se requiere MOFA2"
  )

  assert(
    requireNamespace("msigdbr", quietly = TRUE),
    paste(
      "Se requiere msigdbr.",
      "Instalar con install.packages('msigdbr')"
    )
  )

  assert(
    SET_STATISTIC %in% c("mean.diff", "rank.sum"),
    "MOFA_GSEA_SET_STATISTIC debe ser mean.diff o rank.sum"
  )

  assert(
    STATISTICAL_TEST %in% c(
      "parametric",
      "cor.adj.parametric",
      "permutation"
    ),
    "MOFA_GSEA_TEST inválido"
  )

  assert(
    NPERM >= 1000L || STATISTICAL_TEST != "permutation",
    "Para permutation, usar al menos 1000 permutaciones"
  )

  assert(
    MIN_SET_SIZE >= 5L,
    "MIN_SET_SIZE demasiado pequeño"
  )

  assert(
    is.finite(ALPHA) &&
      ALPHA > 0 &&
      ALPHA < 1,
    "ALPHA inválido"
  )

  upstream <- canonical_path(
    UPSTREAM_DIR,
    TRUE
  )

  selection_run <- if (
    nzchar(SELECTION_RUN)
  ) {
    canonical_path(SELECTION_RUN, TRUE)
  } else {
    find_latest_selection_run(upstream)
  }

  stamp <- paste0(
    format(
      Sys.time(),
      "%Y%m%dT%H%M%SZ",
      tz = "UTC"
    ),
    "_pid",
    Sys.getpid()
  )

  # -------------------------------------------------------------------------
  # Cargar MSigDB una única vez para P1/P2/P3.
  # Esto fija la misma versión y las mismas colecciones en los tres análisis.
  # -------------------------------------------------------------------------
  message("Cargando MSigDB humano mediante msigdbr...")

  msig_cache <- list()

  for (nm in names(LIBRARIES)) {
    message("  ", nm)
    msig_cache[[nm]] <- get_msig_library(
      LIBRARIES[[nm]]
    )
  }

  versions <- unique(unlist(lapply(
    msig_cache,
    function(x) {
      if ("db_version" %in% names(x)) {
        unique(as.character(x$db_version))
      } else {
        NA_character_
      }
    }
  )))

  message(
    "MSigDB version(es): ",
    paste(versions, collapse = ", ")
  )

  statuses <- list()

  for (family in FAMILIES) {
    statuses[[family]] <- tryCatch(
      run_family(
        selection_run = selection_run,
        family = family,
        stamp = stamp,
        msig_cache = msig_cache
      ),
      error = function(e) {
        data.frame(
          family = family,
          n_factors = NA_integer_,
          n_libraries = NA_integer_,
          n_enrichment_rows = NA_integer_,
          n_significant_native_q10 = NA_integer_,
          n_significant_two_sides_q10 = NA_integer_,
          output = NA_character_,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  }

  status <- do.call(rbind, statuses)

  root_status <- file.path(
    selection_run,
    paste0(
      "MOFA_NATIVE_ENRICHMENT_RUN_",
      stamp,
      ".csv"
    )
  )

  write_csv(
    status,
    root_status
  )

  print(
    status,
    row.names = FALSE
  )

  if (
    "error" %in% names(status) &&
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
